use crate::basic::seed::SeedOffset;
use crate::basic::value::BlockId;

pub use super::hit_buffer::{HitBuffer, HitBufferMode, Writer};

/// A seed hit — a match between a query seed and a subject seed.
///
/// This is the fundamental unit flowing through the search pipeline:
/// seeds are matched, hits are created, filtered by ungapped extension,
/// then passed to the gapped alignment phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    /// Query sequence block ID.
    pub query: BlockId,
    /// Subject sequence position (packed).
    pub subject: u64,
    /// Seed offset within the subject sequence.
    pub seed_offset: SeedOffset,
    /// Ungapped extension score (0 if not yet scored).
    pub score: u16,
}

impl Hit {
    pub fn new(query: BlockId, subject: u64, seed_offset: SeedOffset) -> Self {
        Hit {
            query,
            subject,
            seed_offset,
            score: 0,
        }
    }

    pub fn with_score(query: BlockId, subject: u64, seed_offset: SeedOffset, score: u16) -> Self {
        Hit {
            query,
            subject,
            seed_offset,
            score,
        }
    }

    /// Matches C++ `Hit::operator<`.
    pub fn less_than_query(&self, rhs: &Hit) -> bool {
        self.query < rhs.query
    }

    /// Matches C++ `Hit::blank()`.
    pub fn blank(&self) -> bool {
        self.subject == 0
    }

    /// Matches C++ `Hit::operator%(unsigned)`.
    pub fn rem_context_normalized(&self, rhs: u32, query_contexts: u32) -> u32 {
        (self.query / query_contexts) % rhs
    }

    /// Matches C++ `Hit::operator/(size_t)`.
    pub fn div_context_normalized(&self, rhs: usize, query_contexts: u32) -> u32 {
        (self.query / query_contexts) / rhs as u32
    }

    /// Matches C++ `Hit::frame()`.
    pub fn frame(&self, query_contexts: u32) -> u32 {
        self.query % query_contexts
    }

    /// Matches C++ `Hit::global_diagonal()`.
    pub fn global_diagonal(&self) -> i64 {
        self.subject as i64 - self.seed_offset as i64
    }

    /// Matches C++ `Hit::query_id<d>()`.
    pub fn query_id(&self, d: u32) -> u32 {
        self.query / d
    }

    /// Matches C++ `Hit::SourceQuery`.
    pub fn source_query(&self, contexts: u32) -> i32 {
        (self.query / contexts) as i32
    }

    /// Matches C++ `Hit::CmpSubject`.
    pub fn cmp_subject(lhs: &Hit, rhs: &Hit) -> std::cmp::Ordering {
        lhs.subject
            .cmp(&rhs.subject)
            .then_with(|| lhs.query.cmp(&rhs.query))
            .then_with(|| lhs.seed_offset.cmp(&rhs.seed_offset))
    }

    /// Matches C++ `Hit::CmpQueryTarget`.
    pub fn cmp_query_target(lhs: &Hit, rhs: &Hit) -> std::cmp::Ordering {
        lhs.query
            .cmp(&rhs.query)
            .then_with(|| lhs.subject.cmp(&rhs.subject))
    }

    /// Matches C++ `Hit::CmpTargetOffset`.
    pub fn cmp_target_offset(&self, offset: usize) -> std::cmp::Ordering {
        self.subject.cmp(&(offset as u64))
    }

    /// Matches C++ `Hit::cmp_normalized_subject`.
    pub fn cmp_normalized_subject(lhs: &Hit, rhs: &Hit) -> std::cmp::Ordering {
        let x = lhs.subject + rhs.seed_offset as u64;
        let y = rhs.subject + lhs.seed_offset as u64;
        x.cmp(&y)
            .then_with(|| lhs.seed_offset.cmp(&rhs.seed_offset))
    }

    /// Matches C++ `Hit::cmp_frame`.
    pub fn cmp_frame(lhs: &Hit, rhs: &Hit, query_contexts: u32) -> std::cmp::Ordering {
        lhs.frame(query_contexts).cmp(&rhs.frame(query_contexts))
    }

    /// Compute the global diagonal index (query_pos - subject_pos).
    pub fn diagonal(&self, query_pos: i64, subject_pos: i64) -> i64 {
        query_pos - subject_pos
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hit_buffer() {
        let mut buf = HitBuffer::new(4);
        buf.push(0, Hit::new(0, 100, 5));
        buf.push(0, Hit::new(0, 200, 10));
        buf.push(2, Hit::new(2, 300, 15));

        assert_eq!(buf.total_hits(), 3);
        assert_eq!(buf.get_bin(0).len(), 2);
        assert_eq!(buf.get_bin(1).len(), 0);
        assert_eq!(buf.get_bin(2).len(), 1);
    }

    #[test]
    fn test_hit_buffer_sort() {
        let mut buf = HitBuffer::new(1);
        buf.push(0, Hit::new(0, 300, 0));
        buf.push(0, Hit::new(0, 100, 0));
        buf.push(0, Hit::new(0, 200, 0));
        buf.sort_by_subject();
        let hits = buf.get_bin(0);
        assert_eq!(hits[0].subject, 100);
        assert_eq!(hits[1].subject, 200);
        assert_eq!(hits[2].subject, 300);
    }

    #[test]
    fn test_hit_buffer_clear() {
        let mut buf = HitBuffer::new(2);
        buf.push(0, Hit::new(0, 100, 0));
        buf.push(1, Hit::new(1, 200, 0));
        buf.clear();
        assert_eq!(buf.total_hits(), 0);
        assert_eq!(buf.get_bin(0).len(), 0);
    }

    #[test]
    fn test_hit_buffer_cpp_partition_methods() {
        let mut buf = HitBuffer::with_key_partition(vec![10, 20, 35]);
        assert_eq!(buf.begin(0), 0);
        assert_eq!(buf.end(0), 10);
        assert_eq!(buf.begin(2), 20);
        assert_eq!(buf.end(2), 35);
        assert_eq!(buf.bins(), 3);
        assert_eq!(buf.bin(0).unwrap(), 0);
        assert_eq!(buf.bin(19).unwrap(), 1);
        assert_eq!(buf.bin(34).unwrap(), 2);
        assert!(buf.bin(35).is_err());

        buf.push(0, Hit::new(0, 100, 0));
        buf.push(1, Hit::new(1, 200, 0));
        buf.push(1, Hit::new(1, 300, 0));
        assert_eq!(buf.bin_size(1), 2);
        assert_eq!(buf.next_bin_size(), 1);
        buf.set_bins_processed(1);
        assert_eq!(buf.next_bin_size(), 2);
        buf.set_bins_processed(3);
        assert_eq!(buf.next_bin_size(), 0);
    }

    #[test]
    fn test_hit_buffer_writer_new_query_and_write() {
        let mut buf = HitBuffer::with_key_partition_and_contexts(vec![5, 10], 2);
        {
            let mut writer = buf.writer();
            writer.new_query(12, 13).unwrap();
            writer.write(12, 101, 7, 0).unwrap();
            writer.write(12, 111, 9, 0).unwrap();
            assert_eq!(writer.count(0), 0);
            assert_eq!(writer.count(1), 2);

            writer.new_query(2, 3).unwrap();
            writer.write(2, 55, 5, 0).unwrap();
            assert_eq!(writer.count(0), 1);
        }
        assert_eq!(buf.total_hits(), 3);
        assert_eq!(buf.bin_size(0), 1);
        assert_eq!(buf.bin_size(1), 2);
        assert_eq!(buf.get_bin(0)[0], Hit::with_score(2, 55, 3, 5));
        assert_eq!(buf.get_bin(1)[0], Hit::with_score(12, 101, 13, 7));
        assert_eq!(buf.get_bin(1)[1], Hit::with_score(12, 111, 13, 9));
    }

    #[test]
    fn test_hit_buffer_load_retrieve_and_finish_writing() {
        let mut buf = HitBuffer::with_key_partition(vec![10, 20]);
        buf.push(0, Hit::with_score(1, 100, 2, 5));
        buf.push(1, Hit::with_score(11, 200, 4, 6));

        assert_eq!(buf.total_disk_size(), 0);
        assert!(buf.load(0));
        let (hits, begin, end) = buf.retrieve().unwrap();
        assert_eq!((begin, end), (0, 10));
        assert_eq!(hits, &[Hit::with_score(1, 100, 2, 5)]);
        assert_eq!(buf.next_bin_size(), 1);

        assert!(buf.load(1024));
        let (hits, begin, end) = buf.retrieve().unwrap();
        assert_eq!((begin, end), (10, 20));
        assert_eq!(hits, &[Hit::with_score(11, 200, 4, 6)]);
        assert!(!buf.load(1024));
        assert!(buf.retrieve().is_none());

        buf.finish_writing();
        buf.clear();
        assert_eq!(buf.total_hits(), 0);
    }

    #[test]
    fn test_hit_cpp_accessors_and_comparators() {
        let a = Hit::with_score(7, 100, 9, 3);
        let b = Hit::with_score(8, 95, 4, 2);
        assert!(a.less_than_query(&b));
        assert!(!a.blank());
        assert!(Hit::new(0, 0, 0).blank());
        assert_eq!(a.rem_context_normalized(3, 2), 0);
        assert_eq!(a.div_context_normalized(2, 2), 1);
        assert_eq!(a.frame(3), 1);
        assert_eq!(a.global_diagonal(), 91);
        assert_eq!(a.query_id(2), 3);
        assert_eq!(a.source_query(2), 3);
        assert_eq!(Hit::cmp_subject(&a, &b), std::cmp::Ordering::Greater);
        assert_eq!(Hit::cmp_query_target(&a, &b), std::cmp::Ordering::Less);
        assert_eq!(a.cmp_target_offset(101), std::cmp::Ordering::Less);
        assert_eq!(
            Hit::cmp_normalized_subject(&a, &b),
            (100u64 + 4).cmp(&(95u64 + 9)).then_with(|| 9u32.cmp(&4))
        );
        assert_eq!(Hit::cmp_frame(&a, &b, 3), std::cmp::Ordering::Less);
    }
}
