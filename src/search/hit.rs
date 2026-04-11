use crate::basic::seed::SeedOffset;
use crate::basic::value::BlockId;

/// A seed hit — a match between a query seed and a subject seed.
///
/// This is the fundamental unit flowing through the search pipeline:
/// seeds are matched, hits are created, filtered by ungapped extension,
/// then passed to the gapped alignment phase.
#[derive(Debug, Clone, Copy)]
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

    /// Compute the global diagonal index (query_pos - subject_pos).
    pub fn diagonal(&self, query_pos: i64, subject_pos: i64) -> i64 {
        query_pos - subject_pos
    }
}

/// A buffer for collecting seed hits during the search phase.
///
/// Hits are partitioned into bins by query ID for efficient
/// downstream processing.
pub struct HitBuffer {
    /// Hits stored per query bin.
    bins: Vec<Vec<Hit>>,
    /// Total number of hits stored.
    total_hits: u64,
}

impl HitBuffer {
    pub fn new(num_bins: usize) -> Self {
        HitBuffer {
            bins: vec![Vec::new(); num_bins],
            total_hits: 0,
        }
    }

    /// Add a hit to the appropriate bin.
    pub fn push(&mut self, bin: usize, hit: Hit) {
        self.bins[bin].push(hit);
        self.total_hits += 1;
    }

    /// Get hits for a specific bin.
    pub fn get_bin(&self, bin: usize) -> &[Hit] {
        &self.bins[bin]
    }

    /// Get mutable hits for a specific bin.
    pub fn get_bin_mut(&mut self, bin: usize) -> &mut Vec<Hit> {
        &mut self.bins[bin]
    }

    /// Total number of hits across all bins.
    pub fn total_hits(&self) -> u64 {
        self.total_hits
    }

    /// Number of bins.
    pub fn num_bins(&self) -> usize {
        self.bins.len()
    }

    /// Clear all bins.
    pub fn clear(&mut self) {
        for bin in &mut self.bins {
            bin.clear();
        }
        self.total_hits = 0;
    }

    /// Sort hits within each bin by subject position.
    pub fn sort_by_subject(&mut self) {
        for bin in &mut self.bins {
            bin.sort_by_key(|h| h.subject);
        }
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
}
