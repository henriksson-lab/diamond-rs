use crate::basic::value::{DictId, OId};
use crate::output::format::FormatCode;
use crate::output::intermediate::{IntermediateRecord, OutputFormat, FINISHED};
use std::io::{self, Cursor, Read, Write};

/// `SequenceFile` lookup surface used by C++ `JoinRecord`.
pub trait JoinBlocksDb {
    fn oid(&self, target_dict_id: DictId, ref_block: i64) -> OId;
}

/// Rust translation of C++ `JoinWriter`.
#[derive(Debug)]
pub struct JoinWriter<W> {
    pub file: W,
}

impl<W> JoinWriter<W> {
    /// Matches C++ `JoinWriter::JoinWriter(file)`.
    pub fn new(file: W) -> Self {
        Self { file }
    }
}

impl<W: Write> JoinWriter<W> {
    /// Matches C++ `JoinWriter::operator()(buf)`.
    pub fn consume(&mut self, buf: &mut Vec<u8>) -> io::Result<()> {
        self.file.write_all(buf)?;
        buf.clear();
        Ok(())
    }
}

/// Rust translation of C++ `JoinRecord`.
#[derive(Debug, Clone)]
pub struct JoinRecord {
    pub block: i64,
    pub same_subject: bool,
    pub info: IntermediateRecord,
}

impl JoinRecord {
    /// Matches C++ `JoinRecord::cmp_evalue(lhs, rhs)`.
    pub fn cmp_evalue(lhs: &JoinRecord, rhs: &JoinRecord) -> bool {
        rhs.same_subject
            || (!rhs.same_subject
                && (lhs.info.evalue > rhs.info.evalue
                    || (lhs.info.evalue == rhs.info.evalue && Self::cmp_score(lhs, rhs))))
    }

    /// Matches C++ `JoinRecord::cmp_score(lhs, rhs)`.
    pub fn cmp_score(lhs: &JoinRecord, rhs: &JoinRecord) -> bool {
        rhs.same_subject
            || (!rhs.same_subject
                && (lhs.info.score < rhs.info.score
                    || (lhs.info.score == rhs.info.score && rhs.db_precedence(lhs))))
    }

    /// Matches C++ `JoinRecord::db_precedence(rhs)`.
    pub fn db_precedence(&self, rhs: &JoinRecord) -> bool {
        self.info.target_oid < rhs.info.target_oid
    }

    /// Matches C++ `JoinRecord::JoinRecord(ref_block, subject, reader, db)`.
    pub fn new<R: io::Read, D: JoinBlocksDb>(
        ref_block: i64,
        subject: DictId,
        reader: &mut R,
        db: &D,
        output_format: &OutputFormat,
    ) -> io::Result<Self> {
        let mut info = IntermediateRecord::read(reader, output_format)?;
        let same_subject = info.target_dict_id == subject;
        if output_format.code != FormatCode::Daa {
            info.target_oid = db.oid(info.target_dict_id, ref_block);
        }
        Ok(Self {
            block: ref_block,
            same_subject,
            info,
        })
    }

    /// Matches C++ `JoinRecord::push_next(block, subject, reader, records)`.
    pub fn push_next<D: JoinBlocksDb>(
        block: i64,
        subject: DictId,
        reader: &mut Cursor<Vec<u8>>,
        records: &mut Vec<JoinRecord>,
        db: &D,
        output_format: &OutputFormat,
    ) -> io::Result<bool> {
        if reader.position() as usize == reader.get_ref().len() {
            return Ok(false);
        }
        let record = Self::new(block, subject, reader, db, output_format)?;
        records.push(record);
        Ok(true)
    }
}

/// Rust translation of C++ `JoinFetcher`.
#[derive(Debug)]
pub struct JoinFetcher {
    files: Vec<Cursor<Vec<u8>>>,
    query_ids: Vec<u32>,
    query_last: u32,
    pub buf: Vec<Vec<u8>>,
    pub query_id: u32,
    pub unaligned_from: u32,
}

impl JoinFetcher {
    /// Matches C++ `JoinFetcher::init(files)`.
    pub fn init(files: Vec<Vec<u8>>) -> io::Result<Self> {
        let mut cursors = Vec::with_capacity(files.len());
        let mut query_ids = Vec::with_capacity(files.len());
        for file in files {
            let mut cursor = Cursor::new(file);
            query_ids.push(read_u32(&mut cursor)?);
            cursors.push(cursor);
        }
        let blocks = cursors.len();
        Ok(Self {
            files: cursors,
            query_ids,
            query_last: u32::MAX,
            buf: vec![Vec::new(); blocks],
            query_id: FINISHED,
            unaligned_from: 0,
        })
    }

    /// Matches C++ `JoinFetcher::finish()`.
    pub fn finish(&mut self) {
        self.files.clear();
        self.query_ids.clear();
        self.buf.clear();
    }

    /// Matches C++ `JoinFetcher::next()`.
    pub fn next(&self) -> u32 {
        self.query_ids.iter().copied().min().unwrap_or(FINISHED)
    }

    /// Matches C++ `JoinFetcher::block_count()`.
    pub fn block_count(&self) -> usize {
        self.files.len()
    }

    /// Matches C++ `JoinFetcher::fetch(b)`.
    pub fn fetch(&mut self, b: usize) -> io::Result<()> {
        let size = read_u32(&mut self.files[b])? as usize;
        self.buf[b].clear();
        self.buf[b].resize(size, 0);
        self.files[b].read_exact(&mut self.buf[b])?;
        self.query_ids[b] = read_u32(&mut self.files[b])?;
        Ok(())
    }

    /// Matches C++ `JoinFetcher::operator()()`.
    pub fn fetch_next(&mut self) -> io::Result<bool> {
        self.query_id = self.next();
        self.unaligned_from = self.query_last.wrapping_add(1);
        self.query_last = self.query_id;
        for i in 0..self.buf.len() {
            if self.query_ids[i] == self.query_id && self.query_id != FINISHED {
                self.fetch(i)?;
            } else {
                self.buf[i].clear();
            }
        }
        Ok(self.next() != FINISHED)
    }
}

/// Rust translation of C++ `BlockJoiner`.
#[derive(Debug)]
pub struct BlockJoiner {
    pub records: Vec<JoinRecord>,
    cursors: Vec<Cursor<Vec<u8>>>,
    use_score_order: bool,
}

impl BlockJoiner {
    /// Matches C++ `BlockJoiner::BlockJoiner(buf, db, output_format, use_score_order)`.
    ///
    /// The buffers are per-query record payloads. C++ strips the query id and
    /// payload size in `JoinFetcher::fetch` before constructing `BlockJoiner`.
    pub fn new<D: JoinBlocksDb>(
        buf: &[Vec<u8>],
        db: &D,
        output_format: &OutputFormat,
        use_score_order: bool,
    ) -> io::Result<Self> {
        let mut records = Vec::new();
        let mut cursors = Vec::with_capacity(buf.len());
        for (block, bytes) in buf.iter().enumerate() {
            let mut cursor = Cursor::new(bytes.clone());
            JoinRecord::push_next(
                block as i64,
                u32::MAX as DictId,
                &mut cursor,
                &mut records,
                db,
                output_format,
            )?;
            cursors.push(cursor);
        }
        Ok(Self {
            records,
            cursors,
            use_score_order,
        })
    }

    /// Matches C++ `BlockJoiner::get(target_hsp, block_idx, target_oid, db)`.
    pub fn get<D: JoinBlocksDb>(
        &mut self,
        target_hsp: &mut Vec<IntermediateRecord>,
        block_idx: &mut i64,
        target_oid: &mut OId,
        db: &D,
        output_format: &OutputFormat,
    ) -> io::Result<bool> {
        if self.records.is_empty() {
            return Ok(false);
        }

        let first_idx = self.top_record_index();
        let block = self.records[first_idx].block;
        *block_idx = block;
        *target_oid = self.records[first_idx].info.target_oid;
        let subject = self.records[first_idx].info.target_dict_id;
        target_hsp.clear();

        loop {
            if self.records.is_empty() {
                return Ok(true);
            }
            let next_idx = self.top_record_index();
            if self.records[next_idx].block != block
                || self.records[next_idx].info.target_dict_id != subject
            {
                return Ok(true);
            }

            let next = self.records.swap_remove(next_idx);
            target_hsp.push(next.info);
            let cursor = &mut self.cursors[block as usize];
            JoinRecord::push_next(block, subject, cursor, &mut self.records, db, output_format)?;
        }
    }

    fn top_record_index(&self) -> usize {
        let pred = if self.use_score_order {
            JoinRecord::cmp_score
        } else {
            JoinRecord::cmp_evalue
        };
        let mut best = 0;
        for i in 1..self.records.len() {
            if pred(&self.records[best], &self.records[i]) {
                best = i;
            }
        }
        best
    }
}

fn read_u32<R: Read>(reader: &mut R) -> io::Result<u32> {
    let mut bytes = [0u8; std::mem::size_of::<u32>()];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_ne_bytes(bytes))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::dp::swipe::HspValues;
    use crate::util::interval::Interval;

    #[derive(Debug)]
    struct TestDb;

    impl JoinBlocksDb for TestDb {
        fn oid(&self, target_dict_id: DictId, ref_block: i64) -> OId {
            (ref_block as OId) * 1000 + target_dict_id as OId
        }
    }

    fn hsp(score: i32, evalue: f64, q: i32, s: i32) -> Hsp {
        let mut hsp = Hsp::new();
        hsp.score = score;
        hsp.evalue = evalue;
        hsp.query_range = Interval::new(q, q + 3);
        hsp.query_source_range = Interval::new(q, q + 3);
        hsp.subject_range = Interval::new(s, s + 3);
        hsp.identities = 2;
        hsp.length = 3;
        hsp
    }

    fn record_bytes(format: &OutputFormat, dict: DictId, hsp: &Hsp) -> Vec<u8> {
        let mut out = Vec::new();
        IntermediateRecord::write(&mut out, hsp, 0, dict, 0, format).unwrap();
        out
    }

    #[test]
    fn test_join_writer_writes_and_clears_buffer() {
        let mut writer = JoinWriter::new(Vec::new());
        let mut buf = b"abc".to_vec();
        writer.consume(&mut buf).unwrap();
        assert_eq!(writer.file, b"abc");
        assert!(buf.is_empty());
    }

    #[test]
    fn test_join_record_comparators_match_cpp_heap_predicates() {
        let mut low = JoinRecord {
            block: 0,
            same_subject: false,
            info: IntermediateRecord::default(),
        };
        low.info.score = 10;
        low.info.evalue = 1.0e-2;
        low.info.target_oid = 10;
        let mut high = low.clone();
        high.info.score = 20;
        high.info.evalue = 1.0e-5;
        high.info.target_oid = 20;

        assert!(JoinRecord::cmp_evalue(&low, &high));
        assert!(JoinRecord::cmp_score(&low, &high));

        let mut same_subject = low.clone();
        same_subject.same_subject = true;
        assert!(JoinRecord::cmp_evalue(&high, &same_subject));
        assert!(JoinRecord::cmp_score(&high, &same_subject));
    }

    #[test]
    fn test_block_joiner_groups_same_subject_and_refills_block() {
        let format = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let db = TestDb;
        let mut b0 = Vec::new();
        b0.extend(record_bytes(&format, 7, &hsp(80, 1.0e-20, 1, 2)));
        b0.extend(record_bytes(&format, 7, &hsp(70, 1.0e-18, 5, 6)));
        let b1 = record_bytes(&format, 9, &hsp(60, 1.0e-10, 9, 10));

        let mut joiner = BlockJoiner::new(&[b0, b1], &db, &format, false).unwrap();
        let mut hsps = Vec::new();
        let mut block = -1;
        let mut oid = 0;

        assert!(joiner
            .get(&mut hsps, &mut block, &mut oid, &db, &format)
            .unwrap());
        assert_eq!(block, 0);
        assert_eq!(oid, 7);
        assert_eq!(hsps.len(), 2);
        assert_eq!(hsps[0].target_dict_id, 7);
        assert_eq!(hsps[1].target_dict_id, 7);

        assert!(joiner
            .get(&mut hsps, &mut block, &mut oid, &db, &format)
            .unwrap());
        assert_eq!(block, 1);
        assert_eq!(oid, 1009);
        assert_eq!(hsps.len(), 1);

        assert!(!joiner
            .get(&mut hsps, &mut block, &mut oid, &db, &format)
            .unwrap());
    }

    #[test]
    fn test_block_joiner_score_order_uses_score_before_evalue() {
        let format = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let db = TestDb;
        let b0 = record_bytes(&format, 5, &hsp(40, 1.0e-40, 1, 2));
        let b1 = record_bytes(&format, 6, &hsp(90, 1.0e-2, 3, 4));

        let mut joiner = BlockJoiner::new(&[b0, b1], &db, &format, true).unwrap();
        let mut hsps = Vec::new();
        let mut block = -1;
        let mut oid = 0;

        assert!(joiner
            .get(&mut hsps, &mut block, &mut oid, &db, &format)
            .unwrap());
        assert_eq!(block, 1);
        assert_eq!(oid, 1006);
        assert_eq!(hsps[0].score, 90);
    }

    fn query_buffer(query_id: u32, payload: Vec<u8>) -> Vec<u8> {
        let mut out = Vec::new();
        let seek = IntermediateRecord::write_query_intro(&mut out, query_id);
        out.extend(payload);
        IntermediateRecord::finish_query(&mut out, seek);
        out
    }

    fn join_file(records: Vec<(u32, Vec<u8>)>) -> Vec<u8> {
        let mut out = Vec::new();
        for (query_id, payload) in records {
            out.extend(query_buffer(query_id, payload));
        }
        out.extend_from_slice(&FINISHED.to_ne_bytes());
        out
    }

    #[test]
    fn test_join_fetcher_fetches_min_query_and_leaves_other_blocks_empty() {
        let format = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let db = TestDb;
        let f0 = join_file(vec![
            (42, record_bytes(&format, 7, &hsp(80, 1.0e-20, 1, 2))),
            (45, record_bytes(&format, 8, &hsp(70, 1.0e-18, 3, 4))),
        ]);
        let f1 = join_file(vec![(
            43,
            record_bytes(&format, 9, &hsp(60, 1.0e-10, 9, 10)),
        )]);

        let mut fetcher = JoinFetcher::init(vec![f0, f1]).unwrap();
        assert_eq!(fetcher.block_count(), 2);
        assert!(fetcher.fetch_next().unwrap());
        assert_eq!(fetcher.query_id, 42);
        assert_eq!(fetcher.unaligned_from, 0);
        assert!(!fetcher.buf[0].is_empty());
        assert!(fetcher.buf[1].is_empty());

        let mut joiner = BlockJoiner::new(&fetcher.buf, &db, &format, false).unwrap();
        let mut hsps = Vec::new();
        let mut block = -1;
        let mut oid = 0;
        assert!(joiner
            .get(&mut hsps, &mut block, &mut oid, &db, &format)
            .unwrap());
        assert_eq!(block, 0);
        assert_eq!(oid, 7);
        assert_eq!(hsps[0].target_dict_id, 7);

        assert!(fetcher.fetch_next().unwrap());
        assert_eq!(fetcher.query_id, 43);
        assert!(fetcher.buf[0].is_empty());
        assert!(!fetcher.buf[1].is_empty());

        assert!(!fetcher.fetch_next().unwrap());
        assert_eq!(fetcher.query_id, 45);
        assert!(!fetcher.buf[0].is_empty());
        assert!(fetcher.buf[1].is_empty());
        fetcher.finish();
        assert_eq!(fetcher.block_count(), 0);
    }

    #[test]
    fn test_join_fetcher_rejects_truncated_payload() {
        let format = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let mut file = Vec::new();
        file.extend_from_slice(&42u32.to_ne_bytes());
        file.extend_from_slice(&8u32.to_ne_bytes());
        file.extend_from_slice(&[1, 2, 3]);

        let mut fetcher = JoinFetcher::init(vec![file]).unwrap();
        let err = fetcher.fetch_next().unwrap_err().kind();
        assert_eq!(err, io::ErrorKind::UnexpectedEof);

        let file = join_file(vec![(
            42,
            record_bytes(&format, 7, &hsp(80, 1.0e-20, 1, 2)),
        )]);
        let mut fetcher = JoinFetcher::init(vec![file]).unwrap();
        fetcher.fetch_next().unwrap();
        let mut payload = fetcher.buf[0].clone();
        payload.pop();
        let err = BlockJoiner::new(&[payload], &TestDb, &format, false)
            .unwrap_err()
            .kind();
        assert_eq!(err, io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_join_fetcher_finished_only_stops_without_payload() {
        let mut fetcher = JoinFetcher::init(vec![FINISHED.to_ne_bytes().to_vec()]).unwrap();
        assert!(!fetcher.fetch_next().unwrap());
        assert_eq!(fetcher.query_id, FINISHED);
        assert!(fetcher.buf[0].is_empty());
    }
}
