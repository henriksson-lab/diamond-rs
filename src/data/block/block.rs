use crate::basic::sequence::Sequence;
use crate::basic::value::{BlockId, DictId, Letter, Loc, OId, SequenceType, MASK_LETTER};
use crate::data::flags::SeqInfo;
use crate::data::seed_histogram::SeedHistogram;
use crate::data::sequence_file::{DictionaryEntry, SequenceDictionary, SequenceFileFlags};
use crate::data::sequence_set::{SequenceSet, StringSet, TranslatedSequenceView};
use crate::dp::ungapped::self_score;
use crate::masking::{mask_sequence, MaskingAlgo, MaskingTable};
use crate::output::format::OutputFlags;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence::{all_seqids, find_orfs, seqid, translate as translate_sequence};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignModeBlock {
    pub query_translated: bool,
    pub query_contexts: usize,
}

impl AlignModeBlock {
    pub fn blastp() -> Self {
        Self {
            query_translated: false,
            query_contexts: 1,
        }
    }

    pub fn blastx() -> Self {
        Self {
            query_translated: true,
            query_contexts: 6,
        }
    }
}

pub struct Block {
    seqs: SequenceSet,
    source_seqs: SequenceSet,
    unmasked_seqs: SequenceSet,
    ids: StringSet,
    qual: StringSet,
    hst: SeedHistogram,
    block2oid: Vec<OId>,
    masked: Vec<bool>,
    self_aln_score: Vec<f64>,
    soft_masking_table: MaskingTable,
    soft_masked: bool,
    raw_bytes: u64,
}

impl Block {
    pub fn new() -> Self {
        Self {
            seqs: SequenceSet::new(),
            source_seqs: SequenceSet::new(),
            unmasked_seqs: SequenceSet::new(),
            ids: StringSet::new(),
            qual: StringSet::new(),
            hst: SeedHistogram::new(),
            block2oid: Vec::new(),
            masked: Vec::new(),
            self_aln_score: Vec::new(),
            soft_masking_table: MaskingTable::new(),
            soft_masked: false,
            raw_bytes: 0,
        }
    }

    pub fn empty(&self) -> bool {
        self.seqs.len() == 0
    }

    pub fn source_len(&self, block_id: usize, align_mode: AlignModeBlock) -> u32 {
        if align_mode.query_translated {
            self.seqs
                .reverse_translated_len(block_id * align_mode.query_contexts) as u32
        } else {
            self.seqs.length(block_id) as u32
        }
    }

    pub fn translated(
        &self,
        block_id: usize,
        align_mode: AlignModeBlock,
    ) -> TranslatedSequenceView<'_> {
        let source = if align_mode.query_translated {
            self.source_seqs.get(block_id)
        } else {
            self.seqs.get(block_id)
        };
        self.seqs.translated_seq(
            source,
            block_id * align_mode.query_contexts,
            align_mode.query_translated,
        )
    }

    pub fn long_offsets(&self) -> bool {
        self.seqs.data().len() > u32::MAX as usize
    }

    pub fn seqs(&self) -> &SequenceSet {
        &self.seqs
    }

    pub fn seqs_mut(&mut self) -> &mut SequenceSet {
        &mut self.seqs
    }

    pub fn ids(&self) -> Result<&StringSet, String> {
        if self.ids.empty() {
            Err("Block::ids()".to_string())
        } else {
            Ok(&self.ids)
        }
    }

    pub fn ids_mut(&mut self) -> Result<&mut StringSet, String> {
        if self.ids.empty() {
            Err("Block::ids()".to_string())
        } else {
            Ok(&mut self.ids)
        }
    }

    pub fn source_seqs(&self) -> &SequenceSet {
        &self.source_seqs
    }

    pub fn unmasked_seqs(&self) -> &SequenceSet {
        &self.unmasked_seqs
    }

    pub fn unmasked_seqs_mut(&mut self) -> &mut SequenceSet {
        &mut self.unmasked_seqs
    }

    pub fn qual(&self) -> &StringSet {
        &self.qual
    }

    pub fn hst(&mut self) -> &mut SeedHistogram {
        &mut self.hst
    }

    pub fn block_id2oid(&self, i: BlockId) -> OId {
        self.block2oid[i as usize]
    }

    pub fn oid_begin(&self) -> OId {
        *self.block2oid.iter().min().unwrap()
    }

    pub fn oid_end(&self) -> OId {
        self.block2oid.iter().max().unwrap() + 1
    }

    pub fn oid2block_id(&self, i: OId) -> Result<BlockId, String> {
        if self.block2oid.last().unwrap() - self.block2oid.first().unwrap()
            != self.seqs.len() as u64 - 1
        {
            return Err("Block has a sparse OId range.".to_string());
        }
        if i < *self.block2oid.first().unwrap() || i > *self.block2oid.last().unwrap() {
            return Err("OId not contained in block.".to_string());
        }
        Ok((i - self.block2oid.first().unwrap()) as BlockId)
    }

    pub fn fetch_seq_if_unmasked(&self, block_id: usize, seq: &mut Vec<Letter>) -> bool {
        if self.masked.get(block_id).copied().unwrap_or(false) {
            return false;
        }
        seq.clear();
        seq.extend_from_slice(self.seqs.get(block_id));
        true
    }

    pub fn write_masked_seq(&mut self, block_id: usize, seq: &[Letter]) {
        if self.masked.get(block_id).copied().unwrap_or(false) {
            return;
        }
        self.seqs.get_mut(block_id).copy_from_slice(seq);
        if block_id >= self.masked.len() {
            self.masked.resize(block_id + 1, false);
        }
        self.masked[block_id] = true;
    }

    /// Rust ownership adapter for C++ `Block::dict_id`.
    pub fn dict_id<F>(
        &self,
        block: usize,
        block_id: BlockId,
        dictionary: &mut SequenceDictionary,
        db_flags: SequenceFileFlags,
        format_flags: OutputFlags,
        mut database_seqid: F,
    ) -> Result<DictId, String>
    where
        F: FnMut(OId, bool, bool) -> Result<String, String>,
    {
        let block_id = block_id as usize;
        let oid = self.block_id2oid(block_id as BlockId);
        let title = if self.has_ids() {
            let title = std::str::from_utf8(self.ids.get(block_id))
                .map_err(|e| format!("Invalid UTF-8 sequence title: {e}"))?;
            if format_flags.any(OutputFlags::FULL_TITLES) {
                title.to_owned()
            } else if format_flags.any(OutputFlags::ALL_SEQIDS) {
                all_seqids(title)
            } else {
                seqid(title)
            }
        } else if format_flags.any(OutputFlags::SSEQID) {
            database_seqid(oid, format_flags.any(OutputFlags::ALL_SEQIDS), true)?
        } else {
            String::new()
        };
        let self_aln_score = if db_flags.contains(SequenceFileFlags::SELF_ALN_SCORES) {
            if !self.has_self_aln() {
                return Err("Missing self alignment scores in Block.".to_owned());
            }
            self.self_aln_score(block_id as i64)
        } else {
            0.0
        };
        dictionary.dict_id(
            block,
            block_id,
            DictionaryEntry {
                oid,
                len: self.seqs.length(block_id),
                title,
                seq: if self.unmasked_seqs.is_empty() {
                    Vec::new()
                } else {
                    self.unmasked_seqs.get(block_id).to_vec()
                },
                self_aln_score,
            },
        )
    }

    pub fn soft_mask(&mut self, algo: MaskingAlgo) {
        if self.soft_masked {
            return;
        }
        if self.soft_masking_table.blank() {
            for block_id in 0..self.seqs.len() {
                let before = self.seqs.get(block_id).to_vec();
                mask_sequence(self.seqs.get_mut(block_id), algo);
                let after = self.seqs.get(block_id).to_vec();
                let mut begin = 0;
                while begin < before.len() {
                    if before[begin] == after[begin] {
                        begin += 1;
                        continue;
                    }
                    let mut end = begin + 1;
                    while end < before.len() && before[end] != after[end] {
                        end += 1;
                    }
                    self.seqs.get_mut(block_id)[begin..end].copy_from_slice(&before[begin..end]);
                    self.soft_masking_table.add(
                        block_id,
                        begin as i32,
                        end as i32,
                        self.seqs.get_mut(block_id),
                    );
                    begin = end;
                }
            }
        } else {
            self.soft_masking_table.apply(&mut self.seqs);
        }
        self.soft_masked = true;
    }

    pub fn remove_soft_masking(&mut self, template_len: i32, add_bit_mask: bool) {
        if !self.soft_masked {
            return;
        }
        self.soft_masking_table
            .remove(&mut self.seqs, template_len, add_bit_mask);
        self.soft_masked = false;
    }

    pub fn soft_masked(&self) -> bool {
        self.soft_masked
    }

    pub fn soft_masked_letters(&self) -> usize {
        self.soft_masking_table.masked_letters()
    }

    pub fn compute_self_aln(&mut self, score_matrix: &ScoreMatrix) {
        self.self_aln_score.resize(self.seqs.len(), 0.0);
        for i in 0..self.seqs.len() {
            self.self_aln_score[i] =
                score_matrix.bitscore(self_score(self.seqs.get(i), score_matrix) as f64);
        }
    }

    pub fn self_aln_score(&self, block_id: i64) -> f64 {
        self.self_aln_score[block_id as usize]
    }

    pub fn has_self_aln(&self) -> bool {
        self.self_aln_score.len() == self.seqs.len()
    }

    pub fn push_back(
        &mut self,
        seq: &[Letter],
        id: Option<&str>,
        quals: Option<&[u8]>,
        oid: OId,
        seq_type: SequenceType,
        frame_mask: i32,
        dna_translation: bool,
    ) -> Result<i64, String> {
        self.push_back_with_min_orf(
            seq,
            id,
            quals,
            oid,
            seq_type,
            frame_mask,
            dna_translation,
            0,
            false,
        )
    }

    /// Explicit-config equivalent of the C++ global `config.min_orf_len`.
    #[allow(clippy::too_many_arguments)]
    pub fn push_back_with_min_orf(
        &mut self,
        seq: &[Letter],
        id: Option<&str>,
        quals: Option<&[u8]>,
        oid: OId,
        seq_type: SequenceType,
        frame_mask: i32,
        dna_translation: bool,
        run_len: Loc,
        frame_shift: bool,
    ) -> Result<i64, String> {
        const OVERFLOW_ERR: &str = "Sequences in block exceed supported maximum.";
        if self.block2oid.len() == BlockId::MAX as usize {
            return Err(OVERFLOW_ERR.to_string());
        }
        if let Some(id) = id {
            self.ids.push_back(id.as_bytes());
        }
        if let Some(quals) = quals {
            self.qual.push_back(quals);
        }
        self.block2oid.push(oid);
        if seq_type == SequenceType::AminoAcid || !dna_translation {
            self.seqs.push(seq);
            return Ok(seq.len() as i64);
        }
        if self.seqs.len() > BlockId::MAX as usize - 6 {
            return Err(OVERFLOW_ERR.to_string());
        }
        self.source_seqs.push(seq);
        let mut t = translate_sequence(seq);
        let min_len = if run_len != 0 {
            run_len
        } else if t[0].len() < 30 || frame_shift {
            1
        } else if t[0].len() < 100 {
            20
        } else {
            40
        };
        let mut letters = 0i64;
        for (j, frame) in t.iter_mut().enumerate() {
            if (frame_mask & (1 << j)) != 0 {
                letters += find_orfs(frame, min_len as Loc) as i64;
                self.seqs.push(frame);
            } else {
                self.seqs.fill(frame.len(), MASK_LETTER);
            }
        }
        Ok(letters)
    }

    pub fn append(&mut self, b: &Block) {
        for i in 0..b.seqs.len() {
            self.seqs.push(b.seqs.get(i));
        }
        if !b.ids.empty() {
            for i in 0..b.ids.size() as usize {
                self.ids.push_back(b.ids.get(i));
            }
        }
        self.block2oid.extend_from_slice(&b.block2oid);
    }

    pub fn seq_info(&self, id: BlockId, align_mode: AlignModeBlock) -> SeqInfo<'_> {
        let id_usize = id as usize;
        let mate_id = if id % 2 == 0 { id + 1 } else { id - 1 };
        let title = if self.ids.empty() {
            None
        } else {
            Some(std::str::from_utf8(self.ids.get(id_usize)).unwrap_or(""))
        };
        let qual = if self.qual.empty() {
            Some("")
        } else {
            Some(std::str::from_utf8(self.qual.get(id_usize)).unwrap_or(""))
        };
        SeqInfo {
            block_id: id,
            oid: self.block_id2oid(id),
            title,
            qual,
            len: if align_mode.query_translated {
                self.source_seqs.length(id_usize)
            } else {
                self.seqs.length(id_usize)
            },
            source_seq: if align_mode.query_translated {
                Sequence::new(self.source_seqs.get(id_usize))
            } else {
                Sequence::new(self.seqs.get(id_usize))
            },
            mate_seq: if align_mode.query_translated && (mate_id as usize) < self.source_seqs.len()
            {
                Sequence::new(self.source_seqs.get(mate_id as usize))
            } else {
                Sequence::empty()
            },
        }
    }

    pub fn length_sorted(&self, _threads: i32) -> Block {
        let mut lengths = self.seqs.lengths();
        lengths.sort_by(|a, b| b.cmp(a));
        let mut b = Block::new();
        for (_, j) in lengths {
            b.seqs.push(self.seqs.get(j as usize));
            if !self.ids.empty() {
                b.ids.push_back(self.ids.get(j as usize));
            }
            b.block2oid.push(self.block2oid[j as usize]);
        }
        if !self.masked.is_empty() {
            b.masked.resize(self.masked.len(), false);
        }
        b
    }

    pub fn has_ids(&self) -> bool {
        !self.ids.empty()
    }

    pub fn source_seq_count(&self) -> BlockId {
        if self.source_seqs.is_empty() {
            self.seqs.len() as BlockId
        } else {
            self.source_seqs.len() as BlockId
        }
    }

    pub fn mem_size(&self) -> i64 {
        self.seqs.mem_size() as i64
            + self.source_seqs.mem_size() as i64
            + self.unmasked_seqs.mem_size() as i64
            + self.ids.mem_size()
            + self.qual.mem_size()
            + (self.block2oid.len() * std::mem::size_of::<OId>()) as i64
            + self.soft_masking_table.mem_size()
    }

    pub fn raw_bytes(&self) -> u64 {
        self.raw_bytes
    }

    pub fn oid_count(&self) -> BlockId {
        self.block2oid.len() as BlockId
    }

    pub fn offset_oids(&mut self, offset: OId) {
        for oid in &mut self.block2oid {
            *oid += offset;
        }
    }

    pub fn set_raw_bytes(&mut self, raw_bytes: u64) {
        self.raw_bytes = raw_bytes;
    }
}

impl Default for Block {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::STOP_LETTER;
    use crate::data::block::{BlockWrapper, BlockWrapperSeqInfo};
    use crate::stats::score_matrix::ScoreMatrix;

    #[test]
    fn test_block_push_accessors_and_oids() {
        let mut b = Block::new();
        assert!(b.empty());
        assert_eq!(b.ids().unwrap_err(), "Block::ids()");
        assert_eq!(
            b.push_back(
                &[0, 1, 2],
                Some("seq1"),
                Some(b"!!!"),
                10,
                SequenceType::AminoAcid,
                0,
                true
            )
            .unwrap(),
            3
        );
        assert_eq!(
            b.push_back(
                &[3, 4, 5, 6],
                Some("seq2"),
                None,
                11,
                SequenceType::AminoAcid,
                0,
                true
            )
            .unwrap(),
            4
        );
        assert!(!b.empty());
        assert_eq!(b.source_len(1, AlignModeBlock::blastp()), 4);
        assert!(!b.long_offsets());
        assert_eq!(b.ids().unwrap().get(0), b"seq1");
        assert_eq!(b.qual().get(0), b"!!!");
        assert_eq!(b.block_id2oid(1), 11);
        assert_eq!(b.oid_begin(), 10);
        assert_eq!(b.oid_end(), 12);
        assert_eq!(b.oid2block_id(11).unwrap(), 1);
        assert_eq!(
            b.oid2block_id(13).unwrap_err(),
            "OId not contained in block."
        );
        assert_eq!(b.oid_count(), 2);
        assert!(b.has_ids());
    }

    #[test]
    fn test_block_masking_append_sort_and_seq_info() {
        let mut b = Block::new();
        b.push_back(
            &[0, 1, 2],
            Some("short"),
            None,
            20,
            SequenceType::AminoAcid,
            0,
            true,
        )
        .unwrap();
        b.push_back(
            &[3, 4, 5, 6, 7],
            Some("long"),
            None,
            21,
            SequenceType::AminoAcid,
            0,
            true,
        )
        .unwrap();

        let mut seq = Vec::new();
        assert!(b.fetch_seq_if_unmasked(0, &mut seq));
        assert_eq!(seq, vec![0, 1, 2]);
        b.write_masked_seq(0, &[9, 9, 9]);
        assert!(!b.fetch_seq_if_unmasked(0, &mut seq));
        assert_eq!(b.seqs().get(0), &[9, 9, 9]);

        let info = b.seq_info(1, AlignModeBlock::blastp());
        assert_eq!(info.oid, 21);
        assert_eq!(info.title, Some("long"));
        assert_eq!(info.qual, Some(""));
        assert_eq!(info.len, 5);

        let sorted = b.length_sorted(1);
        assert_eq!(sorted.block_id2oid(0), 21);
        assert_eq!(sorted.block_id2oid(1), 20);

        let mut appended = Block::new();
        appended.append(&sorted);
        assert_eq!(appended.oid_count(), 2);
        appended.offset_oids(100);
        assert_eq!(appended.block_id2oid(0), 121);
    }

    #[test]
    fn test_block_translation_soft_mask_and_scores() {
        let mut b = Block::new();
        let letters = b
            .push_back(
                &[0, 1, 2, 3, 0, 1, 2, 3],
                Some("dna"),
                None,
                30,
                SequenceType::Nucleotide,
                0b001111,
                true,
            )
            .unwrap();
        assert!(letters > 0);
        assert_eq!(b.source_seq_count(), 1);
        assert_eq!(b.seqs().len(), 6);
        assert_eq!(b.source_len(0, AlignModeBlock::blastx()), 8);
        let translated = b.translated(0, AlignModeBlock::blastx());
        assert_eq!(translated.source(), &[0, 1, 2, 3, 0, 1, 2, 3]);
        for frame in 0..6 {
            assert_eq!(translated.frame(frame), b.seqs().get(frame));
        }

        assert!(!b.soft_masked());
        b.soft_mask(MaskingAlgo::None);
        assert!(b.soft_masked());
        assert_eq!(b.soft_masked_letters(), 0);
        b.remove_soft_masking(0, false);
        assert!(!b.soft_masked());

        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        b.compute_self_aln(&sm);
        assert!(b.has_self_aln());
        assert!(b.self_aln_score(0).is_finite());
        b.set_raw_bytes(123);
        assert_eq!(b.raw_bytes(), 123);
        assert!(b.mem_size() > 0);
    }

    #[test]
    fn test_block_translation_applies_min_orf_masking() {
        let mut b = Block::new();
        let dna = [3, 0, 0].repeat(30);
        let letters = b
            .push_back(
                &dna,
                Some("stops"),
                None,
                31,
                SequenceType::Nucleotide,
                0b000001,
                true,
            )
            .unwrap();
        assert_eq!(letters, 0);
        assert_eq!(b.seqs().get(0), vec![STOP_LETTER; 30].as_slice());
        for frame in 1..6 {
            assert!(b.seqs().get(frame).iter().all(|&x| x == MASK_LETTER));
        }
    }

    #[test]
    fn test_block_sparse_oid_error_and_dictionary_entry() {
        let mut b = Block::new();
        b.push_back(
            &[0],
            Some("ref|one description"),
            None,
            5,
            SequenceType::AminoAcid,
            0,
            true,
        )
        .unwrap();
        b.push_back(&[0], None, None, 7, SequenceType::AminoAcid, 0, true)
            .unwrap();
        assert_eq!(
            b.oid2block_id(6).unwrap_err(),
            "Block has a sparse OId range."
        );
        let mut dictionary = SequenceDictionary::default();
        dictionary.init_block(0, 2, false);
        let id = b
            .dict_id(
                0,
                0,
                &mut dictionary,
                SequenceFileFlags::NONE,
                OutputFlags::SSEQID,
                |_, _, _| panic!("inline block title must take precedence"),
            )
            .unwrap();
        let entry = dictionary.entry(id).unwrap();
        assert_eq!(entry.oid, 5);
        assert_eq!(entry.len, 1);
        assert_eq!(entry.title, "ref|one");
        assert!(entry.seq.is_empty());
        assert_eq!(entry.self_aln_score, 0.0);
    }

    #[test]
    fn test_block_wrapper_delegates_and_reports_unsupported() {
        let mut b = Block::new();
        b.push_back(
            &[0, 1, 2],
            Some("seq1"),
            None,
            0,
            SequenceType::AminoAcid,
            0,
            true,
        )
        .unwrap();
        b.push_back(
            &[3, 4],
            Some("seq2"),
            None,
            1,
            SequenceType::AminoAcid,
            0,
            true,
        )
        .unwrap();

        let mut wrapper = BlockWrapper::new(&b);
        assert_eq!(wrapper.file_count(), 1);
        assert!(wrapper.files_synced());
        assert_eq!(wrapper.sequence_count(), 2);
        assert_eq!(wrapper.letters(), 5);
        assert_eq!(wrapper.tell_seq(), 0);
        assert!(!wrapper.eof());
        let info = wrapper.read_seqinfo().unwrap();
        assert_eq!(info, BlockWrapperSeqInfo::new(0, 3));
        assert_eq!(wrapper.tell_seq(), 1);
        wrapper.putback_seqinfo();
        assert_eq!(wrapper.tell_seq(), 0);
        assert_eq!(wrapper.seqid(1, false, false).unwrap(), "seq2");
        assert_eq!(
            wrapper
                .id_len(
                    &BlockWrapperSeqInfo::new(1, 2),
                    &BlockWrapperSeqInfo::default()
                )
                .unwrap(),
            4
        );
        let mut pos = 0usize;
        let mut dst = Vec::new();
        wrapper.read_seq_data(&mut dst, 3, &mut pos, false);
        assert_eq!(dst, vec![Sequence::DELIMITER, 0, 1, 2, Sequence::DELIMITER]);
        assert_eq!(pos, 1);
        assert_eq!(wrapper.read_id_data(0, false, false).unwrap(), "seq1");
        wrapper.set_seqinfo_ptr(2);
        assert!(wrapper.eof());
        assert_eq!(wrapper.file_name(), "");
        assert_eq!(wrapper.db_version().unwrap_err(), "Operation not supported");
        assert_eq!(wrapper.read_seq().unwrap_err(), "Operation not supported");
        assert_eq!(
            wrapper.create_partition_balanced(10).unwrap_err(),
            "Operation not supported"
        );
        wrapper.end_random_access(false);
    }
}
