use crate::basic::sequence::Sequence;
use crate::basic::value::{BlockId, DictId, Letter, Loc, OId, SequenceType, MASK_LETTER};
use crate::data::flags::SeqInfo;
use crate::data::seed_histogram::SeedHistogram;
use crate::data::sequence_set::{SequenceSet, StringSet};
use crate::dp::ungapped::self_score;
use crate::masking::{MaskingAlgo, MaskingTable};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence::{find_orfs, translate as translate_sequence};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BlockWrapperSeqInfo {
    pub pos: u64,
    pub seq_len: u32,
}

impl BlockWrapperSeqInfo {
    pub fn new(pos: u64, len: usize) -> Self {
        Self {
            pos,
            seq_len: len as u32,
        }
    }
}

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

    pub fn translated(&self, block_id: usize, align_mode: AlignModeBlock) -> Sequence<'_> {
        if align_mode.query_translated {
            self.seqs.get(block_id * align_mode.query_contexts).into()
        } else {
            self.seqs.get(block_id).into()
        }
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

    pub fn dict_id(
        &self,
        _block: usize,
        _block_id: BlockId,
        _format_flags: u32,
    ) -> Result<DictId, String> {
        Err("Dictionary backing is not available in data::block::Block".to_string())
    }

    pub fn soft_mask(&mut self, _algo: MaskingAlgo) {
        if self.soft_masked {
            return;
        }
        if !self.soft_masking_table.blank() {
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
        let min_len = if t[0].len() < 30 {
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
        (self.seqs.data().len() * std::mem::size_of::<Letter>()) as i64
            + self.source_seqs.letters() as i64 * std::mem::size_of::<Letter>() as i64
            + self.unmasked_seqs.letters() as i64 * std::mem::size_of::<Letter>() as i64
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

pub struct BlockWrapper<'a> {
    block: &'a Block,
    oid: OId,
}

impl<'a> BlockWrapper<'a> {
    pub fn new(block: &'a Block) -> Self {
        Self { block, oid: 0 }
    }

    pub fn file_count(&self) -> i64 {
        1
    }

    pub fn files_synced(&self) -> bool {
        true
    }

    pub fn read_seqinfo(&mut self) -> Result<BlockWrapperSeqInfo, String> {
        if self.oid >= self.block.seqs().len() as OId {
            self.oid += 1;
            return Ok(BlockWrapperSeqInfo::new(0, 0));
        }
        let l = self.block.seqs().length(self.oid as usize);
        if l == 0 {
            return Err("Database with sequence length 0 is not supported".to_string());
        }
        let out = BlockWrapperSeqInfo::new(self.oid, l as usize);
        self.oid += 1;
        Ok(out)
    }

    pub fn putback_seqinfo(&mut self) {
        self.oid -= 1;
    }

    pub fn close(&mut self) {}

    pub fn set_seqinfo_ptr(&mut self, i: OId) {
        self.oid = i;
    }

    pub fn tell_seq(&self) -> OId {
        self.oid
    }

    pub fn eof(&self) -> bool {
        self.oid >= self.block.seqs().len() as OId
    }

    pub fn init_seq_access(&mut self) {
        self.set_seqinfo_ptr(0);
    }

    pub fn read_seq(&mut self) -> Result<(Vec<Letter>, String, Option<Vec<u8>>), String> {
        Err("Operation not supported".to_string())
    }

    pub fn create_partition_balanced(&mut self, _max_letters: i64) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn save_partition(
        &mut self,
        _partition_file_name: &str,
        _annotation: &str,
    ) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn get_n_partition_chunks(&mut self) -> Result<i32, String> {
        Err("Operation not supported".to_string())
    }

    pub fn init_seqinfo_access(&mut self) {}

    pub fn seek_chunk(
        &mut self,
        _chunk_i: i32,
        _offset: usize,
        _n_seqs: i64,
    ) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn seqid(&self, oid: OId, _all: bool, _full_titles: bool) -> Result<String, String> {
        Ok(std::str::from_utf8(self.block.ids()?.get(oid as usize))
            .unwrap_or("")
            .to_string())
    }

    pub fn id_len(
        &self,
        seq_info: &BlockWrapperSeqInfo,
        _seq_info_next: &BlockWrapperSeqInfo,
    ) -> Result<usize, String> {
        Ok(self.block.ids()?.length(seq_info.pos as usize) as usize)
    }

    pub fn seek_offset(&mut self, _p: usize) {}

    pub fn read_seq_data(&self, dst: &mut Vec<Letter>, len: usize, pos: &mut usize, _seek: bool) {
        dst.clear();
        dst.reserve(len + 2);
        dst.push(Sequence::DELIMITER);
        dst.extend_from_slice(self.block.seqs().get(*pos));
        dst.push(Sequence::DELIMITER);
        *pos += 1;
    }

    pub fn read_id_data(&self, oid: i64, _all: bool, _full_titles: bool) -> Result<String, String> {
        Ok(std::str::from_utf8(self.block.ids()?.get(oid as usize))
            .unwrap_or("")
            .to_string())
    }

    pub fn skip_id_data(&mut self) {}

    pub fn sequence_count(&self) -> u64 {
        self.block.seqs().len() as u64
    }

    pub fn letters(&self) -> u64 {
        self.block.seqs().letters()
    }

    pub fn db_version(&self) -> Result<i32, String> {
        Err("Operation not supported".to_string())
    }

    pub fn program_build_version(&self) -> Result<i32, String> {
        Err("Operation not supported".to_string())
    }

    pub fn build_version(&mut self) -> Result<i32, String> {
        Err("Operation not supported".to_string())
    }

    pub fn filter_by_accession(&mut self, _file_name: &str) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn file_name(&self) -> String {
        String::new()
    }

    pub fn taxids(&self, _oid: usize) -> Result<Vec<crate::basic::value::TaxId>, String> {
        Err("Operation not supported".to_string())
    }

    pub fn seq_data(&self, _oid: usize, _dst: &mut Vec<Letter>) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn seq_length(&self, _oid: usize) -> Result<Loc, String> {
        Err("Operation not supported".to_string())
    }

    pub fn end_random_access(&mut self, _dictionary: bool) {}
}

impl<'a> From<&'a [Letter]> for Sequence<'a> {
    fn from(value: &'a [Letter]) -> Self {
        Sequence::new(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::STOP_LETTER;
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
    fn test_block_sparse_oid_error_and_unsupported_dict() {
        let mut b = Block::new();
        b.push_back(&[0], None, None, 5, SequenceType::AminoAcid, 0, true)
            .unwrap();
        b.push_back(&[0], None, None, 7, SequenceType::AminoAcid, 0, true)
            .unwrap();
        assert_eq!(
            b.oid2block_id(6).unwrap_err(),
            "Block has a sparse OId range."
        );
        assert_eq!(
            b.dict_id(0, 0, 0).unwrap_err(),
            "Dictionary backing is not available in data::block::Block"
        );
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
