pub mod blast_seg;
pub mod lambda;
pub mod motifs;
pub mod tantan;
pub mod tantan_simd;

use crate::basic::value::{Letter, LETTER_MASK, MASK_LETTER, SEED_MASK};
use crate::data::sequence_set::SequenceSet;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::enum_utils::FlagBits;
use std::collections::VecDeque;

/// Masking algorithm selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskingAlgo {
    None = 0,
    Tantan = 1,
    Seg = 2,
    Motif = 4,
}

impl MaskingAlgo {
    pub fn to_string(self) -> &'static str {
        match self {
            MaskingAlgo::None => "None",
            MaskingAlgo::Tantan => "tantan",
            MaskingAlgo::Seg => "SEG",
            MaskingAlgo::Motif => panic!("Invalid conversion from enum to string."),
        }
    }

    pub fn parse(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "tantan" | "1" => MaskingAlgo::Tantan,
            "seg" => MaskingAlgo::Seg,
            "0" | "none" => MaskingAlgo::None,
            _ => MaskingAlgo::Tantan,
        }
    }
}

impl FlagBits for MaskingAlgo {
    type Bits = u32;

    fn bits(self) -> Self::Bits {
        self as u32
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskingMode {
    None,
    Tantan,
    BlastSeg,
}

impl MaskingMode {
    pub fn to_string(self) -> &'static str {
        match self {
            MaskingMode::None => "none",
            MaskingMode::Tantan => "tantan",
            MaskingMode::BlastSeg => "seg",
        }
    }

    pub fn parse(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "0" | "none" => Ok(MaskingMode::None),
            "1" | "tantan" => Ok(MaskingMode::Tantan),
            "seg" => Ok(MaskingMode::BlastSeg),
            _ => Err(format!(
                "Invalid value for string field: {}. Permitted values: none, tantan, seg",
                s
            )),
        }
    }
}

/// Matches C++ `MaskingStat::MaskingStat()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MaskingStat {
    pub masked_letters: [u64; 3],
}

impl MaskingStat {
    pub fn new() -> Self {
        Self {
            masked_letters: [0; 3],
        }
    }

    /// Matches C++ `MaskingStat::add(MaskingAlgo, uint64_t)`.
    pub fn add(&mut self, algo: MaskingAlgo, n: u64) {
        let value = algo as u32;
        assert!(value != 0);
        self.masked_letters[value.trailing_zeros() as usize] += n;
    }

    /// Matches C++ `MaskingStat::get(MaskingAlgo)`.
    pub fn get(&self, algo: MaskingAlgo) -> u64 {
        let value = algo as u32;
        assert!(value != 0);
        self.masked_letters[value.trailing_zeros() as usize]
    }
}

impl std::ops::AddAssign for MaskingStat {
    fn add_assign(&mut self, other: Self) {
        for i in 0..self.masked_letters.len() {
            self.masked_letters[i] += other.masked_letters[i];
        }
    }
}

/// Matches C++ `Mask::Ranges::Ranges()`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Ranges {
    ranges: VecDeque<(i32, i32)>,
}

impl Ranges {
    pub fn new() -> Self {
        Self {
            ranges: VecDeque::new(),
        }
    }

    /// Matches C++ `Mask::Ranges::push_back(Loc, Loc)`.
    pub fn push_back(&mut self, begin: i32, end: i32) {
        if self.ranges.is_empty() || begin > self.ranges.back().unwrap().1 {
            self.ranges.push_back((begin, end));
        } else {
            self.ranges.back_mut().unwrap().1 = end;
        }
    }

    /// Matches C++ `Mask::Ranges::push_front(Loc)`.
    pub fn push_front(&mut self, loc: i32) {
        if !self.ranges.is_empty() && self.ranges.front().unwrap().0 == loc + 1 {
            self.ranges.front_mut().unwrap().0 = loc;
        } else {
            self.ranges.push_front((loc, loc + 1));
        }
    }

    pub fn as_slices(&self) -> (&[(i32, i32)], &[(i32, i32)]) {
        self.ranges.as_slices()
    }

    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaskingTable {
    seq_count: usize,
    masked_letters: usize,
    entry: Vec<MaskingTableEntry>,
    seqs: Vec<Vec<Letter>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MaskingTableEntry {
    block_id: usize,
    begin: i32,
}

impl MaskingTable {
    pub fn new() -> Self {
        Self {
            seq_count: 0,
            masked_letters: 0,
            entry: Vec::new(),
            seqs: Vec::new(),
        }
    }

    /// Rust mapping of C++ `MaskingTable::operator=`.
    pub fn assign_from(&mut self, other: &Self) {
        self.seq_count = other.seq_count;
        self.masked_letters = other.masked_letters;
        self.entry.clone_from(&other.entry);
        self.seqs.clone_from(&other.seqs);
    }

    /// Matches C++ `MaskingTable::add(size_t, Loc, Loc, Letter*)`.
    pub fn add(&mut self, block_id: usize, begin: i32, end: i32, seq: &mut [Letter]) {
        self.entry.push(MaskingTableEntry { block_id, begin });
        self.seqs.push(seq[begin as usize..end as usize].to_vec());
        self.seq_count += 1;
        self.masked_letters += (end - begin) as usize;
        seq[begin as usize..end as usize].fill(MASK_LETTER);
    }

    /// Matches C++ `MaskingTable::remove(SequenceSet&, Loc, bool)`.
    pub fn remove(&self, seqs: &mut SequenceSet, template_len: i32, add_bit_mask: bool) {
        for i in 0..self.entry.len() {
            let entry = self.entry[i];
            let ptr = seqs.get_mut(entry.block_id);
            let begin = entry.begin as usize;
            let len = self.seqs[i].len();
            ptr[begin..begin + len].copy_from_slice(&self.seqs[i]);
            if add_bit_mask {
                let i0 = (entry.begin - template_len + 1).max(0) as usize;
                let i1 = begin + len;
                for letter in &mut ptr[i0..i1] {
                    *letter |= SEED_MASK;
                }
            }
        }
    }

    /// Matches C++ `MaskingTable::apply(SequenceSet&)`.
    pub fn apply(&self, seqs: &mut SequenceSet) {
        for i in 0..self.entry.len() {
            let entry = self.entry[i];
            let ptr = seqs.get_mut(entry.block_id);
            let begin = entry.begin as usize;
            ptr[begin..begin + self.seqs[i].len()].fill(MASK_LETTER);
        }
    }

    /// Matches C++ `MaskingTable::blank()`.
    pub fn blank(&self) -> bool {
        self.seq_count == 0
    }

    /// Matches C++ `MaskingTable::masked_letters()`.
    pub fn masked_letters(&self) -> usize {
        self.masked_letters
    }

    /// Matches C++ `MaskingTable::mem_size()`.
    pub fn mem_size(&self) -> i64 {
        (self.entry.len() * std::mem::size_of::<MaskingTableEntry>()
            + self.seqs.iter().map(Vec::len).sum::<usize>() * std::mem::size_of::<Letter>())
            as i64
    }
}

impl Default for MaskingTable {
    fn default() -> Self {
        Self::new()
    }
}

/// Explicit configuration replacing the C++ process-wide masking globals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaskingConfig {
    pub tantan_min_mask_prob: f32,
    pub max_motif_len: usize,
}

impl Default for MaskingConfig {
    fn default() -> Self {
        Self {
            tantan_min_mask_prob: 0.9,
            max_motif_len: 30,
        }
    }
}

/// Precomputed masking state corresponding to C++ `Masking`.
pub struct Masking {
    tantan: tantan::TantanMasker,
    blast_seg: blast_seg::SegParameters,
    config: MaskingConfig,
}

impl Masking {
    /// Matches `Masking::Masking(const ScoreMatrix&)` while making global
    /// configuration explicit.
    pub fn new(score_matrix: &ScoreMatrix, config: MaskingConfig) -> Self {
        Self {
            tantan: tantan::TantanMasker::from_score_matrix(
                score_matrix,
                config.tantan_min_mask_prob,
            ),
            blast_seg: blast_seg::SegParameters::new_aa(),
            config,
        }
    }

    /// Matches `Masking::operator()`. `algo_bits` accepts ORed `MaskingAlgo`
    /// values, as the C++ enum is a flag set despite its scalar spelling.
    pub fn apply_masking(
        &self,
        seq: &mut [Letter],
        algo_bits: u32,
        block_id: usize,
        mut table: Option<&mut MaskingTable>,
    ) -> Result<MaskingStat, String> {
        let mut stats = MaskingStat::new();
        if algo_bits & MaskingAlgo::Tantan.bits() != 0 {
            if let Some(masking_table) = table.as_deref_mut() {
                let ranges = self.tantan.mask_ranges(seq);
                let (first, second) = ranges.as_slices();
                for &(begin, end) in first.iter().chain(second) {
                    masking_table.add(block_id, begin, end, seq);
                    stats.add(MaskingAlgo::Tantan, (end - begin) as u64);
                }
            } else {
                self.tantan.mask_hard(seq);
            }
        }
        if algo_bits & MaskingAlgo::Seg.bits() != 0 {
            let encoded: Vec<u8> = seq.iter().map(|&letter| letter as u8).collect();
            let ranges = blast_seg::seq_buffer_seg(&encoded, 0, Some(&self.blast_seg));
            for (begin, end) in ranges {
                if let Some(masking_table) = table.as_deref_mut() {
                    masking_table.add(block_id, begin, end + 1, seq);
                    stats.add(MaskingAlgo::Seg, (end - begin + 1) as u64);
                } else {
                    seq[begin as usize..=end as usize].fill(MASK_LETTER);
                }
            }
        }
        if algo_bits & MaskingAlgo::Motif.bits() != 0 {
            let masking_table = table
                .as_deref_mut()
                .ok_or_else(|| "Motif masking requires masking table.".to_owned())?;
            let covered = self.mask_motifs(seq, block_id, masking_table);
            stats.add(MaskingAlgo::Motif, covered as u64);
        }
        Ok(stats)
    }

    fn mask_motifs(&self, seq: &mut [Letter], block_id: usize, table: &mut MaskingTable) -> usize {
        use crate::basic::sequence::Sequence;
        use crate::util::kmer::KmerIterator;

        if seq.len() < MOTIF_LEN {
            return 0;
        }
        let motif_table = motifs::init_motif_table();
        let mut ranges = Vec::<(usize, usize)>::new();
        let mut iter = KmerIterator::<MOTIF_LEN>::new(Sequence::new(seq));
        while iter.good() {
            if motif_table.contains_kmer(iter.get()) {
                let begin = iter.offset_from_start() as usize;
                let end = begin + MOTIF_LEN;
                if let Some(last) = ranges.last_mut().filter(|last| begin <= last.1) {
                    last.1 = last.1.max(end);
                } else {
                    ranges.push((begin, end));
                }
            }
            iter.increment();
        }
        let covered = ranges.iter().map(|(begin, end)| end - begin).sum();
        if covered * 2 >= seq.len() {
            return 0;
        }
        for (begin, end) in ranges {
            if end - begin <= self.config.max_motif_len {
                table.add(block_id, begin as i32, end as i32, seq);
            }
        }
        covered
    }

    pub fn mask_bit(&self, seq: &mut [Letter]) {
        self.tantan.mask_bit(seq);
    }

    pub fn bit_to_hard_mask(&self, seq: &mut [Letter], count: &mut usize) {
        for letter in seq {
            if *letter & SEED_MASK != 0 {
                *letter = MASK_LETTER;
                *count += 1;
            }
        }
    }

    pub fn remove_bit_mask(&self, seq: &mut [Letter]) {
        remove_bit_mask(seq);
    }
}

/// Deterministic ownership-safe counterpart of C++ `mask_seqs`; independent
/// sequences retain identical results without relying on global thread state.
pub fn mask_seqs(
    seqs: &mut SequenceSet,
    masking: &Masking,
    hard_mask: bool,
    algo: MaskingAlgo,
    mut table: Option<&mut MaskingTable>,
) -> Result<MaskingStat, String> {
    let mut all = MaskingStat::new();
    if algo == MaskingAlgo::None {
        return Ok(all);
    }
    if algo == MaskingAlgo::Motif && table.is_none() {
        return Err("Motif masking requires masking table.".to_owned());
    }
    for block_id in 0..seqs.len() {
        if hard_mask {
            let stats = masking.apply_masking(
                seqs.get_mut(block_id),
                algo.bits(),
                block_id,
                table.as_deref_mut(),
            )?;
            all += stats;
        } else {
            masking.mask_bit(seqs.get_mut(block_id));
        }
    }
    Ok(all)
}

pub const MOTIF_LEN: usize = 8;

/// Soft-mask a sequence position by setting the high bit.
#[inline]
pub fn soft_mask(letter: &mut Letter) {
    *letter |= SEED_MASK;
}

/// Hard-mask a sequence position by replacing with MASK_LETTER.
#[inline]
pub fn hard_mask(letter: &mut Letter) {
    *letter = MASK_LETTER;
}

/// Check if a letter is soft-masked (high bit set).
#[inline]
pub fn is_soft_masked(letter: Letter) -> bool {
    letter & SEED_MASK != 0
}

/// Convert soft masks to hard masks in a sequence.
pub fn bit_to_hard_mask(seq: &mut [Letter]) {
    for l in seq.iter_mut() {
        if is_soft_masked(*l) {
            *l = MASK_LETTER;
        }
    }
}

/// Remove soft masks (clear high bit).
pub fn remove_bit_mask(seq: &mut [Letter]) {
    for l in seq.iter_mut() {
        *l &= LETTER_MASK;
    }
}

/// Apply masking to a sequence.
pub fn mask_sequence(seq: &mut [Letter], algo: MaskingAlgo) {
    match algo {
        MaskingAlgo::Tantan => tantan::mask_tantan(seq),
        MaskingAlgo::Seg => {
            let encoded: Vec<u8> = seq.iter().map(|&letter| letter as u8).collect();
            for (begin, end) in blast_seg::seq_buffer_seg(&encoded, 0, None) {
                seq[begin as usize..=end as usize].fill(MASK_LETTER);
            }
        }
        MaskingAlgo::Motif => {}
        MaskingAlgo::None => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_soft_mask() {
        let mut l: Letter = 5;
        assert!(!is_soft_masked(l));
        soft_mask(&mut l);
        assert!(is_soft_masked(l));
        assert_eq!(l & LETTER_MASK, 5);
    }

    #[test]
    fn test_hard_mask() {
        let mut l: Letter = 5;
        hard_mask(&mut l);
        assert_eq!(l, MASK_LETTER);
    }

    #[test]
    fn test_bit_to_hard() {
        let mut seq = vec![0i8, 5, SEED_MASK | 3, 10];
        bit_to_hard_mask(&mut seq);
        assert_eq!(seq[0], 0);
        assert_eq!(seq[1], 5);
        assert_eq!(seq[2], MASK_LETTER);
        assert_eq!(seq[3], 10);
    }

    #[test]
    fn test_seg_masks_inclusive_low_complexity_ranges() {
        let mut seq = vec![0; 20];
        mask_sequence(&mut seq, MaskingAlgo::Seg);
        assert_eq!(seq, vec![MASK_LETTER; 20]);
    }

    #[test]
    fn test_masking_algo_parse() {
        assert_eq!(MaskingAlgo::parse("tantan"), MaskingAlgo::Tantan);
        assert_eq!(MaskingAlgo::parse("none"), MaskingAlgo::None);
        assert_eq!(MaskingAlgo::parse("0"), MaskingAlgo::None);
        assert_eq!(MaskingAlgo::Seg.to_string(), "SEG");
    }

    #[test]
    #[should_panic(expected = "Invalid conversion from enum to string.")]
    fn test_masking_algo_motif_to_string_panics() {
        let _ = MaskingAlgo::Motif.to_string();
    }

    #[test]
    fn test_masking_mode_parse() {
        assert_eq!(MaskingMode::parse("0").unwrap(), MaskingMode::None);
        assert_eq!(MaskingMode::parse("none").unwrap(), MaskingMode::None);
        assert_eq!(MaskingMode::parse("1").unwrap(), MaskingMode::Tantan);
        assert_eq!(MaskingMode::parse("tantan").unwrap(), MaskingMode::Tantan);
        assert_eq!(MaskingMode::parse("seg").unwrap(), MaskingMode::BlastSeg);
        assert_eq!(MaskingMode::BlastSeg.to_string(), "seg");
        assert_eq!(
            MaskingMode::parse("x").unwrap_err(),
            "Invalid value for string field: x. Permitted values: none, tantan, seg"
        );
    }

    #[test]
    fn test_masking_stat() {
        let mut a = MaskingStat::new();
        a.add(MaskingAlgo::Tantan, 3);
        a.add(MaskingAlgo::Seg, 5);
        a.add(MaskingAlgo::Motif, 7);
        assert_eq!(a.get(MaskingAlgo::Tantan), 3);
        assert_eq!(a.get(MaskingAlgo::Seg), 5);
        assert_eq!(a.get(MaskingAlgo::Motif), 7);

        let mut b = MaskingStat::new();
        b.add(MaskingAlgo::Seg, 2);
        a += b;
        assert_eq!(a.get(MaskingAlgo::Seg), 7);
    }

    #[test]
    fn masking_seg_table_round_trip_and_stats() {
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let masking = Masking::new(&matrix, MaskingConfig::default());
        let mut seq = vec![0; 20];
        let mut table = MaskingTable::new();
        let stats = masking
            .apply_masking(&mut seq, MaskingAlgo::Seg.bits(), 0, Some(&mut table))
            .unwrap();
        assert_eq!(stats.get(MaskingAlgo::Seg), 20);
        assert_eq!(seq, vec![MASK_LETTER; 20]);

        let mut seqs = SequenceSet::new();
        seqs.push(&seq);
        table.remove(&mut seqs, 0, false);
        assert_eq!(seqs.get(0), &[0; 20]);
    }

    #[test]
    fn masking_motif_requires_table() {
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let masking = Masking::new(&matrix, MaskingConfig::default());
        let mut seq = vec![0; 20];
        assert_eq!(
            masking
                .apply_masking(&mut seq, MaskingAlgo::Motif.bits(), 0, None)
                .unwrap_err(),
            "Motif masking requires masking table."
        );
    }

    #[test]
    fn mask_seqs_soft_path_sets_and_converts_bits() {
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let masking = Masking::new(&matrix, MaskingConfig::default());
        let mut seqs = SequenceSet::new();
        seqs.push(&vec![0; 80]);
        mask_seqs(&mut seqs, &masking, false, MaskingAlgo::Seg, None).unwrap();
        assert!(seqs.get(0).iter().any(|letter| letter & SEED_MASK != 0));
        let mut count = 0;
        masking.bit_to_hard_mask(seqs.get_mut(0), &mut count);
        assert!(count > 0);
        assert_eq!(
            seqs.get(0)
                .iter()
                .filter(|letter| **letter == MASK_LETTER)
                .count(),
            count
        );
    }

    #[test]
    fn test_ranges() {
        let mut ranges = Ranges::new();
        ranges.push_back(5, 8);
        ranges.push_back(8, 10);
        ranges.push_back(12, 13);
        ranges.push_front(4);
        ranges.push_front(2);
        let (front, back) = ranges.as_slices();
        let combined: Vec<_> = front.iter().chain(back.iter()).copied().collect();
        assert_eq!(combined, vec![(2, 3), (4, 10), (12, 13)]);
    }

    #[test]
    fn test_masking_table_add_remove_apply() {
        let mut seqs = SequenceSet::new();
        seqs.push(&[1, 2, 3, 4, 5]);
        seqs.push(&[6, 7, 8, 9]);

        let mut table = MaskingTable::new();
        assert!(table.blank());
        table.add(0, 1, 4, seqs.get_mut(0));
        assert!(!table.blank());
        assert_eq!(table.masked_letters(), 3);
        assert_eq!(seqs.get(0), &[1, MASK_LETTER, MASK_LETTER, MASK_LETTER, 5]);

        table.remove(&mut seqs, 2, true);
        assert_eq!(seqs.get(0)[1] & LETTER_MASK, 2);
        assert_ne!(seqs.get(0)[0] & SEED_MASK, 0);
        assert_ne!(seqs.get(0)[3] & SEED_MASK, 0);
        assert_eq!(seqs.get(0)[4], 5);

        table.apply(&mut seqs);
        assert_eq!(
            seqs.get(0),
            &[SEED_MASK | 1, MASK_LETTER, MASK_LETTER, MASK_LETTER, 5]
        );
        assert!(table.mem_size() > 0);
        assert_eq!(MOTIF_LEN, 8);
    }
}
