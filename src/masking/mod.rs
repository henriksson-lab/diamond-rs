pub mod motifs;
pub mod tantan;
pub mod tantan_simd;

use crate::basic::value::{Letter, LETTER_MASK, MASK_LETTER, SEED_MASK};

/// Masking algorithm selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskingAlgo {
    None,
    Tantan,
    Seg,
}

impl MaskingAlgo {
    pub fn parse(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "tantan" | "1" => MaskingAlgo::Tantan,
            "seg" => MaskingAlgo::Seg,
            "0" | "none" => MaskingAlgo::None,
            _ => MaskingAlgo::Tantan,
        }
    }
}

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
        MaskingAlgo::Seg => {} // SEG not yet implemented
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
    fn test_masking_algo_parse() {
        assert_eq!(MaskingAlgo::parse("tantan"), MaskingAlgo::Tantan);
        assert_eq!(MaskingAlgo::parse("none"), MaskingAlgo::None);
        assert_eq!(MaskingAlgo::parse("0"), MaskingAlgo::None);
    }
}
