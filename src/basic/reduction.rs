//! Alphabet reduction implementation translated from
//! `diamond/src/basic/basic.cpp` (`Reduction::Reduction` and `decode_seed`).

use super::value::{
    Letter, AMINO_ACID_ALPHABET, DELIMITER_LETTER, MASK_LETTER, STOP_LETTER, TRUE_AA,
};

/// `Stats::blosum62.background_freqs` in the original implementation.
/// Order follows DIAMOND's canonical amino-acid alphabet.
const BLOSUM62_BACKGROUND_FREQS: [f64; TRUE_AA as usize] = [
    7.4216205067993410e-02,
    5.1614486141284638e-02,
    4.4645808512757915e-02,
    5.3626000838554413e-02,
    2.4687457167944848e-02,
    3.4259650591416023e-02,
    5.4311925684587502e-02,
    7.4146941452644999e-02,
    2.6212984805266227e-02,
    6.7917367618953756e-02,
    9.8907868497150955e-02,
    5.8155682303079680e-02,
    2.4990197579643110e-02,
    4.7418459742284751e-02,
    3.8538003320306206e-02,
    5.7229029476494421e-02,
    5.0891364550287033e-02,
    1.3029956129972148e-02,
    3.2281512313758580e-02,
    7.2919098205619245e-02,
];

/// Default reduction definition: maps 20 amino acids to 10 Murphy groups.
pub const DEFAULT_REDUCTION: &str = "A KR EDNQ C G H ILVM FYW P ST";

/// Alphabet reduction mapping amino acids to fewer groups for seed computation.
pub struct Reduction {
    map: [u32; 256],
    map8: [Letter; 256],
    map8b: [Letter; 256],
    seed_map8: [Letter; 256],
    size: u32,
    bit_size: i32,
    bit_size_exact: f64,
    freq: [f64; TRUE_AA as usize],
    /// Two Murphy-10 letters packed as one base-10 pair. The 32x32 table is
    /// small enough for L1 and lets the chunked left-most verifier replace ten
    /// independent reduction lookups with five paired lookups.
    pair10: [u16; 32 * 32],
}

impl Reduction {
    /// Create a new reduction from a definition string.
    /// Groups are space-separated; letters within a group map to the same value.
    /// Example: "A KR EDNQ C G H ILVM FYW P ST"
    pub fn new(definition_string: &str, alphabet: &[u8]) -> Self {
        let mut map = [0u32; 256];
        let mut map8 = [0i8; 256];
        let mut map8b = [0i8; 256];

        map[MASK_LETTER as u8 as usize] = MASK_LETTER as u32;
        map[STOP_LETTER as u8 as usize] = MASK_LETTER as u32;

        let tokens: Vec<&str> = definition_string.split_whitespace().collect();
        let size = tokens.len() as u32;
        let bit_size_exact = (size as f64).ln() / 2.0_f64.ln();
        let bit_size = bit_size_exact.ceil() as i32;

        let mut freq = [0.0f64; TRUE_AA as usize];

        // Build a simple char->letter lookup from the alphabet
        let mut char_to_letter = [u8::MAX; 256];
        for (i, &ch) in alphabet.iter().enumerate() {
            char_to_letter[ch as usize] = i as u8;
            char_to_letter[(ch as char).to_ascii_lowercase() as usize] = i as u8;
        }

        for (i, token) in tokens.iter().enumerate() {
            for ch in token.bytes() {
                let letter = char_to_letter[ch as usize];
                assert_ne!(
                    letter,
                    u8::MAX,
                    "Invalid character in sequence: '{}'",
                    ch as char
                );
                let letter = letter as usize;
                map[letter] = i as u32;
                map8[letter] = i as Letter;
                map8b[letter] = i as Letter;
                freq[i] += BLOSUM62_BACKGROUND_FREQS[letter];
            }
        }

        for value in &mut freq {
            *value = value.ln();
        }

        map8[MASK_LETTER as u8 as usize] = size as Letter;
        map8[STOP_LETTER as u8 as usize] = size as Letter;
        map8[DELIMITER_LETTER as u8 as usize] = size as Letter;
        map8b[MASK_LETTER as u8 as usize] = (size + 1) as Letter;
        map8b[STOP_LETTER as u8 as usize] = (size + 1) as Letter;
        map8b[DELIMITER_LETTER as u8 as usize] = (size + 1) as Letter;

        // Seed enumeration treats either a soft-mask bit or any non-amino
        // residue as the reduced MASK sentinel. Precompute the complete byte
        // mapping so whole-sequence reduction is one lookup per residue, as
        // in upstream Reduction::reduce_seq.
        let mut seed_map8 = [MASK_LETTER; 256];
        for byte in 0..256usize {
            if byte & 0x80 == 0 {
                let letter = byte & super::value::LETTER_MASK as usize;
                if letter < TRUE_AA as usize {
                    seed_map8[byte] = map8[letter];
                }
            }
        }

        let mut pair10 = [u16::MAX; 32 * 32];
        if size == 10 {
            for a in 0..32usize {
                for b in 0..32usize {
                    if super::value::is_amino_acid(a as Letter)
                        && super::value::is_amino_acid(b as Letter)
                    {
                        pair10[(a << 5) | b] = (map[a] * 10 + map[b]) as u16;
                    }
                }
            }
        }
        Reduction {
            map,
            map8,
            map8b,
            seed_map8,
            size,
            bit_size,
            bit_size_exact,
            freq,
            pair10,
        }
    }

    /// Create the default reduction.
    pub fn default_reduction() -> Self {
        Self::new(DEFAULT_REDUCTION, AMINO_ACID_ALPHABET)
    }

    pub fn size(&self) -> u32 {
        self.size
    }

    pub fn bit_size(&self) -> i32 {
        self.bit_size
    }

    pub fn bit_size_exact(&self) -> f64 {
        self.bit_size_exact
    }

    /// Map a letter to its reduced value.
    #[inline]
    pub fn reduce(&self, a: Letter) -> u32 {
        self.map[a as u8 as usize]
    }

    #[inline(always)]
    pub(crate) fn reduce_pair10(&self, a: Letter, b: Letter) -> u16 {
        // Callers strip masking bits first, so both values are in 0..32.
        unsafe { *self.pair10.get_unchecked(((a as usize) << 5) | b as usize) }
    }

    /// Get the map8 lookup table (for SIMD use).
    pub fn map8(&self) -> &[Letter; 256] {
        &self.map8
    }

    /// Get the map8b lookup table (for SIMD use).
    pub fn map8b(&self) -> &[Letter; 256] {
        &self.map8b
    }

    #[inline(always)]
    pub(crate) fn reduce_seed_letter(&self, letter: Letter) -> Letter {
        unsafe { *self.seed_map8.get_unchecked(letter as u8 as usize) }
    }

    /// Get the frequency for a reduced bucket.
    pub fn freq(&self, bucket: u32) -> f64 {
        self.freq[bucket as usize]
    }

    /// Reduce an entire sequence.
    pub fn reduce_seq(&self, seq: &[Letter]) -> Vec<Letter> {
        seq.iter().map(|&l| self.reduce(l) as Letter).collect()
    }

    pub fn decode_seed(&self, seed: u64, len: usize) -> String {
        let mut s = vec![b'-'; len];
        let mut c = seed;
        for i in 0..len {
            s[len - i - 1] = AMINO_ACID_ALPHABET[(c % self.size as u64) as usize];
            c /= self.size as u64;
        }
        String::from_utf8(s).unwrap()
    }
}

impl std::fmt::Display for Reduction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.size {
            write!(f, "[")?;
            for j in 0..20 {
                if self.map[j] == i {
                    write!(f, "{}", AMINO_ACID_ALPHABET[j] as char)?;
                }
            }
            write!(f, "]")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_reduction() {
        let r = Reduction::default_reduction();
        // "A KR EDNQ C G H ILVM FYW P ST" = 10 groups
        assert_eq!(r.size(), 10);
        assert_eq!(r.bit_size(), 4);
    }

    #[test]
    fn test_reduction_mapping() {
        let r = Reduction::default_reduction();
        // A=0 maps to group 0
        assert_eq!(r.reduce(0), 0);
        // K=11, R=1 should be in same group
        assert_eq!(r.reduce(11), r.reduce(1));
    }

    #[test]
    fn seed_reduction_table_preserves_masking_semantics() {
        let reduction = Reduction::default_reduction();
        for byte in 0u16..=255 {
            let raw = byte as u8 as Letter;
            let expected = if byte & 0x80 != 0 {
                MASK_LETTER
            } else {
                let letter = raw & super::super::value::LETTER_MASK;
                if letter < TRUE_AA as Letter {
                    reduction.reduce(letter) as Letter
                } else {
                    MASK_LETTER
                }
            };
            assert_eq!(reduction.reduce_seed_letter(raw), expected);
        }
    }

    #[test]
    fn test_reduction_bucket_frequencies_match_blosum62_groups() {
        let r = Reduction::default_reduction();
        assert!((r.freq(0) - BLOSUM62_BACKGROUND_FREQS[0].ln()).abs() < 1e-15);
        let kr = (BLOSUM62_BACKGROUND_FREQS[11] + BLOSUM62_BACKGROUND_FREQS[1]).ln();
        assert!((r.freq(1) - kr).abs() < 1e-15);
        assert!(r.freq(10).is_infinite() && r.freq(10).is_sign_negative());
    }

    #[test]
    fn test_reduction_display_and_decode_seed() {
        let r = Reduction::default_reduction();
        assert_eq!(format!("{}", r), "[A][RK][NDQE][C][G][H][ILMV][FWY][P][ST]");
        assert_eq!(r.decode_seed(0, 3), "AAA");
        assert_eq!(r.decode_seed(1, 3), "AAR");
        assert_eq!(r.decode_seed(10, 3), "ARA");
    }

    #[test]
    #[should_panic(expected = "Invalid character in sequence: '!'")]
    fn test_reduction_rejects_unknown_definition_letters() {
        let _ = Reduction::new("A !", AMINO_ACID_ALPHABET);
    }
}
