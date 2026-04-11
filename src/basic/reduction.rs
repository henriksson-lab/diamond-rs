use super::value::{
    Letter, AMINO_ACID_ALPHABET, DELIMITER_LETTER, MASK_LETTER, STOP_LETTER, TRUE_AA,
};

/// Default reduction definition: maps 20 amino acids to 8 groups.
pub const DEFAULT_REDUCTION: &str = "A KR EDNQ C G H ILVM FYW P ST";

/// Alphabet reduction mapping amino acids to fewer groups for seed computation.
pub struct Reduction {
    map: [u32; 256],
    map8: [Letter; 256],
    map8b: [Letter; 256],
    size: u32,
    bit_size: i32,
    bit_size_exact: f64,
    freq: [f64; TRUE_AA as usize],
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

        let freq = [0.0f64; TRUE_AA as usize];

        // Build a simple char->letter lookup from the alphabet
        let mut char_to_letter = [0u8; 256];
        for (i, &ch) in alphabet.iter().enumerate() {
            char_to_letter[ch as usize] = i as u8;
            char_to_letter[(ch as char).to_ascii_lowercase() as usize] = i as u8;
        }

        for (i, token) in tokens.iter().enumerate() {
            for ch in token.bytes() {
                let letter = char_to_letter[ch as usize] as usize;
                map[letter] = i as u32;
                map8[letter] = i as Letter;
                map8b[letter] = i as Letter;
            }
        }

        map8[MASK_LETTER as u8 as usize] = size as Letter;
        map8[STOP_LETTER as u8 as usize] = size as Letter;
        map8[DELIMITER_LETTER as u8 as usize] = size as Letter;
        map8b[MASK_LETTER as u8 as usize] = (size + 1) as Letter;
        map8b[STOP_LETTER as u8 as usize] = (size + 1) as Letter;
        map8b[DELIMITER_LETTER as u8 as usize] = (size + 1) as Letter;

        Reduction {
            map,
            map8,
            map8b,
            size,
            bit_size,
            bit_size_exact,
            freq,
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

    /// Get the map8 lookup table (for SIMD use).
    pub fn map8(&self) -> &[Letter; 256] {
        &self.map8
    }

    /// Get the map8b lookup table (for SIMD use).
    pub fn map8b(&self) -> &[Letter; 256] {
        &self.map8b
    }

    /// Get the frequency for a reduced bucket.
    pub fn freq(&self, bucket: u32) -> f64 {
        self.freq[bucket as usize]
    }

    /// Reduce an entire sequence.
    pub fn reduce_seq(&self, seq: &[Letter]) -> Vec<Letter> {
        seq.iter().map(|&l| self.reduce(l) as Letter).collect()
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
}
