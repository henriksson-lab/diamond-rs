use super::value::{Letter, Loc, DELIMITER_LETTER, LETTER_MASK, MASK_LETTER};

/// A lightweight view over a letter sequence (analogous to C++ `Sequence`).
/// Does not own its data — it's a borrowed slice reference.
#[derive(Clone, Copy, Debug)]
pub struct Sequence<'a> {
    data: &'a [Letter],
}

impl<'a> Sequence<'a> {
    pub const DELIMITER: Letter = DELIMITER_LETTER;

    pub fn new(data: &'a [Letter]) -> Self {
        Sequence { data }
    }

    pub fn empty() -> Self {
        Sequence { data: &[] }
    }

    #[inline]
    pub fn length(&self) -> Loc {
        self.data.len() as Loc
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    #[inline]
    pub fn data(&self) -> &'a [Letter] {
        self.data
    }

    /// Access a letter with SEQ_MASK applied (matching C++ operator[]).
    #[inline]
    pub fn get(&self, i: usize) -> Letter {
        self.data[i] & LETTER_MASK
    }

    /// Get a subsequence [begin..end).
    pub fn subseq(&self, begin: usize, end: usize) -> Sequence<'a> {
        Sequence {
            data: &self.data[begin..end],
        }
    }

    /// Get a subsequence from begin to end of sequence.
    pub fn subseq_from(&self, begin: usize) -> Sequence<'a> {
        Sequence {
            data: &self.data[begin..],
        }
    }

    /// Convert to a string representation using an alphabet.
    pub fn to_string(&self, alphabet: &[u8]) -> String {
        self.data
            .iter()
            .map(|&l| alphabet[(l & LETTER_MASK) as usize] as char)
            .collect()
    }

    /// Copy the sequence data into a new Vec.
    pub fn to_vec(&self) -> Vec<Letter> {
        self.data.to_vec()
    }

    /// Return a reversed copy.
    pub fn reverse(&self) -> Vec<Letter> {
        let mut v = self.data.to_vec();
        v.reverse();
        v
    }

    /// Count masked letters.
    pub fn masked_letters(&self) -> Loc {
        self.data
            .iter()
            .filter(|&&l| (l & LETTER_MASK) == MASK_LETTER)
            .count() as Loc
    }

    /// Ratio of masked letters.
    pub fn masked_letter_ratio(&self) -> f64 {
        self.masked_letters() as f64 / self.data.len() as f64
    }

    /// Length ratio (smaller/larger) between two sequences.
    pub fn length_ratio(&self, other: &Sequence) -> f64 {
        let a = self.length();
        let b = other.length();
        if a < b {
            a as f64 / b as f64
        } else {
            b as f64 / a as f64
        }
    }

    /// Parse a string into a sequence using the given alphabet traits.
    pub fn from_string(s: &str, alphabet: &[u8], mask_char: Letter) -> Vec<Letter> {
        let mut result = Vec::with_capacity(s.len());
        let mut lookup = [0i8; 256];
        for (i, &ch) in alphabet.iter().enumerate() {
            lookup[ch as usize] = i as i8;
            lookup[(ch as char).to_ascii_lowercase() as usize] = i as i8;
        }
        for ch in s.bytes() {
            let val = lookup[ch as usize];
            if val == 0 && ch != alphabet[0] && (ch as char).to_ascii_lowercase() != alphabet[0] as char {
                result.push(mask_char);
            } else {
                result.push(val);
            }
        }
        result
    }
}

impl<'a> PartialEq for Sequence<'a> {
    fn eq(&self, other: &Self) -> bool {
        if self.data.len() != other.data.len() {
            return false;
        }
        self.data
            .iter()
            .zip(other.data.iter())
            .all(|(&a, &b)| (a & LETTER_MASK) == (b & LETTER_MASK))
    }
}

impl<'a> Eq for Sequence<'a> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::AMINO_ACID_ALPHABET;

    #[test]
    fn test_sequence_basic() {
        let data = vec![0i8, 1, 2, 3, 4];
        let seq = Sequence::new(&data);
        assert_eq!(seq.length(), 5);
        assert_eq!(seq.get(0), 0);
        assert_eq!(seq.get(4), 4);
    }

    #[test]
    fn test_sequence_subseq() {
        let data = vec![0i8, 1, 2, 3, 4, 5];
        let seq = Sequence::new(&data);
        let sub = seq.subseq(1, 4);
        assert_eq!(sub.length(), 3);
        assert_eq!(sub.get(0), 1);
        assert_eq!(sub.get(2), 3);
    }

    #[test]
    fn test_sequence_to_string() {
        let data = vec![0i8, 1, 2, 3];
        let seq = Sequence::new(&data);
        let s = seq.to_string(AMINO_ACID_ALPHABET);
        assert_eq!(s, "ARND");
    }

    #[test]
    fn test_sequence_equality() {
        let a = vec![0i8, 1, 2];
        let b = vec![0i8, 1, 2];
        let c = vec![0i8, 1, 3];
        assert_eq!(Sequence::new(&a), Sequence::new(&b));
        assert_ne!(Sequence::new(&a), Sequence::new(&c));
    }

    #[test]
    fn test_sequence_masked() {
        let data = vec![0i8, MASK_LETTER, 2, MASK_LETTER];
        let seq = Sequence::new(&data);
        assert_eq!(seq.masked_letters(), 2);
        assert!((seq.masked_letter_ratio() - 0.5).abs() < 0.001);
    }
}
