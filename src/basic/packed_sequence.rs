use super::value::{Letter, SequenceType, LETTER_MASK};

/// Check if a nucleotide sequence contains N (value 4).
pub fn has_n(seq: &[Letter]) -> bool {
    seq.iter().any(|&l| (l & LETTER_MASK) == 4)
}

/// Packed sequence representation.
///
/// Amino acids use 5 bits per letter, nucleotides use 2 bits (ACGT)
/// or 3 bits (ACGTN) depending on whether N is present.
pub struct PackedSequence {
    data: Vec<u8>,
    has_n: bool,
}

impl PackedSequence {
    /// Pack a sequence according to its type.
    pub fn new(seq: &[Letter], seq_type: SequenceType) -> Self {
        let has_n = match seq_type {
            SequenceType::Nucleotide => has_n(seq),
            SequenceType::AminoAcid => false,
        };

        let bits_per_letter = match seq_type {
            SequenceType::Nucleotide => {
                if has_n { 3 } else { 2 }
            }
            SequenceType::AminoAcid => 5,
        };

        let data = pack(seq, bits_per_letter);
        PackedSequence { data, has_n }
    }

    /// Create from raw packed data.
    pub fn from_raw(data: Vec<u8>, has_n: bool) -> Self {
        PackedSequence { data, has_n }
    }

    /// Read packed data from a byte slice.
    pub fn from_bytes(bytes: &[u8], len: usize, has_n: bool, bits_per_letter: u32) -> Self {
        let byte_count = (len * bits_per_letter as usize).div_ceil(8);
        let data = bytes[..byte_count].to_vec();
        PackedSequence { data, has_n }
    }

    /// Unpack the sequence back to letter values.
    pub fn unpack(&self, bits_per_letter: u32, len: usize) -> Vec<Letter> {
        let mut dst = Vec::with_capacity(len);
        let mut x: u32 = 0;
        let mut n: u32 = 0;
        let mask: u32 = (1 << bits_per_letter) - 1;
        for &byte in &self.data {
            x |= (byte as u32) << n;
            n += 8;
            while n >= bits_per_letter && dst.len() < len {
                dst.push((x & mask) as Letter);
                n -= bits_per_letter;
                x >>= bits_per_letter;
            }
        }
        dst
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn has_n(&self) -> bool {
        self.has_n
    }
}

fn pack(seq: &[Letter], bits_per_letter: u32) -> Vec<u8> {
    let mut data = Vec::new();
    let mut x: u32 = 0;
    let mut n: u32 = 0;
    for &letter in seq {
        let l = (letter & LETTER_MASK) as u32;
        x |= l << n;
        n += bits_per_letter;
        if n >= 8 {
            data.push((x & 0xff) as u8);
            n -= 8;
            x >>= 8;
        }
    }
    if n > 0 {
        data.push((x & 0xff) as u8);
    }
    data
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_unpack_amino_acid() {
        let seq: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let packed = PackedSequence::new(&seq, SequenceType::AminoAcid);
        let unpacked = packed.unpack(5, seq.len());
        assert_eq!(unpacked, seq);
    }

    #[test]
    fn test_pack_unpack_nucleotide_no_n() {
        let seq: Vec<Letter> = vec![0, 1, 2, 3, 0, 1]; // ACGTAC
        let packed = PackedSequence::new(&seq, SequenceType::Nucleotide);
        assert!(!packed.has_n());
        let unpacked = packed.unpack(2, seq.len());
        assert_eq!(unpacked, seq);
    }

    #[test]
    fn test_pack_unpack_nucleotide_with_n() {
        let seq: Vec<Letter> = vec![0, 1, 2, 3, 4, 0]; // ACGTN A
        let packed = PackedSequence::new(&seq, SequenceType::Nucleotide);
        assert!(packed.has_n());
        let unpacked = packed.unpack(3, seq.len());
        assert_eq!(unpacked, seq);
    }

    #[test]
    fn test_roundtrip_various_lengths() {
        for len in 1..=20 {
            let seq: Vec<Letter> = (0..len).map(|i| (i % 20) as Letter).collect();
            let packed = PackedSequence::new(&seq, SequenceType::AminoAcid);
            let unpacked = packed.unpack(5, seq.len());
            assert_eq!(unpacked, seq, "Failed for length {}", len);
        }
    }
}
