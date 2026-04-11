/// Letter is the fundamental sequence element type, matching the C++ `signed char` (`Letter`).
pub type Letter = i8;

/// Sequence type discriminator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SequenceType {
    AminoAcid = 0,
    Nucleotide = 1,
}

// Type aliases matching the C++ codebase
pub type Loc = i32;
pub type BlockId = u32;
pub type OId = u64;
pub type DictId = i64;
pub type Score = i32;
pub type TaxId = i32;
pub type CentroidId = OId;
pub type SuperBlockId = u32;

// Alphabet strings
pub const AMINO_ACID_ALPHABET: &[u8] = b"ARNDCQEGHILKMFPSTWYVBJZX*_";
pub const AMINO_ACID_COUNT: usize = AMINO_ACID_ALPHABET.len();
pub const NUCLEOTIDE_ALPHABET: &[u8] = b"ACGTN";
pub const NUCLEOTIDE_COUNT: usize = NUCLEOTIDE_ALPHABET.len();

// Special letter values
pub const MASK_LETTER: Letter = 23;
pub const STOP_LETTER: Letter = 24;
pub const SUPER_HARD_MASK: Letter = 25;
pub const DELIMITER_LETTER: Letter = 31;
pub const LETTER_MASK: Letter = 31;
pub const SEED_MASK: Letter = -128;
pub const TRUE_AA: i32 = 20;

#[inline]
pub fn is_amino_acid(x: Letter) -> bool {
    x != MASK_LETTER && x != DELIMITER_LETTER && x != STOP_LETTER
}

#[inline]
pub fn letter_mask(x: Letter) -> Letter {
    // SEQ_MASK is always defined in the default build
    x & LETTER_MASK
}

/// IUPACAA to standard amino acid encoding conversion table.
pub const IUPACAA_TO_STD: [Letter; 32] = [
    -1, 0, 20, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2, MASK_LETTER,
    14, 5, 1, 15, 16, MASK_LETTER, 19, 17, 23, 18, 22, -1, -1, -1, -1, 24,
];

/// NCBI to standard amino acid encoding conversion table.
pub const NCBI_TO_STD: [Letter; 28] = [
    MASK_LETTER, 0, 20, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5,
    1, 15, 16, 19, 17, 23, 18, 22, MASK_LETTER, 24, MASK_LETTER, 21,
];

/// Converts characters to Letter values for a given alphabet.
pub struct CharRepresentation {
    data: [Letter; 256],
}

const INVALID_LETTER: Letter = -1; // 0xff as i8

impl CharRepresentation {
    pub fn new(alphabet: &[u8], mask: Letter, mask_chars: &[u8]) -> Self {
        let mut data = [INVALID_LETTER; 256];
        for (i, &ch) in alphabet.iter().enumerate() {
            data[ch as usize] = i as Letter;
            data[(ch as char).to_ascii_lowercase() as usize] = i as Letter;
        }
        for &ch in mask_chars {
            data[ch as usize] = mask;
            data[(ch as char).to_ascii_lowercase() as usize] = mask;
        }
        CharRepresentation { data }
    }

    pub fn convert(&self, c: u8) -> Result<Letter, String> {
        let val = self.data[c as usize];
        if val == INVALID_LETTER {
            if (32..127).contains(&c) {
                Err(format!("Invalid character in sequence: '{}'", c as char))
            } else {
                Err(format!("Invalid character in sequence: ASCII {}", c))
            }
        } else {
            Ok(val)
        }
    }
}

/// Traits for a particular value/sequence type (amino acid or nucleotide).
pub struct ValueTraits {
    pub alphabet: &'static [u8],
    pub alphabet_size: usize,
    pub mask_char: Letter,
    pub from_char: CharRepresentation,
    pub seq_type: SequenceType,
}

impl ValueTraits {
    pub fn new(alphabet: &'static [u8], mask_char: Letter, ignore: &[u8], seq_type: SequenceType) -> Self {
        ValueTraits {
            alphabet,
            alphabet_size: alphabet.len(),
            mask_char,
            from_char: CharRepresentation::new(alphabet, mask_char, ignore),
            seq_type,
        }
    }

    pub fn to_char(&self, a: Letter) -> char {
        self.alphabet[a as usize] as char
    }
}

/// Alignment mode information.
#[derive(Debug, Clone, Copy)]
pub struct AlignMode {
    pub sequence_type: SequenceType,
    pub input_sequence_type: SequenceType,
    pub mode: i32,
    pub query_contexts: i32,
    pub query_len_factor: i32,
    pub query_translated: bool,
}

impl AlignMode {
    pub const BLASTP: i32 = 2;
    pub const BLASTX: i32 = 3;
    pub const BLASTN: i32 = 4;

    pub fn new(mode: i32) -> Self {
        match mode {
            Self::BLASTP => AlignMode {
                sequence_type: SequenceType::AminoAcid,
                input_sequence_type: SequenceType::AminoAcid,
                mode,
                query_contexts: 1,
                query_len_factor: 1,
                query_translated: false,
            },
            Self::BLASTX => AlignMode {
                sequence_type: SequenceType::AminoAcid,
                input_sequence_type: SequenceType::Nucleotide,
                mode,
                query_contexts: 6,
                query_len_factor: 3,
                query_translated: true,
            },
            Self::BLASTN => AlignMode {
                sequence_type: SequenceType::Nucleotide,
                input_sequence_type: SequenceType::Nucleotide,
                mode,
                query_contexts: 2,
                query_len_factor: 1,
                query_translated: false,
            },
            _ => panic!("Invalid alignment mode: {}", mode),
        }
    }

    pub fn to_string(&self) -> &'static str {
        match self.mode {
            Self::BLASTP => "blastp",
            Self::BLASTX => "blastx",
            Self::BLASTN => "blastn",
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_amino_acid_alphabet() {
        assert_eq!(AMINO_ACID_COUNT, 26);
        assert_eq!(AMINO_ACID_ALPHABET[0], b'A');
        assert_eq!(AMINO_ACID_ALPHABET[1], b'R');
    }

    #[test]
    fn test_nucleotide_alphabet() {
        assert_eq!(NUCLEOTIDE_COUNT, 5);
    }

    #[test]
    fn test_char_representation() {
        let cr = CharRepresentation::new(AMINO_ACID_ALPHABET, MASK_LETTER, b"UO-");
        assert_eq!(cr.convert(b'A').unwrap(), 0);
        assert_eq!(cr.convert(b'a').unwrap(), 0);
        assert_eq!(cr.convert(b'R').unwrap(), 1);
        assert_eq!(cr.convert(b'U').unwrap(), MASK_LETTER);
        assert_eq!(cr.convert(b'O').unwrap(), MASK_LETTER);
        assert!(cr.convert(b'!').is_err());
    }

    #[test]
    fn test_is_amino_acid() {
        assert!(is_amino_acid(0));  // A
        assert!(is_amino_acid(1));  // R
        assert!(!is_amino_acid(MASK_LETTER));
        assert!(!is_amino_acid(STOP_LETTER));
        assert!(!is_amino_acid(DELIMITER_LETTER));
    }

    #[test]
    fn test_letter_mask() {
        assert_eq!(letter_mask(0), 0);
        assert_eq!(letter_mask(SEED_MASK | 5), 5);
    }

    #[test]
    fn test_iupacaa_to_std() {
        // Position 1 => A (0)
        assert_eq!(IUPACAA_TO_STD[1], 0);
        // Position 2 => B (20)
        assert_eq!(IUPACAA_TO_STD[2], 20);
    }

    #[test]
    fn test_align_mode() {
        let mode = AlignMode::new(AlignMode::BLASTP);
        assert_eq!(mode.query_contexts, 1);
        assert!(!mode.query_translated);
        assert_eq!(mode.to_string(), "blastp");

        let mode = AlignMode::new(AlignMode::BLASTX);
        assert_eq!(mode.query_contexts, 6);
        assert!(mode.query_translated);
    }
}
