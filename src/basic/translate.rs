/// DNA strand direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(u8)]
pub enum Strand {
    #[default]
    Forward = 0,
    Reverse = 1,
}

/// Reading frame for translated sequences.
#[derive(Debug, Clone, Copy, Default)]
pub struct Frame {
    pub offset: i32,
    pub strand: Strand,
}


impl Frame {
    pub fn new(strand: Strand, offset: i32) -> Self {
        Frame { offset, strand }
    }

    /// Create from a 0-5 index (0-2 = forward frames, 3-5 = reverse frames).
    pub fn from_index(index: i32) -> Self {
        Frame {
            offset: index % 3,
            strand: if index < 3 {
                Strand::Forward
            } else {
                Strand::Reverse
            },
        }
    }

    /// Convert to a 0-5 index.
    pub fn index(&self) -> i32 {
        self.strand as i32 * 3 + self.offset
    }

    /// The signed frame number (+1/+2/+3 or -1/-2/-3).
    pub fn signed_frame(&self) -> i32 {
        let sign = if self.strand == Strand::Forward {
            1
        } else {
            -1
        };
        (self.offset + 1) * sign
    }

    /// Length of this frame's translation given DNA length.
    pub fn length(&self, dna_len: i32) -> i32 {
        ((dna_len - self.offset) / 3).max(0)
    }
}

/// Standard genetic code (NCBI translation table 1).
/// Maps codon triplets (4^3 = 64 entries) to amino acid letters.
static STANDARD_GENETIC_CODE: [u8; 64] = [
    // TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG
    // CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG
    // ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG
    // GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG
    // K  N  K  N  T  T  T  T  R  S  R  S  I  I  M  I
    // Q  H  Q  H  P  P  P  P  R  R  R  R  L  L  L  L
    // E  D  E  D  A  A  A  A  G  G  G  G  V  V  V  V
    // *  Y  *  Y  S  S  S  S  *  C  W  C  L  F  L  F
    // Using diamond's letter encoding: A=0 R=1 N=2 D=3 C=4 Q=5 E=6 G=7 H=8 I=9
    // L=10 K=11 M=12 F=13 P=14 S=15 T=16 W=17 Y=18 V=19

    // Codon index = base1*16 + base2*4 + base3, where A=0 C=1 G=2 T=3
    11, 2, 11, 2,   16, 16, 16, 16,   1, 15, 1, 15,   9, 9, 12, 9,   // A**
    5, 8, 5, 8,     14, 14, 14, 14,   1, 1, 1, 1,     10, 10, 10, 10, // C**
    6, 3, 6, 3,     0, 0, 0, 0,       7, 7, 7, 7,     19, 19, 19, 19, // G**
    24, 18, 24, 18, 15, 15, 15, 15,   24, 4, 17, 4,   10, 13, 10, 13, // T**
];

/// Translate a DNA codon (3 nucleotide letters) to an amino acid letter.
///
/// Nucleotide encoding: A=0, C=1, G=2, T=3, N=4
/// Returns the amino acid letter, or MASK_LETTER (23) for codons with N.
pub fn translate_codon(a: u8, b: u8, c: u8) -> u8 {
    if a >= 4 || b >= 4 || c >= 4 {
        return 23; // MASK_LETTER for ambiguous
    }
    STANDARD_GENETIC_CODE[(a as usize) * 16 + (b as usize) * 4 + c as usize]
}

/// Translate a DNA sequence into all 6 reading frames.
///
/// Returns 6 amino acid sequences: frames 0-2 are forward, 3-5 are reverse.
pub fn translate_6_frames(dna: &[u8]) -> [Vec<u8>; 6] {
    let len = dna.len();
    let mut frames = [
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ];

    // Reverse complement lookup: A<->T, C<->G
    let complement = |b: u8| -> u8 {
        match b {
            0 => 3, // A -> T
            1 => 2, // C -> G
            2 => 1, // G -> C
            3 => 0, // T -> A
            _ => 4, // N -> N
        }
    };

    // Forward frames
    for (offset, frame) in frames[..3].iter_mut().enumerate() {
        let mut i = offset;
        while i + 2 < len {
            frame.push(translate_codon(dna[i], dna[i + 1], dna[i + 2]));
            i += 3;
        }
    }

    // Reverse complement frames
    for (offset, frame) in frames[3..].iter_mut().enumerate() {
        let mut i = len as i64 - 1 - offset as i64;
        while i >= 2 {
            frame.push(translate_codon(
                complement(dna[i as usize]),
                complement(dna[i as usize - 1]),
                complement(dna[i as usize - 2]),
            ));
            i -= 3;
        }
    }

    frames
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_frame_from_index() {
        let f = Frame::from_index(0);
        assert_eq!(f.strand, Strand::Forward);
        assert_eq!(f.offset, 0);

        let f = Frame::from_index(4);
        assert_eq!(f.strand, Strand::Reverse);
        assert_eq!(f.offset, 1);
    }

    #[test]
    fn test_frame_signed() {
        assert_eq!(Frame::from_index(0).signed_frame(), 1);
        assert_eq!(Frame::from_index(1).signed_frame(), 2);
        assert_eq!(Frame::from_index(3).signed_frame(), -1);
    }

    #[test]
    fn test_translate_codon_met() {
        // ATG = methionine (M=12)
        assert_eq!(translate_codon(0, 3, 2), 12);
    }

    #[test]
    fn test_translate_codon_stop() {
        // TAA = stop (*=24)
        assert_eq!(translate_codon(3, 0, 0), 24);
    }

    #[test]
    fn test_translate_codon_ambiguous() {
        // Codons with N should give mask letter
        assert_eq!(translate_codon(4, 0, 0), 23);
    }

    #[test]
    fn test_translate_6_frames() {
        // Simple ATG sequence
        let dna = vec![0u8, 3, 2]; // ATG
        let frames = translate_6_frames(&dna);
        assert_eq!(frames[0].len(), 1);
        assert_eq!(frames[0][0], 12); // M
    }

    #[test]
    fn test_strand_default() {
        assert_eq!(Strand::default(), Strand::Forward);
    }
}
