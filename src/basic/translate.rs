/// DNA strand direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(u8)]
pub enum Strand {
    #[default]
    Forward = 0,
    Reverse = 1,
}

/// Reading frame for translated sequences.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
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

/// Position in a translated query context.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TranslatedPosition {
    pub frame: Frame,
    pub translated: i32,
}

impl TranslatedPosition {
    pub fn new(translated: i32, frame: Frame) -> Self {
        TranslatedPosition { frame, translated }
    }

    pub fn from_in_strand(in_strand: i32, strand: Strand, query_translated: bool) -> Self {
        TranslatedPosition {
            frame: Frame::new(strand, in_strand % 3),
            translated: Self::in_strand_to_translated(in_strand, query_translated),
        }
    }

    pub fn shift_forward(&mut self) {
        self.frame.offset += 1;
        if self.frame.offset == 3 {
            self.frame.offset = 0;
            self.translated += 1;
        }
    }

    pub fn shift_back(&mut self) {
        self.frame.offset -= 1;
        if self.frame.offset == -1 {
            self.frame.offset = 2;
            self.translated -= 1;
        }
    }

    pub fn shift_forward_by(&mut self, mut k: i32) {
        while k > 0 {
            self.shift_forward();
            k -= 1;
        }
    }

    pub fn frame_shift(&self, x: &TranslatedPosition) -> i32 {
        const FRAME_SHIFT: [[i32; 3]; 3] = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]];
        FRAME_SHIFT[self.frame.offset as usize][x.frame.offset as usize]
    }

    pub fn absolute(&self, dna_len: i32, query_translated: bool, blastn: bool) -> i32 {
        if self.frame.offset == 0 && blastn {
            return dna_len - 1 - self.translated;
        }
        if !query_translated && self.frame.strand == Strand::Forward {
            return self.translated;
        }
        Self::oriented_position(self.in_strand(query_translated), self.frame.strand, dna_len)
    }

    pub fn absolute_interval(
        begin: TranslatedPosition,
        end: TranslatedPosition,
        dna_len: i32,
        query_translated: bool,
    ) -> crate::util::interval::Interval {
        if begin.frame.strand == Strand::Forward {
            crate::util::interval::Interval::new(
                begin.in_strand(query_translated),
                end.in_strand(query_translated),
            )
        } else {
            crate::util::interval::Interval::new(
                Self::oriented_position(
                    end.in_strand(query_translated) - 1,
                    Strand::Reverse,
                    dna_len,
                ),
                Self::oriented_position(
                    begin.in_strand(query_translated) - 1,
                    Strand::Reverse,
                    dna_len,
                ),
            )
        }
    }

    pub fn in_strand_to_translated(in_strand: i32, query_translated: bool) -> i32 {
        if query_translated {
            in_strand / 3
        } else {
            in_strand
        }
    }

    pub fn translated_to_in_strand(translated: i32, frame: Frame, query_translated: bool) -> i32 {
        if query_translated {
            frame.offset + 3 * translated
        } else {
            translated
        }
    }

    pub fn in_strand(&self, query_translated: bool) -> i32 {
        Self::translated_to_in_strand(self.translated, self.frame, query_translated)
    }

    pub fn oriented_position(pos: i32, strand: Strand, dna_len: i32) -> i32 {
        if strand == Strand::Forward {
            pos
        } else {
            dna_len - pos - 1
        }
    }

    pub fn absolute_to_translated(src: i32, frame: Frame, dna_len: i32, translated: bool) -> i32 {
        if !translated {
            src
        } else {
            Self::in_strand_to_translated(Self::oriented_position(src, frame.strand, dna_len), true)
        }
    }

    pub fn translated_to_absolute(translated: i32, frame: Frame, dna_len: i32) -> i32 {
        Self::oriented_position(
            Self::translated_to_in_strand(translated, frame, true),
            frame.strand,
            dna_len,
        )
    }
}

impl std::ops::Add<i32> for TranslatedPosition {
    type Output = TranslatedPosition;

    fn add(self, rhs: i32) -> Self::Output {
        TranslatedPosition::new(self.translated + rhs, self.frame)
    }
}

impl std::ops::Sub<i32> for TranslatedPosition {
    type Output = TranslatedPosition;

    fn sub(self, rhs: i32) -> Self::Output {
        TranslatedPosition::new(self.translated - rhs, self.frame)
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
    11, 2, 11, 2, 16, 16, 16, 16, 1, 15, 1, 15, 9, 9, 12, 9, // A**
    5, 8, 5, 8, 14, 14, 14, 14, 1, 1, 1, 1, 10, 10, 10, 10, // C**
    6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7, 19, 19, 19, 19, // G**
    24, 18, 24, 18, 15, 15, 15, 15, 24, 4, 17, 4, 10, 13, 10, 13, // T**
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

    #[test]
    fn test_translated_position_forward_absolute_interval() {
        let frame = Frame::new(Strand::Forward, 1);
        let begin = TranslatedPosition::new(2, frame);
        let end = TranslatedPosition::new(5, frame);

        assert_eq!(begin.in_strand(true), 7);
        assert_eq!(begin.absolute(30, true, false), 7);
        assert_eq!(
            TranslatedPosition::absolute_interval(begin, end, 30, true),
            crate::util::interval::Interval::new(7, 16)
        );
    }

    #[test]
    fn test_translated_position_reverse_absolute_interval() {
        let frame = Frame::new(Strand::Reverse, 2);
        let begin = TranslatedPosition::new(1, frame);
        let end = TranslatedPosition::new(4, frame);

        assert_eq!(begin.in_strand(true), 5);
        assert_eq!(begin.absolute(30, true, false), 24);
        assert_eq!(
            TranslatedPosition::absolute_interval(begin, end, 30, true),
            crate::util::interval::Interval::new(16, 25)
        );
    }

    #[test]
    fn test_translated_position_shift_and_conversion() {
        let mut pos = TranslatedPosition::new(10, Frame::new(Strand::Forward, 1));
        pos.shift_forward();
        assert_eq!(
            pos,
            TranslatedPosition::new(10, Frame::new(Strand::Forward, 2))
        );
        pos.shift_forward();
        assert_eq!(
            pos,
            TranslatedPosition::new(11, Frame::new(Strand::Forward, 0))
        );
        pos.shift_back();
        assert_eq!(
            pos,
            TranslatedPosition::new(10, Frame::new(Strand::Forward, 2))
        );

        assert_eq!(
            TranslatedPosition::absolute_to_translated(
                24,
                Frame::new(Strand::Reverse, 2),
                30,
                true
            ),
            1
        );
        assert_eq!(
            TranslatedPosition::translated_to_absolute(1, Frame::new(Strand::Reverse, 2), 30),
            24
        );
    }
}
