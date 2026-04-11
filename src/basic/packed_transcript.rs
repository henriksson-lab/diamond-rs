use super::value::{Letter, AMINO_ACID_COUNT};

/// Edit operation types in alignment transcripts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum EditOperation {
    Match = 0,
    Insertion = 1,
    Deletion = 2,
    Substitution = 3,
    FrameshiftForward = 4,
    FrameshiftReverse = 5,
}

/// A single packed operation byte.
///
/// Bits [7-6]: operation type (match, insertion, deletion, substitution)
/// Bits [5-0]: count or letter value
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct PackedOperation {
    pub code: u8,
}

impl PackedOperation {
    const OP_BITS: u32 = 2;
    const COUNT_BITS: u32 = 8 - Self::OP_BITS;
    pub const MAX_COUNT: u32 = (1 << Self::COUNT_BITS) - 1;

    pub fn new(code: u8) -> Self {
        PackedOperation { code }
    }

    pub fn from_op_count(op: EditOperation, count: u32) -> Self {
        PackedOperation {
            code: ((op as u8) << Self::COUNT_BITS) | (count as u8),
        }
    }

    pub fn from_op_letter(op: EditOperation, letter: Letter) -> Self {
        PackedOperation {
            code: ((op as u8) << Self::COUNT_BITS) | (letter as u8),
        }
    }

    /// Get the edit operation, handling frameshift encoding.
    pub fn op(&self) -> EditOperation {
        let o = self.code >> Self::COUNT_BITS;
        match o {
            3 => {
                // Substitution - check for frameshift encoding
                let l = self.letter();
                if l == AMINO_ACID_COUNT as Letter {
                    EditOperation::FrameshiftReverse
                } else if l == (AMINO_ACID_COUNT + 1) as Letter {
                    EditOperation::FrameshiftForward
                } else {
                    EditOperation::Substitution
                }
            }
            0 => EditOperation::Match,
            1 => EditOperation::Insertion,
            2 => EditOperation::Deletion,
            _ => unreachable!(),
        }
    }

    /// Get the count (for match/insertion) or 1 (for deletion/substitution).
    pub fn count(&self) -> u32 {
        match self.op() {
            EditOperation::Match | EditOperation::Insertion => {
                (self.code & Self::MAX_COUNT as u8) as u32
            }
            _ => 1,
        }
    }

    /// Get the letter value (for substitution/deletion).
    pub fn letter(&self) -> Letter {
        (self.code & Self::MAX_COUNT as u8) as Letter
    }

    /// Create a terminator operation (match with count 0).
    pub fn terminator() -> Self {
        Self::from_op_count(EditOperation::Match, 0)
    }

    /// Create a forward frameshift operation.
    pub fn frameshift_forward() -> Self {
        Self::from_op_letter(EditOperation::Substitution, (AMINO_ACID_COUNT + 1) as Letter)
    }

    /// Create a reverse frameshift operation.
    pub fn frameshift_reverse() -> Self {
        Self::from_op_letter(EditOperation::Substitution, AMINO_ACID_COUNT as Letter)
    }

    pub fn is_terminator(&self) -> bool {
        self.code == 0
    }
}

/// Combined operation with accumulated count.
#[derive(Debug, Clone)]
pub struct CombinedOperation {
    pub op: EditOperation,
    pub count: u32,
    pub letter: Letter,
}

/// Packed alignment transcript — a sequence of edit operations.
#[derive(Default, Debug, Clone)]
pub struct PackedTranscript {
    data: Vec<PackedOperation>,
}

impl PackedTranscript {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(n: usize) -> Self {
        PackedTranscript {
            data: Vec::with_capacity(n),
        }
    }

    /// Append a single edit operation, merging with previous if possible.
    pub fn push(&mut self, op: EditOperation) {
        match op {
            EditOperation::FrameshiftForward => {
                self.data.push(PackedOperation::frameshift_forward());
            }
            EditOperation::FrameshiftReverse => {
                self.data.push(PackedOperation::frameshift_reverse());
            }
            _ => {
                if self.data.is_empty()
                    || self.data.last().unwrap().op() != op
                    || self.data.last().unwrap().count() == PackedOperation::MAX_COUNT
                {
                    self.data.push(PackedOperation::from_op_count(op, 1));
                } else {
                    self.data.last_mut().unwrap().code += 1;
                }
            }
        }
    }

    /// Append a deletion or substitution with a specific letter.
    pub fn push_with_letter(&mut self, op: EditOperation, letter: Letter) {
        self.data.push(PackedOperation::from_op_letter(op, letter));
    }

    /// Append an operation with a specific count.
    pub fn push_with_count(&mut self, op: EditOperation, mut count: u32) {
        while count > 0 {
            let n = count.min(PackedOperation::MAX_COUNT);
            self.data.push(PackedOperation::from_op_count(op, n));
            count -= n;
        }
    }

    /// Reverse the transcript.
    pub fn reverse(&mut self) {
        self.data.reverse();
    }

    /// Add a terminator.
    pub fn push_terminator(&mut self) {
        self.data.push(PackedOperation::terminator());
    }

    pub fn clear(&mut self) {
        self.data.clear();
    }

    pub fn raw_length(&self) -> usize {
        self.data.len()
    }

    pub fn data(&self) -> &[PackedOperation] {
        &self.data
    }

    /// Iterate over combined operations (consecutive same-type ops are merged).
    pub fn iter(&self) -> TranscriptIterator<'_> {
        TranscriptIterator::new(&self.data)
    }

    /// Read a transcript from raw bytes.
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let data: Vec<PackedOperation> = bytes.iter().map(|&b| PackedOperation::new(b)).collect();
        PackedTranscript { data }
    }
}

/// Iterator that yields combined operations from a packed transcript.
pub struct TranscriptIterator<'a> {
    ops: &'a [PackedOperation],
    pos: usize,
}

impl<'a> TranscriptIterator<'a> {
    fn new(ops: &'a [PackedOperation]) -> Self {
        TranscriptIterator { ops, pos: 0 }
    }
}

impl<'a> Iterator for TranscriptIterator<'a> {
    type Item = CombinedOperation;

    fn next(&mut self) -> Option<CombinedOperation> {
        if self.pos >= self.ops.len() || self.ops[self.pos].is_terminator() {
            return None;
        }

        let op = self.ops[self.pos].op();
        match op {
            EditOperation::Deletion
            | EditOperation::Substitution
            | EditOperation::FrameshiftForward
            | EditOperation::FrameshiftReverse => {
                let letter = self.ops[self.pos].letter();
                self.pos += 1;
                Some(CombinedOperation {
                    op,
                    count: 1,
                    letter,
                })
            }
            EditOperation::Match | EditOperation::Insertion => {
                let mut count = 0u32;
                while self.pos < self.ops.len()
                    && !self.ops[self.pos].is_terminator()
                    && self.ops[self.pos].op() == op
                {
                    count += self.ops[self.pos].count();
                    self.pos += 1;
                }
                Some(CombinedOperation {
                    op,
                    count,
                    letter: 0,
                })
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_packed_operation_match() {
        let op = PackedOperation::from_op_count(EditOperation::Match, 5);
        assert_eq!(op.op(), EditOperation::Match);
        assert_eq!(op.count(), 5);
    }

    #[test]
    fn test_packed_operation_substitution() {
        let op = PackedOperation::from_op_letter(EditOperation::Substitution, 3);
        assert_eq!(op.op(), EditOperation::Substitution);
        assert_eq!(op.letter(), 3);
        assert_eq!(op.count(), 1);
    }

    #[test]
    fn test_packed_operation_frameshift() {
        let fwd = PackedOperation::frameshift_forward();
        assert_eq!(fwd.op(), EditOperation::FrameshiftForward);
        let rev = PackedOperation::frameshift_reverse();
        assert_eq!(rev.op(), EditOperation::FrameshiftReverse);
    }

    #[test]
    fn test_terminator() {
        let t = PackedOperation::terminator();
        assert!(t.is_terminator());
        assert_eq!(t.code, 0);
    }

    #[test]
    fn test_transcript_merge() {
        let mut t = PackedTranscript::new();
        t.push(EditOperation::Match);
        t.push(EditOperation::Match);
        t.push(EditOperation::Match);
        t.push_terminator();
        // Should merge into one operation with count 3
        let ops: Vec<_> = t.iter().collect();
        assert_eq!(ops.len(), 1);
        assert_eq!(ops[0].op, EditOperation::Match);
        assert_eq!(ops[0].count, 3);
    }

    #[test]
    fn test_transcript_mixed() {
        let mut t = PackedTranscript::new();
        t.push_with_count(EditOperation::Match, 5);
        t.push_with_letter(EditOperation::Substitution, 3);
        t.push_with_count(EditOperation::Match, 3);
        t.push_terminator();

        let ops: Vec<_> = t.iter().collect();
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0].count, 5);
        assert_eq!(ops[1].op, EditOperation::Substitution);
        assert_eq!(ops[2].count, 3);
    }
}
