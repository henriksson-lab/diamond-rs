use crate::basic::value::{Letter, Loc, DELIMITER_LETTER};

/// A collection of sequences stored contiguously in memory.
///
/// Sequences are stored end-to-end separated by DELIMITER_LETTER bytes,
/// with an offset array for O(1) random access to any sequence.
pub struct SequenceSet {
    /// Raw data: all sequences concatenated with delimiter separators.
    data: Vec<Letter>,
    /// Offsets into data where each sequence starts.
    /// offsets[i] is the start of sequence i, offsets[i+1]-1 is the end (exclusive of delimiter).
    offsets: Vec<usize>,
}

impl SequenceSet {
    pub fn new() -> Self {
        SequenceSet {
            data: vec![DELIMITER_LETTER], // sentinel at start
            offsets: vec![1], // first sequence starts after sentinel
        }
    }

    /// Add a sequence to the set.
    pub fn push(&mut self, seq: &[Letter]) {
        self.data.extend_from_slice(seq);
        self.data.push(DELIMITER_LETTER);
        self.offsets.push(self.data.len());
    }

    /// Number of sequences.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get sequence by index.
    pub fn get(&self, i: usize) -> &[Letter] {
        let start = self.offsets[i];
        let end = self.offsets[i + 1] - 1; // -1 for delimiter
        &self.data[start..end]
    }

    /// Length of sequence i.
    pub fn seq_length(&self, i: usize) -> Loc {
        (self.offsets[i + 1] - 1 - self.offsets[i]) as Loc
    }

    /// Total letters across all sequences.
    pub fn letters(&self) -> u64 {
        let mut total = 0u64;
        for i in 0..self.len() {
            total += self.seq_length(i) as u64;
        }
        total
    }

    /// Get the raw data (for SIMD access or direct manipulation).
    pub fn data(&self) -> &[Letter] {
        &self.data
    }

    /// Get mutable access to a sequence's data.
    pub fn get_mut(&mut self, i: usize) -> &mut [Letter] {
        let start = self.offsets[i];
        let end = self.offsets[i + 1] - 1;
        &mut self.data[start..end]
    }
}

impl Default for SequenceSet {
    fn default() -> Self {
        Self::new()
    }
}

/// A Block holds a batch of sequences with their identifiers.
///
/// This is the primary in-memory representation for query and reference
/// sequences during search and alignment.
pub struct Block {
    /// Encoded sequences.
    pub seqs: SequenceSet,
    /// Sequence identifiers.
    pub ids: Vec<String>,
    /// Mapping from block IDs to original database OIDs.
    pub block2oid: Vec<u64>,
}

impl Block {
    pub fn new() -> Self {
        Block {
            seqs: SequenceSet::new(),
            ids: Vec::new(),
            block2oid: Vec::new(),
        }
    }

    /// Add a sequence with its ID.
    pub fn push(&mut self, id: &str, seq: &[Letter]) {
        self.seqs.push(seq);
        self.ids.push(id.to_string());
    }

    /// Number of sequences.
    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }
}

impl Default for Block {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_set() {
        let mut ss = SequenceSet::new();
        ss.push(&[0, 1, 2, 3]); // ARND
        ss.push(&[4, 5, 6]);    // CQE

        assert_eq!(ss.len(), 2);
        assert_eq!(ss.get(0), &[0, 1, 2, 3]);
        assert_eq!(ss.get(1), &[4, 5, 6]);
        assert_eq!(ss.seq_length(0), 4);
        assert_eq!(ss.seq_length(1), 3);
        assert_eq!(ss.letters(), 7);
    }

    #[test]
    fn test_sequence_set_empty() {
        let ss = SequenceSet::new();
        assert_eq!(ss.len(), 0);
        assert!(ss.is_empty());
    }

    #[test]
    fn test_block() {
        let mut block = Block::new();
        block.push("seq1", &[0, 1, 2]);
        block.push("seq2", &[3, 4, 5, 6]);

        assert_eq!(block.len(), 2);
        assert_eq!(block.ids[0], "seq1");
        assert_eq!(block.ids[1], "seq2");
        assert_eq!(block.seqs.get(0), &[0, 1, 2]);
        assert_eq!(block.seqs.get(1), &[3, 4, 5, 6]);
    }

    #[test]
    fn test_sequence_set_mutate() {
        let mut ss = SequenceSet::new();
        ss.push(&[0, 1, 2]);
        ss.get_mut(0)[1] = 10;
        assert_eq!(ss.get(0), &[0, 10, 2]);
    }
}
