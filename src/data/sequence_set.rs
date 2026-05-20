use crate::basic::value::{BlockId, Letter, Loc, DELIMITER_LETTER};

pub trait StringSetValue: Copy + Default {
    fn from_i64(value: i64) -> Self;
}

impl StringSetValue for u8 {
    fn from_i64(value: i64) -> Self {
        value as u8
    }
}

impl StringSetValue for i8 {
    fn from_i64(value: i64) -> Self {
        value as i8
    }
}

/// Contiguous padded string storage matching C++ `StringSetBase`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StringSetBase<T, const PADDING_CHAR: i64, const PADDING_LEN: usize = 1>
where
    T: StringSetValue,
{
    data: Vec<T>,
    limits: Vec<i64>,
}

impl<T, const PADDING_CHAR: i64, const PADDING_LEN: usize>
    StringSetBase<T, PADDING_CHAR, PADDING_LEN>
where
    T: StringSetValue,
{
    pub const PERIMETER_PADDING: usize = 256;
    pub const DELIMITER: i64 = PADDING_CHAR;

    pub fn new() -> Self {
        Self {
            data: vec![T::from_i64(PADDING_CHAR); Self::PERIMETER_PADDING],
            limits: vec![Self::PERIMETER_PADDING as i64],
        }
    }

    pub fn finish_reserve(&mut self) {
        let raw_len = self.raw_len() as usize;
        self.data
            .resize(raw_len + Self::PERIMETER_PADDING, T::from_i64(PADDING_CHAR));
    }

    pub fn reserve(&mut self, n: usize) {
        self.limits
            .push(self.raw_len() + n as i64 + PADDING_LEN as i64);
    }

    pub fn reserve_capacity(&mut self, entries: usize, length: usize) {
        self.limits.reserve(entries + 1);
        self.data
            .reserve(length + 2 * Self::PERIMETER_PADDING + entries * PADDING_LEN);
    }

    pub fn clear(&mut self) {
        self.limits.resize(1, Self::PERIMETER_PADDING as i64);
        self.data
            .resize(Self::PERIMETER_PADDING, T::from_i64(PADDING_CHAR));
    }

    pub fn shrink_to_fit(&mut self) {
        self.limits.shrink_to_fit();
        self.data.shrink_to_fit();
    }

    pub fn push_back(&mut self, s: &[T]) {
        self.limits
            .push(self.raw_len() + s.len() as i64 + PADDING_LEN as i64);
        self.data.extend_from_slice(s);
        self.data
            .extend(std::iter::repeat_n(T::from_i64(PADDING_CHAR), PADDING_LEN));
    }

    pub fn append(&mut self, s: &Self) {
        let n = s.size() as usize;
        if n == 0 {
            return;
        }
        let offset = self.raw_len() - s.limits[0];
        for i in 0..n {
            self.limits.push(s.limits[i + 1] + offset);
        }
        self.data
            .truncate(self.data.len() - Self::PERIMETER_PADDING);
        let begin = s.ptr(0);
        let end = s.end(n - 1) + 1 + Self::PERIMETER_PADDING;
        self.data.extend_from_slice(&s.data[begin..end]);
    }

    pub fn assign(&mut self, i: usize, s: &[T]) {
        let begin = self.ptr(i);
        let end = begin + s.len();
        self.data[begin..end].copy_from_slice(s);
        self.data[end..end + PADDING_LEN].fill(T::from_i64(PADDING_CHAR));
    }

    pub fn fill(&mut self, n: usize, v: T) {
        self.limits
            .push(self.raw_len() + n as i64 + PADDING_LEN as i64);
        self.data.extend(std::iter::repeat_n(v, n));
        self.data
            .extend(std::iter::repeat_n(T::from_i64(PADDING_CHAR), PADDING_LEN));
    }

    pub fn ptr(&self, i: usize) -> usize {
        self.limits[i] as usize
    }

    pub fn end(&self, i: usize) -> usize {
        (self.limits[i + 1] - PADDING_LEN as i64) as usize
    }

    pub fn check_idx(&self, i: usize) -> usize {
        if self.limits.len() < i + 2 {
            panic!("Sequence set index out of bounds.");
        }
        i
    }

    pub fn length(&self, i: usize) -> Loc {
        (self.limits[i + 1] - self.limits[i] - PADDING_LEN as i64) as Loc
    }

    pub fn size(&self) -> BlockId {
        (self.limits.len() - 1) as BlockId
    }

    pub fn empty(&self) -> bool {
        self.limits.len() <= 1
    }

    pub fn raw_len(&self) -> i64 {
        *self.limits.last().unwrap()
    }

    pub fn mem_size(&self) -> i64 {
        (self.data.len() * std::mem::size_of::<T>()
            + self.limits.len() * std::mem::size_of::<i64>()) as i64
    }

    pub fn letters(&self) -> i64 {
        self.raw_len() - self.size() as i64 - Self::PERIMETER_PADDING as i64
    }

    pub fn data(&self, p: u64) -> &[T] {
        &self.data[p as usize..]
    }

    pub fn data_mut(&mut self, p: u64) -> &mut [T] {
        &mut self.data[p as usize..]
    }

    pub fn position(&self, i: BlockId, j: Loc) -> i64 {
        self.limits[i as usize] + j as i64
    }

    pub fn local_position(&self, p: i64) -> (BlockId, Loc) {
        let i = self.limits.partition_point(|&x| x <= p) - 1;
        (i as BlockId, (p - self.limits[i]) as Loc)
    }

    pub fn local_position_batch<Cmp>(&self, positions: &[i64], cmp: Cmp) -> Vec<isize>
    where
        Cmp: Fn(&i64, &i64) -> bool + Copy,
    {
        let mut out = Vec::new();
        crate::util::algo::batch_binary_search(positions, &self.limits, &mut out, cmp, 0);
        out
    }

    pub fn get(&self, i: usize) -> &[T] {
        &self.data[self.ptr(i)..self.end(i)]
    }

    pub fn back(&self) -> &[T] {
        self.get(self.limits.len() - 2)
    }

    pub fn limits(&self) -> &[i64] {
        &self.limits
    }

    pub fn subset(&self, indices: &[usize]) -> Self {
        let mut r = Self::new();
        r.limits.reserve(indices.len());
        for &i in indices {
            r.reserve(self.length(i) as usize);
        }
        r.finish_reserve();
        for (n, &i) in indices.iter().enumerate() {
            r.assign(n, self.get(i));
        }
        r
    }
}

impl<T, const PADDING_CHAR: i64, const PADDING_LEN: usize> Default
    for StringSetBase<T, PADDING_CHAR, PADDING_LEN>
where
    T: StringSetValue,
{
    fn default() -> Self {
        Self::new()
    }
}

pub type StringSet = StringSetBase<u8, 0, 1>;
pub type LetterStringSet = StringSetBase<Letter, { DELIMITER_LETTER as i64 }, 1>;

pub fn max_id_len(ids: &StringSet) -> usize {
    let mut max_len = 0usize;
    for i in 0..ids.size() as usize {
        let id = std::str::from_utf8(ids.get(i)).unwrap_or("");
        let len = id
            .find(|c| crate::util::sequence::ID_DELIMITERS.contains(c))
            .unwrap_or(id.len());
        max_len = max_len.max(len);
    }
    max_len
}

/// A collection of sequences stored contiguously in memory.
///
/// Sequences are stored end-to-end separated by DELIMITER_LETTER bytes,
/// with an offset array for O(1) random access to any sequence.
pub struct SequenceSet {
    /// Raw data: all sequences concatenated with delimiter separators.
    data: Vec<Letter>,
    /// Offsets into data where each sequence starts.
    /// offsets[i] is the start of sequence i, offsets[i+1]-1 is the end (exclusive of delimiter).
    ///
    /// NOTE: C++ `StringSetBase<Letter, DELIMITER, 1>` carries a 256-byte
    /// PERIMETER_PADDING on both ends so SIMD score-profile loaders can
    /// over-read by up to 32 bytes. The simple Rust SequenceSet does NOT
    /// have that padding. Today's live blastp pipeline survives because
    /// every SIMD consumer copies sequences into a fresh `Vec` before
    /// loading (see `swipe_set`, `WorkTarget::new`, `make_profile8` callers
    /// in `gapped_filter`). If any future caller passes `block.seqs().get(i)`
    /// directly into a SIMD profile builder or `_mm256_loadu_si256`, it can
    /// read uninitialised memory. See [[block-set-padding]] in memory.
    offsets: Vec<usize>,
}

impl SequenceSet {
    pub fn new() -> Self {
        SequenceSet {
            data: vec![DELIMITER_LETTER], // sentinel at start
            offsets: vec![1],             // first sequence starts after sentinel
        }
    }

    /// Add a sequence to the set.
    pub fn push(&mut self, seq: &[Letter]) {
        self.data.extend_from_slice(seq);
        self.data.push(DELIMITER_LETTER);
        self.offsets.push(self.data.len());
    }

    /// Matches C++ `StringSetBase::fill(n, value)`.
    pub fn fill(&mut self, n: usize, value: Letter) {
        self.data.extend(std::iter::repeat_n(value, n));
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

    /// Matches C++ `StringSetBase::length(i)`.
    pub fn length(&self, i: usize) -> Loc {
        self.seq_length(i)
    }

    /// Matches C++ `StringSetBase::ptr(i)`.
    pub fn ptr(&self, i: usize) -> usize {
        self.offsets[i]
    }

    /// Matches C++ `StringSetBase::limits_begin()`.
    pub fn offsets(&self) -> &[usize] {
        &self.offsets
    }

    /// Matches C++ `StringSetBase::position(i, j)`.
    pub fn position(&self, i: usize, j: usize) -> usize {
        self.offsets[i] + j
    }

    /// Matches C++ `StringSetBase::local_position(p)`.
    pub fn local_position(&self, p: usize) -> (usize, usize) {
        let i = match self.offsets.binary_search(&p) {
            Ok(i) => i,
            Err(i) => i - 1,
        };
        assert!(i < self.len());
        (i, p - self.offsets[i])
    }

    /// Matches C++ `SequenceSet::len_bounds(min_len)`.
    pub fn len_bounds(&self, min_len: Loc) -> (Loc, Loc) {
        let mut min = Loc::MAX;
        let mut max = 0;
        for i in 0..self.len() {
            let l = self.seq_length(i);
            if l < min_len {
                continue;
            }
            min = min.min(l);
            max = max.max(l);
        }
        (min, max)
    }

    /// Matches C++ `SequenceSet::max_len(begin, end)`.
    pub fn max_len(&self, begin: u32, end: u32) -> Loc {
        let mut max_len = 0;
        for i in begin..end {
            max_len = max_len.max(self.seq_length(i as usize));
        }
        max_len
    }

    /// Matches C++ `SequenceSet::partition(n_part, shortened, context_reduced)`.
    pub fn partition(&self, n_part: u32, shortened: bool, context_reduced: bool) -> Vec<u32> {
        self.partition_with_contexts(n_part, shortened, context_reduced, 1)
    }

    /// Matches C++ `SequenceSet::partition(n_part, shortened, context_reduced)`.
    pub fn partition_with_contexts(
        &self,
        n_part: u32,
        shortened: bool,
        context_reduced: bool,
        query_contexts: usize,
    ) -> Vec<u32> {
        assert!(n_part > 0);
        let target_letters = (self.letters() + n_part as u64 - 1) / n_part as u64;
        let contexts = if context_reduced { query_contexts } else { 1 };
        let mut partitions = Vec::with_capacity(n_part as usize + 1);
        if !shortened {
            partitions.push(0);
        }
        let mut i = 0usize;
        while i < self.len() {
            let mut letters = 0u64;
            while i < self.len() && letters < target_letters {
                for _ in 0..contexts {
                    if i >= self.len() {
                        break;
                    }
                    letters += self.seq_length(i) as u64;
                    i += 1;
                }
            }
            partitions.push((i / contexts) as u32);
        }
        let target_len = n_part as usize + if shortened { 0 } else { 1 };
        while partitions.len() < target_len {
            partitions.push((self.len() / contexts) as u32);
        }
        partitions
    }

    /// Matches C++ `SequenceSet::reverse_translated_len(i)`.
    pub fn reverse_translated_len(&self, i: usize) -> Loc {
        let j = i - i % 6;
        let l = self.seq_length(j);
        if self.seq_length(j + 2) == l {
            l * 3 + 2
        } else if self.seq_length(j + 1) == l {
            l * 3 + 1
        } else {
            l * 3
        }
    }

    /// Matches C++ `SequenceSet::avg_len()`.
    pub fn avg_len(&self) -> usize {
        self.letters() as usize / self.len()
    }

    /// Matches C++ `SequenceSet::lengths()`.
    pub fn lengths(&self) -> Vec<(Loc, u32)> {
        let mut v = Vec::with_capacity(self.len());
        for i in 0..self.len() {
            v.push((self.seq_length(i), i as u32));
        }
        v
    }

    /// Matches C++ `SequenceSet::source_length(i)`.
    pub fn source_length(&self, i: usize) -> Loc {
        self.source_length_with_contexts(i, 1)
    }

    /// Matches C++ `SequenceSet::source_length(i)`.
    pub fn source_length_with_contexts(&self, i: usize, query_contexts: usize) -> Loc {
        if query_contexts == 1 {
            self.seq_length(i)
        } else {
            let j = i - i % query_contexts;
            self.seq_length(j) + self.seq_length(j + 1) + self.seq_length(j + 2) + 2
        }
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

    /// Matches C++ `StringSetBase::data(p)`.
    pub fn data_at(&self, p: u64) -> &[Letter] {
        &self.data[p as usize..]
    }

    /// Matches C++ `StringSetBase::data(p)`.
    pub fn data_mut_at(&mut self, p: u64) -> &mut [Letter] {
        &mut self.data[p as usize..]
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
        ss.push(&[4, 5, 6]); // CQE

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

    #[test]
    fn test_sequence_set_positions() {
        let mut ss = SequenceSet::new();
        ss.push(&[0, 1, 2, 3]);
        ss.push(&[4, 5, 6]);

        assert_eq!(ss.ptr(0), 1);
        assert_eq!(ss.ptr(1), 6);
        assert_eq!(ss.position(0, 2), 3);
        assert_eq!(ss.position(1, 1), 7);
        assert_eq!(ss.local_position(1), (0, 0));
        assert_eq!(ss.local_position(4), (0, 3));
        assert_eq!(ss.local_position(6), (1, 0));
        assert_eq!(ss.local_position(8), (1, 2));
    }

    #[test]
    fn test_sequence_set_lengths_and_partitions() {
        let mut ss = SequenceSet::new();
        ss.push(&[0, 1, 2, 3]);
        ss.push(&[4, 5, 6]);
        ss.push(&[7, 8]);

        assert_eq!(ss.length(1), 3);
        assert_eq!(ss.len_bounds(3), (3, 4));
        assert_eq!(ss.len_bounds(5), (Loc::MAX, 0));
        assert_eq!(ss.max_len(0, 2), 4);
        assert_eq!(ss.avg_len(), 3);
        assert_eq!(ss.lengths(), vec![(4, 0), (3, 1), (2, 2)]);
        assert_eq!(ss.source_length(0), 4);
        assert_eq!(ss.partition(2, false, false), vec![0, 2, 3]);
        assert_eq!(ss.partition(2, false, true), vec![0, 2, 3]);
    }

    #[test]
    fn test_sequence_set_reverse_translated_lengths() {
        let mut ss = SequenceSet::new();
        ss.push(&[0, 1, 2]);
        ss.push(&[0, 1, 2]);
        ss.push(&[0, 1, 2]);
        ss.push(&[0, 1]);
        ss.push(&[0, 1]);
        ss.push(&[0, 1]);
        assert_eq!(ss.reverse_translated_len(0), 11);
        assert_eq!(ss.source_length_with_contexts(4, 6), 11);
    }

    #[test]
    fn test_sequence_set_partition_reduces_translated_contexts() {
        let mut ss = SequenceSet::new();
        for _ in 0..12 {
            ss.push(&[0]);
        }
        assert_eq!(ss.partition_with_contexts(2, false, true, 6), vec![0, 1, 2]);
        assert_eq!(ss.partition_with_contexts(2, true, true, 6), vec![1, 2]);
    }

    #[test]
    fn test_string_set_base_push_and_positions() {
        let mut ids = StringSet::new();
        assert!(ids.empty());
        assert_eq!(ids.raw_len(), 256);

        ids.push_back(b"alpha");
        ids.push_back(b"b longer");

        assert_eq!(ids.size(), 2);
        assert_eq!(ids.get(0), b"alpha");
        assert_eq!(ids.get(1), b"b longer");
        assert_eq!(ids.back(), b"b longer");
        assert_eq!(ids.length(0), 5);
        assert_eq!(ids.length(1), 8);
        assert_eq!(ids.letters(), 13);
        assert_eq!(ids.ptr(0), 256);
        assert_eq!(ids.ptr(1), 262);
        assert_eq!(ids.end(0), 261);
        assert_eq!(ids.position(1, 0), 262);
        assert_eq!(ids.local_position(259), (0, 3));
        assert_eq!(ids.local_position(262), (1, 0));
        assert_eq!(max_id_len(&ids), 5);
    }

    #[test]
    fn test_string_set_base_reserve_assign_subset_append() {
        let mut reserved = StringSet::new();
        reserved.reserve(3);
        reserved.reserve(2);
        reserved.finish_reserve();
        reserved.assign(0, b"cat");
        reserved.assign(1, b"ox");

        assert_eq!(reserved.get(0), b"cat");
        assert_eq!(reserved.get(1), b"ox");

        let subset = reserved.subset(&[1, 0]);
        assert_eq!(subset.get(0), b"ox");
        assert_eq!(subset.get(1), b"cat");

        let mut lhs = StringSet::new();
        lhs.push_back(b"aa");
        lhs.finish_reserve();
        lhs.append(&subset);

        assert_eq!(lhs.size(), 3);
        assert_eq!(lhs.get(0), b"aa");
        assert_eq!(lhs.get(1), b"ox");
        assert_eq!(lhs.get(2), b"cat");
    }

    #[test]
    fn test_letter_string_set_fill_and_batch_position() {
        let mut letters = LetterStringSet::new();
        letters.fill(3, 7);
        letters.push_back(&[1, 2]);

        assert_eq!(letters.get(0), &[7, 7, 7]);
        assert_eq!(letters.get(1), &[1, 2]);
        assert_eq!(
            letters.local_position_batch(&[256, 260], |a, b| a < b),
            vec![0, 1]
        );
    }
}
