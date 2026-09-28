//! Contiguous sequence storage translated from `data/sequence_set.cpp`,
//! `data/sequence_set.h`, and the inherited `data/string_set.h` template.

use std::io::{self, Write};
use std::sync::{Arc, OnceLock};

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
/// with an offset array for O(1) random access to any sequence. As in the C++
/// `StringSetBase`, the first sequence begins after 256 delimiter bytes. The
/// Rust implementation eagerly restores the trailing 256-byte perimeter after
/// mutating operations, which is equivalent to the C++ finalized state and
/// keeps direct SIMD consumers safe.
#[derive(Debug)]
pub struct SequenceSet {
    /// Raw data: all sequences concatenated with delimiter separators.
    data: Vec<Letter>,
    /// Offsets into data where each sequence starts.
    /// offsets[i] is the start of sequence i, offsets[i+1]-1 is the end (exclusive of delimiter).
    offsets: Vec<usize>,
    /// Lazily materialized owned sequence views.  Alignment objects outlive
    /// individual extension stages, so one Arc per database sequence avoids
    /// copying the same letters for every query/target candidate.
    shared: Vec<OnceLock<Arc<[Letter]>>>,
}

impl Clone for SequenceSet {
    fn clone(&self) -> Self {
        Self {
            data: self.data.clone(),
            offsets: self.offsets.clone(),
            shared: std::iter::repeat_with(OnceLock::new)
                .take(self.len())
                .collect(),
        }
    }
}

impl PartialEq for SequenceSet {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data && self.offsets == other.offsets
    }
}

impl Eq for SequenceSet {}

/// Explicit replacement for the C++ `align_mode` fields used by this file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SequenceSetConfig {
    pub query_contexts: usize,
    pub query_translated: bool,
}

impl Default for SequenceSetConfig {
    fn default() -> Self {
        Self {
            query_contexts: 1,
            query_translated: false,
        }
    }
}

/// Borrowed equivalent of C++ `TranslatedSequence` for sequence-set views.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TranslatedSequenceView<'a> {
    source: &'a [Letter],
    translated: [&'a [Letter]; 6],
}

impl<'a> TranslatedSequenceView<'a> {
    pub fn source(&self) -> &'a [Letter] {
        self.source
    }

    pub fn frame(&self, frame: usize) -> &'a [Letter] {
        self.translated[frame]
    }
}

impl SequenceSet {
    pub const PERIMETER_PADDING: usize = 256;

    pub fn new() -> Self {
        SequenceSet {
            data: vec![DELIMITER_LETTER; Self::PERIMETER_PADDING],
            offsets: vec![Self::PERIMETER_PADDING],
            shared: Vec::new(),
        }
    }

    fn truncate_trailing_padding(&mut self) {
        self.data.truncate(self.raw_len());
    }

    /// Matches C++ `StringSetBase::finish_reserve()`.
    pub fn finish_reserve(&mut self) {
        self.data
            .resize(self.raw_len() + Self::PERIMETER_PADDING, DELIMITER_LETTER);
    }

    /// Reserve one sequence of `n` letters, to be populated with `assign`.
    pub fn reserve(&mut self, n: usize) {
        self.truncate_trailing_padding();
        self.offsets.push(self.raw_len() + n + 1);
        self.shared.push(OnceLock::new());
    }

    /// Matches the two-argument C++ `StringSetBase::reserve` overload.
    pub fn reserve_capacity(&mut self, entries: usize, letters: usize) {
        self.offsets.reserve(entries + 1);
        self.shared.reserve(entries);
        self.data
            .reserve(letters + 2 * Self::PERIMETER_PADDING + entries);
    }

    pub fn assign(&mut self, i: usize, seq: &[Letter]) {
        assert_eq!(seq.len(), self.seq_length(i) as usize);
        if self.data.len() < self.raw_len() + Self::PERIMETER_PADDING {
            self.finish_reserve();
        }
        let begin = self.ptr(i);
        let end = begin + seq.len();
        self.shared[i].take();
        self.data[begin..end].copy_from_slice(seq);
        self.data[end] = DELIMITER_LETTER;
    }

    /// Add a sequence to the set.
    pub fn push(&mut self, seq: &[Letter]) {
        self.truncate_trailing_padding();
        self.data.extend_from_slice(seq);
        self.data.push(DELIMITER_LETTER);
        self.offsets.push(self.data.len());
        self.shared.push(OnceLock::new());
        self.finish_reserve();
    }

    /// Matches C++ `StringSetBase::fill(n, value)`.
    pub fn fill(&mut self, n: usize, value: Letter) {
        self.truncate_trailing_padding();
        self.data.extend(std::iter::repeat_n(value, n));
        self.data.push(DELIMITER_LETTER);
        self.offsets.push(self.data.len());
        self.shared.push(OnceLock::new());
        self.finish_reserve();
    }

    pub fn append(&mut self, other: &Self) {
        for i in 0..other.len() {
            self.push(other.get(i));
        }
    }

    pub fn clear(&mut self) {
        self.offsets.clear();
        self.offsets.push(Self::PERIMETER_PADDING);
        self.shared.clear();
        self.data.resize(Self::PERIMETER_PADDING, DELIMITER_LETTER);
    }

    pub fn shrink_to_fit(&mut self) {
        self.offsets.shrink_to_fit();
        self.shared.shrink_to_fit();
        self.data.shrink_to_fit();
    }

    pub fn subset(&self, indices: &[usize]) -> Self {
        let mut result = Self::new();
        result.reserve_capacity(indices.len(), 0);
        for &i in indices {
            result.reserve(self.seq_length(i) as usize);
        }
        result.finish_reserve();
        for (dst, &src) in indices.iter().enumerate() {
            result.assign(dst, self.get(src));
        }
        result
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

    /// Return one shared copy of sequence `i`, reused by all later callers.
    /// Mutating that sequence invalidates the cache for future callers while
    /// existing owners retain the immutable pre-mutation snapshot.
    pub fn shared_get(&self, i: usize) -> Arc<[Letter]> {
        self.shared[i]
            .get_or_init(|| Arc::from(self.get(i)))
            .clone()
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

    pub fn check_idx(&self, i: usize) -> usize {
        assert!(i < self.len(), "Sequence set index out of bounds.");
        i
    }

    pub fn raw_len(&self) -> usize {
        *self.offsets.last().unwrap()
    }

    pub fn mem_size(&self) -> usize {
        self.data.len() * std::mem::size_of::<Letter>()
            + self.offsets.len() * std::mem::size_of::<usize>()
    }

    pub fn back(&self) -> &[Letter] {
        self.get(self.len() - 1)
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
        (i, p - self.offsets[i])
    }

    /// Matches C++ `SequenceSet::len_bounds(min_len)`.
    pub fn len_bounds(&self, min_len: Loc) -> (Loc, Loc) {
        let mut min = Loc::MAX;
        let mut max = 0;
        for i in 0..self.len() {
            let l = self.seq_length(i);
            max = max.max(l);
            if l >= min_len {
                min = min.min(l);
            }
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
        assert!(query_contexts > 0);
        let target_letters = (self.letters() + n_part as u64 - 1) / n_part as u64;
        let contexts = if context_reduced { query_contexts } else { 1 };
        if context_reduced {
            assert_eq!(
                self.len() % contexts,
                0,
                "translated sequence count must be a multiple of query contexts"
            );
        }
        let mut partitions = Vec::with_capacity(n_part as usize + 1);
        if !shortened {
            partitions.push(0);
        }
        let mut i = 0usize;
        while i < self.len() {
            let mut letters = 0u64;
            while i < self.len() && letters < target_letters {
                for _ in 0..contexts {
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

    /// Matches C++ `SequenceSet::translated_seq`, with the former global
    /// `align_mode.query_translated` supplied explicitly.
    pub fn translated_seq<'a>(
        &'a self,
        source: &'a [Letter],
        i: usize,
        query_translated: bool,
    ) -> TranslatedSequenceView<'a> {
        if !query_translated {
            let seq = self.get(i);
            return TranslatedSequenceView {
                source: seq,
                translated: [seq, &[], &[], &[], &[], &[]],
            };
        }
        TranslatedSequenceView {
            source,
            translated: [
                self.get(i),
                self.get(i + 1),
                self.get(i + 2),
                self.get(i + 3),
                self.get(i + 4),
                self.get(i + 5),
            ],
        }
    }

    pub fn translated_seq_with_config<'a>(
        &'a self,
        source: &'a [Letter],
        i: usize,
        config: SequenceSetConfig,
    ) -> TranslatedSequenceView<'a> {
        self.translated_seq(source, i, config.query_translated)
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
        assert!(query_contexts > 0);
        if query_contexts == 1 {
            self.seq_length(i)
        } else {
            let j = i - i % query_contexts;
            self.seq_length(j) + self.seq_length(j + 1) + self.seq_length(j + 2) + 2
        }
    }

    pub fn source_length_with_config(&self, i: usize, config: SequenceSetConfig) -> Loc {
        self.source_length_with_contexts(i, config.query_contexts)
    }

    pub fn partition_with_config(
        &self,
        n_part: u32,
        shortened: bool,
        context_reduced: bool,
        config: SequenceSetConfig,
    ) -> Vec<u32> {
        self.partition_with_contexts(n_part, shortened, context_reduced, config.query_contexts)
    }

    /// C++ `print_stats` rendered without the process-global verbose stream.
    pub fn stats_line(&self) -> String {
        format!(
            "Sequences = {}, letters = {}, average length = {}\n",
            self.len(),
            self.letters(),
            self.avg_len()
        )
    }

    pub fn write_stats<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(self.stats_line().as_bytes())
    }

    /// Total letters across all sequences.
    pub fn letters(&self) -> u64 {
        (self.raw_len() - self.len() - Self::PERIMETER_PADDING) as u64
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
        for cached in &mut self.shared {
            cached.take();
        }
        &mut self.data[p as usize..]
    }

    /// Get mutable access to a sequence's data.
    pub fn get_mut(&mut self, i: usize) -> &mut [Letter] {
        let start = self.offsets[i];
        let end = self.offsets[i + 1] - 1;
        self.shared[i].take();
        &mut self.data[start..end]
    }
}

impl Default for SequenceSet {
    fn default() -> Self {
        Self::new()
    }
}

/// C++ exposes a move constructor from the underlying letter string set.
impl From<LetterStringSet> for SequenceSet {
    fn from(storage: LetterStringSet) -> Self {
        let LetterStringSet { data, limits } = storage;
        let sequence_count = limits.len().saturating_sub(1);
        Self {
            data,
            offsets: limits.into_iter().map(|offset| offset as usize).collect(),
            shared: std::iter::repeat_with(OnceLock::new)
                .take(sequence_count)
                .collect(),
        }
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
    fn shared_sequence_is_reused_and_get_mut_invalidates_only_its_entry() {
        let mut ss = SequenceSet::new();
        ss.push(&[0, 1, 2]);
        ss.push(&[3, 4]);
        let old0 = ss.shared_get(0);
        let old1 = ss.shared_get(1);
        assert!(Arc::ptr_eq(&old0, &ss.shared_get(0)));
        assert!(Arc::ptr_eq(&old1, &ss.shared_get(1)));

        ss.get_mut(0)[1] = 9;
        let new0 = ss.shared_get(0);
        assert_eq!(&*old0, &[0, 1, 2]);
        assert_eq!(&*new0, &[0, 9, 2]);
        assert!(!Arc::ptr_eq(&old0, &new0));
        assert!(Arc::ptr_eq(&old1, &ss.shared_get(1)));
    }

    #[test]
    fn bulk_mutation_clone_clear_and_reserve_maintain_shared_cache() {
        let mut ss = SequenceSet::new();
        ss.reserve(3);
        ss.finish_reserve();
        ss.assign(0, &[1, 2, 3]);
        let old = ss.shared_get(0);
        ss.data_mut_at(ss.ptr(0) as u64)[0] = 7;
        let changed = ss.shared_get(0);
        assert_eq!(&*old, &[1, 2, 3]);
        assert_eq!(&*changed, &[7, 2, 3]);
        assert!(!Arc::ptr_eq(&old, &changed));

        let cloned = ss.clone();
        let cloned_shared = cloned.shared_get(0);
        assert_eq!(&*cloned_shared, &[7, 2, 3]);
        assert!(!Arc::ptr_eq(&changed, &cloned_shared));

        ss.clear();
        ss.push(&[5]);
        assert_eq!(&*ss.shared_get(0), &[5]);
    }

    #[test]
    fn test_sequence_set_positions() {
        let mut ss = SequenceSet::new();
        ss.push(&[0, 1, 2, 3]);
        ss.push(&[4, 5, 6]);

        assert_eq!(ss.ptr(0), 256);
        assert_eq!(ss.ptr(1), 261);
        assert_eq!(ss.position(0, 2), 258);
        assert_eq!(ss.position(1, 1), 262);
        assert_eq!(ss.local_position(256), (0, 0));
        assert_eq!(ss.local_position(259), (0, 3));
        assert_eq!(ss.local_position(261), (1, 0));
        assert_eq!(ss.local_position(263), (1, 2));
    }

    #[test]
    fn test_sequence_set_lengths_and_partitions() {
        let mut ss = SequenceSet::new();
        ss.push(&[0, 1, 2, 3]);
        ss.push(&[4, 5, 6]);
        ss.push(&[7, 8]);

        assert_eq!(ss.length(1), 3);
        assert_eq!(ss.len_bounds(3), (3, 4));
        // C++ applies the threshold only to the minimum; maximum is global.
        assert_eq!(ss.len_bounds(5), (Loc::MAX, 4));
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
    fn sequence_set_has_cpp_perimeter_and_reserve_storage() {
        let mut ss = SequenceSet::new();
        assert_eq!(ss.raw_len(), SequenceSet::PERIMETER_PADDING);
        assert!(ss.data().iter().all(|&x| x == DELIMITER_LETTER));

        ss.reserve(3);
        ss.reserve(2);
        ss.finish_reserve();
        ss.assign(0, &[1, 2, 3]);
        ss.assign(1, &[4, 5]);

        assert_eq!(ss.offsets(), &[256, 260, 263]);
        assert_eq!(ss.get(0), &[1, 2, 3]);
        assert_eq!(ss.get(1), &[4, 5]);
        assert_eq!(ss.data()[255], DELIMITER_LETTER);
        assert_eq!(ss.data()[259], DELIMITER_LETTER);
        assert!(ss.data()[ss.raw_len()..]
            .iter()
            .all(|&x| x == DELIMITER_LETTER));

        let subset = ss.subset(&[1, 0]);
        assert_eq!(subset.get(0), &[4, 5]);
        assert_eq!(subset.get(1), &[1, 2, 3]);

        let mut appended = SequenceSet::new();
        appended.push(&[9]);
        appended.append(&subset);
        assert_eq!(appended.get(0), &[9]);
        assert_eq!(appended.get(1), &[4, 5]);
        assert_eq!(appended.get(2), &[1, 2, 3]);
    }

    #[test]
    fn sequence_set_partition_and_length_sort_match_cpp() {
        let mut ss = SequenceSet::new();
        for n in [2, 2, 3, 3, 1, 1] {
            ss.fill(n, 7);
        }
        let config = SequenceSetConfig {
            query_contexts: 2,
            query_translated: true,
        };
        assert_eq!(
            ss.partition_with_config(3, false, true, config),
            vec![0, 1, 2, 3]
        );
        assert_eq!(
            ss.partition_with_config(5, true, true, config),
            vec![1, 2, 3, 3, 3]
        );

        let mut lengths = ss.lengths();
        lengths.sort_unstable();
        assert_eq!(
            lengths,
            vec![(1, 4), (1, 5), (2, 0), (2, 1), (3, 2), (3, 3)]
        );
    }

    #[test]
    fn translated_view_and_stats_are_explicit() {
        let mut ss = SequenceSet::new();
        for frame in 0..6 {
            ss.push(&[frame as Letter]);
        }
        let source = [10, 11, 12, 13];
        let translated = ss.translated_seq(&source, 0, true);
        assert_eq!(translated.source(), &source);
        for frame in 0..6 {
            assert_eq!(translated.frame(frame), &[frame as Letter]);
        }

        let protein = ss.translated_seq(&source, 2, false);
        assert_eq!(protein.source(), &[2]);
        assert_eq!(protein.frame(0), &[2]);
        assert!(protein.frame(1).is_empty());

        let mut out = Vec::new();
        ss.write_stats(&mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "Sequences = 6, letters = 6, average length = 1\n"
        );
    }

    #[test]
    fn sequence_set_converts_from_cpp_base_storage() {
        let mut storage = LetterStringSet::new();
        storage.push_back(&[3, 4]);
        storage.finish_reserve();
        let ss = SequenceSet::from(storage);
        assert_eq!(ss.get(0), &[3, 4]);
        assert_eq!(ss.ptr(0), 256);
        assert_eq!(ss.data().len(), ss.raw_len() + 256);
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
