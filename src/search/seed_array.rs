//! Partitioned seed array matching C++ DIAMOND's SeedArray.
//!
//! Seeds are partitioned by their low bits (seedp_bits) into buckets.
//! Within each partition, entries are stored contiguously in a flat array.
//! This gives cache-friendly access during the join phase and enables
//! per-partition parallelism.

use crate::basic::reduction::Reduction;
use crate::basic::seed::{self, PackedSeed, SeedPartition};
use crate::basic::shape::Shape;
use crate::basic::value::Letter;
use crate::search::seed_match::SeedMatch;

/// A single entry in the seed array: the seed key (with partition bits stripped)
/// and the location it came from.
#[derive(Clone, Copy, Debug)]
pub struct SeedEntry {
    /// Seed key with low seedp_bits removed (the "offset" part).
    pub key: u32,
    /// Global location in the sequence set, using unpadded cumulative sequence
    /// lengths. This preserves `(seq_id, pos)` ordering while keeping each seed
    /// entry at 8 bytes instead of 12.
    pub loc: u32,
}

/// Partitioned seed array. Seeds are grouped by `seed & seedp_mask` into partitions.
/// Within each partition, entries are stored contiguously.
pub struct SeedArray {
    /// Flat storage of all entries, grouped by partition.
    data: Vec<SeedEntry>,
    /// offsets[p] = start index of partition p in data. Length = num_partitions + 1.
    offsets: Vec<usize>,
    /// Sequence start offsets for decoding `SeedEntry::loc`.
    seq_offsets: Vec<u32>,
    /// Number of partition bits.
    pub seedp_bits: i32,
}

impl SeedArray {
    fn sequence_offsets(seqs: &[&[Letter]]) -> Vec<u32> {
        let mut offsets = Vec::with_capacity(seqs.len() + 1);
        let mut total = 0usize;
        offsets.push(0);
        for seq in seqs {
            total = total
                .checked_add(seq.len())
                .expect("sequence set length overflow");
            offsets.push(u32::try_from(total).expect("seed array compact location overflow"));
        }
        offsets
    }

    #[inline]
    pub fn seq_pos(&self, entry: SeedEntry) -> (u32, u32) {
        decode_seq_pos(&self.seq_offsets, entry)
    }

    pub fn seq_offsets(&self) -> &[u32] {
        &self.seq_offsets
    }

    /// Build a partitioned seed array from sequences.
    ///
    /// Two-pass approach matching C++ SeedArray:
    /// 1. Count seeds per partition (histogram)
    /// 2. Allocate and fill
    pub fn build(
        seqs: &[&[Letter]],
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
    ) -> Self {
        Self::build_with_complexity_cut(seqs, shape, reduction, seedp_bits, 0.0)
    }

    /// Build a seed index, optionally filtering out seeds whose unreduced
    /// composition has entropy below `complexity_cut`. This mirrors
    /// C++ DIAMOND's behavior, which masks low-complexity seeds during the
    /// seed-index construction (per-sensitivity `seed_cut` in
    /// `search/setup.cpp` * ln(2) * shape.weight).
    ///
    /// Parallel across input sequences. Each worker emits seeds into its own
    /// per-partition `Vec<SeedEntry>`; after the parallel phase the per-worker
    /// partitions are concatenated into the final flat array. This matches
    /// C++ DIAMOND's `OnePassBufferedWriter` model and removes the serial
    /// two-pass enumeration that previously dominated `blastp`'s wall time.
    pub fn build_with_complexity_cut(
        seqs: &[&[Letter]],
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
        complexity_cut: f64,
    ) -> Self {
        Self::build_with_complexity_cut_and_min_query_len(
            seqs,
            shape,
            reduction,
            seedp_bits,
            complexity_cut,
            0,
        )
    }

    pub fn build_with_complexity_cut_and_min_query_len(
        seqs: &[&[Letter]],
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
        complexity_cut: f64,
        min_query_len: usize,
    ) -> Self {
        Self::build_with_complexity_cut_and_min_query_len_partition_range(
            seqs,
            shape,
            reduction,
            seedp_bits,
            complexity_cut,
            min_query_len,
            0,
            seed::seedp_count(seedp_bits) as usize,
        )
    }

    pub fn build_with_complexity_cut_and_min_query_len_partition_range(
        seqs: &[&[Letter]],
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
        complexity_cut: f64,
        min_query_len: usize,
        partition_begin: usize,
        partition_end: usize,
    ) -> Self {
        Self::build_with_complexity_cut_and_min_query_len_partition_range_reuse(
            seqs,
            shape,
            reduction,
            seedp_bits,
            complexity_cut,
            min_query_len,
            partition_begin,
            partition_end,
            Vec::new(),
        )
    }

    pub fn build_with_complexity_cut_and_min_query_len_partition_range_reuse(
        seqs: &[&[Letter]],
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
        complexity_cut: f64,
        min_query_len: usize,
        partition_begin: usize,
        partition_end: usize,
        mut data: Vec<SeedEntry>,
    ) -> Self {
        use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

        let num_partitions = seed::seedp_count(seedp_bits) as usize;
        debug_assert!(partition_begin <= partition_end);
        debug_assert!(partition_end <= num_partitions);
        let filter_partitions = partition_begin != 0 || partition_end != num_partitions;
        let mask = seed::seedp_mask(seedp_bits);
        let shape_length = shape.length as usize;
        let seq_offsets = Self::sequence_offsets(seqs);

        if rayon::current_num_threads() == 1 {
            let mut counts = vec![0usize; num_partitions];
            for seq in seqs {
                let slen = seq.len();
                if min_query_len > 0 && slen < min_query_len {
                    continue;
                }
                if slen < shape_length {
                    continue;
                }
                let last = slen - shape_length;
                for pos in 0..=last {
                    let window = &seq[pos..pos + shape_length];
                    let seed = if complexity_cut > 0.0 {
                        shape.set_seed_with_complexity(window, reduction, complexity_cut)
                    } else {
                        shape.set_seed(window, reduction)
                    };
                    if let Some(s) = seed {
                        let p = seed::seed_partition(s, mask) as usize;
                        if filter_partitions && (p < partition_begin || p >= partition_end) {
                            continue;
                        }
                        counts[p] += 1;
                    }
                }
            }

            let mut offsets = vec![0usize; num_partitions + 1];
            for p in 0..num_partitions {
                offsets[p + 1] = offsets[p]
                    .checked_add(counts[p])
                    .expect("seed array partition size overflow");
            }
            let total = offsets[num_partitions];
            let mut cursors = offsets[..num_partitions].to_vec();
            data.clear();
            if data.capacity() < total {
                data.reserve_exact(total);
            }
            let data_ptr = data.as_mut_ptr();
            for (seq_id, seq) in seqs.iter().enumerate() {
                let slen = seq.len();
                if min_query_len > 0 && slen < min_query_len {
                    continue;
                }
                if slen < shape_length {
                    continue;
                }
                let last = slen - shape_length;
                for pos in 0..=last {
                    let window = &seq[pos..pos + shape_length];
                    let seed = if complexity_cut > 0.0 {
                        shape.set_seed_with_complexity(window, reduction, complexity_cut)
                    } else {
                        shape.set_seed(window, reduction)
                    };
                    if let Some(s) = seed {
                        let p = seed::seed_partition(s, mask) as usize;
                        if filter_partitions && (p < partition_begin || p >= partition_end) {
                            continue;
                        }
                        let key = seed::seed_partition_offset(s, seedp_bits as u64) as u32;
                        let idx = cursors[p];
                        cursors[p] += 1;
                        // SAFETY: pass 1 counted the exact number of entries
                        // for every partition, and `idx` advances within the
                        // reserved range `[offsets[p], offsets[p + 1])`.
                        unsafe {
                            std::ptr::write(
                                data_ptr.add(idx),
                                SeedEntry {
                                    key,
                                    loc: seq_offsets[seq_id] + pos as u32,
                                },
                            );
                        }
                    }
                }
            }

            for p in 0..num_partitions {
                assert_eq!(cursors[p], offsets[p + 1]);
            }
            // SAFETY: the cursor checks above confirm that every slot in
            // `0..total` was initialized exactly once.
            unsafe {
                data.set_len(total);
            }
            return SeedArray {
                data,
                offsets,
                seq_offsets,
                seedp_bits,
            };
        }

        // Two-pass build matching C++'s `SeedArray::SeedArray(...with histogram)`:
        //   pass 1: scan all sequences in parallel, count seeds per partition
        //   pass 2: allocate one contiguous Vec sized to the total, then write
        //           each seed directly to its partition's slot via atomic cursor
        // This avoids the gigabyte-scale intermediate `Vec<Raw>` and the
        // sequential bucket-sort copy that the prior implementation needed.
        use std::sync::atomic::{AtomicUsize, Ordering};

        // Pass 1: histogram. Each worker keeps a `Vec<usize>` count per
        // partition; reduce sums them. Keeping counts in `usize` avoids
        // release-mode wrap before allocation on very large inputs.
        let counts: Vec<usize> = seqs
            .into_par_iter()
            .fold(
                || vec![0usize; num_partitions],
                |mut acc, seq| {
                    let slen = seq.len();
                    if min_query_len > 0 && slen < min_query_len {
                        return acc;
                    }
                    if slen < shape_length {
                        return acc;
                    }
                    let last = slen - shape_length;
                    for pos in 0..=last {
                        let window = &seq[pos..pos + shape_length];
                        let seed = if complexity_cut > 0.0 {
                            shape.set_seed_with_complexity(window, reduction, complexity_cut)
                        } else {
                            shape.set_seed(window, reduction)
                        };
                        if let Some(s) = seed {
                            let p = seed::seed_partition(s, mask) as usize;
                            if filter_partitions && (p < partition_begin || p >= partition_end) {
                                continue;
                            }
                            acc[p] += 1;
                        }
                    }
                    acc
                },
            )
            .reduce(
                || vec![0usize; num_partitions],
                |mut a, b| {
                    for (x, y) in a.iter_mut().zip(b.iter()) {
                        *x += *y;
                    }
                    a
                },
            );

        let mut offsets = vec![0usize; num_partitions + 1];
        for p in 0..num_partitions {
            offsets[p + 1] = offsets[p]
                .checked_add(counts[p])
                .expect("seed array partition size overflow");
        }
        let total = offsets[num_partitions];

        // Pass 2: allocate the flat buffer and write each seed directly to its
        // partition slot via an atomic cursor per partition. Using
        // `AtomicUsize::fetch_add(1)` gives each push a unique index — fast on
        // x86 (single LOCK XADD) and avoids serialization through a single
        // global cursor.
        let cursors: Vec<AtomicUsize> = offsets[..num_partitions]
            .iter()
            .map(|&o| AtomicUsize::new(o))
            .collect();
        data.clear();
        if data.capacity() < total {
            data.reserve_exact(total);
        }
        let data_ptr = data.as_mut_ptr() as usize; // smuggle pointer across threads

        seqs.into_par_iter().enumerate().for_each(|(seq_id, seq)| {
            let slen = seq.len();
            if min_query_len > 0 && slen < min_query_len {
                return;
            }
            if slen < shape_length {
                return;
            }
            let last = slen - shape_length;
            for pos in 0..=last {
                let window = &seq[pos..pos + shape_length];
                let seed = if complexity_cut > 0.0 {
                    shape.set_seed_with_complexity(window, reduction, complexity_cut)
                } else {
                    shape.set_seed(window, reduction)
                };
                if let Some(s) = seed {
                    let p = seed::seed_partition(s, mask) as usize;
                    if filter_partitions && (p < partition_begin || p >= partition_end) {
                        continue;
                    }
                    let key = seed::seed_partition_offset(s, seedp_bits as u64) as u32;
                    let idx = cursors[p].fetch_add(1, Ordering::Relaxed);
                    // SAFETY: `idx` is unique per `cursors[p]` and falls in
                    // `[offsets[p], offsets[p+1])`. data_ptr is the start of the
                    // properly-allocated buffer.
                    unsafe {
                        let ptr = (data_ptr as *mut SeedEntry).add(idx);
                        std::ptr::write(
                            ptr,
                            SeedEntry {
                                key,
                                loc: seq_offsets[seq_id] + pos as u32,
                            },
                        );
                    }
                }
            }
        });

        for p in 0..num_partitions {
            assert_eq!(
                cursors[p].load(Ordering::Relaxed),
                offsets[p + 1],
                "seed array partition {p} was not fully initialized"
            );
        }
        // SAFETY: The pass-1 counts determined `offsets`, each write above used
        // a unique cursor slot in `offsets[p]..offsets[p+1]`, and the cursor
        // check confirms every slot in `0..total` was initialized before the
        // vector length becomes visible to safe Rust.
        unsafe {
            data.set_len(total);
        }

        SeedArray {
            data,
            offsets,
            seq_offsets,
            seedp_bits,
        }
    }

    pub fn build_partition_range_one_pass(
        seqs: &[&[Letter]],
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
        min_query_len: usize,
        partition_begin: usize,
        partition_end: usize,
    ) -> Self {
        let num_partitions = seed::seedp_count(seedp_bits) as usize;
        debug_assert!(partition_begin <= partition_end);
        debug_assert!(partition_end <= num_partitions);
        let mask = seed::seedp_mask(seedp_bits);
        let shape_length = shape.length as usize;
        let seq_offsets = Self::sequence_offsets(seqs);
        let range_len = partition_end - partition_begin;
        let mut buckets: Vec<Vec<SeedEntry>> = (0..range_len).map(|_| Vec::new()).collect();

        for (seq_id, seq) in seqs.iter().enumerate() {
            let slen = seq.len();
            if min_query_len > 0 && slen < min_query_len {
                continue;
            }
            if slen < shape_length {
                continue;
            }
            let last = slen - shape_length;
            let seq_offset = seq_offsets[seq_id];
            for pos in 0..=last {
                let window = &seq[pos..pos + shape_length];
                if let Some(s) = shape.set_seed(window, reduction) {
                    let p = seed::seed_partition(s, mask) as usize;
                    if p < partition_begin || p >= partition_end {
                        continue;
                    }
                    let key = seed::seed_partition_offset(s, seedp_bits as u64) as u32;
                    buckets[p - partition_begin].push(SeedEntry {
                        key,
                        loc: seq_offset + pos as u32,
                    });
                }
            }
        }

        let mut offsets = vec![0usize; num_partitions + 1];
        let mut total = 0usize;
        for p in 0..num_partitions {
            if p >= partition_begin && p < partition_end {
                total = total
                    .checked_add(buckets[p - partition_begin].len())
                    .expect("seed array partition size overflow");
            }
            offsets[p + 1] = total;
        }

        let mut data = Vec::with_capacity(total);
        for bucket in buckets {
            data.extend(bucket);
        }

        SeedArray {
            data,
            offsets,
            seq_offsets,
            seedp_bits,
        }
    }

    pub fn count_partitions(
        seqs: &[&[Letter]],
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
        complexity_cut: f64,
        min_query_len: usize,
    ) -> Vec<usize> {
        let num_partitions = seed::seedp_count(seedp_bits) as usize;
        let mask = seed::seedp_mask(seedp_bits);
        let shape_length = shape.length as usize;
        let mut counts = vec![0usize; num_partitions];
        for seq in seqs {
            let slen = seq.len();
            if min_query_len > 0 && slen < min_query_len {
                continue;
            }
            if slen < shape_length {
                continue;
            }
            let last = slen - shape_length;
            for pos in 0..=last {
                let window = &seq[pos..pos + shape_length];
                let seed = if complexity_cut > 0.0 {
                    shape.set_seed_with_complexity(window, reduction, complexity_cut)
                } else {
                    shape.set_seed(window, reduction)
                };
                if let Some(s) = seed {
                    counts[seed::seed_partition(s, mask) as usize] += 1;
                }
            }
        }
        counts
    }

    pub fn build_partition_range_from_counts_reuse(
        seqs: &[&[Letter]],
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
        min_query_len: usize,
        partition_begin: usize,
        partition_end: usize,
        counts: &[usize],
        mut data: Vec<SeedEntry>,
    ) -> Self {
        let num_partitions = seed::seedp_count(seedp_bits) as usize;
        debug_assert_eq!(counts.len(), num_partitions);
        debug_assert!(partition_begin <= partition_end);
        debug_assert!(partition_end <= num_partitions);
        let mask = seed::seedp_mask(seedp_bits);
        let shape_length = shape.length as usize;
        let seq_offsets = Self::sequence_offsets(seqs);

        let mut offsets = vec![0usize; num_partitions + 1];
        let mut total = 0usize;
        for p in 0..num_partitions {
            if p >= partition_begin && p < partition_end {
                total = total
                    .checked_add(counts[p])
                    .expect("seed array partition size overflow");
            }
            offsets[p + 1] = total;
        }

        let mut cursors = offsets[..num_partitions].to_vec();
        data.clear();
        if data.capacity() < total {
            data.reserve_exact(total);
        }
        let data_ptr = data.as_mut_ptr();

        for (seq_id, seq) in seqs.iter().enumerate() {
            let slen = seq.len();
            if min_query_len > 0 && slen < min_query_len {
                continue;
            }
            if slen < shape_length {
                continue;
            }
            let last = slen - shape_length;
            let seq_offset = seq_offsets[seq_id];
            for pos in 0..=last {
                let window = &seq[pos..pos + shape_length];
                if let Some(s) = shape.set_seed(window, reduction) {
                    let p = seed::seed_partition(s, mask) as usize;
                    if p < partition_begin || p >= partition_end {
                        continue;
                    }
                    let key = seed::seed_partition_offset(s, seedp_bits as u64) as u32;
                    let idx = cursors[p];
                    cursors[p] += 1;
                    unsafe {
                        std::ptr::write(
                            data_ptr.add(idx),
                            SeedEntry {
                                key,
                                loc: seq_offset + pos as u32,
                            },
                        );
                    }
                }
            }
        }

        for p in partition_begin..partition_end {
            assert_eq!(cursors[p], offsets[p + 1]);
        }
        unsafe {
            data.set_len(total);
        }

        SeedArray {
            data,
            offsets,
            seq_offsets,
            seedp_bits,
        }
    }

    /// Get the entries for a given partition.
    #[inline]
    pub fn partition(&self, p: SeedPartition) -> &[SeedEntry] {
        let p = p as usize;
        &self.data[self.offsets[p]..self.offsets[p + 1]]
    }

    /// Get mutable entries for a given partition (for sorting).
    #[inline]
    pub fn partition_mut(&mut self, p: SeedPartition) -> &mut [SeedEntry] {
        let p = p as usize;
        let start = self.offsets[p];
        let end = self.offsets[p + 1];
        &mut self.data[start..end]
    }

    /// Borrow the underlying flat buffer mutably; the caller can use
    /// `split_at_mut` walks plus `partition_offset` to obtain disjoint
    /// per-partition mutable slices without copying.
    #[inline]
    pub fn data_mut(&mut self) -> &mut [SeedEntry] {
        &mut self.data
    }

    /// Start offset of partition `p` in the flat buffer. `partition_offset(n)`
    /// equals the total entry count.
    #[inline]
    pub fn partition_offset(&self, p: usize) -> usize {
        self.offsets[p]
    }

    /// Number of partitions.
    pub fn num_partitions(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Total number of entries.
    pub fn total_entries(&self) -> usize {
        self.data.len()
    }

    pub fn into_data(self) -> Vec<SeedEntry> {
        self.data
    }
}

#[inline]
pub fn decode_seq_pos(seq_offsets: &[u32], entry: SeedEntry) -> (u32, u32) {
    let loc = entry.loc;
    let i = seq_offsets.partition_point(|&x| x <= loc) - 1;
    (i as u32, loc - seq_offsets[i])
}

/// Result of joining one partition: matched (query_loc, ref_loc) pairs.
pub struct JoinResult {
    pub query_locs: Vec<(u32, u32)>, // (seq_id, pos)
    pub ref_locs: Vec<(u32, u32)>,   // (seq_id, pos)
}

/// Per-partition matched-key information used to compute frequent-seed
/// statistics. Each entry records the start index and count of matching
/// entries on the query and ref sides for one shared key, into the
/// already-sorted partition slices.
pub struct PartitionBlocks {
    pub blocks: Vec<MatchBlock>,
}

#[derive(Clone, Copy)]
pub struct MatchBlock {
    pub q_start: u32,
    pub q_count: u32,
    pub r_start: u32,
    pub r_count: u32,
}

/// Sort the partitions and find every block of matching keys. Returns the
/// blocks (one per shared key) suitable for accumulating frequent-seed stats
/// and then emitting the join with a threshold filter. The partitions are
/// modified in place — they are left sorted by key.
pub fn match_blocks(query_part: &mut [SeedEntry], ref_part: &mut [SeedEntry]) -> PartitionBlocks {
    if query_part.is_empty() || ref_part.is_empty() {
        return PartitionBlocks { blocks: Vec::new() };
    }
    // C++ radix sort is stable over entries emitted in sequence/position
    // order. The Rust seed-array builder fills partitions in parallel, so an
    // unstable sort by the full logical order is both deterministic and closer
    // to C++ than preserving atomic insertion order for equal keys.
    query_part.sort_unstable_by_key(|e| (e.key, e.loc));
    ref_part.sort_unstable_by_key(|e| (e.key, e.loc));
    let mut blocks = Vec::new();
    let mut qi = 0usize;
    let mut ri = 0usize;
    while qi < query_part.len() && ri < ref_part.len() {
        let qk = query_part[qi].key;
        let rk = ref_part[ri].key;
        if qk < rk {
            qi += 1;
        } else if qk > rk {
            ri += 1;
        } else {
            let qi_start = qi;
            while qi < query_part.len() && query_part[qi].key == qk {
                qi += 1;
            }
            let ri_start = ri;
            while ri < ref_part.len() && ref_part[ri].key == rk {
                ri += 1;
            }
            blocks.push(MatchBlock {
                q_start: qi_start as u32,
                q_count: (qi - qi_start) as u32,
                r_start: ri_start as u32,
                r_count: (ri - ri_start) as u32,
            });
        }
    }
    PartitionBlocks { blocks }
}

/// Emit the join result for an already-sorted partition, skipping blocks
/// whose query or ref count exceeds the given thresholds. Mirrors the
/// frequent-seed mask in `diamond/src/data/frequent_seeds.cpp`.
pub fn emit_join_filtered(
    query_part: &[SeedEntry],
    ref_part: &[SeedEntry],
    query_offsets: &[u32],
    ref_offsets: &[u32],
    blocks: &PartitionBlocks,
    q_max: u32,
    r_max: u32,
) -> JoinResult {
    let mut query_locs = Vec::new();
    let mut ref_locs = Vec::new();
    for b in &blocks.blocks {
        if b.q_count > q_max || b.r_count > r_max {
            continue;
        }
        for q in b.q_start..b.q_start + b.q_count {
            for r in b.r_start..b.r_start + b.r_count {
                query_locs.push(decode_seq_pos(query_offsets, query_part[q as usize]));
                ref_locs.push(decode_seq_pos(ref_offsets, ref_part[r as usize]));
            }
        }
    }
    JoinResult {
        query_locs,
        ref_locs,
    }
}

/// Emit `SeedMatch` values directly for an already-sorted partition, skipping
/// frequent-seed blocks. This preserves the same pair order as
/// `emit_join_filtered` without materializing two intermediate location
/// vectors.
pub fn emit_seed_matches_filtered(
    query_part: &[SeedEntry],
    ref_part: &[SeedEntry],
    query_offsets: &[u32],
    ref_offsets: &[u32],
    blocks: &PartitionBlocks,
    q_max: u32,
    r_max: u32,
) -> Vec<SeedMatch> {
    let capacity = blocks
        .blocks
        .iter()
        .filter(|b| b.q_count <= q_max && b.r_count <= r_max)
        .map(|b| b.q_count as usize * b.r_count as usize)
        .sum();
    let mut matches = Vec::with_capacity(capacity);
    for b in &blocks.blocks {
        if b.q_count > q_max || b.r_count > r_max {
            continue;
        }
        for q in b.q_start..b.q_start + b.q_count {
            let q_entry = query_part[q as usize];
            let (query_id, query_pos) = decode_seq_pos(query_offsets, q_entry);
            for r in b.r_start..b.r_start + b.r_count {
                let r_entry = ref_part[r as usize];
                let (ref_id, ref_pos) = decode_seq_pos(ref_offsets, r_entry);
                matches.push(SeedMatch {
                    query_id,
                    query_pos,
                    ref_id,
                    ref_pos,
                    seed: 0,
                    shape_id: 0,
                });
            }
        }
    }
    matches
}

/// Sort and emit `SeedMatch` values directly, without materializing
/// `MatchBlock`s. This is the normal C++ path when frequent-seed masking is
/// disabled: the hash join streams matching key groups straight to hits.
pub fn sort_merge_seed_matches(
    query_part: &mut [SeedEntry],
    ref_part: &mut [SeedEntry],
    query_offsets: &[u32],
    ref_offsets: &[u32],
    partition: u32,
    seedp_bits: i32,
) -> Vec<SeedMatch> {
    sort_merge_seed_matches_with_complexity(
        query_part,
        ref_part,
        query_offsets,
        ref_offsets,
        partition,
        seedp_bits,
        |_, _, _| true,
    )
}

/// Sort and emit `SeedMatch` values directly, allowing the caller to skip a
/// whole joined key block. The predicate receives the first query entry for the
/// block and the matching seed key.
pub fn sort_merge_seed_matches_with_complexity<F>(
    query_part: &mut [SeedEntry],
    ref_part: &mut [SeedEntry],
    query_offsets: &[u32],
    ref_offsets: &[u32],
    partition: u32,
    seedp_bits: i32,
    mut keep_query_block: F,
) -> Vec<SeedMatch>
where
    F: FnMut(SeedEntry, (u32, u32), u32) -> bool,
{
    if query_part.is_empty() || ref_part.is_empty() {
        return Vec::new();
    }

    query_part.sort_unstable_by_key(|e| (e.key, e.loc));
    ref_part.sort_unstable_by_key(|e| (e.key, e.loc));

    let mut matches = Vec::new();
    let mut qi = 0usize;
    let mut ri = 0usize;
    while qi < query_part.len() && ri < ref_part.len() {
        let qk = query_part[qi].key;
        let rk = ref_part[ri].key;
        if qk < rk {
            qi += 1;
        } else if qk > rk {
            ri += 1;
        } else {
            let qi_start = qi;
            while qi < query_part.len() && query_part[qi].key == qk {
                qi += 1;
            }
            let ri_start = ri;
            while ri < ref_part.len() && ref_part[ri].key == rk {
                ri += 1;
            }
            if !keep_query_block(
                query_part[qi_start],
                decode_seq_pos(query_offsets, query_part[qi_start]),
                qk,
            ) {
                continue;
            }
            let seed = ((qk as PackedSeed) << seedp_bits) | partition as PackedSeed;
            matches.reserve((qi - qi_start) * (ri - ri_start));
            for q_entry in &query_part[qi_start..qi] {
                let (query_id, query_pos) = decode_seq_pos(query_offsets, *q_entry);
                for r_entry in &ref_part[ri_start..ri] {
                    let (ref_id, ref_pos) = decode_seq_pos(ref_offsets, *r_entry);
                    matches.push(SeedMatch {
                        query_id,
                        query_pos,
                        ref_id,
                        ref_pos,
                        seed,
                        shape_id: 0,
                    });
                }
            }
        }
    }
    matches
}

/// Sort-merge join on one partition.
///
/// Sort both sides by key, then walk through matching keys.
/// For each matching key, emit the cross product of query and ref locations.
pub fn sort_merge_join(
    query_part: &mut [SeedEntry],
    ref_part: &mut [SeedEntry],
    query_offsets: &[u32],
    ref_offsets: &[u32],
) -> JoinResult {
    if query_part.is_empty() || ref_part.is_empty() {
        return JoinResult {
            query_locs: Vec::new(),
            ref_locs: Vec::new(),
        };
    }

    // Sort both by key and C++ logical emission order. See `match_blocks`.
    query_part.sort_unstable_by_key(|e| (e.key, e.loc));
    ref_part.sort_unstable_by_key(|e| (e.key, e.loc));

    let mut query_locs = Vec::new();
    let mut ref_locs = Vec::new();

    let mut qi = 0usize;
    let mut ri = 0usize;

    while qi < query_part.len() && ri < ref_part.len() {
        let qk = query_part[qi].key;
        let rk = ref_part[ri].key;

        if qk < rk {
            qi += 1;
        } else if qk > rk {
            ri += 1;
        } else {
            // Keys match — find the extent of equal keys on both sides
            let qi_start = qi;
            while qi < query_part.len() && query_part[qi].key == qk {
                qi += 1;
            }
            let ri_start = ri;
            while ri < ref_part.len() && ref_part[ri].key == rk {
                ri += 1;
            }

            // Emit cross product
            for q in qi_start..qi {
                for r in ri_start..ri {
                    query_locs.push(decode_seq_pos(query_offsets, query_part[q]));
                    ref_locs.push(decode_seq_pos(ref_offsets, ref_part[r]));
                }
            }
        }
    }

    JoinResult {
        query_locs,
        ref_locs,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seed_array_build() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);
        let seq: Vec<Letter> = (0..20).map(|i| (i % 20) as Letter).collect();
        let seqs: Vec<&[Letter]> = vec![&seq];

        let sa = SeedArray::build(&seqs, &shape, &reduction, 4);
        assert_eq!(sa.num_partitions(), 16);
        // Total entries should equal number of valid seed positions
        assert_eq!(sa.total_entries(), 18); // 20 - 3 + 1 = 18
    }

    #[test]
    fn test_seed_array_parallel_build_matches_serial_enumeration() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);
        let seedp_bits = 4;
        let mask = seed::seedp_mask(seedp_bits);
        let seqs_owned: Vec<Vec<Letter>> = (0..96)
            .map(|sid| {
                (0..24)
                    .map(|pos| ((sid + pos * 7) % 20) as Letter)
                    .collect()
            })
            .collect();
        let seqs: Vec<&[Letter]> = seqs_owned.iter().map(Vec::as_slice).collect();

        let sa = SeedArray::build(&seqs, &shape, &reduction, seedp_bits as i32);

        let mut expected = Vec::new();
        for (seq_id, seq) in seqs.iter().enumerate() {
            for pos in 0..=seq.len() - shape.length as usize {
                let window = &seq[pos..pos + shape.length as usize];
                if let Some(packed) = shape.set_seed(window, &reduction) {
                    expected.push((
                        seed::seed_partition(packed, mask) as usize,
                        seed::seed_partition_offset(packed, seedp_bits as u64) as u32,
                        seq_id as u32,
                        pos as u32,
                    ));
                }
            }
        }
        expected.sort_unstable();

        let mut actual = Vec::new();
        for p in 0..sa.num_partitions() {
            for entry in sa.partition(p as SeedPartition) {
                let (seq_id, pos) = sa.seq_pos(*entry);
                actual.push((p, entry.key, seq_id, pos));
            }
        }
        actual.sort_unstable();

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_sort_merge_join() {
        // Create two small partitions with some matching keys
        let mut query = vec![
            SeedEntry { key: 5, loc: 10 },
            SeedEntry { key: 3, loc: 20 },
            SeedEntry { key: 5, loc: 130 },
        ];
        let mut refs = vec![
            SeedEntry { key: 5, loc: 100 },
            SeedEntry { key: 7, loc: 200 },
            SeedEntry { key: 3, loc: 1300 },
        ];

        let result = sort_merge_join(&mut query, &mut refs, &[0, 100], &[0, 1000]);
        // key=3: 1 query × 1 ref = 1 pair
        // key=5: 2 query × 1 ref = 2 pairs
        // key=7: no query match
        assert_eq!(result.query_locs.len(), 3);
        assert_eq!(result.ref_locs.len(), 3);
    }
}
