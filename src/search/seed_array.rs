//! Partitioned seed array matching C++ DIAMOND's SeedArray.
//!
//! Seeds are partitioned by their low bits (seedp_bits) into buckets.
//! Within each partition, entries are stored contiguously in a flat array.
//! This gives cache-friendly access during the join phase and enables
//! per-partition parallelism.

use crate::basic::reduction::Reduction;
use crate::basic::seed::{self, SeedPartition};
use crate::basic::shape::Shape;
use crate::basic::value::Letter;

/// A single entry in the seed array: the seed key (with partition bits stripped)
/// and the location it came from.
#[derive(Clone, Copy, Debug)]
pub struct SeedEntry {
    /// Seed key with low seedp_bits removed (the "offset" part).
    pub key: u32,
    /// Sequence index.
    pub seq_id: u32,
    /// Position within the sequence.
    pub pos: u32,
}

/// Partitioned seed array. Seeds are grouped by `seed & seedp_mask` into partitions.
/// Within each partition, entries are stored contiguously.
pub struct SeedArray {
    /// Flat storage of all entries, grouped by partition.
    data: Vec<SeedEntry>,
    /// offsets[p] = start index of partition p in data. Length = num_partitions + 1.
    offsets: Vec<usize>,
    /// Number of partition bits.
    pub seedp_bits: i32,
}

impl SeedArray {
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
        use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

        let num_partitions = seed::seedp_count(seedp_bits) as usize;
        let mask = seed::seedp_mask(seedp_bits);
        let shape_length = shape.length as usize;

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
                        if let Some(s) =
                            shape.set_seed_with_complexity(window, reduction, complexity_cut)
                        {
                            let p = seed::seed_partition(s, mask) as usize;
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
        let mut data: Vec<SeedEntry> = Vec::with_capacity(total);
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
                if let Some(s) = shape.set_seed_with_complexity(window, reduction, complexity_cut) {
                    let p = seed::seed_partition(s, mask) as usize;
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
                                seq_id: seq_id as u32,
                                pos: pos as u32,
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
    // Stable sort: C++ uses radix sort (`diamond/src/util/algo/radix_sort.h`),
    // which preserves input order for equal keys. Tied keys produce the same
    // cross-product order in C++; unstable sort here would shuffle them and
    // change downstream emission order — affecting tie-breaks in the
    // gapped-score ranking and culling stages.
    query_part.sort_by_key(|e| e.key);
    ref_part.sort_by_key(|e| e.key);
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
                query_locs.push((query_part[q as usize].seq_id, query_part[q as usize].pos));
                ref_locs.push((ref_part[r as usize].seq_id, ref_part[r as usize].pos));
            }
        }
    }
    JoinResult {
        query_locs,
        ref_locs,
    }
}

/// Sort-merge join on one partition.
///
/// Sort both sides by key, then walk through matching keys.
/// For each matching key, emit the cross product of query and ref locations.
pub fn sort_merge_join(query_part: &mut [SeedEntry], ref_part: &mut [SeedEntry]) -> JoinResult {
    if query_part.is_empty() || ref_part.is_empty() {
        return JoinResult {
            query_locs: Vec::new(),
            ref_locs: Vec::new(),
        };
    }

    // Sort both by key. Stable sort: see `match_blocks` above — C++ radix
    // sort preserves input order for equal keys, and downstream culling
    // tie-breaks depend on emission order.
    query_part.sort_by_key(|e| e.key);
    ref_part.sort_by_key(|e| e.key);

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
                    query_locs.push((query_part[q].seq_id, query_part[q].pos));
                    ref_locs.push((ref_part[r].seq_id, ref_part[r].pos));
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
                actual.push((p, entry.key, entry.seq_id, entry.pos));
            }
        }
        actual.sort_unstable();

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_sort_merge_join() {
        // Create two small partitions with some matching keys
        let mut query = vec![
            SeedEntry {
                key: 5,
                seq_id: 0,
                pos: 10,
            },
            SeedEntry {
                key: 3,
                seq_id: 0,
                pos: 20,
            },
            SeedEntry {
                key: 5,
                seq_id: 1,
                pos: 30,
            },
        ];
        let mut refs = vec![
            SeedEntry {
                key: 5,
                seq_id: 0,
                pos: 100,
            },
            SeedEntry {
                key: 7,
                seq_id: 0,
                pos: 200,
            },
            SeedEntry {
                key: 3,
                seq_id: 1,
                pos: 300,
            },
        ];

        let result = sort_merge_join(&mut query, &mut refs);
        // key=3: 1 query × 1 ref = 1 pair
        // key=5: 2 query × 1 ref = 2 pairs
        // key=7: no query match
        assert_eq!(result.query_locs.len(), 3);
        assert_eq!(result.ref_locs.len(), 3);
    }
}
