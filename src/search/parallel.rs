use std::collections::HashMap;

use rayon::prelude::*;

use super::seed_array::{emit_join_filtered, match_blocks, sort_merge_join, SeedArray};
use super::seed_match::{self, SeedLoc, SeedMatch};
use crate::basic::reduction::Reduction;
use crate::basic::seed::PackedSeed;
use crate::basic::shape::Shape;
use crate::basic::value::Letter;

/// Default number of seed partition bits.
/// 10 bits = 1024 partitions (memory `project_seed_array_partitioning`:
/// empirically ~30% faster than 16 partitions on sprot-scale blastp because
/// each partition's join arrays fit in L2). C++ DIAMOND's stage0 always
/// runs with `seedp_bits = SeedPartitionRange::BITS = 10` too — earlier
/// docstring claiming "16 partitions" was stale.
const DEFAULT_SEEDP_BITS: i32 = 10;

/// Build a seed index from sequences in parallel using rayon.
pub fn build_seed_index_parallel(
    seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
) -> HashMap<PackedSeed, Vec<SeedLoc>> {
    // Extract seeds from each sequence in parallel
    let per_seq_seeds: Vec<Vec<(PackedSeed, u32, u32)>> = seqs
        .par_iter()
        .enumerate()
        .map(|(seq_id, seq)| {
            seed_match::extract_seeds(seq, shape, reduction)
                .into_iter()
                .map(|(seed, pos)| (seed, seq_id as u32, pos))
                .collect()
        })
        .collect();

    // Merge into a single index
    let mut index: HashMap<PackedSeed, Vec<SeedLoc>> = HashMap::new();
    for seeds in per_seq_seeds {
        for (seed, seq_id, pos) in seeds {
            index.entry(seed).or_default().push(SeedLoc { seq_id, pos });
        }
    }
    index
}

/// Find seed matches using partitioned seed arrays and per-partition join.
///
/// This replaces the naive HashMap approach with DIAMOND's strategy:
/// 1. Build partitioned SeedArrays for both query and reference
/// 2. For each partition, sort-merge join to find matching seeds
/// 3. Collect results as SeedMatch list
///
/// Much faster than HashMap for large datasets due to:
/// - Pre-allocated contiguous memory (no hash table overhead)
/// - Per-partition parallelism via rayon
/// - Cache-friendly sequential access within each partition
pub fn find_seed_matches_partitioned(
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
) -> Vec<SeedMatch> {
    find_seed_matches_partitioned_with_complexity(query_seqs, ref_seqs, shape, reduction, 0.0)
}

/// Variant that mirrors C++ DIAMOND by skipping seeds whose unreduced
/// composition entropy is below `complexity_cut` during index construction.
pub fn find_seed_matches_partitioned_with_complexity(
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
) -> Vec<SeedMatch> {
    find_seed_matches_partitioned_filtered(
        query_seqs,
        ref_seqs,
        shape,
        reduction,
        complexity_cut,
        0.0,
    )
}

/// Running mean+variance accumulator, mirroring C++ `Sd` in
/// `diamond/src/util/util.h` BYTE-for-BYTE so the `freq_sd` cap matches.
/// C++ keeps `k = n + 1` (initialized to 1, post-increment in `add`) and
/// returns `sqrt(Q / (k - 1))` — i.e. **population** stddev, NOT the sample
/// stddev. Using sample stddev shifts the `mean + freq_sd * sd` threshold
/// and changes which seeds the frequent-seed filter drops.
///
/// The group merge in `Sd::Sd(vector<Sd>&)` also has an idiosyncratic
/// weighting: each group is weighted by its `k` (= n_i + 1), not `n_i`.
/// We replicate that exactly so the cap mirrors C++.
pub(super) struct Sd {
    /// Running mean (`A` in C++).
    pub(super) mean: f64,
    /// Sum of squared deviations (`Q` in C++).
    pub(super) q: f64,
    /// C++'s `k` = sample count + 1 (starts at 1, post-incremented in `add`).
    pub(super) k: f64,
}
impl Sd {
    pub(super) fn new() -> Self {
        Self {
            mean: 0.0,
            q: 0.0,
            k: 1.0,
        }
    }
    pub(super) fn add(&mut self, x: f64) {
        // Welford recurrence with the C++ pre-increment convention:
        //   Q += (k - 1) / k * d * d
        //   A += d / k
        //   ++k
        let d = x - self.mean;
        self.q += (self.k - 1.0) / self.k * d * d;
        self.mean += d / self.k;
        self.k += 1.0;
    }
    pub(super) fn merge_groups(groups: &[Sd]) -> Sd {
        // Direct port of C++ Sd::Sd(const vector<Sd>&) (util.cpp).
        // Each group is weighted by its `k`, including the empty `k = 1`
        // bias — that bias is part of the C++ behavior we mirror.
        let mut k = 0.0f64;
        let mut a = 0.0f64;
        let mut q = 0.0f64;
        for g in groups {
            k += g.k;
            a += g.mean * g.k;
            q += g.q;
        }
        if k > 0.0 {
            a /= k;
        }
        for g in groups {
            let d = g.mean - a;
            q += d * d * g.k;
        }
        Sd { mean: a, q, k }
    }
    pub(super) fn sd(&self) -> f64 {
        if self.k <= 1.0 {
            0.0
        } else {
            (self.q / (self.k - 1.0)).sqrt()
        }
    }
}

/// Full pipeline: build seed arrays with complexity filtering, then apply
/// `freq_sd`-based frequent-seed masking (matches C++ FrequentSeeds::build).
/// When `freq_sd > 0`, drop seed keys whose query- or ref-side occurrence
/// counts exceed `mean + freq_sd * stddev` across all matching keys.
pub fn find_seed_matches_partitioned_filtered(
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    freq_sd: f64,
) -> Vec<SeedMatch> {
    find_seed_matches_partitioned_filtered_min_query_len(
        query_seqs,
        ref_seqs,
        shape,
        reduction,
        complexity_cut,
        freq_sd,
        0,
    )
}

pub fn find_seed_matches_partitioned_filtered_min_query_len(
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    freq_sd: f64,
    min_query_len: usize,
) -> Vec<SeedMatch> {
    let seedp_bits = DEFAULT_SEEDP_BITS;

    // Build partitioned seed arrays (low-complexity seeds filtered on both sides
    // when complexity_cut > 0).
    let query_sa = SeedArray::build_with_complexity_cut_and_min_query_len(
        query_seqs,
        shape,
        reduction,
        seedp_bits,
        complexity_cut,
        min_query_len,
    );
    let ref_sa = SeedArray::build_with_complexity_cut(
        ref_seqs,
        shape,
        reduction,
        seedp_bits,
        complexity_cut,
    );

    let num_partitions = query_sa.num_partitions();

    // Split the seed-array `data` into disjoint mutable slices, one per
    // partition. `split_at_mut` repeatedly walks the buffer without cloning,
    // saving the ~2.3GB allocation + copy that the previous
    // `partition(p).to_vec()` loop incurred on a sprot-scale DB. The borrow
    // checker is satisfied because each call to `split_at_mut` produces two
    // non-overlapping &mut [SeedEntry]s.
    fn split_into_partitions(
        sa: &mut super::seed_array::SeedArray,
        num_partitions: usize,
    ) -> Vec<&mut [super::seed_array::SeedEntry]> {
        let offsets: Vec<usize> = (0..=num_partitions)
            .map(|p| sa.partition_offset(p))
            .collect();
        let mut data_slice: &mut [super::seed_array::SeedEntry] = sa.data_mut();
        let mut out: Vec<&mut [super::seed_array::SeedEntry]> = Vec::with_capacity(num_partitions);
        for p in 0..num_partitions {
            let len = offsets[p + 1] - offsets[p];
            let (left, right) = data_slice.split_at_mut(len);
            out.push(left);
            data_slice = right;
        }
        out
    }
    let mut query_sa = query_sa;
    let mut ref_sa = ref_sa;
    let mut query_parts = split_into_partitions(&mut query_sa, num_partitions);
    let mut ref_parts = split_into_partitions(&mut ref_sa, num_partitions);

    // Phase 1: sort partitions, find match blocks, accumulate q/r count stats.
    let blocks_per_part: Vec<super::seed_array::PartitionBlocks> = query_parts
        .par_iter_mut()
        .zip(ref_parts.par_iter_mut())
        .map(|(q_part, r_part)| match_blocks(q_part, r_part))
        .collect();

    // Compute global q_max, r_max thresholds via parallel SD reduction.
    // C++ `FrequentSeeds::build` (frequent_seeds.cpp:107-116) computes one
    // per-partition `Sd` per worker thread and then constructs a global `Sd`
    // from the vector via `Sd::Sd(vector<Sd>&)`. The group-merge weights each
    // partition by its `k` (= n+1), which is part of the bit-identical
    // threshold computation — pairwise Welford merging would drift.
    let per_part_sds: Vec<(Sd, Sd)> = blocks_per_part
        .par_iter()
        .map(|pb| {
            let mut qs = Sd::new();
            let mut rs = Sd::new();
            for b in &pb.blocks {
                qs.add(b.q_count as f64);
                rs.add(b.r_count as f64);
            }
            (qs, rs)
        })
        .collect();
    let (q_groups, r_groups): (Vec<Sd>, Vec<Sd>) = per_part_sds.into_iter().unzip();
    let q_sd = Sd::merge_groups(&q_groups);
    let r_sd = Sd::merge_groups(&r_groups);
    // Match C++ `(unsigned)(mean + freq_sd*sd)` which is C-style truncation
    // toward zero — equivalent to `as u32` for non-negative values.
    let q_max = if freq_sd > 0.0 {
        (q_sd.mean + freq_sd * q_sd.sd()).max(0.0) as u32
    } else {
        u32::MAX
    };
    let r_max = if freq_sd > 0.0 {
        (r_sd.mean + freq_sd * r_sd.sd()).max(0.0) as u32
    } else {
        u32::MAX
    };

    // Phase 2: emit (q_loc, r_loc) cross product per block, filtering blocks
    // whose counts exceed the frequent-seed threshold.
    let partition_results: Vec<Vec<SeedMatch>> = blocks_per_part
        .par_iter()
        .zip(query_parts.par_iter())
        .zip(ref_parts.par_iter())
        .map(|((blocks, q_part), r_part)| {
            // q_part / r_part are &&mut [SeedEntry] from the zip; the inner
            // ref is what `emit_join_filtered` wants.
            let q_part: &[super::seed_array::SeedEntry] = q_part;
            let r_part: &[super::seed_array::SeedEntry] = r_part;
            let join = emit_join_filtered(q_part, r_part, blocks, q_max, r_max);
            let _ = sort_merge_join; // keep the old function exported
            join.query_locs
                .iter()
                .zip(join.ref_locs.iter())
                .map(|(&(q_seq, q_pos), &(r_seq, r_pos))| SeedMatch {
                    query_id: q_seq,
                    query_pos: q_pos,
                    ref_id: r_seq,
                    ref_pos: r_pos,
                    seed: 0, // seed value not needed downstream
                })
                .collect()
        })
        .collect();

    // Flatten
    partition_results.into_iter().flatten().collect()
}

/// Find seed matches in parallel across multiple queries (legacy HashMap path).
pub fn find_seed_matches_parallel(
    queries: &[&[Letter]],
    ref_index: &HashMap<PackedSeed, Vec<SeedLoc>>,
    shape: &Shape,
    reduction: &Reduction,
) -> Vec<SeedMatch> {
    queries
        .par_iter()
        .enumerate()
        .flat_map(|(query_id, query)| {
            let seeds = seed_match::extract_seeds(query, shape, reduction);
            let mut matches = Vec::new();
            for (seed, query_pos) in seeds {
                if let Some(ref_locs) = ref_index.get(&seed) {
                    for ref_loc in ref_locs {
                        matches.push(SeedMatch {
                            query_id: query_id as u32,
                            query_pos,
                            ref_id: ref_loc.seq_id,
                            ref_pos: ref_loc.pos,
                            seed,
                        });
                    }
                }
            }
            matches
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_seed_index() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);

        let seq1: Vec<Letter> = (0..20).map(|i| (i % 20) as Letter).collect();
        let seq2: Vec<Letter> = (0..20).map(|i| ((i + 5) % 20) as Letter).collect();
        let seqs: Vec<&[Letter]> = vec![&seq1, &seq2];

        let index = build_seed_index_parallel(&seqs, &shape, &reduction);
        assert!(!index.is_empty());
    }

    #[test]
    fn test_parallel_matches_same_as_serial() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);

        let seq: Vec<Letter> = (0..30).map(|i| (i % 20) as Letter).collect();
        let seqs: Vec<&[Letter]> = vec![&seq];

        // Serial
        let serial_index = seed_match::build_seed_index(&seqs, &shape, &reduction);
        let serial_matches =
            seed_match::find_seed_matches(&seqs, &serial_index, &shape, &reduction);

        // Parallel
        let par_index = build_seed_index_parallel(&seqs, &shape, &reduction);
        let par_matches = find_seed_matches_parallel(&seqs, &par_index, &shape, &reduction);

        assert_eq!(
            serial_matches.len(),
            par_matches.len(),
            "Parallel should find same number of matches as serial"
        );
    }

    #[test]
    fn test_partitioned_matches_same_as_hashmap() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);

        let ref_seq: Vec<Letter> = (0..30).map(|i| (i % 20) as Letter).collect();
        let query_seq: Vec<Letter> = (0..30).map(|i| (i % 20) as Letter).collect();
        let ref_seqs: Vec<&[Letter]> = vec![&ref_seq];
        let query_seqs: Vec<&[Letter]> = vec![&query_seq];

        // HashMap path
        let index = build_seed_index_parallel(&ref_seqs, &shape, &reduction);
        let hashmap_matches = find_seed_matches_parallel(&query_seqs, &index, &shape, &reduction);

        // Partitioned path
        let partitioned_matches =
            find_seed_matches_partitioned(&query_seqs, &ref_seqs, &shape, &reduction);

        assert_eq!(
            hashmap_matches.len(),
            partitioned_matches.len(),
            "Partitioned ({}) should find same matches as HashMap ({})",
            partitioned_matches.len(),
            hashmap_matches.len()
        );
    }

    #[test]
    fn test_partitioned_filtered_respects_min_query_len() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);

        let short_query: Vec<Letter> = vec![0, 1, 2, 3];
        let long_query: Vec<Letter> = (0..20).map(|i| (i % 20) as Letter).collect();
        let ref_seq = long_query.clone();
        let query_seqs: Vec<&[Letter]> = vec![&short_query, &long_query];
        let ref_seqs: Vec<&[Letter]> = vec![&ref_seq];

        let all = find_seed_matches_partitioned_filtered_min_query_len(
            &query_seqs,
            &ref_seqs,
            &shape,
            &reduction,
            0.0,
            0.0,
            0,
        );
        assert!(all.iter().any(|m| m.query_id == 0));

        let filtered = find_seed_matches_partitioned_filtered_min_query_len(
            &query_seqs,
            &ref_seqs,
            &shape,
            &reduction,
            0.0,
            0.0,
            10,
        );
        assert!(!filtered.iter().any(|m| m.query_id == 0));
        assert!(filtered.iter().any(|m| m.query_id == 1));
    }
}
