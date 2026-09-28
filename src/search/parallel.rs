use std::collections::HashMap;

use rayon::prelude::*;

use super::seed_array::{
    emit_seed_matches_filtered, match_blocks, sort_merge_join, sort_merge_seed_matches,
    sort_merge_seed_matches_with_complexity, SeedArray,
};
use super::seed_match::{self, SeedLoc, SeedMatch};
use crate::basic::reduction::Reduction;
use crate::basic::seed::PackedSeed;
use crate::basic::shape::Shape;
use crate::basic::value::Letter;
use crate::search::seed_complexity;

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
    find_seed_matches_partitioned_filtered_min_query_len_impl(
        query_seqs,
        ref_seqs,
        shape,
        reduction,
        complexity_cut,
        freq_sd,
        min_query_len,
        None,
        |matches| matches,
    )
}

fn build_search_seed_array(
    seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    seedp_bits: i32,
    min_query_len: usize,
    sketch_size: usize,
) -> SeedArray {
    if sketch_size > 0 {
        SeedArray::build_sketch_with_min_query_len(
            seqs,
            shape,
            reduction,
            seedp_bits,
            sketch_size,
            min_query_len,
        )
    } else {
        SeedArray::build_with_complexity_cut_and_min_query_len(
            seqs,
            shape,
            reduction,
            seedp_bits,
            0.0,
            min_query_len,
        )
    }
}

/// Collect the low-complexity query positions produced by the joined seed
/// groups without retaining the cross product.  Blastp uses this as a masking
/// prepass so its partition consumer can run the left-most filter immediately.
pub fn collect_low_complexity_positions_partitioned_min_query_len(
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    min_query_len: usize,
    sketch_size: usize,
) -> Vec<(u32, u32)> {
    if complexity_cut <= 0.0 {
        return Vec::new();
    }
    let seedp_bits = DEFAULT_SEEDP_BITS;
    let mut query_sa = build_search_seed_array(
        query_seqs,
        shape,
        reduction,
        seedp_bits,
        min_query_len,
        sketch_size,
    );
    let mut ref_sa =
        build_search_seed_array(ref_seqs, shape, reduction, seedp_bits, 0, sketch_size);
    let num_partitions = query_sa.num_partitions();
    let query_offsets = query_sa.seq_offsets().to_vec();

    fn split_into_partitions(
        sa: &mut SeedArray,
        num_partitions: usize,
    ) -> Vec<&mut [super::seed_array::SeedEntry]> {
        let offsets: Vec<usize> = (0..=num_partitions)
            .map(|partition| sa.partition_offset(partition))
            .collect();
        let mut remaining = sa.data_mut();
        let mut partitions = Vec::with_capacity(num_partitions);
        for partition in 0..num_partitions {
            let len = offsets[partition + 1] - offsets[partition];
            let (current, rest) = remaining.split_at_mut(len);
            partitions.push(current);
            remaining = rest;
        }
        partitions
    }

    let mut query_parts = split_into_partitions(&mut query_sa, num_partitions);
    let mut ref_parts = split_into_partitions(&mut ref_sa, num_partitions);
    let collect_partition =
        |query_part: &mut [super::seed_array::SeedEntry],
         ref_part: &mut [super::seed_array::SeedEntry]| {
            let blocks = match_blocks(query_part, ref_part);
            let mut masked = Vec::new();
            for block in blocks.blocks {
                let first = super::seed_array::decode_seq_pos(
                    &query_offsets,
                    query_part[block.q_start as usize],
                );
                let query = query_seqs[first.0 as usize];
                let pos = first.1 as usize;
                if pos >= query.len()
                    || !seed_complexity::seed_is_complex(
                        &query[pos..],
                        shape,
                        complexity_cut,
                        reduction,
                    )
                {
                    masked.extend((block.q_start..block.q_start + block.q_count).map(|index| {
                        super::seed_array::decode_seq_pos(
                            &query_offsets,
                            query_part[index as usize],
                        )
                    }));
                }
            }
            masked
        };
    let per_partition: Vec<Vec<(u32, u32)>> = if rayon::current_num_threads() == 1 {
        query_parts
            .iter_mut()
            .zip(ref_parts.iter_mut())
            .map(|(query_part, ref_part)| collect_partition(query_part, ref_part))
            .collect()
    } else {
        query_parts
            .par_iter_mut()
            .zip(ref_parts.par_iter_mut())
            .map(|(query_part, ref_part)| collect_partition(query_part, ref_part))
            .collect()
    };
    per_partition.into_iter().flatten().collect()
}

/// Join and apply the stage-1 Hamming filter inside each partition worker.
/// Returns filtered matches plus the number of raw joined pairs examined.
/// Keeping raw cross-products partition-local bounds their lifetime and avoids
/// materializing a shape-wide raw `Vec<SeedMatch>`.
pub fn find_seed_matches_partitioned_filtered_hamming_min_query_len(
    query_seqs: &[&[Letter]],
    query_restores: &[Vec<(usize, Letter)>],
    ref_seqs: &[&[Letter]],
    ref_restores: &[Vec<(usize, Letter)>],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    freq_sd: f64,
    min_query_len: usize,
    hamming_filter_id: u32,
) -> (Vec<SeedMatch>, usize) {
    let (matches, raw_count, _) =
        find_seed_matches_partitioned_filtered_hamming_min_query_len_with_masked_positions(
            query_seqs,
            query_restores,
            ref_seqs,
            ref_restores,
            shape,
            reduction,
            complexity_cut,
            freq_sd,
            min_query_len,
            hamming_filter_id,
        );
    (matches, raw_count)
}

/// Stage-1 join plus the query positions that C++ marks with `SEED_MASK` when
/// a joined seed group fails the low-complexity filter. Those bits are consumed
/// later by `left_most_filter`; merely dropping the group loses that state.
pub fn find_seed_matches_partitioned_filtered_hamming_min_query_len_with_masked_positions(
    query_seqs: &[&[Letter]],
    query_restores: &[Vec<(usize, Letter)>],
    ref_seqs: &[&[Letter]],
    ref_restores: &[Vec<(usize, Letter)>],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    freq_sd: f64,
    min_query_len: usize,
    hamming_filter_id: u32,
) -> (Vec<SeedMatch>, usize, Vec<(u32, u32)>) {
    let (matches, raw_count, masked_positions) =
        map_seed_matches_partitioned_filtered_hamming_min_query_len_with_masked_positions(
            query_seqs,
            query_restores,
            ref_seqs,
            ref_restores,
            shape,
            reduction,
            complexity_cut,
            freq_sd,
            min_query_len,
            hamming_filter_id,
            |matches| matches,
        );
    (matches, raw_count, masked_positions)
}

/// Run a partition-local consumer immediately after stage-1 filtering.
///
/// Unlike returning every Hamming survivor, this permits stage 2 to discard
/// candidates before the next partition is joined.  The peak lifetime is then
/// bounded by one joined partition per worker instead of all seed pairs across
/// every shape.
pub fn map_seed_matches_partitioned_filtered_hamming_min_query_len_with_masked_positions<T, F>(
    query_seqs: &[&[Letter]],
    query_restores: &[Vec<(usize, Letter)>],
    ref_seqs: &[&[Letter]],
    ref_restores: &[Vec<(usize, Letter)>],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    freq_sd: f64,
    min_query_len: usize,
    hamming_filter_id: u32,
    map_partition: F,
) -> (Vec<T>, usize, Vec<(u32, u32)>)
where
    T: Send,
    F: Fn(Vec<SeedMatch>) -> Vec<T> + Sync,
{
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Mutex;

    let raw_count = AtomicUsize::new(0);
    let masked_positions = Mutex::new(Vec::new());
    let matches = find_seed_matches_partitioned_filtered_min_query_len_impl(
        query_seqs,
        ref_seqs,
        shape,
        reduction,
        complexity_cut,
        freq_sd,
        min_query_len,
        Some(&masked_positions),
        |mut matches| {
            raw_count.fetch_add(matches.len(), Ordering::Relaxed);
            crate::search::hamming_filter::retain_hamming_filter_sequence_set(
                &mut matches,
                query_seqs,
                query_restores,
                ref_seqs,
                ref_restores,
                hamming_filter_id,
            );
            map_partition(matches)
        },
    );
    (
        matches,
        raw_count.into_inner(),
        masked_positions.into_inner().unwrap(),
    )
}

/// Streaming stage-1 join. Matching-key cross products are evaluated through
/// the fingerprint kernel and delivered to `map_batch` in q-major batches of
/// at most 32; no `Vec<SeedMatch>` proportional to q_count*r_count is built.
pub fn map_seed_matches_partitioned_streaming_hamming_min_query_len<T, F>(
    query_seqs: &[&[Letter]],
    query_fingerprint_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    ref_fingerprint_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    min_query_len: usize,
    hamming_filter_id: u32,
    sketch_size: usize,
    map_batch: F,
) -> (Vec<Vec<T>>, usize)
where
    T: Send,
    F: Fn(&[SeedMatch], &mut Vec<T>) -> usize + Sync,
{
    let mut output = Vec::new();
    let raw_count = visit_seed_matches_partitioned_streaming_hamming_min_query_len(
        query_seqs,
        query_fingerprint_seqs,
        ref_seqs,
        ref_fingerprint_seqs,
        shape,
        reduction,
        complexity_cut,
        min_query_len,
        hamming_filter_id,
        sketch_size,
        usize::MAX,
        map_batch,
        |partitions| output.extend(partitions),
    );
    (output, raw_count)
}

/// Bounded-output variant of
/// [`map_seed_matches_partitioned_streaming_hamming_min_query_len`]. Seed
/// partitions are processed in ordered batches and handed to `consume` before
/// the next batch is built. This lets callers spill retained hits without
/// first materializing a complete shape's output.
#[allow(clippy::too_many_arguments)]
pub fn visit_seed_matches_partitioned_streaming_hamming_min_query_len<T, F, C>(
    query_seqs: &[&[Letter]],
    query_fingerprint_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    ref_fingerprint_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    min_query_len: usize,
    hamming_filter_id: u32,
    sketch_size: usize,
    partition_batch_size: usize,
    map_batch: F,
    consume: C,
) -> usize
where
    T: Send,
    F: Fn(&[SeedMatch], &mut Vec<T>) -> usize + Sync,
    C: FnMut(Vec<Vec<T>>) + Send,
{
    macro_rules! dispatch {
        ($kernel:expr) => {
            return visit_seed_matches_partitioned_streaming_hamming_min_query_len_for::<
                { $kernel },
                T,
                F,
                C,
            >(
                query_seqs,
                query_fingerprint_seqs,
                ref_seqs,
                ref_fingerprint_seqs,
                shape,
                reduction,
                complexity_cut,
                min_query_len,
                hamming_filter_id,
                sketch_size,
                partition_batch_size,
                map_batch,
                consume,
            );
        };
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx512bw") {
            dispatch!(crate::search::hamming_filter::FP_AVX512BW);
        }
        if std::arch::is_x86_feature_detected!("avx2") {
            dispatch!(crate::search::hamming_filter::FP_AVX2);
        }
        if std::arch::is_x86_feature_detected!("sse2") {
            dispatch!(crate::search::hamming_filter::FP_SSE2);
        }
    }
    #[cfg(target_arch = "aarch64")]
    dispatch!(crate::search::hamming_filter::FP_NEON);
    #[cfg(not(target_arch = "aarch64"))]
    dispatch!(crate::search::hamming_filter::FP_SCALAR);
}

#[allow(clippy::too_many_arguments)]
fn visit_seed_matches_partitioned_streaming_hamming_min_query_len_for<
    const HAMMING_KERNEL: u8,
    T,
    F,
    C,
>(
    query_seqs: &[&[Letter]],
    query_fingerprint_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    ref_fingerprint_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    min_query_len: usize,
    hamming_filter_id: u32,
    sketch_size: usize,
    partition_batch_size: usize,
    map_batch: F,
    mut consume: C,
) -> usize
where
    T: Send,
    F: Fn(&[SeedMatch], &mut Vec<T>) -> usize + Sync,
    C: FnMut(Vec<Vec<T>>) + Send,
{
    use std::sync::atomic::{AtomicUsize, Ordering};

    let seedp_bits = DEFAULT_SEEDP_BITS;
    let mut query_sa = build_search_seed_array(
        query_seqs,
        shape,
        reduction,
        seedp_bits,
        min_query_len,
        sketch_size,
    );
    let mut ref_sa =
        build_search_seed_array(ref_seqs, shape, reduction, seedp_bits, 0, sketch_size);
    let num_partitions = query_sa.num_partitions();
    let query_offsets = query_sa.seq_offsets().to_vec();
    let ref_offsets = ref_sa.seq_offsets().to_vec();

    fn split_into_partitions(
        sa: &mut SeedArray,
        num_partitions: usize,
    ) -> Vec<&mut [super::seed_array::SeedEntry]> {
        let offsets: Vec<usize> = (0..=num_partitions)
            .map(|partition| sa.partition_offset(partition))
            .collect();
        let mut remaining = sa.data_mut();
        let mut partitions = Vec::with_capacity(num_partitions);
        for partition in 0..num_partitions {
            let len = offsets[partition + 1] - offsets[partition];
            let (current, rest) = remaining.split_at_mut(len);
            partitions.push(current);
            remaining = rest;
        }
        partitions
    }

    let raw_count = AtomicUsize::new(0);
    // Aggregate once per partition, rather than once per hit, so exposing the
    // same stage boundary as upstream's TENTATIVE_MATCHES1 counter has
    // negligible cost in the hot cross-product loop.
    let hamming_count = AtomicUsize::new(0);
    let ungapped_count = AtomicUsize::new(0);
    let mut query_parts = split_into_partitions(&mut query_sa, num_partitions);
    let mut ref_parts = split_into_partitions(&mut ref_sa, num_partitions);

    let process = |partition: usize,
                   query_part: &mut [super::seed_array::SeedEntry],
                   ref_part: &mut [super::seed_array::SeedEntry]| {
        let blocks = match_blocks(query_part, ref_part);
        let mut output = Vec::new();
        let mut partition_hamming_count = 0usize;
        let mut partition_ungapped_count = 0usize;
        // C++ keeps the decoded locations in its per-worker WorkSet. Reuse
        // these two buffers across joined seed groups in the partition rather
        // than allocating a pair of Vecs for every group.
        let mut query_locs = Vec::new();
        let mut target_locs = Vec::new();
        let mut target_subjects = Vec::new();
        for block in blocks.blocks {
            let first_query = super::seed_array::decode_seq_pos(
                &query_offsets,
                query_part[block.q_start as usize],
            );
            if complexity_cut > 0.0 {
                let query = query_seqs[first_query.0 as usize];
                let pos = first_query.1 as usize;
                if pos >= query.len()
                    || !seed_complexity::seed_is_complex(
                        &query[pos..],
                        shape,
                        complexity_cut,
                        reduction,
                    )
                {
                    continue;
                }
            }
            raw_count.fetch_add(
                block.q_count as usize * block.r_count as usize,
                Ordering::Relaxed,
            );
            query_locs.clear();
            query_locs.extend((block.q_start..block.q_start + block.q_count).map(|index| {
                super::seed_array::decode_seq_pos(&query_offsets, query_part[index as usize])
            }));
            target_locs.clear();
            target_subjects.clear();
            for index in block.r_start..block.r_start + block.r_count {
                let entry = ref_part[index as usize];
                let target = super::seed_array::decode_seq_pos(&ref_offsets, entry);
                // SeedEntry locations omit SequenceSet's perimeter and record
                // delimiters. Convert once per distinct target, then reuse the
                // absolute backing position across the whole query cross
                // product just as upstream's PackedLoc does.
                let subject = entry.loc as u64
                    + target.0 as u64
                    + crate::data::sequence_set::LetterStringSet::PERIMETER_PADDING as u64;
                target_locs.push(target);
                target_subjects.push(subject);
            }
            let mut batch = Vec::with_capacity(32);
            let mut last_query = None;
            crate::search::hamming_filter::visit_hamming_group_for::<HAMMING_KERNEL, _>(
                &query_locs,
                &target_locs,
                query_fingerprint_seqs,
                ref_fingerprint_seqs,
                hamming_filter_id,
                |query, target, target_index| {
                    partition_hamming_count += 1;
                    if last_query.is_some() && last_query != Some(query) && !batch.is_empty() {
                        partition_ungapped_count += map_batch(&batch, &mut output);
                        batch.clear();
                    }
                    last_query = Some(query);
                    // This streaming-only field packs the absolute subject
                    // position above the partition bits. The downstream
                    // left-most filter needs both but never needs the original
                    // seed key; retaining PackedSeed's layout avoids growing
                    // the 32-byte transient match record.
                    let seed =
                        (target_subjects[target_index] << seedp_bits) | partition as PackedSeed;
                    batch.push(SeedMatch {
                        query_id: query.0,
                        query_pos: query.1,
                        ref_id: target.0,
                        ref_pos: target.1,
                        seed,
                        shape_id: 0,
                    });
                    if batch.len() == 32 {
                        partition_ungapped_count += map_batch(&batch, &mut output);
                        batch.clear();
                    }
                },
            );
            if !batch.is_empty() {
                partition_ungapped_count += map_batch(&batch, &mut output);
            }
        }
        hamming_count.fetch_add(partition_hamming_count, Ordering::Relaxed);
        ungapped_count.fetch_add(partition_ungapped_count, Ordering::Relaxed);
        output
    };

    let batch_size = partition_batch_size.max(1).min(num_partitions.max(1));
    if rayon::current_num_threads() == 1 {
        for begin in (0..num_partitions).step_by(batch_size) {
            let end = (begin + batch_size).min(num_partitions);
            let output = query_parts[begin..end]
                .iter_mut()
                .zip(ref_parts[begin..end].iter_mut())
                .enumerate()
                .map(|(offset, (query_part, ref_part))| {
                    process(begin + offset, query_part, ref_part)
                })
                .collect();
            consume(output);
        }
    } else {
        // Maintain one strictly bounded sliding window across all seed
        // partitions. Fixed waves leave every worker waiting for each wave's
        // slowest prefix; here an ordered completed prefix is consumed and
        // replaced immediately. The coordinator yields into Rayon while it
        // waits, so all configured workers remain available for partition
        // jobs. We consume ready outputs before refilling their slots, keeping
        // active jobs + buffered outputs at or below `batch_size`.
        let (sender, receiver) = std::sync::mpsc::channel();
        let receiver = std::sync::Mutex::new(receiver);
        rayon::scope(|scope| {
            let mut pending = query_parts.iter_mut().zip(ref_parts.iter_mut()).enumerate();
            let mut active = 0usize;
            for _ in 0..batch_size {
                let Some((partition, (query_part, ref_part))) = pending.next() else {
                    break;
                };
                let sender = sender.clone();
                let process = &process;
                scope.spawn(move |_| {
                    let output = process(partition, query_part, ref_part);
                    sender
                        .send((partition, output))
                        .expect("stage-1 result receiver dropped");
                });
                active += 1;
            }

            let mut next_partition = 0usize;
            let mut completed = std::collections::BTreeMap::new();
            while active != 0 {
                let (partition, output) = loop {
                    let result = receiver
                        .lock()
                        .expect("stage-1 result receiver poisoned")
                        .try_recv();
                    match result {
                        Ok(result) => break result,
                        Err(std::sync::mpsc::TryRecvError::Empty) => {
                            rayon::yield_now();
                        }
                        Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                            panic!("stage-1 partition worker stopped without a result")
                        }
                    }
                };
                active -= 1;
                completed.insert(partition, output);
                let mut ready = Vec::new();
                while let Some(output) = completed.remove(&next_partition) {
                    ready.push(output);
                    next_partition += 1;
                }
                let refill = ready.len();
                if refill != 0 {
                    consume(ready);
                }
                for _ in 0..refill {
                    let Some((partition, (query_part, ref_part))) = pending.next() else {
                        break;
                    };
                    let sender = sender.clone();
                    let process = &process;
                    scope.spawn(move |_| {
                        let output = process(partition, query_part, ref_part);
                        sender
                            .send((partition, output))
                            .expect("stage-1 result receiver dropped");
                    });
                    active += 1;
                }
            }
        });
    }
    eprintln!(
        "Stage-1 filters: {} -> {} (Hamming) -> {} (ungapped)",
        raw_count.load(Ordering::Relaxed),
        hamming_count.load(Ordering::Relaxed),
        ungapped_count.load(Ordering::Relaxed),
    );
    raw_count.into_inner()
}

fn find_seed_matches_partitioned_filtered_min_query_len_impl<F, T>(
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
    complexity_cut: f64,
    freq_sd: f64,
    min_query_len: usize,
    masked_positions: Option<&std::sync::Mutex<Vec<(u32, u32)>>>,
    process_partition: F,
) -> Vec<T>
where
    T: Send,
    F: Fn(Vec<SeedMatch>) -> Vec<T> + Sync,
{
    let seedp_bits = DEFAULT_SEEDP_BITS;
    let single_thread = rayon::current_num_threads() == 1;

    // C++ default stage0 builds both seed arrays without low-complexity
    // filtering, computes the hash join, then calls `Search::mask_seeds` on
    // matched query seed groups. Doing complexity checks here for every DB seed
    // is both slower and less faithful.
    if single_thread && freq_sd <= 0.0 {
        const INDEX_CHUNKS: usize = 4;
        let num_partitions = crate::basic::seed::seedp_count(seedp_bits) as usize;
        let chunk_size = num_partitions.div_ceil(INDEX_CHUNKS);
        let mut matches = Vec::new();
        let query_counts = SeedArray::count_partitions(
            query_seqs,
            shape,
            reduction,
            seedp_bits,
            0.0,
            min_query_len,
        );
        let ref_counts =
            SeedArray::count_partitions(ref_seqs, shape, reduction, seedp_bits, 0.0, 0);
        let mut query_data = Vec::new();
        let mut ref_data = Vec::new();
        for chunk in 0..INDEX_CHUNKS {
            let partition_begin = chunk * chunk_size;
            let partition_end = ((chunk + 1) * chunk_size).min(num_partitions);
            if partition_begin >= partition_end {
                continue;
            }
            let mut query_sa = SeedArray::build_partition_range_from_counts_reuse(
                query_seqs,
                shape,
                reduction,
                seedp_bits,
                min_query_len,
                partition_begin,
                partition_end,
                &query_counts,
                query_data,
            );
            let mut ref_sa = SeedArray::build_partition_range_from_counts_reuse(
                ref_seqs,
                shape,
                reduction,
                seedp_bits,
                0,
                partition_begin,
                partition_end,
                &ref_counts,
                ref_data,
            );
            let query_offsets = query_sa.seq_offsets().to_vec();
            let ref_offsets = ref_sa.seq_offsets().to_vec();
            for p in partition_begin..partition_end {
                let q_part = query_sa.partition_mut(p as u32);
                let r_part = ref_sa.partition_mut(p as u32);
                let part = if complexity_cut > 0.0 {
                    sort_merge_seed_matches_with_complexity(
                        q_part,
                        r_part,
                        &query_offsets,
                        &ref_offsets,
                        p as u32,
                        seedp_bits,
                        |query_entries, (query_id, query_pos), _| {
                            let query = query_seqs[query_id as usize];
                            let pos = query_pos as usize;
                            let keep = pos < query.len()
                                && seed_complexity::seed_is_complex(
                                    &query[pos..],
                                    shape,
                                    complexity_cut,
                                    reduction,
                                );
                            if !keep {
                                if let Some(masked) = masked_positions {
                                    let mut masked = masked.lock().unwrap();
                                    masked.extend(query_entries.iter().map(|entry| {
                                        super::seed_array::decode_seq_pos(&query_offsets, *entry)
                                    }));
                                }
                            }
                            keep
                        },
                    )
                } else {
                    sort_merge_seed_matches(
                        q_part,
                        r_part,
                        &query_offsets,
                        &ref_offsets,
                        p as u32,
                        seedp_bits,
                    )
                };
                let mut part = process_partition(part);
                matches.append(&mut part);
            }
            query_data = query_sa.into_data();
            ref_data = ref_sa.into_data();
        }
        return matches;
    }

    let query_sa = SeedArray::build_with_complexity_cut_and_min_query_len(
        query_seqs,
        shape,
        reduction,
        seedp_bits,
        0.0,
        min_query_len,
    );
    let ref_sa = SeedArray::build_with_complexity_cut(ref_seqs, shape, reduction, seedp_bits, 0.0);

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
    let query_offsets = query_sa.seq_offsets().to_vec();
    let ref_offsets = ref_sa.seq_offsets().to_vec();
    let mut query_parts = split_into_partitions(&mut query_sa, num_partitions);
    let mut ref_parts = split_into_partitions(&mut ref_sa, num_partitions);

    // Phase 1: sort partitions, find match blocks, accumulate q/r count stats.
    if freq_sd <= 0.0 {
        let emit_part = |p: usize,
                         q_part: &mut [super::seed_array::SeedEntry],
                         r_part: &mut [super::seed_array::SeedEntry]| {
            if complexity_cut > 0.0 {
                sort_merge_seed_matches_with_complexity(
                    q_part,
                    r_part,
                    &query_offsets,
                    &ref_offsets,
                    p as u32,
                    seedp_bits,
                    |query_entries, (query_id, query_pos), _| {
                        let query = query_seqs[query_id as usize];
                        let pos = query_pos as usize;
                        let keep = pos < query.len()
                            && seed_complexity::seed_is_complex(
                                &query[pos..],
                                shape,
                                complexity_cut,
                                reduction,
                            );
                        if !keep {
                            if let Some(masked) = masked_positions {
                                let mut masked = masked.lock().unwrap();
                                masked.extend(query_entries.iter().map(|entry| {
                                    super::seed_array::decode_seq_pos(&query_offsets, *entry)
                                }));
                            }
                        }
                        keep
                    },
                )
            } else {
                sort_merge_seed_matches(
                    q_part,
                    r_part,
                    &query_offsets,
                    &ref_offsets,
                    p as u32,
                    seedp_bits,
                )
            }
        };
        return if single_thread {
            let mut matches = Vec::new();
            for (p, (q_part, r_part)) in
                query_parts.iter_mut().zip(ref_parts.iter_mut()).enumerate()
            {
                let mut part = process_partition(emit_part(p, q_part, r_part));
                matches.append(&mut part);
            }
            matches
        } else {
            query_parts
                .par_iter_mut()
                .zip(ref_parts.par_iter_mut())
                .enumerate()
                .flat_map_iter(|(p, (q_part, r_part))| {
                    process_partition(emit_part(p, q_part, r_part)).into_iter()
                })
                .collect()
        };
    }
    let blocks_per_part: Vec<super::seed_array::PartitionBlocks> = if single_thread {
        query_parts
            .iter_mut()
            .zip(ref_parts.iter_mut())
            .map(|(q_part, r_part)| match_blocks(q_part, r_part))
            .collect()
    } else {
        query_parts
            .par_iter_mut()
            .zip(ref_parts.par_iter_mut())
            .map(|(q_part, r_part)| match_blocks(q_part, r_part))
            .collect()
    };

    // Compute global q_max, r_max thresholds via SD reduction. C++
    // `FrequentSeeds::build` (frequent_seeds.cpp:107-116) computes one
    // per-partition `Sd` per worker thread and then constructs a global `Sd`
    // from the vector via `Sd::Sd(vector<Sd>&)`. The group-merge weights each
    // partition by its `k` (= n+1), which is part of the bit-identical
    // threshold computation — pairwise Welford merging would drift.
    let (q_max, r_max) = if freq_sd > 0.0 {
        let per_part_sds: Vec<(Sd, Sd)> = if single_thread {
            blocks_per_part
                .iter()
                .map(|pb| {
                    let mut qs = Sd::new();
                    let mut rs = Sd::new();
                    for b in &pb.blocks {
                        qs.add(b.q_count as f64);
                        rs.add(b.r_count as f64);
                    }
                    (qs, rs)
                })
                .collect()
        } else {
            blocks_per_part
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
                .collect()
        };
        let (q_groups, r_groups): (Vec<Sd>, Vec<Sd>) = per_part_sds.into_iter().unzip();
        let q_sd = Sd::merge_groups(&q_groups);
        let r_sd = Sd::merge_groups(&r_groups);
        // Match C++ `(unsigned)(mean + freq_sd*sd)` which is C-style truncation
        // toward zero — equivalent to `as u32` for non-negative values.
        (
            (q_sd.mean + freq_sd * q_sd.sd()).max(0.0) as u32,
            (r_sd.mean + freq_sd * r_sd.sd()).max(0.0) as u32,
        )
    } else {
        (u32::MAX, u32::MAX)
    };

    // Phase 2: emit (q_loc, r_loc) cross product per block, filtering blocks
    // whose counts exceed the frequent-seed threshold.
    let emit_part = |p: usize,
                     blocks: &super::seed_array::PartitionBlocks,
                     q_part: &&mut [super::seed_array::SeedEntry],
                     r_part: &&mut [super::seed_array::SeedEntry]| {
        // q_part / r_part are &&mut [SeedEntry] from the zip; the inner
        // ref is what `emit_seed_matches_filtered` wants.
        let q_part: &[super::seed_array::SeedEntry] = q_part;
        let r_part: &[super::seed_array::SeedEntry] = r_part;
        let _ = sort_merge_join; // keep the old function exported
        if complexity_cut > 0.0 {
            let blocks = super::seed_array::PartitionBlocks {
                blocks: blocks
                    .blocks
                    .iter()
                    .copied()
                    .filter(|b| {
                        let (_, (query_id, query_pos)) = {
                            let q_entry = q_part[b.q_start as usize];
                            (
                                q_entry,
                                super::seed_array::decode_seq_pos(&query_offsets, q_entry),
                            )
                        };
                        let query = query_seqs[query_id as usize];
                        let pos = query_pos as usize;
                        pos < query.len()
                            && seed_complexity::seed_is_complex(
                                &query[pos..],
                                shape,
                                complexity_cut,
                                reduction,
                            )
                    })
                    .collect(),
            };
            emit_seed_matches_filtered(
                q_part,
                r_part,
                &query_offsets,
                &ref_offsets,
                &blocks,
                q_max,
                r_max,
                p as u32,
                seedp_bits,
            )
        } else {
            emit_seed_matches_filtered(
                q_part,
                r_part,
                &query_offsets,
                &ref_offsets,
                blocks,
                q_max,
                r_max,
                p as u32,
                seedp_bits,
            )
        }
    };
    if single_thread {
        let mut matches = Vec::new();
        for (p, ((blocks, q_part), r_part)) in blocks_per_part
            .iter()
            .zip(&query_parts)
            .zip(&ref_parts)
            .enumerate()
        {
            let mut part = process_partition(emit_part(p, blocks, q_part, r_part));
            matches.append(&mut part);
        }
        matches
    } else {
        blocks_per_part
            .par_iter()
            .zip(query_parts.par_iter())
            .zip(ref_parts.par_iter())
            .enumerate()
            .flat_map_iter(|(p, ((blocks, q_part), r_part))| {
                process_partition(emit_part(p, blocks, q_part, r_part)).into_iter()
            })
            .collect()
    }
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
                            shape_id: 0,
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

    #[test]
    fn parallel_partition_flatten_preserves_serial_emission_order() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);
        let query_data: Vec<Vec<Letter>> = (0..8)
            .map(|shift| (0..96).map(|i| ((i + shift) % 20) as Letter).collect())
            .collect();
        let ref_data: Vec<Vec<Letter>> = (0..9)
            .map(|shift| (0..104).map(|i| ((i + shift) % 20) as Letter).collect())
            .collect();
        let query_seqs: Vec<&[Letter]> = query_data.iter().map(Vec::as_slice).collect();
        let ref_seqs: Vec<&[Letter]> = ref_data.iter().map(Vec::as_slice).collect();
        let run = || {
            find_seed_matches_partitioned_filtered_min_query_len(
                &query_seqs,
                &ref_seqs,
                &shape,
                &reduction,
                0.0,
                0.0,
                0,
            )
        };

        let serial = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(run);
        let parallel = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(run);

        let signature = |matches: &[SeedMatch]| {
            matches
                .iter()
                .map(|m| (m.query_id, m.query_pos, m.ref_id, m.ref_pos, m.seed))
                .collect::<Vec<_>>()
        };
        assert_eq!(signature(&parallel), signature(&serial));
    }

    #[test]
    fn partition_local_hamming_matches_post_join_filter_exactly() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);
        let query_data: Vec<Vec<Letter>> = (0..5)
            .map(|shift| (0..72).map(|i| ((i + shift) % 20) as Letter).collect())
            .collect();
        let ref_data: Vec<Vec<Letter>> = (0..6)
            .map(|shift| (0..80).map(|i| ((i + shift) % 20) as Letter).collect())
            .collect();
        let query_seqs: Vec<&[Letter]> = query_data.iter().map(Vec::as_slice).collect();
        let ref_seqs: Vec<&[Letter]> = ref_data.iter().map(Vec::as_slice).collect();

        let raw = find_seed_matches_partitioned_filtered_min_query_len(
            &query_seqs,
            &ref_seqs,
            &shape,
            &reduction,
            0.0,
            0.0,
            0,
        );
        let query_restores = vec![Vec::new(); query_seqs.len()];
        let ref_restores = vec![Vec::new(); ref_seqs.len()];
        let mut expected = raw.clone();
        crate::search::hamming_filter::retain_hamming_filter_sequence_set(
            &mut expected,
            &query_seqs,
            &query_restores,
            &ref_seqs,
            &ref_restores,
            40,
        );
        let (fused, raw_count) = find_seed_matches_partitioned_filtered_hamming_min_query_len(
            &query_seqs,
            &query_restores,
            &ref_seqs,
            &ref_restores,
            &shape,
            &reduction,
            0.0,
            0.0,
            0,
            40,
        );

        let signature = |matches: &[SeedMatch]| {
            matches
                .iter()
                .map(|m| (m.query_id, m.query_pos, m.ref_id, m.ref_pos, m.seed))
                .collect::<Vec<_>>()
        };
        assert_eq!(raw_count, raw.len());
        assert_eq!(signature(&fused), signature(&expected));
    }

    #[test]
    fn joined_low_complexity_groups_return_all_query_marker_positions() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);
        let query_data = [vec![0; 12], vec![0; 12]];
        let ref_data = [vec![0; 12]];
        let query_seqs: Vec<&[Letter]> = query_data.iter().map(Vec::as_slice).collect();
        let ref_seqs: Vec<&[Letter]> = ref_data.iter().map(Vec::as_slice).collect();
        let query_restores = vec![Vec::new(); query_seqs.len()];
        let ref_restores = vec![Vec::new(); ref_seqs.len()];
        let (matches, _, mut masked) =
            find_seed_matches_partitioned_filtered_hamming_min_query_len_with_masked_positions(
                &query_seqs,
                &query_restores,
                &ref_seqs,
                &ref_restores,
                &shape,
                &reduction,
                1.0,
                0.0,
                0,
                0,
            );
        masked.sort_unstable();
        assert!(matches.is_empty());
        assert_eq!(masked.len(), 18);
        assert_eq!(masked.first(), Some(&(0, 0)));
        assert_eq!(masked.last(), Some(&(1, 8)));
    }
}
