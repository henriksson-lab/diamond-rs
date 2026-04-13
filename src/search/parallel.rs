use std::collections::HashMap;

use rayon::prelude::*;

use crate::basic::reduction::Reduction;
use crate::basic::seed::PackedSeed;
use crate::basic::shape::Shape;
use crate::basic::value::Letter;
use super::seed_array::{SeedArray, sort_merge_join};
use super::seed_match::{self, SeedLoc, SeedMatch};

/// Default number of seed partition bits.
/// 4 bits = 16 partitions, matching C++ DIAMOND default.
const DEFAULT_SEEDP_BITS: i32 = 4;

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
    let seedp_bits = DEFAULT_SEEDP_BITS;

    // Build partitioned seed arrays
    let mut query_sa = SeedArray::build(query_seqs, shape, reduction, seedp_bits);
    let mut ref_sa = SeedArray::build(ref_seqs, shape, reduction, seedp_bits);

    let num_partitions = query_sa.num_partitions();

    // Collect partition boundaries so we can process in parallel.
    // We need to split the mutable borrows per partition.
    // Extract partition slices into separate owned vecs for parallel processing.
    let mut query_parts: Vec<Vec<super::seed_array::SeedEntry>> = (0..num_partitions as u32)
        .map(|p| query_sa.partition(p).to_vec())
        .collect();
    let mut ref_parts: Vec<Vec<super::seed_array::SeedEntry>> = (0..num_partitions as u32)
        .map(|p| ref_sa.partition(p).to_vec())
        .collect();

    // Join each partition in parallel
    let partition_results: Vec<Vec<SeedMatch>> = query_parts
        .par_iter_mut()
        .zip(ref_parts.par_iter_mut())
        .map(|(q_part, r_part)| {
            let join = sort_merge_join(q_part, r_part);
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
        let par_matches =
            find_seed_matches_parallel(&seqs, &par_index, &shape, &reduction);

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
}
