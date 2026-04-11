use std::collections::HashMap;

use rayon::prelude::*;

use crate::basic::reduction::Reduction;
use crate::basic::seed::PackedSeed;
use crate::basic::shape::Shape;
use crate::basic::value::Letter;
use super::seed_match::{self, SeedLoc, SeedMatch};

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

/// Find seed matches in parallel across multiple queries.
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
}
