use std::collections::HashMap;

use crate::basic::reduction::Reduction;
use crate::basic::seed::PackedSeed;
use crate::basic::shape::Shape;
use crate::basic::value::{Letter, DELIMITER_LETTER};

/// A seed hit location in a sequence set.
#[derive(Debug, Clone, Copy)]
pub struct SeedLoc {
    /// Sequence index.
    pub seq_id: u32,
    /// Position within the sequence.
    pub pos: u32,
}

/// Extract all seeds from a sequence using a given shape.
///
/// Returns a list of (packed_seed, position) pairs for all valid seed positions.
pub fn extract_seeds(
    seq: &[Letter],
    shape: &Shape,
    reduction: &Reduction,
) -> Vec<(PackedSeed, u32)> {
    let mut seeds = Vec::new();
    let slen = seq.len() as i32;

    if slen < shape.length {
        return seeds;
    }

    for pos in 0..=(slen - shape.length) {
        let pos = pos as usize;
        // Check if the window contains a delimiter
        let window = &seq[pos..pos + shape.length as usize];
        if window.contains(&DELIMITER_LETTER) {
            continue;
        }

        if let Some(seed) = shape.set_seed(window, reduction) {
            seeds.push((seed, pos as u32));
        }
    }

    seeds
}

/// Build a seed index (hash table) from a set of sequences.
///
/// Returns a map from packed seed value to list of locations where that seed occurs.
pub fn build_seed_index(
    seqs: &[&[Letter]],
    shape: &Shape,
    reduction: &Reduction,
) -> HashMap<PackedSeed, Vec<SeedLoc>> {
    let mut index: HashMap<PackedSeed, Vec<SeedLoc>> = HashMap::new();

    for (seq_id, seq) in seqs.iter().enumerate() {
        let seeds = extract_seeds(seq, shape, reduction);
        for (seed, pos) in seeds {
            index
                .entry(seed)
                .or_default()
                .push(SeedLoc {
                    seq_id: seq_id as u32,
                    pos,
                });
        }
    }

    index
}

/// Find seed matches between query sequences and a reference seed index.
///
/// This is the core "hash join" operation that identifies potential alignments.
pub fn find_seed_matches(
    queries: &[&[Letter]],
    ref_index: &HashMap<PackedSeed, Vec<SeedLoc>>,
    shape: &Shape,
    reduction: &Reduction,
) -> Vec<SeedMatch> {
    let mut matches = Vec::new();

    for (query_id, query) in queries.iter().enumerate() {
        let seeds = extract_seeds(query, shape, reduction);
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
    }

    matches
}

/// A seed match between a query and reference position.
#[derive(Debug, Clone)]
pub struct SeedMatch {
    pub query_id: u32,
    pub query_pos: u32,
    pub ref_id: u32,
    pub ref_pos: u32,
    pub seed: PackedSeed,
}

impl SeedMatch {
    /// The diagonal of this seed match (query_pos - ref_pos).
    pub fn diagonal(&self) -> i64 {
        self.query_pos as i64 - self.ref_pos as i64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_reduction() -> Reduction {
        Reduction::default_reduction()
    }

    #[test]
    fn test_extract_seeds() {
        let r = default_reduction();
        let shape = Shape::from_code("111", &r);
        // Sequence: A R N D C (0 1 2 3 4)
        let seq = vec![0i8, 1, 2, 3, 4];
        let seeds = extract_seeds(&seq, &shape, &r);
        // Should get 3 seeds for positions 0,1,2
        assert_eq!(seeds.len(), 3);
    }

    #[test]
    fn test_extract_seeds_with_delimiter() {
        let r = default_reduction();
        let shape = Shape::from_code("111", &r);
        let seq = vec![0i8, DELIMITER_LETTER, 2, 3, 4];
        let seeds = extract_seeds(&seq, &shape, &r);
        // Should skip positions that overlap the delimiter
        assert!(seeds.len() < 3);
    }

    #[test]
    fn test_seed_index_and_match() {
        let r = default_reduction();
        let shape = Shape::from_code("11", &r);

        let ref_seq: Vec<Letter> = vec![0, 1, 2, 3, 0, 1];
        let query_seq: Vec<Letter> = vec![0, 1, 4, 5];

        let ref_seqs: Vec<&[Letter]> = vec![&ref_seq];
        let query_seqs: Vec<&[Letter]> = vec![&query_seq];

        let index = build_seed_index(&ref_seqs, &shape, &r);
        let matches = find_seed_matches(&query_seqs, &index, &shape, &r);

        // "01" appears in both query (pos 0) and reference (pos 0, 4)
        assert!(matches.len() >= 2, "Expected at least 2 matches, got {}", matches.len());
    }

    #[test]
    fn test_seed_match_diagonal() {
        let m = SeedMatch {
            query_id: 0,
            query_pos: 10,
            ref_id: 0,
            ref_pos: 5,
            seed: 0,
        };
        assert_eq!(m.diagonal(), 5);
    }

    #[test]
    fn test_self_seed_match() {
        let r = default_reduction();
        let shape = Shape::from_code("111", &r);

        let seq: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let seqs: Vec<&[Letter]> = vec![&seq];

        let index = build_seed_index(&seqs, &shape, &r);
        let matches = find_seed_matches(&seqs, &index, &shape, &r);

        // Self-match should find at least one match per seed position
        assert!(matches.len() >= 8); // 10 - 3 + 1 = 8 positions
    }
}
