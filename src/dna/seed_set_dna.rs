use crate::basic::reduction::Reduction;
use crate::basic::seed_iterator::MinimizerIterator;
use crate::basic::shape::Shape;
use crate::basic::value::{BlockId, Letter, Loc};
use crate::data::sequence_set::LetterStringSet;
use crate::dna::dna_index::Index;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeedMatch {
    i: Loc,
    j: Loc,
    target_id: BlockId,
    score: i32,
}

impl SeedMatch {
    pub fn new(i: Loc, id: BlockId, j: Loc, shape_length: i32) -> Self {
        Self {
            i,
            target_id: id,
            j,
            score: shape_length,
        }
    }

    pub fn ungapped_score(&self) -> i32 {
        self.score
    }

    pub fn score(&mut self, scr: i32) {
        self.score = scr;
    }

    pub fn i(&self) -> Loc {
        self.i
    }

    pub fn j(&self) -> Loc {
        self.j
    }

    pub fn id(&self) -> BlockId {
        self.target_id
    }

    pub fn i_start(&self) -> i32 {
        self.i - self.score
    }

    pub fn j_start(&self) -> i32 {
        self.j - self.score
    }
}

impl PartialOrd for SeedMatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        // C++ defines only `operator>`: hits with the same target and score
        // are incomparable unless every stored field is equal.  In
        // particular, their query/reference coordinates are not tie-breakers.
        if self.id() < other.id()
            || (self.id() == other.id() && self.ungapped_score() > other.ungapped_score())
        {
            Some(std::cmp::Ordering::Greater)
        } else if other.id() < self.id()
            || (self.id() == other.id() && other.ungapped_score() > self.ungapped_score())
        {
            Some(std::cmp::Ordering::Less)
        } else if self == other {
            Some(std::cmp::Ordering::Equal)
        } else {
            None
        }
    }
}

pub fn seed_lookup(
    query: &[Letter],
    target_seqs: &LetterStringSet,
    filter: &Index,
    window_size: Loc,
    seed_shape: &Shape,
    reduction: &Reduction,
) -> Vec<SeedMatch> {
    let query_sequence = query.to_vec();
    let mut seed_matches = Vec::new();
    // The C++ iterator's end pointer precedes its begin pointer for a query
    // shorter than the shape, so it immediately reports `good() == false`.
    // Guard explicitly because Rust's saturating index arithmetic otherwise
    // constructs one invalid candidate position.
    if query_sequence.len() < seed_shape.length.max(0) as usize {
        return seed_matches;
    }
    let mut it = MinimizerIterator::new(&query_sequence, seed_shape, window_size, reduction);

    while it.good() {
        let key = it.get();
        if let Some(ref_position) = filter.contains(key) {
            for pos in ref_position {
                let value = pos.value;
                let (id, loc) = target_seqs.local_position(value.as_i64());
                seed_matches.push(SeedMatch::new(it.pos(), id, loc, seed_shape.length));
            }
        }
        it.increment();
    }

    seed_matches
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_loc::PackedLoc;
    use crate::basic::seed::{seed_partition, seed_partition_offset, seedp_mask};
    use crate::data::seed_array::SeedArrayEntry;

    fn build_index(
        target_seqs: &LetterStringSet,
        seed_shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
    ) -> Index {
        let mut seed_arr = vec![Vec::new(); 1usize << seedp_bits];
        let mask = seedp_mask(seedp_bits);

        for id in 0..target_seqs.size() {
            let seq = target_seqs.get(id as usize);
            if seq.len() < seed_shape.length as usize {
                continue;
            }
            for j in 0..=seq.len() - seed_shape.length as usize {
                if let Some(seed) = seed_shape.set_seed_reduced(&seq[j..], reduction) {
                    let part = seed_partition(seed, mask) as usize;
                    let key = seed_partition_offset(seed, seedp_bits as u64);
                    let pos = target_seqs.position(id, j as Loc);
                    seed_arr[part].push(SeedArrayEntry::new(key, PackedLoc::new(pos as u64)));
                }
            }
        }

        Index::new(seed_arr, seedp_bits, 0.0)
    }

    #[test]
    fn test_seed_lookup_converts_index_hits_to_seed_matches() {
        let reduction = Reduction::default_reduction();
        let seed_shape = Shape::from_code("11", &reduction);
        let mut target_seqs = LetterStringSet::new();
        target_seqs.push_back(&[1, 2, 9, 9]);
        target_seqs.push_back(&[3, 4, 1, 2]);
        let index = build_index(&target_seqs, &seed_shape, &reduction, 2);

        let query = [1, 2, 3, 4];
        let mut matches = seed_lookup(&query, &target_seqs, &index, 1, &seed_shape, &reduction)
            .into_iter()
            .map(|m| (m.i(), m.id(), m.j(), m.ungapped_score()))
            .collect::<Vec<_>>();
        matches.sort();

        assert_eq!(matches, vec![(0, 0, 0, 2), (0, 1, 2, 2), (2, 1, 0, 2)]);
    }

    #[test]
    fn test_seed_lookup_skips_missing_seeds_and_short_queries() {
        let reduction = Reduction::default_reduction();
        let seed_shape = Shape::from_code("11", &reduction);
        let mut target_seqs = LetterStringSet::new();
        target_seqs.push_back(&[1, 2]);
        let index = build_index(&target_seqs, &seed_shape, &reduction, 2);

        assert!(seed_lookup(&[3, 4], &target_seqs, &index, 1, &seed_shape, &reduction).is_empty());
        assert!(seed_lookup(&[1], &target_seqs, &index, 1, &seed_shape, &reduction).is_empty());
        assert!(seed_lookup(&[], &target_seqs, &index, 1, &seed_shape, &reduction).is_empty());
    }

    #[test]
    fn test_seed_match_accessors_score_update_and_starts() {
        let mut hit = SeedMatch::new(17, 3, 29, 4);
        assert_eq!(
            (hit.i(), hit.id(), hit.j(), hit.ungapped_score()),
            (17, 3, 29, 4)
        );
        assert_eq!((hit.i_start(), hit.j_start()), (13, 25));

        hit.score(7);
        assert_eq!(hit.ungapped_score(), 7);
        assert_eq!((hit.i_start(), hit.j_start()), (10, 22));
    }

    #[test]
    fn test_seed_match_greater_than_matches_cpp_partial_relation() {
        let strong = SeedMatch::new(5, 2, 8, 9);
        let weak = SeedMatch::new(5, 2, 8, 4);
        let earlier_target = SeedMatch::new(5, 1, 8, 1);
        assert!(strong > weak);
        assert!(earlier_target > strong);

        let same_rank_different_coordinates = SeedMatch::new(6, 2, 9, 9);
        assert!(!(strong > same_rank_different_coordinates));
        assert!(!(same_rank_different_coordinates > strong));
        assert_eq!(strong.partial_cmp(&same_rank_different_coordinates), None);
    }
}
