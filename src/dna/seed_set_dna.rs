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
        Some(self.cmp(other))
    }
}

impl Ord for SeedMatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.id() < other.id()
            || (self.id() == other.id() && self.ungapped_score() > other.ungapped_score())
        {
            std::cmp::Ordering::Greater
        } else if self == other {
            std::cmp::Ordering::Equal
        } else {
            std::cmp::Ordering::Less
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
}
