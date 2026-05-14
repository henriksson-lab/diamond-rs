use crate::basic::packed_loc::PackedLoc;
use crate::basic::seed::{
    seed_partition, seed_partition_offset, seedp_mask, PackedSeed, SeedOffset,
};
use crate::data::seed_array::SeedArrayEntry;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

#[derive(Debug, Clone)]
pub struct Index {
    seed_arr: Vec<Vec<SeedArrayEntry<PackedLoc>>>,
    dna_index: Vec<HashMap<SeedOffset, usize>>,
    minimizer_counts: Vec<unsigned>,
    n_minimizer: unsigned,
    seedp_bits: i32,
}

#[allow(non_camel_case_types)]
type unsigned = u32;

impl Index {
    pub fn new(
        seed_arr: Vec<Vec<SeedArrayEntry<PackedLoc>>>,
        seedp_bits: i32,
        repetitive_cutoff: f64,
    ) -> Self {
        let mut index = Self {
            dna_index: vec![HashMap::new(); seed_arr.len()],
            minimizer_counts: vec![0; seed_arr.len()],
            n_minimizer: 0,
            seed_arr,
            seedp_bits,
        };
        index.count_minimizers(0..index.seed_arr.len());
        let cutoff = index.filter_repetitive(0..index.seed_arr.len(), repetitive_cutoff);
        index.build_index(0..index.seed_arr.len(), cutoff);
        index
    }

    pub fn contains(&self, seed: PackedSeed) -> Option<&[SeedArrayEntry<PackedLoc>]> {
        let partition = seed_partition(seed, seedp_mask(self.seedp_bits)) as usize;
        let key = seed_partition_offset(seed, self.seedp_bits as u64);

        if partition >= self.dna_index.len() || self.dna_index[partition].is_empty() {
            return None;
        }

        let first = *self.dna_index[partition].get(&key)?;
        let entries = &self.seed_arr[partition];
        let mut end = entries.len();
        for i in first + 1..entries.len() {
            if entries[i].key != entries[first].key {
                end = i;
                break;
            }
        }
        Some(&entries[first..end])
    }

    pub fn count_worker(&mut self, part: usize) {
        self.seed_arr[part].sort();
        let mut count = 0u32;
        let entries = &self.seed_arr[part];
        let mut i = 0usize;
        while i < entries.len() {
            count += 1;
            let key = entries[i].key;
            i += 1;
            while i < entries.len() && entries[i].key == key {
                i += 1;
            }
        }
        self.minimizer_counts[part] = count;
    }

    pub fn count_minimizers(&mut self, range: std::ops::Range<usize>) {
        for part in range {
            self.count_worker(part);
        }
        self.n_minimizer = self.minimizer_counts.iter().copied().sum();
    }

    pub fn index_worker(&mut self, part: usize, cutoff: i32) {
        self.dna_index[part].clear();
        let entries = &self.seed_arr[part];
        let mut i = 0usize;
        while i < entries.len() {
            let key = entries[i].key;
            let first = i;
            i += 1;
            while i < entries.len() && entries[i].key == key {
                i += 1;
            }
            if (i - first) < cutoff as usize {
                self.dna_index[part].insert(key, first);
            }
        }
    }

    pub fn build_index(&mut self, range: std::ops::Range<usize>, repetitive_cutoff: i32) {
        for part in range {
            self.index_worker(part, repetitive_cutoff);
        }
    }

    pub fn filter_worker(&self, part: usize, rep_thread: &mut BinaryHeap<Reverse<i32>>, n: usize) {
        let entries = &self.seed_arr[part];
        let mut i = 0usize;
        while i < entries.len() {
            let key = entries[i].key;
            let first = i;
            i += 1;
            while i < entries.len() && entries[i].key == key {
                i += 1;
            }
            let count = (i - first) as i32;
            if rep_thread.len() < n {
                rep_thread.push(Reverse(count));
            } else if rep_thread.peek().is_some_and(|top| top.0 < count) {
                rep_thread.push(Reverse(count));
                rep_thread.pop();
            }
        }
    }

    pub fn filter_repetitive(&self, range: std::ops::Range<usize>, repetitive_cutoff: f64) -> i32 {
        let n = (self.n_minimizer as f64 * repetitive_cutoff) as usize;

        if n < 1 {
            return i32::MAX;
        }

        let mut repetitive_max = BinaryHeap::new();
        for part in range {
            self.filter_worker(part, &mut repetitive_max, n);
        }

        repetitive_max.peek().map(|x| x.0).unwrap_or(i32::MAX)
    }

    pub fn minimizer_counts(&self) -> &[unsigned] {
        &self.minimizer_counts
    }

    pub fn n_minimizer(&self) -> unsigned {
        self.n_minimizer
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn loc(x: u64) -> PackedLoc {
        PackedLoc::new(x)
    }

    fn entry(key: SeedOffset, value: u64) -> SeedArrayEntry<PackedLoc> {
        SeedArrayEntry::new(key, loc(value))
    }

    fn seed(partition: u64, key: u64, seedp_bits: i32) -> PackedSeed {
        (key << seedp_bits) | partition
    }

    #[test]
    fn test_count_minimizers_sorts_and_counts_groups() {
        let mut index = Index {
            seed_arr: vec![
                vec![entry(3, 0), entry(1, 1), entry(1, 2)],
                vec![entry(2, 3), entry(2, 4), entry(5, 5)],
            ],
            dna_index: vec![HashMap::new(), HashMap::new()],
            minimizer_counts: vec![0, 0],
            n_minimizer: 0,
            seedp_bits: 1,
        };

        index.count_minimizers(0..2);

        assert_eq!(index.minimizer_counts(), &[2, 2]);
        assert_eq!(index.n_minimizer(), 4);
        let first_key = index.seed_arr[0][0].key;
        let third_key = index.seed_arr[0][2].key;
        assert_eq!(first_key, 1);
        assert_eq!(third_key, 3);
    }

    #[test]
    fn test_filter_repetitive_returns_nth_largest_group_count() {
        let mut index = Index {
            seed_arr: vec![
                vec![entry(1, 0), entry(1, 1), entry(2, 2)],
                vec![
                    entry(3, 3),
                    entry(3, 4),
                    entry(3, 5),
                    entry(4, 6),
                    entry(4, 7),
                ],
            ],
            dna_index: vec![HashMap::new(), HashMap::new()],
            minimizer_counts: vec![0, 0],
            n_minimizer: 0,
            seedp_bits: 1,
        };
        index.count_minimizers(0..2);

        assert_eq!(index.filter_repetitive(0..2, 0.5), 2);
        assert_eq!(index.filter_repetitive(0..2, 0.0), i32::MAX);
    }

    #[test]
    fn test_build_index_and_contains_skip_repetitive_groups() {
        let seedp_bits = 2;
        let index = Index::new(
            vec![
                vec![entry(9, 0), entry(7, 1), entry(7, 2)],
                vec![entry(4, 3), entry(4, 4), entry(4, 5), entry(8, 6)],
                vec![],
                vec![entry(1, 7)],
            ],
            seedp_bits,
            0.5,
        );

        assert_eq!(index.minimizer_counts(), &[2, 2, 0, 1]);
        let unique = index.contains(seed(0, 9, seedp_bits)).unwrap();
        assert_eq!(unique.len(), 1);
        let unique_key = unique[0].key;
        assert_eq!(unique_key, 9);

        assert!(index.contains(seed(0, 7, seedp_bits)).is_none());
        assert!(index.contains(seed(1, 4, seedp_bits)).is_none());

        let part3 = index.contains(seed(3, 1, seedp_bits)).unwrap();
        assert_eq!(part3.len(), 1);
        let part3_key = part3[0].key;
        assert_eq!(part3_key, 1);
    }

    #[test]
    fn test_contains_missing_partition_or_key() {
        let index = Index::new(vec![vec![entry(2, 0)]], 1, 0.0);
        assert!(index.contains(seed(0, 3, 1)).is_none());
        assert!(index.contains(seed(1, 2, 1)).is_none());
    }
}
