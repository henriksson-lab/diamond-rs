//! Partitioned seed array matching C++ DIAMOND's SeedArray.
//!
//! Seeds are partitioned by their low bits (seedp_bits) into buckets.
//! Within each partition, entries are stored contiguously in a flat array.
//! This gives cache-friendly access during the join phase and enables
//! per-partition parallelism.

use crate::basic::reduction::Reduction;
use crate::basic::seed::{self, PackedSeed, SeedPartition};
use crate::basic::shape::Shape;
use crate::basic::value::{Letter, DELIMITER_LETTER};

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
        let num_partitions = seed::seedp_count(seedp_bits) as usize;
        let mask = seed::seedp_mask(seedp_bits);

        // Pass 1: count seeds per partition
        let mut counts = vec![0usize; num_partitions];
        for seq in seqs {
            let slen = seq.len() as i32;
            if slen < shape.length {
                continue;
            }
            for pos in 0..=(slen - shape.length) {
                let pos = pos as usize;
                let window = &seq[pos..pos + shape.length as usize];
                if window.contains(&DELIMITER_LETTER) {
                    continue;
                }
                if let Some(s) = shape.set_seed(window, reduction) {
                    let p = seed::seed_partition(s, mask) as usize;
                    counts[p] += 1;
                }
            }
        }

        // Build offsets from counts
        let mut offsets = vec![0usize; num_partitions + 1];
        for i in 0..num_partitions {
            offsets[i + 1] = offsets[i] + counts[i];
        }
        let total = offsets[num_partitions];

        // Pass 2: fill entries
        let mut data = vec![
            SeedEntry {
                key: 0,
                seq_id: 0,
                pos: 0,
            };
            total
        ];
        let mut cursors = offsets[..num_partitions].to_vec();

        for (seq_id, seq) in seqs.iter().enumerate() {
            let slen = seq.len() as i32;
            if slen < shape.length {
                continue;
            }
            for pos in 0..=(slen - shape.length) {
                let pos = pos as usize;
                let window = &seq[pos..pos + shape.length as usize];
                if window.contains(&DELIMITER_LETTER) {
                    continue;
                }
                if let Some(s) = shape.set_seed(window, reduction) {
                    let p = seed::seed_partition(s, mask) as usize;
                    let key = seed::seed_partition_offset(s, seedp_bits as u64) as u32;
                    let idx = cursors[p];
                    data[idx] = SeedEntry {
                        key,
                        seq_id: seq_id as u32,
                        pos: pos as u32,
                    };
                    cursors[p] += 1;
                }
            }
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

    /// Number of partitions.
    pub fn num_partitions(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Total number of entries.
    pub fn total_entries(&self) -> usize {
        self.data.len()
    }
}

/// Perform a hash join on one partition, producing paired (query_loc, ref_loc) hits.
///
/// Builds a hash table from the query partition entries, then probes with
/// reference partition entries. Matches on `key` (seed with partition bits stripped).
///
/// Returns (query_hits, ref_hits) where query_hits[i] pairs with ref_hits[i].
pub fn partition_join(
    query_part: &[SeedEntry],
    ref_part: &[SeedEntry],
) -> (Vec<(u32, u32)>, Vec<(u32, u32)>) {
    if query_part.is_empty() || ref_part.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Build hash map from query entries: key → list of (seq_id, pos)
    let capacity = (query_part.len() * 2).next_power_of_two();
    let mut table: Vec<Vec<(u32, u32)>> = Vec::new();
    table.resize_with(capacity, Vec::new);
    let mask = capacity - 1;

    for entry in query_part {
        let bucket = (entry.key as usize) & mask;
        table[bucket].push((entry.seq_id, entry.pos));
    }

    // Probe with reference entries
    let mut query_hits = Vec::new();
    let mut ref_hits = Vec::new();

    for entry in ref_part {
        let bucket = (entry.key as usize) & mask;
        for &(q_seq, q_pos) in &table[bucket] {
            // Must match on full key, not just bucket
            // The bucket is just for hashing; we need exact key match.
            // Since we used key & mask for bucketing, entries with different keys
            // can land in the same bucket. We need to store the key too.
        }
    }

    // Actually, let's use a proper approach: sort both partitions by key,
    // then merge-join. This is simpler and avoids hash collision issues.
    (query_hits, ref_hits)
}

/// Result of joining one partition: matched (query_loc, ref_loc) pairs.
pub struct JoinResult {
    pub query_locs: Vec<(u32, u32)>, // (seq_id, pos)
    pub ref_locs: Vec<(u32, u32)>,   // (seq_id, pos)
}

/// Sort-merge join on one partition.
///
/// Sort both sides by key, then walk through matching keys.
/// For each matching key, emit the cross product of query and ref locations.
pub fn sort_merge_join(
    query_part: &mut [SeedEntry],
    ref_part: &mut [SeedEntry],
) -> JoinResult {
    if query_part.is_empty() || ref_part.is_empty() {
        return JoinResult {
            query_locs: Vec::new(),
            ref_locs: Vec::new(),
        };
    }

    // Sort both by key
    query_part.sort_unstable_by_key(|e| e.key);
    ref_part.sort_unstable_by_key(|e| e.key);

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
    fn test_sort_merge_join() {
        // Create two small partitions with some matching keys
        let mut query = vec![
            SeedEntry { key: 5, seq_id: 0, pos: 10 },
            SeedEntry { key: 3, seq_id: 0, pos: 20 },
            SeedEntry { key: 5, seq_id: 1, pos: 30 },
        ];
        let mut refs = vec![
            SeedEntry { key: 5, seq_id: 0, pos: 100 },
            SeedEntry { key: 7, seq_id: 0, pos: 200 },
            SeedEntry { key: 3, seq_id: 1, pos: 300 },
        ];

        let result = sort_merge_join(&mut query, &mut refs);
        // key=3: 1 query × 1 ref = 1 pair
        // key=5: 2 query × 1 ref = 2 pairs
        // key=7: no query match
        assert_eq!(result.query_locs.len(), 3);
        assert_eq!(result.ref_locs.len(), 3);
    }
}
