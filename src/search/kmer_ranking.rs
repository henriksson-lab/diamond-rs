use crate::basic::value::BlockId;
use crate::data::flags as data_flags;
use crate::data::sequence_set::SequenceSet;
use crate::util::data_structures::DoubleArray;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PackedLocId {
    pub loc: u64,
    pub block_id: BlockId,
}

impl PackedLocId {
    pub fn new(loc: u64, block_id: BlockId) -> Self {
        Self { loc, block_id }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct KmerRanking {
    pub rank: Vec<f32>,
}

impl KmerRanking {
    pub const PACKED_LOC_ID_SIZE: usize = 9;

    /// Matches unsupported C++ `KmerRanking(... DoubleArray<PackedLoc>*, ...)`.
    pub fn from_packed_loc_seed_hits() -> Result<Self, String> {
        Err("Unsupported".to_string())
    }

    /// Matches C++ `KmerRanking(const SequenceSet&, SeedPartition, DoubleArray<PackedLocId>*, DoubleArray<PackedLocId>*)`.
    pub fn from_packed_loc_id_seed_hits(
        queries: &SequenceSet,
        seedp_count: usize,
        query_seed_hits: &[DoubleArray],
        ref_seed_hits: &[DoubleArray],
    ) -> Self {
        let query_count = queries.len();
        let mut rank = vec![0.0f32; query_count];
        for p in 0..seedp_count {
            let mut q_it = query_seed_hits[p].begin_with_elem_size(Self::PACKED_LOC_ID_SIZE);
            let mut s_it = ref_seed_hits[p].begin_with_elem_size(Self::PACKED_LOC_ID_SIZE);
            while q_it.good() && s_it.good() {
                let contribution = (s_it.count() as f32).sqrt();
                let query_hits = q_it.range_bytes(Self::PACKED_LOC_ID_SIZE);
                for hit in query_hits.chunks_exact(Self::PACKED_LOC_ID_SIZE) {
                    let block_id = u32::from_ne_bytes(hit[5..9].try_into().unwrap()) as usize;
                    rank[block_id] += contribution;
                }
                q_it.increment(Self::PACKED_LOC_ID_SIZE);
                s_it.increment(Self::PACKED_LOC_ID_SIZE);
            }
        }
        Self { rank }
    }

    /// Matches C++ `KmerRanking(const SequenceSet& queries)`.
    pub fn from_queries(queries: &SequenceSet) -> Self {
        let mut rank = Vec::with_capacity(queries.len());
        for i in 0..queries.len() {
            rank.push(queries.seq_length(i) as f32);
        }
        Self { rank }
    }

    /// Matches the supported C++ `KmerRanking` constructor for `PackedLocId` seed hits.
    ///
    /// Each tuple is one joined seed bucket: query hits for the seed and the target hit count.
    /// The contribution to each query block is `sqrt(target_count)`.
    pub fn from_joined_seed_hits(
        query_count: usize,
        joined_hits: &[(&[PackedLocId], usize)],
    ) -> Self {
        let mut rank = vec![0.0f32; query_count];
        for &(query_hits, target_count) in joined_hits {
            let contribution = (target_count as f32).sqrt();
            for hit in query_hits {
                rank[hit.block_id as usize] += contribution;
            }
        }
        Self { rank }
    }

    /// Matches C++ `KmerRanking::highest_ranking`.
    pub fn highest_ranking(&self, hits: &[PackedLocId]) -> i32 {
        let Some(first) = hits.first() else {
            return 0;
        };
        let mut r = 0usize;
        let mut rank = self.rank[first.block_id as usize];
        for (i, hit) in hits.iter().enumerate().skip(1) {
            if self.rank[hit.block_id as usize] > rank {
                rank = self.rank[hit.block_id as usize];
                r = i;
            }
        }
        r as i32
    }

    /// Matches C++ `KmerRanking::highest_ranking` for the packed-layout seed location.
    pub fn highest_ranking_data_flags(&self, hits: &[data_flags::PackedLocId]) -> i32 {
        let Some(first) = hits.first() else {
            return 0;
        };
        let first_block_id = first.block_id;
        let mut r = 0usize;
        let mut rank = self.rank[first_block_id as usize];
        for (i, hit) in hits.iter().enumerate().skip(1) {
            let block_id = hit.block_id;
            if self.rank[block_id as usize] > rank {
                rank = self.rank[block_id as usize];
                r = i;
            }
        }
        r as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_loc::PackedLoc;

    #[test]
    fn test_kmer_ranking_from_queries_and_highest() {
        let mut queries = SequenceSet::new();
        queries.push(&[0, 1, 2]);
        queries.push(&[0, 1, 2, 3, 4]);
        queries.push(&[0, 1]);
        let ranking = KmerRanking::from_queries(&queries);
        assert_eq!(ranking.rank, vec![3.0, 5.0, 2.0]);
        let hits = [
            PackedLocId::new(10, 0),
            PackedLocId::new(20, 2),
            PackedLocId::new(30, 1),
        ];
        assert_eq!(ranking.highest_ranking(&hits), 2);
    }

    #[test]
    fn test_kmer_ranking_from_joined_seed_hits() {
        let q1 = [PackedLocId::new(10, 0), PackedLocId::new(11, 2)];
        let q2 = [PackedLocId::new(20, 1), PackedLocId::new(21, 2)];
        let ranking = KmerRanking::from_joined_seed_hits(3, &[(&q1, 4), (&q2, 9)]);
        assert_eq!(ranking.rank[0], 2.0);
        assert_eq!(ranking.rank[1], 3.0);
        assert_eq!(ranking.rank[2], 5.0);
    }

    #[test]
    fn test_kmer_ranking_from_packed_loc_unsupported() {
        assert_eq!(
            KmerRanking::from_packed_loc_seed_hits().unwrap_err(),
            "Unsupported"
        );
    }

    #[test]
    fn test_kmer_ranking_from_packed_loc_id_seed_hits() {
        let mut queries = SequenceSet::new();
        queries.push(&[0, 1, 2]);
        queries.push(&[0, 1, 2, 3]);
        queries.push(&[0, 1]);

        let mut query_group = Vec::new();
        query_group.extend_from_slice(&2u32.to_ne_bytes());
        data_flags::PackedLocId::with_block_id(PackedLoc::new(10), 0)
            .write_seed_hit_bytes(&mut query_group);
        data_flags::PackedLocId::with_block_id(PackedLoc::new(20), 2)
            .write_seed_hit_bytes(&mut query_group);
        query_group.extend_from_slice(&1u32.to_ne_bytes());
        data_flags::PackedLocId::with_block_id(PackedLoc::new(30), 1)
            .write_seed_hit_bytes(&mut query_group);

        let mut ref_group = Vec::new();
        ref_group.extend_from_slice(&4u32.to_ne_bytes());
        for i in 0..4 {
            data_flags::PackedLocId::with_block_id(PackedLoc::new(100 + i), 0)
                .write_seed_hit_bytes(&mut ref_group);
        }
        ref_group.extend_from_slice(&9u32.to_ne_bytes());
        for i in 0..9 {
            data_flags::PackedLocId::with_block_id(PackedLoc::new(200 + i), 0)
                .write_seed_hit_bytes(&mut ref_group);
        }

        let query_hits = [DoubleArray::from_bytes(
            query_group.clone(),
            query_group.len(),
        )];
        let ref_hits = [DoubleArray::from_bytes(ref_group.clone(), ref_group.len())];
        let ranking =
            KmerRanking::from_packed_loc_id_seed_hits(&queries, 1, &query_hits, &ref_hits);

        assert_eq!(ranking.rank, vec![2.0, 3.0, 2.0]);
    }

    #[test]
    fn test_highest_ranking_data_flags() {
        let ranking = KmerRanking {
            rank: vec![1.0, 5.0, 3.0],
        };
        let hits = [
            data_flags::PackedLocId::with_block_id(PackedLoc::new(10), 0),
            data_flags::PackedLocId::with_block_id(PackedLoc::new(20), 2),
            data_flags::PackedLocId::with_block_id(PackedLoc::new(30), 1),
        ];
        assert_eq!(ranking.highest_ranking_data_flags(&hits), 2);
    }

    trait TestPackedLocIdBytes {
        fn write_seed_hit_bytes(self, out: &mut Vec<u8>);
    }

    impl TestPackedLocIdBytes for data_flags::PackedLocId {
        fn write_seed_hit_bytes(self, out: &mut Vec<u8>) {
            let pos = self.pos;
            out.push(pos.high);
            out.extend_from_slice(&pos.low.to_ne_bytes());
            let block_id = self.block_id;
            out.extend_from_slice(&block_id.to_ne_bytes());
        }
    }
}
