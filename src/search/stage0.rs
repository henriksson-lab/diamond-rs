use crate::data::seed_array::{SeedArray, SeedLocation};
use crate::util::algo::{hash_join, Relation};
use crate::util::data_structures::DoubleArray;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchShapeSeedLocKind {
    PackedLoc,
    PackedLocId,
}

/// Matches C++ `seed_join_worker(...)`.
pub fn seed_join_worker<SeedLoc>(
    query_seeds: &SeedArray<SeedLoc>,
    ref_seeds: &SeedArray<SeedLoc>,
    seedp_begin: usize,
    partition_count: usize,
    query_seed_hits: &mut [DoubleArray],
    ref_seed_hits: &mut [DoubleArray],
) -> Result<(), String>
where
    SeedLoc: SeedLocation + crate::util::algo::HashJoinValue,
{
    let bits = query_seeds.key_bits();
    if bits != ref_seeds.key_bits() {
        return Err("Joining seed arrays with different key lengths.".to_string());
    }
    for p in seedp_begin..partition_count {
        let query = query_seeds.begin(p);
        let reference = ref_seeds.begin(p);
        let join = hash_join(
            Relation::new(query, query.len()),
            Relation::new(reference, reference.len()),
            bits as u32,
        );
        query_seed_hits[p] = join.0;
        ref_seed_hits[p] = join.1;
    }
    Ok(())
}

/// Matches C++ `search_worker(...)`.
pub fn search_worker<SeedLoc, F>(
    stop: bool,
    seedp_begin: usize,
    partition_count: usize,
    shape: u32,
    thread_id: usize,
    query_seed_hits: &[DoubleArray],
    ref_seed_hits: &[DoubleArray],
    mut run_stage1: F,
) -> usize
where
    SeedLoc: SeedLocation,
    F: FnMut(usize, u32, usize, &DoubleArray, &DoubleArray),
{
    let mut processed = 0usize;
    let mut p = seedp_begin;
    while !stop && p < partition_count {
        run_stage1(p, shape, thread_id, &query_seed_hits[p], &ref_seed_hits[p]);
        processed += 1;
        p += 1;
    }
    processed
}

/// Matches C++ `search_shape(unsigned...)`.
pub fn search_shape_seed_loc_kind(keep_target_id: bool) -> SearchShapeSeedLocKind {
    if keep_target_id {
        SearchShapeSeedLocKind::PackedLocId
    } else {
        SearchShapeSeedLocKind::PackedLoc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_loc::PackedLoc;
    use crate::basic::reduction::Reduction;
    use crate::basic::shape_config::ShapeConfig;
    use crate::basic::value::SequenceType;
    use crate::data::block::Block;
    use crate::data::enum_seeds::EnumSeedsContext;
    use crate::data::flags::{EnumCfg, PackedLocId, SeedEncoding, NO_FILTER};
    use crate::data::seed_histogram::SeedPartitionRange;
    use crate::masking::MaskingAlgo;

    fn enum_cfg<'a>(partition: Option<&'a Vec<u32>>) -> EnumCfg<'a> {
        EnumCfg {
            partition,
            shape_begin: 0,
            shape_end: 1,
            code: SeedEncoding::SpacedFactor,
            skip: None,
            filter_masked_seeds: false,
            mask_seeds: false,
            seed_cut: 0.0,
            soft_masking: MaskingAlgo::None,
            minimizer_window: 0,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: 0,
        }
    }

    fn block(seq: &[i8]) -> Block {
        let mut block = Block::new();
        block
            .push_back(seq, Some("s0"), None, 0, SequenceType::AminoAcid, 0, false)
            .unwrap();
        block
    }

    #[test]
    fn test_search_shape_seed_loc_kind() {
        assert_eq!(
            search_shape_seed_loc_kind(false),
            SearchShapeSeedLocKind::PackedLoc
        );
        assert_eq!(
            search_shape_seed_loc_kind(true),
            SearchShapeSeedLocKind::PackedLocId
        );
    }

    #[test]
    fn test_seed_join_worker_joins_matching_partition() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let seedp_bits = 1;
        let range = SeedPartitionRange::with_bounds(0, 2);
        let partition = vec![0, 1];
        let cfg = enum_cfg(Some(&partition));
        let mut query = block(&[0, 1, 2, 3, 4]);
        let mut reference = block(&[9, 0, 1, 2, 8]);

        let query_seeds = SeedArray::<PackedLoc>::one_pass(
            &mut query, &range, seedp_bits, &NO_FILTER, &cfg, &ctx, 1,
        )
        .unwrap();
        let ref_seeds = SeedArray::<PackedLoc>::one_pass(
            &mut reference,
            &range,
            seedp_bits,
            &NO_FILTER,
            &cfg,
            &ctx,
            1,
        )
        .unwrap();
        let mut query_hits = vec![DoubleArray::new(), DoubleArray::new()];
        let mut ref_hits = vec![DoubleArray::new(), DoubleArray::new()];

        seed_join_worker(
            &query_seeds,
            &ref_seeds,
            0,
            range.size() as usize,
            &mut query_hits,
            &mut ref_hits,
        )
        .unwrap();

        let total_query_bytes: u32 = query_hits.iter().map(|x| x.data().len() as u32).sum();
        let total_ref_bytes: u32 = ref_hits.iter().map(|x| x.data().len() as u32).sum();
        assert!(total_query_bytes > 0);
        assert!(total_ref_bytes > 0);
    }

    #[test]
    fn test_seed_join_worker_rejects_key_bit_mismatch() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let range = SeedPartitionRange::with_bounds(0, 2);
        let partition = vec![0, 1];
        let cfg = enum_cfg(Some(&partition));
        let mut query = block(&[0, 1, 2, 3]);
        let mut reference = block(&[0, 1, 2, 3]);
        let query_seeds =
            SeedArray::<PackedLocId>::one_pass(&mut query, &range, 1, &NO_FILTER, &cfg, &ctx, 1)
                .unwrap();
        let ref_seeds = SeedArray::<PackedLocId>::one_pass(
            &mut reference,
            &range,
            2,
            &NO_FILTER,
            &cfg,
            &ctx,
            1,
        )
        .unwrap();
        let mut query_hits = vec![DoubleArray::new(), DoubleArray::new()];
        let mut ref_hits = vec![DoubleArray::new(), DoubleArray::new()];

        assert_eq!(
            seed_join_worker(
                &query_seeds,
                &ref_seeds,
                0,
                range.size() as usize,
                &mut query_hits,
                &mut ref_hits,
            )
            .unwrap_err(),
            "Joining seed arrays with different key lengths."
        );
    }

    #[test]
    fn test_search_worker_stops_or_processes_partitions() {
        let query_hits = vec![DoubleArray::new(), DoubleArray::new(), DoubleArray::new()];
        let ref_hits = vec![DoubleArray::new(), DoubleArray::new(), DoubleArray::new()];
        let mut seen = Vec::new();
        let n = search_worker::<PackedLoc, _>(
            false,
            1,
            3,
            7,
            2,
            &query_hits,
            &ref_hits,
            |p, shape, thread_id, _q, _r| seen.push((p, shape, thread_id)),
        );
        assert_eq!(n, 2);
        assert_eq!(seen, vec![(1, 7, 2), (2, 7, 2)]);

        let n = search_worker::<PackedLoc, _>(
            true,
            0,
            3,
            7,
            2,
            &query_hits,
            &ref_hits,
            |_p, _shape, _thread_id, _q, _r| unreachable!(),
        );
        assert_eq!(n, 0);
    }
}
