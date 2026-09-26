use crate::basic::packed_loc::PackedLoc;
use crate::basic::reduction::Reduction;
use crate::basic::seed::{seedp_count, seedp_mask};
use crate::basic::shape_config::ShapeConfig;
use crate::data::block::Block;
use crate::data::enum_seeds::EnumSeedsContext;
use crate::data::flags::{EnumCfg, PackedLocId, SeedEncoding, NO_FILTER};
use crate::data::frequent_seeds::{FrequentSeeds, FrequentSeedsBuildStats, FrequentSeedsConfig};
use crate::data::seed_array::{SeedArray, SeedLocation};
use crate::data::seed_histogram::{SeedPartitionRange, CURRENT_RANGE};
use crate::data::seed_set::{HashedSeedSet, SeedSet};
use crate::masking::MaskingAlgo;
use crate::search::kmer_ranking::KmerRanking;
use crate::search::seed_complexity::{mask_seeds, MaskSeedsStats, SeedLocBytes};
use crate::util::algo::Partition;
use crate::util::algo::{hash_join, Relation};
use crate::util::data_structures::DoubleArray;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchShapeSeedLocKind {
    PackedLoc,
    PackedLocId,
}

/// The two process-global query-seed filters used by C++ stage 0.
#[derive(Clone, Copy, Default)]
pub enum ReferenceSeedFilter<'a> {
    #[default]
    None,
    BitSet(&'a SeedSet),
    Hashed(&'a HashedSeedSet),
}

/// Explicit inputs for the C++ `search_shape` orchestration.
pub struct SearchShapeConfig<'a> {
    pub shapes: &'a ShapeConfig,
    pub reduction: &'a Reduction,
    pub seedp_bits: i32,
    pub index_chunks: usize,
    pub threads: usize,
    pub seed_encoding: SeedEncoding,
    pub query_skip: Option<&'a Vec<bool>>,
    pub seed_complexity_cut: f64,
    pub soft_masking: MaskingAlgo,
    pub minimizer_window: i32,
    pub sketch_size: i32,
    pub min_query_len: i32,
    pub query_contexts: usize,
    pub reference_filter: ReferenceSeedFilter<'a>,
    pub target_seeds: Option<&'a HashedSeedSet>,
    pub freq_masking: bool,
    pub freq_sd: f64,
    pub linear_stage1_query: bool,
    pub linear_stage1_target: bool,
    pub kmer_ranking: bool,
    pub keep_target_id: bool,
    pub ungapped_raw_score: i32,
}

impl<'a> SearchShapeConfig<'a> {
    pub fn new(shapes: &'a ShapeConfig, reduction: &'a Reduction) -> Self {
        Self {
            shapes,
            reduction,
            seedp_bits: 10,
            index_chunks: 1,
            threads: 1,
            seed_encoding: SeedEncoding::SpacedFactor,
            query_skip: None,
            seed_complexity_cut: 0.0,
            soft_masking: MaskingAlgo::None,
            minimizer_window: 0,
            sketch_size: 0,
            min_query_len: 0,
            query_contexts: 1,
            reference_filter: ReferenceSeedFilter::None,
            target_seeds: None,
            freq_masking: false,
            freq_sd: 0.0,
            linear_stage1_query: false,
            linear_stage1_target: false,
            kmer_ranking: false,
            keep_target_id: false,
            ungapped_raw_score: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SearchShapeChunkReport {
    pub chunk: usize,
    pub range: SeedPartitionRange,
    pub query_seed_count: usize,
    pub reference_seed_count: usize,
    pub mask_stats: MaskSeedsStats,
    pub frequent_seed_stats: Option<FrequentSeedsBuildStats>,
    pub searched_partitions: usize,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SearchShapeReport {
    pub seed_loc_kind: Option<SearchShapeSeedLocKind>,
    pub chunks: Vec<SearchShapeChunkReport>,
}

/// Per-partition state passed to the stage-1 boundary.
pub struct SearchShapePartition<'a> {
    pub shape_id: u32,
    pub relative_partition: usize,
    pub absolute_partition: u32,
    pub query_seed_hits: &'a DoubleArray,
    pub reference_seed_hits: &'a DoubleArray,
    pub previous_shape_patterns: &'a [u32],
    pub shape_patterns: &'a [u32],
    pub ungapped_raw_score: i32,
    pub seedp_mask: u64,
    pub kmer_ranking: Option<&'a KmerRanking>,
}

pub trait Stage0SeedLocation:
    SeedLocation + SeedLocBytes + crate::util::algo::HashJoinValue + Send + Sync
{
    fn build_ranking(
        enabled: bool,
        query: &Block,
        query_hits: &[DoubleArray],
        reference_hits: &[DoubleArray],
    ) -> Option<KmerRanking>;
}

impl Stage0SeedLocation for PackedLoc {
    fn build_ranking(
        _enabled: bool,
        _query: &Block,
        _query_hits: &[DoubleArray],
        _reference_hits: &[DoubleArray],
    ) -> Option<KmerRanking> {
        None
    }
}

impl Stage0SeedLocation for PackedLocId {
    fn build_ranking(
        enabled: bool,
        query: &Block,
        query_hits: &[DoubleArray],
        reference_hits: &[DoubleArray],
    ) -> Option<KmerRanking> {
        Some(if enabled {
            KmerRanking::from_packed_loc_id_seed_hits(
                query.seqs(),
                query_hits.len(),
                query_hits,
                reference_hits,
            )
        } else {
            KmerRanking::from_queries(query.seqs())
        })
    }
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

fn search_shape_with_seed_loc<SeedLoc, F>(
    shape_id: u32,
    query: &mut Block,
    reference: &mut Block,
    config: &SearchShapeConfig<'_>,
    run_stage1: &mut F,
) -> Result<SearchShapeReport, String>
where
    SeedLoc: Stage0SeedLocation,
    F: FnMut(SearchShapePartition<'_>) -> Result<(), String>,
{
    if shape_id >= config.shapes.count() as u32 {
        return Err("Shape id out of range.".to_string());
    }
    if config.threads == 0 || config.threads > u32::MAX as usize {
        return Err("Stage 0 thread count must be positive.".to_string());
    }

    let partition_count = seedp_count(config.seedp_bits) as usize;
    let chunks = Partition::new(partition_count, config.index_chunks);
    let reference_histogram = reference.hst().get(shape_id as usize).clone();
    let reference_partition = reference.hst().partition().clone();
    let query_histogram = if config.target_seeds.is_none() {
        Some(query.hst().get(shape_id as usize).clone())
    } else {
        None
    };
    let query_partition = if config.target_seeds.is_none() {
        Some(query.hst().partition().clone())
    } else {
        None
    };
    let enum_context = EnumSeedsContext {
        shapes: config.shapes,
        reduction: config.reduction,
        min_query_len: config.min_query_len,
        query_contexts: config.query_contexts,
    };
    let patterns = config.shapes.patterns(0, shape_id + 1);
    let mut report = SearchShapeReport {
        seed_loc_kind: Some(search_shape_seed_loc_kind(config.keep_target_id)),
        chunks: Vec::with_capacity(chunks.parts),
    };

    for chunk in 0..chunks.parts {
        let range =
            SeedPartitionRange::with_bounds(chunks.begin(chunk) as u32, chunks.end(chunk) as u32);
        *CURRENT_RANGE.lock().map_err(|e| e.to_string())? = range;

        let reference_cfg = EnumCfg {
            partition: Some(&reference_partition),
            shape_begin: shape_id as i32,
            shape_end: shape_id as i32 + 1,
            code: config.seed_encoding,
            skip: None,
            filter_masked_seeds: false,
            mask_seeds: false,
            seed_cut: config.seed_complexity_cut,
            soft_masking: if matches!(config.reference_filter, ReferenceSeedFilter::None) {
                config.soft_masking
            } else {
                MaskingAlgo::None
            },
            minimizer_window: config.minimizer_window,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: config.sketch_size,
        };
        let reference_seeds = match config.reference_filter {
            ReferenceSeedFilter::None => SeedArray::<SeedLoc>::from_histogram(
                reference,
                &reference_histogram,
                &range,
                config.seedp_bits,
                &NO_FILTER,
                &reference_cfg,
                &enum_context,
            )?,
            ReferenceSeedFilter::BitSet(filter) => SeedArray::<SeedLoc>::from_histogram(
                reference,
                &reference_histogram,
                &range,
                config.seedp_bits,
                filter,
                &reference_cfg,
                &enum_context,
            )?,
            ReferenceSeedFilter::Hashed(filter) => SeedArray::<SeedLoc>::from_histogram(
                reference,
                &reference_histogram,
                &range,
                config.seedp_bits,
                filter,
                &reference_cfg,
                &enum_context,
            )?,
        };

        let query_cfg = EnumCfg {
            partition: query_partition.as_ref(),
            shape_begin: shape_id as i32,
            shape_end: shape_id as i32 + 1,
            code: config.seed_encoding,
            skip: config.query_skip,
            filter_masked_seeds: false,
            mask_seeds: true,
            seed_cut: config.seed_complexity_cut,
            soft_masking: config.soft_masking,
            minimizer_window: config.minimizer_window,
            filter_low_complexity_seeds: matches!(
                config.reference_filter,
                ReferenceSeedFilter::Hashed(_)
            ),
            mask_low_complexity_seeds: matches!(
                config.reference_filter,
                ReferenceSeedFilter::Hashed(_)
            ),
            sketch_size: config.sketch_size,
        };
        let query_seeds = if let Some(filter) = config.target_seeds {
            SeedArray::<SeedLoc>::one_pass(
                query,
                &range,
                config.seedp_bits,
                filter,
                &query_cfg,
                &enum_context,
                config.threads as u32,
            )?
        } else {
            SeedArray::<SeedLoc>::from_histogram(
                query,
                query_histogram.as_ref().unwrap(),
                &range,
                config.seedp_bits,
                &NO_FILTER,
                &query_cfg,
                &enum_context,
            )?
        };

        let relative_count = range.size() as usize;
        let mut query_hits = vec![DoubleArray::new(); relative_count];
        let mut reference_hits = vec![DoubleArray::new(); relative_count];
        seed_join_worker(
            &query_seeds,
            &reference_seeds,
            0,
            relative_count,
            &mut query_hits,
            &mut reference_hits,
        )?;

        let mut mask_stats = MaskSeedsStats::default();
        let mut frequent_seed_stats = None;
        if config.freq_masking && !config.linear_stage1_query && !config.linear_stage1_target {
            frequent_seed_stats = Some(FrequentSeeds::build_with_config::<SeedLoc>(
                shape_id,
                &range,
                &mut query_hits,
                &mut reference_hits,
                query.seqs_mut(),
                FrequentSeedsConfig {
                    freq_sd: config.freq_sd,
                    threads: config.threads,
                },
            )?);
        } else {
            mask_stats = mask_seeds::<SeedLoc>(
                config.shapes.get(shape_id as usize),
                &range,
                &mut query_hits,
                &mut reference_hits,
                query.seqs_mut(),
                config.seed_encoding,
                config.seed_complexity_cut,
                config.reduction,
            );
        }

        let ranking = if config.keep_target_id && config.linear_stage1_query {
            SeedLoc::build_ranking(config.kmer_ranking, query, &query_hits, &reference_hits)
        } else {
            None
        };
        let mut searched_partitions = 0usize;
        for relative_partition in 0..relative_count {
            run_stage1(SearchShapePartition {
                shape_id,
                relative_partition,
                absolute_partition: range.begin() + relative_partition as u32,
                query_seed_hits: &query_hits[relative_partition],
                reference_seed_hits: &reference_hits[relative_partition],
                previous_shape_patterns: &patterns[..patterns.len() - 1],
                shape_patterns: &patterns,
                ungapped_raw_score: config.ungapped_raw_score,
                seedp_mask: seedp_mask(config.seedp_bits),
                kmer_ranking: ranking.as_ref(),
            })?;
            searched_partitions += 1;
        }
        report.chunks.push(SearchShapeChunkReport {
            chunk,
            range,
            query_seed_count: query_seeds.size(),
            reference_seed_count: reference_seeds.size(),
            mask_stats,
            frequent_seed_stats,
            searched_partitions,
        });
    }
    Ok(report)
}

/// Matches both C++ `search_shape` overloads, including packed-location dispatch.
pub fn search_shape<F>(
    shape_id: u32,
    query: &mut Block,
    reference: &mut Block,
    config: &SearchShapeConfig<'_>,
    mut run_stage1: F,
) -> Result<SearchShapeReport, String>
where
    F: FnMut(SearchShapePartition<'_>) -> Result<(), String>,
{
    if config.keep_target_id {
        search_shape_with_seed_loc::<PackedLocId, _>(
            shape_id,
            query,
            reference,
            config,
            &mut run_stage1,
        )
    } else {
        search_shape_with_seed_loc::<PackedLoc, _>(
            shape_id,
            query,
            reference,
            config,
            &mut run_stage1,
        )
    }
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
    use crate::data::seed_histogram::{SeedHistogram, SeedPartitionRange};
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

    fn initialize_histogram(
        block: &mut Block,
        shapes: &ShapeConfig,
        reduction: &Reduction,
        seedp_bits: i32,
    ) {
        let cfg = enum_cfg(None);
        let histogram = SeedHistogram::from_block(
            block, false, &NO_FILTER, &cfg, seedp_bits, 1, shapes, reduction, 0, 1,
        )
        .unwrap();
        *block.hst() = histogram;
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

    #[test]
    fn test_search_shape_runs_every_chunk_and_partition_in_order() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let mut query = block(&[0, 1, 2, 3, 4, 5]);
        let mut reference = block(&[9, 0, 1, 2, 8, 7]);
        initialize_histogram(&mut query, &shapes, &reduction, 2);
        initialize_histogram(&mut reference, &shapes, &reduction, 2);
        let mut config = SearchShapeConfig::new(&shapes, &reduction);
        config.seedp_bits = 2;
        config.index_chunks = 3;
        config.threads = 1;
        config.ungapped_raw_score = 17;
        let mut seen = Vec::new();
        let mut joined_bytes = 0usize;

        let report = search_shape(0, &mut query, &mut reference, &config, |partition| {
            seen.push((
                partition.relative_partition,
                partition.absolute_partition,
                partition.previous_shape_patterns.to_vec(),
                partition.shape_patterns.to_vec(),
                partition.ungapped_raw_score,
                partition.seedp_mask,
            ));
            joined_bytes += partition.query_seed_hits.data().len();
            Ok(())
        })
        .unwrap();

        assert_eq!(
            report.seed_loc_kind,
            Some(SearchShapeSeedLocKind::PackedLoc)
        );
        assert_eq!(report.chunks.len(), 3);
        assert_eq!(
            report.chunks[0].range,
            SeedPartitionRange::with_bounds(0, 2)
        );
        assert_eq!(
            report.chunks[1].range,
            SeedPartitionRange::with_bounds(2, 3)
        );
        assert_eq!(
            report.chunks[2].range,
            SeedPartitionRange::with_bounds(3, 4)
        );
        assert_eq!(
            seen.iter().map(|entry| entry.1).collect::<Vec<_>>(),
            vec![0, 1, 2, 3]
        );
        assert!(seen.iter().all(|entry| {
            entry.2.is_empty() && entry.3 == vec![0b111] && entry.4 == 17 && entry.5 == 3
        }));
        assert!(joined_bytes > 0);
    }

    #[test]
    fn test_search_shape_dispatches_packed_loc_id_and_builds_ranking() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let mut query = block(&[0, 1, 2, 3, 4]);
        let mut reference = block(&[0, 1, 2, 7, 8]);
        initialize_histogram(&mut query, &shapes, &reduction, 1);
        initialize_histogram(&mut reference, &shapes, &reduction, 1);
        let mut config = SearchShapeConfig::new(&shapes, &reduction);
        config.seedp_bits = 1;
        config.keep_target_id = true;
        config.linear_stage1_query = true;
        config.kmer_ranking = true;
        let mut ranking_seen = false;

        let report = search_shape(0, &mut query, &mut reference, &config, |partition| {
            let ranking = partition.kmer_ranking.unwrap();
            assert_eq!(ranking.rank.len(), 1);
            ranking_seen = true;
            Ok(())
        })
        .unwrap();

        assert!(ranking_seen);
        assert_eq!(
            report.seed_loc_kind,
            Some(SearchShapeSeedLocKind::PackedLocId)
        );
        assert_eq!(report.chunks[0].searched_partitions, 2);
    }

    #[test]
    fn test_search_shape_validates_shape_and_thread_count_before_work() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let mut query = block(&[0, 1, 2]);
        let mut reference = block(&[0, 1, 2]);
        let mut config = SearchShapeConfig::new(&shapes, &reduction);
        assert_eq!(
            search_shape(1, &mut query, &mut reference, &config, |_| Ok(())).unwrap_err(),
            "Shape id out of range."
        );
        config.threads = 0;
        assert_eq!(
            search_shape(0, &mut query, &mut reference, &config, |_| Ok(())).unwrap_err(),
            "Stage 0 thread count must be positive."
        );
    }
}
