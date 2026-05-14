use crate::basic::reduction::Reduction;
use crate::basic::seed_iterator::{
    ContiguousSeedIterator, HashedSeedIterator, MinimizerIterator, SeedIterator, SketchIterator,
};
use crate::basic::sequence::Sequence;
use crate::basic::shape_config::ShapeConfig;
use crate::basic::value::{Loc, SEED_MASK};
use crate::data::block::Block;
use crate::data::flags::{EnumCfg, NoFilter, SeedEncoding};
use crate::data::sequence_set::SequenceSet;
use crate::masking::MaskingAlgo;
use crate::search::seed_complexity::{seed_is_complex, SeedStats};

pub trait EnumSeedsCallback {
    fn call(&mut self, key: u64, pos: usize, block_id: usize, shape_id: i32) -> bool;
    fn finish(&mut self) {}
}

pub trait EnumSeedsFilter {
    fn contains(&self, key: u64, shape_id: i32) -> bool;
}

impl EnumSeedsFilter for NoFilter {
    fn contains(&self, key: u64, shape_id: i32) -> bool {
        self.contains(key, shape_id as u64)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MinimizerKind {
    Minimizer,
    Sketch,
}

#[derive(Clone, Copy)]
pub struct EnumSeedsContext<'a> {
    pub shapes: &'a ShapeConfig,
    pub reduction: &'a Reduction,
    pub min_query_len: Loc,
    pub query_contexts: usize,
}

pub fn enum_seeds<F, Filter>(
    seqs: &mut SequenceSet,
    f: &mut F,
    begin: u32,
    end: u32,
    filter: &Filter,
    cfg: &EnumCfg<'_>,
    ctx: &EnumSeedsContext<'_>,
) -> SeedStats
where
    F: EnumSeedsCallback,
    Filter: EnumSeedsFilter,
{
    let stats = SeedStats::default();
    for i in begin as usize..end as usize {
        if cfg
            .skip
            .is_some_and(|skip| skip[i / ctx.query_contexts.max(1)])
        {
            continue;
        }
        if ctx.min_query_len > 0
            && seqs.source_length_with_contexts(i, ctx.query_contexts) < ctx.min_query_len
        {
            continue;
        }
        let seq = seqs.get(i);
        let buf = ctx.reduction.reduce_seq(seq);
        for shape_id in cfg.shape_begin..cfg.shape_end {
            let sh = ctx.shapes.get(shape_id as usize);
            if seq.len() < sh.length as usize {
                continue;
            }
            let mut it = SeedIterator::new(&buf, sh);
            let mut j = 0usize;
            while it.good() {
                if let Some(key) = it.get(sh, ctx.reduction) {
                    if filter.contains(key, shape_id) {
                        f.call(key, seqs.position(i, j), i, shape_id);
                    }
                }
                j += 1;
            }
        }
    }
    f.finish();
    stats
}

pub fn enum_seeds_minimizer<F, Filter>(
    seqs: &mut SequenceSet,
    f: &mut F,
    begin: u32,
    end: u32,
    filter: &Filter,
    cfg: &EnumCfg<'_>,
    it_param: Loc,
    kind: MinimizerKind,
    ctx: &EnumSeedsContext<'_>,
) -> SeedStats
where
    F: EnumSeedsCallback,
    Filter: EnumSeedsFilter,
{
    let stats = SeedStats::default();
    for i in begin as usize..end as usize {
        if cfg
            .skip
            .is_some_and(|skip| skip[i / ctx.query_contexts.max(1)])
        {
            continue;
        }
        if ctx.min_query_len > 0
            && seqs.source_length_with_contexts(i, ctx.query_contexts) < ctx.min_query_len
        {
            continue;
        }
        let seq = seqs.get(i);
        let buf = ctx.reduction.reduce_seq(seq);
        for shape_id in cfg.shape_begin..cfg.shape_end {
            let sh = ctx.shapes.get(shape_id as usize);
            if seq.len() < sh.length as usize {
                continue;
            }
            match kind {
                MinimizerKind::Minimizer => {
                    let mut it = MinimizerIterator::new(&buf, sh, it_param, ctx.reduction);
                    while it.good() {
                        let key = it.get();
                        if filter.contains(key, shape_id) {
                            f.call(key, seqs.position(i, it.pos() as usize), i, shape_id);
                        }
                        it.increment();
                    }
                }
                MinimizerKind::Sketch => {
                    let mut it = SketchIterator::new(&buf, sh, it_param, ctx.reduction);
                    while it.good() {
                        let key = it.get();
                        if filter.contains(key, shape_id) {
                            f.call(key, seqs.position(i, it.pos() as usize), i, shape_id);
                        }
                        it.increment();
                    }
                }
            }
        }
    }
    f.finish();
    stats
}

pub fn enum_seeds_hashed<F, const BITS: u64, Filter>(
    seqs: &mut SequenceSet,
    f: &mut F,
    begin: u32,
    end: u32,
    filter: &Filter,
    cfg: &EnumCfg<'_>,
    ctx: &EnumSeedsContext<'_>,
) where
    F: EnumSeedsCallback,
    Filter: EnumSeedsFilter,
{
    for i in begin as usize..end as usize {
        if cfg
            .skip
            .is_some_and(|skip| skip[i / ctx.query_contexts.max(1)])
        {
            continue;
        }
        if ctx.min_query_len > 0
            && seqs.source_length_with_contexts(i, ctx.query_contexts) < ctx.min_query_len
        {
            continue;
        }
        for shape_id in cfg.shape_begin..cfg.shape_end {
            let sh = ctx.shapes.get(shape_id as usize);
            if seqs.length(i) < sh.length {
                continue;
            }
            let len = seqs.length(i);
            let seq = seqs.get(i).to_vec();
            let mut it = HashedSeedIterator::<BITS>::new(&seq, len, sh, ctx.reduction);
            while it.good() {
                let key = it.get();
                if filter.contains(key, shape_id) {
                    let j = it.seq_ptr(sh);
                    let accept = !cfg.filter_low_complexity_seeds
                        || seed_is_complex(&seq[j..], sh, cfg.seed_cut, ctx.reduction);
                    if accept {
                        f.call(key, seqs.position(i, j), i, shape_id);
                    } else if cfg.mask_low_complexity_seeds {
                        seqs.get_mut(i)[j] |= SEED_MASK;
                    }
                }
                it.increment();
            }
        }
    }
    f.finish();
}

pub fn enum_seeds_contiguous<
    F,
    const L: usize,
    const BITS: u64,
    const FILTER_MASKED: bool,
    Filter,
>(
    seqs: &mut SequenceSet,
    f: &mut F,
    begin: u32,
    end: u32,
    filter: &Filter,
    cfg: &EnumCfg<'_>,
    ctx: &EnumSeedsContext<'_>,
) where
    F: EnumSeedsCallback,
    Filter: EnumSeedsFilter,
{
    for i in begin as usize..end as usize {
        if cfg
            .skip
            .is_some_and(|skip| skip[i / ctx.query_contexts.max(1)])
        {
            continue;
        }
        if ctx.min_query_len > 0
            && seqs.source_length_with_contexts(i, ctx.query_contexts) < ctx.min_query_len
        {
            continue;
        }
        let seq = seqs.get(i);
        if seq.len() < ContiguousSeedIterator::<L, BITS, FILTER_MASKED>::length() as usize {
            continue;
        }
        let mut it = ContiguousSeedIterator::<L, BITS, FILTER_MASKED>::new(
            Sequence::new(seq),
            ctx.reduction,
        );
        let mut j = 0usize;
        while it.good() {
            if let Some(key) = it.get() {
                if filter.contains(key, 0) && !f.call(key, seqs.position(i, j), i, 0) {
                    return;
                }
            }
            j += 1;
        }
    }
    f.finish();
}

pub fn enum_seeds_worker<F, Filter, const FILTER_MASKED: bool>(
    f: &mut F,
    seqs: &mut SequenceSet,
    begin: u32,
    end: u32,
    filter: &Filter,
    stats: &mut SeedStats,
    cfg: &EnumCfg<'_>,
    ctx: &EnumSeedsContext<'_>,
) -> Result<(), String>
where
    F: EnumSeedsCallback,
    Filter: EnumSeedsFilter,
{
    static ERRMSG: &str = "Unsupported contiguous seed.";
    if cfg.code == SeedEncoding::Contiguous {
        let b = ctx.reduction.bit_size() as u64;
        let l = ctx.shapes.get(cfg.shape_begin as usize).length;
        match (l, b) {
            (7, 4) => enum_seeds_contiguous::<F, 7, 4, FILTER_MASKED, Filter>(
                seqs, f, begin, end, filter, cfg, ctx,
            ),
            (6, 4) => enum_seeds_contiguous::<F, 6, 4, FILTER_MASKED, Filter>(
                seqs, f, begin, end, filter, cfg, ctx,
            ),
            (5, 4) => enum_seeds_contiguous::<F, 5, 4, FILTER_MASKED, Filter>(
                seqs, f, begin, end, filter, cfg, ctx,
            ),
            _ => return Err(ERRMSG.to_string()),
        }
    } else if cfg.code == SeedEncoding::Hashed {
        match ctx.reduction.bit_size() {
            4 => enum_seeds_hashed::<F, 4, Filter>(seqs, f, begin, end, filter, cfg, ctx),
            _ => return Err("Unsupported reduction.".to_string()),
        }
    } else if cfg.minimizer_window > 0 {
        *stats = enum_seeds_minimizer(
            seqs,
            f,
            begin,
            end,
            filter,
            cfg,
            cfg.minimizer_window,
            MinimizerKind::Minimizer,
            ctx,
        );
    } else if cfg.sketch_size > 0 {
        *stats = enum_seeds_minimizer(
            seqs,
            f,
            begin,
            end,
            filter,
            cfg,
            cfg.sketch_size,
            MinimizerKind::Sketch,
            ctx,
        );
    } else {
        *stats = enum_seeds(seqs, f, begin, end, filter, cfg, ctx);
    }
    Ok(())
}

pub fn enum_seeds_block<F, Filter>(
    seqs: &mut Block,
    f: &mut [F],
    filter: &Filter,
    cfg: &EnumCfg<'_>,
    ctx: &EnumSeedsContext<'_>,
) -> Result<SeedStats, String>
where
    F: EnumSeedsCallback,
    Filter: EnumSeedsFilter,
{
    if cfg.soft_masking != MaskingAlgo::None {
        seqs.soft_mask(cfg.soft_masking);
    }

    let partition = cfg
        .partition
        .ok_or_else(|| "EnumCfg::partition is required.".to_string())?;
    if partition.len() < f.len() + 1 {
        return Err("EnumCfg::partition does not cover all callbacks.".to_string());
    }

    let mut stats = vec![SeedStats::default(); f.len()];
    for i in 0..f.len() {
        let begin = partition[i];
        let end = partition[i + 1];
        if cfg.filter_masked_seeds {
            enum_seeds_worker::<F, Filter, true>(
                &mut f[i],
                seqs.seqs_mut(),
                begin,
                end,
                filter,
                &mut stats[i],
                cfg,
                ctx,
            )?;
        } else {
            enum_seeds_worker::<F, Filter, false>(
                &mut f[i],
                seqs.seqs_mut(),
                begin,
                end,
                filter,
                &mut stats[i],
                cfg,
                ctx,
            )?;
        }
    }

    for i in 1..stats.len() {
        stats[0].good_seed_positions += stats[i].good_seed_positions;
        stats[0].low_complexity_seeds += stats[i].low_complexity_seeds;
    }

    if cfg.soft_masking != MaskingAlgo::None {
        let mut l = 0;
        for i in cfg.shape_begin..cfg.shape_end {
            l = l.max(ctx.shapes.get(i as usize).length);
        }
        seqs.remove_soft_masking(l, cfg.mask_seeds);
    }

    Ok(stats.into_iter().next().unwrap_or_default())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SequenceType;
    use crate::data::flags::NO_FILTER;

    #[derive(Default)]
    struct CollectSeeds {
        rows: Vec<(u64, usize, usize, i32)>,
        finished: bool,
        stop_after: Option<usize>,
    }

    impl EnumSeedsCallback for CollectSeeds {
        fn call(&mut self, key: u64, pos: usize, block_id: usize, shape_id: i32) -> bool {
            self.rows.push((key, pos, block_id, shape_id));
            self.stop_after.is_none_or(|n| self.rows.len() < n)
        }

        fn finish(&mut self) {
            self.finished = true;
        }
    }

    struct KeyFilter(u64);

    impl EnumSeedsFilter for KeyFilter {
        fn contains(&self, key: u64, _shape_id: i32) -> bool {
            key == self.0
        }
    }

    fn cfg<'a>(partition: Option<&'a Vec<u32>>, code: SeedEncoding, shape_end: i32) -> EnumCfg<'a> {
        EnumCfg {
            partition,
            shape_begin: 0,
            shape_end,
            code,
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

    #[test]
    fn test_enum_seeds_spaced_filter_and_finish() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["101".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 1, 2, 3]);
        seqs.push(&[4, 5, 6]);
        let mut cb = CollectSeeds::default();
        let stats = enum_seeds(
            &mut seqs,
            &mut cb,
            0,
            2,
            &NO_FILTER,
            &cfg(None, SeedEncoding::SpacedFactor, 1),
            &ctx,
        );

        assert_eq!(stats, SeedStats::default());
        assert!(cb.finished);
        assert_eq!(cb.rows, vec![(2, 1, 0, 0), (12, 2, 0, 0), (32, 6, 1, 0)]);
    }

    #[test]
    fn test_enum_seeds_minimizer_and_sketch() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["11".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 1, 2, 3, 4]);
        let base = cfg(None, SeedEncoding::SpacedFactor, 1);
        let mut minimizer = CollectSeeds::default();
        enum_seeds_minimizer(
            &mut seqs,
            &mut minimizer,
            0,
            1,
            &NO_FILTER,
            &base,
            3,
            MinimizerKind::Minimizer,
            &ctx,
        );
        assert!(minimizer.finished);
        assert!(!minimizer.rows.is_empty());

        let mut sketch = CollectSeeds::default();
        enum_seeds_minimizer(
            &mut seqs,
            &mut sketch,
            0,
            1,
            &NO_FILTER,
            &base,
            2,
            MinimizerKind::Sketch,
            &ctx,
        );
        assert_eq!(sketch.rows.len(), 2);
    }

    #[test]
    fn test_enum_seeds_hashed_masks_low_complexity_seed() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 0, 0, 1]);
        let mut c = cfg(None, SeedEncoding::Hashed, 1);
        c.filter_low_complexity_seeds = true;
        c.mask_low_complexity_seeds = true;
        c.seed_cut = 1.0;
        let mut cb = CollectSeeds::default();
        enum_seeds_hashed::<_, 4, _>(&mut seqs, &mut cb, 0, 1, &NO_FILTER, &c, &ctx);
        assert_eq!(cb.rows.len(), 1);
        assert_ne!(seqs.get(0)[0] & SEED_MASK, 0);
    }

    #[test]
    fn test_enum_seeds_contiguous_stop_and_filter_masked() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["11".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 1, 2, 3]);
        let mut cb = CollectSeeds {
            stop_after: Some(1),
            ..Default::default()
        };
        enum_seeds_contiguous::<_, 2, 4, false, _>(
            &mut seqs,
            &mut cb,
            0,
            1,
            &NO_FILTER,
            &cfg(None, SeedEncoding::Contiguous, 1),
            &ctx,
        );
        assert_eq!(cb.rows.len(), 1);
        assert!(!cb.finished);
    }

    #[test]
    fn test_enum_seeds_worker_dispatch_and_errors() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["11111".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 1, 2, 3, 4, 5]);
        let mut stats = SeedStats::default();
        let mut cb = CollectSeeds::default();
        enum_seeds_worker::<_, _, false>(
            &mut cb,
            &mut seqs,
            0,
            1,
            &NO_FILTER,
            &mut stats,
            &cfg(None, SeedEncoding::Contiguous, 1),
            &ctx,
        )
        .unwrap();
        assert!(cb.finished);
        assert_eq!(stats, SeedStats::default());

        let bad_shapes = ShapeConfig::from_codes(&["1111".to_string()], 0, &reduction).unwrap();
        let bad_ctx = EnumSeedsContext {
            shapes: &bad_shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let err = enum_seeds_worker::<_, _, false>(
            &mut cb,
            &mut seqs,
            0,
            1,
            &NO_FILTER,
            &mut stats,
            &cfg(None, SeedEncoding::Contiguous, 1),
            &bad_ctx,
        )
        .unwrap_err();
        assert_eq!(err, "Unsupported contiguous seed.");
    }

    #[test]
    fn test_enum_seeds_block_partitions_and_soft_mask_restore() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["11111".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2, 3, 4],
                Some("id0"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        block
            .push_back(
                &[5, 6, 7, 8, 9],
                Some("id1"),
                None,
                1,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let partition = vec![0, 1, 2];
        let c = cfg(Some(&partition), SeedEncoding::Contiguous, 1);
        let mut callbacks = [CollectSeeds::default(), CollectSeeds::default()];
        let stats = enum_seeds_block(&mut block, &mut callbacks, &NO_FILTER, &c, &ctx).unwrap();
        assert_eq!(stats, SeedStats::default());
        assert_eq!(callbacks[0].rows.len(), 1);
        assert_eq!(callbacks[1].rows.len(), 1);
    }

    #[test]
    fn test_enum_seeds_filter_skip_and_min_len() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["11".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 3,
            query_contexts: 1,
        };
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 1]);
        seqs.push(&[0, 1, 2]);
        let mut cb = CollectSeeds::default();
        let skip = vec![false, false];
        let mut c = cfg(None, SeedEncoding::SpacedFactor, 1);
        c.skip = Some(&skip);
        enum_seeds(&mut seqs, &mut cb, 0, 2, &KeyFilter(1), &c, &ctx);
        assert_eq!(cb.rows, vec![(1, 4, 1, 0)]);
    }
}
