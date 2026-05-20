use crate::align::gapped_filter::SeedHit;
use crate::basic::consts::MAX_CONTEXT;
use crate::basic::statistics::Statistics;
use crate::basic::value::{BlockId, Letter, Loc, TRUE_AA};
use crate::chaining::{hamming_ext, run as chaining_run, HammingExtConfig};
use crate::data::block::Block;
use crate::dp::ungapped::{xdrop_ungapped_with_cbs, DiagonalSegment};
use crate::search::sensitivity::ExtensionMode;
use crate::stats::cbs::{
    adjust_matrix, count_true_aa, CbsMode, CbsThresholds, MatrixAdjustRule, TargetMatrix,
};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::data_structures::FlatArray;
use crate::util::geo::assert_diag_bounds;
use crate::util::hsp::ApproxHsp;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct UngappedStageConfig {
    pub query_contexts: usize,
    pub query_translated: bool,
    pub hamming_ext: HammingExtConfig,
    pub lin_stage1_query: bool,
    pub lin_stage1_target: bool,
    pub comp_based_stats: CbsMode,
    pub anchored_swipe: bool,
    pub log_extend: bool,
    pub no_chaining_merge_hsps: bool,
    pub xdrop: i32,
    pub matrix_adjust_thresholds: CbsThresholds,
}

impl Default for UngappedStageConfig {
    fn default() -> Self {
        Self {
            query_contexts: 1,
            query_translated: false,
            hamming_ext: HammingExtConfig::default(),
            lin_stage1_query: false,
            lin_stage1_target: false,
            comp_based_stats: CbsMode::Disabled,
            anchored_swipe: false,
            log_extend: false,
            no_chaining_merge_hsps: false,
            xdrop: 20,
            matrix_adjust_thresholds: CbsThresholds::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct WorkTarget {
    pub block_id: BlockId,
    pub seq: Vec<Letter>,
    pub ungapped_score: [i32; MAX_CONTEXT as usize],
    pub hsp: [Vec<ApproxHsp>; MAX_CONTEXT as usize],
    pub matrix: Option<Arc<TargetMatrix>>,
    pub done: bool,
}

impl WorkTarget {
    /// Matches C++ `Extension::WorkTarget::WorkTarget`.
    pub fn new(
        block_id: BlockId,
        seq: &[Letter],
        query: &[Letter],
        query_len_true_aa: Loc,
        query_comp: &[f64; TRUE_AA as usize],
        _max_target_len: Loc,
        _stats: &mut Statistics,
        score_matrix: &ScoreMatrix,
        cfg: &UngappedStageConfig,
    ) -> Self {
        let mut target = WorkTarget {
            block_id,
            seq: seq.to_vec(),
            ungapped_score: [0; MAX_CONTEXT as usize],
            hsp: std::array::from_fn(|_| Vec::new()),
            matrix: None,
            done: false,
        };
        if !cfg.anchored_swipe {
            let rule = adjust_matrix(
                query_comp,
                query_len_true_aa,
                cfg.comp_based_stats,
                &target.seq,
                score_matrix.background_scores(),
                cfg.matrix_adjust_thresholds,
            );
            if rule != MatrixAdjustRule::DontAdjustMatrix {
                let target_comp = crate::stats::cbs::compute_composition(&target.seq);
                target.matrix = Some(Arc::new(TargetMatrix::from_hauser_global(
                    query_comp,
                    &target_comp,
                    score_matrix,
                )));
            }
        }
        let _ = query;
        target
    }
}

/// Matches C++ `Extension::ungapped_stage(begin, end, ..., block_id, ...)`.
pub fn ungapped_stage_target(
    hits: &mut [SeedHit],
    query_seq: &[Vec<Letter>],
    query_cbs: &[Vec<i8>],
    query_comp: &[f64; TRUE_AA as usize],
    block_id: BlockId,
    max_target_len: Loc,
    stats: &mut Statistics,
    targets: &Block,
    mode: ExtensionMode,
    cfg: &UngappedStageConfig,
    score_matrix: &ScoreMatrix,
) -> WorkTarget {
    let ref_seqs = targets.seqs();
    let target_seq = ref_seqs.get(block_id as usize);
    let mut target = WorkTarget::new(
        block_id,
        target_seq,
        &query_seq[0],
        count_true_aa(&query_seq[0]),
        query_comp,
        max_target_len,
        stats,
        score_matrix,
        cfg,
    );

    let with_diag_filter = (cfg.hamming_ext.hamming_ext
        || cfg.hamming_ext.diag_filter_cov.is_some()
        || cfg.hamming_ext.diag_filter_id.is_some())
        && cfg.query_contexts == 1;

    if mode == ExtensionMode::Full {
        for hit in hits.iter() {
            let frame = hit.frame as usize;
            target.ungapped_score[frame] = target.ungapped_score[frame].max(hit.score);
        }
        if !with_diag_filter {
            return target;
        }
    }

    if hits.len() == 1 && cfg.query_translated {
        let hit = hits[0];
        let frame = hit.frame as usize;
        target.ungapped_score[frame] = hit.score;
        target.hsp[frame].push(ApproxHsp::from_parts(
            hit.diag(),
            hit.diag(),
            hit.score,
            hit.frame as i32,
            hit.query_range(),
            hit.target_range(),
            crate::util::hsp::Anchor::from_diagonal_segment(hit.diag_segment()),
            f64::MAX,
        ));
        return target;
    }

    hits.sort();
    let use_hauser = cfg.comp_based_stats.hauser();
    let mut diagonal_segments: [Vec<DiagonalSegment>; MAX_CONTEXT as usize] =
        std::array::from_fn(|_| Vec::new());

    for hit in hits.iter() {
        let frame = hit.frame as usize;
        target.ungapped_score[frame] = target.ungapped_score[frame].max(hit.score);
        if let Some(last) = diagonal_segments[frame].last() {
            if last.diag() == hit.diag() && last.subject_end() >= hit.j {
                continue;
            }
        }
        let cbs = if use_hauser {
            query_cbs.get(frame).map(Vec::as_slice)
        } else {
            None
        };
        let d = xdrop_ungapped_with_cbs(
            &query_seq[frame],
            cbs,
            &target.seq,
            hit.i as usize,
            hit.j as usize,
            cfg.xdrop,
            with_diag_filter,
            score_matrix,
        );
        if d.score > 0 {
            diagonal_segments[frame].push(d);
        }
    }

    if with_diag_filter {
        let h = hamming_ext(
            &mut diagonal_segments[0],
            query_seq[0].len() as i32,
            target.seq.len() as i32,
            !cfg.lin_stage1_target && !cfg.lin_stage1_query,
            &cfg.hamming_ext,
            score_matrix,
        );
        if h.score > 0 {
            target.done = true;
            target.hsp[0].push(h);
            return target;
        }
        if h.score < 0 {
            target.ungapped_score[0] = 0;
            return target;
        }
    }

    if mode == ExtensionMode::Full {
        return target;
    }

    for frame in 0..cfg.query_contexts {
        if diagonal_segments[frame].is_empty() {
            continue;
        }
        diagonal_segments[frame]
            .sort_by(|x, y| x.diag().cmp(&y.diag()).then_with(|| x.j.cmp(&y.j)));
        // Chaining `max_shift` (= C++ `config.chaining_maxgap`, default 2000
        // per `diamond/src/basic/config.cpp:546`). This bounds the diagonal
        // shift between two linkable segments inside `backtrace_old` at
        // `chaining.rs:910`. We previously passed `gap_open + 6` (= 17 on
        // BLOSUM62) which silently broke chaining whenever a true HSP had
        // gaps wider than ~17 diagonals — CD209-family (Q8MIS5) alignments
        // have three large gaps with net ~150-letter subject insertions, so
        // 17 truncates the chain to the C-terminal high-density cluster and
        // produces short alignments (q=102..240 instead of q=1..240).
        const CHAINING_MAXGAP: i32 = 2000;
        let (_score, mut hsps) = chaining_run(
            &query_seq[frame],
            &target.seq,
            &diagonal_segments[frame],
            cfg.log_extend,
            frame as u32,
            CHAINING_MAXGAP,
            score_matrix,
            cfg.no_chaining_merge_hsps,
        );
        hsps.sort_by(|x, y| x.frame.cmp(&y.frame).then_with(|| x.d_min.cmp(&y.d_min)));
        target.hsp[frame] = hsps;
    }
    target
}

/// Matches C++ `Extension::ungapped_stage(const Sequence*, ..., FlatArray<SeedHit>::Iterator, ...)`.
pub fn ungapped_stage_seed_hits(
    query_seq: &[Vec<Letter>],
    query_cbs: &[Vec<i8>],
    query_comp: &[f64; TRUE_AA as usize],
    seed_hits: &FlatArray<SeedHit>,
    target_block_ids: &[BlockId],
    stats: &mut Statistics,
    target_block: &Block,
    mode: ExtensionMode,
    cfg: &UngappedStageConfig,
    score_matrix: &ScoreMatrix,
) -> Vec<WorkTarget> {
    let n = seed_hits.size() as usize;
    let mut targets = Vec::new();
    if n == 0 {
        return targets;
    }
    let max_target_len = 0;
    targets.reserve(n);
    for i in 0..n {
        let mut hits = seed_hits.range(i as u64).to_vec();
        let target = ungapped_stage_target(
            &mut hits,
            query_seq,
            query_cbs,
            query_comp,
            target_block_ids[i],
            max_target_len,
            stats,
            target_block,
            mode,
            cfg,
            score_matrix,
        );
        for hsp in &target.hsp[0] {
            assert_diag_bounds(
                hsp.d_max,
                query_seq[0].len() as i32,
                target.seq.len() as i32,
            );
            assert_diag_bounds(
                hsp.d_min,
                query_seq[0].len() as i32,
                target.seq.len() as i32,
            );
            assert!(hsp.score > 0);
            assert!(hsp.max_diag.segment.score > 0);
        }
        targets.push(target);
    }
    targets
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::block::Block;
    use crate::stats::cbs::{compute_composition, hauser_correction};
    use crate::util::data_structures::FlatArray;

    fn sm() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    fn block_with_sequences(seqs: &[&[Letter]]) -> Block {
        let mut block = Block::new();
        for (i, seq) in seqs.iter().enumerate() {
            block
                .push_back(
                    seq,
                    None,
                    None,
                    i as u64,
                    crate::basic::value::SequenceType::AminoAcid,
                    0,
                    false,
                )
                .unwrap();
        }
        block
    }

    #[test]
    fn test_work_target_full_mode_scores_hits_without_chaining() {
        let score_matrix = sm();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let target = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let block = block_with_sequences(&[&target]);
        let comp = compute_composition(&query);
        let mut hits = vec![
            SeedHit::new(1, 1, 12, 0),
            SeedHit::new(3, 3, 27, 0),
            SeedHit::new(2, 1, 5, 0),
        ];
        let mut stats = Statistics::new();
        let cfg = UngappedStageConfig::default();
        let wt = ungapped_stage_target(
            &mut hits,
            &[query],
            &[],
            &comp,
            0,
            0,
            &mut stats,
            &block,
            ExtensionMode::Full,
            &cfg,
            &score_matrix,
        );
        assert_eq!(wt.block_id, 0);
        assert_eq!(wt.seq, target);
        assert_eq!(wt.ungapped_score[0], 27);
        assert!(wt.hsp[0].is_empty());
        assert!(!wt.done);
    }

    #[test]
    fn test_ungapped_stage_translated_single_hit_seed_hsp() {
        let score_matrix = sm();
        let query = vec![0, 1, 2, 3];
        let target = vec![0, 1, 2, 3];
        let block = block_with_sequences(&[&target]);
        let comp = compute_composition(&query);
        let mut hits = vec![SeedHit::new(2, 1, 19, 2)];
        let mut stats = Statistics::new();
        let cfg = UngappedStageConfig {
            query_contexts: 6,
            query_translated: true,
            ..UngappedStageConfig::default()
        };
        let queries = vec![
            query.clone(),
            query.clone(),
            query.clone(),
            query.clone(),
            query.clone(),
            query,
        ];
        let wt = ungapped_stage_target(
            &mut hits,
            &queries,
            &[],
            &comp,
            0,
            0,
            &mut stats,
            &block,
            ExtensionMode::BandedFast,
            &cfg,
            &score_matrix,
        );
        assert_eq!(wt.ungapped_score[2], 19);
        assert_eq!(wt.hsp[2].len(), 1);
        assert_eq!(wt.hsp[2][0].d_min, 1);
        assert_eq!(wt.hsp[2][0].query_range.begin, 2);
        assert_eq!(wt.hsp[2][0].subject_range.begin, 1);
    }

    #[test]
    fn test_ungapped_stage_chains_diagonal_segments() {
        let score_matrix = sm();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let target = query.clone();
        let block = block_with_sequences(&[&target]);
        let comp = compute_composition(&query);
        let cbs = hauser_correction(&query, &score_matrix);
        let mut hits = vec![SeedHit::new(1, 1, 10, 0), SeedHit::new(6, 6, 12, 0)];
        let mut stats = Statistics::new();
        let cfg = UngappedStageConfig {
            comp_based_stats: CbsMode::Hauser,
            ..UngappedStageConfig::default()
        };
        let wt = ungapped_stage_target(
            &mut hits,
            &[query],
            &[cbs],
            &comp,
            0,
            0,
            &mut stats,
            &block,
            ExtensionMode::BandedFast,
            &cfg,
            &score_matrix,
        );
        assert_eq!(wt.ungapped_score[0], 12);
        assert!(!wt.hsp[0].is_empty());
        assert!(wt.hsp[0][0].score > 0);
    }

    #[test]
    fn test_ungapped_stage_seed_hits_groups_targets() {
        let score_matrix = sm();
        let query = vec![0, 1, 2, 3, 4, 5];
        let block = block_with_sequences(&[&query, &[5, 4, 3, 2, 1, 0]]);
        let comp = compute_composition(&query);
        let mut seed_hits = FlatArray::new();
        seed_hits.push_back(&[SeedHit::new(1, 1, 9, 0), SeedHit::new(2, 2, 15, 0)]);
        seed_hits.push_back(&[SeedHit::new(0, 5, 7, 0)]);
        let mut stats = Statistics::new();
        let cfg = UngappedStageConfig::default();
        let out = ungapped_stage_seed_hits(
            &[query],
            &[],
            &comp,
            &seed_hits,
            &[0, 1],
            &mut stats,
            &block,
            ExtensionMode::Full,
            &cfg,
            &score_matrix,
        );
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].block_id, 0);
        assert_eq!(out[0].ungapped_score[0], 15);
        assert_eq!(out[1].block_id, 1);
        assert_eq!(out[1].ungapped_score[0], 7);
    }
}
