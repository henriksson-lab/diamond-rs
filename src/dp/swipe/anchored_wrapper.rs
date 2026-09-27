//! Anchored SWIPE orchestration translated from
//! `diamond/src/dp/swipe/anchored_wrapper.cpp`.
//!
//! SIMD lanes and the C++ thread pool map to the scalar anchored kernel and
//! deterministic task ranges. The range boundaries and statistics retain the
//! original 16-target batching and `swipe_task_size` threshold semantics.

use crate::align::hsp::Hsp;
use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::value::Letter;
use crate::config::Sensitivity;
use crate::dp::anchored::{smith_waterman_simd, Stats, Target};
use crate::dp::score_profile::{make_profile16, LongScoreProfile};
use crate::dp::swipe::{self, DpTarget, Params, Targets};
use crate::stats::cbs::{
    compute_composition, CbsMode, MatrixAdjustRule, TargetMatrix, TargetMatrixAdjustment,
};
use crate::stats::{self, score_matrix::ScoreMatrix};
use crate::util::geo;
use crate::util::interval::Interval;
use std::sync::{Arc, Mutex};
use std::time::Instant;

const DEFAULT_SWIPE_TASK_SIZE: i64 = 100_000_000;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TargetVector {
    pub int16: Vec<Target>,
}

#[derive(Debug, Clone)]
pub struct Profiles {
    pub int16: LongScoreProfile<i16>,
}

impl Profiles {
    pub fn new(seq: &[Letter], cbs: Option<&[i8]>, padding: usize, matrix: &ScoreMatrix) -> Self {
        Self {
            int16: make_profile16(seq, cbs, padding, matrix),
        }
    }

    pub fn reverse(&self) -> Self {
        Self {
            int16: self.int16.reverse(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WrapperConfig {
    pub query_len: i32,
    pub sensitivity: Sensitivity,
    pub score_hint: i32,
}

#[derive(Clone)]
pub struct AnchoredSwipeConfig<'a> {
    pub query: &'a [Letter],
    pub query_cbs: Option<&'a [i8]>,
    pub score_hint: i32,
    pub sensitivity: Sensitivity,
    pub recompute_adjusted: bool,
    pub target_profiles: bool,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub max_evalue: f64,
    pub score_matrix: &'a ScoreMatrix,
    pub statistics: Option<Arc<Mutex<Statistics>>>,
    /// Explicit replacement for C++'s process-global `config.swipe_task_size`.
    pub swipe_task_size: i64,
    /// Scale used by the composition-adjusted recomputation matrix.
    pub cbs_matrix_scale: i32,
}

impl<'a> AnchoredSwipeConfig<'a> {
    pub fn new(query: &'a [Letter], score_matrix: &'a ScoreMatrix) -> Self {
        Self {
            query,
            query_cbs: None,
            score_hint: 0,
            sensitivity: Sensitivity::Default,
            recompute_adjusted: false,
            target_profiles: false,
            query_cover: 0.0,
            subject_cover: 0.0,
            query_or_target_cover: 0.0,
            max_evalue: f64::MAX,
            score_matrix,
            statistics: None,
            swipe_task_size: DEFAULT_SWIPE_TASK_SIZE,
            cbs_matrix_scale: 1,
        }
    }
}

fn inc(cfg: &AnchoredSwipeConfig<'_>, value: StatValue, n: i64) {
    if let Some(statistics) = &cfg.statistics {
        statistics.lock().unwrap().inc(value, n);
    }
}

fn micros(start: Instant) -> i64 {
    start.elapsed().as_micros().min(i64::MAX as u128) as i64
}

pub fn get_band(_query_len: i32, sensitivity: Sensitivity) -> i32 {
    if sensitivity >= Sensitivity::UltraSensitive {
        160
    } else if sensitivity >= Sensitivity::MoreSensitive {
        96
    } else {
        32
    }
}

#[allow(clippy::too_many_arguments)]
pub fn align_right(
    target_seq: &[Letter],
    reverse: bool,
    i: i32,
    j: i32,
    d_begin: i32,
    d_end: i32,
    _prefix_score: i32,
    targets: &mut TargetVector,
    target_idx: i64,
    cfg: &WrapperConfig,
) {
    align_right_with_matrix(
        target_seq, reverse, i, j, d_begin, d_end, targets, target_idx, cfg, None, 1,
    );
}

#[allow(clippy::too_many_arguments)]
fn align_right_with_matrix(
    target_seq: &[Letter],
    reverse: bool,
    i: i32,
    j: i32,
    mut d_begin: i32,
    mut d_end: i32,
    targets: &mut TargetVector,
    target_idx: i64,
    cfg: &WrapperConfig,
    matrix: Option<Arc<TargetMatrix>>,
    matrix_scale: i32,
) {
    let query_len = cfg.query_len - i;
    let mut target_len = target_seq.len() as i32;
    let band = get_band(cfg.query_len, cfg.sensitivity);
    d_begin -= band;
    d_end += band - 1;
    let d0 = geo::clip_diag(geo::diag_sub_matrix(d_begin, i, j), query_len, target_len);
    let d1 = geo::clip_diag(geo::diag_sub_matrix(d_end, i, j), query_len, target_len);
    target_len = target_len.min(geo::j(query_len - 1, d0) + 1);
    assert!(target_len > 0);
    assert!(d1 >= d0);
    let clipped = if reverse {
        target_seq[target_seq.len() - target_len as usize..].to_vec()
    } else {
        target_seq[..target_len as usize].to_vec()
    };
    let mut target = Target::new(clipped, d0, d1 + 1, i, query_len, target_idx, reverse);
    if let Some(matrix) = matrix {
        target = target.with_matrix(matrix, matrix_scale);
    }
    targets.int16.push(target);
}

#[allow(clippy::too_many_arguments)]
pub fn align_left(
    target_seq: &[Letter],
    i: i32,
    j: i32,
    d_begin: i32,
    d_end: i32,
    suffix_score: i32,
    targets: &mut TargetVector,
    target_idx: i64,
    cfg: &WrapperConfig,
) {
    let query_len = cfg.query_len;
    let target_len = target_seq.len() as i32;
    align_right(
        &target_seq[..=j as usize],
        true,
        query_len - 1 - i,
        target_len - 1 - j,
        geo::rev_diag(d_end, query_len, target_len),
        geo::rev_diag(d_begin, query_len, target_len),
        suffix_score,
        targets,
        target_idx,
        cfg,
    );
}

pub fn add_target(
    target: &DpTarget,
    targets: &mut TargetVector,
    target_idx: &mut i64,
    cfg: &WrapperConfig,
) {
    if target.extend_right(cfg.query_len) {
        let i = target.anchor.query_end;
        let j = target.anchor.subject_end;
        align_right_with_matrix(
            &target.seq[j as usize..],
            false,
            i,
            j,
            target.anchor.d_min_right,
            target.anchor.d_max_right,
            targets,
            *target_idx,
            cfg,
            target.matrix.clone(),
            target.matrix_scale,
        );
        *target_idx += 1;
    }
    if target.extend_left() {
        let query_len = cfg.query_len;
        let target_len = target.seq.len() as i32;
        let i = target.anchor.query_begin - 1;
        let j = target.anchor.subject_begin - 1;
        let ir = query_len - 1 - i;
        let jr = target_len - 1 - j;
        align_right_with_matrix(
            &target.seq[..=j as usize],
            true,
            ir,
            jr,
            geo::rev_diag(target.anchor.d_max_left, query_len, target_len),
            geo::rev_diag(target.anchor.d_min_left, query_len, target_len),
            targets,
            *target_idx,
            cfg,
            target.matrix.clone(),
            target.matrix_scale,
        );
        *target_idx += 1;
    }
}

/// C++ `swipe_threads`, executed deterministically with identical task ranges.
pub fn swipe_threads(
    targets: &mut [Target],
    query: &[Letter],
    cfg: &AnchoredSwipeConfig<'_>,
) -> Stats {
    let threshold = cfg.swipe_task_size.max(1);
    let mut ranges = Vec::new();
    let mut begin = 0usize;
    let mut end = 0usize;
    let mut size = 0i64;
    while end < targets.len() {
        let next = (end + 16).min(targets.len());
        size += targets[end..next]
            .iter()
            .map(Target::gross_cells)
            .sum::<i64>();
        end = next;
        if size >= threshold {
            ranges.push((begin, end));
            begin = end;
            size = 0;
        }
    }
    let asynchronous = !ranges.is_empty();
    if !asynchronous {
        ranges.push((begin, end));
    } else if begin < end {
        ranges.push((begin, end));
    }
    let mut out = Stats::default();
    for (begin, end) in ranges.iter().copied() {
        let stats = smith_waterman_simd(query, &mut targets[begin..end], cfg.score_matrix);
        out.gross_cells += stats.gross_cells;
        out.net_cells += stats.net_cells;
    }
    inc(cfg, StatValue::SwipeTasksTotal, ranges.len() as i64);
    if asynchronous {
        inc(cfg, StatValue::SwipeTasksAsync, ranges.len() as i64);
    }
    out
}

pub fn select_matrix<'a>(_query_len: i32, matrix: &'a ScoreMatrix) -> &'a ScoreMatrix {
    matrix
}

pub fn anchored_swipe(targets: &mut Targets, cfg: &AnchoredSwipeConfig<'_>) -> Vec<Hsp> {
    let total = Instant::now();
    let timer = Instant::now();
    let mut target_count = 0i64;
    let mut max_target_len = 0usize;
    for bin in targets.iter() {
        target_count += bin.size();
        if !cfg.target_profiles {
            max_target_len = max_target_len.max(
                bin.as_slice()
                    .iter()
                    .map(|target| target.seq.len())
                    .max()
                    .unwrap_or(0),
            );
        }
    }
    let mut target_vec = TargetVector {
        int16: Vec::with_capacity((target_count * 2).max(0) as usize),
    };
    inc(cfg, StatValue::TimeAnchoredSwipeAlloc, micros(timer));

    let timer = Instant::now();
    let profiles = (!cfg.target_profiles).then(|| {
        Profiles::new(
            cfg.query,
            cfg.query_cbs,
            cfg.query.len() + max_target_len + 32,
            select_matrix(cfg.query.len() as i32, cfg.score_matrix),
        )
    });
    let _profiles_reverse = profiles.as_ref().map(Profiles::reverse);
    inc(cfg, StatValue::TimeProfile, micros(timer));

    let timer = Instant::now();
    let wrapper_cfg = WrapperConfig {
        query_len: cfg.query.len() as i32,
        sensitivity: cfg.sensitivity,
        score_hint: cfg.score_hint,
    };
    let mut target_idx = 0;
    for bin in targets.iter() {
        for target in bin.as_slice() {
            add_target(target, &mut target_vec, &mut target_idx, &wrapper_cfg);
        }
    }
    inc(cfg, StatValue::TimeAnchoredSwipeAdd, micros(timer));

    let timer = Instant::now();
    target_vec.int16.sort_by_key(Target::band);
    inc(cfg, StatValue::TimeAnchoredSwipeSort, micros(timer));

    let timer = Instant::now();
    let cell_stats = swipe_threads(&mut target_vec.int16, cfg.query, cfg);
    inc(cfg, StatValue::GrossDpCells, cell_stats.gross_cells);
    inc(cfg, StatValue::NetDpCells, cell_stats.net_cells);
    inc(cfg, StatValue::TimeSw, micros(timer));

    let timer = Instant::now();
    target_vec.int16.sort_by_key(|target| target.target_idx);
    inc(cfg, StatValue::TimeAnchoredSwipeSort, micros(timer));

    let timer = Instant::now();
    let mut target_it = target_vec.int16.iter();
    let mut out = Vec::new();
    let mut recompute = swipe::targets();
    let query_composition = compute_composition(cfg.query);
    for bin in 0..swipe::BINS {
        for target in targets[bin].as_slice() {
            if target.anchor.score == 0 {
                continue;
            }
            inc(cfg, StatValue::Ext16, 1);
            let mut score = target.anchor.score;
            let mut i0 = target.anchor.query_begin;
            let mut i1 = target.anchor.query_end;
            let mut j0 = target.anchor.subject_begin;
            let mut j1 = target.anchor.subject_end;
            if target.extend_right(cfg.query.len() as i32) {
                let right = target_it.next().expect("missing right anchored extension");
                score += right.score;
                i1 += right.query_end;
                j1 += right.target_end;
            }
            if target.extend_left() {
                let left = target_it.next().expect("missing left anchored extension");
                score += left.score;
                i0 -= left.query_end;
                j0 -= left.target_end;
            }
            let query_cover = (i1 - i0) as f64 / cfg.query.len() as f64 * 100.0;
            let target_cover = (j1 - j0) as f64 / target.seq.len() as f64 * 100.0;
            let coverage_filtered = (cfg.query_or_target_cover > 0.0
                || cfg.query_cover > 0.0
                || cfg.subject_cover > 0.0)
                && (query_cover.max(target_cover) < cfg.query_or_target_cover
                    || query_cover < cfg.query_cover
                    || target_cover < cfg.subject_cover);
            if !cfg.recompute_adjusted && coverage_filtered {
                continue;
            }
            let evalue =
                cfg.score_matrix
                    .evalue(score, cfg.query.len() as u32, target.seq.len() as u32);
            if !cfg.recompute_adjusted && evalue > cfg.max_evalue {
                continue;
            }
            if cfg.recompute_adjusted && (query_cover < 70.0 || target_cover < 70.0) {
                continue;
            }
            let approx_id = stats::approx_id(score, i1 - i0, j1 - j0);
            if cfg.recompute_adjusted && approx_id < 70.0 {
                let matrix = if let (Some(joint_probs), Some(freq_ratios)) = (
                    cfg.score_matrix.joint_probs(),
                    cfg.score_matrix.freq_ratios(),
                ) {
                    let mut local_stats = Statistics::new();
                    let matrix = TargetMatrix::from_composition_adjustment(
                        &query_composition,
                        cfg.query.len() as i32,
                        CbsMode::MatrixAdjust,
                        &target.seq,
                        &mut local_stats,
                        cfg.score_matrix,
                        MatrixAdjustRule::UserSpecifiedRelEntropy,
                        TargetMatrixAdjustment {
                            matrix_scale: cfg.cbs_matrix_scale.max(1),
                            joint_probs,
                            background_freqs: cfg.score_matrix.background_freqs(),
                            freq_ratios: Some(freq_ratios),
                            tolerance: 1.0e-5,
                            max_iterations: 17,
                        },
                    )
                    .expect("BLOSUM composition adjustment inputs are complete");
                    if let Some(statistics) = &cfg.statistics {
                        *statistics.lock().unwrap() += &local_stats;
                    }
                    Some(Arc::new(matrix))
                } else {
                    None
                };
                let mut dp_target = DpTarget::new(
                    target.seq.clone(),
                    target.true_target_len,
                    target.d_begin,
                    target.d_end,
                    target.target_idx,
                    cfg.query.len() as i32,
                    Default::default(),
                    target.anchor,
                );
                if let Some(matrix) = matrix {
                    dp_target = dp_target.with_matrix(matrix, cfg.cbs_matrix_scale.max(1));
                }
                recompute[bin].push_back(dp_target);
                inc(cfg, StatValue::ExtensionsRecompute, 1);
            } else if !coverage_filtered {
                let mut hsp = Hsp::new();
                hsp.score = score;
                hsp.evalue = evalue;
                hsp.bit_score = cfg.score_matrix.bitscore(score as f64);
                hsp.swipe_target = target.target_idx as i32;
                hsp.query_range = Interval::new(i0, i1);
                hsp.query_source_range = hsp.query_range;
                hsp.subject_range = Interval::new(j0, j1);
                hsp.subject_source_range = hsp.subject_range;
                hsp.target_seq = target.seq.clone();
                hsp.approx_id = hsp.approx_id_percent(cfg.query, &target.seq);
                out.push(hsp);
            }
        }
    }
    inc(cfg, StatValue::TimeAnchoredSwipeOutput, micros(timer));
    inc(cfg, StatValue::TimeAnchoredSwipe, micros(total));
    if cfg.recompute_adjusted {
        let mut params = Params::new(cfg.query, cfg.score_matrix);
        params.v = swipe::HspValues::COORDS;
        params.cbs_matrix_scale = cfg.cbs_matrix_scale.max(1);
        params.statistics = cfg.statistics.clone();
        out.extend(swipe::swipe(&recompute, &mut params));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dp::anchored::smith_waterman;
    use crate::dp::swipe::Anchor;

    fn matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn bands_match_sensitivity_thresholds() {
        assert_eq!(get_band(100, Sensitivity::Default), 32);
        assert_eq!(get_band(100, Sensitivity::MoreSensitive), 96);
        assert_eq!(get_band(100, Sensitivity::UltraSensitive), 160);
    }

    #[test]
    fn profiles_reverse_and_matrix_selection_preserve_inputs() {
        let score_matrix = matrix();
        let sequence = vec![0, 1, 2, 3];
        let profiles = Profiles::new(&sequence, None, sequence.len() + 32, &score_matrix);
        assert_eq!(profiles.int16.length(), sequence.len());
        assert_eq!(profiles.reverse().int16.length(), sequence.len());
        assert!(std::ptr::eq(
            select_matrix(sequence.len() as i32, &score_matrix),
            &score_matrix
        ));
    }

    #[test]
    fn task_ranges_follow_sixteen_lane_thresholds() {
        let score_matrix = matrix();
        let query = vec![0; 8];
        let statistics = Arc::new(Mutex::new(Statistics::new()));
        let mut cfg = AnchoredSwipeConfig::new(&query, &score_matrix);
        cfg.statistics = Some(statistics.clone());
        cfg.swipe_task_size = 1;
        let mut targets = (0..17)
            .map(|idx| Target::new(vec![0], 0, 1, 0, 1, idx, false))
            .collect::<Vec<_>>();
        swipe_threads(&mut targets, &query, &cfg);
        let statistics = statistics.lock().unwrap();
        assert_eq!(statistics.get(StatValue::SwipeTasksTotal), 2);
        assert_eq!(statistics.get(StatValue::SwipeTasksAsync), 2);
    }

    #[test]
    fn target_profile_scores_are_used() {
        let score_matrix = matrix();
        let mut scores = vec![-5i8; 26 * 32];
        scores[0] = 20;
        let adjusted = Arc::new(TargetMatrix::new(scores, -5, 20));
        let mut targets = vec![Target::new(vec![0], 0, 1, 0, 1, 0, false).with_matrix(adjusted, 1)];
        smith_waterman(&[0], &mut targets, &score_matrix);
        assert_eq!(targets[0].score, 21);
    }

    #[test]
    fn add_target_preserves_target_matrix() {
        let adjusted = Arc::new(TargetMatrix::new(vec![1; 26 * 32], 1, 1));
        let target = DpTarget {
            seq: vec![0; 10].into(),
            d_begin: -1,
            d_end: 2,
            cols: 10,
            true_target_len: 10,
            target_idx: 3,
            carry_over: Default::default(),
            anchor: Anchor {
                query_begin: 4,
                query_end: 7,
                subject_begin: 4,
                subject_end: 7,
                d_min_left: -1,
                d_max_left: 2,
                d_min_right: -1,
                d_max_right: 2,
                prefix_score: 12,
                score: 20,
            },
            matrix: Some(adjusted.clone()),
            matrix_scale: 2,
        };
        let cfg = WrapperConfig {
            query_len: 20,
            sensitivity: Sensitivity::Default,
            score_hint: 30,
        };
        let mut out = TargetVector::default();
        let mut index = 0;
        add_target(&target, &mut out, &mut index, &cfg);
        assert_eq!(index, 2);
        assert!(out.int16.iter().all(|target| {
            target
                .matrix
                .as_ref()
                .is_some_and(|matrix| Arc::ptr_eq(matrix, &adjusted))
                && target.matrix_scale == 2
        }));
    }

    #[test]
    fn anchored_output_counts_only_nonzero_anchors() {
        let score_matrix = matrix();
        let query = vec![0, 1, 2];
        let statistics = Arc::new(Mutex::new(Statistics::new()));
        let mut targets = swipe::targets();
        for (target_idx, score) in [(17, 30), (18, 0)] {
            let mut target = DpTarget::new(
                query.clone(),
                query.len() as i32,
                -1,
                2,
                target_idx,
                query.len() as i32,
                Default::default(),
                Default::default(),
            );
            target.anchor = Anchor {
                query_begin: 0,
                query_end: 3,
                subject_begin: 0,
                subject_end: 3,
                score,
                ..Anchor::default()
            };
            targets[0].push_back(target);
        }
        let mut cfg = AnchoredSwipeConfig::new(&query, &score_matrix);
        cfg.statistics = Some(statistics.clone());
        let output = anchored_swipe(&mut targets, &cfg);
        assert_eq!(output.len(), 1);
        assert_eq!(output[0].swipe_target, 17);
        assert_eq!(output[0].approx_id, 100.0);
        assert_eq!(statistics.lock().unwrap().get(StatValue::Ext16), 1);
    }
}
