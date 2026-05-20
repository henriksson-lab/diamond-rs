//! Scalar anchored SWIPE scoring.
//!
//! This ports the score-only recurrence from `dp/swipe/anchored.h` into a
//! straightforward scalar implementation. The C++ version batches targets into
//! SIMD lanes and consumes subject sequence in blocks; this version keeps the
//! same anchored band semantics per target.

use crate::align::hsp::Hsp;
use crate::basic::value::{Letter, LETTER_MASK};
use crate::config::Sensitivity;
use crate::dp::score_profile::{make_profile16, LongScoreProfile};
use crate::dp::swipe::{self, DpTarget, Params, Targets};
use crate::stats::{self, score_matrix::ScoreMatrix};
use crate::util::geo;
use crate::util::interval::Interval;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Stats {
    pub gross_cells: i64,
    pub net_cells: i64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Target {
    pub seq: Vec<Letter>,
    pub d_begin: i32,
    pub d_end: i32,
    pub query_start: i32,
    pub query_length: i32,
    pub target_idx: i64,
    pub reverse: bool,
    pub score: i32,
    pub query_end: i32,
    pub target_end: i32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TargetVector {
    pub int16: Vec<Target>,
}

impl Target {
    pub fn new(
        seq: Vec<Letter>,
        d_begin: i32,
        d_end: i32,
        query_start: i32,
        query_length: i32,
        target_idx: i64,
        reverse: bool,
    ) -> Self {
        Target {
            seq,
            d_begin,
            d_end,
            query_start,
            query_length,
            target_idx,
            reverse,
            score: 0,
            query_end: 0,
            target_end: 0,
        }
    }

    pub fn blank(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn reset(&mut self) {
        self.seq.clear();
    }

    pub fn band(&self) -> i32 {
        self.d_end - self.d_begin
    }

    pub fn cells(&self) -> (i64, i64) {
        let mut net = 0i64;
        let mut gross = 0i64;
        for j in 0..self.seq.len() as i32 {
            let i0 = (self.d_begin + j).max(0);
            let i1 = (self.d_end + j).min(self.query_length);
            net += (i1 - i0).max(0) as i64;
            gross += self.band() as i64;
        }
        (gross, net)
    }

    pub fn gross_cells(&self) -> i64 {
        self.band() as i64 * self.seq.len() as i64
    }
}

pub fn limits(targets: &[Target]) -> (i32, i32) {
    let mut band = 0;
    let mut target_len = 0;
    for target in targets {
        band = band.max(target.band());
        target_len = target_len.max(target.seq.len() as i32);
    }
    (band, target_len)
}

pub fn smith_waterman(
    query: &[Letter],
    targets: &mut [Target],
    score_matrix: &ScoreMatrix,
) -> Stats {
    let mut stats = Stats::default();
    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();
    let neg_inf = i32::MIN / 4;

    for target in targets {
        let (gross, net) = target.cells();
        stats.gross_cells += gross;
        stats.net_cells += net;

        if target.blank() || target.band() <= 0 {
            continue;
        }

        let qlen = target.query_length.max(0).min(
            query
                .len()
                .saturating_sub(target.query_start.max(0) as usize) as i32,
        ) as usize;
        let slen = target.seq.len();
        let rows = qlen + 1;
        let idx = |i: usize, j: usize| -> usize { j * rows + i };
        let mut h = vec![neg_inf; rows * (slen + 1)];
        let mut hgap = vec![neg_inf; rows * (slen + 1)];
        let mut vgap = vec![neg_inf; rows * (slen + 1)];

        if target.d_begin <= 0 && target.d_end > 0 {
            h[idx(0, 0)] = 0;
        }

        let mut best_score = -1;
        let mut best_i = 0usize;
        let mut best_j = 0usize;

        for j in 1..=slen {
            let target_pos = if target.reverse { slen - j } else { j - 1 };
            let lower = (target.d_begin + (j as i32 - 1)).max(0);
            let upper = (target.d_end + (j as i32 - 1)).min(qlen as i32 - 1);
            if lower > upper {
                continue;
            }
            for i in (lower as usize + 1)..=(upper as usize + 1) {
                let qpos = target.query_start + i as i32 - 1;
                let qidx = if target.reverse {
                    query.len() - 1 - qpos as usize
                } else {
                    qpos as usize
                };
                let q = query[qidx] & LETTER_MASK;
                let s = target.seq[target_pos] & LETTER_MASK;
                let current = idx(i, j);
                let diag = h[idx(i - 1, j - 1)] + score_matrix.score(q, s);
                let h_open = h[idx(i, j - 1)] - gap_open;
                let h_extend = hgap[idx(i, j - 1)] - gap_extend;
                hgap[current] = h_open.max(h_extend);
                let v_open = h[idx(i - 1, j)] - gap_open;
                let v_extend = vgap[idx(i - 1, j)] - gap_extend;
                vgap[current] = v_open.max(v_extend);
                let score = diag.max(hgap[current]).max(vgap[current]);
                h[current] = score;
                if score > best_score {
                    best_score = score;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        if best_score >= 0 {
            target.score = best_score + 1;
            target.query_end = best_i as i32;
            target.target_end = best_j as i32;
        }
    }

    stats
}

pub fn anchored_swipe_score(
    query: &[Letter],
    targets: &mut [Target],
    score_matrix: &ScoreMatrix,
) -> Stats {
    targets.sort_by_key(|target| target.band());
    smith_waterman(query, targets, score_matrix)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WrapperConfig {
    pub query_len: i32,
    pub sensitivity: Sensitivity,
    pub score_hint: i32,
}

#[derive(Clone, Copy)]
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
}

impl<'a> AnchoredSwipeConfig<'a> {
    pub fn new(query: &'a [Letter], score_matrix: &'a ScoreMatrix) -> Self {
        AnchoredSwipeConfig {
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
        }
    }
}

#[derive(Debug, Clone)]
pub struct Profiles {
    pub int16: LongScoreProfile<i16>,
}

impl Profiles {
    pub fn new(seq: &[Letter], cbs: Option<&[i8]>, padding: usize, matrix: &ScoreMatrix) -> Self {
        Profiles {
            int16: make_profile16(seq, cbs, padding, matrix),
        }
    }

    pub fn reverse(&self) -> Self {
        Profiles {
            int16: self.int16.reverse(),
        }
    }
}

pub fn get_band(_qlen: i32, sensitivity: Sensitivity) -> i32 {
    if sensitivity >= Sensitivity::UltraSensitive {
        160
    } else if sensitivity >= Sensitivity::MoreSensitive {
        96
    } else {
        32
    }
}

pub fn align_right(
    target_seq: &[Letter],
    reverse: bool,
    i: i32,
    j: i32,
    mut d_begin: i32,
    mut d_end: i32,
    _prefix_score: i32,
    targets: &mut TargetVector,
    target_idx: i64,
    cfg: &WrapperConfig,
) {
    let qlen = cfg.query_len - i;
    let mut tlen = target_seq.len() as i32;
    let band = get_band(cfg.query_len, cfg.sensitivity);
    d_begin -= band;
    d_end += band - 1;
    let d0 = geo::clip_diag(geo::diag_sub_matrix(d_begin, i, j), qlen, tlen);
    let d1 = geo::clip_diag(geo::diag_sub_matrix(d_end, i, j), qlen, tlen);
    tlen = tlen.min(geo::j(qlen - 1, d0) + 1);
    assert!(tlen > 0);
    assert!(d1 >= d0);

    let clipped_target = if reverse {
        target_seq[target_seq.len() - tlen as usize..].to_vec()
    } else {
        target_seq[..tlen as usize].to_vec()
    };
    targets.int16.push(Target::new(
        clipped_target,
        d0,
        d1 + 1,
        i,
        qlen,
        target_idx,
        reverse,
    ));
}

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
    let qlen = cfg.query_len;
    let tlen = target_seq.len() as i32;
    let ir = qlen - 1 - i;
    let jr = tlen - 1 - j;
    align_right(
        &target_seq[..=j as usize],
        true,
        ir,
        jr,
        geo::rev_diag(d_end, qlen, tlen),
        geo::rev_diag(d_begin, qlen, tlen),
        suffix_score,
        targets,
        target_idx,
        cfg,
    );
}

pub fn add_target(
    t: &DpTarget,
    targets: &mut TargetVector,
    target_idx: &mut i64,
    cfg: &WrapperConfig,
) {
    if t.extend_right(cfg.query_len) {
        let i = t.anchor.query_end;
        let j = t.anchor.subject_end;
        align_right(
            &t.seq[j as usize..],
            false,
            i,
            j,
            t.anchor.d_min_right,
            t.anchor.d_max_right,
            t.anchor.prefix_score,
            targets,
            *target_idx,
            cfg,
        );
        *target_idx += 1;
    }
    if t.extend_left() {
        let suffix_score = cfg.score_hint - t.anchor.prefix_score + t.anchor.score;
        align_left(
            &t.seq,
            t.anchor.query_begin - 1,
            t.anchor.subject_begin - 1,
            t.anchor.d_min_left,
            t.anchor.d_max_left,
            suffix_score,
            targets,
            *target_idx,
            cfg,
        );
        *target_idx += 1;
    }
}

pub fn anchored_swipe(targets: &mut Targets, cfg: &AnchoredSwipeConfig<'_>) -> Vec<Hsp> {
    let mut target_count = 0i64;
    let mut max_target_len = 0usize;
    if !cfg.target_profiles {
        for bin in targets.iter() {
            for t in bin.as_slice() {
                target_count += 1;
                max_target_len = max_target_len.max(t.seq.len());
            }
        }
    } else {
        for bin in targets.iter() {
            target_count += bin.size();
        }
    }

    let _profiles = if cfg.target_profiles {
        None
    } else {
        Some(Profiles::new(
            cfg.query,
            cfg.query_cbs,
            cfg.query.len() + max_target_len + 32,
            cfg.score_matrix,
        ))
    };
    let _profiles_rev = _profiles.as_ref().map(Profiles::reverse);

    let wrapper_cfg = WrapperConfig {
        query_len: cfg.query.len() as i32,
        sensitivity: cfg.sensitivity,
        score_hint: cfg.score_hint,
    };
    let mut target_vec = TargetVector {
        int16: Vec::with_capacity((target_count * 2).max(0) as usize),
    };
    let mut target_idx = 0i64;
    for bin in targets.iter() {
        for t in bin.as_slice() {
            add_target(t, &mut target_vec, &mut target_idx, &wrapper_cfg);
        }
    }

    target_vec.int16.sort_by_key(|target| target.band());
    smith_waterman(cfg.query, &mut target_vec.int16, cfg.score_matrix);
    target_vec
        .int16
        .sort_by(|a, b| a.target_idx.cmp(&b.target_idx));

    let mut target_it = target_vec.int16.iter();
    let mut out = Vec::new();
    let mut recompute = swipe::targets();

    for bin in 0..swipe::BINS {
        for t in targets[bin].as_slice() {
            if t.anchor.score == 0 {
                continue;
            }
            let mut score = t.anchor.score;
            let mut i0 = t.anchor.query_begin;
            let mut i1 = t.anchor.query_end;
            let mut j0 = t.anchor.subject_begin;
            let mut j1 = t.anchor.subject_end;
            if t.extend_right(cfg.query.len() as i32) {
                let right = target_it
                    .next()
                    .expect("missing right anchored SWIPE extension");
                score += right.score;
                i1 += right.query_end;
                j1 += right.target_end;
            }
            if t.extend_left() {
                let left = target_it
                    .next()
                    .expect("missing left anchored SWIPE extension");
                score += left.score;
                i0 -= left.query_end;
                j0 -= left.target_end;
            }

            let qcov = (i1 - i0) as f64 / cfg.query.len() as f64 * 100.0;
            let tcov = (j1 - j0) as f64 / t.seq.len() as f64 * 100.0;
            let cov_filtered = (cfg.query_or_target_cover > 0.0
                || cfg.query_cover > 0.0
                || cfg.subject_cover > 0.0)
                && (qcov.max(tcov) < cfg.query_or_target_cover
                    || qcov < cfg.query_cover
                    || tcov < cfg.subject_cover);
            if !cfg.recompute_adjusted && cov_filtered {
                continue;
            }
            let evalue = cfg
                .score_matrix
                .evalue(score, cfg.query.len() as u32, t.seq.len() as u32);
            if !cfg.recompute_adjusted && evalue > cfg.max_evalue {
                continue;
            }
            if cfg.recompute_adjusted && (qcov < 70.0 || tcov < 70.0) {
                continue;
            }
            let approx_id = stats::approx_id(score, i1 - i0, j1 - j0);
            if cfg.recompute_adjusted && approx_id < 70.0 {
                recompute[bin].push_back(DpTarget::new(
                    t.seq.clone(),
                    t.true_target_len,
                    t.d_begin,
                    t.d_end,
                    t.target_idx,
                    cfg.query.len() as i32,
                    Default::default(),
                    t.anchor,
                ));
            } else if !cov_filtered {
                let mut hsp = Hsp::new();
                hsp.score = score;
                hsp.evalue = evalue;
                hsp.bit_score = cfg.score_matrix.bitscore(score as f64);
                hsp.swipe_target = t.target_idx as i32;
                hsp.query_range = Interval::new(i0, i1);
                hsp.query_source_range = hsp.query_range;
                hsp.subject_range = Interval::new(j0, j1);
                hsp.subject_source_range = hsp.subject_range;
                hsp.target_seq = t.seq.clone();
                hsp.approx_id = hsp.approx_id_percent(cfg.query, &t.seq);
                out.push(hsp);
            }
        }
    }
    if cfg.recompute_adjusted {
        let mut params = Params::new(cfg.query, cfg.score_matrix);
        params.v = swipe::HspValues::COORDS;
        out.extend(swipe::swipe(&recompute, &mut params));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dp::swipe::Anchor;

    #[test]
    fn test_target_cells() {
        let target = Target::new(vec![0; 5], -1, 2, 0, 4, 7, false);
        assert_eq!(target.band(), 3);
        assert_eq!(target.gross_cells(), 15);
        assert_eq!(target.cells(), (15, 11));
    }

    #[test]
    fn test_anchored_swipe_self_band() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2, 3, 4, 5];
        let mut targets = vec![Target::new(
            query.clone(),
            -1,
            2,
            0,
            query.len() as i32,
            0,
            false,
        )];
        let stats = anchored_swipe_score(&query, &mut targets, &sm);
        assert!(stats.net_cells > 0);
        assert!(targets[0].score > 0);
        assert_eq!(targets[0].query_end, query.len() as i32);
        assert_eq!(targets[0].target_end, query.len() as i32);
    }

    #[test]
    fn test_align_right_adds_clipped_target() {
        let cfg = WrapperConfig {
            query_len: 12,
            sensitivity: Sensitivity::Default,
            score_hint: 0,
        };
        let mut targets = TargetVector::default();
        let seq: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6, 7];
        align_right(&seq, false, 3, 2, -1, 2, 0, &mut targets, 9, &cfg);
        assert_eq!(targets.int16.len(), 1);
        assert_eq!(targets.int16[0].target_idx, 9);
        assert!(!targets.int16[0].seq.is_empty());
        assert!(targets.int16[0].d_end > targets.int16[0].d_begin);
    }

    #[test]
    fn test_add_target_extends_both_sides() {
        let cfg = WrapperConfig {
            query_len: 20,
            sensitivity: Sensitivity::Default,
            score_hint: 30,
        };
        let target = DpTarget {
            seq: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
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
            matrix: None,
            matrix_scale: 1,
        };
        let mut out = TargetVector::default();
        let mut target_idx = 0;
        add_target(&target, &mut out, &mut target_idx, &cfg);
        assert_eq!(target_idx, 2);
        assert_eq!(out.int16.len(), 2);
    }

    #[test]
    fn test_profiles_reverse() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq: Vec<Letter> = vec![0, 1, 2, 3];
        let profiles = Profiles::new(&seq, None, seq.len() + 32, &sm);
        let reversed = profiles.reverse();
        assert_eq!(profiles.int16.length(), seq.len());
        assert_eq!(reversed.int16.length(), seq.len());
    }

    #[test]
    fn test_anchored_swipe_targets_outputs_hsp() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut targets = swipe::targets();
        let mut target = DpTarget::new(
            query.clone(),
            query.len() as i32,
            -1,
            2,
            11,
            query.len() as i32,
            Default::default(),
            Default::default(),
        );
        target.anchor = Anchor {
            query_begin: 2,
            query_end: 5,
            subject_begin: 2,
            subject_end: 5,
            d_min_left: -1,
            d_max_left: 2,
            d_min_right: -1,
            d_max_right: 2,
            prefix_score: 10,
            score: 30,
        };
        targets[0].push_back(target);
        let cfg = AnchoredSwipeConfig::new(&query, &sm);
        let out = anchored_swipe(&mut targets, &cfg);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].swipe_target, 11);
        assert!(out[0].query_range.begin <= 2);
        assert!(out[0].subject_range.begin <= 2);
        assert!(out[0].query_range.end >= 5);
        assert!(out[0].subject_range.end >= 5);
        assert!(out[0].score > 30);
    }

    #[test]
    fn test_anchored_swipe_exact_identity_approx_id_is_percent_identity() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2];
        let mut targets = swipe::targets();
        let mut target = DpTarget::new(
            query.clone(),
            query.len() as i32,
            -1,
            2,
            17,
            query.len() as i32,
            Default::default(),
            Default::default(),
        );
        target.anchor = Anchor {
            query_begin: 0,
            query_end: 3,
            subject_begin: 0,
            subject_end: 3,
            d_min_left: -1,
            d_max_left: 2,
            d_min_right: -1,
            d_max_right: 2,
            prefix_score: 0,
            score: 1,
        };
        targets[0].push_back(target);
        let cfg = AnchoredSwipeConfig::new(&query, &sm);
        let out = anchored_swipe(&mut targets, &cfg);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].approx_id, 100.0);
        assert!(stats::approx_id(out[0].score, 3, 3) < 100.0);
    }
}
