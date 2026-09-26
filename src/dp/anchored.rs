//! Scalar anchored SWIPE scoring.
//!
//! This ports the score-only recurrence from `dp/swipe/anchored.h` into a
//! straightforward scalar implementation. The C++ version batches targets into
//! SIMD lanes and consumes subject sequence in blocks; this version keeps the
//! same anchored band semantics per target.

use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;
use std::sync::Arc;

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
    /// Target-specific composition-adjusted scores (`DpTarget::prof` in C++).
    pub matrix: Option<Arc<TargetMatrix>>,
    pub matrix_scale: i32,
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
            matrix: None,
            matrix_scale: 1,
        }
    }

    pub fn with_matrix(mut self, matrix: Arc<TargetMatrix>, matrix_scale: i32) -> Self {
        self.matrix = Some(matrix);
        self.matrix_scale = matrix_scale.max(1);
        self
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
    let neg_inf = i32::MIN / 4;

    for target in targets {
        let (gross, net) = target.cells();
        stats.gross_cells += gross;
        stats.net_cells += net;

        if target.blank() || target.band() <= 0 {
            continue;
        }

        let scale = target.matrix_scale.max(1);
        let gap_open = (score_matrix.gap_open() + score_matrix.gap_extend()) * scale;
        let gap_extend = score_matrix.gap_extend() * scale;

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
                let substitution = target.matrix.as_ref().map_or_else(
                    || score_matrix.score(q, s),
                    |matrix| matrix.scores[s as usize * 32 + q as usize] as i32,
                );
                let diag = h[idx(i - 1, j - 1)] + substitution;
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

/* The C++ wrapper lives in swipe/anchored_wrapper.cpp. */
pub use crate::dp::swipe::anchored_wrapper::{
    add_target, align_left, align_right, anchored_swipe, get_band, select_matrix, swipe_threads,
    AnchoredSwipeConfig, Profiles, TargetVector, WrapperConfig,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn target_cells_match_banded_geometry() {
        let target = Target::new(vec![0; 5], -1, 2, 0, 4, 7, false);
        assert_eq!(target.band(), 3);
        assert_eq!(target.gross_cells(), 15);
        assert_eq!(target.cells(), (15, 11));
    }

    #[test]
    fn anchored_swipe_score_matches_self_band() {
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
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
        let stats = anchored_swipe_score(&query, &mut targets, &matrix);
        assert!(stats.net_cells > 0);
        assert!(targets[0].score > 0);
        assert_eq!(targets[0].query_end, query.len() as i32);
        assert_eq!(targets[0].target_end, query.len() as i32);
    }
}
