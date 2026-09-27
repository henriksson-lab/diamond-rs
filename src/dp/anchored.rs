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

/// Score anchored extensions in eight AVX2 lanes using the moving diagonal
/// band from C++ `dp/swipe/anchored.h`.
///
/// `i32` lanes intentionally replace DIAMOND's `i16` lanes here.  This keeps
/// the same recurrence while avoiding a second overflow/recompute path.  The
/// scalar implementation remains the runtime fallback on non-AVX2 hosts.
pub fn smith_waterman_simd(
    query: &[Letter],
    targets: &mut [Target],
    score_matrix: &ScoreMatrix,
) -> Stats {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: AVX2 was detected at runtime. The implementation only
            // accesses slices after explicit bounds checks.
            return unsafe { smith_waterman_avx2(query, targets, score_matrix) };
        }
    }
    smith_waterman(query, targets, score_matrix)
}

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn smith_waterman_avx2(
    query: &[Letter],
    targets: &mut [Target],
    score_matrix: &ScoreMatrix,
) -> Stats {
    const LANES: usize = 8;
    // Low enough that subtracting ordinary protein gap penalties stays well
    // away from i32 underflow, but lower than every reachable score.
    const NEG: i32 = i32::MIN / 4;

    let mut stats = Stats::default();
    for target in targets.iter() {
        let (gross, net) = target.cells();
        stats.gross_cells += gross;
        stats.net_cells += net;
    }

    for chunk in targets.chunks_mut(LANES) {
        let band = chunk
            .iter()
            .filter(|target| !target.blank())
            .map(Target::band)
            .max()
            .unwrap_or(0)
            .max(0) as usize;
        if band == 0 {
            continue;
        }
        let max_subject_len = chunk
            .iter()
            .map(|target| target.seq.len())
            .max()
            .unwrap_or(0);
        let neg = arch::_mm256_set1_epi32(NEG);
        let mut prev_h = vec![neg; band];
        let mut curr_h = vec![neg; band];
        let mut prev_gap = vec![neg; band + 1];
        let mut curr_gap = vec![neg; band + 1];

        // C++ Matrix::init_channel_nw places the sole reachable origin at
        // -d_begin in the moving band and leaves all gap states unreachable.
        for (lane, target) in chunk.iter().enumerate() {
            if target.blank() || target.band() <= 0 {
                continue;
            }
            let origin = -target.d_begin;
            if origin >= 0 && (origin as usize) < band {
                let mut values = [NEG; LANES];
                arch::_mm256_storeu_si256(values.as_mut_ptr().cast(), prev_h[origin as usize]);
                values[lane] = 0;
                prev_h[origin as usize] = arch::_mm256_loadu_si256(values.as_ptr().cast());
            }
        }

        let mut best = [-1i32; LANES];
        let mut best_i = [0i32; LANES];
        let mut best_j = [0i32; LANES];

        for j in 0..max_subject_len {
            curr_h.fill(neg);
            curr_gap.fill(neg);
            let mut vertical = neg;
            for r in 0..band {
                let mut subst = [0i32; LANES];
                let mut valid = [0i32; LANES];
                let mut gap_open = [0i32; LANES];
                let mut gap_extend = [0i32; LANES];
                let mut q_positions = [0i32; LANES];
                for (lane, target) in chunk.iter().enumerate() {
                    let qrel = target.d_begin + j as i32 + r as i32;
                    let qlen = target.query_length.max(0).min(
                        query
                            .len()
                            .saturating_sub(target.query_start.max(0) as usize)
                            as i32,
                    );
                    let in_band = !target.blank()
                        && j < target.seq.len()
                        && (r as i32) < target.band()
                        && qrel >= 0
                        && qrel < qlen;
                    if !in_band {
                        continue;
                    }
                    valid[lane] = -1;
                    q_positions[lane] = qrel;
                    let qpos = target.query_start + qrel;
                    let qidx = if target.reverse {
                        query.len() - 1 - qpos as usize
                    } else {
                        qpos as usize
                    };
                    let target_pos = if target.reverse {
                        target.seq.len() - 1 - j
                    } else {
                        j
                    };
                    let q = query[qidx] & LETTER_MASK;
                    let s = target.seq[target_pos] & LETTER_MASK;
                    subst[lane] = target.matrix.as_ref().map_or_else(
                        || score_matrix.score(q, s),
                        |matrix| matrix.scores[s as usize * 32 + q as usize] as i32,
                    );
                    let scale = target.matrix_scale.max(1);
                    gap_open[lane] = (score_matrix.gap_open() + score_matrix.gap_extend()) * scale;
                    gap_extend[lane] = score_matrix.gap_extend() * scale;
                }

                let mask = arch::_mm256_loadu_si256(valid.as_ptr().cast());
                let substitution = arch::_mm256_loadu_si256(subst.as_ptr().cast());
                let go = arch::_mm256_loadu_si256(gap_open.as_ptr().cast());
                let ge = arch::_mm256_loadu_si256(gap_extend.as_ptr().cast());
                let diag = arch::_mm256_add_epi32(prev_h[r], substitution);
                // The horizontal gap state moves one row up when the diagonal
                // band advances to the next subject column. Vertical is the
                // running state within this column.
                let horizontal = prev_gap[r + 1];
                let score =
                    arch::_mm256_max_epi32(diag, arch::_mm256_max_epi32(horizontal, vertical));
                let score = arch::_mm256_blendv_epi8(neg, score, mask);
                curr_h[r] = score;
                let open = arch::_mm256_sub_epi32(score, go);
                let next_horizontal =
                    arch::_mm256_max_epi32(arch::_mm256_sub_epi32(horizontal, ge), open);
                vertical = arch::_mm256_max_epi32(arch::_mm256_sub_epi32(vertical, ge), open);
                curr_gap[r] = arch::_mm256_blendv_epi8(neg, next_horizontal, mask);
                vertical = arch::_mm256_blendv_epi8(neg, vertical, mask);

                let mut values = [NEG; LANES];
                arch::_mm256_storeu_si256(values.as_mut_ptr().cast(), score);
                for lane in 0..chunk.len() {
                    if valid[lane] != 0 && values[lane] > best[lane] {
                        best[lane] = values[lane];
                        best_i[lane] = q_positions[lane] + 1;
                        best_j[lane] = j as i32 + 1;
                    }
                }
            }
            std::mem::swap(&mut prev_h, &mut curr_h);
            std::mem::swap(&mut prev_gap, &mut curr_gap);
        }

        for lane in 0..chunk.len() {
            if best[lane] >= 0 {
                chunk[lane].score = best[lane] + 1;
                chunk[lane].query_end = best_i[lane];
                chunk[lane].target_end = best_j[lane];
            }
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
    smith_waterman_simd(query, targets, score_matrix)
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

    #[test]
    fn anchored_simd_matches_scalar_across_batches_and_bands() {
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let mut state = 0x9e37_79b9u32;
        let mut next_letter = || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            ((state >> 16) % 20) as Letter
        };
        for count in 1..=19 {
            let query: Vec<Letter> = (0..73).map(|_| next_letter()).collect();
            let mut expected = Vec::new();
            for lane in 0..count {
                let len = 9 + lane * 3 % 31;
                let seq = (0..len).map(|_| next_letter()).collect();
                let query_start = (lane * 2 % 13) as i32;
                let query_length = (query.len() as i32 - query_start).min(17 + lane as i32);
                let d_begin = -((lane * 3 % 7) as i32);
                let width = 2 + (lane * 5 % 15) as i32;
                expected.push(Target::new(
                    seq,
                    d_begin,
                    d_begin + width,
                    query_start,
                    query_length,
                    lane as i64,
                    lane % 3 == 0,
                ));
            }
            let mut actual = expected.clone();
            let expected_stats = smith_waterman(&query, &mut expected, &matrix);
            let actual_stats = smith_waterman_simd(&query, &mut actual, &matrix);
            assert_eq!(actual_stats, expected_stats, "count={count}");
            for lane in 0..count {
                assert_eq!(
                    actual[lane].score, expected[lane].score,
                    "count={count} lane={lane}"
                );
                assert_eq!(
                    (actual[lane].query_end, actual[lane].target_end),
                    (expected[lane].query_end, expected[lane].target_end),
                    "count={count} lane={lane}"
                );
            }
        }
    }
}
