//! AVX2 score-only Smith-Waterman over up to sixteen target lanes.
//!
//! The lanes share a query, but may have different lengths and diagonal
//! bands.  Invalid band cells are masked before their horizontal/vertical
//! gap state can escape into a valid cell.  Scores which reach the `i16`
//! ceiling are reported in `overflow_mask` and must be recomputed with the
//! scalar `i32` kernel.

use crate::basic::value::Letter;
use crate::stats::score_matrix::ScoreMatrix;

/// One target lane for [`score_batch_avx2`].
#[derive(Clone, Copy, Debug)]
pub struct ScoreTarget<'a> {
    pub subject: &'a [Letter],
    pub d_begin: i32,
    pub d_end: i32,
}

/// Result from one SIMD batch. Only the first `len` entries are meaningful.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BatchScores {
    pub scores: [i32; 16],
    pub overflow_mask: u16,
    pub len: usize,
}

/// Reusable aligned DP rows for the AVX2 kernel.
///
/// Keeping this at the caller avoids allocating three query-sized rows for
/// every group of sixteen targets.
#[derive(Default)]
pub struct SimdScoreScratch {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_h: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_h: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_e: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_e: Vec<ArchVector>,
}

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type ArchVector = arch::__m256i;

/// Score at most sixteen targets with AVX2.
///
/// Returns `None` when AVX2 is unavailable, the batch is empty or too large,
/// or the penalties cannot be represented safely by this `i16` kernel.
/// Lanes selected by `overflow_mask` have to be recomputed in `i32`.
/// Adjusted target matrices use different score scaling and must stay on the
/// scalar path; this entry point always uses the standard matrix scale.
pub fn score_batch_avx2(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut SimdScoreScratch,
) -> Option<BatchScores> {
    if targets.is_empty()
        || targets.len() > 16
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
    {
        return None;
    }
    let gap_open = score_matrix
        .gap_open()
        .checked_add(score_matrix.gap_extend())?;
    let gap_extend = score_matrix.gap_extend();
    // NEG below is -16384. Keeping penalties below this limit means all
    // deliberately unreachable-state subtraction remains representable.
    if !(0..=16_000).contains(&gap_open) || !(0..=16_000).contains(&gap_extend) {
        return None;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        // SAFETY: guarded by runtime AVX2 detection; all memory accesses in
        // the implementation are bounds checked or use fixed-size arrays.
        return Some(unsafe {
            score_batch_avx2_impl(
                query,
                targets,
                score_matrix.matrix16(),
                query_cbs,
                gap_open as i16,
                gap_extend as i16,
                scratch,
            )
        });
    }

    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, score_matrix, query_cbs, scratch);
        None
    }
}

/// Score up to sixteen complete Smith-Waterman matrices in AVX2 lanes.
///
/// This is the SIMD counterpart of C++ `dp/swipe/full_swipe.h` for the
/// score-only configuration. Unlike [`score_batch_avx2`], it walks query rows
/// directly and therefore does not spend work on the triangular padding of a
/// synthetic full-width diagonal band.
pub fn score_full_batch_avx2(
    query: &[Letter],
    targets: &[&[Letter]],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut SimdScoreScratch,
) -> Option<BatchScores> {
    if targets.is_empty()
        || targets.len() > 16
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
    {
        return None;
    }
    let gap_open = score_matrix
        .gap_open()
        .checked_add(score_matrix.gap_extend())?;
    let gap_extend = score_matrix.gap_extend();
    if !(0..=16_000).contains(&gap_open) || !(0..=16_000).contains(&gap_extend) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        // SAFETY: guarded by runtime AVX2 detection.
        return Some(unsafe {
            score_full_batch_avx2_impl(
                query,
                targets,
                score_matrix.matrix16(),
                query_cbs,
                gap_open as i16,
                gap_extend as i16,
                scratch,
            )
        });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, score_matrix, query_cbs, scratch);
        None
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_full_batch_avx2_impl(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &[i16; 32 * 32],
    query_cbs: &[i8],
    gap_open: i16,
    gap_extend: i16,
    scratch: &mut SimdScoreScratch,
) -> BatchScores {
    const NEG: i16 = -16_384;
    let rows = query.len() + 1;
    let zero = arch::_mm256_setzero_si256();
    let neg = arch::_mm256_set1_epi16(NEG);
    scratch.prev_h.resize(rows, zero);
    scratch.curr_h.resize(rows, zero);
    scratch.prev_e.resize(rows, neg);
    scratch.curr_e.resize(rows, neg);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(neg);
    let go = arch::_mm256_set1_epi16(gap_open);
    let ge = arch::_mm256_set1_epi16(gap_extend);
    let max_i16 = arch::_mm256_set1_epi16(i16::MAX);
    let mut best = zero;
    let mut overflow = zero;
    let max_subject_len = targets.iter().map(|target| target.len()).max().unwrap_or(0);

    for j in 0..max_subject_len {
        scratch.curr_h[0] = zero;
        scratch.curr_e[0] = neg;
        let mut vertical = neg;
        for (qpos, &ql) in query.iter().enumerate() {
            let mut subst = [0i16; 16];
            let mut valid = [0i16; 16];
            for lane in 0..targets.len() {
                if j >= targets[lane].len() {
                    continue;
                }
                valid[lane] = -1;
                let sl = targets[lane][j];
                let base = if sl & crate::basic::value::SEED_MASK != 0 {
                    0
                } else {
                    matrix[((ql & crate::basic::value::LETTER_MASK) as usize) * 32
                        + (sl & crate::basic::value::LETTER_MASK) as usize]
                        as i32
                };
                let value = base + query_cbs.get(qpos).copied().unwrap_or(0) as i32;
                if !(i16::MIN as i32..=i16::MAX as i32).contains(&value) {
                    let mut lane_mask = [0i16; 16];
                    lane_mask[lane] = -1;
                    overflow = arch::_mm256_or_si256(
                        overflow,
                        arch::_mm256_loadu_si256(lane_mask.as_ptr().cast()),
                    );
                }
                subst[lane] = value.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
            }
            let mask = arch::_mm256_loadu_si256(valid.as_ptr().cast());
            let substitution = arch::_mm256_loadu_si256(subst.as_ptr().cast());
            let diag = arch::_mm256_adds_epi16(scratch.prev_h[qpos], substitution);
            let horizontal = scratch.prev_e[qpos + 1];
            let mut score = arch::_mm256_max_epi16(diag, horizontal);
            score = arch::_mm256_max_epi16(score, vertical);
            score = arch::_mm256_max_epi16(score, zero);
            score = arch::_mm256_and_si256(score, mask);
            overflow = arch::_mm256_or_si256(overflow, arch::_mm256_cmpeq_epi16(score, max_i16));
            let open = arch::_mm256_subs_epi16(score, go);
            let next_horizontal =
                arch::_mm256_max_epi16(arch::_mm256_subs_epi16(horizontal, ge), open);
            vertical = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(vertical, ge), open);
            scratch.curr_h[qpos + 1] = score;
            scratch.curr_e[qpos + 1] = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, next_horizontal),
                arch::_mm256_andnot_si256(mask, neg),
            );
            vertical = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, vertical),
                arch::_mm256_andnot_si256(mask, neg),
            );
            best = arch::_mm256_max_epi16(best, score);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }

    let mut raw_scores = [0i16; 16];
    let mut raw_overflow = [0i16; 16];
    arch::_mm256_storeu_si256(raw_scores.as_mut_ptr().cast(), best);
    arch::_mm256_storeu_si256(raw_overflow.as_mut_ptr().cast(), overflow);
    let mut scores = [0i32; 16];
    let mut overflow_mask = 0u16;
    for lane in 0..targets.len() {
        scores[lane] = raw_scores[lane] as i32;
        if raw_overflow[lane] != 0 {
            overflow_mask |= 1 << lane;
        }
    }
    BatchScores {
        scores,
        overflow_mask,
        len: targets.len(),
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_batch_avx2_impl(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &[i16; 32 * 32],
    query_cbs: &[i8],
    gap_open: i16,
    gap_extend: i16,
    scratch: &mut SimdScoreScratch,
) -> BatchScores {
    const NEG: i16 = -16_384;
    // Index rows by diagonal position rather than absolute query position:
    // qpos = d_begin + subject_pos + band_row.  The diagonal predecessor is
    // then the same row in the previous column, and the vertical predecessor
    // is row + 1.  This is the moving-band layout used by DIAMOND and limits
    // work to max(band width), not query length.
    let band_rows = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    if band_rows == 0 {
        return BatchScores {
            scores: [0; 16],
            overflow_mask: 0,
            len: targets.len(),
        };
    }
    let rows = band_rows + 1;
    let zero = arch::_mm256_setzero_si256();
    let neg = arch::_mm256_set1_epi16(NEG);
    scratch.prev_h.resize(rows, zero);
    scratch.curr_h.resize(rows, zero);
    scratch.prev_e.resize(rows, neg);
    scratch.curr_e.resize(rows, neg);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(neg);

    let gap_open_v = arch::_mm256_set1_epi16(gap_open);
    let gap_extend_v = arch::_mm256_set1_epi16(gap_extend);
    let max_i16 = arch::_mm256_set1_epi16(i16::MAX);
    let mut best = zero;
    let mut overflow = zero;
    let max_subject_len = targets.iter().map(|t| t.subject.len()).max().unwrap_or(0);

    for target_pos in 0..max_subject_len {
        // The sentinel row supplies E(row + 1) at the upper band edge.
        scratch.curr_h[band_rows] = zero;
        scratch.curr_e[band_rows] = neg;
        let mut f = neg;

        for band_row in 0..band_rows {
            let mut substitution = [0i16; 16];
            let mut valid = [0i16; 16];
            for lane in 0..targets.len() {
                let target = targets[lane];
                let width = (target.d_end - target.d_begin).max(0) as usize;
                let qpos = target.d_begin + target_pos as i32 + band_row as i32;
                let in_band = target_pos < target.subject.len()
                    && band_row < width
                    && qpos >= 0
                    && qpos < query.len() as i32;
                if in_band {
                    valid[lane] = -1;
                    let qpos = qpos as usize;
                    let sl = target.subject[target_pos];
                    let base = if sl & crate::basic::value::SEED_MASK != 0 {
                        0
                    } else {
                        let qi = (query[qpos] & crate::basic::value::LETTER_MASK) as usize;
                        let si = (sl & crate::basic::value::LETTER_MASK) as usize;
                        matrix[qi * 32 + si] as i32
                    };
                    let cbs = if query_cbs.is_empty() {
                        0
                    } else {
                        query_cbs[qpos] as i32
                    };
                    let value = base + cbs;
                    if value > i16::MAX as i32 || value < i16::MIN as i32 {
                        overflow = arch::_mm256_or_si256(
                            overflow,
                            arch::_mm256_setr_epi16(
                                if lane == 0 { -1 } else { 0 },
                                if lane == 1 { -1 } else { 0 },
                                if lane == 2 { -1 } else { 0 },
                                if lane == 3 { -1 } else { 0 },
                                if lane == 4 { -1 } else { 0 },
                                if lane == 5 { -1 } else { 0 },
                                if lane == 6 { -1 } else { 0 },
                                if lane == 7 { -1 } else { 0 },
                                if lane == 8 { -1 } else { 0 },
                                if lane == 9 { -1 } else { 0 },
                                if lane == 10 { -1 } else { 0 },
                                if lane == 11 { -1 } else { 0 },
                                if lane == 12 { -1 } else { 0 },
                                if lane == 13 { -1 } else { 0 },
                                if lane == 14 { -1 } else { 0 },
                                if lane == 15 { -1 } else { 0 },
                            ),
                        );
                    }
                    substitution[lane] = value.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
                }
            }

            let subst = arch::_mm256_loadu_si256(substitution.as_ptr().cast());
            let mask = arch::_mm256_loadu_si256(valid.as_ptr().cast());
            let diag = arch::_mm256_adds_epi16(scratch.prev_h[band_row], subst);
            let mut h = arch::_mm256_max_epi16(diag, scratch.prev_e[band_row + 1]);
            h = arch::_mm256_max_epi16(h, f);
            h = arch::_mm256_max_epi16(h, zero);
            h = arch::_mm256_and_si256(h, mask);
            overflow = arch::_mm256_or_si256(overflow, arch::_mm256_cmpeq_epi16(h, max_i16));

            let open = arch::_mm256_subs_epi16(h, gap_open_v);
            let e = arch::_mm256_max_epi16(
                arch::_mm256_subs_epi16(scratch.prev_e[band_row + 1], gap_extend_v),
                open,
            );
            f = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(f, gap_extend_v), open);
            scratch.curr_h[band_row] = h;
            scratch.curr_e[band_row] = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, e),
                arch::_mm256_andnot_si256(mask, neg),
            );
            f = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, f),
                arch::_mm256_andnot_si256(mask, neg),
            );
            best = arch::_mm256_max_epi16(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }

    let mut score_lanes = [0i16; 16];
    let mut overflow_lanes = [0i16; 16];
    arch::_mm256_storeu_si256(score_lanes.as_mut_ptr().cast(), best);
    arch::_mm256_storeu_si256(overflow_lanes.as_mut_ptr().cast(), overflow);
    let mut scores = [0i32; 16];
    let mut overflow_mask = 0u16;
    for lane in 0..targets.len() {
        scores[lane] = score_lanes[lane] as i32;
        if overflow_lanes[lane] != 0 {
            overflow_mask |= 1 << lane;
        }
    }
    BatchScores {
        scores,
        overflow_mask,
        len: targets.len(),
    }
}

#[cfg(all(test, any(target_arch = "x86", target_arch = "x86_64")))]
mod tests {
    use super::*;
    use crate::basic::value::{LETTER_MASK, SEED_MASK};

    fn scalar(query: &[Letter], target: ScoreTarget<'_>, matrix: &ScoreMatrix, cbs: &[i8]) -> i32 {
        let qlen = query.len();
        let neg = i32::MIN / 4;
        let go = matrix.gap_open() + matrix.gap_extend();
        let ge = matrix.gap_extend();
        let mut ph = vec![0; qlen + 1];
        let mut pe = vec![neg; qlen + 1];
        let mut best = 0;
        for (spos, &s) in target.subject.iter().enumerate() {
            let mut ch = vec![0; qlen + 1];
            let mut ce = vec![neg; qlen + 1];
            let mut f = neg;
            for i in 1..=qlen {
                let qpos = i - 1;
                let valid = qpos as i32 >= target.d_begin + spos as i32
                    && (qpos as i32) < target.d_end + spos as i32;
                if !valid {
                    continue;
                }
                let subst = if s & SEED_MASK != 0 {
                    0
                } else {
                    matrix.score(query[qpos] & LETTER_MASK, s & LETTER_MASK)
                } + if cbs.is_empty() { 0 } else { cbs[qpos] as i32 };
                let h = (ph[i - 1] + subst).max(pe[i]).max(f).max(0);
                ch[i] = h;
                ce[i] = (pe[i] - ge).max(h - go);
                f = (f - ge).max(h - go);
                best = best.max(h);
            }
            ph = ch;
            pe = ce;
        }
        best
    }

    #[test]
    fn randomized_avx2_matches_i32_scalar_with_strict_bands() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0x5eed_cafe_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        let mut scratch = SimdScoreScratch::default();
        for _ in 0..80 {
            let qlen = 1 + next() as usize % 73;
            let query: Vec<Letter> = (0..qlen).map(|_| (next() % 25) as Letter).collect();
            let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
            let lane_count = 1 + next() as usize % 16;
            let mut subjects = Vec::with_capacity(lane_count);
            let mut bands = Vec::with_capacity(lane_count);
            for _ in 0..lane_count {
                let len = 1 + next() as usize % 81;
                let mut subject: Vec<Letter> = (0..len).map(|_| (next() % 25) as Letter).collect();
                if next() % 3 == 0 {
                    let pos = next() as usize % len;
                    subject[pos] |= SEED_MASK;
                }
                subjects.push(subject);
                let begin = next() as i32 % 31 - 15;
                let width = 1 + next() as i32 % 40;
                bands.push((begin, begin + width));
            }
            let targets: Vec<_> = subjects
                .iter()
                .zip(&bands)
                .map(|(subject, &(d_begin, d_end))| ScoreTarget {
                    subject,
                    d_begin,
                    d_end,
                })
                .collect();
            let got = score_batch_avx2(&query, &targets, &matrix, &cbs, &mut scratch).unwrap();
            assert_eq!(got.overflow_mask, 0);
            for lane in 0..lane_count {
                assert_eq!(
                    got.scores[lane],
                    scalar(&query, targets[lane], &matrix, &cbs),
                    "lane {lane}, qlen {qlen}, band {:?}",
                    bands[lane]
                );
            }
        }
    }

    #[test]
    fn reports_i16_score_overflow_per_lane() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        // W/W scores 11 in BLOSUM62, so this lane exceeds i16::MAX.
        let query = vec![17; 4_000];
        let long = vec![17; 4_000];
        let short = vec![17; 20];
        let targets = [
            ScoreTarget {
                subject: &long,
                d_begin: 0,
                d_end: 1,
            },
            ScoreTarget {
                subject: &short,
                d_begin: 0,
                d_end: 1,
            },
        ];
        let got = score_batch_avx2(
            &query,
            &targets,
            &matrix,
            &[],
            &mut SimdScoreScratch::default(),
        )
        .unwrap();
        assert_eq!(got.overflow_mask, 1);
        assert_eq!(got.scores[1], 220);
    }

    #[test]
    fn randomized_full_matrix_avx2_matches_scalar() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0xd1a0_600du64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        let mut scratch = SimdScoreScratch::default();
        for _ in 0..48 {
            let qlen = 1 + next() as usize % 65;
            let query: Vec<Letter> = (0..qlen).map(|_| (next() % 25) as Letter).collect();
            let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
            let count = 1 + next() as usize % 16;
            let subjects: Vec<Vec<Letter>> = (0..count)
                .map(|_| {
                    let len = 1 + next() as usize % 67;
                    (0..len).map(|_| (next() % 25) as Letter).collect()
                })
                .collect();
            let refs: Vec<&[Letter]> = subjects.iter().map(Vec::as_slice).collect();
            let got = score_full_batch_avx2(&query, &refs, &matrix, &cbs, &mut scratch).unwrap();
            assert_eq!(got.overflow_mask, 0);
            for lane in 0..count {
                let target = ScoreTarget {
                    subject: refs[lane],
                    d_begin: -(refs[lane].len() as i32 - 1),
                    d_end: qlen as i32,
                };
                assert_eq!(got.scores[lane], scalar(&query, target, &matrix, &cbs));
            }
        }
    }
}
