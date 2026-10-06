//! AVX2 banded Smith-Waterman with per-lane traceback.
//!
//! This mirrors the vector-cell/trace-mask configuration in DIAMOND's
//! `dp/swipe/{banded,full}_swipe.h`. Scores use eight `i32` lanes, avoiding
//! the overflow cascade needed by the C++ 8/16-bit dispatch.

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::basic::value::LETTER_MASK;
use crate::dp::smith_waterman::SwResult;
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;

#[derive(Clone, Copy)]
pub struct TraceTarget<'a> {
    pub subject: &'a [Letter],
    pub d_begin: i32,
    pub d_end: i32,
    pub matrix: Option<&'a TargetMatrix>,
    pub matrix_scale: i32,
}

pub struct TraceBatch {
    pub results: Vec<SwResult>,
    pub overflow_mask: u32,
}

/// Dispatch the traceback score width selected by SWIPE's score bin.
///
/// Byte and word lanes use the upstream signed-min score bias, and report
/// saturated lanes for promotion into the next bin. Target-adjusted and
/// ordinary lanes may share a batch; adjusted lanes use their private profile
/// and scaled gap vectors while ordinary lanes retain query CBS.
pub fn trace_batch_tier_avx2(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    score_bin: usize,
    semi_global: bool,
) -> Option<TraceBatch> {
    match score_bin {
        0 => super::simd_trace_narrow::trace_batch_i8(
            query,
            targets,
            score_matrix,
            query_cbs,
            semi_global,
        )
        .map(|batch| TraceBatch {
            results: batch.results,
            overflow_mask: batch.overflow_mask,
        }),
        1 => super::simd_trace_narrow::trace_batch_i16(
            query,
            targets,
            score_matrix,
            query_cbs,
            semi_global,
        )
        .map(|batch| TraceBatch {
            results: batch.results,
            overflow_mask: batch.overflow_mask,
        }),
        _ => trace_batch_avx2(query, targets, score_matrix, query_cbs).map(|results| TraceBatch {
            results,
            overflow_mask: 0,
        }),
    }
}

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

pub fn trace_batch_avx2(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> Option<Vec<SwResult>> {
    let invalid_band = targets.iter().any(|target| target.d_end <= target.d_begin);
    if targets.is_empty()
        || targets.len() > 8
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
        || invalid_band
    {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: guarded by runtime AVX2 detection.
            return Some(unsafe { trace_batch_avx2_impl(query, targets, score_matrix, query_cbs) });
        }
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    let _ = (query, targets, score_matrix, query_cbs);
    None
}

/// Score-only AVX2 path for lanes with composition-adjusted target matrices.
/// These lanes need per-target substitution tables and gap scales, so they
/// cannot use the shared-profile i8/i16 kernels.
pub fn score_adjusted_batch_avx2(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> Option<Vec<i32>> {
    let invalid_band = targets.iter().any(|target| target.d_end <= target.d_begin);
    if targets.is_empty()
        || targets.len() > 8
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
        || invalid_band
    {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        // SAFETY: guarded by runtime AVX2 detection.
        return Some(unsafe { score_adjusted_impl(query, targets, score_matrix, query_cbs) });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, score_matrix, query_cbs);
        None
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_adjusted_impl(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> Vec<i32> {
    const LANES: usize = 8;
    const NEG: i32 = i32::MIN / 4;
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin) as usize)
        .max()
        .unwrap_or(0);
    let max_len = targets
        .iter()
        .map(|target| target.subject.len())
        .max()
        .unwrap_or(0);
    let zero = arch::_mm256_setzero_si256();
    let neg = arch::_mm256_set1_epi32(NEG);
    let mut prev_h = vec![zero; band];
    let mut curr_h = vec![zero; band];
    let mut prev_e = vec![neg; band + 1];
    let mut curr_e = vec![neg; band + 1];
    let mut best = zero;
    for j in 0..max_len {
        curr_h.fill(zero);
        curr_e.fill(neg);
        let mut vertical = neg;
        for r in 0..band {
            let mut valid = [0i32; LANES];
            let mut subst = [0i32; LANES];
            let mut go = [0i32; LANES];
            let mut ge = [0i32; LANES];
            for lane in 0..targets.len() {
                let target = targets[lane];
                let q = target.d_begin + j as i32 + r as i32;
                if j >= target.subject.len()
                    || r >= (target.d_end - target.d_begin) as usize
                    || q < 0
                    || q >= query.len() as i32
                {
                    continue;
                }
                valid[lane] = -1;
                let ql = query[q as usize];
                let sl = target.subject[j];
                subst[lane] = if let Some(matrix) = target.matrix {
                    matrix.scores[(sl & LETTER_MASK) as usize * 32 + (ql & LETTER_MASK) as usize]
                        as i32
                } else {
                    score_matrix.score(ql, sl)
                        + query_cbs.get(q as usize).copied().unwrap_or(0) as i32
                };
                let scale = target.matrix_scale.max(1);
                go[lane] = (score_matrix.gap_open() + score_matrix.gap_extend()) * scale;
                ge[lane] = score_matrix.gap_extend() * scale;
            }
            let mask = arch::_mm256_loadu_si256(valid.as_ptr().cast());
            let diag =
                arch::_mm256_add_epi32(prev_h[r], arch::_mm256_loadu_si256(subst.as_ptr().cast()));
            let horizontal = prev_e[r + 1];
            let mut h = arch::_mm256_max_epi32(diag, horizontal);
            h = arch::_mm256_max_epi32(h, vertical);
            h = arch::_mm256_max_epi32(h, zero);
            h = arch::_mm256_and_si256(h, mask);
            let open = arch::_mm256_sub_epi32(h, arch::_mm256_loadu_si256(go.as_ptr().cast()));
            let ge_v = arch::_mm256_loadu_si256(ge.as_ptr().cast());
            let e = arch::_mm256_max_epi32(arch::_mm256_sub_epi32(horizontal, ge_v), open);
            vertical = arch::_mm256_max_epi32(arch::_mm256_sub_epi32(vertical, ge_v), open);
            curr_h[r] = h;
            curr_e[r] = arch::_mm256_blendv_epi8(neg, e, mask);
            vertical = arch::_mm256_blendv_epi8(neg, vertical, mask);
            best = arch::_mm256_max_epi32(best, h);
        }
        std::mem::swap(&mut prev_h, &mut curr_h);
        std::mem::swap(&mut prev_e, &mut curr_e);
    }
    let mut raw = [0i32; LANES];
    arch::_mm256_storeu_si256(raw.as_mut_ptr().cast(), best);
    raw[..targets.len()].to_vec()
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn trace_batch_avx2_impl(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> Vec<SwResult> {
    const LANES: usize = 8;
    const NEG: i32 = i32::MIN / 4;
    const ACTIVE: u8 = 1;
    const GAP_V: u8 = 1 << 1;
    const GAP_H: u8 = 1 << 2;
    const OPEN_V: u8 = 1 << 3;
    const OPEN_H: u8 = 1 << 4;

    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin) as usize)
        .max()
        .unwrap_or(0);
    let max_subject_len = targets
        .iter()
        .map(|target| target.subject.len())
        .max()
        .unwrap_or(0);
    let zero = arch::_mm256_setzero_si256();
    let neg = arch::_mm256_set1_epi32(NEG);
    let mut prev_h = vec![zero; band];
    let mut curr_h = vec![zero; band];
    let mut prev_e = vec![neg; band + 1];
    let mut curr_e = vec![neg; band + 1];
    let mut trace: Vec<Vec<u8>> = targets
        .iter()
        .map(|target| {
            vec![0; (target.subject.len() + 1) * ((target.d_end - target.d_begin) as usize + 2)]
        })
        .collect();
    let mut best_score = [0i32; LANES];
    let mut best_i = [0usize; LANES];
    let mut best_j = [0usize; LANES];

    for j0 in 0..max_subject_len {
        curr_h.fill(zero);
        curr_e.fill(neg);
        let mut vertical = neg;
        for r in 0..band {
            let mut valid = [0i32; LANES];
            let mut subst = [0i32; LANES];
            let mut go = [0i32; LANES];
            let mut ge = [0i32; LANES];
            let mut qpos_lane = [0usize; LANES];
            let mut trace_idx = [0usize; LANES];
            for lane in 0..targets.len() {
                let target = targets[lane];
                let width = (target.d_end - target.d_begin) as usize;
                let q = target.d_begin + j0 as i32 + r as i32;
                if j0 >= target.subject.len() || r >= width || q < 0 || q >= query.len() as i32 {
                    continue;
                }
                valid[lane] = -1;
                let qpos = q as usize;
                qpos_lane[lane] = qpos;
                let lower = (target.d_begin + j0 as i32).max(0) as usize;
                let band_rows = width + 2;
                trace_idx[lane] = (j0 + 1) * band_rows + (qpos + 1 - (lower + 1));
                let ql = query[qpos];
                let sl = target.subject[j0];
                subst[lane] = if let Some(matrix) = target.matrix {
                    matrix.scores[(sl & LETTER_MASK) as usize * 32 + (ql & LETTER_MASK) as usize]
                        as i32
                } else {
                    score_matrix.score(ql, sl) + query_cbs.get(qpos).copied().unwrap_or(0) as i32
                };
                let scale = target.matrix_scale.max(1);
                go[lane] = (score_matrix.gap_open() + score_matrix.gap_extend()) * scale;
                ge[lane] = score_matrix.gap_extend() * scale;
            }
            let mask = arch::_mm256_loadu_si256(valid.as_ptr().cast());
            let substitution = arch::_mm256_loadu_si256(subst.as_ptr().cast());
            let go_v = arch::_mm256_loadu_si256(go.as_ptr().cast());
            let ge_v = arch::_mm256_loadu_si256(ge.as_ptr().cast());
            let diag = arch::_mm256_add_epi32(prev_h[r], substitution);
            let horizontal = prev_e[r + 1];
            let mut score = arch::_mm256_max_epi32(diag, horizontal);
            score = arch::_mm256_max_epi32(score, vertical);
            score = arch::_mm256_max_epi32(score, zero);
            score = arch::_mm256_and_si256(score, mask);
            let open = arch::_mm256_sub_epi32(score, go_v);
            let next_horizontal =
                arch::_mm256_max_epi32(arch::_mm256_sub_epi32(horizontal, ge_v), open);
            let next_vertical =
                arch::_mm256_max_epi32(arch::_mm256_sub_epi32(vertical, ge_v), open);
            curr_h[r] = score;
            curr_e[r] = arch::_mm256_blendv_epi8(neg, next_horizontal, mask);

            let active_bits = arch::_mm256_movemask_epi8(arch::_mm256_cmpgt_epi32(score, zero));
            let gap_v_bits = arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi32(score, vertical));
            let gap_h_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi32(score, horizontal));
            let open_v_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi32(next_vertical, open));
            let open_h_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi32(next_horizontal, open));

            let mut scores = [0i32; LANES];
            arch::_mm256_storeu_si256(scores.as_mut_ptr().cast(), score);
            for lane in 0..targets.len() {
                if valid[lane] == 0 {
                    continue;
                }
                let s = scores[lane];
                let bit = 1 << (lane * 4);
                trace[lane][trace_idx[lane]] = u8::from(active_bits & bit != 0) * ACTIVE
                    | u8::from(gap_v_bits & bit != 0) * GAP_V
                    | u8::from(gap_h_bits & bit != 0) * GAP_H
                    | u8::from(open_v_bits & bit != 0) * OPEN_V
                    | u8::from(open_h_bits & bit != 0) * OPEN_H;
                let j = j0 + 1;
                if s > best_score[lane] || (s == best_score[lane] && j == best_j[lane]) {
                    best_score[lane] = s;
                    best_i[lane] = qpos_lane[lane] + 1;
                    best_j[lane] = j;
                }
            }
            vertical = arch::_mm256_blendv_epi8(neg, next_vertical, mask);
        }
        std::mem::swap(&mut prev_h, &mut curr_h);
        std::mem::swap(&mut prev_e, &mut curr_e);
    }

    targets
        .iter()
        .enumerate()
        .map(|(lane, target)| {
            if best_score[lane] == 0 {
                return SwResult::default();
            }
            let width = (target.d_end - target.d_begin) as usize;
            let band_rows = width + 2;
            let trace_at = |i: usize, j: usize| -> u8 {
                if j == 0 || j > target.subject.len() || i == 0 {
                    return 0;
                }
                let lower = (target.d_begin + j as i32 - 1).max(0) as usize + 1;
                if i < lower || i - lower >= band_rows {
                    0
                } else {
                    trace[lane][j * band_rows + i - lower]
                }
            };
            let mut i = best_i[lane];
            let mut j = best_j[lane];
            let mut result = SwResult {
                score: best_score[lane],
                query_end: i as i32,
                subject_end: j as i32,
                ..Default::default()
            };
            let mut ops = Vec::new();
            while i > 0 && j > 0 && trace_at(i, j) & ACTIVE != 0 {
                let flags = trace_at(i, j);
                if flags & GAP_V != 0 {
                    let mut n = 0;
                    loop {
                        n += 1;
                        i -= 1;
                        if i == 0 || trace_at(i, j) & OPEN_V != 0 {
                            break;
                        }
                    }
                    ops.push((EditOperation::Insertion, n));
                    result.gap_openings += 1;
                    result.gaps += n;
                    result.length += n;
                } else if flags & GAP_H != 0 {
                    let mut n = 0;
                    loop {
                        n += 1;
                        j -= 1;
                        if j == 0 || trace_at(i, j) & OPEN_H != 0 {
                            break;
                        }
                    }
                    ops.push((EditOperation::Deletion, n));
                    result.gap_openings += 1;
                    result.gaps += n;
                    result.length += n;
                } else {
                    if (query[i - 1] & LETTER_MASK) == (target.subject[j - 1] & LETTER_MASK) {
                        ops.push((EditOperation::Match, 1));
                        result.identities += 1;
                    } else {
                        ops.push((EditOperation::Substitution, 1));
                        result.mismatches += 1;
                    }
                    result.length += 1;
                    i -= 1;
                    j -= 1;
                }
            }
            result.query_begin = i as i32;
            result.subject_begin = j as i32;
            ops.reverse();
            result.operations = ops;
            result
        })
        .collect()
}

#[cfg(all(test, any(target_arch = "x86", target_arch = "x86_64")))]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;
    use crate::dp::swipe::{banded_sw_cbs_range, TracebackScratch};
    use crate::stats::cbs::TargetMatrix;
    use std::sync::Arc;

    fn coalesced_operations(operations: &[(EditOperation, i32)]) -> Vec<(EditOperation, i32)> {
        let mut out = Vec::new();
        for &(op, count) in operations {
            if let Some((last_op, last_count)) = out.last_mut() {
                if *last_op == op {
                    *last_count += count;
                    continue;
                }
            }
            out.push((op, count));
        }
        out
    }

    fn assert_sw_result_eq(actual: &SwResult, expected: &SwResult, context: &str) {
        assert_eq!(actual.score, expected.score, "score: {context}");
        assert_eq!(
            (
                actual.query_begin,
                actual.query_end,
                actual.subject_begin,
                actual.subject_end,
                actual.length,
                actual.identities,
                actual.mismatches,
                actual.gap_openings,
                actual.gaps,
            ),
            (
                expected.query_begin,
                expected.query_end,
                expected.subject_begin,
                expected.subject_end,
                expected.length,
                expected.identities,
                expected.mismatches,
                expected.gap_openings,
                expected.gaps,
            ),
            "coordinates/counts: {context}"
        );
        assert_eq!(
            coalesced_operations(&actual.operations),
            coalesced_operations(&expected.operations),
            "operations: {context}"
        );
    }

    #[test]
    fn randomized_traceback_matches_scalar() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0xa5a5_1234u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        for _ in 0..64 {
            let qlen = 4 + next() as usize % 53;
            let query: Vec<Letter> = (0..qlen).map(|_| (next() % 20) as Letter).collect();
            let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
            let count = 1 + next() as usize % 8;
            let subjects: Vec<Vec<Letter>> = (0..count)
                .map(|_| {
                    let len = 3 + next() as usize % 57;
                    (0..len).map(|_| (next() % 20) as Letter).collect()
                })
                .collect();
            let bands: Vec<(i32, i32)> = (0..count)
                .map(|_| {
                    let begin = (next() % 19) as i32 - 9;
                    (begin, begin + 1 + (next() % 25) as i32)
                })
                .collect();
            let targets: Vec<_> = subjects
                .iter()
                .zip(&bands)
                .map(|(subject, &(d_begin, d_end))| TraceTarget {
                    subject,
                    d_begin,
                    d_end,
                    matrix: None,
                    matrix_scale: 1,
                })
                .collect();
            let got = trace_batch_avx2(&query, &targets, &matrix, &cbs).unwrap();
            for lane in 0..count {
                let mut scratch = TracebackScratch::default();
                let expected = banded_sw_cbs_range(
                    &query,
                    &subjects[lane],
                    bands[lane].0,
                    bands[lane].1,
                    &matrix,
                    &cbs,
                    None,
                    1,
                    &mut scratch,
                );
                assert_eq!(
                    got[lane].score, expected.score,
                    "lane={lane} band={:?}",
                    bands[lane]
                );
                assert_eq!(
                    (
                        got[lane].query_begin,
                        got[lane].query_end,
                        got[lane].subject_begin,
                        got[lane].subject_end
                    ),
                    (
                        expected.query_begin,
                        expected.query_end,
                        expected.subject_begin,
                        expected.subject_end
                    ),
                    "lane={lane} band={:?}",
                    bands[lane]
                );
                assert_eq!(
                    got[lane].operations, expected.operations,
                    "lane={lane} band={:?}",
                    bands[lane]
                );
            }
        }
    }

    #[test]
    fn avx2_trace_treats_seed_mask_as_lookup_only_metadata() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query = vec![17; 12];
        let plain = vec![17; 12];
        let marked: Vec<_> = plain.iter().map(|&letter| letter | SEED_MASK).collect();
        let targets = [
            TraceTarget {
                subject: &plain,
                d_begin: 0,
                d_end: 1,
                matrix: None,
                matrix_scale: 1,
            },
            TraceTarget {
                subject: &marked,
                d_begin: 0,
                d_end: 1,
                matrix: None,
                matrix_scale: 1,
            },
        ];
        let scored = score_adjusted_batch_avx2(&query, &targets, &matrix, &[]).unwrap();
        assert_eq!(scored[0], scored[1]);
        let traced = trace_batch_avx2(&query, &targets, &matrix, &[]).unwrap();
        assert_sw_result_eq(&traced[0], &traced[1], "plain vs SEED_MASK");

        let adjusted = TargetMatrix::new(
            (0..1024)
                .map(|index| if index / 32 == index % 32 { 5 } else { -3 })
                .collect(),
            -3,
            5,
        );
        let adjusted_targets = targets.map(|target| TraceTarget {
            matrix: Some(&adjusted),
            ..target
        });
        let scored = score_adjusted_batch_avx2(&query, &adjusted_targets, &matrix, &[]).unwrap();
        assert_eq!(scored[0], scored[1]);
        let traced = trace_batch_avx2(&query, &adjusted_targets, &matrix, &[]).unwrap();
        assert_sw_result_eq(
            &traced[0],
            &traced[1],
            "plain vs SEED_MASK with adjusted matrix",
        );
    }

    #[test]
    fn randomized_narrow_traceback_matches_scalar() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0x1984_5eed_1234_abcd_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        for score_bin in [0usize, 1] {
            let max_lanes = if score_bin == 0 { 32 } else { 16 };
            for case in 0..48 {
                let qlen = 8 + next() as usize % 37;
                let query: Vec<Letter> = (0..qlen)
                    .map(|i| {
                        (next() % 20) as Letter | if (i + case) % 11 == 0 { i8::MIN } else { 0 }
                    })
                    .collect();
                let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
                let count = 1 + next() as usize % max_lanes;
                let mut subjects: Vec<Vec<Letter>> = (0..count)
                    .map(|lane| {
                        let len = 5 + next() as usize % 39;
                        (0..len)
                            .map(|i| {
                                (next() % 20) as Letter
                                    | if (i + lane + case) % 13 == 0 {
                                        i8::MIN
                                    } else {
                                        0
                                    }
                            })
                            .collect()
                    })
                    .collect();
                if case % 3 == 0 {
                    subjects.iter_mut().for_each(|subject| subject.reverse());
                }
                let bands: Vec<(i32, i32)> = (0..count)
                    .enumerate()
                    .map(|(lane, _)| {
                        if case % 4 == 0 {
                            (-(subjects[lane].len() as i32 - 1), qlen as i32)
                        } else {
                            let begin = (next() % 13) as i32 - 6;
                            (begin, begin + 1 + (next() % 17) as i32)
                        }
                    })
                    .collect();
                let targets: Vec<_> = subjects
                    .iter()
                    .zip(&bands)
                    .map(|(subject, &(d_begin, d_end))| TraceTarget {
                        subject,
                        d_begin,
                        d_end,
                        matrix: None,
                        matrix_scale: 1,
                    })
                    .collect();
                let got = trace_batch_tier_avx2(&query, &targets, &matrix, &cbs, score_bin, false)
                    .expect("AVX2 narrow tier");
                for lane in 0..count {
                    if got.overflow_mask & (1 << lane) != 0 {
                        continue;
                    }
                    let mut scratch = TracebackScratch::default();
                    let expected = banded_sw_cbs_range(
                        &query,
                        &subjects[lane],
                        bands[lane].0,
                        bands[lane].1,
                        &matrix,
                        &cbs,
                        None,
                        1,
                        &mut scratch,
                    );
                    assert_sw_result_eq(
                        &got.results[lane],
                        &expected,
                        &format!("bin={score_bin} case={case} lane={lane}"),
                    );
                }
            }
        }
    }

    #[test]
    fn narrow_overflow_promotes_to_exact_next_tier() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let (letter, score) = (0i8..20)
            .map(|letter| (letter, matrix.score(letter, letter)))
            .max_by_key(|&(_, score)| score)
            .unwrap();
        assert!(score > 0);

        let byte_len = 255usize.div_ceil(score as usize) + 2;
        let query = vec![letter; byte_len];
        let subject = query.clone();
        let target = TraceTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
            matrix: None,
            matrix_scale: 1,
        };
        let byte = trace_batch_tier_avx2(&query, &[target], &matrix, &[], 0, false).unwrap();
        assert_eq!(byte.overflow_mask, 1);
        let word = trace_batch_tier_avx2(&query, &[target], &matrix, &[], 1, false).unwrap();
        assert_eq!(word.overflow_mask, 0);
        let mut scratch = TracebackScratch::default();
        let expected =
            banded_sw_cbs_range(&query, &subject, 0, 1, &matrix, &[], None, 1, &mut scratch);
        assert_sw_result_eq(&word.results[0], &expected, "i8 -> i16 promotion");

        let word_len = 65_535usize.div_ceil(score as usize) + 2;
        let query = vec![letter; word_len];
        let subject = query.clone();
        let target = TraceTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
            matrix: None,
            matrix_scale: 1,
        };
        let word = trace_batch_tier_avx2(&query, &[target], &matrix, &[], 1, false).unwrap();
        assert_eq!(word.overflow_mask, 1);
        let exact = trace_batch_tier_avx2(&query, &[target], &matrix, &[], 2, false).unwrap();
        assert_eq!(exact.overflow_mask, 0);
        let mut scratch = TracebackScratch::default();
        let expected =
            banded_sw_cbs_range(&query, &subject, 0, 1, &matrix, &[], None, 1, &mut scratch);
        assert_sw_result_eq(&exact.results[0], &expected, "i16 -> i32 promotion");
    }

    #[test]
    fn narrow_traceback_preserves_equal_score_tie_coordinates() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let (letter, _) = (0i8..20)
            .map(|letter| (letter, matrix.score(letter, letter)))
            .max_by_key(|&(_, score)| score)
            .unwrap();
        // The isolated identical letters create several equal local maxima.
        // Scalar SWIPE's last-column tie rule must remain byte-exact.
        let query = vec![letter, 1, 1, letter, 1, 1, letter];
        let subjects = [
            vec![letter],
            vec![letter, 2, 2, letter],
            vec![letter | i8::MIN, 3, 3, letter],
        ];
        let targets: Vec<_> = subjects
            .iter()
            .map(|subject| TraceTarget {
                subject,
                d_begin: -(subject.len() as i32 - 1),
                d_end: query.len() as i32,
                matrix: None,
                matrix_scale: 1,
            })
            .collect();
        for score_bin in [0, 1] {
            let got =
                trace_batch_tier_avx2(&query, &targets, &matrix, &[], score_bin, false).unwrap();
            assert_eq!(got.overflow_mask, 0);
            for lane in 0..targets.len() {
                let mut scratch = TracebackScratch::default();
                let expected = banded_sw_cbs_range(
                    &query,
                    &subjects[lane],
                    targets[lane].d_begin,
                    targets[lane].d_end,
                    &matrix,
                    &[],
                    None,
                    1,
                    &mut scratch,
                );
                assert_sw_result_eq(
                    &got.results[lane],
                    &expected,
                    &format!("tie bin={score_bin} lane={lane}"),
                );
            }
        }
    }

    #[test]
    fn adjusted_score_avx2_matches_scalar_with_lane_scales() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query: Vec<Letter> = (0..57).map(|i| (i * 7 % 20) as Letter).collect();
        let matrices: Vec<_> = (0..8)
            .map(|lane| {
                let scores = (0..1024)
                    .map(|idx| {
                        let q = idx % 32;
                        let s = idx / 32;
                        if q == s {
                            3 + lane as i8 % 3
                        } else {
                            -2
                        }
                    })
                    .collect();
                Arc::new(TargetMatrix::new(scores, -2, 5))
            })
            .collect();
        let subjects: Vec<Vec<Letter>> = (0..8)
            .map(|lane| {
                (0..(19 + lane * 3))
                    .map(|i| ((i * 11 + lane) % 20) as Letter)
                    .collect()
            })
            .collect();
        let targets: Vec<_> = subjects
            .iter()
            .enumerate()
            .map(|(lane, subject)| TraceTarget {
                subject,
                d_begin: -(lane as i32 % 5),
                d_end: 12 + lane as i32,
                matrix: Some(matrices[lane].as_ref()),
                matrix_scale: 1 + lane as i32 % 2,
            })
            .collect();
        let got = score_adjusted_batch_avx2(&query, &targets, &matrix, &[]).unwrap();
        for lane in 0..targets.len() {
            let mut scratch = TracebackScratch::default();
            let expected = banded_sw_cbs_range(
                &query,
                &subjects[lane],
                targets[lane].d_begin,
                targets[lane].d_end,
                &matrix,
                &[],
                Some(matrices[lane].as_ref()),
                targets[lane].matrix_scale,
                &mut scratch,
            );
            assert_eq!(got[lane], expected.score, "lane={lane}");
        }
    }

    #[test]
    fn adjusted_narrow_mixed_lanes_match_scalar_score_and_trace() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query: Vec<Letter> = (0..41).map(|i| (i * 7 % 20) as Letter).collect();
        let cbs: Vec<i8> = (0..query.len()).map(|i| (i % 5) as i8 - 2).collect();
        let adjusted: Vec<_> = (0..16)
            .map(|lane| {
                let scores = (0..1024)
                    .map(|idx| {
                        let q = idx % 32;
                        let s = idx / 32;
                        if q == s {
                            3 + lane as i8 % 4
                        } else {
                            -2
                        }
                    })
                    .collect();
                Arc::new(TargetMatrix::new(scores, -2, 6))
            })
            .collect();

        for reverse in [false, true] {
            let mut subjects: Vec<Vec<Letter>> = (0..16)
                .map(|lane| {
                    (0..(11 + lane % 9))
                        .map(|i| ((i * 11 + lane * 3) % 20) as Letter)
                        .collect()
                })
                .collect();
            if reverse {
                subjects.iter_mut().for_each(|s| s.reverse());
            }
            for full in [false, true] {
                let targets: Vec<_> = subjects
                    .iter()
                    .enumerate()
                    .map(|(lane, subject)| TraceTarget {
                        subject,
                        d_begin: if full {
                            -(subject.len() as i32 - 1)
                        } else {
                            -(lane as i32 % 4)
                        },
                        d_end: if full {
                            query.len() as i32
                        } else {
                            9 + lane as i32 % 7
                        },
                        matrix: (lane % 3 != 0).then_some(adjusted[lane].as_ref()),
                        matrix_scale: if lane % 3 != 0 {
                            1 + lane as i32 % 3
                        } else {
                            1
                        },
                    })
                    .collect();
                for score_bin in [0usize, 1] {
                    let scores = if score_bin == 0 {
                        super::super::simd_adjusted_narrow::score_batch_avx2_i8(
                            &query, &targets, &matrix, &cbs, false,
                        )
                    } else {
                        super::super::simd_adjusted_narrow::score_batch_avx2_i16(
                            &query, &targets, &matrix, &cbs, false,
                        )
                    }
                    .unwrap();
                    for lane in 0..targets.len() {
                        if scores.overflow_mask & (1 << lane) != 0 {
                            continue;
                        }
                        let mut scratch = TracebackScratch::default();
                        let expected = banded_sw_cbs_range(
                            &query,
                            &subjects[lane],
                            targets[lane].d_begin,
                            targets[lane].d_end,
                            &matrix,
                            &cbs,
                            targets[lane].matrix,
                            targets[lane].matrix_scale,
                            &mut scratch,
                        );
                        assert_eq!(
                            scores.scores[lane], expected.score,
                            "score bin={score_bin} lane={lane}"
                        );
                    }

                    // Upstream banded SWIPE traces adjusted and ordinary
                    // lanes together; only full-matrix SWIPE rejects an
                    // adjusted traceback.
                    let traced =
                        trace_batch_tier_avx2(&query, &targets, &matrix, &cbs, score_bin, false)
                            .unwrap();
                    for lane in 0..targets.len() {
                        if traced.overflow_mask & (1 << lane) != 0 {
                            continue;
                        }
                        let expected = banded_sw_cbs_range(
                            &query,
                            &subjects[lane],
                            targets[lane].d_begin,
                            targets[lane].d_end,
                            &matrix,
                            &cbs,
                            targets[lane].matrix,
                            targets[lane].matrix_scale,
                            &mut TracebackScratch::default(),
                        );
                        assert_sw_result_eq(
                            &traced.results[lane],
                            &expected,
                            &format!(
                                "trace bin={score_bin} lane={lane} reverse={reverse} full={full}"
                            ),
                        );
                    }

                    // Call the 128-bit backend directly even on an AVX2 host,
                    // so SSE/NEON lane packing and mixed-profile semantics do
                    // not depend only on cross-compilation.
                    let portable_len = if score_bin == 0 { 16 } else { 8 };
                    let portable_targets = &targets[..portable_len];
                    let portable_scores = if score_bin == 0 {
                        super::super::simd_trace_narrow_portable::score_batch_i8(
                            &query,
                            portable_targets,
                            &matrix,
                            &cbs,
                            false,
                        )
                    } else {
                        super::super::simd_trace_narrow_portable::score_batch_i16(
                            &query,
                            portable_targets,
                            &matrix,
                            &cbs,
                            false,
                        )
                    }
                    .unwrap();
                    let portable_trace = if score_bin == 0 {
                        super::super::simd_trace_narrow_portable::trace_batch_i8(
                            &query,
                            portable_targets,
                            &matrix,
                            &cbs,
                            false,
                        )
                    } else {
                        super::super::simd_trace_narrow_portable::trace_batch_i16(
                            &query,
                            portable_targets,
                            &matrix,
                            &cbs,
                            false,
                        )
                    }
                    .unwrap();
                    for lane in 0..portable_len {
                        if portable_scores.overflow_mask & (1 << lane) != 0 {
                            continue;
                        }
                        let mut scratch = TracebackScratch::default();
                        let expected = banded_sw_cbs_range(
                            &query,
                            &subjects[lane],
                            targets[lane].d_begin,
                            targets[lane].d_end,
                            &matrix,
                            &cbs,
                            targets[lane].matrix,
                            targets[lane].matrix_scale,
                            &mut scratch,
                        );
                        assert_eq!(portable_scores.scores[lane], expected.score);
                    }
                    for lane in 0..portable_len {
                        if portable_trace.overflow_mask & (1 << lane) != 0 {
                            continue;
                        }
                        let expected = banded_sw_cbs_range(
                            &query,
                            &subjects[lane],
                            targets[lane].d_begin,
                            targets[lane].d_end,
                            &matrix,
                            &cbs,
                            targets[lane].matrix,
                            targets[lane].matrix_scale,
                            &mut TracebackScratch::default(),
                        );
                        assert_sw_result_eq(
                            &portable_trace.results[lane],
                            &expected,
                            &format!("portable bin={score_bin} lane={lane}"),
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn adjusted_score_saturation_promotes_through_both_bins() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 100_000).unwrap();
        let mut scores = vec![-1i8; 1024];
        scores[0] = 100;
        let adjusted = TargetMatrix::new(scores, -1, 100);

        let query = vec![0; 4];
        let subject = query.clone();
        let target = TraceTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
            matrix: Some(&adjusted),
            matrix_scale: 3,
        };
        let byte_score = super::super::simd_adjusted_narrow::score_batch_avx2_i8(
            &query,
            &[target],
            &matrix,
            &[],
            false,
        )
        .unwrap();
        assert_eq!(byte_score.overflow_mask, 1);

        let query = vec![0; 700];
        let subject = query.clone();
        let target = TraceTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
            matrix: Some(&adjusted),
            matrix_scale: 3,
        };
        let word_score = super::super::simd_adjusted_narrow::score_batch_avx2_i16(
            &query,
            &[target],
            &matrix,
            &[],
            false,
        )
        .unwrap();
        assert_eq!(word_score.overflow_mask, 1);
        let exact = score_adjusted_batch_avx2(&query, &[target], &matrix, &[]).unwrap();
        assert_eq!(exact[0], 70_000);
    }

    #[test]
    fn semi_global_narrow_trace_and_adjusted_score_use_delta_zero() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        // Pinned regression from the score kernels: clamping intermediate
        // cells to zero changes this DELTA=0 alignment's maximum.
        let query: Vec<Letter> = b"011000010110111010111100000101"
            .iter()
            .map(|&x| (x - b'0') as Letter)
            .collect();
        let subject: Vec<Letter> = b"11010110011001100101010101"
            .iter()
            .map(|&x| (x - b'0') as Letter)
            .collect();
        let target = TraceTarget {
            subject: &subject,
            d_begin: -3,
            d_end: 13,
            matrix: None,
            matrix_scale: 1,
        };
        let local = trace_batch_tier_avx2(&query, &[target], &matrix, &[], 0, false).unwrap();
        assert_eq!(local.overflow_mask, 0);
        let expected_local = banded_sw_cbs_range(
            &query,
            &subject,
            target.d_begin,
            target.d_end,
            &matrix,
            &[],
            None,
            1,
            &mut TracebackScratch::default(),
        );
        assert_sw_result_eq(&local.results[0], &expected_local, "local regression");

        let adjusted = TargetMatrix::new(
            (0..1024)
                .map(|index| matrix.matrix8()[(index % 32) * 32 + index / 32])
                .collect(),
            i8::MIN as i32,
            i8::MAX as i32,
        );
        let adjusted_target = TraceTarget {
            matrix: Some(&adjusted),
            ..target
        };
        let adjusted_score = super::super::simd_adjusted_narrow::score_batch_avx2_i8(
            &query,
            &[adjusted_target],
            &matrix,
            &[],
            true,
        )
        .unwrap();
        assert_eq!(adjusted_score.overflow_mask, 0);
        assert_eq!(adjusted_score.scores, vec![69]);

        // A transcript-producing semi-global lane reaches the same narrow
        // DELTA=0 specialization. Use a monotone exact-match traceback, as
        // upstream's vector traceback reconstructs score from the endpoint.
        let trace_query = vec![17; 5];
        let trace_subject = trace_query.clone();
        let trace_target = TraceTarget {
            subject: &trace_subject,
            d_begin: 0,
            d_end: 1,
            matrix: None,
            matrix_scale: 1,
        };
        let semi =
            trace_batch_tier_avx2(&trace_query, &[trace_target], &matrix, &[], 0, true).unwrap();
        assert_eq!(semi.overflow_mask, 0);
        assert_eq!(semi.results[0].score, 55);
        assert_eq!(semi.results[0].operations, vec![(EditOperation::Match, 5)]);
        let adjusted_trace_target = TraceTarget {
            matrix: Some(&adjusted),
            ..trace_target
        };
        let adjusted_trace = trace_batch_tier_avx2(
            &trace_query,
            &[adjusted_trace_target],
            &matrix,
            &[],
            0,
            true,
        )
        .unwrap();
        assert_eq!(adjusted_trace.overflow_mask, 0);
        assert_eq!(adjusted_trace.results[0].score, 55);
        assert_eq!(
            adjusted_trace.results[0].operations,
            semi.results[0].operations
        );

        let portable_adjusted_score = super::super::simd_trace_narrow_portable::score_batch_i8(
            &query,
            &[adjusted_target],
            &matrix,
            &[],
            true,
        )
        .unwrap();
        assert_eq!(portable_adjusted_score.scores, vec![69]);
        let portable_trace = super::super::simd_trace_narrow_portable::trace_batch_i8(
            &trace_query,
            &[adjusted_trace_target],
            &matrix,
            &[],
            true,
        )
        .unwrap();
        assert_eq!(portable_trace.results[0].score, 55);
        assert_eq!(portable_trace.results[0].length, 5);
        assert!(portable_trace.results[0]
            .operations
            .iter()
            .all(|&(op, count)| op == EditOperation::Match && count > 0));
    }

    #[test]
    fn local_narrow_trace_stops_at_nonzero_alignment_start() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query = vec![0, 0, 17, 17, 17, 17, 17];
        let subject = vec![1, 1, 17, 17, 17, 17, 17];
        let target = TraceTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
            matrix: None,
            matrix_scale: 1,
        };
        let expected = banded_sw_cbs_range(
            &query,
            &subject,
            0,
            1,
            &matrix,
            &[],
            None,
            1,
            &mut TracebackScratch::default(),
        );
        assert!(expected.query_begin > 0);
        assert!(expected.subject_begin > 0);

        for score_bin in [0, 1] {
            let got =
                trace_batch_tier_avx2(&query, &[target], &matrix, &[], score_bin, false).unwrap();
            assert_eq!(got.overflow_mask, 0);
            assert_sw_result_eq(
                &got.results[0],
                &expected,
                &format!("AVX2 local bin {score_bin}"),
            );
        }
        let portable = super::super::simd_trace_narrow_portable::trace_batch_i8(
            &query,
            &[target],
            &matrix,
            &[],
            false,
        )
        .unwrap();
        assert_eq!(portable.overflow_mask, 0);
        assert_sw_result_eq(&portable.results[0], &expected, "portable local byte");
    }
}
