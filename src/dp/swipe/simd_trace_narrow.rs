//! Narrow AVX2 traceback tiers used by SWIPE score bins 0 and 1.
//!
//! DIAMOND represents a local-alignment score as an unsigned distance from
//! the signed minimum of the lane type.  This gives the byte and word tiers
//! score ranges 0..=255 and 0..=65535 respectively.  A lane that reaches the
//! signed maximum is promoted to the next score bin and recomputed there.

use super::simd_trace::TraceTarget;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::basic::value::LETTER_MASK;
use crate::dp::smith_waterman::SwResult;
use crate::stats::score_matrix::ScoreMatrix;

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

pub struct NarrowTraceBatch {
    pub results: Vec<SwResult>,
    pub overflow_mask: u32,
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
const ACTIVE: u8 = 1;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
const GAP_V: u8 = 1 << 1;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
const GAP_H: u8 = 1 << 2;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
const OPEN_V: u8 = 1 << 3;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
const OPEN_H: u8 = 1 << 4;

pub fn trace_batch_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> Option<NarrowTraceBatch> {
    if targets.is_empty()
        || targets.len() > 32
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
    {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: AVX2 is runtime-detected and all lane accesses below
            // are bounded by `targets.len()`.
            return Some(unsafe { trace_i8_impl(query, targets, score_matrix, query_cbs) });
        }
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    let _ = score_matrix;
    None
}

pub fn trace_batch_i16(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> Option<NarrowTraceBatch> {
    if targets.is_empty()
        || targets.len() > 16
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
    {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: AVX2 is runtime-detected and all lane accesses below
            // are bounded by `targets.len()`.
            return Some(unsafe { trace_i16_impl(query, targets, score_matrix, query_cbs) });
        }
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    let _ = score_matrix;
    None
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn trace_i8_impl(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> NarrowTraceBatch {
    const LANES: usize = 32;
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    let i1 = targets
        .iter()
        .map(|target| (target.d_end - 1).max(0))
        .min()
        .unwrap_or(0);
    let i0 = i1 + 1 - band as i32;
    let mut subject_start = [0i32; LANES];
    let mut band_offset = [0usize; LANES];
    let mut cols = 0usize;
    for (lane, target) in targets.iter().enumerate() {
        let expanded_begin = target.d_end - band as i32;
        subject_start[lane] = i1 - (target.d_end - 1);
        band_offset[lane] = (target.d_begin - expanded_begin).max(0) as usize;
        let subject_end =
            ((query.len() as i32 - 1 - expanded_begin).min(target.subject.len() as i32 - 1) + 1)
                .max(0);
        cols = cols.max((subject_end - subject_start[lane]).max(0) as usize);
    }
    let zero = arch::_mm256_set1_epi8(i8::MIN);
    let max_value = arch::_mm256_set1_epi8(i8::MAX);
    let mut gap_open = [0i8; LANES];
    let mut gap_extend = [0i8; LANES];
    let mut overflow_mask = 0u32;
    for lane in 0..targets.len() {
        let scale = targets[lane].matrix_scale.max(1);
        let go = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=63).contains(&go) || !(0..=63).contains(&ge) {
            overflow_mask |= 1 << lane;
        }
        gap_open[lane] = go.clamp(0, 63) as i8;
        gap_extend[lane] = ge.clamp(0, 63) as i8;
    }
    let go = arch::_mm256_loadu_si256(gap_open.as_ptr().cast());
    let ge = arch::_mm256_loadu_si256(gap_extend.as_ptr().cast());
    let mut prev_h = vec![zero; band];
    let mut curr_h = vec![zero; band];
    let mut prev_e = vec![zero; band + 1];
    let mut curr_e = vec![zero; band + 1];
    let mut trace = allocate_trace(targets);
    let mut best_score = [i8::MIN; LANES];
    let mut best_i = [0usize; LANES];
    let mut best_j = [0usize; LANES];
    let mut row_masks = Vec::with_capacity(band);
    for r in 0..band {
        let mut lanes = [0i8; LANES];
        for lane in 0..targets.len() {
            if r >= band_offset[lane] {
                lanes[lane] = -1;
            }
        }
        row_masks.push(arch::_mm256_loadu_si256(lanes.as_ptr().cast()));
    }
    for column in 0..cols {
        curr_h.fill(zero);
        curr_e.fill(zero);
        let mut vertical = zero;
        let mut active = [0i8; LANES];
        let mut subject_positions = [0usize; LANES];
        for lane in 0..targets.len() {
            let position = subject_start[lane] + column as i32;
            if position < 0 || position >= targets[lane].subject.len() as i32 {
                continue;
            }
            active[lane] = -1;
            subject_positions[lane] = position as usize;
        }
        let active_v = arch::_mm256_loadu_si256(active.as_ptr().cast());
        for r in 0..band {
            let q = i0 + column as i32 + r as i32;
            if q < 0 || q >= query.len() as i32 {
                vertical = zero;
                continue;
            }
            let qpos = q as usize;
            let mask = arch::_mm256_and_si256(row_masks[r], active_v);
            let mask_bits = arch::_mm256_movemask_epi8(mask) as u32;
            let query_letter = (query[qpos] & LETTER_MASK) as usize;
            let mut subst = [0i8; LANES];
            for lane in 0..targets.len() {
                if mask_bits & (1 << lane) == 0 {
                    continue;
                }
                let sl = targets[lane].subject[subject_positions[lane]];
                let raw = if sl & crate::basic::value::SEED_MASK != 0 {
                    0
                } else if let Some(adjusted) = targets[lane].matrix {
                    adjusted.scores[(sl & LETTER_MASK) as usize * 32 + query_letter]
                } else {
                    matrix.matrix8()[query_letter * 32 + (sl & LETTER_MASK) as usize]
                };
                let value = i16::from(raw)
                    + if targets[lane].matrix.is_none() {
                        i16::from(query_cbs.get(qpos).copied().unwrap_or(0))
                    } else {
                        0
                    };
                if !(i8::MIN as i16..=i8::MAX as i16).contains(&value) {
                    overflow_mask |= 1 << lane;
                }
                subst[lane] = value.clamp(i8::MIN as i16, i8::MAX as i16) as i8;
            }
            let substitution = arch::_mm256_loadu_si256(subst.as_ptr().cast());
            let diag = arch::_mm256_adds_epi8(prev_h[r], substitution);
            let horizontal = prev_e[r + 1];
            let mut score = arch::_mm256_max_epi8(diag, horizontal);
            score = arch::_mm256_max_epi8(score, vertical);
            score = arch::_mm256_max_epi8(score, zero);
            score = arch::_mm256_blendv_epi8(zero, score, mask);
            let open = arch::_mm256_subs_epi8(score, go);
            let next_horizontal =
                arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, ge), open);
            let next_vertical = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical, ge), open);
            curr_h[r] = score;
            curr_e[r] = arch::_mm256_blendv_epi8(zero, next_horizontal, mask);

            let active_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpgt_epi8(score, zero)) as u32;
            let gap_v_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(score, vertical)) as u32;
            let gap_h_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(score, horizontal)) as u32;
            let open_v_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(next_vertical, open)) as u32;
            let open_h_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(next_horizontal, open)) as u32;
            let saturated =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(score, max_value)) as u32;
            overflow_mask |= saturated;

            let mut scores = [i8::MIN; LANES];
            arch::_mm256_storeu_si256(scores.as_mut_ptr().cast(), score);
            for lane in 0..targets.len() {
                let bit = 1u32 << lane;
                if mask_bits & bit == 0 {
                    continue;
                }
                let target = targets[lane];
                let subject_pos = subject_positions[lane];
                let width = (target.d_end - target.d_begin).max(0) as usize;
                let lower = (target.d_begin + subject_pos as i32).max(0) as usize;
                let trace_idx = (subject_pos + 1) * (width + 2) + qpos - lower;
                trace[lane][trace_idx] = u8::from(active_bits & bit != 0) * ACTIVE
                    | u8::from(gap_v_bits & bit != 0) * GAP_V
                    | u8::from(gap_h_bits & bit != 0) * GAP_H
                    | u8::from(open_v_bits & bit != 0) * OPEN_V
                    | u8::from(open_h_bits & bit != 0) * OPEN_H;
                let j = subject_pos + 1;
                if scores[lane] > best_score[lane]
                    || (scores[lane] == best_score[lane] && j == best_j[lane])
                {
                    best_score[lane] = scores[lane];
                    best_i[lane] = qpos + 1;
                    best_j[lane] = j;
                }
            }
            vertical = arch::_mm256_blendv_epi8(zero, next_vertical, mask);
        }
        std::mem::swap(&mut prev_h, &mut curr_h);
        std::mem::swap(&mut prev_e, &mut curr_e);
    }

    let best: Vec<i32> = best_score[..targets.len()]
        .iter()
        .map(|&score| i32::from(score) - i32::from(i8::MIN))
        .collect();
    NarrowTraceBatch {
        results: finish_results(query, targets, &trace, &best, &best_i, &best_j),
        overflow_mask,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn trace_i16_impl(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> NarrowTraceBatch {
    const LANES: usize = 16;
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    let i1 = targets
        .iter()
        .map(|target| (target.d_end - 1).max(0))
        .min()
        .unwrap_or(0);
    let i0 = i1 + 1 - band as i32;
    let mut subject_start = [0i32; LANES];
    let mut band_offset = [0usize; LANES];
    let mut cols = 0usize;
    for (lane, target) in targets.iter().enumerate() {
        let expanded_begin = target.d_end - band as i32;
        subject_start[lane] = i1 - (target.d_end - 1);
        band_offset[lane] = (target.d_begin - expanded_begin).max(0) as usize;
        let subject_end =
            ((query.len() as i32 - 1 - expanded_begin).min(target.subject.len() as i32 - 1) + 1)
                .max(0);
        cols = cols.max((subject_end - subject_start[lane]).max(0) as usize);
    }
    let zero = arch::_mm256_set1_epi16(i16::MIN);
    let max_value = arch::_mm256_set1_epi16(i16::MAX);
    let mut gap_open = [0i16; LANES];
    let mut gap_extend = [0i16; LANES];
    let mut overflow_mask = 0u32;
    for lane in 0..targets.len() {
        let scale = targets[lane].matrix_scale.max(1);
        let go = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=16_000).contains(&go) || !(0..=16_000).contains(&ge) {
            overflow_mask |= 1 << lane;
        }
        gap_open[lane] = go.clamp(0, 16_000) as i16;
        gap_extend[lane] = ge.clamp(0, 16_000) as i16;
    }
    let go = arch::_mm256_loadu_si256(gap_open.as_ptr().cast());
    let ge = arch::_mm256_loadu_si256(gap_extend.as_ptr().cast());
    let mut prev_h = vec![zero; band];
    let mut curr_h = vec![zero; band];
    let mut prev_e = vec![zero; band + 1];
    let mut curr_e = vec![zero; band + 1];
    let mut trace = allocate_trace(targets);
    let mut best_score = [i16::MIN; LANES];
    let mut best_i = [0usize; LANES];
    let mut best_j = [0usize; LANES];
    let mut row_masks = Vec::with_capacity(band);
    for r in 0..band {
        let mut lanes = [0i16; LANES];
        for lane in 0..targets.len() {
            if r >= band_offset[lane] {
                lanes[lane] = -1;
            }
        }
        row_masks.push(arch::_mm256_loadu_si256(lanes.as_ptr().cast()));
    }
    for column in 0..cols {
        curr_h.fill(zero);
        curr_e.fill(zero);
        let mut vertical = zero;
        let mut active = [0i16; LANES];
        let mut subject_positions = [0usize; LANES];
        for lane in 0..targets.len() {
            let position = subject_start[lane] + column as i32;
            if position < 0 || position >= targets[lane].subject.len() as i32 {
                continue;
            }
            active[lane] = -1;
            subject_positions[lane] = position as usize;
        }
        let active_v = arch::_mm256_loadu_si256(active.as_ptr().cast());
        for r in 0..band {
            let q = i0 + column as i32 + r as i32;
            if q < 0 || q >= query.len() as i32 {
                vertical = zero;
                continue;
            }
            let qpos = q as usize;
            let mask = arch::_mm256_and_si256(row_masks[r], active_v);
            let mask_bits = arch::_mm256_movemask_epi8(mask) as u32;
            let query_letter = (query[qpos] & LETTER_MASK) as usize;
            let mut subst = [0i16; LANES];
            for lane in 0..targets.len() {
                if mask_bits & (1 << (2 * lane)) == 0 {
                    continue;
                }
                let sl = targets[lane].subject[subject_positions[lane]];
                let raw = if sl & crate::basic::value::SEED_MASK != 0 {
                    0
                } else if let Some(adjusted) = targets[lane].matrix {
                    adjusted.scores[(sl & LETTER_MASK) as usize * 32 + query_letter]
                } else {
                    matrix.matrix8()[query_letter * 32 + (sl & LETTER_MASK) as usize]
                };
                subst[lane] = i16::from(raw)
                    + if targets[lane].matrix.is_none() {
                        i16::from(query_cbs.get(qpos).copied().unwrap_or(0))
                    } else {
                        0
                    };
            }
            let substitution = arch::_mm256_loadu_si256(subst.as_ptr().cast());
            let diag = arch::_mm256_adds_epi16(prev_h[r], substitution);
            let horizontal = prev_e[r + 1];
            let mut score = arch::_mm256_max_epi16(diag, horizontal);
            score = arch::_mm256_max_epi16(score, vertical);
            score = arch::_mm256_max_epi16(score, zero);
            score = arch::_mm256_blendv_epi8(zero, score, mask);
            let open = arch::_mm256_subs_epi16(score, go);
            let next_horizontal =
                arch::_mm256_max_epi16(arch::_mm256_subs_epi16(horizontal, ge), open);
            let next_vertical = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(vertical, ge), open);
            curr_h[r] = score;
            curr_e[r] = arch::_mm256_blendv_epi8(zero, next_horizontal, mask);

            let active_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpgt_epi16(score, zero)) as u32;
            let gap_v_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(score, vertical)) as u32;
            let gap_h_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(score, horizontal)) as u32;
            let open_v_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(next_vertical, open)) as u32;
            let open_h_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(next_horizontal, open)) as u32;
            let saturated =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(score, max_value)) as u32;

            let mut scores = [i16::MIN; LANES];
            arch::_mm256_storeu_si256(scores.as_mut_ptr().cast(), score);
            for lane in 0..targets.len() {
                let bit = 1u32 << (lane * 2);
                if mask_bits & bit == 0 {
                    continue;
                }
                if saturated & bit != 0 {
                    overflow_mask |= 1 << lane;
                }
                let target = targets[lane];
                let subject_pos = subject_positions[lane];
                let width = (target.d_end - target.d_begin).max(0) as usize;
                let lower = (target.d_begin + subject_pos as i32).max(0) as usize;
                let trace_idx = (subject_pos + 1) * (width + 2) + qpos - lower;
                trace[lane][trace_idx] = u8::from(active_bits & bit != 0) * ACTIVE
                    | u8::from(gap_v_bits & bit != 0) * GAP_V
                    | u8::from(gap_h_bits & bit != 0) * GAP_H
                    | u8::from(open_v_bits & bit != 0) * OPEN_V
                    | u8::from(open_h_bits & bit != 0) * OPEN_H;
                let j = subject_pos + 1;
                if scores[lane] > best_score[lane]
                    || (scores[lane] == best_score[lane] && j == best_j[lane])
                {
                    best_score[lane] = scores[lane];
                    best_i[lane] = qpos + 1;
                    best_j[lane] = j;
                }
            }
            vertical = arch::_mm256_blendv_epi8(zero, next_vertical, mask);
        }
        std::mem::swap(&mut prev_h, &mut curr_h);
        std::mem::swap(&mut prev_e, &mut curr_e);
    }

    let best: Vec<i32> = best_score[..targets.len()]
        .iter()
        .map(|&score| i32::from(score) - i32::from(i16::MIN))
        .collect();
    NarrowTraceBatch {
        results: finish_results(query, targets, &trace, &best, &best_i, &best_j),
        overflow_mask,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn allocate_trace(targets: &[TraceTarget<'_>]) -> Vec<Vec<u8>> {
    targets
        .iter()
        .map(|target| {
            let width = (target.d_end - target.d_begin).max(0) as usize;
            vec![0; (target.subject.len() + 1) * (width + 2)]
        })
        .collect()
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn finish_results<const LANES: usize>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    trace: &[Vec<u8>],
    best_score: &[i32],
    best_i: &[usize; LANES],
    best_j: &[usize; LANES],
) -> Vec<SwResult> {
    targets
        .iter()
        .enumerate()
        .map(|(lane, target)| {
            if best_score[lane] == 0 {
                return SwResult::default();
            }
            let width = (target.d_end - target.d_begin).max(0) as usize;
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
            let mut operations = Vec::new();
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
                    operations.push((EditOperation::Insertion, n));
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
                    operations.push((EditOperation::Deletion, n));
                    result.gap_openings += 1;
                    result.gaps += n;
                    result.length += n;
                } else {
                    if (query[i - 1] & LETTER_MASK) == (target.subject[j - 1] & LETTER_MASK) {
                        operations.push((EditOperation::Match, 1));
                        result.identities += 1;
                    } else {
                        operations.push((EditOperation::Substitution, 1));
                        result.mismatches += 1;
                    }
                    result.length += 1;
                    i -= 1;
                    j -= 1;
                }
            }
            result.query_begin = i as i32;
            result.subject_begin = j as i32;
            operations.reverse();
            result.operations = operations;
            result
        })
        .collect()
}
