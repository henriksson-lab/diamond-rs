//! 128-bit score-only Smith-Waterman kernels for the main SWIPE path.
//!
//! x86 and x86-64 use eight signed-`i16` SSE2 lanes (with an SSSE3 runtime
//! dispatch target), while AArch64 uses eight signed-`i16` NEON lanes.  The
//! public entry points deliberately have the same score and overflow contract
//! as the AVX2 kernel in [`super::simd_score`], so callers can use this module
//! as the next dispatch tier without changing promotion behavior.

use super::simd_score::{BatchScores, ScoreTarget};
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
use crate::basic::value::LETTER_MASK;
use crate::stats::score_matrix::ScoreMatrix;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64 as arch;
#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type ArchVector = arch::__m128i;
#[cfg(target_arch = "aarch64")]
type ArchVector = arch::int16x8_t;

/// Reusable DP rows for the eight-lane kernels.
#[derive(Default)]
pub struct PortableSimdScoreScratch {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
    prev_h: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
    curr_h: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
    prev_e: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
    curr_e: Vec<ArchVector>,
}

fn penalties(matrix: &ScoreMatrix) -> Option<(i16, i16)> {
    let open = matrix.gap_open().checked_add(matrix.gap_extend())?;
    let extend = matrix.gap_extend();
    // NEG in the kernels is -16384. Keep deliberately unreachable-state
    // subtraction representable, matching the AVX2 tier's accepted range.
    if !(0..=16_000).contains(&open) || !(0..=16_000).contains(&extend) {
        return None;
    }
    Some((open as i16, extend as i16))
}

fn valid_common(query: &[Letter], target_count: usize, cbs: &[i8]) -> bool {
    target_count != 0 && target_count <= 8 && (cbs.is_empty() || cbs.len() >= query.len())
}

/// Score up to eight strictly banded targets using the best available
/// non-AVX2 128-bit `i16` kernel.
///
/// `None` means that the current architecture has no supported backend, the
/// arguments are outside the kernel contract, or the gap penalties are not
/// safely representable. Lanes in `overflow_mask` must be recomputed in i32.
/// `semi_global` selects upstream's `DELTA=0` specialization; local alignment
/// uses the `SHRT_MIN`-biased specialization.
pub fn score_batch_portable_i16(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    scratch: &mut PortableSimdScoreScratch,
) -> Option<BatchScores> {
    if !valid_common(query, targets.len(), cbs) {
        return None;
    }
    let (open, extend) = penalties(matrix)?;

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("ssse3") {
            // SAFETY: SSSE3 was detected at runtime. The implementation only
            // accesses slices after checking their bounds.
            return Some(unsafe {
                if semi_global {
                    score_batch_ssse3_impl::<true>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        open,
                        extend,
                        scratch,
                    )
                } else {
                    score_batch_ssse3_impl::<false>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        open,
                        extend,
                        scratch,
                    )
                }
            });
        }
        if std::arch::is_x86_feature_detected!("sse2") {
            // SAFETY: SSE2 was detected at runtime.
            return Some(unsafe {
                if semi_global {
                    score_batch_sse2_impl::<true>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        open,
                        extend,
                        scratch,
                    )
                } else {
                    score_batch_sse2_impl::<false>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        open,
                        extend,
                        scratch,
                    )
                }
            });
        }
        None
    }

    #[cfg(target_arch = "aarch64")]
    {
        if !std::arch::is_aarch64_feature_detected!("neon") {
            return None;
        }
        // SAFETY: NEON was detected at runtime.
        Some(unsafe {
            if semi_global {
                score_batch_neon_impl::<true>(
                    query,
                    targets,
                    matrix.matrix16(),
                    cbs,
                    open,
                    extend,
                    scratch,
                )
            } else {
                score_batch_neon_impl::<false>(
                    query,
                    targets,
                    matrix.matrix16(),
                    cbs,
                    open,
                    extend,
                    scratch,
                )
            }
        })
    }

    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
    {
        let _ = (
            query,
            targets,
            matrix,
            cbs,
            semi_global,
            scratch,
            open,
            extend,
        );
        None
    }
}

/// Score up to eight complete Smith-Waterman matrices using the best
/// available non-AVX2 128-bit `i16` kernel.
pub fn score_full_batch_portable_i16(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    scratch: &mut PortableSimdScoreScratch,
) -> Option<BatchScores> {
    if !valid_common(query, targets.len(), cbs) {
        return None;
    }
    let (open, extend) = penalties(matrix)?;

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("ssse3") {
            // SAFETY: SSSE3 was detected at runtime.
            return Some(unsafe {
                if semi_global {
                    score_full_ssse3_impl::<true>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        open,
                        extend,
                        scratch,
                    )
                } else {
                    score_full_ssse3_impl::<false>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        open,
                        extend,
                        scratch,
                    )
                }
            });
        }
        if std::arch::is_x86_feature_detected!("sse2") {
            // SAFETY: SSE2 was detected at runtime.
            return Some(unsafe {
                if semi_global {
                    score_full_sse2_impl::<true>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        open,
                        extend,
                        scratch,
                    )
                } else {
                    score_full_sse2_impl::<false>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        open,
                        extend,
                        scratch,
                    )
                }
            });
        }
        None
    }

    #[cfg(target_arch = "aarch64")]
    {
        if !std::arch::is_aarch64_feature_detected!("neon") {
            return None;
        }
        // SAFETY: NEON was detected at runtime.
        Some(unsafe {
            if semi_global {
                score_full_neon_impl::<true>(
                    query,
                    targets,
                    matrix.matrix16(),
                    cbs,
                    open,
                    extend,
                    scratch,
                )
            } else {
                score_full_neon_impl::<false>(
                    query,
                    targets,
                    matrix.matrix16(),
                    cbs,
                    open,
                    extend,
                    scratch,
                )
            }
        })
    }

    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
    {
        let _ = (
            query,
            targets,
            matrix,
            cbs,
            semi_global,
            scratch,
            open,
            extend,
        );
        None
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn score_batch_ssse3_impl<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &[i16; 32 * 32],
    cbs: &[i8],
    open: i16,
    extend: i16,
    scratch: &mut PortableSimdScoreScratch,
) -> BatchScores {
    // SSSE3 is the preferred 128-bit dispatch target used by upstream. The
    // recurrence itself needs only SSE2 because banded lanes can have distinct
    // query rows and therefore cannot share a pshufb score-table lookup.
    score_batch_sse2_impl::<SEMI_GLOBAL>(query, targets, matrix, cbs, open, extend, scratch)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn score_full_ssse3_impl<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &[i16; 32 * 32],
    cbs: &[i8],
    open: i16,
    extend: i16,
    scratch: &mut PortableSimdScoreScratch,
) -> BatchScores {
    score_full_sse2_impl::<SEMI_GLOBAL>(query, targets, matrix, cbs, open, extend, scratch)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn score_batch_sse2_impl<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &[i16; 32 * 32],
    cbs: &[i8],
    open: i16,
    extend: i16,
    scratch: &mut PortableSimdScoreScratch,
) -> BatchScores {
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    if band == 0 {
        return empty_result(targets.len());
    }

    let delta = if SEMI_GLOBAL { 0 } else { i16::MIN };
    let zero = arch::_mm_set1_epi16(delta);
    let bits_zero = arch::_mm_setzero_si128();
    prepare_rows(scratch, band + 1, zero, zero);
    let open_v = arch::_mm_set1_epi16(open);
    let extend_v = arch::_mm_set1_epi16(extend);
    let max_v = arch::_mm_set1_epi16(i16::MAX);
    let mut best = zero;
    let mut overflow = bits_zero;
    let max_len = targets
        .iter()
        .map(|target| target.subject.len())
        .max()
        .unwrap_or(0);

    for j in 0..max_len {
        scratch.curr_h[band] = zero;
        scratch.curr_e[band] = zero;
        let mut vertical = zero;
        for row in 0..band {
            let mut substitutions = [0i16; 8];
            let mut valid = [0i16; 8];
            let mut scalar_overflow = 0u8;
            for lane in 0..targets.len() {
                let target = targets[lane];
                let width = (target.d_end - target.d_begin).max(0) as usize;
                let qpos = target.d_begin + j as i32 + row as i32;
                if j >= target.subject.len()
                    || row >= width
                    || qpos < 0
                    || qpos >= query.len() as i32
                {
                    continue;
                }
                valid[lane] = -1;
                let qpos = qpos as usize;
                let value = substitution(query[qpos], target.subject[j], matrix)
                    + cbs.get(qpos).copied().unwrap_or(0) as i32;
                if !(i16::MIN as i32..=i16::MAX as i32).contains(&value) {
                    scalar_overflow |= 1 << lane;
                }
                substitutions[lane] = value.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
            }

            let mask = arch::_mm_loadu_si128(valid.as_ptr().cast());
            let subst = arch::_mm_loadu_si128(substitutions.as_ptr().cast());
            let diag = arch::_mm_adds_epi16(scratch.prev_h[row], subst);
            let horizontal = scratch.prev_e[row + 1];
            let mut h = arch::_mm_max_epi16(diag, horizontal);
            h = arch::_mm_max_epi16(h, vertical);
            h = select_x86(mask, h, zero);
            overflow = arch::_mm_or_si128(overflow, arch::_mm_cmpeq_epi16(h, max_v));
            overflow = mark_x86_overflow(overflow, scalar_overflow);

            let opened = arch::_mm_subs_epi16(h, open_v);
            let e = arch::_mm_max_epi16(arch::_mm_subs_epi16(horizontal, extend_v), opened);
            vertical = arch::_mm_max_epi16(arch::_mm_subs_epi16(vertical, extend_v), opened);
            scratch.curr_h[row] = h;
            scratch.curr_e[row] = select_x86(mask, e, zero);
            vertical = select_x86(mask, vertical, zero);
            best = arch::_mm_max_epi16(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }
    finish_x86(best, overflow, targets.len(), delta)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn score_full_sse2_impl<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &[i16; 32 * 32],
    cbs: &[i8],
    open: i16,
    extend: i16,
    scratch: &mut PortableSimdScoreScratch,
) -> BatchScores {
    let rows = query.len() + 1;
    let delta = if SEMI_GLOBAL { 0 } else { i16::MIN };
    let zero = arch::_mm_set1_epi16(delta);
    let bits_zero = arch::_mm_setzero_si128();
    prepare_rows(scratch, rows, zero, zero);
    let open_v = arch::_mm_set1_epi16(open);
    let extend_v = arch::_mm_set1_epi16(extend);
    let max_v = arch::_mm_set1_epi16(i16::MAX);
    let mut best = zero;
    let mut overflow = bits_zero;
    let max_len = targets.iter().map(|target| target.len()).max().unwrap_or(0);

    for j in 0..max_len {
        scratch.curr_h[0] = zero;
        scratch.curr_e[0] = zero;
        let mut vertical = zero;
        for (qpos, &ql) in query.iter().enumerate() {
            let mut substitutions = [0i16; 8];
            let mut valid = [0i16; 8];
            let mut scalar_overflow = 0u8;
            for lane in 0..targets.len() {
                if j >= targets[lane].len() {
                    continue;
                }
                valid[lane] = -1;
                let value = substitution(ql, targets[lane][j], matrix)
                    + cbs.get(qpos).copied().unwrap_or(0) as i32;
                if !(i16::MIN as i32..=i16::MAX as i32).contains(&value) {
                    scalar_overflow |= 1 << lane;
                }
                substitutions[lane] = value.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
            }
            let mask = arch::_mm_loadu_si128(valid.as_ptr().cast());
            let subst = arch::_mm_loadu_si128(substitutions.as_ptr().cast());
            let diag = arch::_mm_adds_epi16(scratch.prev_h[qpos], subst);
            let horizontal = scratch.prev_e[qpos + 1];
            let mut h = arch::_mm_max_epi16(diag, horizontal);
            h = arch::_mm_max_epi16(h, vertical);
            h = select_x86(mask, h, zero);
            overflow = arch::_mm_or_si128(overflow, arch::_mm_cmpeq_epi16(h, max_v));
            overflow = mark_x86_overflow(overflow, scalar_overflow);
            let opened = arch::_mm_subs_epi16(h, open_v);
            let e = arch::_mm_max_epi16(arch::_mm_subs_epi16(horizontal, extend_v), opened);
            vertical = arch::_mm_max_epi16(arch::_mm_subs_epi16(vertical, extend_v), opened);
            scratch.curr_h[qpos + 1] = h;
            scratch.curr_e[qpos + 1] = select_x86(mask, e, zero);
            vertical = select_x86(mask, vertical, zero);
            best = arch::_mm_max_epi16(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }
    finish_x86(best, overflow, targets.len(), delta)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn prepare_rows(
    scratch: &mut PortableSimdScoreScratch,
    rows: usize,
    zero: ArchVector,
    neg: ArchVector,
) {
    scratch.prev_h.resize(rows, zero);
    scratch.curr_h.resize(rows, zero);
    scratch.prev_e.resize(rows, neg);
    scratch.curr_e.resize(rows, neg);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(neg);
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn select_x86(mask: ArchVector, yes: ArchVector, no: ArchVector) -> ArchVector {
    arch::_mm_or_si128(
        arch::_mm_and_si128(mask, yes),
        arch::_mm_andnot_si128(mask, no),
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn mark_x86_overflow(mut overflow: ArchVector, mask: u8) -> ArchVector {
    if mask != 0 {
        let lanes = [
            if mask & 1 != 0 { -1i16 } else { 0 },
            if mask & 2 != 0 { -1i16 } else { 0 },
            if mask & 4 != 0 { -1i16 } else { 0 },
            if mask & 8 != 0 { -1i16 } else { 0 },
            if mask & 16 != 0 { -1i16 } else { 0 },
            if mask & 32 != 0 { -1i16 } else { 0 },
            if mask & 64 != 0 { -1i16 } else { 0 },
            if mask & 128 != 0 { -1i16 } else { 0 },
        ];
        overflow = arch::_mm_or_si128(overflow, arch::_mm_loadu_si128(lanes.as_ptr().cast()));
    }
    overflow
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn finish_x86(
    best: ArchVector,
    overflow: ArchVector,
    len: usize,
    delta: i16,
) -> BatchScores {
    let mut raw_scores = [0i16; 8];
    let mut raw_overflow = [0i16; 8];
    arch::_mm_storeu_si128(raw_scores.as_mut_ptr().cast(), best);
    arch::_mm_storeu_si128(raw_overflow.as_mut_ptr().cast(), overflow);
    finish_result(raw_scores, raw_overflow, len, delta)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn score_batch_neon_impl<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &[i16; 32 * 32],
    cbs: &[i8],
    open: i16,
    extend: i16,
    scratch: &mut PortableSimdScoreScratch,
) -> BatchScores {
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    if band == 0 {
        return empty_result(targets.len());
    }
    let delta = if SEMI_GLOBAL { 0 } else { i16::MIN };
    let zero = arch::vdupq_n_s16(delta);
    prepare_rows_neon(scratch, band + 1, zero, zero);
    let open_v = arch::vdupq_n_s16(open);
    let extend_v = arch::vdupq_n_s16(extend);
    let max_v = arch::vdupq_n_s16(i16::MAX);
    let mut best = zero;
    let mut overflow = arch::vdupq_n_u16(0);
    let max_len = targets
        .iter()
        .map(|target| target.subject.len())
        .max()
        .unwrap_or(0);

    for j in 0..max_len {
        scratch.curr_h[band] = zero;
        scratch.curr_e[band] = zero;
        let mut vertical = zero;
        for row in 0..band {
            let mut substitutions = [0i16; 8];
            let mut valid = [0i16; 8];
            let mut overflow_lanes = [0u16; 8];
            for lane in 0..targets.len() {
                let target = targets[lane];
                let width = (target.d_end - target.d_begin).max(0) as usize;
                let qpos = target.d_begin + j as i32 + row as i32;
                if j >= target.subject.len()
                    || row >= width
                    || qpos < 0
                    || qpos >= query.len() as i32
                {
                    continue;
                }
                valid[lane] = -1;
                let qpos = qpos as usize;
                let value = substitution(query[qpos], target.subject[j], matrix)
                    + cbs.get(qpos).copied().unwrap_or(0) as i32;
                if !(i16::MIN as i32..=i16::MAX as i32).contains(&value) {
                    overflow_lanes[lane] = u16::MAX;
                }
                substitutions[lane] = value.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
            }
            let mask = arch::vld1q_s16(valid.as_ptr());
            let subst = arch::vld1q_s16(substitutions.as_ptr());
            let diag = arch::vqaddq_s16(scratch.prev_h[row], subst);
            let horizontal = scratch.prev_e[row + 1];
            let mut h = arch::vmaxq_s16(diag, horizontal);
            h = arch::vmaxq_s16(h, vertical);
            h = select_neon(mask, h, zero);
            overflow = arch::vorrq_u16(overflow, arch::vceqq_s16(h, max_v));
            overflow = arch::vorrq_u16(overflow, arch::vld1q_u16(overflow_lanes.as_ptr()));
            let opened = arch::vqsubq_s16(h, open_v);
            let e = arch::vmaxq_s16(arch::vqsubq_s16(horizontal, extend_v), opened);
            vertical = arch::vmaxq_s16(arch::vqsubq_s16(vertical, extend_v), opened);
            scratch.curr_h[row] = h;
            scratch.curr_e[row] = select_neon(mask, e, zero);
            vertical = select_neon(mask, vertical, zero);
            best = arch::vmaxq_s16(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }
    finish_neon(best, overflow, targets.len(), delta)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn score_full_neon_impl<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &[i16; 32 * 32],
    cbs: &[i8],
    open: i16,
    extend: i16,
    scratch: &mut PortableSimdScoreScratch,
) -> BatchScores {
    let rows = query.len() + 1;
    let delta = if SEMI_GLOBAL { 0 } else { i16::MIN };
    let zero = arch::vdupq_n_s16(delta);
    prepare_rows_neon(scratch, rows, zero, zero);
    let open_v = arch::vdupq_n_s16(open);
    let extend_v = arch::vdupq_n_s16(extend);
    let max_v = arch::vdupq_n_s16(i16::MAX);
    let mut best = zero;
    let mut overflow = arch::vdupq_n_u16(0);
    let max_len = targets.iter().map(|target| target.len()).max().unwrap_or(0);

    for j in 0..max_len {
        scratch.curr_h[0] = zero;
        scratch.curr_e[0] = zero;
        let mut vertical = zero;
        for (qpos, &ql) in query.iter().enumerate() {
            let mut substitutions = [0i16; 8];
            let mut valid = [0i16; 8];
            let mut overflow_lanes = [0u16; 8];
            for lane in 0..targets.len() {
                if j >= targets[lane].len() {
                    continue;
                }
                valid[lane] = -1;
                let value = substitution(ql, targets[lane][j], matrix)
                    + cbs.get(qpos).copied().unwrap_or(0) as i32;
                if !(i16::MIN as i32..=i16::MAX as i32).contains(&value) {
                    overflow_lanes[lane] = u16::MAX;
                }
                substitutions[lane] = value.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
            }
            let mask = arch::vld1q_s16(valid.as_ptr());
            let subst = arch::vld1q_s16(substitutions.as_ptr());
            let diag = arch::vqaddq_s16(scratch.prev_h[qpos], subst);
            let horizontal = scratch.prev_e[qpos + 1];
            let mut h = arch::vmaxq_s16(diag, horizontal);
            h = arch::vmaxq_s16(h, vertical);
            h = select_neon(mask, h, zero);
            overflow = arch::vorrq_u16(overflow, arch::vceqq_s16(h, max_v));
            overflow = arch::vorrq_u16(overflow, arch::vld1q_u16(overflow_lanes.as_ptr()));
            let opened = arch::vqsubq_s16(h, open_v);
            let e = arch::vmaxq_s16(arch::vqsubq_s16(horizontal, extend_v), opened);
            vertical = arch::vmaxq_s16(arch::vqsubq_s16(vertical, extend_v), opened);
            scratch.curr_h[qpos + 1] = h;
            scratch.curr_e[qpos + 1] = select_neon(mask, e, zero);
            vertical = select_neon(mask, vertical, zero);
            best = arch::vmaxq_s16(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }
    finish_neon(best, overflow, targets.len(), delta)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn prepare_rows_neon(
    scratch: &mut PortableSimdScoreScratch,
    rows: usize,
    zero: ArchVector,
    neg: ArchVector,
) {
    scratch.prev_h.resize(rows, zero);
    scratch.curr_h.resize(rows, zero);
    scratch.prev_e.resize(rows, neg);
    scratch.curr_e.resize(rows, neg);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(neg);
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn select_neon(mask: ArchVector, yes: ArchVector, no: ArchVector) -> ArchVector {
    arch::vbslq_s16(arch::vreinterpretq_u16_s16(mask), yes, no)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn finish_neon(
    best: ArchVector,
    overflow: arch::uint16x8_t,
    len: usize,
    delta: i16,
) -> BatchScores {
    let mut raw_scores = [0i16; 8];
    let mut raw_overflow = [0u16; 8];
    arch::vst1q_s16(raw_scores.as_mut_ptr(), best);
    arch::vst1q_u16(raw_overflow.as_mut_ptr(), overflow);
    let raw_overflow = raw_overflow.map(|lane| lane as i16);
    finish_result(raw_scores, raw_overflow, len, delta)
}

#[inline]
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
fn substitution(query: Letter, subject: Letter, matrix: &[i16; 32 * 32]) -> i32 {
    matrix[((query & LETTER_MASK) as usize) * 32 + (subject & LETTER_MASK) as usize] as i32
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
fn empty_result(len: usize) -> BatchScores {
    BatchScores {
        scores: [0; 16],
        overflow_mask: 0,
        len,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
fn finish_result(
    raw_scores: [i16; 8],
    raw_overflow: [i16; 8],
    len: usize,
    delta: i16,
) -> BatchScores {
    let mut scores = [0i32; 16];
    let mut overflow_mask = 0u16;
    for lane in 0..len {
        scores[lane] = raw_scores[lane] as i32 - delta as i32;
        if raw_overflow[lane] != 0 {
            overflow_mask |= 1 << lane;
        }
    }
    BatchScores {
        scores,
        overflow_mask,
        len,
    }
}

#[cfg(all(
    test,
    any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")
))]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;

    fn scalar(query: &[Letter], target: ScoreTarget<'_>, matrix: &ScoreMatrix, cbs: &[i8]) -> i32 {
        let neg = i32::MIN / 4;
        let open = matrix.gap_open() + matrix.gap_extend();
        let extend = matrix.gap_extend();
        let mut prev_h = vec![0; query.len() + 1];
        let mut prev_e = vec![neg; query.len() + 1];
        let mut best = 0;
        for (j, &sl) in target.subject.iter().enumerate() {
            let mut curr_h = vec![0; query.len() + 1];
            let mut curr_e = vec![neg; query.len() + 1];
            let mut vertical = neg;
            for i in 1..=query.len() {
                let qpos = i - 1;
                if qpos as i32 >= target.d_begin + j as i32
                    && (qpos as i32) < target.d_end + j as i32
                {
                    let subst = substitution(query[qpos], sl, matrix.matrix16())
                        + cbs.get(qpos).copied().unwrap_or(0) as i32;
                    let h = (prev_h[i - 1] + subst).max(prev_e[i]).max(vertical).max(0);
                    curr_h[i] = h;
                    curr_e[i] = (prev_e[i] - extend).max(h - open);
                    vertical = (vertical - extend).max(h - open);
                    best = best.max(h);
                }
            }
            prev_h = curr_h;
            prev_e = curr_e;
        }
        best
    }

    fn random_cases(mut check: impl FnMut(&[Letter], &[ScoreTarget<'_>], &ScoreMatrix, &[i8])) {
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0x128b_17de_5eed_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        for _ in 0..80 {
            let qlen = 1 + next() as usize % 73;
            let query: Vec<Letter> = (0..qlen).map(|_| (next() % 25) as Letter).collect();
            let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
            let count = 1 + next() as usize % 8;
            let mut subjects = Vec::with_capacity(count);
            let mut bands = Vec::with_capacity(count);
            for _ in 0..count {
                let len = 1 + next() as usize % 81;
                let mut subject: Vec<Letter> = (0..len).map(|_| (next() % 25) as Letter).collect();
                if next() % 3 == 0 {
                    let pos = next() as usize % len;
                    subject[pos] |= SEED_MASK;
                }
                subjects.push(subject);
                let begin = next() as i32 % 31 - 15;
                let width = next() as i32 % 41;
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
            check(&query, &targets, &matrix, &cbs);
        }
    }

    fn assert_batch(
        query: &[Letter],
        targets: &[ScoreTarget<'_>],
        matrix: &ScoreMatrix,
        cbs: &[i8],
        got: BatchScores,
    ) {
        assert_eq!(got.len, targets.len());
        assert_eq!(got.overflow_mask, 0);
        for lane in 0..targets.len() {
            assert_eq!(got.scores[lane], scalar(query, targets[lane], matrix, cbs));
        }
    }

    fn assert_full_batch(
        query: &[Letter],
        targets: &[ScoreTarget<'_>],
        matrix: &ScoreMatrix,
        cbs: &[i8],
        got: BatchScores,
    ) {
        assert_eq!(got.len, targets.len());
        assert_eq!(got.overflow_mask, 0);
        for lane in 0..targets.len() {
            let full_target = ScoreTarget {
                subject: targets[lane].subject,
                d_begin: -(targets[lane].subject.len() as i32 - 1),
                d_end: query.len() as i32,
            };
            assert_eq!(got.scores[lane], scalar(query, full_target, matrix, cbs));
        }
    }

    #[test]
    fn randomized_portable_banded_matches_i32_scalar() {
        random_cases(|query, targets, matrix, cbs| {
            let got = score_batch_portable_i16(
                query,
                targets,
                matrix,
                cbs,
                false,
                &mut PortableSimdScoreScratch::default(),
            )
            .expect("128-bit SIMD is available on supported test targets");
            assert_batch(query, targets, matrix, cbs, got);
        });
    }

    #[test]
    fn portable_i16_treats_seed_mask_as_lookup_only_metadata() {
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query = vec![17; 12];
        let plain = vec![17; 12];
        let marked: Vec<_> = plain.iter().map(|&letter| letter | SEED_MASK).collect();
        let targets = [
            ScoreTarget {
                subject: &plain,
                d_begin: 0,
                d_end: 1,
            },
            ScoreTarget {
                subject: &marked,
                d_begin: 0,
                d_end: 1,
            },
        ];
        let mut scratch = PortableSimdScoreScratch::default();
        let banded =
            score_batch_portable_i16(&query, &targets, &matrix, &[], false, &mut scratch).unwrap();
        assert_eq!(banded.scores[0], banded.scores[1]);
        let full = score_full_batch_portable_i16(
            &query,
            &[plain.as_slice(), marked.as_slice()],
            &matrix,
            &[],
            false,
            &mut scratch,
        )
        .unwrap();
        assert_eq!(full.scores[0], full.scores[1]);
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn randomized_explicit_sse2_and_ssse3_match_i32_scalar() {
        random_cases(|query, targets, matrix, cbs| {
            let subjects: Vec<&[Letter]> = targets.iter().map(|target| target.subject).collect();
            if std::arch::is_x86_feature_detected!("sse2") {
                let got = unsafe {
                    score_batch_sse2_impl::<false>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        12,
                        1,
                        &mut PortableSimdScoreScratch::default(),
                    )
                };
                assert_batch(query, targets, matrix, cbs, got);
                let got = unsafe {
                    score_full_sse2_impl::<false>(
                        query,
                        &subjects,
                        matrix.matrix16(),
                        cbs,
                        12,
                        1,
                        &mut PortableSimdScoreScratch::default(),
                    )
                };
                assert_full_batch(query, targets, matrix, cbs, got);
            }
            if std::arch::is_x86_feature_detected!("ssse3") {
                let got = unsafe {
                    score_batch_ssse3_impl::<false>(
                        query,
                        targets,
                        matrix.matrix16(),
                        cbs,
                        12,
                        1,
                        &mut PortableSimdScoreScratch::default(),
                    )
                };
                assert_batch(query, targets, matrix, cbs, got);
                let got = unsafe {
                    score_full_ssse3_impl::<false>(
                        query,
                        &subjects,
                        matrix.matrix16(),
                        cbs,
                        12,
                        1,
                        &mut PortableSimdScoreScratch::default(),
                    )
                };
                assert_full_batch(query, targets, matrix, cbs, got);
            }
        });
    }

    #[test]
    fn randomized_portable_full_matches_i32_scalar() {
        random_cases(|query, targets, matrix, cbs| {
            let subjects: Vec<&[Letter]> = targets.iter().map(|target| target.subject).collect();
            let got = score_full_batch_portable_i16(
                query,
                &subjects,
                matrix,
                cbs,
                false,
                &mut PortableSimdScoreScratch::default(),
            )
            .expect("128-bit SIMD is available on supported test targets");
            assert_full_batch(query, targets, matrix, cbs, got);
        });
    }

    #[test]
    fn reports_i16_overflow_per_lane() {
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query = vec![17; 7_000];
        let long = vec![17; 7_000];
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
        let got = score_batch_portable_i16(
            &query,
            &targets,
            &matrix,
            &[],
            false,
            &mut PortableSimdScoreScratch::default(),
        )
        .unwrap();
        assert_eq!(got.overflow_mask, 1);
        assert_eq!(got.scores[1], 220);
    }

    #[test]
    fn semi_global_does_not_clamp_each_cell_to_zero() {
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query: Vec<Letter> = b"011000010110111010111100000101"
            .iter()
            .map(|&x| (x - b'0') as Letter)
            .collect();
        let subject: Vec<Letter> = b"11010110011001100101010101"
            .iter()
            .map(|&x| (x - b'0') as Letter)
            .collect();
        let targets = [ScoreTarget {
            subject: &subject,
            d_begin: -3,
            d_end: 13,
        }];
        let got = score_batch_portable_i16(
            &query,
            &targets,
            &matrix,
            &[],
            true,
            &mut PortableSimdScoreScratch::default(),
        )
        .unwrap();
        assert_eq!(got.overflow_mask, 0);
        assert_eq!(got.scores[0], 69);
    }

    #[test]
    fn rejects_invalid_contract_inputs() {
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query = [0];
        let target = ScoreTarget {
            subject: &query,
            d_begin: 0,
            d_end: 1,
        };
        let mut scratch = PortableSimdScoreScratch::default();
        assert!(score_batch_portable_i16(&query, &[], &matrix, &[], false, &mut scratch).is_none());
        assert!(
            score_batch_portable_i16(&query, &[target; 8], &matrix, &[], false, &mut scratch)
                .is_some()
        );
        assert!(
            score_batch_portable_i16(&query, &[target; 9], &matrix, &[], false, &mut scratch)
                .is_none()
        );
        assert!(
            score_batch_portable_i16(&query, &[target], &matrix, &[], false, &mut scratch)
                .is_some()
        );
    }
}
