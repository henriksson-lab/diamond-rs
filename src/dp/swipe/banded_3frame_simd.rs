//! AVX2 score-only banded three-frame SWIPE.
//!
//! This is the target-parallel `i16` specialization of DIAMOND's
//! `banded_3frame_swipe.cpp`: up to sixteen subjects occupy the AVX2 lanes,
//! while the three translated query frames are visited in interleaved order.
//! Traceback intentionally remains on the scalar `i32` path, as it does in
//! the original implementation.

use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
use crate::basic::value::LETTER_MASK;
use crate::dp::swipe::DpTarget;
use crate::stats::score_matrix::ScoreMatrix;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LaneScore {
    pub score: i32,
    pub query_end: i32,
    pub subject_end: i32,
    pub frame_end: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BatchScores {
    pub lanes: [LaneScore; 16],
    pub overflow_mask: u16,
    pub len: usize,
}

impl Default for BatchScores {
    fn default() -> Self {
        Self {
            lanes: [LaneScore::default(); 16],
            overflow_mask: 0,
            len: 0,
        }
    }
}

#[derive(Default)]
pub struct Scratch {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_h256: Vec<Vector256>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_h256: Vec<Vector256>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_e256: Vec<Vector256>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_e256: Vec<Vector256>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_h128: Vec<Vector128>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_h128: Vec<Vector128>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_e128: Vec<Vector128>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_e128: Vec<Vector128>,
    #[cfg(target_arch = "aarch64")]
    prev_h_neon: Vec<NeonVector>,
    #[cfg(target_arch = "aarch64")]
    curr_h_neon: Vec<NeonVector>,
    #[cfg(target_arch = "aarch64")]
    prev_e_neon: Vec<NeonVector>,
    #[cfg(target_arch = "aarch64")]
    curr_e_neon: Vec<NeonVector>,
}

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type Vector256 = arch::__m256i;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type Vector128 = arch::__m128i;
#[cfg(target_arch = "aarch64")]
use std::arch::aarch64 as neon;
#[cfg(target_arch = "aarch64")]
type NeonVector = neon::int16x8_t;

/// Score at most sixteen targets in parallel.
///
/// `None` selects the architecture-neutral scalar fallback. A set bit in
/// `overflow_mask` means that lane reached `i16::MAX` and must be recomputed
/// with the scalar `i32` kernel.
pub fn score_batch_avx2(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    score_matrix: &ScoreMatrix,
    scratch: &mut Scratch,
) -> Option<BatchScores> {
    if targets.is_empty() || targets.len() > 16 {
        return None;
    }
    let gap_open = score_matrix
        .gap_open()
        .checked_add(score_matrix.gap_extend())?;
    let gap_extend = score_matrix.gap_extend();
    let frame_shift = score_matrix.frame_shift();
    if !(0..=16_000).contains(&gap_open)
        || !(0..=16_000).contains(&gap_extend)
        || !(0..=16_000).contains(&frame_shift)
    {
        return None;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        // SAFETY: runtime AVX2 detection above guards the implementation.
        return Some(unsafe {
            score_batch_avx2_impl(
                query,
                targets,
                score_matrix.matrix16(),
                gap_open as i16,
                gap_extend as i16,
                frame_shift as i16,
                scratch,
            )
        });
    }

    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, score_matrix, scratch);
        None
    }
}

/// Best native implementation for the current CPU.
pub fn score_batch_simd(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    score_matrix: &ScoreMatrix,
    scratch: &mut Scratch,
) -> Option<BatchScores> {
    if let Some(scores) = score_batch_avx2(query, targets, score_matrix, scratch) {
        return Some(scores);
    }
    if let Some(scores) = score_batch_sse2(query, targets, score_matrix, scratch) {
        return Some(scores);
    }
    score_batch_neon(query, targets, score_matrix, scratch)
}

/// Number of target lanes accepted by [`score_batch_simd`] on this machine.
pub fn native_lane_count() -> usize {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            return 16;
        }
        if std::arch::is_x86_feature_detected!("sse2") {
            return 8;
        }
        return 1;
    }
    #[cfg(target_arch = "aarch64")]
    {
        8
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
    {
        1
    }
}

pub fn score_batch_sse2(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    score_matrix: &ScoreMatrix,
    scratch: &mut Scratch,
) -> Option<BatchScores> {
    if targets.is_empty() || targets.len() > 8 {
        return None;
    }
    let gap_open = score_matrix
        .gap_open()
        .checked_add(score_matrix.gap_extend())?;
    let gap_extend = score_matrix.gap_extend();
    let frame_shift = score_matrix.frame_shift();
    if !(0..=16_000).contains(&gap_open)
        || !(0..=16_000).contains(&gap_extend)
        || !(0..=16_000).contains(&frame_shift)
    {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("sse2") {
            return None;
        }
        // SAFETY: runtime SSE2 detection above guards the implementation.
        return Some(unsafe {
            score_batch_sse2_impl(
                query,
                targets,
                score_matrix.matrix16(),
                gap_open as i16,
                gap_extend as i16,
                frame_shift as i16,
                scratch,
            )
        });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, score_matrix, scratch);
        None
    }
}

/// Eight-lane AArch64 NEON implementation corresponding to the upstream
/// `ScoreVector<int16_t, SHRT_MIN>` specialization.
pub fn score_batch_neon(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    score_matrix: &ScoreMatrix,
    scratch: &mut Scratch,
) -> Option<BatchScores> {
    if targets.is_empty() || targets.len() > 8 {
        return None;
    }
    let gap_open = score_matrix
        .gap_open()
        .checked_add(score_matrix.gap_extend())?;
    let gap_extend = score_matrix.gap_extend();
    let frame_shift = score_matrix.frame_shift();
    if !(0..=16_000).contains(&gap_open)
        || !(0..=16_000).contains(&gap_extend)
        || !(0..=16_000).contains(&frame_shift)
    {
        return None;
    }
    #[cfg(target_arch = "aarch64")]
    {
        // NEON is mandatory in the AArch64 architecture.
        return Some(unsafe {
            score_batch_neon_impl(
                query,
                targets,
                score_matrix.matrix16(),
                gap_open as i16,
                gap_extend as i16,
                frame_shift as i16,
                scratch,
            )
        });
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        let _ = (query, targets, score_matrix, scratch);
        None
    }
}

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn blend_neon(a: NeonVector, b: NeonVector, mask: neon::uint16x8_t) -> NeonVector {
    neon::vbslq_s16(mask, b, a)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn score_batch_neon_impl(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    matrix: &[i16; 32 * 32],
    gap_open: i16,
    gap_extend: i16,
    frame_shift: i16,
    scratch: &mut Scratch,
) -> BatchScores {
    const NEG: i16 = -16_384;
    let qlen = query[0].len();
    let shared_band = targets.iter().map(DpTarget::band).max().unwrap_or(0);
    let cells = (qlen + 1) * 3;
    let zero = neon::vdupq_n_s16(0);
    let neg = neon::vdupq_n_s16(NEG);
    scratch.prev_h_neon.resize(cells, zero);
    scratch.curr_h_neon.resize(cells, zero);
    scratch.prev_e_neon.resize(cells, neg);
    scratch.curr_e_neon.resize(cells, neg);
    scratch.prev_h_neon.fill(zero);
    scratch.prev_e_neon.fill(neg);
    let go = neon::vdupq_n_s16(gap_open);
    let ge = neon::vdupq_n_s16(gap_extend);
    let fs = neon::vdupq_n_s16(frame_shift);
    let max_i16 = neon::vdupq_n_s16(i16::MAX);
    let mut best = zero;
    let mut overflow = neon::vdupq_n_u16(0);
    let mut result = BatchScores {
        len: targets.len(),
        ..Default::default()
    };
    let max_subject_len = targets.iter().map(|t| t.seq.len()).max().unwrap_or(0);

    for subject_pos in 0..max_subject_len {
        scratch.curr_h_neon.fill(zero);
        scratch.curr_e_neon.fill(neg);
        let mut vgap = [neg; 3];
        for i in 1..=qlen {
            for frame in 0..3 {
                let cell = i * 3 + frame;
                let mut substitutions = [0i16; 8];
                let mut validity = [0i16; 8];
                for lane in 0..targets.len() {
                    let target = &targets[lane];
                    let qpos = i as i32 - 1;
                    let target_pos = subject_pos as i32;
                    if subject_pos < target.seq.len()
                        && i <= query[frame].len()
                        && qpos >= target.d_end - shared_band + target_pos
                        && qpos < target.d_end + target_pos
                    {
                        validity[lane] = -1;
                        let q = (query[frame][i - 1] & LETTER_MASK) as usize;
                        let s = (target.seq[subject_pos] & LETTER_MASK) as usize;
                        substitutions[lane] = matrix[q * 32 + s];
                    }
                }
                let valid = neon::vreinterpretq_u16_s16(neon::vld1q_s16(validity.as_ptr()));
                let subst = neon::vld1q_s16(substitutions.as_ptr());
                let diag = scratch.prev_h_neon[(i - 1) * 3 + frame];
                let fwd = if frame == 0 {
                    if i >= 2 {
                        scratch.prev_h_neon[(i - 2) * 3 + 2]
                    } else {
                        zero
                    }
                } else {
                    scratch.prev_h_neon[(i - 1) * 3 + frame - 1]
                };
                let rev = if frame == 2 {
                    scratch.prev_h_neon[i * 3]
                } else {
                    scratch.prev_h_neon[(i - 1) * 3 + frame + 1]
                };
                let mut current = neon::vqaddq_s16(diag, subst);
                let shifted = neon::vqsubq_s16(subst, fs);
                current = neon::vmaxq_s16(current, neon::vqaddq_s16(fwd, shifted));
                current = neon::vmaxq_s16(current, neon::vqaddq_s16(rev, shifted));
                current = neon::vmaxq_s16(current, scratch.prev_e_neon[cell]);
                current = neon::vmaxq_s16(current, vgap[frame]);
                current = neon::vmaxq_s16(current, zero);
                current = blend_neon(zero, current, valid);
                let mut hgap = neon::vqsubq_s16(scratch.prev_e_neon[cell], ge);
                let mut next_vgap = neon::vqsubq_s16(vgap[frame], ge);
                let opened = neon::vqsubq_s16(current, go);
                hgap = neon::vmaxq_s16(hgap, opened);
                next_vgap = neon::vmaxq_s16(next_vgap, opened);
                scratch.curr_e_neon[cell] = blend_neon(neg, hgap, valid);
                vgap[frame] = blend_neon(neg, next_vgap, valid);
                scratch.curr_h_neon[cell] = current;
                let improves = neon::vcgtq_s16(current, best);
                let mut improved_lanes = [0u16; 8];
                neon::vst1q_u16(improved_lanes.as_mut_ptr(), improves);
                if improved_lanes.iter().any(|&lane| lane != 0) {
                    for lane in 0..targets.len() {
                        if improved_lanes[lane] != 0 {
                            result.lanes[lane].query_end = i as i32;
                            result.lanes[lane].subject_end = subject_pos as i32 + 1;
                            result.lanes[lane].frame_end = frame;
                        }
                    }
                    best = neon::vmaxq_s16(best, current);
                }
                overflow = neon::vorrq_u16(overflow, neon::vceqq_s16(current, max_i16));
            }
        }
        std::mem::swap(&mut scratch.prev_h_neon, &mut scratch.curr_h_neon);
        std::mem::swap(&mut scratch.prev_e_neon, &mut scratch.curr_e_neon);
    }
    let mut best_lanes = [0i16; 8];
    let mut overflow_lanes = [0u16; 8];
    neon::vst1q_s16(best_lanes.as_mut_ptr(), best);
    neon::vst1q_u16(overflow_lanes.as_mut_ptr(), overflow);
    for lane in 0..targets.len() {
        result.lanes[lane].score = best_lanes[lane] as i32;
        if overflow_lanes[lane] != 0 {
            result.overflow_mask |= 1 << lane;
        }
    }
    result
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
unsafe fn blend128(a: Vector128, b: Vector128, mask: Vector128) -> Vector128 {
    arch::_mm_or_si128(
        arch::_mm_and_si128(mask, b),
        arch::_mm_andnot_si128(mask, a),
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn score_batch_sse2_impl(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    matrix: &[i16; 32 * 32],
    gap_open: i16,
    gap_extend: i16,
    frame_shift: i16,
    scratch: &mut Scratch,
) -> BatchScores {
    const NEG: i16 = -16_384;
    let qlen = query[0].len();
    let shared_band = targets.iter().map(DpTarget::band).max().unwrap_or(0);
    let cells = (qlen + 1) * 3;
    let zero = arch::_mm_setzero_si128();
    let neg = arch::_mm_set1_epi16(NEG);
    scratch.prev_h128.resize(cells, zero);
    scratch.curr_h128.resize(cells, zero);
    scratch.prev_e128.resize(cells, neg);
    scratch.curr_e128.resize(cells, neg);
    scratch.prev_h128.fill(zero);
    scratch.prev_e128.fill(neg);
    let go = arch::_mm_set1_epi16(gap_open);
    let ge = arch::_mm_set1_epi16(gap_extend);
    let fs = arch::_mm_set1_epi16(frame_shift);
    let max_i16 = arch::_mm_set1_epi16(i16::MAX);
    let mut best = zero;
    let mut overflow = zero;
    let mut result = BatchScores {
        len: targets.len(),
        ..Default::default()
    };
    let max_subject_len = targets.iter().map(|t| t.seq.len()).max().unwrap_or(0);

    for subject_pos in 0..max_subject_len {
        scratch.curr_h128.fill(zero);
        scratch.curr_e128.fill(neg);
        let mut vgap = [neg; 3];
        for i in 1..=qlen {
            for frame in 0..3 {
                let cell = i * 3 + frame;
                let mut substitutions = [0i16; 8];
                let mut validity = [0i16; 8];
                for lane in 0..targets.len() {
                    let target = &targets[lane];
                    let qpos = i as i32 - 1;
                    let target_pos = subject_pos as i32;
                    if subject_pos < target.seq.len()
                        && i <= query[frame].len()
                        && qpos >= target.d_end - shared_band + target_pos
                        && qpos < target.d_end + target_pos
                    {
                        validity[lane] = -1;
                        let q = (query[frame][i - 1] & LETTER_MASK) as usize;
                        let s = (target.seq[subject_pos] & LETTER_MASK) as usize;
                        substitutions[lane] = matrix[q * 32 + s];
                    }
                }
                let valid = arch::_mm_loadu_si128(validity.as_ptr().cast());
                let subst = arch::_mm_loadu_si128(substitutions.as_ptr().cast());
                let diag = scratch.prev_h128[(i - 1) * 3 + frame];
                let fwd = if frame == 0 {
                    if i >= 2 {
                        scratch.prev_h128[(i - 2) * 3 + 2]
                    } else {
                        zero
                    }
                } else {
                    scratch.prev_h128[(i - 1) * 3 + frame - 1]
                };
                let rev = if frame == 2 {
                    scratch.prev_h128[i * 3]
                } else {
                    scratch.prev_h128[(i - 1) * 3 + frame + 1]
                };
                let mut current = arch::_mm_adds_epi16(diag, subst);
                let shifted = arch::_mm_subs_epi16(subst, fs);
                current = arch::_mm_max_epi16(current, arch::_mm_adds_epi16(fwd, shifted));
                current = arch::_mm_max_epi16(current, arch::_mm_adds_epi16(rev, shifted));
                current = arch::_mm_max_epi16(current, scratch.prev_e128[cell]);
                current = arch::_mm_max_epi16(current, vgap[frame]);
                current = arch::_mm_max_epi16(current, zero);
                current = blend128(zero, current, valid);
                let mut hgap = arch::_mm_subs_epi16(scratch.prev_e128[cell], ge);
                let mut next_vgap = arch::_mm_subs_epi16(vgap[frame], ge);
                let opened = arch::_mm_subs_epi16(current, go);
                hgap = arch::_mm_max_epi16(hgap, opened);
                next_vgap = arch::_mm_max_epi16(next_vgap, opened);
                scratch.curr_e128[cell] = blend128(neg, hgap, valid);
                vgap[frame] = blend128(neg, next_vgap, valid);
                scratch.curr_h128[cell] = current;
                let improves = arch::_mm_cmpgt_epi16(current, best);
                if arch::_mm_movemask_epi8(improves) != 0 {
                    let mut improved_lanes = [0i16; 8];
                    arch::_mm_storeu_si128(improved_lanes.as_mut_ptr().cast(), improves);
                    for lane in 0..targets.len() {
                        if improved_lanes[lane] != 0 {
                            result.lanes[lane].query_end = i as i32;
                            result.lanes[lane].subject_end = subject_pos as i32 + 1;
                            result.lanes[lane].frame_end = frame;
                        }
                    }
                    best = arch::_mm_max_epi16(best, current);
                }
                overflow = arch::_mm_or_si128(overflow, arch::_mm_cmpeq_epi16(current, max_i16));
            }
        }
        std::mem::swap(&mut scratch.prev_h128, &mut scratch.curr_h128);
        std::mem::swap(&mut scratch.prev_e128, &mut scratch.curr_e128);
    }
    let mut best_lanes = [0i16; 8];
    let mut overflow_lanes = [0i16; 8];
    arch::_mm_storeu_si128(best_lanes.as_mut_ptr().cast(), best);
    arch::_mm_storeu_si128(overflow_lanes.as_mut_ptr().cast(), overflow);
    for lane in 0..targets.len() {
        result.lanes[lane].score = best_lanes[lane] as i32;
        if overflow_lanes[lane] != 0 {
            result.overflow_mask |= 1 << lane;
        }
    }
    result
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_batch_avx2_impl(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    matrix: &[i16; 32 * 32],
    gap_open: i16,
    gap_extend: i16,
    frame_shift: i16,
    scratch: &mut Scratch,
) -> BatchScores {
    const NEG: i16 = -16_384;
    let qlen = query[0].len();
    let shared_band = targets.iter().map(DpTarget::band).max().unwrap_or(0);
    let cells = (qlen + 1) * 3;
    let zero = arch::_mm256_setzero_si256();
    let neg = arch::_mm256_set1_epi16(NEG);
    scratch.prev_h256.resize(cells, zero);
    scratch.curr_h256.resize(cells, zero);
    scratch.prev_e256.resize(cells, neg);
    scratch.curr_e256.resize(cells, neg);
    scratch.prev_h256.fill(zero);
    scratch.prev_e256.fill(neg);

    let go = arch::_mm256_set1_epi16(gap_open);
    let ge = arch::_mm256_set1_epi16(gap_extend);
    let fs = arch::_mm256_set1_epi16(frame_shift);
    let max_i16 = arch::_mm256_set1_epi16(i16::MAX);
    let mut best = zero;
    let mut overflow = zero;
    let mut result = BatchScores {
        len: targets.len(),
        ..Default::default()
    };
    let max_subject_len = targets.iter().map(|t| t.seq.len()).max().unwrap_or(0);

    for subject_pos in 0..max_subject_len {
        scratch.curr_h256.fill(zero);
        scratch.curr_e256.fill(neg);
        let mut vgap = [neg; 3];

        for i in 1..=qlen {
            for frame in 0..3 {
                let cell = i * 3 + frame;
                let mut substitutions = [0i16; 16];
                let mut validity = [0i16; 16];
                for lane in 0..targets.len() {
                    let target = &targets[lane];
                    let qpos = i as i32 - 1;
                    let target_pos = subject_pos as i32;
                    let valid = subject_pos < target.seq.len()
                        && i <= query[frame].len()
                        && qpos >= target.d_end - shared_band + target_pos
                        && qpos < target.d_end + target_pos;
                    if valid {
                        validity[lane] = -1;
                        let q = (query[frame][i - 1] & LETTER_MASK) as usize;
                        let s = (target.seq[subject_pos] & LETTER_MASK) as usize;
                        substitutions[lane] = matrix[q * 32 + s];
                    }
                }
                let valid = arch::_mm256_loadu_si256(validity.as_ptr().cast());
                let subst = arch::_mm256_loadu_si256(substitutions.as_ptr().cast());

                let diag = scratch.prev_h256[(i - 1) * 3 + frame];
                let fwd = if frame == 0 {
                    if i >= 2 {
                        scratch.prev_h256[(i - 2) * 3 + 2]
                    } else {
                        zero
                    }
                } else {
                    scratch.prev_h256[(i - 1) * 3 + frame - 1]
                };
                let rev = if frame == 2 {
                    scratch.prev_h256[i * 3]
                } else {
                    scratch.prev_h256[(i - 1) * 3 + frame + 1]
                };

                let mut current = arch::_mm256_adds_epi16(diag, subst);
                let shift_score = arch::_mm256_subs_epi16(subst, fs);
                current =
                    arch::_mm256_max_epi16(current, arch::_mm256_adds_epi16(fwd, shift_score));
                current =
                    arch::_mm256_max_epi16(current, arch::_mm256_adds_epi16(rev, shift_score));
                current = arch::_mm256_max_epi16(current, scratch.prev_e256[cell]);
                current = arch::_mm256_max_epi16(current, vgap[frame]);
                current = arch::_mm256_max_epi16(current, zero);
                current = arch::_mm256_blendv_epi8(zero, current, valid);

                let mut hgap = arch::_mm256_subs_epi16(scratch.prev_e256[cell], ge);
                let mut next_vgap = arch::_mm256_subs_epi16(vgap[frame], ge);
                let opened = arch::_mm256_subs_epi16(current, go);
                hgap = arch::_mm256_max_epi16(hgap, opened);
                next_vgap = arch::_mm256_max_epi16(next_vgap, opened);
                scratch.curr_e256[cell] = arch::_mm256_blendv_epi8(neg, hgap, valid);
                vgap[frame] = arch::_mm256_blendv_epi8(neg, next_vgap, valid);
                scratch.curr_h256[cell] = current;

                let improves = arch::_mm256_cmpgt_epi16(current, best);
                if arch::_mm256_movemask_epi8(improves) != 0 {
                    let mut improved_lanes = [0i16; 16];
                    arch::_mm256_storeu_si256(improved_lanes.as_mut_ptr().cast(), improves);
                    for lane in 0..targets.len() {
                        if improved_lanes[lane] != 0 {
                            result.lanes[lane].query_end = i as i32;
                            result.lanes[lane].subject_end = subject_pos as i32 + 1;
                            result.lanes[lane].frame_end = frame;
                        }
                    }
                    best = arch::_mm256_max_epi16(best, current);
                }
                overflow =
                    arch::_mm256_or_si256(overflow, arch::_mm256_cmpeq_epi16(current, max_i16));
            }
        }
        std::mem::swap(&mut scratch.prev_h256, &mut scratch.curr_h256);
        std::mem::swap(&mut scratch.prev_e256, &mut scratch.curr_e256);
    }

    let mut best_lanes = [0i16; 16];
    let mut overflow_lanes = [0i16; 16];
    arch::_mm256_storeu_si256(best_lanes.as_mut_ptr().cast(), best);
    arch::_mm256_storeu_si256(overflow_lanes.as_mut_ptr().cast(), overflow);
    for lane in 0..targets.len() {
        result.lanes[lane].score = best_lanes[lane] as i32;
        if overflow_lanes[lane] != 0 {
            result.overflow_mask |= 1 << lane;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dp::banded_3frame::banded_3frame_swipe_range;
    use crate::dp::swipe::{Anchor, CarryOver};

    fn next(seed: &mut u64) -> u32 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*seed >> 32) as u32
    }

    #[test]
    fn randomized_avx2_scores_and_endpoints_match_scalar() {
        if !cfg!(any(target_arch = "x86", target_arch = "x86_64")) {
            return;
        }
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }

        let matrix = ScoreMatrix::new("blosum62", 11, 1, 7, 1, 0).unwrap();
        let mut seed = 0x3f84_d5b5_b547_0917;
        let mut scratch = Scratch::default();
        for _case in 0..100 {
            let qlen = 3 + (next(&mut seed) % 29) as usize;
            let mut frames = [Vec::new(), Vec::new(), Vec::new()];
            for (frame, sequence) in frames.iter_mut().enumerate() {
                let len = qlen.saturating_sub((frame + (next(&mut seed) as usize & 1)) / 2);
                sequence.extend((0..len).map(|_| (next(&mut seed) % 20) as Letter));
            }
            let lane_count = 1 + (next(&mut seed) % 16) as usize;
            let mut targets = Vec::with_capacity(lane_count);
            for lane in 0..lane_count {
                let len = 1 + (next(&mut seed) % 27) as usize;
                let subject: Vec<Letter> = (0..len)
                    .map(|_| {
                        let letter = (next(&mut seed) % 20) as Letter;
                        if next(&mut seed) & 7 == 0 {
                            letter | crate::basic::value::SEED_MASK
                        } else {
                            letter
                        }
                    })
                    .collect();
                let d_begin = -((next(&mut seed) % 8) as i32);
                let width = 1 + (next(&mut seed) % 15) as i32;
                targets.push(DpTarget::new(
                    subject,
                    len as i32,
                    d_begin,
                    d_begin + width,
                    lane as i64,
                    qlen as i32,
                    CarryOver::default(),
                    Anchor::default(),
                ));
            }
            let query = [&frames[0][..], &frames[1][..], &frames[2][..]];
            let simd = score_batch_avx2(query, &targets, &matrix, &mut scratch).unwrap();
            assert_eq!(simd.overflow_mask, 0);
            let shared_band = targets.iter().map(DpTarget::band).max().unwrap();
            for (lane, target) in targets.iter().enumerate() {
                let scalar = banded_3frame_swipe_range(
                    query,
                    &target.seq,
                    target.d_end - shared_band,
                    target.d_end,
                    &matrix,
                );
                assert_eq!(
                    simd.lanes[lane],
                    LaneScore {
                        score: scalar.sw.score,
                        query_end: scalar.sw.query_end,
                        subject_end: scalar.sw.subject_end,
                        frame_end: scalar.frame_end,
                    },
                    "case {_case}, lane {lane}, band {}..{}",
                    target.d_begin,
                    target.d_end,
                );
            }

            // The original SSE2 dispatcher uses eight i16 lanes. Exercise it
            // explicitly even on AVX2 hosts so it cannot silently rot.
            let sse_targets = &targets[..targets.len().min(8)];
            let sse = score_batch_sse2(query, sse_targets, &matrix, &mut scratch).unwrap();
            assert_eq!(sse.overflow_mask, 0);
            let shared_band = sse_targets.iter().map(DpTarget::band).max().unwrap();
            for (lane, target) in sse_targets.iter().enumerate() {
                let scalar = banded_3frame_swipe_range(
                    query,
                    &target.seq,
                    target.d_end - shared_band,
                    target.d_end,
                    &matrix,
                );
                assert_eq!(sse.lanes[lane].score, scalar.sw.score);
                assert_eq!(sse.lanes[lane].query_end, scalar.sw.query_end);
                assert_eq!(sse.lanes[lane].subject_end, scalar.sw.subject_end);
                assert_eq!(sse.lanes[lane].frame_end, scalar.frame_end);
            }
        }
    }

    #[test]
    fn rejects_more_than_sixteen_lanes() {
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 7, 1, 0).unwrap();
        let query = [vec![0; 3], vec![1; 3], vec![2; 3]];
        let targets = (0..17)
            .map(|lane| {
                DpTarget::new(
                    vec![lane as Letter % 20; 2],
                    2,
                    -1,
                    2,
                    lane,
                    3,
                    CarryOver::default(),
                    Anchor::default(),
                )
            })
            .collect::<Vec<_>>();
        assert!(score_batch_avx2(
            [&query[0], &query[1], &query[2]],
            &targets,
            &matrix,
            &mut Scratch::default(),
        )
        .is_none());
    }
}
