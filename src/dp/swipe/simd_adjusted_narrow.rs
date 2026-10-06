//! Narrow score-only SWIPE tiers for per-target composition-adjusted matrices.
//!
//! Unlike the ordinary kernels, every lane may use a different substitution
//! table and gap scale.  A 32-entry profile is therefore assembled for each
//! packed subject column; the DP loop itself remains entirely vectorized.

use super::simd_trace::TraceTarget;
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::basic::value::LETTER_MASK;
use crate::stats::score_matrix::ScoreMatrix;

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdjustedScores {
    pub scores: Vec<i32>,
    pub overflow_mask: u32,
}

fn valid(query: &[Letter], targets: &[TraceTarget<'_>], cbs: &[i8], lanes: usize) -> bool {
    !targets.is_empty()
        && targets.len() <= lanes
        && (cbs.is_empty() || cbs.len() >= query.len())
        && targets.iter().all(|t| t.d_end > t.d_begin)
}

/// AVX2 32-lane signed-byte tier. Saturated or unrepresentable lanes are
/// marked for promotion to the word tier.
pub fn score_batch_avx2_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
) -> Option<AdjustedScores> {
    if !valid(query, targets, cbs, 32) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 is runtime detected; fixed arrays cover every lane.
        return Some(unsafe {
            if semi_global {
                score_i8::<true>(query, targets, matrix, cbs)
            } else {
                score_i8::<false>(query, targets, matrix, cbs)
            }
        });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    let _ = matrix;
    None
}

/// AVX2 16-lane signed-word tier. Saturated or unrepresentable lanes are
/// marked for promotion to the exact i32 tier.
pub fn score_batch_avx2_i16(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
) -> Option<AdjustedScores> {
    if !valid(query, targets, cbs, 16) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 is runtime detected; fixed arrays cover every lane.
        return Some(unsafe {
            if semi_global {
                score_i16::<true>(query, targets, matrix, cbs)
            } else {
                score_i16::<false>(query, targets, matrix, cbs)
            }
        });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    let _ = matrix;
    None
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn geometry<const LANES: usize>(
    query_len: usize,
    targets: &[TraceTarget<'_>],
) -> (usize, i32, [i32; LANES], [usize; LANES], usize) {
    let band = targets
        .iter()
        .map(|t| (t.d_end - t.d_begin) as usize)
        .max()
        .unwrap_or(0);
    let i1 = targets
        .iter()
        .map(|t| (t.d_end - 1).max(0))
        .min()
        .unwrap_or(0);
    let i0 = i1 + 1 - band as i32;
    let mut subject_start = [0; LANES];
    let mut band_offset = [0; LANES];
    let mut cols = 0;
    for (lane, target) in targets.iter().enumerate() {
        let expanded_begin = target.d_end - band as i32;
        subject_start[lane] = i1 - (target.d_end - 1);
        band_offset[lane] = (target.d_begin - expanded_begin).max(0) as usize;
        let subject_end =
            ((query_len as i32 - 1 - expanded_begin).min(target.subject.len() as i32 - 1) + 1)
                .max(0);
        cols = cols.max((subject_end - subject_start[lane]).max(0) as usize);
    }
    (band, i0, subject_start, band_offset, cols)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_i8<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> AdjustedScores {
    const LANES: usize = 32;
    let (band, i0, starts, offsets, cols) = geometry::<LANES>(query.len(), targets);
    let delta = if SEMI_GLOBAL { 0 } else { i8::MIN };
    let zero = arch::_mm256_set1_epi8(delta);
    let max = arch::_mm256_set1_epi8(i8::MAX);
    let mut ph = vec![zero; band];
    let mut ch = vec![zero; band];
    let mut pe = vec![zero; band + 1];
    let mut ce = vec![zero; band + 1];
    let mut best = zero;
    let mut overflow = 0u32;
    let mut go = [0i8; LANES];
    let mut ge = [0i8; LANES];
    for lane in 0..targets.len() {
        let scale = targets[lane].matrix_scale.max(1);
        let go32 = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge32 = matrix.gap_extend().saturating_mul(scale);
        if !(0..=63).contains(&go32) || !(0..=63).contains(&ge32) {
            overflow |= 1 << lane;
        }
        go[lane] = go32.clamp(0, 63) as i8;
        ge[lane] = ge32.clamp(0, 63) as i8;
    }
    let gov = arch::_mm256_loadu_si256(go.as_ptr().cast());
    let gev = arch::_mm256_loadu_si256(ge.as_ptr().cast());
    for col in 0..cols {
        ch.fill(zero);
        ce.fill(zero);
        let mut vertical = zero;
        let mut live = [false; LANES];
        let mut pos = [0usize; LANES];
        let mut profile = [[0i8; LANES]; 32];
        for lane in 0..targets.len() {
            let p = starts[lane] + col as i32;
            if p < 0 || p >= targets[lane].subject.len() as i32 {
                continue;
            }
            live[lane] = true;
            pos[lane] = p as usize;
            let sl = targets[lane].subject[p as usize];
            for ql in 0..32 {
                if let Some(adjusted) = targets[lane].matrix {
                    profile[ql][lane] = adjusted.scores[(sl & LETTER_MASK) as usize * 32 + ql];
                } else {
                    profile[ql][lane] = matrix.matrix8()[ql * 32 + (sl & LETTER_MASK) as usize];
                }
            }
        }
        for r in 0..band {
            let q = i0 + col as i32 + r as i32;
            if q < 0 || q >= query.len() as i32 {
                vertical = zero;
                continue;
            }
            let qpos = q as usize;
            let mut mask = [0i8; LANES];
            let mut subst = profile[(query[qpos] & LETTER_MASK) as usize];
            for lane in 0..targets.len() {
                if !live[lane] || r < offsets[lane] {
                    continue;
                }
                mask[lane] = -1;
                if targets[lane].matrix.is_none() {
                    let value =
                        i16::from(subst[lane]) + i16::from(cbs.get(qpos).copied().unwrap_or(0));
                    if !(i8::MIN as i16..=i8::MAX as i16).contains(&value) {
                        overflow |= 1 << lane;
                    }
                    subst[lane] = value.clamp(i8::MIN as i16, i8::MAX as i16) as i8;
                }
            }
            let mv = arch::_mm256_loadu_si256(mask.as_ptr().cast());
            let sv = arch::_mm256_loadu_si256(subst.as_ptr().cast());
            let diag = arch::_mm256_adds_epi8(ph[r], sv);
            let horizontal = pe[r + 1];
            let mut h = arch::_mm256_max_epi8(diag, horizontal);
            h = arch::_mm256_max_epi8(h, vertical);
            if !SEMI_GLOBAL {
                h = arch::_mm256_max_epi8(h, zero);
            }
            h = arch::_mm256_blendv_epi8(zero, h, mv);
            overflow |= arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(h, max)) as u32;
            let open = arch::_mm256_subs_epi8(h, gov);
            let nh = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, gev), open);
            let nv = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical, gev), open);
            ch[r] = h;
            ce[r] = arch::_mm256_blendv_epi8(zero, nh, mv);
            vertical = arch::_mm256_blendv_epi8(zero, nv, mv);
            best = arch::_mm256_max_epi8(best, h);
        }
        std::mem::swap(&mut ph, &mut ch);
        std::mem::swap(&mut pe, &mut ce);
    }
    let mut raw = [i8::MIN; LANES];
    arch::_mm256_storeu_si256(raw.as_mut_ptr().cast(), best);
    AdjustedScores {
        scores: raw[..targets.len()]
            .iter()
            .map(|&x| i32::from(x) - i32::from(delta))
            .collect(),
        overflow_mask: overflow & ((1u64 << targets.len()).wrapping_sub(1) as u32),
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_i16<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> AdjustedScores {
    const LANES: usize = 16;
    let (band, i0, starts, offsets, cols) = geometry::<LANES>(query.len(), targets);
    let delta = if SEMI_GLOBAL { 0 } else { i16::MIN };
    let zero = arch::_mm256_set1_epi16(delta);
    let max = arch::_mm256_set1_epi16(i16::MAX);
    let mut ph = vec![zero; band];
    let mut ch = vec![zero; band];
    let mut pe = vec![zero; band + 1];
    let mut ce = vec![zero; band + 1];
    let mut best = zero;
    let mut overflow = 0u32;
    let mut go = [0i16; LANES];
    let mut ge = [0i16; LANES];
    for lane in 0..targets.len() {
        let scale = targets[lane].matrix_scale.max(1);
        let go32 = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge32 = matrix.gap_extend().saturating_mul(scale);
        if !(0..=16_000).contains(&go32) || !(0..=16_000).contains(&ge32) {
            overflow |= 1 << lane;
        }
        go[lane] = go32.clamp(0, 16_000) as i16;
        ge[lane] = ge32.clamp(0, 16_000) as i16;
    }
    let gov = arch::_mm256_loadu_si256(go.as_ptr().cast());
    let gev = arch::_mm256_loadu_si256(ge.as_ptr().cast());
    for col in 0..cols {
        ch.fill(zero);
        ce.fill(zero);
        let mut vertical = zero;
        let mut live = [false; LANES];
        let mut profile = [[0i16; LANES]; 32];
        for lane in 0..targets.len() {
            let p = starts[lane] + col as i32;
            if p < 0 || p >= targets[lane].subject.len() as i32 {
                continue;
            }
            live[lane] = true;
            let sl = targets[lane].subject[p as usize];
            for ql in 0..32 {
                profile[ql][lane] = if let Some(adjusted) = targets[lane].matrix {
                    i16::from(adjusted.scores[(sl & LETTER_MASK) as usize * 32 + ql])
                } else {
                    i16::from(matrix.matrix8()[ql * 32 + (sl & LETTER_MASK) as usize])
                };
            }
        }
        for r in 0..band {
            let q = i0 + col as i32 + r as i32;
            if q < 0 || q >= query.len() as i32 {
                vertical = zero;
                continue;
            }
            let qpos = q as usize;
            let mut mask = [0i16; LANES];
            let mut subst = profile[(query[qpos] & LETTER_MASK) as usize];
            for lane in 0..targets.len() {
                if live[lane] && r >= offsets[lane] {
                    mask[lane] = -1;
                    if targets[lane].matrix.is_none() {
                        subst[lane] += i16::from(cbs.get(qpos).copied().unwrap_or(0));
                    }
                }
            }
            let mv = arch::_mm256_loadu_si256(mask.as_ptr().cast());
            let sv = arch::_mm256_loadu_si256(subst.as_ptr().cast());
            let diag = arch::_mm256_adds_epi16(ph[r], sv);
            let horizontal = pe[r + 1];
            let mut h = arch::_mm256_max_epi16(diag, horizontal);
            h = arch::_mm256_max_epi16(h, vertical);
            if !SEMI_GLOBAL {
                h = arch::_mm256_max_epi16(h, zero);
            }
            h = arch::_mm256_blendv_epi8(zero, h, mv);
            let saturated = arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(h, max)) as u32;
            for lane in 0..targets.len() {
                if saturated & (1 << (2 * lane)) != 0 {
                    overflow |= 1 << lane;
                }
            }
            let open = arch::_mm256_subs_epi16(h, gov);
            let nh = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(horizontal, gev), open);
            let nv = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(vertical, gev), open);
            ch[r] = h;
            ce[r] = arch::_mm256_blendv_epi8(zero, nh, mv);
            vertical = arch::_mm256_blendv_epi8(zero, nv, mv);
            best = arch::_mm256_max_epi16(best, h);
        }
        std::mem::swap(&mut ph, &mut ch);
        std::mem::swap(&mut pe, &mut ce);
    }
    let mut raw = [i16::MIN; LANES];
    arch::_mm256_storeu_si256(raw.as_mut_ptr().cast(), best);
    AdjustedScores {
        scores: raw[..targets.len()]
            .iter()
            .map(|&x| i32::from(x) - i32::from(delta))
            .collect(),
        overflow_mask: overflow,
    }
}

#[cfg(all(test, any(target_arch = "x86", target_arch = "x86_64")))]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;
    use crate::stats::cbs::TargetMatrix;

    #[test]
    fn adjusted_narrow_tiers_treat_seed_mask_as_lookup_only_metadata() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let adjusted = TargetMatrix::new(
            (0..1024)
                .map(|index| if index / 32 == index % 32 { 5 } else { -3 })
                .collect(),
            -3,
            5,
        );
        let query = vec![7; 12];
        let plain = vec![7; 12];
        let marked: Vec<_> = plain.iter().map(|&letter| letter | SEED_MASK).collect();
        let targets = [
            TraceTarget {
                subject: &plain,
                d_begin: 0,
                d_end: 1,
                matrix: Some(&adjusted),
                matrix_scale: 1,
            },
            TraceTarget {
                subject: &marked,
                d_begin: 0,
                d_end: 1,
                matrix: Some(&adjusted),
                matrix_scale: 1,
            },
        ];
        let byte = score_batch_avx2_i8(&query, &targets, &matrix, &[], false).unwrap();
        assert_eq!(byte.overflow_mask, 0);
        assert_eq!(byte.scores[0], byte.scores[1]);
        let word = score_batch_avx2_i16(&query, &targets, &matrix, &[], false).unwrap();
        assert_eq!(word.overflow_mask, 0);
        assert_eq!(word.scores[0], word.scores[1]);
    }
}
