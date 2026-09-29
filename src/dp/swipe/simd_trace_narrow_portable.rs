//! Narrow traceback tiers for 128-bit SIMD targets.
//!
//! SSE4.1/SSSE3 and AArch64 NEON process sixteen signed-byte lanes.  The
//! promoted signed-word tier processes eight lanes and only requires SSE2 on
//! x86 (substitution scores are gathered before the vector recurrence).

use super::simd_trace::TraceTarget;
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
use crate::basic::{packed_transcript::EditOperation, value::LETTER_MASK};
use crate::dp::smith_waterman::SwResult;
use crate::stats::score_matrix::ScoreMatrix;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64 as arch;
#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type V8 = arch::__m128i;
#[cfg(target_arch = "aarch64")]
type V8 = arch::int8x16_t;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type V16 = arch::__m128i;
#[cfg(target_arch = "aarch64")]
type V16 = arch::int16x8_t;

pub struct NarrowTraceBatch {
    pub results: Vec<SwResult>,
    pub overflow_mask: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NarrowScores {
    pub scores: Vec<i32>,
    pub overflow_mask: u32,
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
const ACTIVE: u8 = 1;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
const GAP_V: u8 = 1 << 1;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
const GAP_H: u8 = 1 << 2;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
const OPEN_V: u8 = 1 << 3;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
const OPEN_H: u8 = 1 << 4;

pub fn lane_width(score_bin: usize) -> usize {
    if !available() {
        return 0;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if score_bin == 0
        && !(std::arch::is_x86_feature_detected!("sse4.1")
            && std::arch::is_x86_feature_detected!("ssse3"))
    {
        // The byte tier is unavailable, but an eight-lane chunk still fits
        // the exact SSE2 i32 fallback later in the dispatch chain.
        return 8;
    }
    match score_bin {
        0 => 16,
        1 => 8,
        _ => 8,
    }
}

pub fn available() -> bool {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        return std::arch::is_x86_feature_detected!("sse2");
    }
    #[cfg(target_arch = "aarch64")]
    {
        return std::arch::is_aarch64_feature_detected!("neon");
    }
    #[allow(unreachable_code)]
    false
}

fn valid(query: &[Letter], targets: &[TraceTarget<'_>], cbs: &[i8], lanes: usize) -> bool {
    !targets.is_empty()
        && targets.len() <= lanes
        && (cbs.is_empty() || cbs.len() >= query.len())
        && targets.iter().all(|t| t.d_end > t.d_begin)
}

#[cfg_attr(
    not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")),
    allow(unused_variables)
)]
pub fn trace_batch_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> Option<NarrowTraceBatch> {
    if !valid(query, targets, cbs, 16) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("sse4.1") && std::arch::is_x86_feature_detected!("ssse3")
    {
        return Some(unsafe { trace_i8_core(query, targets, matrix, cbs) });
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("neon") {
        return Some(unsafe { trace_i8_core(query, targets, matrix, cbs) });
    }
    None
}

#[cfg_attr(
    not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")),
    allow(unused_variables)
)]
pub fn trace_batch_i16(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> Option<NarrowTraceBatch> {
    if !valid(query, targets, cbs, 8) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("sse2") {
        return Some(unsafe { trace_i16_core(query, targets, matrix, cbs) });
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("neon") {
        return Some(unsafe { trace_i16_core(query, targets, matrix, cbs) });
    }
    None
}

/// Score-only adjusted/mixed byte tier for SSE4.1+SSSE3 and AArch64 NEON.
#[cfg_attr(
    not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")),
    allow(unused_variables)
)]
pub fn score_batch_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> Option<NarrowScores> {
    if !valid(query, targets, cbs, 16) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("sse4.1") && std::arch::is_x86_feature_detected!("ssse3")
    {
        return Some(unsafe { score_i8_core(query, targets, matrix, cbs) });
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("neon") {
        return Some(unsafe { score_i8_core(query, targets, matrix, cbs) });
    }
    None
}

/// Score-only adjusted/mixed word tier for SSE2 and AArch64 NEON.
#[cfg_attr(
    not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")),
    allow(unused_variables)
)]
pub fn score_batch_i16(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> Option<NarrowScores> {
    if !valid(query, targets, cbs, 8) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("sse2") {
        return Some(unsafe { score_i16_core(query, targets, matrix, cbs) });
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("neon") {
        return Some(unsafe { score_i16_core(query, targets, matrix, cbs) });
    }
    None
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
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

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse4.1,ssse3")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn score_i8_core(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> NarrowScores {
    const LANES: usize = 16;
    let (band, i0, starts, offsets, cols) = geometry::<LANES>(query.len(), targets);
    let zero = v8_splat(i8::MIN);
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
        let xgo = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let xge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=63).contains(&xgo) || !(0..=63).contains(&xge) {
            overflow |= 1 << lane;
        }
        go[lane] = xgo.clamp(0, 63) as i8;
        ge[lane] = xge.clamp(0, 63) as i8;
    }
    let go = v8_load(go);
    let ge = v8_load(ge);
    for col in 0..cols {
        ch.fill(zero);
        ce.fill(zero);
        let mut vertical = zero;
        let mut live = [false; LANES];
        let mut profile = [[0i8; LANES]; 32];
        for lane in 0..targets.len() {
            let p = starts[lane] + col as i32;
            if p < 0 || p >= targets[lane].subject.len() as i32 {
                continue;
            }
            live[lane] = true;
            let sl = targets[lane].subject[p as usize];
            for ql in 0..32 {
                profile[ql][lane] = if let Some(adjusted) = targets[lane].matrix {
                    adjusted.scores[(sl & LETTER_MASK) as usize * 32 + ql]
                } else {
                    matrix.matrix8()[ql * 32 + (sl & LETTER_MASK) as usize]
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
            let mut mask = [0i8; LANES];
            let mut subst = profile[(query[qpos] & LETTER_MASK) as usize];
            for lane in 0..targets.len() {
                if live[lane] && r >= offsets[lane] {
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
            }
            let maskv = v8_load(mask);
            let diag = v8_adds(ph[r], v8_load(subst));
            let horizontal = pe[r + 1];
            let h = v8_select(
                maskv,
                v8_max(v8_max(v8_max(diag, horizontal), vertical), zero),
                zero,
            );
            let open = v8_subs(h, go);
            let nh = v8_max(v8_subs(horizontal, ge), open);
            let nv = v8_max(v8_subs(vertical, ge), open);
            ch[r] = h;
            ce[r] = v8_select(maskv, nh, zero);
            vertical = v8_select(maskv, nv, zero);
            best = v8_max(best, h);
            for (lane, &x) in v8_store(h)[..targets.len()].iter().enumerate() {
                if x == i8::MAX {
                    overflow |= 1 << lane;
                }
            }
        }
        std::mem::swap(&mut ph, &mut ch);
        std::mem::swap(&mut pe, &mut ce);
    }
    NarrowScores {
        scores: v8_store(best)[..targets.len()]
            .iter()
            .map(|&x| i32::from(x) - i32::from(i8::MIN))
            .collect(),
        overflow_mask: overflow,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn score_i16_core(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> NarrowScores {
    const LANES: usize = 8;
    let (band, i0, starts, offsets, cols) = geometry::<LANES>(query.len(), targets);
    let zero = v16_splat(i16::MIN);
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
        let xgo = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let xge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=16_000).contains(&xgo) || !(0..=16_000).contains(&xge) {
            overflow |= 1 << lane;
        }
        go[lane] = xgo.clamp(0, 16_000) as i16;
        ge[lane] = xge.clamp(0, 16_000) as i16;
    }
    let go = v16_load(go);
    let ge = v16_load(ge);
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
            let maskv = v16_load(mask);
            let diag = v16_adds(ph[r], v16_load(subst));
            let horizontal = pe[r + 1];
            let h = v16_select(
                maskv,
                v16_max(v16_max(v16_max(diag, horizontal), vertical), zero),
                zero,
            );
            let open = v16_subs(h, go);
            let nh = v16_max(v16_subs(horizontal, ge), open);
            let nv = v16_max(v16_subs(vertical, ge), open);
            ch[r] = h;
            ce[r] = v16_select(maskv, nh, zero);
            vertical = v16_select(maskv, nv, zero);
            best = v16_max(best, h);
            for (lane, &x) in v16_store(h)[..targets.len()].iter().enumerate() {
                if x == i16::MAX {
                    overflow |= 1 << lane;
                }
            }
        }
        std::mem::swap(&mut ph, &mut ch);
        std::mem::swap(&mut pe, &mut ce);
    }
    NarrowScores {
        scores: v16_store(best)[..targets.len()]
            .iter()
            .map(|&x| i32::from(x) - i32::from(i16::MIN))
            .collect(),
        overflow_mask: overflow,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse4.1,ssse3")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn trace_i8_core(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> NarrowTraceBatch {
    const LANES: usize = 16;
    let (band, i0, subject_start, band_offset, cols) = geometry::<LANES>(query.len(), targets);
    let zero = v8_splat(i8::MIN);
    let mut go_lanes = [0i8; LANES];
    let mut ge_lanes = [0i8; LANES];
    let mut overflow_mask = 0u32;
    for lane in 0..targets.len() {
        let scale = targets[lane].matrix_scale.max(1);
        let go = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=63).contains(&go) || !(0..=63).contains(&ge) {
            overflow_mask |= 1 << lane;
        }
        go_lanes[lane] = go.clamp(0, 63) as i8;
        ge_lanes[lane] = ge.clamp(0, 63) as i8;
    }
    let go = v8_load(go_lanes);
    let ge = v8_load(ge_lanes);
    let mut ph = vec![zero; band];
    let mut ch = vec![zero; band];
    let mut pe = vec![zero; band + 1];
    let mut ce = vec![zero; band + 1];
    let mut trace = allocate_trace(targets);
    let mut best = [i8::MIN; LANES];
    let mut best_i = [0; LANES];
    let mut best_j = [0; LANES];

    for column in 0..cols {
        ch.fill(zero);
        ce.fill(zero);
        let mut vertical = zero;
        let mut live = [false; LANES];
        let mut subject_pos = [0; LANES];
        for lane in 0..targets.len() {
            let pos = subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                live[lane] = true;
                subject_pos[lane] = pos as usize;
            }
        }
        for r in 0..band {
            let q = i0 + column as i32 + r as i32;
            if q < 0 || q >= query.len() as i32 {
                vertical = zero;
                continue;
            }
            let qpos = q as usize;
            let mut mask = [0i8; LANES];
            let mut subst = [0i8; LANES];
            let qletter = (query[qpos] & LETTER_MASK) as usize;
            for lane in 0..targets.len() {
                if !live[lane] || r < band_offset[lane] {
                    continue;
                }
                mask[lane] = -1;
                let sl = targets[lane].subject[subject_pos[lane]];
                let raw = if let Some(adjusted) = targets[lane].matrix {
                    adjusted.scores[(sl & LETTER_MASK) as usize * 32 + qletter]
                } else {
                    matrix.matrix8()[qletter * 32 + (sl & LETTER_MASK) as usize]
                };
                let value = i16::from(raw)
                    + if targets[lane].matrix.is_none() {
                        i16::from(cbs.get(qpos).copied().unwrap_or(0))
                    } else {
                        0
                    };
                if !(i8::MIN as i16..=i8::MAX as i16).contains(&value) {
                    overflow_mask |= 1 << lane;
                }
                subst[lane] = value.clamp(i8::MIN as i16, i8::MAX as i16) as i8;
            }
            let maskv = v8_load(mask);
            let diag = v8_adds(ph[r], v8_load(subst));
            let horizontal = pe[r + 1];
            let score = v8_select(
                maskv,
                v8_max(v8_max(v8_max(diag, horizontal), vertical), zero),
                zero,
            );
            let open = v8_subs(score, go);
            let next_h = v8_max(v8_subs(horizontal, ge), open);
            let next_v = v8_max(v8_subs(vertical, ge), open);
            ch[r] = score;
            ce[r] = v8_select(maskv, next_h, zero);
            let scores = v8_store(score);
            let verticals = v8_store(vertical);
            let horizontals = v8_store(horizontal);
            let opens = v8_store(open);
            let next_vs = v8_store(next_v);
            let next_hs = v8_store(next_h);
            for lane in 0..targets.len() {
                if mask[lane] == 0 {
                    continue;
                }
                if scores[lane] == i8::MAX {
                    overflow_mask |= 1 << lane;
                }
                let target = targets[lane];
                let width = (target.d_end - target.d_begin) as usize;
                let lower = (target.d_begin + subject_pos[lane] as i32).max(0) as usize;
                let idx = (subject_pos[lane] + 1) * (width + 2) + qpos - lower;
                trace[lane][idx] = u8::from(scores[lane] > i8::MIN) * ACTIVE
                    | u8::from(scores[lane] == verticals[lane]) * GAP_V
                    | u8::from(scores[lane] == horizontals[lane]) * GAP_H
                    | u8::from(next_vs[lane] == opens[lane]) * OPEN_V
                    | u8::from(next_hs[lane] == opens[lane]) * OPEN_H;
                let j = subject_pos[lane] + 1;
                if scores[lane] > best[lane] || (scores[lane] == best[lane] && j == best_j[lane]) {
                    best[lane] = scores[lane];
                    best_i[lane] = qpos + 1;
                    best_j[lane] = j;
                }
            }
            vertical = v8_select(maskv, next_v, zero);
        }
        std::mem::swap(&mut ph, &mut ch);
        std::mem::swap(&mut pe, &mut ce);
    }
    let scores: Vec<i32> = best[..targets.len()]
        .iter()
        .map(|&x| i32::from(x) - i32::from(i8::MIN))
        .collect();
    NarrowTraceBatch {
        results: finish_results(query, targets, &trace, &scores, &best_i, &best_j),
        overflow_mask,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn trace_i16_core(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> NarrowTraceBatch {
    const LANES: usize = 8;
    let (band, i0, subject_start, band_offset, cols) = geometry::<LANES>(query.len(), targets);
    let zero = v16_splat(i16::MIN);
    let mut go_lanes = [0i16; LANES];
    let mut ge_lanes = [0i16; LANES];
    let mut overflow_mask = 0u32;
    for lane in 0..targets.len() {
        let scale = targets[lane].matrix_scale.max(1);
        let go = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=16_000).contains(&go) || !(0..=16_000).contains(&ge) {
            overflow_mask |= 1 << lane;
        }
        go_lanes[lane] = go.clamp(0, 16_000) as i16;
        ge_lanes[lane] = ge.clamp(0, 16_000) as i16;
    }
    let go = v16_load(go_lanes);
    let ge = v16_load(ge_lanes);
    let mut ph = vec![zero; band];
    let mut ch = vec![zero; band];
    let mut pe = vec![zero; band + 1];
    let mut ce = vec![zero; band + 1];
    let mut trace = allocate_trace(targets);
    let mut best = [i16::MIN; LANES];
    let mut best_i = [0; LANES];
    let mut best_j = [0; LANES];

    for column in 0..cols {
        ch.fill(zero);
        ce.fill(zero);
        let mut vertical = zero;
        let mut live = [false; LANES];
        let mut subject_pos = [0; LANES];
        for lane in 0..targets.len() {
            let pos = subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                live[lane] = true;
                subject_pos[lane] = pos as usize;
            }
        }
        for r in 0..band {
            let q = i0 + column as i32 + r as i32;
            if q < 0 || q >= query.len() as i32 {
                vertical = zero;
                continue;
            }
            let qpos = q as usize;
            let qletter = (query[qpos] & LETTER_MASK) as usize;
            let mut mask = [0i16; LANES];
            let mut subst = [0i16; LANES];
            for lane in 0..targets.len() {
                if !live[lane] || r < band_offset[lane] {
                    continue;
                }
                mask[lane] = -1;
                let sl = targets[lane].subject[subject_pos[lane]];
                let raw = if let Some(adjusted) = targets[lane].matrix {
                    adjusted.scores[(sl & LETTER_MASK) as usize * 32 + qletter]
                } else {
                    matrix.matrix8()[qletter * 32 + (sl & LETTER_MASK) as usize]
                };
                subst[lane] = i16::from(raw)
                    + if targets[lane].matrix.is_none() {
                        i16::from(cbs.get(qpos).copied().unwrap_or(0))
                    } else {
                        0
                    };
            }
            let maskv = v16_load(mask);
            let diag = v16_adds(ph[r], v16_load(subst));
            let horizontal = pe[r + 1];
            let score = v16_select(
                maskv,
                v16_max(v16_max(v16_max(diag, horizontal), vertical), zero),
                zero,
            );
            let open = v16_subs(score, go);
            let next_h = v16_max(v16_subs(horizontal, ge), open);
            let next_v = v16_max(v16_subs(vertical, ge), open);
            ch[r] = score;
            ce[r] = v16_select(maskv, next_h, zero);
            let scores = v16_store(score);
            let verticals = v16_store(vertical);
            let horizontals = v16_store(horizontal);
            let opens = v16_store(open);
            let next_vs = v16_store(next_v);
            let next_hs = v16_store(next_h);
            for lane in 0..targets.len() {
                if mask[lane] == 0 {
                    continue;
                }
                if scores[lane] == i16::MAX {
                    overflow_mask |= 1 << lane;
                }
                let target = targets[lane];
                let width = (target.d_end - target.d_begin) as usize;
                let lower = (target.d_begin + subject_pos[lane] as i32).max(0) as usize;
                let idx = (subject_pos[lane] + 1) * (width + 2) + qpos - lower;
                trace[lane][idx] = u8::from(scores[lane] > i16::MIN) * ACTIVE
                    | u8::from(scores[lane] == verticals[lane]) * GAP_V
                    | u8::from(scores[lane] == horizontals[lane]) * GAP_H
                    | u8::from(next_vs[lane] == opens[lane]) * OPEN_V
                    | u8::from(next_hs[lane] == opens[lane]) * OPEN_H;
                let j = subject_pos[lane] + 1;
                if scores[lane] > best[lane] || (scores[lane] == best[lane] && j == best_j[lane]) {
                    best[lane] = scores[lane];
                    best_i[lane] = qpos + 1;
                    best_j[lane] = j;
                }
            }
            vertical = v16_select(maskv, next_v, zero);
        }
        std::mem::swap(&mut ph, &mut ch);
        std::mem::swap(&mut pe, &mut ce);
    }
    let scores: Vec<i32> = best[..targets.len()]
        .iter()
        .map(|&x| i32::from(x) - i32::from(i16::MIN))
        .collect();
    NarrowTraceBatch {
        results: finish_results(query, targets, &trace, &scores, &best_i, &best_j),
        overflow_mask,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
fn allocate_trace(targets: &[TraceTarget<'_>]) -> Vec<Vec<u8>> {
    targets
        .iter()
        .map(|t| {
            let width = (t.d_end - t.d_begin) as usize;
            vec![0; (t.subject.len() + 1) * (width + 2)]
        })
        .collect()
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
fn finish_results<const LANES: usize>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    trace: &[Vec<u8>],
    scores: &[i32],
    best_i: &[usize; LANES],
    best_j: &[usize; LANES],
) -> Vec<SwResult> {
    targets
        .iter()
        .enumerate()
        .map(|(lane, target)| {
            if scores[lane] == 0 {
                return SwResult::default();
            }
            let rows = (target.d_end - target.d_begin) as usize + 2;
            let at = |i: usize, j: usize| {
                if i == 0 || j == 0 || j > target.subject.len() {
                    return 0;
                }
                let lower = (target.d_begin + j as i32 - 1).max(0) as usize + 1;
                if i < lower || i - lower >= rows {
                    0
                } else {
                    trace[lane][j * rows + i - lower]
                }
            };
            let (mut i, mut j) = (best_i[lane], best_j[lane]);
            let mut out = SwResult {
                score: scores[lane],
                query_end: i as i32,
                subject_end: j as i32,
                ..Default::default()
            };
            let mut operations = Vec::new();
            while i > 0 && j > 0 && at(i, j) & ACTIVE != 0 {
                let flags = at(i, j);
                if flags & GAP_V != 0 {
                    let mut n = 0;
                    loop {
                        n += 1;
                        i -= 1;
                        if i == 0 || at(i, j) & OPEN_V != 0 {
                            break;
                        }
                    }
                    operations.push((EditOperation::Insertion, n));
                    out.gap_openings += 1;
                    out.gaps += n;
                    out.length += n;
                } else if flags & GAP_H != 0 {
                    let mut n = 0;
                    loop {
                        n += 1;
                        j -= 1;
                        if j == 0 || at(i, j) & OPEN_H != 0 {
                            break;
                        }
                    }
                    operations.push((EditOperation::Deletion, n));
                    out.gap_openings += 1;
                    out.gaps += n;
                    out.length += n;
                } else {
                    if query[i - 1] & LETTER_MASK == target.subject[j - 1] & LETTER_MASK {
                        operations.push((EditOperation::Match, 1));
                        out.identities += 1;
                    } else {
                        operations.push((EditOperation::Substitution, 1));
                        out.mismatches += 1;
                    }
                    out.length += 1;
                    i -= 1;
                    j -= 1;
                }
            }
            out.query_begin = i as i32;
            out.subject_begin = j as i32;
            operations.reverse();
            out.operations = operations;
            out
        })
        .collect()
}

// The small wrappers below keep the recurrence shared while letting each
// architecture compile it under the exact runtime-checked feature set.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
unsafe fn v8_splat(x: i8) -> V8 {
    arch::_mm_set1_epi8(x)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v8_splat(x: i8) -> V8 {
    arch::vdupq_n_s8(x)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
unsafe fn v8_load(x: [i8; 16]) -> V8 {
    arch::_mm_loadu_si128(x.as_ptr().cast())
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v8_load(x: [i8; 16]) -> V8 {
    arch::vld1q_s8(x.as_ptr())
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
unsafe fn v8_store(v: V8) -> [i8; 16] {
    let mut x = [0; 16];
    arch::_mm_storeu_si128(x.as_mut_ptr().cast(), v);
    x
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v8_store(v: V8) -> [i8; 16] {
    let mut x = [0; 16];
    arch::vst1q_s8(x.as_mut_ptr(), v);
    x
}

macro_rules! x86_v8 {
    ($name:ident,$op:ident) => {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = "sse4.1")]
        unsafe fn $name(a: V8, b: V8) -> V8 {
            arch::$op(a, b)
        }
    };
}
x86_v8!(v8_adds, _mm_adds_epi8);
x86_v8!(v8_subs, _mm_subs_epi8);
x86_v8!(v8_max, _mm_max_epi8);
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
unsafe fn v8_select(m: V8, y: V8, n: V8) -> V8 {
    arch::_mm_blendv_epi8(n, y, m)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v8_adds(a: V8, b: V8) -> V8 {
    arch::vqaddq_s8(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v8_subs(a: V8, b: V8) -> V8 {
    arch::vqsubq_s8(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v8_max(a: V8, b: V8) -> V8 {
    arch::vmaxq_s8(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v8_select(m: V8, y: V8, n: V8) -> V8 {
    arch::vbslq_s8(arch::vreinterpretq_u8_s8(m), y, n)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn v16_splat(x: i16) -> V16 {
    arch::_mm_set1_epi16(x)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v16_splat(x: i16) -> V16 {
    arch::vdupq_n_s16(x)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn v16_load(x: [i16; 8]) -> V16 {
    arch::_mm_loadu_si128(x.as_ptr().cast())
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v16_load(x: [i16; 8]) -> V16 {
    arch::vld1q_s16(x.as_ptr())
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn v16_store(v: V16) -> [i16; 8] {
    let mut x = [0; 8];
    arch::_mm_storeu_si128(x.as_mut_ptr().cast(), v);
    x
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v16_store(v: V16) -> [i16; 8] {
    let mut x = [0; 8];
    arch::vst1q_s16(x.as_mut_ptr(), v);
    x
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn v16_adds(a: V16, b: V16) -> V16 {
    arch::_mm_adds_epi16(a, b)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn v16_subs(a: V16, b: V16) -> V16 {
    arch::_mm_subs_epi16(a, b)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn v16_max(a: V16, b: V16) -> V16 {
    let m = arch::_mm_cmpgt_epi16(a, b);
    arch::_mm_or_si128(arch::_mm_and_si128(m, a), arch::_mm_andnot_si128(m, b))
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn v16_select(m: V16, y: V16, n: V16) -> V16 {
    arch::_mm_or_si128(arch::_mm_and_si128(m, y), arch::_mm_andnot_si128(m, n))
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v16_adds(a: V16, b: V16) -> V16 {
    arch::vqaddq_s16(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v16_subs(a: V16, b: V16) -> V16 {
    arch::vqsubq_s16(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v16_max(a: V16, b: V16) -> V16 {
    arch::vmaxq_s16(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn v16_select(m: V16, y: V16, n: V16) -> V16 {
    arch::vbslq_s16(arch::vreinterpretq_u16_s16(m), y, n)
}

#[cfg(all(test, any(target_arch = "x86", target_arch = "x86_64")))]
mod tests {
    use super::*;
    use crate::dp::swipe::{banded_sw_cbs_range, TracebackScratch};

    fn assert_same(a: &SwResult, b: &SwResult) {
        assert_eq!(a.score, b.score);
        assert_eq!(
            (a.query_begin, a.query_end, a.subject_begin, a.subject_end),
            (b.query_begin, b.query_end, b.subject_begin, b.subject_end)
        );
        assert_eq!(a.operations, b.operations);
        assert_eq!(
            (a.identities, a.mismatches, a.gap_openings, a.gaps, a.length),
            (b.identities, b.mismatches, b.gap_openings, b.gaps, b.length)
        );
    }

    #[test]
    fn forced_sse_narrow_randomized_matches_scalar() {
        if !std::arch::is_x86_feature_detected!("sse4.1") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0x1280_16ee_5eed_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        for bin in [0, 1] {
            for case in 0..48 {
                let qlen = 8 + next() as usize % 45;
                let query: Vec<Letter> = (0..qlen)
                    .map(|i| (next() % 20) as i8 | if (i + case) % 11 == 0 { i8::MIN } else { 0 })
                    .collect();
                let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
                let count = 3 + next() as usize % if bin == 0 { 14 } else { 6 };
                let subjects: Vec<Vec<Letter>> = (0..count)
                    .map(|lane| {
                        let n = 4 + next() as usize % 43;
                        (0..n)
                            .map(|i| {
                                (next() % 20) as i8
                                    | if (i + lane + case) % 13 == 0 {
                                        i8::MIN
                                    } else {
                                        0
                                    }
                            })
                            .collect()
                    })
                    .collect();
                let bands: Vec<_> = subjects
                    .iter()
                    .map(|s| {
                        if case % 4 == 0 {
                            (-(s.len() as i32 - 1), qlen as i32)
                        } else {
                            let b = (next() % 15) as i32 - 7;
                            (b, b + 1 + (next() % 19) as i32)
                        }
                    })
                    .collect();
                let targets: Vec<_> = subjects
                    .iter()
                    .zip(&bands)
                    .map(|(s, &(b, e))| TraceTarget {
                        subject: s,
                        d_begin: b,
                        d_end: e,
                        matrix: None,
                        matrix_scale: 1,
                    })
                    .collect();
                let got = if bin == 0 {
                    trace_batch_i8(&query, &targets, &matrix, &cbs)
                } else {
                    trace_batch_i16(&query, &targets, &matrix, &cbs)
                }
                .unwrap();
                for lane in 0..count {
                    if got.overflow_mask & (1 << lane) != 0 {
                        continue;
                    }
                    let expected = banded_sw_cbs_range(
                        &query,
                        &subjects[lane],
                        bands[lane].0,
                        bands[lane].1,
                        &matrix,
                        &cbs,
                        None,
                        1,
                        &mut TracebackScratch::default(),
                    );
                    assert_same(&got.results[lane], &expected);
                }
            }
        }
    }

    #[test]
    fn forced_sse_overflow_and_ties() {
        if !std::arch::is_x86_feature_detected!("sse4.1") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let (letter, score) = (0i8..20)
            .map(|x| (x, matrix.score(x, x)))
            .max_by_key(|x| x.1)
            .unwrap();
        let n = 255usize.div_ceil(score as usize) + 2;
        let query = vec![letter; n];
        let subject = query.clone();
        let target = TraceTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
            matrix: None,
            matrix_scale: 1,
        };
        assert_eq!(
            trace_batch_i8(&query, &[target], &matrix, &[])
                .unwrap()
                .overflow_mask,
            1
        );
        let word = trace_batch_i16(&query, &[target], &matrix, &[]).unwrap();
        assert_eq!(word.overflow_mask, 0);
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
        assert_same(&word.results[0], &expected);

        let n = 65_535usize.div_ceil(score as usize) + 2;
        let query = vec![letter; n];
        let subject = query.clone();
        let target = TraceTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
            matrix: None,
            matrix_scale: 1,
        };
        assert_eq!(
            trace_batch_i16(&query, &[target], &matrix, &[])
                .unwrap()
                .overflow_mask,
            1
        );
    }
}
