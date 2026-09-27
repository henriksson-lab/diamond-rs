//! 16-lane signed-byte SWIPE tier for 128-bit SIMD targets.

use super::simd_score::ScoreTarget;
use super::simd_score8::BatchScores8;
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
use crate::basic::value::{LETTER_MASK, SEED_MASK};
use crate::stats::score_matrix::ScoreMatrix;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64 as arch;
#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type V = arch::__m128i;
#[cfg(target_arch = "aarch64")]
type V = arch::int8x16_t;
#[derive(Default)]
pub struct PortableScratch8 {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
    prev_h: Vec<V>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
    curr_h: Vec<V>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
    prev_e: Vec<V>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
    curr_e: Vec<V>,
}

pub fn available() -> bool {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        return std::arch::is_x86_feature_detected!("ssse3");
    }
    #[cfg(target_arch = "aarch64")]
    {
        return std::arch::is_aarch64_feature_detected!("neon");
    }
    #[allow(unreachable_code)]
    false
}

fn args(query: &[Letter], len: usize, matrix: &ScoreMatrix, cbs: &[i8]) -> Option<(i8, i8)> {
    if len == 0 || len > 16 || (!cbs.is_empty() && cbs.len() < query.len()) {
        return None;
    }
    let go = matrix.gap_open().checked_add(matrix.gap_extend())?;
    let ge = matrix.gap_extend();
    if !(0..=63).contains(&go) || !(0..=63).contains(&ge) {
        return None;
    }
    Some((go as i8, ge as i8))
}

pub fn score_batch_portable_i8(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    scratch: &mut PortableScratch8,
) -> Option<BatchScores8> {
    let (go, ge) = args(query, targets.len(), matrix, cbs)?;
    dispatch_banded(
        query,
        targets,
        matrix.matrix8(),
        cbs,
        go,
        ge,
        semi_global,
        scratch,
    )
}

pub fn score_full_batch_portable_i8(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    scratch: &mut PortableScratch8,
) -> Option<BatchScores8> {
    let (go, ge) = args(query, targets.len(), matrix, cbs)?;
    dispatch_full(
        query,
        targets,
        matrix.matrix8(),
        cbs,
        go,
        ge,
        semi_global,
        scratch,
    )
}

fn dispatch_banded(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &[i8; 1024],
    cbs: &[i8],
    go: i8,
    ge: i8,
    semi: bool,
    scratch: &mut PortableScratch8,
) -> Option<BatchScores8> {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("sse4.1") {
            return Some(unsafe {
                banded_sse41(query, targets, matrix, cbs, go, ge, semi, scratch)
            });
        }
        if std::arch::is_x86_feature_detected!("ssse3") {
            return Some(unsafe {
                banded_ssse3(query, targets, matrix, cbs, go, ge, semi, scratch)
            });
        }
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("neon") {
        return Some(unsafe { banded_core(query, targets, matrix, cbs, go, ge, semi, scratch) });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
    let _ = (query, targets, matrix, cbs, go, ge, semi, scratch);
    None
}

fn dispatch_full(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &[i8; 1024],
    cbs: &[i8],
    go: i8,
    ge: i8,
    semi: bool,
    scratch: &mut PortableScratch8,
) -> Option<BatchScores8> {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("sse4.1") {
            return Some(unsafe { full_sse41(query, targets, matrix, cbs, go, ge, semi, scratch) });
        }
        if std::arch::is_x86_feature_detected!("ssse3") {
            return Some(unsafe { full_ssse3(query, targets, matrix, cbs, go, ge, semi, scratch) });
        }
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("neon") {
        return Some(unsafe { full_core(query, targets, matrix, cbs, go, ge, semi, scratch) });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
    let _ = (query, targets, matrix, cbs, go, ge, semi, scratch);
    None
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
unsafe fn banded_sse41(
    q: &[Letter],
    t: &[ScoreTarget<'_>],
    m: &[i8; 1024],
    c: &[i8],
    go: i8,
    ge: i8,
    s: bool,
    x: &mut PortableScratch8,
) -> BatchScores8 {
    banded_core(q, t, m, c, go, ge, s, x)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn banded_ssse3(
    q: &[Letter],
    t: &[ScoreTarget<'_>],
    m: &[i8; 1024],
    c: &[i8],
    go: i8,
    ge: i8,
    s: bool,
    x: &mut PortableScratch8,
) -> BatchScores8 {
    banded_core(q, t, m, c, go, ge, s, x)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
unsafe fn full_sse41(
    q: &[Letter],
    t: &[&[Letter]],
    m: &[i8; 1024],
    c: &[i8],
    go: i8,
    ge: i8,
    s: bool,
    x: &mut PortableScratch8,
) -> BatchScores8 {
    full_core(q, t, m, c, go, ge, s, x)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn full_ssse3(
    q: &[Letter],
    t: &[&[Letter]],
    m: &[i8; 1024],
    c: &[i8],
    go: i8,
    ge: i8,
    s: bool,
    x: &mut PortableScratch8,
) -> BatchScores8 {
    full_core(q, t, m, c, go, ge, s, x)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn banded_core(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &[i8; 1024],
    cbs: &[i8],
    go: i8,
    ge: i8,
    semi: bool,
    scratch: &mut PortableScratch8,
) -> BatchScores8 {
    let delta = if semi { 0 } else { i8::MIN };
    let band = targets
        .iter()
        .map(|t| (t.d_end - t.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    if band == 0 {
        return empty(targets.len());
    }
    let zero = splat(delta);
    let neg = splat(i8::MIN);
    prepare(scratch, band + 1, zero, neg);
    let gov = splat(go);
    let gev = splat(ge);
    let maxv = splat(i8::MAX);
    let mut best = zero;
    let mut overflow = splat(0);
    let max_len = targets.iter().map(|t| t.subject.len()).max().unwrap_or(0);
    for j in 0..max_len {
        let (subject, seeded) = pack_banded_column(targets, j);
        let profile_v = build_profile(matrix, subject, seeded);
        let mut profile = [[0i8; 16]; 32];
        for (letter, vector) in profile_v.into_iter().enumerate() {
            profile[letter] = store(vector);
        }
        scratch.curr_h[band] = zero;
        scratch.curr_e[band] = neg;
        let mut vertical = neg;
        for r in 0..band {
            let mut subst = [0i8; 16];
            let mut valid = [0i8; 16];
            let mut ov = [0i8; 16];
            for lane in 0..targets.len() {
                let t = targets[lane];
                let q = t.d_begin + j as i32 + r as i32;
                if j >= t.subject.len()
                    || r >= (t.d_end - t.d_begin).max(0) as usize
                    || q < 0
                    || q >= query.len() as i32
                {
                    continue;
                }
                let q = q as usize;
                let value = profile[(query[q] & LETTER_MASK) as usize][lane] as i32
                    + cbs.get(q).copied().unwrap_or(0) as i32;
                if !(i8::MIN as i32..=i8::MAX as i32).contains(&value) {
                    ov[lane] = -1;
                } else {
                    valid[lane] = -1;
                }
                subst[lane] = value.clamp(i8::MIN as i32, i8::MAX as i32) as i8;
            }
            let mask = load(valid);
            overflow = or(overflow, load(ov));
            let diag = adds(scratch.prev_h[r], load(subst));
            let horizontal = scratch.prev_e[r + 1];
            let mut h = max(diag, horizontal);
            h = max(max(h, vertical), zero);
            h = select(mask, h, zero);
            overflow = or(overflow, eq(h, maxv));
            let opened = subs(h, gov);
            let e = max(subs(horizontal, gev), opened);
            vertical = max(subs(vertical, gev), opened);
            scratch.curr_h[r] = h;
            scratch.curr_e[r] = select(mask, e, neg);
            vertical = select(mask, vertical, neg);
            best = max(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }
    finish(best, overflow, targets.len(), delta)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn full_core(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &[i8; 1024],
    cbs: &[i8],
    go: i8,
    ge: i8,
    semi: bool,
    scratch: &mut PortableScratch8,
) -> BatchScores8 {
    let delta = if semi { 0 } else { i8::MIN };
    let zero = splat(delta);
    let neg = splat(i8::MIN);
    prepare(scratch, query.len() + 1, zero, neg);
    let gov = splat(go);
    let gev = splat(ge);
    let maxv = splat(i8::MAX);
    let mut best = zero;
    let mut overflow = splat(0);
    let max_len = targets.iter().map(|t| t.len()).max().unwrap_or(0);
    for j in 0..max_len {
        let (subject, valid, seeded) = pack_full_column(targets, j);
        let profile = build_profile(matrix, subject, seeded);
        scratch.curr_h[0] = zero;
        scratch.curr_e[0] = neg;
        let mut vertical = neg;
        for (q, &ql) in query.iter().enumerate() {
            let mask = valid;
            let cbs_v = splat(cbs.get(q).copied().unwrap_or(0));
            let base = profile[(ql & LETTER_MASK) as usize];
            overflow = or(overflow, select(mask, add_overflow(base, cbs_v), splat(0)));
            let diag = adds(scratch.prev_h[q], adds(base, cbs_v));
            let horizontal = scratch.prev_e[q + 1];
            let mut h = max(diag, horizontal);
            h = max(max(h, vertical), zero);
            h = select(mask, h, zero);
            overflow = or(overflow, eq(h, maxv));
            let opened = subs(h, gov);
            let e = max(subs(horizontal, gev), opened);
            vertical = max(subs(vertical, gev), opened);
            scratch.curr_h[q + 1] = h;
            scratch.curr_e[q + 1] = select(mask, e, neg);
            vertical = select(mask, vertical, neg);
            best = max(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }
    finish(best, overflow, targets.len(), delta)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn pack_full_column(targets: &[&[Letter]], column: usize) -> (V, V, V) {
    let mut subject = [0i8; 16];
    let mut valid = [0i8; 16];
    let mut seeded = [0i8; 16];
    for lane in 0..targets.len() {
        if let Some(&letter) = targets[lane].get(column) {
            subject[lane] = (letter & LETTER_MASK) as i8;
            valid[lane] = -1;
            seeded[lane] = if letter & SEED_MASK != 0 { -1 } else { 0 };
        }
    }
    (load(subject), load(valid), load(seeded))
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn pack_banded_column(targets: &[ScoreTarget<'_>], column: usize) -> (V, V) {
    let mut subject = [0i8; 16];
    let mut seeded = [0i8; 16];
    for lane in 0..targets.len() {
        if let Some(&letter) = targets[lane].subject.get(column) {
            subject[lane] = (letter & LETTER_MASK) as i8;
            seeded[lane] = if letter & SEED_MASK != 0 { -1 } else { 0 };
        }
    }
    (load(subject), load(seeded))
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn build_profile(matrix: &[i8; 1024], subject: V, seeded: V) -> [V; 32] {
    let index = arch::_mm_and_si128(subject, arch::_mm_set1_epi8(15));
    let high_mask = arch::_mm_cmpgt_epi8(subject, arch::_mm_set1_epi8(15));
    let zero = arch::_mm_setzero_si128();
    let mut profile = [zero; 32];
    for (query_letter, slot) in profile.iter_mut().enumerate() {
        let row = matrix.as_ptr().add(query_letter * 32);
        let low = arch::_mm_loadu_si128(row.cast());
        let high = arch::_mm_loadu_si128(row.add(16).cast());
        let lo_score = arch::_mm_shuffle_epi8(low, index);
        let hi_score = arch::_mm_shuffle_epi8(high, index);
        let score = select(high_mask, hi_score, lo_score);
        *slot = select(seeded, zero, score);
    }
    profile
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn build_profile(matrix: &[i8; 1024], subject: V, seeded: V) -> [V; 32] {
    let index = arch::vandq_s8(subject, arch::vdupq_n_s8(15));
    let high_mask = arch::vcgtq_s8(subject, arch::vdupq_n_s8(15));
    let zero = arch::vdupq_n_s8(0);
    let mut profile = [zero; 32];
    for (query_letter, slot) in profile.iter_mut().enumerate() {
        let row = matrix.as_ptr().add(query_letter * 32);
        let low = arch::vld1q_s8(row);
        let high = arch::vld1q_s8(row.add(16));
        let lo_score = arch::vqtbl1q_s8(low, arch::vreinterpretq_u8_s8(index));
        let hi_score = arch::vqtbl1q_s8(high, arch::vreinterpretq_u8_s8(index));
        let score = arch::vbslq_s8(high_mask, hi_score, lo_score);
        *slot = arch::vbslq_s8(arch::vreinterpretq_u8_s8(seeded), zero, score);
    }
    profile
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn add_overflow(a: V, b: V) -> V {
    let sum = arch::_mm_add_epi8(a, b);
    let bits = arch::_mm_and_si128(arch::_mm_xor_si128(a, sum), arch::_mm_xor_si128(b, sum));
    arch::_mm_cmpgt_epi8(arch::_mm_setzero_si128(), bits)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn add_overflow(a: V, b: V) -> V {
    let sum = arch::vaddq_s8(a, b);
    let bits = arch::vandq_s8(arch::veorq_s8(a, sum), arch::veorq_s8(b, sum));
    arch::vreinterpretq_s8_u8(arch::vcltq_s8(bits, arch::vdupq_n_s8(0)))
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
fn empty(len: usize) -> BatchScores8 {
    BatchScores8 {
        scores: [0; 32],
        overflow_mask: 0,
        len,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn prepare(s: &mut PortableScratch8, n: usize, z: V, neg: V) {
    s.prev_h.resize(n, z);
    s.curr_h.resize(n, z);
    s.prev_e.resize(n, neg);
    s.curr_e.resize(n, neg);
    s.prev_h.fill(z);
    s.prev_e.fill(neg);
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn finish(best: V, overflow: V, len: usize, delta: i8) -> BatchScores8 {
    let raw = store(best);
    let ov = store(overflow);
    let mut scores = [0; 32];
    let mut mask = 0u32;
    for i in 0..len {
        scores[i] = raw[i] as i32 - delta as i32;
        if ov[i] != 0 {
            mask |= 1 << i;
        }
    }
    BatchScores8 {
        scores,
        overflow_mask: mask,
        len,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn splat(x: i8) -> V {
    arch::_mm_set1_epi8(x)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn load(x: [i8; 16]) -> V {
    arch::_mm_loadu_si128(x.as_ptr().cast())
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn store(v: V) -> [i8; 16] {
    let mut x = [0; 16];
    arch::_mm_storeu_si128(x.as_mut_ptr().cast(), v);
    x
}
macro_rules! xb {
    ($n:ident,$op:ident) => {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = "sse2")]
        unsafe fn $n(a: V, b: V) -> V {
            arch::$op(a, b)
        }
    };
}
xb!(adds, _mm_adds_epi8);
xb!(subs, _mm_subs_epi8);
xb!(or, _mm_or_si128);
xb!(eq, _mm_cmpeq_epi8);
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn max(a: V, b: V) -> V {
    let m = arch::_mm_cmpgt_epi8(a, b);
    select(m, a, b)
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn select(m: V, y: V, n: V) -> V {
    arch::_mm_or_si128(arch::_mm_and_si128(m, y), arch::_mm_andnot_si128(m, n))
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn splat(x: i8) -> V {
    arch::vdupq_n_s8(x)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn load(x: [i8; 16]) -> V {
    arch::vld1q_s8(x.as_ptr())
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn store(v: V) -> [i8; 16] {
    let mut x = [0; 16];
    arch::vst1q_s8(x.as_mut_ptr(), v);
    x
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn adds(a: V, b: V) -> V {
    arch::vqaddq_s8(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn subs(a: V, b: V) -> V {
    arch::vqsubq_s8(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn max(a: V, b: V) -> V {
    arch::vmaxq_s8(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn or(a: V, b: V) -> V {
    arch::vorrq_s8(a, b)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn eq(a: V, b: V) -> V {
    arch::vreinterpretq_s8_u8(arch::vceqq_s8(a, b))
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn select(m: V, y: V, n: V) -> V {
    arch::vbslq_s8(arch::vreinterpretq_u8_s8(m), y, n)
}

#[cfg(all(
    test,
    any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")
))]
mod tests {
    use super::*;
    fn scalar(query: &[Letter], target: ScoreTarget<'_>, matrix: &ScoreMatrix, cbs: &[i8]) -> i32 {
        let neg = i32::MIN / 4;
        let go = matrix.gap_open() + matrix.gap_extend();
        let ge = matrix.gap_extend();
        let mut ph = vec![0; query.len() + 1];
        let mut pe = vec![neg; query.len() + 1];
        let mut best = 0;
        for (j, &sl) in target.subject.iter().enumerate() {
            let mut ch = vec![0; query.len() + 1];
            let mut ce = vec![neg; query.len() + 1];
            let mut f = neg;
            for i in 1..=query.len() {
                let q = i - 1;
                if q as i32 >= target.d_begin + j as i32 && (q as i32) < target.d_end + j as i32 {
                    let subst = if sl & SEED_MASK != 0 {
                        0
                    } else {
                        matrix.score(query[q] & LETTER_MASK, sl & LETTER_MASK)
                    } + cbs.get(q).copied().unwrap_or(0) as i32;
                    let h = (ph[i - 1] + subst).max(pe[i]).max(f).max(0);
                    ch[i] = h;
                    ce[i] = (pe[i] - ge).max(h - go);
                    f = (f - ge).max(h - go);
                    best = best.max(h);
                }
            }
            ph = ch;
            pe = ce;
        }
        best
    }
    fn cases(mut check: impl FnMut(&[Letter], &[ScoreTarget<'_>], &ScoreMatrix, &[i8], bool)) {
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0x8bad_f00du64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        for semi in [false, true] {
            for _ in 0..48 {
                let qlen = 4 + next() as usize % 67;
                let query: Vec<_> = (0..qlen).map(|_| (next() % 20) as Letter).collect();
                let cbs: Vec<_> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
                let n = 1 + next() as usize % 16;
                let subjects: Vec<Vec<Letter>> = (0..n)
                    .map(|_| {
                        let l = 3 + next() as usize % 69;
                        (0..l)
                            .map(|_| {
                                let letter = (next() % 20) as Letter;
                                if next() % 11 == 0 {
                                    letter | SEED_MASK
                                } else {
                                    letter
                                }
                            })
                            .collect()
                    })
                    .collect();
                let bands: Vec<_> = (0..n)
                    .map(|_| {
                        let b = (next() % 19) as i32 - 9;
                        (b, b + 1 + (next() % 28) as i32)
                    })
                    .collect();
                let targets: Vec<_> = subjects
                    .iter()
                    .zip(&bands)
                    .map(|(s, &(b, e))| ScoreTarget {
                        subject: s,
                        d_begin: b,
                        d_end: e,
                    })
                    .collect();
                check(&query, &targets, &matrix, &cbs, semi);
            }
        }
    }
    fn assert_scores(
        query: &[Letter],
        targets: &[ScoreTarget<'_>],
        matrix: &ScoreMatrix,
        cbs: &[i8],
        semi: bool,
        got: BatchScores8,
    ) {
        for lane in 0..targets.len() {
            let expected = scalar(query, targets[lane], matrix, cbs);
            let ceiling = if semi { 127 } else { 255 };
            if expected >= ceiling {
                assert_ne!(
                    got.overflow_mask & (1 << lane),
                    0,
                    "lane={lane} score={expected}"
                );
            } else {
                assert_eq!(got.overflow_mask & (1 << lane), 0);
                assert_eq!(got.scores[lane], expected);
            }
        }
    }
    #[test]
    fn randomized_banded_matches_scalar_or_promotes() {
        cases(|q, t, m, c, s| {
            let got =
                score_batch_portable_i8(q, t, m, c, s, &mut PortableScratch8::default()).unwrap();
            assert_scores(q, t, m, c, s, got);
        });
    }
    #[test]
    fn randomized_full_matches_scalar_or_promotes() {
        cases(|q, t, m, c, s| {
            let refs: Vec<_> = t.iter().map(|x| x.subject).collect();
            let got =
                score_full_batch_portable_i8(q, &refs, m, c, s, &mut PortableScratch8::default())
                    .unwrap();
            let full: Vec<_> = refs
                .iter()
                .map(|&subject| ScoreTarget {
                    subject,
                    d_begin: -(subject.len() as i32 - 1),
                    d_end: q.len() as i32,
                })
                .collect();
            assert_scores(q, &full, m, c, s, got);
        });
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn forced_ssse3_and_sse41_backends() {
        let m = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 100).unwrap();
        let q = vec![0, 1, 2, 3, 4];
        let s = vec![0, 1, 2, 3];
        let t = [ScoreTarget {
            subject: &s,
            d_begin: -2,
            d_end: 5,
        }];
        if std::arch::is_x86_feature_detected!("ssse3") {
            let got = unsafe {
                banded_ssse3(
                    &q,
                    &t,
                    m.matrix8(),
                    &[],
                    12,
                    1,
                    false,
                    &mut PortableScratch8::default(),
                )
            };
            assert_scores(&q, &t, &m, &[], false, got);
        }
        if std::arch::is_x86_feature_detected!("sse4.1") {
            let got = unsafe {
                banded_sse41(
                    &q,
                    &t,
                    m.matrix8(),
                    &[],
                    12,
                    1,
                    false,
                    &mut PortableScratch8::default(),
                )
            };
            assert_scores(&q, &t, &m, &[], false, got);
        }
    }

    #[test]
    fn semi_global_overflow_promotes_losslessly_to_portable_i16() {
        use super::super::simd_score_portable::{
            score_batch_portable_i16, PortableSimdScoreScratch,
        };
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 100).unwrap();
        // W/W is 11. This exceeds the semi-global signed-byte ceiling (127)
        // but is well inside i16. The i16 recurrence uses an unshifted zero
        // baseline, which is exactly the decoded DELTA=0 recurrence.
        let query = vec![17; 20];
        let subject = vec![17; 20];
        let targets = [ScoreTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
        }];
        let low = score_batch_portable_i8(
            &query,
            &targets,
            &matrix,
            &[],
            true,
            &mut PortableScratch8::default(),
        )
        .unwrap();
        assert_eq!(low.overflow_mask, 1);
        let promoted = score_batch_portable_i16(
            &query,
            &targets,
            &matrix,
            &[],
            &mut PortableSimdScoreScratch::default(),
        )
        .unwrap();
        assert_eq!(promoted.overflow_mask, 0);
        assert_eq!(promoted.scores[0], scalar(&query, targets[0], &matrix, &[]));
        assert_eq!(promoted.scores[0], 220);
    }
}
