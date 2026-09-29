//! Eight-lane 128-bit SIMD traceback and adjusted-score SWIPE kernels.

use super::simd_trace::TraceTarget;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
use crate::basic::value::LETTER_MASK;
use crate::dp::smith_waterman::SwResult;
use crate::stats::score_matrix::ScoreMatrix;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64 as arch;
#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[derive(Clone, Copy)]
struct V(arch::__m128i, arch::__m128i);
#[cfg(target_arch = "aarch64")]
#[derive(Clone, Copy)]
struct V(arch::int32x4_t, arch::int32x4_t);
fn valid(query: &[Letter], targets: &[TraceTarget<'_>], cbs: &[i8]) -> bool {
    !targets.is_empty()
        && targets.len() <= 8
        && (cbs.is_empty() || cbs.len() >= query.len())
        && targets
            .iter()
            .all(|t| t.d_end > t.d_begin && (0..=i32::MAX / 8).contains(&t.matrix_scale))
}

/// Trace up to eight banded (or full-band) targets with 128-bit SIMD.
pub fn trace_batch_portable(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> Option<Vec<SwResult>> {
    if !valid(query, targets, cbs) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("ssse3") {
            return Some(unsafe { trace_ssse3(query, targets, matrix, cbs) });
        }
        if std::arch::is_x86_feature_detected!("sse2") {
            return Some(unsafe { trace_impl(query, targets, matrix, cbs) });
        }
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("neon") {
        return Some(unsafe { trace_impl(query, targets, matrix, cbs) });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
    let _ = matrix;
    None
}

/// Score lanes with per-target adjusted matrices and gap scales.
pub fn score_adjusted_batch_portable(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> Option<Vec<i32>> {
    trace_batch_portable(query, targets, matrix, cbs)
        .map(|results| results.into_iter().map(|result| result.score).collect())
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn trace_ssse3(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> Vec<SwResult> {
    trace_impl(query, targets, matrix, cbs)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
#[cfg_attr(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature(enable = "sse2")
)]
#[cfg_attr(target_arch = "aarch64", target_feature(enable = "neon"))]
unsafe fn trace_impl(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
) -> Vec<SwResult> {
    const NEG: i32 = i32::MIN / 4;
    const ACTIVE: u8 = 1;
    const GAP_V: u8 = 2;
    const GAP_H: u8 = 4;
    const OPEN_V: u8 = 8;
    const OPEN_H: u8 = 16;
    let band = targets
        .iter()
        .map(|t| (t.d_end - t.d_begin) as usize)
        .max()
        .unwrap();
    let max_len = targets.iter().map(|t| t.subject.len()).max().unwrap_or(0);
    let zero = splat(0);
    let neg = splat(NEG);
    let mut ph = vec![zero; band];
    let mut ch = vec![zero; band];
    let mut pe = vec![neg; band + 1];
    let mut ce = vec![neg; band + 1];
    let mut trace: Vec<Vec<u8>> = targets
        .iter()
        .map(|t| vec![0; (t.subject.len() + 1) * ((t.d_end - t.d_begin) as usize + 2)])
        .collect();
    let mut best = [0i32; 8];
    let mut best_i = [0usize; 8];
    let mut best_j = [0usize; 8];
    for j0 in 0..max_len {
        ch.fill(zero);
        ce.fill(neg);
        let mut vertical = neg;
        for r in 0..band {
            let mut mask = [0i32; 8];
            let mut subst = [0i32; 8];
            let mut go = [0i32; 8];
            let mut ge = [0i32; 8];
            let mut qposes = [0usize; 8];
            let mut indices = [0usize; 8];
            for lane in 0..targets.len() {
                let t = targets[lane];
                let width = (t.d_end - t.d_begin) as usize;
                let q = t.d_begin + j0 as i32 + r as i32;
                if j0 >= t.subject.len() || r >= width || q < 0 || q >= query.len() as i32 {
                    continue;
                }
                mask[lane] = -1;
                let q = q as usize;
                qposes[lane] = q;
                let lower = (t.d_begin + j0 as i32).max(0) as usize;
                indices[lane] = (j0 + 1) * (width + 2) + q + 1 - (lower + 1);
                let ql = query[q];
                let sl = t.subject[j0];
                subst[lane] = if let Some(m) = t.matrix {
                    m.scores[(sl & LETTER_MASK) as usize * 32 + (ql & LETTER_MASK) as usize] as i32
                } else {
                    matrix.score(ql & LETTER_MASK, sl & LETTER_MASK)
                        + cbs.get(q).copied().unwrap_or(0) as i32
                };
                let scale = t.matrix_scale.max(1);
                go[lane] = (matrix.gap_open() + matrix.gap_extend()) * scale;
                ge[lane] = matrix.gap_extend() * scale;
            }
            let maskv = load(mask);
            let diag = add(ph[r], load(subst));
            let horizontal = pe[r + 1];
            let mut score = max(diag, horizontal);
            score = max(score, vertical);
            score = and(max(score, zero), maskv);
            let open = sub(score, load(go));
            let next_h = max(sub(horizontal, load(ge)), open);
            let next_v = max(sub(vertical, load(ge)), open);
            ch[r] = score;
            ce[r] = select(maskv, next_h, neg);
            let active = gt(score, zero);
            let gap_v = eq(score, vertical);
            let gap_h = eq(score, horizontal);
            let open_v = eq(next_v, open);
            let open_h = eq(next_h, open);
            let scores = store(score);
            let ab = store(active);
            let gv = store(gap_v);
            let gh = store(gap_h);
            let ov = store(open_v);
            let oh = store(open_h);
            for lane in 0..targets.len() {
                if mask[lane] == 0 {
                    continue;
                }
                trace[lane][indices[lane]] = u8::from(ab[lane] != 0) * ACTIVE
                    | u8::from(gv[lane] != 0) * GAP_V
                    | u8::from(gh[lane] != 0) * GAP_H
                    | u8::from(ov[lane] != 0) * OPEN_V
                    | u8::from(oh[lane] != 0) * OPEN_H;
                let j = j0 + 1;
                if scores[lane] > best[lane] || (scores[lane] == best[lane] && j == best_j[lane]) {
                    best[lane] = scores[lane];
                    best_i[lane] = qposes[lane] + 1;
                    best_j[lane] = j;
                }
            }
            vertical = select(maskv, next_v, neg);
        }
        std::mem::swap(&mut ph, &mut ch);
        std::mem::swap(&mut pe, &mut ce);
    }
    targets
        .iter()
        .enumerate()
        .map(|(lane, t)| {
            if best[lane] == 0 {
                return SwResult::default();
            }
            let rows = (t.d_end - t.d_begin) as usize + 2;
            let at = |i: usize, j: usize| {
                if j == 0 || j > t.subject.len() || i == 0 {
                    return 0;
                }
                let lower = (t.d_begin + j as i32 - 1).max(0) as usize + 1;
                if i < lower || i - lower >= rows {
                    0
                } else {
                    trace[lane][j * rows + i - lower]
                }
            };
            let (mut i, mut j) = (best_i[lane], best_j[lane]);
            let mut out = SwResult {
                score: best[lane],
                query_end: i as i32,
                subject_end: j as i32,
                ..Default::default()
            };
            let mut ops = Vec::new();
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
                    ops.push((EditOperation::Insertion, n));
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
                    ops.push((EditOperation::Deletion, n));
                    out.gap_openings += 1;
                    out.gaps += n;
                    out.length += n;
                } else {
                    if query[i - 1] & LETTER_MASK == t.subject[j - 1] & LETTER_MASK {
                        ops.push((EditOperation::Match, 1));
                        out.identities += 1;
                    } else {
                        ops.push((EditOperation::Substitution, 1));
                        out.mismatches += 1;
                    }
                    out.length += 1;
                    i -= 1;
                    j -= 1;
                }
            }
            out.query_begin = i as i32;
            out.subject_begin = j as i32;
            ops.reverse();
            out.operations = ops;
            out
        })
        .collect()
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn splat(x: i32) -> V {
    V(arch::_mm_set1_epi32(x), arch::_mm_set1_epi32(x))
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn splat(x: i32) -> V {
    V(arch::vdupq_n_s32(x), arch::vdupq_n_s32(x))
}

macro_rules! x86_bin {
    ($n:ident,$op:ident) => {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = "sse2")]
        unsafe fn $n(a: V, b: V) -> V {
            V(arch::$op(a.0, b.0), arch::$op(a.1, b.1))
        }
    };
}
x86_bin!(add, _mm_add_epi32);
x86_bin!(sub, _mm_sub_epi32);
x86_bin!(and, _mm_and_si128);
x86_bin!(eq, _mm_cmpeq_epi32);
x86_bin!(gt, _mm_cmpgt_epi32);
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn max(a: V, b: V) -> V {
    let m0 = arch::_mm_cmpgt_epi32(a.0, b.0);
    let m1 = arch::_mm_cmpgt_epi32(a.1, b.1);
    V(
        arch::_mm_or_si128(
            arch::_mm_and_si128(m0, a.0),
            arch::_mm_andnot_si128(m0, b.0),
        ),
        arch::_mm_or_si128(
            arch::_mm_and_si128(m1, a.1),
            arch::_mm_andnot_si128(m1, b.1),
        ),
    )
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn select(m: V, y: V, n: V) -> V {
    V(
        arch::_mm_or_si128(
            arch::_mm_and_si128(m.0, y.0),
            arch::_mm_andnot_si128(m.0, n.0),
        ),
        arch::_mm_or_si128(
            arch::_mm_and_si128(m.1, y.1),
            arch::_mm_andnot_si128(m.1, n.1),
        ),
    )
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn load(x: [i32; 8]) -> V {
    V(
        arch::_mm_loadu_si128(x.as_ptr().cast()),
        arch::_mm_loadu_si128(x.as_ptr().add(4).cast()),
    )
}
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn store(v: V) -> [i32; 8] {
    let mut x = [0; 8];
    arch::_mm_storeu_si128(x.as_mut_ptr().cast(), v.0);
    arch::_mm_storeu_si128(x.as_mut_ptr().add(4).cast(), v.1);
    x
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn add(a: V, b: V) -> V {
    V(arch::vaddq_s32(a.0, b.0), arch::vaddq_s32(a.1, b.1))
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn sub(a: V, b: V) -> V {
    V(arch::vsubq_s32(a.0, b.0), arch::vsubq_s32(a.1, b.1))
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn max(a: V, b: V) -> V {
    V(arch::vmaxq_s32(a.0, b.0), arch::vmaxq_s32(a.1, b.1))
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn and(a: V, b: V) -> V {
    V(arch::vandq_s32(a.0, b.0), arch::vandq_s32(a.1, b.1))
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn eq(a: V, b: V) -> V {
    V(
        arch::vreinterpretq_s32_u32(arch::vceqq_s32(a.0, b.0)),
        arch::vreinterpretq_s32_u32(arch::vceqq_s32(a.1, b.1)),
    )
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn gt(a: V, b: V) -> V {
    V(
        arch::vreinterpretq_s32_u32(arch::vcgtq_s32(a.0, b.0)),
        arch::vreinterpretq_s32_u32(arch::vcgtq_s32(a.1, b.1)),
    )
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn select(m: V, y: V, n: V) -> V {
    V(
        arch::vbslq_s32(arch::vreinterpretq_u32_s32(m.0), y.0, n.0),
        arch::vbslq_s32(arch::vreinterpretq_u32_s32(m.1), y.1, n.1),
    )
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn load(x: [i32; 8]) -> V {
    V(
        arch::vld1q_s32(x.as_ptr()),
        arch::vld1q_s32(x.as_ptr().add(4)),
    )
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn store(v: V) -> [i32; 8] {
    let mut x = [0; 8];
    arch::vst1q_s32(x.as_mut_ptr(), v.0);
    arch::vst1q_s32(x.as_mut_ptr().add(4), v.1);
    x
}

#[cfg(all(
    test,
    any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")
))]
mod tests {
    use super::*;
    use crate::dp::swipe::{banded_sw_cbs_range, TracebackScratch};

    #[test]
    fn randomized_banded_and_full_match_scalar() {
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0x51ee_1285_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        for full in [false, true] {
            for _ in 0..64 {
                let qlen = 4 + next() as usize % 53;
                let query: Vec<Letter> = (0..qlen).map(|_| (next() % 20) as Letter).collect();
                let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
                let count = 1 + next() as usize % 8;
                let subjects: Vec<Vec<Letter>> = (0..count)
                    .map(|_| {
                        let n = 3 + next() as usize % 57;
                        (0..n).map(|_| (next() % 20) as Letter).collect()
                    })
                    .collect();
                let bands: Vec<_> = subjects
                    .iter()
                    .map(|s| {
                        if full {
                            (-(s.len() as i32 - 1), qlen as i32)
                        } else {
                            let b = (next() % 19) as i32 - 9;
                            (b, b + 1 + (next() % 25) as i32)
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
                let got = trace_batch_portable(&query, &targets, &matrix, &cbs).unwrap();
                for lane in 0..count {
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
                    assert_eq!(got[lane].score, expected.score);
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
                        "full={full} lane={lane} band={:?}",
                        bands[lane]
                    );
                    assert_eq!(got[lane].operations, expected.operations);
                    assert_eq!(
                        (
                            got[lane].identities,
                            got[lane].mismatches,
                            got[lane].gap_openings,
                            got[lane].gaps,
                            got[lane].length
                        ),
                        (
                            expected.identities,
                            expected.mismatches,
                            expected.gap_openings,
                            expected.gaps,
                            expected.length
                        )
                    );
                }
            }
        }
    }

    #[test]
    fn adjusted_score_matches_scalar_and_ignores_cbs() {
        use crate::stats::cbs::TargetMatrix;
        use std::sync::Arc;
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let subject = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut scores = [-3i8; 1024];
        for x in 0..32 {
            scores[x * 32 + x] = 7;
        }
        let adjusted = Arc::new(TargetMatrix {
            scores: scores.to_vec(),
            score_min: -3,
            score_max: 7,
        });
        let target = TraceTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
            matrix: Some(&adjusted),
            matrix_scale: 2,
        };
        let cbs = vec![50; query.len()];
        let got = score_adjusted_batch_portable(&query, &[target], &matrix, &cbs).unwrap()[0];
        let expected = banded_sw_cbs_range(
            &query,
            &subject,
            0,
            1,
            &matrix,
            &cbs,
            Some(&adjusted),
            2,
            &mut TracebackScratch::default(),
        )
        .score;
        assert_eq!(got, expected);
    }
}
