//! 32-lane AVX2 signed-byte score tier for banded SWIPE.

use super::simd_score::ScoreTarget;
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::basic::value::{LETTER_MASK, SEED_MASK};
use crate::stats::score_matrix::ScoreMatrix;

#[derive(Clone, Copy, Debug)]
pub struct BatchScores8 {
    pub scores: [i32; 32],
    pub overflow_mask: u32,
    pub len: usize,
}

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type V = arch::__m256i;

#[derive(Default)]
pub struct Scratch8 {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_h: Vec<V>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_h: Vec<V>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_e: Vec<V>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_e: Vec<V>,
}

pub fn available() -> bool {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        return std::arch::is_x86_feature_detected!("avx2");
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        false
    }
}

pub fn score_batch_avx2_i8(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    scratch: &mut Scratch8,
) -> Option<BatchScores8> {
    if targets.is_empty() || targets.len() > 32 || (!cbs.is_empty() && cbs.len() < query.len()) {
        return None;
    }
    let go = matrix.gap_open().checked_add(matrix.gap_extend())?;
    let ge = matrix.gap_extend();
    if !(0..=63).contains(&go) || !(0..=63).contains(&ge) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        // SAFETY: AVX2 was detected at runtime.
        Some(unsafe {
            score_impl(
                query,
                targets,
                matrix.matrix8(),
                cbs,
                go as i8,
                ge as i8,
                semi_global,
                scratch,
            )
        })
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, matrix, cbs, semi_global, scratch);
        None
    }
}

pub fn score_full_batch_avx2_i8(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    scratch: &mut Scratch8,
) -> Option<BatchScores8> {
    if targets.is_empty() || targets.len() > 32 || (!cbs.is_empty() && cbs.len() < query.len()) {
        return None;
    }
    let go = matrix.gap_open().checked_add(matrix.gap_extend())?;
    let ge = matrix.gap_extend();
    if !(0..=63).contains(&go) || !(0..=63).contains(&ge) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        // SAFETY: guarded by runtime AVX2 detection.
        Some(unsafe {
            score_full_impl(
                query,
                targets,
                matrix.matrix8(),
                cbs,
                go as i8,
                ge as i8,
                semi_global,
                scratch,
            )
        })
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, matrix, cbs, semi_global, scratch);
        None
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_full_impl(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &[i8; 1024],
    cbs: &[i8],
    go: i8,
    ge: i8,
    semi_global: bool,
    scratch: &mut Scratch8,
) -> BatchScores8 {
    let delta = if semi_global { 0i8 } else { i8::MIN };
    let zero = arch::_mm256_set1_epi8(delta);
    let neg = arch::_mm256_set1_epi8(i8::MIN);
    let rows = query.len() + 1;
    scratch.prev_h.resize(rows, zero);
    scratch.curr_h.resize(rows, zero);
    scratch.prev_e.resize(rows, neg);
    scratch.curr_e.resize(rows, neg);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(neg);
    let go_v = arch::_mm256_set1_epi8(go);
    let ge_v = arch::_mm256_set1_epi8(ge);
    let max_v = arch::_mm256_set1_epi8(i8::MAX);
    let mut best = zero;
    let mut overflow = arch::_mm256_setzero_si256();
    let max_len = targets.iter().map(|target| target.len()).max().unwrap_or(0);
    for j in 0..max_len {
        let (subject, valid, seeded) = pack_full_column(targets, j);
        let profile = build_profile(matrix, subject, seeded);
        scratch.curr_h[0] = zero;
        scratch.curr_e[0] = neg;
        let mut vertical = neg;
        for (q, &ql) in query.iter().enumerate() {
            let mask = valid;
            let cbs_v = arch::_mm256_set1_epi8(cbs.get(q).copied().unwrap_or(0));
            let base = profile[(ql & LETTER_MASK) as usize];
            overflow = arch::_mm256_or_si256(
                overflow,
                arch::_mm256_and_si256(mask, add_overflow(base, cbs_v)),
            );
            let subst = arch::_mm256_adds_epi8(base, cbs_v);
            let diag = arch::_mm256_adds_epi8(scratch.prev_h[q], subst);
            let horizontal = scratch.prev_e[q + 1];
            let mut h = arch::_mm256_max_epi8(diag, horizontal);
            h = arch::_mm256_max_epi8(h, vertical);
            h = arch::_mm256_max_epi8(h, zero);
            h = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, h),
                arch::_mm256_andnot_si256(mask, zero),
            );
            overflow = arch::_mm256_or_si256(overflow, arch::_mm256_cmpeq_epi8(h, max_v));
            let open = arch::_mm256_subs_epi8(h, go_v);
            let e = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, ge_v), open);
            vertical = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical, ge_v), open);
            scratch.curr_h[q + 1] = h;
            scratch.curr_e[q + 1] = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, e),
                arch::_mm256_andnot_si256(mask, neg),
            );
            vertical = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, vertical),
                arch::_mm256_andnot_si256(mask, neg),
            );
            best = arch::_mm256_max_epi8(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }
    finish(best, overflow, targets.len(), delta)
}

/// Build the upstream-style `SwipeProfile<int8_t>` for one target column.
/// Each entry is the score for one query alphabet letter across all lanes.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn pack_full_column(targets: &[&[Letter]], column: usize) -> (V, V, V) {
    let mut subject = [0i8; 32];
    let mut valid = [0i8; 32];
    let mut seeded = [0i8; 32];
    for lane in 0..targets.len() {
        if let Some(&letter) = targets[lane].get(column) {
            subject[lane] = (letter & LETTER_MASK) as i8;
            valid[lane] = -1;
            seeded[lane] = if letter & SEED_MASK != 0 { -1 } else { 0 };
        }
    }
    (
        arch::_mm256_loadu_si256(subject.as_ptr().cast()),
        arch::_mm256_loadu_si256(valid.as_ptr().cast()),
        arch::_mm256_loadu_si256(seeded.as_ptr().cast()),
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn pack_banded_column(targets: &[ScoreTarget<'_>], column: usize) -> (V, V) {
    let mut subject = [0i8; 32];
    let mut seeded = [0i8; 32];
    for lane in 0..targets.len() {
        if let Some(&letter) = targets[lane].subject.get(column) {
            subject[lane] = (letter & LETTER_MASK) as i8;
            seeded[lane] = if letter & SEED_MASK != 0 { -1 } else { 0 };
        }
    }
    (
        arch::_mm256_loadu_si256(subject.as_ptr().cast()),
        arch::_mm256_loadu_si256(seeded.as_ptr().cast()),
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn build_profile(matrix: &[i8; 1024], subject: V, seeded: V) -> [V; 32] {
    let low_index = arch::_mm256_and_si256(subject, arch::_mm256_set1_epi8(15));
    let high_mask = arch::_mm256_cmpgt_epi8(subject, arch::_mm256_set1_epi8(15));
    let zero = arch::_mm256_setzero_si256();
    let mut profile = [zero; 32];
    for (query_letter, slot) in profile.iter_mut().enumerate() {
        let row = matrix.as_ptr().add(query_letter * 32);
        let low = arch::_mm_loadu_si128(row.cast());
        let high = arch::_mm_loadu_si128(row.add(16).cast());
        let low = arch::_mm256_broadcastsi128_si256(low);
        let high = arch::_mm256_broadcastsi128_si256(high);
        let lo_score = arch::_mm256_shuffle_epi8(low, low_index);
        let hi_score = arch::_mm256_shuffle_epi8(high, low_index);
        let score = arch::_mm256_blendv_epi8(lo_score, hi_score, high_mask);
        *slot = arch::_mm256_andnot_si256(seeded, score);
    }
    profile
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn add_overflow(a: V, b: V) -> V {
    // Signed addition overflows iff both inputs have the same sign and the
    // wrapping result has the opposite sign.
    let sum = arch::_mm256_add_epi8(a, b);
    let bits = arch::_mm256_and_si256(
        arch::_mm256_xor_si256(a, sum),
        arch::_mm256_xor_si256(b, sum),
    );
    arch::_mm256_cmpgt_epi8(arch::_mm256_setzero_si256(), bits)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn finish(best: V, overflow: V, len: usize, delta: i8) -> BatchScores8 {
    let mut raw = [0i8; 32];
    arch::_mm256_storeu_si256(raw.as_mut_ptr().cast(), best);
    let mut scores = [0i32; 32];
    for lane in 0..len {
        scores[lane] = raw[lane] as i32 - delta as i32;
    }
    BatchScores8 {
        scores,
        overflow_mask: arch::_mm256_movemask_epi8(overflow) as u32,
        len,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_impl(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix: &[i8; 32 * 32],
    cbs: &[i8],
    go: i8,
    ge: i8,
    semi_global: bool,
    scratch: &mut Scratch8,
) -> BatchScores8 {
    let delta = if semi_global { 0i8 } else { i8::MIN };
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    if band == 0 {
        return BatchScores8 {
            scores: [0; 32],
            overflow_mask: 0,
            len: targets.len(),
        };
    }
    let zero = arch::_mm256_set1_epi8(delta);
    let neg = arch::_mm256_set1_epi8(i8::MIN);
    scratch.prev_h.resize(band + 1, zero);
    scratch.curr_h.resize(band + 1, zero);
    scratch.prev_e.resize(band + 1, neg);
    scratch.curr_e.resize(band + 1, neg);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(neg);
    let go_v = arch::_mm256_set1_epi8(go);
    let ge_v = arch::_mm256_set1_epi8(ge);
    let max_v = arch::_mm256_set1_epi8(i8::MAX);
    let mut best = zero;
    let mut overflow = arch::_mm256_setzero_si256();
    let max_len = targets
        .iter()
        .map(|target| target.subject.len())
        .max()
        .unwrap_or(0);

    for j in 0..max_len {
        let (subject, seeded) = pack_banded_column(targets, j);
        let profile_v = build_profile(matrix, subject, seeded);
        let mut profile = [[0i8; 32]; 32];
        for (letter, vector) in profile_v.into_iter().enumerate() {
            arch::_mm256_storeu_si256(profile[letter].as_mut_ptr().cast(), vector);
        }
        scratch.curr_h[band] = zero;
        scratch.curr_e[band] = neg;
        let mut vertical = neg;
        for r in 0..band {
            let mut scores = [0i8; 32];
            let mut valid = [0i8; 32];
            for lane in 0..targets.len() {
                let target = targets[lane];
                let q = target.d_begin + j as i32 + r as i32;
                if j >= target.subject.len()
                    || r >= (target.d_end - target.d_begin).max(0) as usize
                    || q < 0
                    || q >= query.len() as i32
                {
                    continue;
                }
                valid[lane] = -1;
                let q = q as usize;
                let value = profile[(query[q] & LETTER_MASK) as usize][lane] as i32
                    + cbs.get(q).copied().unwrap_or(0) as i32;
                if !(i8::MIN as i32..=i8::MAX as i32).contains(&value) {
                    valid[lane] = 0;
                    overflow = arch::_mm256_or_si256(
                        overflow,
                        arch::_mm256_cmpeq_epi8(
                            arch::_mm256_set1_epi8(lane as i8),
                            arch::_mm256_setr_epi8(
                                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                                19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                            ),
                        ),
                    );
                }
                scores[lane] = value.clamp(i8::MIN as i32, i8::MAX as i32) as i8;
            }
            let mask = arch::_mm256_loadu_si256(valid.as_ptr().cast());
            let subst = arch::_mm256_loadu_si256(scores.as_ptr().cast());
            let diag = arch::_mm256_adds_epi8(scratch.prev_h[r], subst);
            let horizontal = scratch.prev_e[r + 1];
            let mut h = arch::_mm256_max_epi8(diag, horizontal);
            h = arch::_mm256_max_epi8(h, vertical);
            h = arch::_mm256_max_epi8(h, zero);
            h = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, h),
                arch::_mm256_andnot_si256(mask, zero),
            );
            overflow = arch::_mm256_or_si256(overflow, arch::_mm256_cmpeq_epi8(h, max_v));
            let open = arch::_mm256_subs_epi8(h, go_v);
            let e = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, ge_v), open);
            vertical = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical, ge_v), open);
            scratch.curr_h[r] = h;
            scratch.curr_e[r] = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, e),
                arch::_mm256_andnot_si256(mask, neg),
            );
            vertical = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, vertical),
                arch::_mm256_andnot_si256(mask, neg),
            );
            best = arch::_mm256_max_epi8(best, h);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }
    let mut raw = [0i8; 32];
    arch::_mm256_storeu_si256(raw.as_mut_ptr().cast(), best);
    let mut scores = [0i32; 32];
    for lane in 0..targets.len() {
        scores[lane] = raw[lane] as i32 - delta as i32;
    }
    BatchScores8 {
        scores,
        overflow_mask: arch::_mm256_movemask_epi8(overflow) as u32,
        len: targets.len(),
    }
}

#[cfg(all(test, any(target_arch = "x86", target_arch = "x86_64")))]
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

    #[test]
    fn randomized_i8_matches_scalar_or_marks_overflow() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 17u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        let mut scratch = Scratch8::default();
        for _ in 0..50 {
            let qlen = 8 + next() as usize % 70;
            let query: Vec<_> = (0..qlen).map(|_| (next() % 20) as Letter).collect();
            let cbs: Vec<_> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
            let count = 1 + next() as usize % 32;
            let subjects: Vec<Vec<Letter>> = (0..count)
                .map(|_| {
                    let n = 5 + next() as usize % 70;
                    (0..n)
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
            let bands: Vec<_> = (0..count)
                .map(|_| {
                    let b = (next() % 17) as i32 - 8;
                    (b, b + 1 + (next() % 24) as i32)
                })
                .collect();
            let targets: Vec<_> = subjects
                .iter()
                .zip(&bands)
                .map(|(subject, &(d_begin, d_end))| ScoreTarget {
                    subject,
                    d_begin,
                    d_end,
                })
                .collect();
            let got =
                score_batch_avx2_i8(&query, &targets, &matrix, &cbs, false, &mut scratch).unwrap();
            for lane in 0..count {
                let expected = scalar(&query, targets[lane], &matrix, &cbs);
                if expected >= u8::MAX as i32 {
                    assert_ne!(got.overflow_mask & (1 << lane), 0);
                } else {
                    assert_eq!(
                        got.overflow_mask & (1 << lane),
                        0,
                        "lane={lane} expected={expected} got={}",
                        got.scores[lane]
                    );
                    assert_eq!(got.scores[lane], expected);
                }
            }
        }
    }

    #[test]
    fn randomized_full_i8_matches_scalar_or_marks_overflow() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query: Vec<Letter> = (0..61).map(|i| (i * 7 % 20) as Letter).collect();
        let subjects: Vec<Vec<Letter>> = (1..=32)
            .map(|lane| {
                (0..(7 + lane))
                    .map(|i| {
                        let letter = ((i * 11 + lane) % 20) as Letter;
                        if (i + lane) % 13 == 0 {
                            letter | SEED_MASK
                        } else {
                            letter
                        }
                    })
                    .collect()
            })
            .collect();
        let refs: Vec<&[Letter]> = subjects.iter().map(Vec::as_slice).collect();
        let got =
            score_full_batch_avx2_i8(&query, &refs, &matrix, &[], false, &mut Scratch8::default())
                .unwrap();
        for lane in 0..refs.len() {
            let expected = scalar(
                &query,
                ScoreTarget {
                    subject: refs[lane],
                    d_begin: -(refs[lane].len() as i32 - 1),
                    d_end: query.len() as i32,
                },
                &matrix,
                &[],
            );
            if expected >= 255 {
                assert_ne!(got.overflow_mask & (1 << lane), 0);
            } else {
                assert_eq!(got.scores[lane], expected);
            }
        }
    }

    /// Focused substitution-lookup benchmark. Run with:
    /// `cargo test --release profile_lookup_microbench -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn profile_lookup_microbench() {
        use std::hint::black_box;
        use std::time::Instant;

        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query: Vec<Letter> = (0..256).map(|i| (i * 13 % 20) as Letter).collect();
        let subjects: Vec<Vec<Letter>> = (0..32)
            .map(|lane| vec![((lane * 7) % 20) as Letter])
            .collect();
        let refs: Vec<&[Letter]> = subjects.iter().map(Vec::as_slice).collect();
        let rounds = 50_000;

        let start = Instant::now();
        let mut scalar_sum = 0i64;
        for _ in 0..rounds {
            for &q in &query {
                for target in &refs {
                    scalar_sum += matrix.matrix8()
                        [(q & LETTER_MASK) as usize * 32 + target[0] as usize]
                        as i64;
                }
            }
            black_box(scalar_sum);
        }
        let scalar = start.elapsed();

        let start = Instant::now();
        for _ in 0..rounds {
            let (subject, _, seeded) = unsafe { pack_full_column(&refs, 0) };
            let profile = unsafe { build_profile(matrix.matrix8(), subject, seeded) };
            for &q in &query {
                black_box(profile[(q & LETTER_MASK) as usize]);
            }
        }
        let packed = start.elapsed();
        eprintln!(
            "scalar_lookup={scalar:?} packed_profile={packed:?} speedup={:.2}x checksum={scalar_sum}",
            scalar.as_secs_f64() / packed.as_secs_f64()
        );
    }
}
