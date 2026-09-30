//! 32-lane AVX2 signed-byte score tier for banded SWIPE.

use super::simd_score::ScoreTarget;
use crate::basic::value::Letter;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::basic::value::{AMINO_ACID_COUNT, LETTER_MASK, SUPER_HARD_MASK};
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
    prev_e: Vec<V>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    cbs: Vec<V>,
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
            match (semi_global, cbs.is_empty()) {
                (true, false) => score_impl::<true, true>(
                    query,
                    targets,
                    matrix.matrix8_low(),
                    matrix.matrix8_high(),
                    cbs,
                    go as i8,
                    ge as i8,
                    scratch,
                ),
                (true, true) => score_impl::<true, false>(
                    query,
                    targets,
                    matrix.matrix8_low(),
                    matrix.matrix8_high(),
                    cbs,
                    go as i8,
                    ge as i8,
                    scratch,
                ),
                (false, false) => score_impl::<false, true>(
                    query,
                    targets,
                    matrix.matrix8_low(),
                    matrix.matrix8_high(),
                    cbs,
                    go as i8,
                    ge as i8,
                    scratch,
                ),
                (false, true) => score_impl::<false, false>(
                    query,
                    targets,
                    matrix.matrix8_low(),
                    matrix.matrix8_high(),
                    cbs,
                    go as i8,
                    ge as i8,
                    scratch,
                ),
            }
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
            match (semi_global, cbs.is_empty()) {
                (true, false) => score_full_impl::<true, true>(
                    query,
                    targets,
                    matrix.matrix8_low(),
                    matrix.matrix8_high(),
                    cbs,
                    go as i8,
                    ge as i8,
                    scratch,
                ),
                (true, true) => score_full_impl::<true, false>(
                    query,
                    targets,
                    matrix.matrix8_low(),
                    matrix.matrix8_high(),
                    cbs,
                    go as i8,
                    ge as i8,
                    scratch,
                ),
                (false, false) => score_full_impl::<false, true>(
                    query,
                    targets,
                    matrix.matrix8_low(),
                    matrix.matrix8_high(),
                    cbs,
                    go as i8,
                    ge as i8,
                    scratch,
                ),
                (false, true) => score_full_impl::<false, false>(
                    query,
                    targets,
                    matrix.matrix8_low(),
                    matrix.matrix8_high(),
                    cbs,
                    go as i8,
                    ge as i8,
                    scratch,
                ),
            }
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
unsafe fn score_full_impl<const SEMI_GLOBAL: bool, const HAS_CBS: bool>(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix_low: &[i8; 1024],
    matrix_high: &[i8; 1024],
    cbs: &[i8],
    go: i8,
    ge: i8,
    scratch: &mut Scratch8,
) -> BatchScores8 {
    let delta = if SEMI_GLOBAL { 0i8 } else { i8::MIN };
    let zero = arch::_mm256_set1_epi8(delta);
    let neg = arch::_mm256_set1_epi8(i8::MIN);
    let rows = query.len() + 1;
    scratch.prev_h.resize(rows, zero);
    scratch.prev_e.resize(rows, neg);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(neg);
    if HAS_CBS {
        scratch.cbs.clear();
        scratch.cbs.reserve(query.len());
        for &bias in cbs.get_unchecked(..query.len()) {
            scratch.cbs.push(arch::_mm256_set1_epi8(bias));
        }
    }
    let go_v = arch::_mm256_set1_epi8(go);
    let ge_v = arch::_mm256_set1_epi8(ge);
    let mut best = zero;
    let max_len = targets.iter().map(|target| target.len()).max().unwrap_or(0);
    for j in 0..max_len {
        let (subject, valid) = pack_full_column(targets, j);
        let profile = build_profile(matrix_low, matrix_high, subject);
        // Upstream SWIPE retains one score row and one horizontal-gap row,
        // updating them in place while carrying the overwritten diagonal in
        // a register. Two additional AVX2 rows are costly for long reads.
        let mut diagonal = scratch.prev_h[0];
        scratch.prev_h[0] = zero;
        scratch.prev_e[0] = neg;
        let mut vertical = neg;
        for q in 0..query.len() {
            // `prev_*` have query_len + 1 rows and q is bounded by the query.
            let next_diagonal = *scratch.prev_h.get_unchecked(q + 1);
            let mask = valid;
            let cbs_v = if HAS_CBS {
                *scratch.cbs.get_unchecked(q)
            } else {
                arch::_mm256_setzero_si256()
            };
            let base = *profile.get_unchecked((*query.get_unchecked(q) & LETTER_MASK) as usize);
            let subst = arch::_mm256_adds_epi8(base, cbs_v);
            let diag = arch::_mm256_adds_epi8(diagonal, subst);
            let horizontal = *scratch.prev_e.get_unchecked(q + 1);
            let mut h = arch::_mm256_max_epi8(diag, horizontal);
            h = arch::_mm256_max_epi8(h, vertical);
            if SEMI_GLOBAL {
                h = arch::_mm256_max_epi8(h, zero);
            }
            h = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, h),
                arch::_mm256_andnot_si256(mask, zero),
            );
            let open = arch::_mm256_subs_epi8(h, go_v);
            let e = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, ge_v), open);
            vertical = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical, ge_v), open);
            *scratch.prev_h.get_unchecked_mut(q + 1) = h;
            *scratch.prev_e.get_unchecked_mut(q + 1) = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, e),
                arch::_mm256_andnot_si256(mask, neg),
            );
            vertical = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, vertical),
                arch::_mm256_andnot_si256(mask, neg),
            );
            best = arch::_mm256_max_epi8(best, h);
            diagonal = next_diagonal;
        }
    }
    // Upstream promotes a byte lane exactly when its retained maximum reaches
    // SCHAR_MAX.  Since `best` is monotonic, testing it once is equivalent to
    // OR-ing a comparison into an overflow vector for every DP cell.
    let overflow = arch::_mm256_cmpeq_epi8(best, arch::_mm256_set1_epi8(i8::MAX));
    finish(best, overflow, targets.len(), delta)
}

/// Build the upstream-style `SwipeProfile<int8_t>` for one target column.
/// Each entry is the score for one query alphabet letter across all lanes.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn pack_full_column(targets: &[&[Letter]], column: usize) -> (V, V) {
    let mut subject = [0i8; 32];
    let mut valid = [0i8; 32];
    for lane in 0..targets.len() {
        if let Some(&letter) = targets[lane].get(column) {
            subject[lane] = (letter & LETTER_MASK) as i8;
            valid[lane] = -1;
        }
    }
    (
        arch::_mm256_loadu_si256(subject.as_ptr().cast()),
        arch::_mm256_loadu_si256(valid.as_ptr().cast()),
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn build_profile(matrix_low: &[i8; 1024], matrix_high: &[i8; 1024], subject: V) -> [V; 32] {
    // This is the AVX2 `ScoreVector<int8_t>` constructor from upstream.  The
    // preformatted tables let each query letter use two vector loads
    // and two shuffles, instead of broadcasting four 128-bit halves.  Masked
    // sequence letters retain their substitution score: SEED_MASK affects
    // seeding, while `letter_mask` strips it for DP.
    let subject = arch::_mm256_and_si256(subject, arch::_mm256_set1_epi8(LETTER_MASK));
    let high_mask = arch::_mm256_slli_epi16(
        arch::_mm256_and_si256(subject, arch::_mm256_set1_epi8(16)),
        3,
    );
    let low_index = arch::_mm256_or_si256(subject, high_mask);
    let high_index = arch::_mm256_or_si256(
        subject,
        arch::_mm256_xor_si256(high_mask, arch::_mm256_set1_epi8(i8::MIN)),
    );
    let zero = arch::_mm256_setzero_si256();
    let mut profile = [zero; 32];
    for (query_letter, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
        let low = arch::_mm256_loadu_si256(matrix_low.as_ptr().add(query_letter * 32).cast());
        let high = arch::_mm256_loadu_si256(matrix_high.as_ptr().add(query_letter * 32).cast());
        *slot = arch::_mm256_or_si256(
            arch::_mm256_shuffle_epi8(low, low_index),
            arch::_mm256_shuffle_epi8(high, high_index),
        );
    }
    profile
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
unsafe fn score_impl<const SEMI_GLOBAL: bool, const HAS_CBS: bool>(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix_low: &[i8; 32 * 32],
    matrix_high: &[i8; 32 * 32],
    cbs: &[i8],
    go: i8,
    ge: i8,
    scratch: &mut Scratch8,
) -> BatchScores8 {
    let delta = if SEMI_GLOBAL { 0i8 } else { i8::MIN };
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
    // Mirror upstream banded_swipe.h: expand every lane to a common moving
    // band. Each row then addresses one query letter for all SIMD lanes.
    let i1 = targets
        .iter()
        .map(|target| (target.d_end - 1).max(0))
        .min()
        .unwrap_or(0);
    let i0 = i1 + 1 - band as i32;
    let mut subject_start = [0i32; 32];
    let mut band_offset = [0usize; 32];
    let mut columns = 0usize;
    for lane in 0..targets.len() {
        let target = &targets[lane];
        let expanded_begin = target.d_end - band as i32;
        subject_start[lane] = i1 - (target.d_end - 1);
        band_offset[lane] = (target.d_begin - expanded_begin).max(0) as usize;
        let subject_end =
            ((query.len() as i32 - 1 - expanded_begin).min(target.subject.len() as i32 - 1) + 1)
                .max(0);
        columns = columns.max((subject_end - subject_start[lane]).max(0) as usize);
    }
    scratch.prev_h.resize(band, zero);
    scratch.prev_e.resize(band + 1, zero);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(zero);
    if HAS_CBS {
        // Match upstream CBSBuffer: pay one broadcast per query position,
        // then consume a vector directly in every DP cell. Retaining the
        // allocation in worker-local scratch avoids a per-batch allocation.
        scratch.cbs.clear();
        scratch.cbs.reserve(query.len());
        for &bias in cbs.get_unchecked(..query.len()) {
            scratch.cbs.push(arch::_mm256_set1_epi8(bias));
        }
    }
    let go_v = arch::_mm256_set1_epi8(go);
    let ge_v = arch::_mm256_set1_epi8(ge);
    let mut best = zero;
    // Match upstream STRICT_BAND: partition rows where lanes enter their true
    // bands and hold one mask in a register throughout each partition.
    let mut offsets = band_offset;
    offsets[..targets.len()].sort_unstable();
    let mut part_bounds = [0usize; 34];
    let mut part_masks = [arch::_mm256_setzero_si256(); 33];
    let mut part_count = 0usize;
    let mut part_begin = 0usize;
    for &offset in &offsets[..targets.len()] {
        let offset = offset.min(band);
        if offset <= part_begin || offset == band {
            continue;
        }
        part_bounds[part_count] = part_begin;
        let mut lanes = [i8::MIN; 32];
        for lane in 0..targets.len() {
            if part_begin >= band_offset[lane] {
                lanes[lane] = 0;
            }
        }
        part_masks[part_count] = arch::_mm256_loadu_si256(lanes.as_ptr().cast());
        part_count += 1;
        part_begin = offset;
    }
    part_bounds[part_count] = part_begin;
    let mut lanes = [i8::MIN; 32];
    for lane in 0..targets.len() {
        if part_begin >= band_offset[lane] {
            lanes[lane] = 0;
        }
    }
    part_masks[part_count] = arch::_mm256_loadu_si256(lanes.as_ptr().cast());
    part_count += 1;
    part_bounds[part_count] = band;
    for column in 0..columns {
        let mut subject = [SUPER_HARD_MASK; 32];
        for lane in 0..targets.len() {
            let pos = subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                let letter = targets[lane].subject[pos as usize];
                subject[lane] = (letter & LETTER_MASK) as i8;
            }
        }
        let subject = arch::_mm256_loadu_si256(subject.as_ptr().cast());
        let profile = build_profile(matrix_low, matrix_high, subject);
        let moving_i0 = i0 + column as i32;
        let query_begin = moving_i0.max(0);
        let query_end = (i1 + column as i32).min(query.len() as i32 - 1) + 1;
        let mut vertical = zero;
        let mut col_best = zero;
        let active_r_begin = (query_begin - moving_i0) as usize;
        let active_r_end = (query_end - moving_i0) as usize;
        for part in 0..part_count {
            let r_begin = part_bounds[part].max(active_r_begin);
            let r_end = part_bounds[part + 1].min(active_r_end);
            if r_begin >= r_end {
                continue;
            }
            let band_mask = *part_masks.get_unchecked(part);
            let q_start = (moving_i0 + r_begin as i32) as usize;
            let mut query_ptr = query.as_ptr().add(q_start);
            let mut cbs_ptr = if HAS_CBS {
                scratch.cbs.as_ptr().add(q_start)
            } else {
                std::ptr::NonNull::<V>::dangling().as_ptr()
            };
            let mut r = r_begin;
            macro_rules! update_cell {
                () => {{
                    // The band and query bounds prove all of these indices.
                    let base = arch::_mm256_adds_epi8(
                        *profile.get_unchecked((*query_ptr & LETTER_MASK) as usize),
                        band_mask,
                    );
                    let bias = if HAS_CBS {
                        *cbs_ptr
                    } else {
                        arch::_mm256_setzero_si256()
                    };
                    let subst = arch::_mm256_adds_epi8(base, bias);
                    let diag = arch::_mm256_adds_epi8(*scratch.prev_h.get_unchecked(r), subst);
                    let horizontal =
                        arch::_mm256_adds_epi8(*scratch.prev_e.get_unchecked(r + 1), band_mask);
                    let mut h = arch::_mm256_max_epi8(diag, horizontal);
                    h = arch::_mm256_max_epi8(h, vertical);
                    if SEMI_GLOBAL {
                        h = arch::_mm256_max_epi8(h, zero);
                    }
                    let open = arch::_mm256_subs_epi8(h, go_v);
                    let e = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, ge_v), open);
                    vertical = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical, ge_v), open);
                    *scratch.prev_h.get_unchecked_mut(r) = h;
                    *scratch.prev_e.get_unchecked_mut(r) = e;
                    col_best = arch::_mm256_max_epi8(col_best, h);
                    query_ptr = query_ptr.add(1);
                    if HAS_CBS {
                        cbs_ptr = cbs_ptr.add(1);
                    }
                    r += 1;
                }};
            }
            // Match the native C++ compiler's two-cell unrolling. This keeps
            // the vertical dependency order exact while halving loop-control
            // branches in the dominant recurrence.
            while r + 1 < r_end {
                update_cell!();
                update_cell!();
            }
            if r < r_end {
                update_cell!();
            }
            debug_assert_eq!(r, r_end);
        }
        best = arch::_mm256_max_epi8(best, col_best);
    }
    let mut raw = [0i8; 32];
    arch::_mm256_storeu_si256(raw.as_mut_ptr().cast(), best);
    let mut scores = [0i32; 32];
    for lane in 0..targets.len() {
        scores[lane] = raw[lane] as i32 - delta as i32;
    }
    BatchScores8 {
        scores,
        overflow_mask: arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(
            best,
            arch::_mm256_set1_epi8(i8::MAX),
        )) as u32,
        len: targets.len(),
    }
}

#[cfg(all(test, any(target_arch = "x86", target_arch = "x86_64")))]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;

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
                    let subst = matrix.score(query[q] & LETTER_MASK, sl & LETTER_MASK)
                        + cbs.get(q).copied().unwrap_or(0) as i32;
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
            let semi_global =
                score_batch_avx2_i8(&query, &targets, &matrix, &cbs, true, &mut scratch).unwrap();
            for lane in 0..count {
                let expected = scalar(&query, targets[lane], &matrix, &cbs);
                if expected >= u8::MAX as i32 {
                    assert_ne!(got.overflow_mask & (1 << lane), 0);
                    assert_ne!(semi_global.overflow_mask & (1 << lane), 0);
                } else {
                    assert_eq!(
                        got.overflow_mask & (1 << lane),
                        0,
                        "lane={lane} expected={expected} got={}",
                        got.scores[lane]
                    );
                    assert_eq!(got.scores[lane], expected);
                    assert_eq!(semi_global.scores[lane], expected);
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
        let semi_global =
            score_full_batch_avx2_i8(&query, &refs, &matrix, &[], true, &mut Scratch8::default())
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
                assert_ne!(semi_global.overflow_mask & (1 << lane), 0);
            } else {
                assert_eq!(got.scores[lane], expected);
                assert_eq!(semi_global.scores[lane], expected);
            }
        }
    }

    #[test]
    fn avx2_profile_matches_matrix_for_every_alphabet_letter() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut subject = [SUPER_HARD_MASK; 32];
        for (lane, letter) in subject[..AMINO_ACID_COUNT].iter_mut().enumerate() {
            *letter = lane as Letter | if lane % 2 == 0 { SEED_MASK } else { 0 };
        }
        unsafe {
            let packed = arch::_mm256_loadu_si256(subject.as_ptr().cast());
            let profile = build_profile(matrix.matrix8_low(), matrix.matrix8_high(), packed);
            for (query_letter, scores) in profile.iter().enumerate().take(AMINO_ACID_COUNT) {
                let mut lanes = [0i8; 32];
                arch::_mm256_storeu_si256(lanes.as_mut_ptr().cast(), *scores);
                for subject_letter in 0..AMINO_ACID_COUNT {
                    assert_eq!(
                        lanes[subject_letter],
                        matrix.score(query_letter as Letter, subject_letter as Letter) as i8,
                        "query={query_letter} subject={subject_letter}"
                    );
                }
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
            let (subject, _) = unsafe { pack_full_column(&refs, 0) };
            let profile =
                unsafe { build_profile(matrix.matrix8_low(), matrix.matrix8_high(), subject) };
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
