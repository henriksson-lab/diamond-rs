//! AVX2 score-only Smith-Waterman over up to sixteen target lanes.
//!
//! The lanes share a query, but may have different lengths and diagonal
//! bands.  Invalid band cells are masked before their horizontal/vertical
//! gap state can escape into a valid cell.  Scores which reach the `i16`
//! ceiling are reported in `overflow_mask` and must be recomputed with the
//! scalar `i32` kernel.

use crate::basic::value::{Letter, AMINO_ACID_COUNT, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;
use std::mem::MaybeUninit;

/// One target lane for [`score_batch_avx2`].
#[derive(Clone, Copy, Debug)]
pub struct ScoreTarget<'a> {
    pub subject: &'a [Letter],
    pub d_begin: i32,
    pub d_end: i32,
}

/// Result from one SIMD batch. Only the first `len` entries are meaningful.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BatchScores {
    pub scores: [i32; 16],
    pub overflow_mask: u16,
    pub len: usize,
}

/// Reusable aligned DP rows for the AVX2 kernel.
///
/// Keeping this at the caller avoids allocating three query-sized rows for
/// every group of sixteen targets.
#[derive(Default)]
pub struct SimdScoreScratch {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_h: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_h: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    prev_e: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    curr_e: Vec<ArchVector>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    query_biases: Vec<ArchVector>,
}

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
type ArchVector = arch::__m256i;

/// Score at most sixteen targets with AVX2.
///
/// Returns `None` when AVX2 is unavailable, the batch is empty or too large,
/// or the penalties cannot be represented safely by this `i16` kernel.
/// Lanes selected by `overflow_mask` have to be recomputed in `i32`.
/// Adjusted target matrices use different score scaling and must stay on the
/// scalar path; this entry point always uses the standard matrix scale.
pub fn score_batch_avx2(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut SimdScoreScratch,
) -> Option<BatchScores> {
    if targets.is_empty()
        || targets.len() > 16
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
    {
        return None;
    }
    let query_letters = query.iter().try_fold(0u32, |mask, &letter| {
        let letter = (letter & LETTER_MASK) as usize;
        (letter < AMINO_ACID_COUNT).then_some(mask | (1 << letter))
    })?;
    let gap_open = score_matrix
        .gap_open()
        .checked_add(score_matrix.gap_extend())?;
    let gap_extend = score_matrix.gap_extend();
    // NEG below is -16384. Keeping penalties below this limit means all
    // deliberately unreachable-state subtraction remains representable.
    if !(0..=16_000).contains(&gap_open) || !(0..=16_000).contains(&gap_extend) {
        return None;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        // SAFETY: guarded by runtime AVX2 detection; all memory accesses in
        // the implementation are bounds checked or use fixed-size arrays.
        return Some(unsafe {
            score_batch_avx2_impl(
                query,
                targets,
                score_matrix.matrix8u_low(),
                score_matrix.matrix8u_high(),
                score_matrix.bias(),
                query_cbs,
                gap_open as i16,
                gap_extend as i16,
                query_letters,
                scratch,
            )
        });
    }

    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, score_matrix, query_cbs, scratch);
        None
    }
}

/// Score up to sixteen complete Smith-Waterman matrices in AVX2 lanes.
///
/// This is the SIMD counterpart of C++ `dp/swipe/full_swipe.h` for the
/// score-only configuration. Unlike [`score_batch_avx2`], it walks query rows
/// directly and therefore does not spend work on the triangular padding of a
/// synthetic full-width diagonal band.
pub fn score_full_batch_avx2(
    query: &[Letter],
    targets: &[&[Letter]],
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut SimdScoreScratch,
) -> Option<BatchScores> {
    if targets.is_empty()
        || targets.len() > 16
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
    {
        return None;
    }
    let gap_open = score_matrix
        .gap_open()
        .checked_add(score_matrix.gap_extend())?;
    let gap_extend = score_matrix.gap_extend();
    if !(0..=16_000).contains(&gap_open) || !(0..=16_000).contains(&gap_extend) {
        return None;
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        // SAFETY: guarded by runtime AVX2 detection.
        return Some(unsafe {
            score_full_batch_avx2_impl(
                query,
                targets,
                score_matrix.matrix16(),
                query_cbs,
                gap_open as i16,
                gap_extend as i16,
                scratch,
            )
        });
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, score_matrix, query_cbs, scratch);
        None
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_full_batch_avx2_impl(
    query: &[Letter],
    targets: &[&[Letter]],
    matrix: &[i16; 32 * 32],
    query_cbs: &[i8],
    gap_open: i16,
    gap_extend: i16,
    scratch: &mut SimdScoreScratch,
) -> BatchScores {
    const NEG: i16 = -16_384;
    let rows = query.len() + 1;
    let zero = arch::_mm256_setzero_si256();
    let neg = arch::_mm256_set1_epi16(NEG);
    scratch.prev_h.resize(rows, zero);
    scratch.curr_h.resize(rows, zero);
    scratch.prev_e.resize(rows, neg);
    scratch.curr_e.resize(rows, neg);
    scratch.prev_h.fill(zero);
    scratch.prev_e.fill(neg);
    let go = arch::_mm256_set1_epi16(gap_open);
    let ge = arch::_mm256_set1_epi16(gap_extend);

    let max_i16 = arch::_mm256_set1_epi16(i16::MAX);
    let mut best = zero;
    let mut overflow = zero;
    let max_subject_len = targets.iter().map(|target| target.len()).max().unwrap_or(0);

    for j in 0..max_subject_len {
        scratch.curr_h[0] = zero;
        scratch.curr_e[0] = neg;
        let mut vertical = neg;
        for (qpos, &ql) in query.iter().enumerate() {
            let mut subst = [0i16; 16];
            let mut valid = [0i16; 16];
            for lane in 0..targets.len() {
                if j >= targets[lane].len() {
                    continue;
                }
                valid[lane] = -1;
                let sl = targets[lane][j];
                let base = if sl & crate::basic::value::SEED_MASK != 0 {
                    0
                } else {
                    matrix[((ql & crate::basic::value::LETTER_MASK) as usize) * 32
                        + (sl & crate::basic::value::LETTER_MASK) as usize]
                        as i32
                };
                let value = base + query_cbs.get(qpos).copied().unwrap_or(0) as i32;
                if !(i16::MIN as i32..=i16::MAX as i32).contains(&value) {
                    let mut lane_mask = [0i16; 16];
                    lane_mask[lane] = -1;
                    overflow = arch::_mm256_or_si256(
                        overflow,
                        arch::_mm256_loadu_si256(lane_mask.as_ptr().cast()),
                    );
                }
                subst[lane] = value.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
            }
            let mask = arch::_mm256_loadu_si256(valid.as_ptr().cast());
            let substitution = arch::_mm256_loadu_si256(subst.as_ptr().cast());
            let diag = arch::_mm256_adds_epi16(scratch.prev_h[qpos], substitution);
            let horizontal = scratch.prev_e[qpos + 1];
            let mut score = arch::_mm256_max_epi16(diag, horizontal);
            score = arch::_mm256_max_epi16(score, vertical);
            score = arch::_mm256_max_epi16(score, zero);
            score = arch::_mm256_and_si256(score, mask);
            overflow = arch::_mm256_or_si256(overflow, arch::_mm256_cmpeq_epi16(score, max_i16));
            let open = arch::_mm256_subs_epi16(score, go);
            let next_horizontal =
                arch::_mm256_max_epi16(arch::_mm256_subs_epi16(horizontal, ge), open);
            vertical = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(vertical, ge), open);
            scratch.curr_h[qpos + 1] = score;
            scratch.curr_e[qpos + 1] = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, next_horizontal),
                arch::_mm256_andnot_si256(mask, neg),
            );
            vertical = arch::_mm256_or_si256(
                arch::_mm256_and_si256(mask, vertical),
                arch::_mm256_andnot_si256(mask, neg),
            );
            best = arch::_mm256_max_epi16(best, score);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
    }

    let mut raw_scores = [0i16; 16];
    let mut raw_overflow = [0i16; 16];
    arch::_mm256_storeu_si256(raw_scores.as_mut_ptr().cast(), best);
    arch::_mm256_storeu_si256(raw_overflow.as_mut_ptr().cast(), overflow);
    let mut scores = [0i32; 16];
    let mut overflow_mask = 0u16;
    for lane in 0..targets.len() {
        scores[lane] = raw_scores[lane] as i32;
        if raw_overflow[lane] != 0 {
            overflow_mask |= 1 << lane;
        }
    }
    BatchScores {
        scores,
        overflow_mask,
        len: targets.len(),
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn score_batch_avx2_impl(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix_low: &[i8; 32 * 32],
    matrix_high: &[i8; 32 * 32],
    matrix_bias: i8,
    query_cbs: &[i8],
    gap_open: i16,
    gap_extend: i16,
    query_letters: u32,
    scratch: &mut SimdScoreScratch,
) -> BatchScores {
    score_batch_avx2_upstream_impl(
        query,
        targets,
        matrix_low,
        matrix_high,
        matrix_bias,
        query_cbs,
        gap_open,
        gap_extend,
        query_letters,
        scratch,
    )
}

/// Direct structural port of `dp/swipe/banded_swipe.h` for the AVX2 i16
/// score-only specialization. The common moving band is essential: it keeps
/// the query coordinate uniform across lanes and the substitution profile in
/// vector form throughout the cell loop.
#[target_feature(enable = "avx2")]
unsafe fn score_batch_avx2_upstream_impl(
    query: &[Letter],
    targets: &[ScoreTarget<'_>],
    matrix_low: &[i8; 32 * 32],
    matrix_high: &[i8; 32 * 32],
    matrix_bias: i8,
    query_cbs: &[i8],
    gap_open: i16,
    gap_extend: i16,
    _query_letters: u32,
    scratch: &mut SimdScoreScratch,
) -> BatchScores {
    const LANES: usize = 16;
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    if band == 0 {
        return BatchScores {
            scores: [0; LANES],
            overflow_mask: 0,
            len: targets.len(),
        };
    }
    let i1 = targets
        .iter()
        .map(|target| (target.d_end - 1).max(0))
        .min()
        .unwrap_or(0);
    let i0 = i1 + 1 - band as i32;
    let mut subject_start = [0i32; LANES];
    let mut band_offset = [0usize; LANES];
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

    let zero = arch::_mm256_setzero_si256();
    let dp_zero = arch::_mm256_set1_epi16(i16::MIN);
    scratch.prev_h.resize(band, dp_zero);
    scratch.prev_e.resize(band + 1, dp_zero);
    scratch.prev_h.fill(dp_zero);
    scratch.prev_e.fill(dp_zero);
    scratch.query_biases.resize(query.len(), zero);
    if query_cbs.is_empty() {
        scratch.query_biases.fill(zero);
    } else {
        for (slot, &bias) in scratch.query_biases.iter_mut().zip(query_cbs) {
            *slot = arch::_mm256_set1_epi16(bias as i16);
        }
    }
    let go = arch::_mm256_set1_epi16(gap_open);
    let ge = arch::_mm256_set1_epi16(gap_extend);

    // A query normally contains only the 20 canonical residue codes.  The
    // substitution profile is rebuilt for every target column, so avoiding
    // ambiguity-code rows which this query cannot select is worthwhile.  Do
    // this once per batch, outside the column loop.
    // Direct port of RangePartition: each segment has a constant lane mask,
    // with 0 for lanes whose strict band has begun and i16::MIN otherwise.
    let mut ordered_offsets = [(0usize, 0usize); LANES];
    for lane in 0..targets.len() {
        ordered_offsets[lane] = (band_offset[lane], lane);
    }
    ordered_offsets[..targets.len()].sort_unstable();
    let mut segment_begin = [0usize; LANES];
    let mut segment_end = [0usize; LANES];
    let mut segment_masks = [dp_zero; LANES];
    let mut segment_count = 0usize;
    let mut lanes = [i16::MIN; LANES];
    let mut cursor = 0usize;
    while cursor < targets.len() {
        let begin = ordered_offsets[cursor].0;
        while cursor < targets.len() && ordered_offsets[cursor].0 == begin {
            lanes[ordered_offsets[cursor].1] = 0;
            cursor += 1;
        }
        segment_begin[segment_count] = begin;
        segment_masks[segment_count] = arch::_mm256_loadu_si256(lanes.as_ptr().cast());
        segment_count += 1;
    }
    for segment in 0..segment_count {
        segment_end[segment] = if segment + 1 < segment_count {
            segment_begin[segment + 1]
        } else {
            band
        };
    }
    let mut best = dp_zero;

    for column in 0..columns {
        let mut subject = [0i16; LANES];
        let mut inactive_lanes = [i16::MIN; LANES];
        let mut seeded = [0i16; LANES];
        for lane in 0..targets.len() {
            let pos = subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                let letter = targets[lane].subject[pos as usize];
                subject[lane] = (letter & crate::basic::value::LETTER_MASK) as i16;
                inactive_lanes[lane] = 0;
                seeded[lane] = if letter & crate::basic::value::SEED_MASK != 0 {
                    -1
                } else {
                    0
                };
            }
        }
        let subject = arch::_mm256_loadu_si256(subject.as_ptr().cast());
        let inactive = arch::_mm256_loadu_si256(inactive_lanes.as_ptr().cast());
        let seeded = arch::_mm256_loadu_si256(seeded.as_ptr().cast());
        let mut profile: [MaybeUninit<ArchVector>; AMINO_ACID_COUNT] =
            [const { MaybeUninit::uninit() }; AMINO_ACID_COUNT];
        // Direct port of AVX2 ScoreVector<int16_t>(letter, subject-vector):
        // select the low/high 16-residue table with byte shuffles, expand the
        // biased u8 scores to i16, then remove the matrix bias. This replaces
        // 512 scalar indexed loads per target column.
        let high_mask = arch::_mm256_slli_epi16(
            arch::_mm256_and_si256(subject, arch::_mm256_set1_epi8(0x10)),
            3,
        );
        let seq_low = arch::_mm256_or_si256(subject, high_mask);
        let seq_high = arch::_mm256_or_si256(
            subject,
            arch::_mm256_xor_si256(high_mask, arch::_mm256_set1_epi8(i8::MIN)),
        );
        let byte_mask = arch::_mm256_set1_epi16(255);
        let bias = arch::_mm256_set1_epi16(matrix_bias as i16);
        // Match upstream's fixed 26-row profile loop. Iterating a sparse query
        // mask looked attractive, but tzcnt/blsr and loop-control overhead
        // retired about 0.2% more instructions end-to-end.
        for query_letter in 0..AMINO_ACID_COUNT {
            let row = query_letter * 32;
            // ScoreMatrix mirrors C++ `alignas(32)` and each row is exactly
            // 32 bytes, so these are aligned loads just like upstream.
            let low = arch::_mm256_load_si256(matrix_low.as_ptr().add(row).cast());
            let high = arch::_mm256_load_si256(matrix_high.as_ptr().add(row).cast());
            let lo_score = arch::_mm256_shuffle_epi8(low, seq_low);
            let hi_score = arch::_mm256_shuffle_epi8(high, seq_high);
            let expanded =
                arch::_mm256_and_si256(arch::_mm256_or_si256(lo_score, hi_score), byte_mask);
            profile
                .get_unchecked_mut(query_letter)
                .write(arch::_mm256_andnot_si256(
                    seeded,
                    arch::_mm256_subs_epi16(expanded, bias),
                ));
        }

        let moving_i0 = i0 + column as i32;
        let query_begin = moving_i0.max(0);
        let query_end = (i1 + column as i32).min(query.len() as i32 - 1) + 1;
        let query_row_begin = (query_begin - moving_i0) as usize;
        let query_row_end = (query_end - moving_i0) as usize;
        let mut vertical = dp_zero;
        let mut col_best = dp_zero;
        for segment in 0..segment_count {
            let row_begin = segment_begin[segment].max(query_row_begin);
            let row_end = segment_end[segment].min(query_row_end);
            if row_begin >= row_end {
                continue;
            }
            let target_mask = arch::_mm256_or_si256(segment_masks[segment], inactive);
            vertical = arch::_mm256_adds_epi16(vertical, target_mask);
            for row in row_begin..row_end {
                let q = (moving_i0 + row as i32) as usize;
                // `query_begin/query_end` clamp q to the query, every segment
                // is clamped to the allocated band, and the letter mask is
                // 0..31. Express those already-established invariants here:
                // otherwise LLVM leaves four panic bounds checks in every DP
                // cell, unlike the pointer-based upstream kernel.
                let query_letter = *query.get_unchecked(q) & crate::basic::value::LETTER_MASK;
                let base = *profile
                    .get_unchecked(query_letter as usize)
                    .assume_init_ref();
                let substitution = arch::_mm256_adds_epi16(
                    arch::_mm256_adds_epi16(base, *scratch.query_biases.get_unchecked(q)),
                    target_mask,
                );
                let diagonal =
                    arch::_mm256_adds_epi16(*scratch.prev_h.get_unchecked(row), substitution);
                let horizontal =
                    arch::_mm256_adds_epi16(*scratch.prev_e.get_unchecked(row + 1), target_mask);
                let mut score = arch::_mm256_max_epi16(diagonal, horizontal);
                score = arch::_mm256_max_epi16(score, vertical);
                let open = arch::_mm256_subs_epi16(score, go);
                let next_horizontal =
                    arch::_mm256_max_epi16(arch::_mm256_subs_epi16(horizontal, ge), open);
                vertical = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(vertical, ge), open);
                *scratch.prev_h.get_unchecked_mut(row) = score;
                *scratch.prev_e.get_unchecked_mut(row) = next_horizontal;
                col_best = arch::_mm256_max_epi16(col_best, score);
            }
        }
        best = arch::_mm256_max_epi16(best, col_best);
    }

    let mut raw_scores = [0i16; LANES];
    let mut raw_overflow = [0i16; LANES];
    arch::_mm256_storeu_si256(raw_scores.as_mut_ptr().cast(), best);
    let overflow = arch::_mm256_cmpeq_epi16(best, arch::_mm256_set1_epi16(i16::MAX));
    arch::_mm256_storeu_si256(raw_overflow.as_mut_ptr().cast(), overflow);
    let mut scores = [0i32; LANES];
    let mut overflow_mask = 0u16;
    for lane in 0..targets.len() {
        scores[lane] = raw_scores[lane] as i32 - i16::MIN as i32;
        if raw_overflow[lane] != 0 {
            overflow_mask |= 1 << lane;
        }
    }
    BatchScores {
        scores,
        overflow_mask,
        len: targets.len(),
    }
}

#[cfg(all(test, any(target_arch = "x86", target_arch = "x86_64")))]
mod tests {
    use super::*;
    use crate::basic::value::{LETTER_MASK, SEED_MASK};

    fn scalar(query: &[Letter], target: ScoreTarget<'_>, matrix: &ScoreMatrix, cbs: &[i8]) -> i32 {
        let qlen = query.len();
        let neg = i32::MIN / 4;
        let go = matrix.gap_open() + matrix.gap_extend();
        let ge = matrix.gap_extend();
        let mut ph = vec![0; qlen + 1];
        let mut pe = vec![neg; qlen + 1];
        let mut best = 0;
        for (spos, &s) in target.subject.iter().enumerate() {
            let mut ch = vec![0; qlen + 1];
            let mut ce = vec![neg; qlen + 1];
            let mut f = neg;
            for i in 1..=qlen {
                let qpos = i - 1;
                let valid = qpos as i32 >= target.d_begin + spos as i32
                    && (qpos as i32) < target.d_end + spos as i32;
                if !valid {
                    continue;
                }
                let subst = if s & SEED_MASK != 0 {
                    0
                } else {
                    matrix.score(query[qpos] & LETTER_MASK, s & LETTER_MASK)
                } + if cbs.is_empty() { 0 } else { cbs[qpos] as i32 };
                let h = (ph[i - 1] + subst).max(pe[i]).max(f).max(0);
                ch[i] = h;
                ce[i] = (pe[i] - ge).max(h - go);
                f = (f - ge).max(h - go);
                best = best.max(h);
            }
            ph = ch;
            pe = ce;
        }
        best
    }

    #[test]
    fn randomized_avx2_matches_i32_scalar_with_strict_bands() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0x5eed_cafe_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        let mut scratch = SimdScoreScratch::default();
        for _ in 0..80 {
            let qlen = 1 + next() as usize % 73;
            let query: Vec<Letter> = (0..qlen).map(|_| (next() % 25) as Letter).collect();
            let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
            let lane_count = 1 + next() as usize % 16;
            let mut subjects = Vec::with_capacity(lane_count);
            let mut bands = Vec::with_capacity(lane_count);
            for _ in 0..lane_count {
                let len = 1 + next() as usize % 81;
                let mut subject: Vec<Letter> = (0..len).map(|_| (next() % 25) as Letter).collect();
                if next() % 3 == 0 {
                    let pos = next() as usize % len;
                    subject[pos] |= SEED_MASK;
                }
                subjects.push(subject);
                let begin = next() as i32 % 31 - 15;
                let width = 1 + next() as i32 % 40;
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
            let got = score_batch_avx2(&query, &targets, &matrix, &cbs, &mut scratch).unwrap();
            assert_eq!(got.overflow_mask, 0);
            for lane in 0..lane_count {
                assert_eq!(
                    got.scores[lane],
                    scalar(&query, targets[lane], &matrix, &cbs),
                    "lane {lane}, qlen {qlen}, band {:?}",
                    bands[lane]
                );
            }
        }
    }

    #[test]
    fn reports_i16_score_overflow_per_lane() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        // The banded kernel uses DIAMOND's SHRT_MIN-biased representation,
        // so W/W scores of 11 can accumulate through 65,535 before overflow.
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
        let got = score_batch_avx2(
            &query,
            &targets,
            &matrix,
            &[],
            &mut SimdScoreScratch::default(),
        )
        .unwrap();
        assert_eq!(got.overflow_mask, 1);
        assert_eq!(got.scores[1], 220);
    }

    #[test]
    fn rejects_query_codes_outside_the_amino_acid_profile() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let subject = [0];
        let target = [ScoreTarget {
            subject: &subject,
            d_begin: 0,
            d_end: 1,
        }];
        assert!(score_batch_avx2(
            &[AMINO_ACID_COUNT as Letter],
            &target,
            &matrix,
            &[],
            &mut SimdScoreScratch::default(),
        )
        .is_none());
    }

    #[test]
    fn randomized_full_matrix_avx2_matches_scalar() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0xd1a0_600du64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        let mut scratch = SimdScoreScratch::default();
        for _ in 0..48 {
            let qlen = 1 + next() as usize % 65;
            let query: Vec<Letter> = (0..qlen).map(|_| (next() % 25) as Letter).collect();
            let cbs: Vec<i8> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
            let count = 1 + next() as usize % 16;
            let subjects: Vec<Vec<Letter>> = (0..count)
                .map(|_| {
                    let len = 1 + next() as usize % 67;
                    (0..len).map(|_| (next() % 25) as Letter).collect()
                })
                .collect();
            let refs: Vec<&[Letter]> = subjects.iter().map(Vec::as_slice).collect();
            let got = score_full_batch_avx2(&query, &refs, &matrix, &cbs, &mut scratch).unwrap();
            assert_eq!(got.overflow_mask, 0);
            for lane in 0..count {
                let target = ScoreTarget {
                    subject: refs[lane],
                    d_begin: -(refs[lane].len() as i32 - 1),
                    d_end: qlen as i32,
                };
                assert_eq!(got.scores[lane], scalar(&query, target, &matrix, &cbs));
            }
        }
    }
}
