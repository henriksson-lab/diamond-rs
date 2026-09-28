//! Translation of `diamond/src/dp/ungapped_simd.{h,cpp}`.
//!
//! Each lane is a signed saturating i8 accumulator biased by `SCHAR_MIN`.
//! x86-64 hosts use the 32-lane AVX2 layout and ARM NEON hosts the 16-lane
//! layout from upstream; other hosts retain the architecture-neutral path.

use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

const SCORE_BIAS: i8 = i8::MIN;

/// Number of subjects handled by one production stage-2 batch.
#[inline]
pub fn preferred_lane_count() -> usize {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("avx2") {
        return 32;
    }
    16
}

/// C++ `window_ungapped(...)`, returning an owned Rust result vector.
pub fn window_ungapped(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
) -> Vec<i32> {
    let mut out = vec![0; subjects.len()];
    window_ungapped_into(query, subjects, window, score_matrix, &mut out);
    out
}

/// Direct output-buffer form of C++ `window_ungapped(...)`.
pub fn window_ungapped_into(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    validate_buffers(query, subjects, window, out);

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if subjects.len() <= 32 && std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 was detected above and validate_buffers established
        // that every input has at least `window` accessible letters.
        unsafe {
            window_ungapped_avx2(query, subjects, window, score_matrix, out);
        }
        return;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if subjects.len() <= 16
        && std::arch::is_x86_feature_detected!("ssse3")
        && std::arch::is_x86_feature_detected!("sse4.1")
    {
        // SAFETY: SSSE3 and SSE4.1 were detected above and validate_buffers
        // established that every input has `window` accessible letters.
        unsafe {
            window_ungapped_sse41(query, subjects, window, score_matrix, out);
        }
        return;
    }

    #[cfg(target_arch = "aarch64")]
    if subjects.len() <= 16 {
        // SAFETY: NEON is mandatory on AArch64; validate_buffers established
        // that all reads are in bounds.
        unsafe {
            window_ungapped_neon(query, subjects, window, score_matrix, out);
        }
        return;
    }

    for (score, subject) in out.iter_mut().zip(subjects) {
        *score = saturating_window_score(query, subject, window, score_matrix);
    }
}

/// C++ AArch64 `window_ungapped`: one subject per signed-byte NEON lane.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn window_ungapped_neon(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    use std::arch::aarch64::*;

    debug_assert!(subjects.len() <= 16);
    let mut subject_block = [0i8; 16 * 16];
    let mut subject_rows = [&[][..]; 16];
    let mut best_bytes = [SCORE_BIAS; 16];
    let mut score = vdupq_n_s8(SCORE_BIAS);
    let mut best = score;
    let low_nibble = vdupq_n_s8(0x0f);
    let high_start = vdupq_n_s8(16);

    let lane_offset = 16 - subjects.len();
    for block_begin in (0..window).step_by(16) {
        let block_len = (window - block_begin).min(16);
        if block_len == 16 {
            for (row, subject) in subject_rows.iter_mut().zip(subjects) {
                *row = &subject[block_begin..];
            }
            crate::util::simd::transpose_16(
                &subject_rows[..subjects.len()],
                subjects.len(),
                &mut subject_block,
            );
        } else {
            subject_block.fill(0);
            for position in 0..block_len {
                for (lane, subject) in subjects.iter().enumerate() {
                    subject_block[position * 16 + lane_offset + lane] =
                        subject[block_begin + position];
                }
            }
        }
        for position in 0..block_len {
            let letters = vandq_s8(
                vld1q_s8(subject_block.as_ptr().add(position * 16)),
                vdupq_n_s8(LETTER_MASK),
            );

            let query_letter = (query[block_begin + position] & LETTER_MASK) as usize;
            let row = score_matrix.matrix8().as_ptr().add(query_letter * 32);
            let match_scores = {
                let indices = vreinterpretq_u8_s8(vandq_s8(letters, low_nibble));
                let low_scores = vqtbl1q_s8(vld1q_s8(row), indices);
                let high_scores = vqtbl1q_s8(vld1q_s8(row.add(16)), indices);
                let high_mask = vcgeq_s8(letters, high_start);
                vbslq_s8(high_mask, high_scores, low_scores)
            };
            score = vqaddq_s8(score, match_scores);
            best = vmaxq_s8(best, score);
        }
    }

    vst1q_s8(best_bytes.as_mut_ptr(), best);
    for (dst, &value) in out.iter_mut().zip(&best_bytes[lane_offset..]) {
        *dst = i32::from(value) - i32::from(SCORE_BIAS);
    }
}

/// C++ `ARCH_AVX2::window_ungapped`: one subject per signed-byte lane.
///
/// Upstream transposes subject letters and performs a vector score lookup. A
/// direct gather into a stack vector would leave most of the lookup scalar, so
/// this uses two lane-local `vpshufb` tables for the 32-entry matrix row.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn window_ungapped_avx2(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    debug_assert!(subjects.len() <= 32);
    let mut subject_block = [0i8; 32 * 32];
    let mut tail_block = [[0i8; 32]; 32];
    let mut subject_ptrs = [std::ptr::null(); 32];
    for (pointer, subject) in subject_ptrs.iter_mut().zip(subjects) {
        *pointer = subject.as_ptr();
    }
    let mut best_bytes = [SCORE_BIAS; 32];
    let mut score = _mm256_set1_epi8(SCORE_BIAS);
    let mut best = score;
    let lane_offset = 32 - subjects.len();
    for block_begin in (0..window).step_by(32) {
        let block_len = (window - block_begin).min(32);
        if block_len == 32 {
            crate::util::simd::transpose_32_avx2(
                &subject_ptrs[..subjects.len()],
                subjects.len(),
                block_begin,
                &mut subject_block,
            );
        } else {
            // Upstream can transpose a complete final vector because its
            // SequenceSet allocation is padded. Rust slices do not permit an
            // overread, so make the same full-width input explicitly. Copying
            // one contiguous tail per lane avoids the former position×lane
            // loop and its two bounds checks per byte.
            for (lane, subject) in subjects.iter().enumerate() {
                tail_block[lane].fill(0);
                std::ptr::copy_nonoverlapping(
                    subject.as_ptr().add(block_begin),
                    tail_block[lane].as_mut_ptr(),
                    block_len,
                );
            }
            for lane in 0..subjects.len() {
                subject_ptrs[lane] = tail_block[lane].as_ptr();
            }
            crate::util::simd::transpose_32_avx2(
                &subject_ptrs[..subjects.len()],
                subjects.len(),
                0,
                &mut subject_block,
            );
        }
        for position in 0..block_len {
            let letters = _mm256_and_si256(
                _mm256_loadu_si256(subject_block.as_ptr().add(position * 32).cast()),
                _mm256_set1_epi8(LETTER_MASK),
            );
            let query_letter = (query[block_begin + position] & LETTER_MASK) as usize;
            let row = query_letter * 32;
            let row_low = _mm256_loadu_si256(score_matrix.matrix8_low().as_ptr().add(row).cast());
            let row_high = _mm256_loadu_si256(score_matrix.matrix8_high().as_ptr().add(row).cast());
            // Direct port of ScoreVector<int8_t>(letter, subject-vector): the
            // high bit steers each lane to exactly one of the two shuffle
            // tables, so their results can be ORed without a blend.
            let high_mask = _mm256_slli_epi16(_mm256_and_si256(letters, _mm256_set1_epi8(0x10)), 3);
            let seq_low = _mm256_or_si256(letters, high_mask);
            let seq_high = _mm256_or_si256(
                letters,
                _mm256_xor_si256(high_mask, _mm256_set1_epi8(i8::MIN)),
            );
            let low_scores = _mm256_shuffle_epi8(row_low, seq_low);
            let high_scores = _mm256_shuffle_epi8(row_high, seq_high);
            let match_scores = _mm256_or_si256(low_scores, high_scores);

            score = _mm256_adds_epi8(score, match_scores);
            best = _mm256_max_epi8(best, score);
        }
    }

    _mm256_storeu_si256(best_bytes.as_mut_ptr().cast(), best);
    for (dst, &value) in out.iter_mut().zip(&best_bytes[lane_offset..]) {
        *dst = i32::from(value) - i32::from(SCORE_BIAS);
    }
}

/// C++ `ARCH_SSE4_1::window_ungapped`: one subject per signed-byte lane.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3,sse4.1")]
unsafe fn window_ungapped_sse41(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    debug_assert!(subjects.len() <= 16);
    let mut subject_block = [0i8; 16 * 16];
    let mut subject_rows = [&[][..]; 16];
    let mut best_bytes = [SCORE_BIAS; 16];
    let mut score = _mm_set1_epi8(SCORE_BIAS);
    let mut best = score;
    let low_nibble = _mm_set1_epi8(0x0f);

    let lane_offset = 16 - subjects.len();
    for block_begin in (0..window).step_by(16) {
        let block_len = (window - block_begin).min(16);
        if block_len == 16 {
            for (row, subject) in subject_rows.iter_mut().zip(subjects) {
                *row = &subject[block_begin..];
            }
            crate::util::simd::transpose_16(
                &subject_rows[..subjects.len()],
                subjects.len(),
                &mut subject_block,
            );
        } else {
            subject_block.fill(0);
            for position in 0..block_len {
                for (lane, subject) in subjects.iter().enumerate() {
                    subject_block[position * 16 + lane_offset + lane] =
                        subject[block_begin + position];
                }
            }
        }
        for position in 0..block_len {
            let letters = _mm_and_si128(
                _mm_loadu_si128(subject_block.as_ptr().add(position * 16).cast()),
                _mm_set1_epi8(LETTER_MASK),
            );
            let indices = _mm_and_si128(letters, low_nibble);

            let query_letter = (query[block_begin + position] & LETTER_MASK) as usize;
            let row = score_matrix.matrix8().as_ptr().add(query_letter * 32);
            let row_low = _mm_loadu_si128(row.cast());
            let row_high = _mm_loadu_si128(row.add(16).cast());
            let low_scores = _mm_shuffle_epi8(row_low, indices);
            let high_scores = _mm_shuffle_epi8(row_high, indices);
            let high_mask = _mm_slli_epi16(_mm_and_si128(letters, _mm_set1_epi8(0x10)), 3);
            let match_scores = _mm_blendv_epi8(low_scores, high_scores, high_mask);

            score = _mm_adds_epi8(score, match_scores);
            best = _mm_max_epi8(best, score);
        }
    }

    _mm_storeu_si128(best_bytes.as_mut_ptr().cast(), best);
    for (dst, &value) in out.iter_mut().zip(&best_bytes[lane_offset..]) {
        *dst = i32::from(value) - i32::from(SCORE_BIAS);
    }
}

/// C++ `window_ungapped_best(...)`.
///
/// The widest upstream dispatch uses the scalar implementation for batches
/// smaller than four, avoiding SIMD saturation for those batches, and the
/// saturating lane implementation otherwise.
pub fn window_ungapped_best(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
) -> Vec<i32> {
    let mut out = vec![0; subjects.len()];
    window_ungapped_best_into(query, subjects, window, score_matrix, &mut out);
    out
}

/// Direct output-buffer form of C++ `window_ungapped_best(...)`.
pub fn window_ungapped_best_into(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    validate_buffers(query, subjects, window, out);
    if subjects.len() < 4 {
        for (score, subject) in out.iter_mut().zip(subjects) {
            *score = super::ungapped::ungapped_window(query, subject, window, score_matrix);
        }
    } else {
        window_ungapped_into(query, subjects, window, score_matrix, out);
    }
}

/// AVX2-specialized direct output form used after dispatch at a higher level.
///
/// # Safety
///
/// The caller must establish that AVX2 is available before calling this
/// function. Keeping that check at the stage boundary avoids repeating CPU
/// feature detection for every small subject batch.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
pub(crate) unsafe fn window_ungapped_best_into_avx2(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    debug_assert!(query.len() >= window);
    debug_assert!(subjects.iter().all(|subject| subject.len() >= window));
    debug_assert!(out.len() >= subjects.len());
    if subjects.len() < 4 {
        for (score, subject) in out.iter_mut().zip(subjects) {
            *score = super::ungapped::ungapped_window(query, subject, window, score_matrix);
        }
    } else {
        // SAFETY: inherited from this function's contract; debug builds also
        // verify the buffer invariants above.
        unsafe { window_ungapped_avx2(query, subjects, window, score_matrix, out) };
    }
}

/// SSE4.1-specialized direct output form used after dispatch at a higher
/// level.
///
/// # Safety
///
/// The caller must establish that SSSE3 and SSE4.1 are available before
/// calling this function.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
pub(crate) unsafe fn window_ungapped_best_into_sse41(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    debug_assert!(query.len() >= window);
    debug_assert!(subjects.iter().all(|subject| subject.len() >= window));
    debug_assert!(out.len() >= subjects.len());
    if subjects.len() < 4 {
        for (score, subject) in out.iter_mut().zip(subjects) {
            *score = super::ungapped::ungapped_window(query, subject, window, score_matrix);
        }
    } else {
        // SAFETY: inherited from this function's contract; debug builds also
        // verify the buffer invariants above.
        unsafe { window_ungapped_sse41(query, subjects, window, score_matrix, out) };
    }
}

/// Backwards-compatible name used by existing Rust stage-2 callers.
pub fn window_ungapped_multi(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
) -> Vec<i32> {
    window_ungapped_best(query, subjects, window, score_matrix)
}

fn validate_buffers(query: &[Letter], subjects: &[&[Letter]], window: usize, out: &[i32]) {
    assert!(
        query.len() >= window,
        "query is shorter than ungapped window"
    );
    assert!(
        subjects.iter().all(|subject| subject.len() >= window),
        "subject is shorter than ungapped window"
    );
    assert!(
        out.len() >= subjects.len(),
        "output is shorter than subject count"
    );
}

fn saturating_window_score(
    query: &[Letter],
    subject: &[Letter],
    window: usize,
    score_matrix: &ScoreMatrix,
) -> i32 {
    let mut score = SCORE_BIAS;
    let mut best = SCORE_BIAS;
    for position in 0..window {
        let query_letter = query[position] & LETTER_MASK;
        let subject_letter = subject[position] & LETTER_MASK;
        let match_score =
            score_matrix.matrix8()[(query_letter as usize) * 32 + subject_letter as usize];
        score = score.saturating_add(match_score);
        best = best.max(score);
    }
    i32::from(best) - i32::from(SCORE_BIAS)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    fn scalar_shifted_reference(
        query: &[Letter],
        subject: &[Letter],
        window: usize,
        score_matrix: &ScoreMatrix,
    ) -> i32 {
        let mut running = 0i32;
        let mut best = 0i32;
        for position in 0..window {
            running += score_matrix.score(query[position], subject[position]);
            running = running.clamp(0, 255);
            best = best.max(running);
        }
        best
    }

    #[test]
    fn lane_semantics_match_shifted_saturating_reference() {
        let score_matrix = matrix();
        let query: Vec<Letter> = (0..96).map(|index| (index % 20) as Letter).collect();
        let subject_data: Vec<Vec<Letter>> = (0..9)
            .map(|shift| {
                (0..96)
                    .map(|index| ((index + shift) % 20) as Letter)
                    .collect()
            })
            .collect();
        let subjects: Vec<&[Letter]> = subject_data.iter().map(Vec::as_slice).collect();
        let scores = window_ungapped(&query, &subjects, 96, &score_matrix);
        let expected: Vec<_> = subjects
            .iter()
            .map(|subject| scalar_shifted_reference(&query, subject, 96, &score_matrix))
            .collect();
        assert_eq!(scores, expected);
    }

    #[test]
    fn signed_saturation_resets_negative_runs_and_caps_at_255() {
        let score_matrix = matrix();
        let positive = vec![17; 80];
        let negative = vec![13; 80];
        assert_eq!(
            window_ungapped(&positive, &[&positive], 80, &score_matrix),
            vec![255]
        );

        let mut query = negative.clone();
        let subject = positive.clone();
        query[40..].fill(17);
        assert_eq!(
            window_ungapped(&query, &[&subject], 80, &score_matrix)[0],
            scalar_shifted_reference(&query, &subject, 80, &score_matrix)
        );
    }

    #[test]
    fn best_uses_unbounded_scalar_below_four_and_saturating_lanes_otherwise() {
        let score_matrix = matrix();
        let sequence = vec![17; 80];
        let three = [&sequence[..], &sequence[..], &sequence[..]];
        let scalar =
            super::super::ungapped::ungapped_window(&sequence, &sequence, 80, &score_matrix);
        assert!(scalar > 255);
        assert_eq!(
            window_ungapped_best(&sequence, &three, 80, &score_matrix),
            vec![scalar; 3]
        );

        let four = [&sequence[..], &sequence[..], &sequence[..], &sequence[..]];
        assert_eq!(
            window_ungapped_best(&sequence, &four, 80, &score_matrix),
            vec![255; 4]
        );
    }

    #[test]
    fn masked_letters_score_identically_to_unmasked_letters() {
        let score_matrix = matrix();
        let query = vec![5; 32];
        let subject = vec![5; 32];
        let masked_query = vec![5 | i8::MIN; 32];
        let masked_subject = vec![5 | i8::MIN; 32];
        assert_eq!(
            window_ungapped(&query, &[&subject], 32, &score_matrix),
            window_ungapped(&masked_query, &[&masked_subject], 32, &score_matrix)
        );
    }

    #[test]
    fn output_buffer_writes_only_subject_count_and_empty_window_is_zero() {
        let score_matrix = matrix();
        let sequence = vec![1; 4];
        let mut output = [91, 92, 93];
        window_ungapped_into(
            &sequence,
            &[&sequence, &sequence],
            0,
            &score_matrix,
            &mut output,
        );
        assert_eq!(output, [0, 0, 93]);
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn avx2_matches_scalar_for_all_batch_sizes_and_masked_letters() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let score_matrix = matrix();
        let query: Vec<Letter> = (0..113)
            .map(|i| ((i * 17 + 3) % 25) as Letter | if i % 7 == 0 { i8::MIN } else { 0 })
            .collect();
        let subject_data: Vec<Vec<Letter>> = (0..32)
            .map(|lane| {
                (0..113)
                    .map(|i| {
                        ((i * 11 + lane * 7 + 5) % 25) as Letter
                            | if (i + lane) % 9 == 0 { i8::MIN } else { 0 }
                    })
                    .collect()
            })
            .collect();
        for count in 1..=32 {
            let subjects: Vec<&[Letter]> =
                subject_data[..count].iter().map(Vec::as_slice).collect();
            let expected: Vec<_> = subjects
                .iter()
                .map(|subject| saturating_window_score(&query, subject, 113, &score_matrix))
                .collect();
            let mut actual = vec![0; count];
            unsafe {
                window_ungapped_avx2(&query, &subjects, 113, &score_matrix, &mut actual);
            }
            assert_eq!(actual, expected, "batch size {count}");
        }
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn sse41_matches_scalar_for_all_batch_sizes_and_masked_letters() {
        if !std::arch::is_x86_feature_detected!("ssse3")
            || !std::arch::is_x86_feature_detected!("sse4.1")
        {
            return;
        }
        let score_matrix = matrix();
        let query: Vec<Letter> = (0..113)
            .map(|i| ((i * 17 + 3) % 25) as Letter | if i % 7 == 0 { i8::MIN } else { 0 })
            .collect();
        let subject_data: Vec<Vec<Letter>> = (0..16)
            .map(|lane| {
                (0..113)
                    .map(|i| {
                        ((i * 11 + lane * 7 + 5) % 25) as Letter
                            | if (i + lane) % 9 == 0 { i8::MIN } else { 0 }
                    })
                    .collect()
            })
            .collect();
        for count in 1..=16 {
            let subjects: Vec<&[Letter]> =
                subject_data[..count].iter().map(Vec::as_slice).collect();
            let expected: Vec<_> = subjects
                .iter()
                .map(|subject| saturating_window_score(&query, subject, 113, &score_matrix))
                .collect();
            let mut actual = vec![0; count];
            unsafe {
                window_ungapped_sse41(&query, &subjects, 113, &score_matrix, &mut actual);
            }
            assert_eq!(actual, expected, "batch size {count}");
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn neon_matches_scalar_for_all_batch_sizes_and_masked_letters() {
        let score_matrix = matrix();
        let query: Vec<Letter> = (0..113)
            .map(|i| ((i * 17 + 3) % 25) as Letter | if i % 7 == 0 { i8::MIN } else { 0 })
            .collect();
        let subject_data: Vec<Vec<Letter>> = (0..16)
            .map(|lane| {
                (0..113)
                    .map(|i| {
                        ((i * 11 + lane * 7 + 5) % 25) as Letter
                            | if (i + lane) % 9 == 0 { i8::MIN } else { 0 }
                    })
                    .collect()
            })
            .collect();
        for count in 1..=16 {
            let subjects: Vec<&[Letter]> =
                subject_data[..count].iter().map(Vec::as_slice).collect();
            let expected: Vec<_> = subjects
                .iter()
                .map(|subject| saturating_window_score(&query, subject, 113, &score_matrix))
                .collect();
            let mut actual = vec![0; count];
            unsafe {
                window_ungapped_neon(&query, &subjects, 113, &score_matrix, &mut actual);
            }
            assert_eq!(actual, expected, "batch size {count}");
        }
    }

    #[test]
    #[should_panic(expected = "subject is shorter than ungapped window")]
    fn rejects_short_subject_instead_of_reading_out_of_bounds() {
        let score_matrix = matrix();
        window_ungapped(&[1; 4], &[&[1; 3]], 4, &score_matrix);
    }
}
