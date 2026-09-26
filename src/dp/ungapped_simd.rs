//! Translation of `diamond/src/dp/ungapped_simd.{h,cpp}`.
//!
//! Upstream uses SIMD only as a lane-parallel implementation detail. Each
//! lane is a signed saturating i8 accumulator biased by `SCHAR_MIN`; this
//! architecture-neutral translation preserves those exact lane semantics on
//! every target.

use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

const SCORE_BIAS: i8 = i8::MIN;

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
    for (score, subject) in out.iter_mut().zip(subjects) {
        *score = saturating_window_score(query, subject, window, score_matrix);
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

    #[test]
    #[should_panic(expected = "subject is shorter than ungapped window")]
    fn rejects_short_subject_instead_of_reading_out_of_bounds() {
        let score_matrix = matrix();
        window_ungapped(&[1; 4], &[&[1; 3]], 4, &score_matrix);
    }
}
