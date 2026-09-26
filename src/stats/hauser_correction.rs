//! Translation of `diamond/src/stats/hauser_correction.{h,cpp}`.

use std::ops::{Deref, DerefMut};

use crate::align::hsp::Hsp;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{Letter, AMINO_ACID_COUNT, LETTER_MASK, TRUE_AA};
use crate::dp::ungapped::DiagonalSegment;
use crate::stats::score_matrix::ScoreMatrix;

pub const PADDING: usize = 32;

struct VectorScores {
    scores: [i32; TRUE_AA as usize],
}

impl VectorScores {
    fn new() -> Self {
        Self {
            scores: [0; TRUE_AA as usize],
        }
    }

    fn add(&mut self, letter: Letter, score_matrix: &ScoreMatrix) {
        let letter = letter & LETTER_MASK;
        for i in 0..TRUE_AA as usize {
            self.scores[i] += score_matrix.score(letter, i as Letter);
        }
    }

    fn sub(&mut self, letter: Letter, score_matrix: &ScoreMatrix) {
        let letter = letter & LETTER_MASK;
        for i in 0..TRUE_AA as usize {
            self.scores[i] -= score_matrix.score(letter, i as Letter);
        }
    }
}

/// C++ `HauserCorrection`, including its floating-point vector base and the
/// padded signed-byte representation used by SIMD score profiles.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HauserCorrection {
    pub values: Vec<f32>,
    pub int8: Vec<i8>,
}

impl HauserCorrection {
    /// C++ `HauserCorrection(const Sequence&)`, with the global score matrix
    /// and `config.cbs_window` passed explicitly.
    pub fn new(seq: &[Letter], score_matrix: &ScoreMatrix, window: usize) -> Self {
        let len = seq.len();
        let background_scores = score_matrix.background_scores();
        let mut values = vec![0.0f32; len];
        let mut scores = VectorScores::new();
        let window_half = (window / 2).min(len.wrapping_sub(1));
        let mut n = 0usize;
        let mut h = 0usize;
        let mut m = 0usize;
        let mut t = 0usize;

        while n < window_half && h < len {
            n += 1;
            scores.add(seq[h], score_matrix);
            h += 1;
        }
        while n < window + 1 && h < len {
            n += 1;
            scores.add(seq[h], score_matrix);
            set_correction(
                &mut values,
                m,
                seq[m],
                &scores,
                n,
                background_scores,
                score_matrix,
            );
            h += 1;
            m += 1;
        }
        while h < len {
            scores.add(seq[h], score_matrix);
            scores.sub(seq[t], score_matrix);
            set_correction(
                &mut values,
                m,
                seq[m],
                &scores,
                n,
                background_scores,
                score_matrix,
            );
            h += 1;
            t += 1;
            m += 1;
        }
        while m < len && n > window_half + 1 {
            n -= 1;
            scores.sub(seq[t], score_matrix);
            set_correction(
                &mut values,
                m,
                seq[m],
                &scores,
                n,
                background_scores,
                score_matrix,
            );
            t += 1;
            m += 1;
        }
        while m < len {
            set_correction(
                &mut values,
                m,
                seq[m],
                &scores,
                n,
                background_scores,
                score_matrix,
            );
            m += 1;
        }

        let mut int8 = Vec::with_capacity(len + PADDING);
        int8.extend(values.iter().map(|&value| {
            let rounded = if value < 0.0 {
                value - 0.5
            } else {
                value + 0.5
            };
            rounded as i8
        }));
        int8.resize(len + PADDING, 0);
        Self { values, int8 }
    }

    /// Header-inline `HauserCorrection::operator()(float&, int, int, int)`.
    pub fn apply(&self, score: &mut f32, i: i32, query_anchor: i32, mult: i32) {
        let position = query_anchor + i * mult;
        *score += self.values[position as usize];
    }

    /// C++ `HauserCorrection::operator()(const Hsp&)`.
    pub fn score_hsp(&self, hsp: &Hsp) -> i32 {
        let mut score = 0.0f32;
        let mut iterator = hsp.begin();
        while iterator.good() {
            if matches!(
                iterator.op(),
                EditOperation::Match | EditOperation::Substitution
            ) {
                score += self.values[iterator.query_pos.translated as usize];
            }
            iterator.advance();
        }
        score as i32
    }

    /// C++ `HauserCorrection::operator()(const DiagonalSegment&)`.
    pub fn score_diagonal(&self, diagonal: &DiagonalSegment) -> i32 {
        let mut score = 0.0f32;
        for position in diagonal.i..diagonal.query_end() {
            score += self.values[position as usize];
        }
        score as i32
    }

    /// C++ `HauserCorrection::reverse(const int8_t*, size_t)`. `None`
    /// represents the null-pointer branch.
    pub fn reverse(values: Option<&[i8]>, len: usize) -> Vec<i8> {
        let Some(values) = values else {
            return Vec::new();
        };
        assert!(
            len <= values.len(),
            "Hauser correction reverse length exceeds input"
        );
        values[..len].iter().rev().copied().collect()
    }
}

impl Deref for HauserCorrection {
    type Target = Vec<f32>;

    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

impl DerefMut for HauserCorrection {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.values
    }
}

#[allow(clippy::too_many_arguments)]
fn set_correction(
    values: &mut [f32],
    position: usize,
    letter: Letter,
    scores: &VectorScores,
    n: usize,
    background_scores: &[f64; TRUE_AA as usize],
    score_matrix: &ScoreMatrix,
) {
    let letter = (letter & LETTER_MASK) as usize;
    if letter < TRUE_AA as usize {
        values[position] = background_scores[letter] as f32
            - (scores.scores[letter] - score_matrix.score(letter as Letter, letter as Letter))
                as f32
                / (n - 1) as f32;
    }
}

/// Compatibility helper previously exposed from `stats::cbs`.
pub fn compute_background_scores(score_matrix: &ScoreMatrix) -> [f64; TRUE_AA as usize] {
    *score_matrix.background_scores()
}

/// C++ `Stats::hauser_global`.
pub fn hauser_global(
    query_comp: &[f64; TRUE_AA as usize],
    target_comp: &[f64; TRUE_AA as usize],
    score_matrix: &ScoreMatrix,
) -> Vec<i32> {
    let background_scores = score_matrix.background_scores();
    let mut query_scores = [0.0f64; TRUE_AA as usize];
    let mut target_scores = [0.0f64; TRUE_AA as usize];
    for i in 0..TRUE_AA as usize {
        for j in 0..TRUE_AA as usize {
            query_scores[i] +=
                query_comp[j] * f64::from(score_matrix.score(i as Letter, j as Letter));
            target_scores[i] +=
                target_comp[j] * f64::from(score_matrix.score(i as Letter, j as Letter));
        }
    }
    for i in 0..TRUE_AA as usize {
        query_scores[i] = background_scores[i] - query_scores[i];
        target_scores[i] = background_scores[i] - target_scores[i];
    }

    let mut matrix = vec![0; AMINO_ACID_COUNT * AMINO_ACID_COUNT];
    for i in 0..AMINO_ACID_COUNT {
        for j in 0..AMINO_ACID_COUNT {
            let score = f64::from(score_matrix.score(i as Letter, j as Letter));
            let query = if i < TRUE_AA as usize {
                query_scores[i]
            } else {
                0.0
            };
            let target = if j < TRUE_AA as usize {
                target_scores[j]
            } else {
                0.0
            };
            matrix[i * AMINO_ACID_COUNT + j] = (score + query.min(target)).round() as i32;
        }
    }
    matrix
}

/// Compatibility function for the historical Rust API.
pub fn hauser_correction(seq: &[Letter], score_matrix: &ScoreMatrix) -> Vec<i8> {
    hauser_correction_window(seq, score_matrix, 40)
}

/// Compatibility function with explicit C++ `config.cbs_window`.
pub fn hauser_correction_window(
    seq: &[Letter],
    score_matrix: &ScoreMatrix,
    window: usize,
) -> Vec<i8> {
    HauserCorrection::new(seq, score_matrix, window).int8
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::EditOperation;
    use crate::util::interval::Interval;

    fn matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn constructor_matches_independent_sliding_window_reference() {
        let score_matrix = matrix();
        let seq: Vec<Letter> = (0..73).map(|i| (i % 20) as Letter).collect();
        let correction = HauserCorrection::new(&seq, &score_matrix, 40);
        assert_eq!(correction.values.len(), seq.len());
        assert_eq!(correction.int8.len(), seq.len() + PADDING);
        assert!(correction.int8[seq.len()..].iter().all(|&value| value == 0));

        for position in [0usize, 19, 20, 39, 52, 72] {
            let begin = position.saturating_sub(20);
            let end = (position + 21).min(seq.len());
            let letter = seq[position] as usize;
            let neighbor_sum: i32 = seq[begin..end]
                .iter()
                .enumerate()
                .filter(|(offset, _)| begin + offset != position)
                .map(|(_, &other)| score_matrix.score(letter as Letter, other))
                .sum();
            let expected = score_matrix.background_scores()[letter] as f32
                - neighbor_sum as f32 / (end - begin - 1) as f32;
            assert!((correction.values[position] - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn masked_and_non_true_letters_follow_sequence_and_zero_rules() {
        let score_matrix = matrix();
        let plain = HauserCorrection::new(&[0, 1, 2, 3], &score_matrix, 2);
        let masked = HauserCorrection::new(&[-128, -127, -126, -125], &score_matrix, 2);
        assert_eq!(plain, masked);
        let unknown = HauserCorrection::new(&[23, 23, 23], &score_matrix, 2);
        assert_eq!(unknown.values, [0.0; 3]);
    }

    #[test]
    fn empty_sequence_still_has_cpp_simd_padding() {
        let correction = HauserCorrection::new(&[], &matrix(), 40);
        assert!(correction.values.is_empty());
        assert_eq!(correction.int8, vec![0; PADDING]);
        assert!(HauserCorrection::default().int8.is_empty());
    }

    #[test]
    fn apply_hsp_and_diagonal_score_only_query_consuming_matches() {
        let correction = HauserCorrection {
            values: vec![0.0, 1.25, 2.5, 4.0, 8.0, 16.0],
            int8: Vec::new(),
        };
        let mut score = 3.0;
        correction.apply(&mut score, 2, 5, -1);
        assert_eq!(score, 7.0);

        let mut hsp = Hsp::new();
        hsp.query_range = Interval::new(1, 5);
        hsp.subject_range = Interval::new(0, 3);
        hsp.transcript.push_with_count(EditOperation::Match, 2);
        hsp.transcript.push_with_count(EditOperation::Insertion, 1);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 3);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 4);
        hsp.transcript.push_terminator();
        assert_eq!(correction.score_hsp(&hsp), (1.25f32 + 2.5 + 8.0) as i32);

        let diagonal = DiagonalSegment::new(1, 0, 3, 0);
        assert_eq!(
            correction.score_diagonal(&diagonal),
            (1.25f32 + 2.5 + 4.0) as i32
        );
    }

    #[test]
    fn reverse_handles_null_prefix_and_exact_length() {
        assert!(HauserCorrection::reverse(None, 99).is_empty());
        assert_eq!(HauserCorrection::reverse(Some(&[1, 2, 3, 4]), 3), [3, 2, 1]);
    }

    #[test]
    fn global_matrix_matches_direct_formula_and_extended_cells() {
        let score_matrix = matrix();
        let mut query = [0.0; TRUE_AA as usize];
        let mut target = [0.0; TRUE_AA as usize];
        query[0] = 1.0;
        target[1] = 1.0;
        let adjusted = hauser_global(&query, &target, &score_matrix);
        assert_eq!(adjusted.len(), AMINO_ACID_COUNT * AMINO_ACID_COUNT);
        let i = 2usize;
        let j = 3usize;
        let query_delta =
            score_matrix.background_scores()[i] - f64::from(score_matrix.score(i as Letter, 0));
        let target_delta =
            score_matrix.background_scores()[j] - f64::from(score_matrix.score(j as Letter, 1));
        assert_eq!(
            adjusted[i * AMINO_ACID_COUNT + j],
            (f64::from(score_matrix.score(i as Letter, j as Letter))
                + query_delta.min(target_delta))
            .round() as i32
        );
        let x = AMINO_ACID_COUNT - 1;
        assert_eq!(
            adjusted[x * AMINO_ACID_COUNT + x],
            score_matrix.score(x as Letter, x as Letter)
        );
    }
}
