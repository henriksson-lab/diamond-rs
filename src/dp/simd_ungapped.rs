use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

/// SIMD-accelerated ungapped window scoring for multiple subjects.
///
/// Scores `subject_count` subjects against the same query window simultaneously.
/// Falls back to scalar on platforms without SSE4.1.
pub fn window_ungapped_multi(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
) -> Vec<i32> {
    let subject_count = subjects.len();
    let mut out = vec![0i32; subject_count];

    // Try SIMD on x86
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("sse4.1") && subject_count >= 4 {
            unsafe {
                window_ungapped_sse41(query, subjects, window, score_matrix, &mut out);
            }
            return out;
        }
    }

    // Scalar fallback
    for (i, subject) in subjects.iter().enumerate() {
        out[i] = super::ungapped::ungapped_window(query, subject, window, score_matrix);
    }
    out
}

/// SSE4.1 implementation of multi-subject ungapped scoring.
///
/// Processes 16 subjects simultaneously using 128-bit SIMD vectors of int8.
/// Note: Uses saturating i8 arithmetic (clips at 127), matching C++ DIAMOND behavior.
/// This is intentional — SIMD ungapped scoring is a fast pre-filter; accurate scoring
/// uses wider types in the gapped alignment phase.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1")]
unsafe fn window_ungapped_sse41(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    use std::arch::x86_64::*;

    let subject_count = subjects.len().min(16);
    let matrix8 = score_matrix.matrix8();

    // Initialize score and best vectors to zero
    let mut score = _mm_setzero_si128();
    let mut best = _mm_setzero_si128();
    let zero = _mm_setzero_si128();

    for pos in 0..window.min(query.len()) {
        let ql = (query[pos] & LETTER_MASK) as usize;

        // Build a vector of subject letters at this position
        let mut subject_letters = [0i8; 16];
        for (i, subj) in subjects[..subject_count].iter().enumerate() {
            if pos < subj.len() {
                subject_letters[16 - subject_count + i] =
                    subj[pos] & LETTER_MASK;
            }
        }

        // Look up scores: for each subject letter, get score[query_letter][subject_letter]
        // We need to use the matrix row for the query letter
        let row_offset = ql * 32;
        let mut scores_arr = [0i8; 16];
        for i in 0..16 {
            let sl = subject_letters[i] as usize;
            if sl < 32 {
                scores_arr[i] = matrix8[row_offset + sl];
            }
        }

        let match_scores = _mm_loadu_si128(scores_arr.as_ptr() as *const __m128i);

        // score = max(score + match_scores, 0)
        score = _mm_adds_epi8(score, match_scores);
        score = _mm_max_epi8(score, zero);

        // best = max(best, score)
        best = _mm_max_epi8(best, score);
    }

    // Extract results
    let mut best_arr = [0i8; 16];
    _mm_storeu_si128(best_arr.as_mut_ptr() as *mut __m128i, best);

    let offset = 16 - subject_count;
    for i in 0..subject_count {
        out[i] = best_arr[offset + i] as i32;
    }
}

/// AVX2 implementation of multi-subject ungapped scoring.
///
/// Processes 32 subjects simultaneously using 256-bit SIMD vectors.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(dead_code)] // Will be used when subject_count > 16
unsafe fn window_ungapped_avx2(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    use std::arch::x86_64::*;

    let subject_count = subjects.len().min(32);
    let matrix8 = score_matrix.matrix8();

    let mut score = _mm256_setzero_si256();
    let mut best = _mm256_setzero_si256();
    let zero = _mm256_setzero_si256();

    for pos in 0..window.min(query.len()) {
        let ql = (query[pos] & LETTER_MASK) as usize;
        let row_offset = ql * 32;

        let mut scores_arr = [0i8; 32];
        for (i, subj) in subjects[..subject_count].iter().enumerate() {
            if pos < subj.len() {
                let sl = (subj[pos] & LETTER_MASK) as usize;
                if sl < 32 {
                    scores_arr[32 - subject_count + i] = matrix8[row_offset + sl];
                }
            }
        }

        let match_scores = _mm256_loadu_si256(scores_arr.as_ptr() as *const __m256i);

        score = _mm256_adds_epi8(score, match_scores);
        score = _mm256_max_epi8(score, zero);
        best = _mm256_max_epi8(best, score);
    }

    let mut best_arr = [0i8; 32];
    _mm256_storeu_si256(best_arr.as_mut_ptr() as *mut __m256i, best);

    let offset = 32 - subject_count;
    for i in 0..subject_count {
        out[i] = best_arr[offset + i] as i32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn test_multi_ungapped_self() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..20).map(|i| i as Letter).collect();
        let subjects: Vec<&[Letter]> = vec![&query, &query];
        let scores = window_ungapped_multi(&query, &subjects, 20, &sm);
        assert_eq!(scores.len(), 2);
        // Both subjects are identical to query, so scores should be equal and positive
        assert!(scores[0] > 0);
        assert_eq!(scores[0], scores[1]);
    }

    #[test]
    fn test_multi_ungapped_different() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = vec![0; 20]; // All A
        let good_subject: Vec<Letter> = vec![0; 20]; // All A (matches)
        let bad_subject: Vec<Letter> = vec![13; 20]; // All F (mismatches)
        let subjects: Vec<&[Letter]> = vec![&good_subject, &bad_subject];
        let scores = window_ungapped_multi(&query, &subjects, 20, &sm);
        assert!(scores[0] > scores[1], "Matching subject should score higher");
    }

    #[test]
    fn test_multi_ungapped_many_subjects() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..10).map(|i| i as Letter).collect();
        // Create 8 subjects
        let subject_data: Vec<Vec<Letter>> = (0..8)
            .map(|j| (0..10).map(|i| ((i + j) % 20) as Letter).collect())
            .collect();
        let subjects: Vec<&[Letter]> = subject_data.iter().map(|s| s.as_slice()).collect();
        let scores = window_ungapped_multi(&query, &subjects, 10, &sm);
        assert_eq!(scores.len(), 8);
        // First subject should have highest score (identical)
        assert!(scores[0] >= scores[1]);
    }

    #[test]
    fn test_scalar_matches_simd() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..15).map(|i| (i % 20) as Letter).collect();
        let subject1: Vec<Letter> = (0..15).map(|i| (i % 20) as Letter).collect();
        let subject2: Vec<Letter> = (0..15).map(|i| ((i + 3) % 20) as Letter).collect();

        // Scalar
        let s1 = super::super::ungapped::ungapped_window(&query, &subject1, 15, &sm);
        let s2 = super::super::ungapped::ungapped_window(&query, &subject2, 15, &sm);

        // Multi (may use SIMD)
        let subjects: Vec<&[Letter]> = vec![&subject1, &subject2];
        let multi = window_ungapped_multi(&query, &subjects, 15, &sm);

        assert_eq!(multi[0], s1, "SIMD result should match scalar for subject 1");
        assert_eq!(multi[1], s2, "SIMD result should match scalar for subject 2");
    }
}
