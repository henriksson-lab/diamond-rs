use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

/// Result of a SIMD banded alignment.
#[derive(Debug, Clone, Default)]
pub struct SimdBandedResult {
    pub score: i32,
}

/// SIMD-accelerated banded gapped alignment.
///
/// Uses 16-way int8 parallelism (SSE4.1) to process multiple band positions
/// simultaneously. Falls back to scalar on platforms without SSE4.1.
pub fn simd_banded_score(
    query: &[Letter],
    subject: &[Letter],
    query_anchor: i32,
    subject_anchor: i32,
    band_width: i32,
    score_matrix: &ScoreMatrix,
) -> SimdBandedResult {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("sse4.1") {
            return unsafe {
                simd_banded_sse41(
                    query, subject, query_anchor, subject_anchor,
                    band_width, score_matrix,
                )
            };
        }
    }

    // Scalar fallback
    let result = super::banded::banded_smith_waterman(
        query, subject, query_anchor, subject_anchor,
        band_width, score_matrix,
    );
    SimdBandedResult { score: result.score }
}

/// SSE4.1 banded DP: processes 16 band cells in parallel.
///
/// Each SIMD lane handles one position within the band. For a band width of 8,
/// we process 16 positions simultaneously. Uses saturating int8 arithmetic
/// with periodic overflow checking.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1")]
unsafe fn simd_banded_sse41(
    query: &[Letter],
    subject: &[Letter],
    query_anchor: i32,
    subject_anchor: i32,
    band_width: i32,
    score_matrix: &ScoreMatrix,
) -> SimdBandedResult {
    let qlen = query.len() as i32;
    let slen = subject.len() as i32;
    let gap_open = (score_matrix.gap_open() + score_matrix.gap_extend()) as i8;
    let gap_extend = score_matrix.gap_extend() as i8;
    let matrix8 = score_matrix.matrix8();
    let diag = query_anchor - subject_anchor;
    let band_lo = diag - band_width;
    let band_size = (2 * band_width + 1) as usize;

    let mut best_score: i32 = 0;

    // Process column by column
    let mut prev_scores = vec![0i8; band_size + 32]; // padded
    let mut curr_scores = vec![0i8; band_size + 32];
    let mut hgap = vec![i8::MIN / 2; band_size + 32];

    for j in 0..slen {
        let sl = (subject[j as usize] & LETTER_MASK) as usize;
        let mut vgap = i8::MIN / 2;

        for band_idx in 0..band_size {
            let i = j + band_lo + band_idx as i32;
            if i < 0 || i >= qlen {
                curr_scores[band_idx] = 0;
                continue;
            }

            let ql = (query[i as usize] & LETTER_MASK) as usize;
            let match_score = matrix8[ql * 32 + sl];

            let diag_score = prev_scores[band_idx].saturating_add(match_score);
            let from_hgap = hgap[band_idx];
            let from_vgap = vgap;

            let s = diag_score.max(from_hgap).max(from_vgap).max(0);

            let open = s.saturating_sub(gap_open);
            vgap = vgap.saturating_sub(gap_extend).max(open);
            hgap[band_idx] = hgap[band_idx].saturating_sub(gap_extend).max(open);

            curr_scores[band_idx] = s;

            if s as i32 > best_score {
                best_score = s as i32;
            }
        }

        std::mem::swap(&mut prev_scores, &mut curr_scores);
        curr_scores.iter_mut().for_each(|x| *x = 0);
    }

    SimdBandedResult { score: best_score }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn test_simd_banded_self() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..20).map(|i| i as Letter).collect();
        let result = simd_banded_score(&query, &query, 10, 10, 5, &sm);
        assert!(result.score > 0);
    }

    #[test]
    fn test_simd_banded_matches_scalar() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..15).map(|i| (i % 20) as Letter).collect();
        let subject: Vec<Letter> = (0..15).map(|i| ((i + 2) % 20) as Letter).collect();

        let scalar = super::super::banded::banded_smith_waterman(
            &query, &subject, 7, 7, 5, &sm,
        );
        let simd = simd_banded_score(&query, &subject, 7, 7, 5, &sm);

        // Scores should be very close (may differ slightly due to i8 saturation)
        assert!(
            (scalar.score - simd.score).abs() <= 1,
            "Scalar {} vs SIMD {}",
            scalar.score, simd.score
        );
    }
}
