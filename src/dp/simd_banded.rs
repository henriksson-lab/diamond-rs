use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

/// Result of a SIMD banded alignment.
#[derive(Debug, Clone, Default)]
pub struct SimdBandedResult {
    pub score: i32,
}

/// SIMD-accelerated banded gapped alignment.
///
/// Uses an exact anti-diagonal i32 vector recurrence when SSE4.1 is available.
/// Cells on the same anti-diagonal are independent for the affine-gap
/// recurrence, so this can process four band cells at a time without changing
/// scalar Smith-Waterman semantics. Falls back to scalar on other platforms.
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
        if band_width >= 0 && is_x86_feature_detected!("sse4.1") {
            return unsafe {
                simd_banded_sse41_i32(
                    query,
                    subject,
                    query_anchor,
                    subject_anchor,
                    band_width,
                    score_matrix,
                )
            };
        }
    }

    let result = super::banded::banded_smith_waterman(
        query,
        subject,
        query_anchor,
        subject_anchor,
        band_width,
        score_matrix,
    );
    SimdBandedResult {
        score: result.score,
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1")]
unsafe fn simd_banded_sse41_i32(
    query: &[Letter],
    subject: &[Letter],
    query_anchor: i32,
    subject_anchor: i32,
    band_width: i32,
    score_matrix: &ScoreMatrix,
) -> SimdBandedResult {
    use std::arch::x86_64::*;

    let qlen = query.len() as i32;
    let slen = subject.len() as i32;
    if qlen == 0 || slen == 0 {
        return SimdBandedResult::default();
    }

    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();
    let neg_inf = i32::MIN / 4;
    let diag = query_anchor - subject_anchor;
    let rows = qlen as usize + 1;
    let cols = slen as usize + 1;
    let matrix32 = score_matrix.matrix32();
    let idx = |i: usize, j: usize| -> usize { j * rows + i };
    let in_band = |i: i32, d: i32| -> bool { (2 * i - d - diag).abs() <= band_width };

    let mut h = vec![0i32; rows * cols];
    let mut e = vec![neg_inf; rows * cols];
    let mut f = vec![neg_inf; rows * cols];
    let mut best_score = 0i32;

    let v_gap_open = _mm_set1_epi32(gap_open);
    let v_gap_extend = _mm_set1_epi32(gap_extend);
    let v_zero = _mm_setzero_si128();

    for d in 2..=(qlen + slen) {
        let mut i_start = 1.max(d - slen);
        let mut i_end = qlen.min(d - 1);
        while i_start <= i_end && !in_band(i_start, d) {
            i_start += 1;
        }
        while i_end >= i_start && !in_band(i_end, d) {
            i_end -= 1;
        }
        let mut i = i_start;
        while i <= i_end {
            let mut h_diag_arr = [0i32; 4];
            let mut h_left_arr = [0i32; 4];
            let mut e_left_arr = [neg_inf; 4];
            let mut h_up_arr = [0i32; 4];
            let mut f_up_arr = [neg_inf; 4];
            let mut score_arr = [0i32; 4];
            let mut current_idx = [0usize; 4];
            let mut lane_count = 0usize;

            for lane in 0..4 {
                let ii = i + lane as i32;
                if ii > i_end {
                    break;
                }
                let jj = d - ii;
                let iu = ii as usize;
                let ju = jj as usize;
                let current = idx(iu, ju);
                current_idx[lane] = current;
                h_diag_arr[lane] = h[idx(iu - 1, ju - 1)];
                h_left_arr[lane] = h[idx(iu, ju - 1)];
                e_left_arr[lane] = e[idx(iu, ju - 1)];
                h_up_arr[lane] = h[idx(iu - 1, ju)];
                f_up_arr[lane] = f[idx(iu - 1, ju)];

                let ql = (query[iu - 1] & LETTER_MASK) as usize;
                let sl = (subject[ju - 1] & LETTER_MASK) as usize;
                score_arr[lane] = matrix32[ql * 32 + sl];
                lane_count += 1;
            }

            let h_diag_v = _mm_loadu_si128(h_diag_arr.as_ptr() as *const __m128i);
            let h_left_v = _mm_loadu_si128(h_left_arr.as_ptr() as *const __m128i);
            let e_left_v = _mm_loadu_si128(e_left_arr.as_ptr() as *const __m128i);
            let h_up_v = _mm_loadu_si128(h_up_arr.as_ptr() as *const __m128i);
            let f_up_v = _mm_loadu_si128(f_up_arr.as_ptr() as *const __m128i);
            let score_v = _mm_loadu_si128(score_arr.as_ptr() as *const __m128i);

            let diag_score = _mm_add_epi32(h_diag_v, score_v);
            let e_open = _mm_sub_epi32(h_left_v, v_gap_open);
            let e_extend = _mm_sub_epi32(e_left_v, v_gap_extend);
            let e_v = _mm_max_epi32(e_open, e_extend);
            let f_open = _mm_sub_epi32(h_up_v, v_gap_open);
            let f_extend = _mm_sub_epi32(f_up_v, v_gap_extend);
            let f_v = _mm_max_epi32(f_open, f_extend);
            let h_v = _mm_max_epi32(_mm_max_epi32(diag_score, e_v), _mm_max_epi32(f_v, v_zero));

            let mut e_out = [0i32; 4];
            let mut f_out = [0i32; 4];
            let mut h_out = [0i32; 4];
            _mm_storeu_si128(e_out.as_mut_ptr() as *mut __m128i, e_v);
            _mm_storeu_si128(f_out.as_mut_ptr() as *mut __m128i, f_v);
            _mm_storeu_si128(h_out.as_mut_ptr() as *mut __m128i, h_v);

            for lane in 0..lane_count {
                let current = current_idx[lane];
                e[current] = e_out[lane];
                f[current] = f_out[lane];
                h[current] = h_out[lane];
                best_score = best_score.max(h_out[lane]);
            }

            i += 4;
        }
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

        let scalar = super::super::banded::banded_smith_waterman(&query, &subject, 7, 7, 5, &sm);
        let simd = simd_banded_score(&query, &subject, 7, 7, 5, &sm);

        assert_eq!(simd.score, scalar.score);
    }

    #[test]
    fn test_simd_banded_does_not_saturate_high_scores() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = vec![17; 40];
        let scalar = super::super::banded::banded_smith_waterman(&query, &query, 20, 20, 5, &sm);
        let simd = simd_banded_score(&query, &query, 20, 20, 5, &sm);
        assert_eq!(simd.score, scalar.score);
        assert!(simd.score > i8::MAX as i32);
    }

    #[test]
    fn test_simd_banded_matches_scalar_across_offsets() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..37).map(|i| ((i * 7) % 20) as Letter).collect();
        let subject: Vec<Letter> = (0..31).map(|i| ((i * 5 + 3) % 20) as Letter).collect();
        for query_anchor in [0, 5, 13, 24, 36] {
            for subject_anchor in [0, 4, 11, 20, 30] {
                for band_width in [0, 1, 3, 8, 20] {
                    let scalar = super::super::banded::banded_smith_waterman(
                        &query,
                        &subject,
                        query_anchor,
                        subject_anchor,
                        band_width,
                        &sm,
                    );
                    let simd = simd_banded_score(
                        &query,
                        &subject,
                        query_anchor,
                        subject_anchor,
                        band_width,
                        &sm,
                    );
                    assert_eq!(
                        simd.score, scalar.score,
                        "anchor {query_anchor}/{subject_anchor} band {band_width}"
                    );
                }
            }
        }
    }
}
