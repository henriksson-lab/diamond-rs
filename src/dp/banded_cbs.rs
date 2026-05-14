//! Banded Smith-Waterman with CBS (composition-based statistics) and traceback.
//!
//! Scalar port of C++ `dp/swipe/banded_swipe.h`. The banding constrains the
//! alignment to stay near the seed diagonal, which is critical for correct
//! behaviour with CBS corrections — without banding, CBS penalties at masked
//! positions cause the full SW to introduce spurious gaps.
//!
//! The cell update matches C++ `swipe_cell_update` (cell_update.h):
//!   current = diag + match_score + cbs
//!   current = max(current, hgap, vgap, 0)  [local alignment]
//!   hgap = max(hgap - extend, current - open)
//!   vgap = max(vgap - extend, current - open)

use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{Letter, LETTER_MASK, SEED_MASK};
use crate::dp::smith_waterman::SwResult;
use crate::stats::score_matrix::ScoreMatrix;

/// Banded Smith-Waterman with CBS corrections and traceback.
///
/// Matches C++ banded_swipe.h scalar path:
/// - Band around the seed diagonal prevents spurious gaps
/// - CBS int8 corrections added per query position
/// - Soft-masked query positions (bit 7 set) score 0 (matching SIMD zeroing)
/// - Full traceback to produce alignment operations
pub fn banded_sw_cbs(
    query: &[Letter],
    subject: &[Letter],
    query_anchor: usize,
    subject_anchor: usize,
    band_width: i32,
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> SwResult {
    let qlen = query.len() as i32;
    let slen = subject.len() as i32;
    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();
    let neg_inf = i32::MIN / 4;

    let diag = query_anchor as i32 - subject_anchor as i32;
    let rows = qlen as usize + 1;
    let cols = slen as usize + 1;
    let cell_count = rows * cols;
    let idx = |i: usize, j: usize| -> usize { j * rows + i };

    // H is the local-alignment score. E and F are affine gap states ending
    // at the same cell after consuming subject or query characters.
    let mut h = vec![0i32; cell_count];
    let mut e = vec![neg_inf; cell_count];
    let mut f = vec![neg_inf; cell_count];

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for j in 1..=slen as usize {
        let subject_pos = j as i32 - 1;
        let lower = (subject_pos + diag - band_width).max(0);
        let upper = (subject_pos + diag + band_width).min(qlen - 1);
        if lower > upper {
            continue;
        }
        let i_begin = lower as usize + 1;
        let i_end = upper as usize + 1;

        for i in i_begin..=i_end {
            let query_pos = i - 1;
            let subject_pos = j - 1;

            // Match score: zero if query is soft-masked (matches C++ SIMD zeroing)
            let ql = query[query_pos];
            let sl = subject[subject_pos];
            let match_score = if ql & SEED_MASK != 0 {
                0
            } else {
                score_matrix.score(ql & LETTER_MASK, sl & LETTER_MASK)
            };

            // Cell update (matches C++ swipe_cell_update)
            let diag_score = h[idx(i - 1, j - 1)] + match_score + query_cbs[query_pos] as i32;
            let e_idx = idx(i, j);
            e[e_idx] = (h[idx(i, j - 1)] - gap_open).max(e[idx(i, j - 1)] - gap_extend);
            f[e_idx] = (h[idx(i - 1, j)] - gap_open).max(f[idx(i - 1, j)] - gap_extend);
            let s = diag_score.max(e[e_idx]).max(f[e_idx]).max(0);

            h[e_idx] = s;

            if s > best_score {
                best_score = s;
                best_i = i;
                best_j = j;
            }
        }
    }

    if best_score == 0 {
        return SwResult::default();
    }

    // Traceback from best cell
    let mut i = best_i;
    let mut j = best_j;
    let mut ops: Vec<(EditOperation, i32)> = Vec::new();
    let mut result = SwResult {
        score: best_score,
        query_end: i as i32 + 1,
        subject_end: j as i32 + 1,
        ..Default::default()
    };

    result.query_end = i as i32;
    result.subject_end = j as i32;

    while i > 0 && j > 0 && h[idx(i, j)] > 0 {
        let score = h[idx(i, j)];
        let ql = query[i - 1];
        let sl = subject[j - 1];
        let match_score = if ql & SEED_MASK != 0 {
            0
        } else {
            score_matrix.score(ql & LETTER_MASK, sl & LETTER_MASK)
        };
        let diag_score = h[idx(i - 1, j - 1)] + match_score + query_cbs[i - 1] as i32;

        if score == diag_score {
            // Diagonal move (match/substitution)
            if (ql & LETTER_MASK) == (sl & LETTER_MASK) {
                ops.push((EditOperation::Match, 1));
                result.identities += 1;
            } else {
                ops.push((EditOperation::Substitution, 1));
                result.mismatches += 1;
            }
            result.length += 1;
            i -= 1;
            j -= 1;
        } else if score == e[idx(i, j)] {
            // Gap in the query; consumes subject letters.
            let mut gap_len = 0i32;
            loop {
                gap_len += 1;
                let current = e[idx(i, j)];
                if current == h[idx(i, j - 1)] - gap_open {
                    j -= 1;
                    break;
                }
                j -= 1;
                if j == 0 || current != e[idx(i, j)] - gap_extend {
                    break;
                }
            }
            ops.push((EditOperation::Deletion, gap_len));
            result.gap_openings += 1;
            result.gaps += gap_len;
            result.length += gap_len;
        } else {
            // Gap in the subject; consumes query letters.
            let mut gap_len = 0i32;
            loop {
                gap_len += 1;
                let current = f[idx(i, j)];
                if current == h[idx(i - 1, j)] - gap_open {
                    i -= 1;
                    break;
                }
                i -= 1;
                if i == 0 || current != f[idx(i, j)] - gap_extend {
                    break;
                }
            }
            ops.push((EditOperation::Insertion, gap_len));
            result.gap_openings += 1;
            result.gaps += gap_len;
            result.length += gap_len;
        }
    }

    result.query_begin = i as i32;
    result.subject_begin = j as i32;
    ops.reverse();
    result.operations = ops;
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn test_banded_cbs_self_alignment() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..20).map(|i| i as Letter).collect();
        let cbs = vec![0i8; 20];
        let result = banded_sw_cbs(&query, &query, 10, 10, 15, &sm, &cbs);
        assert!(result.score > 0);
        assert_eq!(result.length, 20);
        assert_eq!(result.identities, 20);
        assert_eq!(result.gaps, 0);
        // Should match the full SW score
        let full = crate::dp::smith_waterman::smith_waterman(&query, &query, &sm);
        assert_eq!(
            result.score, full.score,
            "Banded CBS (zero cbs) should match full SW: {} vs {}",
            result.score, full.score
        );
    }

    #[test]
    fn test_banded_cbs_masked_positions() {
        let sm = make_test_matrix();
        // Create a sequence with some soft-masked positions
        let mut query: Vec<Letter> = (0..20).map(|i| i as Letter).collect();
        query[5] |= SEED_MASK; // mask position 5
        query[6] |= SEED_MASK; // mask position 6
        let cbs = vec![0i8; 20];
        let result = banded_sw_cbs(&query, &query, 10, 10, 15, &sm, &cbs);
        // Masked positions score 0, not the diagonal score
        // So total score should be less than without masking
        let unmasked_query: Vec<Letter> = (0..20).map(|i| i as Letter).collect();
        let unmasked = banded_sw_cbs(&unmasked_query, &unmasked_query, 10, 10, 15, &sm, &cbs);
        assert!(
            result.score < unmasked.score,
            "Masked should score less: {} vs {}",
            result.score,
            unmasked.score
        );
        // But alignment should still be 100% identity (banding prevents gaps)
        assert_eq!(result.identities, 20);
        assert_eq!(result.gaps, 0);
    }

    #[test]
    fn test_banded_cbs_uses_bias() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..8).map(|i| i as Letter).collect();
        let cbs = vec![1i8; 8];
        let biased = banded_sw_cbs(&query, &query, 4, 4, 8, &sm, &cbs);
        let unbiased = banded_sw_cbs(&query, &query, 4, 4, 8, &sm, &[0i8; 8]);
        assert_eq!(biased.score, unbiased.score + 8);
        assert_eq!(biased.identities, 8);
    }
}
