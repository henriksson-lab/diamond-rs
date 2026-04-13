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

    let diag = query_anchor as i32 - subject_anchor as i32;
    let band_size = (2 * band_width + 1) as usize;

    // DP matrix for traceback: dp[j][band_idx]
    // Store full matrix so we can trace back
    let cols = slen as usize;
    let mut dp = vec![vec![0i32; band_size + 2]; cols + 1];
    let mut hgap_store = vec![vec![i32::MIN / 2; band_size + 2]; cols + 1];

    let mut best_score = 0i32;
    let mut best_j = 0i32;
    let mut best_band_idx = 0usize;

    for j in 0..slen {
        let ju = j as usize;
        let mut vgap = i32::MIN / 2;

        for band_idx in 0..band_size {
            let i = j + (diag - band_width) + band_idx as i32;
            if i < 0 || i >= qlen {
                dp[ju + 1][band_idx] = 0;
                continue;
            }
            let iu = i as usize;

            // Match score: zero if query is soft-masked (matches C++ SIMD zeroing)
            let ql = query[iu];
            let sl = subject[ju];
            let match_score = if ql & SEED_MASK != 0 {
                0
            } else {
                score_matrix.score(ql & LETTER_MASK, sl & LETTER_MASK)
            };

            // CBS correction: zero at masked positions
            let cbs = if ql & SEED_MASK != 0 { 0 } else { query_cbs[iu] as i32 };

            // Cell update (matches C++ swipe_cell_update)
            let diag_score = dp[ju][band_idx] + match_score + cbs;
            let from_hgap = hgap_store[ju + 1][band_idx];
            let from_vgap = vgap;

            let s = diag_score.max(from_hgap).max(from_vgap).max(0);

            let open = s - gap_open;
            vgap = (vgap - gap_extend).max(open);
            hgap_store[ju + 1][band_idx + 1] = (hgap_store[ju + 1][band_idx + 1].max(i32::MIN / 2) - gap_extend).max(open);

            dp[ju + 1][band_idx] = s;

            if s > best_score {
                best_score = s;
                best_j = j;
                best_band_idx = band_idx;
            }
        }
    }

    if best_score == 0 {
        return SwResult::default();
    }

    // Traceback from best cell
    let mut j = best_j as usize;
    let mut band_idx = best_band_idx;
    let mut i = (best_j + (diag - band_width) + band_idx as i32) as usize;
    let mut ops: Vec<(EditOperation, i32)> = Vec::new();
    let mut result = SwResult {
        score: best_score,
        query_end: i as i32 + 1,
        subject_end: j as i32 + 1,
        ..Default::default()
    };

    while dp[j + 1][band_idx] > 0 && j < cols && i < query.len() as usize {
        let score = dp[j + 1][band_idx];
        let ql = query[i];
        let sl = subject[j];
        let match_score = if ql & SEED_MASK != 0 {
            0
        } else {
            score_matrix.score(ql & LETTER_MASK, sl & LETTER_MASK)
        };
        let cbs = if ql & SEED_MASK != 0 { 0 } else { query_cbs[i] as i32 };
        let diag_score = dp[j][band_idx] + match_score + cbs;

        if score == diag_score && dp[j][band_idx] >= 0 {
            // Diagonal move (match/substitution)
            if (ql & LETTER_MASK) == (sl & LETTER_MASK) {
                ops.push((EditOperation::Match, 1));
                result.identities += 1;
            } else {
                ops.push((EditOperation::Substitution, 1));
                result.mismatches += 1;
            }
            result.length += 1;
            if i == 0 || j == 0 { break; }
            i -= 1;
            j -= 1;
            // band_idx stays the same on diagonal move
        } else {
            // Gap — check if horizontal or vertical
            // Horizontal gap: move in j (insertion in subject)
            // In banded coords: band_idx increases by 1 when j decreases
            let hgap_score = if band_idx + 1 < band_size + 2 {
                hgap_store[j + 1][band_idx]
            } else {
                i32::MIN / 2
            };

            if score == hgap_score {
                ops.push((EditOperation::Insertion, 1));
                result.gap_openings += 1;
                result.gaps += 1;
                result.length += 1;
                if j == 0 { break; }
                j -= 1;
                band_idx += 1;
                if band_idx >= band_size { break; }
            } else {
                // Vertical gap (deletion in subject)
                ops.push((EditOperation::Deletion, 1));
                result.gap_openings += 1;
                result.gaps += 1;
                result.length += 1;
                if i == 0 { break; }
                i -= 1;
                if band_idx == 0 { break; }
                band_idx -= 1;
            }
        }
    }

    result.query_begin = i as i32;
    result.subject_begin = j as i32;
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
        assert_eq!(result.score, full.score,
            "Banded CBS (zero cbs) should match full SW: {} vs {}", result.score, full.score);
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
        assert!(result.score < unmasked.score,
            "Masked should score less: {} vs {}", result.score, unmasked.score);
        // But alignment should still be 100% identity (banding prevents gaps)
        assert_eq!(result.identities, 20);
        assert_eq!(result.gaps, 0);
    }
}
