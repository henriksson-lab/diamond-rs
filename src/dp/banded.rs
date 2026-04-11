use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

/// Result of a banded gapped alignment.
#[derive(Debug, Clone, Default)]
pub struct BandedResult {
    pub score: i32,
    pub query_begin: i32,
    pub query_end: i32,
    pub subject_begin: i32,
    pub subject_end: i32,
}

/// Banded Smith-Waterman alignment.
///
/// Performs gapped alignment within a diagonal band around a seed hit.
/// The band is defined by [diag - band_width, diag + band_width].
///
/// This is the scalar implementation; SIMD versions process multiple
/// cells in parallel using score vectors.
pub fn banded_smith_waterman(
    query: &[Letter],
    subject: &[Letter],
    query_anchor: i32,
    subject_anchor: i32,
    band_width: i32,
    score_matrix: &ScoreMatrix,
) -> BandedResult {
    let qlen = query.len() as i32;
    let slen = subject.len() as i32;
    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();

    // Diagonal of the anchor
    let diag = query_anchor - subject_anchor;

    // Band boundaries
    let band_lo = diag - band_width;
    // DP arrays for the band
    let band_size = (2 * band_width + 1) as usize;
    let mut prev_row = vec![0i32; band_size + 2];
    let mut curr_row = vec![0i32; band_size + 2];
    let mut hgap = vec![i32::MIN / 2; band_size + 2];

    let mut best_score = 0i32;
    let mut best_i = 0i32;
    let mut best_j = 0i32;

    for j in 0..slen {
        let mut vgap = i32::MIN / 2;

        for band_idx in 0..band_size {
            let i = j + band_lo + band_idx as i32;
            if i < 0 || i >= qlen {
                curr_row[band_idx] = 0;
                continue;
            }

            let match_score = score_matrix.score(
                query[i as usize] & LETTER_MASK,
                subject[j as usize] & LETTER_MASK,
            );

            let diag_score = prev_row[band_idx] + match_score;
            let from_hgap = hgap[band_idx];
            let from_vgap = vgap;

            let s = diag_score.max(from_hgap).max(from_vgap).max(0);

            let open = s - gap_open;
            vgap = (vgap - gap_extend).max(open);
            hgap[band_idx] = (hgap[band_idx] - gap_extend).max(open);

            curr_row[band_idx] = s;

            if s > best_score {
                best_score = s;
                best_i = i;
                best_j = j;
            }
        }

        std::mem::swap(&mut prev_row, &mut curr_row);
        curr_row.fill(0);
    }

    // Note: begin positions are approximate — proper traceback not yet implemented.
    // Use smith_waterman::smith_waterman() for exact coordinates with traceback.
    BandedResult {
        score: best_score,
        query_begin: 0,
        query_end: best_i + 1,
        subject_begin: 0,
        subject_end: best_j + 1,
    }
}

/// Banded Needleman-Wunsch (global within band).
pub fn banded_nw(
    query: &[Letter],
    subject: &[Letter],
    band_width: i32,
    score_matrix: &ScoreMatrix,
) -> i32 {
    let qlen = query.len() as i32;
    let slen = subject.len() as i32;
    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();

    let band_size = (2 * band_width + 1) as usize;
    let mut prev_row = vec![i32::MIN / 2; band_size + 2];
    let mut curr_row = vec![i32::MIN / 2; band_size + 2];
    let mut hgap = vec![i32::MIN / 2; band_size + 2];

    // Initialize center of band
    prev_row[band_width as usize] = 0;

    for j in 0..slen {
        let mut vgap = i32::MIN / 2;

        for band_idx in 0..band_size {
            let i = j - band_width + band_idx as i32;
            if i < 0 || i >= qlen {
                curr_row[band_idx] = i32::MIN / 2;
                continue;
            }

            let match_score = score_matrix.score(
                query[i as usize] & LETTER_MASK,
                subject[j as usize] & LETTER_MASK,
            );

            let diag_score = if band_idx > 0 {
                prev_row[band_idx - 1] + match_score
            } else {
                i32::MIN / 2
            };
            let from_hgap = hgap[band_idx];
            let from_vgap = vgap;

            let s = diag_score.max(from_hgap).max(from_vgap);

            let open = s - gap_open;
            vgap = (vgap - gap_extend).max(open);
            hgap[band_idx] = (hgap[band_idx] - gap_extend).max(open);

            curr_row[band_idx] = s;
        }

        std::mem::swap(&mut prev_row, &mut curr_row);
        curr_row.fill(i32::MIN / 2);
    }

    // Find best score at the end of the subject
    *prev_row.iter().max().unwrap_or(&0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn test_banded_sw_identical() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..20).map(|i| i as Letter).collect();
        let result = banded_smith_waterman(&query, &query, 10, 10, 5, &sm);
        assert!(result.score > 0, "Self-alignment should have positive score");
    }

    #[test]
    fn test_banded_sw_different() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = vec![0; 20]; // All A
        let subject: Vec<Letter> = vec![13; 20]; // All F
        let result = banded_smith_waterman(&query, &subject, 10, 10, 10, &sm);
        // A-F score in BLOSUM62 is -2, so local alignment should be 0 or very low
        assert!(result.score <= 0);
    }

    #[test]
    fn test_banded_sw_wide_band() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = (0..10).map(|i| i as Letter).collect();
        let subject: Vec<Letter> = (0..10).map(|i| i as Letter).collect();
        // Wide band should give similar result to full DP
        let banded = banded_smith_waterman(&query, &subject, 5, 5, 10, &sm);
        let full = super::super::smith_waterman::smith_waterman(&query, &subject, &sm);
        assert_eq!(
            banded.score, full.score,
            "Wide band should match full DP for identical sequences"
        );
    }
}
