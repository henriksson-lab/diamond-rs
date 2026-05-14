use crate::basic::packed_transcript::EditOperation;
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
    pub length: i32,
    pub identities: i32,
    pub mismatches: i32,
    pub gap_openings: i32,
    pub gaps: i32,
    pub operations: Vec<(EditOperation, i32)>,
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
    let neg_inf = i32::MIN / 4;

    let diag = query_anchor - subject_anchor;
    let rows = qlen as usize + 1;
    let cols = slen as usize + 1;
    let idx = |i: usize, j: usize| -> usize { j * rows + i };
    let mut h = vec![0i32; rows * cols];
    let mut e = vec![neg_inf; rows * cols];
    let mut f = vec![neg_inf; rows * cols];

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
        for i in (lower as usize + 1)..=(upper as usize + 1) {
            let match_score =
                score_matrix.score(query[i - 1] & LETTER_MASK, subject[j - 1] & LETTER_MASK);
            let current = idx(i, j);
            let diag_score = h[idx(i - 1, j - 1)] + match_score;
            e[current] = (h[idx(i, j - 1)] - gap_open).max(e[idx(i, j - 1)] - gap_extend);
            f[current] = (h[idx(i - 1, j)] - gap_open).max(f[idx(i - 1, j)] - gap_extend);
            let s = diag_score.max(e[current]).max(f[current]).max(0);
            h[current] = s;

            if s > best_score {
                best_score = s;
                best_i = i;
                best_j = j;
            }
        }
    }

    if best_score == 0 {
        return BandedResult::default();
    }

    let mut result = BandedResult {
        score: best_score,
        query_end: best_i as i32,
        subject_end: best_j as i32,
        ..Default::default()
    };
    let mut i = best_i;
    let mut j = best_j;
    let mut ops = Vec::new();

    while i > 0 && j > 0 && h[idx(i, j)] > 0 {
        let score = h[idx(i, j)];
        let match_score =
            score_matrix.score(query[i - 1] & LETTER_MASK, subject[j - 1] & LETTER_MASK);
        let diag_score = h[idx(i - 1, j - 1)] + match_score;
        if score == diag_score {
            if (query[i - 1] & LETTER_MASK) == (subject[j - 1] & LETTER_MASK) {
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

/// Banded Needleman-Wunsch (global within band).
pub fn banded_nw(
    query: &[Letter],
    subject: &[Letter],
    band_width: i32,
    score_matrix: &ScoreMatrix,
) -> i32 {
    let qlen = query.len();
    let slen = subject.len();
    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();
    let neg_inf = i32::MIN / 4;
    let rows = qlen + 1;
    let idx = |i: usize, j: usize| -> usize { j * rows + i };

    let mut h = vec![neg_inf; rows * (slen + 1)];
    let mut e = vec![neg_inf; rows * (slen + 1)];
    let mut f = vec![neg_inf; rows * (slen + 1)];
    h[idx(0, 0)] = 0;

    for i in 1..=qlen {
        if i as i32 <= band_width {
            f[idx(i, 0)] = -gap_open - (i as i32 - 1) * gap_extend;
            h[idx(i, 0)] = f[idx(i, 0)];
        }
    }
    for j in 1..=slen {
        if j as i32 <= band_width {
            e[idx(0, j)] = -gap_open - (j as i32 - 1) * gap_extend;
            h[idx(0, j)] = e[idx(0, j)];
        }
    }

    for j in 1..=slen {
        let lower = (j as i32 - band_width).max(1) as usize;
        let upper = (j as i32 + band_width).min(qlen as i32) as usize;
        if lower > upper {
            continue;
        }
        for i in lower..=upper {
            let match_score =
                score_matrix.score(query[i - 1] & LETTER_MASK, subject[j - 1] & LETTER_MASK);
            let current = idx(i, j);
            e[current] = (h[idx(i, j - 1)] - gap_open).max(e[idx(i, j - 1)] - gap_extend);
            f[current] = (h[idx(i - 1, j)] - gap_open).max(f[idx(i - 1, j)] - gap_extend);
            h[current] = (h[idx(i - 1, j - 1)] + match_score)
                .max(e[current])
                .max(f[current]);
        }
    }

    h[idx(qlen, slen)]
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
        assert!(
            result.score > 0,
            "Self-alignment should have positive score"
        );
        assert_eq!(result.query_begin, 0);
        assert_eq!(result.subject_begin, 0);
        assert_eq!(result.query_end, query.len() as i32);
        assert_eq!(result.subject_end, query.len() as i32);
        assert_eq!(result.identities, query.len() as i32);
        assert_eq!(result.gaps, 0);
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

    #[test]
    fn test_banded_nw_wide_band_matches_full_nw() {
        let sm = make_test_matrix();
        let query: Vec<Letter> = vec![0, 1, 2, 3, 4, 5];
        let subject: Vec<Letter> = vec![0, 1, 3, 4, 5];
        let banded = banded_nw(&query, &subject, 6, &sm);
        let (full, _) = super::super::smith_waterman::needleman_wunsch(&query, &subject, &sm);
        assert_eq!(banded, full);
    }
}
