use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

/// Result of a Smith-Waterman alignment.
#[derive(Debug, Clone, Default)]
pub struct SwResult {
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
    /// Edit operations (reversed, from end to start).
    pub operations: Vec<(EditOperation, i32)>,
}

/// DP score buffer for Needleman-Wunsch / Smith-Waterman.
struct DpMatrix {
    data: Vec<i32>,
    rows: usize, // query_len + 1
    cols: usize, // subject_len + 1
}

impl DpMatrix {
    fn new(query_len: usize, subject_len: usize) -> Self {
        let rows = query_len + 1;
        let cols = subject_len + 1;
        DpMatrix {
            data: vec![0; rows * cols],
            rows,
            cols,
        }
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> i32 {
        self.data[j * self.rows + i]
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, val: i32) {
        self.data[j * self.rows + i] = val;
    }

    fn find_max(&self) -> (i32, usize, usize) {
        let mut max_score = 0;
        let mut max_i = 0;
        let mut max_j = 0;
        for j in 0..self.cols {
            for i in 0..self.rows {
                let s = self.get(i, j);
                if s > max_score {
                    max_score = s;
                    max_i = i;
                    max_j = j;
                }
            }
        }
        (max_score, max_i, max_j)
    }
}

/// Local Smith-Waterman alignment (scalar, full DP matrix).
///
/// Computes the optimal local alignment between query and subject sequences
/// using affine gap penalties. Returns the alignment result with traceback.
pub fn smith_waterman(
    query: &[Letter],
    subject: &[Letter],
    score_matrix: &ScoreMatrix,
) -> SwResult {
    let qlen = query.len();
    let slen = subject.len();
    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();

    // Fill DP matrix
    let mut dp = DpMatrix::new(qlen, slen);
    let mut hgap = vec![i32::MIN / 2; qlen + 1]; // horizontal gap scores

    for j in 1..=slen {
        let mut vgap = i32::MIN / 2; // vertical gap score
        for i in 1..=qlen {
            let match_score = score_matrix.score(
                query[i - 1] & LETTER_MASK,
                subject[j - 1] & LETTER_MASK,
            );
            let diag = dp.get(i - 1, j - 1) + match_score;
            let s = diag.max(vgap).max(hgap[i]).max(0);

            let open = s - gap_open;
            vgap = (vgap - gap_extend).max(open);
            hgap[i] = (hgap[i] - gap_extend).max(open);

            dp.set(i, j, s);
        }
    }

    // Find maximum score
    let (max_score, max_i, max_j) = dp.find_max();

    if max_score == 0 {
        return SwResult::default();
    }

    // Traceback
    let mut result = SwResult {
        score: max_score,
        query_end: max_i as i32,
        subject_end: max_j as i32,
        ..Default::default()
    };

    let mut i = max_i;
    let mut j = max_j;
    let mut ops = Vec::new();

    while dp.get(i, j) > 0 && i > 0 && j > 0 {
        let score = dp.get(i, j);
        let match_score = score_matrix.score(
            query[i - 1] & LETTER_MASK,
            subject[j - 1] & LETTER_MASK,
        );
        let diag = dp.get(i - 1, j - 1) + match_score;

        if score == diag {
            // Match or substitution
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
        } else if has_hgap(&dp, i, j, gap_open, gap_extend) {
            // Horizontal gap (deletion in query, insertion in subject)
            let mut gap_len = 1;
            let expected = dp.get(i, j - 1) - gap_open;
            if score == expected {
                // gap of length 1
            } else {
                // Find gap length
                for k in 2..=j {
                    if score == dp.get(i, j - k) - gap_open - (k as i32 - 1) * gap_extend {
                        gap_len = k as i32;
                        break;
                    }
                }
            }
            ops.push((EditOperation::Deletion, gap_len));
            result.length += gap_len;
            result.gaps += gap_len;
            result.gap_openings += 1;
            j -= gap_len as usize;
        } else {
            // Vertical gap (insertion in query, deletion in subject)
            let mut gap_len = 1i32;
            for k in 2..=i {
                if score == dp.get(i - k, j) - gap_open - (k as i32 - 1) * gap_extend {
                    gap_len = k as i32;
                    break;
                }
            }
            ops.push((EditOperation::Insertion, gap_len));
            result.length += gap_len;
            result.gaps += gap_len;
            result.gap_openings += 1;
            i -= gap_len as usize;
        }
    }

    result.query_begin = i as i32;
    result.subject_begin = j as i32;
    ops.reverse();
    result.operations = ops;
    result
}

fn has_hgap(dp: &DpMatrix, i: usize, j: usize, gap_open: i32, gap_extend: i32) -> bool {
    let score = dp.get(i, j);
    for k in 1..=j {
        let expected = dp.get(i, j - k) - gap_open - (k as i32 - 1) * gap_extend;
        if score == expected {
            return true;
        }
        if expected < 0 {
            break;
        }
    }
    false
}

/// Global Needleman-Wunsch alignment (scalar).
pub fn needleman_wunsch(
    query: &[Letter],
    subject: &[Letter],
    score_matrix: &ScoreMatrix,
) -> (i32, Vec<(EditOperation, i32)>) {
    let qlen = query.len();
    let slen = subject.len();
    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();

    let mut dp = DpMatrix::new(qlen, slen);

    // Initialize first column and row with gap penalties
    for i in 1..=qlen {
        dp.set(i, 0, -(gap_open + (i as i32 - 1) * gap_extend));
    }
    for j in 1..=slen {
        dp.set(0, j, -(gap_open + (j as i32 - 1) * gap_extend));
    }

    let mut hgap = vec![i32::MIN / 2; qlen + 1];

    for j in 1..=slen {
        let mut vgap = i32::MIN / 2;
        for i in 1..=qlen {
            let match_score = score_matrix.score(
                query[i - 1] & LETTER_MASK,
                subject[j - 1] & LETTER_MASK,
            );
            let diag = dp.get(i - 1, j - 1) + match_score;
            let s = diag.max(vgap).max(hgap[i]);

            let open = s - gap_open;
            vgap = (vgap - gap_extend).max(open);
            hgap[i] = (hgap[i] - gap_extend).max(open);

            dp.set(i, j, s);
        }
    }

    let final_score = dp.get(qlen, slen);
    (final_score, Vec::new()) // NW traceback not implemented — use smith_waterman for traced alignments
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn test_smith_waterman_identical() {
        let sm = make_test_matrix();
        // A R N D = 0 1 2 3
        let query = vec![0i8, 1, 2, 3];
        let subject = vec![0i8, 1, 2, 3];
        let result = smith_waterman(&query, &subject, &sm);
        assert!(result.score > 0);
        assert_eq!(result.identities, 4);
        assert_eq!(result.length, 4);
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.gaps, 0);
    }

    #[test]
    fn test_smith_waterman_mismatch() {
        let sm = make_test_matrix();
        // Query: A A A A (0 0 0 0)
        // Subject: A R A A (0 1 0 0)
        let query = vec![0i8, 0, 0, 0];
        let subject = vec![0i8, 1, 0, 0];
        let result = smith_waterman(&query, &subject, &sm);
        assert!(result.score > 0);
        assert!(result.identities >= 3);
    }

    #[test]
    fn test_smith_waterman_no_match() {
        let sm = make_test_matrix();
        // Very different short sequences
        let query = vec![0i8]; // A
        let subject = vec![17i8]; // W
        let result = smith_waterman(&query, &subject, &sm);
        // A-W score in BLOSUM62 is -3, so no local alignment
        assert_eq!(result.score, 0);
    }

    #[test]
    fn test_smith_waterman_self_alignment() {
        let sm = make_test_matrix();
        // ARNDCQEGHILKMFPST (0-16)
        let seq: Vec<Letter> = (0..17).map(|i| i as Letter).collect();
        let result = smith_waterman(&seq, &seq, &sm);
        assert_eq!(result.identities, 17);
        assert_eq!(result.length, 17);
        assert!(result.score > 50); // Self-score should be high
    }

    #[test]
    fn test_needleman_wunsch() {
        let sm = make_test_matrix();
        let query = vec![0i8, 1, 2, 3];
        let subject = vec![0i8, 1, 2, 3];
        let (score, _) = needleman_wunsch(&query, &subject, &sm);
        // Same as SW for identical sequences with no gaps needed
        assert!(score > 0);
    }
}
