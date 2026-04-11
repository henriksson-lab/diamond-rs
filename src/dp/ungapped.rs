use crate::basic::value::{Letter, Score, DELIMITER_LETTER, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

/// Result of an ungapped extension.
#[derive(Debug, Clone, Default)]
pub struct DiagonalSegment {
    /// Query start position.
    pub i: i32,
    /// Subject start position.
    pub j: i32,
    /// Alignment length.
    pub len: i32,
    /// Alignment score.
    pub score: i32,
    /// Number of identities (optional).
    pub identities: i32,
}

impl DiagonalSegment {
    pub fn new(i: i32, j: i32, len: i32, score: i32) -> Self {
        DiagonalSegment {
            i,
            j,
            len,
            score,
            identities: 0,
        }
    }

    pub fn with_identities(i: i32, j: i32, len: i32, score: i32, identities: i32) -> Self {
        DiagonalSegment {
            i,
            j,
            len,
            score,
            identities,
        }
    }

    pub fn query_end(&self) -> i32 {
        self.i + self.len
    }

    pub fn subject_end(&self) -> i32 {
        self.j + self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// X-drop ungapped alignment extending in both directions from a seed hit.
///
/// This is the core ungapped extension used to filter seed hits before
/// gapped alignment.
pub fn xdrop_ungapped(
    query: &[Letter],
    subject: &[Letter],
    qa: usize,
    sa: usize,
    xdrop: i32,
    score_matrix: &ScoreMatrix,
) -> DiagonalSegment {
    let mut score: i32 = 0;
    let mut st: i32 = 0;
    let mut delta: i32 = 0;
    let mut len: i32 = 0;
    let mut n: i32 = 1;

    // Extend left
    let mut q = qa as i32 - 1;
    let mut s = sa as i32 - 1;
    while q >= 0
        && s >= 0
        && score - st < xdrop
        && query[q as usize] != DELIMITER_LETTER
        && subject[s as usize] != DELIMITER_LETTER
    {
        let ql = query[q as usize] & LETTER_MASK;
        let sl = subject[s as usize] & LETTER_MASK;
        st += score_matrix.score(ql, sl);
        if st > score {
            score = st;
            delta = n;
        }
        q -= 1;
        s -= 1;
        n += 1;
    }

    // Extend right
    let mut q = qa;
    let mut s = sa;
    st = score;
    n = 1;
    while q < query.len()
        && s < subject.len()
        && score - st < xdrop
        && query[q] != DELIMITER_LETTER
        && subject[s] != DELIMITER_LETTER
    {
        let ql = query[q] & LETTER_MASK;
        let sl = subject[s] & LETTER_MASK;
        st += score_matrix.score(ql, sl);
        if st > score {
            score = st;
            len = n;
        }
        q += 1;
        s += 1;
        n += 1;
    }

    DiagonalSegment::new(
        qa as i32 - delta,
        sa as i32 - delta,
        len + delta,
        score,
    )
}

/// X-drop ungapped extension with identity counting.
pub fn xdrop_ungapped_with_identities(
    query: &[Letter],
    subject: &[Letter],
    qa: usize,
    sa: usize,
    xdrop: i32,
    score_matrix: &ScoreMatrix,
) -> DiagonalSegment {
    let mut score: i32 = 0;
    let mut st: i32 = 0;
    let mut delta: i32 = 0;
    let mut len: i32 = 0;
    let mut ident: i32 = 0;
    let mut n: i32 = 1;
    let mut left_ident: i32 = 0;

    // Extend left
    let mut q = qa as i32 - 1;
    let mut s = sa as i32 - 1;
    while q >= 0
        && s >= 0
        && score - st < xdrop
        && query[q as usize] != DELIMITER_LETTER
        && subject[s as usize] != DELIMITER_LETTER
    {
        let ql = query[q as usize] & LETTER_MASK;
        let sl = subject[s as usize] & LETTER_MASK;
        st += score_matrix.score(ql, sl);
        if ql == sl {
            left_ident += 1;
        }
        if st > score {
            score = st;
            delta = n;
            ident += left_ident;
            left_ident = 0;
        }
        q -= 1;
        s -= 1;
        n += 1;
    }

    // Extend right
    let mut q = qa;
    let mut s = sa;
    st = score;
    n = 1;
    let mut right_ident: i32 = 0;
    while q < query.len()
        && s < subject.len()
        && score - st < xdrop
        && query[q] != DELIMITER_LETTER
        && subject[s] != DELIMITER_LETTER
    {
        let ql = query[q] & LETTER_MASK;
        let sl = subject[s] & LETTER_MASK;
        st += score_matrix.score(ql, sl);
        if ql == sl {
            right_ident += 1;
        }
        if st > score {
            score = st;
            len = n;
            ident += right_ident;
            right_ident = 0;
        }
        q += 1;
        s += 1;
        n += 1;
    }

    DiagonalSegment::with_identities(
        qa as i32 - delta,
        sa as i32 - delta,
        len + delta,
        score,
        ident,
    )
}

/// Compute score over a fixed range of a diagonal.
pub fn score_range(
    query: &[Letter],
    subject: &[Letter],
    i_begin: usize,
    j_begin: usize,
    j_end: usize,
    score_matrix: &ScoreMatrix,
) -> i32 {
    let mut score: i32 = 0;
    for (offset, j) in (j_begin..j_end).enumerate() {
        let i = i_begin + offset;
        score += score_matrix.score(
            query[i] & LETTER_MASK,
            subject[j] & LETTER_MASK,
        );
    }
    score
}

/// Compute the self-alignment score of a sequence (best local score of seq vs itself).
pub fn self_score(seq: &[Letter], score_matrix: &ScoreMatrix) -> Score {
    let mut s: Score = 0;
    let mut sl: Score = 0;
    for &l in seq {
        let l = l & LETTER_MASK;
        sl += score_matrix.score(l, l);
        sl = sl.max(0);
        s = s.max(sl);
    }
    s
}

/// X-drop ungapped extension to the right only.
pub fn xdrop_ungapped_right(
    query: &[Letter],
    subject: &[Letter],
    xdrop: i32,
    score_matrix: &ScoreMatrix,
) -> (i32, i32) {
    let mut score: i32 = 0;
    let mut st: i32 = 0;
    let mut len: i32 = 0;
    let mut n: i32 = 1;
    let mut q = 0usize;
    let mut s = 0usize;

    while q < query.len()
        && s < subject.len()
        && score - st < xdrop
        && query[q] != DELIMITER_LETTER
        && subject[s] != DELIMITER_LETTER
    {
        st += score_matrix.score(
            query[q] & LETTER_MASK,
            subject[s] & LETTER_MASK,
        );
        if st > score {
            score = st;
            len = n;
        }
        q += 1;
        s += 1;
        n += 1;
    }
    (score, len)
}

/// Ungapped alignment within a fixed window.
pub fn ungapped_window(
    query: &[Letter],
    subject: &[Letter],
    window: usize,
    score_matrix: &ScoreMatrix,
) -> i32 {
    let mut score: i32 = 0;
    let mut st: i32 = 0;
    for n in 0..window.min(query.len()).min(subject.len()) {
        st += score_matrix.score(
            query[n] & LETTER_MASK,
            subject[n] & LETTER_MASK,
        );
        st = st.max(0);
        score = score.max(st);
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn test_xdrop_ungapped_identical() {
        let sm = make_test_matrix();
        // Two identical sequences: A R N D (letters 0,1,2,3)
        let query = vec![DELIMITER_LETTER, 0, 1, 2, 3, DELIMITER_LETTER];
        let subject = vec![DELIMITER_LETTER, 0, 1, 2, 3, DELIMITER_LETTER];
        let result = xdrop_ungapped(&query, &subject, 1, 1, 20, &sm);
        assert!(result.score > 0);
        assert!(result.len > 0);
    }

    #[test]
    fn test_xdrop_ungapped_different() {
        let sm = make_test_matrix();
        // Completely different sequences
        let query = vec![DELIMITER_LETTER, 0, 0, 0, 0, DELIMITER_LETTER];
        let subject = vec![DELIMITER_LETTER, 13, 13, 13, 13, DELIMITER_LETTER]; // F=13
        let result = xdrop_ungapped(&query, &subject, 1, 1, 5, &sm);
        // Should have low or zero score
        assert!(result.score <= 0 || result.len == 0);
    }

    #[test]
    fn test_self_score() {
        let sm = make_test_matrix();
        let seq = vec![0, 1, 2, 3]; // A R N D
        let s = self_score(&seq, &sm);
        // Self-score should be positive (4+5+6+6 = 21 for BLOSUM62 diagonal)
        assert!(s > 0);
    }

    #[test]
    fn test_score_range() {
        let sm = make_test_matrix();
        let query = vec![0, 0, 0]; // AAA
        let subject = vec![0, 0, 0]; // AAA
        let s = score_range(&query, &subject, 0, 0, 3, &sm);
        assert_eq!(s, 12); // A-A = 4 in BLOSUM62, 4*3=12
    }

    #[test]
    fn test_ungapped_window() {
        let sm = make_test_matrix();
        let query = vec![0, 1, 2, 3];
        let subject = vec![0, 1, 2, 3];
        let s = ungapped_window(&query, &subject, 4, &sm);
        assert!(s > 0);
    }

    #[test]
    fn test_diagonal_segment() {
        let ds = DiagonalSegment::new(5, 10, 20, 100);
        assert_eq!(ds.query_end(), 25);
        assert_eq!(ds.subject_end(), 30);
        assert!(!ds.is_empty());
    }
}
