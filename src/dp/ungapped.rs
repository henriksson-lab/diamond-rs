use crate::basic::value::{Letter, Score, DELIMITER_LETTER, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::interval::Interval;

/// Result of an ungapped extension.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
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

    pub fn query_range(&self) -> Interval {
        Interval::new(self.i, self.i + self.len)
    }

    pub fn subject_range(&self) -> Interval {
        Interval::new(self.j, self.j + self.len)
    }

    pub fn subject_begin(&self) -> i32 {
        self.j
    }

    pub fn subject_last(&self) -> i32 {
        self.j + self.len - 1
    }

    pub fn query_last(&self) -> i32 {
        self.i + self.len - 1
    }

    pub fn subject_end(&self) -> i32 {
        self.j + self.len
    }

    pub fn query_begin(&self) -> i32 {
        self.i
    }

    pub fn diag(&self) -> i32 {
        self.i - self.j
    }

    pub fn id_percent(&self) -> f64 {
        self.identities as f64 / self.len as f64 * 100.0
    }

    pub fn cov_percent(&self, seq_len: i32) -> f64 {
        self.len as f64 / seq_len as f64 * 100.0
    }

    pub fn band_interval(&self, band: i32) -> Interval {
        Interval::new(self.diag() - band, self.diag() + band)
    }

    pub fn set_query_end(&mut self, i: i32) {
        self.len = i - self.i;
    }

    pub fn set_target_end(&mut self, j: i32) {
        self.len = j - self.j;
    }

    pub fn intersect(&self, x: &DiagonalSegment) -> DiagonalSegment {
        if self.diag() != x.diag() {
            DiagonalSegment::default()
        } else {
            let q = self.query_range().intersect(&x.query_range());
            DiagonalSegment::new(
                q.begin,
                self.subject_range().intersect(&x.subject_range()).begin,
                q.length(),
                0,
            )
        }
    }

    pub fn is_enveloped(&self, x: &DiagonalSegment) -> bool {
        self.score <= x.score
            && self.query_range().overlap_factor(&x.query_range()) == 1.0
            && self.subject_range().overlap_factor(&x.subject_range()) == 1.0
    }

    pub fn transpose(&self) -> DiagonalSegment {
        DiagonalSegment::new(self.j, self.i, self.len, self.score)
    }

    pub fn partial_score(&self, diff: i32) -> i32 {
        self.score * (self.len - diff).max(0) / self.len
    }

    pub fn leq_position(&self, rhs: &DiagonalSegment) -> bool {
        self.i + self.len <= rhs.i && self.j + self.len <= rhs.j
    }

    pub fn cmp_subject(x: &DiagonalSegment, y: &DiagonalSegment) -> bool {
        x.j < y.j || (x.j == y.j && x.i < y.i)
    }

    pub fn cmp_score(x: &DiagonalSegment, y: &DiagonalSegment) -> bool {
        x.score > y.score
    }

    pub fn cmp_subject_end(x: &DiagonalSegment, y: &DiagonalSegment) -> bool {
        x.subject_end() < y.subject_end()
    }

    pub fn cmp_heuristic(x: &DiagonalSegment, y: &DiagonalSegment) -> bool {
        (x.subject_end() < y.subject_end() && x.j < y.j)
            || (x.j - y.j < y.subject_end() - x.subject_end())
    }

    pub fn cmp_diag(x: &DiagonalSegment, y: &DiagonalSegment) -> bool {
        x.diag() < y.diag() || (x.diag() == y.diag() && x.j < y.j)
    }

    pub fn cmp_len(x: &DiagonalSegment, y: &DiagonalSegment) -> bool {
        x.len > y.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

pub fn abs_shift(x: &DiagonalSegment, y: &DiagonalSegment) -> i32 {
    (x.diag() - y.diag()).abs()
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

    DiagonalSegment::new(qa as i32 - delta, sa as i32 - delta, len + delta, score)
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
        score += score_matrix.score(query[i] & LETTER_MASK, subject[j] & LETTER_MASK);
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
        st += score_matrix.score(query[q] & LETTER_MASK, subject[s] & LETTER_MASK);
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
        st += score_matrix.score(query[n] & LETTER_MASK, subject[n] & LETTER_MASK);
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

    #[test]
    fn test_diagonal_segment_cpp_methods() {
        let mut ds = DiagonalSegment::with_identities(5, 10, 20, 100, 15);

        assert_eq!(ds.query_range(), Interval::new(5, 25));
        assert_eq!(ds.subject_range(), Interval::new(10, 30));
        assert_eq!(ds.subject_begin(), 10);
        assert_eq!(ds.subject_last(), 29);
        assert_eq!(ds.query_last(), 24);
        assert_eq!(ds.query_begin(), 5);
        assert_eq!(ds.diag(), -5);
        assert_eq!(ds.id_percent(), 75.0);
        assert_eq!(ds.cov_percent(40), 50.0);
        assert_eq!(ds.band_interval(3), Interval::new(-8, -2));
        assert_eq!(ds.partial_score(5), 75);

        ds.set_query_end(22);
        assert_eq!(ds.len, 17);
        ds.set_target_end(24);
        assert_eq!(ds.len, 14);
    }

    #[test]
    fn test_diagonal_segment_intersect_envelope_and_comparators() {
        let a = DiagonalSegment::new(5, 10, 10, 100);
        let b = DiagonalSegment::new(8, 13, 10, 80);
        let c = DiagonalSegment::new(20, 21, 5, 20);
        let env = DiagonalSegment::new(4, 9, 20, 120);

        assert_eq!(a.intersect(&b), DiagonalSegment::new(8, 13, 7, 0));
        assert!(a.intersect(&c).is_empty());
        assert!(a.is_enveloped(&env));
        assert_eq!(a.transpose(), DiagonalSegment::new(10, 5, 10, 100));
        assert!(a.leq_position(&c));
        assert!(DiagonalSegment::cmp_subject(&a, &b));
        assert!(DiagonalSegment::cmp_score(&a, &b));
        assert!(DiagonalSegment::cmp_subject_end(&a, &c));
        assert!(DiagonalSegment::cmp_heuristic(&a, &c));
        assert!(DiagonalSegment::cmp_diag(&a, &c));
        assert!(DiagonalSegment::cmp_len(&a, &c));
        assert_eq!(abs_shift(&a, &c), 4);
    }
}
