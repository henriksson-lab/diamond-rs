use crate::align::hsp::Hsp;
use crate::basic::value::{Letter, Score, DELIMITER_LETTER, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::hsp::Anchor;
use crate::util::interval::Interval;
use crate::util::sequence::window_scores;

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

/// Matches C++ `xdrop_ungapped(const Sequence&, const int8_t*, ..., bool)`.
pub fn xdrop_ungapped_with_cbs(
    query: &[Letter],
    query_cbs: Option<&[i8]>,
    subject: &[Letter],
    qa: usize,
    sa: usize,
    xdrop: i32,
    count_identities: bool,
    score_matrix: &ScoreMatrix,
) -> DiagonalSegment {
    let mut score: i32 = 0;
    let mut st: i32 = 0;
    let mut delta: i32 = 0;
    let mut len: i32 = 0;
    let mut ident: i32 = 0;
    let mut n: i32 = 1;
    let mut pending_ident: i32 = 0;

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
        if let Some(cbs) = query_cbs {
            st += cbs[q as usize] as i32;
        }
        if count_identities && ql == sl {
            pending_ident += 1;
        }
        if st > score {
            score = st;
            delta = n;
            ident += pending_ident;
            pending_ident = 0;
        }
        q -= 1;
        s -= 1;
        n += 1;
    }

    let mut q = qa;
    let mut s = sa;
    st = score;
    n = 1;
    pending_ident = 0;
    while q < query.len()
        && s < subject.len()
        && score - st < xdrop
        && query[q] != DELIMITER_LETTER
        && subject[s] != DELIMITER_LETTER
    {
        let ql = query[q] & LETTER_MASK;
        let sl = subject[s] & LETTER_MASK;
        st += score_matrix.score(ql, sl);
        if let Some(cbs) = query_cbs {
            st += cbs[q] as i32;
        }
        if count_identities && ql == sl {
            pending_ident += 1;
        }
        if st > score {
            score = st;
            len = n;
            ident += pending_ident;
            pending_ident = 0;
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

/// Matches C++ `score_range_s(...)`.
pub fn score_range_s(
    query: &[Letter],
    query_cbs: Option<&[i8]>,
    subject: &[Letter],
    i_begin: usize,
    j_begin: usize,
    j_end: usize,
    score_matrix: &ScoreMatrix,
) -> DiagonalSegment {
    let mut score: i32 = 0;
    let mut i = i_begin;
    for j in j_begin..j_end {
        score += score_matrix.score(query[i] & LETTER_MASK, subject[j] & LETTER_MASK);
        if let Some(cbs) = query_cbs {
            score += cbs[i] as i32;
        }
        i += 1;
    }
    DiagonalSegment::new(
        i_begin as i32,
        j_begin as i32,
        (j_end - j_begin) as i32,
        score,
    )
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

/// Matches C++ `xdrop_anchored(...)`.
pub fn xdrop_anchored_right(
    query: &[Letter],
    subject: &[Letter],
    len: usize,
    xdrop: i32,
    score_matrix: &ScoreMatrix,
) -> (Score, i32) {
    if len == 0 {
        return (0, 0);
    }
    let mut max_score: Score = 0;
    let mut score: Score = 0;
    let mut max_n: i32 = 0;
    let mut n: usize = 0;
    loop {
        score += score_matrix.score(query[n] & LETTER_MASK, subject[n] & LETTER_MASK);
        n += 1;
        if score > max_score {
            max_score = score;
            max_n = n as i32;
        }
        if n >= len || max_score - score >= xdrop {
            break;
        }
    }
    (max_score, max_n)
}

/// Matches C++ `xdrop_anchored(...)`.
pub fn xdrop_anchored_left(
    query: &[Letter],
    subject: &[Letter],
    i_last: usize,
    j_last: usize,
    len: usize,
    xdrop: i32,
    score_matrix: &ScoreMatrix,
) -> (Score, i32) {
    if len == 0 {
        return (0, 0);
    }
    let mut max_score: Score = 0;
    let mut score: Score = 0;
    let mut max_n: i32 = 0;
    let mut n: usize = 0;
    loop {
        score += score_matrix.score(
            query[i_last - n] & LETTER_MASK,
            subject[j_last - n] & LETTER_MASK,
        );
        n += 1;
        if score > max_score {
            max_score = score;
            max_n = n as i32;
        }
        if n >= len || max_score - score >= xdrop {
            break;
        }
    }
    (max_score, max_n)
}

/// Matches C++ `xdrop_ungapped(const Sequence&, const Sequence&, const DiagonalSegment&)`.
pub fn xdrop_ungapped_from_anchor(
    query: &[Letter],
    subject: &[Letter],
    anchor: &DiagonalSegment,
    xdrop: i32,
    score_matrix: &ScoreMatrix,
) -> DiagonalSegment {
    let left_len = (anchor.i.min(anchor.j)).max(0) as usize;
    let left = if left_len == 0 {
        (0, 0)
    } else {
        xdrop_anchored_left(
            query,
            subject,
            anchor.i as usize - 1,
            anchor.j as usize - 1,
            left_len,
            xdrop,
            score_matrix,
        )
    };
    let right = xdrop_anchored_right(
        &query[anchor.query_end() as usize..],
        &subject[anchor.subject_end() as usize..],
        (query.len() as i32 - anchor.query_end())
            .min(subject.len() as i32 - anchor.subject_end())
            .max(0) as usize,
        xdrop,
        score_matrix,
    );
    DiagonalSegment::new(
        anchor.i - left.1,
        anchor.j - left.1,
        anchor.len + left.1 + right.1,
        left.0
            + right.0
            + score_range(
                query,
                subject,
                anchor.i as usize,
                anchor.j as usize,
                anchor.subject_end() as usize,
                score_matrix,
            ),
    )
}

/// Matches C++ `make_clipped_anchor(...)`.
pub fn make_clipped_anchor(
    anchor: &Anchor,
    query: &[Letter],
    query_cbs: Option<&[i8]>,
    target: &[Letter],
    anchor_window: i32,
    anchor_score: f64,
    score_matrix: &ScoreMatrix,
) -> Anchor {
    let q_begin = anchor.segment.query_begin() as usize;
    let q_end = anchor.segment.query_end() as usize;
    let s_begin = anchor.segment.subject_begin() as usize;
    let s_end = anchor.segment.subject_end() as usize;
    let q = &query[q_begin..q_end];
    let t = &target[s_begin..s_end];
    let scores = window_scores(q, t, anchor_window, |a, b| score_matrix.score(a, b));
    if scores.is_empty() {
        return Anchor::from_diagonal_segment(DiagonalSegment::default());
    }
    let cutoff = (anchor_score * anchor_window as f64).round() as Score;
    let mut max_idx = 0usize;
    for idx in 1..scores.len() {
        if scores[idx] > scores[max_idx] {
            max_idx = idx;
        }
    }
    let d1 = ((max_idx + 1)..scores.len())
        .find(|&idx| scores[idx] < cutoff)
        .unwrap_or(scores.len()) as i32;
    let left_idx = (0..max_idx).rev().find(|&idx| scores[idx] < cutoff);
    let mut d0 = left_idx
        .map(|idx| idx as i32 - anchor_window + 2)
        .unwrap_or(1 - anchor_window)
        .max(0);
    let mut d1 = d1;
    while (d0 as usize) < q.len() && q[d0 as usize] != t[d0 as usize] {
        d0 += 1;
    }
    while d1 > 0 && q[(d1 - 1) as usize] != t[(d1 - 1) as usize] {
        d1 -= 1;
    }
    if d1 <= d0 {
        return Anchor::from_diagonal_segment(DiagonalSegment::default());
    }
    let clipped_anchor = score_range_s(
        query,
        query_cbs,
        target,
        (anchor.segment.query_begin() + d0) as usize,
        (anchor.segment.subject_begin() + d0) as usize,
        (anchor.segment.subject_begin() + d1) as usize,
        score_matrix,
    );
    let clipped_score = score_range_s(
        query,
        query_cbs,
        target,
        (anchor.segment.query_begin() + d1) as usize,
        (anchor.segment.subject_begin() + d1) as usize,
        anchor.segment.subject_end() as usize,
        score_matrix,
    )
    .score;
    Anchor::new(
        clipped_anchor,
        anchor.d_min_left,
        anchor.d_max_left,
        anchor.d_min_right,
        anchor.d_max_right,
        anchor.prefix_score - clipped_score,
    )
}

/// Matches C++ `make_null_anchor(const Anchor&)`.
pub fn make_null_anchor(anchor: &Anchor) -> Anchor {
    Anchor::from_diagonal_segment(DiagonalSegment::new(
        anchor.segment.i + anchor.segment.len / 2,
        anchor.segment.j + anchor.segment.len / 2,
        0,
        0,
    ))
}

/// Matches C++ `trivial(Sequence query, Sequence target, Loc dq, Loc dt, ...)`.
pub fn trivial_at(
    query: &[Letter],
    target: &[Letter],
    dq: usize,
    dt: usize,
    query_cbs: Option<&[i8]>,
    max_evalue: f64,
    score_matrix: &ScoreMatrix,
) -> Hsp {
    const WINDOW: usize = 40;
    const ID: u32 = 30;
    let l = (query.len() - dq).min(target.len() - dt);
    let bits = (1u64 << WINDOW) - 1;
    let mut n = 0usize;
    let mut score: Score = 0;
    let mut mask = 0u64;
    for i in 0..l {
        let eq = (query[i + dq] == target[i + dt]) as u64;
        mask = ((mask << 1) | eq) & bits;
        n += 1;
        if n >= WINDOW && mask.count_ones() < ID {
            return Hsp::new();
        }
        score += score_matrix.score(query[i + dq] & LETTER_MASK, target[i + dt] & LETTER_MASK);
        if let Some(cbs) = query_cbs {
            score += cbs[i + dq] as i32;
        }
    }
    let evalue = score_matrix.evalue(score, query.len() as u32, target.len() as u32);
    if evalue > max_evalue {
        return Hsp::new();
    }
    if l < WINDOW && (mask.count_ones() as f64) / (l as f64) < (ID as usize / WINDOW) as f64 {
        return Hsp::new();
    }
    let mut hsp = Hsp::new();
    hsp.score = score;
    hsp.query_range = Interval::new(dq as i32, (dq + l) as i32);
    hsp.query_source_range = hsp.query_range;
    hsp.subject_range = Interval::new(dt as i32, (dt + l) as i32);
    hsp.evalue = evalue;
    hsp.bit_score = score_matrix.bitscore(score as f64);
    hsp
}

/// Matches C++ `trivial(Sequence query, Sequence target, ...)`.
pub fn trivial(
    query: &[Letter],
    target: &[Letter],
    query_cbs: Option<&[i8]>,
    max_evalue: f64,
    score_matrix: &ScoreMatrix,
) -> Hsp {
    if query.len() <= target.len() {
        for i in 0..=target.len() - query.len() {
            let hsp = trivial_at(query, target, 0, i, query_cbs, max_evalue, score_matrix);
            if hsp.score != 0 {
                return hsp;
            }
        }
    } else {
        for i in 0..=query.len() - target.len() {
            let hsp = trivial_at(query, target, i, 0, query_cbs, max_evalue, score_matrix);
            if hsp.score != 0 {
                return hsp;
            }
        }
    }
    Hsp::new()
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
    fn test_xdrop_ungapped_with_cbs_and_identity_flag() {
        let sm = make_test_matrix();
        let query = vec![DELIMITER_LETTER, 0, 1, 2, DELIMITER_LETTER];
        let subject = vec![DELIMITER_LETTER, 0, 1, 2, DELIMITER_LETTER];
        let cbs = vec![0, 1, -2, 3, 0];
        let raw = xdrop_ungapped_with_cbs(&query, None, &subject, 1, 1, 20, true, &sm);
        let adjusted = xdrop_ungapped_with_cbs(&query, Some(&cbs), &subject, 1, 1, 20, true, &sm);
        assert_eq!(raw.i, 1);
        assert_eq!(raw.j, 1);
        assert_eq!(raw.len, 3);
        assert_eq!(raw.identities, 3);
        assert_eq!(adjusted.score, raw.score + 1 - 2 + 3);
        assert_eq!(adjusted.identities, 3);

        let no_id = xdrop_ungapped_with_cbs(&query, Some(&cbs), &subject, 1, 1, 20, false, &sm);
        assert_eq!(no_id.score, adjusted.score);
        assert_eq!(no_id.identities, 0);
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
    fn test_score_range_s_with_and_without_cbs() {
        let sm = make_test_matrix();
        let query = vec![0, 1, 2, 3];
        let subject = vec![0, 1, 2, 3];
        let raw = score_range_s(&query, None, &subject, 1, 1, 4, &sm);
        assert_eq!(raw.i, 1);
        assert_eq!(raw.j, 1);
        assert_eq!(raw.len, 3);
        assert_eq!(raw.score, score_range(&query, &subject, 1, 1, 4, &sm));

        let cbs = vec![0, 2, -1, 3];
        let adjusted = score_range_s(&query, Some(&cbs), &subject, 1, 1, 4, &sm);
        assert_eq!(adjusted.i, raw.i);
        assert_eq!(adjusted.j, raw.j);
        assert_eq!(adjusted.len, raw.len);
        assert_eq!(adjusted.score, raw.score + 2 - 1 + 3);
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
    fn test_xdrop_anchored_and_anchor_extension() {
        let sm = make_test_matrix();
        let query = vec![0, 1, 2, 3, 4];
        let subject = vec![0, 1, 2, 3, 4];
        let right = xdrop_anchored_right(&query[3..], &subject[3..], 2, 20, &sm);
        assert_eq!(right.1, 2);
        assert!(right.0 > 0);
        let left = xdrop_anchored_left(&query, &subject, 1, 1, 2, 20, &sm);
        assert_eq!(left.1, 2);
        assert!(left.0 > 0);

        let anchor = DiagonalSegment::new(2, 2, 1, score_range(&query, &subject, 2, 2, 3, &sm));
        let extended = xdrop_ungapped_from_anchor(&query, &subject, &anchor, 20, &sm);
        assert_eq!(extended.i, 0);
        assert_eq!(extended.j, 0);
        assert_eq!(extended.len, 5);
        assert_eq!(extended.score, score_range(&query, &subject, 0, 0, 5, &sm));
    }

    #[test]
    fn test_make_clipped_and_null_anchor() {
        let sm = make_test_matrix();
        let query = vec![0, 1, 2, 3, 4, 5];
        let target = vec![9, 1, 2, 3, 4, 9];
        let base_segment = DiagonalSegment::new(
            0,
            0,
            query.len() as i32,
            score_range(&query, &target, 0, 0, query.len(), &sm),
        );
        let anchor = Anchor::new(base_segment, -3, 4, -2, 5, 100);
        let clipped = make_clipped_anchor(&anchor, &query, None, &target, 2, 0.0, &sm);
        assert_eq!(clipped.segment.i, 1);
        assert_eq!(clipped.segment.j, 1);
        assert_eq!(clipped.segment.len, 4);
        assert_eq!(clipped.d_min_left, -3);
        assert_eq!(clipped.d_max_right, 5);
        assert_eq!(
            clipped.prefix_score,
            100 - score_range(&query, &target, 5, 5, 6, &sm)
        );

        let cbs = vec![0, 2, 2, 2, 2, 0];
        let adjusted = make_clipped_anchor(&anchor, &query, Some(&cbs), &target, 2, 0.0, &sm);
        assert_eq!(adjusted.segment.score, clipped.segment.score + 8);

        let null_anchor = make_null_anchor(&anchor);
        assert_eq!(null_anchor.segment.i, 3);
        assert_eq!(null_anchor.segment.j, 3);
        assert_eq!(null_anchor.segment.len, 0);
    }

    #[test]
    fn test_trivial_at_and_trivial() {
        let sm = make_test_matrix();
        let query = vec![0, 1, 2];
        let target = vec![0, 1, 2];
        let hsp = trivial(&query, &target, None, f64::MAX, &sm);
        assert!(hsp.score > 0);
        assert_eq!(hsp.query_range, Interval::new(0, 3));
        assert_eq!(hsp.query_source_range, hsp.query_range);
        assert_eq!(hsp.subject_range, Interval::new(0, 3));
        assert_eq!(hsp.bit_score, sm.bitscore(hsp.score as f64));

        let cbs = vec![2, -1, 3];
        let adjusted = trivial(&query, &target, Some(&cbs), f64::MAX, &sm);
        assert_eq!(adjusted.score, hsp.score + 2 - 1 + 3);

        let filtered = trivial(&query, &target, None, 0.0, &sm);
        assert_eq!(filtered.score, 0);
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
