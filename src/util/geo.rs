//! Integer diagonal geometry helpers.
//!
//! Direct Rust counterparts of `diamond/src/util/geo/geo.h`.

use crate::basic::translate::{Frame, TranslatedPosition};
use crate::basic::value::Loc;
use crate::dp::ungapped::DiagonalSegment;
use crate::util::interval::Interval;

pub type Hit = (Loc, Loc);

pub fn i(j: i32, d: i32) -> i32 {
    d + j
}

pub fn j(i: i32, d: i32) -> i32 {
    i - d
}

pub fn diag_sub_matrix(d: i32, i0: i32, j0: i32) -> i32 {
    d + j0 - i0
}

pub fn rev_diag(d: i32, qlen: i32, tlen: i32) -> i32 {
    -d + qlen - tlen
}

pub fn min_diag(_qlen: i32, tlen: i32) -> i32 {
    -(tlen - 1)
}

pub fn max_diag(qlen: i32, _tlen: i32) -> i32 {
    qlen - 1
}

pub fn clip_diag(d: i32, qlen: i32, tlen: i32) -> i32 {
    d.max(min_diag(qlen, tlen)).min(max_diag(qlen, tlen))
}

pub fn assert_diag_bounds(d: i32, qlen: i32, tlen: i32) {
    assert!(d >= min_diag(qlen, tlen) && d <= max_diag(qlen, tlen));
}

/// Direct Rust counterpart of C++ DiagonalSegmentT.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DiagonalSegmentT {
    pub i: TranslatedPosition,
    pub j: i32,
    pub len: i32,
    pub score: i32,
}

impl DiagonalSegmentT {
    pub fn new(i: TranslatedPosition, j: i32, len: i32, score: i32) -> Self {
        DiagonalSegmentT { i, j, len, score }
    }

    pub fn from_diagonal_segment(d: &DiagonalSegment, frame: Frame) -> Self {
        DiagonalSegmentT {
            i: TranslatedPosition::new(d.i, frame),
            j: d.j,
            len: d.len,
            score: d.score,
        }
    }

    pub fn subject_last(&self) -> i32 {
        self.j + self.len - 1
    }

    pub fn query_last(&self) -> TranslatedPosition {
        self.i + self.len - 1
    }

    pub fn subject_end(&self) -> i32 {
        self.j + self.len
    }

    pub fn query_end(&self) -> TranslatedPosition {
        self.i + self.len
    }

    pub fn diag(&self) -> i32 {
        self.i.translated - self.j
    }

    pub fn query_absolute_range(&self, dna_len: i32, query_translated: bool) -> Interval {
        TranslatedPosition::absolute_interval(self.i, self.i + self.len, dna_len, query_translated)
    }

    pub fn query_in_strand_range(&self, query_translated: bool) -> Interval {
        Interval::new(
            self.i.in_strand(query_translated),
            (self.i + self.len).in_strand(query_translated),
        )
    }

    pub fn subject_range(&self) -> Interval {
        Interval::new(self.j, self.j + self.len)
    }

    pub fn partial_score(&self, d: &DiagonalSegmentT, query_translated: bool) -> i32 {
        let overlap = self.subject_range().overlap_factor(&d.subject_range()).max(
            self.query_in_strand_range(query_translated)
                .overlap_factor(&d.query_in_strand_range(query_translated)),
        );
        ((1.0 - overlap) * self.score as f64) as i32
    }

    pub fn cut_out(&mut self, d: &DiagonalSegmentT) {
        let ll = (d.i.translated - self.i.translated).min(d.j - self.j);
        let lr = (self.query_end().translated - d.query_end().translated)
            .min(self.subject_end() - d.subject_end());
        let len2;
        if ll > 0 && ll >= lr {
            len2 = self.len.min(ll);
        } else if lr > 0 && lr >= ll {
            len2 = self.len.min(lr);
            self.i = self.query_end() - len2;
            self.j = self.subject_end() - len2;
        } else {
            len2 = 0;
        }
        self.score = (len2 as f64 / self.len as f64 * self.score as f64) as i32;
        self.len = len2;
    }
}

pub fn diag_count(query_len: i32, target_len: i32) -> i32 {
    query_len + target_len - 1
}

pub fn diag_idx(i: i32, j: i32, target_len: i32) -> i32 {
    i - j + target_len - 1
}

pub fn cmp_diag(a: &Hit, b: &Hit) -> bool {
    let d1 = a.0 - a.1;
    let d2 = b.0 - b.1;
    d1 < d2 || (d1 == d2 && a.0 < b.0)
}

pub fn merge_hits(hits: &[Hit], kmer_size: i32, window: Loc, min_len: Loc) -> Vec<DiagonalSegment> {
    let mut v = Vec::new();
    let Some(first) = hits.first() else {
        return v;
    };
    let mut a = first.0;
    let mut b = first.0 + kmer_size;
    let d = first.0 - first.1;
    for hit in &hits[1..] {
        if hit.0 - b < window {
            b = b.max(hit.0 + kmer_size);
        } else {
            if b - a >= min_len {
                v.push(DiagonalSegment::new(a, a - d, b - a, 0));
            }
            a = hit.0;
            b = hit.0 + kmer_size;
        }
    }
    if b - a >= min_len {
        v.push(DiagonalSegment::new(a, a - d, b - a, 0));
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::translate::{Frame, Strand};

    #[test]
    fn test_geo_diagonals() {
        assert_eq!(i(4, 2), 6);
        assert_eq!(j(6, 2), 4);
        assert_eq!(diag_sub_matrix(5, 2, 7), 10);
        assert_eq!(rev_diag(2, 10, 6), 2);
        assert_eq!(clip_diag(-10, 8, 4), -3);
        assert_eq!(clip_diag(10, 8, 4), 7);
    }

    #[test]
    fn test_diagonal_segment_t_ranges_and_partial_score() {
        let frame = Frame::new(Strand::Forward, 1);
        let a = DiagonalSegmentT::new(TranslatedPosition::new(2, frame), 5, 4, 80);
        let b = DiagonalSegmentT::new(TranslatedPosition::new(4, frame), 7, 3, 50);

        assert_eq!(a.subject_last(), 8);
        assert_eq!(a.query_last(), TranslatedPosition::new(5, frame));
        assert_eq!(a.subject_end(), 9);
        assert_eq!(a.query_end(), TranslatedPosition::new(6, frame));
        assert_eq!(a.diag(), -3);
        assert_eq!(a.subject_range(), Interval::new(5, 9));
        assert_eq!(a.query_in_strand_range(true), Interval::new(7, 19));
        assert_eq!(a.query_absolute_range(100, true), Interval::new(7, 19));
        assert_eq!(a.partial_score(&b, true), 40);
    }

    #[test]
    fn test_diagonal_segment_t_cut_out_left_and_right() {
        let frame = Frame::new(Strand::Forward, 0);
        let mut left = DiagonalSegmentT::new(TranslatedPosition::new(10, frame), 20, 10, 100);
        let d = DiagonalSegmentT::new(TranslatedPosition::new(16, frame), 26, 3, 30);
        left.cut_out(&d);
        assert_eq!(left.i, TranslatedPosition::new(10, frame));
        assert_eq!(left.j, 20);
        assert_eq!(left.len, 6);
        assert_eq!(left.score, 60);

        let mut right = DiagonalSegmentT::new(TranslatedPosition::new(10, frame), 20, 10, 100);
        let d = DiagonalSegmentT::new(TranslatedPosition::new(11, frame), 21, 3, 30);
        right.cut_out(&d);
        assert_eq!(right.i, TranslatedPosition::new(14, frame));
        assert_eq!(right.j, 24);
        assert_eq!(right.len, 6);
        assert_eq!(right.score, 60);
    }

    #[test]
    fn test_diag_count_and_idx() {
        assert_eq!(diag_count(4, 6), 9);
        assert_eq!(diag_idx(3, 1, 6), 7);
    }

    #[test]
    fn test_cmp_diag_and_merge_hits() {
        let mut hits = vec![(8, 5), (0, 0), (3, 3), (20, 17), (26, 23)];
        hits.sort_by(|a, b| {
            if cmp_diag(a, b) {
                std::cmp::Ordering::Less
            } else if cmp_diag(b, a) {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });
        assert_eq!(hits, vec![(0, 0), (3, 3), (8, 5), (20, 17), (26, 23)]);

        let segments = merge_hits(&hits[2..], 4, 4, 4);
        assert_eq!(
            segments,
            vec![
                DiagonalSegment::new(8, 5, 4, 0),
                DiagonalSegment::new(20, 17, 10, 0)
            ]
        );
    }
}
