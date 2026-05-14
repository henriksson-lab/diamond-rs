use crate::basic::value::{Loc, Score};
use crate::dp::ungapped::DiagonalSegment;
use crate::util::interval::Interval;

#[derive(Debug, Clone, PartialEq)]
pub struct Anchor {
    pub segment: DiagonalSegment,
    pub prefix_score: Score,
    pub d_min_left: Loc,
    pub d_max_left: Loc,
    pub d_min_right: Loc,
    pub d_max_right: Loc,
}

impl Default for Anchor {
    fn default() -> Self {
        Self {
            segment: DiagonalSegment::default(),
            prefix_score: 0,
            d_min_left: Loc::MAX,
            d_max_left: Loc::MIN,
            d_min_right: Loc::MAX,
            d_max_right: Loc::MIN,
        }
    }
}

impl Anchor {
    pub fn from_diagonal_segment(d: DiagonalSegment) -> Self {
        Self {
            segment: d,
            ..Self::default()
        }
    }

    pub fn new(
        d: DiagonalSegment,
        d_min_left: Loc,
        d_max_left: Loc,
        d_min_right: Loc,
        d_max_right: Loc,
        prefix_score: Score,
    ) -> Self {
        Self {
            segment: d,
            prefix_score,
            d_min_left,
            d_max_left,
            d_min_right,
            d_max_right,
        }
    }

    pub fn from_coords(query_pos: i32, subject_pos: i32, len: i32, score: i32, ident: Loc) -> Self {
        Self {
            segment: DiagonalSegment::with_identities(query_pos, subject_pos, len, score, ident),
            ..Self::default()
        }
    }

    pub fn assign_diagonal_segment(&mut self, d: DiagonalSegment) -> &mut Self {
        self.segment = d;
        self
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ApproxHsp {
    pub d_min: i32,
    pub d_max: i32,
    pub score: i32,
    pub frame: i32,
    pub query_source_range: Interval,
    pub query_range: Interval,
    pub subject_range: Interval,
    pub evalue: f64,
    pub max_diag: Anchor,
}

impl ApproxHsp {
    pub fn new(frame: u32, score: Score) -> Self {
        Self {
            d_min: i32::MAX,
            d_max: i32::MIN,
            score,
            frame: frame as i32,
            query_source_range: Interval::default(),
            query_range: Interval::default(),
            subject_range: Interval::default(),
            evalue: f64::MAX,
            max_diag: Anchor::default(),
        }
    }

    pub fn from_query_source_range(query_source_range: Interval) -> Self {
        Self {
            query_source_range,
            ..Self::new(0, 0)
        }
    }

    pub fn from_parts(
        d_min: i32,
        d_max: i32,
        score: i32,
        frame: i32,
        query_range: Interval,
        subject_range: Interval,
        max_diag: Anchor,
        evalue: f64,
    ) -> Self {
        Self {
            d_min,
            d_max,
            score,
            frame,
            query_source_range: Interval::default(),
            query_range,
            subject_range,
            evalue,
            max_diag,
        }
    }

    pub fn partial_score_diagonal_segment(&self, d: &DiagonalSegment) -> i32 {
        let overlap = d
            .subject_range()
            .overlap_factor(&self.subject_range)
            .max(d.query_range().overlap_factor(&self.query_range));
        ((1.0 - overlap) * d.score as f64) as i32
    }

    pub fn partial_score_approx_hsp(&self, x: &ApproxHsp) -> i32 {
        let overlap = x
            .subject_range
            .overlap_factor(&self.subject_range)
            .max(x.query_range.overlap_factor(&self.query_range));
        ((1.0 - overlap) * x.score as f64) as i32
    }

    pub fn disjoint_diagonal_segment(&self, d: &DiagonalSegment) -> bool {
        self.query_range.intersect(&d.query_range()).length() == 0
            && self.subject_range.intersect(&d.subject_range()).length() == 0
    }

    pub fn disjoint_approx_hsp(&self, x: &ApproxHsp) -> bool {
        self.query_range.intersect(&x.query_range).length() == 0
            && self.subject_range.intersect(&x.subject_range).length() == 0
    }

    pub fn rel_disjoint_diagonal_segment(&self, d: &DiagonalSegment) -> bool {
        self.query_range.intersect(&d.query_range()).length() == 0
            || self.subject_range.intersect(&d.subject_range()).length() == 0
    }

    pub fn rel_disjoint_approx_hsp(&self, x: &ApproxHsp) -> bool {
        self.query_range.intersect(&x.query_range).length() == 0
            || self.subject_range.intersect(&x.subject_range).length() == 0
    }

    pub fn collinear_approx_hsp(&self, x: &ApproxHsp) -> bool {
        let di = x.query_range.begin - self.query_range.begin;
        let dj = x.subject_range.begin - self.subject_range.begin;
        (di >= 0 && dj >= 0) || (di <= 0 && dj <= 0)
    }

    pub fn collinear_diagonal_segment(&self, d: &DiagonalSegment) -> bool {
        let di = d.i - self.query_range.begin;
        let dj = d.j - self.subject_range.begin;
        (di >= 0 && dj >= 0) || (di <= 0 && dj <= 0)
    }

    pub fn cmp_diag(x: &ApproxHsp, y: &ApproxHsp) -> bool {
        x.frame < y.frame || (x.frame == y.frame && x.d_min < y.d_min)
    }
}

pub struct ApproxHspFrame;

impl ApproxHspFrame {
    pub fn get(x: &ApproxHsp) -> u32 {
        x.frame as u32
    }
}

impl std::fmt::Display for ApproxHsp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "Score={} query_range={} target_range={}",
            self.score, self.query_range, self.subject_range
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_anchor_constructors_and_assign() {
        let d = DiagonalSegment::with_identities(1, 2, 3, 4, 5);
        let mut a = Anchor::from_diagonal_segment(d.clone());
        assert_eq!(a.segment, d);
        assert_eq!(a.prefix_score, 0);
        assert_eq!(a.d_min_left, Loc::MAX);

        let b = Anchor::new(d.clone(), -2, 3, -4, 5, 99);
        assert_eq!(b.prefix_score, 99);
        assert_eq!(b.d_min_right, -4);

        let c = Anchor::from_coords(7, 8, 9, 10, 11);
        assert_eq!(c.segment, DiagonalSegment::with_identities(7, 8, 9, 10, 11));

        a.assign_diagonal_segment(DiagonalSegment::new(10, 20, 5, 6));
        assert_eq!(a.segment, DiagonalSegment::new(10, 20, 5, 6));
    }

    #[test]
    fn test_approx_hsp_ranges_scores_and_relations() {
        let h = ApproxHsp::from_parts(
            -3,
            5,
            100,
            1,
            Interval::new(10, 20),
            Interval::new(30, 40),
            Anchor::default(),
            1.5,
        );
        let d = DiagonalSegment::new(15, 35, 10, 80);
        assert_eq!(h.partial_score_diagonal_segment(&d), 40);
        assert!(!h.disjoint_diagonal_segment(&d));
        assert!(!h.rel_disjoint_diagonal_segment(&d));
        assert!(h.collinear_diagonal_segment(&d));

        let x = ApproxHsp::from_parts(
            1,
            6,
            60,
            1,
            Interval::new(25, 30),
            Interval::new(45, 50),
            Anchor::default(),
            2.0,
        );
        assert_eq!(h.partial_score_approx_hsp(&x), 60);
        assert!(h.disjoint_approx_hsp(&x));
        assert!(h.rel_disjoint_approx_hsp(&x));
        assert!(h.collinear_approx_hsp(&x));

        let y = ApproxHsp::from_parts(
            -5,
            -1,
            10,
            0,
            Interval::new(0, 5),
            Interval::new(100, 105),
            Anchor::default(),
            3.0,
        );
        assert!(ApproxHsp::cmp_diag(&y, &h));
        assert_eq!(ApproxHspFrame::get(&h), 1);
        assert!(format!("{}", h).contains("Score=100"));
    }
}
