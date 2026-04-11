/// A half-open integer interval [begin, end).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Interval {
    pub begin: i32,
    pub end: i32,
}

impl Interval {
    pub fn new(begin: i32, end: i32) -> Self {
        Interval { begin, end }
    }

    #[inline]
    pub fn length(&self) -> i32 {
        if self.end > self.begin {
            self.end - self.begin
        } else {
            0
        }
    }

    pub fn overlap(&self, rhs: &Interval) -> u32 {
        self.intersect(rhs).length() as u32
    }

    pub fn overlap_factor(&self, rhs: &Interval) -> f64 {
        self.overlap(rhs) as f64 / self.length() as f64
    }

    #[inline]
    pub fn includes(&self, p: i32) -> bool {
        p >= self.begin && p < self.end
    }

    pub fn contains(&self, other: &Interval) -> bool {
        self.begin <= other.begin && self.end >= other.end
    }

    pub fn intersect(&self, rhs: &Interval) -> Interval {
        Interval {
            begin: self.begin.max(rhs.begin),
            end: self.end.min(rhs.end),
        }
    }

    pub fn merge(&mut self, other: &Interval) {
        self.begin = self.begin.min(other.begin);
        self.end = self.end.max(other.end);
    }

    pub fn is_empty(&self) -> bool {
        self.end <= self.begin
    }
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.begin.cmp(&other.begin)
    }
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{};{}]", self.begin, self.end)
    }
}

/// Merge overlapping intervals into disjoint intervals.
pub fn make_disjoint(intervals: &mut [Interval]) -> Vec<Interval> {
    if intervals.is_empty() {
        return Vec::new();
    }
    intervals.sort();
    let mut result = Vec::new();
    let mut a = intervals[0].begin;
    let mut b = intervals[0].end;
    for iv in &intervals[1..] {
        if iv.begin <= b {
            b = b.max(iv.end);
        } else {
            result.push(Interval::new(a, b));
            a = iv.begin;
            b = iv.end;
        }
    }
    result.push(Interval::new(a, b));
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_basic() {
        let i = Interval::new(5, 10);
        assert_eq!(i.length(), 5);
        assert!(i.includes(5));
        assert!(i.includes(9));
        assert!(!i.includes(10));
    }

    #[test]
    fn test_interval_overlap() {
        let a = Interval::new(0, 10);
        let b = Interval::new(5, 15);
        assert_eq!(a.overlap(&b), 5);
    }

    #[test]
    fn test_interval_intersect() {
        let a = Interval::new(0, 10);
        let b = Interval::new(5, 15);
        let c = a.intersect(&b);
        assert_eq!(c, Interval::new(5, 10));
    }

    #[test]
    fn test_make_disjoint() {
        let mut ivs = vec![
            Interval::new(0, 5),
            Interval::new(3, 8),
            Interval::new(10, 15),
        ];
        let result = make_disjoint(&mut ivs);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], Interval::new(0, 8));
        assert_eq!(result[1], Interval::new(10, 15));
    }
}
