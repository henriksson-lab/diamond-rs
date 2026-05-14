/// A half-open integer interval [begin, end).
use std::collections::BTreeMap;

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

    pub fn check(&self, len: i32) -> Result<(), &'static str> {
        if self.begin < 0
            || self.end < 0
            || self.end < self.begin
            || self.begin >= len
            || self.end > len
        {
            Err("interval out of range")
        } else {
            Ok(())
        }
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntervalNode {
    pub count: i64,
    pub min_score: i32,
    pub max_score: i32,
}

impl Default for IntervalNode {
    fn default() -> Self {
        Self {
            count: 0,
            min_score: i32::MAX,
            max_score: 0,
        }
    }
}

impl IntervalNode {
    pub fn new(count: i64, min_score: i32, max_score: i32) -> Self {
        Self {
            count,
            min_score,
            max_score,
        }
    }

    pub fn add(self, score: i32, cap: i64) -> Self {
        Self {
            count: self.count + 1,
            min_score: if self.count < cap {
                self.min_score.min(score)
            } else {
                self.min_score
            },
            max_score: self.max_score.max(score),
        }
    }
}

#[derive(Debug, Clone)]
pub struct IntervalPartition {
    pub cap: i64,
    nodes: BTreeMap<i32, IntervalNode>,
}

#[derive(Debug, Clone, Default)]
pub struct IntervalTree {
    intervals: Vec<Interval>,
}

impl IntervalTree {
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
        }
    }

    pub fn is_overlapped(&self, i: Interval, overlap: f64) -> bool {
        let max = (i.length() as f64 * overlap) as i32;
        let mut n = 0;
        for interval in &self.intervals {
            n += interval.intersect(&i).length();
            if n >= max {
                return true;
            }
        }
        false
    }

    pub fn insert(&mut self, i: Interval) {
        self.intervals.push(i);
    }
}

impl IntervalPartition {
    pub fn new(cap: i64) -> Self {
        let mut nodes = BTreeMap::new();
        nodes.insert(0, IntervalNode::default());
        Self { cap, nodes }
    }

    pub fn insert(&mut self, k: Interval, score: i32) {
        let mut i_key = self
            .nodes
            .range(k.begin..)
            .next()
            .map(|(&key, _)| key)
            .unwrap_or(k.begin);
        if !self.nodes.contains_key(&k.begin) {
            let node = self
                .nodes
                .range(..k.begin)
                .next_back()
                .map(|(_, &node)| node)
                .unwrap_or_default();
            self.nodes.insert(k.begin, node);
            i_key = k.begin;
        }

        let keys: Vec<i32> = self
            .nodes
            .range(i_key..k.end)
            .map(|(&key, _)| key)
            .collect();
        let mut last = IntervalNode::default();
        for key in keys {
            if let Some(node) = self.nodes.get_mut(&key) {
                last = *node;
                *node = node.add(score, self.cap);
            }
        }
        self.nodes.entry(k.end).or_insert(last);
    }

    pub fn covered(&self, k: Interval) -> i32 {
        let mut c = 0;
        for (interval, node) in self.segments_from(k.begin) {
            if interval.begin >= k.end {
                break;
            }
            if node.count >= self.cap {
                c += k.overlap(&interval) as i32;
            }
        }
        c
    }

    pub fn covered_max_score(&self, k: Interval, max_score: i32) -> i32 {
        let mut c = 0;
        for (interval, node) in self.segments_from(k.begin) {
            if interval.begin >= k.end {
                break;
            }
            if node.max_score >= max_score {
                c += k.overlap(&interval) as i32;
            }
        }
        c
    }

    pub fn covered_min_score(&self, k: Interval, min_score: i32) -> i32 {
        let mut c = 0;
        for (interval, node) in self.segments_from(k.begin) {
            if interval.begin >= k.end {
                break;
            }
            if node.count >= self.cap && node.min_score >= min_score {
                c += k.overlap(&interval) as i32;
            }
        }
        c
    }

    pub fn min_score(&self, k: Interval) -> i32 {
        let mut s = i32::MAX;
        for (interval, node) in self.segments_from(k.begin) {
            if interval.begin >= k.end {
                break;
            }
            if node.count < self.cap {
                return 0;
            }
            s = s.min(node.min_score);
        }
        s
    }

    pub fn max_score(&self, k: Interval) -> i32 {
        let mut s = i32::MAX;
        for (interval, node) in self.segments_from(k.begin) {
            if interval.begin >= k.end {
                break;
            }
            s = s.min(node.max_score);
        }
        assert!(s != i32::MAX);
        s
    }

    pub fn segments_from(&self, p: i32) -> Vec<(Interval, IntervalNode)> {
        let mut keys: Vec<i32> = self.nodes.keys().copied().collect();
        if keys.is_empty() {
            return Vec::new();
        }
        let mut idx = match keys.binary_search(&p) {
            Ok(idx) => idx,
            Err(0) => 0,
            Err(idx) => idx - 1,
        };
        let mut out = Vec::new();
        while idx < keys.len() {
            let begin = keys[idx];
            let end = keys.get(idx + 1).copied().unwrap_or(i32::MAX);
            out.push((
                Interval::new(begin, end),
                *self.nodes.get(&begin).expect("partition node"),
            ));
            idx += 1;
        }
        keys.clear();
        out
    }
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
        assert_eq!(i.check(10), Ok(()));
        assert!(Interval::new(-1, 2).check(10).is_err());
        assert!(Interval::new(8, 7).check(10).is_err());
        assert!(Interval::new(10, 10).check(10).is_err());
        assert!(Interval::new(9, 11).check(10).is_err());
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

    #[test]
    fn test_interval_node_add() {
        let node = IntervalNode::default().add(9, 2);
        assert_eq!(node, IntervalNode::new(1, 9, 9));
        let node = node.add(7, 2).add(5, 2);
        assert_eq!(node, IntervalNode::new(3, 7, 9));
    }

    #[test]
    fn test_interval_partition_insert_and_covered() {
        let mut p = IntervalPartition::new(2);
        p.insert(Interval::new(0, 10), 100);
        assert_eq!(p.covered(Interval::new(0, 10)), 0);
        p.insert(Interval::new(5, 15), 80);
        assert_eq!(p.covered(Interval::new(0, 20)), 5);
        assert_eq!(p.covered_max_score(Interval::new(0, 20), 90), 10);
        assert_eq!(p.covered_min_score(Interval::new(0, 20), 80), 5);
        assert_eq!(p.min_score(Interval::new(5, 10)), 80);
        assert_eq!(p.max_score(Interval::new(5, 10)), 100);
    }

    #[test]
    fn test_interval_partition_split_preserves_state() {
        let mut p = IntervalPartition::new(1);
        p.insert(Interval::new(10, 20), 50);
        p.insert(Interval::new(12, 18), 70);
        assert_eq!(p.covered(Interval::new(10, 20)), 10);
        assert_eq!(p.covered_max_score(Interval::new(10, 20), 70), 6);
        assert_eq!(p.max_score(Interval::new(10, 12)), 50);
        assert_eq!(p.max_score(Interval::new(12, 18)), 70);
    }

    #[test]
    fn test_interval_tree_overlap_accumulates_until_threshold() {
        let mut tree = IntervalTree::new();
        tree.insert(Interval::new(0, 5));
        tree.insert(Interval::new(8, 12));
        assert!(tree.is_overlapped(Interval::new(0, 10), 0.5));
        assert!(!tree.is_overlapped(Interval::new(0, 10), 0.8));
        tree.insert(Interval::new(5, 8));
        assert!(tree.is_overlapped(Interval::new(0, 10), 0.8));
    }
}
