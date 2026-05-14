use std::cmp::Ordering;
use std::collections::{BTreeSet, VecDeque};

use crate::util::data_structures::{DoubleArray, FlatArray, FlatIndex};
use crate::util::hash::murmur_hash_u64;

pub const MAX_SHAPE_LEN: usize = 19;
pub const JOIN_SPLIT_SIZE: usize = 100_000;
pub const JOIN_SPLIT_KEY_LEN: u32 = 17;
pub const JOIN_HT_FACTOR: f64 = 1.3;
pub const HASH_JOIN_SWAP: bool = false;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Edge<Int> {
    pub node1: Int,
    pub node2: Int,
    pub weight: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Entry<Int> {
    pub node: Int,
    pub depth: Int,
}

impl<Int> Entry<Int> {
    pub fn new(node: Int, depth: Int) -> Self {
        Self { node, depth }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RelPtr {
    pub r: u32,
    pub s: u32,
}

impl RelPtr {
    pub fn new(r: u32) -> Self {
        Self { r, s: 0 }
    }
}

impl From<RelPtr> for u32 {
    fn from(value: RelPtr) -> Self {
        value.r
    }
}

pub trait AlgoInt:
    Copy
    + Ord
    + Default
    + From<u8>
    + TryFrom<usize>
    + TryInto<usize>
    + std::ops::Add<Output = Self>
    + std::fmt::Debug
    + FlatIndex
where
    <Self as TryFrom<usize>>::Error: std::fmt::Debug,
    <Self as TryInto<usize>>::Error: std::fmt::Debug,
{
    const NIL: Self;
}

impl AlgoInt for u32 {
    const NIL: Self = u32::MAX;
}

impl AlgoInt for u64 {
    const NIL: Self = u64::MAX;
}

impl AlgoInt for usize {
    const NIL: Self = usize::MAX;
}

fn int_to_usize<Int: AlgoInt>(x: Int) -> usize
where
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    x.try_into().unwrap()
}

fn usize_to_int<Int: AlgoInt>(x: usize) -> Int
where
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
{
    Int::try_from(x).unwrap()
}

pub fn neighbor_count2<Int>(edges: &[Edge<Int>], centroids: &[Int]) -> Int
where
    Int: AlgoInt,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut n = 0usize;
    let mut last = Int::NIL;
    for edge in edges {
        if centroids[int_to_usize(edge.node2)] == Int::NIL && edge.node2 != last {
            n += 1;
            last = edge.node2;
        }
    }
    usize_to_int(n)
}

fn neighbor_count<Int>(edges: &mut [Edge<Int>], centroids: &[Int]) -> Int
where
    Int: AlgoInt,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut n = 0usize;
    let mut last = Int::NIL;
    let mut w = 0usize;
    let mut i = 0usize;
    while i < edges.len() {
        if edges[i].node1 == Int::NIL {
            break;
        }
        if centroids[int_to_usize(edges[i].node2)] == Int::NIL && edges[i].node2 != last {
            n += 1;
            last = edges[i].node2;
            if w < i {
                edges.swap(w, i);
            }
            w += 1;
        }
        i += 1;
    }
    if w < edges.len() {
        edges[w].node1 = Int::NIL;
    }
    usize_to_int(n)
}

fn neighbor_count_with_member_counts<Int>(
    node: Int,
    edges: &[Edge<Int>],
    centroids: &[Int],
    member_counts: &[Int],
) -> Int
where
    Int: AlgoInt,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut n = int_to_usize(member_counts[int_to_usize(node)]);
    for edge in edges {
        if centroids[int_to_usize(edge.node2)] == Int::NIL {
            n += int_to_usize(member_counts[int_to_usize(edge.node2)]);
        }
    }
    usize_to_int(n)
}

fn fix_assignment<Int>(centroids: &mut [Int])
where
    Int: AlgoInt,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut i = 0usize;
    while i < centroids.len() {
        let ci = int_to_usize(centroids[i]);
        if centroids[ci] != centroids[i] {
            centroids[i] = centroids[ci];
        } else {
            i += 1;
        }
    }
}

pub fn make_cluster_gvc<Int>(
    rep: Int,
    neighbors: &FlatArray<Edge<Int>, Int>,
    centroids: &mut [Int],
    merge_recursive: bool,
) where
    Int: AlgoInt,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    centroids[int_to_usize(rep)] = rep;
    for edge in neighbors.range(rep) {
        let node2 = int_to_usize(edge.node2);
        if centroids[node2] == Int::NIL || (merge_recursive && centroids[node2] == edge.node2) {
            centroids[node2] = rep;
        }
    }
}

pub fn make_cluster_cc<Int>(
    rep: Int,
    neighbors: &FlatArray<Edge<Int>, Int>,
    centroids: &mut [Int],
    depth: Int,
) where
    Int: AlgoInt,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    centroids[int_to_usize(rep)] = rep;
    let mut q = VecDeque::new();
    for edge in neighbors.range(rep) {
        if centroids[int_to_usize(edge.node2)] == Int::NIL {
            q.push_back(Entry::new(edge.node2, Int::from(1)));
        }
    }
    while let Some(node) = q.pop_front() {
        let node_usize = int_to_usize(node.node);
        if centroids[node_usize] != Int::NIL || node.depth > depth {
            continue;
        }
        for edge in neighbors.range(node.node) {
            if centroids[int_to_usize(edge.node2)] == Int::NIL {
                q.push_back(Entry::new(edge.node2, node.depth + Int::from(1)));
            }
        }
        centroids[node_usize] = rep;
    }
}

pub fn greedy_vertex_cover<Int>(
    neighbors: &mut FlatArray<Edge<Int>, Int>,
    member_counts: Option<&[Int]>,
    merge_recursive: bool,
    reassign: bool,
    connected_component_depth: Int,
) -> Vec<Int>
where
    Int: AlgoInt,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let size = int_to_usize(neighbors.size());
    let mut q = std::collections::BinaryHeap::<(Int, Int)>::new();
    let mut centroids = vec![Int::NIL; size];
    for i in 0..size {
        let node = usize_to_int::<Int>(i);
        let count = if let Some(member_counts) = member_counts {
            neighbor_count_with_member_counts(
                node,
                neighbors.range(node),
                &centroids,
                member_counts,
            )
        } else {
            neighbors.count(node)
        };
        q.push((count, node));
    }

    while let Some((_, node)) = q.pop() {
        let node_usize = int_to_usize(node);
        if centroids[node_usize] != Int::NIL {
            continue;
        }
        let count = if let Some(member_counts) = member_counts {
            neighbor_count_with_member_counts(
                node,
                neighbors.range(node),
                &centroids,
                member_counts,
            )
        } else {
            neighbor_count(neighbors.range_mut(node), &centroids)
        };
        if q.peek().is_some_and(|top| count < top.0) {
            q.push((count, node));
        } else if connected_component_depth > Int::from(0) {
            make_cluster_cc(node, neighbors, &mut centroids, connected_component_depth);
        } else {
            make_cluster_gvc(node, neighbors, &mut centroids, merge_recursive);
        }
    }

    if reassign {
        let mut weights = vec![f64::NEG_INFINITY; size];
        for node in 0..size {
            let node_int = usize_to_int::<Int>(node);
            if centroids[node] == node_int {
                for edge in neighbors.range(node_int) {
                    let node2 = int_to_usize(edge.node2);
                    if centroids[node2] != edge.node2 && edge.weight > weights[node2] {
                        weights[node2] = edge.weight;
                        centroids[node2] = node_int;
                    }
                }
            }
        }
    }

    if merge_recursive {
        fix_assignment(&mut centroids);
    }
    centroids
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PatternMatcher {
    min_len: u32,
    suffix_mask: u32,
    table: Vec<u8>,
}

impl PatternMatcher {
    pub const SIZE: usize = 1usize << MAX_SHAPE_LEN;

    pub fn new(patterns: &[u32]) -> Self {
        let mut min_len = 32u32;
        let mut max_len = 0u32;
        for &pattern in patterns {
            assert!(pattern != 0);
            let len = 32 - pattern.leading_zeros();
            max_len = max_len.max(len);
            min_len = min_len.min(len);
        }
        let suffix_mask = (1u32 << max_len) - 1;
        let mut table = vec![0u8; Self::SIZE];
        for s in 0..=suffix_mask {
            for &pattern in patterns {
                if (s & pattern) == pattern {
                    table[s as usize] = 1;
                }
            }
        }
        Self {
            min_len,
            suffix_mask,
            table,
        }
    }

    pub fn hit(&self, mut h: u32, len: u32) -> u32 {
        if len < self.min_len {
            return 0;
        }
        let mask = self.suffix_mask;
        let end = len - self.min_len + 1;
        let mut r = 0u32;
        for i in 0..end {
            r |= (self.table[(h & mask) as usize] as u32) << i;
            h >>= 1;
        }
        r
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JoinArrayIterator<'a, T> {
    data: &'a [T],
    ptr: usize,
    end: usize,
}

impl<'a, T> JoinArrayIterator<'a, T>
where
    T: Copy + Into<usize> + PartialEq + From<u8>,
{
    pub fn new(data: &'a [T]) -> Self {
        Self {
            data,
            ptr: 0,
            end: data.len(),
        }
    }

    pub fn range(&self) -> &[T] {
        let n = self.data[self.ptr].into();
        &self.data[self.ptr + 1..self.ptr + 1 + n]
    }

    pub fn increment(&mut self) -> &mut Self {
        self.ptr += self.data[self.ptr].into() + 1;
        while self.ptr < self.end && self.data[self.ptr] == T::from(0) {
            self.ptr += self.data[self.ptr + 1].into() + 1;
        }
        self
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct JoinArrayIteratorMut<'a, T> {
    data: &'a mut [T],
    ptr: usize,
    end: usize,
}

impl<'a, T> JoinArrayIteratorMut<'a, T>
where
    T: Copy + Into<usize> + PartialEq + From<u8>,
{
    pub fn new(data: &'a mut [T]) -> Self {
        let end = data.len();
        Self { data, ptr: 0, end }
    }

    pub fn range(&self) -> &[T] {
        let n = self.data[self.ptr].into();
        &self.data[self.ptr + 1..self.ptr + 1 + n]
    }

    pub fn increment(&mut self) -> &mut Self {
        self.ptr += self.data[self.ptr].into() + 1;
        while self.ptr < self.end && self.data[self.ptr] == T::from(0) {
            self.ptr += self.data[self.ptr + 1].into() + 1;
        }
        self
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn erase(&mut self) {
        self.data[self.ptr + 1] = self.data[self.ptr];
        self.data[self.ptr] = T::from(0);
    }
}

#[derive(Clone)]
pub struct TransformIterator<'a, T, F> {
    data: &'a [T],
    index: isize,
    f: F,
}

impl<'a, T, F, R> TransformIterator<'a, T, F>
where
    F: Copy + Fn(&T) -> R,
{
    pub fn new(data: &'a [T], f: F) -> Self {
        Self { data, index: 0, f }
    }

    pub fn get(&self) -> R {
        (self.f)(&self.data[self.index as usize])
    }

    pub fn increment(&mut self) -> &mut Self {
        self.index += 1;
        self
    }

    pub fn decrement(&mut self) -> &mut Self {
        self.index -= 1;
        self
    }

    pub fn add(&self, n: isize) -> Self {
        Self {
            data: self.data,
            index: self.index + n,
            f: self.f,
        }
    }

    pub fn distance(&self, it: &Self) -> isize {
        self.index - it.index
    }

    pub fn not_equal(&self, it: &Self) -> bool {
        self.index != it.index
    }
}

pub fn transform<T, F, R>(data: &[T], f: F) -> TransformIterator<'_, T, F>
where
    F: Copy + Fn(&T) -> R,
{
    TransformIterator::new(data, f)
}

pub fn join_sorted_lists<T1, T2, K1, K2, V, R>(
    a: &[T1],
    b: &[T2],
    key1: K1,
    key2: K2,
    value: V,
) -> Result<Vec<R>, String>
where
    K1: Copy + Fn(&T1) -> String,
    K2: Copy + Fn(&T2) -> String,
    V: Copy + Fn(&T1, &T2) -> R,
{
    let mut out = Vec::new();
    let mut i = 0usize;
    let mut j = 0usize;
    while i < a.len() && j < b.len() {
        let ka = key1(&a[i]);
        let kb = key2(&b[j]);
        if ka < kb {
            i += 1;
        } else if kb < ka {
            j += 1;
        } else {
            out.push(value(&a[i], &b[j]));
            i += 1;
            if i == a.len() {
                break;
            }
            if key2(&b[j]) == key1(&a[i]) {
                continue;
            }
            j += 1;
            if j == b.len() {
                break;
            }
            if ka == key2(&b[j]) {
                return Err(format!("Duplicate keys: {}", ka));
            }
        }
    }
    Ok(out)
}

pub fn merge_keys<T, K, V, KeyType, ValueType>(
    items: &[T],
    key: K,
    value: V,
    begin: KeyType,
    count: usize,
) -> Vec<(KeyType, BTreeSet<ValueType>)>
where
    K: Copy + Fn(&T) -> KeyType,
    V: Copy + Fn(&T) -> ValueType,
    KeyType: Copy + Ord + std::ops::AddAssign<i32>,
    ValueType: Ord,
{
    let mut out = Vec::new();
    let mut current = begin;
    let mut i = 0usize;
    for _ in 0..count {
        let mut set = BTreeSet::new();
        while i < items.len() && key(&items[i]) == current {
            set.insert(value(&items[i]));
            i += 1;
        }
        out.push((current, set));
        current += 1;
    }
    out
}

pub fn first<T1: Copy, T2>(p: &(T1, T2)) -> T1 {
    p.0
}

pub fn second<T1, T2: Copy>(p: &(T1, T2)) -> T2 {
    p.1
}

#[derive(Debug, Clone, PartialEq)]
pub struct HyperLogLog {
    p: i32,
    m: i32,
    registers: Vec<u8>,
    alpha: f64,
}

impl HyperLogLog {
    pub fn new(precision: i32) -> Result<Self, &'static str> {
        if !(4..=20).contains(&precision) {
            return Err("Precision must be between 4 and 20");
        }
        let m = 1 << precision;
        let mut hll = Self {
            p: precision,
            m,
            registers: vec![0; m as usize],
            alpha: 0.0,
        };
        hll.compute_alpha();
        Ok(hll)
    }

    pub fn with_default_precision() -> Self {
        Self::new(10).expect("valid default precision")
    }

    pub fn add(&mut self, x: i64) {
        let hash = murmur_hash_u64(x as u64);
        self.process_hash(hash);
    }

    pub fn estimate(&self) -> i64 {
        let mut sum = 0.0;
        let mut zeros = 0;
        let mut all_zero = true;

        for &r in &self.registers {
            sum += 1.0 / ((1u64 << r) as f64);
            if r != 0 {
                all_zero = false;
            }
            if r == 0 {
                zeros += 1;
            }
        }

        if all_zero {
            return 0;
        }

        let z = 1.0 / sum;
        let mut e = self.alpha * self.m as f64 * self.m as f64 * z;

        if e <= 2.5 * self.m as f64 && zeros > 0 {
            e = self.m as f64 * ((self.m as f64) / zeros as f64).ln();
        }

        e.round() as i64
    }

    pub fn merge(&mut self, other: &HyperLogLog) -> Result<(), &'static str> {
        if self.p != other.p {
            return Err("Precision must match for merging");
        }
        for i in 0..self.m as usize {
            if self.registers[i] < other.registers[i] {
                self.registers[i] = other.registers[i];
            }
        }
        Ok(())
    }

    fn compute_alpha(&mut self) {
        self.alpha = if self.m == 16 {
            0.673
        } else if self.m == 32 {
            0.697
        } else if self.m == 64 {
            0.709
        } else {
            0.7213 / (1.0 + 1.079 / self.m as f64)
        };
    }

    fn process_hash(&mut self, hash: u64) {
        let index = (hash >> (64 - self.p)) as u32;
        let w = hash & ((1u64 << (64 - self.p)) - 1);
        let rho_val = if w == 0 {
            (64 - self.p + 1) as u8
        } else {
            (w.leading_zeros() as i32 - self.p + 1) as u8
        };
        if self.registers[index as usize] < rho_val {
            self.registers[index as usize] = rho_val;
        }
    }
}

impl<Int> Edge<Int> {
    pub fn new(node1: Int, node2: Int, weight: f64) -> Self {
        Self {
            node1,
            node2,
            weight,
        }
    }
}

impl<Int: Ord> Eq for Edge<Int> {}

impl<Int: Ord> PartialOrd for Edge<Int> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<Int: Ord> Ord for Edge<Int> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.node1
            .cmp(&other.node1)
            .then_with(|| self.node2.cmp(&other.node2))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Partition {
    pub parts: usize,
    items: usize,
    size: usize,
    remainder: usize,
}

impl Partition {
    pub fn new(items: usize, parts: usize) -> Self {
        let parts = parts.min(items);
        let (size, remainder) = if parts > 0 {
            (items / parts, items % parts)
        } else {
            (0, 0)
        };
        Self {
            parts,
            items,
            size,
            remainder,
        }
    }

    pub fn begin(&self, i: usize) -> usize {
        let b = i.min(self.remainder);
        b * (self.size + 1) + (i - b) * self.size
    }

    pub fn end(&self, i: usize) -> usize {
        self.begin(i) + self.size(i)
    }

    pub fn size(&self, i: usize) -> usize {
        if i < self.remainder {
            self.size + 1
        } else {
            self.size
        }
    }

    pub fn items(&self) -> usize {
        self.items
    }
}

impl Default for Partition {
    fn default() -> Self {
        Self {
            parts: 0,
            items: 0,
            size: 0,
            remainder: 0,
        }
    }
}

pub fn batch_binary_search<T, Cmp>(q: &[T], t: &[T], out: &mut Vec<isize>, cmp: Cmp, ti: isize)
where
    Cmp: Copy + Fn(&T, &T) -> bool,
{
    let nq = q.len();
    if nq == 0 {
        return;
    }
    if nq == 1 {
        let ub = t.partition_point(|x| !cmp(&q[0], x));
        out.push(ub as isize + ti - 1);
        return;
    }

    let nt = t.len();
    if nt == 1 {
        for _ in q {
            out.push(ti);
        }
        return;
    }

    let d = nt / 2;
    let mid_t = &t[d];
    let mid_q = q.partition_point(|x| cmp(x, mid_t));
    batch_binary_search(&q[..mid_q], &t[..d], out, cmp, ti);
    batch_binary_search(&q[mid_q..], &t[d..], out, cmp, ti + d as isize);
}

pub fn write_varuint32(x: u32, out: &mut Vec<u8>) {
    if x < 1 << 7 {
        out.push((x << 1 | 1) as u8);
    } else if x < 1 << 14 {
        out.extend_from_slice(&((x << 2 | 2) as u16).to_le_bytes());
    } else if x < 1 << 21 {
        out.push(((x & 31) << 3 | 4) as u8);
        out.extend_from_slice(&((x >> 5) as u16).to_le_bytes());
    } else if x < 1 << 28 {
        out.extend_from_slice(&(x << 4 | 8).to_le_bytes());
    } else {
        out.push(((x & 7) << 5 | 16) as u8);
        out.extend_from_slice(&(x >> 3).to_le_bytes());
    }
}

pub fn read_varuint32(ptr: &[u8]) -> Result<(u32, usize), String> {
    if ptr.is_empty() {
        return Err("Format error: Invalid varint encoding.".to_string());
    }
    let b0 = ptr[0];
    let c = b0.trailing_zeros();
    match c {
        0 => Ok(((b0 >> 1) as u32, 1)),
        1 => {
            if ptr.len() < 2 {
                return Err("Format error: Invalid varint encoding.".to_string());
            }
            Ok((((ptr[1] as u32) << 6) | ((b0 as u32) >> 2), 2))
        }
        2 => {
            if ptr.len() < 3 {
                return Err("Format error: Invalid varint encoding.".to_string());
            }
            let b2 = u16::from_le_bytes([ptr[1], ptr[2]]) as u32;
            Ok(((b2 << 5) | ((b0 as u32) >> 3), 3))
        }
        3 => {
            if ptr.len() < 4 {
                return Err("Format error: Invalid varint encoding.".to_string());
            }
            let b2 = u16::from_le_bytes([ptr[2], ptr[3]]) as u32;
            Ok(((b2 << 12) | ((ptr[1] as u32) << 4) | ((b0 as u32) >> 4), 4))
        }
        4 => {
            if ptr.len() < 5 {
                return Err("Format error: Invalid varint encoding.".to_string());
            }
            let b3 = u32::from_le_bytes([ptr[1], ptr[2], ptr[3], ptr[4]]);
            Ok(((b3 << 3) | ((b0 as u32) >> 5), 5))
        }
        _ => Err("Format error: Invalid varint encoding.".to_string()),
    }
}

pub fn merge_capped<T: Ord + Copy>(i: &[T], j: &[T], cap: usize, out: &mut Vec<T>) -> usize {
    let m = cap as isize;
    let mut n = 0isize;
    let mut count = 0usize;
    let mut i0 = 0usize;
    let mut j0 = 0usize;

    while n < m {
        if i0 == i.len() {
            let d = ((m - n) as usize).min(j.len() - j0);
            out.extend_from_slice(&j[j0..j0 + d]);
            return count + d;
        }
        if j0 == j.len() {
            let d = ((m - n) as usize).min(i.len() - i0);
            out.extend_from_slice(&i[i0..i0 + d]);
            return count;
        }
        if i[i0] < j[j0] {
            out.push(i[i0]);
            i0 += 1;
        } else {
            out.push(j[j0]);
            j0 += 1;
            count += 1;
        }
        n += 1;
    }
    count
}

pub fn sort_by_value<T>(slice: &[T], _threads: i32) -> Vec<(T, i64)>
where
    T: Copy + Ord,
{
    let mut out = Vec::with_capacity(slice.len());
    for (i, &value) in slice.iter().enumerate() {
        out.push((value, i as i64));
    }
    out.sort();
    out
}

pub fn insertion_sort_by<T, Compare>(slice: &mut [T], comp: Compare)
where
    T: Clone,
    Compare: Copy + Fn(&T, &T) -> bool,
{
    let n = slice.len();
    if n < 2 {
        return;
    }

    let mut min_it = 0usize;
    for it in 1..n {
        if comp(&slice[it], &slice[min_it]) {
            min_it = it;
        }
    }
    if min_it != 0 {
        slice.swap(0, min_it);
    }

    for i in 1..n {
        let val = slice[i].clone();
        let mut j = i;
        while comp(&val, &slice[j - 1]) {
            slice[j] = slice[j - 1].clone();
            j -= 1;
        }
        slice[j] = val;
    }
}

pub fn insertion_sort<T>(slice: &mut [T])
where
    T: Clone + Ord,
{
    insertion_sort_by(slice, |a, b| a < b);
}

pub fn merge_sort_by<T, Compare>(slice: &mut [T], n_threads: u32, cmp: Compare, level: u32)
where
    T: Clone,
    Compare: Copy + Fn(&T, &T) -> bool,
{
    let diff = slice.len();
    if diff <= 1 {
        return;
    }

    if (1u32 << level) >= n_threads {
        slice.sort_by(|a, b| {
            if cmp(a, b) {
                Ordering::Less
            } else if cmp(b, a) {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        });
        return;
    }

    let mid = diff / 2;
    merge_sort_by(&mut slice[..mid], n_threads, cmp, level + 1);
    merge_sort_by(&mut slice[mid..], n_threads, cmp, level + 1);

    let left = slice[..mid].to_vec();
    let right = slice[mid..].to_vec();
    let mut i = 0usize;
    let mut j = 0usize;
    let mut k = 0usize;
    while i < left.len() && j < right.len() {
        if cmp(&right[j], &left[i]) {
            slice[k] = right[j].clone();
            j += 1;
        } else {
            slice[k] = left[i].clone();
            i += 1;
        }
        k += 1;
    }
    while i < left.len() {
        slice[k] = left[i].clone();
        i += 1;
        k += 1;
    }
    while j < right.len() {
        slice[k] = right[j].clone();
        j += 1;
        k += 1;
    }
}

pub fn merge_sort<T>(slice: &mut [T], n_threads: u32)
where
    T: Clone + Ord,
{
    merge_sort_by(slice, n_threads, |a, b| a < b, 0);
}

pub fn partition_table<T, Key, F>(slice: &[T], n: usize, key: F) -> Vec<usize>
where
    Key: Eq,
    F: Copy + Fn(&T) -> Key,
{
    let mut v = Vec::new();
    let count = slice.len();
    if count == 0 {
        return v;
    }
    let p = Partition::new(count, n);
    v.reserve(p.parts);
    let mut e = 0usize;
    v.push(e);
    for i in 0..p.parts {
        let mut it = p.end(i);
        if it <= e {
            continue;
        }
        let k = key(&slice[it - 1]);
        while it < slice.len() && key(&slice[it]) == k {
            it += 1;
        }
        v.push(it);
        e = it;
    }
    v
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Relation<'a, T> {
    pub data: &'a [T],
    pub n: usize,
}

impl<'a, T> Relation<'a, T> {
    pub fn new(data: &'a [T], n: usize) -> Self {
        Self {
            data: &data[..n],
            n,
        }
    }

    pub fn end(&self) -> &'a [T] {
        &self.data[self.n..]
    }

    pub fn part(&self, begin: usize, n: usize) -> Self {
        Relation::new(&self.data[begin..], n)
    }
}

impl<T> std::ops::Index<usize> for Relation<'_, T> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        &self.data[i]
    }
}

pub trait HashJoinValue: Copy {
    fn write_ne_bytes(self, out: &mut Vec<u8>);
}

macro_rules! impl_hash_join_value {
    ($($t:ty),* $(,)?) => {
        $(
            impl HashJoinValue for $t {
                fn write_ne_bytes(self, out: &mut Vec<u8>) {
                    out.extend_from_slice(&self.to_ne_bytes());
                }
            }
        )*
    };
}

impl_hash_join_value!(u8, i8, u16, i16, u32, i32, u64, i64, usize, isize);

pub trait HashJoinRecord: Copy {
    type Key: Copy + Ord + From<u8> + TryFrom<usize> + TryInto<usize>;
    type Value: HashJoinValue;

    fn key(&self) -> Self::Key;
    fn value(&self) -> Self::Value;
}

pub fn hash_table_join<T>(
    r: Relation<'_, T>,
    s: Relation<'_, T>,
    shift: u32,
    dst_r: &mut DoubleArray,
    dst_s: &mut DoubleArray,
) where
    T: HashJoinRecord,
    T::Key: RadixKey,
    <T::Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <T::Key as TryInto<usize>>::Error: std::fmt::Debug,
{
    if r.n == 0 || s.n == 0 {
        return;
    }
    let n = next_pow2_usize((r.n as f64 * JOIN_HT_FACTOR) as usize);
    let key = ExtractBits::new(T::Key::try_from(n).unwrap(), shift);
    let mut table = vec![(T::Key::from(0), RelPtr::default()); n];
    let mut r_slots = Vec::with_capacity(r.n);
    for item in r.data {
        let k = key.get(item.key());
        let mut p = k.try_into().unwrap();
        let mut wrapped = false;
        loop {
            if table[p].0 == k || table[p].1 == RelPtr::default() {
                table[p].0 = k;
                table[p].1.r += 1;
                r_slots.push(p);
                break;
            }
            p += 1;
            if p == n {
                assert!(!wrapped, "Hash table overflow.");
                p = 0;
                wrapped = true;
            }
        }
    }
    let mut s_hits = Vec::new();
    for item in s.data {
        let k = key.get(item.key());
        let mut p = k.try_into().unwrap();
        let mut wrapped = false;
        loop {
            if table[p].1 == RelPtr::default() {
                break;
            }
            if table[p].0 == k {
                table[p].1.s += 1;
                s_hits.push((p, item.value()));
                break;
            }
            p += 1;
            if p == n {
                assert!(!wrapped, "Hash table overflow.");
                p = 0;
                wrapped = true;
            }
        }
    }
    let mut data_r = Vec::new();
    let mut data_s = Vec::new();
    for (p, (_, ptr)) in table.iter().enumerate().filter(|(_, (_, ptr))| ptr.s != 0) {
        data_r.extend_from_slice(&ptr.r.to_ne_bytes());
        data_s.extend_from_slice(&ptr.s.to_ne_bytes());
        for (slot, item) in r_slots.iter().zip(r.data.iter()) {
            if *slot == p {
                item.value().write_ne_bytes(&mut data_r);
            }
        }
        for (slot, value) in &s_hits {
            if *slot == p {
                value.write_ne_bytes(&mut data_s);
            }
        }
    }
    let size_r = data_r.len();
    let size_s = data_s.len();
    *dst_r = DoubleArray::from_bytes(data_r, size_r);
    *dst_s = DoubleArray::from_bytes(data_s, size_s);
}

pub fn table_join<T>(
    r: Relation<'_, T>,
    s: Relation<'_, T>,
    total_bits: u32,
    shift: u32,
    dst_r: &mut DoubleArray,
    dst_s: &mut DoubleArray,
) where
    T: HashJoinRecord,
    T::Key: RadixKey,
    <T::Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <T::Key as TryInto<usize>>::Error: std::fmt::Debug,
{
    let keys = 1usize << (total_bits - shift);
    let key = ExtractBits::new(T::Key::try_from(keys).unwrap(), shift);
    let mut table = vec![RelPtr::default(); keys];
    for item in r.data {
        let p: usize = key.get(item.key()).try_into().unwrap();
        table[p].r += 1;
    }
    let mut s_hits = Vec::new();
    for item in s.data {
        let p: usize = key.get(item.key()).try_into().unwrap();
        if table[p].r != 0 {
            table[p].s += 1;
            s_hits.push((p, item.value()));
        }
    }
    let mut data_r = Vec::new();
    let mut data_s = Vec::new();
    for (p, ptr) in table.iter().enumerate().filter(|(_, ptr)| ptr.s != 0) {
        data_r.extend_from_slice(&ptr.r.to_ne_bytes());
        data_s.extend_from_slice(&ptr.s.to_ne_bytes());
        for item in r.data {
            if key.get(item.key()).try_into().unwrap() == p {
                item.value().write_ne_bytes(&mut data_r);
            }
        }
        for (slot, value) in &s_hits {
            if *slot == p {
                value.write_ne_bytes(&mut data_s);
            }
        }
    }
    let size_r = data_r.len();
    let size_s = data_s.len();
    *dst_r = DoubleArray::from_bytes(data_r, size_r);
    *dst_s = DoubleArray::from_bytes(data_s, size_s);
}

pub fn hash_join_into<T>(
    r: Relation<'_, T>,
    s: Relation<'_, T>,
    out_r: &mut DoubleArray,
    out_s: &mut DoubleArray,
    total_bits: u32,
    shift: u32,
) where
    T: HashJoinRecord,
    T::Key: RadixKey,
    <T::Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <T::Key as TryInto<usize>>::Error: std::fmt::Debug,
{
    if r.n == 0 || s.n == 0 {
        return;
    }
    let key_bits = total_bits - shift;
    let mut tmp_r = DoubleArray::new();
    let mut tmp_s = DoubleArray::new();
    if r.n < JOIN_SPLIT_SIZE || key_bits < JOIN_SPLIT_KEY_LEN {
        if next_pow2_usize((r.n as f64 * JOIN_HT_FACTOR) as usize) < (1usize << key_bits) {
            hash_table_join(r, s, shift, &mut tmp_r, &mut tmp_s);
        } else {
            table_join(r, s, total_bits, shift, &mut tmp_r, &mut tmp_s);
        }
        out_r.append(&tmp_r);
        out_s.append(&tmp_s);
    } else {
        table_join(r, s, total_bits, shift, &mut tmp_r, &mut tmp_s);
        out_r.append(&tmp_r);
        out_s.append(&tmp_s);
    }
}

pub fn hash_join<'a, T>(
    mut r: Relation<'a, T>,
    mut s: Relation<'a, T>,
    total_bits: u32,
) -> (DoubleArray, DoubleArray)
where
    T: HashJoinRecord,
    T::Key: RadixKey,
    <T::Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <T::Key as TryInto<usize>>::Error: std::fmt::Debug,
{
    let swap = HASH_JOIN_SWAP && r.n > s.n;
    if swap {
        std::mem::swap(&mut r, &mut s);
    }
    let mut out_r = DoubleArray::new();
    let mut out_s = DoubleArray::new();
    hash_join_into(r, s, &mut out_r, &mut out_s, total_bits, 0);
    if swap {
        (out_s, out_r)
    } else {
        (out_r, out_s)
    }
}

fn next_pow2_usize(x: usize) -> usize {
    if x <= 1 {
        1
    } else {
        x.next_power_of_two()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExtractBits<Key> {
    pub shift: u32,
    pub mask: Key,
}

impl<Key> ExtractBits<Key>
where
    Key: Copy + From<u8> + std::ops::Sub<Output = Key>,
{
    pub fn new(n: Key, shift: u32) -> Self {
        Self {
            shift,
            mask: n - Key::from(1),
        }
    }
}

impl<Key> ExtractBits<Key>
where
    Key: Copy + std::ops::Shr<u32, Output = Key> + std::ops::BitAnd<Output = Key>,
{
    pub fn get(&self, x: Key) -> Key {
        (x >> self.shift) & self.mask
    }
}

pub trait RadixKey:
    Copy
    + From<u8>
    + TryFrom<usize>
    + TryInto<usize>
    + std::ops::Sub<Output = Self>
    + std::ops::Shr<u32, Output = Self>
    + std::ops::BitAnd<Output = Self>
{
}

impl<T> RadixKey for T where
    T: Copy
        + From<u8>
        + TryFrom<usize>
        + TryInto<usize>
        + std::ops::Sub<Output = T>
        + std::ops::Shr<u32, Output = T>
        + std::ops::BitAnd<Output = T>
{
}

pub fn radix_cluster<T, Key, GetKey>(
    input: Relation<'_, T>,
    shift: u32,
    out: &mut [T],
    hst: &mut [usize],
    get_key: GetKey,
    radix_bits: u32,
    radix_cluster_buffered: bool,
) where
    T: Copy,
    Key: RadixKey,
    <Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <Key as TryInto<usize>>::Error: std::fmt::Debug,
    GetKey: Copy + Fn(&T) -> Key,
{
    const BUF_SIZE: usize = 8;
    let clusters = 1usize << radix_bits;
    assert!(hst.len() >= clusters);
    let radix = ExtractBits::new(Key::try_from(clusters).unwrap(), shift);

    hst[..clusters].fill(0);
    for item in input.data {
        hst[radix.get(get_key(item)).try_into().unwrap()] += 1;
    }
    let mut sum = 0usize;
    for slot in hst.iter_mut().take(clusters) {
        let c = *slot;
        *slot = sum;
        sum += c;
    }

    if radix_cluster_buffered {
        let mut buf: Vec<Vec<T>> = (0..clusters)
            .map(|_| Vec::with_capacity(BUF_SIZE))
            .collect();
        for item in input.data {
            let r: usize = radix.get(get_key(item)).try_into().unwrap();
            buf[r].push(*item);
            if buf[r].len() == BUF_SIZE {
                let begin = hst[r];
                out[begin..begin + BUF_SIZE].copy_from_slice(&buf[r]);
                hst[r] += BUF_SIZE;
                buf[r].clear();
            }
        }
        for r in 0..clusters {
            let begin = hst[r];
            out[begin..begin + buf[r].len()].copy_from_slice(&buf[r]);
            hst[r] += buf[r].len();
        }
    } else {
        for item in input.data {
            let r: usize = radix.get(get_key(item)).try_into().unwrap();
            out[hst[r]] = *item;
            hst[r] += 1;
        }
    }
}

pub fn parallel_radix_cluster_build_hst<T, Key, GetKey>(
    input: Relation<'_, T>,
    shift: u32,
    hst: &mut [usize],
    get_key: GetKey,
    radix_bits: u32,
) where
    Key: RadixKey,
    <Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <Key as TryInto<usize>>::Error: std::fmt::Debug,
    GetKey: Copy + Fn(&T) -> Key,
{
    let clusters = 1usize << radix_bits;
    let radix = ExtractBits::new(Key::try_from(clusters).unwrap(), shift);
    for item in input.data {
        hst[radix.get(get_key(item)).try_into().unwrap()] += 1;
    }
}

pub fn parallel_radix_cluster_scatter<T, Key, GetKey>(
    input: Relation<'_, T>,
    shift: u32,
    hst: &mut [usize],
    out: &mut [T],
    get_key: GetKey,
    radix_bits: u32,
) where
    T: Copy,
    Key: RadixKey,
    <Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <Key as TryInto<usize>>::Error: std::fmt::Debug,
    GetKey: Copy + Fn(&T) -> Key,
{
    let clusters = 1usize << radix_bits;
    let radix = ExtractBits::new(Key::try_from(clusters).unwrap(), shift);
    for item in input.data {
        let r: usize = radix.get(get_key(item)).try_into().unwrap();
        out[hst[r]] = *item;
        hst[r] += 1;
    }
}

pub fn parallel_radix_cluster<T, Key, GetKey>(
    input: Relation<'_, T>,
    shift: u32,
    out: &mut [T],
    thread_count: usize,
    get_key: GetKey,
    radix_bits: u32,
) where
    T: Copy,
    Key: RadixKey,
    <Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <Key as TryInto<usize>>::Error: std::fmt::Debug,
    GetKey: Copy + Fn(&T) -> Key,
{
    let clusters = 1usize << radix_bits;
    let mut hst = vec![0usize; clusters];
    let p = Partition::new(input.n, thread_count);
    let nt = p.parts;

    let mut thread_hst = Vec::with_capacity(nt);
    for i in 0..nt {
        let mut local = vec![0usize; clusters];
        parallel_radix_cluster_build_hst(
            input.part(p.begin(i), p.size(i)),
            shift,
            &mut local,
            get_key,
            radix_bits,
        );
        for j in 0..clusters {
            hst[j] += local[j];
        }
        thread_hst.push(local);
    }

    let mut sum = 0usize;
    for i in 0..clusters {
        let c = hst[i];
        hst[i] = sum;
        sum += c;

        let mut sum2 = 0usize;
        for local in thread_hst.iter_mut().take(nt) {
            let d = local[i];
            local[i] = hst[i] + sum2;
            sum2 += d;
        }
    }

    for (i, local) in thread_hst.iter_mut().enumerate().take(nt) {
        parallel_radix_cluster_scatter(
            input.part(p.begin(i), p.size(i)),
            shift,
            local,
            out,
            get_key,
            radix_bits,
        );
    }
}

pub fn radix_sort<T, Key, GetKey>(
    slice: &mut [T],
    max_key: u32,
    threads: usize,
    get_key: GetKey,
    radix_bits: u32,
) where
    T: Copy,
    Key: RadixKey,
    <Key as TryFrom<usize>>::Error: std::fmt::Debug,
    <Key as TryInto<usize>>::Error: std::fmt::Debug,
    GetKey: Copy + Fn(&T) -> Key,
{
    let n = slice.len();
    if n <= 1 {
        return;
    }
    let bit_len = if max_key == 0 {
        0
    } else {
        u32::BITS - max_key.leading_zeros()
    };
    let rounds = bit_len.div_ceil(radix_bits);
    if rounds == 0 {
        return;
    }

    let mut in_buf = slice.to_vec();
    let mut out_buf = vec![slice[0]; n];
    for i in 0..rounds {
        if threads > 1 {
            parallel_radix_cluster(
                Relation::new(&in_buf, n),
                i * radix_bits,
                &mut out_buf,
                threads,
                get_key,
                radix_bits,
            );
        } else {
            let mut hst = vec![0usize; 1usize << radix_bits];
            radix_cluster(
                Relation::new(&in_buf, n),
                i * radix_bits,
                &mut out_buf,
                &mut hst,
                get_key,
                radix_bits,
                false,
            );
        }
        std::mem::swap(&mut in_buf, &mut out_buf);
    }
    slice.copy_from_slice(&in_buf);
}

pub fn merge_sorted_files<T, I, F>(begin: impl IntoIterator<Item = I>, mut f: F)
where
    T: Ord,
    I: IntoIterator<Item = T>,
    F: FnMut(T),
{
    let mut d: Vec<_> = begin.into_iter().map(|i| i.into_iter()).collect();
    let mut q: Vec<Option<T>> = d.iter_mut().map(|i| i.next()).collect();
    loop {
        let mut idx: Option<usize> = None;
        for i in 0..q.len() {
            if q[i].is_none() {
                continue;
            }
            if let Some(j) = idx {
                if q[i].as_ref().unwrap() < q[j].as_ref().unwrap() {
                    idx = Some(i);
                }
            } else {
                idx = Some(i);
            }
        }
        let Some(i) = idx else {
            break;
        };
        f(q[i].take().unwrap());
        q[i] = d[i].next();
    }
}

pub fn alloc_size<T>(_x: &T) -> usize {
    std::mem::size_of::<T>()
}

pub fn alloc_size_string_u32_pair(x: &(String, u32)) -> usize {
    4 + std::mem::size_of::<String>() + x.0.len()
}

pub fn external_sorter_default_cmp<Type: Ord>(a: &Type, b: &Type) -> bool {
    a < b
}

#[derive(Debug, Clone)]
pub struct ExternalSorter<Type, Cmp = fn(&Type, &Type) -> bool> {
    cmp: Cmp,
    count: usize,
    size: usize,
    files: Vec<Vec<Type>>,
    buf: Vec<Type>,
    queue: Vec<Option<Type>>,
    heads: Vec<usize>,
}

impl<Type> Default for ExternalSorter<Type, fn(&Type, &Type) -> bool>
where
    Type: Clone + Ord,
{
    fn default() -> Self {
        Self {
            cmp: external_sorter_default_cmp::<Type>,
            count: 0,
            size: 0,
            files: Vec::new(),
            buf: Vec::new(),
            queue: Vec::new(),
            heads: Vec::new(),
        }
    }
}

impl<Type, Cmp> ExternalSorter<Type, Cmp>
where
    Type: Clone,
    Cmp: Copy + Fn(&Type, &Type) -> bool,
{
    pub const BUCKET_SIZE: usize = 2 * 1024 * 1024 * 1024;

    pub fn new(cmp: Cmp) -> Self {
        Self {
            cmp,
            count: 0,
            size: 0,
            files: Vec::new(),
            buf: Vec::new(),
            queue: Vec::new(),
            heads: Vec::new(),
        }
    }

    pub fn push(&mut self, x: Type) {
        self.count += 1;
        self.size += alloc_size(&x);
        self.buf.push(x);
        if self.size > Self::BUCKET_SIZE {
            self.flush();
        }
    }

    pub fn init_read(&mut self) {
        self.flush();
        self.queue.clear();
        self.heads = vec![0; self.files.len()];
        for i in 0..self.files.len() {
            let entry = self.get_entry(i);
            self.queue.push(entry);
        }
    }

    pub fn good(&self) -> bool {
        self.queue.iter().any(|x| x.is_some())
    }

    pub fn current(&self) -> Option<&Type> {
        self.queue_index()
            .and_then(|i| self.queue.get(i).and_then(|x| x.as_ref()))
    }

    pub fn increment(&mut self) {
        let Some(b) = self.queue_index() else {
            return;
        };
        self.queue[b] = self.get_entry(b);
    }

    pub fn count(&self) -> usize {
        self.count
    }

    fn queue_index(&self) -> Option<usize> {
        let mut out: Option<usize> = None;
        for i in 0..self.queue.len() {
            if self.queue[i].is_none() {
                continue;
            }
            if let Some(j) = out {
                if (self.cmp)(
                    self.queue[i].as_ref().unwrap(),
                    self.queue[j].as_ref().unwrap(),
                ) {
                    out = Some(i);
                }
            } else {
                out = Some(i);
            }
        }
        out
    }

    fn get_entry(&mut self, bucket: usize) -> Option<Type> {
        if self.heads[bucket] == self.files[bucket].len() {
            return None;
        }
        let out = self.files[bucket][self.heads[bucket]].clone();
        self.heads[bucket] += 1;
        Some(out)
    }

    fn flush(&mut self) {
        if self.buf.is_empty() {
            return;
        }
        let mut idx: Vec<usize> = (0..self.buf.len()).collect();
        idx.sort_by(|&x, &y| {
            if (self.cmp)(&self.buf[x], &self.buf[y]) {
                Ordering::Less
            } else if (self.cmp)(&self.buf[y], &self.buf[x]) {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        });
        let bucket = idx.into_iter().map(|i| self.buf[i].clone()).collect();
        self.files.push(bucket);
        self.buf.clear();
        self.size = 0;
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Dsu {
    pub parent: Vec<usize>,
    pub rank: Vec<usize>,
}

impl Dsu {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    pub fn find(&mut self, mut x: usize) -> usize {
        while x != self.parent[x] {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    pub fn unite(&mut self, a: usize, b: usize) {
        let mut a = self.find(a);
        let mut b = self.find(b);
        if a != b {
            if self.rank[a] < self.rank[b] {
                std::mem::swap(&mut a, &mut b);
            }
            self.parent[b] = a;
            if self.rank[a] == self.rank[b] {
                self.rank[a] += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edge_ordering() {
        let mut edges = vec![
            Edge::new(2, 1, 0.5),
            Edge::new(1, 9, 1.0),
            Edge::new(1, 3, 2.0),
        ];
        edges.sort();
        assert_eq!(edges[0].node2, 3);
        assert_eq!(edges[1].node2, 9);
        assert_eq!(edges[2].node1, 2);
    }

    #[test]
    fn test_partition() {
        let p = Partition::new(10, 3);
        assert_eq!(p.parts, 3);
        assert_eq!(p.items(), 10);
        assert_eq!(p.begin(0), 0);
        assert_eq!(p.end(0), 4);
        assert_eq!(p.begin(1), 4);
        assert_eq!(p.end(1), 7);
        assert_eq!(p.begin(2), 7);
        assert_eq!(p.end(2), 10);
        assert_eq!(Partition::new(0, 5).parts, 0);
    }

    #[test]
    fn test_pattern_matcher() {
        let matcher = PatternMatcher::new(&[0b11, 0b101]);
        assert_eq!(matcher.hit(0b0011, 4), 0b001);
        assert_eq!(matcher.hit(0b1010, 4), 0b010);
        assert_eq!(matcher.hit(0b0001, 1), 0);
    }

    #[test]
    fn test_join_array_iterator() {
        let data = [2usize, 10, 11, 0, 3, 1, 20, 1, 30];
        let mut it = JoinArrayIterator::new(&data);
        assert!(it.good());
        assert_eq!(it.range(), &[10, 11]);
        it.increment();
        assert_eq!(it.range(), &[30]);

        let mut data = [2usize, 10, 11, 1, 20];
        let mut it = JoinArrayIteratorMut::new(&mut data);
        it.erase();
        assert_eq!(it.range(), &[]);
        it.increment();
        assert_eq!(it.range(), &[11, 1]);
        assert_eq!(data, [0, 2, 11, 1, 20]);
    }

    #[test]
    fn test_transform_iterator() {
        let data = [1, 2, 3, 4];
        let mut it = transform(&data, |x| x * 10);
        assert_eq!(it.get(), 10);
        it.increment();
        assert_eq!(it.get(), 20);
        let it2 = it.add(2);
        assert_eq!(it2.get(), 40);
        assert_eq!(it2.distance(&it), 2);
        assert!(it.not_equal(&it2));
        it.increment().decrement();
        assert_eq!(it.get(), 20);
    }

    #[test]
    fn test_sorted_join_merge_keys_first_second() {
        let a = [("a", 1), ("b", 2), ("b", 3), ("d", 4)];
        let b = [("b", 20), ("d", 40)];
        let joined = join_sorted_lists(
            &a,
            &b,
            |x| x.0.to_string(),
            |x| x.0.to_string(),
            |x, y| (x.1, y.1),
        )
        .unwrap();
        assert_eq!(joined, vec![(2, 20), (3, 20), (4, 40)]);

        let duplicate = [("b", 20), ("b", 21)];
        assert!(join_sorted_lists(
            &[("b", 2), ("c", 3)],
            &duplicate,
            |x| x.0.to_string(),
            |x| x.0.to_string(),
            |x, y| (x.1, y.1),
        )
        .is_err());

        let merged = merge_keys(&[(0, 7), (0, 8), (2, 9)], |x| x.0, |x| x.1, 0, 3);
        assert_eq!(merged[0].1.iter().copied().collect::<Vec<_>>(), vec![7, 8]);
        assert!(merged[1].1.is_empty());
        assert_eq!(merged[2].1.iter().copied().collect::<Vec<_>>(), vec![9]);
        assert_eq!(first(&(1, 2)), 1);
        assert_eq!(second(&(1, 2)), 2);
    }

    #[test]
    fn test_hyper_log_log_add_estimate_and_merge() {
        assert!(HyperLogLog::new(3).is_err());
        assert!(HyperLogLog::new(21).is_err());

        let mut h = HyperLogLog::new(10).unwrap();
        assert_eq!(h.estimate(), 0);
        for i in 0..1_000 {
            h.add(i);
        }
        let estimate = h.estimate();
        assert!(estimate > 900 && estimate < 1_100, "{estimate}");

        let mut other = HyperLogLog::new(10).unwrap();
        for i in 1_000..2_000 {
            other.add(i);
        }
        h.merge(&other).unwrap();
        let estimate = h.estimate();
        assert!(estimate > 1_800 && estimate < 2_200, "{estimate}");

        assert!(h.merge(&HyperLogLog::new(9).unwrap()).is_err());
    }

    #[test]
    fn test_batch_binary_search() {
        let q = [1, 4, 7, 10];
        let t = [2, 5, 8];
        let mut out = Vec::new();
        batch_binary_search(&q, &t, &mut out, |a, b| a < b, 0);
        assert_eq!(out, vec![0, 0, 1, 2]);
    }

    #[test]
    fn test_merge_capped() {
        let mut out = Vec::new();
        let count = merge_capped(&[1, 4, 7], &[2, 3, 8], 5, &mut out);
        assert_eq!(out, vec![1, 2, 3, 4, 7]);
        assert_eq!(count, 2);
    }

    #[test]
    fn test_varuint32_roundtrip_boundaries_and_encoding() {
        let values = [
            0,
            1,
            (1 << 7) - 1,
            1 << 7,
            (1 << 14) - 1,
            1 << 14,
            (1 << 21) - 1,
            1 << 21,
            (1 << 28) - 1,
            1 << 28,
            u32::MAX,
        ];
        for value in values {
            let mut out = Vec::new();
            write_varuint32(value, &mut out);
            let (read, consumed) = read_varuint32(&out).unwrap();
            assert_eq!(read, value);
            assert_eq!(consumed, out.len());
        }

        let mut out = Vec::new();
        write_varuint32(0x7f, &mut out);
        assert_eq!(out, vec![0xff]);
        out.clear();
        write_varuint32(0x80, &mut out);
        assert_eq!(out, vec![0x02, 0x02]);
        out.clear();
        write_varuint32(0x123, &mut out);
        assert_eq!(out, vec![0x8e, 0x04]);
        assert!(read_varuint32(&[0]).is_err());
    }

    #[test]
    fn test_partition_table() {
        let rows = [1, 1, 2, 2, 2, 3, 4, 4, 5];
        let parts = partition_table(&rows, 3, |x| *x);
        assert_eq!(parts, vec![0, 5, 6, 9]);
    }

    #[test]
    fn test_sort_by_value() {
        let sorted = sort_by_value(&[30, 10, 20, 10], 4);
        assert_eq!(sorted, vec![(10, 1), (10, 3), (20, 2), (30, 0)]);
    }

    #[test]
    fn test_dsu_find_and_unite() {
        let mut dsu = Dsu::new(5);
        assert_eq!(dsu.find(3), 3);
        dsu.unite(0, 1);
        dsu.unite(1, 2);
        let root = dsu.find(0);
        assert_eq!(dsu.find(1), root);
        assert_eq!(dsu.find(2), root);
        assert_ne!(dsu.find(3), root);
        dsu.unite(3, 4);
        let other = dsu.find(3);
        assert_eq!(dsu.find(4), other);
    }

    #[test]
    fn test_insertion_sort_and_merge_sort() {
        let mut values = [5, 1, 4, 2, 3];
        insertion_sort(&mut values);
        assert_eq!(values, [1, 2, 3, 4, 5]);

        let mut pairs = [(2, 'b'), (1, 'z'), (2, 'a'), (1, 'a')];
        insertion_sort_by(&mut pairs, |a, b| a.1 < b.1);
        assert_eq!(pairs, [(2, 'a'), (1, 'a'), (2, 'b'), (1, 'z')]);

        let mut values = [9, 3, 7, 1, 5, 2, 8, 6, 4];
        merge_sort(&mut values, 4);
        assert_eq!(values, [1, 2, 3, 4, 5, 6, 7, 8, 9]);

        let mut values = [1, 2, 3, 4, 5];
        merge_sort_by(&mut values, 2, |a, b| a > b, 0);
        assert_eq!(values, [5, 4, 3, 2, 1]);
    }

    #[test]
    fn test_radix_cluster_and_parallel_radix_cluster() {
        let input = [(5u32, 'a'), (2, 'b'), (7, 'c'), (4, 'd'), (1, 'e')];
        let relation = Relation::new(&input, input.len());
        assert_eq!(relation[0], (5, 'a'));
        assert_eq!(relation.part(1, 2).data, &[(2, 'b'), (7, 'c')]);
        assert!(relation.end().is_empty());

        let extract = ExtractBits::new(4u32, 1);
        assert_eq!(extract.get(0b1011), 0b0001);

        let mut out = [(0u32, '\0'); 5];
        let mut hst = [0usize; 4];
        radix_cluster(relation, 0, &mut out, &mut hst, |x| x.0, 2, false);
        assert_eq!(out, [(4, 'd'), (5, 'a'), (1, 'e'), (2, 'b'), (7, 'c')]);
        assert_eq!(hst, [1, 3, 4, 5]);

        let mut buffered = [(0u32, '\0'); 5];
        let mut hst = [0usize; 4];
        radix_cluster(relation, 0, &mut buffered, &mut hst, |x| x.0, 2, true);
        assert_eq!(buffered, out);

        let mut parallel = [(0u32, '\0'); 5];
        parallel_radix_cluster(relation, 0, &mut parallel, 2, |x| x.0, 2);
        assert_eq!(parallel, out);
    }

    #[test]
    fn test_hash_join_groups_values_by_shared_key() {
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        struct Item {
            key: u32,
            value: u32,
        }

        impl HashJoinRecord for Item {
            type Key = u32;
            type Value = u32;

            fn key(&self) -> Self::Key {
                self.key
            }

            fn value(&self) -> Self::Value {
                self.value
            }
        }

        let r = [
            Item { key: 1, value: 10 },
            Item { key: 2, value: 20 },
            Item { key: 1, value: 11 },
            Item { key: 4, value: 40 },
        ];
        let s = [
            Item { key: 1, value: 100 },
            Item { key: 3, value: 300 },
            Item { key: 1, value: 101 },
            Item { key: 2, value: 200 },
        ];
        let (out_r, out_s) = hash_join(Relation::new(&r, r.len()), Relation::new(&s, s.len()), 4);
        let mut it_r = out_r.begin_with_elem_size(4);
        let mut it_s = out_s.begin_with_elem_size(4);

        assert!(it_r.good());
        assert!(it_s.good());
        assert_eq!(it_r.count(), 2);
        assert_eq!(it_s.count(), 2);
        assert_eq!(u32_values(it_r.range_bytes(4)), vec![10, 11]);
        assert_eq!(u32_values(it_s.range_bytes(4)), vec![100, 101]);
        it_r.increment(4);
        it_s.increment(4);
        assert!(it_r.good());
        assert!(it_s.good());
        assert_eq!(u32_values(it_r.range_bytes(4)), vec![20]);
        assert_eq!(u32_values(it_s.range_bytes(4)), vec![200]);
        it_r.increment(4);
        it_s.increment(4);
        assert!(!it_r.good());
        assert!(!it_s.good());
    }

    fn u32_values(bytes: &[u8]) -> Vec<u32> {
        bytes
            .chunks_exact(4)
            .map(|chunk| u32::from_ne_bytes(chunk.try_into().unwrap()))
            .collect()
    }

    #[test]
    fn test_radix_sort() {
        let mut input = [(5u32, 'a'), (2, 'b'), (7, 'c'), (4, 'd'), (1, 'e')];
        radix_sort(&mut input, 7, 1, |x| x.0, 2);
        assert_eq!(input, [(1, 'e'), (2, 'b'), (4, 'd'), (5, 'a'), (7, 'c')]);

        let mut input = [(3u32, 0), (1, 1), (3, 2), (0, 3), (1, 4)];
        radix_sort(&mut input, 3, 3, |x| x.0, 2);
        assert_eq!(input, [(0, 3), (1, 1), (1, 4), (3, 0), (3, 2)]);
    }

    #[test]
    fn test_merge_sorted_files() {
        let files = vec![vec![1, 4, 9], vec![2, 3, 10], vec![], vec![0, 5]];
        let mut out = Vec::new();
        merge_sorted_files::<i32, _, _>(files, |x| out.push(x));
        assert_eq!(out, vec![0, 1, 2, 3, 4, 5, 9, 10]);
    }

    #[test]
    fn test_external_sorter() {
        let mut s = ExternalSorter::default();
        for value in [5, 1, 4, 1, 3, 2] {
            s.push(value);
        }
        assert_eq!(s.count(), 6);
        s.init_read();
        let mut out = Vec::new();
        while s.good() {
            out.push(*s.current().unwrap());
            s.increment();
        }
        assert_eq!(out, vec![1, 1, 2, 3, 4, 5]);

        let mut s = ExternalSorter::new(|a: &i32, b: &i32| a > b);
        for value in [1, 3, 2] {
            s.push(value);
        }
        s.init_read();
        let mut out = Vec::new();
        while let Some(&value) = s.current() {
            out.push(value);
            s.increment();
        }
        assert_eq!(out, vec![3, 2, 1]);

        let p = ("abc".to_string(), 7u32);
        assert_eq!(
            alloc_size_string_u32_pair(&p),
            4 + std::mem::size_of::<String>() + 3
        );
        assert_eq!(alloc_size(&7u32), std::mem::size_of::<u32>());
    }

    #[test]
    fn test_dsu() {
        let mut d = Dsu::new(5);
        d.unite(0, 1);
        d.unite(3, 4);
        assert_eq!(d.find(0), d.find(1));
        assert_ne!(d.find(0), d.find(2));
        assert_eq!(d.find(3), d.find(4));
        d.unite(1, 4);
        assert_eq!(d.find(0), d.find(3));
    }

    #[test]
    fn test_neighbor_count_cpp_leaves() {
        let centroids = [u32::MAX, u32::MAX, 7, u32::MAX];
        let edges = [
            Edge::new(0, 1, 0.0),
            Edge::new(0, 1, 0.0),
            Edge::new(0, 2, 0.0),
            Edge::new(0, 3, 0.0),
        ];
        assert_eq!(neighbor_count2(&edges, &centroids), 2);

        let mut mutable_edges = [
            Edge::new(0, 1, 0.0),
            Edge::new(0, 2, 0.0),
            Edge::new(0, 3, 0.0),
            Edge::new(u32::MAX, 0, 0.0),
        ];
        assert_eq!(neighbor_count(&mut mutable_edges, &centroids), 2);
        assert_eq!(mutable_edges[0].node2, 1);
        assert_eq!(mutable_edges[1].node2, 3);
        assert_eq!(mutable_edges[2].node1, u32::MAX);

        let members = [5u32, 2, 11, 4];
        assert_eq!(
            neighbor_count_with_member_counts(0, &edges, &centroids, &members),
            13
        );
    }

    #[test]
    fn test_greedy_vertex_cover_basic_and_reassign() {
        let mut neighbors = FlatArray::<Edge<u32>, u32>::new();
        neighbors.push_back(&[Edge::new(0, 1, 0.5), Edge::new(0, 2, 0.1)]);
        neighbors.push_back(&[Edge::new(1, 0, 0.5), Edge::new(1, 2, 0.9)]);
        neighbors.push_back(&[Edge::new(2, 0, 0.1), Edge::new(2, 1, 0.9)]);
        let centroids = greedy_vertex_cover(&mut neighbors, None, false, true, 0);
        assert_eq!(centroids.len(), 3);
        assert!(centroids.iter().all(|&x| x != u32::MAX));
        assert!(centroids.iter().any(|&x| centroids[x as usize] == x));
    }

    #[test]
    fn test_greedy_vertex_cover_weighted_cc_and_merge() {
        let mut neighbors = FlatArray::<Edge<u64>, u64>::new();
        neighbors.push_back(&[Edge::new(0, 1, 1.0)]);
        neighbors.push_back(&[Edge::new(1, 2, 1.0)]);
        neighbors.push_back(&[Edge::new(2, 3, 1.0)]);
        neighbors.push_back(&[]);
        let members = [1u64, 3, 2, 1];
        let centroids = greedy_vertex_cover(&mut neighbors, Some(&members), true, false, 2);
        assert!(centroids.iter().all(|&x| x != u64::MAX));
        assert!(centroids.iter().any(|&x| centroids[x as usize] == x));
        assert_eq!(centroids[2], centroids[1]);
    }
}
