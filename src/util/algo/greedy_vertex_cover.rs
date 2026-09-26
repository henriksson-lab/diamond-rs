use std::collections::{BinaryHeap, VecDeque};

use crate::util::data_structures::FlatArray;
use crate::util::log_stream::{message_stream, TaskTimer};

use super::{AlgoInt, Edge};

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

fn int_to_usize<Int: AlgoInt>(value: Int) -> usize
where
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    value.try_into().unwrap()
}

fn usize_to_int<Int: AlgoInt>(value: usize) -> Int
where
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
{
    Int::try_from(value).unwrap()
}

/// Translation of the first `neighbor_count` overload in
/// `util/algo/greedy_vertex_cover.cpp`.
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

/// Translation of the compacting `neighbor_count` overload.
///
/// Live, distinct neighbors are moved to the front and the first dead slot is
/// marked with the same `NIL` sentinel used by the C++ implementation.
pub fn neighbor_count<Int>(edges: &mut [Edge<Int>], centroids: &[Int]) -> Int
where
    Int: AlgoInt,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut n = 0usize;
    let mut last = Int::NIL;
    let mut write = 0usize;
    let mut read = 0usize;
    while read < edges.len() {
        if edges[read].node1 == Int::NIL {
            break;
        }
        if centroids[int_to_usize(edges[read].node2)] == Int::NIL && edges[read].node2 != last {
            n += 1;
            last = edges[read].node2;
            if write < read {
                edges.swap(write, read);
            }
            write += 1;
        }
        read += 1;
    }
    if write < edges.len() {
        edges[write].node1 = Int::NIL;
    }
    usize_to_int(n)
}

/// Translation of the member-count-weighted `neighbor_count` overload.
pub fn neighbor_count_with_member_counts<Int>(
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
        let centroid = int_to_usize(centroids[i]);
        if centroids[centroid] != centroids[i] {
            centroids[i] = centroids[centroid];
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
    let mut queue = VecDeque::new();
    for edge in neighbors.range(rep) {
        if centroids[int_to_usize(edge.node2)] == Int::NIL {
            queue.push_back(Entry::new(edge.node2, Int::from(1)));
        }
    }
    while let Some(entry) = queue.pop_front() {
        let node = int_to_usize(entry.node);
        if centroids[node] != Int::NIL || entry.depth > depth {
            continue;
        }
        for edge in neighbors.range(entry.node) {
            if centroids[int_to_usize(edge.node2)] == Int::NIL {
                queue.push_back(Entry::new(edge.node2, entry.depth + Int::from(1)));
            }
        }
        centroids[node] = rep;
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
    let mut timer = TaskTimer::with_msg("Computing edge counts", 1);
    let size = int_to_usize(neighbors.size());
    let mut queue = BinaryHeap::<(Int, Int)>::new();
    let mut centroids = vec![Int::NIL; size];
    for i in 0..size {
        let node = usize_to_int::<Int>(i);
        let count = match member_counts {
            Some(member_counts) => neighbor_count_with_member_counts(
                node,
                neighbors.range(node),
                &centroids,
                member_counts,
            ),
            None => neighbors.count(node),
        };
        queue.push((count, node));
    }

    timer.go_string("Computing vertex cover");
    let mut cluster_count = 0i64;
    while let Some((_, node)) = queue.pop() {
        if centroids[int_to_usize(node)] != Int::NIL {
            continue;
        }
        let count = match member_counts {
            Some(member_counts) => neighbor_count_with_member_counts(
                node,
                neighbors.range(node),
                &centroids,
                member_counts,
            ),
            None => neighbor_count(neighbors.range_mut(node), &centroids),
        };
        if queue.peek().is_some_and(|top| count < top.0) {
            queue.push((count, node));
        } else {
            if connected_component_depth > Int::from(0) {
                make_cluster_cc(node, neighbors, &mut centroids, connected_component_depth);
            } else {
                make_cluster_gvc(node, neighbors, &mut centroids, merge_recursive);
            }
            cluster_count += 1;
        }
    }
    timer.finish();
    if let Ok(mut stream) = message_stream().lock() {
        let _ = stream
            .write("Cluster count = ")
            .and_then(|stream| stream.write(cluster_count))
            .and_then(|stream| stream.endl());
    }

    if reassign {
        timer.go_string("Computing reassignment");
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
        timer.go_string("Computing merges");
        fix_assignment(&mut centroids);
    }

    centroids
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph_u32(rows: &[&[Edge<u32>]]) -> FlatArray<Edge<u32>, u32> {
        let mut graph = FlatArray::new();
        for row in rows {
            graph.push_back(row);
        }
        graph
    }

    #[test]
    fn neighbor_count2_ignores_assigned_nodes_and_adjacent_duplicates() {
        let centroids = [u32::MAX, u32::MAX, 2, u32::MAX];
        let edges = [
            Edge::new(0, 1, 0.0),
            Edge::new(0, 1, 0.0),
            Edge::new(0, 2, 0.0),
            Edge::new(0, 3, 0.0),
        ];
        assert_eq!(neighbor_count2(&edges, &centroids), 2);
    }

    #[test]
    fn neighbor_count_compacts_live_distinct_neighbors_and_sets_sentinel() {
        let centroids = [u32::MAX, u32::MAX, 2, u32::MAX];
        let mut edges = [
            Edge::new(0, 2, 0.0),
            Edge::new(0, 1, 0.0),
            Edge::new(0, 1, 0.0),
            Edge::new(0, 3, 0.0),
        ];
        assert_eq!(neighbor_count(&mut edges, &centroids), 2);
        assert_eq!((edges[0].node2, edges[1].node2), (1, 3));
        assert_eq!(edges[2].node1, u32::MAX);
    }

    #[test]
    fn weighted_neighbor_count_includes_node_and_each_live_edge() {
        let centroids = [u64::MAX, u64::MAX, 2, u64::MAX];
        let members = [5, 2, 11, 4];
        let edges = [
            Edge::new(0, 1, 0.0),
            Edge::new(0, 1, 0.0),
            Edge::new(0, 2, 0.0),
            Edge::new(0, 3, 0.0),
        ];
        assert_eq!(
            neighbor_count_with_member_counts(0, &edges, &centroids, &members),
            13
        );
    }

    #[test]
    fn make_cluster_gvc_can_merge_an_existing_representative() {
        let graph = graph_u32(&[&[], &[Edge::new(1, 2, 1.0)], &[], &[]]);
        let mut centroids = [0, u32::MAX, 2, 2];
        make_cluster_gvc(1, &graph, &mut centroids, false);
        assert_eq!(centroids, [0, 1, 2, 2]);

        let mut centroids = [0, u32::MAX, 2, 2];
        make_cluster_gvc(1, &graph, &mut centroids, true);
        assert_eq!(centroids, [0, 1, 1, 2]);
    }

    #[test]
    fn make_cluster_cc_respects_breadth_first_depth() {
        let graph = graph_u32(&[
            &[Edge::new(0, 1, 1.0)],
            &[Edge::new(1, 2, 1.0)],
            &[Edge::new(2, 3, 1.0)],
            &[],
        ]);
        let mut centroids = [u32::MAX; 4];
        make_cluster_cc(0, &graph, &mut centroids, 2);
        assert_eq!(centroids, [0, 0, 0, u32::MAX]);
    }

    #[test]
    fn greedy_vertex_cover_clusters_a_star_exactly() {
        let mut graph = graph_u32(&[
            &[Edge::new(0, 1, 1.0), Edge::new(0, 2, 1.0)],
            &[Edge::new(1, 0, 1.0)],
            &[Edge::new(2, 0, 1.0)],
        ]);
        assert_eq!(
            greedy_vertex_cover(&mut graph, None, false, false, 0),
            [0, 0, 0]
        );
    }

    #[test]
    fn greedy_vertex_cover_reassigns_to_the_highest_weight_representative() {
        let mut graph = graph_u32(&[
            &[],
            &[],
            &[Edge::new(2, 0, 9.0)],
            &[Edge::new(3, 0, 1.0), Edge::new(3, 1, 1.0)],
        ]);
        assert_eq!(
            greedy_vertex_cover(&mut graph, None, false, true, 0),
            [2, 3, 2, 3]
        );
    }

    #[test]
    fn greedy_vertex_cover_flattens_recursive_merges() {
        let rows: &[&[Edge<u32>]] = &[&[], &[Edge::new(1, 2, 1.0)], &[Edge::new(2, 3, 1.0)], &[]];
        let mut graph = graph_u32(rows);
        assert_eq!(
            greedy_vertex_cover(&mut graph, None, false, false, 0),
            [0, 1, 2, 2]
        );

        let mut graph = graph_u32(rows);
        assert_eq!(
            greedy_vertex_cover(&mut graph, None, true, false, 0),
            [0, 1, 1, 1]
        );
    }

    #[test]
    fn greedy_vertex_cover_supports_u64_member_counts_and_cc_mode() {
        let mut graph = FlatArray::<Edge<u64>, u64>::new();
        graph.push_back(&[Edge::new(0, 1, 1.0)]);
        graph.push_back(&[Edge::new(1, 2, 1.0)]);
        graph.push_back(&[Edge::new(2, 3, 1.0)]);
        graph.push_back(&[]);
        let members = [1u64, 3, 2, 1];
        assert_eq!(
            greedy_vertex_cover(&mut graph, Some(&members), false, false, 2),
            [0, 1, 1, 1]
        );
    }
}
