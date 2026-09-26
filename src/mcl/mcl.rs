//! Markov clustering translated from `diamond/src/contrib/mcl/mcl.cpp` and
//! its inline `dense.h` and `sparse.h` implementations.

use nalgebra::{linalg::Schur, Complex, ComplexField, DMatrix, DVector, Normed};
use std::collections::{BTreeMap, BTreeSet};

pub const MASK_INVERSE: u64 = 0xc000_0000_0000_0000;
pub const MASK_NORMAL_NODE: u64 = 0x4000_0000_0000_0000;
pub const MASK_ATTRACTOR_NODE: u64 = 0x8000_0000_0000_0000;
pub const MASK_SINGLE_NODE: u64 = 0xc000_0000_0000_0000;
pub const DEFAULT_CLUSTERING_THRESHOLD: f64 = 50.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triplet {
    pub row: usize,
    pub col: usize,
    pub value: f32,
}

impl Triplet {
    pub const fn new(row: usize, col: usize, value: f32) -> Self {
        Self { row, col, value }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DenseMatrix {
    n: usize,
    values: Vec<f32>,
}

impl DenseMatrix {
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            values: vec![0.0; n * n],
        }
    }

    pub fn from_triplets(n: usize, triplets: &[Triplet], symmetric: bool) -> Self {
        let mut matrix = Self::zeros(n);
        for triplet in triplets {
            matrix[(triplet.row, triplet.col)] = triplet.value;
            if symmetric && triplet.row != triplet.col {
                matrix[(triplet.col, triplet.row)] = triplet.value;
            }
        }
        matrix
    }

    pub fn size(&self) -> usize {
        self.n
    }

    pub fn norm(&self) -> f32 {
        self.values
            .iter()
            .map(|value| value * value)
            .sum::<f32>()
            .sqrt()
    }

    fn multiply(&self, rhs: &Self) -> Self {
        assert_eq!(self.n, rhs.n);
        let mut result = Self::zeros(self.n);
        for col in 0..self.n {
            for k in 0..self.n {
                let y = rhs[(k, col)];
                for row in 0..self.n {
                    result[(row, col)] += self[(row, k)] * y;
                }
            }
        }
        result
    }
}

impl std::ops::Index<(usize, usize)> for DenseMatrix {
    type Output = f32;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.values[col * self.n + row]
    }
}

impl std::ops::IndexMut<(usize, usize)> for DenseMatrix {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.values[col * self.n + row]
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SparseMatrix {
    n: usize,
    values: BTreeMap<(usize, usize), f32>,
}

impl SparseMatrix {
    pub fn from_triplets(n: usize, triplets: &[Triplet], symmetric: bool) -> Self {
        let mut values = BTreeMap::new();
        for triplet in triplets {
            if triplet.value.abs() > f32::EPSILON {
                values.insert((triplet.row, triplet.col), triplet.value);
                if symmetric && triplet.row != triplet.col {
                    values.insert((triplet.col, triplet.row), triplet.value);
                }
            }
        }
        Self { n, values }
    }

    pub fn size(&self) -> usize {
        self.n
    }

    pub fn triplets(&self) -> impl Iterator<Item = Triplet> + '_ {
        self.values
            .iter()
            .map(|(&(row, col), &value)| Triplet::new(row, col, value))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MclConfig {
    pub inflation: f32,
    pub expansion: f32,
    pub max_iter: u32,
    pub sparsity_switch: f32,
    pub symmetric: bool,
}

impl Default for MclConfig {
    fn default() -> Self {
        Self {
            inflation: 2.0,
            expansion: 2.0,
            max_iter: 100,
            sparsity_switch: 0.8,
            symmetric: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeLabel {
    Singleton,
    Attractor,
    Normal,
    Unassigned,
}

impl NodeLabel {
    pub const fn code(self) -> char {
        match self {
            Self::Singleton => 's',
            Self::Attractor => 'a',
            Self::Normal => 'n',
            Self::Unassigned => 'u',
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Assignment {
    pub sequence_id: String,
    pub cluster_id: Option<u64>,
    pub label: NodeLabel,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MclStatistics {
    pub elements: usize,
    pub components: usize,
    pub non_singleton_components: usize,
    pub clusters: usize,
    pub singleton_clusters: usize,
    pub dense_calculations: usize,
    pub sparse_calculations: usize,
    pub failed_to_converge: usize,
    pub component_sizes: Vec<usize>,
    pub sparsities: Vec<f32>,
    pub average_neighbors: Vec<f32>,
    pub rough_memory_bytes: f32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MclResult {
    pub assignments: Vec<Assignment>,
    pub encoded_assignments: Vec<u64>,
    pub statistics: MclStatistics,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GraphHandle {
    pub sequence_count: usize,
    pub symmetric: bool,
    pub edges: Vec<Triplet>,
    pub similarity: String,
    pub threshold: Option<f64>,
}

/// Explicit-input counterpart of the source graph/search dispatcher.
pub fn get_graph_handle(
    sequence_count: usize,
    edges: &[Triplet],
    symmetric: bool,
    similarity: Option<&str>,
    threshold: Option<f64>,
) -> GraphHandle {
    let similarity = similarity.unwrap_or("normalized_bitscore_global");
    let threshold = if similarity == "normalized_bitscore_global" {
        threshold
            .or(Some(DEFAULT_CLUSTERING_THRESHOLD))
            .filter(|value| *value != 0.0)
    } else {
        threshold
    };
    GraphHandle {
        sequence_count,
        symmetric,
        edges: edges.to_vec(),
        similarity: similarity.to_owned(),
        threshold,
    }
}

#[derive(Debug, Default)]
pub struct Mcl {
    pub failed_to_converge: usize,
}

impl Mcl {
    pub const fn get_description() -> &'static str {
        "Markov clustering according to doi:10.1137/040608635"
    }

    pub const fn get_key() -> &'static str {
        "mcl"
    }

    pub fn get_dense_matrix_and_clear(
        order: &[usize],
        triplets: &mut Vec<Triplet>,
        symmetric: bool,
    ) -> DenseMatrix {
        let matrix = DenseMatrix::from_triplets(order.len(), triplets, symmetric);
        triplets.clear();
        matrix
    }

    pub fn get_sparse_matrix_and_clear(
        order: &[usize],
        triplets: &mut Vec<Triplet>,
        symmetric: bool,
    ) -> SparseMatrix {
        let matrix = SparseMatrix::from_triplets(order.len(), triplets, symmetric);
        triplets.clear();
        matrix
    }

    pub fn get_gamma_dense(input: &DenseMatrix, r: f32) -> DenseMatrix {
        let mut output = DenseMatrix::zeros(input.n);
        for col in 0..input.n {
            let column_sum = (0..input.n)
                .map(|row| input[(row, col)].powf(r))
                .sum::<f32>();
            for row in 0..input.n {
                output[(row, col)] = input[(row, col)].powf(r) / column_sum;
            }
        }
        output
    }

    pub fn sparse_matrix_get_gamma(
        input: &SparseMatrix,
        r: f32,
        thread: usize,
        threads: usize,
    ) -> Vec<Triplet> {
        assert!(threads > 0 && thread < threads);
        let mut result = Vec::new();
        for col in (thread..input.n).step_by(threads) {
            let column: Vec<_> = input
                .values
                .iter()
                .filter(|((_, entry_col), _)| *entry_col == col)
                .map(|(&(row, _), &value)| (row, value))
                .collect();
            let sum = column.iter().map(|(_, value)| value.powf(r)).sum::<f32>();
            for (row, value) in column {
                let value = value.powf(r) / sum;
                if value.abs() > f32::EPSILON {
                    result.push(Triplet::new(row, col, value));
                }
            }
        }
        result
    }

    pub fn get_gamma_sparse(input: &SparseMatrix, r: f32, threads: usize) -> SparseMatrix {
        let mut values = BTreeMap::new();
        for thread in 0..threads.max(1) {
            for triplet in Self::sparse_matrix_get_gamma(input, r, thread, threads.max(1)) {
                values.insert((triplet.row, triplet.col), triplet.value);
            }
        }
        SparseMatrix { n: input.n, values }
    }

    pub fn sparse_matrix_multiply(
        lhs: &SparseMatrix,
        rhs: &SparseMatrix,
        thread: usize,
        threads: usize,
    ) -> Vec<Triplet> {
        assert_eq!(lhs.n, rhs.n);
        assert!(threads > 0 && thread < threads);
        let mut result = Vec::new();
        for col in (thread..rhs.n).step_by(threads) {
            let rhs_column: Vec<_> = rhs
                .values
                .iter()
                .filter(|((_, entry_col), _)| *entry_col == col)
                .map(|(&(row, _), &value)| (row, value))
                .collect();
            let mut column = vec![0.0f32; lhs.n];
            for (k, y) in rhs_column {
                for (&(row, _), &x) in lhs
                    .values
                    .iter()
                    .filter(|((_, entry_col), _)| *entry_col == k)
                {
                    column[row] += x * y;
                }
            }
            for (row, value) in column.into_iter().enumerate() {
                if value.abs() > f32::EPSILON {
                    result.push(Triplet::new(row, col, value));
                }
            }
        }
        result
    }

    pub fn get_exp_sparse(
        input: &SparseMatrix,
        r: f32,
        threads: usize,
    ) -> Result<SparseMatrix, String> {
        if r.fract() != 0.0 {
            return Err("Eigen does not provide an eigenvalue solver for sparse matrices".into());
        }
        let mut output = input.clone();
        for _ in 1..r as u32 {
            let mut values = BTreeMap::new();
            for thread in 0..threads.max(1) {
                for triplet in Self::sparse_matrix_multiply(input, &output, thread, threads.max(1))
                {
                    values.insert((triplet.row, triplet.col), triplet.value);
                }
            }
            output.values = values;
        }
        output.values.retain(|_, value| value.abs() > f32::EPSILON);
        Ok(output)
    }

    pub fn sparse_matrix_get_norm(input: &SparseMatrix, _threads: usize) -> f32 {
        input
            .values
            .values()
            .map(|value| value.powi(2))
            .sum::<f32>()
            .sqrt()
    }

    pub fn get_exp_dense(input: &DenseMatrix, r: f32) -> DenseMatrix {
        if r.fract() == 0.0 {
            let mut output = input.clone();
            for _ in 1..r as u32 {
                output = output.multiply(input);
            }
            return output;
        }

        let matrix = DMatrix::from_fn(input.n, input.n, |row, col| {
            Complex::new(input[(row, col)], 0.0)
        });
        let (q, t) = Schur::new(matrix).unpack();
        let mut triangular_vectors = DMatrix::<Complex<f32>>::zeros(input.n, input.n);
        for col in 0..input.n {
            let lambda = t[(col, col)];
            triangular_vectors[(col, col)] = Complex::new(1.0, 0.0);
            for row in (0..col).rev() {
                let sum = ((row + 1)..=col)
                    .map(|k| t[(row, k)] * triangular_vectors[(k, col)])
                    .sum::<Complex<f32>>();
                let denominator = t[(row, row)] - lambda;
                if denominator.norm() > f32::EPSILON {
                    triangular_vectors[(row, col)] = -sum / denominator;
                }
            }
        }
        let vectors = q * triangular_vectors;
        let Some(inverse) = vectors.clone().try_inverse() else {
            // This is the source's singular-eigenvector path: `out` remains
            // unchanged. The caller allocates it as zero.
            return DenseMatrix::zeros(input.n);
        };
        if vectors.determinant().norm() * 0.5 <= f32::EPSILON {
            return DenseMatrix::zeros(input.n);
        }
        let diagonal =
            DMatrix::from_diagonal(&DVector::from_fn(input.n, |row, _| t[(row, row)].powf(r)));
        let powered = vectors * diagonal * inverse;
        let mut output = DenseMatrix::zeros(input.n);
        for col in 0..input.n {
            for row in 0..input.n {
                output[(row, col)] = powered[(row, col)].re;
            }
        }
        output
    }

    pub fn markov_process_dense(
        &mut self,
        matrix: &mut DenseMatrix,
        inflation: f32,
        expansion: f32,
        max_iter: u32,
    ) {
        *matrix = Self::get_gamma_dense(matrix, 1.0);
        let mut iteration = 0;
        let mut diff_norm = f32::MAX;
        while iteration < max_iter && diff_norm > f32::EPSILON {
            let expanded = Self::get_exp_dense(matrix, expansion);
            let update = Self::get_gamma_dense(&expanded, inflation);
            diff_norm = matrix
                .values
                .iter()
                .zip(&update.values)
                .map(|(left, right)| (left - right).powi(2))
                .sum::<f32>()
                .sqrt();
            *matrix = update;
            iteration += 1;
        }
        if iteration == max_iter {
            self.failed_to_converge += 1;
        }
    }

    pub fn markov_process_sparse(
        &mut self,
        matrix: &mut SparseMatrix,
        inflation: f32,
        expansion: f32,
        max_iter: u32,
        threads: usize,
    ) -> Result<(), String> {
        *matrix = Self::get_gamma_sparse(matrix, 1.0, threads);
        let mut iteration = 0;
        let mut diff_norm = f32::MAX;
        while iteration < max_iter && diff_norm > f32::EPSILON {
            let expanded = Self::get_exp_sparse(matrix, expansion, threads)?;
            let update = Self::get_gamma_sparse(&expanded, inflation, threads);
            let mut difference = matrix.clone();
            for (&key, &value) in &update.values {
                *difference.values.entry(key).or_insert(0.0) -= value;
            }
            difference
                .values
                .retain(|_, value| value.abs() > f32::EPSILON);
            diff_norm = Self::sparse_matrix_get_norm(&difference, threads);
            *matrix = update;
            iteration += 1;
        }
        if iteration == max_iter {
            self.failed_to_converge += 1;
        }
        Ok(())
    }

    pub fn print_stats(
        sequence_count: usize,
        edges: &[Triplet],
        config: MclConfig,
        threads: usize,
    ) -> MclStatistics {
        let components = graph_components(sequence_count, edges);
        let mut component_sizes = Vec::new();
        let mut sparsities = Vec::new();
        let mut average_neighbors = Vec::new();
        let mut memories = Vec::new();
        for component in components.iter().filter(|component| component.len() > 1) {
            let members: BTreeSet<_> = component.iter().copied().collect();
            let local_edges = edges
                .iter()
                .filter(|edge| members.contains(&edge.row) && members.contains(&edge.col))
                .count();
            let size = component.len();
            let non_diagonal = edges
                .iter()
                .filter(|edge| {
                    edge.row != edge.col
                        && members.contains(&edge.row)
                        && members.contains(&edge.col)
                })
                .count();
            let sparsity = 1.0 - local_edges as f32 / (size * size) as f32;
            let neighbors = non_diagonal as f32 / size as f32;
            component_sizes.push(size);
            sparsities.push(sparsity);
            average_neighbors.push(neighbors);
            memories.push(if sparsity >= config.sparsity_switch {
                size as f32 * (1.0 + neighbors.powf(config.expansion)) * 12.0
            } else {
                4.0 * (size * size) as f32
            });
        }
        memories.sort_by(|left, right| right.total_cmp(left));
        let rough_memory_bytes = edges.len() as f32 * 12.0
            + memories
                .into_iter()
                .take(threads)
                .map(|value| 3.0 * value)
                .sum::<f32>();
        MclStatistics {
            elements: edges.len(),
            components: components.len(),
            non_singleton_components: components
                .iter()
                .filter(|component| component.len() > 1)
                .count(),
            component_sizes,
            sparsities,
            average_neighbors,
            rough_memory_bytes,
            ..MclStatistics::default()
        }
    }

    pub fn run(
        &mut self,
        sequence_ids: &[String],
        graph: &GraphHandle,
        config: MclConfig,
    ) -> Result<MclResult, String> {
        if graph.sequence_count != sequence_ids.len() {
            return Err("Graph and sequence counts differ".into());
        }
        for edge in &graph.edges {
            if edge.row >= graph.sequence_count || edge.col >= graph.sequence_count {
                return Err("Graph edge is outside the sequence range".into());
            }
        }

        let mut statistics = Self::print_stats(graph.sequence_count, &graph.edges, config, 1);
        let components = graph_components(graph.sequence_count, &graph.edges);
        let mut encoded = vec![0u64; graph.sequence_count];
        let mut next_cluster = 0u64;

        for order in components {
            if order.len() == 1 {
                encoded[order[0]] = MASK_SINGLE_NODE | next_cluster;
                next_cluster += 1;
                statistics.singleton_clusters += 1;
                continue;
            }

            let remap: BTreeMap<_, _> = order
                .iter()
                .enumerate()
                .map(|(local, &global)| (global, local))
                .collect();
            let mut triplets: Vec<_> = graph
                .edges
                .iter()
                .filter_map(|edge| {
                    Some(Triplet::new(
                        *remap.get(&edge.row)?,
                        *remap.get(&edge.col)?,
                        edge.value,
                    ))
                })
                .collect();
            let sparsity = 1.0 - triplets.len() as f32 / (order.len() * order.len()) as f32;
            let (sets, attractors) =
                if sparsity >= config.sparsity_switch && config.expansion.fract() == 0.0 {
                    statistics.sparse_calculations += 1;
                    let mut matrix =
                        Self::get_sparse_matrix_and_clear(&order, &mut triplets, config.symmetric);
                    self.markov_process_sparse(
                        &mut matrix,
                        config.inflation,
                        config.expansion,
                        config.max_iter,
                        1,
                    )?;
                    sets_from_entries(order.len(), matrix.triplets())
                } else {
                    statistics.dense_calculations += 1;
                    let mut matrix =
                        Self::get_dense_matrix_and_clear(&order, &mut triplets, config.symmetric);
                    self.markov_process_dense(
                        &mut matrix,
                        config.inflation,
                        config.expansion,
                        config.max_iter,
                    );
                    let entries = (0..matrix.n).flat_map(|col| {
                        let matrix = &matrix;
                        (0..matrix.n).filter_map(move |row| {
                            let value = matrix[(row, col)];
                            (value.abs() > f32::EPSILON).then_some(Triplet::new(row, col, value))
                        })
                    });
                    sets_from_entries(order.len(), entries)
                };

            for set in sets {
                for &local in &set {
                    let mask = if attractors.contains(&local) {
                        MASK_ATTRACTOR_NODE
                    } else {
                        MASK_NORMAL_NODE
                    };
                    encoded[order[local]] = mask | next_cluster;
                }
                if set.len() == 1 {
                    statistics.singleton_clusters += 1;
                }
                next_cluster += 1;
            }
        }

        statistics.clusters = next_cluster as usize;
        statistics.failed_to_converge = self.failed_to_converge;
        let assignments = sequence_ids
            .iter()
            .zip(&encoded)
            .map(|(id, &value)| {
                let label = match value & MASK_INVERSE {
                    MASK_SINGLE_NODE => NodeLabel::Singleton,
                    MASK_ATTRACTOR_NODE => NodeLabel::Attractor,
                    MASK_NORMAL_NODE => NodeLabel::Normal,
                    _ => NodeLabel::Unassigned,
                };
                Assignment {
                    sequence_id: crate::util::sequence::seqid(id),
                    cluster_id: (label != NodeLabel::Unassigned)
                        .then_some((value & !MASK_INVERSE) + 1),
                    label,
                }
            })
            .collect();
        Ok(MclResult {
            assignments,
            encoded_assignments: encoded,
            statistics,
        })
    }
}

fn graph_components(sequence_count: usize, edges: &[Triplet]) -> Vec<Vec<usize>> {
    let mut union = UnionFind::new(sequence_count);
    for edge in edges {
        if edge.row < sequence_count && edge.col < sequence_count {
            union.merge(edge.row, edge.col);
        }
    }
    let mut sets = union.sets();
    // `MCL::run` processes the largest disconnected components first.
    sets.sort_unstable_by(|left, right| right.len().cmp(&left.len()));
    sets
}

fn sets_from_entries(
    size: usize,
    entries: impl IntoIterator<Item = Triplet>,
) -> (Vec<Vec<usize>>, BTreeSet<usize>) {
    let mut union = UnionFind::new(size);
    let mut attractors = BTreeSet::new();
    for entry in entries {
        union.merge(entry.row, entry.col);
        if entry.row == entry.col {
            attractors.insert(entry.row);
        }
    }
    (union.sets(), attractors)
}

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u32>,
}

impl UnionFind {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn root(&mut self, value: usize) -> usize {
        if self.parent[value] != value {
            self.parent[value] = self.root(self.parent[value]);
        }
        self.parent[value]
    }

    fn merge(&mut self, left: usize, right: usize) {
        let mut left = self.root(left);
        let mut right = self.root(right);
        if left == right {
            return;
        }
        if self.rank[left] < self.rank[right] {
            std::mem::swap(&mut left, &mut right);
        }
        self.parent[right] = left;
        if self.rank[left] == self.rank[right] {
            self.rank[left] += 1;
        }
    }

    fn sets(&mut self) -> Vec<Vec<usize>> {
        let mut sets = BTreeMap::<usize, Vec<usize>>::new();
        for value in 0..self.parent.len() {
            sets.entry(self.root(value)).or_default().push(value);
        }
        sets.into_values().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(left: f32, right: f32) {
        assert!((left - right).abs() < 1e-4, "{left} != {right}");
    }

    #[test]
    fn dense_and_sparse_expansion_and_inflation_match() {
        let triplets = [
            Triplet::new(0, 0, 0.5),
            Triplet::new(1, 0, 0.5),
            Triplet::new(0, 1, 0.25),
            Triplet::new(1, 1, 0.75),
        ];
        let dense = DenseMatrix::from_triplets(2, &triplets, false);
        let sparse = SparseMatrix::from_triplets(2, &triplets, false);
        let dense_squared = Mcl::get_exp_dense(&dense, 2.0);
        let sparse_squared = Mcl::get_exp_sparse(&sparse, 2.0, 2).unwrap();
        for entry in sparse_squared.triplets() {
            close(dense_squared[(entry.row, entry.col)], entry.value);
        }
        let dense_gamma = Mcl::get_gamma_dense(&dense_squared, 2.0);
        let sparse_gamma = Mcl::get_gamma_sparse(&sparse_squared, 2.0, 2);
        for entry in sparse_gamma.triplets() {
            close(dense_gamma[(entry.row, entry.col)], entry.value);
        }
    }

    #[test]
    fn fractional_dense_power_handles_diagonal_matrix() {
        let matrix = DenseMatrix::from_triplets(
            2,
            &[Triplet::new(0, 0, 4.0), Triplet::new(1, 1, 9.0)],
            false,
        );
        let root = Mcl::get_exp_dense(&matrix, 0.5);
        close(root[(0, 0)], 2.0);
        close(root[(1, 1)], 3.0);
        close(root[(0, 1)], 0.0);
    }

    #[test]
    fn graph_defaults_and_matrix_constructors_preserve_source_rules() {
        let graph = get_graph_handle(2, &[], true, None, None);
        assert_eq!(graph.similarity, "normalized_bitscore_global");
        assert_eq!(graph.threshold, Some(50.0));
        assert_eq!(
            Mcl::get_description(),
            "Markov clustering according to doi:10.1137/040608635"
        );

        let mut entries = vec![Triplet::new(0, 1, 3.0)];
        let dense = Mcl::get_dense_matrix_and_clear(&[4, 8], &mut entries, true);
        assert!(entries.is_empty());
        assert_eq!(dense[(0, 1)], 3.0);
        assert_eq!(dense[(1, 0)], 3.0);
    }

    #[test]
    fn run_clusters_components_and_labels_singletons() {
        let ids = vec![
            "a description".to_owned(),
            "b".to_owned(),
            "c".to_owned(),
            "d".to_owned(),
        ];
        let edges = vec![
            Triplet::new(0, 0, 1.0),
            Triplet::new(0, 1, 0.9),
            Triplet::new(1, 1, 1.0),
            Triplet::new(2, 2, 1.0),
        ];
        let graph = get_graph_handle(ids.len(), &edges, true, None, None);
        let mut mcl = Mcl::default();
        let result = mcl
            .run(
                &ids,
                &graph,
                MclConfig {
                    max_iter: 20,
                    sparsity_switch: 1.0,
                    ..MclConfig::default()
                },
            )
            .unwrap();
        assert_eq!(result.assignments.len(), 4);
        assert_eq!(result.assignments[0].sequence_id, "a");
        assert_eq!(result.assignments[3].label, NodeLabel::Singleton);
        assert_eq!(result.statistics.components, 3);
        assert!(result.statistics.clusters >= 3);
    }

    #[test]
    fn sparse_fractional_expansion_matches_source_error() {
        let matrix = SparseMatrix::from_triplets(1, &[Triplet::new(0, 0, 1.0)], false);
        assert_eq!(
            Mcl::get_exp_sparse(&matrix, 1.5, 1).unwrap_err(),
            "Eigen does not provide an eigenvalue solver for sparse matrices"
        );
    }
}
