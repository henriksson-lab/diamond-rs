//! Safe Rust boundary for NCBI BLAST `nlm_linear_algebra.cpp`.
//!
//! Numerical kernels reuse the audited implementation in
//! [`crate::stats::linear_algebra`]. Matrix allocation uses owned row vectors;
//! `Option` represents the C allocation/null state and Rust ownership replaces
//! the raw contiguous buffers.

use crate::stats::linear_algebra;

pub type DenseMatrix = Vec<Vec<f64>>;
pub type Int4Matrix = Vec<Vec<i32>>;

fn zeroed_matrix<T: Clone + Default>(nrows: usize, ncols: usize) -> Option<Vec<Vec<T>>> {
    let mut matrix = Vec::new();
    matrix.try_reserve_exact(nrows).ok()?;
    for _ in 0..nrows {
        let mut row = Vec::new();
        row.try_reserve_exact(ncols).ok()?;
        row.resize(ncols, T::default());
        matrix.push(row);
    }
    Some(matrix)
}

/// Allocate a dense matrix. Safe initialization replaces the C routine's
/// uninitialized scalar buffer; valid C callers must initialize before read.
pub fn nlm_dense_matrix_new(nrows: usize, ncols: usize) -> Option<DenseMatrix> {
    zeroed_matrix(nrows, ncols)
}

/// Allocate a zero-initialized packed lower-triangular matrix.
pub fn nlm_ltriang_matrix_new(n: usize) -> Option<DenseMatrix> {
    let mut matrix = Vec::new();
    matrix.try_reserve_exact(n).ok()?;
    for width in 1..=n {
        let mut row = Vec::new();
        row.try_reserve_exact(width).ok()?;
        row.resize(width, 0.0);
        matrix.push(row);
    }
    Some(matrix)
}

/// Drop a dense or lower-triangular allocation and leave the handle null.
pub fn nlm_dense_matrix_free(matrix: &mut Option<DenseMatrix>) {
    *matrix = None;
}

/// Allocate the header's `Int4` matrix (`Int4` is a signed 32-bit integer).
pub fn nlm_int4_matrix_new(nrows: usize, ncols: usize) -> Option<Int4Matrix> {
    zeroed_matrix(nrows, ncols)
}

/// Drop an `Int4` allocation and leave the handle null.
pub fn nlm_int4_matrix_free(matrix: &mut Option<Int4Matrix>) {
    *matrix = None;
}

/// In-place Cholesky factorization of a positive-definite lower triangle.
pub fn nlm_factor_ltriang_pos_def(matrix: &mut [Vec<f64>], n: usize) {
    linear_algebra::factor_ltriang_pos_def(matrix, n);
}

/// Solve `L L^T x = b` in place.
pub fn nlm_solve_ltriang_pos_def(x: &mut [f64], n: usize, matrix: &[Vec<f64>]) {
    linear_algebra::solve_ltriang_pos_def(x, n, matrix);
}

/// Stable Euclidean norm using the source routine's scaled accumulation.
pub fn nlm_euclidean_norm(vector: &[f64], n: usize) -> f64 {
    linear_algebra::euclidean_norm(vector, n)
}

/// Apply `y += alpha * x` to the first `n` entries.
pub fn nlm_add_vectors(y: &mut [f64], n: usize, alpha: f64, x: &[f64]) {
    linear_algebra::add_vectors(y, n, alpha, x);
}

/// Largest step in `[0, maximum]` that keeps `x + alpha * step_x` nonnegative.
pub fn nlm_step_bound(x: &[f64], n: usize, step_x: &[f64], maximum: f64) -> f64 {
    linear_algebra::step_bound(x, n, step_x, maximum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allocation_shapes_and_free_nulling_cover_all_matrix_kinds() {
        let mut dense = nlm_dense_matrix_new(3, 4);
        assert_eq!(
            dense
                .as_ref()
                .unwrap()
                .iter()
                .map(Vec::len)
                .collect::<Vec<_>>(),
            [4, 4, 4]
        );
        dense.as_mut().unwrap()[2][3] = 7.5;
        nlm_dense_matrix_free(&mut dense);
        assert!(dense.is_none());

        let mut lower = nlm_ltriang_matrix_new(4);
        assert_eq!(
            lower
                .as_ref()
                .unwrap()
                .iter()
                .map(Vec::len)
                .collect::<Vec<_>>(),
            [1, 2, 3, 4]
        );
        assert!(lower
            .as_ref()
            .unwrap()
            .iter()
            .flatten()
            .all(|&value| value == 0.0));
        nlm_dense_matrix_free(&mut lower);
        assert!(lower.is_none());

        let mut integers = nlm_int4_matrix_new(2, 3);
        assert_eq!(integers.as_ref().unwrap(), &vec![vec![0_i32; 3]; 2]);
        nlm_int4_matrix_free(&mut integers);
        assert!(integers.is_none());

        assert_eq!(nlm_dense_matrix_new(0, 5), Some(Vec::new()));
        assert_eq!(nlm_ltriang_matrix_new(0), Some(Vec::new()));
    }

    #[test]
    fn factor_and_solve_match_a_known_three_by_three_system() {
        let mut lower = vec![vec![25.0], vec![15.0, 18.0], vec![-5.0, 0.0, 11.0]];
        nlm_factor_ltriang_pos_def(&mut lower, 3);
        assert_eq!(lower, vec![vec![5.0], vec![3.0, 3.0], vec![-1.0, 1.0, 3.0]]);

        let mut rhs = [40.0, 51.0, 28.0];
        nlm_solve_ltriang_pos_def(&mut rhs, 3, &lower);
        for (actual, expected) in rhs.into_iter().zip([1.0, 2.0, 3.0]) {
            assert!((actual - expected).abs() < 1.0e-14);
        }
    }

    #[test]
    fn norm_retains_scale_stability_and_prefix_semantics() {
        assert_eq!(nlm_euclidean_norm(&[], 0), 0.0);
        assert_eq!(nlm_euclidean_norm(&[3.0, 4.0, 99.0], 2), 5.0);
        let large = nlm_euclidean_norm(&[3.0e200, 4.0e200], 2);
        assert!(large.is_finite());
        assert!((large / 5.0e200 - 1.0).abs() < 1.0e-15);
        let small = nlm_euclidean_norm(&[3.0e-200, 4.0e-200], 2);
        assert!(small > 0.0);
        assert!((small / 5.0e-200 - 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn vector_add_and_step_bound_cover_prefix_and_boundary_paths() {
        let mut y = [1.0, 2.0, 30.0];
        nlm_add_vectors(&mut y, 2, -0.5, &[4.0, 6.0, 100.0]);
        assert_eq!(y, [-1.0, -1.0, 30.0]);

        assert_eq!(nlm_step_bound(&[2.0, 4.0], 2, &[-1.0, -4.0], 10.0), 1.0);
        assert_eq!(nlm_step_bound(&[2.0, 4.0], 2, &[1.0, 0.0], 3.0), 3.0);
        assert_eq!(nlm_step_bound(&[0.0], 1, &[-2.0], 3.0), 0.0);
    }
}
