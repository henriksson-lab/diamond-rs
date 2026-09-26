//! Mirrored boundary for `diamond/src/stats/blast/matrix_adjust.cpp`.
//!
//! The generic implementations remain available from [`crate::stats::target_freq`]
//! and [`crate::stats::linear_algebra`]. This module supplies the fixed
//! 20-letter entry points used by the optimized upstream implementation.

use crate::stats::linear_algebra;
use crate::stats::target_freq;

pub use crate::stats::target_freq::{blast_optimize_target_frequencies, ReNewtonSystem};

pub const ALPHABET_SIZE: usize = 20;
pub const MATRIX_SIZE: usize = ALPHABET_SIZE * ALPHABET_SIZE;
pub const LINEAR_CONSTRAINTS: usize = 2 * ALPHABET_SIZE - 1;
pub const NEWTON_SIZE: usize = LINEAR_CONSTRAINTS + 1;

/// Safe Rust mapping of C++ `Nlm_DenseMatrixNew`.
pub fn nlm_dense_matrix_new(nrows: usize, ncols: usize) -> Vec<Vec<f64>> {
    vec![vec![0.0; ncols]; nrows]
}

/// Safe Rust mapping of C++ `Nlm_LtriangMatrixNew`.
pub fn nlm_ltriang_matrix_new(size: usize) -> Vec<Vec<f64>> {
    (1..=size).map(|width| vec![0.0; width]).collect()
}

/// Compatibility mapping of C++ `Nlm_DenseMatrixFree`; ownership normally
/// makes this unnecessary, but clearing also mirrors the observable nulling.
pub fn nlm_dense_matrix_free(matrix: &mut Vec<Vec<f64>>) {
    matrix.clear();
}

pub fn nlm_factor_ltriang_pos_def(matrix: &mut [Vec<f64>], size: usize) {
    linear_algebra::factor_ltriang_pos_def(matrix, size);
}

pub fn nlm_solve_ltriang_pos_def(solution: &mut [f64], matrix: &[Vec<f64>], size: usize) {
    linear_algebra::solve_ltriang_pos_def(solution, size, matrix);
}

pub fn nlm_euclidean_norm(vector: &[f64], size: usize) -> f64 {
    linear_algebra::euclidean_norm(vector, size)
}

pub fn nlm_add_vectors(output: &mut [f64], size: usize, alpha: f64, input: &[f64]) {
    linear_algebra::add_vectors(output, size, alpha, input);
}

pub fn nlm_step_bound(values: &[f64], size: usize, step: &[f64], maximum: f64) -> f64 {
    linear_algebra::step_bound(values, size, step, maximum)
}

pub fn scaled_symmetric_product_a20(matrix: &mut [Vec<f64>], diagonal: &[f64]) {
    target_freq::scaled_symmetric_product_a(matrix, diagonal, ALPHABET_SIZE);
}

pub fn multiply_by_a20(beta: f64, output: &mut [f64], alpha: f64, input: &[f64]) {
    target_freq::multiply_by_a(beta, output, ALPHABET_SIZE, alpha, input);
}

pub fn multiply_by_a_transpose20(beta: f64, output: &mut [f64], alpha: f64, input: &[f64]) {
    target_freq::multiply_by_a_transpose(beta, output, ALPHABET_SIZE, alpha, input);
}

pub fn residuals_linear_constraints20(
    residuals: &mut [f64],
    values: &[f64],
    row_sums: &[f64],
    column_sums: &[f64],
) {
    target_freq::residuals_linear_constraints(
        residuals,
        ALPHABET_SIZE,
        values,
        row_sums,
        column_sums,
    );
}

pub fn dual_residuals20(residuals: &mut [f64], gradients: &[Vec<f64>], z: &[f64]) {
    target_freq::dual_residuals(residuals, ALPHABET_SIZE, gradients, z, true);
}

#[allow(clippy::too_many_arguments)]
pub fn calculate_residuals20(
    residuals_x: &mut [f64],
    residuals_z: &mut [f64],
    values: &[f64],
    gradients: &[Vec<f64>],
    row_sums: &[f64],
    column_sums: &[f64],
    x: &[f64],
    z: &[f64],
    target_relative_entropy: f64,
) -> f64 {
    target_freq::calculate_residuals(
        residuals_x,
        ALPHABET_SIZE,
        residuals_z,
        values,
        gradients,
        row_sums,
        column_sums,
        x,
        z,
        true,
        target_relative_entropy,
    )
}

pub fn evaluate_re_functions20(
    values: &mut [f64],
    gradients: &mut [Vec<f64>],
    x: &[f64],
    q: &[f64],
    scores: &[f64],
) {
    target_freq::evaluate_re_functions(values, gradients, ALPHABET_SIZE, x, q, scores, true);
}

pub fn compute_scores_from_probs20(
    scores: &mut [f64],
    target_frequencies: &[f64],
    row_frequencies: &[f64],
    column_frequencies: &[f64],
) {
    target_freq::compute_scores_from_probs(
        scores,
        ALPHABET_SIZE,
        target_frequencies,
        row_frequencies,
        column_frequencies,
    );
}

pub type NewtonSys20 = ReNewtonSystem;

pub fn factor_newton20(
    system: &mut NewtonSys20,
    x: &[f64],
    z: &[f64],
    gradients: &[Vec<f64>],
    workspace: &mut [f64],
) {
    target_freq::factor_re_newton_system(system, x, z, gradients, true, workspace);
}

pub fn solve_newton20(
    step_x: &mut [f64],
    step_z: &mut [f64],
    system: &NewtonSys20,
    workspace: &mut [f64],
) {
    target_freq::solve_re_newton_system(step_x, step_z, system, workspace);
}

#[allow(clippy::too_many_arguments)]
pub fn new_optimize_target_frequencies(
    x: &mut [f64],
    alphabet_size: usize,
    q: &[f64],
    row_sums: &[f64],
    column_sums: &[f64],
    relative_entropy: f64,
    tolerance: f64,
    max_iterations: i32,
) -> (i32, i32) {
    target_freq::new_optimize_target_frequencies(
        x,
        alphabet_size,
        q,
        row_sums,
        column_sums,
        relative_entropy,
        tolerance,
        max_iterations,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allocation_mappings_have_cpp_shapes() {
        let mut dense = nlm_dense_matrix_new(3, 4);
        assert_eq!(dense, vec![vec![0.0; 4]; 3]);
        nlm_dense_matrix_free(&mut dense);
        assert!(dense.is_empty());

        let lower = nlm_ltriang_matrix_new(4);
        assert_eq!(lower.iter().map(Vec::len).collect::<Vec<_>>(), [1, 2, 3, 4]);
        assert!(lower.iter().flatten().all(|&value| value == 0.0));
    }

    #[test]
    fn fixed_twenty_constraint_products_are_exact() {
        let mut input = vec![0.0; MATRIX_SIZE];
        input[0] = 1.25;
        input[ALPHABET_SIZE + 2] = -0.5;
        input[19 * ALPHABET_SIZE + 19] = 2.0;

        let mut output = vec![3.0; LINEAR_CONSTRAINTS];
        multiply_by_a20(0.0, &mut output, 2.0, &input);
        assert_eq!(output[0], 2.5);
        assert_eq!(output[2], -1.0);
        assert_eq!(output[19], 4.0);
        assert_eq!(output[20], -1.0);
        assert_eq!(output[38], 4.0);
        assert!(output
            .iter()
            .enumerate()
            .all(|(index, &value)| matches!(index, 0 | 2 | 19 | 20 | 38) || value == 0.0));

        let mut transpose = vec![0.0; MATRIX_SIZE];
        multiply_by_a_transpose20(0.0, &mut transpose, 1.0, &output);
        assert_eq!(transpose[0], 2.5);
        assert_eq!(transpose[ALPHABET_SIZE + 2], -2.0);
        assert_eq!(transpose[19 * ALPHABET_SIZE + 19], 8.0);
    }

    #[test]
    fn row_constraint_uses_cpp_per_element_rounding_order() {
        let mut input = vec![0.0; MATRIX_SIZE];
        for value in &mut input[ALPHABET_SIZE..2 * ALPHABET_SIZE] {
            *value = 0.1;
        }
        let mut output = vec![0.0; LINEAR_CONSTRAINTS];
        output[ALPHABET_SIZE] = 1.0e16;
        multiply_by_a20(1.0, &mut output, 1.0, &input);
        assert_eq!(output[ALPHABET_SIZE], 1.0e16);
    }

    #[test]
    fn uniform_fixed_point_matches_specialized_entry_point() {
        let row = vec![1.0 / ALPHABET_SIZE as f64; ALPHABET_SIZE];
        let column = row.clone();
        let q = vec![1.0 / MATRIX_SIZE as f64; MATRIX_SIZE];
        let mut x = vec![f64::NAN; MATRIX_SIZE];
        let (status, iterations) = new_optimize_target_frequencies(
            &mut x,
            ALPHABET_SIZE,
            &q,
            &row,
            &column,
            0.0,
            1.0e-12,
            0,
        );
        assert_eq!((status, iterations), (1, 1));
        assert_eq!(x, q);
    }
}
