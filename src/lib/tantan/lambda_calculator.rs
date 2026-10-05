//! Faithful translation of `diamond/src/lib/tantan/LambdaCalculator.cc`.
//!
//! This is the calculator used by DIAMOND's active repeat-masking path.  It is
//! intentionally separate from `masking::lambda`, which mirrors the unrelated
//! standalone helper in `diamond/src/masking/lambda.cpp`.

#[cfg(windows)]
const RAND_MAX_F64: f64 = 32_767.0;
#[cfg(not(windows))]
const RAND_MAX_F64: f64 = 2_147_483_647.0;

unsafe extern "C" {
    fn rand() -> i32;
}

#[derive(Debug, Clone, PartialEq)]
pub struct LambdaCalculator {
    lambda: f64,
    letter_probs1: Vec<f64>,
    letter_probs2: Vec<f64>,
}

impl Default for LambdaCalculator {
    fn default() -> Self {
        Self::new()
    }
}

impl LambdaCalculator {
    pub fn new() -> Self {
        Self {
            lambda: -1.0,
            letter_probs1: Vec::new(),
            letter_probs2: Vec::new(),
        }
    }

    pub fn set_bad(&mut self) {
        self.lambda = -1.0;
        self.letter_probs1.clear();
        self.letter_probs2.clear();
    }

    pub fn is_bad(&self) -> bool {
        self.lambda < 0.0
    }

    pub fn lambda(&self) -> f64 {
        self.lambda
    }

    pub fn letter_probs1(&self) -> Option<&[f64]> {
        (!self.is_bad()).then_some(&self.letter_probs1)
    }

    pub fn letter_probs2(&self) -> Option<&[f64]> {
        (!self.is_bad()).then_some(&self.letter_probs2)
    }

    pub fn calculate(&mut self, matrix: &[Vec<i32>]) {
        self.set_bad();
        let n = matrix.len();
        if matrix.iter().any(|row| row.len() != n) {
            return;
        }
        self.lambda = calculate_lambda(
            matrix,
            &mut self.letter_probs1,
            &mut self.letter_probs2,
            1000,
            100,
            1e-6,
        );
    }

    pub fn calculate_flat_i8(&mut self, scores: &[i8], stride: usize, n: usize) {
        if stride < n || scores.len() < stride.saturating_mul(n) {
            self.set_bad();
            return;
        }
        let matrix = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| i32::from(scores[i * stride + j]))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        self.calculate(&matrix);
    }
}

fn round_to_few_digits(value: f64) -> f64 {
    // `%g` without an explicit precision uses six significant digits.
    if value == 0.0 || !value.is_finite() {
        return value;
    }
    let scale = 10_f64.powf(5.0 - value.abs().log10().floor());
    (value * scale).round() / scale
}

fn max_index(matrix: &[Vec<f64>], column: usize) -> usize {
    let mut maximum = -f64::MAX;
    let mut maximum_index = column;
    for (row, values) in matrix.iter().enumerate().skip(column) {
        let value = values[column].abs();
        // Upstream uses `>` rather than `>=`, so retain the first row on a
        // tie. This can change the floating-point path through the LU solve.
        if value > maximum {
            maximum = value;
            maximum_index = row;
        }
    }
    maximum_index
}

fn lu_pivoting(matrix: &mut [Vec<f64>], indices: &mut [usize]) -> bool {
    for (i, index) in indices.iter_mut().enumerate() {
        *index = i;
    }
    for i in 0..matrix.len() {
        let pivot_row = max_index(matrix, i);
        if matrix[pivot_row][i].abs() < 1e-10 {
            return false;
        }
        matrix.swap(i, pivot_row);
        indices.swap(i, pivot_row);
        matrix[i][i] = 1.0 / matrix[i][i];
        for j in i + 1..matrix.len() {
            matrix[j][i] *= matrix[i][i];
            for k in i + 1..matrix.len() {
                matrix[j][k] -= matrix[j][i] * matrix[i][k];
            }
        }
    }
    true
}

fn solve_p(matrix: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let n = matrix.len();
    let mut y = vec![0.0; n];
    for i in 0..n {
        y[i] = rhs[i];
        for j in 0..i {
            y[i] -= matrix[i][j] * y[j];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = y[i];
        for j in i + 1..n {
            x[i] -= matrix[i][j] * x[j];
        }
        x[i] *= matrix[i][i];
    }
    x
}

fn invert(mut matrix: Vec<Vec<f64>>) -> Option<Vec<Vec<f64>>> {
    let n = matrix.len();
    let mut indices = vec![0; n];
    if !lu_pivoting(&mut matrix, &mut indices) {
        return None;
    }
    let mut transposed_basis = vec![vec![0.0; n]; n];
    for i in 0..n {
        transposed_basis[indices[i]][i] = 1.0;
    }
    let mut inverse = (0..n)
        .map(|i| solve_p(&matrix, &transposed_basis[i]))
        .collect::<Vec<_>>();
    for i in 0..n {
        for j in 0..i {
            let value = inverse[i][j];
            inverse[i][j] = inverse[j][i];
            inverse[j][i] = value;
        }
    }
    Some(inverse)
}

fn inverse_at(matrix: &[Vec<i32>], tau: f64) -> Option<Vec<Vec<f64>>> {
    let exponential = matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|&score| (tau * f64::from(score)).exp())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    invert(exponential)
}

fn calculate_inv_sum(matrix: &[Vec<i32>], tau: f64) -> Option<f64> {
    Some(inverse_at(matrix, tau)?.iter().flatten().sum())
}

fn find_upper_bound(matrix: &[Vec<i32>]) -> Option<f64> {
    let n = matrix.len();
    if n == 0 {
        return None;
    }
    let mut row_max_min = f64::MAX;
    let mut column_max_min = f64::MAX;
    let mut zero_rows = 0;
    let mut zero_columns = 0;

    for row in matrix {
        let row_max = f64::from(*row.iter().max()?);
        let row_min = f64::from(*row.iter().min()?);
        if row_max == 0.0 && row_min == 0.0 {
            zero_rows += 1;
        } else if row_max <= 0.0 || row_min >= 0.0 {
            return None;
        } else {
            row_max_min = row_max_min.min(row_max);
        }
    }
    for j in 0..n {
        let column_max = f64::from((0..n).map(|i| matrix[i][j]).max()?);
        let column_min = f64::from((0..n).map(|i| matrix[i][j]).min()?);
        if column_max == 0.0 && column_min == 0.0 {
            zero_columns += 1;
        } else if column_max <= 0.0 || column_min >= 0.0 {
            return None;
        } else {
            column_max_min = column_max_min.min(column_max);
        }
    }
    if zero_rows == n {
        return None;
    }
    Some(if row_max_min > column_max_min {
        1.1 * ((n - zero_rows) as f64).ln() / row_max_min
    } else {
        1.1 * ((n - zero_columns) as f64).ln() / column_max_min
    })
}

fn random_between(lower: f64, upper: f64) -> f64 {
    // The source uses C `rand()` and intentionally does not seed it here.
    lower + (upper - lower) * f64::from(unsafe { rand() }) / RAND_MAX_F64
}

fn check_lambda(
    matrix: &[Vec<i32>],
    lambda: f64,
    letter_probs1: &mut Vec<f64>,
    letter_probs2: &mut Vec<f64>,
) -> bool {
    let Some(inverse) = inverse_at(matrix, lambda) else {
        letter_probs1.clear();
        letter_probs2.clear();
        return false;
    };
    let n = matrix.len();
    letter_probs1.clear();
    letter_probs2.clear();
    for row in &inverse {
        let probability: f64 = row.iter().sum();
        if !(0.0..=1.0).contains(&probability) {
            letter_probs2.clear();
            return false;
        }
        letter_probs2.push(round_to_few_digits(probability));
    }
    for j in 0..n {
        let probability: f64 = inverse.iter().map(|row| row[j]).sum();
        if !(0.0..=1.0).contains(&probability) {
            letter_probs1.clear();
            letter_probs2.clear();
            return false;
        }
        letter_probs1.push(round_to_few_digits(probability));
    }
    true
}

#[allow(clippy::too_many_arguments)]
fn binary_search(
    matrix: &[Vec<i32>],
    lower_bound: f64,
    upper_bound: f64,
    letter_probs1: &mut Vec<f64>,
    letter_probs2: &mut Vec<f64>,
    lambda: &mut f64,
    max_iterations: usize,
) -> bool {
    let (mut left, mut right, mut left_sum, mut right_sum) = (0.0, 0.0, 0.0, 0.0);
    let mut iteration = 0;
    while iteration < max_iterations
        && (left >= right
            || (left_sum < 1.0 && right_sum < 1.0)
            || (left_sum > 1.0 && right_sum > 1.0))
    {
        left = random_between(lower_bound, upper_bound);
        right = random_between(lower_bound, upper_bound);
        match (
            calculate_inv_sum(matrix, left),
            calculate_inv_sum(matrix, right),
        ) {
            (Some(l), Some(r)) => {
                left_sum = l;
                right_sum = r;
            }
            _ => {
                left = 0.0;
                right = 0.0;
            }
        }
        iteration += 1;
    }
    if iteration >= max_iterations {
        return false;
    }

    while left_sum != 1.0
        && right_sum != 1.0
        && (left + right) / 2.0 != left
        && (left + right) / 2.0 != right
    {
        let middle = (left + right) / 2.0;
        let Some(middle_sum) = calculate_inv_sum(matrix, middle) else {
            return false;
        };
        if middle_sum.abs() >= f64::MAX {
            return false;
        }
        if (left_sum < 1.0 && middle_sum >= 1.0) || (left_sum > 1.0 && middle_sum <= 1.0) {
            right = middle;
            right_sum = middle_sum;
        } else if (right_sum < 1.0 && middle_sum >= 1.0) || (right_sum > 1.0 && middle_sum <= 1.0) {
            left = middle;
            left_sum = middle_sum;
        } else {
            return false;
        }
    }

    let candidate = if (left_sum - 1.0).abs() < (right_sum - 1.0).abs() {
        left
    } else {
        right
    };
    if check_lambda(matrix, candidate, letter_probs1, letter_probs2) {
        *lambda = candidate;
        true
    } else {
        false
    }
}

fn calculate_lambda(
    matrix: &[Vec<i32>],
    letter_probs1: &mut Vec<f64>,
    letter_probs2: &mut Vec<f64>,
    max_iterations: usize,
    max_boundary_search_iterations: usize,
    lower_bound_ratio: f64,
) -> f64 {
    let Some(upper_bound) = find_upper_bound(matrix) else {
        return -1.0;
    };
    let lower_bound = upper_bound * lower_bound_ratio;
    let mut lambda = -1.0;
    for _ in 0..max_iterations {
        if binary_search(
            matrix,
            lower_bound,
            upper_bound,
            letter_probs1,
            letter_probs2,
            &mut lambda,
            max_boundary_search_iterations,
        ) {
            break;
        }
    }
    lambda
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::AMINO_ACID_COUNT;
    use crate::stats::matrices::{BLOSUM62, PAM250};

    #[test]
    fn standard_matrices_match_upstream_states() {
        let mut calculator = LambdaCalculator::new();
        calculator.calculate_flat_i8(&BLOSUM62.scores, AMINO_ACID_COUNT, 20);
        assert!(!calculator.is_bad());
        assert!((calculator.lambda() - 0.324_032).abs() < 1e-5);

        calculator.calculate_flat_i8(&PAM250.scores, AMINO_ACID_COUNT, 20);
        assert!(calculator.is_bad());
        assert_eq!(calculator.lambda(), -1.0);
        assert!(calculator.letter_probs1().is_none());
        assert!(calculator.letter_probs2().is_none());
    }

    #[test]
    fn rejects_invalid_shapes_and_signs() {
        let mut calculator = LambdaCalculator::new();
        calculator.calculate(&[vec![1, -1]]);
        assert!(calculator.is_bad());
        calculator.calculate(&[vec![1, 1], vec![1, 1]]);
        assert!(calculator.is_bad());
    }
}
