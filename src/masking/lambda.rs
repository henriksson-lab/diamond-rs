//! Scale-factor and implicit letter-probability calculation for score matrices.
//!
//! This is the Rust counterpart of `diamond/src/masking/lambda.cpp`.  The
//! implementation intentionally remains self-contained: the original accepts
//! any square integer matrix, rather than a DIAMOND `ScoreMatrix`.

/// Result of [`LambdaCalculator::compute`].
#[derive(Debug, Clone, PartialEq)]
pub struct Result {
    pub ok: bool,
    pub lambda: f64,
    /// Column sums of `inverse(exp(lambda * scores))`.
    pub left_probs: Vec<f64>,
    /// Row sums of `inverse(exp(lambda * scores))`.
    pub right_probs: Vec<f64>,
    pub reason: String,
}

impl Default for Result {
    fn default() -> Self {
        Self {
            ok: false,
            lambda: -1.0,
            left_probs: Vec::new(),
            right_probs: Vec::new(),
            reason: String::new(),
        }
    }
}

/// Computes the scale factor and letter probabilities implicit in a score
/// matrix.
#[derive(Debug, Clone, Copy, Default)]
pub struct LambdaCalculator;

impl LambdaCalculator {
    pub fn new() -> Self {
        Self
    }

    /// Matches C++ `LambdaCalculator::compute` with its default arguments.
    pub fn compute(&self, scores: &[Vec<i32>]) -> Result {
        self.compute_with_options(scores, 1000, 100, 1e-6)
    }

    /// Matches C++ `LambdaCalculator::compute`, exposing the optional C++
    /// arguments explicitly because Rust has no default function arguments.
    pub fn compute_with_options(
        &self,
        scores: &[Vec<i32>],
        max_outer_iters: usize,
        max_bracket_tries: usize,
        lb_ratio: f64,
    ) -> Result {
        let mut out = Result::default();
        let n = scores.len();
        if n == 0 {
            out.reason = "Empty matrix.".to_owned();
            return out;
        }
        if scores.iter().any(|row| row.len() != n) {
            out.reason = "Matrix must be square.".to_owned();
            return out;
        }

        let Some(upper_bound) = Self::find_upper_bound(scores) else {
            out.reason =
                "Failed to find a valid upper bound (score matrix violates sign conditions)."
                    .to_owned();
            return out;
        };
        let lower_bound = lb_ratio * upper_bound;
        let mut rng = Mt19937_64::new(0xC001_D00D);

        for _ in 0..max_outer_iters {
            let mut bracket = None;
            for _ in 0..max_bracket_tries {
                let mut left = rng.uniform(lower_bound, upper_bound);
                let mut right = rng.uniform(lower_bound, upper_bound);
                if left > right {
                    std::mem::swap(&mut left, &mut right);
                }
                let Some(f_left) = Self::inv_sum(scores, left) else {
                    continue;
                };
                let Some(f_right) = Self::inv_sum(scores, right) else {
                    continue;
                };
                if (f_left <= 0.0 && f_right >= 0.0) || (f_left >= 0.0 && f_right <= 0.0) {
                    bracket = Some((left, right, f_left, f_right));
                    break;
                }
            }

            let Some((mut left, mut right, mut f_left, mut f_right)) = bracket else {
                continue;
            };
            let mut safeguard = 200;
            while safeguard > 0 {
                safeguard -= 1;
                let mid = 0.5 * (left + right);
                let Some(f_mid) = Self::inv_sum(scores, mid) else {
                    out.reason = "Singular/unstable matrix during bisection.".to_owned();
                    break;
                };

                if f_mid.abs() <= 1e-12 {
                    if self.finalize(scores, mid, &mut out) {
                        return out;
                    }
                    break;
                }

                if !(mid > left && mid < right) {
                    break;
                }
                if (f_mid < 0.0 && f_left < 0.0) || (f_mid > 0.0 && f_left > 0.0) {
                    left = mid;
                    f_left = f_mid;
                } else {
                    right = mid;
                    f_right = f_mid;
                }

                if (right - left).abs() <= 1.0_f64.max(mid.abs()) * 1e-12 {
                    let candidate = if f_left.abs() < f_right.abs() {
                        left
                    } else {
                        right
                    };
                    if self.finalize(scores, candidate, &mut out) {
                        return out;
                    }
                    break;
                }
            }
        }

        if !out.ok && out.reason.is_empty() {
            out.reason = "Failed to bracket and solve for lambda.".to_owned();
        }
        out
    }

    /// Matches C++ `tidy` (eight significant decimal digits).  Its input here
    /// is a validated probability, so decimal significant-digit rounding is
    /// equivalent to the source's stream-format-and-parse operation.
    fn tidy(x: f64) -> f64 {
        if x == 0.0 {
            return x;
        }
        let scale = 10_f64.powf(7.0 - x.abs().log10().floor());
        (x * scale).round() / scale
    }

    /// Matches C++ `buildExpMatrix`.
    fn build_exp_matrix(scores: &[Vec<i32>], lambda: f64) -> Matrix {
        let n = scores.len();
        let mut matrix = Matrix::new(n);
        for (i, row) in scores.iter().enumerate() {
            for (j, &score) in row.iter().enumerate() {
                matrix[(i, j)] = (lambda * f64::from(score)).exp();
            }
        }
        matrix
    }

    /// Matches C++ `luDecompose`.
    fn lu_decompose(matrix: &mut Matrix, pivots: &mut Vec<usize>, eps: f64) -> bool {
        let n = matrix.n;
        pivots.clear();
        pivots.extend(0..n);

        for k in 0..n {
            let mut pivot = k;
            let mut max_abs = matrix[(k, k)].abs();
            for i in k + 1..n {
                let value = matrix[(i, k)].abs();
                if value > max_abs {
                    max_abs = value;
                    pivot = i;
                }
            }
            if max_abs < eps {
                return false;
            }
            if pivot != k {
                for j in 0..n {
                    matrix.data.swap(k * n + j, pivot * n + j);
                }
                pivots.swap(k, pivot);
            }

            for i in k + 1..n {
                matrix[(i, k)] /= matrix[(k, k)];
                let lower = matrix[(i, k)];
                for j in k + 1..n {
                    matrix[(i, j)] -= lower * matrix[(k, j)];
                }
            }
        }
        true
    }

    /// Matches C++ `luSolve`.
    fn lu_solve(lu: &Matrix, pivots: &[usize], b: &[f64], x: &mut Vec<f64>) {
        let n = lu.n;
        x.clear();
        x.resize(n, 0.0);
        let mut y = vec![0.0; n];
        for i in 0..n {
            y[i] = b[pivots[i]];
        }
        for i in 0..n {
            let mut sum = y[i];
            for j in 0..i {
                sum -= lu[(i, j)] * y[j];
            }
            y[i] = sum;
        }
        for i in (0..n).rev() {
            let mut sum = y[i];
            for j in i + 1..n {
                sum -= lu[(i, j)] * x[j];
            }
            x[i] = sum / lu[(i, i)];
        }
    }

    /// Matches C++ `invert`.
    fn invert(matrix: &Matrix) -> Option<Matrix> {
        let mut lu = matrix.clone();
        let mut pivots = Vec::new();
        if !Self::lu_decompose(&mut lu, &mut pivots, 1e-12) {
            return None;
        }

        let n = matrix.n;
        let mut inverse = Matrix::new(n);
        let mut basis = vec![0.0; n];
        let mut column = Vec::new();
        for j in 0..n {
            basis.fill(0.0);
            basis[j] = 1.0;
            Self::lu_solve(&lu, &pivots, &basis, &mut column);
            for i in 0..n {
                inverse[(i, j)] = column[i];
            }
        }
        Some(inverse)
    }

    /// Matches C++ `invSum`, including subtraction of one for the root
    /// function.
    fn inv_sum(scores: &[Vec<i32>], lambda: f64) -> Option<f64> {
        let matrix = Self::build_exp_matrix(scores, lambda);
        let inverse = Self::invert(&matrix)?;
        let value = inverse.data.iter().copied().sum::<f64>() - 1.0;
        value.is_finite().then_some(value)
    }

    /// Matches C++ `finalize`.
    fn finalize(&self, scores: &[Vec<i32>], lambda: f64, out: &mut Result) -> bool {
        let matrix = Self::build_exp_matrix(scores, lambda);
        let Some(inverse) = Self::invert(&matrix) else {
            out.reason = "Matrix inversion failed at finalization.".to_owned();
            return false;
        };
        let n = matrix.n;
        let mut row_sums = vec![0.0; n];
        let mut column_sums = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                row_sums[i] += inverse[(i, j)];
                column_sums[j] += inverse[(i, j)];
            }
        }
        if row_sums.iter().any(|&p| !(0.0..=1.0).contains(&p)) {
            out.reason = "Row probability outside [0,1].".to_owned();
            return false;
        }
        if column_sums.iter().any(|&p| !(0.0..=1.0).contains(&p)) {
            out.reason = "Column probability outside [0,1].".to_owned();
            return false;
        }
        row_sums.iter_mut().for_each(|p| *p = Self::tidy(*p));
        column_sums.iter_mut().for_each(|p| *p = Self::tidy(*p));

        out.ok = true;
        out.lambda = lambda;
        out.right_probs = row_sums;
        out.left_probs = column_sums;
        out.reason.clear();
        true
    }

    /// Matches C++ `findUpperBound`.
    fn find_upper_bound(scores: &[Vec<i32>]) -> Option<f64> {
        let n = scores.len();
        let mut row_max_min = f64::INFINITY;
        let mut column_max_min = f64::INFINITY;
        let mut zero_rows = 0usize;
        let mut zero_columns = 0usize;

        for row in scores {
            let row_max = *row.iter().max()?;
            let row_min = *row.iter().min()?;
            if row_max == 0 && row_min == 0 {
                zero_rows += 1;
                continue;
            }
            if row_max <= 0 || row_min >= 0 {
                return None;
            }
            row_max_min = row_max_min.min(f64::from(row_max));
        }
        for j in 0..n {
            let mut column_max = i32::MIN;
            let mut column_min = i32::MAX;
            for row in scores {
                column_max = column_max.max(row[j]);
                column_min = column_min.min(row[j]);
            }
            if column_max == 0 && column_min == 0 {
                zero_columns += 1;
                continue;
            }
            if column_max <= 0 || column_min >= 0 {
                return None;
            }
            column_max_min = column_max_min.min(f64::from(column_max));
        }
        if zero_rows == n {
            return None;
        }

        let upper_bound = if row_max_min > column_max_min {
            1.1 * ((n - zero_rows) as f64).ln() / row_max_min
        } else {
            1.1 * ((n - zero_columns) as f64).ln() / column_max_min
        };
        (upper_bound.is_finite() && upper_bound > 0.0).then_some(upper_bound)
    }
}

/// Small row-major matrix matching the source's `Mat` helper.
#[derive(Debug, Clone)]
struct Matrix {
    n: usize,
    data: Vec<f64>,
}

impl Matrix {
    fn new(n: usize) -> Self {
        Self {
            n,
            data: vec![0.0; n * n],
        }
    }
}

impl std::ops::Index<(usize, usize)> for Matrix {
    type Output = f64;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.data[i * self.n + j]
    }
}

impl std::ops::IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        &mut self.data[i * self.n + j]
    }
}

/// Minimal standard MT19937-64 engine used by the C++ source.  Keeping this
/// local avoids adding a crate solely for deterministic bracket sampling.
struct Mt19937_64 {
    state: [u64; 312],
    index: usize,
}

impl Mt19937_64 {
    fn new(seed: u64) -> Self {
        let mut state = [0; 312];
        state[0] = seed;
        for i in 1..312 {
            state[i] = 6_364_136_223_846_793_005u64
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 62))
                .wrapping_add(i as u64);
        }
        Self { state, index: 312 }
    }

    fn next_u64(&mut self) -> u64 {
        if self.index >= 312 {
            for i in 0..312 {
                let x = (self.state[i] & 0xffff_ffff_8000_0000)
                    | (self.state[(i + 1) % 312] & 0x7fff_ffff);
                let mut xa = x >> 1;
                if x & 1 != 0 {
                    xa ^= 0xB502_6F5A_A966_19E9;
                }
                self.state[i] = self.state[(i + 156) % 312] ^ xa;
            }
            self.index = 0;
        }
        let mut x = self.state[self.index];
        self.index += 1;
        x ^= (x >> 29) & 0x5555_5555_5555_5555;
        x ^= (x << 17) & 0x71D6_7FFF_EDA6_0000;
        x ^= (x << 37) & 0xFFF7_EEE0_0000_0000;
        x ^ (x >> 43)
    }

    fn uniform(&mut self, lower: f64, upper: f64) -> f64 {
        // libstdc++'s uniform_real_distribution uses generate_canonical; one
        // 64-bit engine draw supplies all 53 significant bits for `double`.
        let unit = (self.next_u64() as f64) / 18_446_744_073_709_551_616.0;
        lower + (upper - lower) * unit
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} differs from {expected} by more than {tolerance}"
        );
    }

    #[test]
    fn cpp_demo_matrix_matches_lambda_and_uniform_probabilities() {
        // Exact matrix and output shown in lambda.cpp's demo.
        let scores = vec![
            vec![2, -1, -1, -1],
            vec![-1, 2, -1, -1],
            vec![-1, -1, 2, -1],
            vec![-1, -1, -1, 2],
        ];
        let result = LambdaCalculator::new().compute(&scores);
        assert!(result.ok, "{}", result.reason);
        assert_close(result.lambda, 0.264_497_094_314, 2e-12);
        assert_eq!(result.left_probs, vec![0.25; 4]);
        assert_eq!(result.right_probs, vec![0.25; 4]);
        assert!(result.reason.is_empty());
    }

    #[test]
    fn rejects_empty_non_square_and_invalid_sign_matrices() {
        let calculator = LambdaCalculator::new();
        assert_eq!(calculator.compute(&[]).reason, "Empty matrix.");
        assert_eq!(
            calculator.compute(&[vec![1, -1]]).reason,
            "Matrix must be square."
        );
        let invalid = calculator.compute(&[vec![1, 1], vec![-1, 1]]);
        assert!(!invalid.ok);
        assert!(invalid.reason.contains("violates sign conditions"));
        let zero = calculator.compute(&[vec![0, 0], vec![0, 0]]);
        assert!(!zero.ok);
        assert!(zero.reason.contains("violates sign conditions"));
    }

    #[test]
    fn upper_bound_matches_cpp_formula() {
        let scores = vec![
            vec![2, -1, -1, -1],
            vec![-1, 2, -1, -1],
            vec![-1, -1, 2, -1],
            vec![-1, -1, -1, 2],
        ];
        assert_close(
            LambdaCalculator::find_upper_bound(&scores).unwrap(),
            1.1 * 4.0_f64.ln() / 2.0,
            f64::EPSILON,
        );
    }

    #[test]
    fn lu_inverse_uses_partial_pivoting() {
        let mut matrix = Matrix::new(3);
        matrix.data = vec![0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 2.0, 0.0, 1.0];
        let inverse = LambdaCalculator::invert(&matrix).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let product: f64 = (0..3).map(|k| matrix[(i, k)] * inverse[(k, j)]).sum();
                assert_close(product, if i == j { 1.0 } else { 0.0 }, 1e-12);
            }
        }
    }

    #[test]
    fn singular_matrix_is_not_inverted() {
        let mut matrix = Matrix::new(2);
        matrix.data = vec![1.0, 2.0, 2.0, 4.0];
        assert!(LambdaCalculator::invert(&matrix).is_none());
    }

    #[test]
    fn mt19937_64_matches_standard_engine_sequence() {
        // First values from std::mt19937_64 seeded with 0xC001D00D.
        let mut rng = Mt19937_64::new(0xC001_D00D);
        assert_eq!(rng.next_u64(), 2_086_989_319_759_027_617);
        assert_eq!(rng.next_u64(), 12_351_143_575_849_470_403);
    }
}
