//! x86 kernels for the fixed-size (20 amino acid) composition optimizer.
//!
//! DIAMOND carries AVX2 versions of these operations in
//! `stats/blast/matrix_adjust.cpp`.  The vendored implementation declares its
//! public `MatrixFloat` as `double`, while several of those kernels use
//! `float*`; these kernels retain the production `f64` representation.

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

pub(crate) const ALPHABET: usize = 20;
pub(crate) const CELLS: usize = ALPHABET * ALPHABET;
pub(crate) const CONSTRAINTS: usize = 2 * ALPHABET - 1;

#[inline]
pub(crate) fn available() -> bool {
    std::is_x86_feature_detected!("sse2")
}

#[inline]
pub(crate) unsafe fn multiply_by_a(beta: f64, y: &mut [f64], alpha: f64, x: &[f64]) {
    if std::is_x86_feature_detected!("avx2") {
        multiply_by_a_avx2(beta, y, alpha, x);
    } else {
        multiply_by_a_sse2(beta, y, alpha, x);
    }
}

#[inline]
pub(crate) unsafe fn scaled_symmetric_product(w: &mut [Vec<f64>], diagonal: &[f64]) {
    if std::is_x86_feature_detected!("avx2") {
        scaled_symmetric_product_avx2(w, diagonal);
    } else {
        scaled_symmetric_product_sse2(w, diagonal);
    }
}

#[target_feature(enable = "avx2")]
unsafe fn scaled_symmetric_product_avx2(w: &mut [Vec<f64>], diagonal: &[f64]) {
    for (r, row) in w.iter_mut().enumerate().take(CONSTRAINTS) {
        let mut k = 0;
        while k + 4 <= r + 1 {
            _mm256_storeu_pd(row.as_mut_ptr().add(k), _mm256_setzero_pd());
            k += 4;
        }
        row[k..r + 1].fill(0.0);
    }
    let mut diag_sum = [0.0; ALPHABET];
    for i in 0..ALPHABET {
        let source = &diagonal[i * ALPHABET..(i + 1) * ALPHABET];
        // Keep diagonal reductions in scalar order; these values influence the
        // optimizer's convergence boundary.
        for j in 0..ALPHABET {
            diag_sum[j] += source[j];
        }
        if i > 0 {
            let row = &mut w[ALPHABET - 1 + i];
            for j in (0..ALPHABET).step_by(4) {
                let a = _mm256_loadu_pd(row.as_ptr().add(j));
                let b = _mm256_loadu_pd(source.as_ptr().add(j));
                _mm256_storeu_pd(row.as_mut_ptr().add(j), _mm256_add_pd(a, b));
            }
            for &value in source {
                row[ALPHABET - 1 + i] += value;
            }
        }
    }
    for j in 0..ALPHABET {
        w[j][j] += diag_sum[j];
    }
}

#[target_feature(enable = "sse2")]
unsafe fn scaled_symmetric_product_sse2(w: &mut [Vec<f64>], diagonal: &[f64]) {
    for (r, row) in w.iter_mut().enumerate().take(CONSTRAINTS) {
        let mut k = 0;
        while k + 2 <= r + 1 {
            _mm_storeu_pd(row.as_mut_ptr().add(k), _mm_setzero_pd());
            k += 2;
        }
        row[k..r + 1].fill(0.0);
    }
    let mut diag_sum = [0.0; ALPHABET];
    for i in 0..ALPHABET {
        let source = &diagonal[i * ALPHABET..(i + 1) * ALPHABET];
        for j in 0..ALPHABET {
            diag_sum[j] += source[j];
        }
        if i > 0 {
            let row = &mut w[ALPHABET - 1 + i];
            for j in (0..ALPHABET).step_by(2) {
                let a = _mm_loadu_pd(row.as_ptr().add(j));
                let b = _mm_loadu_pd(source.as_ptr().add(j));
                _mm_storeu_pd(row.as_mut_ptr().add(j), _mm_add_pd(a, b));
            }
            for &value in source {
                row[ALPHABET - 1 + i] += value;
            }
        }
    }
    for j in 0..ALPHABET {
        w[j][j] += diag_sum[j];
    }
}

#[target_feature(enable = "avx2")]
unsafe fn multiply_by_a_avx2(beta: f64, y: &mut [f64], alpha: f64, x: &[f64]) {
    debug_assert!(y.len() >= CONSTRAINTS && x.len() >= CELLS);
    let alpha_v = _mm256_set1_pd(alpha);
    let beta_v = _mm256_set1_pd(beta);
    if beta == 0.0 {
        for k in (0..36).step_by(4) {
            _mm256_storeu_pd(y.as_mut_ptr().add(k), _mm256_setzero_pd());
        }
        y[36..CONSTRAINTS].fill(0.0);
    } else if beta != 1.0 {
        for k in (0..36).step_by(4) {
            let v = _mm256_loadu_pd(y.as_ptr().add(k));
            _mm256_storeu_pd(y.as_mut_ptr().add(k), _mm256_mul_pd(v, beta_v));
        }
        for value in &mut y[36..CONSTRAINTS] {
            *value *= beta;
        }
    }
    for i in 0..ALPHABET {
        let row = x.as_ptr().add(i * ALPHABET);
        for j in (0..ALPHABET).step_by(4) {
            let xv = _mm256_loadu_pd(row.add(j));
            let yv = _mm256_loadu_pd(y.as_ptr().add(j));
            _mm256_storeu_pd(
                y.as_mut_ptr().add(j),
                _mm256_add_pd(yv, _mm256_mul_pd(alpha_v, xv)),
            );
        }
        if i > 0 {
            // Preserve the scalar accumulation order used by the generic NCBI
            // implementation; only the independent column updates are packed.
            for j in 0..ALPHABET {
                y[ALPHABET - 1 + i] += alpha * x[i * ALPHABET + j];
            }
        }
    }
}

#[target_feature(enable = "sse2")]
unsafe fn multiply_by_a_sse2(beta: f64, y: &mut [f64], alpha: f64, x: &[f64]) {
    debug_assert!(y.len() >= CONSTRAINTS && x.len() >= CELLS);
    let alpha_v = _mm_set1_pd(alpha);
    let beta_v = _mm_set1_pd(beta);
    if beta == 0.0 {
        for k in (0..38).step_by(2) {
            _mm_storeu_pd(y.as_mut_ptr().add(k), _mm_setzero_pd());
        }
        y[38] = 0.0;
    } else if beta != 1.0 {
        for k in (0..38).step_by(2) {
            let v = _mm_loadu_pd(y.as_ptr().add(k));
            _mm_storeu_pd(y.as_mut_ptr().add(k), _mm_mul_pd(v, beta_v));
        }
        y[38] *= beta;
    }
    for i in 0..ALPHABET {
        let row = x.as_ptr().add(i * ALPHABET);
        for j in (0..ALPHABET).step_by(2) {
            let xv = _mm_loadu_pd(row.add(j));
            let yv = _mm_loadu_pd(y.as_ptr().add(j));
            _mm_storeu_pd(
                y.as_mut_ptr().add(j),
                _mm_add_pd(yv, _mm_mul_pd(alpha_v, xv)),
            );
        }
        if i > 0 {
            for j in 0..ALPHABET {
                y[ALPHABET - 1 + i] += alpha * x[i * ALPHABET + j];
            }
        }
    }
}

#[inline]
pub(crate) unsafe fn multiply_by_a_transpose(beta: f64, y: &mut [f64], alpha: f64, x: &[f64]) {
    if std::is_x86_feature_detected!("avx2") {
        multiply_by_a_transpose_avx2(beta, y, alpha, x);
    } else {
        multiply_by_a_transpose_sse2(beta, y, alpha, x);
    }
}

#[target_feature(enable = "avx2")]
unsafe fn multiply_by_a_transpose_avx2(beta: f64, y: &mut [f64], alpha: f64, x: &[f64]) {
    debug_assert!(y.len() >= CELLS && x.len() >= CONSTRAINTS);
    let alpha_v = _mm256_set1_pd(alpha);
    let beta_v = _mm256_set1_pd(beta);
    for k in (0..CELLS).step_by(4) {
        let value = if beta == 0.0 {
            _mm256_setzero_pd()
        } else {
            let old = _mm256_loadu_pd(y.as_ptr().add(k));
            if beta == 1.0 {
                old
            } else {
                _mm256_mul_pd(old, beta_v)
            }
        };
        _mm256_storeu_pd(y.as_mut_ptr().add(k), value);
    }
    for i in 0..ALPHABET {
        let row_add = if i == 0 { 0.0 } else { x[ALPHABET - 1 + i] };
        let row_add_v = _mm256_set1_pd(row_add);
        for j in (0..ALPHABET).step_by(4) {
            let xv = _mm256_add_pd(_mm256_loadu_pd(x.as_ptr().add(j)), row_add_v);
            let dst = y.as_mut_ptr().add(i * ALPHABET + j);
            let yv = _mm256_loadu_pd(dst);
            _mm256_storeu_pd(dst, _mm256_add_pd(yv, _mm256_mul_pd(alpha_v, xv)));
        }
    }
}

#[target_feature(enable = "sse2")]
unsafe fn multiply_by_a_transpose_sse2(beta: f64, y: &mut [f64], alpha: f64, x: &[f64]) {
    debug_assert!(y.len() >= CELLS && x.len() >= CONSTRAINTS);
    let alpha_v = _mm_set1_pd(alpha);
    let beta_v = _mm_set1_pd(beta);
    for k in (0..CELLS).step_by(2) {
        let value = if beta == 0.0 {
            _mm_setzero_pd()
        } else {
            let old = _mm_loadu_pd(y.as_ptr().add(k));
            if beta == 1.0 {
                old
            } else {
                _mm_mul_pd(old, beta_v)
            }
        };
        _mm_storeu_pd(y.as_mut_ptr().add(k), value);
    }
    for i in 0..ALPHABET {
        let row_add = if i == 0 { 0.0 } else { x[ALPHABET - 1 + i] };
        let row_add_v = _mm_set1_pd(row_add);
        for j in (0..ALPHABET).step_by(2) {
            let xv = _mm_add_pd(_mm_loadu_pd(x.as_ptr().add(j)), row_add_v);
            let dst = y.as_mut_ptr().add(i * ALPHABET + j);
            let yv = _mm_loadu_pd(dst);
            _mm_storeu_pd(dst, _mm_add_pd(yv, _mm_mul_pd(alpha_v, xv)));
        }
    }
}

#[inline]
pub(crate) unsafe fn dual_residuals(resids: &mut [f64], grads0: &[f64], grads1: &[f64], eta: f64) {
    if std::is_x86_feature_detected!("avx2") {
        dual_residuals_avx2(resids, grads0, grads1, eta);
    } else {
        dual_residuals_sse2(resids, grads0, grads1, eta);
    }
}

#[target_feature(enable = "avx2")]
unsafe fn dual_residuals_avx2(resids: &mut [f64], grads0: &[f64], grads1: &[f64], eta: f64) {
    let eta = _mm256_set1_pd(eta);
    for k in (0..CELLS).step_by(4) {
        let g0 = _mm256_loadu_pd(grads0.as_ptr().add(k));
        let g1 = _mm256_loadu_pd(grads1.as_ptr().add(k));
        _mm256_storeu_pd(
            resids.as_mut_ptr().add(k),
            _mm256_sub_pd(_mm256_mul_pd(eta, g1), g0),
        );
    }
}

#[target_feature(enable = "sse2")]
unsafe fn dual_residuals_sse2(resids: &mut [f64], grads0: &[f64], grads1: &[f64], eta: f64) {
    let eta = _mm_set1_pd(eta);
    for k in (0..CELLS).step_by(2) {
        let g0 = _mm_loadu_pd(grads0.as_ptr().add(k));
        let g1 = _mm_loadu_pd(grads1.as_ptr().add(k));
        _mm_storeu_pd(
            resids.as_mut_ptr().add(k),
            _mm_sub_pd(_mm_mul_pd(eta, g1), g0),
        );
    }
}

#[inline]
pub(crate) unsafe fn compute_scores(scores: &mut [f64], target: &[f64], row: &[f64], col: &[f64]) {
    if std::is_x86_feature_detected!("avx2") {
        compute_scores_avx2(scores, target, row, col);
    } else {
        compute_scores_sse2(scores, target, row, col);
    }
}

#[target_feature(enable = "avx2")]
unsafe fn compute_scores_avx2(scores: &mut [f64], target: &[f64], row: &[f64], col: &[f64]) {
    for i in 0..ALPHABET {
        let r = _mm256_set1_pd(row[i]);
        for j in (0..ALPHABET).step_by(4) {
            let denominator = _mm256_mul_pd(r, _mm256_loadu_pd(col.as_ptr().add(j)));
            let ratio = _mm256_div_pd(
                _mm256_loadu_pd(target.as_ptr().add(i * ALPHABET + j)),
                denominator,
            );
            let mut lanes = [0.0; 4];
            _mm256_storeu_pd(lanes.as_mut_ptr(), ratio);
            for lane in 0..4 {
                scores[i * ALPHABET + j + lane] = lanes[lane].ln();
            }
        }
    }
}

#[target_feature(enable = "sse2")]
unsafe fn compute_scores_sse2(scores: &mut [f64], target: &[f64], row: &[f64], col: &[f64]) {
    for i in 0..ALPHABET {
        let r = _mm_set1_pd(row[i]);
        for j in (0..ALPHABET).step_by(2) {
            let denominator = _mm_mul_pd(r, _mm_loadu_pd(col.as_ptr().add(j)));
            let ratio = _mm_div_pd(
                _mm_loadu_pd(target.as_ptr().add(i * ALPHABET + j)),
                denominator,
            );
            let mut lanes = [0.0; 2];
            _mm_storeu_pd(lanes.as_mut_ptr(), ratio);
            scores[i * ALPHABET + j] = lanes[0].ln();
            scores[i * ALPHABET + j + 1] = lanes[1].ln();
        }
    }
}

#[inline]
pub(crate) unsafe fn evaluate_re(
    values: &mut [f64],
    grads0: &mut [f64],
    grads1: &mut [f64],
    x: &[f64],
    q: &[f64],
    scores: &[f64],
) {
    if std::is_x86_feature_detected!("avx2") {
        evaluate_re_avx2(values, grads0, grads1, x, q, scores);
    } else {
        evaluate_re_sse2(values, grads0, grads1, x, q, scores);
    }
}

#[inline]
pub(crate) unsafe fn factor_lower(a: &mut [Vec<f64>]) {
    if std::is_x86_feature_detected!("avx2") {
        factor_lower_avx2(a);
    } else {
        factor_lower_sse2(a);
    }
}

#[target_feature(enable = "avx2")]
unsafe fn factor_lower_avx2(a: &mut [Vec<f64>]) {
    for i in 0..40 {
        for j in 0..i {
            let mut acc = _mm256_setzero_pd();
            let mut k = 0;
            while k + 4 <= j {
                acc = _mm256_add_pd(
                    acc,
                    _mm256_mul_pd(
                        _mm256_loadu_pd(a[i].as_ptr().add(k)),
                        _mm256_loadu_pd(a[j].as_ptr().add(k)),
                    ),
                );
                k += 4;
            }
            let mut lanes = [0.0; 4];
            _mm256_storeu_pd(lanes.as_mut_ptr(), acc);
            let mut dot = lanes.into_iter().sum::<f64>();
            while k < j {
                dot += a[i][k] * a[j][k];
                k += 1;
            }
            a[i][j] = (a[i][j] - dot) / a[j][j];
        }
        let mut acc = _mm256_setzero_pd();
        let mut k = 0;
        while k + 4 <= i {
            let v = _mm256_loadu_pd(a[i].as_ptr().add(k));
            acc = _mm256_add_pd(acc, _mm256_mul_pd(v, v));
            k += 4;
        }
        let mut lanes = [0.0; 4];
        _mm256_storeu_pd(lanes.as_mut_ptr(), acc);
        let mut sum = lanes.into_iter().sum::<f64>();
        while k < i {
            sum += a[i][k] * a[i][k];
            k += 1;
        }
        a[i][i] = (a[i][i] - sum).sqrt();
    }
}

#[target_feature(enable = "sse2")]
unsafe fn factor_lower_sse2(a: &mut [Vec<f64>]) {
    for i in 0..40 {
        for j in 0..i {
            let mut acc = _mm_setzero_pd();
            let mut k = 0;
            while k + 2 <= j {
                acc = _mm_add_pd(
                    acc,
                    _mm_mul_pd(
                        _mm_loadu_pd(a[i].as_ptr().add(k)),
                        _mm_loadu_pd(a[j].as_ptr().add(k)),
                    ),
                );
                k += 2;
            }
            let mut lanes = [0.0; 2];
            _mm_storeu_pd(lanes.as_mut_ptr(), acc);
            let mut dot = lanes[0] + lanes[1];
            while k < j {
                dot += a[i][k] * a[j][k];
                k += 1;
            }
            a[i][j] = (a[i][j] - dot) / a[j][j];
        }
        let mut acc = _mm_setzero_pd();
        let mut k = 0;
        while k + 2 <= i {
            let v = _mm_loadu_pd(a[i].as_ptr().add(k));
            acc = _mm_add_pd(acc, _mm_mul_pd(v, v));
            k += 2;
        }
        let mut lanes = [0.0; 2];
        _mm_storeu_pd(lanes.as_mut_ptr(), acc);
        let mut sum = lanes[0] + lanes[1];
        while k < i {
            sum += a[i][k] * a[i][k];
            k += 1;
        }
        a[i][i] = (a[i][i] - sum).sqrt();
    }
}

#[inline]
pub(crate) unsafe fn solve_lower(x: &mut [f64], l: &[Vec<f64>]) {
    if std::is_x86_feature_detected!("avx2") {
        solve_lower_avx2(x, l);
    } else {
        solve_lower_sse2(x, l);
    }
}

#[target_feature(enable = "avx2")]
unsafe fn solve_lower_avx2(x: &mut [f64], l: &[Vec<f64>]) {
    for i in 0..40 {
        let mut acc = _mm256_setzero_pd();
        let mut j = 0;
        while j + 4 <= i {
            acc = _mm256_add_pd(
                acc,
                _mm256_mul_pd(
                    _mm256_loadu_pd(l[i].as_ptr().add(j)),
                    _mm256_loadu_pd(x.as_ptr().add(j)),
                ),
            );
            j += 4;
        }
        let mut lanes = [0.0; 4];
        _mm256_storeu_pd(lanes.as_mut_ptr(), acc);
        let mut value = x[i] - lanes.into_iter().sum::<f64>();
        while j < i {
            value -= l[i][j] * x[j];
            j += 1;
        }
        x[i] = value / l[i][i];
    }
    for j in (0..40).rev() {
        x[j] /= l[j][j];
        let xv = _mm256_set1_pd(x[j]);
        let mut i = 0;
        while i + 4 <= j {
            let current = _mm256_loadu_pd(x.as_ptr().add(i));
            let coefficient = _mm256_loadu_pd(l[j].as_ptr().add(i));
            _mm256_storeu_pd(
                x.as_mut_ptr().add(i),
                _mm256_sub_pd(current, _mm256_mul_pd(coefficient, xv)),
            );
            i += 4;
        }
        while i < j {
            x[i] -= l[j][i] * x[j];
            i += 1;
        }
    }
}

#[target_feature(enable = "sse2")]
unsafe fn solve_lower_sse2(x: &mut [f64], l: &[Vec<f64>]) {
    for i in 0..40 {
        let mut acc = _mm_setzero_pd();
        let mut j = 0;
        while j + 2 <= i {
            acc = _mm_add_pd(
                acc,
                _mm_mul_pd(
                    _mm_loadu_pd(l[i].as_ptr().add(j)),
                    _mm_loadu_pd(x.as_ptr().add(j)),
                ),
            );
            j += 2;
        }
        let mut lanes = [0.0; 2];
        _mm_storeu_pd(lanes.as_mut_ptr(), acc);
        let mut value = x[i] - lanes[0] - lanes[1];
        while j < i {
            value -= l[i][j] * x[j];
            j += 1;
        }
        x[i] = value / l[i][i];
    }
    for j in (0..40).rev() {
        x[j] /= l[j][j];
        let xv = _mm_set1_pd(x[j]);
        let mut i = 0;
        while i + 2 <= j {
            let current = _mm_loadu_pd(x.as_ptr().add(i));
            let coefficient = _mm_loadu_pd(l[j].as_ptr().add(i));
            _mm_storeu_pd(
                x.as_mut_ptr().add(i),
                _mm_sub_pd(current, _mm_mul_pd(coefficient, xv)),
            );
            i += 2;
        }
        while i < j {
            x[i] -= l[j][i] * x[j];
            i += 1;
        }
    }
}

#[target_feature(enable = "avx2")]
unsafe fn evaluate_re_avx2(
    values: &mut [f64],
    grads0: &mut [f64],
    grads1: &mut [f64],
    x: &[f64],
    q: &[f64],
    scores: &[f64],
) {
    let one = _mm256_set1_pd(1.0);
    values[0] = 0.0;
    values[1] = 0.0;
    for k in (0..CELLS).step_by(4) {
        let xv = _mm256_loadu_pd(x.as_ptr().add(k));
        let ratio = _mm256_div_pd(xv, _mm256_loadu_pd(q.as_ptr().add(k)));
        let mut t = [0.0; 4];
        _mm256_storeu_pd(t.as_mut_ptr(), ratio);
        for value in &mut t {
            *value = value.ln();
        }
        let tv = _mm256_loadu_pd(t.as_ptr());
        let uv = _mm256_add_pd(tv, _mm256_loadu_pd(scores.as_ptr().add(k)));
        _mm256_storeu_pd(grads0.as_mut_ptr().add(k), _mm256_add_pd(tv, one));
        _mm256_storeu_pd(grads1.as_mut_ptr().add(k), _mm256_add_pd(uv, one));
        let mut u = [0.0; 4];
        _mm256_storeu_pd(u.as_mut_ptr(), uv);
        // Scalar accumulation intentionally preserves the production f64
        // convergence path while lookup/arithmetic remains packed.
        for lane in 0..4 {
            values[0] += x[k + lane] * t[lane];
            values[1] += x[k + lane] * u[lane];
        }
    }
}

#[target_feature(enable = "sse2")]
unsafe fn evaluate_re_sse2(
    values: &mut [f64],
    grads0: &mut [f64],
    grads1: &mut [f64],
    x: &[f64],
    q: &[f64],
    scores: &[f64],
) {
    let one = _mm_set1_pd(1.0);
    values[0] = 0.0;
    values[1] = 0.0;
    for k in (0..CELLS).step_by(2) {
        let ratio = _mm_div_pd(
            _mm_loadu_pd(x.as_ptr().add(k)),
            _mm_loadu_pd(q.as_ptr().add(k)),
        );
        let mut t = [0.0; 2];
        _mm_storeu_pd(t.as_mut_ptr(), ratio);
        t[0] = t[0].ln();
        t[1] = t[1].ln();
        let tv = _mm_loadu_pd(t.as_ptr());
        let uv = _mm_add_pd(tv, _mm_loadu_pd(scores.as_ptr().add(k)));
        _mm_storeu_pd(grads0.as_mut_ptr().add(k), _mm_add_pd(tv, one));
        _mm_storeu_pd(grads1.as_mut_ptr().add(k), _mm_add_pd(uv, one));
        let mut u = [0.0; 2];
        _mm_storeu_pd(u.as_mut_ptr(), uv);
        for lane in 0..2 {
            values[0] += x[k + lane] * t[lane];
            values[1] += x[k + lane] * u[lane];
        }
    }
}
