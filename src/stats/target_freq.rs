use super::linear_algebra::{
    add_vectors, euclidean_norm, factor_ltriang_pos_def, solve_ltriang_pos_def, step_bound,
};

#[derive(Debug, Clone)]
pub struct ReNewtonSystem {
    pub alphsize: usize,
    pub constrain_rel_entropy: bool,
    pub w: Vec<Vec<f64>>,
    pub dinv: Vec<f64>,
    pub grad_re: Vec<f64>,
}

impl ReNewtonSystem {
    pub fn new(alphsize: usize) -> Self {
        ReNewtonSystem {
            alphsize,
            constrain_rel_entropy: true,
            w: vec![vec![0.0; 2 * alphsize]; 2 * alphsize],
            dinv: vec![0.0; alphsize * alphsize],
            grad_re: vec![0.0; alphsize * alphsize],
        }
    }
}

pub fn scaled_symmetric_product_a(w: &mut [Vec<f64>], diagonal: &[f64], alphsize: usize) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if alphsize == super::target_freq_simd::ALPHABET
        && w.len() >= super::target_freq_simd::CONSTRAINTS
        && w.iter()
            .take(super::target_freq_simd::CONSTRAINTS)
            .all(|row| row.len() >= super::target_freq_simd::CONSTRAINTS)
        && diagonal.len() >= super::target_freq_simd::CELLS
        && super::target_freq_simd::available()
    {
        // SAFETY: matrix dimensions and runtime feature support are checked.
        unsafe { super::target_freq_simd::scaled_symmetric_product(w, diagonal) };
        return;
    }
    scaled_symmetric_product_a_scalar(w, diagonal, alphsize);
}

fn scaled_symmetric_product_a_scalar(w: &mut [Vec<f64>], diagonal: &[f64], alphsize: usize) {
    let m = 2 * alphsize - 1;
    for (row_w, row) in w.iter_mut().enumerate().take(m) {
        for value in row.iter_mut().take(row_w + 1) {
            *value = 0.0;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let dd = diagonal[i * alphsize + j];
            w[j][j] += dd;
            if i > 0 {
                w[i + alphsize - 1][j] += dd;
                w[i + alphsize - 1][i + alphsize - 1] += dd;
            }
        }
    }
}

pub fn multiply_by_a(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if alphsize == super::target_freq_simd::ALPHABET
        && y.len() >= super::target_freq_simd::CONSTRAINTS
        && x.len() >= super::target_freq_simd::CELLS
        && super::target_freq_simd::available()
    {
        // SAFETY: `available` establishes SSE2 and the callee performs its own
        // AVX2 runtime check. Slice lengths are checked above.
        unsafe { super::target_freq_simd::multiply_by_a(beta, y, alpha, x) };
        return;
    }
    multiply_by_a_scalar(beta, y, alphsize, alpha, x);
}

fn multiply_by_a_scalar(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for yi in y.iter_mut().take(2 * alphsize - 1) {
            *yi = 0.0;
        }
    } else if beta != 1.0 {
        for yi in y.iter_mut().take(2 * alphsize - 1) {
            *yi *= beta;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            y[j] += alpha * x[i * alphsize + j];
        }
    }
    // Keep the second traversal separate, including its floating-point
    // accumulation order. This is the exact order used by NCBI's
    // `MultiplyByA`; combining each row into a temporary sum changes Newton
    // iterations for ill-conditioned composition matrices.
    for i in 1..alphsize {
        for j in 0..alphsize {
            y[i + alphsize - 1] += alpha * x[i * alphsize + j];
        }
    }
}

pub fn multiply_by_a_transpose(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if alphsize == super::target_freq_simd::ALPHABET
        && y.len() >= super::target_freq_simd::CELLS
        && x.len() >= super::target_freq_simd::CONSTRAINTS
        && super::target_freq_simd::available()
    {
        // SAFETY: see `multiply_by_a` above.
        unsafe { super::target_freq_simd::multiply_by_a_transpose(beta, y, alpha, x) };
        return;
    }
    multiply_by_a_transpose_scalar(beta, y, alphsize, alpha, x);
}

fn multiply_by_a_transpose_scalar(
    beta: f64,
    y: &mut [f64],
    alphsize: usize,
    alpha: f64,
    x: &[f64],
) {
    if beta == 0.0 {
        for yk in y.iter_mut().take(alphsize * alphsize) {
            *yk = 0.0;
        }
    } else if beta != 1.0 {
        for yk in y.iter_mut().take(alphsize * alphsize) {
            *yk *= beta;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            y[k] += alpha * x[j];
            if i > 0 {
                y[k] += alpha * x[i + alphsize - 1];
            }
        }
    }
}

pub fn residuals_linear_constraints(
    r_a: &mut [f64],
    alphsize: usize,
    x: &[f64],
    row_sums: &[f64],
    col_sums: &[f64],
) {
    r_a[..alphsize].copy_from_slice(&col_sums[..alphsize]);
    for i in 1..alphsize {
        r_a[i + alphsize - 1] = row_sums[i];
    }
    multiply_by_a(1.0, r_a, alphsize, -1.0, x);
}

pub fn dual_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    grads: &[Vec<f64>],
    z: &[f64],
    constrain_rel_entropy: bool,
) {
    let n = alphsize * alphsize;
    if constrain_rel_entropy {
        let eta = z[2 * alphsize - 1];
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if alphsize == super::target_freq_simd::ALPHABET && super::target_freq_simd::available() {
            // SAFETY: fixed-size buffers are guaranteed by the optimizer and
            // runtime feature detection is complete.
            unsafe { super::target_freq_simd::dual_residuals(resids_x, &grads[0], &grads[1], eta) };
        } else {
            for i in 0..n {
                resids_x[i] = -grads[0][i] + eta * grads[1][i];
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        for i in 0..n {
            resids_x[i] = -grads[0][i] + eta * grads[1][i];
        }
    } else {
        for i in 0..n {
            resids_x[i] = -grads[0][i];
        }
    }
    multiply_by_a_transpose(1.0, resids_x, alphsize, 1.0, z);
}

pub fn calculate_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    resids_z: &mut [f64],
    values: &[f64],
    grads: &[Vec<f64>],
    row_sums: &[f64],
    col_sums: &[f64],
    x: &[f64],
    z: &[f64],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) -> f64 {
    dual_residuals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    let norm_resids_x = euclidean_norm(resids_x, alphsize * alphsize);

    residuals_linear_constraints(resids_z, alphsize, x, row_sums, col_sums);

    let norm_resids_z = if constrain_rel_entropy {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];
        euclidean_norm(resids_z, 2 * alphsize)
    } else {
        euclidean_norm(resids_z, 2 * alphsize - 1)
    };
    (norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z).sqrt()
}

pub fn factor_re_newton_system(
    newton_system: &mut ReNewtonSystem,
    x: &[f64],
    z: &[f64],
    grads: &[Vec<f64>],
    constrain_rel_entropy: bool,
    workspace: &mut [f64],
) {
    let alphsize = newton_system.alphsize;
    let n = alphsize * alphsize;
    let m = if constrain_rel_entropy {
        2 * alphsize
    } else {
        2 * alphsize - 1
    };

    newton_system.constrain_rel_entropy = constrain_rel_entropy;

    if constrain_rel_entropy {
        let eta = z[m - 1];
        for i in 0..n {
            newton_system.dinv[i] = x[i] / (1.0 - eta);
        }
    } else {
        newton_system.dinv[..n].copy_from_slice(&x[..n]);
    }

    scaled_symmetric_product_a(&mut newton_system.w, &newton_system.dinv, alphsize);

    if constrain_rel_entropy {
        newton_system.grad_re[..n].copy_from_slice(&grads[1][..n]);
        newton_system.w[m - 1][m - 1] = 0.0;
        for i in 0..n {
            workspace[i] = newton_system.dinv[i] * newton_system.grad_re[i];
            newton_system.w[m - 1][m - 1] += newton_system.grad_re[i] * workspace[i];
        }
        multiply_by_a(0.0, &mut newton_system.w[m - 1], alphsize, 1.0, workspace);
    }

    factor_ltriang_pos_def(&mut newton_system.w, m);
}

pub fn solve_re_newton_system(
    x: &mut [f64],
    z: &mut [f64],
    newton_system: &ReNewtonSystem,
    workspace: &mut [f64],
) {
    let alphsize = newton_system.alphsize;
    let n = alphsize * alphsize;
    let m_a = 2 * alphsize - 1;
    let m = if newton_system.constrain_rel_entropy {
        m_a + 1
    } else {
        m_a
    };

    for i in 0..n {
        workspace[i] = x[i] * newton_system.dinv[i];
    }
    multiply_by_a(1.0, z, alphsize, -1.0, workspace);

    if newton_system.constrain_rel_entropy {
        for i in 0..n {
            z[m - 1] -= newton_system.grad_re[i] * workspace[i];
        }
    }

    solve_ltriang_pos_def(z, m, &newton_system.w);

    if newton_system.constrain_rel_entropy {
        for (i, xi) in x.iter_mut().enumerate().take(n) {
            *xi += newton_system.grad_re[i] * z[m - 1];
        }
    }
    multiply_by_a_transpose(1.0, x, alphsize, 1.0, z);

    for (i, xi) in x.iter_mut().enumerate().take(n) {
        *xi *= newton_system.dinv[i];
    }
}

pub fn evaluate_re_functions(
    values: &mut [f64],
    grads: &mut [Vec<f64>],
    alphsize: usize,
    x: &[f64],
    q: &[f64],
    scores: &[f64],
    constrain_rel_entropy: bool,
) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if constrain_rel_entropy
        && alphsize == super::target_freq_simd::ALPHABET
        && super::target_freq_simd::available()
    {
        let (first, rest) = grads.split_at_mut(1);
        // SAFETY: fixed-size optimizer buffers and x86 feature support were
        // checked above; split_at_mut proves the gradient rows disjoint.
        unsafe {
            super::target_freq_simd::evaluate_re(values, &mut first[0], &mut rest[0], x, q, scores)
        };
        return;
    }
    values[0] = 0.0;
    values[1] = 0.0;
    for k in 0..alphsize * alphsize {
        let mut temp = (x[k] / q[k]).ln();
        values[0] += x[k] * temp;
        grads[0][k] = temp + 1.0;

        if constrain_rel_entropy {
            temp += scores[k];
            values[1] += x[k] * temp;
            grads[1][k] = temp + 1.0;
        }
    }
}

pub fn compute_scores_from_probs(
    scores: &mut [f64],
    alphsize: usize,
    target_freqs: &[f64],
    row_freqs: &[f64],
    col_freqs: &[f64],
) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if alphsize == super::target_freq_simd::ALPHABET
        && scores.len() >= super::target_freq_simd::CELLS
        && target_freqs.len() >= super::target_freq_simd::CELLS
        && row_freqs.len() >= super::target_freq_simd::ALPHABET
        && col_freqs.len() >= super::target_freq_simd::ALPHABET
        && super::target_freq_simd::available()
    {
        // SAFETY: lengths and runtime CPU support are checked above.
        unsafe {
            super::target_freq_simd::compute_scores(scores, target_freqs, row_freqs, col_freqs)
        };
        return;
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            scores[k] = (target_freqs[k] / (row_freqs[i] * col_freqs[j])).ln();
        }
    }
}

pub fn blast_optimize_target_frequencies(
    x: &mut [f64],
    alphsize: usize,
    q: &[f64],
    row_sums: &[f64],
    col_sums: &[f64],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
    tol: f64,
    maxits: i32,
) -> (i32, i32) {
    let n = alphsize * alphsize;
    let m_a = 2 * alphsize - 1;
    let m = if constrain_rel_entropy { m_a + 1 } else { m_a };

    let mut values = vec![0.0; 2];
    let mut grads = vec![vec![0.0; n], vec![0.0; n]];
    let mut newton_system = ReNewtonSystem::new(alphsize);
    let mut z = vec![0.0; m_a + 1];
    let mut resids_x = vec![0.0; n];
    let mut resids_z = vec![0.0; m_a + 1];
    let mut old_scores = vec![0.0; n];
    let mut workspace = vec![0.0; n];

    compute_scores_from_probs(&mut old_scores, alphsize, q, row_sums, col_sums);
    x[..n].copy_from_slice(&q[..n]);

    let mut its = 0;
    let mut rnorm;
    loop {
        evaluate_re_functions(
            &mut values,
            &mut grads,
            alphsize,
            x,
            q,
            &old_scores,
            constrain_rel_entropy,
        );
        rnorm = calculate_residuals(
            &mut resids_x,
            alphsize,
            &mut resids_z,
            &values,
            &grads,
            row_sums,
            col_sums,
            x,
            &z,
            constrain_rel_entropy,
            relative_entropy,
        );

        if !(rnorm > tol) {
            break;
        }
        its += 1;
        if its > maxits {
            break;
        }

        factor_re_newton_system(
            &mut newton_system,
            x,
            &z,
            &grads,
            constrain_rel_entropy,
            &mut workspace,
        );
        solve_re_newton_system(&mut resids_x, &mut resids_z, &newton_system, &mut workspace);

        let mut alpha = step_bound(x, n, &resids_x, 1.0 / 0.95);
        alpha *= 0.95;

        add_vectors(x, n, alpha, &resids_x);
        add_vectors(&mut z, m, alpha, &resids_z);
    }

    let mut converged = false;
    if its <= maxits && rnorm <= tol {
        if !constrain_rel_entropy || z[m - 1] < 1.0 {
            converged = true;
        }
    }

    (if converged { 0 } else { 1 }, its)
}

pub fn new_optimize_target_frequencies(
    x: &mut [f64],
    alphsize: usize,
    q: &[f64],
    row_sums: &[f64],
    col_sums: &[f64],
    relative_entropy: f64,
    tol: f64,
    maxits: i32,
) -> (i32, i32) {
    if alphsize != 20
        || x.len() < 400
        || q.len() < 400
        || row_sums.len() < 20
        || col_sums.len() < 20
    {
        return (-1, 0);
    }
    blast_optimize_target_frequencies(
        x,
        alphsize,
        q,
        row_sums,
        col_sums,
        true,
        relative_entropy,
        tol,
        maxits,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiply_by_a_and_transpose() {
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let mut y = vec![9.0; 3];
        multiply_by_a(0.0, &mut y, 2, 1.0, &x);
        assert_eq!(y, vec![4.0, 6.0, 7.0]);

        let mut yt = vec![1.0; 4];
        multiply_by_a_transpose(0.0, &mut yt, 2, 1.0, &y);
        assert_eq!(yt, vec![4.0, 6.0, 11.0, 13.0]);
    }

    #[test]
    fn test_scaled_symmetric_product_a() {
        let diagonal = vec![1.0, 2.0, 3.0, 4.0];
        let mut w = vec![vec![99.0; 3]; 3];
        scaled_symmetric_product_a(&mut w, &diagonal, 2);
        assert_eq!(w[0][0], 4.0);
        assert_eq!(w[1][1], 6.0);
        assert_eq!(w[2][0], 3.0);
        assert_eq!(w[2][1], 4.0);
        assert_eq!(w[2][2], 7.0);
    }

    #[test]
    fn test_residuals_linear_constraints() {
        let x = vec![0.1, 0.2, 0.3, 0.4];
        let row = vec![0.3, 0.7];
        let col = vec![0.4, 0.6];
        let mut r = vec![0.0; 3];
        residuals_linear_constraints(&mut r, 2, &x, &row, &col);
        assert!(r.iter().all(|v| v.abs() < 1e-12));
    }

    #[test]
    fn test_evaluate_and_calculate_residuals() {
        let alphsize = 2;
        let x = vec![0.12, 0.18, 0.28, 0.42];
        let q = x.clone();
        let row = vec![0.3, 0.7];
        let col = vec![0.4, 0.6];
        let mut scores = vec![0.0; 4];
        compute_scores_from_probs(&mut scores, alphsize, &q, &row, &col);
        assert!((scores[0] - 0.0).abs() < 1e-12);

        let mut values = vec![0.0; 2];
        let mut grads = vec![vec![0.0; 4], vec![0.0; 4]];
        evaluate_re_functions(&mut values, &mut grads, alphsize, &x, &q, &scores, true);
        assert!(values[0].abs() < 1e-12);
        assert!(values[1].abs() < 1e-12);

        let z = vec![0.0; 4];
        let mut resids_x = vec![0.0; 4];
        let mut resids_z = vec![0.0; 4];
        let rnorm = calculate_residuals(
            &mut resids_x,
            alphsize,
            &mut resids_z,
            &values,
            &grads,
            &row,
            &col,
            &x,
            &z,
            false,
            0.0,
        );
        assert!(rnorm > 0.0);
    }

    #[test]
    fn test_factor_and_solve_re_newton_system() {
        let alphsize = 2;
        let mut sys = ReNewtonSystem::new(alphsize);
        let x_state = vec![0.12, 0.18, 0.28, 0.42];
        let z_state = vec![0.0; 3];
        let grads = vec![vec![1.0; 4], vec![0.0; 4]];
        let mut workspace = vec![0.0; 4];
        factor_re_newton_system(&mut sys, &x_state, &z_state, &grads, false, &mut workspace);
        assert!(!sys.constrain_rel_entropy);
        assert!(sys.w[0][0] > 0.0);

        let mut rx = vec![0.1, -0.2, 0.3, -0.4];
        let mut rz = vec![0.01, -0.02, 0.03];
        solve_re_newton_system(&mut rx, &mut rz, &sys, &mut workspace);
        assert!(rx.iter().all(|v| v.is_finite()));
        assert!(rz.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_blast_optimize_target_frequencies() {
        let alphsize = 2;
        let row = vec![0.3, 0.7];
        let col = vec![0.4, 0.6];
        let q = vec![0.12, 0.18, 0.28, 0.42];
        let mut x = vec![0.0; 4];
        let (status, iterations) = blast_optimize_target_frequencies(
            &mut x, alphsize, &q, &row, &col, false, 0.0, 1e-10, 50,
        );
        assert_eq!(status, 0);
        assert!(iterations >= 0);
        assert!((x.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        assert!(((x[0] + x[1]) - row[0]).abs() < 1e-9);
        assert!(((x[2] + x[3]) - row[1]).abs() < 1e-9);
        assert!(((x[0] + x[2]) - col[0]).abs() < 1e-9);
        assert!(((x[1] + x[3]) - col[1]).abs() < 1e-9);
    }

    #[test]
    fn test_new_optimize_target_frequencies_entry_point() {
        let alphsize = 20;
        let row = vec![1.0 / alphsize as f64; alphsize];
        let col = vec![1.0 / alphsize as f64; alphsize];
        let q = vec![1.0 / (alphsize * alphsize) as f64; alphsize * alphsize];
        let mut x = vec![0.0; alphsize * alphsize];
        let (status, iterations) =
            new_optimize_target_frequencies(&mut x, alphsize, &q, &row, &col, 0.0, 1e-10, 0);
        assert_ne!(status, -1);
        assert!(iterations >= 0);
        assert!(x.iter().all(|v| v.is_finite()));

        let (bad_status, bad_iterations) =
            new_optimize_target_frequencies(&mut x, 19, &q, &row, &col, 0.0, 1e-10, 50);
        assert_eq!(bad_status, -1);
        assert_eq!(bad_iterations, 0);
    }

    #[test]
    fn randomized_fixed_size_simd_linear_maps_match_scalar() {
        let mut state = 0x7a4d_93e2_61c8_05bfu64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((state >> 11) as f64 / ((1u64 << 53) as f64)) * 4.0 - 2.0
        };
        for &(beta, alpha) in &[(0.0, 1.0), (1.0, -1.0), (0.37, 1.91)] {
            let x = (0..400).map(|_| next()).collect::<Vec<_>>();
            let initial_a = (0..39).map(|_| next()).collect::<Vec<_>>();
            let mut expected_a = initial_a.clone();
            multiply_by_a_scalar(beta, &mut expected_a, 20, alpha, &x);
            let mut actual_a = initial_a;
            multiply_by_a(beta, &mut actual_a, 20, alpha, &x);
            for (i, (&actual, &expected)) in actual_a.iter().zip(&expected_a).enumerate() {
                assert!(
                    (actual - expected).abs() <= 2.0e-13,
                    "A[{i}]: {actual} != {expected}"
                );
            }

            let constraints = (0..39).map(|_| next()).collect::<Vec<_>>();
            let initial_at = (0..400).map(|_| next()).collect::<Vec<_>>();
            let mut expected_at = initial_at.clone();
            multiply_by_a_transpose_scalar(beta, &mut expected_at, 20, alpha, &constraints);
            let mut actual_at = initial_at;
            multiply_by_a_transpose(beta, &mut actual_at, 20, alpha, &constraints);
            for (i, (&actual, &expected)) in actual_at.iter().zip(&expected_at).enumerate() {
                assert!(
                    (actual - expected).abs() <= 2.0e-13,
                    "AT[{i}]: {actual} != {expected}"
                );
            }
        }

        let diagonal = (0..400).map(|_| next().abs()).collect::<Vec<_>>();
        let mut expected_w = vec![vec![9.0; 40]; 40];
        scaled_symmetric_product_a_scalar(&mut expected_w, &diagonal, 20);
        let mut actual_w = vec![vec![9.0; 40]; 40];
        scaled_symmetric_product_a(&mut actual_w, &diagonal, 20);
        for i in 0..39 {
            for j in 0..=i {
                assert_eq!(actual_w[i][j], expected_w[i][j], "W[{i}][{j}]");
            }
        }
    }

    #[test]
    fn randomized_dual_residual_simd_matches_scalar() {
        let mut state = 0xd1b5_4a32_d192_ed03u64;
        let mut next = || {
            state = state
                .wrapping_mul(2862933555777941757)
                .wrapping_add(3037000493);
            ((state >> 11) as f64 / ((1u64 << 53) as f64)) * 2.0 - 1.0
        };
        let grads = vec![
            (0..400).map(|_| next()).collect::<Vec<_>>(),
            (0..400).map(|_| next()).collect::<Vec<_>>(),
        ];
        let z = (0..40).map(|_| next()).collect::<Vec<_>>();
        let mut expected = vec![0.0; 400];
        let eta = z[39];
        for i in 0..400 {
            expected[i] = -grads[0][i] + eta * grads[1][i];
        }
        multiply_by_a_transpose_scalar(1.0, &mut expected, 20, 1.0, &z);
        let mut actual = vec![0.0; 400];
        dual_residuals(&mut actual, 20, &grads, &z, true);
        for (i, (&actual, &expected)) in actual.iter().zip(&expected).enumerate() {
            assert!(
                (actual - expected).abs() <= 2.0e-13,
                "dual[{i}]: {actual} != {expected}"
            );
        }
    }

    #[test]
    fn randomized_score_and_re_simd_match_scalar_formulas() {
        let mut state = 0x49c1_7d32_a85e_60fbu64;
        let mut next_positive = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            1.0e-4 + (state >> 11) as f64 / ((1u64 << 53) as f64)
        };
        let row = (0..20).map(|_| next_positive()).collect::<Vec<_>>();
        let col = (0..20).map(|_| next_positive()).collect::<Vec<_>>();
        let target = (0..400).map(|_| next_positive()).collect::<Vec<_>>();
        let mut scores = vec![0.0; 400];
        compute_scores_from_probs(&mut scores, 20, &target, &row, &col);
        for i in 0..20 {
            for j in 0..20 {
                let k = i * 20 + j;
                assert_eq!(scores[k], (target[k] / (row[i] * col[j])).ln());
            }
        }

        let q = (0..400).map(|_| next_positive()).collect::<Vec<_>>();
        let x = (0..400).map(|_| next_positive()).collect::<Vec<_>>();
        let mut values = [0.0; 2];
        let mut grads = vec![vec![0.0; 400], vec![0.0; 400]];
        evaluate_re_functions(&mut values, &mut grads, 20, &x, &q, &scores, true);
        let mut expected_values = [0.0; 2];
        for k in 0..400 {
            let t = (x[k] / q[k]).ln();
            let u = t + scores[k];
            expected_values[0] += x[k] * t;
            expected_values[1] += x[k] * u;
            assert_eq!(grads[0][k], t + 1.0);
            assert_eq!(grads[1][k], u + 1.0);
        }
        assert_eq!(values, expected_values);
    }
}
