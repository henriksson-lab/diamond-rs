pub fn factor_ltriang_pos_def(a: &mut [Vec<f64>], n: usize) {
    for i in 0..n {
        for j in 0..i {
            let mut temp = a[i][j];
            for k in 0..j {
                temp -= a[i][k] * a[j][k];
            }
            a[i][j] = temp / a[j][j];
        }
        let mut temp = a[i][i];
        for k in 0..i {
            temp -= a[i][k] * a[i][k];
        }
        a[i][i] = temp.sqrt();
    }
}

pub fn solve_ltriang_pos_def(x: &mut [f64], n: usize, l: &[Vec<f64>]) {
    for i in 0..n {
        let mut temp = x[i];
        for j in 0..i {
            temp -= l[i][j] * x[j];
        }
        x[i] = temp / l[i][i];
    }
    for j in (0..n).rev() {
        x[j] /= l[j][j];
        for i in 0..j {
            x[i] -= l[j][i] * x[j];
        }
    }
}

pub fn euclidean_norm(v: &[f64], n: usize) -> f64 {
    let mut sum = 1.0;
    let mut scale = 0.0;
    for &vi in v.iter().take(n) {
        if vi != 0.0 {
            let absvi = vi.abs();
            if scale < absvi {
                sum = 1.0 + sum * (scale / absvi) * (scale / absvi);
                scale = absvi;
            } else {
                sum += (absvi / scale) * (absvi / scale);
            }
        }
    }
    scale * sum.sqrt()
}

pub fn add_vectors(y: &mut [f64], n: usize, alpha: f64, x: &[f64]) {
    for i in 0..n {
        y[i] += alpha * x[i];
    }
}

pub fn step_bound(x: &[f64], n: usize, step_x: &[f64], max: f64) -> f64 {
    let mut alpha = max;
    for i in 0..n {
        let alpha_i = -x[i] / step_x[i];
        if alpha_i >= 0.0 && alpha_i < alpha {
            alpha = alpha_i;
        }
    }
    alpha
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factor_and_solve_ltriang_pos_def() {
        let mut a = vec![vec![4.0, 0.0], vec![2.0, 3.0]];
        factor_ltriang_pos_def(&mut a, 2);
        assert!((a[0][0] - 2.0).abs() < 1e-12);
        assert!((a[1][0] - 1.0).abs() < 1e-12);
        assert!((a[1][1] - 2.0_f64.sqrt()).abs() < 1e-12);

        let mut x = vec![6.0, 8.0];
        solve_ltriang_pos_def(&mut x, 2, &a);
        assert!((x[0] - 0.25).abs() < 1e-12);
        assert!((x[1] - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_euclidean_norm() {
        assert_eq!(euclidean_norm(&[0.0, 0.0], 2), 0.0);
        assert!((euclidean_norm(&[3.0, 4.0], 2) - 5.0).abs() < 1e-12);
        assert!(euclidean_norm(&[1.0e200, 1.0e200], 2).is_finite());
    }

    #[test]
    fn test_add_vectors_and_step_bound() {
        let mut y = vec![1.0, 2.0, 3.0];
        add_vectors(&mut y, 3, 2.0, &[4.0, 5.0, 6.0]);
        assert_eq!(y, vec![9.0, 12.0, 15.0]);

        let x = [2.0, 4.0, 6.0];
        let step = [-1.0, -4.0, 3.0];
        assert_eq!(step_bound(&x, 3, &step, 10.0), 1.0);
    }
}
