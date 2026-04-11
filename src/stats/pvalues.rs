use std::f64::consts::PI;

/// Parameters for the ALP finite-size correction.
pub struct AreaParams {
    pub a_i: f64,
    pub b_i: f64,
    pub alpha_i: f64,
    pub beta_i: f64,
    pub a_j: f64,
    pub b_j: f64,
    pub alpha_j: f64,
    pub beta_j: f64,
    pub sigma: f64,
    pub tau: f64,
}

/// Compute the effective search area using the ALP finite-size correction.
///
/// This implements the `get_appr_tail_prob_with_cov_without_errors` function
/// from the ALP library.
pub fn compute_area(p: &AreaParams, score: f64, query_len: f64, subject_len: f64) -> f64 {
    let const_val = 1.0 / (2.0 * PI).sqrt();

    // Query length correction
    let m_li_y = query_len - (p.a_i * score + p.b_i);
    let vi_y = (p.alpha_i * score + p.beta_i).max(0.0);
    let sqrt_vi_y = vi_y.sqrt();

    let m_f = if sqrt_vi_y == 0.0 { 1e100 } else { m_li_y / sqrt_vi_y };

    let p_m_f = normal_cdf(m_f);
    let e_m_f = -const_val * (-0.5 * m_f * m_f).exp();
    let p1 = m_li_y * p_m_f - sqrt_vi_y * e_m_f;

    // Subject length correction
    let n_lj_y = subject_len - (p.a_j * score + p.b_j);
    let vj_y = (p.alpha_j * score + p.beta_j).max(0.0);
    let sqrt_vj_y = vj_y.sqrt();

    let n_f = if sqrt_vj_y == 0.0 { 1e100 } else { n_lj_y / sqrt_vj_y };

    let p_n_f = normal_cdf(n_f);
    let e_n_f = -const_val * (-0.5 * n_f * n_f).exp();
    let p2 = n_lj_y * p_n_f - sqrt_vj_y * e_n_f;

    // Covariance correction
    let c_y = (p.sigma * score + p.tau).max(0.0);

    let area = p1 * p2 + c_y * p_m_f * p_n_f;
    area.max(0.0)
}

/// Compute E-value using the full ALP area computation.
pub fn evalue_with_area(
    lambda: f64,
    k: f64,
    params: &AreaParams,
    score: f64,
    query_len: f64,
    subject_len: f64,
) -> f64 {
    let area = compute_area(params, score, query_len, subject_len);
    area * k * (-lambda * score).exp()
}

/// Standard normal CDF.
fn normal_cdf(x: f64) -> f64 {
    0.5 * erfc(-x / std::f64::consts::SQRT_2)
}

/// Complementary error function (Abramowitz & Stegun approximation).
fn erfc(x: f64) -> f64 {
    if x >= 0.0 {
        erfc_positive(x)
    } else {
        2.0 - erfc_positive(-x)
    }
}

fn erfc_positive(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736
                + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normal_cdf() {
        assert!((normal_cdf(0.0) - 0.5).abs() < 0.001);
        assert!(normal_cdf(3.0) > 0.99);
        assert!(normal_cdf(-3.0) < 0.01);
    }

    #[test]
    fn test_erfc() {
        assert!((erfc(0.0) - 1.0).abs() < 0.001);
        assert!(erfc(3.0) < 0.001);
    }

    #[test]
    fn test_compute_area_basic() {
        let a = 1.9 / 0.267;
        let params = AreaParams {
            a_i: a, b_i: 0.0, alpha_i: 1.9, beta_i: 0.0,
            a_j: a, b_j: 0.0, alpha_j: 1.9, beta_j: 0.0,
            sigma: 43.0, tau: 0.0,
        };
        let area = compute_area(&params, 20.0, 1000.0, 10000.0);
        assert!(area > 0.0, "Area for low score on large seqs: {}", area);
    }

    #[test]
    fn test_evalue_high_score() {
        let params = AreaParams {
            a_i: 1.9, b_i: 0.0, alpha_i: 1.9, beta_i: 0.0,
            a_j: 1.9, b_j: 0.0, alpha_j: 1.9, beta_j: 0.0,
            sigma: 43.0, tau: 0.0,
        };
        let e = evalue_with_area(0.267, 0.041, &params, 2328.0, 426.0, 426.0);
        assert!(e < 1e-200, "E-value for extreme score: {}", e);
    }
}
