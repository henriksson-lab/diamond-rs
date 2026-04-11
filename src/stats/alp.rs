use std::f64::consts::LN_2;

/// Alignment evaluer parameters (matching C++ ALP_set_of_parameters).
#[derive(Debug, Clone)]
pub struct AlpParameters {
    pub lambda: f64,
    pub k: f64,
    pub alpha: f64,
    pub alpha_i: f64,
    pub alpha_j: f64,
    pub beta: f64,
    pub beta_i: f64,
    pub beta_j: f64,
    pub sigma: f64,
    pub tau: f64,
}

/// Alignment E-value calculator with finite-size correction.
///
/// This implements the BLAST-style E-value calculation with
/// edge-effect corrections based on Altschul & Gish (1996).
#[derive(Debug, Clone)]
pub struct AlignmentEvaluer {
    params: AlpParameters,
}

impl AlignmentEvaluer {
    /// Create from pre-computed parameters.
    pub fn new(params: AlpParameters) -> Self {
        AlignmentEvaluer { params }
    }

    /// Create from standard matrix parameters and gap penalties.
    pub fn from_matrix(
        gapped: &super::standard_matrix::Parameters,
        ungapped: &super::standard_matrix::Parameters,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        let g = (gap_open + gap_extend) as f64;
        let b = 2.0 * g * (ungapped.alpha - gapped.alpha);
        let tau = 2.0 * g * (ungapped.alpha_v - gapped.sigma);

        AlignmentEvaluer {
            params: AlpParameters {
                lambda: gapped.lambda,
                k: gapped.k,
                alpha: gapped.alpha,
                alpha_i: gapped.alpha,
                alpha_j: gapped.alpha,
                beta: b,
                beta_i: b,
                beta_j: b,
                sigma: gapped.sigma,
                tau,
            },
        }
    }

    /// Compute E-value for a pairwise alignment.
    pub fn evalue(&self, score: f64, query_len: f64, subject_len: f64) -> f64 {
        self.area(score, query_len, subject_len) * self.evalue_per_area(score)
    }

    /// Compute the effective search area with finite-size correction.
    pub fn area(&self, score: f64, query_len: f64, subject_len: f64) -> f64 {
        let p = &self.params;

        // Effective lengths with edge correction
        let len_correction = p.alpha / p.lambda * (score.ln() + (p.k * query_len * subject_len).ln());
        let eff_query = (query_len - len_correction).max(1.0);
        let eff_subject = (subject_len - len_correction).max(1.0);

        eff_query * eff_subject
    }

    /// Compute log of the effective search area.
    pub fn log_area(&self, score: f64, query_len: f64, subject_len: f64) -> f64 {
        self.area(score, query_len, subject_len).ln()
    }

    /// E-value per unit area.
    pub fn evalue_per_area(&self, score: f64) -> f64 {
        self.params.k * (-self.params.lambda * score).exp()
    }

    /// Bit score from raw score.
    pub fn bit_score(&self, score: f64) -> f64 {
        (self.params.lambda * score - self.params.k.ln()) / LN_2
    }

    /// Access parameters.
    pub fn parameters(&self) -> &AlpParameters {
        &self.params
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_evaluer() -> AlignmentEvaluer {
        use crate::stats::standard_matrix::Parameters;
        // BLOSUM62 with gap_open=11, gap_extend=1
        let gapped = Parameters {
            gap_exist: 11.0, gap_extend: 1.0, reserved: 0.0,
            lambda: 0.267, k: 0.041, h: 0.14,
            alpha: 1.9, beta: -30.0, c: 0.66972,
            alpha_v: 42.6028, sigma: 43.6362,
        };
        let ungapped = Parameters {
            gap_exist: f64::MAX, gap_extend: f64::MAX, reserved: f64::MAX,
            lambda: 0.3176, k: 0.134, h: 0.4012,
            alpha: 0.7916, beta: -3.2, c: 0.623757,
            alpha_v: 4.96466, sigma: 4.96466,
        };
        AlignmentEvaluer::from_matrix(&gapped, &ungapped, 11, 1)
    }

    #[test]
    fn test_evalue_positive_score() {
        let ev = test_evaluer();
        let e = ev.evalue(50.0, 100.0, 1000.0);
        assert!(e > 0.0, "E-value should be positive");
        assert!(e < 1.0, "Score of 50 should give low E-value");
    }

    #[test]
    fn test_evalue_decreases_with_score() {
        let ev = test_evaluer();
        let e1 = ev.evalue(30.0, 100.0, 1000.0);
        let e2 = ev.evalue(50.0, 100.0, 1000.0);
        let e3 = ev.evalue(100.0, 100.0, 1000.0);
        assert!(e1 > e2);
        assert!(e2 > e3);
    }

    #[test]
    fn test_bit_score() {
        let ev = test_evaluer();
        let bs = ev.bit_score(50.0);
        assert!(bs > 0.0);
    }
}
