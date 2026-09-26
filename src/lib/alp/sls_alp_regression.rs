//! Hierarchy-preserving facade for ALP `sls_alp_regression.cpp`.
//!
//! The established numerical implementation remains canonical in
//! [`crate::stats::sls_alp_regression`]. This facade supplies systematic names
//! for the source's `LSM` identifiers without breaking existing consumers.

pub use crate::stats::sls_alp_regression::alp_reg;

pub type AlpReg = alp_reg;

pub fn robust_regression_sum_with_cut_lsm(
    min_length: i64,
    values: &[f64],
    errors: &mut [f64],
    cut_left_tail: bool,
    cut_right_tail: bool,
    y: f64,
) -> Result<Option<RegressionFit>, crate::stats::sls_basic::Error> {
    let mut fit = RegressionFit::default();
    let mut calculated = false;
    AlpReg::robust_regression_sum_with_cut_LSM(
        min_length,
        values.len() as i64,
        values,
        errors,
        cut_left_tail,
        cut_right_tail,
        y,
        &mut fit.beta0,
        &mut fit.beta1,
        &mut fit.beta0_error,
        &mut fit.beta1_error,
        &mut fit.k1,
        &mut fit.k2,
        &mut calculated,
    )?;
    Ok(calculated.then_some(fit))
}

pub fn robust_regression_sum_with_cut_lsm_beta1_is_defined(
    min_length: i64,
    values: &[f64],
    errors: &mut [f64],
    cut_left_tail: bool,
    cut_right_tail: bool,
    y: f64,
    beta1: f64,
    beta1_error: f64,
) -> Result<Option<FixedSlopeFit>, crate::stats::sls_basic::Error> {
    let mut fit = FixedSlopeFit {
        beta1,
        beta1_error,
        ..FixedSlopeFit::default()
    };
    let mut calculated = false;
    AlpReg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
        min_length,
        values.len() as i64,
        values,
        errors,
        cut_left_tail,
        cut_right_tail,
        y,
        &mut fit.beta0,
        beta1,
        &mut fit.beta0_error,
        beta1_error,
        &mut fit.k1,
        &mut fit.k2,
        &mut calculated,
    )?;
    Ok(calculated.then_some(fit))
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct RegressionFit {
    pub beta0: f64,
    pub beta1: f64,
    pub beta0_error: f64,
    pub beta1_error: f64,
    pub k1: i64,
    pub k2: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FixedSlopeFit {
    pub beta0: f64,
    pub beta1: f64,
    pub beta0_error: f64,
    pub beta1_error: f64,
    pub k1: i64,
    pub k2: i64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn root_search_preserves_endpoints_partitions_and_multiple_roots() {
        let root = AlpReg::find_single_tetta_general(|x| x.cos(), 0.0, 2.0, 1.0e-12).unwrap();
        assert!((root - std::f64::consts::FRAC_PI_2).abs() < 2.0e-12);

        let mut roots = Vec::new();
        AlpReg::find_tetta_general(
            |x| (x + 1.0) * x * (x - 1.0),
            -1.0,
            1.0,
            4,
            1.0e-12,
            &mut roots,
        )
        .unwrap();
        assert_eq!(roots, [-1.0, 0.0, 1.0]);
        assert!(AlpReg::find_tetta_general(|x| x, 0.0, 1.0, 0, 1.0e-6, &mut roots).is_err());
    }

    #[test]
    fn free_and_fixed_slope_facades_fit_exact_lines() {
        let values = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0];
        let mut errors = [0.5; 6];
        let fit = robust_regression_sum_with_cut_lsm(0, &values, &mut errors, false, false, 0.0)
            .unwrap()
            .unwrap();
        assert!((fit.beta0 - 3.0).abs() < 1.0e-12);
        assert!((fit.beta1 - 2.0).abs() < 1.0e-12);
        assert_eq!((fit.k1, fit.k2), (0, 5));

        let mut errors = [0.5; 6];
        let fixed = robust_regression_sum_with_cut_lsm_beta1_is_defined(
            0,
            &values,
            &mut errors,
            false,
            false,
            0.0,
            2.0,
            0.0,
        )
        .unwrap()
        .unwrap();
        assert!((fixed.beta0 - 3.0).abs() < 1.0e-12);
        assert_eq!((fixed.beta1, fixed.k1, fixed.k2), (2.0, 0, 5));
    }

    #[test]
    fn degeneracy_error_correction_median_and_robust_sum_match_source() {
        let mut errors = [0.0, 0.0, 0.0];
        AlpReg::correction_of_errors(&mut errors, 3).unwrap();
        assert_eq!(errors, [1.0e-50; 3]);

        let mut beta0 = 99.0;
        let mut beta1 = 99.0;
        let mut beta0_error = 99.0;
        let mut beta1_error = 99.0;
        let mut calculated = true;
        let score = AlpReg::function_for_robust_regression_sum_with_cut_LSM(
            &[4.0],
            &[1.0],
            1,
            7,
            0.0,
            &mut beta0,
            &mut beta1,
            &mut beta0_error,
            &mut beta1_error,
            &mut calculated,
        );
        assert_eq!(score, 0.0);
        assert!(!calculated);
        assert_eq!(AlpReg::median(4, &[9.0, 1.0, 3.0, 5.0]), 4.0);

        let mut removed = Vec::new();
        let robust = AlpReg::robust_sum(&[1.0, 2.0, 100.0, 4.0], 4, 1, &mut removed).unwrap();
        assert_eq!(robust, 7.0 / 3.0);
        assert_eq!(removed, [true, true, false, true]);

        // The C implementation accepts a negative removal count: its loop is
        // empty and the denominator is `dim - N_points`.
        let negative = AlpReg::robust_sum(&[2.0, 4.0], 2, -1, &mut removed).unwrap();
        assert_eq!(negative, 2.0);
        assert_eq!(removed, [true, true]);
        assert_eq!(AlpReg::sqrt_for_errors(-0.0), 0.0);
    }
}
