//! Hierarchy-preserving facade for `diamond/src/lib/alp/sls_pvalues.cpp`.
//!
//! The canonical owned implementation remains in `stats::pvalues`, where the
//! alignment evaluator already consumes it. This module exposes that same
//! implementation under the original `lib/alp` hierarchy.

pub use crate::stats::pvalues::{pvalues, ALP_set_of_parameters, NAT_CUT_OFF_IN_MAX};

/// Idiomatic aliases; the source-compatible names above remain available.
pub type PValues = pvalues;
pub type AlpSetOfParameters = ALP_set_of_parameters;

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_parameters() -> ALP_set_of_parameters {
        ALP_set_of_parameters {
            lambda: 0.5,
            K: 0.2,
            d_params_flag: true,
            ..Default::default()
        }
    }

    #[test]
    fn constructor_and_parameter_validation_match_the_source_surface() {
        let calculator = PValues::new();
        assert!(!calculator.blast);
        assert_eq!(calculator.eps, 0.0001);
        assert_eq!(calculator.a_normal, -10.0);
        assert_eq!(calculator.b_normal, 10.0);
        assert_eq!(calculator.N_normal, 1000);
        assert_eq!(calculator.h_normal, 0.02);
        assert_eq!(calculator.p_normal.len(), 1001);
        assert!(calculator.p_normal[0] < 1.0e-22);
        assert!((calculator.p_normal[500] - 0.5).abs() < 1.0e-15);
        assert_eq!(calculator.p_normal[1000], 1.0);

        let parameters = simple_parameters();
        assert!(PValues::assert_gumbel_parameters(&parameters));
        assert!(!PValues::assert_gumbel_parameters(
            &ALP_set_of_parameters::default()
        ));
    }

    #[test]
    fn zero_variance_tail_has_the_exact_closed_form_area_and_e_value() {
        let parameters = simple_parameters();
        let mut p = 0.0;
        let mut e = 0.0;
        let mut area = 0.0;
        let mut area_is_one = false;
        PValues::get_appr_tail_prob_with_cov_without_errors(
            &parameters,
            true,
            10.0,
            100.0,
            200.0,
            &mut p,
            &mut e,
            &mut area,
            &mut area_is_one,
            false,
        );

        let expected_e = 20_000.0 * 0.2 * (-5.0_f64).exp();
        assert_eq!(area, 20_000.0);
        assert!((e - expected_e).abs() < 1.0e-14);
        assert!((p - (1.0 - (-expected_e).exp())).abs() < 1.0e-15);
        // The implementation deliberately forces blast=false at entry.
        assert!(!area_is_one);
        let log_area = PValues::new().log_area(&parameters, true, 10.0, 100.0, 200.0);
        assert!((log_area - 20_000.0_f64.ln()).abs() < 1.0e-14);
    }

    #[test]
    fn scalar_and_range_entry_points_agree_and_preserve_error_sentinels() {
        let parameters = simple_parameters();
        let calculator = PValues::new();
        let (mut p, mut p_error, mut e, mut e_error) = (0.0, 0.0, 0.0, 0.0);
        calculator
            .calculate_p_values(
                10.0,
                100.0,
                200.0,
                &parameters,
                &mut p,
                &mut p_error,
                &mut e,
                &mut e_error,
                true,
            )
            .unwrap();
        assert_eq!(p_error, -f64::MAX);
        assert_eq!(e_error, -f64::MAX);

        let (mut ps, mut ps_errors, mut es, mut es_errors) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        calculator
            .calculate_p_values_range(
                10,
                10,
                100.0,
                200.0,
                &parameters,
                &mut ps,
                &mut ps_errors,
                &mut es,
                &mut es_errors,
            )
            .unwrap();
        assert_eq!(ps, vec![p]);
        assert_eq!(es, vec![e]);
        assert_eq!(ps_errors, vec![p_error]);
        assert_eq!(es_errors, vec![e_error]);
    }

    #[test]
    fn stream_roundtrip_and_parse_error_semantics_are_preserved() {
        let mut parameters = simple_parameters();
        parameters.m_LambdaSbs = vec![0.5];
        parameters.m_KSbs = vec![0.2];
        parameters.m_CSbs = vec![0.0];
        parameters.m_AJSbs = vec![0.0];
        parameters.m_AISbs = vec![0.0];
        parameters.m_SigmaSbs = vec![0.0];
        parameters.m_AlphaJSbs = vec![0.0];
        parameters.m_AlphaISbs = vec![0.0];
        let text = parameters.to_string();
        let parsed: AlpSetOfParameters = text.parse().unwrap();
        assert_eq!(parsed.lambda, parameters.lambda);
        assert_eq!(parsed.m_LambdaSbs, vec![0.5]);
        assert!(parsed.d_params_flag);

        let error = "header\n1 2".parse::<AlpSetOfParameters>().unwrap_err();
        assert_eq!(error.error_code, 4);
        assert_eq!(error.st, "Error in the input parameters\n");
    }
}
