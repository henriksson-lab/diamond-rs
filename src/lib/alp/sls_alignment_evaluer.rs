//! Hierarchy-preserving facade for NCBI's `sls_alignment_evaluer.cpp`.
//!
//! The numerical implementation remains in `stats`, where its ALP consumers
//! already use it. This module exposes the same owned evaluator at the vendor
//! hierarchy and provides systematic Rust names for the legacy C++ methods.

pub use crate::stats::sls_alignment_evaluer::{
    default_importance_sampling_temperature, gapped_computation_parameters_struct, AlignmentEvaluer,
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::sls_basic::{
        AlignmentEvaluerParameters, AlignmentEvaluerParametersWithErrors,
    };

    fn parameters() -> AlignmentEvaluerParameters {
        AlignmentEvaluerParameters {
            d_lambda: 0.267,
            d_k: 0.041,
            d_a1: 1.9,
            d_b1: 4.0,
            d_a2: 2.1,
            d_b2: 5.0,
            d_alpha1: 1.7,
            d_beta1: 6.0,
            d_alpha2: 1.8,
            d_beta2: 7.0,
            d_sigma: 43.0,
            d_tau: 8.0,
        }
    }

    fn parameters_with_errors() -> AlignmentEvaluerParametersWithErrors {
        AlignmentEvaluerParametersWithErrors {
            d_lambda: 0.267,
            d_lambda_error: 0.001,
            d_k: 0.041,
            d_k_error: 0.001,
            d_a1: 1.9,
            d_a1_error: 0.01,
            d_b1: 4.0,
            d_b1_error: 0.01,
            d_a2: 2.1,
            d_a2_error: 0.01,
            d_b2: 5.0,
            d_b2_error: 0.01,
            d_alpha1: 1.7,
            d_alpha1_error: 0.01,
            d_beta1: 6.0,
            d_beta1_error: 0.01,
            d_alpha2: 1.8,
            d_alpha2_error: 0.01,
            d_beta2: 7.0,
            d_beta2_error: 0.01,
            d_sigma: 43.0,
            d_sigma_error: 0.01,
            d_tau: 8.0,
            d_tau_error: 0.01,
        }
    }

    #[test]
    fn hierarchy_facade_exposes_inline_numerical_surface() {
        let mut evaluator = AlignmentEvaluer::new();
        evaluator.init_parameters(&parameters()).unwrap();

        assert!(evaluator.is_good());
        assert_eq!(
            evaluator.evalue_per_area(50.0),
            0.041 * (-0.267_f64 * 50.0).exp()
        );
        assert_eq!(evaluator.bit_score(50.0), evaluator.bitScore(50.0));
        let area = evaluator.area(50.0, 1000.0, 2000.0).unwrap();
        assert!(
            (evaluator.log_area(50.0, 1000.0, 2000.0).unwrap().exp() - area).abs() / area < 1.0e-12
        );
        assert_eq!(
            evaluator.evalue(50.0, 1000.0, 2000.0).unwrap(),
            area * evaluator.evalue_per_area(50.0)
        );
    }

    #[test]
    fn error_sampling_is_reseeded_and_source_deterministic() {
        let input = parameters_with_errors();
        let mut first = AlignmentEvaluer::new();
        let mut second = AlignmentEvaluer::new();
        first.init_parameters_with_errors(&input).unwrap();
        second.init_parameters_with_errors(&input).unwrap();

        assert_eq!(
            first.parameters().m_LambdaSbs,
            second.parameters().m_LambdaSbs
        );
        assert_eq!(first.parameters().m_TauSbs, second.parameters().m_TauSbs);
        assert_eq!(first.parameters().m_LambdaSbs.len(), 20);
        assert!(first.parameters().m_CalcTime >= 0.0);
    }

    #[test]
    fn invalid_frequency_input_clears_previously_valid_state() {
        let mut evaluator = AlignmentEvaluer::new();
        evaluator.init_parameters(&parameters()).unwrap();
        assert!(evaluator.is_good());

        let error = evaluator
            .assert_gapless_input_parameters(2, &[0.0, 0.0], &[0.5, 0.5], "test")
            .unwrap_err();
        assert_eq!(error.error_code, 1);
        assert!(!evaluator.is_good());
    }

    #[test]
    fn stream_and_clear_paths_preserve_source_state() {
        let mut evaluator = AlignmentEvaluer::new();
        evaluator.init_parameters(&parameters()).unwrap();
        let encoded = evaluator.to_string();
        let mut decoded: AlignmentEvaluer = encoded.parse().unwrap();
        assert!(decoded.is_good());

        decoded
            .set_gapped_computation_parameters_simplified(60.0, 1000, 200)
            .unwrap();
        decoded.gapped_computation_parameters_clear();
        let state = decoded.gapped_computation_parameters();
        assert!(!state.d_parameters_flag);
        assert!(state
            .d_first_stage_preliminary_realizations_numbers_ALP
            .is_empty());
        assert_eq!(state.d_max_time_for_quick_tests, -1.0);
        assert_eq!(state.d_max_time_with_computation_parameters, -1.0);
    }
}
