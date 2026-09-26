//! Hierarchy-preserving facade for `diamond/src/lib/alp/sls_alp_sim.cpp`.
//!
//! Simulation state and kernels are shared with the established `stats`
//! implementation; this path exposes the same owned implementation where the
//! original ALP library places it.

pub use crate::stats::sls_alp_sim::{
    alp_sim, calculate_C_S_constant_flag, quick_tests_trials_number, struct_for_lambda_calculation,
    CALCULATE_C_S_CONSTANT_FLAG, QUICK_TESTS_TRIALS_NUMBER,
};

pub type AlpSim = alp_sim;
pub type LambdaCalculation = struct_for_lambda_calculation;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::sls_alp_data::alp_data;

    fn test_data() -> alp_data {
        let scores = vec![vec![-2, 1], vec![1, -2]];
        let frequencies = vec![0.5, 0.5];
        alp_data::new_from_parameters(
            123,
            None,
            11,
            11,
            11,
            1,
            1,
            1,
            2,
            &scores,
            &frequencies,
            &frequencies,
            1.0,
            10.0,
            128.0,
            0.05,
            0.05,
            false,
            -1.0,
            -1.0,
        )
        .unwrap()
    }

    #[test]
    fn empty_state_adapter_initializes_every_public_estimate() {
        let simulation = AlpSim::new_empty(test_data()).unwrap();
        assert_eq!(simulation.d_n_alp_obj, 0);
        assert!(simulation.d_alp_obj.is_empty());
        assert_eq!(simulation.d_mult_number, 0);
        assert_eq!(simulation.m_Lambda, 0.0);
        assert_eq!(simulation.m_K, 0.0);
        assert_eq!(simulation.m_C, 0.0);
        assert!(simulation.m_LambdaSbs.is_empty());
        assert!(simulation.m_AJSbs.is_empty());
        assert!(CALCULATE_C_S_CONSTANT_FLAG);
        assert_eq!(QUICK_TESTS_TRIALS_NUMBER, 100);
    }

    #[test]
    fn source_constructor_runs_the_workflow_and_rejects_missing_random_state() {
        let mut data = test_data();
        data.d_rand_all = None;
        let error = AlpSim::new(data).unwrap_err();
        assert_eq!(error.error_code, 4);
        assert_eq!(error.st, "Unexpected error\n");
    }

    #[test]
    fn numerical_helpers_preserve_rounding_bounds_and_sentinel_errors() {
        assert_eq!(AlpSim::round_double(12.34, 1), 12.3);
        assert_eq!(AlpSim::round_double(12.36, 1), 12.4);
        assert_eq!(AlpSim::relative_error_in_percents(20.0, 1.0), 5.0);
        assert_eq!(AlpSim::relative_error_in_percents(0.0, 1.0), f64::MAX);
        assert_eq!(AlpSim::get_number_of_subsimulations(6).unwrap(), 3);
        assert_eq!(AlpSim::get_number_of_subsimulations(10_000).unwrap(), 20);
        assert!(AlpSim::get_number_of_subsimulations(5).is_err());
        assert_eq!(AlpSim::lambda_exp(1, &[0.0, 2.0]).unwrap(), 2.0);
        assert_eq!(
            AlpSim::lambda_exp(1, &[0.0, -1.0]).unwrap_err().error_code,
            3
        );
    }

    #[test]
    fn fsc_release_clears_and_deallocates_every_owned_buffer() {
        let mut arrays = (0..23).map(|_| vec![1.0, 2.0]).collect::<Vec<_>>();
        let (head, tail) = arrays.split_at_mut(1);
        let [exp_array] = head else { unreachable!() };
        let [a01, a02, a03, a04, a05, a06, a07, a08, a09, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22] =
            tail
        else {
            unreachable!()
        };
        AlpSim::memory_release_for_calculate_fsc(
            exp_array, a01, a02, a03, a04, a05, a06, a07, a08, a09, a10, a11, a12, a13, a14, a15,
            a16, a17, a18, a19, a20, a21, a22,
        );
        assert!(arrays.iter().all(Vec::is_empty));
        assert!(arrays.iter().all(|values| values.capacity() == 0));
    }
}
