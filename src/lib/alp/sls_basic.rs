//! Hierarchy-preserving facade for `diamond/src/lib/alp/sls_basic.cpp`.

pub use crate::stats::sls_basic::{
    assert_mem, get_current_time, ln_one_minus_val, normal_probability, normal_probability_eps,
    normal_probability_table, one_minus_exp_function, random_seed_from_time, round, tmax, tmax3,
    tmax4, tmin, tmin3, tmin4, AlignmentEvaluerParameters, AlignmentEvaluerParametersWithErrors,
    Error, CONST_VAL, INVSQRTTWO, PI, QUICK_TESTS_TRIALS_NUMBER,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_edge_cases_and_overloads_are_exposed() {
        assert_eq!(round(&-1.5), -1.0);
        assert_eq!(round(&1.5), 2.0);
        assert_eq!(one_minus_exp_function(0.0), -0.0);
        assert_eq!(normal_probability_eps(0.0, 1e-6), 0.5);

        let table = [0.1, 0.4, 0.9];
        assert_eq!(
            normal_probability_table(-1.0, 1.0, 1.0, 2, &table, 1.0, 1e-6),
            0.9
        );
        assert_eq!(tmax4(1, 4, 3, 2), 4);
        assert_eq!(tmin4(1, 4, 3, 2), 1);
    }

    #[test]
    fn clock_and_time_seed_have_the_source_domain() {
        let mut first = 0.0;
        let mut second = 0.0;
        get_current_time(&mut first);
        get_current_time(&mut second);
        assert!(first > 0.0);
        assert!(second >= first);
        assert!(random_seed_from_time() >= 0);
    }
}
