//! Hierarchy-preserving facade for `diamond/src/lib/alp/njn_dynprogproblim.cpp`.

pub use crate::stats::alp_dynprogproblim::DynProgProbLim;

#[cfg(test)]
mod tests {
    use super::*;

    fn add_state(value: i64, state: usize) -> i64 {
        value + state as i64
    }

    #[test]
    fn limits_discard_both_tails_and_preserve_the_window() {
        let probabilities = [0.1, 0.2, 0.4, 0.2, 0.1];
        let mut value = DynProgProbLim::new(
            Some(add_state),
            1,
            Some(&[1.0]),
            -2,
            3,
            Some(&probabilities),
        );
        value.set_limits(-1, 2);
        assert!((value.get_prob_lost() - 0.2).abs() < 1e-15);
        assert_eq!(value.getProb(-1), 0.2);
        assert_eq!(value.getProb(0), 0.4);
        assert_eq!(value.getProb(1), 0.2);
        assert_eq!(value.getProb(-2), 0.0);
        assert_eq!(value.getProb(2), 0.0);
    }

    #[test]
    fn update_accumulates_probability_that_leaves_the_fixed_range() {
        let mut value = DynProgProbLim::new(Some(add_state), 2, Some(&[0.25, 0.75]), 0, 2, None);
        value.set_limits(0, 2);
        value.update();
        value.update();
        assert_eq!(value.getProb(0), 0.0625);
        assert_eq!(value.getProb(1), 0.375);
        assert!((value.get_prob_lost() - 0.5625).abs() < 1e-15);
    }
}
