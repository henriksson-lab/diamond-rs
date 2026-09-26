//! Hierarchy-preserving facade for NCBI ALP `njn_localmaxstatutil.cpp`.
//!
//! The numerical implementation remains canonical in
//! [`crate::stats::alp_localmaxstat_util`]. This module gives its public
//! header surface systematic Rust names while retaining the established API.

pub use crate::stats::alp_localmaxstat_util::{
    delta, flatten, lambda, lambda_matrix, mu, n_bury, n_step, r, REL_TOL,
};
pub use crate::stats::alp_localmaxstat_util::{
    descendingLadderEpoch as descending_ladder_epoch,
    descendingLadderEpochRepeat as descending_ladder_epoch_repeat, isLogarithmic as is_logarithmic,
    isProbDist as is_prob_dist, isScoreIncreasing as is_score_increasing, muAssoc as mu_assoc,
    muPowerAssoc as mu_power_assoc, rMin as r_min, thetaMin as theta_min,
    thetaMinusDelta as theta_minus_delta,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangular_flatten_and_validation_cover_header_surface() {
        let scores = vec![vec![-2, 1, 3], vec![1, -2, 3]];
        let probabilities = vec![vec![0.1, 0.2, 0.0], vec![0.3, 0.4, 0.0]];
        let (dimension, score, probability) = flatten(2, &scores, &probabilities, 3);
        assert_eq!(dimension, 2);
        assert_eq!(score, [-2, 1]);
        assert_eq!(probability, [0.5, 0.5]);
        assert!(is_prob_dist(2, &[0.5, 0.5]));
        assert!(!is_prob_dist(2, &[-0.1, 1.1]));
        assert!(is_score_increasing(3, &[-2, 0, 3]));
        assert!(!is_score_increasing(3, &[-2, -2, 3]));
        assert!(is_logarithmic(2, &[-2, 1], &[0.75, 0.25]));
        assert!(!is_logarithmic(2, &[-2, 1], &[0.25, 0.75]));
    }

    #[test]
    fn scalar_random_walk_parameters_satisfy_analytic_identities() {
        let score = [-2, 1];
        let probability = [0.75, 0.25];
        assert_eq!(mu(2, &score, &probability), -1.25);
        let rate = lambda(2, &score, &probability);
        let analytic = ((3.0_f64 + 21.0_f64.sqrt()) / 2.0).ln();
        assert!((rate - analytic).abs() < 2.0e-6);
        assert!((r(2, &score, &probability, rate) - 1.0).abs() < 2.0e-6);

        let minimum_theta = theta_min(2, &score, &probability, rate);
        assert!((minimum_theta - 6.0_f64.ln() / 3.0).abs() < 2.0e-6);
        let minimum_rate = r_min(2, &score, &probability, rate, minimum_theta);
        assert!(minimum_rate > 0.0 && minimum_rate < 1.0);
        assert!(mu_assoc(2, &score, &probability, rate) > 0.0);
        assert!(mu_power_assoc(2, &score, &probability, rate, 2) > 0.0);
        assert_eq!(delta(2, &score), 1);
        assert_eq!(theta_minus_delta(rate, 2, &score), 1.0 - (-rate).exp());

        let matrix = vec![vec![-2, 1], vec![1, -2]];
        let matrix_rate = lambda_matrix(2, &matrix, &[0.5, 0.5]);
        assert!((r(2, &[-2, 1], &[0.5, 0.5], matrix_rate) - 1.0).abs() < 2.0e-6);
    }

    #[test]
    fn independent_threads_do_not_overwrite_root_finding_parameters() {
        let cases = [([-1, 2], [0.8, 0.2]), ([-3, 1], [0.6, 0.4])];
        std::thread::scope(|scope| {
            let handles = cases.map(|(scores, probabilities)| {
                scope.spawn(move || {
                    for _ in 0..250 {
                        let rate = lambda(2, &scores, &probabilities);
                        assert!((r(2, &scores, &probabilities, rate) - 1.0).abs() < 2.0e-6);
                    }
                })
            });
            for handle in handles {
                handle.join().unwrap();
            }
        });
    }
}
