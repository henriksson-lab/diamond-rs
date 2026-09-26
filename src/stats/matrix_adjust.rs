//! Translation boundary for `diamond/src/stats/matrix_adjust.cpp`.
//!
//! The implementations live in [`super::cbs`] because they also form the
//! composition-based-statistics API.  Re-exporting them here preserves the
//! upstream source hierarchy and its systematic C++-to-Rust names.

pub use super::cbs::{
    apply_pseudocounts, blast_composition_matrix_adj, calc_freq_ratios, composition_matrix_adjust,
    high_pair_either_seq, high_pair_frequencies, relative_entropy, scores_std_alphabet,
    test_to_apply_re_adjustment_conditional, true_aa_to_std_target_freqs, CbsThresholds,
    MatrixAdjustRule, FIXED_RE_BLOSUM62,
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::TRUE_AA;

    #[test]
    fn identical_background_compositions_retain_cpp_nan_angle_semantics() {
        let background = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let permissive = CbsThresholds {
            query_match_distance_threshold: -1.0,
            length_ratio_threshold: -1.0,
            angle: -1.0,
        };

        // All three Jensen-Shannon distances are zero.  C++ evaluates the
        // cosine expression as 0/0, so `angle` is NaN and `angle > -1` is
        // false.  The conditional mode must therefore retain RE adjustment.
        assert_eq!(
            test_to_apply_re_adjustment_conditional(
                100,
                100,
                &background,
                &background,
                &background,
                permissive,
            ),
            MatrixAdjustRule::UserSpecifiedRelEntropy
        );
    }

    #[test]
    fn mirrored_helpers_cover_pseudocount_and_standard_alphabet_paths() {
        let background = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let mut observed = [0.0; TRUE_AA as usize];
        apply_pseudocounts(&mut observed, 0, &background);
        assert_eq!(observed, background);

        let joint = vec![1.0; TRUE_AA as usize * TRUE_AA as usize];
        let mut standard = vec![f64::NAN; 26 * 26];
        true_aa_to_std_target_freqs(&mut standard, 26, 26, &joint);
        assert!((standard[0] - 1.0 / 400.0).abs() < f64::EPSILON);
        assert_eq!(standard[20 * 26], 0.0);
        assert_eq!(standard[20], 0.0);
    }
}
