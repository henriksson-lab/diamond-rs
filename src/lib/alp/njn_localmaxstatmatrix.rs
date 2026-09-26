//! Hierarchy-preserving facade for NCBI's `njn_localmaxstatmatrix.cpp`.
//!
//! The owned implementation lives in `stats` because it is also consumed by
//! the alignment evaluator. Re-exporting the same type here keeps the vendor
//! `lib/alp` path without duplicating its numerical implementation or state.

pub use crate::stats::alp_localmaxstat_matrix::LocalMaxStatMatrix;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::alp_localmaxstat::LocalMaxStat;

    static TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    #[test]
    fn rectangular_matrix_flattens_scores_with_asymmetric_probabilities() {
        let _guard = TEST_LOCK.lock().unwrap();
        let scores = vec![vec![-2, -1, 3], vec![-2, 1, 3]];
        let p = [0.25, 0.75];
        let p2 = [0.7, 0.2, 0.1];
        let matrix = LocalMaxStatMatrix::new(2, Some(&scores), Some(&p), Some(&p2), 3, 0.0);

        assert!(matrix.is_ready());
        assert_eq!(matrix.dim_matrix(), 2);
        assert_eq!(matrix.dim_matrix2(), 3);
        assert_eq!(matrix.score_matrix(), scores);
        assert_eq!(matrix.p(), p);
        assert_eq!(matrix.p2(), p2);
        assert_eq!(matrix.scores(), &[-2, -1, 1, 3]);
        for (actual, expected) in matrix
            .probabilities()
            .iter()
            .zip([0.7_f64, 0.05, 0.15, 0.1])
        {
            assert!((actual - expected).abs() < 1.0e-15);
        }
        assert!(matrix.lambda() > 0.0);
    }

    #[test]
    fn default_second_alphabet_and_copy_from_are_deep_owned() {
        let _guard = TEST_LOCK.lock().unwrap();
        let scores = vec![vec![-1, 2], vec![-1, -1]];
        let p = [0.5, 0.5];
        let original = LocalMaxStatMatrix::new(2, Some(&scores), Some(&p), None, 0, 0.0);
        let mut copied = LocalMaxStatMatrix::default();
        copied.copy_from(&original);

        assert_eq!(copied, original);
        assert_eq!(copied.p2(), p);
        let replacement = vec![vec![-1, 2], vec![-1, -1]];
        copied.copy_matrix(2, &replacement, &[0.6, 0.4], None, 2);
        assert_eq!(original.dim_matrix(), 2);
        assert_eq!(original.score_matrix(), scores);
    }

    #[test]
    fn base_copy_overload_preserves_precomputed_statistics() {
        let mut base = LocalMaxStat::default();
        base.copy_full(
            2,
            &[-1, 2],
            &[0.75, 0.25],
            1.25,
            0.4,
            0.3,
            0.2,
            0.8,
            1,
            -0.1,
            -0.25,
            0.9,
            0.5,
            1.1,
            2.0,
            true,
        );
        let scores = vec![vec![-1, 2]];
        let mut matrix = LocalMaxStatMatrix::default();
        matrix.copy_with_base(base, 1, &scores, &[1.0], None, 1);

        assert_eq!(matrix.lambda(), 1.25);
        assert_eq!(matrix.k(), 0.4);
        assert_eq!(matrix.c(), 0.3);
        assert_eq!(matrix.dimension(), 2);
        assert!(matrix.terminated());
        assert_eq!(matrix.score_matrix(), &[vec![-1]]);
    }
}
