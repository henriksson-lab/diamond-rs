//! Hierarchy-preserving facade for ALP `sls_alp_data.cpp`.
//!
//! The owned implementation remains canonical in [`crate::stats::sls_alp_data`].
//! Rust-style aliases coexist here with the source-compatible names retained by
//! established consumers.

pub use crate::stats::sls_alp_data::{
    alp_data, array, array_positive, data_for_lambda_equation, error_for_single_realization,
    importance_sampling, mb_bytes, q_elem, struct_for_randomization,
};

pub type AlpData = alp_data;
pub type Array<T> = array<T>;
pub type ArrayPositive<T> = array_positive<T>;
pub type DataForLambdaEquation = data_for_lambda_equation;
pub type ErrorForSingleRealization = error_for_single_realization;
pub type ImportanceSampling = importance_sampling;
pub type QElem = q_elem;
pub type StructForRandomization = struct_for_randomization;
pub const MB_BYTES: f64 = mb_bytes;

#[cfg(test)]
mod tests {
    use super::*;

    fn small_data() -> AlpData {
        let matrix = vec![vec![-2, 1], vec![1, -2]];
        AlpData::new_from_parameters(
            17,
            None,
            11,
            11,
            11,
            1,
            1,
            1,
            2,
            &matrix,
            &[0.5, 0.5],
            &[0.5, 0.5],
            1.0,
            -1.0,
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
    fn facade_constructor_and_release_preserve_owned_lifecycle() {
        let mut data = small_data();
        assert_eq!(MB_BYTES, 1_048_576.0);
        assert_eq!(data.d_number_of_AA, 2);
        assert_eq!(data.d_smatr, [[-2, 1], [1, -2]]);
        assert_eq!(data.d_RR1_sum, [0.5, 1.0]);
        assert!(data.d_is.is_some());
        data.release_memory();
        assert!(data.d_smatr.is_empty());
        assert!(data.d_RR1.is_empty());
        assert!(data.d_is.is_none());
        assert!(data.d_rand_all.is_none());
        assert_eq!(data.d_memory_size_in_MB, 0.0);
    }

    #[test]
    fn probability_helpers_normalize_and_preserve_zero_weight_plateaus() {
        let mut probabilities = [0.0, 2.0, 0.0, 3.0];
        let (cumulative, elements) = AlpData::calculate_rr_sum(&mut probabilities, 4).unwrap();
        assert_eq!(probabilities, [0.0, 0.4, 0.0, 0.6]);
        assert_eq!(cumulative, [0.0, 0.4, 0.4, 1.0]);
        assert_eq!(elements, [0, 1, 2, 3]);
        assert_eq!(
            AlpData::random_long_element(0.4, 4, &cumulative, &[10, 20, 30, 40]).unwrap(),
            20
        );
        assert!(AlpData::check_rr_sum(0.0, 4, "").is_err());
    }

    #[test]
    fn output_file_validation_includes_instance_symmetry_contract() {
        let data = small_data();
        // The active vendored initialization leaves this flag false even when
        // the two letter distributions are equal; preserve that behavior.
        assert!(!data.d_smatr_symmetric_flag);
        let path = std::env::temp_dir().join(format!(
            "diamond_rs_alp_data_output_{}.txt",
            std::process::id()
        ));
        std::fs::write(&path, "number of realizations with killing nonsymmetric\n").unwrap();
        data.validate_out_file(path.to_str().unwrap()).unwrap();
        std::fs::write(
            &path,
            "number of realizations with killing 0.5* symmetric\n",
        )
        .unwrap();
        let error = data.validate_out_file(path.to_str().unwrap()).unwrap_err();
        assert_eq!(error.error_code, 3);
        assert!(error.st.contains("corresponds to symmetric case"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn header_array_helpers_grow_both_directions_and_copy_independently() {
        let mut positive = ArrayPositive::<i64>::new(None).unwrap();
        positive.increase_elem_by_x(25, 7);
        assert_eq!(positive.d_dim, 29);
        assert_eq!(positive.d_elem[25], 7);

        let mut source = Array::<i64>::new(None);
        source.set_elem(-12, 4);
        source.set_elem(13, 9);
        let mut copy = Array::<i64>::new(None);
        copy.set_elems(&source);
        source.set_elem(-12, 99);
        assert_eq!(copy.d_elem[(-12 - copy.d_ind0) as usize], 4);
        assert_eq!(copy.d_elem[(13 - copy.d_ind0) as usize], 9);
    }
}
