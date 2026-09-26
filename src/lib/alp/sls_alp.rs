//! Hierarchy-preserving facade for NCBI's `sls_alp.cpp`.
//!
//! The owned recurrence implementation lives in `stats`. This facade supplies
//! systematic Rust names for upstream methods containing capitalized acronyms.

pub use crate::stats::sls_alp::{alp, small_long, state};
pub use crate::stats::sls_basic::Error;

/// Safe equivalent of the header's generic `alp::swap` helper.
pub fn swap<T>(first: &mut T, second: &mut T) {
    std::mem::swap(first, second);
}

pub trait AlpExt {
    fn random_aa1(&self) -> Result<i64, Error>;
    fn random_aa2(&self) -> Result<i64, Error>;
    fn increment_w_matrix(&mut self);
    fn increment_h_matrix(&mut self);
    fn increment_w_weights(&mut self) -> Result<(), Error>;
    fn increment_h_weights(&mut self) -> Result<(), Error>;
    fn increment_h_weights_without_insertions_after_deletions(&mut self) -> Result<(), Error>;
    fn increment_h_weights_with_insertions_after_deletions(&mut self) -> Result<(), Error>;
    fn increment_h_weights_with_sentinels(&mut self, diff_opt: i64) -> Result<(), Error>;
    fn increment_h_weights_with_sentinels_without_insertions_after_deletions(
        &mut self,
        diff_opt: i64,
    ) -> Result<(), Error>;
    fn increment_h_weights_with_sentinels_with_insertions_after_deletions(
        &mut self,
        diff_opt: i64,
    ) -> Result<(), Error>;
    fn john2_weight_calculation(&mut self, length: i64) -> Result<f64, Error>;
}

impl AlpExt for alp {
    fn random_aa1(&self) -> Result<i64, Error> {
        self.random_AA1()
    }
    fn random_aa2(&self) -> Result<i64, Error> {
        self.random_AA2()
    }
    fn increment_w_matrix(&mut self) {
        self.increment_W_matrix();
    }
    fn increment_h_matrix(&mut self) {
        self.increment_H_matrix();
    }
    fn increment_w_weights(&mut self) -> Result<(), Error> {
        self.increment_W_weights()
    }
    fn increment_h_weights(&mut self) -> Result<(), Error> {
        self.increment_H_weights()
    }
    fn increment_h_weights_without_insertions_after_deletions(&mut self) -> Result<(), Error> {
        self.increment_H_weights_without_insertions_after_deletions()
    }
    fn increment_h_weights_with_insertions_after_deletions(&mut self) -> Result<(), Error> {
        self.increment_H_weights_with_insertions_after_deletions()
    }
    fn increment_h_weights_with_sentinels(&mut self, diff_opt: i64) -> Result<(), Error> {
        self.increment_H_weights_with_sentinels(diff_opt)
    }
    fn increment_h_weights_with_sentinels_without_insertions_after_deletions(
        &mut self,
        diff_opt: i64,
    ) -> Result<(), Error> {
        self.increment_H_weights_with_sentinels_without_insertions_after_deletions(diff_opt)
    }
    fn increment_h_weights_with_sentinels_with_insertions_after_deletions(
        &mut self,
        diff_opt: i64,
    ) -> Result<(), Error> {
        self.increment_H_weights_with_sentinels_with_insertions_after_deletions(diff_opt)
    }
    fn john2_weight_calculation(&mut self, length: i64) -> Result<f64, Error> {
        self.John2_weight_calculation(length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::sls_alp_data::{alp_data, mb_bytes};

    fn test_data(insertions: bool) -> alp_data {
        let scores = vec![vec![-2, 1], vec![1, -2]];
        let frequencies = [0.5, 0.5];
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
            insertions,
            -1.0,
            -1.0,
        )
        .unwrap()
    }

    #[test]
    fn snake_case_recurrences_match_known_ladder_points() {
        for insertions in [false, true] {
            let mut simulation = alp::new(test_data(insertions)).unwrap();
            simulation.increment_sequences();
            simulation.d_seqi[..2].copy_from_slice(&[0, 1]);
            simulation.d_seqj[..2].copy_from_slice(&[1, 0]);
            simulation.d_seqi_len = 1;
            simulation.d_seqj_len = 1;
            simulation.increment_h_weights().unwrap();
            assert_eq!(simulation.d_H_ij_next, 1);
            simulation.d_seqi_len = 2;
            simulation.d_seqj_len = 2;
            simulation.increment_h_weights().unwrap();
            assert_eq!(simulation.d_H_ij_next, 2);
            assert_eq!(simulation.d_M, 2);
        }
    }

    #[test]
    fn partial_release_updates_memory_once() {
        let mut simulation = alp::new(test_data(false)).unwrap();
        simulation.increment_sequences();
        simulation.increment_w_matrix();
        simulation.increment_h_matrix();
        simulation.d_H_matr_len = 1;
        let saved = simulation.save_state().unwrap();
        simulation.d_alp_states.push(Some(saved));

        let bytes = (simulation.d_seqi.len() + simulation.d_seqj.len())
            * std::mem::size_of::<i64>()
            + 12 * simulation.d_W_matr_a_len as usize * std::mem::size_of::<f64>()
            + (16 * simulation.d_H_matr_a_len as usize + simulation.d_H_edge_max.len())
                * std::mem::size_of::<i64>()
            + 8 * std::mem::size_of::<i64>();
        let before = simulation.d_alp_data.d_memory_size_in_MB;
        simulation.partially_release_memory();
        let after = simulation.d_alp_data.d_memory_size_in_MB;
        assert!((before - after - bytes as f64 / mb_bytes).abs() < 1.0e-12);
        assert!(simulation.d_seqi.is_empty());
        assert!(simulation.d_H_edge_max.is_empty());
        simulation.partially_release_memory();
        assert_eq!(simulation.d_alp_data.d_memory_size_in_MB, after);
    }

    #[test]
    fn header_swap_and_degree_edges_are_preserved() {
        let (mut left, mut right) = (3, 7);
        swap(&mut left, &mut right);
        assert_eq!((left, right), (7, 3));
        assert_eq!(alp::degree(0.0, 0.0).unwrap(), 1.0);
        assert_eq!(alp::degree(0.0, 2.0).unwrap(), 0.0);
        assert!(alp::degree(-1.0, 2.0).is_err());
    }
}
