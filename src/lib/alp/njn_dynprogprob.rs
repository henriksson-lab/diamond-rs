//! Hierarchy-preserving facade for NCBI ALP `njn_dynprogprob.cpp`.
//!
//! The established implementation remains canonical in
//! [`crate::stats::alp_dynprogprob`]; this path exposes the audited concrete
//! type alongside the neighboring protocol facade.

pub use crate::stats::alp_dynprogprob::{DynProgProb, ValueFct};

/// Definition of `DynProgProb::ARRAY_CAPACITY` from the implementation unit.
pub const ARRAY_CAPACITY: usize = DynProgProb::ARRAY_CAPACITY;

#[cfg(test)]
mod tests {
    use super::*;

    fn add_input(old: i64, state: usize) -> i64 {
        old + state as i64
    }

    fn jump_right(old: i64, state: usize) -> i64 {
        old + 400 + state as i64
    }

    #[test]
    fn default_and_inline_header_surface_match_contract() {
        let mut value = DynProgProb::new(Some(add_input), 2, Some(&[0.25, 0.75]), 0, 0, None);
        assert_eq!(ARRAY_CAPACITY, 256);
        assert!(value.is_ready());
        assert_eq!(value.array_capacity(), ARRAY_CAPACITY);
        assert_eq!(value.value_begin(), -127);
        assert_eq!((value.value_lower(), value.value_upper()), (0, 1));
        assert_eq!(value.step(), 0);
        assert_eq!(value.input_dimension(), 2);
        assert_eq!(value.input_prob(), &[0.25, 0.75]);
        assert!(value.value_fct().is_some());
        assert_eq!(value.get_prob(-128), 0.0);
        assert_eq!(value.get_prob(0), 1.0);
        assert_eq!(value.get_prob(129), 0.0);

        value.update();
        assert_eq!(value.step(), 1);
        assert_eq!(value.get_prob(0), 0.25);
        assert_eq!(value.get_prob(1), 0.75);
        value.clear_default();
        assert_eq!((value.step(), value.get_prob(0)), (0, 1.0));
    }

    #[test]
    fn explicit_copy_overload_preserves_complete_state_and_independence() {
        let arrays = [vec![0.1, 0.9, 0.0], vec![0.3, 0.2, 0.5]];
        let mut value = DynProgProb::new(None, 0, None, 0, 0, None);
        value.copy_state(
            1,
            &arrays,
            3,
            -1,
            -1,
            2,
            Some(add_input),
            2,
            Some(&[0.4, 0.6]),
        );
        assert_eq!(value.step(), 1);
        assert_eq!(value.arrays(), &arrays);
        assert_eq!(value.get_prob(-1), 0.3);
        assert_eq!(value.get_prob(1), 0.5);

        let mut copy = DynProgProb::new(None, 0, None, 0, 0, None);
        copy.copy(&value);
        value.clear_default();
        assert_eq!(copy.get_prob(-1), 0.3);
        assert_eq!(copy.input_prob(), &[0.4, 0.6]);
    }

    #[test]
    fn reserve_right_and_input_resize_preserve_probability_mass() {
        let mut value = DynProgProb::new(Some(jump_right), 2, Some(&[0.5, 0.5]), 0, 0, None);
        value.update();
        assert!(value.array_capacity() >= 1024);
        assert_eq!(value.get_prob(400), 0.5);
        assert_eq!(value.get_prob(401), 0.5);

        value.set_input(3, Some(&[0.2, 0.3, 0.5]));
        assert_eq!(value.input_dimension(), 3);
        assert_eq!(value.input_prob(), &[0.2, 0.3, 0.5]);
        value.set_input(0, None);
        assert!(!value.is_ready());
        assert!(value.input_prob().is_empty());
    }

    #[test]
    fn explicit_probability_range_ignores_trailing_slice_data_like_c_pointer() {
        let value = DynProgProb::new(
            Some(add_input),
            1,
            Some(&[1.0]),
            -1,
            1,
            Some(&[0.25, 0.75, -99.0]),
        );
        assert_eq!(value.array_capacity(), 2);
        assert_eq!(value.get_prob(-1), 0.25);
        assert_eq!(value.get_prob(0), 0.75);
    }
}
