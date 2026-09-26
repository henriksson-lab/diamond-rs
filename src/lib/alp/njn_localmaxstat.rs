//! Hierarchy-preserving facade for `diamond/src/lib/alp/njn_localmaxstat.cpp`.

pub use crate::stats::alp_localmaxstat::LocalMaxStat;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_copy_clear_and_static_time_cover_implementation_state() {
        LocalMaxStat::set_time(0.25);
        assert_eq!(LocalMaxStat::get_time(), 0.25);

        let mut value = LocalMaxStat::default();
        value.copy_full(
            2,
            &[-1, 2],
            &[0.8, 0.2],
            0.4,
            1.0,
            2.0,
            0.1,
            0.9,
            1,
            0.3,
            -0.4,
            1.0,
            0.5,
            0.6,
            3.0,
            true,
        );
        assert_eq!(value.get_r(0.0), 1.0);
        assert_eq!(value.getDimension(), 2);
        assert!(value.getTerminated());

        let copied = value.clone();
        assert_eq!(copied, value);
        value.clear();
        assert!(!value.bool_());
        assert_eq!(value.getDimension(), 0);
        assert_eq!(value.getLambda(), 0.0);
        LocalMaxStat::set_time(0.0);
    }
}
