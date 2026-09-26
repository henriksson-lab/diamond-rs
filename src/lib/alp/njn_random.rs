//! Hierarchy-preserving facade for `diamond/src/lib/alp/njn_random.cpp`.

pub use crate::stats::alp_random::SEED;

pub fn seed(value: i64) {
    crate::stats::alp_random::seed(value);
}

pub fn number() -> i64 {
    crate::stats::alp_random::number()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn facade_exposes_header_seed_and_global_generator() {
        let _guard = crate::stats::alp_random::TEST_RANDOM_LOCK.lock().unwrap();
        seed(SEED);
        assert_eq!(number(), 1_253_310_707);
        assert!((0..=0x7fff_ffff).contains(&number()));
    }
}
