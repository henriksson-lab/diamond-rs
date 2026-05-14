#![allow(non_snake_case)]

use crate::stats::alp_random;

pub trait UniformFloat:
    Copy + PartialOrd + std::ops::Add<Output = Self> + std::ops::Sub<Output = Self>
{
    fn from_f64(value: f64) -> Self;
    fn to_f64(self) -> f64;
}

impl UniformFloat for f64 {
    fn from_f64(value: f64) -> Self {
        value
    }

    fn to_f64(self) -> f64 {
        self
    }
}

impl UniformFloat for f32 {
    fn from_f64(value: f64) -> Self {
        value as f32
    }

    fn to_f64(self) -> f64 {
        self as f64
    }
}

pub fn variate<T: UniformFloat>(mut a_: T, mut b_: T) -> T {
    assert!(a_ != b_);

    if b_ < a_ {
        std::mem::swap(&mut a_, &mut b_);
    }

    while alp_random::number() == 0x7fffffff {}

    T::from_f64(
        a_.to_f64() + ((b_ - a_).to_f64() * alp_random::number() as f64 / 0x7fffffff as f64),
    )
}

pub fn standardVariate<T: UniformFloat>() -> T {
    variate(T::from_f64(0.0), T::from_f64(1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_variate_uses_cpp_second_random_number_behavior() {
        let _guard = alp_random::TEST_RANDOM_LOCK.lock().unwrap();
        alp_random::seed(alp_random::SEED);
        let value = variate(2.0f64, 5.0);
        assert!((2.0..5.0).contains(&value));
        assert!((value - 3.541644313159233).abs() < 1e-15);
    }

    #[test]
    fn test_variate_swaps_bounds_and_supports_f32() {
        let _guard = alp_random::TEST_RANDOM_LOCK.lock().unwrap();
        alp_random::seed(alp_random::SEED);
        let value = variate(5.0f32, 2.0);
        assert!((2.0..5.0).contains(&value));
    }

    #[test]
    fn test_standard_variate_range() {
        let _guard = alp_random::TEST_RANDOM_LOCK.lock().unwrap();
        alp_random::seed(alp_random::SEED);
        let value = standardVariate::<f64>();
        assert!((0.0..1.0).contains(&value));
    }

    #[test]
    #[should_panic]
    fn test_variate_rejects_equal_bounds() {
        let _ = variate(1.0f64, 1.0);
    }
}
