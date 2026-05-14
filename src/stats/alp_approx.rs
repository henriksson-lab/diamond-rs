#![allow(non_snake_case)]

use crate::stats::alp_function::AlpFloat;

pub const FLT_THRESHOLD: f32 = 10.0;
pub const FLT_ROUND: f32 = FLT_THRESHOLD * f32::EPSILON;
pub const DBL_THRESHOLD: f64 = 100.0;
pub const DBL_ROUND: f64 = DBL_THRESHOLD * f64::EPSILON;

pub fn approx<T: AlpFloat>(x_: T, y_: T, eps_: T) -> bool {
    (x_ - y_).abs() <= eps_.abs()
}

pub fn relApprox<T: AlpFloat>(x_: T, y_: T, eps_: T) -> bool {
    approx(x_, y_, eps_ * y_)
}

pub fn absRelApprox<T: AlpFloat>(x_: T, y_: T, tol_: T, rtol_: T) -> bool {
    approx(x_, y_, tol_) || relApprox(x_, y_, rtol_)
}

pub fn eq<T: AlpFloat>(x_: T, y_: T, round_: T) -> bool {
    relApprox(x_, y_, round_)
}

pub fn ge<T: AlpFloat>(x_: T, y_: T, round_: T) -> bool {
    (x_ - y_) >= -(y_.abs()) * round_
}

pub fn gt<T: AlpFloat>(x_: T, y_: T, round_: T) -> bool {
    (x_ - y_) > y_.abs() * round_
}

pub fn ne<T: AlpFloat>(x_: T, y_: T, round_: T) -> bool {
    !eq(x_, y_, round_)
}

pub fn lt<T: AlpFloat>(x_: T, y_: T, round_: T) -> bool {
    !ge(x_, y_, round_)
}

pub fn le<T: AlpFloat>(x_: T, y_: T, round_: T) -> bool {
    !gt(x_, y_, round_)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_approximation_helpers_match_cpp_formulas() {
        assert!(approx(1.0f64, 1.01, 0.02));
        assert!(!approx(1.0f64, 1.03, 0.02));
        assert!(relApprox(100.0f64, 101.0, 0.01));
        assert!(!relApprox(100.0f64, 103.0, 0.01));
        assert!(absRelApprox(100.0f64, 103.0, 0.1, 0.03));
        assert!(absRelApprox(1.0f64, 1.05, 0.1, 0.001));
        assert!(!absRelApprox(1.0f64, 1.2, 0.1, 0.001));
    }

    #[test]
    fn test_rounded_comparisons() {
        assert!(eq(100.0f64, 100.5, 0.01));
        assert!(!eq(100.0f64, 103.0, 0.01));
        assert!(ge(99.5f64, 100.0, 0.01));
        assert!(!ge(98.0f64, 100.0, 0.01));
        assert!(gt(102.0f64, 100.0, 0.01));
        assert!(!gt(100.5f64, 100.0, 0.01));
        assert!(ne(100.0f64, 103.0, 0.01));
        assert!(lt(98.0f64, 100.0, 0.01));
        assert!(le(100.5f64, 100.0, 0.01));
    }

    #[test]
    fn test_constants_and_f32_instantiation() {
        assert_eq!(FLT_THRESHOLD, 10.0);
        assert_eq!(DBL_THRESHOLD, 100.0);
        assert_eq!(FLT_ROUND, 10.0 * f32::EPSILON);
        assert_eq!(DBL_ROUND, 100.0 * f64::EPSILON);
        assert!(approx(1.0f32, 1.001, 0.01));
    }
}
