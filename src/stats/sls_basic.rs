#![allow(non_snake_case)]

use std::time::{SystemTime, UNIX_EPOCH};

use crate::util::log_erfc::erfc;

pub const PI: f64 = 3.1415926535897932384626433832795;
pub const INVSQRTTWO: f64 = std::f64::consts::FRAC_1_SQRT_2;
pub const CONST_VAL: f64 = 0.39894228040143267793994605993438;
pub const QUICK_TESTS_TRIALS_NUMBER: i64 = 100;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error {
    pub st: String,
    pub error_code: i64,
}

impl Error {
    pub fn new(st_: String, error_code_: i64) -> Self {
        Self {
            st: st_,
            error_code: error_code_,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignmentEvaluerParameters {
    pub d_lambda: f64,
    pub d_k: f64,
    pub d_a1: f64,
    pub d_b1: f64,
    pub d_a2: f64,
    pub d_b2: f64,
    pub d_alpha1: f64,
    pub d_beta1: f64,
    pub d_alpha2: f64,
    pub d_beta2: f64,
    pub d_sigma: f64,
    pub d_tau: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignmentEvaluerParametersWithErrors {
    pub d_lambda: f64,
    pub d_lambda_error: f64,
    pub d_k: f64,
    pub d_k_error: f64,
    pub d_a1: f64,
    pub d_a1_error: f64,
    pub d_b1: f64,
    pub d_b1_error: f64,
    pub d_a2: f64,
    pub d_a2_error: f64,
    pub d_b2: f64,
    pub d_b2_error: f64,
    pub d_alpha1: f64,
    pub d_alpha1_error: f64,
    pub d_beta1: f64,
    pub d_beta1_error: f64,
    pub d_alpha2: f64,
    pub d_alpha2_error: f64,
    pub d_beta2: f64,
    pub d_beta2_error: f64,
    pub d_sigma: f64,
    pub d_sigma_error: f64,
    pub d_tau: f64,
    pub d_tau_error: f64,
}

pub fn Tmax<T: PartialOrd + Copy>(i_: T, j_: T) -> T {
    if i_ > j_ {
        i_
    } else {
        j_
    }
}

pub fn Tmin<T: PartialOrd + Copy>(i_: T, j_: T) -> T {
    if i_ < j_ {
        i_
    } else {
        j_
    }
}

pub fn Tmax3<T: PartialOrd + Copy>(x_: T, y_: T, z_: T) -> T {
    Tmax(Tmax(x_, y_), z_)
}

pub fn Tmin3<T: PartialOrd + Copy>(x_: T, y_: T, z_: T) -> T {
    Tmin(Tmin(x_, y_), z_)
}

pub fn Tmax4<T: PartialOrd + Copy>(x_: T, y_: T, z_: T, w_: T) -> T {
    Tmax(Tmax(x_, y_), Tmax(z_, w_))
}

pub fn Tmin4<T: PartialOrd + Copy>(x_: T, y_: T, z_: T, w_: T) -> T {
    Tmin(Tmin(x_, y_), Tmin(z_, w_))
}

pub fn assert_mem<T>(pointer_: *const T) -> Result<(), Error> {
    if pointer_.is_null() {
        Err(Error::new("Memory allocation error\n".to_string(), 41))
    } else {
        Ok(())
    }
}

pub fn round(x_: &f64) -> f64 {
    let x_floor = x_.floor();
    let x_ceil = x_.ceil();
    if (x_ - x_floor).abs() < 0.5 {
        return x_floor;
    }
    x_ceil
}

pub fn get_current_time(seconds_: &mut f64) {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();
    *seconds_ = now.as_secs() as f64 + now.subsec_micros() as f64 / 1_000_000.0;
}

pub fn random_seed_from_time() -> i64 {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();
    let random_factor = now.as_secs() as i64 + now.subsec_micros() as i64 * 10_000_000;
    random_factor.abs()
}

pub fn one_minus_exp_function(y_: f64) -> f64 {
    if y_.abs() > 1e-3 {
        1.0 - y_.exp()
    } else {
        -(y_ * (120.0 + y_ * (60.0 + y_ * (20.0 + y_ * (5.0 + y_)))) / 120.0)
    }
}

pub fn normal_probability(x_: f64) -> f64 {
    0.5 * erfc(-INVSQRTTWO * x_)
}

pub fn normal_probability_eps(x_: f64, mut eps_: f64) -> f64 {
    if x_ == 0.0 {
        return 0.5;
    }

    eps_ = Tmin(1.0, eps_);

    let x_max = 10.0 * eps_ + Tmax(0.0, -2.0 * eps_.ln()).sqrt();

    if x_ >= x_max {
        let x = x_ / 2.0f64.sqrt();
        return 1.0 - 0.5 * (-x * x).exp() / (x * PI.sqrt()) * (1.0 - 1.0 / (2.0 * x * 2.0 * x));
    }

    if x_ <= -x_max {
        let x = x_ / 2.0f64.sqrt();
        return 0.5 * (-x * x).exp() / (-x * PI.sqrt()) * (1.0 - 1.0 / (2.0 * x * 2.0 * x));
    }

    let const_val = 1.0 / (2.0 * PI).sqrt();
    let n = round(&(x_.abs() / eps_)) as i64 + 1;
    let h = x_ / n as f64;

    let mut res = 0.0;
    for i in 0..=n {
        let y = h * i as f64;
        let tmp = (-0.5 * y * y).exp();
        if i == 0 || i == n {
            res += 0.5 * tmp;
        } else {
            res += tmp;
        }
    }

    res *= h;
    0.5 + const_val * res
}

pub fn normal_probability_table(
    a_: f64,
    b_: f64,
    h_: f64,
    n_: i64,
    p_: &[f64],
    x_: f64,
    eps_: f64,
) -> f64 {
    if x_ < a_ || x_ > b_ {
        return normal_probability_eps(x_, eps_);
    }

    let mut x_n = ((x_ - a_) / h_).floor() as i64;
    x_n = Tmin(n_ - 1, x_n);
    let i = x_n as usize;
    p_[i] + (p_[i + 1] - p_[i]) * (x_ - (h_ * x_n as f64 + a_)) / h_
}

pub fn ln_one_minus_val(val_: f64) -> f64 {
    if val_ > 1e-8 {
        return (1.0 - val_).ln();
    }

    -val_ - val_ * val_ / 2.0 - val_ * val_ * val_ / 3.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_min_max_overloads() {
        assert_eq!(Tmax(2, 5), 5);
        assert_eq!(Tmin(2, 5), 2);
        assert_eq!(Tmax3(2, 5, 3), 5);
        assert_eq!(Tmin3(2, 5, 3), 2);
        assert_eq!(Tmax4(2, 5, 3, 7), 7);
        assert_eq!(Tmin4(2, 5, 3, 7), 2);
    }

    #[test]
    fn test_error_and_assert_mem() {
        let error = Error::new("x".to_string(), 7);
        assert_eq!(error.st, "x");
        assert_eq!(error.error_code, 7);
        assert!(assert_mem(&1i32 as *const i32).is_ok());
        let null: *const i32 = std::ptr::null();
        assert_eq!(
            assert_mem(null),
            Err(Error::new("Memory allocation error\n".to_string(), 41))
        );
    }

    #[test]
    fn test_round_matches_cpp_half_behavior() {
        assert_eq!(round(&1.49), 1.0);
        assert_eq!(round(&1.5), 2.0);
        assert_eq!(round(&-1.49), -1.0);
        assert_eq!(round(&-1.5), -1.0);
    }

    #[test]
    fn test_exponential_and_log_stable_forms() {
        assert!((one_minus_exp_function(-0.5) - (1.0 - (-0.5f64).exp())).abs() < 1e-15);
        assert!((one_minus_exp_function(1e-4) - (1.0 - (1e-4f64).exp())).abs() < 1e-16);
        assert!((ln_one_minus_val(0.25) - 0.75f64.ln()).abs() < 1e-15);
        assert!((ln_one_minus_val(1e-10) - (1.0 - 1e-10f64).ln()).abs() < 1e-16);
    }

    #[test]
    fn test_normal_probability_variants() {
        assert!((normal_probability(0.0) - 0.5).abs() < 1e-15);
        assert!((normal_probability(1.0) - 0.8413447460685429).abs() < 1e-15);
        assert!((normal_probability_eps(0.0, 1e-4) - 0.5).abs() < 1e-15);
        assert!((normal_probability_eps(1.0, 1e-4) - normal_probability(1.0)).abs() < 1e-4);

        let p = [0.0, 10.0, 20.0];
        assert_eq!(
            normal_probability_table(0.0, 2.0, 1.0, 2, &p, 0.5, 1e-4),
            5.0
        );
    }
}
