#![allow(non_snake_case)]

pub trait AlpFloat:
    Copy
    + PartialOrd
    + PartialEq
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<Output = Self>
    + std::ops::Div<Output = Self>
    + std::ops::Neg<Output = Self>
{
    const ZERO: Self;
    const ONE: Self;
    const LN_2: Self;
    const LN_10: Self;

    fn exp(self) -> Self;
    fn log(self) -> Self;
    fn sqrt(self) -> Self;
    fn abs(self) -> Self;
}

impl AlpFloat for f64 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const LN_2: Self = std::f64::consts::LN_2;
    const LN_10: Self = std::f64::consts::LN_10;

    fn exp(self) -> Self {
        f64::exp(self)
    }

    fn log(self) -> Self {
        f64::ln(self)
    }

    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }

    fn abs(self) -> Self {
        f64::abs(self)
    }
}

impl AlpFloat for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const LN_2: Self = std::f32::consts::LN_2;
    const LN_10: Self = std::f32::consts::LN_10;

    fn exp(self) -> Self {
        f32::exp(self)
    }

    fn log(self) -> Self {
        f32::ln(self)
    }

    fn sqrt(self) -> Self {
        f32::sqrt(self)
    }

    fn abs(self) -> Self {
        f32::abs(self)
    }
}

pub fn bitsToNats<T: AlpFloat>(x_: T) -> T {
    x_ * T::LN_2
}

pub fn natsToBits<T: AlpFloat>(x_: T) -> T {
    x_ / T::LN_2
}

pub fn hartleysToNats<T: AlpFloat>(x_: T) -> T {
    x_ * T::LN_10
}

pub fn natsToHartleys<T: AlpFloat>(x_: T) -> T {
    x_ / T::LN_10
}

pub fn hartleysToBits<T: AlpFloat>(x_: T) -> T {
    hartleysToNats(natsToBits(x_))
}

pub fn bitsToHartleys<T: AlpFloat>(x_: T) -> T {
    bitsToNats(natsToHartleys(x_))
}

pub fn exp2<T: AlpFloat>(x_: T) -> T {
    (T::LN_2 * x_).exp()
}

pub fn log2<T: AlpFloat>(x_: T) -> T {
    x_.log() / T::LN_2
}

pub fn exp10<T: AlpFloat>(x_: T) -> T {
    (T::LN_10 * x_).exp()
}

pub fn log10<T: AlpFloat>(x_: T) -> T {
    x_.log() / T::LN_10
}

pub fn max<T: PartialOrd + Copy>(x_: T, y_: T) -> T {
    if x_ > y_ {
        x_
    } else {
        y_
    }
}

pub fn max3<T: PartialOrd + Copy>(a_: T, b_: T, c_: T) -> T {
    max(a_, max(b_, c_))
}

pub fn min<T: PartialOrd + Copy>(x_: T, y_: T) -> T {
    if x_ < y_ {
        x_
    } else {
        y_
    }
}

pub fn positive<T: AlpFloat>(x_: T) -> T {
    max(x_, T::ZERO)
}

pub fn negative<T: AlpFloat>(x_: T) -> T {
    max(-x_, T::ZERO)
}

pub fn signum<T: AlpFloat>(x_: T) -> T {
    if x_ > T::ZERO {
        T::ONE
    } else if x_ == T::ZERO {
        T::ZERO
    } else {
        -T::ONE
    }
}

pub fn heaviside<T: AlpFloat>(x_: T) -> T {
    if x_ < T::ZERO {
        T::ZERO
    } else {
        T::ONE
    }
}

pub fn psqrt<T: AlpFloat>(x_: T) -> T {
    positive(x_.sqrt())
}

pub fn square<T: AlpFloat>(x_: T) -> T {
    x_ * x_
}

pub fn bound<T: PartialOrd + Copy>(x_: T, xlo_: T, xhi_: T) -> T {
    if x_ < xlo_ {
        return xlo_;
    }
    if x_ < xhi_ {
        return x_;
    }
    xhi_
}

pub fn probability<T: AlpFloat>(x_: T) -> T {
    bound(x_, T::ZERO, T::ONE)
}

pub fn prob<T: AlpFloat>(x_: T) -> T {
    probability(x_)
}

pub fn isProb<T: AlpFloat>(x_: T) -> bool {
    T::ZERO <= x_ && x_ <= T::ONE
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_base_conversions() {
        assert!((bitsToNats(3.0f64) - 3.0 * std::f64::consts::LN_2).abs() < 1e-12);
        assert!((natsToBits(bitsToNats(3.0f64)) - 3.0).abs() < 1e-12);
        assert!((hartleysToNats(2.0f64) - 2.0 * std::f64::consts::LN_10).abs() < 1e-12);
        assert!((natsToHartleys(hartleysToNats(2.0f64)) - 2.0).abs() < 1e-12);
        assert!((hartleysToBits(2.0f64) - hartleysToNats(natsToBits(2.0))).abs() < 1e-12);
        assert!((bitsToHartleys(2.0f64) - bitsToNats(natsToHartleys(2.0))).abs() < 1e-12);
    }

    #[test]
    fn test_exp_log_and_basic_extrema() {
        assert!((exp2(5.0f64) - 32.0).abs() < 1e-12);
        assert!((log2(32.0f64) - 5.0).abs() < 1e-12);
        assert!((exp10(3.0f64) - 1000.0).abs() < 1e-9);
        assert!((log10(1000.0f64) - 3.0).abs() < 1e-12);
        assert_eq!(max(2, 7), 7);
        assert_eq!(max3(2, 7, 4), 7);
        assert_eq!(min(2, 7), 2);
    }

    #[test]
    fn test_sign_bound_probability_and_square() {
        assert_eq!(positive(-2.0f64), 0.0);
        assert_eq!(positive(2.0f64), 2.0);
        assert_eq!(negative(-2.0f64), 2.0);
        assert_eq!(negative(2.0f64), 0.0);
        assert_eq!(signum(-2.0f64), -1.0);
        assert_eq!(signum(0.0f64), 0.0);
        assert_eq!(signum(2.0f64), 1.0);
        assert_eq!(heaviside(-0.1f64), 0.0);
        assert_eq!(heaviside(0.0f64), 1.0);
        assert_eq!(psqrt(9.0f64), 3.0);
        assert_eq!(square(4.0f64), 16.0);
        assert_eq!(bound(5, 1, 3), 3);
        assert_eq!(probability(-0.5f64), 0.0);
        assert_eq!(prob(1.5f64), 1.0);
        assert!(isProb(0.5f64));
        assert!(!isProb(1.5f64));
    }

    #[test]
    fn test_f32_instantiation() {
        assert!((exp2(4.0f32) - 16.0).abs() < 1e-5);
        assert_eq!(probability(2.0f32), 1.0);
    }
}
