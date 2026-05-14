#![allow(non_snake_case)]

use crate::stats::alp_function::AlpFloat;

pub trait AlpInteger:
    Copy
    + PartialOrd
    + PartialEq
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<Output = Self>
    + std::ops::Div<Output = Self>
    + std::ops::Rem<Output = Self>
    + std::ops::Neg<Output = Self>
{
    const ZERO: Self;
    const ONE: Self;
    const TWO: Self;
}

macro_rules! impl_alp_integer {
    ($($t:ty),* $(,)?) => {
        $(
            impl AlpInteger for $t {
                const ZERO: Self = 0;
                const ONE: Self = 1;
                const TWO: Self = 2;
            }
        )*
    };
}

impl_alp_integer!(i8, i16, i32, i64, isize);

pub fn minimum<Int: PartialOrd + Copy>(i: Int, j: Int) -> Int {
    if i < j {
        i
    } else {
        j
    }
}

pub fn maximum<Int: PartialOrd + Copy>(i: Int, j: Int) -> Int {
    if i > j {
        i
    } else {
        j
    }
}

pub fn mod_<Int: AlpInteger>(i: Int, j: Int) -> Int {
    let abs_j = if j >= Int::ZERO { j } else { -j };

    if j == Int::ZERO {
        panic!("Nks_Mod : j == 0");
    }

    if i < Int::ZERO {
        let k = (if i >= Int::ZERO { i } else { -i }) % abs_j;
        if k == Int::ZERO {
            k
        } else {
            abs_j - k
        }
    } else {
        i % abs_j
    }
}

pub fn euclidAlgorithm<Int: AlpInteger>(i_: Int, j_: Int) -> Int {
    let mut abs_i = maximum(
        if i_ >= Int::ZERO { i_ } else { -i_ },
        if j_ >= Int::ZERO { j_ } else { -j_ },
    );
    let mut abs_j = minimum(
        if i_ >= Int::ZERO { i_ } else { -i_ },
        if j_ >= Int::ZERO { j_ } else { -j_ },
    );
    let mut remainder;

    while Int::ZERO < abs_j {
        remainder = mod_(abs_i, abs_j);
        abs_i = abs_j;
        abs_j = remainder;
    }

    abs_i
}

pub fn euclidAlgorithmVector<Int: AlpInteger>(values: &[Int]) -> Int {
    if values.is_empty() {
        return Int::ZERO;
    }

    let mut gcd = values[0];
    for &x in &values[1..] {
        gcd = euclidAlgorithm(gcd, x);
    }
    gcd
}

pub fn minusOnePower<Int: AlpInteger>(i: Int) -> Int {
    if mod_(i, Int::TWO) == Int::ZERO {
        Int::ONE
    } else {
        -Int::ONE
    }
}

pub fn integerPower<Real: AlpFloat>(mut x: Real, n: i64) -> Real {
    if x == Real::ZERO {
        if n < 0 {
            panic!("Int::integerPower <class Real, class Int> : negative exponent of zero");
        } else if n == 0 {
            return Real::ONE;
        } else {
            return Real::ZERO;
        }
    }

    let mut y = Real::ONE;
    let mut i = if n > 0 { n } else { -n };
    while i > 0 {
        if i % 2 != 0 {
            y = y * x;
        }
        x = x * x;
        i /= 2;
    }

    if n < 0 {
        y = Real::ONE / y;
    }

    y
}

pub fn integerPositivePower<Real: AlpFloat>(mut x: Real, n: u64) -> Real {
    if x == Real::ZERO {
        if n == 0 {
            return Real::ONE;
        } else {
            return Real::ZERO;
        }
    }

    let mut y = Real::ONE;
    let mut i = n;
    while i > 0 {
        if i % 2 != 0 {
            y = y * x;
        }
        x = x * x;
        i /= 2;
    }
    y
}

pub fn intPower<Int: AlpInteger>(mut x: Int, n: i64) -> Int {
    assert!(n >= 0);

    if x == Int::ZERO {
        return if n == 0 { Int::ONE } else { Int::ZERO };
    }

    let mut y = Int::ONE;
    let mut i = if n > 0 { n } else { -n };
    while i > 0 {
        if i % 2 != 0 {
            y = y * x;
        }
        x = x * x;
        i /= 2;
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimum_maximum_and_mod() {
        assert_eq!(minimum(3, 5), 3);
        assert_eq!(maximum(3, 5), 5);
        assert_eq!(mod_(7i32, 3), 1);
        assert_eq!(mod_(-7i32, 3), 2);
        assert_eq!(mod_(7i32, -3), 1);
        assert_eq!(mod_(-6i32, 3), 0);
    }

    #[test]
    #[should_panic(expected = "Nks_Mod : j == 0")]
    fn test_mod_zero_panics_like_cpp_abort() {
        let _ = mod_(1i32, 0);
    }

    #[test]
    fn test_euclid_and_minus_one_power() {
        assert_eq!(euclidAlgorithm(84i32, 30), 6);
        assert_eq!(euclidAlgorithm(-84i32, 30), 6);
        assert_eq!(euclidAlgorithmVector::<i32>(&[]), 0);
        assert_eq!(euclidAlgorithmVector(&[18i32, 30, 42]), 6);
        assert_eq!(minusOnePower(4i32), 1);
        assert_eq!(minusOnePower(5i32), -1);
        assert_eq!(minusOnePower(-1i32), -1);
        assert_eq!(minusOnePower(-2i32), 1);
    }

    #[test]
    fn test_power_functions() {
        assert_eq!(integerPower(2.0f64, 10), 1024.0);
        assert_eq!(integerPower(2.0f64, -3), 0.125);
        assert_eq!(integerPower(0.0f64, 0), 1.0);
        assert_eq!(integerPower(0.0f64, 3), 0.0);
        assert_eq!(integerPositivePower(2.0f32, 5), 32.0);
        assert_eq!(integerPositivePower(0.0f32, 0), 1.0);
        assert_eq!(integerPositivePower(0.0f32, 5), 0.0);
        assert_eq!(intPower(3i32, 4), 81);
        assert_eq!(intPower(0i32, 0), 1);
        assert_eq!(intPower(0i32, 4), 0);
    }

    #[test]
    #[should_panic(expected = "negative exponent of zero")]
    fn test_integer_power_zero_negative_exponent_panics() {
        let _ = integerPower(0.0f64, -1);
    }
}
