use std::ops::{Div, Mul, Rem};

pub trait SaturatedAdd: Copy {
    fn saturated_add(x: Self, y: Self) -> Self;
}

impl SaturatedAdd for i8 {
    fn saturated_add(x: Self, y: Self) -> Self {
        saturated_add_i8(x, y)
    }
}

impl SaturatedAdd for i16 {
    fn saturated_add(x: Self, y: Self) -> Self {
        saturated_add_i16(x, y)
    }
}

impl SaturatedAdd for i32 {
    fn saturated_add(x: Self, y: Self) -> Self {
        saturated_add_i32(x, y)
    }
}

pub fn saturated_add<T: SaturatedAdd>(x: T, y: T) -> T {
    T::saturated_add(x, y)
}

pub fn saturated_add_i8(x: i8, y: i8) -> i8 {
    ((x as i32 + y as i32).max(i8::MIN as i32)) as i8
}

pub fn saturated_add_i16(x: i16, y: i16) -> i16 {
    ((x as i32 + y as i32).max(i16::MIN as i32)) as i16
}

pub fn saturated_add_i32(x: i32, y: i32) -> i32 {
    x + y
}

pub fn bit_length(x: u64) -> i32 {
    assert!(x > 0);
    64 - x.leading_zeros() as i32
}

pub fn next_pow2_usize(x: usize) -> usize {
    if x <= 1 {
        return 1;
    }
    if x > (usize::MAX >> 1) {
        return 0;
    }
    1usize << (usize::BITS - (x - 1).leading_zeros())
}

pub fn next_pow2_f64(x: f64) -> usize {
    next_pow2_usize(x.ceil() as usize)
}

pub trait NextPow2 {
    fn next_pow2(self) -> usize;
}

impl NextPow2 for usize {
    fn next_pow2(self) -> usize {
        next_pow2_usize(self)
    }
}

impl NextPow2 for f64 {
    fn next_pow2(self) -> usize {
        next_pow2_f64(self)
    }
}

pub fn next_pow2<T: NextPow2>(x: T) -> usize {
    x.next_pow2()
}

pub fn next_combination<const N: i32>(values: &mut [i32]) -> bool {
    for value in values {
        if *value < N - 1 {
            *value += 1;
            return true;
        } else {
            *value = 0;
        }
    }
    false
}

pub fn power<I>(x: I, p: I) -> I
where
    I: Copy + From<u8> + PartialEq + Div<Output = I> + Rem<Output = I> + Mul<Output = I>,
{
    if p == I::from(0) {
        return I::from(1);
    }
    if p == I::from(1) {
        return x;
    }

    let t = power(x, p / I::from(2));
    if p % I::from(2) == I::from(0) {
        t * t
    } else {
        x * t * t
    }
}

pub fn digits<Int>(mut x: Int, base: Int) -> i32
where
    Int: Copy + From<u8> + PartialOrd + std::ops::DivAssign,
{
    let mut d = 0;
    while x > Int::from(0) {
        x /= base;
        d += 1;
    }
    d
}

pub fn log2_approximate(x: f32) -> f32 {
    let mut i = x.to_bits();
    let mut log_2 = ((i >> 23) & 255) as f32 - 128.0;
    i &= !(255 << 23);
    i += 127 << 23;
    let f = f32::from_bits(i);
    log_2 += (-0.34484843f32 * f + 2.02466578f32) * f - 0.67487759f32;
    log_2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_saturated_add() {
        assert_eq!(saturated_add_i8(10, 20), 30);
        assert_eq!(saturated_add_i8(-120, -20), i8::MIN);
        assert_eq!(saturated_add_i8(120, 20), -116);
        assert_eq!(saturated_add_i16(-32_000, -1000), i16::MIN);
        assert_eq!(saturated_add_i32(2_000_000_000, 1), 2_000_000_001);
        assert_eq!(saturated_add(-120i8, -20i8), i8::MIN);
        assert_eq!(saturated_add(-32_000i16, -1000i16), i16::MIN);
        assert_eq!(saturated_add(5i32, 6i32), 11);
    }

    #[test]
    fn test_bit_length_and_next_pow2() {
        assert_eq!(bit_length(1), 1);
        assert_eq!(bit_length(2), 2);
        assert_eq!(bit_length(255), 8);
        assert_eq!(next_pow2_usize(0), 1);
        assert_eq!(next_pow2_usize(1), 1);
        assert_eq!(next_pow2_usize(2), 2);
        assert_eq!(next_pow2_usize(3), 4);
        assert_eq!(next_pow2_f64(5.1), 8);
        assert_eq!(next_pow2(5usize), 8);
        assert_eq!(next_pow2(5.1f64), 8);
        assert_eq!(next_pow2_usize(usize::MAX), 0);
    }

    #[test]
    fn test_next_combination() {
        let mut v = [0, 0, 0];
        assert!(next_combination::<2>(&mut v));
        assert_eq!(v, [1, 0, 0]);
        assert!(next_combination::<2>(&mut v));
        assert_eq!(v, [0, 1, 0]);
        for _ in 0..5 {
            next_combination::<2>(&mut v);
        }
        assert_eq!(v, [1, 1, 1]);
        assert!(!next_combination::<2>(&mut v));
        assert_eq!(v, [0, 0, 0]);
    }

    #[test]
    fn test_power_and_digits() {
        assert_eq!(power(3i64, 0i64), 1);
        assert_eq!(power(3i64, 1i64), 3);
        assert_eq!(power(3i64, 4i64), 81);
        assert_eq!(digits(0u32, 10u32), 0);
        assert_eq!(digits(9u32, 10u32), 1);
        assert_eq!(digits(100u32, 10u32), 3);
    }

    #[test]
    fn test_log2_approximate() {
        assert!((log2_approximate(2.0) - 1.0049398).abs() < 1e-6);
        assert!((log2_approximate(8.0) - 3.0049398).abs() < 1e-6);
        assert!((log2_approximate(10.0) - 10.0f32.log2()).abs() < 0.02);
    }
}
