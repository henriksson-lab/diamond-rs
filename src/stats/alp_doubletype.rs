pub const ROOT_2: f64 = 1.4142135623730951e+000;
pub const ROOT_10: f64 = 3.1622776601683795e+000;

pub const PI: f64 = 3.1415926535897931e+000;
pub const ROOT_PI: f64 = 1.7724538509055159e+000;
pub const ROOT_2_PI: f64 = 2.5066282746310002e+000;

pub const E: f64 = 2.7182818284590451e+000;
pub const LN_2: f64 = 6.9314718055994529e-001;
pub const LN_10: f64 = 2.3025850929940459e+000;
pub const LN_LN_2: f64 = -3.6651292058166435e-001;

pub const EULER: f64 = 0.5772156649015325e+000;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constants_match_cpp_literals() {
        assert_eq!(ROOT_2, 1.4142135623730951e+000);
        assert_eq!(ROOT_10, 3.1622776601683795e+000);
        assert_eq!(PI, 3.1415926535897931e+000);
        assert_eq!(ROOT_PI, 1.7724538509055159e+000);
        assert_eq!(ROOT_2_PI, 2.5066282746310002e+000);
        assert_eq!(E, 2.7182818284590451e+000);
        assert_eq!(LN_2, 6.9314718055994529e-001);
        assert_eq!(LN_10, 2.3025850929940459e+000);
        assert_eq!(LN_LN_2, -3.6651292058166435e-001);
        assert_eq!(EULER, 0.5772156649015325e+000);
    }
}
