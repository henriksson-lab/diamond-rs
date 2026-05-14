#![allow(non_snake_case)]

use crate::stats::alp_approx;
use crate::stats::alp_function;

pub const FAILED: f64 = f64::INFINITY;

pub fn newtonRaphson<T, F, DF>(
    y_: f64,
    f_: F,
    df_: DF,
    param_: &T,
    mut p_: f64,
    x_: f64,
    mut q_: f64,
    tol_: f64,
    rtol_: f64,
    itmax_: Option<&mut i64>,
) -> f64
where
    F: Fn(f64, &T) -> f64,
    DF: Fn(f64, &T) -> f64,
{
    assert!(p_ != FAILED && p_ != -FAILED);
    assert!(q_ != FAILED && q_ != -FAILED);

    let fp = f_(p_, param_) - y_;
    let fq = f_(q_, param_) - y_;
    if fp * fq > 0.0 {
        panic!("Root::newtonRaphson : root not bracketed");
    }

    if fp == 0.0 {
        return p_;
    }
    if fq == 0.0 {
        return q_;
    }

    if p_ == q_ {
        panic!("Root::newtonRaphson : p_ == q_");
    }

    let mut x = x_;
    if q_ < p_ {
        std::mem::swap(&mut p_, &mut q_);
    }

    if x_ < p_ || q_ < x_ {
        x = 0.5 * (p_ + q_);
    }

    if fp > 0.0 {
        std::mem::swap(&mut p_, &mut q_);
    }

    let mut dx;
    let mut dxold;
    dxold = p_ - q_;
    dx = dxold;

    let mut iter = 100;
    let itmax = match itmax_ {
        Some(itmax) => itmax,
        None => &mut iter,
    };

    while *itmax > 0 {
        let fx = f_(x, param_) - y_;
        if fx == 0.0 {
            return x;
        } else if fx < 0.0 {
            p_ = x;
        } else {
            q_ = x;
        }
        let dfx = df_(x, param_) - y_;

        if (dfx * (x - p_) - fx) * (dfx * (x - q_) - fx) >= 0.0 || 2.0 * fx.abs() > (dfx * dx).abs()
        {
            dx = dxold;
            dxold = 0.5 * (p_ - q_);
            x = 0.5 * (p_ + q_);

            if dxold.abs() <= tol_ {
                return x;
            }
        } else {
            dx = dxold;
            dxold = fx / dfx;
            x -= dxold;

            if dxold.abs() < tol_ || dxold.abs() < rtol_ * x.abs() {
                if (f_(x - alp_function::signum(dxold) * tol_, param_) - y_) * fx < 0.0 {
                    return x;
                }
            }
        }

        *itmax -= 1;
    }

    FAILED
}

pub fn bisection<T, F>(
    y_: f64,
    f_: F,
    param_: &T,
    mut p_: f64,
    mut q_: f64,
    tol_: f64,
    rtol_: f64,
    itmax_: Option<&mut i64>,
) -> f64
where
    F: Fn(f64, &T) -> f64,
{
    assert!(p_ != FAILED && p_ != -FAILED);
    assert!(q_ != FAILED && q_ != -FAILED);

    let fp = f_(p_, param_) - y_;
    let fq = f_(q_, param_) - y_;
    if fp * fq > 0.0 {
        panic!("Root::bisection : root not bracketed");
    }

    if fp == 0.0 {
        return p_;
    }
    if fq == 0.0 {
        return q_;
    }

    if p_ == q_ {
        panic!("Root::bisection : p_ == q_");
    }

    if fp > 0.0 {
        std::mem::swap(&mut p_, &mut q_);
    }

    let mut iter = 100;
    let itmax = match itmax_ {
        Some(itmax) => itmax,
        None => &mut iter,
    };

    let mut x = 0.5 * (p_ + q_);
    while *itmax > 0 {
        let fx = f_(x, param_) - y_;
        if fx < 0.0 {
            p_ = x;
        } else {
            q_ = x;
        }
        x = 0.5 * (p_ + q_);

        if alp_approx::absRelApprox(p_, x, tol_, rtol_) {
            return x;
        }

        *itmax -= 1;
    }

    FAILED
}

pub fn hunt<T, F>(
    y_: f64,
    f_: F,
    param_: &T,
    mut p_: f64,
    mut q_: f64,
    tol_: f64,
    rtol_: f64,
    itmax_: Option<&mut i64>,
) -> f64
where
    F: Fn(f64, &T) -> f64 + Copy,
{
    assert!(p_ != FAILED && p_ != -FAILED);
    assert!(q_ != FAILED && q_ != -FAILED);

    if p_ == q_ {
        panic!("Root::hunt : p_ == q_");
    }

    let x0 = 0.5 * (p_ + q_);
    let fx0 = f_(x0, param_) - y_;
    if fx0 == 0.0 {
        return x0;
    }

    if q_ < p_ {
        std::mem::swap(&mut p_, &mut q_);
    }

    let mut pts: usize = 2;
    let mut del = 0.5 * (q_ - p_);

    let mut iter = 1000;
    let using_default_itmax = itmax_.is_none();
    let itmax = match itmax_ {
        Some(itmax) => itmax,
        None => &mut iter,
    };

    while tol_ <= del {
        let mut x = p_ + 0.5 * del;
        let mut i = 0usize;
        while i < pts && *itmax > 0 {
            let fx = f_(x, param_) - y_;
            if fx * fx0 < 0.0 {
                return bisection(y_, f_, param_, x, x0, tol_, rtol_, Some(itmax));
            }
            x += del;
            i += 1;
            *itmax -= 1;
        }
        if using_default_itmax && *itmax == 0 {
            return FAILED;
        }

        pts *= 2;
        del *= 0.5;
    }

    FAILED
}

pub fn huntExtreme<T, F>(
    y_: f64,
    f_: F,
    param_: &T,
    mut p_: f64,
    mut q_: f64,
    tol_: f64,
    rtol_: f64,
    itmax_: Option<&mut i64>,
    isLargest_: bool,
) -> f64
where
    F: Fn(f64, &T) -> f64 + Copy,
{
    let mut iter = 1000;
    let itmax = match itmax_ {
        Some(itmax) => itmax,
        None => &mut iter,
    };

    if q_ < p_ {
        std::mem::swap(&mut p_, &mut q_);
    }

    let mut x = hunt(y_, f_, param_, p_, q_, tol_, rtol_, Some(itmax));
    let mut x0 = x;

    if isLargest_ {
        while *itmax > 0 && x != FAILED {
            x0 = x;
            x = hunt(y_, f_, param_, x, q_, tol_, rtol_, Some(itmax));
        }
    } else {
        while *itmax > 0 && x != FAILED {
            x0 = x;
            x = hunt(y_, f_, param_, p_, x, tol_, rtol_, Some(itmax));
        }
    }

    x0
}

pub fn newtonRaphsonNoParam<F, DF>(
    y_: f64,
    f_: F,
    df_: DF,
    p_: f64,
    x_: f64,
    q_: f64,
    tol_: f64,
    rtol_: f64,
    itmax_: Option<&mut i64>,
) -> f64
where
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
{
    newtonRaphson(
        y_,
        |x, _: &()| f_(x),
        |x, _: &()| df_(x),
        &(),
        p_,
        x_,
        q_,
        tol_,
        rtol_,
        itmax_,
    )
}

pub fn bisectionNoParam<F>(
    y_: f64,
    f_: F,
    p_: f64,
    q_: f64,
    tol_: f64,
    rtol_: f64,
    itmax_: Option<&mut i64>,
) -> f64
where
    F: Fn(f64) -> f64,
{
    bisection(y_, |x, _: &()| f_(x), &(), p_, q_, tol_, rtol_, itmax_)
}

pub fn huntNoParam<F>(
    y_: f64,
    f_: F,
    p_: f64,
    q_: f64,
    tol_: f64,
    rtol_: f64,
    itmax_: Option<&mut i64>,
) -> f64
where
    F: Fn(f64) -> f64 + Copy,
{
    hunt(y_, |x, _: &()| f_(x), &(), p_, q_, tol_, rtol_, itmax_)
}

pub fn huntExtremeNoParam<F>(
    y_: f64,
    f_: F,
    p_: f64,
    q_: f64,
    tol_: f64,
    rtol_: f64,
    itmax_: Option<&mut i64>,
    isLargest_: bool,
) -> f64
where
    F: Fn(f64) -> f64 + Copy,
{
    huntExtreme(
        y_,
        |x, _: &()| f_(x),
        &(),
        p_,
        q_,
        tol_,
        rtol_,
        itmax_,
        isLargest_,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bisection_finds_square_root_with_parameter() {
        let param = 2.0;
        let mut itmax = 100;
        let root = bisection(
            0.0,
            |x, p| x * x - *p,
            &param,
            0.0,
            2.0,
            1e-12,
            0.0,
            Some(&mut itmax),
        );
        assert!((root - 2.0f64.sqrt()).abs() < 1e-10);
        assert!(itmax < 100);
    }

    #[test]
    fn test_newton_raphson_finds_square_root() {
        let param = 9.0;
        let root = newtonRaphson(
            0.0,
            |x, p| x * x - *p,
            |x, _| 2.0 * x,
            &param,
            0.0,
            2.0,
            5.0,
            1e-12,
            0.0,
            None,
        );
        assert!((root - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_hunt_finds_root_inside_interval() {
        let root = huntNoParam(0.0, |x| x - 0.375, 0.0, 1.0, 1e-12, 0.0, None);
        assert!((root - 0.375).abs() < 1e-10);
    }

    #[test]
    fn test_hunt_extreme_selects_smallest_or_largest_root() {
        let f = |x: f64| (x - 1.0) * (x - 3.0);
        let smallest = huntExtremeNoParam(0.0, f, 0.0, 4.0, 1e-10, 0.0, None, false);
        let largest = huntExtremeNoParam(0.0, f, 0.0, 4.0, 1e-10, 0.0, None, true);
        assert!((smallest - 1.0).abs() < 1e-8);
        assert!((largest - 3.0).abs() < 1e-8);
    }

    #[test]
    fn test_failed_when_iterations_exhausted() {
        let mut itmax = 0;
        let root = bisectionNoParam(0.0, |x| x - 0.5, 0.0, 1.0, 1e-12, 0.0, Some(&mut itmax));
        assert_eq!(root, FAILED);
    }

    #[test]
    #[should_panic(expected = "Root::bisection : root not bracketed")]
    fn test_bisection_not_bracketed_panics_like_cpp_abort() {
        let _ = bisectionNoParam(0.0, |x| x * x + 1.0, -1.0, 1.0, 1e-12, 0.0, None);
    }
}
