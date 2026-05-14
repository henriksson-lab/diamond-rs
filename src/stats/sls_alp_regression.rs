#![allow(non_snake_case, non_camel_case_types)]

use crate::stats::sls_basic::{round, Error, Tmax};

pub struct alp_reg;

impl alp_reg {
    pub fn find_tetta_general<F>(
        func_: F,
        a_: f64,
        b_: f64,
        n_partition_: i64,
        eps_: f64,
        res_: &mut Vec<f64>,
    ) -> Result<(), Error>
    where
        F: Fn(f64) -> f64,
    {
        res_.clear();
        let mut intervals = Vec::new();

        if n_partition_ <= 0 {
            return Err(Error::new(
                "Error in alp_reg::find_tetta_general\n".to_string(),
                4,
            ));
        }

        let h = (b_ - a_) / n_partition_ as f64;
        let mut x2 = 0.0;
        for i in 0..n_partition_ {
            let x1 = if i == 0 {
                let x1 = func_(a_ + i as f64 * h);
                if x1.abs() < eps_ {
                    res_.push(a_ + i as f64 * h);
                }
                x1
            } else {
                x2
            };

            x2 = func_(a_ + (i + 1) as f64 * h);
            if x2.abs() < eps_ {
                res_.push(a_ + (i + 1) as f64 * h);
            }

            if x1 * x2 < 0.0 && x1.abs() >= eps_ && x2.abs() >= eps_ {
                intervals.push(i);
            }
        }

        for interval in intervals {
            let sol = Self::find_single_tetta_general(
                &func_,
                a_ + interval as f64 * h,
                a_ + (1 + interval) as f64 * h,
                eps_,
            )?;
            res_.push(sol);
        }

        res_.sort_by(|a, b| a.partial_cmp(b).unwrap());
        Ok(())
    }

    pub fn find_single_tetta_general<F>(func_: F, a_: f64, b_: f64, eps_: f64) -> Result<f64, Error>
    where
        F: Fn(f64) -> f64,
    {
        if b_ < a_ {
            return Err(Error::new(
                "Error in alp_reg::find_single_tetta_general\n".to_string(),
                4,
            ));
        }

        let mut x1 = a_;
        let mut x2 = b_;
        let mut precision = (x2 - x1) / 2.0;
        let mut y1 = func_(x1);
        if y1.abs() < eps_ {
            return Ok(x1);
        }

        let y2 = func_(x2);
        if y2.abs() < eps_ {
            return Ok(x2);
        }

        while precision > eps_ {
            let x12 = (x1 + x2) / 2.0;
            let y12 = func_(x12);
            if y12.abs() < eps_ {
                return Ok(x12);
            }
            if y12 * y1 < 0.0 {
                x2 = x12;
            } else {
                x1 = x12;
                y1 = y12;
            }
            precision = (x2 - x1) / 2.0;
        }

        Ok((x1 + x2) / 2.0)
    }

    pub fn correction_of_errors(
        errors_: &mut [f64],
        number_of_elements_: i64,
    ) -> Result<(), Error> {
        if number_of_elements_ <= 0 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        let n = number_of_elements_ as usize;
        let mut average_error = 0.0;
        for &err in &errors_[..n] {
            if err < 0.0 {
                return Err(Error::new(
                    "Error in alp_reg::correction_of_errors: input error in the regression model is less than 0\n".to_string(),
                    4,
                ));
            }
            average_error += err;
        }
        average_error /= number_of_elements_ as f64;

        let error_eps = if average_error <= 0.0 {
            1e-50
        } else {
            average_error
        };

        for err in &mut errors_[..n] {
            if *err == 0.0 {
                *err = error_eps;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn robust_regression_sum_with_cut_LSM(
        min_length_: i64,
        number_of_elements_: i64,
        values_: &[f64],
        errors_: &mut [f64],
        cut_left_tail_: bool,
        cut_right_tail_: bool,
        y_: f64,
        beta0_: &mut f64,
        beta1_: &mut f64,
        beta0_error_: &mut f64,
        beta1_error_: &mut f64,
        k1_opt_: &mut i64,
        k2_opt_: &mut i64,
        res_was_calculated_: &mut bool,
    ) -> Result<(), Error> {
        if number_of_elements_ < 2 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        Self::correction_of_errors(errors_, number_of_elements_)?;
        let c = y_ * y_;

        let (k1_start, k1_end, k2_start, k2_end) = if cut_left_tail_ && cut_right_tail_ {
            (0, number_of_elements_ - 1, 0, number_of_elements_ - 1)
        } else if cut_left_tail_ && !cut_right_tail_ {
            (
                0,
                number_of_elements_ - 1,
                number_of_elements_ - 1,
                number_of_elements_ - 1,
            )
        } else if !cut_left_tail_ && cut_right_tail_ {
            (0, 0, 0, number_of_elements_ - 1)
        } else {
            (0, 0, number_of_elements_ - 1, number_of_elements_ - 1)
        };

        let mut k1_opt = 0;
        let mut k2_opt = 0;
        let mut func_opt = f64::MAX;
        let mut beta0_opt = 0.0;
        let mut beta1_opt = 0.0;
        let mut beta0_opt_error = 0.0;
        let mut beta1_opt_error = 0.0;
        *res_was_calculated_ = false;

        for k1 in k1_start..=k1_end {
            let k2_begin = Tmax(k1 + 1, Tmax(k1, k2_start) + min_length_);
            for k2 in k2_begin..=k2_end {
                let mut beta0_opt_tmp = 0.0;
                let mut beta1_opt_tmp = 0.0;
                let mut beta0_opt_error_tmp = 0.0;
                let mut beta1_opt_error_tmp = 0.0;
                let mut res_was_calculated = false;
                let tmp = Self::function_for_robust_regression_sum_with_cut_LSM(
                    &values_[k1 as usize..],
                    &errors_[k1 as usize..],
                    k2 - k1 + 1,
                    k1,
                    c,
                    &mut beta0_opt_tmp,
                    &mut beta1_opt_tmp,
                    &mut beta0_opt_error_tmp,
                    &mut beta1_opt_error_tmp,
                    &mut res_was_calculated,
                );
                if tmp < func_opt && res_was_calculated {
                    func_opt = tmp;
                    beta0_opt = beta0_opt_tmp;
                    beta1_opt = beta1_opt_tmp;
                    beta0_opt_error = beta0_opt_error_tmp;
                    beta1_opt_error = beta1_opt_error_tmp;
                    k1_opt = k1;
                    k2_opt = k2;
                    *res_was_calculated_ = true;
                }
            }
        }

        if *res_was_calculated_ {
            *beta0_ = beta0_opt;
            *beta1_ = beta1_opt;
            *beta0_error_ = beta0_opt_error;
            *beta1_error_ = beta1_opt_error;
            *k1_opt_ = k1_opt;
            *k2_opt_ = k2_opt;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn function_for_robust_regression_sum_with_cut_LSM(
        values_: &[f64],
        errors_: &[f64],
        number_of_elements_: i64,
        k_start_: i64,
        c_: f64,
        beta0_: &mut f64,
        beta1_: &mut f64,
        beta0_error_: &mut f64,
        beta1_error_: &mut f64,
        res_was_calculated_: &mut bool,
    ) -> f64 {
        let mut a11 = 0.0;
        let mut a12 = 0.0;
        let mut a22 = 0.0;
        let mut y1 = 0.0;
        let mut y2 = 0.0;
        let mut y1_error = 0.0;
        let mut y2_error = 0.0;

        for i in 0..number_of_elements_ as usize {
            if errors_[i] != 0.0 {
                let tmp = 1.0 / (errors_[i] * errors_[i]);
                let k = k_start_ + i as i64;
                a11 += tmp;
                a12 += k as f64 * tmp;
                a22 += (k * k) as f64 * tmp;
                y1 += values_[i] * tmp;
                y1_error += tmp * tmp * errors_[i] * errors_[i];
                y2 += k as f64 * values_[i] * tmp;
                y2_error += k as f64 * k as f64 * tmp * tmp * errors_[i] * errors_[i];
            }
        }

        let a21 = a12;
        y1_error = Self::sqrt_for_errors(y1_error);
        y2_error = Self::sqrt_for_errors(y2_error);

        let eps = 1e-10 * Tmax((a11 * a22).abs(), (a21 * a12).abs());
        let den = a11 * a22 - a21 * a12;
        if den.abs() <= eps {
            *res_was_calculated_ = false;
            return 0.0;
        }
        *res_was_calculated_ = true;

        *beta0_ = (y1 * a22 - a12 * y2) / den;
        *beta1_ = (a11 * y2 - a21 * y1) / den;
        *beta0_error_ =
            (y1_error * y1_error * a22 * a22 + a12 * a12 * y2_error * y2_error).sqrt() / den;
        *beta1_error_ =
            (a11 * a11 * y2_error * y2_error + a21 * a21 * y1_error * y1_error).sqrt() / den;

        let mut res = 0.0;
        for i in 0..number_of_elements_ as usize {
            if errors_[i] != 0.0 {
                let tmp =
                    (*beta0_ + *beta1_ * (i as i64 + k_start_) as f64 - values_[i]) / errors_[i];
                res += tmp * tmp - c_;
            }
        }
        res
    }

    #[allow(clippy::too_many_arguments)]
    pub fn robust_regression_sum_with_cut_LSM_beta1_is_defined(
        min_length_: i64,
        number_of_elements_: i64,
        values_: &[f64],
        errors_: &mut [f64],
        cut_left_tail_: bool,
        cut_right_tail_: bool,
        y_: f64,
        beta0_: &mut f64,
        beta1_: f64,
        beta0_error_: &mut f64,
        beta1_error_: f64,
        k1_opt_: &mut i64,
        k2_opt_: &mut i64,
        res_was_calculated_: &mut bool,
    ) -> Result<(), Error> {
        Self::correction_of_errors(errors_, number_of_elements_)?;
        let c = y_ * y_;

        let (k1_start, k1_end, k2_start, k2_end) = if cut_left_tail_ && cut_right_tail_ {
            (0, number_of_elements_ - 1, 0, number_of_elements_ - 1)
        } else if cut_left_tail_ && !cut_right_tail_ {
            (
                0,
                number_of_elements_ - 1,
                number_of_elements_ - 1,
                number_of_elements_ - 1,
            )
        } else if !cut_left_tail_ && cut_right_tail_ {
            (0, 0, 0, number_of_elements_ - 1)
        } else {
            (0, 0, number_of_elements_ - 1, number_of_elements_ - 1)
        };

        let mut k1_opt = 0;
        let mut k2_opt = 0;
        let mut func_opt = f64::MAX;
        let mut beta0_opt = 0.0;
        let mut beta0_opt_error = 0.0;
        *res_was_calculated_ = false;

        for k1 in k1_start..=k1_end {
            let k2_begin = Tmax(k1, k2_start) + min_length_;
            for k2 in k2_begin..=k2_end {
                let mut beta0_opt_tmp = 0.0;
                let mut beta1_opt_tmp = beta1_;
                let mut beta0_opt_error_tmp = 0.0;
                let mut beta1_opt_error_tmp = beta1_error_;
                let mut res_was_calculated = false;
                let tmp = Self::function_for_robust_regression_sum_with_cut_LSM_beta1_is_defined(
                    &values_[k1 as usize..],
                    &errors_[k1 as usize..],
                    k2 - k1 + 1,
                    k1,
                    c,
                    &mut beta0_opt_tmp,
                    beta1_opt_tmp,
                    &mut beta0_opt_error_tmp,
                    beta1_opt_error_tmp,
                    &mut res_was_calculated,
                );
                if tmp < func_opt && res_was_calculated {
                    func_opt = tmp;
                    beta0_opt = beta0_opt_tmp;
                    beta0_opt_error = beta0_opt_error_tmp;
                    k1_opt = k1;
                    k2_opt = k2;
                    *res_was_calculated_ = true;
                }
                let _ = (&mut beta1_opt_tmp, &mut beta1_opt_error_tmp);
            }
        }

        if *res_was_calculated_ {
            *beta0_ = beta0_opt;
            *beta0_error_ = beta0_opt_error;
            *k1_opt_ = k1_opt;
            *k2_opt_ = k2_opt;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn function_for_robust_regression_sum_with_cut_LSM_beta1_is_defined(
        values_: &[f64],
        errors_: &[f64],
        number_of_elements_: i64,
        k_start_: i64,
        c_: f64,
        beta0_: &mut f64,
        beta1_: f64,
        beta0_error_: &mut f64,
        beta1_error_: f64,
        res_was_calculated_: &mut bool,
    ) -> f64 {
        let mut a11 = 0.0;
        let mut y1 = 0.0;
        let mut y1_error = 0.0;

        for i in 0..number_of_elements_ as usize {
            if errors_[i] != 0.0 {
                let tmp = 1.0 / (errors_[i] * errors_[i]);
                let k = k_start_ + i as i64;
                a11 += tmp;
                y1 += (values_[i] - k as f64 * beta1_) * tmp;
                let error_tmp =
                    errors_[i] * errors_[i] + k as f64 * k as f64 * beta1_error_ * beta1_error_;
                y1_error += tmp * tmp * error_tmp;
            }
        }

        y1_error = y1_error.sqrt();
        let eps = 1e-10 * a11.abs();
        let den = a11;
        if den.abs() <= eps {
            *res_was_calculated_ = false;
            return 0.0;
        }
        *res_was_calculated_ = true;
        *beta0_ = y1 / den;
        *beta0_error_ = y1_error / den;

        let mut res = 0.0;
        for i in 0..number_of_elements_ as usize {
            if errors_[i] != 0.0 {
                let tmp =
                    (*beta0_ + beta1_ * (i as i64 + k_start_) as f64 - values_[i]) / errors_[i];
                res += tmp * tmp - c_;
            }
        }
        res
    }

    pub fn sqrt_for_errors(x_: f64) -> f64 {
        if x_ <= 0.0 {
            0.0
        } else {
            x_.sqrt()
        }
    }

    pub fn median(dim_: i64, array_: &[f64]) -> f64 {
        let mut array_vect = array_[..dim_ as usize].to_vec();
        array_vect.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if dim_ % 2 == 0 {
            let k = round(&(dim_ as f64 / 2.0)) as usize;
            0.5 * (array_vect[k - 1] + array_vect[k])
        } else {
            let k = round(&((dim_ as f64 - 1.0) / 2.0)) as usize;
            array_vect[k]
        }
    }

    pub fn robust_sum(
        values: &[f64],
        dim: i64,
        N_points: i64,
        remove_flag: &mut Vec<bool>,
    ) -> Result<f64, Error> {
        remove_flag.clear();
        if dim <= N_points {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        let dim_usize = dim as usize;
        remove_flag.resize(dim_usize, true);
        let med_val = Self::median(dim, values);
        let mut array_vect = Vec::with_capacity(dim_usize);
        for i in 0..dim_usize {
            array_vect.push((-(values[i] - med_val).abs(), i as i64));
        }
        array_vect.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for i in 0..N_points as usize {
            remove_flag[array_vect[i].1 as usize] = false;
        }

        let mut res = 0.0;
        for i in 0..dim_usize {
            if remove_flag[i] {
                res += values[i];
            }
        }
        res /= (dim - N_points) as f64;
        Ok(res)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sqrt_for_errors() {
        assert_eq!(alp_reg::sqrt_for_errors(-1.0), 0.0);
        assert_eq!(alp_reg::sqrt_for_errors(0.0), 0.0);
        assert_eq!(alp_reg::sqrt_for_errors(9.0), 3.0);
    }

    #[test]
    fn test_find_tetta_general_and_single() {
        let root = alp_reg::find_single_tetta_general(|x| x * x - 4.0, 0.0, 3.0, 1e-9).unwrap();
        assert!((root - 2.0).abs() < 1e-6);

        let mut roots = Vec::new();
        alp_reg::find_tetta_general(|x| x * x - 1.0, -2.0, 2.0, 8, 1e-9, &mut roots).unwrap();
        assert_eq!(roots.len(), 2);
        assert!((roots[0] + 1.0).abs() < 1e-6);
        assert!((roots[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_correction_median_and_robust_sum() {
        let mut errors = vec![0.0, 2.0, 0.0, 4.0];
        alp_reg::correction_of_errors(&mut errors, 4).unwrap();
        assert_eq!(errors, vec![1.5, 2.0, 1.5, 4.0]);
        assert_eq!(alp_reg::median(4, &[4.0, 1.0, 3.0, 2.0]), 2.5);
        assert_eq!(alp_reg::median(3, &[4.0, 1.0, 3.0]), 3.0);

        let mut remove = Vec::new();
        let res = alp_reg::robust_sum(&[1.0, 2.0, 100.0, 4.0, 5.0], 5, 1, &mut remove).unwrap();
        assert_eq!(remove.iter().filter(|&&x| !x).count(), 1);
        assert_eq!(res, (1.0 + 2.0 + 4.0 + 5.0) / 4.0);
    }

    #[test]
    fn test_robust_regression_functions() {
        let values = vec![1.0, 3.0, 5.0, 7.0, 9.0];
        let errors = vec![1.0; 5];
        let mut beta0 = 0.0;
        let mut beta1 = 0.0;
        let mut beta0_error = 0.0;
        let mut beta1_error = 0.0;
        let mut calculated = false;
        let score = alp_reg::function_for_robust_regression_sum_with_cut_LSM(
            &values,
            &errors,
            5,
            0,
            0.0,
            &mut beta0,
            &mut beta1,
            &mut beta0_error,
            &mut beta1_error,
            &mut calculated,
        );
        assert!(calculated);
        assert!(score.abs() < 1e-12);
        assert!((beta0 - 1.0).abs() < 1e-12);
        assert!((beta1 - 2.0).abs() < 1e-12);

        let mut errors2 = vec![1.0; 5];
        let mut k1 = 0;
        let mut k2 = 0;
        alp_reg::robust_regression_sum_with_cut_LSM(
            2,
            5,
            &values,
            &mut errors2,
            false,
            true,
            0.0,
            &mut beta0,
            &mut beta1,
            &mut beta0_error,
            &mut beta1_error,
            &mut k1,
            &mut k2,
            &mut calculated,
        )
        .unwrap();
        assert!(calculated);
        assert!(k1 <= k2);
    }

    #[test]
    fn test_robust_regression_beta1_defined() {
        let values = vec![1.0, 3.0, 5.0, 7.0, 9.0];
        let errors = vec![1.0; 5];
        let mut beta0 = 0.0;
        let mut beta0_error = 0.0;
        let mut calculated = false;
        let score = alp_reg::function_for_robust_regression_sum_with_cut_LSM_beta1_is_defined(
            &values,
            &errors,
            5,
            0,
            0.0,
            &mut beta0,
            2.0,
            &mut beta0_error,
            0.0,
            &mut calculated,
        );
        assert!(calculated);
        assert!(score.abs() < 1e-12);
        assert!((beta0 - 1.0).abs() < 1e-12);

        let mut errors2 = vec![1.0; 5];
        let mut k1 = 0;
        let mut k2 = 0;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            2,
            5,
            &values,
            &mut errors2,
            false,
            true,
            0.0,
            &mut beta0,
            2.0,
            &mut beta0_error,
            0.0,
            &mut k1,
            &mut k2,
            &mut calculated,
        )
        .unwrap();
        assert!(calculated);
        assert!(k1 <= k2);
    }
}
