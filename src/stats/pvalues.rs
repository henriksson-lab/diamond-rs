use std::f64::consts::PI;

use crate::stats::sls_alp_data::alp_data;
use crate::stats::sls_alp_regression::alp_reg;
use crate::stats::sls_basic::{normal_probability, one_minus_exp_function, Error, CONST_VAL};

pub const NAT_CUT_OFF_IN_MAX: f64 = 2.0;

unsafe extern "C" {
    fn rand() -> i32;
}

const RAND_MAX_F64: f64 = 2147483647.0;

#[allow(non_snake_case, non_camel_case_types)]
#[derive(Debug, Clone, PartialEq)]
pub struct ALP_set_of_parameters {
    pub lambda: f64,
    pub lambda_error: f64,
    pub C: f64,
    pub C_error: f64,
    pub K: f64,
    pub K_error: f64,
    pub a_I: f64,
    pub a_I_error: f64,
    pub a_J: f64,
    pub a_J_error: f64,
    pub sigma: f64,
    pub sigma_error: f64,
    pub alpha_I: f64,
    pub alpha_I_error: f64,
    pub alpha_J: f64,
    pub alpha_J_error: f64,
    pub a: f64,
    pub a_error: f64,
    pub alpha: f64,
    pub alpha_error: f64,
    pub gapless_a: f64,
    pub gapless_a_error: f64,
    pub gapless_alpha: f64,
    pub gapless_alpha_error: f64,
    pub G: i64,
    pub G1: i64,
    pub G2: i64,
    pub m_LambdaSbs: Vec<f64>,
    pub m_KSbs: Vec<f64>,
    pub m_CSbs: Vec<f64>,
    pub m_SigmaSbs: Vec<f64>,
    pub m_AlphaISbs: Vec<f64>,
    pub m_AlphaJSbs: Vec<f64>,
    pub m_AISbs: Vec<f64>,
    pub m_AJSbs: Vec<f64>,
    pub m_CalcTime: f64,
    pub d_params_flag: bool,
    pub b_I: f64,
    pub b_I_error: f64,
    pub b_J: f64,
    pub b_J_error: f64,
    pub beta_I: f64,
    pub beta_I_error: f64,
    pub beta_J: f64,
    pub beta_J_error: f64,
    pub tau: f64,
    pub tau_error: f64,
    pub m_BISbs: Vec<f64>,
    pub m_BJSbs: Vec<f64>,
    pub m_BetaISbs: Vec<f64>,
    pub m_BetaJSbs: Vec<f64>,
    pub m_TauSbs: Vec<f64>,
    pub vi_y_thr: f64,
    pub vj_y_thr: f64,
    pub c_y_thr: f64,
}

impl Default for ALP_set_of_parameters {
    fn default() -> Self {
        Self {
            lambda: 0.0,
            lambda_error: 0.0,
            C: 0.0,
            C_error: 0.0,
            K: 0.0,
            K_error: 0.0,
            a_I: 0.0,
            a_I_error: 0.0,
            a_J: 0.0,
            a_J_error: 0.0,
            sigma: 0.0,
            sigma_error: 0.0,
            alpha_I: 0.0,
            alpha_I_error: 0.0,
            alpha_J: 0.0,
            alpha_J_error: 0.0,
            a: 0.0,
            a_error: 0.0,
            alpha: 0.0,
            alpha_error: 0.0,
            gapless_a: 0.0,
            gapless_a_error: 0.0,
            gapless_alpha: 0.0,
            gapless_alpha_error: 0.0,
            G: 0,
            G1: 0,
            G2: 0,
            m_LambdaSbs: Vec::new(),
            m_KSbs: Vec::new(),
            m_CSbs: Vec::new(),
            m_SigmaSbs: Vec::new(),
            m_AlphaISbs: Vec::new(),
            m_AlphaJSbs: Vec::new(),
            m_AISbs: Vec::new(),
            m_AJSbs: Vec::new(),
            m_CalcTime: 0.0,
            d_params_flag: false,
            b_I: 0.0,
            b_I_error: 0.0,
            b_J: 0.0,
            b_J_error: 0.0,
            beta_I: 0.0,
            beta_I_error: 0.0,
            beta_J: 0.0,
            beta_J_error: 0.0,
            tau: 0.0,
            tau_error: 0.0,
            m_BISbs: Vec::new(),
            m_BJSbs: Vec::new(),
            m_BetaISbs: Vec::new(),
            m_BetaJSbs: Vec::new(),
            m_TauSbs: Vec::new(),
            vi_y_thr: 0.0,
            vj_y_thr: 0.0,
            c_y_thr: 0.0,
        }
    }
}

#[allow(non_snake_case, non_camel_case_types)]
#[derive(Debug, Clone)]
pub struct pvalues {
    pub blast: bool,
    pub eps: f64,
    pub a_normal: f64,
    pub b_normal: f64,
    pub N_normal: i64,
    pub h_normal: f64,
}

#[allow(non_snake_case)]
impl pvalues {
    pub fn new() -> Self {
        let a_normal = -10.0;
        let b_normal = 10.0;
        let N_normal = 1000;
        Self {
            blast: false,
            eps: 0.0001,
            a_normal,
            b_normal,
            N_normal,
            h_normal: (b_normal - a_normal) / N_normal as f64,
        }
    }

    pub fn ran3() -> f64 {
        unsafe { rand() as f64 / RAND_MAX_F64 }
    }

    pub fn standard_normal() -> f64 {
        let mut r1 = 0.0;
        while r1 == 0.0 {
            r1 = Self::ran3();
        }
        let mut r2 = 0.0;
        while r2 == 0.0 {
            r2 = Self::ran3();
        }

        let mut v1 = -2.0 * r1.ln();
        if v1 < 0.0 {
            v1 = 0.0;
        }
        v1.sqrt() * (2.0 * crate::stats::sls_basic::PI * r2).cos()
    }

    pub fn get_appr_tail_prob_with_cov(
        par_: &ALP_set_of_parameters,
        blast_: bool,
        y_: f64,
        m_: f64,
        n_: f64,
        P_: &mut f64,
        P_error_: &mut f64,
        E_: &mut f64,
        E_error_: &mut f64,
        area_: &mut f64,
        area_is_1_flag_: &mut bool,
    ) {
        let blast_ = false && blast_;
        let lambda_ = par_.lambda;
        let lambda_error_ = par_.lambda_error;
        let k_ = par_.K;
        let k_error_ = par_.K_error;

        let ai_hat_ = par_.a_I;
        let ai_hat_error_ = par_.a_I_error;
        let mut alphai_hat_ = par_.alpha_I;
        let mut alphai_hat_error_ = par_.alpha_I_error;
        let aj_hat_ = par_.a_J;
        let aj_hat_error_ = par_.a_J_error;
        let mut alphaj_hat_ = par_.alpha_J;
        let mut alphaj_hat_error_ = par_.alpha_J_error;
        let mut sigma_hat_ = par_.sigma;
        let mut sigma_hat_error_ = par_.sigma_error;

        let bi_hat_ = par_.b_I;
        let bi_hat_error_ = par_.b_I_error;
        let mut betai_hat_ = par_.beta_I;
        let mut betai_hat_error_ = par_.beta_I_error;
        let bj_hat_ = par_.b_J;
        let bj_hat_error_ = par_.b_J_error;
        let mut betaj_hat_ = par_.beta_J;
        let mut betaj_hat_error_ = par_.beta_J_error;
        let mut tau_hat_ = par_.tau;
        let mut tau_hat_error_ = par_.tau_error;

        if blast_ {
            alphai_hat_ = 0.0;
            alphai_hat_error_ = 0.0;
            betai_hat_ = 0.0;
            betai_hat_error_ = 0.0;
            alphaj_hat_ = 0.0;
            alphaj_hat_error_ = 0.0;
            betaj_hat_ = 0.0;
            betaj_hat_error_ = 0.0;
            sigma_hat_ = 0.0;
            sigma_hat_error_ = 0.0;
            tau_hat_ = 0.0;
            tau_hat_error_ = 0.0;
        }

        let tmp = ai_hat_ * y_ + bi_hat_;
        let m_li_y_error = alp_data::error_of_the_sum(y_.abs() * ai_hat_error_, bi_hat_error_);
        let m_li_y = m_ - tmp;
        let vi_y_error = alp_data::error_of_the_sum(y_.abs() * alphai_hat_error_, betai_hat_error_);
        let vi_y = par_.vi_y_thr.max(alphai_hat_ * y_ + betai_hat_);
        let sqrt_vi_y_error = alp_data::error_of_the_sqrt(vi_y, vi_y_error);
        let sqrt_vi_y = vi_y.sqrt();
        let (m_F, m_F_error) = if sqrt_vi_y == 0.0 || blast_ {
            (1e100, 0.0)
        } else {
            (
                m_li_y / sqrt_vi_y,
                alp_data::error_of_the_ratio(m_li_y, m_li_y_error, sqrt_vi_y, sqrt_vi_y_error),
            )
        };
        let P_m_F = normal_probability(m_F);
        let P_m_F_error = CONST_VAL * (-0.5 * m_F * m_F).exp() * m_F_error;
        let E_m_F = -CONST_VAL * (-0.5 * m_F * m_F).exp();
        let E_m_F_error = (-E_m_F * m_F).abs() * m_F_error;
        let m_li_y_P_m_F_error =
            alp_data::error_of_the_product(m_li_y, m_li_y_error, P_m_F, P_m_F_error);
        let m_li_y_P_m_F = m_li_y * P_m_F;
        let sqrt_vi_y_E_m_F_error =
            alp_data::error_of_the_product(sqrt_vi_y, sqrt_vi_y_error, E_m_F, E_m_F_error);
        let sqrt_vi_y_E_m_F = sqrt_vi_y * E_m_F;
        let p1_error = alp_data::error_of_the_sum(m_li_y_P_m_F_error, sqrt_vi_y_E_m_F_error);
        let p1 = m_li_y_P_m_F - sqrt_vi_y_E_m_F;

        let tmp = aj_hat_ * y_ + bj_hat_;
        let n_lj_y_error = alp_data::error_of_the_sum(y_.abs() * aj_hat_error_, bj_hat_error_);
        let n_lj_y = n_ - tmp;
        let vj_y_error = alp_data::error_of_the_sum(y_.abs() * alphaj_hat_error_, betaj_hat_error_);
        let vj_y = par_.vj_y_thr.max(alphaj_hat_ * y_ + betaj_hat_);
        let sqrt_vj_y_error = alp_data::error_of_the_sqrt(vj_y, vj_y_error);
        let sqrt_vj_y = vj_y.sqrt();
        let (n_F, n_F_error) = if sqrt_vj_y == 0.0 || blast_ {
            (1e100, 0.0)
        } else {
            (
                n_lj_y / sqrt_vj_y,
                alp_data::error_of_the_ratio(n_lj_y, n_lj_y_error, sqrt_vj_y, sqrt_vj_y_error),
            )
        };
        let P_n_F = normal_probability(n_F);
        let P_n_F_error = CONST_VAL * (-0.5 * n_F * n_F).exp() * n_F_error;
        let E_n_F = -CONST_VAL * (-0.5 * n_F * n_F).exp();
        let E_n_F_error = (-E_n_F * n_F).abs() * n_F_error;
        let n_lj_y_P_n_F_error =
            alp_data::error_of_the_product(n_lj_y, n_lj_y_error, P_n_F, P_n_F_error);
        let n_lj_y_P_n_F = n_lj_y * P_n_F;
        let sqrt_vj_y_E_n_F_error =
            alp_data::error_of_the_product(sqrt_vj_y, sqrt_vj_y_error, E_n_F, E_n_F_error);
        let sqrt_vj_y_E_n_F = sqrt_vj_y * E_n_F;
        let p2_error = alp_data::error_of_the_sum(n_lj_y_P_n_F_error, sqrt_vj_y_E_n_F_error);
        let p2 = n_lj_y_P_n_F - sqrt_vj_y_E_n_F;

        let c_y_error = alp_data::error_of_the_sum(sigma_hat_error_ * y_, tau_hat_error_);
        let c_y = par_.c_y_thr.max(sigma_hat_ * y_ + tau_hat_);
        let P_m_F_P_n_F_error =
            alp_data::error_of_the_product(P_m_F, P_m_F_error, P_n_F, P_n_F_error);
        let P_m_F_P_n_F = P_m_F * P_n_F;
        let c_y_P_m_F_P_n_F_error =
            alp_data::error_of_the_product(c_y, c_y_error, P_m_F_P_n_F, P_m_F_P_n_F_error);
        let c_y_P_m_F_P_n_F = c_y * P_m_F_P_n_F;
        let p1_p2_error = alp_data::error_of_the_product(p1, p1_error, p2, p2_error);
        let p1_p2 = p1 * p2;
        let area_error = alp_data::error_of_the_sum(p1_p2_error, c_y_P_m_F_P_n_F_error);
        let mut area = p1_p2 + c_y_P_m_F_P_n_F;

        if blast_ {
            if area <= 1.0 {
                *area_is_1_flag_ = true;
            }
            if *area_is_1_flag_ {
                area = 1.0;
            }
        }

        let exp_lambda_y_error = (lambda_error_ * y_ * (-lambda_ * y_).exp()).abs();
        let exp_lambda_y = (-lambda_ * y_).exp();
        let k_exp_lambda_y_error =
            alp_data::error_of_the_product(k_, k_error_, exp_lambda_y, exp_lambda_y_error);
        let k_exp_lambda_y = k_ * exp_lambda_y;
        let area_k_exp_lambda_y_error =
            alp_data::error_of_the_product(area, area_error, k_exp_lambda_y, k_exp_lambda_y_error);
        let area_k_exp_lambda_y = -area * k_exp_lambda_y;
        *E_ = -area_k_exp_lambda_y;
        *E_error_ = area_k_exp_lambda_y_error;
        *P_error_ = area_k_exp_lambda_y.exp() * area_k_exp_lambda_y_error;
        *P_ = one_minus_exp_function(area_k_exp_lambda_y);
        *area_ = area;
    }

    pub fn compute_intercepts(par_: &mut ALP_set_of_parameters) -> Result<(), Error> {
        if !par_.d_params_flag {
            return Err(Error::new(
                "Unexpected error: pvalues::compute_intercepts is called for undefined parameters\n"
                    .to_string(),
                1,
            ));
        }

        par_.b_I = 2.0 * par_.G as f64 * (par_.gapless_a - par_.a_I);
        par_.b_I_error =
            2.0 * par_.G as f64 * alp_data::error_of_the_sum(par_.gapless_a_error, par_.a_I_error);
        par_.beta_I = 2.0 * par_.G as f64 * (par_.gapless_alpha - par_.alpha_I);
        par_.beta_I_error = 2.0
            * par_.G as f64
            * alp_data::error_of_the_sum(par_.gapless_alpha_error, par_.alpha_I_error);
        par_.b_J = 2.0 * par_.G as f64 * (par_.gapless_a - par_.a_J);
        par_.b_I_error =
            2.0 * par_.G as f64 * alp_data::error_of_the_sum(par_.gapless_a_error, par_.a_J_error);
        par_.beta_J = 2.0 * par_.G as f64 * (par_.gapless_alpha - par_.alpha_J);
        par_.beta_J_error = 2.0
            * par_.G as f64
            * alp_data::error_of_the_sum(par_.gapless_alpha_error, par_.alpha_J_error);
        par_.tau = 2.0 * par_.G as f64 * (par_.gapless_alpha - par_.sigma);
        par_.tau_error = 2.0
            * par_.G as f64
            * alp_data::error_of_the_sum(par_.gapless_alpha_error, par_.sigma_error);

        let vector_size = par_.m_LambdaSbs.len();
        par_.m_BISbs.resize(vector_size, 0.0);
        par_.m_BJSbs.resize(vector_size, 0.0);
        par_.m_BetaISbs.resize(vector_size, 0.0);
        par_.m_BetaJSbs.resize(vector_size, 0.0);
        par_.m_TauSbs.resize(vector_size, 0.0);

        for i in 0..vector_size {
            par_.m_BISbs[i] = 2.0 * par_.G as f64 * (par_.gapless_a - par_.m_AISbs[i]);
            par_.m_BetaISbs[i] = 2.0 * par_.G as f64 * (par_.gapless_alpha - par_.m_AlphaISbs[i]);
            par_.m_BJSbs[i] = 2.0 * par_.G as f64 * (par_.gapless_a - par_.m_AJSbs[i]);
            par_.m_BetaJSbs[i] = 2.0 * par_.G as f64 * (par_.gapless_alpha - par_.m_AlphaJSbs[i]);
            par_.m_TauSbs[i] = 2.0 * par_.G as f64 * (par_.gapless_alpha - par_.m_SigmaSbs[i]);
        }

        Self::compute_tmp_values(par_)
    }

    pub fn compute_tmp_values(par_: &mut ALP_set_of_parameters) -> Result<(), Error> {
        if !par_.d_params_flag {
            return Err(Error::new(
                "Unexpected call of pvalues::compute_tmp_values\n".to_string(),
                1,
            ));
        }

        if par_.lambda > 0.0 {
            par_.vi_y_thr = (NAT_CUT_OFF_IN_MAX * par_.alpha_I / par_.lambda).max(0.0);
            par_.vj_y_thr = (NAT_CUT_OFF_IN_MAX * par_.alpha_J / par_.lambda).max(0.0);
            par_.c_y_thr = (NAT_CUT_OFF_IN_MAX * par_.sigma / par_.lambda).max(0.0);
        } else {
            par_.vi_y_thr = 0.0;
            par_.vj_y_thr = 0.0;
            par_.c_y_thr = 0.0;
            par_.d_params_flag = false;
        }
        Ok(())
    }

    pub fn get_appr_tail_prob_with_cov_without_errors(
        par_: &ALP_set_of_parameters,
        blast_: bool,
        y_: f64,
        m_: f64,
        n_: f64,
        P_: &mut f64,
        E_: &mut f64,
        area_: &mut f64,
        area_is_1_flag_: &mut bool,
        compute_only_area_: bool,
    ) {
        let blast_ = false && blast_;
        let mut alphai_hat_ = par_.alpha_I;
        let mut alphaj_hat_ = par_.alpha_J;
        let mut sigma_hat_ = par_.sigma;
        let bi_hat_ = par_.b_I;
        let mut betai_hat_ = par_.beta_I;
        let bj_hat_ = par_.b_J;
        let mut betaj_hat_ = par_.beta_J;
        let mut tau_hat_ = par_.tau;

        if blast_ {
            alphai_hat_ = 0.0;
            betai_hat_ = 0.0;
            alphaj_hat_ = 0.0;
            betaj_hat_ = 0.0;
            sigma_hat_ = 0.0;
            tau_hat_ = 0.0;
        }

        let m_li_y = m_ - (par_.a_I * y_ + bi_hat_);
        let vi_y = par_.vi_y_thr.max(alphai_hat_ * y_ + betai_hat_);
        let sqrt_vi_y = vi_y.sqrt();
        let m_F = if sqrt_vi_y == 0.0 || blast_ {
            1e100
        } else {
            m_li_y / sqrt_vi_y
        };
        let P_m_F = normal_probability(m_F);
        let E_m_F = -CONST_VAL * (-0.5 * m_F * m_F).exp();
        let p1 = m_li_y * P_m_F - sqrt_vi_y * E_m_F;

        let n_lj_y = n_ - (par_.a_J * y_ + bj_hat_);
        let vj_y = par_.vj_y_thr.max(alphaj_hat_ * y_ + betaj_hat_);
        let sqrt_vj_y = vj_y.sqrt();
        let n_F = if sqrt_vj_y == 0.0 || blast_ {
            1e100
        } else {
            n_lj_y / sqrt_vj_y
        };
        let P_n_F = normal_probability(n_F);
        let E_n_F = -CONST_VAL * (-0.5 * n_F * n_F).exp();
        let p2 = n_lj_y * P_n_F - sqrt_vj_y * E_n_F;

        let c_y = par_.c_y_thr.max(sigma_hat_ * y_ + tau_hat_);
        let mut area = p1 * p2 + c_y * P_m_F * P_n_F;
        if blast_ {
            if area <= 1.0 {
                *area_is_1_flag_ = true;
            }
            if *area_is_1_flag_ {
                area = 1.0;
            }
        }
        *area_ = area;

        if compute_only_area_ {
            return;
        }

        let area_k_exp_lambda_y = -area * par_.K * (-par_.lambda * y_).exp();
        *E_ = -area_k_exp_lambda_y;
        *P_ = one_minus_exp_function(area_k_exp_lambda_y);
    }

    fn log_sum(mut a: f64, mut b: f64) -> f64 {
        if a < b {
            std::mem::swap(&mut a, &mut b);
        }
        a + (1.0 + (b - a).exp()).ln()
    }

    fn log_diff(a: f64, b: f64) -> f64 {
        if a < b {
            panic!("Taking log of negative number.");
        }
        a + (1.0 - (b - a).exp()).ln()
    }

    pub fn log_area(
        &self,
        par_: &ALP_set_of_parameters,
        blast_: bool,
        y_: f64,
        m_: f64,
        n_: f64,
    ) -> f64 {
        let blast_ = false && blast_;
        let mut alphai_hat_ = par_.alpha_I;
        let mut alphaj_hat_ = par_.alpha_J;
        let mut sigma_hat_ = par_.sigma;
        let bi_hat_ = par_.b_I;
        let mut betai_hat_ = par_.beta_I;
        let bj_hat_ = par_.b_J;
        let mut betaj_hat_ = par_.beta_J;
        let mut tau_hat_ = par_.tau;

        if blast_ {
            alphai_hat_ = 0.0;
            betai_hat_ = 0.0;
            alphaj_hat_ = 0.0;
            betaj_hat_ = 0.0;
            sigma_hat_ = 0.0;
            tau_hat_ = 0.0;
        }

        let m_li_y = m_ - (par_.a_I * y_ + bi_hat_);
        let vi_y = par_.vi_y_thr.max(alphai_hat_ * y_ + betai_hat_);
        let sqrt_vi_y = vi_y.sqrt();
        let m_F = if sqrt_vi_y == 0.0 || blast_ {
            1e100
        } else {
            m_li_y / sqrt_vi_y
        };
        let log_P_m_F = normal_probability(m_F).ln();
        let log_minus_E_m_F = CONST_VAL.ln() + (-0.5 * m_F * m_F);
        let log_minus_sqrt_vi_y_E_m_F = sqrt_vi_y.ln() + log_minus_E_m_F;
        let log_p1 = if m_li_y < 0.0 {
            let log_minus_m_li_y_P_m_F = (-m_li_y).ln() + log_P_m_F;
            Self::log_diff(log_minus_sqrt_vi_y_E_m_F, log_minus_m_li_y_P_m_F)
        } else {
            let log_m_li_y_P_m_F = m_li_y.ln() + log_P_m_F;
            Self::log_sum(log_minus_sqrt_vi_y_E_m_F, log_m_li_y_P_m_F)
        };

        let n_lj_y = n_ - (par_.a_J * y_ + bj_hat_);
        let vj_y = par_.vj_y_thr.max(alphaj_hat_ * y_ + betaj_hat_);
        let sqrt_vj_y = vj_y.sqrt();
        let n_F = if sqrt_vj_y == 0.0 || blast_ {
            1e100
        } else {
            n_lj_y / sqrt_vj_y
        };
        let log_P_n_F = normal_probability(n_F).ln();
        let log_minus_E_n_F = CONST_VAL.ln() + (-0.5 * n_F * n_F);
        let log_minus_sqrt_vj_y_E_n_F = sqrt_vj_y.ln() + log_minus_E_n_F;
        let log_p2 = if n_lj_y < 0.0 {
            let log_minus_n_lj_y_P_n_F = (-n_lj_y).ln() + log_P_n_F;
            Self::log_diff(log_minus_sqrt_vj_y_E_n_F, log_minus_n_lj_y_P_n_F)
        } else {
            let log_n_lj_y_P_n_F = n_lj_y.ln() + log_P_n_F;
            Self::log_sum(log_minus_sqrt_vj_y_E_n_F, log_n_lj_y_P_n_F)
        };

        let log_c_y = par_.c_y_thr.max(sigma_hat_ * y_ + tau_hat_).ln();
        let log_area = Self::log_sum(log_p1 + log_p2, log_c_y + log_P_m_F + log_P_n_F);
        if !log_area.is_finite() {
            panic!("Numerical error in area computation.");
        }
        log_area
    }

    pub fn get_P_error_using_splitting_method(
        par_: &ALP_set_of_parameters,
        blast_: bool,
        y_: f64,
        m_: f64,
        n_: f64,
        P_: &mut f64,
        P_error_: &mut f64,
        E_: &mut f64,
        E_error_: &mut f64,
        area_is_1_flag_: &mut bool,
    ) -> Result<(), Error> {
        let dim = par_.m_LambdaSbs.len();
        if dim == 0 {
            return Err(Error::new(
                "Unexpected error in get_P_error_using_splitting_method\n".to_string(),
                1,
            ));
        }

        *P_ = 0.0;
        *P_error_ = 0.0;
        *E_ = 0.0;
        *E_error_ = 0.0;
        let mut exp_E_values_aver = 0.0;
        let mut exp_E_values_error = 0.0;
        let mut P_values = vec![0.0; dim];
        let mut E_values = vec![0.0; dim];
        let mut exp_E_values = vec![0.0; dim];

        for i in 0..dim {
            let mut par_tmp = ALP_set_of_parameters {
                a_I: par_.m_AISbs[i],
                a_J: par_.m_AJSbs[i],
                gapless_a: par_.gapless_a,
                gapless_a_error: par_.gapless_a_error,
                sigma: par_.m_SigmaSbs[i],
                gapless_alpha: par_.gapless_alpha,
                gapless_alpha_error: par_.gapless_alpha_error,
                C: par_.m_CSbs[i],
                K: par_.m_KSbs[i],
                lambda: par_.m_LambdaSbs[i],
                alpha_I: par_.m_AlphaISbs[i],
                alpha_J: par_.m_AlphaJSbs[i],
                G: par_.G,
                G1: par_.G1,
                G2: par_.G2,
                b_I: par_.m_BISbs[i],
                b_J: par_.m_BJSbs[i],
                beta_I: par_.m_BetaISbs[i],
                beta_J: par_.m_BetaJSbs[i],
                tau: par_.m_TauSbs[i],
                d_params_flag: true,
                ..Default::default()
            };
            par_tmp.a = 0.5 * (par_tmp.a_I + par_tmp.a_J);
            par_tmp.alpha = 0.5 * (par_tmp.alpha_I + par_tmp.alpha_J);
            Self::compute_tmp_values(&mut par_tmp)?;

            let mut P_tmp = 0.0;
            let mut E_tmp = 0.0;
            let mut area_tmp = 0.0;
            Self::get_appr_tail_prob_with_cov_without_errors(
                &par_tmp,
                blast_,
                y_,
                m_,
                n_,
                &mut P_tmp,
                &mut E_tmp,
                &mut area_tmp,
                area_is_1_flag_,
                false,
            );
            P_values[i] = P_tmp;
            *P_ += P_tmp;
            E_values[i] = E_tmp;
            *E_ += E_tmp;
            let exp_E_tmp = (-E_tmp).exp();
            exp_E_values[i] = exp_E_tmp;
            exp_E_values_aver += exp_E_tmp;
        }

        if dim <= 1 || *P_ <= 0.0 || *E_ <= 0.0 {
            return Ok(());
        }

        *P_ /= dim as f64;
        *E_ /= dim as f64;
        exp_E_values_aver /= dim as f64;

        for i in 0..dim {
            if *P_ > 0.0 {
                let tmp = P_values[i] / *P_;
                *P_error_ += tmp * tmp;
            }
            if *E_ > 0.0 {
                let tmp = E_values[i] / *E_;
                *E_error_ += tmp * tmp;
            }
            if exp_E_values_aver > 0.0 {
                let tmp = exp_E_values[i] / exp_E_values_aver;
                exp_E_values_error += tmp * tmp;
            }
        }

        *P_error_ = *P_error_ / dim as f64 - 1.0;
        *E_error_ = *E_error_ / dim as f64 - 1.0;
        exp_E_values_error = exp_E_values_error / dim as f64 - 1.0;
        if *P_ < 1e-4 {
            *P_error_ = *P_ * alp_reg::sqrt_for_errors(*P_error_ / dim as f64);
        } else {
            *P_error_ =
                exp_E_values_aver * alp_reg::sqrt_for_errors(exp_E_values_error / dim as f64);
        }
        *E_error_ = *E_ * alp_reg::sqrt_for_errors(*E_error_ / dim as f64);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn calculate_P_values_range(
        &self,
        Score1: i64,
        Score2: i64,
        Seq1Len: f64,
        Seq2Len: f64,
        ParametersSet: &ALP_set_of_parameters,
        P_values: &mut Vec<f64>,
        P_values_errors: &mut Vec<f64>,
        E_values: &mut Vec<f64>,
        E_values_errors: &mut Vec<f64>,
    ) -> Result<(), Error> {
        if Score2 < Score1 {
            return Err(Error::new("Error - Score2<Score1\n".to_string(), 2));
        }
        if Seq1Len <= 0.0 || Seq2Len <= 0.0 {
            return Err(Error::new(
                "Error - Seq1Len<=0||Seq2Len<=0\n".to_string(),
                2,
            ));
        }

        let len = (Score2 - Score1 + 1) as usize;
        P_values.resize(len, 0.0);
        P_values_errors.resize(len, 0.0);
        E_values.resize(len, 0.0);
        E_values_errors.resize(len, 0.0);
        for y in Score1..=Score2 {
            self.calculate_P_values(
                y as f64,
                Seq1Len,
                Seq2Len,
                ParametersSet,
                &mut P_values[(y - Score1) as usize],
                &mut P_values_errors[(y - Score1) as usize],
                &mut E_values[(y - Score1) as usize],
                &mut E_values_errors[(y - Score1) as usize],
                true,
            )?;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn calculate_P_values(
        &self,
        Score: f64,
        Seq1Len: f64,
        Seq2Len: f64,
        ParametersSet: &ALP_set_of_parameters,
        P_value: &mut f64,
        P_value_error: &mut f64,
        E_value: &mut f64,
        E_value_error: &mut f64,
        read_Sbs_par_flag: bool,
    ) -> Result<(), Error> {
        if Seq1Len <= 0.0 || Seq2Len <= 0.0 {
            return Err(Error::new(
                "Error - Seq1Len<=0||Seq2Len<=0\n".to_string(),
                2,
            ));
        }

        let mut P = 0.0;
        let mut P_error = 0.0;
        let mut E = 0.0;
        let mut E_error = 0.0;
        let mut area = 0.0;
        let mut area_is_1_flag = false;

        if read_Sbs_par_flag {
            Self::get_appr_tail_prob_with_cov_without_errors(
                ParametersSet,
                self.blast,
                Score,
                Seq1Len,
                Seq2Len,
                &mut P,
                &mut E,
                &mut area,
                &mut area_is_1_flag,
                false,
            );

            if !ParametersSet.m_LambdaSbs.is_empty() {
                let mut P_tmp = 0.0;
                let mut E_tmp = 0.0;
                Self::get_P_error_using_splitting_method(
                    ParametersSet,
                    self.blast,
                    Score,
                    Seq1Len,
                    Seq2Len,
                    &mut P_tmp,
                    &mut P_error,
                    &mut E_tmp,
                    &mut E_error,
                    &mut area_is_1_flag,
                )?;
                if P_tmp > 0.0 {
                    P_error = P_error / P_tmp * P;
                }
                *P_value_error = P_error;
                if E_tmp > 0.0 {
                    E_error = E_error / E_tmp * E;
                }
                *E_value_error = E_error;
            } else {
                *P_value_error = -f64::MAX;
                *E_value_error = -f64::MAX;
            }
        } else {
            Self::get_appr_tail_prob_with_cov(
                ParametersSet,
                self.blast,
                Score,
                Seq1Len,
                Seq2Len,
                &mut P,
                &mut P_error,
                &mut E,
                &mut E_error,
                &mut area,
                &mut area_is_1_flag,
            );
            *P_value_error = P_error;
            *E_value_error = E_error;
        }

        *P_value = P;
        *E_value = E;
        Ok(())
    }

    pub fn assert_Gumbel_parameters(par_: &ALP_set_of_parameters) -> bool {
        if !(par_.lambda > 0.0)
            || par_.lambda_error < 0.0
            || !(par_.K > 0.0)
            || par_.K_error < 0.0
            || par_.a_I < 0.0
            || par_.a_I_error < 0.0
            || par_.a_J < 0.0
            || par_.a_J_error < 0.0
            || par_.sigma < 0.0
            || par_.sigma_error < 0.0
            || par_.alpha_I < 0.0
            || par_.alpha_I_error < 0.0
            || par_.alpha_J < 0.0
            || par_.alpha_J_error < 0.0
            || par_.gapless_a < 0.0
            || par_.gapless_a_error < 0.0
            || par_.gapless_alpha < 0.0
            || par_.gapless_alpha_error < 0.0
            || par_.G < 0
            || par_.G1 < 0
            || par_.G2 < 0
            || par_.b_I_error < 0.0
            || par_.b_J_error < 0.0
            || par_.beta_I_error < 0.0
            || par_.beta_J_error < 0.0
            || par_.tau_error < 0.0
        {
            return false;
        }

        let size_tmp = par_.m_LambdaSbs.len();
        par_.m_KSbs.len() == size_tmp
            && par_.m_SigmaSbs.len() == size_tmp
            && par_.m_AlphaISbs.len() == size_tmp
            && par_.m_AlphaJSbs.len() == size_tmp
            && par_.m_AISbs.len() == size_tmp
            && par_.m_AJSbs.len() == size_tmp
            && par_.m_BISbs.len() == size_tmp
            && par_.m_BJSbs.len() == size_tmp
            && par_.m_BetaISbs.len() == size_tmp
            && par_.m_BetaJSbs.len() == size_tmp
            && par_.m_TauSbs.len() == size_tmp
    }
}

impl Default for pvalues {
    fn default() -> Self {
        Self::new()
    }
}

/// Parameters for the ALP finite-size correction.
pub struct AreaParams {
    pub a_i: f64,
    pub b_i: f64,
    pub alpha_i: f64,
    pub beta_i: f64,
    pub a_j: f64,
    pub b_j: f64,
    pub alpha_j: f64,
    pub beta_j: f64,
    pub sigma: f64,
    pub tau: f64,
}

/// Compute the effective search area using the ALP finite-size correction.
///
/// This implements the `get_appr_tail_prob_with_cov_without_errors` function
/// from the ALP library.
pub fn compute_area(p: &AreaParams, score: f64, query_len: f64, subject_len: f64) -> f64 {
    let const_val = 1.0 / (2.0 * PI).sqrt();

    // Query length correction
    let m_li_y = query_len - (p.a_i * score + p.b_i);
    let vi_y = (p.alpha_i * score + p.beta_i).max(0.0);
    let sqrt_vi_y = vi_y.sqrt();

    let m_f = if sqrt_vi_y == 0.0 {
        1e100
    } else {
        m_li_y / sqrt_vi_y
    };

    let p_m_f = normal_cdf(m_f);
    let e_m_f = -const_val * (-0.5 * m_f * m_f).exp();
    let p1 = m_li_y * p_m_f - sqrt_vi_y * e_m_f;

    // Subject length correction
    let n_lj_y = subject_len - (p.a_j * score + p.b_j);
    let vj_y = (p.alpha_j * score + p.beta_j).max(0.0);
    let sqrt_vj_y = vj_y.sqrt();

    let n_f = if sqrt_vj_y == 0.0 {
        1e100
    } else {
        n_lj_y / sqrt_vj_y
    };

    let p_n_f = normal_cdf(n_f);
    let e_n_f = -const_val * (-0.5 * n_f * n_f).exp();
    let p2 = n_lj_y * p_n_f - sqrt_vj_y * e_n_f;

    // Covariance correction
    let c_y = (p.sigma * score + p.tau).max(0.0);

    let area = p1 * p2 + c_y * p_m_f * p_n_f;
    area.max(0.0)
}

/// Compute E-value using the full ALP area computation.
pub fn evalue_with_area(
    lambda: f64,
    k: f64,
    params: &AreaParams,
    score: f64,
    query_len: f64,
    subject_len: f64,
) -> f64 {
    let area = compute_area(params, score, query_len, subject_len);
    area * k * (-lambda * score).exp()
}

/// Standard normal CDF.
fn normal_cdf(x: f64) -> f64 {
    0.5 * erfc(-x / std::f64::consts::SQRT_2)
}

/// Complementary error function (Abramowitz & Stegun approximation).
fn erfc(x: f64) -> f64 {
    if x >= 0.0 {
        erfc_positive(x)
    } else {
        2.0 - erfc_positive(-x)
    }
}

fn erfc_positive(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normal_cdf() {
        assert!((normal_cdf(0.0) - 0.5).abs() < 0.001);
        assert!(normal_cdf(3.0) > 0.99);
        assert!(normal_cdf(-3.0) < 0.01);
    }

    #[test]
    fn test_erfc() {
        assert!((erfc(0.0) - 1.0).abs() < 0.001);
        assert!(erfc(3.0) < 0.001);
    }

    #[test]
    fn test_compute_area_basic() {
        let a = 1.9 / 0.267;
        let params = AreaParams {
            a_i: a,
            b_i: 0.0,
            alpha_i: 1.9,
            beta_i: 0.0,
            a_j: a,
            b_j: 0.0,
            alpha_j: 1.9,
            beta_j: 0.0,
            sigma: 43.0,
            tau: 0.0,
        };
        let area = compute_area(&params, 20.0, 1000.0, 10000.0);
        assert!(area > 0.0, "Area for low score on large seqs: {}", area);
    }

    #[test]
    fn test_evalue_high_score() {
        let params = AreaParams {
            a_i: 1.9,
            b_i: 0.0,
            alpha_i: 1.9,
            beta_i: 0.0,
            a_j: 1.9,
            b_j: 0.0,
            alpha_j: 1.9,
            beta_j: 0.0,
            sigma: 43.0,
            tau: 0.0,
        };
        let e = evalue_with_area(0.267, 0.041, &params, 2328.0, 426.0, 426.0);
        assert!(e < 1e-200, "E-value for extreme score: {}", e);
    }

    fn original_params() -> ALP_set_of_parameters {
        let mut par = ALP_set_of_parameters {
            lambda: 0.267,
            lambda_error: 0.001,
            K: 0.041,
            K_error: 0.002,
            a_I: 1.9,
            a_I_error: 0.03,
            a_J: 2.1,
            a_J_error: 0.04,
            sigma: 43.0,
            sigma_error: 0.5,
            alpha_I: 1.7,
            alpha_I_error: 0.03,
            alpha_J: 1.8,
            alpha_J_error: 0.04,
            gapless_a: 3.0,
            gapless_a_error: 0.1,
            gapless_alpha: 5.0,
            gapless_alpha_error: 0.2,
            G: 11,
            G1: 11,
            G2: 11,
            d_params_flag: true,
            m_LambdaSbs: vec![0.26, 0.27],
            m_KSbs: vec![0.04, 0.042],
            m_CSbs: vec![1.0, 1.1],
            m_SigmaSbs: vec![42.0, 44.0],
            m_AlphaISbs: vec![1.6, 1.8],
            m_AlphaJSbs: vec![1.7, 1.9],
            m_AISbs: vec![1.8, 2.0],
            m_AJSbs: vec![2.0, 2.2],
            ..Default::default()
        };
        pvalues::compute_intercepts(&mut par).unwrap();
        par
    }

    #[test]
    fn test_original_pvalues_compute_intercepts_and_tmp_values() {
        let par = original_params();
        assert!((par.b_I - 24.2).abs() < 1e-12);
        assert!((par.b_J - 19.8).abs() < 1e-12);
        assert!((par.beta_I - 72.6).abs() < 1e-12);
        assert!((par.beta_J - 70.4).abs() < 1e-12);
        assert_eq!(par.tau, -836.0);
        assert_eq!(par.m_BISbs, vec![26.4, 22.0]);
        assert!((par.m_BJSbs[0] - 22.0).abs() < 1e-12);
        assert!((par.m_BJSbs[1] - 17.6).abs() < 1e-12);
        assert_eq!(par.m_TauSbs, vec![-814.0, -858.0]);
        assert!((par.vi_y_thr - 12.734082397003745).abs() < 1e-15);
    }

    #[test]
    fn test_original_pvalues_tail_probability_methods() {
        let par = original_params();
        let mut p = 0.0;
        let mut e = 0.0;
        let mut area = 0.0;
        let mut area_is_1 = false;
        pvalues::get_appr_tail_prob_with_cov_without_errors(
            &par,
            false,
            50.0,
            1000.0,
            2000.0,
            &mut p,
            &mut e,
            &mut area,
            &mut area_is_1,
            false,
        );
        assert!(area > 0.0);
        assert!(e > 0.0);
        assert!(p > 0.0);

        let mut p2 = 0.0;
        let mut p2_error = 0.0;
        let mut e2 = 0.0;
        let mut e2_error = 0.0;
        let mut area2 = 0.0;
        pvalues::get_appr_tail_prob_with_cov(
            &par,
            false,
            50.0,
            1000.0,
            2000.0,
            &mut p2,
            &mut p2_error,
            &mut e2,
            &mut e2_error,
            &mut area2,
            &mut area_is_1,
        );
        assert!((p2 - p).abs() / p < 1e-12);
        assert!((e2 - e).abs() / e < 1e-12);
        assert_eq!(area2, area);
        assert!(p2_error > 0.0);
        assert!(e2_error > 0.0);
    }

    #[test]
    fn test_original_pvalues_calculate_range_and_assert() {
        let par = original_params();
        assert!(pvalues::assert_Gumbel_parameters(&par));
        let pv = pvalues::new();
        let mut p = 0.0;
        let mut p_error = 0.0;
        let mut e = 0.0;
        let mut e_error = 0.0;
        pv.calculate_P_values(
            50.0,
            1000.0,
            2000.0,
            &par,
            &mut p,
            &mut p_error,
            &mut e,
            &mut e_error,
            true,
        )
        .unwrap();
        assert!(p > 0.0 && e > 0.0);
        assert!(p_error >= 0.0 && e_error >= 0.0);

        let mut ps = Vec::new();
        let mut ps_err = Vec::new();
        let mut es = Vec::new();
        let mut es_err = Vec::new();
        pv.calculate_P_values_range(
            49,
            51,
            1000.0,
            2000.0,
            &par,
            &mut ps,
            &mut ps_err,
            &mut es,
            &mut es_err,
        )
        .unwrap();
        assert_eq!(ps.len(), 3);
        assert_eq!(es.len(), 3);
        assert!(ps[0] > ps[1] && ps[1] > ps[2]);
        assert!(es[0] > es[1] && es[1] > es[2]);
    }
}
