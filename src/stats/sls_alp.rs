#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use crate::stats::sls_alp_data::{alp_data, array, array_positive, mb_bytes, q_elem};
use crate::stats::sls_basic::{get_current_time, Error};

pub const small_long: i64 = i64::MIN / 2;

#[derive(Debug, Clone, PartialEq)]
pub struct state {
    pub d_cells_counts: array<i64>,
    pub d_HS_i_const_next: Vec<i64>,
    pub d_HI_i_const_next: Vec<i64>,
    pub d_HD_i_const_next: Vec<i64>,
    pub d_H_i_const_next: Vec<i64>,
    pub d_HS_j_const_next: Vec<i64>,
    pub d_HI_j_const_next: Vec<i64>,
    pub d_HD_j_const_next: Vec<i64>,
    pub d_H_j_const_next: Vec<i64>,
    pub d_HS_ij_next: i64,
    pub d_HI_ij_next: i64,
    pub d_HD_ij_next: i64,
    pub d_H_ij_next: i64,
    pub d_H_matr_len: i64,
    pub d_M: i64,
    pub d_sentinel_i_next: i64,
    pub d_sentinel_j_next: i64,
}

impl Default for state {
    fn default() -> Self {
        Self {
            d_cells_counts: array::new(None),
            d_HS_i_const_next: Vec::new(),
            d_HI_i_const_next: Vec::new(),
            d_HD_i_const_next: Vec::new(),
            d_H_i_const_next: Vec::new(),
            d_HS_j_const_next: Vec::new(),
            d_HI_j_const_next: Vec::new(),
            d_HD_j_const_next: Vec::new(),
            d_H_j_const_next: Vec::new(),
            d_HS_ij_next: 0,
            d_HI_ij_next: 0,
            d_HD_ij_next: 0,
            d_H_ij_next: 0,
            d_H_matr_len: 0,
            d_M: 0,
            d_sentinel_i_next: 0,
            d_sentinel_j_next: 0,
        }
    }
}

impl state {
    pub fn new() -> Self {
        Self::default()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct alp {
    pub d_alp_data: alp_data,
    pub d_a_step: i64,
    pub d_is_now: bool,
    pub d_seqi_len: i64,
    pub d_seqj_len: i64,
    pub d_seq_a_len: i64,
    pub d_H_matr_a_len: i64,
    pub d_W_matr_a_len: i64,
    pub d_seqi: Vec<i64>,
    pub d_seqj: Vec<i64>,
    pub d_H_matr_len: i64,
    pub d_W_matr_len: i64,
    pub d_WS_i_const_pred: Vec<f64>,
    pub d_WI_i_const_pred: Vec<f64>,
    pub d_WD_i_const_pred: Vec<f64>,
    pub d_WS_i_const_next: Vec<f64>,
    pub d_WI_i_const_next: Vec<f64>,
    pub d_WD_i_const_next: Vec<f64>,
    pub d_WS_j_const_pred: Vec<f64>,
    pub d_WI_j_const_pred: Vec<f64>,
    pub d_WD_j_const_pred: Vec<f64>,
    pub d_WS_j_const_next: Vec<f64>,
    pub d_WI_j_const_next: Vec<f64>,
    pub d_WD_j_const_next: Vec<f64>,
    pub d_WS_ij_pred: f64,
    pub d_WI_ij_pred: f64,
    pub d_WD_ij_pred: f64,
    pub d_WS_ij_next: f64,
    pub d_WI_ij_next: f64,
    pub d_WD_ij_next: f64,
    pub d_HS_i_const_pred: Vec<i64>,
    pub d_HI_i_const_pred: Vec<i64>,
    pub d_HD_i_const_pred: Vec<i64>,
    pub d_H_i_const_pred: Vec<i64>,
    pub d_HS_i_const_next: Vec<i64>,
    pub d_HI_i_const_next: Vec<i64>,
    pub d_HD_i_const_next: Vec<i64>,
    pub d_H_i_const_next: Vec<i64>,
    pub d_HS_j_const_pred: Vec<i64>,
    pub d_HI_j_const_pred: Vec<i64>,
    pub d_HD_j_const_pred: Vec<i64>,
    pub d_H_j_const_pred: Vec<i64>,
    pub d_HS_j_const_next: Vec<i64>,
    pub d_HI_j_const_next: Vec<i64>,
    pub d_HD_j_const_next: Vec<i64>,
    pub d_H_j_const_next: Vec<i64>,
    pub d_HS_ij_pred: i64,
    pub d_HI_ij_pred: i64,
    pub d_HD_ij_pred: i64,
    pub d_H_ij_pred: i64,
    pub d_HS_ij_next: i64,
    pub d_HI_ij_next: i64,
    pub d_HD_ij_next: i64,
    pub d_H_ij_next: i64,
    pub d_success: bool,
    pub d_H_edge_max: Vec<i64>,
    pub d_M: i64,
    pub d_nalp: i64,
    pub d_nalp_killing: i64,
    pub d_alp: array_positive<i64>,
    pub d_H_I: array_positive<i64>,
    pub d_H_J: array_positive<i64>,
    pub d_alp_pos: array_positive<i64>,
    pub d_alp_weights: array_positive<f64>,
    pub d_cells_counts: array<i64>,
    pub d_alp_states: Vec<Option<state>>,
    pub d_sentinel_i_next: i64,
    pub d_sentinel_j_next: i64,
    pub d_sentinel_i_pred: i64,
    pub d_sentinel_j_pred: i64,
    pub d_diff_opt: i64,
    pub d_sentinels_flag: bool,
    pub d_check_time_flag: bool,
    pub d_time_error_flag: bool,
    pub d_time_limit_flag: bool,
    pub d_single_realiztion_calculation_flag: bool,
    pub d_IS_state: u8,
}

impl alp {
    pub fn new(alp_data_: alp_data) -> Result<Self, Error> {
        let mut result = Self {
            d_alp_data: alp_data_,
            d_a_step: 30,
            d_is_now: true,
            d_seqi_len: 0,
            d_seqj_len: 0,
            d_seq_a_len: 0,
            d_H_matr_a_len: 0,
            d_W_matr_a_len: 0,
            d_seqi: Vec::new(),
            d_seqj: Vec::new(),
            d_H_matr_len: -1,
            d_W_matr_len: -1,
            d_WS_i_const_pred: Vec::new(),
            d_WI_i_const_pred: Vec::new(),
            d_WD_i_const_pred: Vec::new(),
            d_WS_i_const_next: Vec::new(),
            d_WI_i_const_next: Vec::new(),
            d_WD_i_const_next: Vec::new(),
            d_WS_j_const_pred: Vec::new(),
            d_WI_j_const_pred: Vec::new(),
            d_WD_j_const_pred: Vec::new(),
            d_WS_j_const_next: Vec::new(),
            d_WI_j_const_next: Vec::new(),
            d_WD_j_const_next: Vec::new(),
            d_WS_ij_pred: 0.0,
            d_WI_ij_pred: 0.0,
            d_WD_ij_pred: 0.0,
            d_WS_ij_next: 0.0,
            d_WI_ij_next: 0.0,
            d_WD_ij_next: 0.0,
            d_HS_i_const_pred: Vec::new(),
            d_HI_i_const_pred: Vec::new(),
            d_HD_i_const_pred: Vec::new(),
            d_H_i_const_pred: Vec::new(),
            d_HS_i_const_next: Vec::new(),
            d_HI_i_const_next: Vec::new(),
            d_HD_i_const_next: Vec::new(),
            d_H_i_const_next: Vec::new(),
            d_HS_j_const_pred: Vec::new(),
            d_HI_j_const_pred: Vec::new(),
            d_HD_j_const_pred: Vec::new(),
            d_H_j_const_pred: Vec::new(),
            d_HS_j_const_next: Vec::new(),
            d_HI_j_const_next: Vec::new(),
            d_HD_j_const_next: Vec::new(),
            d_H_j_const_next: Vec::new(),
            d_HS_ij_pred: 0,
            d_HI_ij_pred: 0,
            d_HD_ij_pred: 0,
            d_H_ij_pred: 0,
            d_HS_ij_next: 0,
            d_HI_ij_next: 0,
            d_HD_ij_next: 0,
            d_H_ij_next: 0,
            d_success: true,
            d_H_edge_max: vec![0],
            d_M: 0,
            d_nalp: -1,
            d_nalp_killing: 0,
            d_alp: array_positive::new(None)?,
            d_H_I: array_positive::new(None)?,
            d_H_J: array_positive::new(None)?,
            d_alp_pos: array_positive::new(None)?,
            d_alp_weights: array_positive::new(None)?,
            d_cells_counts: array::new(None),
            d_alp_states: Vec::new(),
            d_sentinel_i_next: 0,
            d_sentinel_j_next: 0,
            d_sentinel_i_pred: 0,
            d_sentinel_j_pred: 0,
            d_diff_opt: 0,
            d_sentinels_flag: false,
            d_check_time_flag: false,
            d_time_error_flag: false,
            d_time_limit_flag: false,
            d_single_realiztion_calculation_flag: false,
            d_IS_state: 0,
        };
        result.d_alp_data.d_memory_size_in_MB += std::mem::size_of::<i64>() as f64 / mb_bytes;
        result.increment_W_weights()?;
        result.increment_H_weights_with_sentinels(0)?;
        Ok(result)
    }

    pub fn random_AA1(&self) -> Result<i64, Error> {
        alp_data::random_long_element(
            self.d_alp_data.ran2(),
            self.d_alp_data.d_number_of_AA,
            &self.d_alp_data.d_RR1_sum,
            &self.d_alp_data.d_RR1_sum_elements,
        )
    }

    pub fn random_AA2(&self) -> Result<i64, Error> {
        alp_data::random_long_element(
            self.d_alp_data.ran2(),
            self.d_alp_data.d_number_of_AA,
            &self.d_alp_data.d_RR2_sum,
            &self.d_alp_data.d_RR2_sum_elements,
        )
    }

    pub fn one_step_of_importance_sampling_without_weight_calculation(
        &mut self,
        d_dim1_: i64,
        d_dim2_: i64,
        random_values: &mut impl Iterator<Item = f64>,
    ) -> Result<bool, Error> {
        let d_is = self.d_alp_data.d_is.as_ref().ok_or_else(|| {
            Error::new(
                "Unexpected error in alp::one_step_of_importance_sampling_without_weight_calculation\n".to_string(),
                4,
            )
        })?;
        if self.d_seqi_len == 0 && self.d_seqj_len == 0 {
            let r = random_values.next().unwrap_or(0.0);
            self.d_IS_state =
                alp_data::random_long_element(r, 3, &d_is.d_for_S, &d_is.d_for_S_states)?;
        }
        if self.d_IS_state == b'D' {
            if self.d_seqi_len == d_dim1_ {
                return Ok(false);
            }
            if self.d_seqi_len > self.d_seq_a_len - 1 {
                self.increment_sequences();
            }
            let r = random_values.next().unwrap_or(0.0);
            self.d_seqi[self.d_seqi_len as usize] = alp_data::random_long_element(
                r,
                self.d_alp_data.d_number_of_AA,
                &self.d_alp_data.d_RR1_sum,
                &self.d_alp_data.d_RR1_sum_elements,
            )?;
            self.d_seqi_len += 1;
            let r = random_values.next().unwrap_or(0.0);
            let d_is = self.d_alp_data.d_is.as_ref().unwrap();
            self.d_IS_state =
                alp_data::random_long_element(r, 3, &d_is.d_for_D, &d_is.d_for_D_states)?;
            return Ok(true);
        }
        if self.d_IS_state == b'I' {
            if self.d_seqj_len == d_dim2_ {
                return Ok(false);
            }
            if self.d_seqj_len > self.d_seq_a_len - 1 {
                self.increment_sequences();
            }
            let r = random_values.next().unwrap_or(0.0);
            self.d_seqj[self.d_seqj_len as usize] = alp_data::random_long_element(
                r,
                self.d_alp_data.d_number_of_AA,
                &self.d_alp_data.d_RR2_sum,
                &self.d_alp_data.d_RR2_sum_elements,
            )?;
            self.d_seqj_len += 1;
            let r = random_values.next().unwrap_or(0.0);
            let d_is = self.d_alp_data.d_is.as_ref().unwrap();
            self.d_IS_state =
                alp_data::random_long_element(r, 2, &d_is.d_for_I, &d_is.d_for_I_states)?;
            return Ok(true);
        }
        if self.d_IS_state == b'S' {
            if self.d_seqi_len == d_dim1_ || self.d_seqj_len == d_dim2_ {
                return Ok(false);
            }
            let r = random_values.next().unwrap_or(0.0);
            let pair: q_elem = {
                let d_is = self.d_alp_data.d_is.as_ref().unwrap();
                alp_data::random_long_element(
                    r,
                    d_is.d_is_number_of_AA * d_is.d_is_number_of_AA,
                    &d_is.d_elements_values,
                    &d_is.d_elements,
                )?
            };
            if self.d_seqi_len > self.d_seq_a_len - 1 || self.d_seqj_len > self.d_seq_a_len - 1 {
                self.increment_sequences();
            }
            self.d_seqi[self.d_seqi_len as usize] = pair.d_a;
            self.d_seqj[self.d_seqj_len as usize] = pair.d_b;
            self.d_seqi_len += 1;
            self.d_seqj_len += 1;
            let r = random_values.next().unwrap_or(0.0);
            let d_is = self.d_alp_data.d_is.as_ref().unwrap();
            self.d_IS_state =
                alp_data::random_long_element(r, 3, &d_is.d_for_S, &d_is.d_for_S_states)?;
            return Ok(true);
        }
        Ok(true)
    }

    pub fn increment_sequences(&mut self) {
        self.d_seq_a_len += self.d_a_step;
        self.d_seqi.resize(self.d_seq_a_len as usize, 0);
        self.d_seqj.resize(self.d_seq_a_len as usize, 0);
        self.d_alp_data.d_memory_size_in_MB +=
            (std::mem::size_of::<i64>() * self.d_a_step as usize * 2) as f64 / mb_bytes;
    }

    pub fn increment_W_matrix(&mut self) {
        let old_len = self.d_W_matr_a_len as usize;
        self.d_W_matr_a_len += self.d_a_step;
        let new_len = self.d_W_matr_a_len as usize;
        let next_copy_len = self.d_W_matr_len.max(0) as usize;
        let pred_copy_len = (self.d_W_matr_len - 1).max(0) as usize;

        let resize_next = |v: &mut Vec<f64>, len: usize| {
            let mut new_v = vec![0.0; new_len];
            new_v[..len].copy_from_slice(&v[..len]);
            *v = new_v;
        };
        let resize_pred = |v: &mut Vec<f64>, len: usize| {
            let mut new_v = vec![0.0; new_len];
            new_v[..len].copy_from_slice(&v[..len]);
            *v = new_v;
        };

        resize_pred(&mut self.d_WS_i_const_pred, pred_copy_len.min(old_len));
        resize_pred(&mut self.d_WI_i_const_pred, pred_copy_len.min(old_len));
        resize_pred(&mut self.d_WD_i_const_pred, pred_copy_len.min(old_len));
        resize_next(&mut self.d_WS_i_const_next, next_copy_len.min(old_len));
        resize_next(&mut self.d_WI_i_const_next, next_copy_len.min(old_len));
        resize_next(&mut self.d_WD_i_const_next, next_copy_len.min(old_len));
        resize_pred(&mut self.d_WS_j_const_pred, pred_copy_len.min(old_len));
        resize_pred(&mut self.d_WI_j_const_pred, pred_copy_len.min(old_len));
        resize_pred(&mut self.d_WD_j_const_pred, pred_copy_len.min(old_len));
        resize_next(&mut self.d_WS_j_const_next, next_copy_len.min(old_len));
        resize_next(&mut self.d_WI_j_const_next, next_copy_len.min(old_len));
        resize_next(&mut self.d_WD_j_const_next, next_copy_len.min(old_len));
        self.d_alp_data.d_memory_size_in_MB +=
            (std::mem::size_of::<f64>() * self.d_a_step as usize * 12) as f64 / mb_bytes;
    }

    pub fn increment_H_matrix(&mut self) {
        let old_len = self.d_H_matr_a_len as usize;
        self.d_H_matr_a_len += self.d_a_step;
        let new_len = self.d_H_matr_a_len as usize;
        let next_copy_len = self.d_H_matr_len.max(0) as usize;
        let pred_copy_len = (self.d_H_matr_len - 1).max(0) as usize;

        let resize_vec = |v: &mut Vec<i64>, len: usize| {
            let mut new_v = vec![0; new_len];
            new_v[..len].copy_from_slice(&v[..len]);
            *v = new_v;
        };
        resize_vec(&mut self.d_HS_i_const_pred, pred_copy_len.min(old_len));
        resize_vec(&mut self.d_HI_i_const_pred, pred_copy_len.min(old_len));
        resize_vec(&mut self.d_HD_i_const_pred, pred_copy_len.min(old_len));
        resize_vec(&mut self.d_H_i_const_pred, pred_copy_len.min(old_len));
        resize_vec(&mut self.d_HS_i_const_next, next_copy_len.min(old_len));
        resize_vec(&mut self.d_HI_i_const_next, next_copy_len.min(old_len));
        resize_vec(&mut self.d_HD_i_const_next, next_copy_len.min(old_len));
        resize_vec(&mut self.d_H_i_const_next, next_copy_len.min(old_len));
        resize_vec(&mut self.d_HS_j_const_pred, pred_copy_len.min(old_len));
        resize_vec(&mut self.d_HI_j_const_pred, pred_copy_len.min(old_len));
        resize_vec(&mut self.d_HD_j_const_pred, pred_copy_len.min(old_len));
        resize_vec(&mut self.d_H_j_const_pred, pred_copy_len.min(old_len));
        resize_vec(&mut self.d_HS_j_const_next, next_copy_len.min(old_len));
        resize_vec(&mut self.d_HI_j_const_next, next_copy_len.min(old_len));
        resize_vec(&mut self.d_HD_j_const_next, next_copy_len.min(old_len));
        resize_vec(&mut self.d_H_j_const_next, next_copy_len.min(old_len));

        let mut new_edge = vec![0; new_len + 1];
        let edge_copy_len = (self.d_H_matr_len + 1).max(0) as usize;
        new_edge[..edge_copy_len].copy_from_slice(&self.d_H_edge_max[..edge_copy_len]);
        self.d_H_edge_max = new_edge;
        self.d_alp_data.d_memory_size_in_MB +=
            (std::mem::size_of::<i64>() * self.d_a_step as usize * 17) as f64 / mb_bytes;
    }

    pub fn increment_W_weights(&mut self) -> Result<(), Error> {
        if self.d_W_matr_len == -1 {
            self.d_WS_ij_next = 1.0;
            self.d_WI_ij_next = 0.0;
            self.d_WD_ij_next = 0.0;
            self.d_W_matr_len += 1;
            self.d_alp_weights.set_elem(0, 1.0);
            return Ok(());
        }
        if self.d_seqi_len < self.d_W_matr_len + 1 || self.d_seqj_len < self.d_W_matr_len + 1 {
            return Err(Error::new(
                "Unexpected error in increment_W_weights\n".to_string(),
                4,
            ));
        }
        if self.d_W_matr_len + 1 > self.d_W_matr_a_len {
            self.increment_W_matrix();
        }
        self.d_W_matr_len += 1;
        std::mem::swap(&mut self.d_WS_i_const_pred, &mut self.d_WS_i_const_next);
        std::mem::swap(&mut self.d_WI_i_const_pred, &mut self.d_WI_i_const_next);
        std::mem::swap(&mut self.d_WD_i_const_pred, &mut self.d_WD_i_const_next);
        std::mem::swap(&mut self.d_WS_j_const_pred, &mut self.d_WS_j_const_next);
        std::mem::swap(&mut self.d_WI_j_const_pred, &mut self.d_WI_j_const_next);
        std::mem::swap(&mut self.d_WD_j_const_pred, &mut self.d_WD_j_const_next);
        self.d_WS_ij_pred = self.d_WS_ij_next;
        self.d_WI_ij_pred = self.d_WI_ij_next;
        self.d_WD_ij_pred = self.d_WD_ij_next;

        let is = self.d_alp_data.d_is.as_ref().ok_or_else(|| {
            Error::new("Unexpected error in increment_W_weights\n".to_string(), 4)
        })?;
        let len1 = self.d_W_matr_len - 1;
        let len2 = self.d_W_matr_len - 2;
        self.d_WS_i_const_next[len1 as usize] = 0.0;
        self.d_WS_j_const_next[len1 as usize] = 0.0;
        self.d_WI_i_const_next[len1 as usize] = 0.0;
        self.d_WD_j_const_next[len1 as usize] = 0.0;
        let deg_tmp = Self::degree(is.d_nu, len1 as f64)?;
        self.d_WD_i_const_next[len1 as usize] = is.d_mu_DS * deg_tmp;
        self.d_WI_j_const_next[len1 as usize] = is.d_mu_IS * deg_tmp;

        for i in (1..=len2).rev() {
            let iu = i as usize;
            self.d_WS_i_const_next[iu] = is.d_exp_s[self.d_seqi[len1 as usize] as usize]
                [self.d_seqj[(len2 - i) as usize] as usize]
                * (is.d_eta * self.d_WS_i_const_pred[iu]
                    + is.d_mu_SI * self.d_WI_i_const_pred[iu]
                    + is.d_mu_SD * self.d_WD_i_const_pred[iu]);
            self.d_WI_i_const_next[iu] = is.d_mu_IS * self.d_WS_i_const_next[(i + 1) as usize]
                + is.d_nu * self.d_WI_i_const_next[(i + 1) as usize]
                + is.d_mu_ID * self.d_WD_i_const_next[(i + 1) as usize];
            self.d_WD_i_const_next[iu] = is.d_mu_DS * self.d_WS_i_const_pred[(i - 1) as usize]
                + is.d_nu * self.d_WD_i_const_pred[(i - 1) as usize];

            self.d_WS_j_const_next[iu] = is.d_exp_s[self.d_seqi[(len2 - i) as usize] as usize]
                [self.d_seqj[len1 as usize] as usize]
                * (is.d_eta * self.d_WS_j_const_pred[iu]
                    + is.d_mu_SI * self.d_WI_j_const_pred[iu]
                    + is.d_mu_SD * self.d_WD_j_const_pred[iu]);
            self.d_WI_j_const_next[iu] = is.d_mu_IS * self.d_WS_j_const_pred[(i - 1) as usize]
                + is.d_nu * self.d_WI_j_const_pred[(i - 1) as usize]
                + is.d_mu_ID * self.d_WD_j_const_pred[(i - 1) as usize];
            self.d_WD_j_const_next[iu] = is.d_mu_DS * self.d_WS_j_const_next[(i + 1) as usize]
                + is.d_nu * self.d_WD_j_const_next[(i + 1) as usize];
        }
        if self.d_W_matr_len > 1 {
            let i = 0usize;
            self.d_WS_i_const_next[i] = is.d_exp_s[self.d_seqi[len1 as usize] as usize]
                [self.d_seqj[len2 as usize] as usize]
                * (is.d_eta * self.d_WS_i_const_pred[i]
                    + is.d_mu_SI * self.d_WI_i_const_pred[i]
                    + is.d_mu_SD * self.d_WD_i_const_pred[i]);
            self.d_WI_i_const_next[i] = is.d_mu_IS * self.d_WS_i_const_next[1]
                + is.d_nu * self.d_WI_i_const_next[1]
                + is.d_mu_ID * self.d_WD_i_const_next[1];
            self.d_WD_i_const_next[i] =
                is.d_mu_DS * self.d_WS_ij_pred + is.d_nu * self.d_WD_ij_pred;
            self.d_WS_j_const_next[i] = is.d_exp_s[self.d_seqi[len2 as usize] as usize]
                [self.d_seqj[len1 as usize] as usize]
                * (is.d_eta * self.d_WS_j_const_pred[i]
                    + is.d_mu_SI * self.d_WI_j_const_pred[i]
                    + is.d_mu_SD * self.d_WD_j_const_pred[i]);
            self.d_WI_j_const_next[i] = is.d_mu_IS * self.d_WS_ij_pred
                + is.d_nu * self.d_WI_ij_pred
                + is.d_mu_ID * self.d_WD_ij_pred;
            self.d_WD_j_const_next[i] =
                is.d_mu_DS * self.d_WS_j_const_next[1] + is.d_nu * self.d_WD_j_const_next[1];
        }
        self.d_WS_ij_next = is.d_exp_s[self.d_seqi[len1 as usize] as usize]
            [self.d_seqj[len1 as usize] as usize]
            * (is.d_eta * self.d_WS_ij_pred
                + is.d_mu_SI * self.d_WI_ij_pred
                + is.d_mu_SD * self.d_WD_ij_pred);
        self.d_WI_ij_next = is.d_mu_IS * self.d_WS_i_const_next[0]
            + is.d_nu * self.d_WI_i_const_next[0]
            + is.d_mu_ID * self.d_WD_i_const_next[0];
        self.d_WD_ij_next =
            is.d_mu_DS * self.d_WS_j_const_next[0] + is.d_nu * self.d_WD_j_const_next[0];
        Ok(())
    }

    pub fn degree(x_: f64, n_: f64) -> Result<f64, Error> {
        if x_ < 0.0 || n_ < 0.0 {
            return Err(Error::new(
                "Error - unexpected parameter in alp::degree\n".to_string(),
                4,
            ));
        }
        if x_ == 0.0 {
            return Ok(if n_ == 0.0 { 1.0 } else { 0.0 });
        }
        Ok((n_ * x_.ln()).exp())
    }

    pub fn increment_H_weights(&mut self) -> Result<(), Error> {
        if self.d_alp_data.d_insertions_after_deletions {
            self.increment_H_weights_with_insertions_after_deletions()
        } else {
            self.increment_H_weights_without_insertions_after_deletions()
        }
    }

    pub fn increment_H_weights_without_insertions_after_deletions(&mut self) -> Result<(), Error> {
        if self.d_H_matr_len == -1 {
            return self.increment_H_weights_with_sentinels_without_insertions_after_deletions(0);
        }
        if self.d_seqi_len < self.d_H_matr_len + 1 || self.d_seqj_len < self.d_H_matr_len + 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        if self.d_H_matr_len + 1 > self.d_H_matr_a_len {
            self.increment_H_matrix();
        }
        self.d_H_matr_len += 1;
        std::mem::swap(&mut self.d_HS_i_const_pred, &mut self.d_HS_i_const_next);
        std::mem::swap(&mut self.d_HI_i_const_pred, &mut self.d_HI_i_const_next);
        std::mem::swap(&mut self.d_HD_i_const_pred, &mut self.d_HD_i_const_next);
        std::mem::swap(&mut self.d_H_i_const_pred, &mut self.d_H_i_const_next);
        std::mem::swap(&mut self.d_HS_j_const_pred, &mut self.d_HS_j_const_next);
        std::mem::swap(&mut self.d_HI_j_const_pred, &mut self.d_HI_j_const_next);
        std::mem::swap(&mut self.d_HD_j_const_pred, &mut self.d_HD_j_const_next);
        std::mem::swap(&mut self.d_H_j_const_pred, &mut self.d_H_j_const_next);

        self.d_HS_ij_pred = self.d_HS_ij_next;
        self.d_HI_ij_pred = self.d_HI_ij_next;
        self.d_HD_ij_pred = self.d_HD_ij_next;
        self.d_H_ij_pred = self.d_H_ij_next;

        let len1 = self.d_H_matr_len - 1;
        let len2 = self.d_H_matr_len - 2;
        let gap_tmp1 = -self.d_alp_data.d_open1 - len1 * self.d_alp_data.d_epen1;
        let gap_tmp2 = -self.d_alp_data.d_open2 - len1 * self.d_alp_data.d_epen2;
        self.d_HS_i_const_next[len1 as usize] = small_long;
        self.d_HS_j_const_next[len1 as usize] = small_long;
        self.d_HI_i_const_next[len1 as usize] = small_long;
        self.d_HD_j_const_next[len1 as usize] = small_long;
        self.d_HD_i_const_next[len1 as usize] = gap_tmp1;
        self.d_HI_j_const_next[len1 as usize] = gap_tmp2;
        self.d_H_i_const_next[len1 as usize] = gap_tmp1;
        self.d_H_j_const_next[len1 as usize] = gap_tmp2;

        for i in (1..=len2).rev() {
            let iu = i as usize;
            self.d_HS_i_const_next[iu] = self.d_alp_data.d_smatr
                [self.d_seqi[len1 as usize] as usize][self.d_seqj[(len2 - i) as usize] as usize]
                + self.d_H_i_const_pred[iu];
            self.d_HI_i_const_next[iu] = (self.d_HS_i_const_next[(i + 1) as usize]
                - self.d_alp_data.d_open2)
                .max(self.d_HI_i_const_next[(i + 1) as usize] - self.d_alp_data.d_epen2);
            self.d_HD_i_const_next[iu] = (self.d_HS_i_const_pred[(i - 1) as usize]
                - self.d_alp_data.d_open1)
                .max(self.d_HD_i_const_pred[(i - 1) as usize] - self.d_alp_data.d_epen1);
            self.d_H_i_const_next[iu] = self.d_HS_i_const_next[iu]
                .max(self.d_HI_i_const_next[iu])
                .max(self.d_HD_i_const_next[iu]);

            self.d_HS_j_const_next[iu] = self.d_alp_data.d_smatr
                [self.d_seqi[(len2 - i) as usize] as usize][self.d_seqj[len1 as usize] as usize]
                + self.d_H_j_const_pred[iu];
            self.d_HI_j_const_next[iu] = (self.d_HS_j_const_pred[(i - 1) as usize]
                - self.d_alp_data.d_open2)
                .max(self.d_HI_j_const_pred[(i - 1) as usize] - self.d_alp_data.d_epen2);
            self.d_HD_j_const_next[iu] = (self.d_HS_j_const_next[(i + 1) as usize]
                - self.d_alp_data.d_open1)
                .max(self.d_HD_j_const_next[(i + 1) as usize] - self.d_alp_data.d_epen1);
            self.d_H_j_const_next[iu] = self.d_HS_j_const_next[iu]
                .max(self.d_HI_j_const_next[iu])
                .max(self.d_HD_j_const_next[iu]);
        }
        if self.d_H_matr_len > 1 {
            let i = 0usize;
            self.d_HS_i_const_next[i] = self.d_alp_data.d_smatr
                [self.d_seqi[len1 as usize] as usize][self.d_seqj[len2 as usize] as usize]
                + self.d_H_i_const_pred[i];
            self.d_HI_i_const_next[i] = (self.d_HS_i_const_next[1] - self.d_alp_data.d_open2)
                .max(self.d_HI_i_const_next[1] - self.d_alp_data.d_epen2);
            self.d_HD_i_const_next[i] = (self.d_HS_ij_pred - self.d_alp_data.d_open1)
                .max(self.d_HD_ij_pred - self.d_alp_data.d_epen1);
            self.d_H_i_const_next[i] = self.d_HS_i_const_next[i]
                .max(self.d_HI_i_const_next[i])
                .max(self.d_HD_i_const_next[i]);

            self.d_HS_j_const_next[i] = self.d_alp_data.d_smatr
                [self.d_seqi[len2 as usize] as usize][self.d_seqj[len1 as usize] as usize]
                + self.d_H_j_const_pred[i];
            self.d_HI_j_const_next[i] = (self.d_HS_ij_pred - self.d_alp_data.d_open2)
                .max(self.d_HI_ij_pred - self.d_alp_data.d_epen2);
            self.d_HD_j_const_next[i] = (self.d_HS_j_const_next[1] - self.d_alp_data.d_open1)
                .max(self.d_HD_j_const_next[1] - self.d_alp_data.d_epen1);
            self.d_H_j_const_next[i] = self.d_HS_j_const_next[i]
                .max(self.d_HI_j_const_next[i])
                .max(self.d_HD_j_const_next[i]);
        }

        self.d_HS_ij_next = self.d_alp_data.d_smatr[self.d_seqi[len1 as usize] as usize]
            [self.d_seqj[len1 as usize] as usize]
            + self.d_H_ij_pred;
        self.d_HI_ij_next = (self.d_HS_i_const_next[0] - self.d_alp_data.d_open2)
            .max(self.d_HI_i_const_next[0] - self.d_alp_data.d_epen2);
        self.d_HD_ij_next = (self.d_HS_j_const_next[0] - self.d_alp_data.d_open1)
            .max(self.d_HD_j_const_next[0] - self.d_alp_data.d_epen1);
        self.d_H_ij_next = self
            .d_HS_ij_next
            .max(self.d_HI_ij_next)
            .max(self.d_HD_ij_next);
        self.d_cells_counts.increase_elem_by_1(self.d_H_ij_next);
        for i in 0..=len1 {
            self.d_cells_counts
                .increase_elem_by_1(self.d_H_i_const_next[i as usize]);
            self.d_cells_counts
                .increase_elem_by_1(self.d_H_j_const_next[i as usize]);
        }
        let mut tmp = self.d_H_ij_next;
        for i in 0..=len1 {
            tmp = tmp.max(self.d_H_i_const_next[i as usize]);
            tmp = tmp.max(self.d_H_j_const_next[i as usize]);
        }
        self.d_H_edge_max[self.d_H_matr_len as usize] = tmp;
        self.d_M = tmp.max(self.d_M);
        self.d_sentinel_i_next = len1;
        self.d_sentinel_j_next = len1;

        if self.d_is_now && tmp > self.d_alp.d_elem[self.d_nalp as usize] {
            self.d_nalp += 1;
            self.d_alp.set_elem(self.d_nalp, tmp);
            self.d_alp_pos.set_elem(self.d_nalp, self.d_H_matr_len);
            let mut I = -1;
            let mut J = -1;
            for i in 0..=len1 {
                if tmp == self.d_H_i_const_next[i as usize] {
                    I = i;
                }
                if tmp == self.d_H_j_const_next[i as usize] {
                    J = i;
                }
            }
            self.d_H_I.set_elem(self.d_nalp, self.d_H_matr_len - I - 1);
            self.d_H_J.set_elem(self.d_nalp, self.d_H_matr_len - J - 1);
            let st = self.save_state()?;
            self.d_alp_states.resize(self.d_nalp as usize + 1, None);
            self.d_alp_states[self.d_nalp as usize] = Some(st);
        }
        self.check_time_function()
    }

    pub fn increment_H_weights_with_insertions_after_deletions(&mut self) -> Result<(), Error> {
        if self.d_H_matr_len == -1 {
            return self.increment_H_weights_with_sentinels_with_insertions_after_deletions(0);
        }
        if self.d_seqi_len < self.d_H_matr_len + 1 || self.d_seqj_len < self.d_H_matr_len + 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        if self.d_H_matr_len + 1 > self.d_H_matr_a_len {
            self.increment_H_matrix();
        }
        self.d_H_matr_len += 1;
        std::mem::swap(&mut self.d_HS_i_const_pred, &mut self.d_HS_i_const_next);
        std::mem::swap(&mut self.d_HI_i_const_pred, &mut self.d_HI_i_const_next);
        std::mem::swap(&mut self.d_HD_i_const_pred, &mut self.d_HD_i_const_next);
        std::mem::swap(&mut self.d_H_i_const_pred, &mut self.d_H_i_const_next);
        std::mem::swap(&mut self.d_HS_j_const_pred, &mut self.d_HS_j_const_next);
        std::mem::swap(&mut self.d_HI_j_const_pred, &mut self.d_HI_j_const_next);
        std::mem::swap(&mut self.d_HD_j_const_pred, &mut self.d_HD_j_const_next);
        std::mem::swap(&mut self.d_H_j_const_pred, &mut self.d_H_j_const_next);
        self.d_HS_ij_pred = self.d_HS_ij_next;
        self.d_HI_ij_pred = self.d_HI_ij_next;
        self.d_HD_ij_pred = self.d_HD_ij_next;
        self.d_H_ij_pred = self.d_H_ij_next;

        let len1 = self.d_H_matr_len - 1;
        let len2 = self.d_H_matr_len - 2;
        let gap_tmp1 = -self.d_alp_data.d_open1 - len1 * self.d_alp_data.d_epen1;
        let gap_tmp2 = -self.d_alp_data.d_open2 - len1 * self.d_alp_data.d_epen2;
        self.d_HS_i_const_next[len1 as usize] = small_long;
        self.d_HS_j_const_next[len1 as usize] = small_long;
        self.d_HI_i_const_next[len1 as usize] = small_long;
        self.d_HD_j_const_next[len1 as usize] = small_long;
        self.d_HD_i_const_next[len1 as usize] = gap_tmp1;
        self.d_HI_j_const_next[len1 as usize] = gap_tmp2;
        self.d_H_i_const_next[len1 as usize] = gap_tmp1;
        self.d_H_j_const_next[len1 as usize] = gap_tmp2;

        for i in (1..=len2).rev() {
            let iu = i as usize;
            self.d_HS_i_const_next[iu] = self.d_alp_data.d_smatr
                [self.d_seqi[len1 as usize] as usize][self.d_seqj[(len2 - i) as usize] as usize]
                + self.d_H_i_const_pred[iu];
            self.d_HI_i_const_next[iu] = (self.d_HS_i_const_next[(i + 1) as usize]
                - self.d_alp_data.d_open2)
                .max(self.d_HI_i_const_next[(i + 1) as usize] - self.d_alp_data.d_epen2)
                .max(self.d_HD_i_const_next[(i + 1) as usize] - self.d_alp_data.d_open2);
            self.d_HD_i_const_next[iu] = (self.d_HS_i_const_pred[(i - 1) as usize]
                - self.d_alp_data.d_open1)
                .max(self.d_HD_i_const_pred[(i - 1) as usize] - self.d_alp_data.d_epen1);
            self.d_H_i_const_next[iu] = self.d_HS_i_const_next[iu]
                .max(self.d_HI_i_const_next[iu])
                .max(self.d_HD_i_const_next[iu]);

            self.d_HS_j_const_next[iu] = self.d_alp_data.d_smatr
                [self.d_seqi[(len2 - i) as usize] as usize][self.d_seqj[len1 as usize] as usize]
                + self.d_H_j_const_pred[iu];
            self.d_HI_j_const_next[iu] = (self.d_HS_j_const_pred[(i - 1) as usize]
                - self.d_alp_data.d_open2)
                .max(self.d_HI_j_const_pred[(i - 1) as usize] - self.d_alp_data.d_epen2)
                .max(self.d_HD_j_const_pred[(i - 1) as usize] - self.d_alp_data.d_open2);
            self.d_HD_j_const_next[iu] = (self.d_HS_j_const_next[(i + 1) as usize]
                - self.d_alp_data.d_open1)
                .max(self.d_HD_j_const_next[(i + 1) as usize] - self.d_alp_data.d_epen1);
            self.d_H_j_const_next[iu] = self.d_HS_j_const_next[iu]
                .max(self.d_HI_j_const_next[iu])
                .max(self.d_HD_j_const_next[iu]);
        }
        if self.d_H_matr_len > 1 {
            let i = 0usize;
            self.d_HS_i_const_next[i] = self.d_alp_data.d_smatr
                [self.d_seqi[len1 as usize] as usize][self.d_seqj[len2 as usize] as usize]
                + self.d_H_i_const_pred[i];
            self.d_HI_i_const_next[i] = (self.d_HS_i_const_next[1] - self.d_alp_data.d_open2)
                .max(self.d_HI_i_const_next[1] - self.d_alp_data.d_epen2)
                .max(self.d_HD_i_const_next[1] - self.d_alp_data.d_open2);
            self.d_HD_i_const_next[i] = (self.d_HS_ij_pred - self.d_alp_data.d_open1)
                .max(self.d_HD_ij_pred - self.d_alp_data.d_epen1);
            self.d_H_i_const_next[i] = self.d_HS_i_const_next[i]
                .max(self.d_HI_i_const_next[i])
                .max(self.d_HD_i_const_next[i]);

            self.d_HS_j_const_next[i] = self.d_alp_data.d_smatr
                [self.d_seqi[len2 as usize] as usize][self.d_seqj[len1 as usize] as usize]
                + self.d_H_j_const_pred[i];
            self.d_HI_j_const_next[i] = (self.d_HS_ij_pred - self.d_alp_data.d_open2)
                .max(self.d_HI_ij_pred - self.d_alp_data.d_epen2)
                .max(self.d_HD_ij_pred - self.d_alp_data.d_open2);
            self.d_HD_j_const_next[i] = (self.d_HS_j_const_next[1] - self.d_alp_data.d_open1)
                .max(self.d_HD_j_const_next[1] - self.d_alp_data.d_epen1);
            self.d_H_j_const_next[i] = self.d_HS_j_const_next[i]
                .max(self.d_HI_j_const_next[i])
                .max(self.d_HD_j_const_next[i]);
        }
        self.d_HS_ij_next = self.d_alp_data.d_smatr[self.d_seqi[len1 as usize] as usize]
            [self.d_seqj[len1 as usize] as usize]
            + self.d_H_ij_pred;
        self.d_HI_ij_next = (self.d_HS_i_const_next[0] - self.d_alp_data.d_open2)
            .max(self.d_HI_i_const_next[0] - self.d_alp_data.d_epen2)
            .max(self.d_HD_i_const_next[0] - self.d_alp_data.d_open2);
        self.d_HD_ij_next = (self.d_HS_j_const_next[0] - self.d_alp_data.d_open1)
            .max(self.d_HD_j_const_next[0] - self.d_alp_data.d_epen1);
        self.d_H_ij_next = self
            .d_HS_ij_next
            .max(self.d_HI_ij_next)
            .max(self.d_HD_ij_next);
        self.d_cells_counts.increase_elem_by_1(self.d_H_ij_next);
        for i in 0..=len1 {
            self.d_cells_counts
                .increase_elem_by_1(self.d_H_i_const_next[i as usize]);
            self.d_cells_counts
                .increase_elem_by_1(self.d_H_j_const_next[i as usize]);
        }
        let mut tmp = self.d_H_ij_next;
        for i in 0..=len1 {
            tmp = tmp.max(self.d_H_i_const_next[i as usize]);
            tmp = tmp.max(self.d_H_j_const_next[i as usize]);
        }
        self.d_H_edge_max[self.d_H_matr_len as usize] = tmp;
        self.d_M = tmp.max(self.d_M);
        self.d_sentinel_i_next = len1;
        self.d_sentinel_j_next = len1;

        if self.d_is_now && tmp > self.d_alp.d_elem[self.d_nalp as usize] {
            self.d_nalp += 1;
            self.d_alp.set_elem(self.d_nalp, tmp);
            self.d_alp_pos.set_elem(self.d_nalp, self.d_H_matr_len);
            let mut I = -1;
            let mut J = -1;
            for i in 0..=len1 {
                if tmp == self.d_H_i_const_next[i as usize] {
                    I = i;
                }
                if tmp == self.d_H_j_const_next[i as usize] {
                    J = i;
                }
            }
            self.d_H_I.set_elem(self.d_nalp, self.d_H_matr_len - I - 1);
            self.d_H_J.set_elem(self.d_nalp, self.d_H_matr_len - J - 1);
            let st = self.save_state()?;
            self.d_alp_states.resize(self.d_nalp as usize + 1, None);
            self.d_alp_states[self.d_nalp as usize] = Some(st);
        }
        self.check_time_function()
    }

    pub fn increment_H_weights_with_sentinels(&mut self, diff_opt_: i64) -> Result<(), Error> {
        if self.d_alp_data.d_insertions_after_deletions {
            self.increment_H_weights_with_sentinels_with_insertions_after_deletions(diff_opt_)
        } else {
            self.increment_H_weights_with_sentinels_without_insertions_after_deletions(diff_opt_)
        }
    }

    pub fn increment_H_weights_with_sentinels_without_insertions_after_deletions(
        &mut self,
        diff_opt_: i64,
    ) -> Result<(), Error> {
        if self.d_H_matr_len == -1 {
            self.d_HS_ij_next = 0;
            self.d_HI_ij_next = 0;
            self.d_HD_ij_next = 0;
            self.d_H_ij_next = 0;
            self.d_M = 0;
            self.d_nalp = 0;
            self.d_alp.set_elem(0, 0);
            self.d_H_I.set_elem(0, 0);
            self.d_H_J.set_elem(0, 0);
            self.d_alp_pos.set_elem(0, 0);
            self.d_cells_counts.increase_elem_by_1(0);
            self.d_H_matr_len += 1;
            self.d_sentinel_i_next = 0;
            self.d_sentinel_j_next = 0;
            let st = self.save_state()?;
            self.d_alp_states.resize(1, None);
            self.d_alp_states[0] = Some(st);
            return Ok(());
        }
        if self.d_seqi_len < self.d_H_matr_len + 1 || self.d_seqj_len < self.d_H_matr_len + 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        if self.d_H_matr_len + 1 > self.d_H_matr_a_len {
            self.increment_H_matrix();
        }
        self.d_H_matr_len += 1;
        std::mem::swap(&mut self.d_HS_i_const_pred, &mut self.d_HS_i_const_next);
        std::mem::swap(&mut self.d_HI_i_const_pred, &mut self.d_HI_i_const_next);
        std::mem::swap(&mut self.d_HD_i_const_pred, &mut self.d_HD_i_const_next);
        std::mem::swap(&mut self.d_H_i_const_pred, &mut self.d_H_i_const_next);
        std::mem::swap(&mut self.d_HS_j_const_pred, &mut self.d_HS_j_const_next);
        std::mem::swap(&mut self.d_HI_j_const_pred, &mut self.d_HI_j_const_next);
        std::mem::swap(&mut self.d_HD_j_const_pred, &mut self.d_HD_j_const_next);
        std::mem::swap(&mut self.d_H_j_const_pred, &mut self.d_H_j_const_next);
        self.d_HS_ij_pred = self.d_HS_ij_next;
        self.d_HI_ij_pred = self.d_HI_ij_next;
        self.d_HD_ij_pred = self.d_HD_ij_next;
        self.d_H_ij_pred = self.d_H_ij_next;
        self.d_sentinel_i_pred = self.d_sentinel_i_next;
        self.d_sentinel_j_pred = self.d_sentinel_j_next;

        let len1 = self.d_H_matr_len - 1;
        let len2 = self.d_H_matr_len - 2;
        let sentinel_i_boundary = (self.d_sentinel_i_pred + 2).min(len1);
        let sentinel_j_boundary = (self.d_sentinel_j_pred + 2).min(len1);
        self.d_HS_i_const_next[sentinel_i_boundary as usize] = small_long;
        self.d_HS_j_const_next[sentinel_j_boundary as usize] = small_long;
        self.d_HI_i_const_next[sentinel_i_boundary as usize] = small_long;
        self.d_HD_j_const_next[sentinel_j_boundary as usize] = small_long;
        self.d_HD_i_const_next[sentinel_i_boundary as usize] = small_long;
        self.d_HI_j_const_next[sentinel_j_boundary as usize] = small_long;
        self.d_H_i_const_next[sentinel_i_boundary as usize] = small_long;
        self.d_H_j_const_next[sentinel_j_boundary as usize] = small_long;

        if sentinel_i_boundary > 1 {
            for i in (1..=sentinel_i_boundary - 1).rev() {
                let iu = i as usize;
                self.d_HS_i_const_next[iu] = self.d_alp_data.d_smatr
                    [self.d_seqi[len1 as usize] as usize]
                    [self.d_seqj[(len2 - i) as usize] as usize]
                    + self.d_H_i_const_pred[iu];
                self.d_HI_i_const_next[iu] = (self.d_HS_i_const_next[(i + 1) as usize]
                    - self.d_alp_data.d_open2)
                    .max(self.d_HI_i_const_next[(i + 1) as usize] - self.d_alp_data.d_epen2);
                self.d_HD_i_const_next[iu] = (self.d_HS_i_const_pred[(i - 1) as usize]
                    - self.d_alp_data.d_open1)
                    .max(self.d_HD_i_const_pred[(i - 1) as usize] - self.d_alp_data.d_epen1);
                self.d_H_i_const_next[iu] = self.d_HS_i_const_next[iu]
                    .max(self.d_HI_i_const_next[iu])
                    .max(self.d_HD_i_const_next[iu]);
            }
        }
        if sentinel_j_boundary > 1 {
            for i in (1..=sentinel_j_boundary - 1).rev() {
                let iu = i as usize;
                self.d_HS_j_const_next[iu] = self.d_alp_data.d_smatr
                    [self.d_seqi[(len2 - i) as usize] as usize]
                    [self.d_seqj[len1 as usize] as usize]
                    + self.d_H_j_const_pred[iu];
                self.d_HI_j_const_next[iu] = (self.d_HS_j_const_pred[(i - 1) as usize]
                    - self.d_alp_data.d_open2)
                    .max(self.d_HI_j_const_pred[(i - 1) as usize] - self.d_alp_data.d_epen2);
                self.d_HD_j_const_next[iu] = (self.d_HS_j_const_next[(i + 1) as usize]
                    - self.d_alp_data.d_open1)
                    .max(self.d_HD_j_const_next[(i + 1) as usize] - self.d_alp_data.d_epen1);
                self.d_H_j_const_next[iu] = self.d_HS_j_const_next[iu]
                    .max(self.d_HI_j_const_next[iu])
                    .max(self.d_HD_j_const_next[iu]);
            }
        }
        if self.d_H_matr_len > 1 {
            let i = 0usize;
            self.d_HS_i_const_next[i] = self.d_alp_data.d_smatr
                [self.d_seqi[len1 as usize] as usize][self.d_seqj[len2 as usize] as usize]
                + self.d_H_i_const_pred[i];
            self.d_HI_i_const_next[i] = (self.d_HS_i_const_next[1] - self.d_alp_data.d_open2)
                .max(self.d_HI_i_const_next[1] - self.d_alp_data.d_epen2);
            self.d_HD_i_const_next[i] = (self.d_HS_ij_pred - self.d_alp_data.d_open1)
                .max(self.d_HD_ij_pred - self.d_alp_data.d_epen1);
            self.d_H_i_const_next[i] = self.d_HS_i_const_next[i]
                .max(self.d_HI_i_const_next[i])
                .max(self.d_HD_i_const_next[i]);

            self.d_HS_j_const_next[i] = self.d_alp_data.d_smatr
                [self.d_seqi[len2 as usize] as usize][self.d_seqj[len1 as usize] as usize]
                + self.d_H_j_const_pred[i];
            self.d_HI_j_const_next[i] = (self.d_HS_ij_pred - self.d_alp_data.d_open2)
                .max(self.d_HI_ij_pred - self.d_alp_data.d_epen2);
            self.d_HD_j_const_next[i] = (self.d_HS_j_const_next[1] - self.d_alp_data.d_open1)
                .max(self.d_HD_j_const_next[1] - self.d_alp_data.d_epen1);
            self.d_H_j_const_next[i] = self.d_HS_j_const_next[i]
                .max(self.d_HI_j_const_next[i])
                .max(self.d_HD_j_const_next[i]);
        }
        self.d_HS_ij_next = self.d_alp_data.d_smatr[self.d_seqi[len1 as usize] as usize]
            [self.d_seqj[len1 as usize] as usize]
            + self.d_H_ij_pred;
        self.d_HI_ij_next = (self.d_HS_i_const_next[0] - self.d_alp_data.d_open2)
            .max(self.d_HI_i_const_next[0] - self.d_alp_data.d_epen2);
        self.d_HD_ij_next = (self.d_HS_j_const_next[0] - self.d_alp_data.d_open1)
            .max(self.d_HD_j_const_next[0] - self.d_alp_data.d_epen1);
        self.d_H_ij_next = self
            .d_HS_ij_next
            .max(self.d_HI_ij_next)
            .max(self.d_HD_ij_next);
        self.d_cells_counts.increase_elem_by_1(self.d_H_ij_next);
        if sentinel_i_boundary >= 1 {
            for i in 0..=sentinel_i_boundary - 1 {
                self.d_cells_counts
                    .increase_elem_by_1(self.d_H_i_const_next[i as usize]);
            }
        }
        if sentinel_j_boundary >= 1 {
            for i in 0..=sentinel_j_boundary - 1 {
                self.d_cells_counts
                    .increase_elem_by_1(self.d_H_j_const_next[i as usize]);
            }
        }
        let mut tmp = self.d_H_ij_next;
        if sentinel_i_boundary >= 1 {
            for i in 0..=sentinel_i_boundary - 1 {
                tmp = tmp.max(self.d_H_i_const_next[i as usize]);
            }
        }
        if sentinel_j_boundary >= 1 {
            for i in 0..=sentinel_j_boundary - 1 {
                tmp = tmp.max(self.d_H_j_const_next[i as usize]);
            }
        }
        self.d_H_edge_max[self.d_H_matr_len as usize] = tmp;
        self.d_M = tmp.max(self.d_M);
        let level = tmp - diff_opt_;
        self.d_sentinel_i_next = 1;
        self.d_sentinel_j_next = 1;
        if sentinel_i_boundary > 1 {
            for i in (1..=sentinel_i_boundary - 1).rev() {
                if self.d_H_i_const_next[i as usize] >= level {
                    self.d_sentinel_i_next = i;
                    break;
                }
            }
        }
        if sentinel_j_boundary > 1 {
            for i in (1..=sentinel_j_boundary - 1).rev() {
                if self.d_H_j_const_next[i as usize] >= level {
                    self.d_sentinel_j_next = i;
                    break;
                }
            }
        }
        if self.d_is_now && tmp > self.d_alp.d_elem[self.d_nalp as usize] {
            self.d_nalp += 1;
            self.d_alp.set_elem(self.d_nalp, tmp);
            self.d_alp_pos.set_elem(self.d_nalp, self.d_H_matr_len);
            let mut I = -1;
            let mut J = -1;
            if sentinel_i_boundary >= 1 {
                for i in 0..=sentinel_i_boundary - 1 {
                    if tmp == self.d_H_i_const_next[i as usize] {
                        I = i;
                    }
                }
            }
            if sentinel_j_boundary >= 1 {
                for i in 0..=sentinel_j_boundary - 1 {
                    if tmp == self.d_H_j_const_next[i as usize] {
                        J = i;
                    }
                }
            }
            self.d_H_I.set_elem(self.d_nalp, self.d_H_matr_len - I - 1);
            self.d_H_J.set_elem(self.d_nalp, self.d_H_matr_len - J - 1);
            let st = self.save_state()?;
            self.d_alp_states.resize(self.d_nalp as usize + 1, None);
            self.d_alp_states[self.d_nalp as usize] = Some(st);
        }
        self.check_time_function()
    }

    pub fn increment_H_weights_with_sentinels_with_insertions_after_deletions(
        &mut self,
        diff_opt_: i64,
    ) -> Result<(), Error> {
        if self.d_H_matr_len == -1 {
            self.d_HS_ij_next = 0;
            self.d_HI_ij_next = 0;
            self.d_HD_ij_next = 0;
            self.d_H_ij_next = 0;
            self.d_M = 0;
            self.d_nalp = 0;
            self.d_alp.set_elem(0, 0);
            self.d_H_I.set_elem(0, 0);
            self.d_H_J.set_elem(0, 0);
            self.d_alp_pos.set_elem(0, 0);
            self.d_cells_counts.increase_elem_by_1(0);
            self.d_H_matr_len += 1;
            self.d_sentinel_i_next = 0;
            self.d_sentinel_j_next = 0;
            let st = self.save_state()?;
            self.d_alp_states.resize(1, None);
            self.d_alp_states[0] = Some(st);
            return Ok(());
        }
        if self.d_seqi_len < self.d_H_matr_len + 1 || self.d_seqj_len < self.d_H_matr_len + 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        if self.d_H_matr_len + 1 > self.d_H_matr_a_len {
            self.increment_H_matrix();
        }
        self.d_H_matr_len += 1;
        std::mem::swap(&mut self.d_HS_i_const_pred, &mut self.d_HS_i_const_next);
        std::mem::swap(&mut self.d_HI_i_const_pred, &mut self.d_HI_i_const_next);
        std::mem::swap(&mut self.d_HD_i_const_pred, &mut self.d_HD_i_const_next);
        std::mem::swap(&mut self.d_H_i_const_pred, &mut self.d_H_i_const_next);
        std::mem::swap(&mut self.d_HS_j_const_pred, &mut self.d_HS_j_const_next);
        std::mem::swap(&mut self.d_HI_j_const_pred, &mut self.d_HI_j_const_next);
        std::mem::swap(&mut self.d_HD_j_const_pred, &mut self.d_HD_j_const_next);
        std::mem::swap(&mut self.d_H_j_const_pred, &mut self.d_H_j_const_next);
        self.d_HS_ij_pred = self.d_HS_ij_next;
        self.d_HI_ij_pred = self.d_HI_ij_next;
        self.d_HD_ij_pred = self.d_HD_ij_next;
        self.d_H_ij_pred = self.d_H_ij_next;
        self.d_sentinel_i_pred = self.d_sentinel_i_next;
        self.d_sentinel_j_pred = self.d_sentinel_j_next;

        let len1 = self.d_H_matr_len - 1;
        let len2 = self.d_H_matr_len - 2;
        let sentinel_i_boundary = (self.d_sentinel_i_pred + 2).min(len1);
        let sentinel_j_boundary = (self.d_sentinel_j_pred + 2).min(len1);
        self.d_HS_i_const_next[sentinel_i_boundary as usize] = small_long;
        self.d_HS_j_const_next[sentinel_j_boundary as usize] = small_long;
        self.d_HI_i_const_next[sentinel_i_boundary as usize] = small_long;
        self.d_HD_j_const_next[sentinel_j_boundary as usize] = small_long;
        self.d_HD_i_const_next[sentinel_i_boundary as usize] = small_long;
        self.d_HI_j_const_next[sentinel_j_boundary as usize] = small_long;
        self.d_H_i_const_next[sentinel_i_boundary as usize] = small_long;
        self.d_H_j_const_next[sentinel_j_boundary as usize] = small_long;

        if sentinel_i_boundary > 1 {
            for i in (1..=sentinel_i_boundary - 1).rev() {
                let iu = i as usize;
                self.d_HS_i_const_next[iu] = self.d_alp_data.d_smatr
                    [self.d_seqi[len1 as usize] as usize]
                    [self.d_seqj[(len2 - i) as usize] as usize]
                    + self.d_H_i_const_pred[iu];
                self.d_HI_i_const_next[iu] = (self.d_HS_i_const_next[(i + 1) as usize]
                    - self.d_alp_data.d_open2)
                    .max(self.d_HI_i_const_next[(i + 1) as usize] - self.d_alp_data.d_epen2)
                    .max(self.d_HD_i_const_next[(i + 1) as usize] - self.d_alp_data.d_open2);
                self.d_HD_i_const_next[iu] = (self.d_HS_i_const_pred[(i - 1) as usize]
                    - self.d_alp_data.d_open1)
                    .max(self.d_HD_i_const_pred[(i - 1) as usize] - self.d_alp_data.d_epen1);
                self.d_H_i_const_next[iu] = self.d_HS_i_const_next[iu]
                    .max(self.d_HI_i_const_next[iu])
                    .max(self.d_HD_i_const_next[iu]);
            }
        }
        if sentinel_j_boundary > 1 {
            for i in (1..=sentinel_j_boundary - 1).rev() {
                let iu = i as usize;
                self.d_HS_j_const_next[iu] = self.d_alp_data.d_smatr
                    [self.d_seqi[(len2 - i) as usize] as usize]
                    [self.d_seqj[len1 as usize] as usize]
                    + self.d_H_j_const_pred[iu];
                self.d_HI_j_const_next[iu] = (self.d_HS_j_const_pred[(i - 1) as usize]
                    - self.d_alp_data.d_open2)
                    .max(self.d_HI_j_const_pred[(i - 1) as usize] - self.d_alp_data.d_epen2)
                    .max(self.d_HD_j_const_pred[(i - 1) as usize] - self.d_alp_data.d_open2);
                self.d_HD_j_const_next[iu] = (self.d_HS_j_const_next[(i + 1) as usize]
                    - self.d_alp_data.d_open1)
                    .max(self.d_HD_j_const_next[(i + 1) as usize] - self.d_alp_data.d_epen1);
                self.d_H_j_const_next[iu] = self.d_HS_j_const_next[iu]
                    .max(self.d_HI_j_const_next[iu])
                    .max(self.d_HD_j_const_next[iu]);
            }
        }
        if self.d_H_matr_len > 1 {
            let i = 0usize;
            self.d_HS_i_const_next[i] = self.d_alp_data.d_smatr
                [self.d_seqi[len1 as usize] as usize][self.d_seqj[len2 as usize] as usize]
                + self.d_H_i_const_pred[i];
            self.d_HI_i_const_next[i] = (self.d_HS_i_const_next[1] - self.d_alp_data.d_open2)
                .max(self.d_HI_i_const_next[1] - self.d_alp_data.d_epen2)
                .max(self.d_HD_i_const_next[1] - self.d_alp_data.d_open2);
            self.d_HD_i_const_next[i] = (self.d_HS_ij_pred - self.d_alp_data.d_open1)
                .max(self.d_HD_ij_pred - self.d_alp_data.d_epen1);
            self.d_H_i_const_next[i] = self.d_HS_i_const_next[i]
                .max(self.d_HI_i_const_next[i])
                .max(self.d_HD_i_const_next[i]);

            self.d_HS_j_const_next[i] = self.d_alp_data.d_smatr
                [self.d_seqi[len2 as usize] as usize][self.d_seqj[len1 as usize] as usize]
                + self.d_H_j_const_pred[i];
            self.d_HI_j_const_next[i] = (self.d_HS_ij_pred - self.d_alp_data.d_open2)
                .max(self.d_HI_ij_pred - self.d_alp_data.d_epen2)
                .max(self.d_HD_ij_pred - self.d_alp_data.d_open2);
            self.d_HD_j_const_next[i] = (self.d_HS_j_const_next[1] - self.d_alp_data.d_open1)
                .max(self.d_HD_j_const_next[1] - self.d_alp_data.d_epen1);
            self.d_H_j_const_next[i] = self.d_HS_j_const_next[i]
                .max(self.d_HI_j_const_next[i])
                .max(self.d_HD_j_const_next[i]);
        }
        self.d_HS_ij_next = self.d_alp_data.d_smatr[self.d_seqi[len1 as usize] as usize]
            [self.d_seqj[len1 as usize] as usize]
            + self.d_H_ij_pred;
        self.d_HI_ij_next = (self.d_HS_i_const_next[0] - self.d_alp_data.d_open2)
            .max(self.d_HI_i_const_next[0] - self.d_alp_data.d_epen2)
            .max(self.d_HD_i_const_next[0] - self.d_alp_data.d_open2);
        self.d_HD_ij_next = (self.d_HS_j_const_next[0] - self.d_alp_data.d_open1)
            .max(self.d_HD_j_const_next[0] - self.d_alp_data.d_epen1);
        self.d_H_ij_next = self
            .d_HS_ij_next
            .max(self.d_HI_ij_next)
            .max(self.d_HD_ij_next);
        self.d_cells_counts.increase_elem_by_1(self.d_H_ij_next);
        if sentinel_i_boundary >= 1 {
            for i in 0..=sentinel_i_boundary - 1 {
                self.d_cells_counts
                    .increase_elem_by_1(self.d_H_i_const_next[i as usize]);
            }
        }
        if sentinel_j_boundary >= 1 {
            for i in 0..=sentinel_j_boundary - 1 {
                self.d_cells_counts
                    .increase_elem_by_1(self.d_H_j_const_next[i as usize]);
            }
        }
        let mut tmp = self.d_H_ij_next;
        if sentinel_i_boundary >= 1 {
            for i in 0..=sentinel_i_boundary - 1 {
                tmp = tmp.max(self.d_H_i_const_next[i as usize]);
            }
        }
        if sentinel_j_boundary >= 1 {
            for i in 0..=sentinel_j_boundary - 1 {
                tmp = tmp.max(self.d_H_j_const_next[i as usize]);
            }
        }
        self.d_H_edge_max[self.d_H_matr_len as usize] = tmp;
        self.d_M = tmp.max(self.d_M);
        let level = tmp - diff_opt_;
        self.d_sentinel_i_next = 1;
        self.d_sentinel_j_next = 1;
        if sentinel_i_boundary > 1 {
            for i in (1..=sentinel_i_boundary - 1).rev() {
                if self.d_H_i_const_next[i as usize] >= level {
                    self.d_sentinel_i_next = i;
                    break;
                }
            }
        }
        if sentinel_j_boundary > 1 {
            for i in (1..=sentinel_j_boundary - 1).rev() {
                if self.d_H_j_const_next[i as usize] >= level {
                    self.d_sentinel_j_next = i;
                    break;
                }
            }
        }
        if self.d_is_now && tmp > self.d_alp.d_elem[self.d_nalp as usize] {
            self.d_nalp += 1;
            self.d_alp.set_elem(self.d_nalp, tmp);
            self.d_alp_pos.set_elem(self.d_nalp, self.d_H_matr_len);
            let mut I = -1;
            let mut J = -1;
            if sentinel_i_boundary >= 1 {
                for i in 0..=sentinel_i_boundary - 1 {
                    if tmp == self.d_H_i_const_next[i as usize] {
                        I = i;
                    }
                }
            }
            if sentinel_j_boundary >= 1 {
                for i in 0..=sentinel_j_boundary - 1 {
                    if tmp == self.d_H_j_const_next[i as usize] {
                        J = i;
                    }
                }
            }
            self.d_H_I.set_elem(self.d_nalp, self.d_H_matr_len - I - 1);
            self.d_H_J.set_elem(self.d_nalp, self.d_H_matr_len - J - 1);
            let st = self.save_state()?;
            self.d_alp_states.resize(self.d_nalp as usize + 1, None);
            self.d_alp_states[self.d_nalp as usize] = Some(st);
        }
        self.check_time_function()
    }

    pub fn check_time_function(&mut self) -> Result<(), Error> {
        if self.d_check_time_flag {
            let mut time_after3 = 0.0;
            get_current_time(&mut time_after3);
            if (time_after3 - self.d_alp_data.d_time_before1) > self.d_alp_data.d_max_time {
                if self.d_time_error_flag {
                    return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
                }
                self.d_time_limit_flag = true;
                if self.d_single_realiztion_calculation_flag {
                    return Err(Error::new(String::new(), 0));
                }
                return Ok(());
            }
        }
        if self.d_alp_data.d_max_time <= 0.0
            && self.d_alp_data.d_max_time_with_computation_parameters > 0.0
        {
            let mut time_after3 = 0.0;
            get_current_time(&mut time_after3);
            if (time_after3 - self.d_alp_data.d_time_before1)
                > self.d_alp_data.d_max_time_with_computation_parameters
            {
                return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
            }
        }
        Ok(())
    }

    pub fn restore_state(&mut self, state_: &state) -> Result<(), Error> {
        self.d_M = state_.d_M;
        self.d_H_matr_len = state_.d_H_matr_len;
        if self.d_H_matr_len < 0 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        self.d_is_now = false;
        self.d_cells_counts = array::new(None);
        self.d_cells_counts.set_elems(&state_.d_cells_counts);
        self.d_HS_ij_next = state_.d_HS_ij_next;
        self.d_HI_ij_next = state_.d_HI_ij_next;
        self.d_HD_ij_next = state_.d_HD_ij_next;
        self.d_H_ij_next = state_.d_H_ij_next;
        for i in 0..self.d_H_matr_len as usize {
            self.d_HS_i_const_next[i] = state_.d_HS_i_const_next[i];
            self.d_HI_i_const_next[i] = state_.d_HI_i_const_next[i];
            self.d_HD_i_const_next[i] = state_.d_HD_i_const_next[i];
            self.d_H_i_const_next[i] = state_.d_H_i_const_next[i];
            self.d_HS_j_const_next[i] = state_.d_HS_j_const_next[i];
            self.d_HI_j_const_next[i] = state_.d_HI_j_const_next[i];
            self.d_HD_j_const_next[i] = state_.d_HD_j_const_next[i];
            self.d_H_j_const_next[i] = state_.d_H_j_const_next[i];
        }
        self.d_sentinel_i_next = state_.d_sentinel_i_next;
        self.d_sentinel_j_next = state_.d_sentinel_j_next;
        Ok(())
    }

    pub fn save_state(&mut self) -> Result<state, Error> {
        if self.d_H_matr_len < 0 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        let mut state_ = state::new();
        self.d_alp_data.d_memory_size_in_MB += std::mem::size_of::<state>() as f64 / mb_bytes;
        state_.d_M = self.d_M;
        state_.d_cells_counts = array::new(None);
        state_.d_cells_counts.set_elems(&self.d_cells_counts);
        state_.d_H_matr_len = self.d_H_matr_len;
        state_.d_HS_ij_next = self.d_HS_ij_next;
        state_.d_HI_ij_next = self.d_HI_ij_next;
        state_.d_HD_ij_next = self.d_HD_ij_next;
        state_.d_H_ij_next = self.d_H_ij_next;
        if self.d_H_matr_len > 0 {
            let n = self.d_H_matr_len as usize;
            state_.d_HS_i_const_next = self.d_HS_i_const_next[..n].to_vec();
            state_.d_HI_i_const_next = self.d_HI_i_const_next[..n].to_vec();
            state_.d_HD_i_const_next = self.d_HD_i_const_next[..n].to_vec();
            state_.d_H_i_const_next = self.d_H_i_const_next[..n].to_vec();
            state_.d_HS_j_const_next = self.d_HS_j_const_next[..n].to_vec();
            state_.d_HI_j_const_next = self.d_HI_j_const_next[..n].to_vec();
            state_.d_HD_j_const_next = self.d_HD_j_const_next[..n].to_vec();
            state_.d_H_j_const_next = self.d_H_j_const_next[..n].to_vec();
            self.d_alp_data.d_memory_size_in_MB +=
                8.0 * (n * std::mem::size_of::<i64>()) as f64 / mb_bytes;
        }
        state_.d_sentinel_i_next = self.d_sentinel_i_next;
        state_.d_sentinel_j_next = self.d_sentinel_j_next;
        Ok(state_)
    }

    pub fn John2_weight_calculation(&mut self, length_: i64) -> Result<f64, Error> {
        if length_ == 0 {
            return Ok(1.0);
        }
        if self.d_W_matr_len > length_ {
            return Err(Error::new(
                "Error - unexpected parameter in alp::John2_weight_calculation\n".to_string(),
                4,
            ));
        }
        while self.d_W_matr_len < length_ {
            self.increment_W_weights()?;
        }
        let is = self.d_alp_data.d_is.as_ref().ok_or_else(|| {
            Error::new(
                "Error - unexpected parameter in alp::John2_weight_calculation\n".to_string(),
                4,
            )
        })?;
        let len1 = self.d_W_matr_len - 1;
        let mut US = 0.0;
        let mut UD = 0.0;
        let mut UI = self.d_WI_j_const_next[len1 as usize] / (1.0 - is.d_nu);
        let mut VS = 0.0;
        let mut VI = 0.0;
        let mut VD = self.d_WD_i_const_next[len1 as usize] / (1.0 - is.d_nu);
        for j in 1..=length_ - 1 {
            let US_next = self.d_alp_data.d_r_i_dot[self.d_seqi[(j - 1) as usize] as usize]
                * (is.d_eta * US + is.d_mu_SI * UI + is.d_mu_SD * UD)
                + self.d_WS_j_const_next[(len1 - j) as usize];
            let UD_next = is.d_mu_DS * US + is.d_nu * UD;
            let UI_next = (is.d_mu_IS * US_next
                + is.d_mu_ID * UD_next
                + self.d_WI_j_const_next[(len1 - j) as usize])
                / (1.0 - is.d_nu);
            let VS_next = self.d_alp_data.d_r_dot_j[self.d_seqj[(j - 1) as usize] as usize]
                * (is.d_eta * VS + is.d_mu_SI * VI + is.d_mu_SD * VD)
                + self.d_WS_i_const_next[(len1 - j) as usize];
            let VI_next = is.d_mu_IS * VS + is.d_mu_ID * VD + is.d_nu * VI;
            let VD_next = (is.d_mu_DS * VS_next + self.d_WD_i_const_next[(len1 - j) as usize])
                / (1.0 - is.d_nu);
            US = US_next;
            UD = UD_next;
            UI = UI_next;
            VS = VS_next;
            VD = VD_next;
            VI = VI_next;
        }
        let j = length_;
        let US_next = self.d_alp_data.d_r_i_dot[self.d_seqi[(j - 1) as usize] as usize]
            * (is.d_eta * US + is.d_mu_SI * UI + is.d_mu_SD * UD)
            + self.d_WS_ij_next;
        let UD_next = is.d_mu_DS * US + is.d_nu * UD;
        let _UI_next =
            (is.d_mu_IS * US_next + is.d_mu_ID * UD_next + self.d_WI_ij_next) / (1.0 - is.d_nu);
        let VS_next = self.d_alp_data.d_r_dot_j[self.d_seqj[(j - 1) as usize] as usize]
            * (is.d_eta * VS + is.d_mu_SI * VI + is.d_mu_SD * VD)
            + self.d_WS_ij_next;
        let VI_next = is.d_mu_IS * VS + is.d_mu_ID * VD + is.d_nu * VI;
        US = US_next;
        UD = UD_next;
        VS = VS_next;
        VI = VI_next;
        let weight = -self.d_WS_ij_next + US + UD + VS + VI;
        if weight == 0.0 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        Ok(1.0 / weight)
    }

    pub fn simulate_next_alp(&mut self) -> Result<(), Error> {
        if !self.d_success {
            return Ok(());
        }
        if !self.d_is_now {
            return Err(Error::new(
                "Unexpected error - ALP can be generated only in the importance sampling mode\n"
                    .to_string(),
                4,
            ));
        }
        let target_nalp = self.d_nalp + 1;
        while self.d_nalp < target_nalp {
            let k = self.d_seqi_len.min(self.d_seqj_len);
            while self.d_seqi_len.min(self.d_seqj_len) != k + 1 {
                let mut rng = std::iter::from_fn(|| Some(crate::stats::pvalues::pvalues::ran3()));
                let success = self.one_step_of_importance_sampling_without_weight_calculation(
                    self.d_alp_data.d_dim1_tmp,
                    self.d_alp_data.d_dim2_tmp,
                    &mut rng,
                )?;
                self.check_time_function()?;
                if !success {
                    self.d_success = false;
                    return Ok(());
                }
            }
            if self.d_sentinels_flag {
                self.increment_H_weights_with_sentinels(self.d_diff_opt)?;
            } else {
                self.increment_H_weights()?;
            }
            if self.d_time_limit_flag {
                self.d_success = false;
                return Ok(());
            }
            self.increment_W_weights()?;
        }
        let weight = self.John2_weight_calculation(self.d_seqi_len.min(self.d_seqj_len))?;
        if weight <= 0.0 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        self.d_alp_weights.set_elem(self.d_nalp, weight);
        Ok(())
    }

    pub fn simulate_alp_upto_the_given_number(&mut self, nalp_: i64) -> Result<(), Error> {
        self.d_sentinels_flag = false;
        while self.d_nalp < nalp_ {
            self.simulate_next_alp()?;
            if !self.d_success {
                return Ok(());
            }
        }
        Ok(())
    }

    pub fn simulate_alp_upto_the_given_level(&mut self, M_min_: i64) -> Result<(), Error> {
        self.d_sentinels_flag = false;
        while self.d_alp.d_elem[self.d_nalp as usize] < M_min_ {
            self.simulate_next_alp()?;
            if !self.d_success {
                return Ok(());
            }
        }
        self.d_nalp_killing = self.d_nalp;
        Ok(())
    }

    pub fn kill_upto_level(
        &mut self,
        M_min_: i64,
        M_level_: i64,
        M_upper_level_: Option<&i64>,
    ) -> Result<(), Error> {
        if self.d_is_now {
            while self.d_alp.d_elem[self.d_nalp as usize] < M_min_ {
                self.simulate_next_alp()?;
                if !self.d_success {
                    return Ok(());
                }
            }
            self.d_is_now = false;
            self.d_nalp_killing = -1;
            for i in 0..=self.d_nalp {
                if self.d_alp.d_elem[i as usize] >= M_min_ {
                    self.d_nalp_killing = i;
                    break;
                }
            }
            if self.d_nalp_killing == -1 {
                return Err(Error::new("Unexpected error\n".to_string(), 4));
            }
            let state_ = self
                .d_alp_states
                .get(self.d_nalp_killing as usize)
                .and_then(|x| x.as_ref())
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?
                .clone();
            self.restore_state(&state_)?;
        }

        while self.d_H_edge_max[self.d_H_matr_len as usize] >= M_level_ {
            if self.d_H_matr_len + 1 >= self.d_alp_data.d_dim1_tmp {
                self.d_success = false;
                return Ok(());
            }
            if let Some(M_upper_level_) = M_upper_level_ {
                if self.d_H_edge_max[self.d_H_matr_len as usize] > *M_upper_level_ {
                    self.d_success = false;
                    return Ok(());
                }
            }
            if self.d_H_matr_len + 1 > self.d_seq_a_len {
                self.increment_sequences();
            }
            self.d_seqi_len = self.d_H_matr_len + 1;
            self.d_seqj_len = self.d_H_matr_len + 1;
            self.d_seqi[(self.d_seqi_len - 1) as usize] = self.random_AA1()?;
            self.d_seqj[(self.d_seqj_len - 1) as usize] = self.random_AA2()?;

            if self.d_sentinels_flag {
                self.increment_H_weights_with_sentinels(self.d_diff_opt)?;
            } else {
                self.increment_H_weights()?;
            }
            if self.d_time_limit_flag {
                self.d_success = false;
                return Ok(());
            }
        }
        self.d_success = true;
        Ok(())
    }

    pub fn partially_release_memory(&mut self) {
        self.d_seqi.clear();
        self.d_seqj.clear();
        self.d_WS_i_const_pred.clear();
        self.d_WI_i_const_pred.clear();
        self.d_WD_i_const_pred.clear();
        self.d_WS_i_const_next.clear();
        self.d_WI_i_const_next.clear();
        self.d_WD_i_const_next.clear();
        self.d_WS_j_const_pred.clear();
        self.d_WI_j_const_pred.clear();
        self.d_WD_j_const_pred.clear();
        self.d_WS_j_const_next.clear();
        self.d_WI_j_const_next.clear();
        self.d_WD_j_const_next.clear();
        self.d_HS_i_const_pred.clear();
        self.d_HI_i_const_pred.clear();
        self.d_HD_i_const_pred.clear();
        self.d_H_i_const_pred.clear();
        self.d_HS_i_const_next.clear();
        self.d_HI_i_const_next.clear();
        self.d_HD_i_const_next.clear();
        self.d_H_i_const_next.clear();
        self.d_HS_j_const_pred.clear();
        self.d_HI_j_const_pred.clear();
        self.d_HD_j_const_pred.clear();
        self.d_H_j_const_pred.clear();
        self.d_HS_j_const_next.clear();
        self.d_HI_j_const_next.clear();
        self.d_HD_j_const_next.clear();
        self.d_H_j_const_next.clear();
        self.d_H_edge_max.clear();
        for state in &mut self.d_alp_states {
            if let Some(state) = state {
                state.d_HS_i_const_next.clear();
                state.d_HI_i_const_next.clear();
                state.d_HD_i_const_next.clear();
                state.d_H_i_const_next.clear();
                state.d_HS_j_const_next.clear();
                state.d_HI_j_const_next.clear();
                state.d_HD_j_const_next.clear();
                state.d_H_j_const_next.clear();
                state.d_cells_counts = array::new(None);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_data() -> alp_data {
        let smatr = vec![vec![-2, 1], vec![1, -2]];
        let rr = vec![0.5, 0.5];
        alp_data::new_from_parameters(
            123, None, 11, 11, 11, 1, 1, 1, 2, &smatr, &rr, &rr, 1.0, 10.0, 128.0, 0.05, 0.05,
            false, -1.0, -1.0,
        )
        .unwrap()
    }

    #[test]
    fn test_degree_matches_cpp_cases() {
        assert_eq!(alp::degree(0.0, 0.0).unwrap(), 1.0);
        assert_eq!(alp::degree(0.0, 2.0).unwrap(), 0.0);
        assert!((alp::degree(2.0, 3.0).unwrap() - 8.0).abs() < 1e-14);
        assert!(alp::degree(-1.0, 2.0).is_err());
    }

    #[test]
    fn test_constructor_initializes_base_alp_and_state() {
        let a = alp::new(test_data()).unwrap();
        assert_eq!(a.d_W_matr_len, 0);
        assert_eq!(a.d_H_matr_len, 0);
        assert_eq!(a.d_nalp, 0);
        assert_eq!(a.d_alp.d_elem[0], 0);
        assert_eq!(a.d_alp_weights.d_elem[0], 1.0);
        assert_eq!(
            a.d_cells_counts.d_elem[(0 - a.d_cells_counts.d_ind0) as usize],
            1
        );
        assert_eq!(a.d_alp_states.len(), 1);
        assert!(a.d_alp_states[0].is_some());
    }

    #[test]
    fn test_sequence_and_matrix_growth() {
        let mut a = alp::new(test_data()).unwrap();
        a.increment_sequences();
        assert_eq!(a.d_seq_a_len, 30);
        assert_eq!(a.d_seqi.len(), 30);
        a.increment_W_matrix();
        assert_eq!(a.d_W_matr_a_len, 30);
        assert_eq!(a.d_WS_i_const_next.len(), 30);
        a.increment_H_matrix();
        assert_eq!(a.d_H_matr_a_len, 30);
        assert_eq!(a.d_HS_i_const_next.len(), 30);
        assert_eq!(a.d_H_edge_max.len(), 31);
    }

    #[test]
    fn test_one_step_importance_sampling_deterministic_stream() {
        let mut a = alp::new(test_data()).unwrap();
        let mut values = [0.0, 0.1, 0.8].into_iter();
        assert!(a
            .one_step_of_importance_sampling_without_weight_calculation(10, 10, &mut values)
            .unwrap());
        assert_eq!(a.d_IS_state, b'S');
        assert_eq!(a.d_seqi_len, 1);
        assert_eq!(a.d_seqj_len, 1);
    }

    #[test]
    fn test_save_restore_state_roundtrip() {
        let mut a = alp::new(test_data()).unwrap();
        a.increment_H_matrix();
        a.d_H_matr_len = 1;
        a.d_M = 7;
        a.d_HS_ij_next = 1;
        a.d_HI_ij_next = 2;
        a.d_HD_ij_next = 3;
        a.d_H_ij_next = 4;
        a.d_HS_i_const_next[0] = 11;
        a.d_HI_i_const_next[0] = 12;
        a.d_HD_i_const_next[0] = 13;
        a.d_H_i_const_next[0] = 14;
        a.d_HS_j_const_next[0] = 21;
        a.d_HI_j_const_next[0] = 22;
        a.d_HD_j_const_next[0] = 23;
        a.d_H_j_const_next[0] = 24;
        let st = a.save_state().unwrap();
        a.d_M = 0;
        a.d_HS_i_const_next[0] = 0;
        a.restore_state(&st).unwrap();
        assert_eq!(a.d_M, 7);
        assert_eq!(a.d_HS_i_const_next[0], 11);
        assert_eq!(a.d_H_j_const_next[0], 24);
    }

    #[test]
    fn test_increment_h_weights_without_insertions_known_sequence() {
        let mut a = alp::new(test_data()).unwrap();
        a.increment_sequences();
        a.d_seqi[0] = 0;
        a.d_seqj[0] = 1;
        a.d_seqi_len = 1;
        a.d_seqj_len = 1;
        a.increment_H_weights_without_insertions_after_deletions()
            .unwrap();
        assert_eq!(a.d_H_matr_len, 1);
        assert_eq!(a.d_HS_ij_next, 1);
        assert_eq!(a.d_H_ij_next, 1);
        assert_eq!(a.d_M, 1);
        assert_eq!(a.d_nalp, 1);

        a.d_seqi[1] = 1;
        a.d_seqj[1] = 0;
        a.d_seqi_len = 2;
        a.d_seqj_len = 2;
        a.increment_H_weights_without_insertions_after_deletions()
            .unwrap();
        assert_eq!(a.d_H_matr_len, 2);
        assert_eq!(a.d_HS_ij_next, 2);
        assert_eq!(a.d_H_ij_next, 2);
        assert_eq!(a.d_M, 2);
        assert_eq!(a.d_alp.d_elem[a.d_nalp as usize], 2);
    }

    #[test]
    fn test_increment_h_weights_with_insertions_known_sequence() {
        let mut data = test_data();
        data.d_insertions_after_deletions = true;
        let mut a = alp::new(data).unwrap();
        a.increment_sequences();
        a.d_seqi[0] = 0;
        a.d_seqj[0] = 1;
        a.d_seqi_len = 1;
        a.d_seqj_len = 1;
        a.increment_H_weights_with_insertions_after_deletions()
            .unwrap();
        assert_eq!(a.d_H_ij_next, 1);

        a.d_seqi[1] = 1;
        a.d_seqj[1] = 0;
        a.d_seqi_len = 2;
        a.d_seqj_len = 2;
        a.increment_H_weights_with_insertions_after_deletions()
            .unwrap();
        assert_eq!(a.d_H_ij_next, 2);
        assert_eq!(a.d_M, 2);
    }

    #[test]
    fn test_increment_h_weights_with_sentinels_known_sequence() {
        let mut a = alp::new(test_data()).unwrap();
        a.increment_sequences();
        a.d_seqi[0] = 0;
        a.d_seqj[0] = 1;
        a.d_seqi_len = 1;
        a.d_seqj_len = 1;
        a.increment_H_weights_with_sentinels_without_insertions_after_deletions(99)
            .unwrap();
        assert_eq!(a.d_H_ij_next, 1);
        assert_eq!(a.d_sentinel_i_next, 1);
        assert_eq!(a.d_sentinel_j_next, 1);

        a.d_seqi[1] = 1;
        a.d_seqj[1] = 0;
        a.d_seqi_len = 2;
        a.d_seqj_len = 2;
        a.increment_H_weights_with_sentinels_without_insertions_after_deletions(99)
            .unwrap();
        assert_eq!(a.d_H_ij_next, 2);
        assert_eq!(a.d_M, 2);
    }

    #[test]
    fn test_kill_upto_level_restores_saved_state() {
        let mut a = alp::new(test_data()).unwrap();
        let mut st = a.d_alp_states[0].clone().unwrap();
        st.d_M = 9;
        st.d_H_ij_next = 7;
        a.d_alp_states[0] = Some(st);
        a.kill_upto_level(0, 99, None).unwrap();
        assert!(!a.d_is_now);
        assert_eq!(a.d_nalp_killing, 0);
        assert_eq!(a.d_M, 9);
        assert_eq!(a.d_H_ij_next, 7);
        assert!(a.d_success);
    }

    #[test]
    fn test_kill_upto_level_upper_level_failure() {
        let mut a = alp::new(test_data()).unwrap();
        a.d_is_now = false;
        a.d_H_matr_len = 0;
        a.d_H_edge_max[0] = 5;
        a.kill_upto_level(0, 1, Some(&4)).unwrap();
        assert!(!a.d_success);
    }
}
