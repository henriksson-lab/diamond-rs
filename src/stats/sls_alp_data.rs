#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use std::fs;
use std::ops::AddAssign;
use std::str::FromStr;

use crate::stats::sls_alp_regression::alp_reg;
use crate::stats::sls_basic::{random_seed_from_time, round, Error, Tmax, Tmax4, Tmin};

pub const mb_bytes: f64 = 1048576.0;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct struct_for_randomization {
    pub d_random_seed: i64,
    pub d_first_stage_preliminary_realizations_numbers_ALP: Vec<i64>,
    pub d_preliminary_realizations_numbers_ALP: Vec<i64>,
    pub d_preliminary_realizations_numbers_killing: Vec<i64>,
    pub d_total_realizations_number_with_ALP: i64,
    pub d_total_realizations_number_with_killing: i64,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct error_for_single_realization {
    pub st: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct data_for_lambda_equation {
    pub d_number_of_AA: i64,
    pub d_smatr: Vec<Vec<i64>>,
    pub d_RR1: Vec<f64>,
    pub d_RR2: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct alp_data {
    pub d_open: i64,
    pub d_open1: i64,
    pub d_open2: i64,
    pub d_epen: i64,
    pub d_epen1: i64,
    pub d_epen2: i64,
    pub d_max_time: f64,
    pub d_max_time_for_quick_tests: f64,
    pub d_max_time_with_computation_parameters: f64,
    pub d_max_mem: f64,
    pub d_eps_lambda: f64,
    pub d_eps_K: f64,
    pub d_insertions_after_deletions: bool,
    pub d_smatr_symmetric_flag: bool,
    pub d_number_of_AA: i64,
    pub d_number_of_AA_smatr: i64,
    pub d_smatr: Vec<Vec<i64>>,
    pub d_RR1: Vec<f64>,
    pub d_RR1_sum: Vec<f64>,
    pub d_RR1_sum_elements: Vec<i64>,
    pub d_RR2: Vec<f64>,
    pub d_RR2_sum: Vec<f64>,
    pub d_RR2_sum_elements: Vec<i64>,
    pub d_random_seed: i64,
    pub d_randout: String,
    pub d_memory_size_in_MB: f64,
    pub d_is: Option<Box<importance_sampling>>,
    pub d_r_i_dot: Vec<f64>,
    pub d_r_dot_j: Vec<f64>,
    pub d_minimum_realizations_number: i64,
    pub d_sentinels_flag: bool,
    pub d_dim1_tmp: i64,
    pub d_dim2_tmp: i64,
    pub d_realizations_number2: i64,
    pub d_time_before1: f64,
    pub d_rand_all: Option<Box<struct_for_randomization>>,
    pub d_rand_flag: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct array_positive<T> {
    pub d_step: i64,
    pub d_dim: i64,
    pub d_elem: Vec<T>,
}

impl<T> Default for array_positive<T> {
    fn default() -> Self {
        Self {
            d_step: 10,
            d_dim: -1,
            d_elem: Vec::new(),
        }
    }
}

impl<T> array_positive<T>
where
    T: Copy + Default + AddAssign + From<u8>,
{
    pub fn new(_alp_data_: Option<&mut alp_data>) -> Result<Self, Error> {
        Ok(Self::default())
    }

    pub fn increment_array(&mut self, ind_: i64) {
        let o_dim = self.d_dim;
        while ind_ > self.d_dim {
            self.d_dim += self.d_step;
        }
        let old_len = if o_dim < 0 { 0 } else { (o_dim + 1) as usize };
        self.d_elem.resize((self.d_dim + 1) as usize, T::default());
        for i in old_len..self.d_elem.len() {
            self.d_elem[i] = T::from(0u8);
        }
    }

    pub fn set_elem(&mut self, ind_: i64, elem_: T) {
        if ind_ > self.d_dim {
            self.increment_array(ind_);
        }
        self.d_elem[ind_ as usize] = elem_;
    }

    pub fn increase_elem_by_1(&mut self, ind_: i64) {
        if ind_ > self.d_dim {
            self.increment_array(ind_);
        }
        self.d_elem[ind_ as usize] += T::from(1u8);
    }

    pub fn increase_elem_by_x(&mut self, ind_: i64, x_: T) {
        if ind_ > self.d_dim {
            self.increment_array(ind_);
        }
        self.d_elem[ind_ as usize] += x_;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct array<T> {
    pub d_step: i64,
    pub d_dim: i64,
    pub d_ind0: i64,
    pub d_dim_plus_d_ind0: i64,
    pub d_elem: Vec<T>,
}

impl<T> Default for array<T> {
    fn default() -> Self {
        Self {
            d_step: 10,
            d_dim: -1,
            d_ind0: 0,
            d_dim_plus_d_ind0: -1,
            d_elem: Vec::new(),
        }
    }
}

impl<T> array<T>
where
    T: Copy + Default + AddAssign + From<u8>,
{
    pub fn new(_alp_data_: Option<&mut alp_data>) -> Self {
        Self::default()
    }

    pub fn increment_array_on_the_right(&mut self, ind_: i64) {
        let o_dim = self.d_dim;
        while ind_ > self.d_dim_plus_d_ind0 {
            self.d_dim += self.d_step;
            self.d_dim_plus_d_ind0 += self.d_step;
        }
        let old_len = if o_dim < 0 { 0 } else { (o_dim + 1) as usize };
        self.d_elem.resize((self.d_dim + 1) as usize, T::default());
        for i in old_len..self.d_elem.len() {
            self.d_elem[i] = T::from(0u8);
        }
    }

    pub fn increment_array_on_the_left(&mut self, ind_: i64) {
        let o_dim = self.d_dim;
        while ind_ < self.d_ind0 {
            self.d_dim += self.d_step;
            self.d_ind0 -= self.d_step;
        }
        let jump = self.d_dim - o_dim;
        if jump > 0 {
            let mut d_elem_new = vec![T::from(0u8); (self.d_dim + 1) as usize];
            for i in 0..=o_dim {
                if i >= 0 {
                    d_elem_new[(i + jump) as usize] = self.d_elem[i as usize];
                }
            }
            self.d_elem = d_elem_new;
            self.d_dim_plus_d_ind0 = self.d_dim + self.d_ind0;
        }
    }

    pub fn set_elems(&mut self, a_: &array<T>) {
        let a0 = a_.d_ind0;
        let a1 = a_.d_dim_plus_d_ind0;
        if a0 > a1 {
            return;
        }
        while a1 > self.d_dim_plus_d_ind0 {
            self.d_dim_plus_d_ind0 += self.d_step;
        }
        while a0 < self.d_ind0 {
            self.d_ind0 -= self.d_step;
        }
        self.d_dim = self.d_dim_plus_d_ind0 - self.d_ind0;
        self.d_elem = vec![T::from(0u8); (self.d_dim + 1) as usize];
        for i in a0..=a1 {
            self.d_elem[(i - self.d_ind0) as usize] = a_.d_elem[(i - a0) as usize];
        }
    }

    pub fn set_elem(&mut self, ind_: i64, elem_: T) {
        if ind_ > self.d_dim_plus_d_ind0 {
            self.increment_array_on_the_right(ind_);
        }
        if ind_ < self.d_ind0 {
            self.increment_array_on_the_left(ind_);
        }
        self.d_elem[(ind_ - self.d_ind0) as usize] = elem_;
    }

    pub fn increase_elem_by_1(&mut self, ind_: i64) {
        if ind_ > self.d_dim_plus_d_ind0 {
            self.increment_array_on_the_right(ind_);
        }
        if ind_ < self.d_ind0 {
            self.increment_array_on_the_left(ind_);
        }
        self.d_elem[(ind_ - self.d_ind0) as usize] += T::from(1u8);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct q_elem {
    pub d_a: i64,
    pub d_b: i64,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct importance_sampling {
    pub d_mu: f64,
    pub d_nu: f64,
    pub d_eta: f64,
    pub d_mu_SI: f64,
    pub d_mu_DS: f64,
    pub d_mu_ID: f64,
    pub d_mu_IS: f64,
    pub d_mu_SD: f64,
    pub d_elements: Vec<q_elem>,
    pub d_elements_values: Vec<f64>,
    pub d_for_D: [f64; 3],
    pub d_for_I: [f64; 2],
    pub d_for_S: [f64; 3],
    pub d_for_D_states: [u8; 3],
    pub d_for_I_states: [u8; 2],
    pub d_for_S_states: [u8; 3],
    pub d_exp_s: Vec<Vec<f64>>,
    pub d_lambda: f64,
    pub d_ungap_lambda: f64,
    pub d_is_number_of_AA: i64,
}

impl importance_sampling {
    pub fn new(
        alp_data_: Option<&mut alp_data>,
        open_: i64,
        epen_: i64,
        temperature_: f64,
        number_of_AA_: i64,
        smatr_: &[Vec<i64>],
        RR1_: &[f64],
        RR2_: &[f64],
    ) -> Result<Self, Error> {
        let mut result = Self::default();
        result.d_is_number_of_AA = number_of_AA_;

        let data = data_for_lambda_equation {
            d_number_of_AA: number_of_AA_,
            d_smatr: smatr_.to_vec(),
            d_RR1: RR1_.to_vec(),
            d_RR2: RR2_.to_vec(),
        };

        let mut smatr_max = smatr_[0][0];
        let mut smatr_max_i = 0usize;
        let mut smatr_max_j = 0usize;
        let mut smatr_min = smatr_[0][0];
        let mut smatr_pos_max = i64::MIN;
        let mut smatr_neg_min = i64::MAX;
        let mut eps = 0.00001;
        let threshold = f64::MIN_POSITIVE * 10.0;
        let mut aver_score = 0.0;
        for i in 0..number_of_AA_ as usize {
            for j in 0..number_of_AA_ as usize {
                if RR1_[i] * RR2_[j] <= threshold {
                    continue;
                }
                aver_score += RR1_[i] * RR2_[j] * smatr_[i][j] as f64;
                if smatr_max < smatr_[i][j] {
                    smatr_max = smatr_[i][j];
                    smatr_max_i = i;
                    smatr_max_j = j;
                }
                smatr_min = Tmin(smatr_min, smatr_[i][j]);
                if smatr_[i][j] > 0 {
                    smatr_pos_max = Tmax(smatr_pos_max, smatr_[i][j]);
                }
                if smatr_[i][j] < 0 {
                    smatr_neg_min = Tmin(smatr_neg_min, smatr_[i][j]);
                }
            }
        }
        let _ = (smatr_min, smatr_pos_max, smatr_neg_min);

        if aver_score >= -threshold {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }
        if smatr_max <= 0 {
            return Err(Error::new(
                "Error - at least one element of the scoring matrix must be positive\n".to_string(),
                3,
            ));
        }

        let mut a = eps;
        while Self::lambda_equation(a, &data) > 0.0 {
            a /= 2.0;
            if a < threshold * 100.0 {
                return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
            }
        }
        if a < threshold * 100.0 {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }
        eps = a / 10.0;

        let tmp_pr = RR1_[smatr_max_i] * RR2_[smatr_max_j];
        let b = ((1.0 + 10.0 * eps).ln() - tmp_pr.ln()) / smatr_max as f64;
        let mut res_lambda = Vec::new();
        alp_reg::find_tetta_general(
            |x| Self::lambda_equation(x, &data),
            a,
            b,
            2,
            eps,
            &mut res_lambda,
        )?;
        res_lambda.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if res_lambda.is_empty() {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        result.d_lambda = res_lambda[res_lambda.len() - 1];
        result.d_ungap_lambda = result.d_lambda;
        result.d_lambda *= temperature_;

        let n2 = (number_of_AA_ * number_of_AA_) as usize;
        result.d_elements = vec![q_elem::default(); n2];
        result.d_elements_values = vec![0.0; n2];
        result.d_exp_s = vec![vec![0.0; number_of_AA_ as usize]; number_of_AA_ as usize];

        let mut ind = 0usize;
        let mut sum = 0.0;
        for a_idx in 0..number_of_AA_ as usize {
            for b_idx in 0..number_of_AA_ as usize {
                result.d_exp_s[a_idx][b_idx] =
                    (result.d_lambda * smatr_[a_idx][b_idx] as f64).exp();
                result.d_elements_values[ind] =
                    RR1_[a_idx] * RR2_[b_idx] * result.d_exp_s[a_idx][b_idx];
                sum += result.d_elements_values[ind];
                ind += 1;
            }
        }

        for a_idx in 0..number_of_AA_ as usize {
            for b_idx in 0..number_of_AA_ as usize {
                result.d_exp_s[a_idx][b_idx] /= sum;
            }
        }
        for ind in 0..n2 {
            result.d_elements_values[ind] /= sum;
        }
        for ind in 1..n2 {
            result.d_elements_values[ind] =
                result.d_elements_values[ind - 1] + result.d_elements_values[ind];
        }

        ind = 0;
        for a_idx in 0..number_of_AA_ {
            for b_idx in 0..number_of_AA_ {
                result.d_elements[ind] = q_elem {
                    d_a: a_idx,
                    d_b: b_idx,
                };
                ind += 1;
            }
        }

        result.d_mu = (-result.d_lambda.abs() * open_ as f64).exp();
        result.d_nu = (-result.d_lambda.abs() * epen_ as f64).exp();
        let tmp = 1.0 + result.d_mu - result.d_nu;
        result.d_eta = (1.0 - result.d_nu) * (1.0 - result.d_nu) / (tmp * tmp);
        result.d_mu_SI = 1.0 - result.d_nu;
        result.d_mu_IS = result.d_mu * (1.0 - result.d_nu) / (tmp * tmp);
        result.d_mu_DS = result.d_mu / tmp;
        result.d_mu_SD = (1.0 - result.d_nu) * (1.0 - result.d_nu) / tmp;
        result.d_mu_ID = result.d_mu * (1.0 - result.d_nu) / tmp;

        result.d_for_D = [
            result.d_nu,
            result.d_nu + result.d_mu_SD,
            result.d_nu + result.d_mu_SD + result.d_mu_ID,
        ];
        result.d_for_D_states = [b'D', b'S', b'I'];
        result.d_for_I = [result.d_nu, result.d_nu + result.d_mu_SI];
        result.d_for_I_states = [b'I', b'S'];
        result.d_for_S = [
            result.d_eta,
            result.d_eta + result.d_mu_DS,
            result.d_eta + result.d_mu_DS + result.d_mu_IS,
        ];
        result.d_for_S_states = [b'S', b'D', b'I'];

        if let Some(alp_data) = alp_data_ {
            alp_data.d_memory_size_in_MB += (std::mem::size_of::<f64>()
                + std::mem::size_of::<q_elem>()) as f64
                * number_of_AA_ as f64
                / mb_bytes;
        }

        Ok(result)
    }

    pub fn lambda_equation(x_: f64, data: &data_for_lambda_equation) -> f64 {
        let mut res = 0.0;
        for i in 0..data.d_number_of_AA as usize {
            for j in 0..data.d_number_of_AA as usize {
                res += data.d_RR1[i] * data.d_RR2[j] * (x_ * data.d_smatr[i][j] as f64).exp();
            }
        }
        res - 1.0
    }
}

fn parse_next<T>(it: &mut std::str::SplitWhitespace<'_>, msg: &str) -> Result<T, Error>
where
    T: FromStr,
{
    it.next()
        .and_then(|s| s.parse::<T>().ok())
        .ok_or_else(|| Error::new(msg.to_string(), 3))
}

impl alp_data {
    #[allow(clippy::too_many_arguments)]
    pub fn new_from_files(
        rand_: i64,
        randout_: &str,
        open_: i64,
        open1_: i64,
        open2_: i64,
        epen_: i64,
        epen1_: i64,
        epen2_: i64,
        smatr_file_name_: &str,
        RR1_file_name_: &str,
        RR2_file_name_: &str,
        temperature_: f64,
        max_time_: f64,
        max_mem_: f64,
        eps_lambda_: f64,
        eps_K_: f64,
        insertions_after_deletions_: bool,
    ) -> Result<Self, Error> {
        let mut data = Self {
            d_rand_flag: true,
            d_rand_all: Some(Box::new(struct_for_randomization::default())),
            ..Self::default()
        };
        data.input_data_for_the_constructor(
            randout_,
            smatr_file_name_,
            RR1_file_name_,
            RR2_file_name_,
            rand_,
        )?;
        data.init_main_class_members(
            data.d_random_seed,
            randout_,
            open_,
            open1_,
            open2_,
            epen_,
            epen1_,
            epen2_,
            temperature_,
            max_time_,
            max_mem_,
            eps_lambda_,
            eps_K_,
            insertions_after_deletions_,
        )?;
        data.d_max_time_with_computation_parameters = -1.0;
        data.d_max_time_for_quick_tests = if max_time_ > 0.0 {
            0.25 * max_time_
        } else {
            1e99
        };
        Ok(data)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_from_parameters(
        rand_: i64,
        randomization_parameters_: Option<&struct_for_randomization>,
        open_: i64,
        open1_: i64,
        open2_: i64,
        epen_: i64,
        epen1_: i64,
        epen2_: i64,
        alphabetSize_: i64,
        substitutionScoreMatrix_: &[Vec<i64>],
        letterFreqs1_: &[f64],
        letterFreqs2_: &[f64],
        temperature_: f64,
        max_time_: f64,
        max_mem_: f64,
        eps_lambda_: f64,
        eps_K_: f64,
        insertions_after_deletions_: bool,
        max_time_for_quick_tests_: f64,
        max_time_with_computation_parameters_: f64,
    ) -> Result<Self, Error> {
        let mut data = Self {
            d_rand_flag: false,
            d_number_of_AA: alphabetSize_,
            d_smatr: substitutionScoreMatrix_.to_vec(),
            d_RR1: letterFreqs1_.to_vec(),
            d_RR2: letterFreqs2_.to_vec(),
            d_rand_all: Some(Box::new(struct_for_randomization::default())),
            ..Self::default()
        };
        if let Some(randomization_parameters) = randomization_parameters_ {
            data.d_rand_flag = true;
            data.d_rand_all = Some(Box::new(randomization_parameters.clone()));
        }
        data.init_main_class_members(
            rand_,
            "",
            open_,
            open1_,
            open2_,
            epen_,
            epen1_,
            epen2_,
            temperature_,
            max_time_,
            max_mem_,
            eps_lambda_,
            eps_K_,
            insertions_after_deletions_,
        )?;
        data.d_max_time_for_quick_tests = if max_time_for_quick_tests_ > 0.0 {
            max_time_for_quick_tests_
        } else if max_time_ > 0.0 {
            0.25 * max_time_
        } else {
            1e99
        };
        data.d_max_time_with_computation_parameters =
            if max_time_with_computation_parameters_ > 0.0 && !(max_time_ > 0.0) {
                max_time_with_computation_parameters_
            } else {
                1e99
            };
        let (rr1_sum, rr1_sum_elements) = Self::calculate_RR_sum(&mut data.d_RR1, alphabetSize_)?;
        data.d_RR1_sum = rr1_sum;
        data.d_RR1_sum_elements = rr1_sum_elements;
        let (rr2_sum, rr2_sum_elements) = Self::calculate_RR_sum(&mut data.d_RR2, alphabetSize_)?;
        data.d_RR2_sum = rr2_sum;
        data.d_RR2_sum_elements = rr2_sum_elements;
        Ok(data)
    }

    pub fn input_data_for_the_constructor(
        &mut self,
        randout_: &str,
        smatr_file_name_: &str,
        RR1_file_name_: &str,
        RR2_file_name_: &str,
        mut rand_: i64,
    ) -> Result<(), Error> {
        let (smatr, number_of_AA_smatr) = Self::read_smatr(smatr_file_name_)?;
        self.d_smatr = smatr;
        self.d_number_of_AA_smatr = number_of_AA_smatr;

        let (rr1, rr1_sum, rr1_sum_elements, number_of_AA_RR1) =
            Self::read_RR_with_sum(RR1_file_name_)?;
        self.d_RR1 = rr1;
        self.d_RR1_sum = rr1_sum;
        self.d_RR1_sum_elements = rr1_sum_elements;

        let (rr2, rr2_sum, rr2_sum_elements, number_of_AA_RR2) =
            Self::read_RR_with_sum(RR2_file_name_)?;
        self.d_RR2 = rr2;
        self.d_RR2_sum = rr2_sum;
        self.d_RR2_sum_elements = rr2_sum_elements;

        if number_of_AA_RR1 == number_of_AA_smatr {
            self.d_number_of_AA = number_of_AA_smatr;
        } else {
            return Err(Error::new(
                format!(
                    "Number of letters is different in the files {} and {}\n",
                    smatr_file_name_, RR1_file_name_
                ),
                3,
            ));
        }
        if number_of_AA_RR2 != number_of_AA_smatr {
            return Err(Error::new(
                format!(
                    "Number of letters is different in the files {} and {}\n",
                    smatr_file_name_, RR2_file_name_
                ),
                3,
            ));
        }

        if !randout_.is_empty() {
            match fs::read_to_string(randout_) {
                Ok(text) => {
                    let mut it = text.split_whitespace();
                    let mut rand_all = struct_for_randomization {
                        d_random_seed: parse_next(&mut it, "File is not correct\n")?,
                        ..Default::default()
                    };
                    if rand_all.d_random_seed < 0 {
                        return Err(Error::new(format!("File {} is not correct\n", randout_), 3));
                    }
                    rand_ = rand_all.d_random_seed;
                    let size: i64 = parse_next(&mut it, "File is not correct\n")?;
                    for _ in 0..size {
                        let tmp: i64 = parse_next(&mut it, "File is not correct\n")?;
                        if tmp < 0 {
                            return Err(Error::new(
                                format!("File {} is not correct\n", randout_),
                                3,
                            ));
                        }
                        rand_all
                            .d_first_stage_preliminary_realizations_numbers_ALP
                            .push(tmp);
                    }
                    let size: i64 = parse_next(&mut it, "File is not correct\n")?;
                    for _ in 0..size {
                        let tmp: i64 = parse_next(&mut it, "File is not correct\n")?;
                        if tmp < 0 {
                            return Err(Error::new(
                                format!("File {} is not correct\n", randout_),
                                3,
                            ));
                        }
                        rand_all.d_preliminary_realizations_numbers_ALP.push(tmp);
                    }
                    let size: i64 = parse_next(&mut it, "File is not correct\n")?;
                    for _ in 0..size {
                        let tmp: i64 = parse_next(&mut it, "File is not correct\n")?;
                        if tmp < 0 {
                            return Err(Error::new(
                                format!("File {} is not correct\n", randout_),
                                3,
                            ));
                        }
                        rand_all
                            .d_preliminary_realizations_numbers_killing
                            .push(tmp);
                    }
                    rand_all.d_total_realizations_number_with_ALP =
                        parse_next(&mut it, "File is not correct\n")?;
                    rand_all.d_total_realizations_number_with_killing =
                        parse_next(&mut it, "File is not correct\n")?;
                    if rand_all.d_total_realizations_number_with_ALP < 0
                        || rand_all.d_total_realizations_number_with_killing < 0
                    {
                        return Err(Error::new(format!("File {} is not correct\n", randout_), 3));
                    }
                    self.d_rand_flag = true;
                    self.d_rand_all = Some(Box::new(rand_all));
                }
                Err(_) => {
                    self.d_rand_flag = false;
                }
            }
        } else {
            self.d_rand_flag = false;
        }

        if !(self.d_number_of_AA > 0
            && !self.d_smatr.is_empty()
            && !self.d_RR1.is_empty()
            && !self.d_RR2.is_empty())
        {
            return Err(Error::new("Incorrect input parameters\n".to_string(), 1));
        }
        self.d_random_seed = rand_;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn init_main_class_members(
        &mut self,
        mut rand_: i64,
        randout_: &str,
        open_: i64,
        open1_: i64,
        open2_: i64,
        epen_: i64,
        epen1_: i64,
        epen2_: i64,
        temperature_: f64,
        max_time_: f64,
        max_mem_: f64,
        eps_lambda_: f64,
        eps_K_: f64,
        insertions_after_deletions_: bool,
    ) -> Result<(), Error> {
        self.d_randout = randout_.to_string();
        if !self.d_rand_flag && rand_ < 0 {
            rand_ = random_seed_from_time();
            self.d_rand_flag = false;
        }
        self.d_random_seed = rand_;
        self.d_number_of_AA_smatr = self.d_number_of_AA;
        self.d_sentinels_flag = false;
        self.d_memory_size_in_MB = 0.0;
        self.d_smatr_symmetric_flag = false;
        for t in 0..self.d_number_of_AA as usize {
            if self.d_RR1[t] != self.d_RR2[t] {
                self.d_smatr_symmetric_flag = false;
                break;
            }
        }
        self.d_insertions_after_deletions = insertions_after_deletions_;
        self.d_open = open_ + epen_;
        self.d_open1 = open1_ + epen1_;
        self.d_open2 = open2_ + epen2_;
        self.d_epen = epen_;
        self.d_epen1 = epen1_;
        self.d_epen2 = epen2_;
        self.d_max_time = max_time_;
        self.d_max_mem = max_mem_;
        self.d_eps_lambda = eps_lambda_;
        self.d_eps_K = eps_K_;
        self.d_minimum_realizations_number = 40;

        let is = importance_sampling::new(
            None,
            self.d_open,
            self.d_epen,
            temperature_,
            self.d_number_of_AA,
            &self.d_smatr,
            &self.d_RR1,
            &self.d_RR2,
        )?;
        self.d_is = Some(Box::new(is));

        self.d_r_i_dot = vec![0.0; self.d_number_of_AA as usize];
        self.d_r_dot_j = vec![0.0; self.d_number_of_AA as usize];
        let is_ref = self.d_is.as_ref().unwrap();
        for k in 0..self.d_number_of_AA as usize {
            self.d_r_i_dot[k] = 0.0;
            if self.d_RR1[k] != 0.0 {
                for i in 0..self.d_number_of_AA as usize {
                    if self.d_RR2[i] != 0.0 {
                        self.d_r_i_dot[k] += is_ref.d_exp_s[k][i] * self.d_RR2[i];
                    }
                }
            }
        }
        for k in 0..self.d_number_of_AA as usize {
            self.d_r_dot_j[k] = 0.0;
            if self.d_RR2[k] != 0.0 {
                for i in 0..self.d_number_of_AA as usize {
                    if self.d_RR1[i] != 0.0 {
                        self.d_r_dot_j[k] += is_ref.d_exp_s[i][k] * self.d_RR1[i];
                    }
                }
            }
        }

        let tmp_size1 = i64::MAX as f64;
        let tmp_size = Tmin(
            tmp_size1,
            (mb_bytes * self.d_max_mem / self.d_minimum_realizations_number as f64)
                / ((std::mem::size_of::<f64>() * 12) as f64
                    + (std::mem::size_of::<i64>() * 17) as f64),
        );
        self.d_dim1_tmp = tmp_size as i64;
        self.d_dim2_tmp = tmp_size as i64;
        Ok(())
    }

    pub fn release_memory(&mut self) {
        self.d_RR1.clear();
        self.d_RR1_sum.clear();
        self.d_RR1_sum_elements.clear();
        self.d_RR2.clear();
        self.d_RR2_sum.clear();
        self.d_RR2_sum_elements.clear();
        self.d_smatr.clear();
        self.d_is = None;
        self.d_r_i_dot.clear();
        self.d_r_dot_j.clear();
        self.d_rand_all = None;
        self.d_memory_size_in_MB = 0.0;
    }

    pub fn check_out_file(out_file_name_: &str) -> Result<Option<bool>, Error> {
        let text = match fs::read_to_string(out_file_name_) {
            Ok(text) => text,
            Err(_) => return Ok(None),
        };
        let first_line = text.lines().next().unwrap_or("");
        if !first_line.contains("number of realizations with killing") {
            return Err(Error::new(
                format!(
                    "The output file {} exists and does not have the correct format;\nplease delete the file and rerun the program\n",
                    out_file_name_
                ),
                3,
            ));
        }
        Ok(Some(first_line.contains("0.5*")))
    }

    pub fn read_smatr(smatr_file_name_: &str) -> Result<(Vec<Vec<i64>>, i64), Error> {
        let text = fs::read_to_string(smatr_file_name_).map_err(|_| {
            Error::new(
                format!("Error - file {} is not found\n", smatr_file_name_),
                3,
            )
        })?;
        let mut it = text.split_whitespace();
        let number_of_AA_smatr_: i64 = parse_next(
            &mut it,
            "Error - number of letters in the scoring matrix file must be greater than 0\n",
        )?;
        if number_of_AA_smatr_ <= 0 {
            return Err(Error::new(
                "Error - number of letters in the scoring matrix file must be greater than 0\n"
                    .to_string(),
                3,
            ));
        }
        let n = number_of_AA_smatr_ as usize;
        let mut smatr_ = vec![vec![0i64; n]; n];
        for row in smatr_.iter_mut().take(n) {
            for cell in row.iter_mut().take(n) {
                *cell = parse_next(&mut it, "Error - scoring matrix file is not correct\n")?;
            }
        }
        Ok((smatr_, number_of_AA_smatr_))
    }

    pub fn read_RR_with_sum(
        RR_file_name_: &str,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<i64>, i64), Error> {
        let (mut RR_, number_of_AA_RR_) = Self::read_RR(RR_file_name_)?;
        let (RR_sum_, RR_sum_elements_) = Self::calculate_RR_sum(&mut RR_, number_of_AA_RR_)?;
        Ok((RR_, RR_sum_, RR_sum_elements_, number_of_AA_RR_))
    }

    pub fn read_RR(RR_file_name_: &str) -> Result<(Vec<f64>, i64), Error> {
        let text = fs::read_to_string(RR_file_name_)
            .map_err(|_| Error::new(format!("Error - file {} is not found\n", RR_file_name_), 3))?;
        let mut it = text.split_whitespace();
        let number_of_AA_RR_: i64 = parse_next(
            &mut it,
            "Error - number of letters in the probabilities file must be greater than 0\n",
        )?;
        if number_of_AA_RR_ <= 0 {
            return Err(Error::new(
                "Error - number of letters in the probabilities file must be greater than 0\n"
                    .to_string(),
                3,
            ));
        }
        let mut RR_ = vec![0.0; number_of_AA_RR_ as usize];
        let mut sum_tmp = 0.0;
        for item in RR_.iter_mut().take(number_of_AA_RR_ as usize) {
            *item = parse_next(&mut it, "Error - probabilities file is not correct\n")?;
            if *item < 0.0 {
                return Err(Error::new(
                    format!(
                        "Error - the frequencies defined in the file {} must be non-negative\n",
                        RR_file_name_
                    ),
                    3,
                ));
            }
            sum_tmp += *item;
        }
        Self::check_RR_sum(sum_tmp, number_of_AA_RR_, RR_file_name_)?;
        Ok((RR_, number_of_AA_RR_))
    }

    pub fn calculate_RR_sum(
        RR_: &mut [f64],
        number_of_AA_RR_: i64,
    ) -> Result<(Vec<f64>, Vec<i64>), Error> {
        if number_of_AA_RR_ <= 0 {
            return Err(Error::new(
                "Error - number of letters in the probabilities file must be greater than 0\n"
                    .to_string(),
                3,
            ));
        }
        let n = number_of_AA_RR_ as usize;
        let mut RR_sum_ = vec![0.0; n];
        let mut RR_sum_elements_ = vec![0; n];
        for i in 0..n {
            if RR_[i] < 0.0 {
                return Err(Error::new(
                    "Error - the frequencies must be non-negative\n".to_string(),
                    3,
                ));
            }
            RR_sum_[i] = if i != 0 {
                RR_sum_[i - 1] + RR_[i]
            } else {
                RR_[i]
            };
            RR_sum_elements_[i] = i as i64;
        }
        let sum_tmp = RR_sum_[n - 1];
        Self::check_RR_sum(sum_tmp, number_of_AA_RR_, "")?;
        if sum_tmp > 0.0 {
            for i in 0..n {
                RR_[i] /= sum_tmp;
                RR_sum_[i] /= sum_tmp;
            }
        }
        Ok((RR_sum_, RR_sum_elements_))
    }

    pub fn check_RR_sum(
        sum_tmp_: f64,
        number_of_AA_RR_: i64,
        RR_file_name_: &str,
    ) -> Result<(), Error> {
        if number_of_AA_RR_ <= 0 {
            return Err(Error::new(
                "Error - number of letters in the probabilities file must be greater than 0\n"
                    .to_string(),
                3,
            ));
        }
        let diff_tmp = (sum_tmp_ - 1.0).abs();
        if diff_tmp > 0.0 {
            let lg_diff = -((diff_tmp.ln() - (number_of_AA_RR_ as f64).ln()) / 10.0f64.ln());
            let lg_eps = -f64::EPSILON.ln() / 10.0f64.ln() - 1.0;
            if lg_diff < lg_eps && sum_tmp_ <= 0.0 {
                if RR_file_name_.is_empty() {
                    return Err(Error::new(
                        "Error: the sum of the probabilities is non-positive\n".to_string(),
                        3,
                    ));
                }
                return Err(Error::new(
                    format!(
                        "Error: the sum of the probabilities from the file {} is non-positive\n",
                        RR_file_name_
                    ),
                    3,
                ));
            }
        }
        Ok(())
    }

    pub fn long_to_string(number_: i64) -> String {
        let mut res_ = String::new();
        let tmp_string = if number_ < 0 { "-" } else { "" };
        let mut number_ = number_.abs();
        loop {
            let reminder = number_ % 10;
            number_ = (number_ - reminder) / 10;
            res_.insert(0, Self::digit_to_string(reminder));
            if number_ == 0 {
                break;
            }
        }
        format!("{}{}", tmp_string, res_)
    }

    pub fn digit_to_string(digit_: i64) -> char {
        match digit_ {
            0 => '0',
            1 => '1',
            2 => '2',
            3 => '3',
            4 => '4',
            5 => '5',
            6 => '6',
            7 => '7',
            8 => '8',
            9 => '9',
            _ => '?',
        }
    }

    pub fn the_value_is_double(str_: &str, val_: &mut f64) -> bool {
        if str_.is_empty() {
            return false;
        }
        match str_.parse::<f64>() {
            Ok(v) => {
                *val_ = v;
                true
            }
            Err(_) => false,
        }
    }

    pub fn the_value_is_long(str_: &str, val_: &mut i64) -> bool {
        if str_.is_empty() {
            return false;
        }
        let bytes = str_.as_bytes();
        if !(bytes[0] == b'+' || bytes[0] == b'-' || bytes[0].is_ascii_digit()) {
            return false;
        }
        let start_digit = if bytes[0] == b'+' || bytes[0] == b'-' {
            1
        } else {
            0
        };
        if str_.len() - start_digit == 0 {
            return false;
        }
        if !str_[start_digit..].bytes().all(|b| b.is_ascii_digit()) {
            return false;
        }
        let mut digits = str_[start_digit..].trim_start_matches('0').to_string();
        if digits.is_empty() {
            digits.push('0');
        }
        if digits.len() > 10 || digits.is_empty() {
            return false;
        }
        if digits.len() == 10 {
            let first = digits.as_bytes()[0];
            if !(first == b'1' || first == b'2') {
                return false;
            }
            if first == b'2' {
                let val2 = digits[1..].parse::<i64>().unwrap_or(i64::MAX);
                let positive = start_digit == 0 || bytes[0] != b'-';
                if positive {
                    if val2 > 147_483_647 {
                        return false;
                    }
                } else if val2 > 147_483_648 {
                    return false;
                }
            }
        }
        match str_.parse::<i64>() {
            Ok(v) => {
                *val_ = v;
                true
            }
            Err(_) => false,
        }
    }

    pub fn get_memory_for_matrix<T: Default + Clone>(dim_: i64) -> Vec<Vec<T>> {
        vec![vec![T::default(); dim_ as usize]; dim_ as usize]
    }

    pub fn delete_memory_for_matrix<T>(matr_: &mut Vec<Vec<T>>) {
        matr_.clear();
    }

    pub fn random_long(value_: f64, dim_: i64) -> Result<i64, Error> {
        if !(0.0..=1.0).contains(&value_) || dim_ <= 0 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        if dim_ == 1 {
            return Ok(0);
        }
        let tmp = (value_ * dim_ as f64).floor() as i64;
        Ok(Tmin(tmp, dim_ - 1))
    }

    pub fn ran2(&self) -> f64 {
        crate::stats::pvalues::pvalues::ran3()
    }

    pub fn random_long_element<T: Copy>(
        value_: f64,
        dim_: i64,
        sum_distr_: &[f64],
        elements_: &[T],
    ) -> Result<T, Error> {
        if value_ < 0.0 || value_ > 1.0 {
            return Err(Error::new(
                "Unexpected error in alp_data::random_long\n".to_string(),
                4,
            ));
        }
        let mut v1 = 0;
        let mut v2 = dim_;
        while v2 - v1 > 1 {
            let v3 = round(&(f64::from((v2 + v1) as i32) / 2.0)) as i64;
            if sum_distr_[(v3 - 1) as usize] == value_ {
                v2 = v3;
                break;
            }
            if sum_distr_[(v3 - 1) as usize] > value_ {
                v2 = v3;
            } else {
                v1 = v3;
            }
        }

        let v2_1 = v2 - 1;
        let mut v2_minus = -1;
        for j in (1..=v2_1).rev() {
            if sum_distr_[j as usize] != sum_distr_[(j - 1) as usize] {
                v2_minus = j;
                break;
            }
        }
        if v2_minus < 0 && sum_distr_[0] > 0.0 {
            v2_minus = 0;
        }
        if v2_minus >= 0 {
            return Ok(elements_[v2_minus as usize]);
        }

        let mut v2_plus = -1;
        for j in v2..dim_ {
            if sum_distr_[j as usize] != sum_distr_[(j - 1) as usize] {
                v2_plus = j;
                break;
            }
        }
        if v2_minus < 0 {
            Err(Error::new(
                "Unexpected error in alp_data::random_long\n".to_string(),
                1,
            ))
        } else {
            Ok(elements_[v2_plus as usize])
        }
    }

    pub fn error_of_the_sum(v1_error_: f64, v2_error_: f64) -> f64 {
        if v1_error_ >= 1e100 || v2_error_ >= 1e100 {
            return 1e100;
        }
        (v1_error_ * v1_error_ + v2_error_ * v2_error_).sqrt()
    }

    pub fn error_of_the_product(v1_: f64, v1_error_: f64, v2_: f64, v2_error_: f64) -> f64 {
        if v1_error_ >= 1e100 || v2_error_ >= 1e100 {
            return 1e100;
        }

        let a1 = (v1_ + v1_error_) * (v2_ + v2_error_);
        let a2 = (v1_ - v1_error_) * (v2_ + v2_error_);
        let a3 = (v1_ + v1_error_) * (v2_ - v2_error_);
        let a4 = (v1_ - v1_error_) * (v2_ - v2_error_);
        let a = v1_ * v2_;

        Tmax4(
            (a1 - a).abs(),
            (a2 - a).abs(),
            (a3 - a).abs(),
            (a4 - a).abs(),
        )
    }

    pub fn error_of_the_sqrt(v1_: f64, v1_error_: f64) -> f64 {
        if v1_error_ >= 1e100 || v1_ < 0.0 {
            return 1e100;
        }

        let s = v1_.sqrt();
        let s1 = Tmax(0.0, v1_ - v1_error_).sqrt();
        let s2 = Tmax(0.0, v1_ + v1_error_).sqrt();
        Tmax((s - s1).abs(), (s - s2).abs())
    }

    pub fn error_of_the_ratio(v1_: f64, v1_error_: f64, v2_: f64, v2_error_: f64) -> f64 {
        if v1_error_ >= 1e100 || v2_error_ >= 1e100 {
            return 1e100;
        }
        if v2_ == 0.0 {
            return 1e100;
        }
        if v1_ == 0.0 && v1_error_ == 0.0 {
            return 0.0;
        }

        let a = v1_ / v2_;
        if (v2_ + v2_error_) * v2_ <= 0.0 {
            let a3 = (v1_ + v1_error_) / (v2_ - v2_error_);
            let a4 = (v1_ - v1_error_) / (v2_ - v2_error_);
            return Tmax((a - a3).abs(), (a - a4).abs());
        }
        if (v2_ - v2_error_) * v2_ <= 0.0 {
            let a1 = (v1_ + v1_error_) / (v2_ + v2_error_);
            let a2 = (v1_ - v1_error_) / (v2_ + v2_error_);
            return Tmax((a - a1).abs(), (a - a2).abs());
        }

        let a1 = (v1_ + v1_error_) / (v2_ + v2_error_);
        let a2 = (v1_ - v1_error_) / (v2_ + v2_error_);
        let a3 = (v1_ + v1_error_) / (v2_ - v2_error_);
        let a4 = (v1_ - v1_error_) / (v2_ - v2_error_);
        Tmax4(
            (a - a1).abs(),
            (a - a2).abs(),
            (a - a3).abs(),
            (a - a4).abs(),
        )
    }

    pub fn error_of_the_lg(v1_: f64, v1_error_: f64) -> f64 {
        if v1_error_ >= 1e100 || v1_ <= 0.0 {
            return 1e100;
        }
        Tmin(
            (v1_.ln() / 10.0f64.ln()).abs(),
            v1_error_ / v1_ / 10.0f64.ln(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_helpers_match_cpp_formulas() {
        assert_eq!(alp_data::error_of_the_sum(1e100, 1.0), 1e100);
        assert!((alp_data::error_of_the_sum(3.0, 4.0) - 5.0).abs() < 1e-15);
        assert!((alp_data::error_of_the_product(2.0, 0.1, 3.0, 0.2) - 0.72).abs() < 1e-15);
        assert!((alp_data::error_of_the_sqrt(4.0, 0.75) - (2.0 - 3.25f64.sqrt())).abs() < 1e-15);
        assert_eq!(alp_data::error_of_the_ratio(0.0, 0.0, 2.0, 0.5), 0.0);
        assert!((alp_data::error_of_the_lg(100.0, 1.0) - (0.01 / 10.0f64.ln())).abs() < 1e-15);
    }

    #[test]
    fn test_array_positive_growth_and_updates() {
        let mut a = array_positive::<i64>::new(None).unwrap();
        a.set_elem(12, 5);
        assert_eq!(a.d_dim, 19);
        assert_eq!(a.d_elem[12], 5);
        a.increase_elem_by_1(12);
        a.increase_elem_by_x(21, 3);
        assert_eq!(a.d_dim, 29);
        assert_eq!(a.d_elem[12], 6);
        assert_eq!(a.d_elem[21], 3);
    }

    #[test]
    fn test_array_growth_left_right_and_copy() {
        let mut a = array::<i64>::new(None);
        a.set_elem(5, 10);
        a.set_elem(-3, 7);
        a.increase_elem_by_1(-3);
        assert!(a.d_ind0 <= -3);
        assert!(a.d_dim_plus_d_ind0 >= 5);
        assert_eq!(a.d_elem[(5 - a.d_ind0) as usize], 10);
        assert_eq!(a.d_elem[(-3 - a.d_ind0) as usize], 8);

        let mut b = array::<i64>::new(None);
        b.set_elems(&a);
        assert_eq!(b.d_ind0, a.d_ind0);
        assert_eq!(b.d_dim_plus_d_ind0, a.d_dim_plus_d_ind0);
        assert_eq!(b.d_elem[(5 - b.d_ind0) as usize], 10);
        assert_eq!(b.d_elem[(-3 - b.d_ind0) as usize], 8);
    }

    #[test]
    fn test_read_matrix_rr_and_random_helpers() {
        let dir = std::env::temp_dir();
        let smatr = dir.join(format!("diamond_rs_smatr_{}.txt", std::process::id()));
        let rr = dir.join(format!("diamond_rs_rr_{}.txt", std::process::id()));
        std::fs::write(&smatr, "2\n1 -1\n-2 3\n").unwrap();
        std::fs::write(&rr, "3\n2\n3\n5\n").unwrap();

        let (matrix, n) = alp_data::read_smatr(smatr.to_str().unwrap()).unwrap();
        assert_eq!(n, 2);
        assert_eq!(matrix, vec![vec![1, -1], vec![-2, 3]]);

        let (mut probs, n_rr) = alp_data::read_RR(rr.to_str().unwrap()).unwrap();
        assert_eq!(n_rr, 3);
        assert_eq!(probs, vec![2.0, 3.0, 5.0]);
        let (sum, elems) = alp_data::calculate_RR_sum(&mut probs, n_rr).unwrap();
        assert_eq!(probs, vec![0.2, 0.3, 0.5]);
        assert_eq!(sum, vec![0.2, 0.5, 1.0]);
        assert_eq!(elems, vec![0, 1, 2]);

        assert_eq!(alp_data::random_long(0.0, 4).unwrap(), 0);
        assert_eq!(alp_data::random_long(1.0, 4).unwrap(), 3);
        let picked =
            alp_data::random_long_element(0.6, 3, &[0.2, 0.5, 1.0], &[10, 20, 30]).unwrap();
        assert_eq!(picked, 30);

        let _ = std::fs::remove_file(smatr);
        let _ = std::fs::remove_file(rr);
    }

    #[test]
    fn test_conversion_and_lambda_equation_helpers() {
        assert_eq!(alp_data::long_to_string(-12345), "-12345");
        assert_eq!(alp_data::long_to_string(0), "0");
        assert_eq!(alp_data::digit_to_string(11), '?');
        let mut d = 0.0;
        assert!(alp_data::the_value_is_double("1.25", &mut d));
        assert_eq!(d, 1.25);
        assert!(!alp_data::the_value_is_double("1.2x", &mut d));
        let mut l = 0;
        assert!(alp_data::the_value_is_long("-2147483648", &mut l));
        assert_eq!(l, -2147483648);
        assert!(!alp_data::the_value_is_long("2147483648", &mut l));

        let data = data_for_lambda_equation {
            d_number_of_AA: 2,
            d_smatr: vec![vec![-1, 2], vec![3, -2]],
            d_RR1: vec![0.25, 0.75],
            d_RR2: vec![0.5, 0.5],
        };
        let v = importance_sampling::lambda_equation(0.0, &data);
        assert!(v.abs() < 1e-15);
    }

    #[test]
    fn test_importance_sampling_constructor_for_small_scheme() {
        let smatr = vec![vec![-2, 1], vec![1, -2]];
        let rr = vec![0.5, 0.5];
        let is = importance_sampling::new(None, 12, 1, 1.0, 2, &smatr, &rr, &rr).unwrap();

        let golden_ratio = (1.0 + 5.0f64.sqrt()) / 2.0;
        assert!((is.d_ungap_lambda - golden_ratio.ln()).abs() < 1e-5);
        assert_eq!(is.d_elements.len(), 4);
        assert_eq!(is.d_elements[0], q_elem { d_a: 0, d_b: 0 });
        assert_eq!(is.d_elements[3], q_elem { d_a: 1, d_b: 1 });
        assert!((is.d_elements_values[3] - 1.0).abs() < 1e-12);
        assert_eq!(is.d_for_D_states, [b'D', b'S', b'I']);
        assert_eq!(is.d_for_I_states, [b'I', b'S']);
        assert_eq!(is.d_for_S_states, [b'S', b'D', b'I']);
        assert!(is.d_for_D[0] <= is.d_for_D[1] && is.d_for_D[1] <= is.d_for_D[2]);
        assert!(is.d_for_S[0] <= is.d_for_S[1] && is.d_for_S[1] <= is.d_for_S[2]);
    }

    #[test]
    fn test_alp_data_from_parameters_initializes_runtime_state() {
        let smatr = vec![vec![-2, 1], vec![1, -2]];
        let rr = vec![0.5, 0.5];
        let data = alp_data::new_from_parameters(
            123, None, 11, 13, 17, 1, 2, 3, 2, &smatr, &rr, &rr, 1.0, 10.0, 128.0, 0.05, 0.05,
            true, -1.0, -1.0,
        )
        .unwrap();

        assert_eq!(data.d_random_seed, 123);
        assert_eq!(data.d_open, 12);
        assert_eq!(data.d_open1, 15);
        assert_eq!(data.d_open2, 20);
        assert_eq!(data.d_epen, 1);
        assert_eq!(data.d_epen1, 2);
        assert_eq!(data.d_epen2, 3);
        assert!(data.d_insertions_after_deletions);
        assert_eq!(data.d_minimum_realizations_number, 40);
        assert_eq!(data.d_RR1_sum, vec![0.5, 1.0]);
        assert_eq!(data.d_RR2_sum_elements, vec![0, 1]);
        assert!(data.d_is.is_some());
        assert_eq!(data.d_r_i_dot.len(), 2);
        assert_eq!(data.d_r_dot_j.len(), 2);
        assert_eq!(data.d_max_time_for_quick_tests, 2.5);
        assert_eq!(data.d_max_time_with_computation_parameters, 1e99);
    }

    #[test]
    fn test_alp_data_from_files_and_check_out_file() {
        let dir = std::env::temp_dir();
        let tag = format!("diamond_rs_alpdata_{}", std::process::id());
        let smatr = dir.join(format!("{}_smatr.txt", tag));
        let rr1 = dir.join(format!("{}_rr1.txt", tag));
        let rr2 = dir.join(format!("{}_rr2.txt", tag));
        let rand = dir.join(format!("{}_rand.txt", tag));
        let out = dir.join(format!("{}_out.txt", tag));
        std::fs::write(&smatr, "2\n-2 1\n1 -2\n").unwrap();
        std::fs::write(&rr1, "2\n0.5\n0.5\n").unwrap();
        std::fs::write(&rr2, "2\n0.5\n0.5\n").unwrap();
        std::fs::write(&rand, "77\n2\n3\n5\n1\n7\n1\n9\n11\n13\n").unwrap();
        std::fs::write(&out, "number of realizations with killing 0.5* symmetric\n").unwrap();

        let data = alp_data::new_from_files(
            -1,
            rand.to_str().unwrap(),
            11,
            11,
            11,
            1,
            1,
            1,
            smatr.to_str().unwrap(),
            rr1.to_str().unwrap(),
            rr2.to_str().unwrap(),
            1.0,
            -1.0,
            128.0,
            0.05,
            0.05,
            false,
        )
        .unwrap();

        assert_eq!(data.d_random_seed, 77);
        assert!(data.d_rand_flag);
        let rand_all = data.d_rand_all.as_ref().unwrap();
        assert_eq!(
            rand_all.d_first_stage_preliminary_realizations_numbers_ALP,
            vec![3, 5]
        );
        assert_eq!(rand_all.d_preliminary_realizations_numbers_ALP, vec![7]);
        assert_eq!(rand_all.d_preliminary_realizations_numbers_killing, vec![9]);
        assert_eq!(rand_all.d_total_realizations_number_with_ALP, 11);
        assert_eq!(rand_all.d_total_realizations_number_with_killing, 13);
        assert_eq!(
            alp_data::check_out_file(out.to_str().unwrap()).unwrap(),
            Some(true)
        );

        let _ = std::fs::remove_file(smatr);
        let _ = std::fs::remove_file(rr1);
        let _ = std::fs::remove_file(rr2);
        let _ = std::fs::remove_file(rand);
        let _ = std::fs::remove_file(out);
    }
}
