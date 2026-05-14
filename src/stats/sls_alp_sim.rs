#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use std::fs;

use crate::stats::sls_alp::alp;
use crate::stats::sls_alp_data::{alp_data, array_positive, mb_bytes};
use crate::stats::sls_alp_regression::alp_reg;
use crate::stats::sls_basic::{get_current_time, round, Error, Tmax, Tmin};

pub static calculate_C_S_constant_flag: bool = true;
pub const quick_tests_trials_number: i64 = 100;

#[derive(Debug, Clone, PartialEq)]
pub struct struct_for_lambda_calculation {
    pub d_alp_distr: Vec<array_positive<f64>>,
    pub d_alp_distr_errors: Vec<array_positive<f64>>,
    pub d_nalp: i64,
    pub d_f_error: f64,
    pub d_last_sum: f64,
    pub d_last_sum_error: f64,
    pub d_calculate_alp_number: bool,
    pub d_alp_number: i64,
}

impl Default for struct_for_lambda_calculation {
    fn default() -> Self {
        Self {
            d_alp_distr: Vec::new(),
            d_alp_distr_errors: Vec::new(),
            d_nalp: 0,
            d_f_error: 0.0,
            d_last_sum: 0.0,
            d_last_sum_error: 0.0,
            d_calculate_alp_number: false,
            d_alp_number: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct alp_sim {
    pub d_alp_data: alp_data,
    pub d_alp_obj: Vec<Option<Box<alp>>>,
    pub d_n_alp_obj: i64,
    pub d_lambda_tmp: array_positive<f64>,
    pub d_lambda_tmp_errors: array_positive<f64>,
    pub d_C_tmp: array_positive<f64>,
    pub d_C_tmp_errors: array_positive<f64>,
    pub d_mult_number: i64,
    pub m_Lambda: f64,
    pub m_LambdaError: f64,
    pub m_K: f64,
    pub m_KError: f64,
    pub m_C: f64,
    pub m_CError: f64,
    pub m_Sigma: f64,
    pub m_SigmaError: f64,
    pub m_GaplessAlpha: f64,
    pub m_GaplessAlphaError: f64,
    pub m_GaplessA: f64,
    pub m_GaplessAError: f64,
    pub m_AlphaI: f64,
    pub m_AlphaIError: f64,
    pub m_AlphaJ: f64,
    pub m_AlphaJError: f64,
    pub m_AI: f64,
    pub m_AIError: f64,
    pub m_AJ: f64,
    pub m_AJError: f64,
    pub m_CalcTime: f64,
    pub m_G: i64,
    pub m_G1: i64,
    pub m_G2: i64,
    pub m_LambdaSbs: Vec<f64>,
    pub m_KSbs: Vec<f64>,
    pub m_CSbs: Vec<f64>,
    pub m_SigmaSbs: Vec<f64>,
    pub m_AlphaISbs: Vec<f64>,
    pub m_AlphaJSbs: Vec<f64>,
    pub m_AISbs: Vec<f64>,
    pub m_AJSbs: Vec<f64>,
}

impl alp_sim {
    pub fn new(alp_data_: alp_data) -> Result<Self, Error> {
        let mut result = Self {
            d_alp_data: alp_data_,
            d_alp_obj: Vec::new(),
            d_n_alp_obj: 0,
            d_lambda_tmp: array_positive::new(None)?,
            d_lambda_tmp_errors: array_positive::new(None)?,
            d_C_tmp: array_positive::new(None)?,
            d_C_tmp_errors: array_positive::new(None)?,
            d_mult_number: 0,
            m_Lambda: 0.0,
            m_LambdaError: 0.0,
            m_K: 0.0,
            m_KError: 0.0,
            m_C: 0.0,
            m_CError: 0.0,
            m_Sigma: 0.0,
            m_SigmaError: 0.0,
            m_GaplessAlpha: 0.0,
            m_GaplessAlphaError: 0.0,
            m_GaplessA: 0.0,
            m_GaplessAError: 0.0,
            m_AlphaI: 0.0,
            m_AlphaIError: 0.0,
            m_AlphaJ: 0.0,
            m_AlphaJError: 0.0,
            m_AI: 0.0,
            m_AIError: 0.0,
            m_AJ: 0.0,
            m_AJError: 0.0,
            m_CalcTime: 0.0,
            m_G: 0,
            m_G1: 0,
            m_G2: 0,
            m_LambdaSbs: Vec::new(),
            m_KSbs: Vec::new(),
            m_CSbs: Vec::new(),
            m_SigmaSbs: Vec::new(),
            m_AlphaISbs: Vec::new(),
            m_AlphaJSbs: Vec::new(),
            m_AISbs: Vec::new(),
            m_AJSbs: Vec::new(),
        };
        result.d_alp_data.d_memory_size_in_MB += (std::mem::size_of::<array_positive<f64>>() * 4
            + std::mem::size_of::<Vec<Option<Box<alp>>>>())
            as f64
            / mb_bytes;
        Ok(result)
    }

    pub fn new_with_simulation(alp_data_: alp_data) -> Result<Self, Error> {
        let mut result = Self::new(alp_data_)?;

        let memory_before1 = result.d_alp_data.d_memory_size_in_MB;
        let mut time_before1 = 0.0;
        get_current_time(&mut time_before1);
        result.d_alp_data.d_time_before1 = time_before1;

        if !result.d_alp_data.d_rand_flag {
            let rand_all = result
                .d_alp_data
                .d_rand_all
                .as_mut()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            rand_all.d_random_seed = result.d_alp_data.d_random_seed;
        }

        let mut M_min = 0;
        let mut nalp = 0;
        let mut time_after_tmp = 0.0;

        result.d_lambda_tmp = array_positive::new(None)?;
        result.d_lambda_tmp_errors = array_positive::new(None)?;
        result.d_C_tmp = array_positive::new(None)?;
        result.d_C_tmp_errors = array_positive::new(None)?;

        result.quick_test(
            quick_tests_trials_number,
            result.d_alp_data.d_max_time_for_quick_tests,
        )?;

        let maximum_number_of_realizations_for_preliminary_simulation = 1000;
        let mut loop_break_flag;
        let mut rand_i = -1;
        let mut lambda_accuracy_flag = true;
        let mut sim_number = 1;

        loop {
            let mut nalp_lambda = 0;
            let C_calculation = false;
            let number_tmp;
            let mut check_time_flag = true;

            if result.d_alp_data.d_rand_flag {
                check_time_flag = false;
                rand_i += 1;
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_ref()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                number_tmp =
                    rand_all.d_first_stage_preliminary_realizations_numbers_ALP[rand_i as usize];
            } else {
                number_tmp = Tmin(
                    maximum_number_of_realizations_for_preliminary_simulation - 1,
                    result.d_n_alp_obj
                        + sim_number * result.d_alp_data.d_minimum_realizations_number
                        - 1,
                );
            }

            result.get_minimal_simulation(
                0,
                number_tmp,
                &mut M_min,
                &mut nalp,
                &mut nalp_lambda,
                C_calculation,
                check_time_flag,
            )?;

            if !result.d_alp_data.d_rand_flag {
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_mut()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                rand_all
                    .d_first_stage_preliminary_realizations_numbers_ALP
                    .push(number_tmp);
            }

            sim_number *= 2;

            if result.d_lambda_tmp.d_elem[nalp as usize] >= 0.0
                && result.d_lambda_tmp_errors.d_elem[nalp as usize]
                    / result.d_lambda_tmp.d_elem[nalp as usize]
                    < result.d_alp_data.d_eps_lambda
            {
                lambda_accuracy_flag = false;
            }

            get_current_time(&mut time_after_tmp);

            if number_tmp >= maximum_number_of_realizations_for_preliminary_simulation - 1
                && !result.d_alp_data.d_rand_flag
            {
                break;
            }

            if result.d_alp_data.d_rand_flag {
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_ref()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                loop_break_flag = rand_i
                    >= rand_all
                        .d_first_stage_preliminary_realizations_numbers_ALP
                        .len() as i64
                        - 1;
            } else {
                loop_break_flag = !(maximum_number_of_realizations_for_preliminary_simulation
                    > result.d_n_alp_obj - 1
                    && lambda_accuracy_flag
                    && ((time_after_tmp - time_before1) <= 0.0
                        || ((time_after_tmp - time_before1)
                            < 0.01 * result.d_alp_data.d_max_time
                            && result.d_alp_data.d_memory_size_in_MB
                                < result.d_alp_data.d_max_mem)));
            }

            if loop_break_flag {
                break;
            }
        }

        let memory_after1 = result.d_alp_data.d_memory_size_in_MB;
        let mut time_after1 = 0.0;
        get_current_time(&mut time_after1);

        if memory_after1 <= memory_before1 {
            result.d_alp_obj.clear();
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        let limit_by_memory = Tmin(
            i64::MAX - 1,
            round(
                &(result.d_alp_data.d_max_mem / 2.0 / (memory_after1 - memory_before1)
                    * result.d_n_alp_obj as f64),
            ) as i64,
        );

        let limit_by_time = if time_after1 <= time_before1 {
            i64::MAX
        } else {
            Tmin(
                i64::MAX - 1,
                round(
                    &(result.d_alp_data.d_max_time / 3.0 / 4.0 / (time_after1 - time_before1)
                        * result.d_n_alp_obj as f64),
                ) as i64,
            )
        };

        let mut realizations_number2 = Tmin(
            round(&(limit_by_time as f64)) as i64 - 1,
            limit_by_memory - 1,
        );
        realizations_number2 = Tmin(
            maximum_number_of_realizations_for_preliminary_simulation - 1,
            realizations_number2,
        );
        realizations_number2 = Tmax(result.d_n_alp_obj - 1, realizations_number2);

        result.d_lambda_tmp = array_positive::new(None)?;
        result.d_lambda_tmp_errors = array_positive::new(None)?;
        result.d_C_tmp = array_positive::new(None)?;
        result.d_C_tmp_errors = array_positive::new(None)?;

        let mut nalp_lambda = 0;
        let C_calculation = false;
        let mut lambda: f64;
        let eps_K = result.d_alp_data.d_eps_K;
        let mut K_C = 0.0;
        let mut K_C_error = 0.0;
        let mut level = 0;
        let mut diff_opt = 0;

        rand_i = -1;
        let mut time_before_ALP = 0.0;
        let mut time_during_ALP = 0.0;
        let mut number_of_realizations_with_ALP = Tmin(
            realizations_number2,
            result.d_n_alp_obj - 1 + result.d_alp_data.d_minimum_realizations_number,
        );
        let mut number_of_realizations_with_ALP_pred;
        get_current_time(&mut time_before_ALP);

        loop {
            let mut check_time_flag = true;
            if result.d_alp_data.d_rand_flag {
                check_time_flag = false;
                rand_i += 1;
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_ref()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                number_of_realizations_with_ALP =
                    rand_all.d_preliminary_realizations_numbers_ALP[rand_i as usize];
            }

            result.get_minimal_simulation(
                0,
                number_of_realizations_with_ALP,
                &mut M_min,
                &mut nalp,
                &mut nalp_lambda,
                C_calculation,
                check_time_flag,
            )?;

            if !result.d_alp_data.d_rand_flag {
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_mut()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                rand_all
                    .d_preliminary_realizations_numbers_ALP
                    .push(number_of_realizations_with_ALP);
            }

            lambda = result.d_lambda_tmp.d_elem[nalp as usize];
            let mut tmp_lambda = 2.0;
            if result.d_lambda_tmp.d_elem[nalp as usize] > 0.0 {
                tmp_lambda = (result.d_lambda_tmp_errors.d_elem[nalp as usize]
                    / result.d_lambda_tmp.d_elem[nalp as usize])
                    / result.d_alp_data.d_eps_lambda;
            }

            number_of_realizations_with_ALP_pred = number_of_realizations_with_ALP;
            get_current_time(&mut time_during_ALP);

            if result.d_alp_data.d_rand_flag {
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_ref()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                loop_break_flag =
                    rand_i >= rand_all.d_preliminary_realizations_numbers_ALP.len() as i64 - 1;
                if loop_break_flag {
                    break;
                }
            }

            if time_during_ALP - time_before1 >= result.d_alp_data.d_max_time * 0.25
                || number_of_realizations_with_ALP >= realizations_number2
                || tmp_lambda <= 1.0
            {
                if !result.d_alp_data.d_rand_flag {
                    break;
                }
            }

            if time_during_ALP <= time_before_ALP {
                number_of_realizations_with_ALP = Tmin(
                    realizations_number2,
                    number_of_realizations_with_ALP
                        + result.d_alp_data.d_minimum_realizations_number,
                );
            } else {
                let max_number_of_realizations = (number_of_realizations_with_ALP as f64
                    * (result.d_alp_data.d_max_time * 0.35 - (time_before_ALP - time_before1))
                    / (time_during_ALP - time_before_ALP))
                    .floor() as i64;
                number_of_realizations_with_ALP = Tmin(
                    realizations_number2,
                    (0.5 * number_of_realizations_with_ALP as f64
                        + 0.5 * max_number_of_realizations as f64)
                        .floor() as i64,
                );
                if number_of_realizations_with_ALP >= max_number_of_realizations {
                    number_of_realizations_with_ALP = Tmin(
                        realizations_number2,
                        number_of_realizations_with_ALP
                            + result.d_alp_data.d_minimum_realizations_number,
                    );
                }

                if ((number_of_realizations_with_ALP - number_of_realizations_with_ALP_pred) as f64
                    / number_of_realizations_with_ALP_pred as f64)
                    < 0.005
                {
                    number_of_realizations_with_ALP = number_of_realizations_with_ALP_pred;
                    if !result.d_alp_data.d_rand_flag {
                        break;
                    }
                }
            }
        }

        realizations_number2 = number_of_realizations_with_ALP;
        let realizations_number2_lambda = number_of_realizations_with_ALP;

        rand_i = -1;
        let mut time_before_kill = 0.0;
        let mut time_during_kill = 0.0;
        let mut number_of_realizations_with_killing = Tmin(
            realizations_number2,
            result.d_alp_data.d_minimum_realizations_number - 1,
        );
        let mut number_of_realizations_with_killing_pred;
        get_current_time(&mut time_before_kill);

        loop {
            let mut check_time_flag = false;
            if result.d_alp_data.d_rand_flag {
                check_time_flag = false;
                rand_i += 1;
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_ref()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                number_of_realizations_with_killing =
                    rand_all.d_preliminary_realizations_numbers_killing[rand_i as usize];
            }

            result.kill(
                check_time_flag,
                0,
                number_of_realizations_with_killing,
                M_min,
                lambda,
                eps_K,
                &mut K_C,
                &mut K_C_error,
                &mut level,
                &mut diff_opt,
            )?;

            if !result.d_alp_data.d_rand_flag {
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_mut()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                rand_all
                    .d_preliminary_realizations_numbers_killing
                    .push(number_of_realizations_with_killing);
            }

            number_of_realizations_with_killing_pred = number_of_realizations_with_killing;
            get_current_time(&mut time_during_kill);

            let mut tmp_K = 2.0;
            if K_C > 0.0 {
                tmp_K = (K_C_error / K_C) / result.d_alp_data.d_eps_K;
            }

            if result.d_alp_data.d_rand_flag {
                let rand_all = result
                    .d_alp_data
                    .d_rand_all
                    .as_ref()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                loop_break_flag =
                    rand_i >= rand_all.d_preliminary_realizations_numbers_killing.len() as i64 - 1;
                if loop_break_flag {
                    break;
                }
            }

            if time_during_kill - time_before1 >= result.d_alp_data.d_max_time
                || number_of_realizations_with_killing >= realizations_number2
                || tmp_K <= 1.0
            {
                if !result.d_alp_data.d_rand_flag {
                    break;
                }
            }

            if time_during_kill <= time_before_kill {
                number_of_realizations_with_killing = Tmin(
                    realizations_number2,
                    number_of_realizations_with_killing
                        + result.d_alp_data.d_minimum_realizations_number,
                );
            } else {
                let max_number_of_realizations = (number_of_realizations_with_killing as f64
                    * (result.d_alp_data.d_max_time - (time_before_kill - time_before1))
                    / (time_during_kill - time_before_kill))
                    .floor() as i64;
                number_of_realizations_with_killing = Tmin(
                    realizations_number2,
                    (0.5 * number_of_realizations_with_killing as f64
                        + 0.5 * max_number_of_realizations as f64)
                        .floor() as i64,
                );
                if number_of_realizations_with_killing >= max_number_of_realizations {
                    number_of_realizations_with_killing = Tmin(
                        realizations_number2,
                        number_of_realizations_with_killing
                            + result.d_alp_data.d_minimum_realizations_number,
                    );
                }

                if ((number_of_realizations_with_killing - number_of_realizations_with_killing_pred)
                    as f64
                    / number_of_realizations_with_killing_pred as f64)
                    < 0.005
                {
                    number_of_realizations_with_killing = number_of_realizations_with_killing_pred;
                    if !result.d_alp_data.d_rand_flag {
                        break;
                    }
                }
            }
        }

        for k in 0..=number_of_realizations_with_killing {
            if let Some(obj) = result.d_alp_obj[k as usize].as_mut() {
                obj.partially_release_memory();
            }
        }

        let realizations_number2_K = number_of_realizations_with_killing;
        result.d_alp_data.d_realizations_number2 = realizations_number2_lambda;

        if K_C <= 0.0 {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        let mut tmp = ((realizations_number2_K + 1) as f64
            * ((K_C_error / K_C) / result.d_alp_data.d_eps_K)
            * ((K_C_error / K_C) / result.d_alp_data.d_eps_K))
            .ceil();
        tmp = Tmin(tmp as i64, i64::MAX) as f64;
        let realizations_number_killing = tmp as i64;
        tmp = ((realizations_number2_lambda + 1) as f64
            * ((result.d_lambda_tmp_errors.d_elem[nalp as usize]
                / result.d_lambda_tmp.d_elem[nalp as usize])
                / result.d_alp_data.d_eps_lambda)
            * ((result.d_lambda_tmp_errors.d_elem[nalp as usize]
                / result.d_lambda_tmp.d_elem[nalp as usize])
                / result.d_alp_data.d_eps_lambda))
            .ceil();
        tmp = Tmin(tmp as i64, i64::MAX) as f64;
        let realizations_number_lambda = tmp as i64;

        let mut _time_after2 = 0.0;
        get_current_time(&mut _time_after2);
        _time_after2 = Tmax(_time_after2, time_after_tmp);

        result.d_lambda_tmp = array_positive::new(None)?;
        result.d_lambda_tmp_errors = array_positive::new(None)?;
        result.d_C_tmp = array_positive::new(None)?;
        result.d_C_tmp_errors = array_positive::new(None)?;

        let mut j = 1;
        let mut kill_j = 0;
        let mut kill_flag = realizations_number_killing > realizations_number2_K + 1 + j;
        let mut lambda_flag = realizations_number_lambda > realizations_number2_lambda + 1 + j;
        let mut nalp_for_simulation = nalp;

        if result.d_alp_data.d_rand_flag {
            let rand_all = result
                .d_alp_data
                .d_rand_all
                .as_ref()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            lambda_flag =
                rand_all.d_total_realizations_number_with_ALP > realizations_number2_K + j - 1;
            kill_flag =
                rand_all.d_total_realizations_number_with_killing > realizations_number2_K + j - 1;
        }

        if kill_flag || lambda_flag {
            let step_for_time = 1;

            while result.d_n_alp_obj < i64::MAX - realizations_number2_lambda - j {
                kill_flag = realizations_number_killing > realizations_number2_K + j;
                lambda_flag = realizations_number_lambda > realizations_number2_lambda + j;

                if result.d_alp_data.d_rand_flag {
                    let rand_all = result
                        .d_alp_data
                        .d_rand_all
                        .as_ref()
                        .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                    lambda_flag = rand_all.d_total_realizations_number_with_ALP
                        > realizations_number2_K + j - 1;
                    kill_flag = rand_all.d_total_realizations_number_with_killing
                        > realizations_number2_K + j - 1;
                }

                if !(kill_flag || lambda_flag) {
                    break;
                }

                if !kill_flag {
                    nalp_for_simulation = Tmin(nalp_lambda, nalp);
                }

                let mut sucess_flag = false;
                let index = realizations_number2_K + j;
                if index > realizations_number2_lambda {
                    if result.d_alp_obj.len() <= index as usize {
                        result.d_alp_obj.resize_with(index as usize + 1, || None);
                    }
                    result.d_alp_obj[index as usize] = None;
                    result.d_n_alp_obj += 1;
                }

                let mut obj = if index as usize >= result.d_alp_obj.len() {
                    None
                } else {
                    result.d_alp_obj[index as usize].take()
                };
                let mut eps_tmp = 0.0;
                while !sucess_flag {
                    let mut check_time_flag = true;
                    if result.d_alp_data.d_rand_flag {
                        check_time_flag = false;
                    }

                    result.get_single_realization(
                        check_time_flag,
                        M_min,
                        nalp_for_simulation,
                        kill_flag,
                        level,
                        diff_opt,
                        &mut obj,
                        &mut sucess_flag,
                        &mut eps_tmp,
                    )?;

                    if !sucess_flag && index > realizations_number2_lambda {
                        result.d_n_alp_obj -= 1;
                    }
                }

                if result.d_alp_obj.len() <= index as usize {
                    result.d_alp_obj.resize_with(index as usize + 1, || None);
                }
                result.d_alp_obj[index as usize] = obj;

                if index > realizations_number2_lambda && kill_flag {
                    kill_j = j;
                }

                if let Some(obj) = result.d_alp_obj[index as usize].as_mut() {
                    obj.partially_release_memory();
                }

                j += 1;

                if result.d_alp_data.d_memory_size_in_MB > result.d_alp_data.d_max_mem
                    && !result.d_alp_data.d_rand_flag
                {
                    break;
                }

                if j % step_for_time == 0 {
                    let mut time_after3 = 0.0;
                    get_current_time(&mut time_after3);
                    if (time_after3 - time_before1) > result.d_alp_data.d_max_time
                        && !result.d_alp_data.d_rand_flag
                    {
                        break;
                    }
                }
            }
        }

        let final_realizations_number_killing = kill_j + realizations_number2_K + 1;
        let final_realizations_number_lambda =
            Tmax(realizations_number2_lambda + 1, j + realizations_number2_K);
        result.d_n_alp_obj = final_realizations_number_lambda;

        let mut time_after100 = 0.0;
        get_current_time(&mut time_after100);

        if !result.d_alp_data.d_rand_flag {
            let rand_all = result
                .d_alp_data
                .d_rand_all
                .as_mut()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            rand_all.d_total_realizations_number_with_ALP = final_realizations_number_lambda - 1;
            rand_all.d_total_realizations_number_with_killing =
                final_realizations_number_killing - 1;
        }

        if !result.d_alp_data.d_rand_flag && !result.d_alp_data.d_randout.is_empty() {
            let rand_all = result
                .d_alp_data
                .d_rand_all
                .as_ref()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            let mut out = String::new();
            out.push_str(&format!("{}\n", rand_all.d_random_seed));
            out.push_str(&format!(
                "{}\t",
                rand_all
                    .d_first_stage_preliminary_realizations_numbers_ALP
                    .len()
            ));
            for val in &rand_all.d_first_stage_preliminary_realizations_numbers_ALP {
                out.push_str(&format!("{}\t", val));
            }
            out.push('\n');
            out.push_str(&format!(
                "{}\t",
                rand_all.d_preliminary_realizations_numbers_ALP.len()
            ));
            for val in &rand_all.d_preliminary_realizations_numbers_ALP {
                out.push_str(&format!("{}\t", val));
            }
            out.push('\n');
            out.push_str(&format!(
                "{}\t",
                rand_all.d_preliminary_realizations_numbers_killing.len()
            ));
            for val in &rand_all.d_preliminary_realizations_numbers_killing {
                out.push_str(&format!("{}\t", val));
            }
            out.push('\n');
            out.push_str(&format!(
                "{}\n{}\n",
                rand_all.d_total_realizations_number_with_ALP,
                rand_all.d_total_realizations_number_with_killing
            ));
            fs::write(&result.d_alp_data.d_randout, out).map_err(|_| {
                Error::new(
                    format!(
                        "Error - file {} cannot be created\n",
                        result.d_alp_data.d_randout
                    ),
                    3,
                )
            })?;
        }

        let mut inside_simulation_flag = false;
        result.m_CalcTime = time_after100 - time_before1;
        result.output_main_parameters2m_new(
            nalp_for_simulation,
            level,
            &mut inside_simulation_flag,
            final_realizations_number_lambda,
            final_realizations_number_killing,
        )?;

        Ok(result)
    }

    pub fn round_double(mut val_: f64, digits_: i64) -> f64 {
        for _ in 0..digits_ {
            val_ *= 10.0;
        }
        val_ = round(&val_);
        for _ in 0..digits_ {
            val_ /= 10.0;
        }
        val_
    }

    pub fn relative_error_in_percents(val_: f64, val_error_: f64) -> f64 {
        if val_ == 0.0 {
            return f64::MAX;
        }
        Self::round_double(val_error_ / val_ * 100.0, 1).abs()
    }

    pub fn get_number_of_subsimulations(number_of_realizations_: i64) -> Result<i64, Error> {
        let max_number_of_subsimulations = 20;
        let min_number_of_subsimulations = 3;
        let min_number_of_realizations_for_subsimulation = 2;
        if number_of_realizations_
            < min_number_of_realizations_for_subsimulation * min_number_of_subsimulations
        {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }
        let mut res_subsimulations = (number_of_realizations_ as f64).sqrt().ceil() as i64;
        res_subsimulations = Tmin(res_subsimulations, max_number_of_subsimulations);
        res_subsimulations = Tmax(res_subsimulations, min_number_of_subsimulations);
        Ok(res_subsimulations)
    }

    pub fn error_in_calculate_main_parameters2m(
        C: f64,
        C_error: &mut f64,
        C_mult2: f64,
        C_mult2_error: f64,
    ) {
        if C != 0.0 && C_mult2 != 0.0 {
            *C_error = (C * C_mult2_error / C_mult2).abs();
        } else {
            *C_error = C_mult2_error;
        }
    }

    pub fn lambda_exp(i_: i64, exp_array_: &[f64]) -> Result<f64, Error> {
        if exp_array_[i_ as usize] == -1.0 {
            return Err(Error::new(
                "The program is not able to calculate the parameters; rescaling penalties and scoring matrix might help\n".to_string(),
                3,
            ));
        }
        Ok(exp_array_[i_ as usize])
    }

    #[allow(clippy::too_many_arguments)]
    pub fn get_single_realization(
        &mut self,
        check_time_: bool,
        M_min_: i64,
        nalp_: i64,
        killing_flag_: bool,
        level_: i64,
        diff_opt_: i64,
        obj_: &mut Option<Box<alp>>,
        sucess_flag_: &mut bool,
        d_eps_: &mut f64,
    ) -> Result<(), Error> {
        if obj_.is_none() {
            *obj_ = Some(Box::new(alp::new(self.d_alp_data.clone())?));
            self.d_alp_data.d_memory_size_in_MB += std::mem::size_of::<alp>() as f64 / mb_bytes;
        }
        let obj = obj_.as_mut().unwrap();
        obj.d_single_realiztion_calculation_flag = true;
        obj.d_check_time_flag = check_time_;
        *d_eps_ = Tmin(self.d_alp_data.d_eps_K, self.d_alp_data.d_eps_lambda);
        obj.d_diff_opt = diff_opt_;
        obj.d_sentinels_flag = self.d_alp_data.d_sentinels_flag;
        *sucess_flag_ = true;

        while obj.d_nalp < nalp_ {
            obj.simulate_next_alp()?;
            if !obj.d_success {
                *sucess_flag_ = false;
                *obj_ = None;
                *d_eps_ = self.d_alp_data.d_eps_lambda;
                self.d_alp_data.d_memory_size_in_MB -= std::mem::size_of::<alp>() as f64 / mb_bytes;
                return Ok(());
            }
        }

        if killing_flag_ {
            obj.kill_upto_level(M_min_, level_, None)?;
            if !obj.d_success {
                *sucess_flag_ = false;
                *obj_ = None;
                *d_eps_ = self.d_alp_data.d_eps_K;
                self.d_alp_data.d_memory_size_in_MB -= std::mem::size_of::<alp>() as f64 / mb_bytes;
                return Ok(());
            }
        }
        Ok(())
    }

    pub fn randomize_realizations(
        &mut self,
        final_realizations_number_lambda_: i64,
        final_realizations_number_killing_: i64,
    ) -> Result<(), Error> {
        self.randomize_realizations_ind(0, final_realizations_number_killing_ - 1)?;
        self.randomize_realizations_ind(
            final_realizations_number_killing_,
            final_realizations_number_lambda_ - 1,
        )
    }

    pub fn randomize_realizations_ind(&mut self, ind1_: i64, ind2_: i64) -> Result<(), Error> {
        if ind1_ >= ind2_ {
            return Ok(());
        }
        if ind2_ > self.d_n_alp_obj - 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        let total_number = ind2_ - ind1_ + 1;
        let perm = self.generate_random_permulation(total_number)?;
        let mut array_ind: Vec<Option<Box<alp>>> = Vec::with_capacity(total_number as usize);
        for i in 0..total_number {
            array_ind.push(self.d_alp_obj[(ind1_ + perm[i as usize]) as usize].take());
        }
        for i in 0..total_number {
            self.d_alp_obj[(ind1_ + i) as usize] = array_ind[i as usize].take();
        }
        Ok(())
    }

    pub fn generate_random_permulation(&self, dim_: i64) -> Result<Vec<i64>, Error> {
        let mut perm = vec![0; dim_ as usize];
        for i in 0..dim_ {
            perm[i as usize] = i;
        }
        for i in 0..dim_ - 1 {
            let ind_swap = i + alp_data::random_long(self.d_alp_data.ran2(), dim_ - i)?;
            let tmp = perm[ind_swap as usize];
            perm[ind_swap as usize] = perm[i as usize];
            perm[i as usize] = tmp;
        }
        Ok(perm)
    }

    pub fn symmetric_parameters_for_symmetric_scheme(&mut self) {
        let mut symmetric_flag = true;
        for i in 0..self.d_alp_data.d_number_of_AA {
            for j in 0..i {
                if self.d_alp_data.d_smatr[i as usize][j as usize]
                    != self.d_alp_data.d_smatr[j as usize][i as usize]
                {
                    symmetric_flag = false;
                    break;
                }
            }
            if !symmetric_flag {
                break;
            }
        }
        if symmetric_flag {
            for i in 0..self.d_alp_data.d_number_of_AA {
                if self.d_alp_data.d_RR1[i as usize] != self.d_alp_data.d_RR2[i as usize] {
                    symmetric_flag = false;
                    break;
                }
            }
        }
        if symmetric_flag
            && (self.d_alp_data.d_epen1 != self.d_alp_data.d_epen2
                || self.d_alp_data.d_open1 != self.d_alp_data.d_open2)
        {
            symmetric_flag = false;
        }
        if symmetric_flag {
            self.m_AI = 0.5 * (self.m_AI + self.m_AJ);
            self.m_AJ = self.m_AI;
            self.m_AIError = 0.5 * (self.m_AIError + self.m_AJError);
            self.m_AJError = self.m_AIError;
            self.m_AlphaI = 0.5 * (self.m_AlphaI + self.m_AlphaJ);
            self.m_AlphaJ = self.m_AlphaI;
            self.m_AlphaIError = 0.5 * (self.m_AlphaIError + self.m_AlphaJError);
            self.m_AlphaJError = self.m_AlphaIError;
        }
    }

    pub fn get_and_allocate_alp_distribution(
        &self,
        ind1_: i64,
        ind2_: i64,
        alp_distr: &mut Vec<array_positive<f64>>,
        alp_distr_errors: &mut Vec<array_positive<f64>>,
        nalp: i64,
    ) -> Result<(), Error> {
        if nalp <= 0 {
            if nalp < 0 {
                return Err(Error::new("Unexpected error\n".to_string(), 4));
            }
            alp_distr.clear();
            alp_distr_errors.clear();
            return Ok(());
        }

        if alp_distr.len() < nalp as usize + 1 {
            alp_distr.resize_with(nalp as usize + 1, || array_positive::new(None).unwrap());
            alp_distr_errors.resize_with(nalp as usize + 1, || array_positive::new(None).unwrap());
        }
        alp_distr[nalp as usize] = array_positive::new(None)?;
        alp_distr_errors[nalp as usize] = array_positive::new(None)?;

        for i in ind1_..=ind2_ {
            let alp_obj_tmp = self.d_alp_obj[i as usize]
                .as_ref()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            let k = nalp;
            let alp_tmp = alp_obj_tmp.d_alp.d_elem[k as usize];
            let weight_tmp = alp_obj_tmp.d_alp_weights.d_elem[k as usize];
            alp_distr[k as usize].increase_elem_by_x(alp_tmp, weight_tmp);
            alp_distr_errors[k as usize].increase_elem_by_x(alp_tmp, weight_tmp * weight_tmp);
        }

        let ind_diff = (ind2_ - ind1_ + 1) as f64;
        let k = nalp as usize;
        for j in 0..=alp_distr[k].d_dim {
            alp_distr[k].d_elem[j as usize] /= ind_diff;
            alp_distr_errors[k].d_elem[j as usize] /= ind_diff;
            alp_distr_errors[k].d_elem[j as usize] -=
                alp_distr[k].d_elem[j as usize] * alp_distr[k].d_elem[j as usize];
            alp_distr_errors[k].d_elem[j as usize] /= ind_diff;
        }
        Ok(())
    }

    pub fn check_K_criterion(
        &self,
        nalp_: i64,
        ind1_: i64,
        ind2_: i64,
        lambda_: f64,
        eps_K_: f64,
        M_min_: &mut i64,
    ) -> Result<bool, Error> {
        if nalp_ <= 0 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        let mut diff = array_positive::new(None)?;
        let mut sum_of_weights = 0.0;
        let mut M_aver = 0.0;
        for i in ind1_..=ind2_ {
            let alp_obj_tmp = self.d_alp_obj[i as usize]
                .as_ref()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            let alp_tmp = alp_obj_tmp.d_alp.d_elem[nalp_ as usize];
            let weight_tmp = alp_obj_tmp.d_alp_weights.d_elem[nalp_ as usize];
            sum_of_weights += weight_tmp;
            M_aver += alp_tmp as f64 * weight_tmp;
            let cells_counts = &alp_obj_tmp.d_cells_counts;
            for k in cells_counts.d_ind0..=Tmin(alp_tmp, cells_counts.d_dim_plus_d_ind0) {
                diff.increase_elem_by_x(
                    alp_tmp - k,
                    cells_counts.d_elem[(k - cells_counts.d_ind0) as usize] as f64 * weight_tmp,
                );
            }
        }
        let mut den = 0.0;
        for i in 0..=diff.d_dim {
            den += (-lambda_ * i as f64).exp() * diff.d_elem[i as usize];
        }
        if den <= 0.0 || sum_of_weights <= 0.0 {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }
        M_aver /= sum_of_weights;
        let delta_val = den * eps_K_ * (1.0 - (-lambda_).exp());
        let mut diff_opt = 1;
        for i in (0..=diff.d_dim).rev() {
            if (-lambda_ * i as f64).exp() * diff.d_elem[i as usize] > delta_val {
                diff_opt = i + 1;
                break;
            }
        }
        *M_min_ = round(&M_aver) as i64;
        if M_aver < diff_opt as f64 {
            return Ok(false);
        }
        Ok(true)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn check_K_criterion_during_killing(
        &self,
        ind1_: i64,
        ind2_: i64,
        lambda_: f64,
        eps_K_: f64,
        current_level_: i64,
        recommended_level_: &mut i64,
        diff_opt_: &mut i64,
        K_C_: &mut f64,
        K_C_error_: &mut f64,
    ) -> Result<bool, Error> {
        if ind1_ > ind2_ {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        let mut diff = array_positive::new(None)?;
        let mut diff_error = array_positive::new(None)?;
        let mut sum_of_weights = 0.0;
        let mut sum_of_weights_error = 0.0;
        let mut M_aver = 0.0;

        for i in ind1_..=ind2_ {
            let alp_obj_tmp = self.d_alp_obj[i as usize]
                .as_ref()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            let alp_tmp = alp_obj_tmp.d_M;
            let weight_tmp = alp_obj_tmp.d_alp_weights.d_elem[alp_obj_tmp.d_nalp_killing as usize];
            sum_of_weights += weight_tmp;
            sum_of_weights_error += weight_tmp * weight_tmp;
            M_aver += alp_tmp as f64 * weight_tmp;
            let cells_counts = &alp_obj_tmp.d_cells_counts;
            for k in cells_counts.d_ind0..=Tmin(alp_tmp, cells_counts.d_dim_plus_d_ind0) {
                let tmp =
                    cells_counts.d_elem[(k - cells_counts.d_ind0) as usize] as f64 * weight_tmp;
                diff.increase_elem_by_x(alp_tmp - k, tmp);
                diff_error.increase_elem_by_x(alp_tmp - k, tmp * tmp);
            }
        }

        let tmp2 = (ind2_ - ind1_ + 1) as f64;
        sum_of_weights /= tmp2;
        sum_of_weights_error /= tmp2;
        sum_of_weights_error -= sum_of_weights * sum_of_weights;
        sum_of_weights_error /= tmp2;
        sum_of_weights_error = alp_reg::sqrt_for_errors(sum_of_weights_error);

        for i in 0..=diff.d_dim {
            diff.d_elem[i as usize] /= tmp2;
            diff_error.d_elem[i as usize] /= tmp2;
            diff_error.d_elem[i as usize] -= diff.d_elem[i as usize] * diff.d_elem[i as usize];
            diff_error.d_elem[i as usize] /= tmp2;
        }

        let mut den = 0.0;
        let mut den_error = 0.0;
        for i in 0..=diff.d_dim {
            let tmp = (-lambda_ * i as f64).exp();
            den += tmp * diff.d_elem[i as usize];
            den_error += tmp * tmp * diff_error.d_elem[i as usize];
        }
        den_error = alp_reg::sqrt_for_errors(den_error);
        if den <= 0.0 || sum_of_weights <= 0.0 {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        *K_C_ = sum_of_weights / den;
        *K_C_error_ =
            alp_data::error_of_the_ratio(sum_of_weights, sum_of_weights_error, den, den_error);
        M_aver /= tmp2;
        M_aver /= sum_of_weights;
        let delta_val = den * eps_K_ * (1.0 - (-lambda_).exp());
        let mut diff_opt = 1;
        for i in (0..=diff.d_dim).rev() {
            if (-lambda_ * i as f64).exp() * diff.d_elem[i as usize] > delta_val {
                diff_opt = i + 1;
                break;
            }
        }
        if M_aver - (diff_opt as f64) < current_level_ as f64 {
            *recommended_level_ = (M_aver - diff_opt as f64 * 1.1).floor() as i64;
            *diff_opt_ = (M_aver - *recommended_level_ as f64).ceil() as i64;
            return Ok(false);
        }
        *recommended_level_ = current_level_;
        *diff_opt_ = (M_aver - *recommended_level_ as f64).ceil() as i64;
        Ok(true)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn the_criterion(
        &mut self,
        upto_nalp_: i64,
        nalp_for_lambda_simulation_: &mut i64,
        ind1_: i64,
        ind2_: i64,
        alp_distr: &mut Vec<array_positive<f64>>,
        alp_distr_errors: &mut Vec<array_positive<f64>>,
        M_min_: &mut i64,
        M_min_flag_: &mut bool,
        nalp_flag_: &mut bool,
        inside_simulation_flag_: &mut bool,
        C_calculation_: bool,
        lambda_: Option<&mut f64>,
        lambda_error_: Option<&mut f64>,
    ) -> Result<bool, Error> {
        *nalp_flag_ = false;
        *M_min_flag_ = false;

        if ind1_ > ind2_ {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        let mut lambda = 0.0;
        let mut lambda_error = 0.0;
        let mut test_difference = 0.0;
        let mut test_difference_error = 0.0;
        let nalp = upto_nalp_;

        if nalp < 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        self.get_and_allocate_alp_distribution(ind1_, ind2_, alp_distr, alp_distr_errors, nalp)?;

        self.calculate_lambda(
            true,
            upto_nalp_,
            nalp_for_lambda_simulation_,
            inside_simulation_flag_,
            alp_distr,
            alp_distr_errors,
            &mut lambda,
            &mut lambda_error,
            &mut test_difference,
            &mut test_difference_error,
        )?;

        if !*inside_simulation_flag_ {
            return Ok(false);
        }

        self.d_lambda_tmp.set_elem(upto_nalp_, lambda);
        self.d_lambda_tmp_errors.set_elem(upto_nalp_, lambda_error);

        if C_calculation_ {
            let mut C = 0.0;
            let mut C_error = 0.0;
            let mut Sc = 0.0;
            let mut Sc_error = 0.0;
            Self::calculate_C(
                0,
                upto_nalp_,
                alp_distr,
                alp_distr_errors,
                lambda,
                lambda_error,
                &mut C,
                &mut C_error,
                &mut Sc,
                &mut Sc_error,
            )?;
            self.d_C_tmp.set_elem(upto_nalp_, C);
            self.d_C_tmp_errors.set_elem(upto_nalp_, C_error);
        }

        if let Some(lambda_ref) = lambda_ {
            *lambda_ref = lambda;
        }
        if let Some(lambda_error_ref) = lambda_error_ {
            *lambda_error_ref = lambda_error;
        }

        if nalp >= 1 && test_difference <= test_difference_error {
            *nalp_flag_ = true;
            *M_min_ = 0;
            return Ok(true);
        }

        Ok(false)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn calculate_lambda(
        &mut self,
        check_the_criteria_: bool,
        nalp_: i64,
        nalp_thr_: &mut i64,
        inside_simulation_flag_: &mut bool,
        alp_distr: &[array_positive<f64>],
        alp_distr_errors: &[array_positive<f64>],
        lambda_: &mut f64,
        lambda_error_: &mut f64,
        test_difference_: &mut f64,
        test_difference_error_: &mut f64,
    ) -> Result<(), Error> {
        let nalp = nalp_;
        if nalp <= 0 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        let is = self
            .d_alp_data
            .d_is
            .as_ref()
            .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
        let mut tmp_struct = struct_for_lambda_calculation {
            d_alp_distr: alp_distr.to_vec(),
            d_alp_distr_errors: alp_distr_errors.to_vec(),
            d_nalp: nalp,
            d_calculate_alp_number: false,
            ..Default::default()
        };

        let a = 0.0;
        let b = is.d_lambda * 2.0;
        let n_partition = 30;
        let eps = 1e-10;
        let h = (b - a) / n_partition as f64;
        let mut res = Vec::new();
        let mut x2 = 0.0;
        let mut intervals = Vec::new();
        for i in 0..n_partition {
            let x1 = if i == 0 {
                let x1 = Self::function_for_lambda_calculation(a + i as f64 * h, &mut tmp_struct)?;
                if x1.abs() < eps {
                    res.push(a + i as f64 * h);
                }
                x1
            } else {
                x2
            };
            x2 = Self::function_for_lambda_calculation(a + (i + 1) as f64 * h, &mut tmp_struct)?;
            if x2.abs() < eps {
                res.push(a + (i + 1) as f64 * h);
            }
            if x1 * x2 < 0.0 && x1.abs() >= eps && x2.abs() >= eps {
                intervals.push(i);
            }
        }
        for interval in intervals {
            let mut left = a + interval as f64 * h;
            let mut right = a + (interval + 1) as f64 * h;
            let mut f_left = Self::function_for_lambda_calculation(left, &mut tmp_struct)?;
            while (right - left).abs() > eps {
                let mid = (left + right) / 2.0;
                let f_mid = Self::function_for_lambda_calculation(mid, &mut tmp_struct)?;
                if f_mid.abs() < eps {
                    left = mid;
                    right = mid;
                    break;
                }
                if f_left * f_mid <= 0.0 {
                    right = mid;
                } else {
                    left = mid;
                    f_left = f_mid;
                }
            }
            res.push((left + right) / 2.0);
        }

        *inside_simulation_flag_ = true;
        if res.is_empty() {
            *inside_simulation_flag_ = false;
            return Ok(());
        }
        *lambda_ = Self::get_root(&res, is.d_lambda)?;

        tmp_struct.d_calculate_alp_number = true;
        let f1 = Self::function_for_lambda_calculation(*lambda_, &mut tmp_struct)?;
        *nalp_thr_ = tmp_struct.d_alp_number;
        tmp_struct.d_calculate_alp_number = false;
        let slope_error = tmp_struct.d_f_error;
        let sum1 = tmp_struct.d_last_sum;
        let sum1_error = tmp_struct.d_last_sum_error;
        let delta_lambda = *lambda_ / 100.0;
        let f2 = Self::function_for_lambda_calculation(*lambda_ + delta_lambda, &mut tmp_struct)?;
        if delta_lambda == 0.0 || f1 == f2 {
            *lambda_error_ = 0.0;
        } else {
            let derivative = (f2 - f1) / delta_lambda;
            *lambda_error_ = (slope_error / derivative).abs();
        }
        if !check_the_criteria_ {
            return Ok(());
        }

        if nalp > 1 {
            Self::function_for_lambda_calculation(
                self.d_lambda_tmp.d_elem[(nalp - 1) as usize],
                &mut tmp_struct,
            )?;
        } else {
            Self::function_for_lambda_calculation(is.d_ungap_lambda, &mut tmp_struct)?;
        }
        let sum2 = tmp_struct.d_last_sum;
        let sum2_error = tmp_struct.d_last_sum_error;
        let max_sum = Tmax(sum1.abs(), sum2.abs());
        if max_sum != 0.0 {
            *test_difference_ = ((sum1 - sum2) / max_sum).abs();
            *test_difference_error_ = 0.5 * (sum1_error + sum2_error) / max_sum;
        } else {
            *test_difference_ = -1.0;
            *test_difference_error_ = 0.0;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn calculate_C(
        starting_point: i64,
        nalp_: i64,
        alp_distr: &[array_positive<f64>],
        alp_distr_errors: &[array_positive<f64>],
        lambda_: f64,
        lambda_error_: f64,
        C_: &mut f64,
        C_error_: &mut f64,
        Sc_: &mut f64,
        Sc_error_: &mut f64,
    ) -> Result<(), Error> {
        let total_number_of_ALP = nalp_;
        if total_number_of_ALP < 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        let mut P = vec![0.0; total_number_of_ALP as usize + 1];
        let mut P_errors = vec![0.0; total_number_of_ALP as usize + 1];
        P[0] = 1.0;
        P_errors[0] = 0.0;
        for j in 1..=total_number_of_ALP {
            let tmp = &alp_distr[j as usize];
            let tmp_errors = &alp_distr_errors[j as usize];
            for i in 0..=tmp.d_dim {
                P[j as usize] += tmp.d_elem[i as usize];
                P_errors[j as usize] += tmp_errors.d_elem[i as usize];
            }
            P_errors[j as usize] = alp_reg::sqrt_for_errors(P_errors[j as usize]);
        }

        let mut values_P_ratio = vec![0.0; total_number_of_ALP as usize];
        let mut errors_P_ratio = vec![0.0; total_number_of_ALP as usize];
        for j in 0..total_number_of_ALP {
            values_P_ratio[j as usize] = P[(j + 1) as usize] / P[j as usize];
            errors_P_ratio[j as usize] = alp_data::error_of_the_ratio(
                P[(j + 1) as usize],
                P_errors[(j + 1) as usize],
                P[j as usize],
                P_errors[j as usize],
            );
        }

        let beta1 = 0.0;
        let beta1_error = 0.0;
        let number_of_elements = total_number_of_ALP;
        let cut_left_tail = true;
        let cut_right_tail = false;
        let y = 2.0;
        let mut k1_opt = 0;
        let mut k2_opt = 0;
        let mut P_beta_inf = 0.0;
        let mut P_beta_inf_error = 0.0;
        let mut res_was_calculated = false;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            0,
            number_of_elements - starting_point,
            &values_P_ratio[starting_point as usize..],
            &mut errors_P_ratio[starting_point as usize..],
            cut_left_tail,
            cut_right_tail,
            y,
            &mut P_beta_inf,
            beta1,
            &mut P_beta_inf_error,
            beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }
        P_beta_inf = 1.0 - P_beta_inf;

        let mut E = vec![0.0; total_number_of_ALP as usize + 1];
        let mut E_errors = vec![0.0; total_number_of_ALP as usize + 1];
        let mut E_T_beta = vec![0.0; total_number_of_ALP as usize + 1];
        let mut E_T_beta_errors = vec![0.0; total_number_of_ALP as usize + 1];
        E[0] = 1.0;
        E_T_beta[0] = 0.0;
        E_errors[0] = 0.0;
        E_T_beta_errors[0] = 0.0;
        for j in 1..=total_number_of_ALP {
            let tmp = &alp_distr[j as usize];
            let tmp_errors = &alp_distr_errors[j as usize];
            for i in 0..=tmp.d_dim {
                let tmp_double = (lambda_ * i as f64).exp();
                E[j as usize] += tmp_double * tmp.d_elem[i as usize];
                E_errors[j as usize] += tmp_double * tmp_double * tmp_errors.d_elem[i as usize];
                let tmp_double = i as f64 * (lambda_ * i as f64).exp();
                E_T_beta[j as usize] += tmp_double * tmp.d_elem[i as usize];
                E_T_beta_errors[j as usize] +=
                    tmp_double * tmp_double * tmp_errors.d_elem[i as usize];
            }
            E_errors[j as usize] = alp_reg::sqrt_for_errors(E_errors[j as usize]);
            E_T_beta_errors[j as usize] = alp_reg::sqrt_for_errors(E_T_beta_errors[j as usize]);
        }

        let E_aver;
        let E_aver_error;
        let E_T_beta_diff_aver;
        let E_T_beta_diff_aver_error;
        if total_number_of_ALP == 1 {
            E_aver = E[1];
            E_aver_error = E_errors[1];
            E_T_beta_diff_aver = E_T_beta[1] - E_T_beta[0];
            E_T_beta_diff_aver_error = E_T_beta_errors[1];
        } else {
            let mut beta0 = 0.0;
            let mut beta1 = 0.0;
            let mut beta0_error = 0.0;
            let mut beta1_error = 0.0;
            res_was_calculated = false;
            let mut E_errors_tail = E_errors[1 + starting_point as usize..].to_vec();
            alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                total_number_of_ALP - starting_point,
                &E[1 + starting_point as usize..],
                &mut E_errors_tail,
                cut_left_tail,
                cut_right_tail,
                y,
                &mut beta0,
                beta1,
                &mut beta0_error,
                beta1_error,
                &mut k1_opt,
                &mut k2_opt,
                &mut res_was_calculated,
            )?;
            if !res_was_calculated {
                return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
            }
            E_aver = beta0;
            E_aver_error = beta0_error;

            res_was_calculated = false;
            let mut E_T_beta_errors_tail = E_T_beta_errors[1 + starting_point as usize..].to_vec();
            alp_reg::robust_regression_sum_with_cut_LSM(
                0,
                total_number_of_ALP - starting_point,
                &E_T_beta[1 + starting_point as usize..],
                &mut E_T_beta_errors_tail,
                cut_left_tail,
                cut_right_tail,
                y,
                &mut beta0,
                &mut beta1,
                &mut beta0_error,
                &mut beta1_error,
                &mut k1_opt,
                &mut k2_opt,
                &mut res_was_calculated,
            )?;
            if !res_was_calculated {
                return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
            }
            E_T_beta_diff_aver = beta1;
            E_T_beta_diff_aver_error = beta1_error;
        }

        let exp_lambda_error = (-lambda_).exp() * lambda_error_;
        let exp_lambda = 1.0 - (-lambda_).exp();
        let den_error = alp_data::error_of_the_product(
            E_T_beta_diff_aver,
            E_T_beta_diff_aver_error,
            exp_lambda,
            exp_lambda_error,
        );
        let den = (1.0 - (-lambda_).exp()) * E_T_beta_diff_aver;
        let (nom, nom_error) = if calculate_C_S_constant_flag {
            *Sc_error_ = E_aver_error;
            *Sc_ = E_aver;
            (
                P_beta_inf * E_aver,
                alp_data::error_of_the_product(P_beta_inf, P_beta_inf_error, E_aver, E_aver_error),
            )
        } else {
            let E_aver_sqr_error =
                alp_data::error_of_the_product(E_aver, E_aver_error, E_aver, E_aver_error);
            let E_aver_sqr = E_aver * E_aver;
            (
                P_beta_inf * E_aver_sqr,
                alp_data::error_of_the_product(
                    P_beta_inf,
                    P_beta_inf_error,
                    E_aver_sqr,
                    E_aver_sqr_error,
                ),
            )
        };
        *C_error_ = alp_data::error_of_the_ratio(nom, nom_error, den, den_error);
        *C_ = nom / den;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn calculate_FSC(
        &self,
        nalp_: i64,
        ind1_: i64,
        ind2_: i64,
        alp_distr: &[array_positive<f64>],
        lambda_: f64,
        Sc_: f64,
        a_I_: &mut f64,
        a_I_error_: &mut f64,
        a_J_: &mut f64,
        a_J_error_: &mut f64,
        sigma_: &mut f64,
        sigma_error_: &mut f64,
        alpha_I_: &mut f64,
        alpha_I_error_: &mut f64,
        alpha_J_: &mut f64,
        alpha_J_error_: &mut f64,
    ) -> Result<(), Error> {
        if nalp_ < 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        let tmp = &alp_distr[nalp_ as usize];
        let dim = tmp.d_dim;
        let dbl_max_log = f64::MAX.ln();
        let mut exp_array = vec![0.0; dim as usize + 1];
        for i in 0..=dim {
            let tmp = i as f64 * lambda_;
            if tmp < dbl_max_log {
                exp_array[i as usize] = tmp.exp();
            } else {
                exp_array[i as usize] = -1.0;
            }
        }

        let n = nalp_ as usize;
        let mut delta_E = vec![0.0; n];
        let mut delta_E_error = vec![0.0; n];
        let mut delta_E_E = vec![0.0; n];
        let mut delta_E_E_error = vec![0.0; n];
        let mut delta_I = vec![0.0; n];
        let mut delta_I_error = vec![0.0; n];
        let mut delta_J = vec![0.0; n];
        let mut delta_J_error = vec![0.0; n];
        let mut delta_I_I = vec![0.0; n];
        let mut delta_I_I_error = vec![0.0; n];
        let mut delta_I_J = vec![0.0; n];
        let mut delta_I_J_error = vec![0.0; n];
        let mut delta_J_J = vec![0.0; n];
        let mut delta_J_J_error = vec![0.0; n];
        let mut cov_J_J = vec![0.0; n];
        let mut cov_J_J_error = vec![0.0; n];
        let mut cov_I_J = vec![0.0; n];
        let mut cov_I_J_error = vec![0.0; n];
        let mut cov_I_I = vec![0.0; n];
        let mut cov_I_I_error = vec![0.0; n];
        let mut cov_E_E = vec![0.0; n];
        let mut cov_E_E_error = vec![0.0; n];

        let mut C_S_constant = 1.0;
        if calculate_C_S_constant_flag && Sc_ > 0.0 {
            C_S_constant = Sc_;
        }
        let one_div_C_S_constant = 1.0 / C_S_constant;

        for i in ind1_..=ind2_ {
            let alp_obj_tmp = self.d_alp_obj[i as usize]
                .as_ref()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            for j in 1..=nalp_ {
                let j_1 = j - 1;
                let E_j_1 = alp_obj_tmp.d_alp.d_elem[j_1 as usize];
                let E_j = alp_obj_tmp.d_alp.d_elem[j as usize];
                let weight_j = alp_obj_tmp.d_alp_weights.d_elem[j as usize];
                let I_j_1 = alp_obj_tmp.d_H_I.d_elem[j_1 as usize];
                let I_j = alp_obj_tmp.d_H_I.d_elem[j as usize];
                let J_j_1 = alp_obj_tmp.d_H_J.d_elem[j_1 as usize];
                let J_j = alp_obj_tmp.d_H_J.d_elem[j as usize];
                let exp_tmp = Self::lambda_exp(E_j, &exp_array)? * one_div_C_S_constant;
                let delta_I_tmp = (I_j - I_j_1) as f64 * exp_tmp * weight_j;
                let delta_J_tmp = (J_j - J_j_1) as f64 * exp_tmp * weight_j;
                let delta_E_tmp = (E_j - E_j_1) as f64 * exp_tmp * weight_j;
                let delta_E_E_tmp =
                    (E_j - E_j_1) as f64 * (E_j - E_j_1) as f64 * exp_tmp * weight_j;
                let delta_I_I_tmp = delta_I_tmp * (I_j - I_j_1) as f64;
                let delta_J_J_tmp = delta_J_tmp * (J_j - J_j_1) as f64;
                let delta_I_J_tmp = delta_I_tmp * (J_j - J_j_1) as f64;
                let idx = j_1 as usize;
                delta_E[idx] += delta_E_tmp;
                delta_E_error[idx] += delta_E_tmp * delta_E_tmp;
                delta_E_E[idx] += delta_E_E_tmp;
                delta_E_E_error[idx] += delta_E_E_tmp * delta_E_E_tmp;
                delta_I[idx] += delta_I_tmp;
                delta_I_error[idx] += delta_I_tmp * delta_I_tmp;
                delta_J[idx] += delta_J_tmp;
                delta_J_error[idx] += delta_J_tmp * delta_J_tmp;
                delta_I_I[idx] += delta_I_I_tmp;
                delta_I_I_error[idx] += delta_I_I_tmp * delta_I_I_tmp;
                delta_I_J[idx] += delta_I_J_tmp;
                delta_I_J_error[idx] += delta_I_J_tmp * delta_I_J_tmp;
                delta_J_J[idx] += delta_J_J_tmp;
                delta_J_J_error[idx] += delta_J_J_tmp * delta_J_J_tmp;
            }
        }

        let ind_diff = (ind2_ - ind1_ + 1) as f64;
        for j in 0..n {
            delta_E[j] /= ind_diff;
            delta_E_error[j] /= ind_diff;
            delta_E_error[j] -= delta_E[j] * delta_E[j];
            delta_E_error[j] /= ind_diff;
            delta_E_error[j] = alp_reg::sqrt_for_errors(delta_E_error[j]);
            delta_E_E[j] /= ind_diff;
            delta_E_E_error[j] /= ind_diff;
            delta_E_E_error[j] -= delta_E_E[j] * delta_E_E[j];
            delta_E_E_error[j] /= ind_diff;

            delta_I[j] /= ind_diff;
            delta_I_error[j] /= ind_diff;
            delta_I_error[j] -= delta_I[j] * delta_I[j];
            delta_I_error[j] /= ind_diff;
            delta_I_error[j] = alp_reg::sqrt_for_errors(delta_I_error[j]);
            delta_J[j] /= ind_diff;
            delta_J_error[j] /= ind_diff;
            delta_J_error[j] -= delta_J[j] * delta_J[j];
            delta_J_error[j] /= ind_diff;
            delta_J_error[j] = alp_reg::sqrt_for_errors(delta_J_error[j]);

            delta_I_J[j] /= ind_diff;
            delta_I_J_error[j] /= ind_diff;
            delta_I_J_error[j] -= delta_I_J[j] * delta_I_J[j];
            delta_I_J_error[j] /= ind_diff;
            delta_I_I[j] /= ind_diff;
            delta_I_I_error[j] /= ind_diff;
            delta_I_I_error[j] -= delta_I_I[j] * delta_I_I[j];
            delta_I_I_error[j] /= ind_diff;
            delta_J_J[j] /= ind_diff;
            delta_J_J_error[j] /= ind_diff;
            delta_J_J_error[j] -= delta_J_J[j] * delta_J_J[j];
            delta_J_J_error[j] /= ind_diff;

            cov_I_J[j] = delta_I_J[j] - delta_I[j] * delta_J[j];
            cov_I_I[j] = delta_I_I[j] - delta_I[j] * delta_I[j];
            cov_J_J[j] = delta_J_J[j] - delta_J[j] * delta_J[j];
            cov_E_E[j] = delta_E_E[j] - delta_E[j] * delta_E[j];

            cov_I_J_error[j] = alp_data::error_of_the_product(
                delta_I[j],
                delta_I_error[j],
                delta_J[j],
                delta_J_error[j],
            );
            cov_I_J_error[j] =
                alp_reg::sqrt_for_errors(delta_I_J_error[j] + cov_I_J_error[j] * cov_I_J_error[j]);
            cov_I_I_error[j] = alp_data::error_of_the_product(
                delta_I[j],
                delta_I_error[j],
                delta_I[j],
                delta_I_error[j],
            );
            cov_I_I_error[j] =
                alp_reg::sqrt_for_errors(delta_I_I_error[j] + cov_I_I_error[j] * cov_I_I_error[j]);
            cov_J_J_error[j] = alp_data::error_of_the_product(
                delta_J[j],
                delta_J_error[j],
                delta_J[j],
                delta_J_error[j],
            );
            cov_J_J_error[j] =
                alp_reg::sqrt_for_errors(delta_J_J_error[j] + cov_J_J_error[j] * cov_J_J_error[j]);
            cov_E_E_error[j] = alp_data::error_of_the_product(
                delta_E[j],
                delta_E_error[j],
                delta_E[j],
                delta_E_error[j],
            );
            cov_E_E_error[j] =
                alp_reg::sqrt_for_errors(delta_E_E_error[j] + cov_E_E_error[j] * cov_E_E_error[j]);
        }

        let beta1 = 0.0;
        let beta1_error = 0.0;
        let number_of_elements = nalp_;
        let cut_left_tail = true;
        let cut_right_tail = false;
        let y = 2.0;
        let mut k1_opt = 0;
        let mut k2_opt = 0;
        let mut res_was_calculated = false;

        let mut delta_I_aver = 0.0;
        let mut delta_I_aver_error = 0.0;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            0,
            number_of_elements,
            &delta_I,
            &mut delta_I_error,
            cut_left_tail,
            cut_right_tail,
            y,
            &mut delta_I_aver,
            beta1,
            &mut delta_I_aver_error,
            beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        let mut delta_J_aver = 0.0;
        let mut delta_J_aver_error = 0.0;
        res_was_calculated = false;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            0,
            number_of_elements,
            &delta_J,
            &mut delta_J_error,
            cut_left_tail,
            cut_right_tail,
            y,
            &mut delta_J_aver,
            beta1,
            &mut delta_J_aver_error,
            beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        let mut delta_E_aver = 0.0;
        let mut delta_E_aver_error = 0.0;
        res_was_calculated = false;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            0,
            number_of_elements,
            &delta_E,
            &mut delta_E_error,
            cut_left_tail,
            cut_right_tail,
            y,
            &mut delta_E_aver,
            beta1,
            &mut delta_E_aver_error,
            beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        let mut cov_I_J_aver = 0.0;
        let mut cov_I_J_aver_error = 0.0;
        res_was_calculated = false;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            0,
            number_of_elements,
            &cov_I_J,
            &mut cov_I_J_error,
            cut_left_tail,
            cut_right_tail,
            y,
            &mut cov_I_J_aver,
            beta1,
            &mut cov_I_J_aver_error,
            beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        let mut cov_I_I_aver = 0.0;
        let mut cov_I_I_aver_error = 0.0;
        res_was_calculated = false;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            0,
            number_of_elements,
            &cov_I_I,
            &mut cov_I_I_error,
            cut_left_tail,
            cut_right_tail,
            y,
            &mut cov_I_I_aver,
            beta1,
            &mut cov_I_I_aver_error,
            beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        let mut cov_J_J_aver = 0.0;
        let mut cov_J_J_aver_error = 0.0;
        res_was_calculated = false;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            0,
            number_of_elements,
            &cov_J_J,
            &mut cov_J_J_error,
            cut_left_tail,
            cut_right_tail,
            y,
            &mut cov_J_J_aver,
            beta1,
            &mut cov_J_J_aver_error,
            beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        let mut cov_E_E_aver = 0.0;
        let mut cov_E_E_aver_error = 0.0;
        res_was_calculated = false;
        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
            0,
            number_of_elements,
            &cov_E_E,
            &mut cov_E_E_error,
            cut_left_tail,
            cut_right_tail,
            y,
            &mut cov_E_E_aver,
            beta1,
            &mut cov_E_E_aver_error,
            beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        if delta_E_aver <= 0.0 {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        *a_I_ = delta_I_aver / delta_E_aver;
        *a_I_error_ = alp_data::error_of_the_ratio(
            delta_I_aver,
            delta_I_aver_error,
            delta_E_aver,
            delta_E_aver_error,
        );
        *a_J_ = delta_J_aver / delta_E_aver;
        *a_J_error_ = alp_data::error_of_the_ratio(
            delta_J_aver,
            delta_J_aver_error,
            delta_E_aver,
            delta_E_aver_error,
        );
        Self::sigma_calculation(
            delta_I_aver,
            delta_I_aver_error,
            delta_J_aver,
            delta_J_aver_error,
            delta_E_aver,
            delta_E_aver_error,
            cov_E_E_aver,
            cov_E_E_aver_error,
            cov_I_J_aver,
            cov_I_J_aver_error,
            sigma_,
            sigma_error_,
        );
        Self::sigma_calculation(
            delta_I_aver,
            delta_I_aver_error,
            delta_I_aver,
            delta_I_aver_error,
            delta_E_aver,
            delta_E_aver_error,
            cov_E_E_aver,
            cov_E_E_aver_error,
            cov_I_I_aver,
            cov_I_I_aver_error,
            alpha_I_,
            alpha_I_error_,
        );
        Self::sigma_calculation(
            delta_J_aver,
            delta_J_aver_error,
            delta_J_aver,
            delta_J_aver_error,
            delta_E_aver,
            delta_E_aver_error,
            cov_E_E_aver,
            cov_E_E_aver_error,
            cov_J_J_aver,
            cov_J_J_aver_error,
            alpha_J_,
            alpha_J_error_,
        );

        *a_I_ = Tmax(*a_I_, 0.0);
        *a_J_ = Tmax(*a_J_, 0.0);
        *sigma_ = Tmax(*sigma_, 0.0);
        *alpha_I_ = Tmax(*alpha_I_, 0.0);
        *alpha_J_ = Tmax(*alpha_J_, 0.0);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn calculate_main_parameters2m(
        &mut self,
        final_realizations_number_lambda_: i64,
        final_realizations_number_killing_: i64,
        nalp_for_lambda_simulation: i64,
        level: i64,
        inside_simulation_flag: &mut bool,
        lambda: &mut f64,
        lambda_error: &mut f64,
        test_difference: &mut f64,
        test_difference_error: &mut f64,
        C: &mut f64,
        C_error: &mut f64,
        K_C: &mut f64,
        K_C_error: &mut f64,
        a_I: &mut f64,
        a_I_error: &mut f64,
        a_J: &mut f64,
        a_J_error: &mut f64,
        sigma: &mut f64,
        sigma_error: &mut f64,
        alpha_I: &mut f64,
        alpha_I_error: &mut f64,
        alpha_J: &mut f64,
        alpha_J_error: &mut f64,
        K: &mut f64,
        K_error: &mut f64,
        flag_: &mut bool,
    ) -> Result<(), Error> {
        *flag_ = false;

        if final_realizations_number_killing_ > final_realizations_number_lambda_ {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        let mult_number_lambda = Self::get_number_of_subsimulations(self.d_n_alp_obj)?;
        let mult_number_K = Self::get_number_of_subsimulations(final_realizations_number_killing_)?;
        self.d_mult_number = Tmin(mult_number_lambda, mult_number_K);

        let mut d_mult_realizations = vec![0i64; self.d_mult_number as usize + 1];
        let mut d_mult_K_realizations = vec![0i64; self.d_mult_number as usize + 1];
        let mut lambda_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut lambda_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut C_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut C_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut a_I_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut a_I_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut a_J_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut a_J_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut sigma_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut sigma_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut alpha_I_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut alpha_I_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut alpha_J_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut alpha_J_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut K_C_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut K_C_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut K_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut K_mult_error = vec![0.0; self.d_mult_number as usize + 1];
        let mut Sc_mult = vec![0.0; self.d_mult_number as usize + 1];
        let mut Sc_mult_error = vec![0.0; self.d_mult_number as usize + 1];

        let mut lambda_mult2 = 0.0;
        let mut C_mult2 = 0.0;
        let mut K_C_mult2 = 0.0;
        let mut a_I_mult2 = 0.0;
        let mut a_J_mult2 = 0.0;
        let mut sigma_mult2 = 0.0;
        let mut alpha_I_mult2 = 0.0;
        let mut alpha_J_mult2 = 0.0;
        let mut K_mult2 = 0.0;

        let mut lambda_mult2_error = 0.0;
        let mut C_mult2_error = 0.0;
        let mut K_C_mult2_error = 0.0;
        let mut a_I_mult2_error = 0.0;
        let mut a_J_mult2_error = 0.0;
        let mut sigma_mult2_error = 0.0;
        let mut alpha_I_mult2_error = 0.0;
        let mut alpha_J_mult2_error = 0.0;
        let mut K_mult2_error = 0.0;

        let mut alp_distr: Vec<array_positive<f64>> = Vec::new();
        let mut alp_distr_errors: Vec<array_positive<f64>> = Vec::new();
        for j in 0..=nalp_for_lambda_simulation {
            self.get_and_allocate_alp_distribution(
                0,
                self.d_n_alp_obj - 1,
                &mut alp_distr,
                &mut alp_distr_errors,
                j,
            )?;
        }

        let mut alp_mult_distr: Vec<Vec<array_positive<f64>>> =
            vec![Vec::new(); self.d_mult_number as usize + 1];
        let mut alp_mult_distr_errors: Vec<Vec<array_positive<f64>>> =
            vec![Vec::new(); self.d_mult_number as usize + 1];
        alp_mult_distr[0] = alp_distr.clone();
        alp_mult_distr_errors[0] = alp_distr_errors.clone();

        let mut real_number =
            (final_realizations_number_lambda_ as f64 / self.d_mult_number as f64).floor() as i64;
        d_mult_realizations[0] = final_realizations_number_lambda_;
        for k in 1..=self.d_mult_number {
            d_mult_realizations[k as usize] = real_number;
        }

        let mut nr_tmp = 0;
        for k in 1..=self.d_mult_number {
            nr_tmp += d_mult_realizations[k as usize];
            for j in 0..=nalp_for_lambda_simulation {
                self.get_and_allocate_alp_distribution(
                    nr_tmp - d_mult_realizations[k as usize],
                    nr_tmp - 1,
                    &mut alp_mult_distr[k as usize],
                    &mut alp_mult_distr_errors[k as usize],
                    j,
                )?;
            }
        }

        for k in 1..=self.d_mult_number {
            let mut nalp_tmp = 0;
            let mut test_difference_tmp = 0.0;
            let mut test_difference_error_tmp = 0.0;
            self.calculate_lambda(
                false,
                nalp_for_lambda_simulation,
                &mut nalp_tmp,
                inside_simulation_flag,
                &alp_mult_distr[k as usize],
                &alp_mult_distr_errors[k as usize],
                &mut lambda_mult[k as usize],
                &mut lambda_mult_error[k as usize],
                &mut test_difference_tmp,
                &mut test_difference_error_tmp,
            )?;

            if !*inside_simulation_flag {
                self.symmetric_parameters_for_symmetric_scheme();
                return Ok(());
            }

            lambda_mult2 += lambda_mult[k as usize];
            lambda_mult2_error += lambda_mult[k as usize] * lambda_mult[k as usize];
        }

        let mut nalp_tmp = 0;
        self.calculate_lambda(
            false,
            nalp_for_lambda_simulation,
            &mut nalp_tmp,
            inside_simulation_flag,
            &alp_distr,
            &alp_distr_errors,
            lambda,
            lambda_error,
            test_difference,
            test_difference_error,
        )?;

        if !*inside_simulation_flag {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        lambda_mult[0] = *lambda;
        lambda_mult_error[0] = *lambda_error;

        for k in 1..=self.d_mult_number {
            Self::calculate_C(
                0,
                nalp_for_lambda_simulation,
                &alp_mult_distr[k as usize],
                &alp_mult_distr_errors[k as usize],
                lambda_mult[k as usize],
                lambda_mult_error[k as usize],
                &mut C_mult[k as usize],
                &mut C_mult_error[k as usize],
                &mut Sc_mult[k as usize],
                &mut Sc_mult_error[k as usize],
            )?;

            C_mult2 += C_mult[k as usize];
            C_mult2_error += C_mult[k as usize] * C_mult[k as usize];
        }

        let mut Sc = 0.0;
        let mut Sc_error = 0.0;
        Self::calculate_C(
            0,
            nalp_for_lambda_simulation,
            &alp_distr,
            &alp_distr_errors,
            *lambda,
            *lambda_error,
            C,
            C_error,
            &mut Sc,
            &mut Sc_error,
        )?;

        C_mult[0] = *C;
        C_mult_error[0] = *C_error;
        Sc_mult[0] = Sc;
        Sc_mult_error[0] = Sc_error;

        nr_tmp = 0;
        for k in 1..=self.d_mult_number {
            nr_tmp += d_mult_realizations[k as usize];
            self.calculate_FSC(
                nalp_for_lambda_simulation,
                nr_tmp - d_mult_realizations[k as usize],
                nr_tmp - 1,
                &alp_mult_distr[k as usize],
                lambda_mult[k as usize],
                Sc_mult[k as usize],
                &mut a_I_mult[k as usize],
                &mut a_I_mult_error[k as usize],
                &mut a_J_mult[k as usize],
                &mut a_J_mult_error[k as usize],
                &mut sigma_mult[k as usize],
                &mut sigma_mult_error[k as usize],
                &mut alpha_I_mult[k as usize],
                &mut alpha_I_mult_error[k as usize],
                &mut alpha_J_mult[k as usize],
                &mut alpha_J_mult_error[k as usize],
            )?;

            a_I_mult2 += a_I_mult[k as usize];
            a_I_mult2_error += a_I_mult[k as usize] * a_I_mult[k as usize];
            a_J_mult2 += a_J_mult[k as usize];
            a_J_mult2_error += a_J_mult[k as usize] * a_J_mult[k as usize];
            sigma_mult2 += sigma_mult[k as usize];
            sigma_mult2_error += sigma_mult[k as usize] * sigma_mult[k as usize];
            alpha_I_mult2 += alpha_I_mult[k as usize];
            alpha_I_mult2_error += alpha_I_mult[k as usize] * alpha_I_mult[k as usize];
            alpha_J_mult2 += alpha_J_mult[k as usize];
            alpha_J_mult2_error += alpha_J_mult[k as usize] * alpha_J_mult[k as usize];
        }

        self.calculate_FSC(
            nalp_for_lambda_simulation,
            0,
            final_realizations_number_lambda_ - 1,
            &alp_distr,
            *lambda,
            Sc,
            a_I,
            a_I_error,
            a_J,
            a_J_error,
            sigma,
            sigma_error,
            alpha_I,
            alpha_I_error,
            alpha_J,
            alpha_J_error,
        )?;

        a_I_mult[0] = *a_I;
        a_I_mult_error[0] = *a_I_error;
        a_J_mult[0] = *a_J;
        a_J_mult_error[0] = *a_J_error;
        sigma_mult[0] = *sigma;
        sigma_mult_error[0] = *sigma_error;
        alpha_I_mult[0] = *alpha_I;
        alpha_I_mult_error[0] = *alpha_I_error;
        alpha_J_mult[0] = *alpha_J;
        alpha_J_mult_error[0] = *alpha_J_error;

        real_number =
            (final_realizations_number_killing_ as f64 / self.d_mult_number as f64).floor() as i64;
        d_mult_K_realizations[0] = final_realizations_number_killing_;
        for k in 1..=self.d_mult_number {
            d_mult_K_realizations[k as usize] = real_number;
        }

        nr_tmp = 0;
        for k in 1..=self.d_mult_number {
            nr_tmp += d_mult_K_realizations[k as usize];
            let mut recommended_level = 0;
            let mut diff_opt = 0;
            self.check_K_criterion_during_killing(
                nr_tmp - d_mult_K_realizations[k as usize],
                nr_tmp - 1,
                lambda_mult[k as usize],
                self.d_alp_data.d_eps_K,
                level,
                &mut recommended_level,
                &mut diff_opt,
                &mut K_C_mult[k as usize],
                &mut K_C_mult_error[k as usize],
            )?;

            K_mult[k as usize] = C_mult[k as usize] * K_C_mult[k as usize];
            K_mult_error[k as usize] = alp_data::error_of_the_product(
                C_mult[k as usize],
                C_mult_error[k as usize],
                K_C_mult[k as usize],
                K_C_mult_error[k as usize],
            );

            K_C_mult2 += K_C_mult[k as usize];
            K_C_mult2_error += K_C_mult[k as usize] * K_C_mult[k as usize];
            K_mult2 += K_mult[k as usize];
            K_mult2_error += K_mult[k as usize] * K_mult[k as usize];
        }

        let mut recommended_level = 0;
        let mut diff_opt = 0;
        self.check_K_criterion_during_killing(
            0,
            final_realizations_number_killing_ - 1,
            *lambda,
            self.d_alp_data.d_eps_K,
            level,
            &mut recommended_level,
            &mut diff_opt,
            K_C,
            K_C_error,
        )?;

        *K = *C * *K_C;
        *K_error = alp_data::error_of_the_product(*C, *C_error, *K_C, *K_C_error);

        K_C_mult[0] = *K_C;
        K_C_mult_error[0] = *K_C_error;
        K_mult[0] = *K;
        K_mult_error[0] = *K_error;

        let d_mult_number_f = self.d_mult_number as f64;
        lambda_mult2 /= d_mult_number_f;
        C_mult2 /= d_mult_number_f;
        K_C_mult2 /= d_mult_number_f;
        a_I_mult2 /= d_mult_number_f;
        a_J_mult2 /= d_mult_number_f;
        sigma_mult2 /= d_mult_number_f;
        alpha_I_mult2 /= d_mult_number_f;
        alpha_J_mult2 /= d_mult_number_f;
        K_mult2 /= d_mult_number_f;

        lambda_mult2_error /= d_mult_number_f;
        C_mult2_error /= d_mult_number_f;
        K_C_mult2_error /= d_mult_number_f;
        a_I_mult2_error /= d_mult_number_f;
        a_J_mult2_error /= d_mult_number_f;
        sigma_mult2_error /= d_mult_number_f;
        alpha_I_mult2_error /= d_mult_number_f;
        alpha_J_mult2_error /= d_mult_number_f;
        K_mult2_error /= d_mult_number_f;

        let mult_number_double_lambda =
            final_realizations_number_lambda_ as f64 / real_number as f64;
        let mult_number_double_K = final_realizations_number_killing_ as f64 / real_number as f64;

        lambda_mult2_error =
            alp_reg::sqrt_for_errors(lambda_mult2_error - lambda_mult2 * lambda_mult2)
                / mult_number_double_lambda.sqrt();
        C_mult2_error = alp_reg::sqrt_for_errors(C_mult2_error - C_mult2 * C_mult2)
            / mult_number_double_lambda.sqrt();
        K_C_mult2_error = alp_reg::sqrt_for_errors(K_C_mult2_error - K_C_mult2 * K_C_mult2)
            / mult_number_double_K.sqrt();
        a_I_mult2_error = alp_reg::sqrt_for_errors(a_I_mult2_error - a_I_mult2 * a_I_mult2)
            / mult_number_double_lambda.sqrt();
        a_J_mult2_error = alp_reg::sqrt_for_errors(a_J_mult2_error - a_J_mult2 * a_J_mult2)
            / mult_number_double_lambda.sqrt();
        sigma_mult2_error = alp_reg::sqrt_for_errors(sigma_mult2_error - sigma_mult2 * sigma_mult2)
            / mult_number_double_lambda.sqrt();
        alpha_I_mult2_error =
            alp_reg::sqrt_for_errors(alpha_I_mult2_error - alpha_I_mult2 * alpha_I_mult2)
                / mult_number_double_lambda.sqrt();
        alpha_J_mult2_error =
            alp_reg::sqrt_for_errors(alpha_J_mult2_error - alpha_J_mult2 * alpha_J_mult2)
                / mult_number_double_lambda.sqrt();
        K_mult2_error = alp_reg::sqrt_for_errors(K_mult2_error - K_mult2 * K_mult2)
            / Tmin(mult_number_double_lambda, mult_number_double_K).sqrt();

        Self::error_in_calculate_main_parameters2m(
            *lambda,
            lambda_error,
            lambda_mult2,
            lambda_mult2_error,
        );
        Self::error_in_calculate_main_parameters2m(*C, C_error, C_mult2, C_mult2_error);
        Self::error_in_calculate_main_parameters2m(*K_C, K_C_error, K_C_mult2, K_C_mult2_error);
        Self::error_in_calculate_main_parameters2m(*a_I, a_I_error, a_I_mult2, a_I_mult2_error);
        Self::error_in_calculate_main_parameters2m(*a_J, a_J_error, a_J_mult2, a_J_mult2_error);
        Self::error_in_calculate_main_parameters2m(
            *sigma,
            sigma_error,
            sigma_mult2,
            sigma_mult2_error,
        );
        Self::error_in_calculate_main_parameters2m(
            *alpha_I,
            alpha_I_error,
            alpha_I_mult2,
            alpha_I_mult2_error,
        );
        Self::error_in_calculate_main_parameters2m(
            *alpha_J,
            alpha_J_error,
            alpha_J_mult2,
            alpha_J_mult2_error,
        );
        Self::error_in_calculate_main_parameters2m(*K, K_error, K_mult2, K_mult2_error);

        *flag_ = true;

        self.m_AI = *a_I;
        self.m_AIError = *a_I_error;
        self.m_AJ = *a_J;
        self.m_AJError = *a_J_error;
        self.m_Sigma = *sigma;
        self.m_SigmaError = *sigma_error;
        self.m_C = *C;
        self.m_CError = *C_error;
        self.m_K = *K;
        self.m_KError = *K_error;
        self.m_Lambda = *lambda;
        self.m_LambdaError = *lambda_error;
        self.m_AlphaI = *alpha_I;
        self.m_AlphaIError = *alpha_I_error;
        self.m_AlphaJ = *alpha_J;
        self.m_AlphaJError = *alpha_J_error;

        self.m_AISbs.resize(self.d_mult_number as usize, 0.0);
        self.m_AJSbs.resize(self.d_mult_number as usize, 0.0);
        self.m_SigmaSbs.resize(self.d_mult_number as usize, 0.0);
        self.m_CSbs.resize(self.d_mult_number as usize, 0.0);
        self.m_KSbs.resize(self.d_mult_number as usize, 0.0);
        self.m_LambdaSbs.resize(self.d_mult_number as usize, 0.0);
        self.m_AlphaISbs.resize(self.d_mult_number as usize, 0.0);
        self.m_AlphaJSbs.resize(self.d_mult_number as usize, 0.0);

        for k in 1..=self.d_mult_number {
            self.m_AISbs[(k - 1) as usize] = a_I_mult[k as usize];
            self.m_AJSbs[(k - 1) as usize] = a_J_mult[k as usize];
            self.m_SigmaSbs[(k - 1) as usize] = sigma_mult[k as usize];
            self.m_CSbs[(k - 1) as usize] = C_mult[k as usize];
            self.m_KSbs[(k - 1) as usize] = K_mult[k as usize];
            self.m_LambdaSbs[(k - 1) as usize] = lambda_mult[k as usize];
            self.m_AlphaISbs[(k - 1) as usize] = alpha_I_mult[k as usize];
            self.m_AlphaJSbs[(k - 1) as usize] = alpha_J_mult[k as usize];
        }

        self.symmetric_parameters_for_symmetric_scheme();
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn memory_release_for_calculate_main_parameters2m(
        _nalp_for_lambda_simulation: i64,
        d_mult_realizations: &mut Vec<i64>,
        d_mult_K_realizations: &mut Vec<i64>,
        lambda_mult: &mut Vec<f64>,
        lambda_mult_error: &mut Vec<f64>,
        C_mult: &mut Vec<f64>,
        C_mult_error: &mut Vec<f64>,
        a_I_mult: &mut Vec<f64>,
        a_I_mult_error: &mut Vec<f64>,
        a_J_mult: &mut Vec<f64>,
        a_J_mult_error: &mut Vec<f64>,
        sigma_mult: &mut Vec<f64>,
        sigma_mult_error: &mut Vec<f64>,
        alpha_I_mult: &mut Vec<f64>,
        alpha_I_mult_error: &mut Vec<f64>,
        alpha_J_mult: &mut Vec<f64>,
        alpha_J_mult_error: &mut Vec<f64>,
        K_C_mult: &mut Vec<f64>,
        K_C_mult_error: &mut Vec<f64>,
        K_mult: &mut Vec<f64>,
        K_mult_error: &mut Vec<f64>,
        Sc_mult: &mut Vec<f64>,
        Sc_mult_error: &mut Vec<f64>,
        alp_distr: &mut Vec<array_positive<f64>>,
        alp_distr_errors: &mut Vec<array_positive<f64>>,
        alp_mult_distr: &mut Vec<Vec<array_positive<f64>>>,
        alp_mult_distr_errors: &mut Vec<Vec<array_positive<f64>>>,
    ) {
        alp_distr.clear();
        alp_distr_errors.clear();
        alp_mult_distr.clear();
        alp_mult_distr_errors.clear();
        d_mult_realizations.clear();
        d_mult_K_realizations.clear();
        lambda_mult.clear();
        lambda_mult_error.clear();
        C_mult.clear();
        C_mult_error.clear();
        a_I_mult.clear();
        a_I_mult_error.clear();
        a_J_mult.clear();
        a_J_mult_error.clear();
        sigma_mult.clear();
        sigma_mult_error.clear();
        alpha_I_mult.clear();
        alpha_I_mult_error.clear();
        alpha_J_mult.clear();
        alpha_J_mult_error.clear();
        K_C_mult.clear();
        K_C_mult_error.clear();
        K_mult.clear();
        K_mult_error.clear();
        Sc_mult.clear();
        Sc_mult_error.clear();
    }

    pub fn output_main_parameters2m_new(
        &mut self,
        nalp_for_lambda_simulation: i64,
        level: i64,
        inside_simulation_flag: &mut bool,
        final_realizations_number_lambda_: i64,
        final_realizations_number_killing_: i64,
    ) -> Result<(), Error> {
        let mut lambda = 0.0;
        let mut lambda_error = 0.0;
        let mut test_difference = 0.0;
        let mut test_difference_error = 0.0;
        let mut C = 0.0;
        let mut C_error = 0.0;
        let mut K_C = 0.0;
        let mut K_C_error = 0.0;
        let mut a_I = 0.0;
        let mut a_I_error = 0.0;
        let mut a_J = 0.0;
        let mut a_J_error = 0.0;
        let mut sigma = 0.0;
        let mut sigma_error = 0.0;
        let mut alpha_I = 0.0;
        let mut alpha_I_error = 0.0;
        let mut alpha_J = 0.0;
        let mut alpha_J_error = 0.0;
        let mut K = 0.0;
        let mut K_error = 0.0;

        let mut flag = false;
        let mut number_of_trials = 0;
        let number_of_trials_threshold = 4;

        while !flag && number_of_trials <= number_of_trials_threshold {
            self.calculate_main_parameters2m(
                final_realizations_number_lambda_,
                final_realizations_number_killing_,
                nalp_for_lambda_simulation,
                level,
                inside_simulation_flag,
                &mut lambda,
                &mut lambda_error,
                &mut test_difference,
                &mut test_difference_error,
                &mut C,
                &mut C_error,
                &mut K_C,
                &mut K_C_error,
                &mut a_I,
                &mut a_I_error,
                &mut a_J,
                &mut a_J_error,
                &mut sigma,
                &mut sigma_error,
                &mut alpha_I,
                &mut alpha_I_error,
                &mut alpha_J,
                &mut alpha_J_error,
                &mut K,
                &mut K_error,
                &mut flag,
            )?;

            number_of_trials += 1;

            if !flag {
                self.randomize_realizations(
                    final_realizations_number_lambda_,
                    final_realizations_number_killing_,
                )?;
            }
        }

        if !flag {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn kill(
        &mut self,
        check_time_: bool,
        ind1_: i64,
        ind2_: i64,
        M_min_: i64,
        lambda_: f64,
        eps_K_: f64,
        K_C_: &mut f64,
        K_C_error_: &mut f64,
        level_: &mut i64,
        diff_opt_: &mut i64,
    ) -> Result<(), Error> {
        let mut flag = false;
        let mut current_level = (M_min_ as f64 * 0.5).floor() as i64;
        let mut recommended_level = 0;
        for i in ind1_..=ind2_ {
            let alp_obj_tmp = self.d_alp_obj[i as usize]
                .as_mut()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            if i - ind1_ + 1 > alp_obj_tmp.d_alp_data.d_minimum_realizations_number {
                alp_obj_tmp.d_check_time_flag = check_time_;
                alp_obj_tmp.d_time_error_flag = check_time_;
            }
        }

        while !flag {
            for i in ind1_..=ind2_ {
                let mut inner_flag = false;
                while !inner_flag {
                    {
                        let alp_obj_tmp = self.d_alp_obj[i as usize]
                            .as_mut()
                            .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                        alp_obj_tmp.d_sentinels_flag = false;
                        alp_obj_tmp.kill_upto_level(M_min_, current_level, None)?;
                    }
                    let success = self.d_alp_obj[i as usize]
                        .as_ref()
                        .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?
                        .d_success;
                    if !success {
                        self.d_alp_obj[i as usize] =
                            Some(Box::new(alp::new(self.d_alp_data.clone())?));
                        {
                            let alp_obj_tmp = self.d_alp_obj[i as usize]
                                .as_mut()
                                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                            if i - ind1_ + 1 > alp_obj_tmp.d_alp_data.d_minimum_realizations_number
                            {
                                alp_obj_tmp.d_check_time_flag = check_time_;
                                alp_obj_tmp.d_time_error_flag = check_time_;
                            }
                            let mut simulation_success = false;
                            while !simulation_success {
                                alp_obj_tmp.simulate_alp_upto_the_given_level(M_min_)?;
                                simulation_success = alp_obj_tmp.d_success;
                            }
                        }
                    }
                    inner_flag = self.d_alp_obj[i as usize]
                        .as_ref()
                        .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?
                        .d_success;
                }
            }

            flag = self.check_K_criterion_during_killing(
                ind1_,
                ind2_,
                lambda_,
                eps_K_,
                current_level,
                &mut recommended_level,
                diff_opt_,
                K_C_,
                K_C_error_,
            )?;
            current_level = recommended_level;
        }
        *level_ = current_level;
        Ok(())
    }

    pub fn quick_test(&mut self, trials_number_: i64, max_time_: f64) -> Result<(), Error> {
        if trials_number_ <= 0 {
            return Err(Error::new(
                "Unexpected error in alp_sim::quick_test\n".to_string(),
                1,
            ));
        }

        let mut check_time_flag = false;
        if max_time_ > 0.0 {
            check_time_flag = true;
        }

        let alp_number = 5;
        let p_thres: f64 = 1e-10;

        let lambda_ungapped = self
            .d_alp_data
            .d_is
            .as_ref()
            .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?
            .d_ungap_lambda;
        if lambda_ungapped <= 0.0 {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }
        let score_diff = round(&(-p_thres.ln() / lambda_ungapped)) as i64;

        let max_number_of_unsuccessful_objects2 = (0.5
            * trials_number_ as f64
            * (self.d_alp_data.d_eps_K + self.d_alp_data.d_eps_lambda))
            .floor() as i64;
        let mut number_of_unsuccessful_objects2 = 0;

        let max_time_store = self.d_alp_data.d_max_time;
        if check_time_flag {
            self.d_alp_data.d_max_time = max_time_;
        }

        for _ in 0..trials_number_ {
            let mut alp_obj_tmp: Option<Box<alp>> = None;
            let mut success3 = false;
            while !success3 {
                alp_obj_tmp = Some(Box::new(alp::new(self.d_alp_data.clone())?));
                self.d_alp_data.d_memory_size_in_MB += std::mem::size_of::<alp>() as f64 / mb_bytes;
                {
                    let alp_ref = alp_obj_tmp.as_mut().unwrap();
                    alp_ref.d_check_time_flag = check_time_flag;
                    alp_ref.d_time_error_flag = check_time_flag;
                    alp_ref.simulate_alp_upto_the_given_number(alp_number + 1)?;
                    success3 = alp_ref.d_success;
                }

                if !success3 {
                    alp_obj_tmp = None;
                    self.d_alp_data.d_memory_size_in_MB -=
                        std::mem::size_of::<alp>() as f64 / mb_bytes;
                    number_of_unsuccessful_objects2 += 1;
                    if number_of_unsuccessful_objects2 > max_number_of_unsuccessful_objects2 {
                        if check_time_flag {
                            self.d_alp_data.d_max_time = max_time_store;
                        }
                        return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
                    }
                }
            }

            let alp_ref = alp_obj_tmp.as_mut().unwrap();
            let last_alp = alp_ref.d_alp.d_elem[alp_number as usize];
            let mut M_upper_level = last_alp + score_diff;
            alp_ref.d_sentinels_flag = false;
            alp_ref.kill_upto_level(last_alp, last_alp - score_diff, Some(&mut M_upper_level))?;
            if !alp_ref.d_success {
                number_of_unsuccessful_objects2 += 1;
                if number_of_unsuccessful_objects2 > max_number_of_unsuccessful_objects2 {
                    if check_time_flag {
                        self.d_alp_data.d_max_time = max_time_store;
                    }
                    return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
                }
            }

            drop(alp_obj_tmp);
            self.d_alp_data.d_memory_size_in_MB -= std::mem::size_of::<alp>() as f64 / mb_bytes;
        }

        if check_time_flag {
            self.d_alp_data.d_max_time = max_time_store;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn get_minimal_simulation(
        &mut self,
        ind1_: i64,
        ind2_: i64,
        M_min_: &mut i64,
        nalp_: &mut i64,
        nalp_lambda_: &mut i64,
        C_calculation_: bool,
        check_time_flag_: bool,
    ) -> Result<(), Error> {
        let mut alp_distr: Vec<array_positive<f64>> = Vec::new();
        let mut alp_distr_errors: Vec<array_positive<f64>> = Vec::new();

        let add_alp_number = 3;
        let mut add_alp_number_count = 0;
        let max_alp_number = 30;

        if self.d_n_alp_obj < ind1_ || self.d_n_alp_obj - 1 > ind2_ {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }

        *nalp_ = 0;

        if self.d_alp_obj.len() < ind2_ as usize + 1 {
            self.d_alp_obj.resize_with(ind2_ as usize + 1, || None);
        }
        for i in self.d_n_alp_obj..=ind2_ {
            self.d_alp_obj[i as usize] = None;
            self.d_alp_obj[i as usize] = Some(Box::new(alp::new(self.d_alp_data.clone())?));
            self.d_alp_data.d_memory_size_in_MB += std::mem::size_of::<alp>() as f64 / mb_bytes;
            let alp_obj_tmp = self.d_alp_obj[i as usize]
                .as_mut()
                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
            alp_obj_tmp.d_check_time_flag = check_time_flag_;
            alp_obj_tmp.d_time_error_flag = check_time_flag_;
        }

        self.d_n_alp_obj = ind2_ + 1;

        let mut M_min_flag = false;
        let mut nalp_flag = false;
        let mut criterion_flag = false;
        let mut number_of_fails = 0;
        let number_of_fails_threshold = 5;

        while !criterion_flag {
            if *nalp_ >= max_alp_number {
                Self::memory_release_for_get_minimal_simulation(
                    *nalp_,
                    &mut alp_distr,
                    &mut alp_distr_errors,
                );
                return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 1));
            }

            for i in ind1_..=ind2_ {
                {
                    let alp_obj_tmp = self.d_alp_obj[i as usize]
                        .as_mut()
                        .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                    alp_obj_tmp.d_check_time_flag = check_time_flag_;
                    alp_obj_tmp.d_time_error_flag = check_time_flag_;

                    if alp_obj_tmp.d_nalp < *nalp_ + 1 {
                        alp_obj_tmp.simulate_alp_upto_the_given_number(*nalp_ + 1)?;
                    }
                }

                let success = self.d_alp_obj[i as usize]
                    .as_ref()
                    .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?
                    .d_success;
                if !success {
                    self.d_alp_obj[i as usize] = None;

                    let mut success2 = false;
                    while !success2 {
                        self.d_alp_obj[i as usize] =
                            Some(Box::new(alp::new(self.d_alp_data.clone())?));
                        {
                            let alp_obj_tmp = self.d_alp_obj[i as usize]
                                .as_mut()
                                .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                            for j in 0..=*nalp_ {
                                alp_obj_tmp.simulate_alp_upto_the_given_number(j + 1)?;
                            }
                            success2 = alp_obj_tmp.d_success;
                        }
                        if !success2 {
                            self.d_alp_obj[i as usize] = None;
                        }
                    }
                }
            }

            *nalp_ += 1;

            let mut inside_simulation_flag = false;
            let mut lambda = 0.0;
            criterion_flag = self.the_criterion(
                *nalp_,
                nalp_lambda_,
                0,
                ind2_,
                &mut alp_distr,
                &mut alp_distr_errors,
                M_min_,
                &mut M_min_flag,
                &mut nalp_flag,
                &mut inside_simulation_flag,
                C_calculation_,
                Some(&mut lambda),
                None,
            )?;

            if inside_simulation_flag {
                if lambda <= 0.0 {
                    criterion_flag = false;
                    inside_simulation_flag = false;
                }
            } else {
                criterion_flag = false;
            }

            if !inside_simulation_flag {
                number_of_fails += 1;

                Self::memory_release_for_get_minimal_simulation(
                    *nalp_,
                    &mut alp_distr,
                    &mut alp_distr_errors,
                );

                M_min_flag = false;
                nalp_flag = false;
                *nalp_ = 0;
                criterion_flag = false;

                for i in ind1_..=ind2_ {
                    self.d_alp_obj[i as usize] = None;
                }

                if number_of_fails > number_of_fails_threshold {
                    return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
                }

                for i in ind1_..=ind2_ {
                    self.d_alp_obj[i as usize] = Some(Box::new(alp::new(self.d_alp_data.clone())?));
                    let alp_obj_tmp = self.d_alp_obj[i as usize]
                        .as_mut()
                        .ok_or_else(|| Error::new("Unexpected error\n".to_string(), 4))?;
                    alp_obj_tmp.d_check_time_flag = check_time_flag_;
                    alp_obj_tmp.d_time_error_flag = check_time_flag_;
                }

                continue;
            }

            if criterion_flag {
                add_alp_number_count += 1;
                if add_alp_number_count < add_alp_number {
                    criterion_flag = false;
                }

                if criterion_flag {
                    criterion_flag = self.check_K_criterion(
                        *nalp_,
                        ind1_,
                        ind2_,
                        lambda,
                        self.d_alp_data.d_eps_K,
                        M_min_,
                    )?;
                }
            } else {
                add_alp_number_count = 0;
            }
        }

        *nalp_lambda_ = *nalp_;

        Self::memory_release_for_get_minimal_simulation(
            *nalp_,
            &mut alp_distr,
            &mut alp_distr_errors,
        );
        Ok(())
    }

    pub fn memory_release_for_get_minimal_simulation(
        _nalp_: i64,
        alp_distr: &mut Vec<array_positive<f64>>,
        alp_distr_errors: &mut Vec<array_positive<f64>>,
    ) {
        alp_distr.clear();
        alp_distr_errors.clear();
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sigma_calculation(
        delta_I_aver_: f64,
        delta_I_aver_error_: f64,
        delta_J_aver_: f64,
        delta_J_aver_error_: f64,
        delta_E_aver_: f64,
        delta_E_aver_error_: f64,
        cov_E_E_aver_: f64,
        cov_E_E_aver_error_: f64,
        cov_I_J_aver_: f64,
        cov_I_J_aver_error_: f64,
        sigma_: &mut f64,
        sigma_error_: &mut f64,
    ) {
        let nom1_1 = delta_I_aver_ * delta_J_aver_;
        let nom2_2 = delta_E_aver_ * delta_E_aver_;
        let den = nom2_2 * delta_E_aver_;
        let nom1 = nom1_1 * cov_E_E_aver_;
        let nom2 = nom2_2 * cov_I_J_aver_;
        *sigma_ = (nom1 + nom2) / den;

        let nom1_sigma_error = alp_data::error_of_the_product(
            delta_I_aver_,
            delta_I_aver_error_,
            delta_J_aver_,
            delta_J_aver_error_,
        );
        let nom1_sigma_error = alp_data::error_of_the_product(
            nom1_1,
            nom1_sigma_error,
            cov_E_E_aver_,
            cov_E_E_aver_error_,
        );
        let nom2_sigma_error_2 = alp_data::error_of_the_product(
            delta_E_aver_,
            delta_E_aver_error_,
            delta_E_aver_,
            delta_E_aver_error_,
        );
        let nom2_sigma_error = alp_data::error_of_the_product(
            nom2_2,
            nom2_sigma_error_2,
            cov_I_J_aver_,
            cov_I_J_aver_error_,
        );
        let den_sigma_error = alp_data::error_of_the_product(
            nom2_2,
            nom2_sigma_error_2,
            delta_E_aver_,
            delta_E_aver_error_,
        );
        let nom_sigma_error = alp_data::error_of_the_sum(nom1_sigma_error, nom2_sigma_error);
        *sigma_error_ =
            alp_data::error_of_the_ratio(nom1 + nom2, nom_sigma_error, den, den_sigma_error);
    }

    pub fn get_root(res_tmp_: &[f64], point_: f64) -> Result<f64, Error> {
        if res_tmp_.is_empty() {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }
        let mut p = 0usize;
        let mut d1 = (point_ - res_tmp_[0]).abs();
        for i in 1..res_tmp_.len() {
            let d2 = (point_ - res_tmp_[i]).abs();
            if d2 < d1 {
                p = i;
                d1 = d2;
            }
        }
        Ok(res_tmp_[p])
    }

    pub fn function_for_lambda_calculation(
        lambda_: f64,
        data_: &mut struct_for_lambda_calculation,
    ) -> Result<f64, Error> {
        let nalp = data_.d_nalp;
        if nalp < 1 {
            return Err(Error::new("Unexpected error\n".to_string(), 4));
        }
        let mut expect = vec![0.0; nalp as usize];
        let mut expect_errors = vec![0.0; nalp as usize];
        for k in 1..=nalp {
            let tmp = &data_.d_alp_distr[k as usize];
            let tmp_errors = &data_.d_alp_distr_errors[k as usize];
            let mut val = 0.0;
            let mut val_error = 0.0;
            for j in 0..=tmp.d_dim {
                if tmp.d_elem[j as usize] <= 0.0 {
                    continue;
                }
                let exp_tmp = (lambda_ * j as f64).exp();
                val += exp_tmp * tmp.d_elem[j as usize];
                val_error += exp_tmp * exp_tmp * tmp_errors.d_elem[j as usize];
            }
            val_error = alp_reg::sqrt_for_errors(val_error);
            expect[(k - 1) as usize] = val;
            expect_errors[(k - 1) as usize] = val_error;
        }
        data_.d_last_sum = expect[(nalp - 1) as usize];
        data_.d_last_sum_error = expect_errors[(nalp - 1) as usize];

        if data_.d_calculate_alp_number {
            let mut tmp = 0.0;
            for err in &expect_errors {
                if *err != 0.0 {
                    tmp += 1.0 / (err * err);
                }
            }
            let mut tmp_alp = nalp;
            let mut tmp1 = 0.0;
            for k in (0..nalp).rev() {
                if expect_errors[k as usize] != 0.0 {
                    tmp1 += 1.0 / (expect_errors[k as usize] * expect_errors[k as usize]);
                }
                if tmp1 > 0.2 * tmp {
                    tmp_alp = k + 1;
                    break;
                }
            }
            data_.d_alp_number = tmp_alp;
        }

        if nalp == 1 {
            data_.d_f_error = expect_errors[0];
            return Ok(expect[0] - 1.0);
        }

        let min_length = 0;
        let number_of_elements = nalp;
        let cut_left_tail = true;
        let cut_right_tail = false;
        let y = 2.0;
        let mut beta0 = 0.0;
        let mut beta1 = 0.0;
        let mut beta0_error = 0.0;
        let mut beta1_error = 0.0;
        let mut k1_opt = 0;
        let mut k2_opt = 0;
        let mut res_was_calculated = false;
        alp_reg::robust_regression_sum_with_cut_LSM(
            min_length,
            number_of_elements,
            &expect,
            &mut expect_errors,
            cut_left_tail,
            cut_right_tail,
            y,
            &mut beta0,
            &mut beta1,
            &mut beta0_error,
            &mut beta1_error,
            &mut k1_opt,
            &mut k2_opt,
            &mut res_was_calculated,
        )?;
        if !res_was_calculated {
            return Err(Error::new("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(), 3));
        }
        data_.d_f_error = beta1_error;
        Ok(beta1)
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
    fn test_round_and_relative_error_match_cpp_cases() {
        assert_eq!(alp_sim::round_double(12.34, 1), 12.3);
        assert_eq!(alp_sim::round_double(12.36, 1), 12.4);
        assert_eq!(alp_sim::relative_error_in_percents(0.0, 1.0), f64::MAX);
        assert_eq!(alp_sim::relative_error_in_percents(20.0, 1.0), 5.0);
    }

    #[test]
    fn test_new_with_simulation_requires_randomization_state() {
        let mut data = test_data();
        data.d_rand_all = None;
        assert!(alp_sim::new_with_simulation(data).is_err());
    }

    #[test]
    fn test_get_number_of_subsimulations_bounds() {
        assert!(alp_sim::get_number_of_subsimulations(5).is_err());
        assert_eq!(alp_sim::get_number_of_subsimulations(6).unwrap(), 3);
        assert_eq!(alp_sim::get_number_of_subsimulations(10_000).unwrap(), 20);
    }

    #[test]
    fn test_get_root_nearest_and_empty_error() {
        assert_eq!(alp_sim::get_root(&[0.1, 0.5, 0.9], 0.6).unwrap(), 0.5);
        assert!(alp_sim::get_root(&[], 0.6).is_err());
    }

    #[test]
    fn test_lambda_exp_rejects_sentinel() {
        assert_eq!(alp_sim::lambda_exp(1, &[0.0, 2.0]).unwrap(), 2.0);
        assert!(alp_sim::lambda_exp(1, &[0.0, -1.0]).is_err());
    }

    #[test]
    fn test_error_in_calculate_main_parameters2m() {
        let mut err = 0.0;
        alp_sim::error_in_calculate_main_parameters2m(10.0, &mut err, 5.0, 2.0);
        assert_eq!(err, 4.0);
        alp_sim::error_in_calculate_main_parameters2m(0.0, &mut err, 5.0, 2.0);
        assert_eq!(err, 2.0);
    }

    #[test]
    fn test_function_for_lambda_calculation_single_alp() {
        let mut distr = array_positive::new(None).unwrap();
        distr.set_elem(0, 0.25);
        distr.set_elem(1, 0.75);
        let mut errors = array_positive::new(None).unwrap();
        errors.set_elem(0, 0.01);
        errors.set_elem(1, 0.04);
        let mut data = struct_for_lambda_calculation {
            d_alp_distr: vec![array_positive::new(None).unwrap(), distr],
            d_alp_distr_errors: vec![array_positive::new(None).unwrap(), errors],
            d_nalp: 1,
            ..Default::default()
        };
        let f = alp_sim::function_for_lambda_calculation(0.0, &mut data).unwrap();
        assert_eq!(f, 0.0);
        assert!((data.d_f_error - alp_reg::sqrt_for_errors(0.05)).abs() < 1e-15);
    }

    #[test]
    fn test_get_single_realization_allocates_and_sets_flags() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut obj = None;
        let mut success = false;
        let mut eps = 0.0;
        sim.get_single_realization(false, 0, 0, false, 0, 3, &mut obj, &mut success, &mut eps)
            .unwrap();
        let obj = obj.unwrap();
        assert!(success);
        assert!(obj.d_single_realiztion_calculation_flag);
        assert!(!obj.d_check_time_flag);
        assert_eq!(obj.d_diff_opt, 3);
        assert_eq!(eps, 0.05);
    }

    #[test]
    fn test_generate_random_permulation_contains_each_index() {
        let sim = alp_sim::new(test_data()).unwrap();
        let mut perm = sim.generate_random_permulation(8).unwrap();
        perm.sort();
        assert_eq!(perm, vec![0, 1, 2, 3, 4, 5, 6, 7]);
    }

    #[test]
    fn test_randomize_realizations_ind_bounds_and_noop() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        sim.d_n_alp_obj = 1;
        sim.d_alp_obj
            .push(Some(Box::new(alp::new(test_data()).unwrap())));
        sim.randomize_realizations_ind(0, 0).unwrap();
        assert!(sim.randomize_realizations_ind(0, 1).is_err());
    }

    #[test]
    fn test_symmetric_parameters_for_symmetric_scheme_averages_pairs() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        sim.m_AI = 2.0;
        sim.m_AJ = 4.0;
        sim.m_AIError = 0.2;
        sim.m_AJError = 0.4;
        sim.m_AlphaI = 6.0;
        sim.m_AlphaJ = 8.0;
        sim.m_AlphaIError = 0.6;
        sim.m_AlphaJError = 0.8;
        sim.symmetric_parameters_for_symmetric_scheme();
        assert_eq!(sim.m_AI, 3.0);
        assert_eq!(sim.m_AJ, 3.0);
        assert_eq!(sim.m_AIError, 0.30000000000000004);
        assert_eq!(sim.m_AJError, 0.30000000000000004);
        assert_eq!(sim.m_AlphaI, 7.0);
        assert_eq!(sim.m_AlphaJ, 7.0);
    }

    #[test]
    fn test_symmetric_parameters_for_asymmetric_scheme_keeps_values() {
        let mut data = test_data();
        data.d_RR2 = vec![0.25, 0.75];
        let mut sim = alp_sim::new(data).unwrap();
        sim.m_AI = 2.0;
        sim.m_AJ = 4.0;
        sim.symmetric_parameters_for_symmetric_scheme();
        assert_eq!(sim.m_AI, 2.0);
        assert_eq!(sim.m_AJ, 4.0);
    }

    #[test]
    fn test_get_and_allocate_alp_distribution_accumulates_weighted_scores() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut a1 = alp::new(test_data()).unwrap();
        let mut a2 = alp::new(test_data()).unwrap();
        a1.d_alp.set_elem(1, 2);
        a1.d_alp_weights.set_elem(1, 0.5);
        a2.d_alp.set_elem(1, 4);
        a2.d_alp_weights.set_elem(1, 1.5);
        sim.d_alp_obj.push(Some(Box::new(a1)));
        sim.d_alp_obj.push(Some(Box::new(a2)));
        sim.d_n_alp_obj = 2;

        let mut distr = Vec::new();
        let mut errors = Vec::new();
        sim.get_and_allocate_alp_distribution(0, 1, &mut distr, &mut errors, 1)
            .unwrap();
        assert!((distr[1].d_elem[2] - 0.25).abs() < 1e-15);
        assert!((distr[1].d_elem[4] - 0.75).abs() < 1e-15);
        assert!((errors[1].d_elem[2] - 0.03125).abs() < 1e-15);
        assert!((errors[1].d_elem[4] - 0.28125).abs() < 1e-15);
    }

    #[test]
    fn test_check_k_criterion_known_counts() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut a = alp::new(test_data()).unwrap();
        a.d_alp.set_elem(1, 4);
        a.d_alp_weights.set_elem(1, 2.0);
        a.d_cells_counts.set_elem(3, 2);
        a.d_cells_counts.set_elem(4, 1);
        sim.d_alp_obj.push(Some(Box::new(a)));
        sim.d_n_alp_obj = 1;
        let mut m_min = 0;
        let ok = sim
            .check_K_criterion(1, 0, 0, 1.0, 0.5, &mut m_min)
            .unwrap();
        assert!(ok);
        assert_eq!(m_min, 4);
    }

    #[test]
    fn test_check_k_criterion_during_killing_known_counts() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut a = alp::new(test_data()).unwrap();
        a.d_M = 4;
        a.d_nalp_killing = 0;
        a.d_alp_weights.set_elem(0, 2.0);
        a.d_cells_counts.set_elem(3, 2);
        a.d_cells_counts.set_elem(4, 1);
        sim.d_alp_obj.push(Some(Box::new(a)));
        sim.d_n_alp_obj = 1;
        let mut recommended = 0;
        let mut diff_opt = 0;
        let mut kc = 0.0;
        let mut kc_error = 0.0;
        let ok = sim
            .check_K_criterion_during_killing(
                0,
                0,
                1.0,
                0.5,
                1,
                &mut recommended,
                &mut diff_opt,
                &mut kc,
                &mut kc_error,
            )
            .unwrap();
        assert!(ok);
        assert_eq!(recommended, 1);
        assert_eq!(diff_opt, 3);
        assert!(kc > 0.0);
        assert_eq!(kc_error, 0.0);
    }

    #[test]
    fn test_the_criterion_sets_lambda_and_flags() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut a = alp::new(test_data()).unwrap();
        a.d_alp.set_elem(1, 0);
        a.d_alp_weights.set_elem(1, 1.0);
        sim.d_alp_obj.push(Some(Box::new(a)));
        sim.d_n_alp_obj = 1;

        let mut alp_distr = Vec::new();
        let mut alp_distr_errors = Vec::new();
        let mut m_min = -1;
        let mut m_min_flag = true;
        let mut nalp_flag = false;
        let mut inside = false;
        let mut nalp_lambda = 0;
        let mut lambda = -1.0;
        let mut lambda_error = -1.0;
        let ok = sim
            .the_criterion(
                1,
                &mut nalp_lambda,
                0,
                0,
                &mut alp_distr,
                &mut alp_distr_errors,
                &mut m_min,
                &mut m_min_flag,
                &mut nalp_flag,
                &mut inside,
                false,
                Some(&mut lambda),
                Some(&mut lambda_error),
            )
            .unwrap();
        assert!(ok);
        assert!(inside);
        assert!(nalp_flag);
        assert!(!m_min_flag);
        assert_eq!(m_min, 0);
        assert_eq!(nalp_lambda, 1);
        assert_eq!(sim.d_lambda_tmp.d_elem[1], lambda);
        assert_eq!(sim.d_lambda_tmp_errors.d_elem[1], lambda_error);
    }

    #[test]
    fn test_quick_test_rejects_nonpositive_trials() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        assert!(sim.quick_test(0, 0.0).is_err());
    }

    #[test]
    fn test_get_minimal_simulation_rejects_inconsistent_object_range() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        sim.d_n_alp_obj = 2;
        let mut m_min = 0;
        let mut nalp = 0;
        let mut nalp_lambda = 0;
        assert!(sim
            .get_minimal_simulation(0, 0, &mut m_min, &mut nalp, &mut nalp_lambda, false, false,)
            .is_err());
    }

    #[test]
    fn test_memory_release_for_get_minimal_simulation_clears_vectors() {
        let mut distr = vec![array_positive::new(None).unwrap()];
        let mut errors = vec![array_positive::new(None).unwrap()];
        alp_sim::memory_release_for_get_minimal_simulation(1, &mut distr, &mut errors);
        assert!(distr.is_empty());
        assert!(errors.is_empty());
    }

    #[test]
    fn test_calculate_main_parameters2m_rejects_killing_more_than_lambda() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut inside = false;
        let mut lambda = 0.0;
        let mut lambda_error = 0.0;
        let mut test_difference = 0.0;
        let mut test_difference_error = 0.0;
        let mut c = 0.0;
        let mut c_error = 0.0;
        let mut k_c = 0.0;
        let mut k_c_error = 0.0;
        let mut a_i = 0.0;
        let mut a_i_error = 0.0;
        let mut a_j = 0.0;
        let mut a_j_error = 0.0;
        let mut sigma = 0.0;
        let mut sigma_error = 0.0;
        let mut alpha_i = 0.0;
        let mut alpha_i_error = 0.0;
        let mut alpha_j = 0.0;
        let mut alpha_j_error = 0.0;
        let mut k = 0.0;
        let mut k_error = 0.0;
        let mut flag = true;
        assert!(sim
            .calculate_main_parameters2m(
                6,
                7,
                1,
                0,
                &mut inside,
                &mut lambda,
                &mut lambda_error,
                &mut test_difference,
                &mut test_difference_error,
                &mut c,
                &mut c_error,
                &mut k_c,
                &mut k_c_error,
                &mut a_i,
                &mut a_i_error,
                &mut a_j,
                &mut a_j_error,
                &mut sigma,
                &mut sigma_error,
                &mut alpha_i,
                &mut alpha_i_error,
                &mut alpha_j,
                &mut alpha_j_error,
                &mut k,
                &mut k_error,
                &mut flag,
            )
            .is_err());
        assert!(!flag);
    }

    #[test]
    fn test_memory_release_for_calculate_main_parameters2m_clears_vectors() {
        let mut ints = vec![1];
        let mut ints2 = vec![2];
        let mut v1 = vec![1.0];
        let mut v2 = vec![2.0];
        let mut v3 = vec![3.0];
        let mut v4 = vec![4.0];
        let mut v5 = vec![5.0];
        let mut v6 = vec![6.0];
        let mut v7 = vec![7.0];
        let mut v8 = vec![8.0];
        let mut v9 = vec![9.0];
        let mut v10 = vec![10.0];
        let mut v11 = vec![11.0];
        let mut v12 = vec![12.0];
        let mut v13 = vec![13.0];
        let mut v14 = vec![14.0];
        let mut v15 = vec![15.0];
        let mut v16 = vec![16.0];
        let mut v17 = vec![17.0];
        let mut v18 = vec![18.0];
        let mut v19 = vec![19.0];
        let mut v20 = vec![20.0];
        let mut alp_distr = vec![array_positive::new(None).unwrap()];
        let mut alp_distr_errors = vec![array_positive::new(None).unwrap()];
        let mut alp_mult_distr = vec![vec![array_positive::new(None).unwrap()]];
        let mut alp_mult_distr_errors = vec![vec![array_positive::new(None).unwrap()]];
        alp_sim::memory_release_for_calculate_main_parameters2m(
            1,
            &mut ints,
            &mut ints2,
            &mut v1,
            &mut v2,
            &mut v3,
            &mut v4,
            &mut v5,
            &mut v6,
            &mut v7,
            &mut v8,
            &mut v9,
            &mut v10,
            &mut v11,
            &mut v12,
            &mut v13,
            &mut v14,
            &mut v15,
            &mut v16,
            &mut v17,
            &mut v18,
            &mut v19,
            &mut v20,
            &mut alp_distr,
            &mut alp_distr_errors,
            &mut alp_mult_distr,
            &mut alp_mult_distr_errors,
        );
        assert!(ints.is_empty());
        assert!(ints2.is_empty());
        assert!(v1.is_empty());
        assert!(v20.is_empty());
        assert!(alp_distr.is_empty());
        assert!(alp_distr_errors.is_empty());
        assert!(alp_mult_distr.is_empty());
        assert!(alp_mult_distr_errors.is_empty());
    }

    #[test]
    fn test_output_main_parameters2m_new_propagates_main_parameter_error() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut inside = false;
        assert!(sim
            .output_main_parameters2m_new(1, 0, &mut inside, 6, 7)
            .is_err());
    }

    #[test]
    fn test_calculate_lambda_single_alp_distribution() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut distr = array_positive::new(None).unwrap();
        distr.set_elem(0, 0.25);
        distr.set_elem(1, 0.75);
        let mut errors = array_positive::new(None).unwrap();
        errors.set_elem(0, 0.01);
        errors.set_elem(1, 0.04);
        let alp_distr = vec![array_positive::new(None).unwrap(), distr];
        let alp_distr_errors = vec![array_positive::new(None).unwrap(), errors];
        let mut nalp_thr = 0;
        let mut inside = false;
        let mut lambda = -1.0;
        let mut lambda_error = -1.0;
        let mut test_difference = -1.0;
        let mut test_difference_error = -1.0;
        sim.calculate_lambda(
            false,
            1,
            &mut nalp_thr,
            &mut inside,
            &alp_distr,
            &alp_distr_errors,
            &mut lambda,
            &mut lambda_error,
            &mut test_difference,
            &mut test_difference_error,
        )
        .unwrap();
        assert!(inside);
        assert_eq!(lambda, 0.0);
        assert_eq!(nalp_thr, 1);
        assert_eq!(lambda_error, 0.0);
        assert_eq!(test_difference, -1.0);
        assert_eq!(test_difference_error, -1.0);
    }

    #[test]
    fn test_calculate_c_single_alp_distribution() {
        let mut distr = array_positive::new(None).unwrap();
        distr.set_elem(0, 0.25);
        distr.set_elem(1, 0.75);
        let mut errors = array_positive::new(None).unwrap();
        errors.set_elem(0, 0.01);
        errors.set_elem(1, 0.04);
        let alp_distr = vec![array_positive::new(None).unwrap(), distr];
        let alp_distr_errors = vec![array_positive::new(None).unwrap(), errors];
        let mut c = -1.0;
        let mut c_error = -1.0;
        let mut sc = -1.0;
        let mut sc_error = -1.0;
        alp_sim::calculate_C(
            0,
            1,
            &alp_distr,
            &alp_distr_errors,
            1.0,
            0.01,
            &mut c,
            &mut c_error,
            &mut sc,
            &mut sc_error,
        )
        .unwrap();
        assert!(c >= 0.0);
        assert!(c_error >= 0.0);
        assert!(sc > 1.0);
        assert!(sc_error > 0.0);
    }

    #[test]
    fn test_calculate_fsc_single_alp_path() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut a = alp::new(test_data()).unwrap();
        a.d_alp.set_elem(0, 0);
        a.d_alp.set_elem(1, 1);
        a.d_alp_weights.set_elem(1, 1.0);
        a.d_H_I.set_elem(0, 0);
        a.d_H_I.set_elem(1, 2);
        a.d_H_J.set_elem(0, 0);
        a.d_H_J.set_elem(1, 3);
        sim.d_alp_obj.push(Some(Box::new(a)));
        sim.d_n_alp_obj = 1;
        let mut distr = array_positive::new(None).unwrap();
        distr.set_elem(1, 1.0);
        let alp_distr = vec![array_positive::new(None).unwrap(), distr];
        let mut a_i = -1.0;
        let mut a_i_error = -1.0;
        let mut a_j = -1.0;
        let mut a_j_error = -1.0;
        let mut sigma = -1.0;
        let mut sigma_error = -1.0;
        let mut alpha_i = -1.0;
        let mut alpha_i_error = -1.0;
        let mut alpha_j = -1.0;
        let mut alpha_j_error = -1.0;
        sim.calculate_FSC(
            1,
            0,
            0,
            &alp_distr,
            0.0,
            1.0,
            &mut a_i,
            &mut a_i_error,
            &mut a_j,
            &mut a_j_error,
            &mut sigma,
            &mut sigma_error,
            &mut alpha_i,
            &mut alpha_i_error,
            &mut alpha_j,
            &mut alpha_j_error,
        )
        .unwrap();
        assert_eq!(a_i, 2.0);
        assert_eq!(a_j, 3.0);
        assert_eq!(sigma, 0.0);
        assert_eq!(alpha_i, 0.0);
        assert_eq!(alpha_j, 0.0);
        assert!(a_i_error >= 0.0);
        assert!(a_j_error >= 0.0);
        assert!(sigma_error >= 0.0);
    }

    #[test]
    fn test_kill_uses_existing_successful_realization() {
        let mut sim = alp_sim::new(test_data()).unwrap();
        let mut a = alp::new(test_data()).unwrap();
        a.d_is_now = false;
        a.d_M = 4;
        a.d_H_matr_len = 0;
        a.d_H_edge_max[0] = 0;
        a.d_nalp_killing = 0;
        a.d_alp_weights.set_elem(0, 2.0);
        a.d_cells_counts.set_elem(3, 2);
        a.d_cells_counts.set_elem(4, 1);
        sim.d_alp_obj.push(Some(Box::new(a)));
        sim.d_n_alp_obj = 1;
        let mut kc = 0.0;
        let mut kc_error = 0.0;
        let mut level = 0;
        let mut diff_opt = 0;
        sim.kill(
            false,
            0,
            0,
            4,
            1.0,
            0.5,
            &mut kc,
            &mut kc_error,
            &mut level,
            &mut diff_opt,
        )
        .unwrap();
        assert_eq!(level, 2);
        assert_eq!(diff_opt, 2);
        assert!(kc > 0.0);
        assert_eq!(kc_error, 0.0);
    }
}
