#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use crate::stats::alp_localmaxstat_matrix::LocalMaxStatMatrix;
use crate::stats::pvalues::{pvalues, ALP_set_of_parameters};
use crate::stats::sls_alp_data::{alp_data, struct_for_randomization};
use crate::stats::sls_alp_sim::alp_sim;
use crate::stats::sls_basic::{
    get_current_time, one_minus_exp_function, AlignmentEvaluerParameters,
    AlignmentEvaluerParametersWithErrors, Error, QUICK_TESTS_TRIALS_NUMBER,
};

pub const default_importance_sampling_temperature: f64 = 1.07;

unsafe extern "C" {
    fn srand(seed: u32);
}

#[derive(Debug, Clone, PartialEq)]
pub struct gapped_computation_parameters_struct {
    pub d_first_stage_preliminary_realizations_numbers_ALP: Vec<i64>,
    pub d_preliminary_realizations_numbers_ALP: Vec<i64>,
    pub d_preliminary_realizations_numbers_killing: Vec<i64>,
    pub d_total_realizations_number_with_ALP: i64,
    pub d_total_realizations_number_with_killing: i64,
    pub d_parameters_flag: bool,
    pub d_max_time_for_quick_tests: f64,
    pub d_max_time_with_computation_parameters: f64,
}

impl Default for gapped_computation_parameters_struct {
    fn default() -> Self {
        Self {
            d_first_stage_preliminary_realizations_numbers_ALP: Vec::new(),
            d_preliminary_realizations_numbers_ALP: Vec::new(),
            d_preliminary_realizations_numbers_killing: Vec::new(),
            d_total_realizations_number_with_ALP: 0,
            d_total_realizations_number_with_killing: 0,
            d_parameters_flag: false,
            d_max_time_for_quick_tests: -1.0,
            d_max_time_with_computation_parameters: -1.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct AlignmentEvaluer {
    d_params: ALP_set_of_parameters,
    d_gapped_computation_parameters: gapped_computation_parameters_struct,
}

impl AlignmentEvaluer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn initGapless(
        &mut self,
        alphabetSize_: i64,
        substitutionScoreMatrix_: &[Vec<i64>],
        letterFreqs1_: &[f64],
        letterFreqs2_: &[f64],
        mut max_time_: f64,
    ) -> Result<(), Error> {
        let mut CurrentTime1 = 0.0;
        get_current_time(&mut CurrentTime1);

        let function_name = "void AlignmentEvaluer::initGapless";
        let (letterFreqs1_normalized, letterFreqs2_normalized) = self
            .assert_Gapless_input_parameters(
                alphabetSize_,
                letterFreqs1_,
                letterFreqs2_,
                function_name,
            )?;

        if max_time_ <= 0.0 {
            max_time_ = 60.0;
        }

        self.d_params.d_params_flag = false;

        let local_max_stat_matrix = LocalMaxStatMatrix::new(
            alphabetSize_ as usize,
            Some(substitutionScoreMatrix_),
            Some(&letterFreqs1_normalized),
            Some(&letterFreqs2_normalized),
            alphabetSize_ as usize,
            max_time_,
        );

        if local_max_stat_matrix.getTerminated() {
            return Err(Error::new(
                "Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(),
                3,
            ));
        }

        let calculation_error = 1e-6;

        self.d_params.gapless_alpha = local_max_stat_matrix.getAlpha().max(0.0);
        self.d_params.gapless_alpha_error = calculation_error;

        self.d_params.gapless_a = local_max_stat_matrix.getA().max(0.0);
        self.d_params.gapless_a_error = calculation_error;

        self.d_params.G = 0;
        self.d_params.G1 = 0;
        self.d_params.G2 = 0;

        self.d_params.lambda = local_max_stat_matrix.getLambda();
        self.d_params.lambda_error = calculation_error;

        self.d_params.K = local_max_stat_matrix.getK();
        self.d_params.K_error = calculation_error;

        self.d_params.C = local_max_stat_matrix.getC();
        self.d_params.C_error = calculation_error;

        self.d_params.sigma = self.d_params.gapless_alpha;
        self.d_params.sigma_error = calculation_error;

        self.d_params.alpha_I = self.d_params.gapless_alpha;
        self.d_params.alpha_I_error = calculation_error;

        self.d_params.alpha_J = self.d_params.gapless_alpha;
        self.d_params.alpha_J_error = calculation_error;

        self.d_params.a_I = self.d_params.gapless_a;
        self.d_params.a_I_error = calculation_error;

        self.d_params.a_J = self.d_params.gapless_a;
        self.d_params.a_J_error = calculation_error;

        self.d_params.m_LambdaSbs = vec![
            self.d_params.lambda,
            self.d_params.lambda + calculation_error,
        ];
        self.d_params.m_KSbs = vec![self.d_params.K, self.d_params.K + calculation_error];
        self.d_params.m_CSbs = vec![self.d_params.C, self.d_params.C + calculation_error];
        self.d_params.m_SigmaSbs =
            vec![self.d_params.sigma, self.d_params.sigma + calculation_error];
        self.d_params.m_AlphaISbs = vec![
            self.d_params.alpha_I,
            self.d_params.alpha_I + calculation_error,
        ];
        self.d_params.m_AlphaJSbs = vec![
            self.d_params.alpha_J,
            self.d_params.alpha_J + calculation_error,
        ];
        self.d_params.m_AISbs = vec![self.d_params.a_I, self.d_params.a_I + calculation_error];
        self.d_params.m_AJSbs = vec![self.d_params.a_J, self.d_params.a_J + calculation_error];

        self.d_params.a = (self.d_params.a_I + self.d_params.a_J) * 0.5;
        self.d_params.a_error = (self.d_params.a_I_error + self.d_params.a_J_error) * 0.5;

        self.d_params.alpha = (self.d_params.alpha_I + self.d_params.alpha_J) * 0.5;
        self.d_params.alpha_error =
            (self.d_params.alpha_I_error + self.d_params.alpha_J_error) * 0.5;

        self.d_params.d_params_flag = true;
        pvalues::compute_intercepts(&mut self.d_params)?;

        if !pvalues::assert_Gumbel_parameters(&self.d_params) || !self.isGood() {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                "Error - computation of the Gumbel parameters is unsuccessful in the function \"void AlignmentEvaluer::initGapless\"\n".to_string(),
                1,
            ));
        }

        let mut CurrentTime2 = 0.0;
        get_current_time(&mut CurrentTime2);
        self.d_params.m_CalcTime = CurrentTime2 - CurrentTime1;

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn initGapped(
        &mut self,
        alphabetSize_: i64,
        substitutionScoreMatrix_: &[Vec<i64>],
        letterFreqs1_: &[f64],
        letterFreqs2_: &[f64],
        gapOpen1_: i64,
        gapEpen1_: i64,
        gapOpen2_: i64,
        gapEpen2_: i64,
        insertions_after_deletions_: bool,
        eps_lambda_: f64,
        eps_K_: f64,
        max_time_: f64,
        max_mem_: f64,
        randomSeed_: i64,
        temperature_: f64,
    ) -> Result<(), Error> {
        let mut CurrentTime1 = 0.0;
        get_current_time(&mut CurrentTime1);

        let function_name = "void AlignmentEvaluer::initGapped";
        let (letterFreqs1_normalized, letterFreqs2_normalized) = self
            .assert_Gapless_input_parameters(
                alphabetSize_,
                letterFreqs1_,
                letterFreqs2_,
                function_name,
            )?;

        if !(gapEpen1_ > 0) {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                format!(
                    "Error - the parameter \"gapEpen1_\" in the function \"{}\" must be positive\n",
                    function_name
                ),
                1,
            ));
        }
        if !(gapEpen2_ > 0) {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                format!(
                    "Error - the parameter \"gapEpen2_\" in the function \"{}\" must be positive\n",
                    function_name
                ),
                1,
            ));
        }
        if !(eps_lambda_ > 0.0) {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                format!(
                    "Error - the parameter \"eps_lambda_\" in the function \"{}\" must be positive\n",
                    function_name
                ),
                1,
            ));
        }
        if !(eps_K_ > 0.0) {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                format!(
                    "Error - the parameter \"eps_K_\" in the function \"{}\" must be positive\n",
                    function_name
                ),
                1,
            ));
        }
        if !(max_mem_ > 0.0) {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                format!(
                    "Error - the parameter \"max_mem_\" in the function \"{}\" must be positive\n",
                    function_name
                ),
                1,
            ));
        }

        self.d_params.d_params_flag = false;

        let GaplessTimePortion = 0.5;
        let mut GaplessCalculationTime = max_time_;
        if max_time_ <= 0.0 {
            GaplessCalculationTime = 120.0;
        }
        GaplessCalculationTime *= GaplessTimePortion;

        let local_max_stat_matrix = LocalMaxStatMatrix::new(
            alphabetSize_ as usize,
            Some(substitutionScoreMatrix_),
            Some(&letterFreqs1_normalized),
            Some(&letterFreqs2_normalized),
            alphabetSize_ as usize,
            GaplessCalculationTime,
        );

        if local_max_stat_matrix.getTerminated() {
            return Err(Error::new(
                "Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n".to_string(),
                3,
            ));
        }

        let calculation_error = 1e-6;
        self.d_params.gapless_alpha = local_max_stat_matrix.getAlpha().max(0.0);
        self.d_params.gapless_alpha_error = calculation_error;
        self.d_params.gapless_a = local_max_stat_matrix.getA().max(0.0);
        self.d_params.gapless_a_error = calculation_error;

        let mut CurrentTimeGaplessPreliminary = 0.0;
        get_current_time(&mut CurrentTimeGaplessPreliminary);
        let GaplessPreliminaryTime = CurrentTimeGaplessPreliminary - CurrentTime1;

        let gapEpen = gapEpen1_.min(gapEpen2_);
        let gapOpen = (gapOpen1_ + gapEpen1_).min(gapOpen2_ + gapEpen2_) - gapEpen;

        let mut randomization_parameters: Option<struct_for_randomization> = None;
        if max_time_ <= 0.0 {
            if self.d_gapped_computation_parameters.d_parameters_flag {
                randomization_parameters = Some(struct_for_randomization {
                    d_first_stage_preliminary_realizations_numbers_ALP: self
                        .d_gapped_computation_parameters
                        .d_first_stage_preliminary_realizations_numbers_ALP
                        .clone(),
                    d_preliminary_realizations_numbers_ALP: self
                        .d_gapped_computation_parameters
                        .d_preliminary_realizations_numbers_ALP
                        .clone(),
                    d_preliminary_realizations_numbers_killing: self
                        .d_gapped_computation_parameters
                        .d_preliminary_realizations_numbers_killing
                        .clone(),
                    d_random_seed: randomSeed_,
                    d_total_realizations_number_with_ALP: self
                        .d_gapped_computation_parameters
                        .d_total_realizations_number_with_ALP,
                    d_total_realizations_number_with_killing: self
                        .d_gapped_computation_parameters
                        .d_total_realizations_number_with_killing,
                });
            } else {
                return Err(Error::new(
                    "Error - d_gapped_computation_parameters must be defined before calling AlignmentEvaluer::initGapped with max_time_<=0\n".to_string(),
                    1,
                ));
            }
        }

        let mut data_obj = alp_data::new_from_parameters(
            randomSeed_,
            randomization_parameters.as_ref(),
            gapOpen,
            gapOpen1_,
            gapOpen2_,
            gapEpen,
            gapEpen1_,
            gapEpen2_,
            alphabetSize_,
            substitutionScoreMatrix_,
            &letterFreqs1_normalized,
            &letterFreqs2_normalized,
            temperature_,
            max_time_,
            max_mem_,
            eps_lambda_,
            eps_K_,
            insertions_after_deletions_,
            self.d_gapped_computation_parameters
                .d_max_time_for_quick_tests,
            self.d_gapped_computation_parameters
                .d_max_time_with_computation_parameters,
        )?;

        data_obj.d_max_time = ((1.0 - GaplessTimePortion) * data_obj.d_max_time)
            .max(data_obj.d_max_time - GaplessPreliminaryTime);

        let mut sim_obj = alp_sim::new_with_simulation(data_obj)?;

        if max_time_ > 0.0 {
            if let Some(rand_all) = sim_obj.d_alp_data.d_rand_all.as_ref() {
                self.d_gapped_computation_parameters.d_parameters_flag = true;
                self.d_gapped_computation_parameters
                    .d_first_stage_preliminary_realizations_numbers_ALP = rand_all
                    .d_first_stage_preliminary_realizations_numbers_ALP
                    .clone();
                self.d_gapped_computation_parameters
                    .d_preliminary_realizations_numbers_ALP =
                    rand_all.d_preliminary_realizations_numbers_ALP.clone();
                self.d_gapped_computation_parameters
                    .d_preliminary_realizations_numbers_killing =
                    rand_all.d_preliminary_realizations_numbers_killing.clone();
                self.d_gapped_computation_parameters
                    .d_total_realizations_number_with_ALP =
                    rand_all.d_total_realizations_number_with_ALP;
                self.d_gapped_computation_parameters
                    .d_total_realizations_number_with_killing =
                    rand_all.d_total_realizations_number_with_killing;
            }
        }

        sim_obj.m_GaplessAlpha = self.d_params.gapless_alpha;
        sim_obj.m_GaplessAlphaError = self.d_params.gapless_alpha_error;
        sim_obj.m_GaplessA = self.d_params.gapless_a;
        sim_obj.m_GaplessAError = self.d_params.gapless_a_error;

        sim_obj.m_G1 = gapOpen1_ + gapEpen1_;
        sim_obj.m_G2 = gapOpen2_ + gapEpen2_;
        sim_obj.m_G = sim_obj.m_G1.min(sim_obj.m_G2);

        self.d_params.G = sim_obj.m_G;
        self.d_params.G1 = sim_obj.m_G1;
        self.d_params.G2 = sim_obj.m_G2;

        self.d_params.lambda = sim_obj.m_Lambda;
        self.d_params.lambda_error = sim_obj.m_LambdaError;
        self.d_params.K = sim_obj.m_K;
        self.d_params.K_error = sim_obj.m_KError;
        self.d_params.C = sim_obj.m_C;
        self.d_params.C_error = sim_obj.m_CError;
        self.d_params.sigma = sim_obj.m_Sigma;
        self.d_params.sigma_error = sim_obj.m_SigmaError;
        self.d_params.alpha_I = sim_obj.m_AlphaI;
        self.d_params.alpha_I_error = sim_obj.m_AlphaIError;
        self.d_params.alpha_J = sim_obj.m_AlphaJ;
        self.d_params.alpha_J_error = sim_obj.m_AlphaJError;
        self.d_params.a_I = sim_obj.m_AI;
        self.d_params.a_I_error = sim_obj.m_AIError;
        self.d_params.a_J = sim_obj.m_AJ;
        self.d_params.a_J_error = sim_obj.m_AJError;

        self.d_params.m_LambdaSbs = sim_obj.m_LambdaSbs.clone();
        self.d_params.m_KSbs = sim_obj.m_KSbs.clone();
        self.d_params.m_CSbs = sim_obj.m_CSbs.clone();
        self.d_params.m_SigmaSbs = sim_obj.m_SigmaSbs.clone();
        self.d_params.m_AlphaISbs = sim_obj.m_AlphaISbs.clone();
        self.d_params.m_AlphaJSbs = sim_obj.m_AlphaJSbs.clone();
        self.d_params.m_AISbs = sim_obj.m_AISbs.clone();
        self.d_params.m_AJSbs = sim_obj.m_AJSbs.clone();

        self.d_params.a = (self.d_params.a_I + self.d_params.a_J) * 0.5;
        self.d_params.a_error = (self.d_params.a_I_error + self.d_params.a_J_error) * 0.5;
        self.d_params.alpha = (self.d_params.alpha_I + self.d_params.alpha_J) * 0.5;
        self.d_params.alpha_error =
            (self.d_params.alpha_I_error + self.d_params.alpha_J_error) * 0.5;

        self.d_params.d_params_flag = true;
        pvalues::compute_intercepts(&mut self.d_params)?;

        let mut CurrentTime2 = 0.0;
        get_current_time(&mut CurrentTime2);
        self.d_params.m_CalcTime = CurrentTime2 - CurrentTime1;

        if !pvalues::assert_Gumbel_parameters(&self.d_params) || !self.isGood() {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                "Error - computation of the Gumbel parameters is unsuccessful in the function \"void AlignmentEvaluer::initGapped\"\n".to_string(),
                1,
            ));
        }

        Ok(())
    }

    pub fn initParameters(
        &mut self,
        parameters_: &AlignmentEvaluerParameters,
    ) -> Result<(), Error> {
        self.d_params.d_params_flag = false;
        let calculation_error = 1e-6;

        self.d_params.lambda = parameters_.d_lambda;
        self.d_params.lambda_error = calculation_error;
        self.d_params.C = 0.0;
        self.d_params.C_error = 0.0;
        self.d_params.K = parameters_.d_k;
        self.d_params.K_error = calculation_error;
        self.d_params.a_I = parameters_.d_a2;
        self.d_params.a_I_error = calculation_error;
        self.d_params.a_J = parameters_.d_a1;
        self.d_params.a_J_error = calculation_error;
        self.d_params.sigma = parameters_.d_sigma;
        self.d_params.sigma_error = calculation_error;
        self.d_params.alpha_I = parameters_.d_alpha2;
        self.d_params.alpha_I_error = calculation_error;
        self.d_params.alpha_J = parameters_.d_alpha1;
        self.d_params.alpha_J_error = calculation_error;
        self.d_params.a = 0.5 * (parameters_.d_a1 + parameters_.d_a2);
        self.d_params.a_error = calculation_error;
        self.d_params.alpha = 0.5 * (parameters_.d_alpha1 + parameters_.d_alpha2);
        self.d_params.alpha_error = calculation_error;
        self.d_params.gapless_a = 0.0;
        self.d_params.gapless_a_error = 0.0;
        self.d_params.gapless_alpha = 0.0;
        self.d_params.gapless_alpha_error = 0.0;
        self.d_params.G = 0;
        self.d_params.G1 = 0;
        self.d_params.G2 = 0;
        self.d_params.b_I = parameters_.d_b2;
        self.d_params.b_I_error = calculation_error;
        self.d_params.b_J = parameters_.d_b1;
        self.d_params.b_J_error = calculation_error;
        self.d_params.beta_I = parameters_.d_beta2;
        self.d_params.beta_I_error = calculation_error;
        self.d_params.beta_J = parameters_.d_beta1;
        self.d_params.beta_J_error = calculation_error;
        self.d_params.tau = parameters_.d_tau;
        self.d_params.tau_error = calculation_error;

        self.d_params.m_LambdaSbs = vec![
            self.d_params.lambda,
            self.d_params.lambda + calculation_error,
        ];
        self.d_params.m_KSbs = vec![self.d_params.K, self.d_params.K + calculation_error];
        self.d_params.m_CSbs = vec![self.d_params.C, self.d_params.C + calculation_error];
        self.d_params.m_SigmaSbs =
            vec![self.d_params.sigma, self.d_params.sigma + calculation_error];
        self.d_params.m_AlphaISbs = vec![
            self.d_params.alpha_I,
            self.d_params.alpha_I + calculation_error,
        ];
        self.d_params.m_AlphaJSbs = vec![
            self.d_params.alpha_J,
            self.d_params.alpha_J + calculation_error,
        ];
        self.d_params.m_AISbs = vec![self.d_params.a_I, self.d_params.a_I + calculation_error];
        self.d_params.m_AJSbs = vec![self.d_params.a_J, self.d_params.a_J + calculation_error];
        self.d_params.m_BJSbs = vec![self.d_params.b_J, self.d_params.b_J + calculation_error];
        self.d_params.m_BISbs = vec![self.d_params.b_I, self.d_params.b_I + calculation_error];
        self.d_params.m_BetaJSbs = vec![
            self.d_params.beta_J,
            self.d_params.beta_J + calculation_error,
        ];
        self.d_params.m_BetaISbs = vec![
            self.d_params.beta_I,
            self.d_params.beta_I + calculation_error,
        ];
        self.d_params.m_TauSbs = vec![self.d_params.tau, self.d_params.tau + calculation_error];

        self.d_params.d_params_flag = true;
        pvalues::compute_tmp_values(&mut self.d_params)?;
        self.d_params.m_CalcTime = 0.0;

        if !pvalues::assert_Gumbel_parameters(&self.d_params) || !self.isGood() {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                "Error - computation of the Gumbel parameters is unsuccessful in the function \"void AlignmentEvaluer::initParameters\"\n".to_string(),
                1,
            ));
        }
        Ok(())
    }

    pub fn initParametersWithErrors(
        &mut self,
        parameters_: &AlignmentEvaluerParametersWithErrors,
    ) -> Result<(), Error> {
        self.d_params.d_params_flag = false;
        let array_dim = 20usize;
        unsafe {
            srand(12345);
        }

        self.d_params.lambda = parameters_.d_lambda;
        self.d_params.lambda_error = parameters_.d_lambda_error;
        self.d_params.C = 0.0;
        self.d_params.C_error = 0.0;
        self.d_params.K = parameters_.d_k;
        self.d_params.K_error = parameters_.d_k_error;
        self.d_params.a_I = parameters_.d_a2;
        self.d_params.a_I_error = parameters_.d_a2_error;
        self.d_params.a_J = parameters_.d_a1;
        self.d_params.a_J_error = parameters_.d_a1_error;
        self.d_params.sigma = parameters_.d_sigma;
        self.d_params.sigma_error = parameters_.d_sigma_error;
        self.d_params.alpha_I = parameters_.d_alpha2;
        self.d_params.alpha_I_error = parameters_.d_alpha2_error;
        self.d_params.alpha_J = parameters_.d_alpha1;
        self.d_params.alpha_J_error = parameters_.d_alpha1_error;
        self.d_params.a = 0.5 * (parameters_.d_a1 + parameters_.d_a2);
        self.d_params.a_error = 0.5 * (parameters_.d_a1_error + parameters_.d_a2_error);
        self.d_params.alpha = 0.5 * (parameters_.d_alpha1 + parameters_.d_alpha2);
        self.d_params.alpha_error = 0.5 * (parameters_.d_alpha1_error + parameters_.d_alpha2_error);
        self.d_params.gapless_a = 0.0;
        self.d_params.gapless_a_error = 0.0;
        self.d_params.gapless_alpha = 0.0;
        self.d_params.gapless_alpha_error = 0.0;
        self.d_params.G = 0;
        self.d_params.G1 = 0;
        self.d_params.G2 = 0;
        self.d_params.m_CalcTime = 0.0;
        self.d_params.b_I = parameters_.d_b2;
        self.d_params.b_I_error = parameters_.d_b2_error;
        self.d_params.b_J = parameters_.d_b1;
        self.d_params.b_J_error = parameters_.d_b1_error;
        self.d_params.beta_I = parameters_.d_beta2;
        self.d_params.beta_I_error = parameters_.d_beta2_error;
        self.d_params.beta_J = parameters_.d_beta1;
        self.d_params.beta_J_error = parameters_.d_beta1_error;
        self.d_params.tau = parameters_.d_tau;
        self.d_params.tau_error = parameters_.d_tau_error;

        self.d_params.m_LambdaSbs.clear();
        self.d_params.m_KSbs.clear();
        self.d_params.m_CSbs.clear();
        self.d_params.m_SigmaSbs.clear();
        self.d_params.m_AlphaISbs.clear();
        self.d_params.m_AlphaJSbs.clear();
        self.d_params.m_AISbs.clear();
        self.d_params.m_AJSbs.clear();
        self.d_params.m_BISbs.clear();
        self.d_params.m_BJSbs.clear();
        self.d_params.m_BetaISbs.clear();
        self.d_params.m_BetaJSbs.clear();
        self.d_params.m_TauSbs.clear();

        let sqrt_array_dim = (array_dim as f64).sqrt();
        for _ in 0..array_dim {
            self.d_params.m_LambdaSbs.push(
                self.d_params.lambda
                    + self.d_params.lambda_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_KSbs.push(
                self.d_params.K
                    + self.d_params.K_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_CSbs.push(0.0);
            self.d_params.m_SigmaSbs.push(
                self.d_params.sigma
                    + self.d_params.sigma_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_AlphaISbs.push(
                self.d_params.alpha_I
                    + self.d_params.alpha_I_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_AlphaJSbs.push(
                self.d_params.alpha_J
                    + self.d_params.alpha_J_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_AISbs.push(
                self.d_params.a_I
                    + self.d_params.a_I_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_AJSbs.push(
                self.d_params.a_J
                    + self.d_params.a_J_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_BISbs.push(
                self.d_params.b_I
                    + self.d_params.b_I_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_BJSbs.push(
                self.d_params.b_J
                    + self.d_params.b_J_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_BetaISbs.push(
                self.d_params.beta_I
                    + self.d_params.beta_I_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_BetaJSbs.push(
                self.d_params.beta_J
                    + self.d_params.beta_J_error * pvalues::standard_normal() * sqrt_array_dim,
            );
            self.d_params.m_TauSbs.push(
                self.d_params.tau
                    + self.d_params.tau_error * pvalues::standard_normal() * sqrt_array_dim,
            );
        }

        self.d_params.d_params_flag = true;
        pvalues::compute_tmp_values(&mut self.d_params)?;
        if !pvalues::assert_Gumbel_parameters(&self.d_params) || !self.isGood() {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                "Error - computation of the Gumbel parameters is unsuccessful in the function \"void AlignmentEvaluer::initParameters\"\n".to_string(),
                1,
            ));
        }
        Ok(())
    }

    pub fn calc(
        &self,
        score_: f64,
        seqlen1_: f64,
        seqlen2_: f64,
        pvalue_: &mut f64,
        pvalueErr_: &mut f64,
        evalue_: &mut f64,
        evalueErr_: &mut f64,
    ) -> Result<(), Error> {
        if seqlen1_ <= 0.0 || seqlen2_ <= 0.0 {
            return Err(Error::new(
                "Error - seqlen1_<=0 or seqlen2_<=0 in \"double AlignmentEvaluer::calc\"\n"
                    .to_string(),
                2,
            ));
        }
        if !self.isGood() {
            return Err(Error::new(
                "Unexpected error - the Gumbel parameters are not defined properly in \"double AlignmentEvaluer::calc\"\n".to_string(),
                1,
            ));
        }

        let pvalues_obj = pvalues::new();
        pvalues_obj.calculate_P_values(
            score_,
            seqlen2_,
            seqlen1_,
            &self.d_params,
            pvalue_,
            pvalueErr_,
            evalue_,
            evalueErr_,
            true,
        )
    }

    pub fn calc_without_errors(
        &self,
        score_: f64,
        seqlen1_: f64,
        seqlen2_: f64,
        pvalue_: &mut f64,
        evalue_: &mut f64,
    ) -> Result<(), Error> {
        if seqlen1_ <= 0.0 || seqlen2_ <= 0.0 {
            return Err(Error::new(
                "Error - seqlen1_<=0 or seqlen2_<=0 in \"double AlignmentEvaluer::calc\"\n"
                    .to_string(),
                2,
            ));
        }
        if !self.isGood() {
            return Err(Error::new(
                "Unexpected error - d_params is not defined in \"double AlignmentEvaluer::calc\"\n"
                    .to_string(),
                1,
            ));
        }

        let pvalues_obj = pvalues::new();
        let mut area = 0.0;
        let mut area_is_1_flag = false;
        pvalues::get_appr_tail_prob_with_cov_without_errors(
            &self.d_params,
            pvalues_obj.blast,
            score_,
            seqlen2_,
            seqlen1_,
            pvalue_,
            evalue_,
            &mut area,
            &mut area_is_1_flag,
            false,
        );
        Ok(())
    }

    pub fn evalue(&self, score_: f64, seqlen1_: f64, seqlen2_: f64) -> Result<f64, Error> {
        Ok(self.area(score_, seqlen1_, seqlen2_)? * self.evaluePerArea(score_))
    }

    pub fn pvalue(evalue_: f64) -> f64 {
        one_minus_exp_function(-evalue_)
    }

    pub fn area(&self, score_: f64, seqlen1_: f64, seqlen2_: f64) -> Result<f64, Error> {
        if seqlen1_ <= 0.0 || seqlen2_ <= 0.0 {
            return Err(Error::new(
                "Error - seqlen1_<=0 or seq2en1_<=0 in \"double AlignmentEvaluer::area\"\n"
                    .to_string(),
                2,
            ));
        }
        if !self.isGood() {
            return Err(Error::new(
                "Unexpected error - the Gumbel parameters are not defined properly in \"double AlignmentEvaluer::area\"\n".to_string(),
                1,
            ));
        }

        let pvalues_obj = pvalues::new();
        let mut P = 0.0;
        let mut E = 0.0;
        let mut area_res = 0.0;
        let mut area_is_1_flag = false;
        pvalues::get_appr_tail_prob_with_cov_without_errors(
            &self.d_params,
            pvalues_obj.blast,
            score_,
            seqlen2_,
            seqlen1_,
            &mut P,
            &mut E,
            &mut area_res,
            &mut area_is_1_flag,
            true,
        );
        Ok(area_res)
    }

    pub fn log_area(&self, score_: f64, seqlen1_: f64, seqlen2_: f64) -> Result<f64, Error> {
        if seqlen1_ <= 0.0 || seqlen2_ <= 0.0 {
            return Err(Error::new(
                "Error - seqlen1_<=0 or seq2en1_<=0 in \"double AlignmentEvaluer::area\"\n"
                    .to_string(),
                2,
            ));
        }
        if !self.isGood() {
            return Err(Error::new(
                "Unexpected error - the Gumbel parameters are not defined properly in \"double AlignmentEvaluer::area\"\n".to_string(),
                1,
            ));
        }
        let pvalues_obj = pvalues::new();
        Ok(pvalues_obj.log_area(
            &self.d_params,
            pvalues_obj.blast,
            score_,
            seqlen2_,
            seqlen1_,
        ))
    }

    pub fn evaluePerArea(&self, score_: f64) -> f64 {
        self.d_params.K * (-self.d_params.lambda * score_).exp()
    }

    pub fn bitScore(&self, score_: f64) -> f64 {
        (self.d_params.lambda * score_ - self.d_params.K.ln()) / 2.0f64.ln()
    }

    pub fn isGood(&self) -> bool {
        self.d_params.d_params_flag
    }

    pub fn parameters(&self) -> &ALP_set_of_parameters {
        &self.d_params
    }

    pub fn gapped_computation_parameters(&self) -> &gapped_computation_parameters_struct {
        &self.d_gapped_computation_parameters
    }

    pub fn set_gapped_computation_parameters_simplified(
        &mut self,
        max_time_: f64,
        number_of_samples_: i64,
        number_of_samples_for_preliminary_stages_: i64,
    ) -> Result<(), Error> {
        if number_of_samples_ <= 0 || number_of_samples_for_preliminary_stages_ <= 0 {
            return Err(Error::new(
                "Error - number_of_samples_<=0 or number_of_samples_for_preliminary_stages_<=0 in \"void AlignmentEvaluer::set_gapped_computation_parameters_simplified\"\n".to_string(),
                2,
            ));
        }

        self.d_gapped_computation_parameters
            .d_first_stage_preliminary_realizations_numbers_ALP =
            vec![number_of_samples_for_preliminary_stages_];
        self.d_gapped_computation_parameters
            .d_preliminary_realizations_numbers_ALP =
            vec![number_of_samples_for_preliminary_stages_];
        self.d_gapped_computation_parameters
            .d_preliminary_realizations_numbers_killing =
            vec![number_of_samples_for_preliminary_stages_];
        self.d_gapped_computation_parameters
            .d_total_realizations_number_with_ALP = number_of_samples_;
        self.d_gapped_computation_parameters
            .d_total_realizations_number_with_killing = number_of_samples_;
        self.d_gapped_computation_parameters.d_parameters_flag = true;
        self.d_gapped_computation_parameters
            .d_max_time_with_computation_parameters = max_time_;
        if max_time_ > 0.0 {
            self.d_gapped_computation_parameters
                .d_max_time_for_quick_tests = 0.5 * max_time_ * QUICK_TESTS_TRIALS_NUMBER as f64
                / (QUICK_TESTS_TRIALS_NUMBER
                    + number_of_samples_
                    + number_of_samples_for_preliminary_stages_) as f64;
        } else {
            self.d_gapped_computation_parameters
                .d_max_time_for_quick_tests = -1.0;
        }
        Ok(())
    }

    pub fn set_gapped_computation_parameters(
        &mut self,
        gapped_computation_parameters_: &gapped_computation_parameters_struct,
    ) {
        self.d_gapped_computation_parameters = gapped_computation_parameters_.clone();
    }

    pub fn gapped_computation_parameters_clear(&mut self) {
        self.d_gapped_computation_parameters
            .d_first_stage_preliminary_realizations_numbers_ALP
            .clear();
        self.d_gapped_computation_parameters
            .d_preliminary_realizations_numbers_ALP
            .clear();
        self.d_gapped_computation_parameters
            .d_preliminary_realizations_numbers_killing
            .clear();
        self.d_gapped_computation_parameters.d_parameters_flag = false;
        self.d_gapped_computation_parameters
            .d_max_time_for_quick_tests = -1.0;
        self.d_gapped_computation_parameters
            .d_max_time_with_computation_parameters = -1.0;
    }

    pub fn assert_Gapless_input_parameters(
        &mut self,
        alphabetSize_: i64,
        letterFreqs1_: &[f64],
        letterFreqs2_: &[f64],
        function_name_: &str,
    ) -> Result<(Vec<f64>, Vec<f64>), Error> {
        if alphabetSize_ <= 0 {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                format!(
                    "Error - the parameter \"alphabetSize_\" in the function \"{}\" must be positive\n",
                    function_name_
                ),
                1,
            ));
        }
        let alphabet_size = alphabetSize_ as usize;
        if letterFreqs1_.len() < alphabet_size || letterFreqs2_.len() < alphabet_size {
            self.d_params.d_params_flag = false;
            return Err(Error::new(
                format!(
                    "Error - frequency arrays are shorter than \"alphabetSize_\" in the function \"{}\"\n",
                    function_name_
                ),
                1,
            ));
        }

        let mut sum1 = 0.0;
        for (i, &freq) in letterFreqs1_.iter().take(alphabet_size).enumerate() {
            if freq < 0.0 {
                self.d_params.d_params_flag = false;
                return Err(Error::new(
                    format!(
                        "Error - the value \"letterFreqs1_[{}]\" in the function \"{}\" must be non-negative\n",
                        i, function_name_
                    ),
                    1,
                ));
            }
            sum1 += freq;
        }
        if sum1 <= 0.0 {
            return Err(Error::new(
                format!(
                    "Error - sum of the frequencies \"letterFreqs1_\" is non-positive in the function \"{}\"\n",
                    function_name_
                ),
                1,
            ));
        }
        let letterFreqs1_normalized = letterFreqs1_[..alphabet_size]
            .iter()
            .map(|&x| x / sum1)
            .collect::<Vec<_>>();

        let mut sum2 = 0.0;
        for (i, &freq) in letterFreqs2_.iter().take(alphabet_size).enumerate() {
            if freq < 0.0 {
                self.d_params.d_params_flag = false;
                return Err(Error::new(
                    format!(
                        "Error - the value \"letterFreqs2_[{}]\" in the function \"{}\" must be non-negative\n",
                        i, function_name_
                    ),
                    1,
                ));
            }
            sum2 += freq;
        }
        if sum2 <= 0.0 {
            return Err(Error::new(
                format!(
                    "Error - sum of the frequencies \"letterFreqs2_\" is non-positive in the function \"{}\"\n",
                    function_name_
                ),
                1,
            ));
        }
        let letterFreqs2_normalized = letterFreqs2_[..alphabet_size]
            .iter()
            .map(|&x| x / sum2)
            .collect::<Vec<_>>();

        Ok((letterFreqs1_normalized, letterFreqs2_normalized))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parameters() -> AlignmentEvaluerParameters {
        AlignmentEvaluerParameters {
            d_lambda: 0.267,
            d_k: 0.041,
            d_a1: 1.9,
            d_b1: 4.0,
            d_a2: 2.1,
            d_b2: 5.0,
            d_alpha1: 1.7,
            d_beta1: 6.0,
            d_alpha2: 1.8,
            d_beta2: 7.0,
            d_sigma: 43.0,
            d_tau: 8.0,
        }
    }

    #[test]
    fn test_init_parameters_and_scoring_methods() {
        let mut ev = AlignmentEvaluer::new();
        ev.initParameters(&parameters()).unwrap();
        assert!(ev.isGood());
        assert_eq!(ev.parameters().a_I, 2.1);
        assert_eq!(ev.parameters().a_J, 1.9);
        assert_eq!(ev.parameters().b_I, 5.0);
        assert_eq!(ev.parameters().b_J, 4.0);
        assert_eq!(ev.parameters().m_LambdaSbs.len(), 2);

        let area = ev.area(50.0, 1000.0, 2000.0).unwrap();
        let log_area = ev.log_area(50.0, 1000.0, 2000.0).unwrap();
        assert!(area > 0.0);
        assert!((log_area.exp() - area).abs() / area < 1e-12);
        assert_eq!(ev.evaluePerArea(50.0), 0.041 * (-0.267f64 * 50.0).exp());
        assert!(ev.bitScore(50.0) > 0.0);
        assert!(AlignmentEvaluer::pvalue(0.5) > 0.0);
    }

    #[test]
    fn test_calc_with_and_without_errors() {
        let mut ev = AlignmentEvaluer::new();
        ev.initParameters(&parameters()).unwrap();
        let mut p = 0.0;
        let mut e = 0.0;
        ev.calc_without_errors(50.0, 1000.0, 2000.0, &mut p, &mut e)
            .unwrap();
        assert!(p > 0.0 && e > 0.0);

        let mut p2 = 0.0;
        let mut p2_err = 0.0;
        let mut e2 = 0.0;
        let mut e2_err = 0.0;
        ev.calc(
            50.0,
            1000.0,
            2000.0,
            &mut p2,
            &mut p2_err,
            &mut e2,
            &mut e2_err,
        )
        .unwrap();
        assert!((p2 - p).abs() / p < 1e-12);
        assert!((e2 - e).abs() / e < 1e-12);
        assert!(p2_err >= 0.0 && e2_err >= 0.0);
        assert_eq!(ev.evalue(50.0, 1000.0, 2000.0).unwrap(), e);
    }

    #[test]
    fn test_init_parameters_with_errors_populates_sbs_arrays() {
        let mut ev = AlignmentEvaluer::new();
        ev.initParametersWithErrors(&AlignmentEvaluerParametersWithErrors {
            d_lambda: 0.267,
            d_lambda_error: 0.001,
            d_k: 0.041,
            d_k_error: 0.001,
            d_a1: 1.9,
            d_a1_error: 0.01,
            d_b1: 4.0,
            d_b1_error: 0.01,
            d_a2: 2.1,
            d_a2_error: 0.01,
            d_b2: 5.0,
            d_b2_error: 0.01,
            d_alpha1: 1.7,
            d_alpha1_error: 0.01,
            d_beta1: 6.0,
            d_beta1_error: 0.01,
            d_alpha2: 1.8,
            d_alpha2_error: 0.01,
            d_beta2: 7.0,
            d_beta2_error: 0.01,
            d_sigma: 43.0,
            d_sigma_error: 0.01,
            d_tau: 8.0,
            d_tau_error: 0.01,
        })
        .unwrap();
        assert_eq!(ev.parameters().m_LambdaSbs.len(), 20);
        assert_eq!(ev.parameters().m_TauSbs.len(), 20);
        assert!(pvalues::assert_Gumbel_parameters(ev.parameters()));
    }

    #[test]
    fn test_gapped_computation_parameters_and_frequency_validation() {
        let mut ev = AlignmentEvaluer::new();
        ev.set_gapped_computation_parameters_simplified(60.0, 1000, 200)
            .unwrap();
        let g = ev.gapped_computation_parameters();
        assert!(g.d_parameters_flag);
        assert_eq!(g.d_total_realizations_number_with_ALP, 1000);
        assert!(g.d_max_time_for_quick_tests > 0.0);

        let (f1, f2) = ev
            .assert_Gapless_input_parameters(2, &[2.0, 6.0], &[1.0, 3.0], "test")
            .unwrap();
        assert_eq!(f1, vec![0.25, 0.75]);
        assert_eq!(f2, vec![0.25, 0.75]);

        ev.gapped_computation_parameters_clear();
        assert!(!ev.gapped_computation_parameters().d_parameters_flag);
    }

    #[test]
    fn test_init_gapless_uses_local_max_stat_matrix() {
        let mut ev = AlignmentEvaluer::new();
        let scores = vec![vec![-1, 2], vec![-1, -1]];
        let freqs = [0.5, 0.5];
        ev.initGapless(2, &scores, &freqs, &freqs, 1.0).unwrap();

        assert!(ev.isGood());
        assert!(ev.parameters().lambda > 0.0);
        assert!(ev.parameters().K > 0.0);
        assert_eq!(ev.parameters().G, 0);
        assert_eq!(ev.parameters().m_LambdaSbs.len(), 2);
        assert_eq!(ev.parameters().m_BISbs.len(), 2);
        assert!(pvalues::assert_Gumbel_parameters(ev.parameters()));
    }

    #[test]
    fn test_init_gapped_validates_inputs_before_simulation() {
        let mut ev = AlignmentEvaluer::new();
        let scores = vec![vec![-1, 2], vec![-1, -1]];
        let freqs = [0.5, 0.5];

        let err = ev
            .initGapped(
                2, &scores, &freqs, &freqs, 5, 0, 5, 1, true, 0.05, 0.05, 1.0, 64.0, 1, 1.07,
            )
            .unwrap_err();
        assert!(err.st.contains("\"gapEpen1_\""));

        let err = ev
            .initGapped(
                2, &scores, &freqs, &freqs, 5, 1, 5, 1, true, 0.05, 0.05, 0.0, 64.0, 1, 1.07,
            )
            .unwrap_err();
        assert!(err.st.contains(
            "d_gapped_computation_parameters must be defined before calling AlignmentEvaluer::initGapped"
        ));
    }
}
