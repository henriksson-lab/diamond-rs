#![allow(non_snake_case)]

use std::fmt;
use std::sync::Mutex;

use crate::stats::alp_function;
use crate::stats::alp_localmaxstat_util as util;

static S_TIME: Mutex<f64> = Mutex::new(0.0);

#[derive(Debug, Clone, PartialEq)]
pub struct LocalMaxStat {
    d_dimension: usize,
    d_score_p: Vec<i64>,
    d_prob_p: Vec<f64>,
    d_lambda: f64,
    d_k: f64,
    d_c: f64,
    d_thetaMin: f64,
    d_rMin: f64,
    d_delta: i64,
    d_thetaMinusDelta: f64,
    d_mu: f64,
    d_sigma: f64,
    d_muAssoc: f64,
    d_sigmaAssoc: f64,
    d_meanWDLE: f64,
    d_terminated: bool,
}

impl Default for LocalMaxStat {
    fn default() -> Self {
        Self {
            d_dimension: 0,
            d_score_p: Vec::new(),
            d_prob_p: Vec::new(),
            d_lambda: 0.0,
            d_k: 0.0,
            d_c: 0.0,
            d_thetaMin: 0.0,
            d_rMin: 0.0,
            d_delta: 0,
            d_thetaMinusDelta: 0.0,
            d_mu: 0.0,
            d_sigma: 0.0,
            d_muAssoc: 0.0,
            d_sigmaAssoc: 0.0,
            d_meanWDLE: 0.0,
            d_terminated: false,
        }
    }
}

impl LocalMaxStat {
    pub fn setTime(time_: f64) {
        assert!(time_ >= 0.0);
        *S_TIME.lock().unwrap() = time_;
    }

    pub fn getTime() -> f64 {
        *S_TIME.lock().unwrap()
    }

    pub fn new(dimension_: usize, score_: Option<&[i64]>, prob_: Option<&[f64]>) -> Self {
        let mut this = Self::default();
        if let (Some(score), Some(prob)) = (score_, prob_) {
            this.copy(dimension_, score, prob);
        } else {
            this.copy(0, &[], &[]);
        }
        this
    }

    pub fn bool_(&self) -> bool {
        self.d_dimension != 0
    }

    pub fn assign(&mut self, localMaxStat_: &Self) -> &mut Self {
        if !std::ptr::eq(self, localMaxStat_) {
            self.copy_local(localMaxStat_);
        }
        self
    }

    pub fn copy_local(&mut self, localMaxStat_: &Self) {
        self.copy_full(
            localMaxStat_.getDimension(),
            localMaxStat_.getScore(),
            localMaxStat_.getProb(),
            localMaxStat_.getLambda(),
            localMaxStat_.getK(),
            localMaxStat_.getC(),
            localMaxStat_.getThetaMin(),
            localMaxStat_.getRMin(),
            localMaxStat_.getDelta(),
            localMaxStat_.getThetaMinusDelta(),
            localMaxStat_.getMu(),
            localMaxStat_.getSigma(),
            localMaxStat_.getMuAssoc(),
            localMaxStat_.getSigmaAssoc(),
            localMaxStat_.getMeanWDLE(),
            localMaxStat_.getTerminated(),
        );
    }

    #[allow(clippy::too_many_arguments)]
    pub fn copy_full(
        &mut self,
        dimension_: usize,
        score_: &[i64],
        prob_: &[f64],
        lambda_: f64,
        k_: f64,
        c_: f64,
        thetaMin_: f64,
        rMin_: f64,
        delta_: i64,
        thetaMinusDelta_: f64,
        mu_: f64,
        sigma_: f64,
        muAssoc_: f64,
        sigmaAssoc_: f64,
        meanLength_: f64,
        terminated_: bool,
    ) {
        self.free2();
        self.init(dimension_);
        self.d_score_p.copy_from_slice(&score_[..dimension_]);
        self.d_prob_p.copy_from_slice(&prob_[..dimension_]);
        self.d_lambda = lambda_;
        self.d_k = k_;
        self.d_c = c_;
        self.d_thetaMin = thetaMin_;
        self.d_rMin = rMin_;
        self.d_delta = delta_;
        self.d_thetaMinusDelta = thetaMinusDelta_;
        self.d_mu = mu_;
        self.d_sigma = sigma_;
        self.d_muAssoc = muAssoc_;
        self.d_sigmaAssoc = sigmaAssoc_;
        self.d_meanWDLE = meanLength_;
        self.d_terminated = terminated_;
    }

    pub fn copy(&mut self, dimension_: usize, score_: &[i64], prob_: &[f64]) {
        if dimension_ == 0 {
            self.clear();
            return;
        }
        if !util::isLogarithmic(dimension_, score_, prob_) {
            panic!("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n");
        }

        self.free2();
        self.init(dimension_);
        self.d_score_p.copy_from_slice(&score_[..dimension_]);
        self.d_prob_p.copy_from_slice(&prob_[..dimension_]);

        self.d_mu = util::mu(self.getDimension(), self.getScore(), self.getProb());
        self.d_sigma = 0.0;
        for i in 0..dimension_ {
            self.d_sigma += score_[i] as f64 * score_[i] as f64 * prob_[i];
        }
        self.d_sigma -= self.getMu() * self.getMu();
        self.d_sigma = alp_function::psqrt(self.getSigma());

        self.d_lambda = util::lambda(self.getDimension(), self.getScore(), self.getProb());
        self.d_muAssoc = util::muAssoc(
            self.getDimension(),
            self.getScore(),
            self.getProb(),
            self.getLambda(),
        );
        self.d_sigmaAssoc = 0.0;
        for i in 0..self.getDimension() {
            self.d_sigmaAssoc += self.getScore()[i] as f64
                * self.getScore()[i] as f64
                * self.getProb()[i]
                * (self.getLambda() * self.getScore()[i] as f64).exp();
        }
        self.d_sigmaAssoc -= self.getMuAssoc() * self.getMuAssoc();
        self.d_sigmaAssoc = alp_function::psqrt(self.d_sigmaAssoc);

        self.d_thetaMin = util::thetaMin(
            self.getDimension(),
            self.getScore(),
            self.getProb(),
            self.getLambda(),
        );
        self.d_rMin = util::rMin(
            self.getDimension(),
            self.getScore(),
            self.getProb(),
            self.getLambda(),
            self.getThetaMin(),
        );
        self.d_delta = util::delta(self.getDimension(), self.getScore());
        self.d_thetaMinusDelta =
            util::thetaMinusDelta(self.getLambda(), self.getDimension(), self.getScore());
        self.dynProgCalc();
    }

    pub fn out(&self) -> String {
        String::new()
    }

    pub fn getR(&self, theta_: f64) -> f64 {
        util::r(self.d_dimension, &self.d_score_p, &self.d_prob_p, theta_)
    }

    pub fn getA(&self) -> f64 {
        if self.getMuAssoc() == 0.0 {
            f64::INFINITY
        } else {
            1.0 / self.getMuAssoc()
        }
    }

    pub fn getAlpha(&self) -> f64 {
        self.getSigmaAssoc() * self.getSigmaAssoc() * self.getA() * self.getA() * self.getA()
    }

    pub fn getDimension(&self) -> usize {
        self.d_dimension
    }

    pub fn getScore(&self) -> &[i64] {
        &self.d_score_p
    }

    pub fn getProb(&self) -> &[f64] {
        &self.d_prob_p
    }

    pub fn getLambda(&self) -> f64 {
        self.d_lambda
    }

    pub fn getK(&self) -> f64 {
        self.d_k
    }

    pub fn getC(&self) -> f64 {
        self.d_c
    }

    pub fn getThetaMin(&self) -> f64 {
        self.d_thetaMin
    }

    pub fn getRMin(&self) -> f64 {
        self.d_rMin
    }

    pub fn getDelta(&self) -> i64 {
        self.d_delta
    }

    pub fn getThetaMinusDelta(&self) -> f64 {
        self.d_thetaMinusDelta
    }

    pub fn getMu(&self) -> f64 {
        self.d_mu
    }

    pub fn getSigma(&self) -> f64 {
        self.d_sigma
    }

    pub fn getMuAssoc(&self) -> f64 {
        self.d_muAssoc
    }

    pub fn getSigmaAssoc(&self) -> f64 {
        self.d_sigmaAssoc
    }

    pub fn getMeanWDLE(&self) -> f64 {
        self.d_meanWDLE
    }

    pub fn getTerminated(&self) -> bool {
        self.d_terminated
    }

    fn init(&mut self, dimension_: usize) {
        self.d_score_p = vec![0; dimension_];
        self.d_prob_p = vec![0.0; dimension_];
        self.d_dimension = dimension_;
    }

    fn free2(&mut self) {
        self.d_score_p.clear();
        self.d_prob_p.clear();
        self.d_dimension = 0;
    }

    pub fn clear(&mut self) {
        self.free2();
        self.init(0);
        self.d_lambda = 0.0;
        self.d_k = 0.0;
        self.d_c = 0.0;
        self.d_thetaMin = 0.0;
        self.d_rMin = 0.0;
        self.d_delta = 0;
        self.d_thetaMinusDelta = 0.0;
        self.d_mu = 0.0;
        self.d_sigma = 0.0;
        self.d_muAssoc = 0.0;
        self.d_sigmaAssoc = 0.0;
        self.d_meanWDLE = 0.0;
        self.d_terminated = false;
    }

    fn dynProgCalc(&mut self) {
        let mut eSumAlpha = 0.0;
        let mut eOneMinusExpSumAlpha = 0.0;
        let dimension = self.getDimension();
        let score = self.getScore().to_vec();
        let prob = self.getProb().to_vec();
        util::descendingLadderEpoch(
            dimension,
            &score,
            &prob,
            Some(&mut eSumAlpha),
            Some(&mut eOneMinusExpSumAlpha),
            false,
            self.getLambda(),
            self.getMu(),
            self.getMuAssoc(),
            self.getThetaMin(),
            self.getRMin(),
            Self::getTime(),
            Some(&mut self.d_terminated),
        );

        if self.getTerminated() {
            return;
        }

        let ratio = eOneMinusExpSumAlpha / eSumAlpha;
        self.d_k = self.getMu() * self.getMu() / self.getThetaMinusDelta() / self.getMuAssoc()
            * ratio
            * ratio;
        self.d_meanWDLE = eSumAlpha / self.getMu();
        self.d_c = self.getK() * self.getMeanWDLE() / eOneMinusExpSumAlpha;
    }
}

impl fmt::Display for LocalMaxStat {
    fn fmt(&self, ostr_: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(ostr_, "{}", self.out())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static TEST_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn test_local_max_stat_empty_and_full_copy() {
        let mut stat = LocalMaxStat::new(0, None, None);
        assert!(!stat.bool_());
        assert_eq!(stat.out(), "");
        assert_eq!(format!("{}", stat), "");
        stat.copy_full(
            2,
            &[-1, 2],
            &[0.8, 0.2],
            0.4,
            1.0,
            2.0,
            0.1,
            0.9,
            1,
            0.3,
            -0.4,
            1.0,
            0.5,
            0.6,
            3.0,
            true,
        );
        assert!(stat.bool_());
        assert_eq!(stat.getDimension(), 2);
        assert_eq!(stat.getScore(), &[-1, 2]);
        assert_eq!(stat.getProb(), &[0.8, 0.2]);
        assert_eq!(stat.getA(), 2.0);
        assert_eq!(stat.getAlpha(), 2.88);
        assert!(stat.getTerminated());

        let mut copied = LocalMaxStat::default();
        copied.copy_local(&stat);
        assert_eq!(copied, stat);
    }

    #[test]
    fn test_local_max_stat_computes_small_distribution() {
        let _guard = TEST_LOCK.lock().unwrap();
        LocalMaxStat::setTime(0.0);
        let stat = LocalMaxStat::new(2, Some(&[-1, 2]), Some(&[0.8, 0.2]));
        assert!(stat.bool_());
        assert!((stat.getMu() + 0.4).abs() < 1e-12);
        assert!((stat.getLambda() - 0.44568071901268186).abs() < 1e-6);
        assert!((stat.getR(stat.getLambda()) - 1.0).abs() < 1e-6);
        assert_eq!(stat.getDelta(), 1);
        assert!(stat.getK().is_finite());
        assert!(stat.getC().is_finite());
        assert!(stat.getMeanWDLE().is_finite());
    }

    #[test]
    #[should_panic(expected = "Error - you have exceeded the calculation time or memory limit.")]
    fn test_local_max_stat_rejects_non_logarithmic_distribution() {
        let _guard = TEST_LOCK.lock().unwrap();
        let _ = LocalMaxStat::new(2, Some(&[-1, 2]), Some(&[0.2, 0.8]));
    }
}
