#![allow(non_snake_case)]

use std::sync::Mutex;

use crate::stats::alp_approx;
use crate::stats::alp_dynprogproblim::DynProgProbLim;
use crate::stats::alp_integer;
use crate::stats::alp_root;
use crate::stats::sls_basic;

pub const REL_TOL: f64 = 1.0e-6;

#[derive(Debug, Clone, Default)]
struct LocalParams {
    dimension: usize,
    score: Vec<i64>,
    prob: Vec<f64>,
    morgue: i64,
    entry: i64,
}

static LOCAL_PARAMS: Mutex<LocalParams> = Mutex::new(LocalParams {
    dimension: 0,
    score: Vec::new(),
    prob: Vec::new(),
    morgue: 0,
    entry: 0,
});

fn n_setParameters(dimension_: usize, score_: &[i64], prob_: &[f64], entry_: i64) {
    let mut params = LOCAL_PARAMS.lock().unwrap();
    params.dimension = dimension_;
    params.score.clear();
    params.score.extend_from_slice(&score_[..dimension_]);
    params.prob.clear();
    params.prob.extend_from_slice(&prob_[..dimension_]);
    params.morgue = score_[0] - 1;
    params.entry = entry_;
}

fn n_totalProbAssoc(x_: f64) -> f64 {
    let params = LOCAL_PARAMS.lock().unwrap();
    let mut sum = 0.0;
    for i in 0..params.dimension {
        sum += params.prob[i] * (x_ * params.score[i] as f64).exp();
    }
    sum
}

fn n_meanPowerAssoc(x_: f64, power_: i64) -> f64 {
    let params = LOCAL_PARAMS.lock().unwrap();
    let mut sum = 0.0;
    for i in 0..params.dimension {
        sum += alp_integer::integerPower(params.score[i] as f64, power_)
            * params.prob[i]
            * (x_ * params.score[i] as f64).exp();
    }
    sum
}

fn n_meanAssoc(x_: f64) -> f64 {
    n_meanPowerAssoc(x_, 1)
}

fn n_bracket() -> (f64, f64) {
    const FACTOR: f64 = 0.5;
    let params = LOCAL_PARAMS.lock().unwrap();
    let mut p = -params.prob[params.dimension - 1].ln() / params.score[params.dimension - 1] as f64;
    drop(params);
    while 1.0 <= n_totalProbAssoc(p) {
        p *= FACTOR;
    }
    let q = p / FACTOR;
    (p, q)
}

pub fn n_step(oldValue_: i64, state_: usize) -> i64 {
    let params = LOCAL_PARAMS.lock().unwrap();
    assert!(state_ < params.dimension);
    if params.morgue < oldValue_ {
        oldValue_ + params.score[state_]
    } else {
        oldValue_
    }
}

pub fn n_bury(oldValue_: i64, state_: usize) -> i64 {
    let params = LOCAL_PARAMS.lock().unwrap();
    assert!(state_ < params.dimension);
    if params.entry < oldValue_ {
        oldValue_
    } else {
        params.morgue
    }
}

pub fn flatten(
    dimension_: usize,
    scoreMatrix_: &[Vec<i64>],
    prob_: &[Vec<f64>],
    dimension2_: usize,
) -> (usize, Vec<i64>, Vec<f64>) {
    let dimension2 = if dimension2_ == 0 {
        dimension_
    } else {
        dimension2_
    };

    let mut sum = 0.0;
    for i in 0..dimension_ {
        for j in 0..dimension2 {
            sum += prob_[i][j];
        }
    }
    const FUDGE: f64 = 20.0;
    assert!(alp_approx::relApprox(sum, 1.0, FUDGE * REL_TOL));

    let mut min = i64::MAX;
    let mut max = i64::MIN;
    for row in scoreMatrix_.iter().take(dimension_) {
        for &score in row.iter().take(dimension2) {
            min = min.min(score);
            max = max.max(score);
        }
    }
    assert!(min <= max);

    let dim = (max - min + 1) as usize;
    let mut p = vec![0.0; dim];
    for i in 0..dimension_ {
        for j in 0..dimension2 {
            p[(scoreMatrix_[i][j] - min) as usize] += prob_[i][j];
        }
    }

    let mut score = Vec::new();
    let mut prob = Vec::new();
    for s in min..=max {
        let ps = p[(s - min) as usize];
        if 0.0 < ps {
            score.push(s);
            prob.push(ps);
        }
    }

    (score.len(), score, prob)
}

pub fn lambda_matrix(dimMatrix_: usize, scoreMatrix_: &[Vec<i64>], q_: &[f64]) -> f64 {
    let mut prob = vec![vec![0.0; dimMatrix_]; dimMatrix_];
    for i in 0..dimMatrix_ {
        for j in 0..dimMatrix_ {
            prob[i][j] = q_[i] * q_[j];
        }
    }
    let (dim, score, p) = flatten(dimMatrix_, scoreMatrix_, &prob, 0);
    lambda(dim, &score, &p)
}

pub fn mu(dimension_: usize, score_: &[i64], prob_: &[f64]) -> f64 {
    let mut mu = 0.0;
    for i in 0..dimension_ {
        mu += score_[i] as f64 * prob_[i];
    }
    mu
}

pub fn lambda(dimension_: usize, score_: &[i64], prob_: &[f64]) -> f64 {
    n_setParameters(dimension_, score_, prob_, 0);
    let (p, q) = n_bracket();
    alp_root::bisectionNoParam(
        1.0,
        n_totalProbAssoc,
        p,
        q,
        REL_TOL * (p - q).abs(),
        0.0,
        None,
    )
}

pub fn muPowerAssoc(
    dimension_: usize,
    score_: &[i64],
    prob_: &[f64],
    mut lambda_: f64,
    power_: i64,
) -> f64 {
    n_setParameters(dimension_, score_, prob_, 0);
    if lambda_ == 0.0 {
        lambda_ = lambda(dimension_, score_, prob_);
    }
    n_meanPowerAssoc(lambda_, power_)
}

pub fn muAssoc(dimension_: usize, score_: &[i64], prob_: &[f64], lambda_: f64) -> f64 {
    muPowerAssoc(dimension_, score_, prob_, lambda_, 1)
}

pub fn thetaMin(dimension_: usize, score_: &[i64], prob_: &[f64], mut lambda_: f64) -> f64 {
    n_setParameters(dimension_, score_, prob_, 0);
    if lambda_ == 0.0 {
        lambda_ = lambda(dimension_, score_, prob_);
    }
    let (p, q) = n_bracket();
    alp_root::bisectionNoParam(
        0.0,
        n_meanAssoc,
        0.0,
        lambda_,
        REL_TOL * (p - q).abs(),
        0.0,
        None,
    )
}

pub fn rMin(
    dimension_: usize,
    score_: &[i64],
    prob_: &[f64],
    lambda_: f64,
    mut thetaMin_: f64,
) -> f64 {
    n_setParameters(dimension_, score_, prob_, 0);
    if thetaMin_ == 0.0 {
        thetaMin_ = thetaMin(dimension_, score_, prob_, lambda_);
    }
    n_totalProbAssoc(thetaMin_)
}

pub fn r(dimension_: usize, score_: &[i64], prob_: &[f64], theta_: f64) -> f64 {
    let mut sum = 0.0;
    for i in 0..dimension_ {
        sum += prob_[i] * (theta_ * score_[i] as f64).exp();
    }
    sum
}

pub fn delta(dimension_: usize, score_: &[i64]) -> i64 {
    let mut delta = 0;
    for &score in score_.iter().take(dimension_) {
        delta = alp_integer::euclidAlgorithm(delta, score);
    }
    delta
}

pub fn thetaMinusDelta(lambda_: f64, dimension_: usize, score_: &[i64]) -> f64 {
    let del = delta(dimension_, score_) as f64;
    (1.0 - (-lambda_ * del).exp()) / del
}

#[allow(clippy::too_many_arguments)]
pub fn descendingLadderEpochRepeat(
    dimension_: usize,
    score_: &[i64],
    prob_: &[f64],
    eSumAlpha_: Option<&mut f64>,
    eOneMinusExpSumAlpha_: Option<&mut f64>,
    isStrict_: bool,
    mut lambda_: f64,
    endW_: usize,
    mut pAlphaW_: Option<&mut [f64]>,
    mut eOneMinusExpSumAlphaW_: Option<&mut [f64]>,
    lambda0_: f64,
    mu0_: f64,
    muAssoc0_: f64,
    thetaMin0_: f64,
    rMin0_: f64,
    time_: f64,
    terminated_: Option<&mut bool>,
) {
    let mu0 = if mu0_ == 0.0 {
        mu(dimension_, score_, prob_)
    } else {
        mu0_
    };
    assert!(mu0 < 0.0);
    let lambda0 = if lambda0_ == 0.0 {
        lambda(dimension_, score_, prob_)
    } else {
        lambda0_
    };
    assert!(0.0 < lambda0);
    if lambda_ == 0.0 {
        lambda_ = lambda0;
    }
    assert!(0.0 < lambda_);
    let muAssoc0 = if muAssoc0_ == 0.0 {
        muAssoc(dimension_, score_, prob_, lambda0)
    } else {
        muAssoc0_
    };
    assert!(0.0 < muAssoc0);
    let thetaMin0 = if thetaMin0_ == 0.0 {
        thetaMin(dimension_, score_, prob_, lambda0)
    } else {
        thetaMin0_
    };
    assert!(0.0 < thetaMin0);
    let rMin0 = if rMin0_ == 0.0 {
        rMin(dimension_, score_, prob_, lambda0, thetaMin0)
    } else {
        rMin0_
    };
    assert!(0.0 < rMin0 && rMin0 < 1.0);

    let ITER_MIN = ((REL_TOL * (1.0 - rMin0)).ln() / rMin0.ln()) as i64;
    assert!(0 < ITER_MIN);
    let ITER = if (endW_ as i64) < ITER_MIN {
        ITER_MIN
    } else {
        endW_ as i64
    };
    assert!(0 < ITER);
    let Y_MAX = (-REL_TOL.ln() / lambda0) as i64;

    let entry = if isStrict_ { -1 } else { 0 };
    n_setParameters(dimension_, score_, prob_, entry);

    let mut time0 = 0.0;
    let mut time1 = 0.0;
    if time_ > 0.0 {
        sls_basic::get_current_time(&mut time0);
    }

    let mut dynProgProb = DynProgProbLim::new(
        Some(n_step),
        dimension_,
        Some(prob_),
        score_[0] - 1,
        Y_MAX,
        None,
    );

    if let Some(pAlphaW) = pAlphaW_.as_deref_mut() {
        pAlphaW[0] = 0.0;
    }
    if let Some(eOneMinusExpSumAlphaW) = eOneMinusExpSumAlphaW_.as_deref_mut() {
        eOneMinusExpSumAlphaW[0] = 0.0;
    }

    dynProgProb.update();

    let mut e_sum_alpha = 0.0;
    let mut e_one_minus_exp_sum_alpha = 0.0;

    for w in 1..ITER as usize {
        if w < endW_ {
            if let Some(pAlphaW) = pAlphaW_.as_deref_mut() {
                pAlphaW[w] = 0.0;
            }
            if let Some(eOneMinusExpSumAlphaW) = eOneMinusExpSumAlphaW_.as_deref_mut() {
                eOneMinusExpSumAlphaW[w] = 0.0;
            }
            for value in score_[0]..=entry {
                if let Some(pAlphaW) = pAlphaW_.as_deref_mut() {
                    pAlphaW[w] += dynProgProb.getProb(value);
                }
                if let Some(eOneMinusExpSumAlphaW) = eOneMinusExpSumAlphaW_.as_deref_mut() {
                    eOneMinusExpSumAlphaW[w] +=
                        dynProgProb.getProb(value) * (1.0 - (lambda_ * value as f64).exp());
                }
            }
        }

        for value in score_[0]..=entry {
            e_sum_alpha += dynProgProb.getProb(value) * value as f64;
            e_one_minus_exp_sum_alpha +=
                dynProgProb.getProb(value) * (1.0 - (lambda_ * value as f64).exp());
        }

        dynProgProb.setValueFct(Some(n_bury));
        dynProgProb.update();
        dynProgProb.setValueFct(Some(n_step));
        dynProgProb.update();

        if time_ > 0.0 {
            sls_basic::get_current_time(&mut time1);
            if time1 - time0 > time_ {
                if let Some(terminated) = terminated_ {
                    *terminated = true;
                }
                return;
            }
        }
    }

    for value in score_[0]..=entry {
        e_sum_alpha += dynProgProb.getProb(value) * value as f64;
        e_one_minus_exp_sum_alpha +=
            dynProgProb.getProb(value) * (1.0 - (lambda_ * value as f64).exp());
    }

    let mut prob = 0.0;
    for value in (entry + 1)..dynProgProb.getValueUpper() {
        prob += dynProgProb.getProb(value);
    }
    prob += dynProgProb.getProbLost();
    const FUDGE: f64 = 100.0;
    assert!(prob <= FUDGE * dimension_ as f64 * REL_TOL);

    if let Some(eSumAlpha) = eSumAlpha_ {
        *eSumAlpha = e_sum_alpha;
    }
    if let Some(eOneMinusExpSumAlpha) = eOneMinusExpSumAlpha_ {
        *eOneMinusExpSumAlpha = e_one_minus_exp_sum_alpha;
    }
}

#[allow(clippy::too_many_arguments)]
pub fn descendingLadderEpoch(
    dimension_: usize,
    score_: &[i64],
    prob_: &[f64],
    eSumAlpha_: Option<&mut f64>,
    eOneMinusExpSumAlpha_: Option<&mut f64>,
    isStrict_: bool,
    lambda0_: f64,
    mu0_: f64,
    muAssoc0_: f64,
    thetaMin0_: f64,
    rMin0_: f64,
    time_: f64,
    terminated_: Option<&mut bool>,
) {
    descendingLadderEpochRepeat(
        dimension_,
        score_,
        prob_,
        eSumAlpha_,
        eOneMinusExpSumAlpha_,
        isStrict_,
        0.0,
        0,
        None,
        None,
        lambda0_,
        mu0_,
        muAssoc0_,
        thetaMin0_,
        rMin0_,
        time_,
        terminated_,
    );
}

pub fn isProbDist(dimension_: usize, prob_: &[f64]) -> bool {
    let mut sum = 0.0;
    for &p in prob_.iter().take(dimension_) {
        if !(0.0..=1.0).contains(&p) {
            return false;
        }
        sum += p;
    }
    alp_approx::relApprox(sum, 1.0, REL_TOL)
}

pub fn isScoreIncreasing(dimension_: usize, score_: &[i64]) -> bool {
    for i in 1..dimension_ {
        if score_[i] <= score_[i - 1] {
            return false;
        }
    }
    true
}

pub fn isLogarithmic(dimension_: usize, score_: &[i64], prob_: &[f64]) -> bool {
    if !isScoreIncreasing(dimension_, score_) {
        return false;
    }
    if !isProbDist(dimension_, prob_) {
        return false;
    }
    if 0.0 <= mu(dimension_, score_, prob_) {
        return false;
    }
    if score_[dimension_ - 1] <= 0 {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    static TEST_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn test_flatten_and_distribution_checks() {
        let scores = vec![vec![-1, 2], vec![2, -1]];
        let probs = vec![vec![0.25, 0.25], vec![0.25, 0.25]];
        let (dim, score, p) = flatten(2, &scores, &probs, 0);
        assert_eq!(dim, 2);
        assert_eq!(score, vec![-1, 2]);
        assert_eq!(p, vec![0.5, 0.5]);
        assert!(isProbDist(2, &[0.5, 0.5]));
        assert!(!isProbDist(2, &[0.5, 0.6]));
        assert!(isScoreIncreasing(2, &[-1, 2]));
        assert!(!isScoreIncreasing(2, &[2, -1]));
    }

    #[test]
    fn test_random_walk_basic_parameters() {
        let _guard = TEST_LOCK.lock().unwrap();
        let score = [-1, 2];
        let prob = [0.8, 0.2];
        assert!(isLogarithmic(2, &score, &prob));
        assert!((mu(2, &score, &prob) + 0.4).abs() < 1e-12);
        let lambda = lambda(2, &score, &prob);
        assert!((lambda - 0.44568071901268186).abs() < 1e-6);
        assert!((r(2, &score, &prob, lambda) - 1.0).abs() < 1e-6);
        assert!(muAssoc(2, &score, &prob, lambda) > 0.0);
        assert_eq!(delta(2, &score), 1);
        assert!((thetaMinusDelta(lambda, 2, &score) - (1.0 - (-lambda).exp())).abs() < 1e-12);
    }

    #[test]
    fn test_descending_ladder_epoch_repeat_small_distribution() {
        let _guard = TEST_LOCK.lock().unwrap();
        let score = [-1, 2];
        let prob = [0.8, 0.2];
        let lambda = lambda(2, &score, &prob);
        let mu0 = mu(2, &score, &prob);
        let mu_assoc = muAssoc(2, &score, &prob, lambda);
        let theta_min = thetaMin(2, &score, &prob, lambda);
        let r_min = rMin(2, &score, &prob, lambda, theta_min);
        let mut e_sum = 0.0;
        let mut e_one_minus = 0.0;
        let mut p_alpha = vec![0.0; 4];
        let mut e_by_w = vec![0.0; 4];
        let mut terminated = false;
        descendingLadderEpochRepeat(
            2,
            &score,
            &prob,
            Some(&mut e_sum),
            Some(&mut e_one_minus),
            false,
            lambda,
            4,
            Some(&mut p_alpha),
            Some(&mut e_by_w),
            lambda,
            mu0,
            mu_assoc,
            theta_min,
            r_min,
            0.0,
            Some(&mut terminated),
        );
        assert!(!terminated);
        assert!(e_sum < 0.0);
        assert!(e_one_minus > 0.0);
        assert_eq!(p_alpha[0], 0.0);
        assert!(p_alpha[1] > 0.0);
        assert!(e_by_w[1] > 0.0);
    }
}
