#![allow(non_snake_case)]

use std::fmt;

use crate::stats::alp_localmaxstat::LocalMaxStat;
use crate::stats::alp_localmaxstat_util;

#[derive(Debug, Clone, PartialEq)]
pub struct LocalMaxStatMatrix {
    base: LocalMaxStat,
    d_dimMatrix: usize,
    d_scoreMatrix_p: Vec<Vec<i64>>,
    d_p_p: Vec<f64>,
    d_p2_p: Vec<f64>,
    d_dimMatrix2: usize,
}

impl Default for LocalMaxStatMatrix {
    fn default() -> Self {
        Self {
            base: LocalMaxStat::default(),
            d_dimMatrix: 0,
            d_scoreMatrix_p: Vec::new(),
            d_p_p: Vec::new(),
            d_p2_p: Vec::new(),
            d_dimMatrix2: 0,
        }
    }
}

impl LocalMaxStatMatrix {
    /// Rust-spelled equivalent of the primary C++ `copy` overload.
    pub fn copy_matrix(
        &mut self,
        dim_matrix: usize,
        score_matrix: &[Vec<i64>],
        p: &[f64],
        p2: Option<&[f64]>,
        dim_matrix2: usize,
    ) {
        self.copy(dim_matrix, score_matrix, p, p2, dim_matrix2);
    }

    /// Rust-spelled equivalent of the C++ `copy(LocalMaxStat, ...)` overload.
    pub fn copy_with_base(
        &mut self,
        local_max_stat: LocalMaxStat,
        dim_matrix: usize,
        score_matrix: &[Vec<i64>],
        p: &[f64],
        p2: Option<&[f64]>,
        dim_matrix2: usize,
    ) {
        self.copy_base(local_max_stat, dim_matrix, score_matrix, p, p2, dim_matrix2);
    }

    /// Rust-spelled equivalent of the inline matrix copy overload.
    pub fn copy_from(&mut self, other: &Self) {
        self.copy_local(other);
    }

    pub fn new(
        dimMatrix_: usize,
        scoreMatrix_: Option<&[Vec<i64>]>,
        p_: Option<&[f64]>,
        p2_: Option<&[f64]>,
        dimMatrix2_: usize,
        time_: f64,
    ) -> Self {
        let mut this = Self::default();
        LocalMaxStat::setTime(time_);
        if let (Some(score_matrix), Some(p)) = (scoreMatrix_, p_) {
            this.copy(dimMatrix_, score_matrix, p, p2_, dimMatrix2_);
        } else {
            this.copy(0, &[], &[], None, 0);
        }
        this
    }

    pub fn assign(&mut self, localMaxStat_: &Self) -> &mut Self {
        if !std::ptr::eq(self, localMaxStat_) {
            self.copy_local(localMaxStat_);
        }
        self
    }

    pub fn copy(
        &mut self,
        dimMatrix_: usize,
        scoreMatrix_: &[Vec<i64>],
        p_: &[f64],
        p2_: Option<&[f64]>,
        dimMatrix2_: usize,
    ) {
        let p2 = p2_.unwrap_or(p_);
        let dimMatrix2 = if dimMatrix2_ == 0 {
            dimMatrix_
        } else {
            dimMatrix2_
        };

        self.free2();
        self.init(dimMatrix_, dimMatrix2);

        if self.getDimMatrix() == 0 {
            self.base.copy_full(
                0,
                &[],
                &[],
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                false,
            );
            return;
        }

        for (dst, src) in self.d_scoreMatrix_p.iter_mut().zip(scoreMatrix_.iter()) {
            dst.copy_from_slice(&src[..dimMatrix2]);
        }
        self.d_p_p.copy_from_slice(&p_[..dimMatrix_]);
        self.d_p2_p.copy_from_slice(&p2[..dimMatrix2]);

        let mut probMatrix = vec![vec![0.0; self.getDimMatrix2()]; self.getDimMatrix()];
        for i in 0..self.getDimMatrix() {
            for j in 0..self.getDimMatrix2() {
                probMatrix[i][j] = p_[i] * p2[j];
            }
        }
        let (dim, score, p) = alp_localmaxstat_util::flatten(
            self.getDimMatrix(),
            self.getScoreMatrix(),
            &probMatrix,
            self.getDimMatrix2(),
        );
        self.base.copy(dim, &score, &p);
    }

    pub fn copy_base(
        &mut self,
        localMaxStat_: LocalMaxStat,
        dimMatrix_: usize,
        scoreMatrix_: &[Vec<i64>],
        p_: &[f64],
        p2_: Option<&[f64]>,
        dimMatrix2_: usize,
    ) {
        let p2 = p2_.unwrap_or(p_);
        let dimMatrix2 = if dimMatrix2_ == 0 {
            dimMatrix_
        } else {
            dimMatrix2_
        };
        self.free2();
        self.init(dimMatrix_, dimMatrix2);
        for (dst, src) in self.d_scoreMatrix_p.iter_mut().zip(scoreMatrix_.iter()) {
            dst.copy_from_slice(&src[..dimMatrix2]);
        }
        self.d_p_p.copy_from_slice(&p_[..dimMatrix_]);
        self.d_p2_p.copy_from_slice(&p2[..dimMatrix2]);
        self.base.copy_local(&localMaxStat_);
    }

    pub fn copy_local(&mut self, localMaxStatMatrix_: &Self) {
        self.copy_base(
            localMaxStatMatrix_.base.clone(),
            localMaxStatMatrix_.getDimMatrix(),
            localMaxStatMatrix_.getScoreMatrix(),
            localMaxStatMatrix_.getP(),
            Some(localMaxStatMatrix_.getP2()),
            localMaxStatMatrix_.getDimMatrix2(),
        );
    }

    pub fn bool_(&self) -> bool {
        self.base.bool_()
    }

    pub fn out(&self) -> String {
        self.base.out()
    }

    pub fn getR(&self, theta_: f64) -> f64 {
        self.base.getR(theta_)
    }

    pub fn getA(&self) -> f64 {
        self.base.getA()
    }

    pub fn getAlpha(&self) -> f64 {
        self.base.getAlpha()
    }

    pub fn getDimension(&self) -> usize {
        self.base.getDimension()
    }

    pub fn getScore(&self) -> &[i64] {
        self.base.getScore()
    }

    pub fn getProb(&self) -> &[f64] {
        self.base.getProb()
    }

    pub fn getLambda(&self) -> f64 {
        self.base.getLambda()
    }

    pub fn getK(&self) -> f64 {
        self.base.getK()
    }

    pub fn getC(&self) -> f64 {
        self.base.getC()
    }

    pub fn getTerminated(&self) -> bool {
        self.base.getTerminated()
    }

    pub fn getThetaMin(&self) -> f64 {
        self.base.getThetaMin()
    }

    pub fn getRMin(&self) -> f64 {
        self.base.getRMin()
    }

    pub fn getDelta(&self) -> i64 {
        self.base.getDelta()
    }

    pub fn getThetaMinusDelta(&self) -> f64 {
        self.base.getThetaMinusDelta()
    }

    pub fn getMu(&self) -> f64 {
        self.base.getMu()
    }

    pub fn getSigma(&self) -> f64 {
        self.base.getSigma()
    }

    pub fn getMuAssoc(&self) -> f64 {
        self.base.getMuAssoc()
    }

    pub fn getSigmaAssoc(&self) -> f64 {
        self.base.getSigmaAssoc()
    }

    pub fn getMeanWDLE(&self) -> f64 {
        self.base.getMeanWDLE()
    }

    pub fn getDimMatrix(&self) -> usize {
        self.d_dimMatrix
    }

    pub fn getScoreMatrix(&self) -> &[Vec<i64>] {
        &self.d_scoreMatrix_p
    }

    pub fn getP(&self) -> &[f64] {
        &self.d_p_p
    }

    pub fn getP2(&self) -> &[f64] {
        &self.d_p2_p
    }

    pub fn getDimMatrix2(&self) -> usize {
        self.d_dimMatrix2
    }

    pub fn is_ready(&self) -> bool {
        self.bool_()
    }

    pub fn get_r(&self, theta: f64) -> f64 {
        self.getR(theta)
    }
    pub fn get_a(&self) -> f64 {
        self.getA()
    }
    pub fn get_alpha(&self) -> f64 {
        self.getAlpha()
    }
    pub fn dimension(&self) -> usize {
        self.getDimension()
    }
    pub fn scores(&self) -> &[i64] {
        self.getScore()
    }
    pub fn probabilities(&self) -> &[f64] {
        self.getProb()
    }
    pub fn lambda(&self) -> f64 {
        self.getLambda()
    }
    pub fn k(&self) -> f64 {
        self.getK()
    }
    pub fn c(&self) -> f64 {
        self.getC()
    }
    pub fn terminated(&self) -> bool {
        self.getTerminated()
    }
    pub fn theta_min(&self) -> f64 {
        self.getThetaMin()
    }
    pub fn r_min(&self) -> f64 {
        self.getRMin()
    }
    pub fn delta(&self) -> i64 {
        self.getDelta()
    }
    pub fn theta_minus_delta(&self) -> f64 {
        self.getThetaMinusDelta()
    }
    pub fn mu(&self) -> f64 {
        self.getMu()
    }
    pub fn sigma(&self) -> f64 {
        self.getSigma()
    }
    pub fn mu_assoc(&self) -> f64 {
        self.getMuAssoc()
    }
    pub fn sigma_assoc(&self) -> f64 {
        self.getSigmaAssoc()
    }
    pub fn mean_wdle(&self) -> f64 {
        self.getMeanWDLE()
    }
    pub fn dim_matrix(&self) -> usize {
        self.getDimMatrix()
    }
    pub fn score_matrix(&self) -> &[Vec<i64>] {
        self.getScoreMatrix()
    }
    pub fn p(&self) -> &[f64] {
        self.getP()
    }
    pub fn p2(&self) -> &[f64] {
        self.getP2()
    }
    pub fn dim_matrix2(&self) -> usize {
        self.getDimMatrix2()
    }

    fn init(&mut self, dimMatrix_: usize, dimMatrix2_: usize) {
        let dimMatrix2 = if dimMatrix2_ == 0 {
            dimMatrix_
        } else {
            dimMatrix2_
        };
        self.d_scoreMatrix_p = vec![vec![0; dimMatrix2]; dimMatrix_];
        self.d_p_p = vec![0.0; dimMatrix_];
        self.d_p2_p = vec![0.0; dimMatrix2];
        self.d_dimMatrix = dimMatrix_;
        self.d_dimMatrix2 = dimMatrix2;
    }

    fn free2(&mut self) {
        self.d_scoreMatrix_p.clear();
        self.d_p_p.clear();
        self.d_p2_p.clear();
        self.d_dimMatrix = 0;
        self.d_dimMatrix2 = 0;
    }
}

impl fmt::Display for LocalMaxStatMatrix {
    fn fmt(&self, ostr_: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(ostr_, "{}", self.out())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    #[test]
    fn test_local_max_stat_matrix_empty_and_copy_base() {
        let matrix = LocalMaxStatMatrix::new(0, None, None, None, 0, 0.0);
        assert!(!matrix.bool_());
        assert_eq!(matrix.out(), "");
        assert_eq!(format!("{}", matrix), "");

        let base = LocalMaxStat::new(0, None, None);
        let scores = vec![vec![-1, 2]];
        let mut copied = LocalMaxStatMatrix::default();
        copied.copy_base(base, 1, &scores, &[1.0], None, 1);
        assert_eq!(copied.getDimMatrix(), 1);
        assert_eq!(copied.getDimMatrix2(), 1);
        assert_eq!(copied.getScoreMatrix(), &[vec![-1]]);
        assert_eq!(copied.getP(), &[1.0]);
        assert_eq!(copied.getP2(), &[1.0]);
    }

    #[test]
    fn test_local_max_stat_matrix_computes_flattened_distribution() {
        let _guard = TEST_LOCK.lock().unwrap();
        let scores = vec![vec![-1, 2], vec![-1, -1]];
        let p = [0.5, 0.5];
        let stat = LocalMaxStatMatrix::new(2, Some(&scores), Some(&p), None, 0, 0.0);
        assert!(stat.bool_());
        assert_eq!(stat.getDimMatrix(), 2);
        assert_eq!(stat.getDimMatrix2(), 2);
        assert_eq!(stat.getScore(), &[-1, 2]);
        assert_eq!(stat.getProb(), &[0.75, 0.25]);
        assert!(stat.getLambda() > 0.0);
    }
}
