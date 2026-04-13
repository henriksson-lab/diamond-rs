use crate::basic::value::{Letter, AMINO_ACID_COUNT};
use super::standard_matrix::StandardMatrix;
use super::matrices;

use std::f64::consts::LN_2;

/// Runtime score matrix used for alignment.
///
/// Holds the scoring data in various representations optimized for
/// different SIMD widths (8-bit, 16-bit, 32-bit), plus statistical
/// parameters for E-value computation.
pub struct ScoreMatrix {
    /// 32x32 score matrix (padded for alignment), 32-bit scores.
    matrix32: [i32; 32 * 32],
    /// 32x32 score matrix, 8-bit scores.
    matrix8: [i8; 32 * 32],
    /// 32x32 unsigned 8-bit scores (biased).
    matrix8u: [u8; 32 * 32],
    /// Bias added to make all scores non-negative (for unsigned representation).
    bias: i8,
    /// Gap open penalty.
    gap_open: i32,
    /// Gap extend penalty.
    gap_extend: i32,
    /// Frame shift penalty (for blastx).
    frame_shift: i32,
    /// Database size in letters (for E-value calculation).
    db_letters: f64,
    /// ln(K) for bit score calculation.
    ln_k: f64,
    /// Scale factor (used for scaled matrix variants).
    #[allow(dead_code)]
    scale: f64,
    /// Lambda parameter.
    lambda: f64,
    /// K parameter.
    k: f64,
    /// Matrix name.
    name: String,
    /// Reference to the underlying standard matrix (used for joint_probs, freq_ratios).
    #[allow(dead_code)]
    standard_matrix: &'static StandardMatrix,
    /// ALP area parameters for E-value computation.
    area_params: super::pvalues::AreaParams,
}

impl ScoreMatrix {
    /// Create a score matrix from a named matrix and gap penalties.
    pub fn new(
        matrix_name: &str,
        gap_open: i32,
        gap_extend: i32,
        frame_shift: i32,
        stop_match_score: i32,
        db_letters: u64,
    ) -> Result<Self, String> {
        let standard_matrix = matrices::get_matrix(matrix_name)
            .ok_or_else(|| format!("Unknown scoring matrix: {}", matrix_name))?;

        let go = if gap_open < 0 {
            standard_matrix.default_gap_open
        } else {
            gap_open
        };
        let ge = if gap_extend < 0 {
            standard_matrix.default_gap_extend
        } else {
            gap_extend
        };

        let params = standard_matrix.constants(go, ge);

        // Build the 32x32 matrices from the AMINO_ACID_COUNT x AMINO_ACID_COUNT scores
        let n = AMINO_ACID_COUNT as i32;
        let mut matrix32 = [i32::MIN; 32 * 32];
        let mut matrix8 = [i8::MIN; 32 * 32];

        for i in 0..32 {
            for j in 0..32 {
                if i < n as usize && j < n as usize {
                    let score = standard_matrix.scores[i * AMINO_ACID_COUNT + j] as i32;
                    matrix32[i * 32 + j] = score;
                    matrix8[i * 32 + j] = score as i8;
                }
            }
        }

        // Override stop-stop score if specified
        if stop_match_score != 1 {
            matrix32[24 * 32 + 24] = stop_match_score;
            matrix8[24 * 32 + 24] = stop_match_score as i8;
        }

        // Compute bias (minimum score, negated)
        let min_score = matrix8.iter().filter(|&&s| s != i8::MIN).copied().min().unwrap_or(0);
        let bias = if min_score < 0 { -min_score } else { 0 };

        // Build unsigned biased matrix
        let mut matrix8u = [0u8; 32 * 32];
        for i in 0..32 {
            for j in 0..32 {
                if i < n as usize && j < n as usize {
                    matrix8u[i * 32 + j] = (matrix8[i * 32 + j] as i16 + bias as i16) as u8;
                }
            }
        }

        let ln_k = params.k.ln();

        // Compute ALP area parameters from actual matrix constants.
        // Matches C++ score_matrix.cpp line 48-49:
        //   b = 2*G*(u.alpha - p.alpha)
        //   beta = 2*G*(u.alpha_v - p.alpha_v)
        //   tau = 2*G*(u.alpha_v - p.sigma)
        //   {lambda, K, p.alpha, b, p.alpha, b, p.alpha_v, beta, p.alpha_v, beta, p.sigma, tau}
        let ungapped = standard_matrix.ungapped_constants();
        let g = (go + ge) as f64;
        let a_val = params.alpha;
        let b_val = 2.0 * g * (ungapped.alpha - params.alpha);
        let beta_val = 2.0 * g * (ungapped.alpha_v - params.alpha_v);
        let tau_val = 2.0 * g * (ungapped.alpha_v - params.sigma);
        let area_params = super::pvalues::AreaParams {
            a_i: a_val, b_i: b_val, alpha_i: params.alpha_v, beta_i: beta_val,
            a_j: a_val, b_j: b_val, alpha_j: params.alpha_v, beta_j: beta_val,
            sigma: params.sigma, tau: tau_val,
        };

        Ok(ScoreMatrix {
            matrix32,
            matrix8,
            matrix8u,
            bias,
            gap_open: go,
            gap_extend: ge,
            frame_shift,
            db_letters: db_letters as f64,
            ln_k,
            scale: 1.0,
            lambda: params.lambda,
            k: params.k,
            name: matrix_name.to_lowercase(),
            standard_matrix,
            area_params,
        })
    }

    /// Score for aligning two letters.
    #[inline]
    pub fn score(&self, a: Letter, b: Letter) -> i32 {
        self.matrix32[(a as usize) * 32 + (b as usize)]
    }

    /// Get a row of the 32-bit score matrix.
    #[inline]
    pub fn row(&self, a: Letter) -> &[i32] {
        let start = (a as usize) * 32;
        &self.matrix32[start..start + 32]
    }

    /// Unsigned biased score.
    #[inline]
    pub fn biased_score(&self, a: Letter, b: Letter) -> u8 {
        self.matrix8u[(a as usize) * 32 + (b as usize)]
    }

    /// Convert raw alignment score to bit score (simple, no length correction).
    pub fn bitscore(&self, raw_score: f64) -> f64 {
        (self.lambda * raw_score - self.ln_k) / LN_2
    }

    /// Convert raw alignment score to bit score with ALP area correction.
    ///
    /// Matches C++ ScoreMatrix::bitscore_corrected():
    ///   (lambda * raw_score - ln(K) - ln(area)) / ln(2)
    pub fn bitscore_corrected(&self, raw_score: i32, query_len: u32, subject_len: u32) -> f64 {
        let area = super::pvalues::compute_area(
            &self.area_params, raw_score as f64, query_len as f64, subject_len as f64,
        );
        let log_area = if area > 0.0 { area.ln() } else { 0.0 };
        (self.lambda * raw_score as f64 - self.ln_k - log_area) / LN_2
    }

    /// Convert bit score to raw score.
    pub fn rawscore(&self, bitscore: f64) -> f64 {
        (bitscore * LN_2 + self.ln_k) / self.lambda
    }

    /// Compute the minimum ungapped raw score for a given query length and e-value threshold.
    ///
    /// Matches C++ CutoffTable: finds the minimum raw score S such that
    /// K * query_len * 1e9 * exp(-lambda * S) <= evalue_threshold.
    /// The factor 1e9 / subject_len normalizes per-residue (C++ evalue_norm), with
    /// subject_len approximated as 1e9 for the table.
    pub fn ungapped_cutoff(&self, query_len: usize, evalue_threshold: f64) -> i32 {
        if evalue_threshold <= 0.0 {
            return 0;
        }
        // E = K * m * n * exp(-lambda * S), with n normalized to 1e9
        // S = -(ln(E / (K * m * 1e9))) / lambda
        let m = query_len as f64;
        let inner = evalue_threshold / (self.k * m * 1e9);
        if inner <= 0.0 {
            return 0;
        }
        let s = -(inner.ln()) / self.lambda;
        s.ceil() as i32
    }

    /// Compute the minimum ungapped raw score using actual database size.
    ///
    /// Uses db_letters for normalization instead of 1e9, so it works correctly
    /// for both small and large databases.
    pub fn ungapped_cutoff_db(&self, query_len: usize, evalue_threshold: f64) -> i32 {
        if evalue_threshold <= 0.0 || self.db_letters <= 0.0 {
            return 0;
        }
        let m = query_len as f64;
        let inner = evalue_threshold / (self.k * m * self.db_letters);
        if inner >= 1.0 {
            return 0; // threshold is so permissive that any positive score passes
        }
        if inner <= 0.0 {
            return 0;
        }
        let s = -(inner.ln()) / self.lambda;
        s.ceil().max(0.0) as i32
    }

    /// Compute E-value from raw score and lengths.
    ///
    /// Uses the ALP pvalues finite-size correction with actual matrix parameters.
    /// Matches C++ ScoreMatrix::evalue():
    ///   evaluer.evalue(raw/scale, qlen, slen) * db_letters / slen
    pub fn evalue(&self, raw_score: i32, query_len: u32, subject_len: u32) -> f64 {
        let pairwise = super::pvalues::evalue_with_area(
            self.lambda, self.k, &self.area_params,
            raw_score as f64, query_len as f64, subject_len as f64,
        );
        // Scale to database level (matches C++ db_letters_ / subject_len)
        if self.db_letters > 0.0 && subject_len > 0 {
            pairwise * self.db_letters / subject_len as f64
        } else {
            pairwise
        }
    }

    pub fn lambda(&self) -> f64 {
        self.lambda
    }

    pub fn k(&self) -> f64 {
        self.k
    }

    pub fn ln_k(&self) -> f64 {
        self.ln_k
    }

    pub fn bias(&self) -> i8 {
        self.bias
    }

    pub fn gap_open(&self) -> i32 {
        self.gap_open
    }

    pub fn gap_extend(&self) -> i32 {
        self.gap_extend
    }

    pub fn frame_shift(&self) -> i32 {
        self.frame_shift
    }

    pub fn db_letters(&self) -> u64 {
        self.db_letters as u64
    }

    pub fn set_db_letters(&mut self, n: u64) {
        self.db_letters = n as f64;
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// Get the raw 8-bit score matrix data (for SIMD kernels).
    pub fn matrix8(&self) -> &[i8; 32 * 32] {
        &self.matrix8
    }

    /// Get the raw 32-bit score matrix data.
    pub fn matrix32(&self) -> &[i32; 32 * 32] {
        &self.matrix32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_score_matrix_creation() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        assert_eq!(sm.gap_open(), 11);
        assert_eq!(sm.gap_extend(), 1);
        assert_eq!(sm.name(), "blosum62");
    }

    #[test]
    fn test_score_matrix_scores() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        // A-A = 4
        assert_eq!(sm.score(0, 0), 4);
        // A-R = -1
        assert_eq!(sm.score(0, 1), -1);
        // Symmetric
        assert_eq!(sm.score(0, 1), sm.score(1, 0));
    }

    #[test]
    fn test_bitscore() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let bs = sm.bitscore(50.0);
        assert!(bs > 0.0);
        // Round-trip
        let rs = sm.rawscore(bs);
        assert!((rs - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_evalue() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let e = sm.evalue(50, 100, 1000);
        assert!(e >= 0.0); // With ALP finite-size correction, may be 0.0 for high scores
        assert!(e < 1.0);
        // Lower scores on larger sequences should produce non-zero E-values
        let e2 = sm.evalue(20, 500, 10000);
        assert!(e2 > 0.0, "Low score on large sequences should give non-zero E-value");
    }

    #[test]
    fn test_default_gap_penalties() {
        let sm = ScoreMatrix::new("blosum62", -1, -1, 0, 1, 0).unwrap();
        assert_eq!(sm.gap_open(), 11);
        assert_eq!(sm.gap_extend(), 1);
    }

    #[test]
    fn test_unknown_matrix() {
        assert!(ScoreMatrix::new("nonexistent", 11, 1, 0, 1, 0).is_err());
    }
}
