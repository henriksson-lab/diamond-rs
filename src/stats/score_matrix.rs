use super::matrices;
use super::standard_matrix::StandardMatrix;
use crate::basic::value::{Letter, AMINO_ACID_COUNT, TRUE_AA};

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
    /// Expected score for each true amino acid against BLOSUM62 background.
    background_scores: [f64; TRUE_AA as usize],
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

        let params = standard_matrix.constants(go, ge)?;

        // Build the 32x32 matrices from the AMINO_ACID_COUNT x AMINO_ACID_COUNT scores.
        // C++ `Scores<T>::Scores` (`diamond/src/stats/score_matrix.h:38-46`) fills
        // out-of-range cells (i >= n || j >= n) with `SCHAR_MIN` cast to T.
        // For T=int8_t that's `i8::MIN = -128` (same as Rust);
        // for T=int32_t that's `-128` (NOT `i32::MIN`);
        // for T=uint8_t that's `(uint8_t)(-128) = 128`.
        // The earlier Rust used `i32::MIN` for matrix32 and `0` for matrix8u,
        // both off vs C++. SIMD profile loads that index by a delimiter / pad
        // letter (e.g. DELIMITER_LETTER=31 in a padded SWIPE lane) saw scores
        // that were ~2 billion / 128 units off and silently corrupted DP.
        let n = AMINO_ACID_COUNT as i32;
        let mut matrix32 = [-128i32; 32 * 32];
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
        let min_score = matrix8
            .iter()
            .filter(|&&s| s != i8::MIN)
            .copied()
            .min()
            .unwrap_or(0);
        let bias = if min_score < 0 { -min_score } else { 0 };

        // Build unsigned biased matrix. Out-of-range cells must match C++:
        // `(uint8_t)SCHAR_MIN = 128` — NOT the biased value, NOT 0.
        let mut matrix8u = [128u8; 32 * 32];
        for i in 0..32 {
            for j in 0..32 {
                if i < n as usize && j < n as usize {
                    matrix8u[i * 32 + j] = (matrix8[i * 32 + j] as i16 + bias as i16) as u8;
                }
            }
        }

        let ln_k = params.k.ln();
        let mut background_scores = [0.0f64; TRUE_AA as usize];
        let bg_freq: [f64; TRUE_AA as usize] = [
            7.4216205067993410e-02,
            5.1614486141284638e-02,
            4.4645808512757915e-02,
            5.3626000838554413e-02,
            2.4687457167944848e-02,
            3.4259650591416023e-02,
            5.4311925684587502e-02,
            7.4146941452644999e-02,
            2.6212984805266227e-02,
            6.7917367618953756e-02,
            9.8907868497150955e-02,
            5.8155682303079680e-02,
            2.4990197579643110e-02,
            4.7418459742284751e-02,
            3.8538003320306206e-02,
            5.7229029476494421e-02,
            5.0891364550287033e-02,
            1.3029956129972148e-02,
            3.2281512313758580e-02,
            7.2919098205619245e-02,
        ];
        for i in 0..TRUE_AA as usize {
            for j in 0..TRUE_AA as usize {
                background_scores[i] += bg_freq[j] * matrix32[i * 32 + j] as f64;
            }
        }

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
        // Thresholds from sls_pvalues.cpp:353-355 (compute_tmp_values).
        // Both alpha_i and alpha_j map to params.alpha_v.
        let nat = super::pvalues::NAT_CUT_OFF_IN_MAX;
        let vi_y_thr = (nat * params.alpha_v / params.lambda).max(0.0);
        let vj_y_thr = vi_y_thr;
        let c_y_thr = (nat * params.sigma / params.lambda).max(0.0);
        let area_params = super::pvalues::AreaParams {
            a_i: a_val,
            b_i: b_val,
            alpha_i: params.alpha_v,
            beta_i: beta_val,
            a_j: a_val,
            b_j: b_val,
            alpha_j: params.alpha_v,
            beta_j: beta_val,
            sigma: params.sigma,
            tau: tau_val,
            vi_y_thr,
            vj_y_thr,
            c_y_thr,
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
            background_scores,
            area_params,
        })
    }

    /// Score for aligning two letters.
    ///
    /// Strip `SEED_MASK` (bit 7) before indexing: a Letter with the high
    /// bit set (i.e. tantan soft-masked) is negative as `i8`, and `as usize`
    /// sign-extends through `i64`, indexing far past the 32×32 = 1024-entry
    /// matrix. C++ `Sequence::operator[]` (`diamond/src/basic/sequence.h:87-94`)
    /// applies the same mask whenever `SEQ_MASK` is defined (the default
    /// build flag), so callers there never see the high bit reach this
    /// lookup; mirror that behavior here.
    #[inline]
    pub fn score(&self, a: Letter, b: Letter) -> i32 {
        let ai = (a & crate::basic::value::LETTER_MASK) as usize;
        let bi = (b & crate::basic::value::LETTER_MASK) as usize;
        self.matrix32[ai * 32 + bi]
    }

    /// Get a row of the 32-bit score matrix.
    #[inline]
    pub fn row(&self, a: Letter) -> &[i32] {
        let start = ((a & crate::basic::value::LETTER_MASK) as usize) * 32;
        &self.matrix32[start..start + 32]
    }

    /// Unsigned biased score.
    #[inline]
    pub fn biased_score(&self, a: Letter, b: Letter) -> u8 {
        let ai = (a & crate::basic::value::LETTER_MASK) as usize;
        let bi = (b & crate::basic::value::LETTER_MASK) as usize;
        self.matrix8u[ai * 32 + bi]
    }

    /// Convert raw alignment score to bit score (simple, no length correction).
    pub fn bitscore(&self, raw_score: f64) -> f64 {
        let s = (raw_score / self.scale).round();
        (self.lambda * s - self.ln_k) / LN_2
    }

    /// Convert raw alignment score to bit score with ALP area correction.
    ///
    /// Matches C++ `ScoreMatrix::bitscore_corrected(int, unsigned, unsigned)`.
    ///   (lambda * raw_score - ln(K) - log_area) / ln(2)
    /// where `log_area` comes from `evaluer.log_area(...)` directly.
    /// Going through `compute_area(...).ln()` (= exp-then-ln) underflows for
    /// large query×subject products (area → inf → log = inf) and snaps to 0
    /// for tiny areas via the `if area > 0.0` guard — both modes silently
    /// lose multiple bits vs C++.
    pub fn bitscore_corrected(&self, raw_score: i32, query_len: u32, subject_len: u32) -> f64 {
        let log_area = super::pvalues::compute_log_area(
            &self.area_params,
            raw_score as f64,
            query_len as f64,
            subject_len as f64,
        );
        (self.lambda * raw_score as f64 - self.ln_k - log_area) / LN_2
    }

    /// Convert bit score to raw score.
    pub fn rawscore(&self, bitscore: f64) -> f64 {
        (bitscore * LN_2 + self.ln_k) / self.lambda
    }

    pub fn rawscore_int(&self, bitscore: f64) -> i32 {
        self.rawscore(bitscore).ceil() as i32
    }

    /// Compute the minimum ungapped raw score for a given query length and e-value threshold.
    ///
    /// Matches C++ `Util::Scores::CutoffTable::operator()(unsigned)`.
    /// the table is precomputed at `m = 1 << (b - 1)` for each bit `b`, and
    /// `operator()` indexes by `b = 32 - clz(query_len)`. The effective `m`
    /// is therefore the **largest power of 2 ≤ query_len** — NOT the actual
    /// query length. For `query_len = 250`, C++ uses `m = 128`; using 250
    /// gives a stricter cutoff that rejects seeds C++ would accept.
    pub fn ungapped_cutoff(&self, query_len: usize, evalue_threshold: f64) -> i32 {
        if evalue_threshold <= 0.0 || query_len == 0 {
            return 0;
        }
        // Snap to floor power of 2 (largest power of 2 ≤ query_len).
        // `bits = 32 - leading_zeros(query_len as u32)` matches C++'s
        // `b = 32 - clz(query_len)`; the table cell is precomputed at
        // `1 << (b - 1)`.
        let bits = 32 - (query_len as u32).leading_zeros();
        let m_binned = (1u32 << (bits - 1)) as f64;
        let inner = evalue_threshold / (self.k * m_binned * 1e9);
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
    /// Matches C++ `ScoreMatrix::evalue(int, unsigned, unsigned)`.
    ///   evaluer.evalue((double)raw_score / scale_, qlen, slen) * db_letters / slen
    pub fn evalue(&self, raw_score: i32, query_len: u32, subject_len: u32) -> f64 {
        // C++ divides raw_score by `scale_` before handing it to the ALP
        // evaluer (`score_matrix.cpp:222`). For `scale = 1.0` (default) this
        // is a no-op, but a scaled matrix (matrix adjustment / CBS rescaling)
        // would otherwise see inflated scores and produce evalues off by a
        // power of `exp(lambda * (scale-1) * raw)`.
        let pairwise = super::pvalues::evalue_with_area(
            self.lambda,
            self.k,
            &self.area_params,
            raw_score as f64 / self.scale,
            query_len as f64,
            subject_len as f64,
        );
        // Scale to database level (matches C++ db_letters_ / subject_len)
        if self.db_letters > 0.0 && subject_len > 0 {
            pairwise * self.db_letters / subject_len as f64
        } else {
            pairwise
        }
    }

    /// Normalized E-value using a database size of 1e9.
    ///
    /// Matches C++ `ScoreMatrix::evalue_norm(int, unsigned, unsigned)`.
    pub fn evalue_norm(&self, raw_score: i32, query_len: u32, subject_len: u32) -> f64 {
        let pairwise = super::pvalues::evalue_with_area(
            self.lambda,
            self.k,
            &self.area_params,
            raw_score as f64 / self.scale,
            query_len as f64,
            subject_len as f64,
        );
        if subject_len > 0 {
            pairwise * 1e9 / subject_len as f64
        } else {
            pairwise
        }
    }

    /// Normalized E-value approximation for a query length only.
    ///
    /// Matches C++ `ScoreMatrix::evalue_norm(int, int)`.
    pub fn evalue_norm_query_len(&self, raw_score: i32, query_len: i32) -> f64 {
        1e9 * query_len as f64 * 2.0f64.powf(-self.bitscore(raw_score as f64 * self.scale))
    }

    pub fn bitscore_norm(&self, evalue: f64, query_len: u32) -> f64 {
        -(evalue / 1e9 / query_len as f64).ln() / LN_2
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

    pub fn low_score(&self) -> i8 {
        let mut low = i8::MAX;
        for i in 0..AMINO_ACID_COUNT {
            for j in i + 1..AMINO_ACID_COUNT {
                low = low.min(self.score(i as Letter, j as Letter) as i8);
            }
        }
        low
    }

    pub fn high_score(&self) -> i8 {
        let mut high = i8::MIN;
        for i in 0..AMINO_ACID_COUNT {
            for j in i..AMINO_ACID_COUNT {
                high = high.max(self.score(i as Letter, j as Letter) as i8);
            }
        }
        high
    }

    pub fn avg_id_score(&self) -> f64 {
        let mut s = 0.0;
        for i in 0..TRUE_AA as usize {
            s += self.score(i as Letter, i as Letter) as f64;
        }
        s / TRUE_AA as f64
    }

    pub fn report_cutoff(
        &self,
        score: i32,
        evalue: f64,
        min_bit_score: f64,
        max_evalue: f64,
    ) -> bool {
        if min_bit_score != 0.0 {
            self.bitscore(score as f64) >= min_bit_score
        } else {
            evalue <= max_evalue
        }
    }

    pub fn ideal_lambda(&self) -> Option<f64> {
        super::cbs::ideal_lambda(self)
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

    pub fn background_scores(&self) -> &[f64; TRUE_AA as usize] {
        &self.background_scores
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CutoffTable {
    data: [i32; Self::MAX_BITS + 1],
}

impl CutoffTable {
    pub const MAX_BITS: usize = 31;

    pub fn new(score_matrix: &ScoreMatrix, evalue: f64) -> Self {
        let mut data = [0; Self::MAX_BITS + 1];
        for b in 1..=Self::MAX_BITS {
            data[b] = score_matrix.rawscore_int(score_matrix.bitscore_norm(evalue, 1 << (b - 1)));
        }
        Self { data }
    }

    pub fn get(&self, query_len: i32) -> i32 {
        let b = 32 - (query_len as u32).leading_zeros() as usize;
        self.data[b]
    }

    pub fn call(&self, query_len: i32) -> i32 {
        self.get(query_len)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CutoffTable2D {
    data: [[i32; Self::MAX_BITS + 1]; Self::MAX_BITS + 1],
}

impl CutoffTable2D {
    pub const MAX_BITS: usize = 31;

    pub fn new(score_matrix: &ScoreMatrix, evalue: f64) -> Self {
        let mut data = [[0; Self::MAX_BITS + 1]; Self::MAX_BITS + 1];
        for b1 in 1..=Self::MAX_BITS {
            for b2 in 1..=Self::MAX_BITS {
                data[b1][b2] =
                    Self::calc_min_score(score_matrix, 1 << (b1 - 1), 1 << (b2 - 1), evalue);
            }
        }
        Self { data }
    }

    pub fn get(&self, query_len: i32, target_len: i32) -> i32 {
        let b1 = 32 - (query_len as u32).leading_zeros() as usize;
        let b2 = 32 - (target_len as u32).leading_zeros() as usize;
        self.data[b1][b2]
    }

    pub fn call(&self, query_len: i32, target_len: i32) -> i32 {
        self.get(query_len, target_len)
    }

    fn calc_min_score(score_matrix: &ScoreMatrix, qlen: u32, slen: u32, evalue: f64) -> i32 {
        for i in 10..1000 {
            if score_matrix.evalue_norm(i, qlen, slen) <= evalue {
                return i;
            }
        }
        1000
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
        assert_eq!(sm.background_scores().len(), TRUE_AA as usize);
        assert!(sm.background_scores()[0] < 1.0);
        assert!(sm.ideal_lambda().unwrap() > 0.0);
        assert_eq!(sm.low_score(), -12);
        assert_eq!(sm.high_score(), 11);
        assert!((sm.avg_id_score() - 5.8).abs() < 1e-12);
        assert!(sm.report_cutoff(50, 1e-10, 0.0, 1e-5));
        assert!(!sm.report_cutoff(50, 1e-4, 0.0, 1e-5));
        assert!(sm.report_cutoff(50, 1e-4, 10.0, 1e-5));
    }

    #[test]
    fn test_bitscore() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let bs = sm.bitscore(50.0);
        assert!(bs > 0.0);
        assert_eq!(sm.bitscore(50.4), bs);
        // Round-trip
        let rs = sm.rawscore(bs);
        assert!((rs - 50.0).abs() < 0.01);
        assert_eq!(sm.rawscore_int(bs), 50);
    }

    #[test]
    fn test_default_ungapped_xdrop_bits_to_raw_score() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        assert_eq!(sm.rawscore_int(12.3), 20);
    }

    #[test]
    fn test_evalue() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let e = sm.evalue(50, 100, 1000);
        assert!(e >= 0.0); // With ALP finite-size correction, may be 0.0 for high scores
        assert!(e < 1.0);
        // Lower scores on larger sequences should produce non-zero E-values
        let e2 = sm.evalue(20, 500, 10000);
        assert!(
            e2 > 0.0,
            "Low score on large sequences should give non-zero E-value"
        );
        let norm = sm.evalue_norm(20, 500, 10000);
        assert!(norm > 0.0);
        let short_norm = sm.evalue_norm_query_len(20, 500);
        assert!(short_norm > 0.0);
        assert!((sm.bitscore_norm(1.0, 500) - 38.86313713864835).abs() < 1e-12);
    }

    #[test]
    fn test_cutoff_tables() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let table = CutoffTable::new(&sm, 1e-3);
        assert_eq!(table.get(1), sm.rawscore_int(sm.bitscore_norm(1e-3, 1)));
        assert_eq!(table.get(513), sm.rawscore_int(sm.bitscore_norm(1e-3, 512)));

        let table2 = CutoffTable2D::new(&sm, 1e-3);
        let cutoff = table2.get(200, 700);
        assert!(cutoff >= 10);
        assert!(sm.evalue_norm(cutoff, 128, 512) <= 1e-3);
        if cutoff > 10 {
            assert!(sm.evalue_norm(cutoff - 1, 128, 512) > 1e-3);
        }
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

    #[test]
    fn test_unsupported_gap_penalties() {
        assert!(ScoreMatrix::new("blosum62", 3, 7, 0, 1, 0).is_err());
    }
}
