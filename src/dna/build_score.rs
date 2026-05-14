use std::f64::consts::LN_2;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Blastn_Score {
    lambda: f64,
    k: f64,
    h: f64,
    log_k: f64,
    reward_: i32,
    penalty_: i32,
    gap_open_: i32,
    gap_extend_: i32,
    target_length_: u64,
    db_size_: i64,
}

#[allow(non_snake_case)]
impl Blastn_Score {
    pub fn new(
        reward: i32,
        penalty: i32,
        gapopen: i32,
        gapextend: i32,
        db_letters: u64,
        sequence_count: i64,
        lambda: f64,
        k: f64,
        h: f64,
    ) -> Self {
        Self {
            lambda,
            k,
            h,
            log_k: k.ln(),
            reward_: reward,
            penalty_: penalty,
            gap_open_: gapopen,
            gap_extend_: gapextend,
            target_length_: db_letters,
            db_size_: sequence_count,
        }
    }

    pub fn blast_bit_Score(&self, raw_score: i32) -> f64 {
        ((raw_score as f64 * self.lambda) - self.log_k) / LN_2
    }

    pub fn blast_eValue(&self, raw_score: i32, query_length: i32) -> f64 {
        let searchspace = self.calculate_length_adjustment(
            query_length as u64,
            query_length,
            self.target_length_,
        ) * self.calculate_length_adjustment_db(
            self.target_length_,
            query_length,
            self.target_length_,
            self.db_size_,
        );
        self.blast_karlin_stoe_simple(raw_score, searchspace as i64)
    }

    pub fn reward(&self) -> i32 {
        self.reward_
    }

    pub fn penalty(&self) -> i32 {
        self.penalty_
    }

    pub fn gap_open(&self) -> i32 {
        self.gap_open_
    }

    pub fn gap_extend(&self) -> i32 {
        self.gap_extend_
    }

    pub fn calculate_length_adjustment_db(
        &self,
        length: u64,
        query_length: i32,
        target_length: u64,
        db_size: i64,
    ) -> f64 {
        length as f64 - (self.expected_hsp_value(query_length, target_length) * db_size as f64)
    }

    pub fn calculate_length_adjustment(
        &self,
        length: u64,
        query_length: i32,
        target_length: u64,
    ) -> f64 {
        if self.expected_hsp_value(query_length, target_length) < (1.0 / self.k) {
            1.0 / self.k
        } else {
            length as f64 - self.expected_hsp_value(query_length, target_length)
        }
    }

    pub fn expected_hsp_value(&self, query_length: i32, target_length: u64) -> f64 {
        (self.k * query_length as f64 * target_length as f64).ln() / self.h
    }

    pub fn blast_karlin_stoe_simple(&self, raw_score: i32, searchsp: i64) -> f64 {
        if self.lambda < 0.0 || self.k < 0.0 || self.h < 0.0 {
            return -1.0;
        }
        searchsp as f64 * (-(self.lambda * raw_score as f64) + self.log_k).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn score() -> Blastn_Score {
        Blastn_Score::new(2, -3, 5, 2, 1_000_000, 2_000, 0.625, 0.41, 1.1)
    }

    #[test]
    fn test_blastn_score_accessors() {
        let s = score();
        assert_eq!(s.reward(), 2);
        assert_eq!(s.penalty(), -3);
        assert_eq!(s.gap_open(), 5);
        assert_eq!(s.gap_extend(), 2);
    }

    #[test]
    fn test_blast_bit_score_matches_cpp_formula() {
        let s = score();
        let expected = ((80.0 * 0.625) - 0.41_f64.ln()) / LN_2;
        assert!((s.blast_bit_Score(80) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_length_adjustments_and_evalue_match_cpp_formula() {
        let s = score();
        let expected_hsp = (0.41_f64 * 100.0 * 1_000_000.0).ln() / 1.1;
        assert!((s.expected_hsp_value(100, 1_000_000) - expected_hsp).abs() < 1e-12);

        let q_eff = if expected_hsp < 1.0 / 0.41 {
            1.0 / 0.41
        } else {
            100.0 - expected_hsp
        };
        let db_eff = 1_000_000.0 - expected_hsp * 2_000.0;
        assert!((s.calculate_length_adjustment(100, 100, 1_000_000) - q_eff).abs() < 1e-12);
        assert!(
            (s.calculate_length_adjustment_db(1_000_000, 100, 1_000_000, 2_000) - db_eff).abs()
                < 1e-12
        );

        let searchspace = (q_eff * db_eff) as i64;
        let expected_e = searchspace as f64 * (-(0.625 * 80.0) + 0.41_f64.ln()).exp();
        assert!((s.blast_eValue(80, 100) - expected_e).abs() < expected_e * 1e-12);
    }
}
