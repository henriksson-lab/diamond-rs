use crate::basic::value::{Letter, LETTER_MASK, TRUE_AA};
use super::score_matrix::ScoreMatrix;

/// Composition-based statistics mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbsMode {
    /// No composition adjustment.
    Disabled = 0,
    /// Hauser correction only.
    Hauser = 1,
    /// Legacy mode (deprecated).
    Deprecated1 = 2,
    /// Hauser + matrix adjustment.
    HauserAndMatrixAdjust = 3,
    /// Matrix adjustment only.
    MatrixAdjust = 4,
}

impl CbsMode {
    pub fn parse(s: &str) -> Self {
        match s {
            "0" => CbsMode::Disabled,
            "1" => CbsMode::Hauser,
            "2" => CbsMode::Deprecated1,
            "3" => CbsMode::HauserAndMatrixAdjust,
            "4" => CbsMode::MatrixAdjust,
            _ => CbsMode::Hauser, // default
        }
    }

    pub fn uses_hauser(&self) -> bool {
        matches!(self, CbsMode::Hauser | CbsMode::HauserAndMatrixAdjust)
    }
}

/// Compute amino acid composition frequencies for a sequence.
pub fn compute_composition(seq: &[Letter]) -> [f64; TRUE_AA as usize] {
    let mut counts = [0u32; TRUE_AA as usize];
    let mut total = 0u32;
    for &l in seq {
        let l = (l & LETTER_MASK) as usize;
        if l < TRUE_AA as usize {
            counts[l] += 1;
            total += 1;
        }
    }
    let mut freq = [0.0f64; TRUE_AA as usize];
    if total > 0 {
        for i in 0..TRUE_AA as usize {
            freq[i] = counts[i] as f64 / total as f64;
        }
    }
    freq
}

/// Compute Hauser correction values for each position in a sequence.
///
/// Returns a vector of per-position score adjustments (as i8).
pub fn hauser_correction(
    seq: &[Letter],
    score_matrix: &ScoreMatrix,
) -> Vec<i8> {
    let len = seq.len();
    let composition = compute_composition(seq);

    // Compute background score for each letter based on composition
    let mut bg_score = [0.0f64; TRUE_AA as usize];
    for (i, bg) in bg_score.iter_mut().enumerate() {
        for (j, &comp) in composition.iter().enumerate() {
            *bg += comp * score_matrix.score(i as Letter, j as Letter) as f64;
        }
    }

    // Per-position correction
    let mut correction = vec![0i8; len];
    for (pos, &letter) in seq.iter().enumerate() {
        let l = (letter & LETTER_MASK) as usize;
        if l < TRUE_AA as usize {
            // Correction = -(background_score for this letter)
            correction[pos] = (-bg_score[l].round()) as i8;
        }
    }

    correction
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_composition() {
        // Sequence of all A's (letter 0)
        let seq = vec![0i8; 100];
        let comp = compute_composition(&seq);
        assert!((comp[0] - 1.0).abs() < 0.001);
        assert!((comp[1] - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_hauser_correction() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq: Vec<Letter> = (0..20).map(|i| i as Letter).collect();
        let corr = hauser_correction(&seq, &sm);
        assert_eq!(corr.len(), 20);
    }

    #[test]
    fn test_cbs_mode_parse() {
        assert_eq!(CbsMode::parse("0"), CbsMode::Disabled);
        assert_eq!(CbsMode::parse("1"), CbsMode::Hauser);
        assert!(CbsMode::Hauser.uses_hauser());
        assert!(!CbsMode::Disabled.uses_hauser());
    }
}
