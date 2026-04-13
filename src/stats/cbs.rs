use crate::basic::value::{Letter, LETTER_MASK, TRUE_AA};
use super::score_matrix::ScoreMatrix;

/// Composition-based statistics mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbsMode {
    Disabled = 0,
    Hauser = 1,
    Deprecated1 = 2,
    HauserAndMatrixAdjust = 3,
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
            _ => CbsMode::Hauser,
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

/// Accumulator for vector scores in a sliding window.
struct VectorScores {
    scores: [i32; 20],
}

impl VectorScores {
    fn new() -> Self {
        VectorScores { scores: [0; 20] }
    }

    fn add(&mut self, letter: Letter, score_matrix: &ScoreMatrix) {
        let l = (letter & LETTER_MASK) as usize;
        if l < 20 {
            for i in 0..20 {
                self.scores[i] += score_matrix.score(l as Letter, i as Letter);
            }
        }
    }

    fn sub(&mut self, letter: Letter, score_matrix: &ScoreMatrix) {
        let l = (letter & LETTER_MASK) as usize;
        if l < 20 {
            for i in 0..20 {
                self.scores[i] -= score_matrix.score(l as Letter, i as Letter);
            }
        }
    }
}

/// Compute background scores: for each letter, the expected score
/// against the BLOSUM62 background frequencies.
pub fn compute_background_scores(score_matrix: &ScoreMatrix) -> [f64; 20] {
    // BLOSUM62 background frequencies (exact values from C++ blosum62.h lines 240-246).
    // Order: A R N D C Q E G H I L K M F P S T W Y V
    let bg_freq: [f64; 20] = [
        7.4216205067993410e-02, 5.1614486141284638e-02, 4.4645808512757915e-02,
        5.3626000838554413e-02, 2.4687457167944848e-02, 3.4259650591416023e-02,
        5.4311925684587502e-02, 7.4146941452644999e-02, 2.6212984805266227e-02,
        6.7917367618953756e-02, 9.8907868497150955e-02, 5.8155682303079680e-02,
        2.4990197579643110e-02, 4.7418459742284751e-02, 3.8538003320306206e-02,
        5.7229029476494421e-02, 5.0891364550287033e-02, 1.3029956129972148e-02,
        3.2281512313758580e-02, 7.2919098205619245e-02,
    ];

    let mut bg_scores = [0.0f64; 20];
    for (i, bg) in bg_scores.iter_mut().enumerate() {
        for (j, &freq) in bg_freq.iter().enumerate() {
            *bg += freq * score_matrix.score(i as Letter, j as Letter) as f64;
        }
    }
    bg_scores
}

/// Compute sliding-window Hauser correction for each position in a sequence.
///
/// This matches the C++ `HauserCorrection` constructor which uses a sliding window
/// of composition around each position to compute per-position score adjustments.
///
/// The default window size is 40 (matching C++ `config.cbs_window`).
pub fn hauser_correction(seq: &[Letter], score_matrix: &ScoreMatrix) -> Vec<i8> {
    hauser_correction_window(seq, score_matrix, 40)
}

/// Hauser correction with configurable window size.
pub fn hauser_correction_window(
    seq: &[Letter],
    score_matrix: &ScoreMatrix,
    window: usize,
) -> Vec<i8> {
    let len = seq.len();
    if len == 0 {
        return Vec::new();
    }

    let bg_scores = compute_background_scores(score_matrix);
    let mut corrections = vec![0.0f32; len];
    let mut vs = VectorScores::new();

    let window_half = window.min(2 * (len - 1)) / 2;
    let mut n: usize = 0;
    let mut h: usize = 0; // head (right edge of window)
    let mut m: usize = 0; // current position
    let mut t: usize = 0; // tail (left edge of window)

    // Phase 1: build initial half-window
    while n < window_half && h < len {
        n += 1;
        vs.add(seq[h], score_matrix);
        h += 1;
    }

    // Phase 2: expand to full window, start producing corrections
    while n < window + 1 && h < len {
        n += 1;
        vs.add(seq[h], score_matrix);
        let r = (seq[m] & LETTER_MASK) as usize;
        if r < 20 {
            let self_score = score_matrix.score(r as Letter, r as Letter);
            corrections[m] =
                bg_scores[r] as f32 - (vs.scores[r] - self_score) as f32 / (n - 1) as f32;
        }
        h += 1;
        m += 1;
    }

    // Phase 3: slide the full window
    while h < len {
        vs.add(seq[h], score_matrix);
        vs.sub(seq[t], score_matrix);
        let r = (seq[m] & LETTER_MASK) as usize;
        if r < 20 {
            let self_score = score_matrix.score(r as Letter, r as Letter);
            corrections[m] =
                bg_scores[r] as f32 - (vs.scores[r] - self_score) as f32 / (n - 1) as f32;
        }
        h += 1;
        t += 1;
        m += 1;
    }

    // Phase 4: shrink window from the right
    while m < len && n > window_half + 1 {
        n -= 1;
        vs.sub(seq[t], score_matrix);
        let r = (seq[m] & LETTER_MASK) as usize;
        if r < 20 {
            let self_score = score_matrix.score(r as Letter, r as Letter);
            corrections[m] =
                bg_scores[r] as f32 - (vs.scores[r] - self_score) as f32 / (n - 1) as f32;
        }
        t += 1;
        m += 1;
    }

    // Phase 5: remaining positions with shrinking window
    while m < len {
        let r = (seq[m] & LETTER_MASK) as usize;
        if r < 20 && n > 1 {
            let self_score = score_matrix.score(r as Letter, r as Letter);
            corrections[m] =
                bg_scores[r] as f32 - (vs.scores[r] - self_score) as f32 / (n - 1) as f32;
        }
        m += 1;
    }

    // Convert to i8 with rounding (matching C++)
    corrections
        .iter()
        .map(|&f| {
            if f < 0.0 {
                (f - 0.5) as i8
            } else {
                (f + 0.5) as i8
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_composition() {
        let seq = vec![0i8, 0, 1, 1, 1, 2];
        let comp = compute_composition(&seq);
        assert!((comp[0] - 2.0 / 6.0).abs() < 0.001);
        assert!((comp[1] - 3.0 / 6.0).abs() < 0.001);
    }

    #[test]
    fn test_hauser_correction_length() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq: Vec<Letter> = (0..100).map(|i| (i % 20) as Letter).collect();
        let corr = hauser_correction(&seq, &sm);
        assert_eq!(corr.len(), 100);
    }

    #[test]
    fn test_hauser_correction_empty() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let corr = hauser_correction(&[], &sm);
        assert!(corr.is_empty());
    }

    #[test]
    fn test_hauser_correction_short() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq = vec![0i8, 1, 2];
        let corr = hauser_correction(&seq, &sm);
        assert_eq!(corr.len(), 3);
    }

    #[test]
    fn test_cbs_mode_parse() {
        assert_eq!(CbsMode::parse("0"), CbsMode::Disabled);
        assert_eq!(CbsMode::parse("1"), CbsMode::Hauser);
        assert!(CbsMode::Hauser.uses_hauser());
        assert!(!CbsMode::Disabled.uses_hauser());
    }

    #[test]
    fn test_background_scores() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let bg = compute_background_scores(&sm);
        // Background scores should be negative (expected score against random)
        for &s in &bg {
            assert!(s < 1.0, "Background score should be small: {}", s);
        }
    }
}
