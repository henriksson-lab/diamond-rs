use super::score_matrix::ScoreMatrix;
use super::target_freq::blast_optimize_target_frequencies;
use crate::basic::value::{Letter, AMINO_ACID_COUNT, LETTER_MASK, MASK_LETTER, TRUE_AA};

const BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT: f64 = 1.0e-5;
const BLAST_KARLIN_LAMBDA_ITER_DEFAULT: i32 = 17;
const COMPO_SCORE_MIN: f64 = -128.0;
const LAMBDA_RATIO_LOWER_BOUND: f64 = 0.5;
const K_MAXIMUM_X_SCORE: f64 = -1.0;
const RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS: f64 = 20.0;
const HIGH_PAIR_THRESHOLD: f64 = 0.4;
const LENGTH_LOWER_THRESHOLD: i32 = 50;
const HALF_CIRCLE_DEGREES: f64 = 180.0;
const PI: f64 = 3.1415926543;
pub const FIXED_RE_BLOSUM62: f64 = 0.44;

pub const NCBI_ALPH: usize = 28;
pub const ALPH_TO_NCBI: [usize; TRUE_AA as usize] = [
    1, 16, 13, 4, 3, 15, 5, 7, 8, 9, 11, 10, 12, 6, 14, 17, 18, 20, 22, 19,
];

/// Composition-based statistics mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbsMode {
    Disabled = 0,
    Hauser = 1,
    Deprecated1 = 2,
    HauserAndMatrixAdjust = 3,
    MatrixAdjust = 4,
    CompBasedStatsAndMatrixAdjust = 5,
    ConditionalMatrixAdjust = 6,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixAdjustRule {
    DontAdjustMatrix = -1,
    CompoScaleOldMatrix = 0,
    UnconstrainedRelEntropy = 1,
    RelEntropyOldMatrixNewContext = 2,
    RelEntropyOldMatrixOldContext = 3,
    UserSpecifiedRelEntropy = 4,
}

#[derive(Debug, Clone, Copy)]
pub struct CbsThresholds {
    pub query_match_distance_threshold: f64,
    pub length_ratio_threshold: f64,
    pub angle: f64,
}

impl CbsThresholds {
    pub fn new(
        query_match_distance_threshold: f64,
        length_ratio_threshold: f64,
        angle: f64,
    ) -> Self {
        let mut thresholds = CbsThresholds {
            query_match_distance_threshold: -1.0,
            length_ratio_threshold: -1.0,
            angle: 50.0,
        };
        if angle != -1.0 {
            thresholds.angle = angle;
        }
        if query_match_distance_threshold != 1.0 {
            thresholds.query_match_distance_threshold = query_match_distance_threshold;
        }
        if length_ratio_threshold != -1.0 {
            thresholds.length_ratio_threshold = length_ratio_threshold;
        }
        thresholds
    }
}

impl Default for CbsThresholds {
    fn default() -> Self {
        CbsThresholds::new(-1.0, -1.0, -1.0)
    }
}

impl CbsMode {
    pub fn parse(s: &str) -> Self {
        match s {
            "0" => CbsMode::Disabled,
            "1" => CbsMode::Hauser,
            "2" => CbsMode::Deprecated1,
            "3" => CbsMode::HauserAndMatrixAdjust,
            "4" => CbsMode::MatrixAdjust,
            "5" => CbsMode::CompBasedStatsAndMatrixAdjust,
            "6" => CbsMode::ConditionalMatrixAdjust,
            _ => CbsMode::Hauser,
        }
    }

    pub fn uses_hauser(&self) -> bool {
        self.hauser()
    }

    pub fn hauser(self) -> bool {
        matches!(
            self,
            CbsMode::Hauser | CbsMode::Deprecated1 | CbsMode::HauserAndMatrixAdjust
        )
    }

    pub fn matrix_adjust(self) -> bool {
        matches!(
            self,
            CbsMode::Deprecated1
                | CbsMode::HauserAndMatrixAdjust
                | CbsMode::MatrixAdjust
                | CbsMode::ConditionalMatrixAdjust
                | CbsMode::CompBasedStatsAndMatrixAdjust
        )
    }

    pub fn support_translated(self) -> bool {
        matches!(self, CbsMode::Disabled | CbsMode::Hauser)
    }

    pub fn conditioned(self) -> bool {
        matches!(
            self,
            CbsMode::Deprecated1
                | CbsMode::HauserAndMatrixAdjust
                | CbsMode::ConditionalMatrixAdjust
                | CbsMode::CompBasedStatsAndMatrixAdjust
        )
    }

    pub fn tantan(self) -> i32 {
        match self {
            CbsMode::Disabled | CbsMode::Hauser => 1,
            _ => 0,
        }
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

pub fn count_true_aa(seq: &[Letter]) -> i32 {
    seq.iter()
        .filter(|&&letter| i32::from(letter & LETTER_MASK) < TRUE_AA)
        .count() as i32
}

pub fn use_seg_masking(_a: &[Letter], _b: &[Letter]) -> bool {
    true
}

pub fn high_pair_frequencies(letter_probs: &[f64], length: i32) -> bool {
    if length <= LENGTH_LOWER_THRESHOLD {
        return false;
    }
    let mut max = 0.0;
    let mut second = 0.0;
    for &p in letter_probs.iter().take(TRUE_AA as usize) {
        if p > second {
            second = p;
            if p > max {
                second = max;
                max = p;
            }
        }
    }
    max + second > HIGH_PAIR_THRESHOLD
}

pub fn high_pair_either_seq(
    query_probs: &[f64],
    query_len: i32,
    match_probs: &[f64],
    match_len: i32,
) -> bool {
    high_pair_frequencies(query_probs, query_len) || high_pair_frequencies(match_probs, match_len)
}

pub fn relative_entropy(a: &[f64], b: &[f64]) -> f64 {
    let mut value = 0.0;
    for i in 0..TRUE_AA as usize {
        let temp = (a[i] + b[i]) / 2.0;
        if temp > 0.0 {
            if a[i] > 0.0 {
                value += a[i] * (a[i] / temp).ln() / 2.0;
            }
            if b[i] > 0.0 {
                value += b[i] * (b[i] / temp).ln() / 2.0;
            }
        }
    }
    if value < 0.0 {
        value = 0.0;
    }
    value.sqrt()
}

pub fn test_to_apply_re_adjustment_conditional(
    query_len: i32,
    match_len: i32,
    query_probs: &[f64],
    match_probs: &[f64],
    background_freqs: &[f64],
    thresholds: CbsThresholds,
) -> MatrixAdjustRule {
    let p_matrix = background_freqs;
    let d_m_mat = relative_entropy(match_probs, p_matrix);
    let d_q_mat = relative_entropy(query_probs, p_matrix);
    let d_m_q = relative_entropy(match_probs, query_probs);

    let denom = 2.0 * d_m_mat * d_q_mat;
    let mut angle = if denom == 0.0 {
        0.0
    } else {
        ((d_m_mat * d_m_mat + d_q_mat * d_q_mat - d_m_q * d_m_q) / denom).acos()
    };
    angle = angle * HALF_CIRCLE_DEGREES / PI;

    let len_q = query_len as f64;
    let len_m = match_len as f64;
    let len_large = len_q.max(len_m);
    let len_small = len_q.min(len_m);

    if high_pair_either_seq(query_probs, query_len, match_probs, match_len) {
        MatrixAdjustRule::UserSpecifiedRelEntropy
    } else if d_m_q > thresholds.query_match_distance_threshold
        && len_large / len_small > thresholds.length_ratio_threshold
        && angle > thresholds.angle
    {
        MatrixAdjustRule::CompoScaleOldMatrix
    } else {
        MatrixAdjustRule::UserSpecifiedRelEntropy
    }
}

pub fn adjust_matrix(
    query_comp: &[f64; TRUE_AA as usize],
    query_len: i32,
    cbs: CbsMode,
    target: &[Letter],
    background_freqs: &[f64],
    thresholds: CbsThresholds,
) -> MatrixAdjustRule {
    if !cbs.matrix_adjust() || target.is_empty() || query_len == 0 {
        return MatrixAdjustRule::DontAdjustMatrix;
    }

    let target_comp = compute_composition(target);
    if cbs.conditioned() {
        let rule = test_to_apply_re_adjustment_conditional(
            query_len,
            target.len() as i32,
            query_comp,
            &target_comp,
            background_freqs,
            thresholds,
        );
        if cbs == CbsMode::CompBasedStatsAndMatrixAdjust {
            rule
        } else if rule == MatrixAdjustRule::UserSpecifiedRelEntropy {
            MatrixAdjustRule::UserSpecifiedRelEntropy
        } else {
            MatrixAdjustRule::DontAdjustMatrix
        }
    } else {
        MatrixAdjustRule::UserSpecifiedRelEntropy
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetMatrix {
    pub scores: Vec<i8>,
    pub score_min: i32,
    pub score_max: i32,
}

impl TargetMatrix {
    pub fn new(scores: Vec<i8>, score_min: i32, score_max: i32) -> Self {
        TargetMatrix {
            scores,
            score_min,
            score_max,
        }
    }

    pub fn from_hauser_global(
        query_comp: &[f64; TRUE_AA as usize],
        target_comp: &[f64; TRUE_AA as usize],
        score_matrix: &ScoreMatrix,
    ) -> Self {
        let adjusted = hauser_global(query_comp, target_comp, score_matrix);
        let mut scores = vec![0i8; 32 * AMINO_ACID_COUNT];
        let mut score_min = i32::MAX;
        let mut score_max = i32::MIN;
        for i in 0..AMINO_ACID_COUNT {
            for j in 0..AMINO_ACID_COUNT {
                let s = adjusted[i * AMINO_ACID_COUNT + j];
                scores[i * 32 + j] = s.clamp(i8::MIN as i32, i8::MAX as i32) as i8;
                score_min = score_min.min(s);
                score_max = score_max.max(s);
            }
        }
        TargetMatrix {
            scores,
            score_min,
            score_max,
        }
    }

    pub fn score_width(&self) -> usize {
        if self.score_max > i8::MAX as i32 || self.score_min < i8::MIN as i32 {
            1
        } else {
            0
        }
    }
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
        // C++ VectorScores::operator+= does NOT guard for l < 20 — all letters
        // including X (23) contribute their scores to the window sum.
        // C++ derives `l` via `Sequence::operator[]` which applies `LETTER_MASK`
        // under SEQ_MASK (`diamond/CMakeLists.txt:71` makes SEQ_MASK on by
        // default). Mask here too so a raw SEED_MASK-bearing byte (0x80 / -128)
        // doesn't sign-extend through `as usize` into a wild index.
        let l = (letter & crate::basic::value::LETTER_MASK) as usize;
        for i in 0..20 {
            self.scores[i] += score_matrix.score(l as Letter, i as Letter);
        }
    }

    fn sub(&mut self, letter: Letter, score_matrix: &ScoreMatrix) {
        // See `add` for the LETTER_MASK rationale.
        let l = (letter & crate::basic::value::LETTER_MASK) as usize;
        for i in 0..20 {
            self.scores[i] -= score_matrix.score(l as Letter, i as Letter);
        }
    }
}

/// Compute background scores: for each letter, the expected score
/// against the BLOSUM62 background frequencies.
pub fn compute_background_scores(score_matrix: &ScoreMatrix) -> [f64; 20] {
    *score_matrix.background_scores()
}

pub fn hauser_global(
    query_comp: &[f64; TRUE_AA as usize],
    target_comp: &[f64; TRUE_AA as usize],
    score_matrix: &ScoreMatrix,
) -> Vec<i32> {
    let background_scores = compute_background_scores(score_matrix);
    let mut qscores = [0.0f64; TRUE_AA as usize];
    let mut tscores = [0.0f64; TRUE_AA as usize];

    for i in 0..TRUE_AA as usize {
        for j in 0..TRUE_AA as usize {
            qscores[i] += query_comp[j] * score_matrix.score(i as Letter, j as Letter) as f64;
            tscores[i] += target_comp[j] * score_matrix.score(i as Letter, j as Letter) as f64;
        }
    }
    for i in 0..TRUE_AA as usize {
        qscores[i] = background_scores[i] - qscores[i];
        tscores[i] = background_scores[i] - tscores[i];
    }

    let mut matrix = vec![0i32; AMINO_ACID_COUNT * AMINO_ACID_COUNT];
    for i in 0..AMINO_ACID_COUNT {
        for j in 0..AMINO_ACID_COUNT {
            let s = score_matrix.score(i as Letter, j as Letter) as f64;
            let q = if i < TRUE_AA as usize {
                qscores[i]
            } else {
                0.0
            };
            let t = if j < TRUE_AA as usize {
                tscores[j]
            } else {
                0.0
            };
            matrix[i * AMINO_ACID_COUNT + j] = (s + q.min(t)).round() as i32;
        }
    }
    matrix
}

fn blast_gcd(mut a: i32, mut b: i32) -> i32 {
    b = b.abs();
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

fn prob_at(probs: &[f64], min_score: i32, score: i32) -> f64 {
    probs[(score - min_score) as usize]
}

fn nlm_karlin_lambda_nr(
    probs: &[f64],
    d: i32,
    low: i32,
    high: i32,
    lambda0: f64,
    tolx: f64,
    itmax: i32,
    max_newton: i32,
) -> f64 {
    debug_assert!(d > 0);
    let x0 = (-lambda0).exp();
    let mut x = if 0.0 < x0 && x0 < 1.0 { x0 } else { 0.5 };
    let mut a = 0.0;
    let mut b = 1.0;
    let mut f = 4.0;
    let mut is_newton = false;

    for k in 0..itmax {
        let fold = f;
        let was_newton = is_newton;
        is_newton = false;

        let mut g = 0.0;
        f = prob_at(probs, low, low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + prob_at(probs, low, i);
            i += d;
        }
        g = x * g + f;
        f = f * x + prob_at(probs, low, 0) - 1.0;
        i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + prob_at(probs, low, i);
            i += d;
        }

        if f > 0.0 {
            a = x;
        } else if f < 0.0 {
            b = x;
        } else {
            break;
        }
        if b - a < 2.0 * a * (1.0 - b) * tolx {
            x = (a + b) / 2.0;
            break;
        }

        if k >= max_newton || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
            x = (a + b) / 2.0;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() < tolx * x * (1.0 - x) {
                    break;
                }
            }
        }
    }

    -x.ln() / d as f64
}

pub fn calc_lambda(probs: &[f64], min_score: i32, max_score: i32, lambda0: f64) -> f64 {
    let score_range = max_score - min_score + 1;
    debug_assert_eq!(probs.len(), score_range as usize);
    let mut avg = 0.0;
    for i in 0..score_range {
        avg += (min_score + i) as f64 * probs[i as usize];
    }
    if avg >= 0.0 {
        return -1.0;
    }

    let mut d = -min_score;
    let mut i = 1;
    while i <= max_score - min_score && d > 1 {
        if probs[i as usize] != 0.0 {
            d = blast_gcd(d, i);
        }
        i += 1;
    }

    nlm_karlin_lambda_nr(
        probs,
        d,
        min_score,
        max_score,
        lambda0,
        BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
        20,
        20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT,
    )
}

pub fn matrix_score_probs(
    matrix: &[i32],
    row_stride: usize,
    alphsize: usize,
    subject_probs: &[f64],
    query_probs: &[f64],
) -> (Vec<f64>, i32, i32) {
    let mut obs_min = 0;
    let mut obs_max = 0;
    for irow in 0..alphsize {
        for aa in 0..TRUE_AA as usize {
            let score = matrix[irow * row_stride + aa];
            obs_min = obs_min.min(score);
            obs_max = obs_max.max(score);
        }
    }

    let mut probs = vec![0.0f64; (obs_max - obs_min + 1) as usize];
    for irow in 0..alphsize {
        for aa in 0..TRUE_AA as usize {
            let score = matrix[irow * row_stride + aa];
            if score >= obs_min {
                probs[(score - obs_min) as usize] += query_probs[irow] * subject_probs[aa];
            }
        }
    }
    (probs, obs_min, obs_max)
}

pub fn freq_ratio_to_score(matrix: &mut [Vec<f64>], lambda: f64) {
    for row in matrix {
        for score in row {
            if *score == 0.0 {
                *score = COMPO_SCORE_MIN;
            } else {
                *score = score.ln() / lambda;
            }
        }
    }
}

pub fn round_score_matrix(float_scores: &[Vec<f64>]) -> Vec<Vec<i32>> {
    float_scores
        .iter()
        .map(|row| {
            row.iter()
                .map(|&score| {
                    if score < i32::MIN as f64 {
                        i32::MIN
                    } else {
                        score.round() as i32
                    }
                })
                .collect()
        })
        .collect()
}

pub fn apply_pseudocounts(
    probs20: &mut [f64],
    number_of_observations: i32,
    background_probs20: &[f64],
) {
    let mut sum = 0.0;
    for &p in probs20.iter().take(TRUE_AA as usize) {
        sum += p;
    }
    if sum == 0.0 {
        sum = 1.0;
    }
    let weight = RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS
        / (number_of_observations as f64 + RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS);
    for i in 0..TRUE_AA as usize {
        probs20[i] = (1.0 - weight) * probs20[i] / sum + weight * background_probs20[i];
    }
}

pub fn true_aa_to_std_target_freqs(
    std_freq: &mut [f64],
    stride: usize,
    std_alphsize: usize,
    freq: &[f64],
) {
    let mut sum = 0.0;
    for a in 0..TRUE_AA as usize {
        for b in 0..TRUE_AA as usize {
            sum += freq[a * TRUE_AA as usize + b];
        }
    }
    for a_big in 0..std_alphsize {
        if a_big >= TRUE_AA as usize {
            for b_big in 0..std_alphsize {
                std_freq[a_big * stride + b_big] = 0.0;
            }
        } else {
            for b_big in 0..std_alphsize {
                if b_big >= TRUE_AA as usize {
                    std_freq[a_big * stride + b_big] = 0.0;
                } else {
                    std_freq[a_big * stride + b_big] = freq[a_big * TRUE_AA as usize + b_big] / sum;
                }
            }
        }
    }
}

pub fn calc_freq_ratios(
    ratios: &mut [f64],
    stride: usize,
    alphsize: usize,
    row_prob: &[f64],
    col_prob: &[f64],
) {
    for i in 0..alphsize {
        if row_prob[i] > 0.0 {
            for j in 0..alphsize {
                if col_prob[j] > 0.0 {
                    ratios[i * stride + j] /= row_prob[i] * col_prob[j];
                }
            }
        }
    }
}

pub fn scores_std_alphabet(
    matrix: &mut [i32],
    matrix_stride: usize,
    alphsize: usize,
    target_freq: &[f64],
    row_prob: &[f64],
    col_prob: &[f64],
    lambda: f64,
) -> Result<(), ()> {
    let mut scores = vec![0.0f64; alphsize * alphsize];
    true_aa_to_std_target_freqs(&mut scores, alphsize, alphsize, target_freq);
    calc_freq_ratios(&mut scores, alphsize, TRUE_AA as usize, row_prob, col_prob);
    let mut rows = (0..alphsize)
        .map(|i| scores[i * alphsize..i * alphsize + alphsize].to_vec())
        .collect::<Vec<_>>();
    freq_ratio_to_score(&mut rows, lambda);
    for i in 0..alphsize {
        for j in 0..alphsize {
            scores[i * alphsize + j] = rows[i][j];
        }
    }
    set_xuo_scores(&mut scores, alphsize, TRUE_AA as usize, row_prob, col_prob);
    for i in 0..alphsize {
        for j in 0..alphsize {
            let score = scores[i * alphsize + j];
            matrix[i * matrix_stride + j] = if score < i32::MIN as f64 {
                i32::MIN
            } else {
                score.round() as i32
            };
        }
    }
    Ok(())
}

pub fn blast_composition_matrix_adj(
    matrix: &mut [i32],
    matrix_stride: usize,
    matrix_adjust_rule: MatrixAdjustRule,
    length1: i32,
    length2: i32,
    stdaa_row_probs: &[f64],
    stdaa_col_probs: &[f64],
    lambda: f64,
    joint_probs: &[f64],
    background_freqs: &[f64],
    tol: f64,
    maxits: i32,
) -> Result<i32, ()> {
    let desired_re = match matrix_adjust_rule {
        MatrixAdjustRule::UserSpecifiedRelEntropy => FIXED_RE_BLOSUM62,
        _ => return Err(()),
    };
    let mut row_probs = [0.0f64; TRUE_AA as usize];
    let mut col_probs = [0.0f64; TRUE_AA as usize];
    row_probs.copy_from_slice(&stdaa_row_probs[..TRUE_AA as usize]);
    col_probs.copy_from_slice(&stdaa_col_probs[..TRUE_AA as usize]);

    apply_pseudocounts(&mut row_probs, length1, background_freqs);
    apply_pseudocounts(&mut col_probs, length2, background_freqs);

    let mut mat_final = vec![0.0f64; TRUE_AA as usize * TRUE_AA as usize];
    let (status, iteration_count) = blast_optimize_target_frequencies(
        &mut mat_final,
        TRUE_AA as usize,
        joint_probs,
        &row_probs,
        &col_probs,
        desired_re > 0.0,
        desired_re,
        tol,
        maxits,
    );
    if status != 0 {
        return Err(());
    }

    scores_std_alphabet(
        matrix,
        matrix_stride,
        AMINO_ACID_COUNT,
        &mat_final,
        &row_probs,
        &col_probs,
        lambda,
    )?;
    Ok(iteration_count)
}

pub fn composition_matrix_adjust(
    query_len: i32,
    target_len: i32,
    query_comp: &[f64],
    target_comp: &[f64],
    scale: i32,
    ungapped_lambda: f64,
    joint_probs: &[f64],
    background_freqs: &[f64],
    score_matrix: &ScoreMatrix,
    tol: f64,
    maxits: i32,
) -> Vec<i32> {
    let mut out = vec![0i32; AMINO_ACID_COUNT * AMINO_ACID_COUNT];
    let result = blast_composition_matrix_adj(
        &mut out,
        AMINO_ACID_COUNT,
        MatrixAdjustRule::UserSpecifiedRelEntropy,
        query_len,
        target_len,
        query_comp,
        target_comp,
        ungapped_lambda / scale as f64,
        joint_probs,
        background_freqs,
        tol,
        maxits,
    );
    if result.is_err() {
        for i in 0..AMINO_ACID_COUNT {
            for j in 0..AMINO_ACID_COUNT {
                out[i * AMINO_ACID_COUNT + j] =
                    score_matrix.score(i as Letter, j as Letter) * scale;
            }
        }
    }
    out
}

pub fn calc_avg_score(
    scores: &[f64],
    offset: usize,
    alphsize: usize,
    inc: usize,
    probs: &[f64],
) -> f64 {
    let mut score_ix = 0.0;
    for j in 0..alphsize {
        score_ix += scores[offset + j * inc] * probs[j];
    }
    score_ix
}

pub fn calc_x_score(
    scores: &[f64],
    offset: usize,
    alphsize: usize,
    inc: usize,
    probs: &[f64],
) -> f64 {
    calc_avg_score(scores, offset, alphsize, inc, probs).min(K_MAXIMUM_X_SCORE)
}

pub fn set_xuo_scores(
    scores: &mut [f64],
    stride: usize,
    alphsize: usize,
    row_probs: &[f64],
    col_probs: &[f64],
) {
    let mask = MASK_LETTER as usize;
    let mut score_xx = 0.0;
    for i in 0..alphsize {
        let avg_ix = calc_avg_score(scores, i * stride, alphsize, 1, col_probs);
        scores[i * stride + mask] = avg_ix.min(K_MAXIMUM_X_SCORE);
        score_xx += avg_ix * row_probs[i];

        scores[mask * stride + i] = calc_x_score(scores, i, alphsize, stride, row_probs);
    }
    scores[mask * stride + mask] = score_xx.min(K_MAXIMUM_X_SCORE);
}

pub fn scale_square_matrix(
    matrix: &mut [i32],
    stride: usize,
    alphsize: usize,
    row_probs: &[f64],
    col_probs: &[f64],
    lambda: f64,
    freq_ratios: &[[f64; NCBI_ALPH]; NCBI_ALPH],
) -> Result<(), ()> {
    let mut scores = vec![0.0f64; alphsize * stride];
    for i in 0..TRUE_AA as usize {
        for j in 0..TRUE_AA as usize {
            scores[i * stride + j] = freq_ratios[ALPH_TO_NCBI[i]][ALPH_TO_NCBI[j]];
        }
    }
    let mut rows = (0..alphsize)
        .map(|i| scores[i * stride..i * stride + alphsize].to_vec())
        .collect::<Vec<_>>();
    freq_ratio_to_score(&mut rows, lambda);
    for i in 0..alphsize {
        for j in 0..alphsize {
            scores[i * stride + j] = rows[i][j];
        }
    }
    set_xuo_scores(&mut scores, stride, TRUE_AA as usize, row_probs, col_probs);
    for i in 0..alphsize {
        for j in 0..alphsize {
            let score = scores[i * stride + j];
            matrix[i * stride + j] = if score < i32::MIN as f64 {
                i32::MIN
            } else {
                score.round() as i32
            };
        }
    }
    Ok(())
}

pub fn blast_composition_based_stats(
    matrix: &mut [i32],
    stride: usize,
    matrix_in: &[i32],
    matrix_in_stride: usize,
    query_probs: &[f64],
    res_probs: &[f64],
    lambda: f64,
    cbs_matrix_scale: f64,
    freq_ratios: &[[f64; NCBI_ALPH]; NCBI_ALPH],
) -> Result<f64, ()> {
    let (score_array, obs_min, obs_max) = matrix_score_probs(
        matrix_in,
        matrix_in_stride,
        TRUE_AA as usize,
        res_probs,
        query_probs,
    );
    let ungapped_lambda = lambda / cbs_matrix_scale;
    let correct_ungapped_lambda = calc_lambda(&score_array, obs_min, obs_max, ungapped_lambda);
    if correct_ungapped_lambda < 0.0 {
        return Err(());
    }

    let mut lambda_ratio = correct_ungapped_lambda / ungapped_lambda;
    lambda_ratio = lambda_ratio.min(1.0);
    lambda_ratio = lambda_ratio.max(LAMBDA_RATIO_LOWER_BOUND);

    if lambda_ratio > 0.0 {
        let scaled_lambda = ungapped_lambda / lambda_ratio;
        scale_square_matrix(
            matrix,
            stride,
            AMINO_ACID_COUNT,
            query_probs,
            res_probs,
            scaled_lambda,
            freq_ratios,
        )?;
    }
    Ok(lambda_ratio)
}

pub fn ideal_lambda(score_matrix: &ScoreMatrix) -> Option<f64> {
    let robinson = [
        (0usize, 78.05),
        (4, 19.25),
        (3, 53.64),
        (6, 62.95),
        (13, 38.56),
        (7, 73.77),
        (8, 21.99),
        (9, 51.42),
        (11, 57.44),
        (10, 90.19),
        (12, 22.43),
        (2, 44.87),
        (14, 52.03),
        (5, 42.64),
        (1, 51.29),
        (15, 71.20),
        (16, 58.41),
        (19, 64.41),
        (17, 13.30),
        (18, 32.16),
    ];
    let mut bg = [0.0f64; TRUE_AA as usize];
    let mut sum = 0.0;
    for &(letter, p) in &robinson {
        bg[letter] = p;
        sum += p;
    }
    for p in &mut bg {
        *p /= sum;
    }
    let (probs, obs_min, obs_max) =
        matrix_score_probs(score_matrix.matrix32(), 32, TRUE_AA as usize, &bg, &bg);
    let lambda = calc_lambda(&probs, obs_min, obs_max, 0.5);
    if lambda < 0.0 {
        None
    } else {
        Some(lambda)
    }
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

    // Convert to i8 with rounding (matching C++). Append 32 zero bytes of
    // trailing padding so SIMD score-profile loaders downstream can over-read
    // up to one AVX2 vector past `len` without UB. C++ does this explicitly
    // (`diamond/src/stats/hauser_correction.cpp:107-110`: `reserve(len + PADDING)`
    // followed by `for (Loc i=0; i<PADDING; ++i) int8.push_back(0)`, with
    // `PADDING = 32`). Without it, SIMD `bias` lane reads past `int8.size()`
    // hit either uninitialised heap memory or panic on bounds-checked indexing.
    //
    // Callers that index the slice positionally (the current `query_cbs[qi]`
    // pattern in blastp.rs) get the same numerical result either way since
    // they only touch `0..len`; the padding is purely for over-read safety.
    const PADDING: usize = 32;
    let mut out: Vec<i8> = Vec::with_capacity(len + PADDING);
    for &f in &corrections {
        let v = if f < 0.0 { f - 0.5 } else { f + 0.5 };
        out.push(v as i8);
    }
    for _ in 0..PADDING {
        out.push(0);
    }
    out
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
        assert_eq!(count_true_aa(&seq), 6);
        assert_eq!(count_true_aa(&[0, 1, 23, 30]), 2);
    }

    #[test]
    fn test_hauser_correction_length() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq: Vec<Letter> = (0..100).map(|i| (i % 20) as Letter).collect();
        let corr = hauser_correction(&seq, &sm);
        // C++ appends 32 zero bytes of SIMD over-read padding after `len`
        // (`hauser_correction.cpp:107-110`); Rust mirrors that, so the
        // returned slice is `len + 32` long. The numerically meaningful
        // values are in `corr[..seq.len()]`; the padding is all zero.
        assert_eq!(corr.len(), 100 + 32);
        assert!(corr[100..].iter().all(|&b| b == 0));
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
        assert_eq!(corr.len(), 3 + 32);
        assert!(corr[3..].iter().all(|&b| b == 0));
    }

    #[test]
    fn test_cbs_mode_parse() {
        assert_eq!(CbsMode::parse("0"), CbsMode::Disabled);
        assert_eq!(CbsMode::parse("1"), CbsMode::Hauser);
        assert!(CbsMode::Hauser.uses_hauser());
        assert!(!CbsMode::Disabled.uses_hauser());
        assert!(CbsMode::MatrixAdjust.matrix_adjust());
        assert!(CbsMode::ConditionalMatrixAdjust.conditioned());
        assert!(CbsMode::Hauser.support_translated());
        assert_eq!(CbsMode::MatrixAdjust.tantan(), 0);
    }

    #[test]
    fn test_matrix_adjust_rule_selection() {
        let q = compute_composition(&[0, 1, 2, 3, 4, 5, 6, 7]);
        let target: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let background = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        assert_eq!(
            adjust_matrix(
                &q,
                q.len() as i32,
                CbsMode::Disabled,
                &target,
                &background,
                CbsThresholds::default()
            ),
            MatrixAdjustRule::DontAdjustMatrix
        );
        assert_eq!(
            adjust_matrix(
                &q,
                q.len() as i32,
                CbsMode::MatrixAdjust,
                &target,
                &background,
                CbsThresholds::default()
            ),
            MatrixAdjustRule::UserSpecifiedRelEntropy
        );
        assert!(!high_pair_frequencies(&q, 8));

        let mut skewed = [0.0; TRUE_AA as usize];
        skewed[0] = 0.25;
        skewed[1] = 0.20;
        assert!(high_pair_frequencies(&skewed, 100));
        assert!(relative_entropy(&q, &background) >= 0.0);
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

    #[test]
    fn test_hauser_global_and_target_matrix() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let q = compute_composition(&[0, 0, 1, 2, 3, 4]);
        let t = compute_composition(&[0, 1, 1, 2, 2, 3]);
        let adjusted = hauser_global(&q, &t, &sm);
        assert_eq!(adjusted.len(), AMINO_ACID_COUNT * AMINO_ACID_COUNT);
        let tm = TargetMatrix::from_hauser_global(&q, &t, &sm);
        assert_eq!(tm.scores.len(), 32 * AMINO_ACID_COUNT);
        assert_eq!(tm.score_width(), 0);
        assert!(use_seg_masking(&[], &[]));
    }

    #[test]
    fn test_lambda_and_score_prob_helpers() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let uniform = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let (probs, obs_min, obs_max) =
            matrix_score_probs(sm.matrix32(), 32, TRUE_AA as usize, &uniform, &uniform);
        assert_eq!(probs.len(), (obs_max - obs_min + 1) as usize);
        assert!((probs.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        assert!(calc_lambda(&probs, obs_min, obs_max, 0.5) > 0.0);
        assert!(ideal_lambda(&sm).unwrap() > 0.0);
    }

    #[test]
    fn test_freq_ratio_round_score_matrix() {
        let mut scores = vec![vec![1.0, std::f64::consts::E, 0.0]];
        freq_ratio_to_score(&mut scores, 1.0);
        assert_eq!(scores[0][0], 0.0);
        assert!((scores[0][1] - 1.0).abs() < 1e-12);
        assert_eq!(scores[0][2], COMPO_SCORE_MIN);
        let rounded = round_score_matrix(&scores);
        assert_eq!(rounded[0], vec![0, 1, COMPO_SCORE_MIN as i32]);
    }

    #[test]
    fn test_xuo_and_scale_square_matrix() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let mut freq_ratios = [[1.0f64; NCBI_ALPH]; NCBI_ALPH];
        for i in 0..TRUE_AA as usize {
            for j in 0..TRUE_AA as usize {
                freq_ratios[ALPH_TO_NCBI[i]][ALPH_TO_NCBI[j]] =
                    (sm.score(i as Letter, j as Letter) as f64).exp();
            }
        }

        let uniform = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let mut matrix = vec![0i32; AMINO_ACID_COUNT * 32];
        scale_square_matrix(
            &mut matrix,
            32,
            AMINO_ACID_COUNT,
            &uniform,
            &uniform,
            1.0,
            &freq_ratios,
        )
        .unwrap();
        assert_eq!(matrix[0], sm.score(0, 0));
        assert_eq!(matrix[0 * 32 + 1], sm.score(0, 1));
        assert!(matrix[0 * 32 + MASK_LETTER as usize] <= K_MAXIMUM_X_SCORE as i32);
        assert!(matrix[MASK_LETTER as usize * 32] <= K_MAXIMUM_X_SCORE as i32);
    }

    #[test]
    fn test_blast_composition_based_stats_wrapper() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let mut freq_ratios = [[1.0f64; NCBI_ALPH]; NCBI_ALPH];
        let lambda = sm.lambda();
        for i in 0..TRUE_AA as usize {
            for j in 0..TRUE_AA as usize {
                freq_ratios[ALPH_TO_NCBI[i]][ALPH_TO_NCBI[j]] =
                    (sm.score(i as Letter, j as Letter) as f64 * lambda).exp();
            }
        }

        let uniform = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let mut out = vec![0i32; AMINO_ACID_COUNT * 32];
        let ratio = blast_composition_based_stats(
            &mut out,
            32,
            sm.matrix32(),
            32,
            &uniform,
            &uniform,
            lambda,
            1.0,
            &freq_ratios,
        )
        .unwrap();
        assert!((LAMBDA_RATIO_LOWER_BOUND..=1.0).contains(&ratio));
        assert!(out[0] > 0);
    }

    #[test]
    fn test_pseudocounts_and_target_freqs() {
        let background = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let mut probs = [0.0f64; TRUE_AA as usize];
        probs[0] = 10.0;
        apply_pseudocounts(&mut probs, 10, &background);
        assert!((probs.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(probs[0] > background[0]);
        assert!(probs[1] > 0.0);

        let mut target = vec![0.0f64; TRUE_AA as usize * TRUE_AA as usize];
        target[0] = 2.0;
        target[1] = 1.0;
        let mut std = vec![1.0f64; AMINO_ACID_COUNT * AMINO_ACID_COUNT];
        true_aa_to_std_target_freqs(&mut std, AMINO_ACID_COUNT, AMINO_ACID_COUNT, &target);
        assert!((std[0] - 2.0 / 3.0).abs() < 1e-12);
        assert!((std[1] - 1.0 / 3.0).abs() < 1e-12);
        assert_eq!(std[MASK_LETTER as usize * AMINO_ACID_COUNT], 0.0);
        assert_eq!(std[MASK_LETTER as usize], 0.0);
    }

    #[test]
    fn test_calc_freq_ratios_and_scores_std_alphabet() {
        let row = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let col = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let mut target = vec![0.0f64; TRUE_AA as usize * TRUE_AA as usize];
        for i in 0..TRUE_AA as usize {
            for j in 0..TRUE_AA as usize {
                target[i * TRUE_AA as usize + j] = row[i] * col[j];
            }
        }

        let mut ratios = target.clone();
        calc_freq_ratios(&mut ratios, TRUE_AA as usize, TRUE_AA as usize, &row, &col);
        for &ratio in &ratios {
            assert!((ratio - 1.0).abs() < 1e-12);
        }

        let mut matrix = vec![0i32; AMINO_ACID_COUNT * 32];
        scores_std_alphabet(&mut matrix, 32, AMINO_ACID_COUNT, &target, &row, &col, 1.0).unwrap();
        assert_eq!(matrix[0], 0);
        assert_eq!(matrix[1], 0);
        assert_eq!(
            matrix[0 * 32 + MASK_LETTER as usize],
            K_MAXIMUM_X_SCORE as i32
        );
        assert_eq!(matrix[MASK_LETTER as usize * 32], K_MAXIMUM_X_SCORE as i32);
        assert_eq!(FIXED_RE_BLOSUM62, 0.44);
    }

    #[test]
    fn test_composition_matrix_adjust_path_and_fallback() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let row = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let col = [1.0 / TRUE_AA as f64; TRUE_AA as usize];
        let mut joint = vec![0.0f64; TRUE_AA as usize * TRUE_AA as usize];
        for i in 0..TRUE_AA as usize {
            for j in 0..TRUE_AA as usize {
                joint[i * TRUE_AA as usize + j] = row[i] * col[j];
            }
        }

        let mut out = vec![0i32; AMINO_ACID_COUNT * AMINO_ACID_COUNT];
        let iterations = blast_composition_matrix_adj(
            &mut out,
            AMINO_ACID_COUNT,
            MatrixAdjustRule::UserSpecifiedRelEntropy,
            100,
            100,
            &row,
            &col,
            1.0,
            &joint,
            &row,
            1.0e9,
            10,
        )
        .unwrap();
        assert_eq!(iterations, 0);
        assert_eq!(out[0], 0);
        assert_eq!(out[MASK_LETTER as usize], K_MAXIMUM_X_SCORE as i32);

        let fallback =
            composition_matrix_adjust(100, 100, &row, &col, 2, 1.0, &joint, &row, &sm, 0.0, -1);
        assert_eq!(fallback[0], sm.score(0, 0) * 2);
        assert_eq!(fallback[1], sm.score(0, 1) * 2);
    }
}
