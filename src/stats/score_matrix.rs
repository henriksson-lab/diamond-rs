use super::matrices;
use super::standard_matrix::StandardMatrix;
use crate::basic::value::{Letter, AMINO_ACID_COUNT, TRUE_AA};

use std::f64::consts::LN_2;
use std::fmt;
#[cfg(not(windows))]
use std::fs;
use std::path::Path;

const BLOSUM62_BACKGROUND_FREQS: [f64; TRUE_AA as usize] = [
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

fn build_i8_scores(
    scores: &[i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT],
    stop_match_score: i32,
    bias: i8,
    modulo: usize,
    offset: usize,
) -> [i8; 32 * 32] {
    let mut data = [i8::MIN; 32 * 32];
    for i in 0..32 {
        for j in 0..32 {
            let j2 = j % modulo + offset;
            if i < AMINO_ACID_COUNT && j2 < AMINO_ACID_COUNT {
                data[i * 32 + j] = (scores[i * AMINO_ACID_COUNT + j2] as i16 + bias as i16) as i8;
            }
        }
    }
    if stop_match_score != 1 {
        data[24 * 32 + 24] = stop_match_score as i8;
    }
    data
}

fn build_i16_scores(
    scores: &[i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT],
    stop_match_score: i32,
) -> [i16; 32 * 32] {
    let mut data = [-128i16; 32 * 32];
    for i in 0..AMINO_ACID_COUNT {
        for j in 0..AMINO_ACID_COUNT {
            data[i * 32 + j] = scores[i * AMINO_ACID_COUNT + j] as i16;
        }
    }
    if stop_match_score != 1 {
        data[24 * 32 + 24] = stop_match_score as i16;
    }
    data
}

fn matrix_low_score(matrix: &[i32; 32 * 32]) -> i8 {
    let mut low = i8::MAX;
    for i in 0..AMINO_ACID_COUNT {
        for j in i + 1..AMINO_ACID_COUNT {
            low = low.min(matrix[i * 32 + j] as i8);
        }
    }
    low
}

/// Checked String compatibility for the Source byte-oriented matrix parser.
pub fn custom_scores_from_str(
    contents: &str,
    mask_score: i32,
) -> Result<[i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT], String> {
    custom_scores_from_bytes(contents.as_bytes(), mask_score).map_err(|error| error.to_string())
}

/// `stats/score_matrix.cpp::custom_scores`: getline preserves every byte,
/// comment recognition occurs before alphabet conversion, and the row
/// stringstream retains its read/write cursors between appended lines.
pub fn custom_scores_from_bytes(
    contents: &[u8],
    mask_score: i32,
) -> std::io::Result<[i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT]> {
    let alphabet = crate::basic::value::CharRepresentation::new(
        crate::basic::value::AMINO_ACID_ALPHABET,
        crate::basic::value::MASK_LETTER,
        b"UO-",
    );
    let mut scores = [mask_score as i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT];
    let mut positions = Vec::new();
    let mut row = 0usize;
    let mut stream = Vec::new();
    let mut cursor = 0usize;
    let mut min_score = i32::MAX;
    for line in contents.split(|&byte| byte == b'\n') {
        if line.first() == Some(&b'#') {
            continue;
        }
        if positions.is_empty() {
            for &byte in line {
                if byte != b' ' && byte != b'\t' {
                    positions.push(alphabet.convert(byte).map_err(std::io::Error::other)? as usize);
                }
            }
        } else {
            if row >= positions.len() {
                break;
            }
            stream.extend_from_slice(line);
            let first = *stream.get(cursor).ok_or_else(|| {
                std::io::Error::other(
                    "Source custom matrix row reads beyond its actual stream at end of input.",
                )
            })?;
            cursor += 1;
            if alphabet.convert(first).map_err(std::io::Error::other)? as usize != positions[row] {
                return Err(std::io::Error::other(
                    "Invalid custom scoring matrix file format.",
                ));
            }
            let mut failed = false;
            for &column in &positions {
                // A failed sentry leaves Source's fresh automatic int
                // indeterminate. A num_get failure itself assigns zero (or
                // the saturated limit), then leaves failbit for later cells.
                if failed {
                    return Err(std::io::Error::other("Source custom matrix score extraction leaves an indeterminate value after formatted-input failure."));
                }
                while stream.get(cursor).is_some_and(u8::is_ascii_whitespace) {
                    cursor += 1;
                }
                if cursor == stream.len() {
                    return Err(std::io::Error::other("Source custom matrix score extraction leaves an indeterminate value after formatted-input failure."));
                }
                let negative = stream.get(cursor) == Some(&b'-');
                if matches!(stream.get(cursor), Some(b'+' | b'-')) {
                    cursor += 1;
                }
                let digits = cursor;
                let mut magnitude = 0u64;
                while let Some(digit) = stream.get(cursor).filter(|byte| byte.is_ascii_digit()) {
                    magnitude = magnitude
                        .saturating_mul(10)
                        .saturating_add(u64::from(*digit - b'0'));
                    cursor += 1;
                }
                let limit = if negative {
                    (i32::MAX as u64) + 1
                } else {
                    i32::MAX as u64
                };
                failed = cursor == digits || magnitude > limit;
                let score = if cursor == digits {
                    0
                } else if negative {
                    -(magnitude.min(limit) as i64) as i32
                } else {
                    magnitude.min(limit) as i32
                };
                scores[positions[row] * AMINO_ACID_COUNT + column] = score as i8;
                min_score = min_score.min(score);
            }
            row += 1;
        }
    }
    let hard_mask = crate::basic::value::SUPER_HARD_MASK as usize;
    for i in 0..AMINO_ACID_COUNT {
        scores[i * AMINO_ACID_COUNT + hard_mask] = min_score as i8;
        scores[hard_mask * AMINO_ACID_COUNT + i] = min_score as i8;
    }
    Ok(scores)
}

fn area_params_from_alp(p: &super::pvalues::ALP_set_of_parameters) -> super::pvalues::AreaParams {
    super::pvalues::AreaParams {
        a_i: p.a_I,
        b_i: p.b_I,
        alpha_i: p.alpha_I,
        beta_i: p.beta_I,
        a_j: p.a_J,
        b_j: p.b_J,
        alpha_j: p.alpha_J,
        beta_j: p.beta_J,
        sigma: p.sigma,
        tau: p.tau,
        vi_y_thr: p.vi_y_thr,
        vj_y_thr: p.vj_y_thr,
        c_y_thr: p.c_y_thr,
    }
}

/// Runtime score matrix used for alignment.
///
/// Holds the scoring data in various representations optimized for
/// different SIMD widths (8-bit, 16-bit, 32-bit), plus statistical
/// parameters for E-value computation.
// C++ `Scores<T>::data` is `alignas(32)`.  Keep the tables in declaration
// order and align the containing object so each 32-wide row has the same
// alignment; the preceding table sizes are all multiples of 32 bytes.
#[repr(C, align(32))]
#[derive(Clone)]
pub struct ScoreMatrix {
    /// 32x32 score matrix (padded for alignment), 32-bit scores.
    matrix32: [i32; 32 * 32],
    /// 32x32 score matrix, 8-bit scores.
    matrix8: [i8; 32 * 32],
    /// 32x32 unsigned 8-bit scores (biased).
    matrix8u: [u8; 32 * 32],
    matrix8_low: [i8; 32 * 32],
    matrix8_high: [i8; 32 * 32],
    matrix8u_low: [i8; 32 * 32],
    matrix8u_high: [i8; 32 * 32],
    matrix16: [i16; 32 * 32],
    matrix32_scaled: [i32; 32 * 32],
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
    standard_matrix: Option<&'static StandardMatrix>,
    /// Expected score for each true amino acid against BLOSUM62 background.
    background_scores: [f64; TRUE_AA as usize],
    /// ALP area parameters for E-value computation.
    area_params: super::pvalues::AreaParams,
    dna_statistics: Option<crate::dna::build_score::Blastn_Score>,
    source_static_dna: bool,
}

impl ScoreMatrix {
    /// Reporting/scoring metadata for an already validated nucleotide model.
    ///
    /// This typed constructor explicitly supplies a usable DNA model. Source's
    /// static blastn global has different zero-initialized reporting metadata;
    /// native CLI uses source_blastn_metadata instead. Neither constructor
    /// changes new("DNA") or chooses replacement gap defaults.
    pub fn from_dna_statistics(
        stats: crate::dna::build_score::Blastn_Score,
    ) -> Result<Self, String> {
        if !stats.lambda().is_finite()
            || stats.lambda() <= 0.0
            || !stats.k().is_finite()
            || stats.k() <= 0.0
            || !stats.h().is_finite()
            || stats.h() <= 0.0
        {
            return Err("Invalid nucleotide Karlin statistics".to_owned());
        }
        // KSW explicitly narrows its substitution bytes; WFA and nucleotide
        // statistics retain full i32 costs. A reporting carrier must not reject
        // valid gcd-scaled Karlin/WFA configurations solely for byte tables.
        let reward = stats.reward() as i8;
        let penalty = stats.penalty() as i8;
        let mut scores = [0i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT];
        for i in 0..4 {
            for j in 0..4 {
                scores[i * AMINO_ACID_COUNT + j] = if i == j { reward } else { penalty };
            }
        }
        let matrix8 = build_i8_scores(&scores, 1, 0, 32, 0);
        let bias = penalty.wrapping_neg().max(0);
        let matrix8u = build_i8_scores(&scores, 1, bias, 32, 0).map(|v| v as u8);
        let mut matrix16 = build_i16_scores(&scores, 1);
        let mut matrix32 = matrix16.map(i32::from);
        for i in 0..4 {
            for j in 0..4 {
                let score = if i == j {
                    stats.reward()
                } else {
                    stats.penalty()
                };
                matrix16[i * 32 + j] = score as i16;
                matrix32[i * 32 + j] = score;
            }
        }
        Ok(Self {
            matrix32,
            matrix8,
            matrix8u,
            matrix8_low: build_i8_scores(&scores, 1, 0, 16, 0),
            matrix8_high: build_i8_scores(&scores, 1, 0, 16, 16),
            matrix8u_low: build_i8_scores(&scores, 1, bias, 16, 0),
            matrix8u_high: build_i8_scores(&scores, 1, bias, 16, 16),
            matrix16,
            matrix32_scaled: matrix32,
            bias,
            gap_open: stats.gap_open(),
            gap_extend: stats.gap_extend(),
            frame_shift: 0,
            db_letters: stats.database_letters() as f64,
            ln_k: stats.k().ln(),
            scale: 1.0,
            lambda: stats.lambda(),
            k: stats.k(),
            name: "DNA".to_owned(),
            standard_matrix: None,
            background_scores: [0.0; TRUE_AA as usize],
            // No ALP model is invented for DNA. DNA entry points below bypass it.
            area_params: super::pvalues::AreaParams {
                a_i: 0.0,
                b_i: 0.0,
                alpha_i: 0.0,
                beta_i: 0.0,
                a_j: 0.0,
                b_j: 0.0,
                alpha_j: 0.0,
                beta_j: 0.0,
                sigma: 0.0,
                tau: 0.0,
                vi_y_thr: 0.0,
                vj_y_thr: 0.0,
                c_y_thr: 0.0,
            },
            dna_statistics: Some(stats),
            source_static_dna: false,
        })
    }

    /// Exact defined state of Source's static global for blastn. Static zero
    /// initialization precedes ScoreMatrix's ln_k=0 and ALP's flag=false
    /// constructors. Keep the real DNA HSP model solely for the separate
    /// WITH_DNA IntermediateRecord constructor, not global metadata or scores.
    pub(crate) fn source_blastn_metadata(
        stats: crate::dna::build_score::Blastn_Score,
    ) -> Result<Self, String> {
        let mut matrix = Self::from_dna_statistics(stats)?;
        matrix.matrix32.fill(0);
        matrix.matrix8.fill(0);
        matrix.matrix8u.fill(0);
        matrix.matrix8_low.fill(0);
        matrix.matrix8_high.fill(0);
        matrix.matrix8u_low.fill(0);
        matrix.matrix8u_high.fill(0);
        matrix.matrix16.fill(0);
        matrix.matrix32_scaled.fill(0);
        matrix.bias = 0;
        matrix.gap_open = 0;
        matrix.gap_extend = 0;
        matrix.frame_shift = 0;
        matrix.db_letters = 0.0;
        matrix.ln_k = 0.0;
        matrix.scale = 0.0;
        matrix.lambda = 0.0;
        matrix.k = 0.0;
        matrix.name.clear();
        matrix.source_static_dna = true;
        Ok(matrix)
    }

    pub(crate) fn source_dna_statistics(&self) -> Option<crate::dna::build_score::Blastn_Score> {
        self.source_static_dna
            .then_some(self.dna_statistics)
            .flatten()
    }

    fn source_evaluator_error(&self, query_len: u32, subject_len: u32) -> Option<&'static str> {
        if !self.source_static_dna {
            return None;
        }
        // Source checks lengths before testing the static ALP model's flag.
        Some(if query_len == 0 || subject_len == 0 {
            "Error - seqlen1_<=0 or seq2en1_<=0 in \"double AlignmentEvaluer::area\"\n"
        } else {
            "Unexpected error - the Gumbel parameters are not defined properly in \"double AlignmentEvaluer::area\"\n"
        })
    }

    /// Create a score matrix from a named matrix and gap penalties.
    pub fn new(
        matrix_name: &str,
        gap_open: i32,
        gap_extend: i32,
        frame_shift: i32,
        stop_match_score: i32,
        db_letters: u64,
    ) -> Result<Self, String> {
        Self::new_scaled(
            matrix_name,
            gap_open,
            gap_extend,
            frame_shift,
            stop_match_score,
            db_letters,
            1,
        )
    }

    /// Named-matrix constructor with the C++ `scale` argument exposed.
    #[allow(clippy::too_many_arguments)]
    pub fn new_scaled(
        matrix_name: &str,
        gap_open: i32,
        gap_extend: i32,
        frame_shift: i32,
        stop_match_score: i32,
        db_letters: u64,
        scale: i32,
    ) -> Result<Self, String> {
        if scale <= 0 {
            return Err("Score matrix scale must be positive.".to_string());
        }
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
        let matrix_statistics = super::matrix_data::get(standard_matrix);

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
        let min_score = matrix_low_score(&matrix32);
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

        let matrix8_low = build_i8_scores(&standard_matrix.scores, stop_match_score, 0, 16, 0);
        let matrix8_high = build_i8_scores(&standard_matrix.scores, stop_match_score, 0, 16, 16);
        let matrix8u_low = build_i8_scores(&standard_matrix.scores, stop_match_score, bias, 16, 0);
        let matrix8u_high =
            build_i8_scores(&standard_matrix.scores, stop_match_score, bias, 16, 16);
        let matrix16 = build_i16_scores(&standard_matrix.scores, stop_match_score);
        let mut matrix32_scaled = [-128i32; 32 * 32];
        for i in 0..32 {
            for j in 0..32 {
                if i < AMINO_ACID_COUNT && j < AMINO_ACID_COUNT {
                    matrix32_scaled[i * 32 + j] = if i < TRUE_AA as usize && j < TRUE_AA as usize {
                        let ni = super::cbs::ALPH_TO_NCBI[i];
                        let nj = super::cbs::ALPH_TO_NCBI[j];
                        (matrix_statistics.freq_ratios[ni][nj].ln()
                            / standard_matrix.ungapped_constants().lambda
                            * scale as f64)
                            .round() as i32
                    } else {
                        matrix32[i * 32 + j] * scale
                    };
                }
            }
        }

        let ln_k = params.k.ln();
        let mut background_scores = [0.0f64; TRUE_AA as usize];
        for i in 0..TRUE_AA as usize {
            for j in 0..TRUE_AA as usize {
                background_scores[i] += BLOSUM62_BACKGROUND_FREQS[j] * matrix32[i * 32 + j] as f64;
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
            matrix8_low,
            matrix8_high,
            matrix8u_low,
            matrix8u_high,
            matrix16,
            matrix32_scaled,
            bias,
            gap_open: go,
            gap_extend: ge,
            frame_shift,
            db_letters: db_letters as f64,
            ln_k,
            scale: scale as f64,
            lambda: params.lambda,
            k: params.k,
            name: matrix_name.to_lowercase(),
            standard_matrix: Some(standard_matrix),
            background_scores,
            area_params,
            dna_statistics: None,
            source_static_dna: false,
        })
    }

    /// Construct a score matrix from the custom matrix-file format.
    pub fn from_custom_file(
        matrix_file: impl AsRef<Path>,
        gap_open: i32,
        gap_extend: i32,
        stop_match_score: i32,
        db_letters: u64,
    ) -> Result<Self, String> {
        Self::from_custom_file_path_io(
            matrix_file.as_ref(),
            gap_open,
            gap_extend,
            stop_match_score,
            db_letters,
        )
        .map_err(|error| error.to_string())
    }

    /// Source-valid files are byte streams, not UTF-8 documents. The Source
    /// unchecked ifstream loops forever after a failed open; this fallible
    /// native boundary is an explicitly bounded extension for that domain.
    pub fn from_custom_file_path_io(
        matrix_file: &Path,
        gap_open: i32,
        gap_extend: i32,
        stop_match_score: i32,
        db_letters: u64,
    ) -> std::io::Result<Self> {
        if matrix_file.as_os_str().is_empty() {
            return Self::from_custom_scores(
                [-gap_extend as i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT],
                gap_open,
                gap_extend,
                stop_match_score,
                db_letters,
            )
            .map_err(std::io::Error::other);
        }
        #[cfg(not(windows))]
        let bytes = fs::read(crate::util::io::source_c_path(matrix_file).as_ref());
        #[cfg(windows)]
        let bytes = (|| {
            // Source ifstream uses CRT text mode on Windows. Read through the
            // maintained CRT owner rather than interpreting UTF-8 or manually
            // rewriting CR/LF (which would miss the CRT's other text rules).
            let mut file =
                crate::util::io::cstdio_file::OwnedCstdio::open_file_role_path(matrix_file, "r")
                    .map_err(std::io::Error::other)?;
            let mut bytes = Vec::new();
            let mut chunk = [0u8; 8192];
            loop {
                let count = file.read_bytes(&mut chunk).map_err(std::io::Error::other)?;
                bytes.extend_from_slice(&chunk[..count]);
                if count < chunk.len() {
                    break;
                }
            }
            file.close_file_role().map_err(std::io::Error::other)?;
            Ok::<_, std::io::Error>(bytes)
        })();
        let contents = bytes.map_err(|error| {
            // Preserve the existing native bounded error, retaining raw path
            // bytes rather than claiming Source's nontermination as parity.
            let path = match crate::util::io::path_diagnostic_bytes(matrix_file) {
                Ok(path) => path,
                Err(error) => return std::io::Error::other(error),
            };
            let mut message = b"Failed to read custom scoring matrix '".to_vec();
            message.extend_from_slice(&path);
            message.extend_from_slice(b"': ");
            message.extend_from_slice(error.to_string().as_bytes());
            crate::util::io::RawDiagnosticError(message).into_io()
        })?;
        Self::from_custom_bytes(
            &contents,
            gap_open,
            gap_extend,
            stop_match_score,
            db_letters,
        )
    }

    /// In-memory counterpart of the custom-file constructor.
    pub fn from_custom_str(
        contents: &str,
        gap_open: i32,
        gap_extend: i32,
        stop_match_score: i32,
        db_letters: u64,
    ) -> Result<Self, String> {
        Self::from_custom_bytes(
            contents.as_bytes(),
            gap_open,
            gap_extend,
            stop_match_score,
            db_letters,
        )
        .map_err(|error| error.to_string())
    }
    pub fn from_custom_bytes(
        contents: &[u8],
        gap_open: i32,
        gap_extend: i32,
        stop_match_score: i32,
        db_letters: u64,
    ) -> std::io::Result<Self> {
        let scores = custom_scores_from_bytes(contents, -gap_extend)?;
        Self::from_custom_scores(scores, gap_open, gap_extend, stop_match_score, db_letters)
            .map_err(std::io::Error::other)
    }
    fn from_custom_scores(
        score_array: [i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT],
        gap_open: i32,
        gap_extend: i32,
        stop_match_score: i32,
        db_letters: u64,
    ) -> Result<Self, String> {
        let mut matrix32 = [-128i32; 32 * 32];
        let mut matrix8 = [i8::MIN; 32 * 32];
        for i in 0..AMINO_ACID_COUNT {
            for j in 0..AMINO_ACID_COUNT {
                let score = score_array[i * AMINO_ACID_COUNT + j];
                matrix8[i * 32 + j] = score;
                matrix32[i * 32 + j] = score as i32;
            }
        }
        if stop_match_score != 1 {
            matrix8[24 * 32 + 24] = stop_match_score as i8;
            matrix32[24 * 32 + 24] = stop_match_score;
        }
        let bias = -matrix_low_score(&matrix32);
        let mut matrix8u = [128u8; 32 * 32];
        for i in 0..AMINO_ACID_COUNT {
            for j in 0..AMINO_ACID_COUNT {
                matrix8u[i * 32 + j] = (matrix8[i * 32 + j] as i16 + bias as i16) as u8;
            }
        }
        let matrix8_low = build_i8_scores(&score_array, stop_match_score, 0, 16, 0);
        let matrix8_high = build_i8_scores(&score_array, stop_match_score, 0, 16, 16);
        let matrix8u_low = build_i8_scores(&score_array, stop_match_score, bias, 16, 0);
        let matrix8u_high = build_i8_scores(&score_array, stop_match_score, bias, 16, 16);
        let matrix16 = build_i16_scores(&score_array, stop_match_score);
        let matrix32_scaled = matrix32;

        let substitution_scores = (0..TRUE_AA as usize)
            .map(|i| {
                (0..TRUE_AA as usize)
                    .map(|j| score_array[i * AMINO_ACID_COUNT + j] as i64)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut evaluer = super::sls_alignment_evaluer::AlignmentEvaluer::new();
        evaluer
            .initGapped(
                TRUE_AA as i64,
                &substitution_scores,
                &BLOSUM62_BACKGROUND_FREQS,
                &BLOSUM62_BACKGROUND_FREQS,
                gap_open as i64,
                gap_extend as i64,
                gap_open as i64,
                gap_extend as i64,
                false,
                0.01,
                0.05,
                120.0,
                1024.0,
                1,
                1.07,
            )
            .map_err(|_| {
                "The ALP library failed to compute the statistical parameters for this matrix. It may help to adjust the gap penalty settings."
                    .to_string()
            })?;
        let params = evaluer.parameters();
        let mut background_scores = [0.0; TRUE_AA as usize];
        for i in 0..TRUE_AA as usize {
            for j in 0..TRUE_AA as usize {
                background_scores[i] += BLOSUM62_BACKGROUND_FREQS[j] * matrix32[i * 32 + j] as f64;
            }
        }
        Ok(Self {
            matrix32,
            matrix8,
            matrix8u,
            matrix8_low,
            matrix8_high,
            matrix8u_low,
            matrix8u_high,
            matrix16,
            matrix32_scaled,
            bias,
            gap_open,
            gap_extend,
            frame_shift: 0,
            db_letters: db_letters as f64,
            ln_k: params.K.ln(),
            scale: 1.0,
            lambda: params.lambda,
            k: params.K,
            name: "custom".to_string(),
            standard_matrix: None,
            background_scores,
            area_params: area_params_from_alp(params),
            dna_statistics: None,
            source_static_dna: false,
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
        if let Some(message) = self.source_evaluator_error(query_len, subject_len) {
            panic!("{message}");
        }
        if let Some(stats) = self.dna_statistics {
            return stats.blast_bit_score(raw_score);
        }
        let log_area = super::pvalues::compute_log_area(
            &self.area_params,
            raw_score as f64,
            subject_len as f64,
            query_len as f64,
        );
        // AlignmentEvaluer::log_area passes subject/query lengths in that order.
        // The pinned -march=native oracle contracts lambda*score-ln(K) into
        // one FMA. Make its rounding explicit for binary edge records.
        (self.lambda.mul_add(raw_score as f64, -self.ln_k) - log_area) / LN_2
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
        if let Some(message) = self.source_evaluator_error(query_len, subject_len) {
            panic!("{message}");
        }
        if let Some(stats) = self.dna_statistics {
            return stats.blast_e_value(raw_score, query_len as i32);
        }
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
        if let Some(message) = self.source_evaluator_error(query_len, subject_len) {
            panic!("{message}");
        }
        if self.dna_statistics.is_some() {
            return self.evalue_norm_query_len(raw_score, query_len as i32);
        }
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

    /// Update the global reporting/scoring carrier. Source keeps its separate
    /// Blastn_Score model at full database size even with --dbsize/filtering;
    /// therefore this deliberately does not rewrite DNA model statistics.
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

    pub fn matrix8_low(&self) -> &[i8; 32 * 32] {
        &self.matrix8_low
    }

    pub fn matrix8_high(&self) -> &[i8; 32 * 32] {
        &self.matrix8_high
    }

    pub fn matrix8u(&self) -> &[u8; 32 * 32] {
        &self.matrix8u
    }

    pub fn matrix8u_low(&self) -> &[i8; 32 * 32] {
        &self.matrix8u_low
    }

    pub fn matrix8u_high(&self) -> &[i8; 32 * 32] {
        &self.matrix8u_high
    }

    pub fn matrix16(&self) -> &[i16; 32 * 32] {
        &self.matrix16
    }

    /// Get the raw 32-bit score matrix data.
    pub fn matrix32(&self) -> &[i32; 32 * 32] {
        &self.matrix32
    }

    pub fn matrix32_scaled(&self) -> &[i32; 32 * 32] {
        &self.matrix32_scaled
    }

    pub fn matrix32_scaled_rows(&self) -> Vec<&[i32]> {
        self.matrix32_scaled.chunks_exact(32).collect()
    }

    pub fn ungapped_lambda(&self) -> Option<f64> {
        self.standard_matrix
            .map(|matrix| matrix.ungapped_constants().lambda)
    }

    pub fn standard_matrix(&self) -> Option<&'static StandardMatrix> {
        self.standard_matrix
    }

    pub fn joint_probs(&self) -> Option<&'static [f64; 20 * 20]> {
        self.standard_matrix
            .map(|matrix| super::matrix_data::get(matrix).joint_probs)
    }

    pub fn background_freqs(&self) -> &'static [f64; 20] {
        self.standard_matrix
            .map(|matrix| super::matrix_data::get(matrix).background_freqs)
            .unwrap_or(&BLOSUM62_BACKGROUND_FREQS)
    }

    pub fn freq_ratios(&self) -> Option<&'static [[f64; 28]; 28]> {
        self.standard_matrix
            .map(|matrix| super::matrix_data::get(matrix).freq_ratios)
    }

    pub fn background_scores(&self) -> &[f64; TRUE_AA as usize] {
        &self.background_scores
    }
}

impl fmt::Display for ScoreMatrix {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "(Matrix={} Lambda={} K={} Penalties={}/{})",
            self.name, self.lambda, self.k, self.gap_open, self.gap_extend
        )
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
        assert_eq!(sm.matrix8_low()[0], 4);
        assert_eq!(sm.matrix8_high()[0], 0);
        assert_eq!(sm.matrix8u_low()[0], 16);
        assert_eq!(sm.matrix16()[0], 4);
        assert_eq!(sm.matrix32_scaled()[0], 4);
        assert_eq!(sm.matrix32_scaled_rows().len(), 32);
        assert_eq!(sm.ungapped_lambda(), Some(0.3176));
        assert_eq!(
            sm.to_string(),
            "(Matrix=blosum62 Lambda=0.267 K=0.041 Penalties=11/1)"
        );
    }

    #[test]
    fn scaled_matrix_and_score_views_match_cpp_layout() {
        let sm = ScoreMatrix::new_scaled("blosum62", 11, 1, 0, 1, 0, 4).unwrap();
        assert_eq!(sm.matrix32_scaled()[0], 17);
        assert_eq!(sm.matrix32_scaled()[31], -128);
        assert_eq!(sm.matrix8_high()[0], 0); // A-T: high half starts at index 16
        assert_eq!(sm.matrix8u_high()[0], 12); // A-T (0) plus BLOSUM62 bias 12
        assert_eq!(sm.matrix8u()[31], 128);
        assert_eq!(sm.joint_probs().unwrap().len(), 400);
        assert_eq!(sm.background_freqs(), &BLOSUM62_BACKGROUND_FREQS);
        assert!((sm.freq_ratios().unwrap()[1][1] - 3.90294070).abs() < 1e-8);
    }

    #[test]
    fn simd_score_tables_are_32_byte_aligned() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        assert_eq!(sm.matrix8u_low().as_ptr() as usize % 32, 0);
        assert_eq!(sm.matrix8u_high().as_ptr() as usize % 32, 0);
        assert_eq!(sm.matrix16().as_ptr() as usize % 32, 0);
    }

    #[test]
    fn custom_matrix_byte_comments_compact_headers_and_source_aliases() {
        let scores =
            custom_scores_from_bytes(b"# raw comment \xff\x80\nAR\nA5 -2\nR-3 7\n", -1).unwrap();
        assert_eq!(
            (
                scores[0],
                scores[1],
                scores[AMINO_ACID_COUNT],
                scores[AMINO_ACID_COUNT + 1]
            ),
            (5, -2, -3, 7)
        );
        let aliases = custom_scores_from_bytes(b"AU\nA5 -2\nU-3 7", -1).unwrap();
        let mask = crate::basic::value::MASK_LETTER as usize;
        assert_eq!(aliases[mask], -2);
        assert_eq!(aliases[mask * AMINO_ACID_COUNT], -3);
        assert_eq!(aliases[mask * AMINO_ACID_COUNT + mask], 7);
        // The FIRST formatted extraction failure assigns 0; no following
        // extraction observes failbit in this one-column, one-row matrix.
        assert_eq!(custom_scores_from_bytes(b"A\nA X", -1).unwrap()[0], 0);
        assert_eq!(
            custom_scores_from_bytes(b"A\nA 2147483648", -1).unwrap()[0],
            -1
        );
    }

    #[test]
    fn custom_score_parser_maps_header_rows_and_hard_mask() {
        let scores = custom_scores_from_str("# compact matrix\nA R\nA 5 -2\nR -3 7\n", -1).unwrap();
        assert_eq!(scores[0], 5);
        assert_eq!(scores[1], -2);
        assert_eq!(scores[AMINO_ACID_COUNT], -3);
        assert_eq!(scores[AMINO_ACID_COUNT + 1], 7);
        let hard_mask = crate::basic::value::SUPER_HARD_MASK as usize;
        assert_eq!(scores[hard_mask], -3);
        assert_eq!(scores[hard_mask * AMINO_ACID_COUNT], -3);
        assert!(custom_scores_from_str("A R\nR 1 2\n", -1).is_err());
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
