use crate::basic::reduction::Reduction;
use crate::basic::seed::SeedOffset;
use crate::basic::shape_config::ShapeConfig;
use crate::config::Sensitivity;
use crate::masking::MaskingAlgo;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::enum_utils::FlagBits;
use crate::util::math::{bit_length, power};

pub const SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE: f64 = 0.15;

pub const DEFAULT_CONTIGUOUS_REDUCTION: &str = "KR EQ D N C G H F Y IV LM W P S T A";

/// Matches C++ `Extension::Mode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtensionMode {
    BandedFast,
    BandedSlow,
    Full,
    None,
    Global,
}

impl ExtensionMode {
    /// Matches C++ `EnumTraits<Extension::Mode>::from_string`.
    pub fn from_string(s: &str) -> Result<Self, String> {
        match s {
            "banded-fast" => Ok(Self::BandedFast),
            "banded-slow" => Ok(Self::BandedSlow),
            "full" => Ok(Self::Full),
            "none" => Ok(Self::None),
            "global" => Ok(Self::Global),
            _ => Err(format!(
                "Invalid value for string field: {}. Permitted values: banded-fast, banded-slow, full, none, global",
                s
            )),
        }
    }
}

/// Matches C++ `Extension::default_ext_mode`.
pub fn default_ext_mode(sens: Sensitivity) -> ExtensionMode {
    match sens {
        Sensitivity::Faster
        | Sensitivity::Fast
        | Sensitivity::Shapes6x10
        | Sensitivity::Shapes30x10
        | Sensitivity::Linclust20
        | Sensitivity::Linclust40
        | Sensitivity::Default
        | Sensitivity::MidSensitive
        | Sensitivity::Sensitive => ExtensionMode::BandedFast,
        Sensitivity::MoreSensitive | Sensitivity::VerySensitive | Sensitivity::UltraSensitive => {
            ExtensionMode::BandedSlow
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SetupSearchInput {
    pub freq_sd: Option<f64>,
    pub min_identities: Option<u32>,
    pub ungapped_evalue: Option<f64>,
    pub ungapped_evalue_short: Option<f64>,
    pub gapped_filter_evalue: Option<f64>,
    pub query_bins: Option<u32>,
    pub threads: i32,
    pub minimizer_window: Option<i32>,
    pub sketch_size: Option<i32>,
    pub contiguous_seed_mode: bool,
    pub shape_mask: Vec<String>,
    pub shapes: u32,
    pub lin_stage1_query: bool,
    pub lin_stage1_target: bool,
    pub lin_stage1_combo: bool,
    pub gapped_filter_diag_bit_score: f64,
    pub seed_cut: f64,
    pub motif_masking: String,
    pub swipe_all: bool,
    pub freq_masking: bool,
    pub soft_masking: String,
    pub ext: String,
    pub global_ranking_targets: bool,
    pub linclust_banded_ext: bool,
    pub frame_shift: i32,
    pub anchored_swipe: bool,
    pub aln_out: String,
    pub parallel_tmpdir: String,
    pub approx_min_id: f64,
}

impl Default for SetupSearchInput {
    fn default() -> Self {
        Self {
            freq_sd: None,
            min_identities: None,
            ungapped_evalue: None,
            ungapped_evalue_short: None,
            gapped_filter_evalue: None,
            query_bins: None,
            threads: 1,
            minimizer_window: None,
            sketch_size: None,
            contiguous_seed_mode: false,
            shape_mask: Vec::new(),
            shapes: 0,
            lin_stage1_query: false,
            lin_stage1_target: false,
            lin_stage1_combo: false,
            gapped_filter_diag_bit_score: 0.0,
            seed_cut: 0.0,
            motif_masking: String::new(),
            swipe_all: false,
            freq_masking: false,
            soft_masking: String::new(),
            ext: String::new(),
            global_ranking_targets: false,
            linclust_banded_ext: false,
            frame_shift: 0,
            anchored_swipe: false,
            aln_out: String::new(),
            parallel_tmpdir: String::new(),
            approx_min_id: 0.0,
        }
    }
}

pub struct SetupSearchResult {
    pub sensitivity: Sensitivity,
    pub freq_sd: f64,
    pub hamming_filter_id: u32,
    pub ungapped_evalue: f64,
    pub ungapped_evalue_short: f64,
    pub gapped_filter_evalue: f64,
    pub query_bins: u32,
    pub minimizer_window: i32,
    pub sketch_size: i32,
    pub reduction: Reduction,
    pub shapes: ShapeConfig,
    pub gapped_filter_diag_score: i32,
    pub seed_complexity_cut: f64,
    pub soft_masking_bits: u32,
    pub extension_mode: ExtensionMode,
}

/// Matches C++ `Round`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Round {
    pub sensitivity: Sensitivity,
    pub query_indexed: bool,
}

impl Round {
    pub const fn new(sensitivity: Sensitivity, query_indexed: bool) -> Self {
        Self {
            sensitivity,
            query_indexed,
        }
    }
}

pub const ITERATED_FASTER: &[Round] = &[];
pub const ITERATED_FAST: &[Round] = &[Round::new(Sensitivity::Fast, true)];
pub const ITERATED_DEFAULT: &[Round] = &[
    Round::new(Sensitivity::Fast, true),
    Round::new(Sensitivity::Linclust40, true),
];
pub const ITERATED_LINCLUST40: &[Round] = &[
    Round::new(Sensitivity::Fast, true),
    Round::new(Sensitivity::Linclust40, true),
];
pub const ITERATED_LINCLUST20: &[Round] = &[
    Round::new(Sensitivity::Fast, true),
    Round::new(Sensitivity::Linclust20, true),
];
pub const ITERATED_SHAPES30X10: &[Round] = &[
    Round::new(Sensitivity::Fast, true),
    Round::new(Sensitivity::Shapes30x10, true),
];
pub const ITERATED_SHAPES6X10: &[Round] = &[Round::new(Sensitivity::Shapes6x10, true)];
pub const ITERATED_MID_SENSITIVE: &[Round] = &[
    Round::new(Sensitivity::Fast, true),
    Round::new(Sensitivity::Linclust40, true),
    Round::new(Sensitivity::Default, false),
];
pub const ITERATED_SENSITIVE: &[Round] = &[
    Round::new(Sensitivity::Fast, true),
    Round::new(Sensitivity::Linclust20, true),
    Round::new(Sensitivity::Default, false),
];
pub const ITERATED_MORE_SENSITIVE: &[Round] = &[
    Round::new(Sensitivity::Fast, true),
    Round::new(Sensitivity::Linclust20, true),
    Round::new(Sensitivity::Default, false),
];
pub const ITERATED_VERY_ULTRA_SENSITIVE: &[Round] = &[
    Round::new(Sensitivity::Fast, true),
    Round::new(Sensitivity::Linclust20, true),
    Round::new(Sensitivity::Default, false),
    Round::new(Sensitivity::MoreSensitive, false),
];

/// Sensitivity-specific search parameters.
#[derive(Debug, Clone)]
pub struct SensitivityTraits {
    pub query_indexed: bool,
    pub motif_masking: bool,
    pub freq_sd: f64,
    pub min_identities: u32,
    pub ungapped_evalue: u64,
    pub ungapped_evalue_short: u64,
    pub gapped_filter_evalue: u64,
    pub index_chunks: u32,
    pub query_bins: u32,
    pub contiguous_seed: Option<&'static str>,
    pub seed_cut: f64,
    pub block_size: f64,
    pub minimizer_window: i32,
    pub sketch_size: i32,
}

/// Get sensitivity traits for a given level.
pub fn get_traits(sens: Sensitivity) -> SensitivityTraits {
    match sens {
        Sensitivity::Faster => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 50.0,
            min_identities: 11,
            ungapped_evalue: 0,
            ungapped_evalue_short: 0,
            gapped_filter_evalue: 0,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: None,
            seed_cut: 0.9,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 21,
        },
        Sensitivity::Fast => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 50.0,
            min_identities: 11,
            ungapped_evalue: 0,
            ungapped_evalue_short: 0,
            gapped_filter_evalue: 0,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: None,
            seed_cut: 0.9,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::Shapes6x10 => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 50.0,
            min_identities: 11,
            ungapped_evalue: 0,
            ungapped_evalue_short: 0,
            gapped_filter_evalue: 0,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: None,
            seed_cut: 0.9,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::Shapes30x10 => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 50.0,
            min_identities: 11,
            ungapped_evalue: 0,
            ungapped_evalue_short: 0,
            gapped_filter_evalue: 0,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: None,
            seed_cut: 0.9,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::Default => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 50.0,
            min_identities: 11,
            ungapped_evalue: 10000,
            ungapped_evalue_short: 10000,
            gapped_filter_evalue: 0,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: Some("111111"),
            seed_cut: 0.8,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::Linclust40 => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 50.0,
            min_identities: 11,
            ungapped_evalue: 0,
            ungapped_evalue_short: 0,
            gapped_filter_evalue: 0,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: None,
            seed_cut: 0.9,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::Linclust20 => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 50.0,
            min_identities: 11,
            ungapped_evalue: 0,
            ungapped_evalue_short: 0,
            gapped_filter_evalue: 0,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: None,
            seed_cut: 0.9,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::MidSensitive => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 20.0,
            min_identities: 11,
            ungapped_evalue: 10000,
            ungapped_evalue_short: 10000,
            gapped_filter_evalue: 0,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: None,
            seed_cut: 1.0,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::Sensitive => SensitivityTraits {
            query_indexed: true,
            motif_masking: true,
            freq_sd: 20.0,
            min_identities: 11,
            ungapped_evalue: 10000,
            ungapped_evalue_short: 10000,
            gapped_filter_evalue: 1,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: Some("11111"),
            seed_cut: 1.0,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::MoreSensitive => SensitivityTraits {
            query_indexed: true,
            motif_masking: false,
            freq_sd: 200.0,
            min_identities: 11,
            ungapped_evalue: 10000,
            ungapped_evalue_short: 10000,
            gapped_filter_evalue: 1,
            index_chunks: 4,
            query_bins: 16,
            contiguous_seed: Some("11111"),
            seed_cut: 1.0,
            block_size: 2.0,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::VerySensitive => SensitivityTraits {
            query_indexed: true,
            motif_masking: false,
            freq_sd: 15.0,
            min_identities: 9,
            ungapped_evalue: 100000,
            ungapped_evalue_short: 30000,
            gapped_filter_evalue: 1,
            index_chunks: 1,
            query_bins: 16,
            contiguous_seed: None,
            seed_cut: 1.0,
            block_size: 0.4,
            minimizer_window: 0,
            sketch_size: 0,
        },
        Sensitivity::UltraSensitive => SensitivityTraits {
            query_indexed: true,
            motif_masking: false,
            freq_sd: 20.0,
            min_identities: 9,
            ungapped_evalue: 300000,
            ungapped_evalue_short: 30000,
            gapped_filter_evalue: 1,
            index_chunks: 1,
            query_bins: 64,
            contiguous_seed: None,
            seed_cut: 1.0,
            block_size: 0.4,
            minimizer_window: 0,
            sketch_size: 0,
        },
    }
}

/// Shape codes for each sensitivity level.
pub fn get_shape_codes(sens: Sensitivity) -> &'static [&'static str] {
    match sens {
        Sensitivity::Faster | Sensitivity::Fast => &["1101110101101111"],
        Sensitivity::Default => &["111101110111", "111011010010111"],
        Sensitivity::Linclust20 => &[
            "111111111111",
            "1111111011111",
            "1111110111111",
            "11111111010111",
            "11011101111111",
            "11111011110111",
            "11110011111111",
            "11101111101111",
            "11110111111011",
            "110111110110111",
            "111101111011011",
            "111101100111111",
            "111010111110111",
            "111101011101111",
            "111110110011111",
            "111011101011111",
            "111111010011111",
            "111111001111011",
            "111110101101111",
            "111011110101111",
            "1110101110011111",
            "1111100110110111",
            "1110111001101111",
            "1111110010101111",
            "1111001010111111",
            "1110101101110111",
            "1110110111001111",
            "1110110101110111",
            "1111010101101111",
            "1111011011010111",
        ],
        Sensitivity::Linclust40 => &[
            "111111111111",
            "1111111011111",
            "1111110111111",
            "11111111010111",
            "11011101111111",
            "11111011110111",
            "11110011111111",
            "11101111101111",
            "11110111111011",
            "110111110110111",
            "111101111011011",
            "111101100111111",
            "111010111110111",
            "111101011101111",
            "111110110011111",
        ],
        Sensitivity::Shapes6x10 => &[
            "10111111111",
            "111110110111",
            "1101110111011",
            "111111101011",
            "1111011110011",
            "111111100100011",
        ],
        Sensitivity::Shapes30x10 => &[
            "10111111111",
            "111110110111",
            "1101110111011",
            "111111101011",
            "1111011110011",
            "111111100100011",
            "110111010011011",
            "1111100110010011",
            "11101100111101",
            "111011011010101",
            "11011010101111",
            "11111110000010011",
            "11011001100110011",
            "101011100011111",
            "111011111101",
            "111110101100101",
            "1111010101001011",
            "11100111011001001",
            "1110110001111001",
            "110111011000010011",
            "11001100101100111",
            "11111000000111101",
            "11011110011010001",
            "110101101010011001",
            "111010111000010101",
            "1111101000100010011",
            "11010100100111011",
            "101001111100111",
            "101110010001010111",
            "11001101001011011",
        ],
        Sensitivity::MidSensitive => &[
            "11110110111",
            "1101100111101",
            "1110010101111",
            "11010101100111",
            "11101110001011",
            "1110100100010111",
            "1101000011010111",
            "1110011000011011",
        ],
        Sensitivity::Sensitive | Sensitivity::MoreSensitive => &[
            "1011110111",
            "110100100010111",
            "11001011111",
            "101110001111",
            "11011101100001",
            "1111010010101",
            "111001001001011",
            "10101001101011",
            "111101010011",
            "1111000010000111",
            "1100011011011",
            "1101010000011011",
            "1110001010101001",
            "110011000110011",
            "11011010001101",
            "1101001100010011",
        ],
        Sensitivity::VerySensitive => &[
            "11101111",
            "110110111",
            "111111001",
            "1010111011",
            "11110001011",
            "110100101011",
            "110110001101",
            "1010101000111",
            "1100101001011",
            "1101010101001",
            "1110010010011",
            "110110000010011",
            "111001000100011",
            "1101000100010011",
        ],
        Sensitivity::UltraSensitive => &[
            "1111111",
            "11101111",
            "110011111",
            "110110111",
            "111111001",
            "1010111011",
            "1011110101",
            "1111000111",
            "10011110011",
            "10101101101",
            "10111010101",
            "11001010111",
            "11001100111",
            "11010101101",
            "11110001011",
            "100111010011",
            "101100110101",
            "101110000111",
            "110100101011",
            "110110001101",
            "111000110011",
            "1010001011011",
            "1010101000111",
            "1010110100011",
            "1100100110011",
            "1100101001011",
            "1101001100101",
            "1101010101001",
            "1110001010101",
            "1110010010011",
            "10100001101101",
            "11000100010111",
            "11010000100111",
            "11010100110001",
            "11101000011001",
            "11110000001101",
            "11110100000011",
            "101001000001111",
            "110000100101011",
            "110010010000111",
            "110101100001001",
            "110110000010011",
            "111001000100011",
            "111100000100101",
            "1000110010010101",
            "1001000100101101",
            "1001000110011001",
            "1010001001001011",
            "1010001010010011",
            "1010010001010101",
            "1010010100010011",
            "1010010101001001",
            "1010100000101011",
            "1010100011000101",
            "1011000010001011",
            "1100010000111001",
            "1100010010001011",
            "1100100001001011",
            "1100100100100011",
            "1100110000001101",
            "1101000100010011",
            "1101000110000101",
            "1110000001010011",
            "1110100000010101",
        ],
    }
}

/// Matches C++ `approx_id_to_hamming_id`.
pub fn approx_id_to_hamming_id() -> &'static [(f64, u32)] {
    &[(50.0, 20), (90.0, 30)]
}

/// Matches C++ `hamming_id_cutoff`.
pub fn hamming_id_cutoff(approx_id: f64) -> u32 {
    let mut cutoff = 0;
    for &(threshold, value) in approx_id_to_hamming_id() {
        if threshold > approx_id {
            break;
        }
        cutoff = value;
    }
    cutoff
}

/// Matches C++ `iterated_sens`.
pub fn iterated_sens(sens: Sensitivity) -> &'static [Round] {
    match sens {
        Sensitivity::Faster => ITERATED_FASTER,
        Sensitivity::Fast => ITERATED_FAST,
        Sensitivity::Default => ITERATED_DEFAULT,
        Sensitivity::Linclust40 => ITERATED_LINCLUST40,
        Sensitivity::Linclust20 => ITERATED_LINCLUST20,
        Sensitivity::Shapes30x10 => ITERATED_SHAPES30X10,
        Sensitivity::Shapes6x10 => ITERATED_SHAPES6X10,
        Sensitivity::MidSensitive => ITERATED_MID_SENSITIVE,
        Sensitivity::Sensitive => ITERATED_SENSITIVE,
        Sensitivity::MoreSensitive => ITERATED_MORE_SENSITIVE,
        Sensitivity::VerySensitive | Sensitivity::UltraSensitive => ITERATED_VERY_ULTRA_SENSITIVE,
    }
}

/// Matches C++ `Search::soft_masking_algo`.
pub fn soft_masking_algo(
    traits: &SensitivityTraits,
    motif_masking: &str,
    swipe_all: bool,
    freq_masking: bool,
) -> Result<MaskingAlgo, String> {
    if motif_masking.is_empty() {
        Ok(if !swipe_all && !freq_masking && traits.motif_masking {
            MaskingAlgo::Motif
        } else {
            MaskingAlgo::None
        })
    } else if motif_masking == "0" {
        Ok(MaskingAlgo::None)
    } else if motif_masking == "1" {
        if swipe_all {
            Err("Soft masking is not supported for --swipe.".to_string())
        } else {
            Ok(MaskingAlgo::Motif)
        }
    } else {
        Err("Permitted values for --motif-masking: 0, 1".to_string())
    }
}

/// Matches C++ `Search::setup_search`.
pub fn setup_search(
    sens: Sensitivity,
    input: &SetupSearchInput,
    score_matrix: &ScoreMatrix,
) -> Result<SetupSearchResult, String> {
    let traits = get_traits(sens);
    let freq_sd = input.freq_sd.unwrap_or(traits.freq_sd);
    let hamming_filter_id = input
        .min_identities
        .unwrap_or(traits.min_identities)
        .max(hamming_id_cutoff(input.approx_min_id));
    let ungapped_evalue = input
        .ungapped_evalue
        .unwrap_or(traits.ungapped_evalue as f64);
    let ungapped_evalue_short = input
        .ungapped_evalue_short
        .unwrap_or(traits.ungapped_evalue_short as f64);
    let gapped_filter_evalue = input
        .gapped_filter_evalue
        .unwrap_or(traits.gapped_filter_evalue as f64);
    let query_bins = input
        .query_bins
        .unwrap_or(((input.threads as f64 / 8.0).round() as u32).max(traits.query_bins));
    let minimizer_window = input.minimizer_window.unwrap_or(traits.minimizer_window);
    let sketch_size = input.sketch_size.unwrap_or(traits.sketch_size);

    let reduction = if input.contiguous_seed_mode && sens == Sensitivity::Default {
        Reduction::new(
            DEFAULT_CONTIGUOUS_REDUCTION,
            crate::basic::value::AMINO_ACID_ALPHABET,
        )
    } else {
        Reduction::default_reduction()
    };
    let shape_codes: Vec<String> = if input.contiguous_seed_mode {
        let contiguous_seed = traits.contiguous_seed.ok_or_else(|| {
            "Contiguous seed mode is not supported for this sensitivity setting.".to_string()
        })?;
        vec![contiguous_seed.to_string()]
    } else if input.shape_mask.is_empty() {
        get_shape_codes(sens)
            .iter()
            .map(|x| x.to_string())
            .collect()
    } else {
        input.shape_mask.clone()
    };
    let shapes = ShapeConfig::from_codes(&shape_codes, input.shapes, &reduction)?;
    if (input.lin_stage1_target || input.lin_stage1_query || input.lin_stage1_combo)
        && shapes[0].weight < 10
    {
        return Err("Linearization is only supported for seed shapes of weight >= 10.".to_string());
    }

    let seed_cut = if input.seed_cut == 0.0 {
        traits.seed_cut
    } else {
        input.seed_cut
    };
    let seed_complexity_cut = seed_cut * std::f64::consts::LN_2 * shapes[0].weight as f64;
    let mut soft_masking_bits = soft_masking_algo(
        &traits,
        &input.motif_masking,
        input.swipe_all,
        input.freq_masking,
    )?
    .bits();
    if !input.soft_masking.is_empty() {
        soft_masking_bits |= match input.soft_masking.to_lowercase().as_str() {
            "0" | "none" => MaskingAlgo::None.bits(),
            "1" | "tantan" => MaskingAlgo::Tantan.bits(),
            "seg" => MaskingAlgo::Seg.bits(),
            _ => {
                return Err(format!(
                    "Invalid value for string field: {}. Permitted values: none, tantan, seg",
                    input.soft_masking
                ))
            }
        };
    }

    let extension_mode = if input.ext.is_empty() {
        if input.global_ranking_targets
            || input.swipe_all
            || ((input.lin_stage1_query || input.lin_stage1_target || input.lin_stage1_combo)
                && !input.linclust_banded_ext)
        {
            ExtensionMode::Full
        } else {
            default_ext_mode(sens)
        }
    } else {
        let mode = ExtensionMode::from_string(&input.ext)?;
        if mode != ExtensionMode::Full {
            if input.global_ranking_targets {
                return Err("Global ranking only supports full matrix extension.".to_string());
            }
            if input.swipe_all {
                return Err("--swipe only supports full matrix extension.".to_string());
            }
        }
        mode
    };
    if extension_mode == ExtensionMode::Full {
        if input.frame_shift > 0 {
            return Err("Frameshift alignment does not support full matrix extension.".to_string());
        }
        if input.anchored_swipe {
            return Err("Anchored swipe does not support full matrix extension.".to_string());
        }
    }
    if !input.aln_out.is_empty() && input.parallel_tmpdir.is_empty() {
        return Err(
            "Alignment output to file is not supported without --parallel-tmpdir.".to_string(),
        );
    }

    Ok(SetupSearchResult {
        sensitivity: sens,
        freq_sd,
        hamming_filter_id,
        ungapped_evalue,
        ungapped_evalue_short,
        gapped_filter_evalue,
        query_bins,
        minimizer_window,
        sketch_size,
        reduction,
        shapes,
        gapped_filter_diag_score: score_matrix.rawscore_int(input.gapped_filter_diag_bit_score),
        seed_complexity_cut,
        soft_masking_bits,
        extension_mode,
    })
}

/// Matches C++ `Search::seedp_bits`.
pub fn seedp_bits(
    shape_weight: i32,
    threads: i32,
    index_chunks: i32,
    reduction: &Reduction,
) -> i32 {
    let seed_space = power(reduction.size() as i64, shape_weight as i64).saturating_sub(1) as u64;
    let seed_offset_bits = std::mem::size_of::<SeedOffset>() as i32 * 8;
    let by_seed_space = bit_length(seed_space) - seed_offset_bits;
    let by_threads = bit_length((threads as i64 * 4 * index_chunks as i64 - 1) as u64);
    by_seed_space.max(by_threads).max(10)
}

/// Matches C++ `Search::use_single_indexed`.
pub fn use_single_indexed(
    coverage: f64,
    query_letters: usize,
    ref_letters: usize,
    sensitivity: Sensitivity,
) -> bool {
    if coverage >= SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE {
        return false;
    }
    if sensitivity >= Sensitivity::Sensitive {
        query_letters < 300_000 && query_letters * 20_000 < ref_letters
    } else {
        query_letters < 3_000_000 && query_letters * 2_000 < ref_letters
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sensitivity_traits() {
        let t = get_traits(Sensitivity::Default);
        assert_eq!(t.index_chunks, 4);
        assert_eq!(t.contiguous_seed, Some("111111"));
        assert_eq!(get_traits(Sensitivity::Faster).sketch_size, 21);
        assert_eq!(get_traits(Sensitivity::Shapes6x10).seed_cut, 0.9);
    }

    #[test]
    fn test_shape_codes() {
        let codes = get_shape_codes(Sensitivity::Default);
        assert_eq!(codes.len(), 2);
        assert_eq!(codes[0], "111101110111");

        let codes = get_shape_codes(Sensitivity::UltraSensitive);
        assert_eq!(codes.len(), 64);

        let codes = get_shape_codes(Sensitivity::Linclust20);
        assert_eq!(codes.len(), 30);
        assert_eq!(codes[0], "111111111111");

        let codes = get_shape_codes(Sensitivity::Linclust40);
        assert_eq!(codes.len(), 15);

        let codes = get_shape_codes(Sensitivity::Shapes6x10);
        assert_eq!(codes.len(), 6);

        let codes = get_shape_codes(Sensitivity::Shapes30x10);
        assert_eq!(codes.len(), 30);
    }

    #[test]
    fn test_fast_shape() {
        let codes = get_shape_codes(Sensitivity::Fast);
        assert_eq!(codes.len(), 1);
        assert_eq!(codes[0], "1101110101101111");
    }

    #[test]
    fn test_seedp_bits() {
        let reduction = Reduction::default_reduction();
        assert_eq!(seedp_bits(16, 1, 4, &reduction), 22);
        assert_eq!(seedp_bits(4, 1, 1, &reduction), 10);
        assert_eq!(seedp_bits(4, 64, 4, &reduction), 10);
        assert_eq!(seedp_bits(4, 1024, 4, &reduction), 14);
    }

    #[test]
    fn test_use_single_indexed() {
        assert!(!use_single_indexed(
            0.15,
            1,
            1_000_000_000,
            Sensitivity::Fast
        ));
        assert!(use_single_indexed(
            0.1,
            100_000,
            300_000_000,
            Sensitivity::Fast
        ));
        assert!(!use_single_indexed(
            0.1,
            3_000_000,
            10_000_000_000,
            Sensitivity::Fast
        ));
        assert!(use_single_indexed(
            0.1,
            10_000,
            300_000_000,
            Sensitivity::Sensitive
        ));
        assert!(!use_single_indexed(
            0.1,
            300_000,
            10_000_000_000,
            Sensitivity::Sensitive
        ));
    }

    #[test]
    fn test_hamming_id_cutoff() {
        assert_eq!(hamming_id_cutoff(49.9), 0);
        assert_eq!(hamming_id_cutoff(50.0), 20);
        assert_eq!(hamming_id_cutoff(89.9), 20);
        assert_eq!(hamming_id_cutoff(90.0), 30);
    }

    #[test]
    fn test_iterated_sens() {
        assert_eq!(iterated_sens(Sensitivity::Faster), &[]);
        assert_eq!(
            iterated_sens(Sensitivity::Default),
            &[
                Round::new(Sensitivity::Fast, true),
                Round::new(Sensitivity::Linclust40, true)
            ]
        );
        assert_eq!(
            iterated_sens(Sensitivity::VerySensitive),
            &[
                Round::new(Sensitivity::Fast, true),
                Round::new(Sensitivity::Linclust20, true),
                Round::new(Sensitivity::Default, false),
                Round::new(Sensitivity::MoreSensitive, false)
            ]
        );
    }

    #[test]
    fn test_soft_masking_algo() {
        let motif_traits = get_traits(Sensitivity::Fast);
        assert_eq!(
            soft_masking_algo(&motif_traits, "", false, false).unwrap(),
            MaskingAlgo::Motif
        );
        assert_eq!(
            soft_masking_algo(&motif_traits, "", false, true).unwrap(),
            MaskingAlgo::None
        );
        assert_eq!(
            soft_masking_algo(&motif_traits, "0", false, false).unwrap(),
            MaskingAlgo::None
        );
        assert_eq!(
            soft_masking_algo(&motif_traits, "1", false, false).unwrap(),
            MaskingAlgo::Motif
        );
        assert!(soft_masking_algo(&motif_traits, "1", true, false).is_err());
        assert!(soft_masking_algo(&motif_traits, "bad", false, false).is_err());
    }

    #[test]
    fn test_extension_mode_from_string_and_defaults() {
        assert_eq!(
            ExtensionMode::from_string("banded-fast").unwrap(),
            ExtensionMode::BandedFast
        );
        assert_eq!(
            ExtensionMode::from_string("full").unwrap(),
            ExtensionMode::Full
        );
        assert!(ExtensionMode::from_string("bad").is_err());
        assert_eq!(
            default_ext_mode(Sensitivity::Sensitive),
            ExtensionMode::BandedFast
        );
        assert_eq!(
            default_ext_mode(Sensitivity::MoreSensitive),
            ExtensionMode::BandedSlow
        );
    }

    #[test]
    fn test_setup_search_defaults_and_shape_setup() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let input = SetupSearchInput::default();
        let result = setup_search(Sensitivity::Default, &input, &score_matrix).unwrap();
        assert_eq!(result.sensitivity, Sensitivity::Default);
        assert_eq!(result.freq_sd, 50.0);
        assert_eq!(result.hamming_filter_id, 11);
        assert_eq!(result.ungapped_evalue, 10000.0);
        assert_eq!(result.query_bins, 16);
        assert_eq!(result.extension_mode, ExtensionMode::BandedFast);
        assert_eq!(result.shapes.count(), 2);
        assert_eq!(format!("{}", result.shapes), "111101110111,111011010010111");
        assert_eq!(
            result.seed_complexity_cut,
            0.8 * std::f64::consts::LN_2 * 10.0
        );
        assert_ne!(result.soft_masking_bits & MaskingAlgo::Motif.bits(), 0);
    }

    #[test]
    fn test_setup_search_overrides_and_full_extension_guards() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let mut input = SetupSearchInput {
            freq_sd: Some(7.0),
            min_identities: Some(3),
            approx_min_id: 90.0,
            query_bins: Some(2),
            shape_mask: vec!["1111111111".to_string()],
            lin_stage1_query: true,
            soft_masking: "tantan".to_string(),
            ext: "full".to_string(),
            ..SetupSearchInput::default()
        };
        let result = setup_search(Sensitivity::Fast, &input, &score_matrix).unwrap();
        assert_eq!(result.freq_sd, 7.0);
        assert_eq!(result.hamming_filter_id, 30);
        assert_eq!(result.query_bins, 2);
        assert_eq!(result.shapes.count(), 1);
        assert_ne!(result.soft_masking_bits & MaskingAlgo::Tantan.bits(), 0);
        assert_eq!(result.extension_mode, ExtensionMode::Full);

        input.frame_shift = 1;
        assert_eq!(
            setup_search(Sensitivity::Fast, &input, &score_matrix)
                .err()
                .unwrap(),
            "Frameshift alignment does not support full matrix extension."
        );
    }

    #[test]
    fn test_setup_search_contiguous_and_error_paths() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let input = SetupSearchInput {
            contiguous_seed_mode: true,
            ..SetupSearchInput::default()
        };
        let result = setup_search(Sensitivity::Default, &input, &score_matrix).unwrap();
        assert_eq!(format!("{}", result.shapes), "111111");
        assert_eq!(
            format!("{}", result.reduction),
            "[RK][QE][D][N][C][G][H][F][Y][IV][LM][W][P][S][T][A]"
        );
        assert_eq!(
            setup_search(Sensitivity::Fast, &input, &score_matrix)
                .err()
                .unwrap(),
            "Contiguous seed mode is not supported for this sensitivity setting."
        );

        let input = SetupSearchInput {
            shape_mask: vec!["111".to_string()],
            lin_stage1_target: true,
            ..SetupSearchInput::default()
        };
        assert_eq!(
            setup_search(Sensitivity::Fast, &input, &score_matrix)
                .err()
                .unwrap(),
            "Linearization is only supported for seed shapes of weight >= 10."
        );
    }

    #[test]
    fn test_setup_search_extension_and_output_errors() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let input = SetupSearchInput {
            ext: "banded-fast".to_string(),
            global_ranking_targets: true,
            ..SetupSearchInput::default()
        };
        assert_eq!(
            setup_search(Sensitivity::Fast, &input, &score_matrix)
                .err()
                .unwrap(),
            "Global ranking only supports full matrix extension."
        );

        let input = SetupSearchInput {
            ext: "banded-fast".to_string(),
            swipe_all: true,
            ..SetupSearchInput::default()
        };
        assert_eq!(
            setup_search(Sensitivity::Fast, &input, &score_matrix)
                .err()
                .unwrap(),
            "--swipe only supports full matrix extension."
        );

        let input = SetupSearchInput {
            aln_out: "x.daa".to_string(),
            ..SetupSearchInput::default()
        };
        assert_eq!(
            setup_search(Sensitivity::Fast, &input, &score_matrix)
                .err()
                .unwrap(),
            "Alignment output to file is not supported without --parallel-tmpdir."
        );
    }
}
