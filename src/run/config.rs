//! Translation of `diamond/src/run/config.{h,cpp}`.
//!
//! C++ obtains most constructor inputs from mutable process-wide globals. This
//! module makes that dependency explicit through [`ConfigInput`] and models
//! `init_output(max_target_seqs)` as a callback.

use std::cmp::Ordering;

use crate::config::Sensitivity;
use crate::data::flags::SeedEncoding;
use crate::masking::{MaskingAlgo, MaskingMode};
use crate::search::sensitivity::iterated_sens;

/// C++ `Search::Round`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Round {
    pub sensitivity: Sensitivity,
    pub linearize: bool,
}

impl Round {
    /// C++ `Round::Round(Sensitivity, bool)`.
    pub const fn new(sensitivity: Sensitivity, linearize: bool) -> Self {
        Self {
            sensitivity,
            linearize,
        }
    }

    /// C++ `Round::operator<`.
    pub fn less_than(self, other: Self) -> bool {
        (self.linearize && !other.linearize)
            || (self.linearize == other.linearize && self.sensitivity < other.sensitivity)
    }

    /// C++ `Round::operator==`.
    pub fn equals(self, other: Self) -> bool {
        self.sensitivity == other.sensitivity && self.linearize == other.linearize
    }

    /// C++ `Round::operator!=`.
    pub fn not_equal(self, other: Self) -> bool {
        self.sensitivity != other.sensitivity || self.linearize != other.linearize
    }
}

impl Ord for Round {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.less_than(*other) {
            Ordering::Less
        } else if other.less_than(*self) {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

impl PartialOrd for Round {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// C++ `::Config::Algo` values used by this constructor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algo {
    Auto,
    DoubleIndexed,
    QueryIndexed,
    CtgSeed,
}

/// Explicit replacement for the globals read by C++ `Search::Config::Config`.
#[derive(Debug, Clone, PartialEq)]
pub struct ConfigInput {
    pub self_search: bool,
    pub target_indexed: bool,
    /// `None` means `--iterate` was absent; `Some([])` means it was supplied
    /// without explicit sensitivity levels.
    pub iterate: Option<Vec<String>>,
    pub multiprocessing: bool,
    pub lin_stage1_query: bool,
    pub lin_stage1_target: bool,
    pub sensitivity: Sensitivity,
    pub unaligned: String,
    pub aligned_file: String,
    pub taxonlist: String,
    pub taxon_exclude: String,
    pub global_ranking_targets: i64,
    pub frame_shift: i32,
    pub algo: Algo,
    pub blastn: bool,
    pub masking: MaskingMode,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub freq_masking: bool,
    pub seed_cut: f64,
    pub freq_sd: f64,
    pub minimizer_window: i32,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub min_length_ratio: f64,
    pub query_translated: bool,
}

impl Default for ConfigInput {
    fn default() -> Self {
        Self {
            self_search: false,
            target_indexed: false,
            iterate: None,
            multiprocessing: false,
            lin_stage1_query: false,
            lin_stage1_target: false,
            sensitivity: Sensitivity::Default,
            unaligned: String::new(),
            aligned_file: String::new(),
            taxonlist: String::new(),
            taxon_exclude: String::new(),
            global_ranking_targets: 0,
            frame_shift: 0,
            algo: Algo::Auto,
            blastn: false,
            masking: MaskingMode::Tantan,
            gap_open: -1,
            gap_extend: -1,
            freq_masking: false,
            seed_cut: 0.0,
            freq_sd: 0.0,
            minimizer_window: 0,
            query_cover: 0.0,
            subject_cover: 0.0,
            min_length_ratio: 0.0,
            query_translated: false,
        }
    }
}

/// Constructor-owned surface of C++ `Search::Config`.
///
/// Fields populated later by `setup_search` or workflow setup are deliberately
/// not fabricated here. `O` is the owned output-format object returned by the
/// explicit output initializer.
#[derive(Debug)]
pub struct Config<O = ()> {
    pub self_search: bool,
    pub sensitivity: Vec<Round>,
    pub seed_encoding: SeedEncoding,
    pub query_masking: MaskingAlgo,
    pub target_masking: MaskingAlgo,
    pub soft_masking: MaskingAlgo,
    pub lazy_masking: bool,
    pub track_aligned_queries: bool,
    pub lin_stage1_target: bool,
    pub min_length_ratio: f64,
    pub max_target_seqs: i64,
    pub output_format: O,
    pub iteration_query_aligned: u32,
    /// Effective values after the BLASTN `-1` defaults are applied.
    pub gap_open: i32,
    pub gap_extend: i32,
    /// Text written to C++ `message_stream`, without the trailing newline.
    pub iteration_message: Option<String>,
}

impl Config<()> {
    /// Convenience construction when no concrete output object is required.
    pub fn from_input(input: &ConfigInput) -> Result<Self, String> {
        Self::new(input, |_| Ok(()))
    }
}

impl<O> Config<O> {
    /// C++ `Search::Config::Config()` with globals and output creation made
    /// explicit. The callback receives the zero-initialized
    /// `max_target_seqs` by mutable reference, matching `init_output`.
    pub fn new<F>(input: &ConfigInput, init_output: F) -> Result<Self, String>
    where
        F: FnOnce(&mut i64) -> Result<O, String>,
    {
        let mut sensitivity = build_rounds(input)?;
        sensitivity.sort();
        if sensitivity.windows(2).any(|rounds| rounds[0] == rounds[1]) {
            return Err(
                "The same sensitivity level was specified multiple times for --iterate."
                    .to_string(),
            );
        }

        let iteration_message = (sensitivity.len() > 1).then(|| {
            let steps = sensitivity
                .iter()
                .map(|round| {
                    let name = sensitivity_name(round.sensitivity);
                    if round.linearize {
                        format!("{name} (linear)")
                    } else {
                        name.to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join(", ");
            format!("Running iterated search mode with sensitivity steps: {steps}")
        });
        let track_aligned_queries = iteration_message.is_some()
            || !input.unaligned.is_empty()
            || !input.aligned_file.is_empty();

        validate_modes(input)?;

        let (query_masking, target_masking) = if input.blastn {
            (MaskingAlgo::None, MaskingAlgo::None)
        } else {
            match input.masking {
                MaskingMode::BlastSeg => (MaskingAlgo::None, MaskingAlgo::Seg),
                MaskingMode::Tantan => (MaskingAlgo::Tantan, MaskingAlgo::Tantan),
                MaskingMode::None => (MaskingAlgo::None, MaskingAlgo::None),
            }
        };
        let gap_open = if input.blastn && input.gap_open == -1 {
            5
        } else {
            input.gap_open
        };
        let gap_extend = if input.blastn && input.gap_extend == -1 {
            2
        } else {
            input.gap_extend
        };

        let min_length_ratio = derive_min_length_ratio(input, &sensitivity)?;
        let mut max_target_seqs = 0;
        let output_format = init_output(&mut max_target_seqs)?;

        Ok(Self {
            self_search: input.self_search,
            sensitivity,
            seed_encoding: if input.target_indexed {
                SeedEncoding::Hashed
            } else {
                SeedEncoding::SpacedFactor
            },
            query_masking,
            target_masking,
            soft_masking: MaskingAlgo::None,
            lazy_masking: false,
            track_aligned_queries,
            // C++ initializes this runtime field to false; setup_search fills
            // it later from the selected sensitivity traits.
            lin_stage1_target: false,
            min_length_ratio,
            max_target_seqs,
            output_format,
            iteration_query_aligned: 0,
            gap_open,
            gap_extend,
            iteration_message,
        })
    }

    /// C++ inline `Search::Config::iterated()`.
    pub fn iterated(&self) -> bool {
        self.sensitivity.len() > 1
    }

    /// C++ `Search::Config::free()` is intentionally empty.
    pub fn free(&mut self) {}
}

fn build_rounds(input: &ConfigInput) -> Result<Vec<Round>, String> {
    let mut rounds = Vec::new();
    if let Some(iterate) = &input.iterate {
        if input.multiprocessing {
            return Err("Iterated search is not compatible with --multiprocessing.".to_string());
        }
        if input.target_indexed {
            return Err("Iterated search is not compatible with --target-indexed.".to_string());
        }
        if input.self_search {
            return Err("Iterated search is not compatible with --self.".to_string());
        }
        if input.lin_stage1_query {
            return Err("Iterated search is not compatible with --lin-stage1.".to_string());
        }

        if iterate.is_empty() {
            rounds.push(Round::new(Sensitivity::Faster, true));
            rounds.extend(
                iterated_sens(input.sensitivity)
                    .iter()
                    .map(|round| Round::new(round.sensitivity, round.query_indexed)),
            );
        } else {
            let target = Round::new(
                input.sensitivity,
                input.lin_stage1_query || input.lin_stage1_target,
            );
            for value in iterate {
                let (name, linearize) = match value.strip_suffix("_lin") {
                    Some(name) => (name, true),
                    None => (value.as_str(), false),
                };
                let round = Round::new(parse_sensitivity(name)?, linearize);
                if round >= target {
                    return Err(
                        "Sensitivity levels set for --iterate must be below target sensitivity."
                            .to_string(),
                    );
                }
                rounds.push(round);
            }
        }
    }

    let target = Round::new(input.sensitivity, input.lin_stage1_target);
    if rounds.last().copied() != Some(target) {
        rounds.push(target);
    }
    Ok(rounds)
}

fn validate_modes(input: &ConfigInput) -> Result<(), String> {
    if input.multiprocessing && (!input.taxonlist.is_empty() || !input.taxon_exclude.is_empty()) {
        return Err("Multiprocessing mode is not compatible with database filtering.".to_string());
    }
    if input.global_ranking_targets != 0 {
        if input.frame_shift != 0 {
            return Err(
                "Global ranking mode is not compatible with frameshift alignments.".to_string(),
            );
        }
        if input.multiprocessing {
            return Err(
                "Global ranking mode is not compatible with --multiprocessing.".to_string(),
            );
        }
    }
    if input.target_indexed && !matches!(input.algo, Algo::Auto | Algo::DoubleIndexed) {
        return Err("--target-indexed requires --algo 0".to_string());
    }
    if input.freq_masking && input.seed_cut != 0.0 {
        return Err("Incompatible options: --freq-masking, --seed-cut.".to_string());
    }
    if input.freq_sd != 0.0 && !input.freq_masking {
        return Err("--freq-sd requires --freq-masking.".to_string());
    }
    if input.minimizer_window != 0 && input.algo == Algo::CtgSeed {
        return Err("Minimizer setting is not compatible with contiguous seed mode.".to_string());
    }
    Ok(())
}

fn derive_min_length_ratio(input: &ConfigInput, rounds: &[Round]) -> Result<f64, String> {
    if input.query_cover >= 50.0
        && input.query_cover == input.subject_cover
        && input.min_length_ratio == 0.0
        && !input.query_translated
    {
        let coverage = input.query_cover / 100.0;
        Ok(
            if input.lin_stage1_query
                && rounds
                    .last()
                    .is_some_and(|round| round.sensitivity < Sensitivity::Linclust40)
            {
                (coverage + 0.05).min(0.92)
            } else {
                (coverage - 0.05).max(0.0)
            },
        )
    } else if input.query_translated && input.min_length_ratio != 0.0 {
        Err("--min-len-ratio is not supported for translated searches".to_string())
    } else {
        Ok(input.min_length_ratio)
    }
}

fn parse_sensitivity(value: &str) -> Result<Sensitivity, String> {
    match value {
        "faster" => Ok(Sensitivity::Faster),
        "fast" => Ok(Sensitivity::Fast),
        "default" => Ok(Sensitivity::Default),
        "linclust-40" => Ok(Sensitivity::Linclust40),
        "linclust-20" => Ok(Sensitivity::Linclust20),
        "shapes-6x10" => Ok(Sensitivity::Shapes6x10),
        "shapes-30x10" => Ok(Sensitivity::Shapes30x10),
        "mid-sensitive" => Ok(Sensitivity::MidSensitive),
        "sensitive" => Ok(Sensitivity::Sensitive),
        "more-sensitive" => Ok(Sensitivity::MoreSensitive),
        "very-sensitive" => Ok(Sensitivity::VerySensitive),
        "ultra-sensitive" => Ok(Sensitivity::UltraSensitive),
        _ => Err(format!("Invalid sensitivity level: {value}")),
    }
}

fn sensitivity_name(value: Sensitivity) -> &'static str {
    match value {
        Sensitivity::Faster => "faster",
        Sensitivity::Fast => "fast",
        Sensitivity::Default => "default",
        Sensitivity::Linclust40 => "linclust-40",
        Sensitivity::Linclust20 => "linclust-20",
        Sensitivity::Shapes6x10 => "shapes-6x10",
        Sensitivity::Shapes30x10 => "shapes-30x10",
        Sensitivity::MidSensitive => "mid-sensitive",
        Sensitivity::Sensitive => "sensitive",
        Sensitivity::MoreSensitive => "more-sensitive",
        Sensitivity::VerySensitive => "very-sensitive",
        Sensitivity::UltraSensitive => "ultra-sensitive",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_inline_operators_match_cpp_ordering() {
        let linear_default = Round::new(Sensitivity::Default, true);
        let normal_fast = Round::new(Sensitivity::Fast, false);
        let normal_default = Round::new(Sensitivity::Default, false);
        assert!(linear_default.less_than(normal_fast));
        assert!(normal_fast.less_than(normal_default));
        assert!(normal_default.equals(normal_default));
        assert!(normal_default.not_equal(normal_fast));
        assert!(!normal_default.not_equal(normal_default));
    }

    #[test]
    fn default_constructor_maps_defaults_and_output_reference() {
        let config = Config::new(&ConfigInput::default(), |max_targets| {
            assert_eq!(*max_targets, 0);
            *max_targets = 25;
            Ok("format")
        })
        .unwrap();
        assert_eq!(
            config.sensitivity,
            [Round::new(Sensitivity::Default, false)]
        );
        assert_eq!(config.seed_encoding, SeedEncoding::SpacedFactor);
        assert_eq!(config.query_masking, MaskingAlgo::Tantan);
        assert_eq!(config.target_masking, MaskingAlgo::Tantan);
        assert_eq!(config.max_target_seqs, 25);
        assert_eq!(config.output_format, "format");
        assert!(!config.iterated());
    }

    #[test]
    fn default_iteration_builds_sorts_and_formats_rounds() {
        let input = ConfigInput {
            iterate: Some(Vec::new()),
            sensitivity: Sensitivity::VerySensitive,
            ..ConfigInput::default()
        };
        let config = Config::from_input(&input).unwrap();
        assert_eq!(
            config.sensitivity,
            [
                Round::new(Sensitivity::Faster, true),
                Round::new(Sensitivity::Fast, true),
                Round::new(Sensitivity::Linclust20, true),
                Round::new(Sensitivity::Default, false),
                Round::new(Sensitivity::MoreSensitive, false),
                Round::new(Sensitivity::VerySensitive, false),
            ]
        );
        assert!(config.track_aligned_queries);
        assert_eq!(
            config.iteration_message.as_deref(),
            Some("Running iterated search mode with sensitivity steps: faster (linear), fast (linear), linclust-20 (linear), default, more-sensitive, very-sensitive")
        );
    }

    #[test]
    fn explicit_iteration_suffix_and_duplicate_rules_match_cpp() {
        let input = ConfigInput {
            iterate: Some(vec!["fast_lin".into(), "default".into()]),
            sensitivity: Sensitivity::Sensitive,
            ..ConfigInput::default()
        };
        let config = Config::from_input(&input).unwrap();
        assert_eq!(config.sensitivity[0], Round::new(Sensitivity::Fast, true));
        assert_eq!(
            config.sensitivity[1],
            Round::new(Sensitivity::Default, false)
        );

        let duplicate = ConfigInput {
            iterate: Some(vec!["fast".into(), "fast".into()]),
            sensitivity: Sensitivity::Sensitive,
            ..ConfigInput::default()
        };
        assert_eq!(
            Config::from_input(&duplicate).unwrap_err(),
            "The same sensitivity level was specified multiple times for --iterate."
        );
    }

    #[test]
    fn iteration_and_mode_incompatibilities_are_reported() {
        let iterated_mp = ConfigInput {
            iterate: Some(Vec::new()),
            multiprocessing: true,
            ..ConfigInput::default()
        };
        assert_eq!(
            Config::from_input(&iterated_mp).unwrap_err(),
            "Iterated search is not compatible with --multiprocessing."
        );

        let filtered_mp = ConfigInput {
            multiprocessing: true,
            taxonlist: "2".into(),
            ..ConfigInput::default()
        };
        assert_eq!(
            Config::from_input(&filtered_mp).unwrap_err(),
            "Multiprocessing mode is not compatible with database filtering."
        );

        let bad_index = ConfigInput {
            target_indexed: true,
            algo: Algo::QueryIndexed,
            ..ConfigInput::default()
        };
        assert_eq!(
            Config::from_input(&bad_index).unwrap_err(),
            "--target-indexed requires --algo 0"
        );
    }

    #[test]
    fn masking_blastn_gaps_and_seed_encoding_match_cpp() {
        let protein = Config::from_input(&ConfigInput {
            masking: MaskingMode::BlastSeg,
            target_indexed: true,
            algo: Algo::DoubleIndexed,
            ..ConfigInput::default()
        })
        .unwrap();
        assert_eq!(protein.query_masking, MaskingAlgo::None);
        assert_eq!(protein.target_masking, MaskingAlgo::Seg);
        assert_eq!(protein.seed_encoding, SeedEncoding::Hashed);

        let blastn = Config::from_input(&ConfigInput {
            blastn: true,
            masking: MaskingMode::Tantan,
            ..ConfigInput::default()
        })
        .unwrap();
        assert_eq!(blastn.query_masking, MaskingAlgo::None);
        assert_eq!(blastn.target_masking, MaskingAlgo::None);
        assert_eq!((blastn.gap_open, blastn.gap_extend), (5, 2));
    }

    #[test]
    fn min_length_ratio_derivation_and_translated_rejection_match_cpp() {
        let symmetric = Config::from_input(&ConfigInput {
            query_cover: 80.0,
            subject_cover: 80.0,
            ..ConfigInput::default()
        })
        .unwrap();
        assert!((symmetric.min_length_ratio - 0.75).abs() < f64::EPSILON);

        let linear = Config::from_input(&ConfigInput {
            sensitivity: Sensitivity::Fast,
            lin_stage1_query: true,
            query_cover: 90.0,
            subject_cover: 90.0,
            ..ConfigInput::default()
        })
        .unwrap();
        assert!((linear.min_length_ratio - 0.92).abs() < f64::EPSILON);

        let translated = ConfigInput {
            query_translated: true,
            min_length_ratio: 0.5,
            ..ConfigInput::default()
        };
        assert_eq!(
            Config::from_input(&translated).unwrap_err(),
            "--min-len-ratio is not supported for translated searches"
        );
    }

    #[test]
    fn frequency_minimizer_and_global_ranking_checks_match_cpp() {
        let cases = [
            (
                ConfigInput {
                    freq_masking: true,
                    seed_cut: 1.0,
                    ..ConfigInput::default()
                },
                "Incompatible options: --freq-masking, --seed-cut.",
            ),
            (
                ConfigInput {
                    freq_sd: 1.0,
                    ..ConfigInput::default()
                },
                "--freq-sd requires --freq-masking.",
            ),
            (
                ConfigInput {
                    minimizer_window: 10,
                    algo: Algo::CtgSeed,
                    ..ConfigInput::default()
                },
                "Minimizer setting is not compatible with contiguous seed mode.",
            ),
            (
                ConfigInput {
                    global_ranking_targets: 10,
                    frame_shift: 15,
                    ..ConfigInput::default()
                },
                "Global ranking mode is not compatible with frameshift alignments.",
            ),
        ];
        for (input, message) in cases {
            assert_eq!(Config::from_input(&input).unwrap_err(), message);
        }
    }

    #[test]
    fn output_initializer_is_not_called_after_validation_failure() {
        let input = ConfigInput {
            freq_sd: 1.0,
            ..ConfigInput::default()
        };
        let mut called = false;
        let error = Config::new(&input, |_| {
            called = true;
            Ok(())
        })
        .unwrap_err();
        assert_eq!(error, "--freq-sd requires --freq-masking.");
        assert!(!called);
    }
}
