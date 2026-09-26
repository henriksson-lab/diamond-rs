//! Translation boundary for `diamond/src/search/setup.cpp`.
//!
//! Mutable C++ globals (configuration, score matrix, reduction, and shape
//! table) are explicit inputs and outputs here. The implementation is shared
//! with the older [`super::sensitivity`] compatibility module.

use crate::basic::reduction::Reduction;
use crate::config::Sensitivity;
use crate::masking::MaskingAlgo;
use crate::stats::score_matrix::ScoreMatrix;

pub use super::sensitivity::{
    ExtensionMode, Round, SensitivityTraits, SetupSearchInput, SetupSearchResult,
    DEFAULT_CONTIGUOUS_REDUCTION, SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE,
};

pub fn sensitivity_traits(sensitivity: Sensitivity) -> SensitivityTraits {
    super::sensitivity::get_traits(sensitivity)
}

pub fn shape_codes(sensitivity: Sensitivity) -> &'static [&'static str] {
    super::sensitivity::get_shape_codes(sensitivity)
}

pub fn iterated_sens(sensitivity: Sensitivity) -> &'static [Round] {
    super::sensitivity::iterated_sens(sensitivity)
}

pub fn approx_id_to_hamming_id() -> &'static [(f64, u32)] {
    super::sensitivity::approx_id_to_hamming_id()
}

pub fn hamming_id_cutoff(approx_id: f64) -> u32 {
    super::sensitivity::hamming_id_cutoff(approx_id)
}

pub fn seedp_bits(
    shape_weight: i32,
    threads: i32,
    index_chunks: i32,
    reduction: &Reduction,
) -> i32 {
    super::sensitivity::seedp_bits(shape_weight, threads, index_chunks, reduction)
}

pub fn use_single_indexed(
    coverage: f64,
    query_letters: usize,
    reference_letters: usize,
    sensitivity: Sensitivity,
) -> bool {
    super::sensitivity::use_single_indexed(coverage, query_letters, reference_letters, sensitivity)
}

pub fn soft_masking_algo(
    traits: &SensitivityTraits,
    motif_masking: &str,
    swipe_all: bool,
    freq_masking: bool,
) -> Result<MaskingAlgo, String> {
    super::sensitivity::soft_masking_algo(traits, motif_masking, swipe_all, freq_masking)
}

pub fn setup_search(
    sensitivity: Sensitivity,
    input: &SetupSearchInput,
    score_matrix: &ScoreMatrix,
) -> Result<SetupSearchResult, String> {
    super::sensitivity::setup_search(sensitivity, input, score_matrix)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn score_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn explicit_override_and_sentinel_paths_match_set_option() {
        let matrix = score_matrix();
        let explicit = SetupSearchInput {
            min_identities: Some(3),
            approx_min_id: 90.0,
            ..SetupSearchInput::default()
        };
        assert_eq!(
            setup_search(Sensitivity::Fast, &explicit, &matrix)
                .unwrap()
                .hamming_filter_id,
            3
        );

        let sentinels = SetupSearchInput {
            freq_sd: Some(0.0),
            min_identities: Some(0),
            ungapped_evalue: Some(-1.0),
            ungapped_evalue_short: Some(-1.0),
            gapped_filter_evalue: Some(-1.0),
            minimizer_window: Some(0),
            sketch_size: Some(0),
            ..SetupSearchInput::default()
        };
        let result = setup_search(Sensitivity::Faster, &sentinels, &matrix).unwrap();
        assert_eq!(result.freq_sd, 50.0);
        assert_eq!(result.hamming_filter_id, 11);
        assert_eq!(result.ungapped_evalue, 0.0);
        assert_eq!(result.ungapped_evalue_short, 0.0);
        assert_eq!(result.gapped_filter_evalue, 0.0);
        assert_eq!(result.minimizer_window, 0);
        assert_eq!(result.sketch_size, 21);
    }

    #[test]
    fn setup_builds_cutoff_tables_and_restores_murphy_reduction() {
        let matrix = score_matrix();
        let input = SetupSearchInput {
            contiguous_seed_mode: true,
            ungapped_evalue: Some(10_000.0),
            ungapped_evalue_short: Some(30_000.0),
            ..SetupSearchInput::default()
        };
        let result = setup_search(Sensitivity::Default, &input, &matrix).unwrap();
        assert_eq!(format!("{}", result.shapes), "111111");
        assert_eq!(result.reduction.size(), 10);
        assert_eq!(
            result.cutoff_table.get(128),
            matrix.rawscore_int(matrix.bitscore_norm(10_000.0, 128))
        );
        assert_eq!(
            result.cutoff_table_short.get(128),
            matrix.rawscore_int(matrix.bitscore_norm(30_000.0, 128))
        );
    }

    #[test]
    fn soft_masking_string_parser_is_exact_and_case_sensitive() {
        let matrix = score_matrix();
        for invalid in ["1", "SEG", "TANTAN"] {
            let input = SetupSearchInput {
                soft_masking: invalid.to_string(),
                ..SetupSearchInput::default()
            };
            assert_eq!(
                setup_search(Sensitivity::Fast, &input, &matrix)
                    .err()
                    .unwrap(),
                format!(
                    "Invalid value for string field: {invalid}. Permitted values: 0, none, seg, tantan"
                )
            );
        }
    }

    #[test]
    #[should_panic(expected = "no iterated-search schedule for shapes-6x10")]
    fn shapes_6x10_has_no_upstream_iterated_schedule() {
        let _ = iterated_sens(Sensitivity::Shapes6x10);
    }
}
