use crate::config::Sensitivity;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SensitivityTraits {
    pub query_indexed: bool,
    pub seed_freq_masking: bool,
    pub freq_sd: f64,
    pub shapes: i32,
    pub index_mode: i32,
    pub freq_masking: i32,
    pub seed_cut: i32,
    pub block_size: i32,
    pub minimizer_window: i32,
    pub shape_mask: Option<&'static str>,
    pub ungapped_evalue: f64,
    pub gapped_filter_evalue: f64,
    pub reduction: &'static str,
    pub min_chain_score: i32,
    pub chain_fraction_align: f64,
    pub band_extension: i32,
    pub max_overlap_extension: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DnaSearchConfig {
    pub sensitivity: Sensitivity,
    pub chain_fraction_align: f64,
    pub min_chain_score: i32,
    pub max_overlap_extension: f64,
    pub chain_pen_gap: f64,
    pub chain_pen_skip: f64,
    pub chain_fraction_align_option: Option<f64>,
    pub min_chain_score_option: Option<i32>,
    pub max_overlap_extension_option: Option<f64>,
    pub chain_pen_gap_scale: f64,
    pub chain_pen_skip_scale: f64,
}

impl Default for DnaSearchConfig {
    fn default() -> Self {
        Self {
            sensitivity: Sensitivity::Default,
            chain_fraction_align: 0.0,
            min_chain_score: 0,
            max_overlap_extension: 0.0,
            chain_pen_gap: 0.0,
            chain_pen_skip: 0.0,
            chain_fraction_align_option: None,
            min_chain_score_option: None,
            max_overlap_extension_option: None,
            chain_pen_gap_scale: 1.0,
            chain_pen_skip_scale: 1.0,
        }
    }
}

pub fn sensitivity_traits(sens: Sensitivity) -> Option<SensitivityTraits> {
    match sens {
        Sensitivity::Faster => Some(SensitivityTraits {
            query_indexed: true,
            seed_freq_masking: false,
            freq_sd: 20.0,
            shapes: 9,
            index_mode: 0,
            freq_masking: 0,
            seed_cut: 0,
            block_size: 1,
            minimizer_window: 16,
            shape_mask: None,
            ungapped_evalue: 1.0,
            gapped_filter_evalue: 2.0,
            reduction: "dna",
            min_chain_score: 12,
            chain_fraction_align: 0.0,
            band_extension: 40,
            max_overlap_extension: 0.1,
        }),
        Sensitivity::Fast => Some(SensitivityTraits {
            query_indexed: true,
            seed_freq_masking: false,
            freq_sd: 20.0,
            shapes: 9,
            index_mode: 0,
            freq_masking: 0,
            seed_cut: 0,
            block_size: 1,
            minimizer_window: 16,
            shape_mask: None,
            ungapped_evalue: 1.0,
            gapped_filter_evalue: 2.0,
            reduction: "dna",
            min_chain_score: 10,
            chain_fraction_align: 0.0,
            band_extension: 40,
            max_overlap_extension: 0.5,
        }),
        Sensitivity::Default => Some(SensitivityTraits {
            query_indexed: true,
            seed_freq_masking: false,
            freq_sd: 20.0,
            shapes: 9,
            index_mode: 0,
            freq_masking: 0,
            seed_cut: 0,
            block_size: 1,
            minimizer_window: 16,
            shape_mask: None,
            ungapped_evalue: 1.0,
            gapped_filter_evalue: 2.0,
            reduction: "dna",
            min_chain_score: 10,
            chain_fraction_align: 0.0,
            band_extension: 20,
            max_overlap_extension: 0.5,
        }),
        Sensitivity::Sensitive => Some(SensitivityTraits {
            query_indexed: true,
            seed_freq_masking: false,
            freq_sd: 20.0,
            shapes: 9,
            index_mode: 0,
            freq_masking: 0,
            seed_cut: 0,
            block_size: 1,
            minimizer_window: 16,
            shape_mask: None,
            ungapped_evalue: 1.0,
            gapped_filter_evalue: 2.0,
            reduction: "dna",
            min_chain_score: 6,
            chain_fraction_align: 0.0,
            band_extension: 20,
            max_overlap_extension: 0.5,
        }),
        Sensitivity::VerySensitive => Some(SensitivityTraits {
            query_indexed: true,
            seed_freq_masking: false,
            freq_sd: 20.0,
            shapes: 9,
            index_mode: 0,
            freq_masking: 0,
            seed_cut: 0,
            block_size: 1,
            minimizer_window: 16,
            shape_mask: None,
            ungapped_evalue: 1.0,
            gapped_filter_evalue: 2.0,
            reduction: "dna",
            min_chain_score: 5,
            chain_fraction_align: 0.0,
            band_extension: 17,
            max_overlap_extension: 0.5,
        }),
        Sensitivity::UltraSensitive => Some(SensitivityTraits {
            query_indexed: true,
            seed_freq_masking: false,
            freq_sd: 20.0,
            shapes: 9,
            index_mode: 0,
            freq_masking: 0,
            seed_cut: 0,
            block_size: 1,
            minimizer_window: 16,
            shape_mask: None,
            ungapped_evalue: 1.0,
            gapped_filter_evalue: 2.0,
            reduction: "dna",
            min_chain_score: 4,
            chain_fraction_align: 0.0,
            band_extension: 15,
            max_overlap_extension: 0.5,
        }),
        _ => None,
    }
}

pub fn shape_codes(sens: Sensitivity) -> Option<&'static [&'static str]> {
    match sens {
        Sensitivity::UltraSensitive => Some(&["111111111111"]),
        Sensitivity::VerySensitive => Some(&["1111111111111"]),
        Sensitivity::Sensitive => Some(&["11111111111111"]),
        Sensitivity::Default => Some(&["111111111111111"]),
        Sensitivity::Fast => Some(&["111111111111111"]),
        Sensitivity::Faster => Some(&["111111111111111111"]),
        _ => None,
    }
}

pub fn setup_search(sens: Sensitivity, cfg: &mut DnaSearchConfig, shape_length: i32) -> bool {
    let Some(traits) = sensitivity_traits(sens) else {
        return false;
    };

    cfg.sensitivity = sens;
    cfg.chain_fraction_align = cfg
        .chain_fraction_align_option
        .unwrap_or(traits.chain_fraction_align);
    cfg.min_chain_score = cfg.min_chain_score_option.unwrap_or(traits.min_chain_score);
    cfg.max_overlap_extension = cfg
        .max_overlap_extension_option
        .unwrap_or(traits.max_overlap_extension);
    cfg.chain_pen_gap = cfg.chain_pen_gap_scale * 0.01 * f64::from(shape_length);
    cfg.chain_pen_skip = cfg.chain_pen_skip_scale * 0.01 * f64::from(shape_length);
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sensitivity_traits_match_dna_table() {
        let faster = sensitivity_traits(Sensitivity::Faster).unwrap();
        assert_eq!(faster.min_chain_score, 12);
        assert_eq!(faster.band_extension, 40);
        assert_eq!(faster.max_overlap_extension, 0.1);

        let ultra = sensitivity_traits(Sensitivity::UltraSensitive).unwrap();
        assert_eq!(ultra.min_chain_score, 4);
        assert_eq!(ultra.band_extension, 15);
        assert_eq!(ultra.max_overlap_extension, 0.5);

        assert!(sensitivity_traits(Sensitivity::MoreSensitive).is_none());
    }

    #[test]
    fn test_shape_codes_match_dna_table() {
        assert_eq!(
            shape_codes(Sensitivity::Faster).unwrap(),
            &["111111111111111111"]
        );
        assert_eq!(
            shape_codes(Sensitivity::Fast).unwrap(),
            &["111111111111111"]
        );
        assert_eq!(
            shape_codes(Sensitivity::Default).unwrap(),
            &["111111111111111"]
        );
        assert_eq!(
            shape_codes(Sensitivity::Sensitive).unwrap(),
            &["11111111111111"]
        );
        assert_eq!(
            shape_codes(Sensitivity::VerySensitive).unwrap(),
            &["1111111111111"]
        );
        assert_eq!(
            shape_codes(Sensitivity::UltraSensitive).unwrap(),
            &["111111111111"]
        );
        assert!(shape_codes(Sensitivity::MoreSensitive).is_none());
    }

    #[test]
    fn test_setup_search_applies_overrides_and_chain_penalties() {
        let mut cfg = DnaSearchConfig {
            chain_fraction_align_option: Some(0.25),
            min_chain_score_option: Some(99),
            max_overlap_extension_option: Some(0.75),
            chain_pen_gap_scale: 2.0,
            chain_pen_skip_scale: 3.0,
            ..DnaSearchConfig::default()
        };

        assert!(setup_search(Sensitivity::Sensitive, &mut cfg, 15));

        assert_eq!(cfg.sensitivity, Sensitivity::Sensitive);
        assert_eq!(cfg.chain_fraction_align, 0.25);
        assert_eq!(cfg.min_chain_score, 99);
        assert_eq!(cfg.max_overlap_extension, 0.75);
        assert_eq!(cfg.chain_pen_gap, 0.3);
        assert!((cfg.chain_pen_skip - 0.45).abs() < f64::EPSILON);
    }

    #[test]
    fn test_setup_search_uses_trait_defaults() {
        let mut cfg = DnaSearchConfig::default();

        assert!(setup_search(Sensitivity::VerySensitive, &mut cfg, 13));

        assert_eq!(cfg.min_chain_score, 5);
        assert_eq!(cfg.chain_fraction_align, 0.0);
        assert_eq!(cfg.max_overlap_extension, 0.5);
        assert_eq!(cfg.chain_pen_gap, 0.13);
        assert_eq!(cfg.chain_pen_skip, 0.13);
        assert!(!setup_search(Sensitivity::MidSensitive, &mut cfg, 13));
    }
}
