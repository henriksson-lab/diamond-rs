//! Hierarchy-preserving facade for `diamond/src/align/extend.cpp`.
//!
//! The implementation was established in `align::target`, `align::hsp`, and
//! `align::gapped_filter`; this module exposes the source translation-unit
//! surface without duplicating its alignment machinery.

pub use super::gapped_filter::{seed_only_hsp, seed_only_matches};
pub use super::hsp::Match;
pub use super::target::{
    add_self_aln, extend, extend_chunk, extend_mut, extend_seed_hit_list, extend_seed_hit_list_mut,
    first_round_hspv, have_filters, lazy_masking, ranking_chunk_size, ranking_terminate,
    GappedScoreConfig,
};
pub use crate::search::sensitivity::ExtensionMode as Mode;

/// Source-order table behind C++ `Extension::default_ext_mode`.
pub const DEFAULT_EXT_MODE: &[(crate::config::Sensitivity, Mode)] = &[
    (crate::config::Sensitivity::Faster, Mode::BandedFast),
    (crate::config::Sensitivity::Fast, Mode::BandedFast),
    (crate::config::Sensitivity::Shapes6x10, Mode::BandedFast),
    (crate::config::Sensitivity::Shapes30x10, Mode::BandedFast),
    (crate::config::Sensitivity::Linclust20, Mode::BandedFast),
    (crate::config::Sensitivity::Linclust40, Mode::BandedFast),
    (crate::config::Sensitivity::Default, Mode::BandedFast),
    (crate::config::Sensitivity::MidSensitive, Mode::BandedFast),
    (crate::config::Sensitivity::Sensitive, Mode::BandedFast),
    (crate::config::Sensitivity::MoreSensitive, Mode::BandedSlow),
    (crate::config::Sensitivity::VerySensitive, Mode::BandedSlow),
    (crate::config::Sensitivity::UltraSensitive, Mode::BandedSlow),
];

pub fn default_ext_mode(sensitivity: crate::config::Sensitivity) -> Option<Mode> {
    DEFAULT_EXT_MODE
        .iter()
        .find_map(|&(key, value)| (key == sensitivity).then_some(value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Sensitivity;

    #[test]
    fn extension_mode_table_matches_cpp() {
        assert_eq!(DEFAULT_EXT_MODE.len(), 12);
        assert_eq!(
            default_ext_mode(Sensitivity::Faster),
            Some(Mode::BandedFast)
        );
        assert_eq!(
            default_ext_mode(Sensitivity::MoreSensitive),
            Some(Mode::BandedSlow)
        );
        assert_eq!(
            default_ext_mode(Sensitivity::UltraSensitive),
            Some(Mode::BandedSlow)
        );
    }
}
