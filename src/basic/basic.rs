//! Compatibility surface for `diamond/src/basic/basic.cpp`.
//!
//! Upstream collects unrelated out-of-line implementations and mutable process
//! globals in one translation unit. Rust keeps the owning types in their
//! established modules and re-exports that surface here. [`BasicState`] is the
//! explicit replacement for the mutable `align_mode`, `shapes`, `shape_from`,
//! and `shape_to` globals.

pub use super::consts::{PROGRAM_NAME, VERSION_STRING};
pub use super::reduction::Reduction;
pub use super::seed::Seed;
pub use super::sequence::Sequence;
pub use super::shape_config::ShapeConfig;
pub use super::statistics::{StatValue, Statistics, StatisticsReport, STATISTICS};
pub use super::translate::{translate_6_frames, translate_6_frames_with_genetic_code, Translator};
pub use super::value::{AlignMode, Letter, SequenceType, ValueTraits};

#[derive(Debug)]
pub struct BasicState {
    pub align_mode: AlignMode,
    pub statistics: Statistics,
    pub shapes: ShapeConfig,
    pub shape_from: u32,
    pub shape_to: u32,
}

impl Default for BasicState {
    fn default() -> Self {
        Self {
            align_mode: AlignMode::new(AlignMode::BLASTP),
            statistics: Statistics::new(),
            shapes: ShapeConfig::new(),
            shape_from: 0,
            shape_to: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::AMINO_ACID_ALPHABET;

    #[test]
    fn compatibility_surface_routes_basic_cpp_implementations() {
        assert_eq!(PROGRAM_NAME, "diamond");
        assert_eq!(VERSION_STRING, "2.1.24");
        assert_eq!(
            AlignMode::from_command(AlignMode::COMMAND_BLASTX),
            AlignMode::BLASTX
        );

        let sequence_data = [1, 2, 3];
        assert_eq!(Sequence::new(&sequence_data).reverse(), vec![3, 2, 1]);
        assert_eq!(Reduction::default_reduction().decode_seed(1, 3), "AAR");

        let traits = ValueTraits::new(AMINO_ACID_ALPHABET, 23, b"-U", SequenceType::AminoAcid);
        assert_eq!(
            Sequence::from_string("AR", &traits, 17).unwrap(),
            vec![0, 1]
        );
    }

    #[test]
    fn basic_state_replaces_mutable_process_globals() {
        let mut first = BasicState::default();
        let second = BasicState::default();
        first.align_mode = AlignMode::new(AlignMode::BLASTX);
        first.shape_from = 2;
        first.statistics.inc(StatValue::Matches, 3);

        assert_eq!(first.align_mode.query_contexts, 6);
        assert_eq!(first.shape_from, 2);
        assert_eq!(first.statistics.get(StatValue::Matches), 3);
        assert_eq!(second.align_mode.mode, AlignMode::BLASTP);
        assert_eq!(second.shape_from, 0);
        assert_eq!(second.statistics.get(StatValue::Matches), 0);
    }

    #[test]
    fn invalid_genetic_code_preserves_cpp_error() {
        assert_eq!(
            translate_6_frames_with_genetic_code(&[0, 3, 2], 7).unwrap_err(),
            "Invalid genetic code id."
        );
    }
}
