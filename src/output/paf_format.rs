//! Mirrored facade for `diamond/src/output/paf_format.cpp`.

use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::output::format::OutputFormatSpec;
use crate::stats::score_matrix::ScoreMatrix;

pub use super::paf::{
    print_match_context, print_match_context_with_command, print_query_intro, write_paf_record,
    PafRecord,
};

/// Rust counterpart of C++ `PAFFormat`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PafFormat;

impl PafFormat {
    /// Header-declared constructor state: PAF code, transcript DP fields,
    /// subject IDs, and default reporting of unaligned queries.
    pub fn output_format_spec(&self) -> OutputFormatSpec {
        OutputFormatSpec::paf_format()
    }

    pub fn print_query_intro<W: Write>(
        &self,
        writer: &mut W,
        query_title: &str,
        unaligned: bool,
    ) -> io::Result<()> {
        super::paf::print_query_intro(writer, query_title, unaligned)
    }

    pub fn print_match<W: Write>(
        &self,
        writer: &mut W,
        record: &HspContext,
        score_matrix: &ScoreMatrix,
    ) -> io::Result<()> {
        super::paf::print_match_context(writer, record, score_matrix)
    }

    /// Full C++ mapping with the global blastn command state made explicit.
    pub fn print_match_with_command<W: Write>(
        &self,
        writer: &mut W,
        record: &HspContext,
        score_matrix: &ScoreMatrix,
        blastn_command: bool,
    ) -> io::Result<()> {
        super::paf::print_match_context_with_command(writer, record, score_matrix, blastn_command)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dp::swipe::HspValues;
    use crate::output::format::{OutputFlags, OutputFormatKind};

    #[test]
    fn paf_format_header_state_matches_cpp_constructor() {
        let spec = PafFormat.output_format_spec();
        assert_eq!(spec.code, OutputFormatKind::Paf);
        assert_eq!(spec.hsp_values, HspValues::TRANSCRIPT);
        assert_eq!(
            spec.flags,
            OutputFlags::SSEQID | OutputFlags::DEFAULT_REPORT_UNALIGNED
        );
    }

    #[test]
    fn unaligned_intro_delegates_exactly() {
        let mut direct = Vec::new();
        let mut facade = Vec::new();
        super::super::paf::print_query_intro(&mut direct, "query one", true).unwrap();
        PafFormat
            .print_query_intro(&mut facade, "query one", true)
            .unwrap();
        assert_eq!(facade, direct);
    }
}
