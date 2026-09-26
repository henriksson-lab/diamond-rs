//! Compatibility facade mirroring `diamond/src/output/sam_format.cpp` and the
//! `SamFormat` declaration in `output_format.h`.
//!
//! Existing callers may continue to use [`super::sam`]; this module supplies
//! the upstream file boundary and its concrete format type.

pub use super::sam::{
    print_cigar, print_header, print_match_context, print_md, print_query_intro, write_sam_header,
    write_sam_record, SamRecord,
};

use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::output::format::OutputFormatSpec;
use crate::stats::score_matrix::ScoreMatrix;

/// Rust counterpart of C++ `SamFormat`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SamFormat;

impl SamFormat {
    /// Header-declared constructor state: SAM code, transcript DP fields,
    /// subject IDs, and default reporting of unaligned queries.
    pub fn output_format_spec(&self) -> OutputFormatSpec {
        OutputFormatSpec::sam_format()
    }

    /// Matches `SamFormat::print_query_intro` with `Output::Info` flattened to
    /// the two values this implementation consumes.
    pub fn print_query_intro<W: Write>(
        &self,
        writer: &mut W,
        query_title: &str,
        unaligned: bool,
    ) -> io::Result<()> {
        super::sam::print_query_intro(writer, query_title, unaligned)
    }

    /// Matches `SamFormat::print_match`. C++ derives one title-policy boolean
    /// from `salltitles || command == view`; accepting those inputs here avoids
    /// exposing combinations that the original method cannot produce.
    pub fn print_match<W: Write>(
        &self,
        writer: &mut W,
        record: &HspContext,
        score_matrix: &ScoreMatrix,
        salltitles: bool,
        view_command: bool,
        sam_qlen_field: bool,
    ) -> io::Result<()> {
        let long_titles = salltitles || view_command;
        super::sam::print_match_context(
            writer,
            record,
            score_matrix,
            long_titles,
            long_titles,
            sam_qlen_field,
        )
    }

    /// Matches `SamFormat::print_header`. Matrix, gap and first-query
    /// parameters from the virtual interface are intentionally absent because
    /// the C++ implementation does not read them.
    pub fn print_header<W: Write>(
        &self,
        writer: &mut W,
        mode: u32,
        invocation: &str,
    ) -> io::Result<()> {
        super::sam::print_header(writer, mode, invocation)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dp::swipe::HspValues;
    use crate::output::format::{OutputFlags, OutputFormatKind};

    #[test]
    fn sam_format_header_state_matches_cpp_constructor() {
        let spec = SamFormat.output_format_spec();
        assert_eq!(spec.code, OutputFormatKind::Sam);
        assert_eq!(spec.hsp_values, HspValues::TRANSCRIPT);
        assert_eq!(
            spec.flags,
            OutputFlags::SSEQID | OutputFlags::DEFAULT_REPORT_UNALIGNED
        );
    }

    #[test]
    fn sam_format_methods_preserve_free_function_output() {
        let mut via_type = Vec::new();
        SamFormat
            .print_query_intro(&mut via_type, "query title", true)
            .unwrap();
        let mut via_compatibility = Vec::new();
        print_query_intro(&mut via_compatibility, "query title", true).unwrap();
        assert_eq!(via_type, via_compatibility);

        via_type.clear();
        SamFormat
            .print_header(&mut via_type, 3, "diamond blastx")
            .unwrap();
        via_compatibility.clear();
        print_header(&mut via_compatibility, 3, "diamond blastx").unwrap();
        assert_eq!(via_type, via_compatibility);
    }
}
