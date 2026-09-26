//! Compatibility facade for `diamond/src/output/blast_pairwise_format.cpp`
//! and the `PairwiseFormat` declaration in `output_format.h`.

use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::output::format::OutputFormatSpec;
use crate::stats::score_matrix::ScoreMatrix;

pub use super::pairwise::{
    print_footer, print_header, print_match_context, print_query_epilog, print_query_intro,
    write_header, write_pairwise, write_query_header, PairwiseData,
};

/// Rust counterpart of C++ `PairwiseFormat`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PairwiseFormat;

impl PairwiseFormat {
    /// Header-declared constructor state: pairwise code, transcript fields,
    /// full/all subject titles, and default reporting of unaligned queries.
    pub fn output_format_spec(&self) -> OutputFormatSpec {
        OutputFormatSpec::pairwise_format()
    }

    /// C++ `PairwiseFormat::print_match` with process globals explicit.
    pub fn print_match<W: Write>(
        &self,
        writer: &mut W,
        record: &HspContext,
        score_matrix: &ScoreMatrix,
        query_translated: bool,
        blastn_command: bool,
    ) -> io::Result<()> {
        super::pairwise::print_match_context(
            writer,
            record,
            score_matrix,
            query_translated,
            blastn_command,
        )
    }

    /// C++ `PairwiseFormat::print_header`. Its virtual-interface metadata
    /// arguments are omitted because the implementation ignores every one.
    pub fn print_header<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        super::pairwise::print_header(writer)
    }

    /// C++ `PairwiseFormat::print_query_intro`.
    pub fn print_query_intro<W: Write>(
        &self,
        writer: &mut W,
        query_title: &str,
        query_len: i32,
        unaligned: bool,
    ) -> io::Result<()> {
        super::pairwise::print_query_intro(writer, query_title, query_len, unaligned)
    }

    /// C++ `PairwiseFormat::print_query_epilog` (intentionally empty).
    pub fn print_query_epilog<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        super::pairwise::print_query_epilog(writer)
    }

    /// C++ `PairwiseFormat::print_footer` (intentionally empty).
    pub fn print_footer<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        super::pairwise::print_footer(writer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::{Hsp, HspContext};
    use crate::basic::packed_transcript::EditOperation;
    use crate::dp::swipe::HspValues;
    use crate::output::format::{OutputFlags, OutputFormatKind};
    use crate::util::interval::Interval;

    fn score_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1000).unwrap()
    }

    fn blastn_context() -> HspContext {
        let mut hsp = Hsp::new();
        hsp.score = 20;
        hsp.bit_score = 10.0;
        hsp.evalue = 0.0;
        hsp.frame = 0;
        hsp.identities = 2;
        hsp.positives = 2;
        hsp.length = 2;
        hsp.query_range = Interval::new(2, 4);
        hsp.query_source_range = Interval::new(6, 8);
        hsp.subject_range = Interval::new(4, 6);
        hsp.transcript.push_with_count(EditOperation::Match, 2);
        HspContext::new(
            hsp,
            0,
            0,
            vec![vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1]],
            10,
            "query",
            0,
            10,
            "subject",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        )
    }

    #[test]
    fn constructor_state_matches_output_format_header() {
        let spec = PairwiseFormat.output_format_spec();
        assert_eq!(spec.code, OutputFormatKind::BlastPairwise);
        assert_eq!(spec.hsp_values, HspValues::TRANSCRIPT);
        assert_eq!(
            spec.flags,
            OutputFlags::FULL_TITLES
                | OutputFlags::SSEQID
                | OutputFlags::DEFAULT_REPORT_UNALIGNED
                | OutputFlags::ALL_SEQIDS
        );
    }

    #[test]
    fn header_and_unaligned_intro_are_byte_exact() {
        let mut output = Vec::new();
        PairwiseFormat.print_header(&mut output).unwrap();
        PairwiseFormat
            .print_query_intro(&mut output, "query one", 12, true)
            .unwrap();
        PairwiseFormat.print_query_epilog(&mut output).unwrap();
        PairwiseFormat.print_footer(&mut output).unwrap();
        assert_eq!(
            output,
            b"BLASTP 2.3.0+\n\n\nQuery= query one\n\nLength=12\n\n\n***** No hits found *****\n\n\n"
        );
    }

    #[test]
    fn blastn_match_uses_reverse_absolute_query_coordinates() {
        let mut output = Vec::new();
        PairwiseFormat
            .print_match(&mut output, &blastn_context(), &score_matrix(), false, true)
            .unwrap();
        let output = String::from_utf8(output).unwrap();
        assert_eq!(
            output,
            concat!(
                ">subject\nLength=10\n\n",
                " Score = 10.0 bits (20),  Expect = 0.0\n",
                " Identities = 2/2 (100%), Positives = 2/2 (100%), Gaps = 0/2 (0%)\n",
                " Strand = Minus/Plus\n\n",
                "Query  8  ND 7\n",
                "          ND\n",
                "Sbjct  5  ND 6\n\n",
            )
        );
    }

    #[test]
    fn facade_match_is_byte_identical_to_compatibility_function() {
        let context = blastn_context();
        let matrix = score_matrix();
        let mut facade = Vec::new();
        PairwiseFormat
            .print_match(&mut facade, &context, &matrix, false, true)
            .unwrap();
        let mut direct = Vec::new();
        super::super::pairwise::print_match_context(&mut direct, &context, &matrix, false, true)
            .unwrap();
        assert_eq!(facade, direct);
    }
}
