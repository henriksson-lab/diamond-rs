//! Compatibility facade mirroring `diamond/src/output/xml_format.cpp` and
//! the C++ `XMLFormat` declaration in `output_format.h`.

use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::output::format::OutputFormatSpec;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence::AccessionParsing;

pub use super::xml::{
    print_footer, print_header, print_match_context, print_match_context_with_accession_stats,
    print_query_epilog, print_query_intro, write_hit, write_iteration_end, write_iteration_start,
    write_xml_footer, write_xml_header, XmlHit, XmlHsp,
};

/// Rust counterpart of C++ `XMLFormat`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct XmlFormat;

/// Descriptive alias matching the Rust module/file name.
pub type BlastXmlFormat = XmlFormat;

impl XmlFormat {
    /// Header-declared constructor state.
    pub fn output_format_spec(&self) -> OutputFormatSpec {
        OutputFormatSpec::xml_format()
    }

    pub fn print_match<W: Write>(
        &self,
        writer: &mut W,
        record: &HspContext,
        score_matrix: &ScoreMatrix,
        query_translated: bool,
        xml_blord_format: bool,
        no_parse_seqids: bool,
        accession_stats: &mut AccessionParsing,
    ) -> io::Result<()> {
        super::xml::print_match_context_with_accession_stats(
            writer,
            record,
            score_matrix,
            query_translated,
            xml_blord_format,
            no_parse_seqids,
            accession_stats,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn print_header<W: Write>(
        &self,
        writer: &mut W,
        mode: u32,
        matrix: &str,
        gap_open: i32,
        gap_extend: i32,
        evalue: f64,
        first_query_name: &str,
        first_query_len: u32,
        database: &str,
    ) -> io::Result<()> {
        super::xml::print_header(
            writer,
            mode,
            matrix,
            gap_open,
            gap_extend,
            evalue,
            first_query_name,
            first_query_len,
            database,
        )
    }

    pub fn print_query_intro<W: Write>(
        &self,
        writer: &mut W,
        query_oid: u64,
        query_title: &str,
        query_len: i32,
    ) -> io::Result<()> {
        super::xml::print_query_intro(writer, query_oid, query_title, query_len)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn print_query_epilog<W: Write>(
        &self,
        writer: &mut W,
        unaligned: bool,
        db_seqs: u64,
        db_letters: u64,
        kappa: f64,
        lambda: f64,
    ) -> io::Result<()> {
        super::xml::print_query_epilog(writer, unaligned, db_seqs, db_letters, kappa, lambda)
    }

    pub fn print_footer<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        super::xml::print_footer(writer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dp::swipe::HspValues;
    use crate::output::format::{OutputFlags, OutputFormatKind};

    #[test]
    fn xml_format_header_state_matches_cpp_constructor() {
        let spec = XmlFormat.output_format_spec();
        assert_eq!(spec.code, OutputFormatKind::BlastXml);
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
    fn facade_intro_and_footer_are_byte_identical() {
        let mut direct = Vec::new();
        print_query_intro(&mut direct, 4, "query <five>", 17).unwrap();
        let mut facade = Vec::new();
        XmlFormat
            .print_query_intro(&mut facade, 4, "query <five>", 17)
            .unwrap();
        assert_eq!(facade, direct);

        direct.clear();
        facade.clear();
        print_footer(&mut direct).unwrap();
        XmlFormat.print_footer(&mut facade).unwrap();
        assert_eq!(facade, direct);
    }
}
