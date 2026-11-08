//! Compatibility facade mirroring `diamond/src/output/output_format.{h,cpp}`.
//!
//! The Rust implementation keeps the substantial tabular-format machinery in
//! [`super::format`], binary intermediate records in [`super::intermediate`],
//! and the binary edge format in [`super::edge`]. Re-exporting them here keeps
//! the upstream source-file boundary available without breaking established
//! Rust module paths.

pub use super::edge::{print_match as print_edge_match, Edge, EdgeData};
pub use super::format::*;
pub use super::intermediate::{IntermediateRecord, OutputFormat, FINISHED};

use crate::align::hsp::Hsp as AlignHsp;
use crate::basic::value::Letter;

/// C++ `Null_format`: the ordinary transcript pipeline with inherited no-op
/// output callbacks. Cloning preserves its stateless constructor state.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NullFormat;

impl NullFormat {
    pub fn output_format_spec(&self) -> OutputFormatSpec {
        OutputFormatSpec::null_format()
    }

    pub fn print_match<W: std::io::Write>(
        &self,
        _writer: &mut W,
        _context: &crate::align::hsp::HspContext,
    ) -> std::io::Result<()> {
        Ok(())
    }
    pub fn print_match_bytes<W: std::io::Write>(
        &self,
        _writer: &mut W,
        _context: &crate::align::hsp::ByteHspContext,
    ) -> std::io::Result<()> {
        Ok(())
    }
}

/// Matches C++ `print_hsp` as it exists in `output_format.cpp`.
///
/// The upstream function's pairwise rendering and stdout write are commented
/// out; its only observable behavior is therefore a no-op.
pub fn print_hsp(_hsp: &mut AlignHsp, _query: &[Letter]) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn print_hsp_matches_upstream_noop() {
        let mut hsp = AlignHsp::new();
        hsp.score = 42;
        print_hsp(&mut hsp, &[0, 1, 2]);
        assert_eq!(hsp.score, 42);
    }

    #[test]
    fn compatibility_exports_are_available() {
        let format = OutputFormat::edge();
        assert_eq!(format.hsp_values, crate::dp::swipe::HspValues::COORDS);
        assert_eq!(OutputFormatSpec::edge_format().code, OutputFormatKind::Edge);
        let _edge = Edge;
    }
}
