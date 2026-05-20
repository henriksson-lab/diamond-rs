use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::basic::value::OId;

/// Rust translation of C++ `Output::Format::Edge::Data`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeData {
    pub query: OId,
    pub target: OId,
    pub qcovhsp: f32,
    pub scovhsp: f32,
    pub evalue: f64,
}

impl EdgeData {
    pub const SIZE: usize = 32;

    pub fn from_context(r: &HspContext) -> Self {
        EdgeData {
            query: r.query_oid,
            target: r.subject_oid,
            qcovhsp: r.qcovhsp() as f32,
            scovhsp: r.scovhsp() as f32,
            evalue: r.corrected_bit_score(),
        }
    }

    pub fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.query.to_ne_bytes())?;
        writer.write_all(&self.target.to_ne_bytes())?;
        writer.write_all(&self.qcovhsp.to_ne_bytes())?;
        writer.write_all(&self.scovhsp.to_ne_bytes())?;
        writer.write_all(&self.evalue.to_ne_bytes())?;
        Ok(())
    }
}

/// Matches C++ `Output::Format::Edge::print_match(r)`.
pub fn print_match<W: Write>(writer: &mut W, r: &HspContext) -> io::Result<()> {
    EdgeData::from_context(r).write(writer)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::{Hsp, HspContext};
    use crate::util::interval::Interval;

    #[test]
    fn test_edge_data_from_context() {
        let mut hsp = Hsp::new();
        hsp.query_source_range = Interval::new(10, 30);
        hsp.subject_range = Interval::new(20, 50);
        hsp.corrected_bit_score = 47.25;
        let ctx = HspContext::new(
            hsp,
            0,
            11,
            Vec::new(),
            100,
            "",
            22,
            120,
            "",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );
        let data = EdgeData::from_context(&ctx);
        assert_eq!(data.query, 11);
        assert_eq!(data.target, 22);
        assert_eq!(data.qcovhsp, 20.0);
        assert_eq!(data.scovhsp, 25.0);
        assert_eq!(data.evalue, 47.25);
    }

    #[test]
    fn test_print_match_binary_layout() {
        let mut hsp = Hsp::new();
        hsp.query_source_range = Interval::new(0, 4);
        hsp.subject_range = Interval::new(2, 7);
        hsp.corrected_bit_score = 19.5;
        let ctx = HspContext::new(
            hsp,
            0,
            3,
            Vec::new(),
            8,
            "",
            5,
            10,
            "",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );
        let mut buf = Vec::new();
        print_match(&mut buf, &ctx).unwrap();
        assert_eq!(buf.len(), EdgeData::SIZE);
        assert_eq!(u64::from_ne_bytes(buf[0..8].try_into().unwrap()), 3);
        assert_eq!(u64::from_ne_bytes(buf[8..16].try_into().unwrap()), 5);
        assert_eq!(f32::from_ne_bytes(buf[16..20].try_into().unwrap()), 50.0);
        assert_eq!(f32::from_ne_bytes(buf[20..24].try_into().unwrap()), 50.0);
        assert_eq!(f64::from_ne_bytes(buf[24..32].try_into().unwrap()), 19.5);
    }
}
