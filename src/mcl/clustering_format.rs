//! Binary clustering output formats from `clustering_format.cpp`.

use super::recursive_parser::{ParseError, RecursiveParser};
use crate::align::hsp::HspContext;
use crate::dp::swipe::HspValues;
use crate::output::format::OutputFlags;

#[derive(Debug, Clone, PartialEq)]
pub struct ClusteringFormat {
    pub format: String,
    pub hsp_values: HspValues,
    pub flags: OutputFlags,
}

impl ClusteringFormat {
    pub fn new(format: &str) -> Result<Self, ParseError> {
        let format = RecursiveParser::clean_expression(format);
        let mut parser = RecursiveParser::new(None, &format);
        parser.evaluate()?;
        let mut hsp_values = HspValues::NONE;
        let mut flags = OutputFlags::NONE;
        for variable in parser.variables() {
            hsp_values = hsp_values | variable.hsp_values;
            flags |= variable.flags;
        }
        Ok(Self {
            format,
            hsp_values,
            flags,
        })
    }

    pub fn print_match(
        &self,
        context: &HspContext,
        query_translated: bool,
        output: &mut Vec<u8>,
    ) -> Result<(), ParseError> {
        output.extend_from_slice(&(context.query_oid as u32).to_ne_bytes());
        output.extend_from_slice(&(context.subject_oid as u32).to_ne_bytes());
        let mut parser =
            RecursiveParser::new_with_mode(Some(context), &self.format, query_translated);
        output.extend_from_slice(&parser.evaluate()?.to_ne_bytes());
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Bin1Format;

impl Bin1Format {
    pub const HSP_VALUES: HspValues = HspValues::TRANSCRIPT;

    pub fn print_query_intro(block_id: u32, output: &mut Vec<u8>) {
        output.extend_from_slice(&u32::MAX.to_ne_bytes());
        output.extend_from_slice(&block_id.to_ne_bytes());
    }

    pub fn print_match(context: &HspContext, output: &mut Vec<u8>) {
        if u64::from(context.query_id) < context.subject_oid {
            output.extend_from_slice(&(context.subject_oid as u32).to_ne_bytes());
            let denominator = context.query_source_len.max(context.subject_len) as f64;
            output.extend_from_slice(&(context.bit_score() / denominator).to_ne_bytes());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;

    fn context() -> HspContext {
        let mut hsp = Hsp::default();
        hsp.bit_score = 80.0;
        hsp.length = 10;
        hsp.identities = 7;
        HspContext::new(
            hsp,
            2,
            7,
            vec![vec![0; 20]],
            20,
            "q",
            9,
            40,
            "s",
            0,
            0,
            vec![0; 40],
            100.0,
            80.0,
        )
    }

    #[test]
    fn clustering_format_discovers_requirements_and_writes_native_record() {
        let format = ClusteringFormat::new(" pident + normalized_bitscore_global ").unwrap();
        assert_eq!(format.format, "pident+normalized_bitscore_global");
        assert!(format.hsp_values.0 != 0);
        assert!(format.flags.any(OutputFlags::SELF_ALN_SCORES));

        let mut output = Vec::new();
        format.print_match(&context(), false, &mut output).unwrap();
        assert_eq!(u32::from_ne_bytes(output[0..4].try_into().unwrap()), 7);
        assert_eq!(u32::from_ne_bytes(output[4..8].try_into().unwrap()), 9);
        assert_eq!(f64::from_ne_bytes(output[8..16].try_into().unwrap()), 150.0);
    }

    #[test]
    fn bin1_writes_intro_and_only_upper_triangle_matches() {
        let mut output = Vec::new();
        Bin1Format::print_query_intro(12, &mut output);
        Bin1Format::print_match(&context(), &mut output);
        assert_eq!(
            u32::from_ne_bytes(output[0..4].try_into().unwrap()),
            u32::MAX
        );
        assert_eq!(u32::from_ne_bytes(output[4..8].try_into().unwrap()), 12);
        assert_eq!(u32::from_ne_bytes(output[8..12].try_into().unwrap()), 9);
        assert_eq!(f64::from_ne_bytes(output[12..20].try_into().unwrap()), 2.0);

        let mut reverse = context();
        reverse.query_id = 10;
        Bin1Format::print_match(&reverse, &mut output);
        assert_eq!(output.len(), 20);
    }
}
