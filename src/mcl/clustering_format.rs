//! Binary clustering output formats from `clustering_format.cpp`.

use super::recursive_parser::{ParseError, RecursiveParser};
use crate::align::hsp::HspContext;
use crate::dp::swipe::HspValues;
use crate::output::format::OutputFlags;
use std::sync::Mutex;

/// Validated Source clustering expression shared by extension workers.
/// Evaluation errors are retained until the command boundary can return them.
#[derive(Debug)]
pub struct ClusteringFilter {
    pub format: ClusteringFormat,
    pub minimum: Option<f64>,
    invalid_output_class: bool,
    error: Mutex<Option<ParseError>>,
}

impl PartialEq for ClusteringFilter {
    fn eq(&self, other: &Self) -> bool {
        self.format == other.format
            && self.minimum == other.minimum
            && self.invalid_output_class == other.invalid_output_class
    }
}

impl ClusteringFilter {
    pub fn new(expression: &str, minimum: Option<f64>) -> Result<Self, ParseError> {
        Self::new_bytes(expression.as_bytes(), minimum)
    }
    pub fn new_bytes(expression: &[u8], minimum: Option<f64>) -> Result<Self, ParseError> {
        Ok(Self {
            format: ClusteringFormat::new_bytes(expression)?,
            minimum,
            invalid_output_class: false,
            error: Mutex::new(None),
        })
    }

    /// Source WITH_MCL dynamic_casts the selected output to Clustering_format
    /// at an actual threshold lookup. Other output classes dereference null;
    /// retain a deferred defined diagnostic without inventing an expression.
    pub(crate) fn invalid_output_class(minimum: f64) -> Self {
        Self {
            format: ClusteringFormat {
                format: Vec::new(),
                hsp_values: HspValues::NONE,
                flags: OutputFlags::NONE,
            },
            minimum: Some(minimum),
            invalid_output_class: true,
            error: Mutex::new(None),
        }
    }

    pub fn evaluate<T>(
        &self,
        context: &HspContext<T>,
        translated: bool,
    ) -> Result<f64, ParseError> {
        if self.invalid_output_class {
            return Err(ParseError::new("Source WITH_MCL cluster-threshold lookup dereferences a null Clustering_format cast for the selected output format."));
        }
        RecursiveParser::new_with_context_mode_bytes(Some(context), &self.format.format, translated)
            .evaluate()
    }

    pub fn accepts<T>(&self, context: &HspContext<T>, translated: bool) -> bool {
        let Some(minimum) = self.minimum else {
            return true;
        };
        // Without WITH_MCL the Source threshold evaluation is compiled out.
        if self.invalid_output_class && !cfg!(feature = "mcl") {
            return true;
        }

        match self.evaluate(context, translated) {
            Ok(value) => value >= minimum,
            Err(error) => {
                *self
                    .error
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner()) = Some(error);
                false
            }
        }
    }

    pub fn check(&self) -> Result<(), ParseError> {
        match self
            .error
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .clone()
        {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ClusteringFormat {
    pub format: Vec<u8>,
    pub hsp_values: HspValues,
    pub flags: OutputFlags,
}

impl ClusteringFormat {
    pub fn new(format: &str) -> Result<Self, ParseError> {
        Self::new_bytes(format.as_bytes())
    }
    pub fn new_bytes(format: &[u8]) -> Result<Self, ParseError> {
        let format = RecursiveParser::clean_expression_bytes(format);
        let mut parser = RecursiveParser::new_bytes(None, &format);
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

    pub fn print_match<T>(
        &self,
        context: &HspContext<T>,
        query_translated: bool,
        output: &mut Vec<u8>,
    ) -> Result<(), ParseError> {
        output.extend_from_slice(&(context.query_oid as u32).to_ne_bytes());
        output.extend_from_slice(&(context.subject_oid as u32).to_ne_bytes());
        let mut parser = RecursiveParser::new_with_context_mode_bytes(
            Some(context),
            &self.format,
            query_translated,
        );
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

    pub fn print_match<T>(context: &HspContext<T>, output: &mut Vec<u8>) {
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
    fn byte_titles_do_not_change_numeric_clustering_expression_or_record() {
        let text = context();
        let bytes = crate::align::hsp::ByteHspContext::new_bytes(
            text.hsp.clone(),
            text.query_id,
            text.query_oid,
            text.query.clone(),
            text.query_source_len,
            b"q\xff".to_vec(),
            text.subject_oid,
            text.subject_len,
            b"s\xfe".to_vec(),
            text.hit_num,
            text.hsp_num,
            text.subject_seq.clone(),
            text.query_self_aln_score,
            text.target_self_aln_score,
        );
        let format = ClusteringFormat::new("pident+normalized_bitscore_global").unwrap();
        let mut expected = Vec::new();
        let mut actual = Vec::new();
        format.print_match(&text, false, &mut expected).unwrap();
        format.print_match(&bytes, false, &mut actual).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(
            ClusteringFilter::new_bytes(&format.format, Some(150.0))
                .unwrap()
                .evaluate(&bytes, false)
                .unwrap(),
            150.0
        );
    }

    #[test]
    fn clustering_format_discovers_requirements_and_writes_native_record() {
        let format = ClusteringFormat::new(" pident + normalized_bitscore_global ").unwrap();
        assert_eq!(format.format, b"pident+normalized_bitscore_global");
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
