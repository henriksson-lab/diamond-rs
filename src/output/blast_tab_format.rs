//! Source-mirrored facade for `diamond/src/output/blast_tab_format.cpp`.
//!
//! Field definitions and callbacks historically lived in [`super::format`].
//! They remain the single implementation and are re-exported here so existing
//! callers keep their API while new code follows the upstream source boundary.

use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::basic::value::{Letter, TaxId};
use crate::data::taxonomy::TaxonomyTree;

pub use super::format::{
    output_header, parse_tabular_fields, print_full_qqual, print_full_qseq_mate, print_lineage,
    print_qqual, print_rank_taxon_names, print_staxids, print_tabular_footer, print_tabular_header,
    print_taxon_names, print_title, print_title_escaped, FieldId, Header, OutputField,
    TabularFormat, Workflow, DEFAULT_TABULAR_FIELDS,
};

/// All former process-global inputs used by C++ `TabularFormat`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastTabConfig {
    pub output_format: Vec<String>,
    pub json: bool,
    pub frame_shift: i32,
    pub query_translated: bool,
    pub score_matrix_name: String,
    pub workflow: Workflow,
    /// `None` means the option was absent; `Some([])` requests the workflow
    /// default header, matching C++ `config.output_header.present()/empty()`.
    pub output_header: Option<Vec<String>>,
    pub version: String,
    pub invocation: String,
    pub report_unaligned: i32,
    pub qnum_offset: u64,
    pub snum_offset: u64,
}

impl Default for BlastTabConfig {
    fn default() -> Self {
        Self {
            output_format: vec!["6".into()],
            json: false,
            frame_shift: 0,
            query_translated: false,
            score_matrix_name: "BLOSUM62".into(),
            workflow: Workflow::BlastP,
            output_header: None,
            version: env!("CARGO_PKG_VERSION").into(),
            invocation: String::new(),
            report_unaligned: 0,
            qnum_offset: 0,
            snum_offset: 0,
        }
    }
}

/// Configured BLAST tabular/JSON formatter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastTabFormat {
    pub format: TabularFormat,
    pub config: BlastTabConfig,
}

impl BlastTabFormat {
    pub fn new(config: BlastTabConfig) -> Result<Self, String> {
        let fields = config
            .output_format
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        let format = TabularFormat::new(
            &fields,
            config.json,
            config.frame_shift,
            config.query_translated,
            &config.score_matrix_name,
        )?;
        // Validate header configuration at construction, as the C++ header
        // path otherwise defers this error until output starts.
        let header = config
            .output_header
            .as_ref()
            .map(|values| values.iter().map(String::as_str).collect::<Vec<_>>());
        TabularFormat::header_format(config.workflow, header.as_deref())?;
        Ok(Self { format, config })
    }

    pub fn header(&self) -> Result<String, String> {
        let values = self
            .config
            .output_header
            .as_ref()
            .map(|values| values.iter().map(String::as_str).collect::<Vec<_>>());
        let header = TabularFormat::header_format(self.config.workflow, values.as_deref())?;
        print_tabular_header(
            &self.format.fields,
            header,
            self.config.json,
            &self.config.version,
            &self.config.invocation,
        )
    }

    pub fn footer(&self) -> &'static str {
        print_tabular_footer(self.config.json)
    }

    pub fn field_header(&self, cluster: bool) -> Result<String, String> {
        output_header(&self.format.fields, cluster)
    }

    /// C++ `print_match`, with taxonomy and per-query data supplied explicitly.
    #[allow(clippy::too_many_arguments)]
    pub fn write_match<W: Write>(
        &self,
        writer: &mut W,
        context: &HspContext,
        subject_taxids: &[TaxId],
        taxonomy: Option<&TaxonomyTree>,
        query_quality: &str,
        mate_sequence: Option<&[Letter]>,
    ) -> io::Result<()> {
        if self.config.json {
            return super::format::write_tabular_context_row_json(
                writer,
                context,
                &self.format.fields,
                self.config.query_translated,
                self.config.frame_shift != 0,
                self.config.qnum_offset,
                self.config.snum_offset,
                subject_taxids,
                taxonomy,
                query_quality,
                mate_sequence,
            );
        }

        for (index, field) in self.format.fields.iter().enumerate() {
            if index > 0 {
                writer.write_all(b"\t")?;
            }
            let mut cell = Vec::new();
            if matches!(
                field,
                FieldId::STaxIds
                    | FieldId::SSciNames
                    | FieldId::SSKingdoms
                    | FieldId::SKingdoms
                    | FieldId::SPhylums
                    | FieldId::SLineages
            ) {
                if let Some(tree) = taxonomy {
                    super::format::write_tabular_context_row_with_taxonomy(
                        &mut cell,
                        context,
                        std::slice::from_ref(field),
                        self.config.query_translated,
                        self.config.frame_shift != 0,
                        self.config.qnum_offset,
                        self.config.snum_offset,
                        subject_taxids,
                        tree,
                        false,
                    )?;
                } else if *field == FieldId::STaxIds {
                    write!(cell, "{}\n", print_staxids(subject_taxids, false))?;
                } else {
                    cell.extend_from_slice(b"N/A\n");
                }
            } else if matches!(
                field,
                FieldId::QQual | FieldId::FullQQual | FieldId::FullQSeqMate
            ) {
                super::format::write_tabular_context_row_with_query_info(
                    &mut cell,
                    context,
                    std::slice::from_ref(field),
                    self.config.query_translated,
                    self.config.frame_shift != 0,
                    self.config.qnum_offset,
                    self.config.snum_offset,
                    query_quality,
                    mate_sequence,
                )?;
            } else {
                super::format::write_tabular_context_row(
                    &mut cell,
                    context,
                    std::slice::from_ref(field),
                    self.config.query_translated,
                    self.config.frame_shift != 0,
                    self.config.qnum_offset,
                    self.config.snum_offset,
                )?;
            }
            if cell.last() == Some(&b'\n') {
                cell.pop();
            }
            writer.write_all(&cell)?;
        }
        writer.write_all(b"\n")
    }

    /// C++ emits unaligned query rows only for the exact value `1` here.
    pub fn write_query_intro<W: Write>(
        &self,
        writer: &mut W,
        query_title: &str,
        query_len: i32,
        query_sequence: &[Letter],
        query_quality: &str,
        unaligned: bool,
    ) -> io::Result<()> {
        let enabled = unaligned && self.config.report_unaligned == 1;
        if enabled {
            // These C++ callbacks intentionally define match output only. The
            // inherited invalid query-intro handler throws instead of printing
            // the generic numeric `-1` sentinel.
            if let Some(field) = self.format.fields.iter().find(|field| {
                matches!(
                    field,
                    FieldId::HspNum
                        | FieldId::NormalizedBitscore
                        | FieldId::NormalizedNident
                        | FieldId::CorrectedBitScore
                        | FieldId::Reserved1
                        | FieldId::Reserved2
                )
            }) {
                let key = field.output_field().map_or("reserved", |def| def.key);
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Invalid output field: {key}"),
                ));
            }
        }
        super::format::write_tabular_query_intro_with_quality(
            writer,
            query_title,
            query_len,
            query_sequence,
            query_quality,
            &self.format.fields,
            enabled,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::data::taxonomy::TaxonomyNode;
    use crate::util::interval::Interval;

    fn context() -> HspContext {
        let mut hsp = Hsp::new();
        hsp.score = 42;
        hsp.query_source_range = Interval::new(1, 4);
        hsp.subject_range = Interval::new(0, 3);
        hsp.subject_source_range = Interval::new(0, 3);
        HspContext::new(
            hsp,
            0,
            7,
            vec![vec![0, 1, 2, 3, 4]],
            5,
            "query title",
            11,
            3,
            "subject title\x01second title",
            0,
            0,
            vec![0, 1, 2],
            0.0,
            0.0,
        )
    }

    fn taxonomy() -> TaxonomyTree {
        let mut tree = TaxonomyTree::new();
        tree.add_node(TaxonomyNode {
            taxid: 1,
            parent: 1,
            rank: "no rank".into(),
            name: "root".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 2,
            parent: 1,
            rank: "superkingdom".into(),
            name: "Bacteria".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 20,
            parent: 2,
            rank: "phylum".into(),
            name: "Proteobacteria".into(),
        });
        tree
    }

    #[test]
    fn tabular_facade_is_byte_exact_for_fields_and_taxonomy() {
        let formatter = BlastTabFormat::new(BlastTabConfig {
            output_format: [
                "6",
                "qseqid",
                "sallseqid",
                "staxids",
                "sscinames",
                "qqual",
                "score",
            ]
            .into_iter()
            .map(str::to_string)
            .collect(),
            ..Default::default()
        })
        .unwrap();
        let mut output = Vec::new();
        formatter
            .write_match(
                &mut output,
                &context(),
                &[20],
                Some(&taxonomy()),
                "abcdef",
                None,
            )
            .unwrap();
        assert_eq!(
            output,
            b"query\tsubject;second\t20\tProteobacteria\tbcd\t42\n"
        );
    }

    #[test]
    fn json_header_match_footer_are_byte_exact() {
        let formatter = BlastTabFormat::new(BlastTabConfig {
            output_format: ["104", "qseqid", "sallseqid", "staxids", "slineages"]
                .into_iter()
                .map(str::to_string)
                .collect(),
            json: true,
            ..Default::default()
        })
        .unwrap();
        let mut output = formatter.header().unwrap().into_bytes();
        formatter
            .write_match(&mut output, &context(), &[20], Some(&taxonomy()), "", None)
            .unwrap();
        output.extend_from_slice(formatter.footer().as_bytes());
        assert_eq!(String::from_utf8(output).unwrap(), "[\n\t{\n\t\"qseqid\":\"query\",\n\t\"sallseqid\":[\"subject\",\"second\"],\n\t\"staxids\":[20],\n\t\"slineages\": [\n\t\t[\"Bacteria\", \"Proteobacteria\"]\n\t]\n\t}\n]");
    }

    #[test]
    fn unaligned_intro_and_verbose_header_use_explicit_config() {
        let formatter = BlastTabFormat::new(BlastTabConfig {
            output_format: ["6", "qseqid", "qlen", "sseqid", "full_qqual"]
                .into_iter()
                .map(str::to_string)
                .collect(),
            output_header: Some(Vec::new()),
            version: "9.9".into(),
            invocation: "diamond test".into(),
            report_unaligned: 1,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(formatter.header().unwrap(), "# DIAMOND v9.9. http://github.com/bbuchfink/diamond\n# Invocation: diamond test\n# Fields: Query Seq - id, Query sequence length, Subject Seq - id, Query quality values\n");
        let mut output = Vec::new();
        formatter
            .write_query_intro(&mut output, "query title", 4, &[0, 1, 2, 3], "!!!!", true)
            .unwrap();
        assert_eq!(output, b"query\t4\t*\t!!!!\n");

        let invalid = BlastTabFormat::new(BlastTabConfig {
            output_format: ["6", "hspnum"].into_iter().map(str::to_string).collect(),
            report_unaligned: 1,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(
            invalid
                .write_query_intro(&mut Vec::new(), "q", 1, &[0], "", true)
                .unwrap_err()
                .to_string(),
            "Invalid output field: hspnum"
        );
    }
}
