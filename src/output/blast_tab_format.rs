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
pub struct BlastTabConfig<T = String> {
    pub output_format: Vec<T>,
    pub json: bool,
    pub frame_shift: i32,
    pub query_translated: bool,
    /// Source `value_traits` is nucleotide for blastn, while query_translated
    /// remains false. Keep the alphabet separate from coordinate translation.
    pub nucleotide_query: bool,
    pub score_matrix_name: T,
    pub workflow: Workflow,
    /// `None` means the option was absent; `Some([])` requests the workflow
    /// default header, matching C++ `config.output_header.present()/empty()`.
    pub output_header: Option<Vec<T>>,
    pub version: String,
    pub invocation: Vec<u8>,
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
            nucleotide_query: false,
            score_matrix_name: "BLOSUM62".into(),
            workflow: Workflow::BlastP,
            output_header: None,
            version: env!("CARGO_PKG_VERSION").into(),
            invocation: Vec::new(),
            report_unaligned: 0,
            qnum_offset: 0,
            snum_offset: 0,
        }
    }
}

/// Configured BLAST tabular/JSON formatter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastTabFormat<T = String> {
    pub format: TabularFormat,
    pub config: BlastTabConfig<T>,
}

pub type ByteBlastTabConfig = BlastTabConfig<Vec<u8>>;
pub type ByteBlastTabFormat = BlastTabFormat<Vec<u8>>;

impl Default for ByteBlastTabConfig {
    fn default() -> Self {
        let config = BlastTabConfig::<String>::default();
        Self {
            output_format: config
                .output_format
                .into_iter()
                .map(String::into_bytes)
                .collect(),
            score_matrix_name: config.score_matrix_name.into_bytes(),
            output_header: config
                .output_header
                .map(|values| values.into_iter().map(String::into_bytes).collect()),
            json: config.json,
            frame_shift: config.frame_shift,
            query_translated: config.query_translated,
            nucleotide_query: config.nucleotide_query,
            workflow: config.workflow,
            version: config.version,
            invocation: config.invocation,
            report_unaligned: config.report_unaligned,
            qnum_offset: config.qnum_offset,
            snum_offset: config.snum_offset,
        }
    }
}
impl BlastTabFormat<String> {
    pub fn new(config: BlastTabConfig) -> Result<Self, String> {
        Self::new_core(config).map_err(|error| error.to_string())
    }
}
impl ByteBlastTabFormat {
    pub fn new_bytes(config: ByteBlastTabConfig) -> io::Result<Self> {
        Self::new_core(config)
    }
}
impl<N: AsRef<[u8]>> BlastTabFormat<N> {
    fn new_core(config: BlastTabConfig<N>) -> io::Result<Self> {
        let fields = config
            .output_format
            .iter()
            .map(AsRef::as_ref)
            .collect::<Vec<_>>();
        let format = TabularFormat::new_bytes(
            &fields,
            config.json,
            config.frame_shift,
            config.query_translated,
            config.score_matrix_name.as_ref(),
        )?;
        let header = config
            .output_header
            .as_ref()
            .map(|values| values.iter().map(AsRef::as_ref).collect::<Vec<_>>());
        TabularFormat::header_format_bytes(config.workflow, header.as_deref())?;
        Ok(Self { format, config })
    }

    pub fn header(&self) -> Result<String, String> {
        String::from_utf8(self.header_bytes().map_err(|error| error.to_string())?)
            .map_err(|error| error.to_string())
    }
    pub fn header_bytes(&self) -> io::Result<Vec<u8>> {
        let values = self
            .config
            .output_header
            .as_ref()
            .map(|values| values.iter().map(AsRef::as_ref).collect::<Vec<_>>());
        let header = TabularFormat::header_format_bytes(self.config.workflow, values.as_deref())?;
        super::format::print_tabular_header_bytes(
            &self.format.fields,
            header,
            self.config.json,
            &self.config.version,
            &self.config.invocation,
        )
        .map_err(io::Error::other)
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
        self.write_match_bytes(
            writer,
            context,
            subject_taxids,
            taxonomy,
            query_quality.as_bytes(),
            mate_sequence,
        )
    }

    pub fn write_match_bytes<W: Write, T: AsRef<[u8]>>(
        &self,
        writer: &mut W,
        context: &HspContext<T>,
        subject_taxids: &[TaxId],
        taxonomy: Option<&TaxonomyTree>,
        query_quality: &[u8],
        mate_sequence: Option<&[Letter]>,
    ) -> io::Result<()> {
        if self.config.json {
            return super::format::write_tabular_context_row_json_bytes(
                writer,
                context,
                context.hit_num != 0,
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
            if field.taxonomy_rank().is_some()
                || matches!(
                    field,
                    FieldId::STaxIds
                        | FieldId::SSciNames
                        | FieldId::SSKingdoms
                        | FieldId::SKingdoms
                        | FieldId::SPhylums
                        | FieldId::SLineages
                )
            {
                if let Some(tree) = taxonomy {
                    super::format::write_tabular_context_row_with_taxonomy_bytes(
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
                super::format::write_tabular_context_row_with_query_info_bytes(
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
                super::format::write_tabular_context_row_bytes(
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
        self.write_query_intro_bytes(
            writer,
            query_title.as_bytes(),
            query_len,
            query_sequence,
            query_quality.as_bytes(),
            unaligned,
        )
    }

    pub fn write_query_intro_bytes<W: Write>(
        &self,
        writer: &mut W,
        query_title: &[u8],
        query_len: i32,
        query_sequence: &[Letter],
        query_quality: &[u8],
        unaligned: bool,
    ) -> io::Result<()> {
        let enabled = unaligned && self.config.report_unaligned == 1;
        super::format::write_tabular_query_intro_with_quality_for_type_bytes(
            writer,
            query_title,
            query_len,
            query_sequence,
            query_quality,
            &self.format.fields,
            enabled,
            self.config.query_translated || self.config.nucleotide_query,
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
    fn byte_config_preserves_invocation_and_invalid_field_diagnostic() {
        let formatter = ByteBlastTabFormat::new_bytes(ByteBlastTabConfig {
            output_header: Some(vec![b"verbose".to_vec()]),
            invocation: b"diamond_\xff --view".to_vec(),
            ..Default::default()
        })
        .unwrap();
        assert!(formatter
            .header_bytes()
            .unwrap()
            .windows(b"diamond_\xff --view".len())
            .any(|part| part == b"diamond_\xff --view"));
        let error = ByteBlastTabFormat::new_bytes(ByteBlastTabConfig {
            output_format: vec![b"6".to_vec(), b"field_\xff".to_vec()],
            ..Default::default()
        })
        .unwrap_err();
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error),
            Some(b"Invalid output field: field_\xff".as_slice())
        );
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
            invocation: b"diamond test".to_vec(),
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
