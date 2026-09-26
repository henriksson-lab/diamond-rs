//! Cluster realignment workflow from `diamond/src/cluster/realign.cpp`.
//!
//! The lower-level `Cluster::realign` overloads are implemented upstream in
//! `cluster/output.cpp`, not in `realign.cpp`.  [`RealignBackend`] keeps that
//! source boundary explicit while this module owns all workflow policy from
//! the mirrored source file.

use std::io::Write;

use crate::align::hsp::HspContext;
use crate::basic::value::OId;
use crate::data::sequence_file::SequenceFileFlags;
use crate::dp::swipe::HspValues;
use crate::output::format::{
    output_header, FormatCode, Header, OutputFlags, TabularFormat, Workflow,
};
use crate::util::data_structures::FlatArray;

/// C++ `DEFAULT_FORMAT`, including its leading output-format code.
pub const DEFAULT_REALIGN_FORMAT: [&str; 10] = [
    "6",
    "qseqid",
    "sseqid",
    "approx_pident",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
];

/// Explicit counterpart of the global C++ configuration consumed by
/// `Cluster::realign()`.
#[derive(Debug, Clone, PartialEq)]
pub struct RealignConfig {
    pub database: String,
    pub clustering: String,
    pub output_format: Vec<String>,
    pub output_header: Option<Vec<String>>,
    pub db_size: Option<u64>,
    pub max_evalue: f64,
    pub frame_shift: i32,
    pub query_translated: bool,
    pub score_matrix_name: String,
}

impl Default for RealignConfig {
    fn default() -> Self {
        Self {
            database: String::new(),
            clustering: String::new(),
            output_format: Vec::new(),
            output_header: None,
            db_size: None,
            max_evalue: 0.0,
            frame_shift: 0,
            query_translated: false,
            score_matrix_name: "BLOSUM62".to_owned(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DatabaseMetadata {
    pub sequence_count: u64,
    pub letters: u64,
    pub titles_lazy: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealignSummary {
    pub database: DatabaseMetadata,
    pub score_matrix_db_letters: u64,
    pub centroid_count: usize,
    pub mapping_count: usize,
}

/// Adapter for facilities whose upstream implementations live outside
/// `cluster/realign.cpp` (sequence-file construction, TSV cluster loading,
/// and the block-pair alignment overload in `cluster/output.cpp`).
pub trait RealignBackend {
    fn open_database(
        &mut self,
        path: &str,
        flags: SequenceFileFlags,
    ) -> Result<DatabaseMetadata, String>;

    fn set_score_matrix_db_letters(&mut self, letters: u64) -> Result<(), String>;

    fn init_random_access(&mut self) -> Result<(), String>;

    fn read_centroid_sorted(&mut self, path: &str) -> Result<(FlatArray<OId>, Vec<OId>), String>;

    fn realign_clusters(
        &mut self,
        clusters: &FlatArray<OId>,
        centroids: &[OId],
        hsp_values: HspValues,
        emit: &mut dyn FnMut(&HspContext) -> Result<(), String>,
    ) -> Result<(), String>;

    fn close_database(&mut self) -> Result<(), String>;
}

/// Execute C++ `Cluster::realign()` with explicit ownership adapters.
///
/// `format_match` corresponds to `OutputFormat::print_match`; it receives the
/// reused C++ `TextBuffer` equivalent and must append exactly one formatted
/// match. The buffer is written and cleared after every callback invocation.
pub fn realign<B, W, F>(
    config: &mut RealignConfig,
    backend: &mut B,
    output: &mut W,
    mut format_match: F,
) -> Result<RealignSummary, String>
where
    B: RealignBackend,
    W: Write,
    F: FnMut(&HspContext, &mut Vec<u8>) -> Result<(), String>,
{
    if config.database.is_empty() {
        return Err("Database file is required for the realign workflow.".to_owned());
    }
    if config.clustering.is_empty() {
        return Err("Clustering file is required for the realign workflow.".to_owned());
    }

    if config.output_format.is_empty() {
        config.output_format = DEFAULT_REALIGN_FORMAT
            .iter()
            .map(|value| (*value).to_owned())
            .collect();
    }
    let format_code = FormatCode::parse(&config.output_format[0])
        .ok_or_else(|| format!("Invalid output format: {}", config.output_format[0]))?;
    if format_code != FormatCode::Tabular {
        return Err("The realign workflow only supports tabular output format.".to_owned());
    }

    let format_args: Vec<&str> = config.output_format.iter().map(String::as_str).collect();
    let tabular = TabularFormat::new(
        &format_args,
        false,
        config.frame_shift,
        config.query_translated,
        &config.score_matrix_name,
    )?;
    for field in &tabular.fields {
        if field.output_flags().any(OutputFlags::NO_REALIGN) {
            let key = field
                .output_field()
                .map(|definition| definition.key)
                .unwrap_or("unknown");
            return Err(format!(
                "Unsupported output field for the realign workflow: {key}"
            ));
        }
    }

    let header_args = config
        .output_header
        .as_ref()
        .map(|values| values.iter().map(String::as_str).collect::<Vec<&str>>());
    let header = TabularFormat::header_format(Workflow::Cluster, header_args.as_deref())?;
    if header == Header::Simple {
        output
            .write_all(output_header(&tabular.fields, true)?.as_bytes())
            .map_err(|error| error.to_string())?;
    }

    let flags = SequenceFileFlags::NEED_LETTER_COUNT | SequenceFileFlags::ACC_TO_OID_MAPPING;
    let database = backend.open_database(&config.database, flags)?;
    let execution = (|| {
        let score_matrix_db_letters = config.db_size.unwrap_or(database.letters);
        backend.set_score_matrix_db_letters(score_matrix_db_letters)?;
        config.max_evalue = f64::MAX;
        if database.titles_lazy {
            backend.init_random_access()?;
        }
        let (clusters, centroids) = backend.read_centroid_sorted(&config.clustering)?;
        let mapping_count = usize::try_from(clusters.data_size())
            .map_err(|_| "Cluster mapping count does not fit usize.".to_owned())?;
        let mut buffer = Vec::new();
        let mut emit = |hsp: &HspContext| {
            format_match(hsp, &mut buffer)?;
            output
                .write_all(&buffer)
                .map_err(|error| error.to_string())?;
            buffer.clear();
            Ok(())
        };
        backend.realign_clusters(&clusters, &centroids, tabular.hsp_values, &mut emit)?;
        output.flush().map_err(|error| error.to_string())?;
        Ok(RealignSummary {
            database,
            score_matrix_db_letters,
            centroid_count: centroids.len(),
            mapping_count,
        })
    })();

    let close = backend.close_database();
    match (execution, close) {
        (Err(error), _) => Err(error),
        (Ok(_), Err(error)) => Err(error),
        (Ok(summary), Ok(())) => Ok(summary),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::format::FieldId;

    #[derive(Default)]
    struct FakeBackend {
        calls: Vec<String>,
        flags: Option<SequenceFileFlags>,
        hsp_values: Option<HspValues>,
        score_matrix_db_letters: Option<u64>,
        fail_realign: bool,
    }

    impl RealignBackend for FakeBackend {
        fn open_database(
            &mut self,
            path: &str,
            flags: SequenceFileFlags,
        ) -> Result<DatabaseMetadata, String> {
            self.calls.push(format!("open:{path}"));
            self.flags = Some(flags);
            Ok(DatabaseMetadata {
                sequence_count: 4,
                letters: 123,
                titles_lazy: true,
            })
        }

        fn set_score_matrix_db_letters(&mut self, letters: u64) -> Result<(), String> {
            self.calls.push(format!("db_letters:{letters}"));
            self.score_matrix_db_letters = Some(letters);
            Ok(())
        }

        fn init_random_access(&mut self) -> Result<(), String> {
            self.calls.push("random_access".to_owned());
            Ok(())
        }

        fn read_centroid_sorted(
            &mut self,
            path: &str,
        ) -> Result<(FlatArray<OId>, Vec<OId>), String> {
            self.calls.push(format!("read:{path}"));
            Ok((
                FlatArray::from_limits_data(vec![0_u64, 2, 3], vec![0, 1, 3]),
                vec![0, 3],
            ))
        }

        fn realign_clusters(
            &mut self,
            _: &FlatArray<OId>,
            _: &[OId],
            hsp_values: HspValues,
            emit: &mut dyn FnMut(&HspContext) -> Result<(), String>,
        ) -> Result<(), String> {
            self.calls.push("realign".to_owned());
            self.hsp_values = Some(hsp_values);
            if self.fail_realign {
                return Err("alignment failed".to_owned());
            }
            for title in ["first", "second"] {
                emit(&HspContext {
                    query_title: title.to_owned(),
                    ..HspContext::default()
                })?;
            }
            Ok(())
        }

        fn close_database(&mut self) -> Result<(), String> {
            self.calls.push("close".to_owned());
            Ok(())
        }
    }

    fn configured() -> RealignConfig {
        RealignConfig {
            database: "db.dmnd".to_owned(),
            clustering: "clusters.tsv".to_owned(),
            ..RealignConfig::default()
        }
    }

    #[test]
    fn default_workflow_preserves_order_flags_and_callback_output() {
        let mut config = configured();
        config.output_header = Some(Vec::new());
        let mut backend = FakeBackend::default();
        let mut output = Vec::new();
        let summary = realign(&mut config, &mut backend, &mut output, |hsp, buffer| {
            buffer.extend_from_slice(hsp.query_title.as_bytes());
            buffer.push(b'\n');
            Ok(())
        })
        .unwrap();

        assert_eq!(config.output_format, DEFAULT_REALIGN_FORMAT);
        assert_eq!(config.max_evalue, f64::MAX);
        assert_eq!(summary.score_matrix_db_letters, 123);
        assert_eq!((summary.centroid_count, summary.mapping_count), (2, 3));
        assert_eq!(
            backend.calls,
            [
                "open:db.dmnd",
                "db_letters:123",
                "random_access",
                "read:clusters.tsv",
                "realign",
                "close"
            ]
        );
        let flags = backend.flags.unwrap();
        assert!(flags.contains(SequenceFileFlags::NEED_LETTER_COUNT));
        assert!(flags.contains(SequenceFileFlags::ACC_TO_OID_MAPPING));
        let expected_hsp_values = [
            FieldId::QSeqId,
            FieldId::SSeqId,
            FieldId::ApproxPIdent,
            FieldId::QStart,
            FieldId::QEnd,
            FieldId::SStart,
            FieldId::SEnd,
            FieldId::EValue,
            FieldId::BitScore,
        ]
        .into_iter()
        .fold(HspValues::NONE, |values, field| values | field.hsp_values());
        assert_eq!(backend.hsp_values, Some(expected_hsp_values));
        assert_eq!(backend.score_matrix_db_letters, Some(123));
        assert_eq!(
            String::from_utf8(output).unwrap(),
            "cseqid\tmseqid\tapprox_pident\tcstart\tcend\tmstart\tmend\tevalue\tBitscore\nfirst\nsecond\n"
        );
    }

    #[test]
    fn configured_db_size_overrides_database_letters() {
        let mut config = configured();
        config.db_size = Some(999);
        let mut backend = FakeBackend::default();
        let summary = realign(&mut config, &mut backend, &mut Vec::new(), |_, _| Ok(())).unwrap();
        assert_eq!(summary.score_matrix_db_letters, 999);
        assert_eq!(backend.score_matrix_db_letters, Some(999));
    }

    #[test]
    fn rejects_non_tabular_and_no_realign_fields_verbatim() {
        let mut backend = FakeBackend::default();
        let mut config = configured();
        config.output_format = vec!["101".to_owned()];
        assert_eq!(
            realign(&mut config, &mut backend, &mut Vec::new(), |_, _| Ok(())).unwrap_err(),
            "The realign workflow only supports tabular output format."
        );

        config.output_format = vec!["6".to_owned(), "qseq".to_owned()];
        assert_eq!(
            realign(&mut config, &mut backend, &mut Vec::new(), |_, _| Ok(())).unwrap_err(),
            "Unsupported output field for the realign workflow: qseq"
        );
        assert!(backend.calls.is_empty());
    }

    #[test]
    fn closes_database_after_alignment_error() {
        let mut config = configured();
        let mut backend = FakeBackend {
            fail_realign: true,
            ..FakeBackend::default()
        };
        let error = realign(&mut config, &mut backend, &mut Vec::new(), |_, _| Ok(())).unwrap_err();
        assert_eq!(error, "alignment failed");
        assert_eq!(backend.calls.last().map(String::as_str), Some("close"));
    }

    #[test]
    fn requires_both_input_paths_before_side_effects() {
        let mut backend = FakeBackend::default();
        let mut config = RealignConfig::default();
        assert_eq!(
            realign(&mut config, &mut backend, &mut Vec::new(), |_, _| Ok(())).unwrap_err(),
            "Database file is required for the realign workflow."
        );
        config.database = "db.dmnd".to_owned();
        assert_eq!(
            realign(&mut config, &mut backend, &mut Vec::new(), |_, _| Ok(())).unwrap_err(),
            "Clustering file is required for the realign workflow."
        );
        assert!(backend.calls.is_empty());
    }
}
