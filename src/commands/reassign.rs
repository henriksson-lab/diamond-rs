//! Concrete database, search and output adapters for `Cluster::reassign`.
use std::cell::RefCell;
use std::ffi::OsString;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use super::blastp::{self, BlastpConfig};
use super::cluster_database::NativeDatabase;
use crate::basic::value::OId;
use crate::cluster::helpers;
use crate::cluster::realign::DatabaseMetadata;
use crate::cluster::reassign::{
    self, ReassignBackend, ReassignConfig, ReassignSearchConfig, SubDatabase,
};
use crate::data::fasta::ByteFastaRecord as FastaRecord;
use crate::data::sequence_file::SequenceFileFlags;
use crate::masking::{MaskingAlgo, MaskingMode};
use crate::output::format::{Header, TabularFormat, Workflow};
use crate::util::command_line_parser::CommandLineParser;
use crate::util::options::OptionValue;
use crate::util::text_buffer::TextBuffer;

/// Actual Source OutputFile owner behind workflow std::io::Write callbacks.
/// Explicit close follows the CLI workflow; public generic APIs still receive
/// their caller-owned writers and do not construct this adapter.
pub(super) struct SourceWorkflowOutput(crate::util::io::OutputFile);
impl SourceWorkflowOutput {
    pub(super) fn new(path: &str, buffer_size: usize) -> Result<Self, String> {
        Self::new_io(path, buffer_size).map_err(|error| error.to_string())
    }
    pub(super) fn new_io(path: &str, buffer_size: usize) -> io::Result<Self> {
        Self::new_io_path(Path::new(path), buffer_size)
    }
    pub(super) fn new_io_path(path: &Path, buffer_size: usize) -> io::Result<Self> {
        crate::util::io::OutputFile::with_source_profile_path(
            path,
            crate::util::io::Compressor::None,
            "wb",
            buffer_size,
        )
        .map(Self)
        .map_err(io::Error::other)
    }
    pub(super) fn close(&mut self) -> Result<(), String> {
        self.0.close().map_err(|error| error.to_string())
    }
}
impl Write for SourceWorkflowOutput {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.write_raw(bytes).map_err(io::Error::other)?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        self.0.flush().map_err(io::Error::other)
    }
}

struct Backend {
    errors: crate::util::io::IoErrorOwner,
    database: Option<NativeDatabase>,
    members: Vec<FastaRecord>,
    centroids: Vec<FastaRecord>,
    destination: Option<SourceWorkflowOutput>,
    output: PathBuf,
    simple_header: bool,
    oid_output: bool,
    search: BlastpConfig,
    extension: Option<crate::search::sensitivity::ExtensionMode>,
    tantan_min_mask_prob: f64,
}

impl ReassignBackend for Backend {
    fn io_errors(&self) -> crate::util::io::IoErrorOwner {
        self.errors.clone()
    }

    fn open_database(
        &mut self,
        path: &str,
        flags: SequenceFileFlags,
    ) -> Result<DatabaseMetadata, String> {
        self.open_database_path(Path::new(path), flags)
    }
    fn read_clustering(&mut self, path: &str) -> Result<Vec<OId>, String> {
        self.read_clustering_path(Path::new(path))
    }

    fn open_database_path(
        &mut self,
        path: &Path,
        flags: SequenceFileFlags,
    ) -> Result<DatabaseMetadata, String> {
        let database = self
            .errors
            .capture(NativeDatabase::open_source_with_buffer_owned_path(
                path,
                flags,
                Some(self.search.file_buffer_size),
                &self.search.tmpdir,
                self.search.tuning.no_unlink,
                &self.errors,
            ))?;
        let metadata = DatabaseMetadata {
            sequence_count: database.count,
            letters: database.letters,
            titles_lazy: false,
        };
        self.search.tuning.db_size = database.letters;
        self.database = Some(database);
        Ok(metadata)
    }
    fn open_output(&mut self) -> Result<(), String> {
        let mut output = self.errors.capture(SourceWorkflowOutput::new_io_path(
            &self.output,
            self.search.file_buffer_size,
        ))?;
        if self.simple_header {
            writeln!(output, "centroid\tmember").map_err(|e| e.to_string())?;
        }
        self.destination = Some(output);
        Ok(())
    }
    fn read_clustering_path(&mut self, path: &Path) -> Result<Vec<OId>, String> {
        helpers::read_mapping_with_source_buffer_owned(
            path,
            self.database.as_ref().unwrap(),
            self.simple_header,
            self.search.file_buffer_size,
            &self.errors,
        )
    }
    fn create_subdatabase(&mut self, kind: SubDatabase, ids: &[OId]) -> Result<(), String> {
        let records = self.database.as_mut().unwrap().selected_records(ids)?;
        if self.search.tuning.no_unlink {
            // The optimized native adapter keeps subdatabases in memory. When
            // Source explicitly requests retained temporaries, materialize the
            // same FASTA sequence payload rather than silently losing that effect.
            let directory = crate::util::io::temp_file::source_temp_directory_path(
                &self.search.tmpdir,
                self.search.output.as_deref().unwrap_or(Path::new("")),
            );
            let data = self.errors.capture(
                crate::util::io::TempFileData::init_in(false, directory, true)
                    .map_err(io::Error::other),
            )?;
            let mut output = self.errors.capture(
                crate::util::io::OutputFile::with_source_profile_path(
                    &data.name,
                    crate::util::io::Compressor::None,
                    "wb",
                    self.search.file_buffer_size,
                )
                .map_err(io::Error::other),
            )?;
            let traits = protein_traits();
            let mut bytes = TextBuffer::new();
            for record in &records {
                bytes.clear();
                crate::util::sequence::format_bytes(
                    &record.sequence,
                    &record.id,
                    None,
                    &mut bytes,
                    "fasta",
                    &traits,
                    160,
                )?;
                output
                    .write_raw(bytes.data())
                    .map_err(|error| error.to_string())?;
            }
            output.close().map_err(|error| error.to_string())?;
        }
        match kind {
            SubDatabase::Members => self.members = records,
            SubDatabase::Centroids => self.centroids = records,
        }
        Ok(())
    }
    fn reset_statistics(&mut self) {}
    fn run_search(
        &mut self,
        config: &ReassignSearchConfig,
        consume: &mut dyn FnMut(&[u8]) -> Result<(), String>,
    ) -> Result<(), String> {
        run_native_search(
            &self.search,
            self.extension,
            self.tantan_min_mask_prob,
            config,
            None,
            std::mem::take(&mut self.centroids),
            std::mem::take(&mut self.members),
            self.database.as_ref().unwrap().letters,
            consume,
            &self.errors,
        )
    }

    fn init_random_access(&mut self) -> Result<(), String> {
        Ok(())
    }
    fn report_reassigned(&mut self, reassigned: usize, members: usize) -> Result<(), String> {
        crate::util::log_stream::message_stream()
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .write(format!("Reassigned members: {reassigned}/{members}\n"))
            .map(|_| ())
            .map_err(|error| error.to_string())
    }
    fn write_clustering(&mut self, clustering: &[OId]) -> Result<(), String> {
        let database = self.database.as_mut().unwrap();
        if !self.oid_output {
            self.errors
                .capture(database.prepare_output_accessions(self.search.threads.max(1) as usize))?;
        }
        let output = self.destination.as_mut().unwrap();
        helpers::output_mem_mapping_stream(output, database, clustering, self.oid_output)?;
        output.flush().map_err(|error| error.to_string())
    }

    fn close_database(&mut self) -> Result<(), String> {
        self.database = None;
        if let Some(output) = self.destination.as_mut() {
            output.close()?;
        }
        Ok(())
    }
}

pub(crate) fn parse_masking(value: &str) -> Result<MaskingMode, String> {
    parse_masking_bytes(value.as_bytes()).map_err(|error| error.to_string())
}
pub(crate) fn parse_masking_bytes(value: &[u8]) -> io::Result<MaskingMode> {
    match value {
        b"0" | b"none" => Ok(MaskingMode::None),
        b"1" | b"tantan" => Ok(MaskingMode::Tantan),
        b"seg" => Ok(MaskingMode::BlastSeg),
        _ => Err(enum_error_bytes(value, b"none, seg, tantan")),
    }
}
pub(crate) fn parse_soft_masking(value: &str) -> Result<MaskingAlgo, String> {
    parse_soft_masking_bytes(value.as_bytes()).map_err(|error| error.to_string())
}
pub(crate) fn parse_soft_masking_bytes(value: &[u8]) -> io::Result<MaskingAlgo> {
    match value {
        b"0" | b"none" => Ok(MaskingAlgo::None),
        b"seg" => Ok(MaskingAlgo::Seg),
        b"tantan" => Ok(MaskingAlgo::Tantan),
        _ => Err(enum_error_bytes(value, b"0, none, seg, tantan")),
    }
}
fn enum_error_bytes(value: &[u8], permitted: &[u8]) -> io::Error {
    let mut message = b"Invalid value for string field: ".to_vec();
    message.extend_from_slice(value);
    message.extend_from_slice(b". Permitted values: ");
    message.extend_from_slice(permitted);
    crate::util::io::RawDiagnosticError(message).into_io()
}

/// Shared source `Search::run(centroid_db, member_db, Mapback)` adapter.
pub(crate) fn run_native_search(
    base: &BlastpConfig,
    extension: Option<crate::search::sensitivity::ExtensionMode>,
    tantan_min_mask_prob: f64,
    config: &ReassignSearchConfig,
    iterate: Option<&[Vec<u8>]>,
    centroids: Vec<FastaRecord>,
    members: Vec<FastaRecord>,
    db_letters: u64,
    consume: &mut dyn FnMut(&[u8]) -> Result<(), String>,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<(), String> {
    let mut search = base.clone();
    search.tuning.self_search = false;
    search.sensitivity = config.sensitivity;
    search.tuning.db_size = db_letters;
    search.block_size = Some(search.block_size.unwrap_or(0.0));
    errors.capture(blastp::begin_algorithm_scope(&mut search))?;
    search.query_cover = config.query_cover;
    search.subject_cover = config.subject_cover;
    search.masking = errors.capture(parse_masking_bytes(&config.masking))?;
    let soft_masking = errors.capture(parse_soft_masking_bytes(&config.soft_masking))?;
    crate::data::sequence_file::check_protein_input(
        centroids.iter().map(|record| record.sequence.as_slice()),
        search.tuning.ignore_warnings,
    )?;
    crate::data::sequence_file::check_protein_input(
        members.iter().map(|record| record.sequence.as_slice()),
        search.tuning.ignore_warnings,
    )?;
    let matrix = errors.capture(blastp::native_score_matrix(&search, 0))?;
    search.scoring_matrix = Some(std::sync::Arc::new(matrix.clone()));
    search.tuning.tantan_min_mask_prob = tantan_min_mask_prob;
    let masker = blastp::native_tantan_masker(&search, &matrix);
    // The concrete edge driver owns all query/reference blocks and rounds.
    // Preserve Source's aligned-before-join skip state rather than rebuilding
    // the pending query set from the final formatted edge records.
    search.iterate = iterate
        .map(<[Vec<u8>]>::to_vec)
        .or_else(|| search.iterate.clone());
    let mut edges = errors.capture(blastp::run_edges_in_memory_bytes(
        &search,
        centroids,
        Some(members),
        config.approx_min_id,
        false,
        false,
        config.self_search,
        search.query_or_subject_cover,
        false,
        crate::chaining::HammingExtConfig {
            approx_min_id: config.approx_min_id,
            query_cover: config.query_cover,
            subject_cover: config.subject_cover,
            max_evalue: search.max_evalue,
            diag_filter_id: config.diag_filter_id,
            diag_filter_cov: config.diag_filter_cov,
            ..Default::default()
        },
        Some(masker.as_ref()),
        soft_masking,
        false,
        Some(db_letters),
        extension,
    ))?;
    edges.sort_by_key(|edge| edge.edge.query);
    let mut begin = 0;
    while begin < edges.len() {
        let end = begin
            + edges[begin..].partition_point(|edge| edge.edge.query == edges[begin].edge.query);
        let mut bytes = Vec::with_capacity((end - begin) * crate::output::edge::EdgeData::SIZE);
        for edge in &edges[begin..end] {
            edge.edge.write(&mut bytes).map_err(|e| e.to_string())?;
        }
        consume(&bytes)?;
        begin = end;
    }
    Ok(())
}

fn protein_traits() -> crate::basic::value::ValueTraits {
    crate::basic::value::ValueTraits::new(
        crate::basic::value::AMINO_ACID_ALPHABET,
        crate::basic::value::MASK_LETTER,
        b" \t\r\n",
        crate::basic::value::SequenceType::AminoAcid,
    )
}

/// Compatibility hook retained for callers; these workflows are fully native.
pub fn requires_legacy(args: &[String]) -> bool {
    requires_legacy_for(args, "reassign")
}

pub(crate) fn requires_legacy_for(_args: &[String], _workflow: &str) -> bool {
    // All consumed controls reach the native search/cascade adapters. Source
    // linclust-chunk-size is read only by external linclust, which recursive
    // reclustering and reassignment never enter (recluster.cpp:128).
    false
}

pub(crate) struct ParsedOptions {
    pub config: ReassignConfig,
    pub search: BlastpConfig,
    pub output: PathBuf,
    pub simple_header: bool,
    pub oid_output: bool,
    pub quiet: bool,
    pub extension: Option<crate::search::sensitivity::ExtensionMode>,
    pub tantan_min_mask_prob: f64,
    pub cascaded: crate::cluster::cascaded::CascadedConfig,
}

/// The full Source registry owns permissions, aliases, and parameter counts.
/// A workflow parser stores only its own fields; common search controls are
/// read once by the shared tuning factory from the original argument list.
pub(crate) fn workflow_arguments(
    args: &[String],
    owned: &[(&str, char)],
) -> Result<Vec<String>, String> {
    let args: Vec<OsString> = args.iter().map(Into::into).collect();
    workflow_arguments_os(&args, owned)
        .map_err(|error| error.to_string())?
        .into_iter()
        .map(|arg| {
            arg.into_string()
                .map_err(|_| "OS argument is not UTF-8; use workflow_arguments_os".to_owned())
        })
        .collect()
}
pub(crate) fn workflow_arguments_os(
    args: &[OsString],
    owned: &[(&str, char)],
) -> io::Result<Vec<OsString>> {
    super::argument_registry::validate_os(args)?;
    let mut selected = args.iter().take(2).cloned().collect::<Vec<_>>();
    let mut at = 2;
    while at < args.len() {
        let token =
            crate::util::io::source_os_bytes(args[at].as_os_str()).map_err(io::Error::other)?;
        let keep = if token.starts_with(b"--") {
            owned
                .iter()
                .any(|&(long, _)| long.as_bytes() == &token[2..])
        } else {
            owned
                .iter()
                .any(|&(_, short)| short != '\0' && token.get(1).copied() == Some(short as u8))
        };
        let begin = at;
        at += 1;
        while at < args.len() {
            let bytes =
                crate::util::io::source_os_bytes(args[at].as_os_str()).map_err(io::Error::other)?;
            if bytes.len() > 1 && bytes[0] == b'-' && !bytes[1].is_ascii_digit() {
                break;
            }
            at += 1;
        }
        if keep {
            selected.extend_from_slice(&args[begin..at]);
        }
    }
    Ok(selected)
}
pub(crate) fn invocation_bytes_os(args: &[OsString]) -> io::Result<Vec<u8>> {
    let mut bytes = Vec::new();
    for (i, arg) in args.iter().enumerate() {
        if i != 0 {
            bytes.push(b' ');
        }
        bytes.extend_from_slice(
            &crate::util::io::source_os_bytes(arg.as_os_str()).map_err(io::Error::other)?,
        );
    }
    Ok(bytes)
}

pub(crate) fn parse_options(args: &[String], workflow: &str) -> Result<ParsedOptions, String> {
    parse_options_with_startup(args, workflow, None)
}

pub(crate) fn parse_options_with_startup(
    args: &[String],
    workflow: &str,
    startup: Option<&super::source_config::Startup>,
) -> Result<ParsedOptions, String> {
    let args: Vec<OsString> = args.iter().map(Into::into).collect();
    parse_options_with_startup_os(&args, workflow, startup).map_err(|error| error.to_string())
}
pub(crate) fn parse_options_with_startup_os(
    args: &[OsString],
    workflow: &str,
    startup: Option<&super::source_config::Startup>,
) -> io::Result<ParsedOptions> {
    let errors = crate::util::io::IoErrorOwner::default();
    let execute = || -> Result<ParsedOptions, String> {
        let mut parser = CommandLineParser::new();
        parser.add_command(workflow, "Refine clustered sequences", 1);
        let group = parser.add_group("Reassign options", vec![1], false);
        let mut allowed = Vec::new();
        macro_rules! option {
            ($name:ident, $id:literal, $short:expr, $default:expr, $min:expr) => {
                let $name = Rc::new(RefCell::new($default));
                allowed.push(($id, $short));
                parser.add_option(group, $id, $short, "", $name.clone(), $default, $min);
            };
        }
        option!(database, "db", 'd', PathBuf::new(), 1);
        option!(clustering, "clusters", '\0', PathBuf::new(), 1);
        option!(output, "out", 'o', PathBuf::new(), 1);
        option!(verbose, "verbose", 'v', false, 0);
        option!(quiet, "quiet", '\0', false, 0);
        option!(log, "log", '\0', false, 0);
        option!(memory, "memory-limit", 'M', Vec::<u8>::new(), 1);
        option!(
            header,
            "header",
            '\0',
            OptionValue::<Vec<Vec<u8>>>::new(),
            0
        );
        option!(
            member_cover,
            "member-cover",
            '\0',
            OptionValue::<f64>::new(),
            1
        );
        option!(
            mutual_cover,
            "mutual-cover",
            '\0',
            OptionValue::<f64>::new(),
            1
        );
        option!(approx_id, "approx-id", '\0', OptionValue::<f64>::new(), 1);
        option!(masking, "masking", '\0', OptionValue::<Vec<u8>>::new(), 1);
        option!(soft_masking, "soft-masking", '\0', Vec::<u8>::new(), 1);
        option!(extension, "ext", '\0', Vec::<u8>::new(), 1);
        option!(oid_output, "oid-output", '\0', false, 0);
        option!(tantan_prob, "tantan-minMaskProb", '\0', 0.9f64, 1);
        option!(no_auto_append, "no-auto-append", '\0', false, 0);
        option!(
            linclust_chunk_size,
            "linclust-chunk-size",
            '\0',
            b"10G".to_vec(),
            1
        );
        let source_general = parser.add_group(
            "Advanced/general",
            if workflow == "reassign" {
                vec![1]
            } else {
                vec![2]
            },
            false,
        );
        macro_rules! source_option {
            ($name:ident, $id:literal, $default:expr, $min:expr) => {
                let $name = Rc::new(RefCell::new($default));
                allowed.push(($id, '\0'));
                parser.add_option(source_general, $id, '\0', "", $name.clone(), $default, $min);
            };
        }
        source_option!(no_parse_seqids, "no-parse-seqids", false, 0);

        option!(
            depth,
            "connected-component-depth",
            '\0',
            Vec::<Vec<u8>>::new(),
            1
        );
        option!(no_reassign, "no-reassign", '\0', false, 0);
        let cluster_group = parser.add_group(
            "Clustering options",
            if workflow == "recluster" {
                vec![1]
            } else {
                vec![2]
            },
            false,
        );
        macro_rules! cluster_option {
            ($name:ident, $id:literal, $default:expr, $min:expr) => {
                let $name = Rc::new(RefCell::new($default));
                allowed.push(($id, '\0'));
                parser.add_option(cluster_group, $id, '\0', "", $name.clone(), $default, $min);
            };
        }
        cluster_option!(cluster_steps, "cluster-steps", Vec::<Vec<u8>>::new(), 1);
        cluster_option!(round_coverage, "round-coverage", Vec::<Vec<u8>>::new(), 1);
        cluster_option!(round_approx, "round-approx-id", Vec::<Vec<u8>>::new(), 1);
        cluster_option!(kmer_ranking, "kmer-ranking", false, 0);
        cluster_option!(alignment_output, "aln-out", PathBuf::new(), 1);
        let workflow_args = errors.capture(workflow_arguments_os(args, &allowed))?;
        let argv = workflow_args
            .iter()
            .map(OsString::as_os_str)
            .collect::<Vec<_>>();
        errors.capture(parser.store_os(&argv, &mut 0))?;
        if startup.is_none() && args.len() == 2 {
            print!("{}", parser.print_documentation(1));
        }
        parser.require("db")?;
        parser.require("clusters")?;
        let header = header.borrow().get_present().ok();
        let header_args = header
            .as_ref()
            .map(|h| h.iter().map(Vec::as_slice).collect::<Vec<_>>());
        let simple_header = errors.capture(TabularFormat::header_format_bytes(
            Workflow::Cluster,
            header_args.as_deref(),
        ))? == Header::Simple;
        let config = ReassignConfig {
            database: database.borrow().clone(),
            clustering: clustering.borrow().clone(),
            member_cover: member_cover.borrow().get_present().ok(),
            mutual_cover: mutual_cover.borrow().get_present().ok(),
            approx_min_id: approx_id.borrow().get_present().ok(),
            cluster_steps: cluster_steps.borrow().clone(),
            masking: masking.borrow().get_present().ok(),
            soft_masking: (!soft_masking.borrow().is_empty())
                .then(|| soft_masking.borrow().clone()),
            ..Default::default()
        };
        let memory = if memory.borrow().is_empty() {
            None
        } else {
            Some(errors.capture(crate::util::string::interpret_number_bytes(
                &memory.borrow(),
            ))? as u64)
        };
        let extension_name = extension.borrow().clone();
        let extension = if extension.borrow().is_empty() {
            None
        } else {
            Some(
                errors.capture(crate::search::sensitivity::ExtensionMode::from_bytes(
                    &extension.borrow(),
                ))?,
            )
        };
        // Source switches the existing common configuration to edge output; it
        // does not construct a second search configuration with fresh defaults.
        // Project only its forced output format before the shared parser, so an
        // unrelated requested DAA/XML format cannot affect this internal stream.
        let mut search_args = args.to_vec();
        search_args.extend([OsString::from("--outfmt"), OsString::from("edge")]);
        let invocation = errors.capture(invocation_bytes_os(args))?;
        let mut search = errors.capture(super::dispatch::parse_native_blastp_config_os(
            &search_args,
            &invocation,
        ))?;
        if let Some(startup) = startup {
            startup.apply_search(&mut search);
        }
        search.output = None;
        search.max_target_seqs = 1;
        search.tuning.setup.aln_out = alignment_output.borrow().clone();
        // These options affect other workflows, not reassignment's edge stream.
        let _ = (no_auto_append, no_parse_seqids, linclust_chunk_size);
        let cascaded = crate::cluster::cascaded::CascadedConfig {
            matrix: search.matrix.clone(),
            gap_open: search.gap_open,
            gap_extend: search.gap_extend,
            max_evalue: search.max_evalue,
            comp_based_stats: search.comp_based_stats as i32,
            motif_masking: search.motif_masking.clone(),
            extension_mode: extension_name,
            memory_limit: memory.unwrap_or(16_000_000_000),
            threads: if search.threads > 0 {
                search.threads
            } else {
                crate::util::system::hardware_concurrency() as i32
            },
            helpers: crate::cluster::cascaded::helpers::CascadedHelpersConfig {
                cluster_steps: cluster_steps.borrow().clone(),
                connected_component_depth: depth.borrow().clone(),
            },
            round_coverage: round_coverage.borrow().clone(),
            round_approx_id: round_approx.borrow().clone(),
            kmer_ranking: *kmer_ranking.borrow(),
            no_gvc_reassign: *no_reassign.borrow(),
            alignment_output: (!alignment_output.borrow().as_os_str().is_empty())
                .then(|| alignment_output.borrow().clone()),
            ..Default::default()
        };
        let parsed = ParsedOptions {
            config,
            search,
            output: startup.map_or_else(
                || output.borrow().clone(),
                |startup| startup.config.output_file.clone(),
            ),
            simple_header,
            oid_output: *oid_output.borrow(),
            extension,
            tantan_min_mask_prob: *tantan_prob.borrow(),
            quiet: *quiet.borrow(),
            cascaded,
        };
        let _ = (verbose, log);
        Ok(parsed)
    };
    errors.finish(execute())
}

pub fn run(args: &[String]) -> io::Result<()> {
    let args: Vec<OsString> = args.iter().map(Into::into).collect();
    run_os(&args)
}
pub fn run_os(args: &[OsString]) -> io::Result<()> {
    run_with_startup_os(args, None)
}

pub(crate) fn run_with_startup(
    args: &[String],
    startup: Option<std::sync::Arc<super::source_config::Startup>>,
) -> io::Result<()> {
    let args: Vec<OsString> = args.iter().map(Into::into).collect();
    run_with_startup_os(&args, startup)
}
pub(crate) fn run_with_startup_os(
    args: &[OsString],
    startup: Option<std::sync::Arc<super::source_config::Startup>>,
) -> io::Result<()> {
    let startup = match startup {
        Some(startup) => startup,
        None => super::source_config::prepare_os(args)?,
    };
    let mut options = parse_options_with_startup_os(args, "reassign", Some(&startup))?;
    let errors = crate::util::io::IoErrorOwner::default();
    let mut backend = Backend {
        errors: errors.clone(),
        database: None,
        members: Vec::new(),
        centroids: Vec::new(),
        destination: None,
        output: options.output,
        simple_header: options.simple_header,
        oid_output: options.oid_output,
        extension: options.extension,
        tantan_min_mask_prob: options.tantan_min_mask_prob,
        search: options.search,
    };
    errors.finish(reassign::reassign(&mut options.config, &mut backend))?;
    Ok(())
}
