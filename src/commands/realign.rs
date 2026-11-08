//! Concrete adapters for the translated cluster realignment workflow.

use std::cell::RefCell;
use std::ffi::OsString;
use std::io;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use super::cluster_database::NativeDatabase;
use crate::align::hsp::ByteHspContext as HspContext;
use crate::basic::value::OId;
use crate::cluster::helpers;
use crate::cluster::output::{self, ClusterOutputConfig};
use crate::cluster::realign::{self, DatabaseMetadata, RealignBackend, RealignConfig};
use crate::data::sequence_file::SequenceFileFlags;
use crate::dp::swipe::HspValues;
use crate::output::format::{write_tabular_context_row_bytes, TabularFormat, Workflow};
use crate::stats::cbs::CbsMode;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::command_line_parser::CommandLineParser;
use crate::util::data_structures::FlatArray;
use crate::util::options::OptionValue;

struct Backend {
    errors: crate::util::io::IoErrorOwner,
    database: Option<NativeDatabase>,
    file_buffer_size: usize,
    matrix: ScoreMatrix,
    output_config: ClusterOutputConfig,
    io_config: output::ClusterOutputIo,
    simple_header: bool,
}

impl RealignBackend<Vec<u8>> for Backend {
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
    fn read_centroid_sorted(&mut self, path: &str) -> Result<(FlatArray<OId>, Vec<OId>), String> {
        self.read_centroid_sorted_path(Path::new(path))
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
                Some(self.file_buffer_size),
                &self.io_config.tmpdir,
                self.io_config.no_unlink,
                &self.errors,
            ))?;
        let metadata = DatabaseMetadata {
            sequence_count: database.count,
            letters: database.letters,
            titles_lazy: false,
        };
        self.database = Some(database);
        Ok(metadata)
    }
    fn set_score_matrix_db_letters(&mut self, letters: u64) -> Result<(), String> {
        self.matrix
            .set_db_letters(self.output_config.db_size.unwrap_or(letters));
        Ok(())
    }
    fn init_random_access(&mut self) -> Result<(), String> {
        Ok(())
    }
    fn read_centroid_sorted_path(
        &mut self,
        path: &Path,
    ) -> Result<(FlatArray<OId>, Vec<OId>), String> {
        helpers::read_centroid_sorted_with_source_buffer_owned(
            path,
            self.database.as_ref().unwrap(),
            self.simple_header,
            self.file_buffer_size,
            &self.errors,
        )
    }
    fn realign_clusters(
        &mut self,
        clusters: &FlatArray<OId>,
        centroids: &[OId],
        hsp_values: HspValues,
        emit: &mut dyn FnMut(&HspContext) -> Result<(), String>,
    ) -> Result<(), String> {
        let statistics = output::realign_clusters_with_io_errors_bytes(
            clusters,
            centroids,
            self.database.as_mut().unwrap(),
            emit,
            hsp_values,
            &self.output_config,
            &mut self.matrix,
            &self.io_config,
            &self.errors,
        )?;
        crate::basic::statistics::merge_statistics(&statistics);
        // cluster/output.cpp reports its accumulated statistics without a
        // Search::run reset or TotalTime update.
        let _ = crate::basic::statistics::report_statistics();
        Ok(())
    }
    fn close_database(&mut self) -> Result<(), String> {
        self.database = None;
        Ok(())
    }
}

/// Standard upstream command options, parsed by its translated parser.
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
    let errors = crate::util::io::IoErrorOwner::default();
    let execute = || -> Result<(), String> {
        let mut parser = CommandLineParser::new();
        parser.add_command(
            "realign",
            "Realign clustered sequences against their centroids",
            1,
        );
        let group = parser.add_group("Realign options", vec![1], false);
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
        option!(threads, "threads", 'p', 0i32, 1);
        option!(verbose, "verbose", 'v', false, 0);
        option!(quiet, "quiet", '\0', false, 0);
        option!(log, "log", '\0', false, 0);
        option!(tmpdir, "tmpdir", 't', PathBuf::new(), 1);
        option!(memory, "memory-limit", 'M', Vec::<u8>::new(), 1);
        option!(file_buffer, "file-buffer-size", '\0', 67_108_864i64, 1);
        option!(
            header,
            "header",
            '\0',
            OptionValue::<Vec<Vec<u8>>>::new(),
            0
        );
        option!(format, "outfmt", 'f', Vec::<Vec<u8>>::new(), 1);
        option!(qoffset, "qnum-offset", '\0', 0i64, 1);
        option!(soffset, "snum-offset", '\0', 0i64, 1);
        option!(cbs, "comp-based-stats", '\0', 1u32, 1);
        option!(masking, "masking", '\0', Vec::<u8>::new(), 1);
        option!(soft_masking, "soft-masking", '\0', Vec::<u8>::new(), 1);
        option!(no_block_limit, "no-block-size-limit", '\0', false, 0);
        let workflow_args =
            errors.capture(super::reassign::workflow_arguments_os(args, &allowed))?;
        let argv = workflow_args
            .iter()
            .map(OsString::as_os_str)
            .collect::<Vec<_>>();
        errors.capture(parser.store_os(&argv, &mut 0))?;
        CbsMode::from_search_code(*cbs.borrow(), false, 0)?;
        parser.require("db")?;
        parser.require("clusters")?;
        let header = if header.borrow().present() {
            Some(header.borrow().get_present()?)
        } else {
            None
        };
        let mut config = RealignConfig {
            database: database.borrow().clone(),
            clustering: clustering.borrow().clone(),
            output_format: format.borrow().clone(),
            output_header: header,
            ..RealignConfig::default()
        };
        let tabular = errors.capture(realign::tabular_format_io(&mut config))?;
        let header_args = config
            .output_header
            .as_ref()
            .map(|h| h.iter().map(Vec::as_slice).collect::<Vec<_>>());
        let simple_header = errors.capture(TabularFormat::header_format_bytes(
            Workflow::Cluster,
            header_args.as_deref(),
        ))? == crate::output::format::Header::Simple;
        let memory = if memory.borrow().is_empty() {
            b"16G".to_vec()
        } else {
            memory.borrow().clone()
        };
        let block_size = errors.capture(crate::util::string::interpret_number_bytes(&memory))? / 2;
        let mut tuning = errors.capture(super::dispatch::parse_native_search_tuning_os(args))?;
        startup.apply_tuning(&mut tuning);
        let values = errors.capture(super::argument_registry::values_os(args))?;
        let matrix_name = if let Some(value) = values.get("matrix").and_then(|values| values.last())
        {
            errors.capture(
                crate::util::io::source_os_bytes(value)
                    .map(|value| value.into_owned())
                    .map_err(io::Error::other),
            )?
        } else {
            b"BLOSUM62".to_vec()
        };
        let gap_open = errors.capture(super::dispatch::parse_signed_search_option_os(
            args, "gapopen", -1,
        ))?;
        let gap_extend = errors.capture(super::dispatch::parse_signed_search_option_os(
            args,
            "gapextend",
            -1,
        ))?;
        let file_buffer_size = usize::try_from(*file_buffer.borrow()).map_err(|_| {
            "Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage.".to_owned()
        })?;
        super::blastp::preflight_file_buffer_size(file_buffer_size)
            .map_err(|error| error.to_string())?;
        let mut backend =
            Backend {
                errors: errors.clone(),
                database: None,
                file_buffer_size,
                matrix: errors.capture(super::blastp::native_score_matrix_with_tuning_bytes(
                    &matrix_name,
                    gap_open,
                    gap_extend,
                    0,
                    &tuning,
                    0,
                ))?,
                io_config: output::ClusterOutputIo {
                    tmpdir: startup.config.tmpdir.clone().into(),
                    file_buffer_size,
                    no_unlink: tuning.no_unlink,
                },
                output_config: ClusterOutputConfig {
                    cbs_window: tuning.cbs_window,
                    block_size: block_size as u64,
                    comp_based_stats_hauser: *cbs.borrow() == 1,
                    db_size: (tuning.db_size != 0).then_some(tuning.db_size),
                    swipe_task_size: tuning.swipe_task_size,
                    transcript_len_estimate: tuning.transcript_len_estimate,
                    threads_align: tuning.threads_align,
                    target_parallel_verbosity: tuning.target_parallel_verbosity,
                    cutoff_score_8bit: tuning.cutoff_score_8bit,
                    max_swipe_dp: tuning.max_swipe_dp,
                    band_bin: tuning.band_bin,
                    col_bin: tuning.col_bin,
                    approx_backtrace: tuning.approx_backtrace,
                    // Source realign forces max_evalue=DBL_MAX, while the
                    // shared SWIPE reverse pass still consumes these globals.
                    query_cover: errors.capture(super::dispatch::parse_float_search_option_os(
                        args,
                        "query-cover",
                        0.0,
                    ))?,
                    subject_cover: errors.capture(
                        super::dispatch::parse_float_search_option_os(args, "subject-cover", 0.0),
                    )?,
                    query_or_target_cover: errors.capture(
                        super::dispatch::parse_float_search_option_os(
                            args,
                            "query-or-subject-cover",
                            0.0,
                        ),
                    )?,
                    approx_min_id: errors.capture(
                        super::dispatch::parse_float_search_option_os(args, "approx-id", 0.0),
                    )?,
                    cbs_matrix_scale: tuning.matrix_scale,
                    ..ClusterOutputConfig::default()
                },
                simple_header,
            };
        let mut destination =
            errors.capture(super::reassign::SourceWorkflowOutput::new_io_path(
                &startup.config.output_file,
                file_buffer_size,
            ))?;
        let qoffset = *qoffset.borrow() as u64;
        let soffset = *soffset.borrow() as u64;
        // Source resolves omitted/-p0 to hardware_concurrency before this
        // workflow; do not let RAYON_NUM_THREADS replace that resolved budget.
        let threads = startup.config.threads;
        let mut pool = rayon::ThreadPoolBuilder::new();
        if threads > 0 {
            pool = pool.num_threads(threads as usize);
        }
        pool.build().map_err(|e| e.to_string())?.install(|| {
            realign::realign_bytes(
                &mut config,
                &mut backend,
                &mut destination,
                |context, buffer| {
                    write_tabular_context_row_bytes(
                        buffer,
                        context,
                        &tabular.fields,
                        false,
                        false,
                        qoffset,
                        soffset,
                    )
                    .map_err(|e| e.to_string())
                },
            )
        })?;
        destination.close()?;
        // Upstream accepts these general/search settings but realign never
        // applies its search masking pipeline to the supplied cluster pairs.
        let _ = (
            output,
            verbose,
            quiet,
            log,
            tmpdir,
            no_block_limit,
            masking,
            soft_masking,
        );
        Ok(())
    };
    errors.finish(execute())
}
