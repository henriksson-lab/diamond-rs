//! Concrete adapters for the recursive source-mirrored reclustering workflow.
use std::ffi::OsString;
use std::io::{self, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use super::cluster_database::NativeDatabase;
use super::reassign::{self, ParsedOptions};
use crate::align::hsp::ByteHspContext as HspContext;
use crate::basic::value::{OId, SuperBlockId};
use crate::cluster::cascaded::recluster::{
    self, ReclusterBackend, ReclusterConfig, ReclusterSearchConfig,
};
use crate::cluster::helpers;
use crate::cluster::output::{self, ClusterOutputConfig};
use crate::cluster::realign::DatabaseMetadata;
use crate::cluster::reassign::ReassignSearchConfig;
use crate::data::sequence_file::SequenceFileFlags;
use crate::dp::swipe::HspValues;
use crate::stats::cbs::CbsMode;
use crate::util::data_structures::FlatArray;

type Database = Arc<Mutex<NativeDatabase>>;

struct Backend {
    errors: crate::util::io::IoErrorOwner,
    options: ParsedOptions,
    destination: Option<reassign::SourceWorkflowOutput>,
}

impl Backend {
    fn sync_thresholds(&mut self, config: &ReclusterConfig) {
        let core = &mut self.options.cascaded;
        core.member_cover = config.member_cover.unwrap_or(80.0);
        core.mutual_cover = config.mutual_cover;
        core.approx_min_id = config.approx_min_id.unwrap();
        core.diag_filter_id = config.diag_filter_id;
        core.diag_filter_cov = config.diag_filter_cov;
        core.db_size = config.db_size.unwrap_or(0);
    }
}

impl ReclusterBackend<Vec<u8>> for Backend {
    fn io_errors(&self) -> crate::util::io::IoErrorOwner {
        self.errors.clone()
    }

    fn open_database(
        &mut self,
        path: &str,
        flags: SequenceFileFlags,
    ) -> Result<(Self::Database, DatabaseMetadata), String> {
        self.open_database_path(Path::new(path), flags)
    }
    fn read_clustering(
        &mut self,
        path: &str,
        database: &Self::Database,
    ) -> Result<Vec<OId>, String> {
        self.read_clustering_path(Path::new(path), database)
    }

    type Database = Database;

    fn open_database_path(
        &mut self,
        path: &Path,
        flags: SequenceFileFlags,
    ) -> Result<(Database, DatabaseMetadata), String> {
        let database = self
            .errors
            .capture(NativeDatabase::open_source_with_buffer_owned_path(
                path,
                flags,
                Some(self.options.search.file_buffer_size),
                &self.options.search.tmpdir,
                self.options.search.tuning.no_unlink,
                &self.errors,
            ))?;
        let metadata = DatabaseMetadata {
            sequence_count: database.count,
            letters: database.letters,
            titles_lazy: false,
        };
        self.options.search.tuning.db_size = database.letters;
        Ok((Arc::new(Mutex::new(database)), metadata))
    }
    fn close_database(&mut self, database: Database) -> Result<(), String> {
        drop(database);
        if let Some(output) = self.destination.as_mut() {
            output.close()?;
        }
        Ok(())
    }
    fn open_output(&mut self) -> Result<(), String> {
        let mut output = self
            .errors
            .capture(reassign::SourceWorkflowOutput::new_io_path(
                &self.options.output,
                self.options.search.file_buffer_size,
            ))?;
        if self.options.simple_header {
            writeln!(output, "centroid\tmember").map_err(|e| e.to_string())?;
        }
        self.destination = Some(output);
        Ok(())
    }
    fn read_clustering_path(
        &mut self,
        path: &Path,
        database: &Database,
    ) -> Result<Vec<OId>, String> {
        helpers::read_mapping_with_source_buffer_owned(
            path,
            &*database.lock().unwrap(),
            self.options.simple_header,
            self.options.search.file_buffer_size,
            &self.errors,
        )
    }

    fn sequence_count(&self, database: &Database) -> usize {
        database.lock().unwrap().count as usize
    }
    fn create_subdatabase(&mut self, database: &Database, ids: &[OId]) -> Result<Database, String> {
        let records = database.lock().unwrap().selected_records(ids)?;
        Ok(Arc::new(Mutex::new(NativeDatabase::from_records(records)?)))
    }
    fn realign(
        &mut self,
        database: &Database,
        clusters: &FlatArray<OId>,
        centroids: &[OId],
        hsp_values: HspValues,
        config: &ReclusterConfig,
        callback: &mut dyn FnMut(&HspContext),
    ) -> Result<(), String> {
        self.sync_thresholds(config);
        let core = &self.options.cascaded;
        let mut search = self.options.search.clone();
        search.tuning.self_search = false;
        search.matrix = core.matrix.clone();
        search.gap_open = core.gap_open;
        search.gap_extend = core.gap_extend;
        search.tuning.db_size = core.db_size;
        let mut matrix =
            super::blastp::native_score_matrix(&search, database.lock().unwrap().letters)
                .map_err(|error| error.to_string())?;
        let output = ClusterOutputConfig {
            cbs_window: search.tuning.cbs_window,
            block_size: core.memory_limit / 2,
            db_size: (core.db_size != 0).then_some(core.db_size),
            comp_based_stats_hauser: core.comp_based_stats == 1,
            max_evalue: core.max_evalue,
            query_cover: core.query_cover,
            subject_cover: core.subject_cover,
            query_or_target_cover: core.query_or_target_cover,
            approx_min_id: core.approx_min_id,
            swipe_task_size: search.tuning.swipe_task_size,
            transcript_len_estimate: search.tuning.transcript_len_estimate,
            threads_align: search.tuning.threads_align,
            target_parallel_verbosity: search.tuning.target_parallel_verbosity,
            cutoff_score_8bit: search.tuning.cutoff_score_8bit,
            max_swipe_dp: search.tuning.max_swipe_dp,
            band_bin: search.tuning.band_bin,
            col_bin: search.tuning.col_bin,
            approx_backtrace: search.tuning.approx_backtrace,
            cbs_matrix_scale: search.tuning.matrix_scale,
            ..Default::default()
        };
        let io_config = output::ClusterOutputIo {
            tmpdir: search.tmpdir.clone(),
            file_buffer_size: search.file_buffer_size,
            no_unlink: search.tuning.no_unlink,
        };
        let statistics = output::realign_clusters_with_io_errors_bytes(
            clusters,
            centroids,
            &mut *database.lock().unwrap(),
            |context| {
                callback(context);
                Ok(())
            },
            hsp_values,
            &output,
            &mut matrix,
            &io_config,
            &self.errors,
        )?;
        crate::basic::statistics::merge_statistics(&statistics);
        // cluster/output.cpp reports its accumulated statistics without a
        // Search::run reset or TotalTime update.
        let _ = crate::basic::statistics::report_statistics();
        Ok(())
    }
    fn reset_statistics(&mut self) {}
    fn run_search(
        &mut self,
        centroids: &Database,
        queries: &Database,
        config: &ReclusterSearchConfig,
        consume: &mut dyn FnMut(&[u8]) -> Result<(), String>,
    ) -> Result<(), String> {
        let records = |database: &Database| {
            let mut database = database.lock().unwrap();
            let ids = (0..database.count).collect::<Vec<_>>();
            database.selected_records(&ids)
        };
        let core = &mut self.options.cascaded;
        core.query_cover = config.query_cover;
        core.subject_cover = config.subject_cover;
        core.query_or_target_cover = config.query_or_target_cover;
        core.sensitivity = config.sensitivity;
        core.lin_stage1_query = false;
        core.chunk_size = config.chunk_size;
        core.lowmem = config.lowmem as i32;
        // Source resets these controls immediately before the reassignment
        // search, independently of the earlier cluster realignment phase.
        let mut native_search = self.options.search.clone();
        native_search.tuning.self_search = false;
        native_search.index_chunks = config.lowmem as u32;
        native_search.block_size = Some(config.chunk_size);
        native_search.query_or_subject_cover = config.query_or_target_cover;
        native_search.tuning.setup.lin_stage1_query = config.lin_stage1_query;
        native_search.tuning.setup.lin_stage1_target = config.lin_stage1_target;
        let search = ReassignSearchConfig {
            command: config.command,
            max_target_seqs: config.max_target_seqs,
            output_format: config.output_format,
            self_search: config.self_search,
            query_cover: config.query_cover,
            subject_cover: config.subject_cover,
            sensitivity: config.sensitivity,
            approx_min_id: config.approx_min_id,
            masking: config.masking.clone(),
            soft_masking: config.soft_masking.clone(),
            diag_filter_id: config.diag_filter_id,
            diag_filter_cov: config.diag_filter_cov,
        };
        reassign::run_native_search(
            &native_search,
            self.options.extension,
            self.options.tantan_min_mask_prob,
            &search,
            Some(&config.iterate),
            records(centroids)?,
            records(queries)?,
            core.db_size,
            consume,
            &self.errors,
        )
    }
    fn cascaded(
        &mut self,
        database: &Database,
        linear: bool,
        config: &mut ReclusterConfig,
    ) -> Result<Vec<SuperBlockId>, String> {
        self.sync_thresholds(config);
        let core = &mut self.options.cascaded;
        core.masking = self.errors.capture(reassign::parse_masking_bytes(
            config.masking.as_ref().unwrap(),
        ))?;
        core.soft_masking = self.errors.capture(reassign::parse_soft_masking_bytes(
            config.soft_masking.as_ref().unwrap(),
        ))?;
        let mut database = database.lock().unwrap();
        let ids = (0..database.count).collect::<Vec<_>>();
        let records = database.selected_records(&ids)?;
        // Source cascaded/recluster.cpp calls cascaded(unmapped, false)
        // directly here; its surrounding phase messages belong to the driver.
        let mapping = super::cluster_cmd::cascaded_records_with_search_owned(
            records,
            core,
            linear,
            self.options.tantan_min_mask_prob as f32,
            &self.options.search,
            &self.errors,
        )?;
        config.db_size = Some(core.db_size);
        config.approx_min_id = Some(core.approx_min_id);
        self.options.search.comp_based_stats =
            CbsMode::from_search_code(core.comp_based_stats as u32, false, 0)?;
        self.options.search.max_evalue = core.max_evalue;
        self.options.extension = if core.extension_mode.is_empty() {
            None
        } else {
            Some(
                self.errors
                    .capture(crate::search::sensitivity::ExtensionMode::from_bytes(
                        &core.extension_mode,
                    ))?,
            )
        };
        Ok(mapping)
    }
    fn init_random_access(&mut self, _database: &mut Database) -> Result<(), String> {
        Ok(())
    }
    fn write_clustering(&mut self, database: &Database, clustering: &[OId]) -> Result<(), String> {
        let mut database = database.lock().unwrap();
        if !self.options.oid_output {
            self.errors.capture(
                database.prepare_output_accessions(self.options.search.threads.max(1) as usize),
            )?;
        }
        helpers::output_mem_mapping_stream(
            self.destination.as_mut().unwrap(),
            &*database,
            clustering,
            self.options.oid_output,
        )
    }
}

pub fn requires_legacy(args: &[String]) -> bool {
    reassign::requires_legacy_for(args, "recluster")
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
    let mut options = reassign::parse_options_with_startup_os(args, "recluster", Some(&startup))?;
    let parsed = std::mem::take(&mut options.config);
    let mut config = ReclusterConfig {
        database: parsed.database,
        clustering: parsed.clustering,
        member_cover: parsed.member_cover,
        mutual_cover: parsed.mutual_cover,
        approx_min_id: parsed.approx_min_id,
        soft_masking: parsed.soft_masking,
        masking: parsed.masking,
        diag_filter_id: parsed.diag_filter_id,
        diag_filter_cov: parsed.diag_filter_cov,
        cluster_steps: parsed.cluster_steps,
        db_size: parsed.db_size,
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.cascaded.threads as usize)
        .build()
        .map_err(io::Error::other)?;
    let errors = crate::util::io::IoErrorOwner::default();
    let mut backend = Backend {
        errors: errors.clone(),
        options,
        destination: None,
    };
    errors.finish(pool.install(|| recluster::recluster_bytes(&mut config, &mut backend)))?;
    Ok(())
}
