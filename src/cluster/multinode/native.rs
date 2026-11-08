//! Concrete retained-library backend for Source's historical multinode driver.
//!
//! Source's command wrapper comments out its multinode call. These entry points
//! connect the retained library algorithm to real native search and graph
//! owners without adding a CLI command or rebuilding argv between searches.
use std::io;
use std::path::Path;

use super::data::{
    ByteSequenceRecord, FastaSequenceSource, IndexedSequenceReader, SequenceReader, SequenceRecord,
    SequenceSource,
};
use super::multinode::{
    multinode_io, BlockSearchRequest, MultinodeBackend, MultinodeConfig, VolumedFile,
};
use crate::basic::value::OId;
use crate::commands::blastp::{self, BlastpConfig};
use crate::data::sequence_file::SequenceFileFlags;
use crate::run::config::Algo;
use crate::stats::cbs::CbsMode;
use crate::tools::greedy_vertex_cover::{
    greedy_vertex_cover_from_source_files_with_threads, GreedyVertexCoverConfig,
};
use crate::util::io::IoErrorOwner;

/// Retains the actual initialized search configuration and the untouched graph
/// options that Source obtains from its shared Config. Only fields assigned by
/// multinode.cpp are changed; scoring/masking kernels and raw path owners are
/// inherited from the supplied native configuration.
pub struct NativeMultinodeBackend {
    search: BlastpConfig,
    graph: GreedyVertexCoverConfig<Vec<u8>>,
    source: FastaSequenceSource,
    errors: IoErrorOwner,
}

impl NativeMultinodeBackend {
    pub fn new(search: BlastpConfig, graph: GreedyVertexCoverConfig<Vec<u8>>) -> Self {
        Self::with_io_errors(search, graph, IoErrorOwner::default())
    }

    /// Retain the invocation's existing error policy and primary-error slot.
    /// No input/output stream is opened by constructing this backend.
    pub fn with_io_errors(
        search: BlastpConfig,
        graph: GreedyVertexCoverConfig<Vec<u8>>,
        errors: IoErrorOwner,
    ) -> Self {
        Self {
            search,
            graph,
            source: FastaSequenceSource,
            errors,
        }
    }

    pub fn search_config(&self) -> &BlastpConfig {
        &self.search
    }

    pub fn graph_config(&self) -> &GreedyVertexCoverConfig<Vec<u8>> {
        &self.graph
    }

    /// Project the same assignments made immediately before Source Search::run.
    /// This is a trusted typed boundary, not a second CLI validation phase.
    fn apply_search(&mut self, request: &BlockSearchRequest) -> io::Result<()> {
        if request.command != "blastp" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Multinode block searches require blastp.",
            ));
        }
        self.search.query_files = request.query_file.iter().cloned().collect();
        self.search.database = request.database.clone();
        self.search.output = Some(request.output_file.clone());
        self.search.outfmt = request
            .output_fields
            .iter()
            .map(|field| field.as_bytes().to_vec())
            .collect();
        self.search.tuning.output.qnum_offset = request.query_offset as i64;
        self.search.tuning.output.snum_offset = request.subject_offset as i64;
        self.search.query_cover = request.query_cover;
        self.search.subject_cover = request.subject_cover;
        self.search.query_or_subject_cover = request.query_or_target_cover;
        self.search.tuning.self_search = request.self_search;
        self.search.tuning.setup.lin_stage1_query = request.lin_stage1_query;
        self.search.tuning.setup.lin_stage1_combo = request.lin_stage1_combo;
        self.search.tuning.setup.lin_stage1_target = request.lin_stage1_target;
        self.search.round_linear_target = request.lin_stage1_target;
        self.search.max_target_seqs = request.max_target_seqs;
        self.search.toppercent = request.top_percent;
        if request.iterate_empty {
            self.search.iterate = None;
        }
        if request.double_indexed {
            self.search.query_indexed = false;
            self.search.contiguous_seed_mode = false;
            self.search.tuning.setup.contiguous_seed_mode = false;
            self.search
                .tuning
                .algorithm
                .force_source_algorithm(request.sensitivity, Algo::DoubleIndexed)?;
        }
        self.search.tuning.mapany = request.map_any;
        self.search.index_chunks = request.lowmem as u32;
        self.search.block_size = Some(request.chunk_size as f64);
        self.search.tuning.db_size = request.db_size;
        self.search.comp_based_stats = CbsMode::from_search_code(
            request.comp_based_stats as u32,
            false,
            self.search.global_ranking_targets,
        )
        .map_err(io::Error::other)?;
        self.search.tuning.min_length_ratio = request.min_length_ratio;
        self.search.tuning.search_min_length_ratio = None;
        self.search.sensitivity = request.sensitivity;
        self.search.approx_min_id = request.approx_min_id;
        self.search.tuning.setup.approx_min_id = request.approx_min_id;
        self.search.max_evalue = request.max_evalue;
        // The subsequent GVC call observes the same mutated global coverage.
        self.graph.query_or_target_cover = request.query_or_target_cover;
        Ok(())
    }

    pub fn run_search_io(&mut self, request: &BlockSearchRequest) -> io::Result<()> {
        self.apply_search(request)?;
        let _ = crate::util::system::log_rss();
        blastp::run(&self.search)
    }

    pub fn greedy_vertex_cover_io(
        &mut self,
        database: &Path,
        edges: &Path,
        output: &Path,
    ) -> io::Result<()> {
        // Source mutates the same Config object used by the preceding search.
        self.search.database = database.to_path_buf();
        self.search.output = Some(output.to_path_buf());
        self.graph.database = database.to_path_buf();
        self.graph.edges = edges.to_path_buf();
        self.graph.output_file = Some(output.to_path_buf());
        self.graph.edge_format.clear();
        greedy_vertex_cover_from_source_files_with_threads(
            &self.graph,
            self.search.file_buffer_size,
            self.search.threads as usize,
        )
        .map(|_| ())
    }
}

impl SequenceSource for NativeMultinodeBackend {
    fn index_volume(&self, path: &Path, index: &Path, errors: &IoErrorOwner) -> Result<(), String> {
        self.source.index_volume(path, index, errors)
    }
    fn open_indexed_records(
        &self,
        path: &Path,
        index: &Path,
        flags: SequenceFileFlags,
        errors: &IoErrorOwner,
    ) -> Result<Box<dyn IndexedSequenceReader>, String> {
        self.source.open_indexed_records(path, index, flags, errors)
    }
    fn open_records(
        &self,
        path: &Path,
        errors: &IoErrorOwner,
    ) -> Result<Box<dyn SequenceReader>, String> {
        self.source.open_records(path, errors)
    }
    fn read_records_bytes_with_errors(
        &self,
        path: &Path,
        errors: &IoErrorOwner,
    ) -> Result<Vec<ByteSequenceRecord>, String> {
        self.source.read_records_bytes_with_errors(path, errors)
    }
    fn read_records(&self, path: &Path) -> Result<Vec<SequenceRecord>, String> {
        self.source.read_records(path)
    }
    fn read_records_bytes(&self, path: &Path) -> Result<Vec<ByteSequenceRecord>, String> {
        self.source
            .read_records_bytes_with_errors(path, &self.errors)
    }
}

impl MultinodeBackend for NativeMultinodeBackend {
    fn io_errors(&self) -> IoErrorOwner {
        self.errors.clone()
    }
    fn begin_workflow(&mut self, config: &MultinodeConfig) -> Result<(), String> {
        self.search.threads = i32::try_from(config.threads)
            .map_err(|_| "Multinode thread count exceeds signed32 range.".to_owned())?;
        self.search.tuning.setup.threads = self.search.threads;
        self.search.file_buffer_size = config.file_buffer_size;
        self.search.tuning.ignore_warnings = config.ignore_warnings;
        self.search.tuning.mutual_cover = config.mutual_cover;
        self.graph.member_cover = Some(config.member_cover);
        Ok(())
    }
    fn restore_output(&mut self, output: &Path) -> Result<(), String> {
        self.search.output = Some(output.to_path_buf());
        self.graph.output_file = Some(output.to_path_buf());
        Ok(())
    }
    fn run_search(&mut self, request: &BlockSearchRequest) -> Result<(), String> {
        let result = self.run_search_io(request);
        self.errors.capture(result)
    }
    fn greedy_vertex_cover(
        &mut self,
        database: &Path,
        edges: &Path,
        output: &Path,
    ) -> Result<(), String> {
        let result = self.greedy_vertex_cover_io(database, edges, output);
        self.errors.capture(result)
    }
}

/// Concrete native library entry for the retained Source workflow. The Source
/// CLI's commented multinode call remains inactive.
pub fn multinode_native_io(
    config: &mut MultinodeConfig,
    input: &VolumedFile,
    search: BlastpConfig,
    graph: GreedyVertexCoverConfig<Vec<u8>>,
) -> io::Result<OId> {
    let mut backend = NativeMultinodeBackend::new(search, graph);
    multinode_io(config, input, &mut backend)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::multinode::multinode::{run_block_combo, Job};
    use crate::cluster::multinode::output::Volume;
    use crate::run::main::RunBackend;
    use std::ffi::OsString;
    use std::fs;
    use std::path::PathBuf;

    fn directory(label: &str) -> PathBuf {
        let root = std::env::temp_dir().join(format!(
            "diamond-native-multinode-{label}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
        ));
        fs::create_dir_all(&root).unwrap();
        root
    }

    fn graph() -> GreedyVertexCoverConfig<Vec<u8>> {
        GreedyVertexCoverConfig {
            database: PathBuf::new(),
            edges: PathBuf::new(),
            query_or_target_cover: 99.0,
            member_cover: None,
            edge_format: b"triplet".to_vec(),
            symmetric: false,
            strict_gvc: false,
            no_gvc_reassign: false,
            connected_component_depth: Vec::new(),
            centroid_out: None,
            output_file: None,
        }
    }

    fn args(query: &Path, database: &Path, temporary: &Path) -> Vec<OsString> {
        let mut args: Vec<OsString> = [
            "diamond",
            "blastp",
            "--quiet",
            "--threads",
            "1",
            "--masking",
            "0",
            "--soft-masking",
            "none",
            "--motif-masking",
            "0",
            "--comp-based-stats",
            "0",
            "--algo",
            "1",
            "--shapes",
            "1",
        ]
        .into_iter()
        .map(Into::into)
        .collect();
        for (flag, path) in [
            ("--query", query),
            ("--db", database),
            ("--tmpdir", temporary),
        ] {
            args.push(flag.into());
            args.push(path.as_os_str().to_owned());
        }
        args
    }

    #[test]
    fn native_backend_indexes_raw_fasta_runs_search_and_consumes_real_gvc_output() {
        // Use the actual command execution guard and startup, including motif,
        // scoring/masking initialization and real temporary/input checks.
        let mut command_guard = crate::commands::dispatch::NativeBackend::new();
        command_guard.init_motif_table().unwrap();
        let root = directory("connected");
        let query = root.join("query.faa");
        let target = root.join("target.faa");
        fs::write(
            &query,
            b">query\xff\nARNDCQEGHILKMFPSTWYVARNDCQEGHILKMFPSTWYV\n",
        )
        .unwrap();
        fs::write(
            &target,
            b">target\xfe\nARNDCQEGHILKMFPSTWYVARNDCQEGHILKMFPSTWYV\n",
        )
        .unwrap();
        let argv = args(&query, &target, &root);
        let startup = crate::commands::source_config::prepare_os(&argv).unwrap();
        let mut search = crate::commands::dispatch::parse_native_blastp_config_os(
            &argv,
            b"native multinode integration",
        )
        .unwrap();
        startup.apply_search(&mut search);
        let mut backend = NativeMultinodeBackend::new(search, graph());
        let mut config = MultinodeConfig {
            root_dir: root.join("job"),
            output_file: root.join("final.tsv"),
            threads: 1,
            member_cover: 80.0,
            approx_min_id: 0.0,
            max_evalue: 1000.0,
            ..Default::default()
        };
        // Same outer Source mutation performed before constructing the Job.
        config.file_buffer_size = 65_536;
        backend.begin_workflow(&config).unwrap();
        assert_eq!(backend.search_config().file_buffer_size, 65_536);
        let errors = backend.io_errors();
        let index = root.join("query.index");
        backend.index_volume(&query, &index, &errors).unwrap();
        assert!(fs::metadata(&index).unwrap().len() > 0);
        let mut indexed = backend
            .open_indexed_records(
                &query,
                &index,
                SequenceFileFlags::NEED_LETTER_COUNT
                    | SequenceFileFlags::SEQS
                    | SequenceFileFlags::TITLES,
                &errors,
            )
            .unwrap();
        assert_eq!(indexed.letters(), 40);
        assert_eq!(indexed.read_record().unwrap().unwrap().id, b"query\xff");
        indexed.close().unwrap();
        drop(indexed);
        let mut job = Job::new(1, 1 << 20, &config.root_dir).unwrap();
        job.set_round_count(1);
        let volumes = VolumedFile(vec![
            Volume::new(query.clone(), 0, 1, 1),
            Volume::new(target.clone(), 1, 2, 1),
        ]);
        let combos = root.join("combos");
        fs::create_dir_all(&combos).unwrap();
        run_block_combo(&job, &volumes, 0, 1, &combos, &config, &mut backend).unwrap();
        let alignment = combos.join("0_1.tsv");
        let bytes = fs::read(&alignment).unwrap();
        let rows: Vec<_> = bytes
            .split(|byte| *byte == b'\n')
            .filter(|row| !row.is_empty())
            .collect();
        assert_eq!(rows.len(), 1);
        let fields: Vec<_> = rows[0].split(|byte| *byte == b'\t').collect();
        assert_eq!(fields.len(), 5);
        assert_eq!(&fields[..4], [b"0".as_slice(), b"1", b"100", b"100"]);
        assert!(std::str::from_utf8(fields[4])
            .unwrap()
            .parse::<f64>()
            .unwrap()
            .is_finite());
        assert_eq!(
            backend.search_config().tuning.algorithm.selected().unwrap(),
            Algo::DoubleIndexed
        );
        assert_eq!(backend.search_config().tuning.setup.shapes, 1);
        assert!(std::sync::Arc::ptr_eq(
            backend.search_config().scoring_matrix.as_ref().unwrap(),
            startup.matrix.as_ref().unwrap(),
        ));
        let mapping = root.join("oids.txt");
        fs::write(&mapping, b"0\n1\n").unwrap();
        let clusters = root.join("clusters.tsv");
        backend
            .greedy_vertex_cover_io(&mapping, &alignment, &clusters)
            .unwrap();
        assert_eq!(backend.graph_config().query_or_target_cover, 0.0);
        assert!(backend.graph_config().edge_format.is_empty());
        let expected: &[u8] = if cfg!(windows) {
            b"1\t0\r\n1\t1\r\n"
        } else {
            b"1\t0\n1\t1\n"
        };
        assert_eq!(fs::read(&clusters).unwrap(), expected);
        backend.restore_output(&config.output_file).unwrap();
        assert_eq!(
            backend.search_config().output.as_deref(),
            Some(config.output_file.as_path())
        );
        drop(job);
        drop(backend);
        drop(startup);
        command_guard.finish_execution();
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn native_workflow_retains_primary_index_input_error() {
        let mut command_guard = crate::commands::dispatch::NativeBackend::new();
        command_guard.init_motif_table().unwrap();
        let root = directory("typed-error");
        let query = root.join("missing.faa");
        let argv = args(&query, &query, &root);
        let startup = crate::commands::source_config::prepare_os(&argv).unwrap();
        let mut search = crate::commands::dispatch::parse_native_blastp_config_os(
            &argv,
            b"native multinode missing input",
        )
        .unwrap();
        startup.apply_search(&mut search);
        let mut config = MultinodeConfig {
            root_dir: root.join("job"),
            output_file: root.join("clusters.tsv"),
            threads: 1,
            ..Default::default()
        };
        let input = VolumedFile(vec![Volume::new(query.clone(), 0, 1, 1)]);
        let error = multinode_native_io(&mut config, &input, search, graph()).unwrap_err();
        let mut expected = b"Failed to open file: ".to_vec();
        expected.extend_from_slice(&crate::util::io::source_os_bytes(query.as_os_str()).unwrap());
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error),
            Some(expected.as_slice())
        );
        assert!(!config.output_file.exists());
        drop(startup);
        command_guard.finish_execution();
        fs::remove_dir_all(root).unwrap();
    }
}
