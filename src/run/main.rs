//! Translation of `diamond/src/run/main.cpp`.
//!
//! The C++ entry point mutates process-global configuration and calls command
//! implementations spread across the program.  [`RunBackend`] makes those
//! boundaries explicit while [`main`] retains the dispatch, diagnostics, and
//! exit-status behavior of the original entry point.

use std::io::Write;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunCommand {
    Help,
    Version,
    MakeDb,
    BlastP,
    BlastX,
    View,
    GetSeq,
    RandomSeqs,
    Mask,
    Fastq2Fasta,
    DbInfo,
    Info,
    SmithWaterman,
    Cluster,
    DeepClust,
    LinClust,
    Benchmark,
    Split,
    RegressionTest,
    ReverseSeqs,
    Roc,
    RocId,
    MakeIdx,
    FindShapes,
    HashSeqs,
    PrepDb,
    ListSeeds,
    ClusterRealign,
    GreedyVertexCover,
    ClusterReassign,
    Recluster,
    MergeDaa,
    BlastN,
    WordCount,
    Cut,
    ProfileRecluster,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MainConfig {
    pub command: RunCommand,
    pub program_name: String,
    pub version_string: String,
    pub daa_file: String,
    pub cluster_similarity: String,
    pub cluster_algorithm: Option<String>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RunFeatures {
    pub with_mcl: bool,
    pub extra: bool,
    pub with_famsa: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RunError {
    BadAlloc(String),
    FileOpen,
    Standard(String),
    Unknown,
}

pub trait RunBackend {
    fn init_motif_table(&mut self) -> Result<(), RunError>;
    fn parse_config(&mut self, args: &[String]) -> Result<MainConfig, RunError>;
    fn make_db(&mut self) -> Result<(), RunError>;
    fn search(&mut self) -> Result<(), RunError>;
    fn view_daa(&mut self) -> Result<(), RunError>;
    fn get_seq(&mut self) -> Result<(), RunError>;
    fn random_seqs(&mut self) -> Result<(), RunError>;
    fn run_masker(&mut self) -> Result<(), RunError>;
    fn fastq2fasta(&mut self) -> Result<(), RunError>;
    fn db_info(&mut self) -> Result<(), RunError>;
    fn info(&mut self) -> Result<(), RunError>;
    fn pairwise(&mut self) -> Result<(), RunError>;
    fn evaluate_cluster_similarity(&mut self, expression: &str) -> Result<(), RunError>;
    fn run_cluster(&mut self, algorithm: &str) -> Result<(), RunError>;
    fn benchmark(&mut self) -> Result<(), RunError>;
    fn split(&mut self) -> Result<(), RunError>;
    fn regression_test(&mut self) -> Result<i32, RunError>;
    fn reverse(&mut self) -> Result<(), RunError>;
    fn makeindex(&mut self) -> Result<(), RunError>;
    fn find_shapes(&mut self) -> Result<(), RunError>;
    fn hash_seqs(&mut self) -> Result<(), RunError>;
    fn set_warning_color(&mut self) -> Result<(), RunError>;
    fn reset_warning_color(&mut self) -> Result<(), RunError>;
    fn list_seeds(&mut self) -> Result<(), RunError>;
    fn cluster_realign(&mut self) -> Result<(), RunError>;
    fn greedy_vertex_cover(&mut self) -> Result<(), RunError>;
    fn cluster_reassign(&mut self) -> Result<(), RunError>;
    fn recluster(&mut self) -> Result<(), RunError>;
    fn merge_daa(&mut self) -> Result<(), RunError>;
    fn word_count(&mut self) -> Result<(), RunError>;
    fn cut(&mut self) -> Result<(), RunError>;
    fn profile_recluster(&mut self) -> Result<(), RunError>;
    fn message(&mut self, text: &str);
}

fn write_line<W: Write>(out: &mut W, text: &str) -> Result<(), RunError> {
    writeln!(out, "{text}").map_err(|error| RunError::Standard(error.to_string()))
}

fn dispatch<B: RunBackend, O: Write, E: Write>(
    backend: &mut B,
    config: &MainConfig,
    features: RunFeatures,
    out: &mut O,
    err: &mut E,
) -> Result<i32, RunError> {
    match config.command {
        RunCommand::Help => {}
        RunCommand::Version => write_line(
            out,
            &format!("{} version {}", config.program_name, config.version_string),
        )?,
        RunCommand::MakeDb => backend.make_db()?,
        RunCommand::BlastP | RunCommand::BlastX => backend.search()?,
        RunCommand::View => {
            if config.daa_file.is_empty() {
                return Err(RunError::Standard(
                    "The view command requires a DAA (option -a) input file.".to_owned(),
                ));
            }
            backend.view_daa()?;
        }
        RunCommand::GetSeq => backend.get_seq()?,
        RunCommand::RandomSeqs => backend.random_seqs()?,
        RunCommand::Mask => backend.run_masker()?,
        RunCommand::Fastq2Fasta => backend.fastq2fasta()?,
        RunCommand::DbInfo => backend.db_info()?,
        RunCommand::Info => backend.info()?,
        RunCommand::SmithWaterman => backend.pairwise()?,
        RunCommand::Cluster | RunCommand::DeepClust | RunCommand::LinClust => {
            if features.with_mcl && !config.cluster_similarity.is_empty() {
                if let Err(error) = backend.evaluate_cluster_similarity(&config.cluster_similarity)
                {
                    backend.message(&format!(
                        "Could not evaluate the expression: {}\n",
                        config.cluster_similarity
                    ));
                    return Err(error);
                }
            }
            backend.run_cluster(config.cluster_algorithm.as_deref().unwrap_or("cascaded"))?;
        }
        RunCommand::Benchmark => backend.benchmark()?,
        RunCommand::Split => backend.split()?,
        RunCommand::RegressionTest => return backend.regression_test(),
        RunCommand::ReverseSeqs => backend.reverse()?,
        RunCommand::Roc => return Err(RunError::Standard("Deprecated command: roc".to_owned())),
        RunCommand::RocId => {
            return Err(RunError::Standard("Deprecated command: rocid".to_owned()));
        }
        RunCommand::MakeIdx => backend.makeindex()?,
        RunCommand::FindShapes => backend.find_shapes()?,
        RunCommand::HashSeqs => backend.hash_seqs()?,
        RunCommand::PrepDb => {
            backend.set_warning_color()?;
            write_line(err, "Warning: prepdb is deprecated since v2.1.14 and no longer needed to use BLAST databases. No action was taken.")?;
            backend.reset_warning_color()?;
        }
        RunCommand::ListSeeds => backend.list_seeds()?,
        RunCommand::ClusterRealign => backend.cluster_realign()?,
        RunCommand::GreedyVertexCover => backend.greedy_vertex_cover()?,
        RunCommand::ClusterReassign => backend.cluster_reassign()?,
        RunCommand::Recluster => backend.recluster()?,
        RunCommand::MergeDaa => backend.merge_daa()?,
        RunCommand::BlastN if features.extra => backend.search()?,
        RunCommand::WordCount if features.extra => backend.word_count()?,
        RunCommand::Cut if features.extra => backend.cut()?,
        RunCommand::ProfileRecluster if features.extra && features.with_famsa => {
            backend.profile_recluster()?
        }
        RunCommand::BlastN
        | RunCommand::WordCount
        | RunCommand::Cut
        | RunCommand::ProfileRecluster
        | RunCommand::Unknown => return Ok(1),
    }
    Ok(0)
}

/// Execute the translated DIAMOND entry point and return its process status.
///
/// The C++ terminate handler calls `abort`; Rust's process-level panic/abort
/// policy remains outside this library-facing dispatcher.
pub fn main<B: RunBackend, O: Write, E: Write, L: Write>(
    backend: &mut B,
    args: &[String],
    features: RunFeatures,
    out: &mut O,
    err: &mut E,
    log: &mut L,
) -> i32 {
    let result = backend
        .init_motif_table()
        .and_then(|()| backend.parse_config(args))
        .and_then(|config| dispatch(backend, &config, features, out, err));
    match result {
        Ok(status) => status,
        Err(RunError::BadAlloc(message)) => {
            let _ = writeln!(err, "Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage.");
            let _ = writeln!(log, "Error: {message}");
            1
        }
        Err(RunError::FileOpen) => 1,
        Err(RunError::Standard(message)) => {
            let _ = writeln!(err, "Error: {message}");
            let _ = writeln!(log, "Error: {message}");
            1
        }
        Err(RunError::Unknown) => {
            let _ = writeln!(err, "Exception of unknown type!");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Backend {
        config: Result<MainConfig, RunError>,
        fail: Option<RunError>,
        regression_status: i32,
        calls: Vec<String>,
        messages: Vec<String>,
    }

    impl Backend {
        fn new(command: RunCommand) -> Self {
            Self {
                config: Ok(MainConfig {
                    command,
                    program_name: "diamond".to_owned(),
                    version_string: "2.1".to_owned(),
                    daa_file: String::new(),
                    cluster_similarity: String::new(),
                    cluster_algorithm: None,
                }),
                fail: None,
                regression_status: 0,
                calls: Vec::new(),
                messages: Vec::new(),
            }
        }

        fn call(&mut self, name: &str) -> Result<(), RunError> {
            self.calls.push(name.to_owned());
            match self.fail.take() {
                Some(error) => Err(error),
                None => Ok(()),
            }
        }
    }

    macro_rules! command_methods {
        ($($method:ident),+ $(,)?) => {$(
            fn $method(&mut self) -> Result<(), RunError> {
                self.call(stringify!($method))
            }
        )+};
    }

    impl RunBackend for Backend {
        fn init_motif_table(&mut self) -> Result<(), RunError> {
            self.call("init_motif_table")
        }
        fn parse_config(&mut self, _: &[String]) -> Result<MainConfig, RunError> {
            self.calls.push("parse_config".to_owned());
            self.config.clone()
        }
        command_methods!(
            make_db,
            search,
            view_daa,
            get_seq,
            random_seqs,
            run_masker,
            fastq2fasta,
            db_info,
            info,
            pairwise,
            benchmark,
            split,
            reverse,
            makeindex,
            find_shapes,
            hash_seqs,
            set_warning_color,
            reset_warning_color,
            list_seeds,
            cluster_realign,
            greedy_vertex_cover,
            cluster_reassign,
            recluster,
            merge_daa,
            word_count,
            cut,
            profile_recluster,
        );
        fn evaluate_cluster_similarity(&mut self, expression: &str) -> Result<(), RunError> {
            self.call(&format!("evaluate:{expression}"))
        }
        fn run_cluster(&mut self, algorithm: &str) -> Result<(), RunError> {
            self.call(&format!("cluster:{algorithm}"))
        }
        fn regression_test(&mut self) -> Result<i32, RunError> {
            self.calls.push("regression_test".to_owned());
            Ok(self.regression_status)
        }
        fn message(&mut self, text: &str) {
            self.messages.push(text.to_owned());
        }
    }

    fn run(backend: &mut Backend, features: RunFeatures) -> (i32, String, String, String) {
        let mut out = Vec::new();
        let mut err = Vec::new();
        let mut log = Vec::new();
        let status = main(
            backend,
            &["diamond".to_owned()],
            features,
            &mut out,
            &mut err,
            &mut log,
        );
        (
            status,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
            String::from_utf8(log).unwrap(),
        )
    }

    #[test]
    fn version_view_and_deprecated_paths_preserve_output_and_errors() {
        let mut backend = Backend::new(RunCommand::Version);
        let result = run(&mut backend, RunFeatures::default());
        assert_eq!(
            result,
            (0, "diamond version 2.1\n".into(), "".into(), "".into())
        );

        let mut backend = Backend::new(RunCommand::View);
        let result = run(&mut backend, RunFeatures::default());
        assert_eq!(result.0, 1);
        assert_eq!(
            result.2,
            "Error: The view command requires a DAA (option -a) input file.\n"
        );
        assert_eq!(result.2, result.3);

        let mut backend = Backend::new(RunCommand::RocId);
        let result = run(&mut backend, RunFeatures::default());
        assert_eq!(result.2, "Error: Deprecated command: rocid\n");
    }

    #[test]
    fn search_cluster_and_optional_commands_follow_compile_feature_branches() {
        let mut backend = Backend::new(RunCommand::BlastX);
        assert_eq!(run(&mut backend, RunFeatures::default()).0, 0);
        assert!(backend.calls.contains(&"search".to_owned()));

        let mut backend = Backend::new(RunCommand::Cluster);
        if let Ok(config) = &mut backend.config {
            config.cluster_similarity = "identity >= 50".to_owned();
            config.cluster_algorithm = Some("linclust".to_owned());
        }
        assert_eq!(
            run(
                &mut backend,
                RunFeatures {
                    with_mcl: true,
                    ..Default::default()
                }
            )
            .0,
            0
        );
        assert!(backend
            .calls
            .contains(&"evaluate:identity >= 50".to_owned()));
        assert!(backend.calls.contains(&"cluster:linclust".to_owned()));

        let mut backend = Backend::new(RunCommand::BlastN);
        assert_eq!(run(&mut backend, RunFeatures::default()).0, 1);
        let mut backend = Backend::new(RunCommand::BlastN);
        assert_eq!(
            run(
                &mut backend,
                RunFeatures {
                    extra: true,
                    ..Default::default()
                }
            )
            .0,
            0
        );
    }

    #[test]
    fn prepdb_and_regression_return_preserve_special_behavior() {
        let mut backend = Backend::new(RunCommand::PrepDb);
        let result = run(&mut backend, RunFeatures::default());
        assert_eq!(result.0, 0);
        assert_eq!(result.2, "Warning: prepdb is deprecated since v2.1.14 and no longer needed to use BLAST databases. No action was taken.\n");
        assert_eq!(
            backend.calls[2..],
            ["set_warning_color", "reset_warning_color"]
        );

        let mut backend = Backend::new(RunCommand::RegressionTest);
        backend.regression_status = 17;
        assert_eq!(run(&mut backend, RunFeatures::default()).0, 17);
    }

    #[test]
    fn exception_categories_match_cpp_catch_diagnostics() {
        let mut backend = Backend::new(RunCommand::MakeDb);
        backend.fail = Some(RunError::BadAlloc("oom detail".to_owned()));
        let result = run(&mut backend, RunFeatures::default());
        assert_eq!(result.0, 1);
        assert!(result
            .2
            .starts_with("Failed to allocate sufficient memory."));
        assert_eq!(result.3, "Error: oom detail\n");

        let mut backend = Backend::new(RunCommand::MakeDb);
        backend.config = Err(RunError::FileOpen);
        let result = run(&mut backend, RunFeatures::default());
        assert_eq!(result, (1, "".into(), "".into(), "".into()));

        let mut backend = Backend::new(RunCommand::MakeDb);
        backend.config = Err(RunError::Unknown);
        let result = run(&mut backend, RunFeatures::default());
        assert_eq!(result.2, "Exception of unknown type!\n");
        assert!(result.3.is_empty());
    }
}
