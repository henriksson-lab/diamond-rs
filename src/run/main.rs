//! Translation of `diamond/src/run/main.cpp`.
//!
//! The C++ entry point mutates process-global configuration and calls command
//! implementations spread across the program.  [`RunBackend`] makes those
//! boundaries explicit while [`main`] retains the dispatch, diagnostics, and
//! exit-status behavior of the original entry point.

use std::ffi::OsString;
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
    pub daa_file: std::path::PathBuf,
    pub cluster_similarity: Vec<u8>,
    pub cluster_algorithm: Option<Vec<u8>>,
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
    /// Historical Source exception: its failing opener already emitted perror.
    FileOpen,
    /// Native deferred counterpart, retaining the failing operation's captured
    /// OS diagnostic across this cloneable execution boundary.
    FileOpenDetails {
        kind: std::io::ErrorKind,
        error: crate::util::io::FileOpenError,
    },
    /// TempFile's constructor emits perror and then throws a standard error.
    /// Keep both parts across the deferred native execution boundary.
    TemporaryFileOpenDetails {
        kind: std::io::ErrorKind,
        error: crate::util::io::temp_file::TempFileOpenError,
    },
    StdioDetails {
        kind: std::io::ErrorKind,
        error: crate::util::io::cstdio_file::StdioError,
    },
    /// The Source perror part was emitted before a destructor ran. The
    /// ordinary exception catch and log remain required for standard errors.
    PerrorReported(Box<RunError>),
    Standard(String),
    Io {
        kind: std::io::ErrorKind,
        message: String,
    },
    RawIo {
        kind: std::io::ErrorKind,
        message: Vec<u8>,
    },
    Untranslated,
    /// Configuration already emitted its diagnostic or completed early.
    Status(i32),
    Unknown,
}

impl RunError {
    pub(crate) fn without_perror_marker(&self) -> &Self {
        let mut error = self;
        while let Self::PerrorReported(inner) = error {
            error = inner;
        }
        error
    }
}

pub trait RunBackend {
    /// Concrete native execution can report an exception before releasing its
    /// log, preserving Source's terminal-then-log order. Generic adapters leave
    /// reporting to this module's supplied streams.
    fn error_was_reported(&self) -> bool {
        false
    }
    /// Release invocation-owned resources after success or a reported error.
    /// Concrete native backends may be reused without retaining command locks.
    fn finish_execution(&mut self) {}
    /// Observe the real result while invocation-owned log streams still live.
    fn finish_execution_result(&mut self, _: &Result<i32, RunError>) {
        self.finish_execution();
    }
    fn init_motif_table(&mut self) -> Result<(), RunError>;
    fn parse_config(&mut self, args: &[String]) -> Result<MainConfig, RunError>;
    fn parse_config_os(&mut self, args: &[OsString]) -> Result<MainConfig, RunError> {
        let args = args
            .iter()
            .map(|s| {
                s.to_str().map(str::to_owned).ok_or_else(|| {
                    RunError::Standard(
                        "OS arguments are not UTF-8; implement parse_config_os".into(),
                    )
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.parse_config(&args)
    }
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
    fn evaluate_cluster_similarity_bytes(&mut self, expression: &[u8]) -> Result<(), RunError> {
        self.evaluate_cluster_similarity(
            std::str::from_utf8(expression).map_err(|e| RunError::Standard(e.to_string()))?,
        )
    }
    fn run_cluster(&mut self, algorithm: &str) -> Result<(), RunError>;
    fn run_cluster_bytes(&mut self, algorithm: &[u8]) -> Result<(), RunError> {
        self.run_cluster(
            std::str::from_utf8(algorithm).map_err(|e| RunError::Standard(e.to_string()))?,
        )
    }
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
    fn message_bytes(&mut self, text: &[u8]) {
        self.message(
            std::str::from_utf8(text).expect("Message is not UTF-8; implement message_bytes"),
        )
    }
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
            if config.daa_file.as_os_str().is_empty() {
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
                if let Err(error) =
                    backend.evaluate_cluster_similarity_bytes(&config.cluster_similarity)
                {
                    let mut message = b"Could not evaluate the expression: ".to_vec();
                    message.extend_from_slice(&config.cluster_similarity);
                    message.push(b'\n');
                    backend.message_bytes(&message);
                    return Err(error);
                }
            }
            backend
                .run_cluster_bytes(config.cluster_algorithm.as_deref().unwrap_or(b"cascaded"))?;
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
    main_os(
        backend,
        &args.iter().map(OsString::from).collect::<Vec<_>>(),
        features,
        out,
        err,
        log,
    )
}
pub fn main_os<B: RunBackend, O: Write, E: Write, L: Write>(
    backend: &mut B,
    args: &[OsString],
    features: RunFeatures,
    out: &mut O,
    err: &mut E,
    log: &mut L,
) -> i32 {
    let result = execute_os(backend, args, features, out, err);
    if backend.error_was_reported() {
        return match result {
            Ok(status) | Err(RunError::Status(status)) => status,
            _ => 1,
        };
    }
    let (result, perror_reported) = match result {
        Err(mut error) => {
            let mut reported = false;
            while let RunError::PerrorReported(inner) = error {
                error = *inner;
                reported = true;
            }
            (Err(error), reported)
        }
        ok => (ok, false),
    };
    match result {
        Ok(status) => status,
        Err(RunError::Status(status)) => status,
        Err(RunError::BadAlloc(message)) => {
            let _ = writeln!(err, "Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage.");
            let _ = writeln!(log, "Error: {message}");
            1
        }
        Err(RunError::FileOpen) => 1,
        Err(RunError::FileOpenDetails { error, .. }) => {
            // Source prints this in the file constructor and its dedicated
            // main catch returns silently. Native openers defer that one
            // diagnostic to the supplied terminal stream; never the log.
            if !perror_reported {
                match error.try_perror_bytes() {
                    Ok(bytes) => {
                        let _ = err.write_all(&bytes);
                    }
                    Err(error) => {
                        let _ = writeln!(err, "Error: {error}");
                    }
                }
            }
            1
        }
        Err(RunError::TemporaryFileOpenDetails { error, .. }) => {
            if !perror_reported {
                let _ = err.write_all(&error.reason);
                let _ = err.write_all(b"\n");
            }
            match error.try_message_bytes() {
                Ok(bytes) => {
                    for stream in [&mut *err as &mut dyn Write, &mut *log as &mut dyn Write] {
                        let _ = stream.write_all(b"Error: ");
                        let _ = stream.write_all(&bytes);
                        let _ = stream.write_all(b"\n");
                    }
                }
                Err(error) => {
                    let _ = writeln!(err, "Error: {error}");
                    let _ = writeln!(log, "Error: {error}");
                }
            }
            1
        }
        Err(RunError::StdioDetails { error, .. }) => {
            if !perror_reported {
                if let Some(perror) = &error.perror {
                    let _ = err.write_all(perror);
                    let _ = err.write_all(b"\n");
                }
            }
            let _ = out.flush();
            for stream in [&mut *err as &mut dyn Write, &mut *log as &mut dyn Write] {
                let _ = stream.write_all(b"Error: ");
                let _ = stream.write_all(error.message_bytes());
                let _ = stream.write_all(b"\n");
            }
            1
        }
        Err(RunError::PerrorReported(_)) => unreachable!("perror wrappers were removed"),
        Err(RunError::Standard(message) | RunError::Io { message, .. }) => {
            let _ = writeln!(err, "Error: {message}");
            let _ = writeln!(log, "Error: {message}");
            1
        }
        Err(RunError::RawIo { message, .. }) => {
            let message = message.split(|&byte| byte == 0).next().unwrap_or_default();
            for stream in [&mut *err as &mut dyn Write, &mut *log as &mut dyn Write] {
                let _ = stream.write_all(b"Error: ");
                let _ = stream.write_all(message);
                let _ = stream.write_all(b"\n");
            }
            1
        }
        Err(RunError::Untranslated) => {
            let _ = writeln!(
                err,
                "Error: Native execution of this command or option is not yet connected."
            );
            1
        }
        Err(RunError::Unknown) => {
            let _ = writeln!(err, "Exception of unknown type!");
            1
        }
    }
}

/// Prepare and execute one translated entry point, leaving error reporting to
/// its caller. The executable and safe native facade retain the same concrete
/// backend without recursively entering command setup or printing errors twice.
pub fn execute<B: RunBackend, O: Write, E: Write>(
    backend: &mut B,
    args: &[String],
    features: RunFeatures,
    out: &mut O,
    err: &mut E,
) -> Result<i32, RunError> {
    execute_os(
        backend,
        &args.iter().map(OsString::from).collect::<Vec<_>>(),
        features,
        out,
        err,
    )
}
pub fn execute_os<B: RunBackend, O: Write, E: Write>(
    backend: &mut B,
    args: &[OsString],
    features: RunFeatures,
    out: &mut O,
    err: &mut E,
) -> Result<i32, RunError> {
    let result = (|| {
        backend.init_motif_table()?;
        let config = backend.parse_config_os(args)?;
        dispatch(backend, &config, features, out, err)
    })();
    backend.finish_execution_result(&result);
    result
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
                    daa_file: std::path::PathBuf::new(),
                    cluster_similarity: Vec::new(),
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
            config.cluster_similarity = b"identity >= 50".to_vec();
            config.cluster_algorithm = Some(b"linclust".to_vec());
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

    #[test]
    fn exception_diagnostic_preserves_raw_accession_bytes_and_what_terminator() {
        let mut backend = Backend::new(RunCommand::MakeDb);
        backend.config = Err(RunError::RawIo {
            kind: std::io::ErrorKind::Other,
            message: b"Missing accession: \xff\0ignored".to_vec(),
        });
        let (mut out, mut err, mut log) = (Vec::new(), Vec::new(), Vec::new());
        assert_eq!(
            super::main(
                &mut backend,
                &[],
                RunFeatures::default(),
                &mut out,
                &mut err,
                &mut log,
            ),
            1
        );
        assert!(out.is_empty());
        assert_eq!(err, b"Error: Missing accession: \xff\n");
        assert_eq!(log, err);
    }
}
