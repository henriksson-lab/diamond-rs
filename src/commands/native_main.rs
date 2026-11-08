//! Concrete owners for the translated Source main dispatcher.
//!
//! Preparation parses argv once. The generic runtime selects a backend method;
//! that method consumes its prepared command and invokes the actual native
//! dependency, without re-entering command routing or its execution lock.
use super::*;
use crate::run::main::{MainConfig, RunBackend, RunCommand, RunError, RunFeatures};
use std::ffi::OsString;
use std::io;
use std::path::PathBuf;

pub(super) enum NativeCommand {
    Help,
    Version,
    DbInfo {
        database: PathBuf,
        file_buffer_size: usize,
    },
    MakeDb {
        inputs: Vec<PathBuf>,
        database: PathBuf,
        threads: i32,
        options: crate::commands::makedb::MakedbOptions,
        runtime: crate::commands::makedb::RuntimeOptions,
    },
    BlastP(Box<BlastpConfig>),
    BlastX(Box<BlastxConfig>),
    MakeIndex(crate::data::index::MakeIndexConfig),
    View(ViewConfig, crate::commands::view::RuntimeOptions),
    Cluster(Box<ClusterConfig>),
    #[cfg(feature = "mcl")]
    Mcl(Box<BlastpConfig>, crate::commands::mcl_cmd::MclOptions),
    MergeDaa {
        inputs: Vec<PathBuf>,
        output: PathBuf,
        buffer: usize,
    },
    /// These workflows own their parser and receive registry-validated argv.
    Workflow {
        kind: RunCommand,
        args: Vec<OsString>,
        startup: std::sync::Arc<crate::commands::source_config::Startup>,
    },
    #[cfg(feature = "extra")]
    Extra(crate::commands::extra::PreparedExtra),
    PrepDb,
    Roc,
}

pub(super) struct PreparedNative {
    pub(super) command: Option<NativeCommand>,
    pub(super) config: MainConfig,
    pub(super) threads: Option<usize>,
    pub(super) _startup: std::sync::Arc<crate::commands::source_config::Startup>,
}
pub(super) enum Preparation {
    Ready(PreparedNative),
    Finished(i32),
}

pub(super) fn ready(
    command: NativeCommand,
    kind: RunCommand,
    args: &[String],
    startup: std::sync::Arc<crate::commands::source_config::Startup>,
) -> NativeResult<Preparation> {
    ready_os(
        command,
        kind,
        &args.iter().map(OsString::from).collect::<Vec<_>>(),
        startup,
    )
}
pub(super) fn ready_os(
    mut command: NativeCommand,
    kind: RunCommand,
    args: &[OsString],
    startup: std::sync::Arc<crate::commands::source_config::Startup>,
) -> NativeResult<Preparation> {
    match &mut command {
        NativeCommand::BlastP(search) => startup.apply_search(search),
        NativeCommand::BlastX(search) => startup.apply_translated(search),
        NativeCommand::MakeDb {
            database,
            threads,
            runtime,
            ..
        } => {
            *database = startup.config.database.clone();
            *threads = startup.config.threads;
            startup.apply_tuning(&mut runtime.tuning);
        }
        NativeCommand::MakeIndex(config) => {
            config.sensitivity = startup.config.sensitivity;
            config.threads = startup.config.threads.max(0) as usize;
        }
        NativeCommand::View(config, runtime) => {
            config.daa_file = startup.config.daa_file.clone();
            config.output = startup.config.output_file.clone();
            config.threads = startup.config.threads.max(0) as usize;
            config.verbosity = startup.config.verbosity;
            config.invocation = startup.config.invocation.clone();
            runtime.database = startup.config.database.clone();
            if runtime.source_extra && !runtime.top_present {
                config.toppercent = startup.config.top;
            }
            runtime.frame_shift = startup.config.frame_shift;
            runtime.query_range_culling = startup.config.range_culling;
            runtime.search_input.sensitivity = startup.config.sensitivity;
            runtime.search_input.algo = match startup.config.algo {
                crate::basic::config::Algo::Auto => crate::run::config::Algo::Auto,
                crate::basic::config::Algo::DoubleIndexed => {
                    crate::run::config::Algo::DoubleIndexed
                }
                crate::basic::config::Algo::QueryIndexed => crate::run::config::Algo::QueryIndexed,
                crate::basic::config::Algo::CtgSeed => crate::run::config::Algo::CtgSeed,
            };
        }
        NativeCommand::Cluster(config) => {
            startup.apply_tuning(&mut config.tuning);
            config.search_args = Some(startup.argv.clone());
            config.database = startup.config.database.clone();
            config.output = startup.config.output_file.clone();
            config.threads = startup.config.threads;
            config.tmpdir = startup.config.tmpdir.clone().into();
        }
        NativeCommand::MergeDaa { output, .. } => *output = startup.config.output_file.clone(),
        #[cfg(feature = "mcl")]
        NativeCommand::Mcl(search, _) => startup.apply_search(search),
        _ => (),
    }
    let daa_file = match &command {
        NativeCommand::View(config, _) => config.daa_file.clone(),
        _ => PathBuf::new(),
    };
    let options = super::dispatch_os::Options::new(args)?;
    let cluster = matches!(
        kind,
        RunCommand::Cluster | RunCommand::DeepClust | RunCommand::LinClust
    );
    let config = MainConfig {
        command: kind,
        program_name: "diamond".into(),
        version_string: "2.1.24".into(),
        daa_file,
        cluster_similarity: if cluster {
            options.raw("cluster-similarity")?.unwrap_or_default()
        } else {
            Vec::new()
        },
        cluster_algorithm: if cluster {
            options.raw("cluster-algo")?
        } else {
            None
        },
    };
    let threads = matches!(
        kind,
        RunCommand::BlastP
            | RunCommand::BlastX
            | RunCommand::MakeDb
            | RunCommand::MakeIdx
            | RunCommand::Cluster
            | RunCommand::LinClust
            | RunCommand::DeepClust
            | RunCommand::ClusterReassign
    )
    .then_some(startup.config.threads.max(0) as usize);
    Ok(Preparation::Ready(PreparedNative {
        command: Some(command),
        config,
        threads,
        _startup: startup,
    }))
}

fn runtime_error(error: NativeDispatchError) -> RunError {
    match error {
        NativeDispatchError::Io(error) => {
            let perror_reported = crate::util::io::PerrorEmittedError::contains(&error);
            let result = if let Some(temp) = crate::util::io::temp_file::find_temp_open(&error) {
                RunError::TemporaryFileOpenDetails {
                    kind: error.kind(),
                    error: temp.clone(),
                }
            } else if let Some(open) = crate::util::io::file_open_error::find(&error) {
                RunError::FileOpenDetails {
                    kind: error.kind(),
                    error: open.clone(),
                }
            } else if let Some(stdio) = crate::util::io::cstdio_file::StdioError::find(&error) {
                RunError::StdioDetails {
                    kind: error.kind(),
                    error: stdio.clone(),
                }
            } else if let Some(allocation) = crate::util::io::BadAllocError::find(&error) {
                RunError::BadAlloc(allocation.to_string())
            } else if let Some(message) = crate::util::io::RawDiagnosticError::find(&error) {
                RunError::RawIo {
                    kind: error.kind(),
                    message: message.to_vec(),
                }
            } else {
                RunError::Io {
                    kind: error.kind(),
                    message: error.to_string(),
                }
            };
            if perror_reported {
                RunError::PerrorReported(Box::new(result))
            } else {
                result
            }
        }
        NativeDispatchError::Untranslated => RunError::Untranslated,
    }
}

/// A fresh backend owns one command and its Source logging lifetime.
pub struct NativeBackend {
    prepared: Option<PreparedNative>,
    session: crate::commands::source_config::Session,
    execution: Option<std::sync::MutexGuard<'static, ()>>,
    error_reported: bool,
}
static COMMAND_EXECUTION: std::sync::Mutex<()> = std::sync::Mutex::new(());
impl Default for NativeBackend {
    fn default() -> Self {
        Self::new()
    }
}
impl NativeBackend {
    pub fn new() -> Self {
        Self {
            prepared: None,
            session: Default::default(),
            execution: None,
            error_reported: false,
        }
    }

    fn perform(&mut self, expected: &[RunCommand]) -> Result<i32, RunError> {
        let prepared = self
            .prepared
            .as_mut()
            .ok_or_else(|| RunError::Standard("Command has not been prepared.".into()))?;
        if !expected.contains(&prepared.config.command) {
            return Err(RunError::Standard(
                "Prepared command does not match runtime dispatch.".into(),
            ));
        }
        // Leave the guard in the backend until generic dispatch has completed.
        let command = prepared.command.take().ok_or_else(|| {
            RunError::Standard("Prepared command has already been consumed.".into())
        })?;
        let session = &mut self.session;
        let execute = move || command.execute(session);
        let result = if let Some(threads) = prepared.threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(io::Error::other)
                .map_err(NativeDispatchError::from)
                .map_err(runtime_error)?
                .install(execute)
        } else {
            execute()
        };
        result.map_err(runtime_error)
    }
    fn perform_unit(&mut self, expected: &[RunCommand]) -> Result<(), RunError> {
        match self.perform(expected)? {
            0 => Ok(()),
            status => Err(RunError::Status(status)),
        }
    }
}

impl NativeCommand {
    fn execute(self, session: &mut crate::commands::source_config::Session) -> NativeResult<i32> {
        match self {
            Self::DbInfo {
                database,
                file_buffer_size,
            } => crate::commands::dbinfo::run_path_with_buffer_size(&database, file_buffer_size)?,
            Self::MakeDb {
                inputs,
                database,
                threads,
                options,
                runtime,
            } => crate::commands::makedb::run_with_runtime_paths(
                &inputs, &database, threads, &options, runtime,
            )?,
            Self::BlastP(config) => crate::commands::blastp::run(&config)?,
            Self::BlastX(config) => crate::commands::blastx::run(&config)?,
            Self::MakeIndex(config) => {
                let path = crate::data::index::make_index(&config)?;
                let mut message = b"Wrote seed index: ".to_vec();
                message.extend_from_slice(
                    &crate::util::io::path_diagnostic_bytes(&path).map_err(io::Error::other)?,
                );
                message.push(b'\n');
                crate::util::log_stream::message_stream()
                    .lock()
                    .unwrap()
                    .write_bytes(&message)?;
            }
            Self::View(config, runtime) => {
                crate::commands::view::run_with_options(&config, &runtime)?;
            }
            Self::Cluster(config) => crate::commands::cluster_cmd::run(&config)?,
            #[cfg(feature = "mcl")]
            Self::Mcl(search, options) => crate::commands::mcl_cmd::run(&search, &options)?,
            Self::MergeDaa {
                inputs,
                output,
                buffer,
            } => crate::commands::merge_daa::run_with_buffer_size_paths(&inputs, &output, buffer)?,
            Self::Workflow {
                kind,
                args,
                startup,
            } => match kind {
                RunCommand::GetSeq => {
                    crate::commands::getseq::run_with_startup_os(&args, Some(startup))?
                }
                RunCommand::GreedyVertexCover => {
                    crate::commands::greedy_vertex_cover::run_os(&args)?
                }
                RunCommand::ClusterRealign => {
                    crate::commands::realign::run_with_startup_os(&args, Some(startup))?
                }
                RunCommand::ClusterReassign => {
                    crate::commands::reassign::run_with_startup_os(&args, Some(startup))?
                }
                RunCommand::Recluster => {
                    crate::commands::recluster::run_with_startup_os(&args, Some(startup))?
                }
                RunCommand::Benchmark => {
                    crate::commands::benchmark::run_with_startup_os(&args, &startup)?
                }
                RunCommand::RegressionTest => {
                    return crate::commands::test_cmd::run_with_session_os(&args, session)
                        .map_err(Into::into);
                }
                _ => return Err(NativeDispatchError::Untranslated),
            },
            #[cfg(feature = "extra")]
            Self::Extra(prepared) => return prepared.execute().map_err(Into::into),
            Self::Help | Self::Version | Self::PrepDb | Self::Roc => {
                return Err(argument_error("Command has no backend execution phase.").into())
            }
        }
        Ok(0)
    }
}

macro_rules! methods {
    ($($method:ident => $kind:ident),* $(,)?) => {$ (
        fn $method(&mut self) -> Result<(), RunError> { self.perform_unit(&[RunCommand::$kind]) }
    )*};
}
impl RunBackend for NativeBackend {
    fn finish_execution(&mut self) {
        // Restore Source's command-global streams before admitting the next
        // caller. Both lifetimes also unwind through ordinary field drops.
        self.prepared.take();
        self.session.clear();
        self.execution.take();
    }
    fn finish_execution_result(&mut self, result: &Result<i32, RunError>) {
        // Source's catch writes the terminal diagnostic before writing its log.
        // Both effects must occur before this invocation releases its streams.
        self.error_reported = false;
        if let Err(NativeDispatchError::Io(error)) = runtime_result(result.clone()) {
            super::report_io_error(&error);
            self.error_reported = true;
        }
        let result = result.as_ref().map_err(RunError::without_perror_marker);
        let temporary_message = match result {
            Err(RunError::TemporaryFileOpenDetails { error, .. }) => Some(
                error
                    .try_message_bytes()
                    .unwrap_or_else(|conversion| conversion.to_string().into_bytes()),
            ),
            _ => None,
        };
        let message = match result {
            Err(
                RunError::Standard(message)
                | RunError::Io { message, .. }
                | RunError::BadAlloc(message),
            ) => Some(message.as_str()),
            Err(RunError::Untranslated) => {
                Some("Native execution of this command or option is not yet connected.")
            }
            Err(RunError::TemporaryFileOpenDetails { .. }) => None,
            // FileOpenException is never added to Source's log stream; its
            // single perror diagnostic belongs only to the terminal reporter.
            Err(RunError::FileOpen | RunError::FileOpenDetails { .. }) => None,
            _ => None,
        };
        if let Some(message) = message {
            let _ = crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .write(format!("Error: {message}\n"));
        }
        if let Some(message) = temporary_message {
            let mut bytes = b"Error: ".to_vec();
            bytes.extend_from_slice(&message);
            bytes.push(b'\n');
            let _ = crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .write_bytes(&bytes);
        }
        if let Err(RunError::RawIo { message, .. }) = result {
            let mut bytes = b"Error: ".to_vec();
            bytes.extend_from_slice(message.split(|&byte| byte == 0).next().unwrap_or_default());
            bytes.push(b'\n');
            let _ = crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .write_bytes(&bytes);
        }
        if let Err(RunError::StdioDetails { error, .. }) = result {
            let mut bytes = b"Error: ".to_vec();
            bytes.extend_from_slice(error.message_bytes());
            bytes.push(b'\n');
            let _ = crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .write_bytes(&bytes);
        }
        self.finish_execution();
    }
    fn error_was_reported(&self) -> bool {
        self.error_reported
    }
    fn init_motif_table(&mut self) -> Result<(), RunError> {
        self.error_reported = false;
        if self.execution.is_some() {
            return Err(RunError::Standard(
                "A native backend is already executing a command.".into(),
            ));
        }
        self.execution = Some(
            COMMAND_EXECUTION
                .lock()
                .map_err(|_| RunError::Standard("Command execution lock poisoned.".into()))?,
        );
        crate::masking::motifs::init_motif_table();
        Ok(())
    }
    fn parse_config(&mut self, args: &[String]) -> Result<MainConfig, RunError> {
        self.parse_config_os(&args.iter().map(OsString::from).collect::<Vec<_>>())
    }
    fn parse_config_os(&mut self, args: &[OsString]) -> Result<MainConfig, RunError> {
        if self.prepared.is_some() {
            return Err(RunError::Standard(
                "A native backend can execute only one command.".into(),
            ));
        }
        match super::prepare_native_os(args, &mut self.session).map_err(runtime_error)? {
            Preparation::Finished(status) => Err(RunError::Status(status)),
            Preparation::Ready(prepared) => {
                let config = prepared.config.clone();
                self.prepared = Some(prepared);
                Ok(config)
            }
        }
    }
    methods!(make_db=>MakeDb, view_daa=>View, get_seq=>GetSeq,
        random_seqs=>RandomSeqs, run_masker=>Mask, fastq2fasta=>Fastq2Fasta,
        db_info=>DbInfo, info=>Info, pairwise=>SmithWaterman,
        benchmark=>Benchmark, split=>Split, reverse=>ReverseSeqs,
        makeindex=>MakeIdx, find_shapes=>FindShapes, hash_seqs=>HashSeqs,
        list_seeds=>ListSeeds, cluster_realign=>ClusterRealign,
        greedy_vertex_cover=>GreedyVertexCover, cluster_reassign=>ClusterReassign,
        recluster=>Recluster, merge_daa=>MergeDaa, word_count=>WordCount,
        cut=>Cut, profile_recluster=>ProfileRecluster);
    fn search(&mut self) -> Result<(), RunError> {
        self.perform_unit(&[RunCommand::BlastP, RunCommand::BlastX, RunCommand::BlastN])
    }
    fn run_cluster(&mut self, algorithm: &str) -> Result<(), RunError> {
        self.run_cluster_bytes(algorithm.as_bytes())
    }
    fn run_cluster_bytes(&mut self, algorithm: &[u8]) -> Result<(), RunError> {
        let selected = self
            .prepared
            .as_ref()
            .and_then(|p| p.config.cluster_algorithm.as_deref())
            .unwrap_or(b"cascaded");
        if selected != algorithm {
            return Err(RunError::Standard(
                "Prepared clustering algorithm does not match runtime dispatch.".into(),
            ));
        }
        if algorithm != b"cascaded" && !(cfg!(feature = "mcl") && algorithm == b"mcl") {
            return Err(RunError::Standard("Clustering algorithm not found.".into()));
        }
        self.perform_unit(&[
            RunCommand::Cluster,
            RunCommand::DeepClust,
            RunCommand::LinClust,
        ])
    }
    fn regression_test(&mut self) -> Result<i32, RunError> {
        self.perform(&[RunCommand::RegressionTest])
    }
    fn evaluate_cluster_similarity(&mut self, expression: &str) -> Result<(), RunError> {
        self.evaluate_cluster_similarity_bytes(expression.as_bytes())
    }
    fn evaluate_cluster_similarity_bytes(&mut self, expression: &[u8]) -> Result<(), RunError> {
        #[cfg(feature = "mcl")]
        {
            crate::mcl::recursive_parser::evaluate_bytes(expression)
                .map(|_| ())
                .map_err(|e| runtime_error(NativeDispatchError::Io(e)))
        }
        #[cfg(not(feature = "mcl"))]
        {
            let _ = expression;
            Err(RunError::Standard("MCL support is not enabled.".into()))
        }
    }
    fn set_warning_color(&mut self) -> Result<(), RunError> {
        crate::util::system::set_color(crate::util::system::Color::Yellow, true);
        Ok(())
    }
    fn reset_warning_color(&mut self) -> Result<(), RunError> {
        crate::util::system::reset_color(true);
        Ok(())
    }
    fn message(&mut self, text: &str) {
        let _ = crate::util::log_stream::message_stream()
            .lock()
            .unwrap()
            .write(text);
    }
    fn message_bytes(&mut self, text: &[u8]) {
        let _ = crate::util::log_stream::message_stream()
            .lock()
            .unwrap()
            .write_bytes(text);
    }
}

pub(super) fn try_run(args: &[String]) -> NativeResult<i32> {
    try_run_os(&args.iter().map(OsString::from).collect::<Vec<_>>())
}
pub(super) fn try_run_os(args: &[OsString]) -> NativeResult<i32> {
    let mut backend = NativeBackend::new();
    let result = crate::run::main::execute_os(
        &mut backend,
        args,
        RunFeatures {
            with_mcl: cfg!(feature = "mcl"),
            extra: cfg!(feature = "extra"),
            with_famsa: false,
        },
        &mut crate::util::io::console::stdout_version(),
        &mut crate::util::io::console::stderr(),
    );
    let result = runtime_result(result);
    if backend.error_was_reported() {
        result.map_err(|error| match error {
            NativeDispatchError::Io(error) => {
                NativeDispatchError::Io(crate::util::io::ReportedIoError::wrap(error))
            }
            other => other,
        })
    } else {
        result
    }
}

fn runtime_result(result: Result<i32, RunError>) -> NativeResult<i32> {
    match result {
        Ok(status) | Err(RunError::Status(status)) => Ok(status),
        Err(RunError::Untranslated) => Err(NativeDispatchError::Untranslated),
        Err(RunError::TemporaryFileOpenDetails { kind, error }) => {
            Err(io::Error::new(kind, error).into())
        }
        Err(RunError::FileOpenDetails { kind, error }) => Err(io::Error::new(kind, error).into()),
        Err(RunError::StdioDetails { kind, error }) => Err(io::Error::new(kind, error).into()),
        Err(RunError::PerrorReported(error)) => match runtime_result(Err(*error)) {
            Err(NativeDispatchError::Io(error)) => {
                Err(crate::util::io::PerrorEmittedError::wrap(error).into())
            }
            other => other,
        },
        Err(RunError::Io { kind, message }) => Err(io::Error::new(kind, message).into()),
        Err(RunError::RawIo { kind, message }) => {
            Err(io::Error::new(kind, crate::util::io::RawDiagnosticError(message)).into())
        }
        Err(RunError::Standard(message)) => Err(io::Error::other(message).into()),
        Err(RunError::BadAlloc(message)) => {
            Err(crate::util::io::BadAllocError(message).into_io().into())
        }
        Err(RunError::FileOpen) => Ok(1),
        Err(RunError::Unknown) => Err(io::Error::other("Exception of unknown type!").into()),
    }
}
