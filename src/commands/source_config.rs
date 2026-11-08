//! Concrete optional-command startup for the translated basic/config.cpp.
//! Configuration owns its scoring/masking/translation state rather than
//! installing mutable process globals in a library caller.
use crate::util::io::console::{source_eprint, source_eprintln, source_print};

use crate::basic::config::{Algo, Command, Config, ConfigBackend, ConfigInput};
use crate::basic::translate::Translator;
use crate::config::Sensitivity;
use crate::masking::{Masking, MaskingConfig};
use crate::stats::{cbs::CbsMode, score_matrix::ScoreMatrix};
use crate::util::command_line_parser::ReadOption;
use std::ffi::OsString;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Actual Search::run directory resolution/diagnostic. Source evaluates
/// TempFile::get_temp_dir even when its message stream is disabled, before
/// opening the query/reference/output owners.
pub(crate) fn begin_search_temp(config: &super::blastp::BlastpConfig) -> io::Result<()> {
    // Source inserts the prefix before evaluating get_temp_dir. Retain it on
    // the console and in diamond.log even when that real owner fails to open.
    crate::util::log_stream::message_stream()
        .lock()
        .unwrap()
        .write("Temporary directory: ")?;
    let directory = crate::util::io::temp_file::source_temp_directory_path(
        &config.tmpdir,
        config.output.as_deref().unwrap_or_else(|| Path::new("")),
    );
    let resolved = crate::util::io::TempFile::get_temp_dir_path_in_with_no_unlink(
        &directory,
        config.tuning.no_unlink,
    )
    .map_err(io::Error::other)?;
    crate::util::log_stream::message_stream()
        .lock()
        .unwrap()
        .write_bytes(
            &crate::util::io::source_os_bytes(resolved.as_os_str()).map_err(io::Error::other)?,
        )?
        .write("\n")?;
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MatrixKey {
    name: Vec<u8>,
    file: PathBuf,
    gap_open: i32,
    gap_extend: i32,
    frame_shift: i32,
    stop: i32,
    scale: i32,
}
impl MatrixKey {
    fn new(
        name: &[u8],
        gap_open: i32,
        gap_extend: i32,
        frame_shift: i32,
        tuning: &super::blastp::SearchTuning,
    ) -> Self {
        Self {
            name: name.to_ascii_uppercase(),
            file: tuning.matrix_file.clone(),
            gap_open,
            gap_extend,
            frame_shift,
            stop: tuning.stop_match_score,
            scale: tuning.matrix_scale,
        }
    }
}
struct Resources {
    key: MatrixKey,
    matrix: Option<Arc<ScoreMatrix>>,
    masking: Option<Arc<Masking>>,
    translator: Translator,
}

/// Initialized immutable dependencies owned by one Source configuration.
/// Cloned search configurations share these kernels without process globals.
#[derive(Clone, Default)]
pub struct Runtime(Option<Arc<Resources>>);
impl std::fmt::Debug for Runtime {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SourceRuntime")
            .field("initialized", &self.0.is_some())
            .finish()
    }
}
impl PartialEq for Runtime {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}
impl Runtime {
    pub fn tantan(&self) -> Option<Arc<crate::masking::tantan::TantanMasker>> {
        self.masking().map(Masking::tantan_shared)
    }
    pub fn masking(&self) -> Option<&Masking> {
        self.0.as_ref()?.masking.as_deref()
    }
    pub fn translator(&self) -> Option<&Translator> {
        Some(&self.0.as_ref()?.translator)
    }
    pub(crate) fn matrix_for(
        &self,
        name: &[u8],
        gap_open: i32,
        gap_extend: i32,
        frame_shift: i32,
        tuning: &super::blastp::SearchTuning,
    ) -> Option<Arc<ScoreMatrix>> {
        let resources = self.0.as_ref()?;
        (resources.key == MatrixKey::new(name, gap_open, gap_extend, frame_shift, tuning))
            .then(|| resources.matrix.clone())
            .flatten()
    }
}

#[derive(Default)]
pub(crate) struct Session {
    logging: Option<super::search_logging::SearchLogging>,
}
impl Session {
    pub(crate) fn clear(&mut self) {
        self.logging.take();
    }
}

pub(crate) struct Startup {
    pub argv: Vec<OsString>,
    pub config: Config,
    pub matrix: Option<Arc<ScoreMatrix>>,
    pub masking: Option<Arc<Masking>>,
    pub _translator: Translator,
    pub runtime: Runtime,
    _owned_session: Option<Session>,
}
impl Startup {
    pub(crate) fn apply_tuning(&self, tuning: &mut super::blastp::SearchTuning) {
        tuning.runtime = self.runtime.clone();
        tuning.store_query_quality = self.config.store_query_quality;
        let requested = match self.config.algo {
            Algo::Auto => crate::run::config::Algo::Auto,
            Algo::DoubleIndexed => crate::run::config::Algo::DoubleIndexed,
            Algo::QueryIndexed => crate::run::config::Algo::QueryIndexed,
            Algo::CtgSeed => crate::run::config::Algo::CtgSeed,
        };
        tuning.algorithm = super::search_algorithm::AlgorithmControl::default();
        tuning.algorithm.requested = requested;
        tuning
            .algorithm
            .begin_search(self.config.sensitivity, requested);
    }
    pub(crate) fn apply_search(&self, search: &mut super::blastp::BlastpConfig) {
        self.apply_tuning(&mut search.tuning);
        search.scoring_matrix = self.matrix.clone();
        search.database = self.config.database.clone();
        search.output = (!self.config.output_file.as_os_str().is_empty())
            .then(|| self.config.output_file.clone());
        search.frame_shift = self.config.frame_shift;
        search.range_culling = self.config.range_culling;
        search.toppercent = self.config.top;
        search.sensitivity = self.config.sensitivity;
        search.threads = self.config.threads;
        search.tmpdir = self.config.tmpdir.clone().into();
        search.query_indexed = self.config.algo == Algo::QueryIndexed;
        search.contiguous_seed_mode = self.config.algo == Algo::CtgSeed;
        search.sam.invocation = self.config.invocation.clone();
    }
    pub(crate) fn apply_translated(&self, search: &mut super::blastx::BlastxConfig) {
        self.apply_tuning(&mut search.tuning);
        search.scoring_matrix = self.matrix.clone();
        search.database = self.config.database.clone();
        search.output = (!self.config.output_file.as_os_str().is_empty())
            .then(|| self.config.output_file.clone());
        search.frame_shift = self.config.frame_shift;
        search.range_culling = self.config.range_culling;
        search.toppercent = self.config.top;
        search.sensitivity = self.config.sensitivity;
        search.threads = self.config.threads;
        search.tmpdir = self.config.tmpdir.clone().into();
        search.query_indexed = self.config.algo == Algo::QueryIndexed;
        search.contiguous_seed_mode = self.config.algo == Algo::CtgSeed;
        search.sam.invocation = self.config.invocation.clone();
    }
}
struct Backend<'a> {
    session: &'a mut Session,
    errors: &'a crate::util::io::IoErrorOwner,
    no_logfile: bool,
    matrix_name: Vec<u8>,
    tuning: super::blastp::SearchTuning,
    masking_config: MaskingConfig,
    matrix: Option<ScoreMatrix>,
    masking: Option<Masking>,
    translator: Option<Translator>,
    cbs: CbsMode,
    raw_algo: Vec<u8>,
    raw_dbtype: Vec<u8>,
    #[cfg(feature = "dna")]
    raw_dna_extension: Vec<u8>,
    algo: Option<Algo>,
}
impl ConfigBackend for Backend<'_> {
    fn invocation(&mut self, args: &[OsString]) -> Result<Vec<u8>, String> {
        self.errors
            .capture(crate::basic::config::invocation_bytes_io(args))
    }
    fn auxiliary_formats(&mut self, input: &ConfigInput) -> Result<(), String> {
        if !input.aligned_file.as_os_str().is_empty() {
            let mut message = b"Aligned file format: ".to_vec();
            message.extend_from_slice(&input.alfmt);
            message.push(b'\n');
            self.errors.capture(
                crate::util::log_stream::log_stream()
                    .lock()
                    .unwrap()
                    .write_bytes(&message)
                    .map(|_| ()),
            )?;
        }
        Ok(())
    }

    fn presentation(&mut self, input: &ConfigInput, verbosity: u8) -> Result<(), String> {
        let command = self.errors.capture(command_name(&input.argv))?;
        if verbosity >= 1 || input.command == Command::RegressionTest {
            let header = format!("diamond v{}.{} (C) Max Planck Society for the Advancement of Science, Benjamin J. Buchfink, University of Tuebingen\nDocumentation, support and updates available at http://www.diamondsearch.org\nPlease cite: http://dx.doi.org/10.1038/s41592-021-01101-x Nature Methods (2021)\n\n",
                crate::basic::consts::VERSION_STRING, crate::basic::consts::BUILD_VERSION);
            if command == "help" {
                source_print!("{header}");
            } else {
                source_eprint!("{header}");
            }
        }
        self.errors.capture(
            crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .write(format!(
                    "diamond v{}.{}\n",
                    crate::basic::consts::VERSION_STRING,
                    crate::basic::consts::BUILD_VERSION
                ))
                .map(|_| ()),
        )?;
        if input.argc == 2 && !matches!(input.command, Command::Version | Command::RegressionTest) {
            crate::util::io::console::stdout_endls_fmt(format_args!(
                "{}",
                super::argument_registry::documentation(command)
            ));
        }
        Ok(())
    }
    fn streams(&mut self, input: &ConfigInput, verbosity: u8) -> Result<(), String> {
        if self.session.logging.is_some() {
            super::search_logging::SearchLogging::reconfigure_source(verbosity, !self.no_logfile);
        } else {
            self.session.logging = Some(self.errors.capture(
                super::search_logging::SearchLogging::establish_source(verbosity, !self.no_logfile),
            )?);
        }
        let invocation = self.invocation(&input.argv)?;
        self.errors.capture(
            crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .write_bytes(&invocation)
                .and_then(|stream| stream.write("\n"))
                .map(|_| ()),
        )?;
        if cfg!(debug_assertions) {
            self.errors.capture(
                crate::util::log_stream::verbose_stream()
                    .lock()
                    .unwrap()
                    .write("Assertions enabled.\n")
                    .map(|_| ()),
            )?;
        }
        let command = self.errors.capture(command_name(&input.argv))?;
        if input.argc != 2
            && matches!(
                command,
                "makedb"
                    | "blastp"
                    | "blastx"
                    | "blastn"
                    | "view"
                    | "cluster"
                    | "deepclust"
                    | "linclust"
                    | "test"
                    | "compute-medoids"
                    | "reassign"
                    | "greedy-vertex-cover"
                    | "recluster"
            )
        {
            let threads = input
                .threads
                .filter(|threads| *threads != 0)
                .unwrap_or(input.hardware_threads);
            self.errors.capture(
                crate::util::log_stream::message_stream()
                    .lock()
                    .unwrap()
                    .write(format!("#CPU threads: {threads}\n"))
                    .map(|_| ()),
            )?;
        }
        Ok(())
    }
    fn completed(&mut self, config: &Config, _: &ConfigInput) -> Result<(), String> {
        if config.old_warning {
            crate::util::system::set_color(crate::util::system::Color::Yellow, true);
            crate::util::log_stream::flush_tied_stdout();
            source_eprintln!("Warning: This version of DIAMOND is more than 180 days old. It is recommended to always use the latest version.");
            crate::util::system::reset_color(false);
        }
        self.errors.capture(
            crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .write(format!(
                    "TRACE_PT_MEMBUF={}\nMAX_SHAPE_LEN=19 SEQ_MASK STRICT_BAND\n",
                    if config.trace_pt_membuf { "ON" } else { "OFF" }
                ))
                .map(|_| ()),
        )?;
        Ok(())
    }
    fn cbs(&mut self, code: u32) -> Result<(), String> {
        self.cbs = CbsMode::parse(&code.to_string());
        Ok(())
    }
    fn scoring(&mut self, input: &ConfigInput) -> Result<(), String> {
        // BasicConfig passes to_upper_case(config.matrix) to the builtin
        // constructor, whose diagnostic name preserves that argument. The
        // config literal itself remains unchanged for archive metadata.
        let constructor_name = self.matrix_name.to_ascii_uppercase();
        self.matrix = Some(self.errors.capture(
            super::blastp::native_score_matrix_with_tuning_bytes(
                &constructor_name,
                input.gap_open,
                input.gap_extend,
                input.frame_shift,
                &self.tuning,
                0,
            ),
        )?);
        if input.argc != 2 {
            let matrix = self.matrix.as_ref().unwrap();
            let number = crate::output::format::format_default_float;
            let mut message = b"Scoring parameters: (Matrix=".to_vec();
            if self.tuning.matrix_file.as_os_str().is_empty() {
                message.extend_from_slice(&constructor_name);
            } else {
                message.extend_from_slice(matrix.name().as_bytes());
            }
            message.extend_from_slice(
                format!(
                    " Lambda={} K={} Penalties={}/{})\n",
                    number(matrix.lambda()),
                    number(matrix.k()),
                    matrix.gap_open(),
                    matrix.gap_extend(),
                )
                .as_bytes(),
            );
            self.errors.capture(
                crate::util::log_stream::message_stream()
                    .lock()
                    .unwrap()
                    .write_bytes(&message)
                    .map(|_| ()),
            )?;
        }
        Ok(())
    }
    fn masking(&mut self) -> Result<(), String> {
        self.masking = Some(Masking::new(
            self.matrix
                .as_ref()
                .ok_or("Scoring matrix was not initialized")?,
            self.masking_config.clone(),
        ));
        Ok(())
    }
    fn enums(&mut self) -> Result<(), String> {
        self.algo = Some(match self.raw_algo.as_slice() {
            b"" => Algo::Auto,
            b"0" => Algo::DoubleIndexed,
            b"1" => Algo::QueryIndexed,
            b"ctg" => Algo::CtgSeed,
            value => return self.errors.capture(Err(enum_error(value, b", 0, 1, ctg"))),
        });
        if !matches!(self.raw_dbtype.as_slice(), b"prot" | b"nucl") {
            return self
                .errors
                .capture(Err(enum_error(&self.raw_dbtype, b"nucl, prot")));
        }
        #[cfg(feature = "dna")]
        if !matches!(self.raw_dna_extension.as_slice(), b"ksw" | b"wfa") {
            return self
                .errors
                .capture(Err(enum_error(&self.raw_dna_extension, b"ksw, wfa")));
        }
        Ok(())
    }
    fn translator(&mut self, code: u32) -> Result<(), String> {
        self.translator = Some(Translator::init(code)?);
        Ok(())
    }
    fn temp(&mut self, path: &str) -> Result<(), String> {
        self.temp_path(Path::new(path))
    }
    fn temp_path(&mut self, path: &Path) -> Result<(), String> {
        self.errors.capture(
            crate::util::io::TempFileHandler::new()
                .init_path(path)
                .map_err(io::Error::other),
        )?;
        self.errors.capture(
            crate::util::log_stream::verbose_stream()
                .lock()
                .unwrap()
                .write(format!(
                    "CPU features detected: {}\n",
                    crate::util::simd::features()
                ))
                .map(|_| ()),
        )?;
        self.errors.capture(
            crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .write(format!(
                    "L3 cache size: {}\n",
                    crate::util::system::l3_cache_size()
                ))
                .map(|_| ()),
        )?;
        Ok(())
    }
    fn mkdir(&mut self, path: &str) -> Result<(), String> {
        self.mkdir_path(Path::new(path))
    }
    fn mkdir_path(&mut self, path: &Path) -> Result<(), String> {
        self.errors.capture(crate::util::system::mkdir_path(path))
    }
}
fn enum_error(value: &[u8], permitted: &[u8]) -> io::Error {
    let mut message = b"Invalid value for string field: ".to_vec();
    message.extend_from_slice(value);
    message.extend_from_slice(b". Permitted values: ");
    message.extend_from_slice(permitted);
    crate::util::io::RawDiagnosticError(message).into_io()
}
struct Values {
    items: super::argument_registry::OsValues,
    errors: crate::util::io::IoErrorOwner,
}
impl std::ops::Deref for Values {
    type Target = super::argument_registry::OsValues;
    fn deref(&self) -> &Self::Target {
        &self.items
    }
}
fn scalar<T: ReadOption>(values: &Values, name: &str, mut default: T) -> Result<T, String> {
    if let Some(value) = values.get(name) {
        values
            .errors
            .capture(T::read_option_os(&mut default, value))?;
    }
    Ok(default)
}

fn command_name(args: &[OsString]) -> io::Result<&'static str> {
    Ok(
        super::argument_registry::command_name_os(args)?.unwrap_or(if args.len() < 2 {
            "help"
        } else {
            ""
        }),
    )
}
pub(super) fn prepare(args: &[String]) -> io::Result<Arc<Startup>> {
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    prepare_os(&args)
}
pub(crate) fn prepare_os(args: &[OsString]) -> io::Result<Arc<Startup>> {
    prepare_with_io_os(args, true)
}
pub(super) fn prepare_without_io(args: &[String]) -> io::Result<Arc<Startup>> {
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    prepare_without_io_os(&args)
}
pub(super) fn prepare_without_io_os(args: &[OsString]) -> io::Result<Arc<Startup>> {
    prepare_with_io_os(args, false)
}
fn prepare_with_io_os(args: &[OsString], check_io: bool) -> io::Result<Arc<Startup>> {
    let mut session = Session::default();
    let mut startup = prepare_inner_os(args, &mut session, check_io)?;
    startup._owned_session = Some(session);
    Ok(Arc::new(startup))
}
pub(super) fn prepare_in(args: &[String], session: &mut Session) -> io::Result<Arc<Startup>> {
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    prepare_os_in(&args, session)
}
pub(super) fn prepare_os_in(args: &[OsString], session: &mut Session) -> io::Result<Arc<Startup>> {
    prepare_inner_os(args, session, true).map(Arc::new)
}
pub(super) fn prepare_without_io_in(
    args: &[String],
    session: &mut Session,
) -> io::Result<Arc<Startup>> {
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    prepare_without_io_os_in(&args, session)
}
pub(super) fn prepare_without_io_os_in(
    args: &[OsString],
    session: &mut Session,
) -> io::Result<Arc<Startup>> {
    prepare_inner_os(args, session, false).map(Arc::new)
}
fn prepare_inner_os(
    args: &[OsString],
    session: &mut Session,
    check_io: bool,
) -> io::Result<Startup> {
    let errors = crate::util::io::IoErrorOwner::default();
    let values = Values {
        items: super::argument_registry::values_os(args)?,
        errors: errors.clone(),
    };
    errors.finish(prepare_values_os(args, session, check_io, &values))
}
fn prepare_values_os(
    args: &[OsString],
    session: &mut Session,
    check_io: bool,
    values: &Values,
) -> Result<Startup, String> {
    let command_name = values.errors.capture(command_name(args))?;
    let command = match command_name {
        "makedb" => Command::MakeDb,
        "blastp" => Command::BlastP,
        "blastx" => Command::BlastX,
        "view" => Command::View,
        "dbinfo" => Command::DbInfo,
        "version" => Command::Version,
        "test" => Command::RegressionTest,
        "cluster" => Command::Cluster,
        "deepclust" => Command::DeepClust,
        "linclust" => Command::LinClust,
        "recluster" => Command::Recluster,
        "reassign" => Command::ClusterReassign,
        "realign" => Command::ClusterRealign,
        "benchmark" => Command::Benchmark,
        "getseq" => Command::GetSeq,
        "mask" => Command::Mask,
        "fastq2fasta" => Command::Fastq2Fasta,
        "listseeds" => Command::ListSeeds,
        "compute-medoids" => Command::ComputeMedoids,
        "model-seqs" => Command::ModelSeqs,
        "make-seed-table" => Command::MakeSeedTable,
        "blastn" => Command::BlastN,
        _ => Command::Other,
    };
    let mut input = ConfigInput {
        argc: args.len(),
        argv: args.to_vec(),
        command,
        check_io,
        hardware_threads: crate::util::system::hardware_concurrency() as i32,
        now: std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |time| time.as_secs() as i64),
        ..ConfigInput::default()
    };
    macro_rules! bind {
        ($field:ident, $name:literal) => {
            input.$field = scalar(&values, $name, input.$field)?;
        };
    }
    macro_rules! flag {
        ($field:ident, $name:literal) => {
            input.$field = values.contains_key($name);
        };
    }
    bind!(database, "db");
    bind!(output_file, "out");
    bind!(daa_file, "daa");
    bind!(compression, "compress");
    bind!(chunk_size, "block-size");
    bind!(frame_shift, "frameshift");
    bind!(global_ranking, "global-ranking");
    bind!(taxon_k, "taxon-k");
    bind!(comp_based_stats, "comp-based-stats");
    bind!(matrix_file, "custom-matrix");
    bind!(gap_open, "gapopen");
    bind!(gap_extend, "gapextend");
    bind!(query_gencode, "query-gencode");
    bind!(strand, "strand");
    bind!(unfmt, "unfmt");
    bind!(alfmt, "alfmt");
    bind!(aligned_file, "al");
    bind!(lowmem, "index-chunks");
    bind!(parallel_tmpdir, "parallel-tmpdir");
    bind!(tmpdir, "tmpdir");
    flag!(debug_log, "log");
    flag!(quiet, "quiet");
    flag!(verbose, "verbose");
    flag!(no_auto_append, "no-auto-append");
    flag!(no_self_hits, "no-self-hits");
    flag!(long_reads, "long-reads");
    flag!(range_culling, "range-culling");
    flag!(multiprocessing, "multiprocessing");
    flag!(mp_init, "mp-init");
    flag!(mp_recover, "mp-recover");
    flag!(store_query_quality, "store-query-quality");
    flag!(target_indexed, "target-indexed");
    flag!(swipe_all, "swipe");
    flag!(hit_membuf, "hit-membuf");
    input.translated_cbs_supported =
        CbsMode::parse(&input.comp_based_stats.to_string()).support_translated();
    input.top_present = values.contains_key("top");
    input.max_target_present = values.contains_key("max-target-seqs");
    if input.top_present {
        input.top = Some(scalar(&values, "top", 0.0f64)?);
    }
    if values.contains_key("threads") {
        input.threads = Some(scalar(&values, "threads", 0i32)?);
    }
    input.query_file = scalar(values, "query", Vec::<PathBuf>::new())?;
    input.output_format = scalar(values, "outfmt", Vec::<Vec<u8>>::new())?;
    for (name, sensitivity) in [
        ("faster", Sensitivity::Faster),
        ("fast", Sensitivity::Fast),
        ("mid-sensitive", Sensitivity::MidSensitive),
        ("sensitive", Sensitivity::Sensitive),
        ("more-sensitive", Sensitivity::MoreSensitive),
        ("very-sensitive", Sensitivity::VerySensitive),
        ("ultra-sensitive", Sensitivity::UltraSensitive),
        ("shapes-30x10", Sensitivity::Shapes30x10),
        ("shapes-6x10", Sensitivity::Shapes6x10),
        ("linclust-20", Sensitivity::Linclust20),
    ] {
        if values.contains_key(name) {
            input.sensitivities.push(sensitivity);
        }
    }
    let mut backend = Backend {
        session,
        errors: &values.errors,
        no_logfile: values.contains_key("no-logfile"),
        matrix_name: scalar(&values, "matrix", b"blosum62".to_vec())?,
        tuning: values
            .errors
            .capture(super::dispatch::parse_source_config_tuning_os(args))?,
        masking_config: MaskingConfig {
            tantan_min_mask_prob: scalar(&values, "tantan-minMaskProb", 0.9f64)? as f32,
            max_motif_len: scalar(&values, "max-motif-len", 30i32)?.max(0) as usize,
        },
        matrix: None,
        masking: None,
        translator: None,
        cbs: CbsMode::Hauser,
        raw_algo: scalar(&values, "algo", Vec::<u8>::new())?,
        raw_dbtype: scalar(&values, "dbtype", b"prot".to_vec())?,
        #[cfg(feature = "dna")]
        raw_dna_extension: scalar(&values, "extension", b"ksw".to_vec())?,
        algo: None,
    };
    let mut config = if cfg!(feature = "extra") {
        Config::new_extra(&input, &mut backend)?
    } else {
        Config::new(&input, &mut backend)?
    };
    if !input.swipe_all {
        config.algo = backend.algo.ok_or("Algorithm was not initialized")?;
    }
    let key = MatrixKey::new(
        &backend.matrix_name,
        input.gap_open,
        input.gap_extend,
        config.frame_shift,
        &backend.tuning,
    );
    let matrix = backend.matrix.map(Arc::new);
    let masking = backend.masking.map(Arc::new);
    let translator = backend.translator.ok_or("Translator was not initialized")?;
    let runtime = Runtime(Some(Arc::new(Resources {
        key,
        matrix: matrix.clone(),
        masking: masking.clone(),
        translator: translator.clone(),
    })));
    Ok(Startup {
        argv: input.argv.clone(),
        config,
        matrix,
        masking,
        _translator: translator,
        runtime,
        _owned_session: None,
    })
}
