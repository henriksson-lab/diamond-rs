//! Native command parsing and concrete command dispatch shared by the CLI and library.
use crate::util::io::console::{source_eprintln, source_println};

#[derive(Debug)]
pub enum NativeDispatchError {
    Io(std::io::Error),
    Untranslated,
}
type NativeResult<T> = Result<T, NativeDispatchError>;
impl From<std::io::Error> for NativeDispatchError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}
fn argument_error(error: impl std::fmt::Display) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidInput, error.to_string())
}

/// Reject capacities that cannot fit the target or cannot be reserved now.
/// This does not guarantee that later allocations succeed: configured input,
/// output and codec buffers have independent lifetimes and allocations.
fn preflight_file_buffer_size(value: i64) -> Option<usize> {
    let capacity = usize::try_from(value).ok()?;
    let mut allocation = Vec::<u8>::new();
    allocation.try_reserve_exact(capacity).ok()?;
    Some(capacity)
}

/// Run connected commands natively and return a status without exiting the process.
/// Untranslated commands/options report a native error; they never invoke C++.
pub fn run(args: &[&str]) -> i32 {
    let args: Vec<String> = args.iter().map(|value| (*value).to_owned()).collect();
    finish(try_run(&args))
}

/// Report captured Source I/O diagnostics through the process or embedded CLI.
/// FileOpen has its own catch; TempFile combines perror with a standard error.
pub fn report_io_error(error: &std::io::Error) -> i32 {
    if crate::util::io::ReportedIoError::contains(error) {
        return 1;
    }
    if let Some(temp) = crate::util::io::temp_file::find_temp_open(error) {
        // TempFile emits perror without a cout tie, followed by the standard
        // exception catch's tied cerr diagnostic. The OS reason is captured.
        if let Err(conversion) = crate::util::io::emit_source_perror(error) {
            return report_io_error(&std::io::Error::other(conversion));
        }
        crate::util::log_stream::flush_tied_stdout();
        let message = match temp.try_message_bytes() {
            Ok(message) => message,
            Err(conversion) => return report_io_error(&std::io::Error::other(conversion)),
        };
        use std::io::Write;
        let mut stderr = crate::util::io::console::stderr();
        let _ = stderr.write_all(b"Error: ");
        let _ = stderr.write_all(&message);
        let _ = stderr.write_all(b"\n");
    } else if crate::util::io::file_open_error::find(error).is_some() {
        // Source emits perror at the failing constructor, then its dedicated
        // main catch returns without another message. perror has no cerr tie.
        if let Err(conversion) = crate::util::io::emit_source_perror(error) {
            return report_io_error(&std::io::Error::other(conversion));
        }
    } else if let Some(stdio) = crate::util::io::cstdio_file::StdioError::find(error) {
        use std::io::Write;
        if let Err(conversion) = crate::util::io::emit_source_perror(error) {
            return report_io_error(&std::io::Error::other(conversion));
        }
        crate::util::log_stream::flush_tied_stdout();
        let mut stderr = crate::util::io::console::stderr();
        let _ = stderr.write_all(b"Error: ");
        let _ = stderr.write_all(stdio.message_bytes());
        let _ = stderr.write_all(b"\n");
    } else if crate::util::io::BadAllocError::find(error).is_some() {
        crate::util::log_stream::flush_tied_stdout();
        source_eprintln!("Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage.");
    } else if let Some(message) = crate::util::io::RawDiagnosticError::find(error) {
        use std::io::Write;
        crate::util::log_stream::flush_tied_stdout();
        let mut stderr = crate::util::io::console::stderr();
        let _ = stderr.write_all(b"Error: ");
        let _ = stderr.write_all(message);
        let _ = stderr.write_all(b"\n");
    } else {
        crate::util::log_stream::flush_tied_stdout();
        source_eprintln!("Error: {error}");
    }
    1
}

fn finish(result: NativeResult<i32>) -> i32 {
    match result {
        Ok(code) => code,
        Err(NativeDispatchError::Io(error)) => report_io_error(&error),
        Err(NativeDispatchError::Untranslated) => {
            crate::util::log_stream::flush_tied_stdout();
            source_eprintln!(
                "Error: Native execution of this command or option is not yet connected."
            );
            1
        }
    }
}

use std::path::PathBuf;

use crate::cluster::cascaded::GraphAlgo;
use crate::commands::blastp::BlastpConfig;
use crate::commands::blastx::BlastxConfig;
use crate::commands::cluster_cmd::{ClusterConfig, ClusterWorkflow};
use crate::commands::view::ViewConfig;
use crate::config::Sensitivity;
use crate::output::format::FieldId;
use crate::stats::cbs::CbsMode;

/// Execute native command routing, leaving untranslated routes to the caller.
/// Only an explicit oracle request can use the executable's optional C++
/// bridge; this function and the public `diamond::run` API never call C++.
#[path = "dispatch_os.rs"]
mod dispatch_os;
pub(crate) use dispatch_os::{
    parse_float_search_option_os, parse_native_blastp_config_os, parse_native_search_tuning_os,
    parse_search_string_option_os, parse_signed_search_option_os, parse_source_config_tuning_os,
};

#[path = "native_main.rs"]
mod native_main;
pub use native_main::NativeBackend;

pub fn try_run(args: &[String]) -> NativeResult<i32> {
    native_main::try_run(args)
}
pub fn try_run_os(args: &[std::ffi::OsString]) -> NativeResult<i32> {
    native_main::try_run_os(args)
}
pub fn run_os(args: &[std::ffi::OsString]) -> i32 {
    finish(try_run_os(args))
}

fn prepare_native(
    args: &[String],
    session: &mut super::source_config::Session,
) -> NativeResult<native_main::Preparation> {
    prepare_native_os(
        &args
            .iter()
            .map(std::ffi::OsString::from)
            .collect::<Vec<_>>(),
        session,
    )
}
fn prepare_native_os(
    args: &[std::ffi::OsString],
    session: &mut super::source_config::Session,
) -> NativeResult<native_main::Preparation> {
    dispatch_os::prepare(args, session)
}

pub(crate) fn print_banner() {
    source_eprintln!(
        "diamond v2.1.24.178 (C) Max Planck Society for the Advancement of Science, Benjamin J. Buchfink, University of Tuebingen"
    );
    source_eprintln!(
        "Documentation, support and updates available at http://www.diamondsearch.org"
    );
    source_eprintln!(
        "Please cite: http://dx.doi.org/10.1038/s41592-021-01101-x Nature Methods (2021)"
    );
    source_eprintln!();
}

fn print_usage() {
    source_println!("diamond v2.1.24 — Rust port");
    source_println!();
    source_println!("Commands:");
    source_println!("  makedb     Build DIAMOND database from FASTA");
    source_println!("  blastp     Protein-protein alignment");
    source_println!("  blastx     Translated DNA-protein alignment");
    source_println!("  view       View DAA file");
    source_println!("  dbinfo     Print database info");
    source_println!("  getseq     Retrieve sequences from database");
    source_println!("  makeidx    Build persistent seed index");
    source_println!("  cluster    Cluster sequences");
    source_println!("  greedy-vertex-cover  Cluster a supplied edge graph");
    source_println!("  merge-daa  Merge DAA files");
    source_println!("  version    Show version");
    source_println!("  test       Run regression tests");
    source_println!();
    source_println!("Use 'diamond help' for options, or 'diamond blastp --help' for search help.");
    source_println!(
        "Native blastp/blastx: --memory-limit 16G is the default soft RSS ceiling; 0G forces disk mode."
    );
}

pub(crate) fn source_output_format_description() -> String {
    let base_fields = [
        "qseqid",
        "qlen",
        "sseqid",
        "sallseqid",
        "slen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "qseq",
        "sseq",
        "evalue",
        "bitscore",
        "score",
        "length",
        "pident",
        "nident",
        "mismatch",
        "positive",
        "gapopen",
        "gaps",
        "ppos",
        "qframe",
        "btop",
        "staxids",
        "sscinames",
        "sskingdoms",
        "stitle",
        "salltitles",
        "qcovhsp",
        "qtitle",
        "full_sseq",
        "qqual",
        "qnum",
        "snum",
        "scovhsp",
        "full_qqual",
        "full_qseq",
        "qseq_gapped",
        "sseq_gapped",
        "qstrand",
        "cigar",
        "skingdoms",
        "sphylums",
        "full_qseq_mate",
        "qseq_translated",
        "hspnum",
        "normalized_bitscore",
        "normalized_nident",
        "approx_pident",
        "corrected_bitscore",
        "slineages",
    ];
    let mut fields = base_fields
        .iter()
        .map(|key| {
            let field = FieldId::from_name(key)
                .and_then(FieldId::output_field)
                .expect("documented output field must be registered");
            (field.key.to_string(), field.description.to_string())
        })
        .collect::<Vec<_>>();
    for rank in crate::data::taxonomy_nodes::Rank::NAMES.iter().skip(1) {
        fields.push((
            format!("s{}", rank.replace(' ', "_")),
            format!("Unique subject {rank}(s), separated by a ';'"),
        ));
    }
    let width = fields.iter().map(|(key, _)| key.len()).max().unwrap_or(0);
    let mut format_description = String::from(
        "output format\n\
\t0   = BLAST pairwise\n\
\t5   = BLAST XML\n\
\t6   = BLAST tabular\n\
\t100 = DIAMOND alignment archive (DAA)\n\
\t101 = SAM\n\
\t102 = Taxonomic classification\n\
\t103 = PAF\n\
\t104 = JSON (flat)\n\n\
\tValues 6 and 104 may be followed by a space-separated list of these keywords:\n\n",
    );
    for (key, description) in fields {
        format_description.push_str(&format!("\t{key:>width$} means {description}\n"));
    }
    format_description.push_str(
        "\n\tDefault: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\n",
    );
    format_description
}

fn get_arg(args: &[String], flags: &[&str]) -> Option<String> {
    for (i, arg) in args.iter().enumerate().rev() {
        if flags.contains(&arg.as_str()) {
            return args.get(i + 1).cloned();
        }
        for flag in flags {
            if flag.starts_with("--") {
                let Some(value) = arg.strip_prefix(&format!("{flag}=")) else {
                    continue;
                };
                return Some(value.to_string());
            }
            if flag.starts_with('-') && flag.len() == 2 {
                let Some(value) = arg.strip_prefix(flag) else {
                    continue;
                };
                if !value.is_empty() {
                    return Some(value.to_string());
                }
            }
        }
    }
    None
}

#[cfg(test)]
fn get_all_args(args: &[String], flags: &[&str]) -> Vec<String> {
    let mut values = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let mut joined = None;
        if flags.contains(&args[i].as_str()) {
            i += 1;
        } else {
            for flag in flags {
                if flag.starts_with("--") {
                    joined = args[i]
                        .strip_prefix(&format!("{flag}="))
                        .map(ToOwned::to_owned);
                } else if flag.starts_with('-') && flag.len() == 2 {
                    joined = args[i]
                        .strip_prefix(flag)
                        .filter(|value| !value.is_empty())
                        .map(ToOwned::to_owned);
                }
                if joined.is_some() {
                    break;
                }
            }
            if joined.is_none() {
                i += 1;
                continue;
            }
            i += 1;
        }
        if let Some(value) = joined {
            values.push(value);
        }
        while i < args.len() && !is_cpp_option_start(&args[i]) {
            values.push(args[i].clone());
            i += 1;
        }
    }
    values
}

/// Return the values belonging to the last occurrence of an option.
///
/// DIAMOND's parser assigns repeated vector-valued options instead of
/// appending them, including query inputs and output-field lists.
fn get_last_arg_group(args: &[String], flags: &[&str]) -> Vec<String> {
    let mut last = None;
    let mut i = 0;
    while i < args.len() {
        let mut joined = None;
        if flags.contains(&args[i].as_str()) {
            i += 1;
        } else {
            for flag in flags {
                if flag.starts_with("--") {
                    joined = args[i]
                        .strip_prefix(&format!("{flag}="))
                        .map(ToOwned::to_owned);
                } else if flag.starts_with('-') && flag.len() == 2 {
                    joined = args[i]
                        .strip_prefix(flag)
                        .filter(|value| !value.is_empty())
                        .map(ToOwned::to_owned);
                }
                if joined.is_some() {
                    break;
                }
            }
            if joined.is_none() {
                i += 1;
                continue;
            }
            i += 1;
        }

        let mut values = Vec::new();
        if let Some(value) = joined {
            values.push(value);
        }
        while i < args.len() && !is_cpp_option_start(&args[i]) {
            values.push(args[i].clone());
            i += 1;
        }
        last = Some(values);
    }
    last.unwrap_or_default()
}

/// Parse a vector-valued option using DIAMOND's last-occurrence-wins rule,
/// while rejecting an occurrence without any values as the upstream parser
/// does for options whose minimum argument count is one.
fn get_last_arg_group_checked(args: &[String], flags: &[&str]) -> Result<Vec<String>, String> {
    let values = get_last_arg_group(args, flags);
    let present = args.iter().any(|arg| {
        flags.iter().any(|flag| {
            arg == flag
                || (flag.starts_with("--") && arg.starts_with(&format!("{flag}=")))
                || (flag.starts_with('-')
                    && flag.len() == 2
                    && arg.starts_with(flag)
                    && arg.len() > flag.len())
        })
    });
    if present && values.is_empty() {
        Err(format!("Option {} requires at least one value.", flags[0]))
    } else {
        Ok(values)
    }
}

fn parse_search_iterate(args: &[String]) -> Result<Option<Vec<String>>, String> {
    parse_search_vector(args, "iterate")
}

fn parse_search_header(args: &[String]) -> Result<Option<Vec<String>>, String> {
    parse_search_vector(args, "header")
}

fn parse_search_vector(args: &[String], name: &str) -> Result<Option<Vec<String>>, String> {
    let mut result = None;
    for (index, argument) in args.iter().enumerate().skip(2) {
        if argument.starts_with(&format!("--{name}=")) {
            return Err(format!("Invalid option: {}", &argument[2..]));
        }
        if argument == &format!("--{name}") {
            result = Some(
                args[index + 1..]
                    .iter()
                    .take_while(|value| !is_cpp_option_start(value))
                    .cloned()
                    .collect(),
            );
        }
    }
    Ok(result)
}

pub(crate) fn parse_search_string_option(
    args: &[String],
    name: &str,
) -> Result<Option<String>, String> {
    let flag = format!("--{name}");
    let mut result = None;
    for (index, argument) in args.iter().enumerate().skip(2) {
        if argument.starts_with(&format!("{flag}=")) {
            return Err(format!("Invalid option: {}", &argument[2..]));
        }
        if argument != &flag {
            continue;
        }
        let values = args[index + 1..]
            .iter()
            .take_while(|value| !is_cpp_option_start(value))
            .collect::<Vec<_>>();
        if values.len() != 1 {
            return Err(format!("Invalid parameter count for option '{flag}'"));
        }
        result = Some(values[0].clone());
    }
    Ok(result)
}

/// Common native search options, also used by the explicit DNA command.
/// The caller retains workflow-specific Source configuration validation.
pub(crate) fn parse_native_blastp_config(
    args: &[String],
    invocation: &str,
) -> std::io::Result<BlastpConfig> {
    parse_native_blastp_config_os(
        &args
            .iter()
            .map(std::ffi::OsString::from)
            .collect::<Vec<_>>(),
        invocation.as_bytes(),
    )
}

fn parse_native_search_output(args: &[String]) -> Result<(Option<String>, Vec<String>), String> {
    let daa_file = get_arg(args, &["-a", "--daa"]);
    let mut output = get_arg(args, &["-o", "--out"]);
    let mut outfmt = get_last_arg_group(args, &["-f", "--outfmt"]);

    if let Some(format) = outfmt.first() {
        if !matches!(format.as_str(), "null" | "edge")
            && crate::output::format::FormatCode::parse(format).is_none()
        {
            return Err(format!("Invalid output format: {format}\nAllowed values: 0,5,xml,6,tab,100,daa,101,sam,102,103,104,paf"));
        }
    }

    if let Some(mut daa_file) = daa_file {
        if output.is_some() {
            return Err("Options --daa and --out cannot be used together.".to_string());
        }
        if outfmt
            .first()
            .is_some_and(|format| !format.eq_ignore_ascii_case("daa"))
        {
            return Err(
                "Invalid parameter: --daa/-a. Output file is specified with the --out/-o parameter."
                    .to_string(),
            );
        }
        if !has_flag(args, "--no-auto-append")
            && !daa_file.is_empty()
            && !daa_file.ends_with(".daa")
        {
            daa_file.push_str(".daa");
        }
        output = Some(daa_file);
        if outfmt.is_empty() {
            outfmt.push("daa".to_string());
        }
    } else if outfmt
        .first()
        .is_some_and(|format| format == "100" || format.eq_ignore_ascii_case("daa"))
    {
        let path = output
            .as_mut()
            .ok_or_else(|| "DAA output requires -a/--daa or -o/--out.".to_string())?;
        if !has_flag(args, "--no-auto-append") && !path.is_empty() && !path.ends_with(".daa") {
            path.push_str(".daa");
        }
    }

    let compression = parse_search_string_option(args, "compress")?.unwrap_or_default();
    if !compression.is_empty()
        && outfmt
            .first()
            .is_some_and(|format| matches!(format.as_str(), "100" | "daa"))
    {
        return Err("Compression is not supported for DAA format.".into());
    }
    if !has_flag(args, "--no-auto-append") && matches!(compression.as_str(), "1" | "zstd") {
        if let Some(path) = output.as_mut().filter(|path| !path.is_empty()) {
            let suffix = if compression == "1" { ".gz" } else { ".zst" };
            if !path.ends_with(suffix) {
                path.push_str(suffix);
            }
        }
    }

    if outfmt
        .first()
        .is_some_and(|format| matches!(format.as_str(), "6" | "tab" | "104" | "json-flat"))
    {
        for name in outfmt.iter().skip(1) {
            if FieldId::from_name(name)
                .and_then(FieldId::output_field)
                .is_none()
            {
                return Err(format!("Invalid output field: {name}"));
            }
        }
        crate::output::format::TabularFormat::new(
            &outfmt.iter().map(String::as_str).collect::<Vec<_>>(),
            matches!(outfmt[0].as_str(), "104" | "json-flat"),
            parse_frame_shift(args)?,
            args[1] == "blastx",
            &get_arg(args, &["--matrix"])
                .unwrap_or_else(|| "BLOSUM62".into())
                .to_ascii_uppercase(),
        )?;
    }
    Ok((output, outfmt))
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|a| a == flag)
}

/// Whether an option occurs in separated (`--foo value`), long joined
/// (`--foo=value`), or short joined (`-b2`) form. This deliberately checks
/// presence rather than successfully parsing a value: an option unsupported by
/// the native pipeline must never be silently ignored just because it is
/// malformed.
fn has_option(args: &[String], flags: &[&str]) -> bool {
    args.iter().any(|arg| {
        flags.iter().any(|flag| {
            arg == flag
                || (flag.starts_with("--")
                    && arg
                        .strip_prefix(flag)
                        .is_some_and(|suffix| suffix.starts_with('=')))
                || (flag.starts_with('-')
                    && flag.len() == 2
                    && arg
                        .strip_prefix(flag)
                        .is_some_and(|suffix| !suffix.is_empty()))
        })
    })
}

#[derive(Debug, Clone, PartialEq)]
struct NativeViewCli {
    daa_file: PathBuf,
    output: PathBuf,
    outfmt: Vec<Vec<u8>>,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    forward_only: bool,
    header: Vec<Vec<u8>>,
    threads: usize,
    qnum_offset: u64,
    snum_offset: u64,
    verbosity: u8,
    file_buffer_size: usize,
    runtime: super::view::RuntimeOptions,
}

impl Default for NativeViewCli {
    fn default() -> Self {
        Self {
            daa_file: PathBuf::new(),
            output: PathBuf::new(),
            outfmt: Vec::new(),
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            header: Vec::new(),
            threads: 0,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
            file_buffer_size: 67_108_864,
            runtime: super::view::RuntimeOptions::default(),
        }
    }
}

/// Match upstream's `atoll`-based parsing for signed ordinal offsets, but make
/// its result deterministic across Rust targets: leading whitespace and a
/// sign are accepted, parsing stops at the first non-digit, no digits means
/// zero, and overflow saturates to the signed 64-bit limits.
fn parse_cpp_atoll(value: &str) -> i64 {
    let mut value = value.as_bytes().to_vec();
    value.push(0);
    unsafe { libc::atoll(value.as_ptr().cast()) as i64 }
}
fn parse_cpp_atof(value: &str) -> f64 {
    let mut value = value.as_bytes().to_vec();
    value.push(0);
    unsafe { libc::atof(value.as_ptr().cast()) }
}

fn is_cpp_option_start(value: &str) -> bool {
    let bytes = value.as_bytes();
    bytes.len() > 1 && bytes[0] == b'-' && !bytes[1].is_ascii_digit()
}

fn unsupported_view_option(id: &str) -> bool {
    matches!(
        id,
        "query"
            | "q"
            | "db"
            | "d"
            | "evalue"
            | "e"
            | "matrix"
            | "gapopen"
            | "gapextend"
            | "unal"
            | "compress"
            | "sam-query-len"
            | "masking"
            | "comp-based-stats"
            | "id"
            | "query-cover"
            | "subject-cover"
            | "block-size"
            | "b"
            | "index-chunks"
            | "c"
            | "sensitive"
            | "very-sensitive"
            | "ultra-sensitive"
    )
}

/// Parse the native view surface with the same grouping rules as DIAMOND's
/// CommandLineParser: a switch owns every following token until the next
/// non-numeric `-` token, scalar options require exactly one value, vector
/// options are replaced by their last occurrence, and joined long options are
/// invalid while joined short options are accepted.
fn parse_native_view_cli(args: &[String]) -> Result<NativeViewCli, String> {
    let args = args
        .iter()
        .map(std::ffi::OsString::from)
        .collect::<Vec<_>>();
    dispatch_os::view_cli(
        &dispatch_os::Options::new(&args).map_err(|e| e.to_string())?,
        args.len(),
    )
    .map_err(|e| e.to_string())
}
#[cfg(feature = "extra")]
fn parse_extra_native_view_cli(args: &[String]) -> Result<NativeViewCli, String> {
    parse_native_view_cli(args)
}
fn parse_baseline_native_view_cli(args: &[String]) -> Result<NativeViewCli, String> {
    parse_native_view_cli(args)
}

fn validate_native_cluster_options(args: &[String]) -> Result<(), String> {
    let parallel = get_arg(args, &["--parallel-tmpdir"]).unwrap_or_default();
    if !parallel.is_empty() && !args.get(1).is_some_and(|command| command == "linclust") {
        return Err("Option is not permitted for this workflow: --parallel-tmpdir".to_string());
    }
    if parallel.is_empty() && get_arg(args, &["--aln-out"]).is_some_and(|path| !path.is_empty()) {
        return Err(
            "Alignment output to file is not supported without --parallel-tmpdir.".to_string(),
        );
    }
    Ok(())
}

fn parse_optional_cluster_number<T: crate::util::command_line_parser::ReadOption + Default>(
    args: &[String],
    flags: &[&str],
) -> Result<Option<T>, String> {
    match get_arg(args, flags) {
        Some(value) => {
            // Source scalar numeric readers use atoi/atoll/atof, including
            // prefix conversion and zero for nonnumeric values.
            let mut parsed = T::default();
            T::read_option(&mut parsed, &[value])?;
            Ok(Some(parsed))
        }
        None if has_option(args, flags) => Err(format!("Option {} requires a value.", flags[0])),
        None => Ok(None),
    }
}

fn parse_optional_cluster_string(
    args: &[String],
    flags: &[&str],
) -> Result<Option<String>, String> {
    match get_arg(args, flags) {
        Some(value) => Ok(Some(value)),
        None if has_option(args, flags) => Err(format!("Option {} requires a value.", flags[0])),
        None => Ok(None),
    }
}

#[cfg(feature = "mcl")]
fn parse_native_mcl_options(args: &[String]) -> Result<super::mcl_cmd::MclOptions, String> {
    let args = args
        .iter()
        .map(std::ffi::OsString::from)
        .collect::<Vec<_>>();
    dispatch_os::mcl_options(&dispatch_os::Options::new(&args).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())
}
fn parse_native_cluster_config(args: &[String]) -> Result<ClusterConfig, String> {
    let args = args
        .iter()
        .map(std::ffi::OsString::from)
        .collect::<Vec<_>>();
    dispatch_os::cluster_config(
        &dispatch_os::Options::new(&args).map_err(|e| e.to_string())?,
        &args,
    )
    .map_err(|e| e.to_string())
}

#[path = "source_no_effect.rs"]
mod source_no_effect;

/// Native search accepts only connected controls and proven Source no-effect
/// registrations. Other requests produce an explicit native error.
fn has_unknown_search_option(args: &[String], blastx: bool) -> bool {
    const COMMON: &[&str] = &[
        "--cbs-window",
        "--max-motif-len",
        // Source parses these settings but only consumes them in other
        // workflows, or stores them in constructor-local deprecated scalars.
        "--K",
        "--lambda",
        "--bootstrap",
        "--check-multi-target",
        "--cut-bar",
        "--family-map",
        "--family-map-query",
        "--gapped-xdrop",
        "--hit-band",
        "--hit-score",
        "--linclust-chunk-size",
        "--log-evalue-scale",
        "--oid-output",
        "--rank-ratio2",
        "--roc-file",
        "--ungapped-score",
        "--window",
        "--tantan-minMaskProb",
        "--qnum-offset",
        "--snum-offset",
        "--xml-blord-format",
        "--no-parse-seqids",
        "--query",
        "--query-gencode",
        "--strand",
        "--min-orf",
        "--swipe",
        "--sam-query-len",
        "--include-lineage",
        "--no-auto-append",
        "--salltitles",
        "--sallseqid",
        "--db",
        "--out",
        "--daa",
        "--matrix",
        "--gapopen",
        "--gapextend",
        "--evalue",
        "--max-target-seqs",
        "--ext-chunk-size",
        "--top",
        "--global-ranking",
        "--frameshift",
        "--long-reads",
        "--quiet",
        "--verbose",
        "--log",
        "--no-logfile",
        "--ext",
        "--id",
        "--threads",
        "--outfmt",
        "--faster",
        "--fast",
        "--mid-sensitive",
        "--linclust-20",
        "--shapes-6x10",
        "--shapes-30x10",
        "--sensitive",
        "--more-sensitive",
        "--very-sensitive",
        "--ultra-sensitive",
        "--masking",
        "--motif-masking",
        "--min-query-len",
        "--query-cover",
        "--subject-cover",
        "--comp-based-stats",
        "--no-self-hits",
        "--xdrop",
        "--memory-limit",
        "--target-indexed",
        "--tmpdir",
        "--legacy",
        "--help",
        // Recognised incompatibilities are included here because the explicit
        // routing checks below decide whether a value is natively supported.
        "--max-hsps",
        "--unal",
        "--block-size",
        "--index-chunks",
        "--algo",
        "--custom-matrix",
        "--stop-match-score",
        "--dbsize",
        "--no-ranking",
        "--no-reextend",
        "--mapany",
        "--self",
        "--add-self-aln",
        "--target_hard_cap",
        "--score-drop-factor",
        "--ranking-cutoff-bitscore",
        "--load-balancing",
        "--kmer-ranking",
        "--log-query",
        "--log-extend",
        "--target-parallel-verbosity",
        "--band-bin",
        "--col-bin",
        "--palign",
        "--min_task_trace_pts",
        "--trace-pt-fetch-size",
        "--chaining-maxgap",
        "--chaining-range-cover",
        "--chaining-maxnodes",
        "--chaining-len-cap",
        "--chaining-min-nodes",
        "--no_chaining_merge_hsps",
        "--chaining-stacked-hsp-ratio",
        "--aln-out",
        "--cluster-threshold",
        "--hamming-ext",
        "--diag-filter-id",
        "--diag-filter-cov",
        "--mutual-cover",
        "--radix-cluster-buffered",
        "--join-split-size",
        "--join-split-key-len",
        "--radix-bits",
        "--join-ht-factor",
        "--hash_join_swap",
        "--min-len-ratio",
        "--rank-ratio",
        "--rank-factor",
        "--cutoff-score-8bit",
        "--max-swipe-dp",
        "--approx-backtrace",
        "--transcript-len-estimate",
        "--swipe-chunk-size",
        "--cbs-err-tolerance",
        "--cbs-it-limit",
        "--anchor-window",
        "--anchor-score",
        "--cbs-matrix-scale",
        "--heartbeat",
        "--unaligned-targets",
        "--load-threads",
        "--minichunk",
        "--linsearch",
        "--lin-stage1",
        "--lin-combo",
        "--linclust-banded-ext",
        "--tile-size",
        "--swipe-task-size",
        "--anchored-swipe",
        "--band",
        "--shape-mask",
        "--shapes",
        "--seed-cut",
        "--freq-masking",
        "--freq-sd",
        "--id2",
        "--bin",
        "--sketch-size",
        "--ungapped-evalue",
        "--ungapped-evalue-short",
        "--gapped-filter-evalue",
        "--gapped-filter-diag-score",
        "--short-query-ungapped-bitscore",
        "--short-query-max-len",
        "--store-query-quality",
        "--ungapped-window",
        "--left-most-interval",
        "--gapped-filter-evalue1",
        "--gapped-filter-window",
        "--query-match-distance-threshold",
        "--length-ratio-threshold",
        "--cbs-angle",
        "--un",
        "--al",
        "--unfmt",
        "--alfmt",
        "--taxonlist",
        "--taxon-exclude",
        "--seqidlist",
        "--skip-missing-seqids",
        "--range-culling",
        "--range-cover",
        "--culling-overlap",
        "--taxon-k",
        "--query-or-subject-cover",
        "--no-unlink",
        "--ignore-warnings",
        "--no-block-size-limit",
        "--hit-membuf",
        "--multiprocessing",
        "--mp-init",
        "--mp-recover",
        "--mp-query-chunk",
        "--mp-self",
        "--parallel-tmpdir",
        "--minimizer-window",
        "--soft-masking",
        "--approx-id",
        "--file-buffer-size",
        "--query-parallel-limit",
        "--min-score",
        "--header",
        "--compress",
        "--iterate",
        "--narrow-band-cov",
        "--narrow-band-factor",
        "--min-band-overlap",
    ];
    const BLASTX: &[&str] = &[];
    const SHORT_COMMON: &[&str] = &[
        "-q", "-d", "-o", "-a", "-e", "-k", "-p", "-f", "-x", "-t", "-b", "-c", "-g", "-h", "-F",
        "-s", "-v", "-l", "-w", "-X", "-M", "-z",
    ];

    args.iter().skip(2).any(|arg| {
        if let Some(long) = arg.strip_prefix("--") {
            let name = &arg[..2 + long.find('=').unwrap_or(long.len())];
            !COMMON.contains(&name)
                && !source_no_effect::contains(name)
                && !(blastx && BLASTX.contains(&name))
        } else if arg.starts_with('-') && arg.as_bytes().get(1).is_some_and(u8::is_ascii_alphabetic)
        {
            let name = &arg[..2];
            (name == "-h" && arg.len() != 2)
                || (!SHORT_COMMON.contains(&name) && !(blastx && name == "-F"))
        } else {
            false
        }
    })
}

/// Source locals in basic/config.cpp:489-500 have no consumer. Other listed
/// fields are read only in test/ROC/clustering, never the search workflow.
/// Preserve their typed reads instead of sending legitimate ignored flags to
/// the optional C++ diagnostic adapter.
fn parse_search_ignored_options(args: &[String]) -> Result<(), String> {
    for name in [
        "K",
        "lambda",
        "gapped-xdrop",
        "log-evalue-scale",
        "rank-ratio2",
    ] {
        parse_float_search_option(args, name, 0.0)?;
    }
    for name in ["hit-band", "hit-score", "ungapped-score", "window"] {
        parse_unsigned_search_option(args, name, 0)?;
    }
    for name in [
        "family-map",
        "family-map-query",
        "linclust-chunk-size",
        "roc-file",
    ] {
        parse_search_string_option(args, name)?;
    }
    Ok(())
}

fn route_common_search_to_legacy(_args: &[String], _query_translated: bool) -> bool {
    // All Source output formats and algorithms have native implementations.
    // Invalid values and incompatible combinations belong to their native
    // parsers/configuration checks, not the optional conformance adapter.
    false
}

fn parse_search_soft_masking(args: &[String]) -> Result<crate::masking::MaskingAlgo, String> {
    let value = parse_search_string_option(args, "soft-masking")?.unwrap_or_default();
    if value.is_empty() {
        return Ok(crate::masking::MaskingAlgo::None);
    }
    crate::commands::reassign::parse_soft_masking(&value)
}

fn route_blastp_to_legacy(args: &[String]) -> bool {
    has_unknown_search_option(args, false) || route_common_search_to_legacy(args, false)
}

fn route_blastx_to_legacy(args: &[String]) -> bool {
    has_unknown_search_option(args, true) || route_common_search_to_legacy(args, true)
}

fn parse_arg_or<T: crate::util::command_line_parser::ReadOption>(
    args: &[String],
    flags: &[&str],
    mut default: T,
) -> T {
    if let Some(value) = get_arg(args, flags) {
        // Use the translated scalar readers, including Source's C numeric
        // prefix conversions, instead of replacing malformed numbers with
        // the option's default value.
        let _ = T::read_option(&mut default, &[value]);
    }
    default
}

fn parse_block_size(args: &[String]) -> Result<Option<f64>, String> {
    let mut size = None;
    for (index, argument) in args.iter().enumerate().skip(2) {
        let attached = if argument == "--block-size" || argument == "-b" {
            None
        } else if let Some(value) = argument.strip_prefix("-b") {
            Some(value)
        } else if argument.starts_with("--block-size=") {
            return Err(format!("Invalid option: {}", &argument[2..]));
        } else {
            continue;
        };
        let mut values: Vec<&str> = attached.into_iter().collect();
        values.extend(
            args[index + 1..]
                .iter()
                .take_while(|value| !is_cpp_option_start(value))
                .map(String::as_str),
        );
        if values.len() != 1 {
            return Err("Invalid parameter count for option '-b/--block-size'".into());
        }
        size = Some(parse_cpp_atof(values[0]));
    }
    Ok(size)
}

fn parse_search_algorithm(args: &[String]) -> Result<Option<String>, String> {
    let mut algorithm = None;
    for (index, argument) in args.iter().enumerate().skip(2) {
        if argument.starts_with("--algo=") {
            return Err(format!("Invalid option: {}", &argument[2..]));
        }
        if argument != "--algo" {
            continue;
        }
        let values: Vec<_> = args[index + 1..]
            .iter()
            .take_while(|value| !is_cpp_option_start(value))
            .collect();
        if values.len() != 1 {
            return Err("Invalid parameter count for option '--algo'".into());
        }
        algorithm = Some(values[0].clone());
    }
    if let Some(value) = algorithm.as_deref() {
        if !matches!(value, "" | "0" | "1" | "ctg") {
            return Err(format!(
                "Invalid value for string field: {value}. Permitted values: , 0, 1, ctg"
            ));
        }
    }
    Ok(algorithm)
}

fn parse_index_chunks(args: &[String]) -> Result<u32, String> {
    let mut chunks = 0;
    for (index, argument) in args.iter().enumerate().skip(2) {
        let attached = if argument == "--index-chunks" || argument == "-c" {
            None
        } else if let Some(value) = argument.strip_prefix("-c") {
            Some(value)
        } else if argument.starts_with("--index-chunks=") {
            return Err(format!("Invalid option: {}", &argument[2..]));
        } else {
            continue;
        };
        let mut values: Vec<&str> = attached.into_iter().collect();
        values.extend(
            args[index + 1..]
                .iter()
                .take_while(|value| !is_cpp_option_start(value))
                .map(String::as_str),
        );
        if values.len() != 1 {
            return Err("Invalid parameter count for option '-c/--index-chunks'".into());
        }
        chunks = parse_cpp_atoll(values[0]) as u32;
    }
    Ok(chunks)
}

fn parse_global_ranking(args: &[String]) -> Result<i64, String> {
    let mut count = 0;
    for (index, argument) in args.iter().enumerate().skip(2) {
        let attached = if argument == "--global-ranking" || argument == "-g" {
            None
        } else if let Some(value) = argument.strip_prefix("-g") {
            Some(value)
        } else if argument.starts_with("--global-ranking=") {
            return Err(format!("Invalid option: {}", &argument[2..]));
        } else {
            continue;
        };
        let mut values: Vec<&str> = attached.into_iter().collect();
        values.extend(
            args[index + 1..]
                .iter()
                .take_while(|value| !is_cpp_option_start(value))
                .map(String::as_str),
        );
        if values.len() != 1 {
            return Err("Invalid parameter count for option '-g/--global-ranking'".into());
        }
        count = parse_cpp_atoll(values[0]);
    }
    Ok(count)
}

fn parse_frame_shift(args: &[String]) -> Result<i32, String> {
    let mut count = 0;
    for (index, argument) in args.iter().enumerate().skip(2) {
        let attached = if argument == "--frameshift" || argument == "-F" {
            None
        } else if let Some(value) = argument.strip_prefix("-F") {
            Some(value)
        } else if argument.starts_with("--frameshift=") {
            return Err(format!("Invalid option: {}", &argument[2..]));
        } else {
            continue;
        };
        let mut values: Vec<&str> = attached.into_iter().collect();
        values.extend(
            args[index + 1..]
                .iter()
                .take_while(|value| !is_cpp_option_start(value))
                .map(String::as_str),
        );
        if values.len() != 1 {
            return Err("Invalid parameter count for option '-F/--frameshift'".into());
        }
        count = parse_cpp_atoll(values[0]) as i32;
    }
    Ok(if count == 0 && has_flag(args, "--long-reads") {
        15
    } else {
        count
    })
}

pub(crate) fn parse_signed_search_option(
    args: &[String],
    name: &str,
    default: i32,
) -> Result<i32, String> {
    Ok(parse_search_string_option(args, name)?
        .map_or(default, |value| parse_cpp_atoll(&value) as i32))
}

pub(crate) fn parse_native_search_tuning(
    args: &[String],
) -> Result<crate::commands::blastp::SearchTuning, String> {
    parse_native_search_tuning_inner(args, true)
}

pub(crate) fn parse_source_config_tuning(
    args: &[String],
) -> Result<crate::commands::blastp::SearchTuning, String> {
    parse_native_search_tuning_inner(args, false)
}

fn parse_native_search_tuning_inner(
    args: &[String],
    resolve: bool,
) -> Result<crate::commands::blastp::SearchTuning, String> {
    let args = args
        .iter()
        .map(std::ffi::OsString::from)
        .collect::<Vec<_>>();
    if resolve {
        parse_native_search_tuning_os(&args)
    } else {
        parse_source_config_tuning_os(&args)
    }
    .map_err(|e| e.to_string())
}

fn parse_auxiliary_output(
    args: &[String],
) -> Result<crate::commands::query_output::AuxiliaryOutputOptions, String> {
    Ok(crate::commands::query_output::AuxiliaryOutputOptions {
        unaligned: PathBuf::from(parse_search_string_option(args, "un")?.unwrap_or_default()),
        aligned: PathBuf::from(parse_search_string_option(args, "al")?.unwrap_or_default()),
        unaligned_format: parse_search_string_option(args, "unfmt")?
            .unwrap_or_else(|| "fasta".into())
            .into_bytes(),
        aligned_format: parse_search_string_option(args, "alfmt")?
            .unwrap_or_else(|| "fasta".into())
            .into_bytes(),
    })
}

fn parse_search_shapes(args: &[String]) -> Result<u32, String> {
    let mut count = 0;
    for (index, arg) in args.iter().enumerate().skip(2) {
        let attached = if arg == "--shapes" || arg == "-s" {
            None
        } else if let Some(value) = arg.strip_prefix("-s") {
            Some(value)
        } else if arg.starts_with("--shapes=") {
            return Err(format!("Invalid option: {}", &arg[2..]));
        } else {
            continue;
        };
        let mut values: Vec<&str> = attached.into_iter().collect();
        values.extend(
            args[index + 1..]
                .iter()
                .take_while(|value| !is_cpp_option_start(value))
                .map(String::as_str),
        );
        if values.len() != 1 {
            return Err("Invalid parameter count for option '-s/--shapes'".into());
        }
        count = parse_cpp_atoll(values[0]) as u32;
    }
    Ok(count)
}

fn parse_unsigned_search_option(args: &[String], name: &str, default: u32) -> Result<u32, String> {
    let mut maximum = default;
    let flag = format!("--{name}");
    let joined = format!("{flag}=");
    for (index, argument) in args.iter().enumerate().skip(2) {
        if argument.starts_with(&joined) {
            return Err(format!("Invalid option: {}", &argument[2..]));
        }
        if argument != &flag {
            continue;
        }
        let values: Vec<_> = args[index + 1..]
            .iter()
            .take_while(|value| !is_cpp_option_start(value))
            .collect();
        if values.len() != 1 {
            return Err(format!("Invalid parameter count for option '{flag}'"));
        }
        // Upstream reads unsigned options via atoll followed by a cast and
        // replaces repeated scalar occurrences.
        maximum = parse_cpp_atoll(values[0]) as u32;
    }
    Ok(maximum)
}

pub(crate) fn parse_float_search_option(
    args: &[String],
    name: &str,
    default: f64,
) -> Result<f64, String> {
    let mut value = default;
    let flag = format!("--{name}");
    let joined = format!("{flag}=");
    for (index, argument) in args.iter().enumerate().skip(2) {
        if argument.starts_with(&joined) {
            return Err(format!("Invalid option: {}", &argument[2..]));
        }
        if argument != &flag {
            continue;
        }
        let values: Vec<_> = args[index + 1..]
            .iter()
            .take_while(|value| !is_cpp_option_start(value))
            .collect();
        if values.len() != 1 {
            return Err(format!("Invalid parameter count for option '{flag}'"));
        }
        value = parse_cpp_atof(values[0]);
    }
    Ok(value)
}

fn parse_max_target_seqs(args: &[String]) -> i64 {
    match parse_arg_or(args, &["-k", "--max-target-seqs"], 25_i64) {
        // C++ normalizes the documented unlimited sentinel before search or
        // view code reaches `Config::output_range`.
        0 => i64::MAX,
        value => value,
    }
}

#[cfg(test)]
fn native_getseq_option_error(args: &[String]) -> Option<String> {
    args.iter().skip(2).find_map(|arg| {
        let error = match arg.as_str() {
            // Hidden and disabled in standard upstream builds.
            "--reverse" | "--hard-masked" | "--oid-list" | "--chunk-size" => {
                format!("Invalid option: {}", &arg[2..])
            }
            // Registered for other workflows, but not getseq.
            "--query" | "--block-size" => {
                format!("Option is not permitted for this workflow: {}", &arg[2..])
            }
            "-q" | "-b" => {
                format!("Option is not permitted for this workflow: {}", &arg[1..])
            }
            value if value.starts_with("-q") || value.starts_with("-b") => {
                format!(
                    "Option is not permitted for this workflow: {}",
                    &value[1..2]
                )
            }
            value
                if value.starts_with("--reverse=")
                    || value.starts_with("--hard-masked=")
                    || value.starts_with("--oid-list=")
                    || value.starts_with("--chunk-size=")
                    || value.starts_with("--query=")
                    || value.starts_with("--block-size=") =>
            {
                format!("Invalid option: {}", &value[2..])
            }
            _ => return None,
        };
        Some(error)
    })
}

fn parse_target_indexed(args: &[String], index_chunks: u32) -> std::io::Result<bool> {
    if !has_flag(args, "--target-indexed") {
        return Ok(false);
    }
    if index_chunks != 1 {
        print_banner();
        return Err(argument_error("--target-indexed requires -c1."));
    }
    Ok(true)
}

const DEFAULT_NATIVE_MEMORY_LIMIT: u64 = 16_000_000_000;

fn parse_memory_limit_arg(args: &[String]) -> std::io::Result<Option<u64>> {
    match get_arg(args, &["-M", "--memory-limit"]) {
        Some(value) => parse_byte_size(&value)
            .map(Some)
            .map_err(|error| argument_error(format!("invalid --memory-limit '{value}': {error}"))),
        None => Ok(Some(DEFAULT_NATIVE_MEMORY_LIMIT)),
    }
}

fn parse_byte_size(value: &str) -> Result<u64, &'static str> {
    let bytes = crate::util::string::interpret_number(value)
        .map_err(|_| "use K, M, G, or T (for example 4G)")?;
    u64::try_from(bytes).map_err(|_| "value is out of range")
}

fn parse_sensitivity_checked(args: &[String]) -> Result<Sensitivity, String> {
    let mut selected = None;
    for (flag, level) in [
        ("--faster", Sensitivity::Faster),
        ("--fast", Sensitivity::Fast),
        ("--mid-sensitive", Sensitivity::MidSensitive),
        ("--sensitive", Sensitivity::Sensitive),
        ("--more-sensitive", Sensitivity::MoreSensitive),
        ("--very-sensitive", Sensitivity::VerySensitive),
        ("--ultra-sensitive", Sensitivity::UltraSensitive),
        ("--shapes-30x10", Sensitivity::Shapes30x10),
        ("--shapes-6x10", Sensitivity::Shapes6x10),
        ("--linclust-20", Sensitivity::Linclust20),
    ] {
        for (index, argument) in args.iter().enumerate().skip(2) {
            if argument == flag
                && args
                    .get(index + 1)
                    .is_some_and(|value| !is_cpp_option_start(value))
            {
                return Err(format!("Invalid parameter count for option '{flag}'"));
            }
        }
        if has_flag(args, flag) {
            if selected.is_some() {
                return Err("Sensitivity switches are mutually exclusive.".into());
            }
            selected = Some(level);
        }
    }
    Ok(selected.unwrap_or(Sensitivity::Default))
}

fn parse_sensitivity(args: &[String]) -> std::io::Result<Sensitivity> {
    parse_sensitivity_checked(args).map_err(argument_error)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_all_args_collects_outfmt_field_list() {
        let args = vec![
            "diamond".to_string(),
            "view".to_string(),
            "-f".to_string(),
            "6".to_string(),
            "qseqid".to_string(),
            "sseqid".to_string(),
            "--top".to_string(),
            "10".to_string(),
            "--outfmt".to_string(),
            "104".to_string(),
            "score".to_string(),
        ];

        assert_eq!(
            get_all_args(&args, &["-f", "--outfmt"]),
            ["6", "qseqid", "sseqid", "104", "score"]
        );
    }

    #[test]
    fn native_getseq_rejects_ignored_options_like_standard_upstream() {
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "--hard-masked", "--reverse"])),
            Some("Invalid option: hard-masked".to_string())
        );
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "--query", "ids.tsv"])),
            Some("Option is not permitted for this workflow: query".to_string())
        );
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "-b2"])),
            Some("Option is not permitted for this workflow: b".to_string())
        );
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "--oid-list=ids.tsv"])),
            Some("Invalid option: oid-list=ids.tsv".to_string())
        );
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "--seq", "1"])),
            None
        );
    }

    #[test]
    fn test_get_all_args_accepts_joined_first_value() {
        assert_eq!(
            get_all_args(
                &args(&["diamond", "blastp", "--outfmt=6", "qseqid", "sseqid"]),
                &["-f", "--outfmt"]
            ),
            ["6", "qseqid", "sseqid"]
        );
        assert_eq!(
            get_all_args(
                &args(&["diamond", "blastp", "-f6", "qseqid"]),
                &["-f", "--outfmt"]
            ),
            ["6", "qseqid"]
        );
    }

    #[test]
    fn last_arg_group_replaces_repeated_occurrences() {
        assert_eq!(
            get_last_arg_group(
                &args(&[
                    "diamond",
                    "makedb",
                    "--in",
                    "first.faa",
                    "--in",
                    "second.faa"
                ]),
                &["--in"]
            ),
            ["second.faa"]
        );
        assert_eq!(
            get_last_arg_group(
                &args(&["diamond", "makedb", "--in", "one.faa", "two.faa"]),
                &["--in"]
            ),
            ["one.faa", "two.faa"]
        );
        assert_eq!(
            get_last_arg_group(
                &args(&["diamond", "makedb", "--in=first.faa", "--in=last.faa"]),
                &["--in"]
            ),
            ["last.faa"]
        );
    }

    #[test]
    fn checked_last_arg_group_rejects_missing_values() {
        assert_eq!(
            get_last_arg_group_checked(
                &args(&["diamond", "cluster", "--cluster-steps"]),
                &["--cluster-steps"]
            ),
            Err("Option --cluster-steps requires at least one value.".to_string())
        );
        assert_eq!(
            get_last_arg_group_checked(
                &args(&[
                    "diamond",
                    "cluster",
                    "--cluster-steps",
                    "fast",
                    "--cluster-steps",
                    "default"
                ]),
                &["--cluster-steps"]
            ),
            Ok(vec!["default".to_string()])
        );
    }

    #[test]
    fn native_cluster_checks_external_workflow_requirements() {
        assert_eq!(
            validate_native_cluster_options(&args(&[
                "diamond",
                "cluster",
                "--aln-out",
                "edges.tsv"
            ])),
            Err("Alignment output to file is not supported without --parallel-tmpdir.".to_string())
        );
        assert_eq!(
            validate_native_cluster_options(&args(&[
                "diamond",
                "cluster",
                "--parallel-tmpdir=/tmp/work"
            ])),
            Err("Option is not permitted for this workflow: --parallel-tmpdir".to_string())
        );
        assert_eq!(
            validate_native_cluster_options(&args(&[
                "diamond",
                "linclust",
                "--parallel-tmpdir",
                "/tmp/work"
            ])),
            Ok(())
        );
    }

    #[test]
    fn native_cluster_parser_forwards_vector_threshold_and_memory_controls() {
        let mut arguments = args(&[
            "diamond",
            "cluster",
            "-d",
            "input.dmnd",
            "-o",
            "clusters.tsv",
            "--member-cover",
            "75",
            "--mutual-cover",
            "80",
            "--cluster-steps",
            "fast",
            "default",
            "--cluster-steps",
            "faster_lin",
            "sensitive",
            "--connected-component-depth",
            "2",
            "--round-coverage",
            "81",
            "82",
            "--round-approx-id",
            "55",
            "60",
            "--kmer-ranking",
            "--no-reassign",
            "--matrix",
            "pam250",
            "--gapopen",
            "15",
            "--gapextend",
            "3",
            "--evalue",
            "1e-20",
            "--comp-based-stats",
            "4",
            "--masking",
            "tantan",
            "--soft-masking",
            "seg",
            "--motif-masking",
            "0",
            "--ext",
            "banded-slow",
            "--oid-output",
            "--header",
            "simple",
            "-M3G",
        ]);
        if cfg!(feature = "extra") {
            arguments.extend(args(&[
                "--weighted-gvc",
                "--strict-gvc",
                "--graph-algo",
                "len",
            ]));
        }
        let parsed = parse_native_cluster_config(&arguments).unwrap();
        assert_eq!(parsed.member_cover, Some(75.0));
        assert_eq!(parsed.mutual_cover, Some(80.0));
        assert_eq!(
            parsed.cluster_steps,
            vec![b"faster_lin".to_vec(), b"sensitive".to_vec()]
        );
        assert_eq!(parsed.connected_component_depth, vec![b"2".to_vec()]);
        assert_eq!(parsed.round_coverage, vec![b"81".to_vec(), b"82".to_vec()]);
        assert_eq!(parsed.round_approx_id, vec![b"55".to_vec(), b"60".to_vec()]);
        assert!(parsed.kmer_ranking && parsed.no_reassign);
        #[cfg(feature = "extra")]
        {
            assert!(parsed.weighted_gvc && parsed.strict_gvc);
            assert_eq!(parsed.graph_algo, GraphAlgo::LengthSorted);
        }
        #[cfg(not(feature = "extra"))]
        {
            assert!(!parsed.weighted_gvc && !parsed.strict_gvc);
            assert_eq!(parsed.graph_algo, GraphAlgo::GreedyVertexCover);
        }
        assert_eq!(parsed.memory_limit, 3_000_000_000);
        assert_eq!(parsed.matrix, b"pam250");
        assert_eq!(parsed.gap_open, 15);
        assert_eq!(parsed.gap_extend, 3);
        assert_eq!(parsed.max_evalue, 1e-20);
        assert_eq!(parsed.comp_based_stats, 4);
        assert_eq!(parsed.masking, crate::masking::MaskingMode::Tantan);
        assert_eq!(parsed.soft_masking, crate::masking::MaskingAlgo::Seg);
        assert_eq!(parsed.motif_masking, b"0");
        assert_eq!(parsed.extension_mode, b"banded-slow");
        assert!(parsed.oid_output);
        assert!(parsed.simple_header);
    }

    #[test]
    fn native_cluster_search_control_defaults_match_init_thresholds() {
        let parsed = parse_native_cluster_config(&args(&["diamond", "cluster"])).unwrap();
        assert_eq!(parsed.max_evalue, 0.001);
        assert_eq!(parsed.comp_based_stats, 1);
        assert_eq!(parsed.masking, crate::masking::MaskingMode::None);
        assert_eq!(parsed.soft_masking, crate::masking::MaskingAlgo::Tantan);
        assert!(parsed.motif_masking.is_empty());
        assert!(parsed.extension_mode.is_empty());
        assert!(parse_native_cluster_config(&args(&[
            "diamond",
            "cluster",
            "--comp-based-stats=6"
        ]))
        .is_err());
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--comp-based-stats=5"]))
                .unwrap()
                .comp_based_stats,
            5
        );
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--comp-based-stats=7"]))
                .unwrap_err(),
            "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5."
        );
        assert!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--ext=invalid"])).is_err()
        );
    }

    #[test]
    fn native_cluster_header_matches_upstream_forms_and_errors() {
        assert!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--header"]))
                .unwrap()
                .simple_header
        );
        assert!(
            !parse_native_cluster_config(&args(&["diamond", "cluster", "--header", "0"]))
                .unwrap()
                .simple_header
        );
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--header", "verbose"]))
                .unwrap_err(),
            "Verbose header format is not supported for cluster workflow."
        );
    }

    #[test]
    fn native_cluster_scoring_controls_require_values() {
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--matrix"])).unwrap_err(),
            "Option --matrix requires a value."
        );
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--gapopen", "bad"]))
                .unwrap()
                .gap_open,
            0
        );
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--gapextend"])).unwrap_err(),
            "Option --gapextend requires a value."
        );
    }

    #[test]
    fn native_search_output_parses_daa_aliases_and_suffixes() {
        assert_eq!(
            parse_native_search_output(&args(&["diamond", "blastp", "-a", "hits"])),
            Ok((Some("hits.daa".to_string()), vec!["daa".to_string()]))
        );
        assert_eq!(
            parse_native_search_output(&args(&[
                "diamond", "blastx", "--outfmt", "100", "--out", "hits"
            ])),
            Ok((Some("hits.daa".to_string()), vec!["100".to_string()]))
        );
        assert_eq!(
            parse_native_search_output(&args(&[
                "diamond",
                "blastp",
                "--daa=hits.daa",
                "--outfmt",
                "daa"
            ])),
            Ok((Some("hits.daa".to_string()), vec!["daa".to_string()]))
        );
    }

    #[test]
    fn native_search_output_rejects_conflicting_or_invalid_daa_options() {
        assert!(parse_native_search_output(&args(&[
            "diamond", "blastp", "-a", "hits", "-o", "other"
        ]))
        .is_err());
        assert!(parse_native_search_output(&args(&[
            "diamond", "blastp", "-a", "hits", "--outfmt", "6"
        ]))
        .is_err());
        assert!(
            parse_native_search_output(&args(&["diamond", "blastp", "--outfmt", "100"])).is_err()
        );
    }

    #[test]
    fn parse_memory_sizes() {
        assert_eq!(parse_byte_size("512M"), Ok(512_000_000));
        assert_eq!(parse_byte_size("1.5G"), Ok(1_500_000_000));
        assert!(parse_byte_size("lots").is_err());
        assert!(parse_byte_size("512").is_err());

        let args = vec!["diamond".to_string(), "blastp".to_string()];
        assert_eq!(parse_memory_limit_arg(&args).unwrap(), Some(16_000_000_000));
        let args = vec![
            "diamond".to_string(),
            "blastp".to_string(),
            "--memory-limit".to_string(),
            "0G".to_string(),
        ];
        assert_eq!(parse_memory_limit_arg(&args).unwrap(), Some(0));
    }

    #[test]
    fn zero_max_target_sequences_means_unlimited() {
        assert_eq!(parse_max_target_seqs(&args(&["diamond", "view"])), 25);
        assert_eq!(
            parse_max_target_seqs(&args(&["diamond", "view", "-k", "0"])),
            i64::MAX
        );
        assert_eq!(
            parse_max_target_seqs(&args(&["diamond", "blastp", "--max-target-seqs=7"])),
            7
        );
    }

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
    }

    #[test]
    fn native_view_rejects_search_only_options_like_upstream() {
        #[cfg(not(feature = "extra"))]
        {
            assert_eq!(
                parse_native_view_cli(&args(&["diamond", "view", "--unal", "0"])).unwrap_err(),
                "Option is not permitted for this workflow: unal"
            );
            assert_eq!(
                parse_native_view_cli(&args(&["diamond", "view", "--compress", "1"])).unwrap_err(),
                "Option is not permitted for this workflow: compress"
            );
        }
        #[cfg(feature = "extra")]
        {
            assert!(
                parse_native_view_cli(&args(&["diamond", "view", "--unal", "0"]))
                    .unwrap()
                    .runtime
                    .source_extra
            );
            assert_eq!(
                parse_native_view_cli(&args(&["diamond", "view", "--compress", "1"]))
                    .unwrap()
                    .runtime
                    .compression,
                b"1"
            );
            assert_eq!(
                parse_native_view_cli(&args(&["diamond", "view", "--compress"])).unwrap_err(),
                "Invalid parameter count for option '--compress'"
            );
        }
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--wat"])).unwrap_err(),
            "Invalid option: wat"
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--top=10"])).unwrap_err(),
            "Invalid option: top=10"
        );
    }

    #[test]
    fn native_view_offsets_match_upstream_atoll_and_unsigned_conversion() {
        assert_eq!(parse_cpp_atoll("nope"), 0);
        assert_eq!(parse_cpp_atoll("  +12tail"), 12);
        assert_eq!(parse_cpp_atoll("9223372036854775808"), i64::MAX);
        assert_eq!(parse_cpp_atoll("-9223372036854775809"), i64::MIN);
        assert_eq!(
            parse_native_view_cli(&args(&[
                "diamond",
                "view",
                "--qnum-offset",
                "-1",
                "--snum-offset",
                "20",
            ]))
            .map(|config| (config.qnum_offset, config.snum_offset)),
            Ok((u64::MAX, 20))
        );
        // Like the C++ option parser, the last repeated occurrence wins.
        assert_eq!(
            parse_native_view_cli(&args(&[
                "diamond",
                "view",
                "--qnum-offset",
                "3",
                "--qnum-offset",
                "-2",
            ]))
            .map(|config| (config.qnum_offset, config.snum_offset)),
            Ok((u64::MAX - 1, 0))
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--qnum-offset"])).map(|_| ()),
            Err("Invalid parameter count for option '--qnum-offset'".to_string())
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--qnum-offset=-1",])).map(|_| ()),
            Err("Invalid option: qnum-offset=-1".to_string())
        );
    }

    #[test]
    fn native_view_top_and_vector_options_follow_cpp_rules() {
        let top_zero =
            parse_native_view_cli(&args(&["diamond", "view", "--top", "nonsense"])).unwrap();
        assert_eq!(top_zero.toppercent, Some(0.0));

        let unlimited = parse_native_view_cli(&args(&["diamond", "view", "--top", "100"])).unwrap();
        assert_eq!(unlimited.toppercent, None);
        assert_eq!(unlimited.max_target_seqs, i64::MAX);

        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--top", "-1"])).unwrap_err(),
            "Allowed value range for --top is between 0.0 and 100.0"
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--top", "10", "-k5"])).unwrap_err(),
            "--top and -k/--max-target-seqs are mutually exclusive."
        );

        let repeated = parse_native_view_cli(&args(&[
            "diamond", "view", "-f", "6", "qseqid", "-f104", "score", "--header",
        ]))
        .unwrap();
        assert_eq!(repeated.outfmt, vec![b"104".to_vec(), b"score".to_vec()]);
        assert_eq!(repeated.header, vec![Vec::<u8>::new()]);

        // Source view.cpp opens DAA and initializes its score matrix before
        // Search::Config constructs/validates the output format. The parser
        // must retain these controls for that actual owner boundary.
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "-f", "wat"]))
                .unwrap()
                .outfmt,
            vec![b"wat".to_vec()]
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "-f102"]))
                .unwrap()
                .outfmt,
            vec![b"102".to_vec()]
        );
    }

    #[test]
    fn native_source_sensitivity_flags_are_complete_and_exclusive() {
        for (flag, level) in [
            ("--faster", Sensitivity::Faster),
            ("--fast", Sensitivity::Fast),
            ("--mid-sensitive", Sensitivity::MidSensitive),
            ("--sensitive", Sensitivity::Sensitive),
            ("--more-sensitive", Sensitivity::MoreSensitive),
            ("--very-sensitive", Sensitivity::VerySensitive),
            ("--ultra-sensitive", Sensitivity::UltraSensitive),
            ("--shapes-30x10", Sensitivity::Shapes30x10),
            ("--shapes-6x10", Sensitivity::Shapes6x10),
            ("--linclust-20", Sensitivity::Linclust20),
        ] {
            for mode in ["blastp", "blastx", "makeidx"] {
                let command = args(&["diamond", mode, flag]);
                assert_eq!(parse_sensitivity_checked(&command), Ok(level));
                if mode != "makeidx" {
                    assert!(!has_unknown_search_option(&command, mode == "blastx"));
                }
                assert_eq!(
                    parse_sensitivity_checked(&args(&["diamond", mode, flag, flag])),
                    Ok(level)
                );
                assert!(
                    parse_sensitivity_checked(&args(&["diamond", mode, flag, "1"]))
                        .unwrap_err()
                        .contains("Invalid parameter count")
                );
                if flag != "--sensitive" {
                    assert_eq!(
                        parse_sensitivity_checked(&args(&["diamond", mode, flag, "--sensitive"]))
                            .unwrap_err(),
                        "Sensitivity switches are mutually exclusive."
                    );
                }
            }
        }
    }

    #[test]
    fn native_query_indexed_algorithm_parser_and_route() {
        for options in [
            vec!["--algo", "0"],
            vec!["--algo", "1"],
            vec!["--algo", ""],
            vec!["--algo", "0", "--algo", "1"],
        ] {
            for mode in ["blastp", "blastx"] {
                let mut command = args(&["diamond", mode]);
                command.extend(options.iter().map(|value| value.to_string()));
                assert_eq!(
                    parse_search_algorithm(&command).unwrap().as_deref(),
                    options.last().copied()
                );
                assert!(!route_common_search_to_legacy(&command, mode == "blastx"));
            }
        }
        for options in [vec!["--algo"], vec!["--algo", "1", "2"]] {
            let mut command = args(&["diamond", "blastp"]);
            command.extend(options.iter().map(|value| value.to_string()));
            assert_eq!(
                parse_search_algorithm(&command).unwrap_err(),
                "Invalid parameter count for option '--algo'"
            );
        }
        assert_eq!(
            parse_search_algorithm(&args(&["diamond", "blastp", "--algo=1"])).unwrap_err(),
            "Invalid option: algo=1"
        );
        assert!(!route_common_search_to_legacy(
            &args(&["diamond", "blastp", "--algo", "ctg"]),
            false
        ));
    }

    #[test]
    fn native_index_chunks_preserves_source_scalar_parsing() {
        for (options, expected) in [
            (vec![], 0),
            (vec!["-c1"], 1),
            (vec!["-c", "2"], 2),
            (vec!["--index-chunks", "3"], 3),
            (vec!["-c1", "--index-chunks", "4"], 4),
            (vec!["-c0"], 0),
            (vec!["-cno-number"], 0),
            (vec!["-c-1"], u32::MAX),
        ] {
            for mode in ["blastp", "blastx"] {
                let mut command = args(&["diamond", mode]);
                command.extend(options.iter().map(|value| value.to_string()));
                assert_eq!(parse_index_chunks(&command), Ok(expected));
                assert!(!route_common_search_to_legacy(&command, mode == "blastx"));
            }
        }
        for options in [
            vec!["-c"],
            vec!["-c1", "2"],
            vec!["--index-chunks", "1", "2"],
        ] {
            let mut command = args(&["diamond", "blastp"]);
            command.extend(options.iter().map(|value| value.to_string()));
            assert_eq!(
                parse_index_chunks(&command).unwrap_err(),
                "Invalid parameter count for option '-c/--index-chunks'"
            );
        }
        assert_eq!(
            parse_index_chunks(&args(&["diamond", "blastp", "--index-chunks=2"])).unwrap_err(),
            "Invalid option: index-chunks=2"
        );
    }

    #[test]
    fn native_block_routes_preserve_parser_values_and_supported_modes() {
        for mode in ["blastp", "blastx"] {
            for (options, expected) in [
                (vec!["-b0.0000001"], Some(0.0000001)),
                (vec!["--block-size", "0"], Some(0.0)),
                (vec!["-b1", "-b", "2tail"], Some(2.0)),
                (vec![], None),
            ] {
                let mut command = args(&["diamond", mode]);
                command.extend(options.into_iter().map(str::to_owned));
                assert_eq!(parse_block_size(&command), Ok(expected));
                assert!(!route_common_search_to_legacy(&command, mode == "blastx"));
            }
            let command = args(&["diamond", mode, "-b1", "-f", "6", "staxids"]);
            assert!(!route_common_search_to_legacy(&command, mode == "blastx"));
        }
        for options in [vec!["-b"], vec!["-b1", "2"], vec!["--block-size", "1", "2"]] {
            let mut command = args(&["diamond", "blastp"]);
            command.extend(options.into_iter().map(str::to_owned));
            assert_eq!(
                parse_block_size(&command).unwrap_err(),
                "Invalid parameter count for option '-b/--block-size'"
            );
        }
        assert_eq!(
            parse_block_size(&args(&["diamond", "blastp", "--block-size=1"])).unwrap_err(),
            "Invalid option: block-size=1"
        );
    }

    #[test]
    fn native_route_never_silently_ignores_unsupported_options() {
        for mode in ["blastp", "blastx"] {
            for option in [
                ["--soft-masking", "0"],
                ["--soft-masking", "seg"],
                ["--soft-masking", "tantan"],
                ["--approx-id", "80"],
            ] {
                let mut command = args(&["diamond", mode]);
                command.extend(option.into_iter().map(str::to_owned));
                assert!(!route_common_search_to_legacy(&command, mode == "blastx"));
            }
        }

        for mode in ["blastp", "blastx"] {
            for format in ["100", "6", "104", "sam", "null"] {
                for flag in ["--salltitles", "--sallseqid"] {
                    let command = args(&["diamond", mode, "-f", format, flag]);
                    assert!(!has_unknown_search_option(&command, mode == "blastx"));
                    assert!(!route_common_search_to_legacy(&command, mode == "blastx"));
                }
            }
            for value in ["0", "1"] {
                let command = args(&["diamond", mode, "-f", "100", "--unal", value]);
                assert!(!route_common_search_to_legacy(&command, mode == "blastx"));
            }
        }

        for format in [
            "0", "5", "xml", "103", "paf", "101", "sam", "null", "edge", "102",
        ] {
            assert!(!route_blastp_to_legacy(&args(&[
                "diamond", "blastp", "--outfmt", format
            ])));
            assert!(!route_blastx_to_legacy(&args(&[
                "diamond", "blastx", "--outfmt", format
            ])));
        }
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond", "blastp", "--outfmt", "6", "qseqid", "sseqid"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt=6",
            "qlen",
            "nident",
            "score"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond", "blastp", "--outfmt", "6", "qseq"
        ])));
        // Native DMND databases carry the same mapping, nodes, ranks, and
        // scientific-name payloads as upstream DatabaseFile.
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt",
            "6",
            "staxids",
            "sscinames",
            "sskingdoms",
            "skingdoms",
            "sphylums",
            "slineages"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt",
            "6",
            "full_qseq_mate"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt",
            "6",
            "not_a_field"
        ])));
        // Pinned DIAMOND does not register the BLAST accession aliases in its
        // tabular field map and rejects them as invalid output fields.
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond", "blastp", "--outfmt", "6", "qacc"
        ])));
        assert!(route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--definitely-not-a-diamond-option"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "-qquery.faa"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--target-indexed",
            "-c1"
        ])));
    }

    #[test]
    fn blastx_specific_unsupported_modes_route_to_legacy() {
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--frameshift",
            "15"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond", "blastx", "--swipe"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--min-query-len",
            "30"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--query-gencode",
            "2",
            "--strand",
            "plus",
            "--frameshift",
            "0",
            "--outfmt",
            "6"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond", "blastx", "--outfmt", "6", "qseqid", "qlen", "sseqid"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--outfmt",
            "6",
            "qseqid",
            "full_qseq_mate",
            "sseqid"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--no-self-hits"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "-q",
            "query.fna",
            "-d",
            "db.faa",
            "--no-self-hits",
            "-o",
            "out.tsv"
        ])));
    }
}
