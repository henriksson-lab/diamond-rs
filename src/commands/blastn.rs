//! Optional nucleotide search with native minimizer, KSW/WFA and output engines.
use crate::dna::native::DnaOptions;
use crate::util::command_line_parser::ReadOption;
use std::ffi::OsString;
use std::io;

#[derive(Clone)]
pub struct BlastnConfig {
    pub search: super::blastp::BlastpConfig,
    pub dna: DnaOptions,
}

pub(crate) fn run_cli(args: &[String]) -> io::Result<()> {
    run_cli_with_startup(args, None)
}

pub(crate) fn run_cli_with_startup(
    args: &[String],
    startup: Option<&super::source_config::Startup>,
) -> io::Result<()> {
    let args = args.iter().map(OsString::from).collect::<Vec<_>>();
    run_with_startup_os(&args, startup)
}

pub(crate) fn run_with_startup_os(
    args: &[OsString],
    startup: Option<&super::source_config::Startup>,
) -> io::Result<()> {
    let owned_startup;
    let startup = match startup {
        Some(startup) => startup,
        None => {
            owned_startup = super::source_config::prepare_os(args)?;
            owned_startup.as_ref()
        }
    };
    let values = super::argument_registry::extra_values_os(args)?;
    let mut search = super::dispatch::parse_native_blastp_config_os(
        args,
        &crate::basic::config::invocation_bytes_io(args)?,
    )?;
    startup.apply_search(&mut search);
    let mut dna = DnaOptions::default();
    macro_rules! scalar {
        ($name:literal,$field:ident,$ty:ty) => {
            if let Some(items) = values.get($name) {
                <$ty>::read_option_os(&mut dna.$field, items)?;
            }
        };
    }
    scalar!("zdrop", zdrop, i32);
    scalar!("repetition-cutoff", repetitive_cutoff, f64);
    scalar!("chain-pen-gap-scale", chain_pen_gap_scale, f64);
    scalar!("chain-pen-skip-scale", chain_pen_skip_scale, f64);
    scalar!("penalty", penalty, i32);
    scalar!("reward", reward, i32);
    scalar!("zdrop-extension", zdrop_extension, i32);
    scalar!("zdrop-global", zdrop_global, i32);
    scalar!("band-extension", band_extension, i32);
    scalar!("band-global", band_global, i32);
    if let Some(items) = values.get("extension") {
        let mut name = Vec::<u8>::new();
        Vec::<u8>::read_option_os(&mut name, items)?;
        dna.algorithm = match name.as_slice() {
            b"ksw" => crate::dna::extension::DnaExtensionAlgo::Ksw,
            b"wfa" => crate::dna::extension::DnaExtensionAlgo::Wfa,
            _ => {
                let mut message = b"Invalid DNA extension algorithm: ".to_vec();
                message.extend_from_slice(&name);
                return Err(crate::util::io::RawDiagnosticError(message).into_io());
            }
        };
    }
    macro_rules! optional {
        ($name:literal,$field:ident,$ty:ty) => {
            if let Some(items) = values.get($name) {
                let mut value = <$ty>::default();
                <$ty>::read_option_os(&mut value, items)?;
                dna.$field = Some(value);
            }
        };
    }
    optional!("chain-align-cutoff", chain_align_cutoff, f64);
    optional!("min-chain-score", min_chain_score, i32);
    optional!("max-overlap-extension", max_overlap_extension, f64);
    dna.chaining_out = values.contains_key("chaining-out");
    dna.align_long_reads = values.contains_key("align-long-reads");
    dna.best_hsp_only = values.contains_key("best-hsp-only");
    run(&BlastnConfig { search, dna })
}

pub fn run(config: &BlastnConfig) -> io::Result<()> {
    if config.search.database.as_os_str().is_empty() {
        return Err(io::Error::other("Missing required option: --db/-d"));
    }
    let size = config.search.block_size.unwrap_or(1.0);
    if !size.is_finite() || size <= 0.0 {
        return Err(io::Error::other("Invalid block size"));
    }
    let budget = (size * 1e9) as u64;
    super::search_blocks::run_nucleotide_files(config, budget.max(1))
}
