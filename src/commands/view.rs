use std::collections::BTreeMap;
#[cfg(test)]
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::{mpsc, Arc, Mutex};

use crate::basic::consts::VERSION_STRING;
use crate::data::daa::{
    finish_daa_from_file, init_daa, view_query_daa, view_query_edge, view_query_null,
    view_query_paf, view_query_pairwise, view_query_sam, ByteDaaFile as DaaFile,
    ByteDaaQueryRecord as DaaQueryRecord, DaaViewConfig,
};
use crate::dp::swipe::HspValues;
use crate::output::format::{
    print_tabular_footer, print_tabular_header_bytes, Header, TabularFormat,
};
use crate::output::{pairwise, sam, xml};
use crate::stats::score_matrix::ScoreMatrix;

macro_rules! view_status {
    ($verbosity:expr, $($argument:tt)*) => {
        crate::util::log_stream::MessageStream::new($verbosity > 0, $verbosity == 3)
            .write(format!("{}\n", format_args!($($argument)*)))?;
    };
}

#[derive(Debug, Clone)]
pub struct ViewConfig {
    pub daa_file: PathBuf,
    pub output: PathBuf,
    pub outfmt: Vec<Vec<u8>>,
    pub max_target_seqs: i64,
    pub toppercent: Option<f64>,
    pub forward_only: bool,
    pub report_unaligned: bool,
    pub sam_qlen_field: bool,
    pub invocation: Vec<u8>,
    /// --header value (`simple` / `verbose` / unset).
    pub header: Vec<Vec<u8>>,
    /// Worker count (`0` selects the host's available parallelism, as upstream).
    pub threads: usize,
    /// Unsigned representation of upstream's signed int64 ordinal offsets.
    pub qnum_offset: u64,
    pub snum_offset: u64,
    /// Upstream console verbosity: 0=quiet, 1=normal, 2=verbose, 3=debug log.
    pub verbosity: u8,
    pub file_buffer_size: usize,
}
impl Default for ViewConfig {
    fn default() -> Self {
        Self {
            daa_file: PathBuf::new(),
            output: PathBuf::new(),
            outfmt: Vec::new(),
            max_target_seqs: 0,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: Vec::new(),
            header: Vec::new(),
            threads: 0,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
            file_buffer_size: 67_108_864,
        }
    }
}

/// Actual Source view controls admitted by EXTRA, separate from the existing
/// typed library configuration. Compression is resolved at writer construction
/// after archive/model/output initialization, preserving Source error order.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct RuntimeOptions {
    pub compression: Vec<u8>,
    pub source_extra: bool,
    pub xml_blord_format: bool,
    pub no_parse_seqids: bool,
    pub database: PathBuf,
    pub taxonomy_names: PathBuf,
    pub frame_shift: i32,
    pub extension: Vec<u8>,
    pub multiprocessing: bool,
    pub global_ranking: i64,
    pub query_range_culling: bool,
    pub search_input: crate::run::config::ConfigInput,
    pub masking: Vec<u8>,
    pub include_lineage: bool,
    pub min_bit_score: f64,
    pub top_present: bool,
}

pub(crate) fn run_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    match config.outfmt.first().map(Vec::as_slice).unwrap_or(b"6") {
        b"100" | b"daa" => run_daa_with_options(config, runtime),
        b"6" | b"tab" | b"104" | b"json-flat" => run_tabular_with_options(config, runtime),
        b"103" | b"paf" => run_paf_with_options(config, runtime),
        b"101" | b"sam" => run_sam_with_options(config, runtime),
        b"0" => run_pairwise_with_options(config, runtime),
        b"5" | b"xml" => run_xml_with_options(config, runtime),
        b"null" => run_null_with_options(config, runtime),
        b"edge" => run_edge_with_options(config, runtime),
        value => Err(raw_error(b"Invalid output format: ", value, b"")),
    }
}

fn raw_error(prefix: &[u8], value: &[u8], suffix: &[u8]) -> io::Error {
    let mut message = prefix.to_vec();
    message.extend_from_slice(value);
    message.extend_from_slice(suffix);
    io::Error::other(crate::util::io::RawDiagnosticError(message))
}

fn prepare_source_output(
    config: &ViewConfig,
    runtime: &RuntimeOptions,
    daa: &DaaFile,
) -> io::Result<Option<crate::output::format::OutputInit>> {
    if !runtime.source_extra {
        return Ok(None);
    }
    let mut input = runtime.search_input.clone();
    input.query_translated = daa.mode() == 3;
    input.frame_shift = runtime.frame_shift;
    input.global_ranking_targets = runtime.global_ranking;
    input.multiprocessing = runtime.multiprocessing;
    input.masking = match runtime.masking.to_ascii_lowercase().as_slice() {
        b"" | b"1" | b"tantan" => crate::masking::MaskingMode::Tantan,
        b"0" | b"none" => crate::masking::MaskingMode::None,
        b"seg" => crate::masking::MaskingMode::BlastSeg,
        _ => {
            return Err(raw_error(
                b"Invalid value for string field: ",
                &runtime.masking,
                b". Permitted values: none, tantan, seg",
            ))
        }
    };
    let fields = config.outfmt.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let archive_matrix = daa.score_matrix_for_alignment()?;
    let errors = crate::util::io::IoErrorOwner::default();
    let source = errors.finish(crate::run::config::Config::new_with_io_owner(
        &input,
        |maximum| {
            let output = errors.capture(crate::output::format::init_output_bytes(
                &crate::output::format::OutputInitConfigBytes {
                    output_format: &fields,
                    daa_file: &config.daa_file,
                    command: crate::output::format::OutputCommand::View,
                    include_lineage: runtime.include_lineage,
                    frame_shift: runtime.frame_shift,
                    query_translated: daa.mode() == 3,
                    score_matrix_name: archive_matrix.as_bytes(),
                    ext: &runtime.extension,
                    multiprocessing: runtime.multiprocessing,
                    global_ranking_targets: u32::from(runtime.global_ranking != 0),
                    toppercent: config.toppercent,
                    max_target_seqs: Some(config.max_target_seqs),
                    min_bit_score: runtime.min_bit_score,
                    query_range_culling: runtime.query_range_culling,
                },
            ))?;
            *maximum = output.max_target_seqs;
            Ok(output)
        },
        &errors,
    ))?;
    if let Some(message) = source.iteration_message {
        view_status!(config.verbosity, "{}", message);
    }
    if config.verbosity >= 3 {
        view_status!(
            config.verbosity,
            "Min length ratio: {}",
            crate::output::format::format_default_float(source.min_length_ratio)
        );
    }
    Ok(Some(source.output_format))
}

/// Source calls taxonomy.init() after Search::Config and before ViewWriter.
/// Keep its actual names owner resident through rendering and output closure.
fn load_view_taxonomy(
    config: &ViewConfig,
    runtime: &RuntimeOptions,
) -> io::Result<crate::data::taxonomy::ByteTaxonomy> {
    let mut taxonomy = crate::data::taxonomy::ByteTaxonomy::default();
    if runtime.source_extra && !runtime.taxonomy_names.is_empty() {
        let mut timer = crate::util::log_stream::TaskTimer::with_msg("Loading taxonomy names", 1);
        let count =
            taxonomy.init_source_bytes_io_path(&runtime.taxonomy_names, config.file_buffer_size)?;
        timer.finish();
        view_status!(
            config.verbosity,
            "Loaded taxonomy names for {} taxon ids.",
            count
        );
    }
    Ok(taxonomy)
}

fn print_daa_metadata(daa: &DaaFile, score_matrix: &ScoreMatrix, verbosity: u8) -> io::Result<()> {
    if verbosity == 0 {
        return Ok(());
    }
    // Upstream constructs `ScoreMatrix` from the archive's matrix name, gap
    // penalties, and database size, then prints that derived object.  It does
    // not echo the archived K/Lambda header fields directly.
    view_status!(
        verbosity,
        "Scoring parameters: (Matrix={} Lambda={} K={} Penalties={}/{})",
        daa.score_matrix_for_alignment()?,
        crate::output::format::format_default_float(score_matrix.lambda()),
        crate::output::format::format_default_float(score_matrix.k()),
        score_matrix.gap_open(),
        score_matrix.gap_extend()
    );
    if verbosity >= 2 {
        view_status!(verbosity, "Build version = {}", daa.diamond_build());
    }
    view_status!(verbosity, "DB sequences = {}", daa.db_seqs());
    view_status!(verbosity, "DB sequences used = {}", daa.db_seqs_used());
    view_status!(verbosity, "DB letters = {}", daa.db_letters());
    Ok(())
}

/// Print the messages produced when upstream constructs `Search::Config`
/// inside `view_daa`.  This happens after the archive metadata has been
/// printed, so keeping it here (rather than in the CLI wrapper) also preserves
/// the observable ordering of verbose/debug output.
fn print_view_selection(
    config: &ViewConfig,
    dp_fields: HspValues,
    runtime: &RuntimeOptions,
) -> io::Result<()> {
    if config.verbosity == 0 {
        return Ok(());
    }
    if config.verbosity >= 3 && !runtime.source_extra {
        view_status!(config.verbosity, "Min length ratio: 0");
    }
    if let Some(top) = config.toppercent {
        view_status!(
            config.verbosity,
            "Percentage range of top alignment score to report hits: {}",
            crate::output::format::format_default_float(top)
        );
    } else if config.max_target_seqs == i64::MAX {
        view_status!(
            config.verbosity,
            "#Target sequences to report alignments for: unlimited"
        );
    } else {
        view_status!(
            config.verbosity,
            "#Target sequences to report alignments for: {}",
            config.max_target_seqs
        );
    }
    if config.verbosity >= 3 {
        view_status!(config.verbosity, "DP fields: {}", dp_fields.0);
    }
    Ok(())
}

fn view_threads(requested: usize) -> usize {
    if requested == 0 {
        crate::util::system::hardware_concurrency().max(1) as usize
    } else {
        requested
    }
}

/// Parse the `--header` flag value (as a Vec<String> coming from CLI parsing).
/// Mirrors `TabularFormat::header_format(Config::blastp)`: missing flag means
/// no header, bare `--header` means verbose, and an explicit value must be one
/// of `0`, `simple`, or `verbose`.
fn parse_header_arg(header: &[Vec<u8>]) -> io::Result<Header> {
    if header.is_empty() {
        return Ok(Header::None);
    }
    if header.len() > 1 {
        return Err(raw_error(
            b"Invalid header format: ",
            &header.join(&b' '),
            b"",
        ));
    }
    match header[0].as_slice() {
        b"" | b"verbose" => Ok(Header::Verbose),
        b"0" => Ok(Header::None),
        b"simple" => Ok(Header::Simple),
        value => Err(raw_error(
            b"Invalid value for string field: ",
            value,
            b". Permitted values: 0, simple, verbose",
        )),
    }
}

fn reject_unsupported_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    if config.report_unaligned && !runtime.source_extra {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Option is not permitted for this workflow: unal",
        ));
    }
    Ok(())
}

/// Source ViewWriter owns an OutputFile and closes it even during unwind.
/// Actual descriptor adoption follows standalone CLI versus embedded policy.
fn open_output_writer(
    path: &Path,
    file_buffer_size: usize,
    runtime: &RuntimeOptions,
) -> io::Result<crate::output::daa::buffered_io::DaaOutput> {
    let compressor = crate::basic::config::compressor_bytes(&runtime.compression)?;
    crate::output::daa::buffered_io::DaaOutput::create_view_with_compressor(
        path,
        file_buffer_size,
        compressor,
    )
}

const VIEW_QUERY_BATCH_SIZE: usize = 32;

type ViewQueryBuffer = (Vec<u8>, usize);

fn read_view_batch(daa: &mut DaaFile) -> io::Result<(Option<Vec<ViewQueryBuffer>>, bool)> {
    let mut batch = Vec::with_capacity(VIEW_QUERY_BATCH_SIZE);
    let mut at_end = false;
    for _ in 0..VIEW_QUERY_BATCH_SIZE {
        match daa.read_query_buffer()? {
            Some(query) => batch.push(query),
            None => {
                at_end = true;
                break;
            }
        }
    }
    Ok(((!batch.is_empty()).then_some(batch), at_end))
}

/// Source renders the first record before printing the format header or
/// starting workers. A malformed first alignment must not leave a header.
fn render_first_view_query<Render>(
    daa: &mut DaaFile,
    render: &Render,
) -> io::Result<Option<(DaaQueryRecord, Vec<u8>)>>
where
    Render: Fn(&DaaQueryRecord, &DaaFile, &mut Vec<u8>) -> Result<(), String>,
{
    let Some((buffer, ordinal)) = daa.read_query_buffer()? else {
        return Ok(None);
    };
    let record = DaaQueryRecord::from_buffer_bytes(daa, &buffer, ordinal)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    let mut output = Vec::new();
    render(&record, daa, &mut output)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    Ok(Some((record, output)))
}

/// Headerless formats still execute the first record synchronously.
fn run_ordered_view_workers<Render, WriteBuffer>(
    daa: &mut DaaFile,
    daa_path: &Path,
    requested_threads: usize,
    file_buffer_size: usize,
    render: Render,
    mut write_buffer: WriteBuffer,
) -> io::Result<()>
where
    Render: Fn(&DaaQueryRecord, &DaaFile, &mut Vec<u8>) -> Result<(), String> + Sync,
    WriteBuffer: FnMut(&[u8]) -> io::Result<()>,
{
    if let Some((_, output)) = render_first_view_query(daa, &render)? {
        write_buffer(&output)?;
        run_remaining_view_workers(
            daa,
            daa_path,
            requested_threads,
            file_buffer_size,
            render,
            write_buffer,
        )?;
    }
    Ok(())
}

/// Run the pinned `view_worker` boundary: 32 query records per task and at
/// most `3 * threads` ordered tasks resident.  Each worker renders a whole
/// task into one buffer, and the caller writes completed buffers in input
/// order.  The mutable archive reader remains on this thread; workers use a
/// second handle only for immutable header/reference metadata.
fn run_remaining_view_workers<Render, WriteBuffer>(
    daa: &mut DaaFile,
    daa_path: &Path,
    requested_threads: usize,
    file_buffer_size: usize,
    render: Render,
    mut write_buffer: WriteBuffer,
) -> io::Result<()>
where
    Render: Fn(&DaaQueryRecord, &DaaFile, &mut Vec<u8>) -> Result<(), String> + Sync,
    WriteBuffer: FnMut(&[u8]) -> io::Result<()>,
{
    let threads = view_threads(requested_threads);
    let queue_limit = 3usize.saturating_mul(threads).max(1);
    let worker_daa = DaaFile::open_bytes_with_buffer_size(daa_path, file_buffer_size)?;
    let (job_tx, job_rx) = mpsc::sync_channel::<(usize, Vec<ViewQueryBuffer>)>(queue_limit);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (result_tx, result_rx) = mpsc::channel::<(usize, Result<Vec<u8>, String>)>();

    std::thread::scope(|scope| -> io::Result<()> {
        for _ in 0..threads {
            let job_rx = Arc::clone(&job_rx);
            let result_tx = result_tx.clone();
            let worker_daa = &worker_daa;
            let render = &render;
            let worker = std::thread::Builder::new().spawn_scoped(scope, move || loop {
                let job = match job_rx.lock() {
                    Ok(rx) => rx.recv(),
                    Err(_) => return,
                };
                let Ok((sequence, batch)) = job else {
                    return;
                };
                let rendered = (|| {
                    let mut output = Vec::new();
                    for (buf, query_num) in batch {
                        let record =
                            DaaQueryRecord::from_buffer_bytes(worker_daa, &buf, query_num)?;
                        render(&record, worker_daa, &mut output)?;
                    }
                    Ok(output)
                })();
                if result_tx.send((sequence, rendered)).is_err() {
                    return;
                }
            });
            if let Err(error) = worker {
                // Wake already-created workers before scope waits for them.
                // C++ thread construction throws; surface a fallible native
                // resource error instead of panicking or stranding receivers.
                drop(job_tx);
                return Err(error);
            }
        }
        drop(result_tx);

        let mut pending = BTreeMap::new();
        let mut next_submit = 0usize;
        let mut next_write = 0usize;
        let mut resident = 0usize;
        let mut at_end = false;

        let result = (|| -> io::Result<()> {
            while resident > 0 || !at_end {
                while resident < queue_limit && !at_end {
                    let (batch, reached_end) = read_view_batch(daa)?;
                    at_end = reached_end;
                    match batch {
                        Some(batch) => {
                            job_tx.send((next_submit, batch)).map_err(|_| {
                                io::Error::new(
                                    io::ErrorKind::BrokenPipe,
                                    "view worker queue closed",
                                )
                            })?;
                            next_submit += 1;
                            resident += 1;
                        }
                        None => {}
                    }
                }
                if resident == 0 {
                    break;
                }

                let (sequence, rendered) = result_rx.recv().map_err(|_| {
                    io::Error::new(io::ErrorKind::BrokenPipe, "view result queue closed")
                })?;
                pending.insert(sequence, rendered);
                while let Some(rendered) = pending.remove(&next_write) {
                    let output =
                        rendered.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                    write_buffer(&output)?;
                    next_write += 1;
                    resident -= 1;
                }
            }
            Ok(())
        })();
        drop(job_tx);
        result
    })
}

/// Run the view command — show DAA file metadata.
///
/// Record renderers below use the native DAA decoding and output modules.
pub fn run(daa_file: &str) -> io::Result<()> {
    run_path(Path::new(daa_file))
}

pub fn run_path(daa_file: &Path) -> io::Result<()> {
    let daa = DaaFile::open_bytes(daa_file)?;

    let mode_str = match daa.mode() {
        2 => "blastp",
        3 => "blastx",
        4 => "blastn",
        _ => "unknown",
    };

    // This legacy inspection helper is a native library convenience, not
    // Source view_daa's format writer. Retain its report, including raw names.
    let mut diagnostics = crate::util::io::console::stderr();
    diagnostics.write_all(b"DAA file: ")?;
    diagnostics.write_all(
        &crate::util::io::source_os_bytes(daa_file.as_os_str()).map_err(io::Error::other)?,
    )?;
    diagnostics.write_all(b"\n")?;
    writeln!(diagnostics, "Mode: {}", mode_str)?;
    writeln!(diagnostics, "Database sequences: {}", daa.db_seqs())?;
    writeln!(
        diagnostics,
        "Database sequences used: {}",
        daa.db_seqs_used()
    )?;
    writeln!(diagnostics, "Database letters: {}", daa.db_letters())?;
    writeln!(diagnostics, "Query records: {}", daa.query_records())?;
    diagnostics.write_all(b"Matrix: ")?;
    diagnostics.write_all(daa.score_matrix_bytes())?;
    diagnostics.write_all(b"\n")?;
    writeln!(
        diagnostics,
        "Gap penalties: {}/{}",
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    )?;
    if daa.db_seqs_used() > 0 {
        diagnostics.write_all(b"First reference: ")?;
        diagnostics.write_all(daa.ref_name_bytes(0))?;
        writeln!(diagnostics, " ({})", daa.ref_len_at(0))?;
    }

    Ok(())
}

/// Matches the DAA-output branch of C++ `view_daa`.
pub fn run_daa(config: &ViewConfig) -> io::Result<()> {
    run_daa_with_options(config, &RuntimeOptions::default())
}
fn run_daa_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    reject_unsupported_options(config, runtime)?;
    let mut daa = DaaFile::open_bytes_with_buffer_size(&config.daa_file, config.file_buffer_size)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix_for_alignment()?,
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_daa_metadata(&daa, &score_matrix, config.verbosity)?;
    let source_output = prepare_source_output(config, runtime, &daa)?;
    print_view_selection(
        config,
        source_output
            .as_ref()
            .map_or(HspValues::TRANSCRIPT, |o| o.dp_fields),
        runtime,
    )?;

    let _taxonomy_owner = load_view_taxonomy(config, runtime)?;
    let mut out = open_output_writer(&config.output, config.file_buffer_size, runtime)?;
    init_daa(&mut out)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        config.file_buffer_size,
        |record, worker_daa, query_out| {
            view_query_daa(record, worker_daa, query_out, &score_matrix, &view_config)
        },
        |query_out| out.write_all(query_out),
    )?;
    finish_daa_from_file(&mut out, &daa)?;
    out.close()
}

/// Matches the tabular-output branch of C++ `view_daa`.
pub fn run_tabular(config: &ViewConfig) -> io::Result<()> {
    run_tabular_with_options(config, &RuntimeOptions::default())
}
fn run_tabular_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    reject_unsupported_options(config, runtime)?;
    let mut daa = DaaFile::open_bytes_with_buffer_size(&config.daa_file, config.file_buffer_size)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix_for_alignment()?,
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    print_daa_metadata(&daa, &score_matrix, config.verbosity)?;
    let outfmt = if config.outfmt.is_empty() {
        vec![b"tab".to_vec()]
    } else {
        config.outfmt.clone()
    };
    let is_json = outfmt
        .first()
        .is_some_and(|f| f == b"json-flat" || f == b"104");
    let outfmt_refs = outfmt.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let source_output = prepare_source_output(config, runtime, &daa)?;
    let format = if let Some(output) = source_output {
        match output.output_format {
            crate::output::format::ResolvedOutputFormat::Tabular(format) => format,
            _ => return Err(io::Error::other("Expected tabular view format")),
        }
    } else {
        TabularFormat::new_bytes(
            &outfmt_refs,
            is_json,
            runtime.frame_shift,
            daa.mode() == 3,
            daa.score_matrix_for_alignment()?.as_bytes(),
        )?
    };
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_view_selection(config, format.hsp_values, runtime)?;

    let _taxonomy_owner = load_view_taxonomy(config, runtime)?;
    let mut out = open_output_writer(&config.output, config.file_buffer_size, runtime)?;
    let render = |record: &DaaQueryRecord, worker_daa: &DaaFile, query_out: &mut Vec<u8>| {
        crate::data::daa::view_query_tabular_with_options(
            record,
            worker_daa,
            query_out,
            &score_matrix,
            &format,
            &view_config,
            config.report_unaligned && !runtime.source_extra,
            runtime.frame_shift != 0,
        )
    };
    let first = render_first_view_query(&mut daa, &render)?;
    let header_mode = parse_header_arg(&config.header)?;
    let header = print_tabular_header_bytes(
        &format.fields,
        header_mode,
        format.is_json,
        VERSION_STRING,
        &config.invocation,
    )
    .map_err(io::Error::other)?;
    out.write_all(&header)?;
    if let Some((_, query_out)) = first {
        out.write_all(&query_out)?;
        run_remaining_view_workers(
            &mut daa,
            &config.daa_file,
            config.threads,
            config.file_buffer_size,
            render,
            |query_out| out.write_all(query_out),
        )?;
    }
    out.write_all(print_tabular_footer(format.is_json).as_bytes())?;
    out.close()
}

/// Matches the PAF-output branch of C++ `view_daa`.
pub fn run_paf(config: &ViewConfig) -> io::Result<()> {
    run_paf_with_options(config, &RuntimeOptions::default())
}
fn run_paf_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    reject_unsupported_options(config, runtime)?;
    let mut daa = DaaFile::open_bytes_with_buffer_size(&config.daa_file, config.file_buffer_size)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix_for_alignment()?,
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_daa_metadata(&daa, &score_matrix, config.verbosity)?;
    let source_output = prepare_source_output(config, runtime, &daa)?;
    print_view_selection(
        config,
        source_output
            .as_ref()
            .map_or(HspValues::TRANSCRIPT, |o| o.dp_fields),
        runtime,
    )?;

    let _taxonomy_owner = load_view_taxonomy(config, runtime)?;
    let mut out = open_output_writer(&config.output, config.file_buffer_size, runtime)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        config.file_buffer_size,
        |record, worker_daa, query_out| {
            view_query_paf(record, worker_daa, query_out, &score_matrix, &view_config)
        },
        |query_out| out.write_all(query_out),
    )?;
    out.close()
}

/// Matches the SAM-output branch of C++ `view_daa`.
pub fn run_sam(config: &ViewConfig) -> io::Result<()> {
    run_sam_with_options(config, &RuntimeOptions::default())
}
fn run_sam_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    reject_unsupported_options(config, runtime)?;
    let mut daa = DaaFile::open_bytes_with_buffer_size(&config.daa_file, config.file_buffer_size)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix_for_alignment()?,
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_daa_metadata(&daa, &score_matrix, config.verbosity)?;
    let source_output = prepare_source_output(config, runtime, &daa)?;
    print_view_selection(
        config,
        source_output
            .as_ref()
            .map_or(HspValues::TRANSCRIPT, |o| o.dp_fields),
        runtime,
    )?;

    let _taxonomy_owner = load_view_taxonomy(config, runtime)?;
    let mut out = open_output_writer(&config.output, config.file_buffer_size, runtime)?;
    let render = |record: &DaaQueryRecord, worker_daa: &DaaFile, query_out: &mut Vec<u8>| {
        view_query_sam(
            record,
            worker_daa,
            query_out,
            &score_matrix,
            &view_config,
            config.sam_qlen_field,
        )
    };
    let first = render_first_view_query(&mut daa, &render)?;
    sam::print_header_bytes(&mut out, daa.mode(), &config.invocation)?;
    if let Some((_, query_out)) = first {
        out.write_all(&query_out)?;
        run_remaining_view_workers(
            &mut daa,
            &config.daa_file,
            config.threads,
            config.file_buffer_size,
            render,
            |query_out| out.write_all(query_out),
        )?;
    }

    out.close()
}

/// Matches the pairwise-output branch of C++ `view_daa`.
pub fn run_pairwise(config: &ViewConfig) -> io::Result<()> {
    run_pairwise_with_options(config, &RuntimeOptions::default())
}
fn run_pairwise_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    reject_unsupported_options(config, runtime)?;
    let mut daa = DaaFile::open_bytes_with_buffer_size(&config.daa_file, config.file_buffer_size)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix_for_alignment()?,
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_daa_metadata(&daa, &score_matrix, config.verbosity)?;
    let source_output = prepare_source_output(config, runtime, &daa)?;
    print_view_selection(
        config,
        source_output
            .as_ref()
            .map_or(HspValues::TRANSCRIPT, |o| o.dp_fields),
        runtime,
    )?;

    let _taxonomy_owner = load_view_taxonomy(config, runtime)?;
    let mut out = open_output_writer(&config.output, config.file_buffer_size, runtime)?;
    let render = |record: &DaaQueryRecord, worker_daa: &DaaFile, query_out: &mut Vec<u8>| {
        view_query_pairwise(record, worker_daa, query_out, &score_matrix, &view_config)
    };
    let first = render_first_view_query(&mut daa, &render)?;
    pairwise::write_header(&mut out)?;
    if let Some((_, query_out)) = first {
        out.write_all(&query_out)?;
        run_remaining_view_workers(
            &mut daa,
            &config.daa_file,
            config.threads,
            config.file_buffer_size,
            render,
            |query_out| out.write_all(query_out),
        )?;
    }
    pairwise::print_footer(&mut out)?;
    out.close()
}

/// Matches the XML-output branch of C++ `view_daa`.
pub fn run_xml(config: &ViewConfig) -> io::Result<()> {
    run_xml_with_options(config, &RuntimeOptions::default())
}
fn run_xml_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    reject_unsupported_options(config, runtime)?;
    let mut daa = DaaFile::open_bytes_with_buffer_size(&config.daa_file, config.file_buffer_size)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix_for_alignment()?,
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_daa_metadata(&daa, &score_matrix, config.verbosity)?;
    let source_output = prepare_source_output(config, runtime, &daa)?;
    print_view_selection(
        config,
        source_output
            .as_ref()
            .map_or(HspValues::TRANSCRIPT, |o| o.dp_fields),
        runtime,
    )?;

    let _taxonomy_owner = load_view_taxonomy(config, runtime)?;
    let mut out = open_output_writer(&config.output, config.file_buffer_size, runtime)?;
    if let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer_bytes(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut query_out = Vec::new();
        crate::data::daa::view_query_xml_with_options(
            &r,
            &daa,
            &mut query_out,
            &score_matrix,
            &view_config,
            runtime.xml_blord_format,
            runtime.no_parse_seqids,
        )
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        xml::print_header_bytes(
            &mut out,
            daa.mode(),
            &daa.score_matrix_for_alignment()?,
            daa.gap_open_penalty(),
            daa.gap_extension_penalty(),
            daa.evalue(),
            &r.query_name,
            r.query_len() as u32,
            &crate::util::io::source_os_bytes(runtime.database.as_os_str())
                .map_err(io::Error::other)?,
        )?;
        out.write_all(&query_out)?;
        run_remaining_view_workers(
            &mut daa,
            &config.daa_file,
            config.threads,
            config.file_buffer_size,
            |record, worker_daa, query_out| {
                crate::data::daa::view_query_xml_with_options(
                    record,
                    worker_daa,
                    query_out,
                    &score_matrix,
                    &view_config,
                    runtime.xml_blord_format,
                    runtime.no_parse_seqids,
                )
            },
            |query_out| out.write_all(query_out),
        )?;
    } else {
        xml::print_header_bytes(
            &mut out,
            daa.mode(),
            &daa.score_matrix_for_alignment()?,
            daa.gap_open_penalty(),
            daa.gap_extension_penalty(),
            daa.evalue(),
            b"",
            0,
            &crate::util::io::source_os_bytes(runtime.database.as_os_str())
                .map_err(io::Error::other)?,
        )?;
    }
    xml::print_footer(&mut out)?;
    out.close()
}

/// Matches the null-output branch of C++ `view_daa`.
pub fn run_null(config: &ViewConfig) -> io::Result<()> {
    run_null_with_options(config, &RuntimeOptions::default())
}
fn run_null_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    reject_unsupported_options(config, runtime)?;
    let mut daa = DaaFile::open_bytes_with_buffer_size(&config.daa_file, config.file_buffer_size)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix_for_alignment()?,
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_daa_metadata(&daa, &score_matrix, config.verbosity)?;
    let source_output = prepare_source_output(config, runtime, &daa)?;
    print_view_selection(
        config,
        source_output
            .as_ref()
            .map_or(HspValues::TRANSCRIPT, |o| o.dp_fields),
        runtime,
    )?;

    let _taxonomy_owner = load_view_taxonomy(config, runtime)?;
    let mut out = open_output_writer(&config.output, config.file_buffer_size, runtime)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        config.file_buffer_size,
        |record, worker_daa, _| view_query_null(record, worker_daa, &score_matrix, &view_config),
        |_| Ok(()),
    )?;
    out.close()
}

/// Matches the edge-output branch of C++ `view_daa`.
pub fn run_edge(config: &ViewConfig) -> io::Result<()> {
    run_edge_with_options(config, &RuntimeOptions::default())
}
fn run_edge_with_options(config: &ViewConfig, runtime: &RuntimeOptions) -> io::Result<()> {
    reject_unsupported_options(config, runtime)?;
    let mut daa = DaaFile::open_bytes_with_buffer_size(&config.daa_file, config.file_buffer_size)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix_for_alignment()?,
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_daa_metadata(&daa, &score_matrix, config.verbosity)?;
    let source_output = prepare_source_output(config, runtime, &daa)?;
    print_view_selection(
        config,
        source_output
            .as_ref()
            .map_or(HspValues::COORDS, |o| o.dp_fields),
        runtime,
    )?;

    let _taxonomy_owner = load_view_taxonomy(config, runtime)?;
    let mut out = open_output_writer(&config.output, config.file_buffer_size, runtime)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        config.file_buffer_size,
        |record, worker_daa, query_out| {
            view_query_edge(record, worker_daa, query_out, &score_matrix, &view_config)
        },
        |query_out| out.write_all(query_out),
    )?;
    out.close()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::{EditOperation, PackedOperation};
    use crate::basic::value::SequenceType;
    use crate::data::daa::{
        compute_flag, write_daa_query_record, BlockType, DaaHeader1, DaaHeader2,
    };
    use std::io::Write;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static VIEW_TEST_FILE_ID: AtomicUsize = AtomicUsize::new(0);

    #[test]
    fn test_view_stub() {
        // Historical metadata-helper placeholder; actual native archive and
        // record rendering is exercised by the owner-based tests below.
    }

    #[cfg(feature = "extra")]
    #[test]
    fn source_extra_view_retains_names_and_typed_input_error() {
        let id = VIEW_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "diamond_rs_view_names_{}_{}.dmp",
            std::process::id(),
            id
        ));
        std::fs::write(
            &path,
            b"1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        let runtime = RuntimeOptions {
            source_extra: true,
            taxonomy_names: path.clone(),
            ..RuntimeOptions::default()
        };
        let owner = load_view_taxonomy(&ViewConfig::default(), &runtime).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(
            owner.name,
            vec![b"".to_vec(), b"root".to_vec(), b"Bacteria".to_vec()]
        );
        let error = load_view_taxonomy(&ViewConfig::default(), &runtime).unwrap_err();
        assert!(
            crate::util::io::file_open_error::find(&error).is_some(),
            "{error}"
        );
        // The returned owner remains valid after its source descriptor closes.
        assert_eq!(owner.name[2], b"Bacteria");
    }

    #[test]
    fn source_view_translated_mate_field_is_absent_not_unavailable() {
        let id = VIEW_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let directory =
            std::env::temp_dir().join(format!("diamond_view_mate_{}_{}", std::process::id(), id));
        std::fs::create_dir(&directory).unwrap();
        let input = directory.join("input.daa");
        let output = directory.join("output.tsv");
        std::fs::write(
            &input,
            include_bytes!("../../tests/data/min_score/blastx.cbs1.score100.daa"),
        )
        .unwrap();
        run_tabular(&ViewConfig {
            daa_file: input.clone(),
            output: output.clone(),
            outfmt: vec![
                "6".into(),
                "qseqid".into(),
                "full_qseq_mate".into(),
                "qseq_translated".into(),
            ],
            max_target_seqs: 25,
            threads: 1,
            ..ViewConfig::default()
        })
        .unwrap();
        assert_eq!(
            std::fs::read(&output).unwrap(),
            include_bytes!("../../tests/data/view_extra/unpaired-translated.tsv")
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn test_run_daa_rewrites_archive() {
        let id = VIEW_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let input = std::env::temp_dir().join(format!(
            "diamond_rs_view_input_{}_{}.daa",
            std::process::id(),
            id
        ));
        let output = std::env::temp_dir().join(format!(
            "diamond_rs_view_output_{}_{}.daa",
            std::process::id(),
            id
        ));

        let mut query = Vec::new();
        let seek_pos =
            write_daa_query_record(&mut query, "query0", &[0, 1, 2], SequenceType::AminoAcid);
        crate::data::daa::finish_daa_query_record(&mut query, seek_pos);
        query.extend_from_slice(&0u32.to_ne_bytes());

        let mut h2 = DaaHeader2::with_params(1, 100, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);
        h2.db_seqs_used = 1;
        h2.query_records = 1;
        h2.block_type[0] = BlockType::Alignments as u8;
        h2.block_type[1] = BlockType::RefNames as u8;
        h2.block_type[2] = BlockType::RefLengths as u8;
        h2.block_size[0] = query.len() as u64;
        h2.block_size[1] = 5;
        h2.block_size[2] = 4;

        let mut f = File::create(&input).unwrap();
        DaaHeader1::new().write_to(&mut f).unwrap();
        h2.write_to(&mut f).unwrap();
        f.write_all(&query).unwrap();
        f.write_all(b"ref0\0").unwrap();
        f.write_all(&100u32.to_ne_bytes()).unwrap();

        run_daa(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output.clone(),
            outfmt: vec![b"daa".to_vec()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: Vec::new(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        let mut out_daa = DaaFile::open_bytes(&output).unwrap();
        assert_eq!(out_daa.query_records(), 1);
        assert_eq!(out_daa.db_seqs_used(), 1);
        assert_eq!(out_daa.ref_name_bytes(0), b"ref0");
        let (buf, query_num) = out_daa.read_query_buffer().unwrap().unwrap();
        let r = DaaQueryRecord::from_buffer_bytes(&out_daa, &buf, query_num).unwrap();
        assert_eq!(r.query_name, b"query0");
        assert!(out_daa.read_query_buffer().unwrap().is_none());

        let _ = std::fs::remove_file(input);
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn test_run_tabular_writes_alignment_rows() {
        let id = VIEW_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let input = std::env::temp_dir().join(format!(
            "diamond_rs_view_tab_input_{}_{}.daa",
            std::process::id(),
            id
        ));
        let output = std::env::temp_dir().join(format!(
            "diamond_rs_view_tab_output_{}_{}.txt",
            std::process::id(),
            id
        ));
        let output_json = std::env::temp_dir().join(format!(
            "diamond_rs_view_json_output_{}_{}.json",
            std::process::id(),
            id
        ));
        let output_paf = std::env::temp_dir().join(format!(
            "diamond_rs_view_paf_output_{}_{}.paf",
            std::process::id(),
            id
        ));
        let output_sam = std::env::temp_dir().join(format!(
            "diamond_rs_view_sam_output_{}_{}.sam",
            std::process::id(),
            id
        ));
        let output_pairwise = std::env::temp_dir().join(format!(
            "diamond_rs_view_pairwise_output_{}_{}.txt",
            std::process::id(),
            id
        ));
        let output_xml = std::env::temp_dir().join(format!(
            "diamond_rs_view_xml_output_{}_{}.xml",
            std::process::id(),
            id
        ));
        let output_null = std::env::temp_dir().join(format!(
            "diamond_rs_view_null_output_{}_{}.txt",
            std::process::id(),
            id
        ));
        let output_edge = std::env::temp_dir().join(format!(
            "diamond_rs_view_edge_output_{}_{}.bin",
            std::process::id(),
            id
        ));

        let mut query = Vec::new();
        let seek_pos =
            write_daa_query_record(&mut query, "query0", &[0, 1, 2], SequenceType::AminoAcid);
        let flag = compute_flag(30, 0, 0, false);
        query.extend_from_slice(&0u32.to_ne_bytes());
        query.push(flag);
        query.push(30);
        query.push(0);
        query.push(0);
        query.push(PackedOperation::from_op_count(EditOperation::Match, 3).code);
        query.push(PackedOperation::terminator().code);
        crate::data::daa::finish_daa_query_record(&mut query, seek_pos);
        query.extend_from_slice(&0u32.to_ne_bytes());

        let mut h2 = DaaHeader2::with_params(1, 100, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);
        h2.db_seqs_used = 1;
        h2.query_records = 1;
        h2.block_type[0] = BlockType::Alignments as u8;
        h2.block_type[1] = BlockType::RefNames as u8;
        h2.block_type[2] = BlockType::RefLengths as u8;
        h2.block_size[0] = query.len() as u64;
        h2.block_size[1] = 5;
        h2.block_size[2] = 4;

        let mut f = File::create(&input).unwrap();
        DaaHeader1::new().write_to(&mut f).unwrap();
        h2.write_to(&mut f).unwrap();
        f.write_all(&query).unwrap();
        f.write_all(b"ref0\0").unwrap();
        f.write_all(&100u32.to_ne_bytes()).unwrap();

        run_tabular(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output.clone(),
            outfmt: vec![
                b"tab".to_vec(),
                b"qseqid".to_vec(),
                b"sseqid".to_vec(),
                b"score".to_vec(),
                b"length".to_vec(),
                b"pident".to_vec(),
                b"hspnum".to_vec(),
                b"qnum".to_vec(),
                b"snum".to_vec(),
            ],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: Vec::new(),
            // Exercise the ordered parallel tabular path; the remaining
            // format cases below retain the direct one-thread path.
            threads: 4,
            qnum_offset: (-1_i64) as u64,
            snum_offset: 20,
            verbosity: 0,
        })
        .unwrap();

        let out = std::fs::read_to_string(&output).unwrap();
        assert_eq!(
            out,
            "query0\tref0\t30\t3\t100\t0\t18446744073709551615\t20\n"
        );

        run_tabular(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output_json.clone(),
            outfmt: vec![
                b"104".to_vec(),
                b"qseqid".to_vec(),
                b"sseqid".to_vec(),
                b"score".to_vec(),
                b"qnum".to_vec(),
                b"snum".to_vec(),
            ],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: Vec::new(),
            threads: 1,
            qnum_offset: 10,
            snum_offset: (-2_i64) as u64,
            verbosity: 0,
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_json).unwrap();
        assert_eq!(
            out,
            concat!(
                "[\n\t{\n",
                "\t\"qseqid\":\"query0\",\n",
                "\t\"sseqid\":\"ref0\",\n",
                "\t\"score\":30,\n",
                "\t\"qnum\":10,\n",
                "\t\"snum\":18446744073709551614\n",
                "\t}\n]"
            )
        );

        run_paf(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output_paf.clone(),
            outfmt: vec![b"103".to_vec()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: Vec::new(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_paf).unwrap();
        assert!(out.starts_with("query0\t3\t0\t2\t+\tref0\t100\t0\t2\t3\t3\t255\tAS:i:"));
        assert!(out.contains("\tZR:i:30\tZE:f:"));

        run_sam(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output_sam.clone(),
            outfmt: vec![b"101".to_vec()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: true,
            header: vec![],
            invocation: b"diamond view -f sam".to_vec(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_sam).unwrap();
        assert!(out.starts_with("@HD\tVN:1.5\tSO:query\n"));
        assert!(out.contains("@PG\tPN:DIAMOND\tVN:2.1.24\tCL:diamond view -f sam\n"));
        assert!(out.contains("\nquery0\t0\tref0\t1\t255\t3M\t*\t0\t0\tARN\t*"));
        assert!(out.contains("\tNM:i:0\tZL:i:100\tZR:i:30\t"));
        assert!(out.contains("\tZI:i:100\tZF:i:1\tZS:i:1\tMD:Z:3\tZQ:i:3\n"));

        run_pairwise(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output_pairwise.clone(),
            outfmt: vec![b"0".to_vec()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: Vec::new(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_pairwise).unwrap();
        assert!(out.starts_with("BLASTP 2.3.0+\n\n\nQuery= query0\n\nLength=3\n\n"));
        assert!(out.contains(">ref0\nLength=100\n\n"));
        assert!(out.contains("Score = "));
        assert!(out.contains("Identities = 3/3 (100%), Positives = 3/3 (100%), Gaps = 0/3 (0%)"));
        assert!(out.contains("Query  1  ARN 3\n"));
        assert!(out.contains("Sbjct  1  ARN 3\n"));

        run_xml(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output_xml.clone(),
            outfmt: vec![b"5".to_vec()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: Vec::new(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_xml).unwrap();
        assert!(out.starts_with("<?xml version=\"1.0\"?>\n"));
        assert!(out.contains("<BlastOutput_program>blastp</BlastOutput_program>"));
        assert!(out.contains("<BlastOutput_query-def>query0</BlastOutput_query-def>"));
        assert!(out.contains("<Iteration_query-ID>Query_1</Iteration_query-ID>"));
        assert!(out.contains("<Hit_id>ref0</Hit_id>"));
        assert!(out.contains("<Hsp_qseq>ARN</Hsp_qseq>"));
        assert!(out.contains("<Hsp_hseq>ARN</Hsp_hseq>"));
        assert!(out.ends_with("</BlastOutput_iterations>\n</BlastOutput>"));

        run_null(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output_null.clone(),
            outfmt: vec![b"null".to_vec()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: Vec::new(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        assert!(std::fs::read(&output_null).unwrap().is_empty());

        run_edge(&ViewConfig {
            file_buffer_size: 67_108_864,
            daa_file: input.clone(),
            output: output_edge.clone(),
            outfmt: vec![b"edge".to_vec()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: Vec::new(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        let out = std::fs::read(&output_edge).unwrap();
        assert_eq!(out.len(), crate::output::edge::EdgeData::SIZE);
        assert_eq!(u64::from_ne_bytes(out[0..8].try_into().unwrap()), 0);
        assert_eq!(u64::from_ne_bytes(out[8..16].try_into().unwrap()), 0);
        assert_eq!(f32::from_ne_bytes(out[16..20].try_into().unwrap()), 100.0);
        assert_eq!(f32::from_ne_bytes(out[20..24].try_into().unwrap()), 3.0);

        let _ = std::fs::remove_file(input);
        let _ = std::fs::remove_file(output);
        let _ = std::fs::remove_file(output_json);
        let _ = std::fs::remove_file(output_paf);
        let _ = std::fs::remove_file(output_sam);
        let _ = std::fs::remove_file(output_pairwise);
        let _ = std::fs::remove_file(output_xml);
        let _ = std::fs::remove_file(output_null);
        let _ = std::fs::remove_file(output_edge);
    }

    #[test]
    fn test_run_tabular_json_preserves_upstream_query_boundary_bytes() {
        let id = VIEW_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let input = std::env::temp_dir().join(format!(
            "diamond_rs_view_multi_json_input_{}_{}.daa",
            std::process::id(),
            id
        ));
        let output = std::env::temp_dir().join(format!(
            "diamond_rs_view_multi_json_output_{}_{}.json",
            std::process::id(),
            id
        ));

        let mut alignments = Vec::new();
        for (name, score) in [("query0", 30u8), ("query1", 29u8)] {
            let seek_pos =
                write_daa_query_record(&mut alignments, name, &[0, 1, 2], SequenceType::AminoAcid);
            let flag = compute_flag(u32::from(score), 0, 0, false);
            alignments.extend_from_slice(&0u32.to_ne_bytes());
            alignments.push(flag);
            alignments.push(score);
            alignments.push(0);
            alignments.push(0);
            alignments.push(PackedOperation::from_op_count(EditOperation::Match, 3).code);
            alignments.push(PackedOperation::terminator().code);
            crate::data::daa::finish_daa_query_record(&mut alignments, seek_pos);
        }
        alignments.extend_from_slice(&0u32.to_ne_bytes());

        let mut h2 = DaaHeader2::with_params(1, 100, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);
        h2.db_seqs_used = 1;
        h2.query_records = 2;
        h2.block_type[0] = BlockType::Alignments as u8;
        h2.block_type[1] = BlockType::RefNames as u8;
        h2.block_type[2] = BlockType::RefLengths as u8;
        h2.block_size[0] = alignments.len() as u64;
        h2.block_size[1] = 5;
        h2.block_size[2] = 4;

        let mut file = File::create(&input).unwrap();
        DaaHeader1::new().write_to(&mut file).unwrap();
        h2.write_to(&mut file).unwrap();
        file.write_all(&alignments).unwrap();
        file.write_all(b"ref0\0").unwrap();
        file.write_all(&100u32.to_ne_bytes()).unwrap();
        drop(file);

        run_tabular(&ViewConfig {
            daa_file: input.clone(),
            output: output.clone(),
            outfmt: vec![
                b"104".to_vec(),
                b"qseqid".to_vec(),
                b"sseqid".to_vec(),
                b"score".to_vec(),
            ],
            max_target_seqs: 25,
            ..ViewConfig::default()
        })
        .unwrap();

        let actual = std::fs::read_to_string(&output).unwrap();
        assert_eq!(
            actual,
            concat!(
                "[\n\t{\n",
                "\t\"qseqid\":\"query0\",\n",
                "\t\"sseqid\":\"ref0\",\n",
                "\t\"score\":30\n",
                "\t}\n\t{\n",
                "\t\"qseqid\":\"query1\",\n",
                "\t\"sseqid\":\"ref0\",\n",
                "\t\"score\":29\n",
                "\t}\n]"
            )
        );
        assert!(actual.contains("\t}\n\t{\n"));

        let _ = std::fs::remove_file(input);
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn test_view_rejects_unaligned_reporting_before_opening_archive() {
        let error = run_tabular(&ViewConfig {
            daa_file: "does-not-exist.daa".into(),
            report_unaligned: true,
            ..ViewConfig::default()
        })
        .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(
            error.to_string(),
            "Option is not permitted for this workflow: unal"
        );
    }

    #[test]
    fn raw_view_header_error_retains_source_option_bytes() {
        let error = parse_header_arg(&[b"bad_\xff".to_vec()]).unwrap_err();
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error),
            Some(
                b"Invalid value for string field: bad_\xff. Permitted values: 0, simple, verbose"
                    .as_slice()
            )
        );
    }

    #[test]
    fn test_view_header_modes_match_upstream_tabular_rules() {
        assert_eq!(parse_header_arg(&[]).unwrap(), Header::None);
        assert_eq!(parse_header_arg(&[Vec::new()]).unwrap(), Header::Verbose);
        assert_eq!(
            parse_header_arg(&[b"simple".to_vec()]).unwrap(),
            Header::Simple
        );
        assert_eq!(parse_header_arg(&[b"0".to_vec()]).unwrap(), Header::None);
        assert_eq!(
            parse_header_arg(&[b"bad".to_vec()])
                .unwrap_err()
                .to_string(),
            "Invalid value for string field: bad. Permitted values: 0, simple, verbose"
        );
        assert_eq!(
            parse_header_arg(&[b"simple".to_vec(), b"verbose".to_vec()])
                .unwrap_err()
                .to_string(),
            "Invalid header format: simple verbose"
        );
    }
}
