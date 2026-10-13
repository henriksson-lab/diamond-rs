use std::collections::BTreeMap;
use std::fs::File;
#[cfg(not(windows))]
use std::fs::OpenOptions;
use std::io::{self, BufWriter, Write};
use std::sync::{mpsc, Arc, Mutex};

use crate::basic::consts::VERSION_STRING;
use crate::data::daa::{
    finish_daa_from_file, init_daa, view_query_daa, view_query_edge, view_query_null,
    view_query_paf, view_query_pairwise, view_query_sam, view_query_tabular, view_query_xml,
    DaaFile, DaaQueryRecord, DaaViewConfig,
};
use crate::dp::swipe::HspValues;
use crate::output::format::{print_tabular_footer, print_tabular_header, Header, TabularFormat};
use crate::output::{pairwise, sam, xml};
use crate::stats::score_matrix::ScoreMatrix;

#[derive(Debug, Clone, Default)]
pub struct ViewConfig {
    pub daa_file: String,
    pub output: String,
    pub outfmt: Vec<String>,
    pub max_target_seqs: i64,
    pub toppercent: Option<f64>,
    pub forward_only: bool,
    pub report_unaligned: bool,
    pub sam_qlen_field: bool,
    pub invocation: String,
    /// --header value (`simple` / `verbose` / unset).
    pub header: Vec<String>,
    /// Worker count (`0` selects the host's available parallelism, as upstream).
    pub threads: usize,
    /// Unsigned representation of upstream's signed int64 ordinal offsets.
    pub qnum_offset: u64,
    pub snum_offset: u64,
    /// Upstream console verbosity: 0=quiet, 1=normal, 2=verbose, 3=debug log.
    pub verbosity: u8,
}

fn print_daa_metadata(daa: &DaaFile, score_matrix: &ScoreMatrix, verbosity: u8) {
    if verbosity == 0 {
        return;
    }
    // Upstream constructs `ScoreMatrix` from the archive's matrix name, gap
    // penalties, and database size, then prints that derived object.  It does
    // not echo the archived K/Lambda header fields directly.
    eprintln!(
        "Scoring parameters: (Matrix={} Lambda={} K={} Penalties={}/{})",
        daa.score_matrix(),
        score_matrix.lambda(),
        score_matrix.k(),
        score_matrix.gap_open(),
        score_matrix.gap_extend()
    );
    if verbosity >= 2 {
        eprintln!("Build version = {}", daa.diamond_build());
    }
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());
}

/// Print the messages produced when upstream constructs `Search::Config`
/// inside `view_daa`.  This happens after the archive metadata has been
/// printed, so keeping it here (rather than in the CLI wrapper) also preserves
/// the observable ordering of verbose/debug output.
fn print_view_selection(config: &ViewConfig, dp_fields: HspValues) {
    if config.verbosity == 0 {
        return;
    }
    if config.verbosity >= 3 {
        eprintln!("Min length ratio: 0");
    }
    if let Some(top) = config.toppercent {
        eprintln!("Percentage range of top alignment score to report hits: {top}");
    } else if config.max_target_seqs == i64::MAX {
        eprintln!("#Target sequences to report alignments for: unlimited");
    } else {
        eprintln!(
            "#Target sequences to report alignments for: {}",
            config.max_target_seqs
        );
    }
    if config.verbosity >= 3 {
        eprintln!("DP fields: {}", dp_fields.0);
    }
}

fn view_threads(requested: usize) -> usize {
    if requested == 0 {
        std::thread::available_parallelism().map_or(1, usize::from)
    } else {
        requested
    }
}

/// Parse the `--header` flag value (as a Vec<String> coming from CLI parsing).
/// Mirrors `TabularFormat::header_format(Config::blastp)`: missing flag means
/// no header, bare `--header` means verbose, and an explicit value must be one
/// of `0`, `simple`, or `verbose`.
fn parse_header_arg(header: &[String]) -> io::Result<Header> {
    if header.is_empty() {
        return Ok(Header::None);
    }
    if header.len() > 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Invalid header format: {}", header.join(" ")),
        ));
    }
    let value = header[0].trim();
    if value.is_empty() {
        return Ok(Header::Verbose);
    }
    Header::parse(value).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Invalid value for string field: {value}. Permitted values: 0, simple, verbose"
            ),
        )
    })
}

fn reject_unsupported_options(config: &ViewConfig) -> io::Result<()> {
    if config.report_unaligned {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Option is not permitted for this workflow: unal",
        ));
    }
    Ok(())
}

/// Upstream's empty output filename addresses the inherited stdout stream and
/// leaves that process-global stream open. Owning a `Stdout` handle has the
/// same lifetime behavior; the buffer is explicitly flushed by each command.
fn open_output_writer(path: &str) -> io::Result<Box<dyn Write>> {
    if path.is_empty() {
        Ok(Box::new(BufWriter::new(io::stdout())))
    } else {
        Ok(Box::new(File::create(path)?))
    }
}

/// DAA finalization rewrites header 2, so its stdout sink must retain seek
/// support when stdout is redirected to a regular file. This matches C++'s
/// `fdopen(STDOUT_FILENO, "wb")`: a pipe/terminal is accepted initially and
/// reports the seek failure only when finalization attempts `fseek`.
fn open_seekable_output(path: &str) -> io::Result<File> {
    if !path.is_empty() {
        return File::create(path);
    }
    #[cfg(windows)]
    {
        crate::util::io::duplicate_standard_handle(
            windows_sys::Win32::System::Console::STD_OUTPUT_HANDLE,
        )
    }
    #[cfg(not(windows))]
    {
        OpenOptions::new().write(true).open("/dev/stdout")
    }
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

/// Run the pinned `view_worker` boundary: 32 query records per task and at
/// most `3 * threads` ordered tasks resident.  Each worker renders a whole
/// task into one buffer, and the caller writes completed buffers in input
/// order.  The mutable archive reader remains on this thread; workers use a
/// second handle only for immutable header/reference metadata.
fn run_ordered_view_workers<Render, WriteBuffer>(
    daa: &mut DaaFile,
    daa_path: &str,
    requested_threads: usize,
    render: Render,
    mut write_buffer: WriteBuffer,
) -> io::Result<()>
where
    Render: Fn(&DaaQueryRecord, &DaaFile, &mut Vec<u8>) -> Result<(), String> + Sync,
    WriteBuffer: FnMut(&[u8]) -> io::Result<()>,
{
    let threads = view_threads(requested_threads);
    let queue_limit = 3usize.saturating_mul(threads).max(1);
    let worker_daa = DaaFile::open(daa_path)?;
    let (job_tx, job_rx) = mpsc::sync_channel::<(usize, Vec<ViewQueryBuffer>)>(queue_limit);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (result_tx, result_rx) = mpsc::channel::<(usize, Result<Vec<u8>, String>)>();

    std::thread::scope(|scope| -> io::Result<()> {
        for _ in 0..threads {
            let job_rx = Arc::clone(&job_rx);
            let result_tx = result_tx.clone();
            let worker_daa = &worker_daa;
            let render = &render;
            scope.spawn(move || loop {
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
                        let record = DaaQueryRecord::from_buffer(worker_daa, &buf, query_num)?;
                        render(&record, worker_daa, &mut output)?;
                    }
                    Ok(output)
                })();
                if result_tx.send((sequence, rendered)).is_err() {
                    return;
                }
            });
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
/// Non-DAA record rendering still uses FFI via the binary fallback.
pub fn run(daa_file: &str) -> io::Result<()> {
    let daa = DaaFile::open(daa_file)?;

    let mode_str = match daa.mode() {
        2 => "blastp",
        3 => "blastx",
        4 => "blastn",
        _ => "unknown",
    };

    eprintln!("DAA file: {}", daa_file);
    eprintln!("Mode: {}", mode_str);
    eprintln!("Database sequences: {}", daa.db_seqs());
    eprintln!("Database sequences used: {}", daa.db_seqs_used());
    eprintln!("Database letters: {}", daa.db_letters());
    eprintln!("Query records: {}", daa.query_records());
    eprintln!("Matrix: {}", daa.score_matrix());
    eprintln!(
        "Gap penalties: {}/{}",
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );

    if daa.db_seqs_used() > 0 {
        eprintln!(
            "First reference: {} ({})",
            daa.ref_name(0),
            daa.ref_len_at(0)
        );
    }

    Ok(())
}

/// Matches the DAA-output branch of C++ `view_daa`.
pub fn run_daa(config: &ViewConfig) -> io::Result<()> {
    reject_unsupported_options(config)?;
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
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

    print_daa_metadata(&daa, &score_matrix, config.verbosity);
    print_view_selection(config, HspValues::TRANSCRIPT);

    let mut out = open_seekable_output(&config.output)?;
    init_daa(&mut out)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        |record, worker_daa, query_out| {
            view_query_daa(record, worker_daa, query_out, &score_matrix, &view_config)
        },
        |query_out| out.write_all(query_out),
    )?;
    finish_daa_from_file(&mut out, &daa)
}

/// Matches the tabular-output branch of C++ `view_daa`.
pub fn run_tabular(config: &ViewConfig) -> io::Result<()> {
    reject_unsupported_options(config)?;
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let outfmt = if config.outfmt.is_empty() {
        vec!["tab".to_string()]
    } else {
        config.outfmt.clone()
    };
    let is_json = outfmt
        .first()
        .is_some_and(|f| f == "json-flat" || f == "104");
    let outfmt_refs = outfmt.iter().map(String::as_str).collect::<Vec<_>>();
    let format = TabularFormat::new(
        &outfmt_refs,
        is_json,
        0,
        daa.mode() == 3,
        &daa.score_matrix(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
        qnum_offset: config.qnum_offset,
        snum_offset: config.snum_offset,
    };

    print_daa_metadata(&daa, &score_matrix, config.verbosity);
    print_view_selection(config, format.hsp_values);

    let mut out = open_output_writer(&config.output)?;
    let header_mode = parse_header_arg(&config.header)?;
    let header = print_tabular_header(
        &format.fields,
        header_mode,
        format.is_json,
        VERSION_STRING,
        &config.invocation,
    )
    .map_err(io::Error::other)?;
    out.write_all(header.as_bytes())?;
    let threads = view_threads(config.threads);
    if threads == 1 {
        while let Some((buf, query_num)) = daa.read_query_buffer()? {
            let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let mut query_out = Vec::new();
            view_query_tabular(
                &r,
                &daa,
                &mut query_out,
                &score_matrix,
                &format,
                &view_config,
                config.report_unaligned,
            )
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            out.write_all(&query_out)?;
        }
        out.write_all(print_tabular_footer(format.is_json).as_bytes())?;
        out.flush()?;
        return Ok(());
    }

    // Upstream renders through a bounded ordered worker queue. Keep at most
    // one 32-record window resident, but refill slots as soon as their ordered
    // output has been written instead of imposing a barrier between windows.
    // A second DAA handle supplies immutable archive metadata to workers while
    // the original handle advances the alignment-block reader.
    const VIEW_BATCH_SIZE: usize = 32;
    let worker_daa = DaaFile::open(&config.daa_file)?;
    let (job_tx, job_rx) = mpsc::sync_channel::<(usize, Vec<u8>, usize)>(VIEW_BATCH_SIZE);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (result_tx, result_rx) = mpsc::channel();

    std::thread::scope(|scope| -> io::Result<()> {
        for _ in 0..threads {
            let job_rx = Arc::clone(&job_rx);
            let result_tx = result_tx.clone();
            let worker_daa = &worker_daa;
            let score_matrix = &score_matrix;
            let format = &format;
            let view_config = &view_config;
            scope.spawn(move || loop {
                let job = match job_rx.lock() {
                    Ok(rx) => rx.recv(),
                    Err(_) => return,
                };
                let Ok((sequence, buf, query_num)) = job else {
                    return;
                };
                let rendered = (|| {
                    let r = DaaQueryRecord::from_buffer(worker_daa, &buf, query_num)?;
                    let mut query_out = Vec::new();
                    view_query_tabular(
                        &r,
                        worker_daa,
                        &mut query_out,
                        score_matrix,
                        format,
                        view_config,
                        config.report_unaligned,
                    )?;
                    Ok::<_, String>(query_out)
                })();
                if result_tx.send((sequence, rendered)).is_err() {
                    return;
                }
            });
        }
        drop(result_tx);

        let mut pending = BTreeMap::new();
        let mut next_submit = 0usize;
        let mut next_write = 0usize;
        let mut resident = 0usize;
        let mut at_end = false;

        // Always close the job queue before leaving the scope, including on
        // malformed input or output failure. Otherwise scoped workers could
        // remain blocked in `recv` while the scope waits to join them.
        let result = (|| -> io::Result<()> {
            while resident > 0 || !at_end {
                while resident < VIEW_BATCH_SIZE && !at_end {
                    match daa.read_query_buffer()? {
                        Some((buf, query_num)) => {
                            job_tx.send((next_submit, buf, query_num)).map_err(|_| {
                                io::Error::new(
                                    io::ErrorKind::BrokenPipe,
                                    "view worker queue closed",
                                )
                            })?;
                            next_submit += 1;
                            resident += 1;
                        }
                        None => at_end = true,
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
                    let query_out =
                        rendered.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                    // Pinned upstream writes each ordered per-query buffer verbatim.
                    out.write_all(&query_out)?;
                    next_write += 1;
                    resident -= 1;
                }
            }
            Ok(())
        })();
        drop(job_tx);
        result
    })?;
    out.write_all(print_tabular_footer(format.is_json).as_bytes())?;
    out.flush()
}

/// Matches the PAF-output branch of C++ `view_daa`.
pub fn run_paf(config: &ViewConfig) -> io::Result<()> {
    reject_unsupported_options(config)?;
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
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

    print_daa_metadata(&daa, &score_matrix, config.verbosity);
    print_view_selection(config, HspValues::TRANSCRIPT);

    let mut out = open_output_writer(&config.output)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        |record, worker_daa, query_out| {
            view_query_paf(record, worker_daa, query_out, &score_matrix, &view_config)
        },
        |query_out| out.write_all(query_out),
    )?;
    out.flush()
}

/// Matches the SAM-output branch of C++ `view_daa`.
pub fn run_sam(config: &ViewConfig) -> io::Result<()> {
    reject_unsupported_options(config)?;
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
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

    print_daa_metadata(&daa, &score_matrix, config.verbosity);
    print_view_selection(config, HspValues::TRANSCRIPT);

    let mut out = open_output_writer(&config.output)?;
    sam::print_header(&mut out, daa.mode(), &config.invocation)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        |record, worker_daa, query_out| {
            view_query_sam(
                record,
                worker_daa,
                query_out,
                &score_matrix,
                &view_config,
                config.sam_qlen_field,
            )
        },
        |query_out| out.write_all(query_out),
    )?;
    out.flush()
}

/// Matches the pairwise-output branch of C++ `view_daa`.
pub fn run_pairwise(config: &ViewConfig) -> io::Result<()> {
    reject_unsupported_options(config)?;
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
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

    print_daa_metadata(&daa, &score_matrix, config.verbosity);
    print_view_selection(config, HspValues::TRANSCRIPT);

    let mut out = open_output_writer(&config.output)?;
    pairwise::write_header(&mut out)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        |record, worker_daa, query_out| {
            view_query_pairwise(record, worker_daa, query_out, &score_matrix, &view_config)
        },
        |query_out| out.write_all(query_out),
    )?;
    pairwise::print_footer(&mut out)?;
    out.flush()
}

/// Matches the XML-output branch of C++ `view_daa`.
pub fn run_xml(config: &ViewConfig) -> io::Result<()> {
    reject_unsupported_options(config)?;
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
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

    print_daa_metadata(&daa, &score_matrix, config.verbosity);
    print_view_selection(config, HspValues::TRANSCRIPT);

    let mut out = open_output_writer(&config.output)?;
    if let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        xml::print_header(
            &mut out,
            daa.mode(),
            &daa.score_matrix(),
            daa.gap_open_penalty(),
            daa.gap_extension_penalty(),
            daa.evalue(),
            &r.query_name,
            r.query_len() as u32,
            "",
        )?;
        let mut query_out = Vec::new();
        view_query_xml(&r, &daa, &mut query_out, &score_matrix, &view_config)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        out.write_all(&query_out)?;
        run_ordered_view_workers(
            &mut daa,
            &config.daa_file,
            config.threads,
            |record, worker_daa, query_out| {
                view_query_xml(record, worker_daa, query_out, &score_matrix, &view_config)
            },
            |query_out| out.write_all(query_out),
        )?;
    } else {
        xml::print_header(
            &mut out,
            daa.mode(),
            &daa.score_matrix(),
            daa.gap_open_penalty(),
            daa.gap_extension_penalty(),
            daa.evalue(),
            "",
            0,
            "",
        )?;
    }
    xml::print_footer(&mut out)?;
    out.flush()
}

/// Matches the null-output branch of C++ `view_daa`.
pub fn run_null(config: &ViewConfig) -> io::Result<()> {
    reject_unsupported_options(config)?;
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
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

    print_daa_metadata(&daa, &score_matrix, config.verbosity);
    print_view_selection(config, HspValues::TRANSCRIPT);

    let mut out = open_output_writer(&config.output)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        |record, worker_daa, _| view_query_null(record, worker_daa, &score_matrix, &view_config),
        |_| Ok(()),
    )?;
    out.flush()
}

/// Matches the edge-output branch of C++ `view_daa`.
pub fn run_edge(config: &ViewConfig) -> io::Result<()> {
    reject_unsupported_options(config)?;
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
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

    print_daa_metadata(&daa, &score_matrix, config.verbosity);
    print_view_selection(config, HspValues::COORDS);

    let mut out = open_output_writer(&config.output)?;
    run_ordered_view_workers(
        &mut daa,
        &config.daa_file,
        config.threads,
        |record, worker_daa, query_out| {
            view_query_edge(record, worker_daa, query_out, &score_matrix, &view_config)
        },
        |query_out| out.write_all(query_out),
    )?;
    out.flush()
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
        // View command shows DAA metadata — full record parsing uses FFI
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
            daa_file: input.to_string_lossy().to_string(),
            output: output.to_string_lossy().to_string(),
            outfmt: vec!["daa".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: String::new(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        let mut out_daa = DaaFile::open(&output).unwrap();
        assert_eq!(out_daa.query_records(), 1);
        assert_eq!(out_daa.db_seqs_used(), 1);
        assert_eq!(out_daa.ref_name(0), "ref0");
        let (buf, query_num) = out_daa.read_query_buffer().unwrap().unwrap();
        let r = DaaQueryRecord::from_buffer(&out_daa, &buf, query_num).unwrap();
        assert_eq!(r.query_name, "query0");
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
            daa_file: input.to_string_lossy().to_string(),
            output: output.to_string_lossy().to_string(),
            outfmt: vec![
                "tab".to_string(),
                "qseqid".to_string(),
                "sseqid".to_string(),
                "score".to_string(),
                "length".to_string(),
                "pident".to_string(),
                "hspnum".to_string(),
                "qnum".to_string(),
                "snum".to_string(),
            ],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: String::new(),
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
            daa_file: input.to_string_lossy().to_string(),
            output: output_json.to_string_lossy().to_string(),
            outfmt: vec![
                "104".to_string(),
                "qseqid".to_string(),
                "sseqid".to_string(),
                "score".to_string(),
                "qnum".to_string(),
                "snum".to_string(),
            ],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: String::new(),
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
            daa_file: input.to_string_lossy().to_string(),
            output: output_paf.to_string_lossy().to_string(),
            outfmt: vec!["103".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: String::new(),
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
            daa_file: input.to_string_lossy().to_string(),
            output: output_sam.to_string_lossy().to_string(),
            outfmt: vec!["101".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: true,
            header: vec![],
            invocation: "diamond view -f sam".to_string(),
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
            daa_file: input.to_string_lossy().to_string(),
            output: output_pairwise.to_string_lossy().to_string(),
            outfmt: vec!["0".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: String::new(),
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
            daa_file: input.to_string_lossy().to_string(),
            output: output_xml.to_string_lossy().to_string(),
            outfmt: vec!["5".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: String::new(),
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
            daa_file: input.to_string_lossy().to_string(),
            output: output_null.to_string_lossy().to_string(),
            outfmt: vec!["null".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: String::new(),
            threads: 1,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        })
        .unwrap();

        assert!(std::fs::read(&output_null).unwrap().is_empty());

        run_edge(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output_edge.to_string_lossy().to_string(),
            outfmt: vec!["edge".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            header: vec![],
            invocation: String::new(),
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
            daa_file: input.to_string_lossy().to_string(),
            output: output.to_string_lossy().to_string(),
            outfmt: vec![
                "104".to_string(),
                "qseqid".to_string(),
                "sseqid".to_string(),
                "score".to_string(),
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
            daa_file: "does-not-exist.daa".to_string(),
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
    fn test_view_header_modes_match_upstream_tabular_rules() {
        assert_eq!(parse_header_arg(&[]).unwrap(), Header::None);
        assert_eq!(parse_header_arg(&[String::new()]).unwrap(), Header::Verbose);
        assert_eq!(
            parse_header_arg(&["simple".to_string()]).unwrap(),
            Header::Simple
        );
        assert_eq!(parse_header_arg(&["0".to_string()]).unwrap(), Header::None);
        assert_eq!(
            parse_header_arg(&["bad".to_string()])
                .unwrap_err()
                .to_string(),
            "Invalid value for string field: bad. Permitted values: 0, simple, verbose"
        );
        assert_eq!(
            parse_header_arg(&["simple".to_string(), "verbose".to_string()])
                .unwrap_err()
                .to_string(),
            "Invalid header format: simple verbose"
        );
    }
}
