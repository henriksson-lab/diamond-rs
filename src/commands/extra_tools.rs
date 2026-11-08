//! Concrete optional tools from `tools/tools.cpp` and `tools/tsv.cpp`.
//! SequenceFile consumers share the native generic reader. Hashseqs supplies
//! the commented-out Source record reader with finite, validated EOF handling.

use std::collections::BTreeMap;
use std::ffi::OsString;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use crate::basic::reduction::Reduction;
use crate::basic::shape_config::ShapeConfig;
use crate::basic::value::{SequenceType, ValueTraits, AMINO_ACID_ALPHABET, MASK_LETTER};
use crate::data::enum_seeds::{enum_seeds, EnumSeedsCallback, EnumSeedsContext};
use crate::data::flags::{EnumCfg, SeedEncoding, NO_FILTER};
use crate::data::sequence_set::SequenceSet;
use crate::masking::tantan::TantanMasker;
use crate::masking::MaskingAlgo;
use crate::util::command_line_parser::ReadOption;
use crate::util::io::{Compressor, OutputFile};
use crate::util::sequence::{self, seqid_bytes};
use crate::util::text_buffer::TextBuffer;

type Values = BTreeMap<&'static str, Vec<OsString>>;

fn scalar<T: ReadOption>(values: &Values, name: &str, mut default: T) -> io::Result<T> {
    if let Some(items) = values.get(name) {
        T::read_option_os(&mut default, items)?;
    }
    Ok(default)
}

fn input(values: &Values) -> io::Result<&Path> {
    match values.get("query") {
        Some(paths) if paths.len() > 1 => {
            Err(io::Error::other("--query/-q has more than one argument."))
        }
        Some(paths) => Ok(paths.first().map(|s| Path::new(s)).unwrap_or(Path::new(""))),
        None => Ok(Path::new("")),
    }
}

fn file_buffer(values: &Values) -> io::Result<usize> {
    let size = scalar(values, "file-buffer-size", 67_108_864usize)?;
    super::blastp::preflight_file_buffer_size(size)?;
    Ok(size)
}

fn open_output(path: &Path, compressor: Compressor, buffer: usize) -> io::Result<OutputFile> {
    OutputFile::with_source_profile_path(path, compressor, "wb", buffer).map_err(io::Error::other)
}

/// Source TSV's stdout OutputFile is opened after its input and verbose
/// size diagnostic. Keep that descriptor owner separate from Rust's buffered
/// stdout, whose pending reset is discarded when the standalone owner closes.
struct CutOutput<'a>(&'a std::cell::RefCell<Option<OutputFile>>);
impl Write for CutOutput<'_> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let mut owner = self.0.borrow_mut();
        let file = owner
            .as_mut()
            .ok_or_else(|| io::Error::other("TSV output is not initialized"))?;
        file.write_raw(bytes).map_err(io::Error::other)?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        if let Some(file) = self.0.borrow_mut().as_mut() {
            file.flush().map_err(io::Error::other)?;
        }
        Ok(())
    }
}

/// Handle only concrete tools; remaining registered EXTRA modes retain their
/// untranslated status. Source wc/cut read a two-column TSV, not sequence data.
pub(super) fn run(args: &[String], startup: &super::source_config::Startup) -> io::Result<bool> {
    run_os(
        &args.iter().map(OsString::from).collect::<Vec<_>>(),
        startup,
    )
}
pub(super) fn run_os(
    args: &[OsString],
    startup: &super::source_config::Startup,
) -> io::Result<bool> {
    let command = args
        .get(1)
        .and_then(|s| s.to_str())
        .map(|s| s.strip_prefix("--").unwrap_or(s));
    if !matches!(
        command,
        Some("split" | "hashseqs" | "listseeds" | "wc" | "cut" | "roc-id" | "find-shapes")
    ) {
        return Ok(false);
    }
    let values = super::argument_registry::extra_values_os(args)?;
    if command == Some("roc-id") {
        return Err(io::Error::other("Deprecated command: rocid"));
    }
    let path = input(&values)?;
    match command.unwrap() {
        "find-shapes" => crate::tools::find_shapes::find_shapes_file_path_with_buffer_size(
            path,
            &mut crate::util::io::console::stdout_endl(),
            file_buffer(&values)?,
        )?,
        "split" => split(path, &values)?,
        "hashseqs" => hash_sequences(path, file_buffer(&values)?)?,
        "listseeds" => list_seeds(&values, startup)?,
        "wc" | "cut" => {
            // Source front() on an absent --query is undefined. Report the
            // missing TSV path rather than interpreting stdin as a TSV file.
            if path.as_os_str().is_empty() {
                return Err(io::Error::other("Missing required parameter: query"));
            }
            let size = scalar(&values, "tsv-read-size", 1i64 << 30)?;
            let threads = scalar(&values, "threads", 0i32)?;
            let threads = if threads == 0 {
                crate::util::system::hardware_concurrency() as usize
            } else {
                threads.max(0) as usize
            };
            let quiet = values.contains_key("quiet");
            let log = values.contains_key("log");
            let mut verbose = crate::util::log_stream::MessageStream::new(
                !quiet && values.contains_key("verbose"),
                log,
            );
            let mut opened = |file_size| {
                verbose
                    .write(format!("File size: {file_size} bytes\n"))
                    .map(|_| ())
                    .map_err(io::Error::other)
            };
            let (_file_size, count) = if command == Some("wc") {
                let result = crate::tools::tsv::word_count_with_source_open_callback_io_typed(
                    path,
                    size,
                    threads,
                    file_buffer(&values)?,
                    &mut opened,
                )
                .map_err(io::Error::other)?;
                (result.file_size, Some(result.records))
            } else {
                let buffer_size = file_buffer(&values)?;
                let output = std::cell::RefCell::new(None);
                let mut writer = CutOutput(&output);
                let mut cut_opened = |file_size| {
                    opened(file_size)?;
                    *output.borrow_mut() =
                        Some(open_output(Path::new(""), Compressor::None, buffer_size)?);
                    Ok(())
                };
                let result = crate::tools::tsv::cut_with_source_open_callback_io_typed(
                    path,
                    size,
                    threads,
                    &mut writer,
                    buffer_size,
                    &mut cut_opened,
                )
                .map_err(io::Error::other)?;
                output
                    .into_inner()
                    .ok_or_else(|| io::Error::other("TSV output is not initialized"))?
                    .close()
                    .map_err(io::Error::other)?;
                (result.file_size, None)
            };
            if let Some(count) = count {
                crate::util::log_stream::MessageStream::new(!quiet, log)
                    .write(format!("{count}\n"))?;
            }
        }
        _ => unreachable!(),
    }
    Ok(true)
}

fn split(path: &Path, values: &Values) -> io::Result<()> {
    let buffer_size = file_buffer(values)?;
    let mut source = super::getseq::Source::open_path_with_buffer(path, Some(buffer_size))?;
    let size = scalar(values, "block-size", 0.0f64)?;
    let letters_per_file = (size * 1e9) as usize;
    let traits = ValueTraits::new(
        AMINO_ACID_ALPHABET,
        MASK_LETTER,
        b"UO-",
        SequenceType::AminoAcid,
    );
    let mut number = 0usize;
    // Source creates 0.faa.zst even for empty input and creates an empty first
    // file before the first sequence when the configured letter budget is 0.
    let mut output = open_output(Path::new("0.faa.zst"), Compressor::Zstd, buffer_size)?;
    let mut letters = 0usize;
    let mut processed = 0usize;
    let mut buffer = TextBuffer::new();
    let mut emit = |sequence: &[crate::basic::value::Letter], id: &[u8]| -> io::Result<()> {
        if letters >= letters_per_file {
            output.close().map_err(io::Error::other)?;
            number = number
                .checked_add(1)
                .ok_or_else(|| io::Error::other("Split file count overflow"))?;
            output = open_output(
                Path::new(&format!("{number}.faa.zst")),
                Compressor::Zstd,
                buffer_size,
            )?;
            letters = 0;
        }
        buffer.clear();
        sequence::format_bytes(
            sequence,
            &seqid_bytes(id),
            None,
            &mut buffer,
            "fasta",
            &traits,
            160,
        )
        .map_err(io::Error::other)?;
        output.write_raw(buffer.data()).map_err(io::Error::other)?;
        letters = letters
            .checked_add(sequence.len())
            .ok_or_else(|| io::Error::other("Split letter count overflow"))?;
        processed += 1;
        if processed % 1_000_000 == 0 {
            writeln!(
                crate::util::io::console::stderr(),
                "#Sequences processed: {processed} #letters:{letters}"
            )?;
        }
        Ok(())
    };
    while let Some((sequence, id)) = source.read_bytes()? {
        emit(&sequence, &id)?;
    }
    output.close().map_err(io::Error::other)
}

fn hash_sequences(path: &Path, buffer_size: usize) -> io::Result<()> {
    let mut out = crate::util::io::console::stdout_endl();
    let emit = |record: crate::data::fasta::ByteFastaRecord| -> io::Result<()> {
        let records = [crate::tools::tools::SequenceRecord {
            title: record.id,
            sequence: record.sequence,
        }];
        for line in crate::tools::tools::hash_seqs_bytes(&records) {
            out.write_all(&line)?;
            writeln!(out)?;
        }
        Ok(())
    };
    crate::data::fasta::for_each_fasta_file_with_buffer_size_bytes(
        path,
        SequenceType::AminoAcid,
        buffer_size,
        emit,
    )?;
    Ok(())
}

fn list_seeds(values: &Values, startup: &super::source_config::Startup) -> io::Result<()> {
    let database = scalar(values, "db", PathBuf::new())?;
    let mut source =
        super::getseq::Source::open_path_with_buffer(&database, Some(file_buffer(values)?))?;
    let threshold = scalar(values, "tantan-minMaskProb", 0.9f64)? as f32;
    let masker = TantanMasker::from_score_matrix(
        startup
            .matrix
            .as_ref()
            .ok_or_else(|| io::Error::other("Scoring matrix was not initialized"))?,
        threshold,
    );
    let mut sequences = SequenceSet::new();
    while let Some((mut sequence, _)) = source.read_bytes()? {
        masker.mask_hard(&mut sequence);
        sequences.push(&sequence);
    }
    let reduction = Reduction::new(
        "A R N D C Q E G H I L K M F P S T W Y V",
        AMINO_ACID_ALPHABET,
    );
    let codes = if let Some(codes) = values.get("shape-mask") {
        codes
            .iter()
            .map(|s| {
                crate::util::io::source_os_bytes(s)
                    .map(|s| s.into_owned())
                    .map_err(io::Error::other)
            })
            .collect::<io::Result<Vec<_>>>()?
    } else {
        crate::search::sensitivity::get_shape_codes(crate::config::Sensitivity::Default)
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect()
    };
    let shapes = ShapeConfig::from_codes_bytes(&codes, scalar(values, "shapes", 0u32)?, &reduction)
        .map_err(io::Error::other)?;
    if shapes.count() == 0 {
        return Err(io::Error::other("No seed shape configured."));
    }
    let cfg = EnumCfg {
        partition: None,
        shape_begin: 0,
        shape_end: 1,
        code: SeedEncoding::SpacedFactor,
        skip: None,
        filter_masked_seeds: false,
        mask_seeds: false,
        seed_cut: scalar(values, "seed-cut", 0.0f64)?,
        soft_masking: MaskingAlgo::None,
        minimizer_window: 0,
        filter_low_complexity_seeds: false,
        mask_low_complexity_seeds: false,
        sketch_size: 0,
    };
    let context = EnumSeedsContext {
        shapes: &shapes,
        reduction: &reduction,
        min_query_len: scalar(values, "min-query-len", 0i32)?,
        query_contexts: 1,
    };
    struct Collector(Vec<u64>);
    impl EnumSeedsCallback for Collector {
        fn call(&mut self, seed: u64, _: usize, _: usize, _: i32) -> bool {
            self.0.push(seed);
            true
        }
    }
    let mut collector = Collector(Vec::new());
    let end = u32::try_from(sequences.len()).map_err(io::Error::other)?;
    enum_seeds(
        &mut sequences,
        &mut collector,
        0,
        end,
        &NO_FILTER,
        &cfg,
        &context,
    );
    let top = scalar(values, "query-count", 1usize)?;
    let mut output = crate::util::io::console::stdout_endl();
    for (count, seed) in crate::tools::tools::list_seeds(
        &collector.0,
        &reduction,
        shapes.get(0).weight as usize,
        top,
    ) {
        writeln!(output, "{count}\t{seed}")?;
    }
    Ok(())
}

/// These commands are registered by EXTRA but have no switch arm in the
/// pinned run/main.cpp. Their native command result is the Source default 1;
/// calling an unrelated translated algorithm would invent command behavior.
pub(super) fn source_fallthrough(args: &[String]) -> io::Result<bool> {
    source_fallthrough_os(&args.iter().map(OsString::from).collect::<Vec<_>>())
}
pub(super) fn source_fallthrough_os(args: &[OsString]) -> io::Result<bool> {
    let command = args
        .get(1)
        .and_then(|s| s.to_str())
        .map(|s| s.strip_prefix("--").unwrap_or(s))
        .unwrap_or("");
    if !matches!(
        command,
        "sort"
            | "dbstat"
            | "seed-stat"
            | "simulate-seqs"
            | "upgma"
            | "upgmamc"
            | "compute-medoids"
            | "mutate"
            | "model-seqs"
            | "make-seed-table"
    ) {
        return Ok(false);
    }
    // Source startup was already executed by the shared optional adapter.
    Ok(true)
}
