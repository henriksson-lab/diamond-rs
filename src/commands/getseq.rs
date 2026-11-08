use std::cell::RefCell;
use std::collections::BTreeSet;
use std::ffi::OsString;
use std::io;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use crate::basic::value::Letter;
use crate::basic::value::AMINO_ACID_ALPHABET;
use crate::data::blastdb::blastdb::BlastDB;
use crate::data::dmnd::{DatabaseFile, DatabaseOpenOptions, SourceDatabaseFile};
use crate::data::sequence_file::{
    detect_type_source_with_buffer, SequenceFileConfig, SequenceFileFlags, SequenceFileType,
};
use crate::util::command_line_parser::CommandLineParser;
use crate::util::io::{Compressor, OutputFile};
use crate::util::log_stream::MessageStream;
use crate::util::options::OptionValue;

/// Baseline command-line protocol. Keep this separate from the typed API's
/// comma-list convenience: source --seq is a vector of individual atoi inputs.
pub(crate) fn run_cli(args: &[String]) -> io::Result<()> {
    run_os(&args.iter().map(OsString::from).collect::<Vec<_>>())
}
pub(crate) fn run_cli_with_startup(
    args: &[String],
    startup: Option<std::sync::Arc<super::source_config::Startup>>,
) -> io::Result<()> {
    run_with_startup_os(
        &args.iter().map(OsString::from).collect::<Vec<_>>(),
        startup,
    )
}
pub(crate) fn run_os(args: &[OsString]) -> io::Result<()> {
    run_with_startup_os(args, None)
}
pub(crate) fn run_with_startup_os(
    args: &[OsString],
    startup: Option<std::sync::Arc<super::source_config::Startup>>,
) -> io::Result<()> {
    let startup = match startup {
        Some(startup) => startup,
        None => super::source_config::prepare_os(args)?,
    };
    #[cfg(feature = "extra")]
    if args.len() > 2 {
        return run_extra_os(args, &startup);
    }
    run_baseline_os(args, &startup)
}
fn run_baseline_os(args: &[OsString], startup: &super::source_config::Startup) -> io::Result<()> {
    let mut parser = CommandLineParser::new();
    parser.add_command(
        "getseq",
        "Retrieve sequences from a DIAMOND database file",
        1,
    );
    let general = parser.add_group("General options", vec![1], false);
    let db_group = parser.add_group("General/database options", vec![1], false);
    let out_group = parser.add_group("General/output", vec![1], false);
    let seq_group = parser.add_group("Getseq options", vec![1], false);
    let forbidden = parser.add_group("Other workflows", vec![], false);
    let disabled = parser.add_group("", vec![], true);
    let database = Rc::new(RefCell::new(OptionValue::<PathBuf>::new()));
    let output = Rc::new(RefCell::new(PathBuf::new()));
    let numbers = Rc::new(RefCell::new(Vec::<Vec<u8>>::new()));
    let threads = Rc::new(RefCell::new(0i32));
    let verbose = Rc::new(RefCell::new(false));
    let log = Rc::new(RefCell::new(false));
    let quiet = Rc::new(RefCell::new(false));
    let tmpdir = Rc::new(RefCell::new(PathBuf::new()));
    let rejected = Rc::new(RefCell::new(false));
    let mut hidden = false;
    // REGISTERED_OPTIONS preserves constructor order. The hidden group starts
    // at match1; DNA-only names are absent from the baseline. Register rejected
    // fields too so permission errors and short-ID shadowing stay source-exact.
    for id in crate::basic::config::registered_options(false) {
        if id == "match1" {
            hidden = true;
        }
        let short = match id {
            "threads" => 'p',
            "verbose" => 'v',
            "tmpdir" => 't',
            "db" => 'd',
            "out" => 'o',
            "evalue" => 'e',
            "max-target-seqs" => 'k',
            "shapes" => 's',
            "query" => 'q',
            "global-ranking" => 'g',
            "block-size" => 'b',
            "index-chunks" => 'c',
            "frameshift" => 'F',
            "outfmt" => 'f',
            "memory-limit" => 'M',
            "min-orf" => 'l',
            "xdrop" => 'x',
            "daa" => 'a',
            "window" | "seed-weight" => 'w',
            "gapped-xdrop" => 'X',
            _ => '\0',
        };
        match id {
            "threads" => {
                parser.add_option(
                    general,
                    id,
                    short,
                    "number of CPU threads",
                    threads.clone(),
                    0i32,
                    1,
                );
            }
            "verbose" => {
                parser.add_option(
                    general,
                    id,
                    short,
                    "verbose console output",
                    verbose.clone(),
                    false,
                    1,
                );
            }
            "log" => {
                parser.add_option(
                    general,
                    id,
                    short,
                    "enable debug log",
                    log.clone(),
                    false,
                    1,
                );
            }
            "quiet" => {
                parser.add_option(
                    general,
                    id,
                    short,
                    "disable console output",
                    quiet.clone(),
                    false,
                    1,
                );
            }
            "tmpdir" => {
                parser.add_option(
                    general,
                    id,
                    short,
                    "directory for temporary files",
                    tmpdir.clone(),
                    PathBuf::new(),
                    1,
                );
            }
            "db" => {
                parser.add_option(
                    db_group,
                    id,
                    short,
                    "database file",
                    database.clone(),
                    OptionValue::<PathBuf>::new(),
                    1,
                );
            }
            "out" => {
                parser.add_option(
                    out_group,
                    id,
                    short,
                    "output file",
                    output.clone(),
                    PathBuf::new(),
                    1,
                );
            }
            "seq" => {
                parser.add_option(
                    seq_group,
                    id,
                    short,
                    "Space-separated list of sequence numbers to display.",
                    numbers.clone(),
                    Vec::<Vec<u8>>::new(),
                    1,
                );
            }
            _ => {
                parser.add_option(
                    if hidden { disabled } else { forbidden },
                    id,
                    short,
                    "",
                    rejected.clone(),
                    false,
                    1,
                );
            }
        }
    }
    let mut command = 0;
    parser.store_os(
        &args.iter().map(OsString::as_os_str).collect::<Vec<_>>(),
        &mut command,
    )?;
    let mut messages = crate::util::log_stream::message_stream()
        .lock()
        .unwrap()
        .clone();
    database.borrow().require().map_err(argument_error)?;
    let _ = (verbose, log, quiet);
    let database = startup.config.database.clone();
    let output = output.borrow();
    let selected = numbers.borrow();
    run_selected_path(
        &database,
        &selected,
        (!output.as_os_str().is_empty()).then_some(output.as_path()),
        Some(&mut messages),
        &GetSeqRuntime {
            source_cursor: true,
            ..Default::default()
        },
    )
}

#[derive(Default)]
struct GetSeqRuntime {
    source_cursor: bool,
    file_buffer_size: Option<usize>,
    titles_file: PathBuf,
    oid_list: PathBuf,
    reverse: bool,
    hardmasked: bool,
    max_letters: Option<usize>,
}

#[cfg(feature = "extra")]
fn run_extra_os(args: &[OsString], startup: &super::source_config::Startup) -> io::Result<()> {
    use crate::util::command_line_parser::ReadOption;
    let values = super::argument_registry::values_os(args)?;
    let scalar = |name: &str| -> PathBuf {
        values
            .get(name)
            .and_then(|items| items.first())
            .cloned()
            .unwrap_or_default()
            .into()
    };
    let mut size = 67_108_864i64;
    if let Some(items) = values.get("file-buffer-size") {
        i64::read_option_os(&mut size, items)?;
    }
    let allocation_error = || {
        io::Error::other("Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage.")
    };
    let size = usize::try_from(size).map_err(|_| allocation_error())?;
    Vec::<u8>::new()
        .try_reserve_exact(size)
        .map_err(|_| allocation_error())?;
    let mut chunk = 0.0f64;
    if let Some(items) = values.get("block-size") {
        f64::read_option_os(&mut chunk, items)?;
    }
    let mut messages = crate::util::log_stream::message_stream()
        .lock()
        .unwrap()
        .clone();
    if startup.config.database.as_os_str().is_empty() {
        return Err(argument_error("Missing parameter: --db/-d".into()));
    }
    let titles = values
        .get("query")
        .and_then(|files| files.first())
        .cloned()
        .unwrap_or_default();
    let runtime = GetSeqRuntime {
        source_cursor: true,
        file_buffer_size: Some(size),
        titles_file: titles.into(),
        oid_list: scalar("oid-list"),
        reverse: values.contains_key("reverse"),
        hardmasked: values.contains_key("hard-masked"),
        max_letters: (chunk != 0.0).then_some((chunk * 1e9) as usize),
    };
    let numbers = values
        .get("seq")
        .map(Vec::as_slice)
        .unwrap_or(&[])
        .iter()
        .map(|v| {
            crate::util::io::source_os_bytes(v)
                .map(|v| v.into_owned())
                .map_err(io::Error::other)
        })
        .collect::<io::Result<Vec<_>>>()?;
    run_selected_path(
        &startup.config.database,
        &numbers,
        (!startup.config.output_file.as_os_str().is_empty())
            .then_some(startup.config.output_file.as_path()),
        Some(&mut messages),
        &runtime,
    )
}

fn argument_error(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message)
}

pub(crate) enum Source {
    Dmnd(Box<DatabaseFile>, u64),
    SourceDmnd(Box<SourceDatabaseFile>, u64),
    Fasta(Box<super::search_blocks::RecordStream>),
    Blast(Box<BlastDB>),
}
impl Source {
    pub(crate) fn open(database: &str) -> io::Result<Self> {
        Self::open_path(Path::new(database))
    }
    pub(crate) fn open_path(database: &Path) -> io::Result<Self> {
        Self::open_path_with_buffer(database, None)
    }
    pub(crate) fn open_with_buffer(database: &str, buffer: Option<usize>) -> io::Result<Self> {
        Self::open_path_with_buffer(Path::new(database), buffer)
    }
    pub(crate) fn open_path_with_buffer(
        database: &Path,
        buffer: Option<usize>,
    ) -> io::Result<Self> {
        Self::open_with_cursor_path(database, buffer, false)
    }
    fn open_with_cursor_path(
        database: &Path,
        buffer: Option<usize>,
        source_cursor: bool,
    ) -> io::Result<Self> {
        if let Some(size) = buffer {
            super::blastp::preflight_file_buffer_size(size)?;
        }
        let (kind, path) = detect_type_source_with_buffer(
            &[database.to_owned()],
            SequenceFileFlags::NONE,
            SequenceFileConfig::default(),
            buffer.unwrap_or(67_108_864),
        )?;
        match kind {
            SequenceFileType::Dmnd if source_cursor => Ok(Self::SourceDmnd(
                Box::new(SourceDatabaseFile::open(
                    &path,
                    DatabaseOpenOptions {
                        file_buffer_size: buffer,
                        ..Default::default()
                    },
                )?),
                0,
            )),
            SequenceFileType::Dmnd => Ok(Self::Dmnd(
                Box::new(
                    DatabaseFile::open_with_options_io(
                        path,
                        DatabaseOpenOptions {
                            file_buffer_size: buffer,
                            ..Default::default()
                        },
                    )
                    .map_err(io::Error::other)?,
                ),
                0,
            )),
            // Upstream accepts FASTA here, but its NONE-flags constructor can
            // leave the sequence-count member uninitialized. Iterate the
            // concrete reader to EOF instead of reproducing that undefined bound.
            SequenceFileType::Fasta => Ok(Self::Fasta(Box::new(
                super::search_blocks::RecordStream::open_ready(
                    path,
                    crate::basic::value::SequenceType::AminoAcid,
                    false,
                    buffer.unwrap_or(67_108_864),
                )?,
            ))),
            SequenceFileType::Blast => Ok(Self::Blast(Box::new(BlastDB::new_path(
                &path,
                SequenceFileFlags::NONE,
            )?))),
            SequenceFileType::Block => {
                unreachable!("auto_create does not infer an in-memory block")
            }
        }
    }
    /// Checked legacy String access; the actual command consumes raw bytes.
    pub(crate) fn read(&mut self) -> io::Result<Option<(Vec<Letter>, String)>> {
        self.read_bytes()?
            .map(|(sequence, title)| {
                String::from_utf8(title)
                    .map(|title| (sequence, title))
                    .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
            })
            .transpose()
    }
    pub(crate) fn read_bytes(&mut self) -> io::Result<Option<(Vec<Letter>, Vec<u8>)>> {
        match self {
            Self::SourceDmnd(database, cursor) => {
                if *cursor == database.ref_header.sequences {
                    return Ok(None);
                }
                let record = database.read_seq_bytes()?;
                *cursor += 1;
                Ok(Some(record))
            }
            Self::Dmnd(database, cursor) => {
                if *cursor == database.sequence_count() {
                    return Ok(None);
                }
                let record = database.record_bytes_io(*cursor)?;
                *cursor += 1;
                Ok(Some(record))
            }
            Self::Fasta(database) => {
                Ok(database.next()?.map(|record| (record.sequence, record.id)))
            }
            Self::Blast(database) => {
                if database.eof() {
                    Ok(None)
                } else {
                    database.source_io(|db| db.read_seq_bytes()).map(Some)
                }
            }
        }
    }
}

/// Run the getseq command — retrieve sequences from a DIAMOND database.
pub fn run(database: &str, seq_numbers: Option<&str>) -> io::Result<()> {
    run_with_output(database, seq_numbers, None)
}

/// Retrieve sequences and optionally write them to `-o/--out`. Keeping stdout
/// as the `None` case preserves the library API while matching the CLI's
/// upstream output-file behavior.
pub fn run_with_output(
    database: &str,
    seq_numbers: Option<&str>,
    output: Option<&Path>,
) -> io::Result<()> {
    let seq_numbers: Vec<String> = seq_numbers
        .map(|numbers| {
            numbers
                .split(|c: char| c == ',' || c.is_ascii_whitespace())
                .filter(|number| !number.is_empty())
                .map(ToOwned::to_owned)
                .collect()
        })
        .unwrap_or_default();
    run_with_output_numbers(database, &seq_numbers, output)
}

/// Retrieve sequences by DIAMOND's 1-based numeric sequence numbers.
///
/// The CLI accepts a space-separated list after `--seq`. The compatibility
/// wrapper above also accepts comma-separated numbers for existing Rust API
/// callers.
pub fn run_with_output_numbers(
    database: &str,
    seq_numbers: &[String],
    output: Option<&Path>,
) -> io::Result<()> {
    run_numbers(database, seq_numbers, output, None)
}

pub(crate) fn run_path_with_output_numbers(
    database: &Path,
    numbers: &[Vec<u8>],
    output: Option<&Path>,
) -> io::Result<()> {
    run_selected_path(database, numbers, output, None, &GetSeqRuntime::default())
}
fn run_numbers(
    database: &str,
    seq_numbers: &[String],
    output: Option<&Path>,
    messages: Option<&mut MessageStream>,
) -> io::Result<()> {
    let numbers = seq_numbers
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect::<Vec<_>>();
    run_selected_path(
        Path::new(database),
        &numbers,
        output,
        messages,
        &GetSeqRuntime::default(),
    )
}
fn run_selected_path(
    database: &Path,
    seq_numbers: &[Vec<u8>],
    output: Option<&Path>,
    mut messages: Option<&mut MessageStream>,
    runtime: &GetSeqRuntime,
) -> io::Result<()> {
    // Source opens the database before reading either selection list and
    // before opening/truncating the output.
    let mut source =
        Source::open_with_cursor_path(database, runtime.file_buffer_size, runtime.source_cursor)?;
    let mut titles = std::collections::BTreeMap::new();
    if !runtime.titles_file.as_os_str().is_empty() {
        let mut input = crate::util::io::text_input_file::ByteTextInputFile::from_source_path_bytes_with_buffer_size(
            &runtime.titles_file,
            runtime.file_buffer_size.unwrap_or(67_108_864),
        )
        .map_err(io::Error::other)?;
        loop {
            input.getline().map_err(io::Error::other)?;
            // This particular Source loop ignores an unterminated final line.
            if input.eof() {
                break;
            }
            let fields = crate::util::string::tokenize_bytes(&input.line, b"\t");
            if fields.len() != 2 {
                return Err(io::Error::other("Query file format error."));
            }
            titles.insert(fields[0].clone(), fields[1].clone());
        }
        input.into_inner().close().map_err(io::Error::other)?;
    }
    let all =
        seq_numbers.is_empty() && titles.is_empty() && runtime.oid_list.as_os_str().is_empty();
    let mut requested = seq_numbers
        .iter()
        .map(|number| c_atoi(number).wrapping_sub(1) as usize)
        .collect::<BTreeSet<_>>();
    if !runtime.oid_list.as_os_str().is_empty() {
        let mut input = crate::util::io::text_input_file::ByteTextInputFile::from_source_path_bytes_with_buffer_size(
            &runtime.oid_list,
            runtime.file_buffer_size.unwrap_or(67_108_864),
        )
        .map_err(io::Error::other)?;
        loop {
            input.getline().map_err(io::Error::other)?;
            if input.eof() && input.line.is_empty() {
                break;
            }
            let mut tokenizer = crate::util::string::ByteTokenizer::new(&input.line, b'\t');
            requested.insert(tokenizer.read_u64().map_err(io::Error::other)? as usize);
        }
        input.into_inner().close().map_err(io::Error::other)?;
    }
    if !requested.is_empty() {
        if let Some(messages) = messages.as_mut() {
            messages.write(format!("#Selected sequences: {}\n", requested.len()))?;
        }
    }
    let output = output.unwrap_or(Path::new(""));
    let mut writer = OutputFile::with_source_profile_path(
        output,
        Compressor::None,
        "wb",
        runtime.file_buffer_size.unwrap_or(67_108_864),
    )
    .map_err(io::Error::other)?;
    let mut record = crate::util::text_buffer::TextBuffer::new();
    let value_traits = crate::basic::value::ValueTraits::new(
        AMINO_ACID_ALPHABET,
        crate::basic::value::MASK_LETTER,
        b"UO-",
        crate::basic::value::SequenceType::AminoAcid,
    );
    let mut ordinal = 0usize;
    let mut letters = 0usize;
    while let Some((sequence, id)) = source.read_bytes()? {
        let mapped = titles.get(&crate::util::sequence::seqid_bytes(&id));
        if all || requested.contains(&ordinal) || mapped.is_some() {
            if runtime.reverse || runtime.hardmasked {
                record.clear();
                record.push_back(b'>');
                record.write_raw(mapped.map_or(id.as_slice(), Vec::as_slice));
                record.push_back(b'\n');
                if runtime.reverse {
                    for &letter in sequence.iter().rev() {
                        record.push_back(
                            AMINO_ACID_ALPHABET
                                .get((letter as u8 & 127) as usize)
                                .copied()
                                .unwrap_or(b'X'),
                        );
                    }
                } else {
                    for &letter in &sequence {
                        record.push_back(if letter as u8 & 128 != 0 {
                            b'X'
                        } else {
                            AMINO_ACID_ALPHABET
                                .get(letter as usize)
                                .copied()
                                .unwrap_or(b'X')
                        });
                    }
                }
                record.push_back(b'\n');
            } else {
                // Source ordinary formatting uses the original title even
                // when the selection map supplied a replacement title.
                record.clear();
                crate::util::sequence::format_bytes(
                    &sequence,
                    &id,
                    None,
                    &mut record,
                    "fasta",
                    &value_traits,
                    80,
                )
                .map_err(io::Error::other)?;
            }
            writer.write_raw(record.data()).map_err(io::Error::other)?;
        }
        letters = letters.wrapping_add(sequence.len());
        if runtime.max_letters.is_some_and(|limit| letters >= limit) {
            break;
        }
        ordinal += 1;
    }
    writer.close().map_err(io::Error::other)
}

// Actual selected CRT atoi semantics, including C-locale whitespace/prefixes.
fn c_atoi(value: impl AsRef<[u8]>) -> i32 {
    let mut value = value.as_ref().to_vec();
    value.push(0);
    unsafe { libc::atoi(value.as_ptr().cast()) }
}

#[cfg(test)]
fn format_fasta_record(record: &mut Vec<u8>, id: &str, sequence: &[Letter]) {
    record.clear();
    let line_count = sequence.len().div_ceil(80);
    record.reserve(1 + id.len() + 1 + sequence.len() + line_count);
    record.push(b'>');
    record.extend_from_slice(id.as_bytes());
    record.push(b'\n');

    // C++ getseq writes 80-char lines (sequence_file.cpp:get_seq).
    for chunk in sequence.chunks(80) {
        for &letter in chunk {
            let idx = (letter & 0x1F) as usize;
            record.push(AMINO_ACID_ALPHABET.get(idx).copied().unwrap_or(b'X'));
        }
        record.push(b'\n');
    }
}

#[cfg(test)]
mod tests {
    use super::{c_atoi, format_fasta_record};

    #[test]
    fn atoi_matches_upstream_sequence_number_inputs() {
        assert_eq!(c_atoi("1"), 1);
        assert_eq!(c_atoi("  +3tail"), 3);
        assert_eq!(c_atoi("-2"), -2);
        assert_eq!(c_atoi("not-a-number"), 0);
        assert_eq!(c_atoi(b" +3\xff"), 3);
        assert_eq!(c_atoi(b"\xc2\xa03"), 0);
    }

    #[test]
    fn reusable_record_buffer_formats_wrapping_mask_bits_and_unknown_letters() {
        let mut sequence = vec![0i8; 80];
        sequence.extend_from_slice(&[-127, 31]);
        let mut record = b"stale bytes".to_vec();
        format_fasta_record(&mut record, "seq", &sequence);

        let mut expected = b">seq\n".to_vec();
        expected.extend_from_slice(&[b'A'; 80]);
        expected.extend_from_slice(b"\nRX\n");
        assert_eq!(record, expected);
    }
}
