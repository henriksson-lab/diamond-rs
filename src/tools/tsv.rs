//! Small TSV command helpers.
//!
//! This mirrors `diamond/src/tools/tsv.cpp`. The C++ entry points read their
//! input, chunk size, thread count, and output stream from global state. The
//! Rust API makes those dependencies explicit and returns the values that the
//! command layer can use for verbose and message output.

use std::io::Write;
use std::path::Path;

use crate::util::tsv::{streaming::StreamingFile, Type};
use rayon::prelude::*;

/// Result of [`word_count`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WordCount {
    /// Size of the input file in bytes, matching the C++ verbose message.
    pub file_size: i64,
    /// Number of parsed TSV records, matching the C++ message output.
    pub records: i64,
}

/// Result of [`cut`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CutStats {
    /// Size of the input file in bytes, matching the C++ verbose message.
    pub file_size: i64,
    /// Number of projected records written.
    pub records: i64,
}

/// Count two-column TSV records in chunked reads.
///
/// `read_size` and `threads` correspond to `config.tsv_read_size` and
/// `config.threads_`. As in the source, `threads` is forwarded to the TSV
/// reader and a final empty chunk terminates the loop.
pub fn word_count(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
) -> Result<WordCount, String> {
    word_count_with_buffer_size(input, read_size, threads, 1 << 20)
}

pub fn word_count_with_buffer_size(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    buffer_size: usize,
) -> Result<WordCount, String> {
    word_count_with_open_callback(input, read_size, threads, buffer_size, |_| Ok(()))
}

pub fn word_count_with_open_callback(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    buffer_size: usize,
    opened: impl FnMut(i64) -> Result<(), String>,
) -> Result<WordCount, String> {
    word_count_with_profile(input, read_size, threads, buffer_size, false, opened)
}

pub fn word_count_with_source_open_callback(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    buffer_size: usize,
    opened: impl FnMut(i64) -> Result<(), String>,
) -> Result<WordCount, String> {
    word_count_with_profile(input, read_size, threads, buffer_size, true, opened)
}

fn word_count_with_profile(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    buffer_size: usize,
    source_profile: bool,
    opened: impl FnMut(i64) -> Result<(), String>,
) -> Result<WordCount, String> {
    let input = open_two_column_input(input.as_ref(), buffer_size, source_profile)?;
    word_count_from_input(input, read_size, threads, opened)
}

fn word_count_from_input(
    input: StreamingFile,
    read_size: i64,
    threads: usize,
    mut opened: impl FnMut(i64) -> Result<(), String>,
) -> Result<WordCount, String> {
    word_count_from_input_io(input, read_size, threads, |size| {
        opened(size).map_err(std::io::Error::other)
    })
    .map_err(|error| error.to_string())
}
fn word_count_from_input_io(
    mut input: StreamingFile,
    read_size: i64,
    threads: usize,
    mut opened: impl FnMut(i64) -> std::io::Result<()>,
) -> std::io::Result<WordCount> {
    let file_size = input.size();
    opened(file_size)?;
    let mut records = 0i64;

    loop {
        let table = input.read_io(read_size, threads)?;
        records = records
            .checked_add(table.size())
            .ok_or_else(|| std::io::Error::other("TSV record count overflow"))?;
        if table.empty() {
            break;
        }
    }

    Ok(WordCount { file_size, records })
}

/// Project the first column of a two-column TSV input to `output`.
///
/// Records retain their input order, including empty first fields. Supplying
/// the writer explicitly is the safe counterpart of the C++ write-only `File`
/// constructed with an empty path (standard output).
pub fn cut(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
) -> Result<CutStats, String> {
    cut_with_buffer_size(input, read_size, threads, output, 1 << 20)
}

pub fn cut_with_buffer_size(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
    buffer_size: usize,
) -> Result<CutStats, String> {
    cut_with_open_callback(input, read_size, threads, output, buffer_size, |_| Ok(()))
}

pub fn cut_with_open_callback(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
    buffer_size: usize,
    opened: impl FnMut(i64) -> Result<(), String>,
) -> Result<CutStats, String> {
    cut_with_profile(
        input,
        read_size,
        threads,
        output,
        buffer_size,
        false,
        opened,
    )
}

pub fn cut_with_source_open_callback(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
    buffer_size: usize,
    opened: impl FnMut(i64) -> Result<(), String>,
) -> Result<CutStats, String> {
    cut_with_profile(input, read_size, threads, output, buffer_size, true, opened)
}

fn cut_with_profile(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
    buffer_size: usize,
    source_profile: bool,
    opened: impl FnMut(i64) -> Result<(), String>,
) -> Result<CutStats, String> {
    let input = open_two_column_input(input.as_ref(), buffer_size, source_profile)?;
    cut_from_input(input, read_size, threads, output, opened)
}

fn cut_from_input(
    input: StreamingFile,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
    mut opened: impl FnMut(i64) -> Result<(), String>,
) -> Result<CutStats, String> {
    cut_from_input_io(input, read_size, threads, output, |size| {
        opened(size).map_err(std::io::Error::other)
    })
    .map_err(|error| error.to_string())
}
fn cut_from_input_io(
    mut input: StreamingFile,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
    mut opened: impl FnMut(i64) -> std::io::Result<()>,
) -> std::io::Result<CutStats> {
    let file_size = input.size();
    opened(file_size)?;
    let mut records = 0i64;
    loop {
        let table = input.read_io(read_size, threads)?;
        if table.empty() {
            break;
        }
        // Source Table::map claims 1024-row batches in parallel and sends
        // the byte buffers through an ordered output queue.
        if threads == 0 {
            continue;
        }
        let batches = (table.size() as usize).div_ceil(1024);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.min(batches))
            .build()
            .map_err(std::io::Error::other)?;
        let buffers = pool.install(|| {
            (0..batches)
                .into_par_iter()
                .map(|batch| {
                    let mut bytes = Vec::new();
                    for row in batch * 1024..((batch + 1) * 1024).min(table.size() as usize) {
                        bytes.extend_from_slice(table.record(row as i64).string_bytes(0));
                        bytes.push(b'\n');
                    }
                    bytes
                })
                .collect::<Vec<_>>()
        });
        for bytes in buffers {
            output.write_all(&bytes).map_err(std::io::Error::other)?;
        }
        records = records
            .checked_add(table.size())
            .ok_or_else(|| std::io::Error::other("TSV record count overflow"))?;
    }

    Ok(CutStats { file_size, records })
}

/// Source command API preserving actual file-open errors. The established
/// string-returning wrappers above retain their library signatures.
pub fn word_count_with_source_open_callback_io(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    buffer_size: usize,
    mut opened: impl FnMut(i64) -> Result<(), String>,
) -> std::io::Result<WordCount> {
    let input = source_input_io(input.as_ref(), buffer_size)?;
    word_count_from_input_io(input, read_size, threads, |size| {
        opened(size).map_err(std::io::Error::other)
    })
}

pub fn cut_with_source_open_callback_io(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
    buffer_size: usize,
    mut opened: impl FnMut(i64) -> Result<(), String>,
) -> std::io::Result<CutStats> {
    let input = source_input_io(input.as_ref(), buffer_size)?;
    cut_from_input_io(input, read_size, threads, output, |size| {
        opened(size).map_err(std::io::Error::other)
    })
}

pub fn word_count_with_source_open_callback_io_typed(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    buffer_size: usize,
    opened: impl FnMut(i64) -> std::io::Result<()>,
) -> std::io::Result<WordCount> {
    word_count_from_input_io(
        source_input_io(input.as_ref(), buffer_size)?,
        read_size,
        threads,
        opened,
    )
}
pub fn cut_with_source_open_callback_io_typed(
    input: impl AsRef<Path>,
    read_size: i64,
    threads: usize,
    output: &mut impl Write,
    buffer_size: usize,
    opened: impl FnMut(i64) -> std::io::Result<()>,
) -> std::io::Result<CutStats> {
    cut_from_input_io(
        source_input_io(input.as_ref(), buffer_size)?,
        read_size,
        threads,
        output,
        opened,
    )
}

fn source_input_io(path: &Path, buffer_size: usize) -> std::io::Result<StreamingFile> {
    StreamingFile::open_source_path_io(vec![Type::String, Type::String], path, buffer_size)
}

fn open_two_column_input(
    path: &Path,
    buffer_size: usize,
    source_profile: bool,
) -> Result<StreamingFile, String> {
    let path = path
        .to_str()
        .ok_or_else(|| "TSV input path is not valid UTF-8".to_string())?;
    if source_profile {
        StreamingFile::open_source(vec![Type::String, Type::String], path, buffer_size)
    } else {
        StreamingFile::open(vec![Type::String, Type::String], path, buffer_size)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_path(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-tools-tsv-{label}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[cfg(unix)]
    #[test]
    fn source_tsv_raw_filename_keeps_byte_rows_and_typed_open_errors() {
        let root = temp_path("raw-name");
        std::fs::create_dir(&root).unwrap();
        let (path, _) = crate::util::test_os_path_or_unicode(&root, b"rows\xff.tsv", "rows-é.tsv");
        std::fs::write(&path, b"a\xfe\t1\nb\t2\n").unwrap();
        let mut output = Vec::new();
        let result =
            cut_with_source_open_callback_io_typed(&path, 4096, 1, &mut output, 64, |_| Ok(()))
                .unwrap();
        assert_eq!(result.records, 2);
        assert_eq!(output, b"a\xfe\nb\n");
        std::fs::remove_file(&path).unwrap();
        let error = word_count_with_source_open_callback_io_typed(&path, 4096, 1, 64, |_| Ok(()))
            .unwrap_err();
        assert_eq!(
            crate::util::io::file_open_error::find(&error).unwrap().path,
            path
        );
        std::fs::remove_dir(root).unwrap();
    }

    #[test]
    fn word_count_reports_source_file_size_and_all_records() {
        let path = temp_path("word-count");
        let input = "alpha\t1\n\t2\ngamma\t3";
        std::fs::write(&path, input).unwrap();

        let count = word_count(&path, 1, 3).unwrap();

        assert_eq!(count.file_size, input.len() as i64);
        assert_eq!(count.records, 3);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn cut_projects_first_column_in_order() {
        let path = temp_path("cut");
        let input = "alpha\t1\n\t2\ngamma\t3\n";
        std::fs::write(&path, input).unwrap();
        let mut output = Vec::new();

        let stats = cut(&path, 1, 4, &mut output).unwrap();

        assert_eq!(stats.file_size, input.len() as i64);
        assert_eq!(stats.records, 3);
        assert_eq!(output, b"alpha\n\ngamma\n");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn both_tools_reject_rows_missing_the_second_column() {
        let path = temp_path("missing-column");
        std::fs::write(&path, "alpha\t1\nbroken\n").unwrap();

        assert_eq!(
            word_count(&path, 1024, 1).unwrap_err(),
            "Missing fields in input line"
        );
        let mut output = Vec::new();
        assert_eq!(
            cut(&path, 1024, 1, &mut output).unwrap_err(),
            "Missing fields in input line"
        );
        assert!(output.is_empty());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn cut_propagates_writer_failure() {
        struct FailingWriter;

        impl Write for FailingWriter {
            fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::other("closed"))
            }

            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }

        let path = temp_path("writer-error");
        std::fs::write(&path, "alpha\t1\n").unwrap();

        let error = cut(&path, 1024, 1, &mut FailingWriter).unwrap_err();

        assert_eq!(error, "closed");
        std::fs::remove_file(path).unwrap();
    }
}
