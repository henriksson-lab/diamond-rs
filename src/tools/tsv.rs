//! Small TSV command helpers.
//!
//! This mirrors `diamond/src/tools/tsv.cpp`. The C++ entry points read their
//! input, chunk size, thread count, and output stream from global state. The
//! Rust API makes those dependencies explicit and returns the values that the
//! command layer can use for verbose and message output.

use std::io::Write;
use std::path::Path;

use crate::util::tsv::{File, Flags, Type};

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
    let mut input = open_two_column_input(input.as_ref())?;
    let file_size = input.size();
    let mut records = 0i64;

    loop {
        let table = input.read_chunk(read_size, threads);
        records = records
            .checked_add(table.size())
            .ok_or_else(|| "TSV record count overflow".to_string())?;
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
    let mut input = open_two_column_input(input.as_ref())?;
    let file_size = input.size();
    let mut records = 0i64;

    loop {
        let table = input.read_chunk(read_size, threads);
        if table.empty() {
            break;
        }
        for i in 0..table.size() {
            let first = table.record(i).get(0);
            output
                .write_all(first.as_bytes())
                .and_then(|()| output.write_all(b"\n"))
                .map_err(|error| error.to_string())?;
            records = records
                .checked_add(1)
                .ok_or_else(|| "TSV record count overflow".to_string())?;
        }
    }

    Ok(CutStats { file_size, records })
}

fn open_two_column_input(path: &Path) -> Result<File, String> {
    let path = path
        .to_str()
        .ok_or_else(|| "TSV input path is not valid UTF-8".to_string())?;
    File::new(
        vec![Type::String, Type::String],
        path,
        Flags::READ,
        Default::default(),
    )
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
