//! Top-level TSV utilities.
//!
//! This mirrors `diamond/src/util/tsv/tsv.{h,cpp}`.

use crate::util::io::InputFile;

/// Count newline-delimited records using DIAMOND's input stream abstraction.
///
/// A final unterminated non-empty line counts as one record. Empty input does
/// not. `InputFile` also preserves the C++ path's gzip/zlib and zstd
/// autodetection.
pub fn count_lines(file_name: &str) -> Result<i64, String> {
    count_lines_with_buffer(file_name, 1 << 20)
}

pub fn count_lines_with_buffer(file_name: &str, buffer_size: usize) -> Result<i64, String> {
    count_lines_with_buffer_path(std::path::Path::new(file_name), buffer_size)
}

pub fn count_lines_with_buffer_path(
    file_name: &std::path::Path,
    buffer_size: usize,
) -> Result<i64, String> {
    let file = InputFile::with_buffer_size_path(file_name, 0, buffer_size)
        .map_err(|error| error.to_string())?;
    count_lines_input(file).map_err(|error| error.to_string())
}

pub fn count_lines_with_source_buffer(file_name: &str, buffer_size: usize) -> Result<i64, String> {
    count_lines_with_source_buffer_io(file_name, buffer_size).map_err(|error| error.to_string())
}
pub fn count_lines_with_source_buffer_io(
    file_name: &str,
    buffer_size: usize,
) -> std::io::Result<i64> {
    count_lines_with_source_buffer_path_io(std::path::Path::new(file_name), buffer_size)
}
pub fn count_lines_with_source_buffer_path_io(
    file_name: &std::path::Path,
    buffer_size: usize,
) -> std::io::Result<i64> {
    let file = InputFile::with_source_profile_path(file_name, 0, buffer_size)
        .map_err(std::io::Error::other)?;
    count_lines_input(file)
}
fn count_lines_input(mut file: InputFile) -> std::io::Result<i64> {
    let mut count = 0i64;
    loop {
        let mut line = Vec::new();
        let delimiter_found = file
            .read_to(&mut line, b'\n')
            .map_err(std::io::Error::other)?;
        if !line.is_empty() || delimiter_found {
            count += 1;
        }
        if !delimiter_found {
            break;
        }
    }
    file.close().map_err(std::io::Error::other)?;
    Ok(count)
}

/// In-memory counterpart with the same final-line semantics.
pub fn count_lines_str(text: &str) -> i64 {
    text.split_terminator('\n').count() as i64
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn temp_path(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-tsv-count-lines-{label}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[test]
    fn in_memory_line_count_matches_getline_eof_semantics() {
        assert_eq!(count_lines_str(""), 0);
        assert_eq!(count_lines_str("a"), 1);
        assert_eq!(count_lines_str("a\n"), 1);
        assert_eq!(count_lines_str("\n"), 1);
        assert_eq!(count_lines_str("\n\n"), 2);
        assert_eq!(count_lines_str("a\n\nb"), 3);
    }

    #[test]
    fn file_counter_accepts_non_utf8_and_unterminated_final_line() {
        let path = temp_path("binary");
        std::fs::write(&path, [0xff, b'\n', b'x']).unwrap();

        assert_eq!(count_lines(path.to_str().unwrap()).unwrap(), 2);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn file_counter_autodetects_gzip_and_zstd() {
        let input = b"first\n\nthird";
        let gzip_path = temp_path("gzip");
        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip.write_all(input).unwrap();
        std::fs::write(&gzip_path, gzip.finish().unwrap()).unwrap();

        let zstd_path = temp_path("zstd");
        std::fs::write(&zstd_path, crate::util::io::zstd_compress_for_test(input)).unwrap();

        assert_eq!(count_lines(gzip_path.to_str().unwrap()).unwrap(), 3);
        assert_eq!(count_lines(zstd_path.to_str().unwrap()).unwrap(), 3);

        std::fs::remove_file(gzip_path).unwrap();
        std::fs::remove_file(zstd_path).unwrap();
    }
}
