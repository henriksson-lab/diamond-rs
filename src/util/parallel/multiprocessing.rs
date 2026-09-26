//! Miscellaneous helpers used by DIAMOND's multiprocessing implementation.
//!
//! This mirrors `diamond/src/util/parallel/multiprocessing.{h,cpp}`. The
//! binary serialization deliberately uses native byte order and native
//! `usize`, just like the original writes `size_t` directly.

use std::fmt::Display;
use std::fs::File;
use std::io::{Read, Write};

/// Split `text` on `delim`, matching `std::getline`'s treatment of trailing
/// delimiters (there is no final empty field).
pub fn split(text: &str, delim: char) -> Vec<String> {
    let mut segments = Vec::new();
    let mut segment = String::new();
    for c in text.chars() {
        if c == delim {
            segments.push(std::mem::take(&mut segment));
        } else {
            segment.push(c);
        }
    }
    if !text.is_empty() && !text.ends_with(delim) {
        segments.push(segment);
    }
    segments
}

pub fn join(tokens: &[String], delim: char) -> String {
    let mut result = String::new();
    for (i, token) in tokens.iter().enumerate() {
        result.push_str(token);
        if i + 1 < tokens.len() {
            result.push(delim);
        }
    }
    result
}

pub fn quote(text: &str) -> String {
    format!("\"{text}\"")
}

pub fn unquote(text: &str) -> String {
    if text.len() >= 2 && text.starts_with('"') && text.ends_with('"') {
        text[1..text.len() - 1].to_owned()
    } else {
        text.to_owned()
    }
}

/// Copy a file byte-for-byte.
///
/// The destination is opened even if the source cannot be opened. This
/// preserves the original's construction order (`ifstream src; ofstream
/// dst;`), including its truncation/creation side effect on `dst`.
pub fn copy(src_file_name: &str, dst_file_name: &str) -> std::io::Result<()> {
    let src = File::open(src_file_name);
    let mut dst = File::create(dst_file_name)?;
    let mut src = src?;
    std::io::copy(&mut src, &mut dst)?;
    Ok(())
}

/// Join two path components using the platform separator without normalizing
/// either component, as in the original string-concatenating helper.
pub fn join_path(path_1: &str, path_2: &str) -> String {
    #[cfg(unix)]
    const SEP: &str = "/";
    #[cfg(not(unix))]
    const SEP: &str = "\\";
    format!("{path_1}{SEP}{path_2}")
}

pub fn file_exists(file_name: &str) -> bool {
    File::open(file_name).is_ok()
}

/// Types that can safely use the original native-representation scalar and
/// vector serialization format.
pub trait ScalarBytes: Copy + Default + Sized {
    const SIZE: usize;
    fn from_ne_bytes(bytes: &[u8]) -> Self;
    fn write_ne_bytes(self, out: &mut Vec<u8>);
}

macro_rules! impl_scalar_bytes {
    ($($t:ty),* $(,)?) => {
        $(
            impl ScalarBytes for $t {
                const SIZE: usize = std::mem::size_of::<$t>();

                fn from_ne_bytes(bytes: &[u8]) -> Self {
                    <$t>::from_ne_bytes(bytes.try_into().unwrap())
                }

                fn write_ne_bytes(self, out: &mut Vec<u8>) {
                    out.extend_from_slice(&self.to_ne_bytes());
                }
            }
        )*
    };
}

impl_scalar_bytes!(u8, i8, u16, i16, u32, i32, u64, i64, usize, isize, f32, f64);

pub fn load_scalar<R: Read, T: ScalarBytes>(reader: &mut R, value: &mut T) -> std::io::Result<()> {
    let mut bytes = vec![0; T::SIZE];
    reader.read_exact(&mut bytes)?;
    *value = T::from_ne_bytes(&bytes);
    Ok(())
}

pub fn save_scalar<W: Write, T: ScalarBytes>(writer: &mut W, value: T) -> std::io::Result<()> {
    let mut bytes = Vec::with_capacity(T::SIZE);
    value.write_ne_bytes(&mut bytes);
    writer.write_all(&bytes)
}

pub fn load_string<R: Read>(reader: &mut R, value: &mut String) -> std::io::Result<()> {
    let mut size = 0usize;
    load_scalar(reader, &mut size)?;
    let mut bytes = vec![0; size];
    reader.read_exact(&mut bytes)?;

    // The C++ implementation appends a NUL then assigns from `char*`, so an
    // embedded NUL terminates the resulting std::string.
    let logical_end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    *value = String::from_utf8(bytes[..logical_end].to_vec())
        .map_err(|error| std::io::Error::new(std::io::ErrorKind::InvalidData, error))?;
    Ok(())
}

/// Companion writer for `load_string` (the upstream header only supplies the
/// reader, but this preserves the pre-existing Rust API).
pub fn save_string<W: Write>(writer: &mut W, value: &str) -> std::io::Result<()> {
    save_scalar(writer, value.len())?;
    writer.write_all(value.as_bytes())
}

pub fn load_vector<R: Read, T: ScalarBytes>(
    reader: &mut R,
    value: &mut Vec<T>,
) -> std::io::Result<()> {
    let mut size = 0usize;
    load_scalar(reader, &mut size)?;
    value.clear();
    value.reserve(size);
    for _ in 0..size {
        let mut element = T::default();
        load_scalar(reader, &mut element)?;
        value.push(element);
    }
    Ok(())
}

pub fn save_vector<W: Write, T: ScalarBytes>(writer: &mut W, value: &[T]) -> std::io::Result<()> {
    save_scalar(writer, value.len())?;
    for &element in value {
        save_scalar(writer, element)?;
    }
    Ok(())
}

pub const DEFAULT_LABEL_WIDTH: usize = 6;

/// Append a zero-filled label with a minimum field width.
///
/// `std::setw` uses right alignment, so the fill precedes a minus sign (for
/// example, `-1` at width six becomes `0000-1`). Formatting the value first
/// and padding that complete representation reproduces that behavior.
pub fn append_label<T: Display>(text: &str, label: T, width: usize) -> String {
    let label = label.to_string();
    let padding = width.saturating_sub(label.len());
    let mut result = String::with_capacity(text.len() + padding + label.len());
    result.push_str(text);
    result.extend(std::iter::repeat('0').take(padding));
    result.push_str(&label);
    result
}

/// Rust spelling of the C++ overload's default `width = 6` behavior.
pub fn append_label_default<T: Display>(text: &str, label: T) -> String {
    append_label(text, label, DEFAULT_LABEL_WIDTH)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn split_matches_getline_edges() {
        assert_eq!(split("", ','), Vec::<String>::new());
        assert_eq!(split(",", ','), vec![""]);
        assert_eq!(split(",,", ','), vec!["", ""]);
        assert_eq!(split("a,,b,", ','), vec!["a", "", "b"]);
        assert_eq!(split("å:β", ':'), vec!["å", "β"]);
        assert_eq!(join(&["a".into(), String::new(), "b".into()], ':'), "a::b");
    }

    #[test]
    fn quote_unquote_only_strip_a_matching_outer_pair() {
        assert_eq!(quote("a\"b"), "\"a\"b\"");
        assert_eq!(unquote("\"a\"b\""), "a\"b");
        assert_eq!(unquote("\""), "\"");
        assert_eq!(unquote("\"\""), "");
        assert_eq!(unquote("plain"), "plain");
    }

    #[test]
    fn serialization_is_native_and_round_trips() {
        let mut encoded = Vec::new();
        save_scalar(&mut encoded, 0x1234u16).unwrap();
        save_vector(&mut encoded, &[-2i32, 7]).unwrap();
        assert_eq!(&encoded[..2], &0x1234u16.to_ne_bytes());

        let mut reader = Cursor::new(encoded);
        let mut scalar = 0u16;
        let mut vector = Vec::<i32>::new();
        load_scalar(&mut reader, &mut scalar).unwrap();
        load_vector(&mut reader, &mut vector).unwrap();
        assert_eq!(scalar, 0x1234);
        assert_eq!(vector, vec![-2, 7]);
    }

    #[test]
    fn load_string_matches_cpp_c_string_assignment() {
        let mut encoded = Vec::new();
        save_scalar(&mut encoded, 5usize).unwrap();
        encoded.extend_from_slice(b"ab\0cd");
        let mut decoded = String::new();
        load_string(&mut Cursor::new(encoded), &mut decoded).unwrap();
        assert_eq!(decoded, "ab");
    }

    #[test]
    fn label_padding_matches_setw_right_alignment() {
        assert_eq!(append_label("x", 7, 4), "x0007");
        assert_eq!(append_label("x", 12345, 3), "x12345");
        assert_eq!(append_label("x", -1, 6), "x0000-1");
        assert_eq!(append_label_default("block_", 42), "block_000042");
    }

    #[test]
    fn copy_and_file_helpers_preserve_cpp_edges() {
        let root = std::env::temp_dir().join(format!(
            "diamond-rs-multiprocessing-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&root).unwrap();
        let src = root.join("src");
        let dst = root.join("dst");
        std::fs::write(&src, b"abc\0def").unwrap();
        copy(src.to_str().unwrap(), dst.to_str().unwrap()).unwrap();
        assert_eq!(std::fs::read(&dst).unwrap(), b"abc\0def");
        assert!(file_exists(dst.to_str().unwrap()));

        std::fs::write(&dst, b"old contents").unwrap();
        let missing = root.join("missing");
        assert!(copy(missing.to_str().unwrap(), dst.to_str().unwrap()).is_err());
        assert_eq!(std::fs::read(&dst).unwrap(), b"");

        let expected = if cfg!(unix) {
            "left/right"
        } else {
            "left\\right"
        };
        assert_eq!(join_path("left", "right"), expected);
        let doubled = if cfg!(unix) {
            "left///right"
        } else {
            "left/\\/right"
        };
        assert_eq!(join_path("left/", "/right"), doubled);
        std::fs::remove_dir_all(root).unwrap();
    }
}
