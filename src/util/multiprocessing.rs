use std::fmt::Display;
use std::io::{Read, Write};

pub fn split(str_: &str, delim: char) -> Vec<String> {
    let mut segments = Vec::new();
    let mut segment = String::new();
    for c in str_.chars() {
        if c == delim {
            segments.push(std::mem::take(&mut segment));
        } else {
            segment.push(c);
        }
    }
    if !str_.is_empty() && !str_.ends_with(delim) {
        segments.push(segment);
    }
    segments
}

pub fn join(tokens: &[String], delim: char) -> String {
    let mut buf = String::new();
    for (i, token) in tokens.iter().enumerate() {
        buf.push_str(token);
        if i + 1 < tokens.len() {
            buf.push(delim);
        }
    }
    buf
}

pub fn quote(str_: &str) -> String {
    format!("\"{str_}\"")
}

pub fn unquote(str_: &str) -> String {
    if str_.len() >= 2 && str_.starts_with('"') && str_.ends_with('"') {
        str_[1..str_.len() - 1].to_string()
    } else {
        str_.to_string()
    }
}

pub fn copy(src_file_name: &str, dst_file_name: &str) -> std::io::Result<()> {
    let mut src = std::fs::File::open(src_file_name)?;
    let mut dst = std::fs::File::create(dst_file_name)?;
    std::io::copy(&mut src, &mut dst)?;
    Ok(())
}

pub fn join_path(path_1: &str, path_2: &str) -> String {
    #[cfg(unix)]
    const SEP: &str = "/";
    #[cfg(not(unix))]
    const SEP: &str = "\\";
    format!("{path_1}{SEP}{path_2}")
}

pub fn file_exists(file_name: &str) -> bool {
    std::fs::File::open(file_name).is_ok()
}

pub trait ScalarBytes: Copy + Sized {
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

pub fn load_scalar<R: Read, T: ScalarBytes>(is: &mut R, v: &mut T) -> std::io::Result<()> {
    let mut buf = vec![0u8; T::SIZE];
    is.read_exact(&mut buf)?;
    *v = T::from_ne_bytes(&buf);
    Ok(())
}

pub fn save_scalar<W: Write, T: ScalarBytes>(os: &mut W, v: T) -> std::io::Result<()> {
    let mut buf = Vec::with_capacity(T::SIZE);
    v.write_ne_bytes(&mut buf);
    os.write_all(&buf)
}

pub fn load_string<R: Read>(is: &mut R, s: &mut String) -> std::io::Result<()> {
    let mut size = 0usize;
    load_scalar(is, &mut size)?;
    let mut v = vec![0u8; size];
    is.read_exact(&mut v)?;
    *s = String::from_utf8(v)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    Ok(())
}

pub fn save_string<W: Write>(os: &mut W, s: &str) -> std::io::Result<()> {
    save_scalar(os, s.len())?;
    os.write_all(s.as_bytes())
}

pub fn load_vector<R: Read, T: ScalarBytes>(is: &mut R, v: &mut Vec<T>) -> std::io::Result<()> {
    let mut size = 0usize;
    load_scalar(is, &mut size)?;
    v.clear();
    v.reserve(size);
    for _ in 0..size {
        let mut x = T::from_ne_bytes(&vec![0u8; T::SIZE]);
        load_scalar(is, &mut x)?;
        v.push(x);
    }
    Ok(())
}

pub fn save_vector<W: Write, T: ScalarBytes>(os: &mut W, v: &[T]) -> std::io::Result<()> {
    save_scalar(os, v.len())?;
    for &x in v {
        save_scalar(os, x)?;
    }
    Ok(())
}

pub fn append_label<T: Display>(str_: &str, label: T, width: usize) -> String {
    format!("{str_}{label:0width$}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_string_path_helpers() {
        assert_eq!(split("a,b,c", ','), vec!["a", "b", "c"]);
        assert_eq!(split("a,b,", ','), vec!["a", "b"]);
        assert_eq!(split("a,,b", ','), vec!["a", "", "b"]);
        assert_eq!(split(",", ','), vec![""]);
        assert_eq!(split("", ','), Vec::<String>::new());
        assert_eq!(join(&["a".to_string(), "b".to_string()], ':'), "a:b");
        assert_eq!(quote("x"), "\"x\"");
        assert_eq!(unquote("\"x\""), "x");
        assert_eq!(unquote("x"), "x");
        assert!(join_path("a", "b").ends_with("a/b") || join_path("a", "b").ends_with("a\\b"));
        assert_eq!(append_label("x", 7, 4), "x0007");
    }

    #[test]
    fn test_scalar_string_vector_io() {
        let mut buf = Vec::new();
        save_scalar(&mut buf, 7u32).unwrap();
        save_string(&mut buf, "abc").unwrap();
        save_vector(&mut buf, &[1i16, 2, 3]).unwrap();
        let mut cursor = std::io::Cursor::new(buf);

        let mut x = 0u32;
        load_scalar(&mut cursor, &mut x).unwrap();
        assert_eq!(x, 7);
        let mut s = String::new();
        load_string(&mut cursor, &mut s).unwrap();
        assert_eq!(s, "abc");
        let mut v = Vec::<i16>::new();
        load_vector(&mut cursor, &mut v).unwrap();
        assert_eq!(v, vec![1, 2, 3]);
    }

    #[test]
    fn test_copy_and_exists() {
        let src = std::env::temp_dir().join(format!("diamond-rs-mp-src-{}", std::process::id()));
        let dst = std::env::temp_dir().join(format!("diamond-rs-mp-dst-{}", std::process::id()));
        std::fs::write(&src, b"abc").unwrap();
        assert!(file_exists(src.to_str().unwrap()));
        copy(src.to_str().unwrap(), dst.to_str().unwrap()).unwrap();
        assert_eq!(std::fs::read(&dst).unwrap(), b"abc");
        let _ = std::fs::remove_file(src);
        let _ = std::fs::remove_file(dst);
    }
}
