//! BLAST protein alias (.pal) parsing, with opaque Source path/value bytes.
use super::volume::BlastVolume;
use crate::util::io::source_c_path;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

pub(crate) fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    name.into()
}
pub(crate) fn path_bytes(path: &Path) -> io::Result<Vec<u8>> {
    crate::util::io::path_diagnostic_bytes(path).map_err(io::Error::other)
}
pub(crate) fn raw_error(message: Vec<u8>) -> io::Error {
    crate::util::io::RawDiagnosticError(message).into_io()
}
fn bytes_path(bytes: &[u8]) -> io::Result<PathBuf> {
    crate::util::io::source_bytes_os(bytes)
        .map(Into::into)
        .map_err(io::Error::other)
}
pub(crate) fn absolute_parts(path: &Path) -> io::Result<(PathBuf, std::ffi::OsString)> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::{OsStrExt, OsStringExt};
        let original = path.as_os_str().as_bytes();
        let original = if original.is_empty() { b"." } else { original };
        let base = original
            .split(|b| *b == b'/')
            .filter(|s| !s.is_empty())
            .next_back()
            .unwrap_or_default();
        let directory = original.ends_with(b"/") || base == b"." || base == b"..";
        let joined = if original.starts_with(b"/") {
            original.to_vec()
        } else {
            let mut v = std::env::current_dir()?.as_os_str().as_bytes().to_vec();
            v.push(b'/');
            v.extend_from_slice(original);
            v
        };
        let mut components: Vec<&[u8]> = Vec::new();
        for part in joined.split(|b| *b == b'/') {
            match part {
                b"" | b"." => {}
                b".." => {
                    components.pop();
                }
                _ => components.push(part),
            }
        }
        let mut normalized = vec![b'/'];
        for (i, part) in components.iter().enumerate() {
            if i > 0 {
                normalized.push(b'/');
            }
            normalized.extend_from_slice(part);
        }
        let full = PathBuf::from(std::ffi::OsString::from_vec(normalized));
        if directory {
            Ok((full, std::ffi::OsString::new()))
        } else {
            Ok((
                full.parent().unwrap_or(Path::new("/")).to_path_buf(),
                std::ffi::OsString::from_vec(base.to_vec()),
            ))
        }
    }
    #[cfg(not(unix))]
    {
        let path = if path.as_os_str().is_empty() {
            Path::new(".")
        } else {
            path
        };
        let full = std::path::absolute(path)?;
        let raw = path_bytes(path)?;
        let directory = raw.ends_with(b"/") || raw.ends_with(b"\\") || raw == b"." || raw == b"..";
        if directory {
            Ok((full, std::ffi::OsString::new()))
        } else {
            Ok((
                full.parent().unwrap_or(Path::new("")).to_path_buf(),
                path.file_name().unwrap_or_default().to_owned(),
            ))
        }
    }
}
fn trim_bytes(text: &[u8]) -> &[u8] {
    let start = text
        .iter()
        .position(|b| !b" \t\r\n".contains(b))
        .unwrap_or(text.len());
    let end = text
        .iter()
        .rposition(|b| !b" \t\r\n".contains(b))
        .map_or(start, |i| i + 1);
    &text[start..end]
}
fn is_absolute(bytes: &[u8]) -> bool {
    bytes.first().is_some_and(|b| *b == b'/' || *b == b'\\')
        || (bytes.len() >= 3
            && bytes[0].is_ascii_alphabetic()
            && bytes[1] == b':'
            && matches!(bytes[2], b'/' | b'\\'))
}
fn relative_to(dir: &Path, bytes: &[u8]) -> io::Result<PathBuf> {
    let path = bytes_path(bytes)?;
    if is_absolute(bytes) {
        Ok(path)
    } else {
        Ok(dir.join(path))
    }
}
fn line_error(prefix: String, line: &[u8]) -> io::Error {
    let mut v = prefix.into_bytes();
    v.extend_from_slice(line);
    raw_error(v)
}
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Pal {
    pub volumes: Vec<PathBuf>,
    pub metadata: BTreeMap<String, Vec<u8>>,
    pub oid_index: Vec<u64>,
    pub sequence_count: u64,
    pub letters: u64,
    pub version: i32,
}
impl Pal {
    pub fn new(path: &str) -> Result<Self, String> {
        Self::new_path(Path::new(path)).map_err(|e| e.to_string())
    }
    pub fn new_path(path: &Path) -> io::Result<Self> {
        let supported: BTreeSet<&[u8]> = [
            b"TITLE".as_slice(),
            b"MEMB_BIT",
            b"SEQIDLIST",
            b"NSEQ",
            b"LENGTH",
            b"TAXIDLIST",
        ]
        .into_iter()
        .collect();
        let (db_dir, file) = absolute_parts(path)?;
        let mut pal = Self::default();
        let raw = path_bytes(path)?;
        if !source_c_path(&append_suffix(path, ".pal")).exists() && !raw.ends_with(b".pal") {
            pal.volumes.push(db_dir.join(file));
        } else {
            let pal_path = if raw.ends_with(b".pal") {
                path.to_owned()
            } else {
                append_suffix(path, ".pal")
            };
            let file = std::fs::File::open(source_c_path(&pal_path).as_ref()).map_err(|_| {
                let mut v = b"Unable to open PAL file: ".to_vec();
                let bytes = match path_bytes(&pal_path) {
                    Ok(bytes) => bytes,
                    Err(error) => return error,
                };
                v.extend_from_slice(&bytes);
                raw_error(v)
            })?;
            let mut reader = BufReader::new(file);
            let mut line_number = 0;
            loop {
                let mut raw = Vec::new();
                match reader.read_until(b'\n', &mut raw) {
                    Ok(0) | Err(_) => break,
                    Ok(_) => {}
                }
                line_number += 1;
                let raw = raw.split(|b| *b == b'#').next().unwrap_or_default();
                let line = trim_bytes(raw);
                if line.is_empty() {
                    continue;
                }
                let end = line
                    .iter()
                    .position(|b| matches!(b, b' ' | b'\t'))
                    .ok_or_else(|| {
                        line_error(
                            format!(
                                "Error parsing PAL file: line {line_number} is missing a value: "
                            ),
                            line,
                        )
                    })?;
                let key = &line[..end];
                let value = trim_bytes(&line[end + 1..]);
                if value.is_empty() {
                    return Err(line_error(
                        format!("Error parsing PAL file: line {line_number} has an empty value: "),
                        line,
                    ));
                }
                if key == b"DBLIST" {
                    let values = value
                        .split(|b| b" \t\r\n\x0b\x0c".contains(b))
                        .filter(|v| !v.is_empty())
                        .map(bytes_path)
                        .collect::<io::Result<Vec<_>>>()?;
                    if values.is_empty() {
                        return Err(raw_error(format!("Error parsing PAL file: DBLIST on line {line_number} does not list any volumes").into_bytes()));
                    }
                    pal.volumes.extend(values);
                    // Source deliberately revisits all accumulated volume entries.
                    for path in &mut pal.volumes {
                        let raw = path_bytes(path)?;
                        if !is_absolute(&raw) && !raw.is_empty() && !raw.starts_with(b"\"") {
                            *path = db_dir.join(&*path);
                        }
                    }
                    continue;
                }
                if !supported.contains(key) {
                    let mut v = b"Error parsing PAL file: Unsupported PAL key '".to_vec();
                    v.extend_from_slice(key);
                    v.extend_from_slice(format!("' on line {line_number}").as_bytes());
                    return Err(raw_error(v));
                }
                let key = std::str::from_utf8(key).expect("validated ASCII PAL key");
                if pal.metadata.contains_key(key) {
                    return Err(raw_error(
                        format!(
                            "Error parsing PAL file: Duplicate key '{key}' on line {line_number}"
                        )
                        .into_bytes(),
                    ));
                }
                pal.metadata.insert(key.to_owned(), value.to_vec());
            }
        }
        pal.oid_index.push(0);
        let mut it = 0;
        while it < pal.volumes.len() {
            let volume = pal.volumes[it].clone();
            let raw = path_bytes(&volume)?;
            if raw.len() >= 2 && raw.starts_with(b"\"") && raw.ends_with(b"\"") {
                let nested = relative_to(&db_dir, &raw[1..raw.len() - 1])?;
                pal.volumes.remove(it);
                it += pal.recurse(&nested, it)?;
            } else {
                let volume = BlastVolume::new_path(&volume, 0, 0, 0, false)?;
                pal.sequence_count += u64::from(volume.index().num_oids);
                pal.oid_index
                    .push(u64::from(volume.index().num_oids) + pal.oid_index.last().unwrap());
                pal.letters += volume.index().total_length;
                pal.version = volume.index().version as i32;
                it += 1;
            }
        }
        if let Some(value) = pal.metadata.get_mut("SEQIDLIST") {
            if value.ends_with(b".bsl") {
                let mut v =
                    b"Binary SEQIDLIST files(.bsl) are not supported, use text file instead : "
                        .to_vec();
                v.extend_from_slice(value);
                return Err(raw_error(v));
            }
            if !is_absolute(value) {
                *value = path_bytes(&relative_to(&db_dir, value)?)?;
            }
        }
        if let Some(value) = pal.metadata.get_mut("TAXIDLIST") {
            if !is_absolute(value) {
                *value = path_bytes(&relative_to(&db_dir, value)?)?;
            }
        }
        assert!(pal.sequence_count > 0);
        Ok(pal)
    }
    fn recurse(&mut self, path: &Path, volume_it: usize) -> io::Result<usize> {
        let pal = Self::new_path(path)?;
        let inserted = pal.volumes.len();
        self.volumes.splice(volume_it..volume_it, pal.volumes);
        for (key, value) in pal.metadata {
            if self.metadata.contains_key(&key) {
                if key == "TITLE" || key == "NSEQ" || key == "LENGTH" {
                    continue;
                }
                let mut v = format!("Duplicate key '{key}' in nested PAL file: ").into_bytes();
                v.extend_from_slice(&path_bytes(path)?);
                return Err(raw_error(v));
            }
            self.metadata.insert(key, value);
        }
        let base = *self.oid_index.last().unwrap();
        self.oid_index
            .extend(pal.oid_index.iter().skip(1).map(|oid| oid + base));
        self.sequence_count += pal.sequence_count;
        self.letters += pal.letters;
        self.version = pal.version;
        Ok(inserted)
    }
    pub fn volume(&self, oid: u64) -> i32 {
        assert!(oid < self.sequence_count);
        self.oid_index.partition_point(|&x| x <= oid) as i32 - 1
    }
}
#[cfg(test)]
fn trim(text: &str) -> String {
    String::from_utf8(trim_bytes(text.as_bytes()).to_vec()).unwrap()
}
#[cfg(test)]
fn split_whitespace(text: &str) -> Vec<String> {
    text.split_ascii_whitespace().map(str::to_owned).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_pal(contents: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-pal-errors-{}-{}.pal",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn trim_and_split_match_stream_whitespace_semantics() {
        assert_eq!(trim(" \t alpha beta\r\n"), "alpha beta");
        assert_eq!(
            split_whitespace(" one\t two\r\nthree "),
            ["one", "two", "three"]
        );
    }

    #[test]
    fn rejects_missing_unsupported_and_duplicate_values() {
        let cases = [
            ("TITLE\n", "line 1 is missing a value"),
            ("TITLE  \t\n", "line 1 is missing a value"),
            ("UNKNOWN value\n", "Unsupported PAL key 'UNKNOWN' on line 1"),
            ("TITLE one\nTITLE two\n", "Duplicate key 'TITLE' on line 2"),
        ];
        for (text, expected) in cases {
            let path = temp_pal(text);
            let error = Pal::new(path.to_str().unwrap()).unwrap_err();
            assert!(error.contains(expected), "{error:?}");
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn comments_are_removed_before_validation() {
        let path = temp_pal("# comment\nTITLE # missing after comment\n");
        let error = Pal::new(path.to_str().unwrap()).unwrap_err();
        assert!(error.contains("line 2 is missing a value"), "{error:?}");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn volume_lookup_uses_upper_bound_at_boundaries() {
        let pal = Pal {
            oid_index: vec![0, 2, 5, 9],
            sequence_count: 9,
            ..Pal::default()
        };
        assert_eq!(pal.volume(0), 0);
        assert_eq!(pal.volume(1), 0);
        assert_eq!(pal.volume(2), 1);
        assert_eq!(pal.volume(8), 2);
    }
}
