//! Source's exception-disabled text streams and classic formatted extraction.
//! Shared Source ifstream/ofstream contracts for concrete workflow owners.
use crate::util::io::{cstdio_file::OwnedCstdio, IoErrorOwner, RawDiagnosticError};
use std::path::Path;

pub(crate) fn path_error(owner: &IoErrorOwner, prefix: &[u8], path: &Path) -> String {
    let result = (|| -> std::io::Result<()> {
        let mut message = prefix.to_vec();
        message.extend_from_slice(
            &crate::util::io::source_os_bytes(path.as_os_str()).map_err(std::io::Error::other)?,
        );
        Err(RawDiagnosticError(message).into_io())
    })();
    owner.capture(result).unwrap_err()
}

/// Retained Source ifstream. Construction opens only; extraction advances the
/// real text-mode owner and preserves eofbit/failbit across calls.
pub(crate) struct SourceTextInput {
    stream: Option<OwnedCstdio>,
    pending: Option<u8>,
    terminal: bool,
    eofbit: bool,
    failed: bool,
    read_failed: bool,
}
impl SourceTextInput {
    pub(crate) fn open(path: &Path, owner: &IoErrorOwner, prefix: &[u8]) -> Result<Self, String> {
        let stream = OwnedCstdio::open_file_role_path(path, "r")
            .map_err(|_| path_error(owner, prefix, path))?;
        Ok(Self {
            stream: Some(stream),
            pending: None,
            terminal: false,
            eofbit: false,
            failed: false,
            read_failed: false,
        })
    }
    pub(crate) fn open_unchecked(path: &Path) -> Self {
        let stream = OwnedCstdio::open_file_role_path(path, "r").ok();
        let failed = stream.is_none();
        Self {
            stream,
            pending: None,
            terminal: failed,
            eofbit: false,
            failed,
            read_failed: false,
        }
    }
    fn peek(&mut self) -> Option<u8> {
        if self.pending.is_some() {
            return self.pending;
        }
        if self.terminal {
            return None;
        }
        let mut byte = [0];
        let Some(stream) = self.stream.as_mut() else {
            self.failed = true;
            return None;
        };
        if stream
            .read_file_role(&mut byte)
            .map(|(n, _)| n)
            .unwrap_or(0)
            == 1
        {
            self.pending = Some(byte[0]);
        } else {
            self.terminal = true;
            self.eofbit = stream.eof().unwrap_or(false);
            self.read_failed = !self.eofbit;
        }
        self.pending
    }
    fn take(&mut self) -> Option<u8> {
        let value = self.peek();
        self.pending = None;
        value
    }
    pub(crate) fn good(&self) -> bool {
        !self.failed && !self.read_failed
    }
    pub(crate) fn line(&mut self) -> Option<Vec<u8>> {
        if !self.good() || self.eofbit {
            self.failed = true;
            return None;
        }
        let mut bytes = Vec::new();
        loop {
            match self.take() {
                Some(b'\n') => return Some(bytes),
                Some(byte) => bytes.push(byte),
                None if self.read_failed || bytes.is_empty() => {
                    self.failed = true;
                    return None;
                }
                None => return Some(bytes),
            }
        }
    }
    /// A sentry failure leaves the previous value intact. Once num_get runs,
    /// an invalid conversion stores zero, and overflow stores UINT64_MAX,
    /// both with failbit. Only the initial sentry/EOF leaves None undefined.
    pub(crate) fn number_into(&mut self, value: &mut Option<u64>) {
        if !self.good() || self.eofbit {
            self.failed = true;
            return;
        }
        while self.peek().is_some_and(classic_space) {
            self.take();
        }
        if self.peek().is_none() {
            self.failed = true;
            return;
        }
        let negative = self.peek() == Some(b'-');
        if matches!(self.peek(), Some(b'-' | b'+')) {
            self.take();
        }
        let mut magnitude = 0u64;
        let mut digits = false;
        let mut overflow = false;
        while let Some(byte) = self.peek().filter(u8::is_ascii_digit) {
            self.take();
            digits = true;
            match magnitude
                .checked_mul(10)
                .and_then(|n| n.checked_add((byte - b'0') as u64))
            {
                Some(next) => magnitude = next,
                None => overflow = true,
            }
        }
        if self.read_failed {
            self.failed = true;
            return;
        }
        if !digits {
            *value = Some(0);
            self.failed = true;
        } else if overflow {
            *value = Some(u64::MAX);
            self.failed = true;
        } else {
            *value = Some(if negative {
                magnitude.wrapping_neg()
            } else {
                magnitude
            });
        }
    }
    pub(crate) fn number(&mut self) -> Option<u64> {
        let mut value = None;
        self.number_into(&mut value);
        if self.good() {
            value
        } else {
            None
        }
    }
    pub(crate) fn close(&mut self) {
        if let Some(mut stream) = self.stream.take() {
            let _ = stream.close_file_role();
        }
        self.terminal = true;
        self.pending = None;
    }
}
impl Drop for SourceTextInput {
    fn drop(&mut self) {
        self.close();
    }
}
fn classic_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\r' | b'\n' | 0x0b | 0x0c)
}

pub(crate) struct TextOutput {
    stream: Option<OwnedCstdio>,
    failed: bool,
    sticky: bool,
}
impl TextOutput {
    pub(crate) fn open(
        path: &Path,
        owner: &IoErrorOwner,
        prefix: &[u8],
        sticky: bool,
    ) -> Result<Self, String> {
        let stream = OwnedCstdio::open_file_role_path(path, "w")
            .map_err(|_| path_error(owner, prefix, path))?;
        Ok(Self {
            stream: Some(stream),
            failed: false,
            sticky,
        })
    }
    pub(crate) fn open_unchecked(path: &Path, sticky: bool) -> Self {
        let stream = OwnedCstdio::open_file_role_path(path, "w").ok();
        let failed = stream.is_none();
        Self {
            stream,
            failed,
            sticky,
        }
    }
    pub(crate) fn endl(&mut self) {
        self.write(b"\n");
        if !self.failed {
            if let Some(stream) = self.stream.as_mut() {
                self.failed = !stream.flush_text_output().unwrap_or(false);
            }
        }
    }
    pub(crate) fn write(&mut self, bytes: &[u8]) {
        if self.sticky && self.failed {
            return;
        }
        if let Some(stream) = self.stream.as_mut() {
            self.failed |= !stream.write_text_output(bytes).unwrap_or(false);
        }
    }
}
impl Drop for TextOutput {
    fn drop(&mut self) {
        if let Some(mut stream) = self.stream.take() {
            let _ = stream.close_file_role();
        }
    }
}

/// `operator>> uint64_t` advances through the numeric prefix, preserves its
/// failbit, and wraps a representable negative magnitude. No UTF-8 conversion.
pub(crate) struct Fields<'a> {
    bytes: &'a [u8],
    at: usize,
    failed: bool,
    eofbit: bool,
}
impl<'a> Fields<'a> {
    pub(crate) fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            at: 0,
            failed: false,
            eofbit: false,
        }
    }
    fn space(&mut self) {
        while self.at < self.bytes.len()
            && matches!(
                self.bytes[self.at],
                b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c
            )
        {
            self.at += 1;
        }
    }
    pub(crate) fn word(&mut self) -> Option<&'a [u8]> {
        if self.failed || self.eofbit {
            self.failed = true;
            return None;
        }
        self.space();
        let begin = self.at;
        while self.at < self.bytes.len()
            && !matches!(
                self.bytes[self.at],
                b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c
            )
        {
            self.at += 1;
        }
        if begin == self.at {
            self.failed = true;
            return None;
        }
        self.eofbit = self.at == self.bytes.len();
        Some(&self.bytes[begin..self.at])
    }
    pub(crate) fn number_into(&mut self, value: &mut Option<u64>) {
        if self.failed || self.eofbit {
            self.failed = true;
            return;
        }
        self.space();
        if self.at == self.bytes.len() {
            self.eofbit = true;
            self.failed = true;
            return;
        }
        let negative = self.bytes.get(self.at) == Some(&b'-');
        if matches!(self.bytes.get(self.at), Some(b'-' | b'+')) {
            self.at += 1;
        }
        let begin = self.at;
        let mut magnitude = 0u64;
        let mut overflow = false;
        while let Some(&digit) = self.bytes.get(self.at).filter(|b| b.is_ascii_digit()) {
            match magnitude
                .checked_mul(10)
                .and_then(|n| n.checked_add((digit - b'0') as u64))
            {
                Some(next) => magnitude = next,
                None => overflow = true,
            }
            self.at += 1;
        }
        self.eofbit = self.at == self.bytes.len();
        if self.at == begin {
            *value = Some(0);
            self.failed = true;
        } else if overflow {
            *value = Some(u64::MAX);
            self.failed = true;
        } else {
            *value = Some(if negative {
                magnitude.wrapping_neg()
            } else {
                magnitude
            });
        }
    }
    pub(crate) fn number(&mut self) -> Option<u64> {
        let mut value = None;
        self.number_into(&mut value);
        if self.failed {
            None
        } else {
            value
        }
    }
}

/// Expected bytes for Source ofstream/fprintf text files, never for binary
/// manifests or captured data. Project only fixtures through CRT text mode.
#[cfg(test)]
pub(crate) fn expected_text(bytes: &[u8]) -> Vec<u8> {
    let mut expected = Vec::with_capacity(bytes.len());
    for &byte in bytes {
        if cfg!(windows) && byte == b'\n' {
            expected.push(b'\r');
        }
        expected.push(byte);
    }
    expected
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_numeric_extraction_keeps_values_after_sentry_failure() {
        let root = std::env::temp_dir().join(format!(
            "diamond-multinode-number-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&root, b"-1+2 17x").unwrap();
        let mut input =
            SourceTextInput::open(&root, &IoErrorOwner::default(), b"Error opening file ").unwrap();
        let mut value = None;
        input.number_into(&mut value);
        assert_eq!(value, Some(u64::MAX));
        input.number_into(&mut value);
        assert_eq!(value, Some(2));
        input.number_into(&mut value);
        assert_eq!(value, Some(17));
        input.number_into(&mut value);
        assert_eq!(value, Some(0));
        assert!(!input.good());
        value = Some(99);
        input.number_into(&mut value);
        assert_eq!(value, Some(99));
        drop(input);
        std::fs::write(&root, b"8").unwrap();
        let mut input =
            SourceTextInput::open(&root, &IoErrorOwner::default(), b"Error opening file ").unwrap();
        value = None;
        input.number_into(&mut value);
        assert_eq!(value, Some(8));
        assert!(input.good());
        input.number_into(&mut value);
        assert_eq!(value, Some(8));
        assert!(!input.good());
        drop(input);
        std::fs::remove_file(root).unwrap();
    }
    #[cfg(unix)]
    #[test]
    fn retained_line_and_number_do_not_commit_pending_bytes_on_read_failure() {
        let directory = std::env::temp_dir().join(format!(
            "diamond-multinode-partial-read-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&directory).unwrap();
        let empty = directory.join("empty");
        std::fs::write(&empty, b"").unwrap();
        let mut clean =
            SourceTextInput::open(&empty, &IoErrorOwner::default(), b"Error opening file ")
                .unwrap();
        clean.pending = Some(b'x');
        assert_eq!(clean.line(), Some(b"x".to_vec()));
        drop(clean);
        let mut failed =
            SourceTextInput::open(&directory, &IoErrorOwner::default(), b"Error opening file ")
                .unwrap();
        failed.pending = Some(b'x');
        assert_eq!(failed.line(), None);
        assert!(failed.read_failed);
        drop(failed);
        let mut failed =
            SourceTextInput::open(&directory, &IoErrorOwner::default(), b"Error opening file ")
                .unwrap();
        failed.pending = Some(b'3');
        let mut value = Some(7);
        failed.number_into(&mut value);
        assert_eq!(value, Some(7));
        assert!(!failed.good());
        drop(failed);
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn numeric_extraction_uses_prefix_cursor_and_sticky_failure() {
        let mut fields = Fields::new(b"-1+2 17\xff 99");
        assert_eq!(fields.number(), Some(u64::MAX));
        assert_eq!(fields.number(), Some(2));
        assert_eq!(fields.number(), Some(17));
        assert_eq!(fields.number(), None);
        assert_eq!(fields.number(), None);
        let mut overflow = Fields::new(b"18446744073709551616 2");
        let mut value = Some(0);
        overflow.number_into(&mut value);
        assert_eq!(value, Some(u64::MAX));
        overflow.number_into(&mut value);
        assert_eq!(value, Some(u64::MAX));
    }
    #[cfg(unix)]
    #[test]
    fn directory_read_failure_is_not_clean_eof() {
        let path = std::env::temp_dir().join(format!(
            "diamond-multinode-read-state-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&path).unwrap();
        let mut input =
            SourceTextInput::open(&path, &IoErrorOwner::default(), b"Error opening file ").unwrap();
        assert!(input.line().is_none());
        assert!(input.read_failed);
        assert!(!input.eofbit);
        drop(input);
        std::fs::remove_dir(path).unwrap();
    }
}
