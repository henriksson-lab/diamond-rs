//! Source's exception-disabled text streams and classic formatted extraction.
//! These retained library helpers do not install a CLI route.
use crate::util::io::{cstdio_file::OwnedCstdio, IoErrorOwner, RawDiagnosticError};
use std::path::Path;

pub(super) fn path_error(owner: &IoErrorOwner, prefix: &[u8], path: &Path) -> String {
    let result = (|| -> std::io::Result<()> {
        let mut message = prefix.to_vec();
        message.extend_from_slice(
            &crate::util::io::source_os_bytes(path.as_os_str()).map_err(std::io::Error::other)?,
        );
        Err(RawDiagnosticError(message).into_io())
    })();
    owner.capture(result).unwrap_err()
}

/// Bytes read by an exception-disabled ifstream, together with the terminal
/// state observed on the actual CRT owner. A short ferror is not clean EOF.
pub(super) struct ReadText {
    bytes: Vec<u8>,
    clean_eof: bool,
}
impl ReadText {
    pub(super) fn lines(&self) -> impl Iterator<Item = &[u8]> {
        self.bytes
            .split_inclusive(|&byte| byte == b'\n')
            .filter_map(move |part| {
                if part.last() == Some(&b'\n') {
                    Some(&part[..part.len() - 1])
                } else if self.clean_eof {
                    Some(part)
                } else {
                    None
                }
            })
    }
}

pub(super) fn read_text(
    path: &Path,
    owner: &IoErrorOwner,
    prefix: &[u8],
) -> Result<ReadText, String> {
    let mut input =
        OwnedCstdio::open_file_role_path(path, "r").map_err(|_| path_error(owner, prefix, path))?;
    let mut bytes = Vec::new();
    let mut chunk = [0; 8192];
    let clean_eof = loop {
        // ifstream has no exception mask: read failure sets fail/eof state.
        let count = input
            .read_file_role(&mut chunk)
            .map(|(count, _)| count)
            .unwrap_or(0);
        bytes.extend_from_slice(&chunk[..count]);
        if count < chunk.len() {
            break input.eof().unwrap_or(false);
        }
    };
    let _ = input.close_file_role();
    Ok(ReadText { bytes, clean_eof })
}

pub(super) struct TextOutput {
    stream: OwnedCstdio,
    failed: bool,
    sticky: bool,
}
impl TextOutput {
    pub(super) fn open(
        path: &Path,
        owner: &IoErrorOwner,
        prefix: &[u8],
        sticky: bool,
    ) -> Result<Self, String> {
        let stream = OwnedCstdio::open_file_role_path(path, "w")
            .map_err(|_| path_error(owner, prefix, path))?;
        Ok(Self {
            stream,
            failed: false,
            sticky,
        })
    }
    pub(super) fn endl(&mut self) {
        self.write(b"\n");
        if !self.failed {
            self.failed = !self.stream.flush_text_output().unwrap_or(false);
        }
    }
    pub(super) fn write(&mut self, bytes: &[u8]) {
        if self.sticky && self.failed {
            return;
        }
        self.failed |= !self.stream.write_text_output(bytes).unwrap_or(false);
    }
}
impl Drop for TextOutput {
    fn drop(&mut self) {
        let _ = self.stream.close_file_role();
    }
}

/// `operator>> uint64_t` advances through the numeric prefix, preserves its
/// failbit, and wraps a representable negative magnitude. No UTF-8 conversion.
pub(super) struct Fields<'a> {
    bytes: &'a [u8],
    at: usize,
    failed: bool,
    clean_end: bool,
}
impl<'a> Fields<'a> {
    pub(super) fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            at: 0,
            failed: false,
            clean_end: true,
        }
    }
    pub(super) fn from_text(text: &'a ReadText) -> Self {
        Self {
            bytes: &text.bytes,
            at: 0,
            failed: false,
            clean_end: text.clean_eof,
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
    pub(super) fn word(&mut self) -> Option<&'a [u8]> {
        if self.failed {
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
        if begin == self.at || (self.at == self.bytes.len() && !self.clean_end) {
            self.failed = true;
            return None;
        }
        Some(&self.bytes[begin..self.at])
    }
    pub(super) fn number(&mut self) -> Option<u64> {
        if self.failed {
            return None;
        }
        self.space();
        let negative = self.bytes.get(self.at) == Some(&b'-');
        if matches!(self.bytes.get(self.at), Some(b'-' | b'+')) {
            self.at += 1;
        }
        let begin = self.at;
        let mut value = 0u64;
        let mut overflow = false;
        while let Some(&digit) = self.bytes.get(self.at).filter(|b| b.is_ascii_digit()) {
            match value
                .checked_mul(10)
                .and_then(|n| n.checked_add((digit - b'0') as u64))
            {
                Some(next) => value = next,
                None => overflow = true,
            }
            self.at += 1;
        }
        if self.at == begin || overflow || (self.at == self.bytes.len() && !self.clean_end) {
            self.failed = true;
            return None;
        }
        Some(if negative {
            value.wrapping_neg()
        } else {
            value
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let input = read_text(&path, &IoErrorOwner::default(), b"Error opening file ").unwrap();
        assert!(!input.clean_eof);
        assert!(input.bytes.is_empty());
        assert_eq!(input.lines().count(), 0);
        std::fs::remove_dir(path).unwrap();
    }

    #[test]
    fn terminal_read_failure_does_not_commit_partial_fields_or_lines() {
        let input = ReadText {
            bytes: b"1 2\n3 4".to_vec(),
            clean_eof: false,
        };
        let mut fields = Fields::from_text(&input);
        assert_eq!(fields.number(), Some(1));
        assert_eq!(fields.number(), Some(2));
        assert_eq!(fields.number(), Some(3));
        assert_eq!(fields.number(), None);
        assert_eq!(input.lines().collect::<Vec<_>>(), [b"1 2".as_slice()]);
        let token_input = ReadText {
            bytes: b"path".to_vec(),
            clean_eof: false,
        };
        let mut token = Fields::from_text(&token_input);
        assert_eq!(token.word(), None);
    }

    #[test]
    fn clean_eof_accepts_a_final_field_and_unterminated_line() {
        let input = ReadText {
            bytes: b"1 2\n3 4".to_vec(),
            clean_eof: true,
        };
        let mut fields = Fields::from_text(&input);
        assert_eq!((fields.number(), fields.number()), (Some(1), Some(2)));
        assert_eq!((fields.number(), fields.number()), (Some(3), Some(4)));
        assert_eq!(
            input.lines().collect::<Vec<_>>(),
            [b"1 2".as_slice(), b"3 4".as_slice()]
        );
        let completed = ReadText {
            bytes: b"5 6 ".to_vec(),
            clean_eof: false,
        };
        let mut fields = Fields::from_text(&completed);
        assert_eq!((fields.number(), fields.number()), (Some(5), Some(6)));
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
        assert_eq!(overflow.number(), None);
        assert_eq!(overflow.number(), None);
    }
}
