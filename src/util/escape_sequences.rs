use crate::util::text_buffer::{find_first_of, TextBuffer};
use std::sync::LazyLock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EscapeSequence<'a> {
    pub c: u8,
    pub seq: &'a str,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EscapeSequences {
    seq: [String; 256],
    bytes: [Vec<u8>; 256],
}

pub static XML: LazyLock<EscapeSequences> = LazyLock::new(EscapeSequences::xml);

impl EscapeSequences {
    pub fn new(seqs: &[EscapeSequence<'_>]) -> Self {
        let seq = std::array::from_fn(|i| (i as u8 as char).to_string());
        let bytes = std::array::from_fn(|i| vec![i as u8]);
        let mut out = Self { seq, bytes };
        for s in seqs {
            out.seq[s.c as usize] = s.seq.to_string();
            out.bytes[s.c as usize] = s.seq.as_bytes().to_vec();
        }
        out
    }

    pub fn xml() -> Self {
        Self::new(&[
            EscapeSequence {
                c: b'<',
                seq: "&lt;",
            },
            EscapeSequence {
                c: b'>',
                seq: "&gt;",
            },
            EscapeSequence {
                c: b'&',
                seq: "&amp;",
            },
            EscapeSequence {
                c: b'\'',
                seq: "&apos;",
            },
            EscapeSequence {
                c: b'"',
                seq: "&quot;",
            },
        ])
    }

    pub fn escape_bytes_into(&self, input: &[u8], out: &mut Vec<u8>) {
        for &byte in input {
            out.extend_from_slice(&self.bytes[byte as usize]);
        }
    }
    pub fn escape_bytes(&self, input: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(input.len());
        self.escape_bytes_into(input, &mut out);
        out
    }
    pub fn escape_char_bytes(&self, byte: u8) -> &[u8] {
        &self.bytes[byte as usize]
    }
    /// Legacy per-character String view. Actual Source bytes use
    /// escape_char_bytes; a lone non-UTF8 default byte has no &str view.
    pub fn escape_char(&self, c: u8) -> &str {
        &self.seq[c as usize]
    }

    pub fn escape_len(&self, s: &str, len: usize, out: &mut String) {
        let bytes = self.escape_bytes(&s.as_bytes()[..len]);
        // Input String and configured replacement Strings are UTF-8. The
        // canonical byte table preserves unchanged multibyte sequences.
        out.push_str(std::str::from_utf8(&bytes).expect("String escape output must remain UTF-8"));
    }

    pub fn escape_str_into(&self, s: &str, out: &mut String) {
        self.escape_len(s, s.len(), out);
    }

    pub fn try_escape_string(&self, s: &str) -> Result<String, std::string::FromUtf8Error> {
        String::from_utf8(self.escape_bytes(s.as_bytes()))
    }
    pub fn escape_string(&self, s: &str) -> String {
        self.try_escape_string(s)
            .expect("String escape output is not UTF-8; use escape_bytes or try_escape_string")
    }
}

pub fn print_escaped_until(
    buf: &mut TextBuffer,
    s: &str,
    delimiters: &str,
    esc: Option<&EscapeSequences>,
) {
    if let Some(esc) = esc {
        let mut tmp = String::new();
        esc.escape_len(s, find_first_of(s, delimiters), &mut tmp);
        buf.append_str(&tmp);
    } else {
        buf.write_until(s, delimiters);
    }
}

pub fn print_escaped(buf: &mut TextBuffer, s: &str, esc: Option<&EscapeSequences>) {
    if let Some(esc) = esc {
        let mut tmp = String::new();
        esc.escape_str_into(s, &mut tmp);
        buf.append_str(&tmp);
    } else {
        buf.append_str(s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_escape_sequences_xml() {
        let xml = EscapeSequences::xml();
        assert_eq!(xml.escape_char(b'<'), "&lt;");
        assert_eq!(xml.escape_char(b'A'), "A");
        assert_eq!(
            xml.escape_string("<tag a=\"b&c\">"),
            "&lt;tag a=&quot;b&amp;c&quot;&gt;"
        );
        assert_eq!(XML.escape_char(b'&'), "&amp;");
    }

    #[test]
    fn test_print_escaped() {
        let xml = EscapeSequences::xml();
        let mut b = TextBuffer::new();
        print_escaped_until(&mut b, "a<b,c", ",", Some(&xml));
        assert_eq!(String::from_utf8(b.data().to_vec()).unwrap(), "a&lt;b");

        b.clear();
        print_escaped_until(&mut b, "a<b,c", ",", None);
        assert_eq!(String::from_utf8(b.data().to_vec()).unwrap(), "a<b");

        b.clear();
        print_escaped(&mut b, "\"x\"", Some(&xml));
        assert_eq!(
            String::from_utf8(b.data().to_vec()).unwrap(),
            "&quot;x&quot;"
        );
    }
}
