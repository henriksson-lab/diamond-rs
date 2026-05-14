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
}

pub static XML: LazyLock<EscapeSequences> = LazyLock::new(EscapeSequences::xml);

impl EscapeSequences {
    pub fn new(seqs: &[EscapeSequence<'_>]) -> Self {
        let seq = std::array::from_fn(|i| (i as u8 as char).to_string());
        let mut out = Self { seq };
        for s in seqs {
            out.seq[s.c as usize] = s.seq.to_string();
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

    pub fn escape_char(&self, c: u8) -> &str {
        &self.seq[c as usize]
    }

    pub fn escape_len(&self, s: &str, len: usize, out: &mut String) {
        out.reserve(len);
        for &b in &s.as_bytes()[..len] {
            out.push_str(self.escape_char(b));
        }
    }

    pub fn escape_str_into(&self, s: &str, out: &mut String) {
        self.escape_len(s, s.len(), out);
    }

    pub fn escape_string(&self, s: &str) -> String {
        let mut out = String::new();
        self.escape_str_into(s, &mut out);
        out
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
