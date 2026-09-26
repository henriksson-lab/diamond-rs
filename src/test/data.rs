//! Test sequence corpus from `diamond/src/test/data.cpp`.
//!
//! This is a test-only data translation unit with no executable routines. The
//! parser decodes its C++ string-literal initializer once, preserving the
//! source corpus without maintaining a second 100-KiB hand-copied table.

use std::sync::OnceLock;

static SEQS: OnceLock<Vec<(String, String)>> = OnceLock::new();

/// Returns the upstream `(identifier, amino-acid sequence)` corpus in source
/// order. This is the Rust equivalent of `Test::seqs`.
pub fn seqs() -> &'static [(String, String)] {
    SEQS.get_or_init(|| {
        parse_cpp_pairs(include_bytes!("../../diamond/src/test/data.cpp"))
            .expect("vendored test sequence initializer must remain valid")
    })
}

fn parse_cpp_pairs(source: &[u8]) -> Result<Vec<(String, String)>, String> {
    let declaration = b"const vector<pair<string, string>> seqs = {";
    let start = source
        .windows(declaration.len())
        .position(|window| window == declaration)
        .ok_or_else(|| "missing Test::seqs declaration".to_string())?
        + declaration.len();
    let end = source[start..]
        .windows(3)
        .position(|window| window == b"};\n")
        .ok_or_else(|| "missing Test::seqs terminator".to_string())?
        + start;
    let source = &source[start..end];

    let mut literals = Vec::new();
    let mut index = 0;
    while index < source.len() {
        if source[index] != b'"' {
            index += 1;
            continue;
        }
        index += 1;
        let mut value = Vec::new();
        while index < source.len() && source[index] != b'"' {
            if source[index] == b'\\' {
                if source.get(index + 1) == Some(&b'\n') {
                    index += 2;
                    continue;
                }
                if source.get(index + 1) == Some(&b'\r') && source.get(index + 2) == Some(&b'\n') {
                    index += 3;
                    continue;
                }
                let escaped = *source
                    .get(index + 1)
                    .ok_or_else(|| "unterminated string escape".to_string())?;
                value.push(match escaped {
                    b'n' => b'\n',
                    b'r' => b'\r',
                    b't' => b'\t',
                    b'\\' => b'\\',
                    b'"' => b'"',
                    other => other,
                });
                index += 2;
            } else {
                value.push(source[index]);
                index += 1;
            }
        }
        if source.get(index) != Some(&b'"') {
            return Err("unterminated string literal".to_string());
        }
        index += 1;
        literals.push(String::from_utf8(value).map_err(|error| error.to_string())?);
    }

    if literals.len() % 2 != 0 {
        return Err("unpaired identifier/sequence literal".to_string());
    }
    Ok(literals
        .chunks_exact(2)
        .map(|pair| (pair[0].clone(), pair[1].clone()))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entire_cpp_fixture_is_decoded_as_identifier_sequence_pairs() {
        let data = seqs();
        assert!(
            data.len() > 100,
            "unexpectedly small corpus: {}",
            data.len()
        );
        assert_eq!(data.first().unwrap().0, "d2dc3a_");
        assert!(data.first().unwrap().1.starts_with("eelseaerkavqamwarly"));
        assert!(data.iter().all(|(id, sequence)| !id.is_empty()
            && !sequence.is_empty()
            && sequence.bytes().all(|byte| byte.is_ascii_alphabetic())));
        let unique: std::collections::HashSet<_> = data.iter().map(|pair| &pair.0).collect();
        assert_eq!(unique.len(), data.len());
    }

    #[test]
    fn parser_handles_cpp_line_splicing_and_escapes() {
        let fixture = b"const vector<pair<string, string>> seqs = {{\"id\", \"AB\\\nCD\"}, {\"x\", \"Y\\tZ\"}};\n";
        assert_eq!(
            parse_cpp_pairs(fixture).unwrap(),
            vec![
                ("id".to_string(), "ABCD".to_string()),
                ("x".to_string(), "Y\tZ".to_string())
            ]
        );
    }
}
