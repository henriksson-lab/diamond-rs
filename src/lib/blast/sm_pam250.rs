//! Rust representation of the packed PAM250 table from
//! `diamond/src/lib/blast/sm_pam250.c`.
//!
//! DIAMOND already stores the same scores in its canonical 26-symbol
//! [`StandardMatrix`](crate::stats::standard_matrix::StandardMatrix). This
//! facade preserves the vendored NCBI table's 25-symbol view without keeping
//! a second 625-entry copy of the data.

use super::sm_blosum45::PackedScoreMatrix;
use crate::stats::matrices::PAM250;

/// Residue order used by the packed NCBI score matrix.
pub const PAM250_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Translation of the exported C descriptor `NCBISM_Pam250`.
pub const NCBISM_PAM250: PackedScoreMatrix = PackedScoreMatrix {
    symbols: PAM250_SYMBOLS,
    matrix: &PAM250,
    default_score: -8,
};

#[cfg(test)]
mod tests {
    use super::*;

    const SOURCE: &str = include_str!("../../../diamond/src/lib/blast/sm_pam250.c");

    fn upstream_scores() -> Vec<i8> {
        let declaration = SOURCE
            .split("s_Pam250PSM[25 * 25] = {")
            .nth(1)
            .expect("upstream score declaration");
        let initializer = declaration
            .split("};")
            .next()
            .expect("upstream score initializer");

        let mut without_comments = String::with_capacity(initializer.len());
        let mut rest = initializer;
        while let Some(start) = rest.find("/*") {
            without_comments.push_str(&rest[..start]);
            let after_start = &rest[start + 2..];
            let end = after_start.find("*/").expect("closed block comment");
            rest = &after_start[end + 2..];
        }
        without_comments.push_str(rest);

        without_comments
            .split(',')
            .filter_map(|value| {
                let value = value.trim();
                (!value.is_empty()).then(|| value.parse::<i8>().expect("integer score"))
            })
            .collect()
    }

    fn upstream_descriptor() -> (&'static [u8], i8) {
        let initializer = SOURCE
            .split("const SNCBIPackedScoreMatrix NCBISM_Pam250 = {")
            .nth(1)
            .expect("upstream descriptor")
            .split("};")
            .next()
            .expect("upstream descriptor terminator");
        let quote = initializer.find('"').expect("symbol string");
        let after_quote = &initializer[quote + 1..];
        let end_quote = after_quote.find('"').expect("closed symbol string");
        let symbols = &after_quote.as_bytes()[..end_quote];
        let default_score = initializer
            .rsplit(',')
            .next()
            .expect("default score")
            .trim()
            .parse::<i8>()
            .expect("integer default score");
        (symbols, default_score)
    }

    #[test]
    fn descriptor_metadata_matches_upstream() {
        let (symbols, default_score) = upstream_descriptor();
        assert_eq!(NCBISM_PAM250.symbols, symbols);
        assert_eq!(NCBISM_PAM250.default_score, default_score);
    }

    #[test]
    fn every_packed_score_matches_upstream() {
        let expected = upstream_scores();
        assert_eq!(expected.len(), 25 * 25);
        for row in 0..25 {
            for column in 0..25 {
                assert_eq!(
                    NCBISM_PAM250.score(row, column),
                    Some(expected[row * 25 + column]),
                    "score mismatch at ({row}, {column})"
                );
            }
        }
        assert_eq!(NCBISM_PAM250.score(25, 0), None);
        assert_eq!(NCBISM_PAM250.score(0, 25), None);
    }
}
