//! Rust representation of `diamond/src/lib/blast/sm_blosum80.c`.

use super::sm_blosum45::PackedScoreMatrix;
use crate::stats::matrices::BLOSUM80;

pub const BLOSUM80_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Translation of the exported C descriptor `NCBISM_Blosum80`.
pub const NCBISM_BLOSUM80: PackedScoreMatrix = PackedScoreMatrix {
    symbols: BLOSUM80_SYMBOLS,
    matrix: &BLOSUM80,
    default_score: -6,
};

#[cfg(test)]
mod tests {
    use super::*;

    fn upstream_scores() -> Vec<i8> {
        let source = include_str!("../../../diamond/src/lib/blast/sm_blosum80.c");
        let initializer = source
            .split("s_Blosum80PSM[25 * 25] = {")
            .nth(1)
            .expect("upstream score declaration")
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

    #[test]
    fn descriptor_and_every_packed_score_match_upstream() {
        assert_eq!(NCBISM_BLOSUM80.symbols, b"ARNDCQEGHILKMFPSTWYVBJZX*");
        assert_eq!(NCBISM_BLOSUM80.default_score, -6);
        let expected = upstream_scores();
        assert_eq!(expected.len(), 25 * 25);
        for row in 0..25 {
            for column in 0..25 {
                assert_eq!(
                    NCBISM_BLOSUM80.score(row, column),
                    Some(expected[row * 25 + column]),
                    "score mismatch at ({row}, {column})"
                );
            }
        }
    }
}
