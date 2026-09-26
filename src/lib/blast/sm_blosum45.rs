//! Rust representation of the packed BLOSUM45 table from
//! `diamond/src/lib/blast/sm_blosum45.c`.
//!
//! DIAMOND already stores the same scores in its canonical 26-symbol
//! [`StandardMatrix`](crate::stats::standard_matrix::StandardMatrix).  This
//! facade preserves the vendored NCBI table's 25-symbol view without keeping
//! a second 625-entry copy of the data.

use crate::stats::{matrices::BLOSUM45, standard_matrix::StandardMatrix};

/// Residue order used by the packed NCBI score matrix.
pub const BLOSUM45_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Safe equivalent of `SNCBIPackedScoreMatrix` for an existing DIAMOND
/// standard matrix.
#[derive(Clone, Copy)]
pub struct PackedScoreMatrix {
    pub symbols: &'static [u8],
    pub matrix: &'static StandardMatrix,
    pub default_score: i8,
}

impl PackedScoreMatrix {
    /// Return an entry by packed-table index.
    pub fn score(&self, row: usize, column: usize) -> Option<i8> {
        if row >= self.symbols.len() || column >= self.symbols.len() {
            return None;
        }
        Some(self.matrix.scores[row * 26 + column])
    }
}

/// Translation of the exported C descriptor `NCBISM_Blosum45`.
pub const NCBISM_BLOSUM45: PackedScoreMatrix = PackedScoreMatrix {
    symbols: BLOSUM45_SYMBOLS,
    matrix: &BLOSUM45,
    default_score: -5,
};

#[cfg(test)]
mod tests {
    use super::*;

    fn upstream_scores() -> Vec<i8> {
        let source = include_str!("../../../diamond/src/lib/blast/sm_blosum45.c");
        let declaration = source
            .split("s_Blosum45PSM[25 * 25] = {")
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

    #[test]
    fn descriptor_matches_upstream() {
        assert_eq!(NCBISM_BLOSUM45.symbols, b"ARNDCQEGHILKMFPSTWYVBJZX*");
        assert_eq!(NCBISM_BLOSUM45.default_score, -5);
    }

    #[test]
    fn every_packed_score_matches_upstream() {
        let expected = upstream_scores();
        assert_eq!(expected.len(), 25 * 25);
        for row in 0..25 {
            for column in 0..25 {
                assert_eq!(
                    NCBISM_BLOSUM45.score(row, column),
                    Some(expected[row * 25 + column]),
                    "score mismatch at ({row}, {column})"
                );
            }
        }
        assert_eq!(NCBISM_BLOSUM45.score(25, 0), None);
        assert_eq!(NCBISM_BLOSUM45.score(0, 25), None);
    }
}
