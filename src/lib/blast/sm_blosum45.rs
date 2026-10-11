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
