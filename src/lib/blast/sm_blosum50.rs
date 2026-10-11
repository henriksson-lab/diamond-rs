//! Rust representation of the packed BLOSUM50 table from
//! `diamond/src/lib/blast/sm_blosum50.c`.
//!
//! DIAMOND already stores the same scores in its canonical 26-symbol
//! [`StandardMatrix`](crate::stats::standard_matrix::StandardMatrix). This
//! facade preserves the vendored NCBI table's 25-symbol view without keeping
//! a second 625-entry copy of the data.

use super::sm_blosum45::PackedScoreMatrix;
use crate::stats::matrices::BLOSUM50;

/// Residue order used by the packed NCBI score matrix.
pub const BLOSUM50_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Translation of the exported C descriptor `NCBISM_Blosum50`.
pub const NCBISM_BLOSUM50: PackedScoreMatrix = PackedScoreMatrix {
    symbols: BLOSUM50_SYMBOLS,
    matrix: &BLOSUM50,
    default_score: -5,
};
