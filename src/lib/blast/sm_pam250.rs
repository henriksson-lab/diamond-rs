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
