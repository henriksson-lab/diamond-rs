//! Rust representation of `diamond/src/lib/blast/sm_blosum62.c`.

use super::sm_blosum45::PackedScoreMatrix;
use crate::stats::matrices::BLOSUM62;

/// Residue order used by the packed NCBI BLOSUM62 table.
pub const BLOSUM62_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Translation of the exported C descriptor `NCBISM_Blosum62`.
pub const NCBISM_BLOSUM62: PackedScoreMatrix = PackedScoreMatrix {
    symbols: BLOSUM62_SYMBOLS,
    matrix: &BLOSUM62,
    default_score: -4,
};
