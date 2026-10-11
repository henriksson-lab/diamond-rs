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
