//! Rust representation of `diamond/src/lib/blast/sm_blosum90.c`.

use super::sm_blosum45::PackedScoreMatrix;
use crate::stats::matrices::BLOSUM90;

pub const BLOSUM90_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Translation of the exported C descriptor `NCBISM_Blosum90`.
pub const NCBISM_BLOSUM90: PackedScoreMatrix = PackedScoreMatrix {
    symbols: BLOSUM90_SYMBOLS,
    matrix: &BLOSUM90,
    default_score: -6,
};
