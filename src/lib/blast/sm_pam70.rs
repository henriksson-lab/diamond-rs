//! Rust representation of `diamond/src/lib/blast/sm_pam70.c`.

use super::sm_blosum45::PackedScoreMatrix;
use crate::stats::matrices::PAM70;

pub const PAM70_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Translation of the exported C descriptor `NCBISM_Pam70`.
pub const NCBISM_PAM70: PackedScoreMatrix = PackedScoreMatrix {
    symbols: PAM70_SYMBOLS,
    matrix: &PAM70,
    default_score: -11,
};
