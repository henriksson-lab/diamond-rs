//! Rust representation of `diamond/src/lib/blast/sm_pam30.c`.

use super::sm_blosum45::PackedScoreMatrix;
use crate::stats::matrices::PAM30;

pub const PAM30_SYMBOLS: &[u8; 25] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Translation of the exported C descriptor `NCBISM_Pam30`.
pub const NCBISM_PAM30: PackedScoreMatrix = PackedScoreMatrix {
    symbols: PAM30_SYMBOLS,
    matrix: &PAM30,
    default_score: -17,
};
