//! Mirrored translation facade for `diamond/src/stats/comp_based_stats.cpp`.
//!
//! The implementations predate this hierarchy-preserving module and remain in
//! [`super::cbs`], where their callers already import them. Re-exporting the
//! audited surface here preserves those callers while matching the upstream
//! source layout.

pub use super::cbs::{
    blast_composition_based_stats, blast_karlin_lambda_nr, calc_avg_score, calc_lambda,
    calc_x_score, composition_based_stats, freq_ratio_to_score, get_score_range, ideal_lambda,
    matrix_score_probs, round_score_matrix, scale_square_matrix, set_xuo_scores, BlastScoreFreq,
    ALPH_TO_NCBI, NCBI_ALPH,
};
