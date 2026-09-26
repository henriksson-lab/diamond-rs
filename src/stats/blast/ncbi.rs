//! NCBI target-frequency optimization translated from
//! `diamond/src/stats/blast/ncbi.cpp`.
//!
//! The implementation is generic over alphabet size and is shared with the
//! composition-adjustment layer. Rust ownership replaces the source file's
//! explicit `ReNewtonSystemNew`/`ReNewtonSystemFree` allocation pair.

pub use crate::stats::target_freq::{
    blast_optimize_target_frequencies, calculate_residuals, compute_scores_from_probs,
    dual_residuals, evaluate_re_functions, factor_re_newton_system, multiply_by_a,
    multiply_by_a_transpose, residuals_linear_constraints, scaled_symmetric_product_a,
    solve_re_newton_system, ReNewtonSystem,
};
