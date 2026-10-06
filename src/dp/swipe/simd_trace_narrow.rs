//! Narrow AVX2 traceback tiers used by SWIPE score bins 0 and 1.
//!
//! The implementation lives in `simd_trace_upstream`: it is a structural
//! port of upstream's moving band, vector profile, rolling DP rows, and packed
//! traceback masks. This module retains the dispatch surface used by SWIPE.

use super::simd_trace::TraceTarget;
use crate::basic::value::Letter;
use crate::dp::smith_waterman::SwResult;
use crate::stats::score_matrix::ScoreMatrix;

pub struct NarrowTraceBatch {
    pub results: Vec<SwResult>,
    pub overflow_mask: u32,
}

pub fn trace_batch_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    _score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    _semi_global: bool,
) -> Option<NarrowTraceBatch> {
    if targets.is_empty()
        || targets.len() > 32
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
    {
        return None;
    }
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 is runtime-detected and the kernel bounds every lane.
        let (results, overflow_mask) = unsafe {
            super::simd_trace_upstream::trace_i8(
                query,
                targets,
                _score_matrix,
                query_cbs,
                _semi_global,
            )
        };
        return Some(NarrowTraceBatch {
            results,
            overflow_mask,
        });
    }
    None
}

pub fn trace_batch_i16(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    _score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    _semi_global: bool,
) -> Option<NarrowTraceBatch> {
    if targets.is_empty()
        || targets.len() > 16
        || (!query_cbs.is_empty() && query_cbs.len() < query.len())
    {
        return None;
    }
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 is runtime-detected and the kernel bounds every lane.
        let (results, overflow_mask) = unsafe {
            super::simd_trace_upstream::trace_i16(
                query,
                targets,
                _score_matrix,
                query_cbs,
                _semi_global,
            )
        };
        return Some(NarrowTraceBatch {
            results,
            overflow_mask,
        });
    }
    None
}
