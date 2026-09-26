//! Mirrored module for `diamond/src/align/gapped_score.cpp`.
//!
//! The public implementations remain in `target` for API compatibility; this
//! module restores the upstream source hierarchy and its systematic names.

pub use super::target::align_work_targets as align;
pub use super::target::{add_dp_targets, align_work_targets, band, hsp_band};
