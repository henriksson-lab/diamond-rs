//! DIAMOND database format.

pub mod dmnd;

// Preserve the historical `crate::data::dmnd::ReferenceHeader` API after
// introducing the source-mirrored `data::dmnd::dmnd` hierarchy.
pub use dmnd::*;
