//! HSP operations translated from `diamond/src/basic/hssp.cpp`.
//!
//! The core types are owned by `align::hsp` because alignment and output code
//! share them. This mirrored module preserves the upstream source hierarchy
//! and provides the canonical compatibility surface for the translated file.

pub use crate::align::hsp::{normalized_range, Hsp, HspContext, HspContextIterator, HspIterator};
