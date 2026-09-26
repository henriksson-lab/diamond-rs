//! DAA output writers.
//!
//! This mirrors the upstream `src/output/daa` source hierarchy.  The wire
//! format itself remains owned by [`crate::data::daa`], while this module is
//! the output-facing compatibility layer corresponding to DIAMOND's
//! `daa_write.cpp`.

pub mod daa_record;
pub mod daa_write;
pub mod merge;
pub mod view;
