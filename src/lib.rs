//! # diamond-rs
//!
//! Rust port of the DIAMOND protein sequence aligner.
//!
//! DIAMOND is a high-performance sequence aligner for protein and translated
//! DNA searches, designed for big sequence data analysis. This crate provides
//! both a library API and a CLI binary.
//!
//! ## Library Usage
//!
//! ```rust,no_run
//! use diamond::prelude::*;
//!
//! // Parse FASTA sequences
//! let records = diamond::data::fasta::read_fasta_amino_acid(
//!     b">seq1\nARNDCQEGHILKMFPSTWYV\n" as &[u8]
//! ).unwrap();
//!
//! // Create a scoring matrix
//! let score_matrix = diamond::stats::score_matrix::ScoreMatrix::new(
//!     "blosum62", 11, 1, 0, 1, 0
//! ).unwrap();
//!
//! // Run Smith-Waterman alignment
//! let result = diamond::dp::smith_waterman::smith_waterman(
//!     &records[0].sequence,
//!     &records[0].sequence,
//!     &score_matrix,
//! );
//! assert!(result.score > 0);
//! ```

pub mod align;
#[path = "lib/alp/mod.rs"]
pub mod alp;
pub mod basic;
#[path = "lib/blast/mod.rs"]
pub mod blast;
pub mod chaining;
pub mod cluster;
pub mod commands;
pub mod config;
pub mod data;
pub mod dna;
pub mod dp;
pub mod ffi;
pub mod masking;
pub mod mcl;
#[path = "lib/murmurhash/mod.rs"]
pub mod murmurhash;
pub mod output;
pub mod run;
pub mod search;
pub mod stats;
#[path = "lib/tantan/mod.rs"]
pub mod tantan;
pub mod tools;
// Translation-audit fixtures read the pinned upstream C/C++ sources. Keep
// them out of normal library builds so a published crate remains buildable
// without a sibling `diamond/` checkout; CI checks them in unit-test builds.
#[cfg(test)]
#[path = "test/mod.rs"]
mod upstream_test;
pub mod util;

/// Convenient re-exports for common types.
pub mod prelude {
    pub use crate::align::hsp::{Hsp, Match};
    pub use crate::basic::sequence::Sequence;
    pub use crate::basic::value::{Letter, Score, SequenceType};
    pub use crate::data::fasta::FastaRecord;
    pub use crate::dp::smith_waterman::SwResult;
    pub use crate::dp::ungapped::DiagonalSegment;
    pub use crate::output::format::FieldId;
    pub use crate::stats::score_matrix::ScoreMatrix;
}

/// Run DIAMOND with command-line arguments through the C++ FFI implementation.
///
/// The FFI backend is intended for non-Windows conformance testing and is only
/// available when the crate is built with `--features ffi`.
///
/// # Example
/// ```rust,no_run,ignore
/// let code = diamond::run(&["diamond", "version"]);
/// assert_eq!(code, 0);
/// ```
#[cfg(all(feature = "ffi", not(windows)))]
pub fn run(args: &[&str]) -> i32 {
    crate::ffi::run(args)
}

/// Fallback when the C++ FFI backend is not compiled.
#[cfg(not(all(feature = "ffi", not(windows))))]
pub fn run(_args: &[&str]) -> i32 {
    eprintln!("C++ FFI backend is not available in this build.");
    1
}
