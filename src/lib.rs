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
pub mod basic;
pub mod cluster;
pub mod commands;
pub mod config;
pub mod data;
pub mod dp;
pub mod ffi;
pub mod masking;
pub mod output;
pub mod search;
pub mod stats;
pub mod util;

/// Convenient re-exports for common types.
pub mod prelude {
    pub use crate::basic::value::{Letter, Score, SequenceType};
    pub use crate::basic::sequence::Sequence;
    pub use crate::data::fasta::FastaRecord;
    pub use crate::stats::score_matrix::ScoreMatrix;
    pub use crate::dp::smith_waterman::SwResult;
    pub use crate::dp::ungapped::DiagonalSegment;
    pub use crate::align::hsp::{Hsp, Match};
    pub use crate::output::format::FieldId;
}

/// Run DIAMOND with command-line arguments (FFI to C++ implementation).
///
/// This is the primary entry point that delegates to the C++ DIAMOND
/// implementation for full feature support and bit-exact output.
///
/// # Example
/// ```rust,no_run
/// let code = diamond::run(&["diamond", "version"]);
/// assert_eq!(code, 0);
/// ```
pub fn run(args: &[&str]) -> i32 {
    ffi::run(args)
}
