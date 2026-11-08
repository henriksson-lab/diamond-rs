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
#[path = "test/mod.rs"]
pub mod upstream_test;
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

/// Run connected DIAMOND commands through the native Rust dispatcher.
///
/// Returns the command status without terminating the calling process. The
/// first argument is the program name. Untranslated commands/options return
/// an explicit error; C++ conformance is available separately via `ffi::run_cpp`
/// only in a feature-enabled source checkout.
///
/// ```no_run
/// assert_eq!(diamond::run(&["diamond", "version"]), 0);
/// ```
pub fn run(args: &[&str]) -> i32 {
    crate::commands::dispatch::run(args)
}

/// Run the native CLI with OS arguments, preserving Unix filename bytes and
/// Windows UTF-16 paths. No production route invokes the C++ oracle.
pub fn run_os(args: &[&std::ffi::OsStr]) -> i32 {
    let args: Vec<std::ffi::OsString> = args.iter().map(|value| (*value).to_owned()).collect();
    crate::commands::dispatch::run_os(&args)
}
