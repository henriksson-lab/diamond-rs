//! Compatibility path for the original pre-hierarchy Rust port.
//!
//! The upstream source lives at `util/parallel/multiprocessing.{h,cpp}`, so
//! new code should use [`crate::util::parallel::multiprocessing`].

pub use crate::util::parallel::multiprocessing::*;
