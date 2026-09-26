//! Compatibility facade for declarations and definitions in `util/util.cpp`
//! and `util/util.h`.
//!
//! Implementations remain grouped in the established Rust utility modules;
//! this mirrored path exposes the upstream file's snake-case surface.

pub use super::log_stream::exit_with_error;
pub use super::misc::*;
