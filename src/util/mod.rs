pub mod algo;
pub mod binary_buffer;
pub mod command_line_parser;
#[doc(hidden)]
pub mod compat_rng;
pub mod data_structures;
pub mod endianness;
pub mod enum_utils;
pub mod escape_sequences;
pub mod geo;
pub mod hash;
pub mod hsp;
pub mod interval;
pub mod intrin;
pub mod io;
pub mod kmer;
pub mod log_erfc;
pub mod log_stream;
pub mod math;
pub mod memory;
pub mod misc;
pub mod multiprocessing;
pub mod optional;
pub mod options;
pub mod parallel;
pub mod profiler;
pub mod range;
pub mod scores;
pub mod search;
pub mod seed_table;
pub mod sequence;
pub mod simd;
pub mod string;
pub mod system;
pub mod table;
pub mod text_buffer;
pub mod tsv;
pub mod util;

/// Filesystem-safe, deterministic discriminator for paths created by tests.
///
/// Rust names test threads after the fully qualified test path. Those names
/// contain `::`, which is not a valid Windows filename component. Hashing the
/// name also keeps paths comfortably below the legacy Windows path limit.
#[cfg(test)]
pub(crate) fn test_path_component() -> String {
    use std::hash::{Hash, Hasher};

    let thread = std::thread::current();
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    match thread.name() {
        Some(name) => name.hash(&mut hasher),
        None => format!("{:?}", thread.id()).hash(&mut hasher),
    }
    format!("{:016x}", hasher.finish())
}
