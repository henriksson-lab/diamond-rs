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

/// Pick a filesystem-supported test name without treating opaque path bytes
/// as universally creatable on Unix (APFS rejects invalid UTF-8 with EILSEQ).
/// The rejection path still exercises the real Source opener and its saved
/// raw-path diagnostic; callers then run their complete payload workflow with
/// an explicit valid-Unicode alternative.
#[cfg(all(test, unix))]
pub(crate) fn test_os_path_or_unicode(
    directory: &std::path::Path,
    raw_component: &[u8],
    unicode_component: &str,
) -> (std::path::PathBuf, bool) {
    use std::os::unix::ffi::{OsStrExt, OsStringExt};
    let path = directory.join(std::ffi::OsString::from_vec(raw_component.to_vec()));
    match std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
    {
        Ok(file) => {
            drop(file);
            std::fs::remove_file(&path).unwrap();
            (path, true)
        }
        Err(error) if error.raw_os_error() == Some(libc::EILSEQ) => {
            let source_error = io::FileSink::new_path(&path, "wb", false, 64).unwrap_err();
            let saved = match source_error {
                io::IoError::FileOpenDetails(saved) => saved,
                other => panic!("Expected saved Source open failure, got {other:?}"),
            };
            assert_eq!(saved.path, path);
            assert_eq!(saved.raw_os_error, Some(libc::EILSEQ));
            assert!(!saved.reason.is_empty());
            let mut expected = b"Error opening file ".to_vec();
            expected.extend_from_slice(path.as_os_str().as_bytes());
            assert_eq!(saved.try_message_bytes().unwrap(), expected);
            assert!(!path.exists());
            (directory.join(unicode_component), false)
        }
        Err(error) => panic!("Unexpected raw-filename capability error for {path:?}: {error}"),
    }
}
