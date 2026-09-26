pub mod ffi;

#[cfg(all(feature = "ffi", not(windows)))]
use std::ffi::CString;
#[cfg(all(feature = "ffi", not(windows)))]
use std::os::raw::{c_char, c_int};

#[cfg(all(feature = "ffi", not(windows)))]
extern "C" {
    fn diamond_main(argc: c_int, argv: *const *const c_char) -> c_int;
}

/// Run DIAMOND with the given command-line arguments.
///
/// The first argument should be the program name (e.g., "diamond").
/// Returns the exit code from DIAMOND.
#[cfg(all(feature = "ffi", not(windows)))]
pub fn run_cpp(args: &[&str]) -> i32 {
    let c_strings: Vec<CString> = args
        .iter()
        .map(|s| CString::new(*s).expect("argument contains null byte"))
        .collect();
    let c_ptrs: Vec<*const c_char> = c_strings.iter().map(|s| s.as_ptr()).collect();

    unsafe { diamond_main(c_ptrs.len() as c_int, c_ptrs.as_ptr()) }
}

/// Backward-compatible name for the optional C++ conformance adapter.
#[cfg(all(feature = "ffi", not(windows)))]
pub fn run(args: &[&str]) -> i32 {
    run_cpp(args)
}
