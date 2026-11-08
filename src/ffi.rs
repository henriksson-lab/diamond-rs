pub mod ffi;

#[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
use std::ffi::CString;
#[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
use std::os::raw::{c_char, c_int};

#[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
extern "C" {
    fn diamond_main(argc: c_int, argv: *const *const c_char) -> c_int;
}

/// Run DIAMOND with the given command-line arguments.
///
/// The first argument should be the program name (e.g., "diamond").
/// Returns the exit code from DIAMOND.
#[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
pub fn run_cpp(args: &[&str]) -> i32 {
    let args: Vec<&[u8]> = args.iter().map(|value| value.as_bytes()).collect();
    run_cpp_bytes(&args)
}

/// Raw argument entry for the explicit, default-off conformance oracle.
#[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
pub fn run_cpp_bytes(args: &[&[u8]]) -> i32 {
    let c_strings: Vec<CString> = args
        .iter()
        .map(|value| CString::new(*value).expect("argument contains null byte"))
        .collect();
    let argc = c_int::try_from(c_strings.len()).expect("too many arguments for C++ adapter");
    let mut c_ptrs: Vec<*const c_char> = c_strings.iter().map(|value| value.as_ptr()).collect();
    // C main's argv contract includes argv[argc] == nullptr.
    c_ptrs.push(std::ptr::null());
    unsafe { diamond_main(argc, c_ptrs.as_ptr()) }
}

/// Backward-compatible name for the optional C++ conformance adapter.
#[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
pub fn run(args: &[&str]) -> i32 {
    run_cpp(args)
}
