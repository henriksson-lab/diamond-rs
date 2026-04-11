use std::ffi::CString;
use std::os::raw::{c_char, c_int};

extern "C" {
    fn diamond_main(argc: c_int, argv: *const *const c_char) -> c_int;
}

/// Run DIAMOND with the given command-line arguments.
///
/// The first argument should be the program name (e.g., "diamond").
/// Returns the exit code from DIAMOND.
pub fn run(args: &[&str]) -> i32 {
    let c_strings: Vec<CString> = args
        .iter()
        .map(|s| CString::new(*s).expect("argument contains null byte"))
        .collect();
    let c_ptrs: Vec<*const c_char> = c_strings.iter().map(|s| s.as_ptr()).collect();

    unsafe { diamond_main(c_ptrs.len() as c_int, c_ptrs.as_ptr()) }
}
