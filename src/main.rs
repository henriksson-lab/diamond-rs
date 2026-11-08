use diamond::commands::dispatch::{self, NativeDispatchError};

#[cfg(unix)]
mod source_signals;

fn main() {
    // Source inherits its parent's SIGPIPE policy. Rust's runtime changes it
    // before main, so the executable captures and restores the original value.
    #[cfg(unix)]
    source_signals::restore();
    let args: Vec<std::ffi::OsString> = std::env::args_os().collect();
    // SAFETY: this executable owns the process standard descriptors. Library
    // run APIs never enable this policy and preserve their host's descriptors.
    unsafe {
        diamond::util::io::file_sink::enable_process_stdout_ownership();
        diamond::util::io::file_source::enable_process_stdin_ownership();
    }
    let code = match dispatch::try_run_os(&args) {
        Ok(code) => code,
        Err(NativeDispatchError::Io(error)) => dispatch::report_io_error(&error),
        // The optional C++ build is an explicitly selected comparison oracle.
        // Untranslated production requests must never silently execute it.
        Err(NativeDispatchError::Untranslated) if args.iter().any(|arg| arg == "--legacy") => {
            let filtered: Vec<&std::ffi::OsStr> = args
                .iter()
                .map(std::ffi::OsString::as_os_str)
                .filter(|value| *value != "--legacy")
                .collect();
            run_legacy_os(&filtered)
        }
        Err(NativeDispatchError::Untranslated) => {
            diamond::util::log_stream::flush_tied_stdout();
            diamond::util::log_stream::write_source_error_line(format_args!(
                "Error: Native execution of this command or option is not yet connected."
            ));
            1
        }
    };
    #[cfg(windows)]
    {
        use std::io::Write;
        // Rust exits Windows through ExitProcess, which skips CRT exit
        // handlers. Source returns through the C runtime, flushing its actual
        // stdout FILE (which FileSink deliberately does not close on MSVC).
        // Flush Rust's separate writers before entering that same CRT exit.
        let _ = std::io::stdout().flush();
        let _ = std::io::stderr().flush();
        unsafe { libc::exit(code) }
    }
    #[cfg(not(windows))]
    std::process::exit(code);
}

#[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
fn run_legacy_os(args: &[&std::ffi::OsStr]) -> i32 {
    let values = args
        .iter()
        .map(|value| diamond::util::io::source_os_bytes(value))
        .collect::<diamond::util::io::IoResult<Vec<_>>>();
    match values {
        Ok(values) => {
            let views: Vec<&[u8]> = values.iter().map(|value| value.as_ref()).collect();
            diamond::ffi::run_cpp_bytes(&views)
        }
        Err(error) => dispatch::report_io_error(&std::io::Error::other(error)),
    }
}

#[cfg(not(all(feature = "ffi", diamond_ffi_available, not(windows))))]
fn run_legacy_os(_args: &[&std::ffi::OsStr]) -> i32 {
    diamond::util::log_stream::flush_tied_stdout();
    diamond::util::log_stream::write_source_error_line(format_args!(
        "Error: C++ FFI fallback is not available in this build."
    ));
    diamond::util::log_stream::write_source_error_line(format_args!(
        "A non-Windows build with the `ffi` feature and C++ oracle sources is required to enable it."
    ));
    1
}
