//! Source console transport, distinct from FileSink and arbitrary user writers.
//! Supported platforms borrow their actual standard FILEs, preserving CRT
//! text mode, cout ties and cerr unit-buffering without rewriting bytes.
//! Other platforms retain the Rust stream fallback pending audited bindings.
use std::fmt;
use std::io::{self, Write};

#[derive(Clone, Copy)]
enum Kind {
    Out,
    Err,
}

pub(crate) struct ConsoleWriter {
    kind: Kind,
    // NativeMain's generic stdout consumer is Source's version << endl only.
    // Other stdout text helpers deliberately do not flush on newline bytes.
    #[cfg_attr(not(windows), allow(dead_code))]
    format_endl: bool,
    // Only bind this adapter where every formatted newline is a proven
    // Source std::endl (data write_all calls are deliberately unaffected).
    newline_endls: bool,
}

pub(crate) fn stdout_version() -> ConsoleWriter {
    ConsoleWriter {
        kind: Kind::Out,
        format_endl: true,
        newline_endls: false,
    }
}
pub(crate) fn stdout() -> ConsoleWriter {
    ConsoleWriter {
        kind: Kind::Out,
        format_endl: false,
        newline_endls: false,
    }
}
pub(crate) fn stderr() -> ConsoleWriter {
    ConsoleWriter {
        kind: Kind::Err,
        format_endl: false,
        newline_endls: false,
    }
}

/// Adapter for generic row renderers whose formatted newlines correspond
/// exactly to Source std::endl. Non-newline formatted fragments do not flush.
pub(crate) fn stdout_endl() -> ConsoleWriter {
    ConsoleWriter {
        kind: Kind::Out,
        format_endl: false,
        newline_endls: true,
    }
}

/// C printf output has no ostream failbit or cout tie. Its callers ignore the
/// printf return value, and never request an implicit per-newline flush.
pub(crate) struct PrintfStdout;
pub(crate) fn printf_stdout() -> PrintfStdout {
    PrintfStdout
}
impl Write for PrintfStdout {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        // These adapters carry formatted C text (the benchmark's numeric
        // fragments), not binary output. Precision bounds the non-NUL-terminated
        // Rust slice; printf retains its actual C stdout/text-mode buffering.
        for chunk in bytes.chunks(i32::MAX as usize) {
            unsafe {
                libc::printf(
                    c"%.*s".as_ptr(),
                    chunk.len() as libc::c_int,
                    chunk.as_ptr().cast::<libc::c_char>(),
                );
            }
        }
        Ok(bytes.len())
    }
    fn write_fmt(&mut self, args: fmt::Arguments<'_>) -> io::Result<()> {
        let text = args.to_string();
        // Keep one formatted Source printf invocation together, rather than
        // letting Rust's formatter split a value and its trailing separator.
        self.write_all(text.as_bytes())
    }
    fn flush(&mut self) -> io::Result<()> {
        // Source's printf call sites do not call fflush. In particular do not
        // flush unrelated streams or turn a C return value into an exception.
        Ok(())
    }
}

mod cstdio {
    use super::Kind;
    #[cfg(windows)]
    use crate::util::io::windows_crt_bindings as crt;
    use std::sync::Mutex;

    struct StreamState {
        failed: bool,
    }
    static COUT: Mutex<StreamState> = Mutex::new(StreamState { failed: false });
    static CERR: Mutex<StreamState> = Mutex::new(StreamState { failed: false });

    fn state(kind: Kind) -> &'static Mutex<StreamState> {
        match kind {
            Kind::Out => &COUT,
            Kind::Err => &CERR,
        }
    }
    pub(super) const SUPPORTED: bool = cfg!(any(
        windows,
        all(
            target_os = "linux",
            any(target_env = "gnu", target_env = "musl")
        ),
        target_vendor = "apple"
    ));
    fn pointer(kind: Kind) -> *mut libc::FILE {
        // Borrow the same FILE used by printf/std::cout. Never fdopen an
        // alias, close it, change mode or flush unrelated streams.
        #[cfg(windows)]
        {
            unsafe {
                crt::__acrt_iob_func(match kind {
                    Kind::Out => 1,
                    Kind::Err => 2,
                })
            }
        }
        #[cfg(any(
            all(target_os = "linux", any(target_env = "gnu", target_env = "musl")),
            target_vendor = "apple"
        ))]
        {
            unsafe {
                match kind {
                    Kind::Out => crate::util::io::posix_stdio_bindings::stdout(),
                    Kind::Err => crate::util::io::posix_stdio_bindings::stderr(),
                }
            }
        }
        #[cfg(not(any(
            windows,
            all(target_os = "linux", any(target_env = "gnu", target_env = "musl")),
            target_vendor = "apple"
        )))]
        {
            let _ = kind;
            std::ptr::null_mut()
        }
    }
    pub(super) fn flush(kind: Kind) {
        let mut state = state(kind).lock().unwrap_or_else(|p| p.into_inner());
        if state.failed {
            return;
        }
        let stream = pointer(kind);
        if stream.is_null() || unsafe { libc::fflush(stream) } != 0 {
            // ostream exceptions are off. Retain the fail state, do not turn
            // this into a new Source command exception or retry the output.
            state.failed = true;
        }
    }
    pub(super) fn write(kind: Kind, bytes: &[u8]) {
        if matches!(kind, Kind::Err) {
            super::flush_tied_stdout();
        }
        let mut state = state(kind).lock().unwrap_or_else(|p| p.into_inner());
        if state.failed {
            return;
        }
        let stream = pointer(kind);
        if stream.is_null()
            || unsafe { libc::fwrite(bytes.as_ptr().cast(), 1, bytes.len(), stream) } != bytes.len()
        {
            state.failed = true;
        } else if matches!(kind, Kind::Err) && unsafe { libc::fflush(stream) } != 0 {
            // std::cerr is unitbuf; only its own FILE is synchronized.
            state.failed = true;
        }
    }
    pub(super) fn perror(bytes: &[u8]) {
        // C perror is independent of cerr's failbit and has no cout tie.
        let stream = pointer(Kind::Err);
        if !stream.is_null() {
            unsafe { libc::fwrite(bytes.as_ptr().cast(), 1, bytes.len(), stream) };
        }
    }
    #[cfg(windows)]
    pub(super) fn append_log(bytes: &[u8]) {
        // MessageStream creates a fresh std::ofstream(out|app) for each
        // insertion. In the MSVC filebuf that is an unchecked text-mode FILE;
        // default exceptions remain off, including open/write/close failures.
        let stream = unsafe { libc::fopen(c"diamond.log".as_ptr(), c"a".as_ptr()) };
        if !stream.is_null() {
            unsafe {
                libc::fwrite(bytes.as_ptr().cast(), 1, bytes.len(), stream);
                libc::fclose(stream);
            }
        }
    }
}

impl Write for ConsoleWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if cstdio::SUPPORTED {
            // Drain earlier Rust-host output before entering Source's stream.
            // Source-owned console writes themselves stay in the actual FILE.
            if matches!(self.kind, Kind::Out) {
                let _ = io::stdout().flush();
            }
            cstdio::write(self.kind, bytes);
            Ok(bytes.len())
        } else {
            match self.kind {
                Kind::Out => io::stdout().write(bytes),
                Kind::Err => io::stderr().write(bytes),
            }
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        if cstdio::SUPPORTED {
            if matches!(self.kind, Kind::Err) {
                flush_tied_stdout();
            }
            cstdio::flush(self.kind);
            Ok(())
        } else {
            match self.kind {
                Kind::Out => io::stdout().flush(),
                Kind::Err => io::stderr().flush(),
            }
        }
    }
    fn write_fmt(&mut self, args: fmt::Arguments<'_>) -> io::Result<()> {
        if self.newline_endls {
            let text = args.to_string();
            for part in text.split_inclusive('\n') {
                self.write_all(part.as_bytes())?;
                if part.ends_with('\n') {
                    self.flush()?;
                }
            }
            return Ok(());
        }
        if cstdio::SUPPORTED {
            let text = args.to_string();
            self.write_all(text.as_bytes())?;
            if self.format_endl {
                self.flush()?;
            }
            Ok(())
        } else {
            match self.kind {
                Kind::Out => io::stdout().write_fmt(args),
                Kind::Err => io::stderr().write_fmt(args),
            }
        }
    }
}

pub(crate) fn flush_tied_stdout() {
    // Retain compatibility with still-connected Rust console writers, then
    // flush precisely cout's FILE. On MSVC FileSink stdout shares this FILE;
    // GNU FileSink fdopen(1) owns a different FILE and is not flushed here.
    let _ = io::stdout().flush();
    if cstdio::SUPPORTED {
        cstdio::flush(Kind::Out);
    }
}

pub(crate) fn write_perror(bytes: &[u8]) {
    if cstdio::SUPPORTED {
        cstdio::perror(bytes);
    } else {
        let _ = io::stderr().lock().write_all(bytes);
    }
}

#[cfg(windows)]
pub(crate) fn append_log(bytes: &[u8]) {
    cstdio::append_log(bytes);
}

pub(crate) fn stdout_fmt(args: fmt::Arguments<'_>) {
    let _ = stdout().write_fmt(args);
}
pub(crate) fn stdout_endls_fmt(args: fmt::Arguments<'_>) {
    let _ = stdout_endl().write_fmt(args);
}
pub(crate) fn stderr_fmt(args: fmt::Arguments<'_>) {
    let _ = stderr().write_fmt(args);
}
pub(crate) fn stderr_line(args: fmt::Arguments<'_>) {
    let _ = writeln!(stderr(), "{args}");
}
pub(crate) fn stdout_line(args: fmt::Arguments<'_>) {
    let mut out = stdout();
    let _ = writeln!(out, "{args}");
    let _ = out.flush();
}

macro_rules! source_print { ($($arg:tt)*) => { $crate::util::io::console::stdout_fmt(format_args!($($arg)*)) }; }
macro_rules! source_eprint { ($($arg:tt)*) => { $crate::util::io::console::stderr_fmt(format_args!($($arg)*)) }; }
macro_rules! source_eprintln {
    () => { $crate::util::io::console::stderr_line(format_args!("")) };
    ($($arg:tt)+) => { $crate::util::io::console::stderr_line(format_args!($($arg)*)) };
}
macro_rules! source_println {
    () => { $crate::util::io::console::stdout_line(format_args!("")) };
    ($($arg:tt)+) => { $crate::util::io::console::stdout_line(format_args!($($arg)*)) };
}
pub(crate) use {source_eprint, source_eprintln, source_print, source_println};
