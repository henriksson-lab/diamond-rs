use std::fs::File as StdFile;
use std::io::SeekFrom;
#[cfg(test)]
use std::io::Write;

/// Source join_worker handles its own std::exception and exits the standalone
/// process immediately. Embedded calls retain their typed error and host.
/// Runtime initialization alone never grants this process ownership policy.
pub(crate) fn exit_source_worker_error(error: &std::io::Error) {
    #[cfg(any(unix, windows))]
    let process_owned = file_sink::owns_process_stdio();
    #[cfg(not(any(unix, windows)))]
    let process_owned = false;
    if !process_owned {
        return;
    }
    let what = if let Some(raw) = RawDiagnosticError::find(error) {
        Ok(raw.to_vec())
    } else if let Some(open) = file_open_error::find(error) {
        open.try_message_bytes()
    } else if let Some(temp) = temp_file::find_temp_open(error) {
        temp.try_message_bytes()
    } else if let Some(stdio) = cstdio_file::StdioError::find(error) {
        Ok(stdio.message_bytes().to_vec())
    } else if BadAllocError::find(error).is_some() {
        Ok(if cfg!(target_env = "msvc") {
            b"bad allocation".to_vec()
        } else {
            b"std::bad_alloc".to_vec()
        })
    } else {
        Ok(error.to_string().into_bytes())
    };
    let what = what.unwrap_or_else(|conversion: IoError| conversion.to_string().into_bytes());
    // A constructor's saved perror precedes the worker's tied cerr catch.
    let _ = emit_source_perror(error);
    crate::util::log_stream::flush_tied_stdout();
    let mut message = b"Error: ".to_vec();
    message.extend_from_slice(what.split(|byte| *byte == 0).next().unwrap_or_default());
    message.push(b'\n');
    {
        use std::io::Write;
        let _ = console::stderr().write_all(&message);
        let mut log = crate::util::log_stream::log_stream()
            .lock()
            .unwrap_or_else(|poison| poison.into_inner());
        let _ = log.write_bytes(&message).and_then(|stream| stream.flush());
        let _ = std::io::stdout().flush();
        let _ = std::io::stderr().flush();
    }
    // Source exit() flushes CRT owners without unwinding worker/file objects.
    unsafe { libc::exit(1) }
}

pub mod deserializer;
pub use deserializer::Deserializer;
pub mod serializer;
pub use serializer::Serializer;
pub mod file;
pub use file::{File, FilePrimitive, Temporary};
pub mod file_open_error;
pub use file_open_error::{path_diagnostic_bytes, FileOpenError};
pub mod file_source;
pub use file_source::FileSource;
pub mod file_sink;
pub use file_sink::FileSink;
pub(crate) mod console;
pub mod cstdio_file;
#[cfg(any(
    all(target_os = "linux", any(target_env = "gnu", target_env = "musl")),
    target_vendor = "apple"
))]
mod posix_stdio_bindings;
pub mod temp_file;
#[cfg(windows)]
mod windows_crt_bindings;
#[cfg(windows)]
mod windows_cstdio;
pub use temp_file::{TempFile, TempFileData, TempFileHandler, TEMP_FILE_HANDLER};
pub mod output_file;
pub use output_file::{decompress, decompress_file, Compressor, OutputFile};
pub mod zstd_stream;
#[cfg(test)]
pub(crate) use zstd_stream::zstd_compress_for_test;
pub use zstd_stream::{zstd_decompress, zstd_decompress_file, ZstdSink, ZstdSource};
pub mod compressed_stream;
pub use compressed_stream::{zlib_decompress, zlib_decompress_file, ZlibSink, ZlibSource};
pub mod output_stream_buffer;
pub use output_stream_buffer::OutputStreamBuffer;
pub mod input_stream_buffer;
pub use input_stream_buffer::InputStreamBuffer;
pub(crate) mod source_text;
pub mod text_input_file;
pub use text_input_file::{ByteTextInputFile, TextInputFile};
pub mod input_file;
pub use input_file::{detect_compressor, InputFile, InputFileStream};

pub(crate) const DEFAULT_FILE_BUFFER_SIZE: usize = 1 << 20;
pub const MEGABYTES: usize = 1 << 20;
pub const GIGABYTES: usize = 1 << 30;
pub const KILOBYTES: usize = 1 << 10;

/// Duplicate an inherited Windows standard handle into an owned `File`.
///
/// Opening `CONIN$`/`CONOUT$` addresses the console device rather than the
/// process's inherited stream and therefore breaks pipes, redirection, and
/// headless CI jobs. Duplicating the inherited handle matches Unix's
/// `/dev/stdin` and `/dev/stdout` behavior while giving `FileSource` and
/// `FileSink` an independently owned RAII handle.
#[cfg(windows)]
pub(crate) fn duplicate_standard_handle(which: u32) -> std::io::Result<StdFile> {
    use std::os::windows::io::FromRawHandle;
    use windows_sys::Win32::Foundation::{
        DuplicateHandle, DUPLICATE_SAME_ACCESS, INVALID_HANDLE_VALUE,
    };
    use windows_sys::Win32::System::Console::GetStdHandle;
    use windows_sys::Win32::System::Threading::GetCurrentProcess;

    // SAFETY: `GetStdHandle` returns a borrowed process handle. DuplicateHandle
    // creates a distinct owned handle, which is transferred exactly once to
    // `StdFile`; dropping the file cannot close the original standard stream.
    unsafe {
        let source = GetStdHandle(which);
        if source.is_null() || source == INVALID_HANDLE_VALUE {
            return Err(std::io::Error::last_os_error());
        }
        let process = GetCurrentProcess();
        let mut duplicate = std::ptr::null_mut();
        if DuplicateHandle(
            process,
            source,
            process,
            &mut duplicate,
            0,
            0,
            DUPLICATE_SAME_ACCESS,
        ) == 0
        {
            return Err(std::io::Error::last_os_error());
        }
        Ok(StdFile::from_raw_handle(duplicate))
    }
}

/// A command exception whose terminal diagnostic was emitted while its
/// invocation-owned Source log was still alive. Retain the original cause for
/// callers while preventing the outer CLI facade from reporting it twice.
#[derive(Debug)]
pub(crate) struct ReportedIoError(std::io::Error);
impl std::fmt::Display for ReportedIoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}
impl std::error::Error for ReportedIoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.0)
    }
}
impl ReportedIoError {
    pub(crate) fn wrap(error: std::io::Error) -> std::io::Error {
        std::io::Error::new(error.kind(), Self(error))
    }
    pub(crate) fn contains(error: &(dyn std::error::Error + 'static)) -> bool {
        let mut current = Some(error);
        while let Some(error) = current {
            if error.is::<Self>() {
                return true;
            }
            if let Some(io) = error.downcast_ref::<std::io::Error>() {
                if let Some(inner) = io.get_ref() {
                    current = Some(inner);
                    continue;
                }
            }
            current = error.source();
        }
        false
    }
}

/// Only the perror portion of an exception has already been emitted. Its
/// ordinary Source catch (and log) must still run unless it is FileOpen.
#[derive(Debug)]
pub(crate) struct PerrorEmittedError(std::io::Error);
impl std::fmt::Display for PerrorEmittedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}
impl std::error::Error for PerrorEmittedError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.0)
    }
}
impl PerrorEmittedError {
    pub(crate) fn wrap(error: std::io::Error) -> std::io::Error {
        if Self::contains(&error) {
            error
        } else {
            std::io::Error::new(error.kind(), Self(error))
        }
    }
    pub(crate) fn contains(error: &(dyn std::error::Error + 'static)) -> bool {
        let mut current = Some(error);
        while let Some(error) = current {
            if error.is::<Self>() {
                return true;
            }
            if let Some(io) = error.downcast_ref::<std::io::Error>() {
                if let Some(inner) = io.get_ref() {
                    current = Some(inner);
                    continue;
                }
            }
            current = error.source();
        }
        false
    }
}

/// Captured Source perror bytes, before any generic exception catch. There
/// is no cout tie here; a later cerr/what diagnostic owns that separate flush.
pub(crate) fn source_perror_bytes(
    error: &(dyn std::error::Error + 'static),
) -> IoResult<Option<Vec<u8>>> {
    if PerrorEmittedError::contains(error) || ReportedIoError::contains(error) {
        return Ok(None);
    }
    let mut bytes = if let Some(temp) = temp_file::find_temp_open(error) {
        temp.reason.clone()
    } else if let Some(open) = file_open_error::find(error) {
        let mut bytes = b"\n".to_vec();
        bytes.extend_from_slice(&open.try_message_bytes()?);
        bytes.extend_from_slice(b": ");
        bytes.extend_from_slice(&open.reason);
        bytes
    } else {
        let Some(stdio) = cstdio_file::StdioError::find(error) else {
            return Ok(None);
        };
        let Some(perror) = &stdio.perror else {
            return Ok(None);
        };
        perror.clone()
    };
    if let Some(nul) = bytes.iter().position(|&byte| byte == 0) {
        bytes.truncate(nul);
    }
    bytes.push(b'\n');
    Ok(Some(bytes))
}

pub(crate) fn emit_source_perror(error: &(dyn std::error::Error + 'static)) -> IoResult<bool> {
    if let Some(bytes) = source_perror_bytes(error)? {
        console::write_perror(&bytes);
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Report only a diagnostic-codec failure at a Source noexcept boundary.
/// The argument must be the conversion error, never the original user-path
/// error; raw typed messages are retained and fallback text is ASCII only.
pub(crate) fn abort_source_diagnostic_error(error: &(dyn std::error::Error + 'static)) -> ! {
    let text;
    let message = if let Some(bytes) = RawDiagnosticError::find(error) {
        bytes
    } else {
        text = error.to_string();
        if text.is_ascii() {
            text.as_bytes()
        } else {
            b"Source diagnostic encoding failed."
        }
    };
    let mut bytes = b"Error: ".to_vec();
    bytes.extend_from_slice(message);
    bytes.push(b'\n');
    let _ = std::io::Write::write_all(&mut console::stderr(), &bytes);
    let _ = std::io::Write::flush(&mut console::stderr());
    abort_source_process()
}

/// Source's terminate handler calls std::abort in the selected C runtime.
/// In particular, Windows must use that CRT's abort behavior and SIGABRT path.
pub(crate) fn abort_source_process() -> ! {
    #[cfg(any(unix, windows))]
    unsafe {
        libc::abort()
    }
    #[cfg(not(any(unix, windows)))]
    std::process::abort()
}

/// Concrete native counterpart of a failed Source allocation. Keeping its
/// class through io::Error distinguishes it from an OS I/O failure with ENOMEM.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BadAllocError(pub String);
impl Default for BadAllocError {
    fn default() -> Self {
        Self("std::bad_alloc".into())
    }
}
impl std::fmt::Display for BadAllocError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for BadAllocError {}
impl BadAllocError {
    pub fn into_io(self) -> std::io::Error {
        std::io::Error::new(std::io::ErrorKind::OutOfMemory, self)
    }
    pub fn find<'a>(error: &'a (dyn std::error::Error + 'static)) -> Option<&'a Self> {
        let mut current = Some(error);
        while let Some(error) = current {
            if let Some(allocation) = error.downcast_ref::<Self>() {
                return Some(allocation);
            }
            if let Some(io) = error.downcast_ref::<std::io::Error>() {
                if let Some(inner) = io.get_ref() {
                    current = Some(inner);
                    continue;
                }
            }
            current = error.source();
        }
        None
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IoError {
    UnsupportedOperation,
    FileOpen(String),
    RawDetails(RawDiagnosticError),
    FileOpenDetails(FileOpenError),
    TemporaryFileOpenDetails(temp_file::TempFileOpenError),
    StdioDetails(cstdio_file::StdioError),
    FileRead(String),
    FileWrite(String),
    EndOfStream,
    StreamRead { line_count: usize, msg: String },
    Other(String),
}

impl std::fmt::Display for IoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IoError::RawDetails(error) => std::fmt::Display::fmt(error, f),
            IoError::UnsupportedOperation => f.write_str("Unsupported I/O operation."),
            IoError::FileOpen(file_name) => write!(f, "Error opening file {file_name}"),
            IoError::FileOpenDetails(error) => std::fmt::Display::fmt(error, f),
            IoError::TemporaryFileOpenDetails(error) => std::fmt::Display::fmt(error, f),
            IoError::StdioDetails(error) => std::fmt::Display::fmt(error, f),
            IoError::FileRead(file_name) => write!(f, "Error reading file {file_name}"),
            IoError::FileWrite(file_name) => write!(f, "Error writing file {file_name}"),
            IoError::EndOfStream => f.write_str("Unexpected end of input."),
            IoError::StreamRead { line_count, msg } => {
                write!(f, "Error reading input stream at line {line_count}: {msg}")
            }
            IoError::Other(msg) => f.write_str(msg),
        }
    }
}

impl IoError {
    pub fn path_message(prefix: &[u8], path: &std::path::Path) -> Self {
        let path = match path_diagnostic_bytes(path) {
            Ok(bytes) => bytes,
            Err(error) => return error,
        };
        let mut bytes = prefix.to_vec();
        bytes.extend_from_slice(&path);
        Self::RawDetails(RawDiagnosticError(bytes))
    }
    pub fn file_open(path: impl AsRef<std::path::Path>, error: std::io::Error) -> Self {
        Self::FileOpenDetails(FileOpenError::new(path, error))
    }
}
impl std::error::Error for IoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::RawDetails(error) => Some(error),
            Self::FileOpenDetails(error) => Some(error),
            Self::TemporaryFileOpenDetails(error) => Some(error),
            Self::StdioDetails(error) => Some(error),
            _ => None,
        }
    }
}

/// The C pathname boundary used by Source fopen/stat/remove calls.
/// Retained owners keep the complete path; only the individual OS operation
/// consumes this prefix, just as std::string::c_str() does in Source.
pub fn source_c_path(path: &std::path::Path) -> std::borrow::Cow<'_, std::path::Path> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        let bytes = path.as_os_str().as_bytes();
        return match bytes.iter().position(|&byte| byte == 0) {
            Some(end) => std::borrow::Cow::Borrowed(std::path::Path::new(
                std::ffi::OsStr::from_bytes(&bytes[..end]),
            )),
            None => std::borrow::Cow::Borrowed(path),
        };
    }
    #[cfg(windows)]
    {
        use std::os::windows::ffi::{OsStrExt, OsStringExt};
        let units = path.as_os_str().encode_wide().collect::<Vec<_>>();
        return match units.iter().position(|&unit| unit == 0) {
            Some(end) => std::borrow::Cow::Owned(std::path::PathBuf::from(
                std::ffi::OsString::from_wide(&units[..end]),
            )),
            None => std::borrow::Cow::Borrowed(path),
        };
    }
    #[cfg(not(any(unix, windows)))]
    {
        match path
            .to_str()
            .and_then(|text| text.find('\0').map(|end| &text[..end]))
        {
            Some(prefix) => std::borrow::Cow::Borrowed(std::path::Path::new(prefix)),
            None => std::borrow::Cow::Borrowed(path),
        }
    }
}

/// Convert an OS argument to the pinned Source narrow-string domain.
/// Unix borrows the original bytes; Windows uses the CRT ANSI code page,
/// independently of the retained wide path used by native filename opens.
pub fn source_os_bytes(value: &std::ffi::OsStr) -> IoResult<std::borrow::Cow<'_, [u8]>> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        Ok(std::borrow::Cow::Borrowed(value.as_bytes()))
    }
    #[cfg(windows)]
    {
        windows_cstdio::source_os_bytes(value).map(std::borrow::Cow::Owned)
    }
    #[cfg(not(any(unix, windows)))]
    {
        value
            .to_str()
            .map(|text| std::borrow::Cow::Borrowed(text.as_bytes()))
            .ok_or(IoError::UnsupportedOperation)
    }
}
pub fn source_bytes_os(bytes: &[u8]) -> IoResult<std::ffi::OsString> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStringExt;
        Ok(std::ffi::OsString::from_vec(bytes.to_vec()))
    }
    #[cfg(windows)]
    {
        windows_cstdio::source_bytes_os(bytes)
    }
    #[cfg(not(any(unix, windows)))]
    {
        String::from_utf8(bytes.to_vec())
            .map(Into::into)
            .map_err(|error| IoError::Other(error.to_string()))
    }
}

pub type IoResult<T> = Result<T, IoError>;

pub trait Consumer {
    fn consume(&mut self, ptr: &[u8]) -> IoResult<()>;
    fn finalize(&mut self) -> IoResult<()> {
        Ok(())
    }
}

pub trait StreamEntity {
    fn rewind(&mut self) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn seek(&mut self, _p: i64, _origin: SeekFrom) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn tell(&mut self) -> IoResult<i64> {
        Err(IoError::UnsupportedOperation)
    }

    fn read(&mut self, _ptr: &mut [u8]) -> IoResult<usize> {
        Err(IoError::UnsupportedOperation)
    }

    fn fetch(&mut self) -> IoResult<bool> {
        Err(IoError::UnsupportedOperation)
    }

    fn close(&mut self) -> IoResult<()> {
        Ok(())
    }

    /// Canonical OS path. Forwarding streams delegate this directly.
    fn file_path(&self) -> &std::path::Path {
        std::path::Path::new(self.file_name())
    }
    /// Legacy UTF-8 accessor. Actual Source file operations use file_path.
    fn file_name(&self) -> &str {
        ""
    }
    #[cfg(unix)]
    fn file_name_bytes(&self) -> &[u8] {
        use std::os::unix::ffi::OsStrExt;
        self.file_path().as_os_str().as_bytes()
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        Err(IoError::UnsupportedOperation)
    }

    /// Source FILE alias service. Ordinary streams retain their existing
    /// independent descriptor constructors unless they expose a shared owner.
    #[cfg(any(unix, windows))]
    fn shared_stdio(&mut self) -> IoResult<Option<cstdio_file::SharedCstdio>> {
        Ok(None)
    }

    fn clone_file(&mut self) -> IoResult<StdFile> {
        let name = self.file_path().to_path_buf();
        self.file()?
            .try_clone()
            .map_err(|error| IoError::file_open(name, error))
    }

    fn data(&self) -> &[u8] {
        &[]
    }

    fn write(&mut self, _ptr: &[u8]) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn flush(&mut self) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn eof(&mut self) -> IoResult<bool> {
        Err(IoError::UnsupportedOperation)
    }

    fn file_size(&mut self) -> IoResult<i64> {
        Err(IoError::UnsupportedOperation)
    }

    fn seekable(&self) -> bool {
        false
    }
}

pub mod compressed_buffer;
pub use compressed_buffer::CompressedBuffer;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VecStream {
    data: Vec<u8>,
    pos: usize,
    name: String,
}

impl VecStream {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_vec(data: Vec<u8>) -> Self {
        Self {
            data,
            pos: 0,
            name: String::new(),
        }
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }
}

impl StreamEntity for VecStream {
    fn rewind(&mut self) -> IoResult<()> {
        self.pos = 0;
        Ok(())
    }

    fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        let pos = match origin {
            SeekFrom::Start(_) => p,
            SeekFrom::End(_) => self.data.len() as i64 + p,
            SeekFrom::Current(_) => self.pos as i64 + p,
        };
        if pos < 0 {
            return Err(IoError::UnsupportedOperation);
        }
        self.pos = pos as usize;
        if self.pos > self.data.len() {
            self.data.resize(self.pos, 0);
        }
        Ok(())
    }

    fn tell(&mut self) -> IoResult<i64> {
        Ok(self.pos as i64)
    }

    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let n = ptr.len().min(self.data.len().saturating_sub(self.pos));
        ptr[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }

    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        if self.pos + ptr.len() > self.data.len() {
            self.data.resize(self.pos + ptr.len(), 0);
        }
        self.data[self.pos..self.pos + ptr.len()].copy_from_slice(ptr);
        self.pos += ptr.len();
        Ok(())
    }

    fn flush(&mut self) -> IoResult<()> {
        Ok(())
    }

    fn eof(&mut self) -> IoResult<bool> {
        Ok(self.pos >= self.data.len())
    }

    fn file_size(&mut self) -> IoResult<i64> {
        Ok(self.data.len() as i64)
    }

    fn seekable(&self) -> bool {
        true
    }

    fn file_name(&self) -> &str {
        &self.name
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_worker_exit_probe() {
        let Some(marker) = std::env::var_os("DIAMOND_RS_WORKER_EXIT_PROBE_MARKER") else {
            return;
        };
        struct UnwindMarker(std::path::PathBuf);
        impl Drop for UnwindMarker {
            fn drop(&mut self) {
                let _ = std::fs::write(&self.0, b"unwound");
            }
        }
        let _marker = UnwindMarker(marker.into());
        *crate::util::log_stream::log_stream().lock().unwrap() =
            crate::util::log_stream::MessageStream::new(false, true);
        // This isolated subprocess owns its descriptors just like main().
        unsafe {
            file_sink::enable_process_stdout_ownership();
        }
        exit_source_worker_error(
            &RawDiagnosticError(b"worker \xff\0not printed".to_vec()).into_io(),
        );
        panic!("Source process worker catch returned");
    }

    #[test]
    fn source_worker_failure_preserves_embedded_host_and_exits_owned_child_without_unwind() {
        let error = RawDiagnosticError(b"embedded failure".to_vec()).into_io();
        exit_source_worker_error(&error); // No process ownership in the host.
        assert_eq!(
            RawDiagnosticError::find(&error),
            Some(b"embedded failure".as_slice())
        );
        let root = std::env::temp_dir().join(format!(
            "diamond-worker-error-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&root).unwrap();
        let marker = root.join("unwind-marker");
        let child = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "util::io::tests::source_worker_exit_probe",
                "--nocapture",
            ])
            .env("DIAMOND_RS_WORKER_EXIT_PROBE_MARKER", &marker)
            .current_dir(&root)
            .output()
            .unwrap();
        assert_eq!(child.status.code(), Some(1));
        let expected = source_text::expected_text(b"Error: worker \xff\n");
        assert_eq!(child.stderr, expected);
        assert_eq!(std::fs::read(root.join("diamond.log")).unwrap(), expected);
        assert!(
            !marker.exists(),
            "Source exit must bypass worker object destruction"
        );
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn typed_io_capture_retains_primary_and_returned_stdio_error() {
        let details = cstdio_file::StdioError {
            operation: cstdio_file::StdioOperation::Read,
            path: "primary".into(),
            raw_os_error: Some(5),
            perror: Some(b"Input/output error".to_vec()),
            message: b"Error reading file primary".to_vec(),
        };
        let owner = IoErrorOwner::default().with_source_destructors();
        let returned = owner
            .capture_io::<()>(Err(std::io::Error::other(IoError::StdioDetails(
                details.clone(),
            ))))
            .unwrap_err();
        assert_eq!(cstdio_file::StdioError::find(&returned), Some(&details));
        let restored = owner
            .finish::<()>(Err("transport String".into()))
            .unwrap_err();
        assert_eq!(cstdio_file::StdioError::find(&restored), Some(&details));
        // Source-mode finish leaves an owned copy for later noexcept cleanup.
        let restored_again = owner
            .finish::<()>(Err("transport String".into()))
            .unwrap_err();
        assert_eq!(
            cstdio_file::StdioError::find(&restored_again),
            Some(&details)
        );
        assert_eq!(returned.to_string(), "Error reading file primary");
    }

    #[test]
    fn test_serializer_deserializer_roundtrip() {
        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        ser.write_i16(-2).unwrap();
        ser.write_u16(0x0102).unwrap();
        ser.write_i32(-7).unwrap();
        ser.write_u64(0x0102_0304_0506_0708).unwrap();
        ser.write_f64(1.5).unwrap();
        ser.write_value(-9i32).unwrap();
        ser.write_values(&[11u16, 12u16, 13u16]).unwrap();
        ser.write_string("abc").unwrap();
        let mut stream = ser.into_inner().unwrap();
        stream.rewind().unwrap();

        let mut de = Deserializer::new(stream);
        assert_eq!(de.read_i16().unwrap(), -2);
        assert_eq!(de.read_u16().unwrap(), 0x0102);
        assert_eq!(de.read_i32().unwrap(), -7);
        assert_eq!(de.read_u64().unwrap(), 0x0102_0304_0506_0708);
        assert_eq!(de.read_f64().unwrap(), 1.5);
        assert_eq!(de.read_value::<i32>().unwrap(), -9);
        let mut values = [0u16; 3];
        de.read_values(&mut values).unwrap();
        assert_eq!(values, [11, 12, 13]);
        assert_eq!(de.read_string().unwrap(), "abc");
    }

    #[test]
    fn test_serializer_buffered_write_seek_and_reset() {
        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        let data = vec![b'a'; DEFAULT_FILE_BUFFER_SIZE + 17];
        ser.write_raw(&data).unwrap();
        assert_eq!(ser.tell().unwrap(), DEFAULT_FILE_BUFFER_SIZE + 17);
        ser.seek(4, SeekFrom::Start(0)).unwrap();
        ser.write_raw(b"BC").unwrap();
        let stream = ser.into_inner().unwrap();
        assert_eq!(&stream.data()[..4], b"aaaa");
        assert_eq!(&stream.data()[4..6], b"BC");
        assert_eq!(stream.data().len(), DEFAULT_FILE_BUFFER_SIZE + 17);

        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        ser.reset_buffer();
        ser.write_raw(b"xy").unwrap();
        assert_eq!(ser.into_inner().unwrap().data(), b"xy");
    }

    #[test]
    fn test_deserializer_read_to_seek_and_peek() {
        let stream = VecStream::from_vec(b"abc,def\n>rec\nx".to_vec());
        let mut de = Deserializer::new(stream);
        assert!(de.data().is_empty());
        assert_eq!(de.peek(3).unwrap(), "abc");
        let mut out = Vec::new();
        assert!(de.read_to(&mut out, b',').unwrap());
        assert_eq!(out, b"abc");
        de.seek_forward(3).unwrap();
        let mut out = Vec::new();
        assert_eq!(
            de.read_to_record_start(&mut out, b'\n', b'>').unwrap(),
            (true, 1)
        );
        assert!(out.is_empty());
        assert_eq!(de.read_string().unwrap_err(), IoError::EndOfStream);

        let mut input = InputStreamBuffer::new(VecStream::from_vec(b"peek-pop".to_vec()), 0);
        input.fetch().unwrap();
        let mut de = Deserializer::new(input);
        assert_eq!(de.peek(4).unwrap(), "peek");
        let mut popped = [0u8; 4];
        de.pop(&mut popped).unwrap();
        assert_eq!(&popped, b"peek");
    }

    #[test]
    fn test_file_temporary_write_read_size() {
        let mut file = File::new_temporary(Temporary).unwrap();
        file.write(b"abcd").unwrap();
        file.write_value(7u16).unwrap();
        assert_eq!(file.size().unwrap(), 6);
        file.seek(0, SeekFrom::Start(0)).unwrap();
        let mut buf = [0u8; 4];
        file.read_exact(&mut buf).unwrap();
        assert_eq!(&buf, b"abcd");
        assert!(!file.eof().unwrap());
        assert_eq!(file.read(2).unwrap(), &7u16.to_ne_bytes());
        assert!(file.eof().unwrap());
        file.seek(4, SeekFrom::Start(0)).unwrap();
        assert_eq!(file.read_value::<u16>().unwrap(), 7);
        assert!(!file.file_name().is_empty());
        #[cfg(unix)]
        assert!(!std::path::Path::new(file.file_name()).exists());
        assert!(file.file().is_ok());
        file.close().unwrap();

        let path = std::env::temp_dir().join(format!(
            "diamond-rs-file-new-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut named = File::new(&name, "w+b").unwrap();
            named.write(b"x").unwrap();
            assert_eq!(named.file_name(), name);
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_stream_entity_file_delegation() {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "diamond-rs-stream-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut serializer = Serializer::new(OutputStreamBuffer::new(
                FileSink::new(&name, "w+b", false, 0).unwrap(),
            ));
            serializer.write_raw(b"abcdef").unwrap();
            assert!(serializer.file().is_ok());
            serializer.close().unwrap();
        }
        {
            let mut de =
                Deserializer::new(InputStreamBuffer::new(FileSource::new(&name).unwrap(), 0));
            assert!(de.file().is_ok());
            let mut b = [0u8; 2];
            de.read_exact(&mut b).unwrap();
            assert_eq!(&b, b"ab");
            de.close().unwrap();
        }
        let _ = std::fs::remove_file(&name);
    }

    #[test]
    fn test_errors_display() {
        assert_eq!(
            IoError::UnsupportedOperation.to_string(),
            "Unsupported I/O operation."
        );
        assert_eq!(IoError::EndOfStream.to_string(), "Unexpected end of input.");
        assert!(IoError::StreamRead {
            line_count: 3,
            msg: "bad".to_string()
        }
        .to_string()
        .contains("line 3"));
    }

    #[test]
    fn test_detect_compressor_and_decompress() {
        assert_eq!(detect_compressor(&[0x1f, 0x8b, 0, 0]), Compressor::Zlib);
        assert_eq!(detect_compressor(&[0x78, 0x9c, 0, 0]), Compressor::Zlib);
        assert_eq!(
            detect_compressor(&[0x28, 0xb5, 0x2f, 0xfd]),
            Compressor::Zstd
        );
        assert_eq!(detect_compressor(b"plain"), Compressor::None);

        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"abcabc").unwrap();
        let compressed = zlib_encoder.finish().unwrap();
        let mut out = [0u8; 16];
        let n = decompress(&compressed, &mut out, Compressor::Zlib).unwrap();
        assert_eq!(&out[..n], b"abcabc");

        let mut gzip_encoder =
            flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip_encoder.write_all(b"xyz").unwrap();
        let compressed = gzip_encoder.finish().unwrap();
        let n = decompress(&compressed, &mut out, Compressor::Zlib).unwrap();
        assert_eq!(&out[..n], b"xyz");
        let n = zlib_decompress(&compressed, &mut out).unwrap();
        assert_eq!(&out[..n], b"xyz");

        let compressed = zstd_compress_for_test(b"zstd-data");
        let n = decompress(&compressed, &mut out, Compressor::Zstd).unwrap();
        assert_eq!(&out[..n], b"zstd-data");
        let n = zstd_decompress(&compressed, &mut out).unwrap();
        assert_eq!(&out[..n], b"zstd-data");
        assert!(zstd_decompress(&compressed, &mut [0u8; 2]).is_err());

        assert!(decompress(&compressed, &mut [0u8; 2], Compressor::Zlib).is_err());

        let mut first = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        first.write_all(b"one").unwrap();
        let mut concatenated = first.finish().unwrap();
        let mut second =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        second.write_all(b"two").unwrap();
        concatenated.extend_from_slice(&second.finish().unwrap());
        let mut out = [0u8; 6];
        let n = zlib_decompress(&concatenated, &mut out).unwrap();
        assert_eq!(&out[..n], b"onetwo");
    }

    #[test]
    fn test_compressed_buffer_roundtrip_and_clear() {
        let compressor = if cfg!(feature = "zstd") {
            Compressor::Zstd
        } else {
            Compressor::Zlib
        };
        let mut compressed = CompressedBuffer::new();
        compressed.write(b"abc").unwrap();
        compressed.write_value(7u16).unwrap();
        compressed.write(b"xyzxyzxyz").unwrap();
        compressed.finish().unwrap();
        assert!(compressed.size() > 0);
        assert_eq!(detect_compressor(compressed.data()), compressor);

        let mut out = [0u8; 32];
        let n = decompress(compressed.data(), &mut out, compressor).unwrap();
        let mut expected = Vec::new();
        expected.extend_from_slice(b"abc");
        expected.extend_from_slice(&7u16.to_ne_bytes());
        expected.extend_from_slice(b"xyzxyzxyz");
        assert_eq!(&out[..n], expected.as_slice());

        compressed.clear();
        assert_eq!(compressed.size(), 0);
        compressed.write(b"second").unwrap();
        compressed.finish().unwrap();
        let n = decompress(compressed.data(), &mut out, compressor).unwrap();
        assert_eq!(&out[..n], b"second");
    }

    #[test]
    fn test_file_source_sink_and_stream_buffers() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let mut out = OutputStreamBuffer::new(sink);
            out.write(b"abc").unwrap();
            out.write(b"def").unwrap();
            assert_eq!(out.tell().unwrap(), 6);
            out.close().unwrap();
        }
        {
            let source = FileSource::new(&name).unwrap();
            assert!(source.seekable());
            let mut input = InputStreamBuffer::new(source, 0);
            let mut buf = [0u8; 6];
            assert_eq!(input.read(&mut buf).unwrap(), 6);
            assert_eq!(&buf, b"abcdef");
            assert_eq!(input.file_size().unwrap(), 6);
        }
        let _ = std::fs::remove_file(path);
    }

    #[cfg(unix)]
    #[test]
    fn test_file_source_sink_standard_stream_names() {
        let stdin_empty = FileSource::new("").unwrap();
        assert_eq!(stdin_empty.file_name(), "");
        assert!(!stdin_empty.seekable());

        let stdin_dash = FileSource::new("-").unwrap();
        assert_eq!(stdin_dash.file_name(), "-");
        assert!(!stdin_dash.seekable());

        let mut stdout_sink = FileSink::new("", "wb", false, 0).unwrap();
        assert_eq!(stdout_sink.file_name(), "");
        stdout_sink.close().unwrap();
    }

    #[test]
    fn test_file_sink_modes_and_from_file() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-fd-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        let file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&path)
            .unwrap();
        assert_eq!(
            FileSink::new(&name, "bad", false, 0)
                .unwrap_err()
                .to_string(),
            "Invalid fopen mode."
        );
        assert_eq!(
            FileSink::from_file(&name, file.try_clone().unwrap(), "bad", false, 0)
                .unwrap_err()
                .to_string(),
            "Invalid fopen mode."
        );
        let mut sink = FileSink::from_file(&name, file, "w+b", false, 0).unwrap();
        sink.write(b"fd").unwrap();
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"fd");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_output_stream_buffer_write_buffer_and_flush_count() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-output-buffer-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let mut out = OutputStreamBuffer::new(sink);
            let (buf, size) = out.write_buffer_range();
            assert_eq!(size, DEFAULT_FILE_BUFFER_SIZE);
            buf[..4].copy_from_slice(b"abcd");
            out.flush_count(4).unwrap();
            out.write(b"ef").unwrap();
            out.close().unwrap();
        }
        assert_eq!(std::fs::read(&path).unwrap(), b"abcdef");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_compressed_stream_entities_roundtrip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-compressed-stream-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let out = OutputStreamBuffer::new(sink);
            let mut zlib = ZlibSink::new(out);
            zlib.write(b"stream-zlib").unwrap();
            zlib.close().unwrap();
        }
        assert_eq!(
            detect_compressor(&std::fs::read(&path).unwrap()),
            Compressor::Zlib
        );
        {
            let source = FileSource::new(&name).unwrap();
            let input = InputStreamBuffer::new(source, 0);
            let mut zlib = ZlibSource::new(input).unwrap();
            let mut buf = [0u8; 11];
            assert_eq!(zlib.read(&mut buf).unwrap(), 11);
            assert_eq!(&buf, b"stream-zlib");
            let mut tail = [0u8; 1];
            assert_eq!(zlib.read(&mut tail).unwrap(), 0);
            assert!(zlib.eof().unwrap());
            zlib.rewind().unwrap();
            let mut head = [0u8; 6];
            assert_eq!(zlib.read(&mut head).unwrap(), 6);
            assert_eq!(&head, b"stream");
        }

        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let out = OutputStreamBuffer::new(sink);
            let mut zstd = ZstdSink::new(out);
            zstd.write(b"stream-zstd").unwrap();
            zstd.close().unwrap();
        }
        assert_eq!(
            detect_compressor(&std::fs::read(&path).unwrap()),
            Compressor::Zstd
        );
        {
            let source = FileSource::new(&name).unwrap();
            let input = InputStreamBuffer::new(source, 0);
            let mut zstd = ZstdSource::new(input).unwrap();
            let mut buf = [0u8; 11];
            assert_eq!(zstd.read(&mut buf).unwrap(), 11);
            assert_eq!(&buf, b"stream-zstd");
            let mut tail = [0u8; 1];
            assert_eq!(zstd.read(&mut tail).unwrap(), 0);
            assert!(zstd.eof().unwrap());
        }

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_temp_file_data_and_temp_file() {
        let data = TempFileData::init(false).unwrap();
        {
            let mut out = OutputFile::from_temp_file_data(&data, Compressor::None, "w+b").unwrap();
            out.write_raw(b"temp").unwrap();
            out.close().unwrap();
        }
        assert_eq!(std::fs::read(&data.name).unwrap(), b"temp");
        std::fs::remove_file(&data.name).unwrap();

        let dir = TempFile::get_temp_dir().unwrap();
        // Source's empty configured directory names the current directory.
        assert_eq!(dir, "");

        let mut tmp = TempFile::new(true).unwrap();
        #[cfg(not(unix))]
        let tmp_name = tmp.file_name().to_string();
        assert_eq!(tmp.unlinked(), cfg!(unix));
        tmp.write_raw(b"x").unwrap();
        tmp.close().unwrap();
        #[cfg(not(unix))]
        std::fs::remove_file(tmp_name).unwrap();
    }

    #[test]
    fn test_output_file_input_file_and_hash() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-file-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
            out.write_i16(-123).unwrap();
            out.write_u16(0xCAFE).unwrap();
            out.write_i32(-7).unwrap();
            out.write_u64(0x0102_0304_0506_0708).unwrap();
            out.write_value(99u32).unwrap();
            out.write_values(&[3i16, 4i16]).unwrap();
            out.write_string("abc").unwrap();
            assert!(out.file().is_ok());
            assert!(out.tell().unwrap() > 0);
            out.close().unwrap();
        }
        {
            let mut input = InputFile::new(&name, 0).unwrap();
            assert_eq!(input.read_i16().unwrap(), -123);
            assert_eq!(input.read_u16().unwrap(), 0xCAFE);
            assert_eq!(input.read_i32().unwrap(), -7);
            assert_eq!(input.read_u64().unwrap(), 0x0102_0304_0506_0708);
            assert_eq!(input.read_value::<u32>().unwrap(), 99);
            let mut values = [0i16; 2];
            input.read_values(&mut values).unwrap();
            assert_eq!(values, [3, 4]);
            assert_eq!(input.read_string().unwrap(), "abc");
            assert!(input.file().is_ok());
            input.rewind().unwrap();
            let h = input.hash().unwrap();
            let bytes = std::fs::read(&name).unwrap();
            assert_eq!(h, crate::util::hash::file_hash(&bytes));
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_output_file_compressors_roundtrip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-compressor-out-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        {
            let mut out = OutputFile::new(&name, Compressor::Zlib, "w+b").unwrap();
            out.write_raw(b"gzip-output").unwrap();
            out.close().unwrap();
        }
        let compressed = std::fs::read(&path).unwrap();
        assert_eq!(detect_compressor(&compressed), Compressor::Zlib);
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"gzip-output");
        input.close().unwrap();

        {
            let mut out = OutputFile::new(&name, Compressor::Zstd, "w+b").unwrap();
            out.write_raw(b"zstd-output").unwrap();
            out.close().unwrap();
        }
        let compressed = std::fs::read(&path).unwrap();
        assert_eq!(detect_compressor(&compressed), Compressor::Zstd);
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zstd-output");
        input.close().unwrap();

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_input_file_temp_constructors_with_compressors() {
        let mut tmp = TempFile::new(true).unwrap();
        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"temp-zlib").unwrap();
        tmp.write_raw(&zlib_encoder.finish().unwrap()).unwrap();
        let mut input = InputFile::from_temp_file(&mut tmp, 0, Compressor::Zlib).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"temp-zlib");
        input.close().unwrap();
        tmp.close().unwrap();

        let data = TempFileData::init(false).unwrap();
        {
            let mut out = OutputFile::from_temp_file_data(&data, Compressor::Zstd, "w+b").unwrap();
            out.write_raw(b"temp-zstd").unwrap();
            out.close().unwrap();
        }
        let mut input = InputFile::from_temp_file_data(&data, 0, Compressor::Zstd).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"temp-zstd");
        input.close().unwrap();
        std::fs::remove_file(&data.name).unwrap();
    }

    #[test]
    fn test_input_file_autodetects_zlib_and_gzip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-compressed-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"zlib-data").unwrap();
        std::fs::write(&path, zlib_encoder.finish().unwrap()).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zlib-data");
        input.rewind().unwrap();
        let mut buf = [0u8; 4];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zlib");
        input.close().unwrap();

        let mut first = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        first.write_all(b"one").unwrap();
        let mut concatenated = first.finish().unwrap();
        let mut second =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        second.write_all(b"two").unwrap();
        concatenated.extend_from_slice(&second.finish().unwrap());
        std::fs::write(&path, &concatenated).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 6];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"onetwo");
        input.close().unwrap();

        let mut gzip_encoder =
            flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip_encoder.write_all(b"gzip-data").unwrap();
        let gzip_compressed = gzip_encoder.finish().unwrap();
        std::fs::write(&path, &gzip_compressed).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"gzip-data");
        input.close().unwrap();

        std::fs::write(&path, zstd_compress_for_test(b"zstd-data")).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zstd-data");
        input.close().unwrap();

        std::fs::write(&path, &gzip_compressed).unwrap();
        let mut input = InputFile::new(&name, InputFile::NO_AUTODETECT).unwrap();
        let mut header = [0u8; 2];
        input.read_raw(&mut header).unwrap();
        assert_eq!(header, [0x1f, 0x8b]);
        input.close().unwrap();

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_input_file_from_output_file() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-from-output-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
            out.write_raw(b"from-output").unwrap();
            let mut input = InputFile::from_output_file(&mut out, 0).unwrap();
            let mut buf = [0u8; 11];
            input.read_raw(&mut buf).unwrap();
            assert_eq!(&buf, b"from-output");
            assert!(input.temp_file);
            input.close().unwrap();
            out.close().unwrap();
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_input_file_from_unlinked_temp_file() {
        let mut tmp = TempFile::new(true).unwrap();
        #[cfg(not(unix))]
        let tmp_name = tmp.file_name().to_string();
        assert_eq!(tmp.unlinked(), cfg!(unix));
        tmp.write_raw(b"hidden-temp").unwrap();
        let mut input = InputFile::from_temp_file(&mut tmp, 0, Compressor::None).unwrap();
        assert!(input.temp_file);
        assert_eq!(input.unlinked, cfg!(unix));
        #[cfg(unix)]
        assert!(!std::path::Path::new(&input.file_name).exists());
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"hidden-temp");
        input.close().unwrap();
        tmp.close().unwrap();
        #[cfg(not(unix))]
        std::fs::remove_file(tmp_name).unwrap();
    }

    #[test]
    fn test_text_input_file_getline_putback() {
        let stream = VecStream::from_vec(b"a\r\nb\n".to_vec());
        let de = Deserializer::new(stream);
        let mut text = TextInputFile::new(de, b'\n');
        text.getline().unwrap();
        assert_eq!(text.line, "a");
        assert_eq!(text.line_count, 1);
        text.putback_line();
        text.getline().unwrap();
        assert_eq!(text.line, "a");
        text.getline().unwrap();
        assert_eq!(text.line, "b");
        assert!(!text.eof());
        text.getline().unwrap();
        assert_eq!(text.line, "");
        assert!(text.eof());
        text.rewind().unwrap();
        assert_eq!(text.line_count, 0);
    }

    #[test]
    fn test_text_input_file_constructors() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-text-input-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        std::fs::write(&path, b"one\r\ntwo\n").unwrap();

        let mut text = TextInputFile::from_file_name(&name, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "one");
        text.getline().unwrap();
        assert_eq!(text.line, "two");
        text.into_inner().close().unwrap();

        let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
        out.write_raw(b"out\n").unwrap();
        let mut text = TextInputFile::from_output_file(&mut out, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "out");
        text.into_inner().close().unwrap();
        out.close().unwrap();

        let data = TempFileData::init(false).unwrap();
        let mut tmp = TempFile::from_temp_file_data(&data).unwrap();
        tmp.write_raw(b"tmp\n").unwrap();
        let mut text = TextInputFile::from_temp_file(&mut tmp, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "tmp");
        text.into_inner().close().unwrap();
        tmp.close().unwrap();
        std::fs::remove_file(&data.name).unwrap();

        let mut unlinked = TempFile::new(true).unwrap();
        unlinked.write_raw(b"unlinked\n").unwrap();
        let mut text = TextInputFile::from_temp_file(&mut unlinked, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "unlinked");
        text.into_inner().close().unwrap();
        unlinked.close().unwrap();

        std::fs::remove_file(path).unwrap();
    }
}

pub mod std_stream;

/// A Source std::exception diagnostic can contain arbitrary accession bytes.
/// Display is a convenience for String APIs; native reporting uses the
/// retained bytes and never round-trips the diagnostic through UTF-8.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawDiagnosticError(pub Vec<u8>);

impl std::fmt::Display for RawDiagnosticError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(self.message_bytes()))
    }
}

impl std::error::Error for RawDiagnosticError {}

impl RawDiagnosticError {
    pub fn message_bytes(&self) -> &[u8] {
        self.0.split(|&byte| byte == 0).next().unwrap_or_default()
    }
    pub fn into_io(self) -> std::io::Error {
        std::io::Error::other(self)
    }

    pub fn find<'a>(error: &'a (dyn std::error::Error + 'static)) -> Option<&'a [u8]> {
        let mut current = Some(error);
        while let Some(error) = current {
            if let Some(raw) = error.downcast_ref::<Self>() {
                return Some(raw.message_bytes());
            }
            if let Some(io) = error.downcast_ref::<std::io::Error>() {
                if let Some(inner) = io.get_ref() {
                    current = Some(inner);
                    continue;
                }
            }
            current = error.source();
        }
        None
    }
}

/// Two owners of the same concrete error: its source chain and raw diagnostic
/// remain available both to ordinary propagation and noexcept owner cleanup.
#[derive(Debug)]
struct SharedIoError(std::sync::Arc<std::io::Error>);
impl std::fmt::Display for SharedIoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}
impl std::error::Error for SharedIoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(self.0.as_ref())
    }
}

/// Invocation-owned transport for fatal concrete I/O through translated
/// backend traits whose public error type remains String. Call capture only
/// when that failure is immediately propagated; finish restores its type.
#[derive(Clone, Default, Debug)]
pub struct IoErrorOwner {
    primary: std::sync::Arc<std::sync::Mutex<Option<std::io::Error>>>,
    source_destructors: bool,
}
impl IoErrorOwner {
    /// Opt the actual Source command owner into noexcept destructor behavior.
    /// Ordinary reusable library owners retain fallible extension methods.
    pub(crate) fn with_source_destructors(&self) -> Self {
        Self {
            primary: self.primary.clone(),
            source_destructors: true,
        }
    }
    pub(crate) fn source_destructors(&self) -> bool {
        self.source_destructors
    }
    pub(crate) fn cleanup_owner(&self) -> Self {
        Self {
            primary: Default::default(),
            source_destructors: self.source_destructors,
        }
    }
    pub(crate) fn emit_pending_perror(&self) -> IoResult<()> {
        let mut pending = self
            .primary
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if let Some(error) = pending.take() {
            let emitted = match emit_source_perror(&error) {
                Ok(emitted) => emitted,
                Err(conversion) => {
                    *pending = Some(error);
                    return Err(conversion);
                }
            };
            *pending = Some(if emitted {
                PerrorEmittedError::wrap(error)
            } else {
                error
            });
        }
        Ok(())
    }
    /// Save the actual typed failure before a Source owner can be destroyed.
    /// Returning another shared handle preserves the public io::Result API.
    pub(crate) fn capture_io<T>(&self, result: std::io::Result<T>) -> std::io::Result<T> {
        result.map_err(|error| {
            let kind = error.kind();
            let error = std::sync::Arc::new(error);
            *self
                .primary
                .lock()
                .unwrap_or_else(|error| error.into_inner()) =
                Some(std::io::Error::new(kind, SharedIoError(error.clone())));
            std::io::Error::new(kind, SharedIoError(error))
        })
    }
    pub fn capture<T>(&self, result: std::io::Result<T>) -> Result<T, String> {
        result.map_err(|error| {
            let message = error.to_string();
            *self
                .primary
                .lock()
                .unwrap_or_else(|error| error.into_inner()) = Some(error);
            message
        })
    }
    pub fn finish<T>(&self, result: Result<T, String>) -> std::io::Result<T> {
        self.finish_with_kind(result, std::io::ErrorKind::Other)
    }
    /// Restore captured I/O errors while retaining the caller's existing kind
    /// for validation errors that originated in a String-returning adapter.
    pub fn finish_with_kind<T>(
        &self,
        result: Result<T, String>,
        fallback_kind: std::io::ErrorKind,
    ) -> std::io::Result<T> {
        let error = self
            .primary
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .take();
        let result = result.map_err(|message| {
            error.unwrap_or_else(|| std::io::Error::new(fallback_kind, message))
        });
        if self.source_destructors {
            self.capture_io(result)
        } else {
            result
        }
    }
}
