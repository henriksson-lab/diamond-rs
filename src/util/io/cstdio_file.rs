//! Owned standard-C stdio operations corresponding to Source FileSink.
//!
//! The FILE layout and ABI come exclusively from the maintained libc crate.
//! Named and transferred descriptors have one owner; legacy Rust File access
//! uses an on-demand duplicate, never an alias owning the original descriptor.

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StdioOperation {
    Write,
    Read,
    Close,
    Seek,
    Tell,
}

/// A diagnostic captured at the actual failed C operation, before cleanup can
/// change errno. Reusable I/O does not print; the Source command boundary owns
/// perror/main reporting, including the separate noexcept-destructor policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StdioError {
    pub operation: StdioOperation,
    pub path: std::path::PathBuf,
    pub raw_os_error: Option<i32>,
    pub perror: Option<Vec<u8>>,
    pub message: Vec<u8>,
}

impl StdioError {
    pub fn message_bytes(&self) -> &[u8] {
        self.message
            .split(|&byte| byte == 0)
            .next()
            .unwrap_or_default()
    }

    pub fn find<'a>(error: &'a (dyn std::error::Error + 'static)) -> Option<&'a Self> {
        let mut current = Some(error);
        while let Some(error) = current {
            if let Some(error) = error.downcast_ref::<Self>() {
                return Some(error);
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

impl fmt::Display for StdioError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&String::from_utf8_lossy(self.message_bytes()))
    }
}
impl std::error::Error for StdioError {}

#[cfg(unix)]
mod posix {
    use super::{StdioError, StdioOperation};
    use crate::util::io::{IoError, IoResult};
    use std::ffi::{CStr, CString};
    use std::fs::File;
    use std::io::SeekFrom;
    use std::os::fd::{FromRawFd, IntoRawFd};
    use std::os::unix::ffi::OsStrExt;
    use std::path::Path;
    use std::ptr::NonNull;

    /// The C FILE is the sole owner of the transferred descriptor. No raw FILE
    /// pointer is exposed. All stream operations require exclusive Rust access.
    #[derive(Debug)]
    pub(crate) struct OwnedCstdio {
        stream: Option<NonNull<libc::FILE>>,
        name: std::path::PathBuf,
        file_view: Option<File>,
        view_dirty: bool,
    }

    // SAFETY: FILE ownership is unique, no borrowed FILE escapes this type,
    // and every operation on it requires &mut self. Moving or sharing an
    // immutable owner cannot race a stream operation. libc's FILE additionally
    // supplies its normal internal stdio synchronization.
    unsafe impl Send for OwnedCstdio {}
    unsafe impl Sync for OwnedCstdio {}

    impl OwnedCstdio {
        /// Distinct Source File facade uses fopen, not FileSource's open+fdopen.
        pub(crate) fn open_file_role_path(name: &Path, mode: &str) -> IoResult<Self> {
            let path = CString::new(crate::util::io::source_c_path(name).as_os_str().as_bytes())
                .map_err(|_| {
                    IoError::file_open(name, std::io::Error::from_raw_os_error(libc::EINVAL))
                })?;
            let mode = CString::new(mode).map_err(|_| {
                IoError::file_open(name, std::io::Error::from_raw_os_error(libc::EINVAL))
            })?;
            let stream = unsafe { libc::fopen(path.as_ptr(), mode.as_ptr()) };
            let Some(stream) = NonNull::new(stream) else {
                return Err(IoError::file_open(name, std::io::Error::last_os_error()));
            };
            Ok(Self {
                stream: Some(stream),
                name: name.to_path_buf(),
                file_view: None,
                view_dirty: false,
            })
        }
        /// Source File::read calls fread once and inspects only its count.
        pub(crate) fn close_file_role(&mut self) -> IoResult<()> {
            self.file_view.take();
            if let Some(stream) = self.stream.take() {
                unsafe {
                    libc::fclose(stream.as_ptr());
                }
            }
            self.view_dirty = false;
            Ok(())
        }
        pub(crate) fn read_file_role(&mut self, bytes: &mut [u8]) -> IoResult<(usize, i32)> {
            self.sync_view()?;
            let pointer = self.pointer()?;
            let count =
                unsafe { libc::fread(bytes.as_mut_ptr().cast(), 1, bytes.len(), pointer.as_ptr()) };
            let saved = std::io::Error::last_os_error().raw_os_error().unwrap_or(0);
            Ok((count, saved))
        }
        pub(crate) fn set_file_role_buffer(&mut self, size: usize) -> IoResult<()> {
            let pointer = self.pointer()?;
            if unsafe { libc::setvbuf(pointer.as_ptr(), std::ptr::null_mut(), libc::_IOFBF, size) }
                != 0
            {
                let saved =
                    super::super::FileOpenError::new(&self.name, std::io::Error::last_os_error());
                let mut message = b"Error setting buffer size for temporary file ".to_vec();
                message.extend_from_slice(self.name.as_os_str().as_bytes());
                message.extend_from_slice(b". ");
                message.extend_from_slice(&saved.reason);
                return Err(IoError::RawDetails(super::super::RawDiagnosticError(
                    message,
                )));
            }
            Ok(())
        }

        pub(crate) fn open(name: &str, mode: &str, buffer_size: usize) -> IoResult<Self> {
            Self::open_path(Path::new(name), mode, buffer_size)
        }
        pub(crate) fn open_path(name: &Path, mode: &str, buffer_size: usize) -> IoResult<Self> {
            let flags = match mode {
                "rb" => libc::O_RDONLY,
                "wb" => libc::O_WRONLY | libc::O_CREAT | libc::O_TRUNC,
                "r+b" => libc::O_RDWR,
                "w+b" => libc::O_RDWR | libc::O_CREAT | libc::O_TRUNC,
                _ => return Err(IoError::Other("Invalid fopen mode.".into())),
            };
            let path = CString::new(crate::util::io::source_c_path(name).as_os_str().as_bytes())
                .map_err(|_| {
                    IoError::file_open(name, std::io::Error::from_raw_os_error(libc::EINVAL))
                })?;
            // SAFETY: path is NUL terminated. C promotes the variadic mode
            // argument to an int-width type (Darwin's mode_t is only u16).
            // Supply it even without O_CREAT, as Source's POSIX_OPEN does.
            let descriptor = unsafe { libc::open(path.as_ptr(), flags, 0o666 as libc::c_uint) };
            if descriptor < 0 {
                return Err(IoError::file_open(name, std::io::Error::last_os_error()));
            }
            // SAFETY: successful open returned a uniquely owned descriptor.
            let file = unsafe { File::from_raw_fd(descriptor) };
            Self::from_file_path(name, file, mode, buffer_size)
        }

        pub(crate) fn from_file(
            name: &str,
            file: File,
            mode: &str,
            buffer_size: usize,
        ) -> IoResult<Self> {
            Self::from_file_path(Path::new(name), file, mode, buffer_size)
        }
        pub(crate) fn from_file_path(
            name: &Path,
            file: File,
            mode: &str,
            buffer_size: usize,
        ) -> IoResult<Self> {
            if !matches!(mode, "rb" | "wb" | "r+b" | "w+b") {
                return Err(IoError::Other("Invalid fopen mode.".into()));
            }
            let mode = CString::new(mode).expect("validated stdio mode");
            let descriptor = file.into_raw_fd();
            // SAFETY: fdopen receives ownership of the unique open descriptor
            // on success and the maintained libc declaration supplies its ABI.
            let stream = unsafe { libc::fdopen(descriptor, mode.as_ptr()) };
            let Some(stream) = NonNull::new(stream) else {
                let error = std::io::Error::last_os_error();
                // fdopen does not own the descriptor when it fails.
                unsafe { libc::close(descriptor) };
                return Err(IoError::file_open(name, error));
            };
            let result = Self {
                stream: Some(stream),
                name: name.to_path_buf(),
                file_view: None,
                view_dirty: false,
            };
            if buffer_size != 0 {
                // SAFETY: the fresh stream has performed no I/O. A null buffer
                // asks libc to own allocation, exactly as Source setvbuf does.
                if unsafe {
                    libc::setvbuf(
                        stream.as_ptr(),
                        std::ptr::null_mut(),
                        libc::_IOFBF,
                        buffer_size,
                    )
                } != 0
                {
                    return Err(IoError::Other("Error calling setvbuf.".into()));
                }
            }
            Ok(result)
        }

        fn pointer(&self) -> IoResult<NonNull<libc::FILE>> {
            self.stream
                .ok_or_else(|| IoError::path_message(b"Error writing file ", &self.name))
        }

        fn failure(&self, operation: StdioOperation, saved: std::io::Error) -> IoError {
            let raw_os_error = saved.raw_os_error();
            let reason = if let Some(code) = raw_os_error {
                // strerror consumes the saved code, not a subsequently sampled
                // errno. Copy its bytes before another operation can reuse it.
                let ptr = unsafe { libc::strerror(code) };
                if ptr.is_null() {
                    saved.to_string().into_bytes()
                } else {
                    unsafe { CStr::from_ptr(ptr) }.to_bytes().to_vec()
                }
            } else {
                saved.to_string().into_bytes()
            };
            let mut message = match operation {
                StdioOperation::Write => b"Error writing file ".to_vec(),
                StdioOperation::Read => b"Error reading file ".to_vec(),
                StdioOperation::Close => b"Error closing file ".to_vec(),
                StdioOperation::Seek => b"Error calling fseek.".to_vec(),
                StdioOperation::Tell => b"Error calling ftell.".to_vec(),
            };
            if matches!(
                operation,
                StdioOperation::Write | StdioOperation::Read | StdioOperation::Close
            ) {
                message.extend_from_slice(self.name.as_os_str().as_bytes());
            }
            IoError::StdioDetails(StdioError {
                operation,
                path: self.name.clone(),
                raw_os_error,
                perror: Some(reason),
                message,
            })
        }

        /// A legacy raw descriptor view may have moved the shared OS offset.
        /// Reposition FILE absolutely from that actual offset, not from a
        /// potentially cached FILE logical position.
        fn sync_view(&mut self) -> IoResult<()> {
            if !self.view_dirty {
                return Ok(());
            }
            let pointer = self.pointer()?;
            // A transferred input reader owns a separate descriptor, but dup
            // preserves this same open-file description. Query its real offset
            // without creating another owned File or trusting FILE's cache.
            let descriptor = unsafe { libc::fileno(pointer.as_ptr()) };
            let position = unsafe { libc::lseek(descriptor, 0, libc::SEEK_CUR) };
            if position < 0 {
                return Err(self.failure(StdioOperation::Tell, std::io::Error::last_os_error()));
            }
            if unsafe { libc::fseek(pointer.as_ptr(), position as libc::c_long, libc::SEEK_SET) }
                != 0
            {
                return Err(self.failure(StdioOperation::Seek, std::io::Error::last_os_error()));
            }
            self.view_dirty = false;
            Ok(())
        }

        /// Exception-disabled Source ofstream transport. The caller owns its
        /// sticky failbit; short writes do not materialize errno diagnostics.
        pub(crate) fn write_text_output(&mut self, bytes: &[u8]) -> IoResult<bool> {
            self.sync_view()?;
            let pointer = self.pointer()?;
            Ok(
                unsafe { libc::fwrite(bytes.as_ptr().cast(), 1, bytes.len(), pointer.as_ptr()) }
                    == bytes.len(),
            )
        }
        /// Flush exactly this text stream at a Source endl boundary.
        pub(crate) fn flush_text_output(&mut self) -> IoResult<bool> {
            self.sync_view()?;
            let pointer = self.pointer()?;
            Ok(unsafe { libc::fflush(pointer.as_ptr()) } == 0)
        }

        pub(crate) fn write_bytes(&mut self, bytes: &[u8]) -> IoResult<()> {
            self.sync_view()?;
            let pointer = self.pointer()?;
            // Source makes one fwrite, including for a zero-sized write. Do
            // not retry a short result through Rust's write_all abstraction.
            let written =
                unsafe { libc::fwrite(bytes.as_ptr().cast(), 1, bytes.len(), pointer.as_ptr()) };
            if written != bytes.len() {
                return Err(self.failure(StdioOperation::Write, std::io::Error::last_os_error()));
            }
            Ok(())
        }

        pub(crate) fn read_bytes(&mut self, bytes: &mut [u8]) -> IoResult<usize> {
            self.sync_view()?;
            let pointer = self
                .stream
                .ok_or_else(|| IoError::path_message(b"Error reading file ", &self.name))?;
            let count =
                unsafe { libc::fread(bytes.as_mut_ptr().cast(), 1, bytes.len(), pointer.as_ptr()) };
            if count != bytes.len() && unsafe { libc::feof(pointer.as_ptr()) } == 0 {
                return Err(self.failure(StdioOperation::Read, std::io::Error::last_os_error()));
            }
            Ok(count)
        }

        pub(crate) fn eof(&mut self) -> IoResult<bool> {
            let pointer = self
                .stream
                .ok_or_else(|| IoError::path_message(b"Error reading file ", &self.name))?;
            Ok(unsafe { libc::feof(pointer.as_ptr()) } != 0)
        }

        fn input_error(&self, error: IoError, operation: StdioOperation) -> IoError {
            match error {
                IoError::StdioDetails(mut error) => {
                    error.message = match operation {
                        StdioOperation::Seek => {
                            let mut bytes = b"Error executing seek on file ".to_vec();
                            bytes.extend_from_slice(self.name.as_os_str().as_bytes());
                            bytes
                        }
                        StdioOperation::Tell => {
                            let mut bytes = b"Error executing ftell on stream ".to_vec();
                            bytes.extend_from_slice(self.name.as_os_str().as_bytes());
                            bytes
                        }
                        _ => error.message,
                    };
                    IoError::StdioDetails(error)
                }
                error => error,
            }
        }

        pub(crate) fn input_seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
            self.seek(position, origin)
                .map_err(|error| self.input_error(error, StdioOperation::Seek))
        }

        pub(crate) fn input_position(&mut self) -> IoResult<i64> {
            self.position()
                .map_err(|error| self.input_error(error, StdioOperation::Tell))
        }

        pub(crate) fn close(&mut self) -> IoResult<()> {
            // The duplicate must not keep the stream descriptor alive after
            // Source's explicit close; clients' independent clones stay theirs.
            self.file_view.take();
            self.view_dirty = false;
            let Some(pointer) = self.stream.take() else {
                return Ok(());
            };
            // fclose invalidates FILE even on failure. Taking it first prevents
            // a later close/drop from retrying an invalid pointer.
            if unsafe { libc::fclose(pointer.as_ptr()) } != 0 {
                return Err(self.failure(StdioOperation::Close, std::io::Error::last_os_error()));
            }
            Ok(())
        }

        pub(crate) fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
            if matches!(origin, SeekFrom::Current(_)) {
                self.sync_view()?;
            }
            let pointer = self.pointer()?;
            let whence = match origin {
                SeekFrom::Start(_) => libc::SEEK_SET,
                SeekFrom::Current(_) => libc::SEEK_CUR,
                SeekFrom::End(_) => libc::SEEK_END,
            };
            if unsafe { libc::fseek(pointer.as_ptr(), position as libc::c_long, whence) } != 0 {
                return Err(self.failure(StdioOperation::Seek, std::io::Error::last_os_error()));
            }
            self.view_dirty = false;
            Ok(())
        }

        pub(crate) fn position(&mut self) -> IoResult<i64> {
            self.sync_view()?;
            let pointer = self.pointer()?;
            let result = unsafe { libc::ftell(pointer.as_ptr()) };
            if result < 0 {
                return Err(self.failure(StdioOperation::Tell, std::io::Error::last_os_error()));
            }
            Ok(result as i64)
        }

        pub(crate) fn rewind(&mut self) -> IoResult<()> {
            let pointer = self.pointer()?;
            // Source's void rewind deliberately does not throw on seek errors.
            unsafe { libc::rewind(pointer.as_ptr()) };
            self.view_dirty = false;
            Ok(())
        }

        pub(crate) fn clone_file(&mut self) -> IoResult<File> {
            self.sync_view()?;
            let pointer = self.pointer()?;
            let descriptor = unsafe { libc::fileno(pointer.as_ptr()) };
            let cloned = unsafe { libc::dup(descriptor) };
            if cloned < 0 {
                return Err(IoError::file_open(
                    &self.name,
                    std::io::Error::last_os_error(),
                ));
            }
            // Reads/seeks through this duplicate can advance the shared kernel
            // offset while libc retains a cached FILE position. Reconcile it
            // before the next owner operation, just as for the legacy view.
            self.view_dirty = true;
            // SAFETY: dup returned a new descriptor with independent ownership.
            Ok(unsafe { File::from_raw_fd(cloned) })
        }

        pub(crate) fn prepare_legacy(&mut self) -> IoResult<()> {
            self.sync_view()?;
            let pointer = self.pointer()?;
            if unsafe { libc::fflush(pointer.as_ptr()) } != 0 {
                return Err(self.failure(StdioOperation::Write, std::io::Error::last_os_error()));
            }
            self.view_dirty = true;
            Ok(())
        }

        pub(crate) fn legacy_file(&mut self) -> IoResult<&mut File> {
            let pointer = self.pointer()?;
            // Synchronize the C stream before giving an unbuffered Rust handle
            // to the same open-file description. Outer Serializer data remains
            // untouched: callers must still request its real flush/rewind.
            if unsafe { libc::fflush(pointer.as_ptr()) } != 0 {
                return Err(self.failure(StdioOperation::Write, std::io::Error::last_os_error()));
            }
            if self.file_view.is_none() {
                self.file_view = Some(self.clone_file()?);
            }
            self.view_dirty = true;
            Ok(self.file_view.as_mut().unwrap())
        }
    }

    impl Drop for OwnedCstdio {
        fn drop(&mut self) {
            // Resource-safe cleanup of an abandoned owner. Source's low-level
            // FileSink destructor does not enforce close/fatal behavior; only
            // actual higher-level noexcept destructor sites own that policy.
            self.file_view.take();
            if let Some(pointer) = self.stream.take() {
                unsafe { libc::fclose(pointer.as_ptr()) };
            }
        }
    }
}

#[cfg(windows)]
pub(crate) use super::windows_cstdio::{source_external_stdio_limit, OwnedCstdio};
#[cfg(unix)]
pub(crate) use posix::OwnedCstdio;

/// Opaque shared ownership for Source's InputFile(OutputFile&) FILE alias.
/// Every operation locks the one stream; explicit close invalidates all aliases.
#[cfg(any(unix, windows))]
#[derive(Debug, Clone)]
pub struct SharedCstdio(std::sync::Arc<std::sync::Mutex<OwnedCstdio>>);

#[cfg(any(unix, windows))]
impl SharedCstdio {
    pub(crate) fn new(owner: OwnedCstdio) -> Self {
        Self(std::sync::Arc::new(std::sync::Mutex::new(owner)))
    }

    pub(crate) fn with<T>(
        &self,
        operation: impl FnOnce(&mut OwnedCstdio) -> super::IoResult<T>,
    ) -> super::IoResult<T> {
        let mut owner = self
            .0
            .lock()
            .map_err(|_| super::IoError::Other("Poisoned stdio owner.".into()))?;
        operation(&mut owner)
    }
}

#[cfg(all(test, unix))]
mod tests {
    use super::OwnedCstdio;
    use super::StdioOperation;
    use crate::util::io::IoError;
    use std::io::{Read, Seek, SeekFrom, Write};
    use std::os::fd::AsRawFd;
    use std::os::unix::fs::MetadataExt;

    fn path(label: &str) -> std::path::PathBuf {
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        std::env::temp_dir().join(format!(
            "diamond-cstdio-{label}-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ))
    }

    #[test]
    fn transferred_descriptor_survives_unlink_and_closes_once() {
        let path = path("transfer");
        let file = std::fs::OpenOptions::new()
            .create_new(true)
            .read(true)
            .write(true)
            .open(&path)
            .unwrap();
        let descriptor = file.as_raw_fd();
        let inode = file.metadata().unwrap().ino();
        std::fs::remove_file(&path).unwrap();
        let mut owner = OwnedCstdio::from_file(path.to_str().unwrap(), file, "w+b", 0).unwrap();
        owner.write_bytes(b"same descriptor").unwrap();
        assert_eq!(owner.position().unwrap(), 15);
        owner.rewind().unwrap();
        let mut reader = owner.clone_file().unwrap();
        assert_ne!(reader.as_raw_fd(), descriptor);
        assert_eq!(reader.metadata().unwrap().ino(), inode);
        assert_eq!(reader.metadata().unwrap().nlink(), 0);
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).unwrap();
        assert_eq!(bytes, b"same descriptor");
        assert_eq!(owner.position().unwrap(), 15);
        owner.close().unwrap();
        assert_eq!(unsafe { libc::fcntl(descriptor, libc::F_GETFD) }, -1);
        // The independent owned reader remains valid after fclose.
        reader.seek(SeekFrom::Start(0)).unwrap();
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).unwrap();
        assert_eq!(bytes, b"same descriptor");
        owner.close().unwrap();
    }

    #[test]
    fn legacy_duplicate_offset_is_resynchronized_before_stdio_write() {
        let path = path("legacy-offset");
        let mut owner = OwnedCstdio::open(path.to_str().unwrap(), "w+b", 0).unwrap();
        owner.write_bytes(b"abcdef").unwrap();
        let file = owner.legacy_file().unwrap();
        file.seek(SeekFrom::Start(2)).unwrap();
        file.write_all(b"X").unwrap();
        owner.write_bytes(b"Y").unwrap();
        assert_eq!(owner.position().unwrap(), 4);
        owner.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"abXYef");
        std::fs::remove_file(path).unwrap();
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn buffered_dev_full_small_write_fails_at_close_large_write_at_write() {
        let mut small = OwnedCstdio::open("/dev/full", "wb", 0).unwrap();
        small.write_bytes(b"x").unwrap();
        assert_eq!(small.position().unwrap(), 1);
        let IoError::StdioDetails(error) = small.close().unwrap_err() else {
            panic!("close must retain actual stdio metadata")
        };
        assert_eq!(error.operation, StdioOperation::Close);
        assert_eq!(error.raw_os_error, Some(libc::ENOSPC));
        assert_eq!(error.path, std::path::Path::new("/dev/full"));
        assert_eq!(error.message_bytes(), b"Error closing file /dev/full");
        small.close().unwrap();

        let mut large = OwnedCstdio::open("/dev/full", "wb", 0).unwrap();
        let IoError::StdioDetails(error) = large.write_bytes(&[b'x'; 8192]).unwrap_err() else {
            panic!("write must retain actual stdio metadata")
        };
        assert_eq!(error.operation, StdioOperation::Write);
        assert_eq!(error.raw_os_error, Some(libc::ENOSPC));
        assert_eq!(error.message_bytes(), b"Error writing file /dev/full");
        // No retry is performed by owned cleanup, and errors are not printed
        // by this reusable primitive.
        let _ = large.close();
    }
}
