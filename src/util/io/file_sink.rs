//! Translation of `diamond/src/util/io/file_sink.{h,cpp}`.
//!
//! POSIX/macOS use maintained libc bindings; Windows uses the generated CRT
//! declarations selected for its target ABI. Both retain actual FILE buffering,
//! with no C++ dependency or raw FILE pointer exposed to library callers.

use super::{IoError, IoResult, StreamEntity};
use std::fs::{File as StdFile, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
#[cfg(any(unix, windows))]
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

#[cfg(any(unix, windows))]
static PROCESS_STDOUT: AtomicBool = AtomicBool::new(false);
#[cfg(unix)]
static STDOUT_CLAIMED: AtomicBool = AtomicBool::new(false);

/// Opt the standalone executable into Source's stdout ownership.
///
/// # Safety
/// Call before command work. This grants FileSink the Source CRT stdout
/// ownership policy; embedded callers keep duplicate descriptors.
#[doc(hidden)]
pub unsafe fn enable_process_stdout_ownership() {
    #[cfg(any(unix, windows))]
    PROCESS_STDOUT.store(true, Ordering::Release);
}

#[cfg(any(unix, windows))]
pub(crate) fn owns_process_stdio() -> bool {
    PROCESS_STDOUT.load(Ordering::Acquire)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FileMode {
    Write,
    ReadWrite,
    WriteRead,
}

impl FileMode {
    fn parse(mode: &str) -> IoResult<Self> {
        match mode {
            "wb" => Ok(Self::Write),
            "r+b" => Ok(Self::ReadWrite),
            "w+b" => Ok(Self::WriteRead),
            _ => Err(IoError::Other("Invalid fopen mode.".to_string())),
        }
    }
}

#[cfg(any(unix, windows))]
#[derive(Debug)]
pub struct FileSink {
    file: super::cstdio_file::SharedCstdio,
    file_view: Option<StdFile>,
    file_name: PathBuf,
    async_: bool,
    mtx: Mutex<()>,
}

#[cfg(any(unix, windows))]
impl FileSink {
    pub fn file_path(&self) -> &Path {
        &self.file_name
    }

    pub fn new(file_name: &str, mode: &str, async_: bool, buffer_size: usize) -> IoResult<Self> {
        Self::new_path(Path::new(file_name), mode, async_, buffer_size)
    }
    pub fn new_path(
        file_name: &Path,
        mode: &str,
        async_: bool,
        buffer_size: usize,
    ) -> IoResult<Self> {
        let _ = FileMode::parse(mode)?;
        let file = if file_name.as_os_str().is_empty() {
            #[cfg(windows)]
            {
                if owns_process_stdio() {
                    super::cstdio_file::OwnedCstdio::standard_output("", mode, buffer_size)?
                } else {
                    let file = crate::util::io::duplicate_standard_handle(
                        windows_sys::Win32::System::Console::STD_OUTPUT_HANDLE,
                    )
                    .map_err(|error| IoError::file_open(file_name, error))?;
                    super::cstdio_file::OwnedCstdio::from_file_path(
                        file_name,
                        file,
                        mode,
                        buffer_size,
                    )?
                }
            }
            #[cfg(unix)]
            {
                use std::os::fd::{AsFd, FromRawFd};
                let file = if PROCESS_STDOUT.load(Ordering::Acquire) {
                    if unsafe { libc::fcntl(libc::STDOUT_FILENO, libc::F_GETFD) } < 0 {
                        return Err(IoError::file_open(
                            file_name,
                            std::io::Error::last_os_error(),
                        ));
                    }
                    if STDOUT_CLAIMED
                        .compare_exchange(false, true, Ordering::AcqRel, Ordering::Acquire)
                        .is_err()
                    {
                        return Err(IoError::FileOpen(String::new()));
                    }
                    // SAFETY: process setup granted the unique valid FD1 owner.
                    unsafe { StdFile::from_raw_fd(libc::STDOUT_FILENO) }
                } else {
                    std::io::stdout()
                        .as_fd()
                        .try_clone_to_owned()
                        .map(StdFile::from)
                        .map_err(|error| IoError::file_open(file_name, error))?
                };
                super::cstdio_file::OwnedCstdio::from_file_path(file_name, file, mode, buffer_size)?
            }
        } else {
            super::cstdio_file::OwnedCstdio::open_path(file_name, mode, buffer_size)?
        };
        Ok(Self {
            file: super::cstdio_file::SharedCstdio::new(file),
            file_view: None,
            file_name: file_name.into(),
            async_,
            mtx: Mutex::new(()),
        })
    }

    pub fn from_file(
        file_name: &str,
        file: StdFile,
        mode: &str,
        async_: bool,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Self::from_file_path(Path::new(file_name), file, mode, async_, buffer_size)
    }
    pub fn from_file_path(
        file_name: &Path,
        file: StdFile,
        mode: &str,
        async_: bool,
        buffer_size: usize,
    ) -> IoResult<Self> {
        let _ = FileMode::parse(mode)?;
        Ok(Self {
            file: super::cstdio_file::SharedCstdio::new(
                super::cstdio_file::OwnedCstdio::from_file_path(
                    file_name,
                    file,
                    mode,
                    buffer_size,
                )?,
            ),
            file_view: None,
            file_name: file_name.into(),
            async_,
            mtx: Mutex::new(()),
        })
    }
}

#[cfg(any(unix, windows))]
impl StreamEntity for FileSink {
    fn close(&mut self) -> IoResult<()> {
        self.file_view.take();
        self.file.with(|file| file.close())
    }

    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        if self.async_ {
            let _guard = self
                .mtx
                .lock()
                .map_err(|_| IoError::path_message(b"Error writing file ", &self.file_name))?;
            self.file.with(|file| file.write_bytes(ptr))
        } else {
            self.file.with(|file| file.write_bytes(ptr))
        }
    }
    fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        self.file.with(|file| file.seek(position, origin))
    }
    fn rewind(&mut self) -> IoResult<()> {
        self.file.with(|file| file.rewind())
    }
    fn tell(&mut self) -> IoResult<i64> {
        self.file.with(|file| file.position())
    }
    fn file_path(&self) -> &Path {
        &self.file_name
    }
    fn file_name(&self) -> &str {
        self.file_name
            .to_str()
            .expect("Non-UTF8 OS path; use file_path")
    }
    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.file.with(|file| file.prepare_legacy())?;
        if self.file_view.is_none() {
            self.file_view = Some(self.file.with(|file| file.clone_file())?);
        }
        Ok(self.file_view.as_mut().unwrap())
    }
    fn shared_stdio(&mut self) -> IoResult<Option<super::cstdio_file::SharedCstdio>> {
        Ok(Some(self.file.clone()))
    }
    fn clone_file(&mut self) -> IoResult<StdFile> {
        self.file.with(|file| file.clone_file())
    }
    fn flush(&mut self) -> IoResult<()> {
        // Source's serializer flush calls fwrite, not fflush. In particular,
        // a tiny /dev/full write may succeed and fail only on explicit fclose.
        Ok(())
    }
    fn file_size(&mut self) -> IoResult<i64> {
        let position = self.tell()?;
        self.seek(0, SeekFrom::End(0))?;
        let size = self.tell()?;
        self.seek(position, SeekFrom::Start(0))?;
        Ok(size)
    }
    fn seekable(&self) -> bool {
        false
    }
}

#[cfg(not(any(unix, windows)))]
mod fallback {
    use super::*;
    #[derive(Debug)]
    pub struct FileSink {
        file: Option<StdFile>,
        file_name: PathBuf,
        async_: bool,
        is_stdout: bool,
        mtx: Mutex<()>,
    }

    impl FileSink {
        /// C++ `FileSink(file_name, mode, async, buffer_size)`.
        pub fn new(
            file_name: &str,
            mode: &str,
            async_: bool,
            _buffer_size: usize,
        ) -> IoResult<Self> {
            let file_name = Path::new(file_name);
            let mode = FileMode::parse(mode)?;

            let mut options = OpenOptions::new();
            let is_stdout = file_name.as_os_str().is_empty();
            if is_stdout {
                options.write(true);
                if mode != FileMode::Write {
                    options.read(true);
                }
            } else {
                match mode {
                    FileMode::Write => {
                        options.write(true).create(true).truncate(true);
                    }
                    FileMode::ReadWrite => {
                        options.read(true).write(true);
                    }
                    FileMode::WriteRead => {
                        options.read(true).write(true).create(true).truncate(true);
                    }
                }
            }

            let file = if is_stdout {
                #[cfg(windows)]
                {
                    crate::util::io::duplicate_standard_handle(
                        windows_sys::Win32::System::Console::STD_OUTPUT_HANDLE,
                    )
                    .map_err(|error| IoError::file_open(file_name, error))?
                }
                #[cfg(unix)]
                {
                    use std::os::fd::{AsFd, FromRawFd};
                    if PROCESS_STDOUT.load(Ordering::Acquire) {
                        // Source fdopen(1) creates a new FILE*, distinct from
                        // stdout; fclose therefore closes the actual descriptor.
                        // Validate before constructing the unique Rust owner.
                        if unsafe { libc::fcntl(libc::STDOUT_FILENO, libc::F_GETFD) } < 0
                            || STDOUT_CLAIMED
                                .compare_exchange(false, true, Ordering::AcqRel, Ordering::Acquire)
                                .is_err()
                        {
                            return Err(IoError::path_message(b"Error opening file ", file_name));
                        }
                        // SAFETY: the process entry point granted ownership, the
                        // descriptor is valid, and the atomic claim prevents a
                        // second File owner from adopting it.
                        unsafe { StdFile::from_raw_fd(libc::STDOUT_FILENO) }
                    } else {
                        // Embedded/library calls cannot take the host's FD1.
                        // Clone the descriptor itself, including pipe redirections,
                        // rather than reopening a platform-specific /dev path.
                        std::io::stdout()
                            .as_fd()
                            .try_clone_to_owned()
                            .map(StdFile::from)
                            .map_err(|error| IoError::file_open(file_name, error))?
                    }
                }
                #[cfg(not(any(unix, windows)))]
                {
                    options
                        .open("/dev/stdout")
                        .map_err(|error| IoError::file_open(file_name, error))?
                }
            } else {
                options
                    .open(super::super::source_c_path(file_name).as_ref())
                    .map_err(|error| IoError::file_open(file_name, error))?
            };
            Ok(Self {
                file: Some(file),
                file_name: file_name.to_path_buf(),
                async_,
                is_stdout,
                mtx: Mutex::new(()),
            })
        }

        /// Safe equivalent of Unix C++ `FileSink(file_name, fd, mode, async,
        /// buffer_size)`, taking ownership of an already-open file.
        pub fn from_file(
            file_name: &str,
            file: StdFile,
            mode: &str,
            async_: bool,
            _buffer_size: usize,
        ) -> IoResult<Self> {
            let file_name = Path::new(file_name);
            let _ = FileMode::parse(mode)?;
            Ok(Self {
                file: Some(file),
                file_name: file_name.to_path_buf(),
                async_,
                is_stdout: false,
                mtx: Mutex::new(()),
            })
        }

        fn open_file(&mut self) -> IoResult<&mut StdFile> {
            self.file
                .as_mut()
                .ok_or_else(|| IoError::path_message(b"Error writing file ", &self.file_name))
        }
    }

    impl StreamEntity for FileSink {
        fn close(&mut self) -> IoResult<()> {
            // On POSIX Source fdopen(1) != stdout and fclose owns descriptor 1.
            // A library sink closes only its duplicate; CLI sinks close FD1 itself.
            // Windows instead uses the existing stdout FILE* and leaves it open.
            if !self.is_stdout || !cfg!(windows) {
                self.file.take();
            }
            Ok(())
        }

        fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
            if self.async_ {
                let _guard = self
                    .mtx
                    .lock()
                    .map_err(|_| IoError::path_message(b"Error writing file ", &self.file_name))?;
                self.file
                    .as_mut()
                    .ok_or_else(|| IoError::path_message(b"Error writing file ", &self.file_name))?
                    .write_all(ptr)
                    .map_err(|_| IoError::path_message(b"Error writing file ", &self.file_name))
            } else {
                self.open_file()?
                    .write_all(ptr)
                    .map_err(|_| IoError::path_message(b"Error writing file ", &self.file_name))
            }
        }

        fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
            let target = match origin {
                SeekFrom::Start(_) if position < 0 => {
                    return Err(IoError::Other("Error calling fseek.".to_string()));
                }
                SeekFrom::Start(_) => SeekFrom::Start(position as u64),
                SeekFrom::End(_) => SeekFrom::End(position),
                SeekFrom::Current(_) => SeekFrom::Current(position),
            };
            self.file
                .as_mut()
                .ok_or_else(|| IoError::Other("Error calling fseek.".to_string()))?
                .seek(target)
                .map(|_| ())
                .map_err(|_| IoError::Other("Error calling fseek.".to_string()))
        }

        fn rewind(&mut self) -> IoResult<()> {
            self.seek(0, SeekFrom::Start(0))
        }

        fn tell(&mut self) -> IoResult<i64> {
            let file_name = self.file_name.clone();
            self.file
                .as_mut()
                .ok_or_else(|| {
                    IoError::Other(format!("Error executing ftell on stream {file_name}"))
                })?
                .stream_position()
                .map(|position| position as i64)
                .map_err(|_| IoError::Other("Error calling ftell.".to_string()))
        }

        fn file_path(&self) -> &Path {
            &self.file_name
        }
        fn file_name(&self) -> &str {
            self.file_name
                .to_str()
                .expect("Non-UTF8 OS path; use file_path")
        }

        fn file(&mut self) -> IoResult<&mut StdFile> {
            self.open_file()
        }

        fn flush(&mut self) -> IoResult<()> {
            let file_name = self.file_name.clone();
            self.file
                .as_mut()
                .map(|file| {
                    file.flush()
                        .map_err(|_| IoError::path_message(b"Error writing file ", &file_name))
                })
                .unwrap_or(Ok(()))
        }

        fn file_size(&mut self) -> IoResult<i64> {
            let position = self.tell()?;
            self.seek(0, SeekFrom::End(0))?;
            let size = self.tell()?;
            self.seek(position, SeekFrom::Start(0))?;
            Ok(size)
        }

        // FileSink does not opt into `StreamEntity(true)` in C++.
        fn seekable(&self) -> bool {
            false
        }
    }
}
#[cfg(not(any(unix, windows)))]
pub use fallback::FileSink;

#[cfg(test)]
mod tests {
    use super::*;

    fn path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-file-sink-{name}-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ))
    }

    #[test]
    fn c_path_output_writes_prefix_while_retaining_full_owner() {
        let path = path("c-prefix");
        let mut full = path.clone();
        full.as_mut_os_string().push("\0hidden-suffix");
        let mut sink = FileSink::new_path(&full, "wb", false, 4096).unwrap();
        assert_eq!(sink.file_path(), full.as_path());
        sink.write(b"source-prefix").unwrap();
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"source-prefix");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn writes_seeks_rewinds_and_preserves_position_when_sizing() {
        let path = path("binary");
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        let mut sink = FileSink::new(&name, "w+b", false, 4096).unwrap();
        assert!(!sink.seekable());
        sink.write(b"abcdef").unwrap();
        sink.seek(2, SeekFrom::Start(0)).unwrap();
        sink.write(b"XY").unwrap();
        assert_eq!(sink.tell().unwrap(), 4);
        assert_eq!(sink.file_size().unwrap(), 6);
        assert_eq!(sink.tell().unwrap(), 4);
        sink.rewind().unwrap();
        assert_eq!(sink.tell().unwrap(), 0);
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"abXYef");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn modes_match_fopen_create_and_truncate_behavior() {
        let path = path("modes");
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, b"old").unwrap();
        let name = path.to_string_lossy();
        let mut update = FileSink::new(&name, "r+b", false, 0).unwrap();
        update.write(b"N").unwrap();
        update.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"Nld");

        let mut truncate = FileSink::new(&name, "wb", false, 0).unwrap();
        truncate.write(b"x").unwrap();
        truncate.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"x");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn async_write_and_closed_handle_report_expected_state() {
        let path = path("close");
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy().into_owned();
        let mut sink = FileSink::new(&name, "wb", true, 0).unwrap();
        sink.write(b"data").unwrap();
        sink.close().unwrap();
        sink.close().unwrap();
        let error = sink.write(b"x").unwrap_err();
        let mut expected = b"Error writing file ".to_vec();
        expected.extend_from_slice(name.as_bytes());
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error),
            Some(expected.as_slice())
        );
        assert_eq!(std::fs::read(&path).unwrap(), b"data");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn invalid_mode_and_negative_absolute_seek_are_rejected() {
        let path = path("errors");
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        assert_eq!(
            FileSink::new(&name, "append", false, 0).unwrap_err(),
            IoError::Other("Invalid fopen mode.".to_string())
        );
        let mut sink = FileSink::new(&name, "w+b", false, 0).unwrap();
        assert_eq!(
            sink.seek(-1, SeekFrom::Start(0)).unwrap_err().to_string(),
            "Error calling fseek."
        );
        sink.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn file_mode_accepts_only_upstream_spellings() {
        assert!(matches!(FileMode::parse("wb"), Ok(FileMode::Write)));
        assert!(matches!(FileMode::parse("r+b"), Ok(FileMode::ReadWrite)));
        assert!(matches!(FileMode::parse("w+b"), Ok(FileMode::WriteRead)));
        for rejected in ["", "rb", "ab", "append"] {
            assert_eq!(
                FileMode::parse(rejected).unwrap_err(),
                IoError::Other("Invalid fopen mode.".to_string())
            );
        }
    }

    #[test]
    fn library_stdout_close_preserves_the_process_handle() {
        let mut sink = FileSink::new("", "wb", false, 0).unwrap();
        sink.close().unwrap();
        assert_eq!(sink.file_name(), "");
        #[cfg(unix)]
        {
            use std::os::fd::AsFd;
            assert!(sink.file().is_err());
            assert!(std::io::stdout().as_fd().try_clone_to_owned().is_ok());
        }
        #[cfg(windows)]
        {
            // Embedded sinks own a duplicate. Actual MSVC stdout is adopted
            // only by the standalone startup profile and remains nonclosing.
            assert!(sink.file().is_err());
            assert!(crate::util::io::duplicate_standard_handle(
                windows_sys::Win32::System::Console::STD_OUTPUT_HANDLE,
            )
            .is_ok());
        }
    }

    #[cfg(unix)]
    #[test]
    fn from_file_takes_ownership_and_writes_exact_bytes() {
        let path = path("fd");
        let _ = std::fs::remove_file(&path);
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        let name = path.to_string_lossy();
        let mut sink = FileSink::from_file(&name, file, "w+b", false, 0).unwrap();
        sink.write(b"fd").unwrap();
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"fd");
        std::fs::remove_file(path).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn from_file_validates_write_mode_without_retruncating_open_file() {
        let path = path("from-file-no-truncate");
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, b"existing").unwrap();
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .open(&path)
            .unwrap();
        file.seek(SeekFrom::End(0)).unwrap();
        let name = path.to_string_lossy();
        let mut sink = FileSink::from_file(&name, file, "w+b", false, 0).unwrap();
        sink.write(b"+").unwrap();
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"existing+");
        std::fs::remove_file(path).unwrap();
    }
}
