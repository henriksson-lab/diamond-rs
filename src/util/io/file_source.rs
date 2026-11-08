//! Translation of `diamond/src/util/io/file_source.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use super::{IoError, IoResult, StreamEntity};

#[cfg(any(unix, windows))]
static PROCESS_STDIN: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Grant Source-profile standalone input constructors ownership of CRT stdin.
/// Embedded callers instead retain a duplicate of their host input descriptor.
///
/// # Safety
/// Call only at standalone process startup before input is read.
#[doc(hidden)]
pub unsafe fn enable_process_stdin_ownership() {
    #[cfg(any(unix, windows))]
    PROCESS_STDIN.store(true, std::sync::atomic::Ordering::Release);
}

#[cfg(any(unix, windows))]
fn owns_process_stdin() -> bool {
    PROCESS_STDIN.load(std::sync::atomic::Ordering::Acquire)
}

#[cfg(any(unix, windows))]
#[derive(Debug)]
pub struct FileSource {
    owner: Option<super::cstdio_file::SharedCstdio>,
    raw: Option<fallback::RawFileSource>,
    file_view: Option<StdFile>,
    file_name: PathBuf,
    seekable: bool,
}

#[cfg(any(unix, windows))]
impl FileSource {
    pub fn file_path(&self) -> &Path {
        &self.file_name
    }

    pub fn new(file_name: &str) -> IoResult<Self> {
        Self::new_path(Path::new(file_name))
    }
    pub fn new_path(file_name: &Path) -> IoResult<Self> {
        let is_stdin =
            file_name.as_os_str().is_empty() || file_name.as_os_str() == std::ffi::OsStr::new("-");
        // Source POSIX stats before open. Its MSVC constructor uses fopen and
        // leaves StreamEntity(true), even for redirected stdin/special files.
        #[cfg(unix)]
        let seekable = if is_stdin {
            false
        } else {
            std::fs::metadata(super::source_c_path(file_name).as_ref())
                .map_err(|error| IoError::file_open(file_name, error))?
                .is_file()
        };
        #[cfg(all(windows, not(target_env = "msvc")))]
        let seekable = if is_stdin {
            false
        } else {
            super::cstdio_file::OwnedCstdio::source_seekable_path(file_name)?
        };
        #[cfg(all(windows, target_env = "msvc"))]
        let seekable = true;
        let owner = if is_stdin {
            #[cfg(windows)]
            {
                if owns_process_stdin() {
                    super::cstdio_file::OwnedCstdio::standard_input(
                        file_name.to_str().expect("stdin sentinel is UTF-8"),
                    )?
                } else {
                    let file = crate::util::io::duplicate_standard_handle(
                        windows_sys::Win32::System::Console::STD_INPUT_HANDLE,
                    )
                    .map_err(|error| IoError::file_open(file_name, error))?;
                    super::cstdio_file::OwnedCstdio::from_file_path(file_name, file, "rb", 0)?
                }
            }
            #[cfg(unix)]
            {
                use std::os::fd::{AsFd, FromRawFd};
                let file = if owns_process_stdin() {
                    // Standalone startup grants ownership of actual fd0.
                    unsafe { StdFile::from_raw_fd(libc::STDIN_FILENO) }
                } else {
                    std::io::stdin()
                        .as_fd()
                        .try_clone_to_owned()
                        .map(StdFile::from)
                        .map_err(|error| IoError::file_open(file_name, error))?
                };
                super::cstdio_file::OwnedCstdio::from_file_path(file_name, file, "rb", 0)?
            }
        } else {
            super::cstdio_file::OwnedCstdio::open_path(file_name, "rb", 0)?
        };
        Ok(Self {
            owner: Some(super::cstdio_file::SharedCstdio::new(owner)),
            raw: None,
            file_view: None,
            file_name: file_name.into(),
            seekable,
        })
    }

    /// Independent supplied descriptor constructor. Source FILE aliases use
    /// from_shared instead, preserving one libc buffer/EOF/close state.
    pub fn from_file(file_name: &str, file: StdFile) -> Self {
        Self::from_file_path(Path::new(file_name), file)
    }
    pub fn from_file_path(file_name: &Path, file: StdFile) -> Self {
        Self {
            owner: None,
            raw: Some(fallback::RawFileSource::from_file(file_name, file)),
            file_view: None,
            file_name: file_name.into(),
            seekable: false,
        }
    }

    pub fn try_from_file(file_name: &str, file: StdFile) -> IoResult<Self> {
        Self::try_from_file_path(Path::new(file_name), file)
    }
    pub fn try_from_file_path(file_name: &Path, file: StdFile) -> IoResult<Self> {
        let owner = super::cstdio_file::OwnedCstdio::from_file_path(file_name, file, "rb", 0)?;
        Ok(Self::from_shared_path(
            file_name,
            super::cstdio_file::SharedCstdio::new(owner),
        ))
    }

    /// Source's `FileSource(name, FILE*)` over an opaque shared stream. Reads,
    /// EOF, rewind and explicit close act on the same FILE as the output owner.
    /// Closing either alias invalidates it for every owner; unlike Source's
    /// dangling-pointer reuse, subsequent Rust operations return an error.
    pub fn from_shared(file_name: &str, owner: super::cstdio_file::SharedCstdio) -> Self {
        Self::from_shared_path(Path::new(file_name), owner)
    }
    pub fn from_shared_path(file_name: &Path, owner: super::cstdio_file::SharedCstdio) -> Self {
        Self {
            owner: Some(owner),
            raw: None,
            file_view: None,
            file_name: file_name.into(),
            seekable: false,
        }
    }
}

#[cfg(any(unix, windows))]
impl StreamEntity for FileSource {
    fn rewind(&mut self) -> IoResult<()> {
        if let Some(raw) = &mut self.raw {
            return raw.rewind();
        }
        self.owner.as_ref().unwrap().with(|owner| owner.rewind())
    }
    fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        if let Some(raw) = &mut self.raw {
            return raw.seek(position, origin);
        }
        self.owner
            .as_ref()
            .unwrap()
            .with(|owner| owner.input_seek(position, origin))
    }
    fn tell(&mut self) -> IoResult<i64> {
        if let Some(raw) = &mut self.raw {
            return raw.tell();
        }
        self.owner
            .as_ref()
            .unwrap()
            .with(|owner| owner.input_position())
    }
    fn read(&mut self, bytes: &mut [u8]) -> IoResult<usize> {
        if let Some(raw) = &mut self.raw {
            return raw.read(bytes);
        }
        self.owner
            .as_ref()
            .unwrap()
            .with(|owner| owner.read_bytes(bytes))
    }
    fn eof(&mut self) -> IoResult<bool> {
        if let Some(raw) = &mut self.raw {
            return raw.eof();
        }
        self.owner.as_ref().unwrap().with(|owner| owner.eof())
    }
    fn close(&mut self) -> IoResult<()> {
        if let Some(raw) = &mut self.raw {
            return raw.close();
        }
        self.file_view.take();
        self.owner.as_ref().unwrap().with(|owner| owner.close())
    }
    fn file_size(&mut self) -> IoResult<i64> {
        let position = self.tell()?;
        self.seek(0, SeekFrom::End(0))?;
        let size = self.tell()?;
        self.seek(position, SeekFrom::Start(0))?;
        Ok(size)
    }
    fn file_path(&self) -> &Path {
        &self.file_name
    }
    fn file_name(&self) -> &str {
        self.file_name
            .to_str()
            .expect("Non-UTF8 OS path; use file_path")
    }
    fn seekable(&self) -> bool {
        self.seekable
    }
    fn file(&mut self) -> IoResult<&mut StdFile> {
        if let Some(raw) = &mut self.raw {
            return raw.file();
        }
        let owner = self.owner.as_ref().unwrap();
        owner.with(|owner| owner.prepare_legacy())?;
        if self.file_view.is_none() {
            self.file_view = Some(owner.with(|owner| owner.clone_file())?);
        }
        Ok(self.file_view.as_mut().unwrap())
    }
    fn clone_file(&mut self) -> IoResult<StdFile> {
        if let Some(raw) = &mut self.raw {
            return raw.clone_file();
        }
        self.owner
            .as_ref()
            .unwrap()
            .with(|owner| owner.clone_file())
    }
    fn shared_stdio(&mut self) -> IoResult<Option<super::cstdio_file::SharedCstdio>> {
        Ok(self.owner.clone())
    }
}

mod fallback {
    use super::*;
    #[derive(Debug)]
    pub struct RawFileSource {
        file: Option<StdFile>,
        file_name: PathBuf,
        seekable: bool,
        eof: bool,
    }

    impl RawFileSource {
        /// C++ `RawFileSource(const string&)`.
        pub fn new(file_name: impl AsRef<Path>) -> IoResult<Self> {
            let file_name = file_name.as_ref();
            let is_stdin = file_name.as_os_str().is_empty()
                || file_name.as_os_str() == std::ffi::OsStr::new("-");
            let file = if is_stdin {
                #[cfg(windows)]
                {
                    crate::util::io::duplicate_standard_handle(
                        windows_sys::Win32::System::Console::STD_INPUT_HANDLE,
                    )
                    .map_err(|error| IoError::file_open(file_name, error))?
                }
                #[cfg(not(windows))]
                {
                    StdFile::open("/dev/stdin")
                        .map_err(|error| IoError::file_open(file_name, error))?
                }
            } else {
                StdFile::open(super::super::source_c_path(file_name).as_ref())
                    .map_err(|error| IoError::file_open(file_name, error))?
            };
            let seekable = !is_stdin && file.metadata().map(|meta| meta.is_file()).unwrap_or(false);
            Ok(Self {
                file: Some(file),
                file_name: file_name.to_path_buf(),
                seekable,
                eof: false,
            })
        }

        /// C++ `RawFileSource(const string&, FILE*)`.
        ///
        /// Rust takes ownership of a cloned or otherwise owned file handle rather
        /// than borrowing a raw `FILE*`.
        pub fn from_file(file_name: impl AsRef<Path>, file: StdFile) -> Self {
            let file_name = file_name.as_ref();
            Self {
                file: Some(file),
                file_name: file_name.to_path_buf(),
                seekable: false,
                eof: false,
            }
        }

        fn open_file(&mut self) -> IoResult<&mut StdFile> {
            self.file
                .as_mut()
                .ok_or_else(|| IoError::path_message(b"Error reading file ", &self.file_name))
        }
    }

    impl StreamEntity for RawFileSource {
        fn rewind(&mut self) -> IoResult<()> {
            let name = self.file_name.clone();
            self.open_file()?
                .seek(SeekFrom::Start(0))
                .map_err(|_| IoError::path_message(b"Error executing seek on file ", &name))?;
            self.eof = false;
            Ok(())
        }

        fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
            let target = match origin {
                SeekFrom::Start(_) if position < 0 => {
                    return Err(IoError::path_message(
                        b"Error executing seek on file ",
                        &self.file_name,
                    ));
                }
                SeekFrom::Start(_) => SeekFrom::Start(position as u64),
                SeekFrom::End(_) => SeekFrom::End(position),
                SeekFrom::Current(_) => SeekFrom::Current(position),
            };
            let name = self.file_name.clone();
            self.open_file()?
                .seek(target)
                .map_err(|_| IoError::path_message(b"Error executing seek on file ", &name))?;
            self.eof = false;
            Ok(())
        }

        fn tell(&mut self) -> IoResult<i64> {
            let name = self.file_name.clone();
            self.open_file()?
                .stream_position()
                .map(|position| position as i64)
                .map_err(|_| IoError::path_message(b"Error executing ftell on stream ", &name))
        }

        fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
            if ptr.is_empty() {
                return Ok(0);
            }

            let mut read = 0;
            while read < ptr.len() {
                let name = self.file_name.clone();
                match self.open_file()?.read(&mut ptr[read..]) {
                    Ok(0) => {
                        self.eof = true;
                        break;
                    }
                    Ok(count) => read += count,
                    Err(_) => return Err(IoError::path_message(b"Error reading file ", &name)),
                }
            }
            Ok(read)
        }

        fn close(&mut self) -> IoResult<()> {
            self.file.take();
            Ok(())
        }

        fn file_size(&mut self) -> IoResult<i64> {
            let position = self.tell()?;
            self.seek(0, SeekFrom::End(0))?;
            let size = self.tell()?;
            self.seek(position, SeekFrom::Start(0))?;
            Ok(size)
        }

        fn eof(&mut self) -> IoResult<bool> {
            Ok(self.eof)
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

        fn seekable(&self) -> bool {
            self.seekable
        }
    }
}
#[cfg(not(any(unix, windows)))]
pub use fallback::RawFileSource as FileSource;

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str, bytes: &[u8]) -> (std::path::PathBuf, String) {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-file-source-{name}-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ));
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, bytes).unwrap();
        let file_name = path.to_string_lossy().into_owned();
        (path, file_name)
    }

    #[test]
    fn c_path_open_retains_full_owner_and_perror_clips_only_its_label() {
        let (path, _) = fixture("c-prefix", b"abc");
        let mut full = path.clone();
        full.as_mut_os_string().push("\0hidden-suffix");
        let mut source = FileSource::new_path(&full).unwrap();
        assert_eq!(source.file_path(), full.as_path());
        let mut bytes = [0; 3];
        source.read(&mut bytes).unwrap();
        assert_eq!(&bytes, b"abc");
        source.close().unwrap();
        std::fs::remove_file(&path).unwrap();
        let error = FileSource::new_path(&full).unwrap_err();
        let saved = match error {
            IoError::FileOpenDetails(saved) => saved,
            other => panic!("{other}"),
        };
        assert_eq!(saved.path, full);
        let mut label = b"\nError opening file ".to_vec();
        label.extend_from_slice(&crate::util::io::path_diagnostic_bytes(&path).unwrap());
        label.extend_from_slice(b": ");
        label.extend_from_slice(&saved.reason);
        label.push(b'\n');
        assert_eq!(saved.try_perror_bytes().unwrap(), label);
        assert!(saved
            .try_message_bytes()
            .unwrap()
            .ends_with(b"\0hidden-suffix"));
    }

    #[test]
    fn reads_requested_count_and_sets_eof_only_after_eof_is_observed() {
        let (path, name) = fixture("eof", b"abc");
        let mut source = FileSource::new(&name).unwrap();
        assert!(source.seekable());

        let mut exact = [0; 3];
        assert_eq!(source.read(&mut exact).unwrap(), 3);
        assert_eq!(&exact, b"abc");
        assert!(!source.eof().unwrap());

        let mut tail = [0; 2];
        assert_eq!(source.read(&mut tail).unwrap(), 0);
        assert!(source.eof().unwrap());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn partial_read_marks_eof_and_empty_read_preserves_it() {
        let (path, name) = fixture("partial", b"xy");
        let mut source = FileSource::new(&name).unwrap();
        let mut bytes = [0; 4];
        assert_eq!(source.read(&mut bytes).unwrap(), 2);
        assert_eq!(&bytes[..2], b"xy");
        assert!(source.eof().unwrap());
        assert_eq!(source.read(&mut []).unwrap(), 0);
        assert!(source.eof().unwrap());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn seek_rewind_and_file_size_match_stdio_state_changes() {
        let (path, name) = fixture("seek", b"abcdef");
        let mut source = FileSource::new(&name).unwrap();
        source.seek(2, SeekFrom::Start(0)).unwrap();
        assert_eq!(source.tell().unwrap(), 2);
        assert_eq!(source.file_size().unwrap(), 6);
        assert_eq!(source.tell().unwrap(), 2);
        source.seek(0, SeekFrom::End(0)).unwrap();
        let mut byte = [0];
        assert_eq!(source.read(&mut byte).unwrap(), 0);
        assert!(source.eof().unwrap());
        source.rewind().unwrap();
        assert_eq!(source.tell().unwrap(), 0);
        assert!(!source.eof().unwrap());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn supplied_file_is_non_seekable_by_contract_and_close_invalidates_handle() {
        let (path, name) = fixture("owned", b"data");
        let file = StdFile::open(&path).unwrap();
        let mut source = FileSource::from_file(&name, file);
        assert!(!source.seekable());
        assert_eq!(source.file_name(), name);
        source.close().unwrap();
        assert_eq!(
            source.read(&mut [0]).unwrap_err(),
            IoError::path_message(b"Error reading file ", &name)
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn raw_file_access_shares_position() {
        let (path, name) = fixture("raw", b"12");
        let mut source = FileSource::new(&name).unwrap();
        source.file().unwrap().seek(SeekFrom::Start(1)).unwrap();
        assert_eq!(source.tell().unwrap(), 1);
        std::fs::remove_file(path).unwrap();
    }

    #[cfg(any(unix, windows))]
    #[test]
    fn shared_output_input_uses_one_stdio_cursor_eof_and_close_state() {
        let (path, name) = fixture("shared-stdio", b"");
        let mut output = crate::util::io::FileSink::new(&name, "w+b", false, 0).unwrap();
        output.write(b"abcdef").unwrap();
        let shared = output.shared_stdio().unwrap().unwrap();
        output.rewind().unwrap();
        let mut input = FileSource::from_shared(&name, shared);
        let mut byte = [0];
        assert_eq!(input.read(&mut byte).unwrap(), 1);
        assert_eq!(byte, [b'a']);
        // fread can prefetch beyond one byte into libc's input buffer. Both
        // aliases must see its logical position, not the kernel offset.
        assert_eq!(output.tell().unwrap(), 1);
        assert_eq!(input.tell().unwrap(), 1);
        assert_eq!(input.read(&mut [0; 8]).unwrap(), 5);
        assert!(input.eof().unwrap());
        output.rewind().unwrap();
        assert!(!input.eof().unwrap());
        assert_eq!(input.tell().unwrap(), 0);
        input.close().unwrap();
        // Source reusing a dangling FILE after explicit close is undefined;
        // Rust keeps a deterministic invalidated-owner boundary instead.
        assert!(output.tell().is_err());
        output.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn stdio_eof_is_sticky_until_rewind_even_if_file_grows() {
        use std::io::Write;
        let (path, name) = fixture("sticky-eof", b"a");
        let mut source = FileSource::new(&name).unwrap();
        assert_eq!(source.read(&mut [0; 2]).unwrap(), 1);
        assert!(source.eof().unwrap());
        std::fs::OpenOptions::new()
            .append(true)
            .open(&path)
            .unwrap()
            .write_all(b"b")
            .unwrap();
        assert_eq!(source.read(&mut [0]).unwrap(), 0);
        source.rewind().unwrap();
        let mut data = [0; 2];
        assert_eq!(source.read(&mut data).unwrap(), 2);
        assert_eq!(&data, b"ab");
        source.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn directory_read_preserves_actual_stdio_error() {
        let path = std::env::temp_dir();
        let name = path.to_string_lossy();
        let mut source = FileSource::new(&name).unwrap();
        let IoError::StdioDetails(error) = source.read(&mut [0]).unwrap_err() else {
            panic!("fread failure must preserve its saved errno")
        };
        assert_eq!(
            error.operation,
            crate::util::io::cstdio_file::StdioOperation::Read
        );
        assert_eq!(error.raw_os_error, Some(libc::EISDIR));
        assert_eq!(
            error.message_bytes(),
            format!("Error reading file {name}").as_bytes()
        );
        source.close().unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn nonseekable_stdio_rewind_ignores_seek_error_and_retains_read_buffer() {
        use std::io::Write;
        use std::os::fd::FromRawFd;
        let mut descriptors = [-1; 2];
        // A real pipe matches Source's stdin control and needs no socket
        // permission in restricted build/test environments.
        assert_eq!(unsafe { libc::pipe(descriptors.as_mut_ptr()) }, 0);
        let consumer = unsafe { StdFile::from_raw_fd(descriptors[0]) };
        let mut producer = unsafe { StdFile::from_raw_fd(descriptors[1]) };
        producer.write_all(b"abc").unwrap();
        drop(producer);
        let mut source = FileSource::try_from_file("pipe", consumer).unwrap();
        let mut first = [0];
        assert_eq!(source.read(&mut first).unwrap(), 1);
        assert_eq!(first, [b'a']);
        let IoError::StdioDetails(error) = source.tell().unwrap_err() else {
            panic!("nonseekable ftell must preserve saved errno")
        };
        assert_eq!(error.raw_os_error, Some(libc::ESPIPE));
        assert_eq!(
            error.message_bytes(),
            b"Error executing ftell on stream pipe"
        );
        source.rewind().unwrap();
        let mut tail = [0; 8];
        assert_eq!(source.read(&mut tail).unwrap(), 2);
        assert_eq!(&tail[..2], b"bc");
        assert!(source.eof().unwrap());
        source.close().unwrap();
    }
}
