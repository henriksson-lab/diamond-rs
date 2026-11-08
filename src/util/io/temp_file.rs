//! Translation of `diamond/src/util/io/temp_file.{h,cpp}`.

#[cfg(unix)]
use std::fs::File as StdFile;
#[cfg(all(not(unix), not(windows)))]
use std::fs::OpenOptions;
use std::path::Path;
use std::path::PathBuf;
#[cfg(windows)]
use std::sync::atomic::AtomicU32;
#[cfg(all(not(unix), not(windows)))]
use std::sync::atomic::AtomicU64;
#[cfg(not(unix))]
use std::sync::atomic::Ordering;
use std::sync::Mutex;

use super::{Compressor, IoError, IoResult, OutputFile};

/// Saved diagnostic for Source mkstemp failure. Source prints perror(0), then
/// throws a standard runtime_error; this is distinct from FileOpenException.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TempFileOpenError {
    pub path: PathBuf,
    pub reason: Vec<u8>,
    pub raw_os_error: Option<i32>,
}
impl TempFileOpenError {
    pub fn try_message_bytes(&self) -> IoResult<Vec<u8>> {
        let mut bytes = b"Error opening temporary file ".to_vec();
        bytes.extend_from_slice(&super::path_diagnostic_bytes(&self.path)?);
        Ok(bytes)
    }

    fn new(path: &Path, error: std::io::Error) -> Self {
        let saved = super::FileOpenError::new(path, error);
        Self {
            path: saved.path,
            reason: saved.reason,
            raw_os_error: saved.raw_os_error,
        }
    }
}
impl std::fmt::Display for TempFileOpenError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Error opening temporary file {}", self.path.display())
    }
}
impl std::error::Error for TempFileOpenError {}
pub fn find_temp_open<'a>(
    error: &'a (dyn std::error::Error + 'static),
) -> Option<&'a TempFileOpenError> {
    let mut current = Some(error);
    while let Some(error) = current {
        if let Some(temp) = error.downcast_ref::<TempFileOpenError>() {
            return Some(temp);
        }
        if let Some(super::IoError::TemporaryFileOpenDetails(temp)) =
            error.downcast_ref::<super::IoError>()
        {
            return Some(temp);
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

/// Source Config::tmpdir: explicit path wins; otherwise extract the output
/// directory with the platform Source separator. Empty means current directory,
/// not the environment's system temporary directory.
pub fn source_temp_directory(explicit: &Path, output: &str) -> PathBuf {
    source_temp_directory_path(explicit, Path::new(output))
}
pub fn source_temp_directory_path(explicit: &Path, output: &Path) -> PathBuf {
    if !explicit.as_os_str().is_empty() {
        return explicit.to_path_buf();
    }
    #[cfg(unix)]
    {
        use std::os::unix::ffi::{OsStrExt, OsStringExt};
        let bytes = output.as_os_str().as_bytes();
        return bytes
            .iter()
            .rposition(|&byte| byte == b'/')
            .map(|i| PathBuf::from(std::ffi::OsString::from_vec(bytes[..i].to_vec())))
            .unwrap_or_default();
    }
    #[cfg(windows)]
    {
        use std::os::windows::ffi::{OsStrExt, OsStringExt};
        let wide: Vec<u16> = output.as_os_str().encode_wide().collect();
        return wide
            .iter()
            .rposition(|&unit| unit == u16::from(b'\\'))
            .map(|i| PathBuf::from(std::ffi::OsString::from_wide(&wide[..i])))
            .unwrap_or_default();
    }
    #[cfg(not(any(unix, windows)))]
    {
        output
            .parent()
            .unwrap_or_else(|| Path::new(""))
            .to_path_buf()
    }
}

#[cfg(all(not(unix), not(windows)))]
static NEXT_TEMP_FILE: AtomicU64 = AtomicU64::new(0);

#[cfg(windows)]
static NEXT_TEMP_FILE: AtomicU32 = AtomicU32::new(0);
#[cfg(windows)]
static TEMP_FILE_HASH: std::sync::OnceLock<u64> = std::sync::OnceLock::new();

#[derive(Debug)]
pub struct TempFileData {
    pub name: PathBuf,
    pub unlinked: bool,
    #[cfg(unix)]
    file: Option<StdFile>,
}

impl TempFileData {
    pub fn init(unlink: bool) -> IoResult<Self> {
        Self::init_in(unlink, "", false)
    }

    /// Explicit-config equivalent of C++ `TempFile::init`, replacing accesses
    /// to the global `config.tmpdir` and `config.no_unlink`.
    pub fn init_in(unlink: bool, temp_dir: impl AsRef<Path>, no_unlink: bool) -> IoResult<Self> {
        #[cfg(not(unix))]
        let _ = (unlink, no_unlink);
        let temp_dir = temp_dir.as_ref();
        #[cfg(unix)]
        {
            use std::os::fd::FromRawFd;
            use std::os::unix::ffi::{OsStrExt, OsStringExt};
            let path = temp_dir.join("diamond-tmp-XXXXXX");
            let template =
                std::ffi::CString::new(super::source_c_path(&path).as_os_str().as_bytes())
                    .map_err(|_| {
                        IoError::TemporaryFileOpenDetails(TempFileOpenError::new(
                            &path,
                            std::io::Error::from_raw_os_error(libc::EINVAL),
                        ))
                    })?;
            let mut bytes = template.into_bytes_with_nul();
            let descriptor = unsafe { libc::mkstemp(bytes.as_mut_ptr().cast()) };
            if descriptor < 0 {
                let error = std::io::Error::last_os_error();
                let name = PathBuf::from(std::ffi::OsString::from_vec(
                    std::ffi::CStr::from_bytes_with_nul(&bytes)
                        .unwrap()
                        .to_bytes()
                        .to_vec(),
                ));
                return Err(IoError::TemporaryFileOpenDetails(TempFileOpenError::new(
                    &name, error,
                )));
            }
            // SAFETY: successful mkstemp returns a new exclusively-owned FD.
            let file = unsafe { StdFile::from_raw_fd(descriptor) };
            let name = PathBuf::from(std::ffi::OsString::from_vec(
                std::ffi::CStr::from_bytes_with_nul(&bytes)
                    .unwrap()
                    .to_bytes()
                    .to_vec(),
            ));
            let unlinked = if no_unlink || !unlink {
                false
            } else {
                std::fs::remove_file(super::source_c_path(&name).as_ref()).is_ok()
            };
            return Ok(Self {
                name,
                unlinked,
                file: Some(file),
            });
        }
        #[cfg(windows)]
        {
            // Source init only generates a name on Windows. The later
            // OutputFile/FileSink constructor performs the actual creating
            // fopen and owns any FileOpenException at that boundary.
            let hash = TEMP_FILE_HASH.get_or_init(|| {
                let mut count = 0i64;
                // SAFETY: QueryPerformanceCounter writes one initialized i64.
                unsafe {
                    windows_sys::Win32::System::Performance::QueryPerformanceCounter(&mut count);
                }
                let high = (count >> 32) as u32;
                let low = count as u32;
                // LARGE_INTEGER HighPart + LowPart first uses unsigned32
                // arithmetic, then promotes for addition of QuadPart and PID.
                u64::from(high.wrapping_add(low))
                    .wrapping_add(count as u64)
                    .wrapping_add(u64::from(std::process::id()))
            });
            let serial = NEXT_TEMP_FILE.fetch_add(1, Ordering::Relaxed);
            return Ok(Self {
                name: temp_dir.join(format!("diamond-{hash:x}-{serial}.tmp")),
                unlinked: false,
            });
        }
        #[cfg(all(not(unix), not(windows)))]
        for _ in 0..128 {
            let serial = NEXT_TEMP_FILE.fetch_add(1, Ordering::Relaxed);
            let path = temp_dir.join(format!(
                "diamond-tmp-{:x}-{:x}.tmp",
                std::process::id(),
                serial
            ));
            let mut options = OpenOptions::new();
            options.read(true).write(true).create_new(true);
            #[cfg(unix)]
            options.mode(0o600);
            match options.open(super::source_c_path(&path).as_ref()) {
                Ok(file) => {
                    let name = path.to_path_buf();
                    #[cfg(unix)]
                    {
                        let unlinked = if no_unlink || !unlink {
                            false
                        } else {
                            std::fs::remove_file(super::source_c_path(&path).as_ref()).is_ok()
                        };
                        return Ok(Self {
                            name,
                            unlinked,
                            file: Some(file),
                        });
                    }
                    #[cfg(not(unix))]
                    {
                        // Fallback reservation for platforms without either
                        // the Source POSIX or Windows temporary implementation.
                        drop(file);
                        return Ok(Self {
                            name,
                            unlinked: false,
                        });
                    }
                }
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => {
                    return Err(IoError::Other(format!(
                        "Error opening temporary file {}. {error}",
                        path.display()
                    )));
                }
            }
        }
        #[cfg(all(not(unix), not(windows)))]
        {
            Err(IoError::Other(format!(
                "Error opening temporary file in {}",
                temp_dir.display()
            )))
        }
    }

    #[cfg(unix)]
    pub(crate) fn try_clone_file(&self) -> IoResult<StdFile> {
        self.file
            .as_ref()
            .ok_or_else(|| IoError::path_message(b"Error opening temporary file ", &self.name))?
            .try_clone()
            .map_err(|_| IoError::path_message(b"Error opening temporary file ", &self.name))
    }

    #[cfg(unix)]
    pub(crate) fn take_file(&mut self) -> IoResult<StdFile> {
        self.file
            .take()
            .ok_or_else(|| IoError::path_message(b"Error opening temporary file ", &self.name))
    }
}

#[derive(Debug, Default)]
pub struct TempFileHandler;

impl TempFileHandler {
    pub const fn new() -> Self {
        Self
    }

    pub fn init(&mut self, path: &str) -> IoResult<()> {
        self.init_path(Path::new(path))
    }
    pub fn init_path(&mut self, _path: &Path) -> IoResult<()> {
        // The C++ `path_` is a default-constructed const string and `init`
        // never assigns it, so repeated calls remain no-ops on Unix.
        Ok(())
    }
}

pub static TEMP_FILE_HANDLER: Mutex<TempFileHandler> = Mutex::new(TempFileHandler::new());

#[derive(Debug)]
pub struct TempFile {
    output: OutputFile,
    unlinked: bool,
}

impl TempFile {
    /// C++ `TempFile(bool unlink = true)` with the argument explicit in Rust.
    pub fn new(unlink: bool) -> IoResult<Self> {
        let data = Self::init(unlink)?;
        Self::from_temp_file_data(&data)
    }

    pub fn temporary() -> IoResult<Self> {
        Self::new(true)
    }

    /// C++ `TempFile(const string&)`.
    pub fn from_file_name(file_name: &str) -> IoResult<Self> {
        Self::from_file_path(Path::new(file_name))
    }
    pub fn from_file_path(file_name: &Path) -> IoResult<Self> {
        Ok(Self {
            output: OutputFile::new_path(file_name, Compressor::None, "wb")?,
            unlinked: false,
        })
    }

    /// C++ `TempFile(const TempFileData&)`.
    pub fn from_temp_file_data(data: &TempFileData) -> IoResult<Self> {
        Self::from_temp_file_data_with_buffer_size(data, super::DEFAULT_FILE_BUFFER_SIZE)
    }
    pub fn from_temp_file_data_with_buffer_size(
        data: &TempFileData,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Ok(Self {
            output: OutputFile::from_temp_file_data_with_buffer_size(
                data,
                Compressor::None,
                "w+b",
                buffer_size,
            )?,
            unlinked: data.unlinked,
        })
    }

    pub fn init(unlink: bool) -> IoResult<TempFileData> {
        TempFileData::init(unlink)
    }

    pub fn init_in(
        unlink: bool,
        temp_dir: impl AsRef<Path>,
        no_unlink: bool,
    ) -> IoResult<TempFileData> {
        TempFileData::init_in(unlink, temp_dir, no_unlink)
    }

    pub fn finalize(&mut self) -> IoResult<()> {
        Ok(())
    }

    pub fn get_temp_dir_path() -> IoResult<PathBuf> {
        Self::get_temp_dir_path_in_with_no_unlink(Path::new(""), false)
    }
    pub fn get_temp_dir_path_in(temp_dir: impl AsRef<Path>) -> IoResult<PathBuf> {
        Self::get_temp_dir_path_in_with_no_unlink(temp_dir, false)
    }
    pub fn get_temp_dir() -> IoResult<String> {
        Self::get_temp_dir_in("")
    }

    pub fn get_temp_dir_in(temp_dir: impl AsRef<Path>) -> IoResult<String> {
        Self::get_temp_dir_in_with_no_unlink(temp_dir, false)
    }

    pub fn get_temp_dir_in_with_no_unlink(
        temp_dir: impl AsRef<Path>,
        no_unlink: bool,
    ) -> IoResult<String> {
        Self::get_temp_dir_path_in_with_no_unlink(temp_dir, no_unlink)?
            .into_os_string()
            .into_string()
            .map_err(|_| {
                IoError::Other(
                    "OS path is not UTF-8; use get_temp_dir_path_in_with_no_unlink".into(),
                )
            })
    }
    pub fn get_temp_dir_path_in_with_no_unlink(
        temp_dir: impl AsRef<Path>,
        no_unlink: bool,
    ) -> IoResult<PathBuf> {
        let data = Self::init_in(true, temp_dir, no_unlink)?;
        let directory = data
            .name
            .parent()
            .unwrap_or_else(|| Path::new(""))
            .to_path_buf();
        let mut temp = Self::from_temp_file_data(&data)?;
        temp.close()?;
        if !data.unlinked {
            let _ = std::fs::remove_file(super::source_c_path(&data.name).as_ref());
        }
        Ok(directory)
    }

    pub fn output_file(&mut self) -> &mut OutputFile {
        &mut self.output
    }

    pub fn file_path(&self) -> &Path {
        self.output.file_path()
    }
    pub fn file_name(&self) -> String {
        self.output.file_name()
    }

    pub fn unlinked(&self) -> bool {
        self.unlinked
    }
}

impl std::ops::Deref for TempFile {
    type Target = OutputFile;

    fn deref(&self) -> &Self::Target {
        &self.output
    }
}

impl std::ops::DerefMut for TempFile {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.output
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(unix)]
    use std::io::Write;
    use std::io::{Read, Seek, SeekFrom};

    fn test_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "diamond-rs-temp-file-{name}-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir(&dir).unwrap();
        dir
    }

    #[test]
    fn init_in_honors_directory_unlink_and_no_unlink() {
        let dir = test_dir("init");
        let linked = TempFile::init_in(true, &dir, true).unwrap();
        assert_eq!(Path::new(&linked.name).parent(), Some(dir.as_path()));
        assert!(!linked.unlinked);
        #[cfg(not(windows))]
        {
            assert!(Path::new(&linked.name).exists());
            std::fs::remove_file(&linked.name).unwrap();
        }
        #[cfg(windows)]
        {
            assert!(!Path::new(&linked.name).exists());
            let mut opened = TempFile::from_temp_file_data(&linked).unwrap();
            assert!(Path::new(&linked.name).exists());
            opened.close().unwrap();
            std::fs::remove_file(&linked.name).unwrap();
        }

        let unlinked = TempFile::init_in(true, &dir, false).unwrap();
        #[cfg(all(not(unix), not(windows)))]
        let unlinked_name = unlinked.name.clone();
        #[cfg(unix)]
        {
            assert!(unlinked.unlinked);
            assert!(!Path::new(&unlinked.name).exists());
        }
        drop(unlinked);
        #[cfg(all(not(unix), not(windows)))]
        std::fs::remove_file(unlinked_name).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn temp_file_data_constructor_keeps_unlinked_descriptor_usable() {
        let dir = test_dir("descriptor");
        let data = TempFile::init_in(true, &dir, false).unwrap();
        #[cfg(not(unix))]
        let name = data.name.clone();
        let mut temp = TempFile::from_temp_file_data(&data).unwrap();
        drop(data);
        temp.write_raw(b"abc").unwrap();
        temp.flush().unwrap();
        temp.file().unwrap().seek(SeekFrom::Start(0)).unwrap();
        let mut bytes = [0; 3];
        temp.file().unwrap().read_exact(&mut bytes).unwrap();
        assert_eq!(&bytes, b"abc");
        temp.close().unwrap();
        #[cfg(not(unix))]
        std::fs::remove_file(name).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn cloned_temp_files_share_cursor_and_preserve_closed_error_text() {
        let dir = test_dir("clone-cursor");
        let mut data = TempFile::init_in(true, &dir, false).unwrap();
        let mut first = data.try_clone_file().unwrap();
        let mut second = data.try_clone_file().unwrap();

        first.write_all(b"abcdef").unwrap();
        assert_eq!(second.stream_position().unwrap(), 6);
        second.seek(SeekFrom::Start(2)).unwrap();
        assert_eq!(first.stream_position().unwrap(), 2);

        drop(first);
        drop(second);
        drop(data.take_file().unwrap());
        assert_eq!(
            data.try_clone_file().unwrap_err().to_string(),
            format!("Error opening temporary file {}", data.name.display())
        );
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn named_constructor_uses_write_only_truncating_mode() {
        let dir = test_dir("named");
        let path = dir.join("named.tmp");
        std::fs::write(&path, b"old").unwrap();
        let name = path.to_string_lossy();
        let mut temp = TempFile::from_file_name(&name).unwrap();
        assert!(!temp.unlinked());
        temp.write_raw(b"new").unwrap();
        temp.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"new");
        std::fs::remove_file(path).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn finalize_is_noop_and_get_temp_dir_uses_requested_directory() {
        let dir = test_dir("directory");
        let mut temp = TempFile::new(false).unwrap();
        temp.write_raw(b"pending").unwrap();
        let position = temp.tell().unwrap();
        temp.finalize().unwrap();
        assert_eq!(temp.tell().unwrap(), position);
        temp.close().unwrap();
        std::fs::remove_file(temp.file_name()).unwrap();

        assert_eq!(
            TempFile::get_temp_dir_in(&dir).unwrap(),
            dir.to_string_lossy()
        );
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn handler_init_is_repeatable_like_the_const_cpp_state() {
        let mut handler = TempFileHandler::new();
        handler.init("first").unwrap();
        handler.init("second").unwrap();
    }

    #[test]
    fn source_directory_prefers_explicit_then_output_parent_then_cwd() {
        let sep = if cfg!(windows) { '\\' } else { '/' };
        let output = format!("parent{sep}result.tsv");
        assert_eq!(
            source_temp_directory(Path::new("explicit"), &output),
            PathBuf::from("explicit")
        );
        assert_eq!(
            source_temp_directory(Path::new(""), &output),
            PathBuf::from("parent")
        );
        assert_eq!(
            source_temp_directory(Path::new(""), "result.tsv"),
            PathBuf::new()
        );
        assert_eq!(source_temp_directory(Path::new(""), ""), PathBuf::new());
    }

    #[cfg(windows)]
    #[test]
    fn windows_temp_creation_is_deferred_to_the_file_sink() {
        let dir = test_dir("missing-parent");
        std::fs::remove_dir(&dir).unwrap();
        let data = TempFile::init_in(true, &dir, false).unwrap();
        assert!(!Path::new(&data.name).exists());
        let error = TempFile::from_temp_file_data(&data).unwrap_err();
        assert!(super::super::file_open_error::find(&error).is_some());
        assert!(find_temp_open(&error).is_none());
        assert!(!dir.exists());
    }

    #[cfg(unix)]
    #[test]
    fn failed_mkstemp_retains_its_own_error_kind_and_path() {
        let dir = test_dir("missing-parent");
        std::fs::remove_dir(&dir).unwrap();
        let error = TempFile::init_in(true, &dir, false).unwrap_err();
        let saved = find_temp_open(&error).unwrap();
        assert_eq!(saved.raw_os_error, Some(libc::ENOENT));
        assert_eq!(Path::new(&saved.path).parent(), Some(dir.as_path()));
        // TempFile throws standard runtime_error in Source, never FileOpenException.
        assert!(super::super::file_open_error::find(&error).is_none());
        assert!(!dir.exists());
    }
}
