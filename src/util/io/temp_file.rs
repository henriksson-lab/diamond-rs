//! Translation of `diamond/src/util/io/temp_file.{h,cpp}`.

use std::fs::OpenOptions;
use std::path::Path;
#[cfg(test)]
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

#[cfg(unix)]
use std::os::fd::IntoRawFd;
#[cfg(unix)]
use std::os::unix::fs::OpenOptionsExt;

use super::{Compressor, IoError, IoResult, OutputFile};

static NEXT_TEMP_FILE: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, PartialEq, Eq)]
pub struct TempFileData {
    pub name: String,
    pub fd: i32,
    pub unlinked: bool,
}

impl TempFileData {
    pub fn init(unlink: bool) -> IoResult<Self> {
        Self::init_in(unlink, std::env::temp_dir(), false)
    }

    /// Explicit-config equivalent of C++ `TempFile::init`, replacing accesses
    /// to the global `config.tmpdir` and `config.no_unlink`.
    pub fn init_in(unlink: bool, temp_dir: impl AsRef<Path>, no_unlink: bool) -> IoResult<Self> {
        #[cfg(not(unix))]
        let _ = (unlink, no_unlink);
        let temp_dir = temp_dir.as_ref();
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
            match options.open(&path) {
                Ok(file) => {
                    let name = path.to_string_lossy().into_owned();
                    #[cfg(unix)]
                    {
                        let fd = file.into_raw_fd();
                        let unlinked = if no_unlink || !unlink {
                            false
                        } else {
                            std::fs::remove_file(&path).is_ok()
                        };
                        return Ok(Self { name, fd, unlinked });
                    }
                    #[cfg(not(unix))]
                    {
                        // Windows C++ returns a name and opens it later.
                        drop(file);
                        return Ok(Self {
                            name,
                            fd: -1,
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
        Err(IoError::Other(format!(
            "Error opening temporary file in {}",
            temp_dir.display()
        )))
    }

    #[cfg(unix)]
    pub(crate) fn take_fd(&mut self) -> i32 {
        std::mem::replace(&mut self.fd, -1)
    }
}

#[cfg(unix)]
impl Drop for TempFileData {
    fn drop(&mut self) {
        if self.fd >= 0 {
            unsafe extern "C" {
                fn close(fd: i32) -> i32;
            }
            let _ = unsafe { close(self.fd) };
            self.fd = -1;
        }
    }
}

#[derive(Debug, Default)]
pub struct TempFileHandler;

impl TempFileHandler {
    pub const fn new() -> Self {
        Self
    }

    pub fn init(&mut self, _path: &str) -> IoResult<()> {
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
        Ok(Self {
            output: OutputFile::new(file_name, Compressor::None, "wb")?,
            unlinked: false,
        })
    }

    /// C++ `TempFile(const TempFileData&)`.
    pub fn from_temp_file_data(data: &TempFileData) -> IoResult<Self> {
        Ok(Self {
            output: OutputFile::from_temp_file_data(data, Compressor::None, "w+b")?,
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

    pub fn get_temp_dir() -> IoResult<String> {
        Self::get_temp_dir_in(std::env::temp_dir())
    }

    pub fn get_temp_dir_in(temp_dir: impl AsRef<Path>) -> IoResult<String> {
        let data = Self::init_in(true, temp_dir, false)?;
        let directory = Path::new(&data.name)
            .parent()
            .map(|path| path.to_string_lossy().into_owned())
            .unwrap_or_default();
        let mut temp = Self::from_temp_file_data(&data)?;
        temp.close()?;
        if !data.unlinked {
            let _ = std::fs::remove_file(&data.name);
        }
        Ok(directory)
    }

    pub fn output_file(&mut self) -> &mut OutputFile {
        &mut self.output
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
    use std::io::{Read, Seek, SeekFrom};

    fn test_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "diamond-rs-temp-file-{name}-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
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
        assert!(Path::new(&linked.name).exists());
        std::fs::remove_file(&linked.name).unwrap();

        let unlinked = TempFile::init_in(true, &dir, false).unwrap();
        #[cfg(unix)]
        {
            assert!(unlinked.unlinked);
            assert!(!Path::new(&unlinked.name).exists());
        }
        drop(unlinked);
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn temp_file_data_constructor_keeps_unlinked_descriptor_usable() {
        let dir = test_dir("descriptor");
        let data = TempFile::init_in(true, &dir, false).unwrap();
        let mut temp = TempFile::from_temp_file_data(&data).unwrap();
        temp.write_raw(b"abc").unwrap();
        temp.flush().unwrap();
        temp.file().unwrap().seek(SeekFrom::Start(0)).unwrap();
        let mut bytes = [0; 3];
        temp.file().unwrap().read_exact(&mut bytes).unwrap();
        assert_eq!(&bytes, b"abc");
        temp.close().unwrap();
        drop(data);
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
}
