//! Translation of `diamond/src/util/io/file_sink.{h,cpp}`.

use std::fs::{File as StdFile, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
#[cfg(unix)]
use std::os::fd::FromRawFd;
use std::sync::Mutex;

use super::{IoError, IoResult, StreamEntity};

#[cfg(unix)]
pub fn posix_flags(mode: &str) -> IoResult<i32> {
    const O_WRONLY: i32 = 1;
    const O_RDWR: i32 = 2;
    const O_CREAT: i32 = 64;
    const O_TRUNC: i32 = 512;
    match mode {
        "wb" => Ok(O_WRONLY | O_CREAT | O_TRUNC),
        "r+b" => Ok(O_RDWR),
        "w+b" => Ok(O_RDWR | O_CREAT | O_TRUNC),
        _ => Err(IoError::Other("Invalid fopen mode.".to_string())),
    }
}

#[derive(Debug)]
pub struct FileSink {
    file: Option<StdFile>,
    file_name: String,
    async_: bool,
    is_stdout: bool,
    mtx: Mutex<()>,
}

impl FileSink {
    /// C++ `FileSink(file_name, mode, async, buffer_size)`.
    pub fn new(file_name: &str, mode: &str, async_: bool, _buffer_size: usize) -> IoResult<Self> {
        #[cfg(unix)]
        let _ = posix_flags(mode)?;

        let mut options = OpenOptions::new();
        let is_stdout = file_name.is_empty();
        if is_stdout {
            options.write(true);
            if mode.contains('+') {
                options.read(true);
            }
        } else {
            match mode {
                "wb" => {
                    options.write(true).create(true).truncate(true);
                }
                "r+b" => {
                    options.read(true).write(true);
                }
                "w+b" => {
                    options.read(true).write(true).create(true).truncate(true);
                }
                _ => return Err(IoError::Other("Invalid fopen mode.".to_string())),
            }
        }

        let stdout_name = if cfg!(windows) {
            "CONOUT$"
        } else {
            "/dev/stdout"
        };
        let open_name = if is_stdout { stdout_name } else { file_name };
        let file = options
            .open(open_name)
            .map_err(|_| IoError::FileOpen(file_name.to_string()))?;
        Ok(Self {
            file: Some(file),
            file_name: file_name.to_string(),
            async_,
            is_stdout,
            mtx: Mutex::new(()),
        })
    }

    /// Unix C++ `FileSink(file_name, fd, mode, async, buffer_size)`.
    ///
    /// # Safety
    ///
    /// `fd` must be a valid, open descriptor whose ownership is transferred to
    /// the returned sink. It must not be used elsewhere after this call.
    #[cfg(unix)]
    pub unsafe fn from_fd(
        file_name: &str,
        fd: i32,
        mode: &str,
        async_: bool,
        _buffer_size: usize,
    ) -> IoResult<Self> {
        let _ = posix_flags(mode)?;
        Ok(Self {
            file: Some(unsafe { StdFile::from_raw_fd(fd) }),
            file_name: file_name.to_string(),
            async_,
            is_stdout: false,
            mtx: Mutex::new(()),
        })
    }

    fn open_file(&mut self) -> IoResult<&mut StdFile> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileWrite(self.file_name.clone()))
    }
}

impl StreamEntity for FileSink {
    fn close(&mut self) -> IoResult<()> {
        // The C++ implementation deliberately leaves the process stdout stream
        // open. Our `/dev/stdout`/`CONOUT$` handle follows the same observable
        // behavior by remaining usable after close.
        if !self.is_stdout {
            self.file.take();
        }
        Ok(())
    }

    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        if self.async_ {
            let _guard = self
                .mtx
                .lock()
                .map_err(|_| IoError::FileWrite(self.file_name.clone()))?;
            self.file
                .as_mut()
                .ok_or_else(|| IoError::FileWrite(self.file_name.clone()))?
                .write_all(ptr)
                .map_err(|_| IoError::FileWrite(self.file_name.clone()))
        } else {
            self.open_file()?
                .write_all(ptr)
                .map_err(|_| IoError::FileWrite(self.file_name.clone()))
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
            .ok_or_else(|| IoError::Other(format!("Error executing ftell on stream {file_name}")))?
            .stream_position()
            .map(|position| position as i64)
            .map_err(|_| IoError::Other("Error calling ftell.".to_string()))
    }

    fn file_name(&self) -> &str {
        &self.file_name
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
                    .map_err(|_| IoError::FileWrite(file_name.clone()))
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

#[cfg(test)]
mod tests {
    use super::*;

    fn path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-file-sink-{name}-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
        ))
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
        assert_eq!(sink.write(b"x").unwrap_err(), IoError::FileWrite(name));
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
            sink.seek(-1, SeekFrom::Start(0)).unwrap_err(),
            IoError::Other("Error calling fseek.".to_string())
        );
        sink.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn closing_standard_output_preserves_its_handle() {
        let mut sink = FileSink::new("", "wb", false, 0).unwrap();
        sink.close().unwrap();
        assert_eq!(sink.file_name(), "");
        assert!(sink.file().is_ok());
    }

    #[cfg(unix)]
    #[test]
    fn from_fd_takes_ownership_and_writes_exact_bytes() {
        use std::os::fd::IntoRawFd;

        let path = path("fd");
        let _ = std::fs::remove_file(&path);
        let descriptor = OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap()
            .into_raw_fd();
        let name = path.to_string_lossy();
        let mut sink = unsafe { FileSink::from_fd(&name, descriptor, "w+b", false, 0) }.unwrap();
        sink.write(b"fd").unwrap();
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"fd");
        std::fs::remove_file(path).unwrap();
    }
}
