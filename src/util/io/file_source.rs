//! Translation of `diamond/src/util/io/file_source.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::{Read, Seek, SeekFrom};

use super::{IoError, IoResult, StreamEntity};

#[derive(Debug)]
pub struct FileSource {
    file: Option<StdFile>,
    file_name: String,
    seekable: bool,
    eof: bool,
}

impl FileSource {
    /// C++ `FileSource(const string&)`.
    pub fn new(file_name: &str) -> IoResult<Self> {
        let is_stdin = file_name.is_empty() || file_name == "-";
        let stdin_name = if cfg!(windows) {
            "CONIN$"
        } else {
            "/dev/stdin"
        };
        let open_name = if is_stdin { stdin_name } else { file_name };

        let file =
            StdFile::open(open_name).map_err(|_| IoError::FileOpen(file_name.to_string()))?;
        let seekable = !is_stdin && file.metadata().map(|meta| meta.is_file()).unwrap_or(false);
        Ok(Self {
            file: Some(file),
            file_name: file_name.to_string(),
            seekable,
            eof: false,
        })
    }

    /// C++ `FileSource(const string&, FILE*)`.
    ///
    /// Rust takes ownership of a cloned or otherwise owned file handle rather
    /// than borrowing a raw `FILE*`.
    pub fn from_file(file_name: &str, file: StdFile) -> Self {
        Self {
            file: Some(file),
            file_name: file_name.to_string(),
            seekable: false,
            eof: false,
        }
    }

    fn open_file(&mut self) -> IoResult<&mut StdFile> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileRead(self.file_name.clone()))
    }
}

impl StreamEntity for FileSource {
    fn rewind(&mut self) -> IoResult<()> {
        let name = self.file_name.clone();
        self.open_file()?
            .seek(SeekFrom::Start(0))
            .map_err(|_| IoError::Other(format!("Error executing seek on file {name}")))?;
        self.eof = false;
        Ok(())
    }

    fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        let target = match origin {
            SeekFrom::Start(_) if position < 0 => {
                return Err(IoError::Other(format!(
                    "Error executing seek on file {}",
                    self.file_name
                )));
            }
            SeekFrom::Start(_) => SeekFrom::Start(position as u64),
            SeekFrom::End(_) => SeekFrom::End(position),
            SeekFrom::Current(_) => SeekFrom::Current(position),
        };
        let name = self.file_name.clone();
        self.open_file()?
            .seek(target)
            .map_err(|_| IoError::Other(format!("Error executing seek on file {name}")))?;
        self.eof = false;
        Ok(())
    }

    fn tell(&mut self) -> IoResult<i64> {
        let name = self.file_name.clone();
        self.open_file()?
            .stream_position()
            .map(|position| position as i64)
            .map_err(|_| IoError::Other(format!("Error executing ftell on stream {name}")))
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
                Err(_) => return Err(IoError::FileRead(name)),
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

    fn file_name(&self) -> &str {
        &self.file_name
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.open_file()
    }

    fn seekable(&self) -> bool {
        self.seekable
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str, bytes: &[u8]) -> (std::path::PathBuf, String) {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-file-source-{name}-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
        ));
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, bytes).unwrap();
        let file_name = path.to_string_lossy().into_owned();
        (path, file_name)
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
        assert_eq!(source.read(&mut [0]).unwrap_err(), IoError::FileRead(name));
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
}
