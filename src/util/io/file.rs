//! Translation of `diamond/src/util/io/file.{h,cpp}`.

use super::{IoError, IoResult};
use std::fs::{File as StdFile, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};

#[cfg(unix)]
use super::TempFileData;

#[derive(Debug)]
pub struct File {
    file: Option<StdFile>,
    auto_delete: bool,
    unlinked: bool,
    file_name: String,
    // Safe instance-owned equivalent of C++ `read(size_t)`'s static MemBuffer.
    scratch: Vec<u8>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Temporary;

impl File {
    /// C++ `File(const string&, const char*)`.
    pub fn new(name: &str, mode: &str) -> IoResult<Self> {
        Self::open(name, mode)
    }

    /// C++ overload `File(Temporary)`.
    pub fn new_temporary(_temporary: Temporary) -> IoResult<Self> {
        Self::temporary()
    }

    pub fn open(name: &str, mode: &str) -> IoResult<Self> {
        let mut options = OpenOptions::new();
        if mode.contains('r') {
            options.read(true);
        }
        if mode.contains('w') {
            options.write(true).create(true).truncate(true);
        }
        if mode.contains('a') {
            options.append(true).create(true);
        }
        if mode.contains('+') {
            options.read(true).write(true);
        }
        let file = options
            .open(name)
            .map_err(|error| IoError::Other(format!("Error opening file {name}. {error}")))?;
        Ok(Self {
            file: Some(file),
            auto_delete: false,
            unlinked: false,
            file_name: name.to_string(),
            scratch: Vec::new(),
        })
    }

    pub fn temporary() -> IoResult<Self> {
        #[cfg(unix)]
        {
            let mut data = TempFileData::init(true)?;
            let file = data.take_file()?;
            return Ok(Self {
                file: Some(file),
                auto_delete: true,
                unlinked: data.unlinked,
                file_name: data.name.clone(),
                scratch: Vec::new(),
            });
        }
        #[cfg(not(unix))]
        {
            let mut path = std::env::temp_dir();
            path.push(format!(
                "diamond-rs-{}-{}.tmp",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            let name = path.to_string_lossy().into_owned();
            let file = OpenOptions::new()
                .read(true)
                .write(true)
                .create_new(true)
                .open(&path)
                .map_err(|error| {
                    IoError::Other(format!("Error opening temporary file {name}. {error}"))
                })?;
            Ok(Self {
                file: Some(file),
                auto_delete: true,
                unlinked: false,
                file_name: name,
                scratch: Vec::new(),
            })
        }
    }

    pub fn close(&mut self) -> IoResult<()> {
        let was_open = self.file.take().is_some();
        if was_open && self.auto_delete && !self.unlinked {
            if let Err(error) = std::fs::remove_file(&self.file_name) {
                eprintln!(
                    "Warning: Failed to delete temporary file {}: {error}",
                    self.file_name
                );
            }
        }
        Ok(())
    }

    pub fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileWrite(self.file_name.clone()))?
            .write_all(ptr)
            .map_err(|_| IoError::FileWrite(self.file_name.clone()))
    }

    /// Header template `write(const T&)`.
    pub fn write_value<T: FilePrimitive>(&mut self, value: T) -> IoResult<()> {
        self.write(&value.to_ne_bytes_vec())
    }

    /// Header template `read(T&)`, returned by value for safe Rust ownership.
    pub fn read_value<T: FilePrimitive>(&mut self) -> IoResult<T> {
        let mut bytes = vec![0; T::SIZE];
        self.read_exact(&mut bytes)?;
        T::from_ne_bytes_slice(&bytes)
    }

    pub fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
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

    pub fn tell(&mut self) -> IoResult<i64> {
        self.file
            .as_mut()
            .ok_or_else(|| {
                IoError::Other(format!(
                    "Error executing ftell on stream {}",
                    self.file_name
                ))
            })?
            .stream_position()
            .map(|position| position as i64)
            .map_err(|_| IoError::Other("Error calling ftell.".to_string()))
    }

    pub fn size(&mut self) -> IoResult<i64> {
        let position = self.tell()?;
        self.seek(0, SeekFrom::End(0))?;
        let size = self.tell()?;
        self.seek(position, SeekFrom::Start(0))?;
        Ok(size)
    }

    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::Other(format!("Closed file {}", self.file_name)))
    }

    pub fn read_exact(&mut self, ptr: &mut [u8]) -> IoResult<()> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileRead(self.file_name.clone()))?
            .read_exact(ptr)
            .map_err(|_| IoError::FileRead(self.file_name.clone()))
    }

    pub fn read_max(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileRead(self.file_name.clone()))?
            .read(ptr)
            .map_err(|_| IoError::FileRead(self.file_name.clone()))
    }

    pub fn read(&mut self, count: usize) -> IoResult<&[u8]> {
        let mut bytes = std::mem::take(&mut self.scratch);
        bytes.resize(count, 0);
        let result = self.read_exact(&mut bytes);
        self.scratch = bytes;
        result?;
        Ok(&self.scratch)
    }

    /// Rust compatibility helper; equivalent to testing `tell() == size()`.
    pub fn eof(&mut self) -> IoResult<bool> {
        let file = self
            .file
            .as_mut()
            .ok_or_else(|| IoError::FileRead(self.file_name.clone()))?;
        let position = file
            .stream_position()
            .map_err(|_| IoError::Other("Error calling ftell.".to_string()))?;
        let length = file
            .metadata()
            .map_err(|_| IoError::FileRead(self.file_name.clone()))?
            .len();
        Ok(position >= length)
    }

    pub fn file_name(&self) -> &str {
        &self.file_name
    }
}

impl Drop for File {
    fn drop(&mut self) {
        let _ = self.close();
    }
}

pub trait FilePrimitive: Copy {
    const SIZE: usize;
    fn to_ne_bytes_vec(self) -> Vec<u8>;
    fn from_ne_bytes_slice(bytes: &[u8]) -> IoResult<Self>;
}

macro_rules! impl_file_primitive {
    ($($type:ty),* $(,)?) => {
        $(
            impl FilePrimitive for $type {
                const SIZE: usize = std::mem::size_of::<$type>();

                fn to_ne_bytes_vec(self) -> Vec<u8> {
                    self.to_ne_bytes().to_vec()
                }

                fn from_ne_bytes_slice(bytes: &[u8]) -> IoResult<Self> {
                    Ok(<$type>::from_ne_bytes(bytes.try_into().unwrap()))
                }
            }
        )*
    };
}

impl_file_primitive!(u8, i8, u16, i16, u32, i32, u64, i64, f32, f64);

#[cfg(test)]
mod tests {
    use super::*;

    fn named_path(test: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-file-{test}-{}-{}.tmp",
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
        ))
    }

    #[test]
    fn named_file_binary_seek_size_and_position_restore() {
        let path = named_path("binary");
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        let mut file = File::new(&name, "w+b").unwrap();
        file.write(b"abc").unwrap();
        file.write_value(0x0102u16).unwrap();
        assert_eq!(file.tell().unwrap(), 5);
        file.seek(1, SeekFrom::Start(0)).unwrap();
        assert_eq!(file.size().unwrap(), 5);
        assert_eq!(file.tell().unwrap(), 1);
        assert_eq!(file.read(2).unwrap(), b"bc");
        file.seek(3, SeekFrom::Start(0)).unwrap();
        assert_eq!(file.read_value::<u16>().unwrap(), 0x0102);
        assert!(file.eof().unwrap());
        file.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn read_max_returns_partial_count_and_exact_read_errors() {
        let path = named_path("partial");
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, b"xy").unwrap();
        let name = path.to_string_lossy();
        let mut file = File::new(&name, "rb").unwrap();
        let mut bytes = [0; 4];
        assert_eq!(file.read_max(&mut bytes).unwrap(), 2);
        assert_eq!(&bytes[..2], b"xy");
        file.seek(0, SeekFrom::Start(0)).unwrap();
        assert_eq!(
            file.read_exact(&mut bytes).unwrap_err(),
            IoError::FileRead(name.into_owned())
        );
        drop(file);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn negative_absolute_seek_and_closed_file_report_errors() {
        let path = named_path("errors");
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        let mut file = File::new(&name, "w+b").unwrap();
        assert_eq!(
            file.seek(-1, SeekFrom::Start(0)).unwrap_err(),
            IoError::Other("Error calling fseek.".to_string())
        );
        file.close().unwrap();
        assert_eq!(
            file.write(b"x").unwrap_err(),
            IoError::FileWrite(name.into_owned())
        );
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn temporary_file_owns_descriptor_and_is_move_safe() {
        let mut original = File::new_temporary(Temporary).unwrap();
        #[cfg(unix)]
        let name = original.file_name().to_string();
        original.write(b"data").unwrap();
        let mut moved = original;
        moved.seek(0, SeekFrom::Start(0)).unwrap();
        assert_eq!(moved.read(4).unwrap(), b"data");
        #[cfg(unix)]
        assert!(!std::path::Path::new(&name).exists());
        moved.close().unwrap();
    }
}
