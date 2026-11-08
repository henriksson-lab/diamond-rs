//! Translation of `diamond/src/util/io/file.{h,cpp}`.
//!
//! This Source class has its own stdio contract: reads and fopen failures are
//! ordinary exceptions, writes/seeks retain perror, and fclose is unchecked.

use super::cstdio_file::OwnedCstdio;
use super::TempFileData;
use super::{IoError, IoResult};
use std::fs::File as StdFile;
use std::io::SeekFrom;
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub struct File {
    file: Option<OwnedCstdio>,
    auto_delete: bool,
    unlinked: bool,
    file_name: PathBuf,
    // Safe instance-owned equivalent of C++ read(size_t)'s static MemBuffer.
    scratch: Vec<u8>,
}
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Temporary;

fn path_bytes(path: &Path) -> IoResult<Vec<u8>> {
    super::path_diagnostic_bytes(path)
}
fn ordinary_error(prefix: &[u8], path: &Path, reason: &[u8]) -> IoError {
    let mut message = prefix.to_vec();
    let bytes = match path_bytes(path) {
        Ok(bytes) => bytes,
        Err(error) => return error,
    };
    message.extend_from_slice(&bytes);
    message.extend_from_slice(b". ");
    message.extend_from_slice(reason);
    IoError::RawDetails(super::RawDiagnosticError(message))
}
fn open_error(error: IoError, path: &Path, temporary: bool) -> IoError {
    match error {
        IoError::FileOpenDetails(error) => ordinary_error(
            if temporary {
                b"Error opening temporary file "
            } else {
                b"Error opening file "
            },
            path,
            &error.reason,
        ),
        error => error,
    }
}
impl File {
    pub fn new(name: &str, mode: &str) -> IoResult<Self> {
        Self::open_path(Path::new(name), mode)
    }
    pub fn new_path(name: &Path, mode: &str) -> IoResult<Self> {
        Self::open_path(name, mode)
    }
    pub fn open(name: &str, mode: &str) -> IoResult<Self> {
        Self::open_path(Path::new(name), mode)
    }
    pub fn open_path(name: &Path, mode: &str) -> IoResult<Self> {
        let file =
            OwnedCstdio::open_file_role_path(name, mode).map_err(|e| open_error(e, name, false))?;
        Ok(Self {
            file: Some(file),
            auto_delete: false,
            unlinked: false,
            file_name: name.to_owned(),
            scratch: Vec::new(),
        })
    }
    pub fn new_temporary(_: Temporary) -> IoResult<Self> {
        Self::temporary()
    }
    pub fn temporary() -> IoResult<Self> {
        #[cfg(unix)]
        let mut data = TempFileData::init(true)?;
        #[cfg(not(unix))]
        let data = TempFileData::init(true)?;
        #[cfg(unix)]
        let descriptor = data.take_file()?;
        #[cfg(unix)]
        let mut file = OwnedCstdio::from_file_path(&data.name, descriptor, "w+b", 0)
            .map_err(|e| open_error(e, &data.name, true))?;
        #[cfg(not(unix))]
        let mut file = OwnedCstdio::open_file_role_path(&data.name, "w+b")
            .map_err(|e| open_error(e, &data.name, true))?;
        file.set_file_role_buffer(64 * super::MEGABYTES)?;
        Ok(Self {
            file: Some(file),
            auto_delete: true,
            unlinked: data.unlinked,
            file_name: data.name,
            scratch: Vec::new(),
        })
    }
    pub fn close(&mut self) -> IoResult<()> {
        if let Some(mut file) = self.file.take() {
            // Source File::close explicitly ignores fclose's return value.
            let _ = file.close_file_role();
            if self.auto_delete
                && !self.unlinked
                && std::fs::remove_file(super::source_c_path(&self.file_name).as_ref()).is_err()
            {
                let mut message = b"Warning: Failed to delete temporary file ".to_vec();
                message.extend_from_slice(&path_bytes(
                    super::source_c_path(&self.file_name).as_ref(),
                )?);
                message.push(b'\n');
                // Source uses fprintf(stderr), without a cout tie or errno suffix.
                super::console::write_perror(&message);
            }
        }
        Ok(())
    }
    pub fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        let file = self
            .file
            .as_mut()
            .ok_or_else(|| IoError::path_message(b"Error writing file ", &self.file_name))?;
        file.write_bytes(ptr).map_err(|error| match error {
            IoError::StdioDetails(mut error) => {
                error.message = b"Error writing to temporary file ".to_vec();
                let bytes = match path_bytes(&self.file_name) {
                    Ok(bytes) => bytes,
                    Err(conversion) => return conversion,
                };
                error.message.extend_from_slice(&bytes);
                IoError::StdioDetails(error)
            }
            error => error,
        })
    }
    pub fn write_value<T: FilePrimitive>(&mut self, value: T) -> IoResult<()> {
        self.write(&value.to_ne_bytes_vec())
    }
    pub fn read_value<T: FilePrimitive>(&mut self) -> IoResult<T> {
        let mut bytes = vec![0; T::SIZE];
        self.read_exact(&mut bytes)?;
        T::from_ne_bytes_slice(&bytes)
    }
    pub fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::Other("Error calling fseek.".into()))?
            .seek(position, origin)
    }
    pub fn tell(&mut self) -> IoResult<i64> {
        let file = self.file.as_mut().ok_or_else(|| {
            IoError::path_message(b"Error executing ftell on stream ", &self.file_name)
        })?;
        file.position().map_err(|error| {
            #[cfg(windows)]
            if let IoError::StdioDetails(mut error) = error {
                // Source WIN32 _ftelli64 throws without calling perror.
                error.perror = None;
                error.message = b"Error executing ftell on stream ".to_vec();
                let bytes = match path_bytes(&self.file_name) {
                    Ok(bytes) => bytes,
                    Err(conversion) => return conversion,
                };
                error.message.extend_from_slice(&bytes);
                return IoError::StdioDetails(error);
            }
            error
        })
    }
    pub fn size(&mut self) -> IoResult<i64> {
        let position = self.tell()?;
        self.seek(0, SeekFrom::End(0))?;
        let size = self.tell()?;
        self.seek(position, SeekFrom::Start(0))?;
        Ok(size)
    }
    /// Rust-only descriptor compatibility view. Actual Source operations retain FILE.
    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::path_message(b"Closed file ", &self.file_name))?
            .legacy_file()
    }
    pub fn read_exact(&mut self, ptr: &mut [u8]) -> IoResult<()> {
        let file = self
            .file
            .as_mut()
            .ok_or_else(|| IoError::path_message(b"Error reading file ", &self.file_name))?;
        let (count, saved_errno) = file.read_file_role(ptr)?;
        if count == ptr.len() {
            Ok(())
        } else {
            // Only Source's short exact-read branch evaluates strerror(errno).
            let reason = unsafe { libc::strerror(saved_errno) };
            let reason = if reason.is_null() {
                Vec::new()
            } else {
                unsafe { std::ffi::CStr::from_ptr(reason) }
                    .to_bytes()
                    .to_vec()
            };
            Err(ordinary_error(
                b"Error reading file ",
                &self.file_name,
                &reason,
            ))
        }
    }
    pub fn read_max(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let file = self
            .file
            .as_mut()
            .ok_or_else(|| IoError::path_message(b"Error reading file ", &self.file_name))?;
        // Source returns fread's count even when ferror is set.
        file.read_file_role(ptr).map(|(count, _)| count)
    }
    pub fn read(&mut self, count: usize) -> IoResult<&[u8]> {
        let mut bytes = std::mem::take(&mut self.scratch);
        bytes.resize(count, 0);
        let result = self.read_exact(&mut bytes);
        self.scratch = bytes;
        result?;
        Ok(&self.scratch)
    }
    /// Rust extension; Source File exposes no eof member.
    pub fn eof(&mut self) -> IoResult<bool> {
        Ok(self.tell()? >= self.size()?)
    }
    pub fn file_path(&self) -> &Path {
        &self.file_name
    }
    pub fn try_file_name(&self) -> IoResult<&str> {
        self.file_name
            .to_str()
            .ok_or_else(|| IoError::Other("File name is not UTF-8; use file_path.".into()))
    }
    /// Checked legacy String API; actual raw-path callers use file_path.
    pub fn file_name(&self) -> &str {
        self.try_file_name()
            .expect("File name is not UTF-8; use file_path or try_file_name")
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
            crate::util::test_path_component()
        ))
    }

    #[test]
    fn file_role_retains_stdio_buffer_until_seek_or_unchecked_close() {
        let path = named_path("buffered");
        let mut file = File::open_path(&path, "w+b").unwrap();
        file.write(b"x").unwrap();
        assert_eq!(file.tell().unwrap(), 1);
        assert_eq!(std::fs::metadata(&path).unwrap().len(), 0);
        assert_eq!(file.size().unwrap(), 1);
        assert_eq!(std::fs::metadata(&path).unwrap().len(), 1);
        file.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn file_role_read_max_ignores_directory_fread_error() {
        let mut file = File::open_path(&std::env::temp_dir(), "rb").unwrap();
        assert_eq!(file.read_max(&mut [0; 1]).unwrap(), 0);
        let IoError::RawDetails(error) = file.read_exact(&mut [0; 1]).unwrap_err() else {
            panic!("read is ordinary exception")
        };
        assert!(error.message_bytes().starts_with(b"Error reading file "));
        file.close().unwrap();
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
        let error = file.read_exact(&mut bytes).unwrap_err();
        let IoError::RawDetails(error) = error else {
            panic!("Source File read must be an ordinary raw exception")
        };
        let mut prefix = b"Error reading file ".to_vec();
        prefix.extend_from_slice(name.as_bytes());
        prefix.extend_from_slice(b". ");
        assert!(error.message_bytes().starts_with(&prefix));
        drop(file);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn negative_absolute_seek_and_closed_file_report_errors() {
        let path = named_path("errors");
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        let mut file = File::new(&name, "w+b").unwrap();
        let IoError::StdioDetails(error) = file.seek(-1, SeekFrom::Start(0)).unwrap_err() else {
            panic!("fseek must preserve saved perror")
        };
        assert_eq!(error.message, b"Error calling fseek.");
        assert!(error.perror.is_some());
        file.close().unwrap();
        assert_eq!(
            file.write(b"x").unwrap_err(),
            IoError::path_message(b"Error writing file ", Path::new(name.as_ref()))
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
