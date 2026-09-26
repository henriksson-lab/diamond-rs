//! Translation of `diamond/src/util/io/output_file.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::{Read, SeekFrom};

use super::{
    zlib_decompress, zstd_decompress, FilePrimitive, FileSink, IoError, IoResult,
    OutputStreamBuffer, Serializer, StreamEntity, TempFileData, ZlibSink, ZstdSink,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compressor {
    None,
    Zlib,
    Zstd,
}

#[derive(Debug)]
enum OutputFileStream {
    Plain(OutputStreamBuffer<FileSink>),
    Zlib(OutputStreamBuffer<ZlibSink<FileSink>>),
    Zstd(OutputStreamBuffer<ZstdSink<FileSink>>),
}

fn make_compressor(sink: FileSink, compressor: Compressor) -> OutputFileStream {
    let buffer = OutputStreamBuffer::new(sink);
    match compressor {
        Compressor::None => OutputFileStream::Plain(buffer),
        Compressor::Zlib => OutputFileStream::Zlib(OutputStreamBuffer::new(ZlibSink::new(buffer))),
        Compressor::Zstd => OutputFileStream::Zstd(OutputStreamBuffer::new(ZstdSink::new(buffer))),
    }
}

macro_rules! dispatch_stream {
    ($self:expr, $stream:ident => $expression:expr) => {
        match $self {
            Self::Plain($stream) => $expression,
            Self::Zlib($stream) => $expression,
            Self::Zstd($stream) => $expression,
        }
    };
}

impl StreamEntity for OutputFileStream {
    fn rewind(&mut self) -> IoResult<()> {
        dispatch_stream!(self, stream => stream.rewind())
    }

    fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        dispatch_stream!(self, stream => stream.seek(position, origin))
    }

    fn tell(&mut self) -> IoResult<i64> {
        dispatch_stream!(self, stream => stream.tell())
    }

    fn close(&mut self) -> IoResult<()> {
        dispatch_stream!(self, stream => stream.close())
    }

    fn file_name(&self) -> &str {
        dispatch_stream!(self, stream => stream.file_name())
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        dispatch_stream!(self, stream => stream.file())
    }

    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        dispatch_stream!(self, stream => stream.write(ptr))
    }

    fn flush(&mut self) -> IoResult<()> {
        dispatch_stream!(self, stream => stream.flush())
    }

    fn file_size(&mut self) -> IoResult<i64> {
        dispatch_stream!(self, stream => stream.file_size())
    }

    fn seekable(&self) -> bool {
        dispatch_stream!(self, stream => stream.seekable())
    }
}

#[derive(Debug)]
pub struct OutputFile {
    serializer: Serializer<OutputFileStream>,
    file_name: String,
}

impl OutputFile {
    /// C++ `OutputFile(const string&, Compressor, const char*)`.
    pub fn new(file_name: &str, compressor: Compressor, mode: &str) -> IoResult<Self> {
        let sink = FileSink::new(file_name, mode, false, 0)?;
        Ok(Self {
            serializer: Serializer::new(make_compressor(sink, compressor)),
            file_name: file_name.to_string(),
        })
    }

    /// C++ `OutputFile(const TempFileData&, Compressor, const char*)`.
    pub fn from_temp_file_data(
        data: &TempFileData,
        compressor: Compressor,
        mode: &str,
    ) -> IoResult<Self> {
        #[cfg(unix)]
        {
            let fd = unsafe { super::dup(data.fd) };
            if fd < 0 {
                return Err(IoError::Other(format!(
                    "Error opening temporary file {}",
                    data.name
                )));
            }
            let sink = unsafe { FileSink::from_fd(&data.name, fd, mode, false, 0)? };
            return Ok(Self {
                serializer: Serializer::new(make_compressor(sink, compressor)),
                file_name: data.name.clone(),
            });
        }
        #[cfg(not(unix))]
        {
            Self::new(&data.name, compressor, mode)
        }
    }

    pub fn remove(&mut self) {
        if std::fs::remove_file(&self.file_name).is_err() {
            eprintln!("Warning: Failed to delete file {}", self.file_name);
        }
    }

    pub fn advise_need(&mut self) -> IoResult<()> {
        #[cfg(all(unix, not(target_os = "macos")))]
        {
            use std::os::fd::AsRawFd;

            unsafe extern "C" {
                fn posix_fadvise(fd: i32, offset: i64, len: i64, advice: i32) -> i32;
            }
            let size = self.serializer.file_size()?;
            let fd = self.serializer.file()?.as_raw_fd();
            // Preserve the upstream bitwise expression. On POSIX this resolves
            // to the WILLNEED advice value.
            const POSIX_FADV_SEQUENTIAL: i32 = 2;
            const POSIX_FADV_WILLNEED: i32 = 3;
            let _ =
                unsafe { posix_fadvise(fd, 0, size, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED) };
        }
        Ok(())
    }

    pub fn file_name(&self) -> String {
        self.file_name.clone()
    }

    pub fn write_raw(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.serializer.write_raw(ptr)
    }

    pub fn write_i32(&mut self, value: i32) -> IoResult<&mut Self> {
        self.serializer.write_i32(value)?;
        Ok(self)
    }

    pub fn write_i16(&mut self, value: i16) -> IoResult<&mut Self> {
        self.serializer.write_i16(value)?;
        Ok(self)
    }

    pub fn write_u16(&mut self, value: u16) -> IoResult<&mut Self> {
        self.serializer.write_u16(value)?;
        Ok(self)
    }

    pub fn write_i64(&mut self, value: i64) -> IoResult<&mut Self> {
        self.serializer.write_i64(value)?;
        Ok(self)
    }

    pub fn write_u32(&mut self, value: u32) -> IoResult<&mut Self> {
        self.serializer.write_u32(value)?;
        Ok(self)
    }

    pub fn write_u64(&mut self, value: u64) -> IoResult<&mut Self> {
        self.serializer.write_u64(value)?;
        Ok(self)
    }

    pub fn write_f64(&mut self, value: f64) -> IoResult<&mut Self> {
        self.serializer.write_f64(value)?;
        Ok(self)
    }

    pub fn write_value<T: FilePrimitive>(&mut self, value: T) -> IoResult<&mut Self> {
        self.serializer.write_value(value)?;
        Ok(self)
    }

    pub fn write_values<T: FilePrimitive>(&mut self, values: &[T]) -> IoResult<&mut Self> {
        self.serializer.write_values(values)?;
        Ok(self)
    }

    pub fn write_string(&mut self, value: &str) -> IoResult<&mut Self> {
        self.serializer.write_string(value)?;
        Ok(self)
    }

    pub fn file_size(&mut self) -> IoResult<i64> {
        self.serializer.file_size()
    }

    pub fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        self.serializer.seek(position, origin)
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.serializer.rewind()
    }

    pub fn tell(&mut self) -> IoResult<usize> {
        self.serializer.tell()
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.serializer.close()
    }

    pub fn flush(&mut self) -> IoResult<()> {
        self.serializer.flush()
    }

    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.serializer.file()
    }
}

/// Safe slice-based compatibility form of C++ `decompress(FILE*, ...)`.
pub fn decompress(src: &[u8], dst: &mut [u8], compressor: Compressor) -> IoResult<usize> {
    match compressor {
        Compressor::Zlib => zlib_decompress(src, dst),
        Compressor::Zstd => zstd_decompress(src, dst),
        Compressor::None => Err(IoError::Other(
            "Invalid compressor in decompress".to_string(),
        )),
    }
}

/// Direct mapping of C++ `decompress(FILE*, ...)`, reading from the current
/// file position through EOF.
pub fn decompress_file(
    src: &mut StdFile,
    dst: &mut [u8],
    compressor: Compressor,
) -> IoResult<usize> {
    let mut compressed = Vec::new();
    src.read_to_end(&mut compressed)
        .map_err(|error| IoError::Other(format!("Error reading file: {error}")))?;
    decompress(&compressed, dst, compressor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Seek, SeekFrom};

    fn path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-output-file-{name}-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
        ))
    }

    #[test]
    fn plain_output_preserves_exact_serializer_bytes_and_position() {
        let path = path("plain");
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        let mut output = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
        output.write_u32(0x0102_0304).unwrap();
        output.write_string("x").unwrap();
        assert_eq!(output.tell().unwrap(), 6);
        assert_eq!(output.file_size().unwrap(), 6);
        output.close().unwrap();
        assert_eq!(
            std::fs::read(&path).unwrap(),
            [0x04, 0x03, 0x02, 0x01, b'x', 0]
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn compressed_output_works_with_upstream_default_write_only_mode() {
        for (compressor, expected) in [
            (Compressor::Zlib, Compressor::Zlib),
            (Compressor::Zstd, Compressor::Zstd),
        ] {
            let path = path(match compressor {
                Compressor::Zlib => "gzip",
                Compressor::Zstd => "zstd",
                Compressor::None => unreachable!(),
            });
            let _ = std::fs::remove_file(&path);
            let name = path.to_string_lossy();
            let mut output = OutputFile::new(&name, compressor, "wb").unwrap();
            output.write_raw(b"compressed-output").unwrap();
            output.close().unwrap();
            let bytes = std::fs::read(&path).unwrap();
            assert_eq!(super::super::detect_compressor(&bytes), expected);
            let mut decoded = [0; 17];
            assert_eq!(decompress(&bytes, &mut decoded, compressor).unwrap(), 17);
            assert_eq!(&decoded, b"compressed-output");
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn temp_data_constructor_and_file_decompress_preserve_bytes() {
        let data = TempFileData::init(false).unwrap();
        {
            let mut output =
                OutputFile::from_temp_file_data(&data, Compressor::Zlib, "w+b").unwrap();
            output.write_raw(b"temporary").unwrap();
            output.close().unwrap();
        }
        let mut file = StdFile::open(&data.name).unwrap();
        let mut decoded = [0; 9];
        assert_eq!(
            decompress_file(&mut file, &mut decoded, Compressor::Zlib).unwrap(),
            9
        );
        assert_eq!(&decoded, b"temporary");
        let name = data.name.clone();
        drop(data);
        std::fs::remove_file(name).unwrap();
    }

    #[test]
    fn remove_and_advise_need_map_header_helpers() {
        let path = path("helpers");
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        let mut output = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
        output.write_raw(b"data").unwrap();
        output.flush().unwrap();
        output.file().unwrap().seek(SeekFrom::Start(0)).unwrap();
        output.advise_need().unwrap();
        output.close().unwrap();
        output.remove();
        assert!(!path.exists());
    }

    #[test]
    fn invalid_decompressor_is_rejected() {
        assert_eq!(
            decompress(b"data", &mut [0; 4], Compressor::None).unwrap_err(),
            IoError::Other("Invalid compressor in decompress".to_string())
        );
    }
}
