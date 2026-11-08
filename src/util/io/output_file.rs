//! Translation of `diamond/src/util/io/output_file.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::{Read, SeekFrom};
use std::path::{Path, PathBuf};

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
    make_compressor_with_buffer_size(sink, compressor, super::DEFAULT_FILE_BUFFER_SIZE)
}

fn make_compressor_with_buffer_size(
    sink: FileSink,
    compressor: Compressor,
    buffer_size: usize,
) -> OutputFileStream {
    let buffer = OutputStreamBuffer::with_source_buffer_size(sink, buffer_size);
    match compressor {
        Compressor::None => OutputFileStream::Plain(buffer),
        Compressor::Zlib => OutputFileStream::Zlib(OutputStreamBuffer::with_source_buffer_size(
            ZlibSink::new(buffer),
            buffer_size,
        )),
        Compressor::Zstd => OutputFileStream::Zstd(OutputStreamBuffer::with_source_buffer_size(
            ZstdSink::new(buffer),
            buffer_size,
        )),
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

    fn file_path(&self) -> &Path {
        dispatch_stream!(self, stream => stream.file_path())
    }
    fn file_name(&self) -> &str {
        self.file_path()
            .to_str()
            .expect("Non-UTF8 OS path; use file_path")
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        dispatch_stream!(self, stream => stream.file())
    }

    #[cfg(any(unix, windows))]
    fn shared_stdio(&mut self) -> IoResult<Option<super::cstdio_file::SharedCstdio>> {
        dispatch_stream!(self, stream => stream.shared_stdio())
    }

    fn clone_file(&mut self) -> IoResult<StdFile> {
        dispatch_stream!(self, stream => stream.clone_file())
    }

    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        dispatch_stream!(self, stream => stream.forward_write(ptr))
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
    file_name: PathBuf,
}

impl OutputFile {
    /// C++ `OutputFile(const string&, Compressor, const char*)`.
    pub fn new(file_name: &str, compressor: Compressor, mode: &str) -> IoResult<Self> {
        Self::new_path(Path::new(file_name), compressor, mode)
    }
    pub fn new_path(file_name: &Path, compressor: Compressor, mode: &str) -> IoResult<Self> {
        let sink = FileSink::new_path(file_name, mode, false, 0)?;
        Ok(Self {
            serializer: Serializer::new(make_compressor(sink, compressor)),
            file_name: file_name.to_path_buf(),
        })
    }

    /// Source construction with explicit `config.file_buffer_size`.
    pub fn with_buffer_size(
        file_name: &str,
        compressor: Compressor,
        mode: &str,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Self::with_buffer_size_path(Path::new(file_name), compressor, mode, buffer_size)
    }
    pub fn with_buffer_size_path(
        file_name: &Path,
        compressor: Compressor,
        mode: &str,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Self::with_codec_profile(file_name, compressor, mode, buffer_size, true, false)
    }
    /// Source command build profile; low-level typed constructors stay capable.
    pub fn with_source_profile(
        file_name: &str,
        compressor: Compressor,
        mode: &str,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Self::with_source_profile_path(Path::new(file_name), compressor, mode, buffer_size)
    }
    pub fn with_source_profile_path(
        file_name: &Path,
        compressor: Compressor,
        mode: &str,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Self::with_codec_profile(
            file_name,
            compressor,
            mode,
            buffer_size,
            cfg!(feature = "zstd"),
            true,
        )
    }
    pub fn source_new(file_name: &str, compressor: Compressor, mode: &str) -> IoResult<Self> {
        Self::source_new_path(Path::new(file_name), compressor, mode)
    }
    pub fn source_new_path(file_name: &Path, compressor: Compressor, mode: &str) -> IoResult<Self> {
        Self::with_source_profile_path(file_name, compressor, mode, super::DEFAULT_FILE_BUFFER_SIZE)
    }
    fn with_codec_profile(
        file_name: &Path,
        compressor: Compressor,
        mode: &str,
        buffer_size: usize,
        allow_zstd: bool,
        source_profile: bool,
    ) -> IoResult<Self> {
        if buffer_size == 0 && !source_profile {
            return Err(IoError::Other("File buffer size must be positive.".into()));
        }
        let sink = FileSink::new_path(file_name, mode, false, 0)?;
        // Source creates/truncates the sink before make_compressor rejects
        // an unavailable build feature. Keep that observable ordering.
        if compressor == Compressor::Zstd && !allow_zstd {
            return Err(IoError::Other(
                "Executable was not compiled with ZStd support.".into(),
            ));
        }
        let serializer_size = if file_name.as_os_str().is_empty() {
            OutputStreamBuffer::<FileSink>::STDOUT_BUF_SIZE
        } else {
            buffer_size
        };
        Ok(Self {
            serializer: Serializer::with_source_buffer_size(
                make_compressor_with_buffer_size(sink, compressor, buffer_size),
                serializer_size,
            ),
            file_name: file_name.to_path_buf(),
        })
    }

    /// C++ `OutputFile(const TempFileData&, Compressor, const char*)`.
    pub fn from_temp_file_data(
        data: &TempFileData,
        compressor: Compressor,
        mode: &str,
    ) -> IoResult<Self> {
        Self::from_temp_file_data_with_buffer_size(
            data,
            compressor,
            mode,
            super::DEFAULT_FILE_BUFFER_SIZE,
        )
    }
    pub fn from_temp_file_data_with_buffer_size(
        data: &TempFileData,
        compressor: Compressor,
        mode: &str,
        buffer_size: usize,
    ) -> IoResult<Self> {
        #[cfg(unix)]
        {
            let file = data.try_clone_file()?;
            let sink = FileSink::from_file_path(&data.name, file, mode, false, 0)?;
            Ok(Self {
                serializer: Serializer::with_source_buffer_size(
                    make_compressor_with_buffer_size(sink, compressor, buffer_size),
                    buffer_size,
                ),
                file_name: data.name.clone(),
            })
        }
        #[cfg(not(unix))]
        {
            Self::with_buffer_size_path(&data.name, compressor, mode, buffer_size)
        }
    }

    pub fn remove(&mut self) {
        let _ = self.try_remove();
    }
    pub fn try_remove(&mut self) -> IoResult<()> {
        if std::fs::remove_file(super::source_c_path(&self.file_name).as_ref()).is_err() {
            let mut bytes = b"Warning: Failed to delete file ".to_vec();
            bytes.extend_from_slice(&super::path_diagnostic_bytes(&self.file_name)?);
            bytes.push(b'\n');
            let _ = std::io::Write::write_all(&mut super::console::stderr(), &bytes);
        }
        Ok(())
    }

    pub fn advise_need(&mut self) -> IoResult<()> {
        #[cfg(all(unix, not(target_os = "macos")))]
        {
            use std::num::NonZeroU64;

            let size = u64::try_from(self.serializer.file_size()?).unwrap_or(0);
            let file = self.serializer.clone_file()?;
            // Upstream passes SEQUENTIAL | WILLNEED. POSIX advice values are
            // an enum rather than bitflags, and on the supported Unix targets
            // that expression resolves to WILLNEED. Preserve that one-call
            // behavior while using a typed descriptor API.
            let _ = rustix::fs::fadvise(
                &file,
                0,
                NonZeroU64::new(size),
                rustix::fs::Advice::WillNeed,
            );
        }
        Ok(())
    }

    pub fn file_path(&self) -> &Path {
        &self.file_name
    }
    #[cfg(unix)]
    pub fn file_name_bytes(&self) -> &[u8] {
        use std::os::unix::ffi::OsStrExt;
        self.file_name.as_os_str().as_bytes()
    }
    pub fn try_file_name(&self) -> IoResult<&str> {
        self.file_name
            .to_str()
            .ok_or_else(|| IoError::Other("OS path is not UTF-8; use file_path".into()))
    }
    pub fn file_name(&self) -> String {
        self.try_file_name()
            .expect("Non-UTF8 OS path; use file_path")
            .to_owned()
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

    /// Owned descriptor view, without reopening the pathname. Rewind first
    /// when transferring a completed output to an input reader.
    #[cfg(any(unix, windows))]
    pub(crate) fn shared_stdio(&mut self) -> IoResult<Option<super::cstdio_file::SharedCstdio>> {
        self.serializer.shared_stdio()
    }

    pub fn clone_file(&mut self) -> IoResult<StdFile> {
        self.serializer.clone_file()
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
            crate::util::test_path_component()
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
