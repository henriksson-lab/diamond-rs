use std::fs::File as StdFile;
use std::io::SeekFrom;
#[cfg(test)]
use std::io::Write;

pub mod deserializer;
pub use deserializer::Deserializer;
pub mod serializer;
pub use serializer::Serializer;
pub mod file;
pub use file::{File, FilePrimitive, Temporary};
pub mod file_source;
pub use file_source::FileSource;
pub mod file_sink;
pub use file_sink::FileSink;
pub mod temp_file;
pub use temp_file::{TempFile, TempFileData, TempFileHandler, TEMP_FILE_HANDLER};
pub mod output_file;
pub use output_file::{decompress, decompress_file, Compressor, OutputFile};
pub mod zstd_stream;
#[cfg(test)]
pub(crate) use zstd_stream::zstd_compress_for_test;
pub use zstd_stream::{zstd_decompress, zstd_decompress_file, ZstdSink, ZstdSource};
pub mod compressed_stream;
pub use compressed_stream::{zlib_decompress, zlib_decompress_file, ZlibSink, ZlibSource};
pub mod output_stream_buffer;
pub use output_stream_buffer::OutputStreamBuffer;
pub mod input_stream_buffer;
pub use input_stream_buffer::InputStreamBuffer;
pub mod text_input_file;
pub use text_input_file::TextInputFile;
pub mod input_file;
pub use input_file::{detect_compressor, InputFile, InputFileStream};

const DEFAULT_FILE_BUFFER_SIZE: usize = 1 << 20;
pub const MEGABYTES: usize = 1 << 20;
pub const GIGABYTES: usize = 1 << 30;
pub const KILOBYTES: usize = 1 << 10;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IoError {
    UnsupportedOperation,
    FileOpen(String),
    FileRead(String),
    FileWrite(String),
    EndOfStream,
    StreamRead { line_count: usize, msg: String },
    Other(String),
}

impl std::fmt::Display for IoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IoError::UnsupportedOperation => f.write_str("Unsupported I/O operation."),
            IoError::FileOpen(file_name) => write!(f, "Error opening file {file_name}"),
            IoError::FileRead(file_name) => write!(f, "Error reading file {file_name}"),
            IoError::FileWrite(file_name) => write!(f, "Error writing file {file_name}"),
            IoError::EndOfStream => f.write_str("Unexpected end of input."),
            IoError::StreamRead { line_count, msg } => {
                write!(f, "Error reading input stream at line {line_count}: {msg}")
            }
            IoError::Other(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for IoError {}

pub type IoResult<T> = Result<T, IoError>;

pub trait Consumer {
    fn consume(&mut self, ptr: &[u8]) -> IoResult<()>;
    fn finalize(&mut self) -> IoResult<()> {
        Ok(())
    }
}

pub trait StreamEntity {
    fn rewind(&mut self) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn seek(&mut self, _p: i64, _origin: SeekFrom) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn tell(&mut self) -> IoResult<i64> {
        Err(IoError::UnsupportedOperation)
    }

    fn read(&mut self, _ptr: &mut [u8]) -> IoResult<usize> {
        Err(IoError::UnsupportedOperation)
    }

    fn fetch(&mut self) -> IoResult<bool> {
        Err(IoError::UnsupportedOperation)
    }

    fn close(&mut self) -> IoResult<()> {
        Ok(())
    }

    fn file_name(&self) -> &str {
        ""
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        Err(IoError::UnsupportedOperation)
    }

    fn data(&self) -> &[u8] {
        &[]
    }

    fn write(&mut self, _ptr: &[u8]) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn flush(&mut self) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn eof(&mut self) -> IoResult<bool> {
        Err(IoError::UnsupportedOperation)
    }

    fn file_size(&mut self) -> IoResult<i64> {
        Err(IoError::UnsupportedOperation)
    }

    fn seekable(&self) -> bool {
        false
    }
}

pub mod compressed_buffer;
pub use compressed_buffer::CompressedBuffer;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VecStream {
    data: Vec<u8>,
    pos: usize,
    name: String,
}

impl VecStream {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_vec(data: Vec<u8>) -> Self {
        Self {
            data,
            pos: 0,
            name: String::new(),
        }
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }
}

impl StreamEntity for VecStream {
    fn rewind(&mut self) -> IoResult<()> {
        self.pos = 0;
        Ok(())
    }

    fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        let pos = match origin {
            SeekFrom::Start(_) => p,
            SeekFrom::End(_) => self.data.len() as i64 + p,
            SeekFrom::Current(_) => self.pos as i64 + p,
        };
        if pos < 0 {
            return Err(IoError::UnsupportedOperation);
        }
        self.pos = pos as usize;
        if self.pos > self.data.len() {
            self.data.resize(self.pos, 0);
        }
        Ok(())
    }

    fn tell(&mut self) -> IoResult<i64> {
        Ok(self.pos as i64)
    }

    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let n = ptr.len().min(self.data.len().saturating_sub(self.pos));
        ptr[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }

    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        if self.pos + ptr.len() > self.data.len() {
            self.data.resize(self.pos + ptr.len(), 0);
        }
        self.data[self.pos..self.pos + ptr.len()].copy_from_slice(ptr);
        self.pos += ptr.len();
        Ok(())
    }

    fn flush(&mut self) -> IoResult<()> {
        Ok(())
    }

    fn eof(&mut self) -> IoResult<bool> {
        Ok(self.pos >= self.data.len())
    }

    fn file_size(&mut self) -> IoResult<i64> {
        Ok(self.data.len() as i64)
    }

    fn seekable(&self) -> bool {
        true
    }

    fn file_name(&self) -> &str {
        &self.name
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serializer_deserializer_roundtrip() {
        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        ser.write_i16(-2).unwrap();
        ser.write_u16(0x0102).unwrap();
        ser.write_i32(-7).unwrap();
        ser.write_u64(0x0102_0304_0506_0708).unwrap();
        ser.write_f64(1.5).unwrap();
        ser.write_value(-9i32).unwrap();
        ser.write_values(&[11u16, 12u16, 13u16]).unwrap();
        ser.write_string("abc").unwrap();
        let mut stream = ser.into_inner().unwrap();
        stream.rewind().unwrap();

        let mut de = Deserializer::new(stream);
        assert_eq!(de.read_i16().unwrap(), -2);
        assert_eq!(de.read_u16().unwrap(), 0x0102);
        assert_eq!(de.read_i32().unwrap(), -7);
        assert_eq!(de.read_u64().unwrap(), 0x0102_0304_0506_0708);
        assert_eq!(de.read_f64().unwrap(), 1.5);
        assert_eq!(de.read_value::<i32>().unwrap(), -9);
        let mut values = [0u16; 3];
        de.read_values(&mut values).unwrap();
        assert_eq!(values, [11, 12, 13]);
        assert_eq!(de.read_string().unwrap(), "abc");
    }

    #[test]
    fn test_serializer_buffered_write_seek_and_reset() {
        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        let data = vec![b'a'; DEFAULT_FILE_BUFFER_SIZE + 17];
        ser.write_raw(&data).unwrap();
        assert_eq!(ser.tell().unwrap(), DEFAULT_FILE_BUFFER_SIZE + 17);
        ser.seek(4, SeekFrom::Start(0)).unwrap();
        ser.write_raw(b"BC").unwrap();
        let stream = ser.into_inner().unwrap();
        assert_eq!(&stream.data()[..4], b"aaaa");
        assert_eq!(&stream.data()[4..6], b"BC");
        assert_eq!(stream.data().len(), DEFAULT_FILE_BUFFER_SIZE + 17);

        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        ser.reset_buffer();
        ser.write_raw(b"xy").unwrap();
        assert_eq!(ser.into_inner().unwrap().data(), b"xy");
    }

    #[test]
    fn test_deserializer_read_to_seek_and_peek() {
        let stream = VecStream::from_vec(b"abc,def\n>rec\nx".to_vec());
        let mut de = Deserializer::new(stream);
        assert!(de.data().is_empty());
        assert_eq!(de.peek(3).unwrap(), "abc");
        let mut out = Vec::new();
        assert!(de.read_to(&mut out, b',').unwrap());
        assert_eq!(out, b"abc");
        de.seek_forward(3).unwrap();
        let mut out = Vec::new();
        assert_eq!(
            de.read_to_record_start(&mut out, b'\n', b'>').unwrap(),
            (true, 1)
        );
        assert!(out.is_empty());
        assert_eq!(de.read_string().unwrap_err(), IoError::EndOfStream);

        let mut input = InputStreamBuffer::new(VecStream::from_vec(b"peek-pop".to_vec()), 0);
        input.fetch().unwrap();
        let mut de = Deserializer::new(input);
        assert_eq!(de.peek(4).unwrap(), "peek");
        let mut popped = [0u8; 4];
        de.pop(&mut popped).unwrap();
        assert_eq!(&popped, b"peek");
    }

    #[test]
    fn test_file_temporary_write_read_size() {
        let mut file = File::new_temporary(Temporary).unwrap();
        file.write(b"abcd").unwrap();
        file.write_value(7u16).unwrap();
        assert_eq!(file.size().unwrap(), 6);
        file.seek(0, SeekFrom::Start(0)).unwrap();
        let mut buf = [0u8; 4];
        file.read_exact(&mut buf).unwrap();
        assert_eq!(&buf, b"abcd");
        assert!(!file.eof().unwrap());
        assert_eq!(file.read(2).unwrap(), &7u16.to_ne_bytes());
        assert!(file.eof().unwrap());
        file.seek(4, SeekFrom::Start(0)).unwrap();
        assert_eq!(file.read_value::<u16>().unwrap(), 7);
        assert!(!file.file_name().is_empty());
        #[cfg(unix)]
        assert!(!std::path::Path::new(file.file_name()).exists());
        assert!(file.file().is_ok());
        file.close().unwrap();

        let path = std::env::temp_dir().join(format!(
            "diamond-rs-file-new-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut named = File::new(&name, "w+b").unwrap();
            named.write(b"x").unwrap();
            assert_eq!(named.file_name(), name);
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_stream_entity_file_delegation() {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "diamond-rs-stream-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut serializer = Serializer::new(OutputStreamBuffer::new(
                FileSink::new(&name, "w+b", false, 0).unwrap(),
            ));
            serializer.write_raw(b"abcdef").unwrap();
            assert!(serializer.file().is_ok());
            serializer.close().unwrap();
        }
        {
            let mut de =
                Deserializer::new(InputStreamBuffer::new(FileSource::new(&name).unwrap(), 0));
            assert!(de.file().is_ok());
            let mut b = [0u8; 2];
            de.read_exact(&mut b).unwrap();
            assert_eq!(&b, b"ab");
            de.close().unwrap();
        }
        let _ = std::fs::remove_file(&name);
    }

    #[test]
    fn test_errors_display() {
        assert_eq!(
            IoError::UnsupportedOperation.to_string(),
            "Unsupported I/O operation."
        );
        assert_eq!(IoError::EndOfStream.to_string(), "Unexpected end of input.");
        assert!(IoError::StreamRead {
            line_count: 3,
            msg: "bad".to_string()
        }
        .to_string()
        .contains("line 3"));
    }

    #[test]
    fn test_detect_compressor_and_decompress() {
        assert_eq!(detect_compressor(&[0x1f, 0x8b, 0, 0]), Compressor::Zlib);
        assert_eq!(detect_compressor(&[0x78, 0x9c, 0, 0]), Compressor::Zlib);
        assert_eq!(
            detect_compressor(&[0x28, 0xb5, 0x2f, 0xfd]),
            Compressor::Zstd
        );
        assert_eq!(detect_compressor(b"plain"), Compressor::None);

        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"abcabc").unwrap();
        let compressed = zlib_encoder.finish().unwrap();
        let mut out = [0u8; 16];
        let n = decompress(&compressed, &mut out, Compressor::Zlib).unwrap();
        assert_eq!(&out[..n], b"abcabc");

        let mut gzip_encoder =
            flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip_encoder.write_all(b"xyz").unwrap();
        let compressed = gzip_encoder.finish().unwrap();
        let n = decompress(&compressed, &mut out, Compressor::Zlib).unwrap();
        assert_eq!(&out[..n], b"xyz");
        let n = zlib_decompress(&compressed, &mut out).unwrap();
        assert_eq!(&out[..n], b"xyz");

        let compressed = zstd_compress_for_test(b"zstd-data");
        let n = decompress(&compressed, &mut out, Compressor::Zstd).unwrap();
        assert_eq!(&out[..n], b"zstd-data");
        let n = zstd_decompress(&compressed, &mut out).unwrap();
        assert_eq!(&out[..n], b"zstd-data");
        assert!(zstd_decompress(&compressed, &mut [0u8; 2]).is_err());

        assert!(decompress(&compressed, &mut [0u8; 2], Compressor::Zlib).is_err());

        let mut first = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        first.write_all(b"one").unwrap();
        let mut concatenated = first.finish().unwrap();
        let mut second =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        second.write_all(b"two").unwrap();
        concatenated.extend_from_slice(&second.finish().unwrap());
        let mut out = [0u8; 6];
        let n = zlib_decompress(&concatenated, &mut out).unwrap();
        assert_eq!(&out[..n], b"onetwo");
    }

    #[test]
    fn test_compressed_buffer_roundtrip_and_clear() {
        let mut compressed = CompressedBuffer::new();
        compressed.write(b"abc").unwrap();
        compressed.write_value(7u16).unwrap();
        compressed.write(b"xyzxyzxyz").unwrap();
        compressed.finish().unwrap();
        assert!(compressed.size() > 0);
        assert_eq!(detect_compressor(compressed.data()), Compressor::Zlib);

        let mut out = [0u8; 32];
        let n = decompress(compressed.data(), &mut out, Compressor::Zlib).unwrap();
        let mut expected = Vec::new();
        expected.extend_from_slice(b"abc");
        expected.extend_from_slice(&7u16.to_ne_bytes());
        expected.extend_from_slice(b"xyzxyzxyz");
        assert_eq!(&out[..n], expected.as_slice());

        compressed.clear();
        assert_eq!(compressed.size(), 0);
        compressed.write(b"second").unwrap();
        compressed.finish().unwrap();
        let n = decompress(compressed.data(), &mut out, Compressor::Zlib).unwrap();
        assert_eq!(&out[..n], b"second");
    }

    #[test]
    fn test_file_source_sink_and_stream_buffers() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let mut out = OutputStreamBuffer::new(sink);
            out.write(b"abc").unwrap();
            out.write(b"def").unwrap();
            assert_eq!(out.tell().unwrap(), 6);
            out.close().unwrap();
        }
        {
            let source = FileSource::new(&name).unwrap();
            assert!(source.seekable());
            let mut input = InputStreamBuffer::new(source, 0);
            let mut buf = [0u8; 6];
            assert_eq!(input.read(&mut buf).unwrap(), 6);
            assert_eq!(&buf, b"abcdef");
            assert_eq!(input.file_size().unwrap(), 6);
        }
        let _ = std::fs::remove_file(path);
    }

    #[cfg(unix)]
    #[test]
    fn test_file_source_sink_standard_stream_names() {
        let stdin_empty = FileSource::new("").unwrap();
        assert_eq!(stdin_empty.file_name(), "");
        assert!(!stdin_empty.seekable());

        let stdin_dash = FileSource::new("-").unwrap();
        assert_eq!(stdin_dash.file_name(), "-");
        assert!(!stdin_dash.seekable());

        let mut stdout_sink = FileSink::new("", "wb", false, 0).unwrap();
        assert_eq!(stdout_sink.file_name(), "");
        stdout_sink.close().unwrap();
    }

    #[test]
    fn test_file_sink_modes_and_from_file() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-fd-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        let file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&path)
            .unwrap();
        assert_eq!(
            FileSink::new(&name, "bad", false, 0)
                .unwrap_err()
                .to_string(),
            "Invalid fopen mode."
        );
        assert_eq!(
            FileSink::from_file(&name, file.try_clone().unwrap(), "bad", false, 0)
                .unwrap_err()
                .to_string(),
            "Invalid fopen mode."
        );
        let mut sink = FileSink::from_file(&name, file, "w+b", false, 0).unwrap();
        sink.write(b"fd").unwrap();
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"fd");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_output_stream_buffer_write_buffer_and_flush_count() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-output-buffer-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let mut out = OutputStreamBuffer::new(sink);
            let (buf, size) = out.write_buffer_range();
            assert_eq!(size, DEFAULT_FILE_BUFFER_SIZE);
            buf[..4].copy_from_slice(b"abcd");
            out.flush_count(4).unwrap();
            out.write(b"ef").unwrap();
            out.close().unwrap();
        }
        assert_eq!(std::fs::read(&path).unwrap(), b"abcdef");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_compressed_stream_entities_roundtrip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-compressed-stream-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let out = OutputStreamBuffer::new(sink);
            let mut zlib = ZlibSink::new(out);
            zlib.write(b"stream-zlib").unwrap();
            zlib.close().unwrap();
        }
        assert_eq!(
            detect_compressor(&std::fs::read(&path).unwrap()),
            Compressor::Zlib
        );
        {
            let source = FileSource::new(&name).unwrap();
            let input = InputStreamBuffer::new(source, 0);
            let mut zlib = ZlibSource::new(input).unwrap();
            let mut buf = [0u8; 11];
            assert_eq!(zlib.read(&mut buf).unwrap(), 11);
            assert_eq!(&buf, b"stream-zlib");
            let mut tail = [0u8; 1];
            assert_eq!(zlib.read(&mut tail).unwrap(), 0);
            assert!(zlib.eof().unwrap());
            zlib.rewind().unwrap();
            let mut head = [0u8; 6];
            assert_eq!(zlib.read(&mut head).unwrap(), 6);
            assert_eq!(&head, b"stream");
        }

        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let out = OutputStreamBuffer::new(sink);
            let mut zstd = ZstdSink::new(out);
            zstd.write(b"stream-zstd").unwrap();
            zstd.close().unwrap();
        }
        assert_eq!(
            detect_compressor(&std::fs::read(&path).unwrap()),
            Compressor::Zstd
        );
        {
            let source = FileSource::new(&name).unwrap();
            let input = InputStreamBuffer::new(source, 0);
            let mut zstd = ZstdSource::new(input).unwrap();
            let mut buf = [0u8; 11];
            assert_eq!(zstd.read(&mut buf).unwrap(), 11);
            assert_eq!(&buf, b"stream-zstd");
            let mut tail = [0u8; 1];
            assert_eq!(zstd.read(&mut tail).unwrap(), 0);
            assert!(zstd.eof().unwrap());
        }

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_temp_file_data_and_temp_file() {
        let data = TempFileData::init(false).unwrap();
        {
            let mut out = OutputFile::from_temp_file_data(&data, Compressor::None, "w+b").unwrap();
            out.write_raw(b"temp").unwrap();
            out.close().unwrap();
        }
        assert_eq!(std::fs::read(&data.name).unwrap(), b"temp");
        std::fs::remove_file(&data.name).unwrap();

        let dir = TempFile::get_temp_dir().unwrap();
        assert!(!dir.is_empty());

        let mut tmp = TempFile::new(true).unwrap();
        assert!(tmp.unlinked());
        tmp.write_raw(b"x").unwrap();
        tmp.close().unwrap();
    }

    #[test]
    fn test_output_file_input_file_and_hash() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-file-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
            out.write_i16(-123).unwrap();
            out.write_u16(0xCAFE).unwrap();
            out.write_i32(-7).unwrap();
            out.write_u64(0x0102_0304_0506_0708).unwrap();
            out.write_value(99u32).unwrap();
            out.write_values(&[3i16, 4i16]).unwrap();
            out.write_string("abc").unwrap();
            assert!(out.file().is_ok());
            assert!(out.tell().unwrap() > 0);
            out.close().unwrap();
        }
        {
            let mut input = InputFile::new(&name, 0).unwrap();
            assert_eq!(input.read_i16().unwrap(), -123);
            assert_eq!(input.read_u16().unwrap(), 0xCAFE);
            assert_eq!(input.read_i32().unwrap(), -7);
            assert_eq!(input.read_u64().unwrap(), 0x0102_0304_0506_0708);
            assert_eq!(input.read_value::<u32>().unwrap(), 99);
            let mut values = [0i16; 2];
            input.read_values(&mut values).unwrap();
            assert_eq!(values, [3, 4]);
            assert_eq!(input.read_string().unwrap(), "abc");
            assert!(input.file().is_ok());
            input.rewind().unwrap();
            let h = input.hash().unwrap();
            let bytes = std::fs::read(&name).unwrap();
            assert_eq!(h, crate::util::hash::file_hash(&bytes));
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_output_file_compressors_roundtrip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-compressor-out-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        {
            let mut out = OutputFile::new(&name, Compressor::Zlib, "w+b").unwrap();
            out.write_raw(b"gzip-output").unwrap();
            out.close().unwrap();
        }
        let compressed = std::fs::read(&path).unwrap();
        assert_eq!(detect_compressor(&compressed), Compressor::Zlib);
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"gzip-output");

        {
            let mut out = OutputFile::new(&name, Compressor::Zstd, "w+b").unwrap();
            out.write_raw(b"zstd-output").unwrap();
            out.close().unwrap();
        }
        let compressed = std::fs::read(&path).unwrap();
        assert_eq!(detect_compressor(&compressed), Compressor::Zstd);
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zstd-output");

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_input_file_temp_constructors_with_compressors() {
        let mut tmp = TempFile::new(true).unwrap();
        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"temp-zlib").unwrap();
        tmp.write_raw(&zlib_encoder.finish().unwrap()).unwrap();
        let mut input = InputFile::from_temp_file(&mut tmp, 0, Compressor::Zlib).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"temp-zlib");
        input.close().unwrap();
        tmp.close().unwrap();

        let data = TempFileData::init(false).unwrap();
        {
            let mut out = OutputFile::from_temp_file_data(&data, Compressor::Zstd, "w+b").unwrap();
            out.write_raw(b"temp-zstd").unwrap();
            out.close().unwrap();
        }
        let mut input = InputFile::from_temp_file_data(&data, 0, Compressor::Zstd).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"temp-zstd");
        input.close().unwrap();
        std::fs::remove_file(&data.name).unwrap();
    }

    #[test]
    fn test_input_file_autodetects_zlib_and_gzip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-compressed-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"zlib-data").unwrap();
        std::fs::write(&path, zlib_encoder.finish().unwrap()).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zlib-data");
        input.rewind().unwrap();
        let mut buf = [0u8; 4];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zlib");

        let mut first = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        first.write_all(b"one").unwrap();
        let mut concatenated = first.finish().unwrap();
        let mut second =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        second.write_all(b"two").unwrap();
        concatenated.extend_from_slice(&second.finish().unwrap());
        std::fs::write(&path, &concatenated).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 6];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"onetwo");

        let mut gzip_encoder =
            flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip_encoder.write_all(b"gzip-data").unwrap();
        let gzip_compressed = gzip_encoder.finish().unwrap();
        std::fs::write(&path, &gzip_compressed).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"gzip-data");

        std::fs::write(&path, zstd_compress_for_test(b"zstd-data")).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zstd-data");

        std::fs::write(&path, &gzip_compressed).unwrap();
        let mut input = InputFile::new(&name, InputFile::NO_AUTODETECT).unwrap();
        let mut header = [0u8; 2];
        input.read_raw(&mut header).unwrap();
        assert_eq!(header, [0x1f, 0x8b]);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_input_file_from_output_file() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-from-output-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
            out.write_raw(b"from-output").unwrap();
            let mut input = InputFile::from_output_file(&mut out, 0).unwrap();
            let mut buf = [0u8; 11];
            input.read_raw(&mut buf).unwrap();
            assert_eq!(&buf, b"from-output");
            assert!(input.temp_file);
            input.close().unwrap();
            out.close().unwrap();
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_input_file_from_unlinked_temp_file() {
        let mut tmp = TempFile::new(true).unwrap();
        assert!(tmp.unlinked());
        tmp.write_raw(b"hidden-temp").unwrap();
        let mut input = InputFile::from_temp_file(&mut tmp, 0, Compressor::None).unwrap();
        assert!(input.temp_file);
        assert!(input.unlinked);
        #[cfg(unix)]
        assert!(!std::path::Path::new(&input.file_name).exists());
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"hidden-temp");
        input.close().unwrap();
        tmp.close().unwrap();
    }

    #[test]
    fn test_text_input_file_getline_putback() {
        let stream = VecStream::from_vec(b"a\r\nb\n".to_vec());
        let de = Deserializer::new(stream);
        let mut text = TextInputFile::new(de, b'\n');
        text.getline().unwrap();
        assert_eq!(text.line, "a");
        assert_eq!(text.line_count, 1);
        text.putback_line();
        text.getline().unwrap();
        assert_eq!(text.line, "a");
        text.getline().unwrap();
        assert_eq!(text.line, "b");
        assert!(!text.eof());
        text.getline().unwrap();
        assert_eq!(text.line, "");
        assert!(text.eof());
        text.rewind().unwrap();
        assert_eq!(text.line_count, 0);
    }

    #[test]
    fn test_text_input_file_constructors() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-text-input-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        std::fs::write(&path, b"one\r\ntwo\n").unwrap();

        let mut text = TextInputFile::from_file_name(&name, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "one");
        text.getline().unwrap();
        assert_eq!(text.line, "two");

        let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
        out.write_raw(b"out\n").unwrap();
        let mut text = TextInputFile::from_output_file(&mut out, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "out");
        out.close().unwrap();

        let data = TempFileData::init(false).unwrap();
        let mut tmp = TempFile::from_temp_file_data(&data).unwrap();
        tmp.write_raw(b"tmp\n").unwrap();
        let mut text = TextInputFile::from_temp_file(&mut tmp, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "tmp");
        tmp.close().unwrap();
        std::fs::remove_file(&data.name).unwrap();

        let mut unlinked = TempFile::new(true).unwrap();
        unlinked.write_raw(b"unlinked\n").unwrap();
        let mut text = TextInputFile::from_temp_file(&mut unlinked, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "unlinked");
        unlinked.close().unwrap();

        std::fs::remove_file(path).unwrap();
    }
}
