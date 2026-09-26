//! Translation of `diamond/src/util/io/input_file.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::{Seek, SeekFrom};
#[cfg(unix)]
use std::os::fd::FromRawFd;

use super::{
    Compressor, Deserializer, FilePrimitive, FileSource, InputStreamBuffer, IoError, IoResult,
    OutputFile, StreamEntity, TempFile, TempFileData, VecStream, ZlibSource, ZstdSource,
};

pub fn detect_compressor(bytes: &[u8]) -> Compressor {
    if bytes.len() >= 2
        && ((bytes[0] == 0x1f && bytes[1] == 0x8b)
            || (bytes[0] == 0x78 && matches!(bytes[1], 0x01 | 0x9c | 0xda)))
    {
        Compressor::Zlib
    } else if bytes.len() >= 4 && bytes[..4] == [0x28, 0xb5, 0x2f, 0xfd] {
        Compressor::Zstd
    } else {
        Compressor::None
    }
}

#[derive(Debug)]
pub enum InputFileStream {
    File(InputStreamBuffer<FileSource>),
    Zlib(InputStreamBuffer<ZlibSource<FileSource>>),
    Zstd(InputStreamBuffer<ZstdSource<FileSource>>),
    /// Compatibility variant for existing callers that construct an input
    /// file over owned memory. C++ constructors use the three variants above.
    Memory(VecStream),
}

fn make_decompressor(
    compressor: Compressor,
    buffer: InputStreamBuffer<FileSource>,
) -> IoResult<InputFileStream> {
    match compressor {
        Compressor::Zlib => Ok(InputFileStream::Zlib(InputStreamBuffer::new(
            ZlibSource::new(buffer)?,
            0,
        ))),
        Compressor::Zstd => Ok(InputFileStream::Zstd(InputStreamBuffer::new(
            ZstdSource::new(buffer)?,
            0,
        ))),
        Compressor::None => Err(IoError::Other(String::new())),
    }
}

macro_rules! dispatch_stream {
    ($self:expr, $stream:ident => $expression:expr) => {
        match $self {
            Self::File($stream) => $expression,
            Self::Zlib($stream) => $expression,
            Self::Zstd($stream) => $expression,
            Self::Memory($stream) => $expression,
        }
    };
}

impl StreamEntity for InputFileStream {
    fn rewind(&mut self) -> IoResult<()> {
        dispatch_stream!(self, stream => stream.rewind())
    }

    fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        dispatch_stream!(self, stream => stream.seek(position, origin))
    }

    fn tell(&mut self) -> IoResult<i64> {
        dispatch_stream!(self, stream => stream.tell())
    }

    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        dispatch_stream!(self, stream => stream.read(ptr))
    }

    fn fetch(&mut self) -> IoResult<bool> {
        dispatch_stream!(self, stream => stream.fetch())
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

    fn data(&self) -> &[u8] {
        dispatch_stream!(self, stream => stream.data())
    }

    fn eof(&mut self) -> IoResult<bool> {
        dispatch_stream!(self, stream => stream.eof())
    }

    fn file_size(&mut self) -> IoResult<i64> {
        dispatch_stream!(self, stream => stream.file_size())
    }

    fn seekable(&self) -> bool {
        dispatch_stream!(self, stream => stream.seekable())
    }
}

#[derive(Debug)]
pub struct InputFile {
    pub(crate) deserializer: Deserializer<InputFileStream>,
    pub file_name: String,
    pub unlinked: bool,
    pub temp_file: bool,
}

impl InputFile {
    pub const BUFFERED: i32 = 1;
    pub const NO_AUTODETECT: i32 = 2;

    /// C++ `InputFile(const string&, int)`.
    pub fn new(file_name: &str, flags: i32) -> IoResult<Self> {
        let mut buffer = InputStreamBuffer::new(FileSource::new(file_name)?, flags);
        let mut compressor = None;

        if !file_name.is_empty() && file_name != "-" {
            #[cfg(not(windows))]
            {
                let metadata = std::fs::metadata(file_name).map_err(|_| {
                    IoError::Other(format!("Error calling stat on file {file_name}"))
                })?;
                if metadata.is_file() && flags & Self::NO_AUTODETECT == 0 {
                    buffer.fetch()?;
                    if buffer.data().len() >= 4 {
                        let detected = detect_compressor(buffer.data());
                        if detected != Compressor::None {
                            compressor = Some(detected);
                        }
                    }
                }
            }
            #[cfg(windows)]
            if flags & Self::NO_AUTODETECT == 0 {
                buffer.fetch()?;
                if buffer.data().len() >= 4 {
                    let detected = detect_compressor(buffer.data());
                    if detected != Compressor::None {
                        compressor = Some(detected);
                    }
                }
            }
        }

        let stream = match compressor {
            Some(compressor) => make_decompressor(compressor, buffer)?,
            None => InputFileStream::File(buffer),
        };

        Ok(Self {
            deserializer: Deserializer::new(stream),
            file_name: file_name.to_string(),
            unlinked: false,
            temp_file: false,
        })
    }

    /// Rust compatibility constructor over the descriptor container used by
    /// `TempFile`. The C++ constructor itself accepts `TempFile&`.
    pub fn from_temp_file_data(
        data: &TempFileData,
        flags: i32,
        compressor: Compressor,
    ) -> IoResult<Self> {
        #[cfg(unix)]
        let source = {
            let fd = unsafe { super::dup(data.fd) };
            if fd < 0 {
                return Err(IoError::Other(format!(
                    "Error opening temporary file {}",
                    data.name
                )));
            }
            let mut file = unsafe { StdFile::from_raw_fd(fd) };
            file.seek(SeekFrom::Start(0))
                .map_err(|_| IoError::Other("Error calling fseek.".to_string()))?;
            FileSource::from_file(&data.name, file)
        };
        #[cfg(not(unix))]
        let source = FileSource::new(&data.name)?;

        Self::from_source(source, data.name.clone(), data.unlinked, flags, compressor)
    }

    /// C++ `InputFile(TempFile&, int, Compressor)`.
    pub fn from_temp_file(
        tmp_file: &mut TempFile,
        flags: i32,
        compressor: Compressor,
    ) -> IoResult<Self> {
        let file_name = tmp_file.file_name();
        let file = tmp_file
            .file()?
            .try_clone()
            .map_err(|_| IoError::FileOpen(file_name.clone()))?;
        tmp_file.rewind()?;
        Self::from_source(
            FileSource::from_file(&file_name, file),
            file_name,
            tmp_file.unlinked(),
            flags,
            compressor,
        )
    }

    /// C++ `InputFile(OutputFile&, int)`.
    pub fn from_output_file(out_file: &mut OutputFile, flags: i32) -> IoResult<Self> {
        let file_name = out_file.file_name();
        let file = out_file
            .file()?
            .try_clone()
            .map_err(|_| IoError::FileOpen(file_name.clone()))?;
        out_file.rewind()?;
        Self::from_source(
            FileSource::from_file(&file_name, file),
            file_name,
            false,
            flags,
            Compressor::None,
        )
    }

    fn from_source(
        source: FileSource,
        file_name: String,
        unlinked: bool,
        flags: i32,
        compressor: Compressor,
    ) -> IoResult<Self> {
        let buffer = InputStreamBuffer::new(source, flags);
        let stream = if compressor == Compressor::None {
            InputFileStream::File(buffer)
        } else {
            make_decompressor(compressor, buffer)?
        };
        Ok(Self {
            deserializer: Deserializer::new(stream),
            file_name,
            unlinked,
            temp_file: true,
        })
    }

    pub fn close_and_delete(&mut self) -> IoResult<()> {
        self.close()?;
        if !self.unlinked && std::fs::remove_file(&self.file_name).is_err() {
            eprintln!(
                "Warning: Failed to delete temporary file {}",
                self.file_name
            );
        }
        Ok(())
    }

    pub fn hash(&mut self) -> IoResult<u64> {
        let mut seed = [0u8; 16];
        let mut buffer = [0u8; 4096];
        loop {
            let count = self.read_raw(&mut buffer)?;
            if count == 0 {
                break;
            }
            seed = crate::util::hash::murmurhash3_x64_128(&buffer[..count], &seed);
        }
        Ok(u64::from_ne_bytes(seed[..8].try_into().unwrap()))
    }

    pub(crate) fn into_deserializer(self) -> Deserializer<InputFileStream> {
        self.deserializer
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.deserializer.rewind()
    }

    pub fn seek(&mut self, position: i64) -> IoResult<&mut Self> {
        self.deserializer.seek(position)?;
        Ok(self)
    }

    pub fn seek_forward(&mut self, count: usize) -> IoResult<()> {
        self.deserializer.seek_forward(count)
    }

    pub fn seek_forward_delim(&mut self, delimiter: u8) -> IoResult<bool> {
        self.deserializer.seek_forward_delim(delimiter)
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.deserializer.close()
    }

    pub fn read_u32(&mut self) -> IoResult<u32> {
        self.deserializer.read_u32()
    }

    pub fn read_i32(&mut self) -> IoResult<i32> {
        self.deserializer.read_i32()
    }

    pub fn read_i16(&mut self) -> IoResult<i16> {
        self.deserializer.read_i16()
    }

    pub fn read_u16(&mut self) -> IoResult<u16> {
        self.deserializer.read_u16()
    }

    pub fn read_i64(&mut self) -> IoResult<i64> {
        self.deserializer.read_i64()
    }

    pub fn read_u64(&mut self) -> IoResult<u64> {
        self.deserializer.read_u64()
    }

    pub fn read_f64(&mut self) -> IoResult<f64> {
        self.deserializer.read_f64()
    }

    pub fn read_value<T: FilePrimitive>(&mut self) -> IoResult<T> {
        self.deserializer.read_value()
    }

    pub fn read_values<T: FilePrimitive>(&mut self, ptr: &mut [T]) -> IoResult<()> {
        self.deserializer.read_values(ptr)
    }

    pub fn read_string(&mut self) -> IoResult<String> {
        self.deserializer.read_string()
    }

    pub fn read_raw(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        self.deserializer.read_raw(ptr)
    }

    pub fn read_to(&mut self, dst: &mut Vec<u8>, delimiter: u8) -> IoResult<bool> {
        self.deserializer.read_to(dst, delimiter)
    }

    pub fn peek(&mut self, count: usize) -> IoResult<String> {
        self.deserializer.peek(count)
    }

    pub fn file_size(&mut self) -> IoResult<i64> {
        self.deserializer.file_size()
    }

    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.deserializer.file()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn fixture(name: &str, bytes: &[u8]) -> (std::path::PathBuf, String) {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-input-file-{name}-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
        ));
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, bytes).unwrap();
        let file_name = path.to_string_lossy().into_owned();
        (path, file_name)
    }

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn detects_exact_magic_values_and_short_inputs_safely() {
        assert_eq!(detect_compressor(&[0x1f, 0x8b]), Compressor::Zlib);
        for second in [0x01, 0x9c, 0xda] {
            assert_eq!(detect_compressor(&[0x78, second]), Compressor::Zlib);
        }
        assert_eq!(
            detect_compressor(&[0x28, 0xb5, 0x2f, 0xfd]),
            Compressor::Zstd
        );
        assert_eq!(detect_compressor(&[]), Compressor::None);
        assert_eq!(detect_compressor(&[0x28, 0xb5, 0x2f]), Compressor::None);
    }

    #[test]
    fn named_constructor_streams_autodetected_formats_and_honors_disable_flag() {
        let plain = b"streamed gzip payload";
        let encoded = gzip(plain);
        let (path, name) = fixture("autodetect", &encoded);

        let mut input = InputFile::new(&name, 0).unwrap();
        // A streaming source retains the underlying file delegation; the old
        // eager in-memory implementation did not.
        assert!(input.file().is_ok());
        assert!(matches!(input.deserializer.data(), []));
        let mut decoded = vec![0; plain.len()];
        assert_eq!(input.read_raw(&mut decoded).unwrap(), plain.len());
        assert_eq!(decoded, plain);
        input.rewind().unwrap();
        assert_eq!(input.read_raw(&mut decoded[..8]).unwrap(), 8);
        assert_eq!(&decoded[..8], &plain[..8]);

        let mut raw = InputFile::new(&name, InputFile::NO_AUTODETECT).unwrap();
        let mut magic = [0; 2];
        raw.read_raw(&mut magic).unwrap();
        assert_eq!(magic, [0x1f, 0x8b]);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn short_and_plain_files_preserve_bytes_already_fetched_for_detection() {
        for (name, bytes) in [("short", &b"abc"[..]), ("plain", &b"plain data"[..])] {
            let (path, file_name) = fixture(name, bytes);
            let mut input = InputFile::new(&file_name, 0).unwrap();
            let mut output = vec![0; bytes.len()];
            assert_eq!(input.read_raw(&mut output).unwrap(), bytes.len());
            assert_eq!(output, bytes);
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn temp_and_output_constructors_rewind_and_set_public_state() {
        let mut temp = TempFile::new(true).unwrap();
        temp.write_raw(b"temp").unwrap();
        let mut from_temp = InputFile::from_temp_file(&mut temp, 0, Compressor::None).unwrap();
        assert!(from_temp.temp_file);
        assert!(from_temp.unlinked);
        let mut bytes = [0; 4];
        from_temp.read_raw(&mut bytes).unwrap();
        assert_eq!(&bytes, b"temp");

        let (path, name) = fixture("output", b"");
        let mut output = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
        output.write_raw(b"output").unwrap();
        let mut from_output = InputFile::from_output_file(&mut output, 0).unwrap();
        assert!(from_output.temp_file);
        assert!(!from_output.unlinked);
        let mut bytes = [0; 6];
        from_output.read_raw(&mut bytes).unwrap();
        assert_eq!(&bytes, b"output");
        output.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn explicit_temp_compressor_is_streamed_and_rewindable() {
        let mut temp = TempFile::new(true).unwrap();
        temp.write_raw(&gzip(b"compressed temp")).unwrap();
        let mut input = InputFile::from_temp_file(&mut temp, 0, Compressor::Zlib).unwrap();
        let mut bytes = [0; 15];
        input.read_raw(&mut bytes).unwrap();
        assert_eq!(&bytes, b"compressed temp");
        input.rewind().unwrap();
        input.read_raw(&mut bytes[..10]).unwrap();
        assert_eq!(&bytes[..10], b"compressed");
    }

    #[test]
    fn hash_is_chained_in_4096_byte_chunks_from_current_position() {
        let bytes: Vec<_> = (0..9000).map(|index| (index % 251) as u8).collect();
        let (path, name) = fixture("hash", &bytes);
        let mut input = InputFile::new(&name, InputFile::NO_AUTODETECT).unwrap();
        input.seek(17).unwrap();
        assert_eq!(
            input.hash().unwrap(),
            crate::util::hash::file_hash(&bytes[17..])
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn close_and_delete_respects_linked_and_unlinked_state() {
        let (path, name) = fixture("delete", b"delete me");
        let mut input = InputFile::new(&name, InputFile::NO_AUTODETECT).unwrap();
        input.close_and_delete().unwrap();
        assert!(!path.exists());

        let mut temp = TempFile::new(true).unwrap();
        temp.write_raw(b"unlinked").unwrap();
        let mut input = InputFile::from_temp_file(&mut temp, 0, Compressor::None).unwrap();
        input.close_and_delete().unwrap();
    }
}
