//! Translation of `diamond/src/util/io/input_file.{h,cpp}`.

use std::fs::File as StdFile;
#[cfg(unix)]
use std::io::Seek;
use std::io::SeekFrom;
use std::path::{Path, PathBuf};

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

fn make_decompressor_with_buffer_size(
    compressor: Compressor,
    buffer: InputStreamBuffer<FileSource>,
    buffer_size: usize,
) -> IoResult<InputFileStream> {
    match compressor {
        Compressor::Zlib => Ok(InputFileStream::Zlib(InputStreamBuffer::with_buffer_size(
            ZlibSource::new(buffer)?,
            0,
            buffer_size,
        ))),
        Compressor::Zstd => Ok(InputFileStream::Zstd(InputStreamBuffer::with_buffer_size(
            ZstdSource::new(buffer)?,
            0,
            buffer_size,
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
    pub file_name: PathBuf,
    pub unlinked: bool,
    pub temp_file: bool,
}

impl std::io::Read for InputFile {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        self.read_raw(bytes).map_err(std::io::Error::other)
    }
}

impl InputFile {
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
    pub const BUFFERED: i32 = 1;
    pub const NO_AUTODETECT: i32 = 2;

    /// C++ `InputFile(const string&, int)`.
    pub fn new(file_name: &str, flags: i32) -> IoResult<Self> {
        Self::new_path(Path::new(file_name), flags)
    }
    pub fn new_path(file_name: &Path, flags: i32) -> IoResult<Self> {
        Self::with_buffer_size_path(file_name, flags, super::DEFAULT_FILE_BUFFER_SIZE)
    }

    /// Explicit replacement for Source's process-global file buffer setting.
    pub fn with_buffer_size(file_name: &str, flags: i32, buffer_size: usize) -> IoResult<Self> {
        Self::with_buffer_size_path(Path::new(file_name), flags, buffer_size)
    }
    pub fn with_buffer_size_path(
        file_name: &Path,
        flags: i32,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Self::with_codec_profile(file_name, flags, buffer_size, true, false)
    }

    /// Source command build profile, preserving the same stream owner/type.
    pub fn with_source_profile(file_name: &str, flags: i32, buffer_size: usize) -> IoResult<Self> {
        Self::with_source_profile_path(Path::new(file_name), flags, buffer_size)
    }
    pub fn with_source_profile_path(
        file_name: &Path,
        flags: i32,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Self::with_codec_profile(file_name, flags, buffer_size, cfg!(feature = "zstd"), true)
    }
    pub fn source_new(file_name: &str, flags: i32) -> IoResult<Self> {
        Self::source_new_path(Path::new(file_name), flags)
    }
    pub fn source_new_path(file_name: &Path, flags: i32) -> IoResult<Self> {
        Self::with_source_profile_path(file_name, flags, super::DEFAULT_FILE_BUFFER_SIZE)
    }
    fn with_codec_profile(
        file_name: &Path,
        flags: i32,
        buffer_size: usize,
        allow_zstd: bool,
        source_profile: bool,
    ) -> IoResult<Self> {
        let source = FileSource::new_path(file_name)?;
        if buffer_size == 0 && !source_profile {
            // Source's zero-capacity regular-file compressor peek fails
            // without setting EOF (fread(..., 0) performs no read).
            return Err(IoError::Other("Invalid peek".into()));
        }
        Vec::<u8>::new().try_reserve_exact(buffer_size).map_err(|_| IoError::Other(
            "Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage.".into()
        ))?;
        let mut buffer = InputStreamBuffer::with_source_buffer_size(source, flags, buffer_size);
        let mut compressor = None;

        if !file_name.as_os_str().is_empty() && file_name.as_os_str() != std::ffi::OsStr::new("-") {
            #[cfg(not(windows))]
            {
                let metadata = std::fs::metadata(super::source_c_path(file_name).as_ref())
                    .map_err(|_| {
                        IoError::path_message(b"Error calling stat on file ", file_name)
                    })?;
                if metadata.is_file() && flags & Self::NO_AUTODETECT == 0 {
                    buffer.fetch()?;
                    if buffer.data().len() < 4 && !buffer.eof()? {
                        return Err(IoError::Other("Invalid peek".into()));
                    }
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
                if buffer.data().len() < 4 && !buffer.eof()? {
                    return Err(IoError::Other("Invalid peek".into()));
                }
                if buffer.data().len() >= 4 {
                    let detected = detect_compressor(buffer.data());
                    if detected != Compressor::None {
                        compressor = Some(detected);
                    }
                }
            }
        }

        let stream = match compressor {
            Some(compressor) => {
                if compressor == Compressor::Zstd && !allow_zstd {
                    return Err(IoError::Other(
                        "Executable was not compiled with ZStd support.".into(),
                    ));
                }
                make_decompressor_with_buffer_size(compressor, buffer, buffer_size)?
            }
            None => InputFileStream::File(buffer),
        };

        Ok(Self {
            deserializer: Deserializer::new(stream),
            file_name: file_name.to_path_buf(),
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
            let mut file = data.try_clone_file()?;
            file.seek(SeekFrom::Start(0))
                .map_err(|_| IoError::Other("Error calling fseek.".to_string()))?;
            FileSource::try_from_file_path(&data.name, file)?
        };
        #[cfg(not(unix))]
        let source = FileSource::new_path(&data.name)?;

        Self::from_source(source, data.name.clone(), data.unlinked, flags, compressor)
    }

    /// C++ `InputFile(TempFile&, int, Compressor)`.
    pub fn from_temp_file(
        tmp_file: &mut TempFile,
        flags: i32,
        compressor: Compressor,
    ) -> IoResult<Self> {
        Self::from_temp_file_with_buffer_size(
            tmp_file,
            flags,
            compressor,
            super::DEFAULT_FILE_BUFFER_SIZE,
        )
    }
    pub fn from_temp_file_with_buffer_size(
        tmp_file: &mut TempFile,
        flags: i32,
        compressor: Compressor,
        buffer_size: usize,
    ) -> IoResult<Self> {
        let file_name = tmp_file.file_path().to_path_buf();
        #[cfg(any(unix, windows))]
        let shared = tmp_file.shared_stdio()?;
        tmp_file.rewind()?;
        #[cfg(any(unix, windows))]
        let source = if let Some(shared) = shared {
            FileSource::from_shared_path(&file_name, shared)
        } else {
            FileSource::try_from_file_path(&file_name, tmp_file.clone_file()?)?
        };
        #[cfg(not(any(unix, windows)))]
        let source = FileSource::from_file_path(&file_name, tmp_file.clone_file()?);
        Self::from_source_with_buffer_size(
            source,
            file_name,
            tmp_file.unlinked(),
            flags,
            compressor,
            buffer_size,
        )
    }

    /// C++ `InputFile(OutputFile&, int)`.
    pub fn from_output_file(out_file: &mut OutputFile, flags: i32) -> IoResult<Self> {
        Self::from_output_file_with_buffer_size(out_file, flags, super::DEFAULT_FILE_BUFFER_SIZE)
    }

    pub fn from_output_file_with_buffer_size(
        out_file: &mut OutputFile,
        flags: i32,
        buffer_size: usize,
    ) -> IoResult<Self> {
        let file_name = out_file.file_path().to_path_buf();
        #[cfg(any(unix, windows))]
        let shared = out_file.shared_stdio()?;
        out_file.rewind()?;
        #[cfg(any(unix, windows))]
        let source = if let Some(shared) = shared {
            FileSource::from_shared_path(&file_name, shared)
        } else {
            FileSource::try_from_file_path(&file_name, out_file.clone_file()?)?
        };
        #[cfg(not(any(unix, windows)))]
        let source = FileSource::from_file_path(&file_name, out_file.clone_file()?);
        Self::from_source_with_buffer_size(
            source,
            file_name,
            false,
            flags,
            Compressor::None,
            buffer_size,
        )
    }

    fn from_source(
        source: FileSource,
        file_name: PathBuf,
        unlinked: bool,
        flags: i32,
        compressor: Compressor,
    ) -> IoResult<Self> {
        Self::from_source_with_buffer_size(
            source,
            file_name,
            unlinked,
            flags,
            compressor,
            super::DEFAULT_FILE_BUFFER_SIZE,
        )
    }
    fn from_source_with_buffer_size(
        source: FileSource,
        file_name: PathBuf,
        unlinked: bool,
        flags: i32,
        compressor: Compressor,
        buffer_size: usize,
    ) -> IoResult<Self> {
        let buffer = InputStreamBuffer::with_source_buffer_size(source, flags, buffer_size);
        let stream = if compressor == Compressor::None {
            InputFileStream::File(buffer)
        } else {
            make_decompressor_with_buffer_size(compressor, buffer, buffer_size)?
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
        if !self.unlinked
            && std::fs::remove_file(super::source_c_path(&self.file_name).as_ref()).is_err()
        {
            let mut bytes = b"Warning: Failed to delete temporary file ".to_vec();
            bytes.extend_from_slice(&super::path_diagnostic_bytes(&self.file_name)?);
            bytes.push(b'\n');
            let _ = std::io::Write::write_all(&mut super::console::stderr(), &bytes);
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

    /// Codec-backed inputs remain non-seekable even with prefetched data.
    pub fn seekable(&self) -> bool {
        self.deserializer.seekable()
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

    pub fn peek_bytes(&mut self, count: usize) -> IoResult<Vec<u8>> {
        self.deserializer.peek_bytes(count)
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
    #[cfg(unix)]
    #[test]
    fn os_paths_preserve_non_utf8_through_codecs_errors_and_temp_handoff() {
        use std::os::unix::ffi::{OsStrExt, OsStringExt};
        let root = std::env::temp_dir().join(format!(
            "diamond-raw-path-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ));
        std::fs::create_dir(&root).unwrap();
        let (dir, raw_names) =
            crate::util::test_os_path_or_unicode(&root, b"directory\xff", "directory-é");
        std::fs::create_dir(&dir).unwrap();
        let path = if raw_names {
            dir.join(std::ffi::OsString::from_vec(b"file\x80.gz".to_vec()))
        } else {
            dir.join("file-λ.gz")
        };
        let mut output =
            super::OutputFile::with_source_profile_path(&path, super::Compressor::Zlib, "wb", 32)
                .unwrap();
        assert_eq!(output.file_path(), path);
        if raw_names {
            assert!(output.try_file_name().is_err());
        } else {
            assert_eq!(
                output.try_file_name().unwrap().as_bytes(),
                path.as_os_str().as_bytes()
            );
        }
        output.write_raw(b"actual raw path").unwrap();
        output.close().unwrap();
        let mut input = super::InputFile::with_source_profile_path(&path, 0, 32).unwrap();
        assert_eq!(input.file_name_bytes(), path.as_os_str().as_bytes());
        let mut bytes = [0; 15];
        assert_eq!(input.read_raw(&mut bytes).unwrap(), 15);
        assert_eq!(&bytes, b"actual raw path");
        input.close().unwrap();
        let missing = dir.join(std::ffi::OsString::from_vec(b"missing\xfe".to_vec()));
        let error = super::InputFile::new_path(&missing, 0).unwrap_err();
        let saved = crate::util::io::file_open_error::find(&error).unwrap();
        assert_eq!(saved.path, missing);
        let mut expected = b"Error opening file ".to_vec();
        expected.extend_from_slice(missing.as_os_str().as_bytes());
        assert_eq!(saved.try_message_bytes().unwrap(), expected);
        let data = crate::util::io::TempFileData::init_in(true, &dir, false).unwrap();
        assert_eq!(data.name.parent(), Some(dir.as_path()));
        assert!(data.unlinked);
        assert!(!data.name.exists());
        let mut temp = crate::util::io::TempFile::from_temp_file_data(&data).unwrap();
        temp.write_raw(b"abc").unwrap();
        let mut alias = super::InputFile::from_temp_file_with_buffer_size(
            &mut temp,
            0,
            super::Compressor::None,
            32,
        )
        .unwrap();
        assert_eq!(alias.file_path(), data.name);
        let mut bytes = [0; 3];
        alias.read_raw(&mut bytes).unwrap();
        assert_eq!(&bytes, b"abc");
        alias.close().unwrap();
        drop(alias);
        drop(temp);
        drop(data);
        drop(input);
        drop(output);
        std::fs::remove_dir_all(root).unwrap();
    }

    use super::*;
    use std::io::Write;

    fn fixture(name: &str, bytes: &[u8]) -> (std::path::PathBuf, String) {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-input-file-{name}-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
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
        #[cfg(not(unix))]
        let temp_name = temp.file_name().to_string();
        temp.write_raw(b"temp").unwrap();
        let mut from_temp = InputFile::from_temp_file(&mut temp, 0, Compressor::None).unwrap();
        assert!(from_temp.temp_file);
        assert_eq!(from_temp.unlinked, cfg!(unix));
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
        from_output.close().unwrap();
        output.close().unwrap();
        std::fs::remove_file(path).unwrap();

        from_temp.close().unwrap();
        temp.close().unwrap();
        #[cfg(not(unix))]
        std::fs::remove_file(temp_name).unwrap();
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
