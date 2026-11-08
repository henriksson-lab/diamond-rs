//! Bridges standard Rust readers to the translated Source stream/codec stack.
//! No decoder implementation is duplicated here.
use super::{
    Compressor, InputStreamBuffer, IoError, IoResult, StreamEntity, ZlibSource, ZstdSource,
};
use std::io::{self, Read};
use std::path::Path;

pub struct StdReadStream<R> {
    reader: R,
    label: String,
    eof: bool,
}
impl<R: Read> StdReadStream<R> {
    pub fn new(reader: R, label: impl Into<String>) -> Self {
        Self {
            reader,
            label: label.into(),
            eof: false,
        }
    }
}
impl<R: Read> StreamEntity for StdReadStream<R> {
    fn read(&mut self, buffer: &mut [u8]) -> IoResult<usize> {
        let mut used = 0;
        while used < buffer.len() {
            match self.reader.read(&mut buffer[used..]) {
                Ok(0) => {
                    self.eof = true;
                    break;
                }
                Ok(count) => used += count,
                Err(error) => {
                    // An actual Source InputFile can reach this adapter after
                    // codec/format preparation. Preserve the errno and raw
                    // perror captured by fread instead of reducing it to text.
                    if let Some(saved) = super::cstdio_file::StdioError::find(&error) {
                        return Err(IoError::StdioDetails(saved.clone()));
                    }
                    return Err(if self.label.is_empty() {
                        IoError::Other(error.to_string())
                    } else {
                        IoError::FileRead(self.label.clone())
                    });
                }
            }
        }
        Ok(used)
    }
    fn eof(&mut self) -> IoResult<bool> {
        Ok(self.eof)
    }
    fn file_name(&self) -> &str {
        &self.label
    }
}

pub struct StreamRead<S: StreamEntity>(pub S);
impl<S: StreamEntity> Read for StreamRead<S> {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        self.0.read(bytes).map_err(io::Error::other)
    }
}

#[derive(Clone, Copy, Debug)]
pub enum DecodeMode {
    Plain,
    /// Library convenience supports existing translated low-level codecs.
    Auto {
        allow_zstd: bool,
    },
}
impl DecodeMode {
    pub const NATIVE_AUTO: Self = Self::Auto { allow_zstd: true };
}

/// Source InputFile skips compressor probing for stdin and, on Unix, all
/// nonregular files. Cargo `zstd` selects Source's optional WITH_ZSTD profile.
/// This policy does not change InputFile/ZstdSource's typed library APIs.
pub fn source_path_mode(path: &Path) -> io::Result<DecodeMode> {
    if path.as_os_str().is_empty() || path == Path::new("-") {
        return Ok(DecodeMode::Plain);
    }
    #[cfg(unix)]
    if !std::fs::metadata(path)?.is_file() {
        return Ok(DecodeMode::Plain);
    }
    Ok(DecodeMode::Auto {
        allow_zstd: cfg!(feature = "zstd"),
    })
}

type Raw<R> = InputStreamBuffer<StdReadStream<R>>;
pub enum DecodedRead<R: Read> {
    Plain(StreamRead<Raw<R>>),
    /// Source stdin/nonregular input skips the constructor's compressor peek.
    /// FASTA's subsequent format peek(1) is the first operation to fail.
    ZeroFormatPeek(StreamRead<Raw<R>>),
    Zlib(StreamRead<InputStreamBuffer<ZlibSource<StdReadStream<R>>>>),
    Zstd(StreamRead<InputStreamBuffer<ZstdSource<StdReadStream<R>>>>),
}
impl<R: Read> Read for DecodedRead<R> {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Plain(reader) => reader.read(bytes),
            Self::ZeroFormatPeek(_) => {
                if bytes.is_empty() {
                    Ok(0)
                } else {
                    Err(io::Error::other("Invalid peek"))
                }
            }
            Self::Zlib(reader) => reader.read(bytes),
            Self::Zstd(reader) => reader.read(bytes),
        }
    }
}

pub fn decode_reader<R: Read>(
    reader: R,
    capacity: usize,
    mode: DecodeMode,
    label: impl Into<String>,
) -> io::Result<DecodedRead<R>> {
    if capacity == 0 {
        return Err(io::Error::other("File buffer size must be positive."));
    }
    decode_source_format_reader(reader, capacity, mode, label)
}

/// Source sequence-reader construction followed by its format peek. In Plain
/// mode (stdin/pipes), failure is deferred until the first nonempty format read.
/// General InputFile raw reads intentionally do not use this format adapter.
pub fn decode_source_format_reader<R: Read>(
    reader: R,
    capacity: usize,
    mode: DecodeMode,
    label: impl Into<String>,
) -> io::Result<DecodedRead<R>> {
    Vec::<u8>::new().try_reserve_exact(capacity).map_err(|_|io::Error::other(
        "Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage."))?;
    let mut input =
        InputStreamBuffer::with_source_buffer_size(StdReadStream::new(reader, label), 0, capacity);
    let compressor = if let DecodeMode::Auto { .. } = mode {
        input.fetch().map_err(io::Error::other)?;
        // Source Deserializer::peek(4) rejects insufficient staging capacity
        // unless the underlying stream has already reached EOF.
        if input.data().len() < 4 && !input.eof().map_err(io::Error::other)? {
            return Err(io::Error::other("Invalid peek"));
        }
        if input.data().len() >= 4 {
            super::detect_compressor(input.data())
        } else {
            Compressor::None
        }
    } else {
        Compressor::None
    };
    match compressor {
        Compressor::None if capacity == 0 => Ok(DecodedRead::ZeroFormatPeek(StreamRead(input))),
        Compressor::None => Ok(DecodedRead::Plain(StreamRead(input))),
        Compressor::Zlib => Ok(DecodedRead::Zlib(StreamRead(
            InputStreamBuffer::with_buffer_size(
                ZlibSource::new(input).map_err(io::Error::other)?,
                0,
                capacity,
            ),
        ))),
        Compressor::Zstd => {
            if matches!(mode, DecodeMode::Auto { allow_zstd: false }) {
                return Err(io::Error::other(
                    "Executable was not compiled with ZStd support.",
                ));
            }
            Ok(DecodedRead::Zstd(StreamRead(
                InputStreamBuffer::with_buffer_size(
                    ZstdSource::new(input).map_err(io::Error::other)?,
                    0,
                    capacity,
                ),
            )))
        }
    }
}
