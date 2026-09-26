//! Translation of `diamond/src/util/io/zstd_stream.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::Read;

use zstd::zstd_safe::{CCtx, DCtx, InBuffer, OutBuffer};

use super::{InputStreamBuffer, IoError, IoResult, OutputStreamBuffer, StreamEntity};

pub struct ZstdSink<S: StreamEntity> {
    prev: OutputStreamBuffer<S>,
    stream: Option<CCtx<'static>>,
}

impl<S: StreamEntity> std::fmt::Debug for ZstdSink<S> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ZstdSink")
            .field("closed", &self.stream.is_none())
            .finish_non_exhaustive()
    }
}

impl<S: StreamEntity> ZstdSink<S> {
    pub fn new(prev: OutputStreamBuffer<S>) -> Self {
        let mut stream = CCtx::create();
        // `ZSTD_createCStream` uses the default compression level.
        stream
            .init(zstd::DEFAULT_COMPRESSION_LEVEL)
            .expect("ZSTD_createCStream error");
        Self {
            prev,
            stream: Some(stream),
        }
    }
}

impl<S: StreamEntity> StreamEntity for ZstdSink<S> {
    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        let stream = self
            .stream
            .as_mut()
            .ok_or_else(|| IoError::Other("ZSTD_compressStream".to_string()))?;
        let mut input = InBuffer::around(ptr);
        loop {
            let output_position = {
                let (buffer, _) = self.prev.write_buffer_range();
                let mut output = OutBuffer::around(buffer);
                stream
                    .compress_stream(&mut output, &mut input)
                    .map_err(|_| IoError::Other("ZSTD_compressStream".to_string()))?;
                output.pos()
            };
            self.prev.flush_count(output_position)?;
            if input.pos() == ptr.len() {
                break;
            }
        }
        Ok(())
    }

    fn close(&mut self) -> IoResult<()> {
        let Some(stream) = self.stream.as_mut() else {
            return Ok(());
        };
        loop {
            let (remaining, output_position) = {
                let (buffer, _) = self.prev.write_buffer_range();
                let mut output = OutBuffer::around(buffer);
                let remaining = stream
                    .end_stream(&mut output)
                    .map_err(|_| IoError::Other("ZSTD_endStream".to_string()))?;
                (remaining, output.pos())
            };
            self.prev.flush_count(output_position)?;
            if remaining == 0 {
                break;
            }
        }
        // Match C++ ordering: release/clear the stream after a successful end,
        // but before closing the predecessor.
        self.stream.take();
        self.prev.close()
    }

    fn flush(&mut self) -> IoResult<()> {
        self.prev.flush()
    }

    fn file_name(&self) -> &str {
        self.prev.file_name()
    }

    fn file_size(&mut self) -> IoResult<i64> {
        self.prev.file_size()
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.prev.file()
    }
}

pub struct ZstdSource<S: StreamEntity> {
    prev: InputStreamBuffer<S>,
    stream: Option<DCtx<'static>>,
    eos: bool,
}

impl<S: StreamEntity> std::fmt::Debug for ZstdSource<S> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ZstdSource")
            .field("closed", &self.stream.is_none())
            .field("eos", &self.eos)
            .finish_non_exhaustive()
    }
}

impl<S: StreamEntity> ZstdSource<S> {
    pub fn new(prev: InputStreamBuffer<S>) -> IoResult<Self> {
        Ok(Self {
            prev,
            stream: Some(Self::init()?),
            eos: false,
        })
    }

    fn init() -> IoResult<DCtx<'static>> {
        let mut stream = DCtx::create();
        stream
            .init()
            .map_err(|error| IoError::Other(format!("ZSTD_initDStream: {error}")))?;
        Ok(stream)
    }
}

impl<S: StreamEntity> StreamEntity for ZstdSource<S> {
    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let stream = self
            .stream
            .as_mut()
            .ok_or_else(|| IoError::Other("ZstdSource closed".to_string()))?;
        let mut output_position = 0;

        while output_position < ptr.len() {
            if self.prev.begin == self.prev.end && !self.prev.fetch()? {
                self.eos = true;
                break;
            }

            let output_before = output_position;
            let consumed = {
                let input_slice = &self.prev.buf[self.prev.begin..self.prev.end];
                let mut input = InBuffer::around(input_slice);
                let mut output = OutBuffer::around_pos(ptr, output_position);
                stream
                    .decompress_stream(&mut output, &mut input)
                    .map_err(|error| IoError::Other(format!("ZSTD_decompressStream: {error}")))?;
                output_position = output.pos();
                input.pos()
            };
            self.prev.begin += consumed;

            if consumed == 0 && output_position == output_before && output_position < ptr.len() {
                return Err(IoError::Other(
                    "ZSTD_decompressStream: no progress".to_string(),
                ));
            }
        }
        Ok(output_position)
    }

    fn close(&mut self) -> IoResult<()> {
        if self.stream.take().is_none() {
            return Ok(());
        }
        self.prev.close()
    }

    fn rewind(&mut self) -> IoResult<()> {
        self.prev.rewind()?;
        self.stream = Some(Self::init()?);
        self.eos = false;
        Ok(())
    }

    fn eof(&mut self) -> IoResult<bool> {
        Ok(self.eos)
    }

    fn file_name(&self) -> &str {
        self.prev.file_name()
    }

    fn file_size(&mut self) -> IoResult<i64> {
        self.prev.file_size()
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.prev.file()
    }
}

pub fn zstd_decompress(src: &[u8], dst: &mut [u8]) -> IoResult<usize> {
    let mut stream = DCtx::create();
    stream
        .init()
        .map_err(|error| IoError::Other(format!("Failed decompressing zstd stream: {error}")))?;
    let mut input = InBuffer::around(src);
    let mut output = OutBuffer::around(dst);
    let mut last_return = 1;

    while input.pos() < src.len() {
        let input_before = input.pos();
        let output_before = output.pos();
        last_return = stream
            .decompress_stream(&mut output, &mut input)
            .map_err(|error| {
                IoError::Other(format!("Failed decompressing zstd stream: {error}"))
            })?;
        if input.pos() == input_before && output.pos() == output_before {
            return Err(IoError::Other(
                "Failed decompressing zstd stream: output buffer too small".to_string(),
            ));
        }
    }
    if last_return != 0 {
        return Err(IoError::Other(
            "Failed decompressing zstd stream".to_string(),
        ));
    }
    Ok(output.pos())
}

pub fn zstd_decompress_file(src: &mut StdFile, dst: &mut [u8]) -> IoResult<usize> {
    let mut compressed = Vec::new();
    src.read_to_end(&mut compressed)
        .map_err(|error| IoError::Other(format!("Error reading file: {error}")))?;
    zstd_decompress(&compressed, dst)
}

#[cfg(test)]
mod tests {
    use super::super::VecStream;
    use super::*;

    fn compressed(bytes: &[u8]) -> Vec<u8> {
        zstd::stream::encode_all(bytes, 0).unwrap()
    }

    #[test]
    fn source_reads_across_boundaries_and_observes_eof_lazily() {
        let plain = vec![b'x'; (1 << 20) + 37];
        let encoded = compressed(&plain);
        let input = InputStreamBuffer::new(VecStream::from_vec(encoded), 0);
        let mut source = ZstdSource::new(input).unwrap();
        let mut decoded = vec![0; plain.len()];
        assert_eq!(source.read(&mut decoded).unwrap(), plain.len());
        assert_eq!(decoded, plain);
        assert!(!source.eof().unwrap());
        assert_eq!(source.read(&mut [0]).unwrap(), 0);
        assert!(source.eof().unwrap());
    }

    #[test]
    fn source_matches_cpp_truncated_input_behavior() {
        let mut encoded = compressed(b"truncated-payload");
        encoded.truncate(encoded.len() - 3);
        let input = InputStreamBuffer::new(VecStream::from_vec(encoded), 0);
        let mut source = ZstdSource::new(input).unwrap();
        let mut decoded = [0; 64];
        let count = source.read(&mut decoded).unwrap();
        assert!(count <= b"truncated-payload".len());
        assert!(source.eof().unwrap());
    }

    #[test]
    fn source_rewind_reinitializes_decoder() {
        let input = InputStreamBuffer::new(VecStream::from_vec(compressed(b"rewind")), 0);
        let mut source = ZstdSource::new(input).unwrap();
        let mut first = [0; 6];
        source.read(&mut first).unwrap();
        source.rewind().unwrap();
        let mut second = [0; 6];
        source.read(&mut second).unwrap();
        assert_eq!(first, second);
        assert_eq!(&second, b"rewind");
    }

    #[test]
    fn source_accepts_decoder_output_without_new_input_consumption() {
        // Highly compressible input makes the decoder retain substantial
        // output internally; uneven small reads exercise calls that can emit
        // bytes without consuming additional compressed input.
        let plain = vec![b'q'; (2 << 20) + 123];
        let input = InputStreamBuffer::new(VecStream::from_vec(compressed(&plain)), 0);
        let mut source = ZstdSource::new(input).unwrap();
        let mut decoded = Vec::new();
        let mut chunk = [0; 997];
        loop {
            let count = source.read(&mut chunk).unwrap();
            if count == 0 {
                break;
            }
            decoded.extend_from_slice(&chunk[..count]);
        }
        assert_eq!(decoded, plain);
    }

    #[test]
    fn sink_streams_and_close_finishes_exact_frame() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-zstd-stream-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
        ));
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        let file = super::super::FileSink::new(&name, "wb", false, 0).unwrap();
        let mut sink = ZstdSink::new(OutputStreamBuffer::new(file));
        let first = vec![b'a'; (1 << 20) + 11];
        sink.write(&first).unwrap();
        sink.write(b"tail").unwrap();
        sink.close().unwrap();
        sink.close().unwrap();
        assert_eq!(
            zstd::stream::decode_all(std::fs::read(&path).unwrap().as_slice()).unwrap(),
            [first, b"tail".to_vec()].concat()
        );
        assert!(sink.write(b"closed").is_err());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn whole_file_decompress_rejects_truncation_and_small_output() {
        let encoded = compressed(b"complete");
        let mut output = [0; 8];
        assert_eq!(zstd_decompress(&encoded, &mut output).unwrap(), 8);
        assert_eq!(&output, b"complete");

        let mut truncated = encoded.clone();
        truncated.pop();
        assert!(zstd_decompress(&truncated, &mut output).is_err());
        assert!(zstd_decompress(&encoded, &mut [0; 2]).is_err());
    }
}
