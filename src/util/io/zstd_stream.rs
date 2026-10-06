//! Translation of `diamond/src/util/io/zstd_stream.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::Read;

use zstd_pure_rs::decompress::zstd_decompress::ZSTD_decompressContinue_into_history;
use zstd_pure_rs::prelude::{
    ERR_getErrorName, ERR_isError, ZSTD_CStream, ZSTD_DStream, ZSTD_compressStream,
    ZSTD_createCStream, ZSTD_createDStream, ZSTD_decompress, ZSTD_endStream, ZSTD_initCStream,
    ZSTD_initDStream, ZSTD_nextSrcSizeToDecompress, ZSTD_resetDStream, ZSTD_CLEVEL_DEFAULT,
};

use super::{InputStreamBuffer, IoError, IoResult, OutputStreamBuffer, StreamEntity};

pub struct ZstdSink<S: StreamEntity> {
    prev: OutputStreamBuffer<S>,
    stream: Option<Box<ZSTD_CStream>>,
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
        let mut stream = ZSTD_createCStream().expect("ZSTD_createCStream error");
        // `ZSTD_createCStream` uses the default compression level.
        let result = ZSTD_initCStream(&mut stream, ZSTD_CLEVEL_DEFAULT);
        assert!(!ERR_isError(result), "ZSTD_initCStream error");
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
        let mut input_position = 0;
        while input_position < ptr.len() {
            let input_before = input_position;
            let output_position = {
                let (buffer, _) = self.prev.write_buffer_range();
                let mut output_position = 0;
                let result = ZSTD_compressStream(
                    stream,
                    buffer,
                    &mut output_position,
                    ptr,
                    &mut input_position,
                );
                if ERR_isError(result) {
                    return Err(IoError::Other(format!(
                        "ZSTD_compressStream: {}",
                        ERR_getErrorName(result)
                    )));
                }
                output_position
            };
            self.prev.flush_count(output_position)?;
            if input_position == input_before && output_position == 0 {
                return Err(IoError::Other(
                    "ZSTD_compressStream: no progress".to_string(),
                ));
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
                let mut output_position = 0;
                let remaining = ZSTD_endStream(stream, buffer, &mut output_position);
                if ERR_isError(remaining) {
                    return Err(IoError::Other(format!(
                        "ZSTD_endStream: {}",
                        ERR_getErrorName(remaining)
                    )));
                }
                (remaining, output_position)
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
    stream: Option<Box<ZSTD_DStream>>,
    eos: bool,
    compressed_chunk: Vec<u8>,
    pending_begin: usize,
    pending_end: usize,
}

impl<S: StreamEntity> std::fmt::Debug for ZstdSource<S> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ZstdSource")
            .field("closed", &self.stream.is_none())
            .field("eos", &self.eos)
            .field("pending", &(self.pending_end - self.pending_begin))
            .finish_non_exhaustive()
    }
}

impl<S: StreamEntity> ZstdSource<S> {
    pub fn new(prev: InputStreamBuffer<S>) -> IoResult<Self> {
        Ok(Self {
            prev,
            stream: Some(Self::init()?),
            eos: false,
            compressed_chunk: Vec::new(),
            pending_begin: 0,
            pending_end: 0,
        })
    }

    fn init() -> IoResult<Box<ZSTD_DStream>> {
        let mut stream =
            ZSTD_createDStream().ok_or_else(|| IoError::Other("ZSTD_createDStream".to_string()))?;
        let result = ZSTD_initDStream(&mut stream);
        if ERR_isError(result) {
            return Err(IoError::Other(format!(
                "ZSTD_initDStream: {}",
                ERR_getErrorName(result)
            )));
        }
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
            if self.pending_begin < self.pending_end {
                let count =
                    (ptr.len() - output_position).min(self.pending_end - self.pending_begin);
                ptr[output_position..output_position + count].copy_from_slice(
                    &stream.historyBuffer[self.pending_begin..self.pending_begin + count],
                );
                self.pending_begin += count;
                output_position += count;
                continue;
            }

            let expected = ZSTD_nextSrcSizeToDecompress(stream);
            if expected == 0 {
                if self.prev.begin == self.prev.end && !self.prev.fetch()? {
                    self.eos = true;
                    break;
                }
                let result = ZSTD_resetDStream(stream);
                if ERR_isError(result) {
                    return Err(IoError::Other(format!(
                        "ZSTD_resetDStream: {}",
                        ERR_getErrorName(result)
                    )));
                }
                continue;
            }

            self.compressed_chunk.clear();
            while self.compressed_chunk.len() < expected {
                if self.prev.begin == self.prev.end && !self.prev.fetch()? {
                    self.eos = true;
                    return Ok(output_position);
                }
                let count =
                    (expected - self.compressed_chunk.len()).min(self.prev.end - self.prev.begin);
                self.compressed_chunk
                    .extend_from_slice(&self.prev.buf[self.prev.begin..self.prev.begin + count]);
                self.prev.begin += count;
            }

            let produced_len = ZSTD_decompressContinue_into_history(stream, &self.compressed_chunk)
                .map_err(|error| {
                    IoError::Other(format!(
                        "ZSTD_decompressContinue: {}",
                        ERR_getErrorName(error)
                    ))
                })?
                .len();
            self.pending_end = stream.historyBuffer.len();
            self.pending_begin = self.pending_end - produced_len;
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
        self.compressed_chunk.clear();
        self.pending_begin = 0;
        self.pending_end = 0;
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
    let result = ZSTD_decompress(dst, src);
    if ERR_isError(result) {
        return Err(IoError::Other(format!(
            "Failed decompressing zstd stream: {}",
            ERR_getErrorName(result)
        )));
    }
    Ok(result)
}

#[cfg(test)]
pub(crate) fn zstd_compress_for_test(src: &[u8]) -> Vec<u8> {
    use zstd_pure_rs::prelude::{ZSTD_compress, ZSTD_compressBound};

    let mut dst = vec![0; ZSTD_compressBound(src.len())];
    let result = ZSTD_compress(&mut dst, src, ZSTD_CLEVEL_DEFAULT);
    assert!(!ERR_isError(result), "{}", ERR_getErrorName(result));
    dst.truncate(result);
    dst
}

pub fn zstd_decompress_file(src: &mut StdFile, dst: &mut [u8]) -> IoResult<usize> {
    let mut compressed = Vec::new();
    src.read_to_end(&mut compressed)
        .map_err(|error| IoError::Other(format!("Error reading file: {error}")))?;
    zstd_decompress(&compressed, dst)
}

#[cfg(test)]
mod tests {
    use super::super::{Compressor, InputFile, OutputFile, TempFileData, VecStream};
    use super::*;

    fn compressed(bytes: &[u8]) -> Vec<u8> {
        zstd_compress_for_test(bytes)
    }

    #[test]
    fn decodes_frame_emitted_by_upstream_cpp_stream() {
        // Produced by the checked-in `diamond/src/util/io/zstd_stream.cpp`
        // through ZstdSink, linked to libzstd 1.4.8. Keeping the encoded bytes
        // here makes C++ -> Rust compatibility a normal, toolchain-independent
        // test instead of relying only on frames produced by our own backend.
        const CPP_FRAME: &[u8] = &[
            0x28, 0xb5, 0x2f, 0xfd, 0x00, 0x58, 0xd1, 0x00, 0x00, 0x75, 0x70, 0x73, 0x74, 0x72,
            0x65, 0x61, 0x6d, 0x2d, 0x63, 0x70, 0x70, 0x2d, 0x7a, 0x73, 0x74, 0x64, 0x2d, 0x66,
            0x69, 0x78, 0x74, 0x75, 0x72, 0x65, 0x0a,
        ];
        const PLAIN: &[u8] = b"upstream-cpp-zstd-fixture\n";

        let mut whole = vec![0; PLAIN.len()];
        assert_eq!(zstd_decompress(CPP_FRAME, &mut whole).unwrap(), PLAIN.len());
        assert_eq!(whole, PLAIN);

        let input =
            InputStreamBuffer::with_buffer_size(VecStream::from_vec(CPP_FRAME.to_vec()), 0, 3);
        let mut source = ZstdSource::new(input).unwrap();
        let mut streamed = Vec::new();
        let mut chunk = [0; 7];
        loop {
            let count = source.read(&mut chunk).unwrap();
            if count == 0 {
                break;
            }
            streamed.extend_from_slice(&chunk[..count]);
        }
        assert_eq!(streamed, PLAIN);
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
    fn empty_frame_has_standard_magic_and_roundtrips() {
        let encoded = compressed(b"");
        // Zstandard frames use this little-endian magic in the upstream C++
        // implementation and in the Rust zstd backend.
        assert!(encoded.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]));

        let mut whole_output = [];
        assert_eq!(zstd_decompress(&encoded, &mut whole_output).unwrap(), 0);

        let input = InputStreamBuffer::new(VecStream::from_vec(encoded), 0);
        let mut source = ZstdSource::new(input).unwrap();
        let mut byte = [0; 1];
        assert_eq!(source.read(&mut byte).unwrap(), 0);
        assert!(source.eof().unwrap());
    }

    #[test]
    fn concatenated_frames_roundtrip_with_whole_and_streaming_decoders() {
        let expected = b"first-framesecond-frame";
        let mut encoded = compressed(b"first-frame");
        encoded.extend_from_slice(&compressed(b"second-frame"));

        let mut whole_output = vec![0; expected.len()];
        assert_eq!(
            zstd_decompress(&encoded, &mut whole_output).unwrap(),
            expected.len()
        );
        assert_eq!(whole_output, expected);

        let input = InputStreamBuffer::new(VecStream::from_vec(encoded), 0);
        let mut source = ZstdSource::new(input).unwrap();
        let mut streamed = Vec::new();
        let mut chunk = [0; 5];
        loop {
            let count = source.read(&mut chunk).unwrap();
            if count == 0 {
                break;
            }
            streamed.extend_from_slice(&chunk[..count]);
        }
        assert_eq!(streamed, expected);
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
        let stream = source.stream.as_ref().unwrap();
        assert!(stream.stream_in_buffer.is_empty());
        assert!(stream.stream_out_buffer.is_empty());
        assert!(source.compressed_chunk.capacity() < plain.len());
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
        let expected = [first, b"tail".to_vec()].concat();
        let mut decoded = vec![0; expected.len()];
        let count = zstd_decompress(&std::fs::read(&path).unwrap(), &mut decoded).unwrap();
        decoded.truncate(count);
        assert_eq!(decoded, expected);
        assert!(sink.write(b"closed").is_err());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn large_temp_file_roundtrips_through_output_and_input_streams() {
        let data = TempFileData::init(true).unwrap();
        let name = data.name.clone();
        let unlinked = data.unlinked;
        let plain = (0..(3 * (1 << 20) + 137))
            .map(|index| ((index * 31 + index / 251) & 0xff) as u8)
            .collect::<Vec<_>>();

        {
            let mut output =
                OutputFile::from_temp_file_data(&data, Compressor::Zstd, "w+b").unwrap();
            for chunk in plain.chunks(65_521) {
                output.write_raw(chunk).unwrap();
            }
            output.close().unwrap();
        }

        let mut input = InputFile::from_temp_file_data(&data, 0, Compressor::Zstd).unwrap();
        let mut decoded = Vec::with_capacity(plain.len());
        let mut chunk = vec![0; 37_117];
        loop {
            let count = input.read_raw(&mut chunk).unwrap();
            if count == 0 {
                break;
            }
            decoded.extend_from_slice(&chunk[..count]);
        }
        input.close().unwrap();
        assert_eq!(decoded, plain);

        drop(data);
        if !unlinked {
            std::fs::remove_file(name).unwrap();
        }
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
