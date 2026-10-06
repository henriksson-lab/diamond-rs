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
    use super::super::{Compressor, InputFile, OutputFile, TempFileData, VecStream};
    use super::*;

    fn compressed(bytes: &[u8]) -> Vec<u8> {
        zstd::stream::encode_all(bytes, 0).unwrap()
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

        let input = InputStreamBuffer::new(VecStream::from_vec(CPP_FRAME.to_vec()), 0);
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
