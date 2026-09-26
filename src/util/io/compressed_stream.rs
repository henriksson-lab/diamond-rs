//! Translation of `diamond/src/util/io/compressed_stream.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::Read;

use flate2::{Compress, Compression, Decompress, FlushCompress, FlushDecompress, Status};

use super::{InputStreamBuffer, IoError, IoResult, OutputStreamBuffer, StreamEntity};

enum Decoder {
    Gzip(Decompress),
    Zlib(Decompress),
}

impl Decoder {
    fn detect(input: &[u8]) -> Self {
        if input.len() >= 2 && input[0] == 0x1f && input[1] == 0x8b {
            Self::Gzip(Decompress::new_gzip(15))
        } else {
            Self::Zlib(Decompress::new(true))
        }
    }

    fn stream(&mut self) -> &mut Decompress {
        match self {
            Self::Gzip(stream) | Self::Zlib(stream) => stream,
        }
    }
}

pub struct ZlibSource<S: StreamEntity> {
    prev: InputStreamBuffer<S>,
    stream: Option<Decoder>,
    eos: bool,
    closed: bool,
}

impl<S: StreamEntity> std::fmt::Debug for ZlibSource<S> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ZlibSource")
            .field("closed", &self.closed)
            .field("eos", &self.eos)
            .finish_non_exhaustive()
    }
}

impl<S: StreamEntity> ZlibSource<S> {
    pub fn new(prev: InputStreamBuffer<S>) -> IoResult<Self> {
        let mut source = Self {
            prev,
            stream: None,
            eos: false,
            closed: false,
        };
        source.init();
        Ok(source)
    }

    fn init(&mut self) {
        self.eos = false;
        self.closed = false;
        self.stream = None;
    }
}

impl<S: StreamEntity> StreamEntity for ZlibSource<S> {
    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        if self.closed {
            return Err(IoError::Other("ZlibSource closed".to_string()));
        }
        let mut output_position = 0;
        while output_position < ptr.len() && !self.eos {
            if self.prev.begin == self.prev.end && !self.prev.fetch()? {
                self.eos = true;
                break;
            }
            if self.stream.is_none() {
                self.stream = Some(Decoder::detect(
                    &self.prev.buf[self.prev.begin..self.prev.end],
                ));
            }

            let stream = self.stream.as_mut().unwrap().stream();
            let input_before = stream.total_in();
            let output_before_total = stream.total_out();
            let status = stream
                .decompress(
                    &self.prev.buf[self.prev.begin..self.prev.end],
                    &mut ptr[output_position..],
                    FlushDecompress::None,
                )
                .map_err(|_| {
                    IoError::Other(format!(
                        "Error reading gzip-compressed input file. The file may be corrupted: {}",
                        self.prev.file_name()
                    ))
                })?;
            let consumed = (stream.total_in() - input_before) as usize;
            let produced = (stream.total_out() - output_before_total) as usize;
            self.prev.begin += consumed;
            output_position += produced;

            if status == Status::StreamEnd {
                self.stream = None;
            } else if consumed == 0 && produced == 0 {
                return Err(IoError::Other(format!(
                    "Error reading gzip-compressed input file. The file may be corrupted: {}",
                    self.prev.file_name()
                )));
            }
        }
        Ok(output_position)
    }

    fn close(&mut self) -> IoResult<()> {
        // `inflateEnd` is harmless even between concatenated members; use eos
        // to distinguish a closed source from a member boundary.
        if self.closed {
            return Ok(());
        }
        self.stream = None;
        self.eos = true;
        self.closed = true;
        self.prev.close()
    }

    fn rewind(&mut self) -> IoResult<()> {
        self.prev.rewind()?;
        self.init();
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

pub struct ZlibSink<S: StreamEntity> {
    prev: OutputStreamBuffer<S>,
    stream: Option<Compress>,
}

impl<S: StreamEntity> std::fmt::Debug for ZlibSink<S> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ZlibSink")
            .field("closed", &self.stream.is_none())
            .finish_non_exhaustive()
    }
}

impl<S: StreamEntity> ZlibSink<S> {
    pub const CHUNK_SIZE: usize = 1 << 20;

    pub fn new(prev: OutputStreamBuffer<S>) -> Self {
        Self {
            prev,
            stream: Some(Compress::new_gzip(Compression::default(), 15)),
        }
    }

    pub fn deflate_loop(&mut self, ptr: &[u8], finish: bool) -> IoResult<()> {
        let stream = self
            .stream
            .as_mut()
            .ok_or_else(|| IoError::Other("deflate error".to_string()))?;
        let flush = if finish {
            FlushCompress::Finish
        } else {
            FlushCompress::None
        };
        let mut input_position = 0;
        loop {
            let input_before = stream.total_in();
            let output_before = stream.total_out();
            let (status, output_capacity) = {
                let (output, capacity) = self.prev.write_buffer_range();
                let status = stream
                    .compress(&ptr[input_position..], output, flush)
                    .map_err(|_| IoError::Other("deflate error".to_string()))?;
                (status, capacity)
            };
            let consumed = (stream.total_in() - input_before) as usize;
            let produced = (stream.total_out() - output_before) as usize;
            input_position += consumed;
            self.prev.flush_count(produced)?;

            if finish {
                if status == Status::StreamEnd {
                    break;
                }
            } else if input_position == ptr.len() && produced < output_capacity {
                break;
            }
            if consumed == 0 && produced == 0 {
                return Err(IoError::Other("deflate error".to_string()));
            }
        }
        Ok(())
    }
}

impl<S: StreamEntity> StreamEntity for ZlibSink<S> {
    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.deflate_loop(ptr, false)
    }

    fn close(&mut self) -> IoResult<()> {
        if self.stream.is_none() {
            return Ok(());
        }
        self.deflate_loop(&[], true)?;
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

pub fn zlib_decompress(src: &[u8], dst: &mut [u8]) -> IoResult<usize> {
    let mut input_position = 0;
    let mut output_position = 0;
    let mut ended_last_stream = false;

    while input_position < src.len() {
        let mut decoder = Decoder::detect(&src[input_position..]);
        loop {
            let stream = decoder.stream();
            let input_before = stream.total_in();
            let output_before = stream.total_out();
            let mut scratch = [0u8; 1];
            let output = if output_position < dst.len() {
                &mut dst[output_position..]
            } else {
                &mut scratch[..]
            };
            let status = stream
                .decompress(&src[input_position..], output, FlushDecompress::None)
                .map_err(|error| {
                    IoError::Other(format!("Error during zlib decompression: {error}"))
                })?;
            let consumed = (stream.total_in() - input_before) as usize;
            let produced = (stream.total_out() - output_before) as usize;
            input_position += consumed;
            if output_position == dst.len() && produced != 0 {
                return Err(IoError::Other(
                    "zlib_decompress: output buffer too small".to_string(),
                ));
            }
            output_position += produced;

            if status == Status::StreamEnd {
                ended_last_stream = true;
                break;
            }
            ended_last_stream = false;
            if consumed == 0 && produced == 0 {
                return Err(IoError::Other("Unexpected end of zlib stream".to_string()));
            }
            if input_position == src.len() {
                break;
            }
        }
    }
    if !ended_last_stream {
        return Err(IoError::Other("Unexpected end of zlib stream".to_string()));
    }
    Ok(output_position)
}

pub fn zlib_decompress_file(src: &mut StdFile, dst: &mut [u8]) -> IoResult<usize> {
    let mut compressed = Vec::new();
    src.read_to_end(&mut compressed)
        .map_err(|error| IoError::Other(format!("Error reading file: {error}")))?;
    zlib_decompress(&compressed, dst)
}

#[cfg(test)]
mod tests {
    use super::super::{FileSink, VecStream};
    use super::*;
    use std::io::Write;

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    fn zlib(bytes: &[u8]) -> Vec<u8> {
        let mut encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn source_detects_gzip_and_zlib_and_reads_concatenated_members() {
        for encoded in [gzip(b"gzip"), zlib(b"zlib")] {
            let input = InputStreamBuffer::new(VecStream::from_vec(encoded), 0);
            let mut source = ZlibSource::new(input).unwrap();
            let mut decoded = [0; 4];
            assert_eq!(source.read(&mut decoded).unwrap(), 4);
            assert!(!source.eof().unwrap());
            assert_eq!(source.read(&mut [0]).unwrap(), 0);
            assert!(source.eof().unwrap());
        }

        let encoded = [gzip(b"one").as_slice(), zlib(b"two").as_slice()].concat();
        let input = InputStreamBuffer::new(VecStream::from_vec(encoded), 0);
        let mut source = ZlibSource::new(input).unwrap();
        let mut decoded = [0; 6];
        assert_eq!(source.read(&mut decoded).unwrap(), 6);
        assert_eq!(&decoded, b"onetwo");
    }

    #[test]
    fn streaming_source_allows_truncated_input_but_whole_file_rejects_it() {
        let mut encoded = gzip(b"truncated-data");
        encoded.truncate(encoded.len() - 4);
        let input = InputStreamBuffer::new(VecStream::from_vec(encoded.clone()), 0);
        let mut source = ZlibSource::new(input).unwrap();
        let mut decoded = [0; 64];
        assert!(source.read(&mut decoded).is_ok());
        assert!(source.eof().unwrap());
        assert!(zlib_decompress(&encoded, &mut decoded).is_err());
    }

    #[test]
    fn source_rewind_resets_decoder() {
        let input = InputStreamBuffer::new(VecStream::from_vec(gzip(b"rewind")), 0);
        let mut source = ZlibSource::new(input).unwrap();
        let mut first = [0; 6];
        source.read(&mut first).unwrap();
        source.rewind().unwrap();
        let mut second = [0; 6];
        source.read(&mut second).unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn sink_streams_large_boundaries_and_finalizes_once() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-zlib-stream-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
        ));
        let _ = std::fs::remove_file(&path);
        let name = path.to_string_lossy();
        let file = FileSink::new(&name, "wb", false, 0).unwrap();
        let mut sink = ZlibSink::new(OutputStreamBuffer::new(file));
        let first = vec![b'a'; (1 << 20) + 17];
        sink.write(&first).unwrap();
        sink.write(b"tail").unwrap();
        sink.close().unwrap();
        sink.close().unwrap();
        assert_eq!(
            flate2::read::GzDecoder::new(std::fs::File::open(&path).unwrap())
                .bytes()
                .collect::<Result<Vec<_>, _>>()
                .unwrap(),
            [first, b"tail".to_vec()].concat()
        );
        assert!(sink.write(b"closed").is_err());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn whole_file_handles_concatenation_capacity_and_corruption() {
        let encoded = [gzip(b"abc").as_slice(), zlib(b"def").as_slice()].concat();
        let mut decoded = [0; 6];
        assert_eq!(zlib_decompress(&encoded, &mut decoded).unwrap(), 6);
        assert_eq!(&decoded, b"abcdef");
        assert!(zlib_decompress(&encoded, &mut [0; 5]).is_err());
        assert!(zlib_decompress(b"not compressed", &mut decoded).is_err());
    }
}
