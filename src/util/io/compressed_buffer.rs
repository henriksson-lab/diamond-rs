//! Streaming compressed byte buffer.
//!
//! This mirrors `diamond/src/util/io/compressed_buffer.{h,cpp}`. DIAMOND's
//! default build uses a zlib stream; the `zstd` feature selects Source's
//! `WITH_ZSTD` branch for these internal streams as well.

#[cfg(not(feature = "zstd"))]
use flate2::{Compress, Compression, FlushCompress, Status};

use super::{FilePrimitive, IoError, IoResult};
#[cfg(feature = "zstd")]
use zstd_pure_rs::compress::zstd_compress::ZSTD_CCtx_setPledgedSrcSize;
#[cfg(feature = "zstd")]
use zstd_pure_rs::prelude::{
    ERR_getErrorName, ERR_isError, ZSTD_CStream, ZSTD_compressStream, ZSTD_createCStream,
    ZSTD_endStream, ZSTD_initCStream,
};

/// An in-memory Source-profile compression stream that can be finalized and reused.
pub struct CompressedBuffer {
    #[cfg(not(feature = "zstd"))]
    encoder: Option<Compress>,
    #[cfg(feature = "zstd")]
    encoder: Option<Box<ZSTD_CStream>>,
    buf: Vec<u8>,
    size: usize,
    #[cfg(feature = "zstd")]
    input_seen: bool,
}

impl std::fmt::Debug for CompressedBuffer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CompressedBuffer")
            .field("closed", &self.encoder.is_none())
            .field("size", &self.size())
            .finish_non_exhaustive()
    }
}

impl Default for CompressedBuffer {
    fn default() -> Self {
        Self::new()
    }
}

impl CompressedBuffer {
    /// Growth quantum used by the C++ implementation.
    pub const BUF_SIZE: usize = 32_768;

    /// Construct an empty, active compression stream.
    pub fn new() -> Self {
        let mut out = Self {
            encoder: None,
            buf: Vec::with_capacity(Self::BUF_SIZE),
            size: 0,
            #[cfg(feature = "zstd")]
            input_seen: false,
        };
        out.clear();
        out
    }

    /// Add uncompressed bytes to the active stream.
    #[cfg(not(feature = "zstd"))]
    pub fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        let stream = self
            .encoder
            .as_mut()
            .ok_or_else(|| IoError::Other("CompressedBuffer stream is closed.".to_string()))?;
        // Source caps avail_in to uInt, then keeps calling Z_NO_FLUSH until
        // that input chunk is consumed, growing output in 32-KiB increments.
        for chunk in ptr.chunks(u32::MAX as usize) {
            let mut input_position = 0;
            while input_position < chunk.len() {
                if self.size == self.buf.len() {
                    self.buf.resize(self.buf.len() + Self::BUF_SIZE, 0);
                }
                let input_before = stream.total_in();
                let output_before = stream.total_out();
                stream
                    .compress(
                        &chunk[input_position..],
                        &mut self.buf[self.size..],
                        FlushCompress::None,
                    )
                    .map_err(|_| IoError::Other("deflate(Z_NO_FLUSH) failed".into()))?;
                input_position += (stream.total_in() - input_before) as usize;
                self.size += (stream.total_out() - output_before) as usize;
                if input_position < chunk.len() && self.size == self.buf.len() {
                    self.buf.resize(self.buf.len() + Self::BUF_SIZE, 0);
                }
            }
        }
        Ok(())
    }

    #[cfg(feature = "zstd")]
    pub fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        let stream = self
            .encoder
            .as_mut()
            .ok_or_else(|| IoError::Other("CompressedBuffer stream is closed.".to_string()))?;
        let mut input_position = 0;
        while input_position < ptr.len() {
            if self.size == self.buf.len() {
                self.buf.resize(self.buf.len() + Self::BUF_SIZE, 0);
            }
            let mut output_position = 0;
            let result = ZSTD_compressStream(
                stream,
                &mut self.buf[self.size..],
                &mut output_position,
                ptr,
                &mut input_position,
            );
            if ERR_isError(result) {
                return Err(IoError::Other(ERR_getErrorName(result).to_string()));
            }
            self.input_seen |= input_position != 0;
            self.size += output_position;
            if output_position == 0 && input_position < ptr.len() {
                self.buf.resize(self.buf.len() + Self::BUF_SIZE, 0);
            }
        }
        Ok(())
    }

    /// Write a primitive in native byte order, matching C++ object-byte output.
    pub fn write_value<T: FilePrimitive>(&mut self, x: T) -> IoResult<()> {
        self.write(&x.to_ne_bytes_vec())
    }

    /// Finish the stream and make all compressed bytes available through [`Self::data`].
    #[cfg(not(feature = "zstd"))]
    pub fn finish(&mut self) -> IoResult<()> {
        let Some(stream) = self.encoder.as_mut() else {
            return Ok(());
        };
        loop {
            if self.size == self.buf.len() {
                self.buf.resize(self.buf.len() + Self::BUF_SIZE, 0);
            }
            let before = stream.total_out();
            let status = stream
                .compress(&[], &mut self.buf[self.size..], FlushCompress::Finish)
                .map_err(|_| IoError::Other("deflate(Z_FINISH) failed".into()))?;
            self.size += (stream.total_out() - before) as usize;
            match status {
                Status::StreamEnd => break,
                Status::Ok => {}
                Status::BufError => return Err(IoError::Other("deflate(Z_FINISH) failed".into())),
            }
        }
        self.encoder.take();
        Ok(())
    }

    #[cfg(feature = "zstd")]
    pub fn finish(&mut self) -> IoResult<()> {
        let Some(stream) = self.encoder.as_mut() else {
            return Ok(());
        };
        // Upstream endStream performs first-call e_end initialization and
        // implicitly pledges zero for a never-fed stream. The pure backend's
        // endStream bypasses that generic initialization; restore its actual
        // parameter transition rather than rewriting emitted frame bytes.
        if !self.input_seen {
            let result = ZSTD_CCtx_setPledgedSrcSize(stream, 0);
            if ERR_isError(result) {
                return Err(IoError::Other(ERR_getErrorName(result).to_string()));
            }
        }
        loop {
            if self.size == self.buf.len() {
                self.buf.resize(self.buf.len() + Self::BUF_SIZE, 0);
            }
            let mut output_position = 0;
            let remaining =
                ZSTD_endStream(stream, &mut self.buf[self.size..], &mut output_position);
            if ERR_isError(remaining) {
                return Err(IoError::Other(ERR_getErrorName(remaining).to_string()));
            }
            self.size += output_position;
            if remaining == 0 {
                break;
            }
        }
        self.encoder.take();
        Ok(())
    }

    /// Discard the current contents and start a fresh Source-profile stream.
    #[cfg(not(feature = "zstd"))]
    pub fn clear(&mut self) {
        self.encoder = Some(Compress::new(Compression::default(), true));
        self.size = 0;
        if self.buf.is_empty() {
            self.buf.resize(Self::BUF_SIZE, 0);
        }
    }

    #[cfg(feature = "zstd")]
    pub fn clear(&mut self) {
        let mut stream = ZSTD_createCStream().expect("ZSTD_createCStream error");
        // Source passes level zero, selecting the library's default level.
        let result = ZSTD_initCStream(&mut stream, 0);
        assert!(!ERR_isError(result), "{}", ERR_getErrorName(result));
        self.encoder = Some(stream);
        self.size = 0;
        self.input_seen = false;
        if self.buf.is_empty() {
            self.buf.resize(Self::BUF_SIZE, 0);
        }
    }

    /// Return the compressed prefix produced so far.
    pub fn data(&self) -> &[u8] {
        &self.buf[..self.size]
    }

    /// Return the number of compressed bytes produced so far.
    pub fn size(&self) -> usize {
        self.data().len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::{decompress, detect_compressor, Compressor};

    fn profile_compressor() -> Compressor {
        if cfg!(feature = "zstd") {
            Compressor::Zstd
        } else {
            Compressor::Zlib
        }
    }

    #[test]
    #[cfg(not(feature = "zstd"))]
    fn empty_stream_has_canonical_zlib_bytes() {
        let mut buffer = CompressedBuffer::new();
        assert_eq!(buffer.size(), 0);
        buffer.finish().unwrap();

        assert_eq!(
            buffer.data(),
            &[0x78, 0x9c, 0x03, 0x00, 0x00, 0x00, 0x00, 0x01]
        );
        assert_eq!(detect_compressor(buffer.data()), Compressor::Zlib);
    }

    #[test]
    #[cfg(not(feature = "zstd"))]
    fn source_z_no_flush_exposes_header_on_first_byte() {
        let mut buffer = CompressedBuffer::new();
        buffer.write(b"A").unwrap();
        // Actual immutable Source object, libz1.2.11: this prefix is visible
        // immediately, before either finish or a full input block.
        assert_eq!(buffer.data(), &[0x78, 0x9c]);
        buffer.finish().unwrap();
        buffer.clear();
        buffer.write(b"second").unwrap();
        assert_eq!(buffer.data(), &[0x78, 0x9c]);
    }

    #[test]
    #[cfg(feature = "zstd")]
    fn source_with_zstd_internal_buffer_emits_zstd_frame() {
        let mut buffer = CompressedBuffer::new();
        assert_eq!(buffer.size(), 0);
        buffer.finish().unwrap();
        assert_eq!(detect_compressor(buffer.data()), Compressor::Zstd);
        // Unchanged Source wrapper/libzstd1.4.8 and actual libzstd1.5.7
        // endStream both initialize the first empty frame with known size0.
        assert_eq!(buffer.data(), &[0x28, 0xb5, 0x2f, 0xfd, 0x20, 0, 1, 0, 0]);
        assert_eq!(
            decompress(buffer.data(), &mut [0; 1], Compressor::Zstd).unwrap(),
            0
        );
    }

    #[test]
    fn preserves_native_binary_values_across_buffer_growth() {
        let payload: Vec<u8> = (0..CompressedBuffer::BUF_SIZE * 3 + 17)
            .map(|i| (i.wrapping_mul(31) & 0xff) as u8)
            .collect();
        let mut buffer = CompressedBuffer::new();
        buffer.write_value(0x1234_5678u32).unwrap();
        buffer.write_value(-12_345i16).unwrap();
        buffer.write(&payload).unwrap();
        buffer.finish().unwrap();

        let mut decoded = vec![0; 6 + payload.len()];
        assert_eq!(detect_compressor(buffer.data()), profile_compressor());
        let decoded_len = decompress(buffer.data(), &mut decoded, profile_compressor()).unwrap();
        let mut expected = Vec::with_capacity(decoded.len());
        expected.extend_from_slice(&0x1234_5678u32.to_ne_bytes());
        expected.extend_from_slice(&(-12_345i16).to_ne_bytes());
        expected.extend_from_slice(&payload);
        assert_eq!(decoded_len, expected.len());
        assert_eq!(decoded, expected);
    }

    #[test]
    fn closed_stream_rejects_writes_and_clear_reopens_it() {
        let mut buffer = CompressedBuffer::new();
        buffer.write(b"first").unwrap();
        buffer.finish().unwrap();
        let finished = buffer.data().to_vec();

        assert_eq!(
            buffer.write(b"late"),
            Err(IoError::Other(
                "CompressedBuffer stream is closed.".to_string()
            ))
        );
        assert_eq!(buffer.data(), finished);

        buffer.clear();
        assert_eq!(buffer.size(), 0);
        buffer.write(b"second").unwrap();
        buffer.finish().unwrap();
        let mut decoded = [0; 6];
        assert_eq!(
            decompress(buffer.data(), &mut decoded, profile_compressor()).unwrap(),
            6
        );
        assert_eq!(&decoded, b"second");
    }

    #[test]
    fn streaming_prefix_and_final_output_survive_multiple_output_growths() {
        let mut state = 0x1234_5678u32;
        let payload: Vec<u8> = (0..CompressedBuffer::BUF_SIZE * 17)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                state as u8
            })
            .collect();
        let mut buffer = CompressedBuffer::new();
        for chunk in payload.chunks(17_003) {
            buffer.write(chunk).unwrap();
        }
        assert!(buffer.size() > CompressedBuffer::BUF_SIZE);
        let prefix = buffer.data().to_vec();
        buffer.finish().unwrap();
        assert!(buffer.data().starts_with(&prefix));
        assert_eq!(detect_compressor(buffer.data()), profile_compressor());
        let mut decoded = vec![0; payload.len()];
        assert_eq!(
            decompress(buffer.data(), &mut decoded, profile_compressor()).unwrap(),
            payload.len()
        );
        assert_eq!(decoded, payload);
    }
}
