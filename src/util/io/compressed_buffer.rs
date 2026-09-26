//! Streaming compressed byte buffer.
//!
//! This mirrors `diamond/src/util/io/compressed_buffer.{h,cpp}`. DIAMOND's
//! default build uses a zlib stream; its optional `WITH_ZSTD` compile-time
//! branch is represented elsewhere by the explicit zstd output facilities.

use std::io::Write;

use super::{FilePrimitive, IoError, IoResult};

/// An in-memory zlib stream that can be finalized and reused.
pub struct CompressedBuffer {
    encoder: Option<flate2::write::ZlibEncoder<Vec<u8>>>,
    buf: Vec<u8>,
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
        };
        out.clear();
        out
    }

    /// Add uncompressed bytes to the active stream.
    pub fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.encoder
            .as_mut()
            .ok_or_else(|| IoError::Other("CompressedBuffer stream is closed.".to_string()))?
            .write_all(ptr)
            .map_err(|e| IoError::Other(e.to_string()))
    }

    /// Write a primitive in native byte order, matching C++ object-byte output.
    pub fn write_value<T: FilePrimitive>(&mut self, x: T) -> IoResult<()> {
        self.write(&x.to_ne_bytes_vec())
    }

    /// Finish the stream and make all compressed bytes available through [`Self::data`].
    pub fn finish(&mut self) -> IoResult<()> {
        if let Some(encoder) = self.encoder.take() {
            self.buf = encoder
                .finish()
                .map_err(|e| IoError::Other(e.to_string()))?;
        }
        Ok(())
    }

    /// Discard the current contents and start a fresh zlib stream.
    pub fn clear(&mut self) {
        self.buf.clear();
        self.encoder = Some(flate2::write::ZlibEncoder::new(
            Vec::with_capacity(Self::BUF_SIZE),
            flate2::Compression::default(),
        ));
    }

    /// Return the compressed prefix produced so far.
    pub fn data(&self) -> &[u8] {
        if let Some(encoder) = &self.encoder {
            encoder.get_ref()
        } else {
            &self.buf
        }
    }

    /// Return the number of compressed bytes produced so far.
    pub fn size(&self) -> usize {
        self.data().len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::{detect_compressor, zlib_decompress, Compressor};

    #[test]
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
        let decoded_len = zlib_decompress(buffer.data(), &mut decoded).unwrap();
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
        assert_eq!(zlib_decompress(buffer.data(), &mut decoded).unwrap(), 6);
        assert_eq!(&decoded, b"second");
    }
}
