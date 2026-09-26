//! Translation of `diamond/src/util/io/serializer.{h,cpp}`.

use std::fs::File;
use std::io::SeekFrom;

use super::{Consumer, FilePrimitive, IoResult, StreamEntity, DEFAULT_FILE_BUFFER_SIZE};

#[derive(Debug, Clone)]
pub struct Serializer<S: StreamEntity> {
    buffer: S,
    pending: Vec<u8>,
    pending_pos: usize,
}

impl<S: StreamEntity> Serializer<S> {
    /// C++ `Serializer::Serializer(StreamEntity*)`, including initial buffer
    /// acquisition through `reset_buffer()`.
    pub fn new(buffer: S) -> Self {
        Self {
            buffer,
            pending: vec![0; DEFAULT_FILE_BUFFER_SIZE],
            pending_pos: 0,
        }
    }

    pub fn write_i32(&mut self, value: i32) -> IoResult<&mut Self> {
        self.write_raw(&value.to_le_bytes())?;
        Ok(self)
    }

    /// Compatibility extension corresponding to the deserializer's C++
    /// `operator>>(short&)` little-endian representation.
    pub fn write_i16(&mut self, value: i16) -> IoResult<&mut Self> {
        self.write_raw(&value.to_le_bytes())?;
        Ok(self)
    }

    pub fn write_u16(&mut self, value: u16) -> IoResult<&mut Self> {
        self.write_raw(&value.to_le_bytes())?;
        Ok(self)
    }

    pub fn write_i64(&mut self, value: i64) -> IoResult<&mut Self> {
        self.write_raw(&value.to_le_bytes())?;
        Ok(self)
    }

    pub fn write_u32(&mut self, value: u32) -> IoResult<&mut Self> {
        self.write_raw(&value.to_le_bytes())?;
        Ok(self)
    }

    pub fn write_u64(&mut self, value: u64) -> IoResult<&mut Self> {
        self.write_raw(&value.to_le_bytes())?;
        Ok(self)
    }

    /// C++ deliberately writes double without `big_endian_byteswap`.
    pub fn write_f64(&mut self, value: f64) -> IoResult<&mut Self> {
        self.write_raw(&value.to_ne_bytes())?;
        Ok(self)
    }

    /// C++ template `write(const T&)` (raw native representation).
    pub fn write_value<T: FilePrimitive>(&mut self, value: T) -> IoResult<&mut Self> {
        self.write_raw(&value.to_ne_bytes_vec())?;
        Ok(self)
    }

    /// C++ template `write(const T*, size_t)`.
    pub fn write_values<T: FilePrimitive>(&mut self, values: &[T]) -> IoResult<&mut Self> {
        for &value in values {
            self.write_raw(&value.to_ne_bytes_vec())?;
        }
        Ok(self)
    }

    pub fn write_string(&mut self, value: &str) -> IoResult<&mut Self> {
        self.write_raw(value.as_bytes())?;
        self.write_raw(&[0])?;
        Ok(self)
    }

    pub fn write_raw(&mut self, ptr: &[u8]) -> IoResult<()> {
        let mut offset = 0;
        while offset < ptr.len() {
            if self.pending.is_empty() || self.pending_pos == self.pending.len() {
                self.reset_buffer();
            }
            let count = (self.pending.len() - self.pending_pos).min(ptr.len() - offset);
            self.pending[self.pending_pos..self.pending_pos + count]
                .copy_from_slice(&ptr[offset..offset + count]);
            self.pending_pos += count;
            offset += count;
            if self.pending_pos == self.pending.len() {
                self.flush()?;
                self.reset_buffer();
            }
        }
        Ok(())
    }

    /// Matches the header inline `file_size`: pending serializer bytes are not
    /// flushed implicitly.
    pub fn file_size(&mut self) -> IoResult<i64> {
        self.buffer.file_size()
    }

    pub fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        self.flush()?;
        self.buffer.seek(position, origin)?;
        self.reset_buffer();
        Ok(())
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.flush()?;
        self.buffer.rewind()?;
        self.reset_buffer();
        Ok(())
    }

    pub fn tell(&mut self) -> IoResult<usize> {
        self.flush()?;
        self.reset_buffer();
        self.buffer.tell().map(|position| position as usize)
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.flush()?;
        self.buffer.close()
    }

    pub fn file_name(&self) -> String {
        self.buffer.file_name().to_string()
    }

    /// Matches C++ `file()` without an implicit serializer flush.
    pub fn file(&mut self) -> IoResult<&mut File> {
        self.buffer.file()
    }

    pub fn consume(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.write_raw(ptr)
    }

    pub fn finalize(&mut self) -> IoResult<()> {
        self.close()
    }

    pub fn flush(&mut self) -> IoResult<()> {
        if self.pending_pos != 0 {
            self.buffer.write(&self.pending[..self.pending_pos])?;
            self.pending_pos = 0;
        }
        self.buffer.flush()
    }

    pub fn reset_buffer(&mut self) {
        if self.pending.is_empty() {
            self.pending = vec![0; DEFAULT_FILE_BUFFER_SIZE];
        }
        self.pending_pos = 0;
    }

    pub fn avail(&self) -> usize {
        self.pending.len() - self.pending_pos
    }

    pub fn into_inner(mut self) -> IoResult<S> {
        self.flush()?;
        Ok(self.buffer)
    }
}

impl<S: StreamEntity> Consumer for Serializer<S> {
    fn consume(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.write_raw(ptr)
    }

    fn finalize(&mut self) -> IoResult<()> {
        self.close()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::{IoError, VecStream};

    #[test]
    fn scalar_operators_emit_exact_cpp_bytes() {
        let mut serializer = Serializer::new(VecStream::new());
        serializer.write_i32(-7).unwrap();
        serializer.write_i64(-9).unwrap();
        serializer.write_u32(0x0102_0304).unwrap();
        serializer.write_u64(0x0102_0304_0506_0708).unwrap();
        serializer.write_f64(1.5).unwrap();
        serializer.write_string("abc").unwrap();
        let stream = serializer.into_inner().unwrap();

        let mut expected = Vec::new();
        expected.extend_from_slice(&(-7i32).to_le_bytes());
        expected.extend_from_slice(&(-9i64).to_le_bytes());
        expected.extend_from_slice(&0x0102_0304u32.to_le_bytes());
        expected.extend_from_slice(&0x0102_0304_0506_0708u64.to_le_bytes());
        expected.extend_from_slice(&1.5f64.to_ne_bytes());
        expected.extend_from_slice(b"abc\0");
        assert_eq!(stream.data(), expected);
    }

    #[test]
    fn file_size_excludes_pending_bytes_until_flush() {
        let mut serializer = Serializer::new(VecStream::new());
        serializer.write_raw(b"pending").unwrap();
        assert_eq!(serializer.file_size().unwrap(), 0);
        serializer.flush().unwrap();
        assert_eq!(serializer.file_size().unwrap(), 7);
    }

    #[test]
    fn seek_tell_rewind_and_buffer_boundary_preserve_order() {
        let mut serializer = Serializer::new(VecStream::new());
        let bytes = vec![b'a'; DEFAULT_FILE_BUFFER_SIZE + 3];
        serializer.write_raw(&bytes).unwrap();
        assert_eq!(serializer.tell().unwrap(), bytes.len());
        serializer.seek(1, SeekFrom::Start(0)).unwrap();
        serializer.write_raw(b"BC").unwrap();
        serializer.rewind().unwrap();
        serializer.write_raw(b"Z").unwrap();
        let stream = serializer.into_inner().unwrap();
        assert_eq!(&stream.data()[..4], b"ZBCa");
        assert_eq!(stream.data().len(), bytes.len());
    }

    #[derive(Debug, Clone)]
    struct FailingSink;

    impl StreamEntity for FailingSink {
        fn write(&mut self, _ptr: &[u8]) -> IoResult<()> {
            Err(IoError::Other("write failed".to_string()))
        }

        fn flush(&mut self) -> IoResult<()> {
            Ok(())
        }
    }

    #[test]
    fn flush_and_finalize_propagate_write_errors() {
        let mut serializer = Serializer::new(FailingSink);
        serializer.consume(b"x").unwrap();
        assert_eq!(
            serializer.flush().unwrap_err(),
            IoError::Other("write failed".to_string())
        );

        let mut serializer = Serializer::new(FailingSink);
        serializer.consume(b"x").unwrap();
        assert_eq!(
            serializer.finalize().unwrap_err(),
            IoError::Other("write failed".to_string())
        );
    }
}
