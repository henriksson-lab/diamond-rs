//! Translation of `diamond/src/util/io/output_stream_buffer.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::SeekFrom;

use super::{IoError, IoResult, StreamEntity, DEFAULT_FILE_BUFFER_SIZE};

#[derive(Debug, Clone)]
pub struct OutputStreamBuffer<S: StreamEntity> {
    prev: S,
    buf_size: usize,
    buf: Vec<u8>,
}

impl<S: StreamEntity> OutputStreamBuffer<S> {
    pub const STDOUT_BUF_SIZE: usize = 4096;

    /// C++ constructor using the configured default file-buffer size.
    pub fn new(prev: S) -> Self {
        Self::with_buffer_size(prev, DEFAULT_FILE_BUFFER_SIZE)
    }

    /// Explicit-config form replacing C++ `config.file_buffer_size`.
    pub fn with_buffer_size(prev: S, file_buffer_size: usize) -> Self {
        let buf_size = if prev.file_name().is_empty() {
            Self::STDOUT_BUF_SIZE
        } else {
            file_buffer_size
        };
        assert!(buf_size != 0, "output stream buffer size must be nonzero");
        Self {
            prev,
            buf_size,
            buf: Vec::with_capacity(buf_size),
        }
    }

    /// C++ `write_buffer()` pair represented as a mutable slice and length.
    pub fn write_buffer(&mut self) -> (&mut [u8], usize) {
        if self.buf.len() != self.buf_size {
            self.buf.resize(self.buf_size, 0);
        }
        (self.buf.as_mut_slice(), self.buf_size)
    }

    pub fn write_buffer_range(&mut self) -> (&mut [u8], usize) {
        self.write_buffer()
    }

    /// C++ `flush(size_t)`; kept distinct from the no-argument Rust stream
    /// flush because Rust does not support overloads.
    pub fn flush_count(&mut self, count: usize) -> IoResult<()> {
        if count > self.buf_size || count > self.buf.len() {
            return Err(IoError::Other(
                "OutputStreamBuffer::flush count exceeds buffer size.".to_string(),
            ));
        }
        self.prev.write(&self.buf[..count])?;
        self.buf.clear();
        Ok(())
    }

    pub fn buffer_size(&self) -> usize {
        self.buf_size
    }

    pub fn into_inner(mut self) -> IoResult<S> {
        self.flush()?;
        Ok(self.prev)
    }
}

impl<S: StreamEntity> StreamEntity for OutputStreamBuffer<S> {
    /// Compatibility extension used by the safe Rust serializer. The C++
    /// serializer writes through `write_buffer` and `flush(count)` directly.
    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        let mut offset = 0;
        while offset < ptr.len() {
            let available = self.buf_size - self.buf.len();
            let count = available.min(ptr.len() - offset);
            self.buf.extend_from_slice(&ptr[offset..offset + count]);
            offset += count;
            if self.buf.len() == self.buf_size {
                self.flush()?;
            }
        }
        Ok(())
    }

    fn flush(&mut self) -> IoResult<()> {
        if !self.buf.is_empty() {
            self.prev.write(&self.buf)?;
            self.buf.clear();
        }
        self.prev.flush()
    }

    fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        self.flush()?;
        self.prev.seek(position, origin)
    }

    fn rewind(&mut self) -> IoResult<()> {
        self.flush()?;
        self.prev.rewind()
    }

    fn tell(&mut self) -> IoResult<i64> {
        self.flush()?;
        self.prev.tell()
    }

    fn close(&mut self) -> IoResult<()> {
        self.flush()?;
        self.prev.close()
    }

    fn file_name(&self) -> &str {
        self.prev.file_name()
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.flush()?;
        self.prev.file()
    }

    fn file_size(&mut self) -> IoResult<i64> {
        self.flush()?;
        self.prev.file_size()
    }

    // `OutputStreamBuffer(StreamEntity*)` does not opt into seekable state.
    fn seekable(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Clone, Default)]
    struct RecordingStream {
        data: Vec<u8>,
        position: usize,
        flushes: usize,
        closes: usize,
        name: String,
        fail_write: bool,
    }

    impl StreamEntity for RecordingStream {
        fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
            if self.fail_write {
                return Err(IoError::FileWrite(self.name.clone()));
            }
            if self.position + ptr.len() > self.data.len() {
                self.data.resize(self.position + ptr.len(), 0);
            }
            self.data[self.position..self.position + ptr.len()].copy_from_slice(ptr);
            self.position += ptr.len();
            Ok(())
        }

        fn flush(&mut self) -> IoResult<()> {
            self.flushes += 1;
            Ok(())
        }

        fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
            let target = match origin {
                SeekFrom::Start(_) => position,
                SeekFrom::End(_) => self.data.len() as i64 + position,
                SeekFrom::Current(_) => self.position as i64 + position,
            };
            if target < 0 {
                return Err(IoError::UnsupportedOperation);
            }
            self.position = target as usize;
            Ok(())
        }

        fn rewind(&mut self) -> IoResult<()> {
            self.position = 0;
            Ok(())
        }

        fn tell(&mut self) -> IoResult<i64> {
            Ok(self.position as i64)
        }

        fn close(&mut self) -> IoResult<()> {
            self.closes += 1;
            Ok(())
        }

        fn file_name(&self) -> &str {
            &self.name
        }

        fn file_size(&mut self) -> IoResult<i64> {
            Ok(self.data.len() as i64)
        }

        fn seekable(&self) -> bool {
            true
        }
    }

    #[test]
    fn explicit_buffer_flushes_at_boundaries_and_preserves_order() {
        let stream = RecordingStream {
            name: "named".to_string(),
            ..RecordingStream::default()
        };
        let mut buffer = OutputStreamBuffer::with_buffer_size(stream, 4);
        assert_eq!(buffer.buffer_size(), 4);
        assert!(!buffer.seekable());
        buffer.write(b"abc").unwrap();
        assert!(buffer.prev.data.is_empty());
        buffer.write(b"defghi").unwrap();
        assert_eq!(buffer.prev.data, b"abcdefgh");
        buffer.flush().unwrap();
        assert_eq!(buffer.prev.data, b"abcdefghi");
    }

    #[test]
    fn write_buffer_and_flush_count_write_only_requested_prefix() {
        let stream = RecordingStream {
            name: "named".to_string(),
            ..RecordingStream::default()
        };
        let mut buffer = OutputStreamBuffer::with_buffer_size(stream, 8);
        let (bytes, size) = buffer.write_buffer();
        assert_eq!(size, 8);
        bytes.copy_from_slice(b"abcdefgh");
        buffer.flush_count(3).unwrap();
        assert_eq!(buffer.prev.data, b"abc");
        assert_eq!(
            buffer.flush_count(9).unwrap_err(),
            IoError::Other("OutputStreamBuffer::flush count exceeds buffer size.".to_string())
        );
    }

    #[test]
    fn stdout_name_forces_small_buffer() {
        let stdout = OutputStreamBuffer::with_buffer_size(RecordingStream::default(), 1 << 20);
        assert_eq!(
            stdout.buffer_size(),
            OutputStreamBuffer::<RecordingStream>::STDOUT_BUF_SIZE
        );
    }

    #[test]
    fn seek_rewind_tell_size_and_close_flush_pending_bytes_first() {
        let stream = RecordingStream {
            name: "named".to_string(),
            ..RecordingStream::default()
        };
        let mut buffer = OutputStreamBuffer::with_buffer_size(stream, 8);
        buffer.write(b"abcd").unwrap();
        assert_eq!(buffer.tell().unwrap(), 4);
        buffer.seek(1, SeekFrom::Start(0)).unwrap();
        buffer.write(b"Z").unwrap();
        assert_eq!(buffer.file_size().unwrap(), 4);
        buffer.rewind().unwrap();
        buffer.write(b"Q").unwrap();
        buffer.close().unwrap();
        assert_eq!(buffer.prev.data, b"QZcd");
        assert_eq!(buffer.prev.closes, 1);
    }

    #[test]
    fn write_and_close_propagate_predecessor_errors_without_losing_bytes() {
        let stream = RecordingStream {
            name: "failure".to_string(),
            fail_write: true,
            ..RecordingStream::default()
        };
        let mut buffer = OutputStreamBuffer::with_buffer_size(stream, 4);
        assert_eq!(
            buffer.write(b"four").unwrap_err(),
            IoError::FileWrite("failure".to_string())
        );
        assert_eq!(buffer.buf, b"four");
        assert_eq!(
            buffer.close().unwrap_err(),
            IoError::FileWrite("failure".to_string())
        );
        assert_eq!(buffer.prev.closes, 0);
    }
}
