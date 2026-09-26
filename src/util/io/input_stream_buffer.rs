//! Translation of `diamond/src/util/io/input_stream_buffer.{h,cpp}`.

use std::fs::File as StdFile;
use std::io::SeekFrom;

use super::{IoError, IoResult, StreamEntity, DEFAULT_FILE_BUFFER_SIZE};

#[derive(Debug, Clone)]
pub struct InputStreamBuffer<S: StreamEntity> {
    prev: S,
    buf_size: usize,
    pub(crate) buf: Vec<u8>,
    pub(crate) begin: usize,
    pub(crate) end: usize,
    load_buf: Option<Vec<u8>>,
    load_result: Option<IoResult<usize>>,
    file_offset: usize,
    async_: bool,
}

impl<S: StreamEntity> InputStreamBuffer<S> {
    pub const ASYNC: i32 = 4;

    /// C++ constructor using the configured default file-buffer size.
    pub fn new(prev: S, flags: i32) -> Self {
        Self::with_buffer_size(prev, flags, DEFAULT_FILE_BUFFER_SIZE)
    }

    /// Explicit-config form replacing C++ `config.file_buffer_size`.
    pub fn with_buffer_size(prev: S, flags: i32, file_buffer_size: usize) -> Self {
        assert!(
            file_buffer_size != 0,
            "input stream buffer size must be nonzero"
        );
        let async_ = flags & Self::ASYNC != 0;
        Self {
            prev,
            buf_size: file_buffer_size,
            buf: vec![0; file_buffer_size],
            begin: 0,
            end: 0,
            load_buf: async_.then(|| vec![0; file_buffer_size]),
            load_result: None,
            file_offset: 0,
            async_,
        }
    }

    /// Bytes in the public C++ `[begin, end)` range.
    pub fn data(&self) -> &[u8] {
        &self.buf[self.begin..self.end]
    }

    /// Safe equivalent of advancing C++'s public `begin` pointer.
    pub fn consume(&mut self, count: usize) -> IoResult<()> {
        if count > self.end - self.begin {
            return Err(IoError::EndOfStream);
        }
        self.begin += count;
        Ok(())
    }

    /// Compatibility helper used by Rust callers of the old embedded type.
    pub fn pop(&mut self, dst: &mut [u8]) -> usize {
        let count = dst.len().min(self.end - self.begin);
        dst[..count].copy_from_slice(&self.buf[self.begin..self.begin + count]);
        self.begin += count;
        count
    }

    fn load_worker(&mut self) {
        let load_buf = self
            .load_buf
            .as_mut()
            .expect("async input stream has a load buffer");
        // C++ performs this read on `load_worker_`. Keeping the completed read
        // pending preserves its double-buffer/read-ahead behavior without
        // imposing `Send + 'static` on every StreamEntity implementation.
        self.load_result = Some(self.prev.read(&mut load_buf[..self.buf_size]));
    }
}

impl<S: StreamEntity> StreamEntity for InputStreamBuffer<S> {
    fn rewind(&mut self) -> IoResult<()> {
        self.prev.rewind()?;
        self.file_offset = 0;
        self.begin = 0;
        self.end = 0;
        self.load_result = None;
        Ok(())
    }

    fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        self.prev.seek(position, origin)?;
        self.file_offset = 0;
        self.begin = 0;
        self.end = 0;
        self.load_result = None;
        Ok(())
    }

    fn fetch(&mut self) -> IoResult<bool> {
        if let Some(result) = self.load_result.take() {
            let count = result?;
            std::mem::swap(
                &mut self.buf,
                self.load_buf
                    .as_mut()
                    .expect("async input stream has a load buffer"),
            );
            self.begin = 0;
            self.end = count;
        } else {
            let count = self.prev.read(&mut self.buf[..self.buf_size])?;
            if self.prev.seekable() {
                self.file_offset = self.prev.tell()? as usize;
            }
            self.begin = 0;
            self.end = count;
        }

        if self.async_ {
            self.load_worker();
        }
        Ok(self.end > self.begin)
    }

    /// Compatibility extension: C++ consumers operate directly on
    /// `begin`/`end`, while generic Rust consumers use `read`.
    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let mut total = 0;
        while total < ptr.len() {
            if self.begin == self.end && !self.fetch()? {
                break;
            }
            total += self.pop(&mut ptr[total..]);
        }
        Ok(total)
    }

    fn close(&mut self) -> IoResult<()> {
        // The C++ close joins an outstanding worker but deliberately ignores
        // already-prefetched bytes.
        self.load_result = None;
        self.prev.close()
    }

    fn tell(&mut self) -> IoResult<i64> {
        if !self.seekable() {
            return Err(IoError::Other(
                "Calling tell on non seekable stream.".to_string(),
            ));
        }
        // Exact C++ behavior: this is the predecessor's position immediately
        // after the synchronous fetch, not `begin`'s logical position.
        Ok(self.file_offset as i64)
    }

    fn eof(&mut self) -> IoResult<bool> {
        self.prev.eof()
    }

    fn file_size(&mut self) -> IoResult<i64> {
        self.prev.file_size()
    }

    fn seekable(&self) -> bool {
        self.prev.seekable()
    }

    fn file_name(&self) -> &str {
        self.prev.file_name()
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.prev.file()
    }

    fn data(&self) -> &[u8] {
        self.data()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::VecStream;
    use std::cell::RefCell;
    use std::rc::Rc;

    #[derive(Clone)]
    struct ObservedStream {
        bytes: Vec<u8>,
        position: usize,
        reads: Rc<RefCell<Vec<usize>>>,
        fail_read: Option<usize>,
        closed: Rc<RefCell<bool>>,
        seekable: bool,
    }

    impl ObservedStream {
        fn new(bytes: &[u8], seekable: bool) -> Self {
            Self {
                bytes: bytes.to_vec(),
                position: 0,
                reads: Rc::new(RefCell::new(Vec::new())),
                fail_read: None,
                closed: Rc::new(RefCell::new(false)),
                seekable,
            }
        }
    }

    impl StreamEntity for ObservedStream {
        fn rewind(&mut self) -> IoResult<()> {
            self.position = 0;
            Ok(())
        }

        fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
            let next = match origin {
                SeekFrom::Start(_) => position,
                SeekFrom::Current(_) => self.position as i64 + position,
                SeekFrom::End(_) => self.bytes.len() as i64 + position,
            };
            if next < 0 {
                return Err(IoError::UnsupportedOperation);
            }
            self.position = next as usize;
            Ok(())
        }

        fn tell(&mut self) -> IoResult<i64> {
            Ok(self.position as i64)
        }

        fn read(&mut self, dst: &mut [u8]) -> IoResult<usize> {
            let call = self.reads.borrow().len();
            self.reads.borrow_mut().push(dst.len());
            if self.fail_read == Some(call) {
                return Err(IoError::FileRead("observed".to_string()));
            }
            let count = dst
                .len()
                .min(self.bytes.len().saturating_sub(self.position));
            dst[..count].copy_from_slice(&self.bytes[self.position..self.position + count]);
            self.position += count;
            Ok(count)
        }

        fn close(&mut self) -> IoResult<()> {
            *self.closed.borrow_mut() = true;
            Ok(())
        }

        fn eof(&mut self) -> IoResult<bool> {
            Ok(self.position == self.bytes.len())
        }

        fn file_size(&mut self) -> IoResult<i64> {
            Ok(self.bytes.len() as i64)
        }

        fn seekable(&self) -> bool {
            self.seekable
        }
    }

    #[test]
    fn fetch_exposes_exact_chunks_and_tell_is_post_fetch_offset() {
        let mut input =
            InputStreamBuffer::with_buffer_size(VecStream::from_vec(b"abcdefg".to_vec()), 0, 3);
        assert!(input.fetch().unwrap());
        assert_eq!(input.data(), b"abc");
        assert_eq!(input.tell().unwrap(), 3);
        assert_eq!(input.pop(&mut [0; 2]), 2);
        assert_eq!(input.data(), b"c");
        assert_eq!(input.tell().unwrap(), 3);

        // Like C++, fetch replaces even an incompletely consumed range.
        assert!(input.fetch().unwrap());
        assert_eq!(input.data(), b"def");
        assert_eq!(input.tell().unwrap(), 6);
        assert!(input.fetch().unwrap());
        assert_eq!(input.data(), b"g");
        assert_eq!(input.tell().unwrap(), 7);
        assert!(!input.fetch().unwrap());
    }

    #[test]
    fn read_and_consume_cross_boundaries_without_overrun() {
        let mut input =
            InputStreamBuffer::with_buffer_size(VecStream::from_vec(b"abcdefgh".to_vec()), 0, 3);
        let mut dst = [0; 7];
        assert_eq!(input.read(&mut dst).unwrap(), 7);
        assert_eq!(&dst, b"abcdefg");
        assert_eq!(input.data(), b"h");
        assert_eq!(input.consume(1), Ok(()));
        assert_eq!(input.consume(1), Err(IoError::EndOfStream));
        assert_eq!(input.read(&mut [0; 1]).unwrap(), 0);
    }

    #[test]
    fn async_flag_prefetches_and_swaps_double_buffers() {
        let source = ObservedStream::new(b"abcdefghi", true);
        let reads = Rc::clone(&source.reads);
        let mut input = InputStreamBuffer::with_buffer_size(
            source,
            InputStreamBuffer::<ObservedStream>::ASYNC,
            3,
        );

        assert!(input.fetch().unwrap());
        assert_eq!(input.data(), b"abc");
        assert_eq!(&*reads.borrow(), &[3, 3]);
        // The async worker advances the predecessor, but C++ file_offset_ is
        // updated only by the initial synchronous fetch.
        assert_eq!(input.tell().unwrap(), 3);

        assert!(input.fetch().unwrap());
        assert_eq!(input.data(), b"def");
        assert_eq!(&*reads.borrow(), &[3, 3, 3]);
        assert_eq!(input.tell().unwrap(), 3);
        assert!(input.eof().unwrap());
    }

    #[test]
    fn async_read_error_is_delivered_when_prefetched_buffer_is_fetched() {
        let mut source = ObservedStream::new(b"abcdef", true);
        source.fail_read = Some(1);
        let mut input = InputStreamBuffer::with_buffer_size(
            source,
            InputStreamBuffer::<ObservedStream>::ASYNC,
            3,
        );
        assert!(input.fetch().unwrap());
        assert_eq!(input.data(), b"abc");
        assert_eq!(
            input.fetch(),
            Err(IoError::FileRead("observed".to_string()))
        );
    }

    #[test]
    fn rewind_seek_close_and_nonseekable_tell_match_wrapped_stream() {
        let source = ObservedStream::new(b"abcdef", true);
        let closed = Rc::clone(&source.closed);
        let mut input = InputStreamBuffer::with_buffer_size(source, 0, 2);
        input.fetch().unwrap();
        input.seek(3, SeekFrom::Start(0)).unwrap();
        assert!(input.data().is_empty());
        assert_eq!(input.tell().unwrap(), 0);
        input.fetch().unwrap();
        assert_eq!(input.data(), b"de");
        input.rewind().unwrap();
        assert!(input.data().is_empty());
        assert_eq!(input.tell().unwrap(), 0);
        input.fetch().unwrap();
        assert_eq!(input.data(), b"ab");
        assert_eq!(input.file_size().unwrap(), 6);
        input.close().unwrap();
        assert!(*closed.borrow());

        let mut nonseekable =
            InputStreamBuffer::with_buffer_size(ObservedStream::new(b"x", false), 0, 1);
        assert_eq!(
            nonseekable.tell(),
            Err(IoError::Other(
                "Calling tell on non seekable stream.".to_string()
            ))
        );
    }
}
