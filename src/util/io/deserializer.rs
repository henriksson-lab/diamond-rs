//! Translation of `diamond/src/util/io/deserializer.{h,cpp}`.

use std::fs::File;
use std::io::SeekFrom;

use super::{FilePrimitive, IoError, IoResult, StreamEntity};

#[derive(Debug, Clone)]
pub struct Deserializer<S: StreamEntity> {
    buffer: S,
    // C++ leaves `record_start` at `buffer_->begin`. A one-byte lookahead
    // provides the same behavior for generic and non-seekable Rust streams.
    lookahead: Option<u8>,
    lookahead_data: [u8; 1],
}

impl<S: StreamEntity> Deserializer<S> {
    pub fn new(buffer: S) -> Self {
        Self {
            buffer,
            lookahead: None,
            lookahead_data: [0],
        }
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.lookahead = None;
        self.buffer.rewind()
    }

    pub fn seek(&mut self, pos: i64) -> IoResult<&mut Self> {
        self.lookahead = None;
        self.buffer.seek(pos, SeekFrom::Start(0))?;
        Ok(self)
    }

    pub fn seek_forward(&mut self, n: usize) -> IoResult<()> {
        if n == 0 {
            return Ok(());
        }
        let remaining = n - usize::from(self.lookahead.take().is_some());
        if remaining == 0 {
            Ok(())
        } else {
            self.buffer.seek(remaining as i64, SeekFrom::Current(0))
        }
    }

    pub fn seek_forward_delim(&mut self, delimiter: u8) -> IoResult<bool> {
        while let Some(byte) = self.read_one()? {
            if byte == delimiter {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.lookahead = None;
        self.buffer.close()
    }

    pub fn read_u32(&mut self) -> IoResult<u32> {
        let mut bytes = [0; 4];
        self.read_exact(&mut bytes)?;
        Ok(u32::from_le_bytes(bytes))
    }

    pub fn read_i32(&mut self) -> IoResult<i32> {
        let mut bytes = [0; 4];
        self.read_exact(&mut bytes)?;
        Ok(i32::from_le_bytes(bytes))
    }

    pub fn read_i16(&mut self) -> IoResult<i16> {
        let mut bytes = [0; 2];
        self.read_exact(&mut bytes)?;
        Ok(i16::from_le_bytes(bytes))
    }

    pub fn read_u16(&mut self) -> IoResult<u16> {
        let mut bytes = [0; 2];
        self.read_exact(&mut bytes)?;
        Ok(u16::from_le_bytes(bytes))
    }

    pub fn read_i64(&mut self) -> IoResult<i64> {
        let mut bytes = [0; 8];
        self.read_exact(&mut bytes)?;
        Ok(i64::from_le_bytes(bytes))
    }

    pub fn read_u64(&mut self) -> IoResult<u64> {
        let mut bytes = [0; 8];
        self.read_exact(&mut bytes)?;
        Ok(u64::from_le_bytes(bytes))
    }

    /// C++ deliberately does not apply `big_endian_byteswap` to doubles.
    pub fn read_f64(&mut self) -> IoResult<f64> {
        let mut bytes = [0; 8];
        self.read_exact(&mut bytes)?;
        Ok(f64::from_ne_bytes(bytes))
    }

    pub fn read_value<T: FilePrimitive>(&mut self) -> IoResult<T> {
        let mut bytes = vec![0; T::SIZE];
        self.read_exact(&mut bytes)?;
        T::from_ne_bytes_slice(&bytes)
    }

    pub fn read_values<T: FilePrimitive>(&mut self, ptr: &mut [T]) -> IoResult<()> {
        for value in ptr {
            *value = self.read_value()?;
        }
        Ok(())
    }

    /// C++ `read(T* ptr, size_t count)`, returning the number of complete
    /// elements read rather than requiring the whole slice.
    pub fn read_items<T: FilePrimitive>(&mut self, ptr: &mut [T]) -> IoResult<usize> {
        // Initialize from the destination because C++ writes directly into
        // it: a trailing partial element changes only the bytes actually read.
        let mut bytes = Vec::with_capacity(ptr.len() * T::SIZE);
        for &value in ptr.iter() {
            bytes.extend_from_slice(&value.to_ne_bytes_vec());
        }
        let byte_count = self.read_raw(&mut bytes)?;
        let item_count = byte_count / T::SIZE;
        for (value, chunk) in ptr.iter_mut().zip(bytes.chunks_exact(T::SIZE)) {
            *value = T::from_ne_bytes_slice(chunk)?;
        }
        Ok(item_count)
    }

    pub fn read_string(&mut self) -> IoResult<String> {
        let mut out = Vec::new();
        if !self.read_to(&mut out, 0)? {
            return Err(IoError::EndOfStream);
        }
        String::from_utf8(out).map_err(|e| IoError::Other(e.to_string()))
    }

    pub fn read_exact(&mut self, ptr: &mut [u8]) -> IoResult<()> {
        if self.read_raw(ptr)? == ptr.len() {
            Ok(())
        } else {
            Err(IoError::EndOfStream)
        }
    }

    pub fn read_raw(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let mut total = 0;
        if !ptr.is_empty() {
            if let Some(byte) = self.lookahead.take() {
                ptr[0] = byte;
                total = 1;
            }
        }
        while total < ptr.len() {
            let n = self.buffer.read(&mut ptr[total..])?;
            if n == 0 {
                break;
            }
            total += n;
        }
        Ok(total)
    }

    pub fn pop(&mut self, dst: &mut [u8]) -> IoResult<()> {
        self.read_exact(dst)
    }

    pub fn read_to(&mut self, dst: &mut Vec<u8>, delimiter: u8) -> IoResult<bool> {
        while let Some(byte) = self.read_one()? {
            if byte == delimiter {
                return Ok(true);
            }
            dst.push(byte);
        }
        Ok(false)
    }

    pub fn data(&self) -> &[u8] {
        if self.lookahead.is_some() {
            &self.lookahead_data
        } else {
            self.buffer.data()
        }
    }

    /// C++ protected `avail()`; public for compatibility with Rust wrappers
    /// that cannot use inheritance.
    pub fn avail(&self) -> usize {
        usize::from(self.lookahead.is_some()) + self.buffer.data().len()
    }

    pub fn read_to_record_start(
        &mut self,
        dst: &mut Vec<u8>,
        line_delimiter: u8,
        record_start: u8,
    ) -> IoResult<(bool, i64)> {
        let mut newlines = 0;
        while let Some(byte) = self.read_one()? {
            if byte != line_delimiter {
                dst.push(byte);
                continue;
            }
            newlines += 1;
            match self.read_one()? {
                None => return Ok((false, newlines)),
                Some(next) => {
                    self.lookahead = Some(next);
                    self.lookahead_data[0] = next;
                    if next == record_start {
                        return Ok((true, newlines));
                    }
                    dst.push(line_delimiter);
                }
            }
        }
        Ok((false, newlines))
    }

    pub fn peek(&mut self, n: usize) -> IoResult<String> {
        if n == 0 {
            return Ok(String::new());
        }
        if self.lookahead.is_none() && self.buffer.data().is_empty() {
            match self.buffer.fetch() {
                Ok(_) => {}
                Err(IoError::UnsupportedOperation) => return self.peek_via_seek(n),
                Err(error) => return Err(error),
            }
        }
        let available = self.available_prefix(n);
        if available.len() < n && !self.buffer.eof()? {
            return Err(IoError::Other("Invalid peek".to_string()));
        }
        String::from_utf8(available).map_err(|e| IoError::Other(e.to_string()))
    }

    pub fn file_size(&mut self) -> IoResult<i64> {
        self.buffer.file_size()
    }

    pub fn file(&mut self) -> IoResult<&mut File> {
        self.buffer.file()
    }

    pub fn into_inner(self) -> S {
        self.buffer
    }

    fn read_one(&mut self) -> IoResult<Option<u8>> {
        if let Some(byte) = self.lookahead.take() {
            return Ok(Some(byte));
        }
        let mut byte = [0];
        Ok(match self.buffer.read(&mut byte)? {
            0 => None,
            _ => Some(byte[0]),
        })
    }

    fn available_prefix(&self, n: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(n);
        if let Some(byte) = self.lookahead {
            out.push(byte);
        }
        let take = n.saturating_sub(out.len()).min(self.buffer.data().len());
        out.extend_from_slice(&self.buffer.data()[..take]);
        out
    }

    /// Compatibility for Rust callers that wrap a seekable `StreamEntity`
    /// directly rather than the C++-equivalent `InputStreamBuffer`.
    fn peek_via_seek(&mut self, n: usize) -> IoResult<String> {
        let pos = self.buffer.tell()?;
        let mut out = Vec::with_capacity(n);
        if let Some(byte) = self.lookahead {
            out.push(byte);
        }
        let mut tail = vec![0; n - out.len()];
        let read = self.buffer.read(&mut tail)?;
        let eof = self.buffer.eof()?;
        self.buffer.seek(pos, SeekFrom::Start(0))?;
        tail.truncate(read);
        out.extend_from_slice(&tail);
        if out.len() < n && !eof {
            return Err(IoError::Other("Invalid peek".to_string()));
        }
        String::from_utf8(out).map_err(|e| IoError::Other(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::{InputStreamBuffer, VecStream};

    fn buffered(bytes: &[u8]) -> Deserializer<InputStreamBuffer<VecStream>> {
        Deserializer::new(InputStreamBuffer::new(
            VecStream::from_vec(bytes.to_vec()),
            0,
        ))
    }

    #[test]
    fn reads_little_endian_scalars_and_reports_short_binary_input() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&0x0102_0304u32.to_le_bytes());
        bytes.extend_from_slice(&(-7i16).to_le_bytes());
        bytes.extend_from_slice(&1.5f64.to_ne_bytes());
        let mut input = buffered(&bytes);
        assert_eq!(input.read_u32().unwrap(), 0x0102_0304);
        assert_eq!(input.read_i16().unwrap(), -7);
        assert_eq!(input.read_f64().unwrap(), 1.5);
        assert_eq!(input.read_u32().unwrap_err(), IoError::EndOfStream);

        let mut input = buffered(&[1, 0, 2]);
        let mut values = [99u16; 2];
        assert_eq!(input.read_items(&mut values).unwrap(), 1);
        let mut second = 99u16.to_ne_bytes();
        second[0] = 2;
        assert_eq!(values, [1, u16::from_ne_bytes(second)]);
    }

    #[test]
    fn seek_rewind_delimiter_and_partial_raw_read() {
        let mut input = buffered(b"zero,one");
        assert!(input.seek_forward_delim(b',').unwrap());
        let mut tail = [0; 8];
        assert_eq!(input.read_raw(&mut tail).unwrap(), 3);
        assert_eq!(&tail[..3], b"one");
        input.rewind().unwrap();
        input.seek(5).unwrap();
        let mut value = [0; 3];
        input.pop(&mut value).unwrap();
        assert_eq!(&value, b"one");
    }

    #[test]
    fn record_start_is_not_consumed_on_non_seekable_stream() {
        #[derive(Debug, Clone)]
        struct NonSeekable(VecStream);
        impl StreamEntity for NonSeekable {
            fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
                self.0.read(ptr)
            }
            fn eof(&mut self) -> IoResult<bool> {
                self.0.eof()
            }
        }

        let mut input = Deserializer::new(NonSeekable(VecStream::from_vec(b"abc\n>next".to_vec())));
        let mut record = Vec::new();
        assert_eq!(
            input
                .read_to_record_start(&mut record, b'\n', b'>')
                .unwrap(),
            (true, 1)
        );
        assert_eq!(record, b"abc");
        input.seek_forward(1).unwrap();
        let mut marker = [0];
        input.read_exact(&mut marker).unwrap();
        assert_eq!(marker[0], b'n');
    }

    #[test]
    fn peek_matches_eof_and_invalid_peek_rules() {
        let mut input = buffered(b"abc");
        assert_eq!(input.peek(3).unwrap(), "abc");
        assert_eq!(input.peek(5).unwrap(), "abc");

        #[derive(Debug, Clone)]
        struct NotEofBuffer {
            data: Vec<u8>,
            fetched: bool,
        }
        impl StreamEntity for NotEofBuffer {
            fn fetch(&mut self) -> IoResult<bool> {
                self.fetched = true;
                Ok(true)
            }
            fn data(&self) -> &[u8] {
                if self.fetched {
                    &self.data
                } else {
                    &[]
                }
            }
            fn eof(&mut self) -> IoResult<bool> {
                Ok(false)
            }
        }
        let mut input = Deserializer::new(NotEofBuffer {
            data: b"ab".to_vec(),
            fetched: false,
        });
        assert_eq!(
            input.peek(3).unwrap_err(),
            IoError::Other("Invalid peek".to_string())
        );
    }
}
