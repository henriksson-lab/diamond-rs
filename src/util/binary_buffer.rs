use crate::util::algo::read_varuint32;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BinaryBuffer {
    data: Vec<u8>,
}

impl BinaryBuffer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_vec(data: Vec<u8>) -> Self {
        Self { data }
    }

    pub fn push(&mut self, b: u8) {
        self.data.push(b);
    }

    pub fn extend_from_slice(&mut self, data: &[u8]) {
        self.data.extend_from_slice(data);
    }

    pub fn clear(&mut self) {
        self.data.clear();
    }

    pub fn resize(&mut self, new_len: usize, value: u8) {
        self.data.resize(new_len, value);
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [u8] {
        &mut self.data
    }

    pub fn as_vec(&self) -> &Vec<u8> {
        &self.data
    }

    pub fn as_mut_vec(&mut self) -> &mut Vec<u8> {
        &mut self.data
    }

    pub fn begin(&self) -> BinaryBufferIterator<'_> {
        BinaryBufferIterator::new(&self.data)
    }
}

impl From<Vec<u8>> for BinaryBuffer {
    fn from(data: Vec<u8>) -> Self {
        Self::from_vec(data)
    }
}

impl std::ops::Deref for BinaryBuffer {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl std::ops::DerefMut for BinaryBuffer {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl Extend<u8> for BinaryBuffer {
    fn extend<T: IntoIterator<Item = u8>>(&mut self, iter: T) {
        self.data.extend(iter);
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BinaryBufferIterator<'a> {
    data: &'a [u8],
    ptr: usize,
    end: usize,
}

impl<'a> BinaryBufferIterator<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            ptr: 0,
            end: data.len(),
        }
    }

    pub fn read<T: BinaryReadable>(&mut self) -> Result<T, String> {
        self.check(T::SIZE)?;
        let x = T::read_ne(&self.data[self.ptr..self.ptr + T::SIZE]);
        self.ptr += T::SIZE;
        Ok(x)
    }

    pub fn read_vec<T: BinaryReadable>(&mut self, count: usize) -> Result<Vec<T>, String> {
        let mut v = Vec::with_capacity(count);
        for _ in 0..count {
            v.push(self.read()?);
        }
        Ok(v)
    }

    pub fn read_packed_u32(&mut self, length: u8) -> Result<u32, String> {
        match length {
            0 => Ok(self.read::<u8>()? as u32),
            1 => Ok(self.read::<u16>()? as u32),
            2 => self.read::<u32>(),
            _ => Err("invalid packed integer length".to_string()),
        }
    }

    pub fn read_packed_i32(&mut self, length: u8) -> Result<i32, String> {
        match length {
            0 => Ok(self.read::<u8>()? as i32),
            1 => Ok(self.read::<u16>()? as i32),
            2 => self.read::<i32>(),
            _ => Err("invalid packed integer length".to_string()),
        }
    }

    pub fn read_varint(&mut self) -> Result<u32, String> {
        self.check(1)?;
        let n = (self.end - self.ptr).min(5);
        let (value, consumed) = read_varuint32(&self.data[self.ptr..self.ptr + n])?;
        self.check(consumed)?;
        self.ptr += consumed;
        Ok(value)
    }

    pub fn read_string(&mut self) -> Result<String, String> {
        let start = self.ptr;
        while {
            self.check(1)?;
            let c = self.data[self.ptr];
            self.ptr += 1;
            c != 0
        } {}
        String::from_utf8(self.data[start..self.ptr - 1].to_vec()).map_err(|e| e.to_string())
    }

    pub fn read_bytes(&mut self, size: usize) -> Result<&'a [u8], String> {
        self.check(size)?;
        let begin = self.ptr;
        self.ptr += size;
        Ok(&self.data[begin..self.ptr])
    }

    pub fn remaining(&self) -> &'a [u8] {
        &self.data[self.ptr..self.end]
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn position(&self) -> usize {
        self.ptr
    }

    fn check(&self, size: usize) -> Result<(), String> {
        if self.ptr + size > self.end {
            Err("Unexpected end of file.".to_string())
        } else {
            Ok(())
        }
    }
}

pub trait BinaryReadable: Sized {
    const SIZE: usize;
    fn read_ne(bytes: &[u8]) -> Self;
}

macro_rules! impl_binary_readable {
    ($($t:ty),* $(,)?) => {
        $(
            impl BinaryReadable for $t {
                const SIZE: usize = std::mem::size_of::<$t>();

                fn read_ne(bytes: &[u8]) -> Self {
                    <$t>::from_ne_bytes(bytes.try_into().unwrap())
                }
            }
        )*
    };
}

impl_binary_readable!(u8, i8, u16, i16, u32, i32, u64, i64, f32, f64);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::algo::write_varuint32;

    #[test]
    fn test_binary_buffer_read_scalars_vec_and_string() {
        let mut data = Vec::new();
        data.extend_from_slice(&7u8.to_ne_bytes());
        data.extend_from_slice(&0x1234u16.to_ne_bytes());
        data.extend_from_slice(&[1u32.to_ne_bytes(), 2u32.to_ne_bytes()].concat());
        data.extend_from_slice(b"abc\0");
        let b = BinaryBuffer::from(data);
        let mut it = b.begin();
        assert_eq!(it.read::<u8>().unwrap(), 7);
        assert_eq!(it.read::<u16>().unwrap(), 0x1234);
        assert_eq!(it.read_vec::<u32>(2).unwrap(), vec![1, 2]);
        assert_eq!(it.read_string().unwrap(), "abc");
        assert_eq!(it.remaining(), &[] as &[u8]);
        assert!(!it.good());
    }

    #[test]
    fn test_binary_buffer_read_bytes_and_remaining() {
        let b = BinaryBuffer::from(vec![1, 2, 3, 4]);
        let mut it = b.begin();
        assert_eq!(it.read_bytes(2).unwrap(), &[1, 2]);
        assert_eq!(it.remaining(), &[3, 4]);
        assert_eq!(it.read_bytes(2).unwrap(), &[3, 4]);
        assert!(it.read_bytes(1).is_err());
    }

    #[test]
    fn test_binary_buffer_vector_surface() {
        let mut b = BinaryBuffer::new();
        assert!(b.is_empty());
        b.push(1);
        b.extend_from_slice(&[2, 3]);
        b.extend([4, 5]);
        assert_eq!(b.len(), 5);
        assert_eq!(&b[..], &[1, 2, 3, 4, 5]);
        b.data_mut()[1] = 9;
        assert_eq!(b.as_vec(), &vec![1, 9, 3, 4, 5]);
        b.as_mut_vec().push(6);
        assert_eq!(b.data(), &[1, 9, 3, 4, 5, 6]);
        b.resize(3, 0);
        assert_eq!(b.data(), &[1, 9, 3]);
        b.clear();
        assert!(b.is_empty());
    }

    #[test]
    fn test_binary_buffer_read_packed_and_varint() {
        let mut data = Vec::new();
        data.extend_from_slice(&255u8.to_ne_bytes());
        data.extend_from_slice(&256u16.to_ne_bytes());
        data.extend_from_slice(&70_000u32.to_ne_bytes());
        write_varuint32(16_384, &mut data);
        let mut it = BinaryBufferIterator::new(&data);
        assert_eq!(it.read_packed_u32(0).unwrap(), 255);
        assert_eq!(it.read_packed_i32(1).unwrap(), 256);
        assert_eq!(it.read_packed_u32(2).unwrap(), 70_000);
        assert_eq!(it.read_varint().unwrap(), 16_384);
        assert!(it.read::<u8>().is_err());
    }
}
