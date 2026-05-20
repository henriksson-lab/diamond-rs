use crate::util::algo::write_varuint32;
use crate::util::string::format_double;

pub fn find_first_of(s: &str, delimiters: &str) -> usize {
    let bytes = s.as_bytes();
    let del = delimiters.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() && bytes[i] != 0 && !del.contains(&bytes[i]) {
        i += 1;
    }
    i
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TextBuffer {
    data: Vec<u8>,
    alloc_size: usize,
}

impl TextBuffer {
    const BLOCK_SIZE: usize = 4096;

    pub fn new() -> Self {
        Self::default()
    }

    pub fn reserve(&mut self, n: usize) {
        let s = self.data.len();
        if s + n < self.alloc_size {
            return;
        }
        self.alloc_size = s + n + Self::BLOCK_SIZE - ((s + n) & (Self::BLOCK_SIZE - 1));
        self.data.reserve(self.alloc_size.saturating_sub(s));
    }

    pub fn advance(&mut self, n: usize) {
        self.reserve(n);
        self.data.resize(self.data.len() + n, 0);
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [u8] {
        &mut self.data
    }

    pub fn clear(&mut self) {
        self.data.clear();
    }

    pub fn write<T: TextBufferPrimitive>(&mut self, data: T) -> &mut Self {
        let bytes = data.to_ne_bytes_vec();
        self.reserve(bytes.len());
        self.data.extend_from_slice(&bytes);
        self
    }

    pub fn write_c_str(&mut self, s: &str) {
        self.reserve(s.len() + 1);
        self.data.extend_from_slice(s.as_bytes());
        self.data.push(0);
    }

    pub fn write_c_str_len(&mut self, s: &str, len: usize) {
        self.reserve(len + 1);
        self.data.extend_from_slice(&s.as_bytes()[..len]);
        self.data.push(0);
    }

    pub fn write_raw(&mut self, s: &[u8]) {
        self.reserve(s.len());
        self.data.extend_from_slice(s);
    }

    pub fn push_back(&mut self, c: u8) {
        self.data.push(c);
    }

    pub fn write_until(&mut self, s: &str, delimiters: &str) {
        let n = find_first_of(s, delimiters);
        self.write_raw(&s.as_bytes()[..n]);
    }

    pub fn write_packed(&mut self, x: u32) -> &mut Self {
        if x <= u8::MAX as u32 {
            self.write(x as u8)
        } else if x <= u16::MAX as u32 {
            self.write(x as u16)
        } else {
            self.write(x)
        }
    }

    pub fn write_varint(&mut self, x: u32) -> &mut Self {
        let mut buf = Vec::new();
        write_varuint32(x, &mut buf);
        self.write_raw(&buf);
        self
    }

    pub fn size(&self) -> usize {
        self.data.len()
    }

    pub fn alloc_size(&self) -> usize {
        self.alloc_size
    }

    pub fn append_str(&mut self, s: &str) -> &mut Self {
        self.reserve(s.len());
        self.data.extend_from_slice(s.as_bytes());
        self
    }

    pub fn append_char(&mut self, c: char) -> &mut Self {
        let mut buf = [0u8; 4];
        let s = c.encode_utf8(&mut buf);
        self.append_str(s)
    }

    pub fn append_display<T: std::fmt::Display>(&mut self, x: T) -> &mut Self {
        self.append_str(&x.to_string())
    }

    pub fn append_double(&mut self, x: f64) -> &mut Self {
        self.append_str(&format_double(x))
    }

    pub fn print_d(&mut self, x: f64) -> &mut Self {
        self.append_str(&format!("{:.6}", x))
    }

    pub fn print_e(&mut self, x: f64) -> &mut Self {
        // C++ `TextBuffer::print_e` (`util/text_buffer.h:249-257`) emits
        // `snprintf("%.2e", x)` which has an explicit `+`/`-` exponent sign
        // and pads single-digit exponents to two digits (`1.23e-05`, not
        // `1.23e-5`). Rust's `{:.2e}` omits the sign on positive exponents
        // and the leading zero — fix it up here to match C printf.
        if x == 0.0 {
            return self.append_str("0.0");
        }
        let s = format!("{:.2e}", x);
        let bytes = s.as_bytes();
        let Some(epos) = bytes.iter().position(|&b| b == b'e') else {
            return self.append_str(&s);
        };
        let (mantissa, exp_part) = s.split_at(epos);
        let exp = &exp_part[1..];
        let (sign, digits) = match exp.as_bytes().first() {
            Some(b'-') => ("-", &exp[1..]),
            Some(b'+') => ("+", &exp[1..]),
            _ => ("+", exp),
        };
        let formatted = if digits.len() == 1 {
            format!("{}e{}0{}", mantissa, sign, digits)
        } else {
            format!("{}e{}{}", mantissa, sign, digits)
        };
        self.append_str(&formatted)
    }

    pub fn print(&mut self, i: u32, width: usize) -> &mut Self {
        let s = format!("{:>width$}", i, width = width);
        self.reserve(16);
        self.data
            .extend_from_slice(&s.as_bytes()[..width.min(s.len())]);
        self
    }

    pub fn print_vec<T: std::fmt::Display>(&mut self, v: &[T], _separator: char) -> &mut Self {
        if v.is_empty() {
            return self;
        }
        self.append_display(&v[0]);
        for item in &v[1..] {
            self.append_char(';');
            self.append_display(item);
        }
        self
    }

    pub fn write_vec<T: TextBufferPrimitive>(&mut self, v: &[T]) -> &mut Self {
        for item in v {
            self.write(*item);
        }
        self
    }

    pub fn get_mut(&mut self, pos: usize) -> &mut u8 {
        &mut self.data[pos]
    }
}

pub trait TextBufferPrimitive: Copy {
    fn to_ne_bytes_vec(self) -> Vec<u8>;
}

macro_rules! impl_text_buffer_primitive {
    ($($t:ty),* $(,)?) => {
        $(
            impl TextBufferPrimitive for $t {
                fn to_ne_bytes_vec(self) -> Vec<u8> {
                    self.to_ne_bytes().to_vec()
                }
            }
        )*
    };
}

impl_text_buffer_primitive!(u8, i8, u16, i16, u32, i32, u64, i64, usize, isize, f32, f64);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::algo::read_varuint32;

    #[test]
    fn test_find_first_of() {
        assert_eq!(find_first_of("abc,def", ","), 3);
        assert_eq!(find_first_of("abcdef", ","), 6);
        assert_eq!(find_first_of("abc\0def,ghi", ","), 3);
    }

    #[test]
    fn test_text_buffer_write_methods() {
        let mut b = TextBuffer::new();
        b.reserve(1);
        assert!(b.alloc_size() >= 4096);
        b.write(7u8).write(0x1234u16);
        b.write_c_str("abc");
        b.write_c_str_len("defghi", 3);
        assert_eq!(b.data()[0], 7);
        assert_eq!(&b.data()[3..7], b"abc\0");
        assert_eq!(&b.data()[7..11], b"def\0");
        b.clear();
        assert_eq!(b.size(), 0);
    }

    #[test]
    fn test_text_buffer_print_and_packed() {
        let mut b = TextBuffer::new();
        b.append_str("x")
            .append_char(':')
            .append_display(12u32)
            .append_double(9.94)
            .print_d(1.5)
            .print_e(0.0)
            .print(7, 3)
            .print_vec(&[1, 2, 3], ',');
        assert_eq!(
            String::from_utf8(b.data().to_vec()).unwrap(),
            "x:129.91.5000000.0  71;2;3"
        );

        b.clear();
        b.write_packed(255).write_packed(256);
        assert_eq!(&b.data()[..1], &[255]);
        assert_eq!(u16::from_ne_bytes([b.data()[1], b.data()[2]]), 256);

        b.clear();
        b.write_varint(16_384);
        assert_eq!(read_varuint32(b.data()).unwrap().0, 16_384);

        b.clear();
        b.print(1234, 2);
        assert_eq!(b.data(), b"12");
    }

    #[test]
    fn test_text_buffer_write_until_and_vec() {
        let mut b = TextBuffer::new();
        b.write_until("abc,def", ",");
        assert_eq!(b.data(), b"abc");
        b.push_back(b'!');
        assert_eq!(b.data(), b"abc!");
        *b.get_mut(3) = b'?';
        assert_eq!(b.data(), b"abc?");
        b.clear();
        b.write_vec(&[1u16, 2u16]);
        assert_eq!(b.size(), 4);
        b.advance(2);
        assert_eq!(b.size(), 6);
    }
}
