pub fn is_little_endian() -> bool {
    cfg!(target_endian = "little")
}

pub trait BigEndianByteswap {
    fn big_endian_byteswap(self) -> Self;
}

macro_rules! impl_big_endian_byteswap {
    ($($t:ty),* $(,)?) => {
        $(
            impl BigEndianByteswap for $t {
                fn big_endian_byteswap(self) -> Self {
                    <$t>::from_be(self)
                }
            }
        )*
    };
}

impl_big_endian_byteswap!(u16, u32, u64, i16, i32, i64);

pub fn big_endian_byteswap<T: BigEndianByteswap>(x: T) -> T {
    x.big_endian_byteswap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_endianness_and_byteswap() {
        assert_eq!(is_little_endian(), cfg!(target_endian = "little"));
        if cfg!(target_endian = "little") {
            assert_eq!(big_endian_byteswap(0x1234u16), 0x3412);
            assert_eq!(big_endian_byteswap(0x0102_0304u32), 0x0403_0201);
            assert_eq!(
                big_endian_byteswap(0x0102_0304_0506_0708u64),
                0x0807_0605_0403_0201
            );
        } else {
            assert_eq!(big_endian_byteswap(0x1234u16), 0x1234);
        }
    }
}
