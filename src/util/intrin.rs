pub fn popcount32(x: u32) -> u32 {
    x.count_ones()
}

pub fn popcount64(x: u64) -> u32 {
    x.count_ones()
}

pub fn ctz32(x: u32) -> i32 {
    x.trailing_zeros() as i32
}

pub fn ctz64(x: u64) -> i32 {
    x.trailing_zeros() as i32
}

pub fn clz32(x: u32) -> i32 {
    x.leading_zeros() as i32
}

pub fn clz64(x: u64) -> i32 {
    x.leading_zeros() as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_popcount() {
        assert_eq!(popcount32(0), 0);
        assert_eq!(popcount32(0b1011), 3);
        assert_eq!(popcount64(0xffff_0000_ffff_0000), 32);
    }

    #[test]
    fn test_ctz_and_clz() {
        assert_eq!(ctz32(0b1000), 3);
        assert_eq!(ctz64(0x1_0000_0000), 32);
        assert_eq!(clz32(1), 31);
        assert_eq!(clz64(1), 63);
    }
}
