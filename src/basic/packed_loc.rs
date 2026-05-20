#[repr(C, packed)]
#[derive(Clone, Copy, Default, Eq)]
pub struct PackedUint40 {
    pub high: u8,
    pub low: u32,
}

pub type PackedLoc = PackedUint40;

impl PackedUint40 {
    pub fn new(v: u64) -> Self {
        Self {
            high: (v >> 32) as u8,
            low: (v & 0xffff_ffff) as u32,
        }
    }

    pub fn assign_u32(&mut self, x: u32) {
        self.high = 0;
        self.low = x;
    }

    pub fn as_u64(self) -> u64 {
        ((self.high as u64) << 32) | self.low as u64
    }

    pub fn as_i64(self) -> i64 {
        ((self.high as i64) << 32) | self.low as i64
    }

    pub fn as_u32(self) -> u32 {
        self.low
    }

    pub fn as_i32(self) -> i32 {
        self.low as i32
    }
}

impl std::fmt::Debug for PackedUint40 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let high = self.high;
        let low = self.low;
        f.debug_struct("PackedUint40")
            .field("high", &high)
            .field("low", &low)
            .finish()
    }
}

impl PartialEq for PackedUint40 {
    fn eq(&self, rhs: &Self) -> bool {
        self.high == rhs.high && self.low == rhs.low
    }

    // NOTE: this `ne` is intentionally *not* `!self.eq(rhs)`. It mirrors the
    // C++ bug in `diamond/src/basic/packed_loc.h:61` which writes
    // `high != rhs.high && low != rhs.low` — so for `(1, 100)` vs `(2, 100)`
    // BOTH `==` and `!=` return false. We keep the buggy semantics so any
    // C++ code path relying on it gives the same answer in Rust.
    #[allow(clippy::partialeq_ne_impl)]
    fn ne(&self, rhs: &Self) -> bool {
        self.high != rhs.high && self.low != rhs.low
    }
}

impl From<u64> for PackedUint40 {
    fn from(value: u64) -> Self {
        Self::new(value)
    }
}

impl From<u32> for PackedUint40 {
    fn from(value: u32) -> Self {
        Self {
            high: 0,
            low: value,
        }
    }
}

impl From<PackedUint40> for u64 {
    fn from(value: PackedUint40) -> Self {
        value.as_u64()
    }
}

impl PartialOrd for PackedUint40 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PackedUint40 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let high = self.high;
        let other_high = other.high;
        let low = self.low;
        let other_low = other.low;
        high.cmp(&other_high).then_with(|| low.cmp(&other_low))
    }
}

impl std::ops::Sub for PackedUint40 {
    type Output = u64;

    fn sub(self, rhs: PackedUint40) -> Self::Output {
        self.as_u64() - rhs.as_u64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_packed_uint40_roundtrip_and_casts() {
        let p = PackedUint40::new(0x12_3456_789a);
        assert_eq!(p.high, 0x12);
        let low = p.low;
        assert_eq!(low, 0x3456_789a);
        assert_eq!(p.as_u64(), 0x12_3456_789a);
        assert_eq!(p.as_i64(), 0x12_3456_789a);
        assert_eq!(p.as_u32(), 0x3456_789a);
        assert_eq!(p.as_i32(), 0x3456_789a_u32 as i32);
    }

    #[test]
    fn test_packed_uint40_assign_compare_and_sub() {
        let mut p = PackedUint40::new(0x01_0000_0005);
        p.assign_u32(7);
        assert_eq!(p, PackedUint40::from(7u32));
        assert!(PackedUint40::new(0x01_0000_0000) > PackedUint40::new(0xffff_ffff));
        assert_eq!(
            PackedUint40::new(0x01_0000_0005) - PackedUint40::new(0x01_0000_0000),
            5
        );
    }

    #[test]
    fn test_packed_uint40_layout_and_cpp_ne() {
        assert_eq!(std::mem::size_of::<PackedUint40>(), 5);
        assert_eq!(std::mem::align_of::<PackedUint40>(), 1);

        let a = PackedUint40::new(0x01_0000_0001);
        let same_high = PackedUint40::new(0x01_0000_0002);
        let same_low = PackedUint40::new(0x02_0000_0001);
        assert!(!(a != same_high));
        assert!(!(a != same_low));
        assert!(a != PackedUint40::new(0x02_0000_0002));
    }
}
