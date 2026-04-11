use super::value::Letter;

/// Packed seed representation (uint64_t in C++).
pub type PackedSeed = u64;

/// Seed offset type for indexing.
pub type SeedOffset = u32;

/// Seed partition type.
pub type SeedPartition = u32;

/// Maximum seed weight (number of positions in a spaced seed pattern).
pub const MAX_SEED_WEIGHT: usize = 32;

#[inline]
pub fn seedp_mask(seedp_bits: i32) -> PackedSeed {
    (1u64 << seedp_bits) - 1
}

#[inline]
pub fn seedp_count(seedp_bits: i32) -> PackedSeed {
    1u64 << seedp_bits
}

#[inline]
pub fn seed_partition(s: PackedSeed, mask: PackedSeed) -> SeedPartition {
    (s & mask) as SeedPartition
}

#[inline]
pub fn seed_partition_offset(s: PackedSeed, seedp_bits: PackedSeed) -> SeedOffset {
    (s >> seedp_bits) as SeedOffset
}

/// A seed as an array of letter values (expanded form).
#[derive(Clone, Default)]
pub struct Seed {
    data: [Letter; MAX_SEED_WEIGHT],
}

impl Seed {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn get(&self, i: usize) -> Letter {
        self.data[i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, val: Letter) {
        self.data[i] = val;
    }

    /// Convert to packed representation using alphabet size 20.
    pub fn packed(&self, weight: usize) -> u64 {
        let mut s = 0u64;
        for i in 0..weight {
            s *= 20;
            s += self.data[i] as u64;
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seedp_mask() {
        assert_eq!(seedp_mask(4), 0xF);
        assert_eq!(seedp_mask(8), 0xFF);
    }

    #[test]
    fn test_seed_partition() {
        let s: PackedSeed = 0xABCD;
        assert_eq!(seed_partition(s, 0xFF), 0xCD);
    }

    #[test]
    fn test_seed_packed() {
        let mut seed = Seed::new();
        seed.set(0, 1);
        seed.set(1, 2);
        // packed = 1*20 + 2 = 22
        assert_eq!(seed.packed(2), 22);
    }
}
