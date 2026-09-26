// MurmurHash3 implementation (public domain, by Austin Appleby).
// This is the x64_128 variant used by DIAMOND for database hashing
// and test output verification.

pub fn hash64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^= x >> 31;
    x
}

pub fn murmur_hash_u64(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    h ^= h >> 33;
    h
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MurmurHash;

impl MurmurHash {
    pub fn call(&self, h: u64) -> u64 {
        murmur_hash_u64(h)
    }
}

/// MurmurHash3_x64_128 with a 16-byte seed.
///
/// This matches the C++ signature: `MurmurHash3_x64_128(key, len, seed, out)`
/// where seed is a `char[16]` used as two u64 values.
pub fn murmurhash3_x64_128(data: &[u8], seed: &[u8; 16]) -> [u8; 16] {
    crate::murmurhash::murmur_hash3::murmur_hash3_x64_128_bytes(data, seed)
}

/// Compute the iterative hash used by DIAMOND for output file verification.
///
/// This matches the C++ `InputFile::hash()` function which processes
/// data in 4096-byte chunks, chaining the hash as a seed.
pub fn file_hash(data: &[u8]) -> u64 {
    let mut seed = [0u8; 16];
    let chunk_size = 4096;

    let mut offset = 0;
    while offset < data.len() {
        let end = (offset + chunk_size).min(data.len());
        let chunk = &data[offset..end];
        seed = murmurhash3_x64_128(chunk, &seed);
        offset = end;
    }

    u64::from_ne_bytes(seed[0..8].try_into().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_murmurhash3_nonempty() {
        let seed = [0u8; 16];
        let result = murmurhash3_x64_128(b"test", &seed);
        assert_ne!(result, [0u8; 16]);
    }

    #[test]
    fn test_hash64() {
        assert_eq!(hash64(0), 0xe220_a839_7b1d_cdaf);
        assert_eq!(hash64(1), 0x910a_2dec_8902_5cc1);
        assert_eq!(hash64(u64::MAX), 0xe4d9_7177_1b65_2c20);
    }

    #[test]
    fn test_murmur_hash_u64() {
        assert_eq!(murmur_hash_u64(0), 0);
        assert_eq!(murmur_hash_u64(1), 0xb456_bcfc_34c2_cb2c);
        assert_eq!(murmur_hash_u64(u64::MAX), 0x64b5_720b_4b82_5f21);
        let h = MurmurHash;
        assert_eq!(h.call(1), murmur_hash_u64(1));
    }

    #[test]
    fn test_murmurhash3_deterministic() {
        let seed = [0u8; 16];
        let data = b"Hello, World!";
        let h1 = murmurhash3_x64_128(data, &seed);
        let h2 = murmurhash3_x64_128(data, &seed);
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_murmurhash3_different_inputs() {
        let seed = [0u8; 16];
        let h1 = murmurhash3_x64_128(b"abc", &seed);
        let h2 = murmurhash3_x64_128(b"abd", &seed);
        assert_ne!(h1, h2);
    }

    #[test]
    fn test_file_hash() {
        let data = vec![0u8; 8192]; // Two full chunks
        let h = file_hash(&data);
        assert_ne!(h, 0);
    }

    #[test]
    fn test_file_hash_small() {
        let data = b"test data";
        let h = file_hash(data);
        assert_ne!(h, 0);
    }
}
