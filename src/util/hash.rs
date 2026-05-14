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

#[inline]
fn rotl64(x: u64, r: u32) -> u64 {
    x.rotate_left(r)
}

#[inline]
fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;
    k
}

/// MurmurHash3_x64_128 with a 16-byte seed.
///
/// This matches the C++ signature: `MurmurHash3_x64_128(key, len, seed, out)`
/// where seed is a `char[16]` used as two u64 values.
pub fn murmurhash3_x64_128(data: &[u8], seed: &[u8; 16]) -> [u8; 16] {
    let len = data.len();
    let nblocks = len / 16;

    let mut h1 = u64::from_le_bytes(seed[0..8].try_into().unwrap());
    let mut h2 = u64::from_le_bytes(seed[8..16].try_into().unwrap());

    let c1: u64 = 0x87c37b91114253d5;
    let c2: u64 = 0x4cf5ad432745937f;

    // Body
    for i in 0..nblocks {
        let offset = i * 16;
        let mut k1 = u64::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
        let mut k2 = u64::from_le_bytes(data[offset + 8..offset + 16].try_into().unwrap());

        k1 = k1.wrapping_mul(c1);
        k1 = rotl64(k1, 31);
        k1 = k1.wrapping_mul(c2);
        h1 ^= k1;

        h1 = rotl64(h1, 27);
        h1 = h1.wrapping_add(h2);
        h1 = h1.wrapping_mul(5).wrapping_add(0x52dce729);

        k2 = k2.wrapping_mul(c2);
        k2 = rotl64(k2, 33);
        k2 = k2.wrapping_mul(c1);
        h2 ^= k2;

        h2 = rotl64(h2, 31);
        h2 = h2.wrapping_add(h1);
        h2 = h2.wrapping_mul(5).wrapping_add(0x38495ab5);
    }

    // Tail
    let tail = &data[nblocks * 16..];
    let mut k1: u64 = 0;
    let mut k2: u64 = 0;

    // Fall-through switch emulation
    let tail_len = len & 15;
    if tail_len >= 15 {
        k2 ^= (tail[14] as u64) << 48;
    }
    if tail_len >= 14 {
        k2 ^= (tail[13] as u64) << 40;
    }
    if tail_len >= 13 {
        k2 ^= (tail[12] as u64) << 32;
    }
    if tail_len >= 12 {
        k2 ^= (tail[11] as u64) << 24;
    }
    if tail_len >= 11 {
        k2 ^= (tail[10] as u64) << 16;
    }
    if tail_len >= 10 {
        k2 ^= (tail[9] as u64) << 8;
    }
    if tail_len >= 9 {
        k2 ^= tail[8] as u64;
        k2 = k2.wrapping_mul(c2);
        k2 = rotl64(k2, 33);
        k2 = k2.wrapping_mul(c1);
        h2 ^= k2;
    }
    if tail_len >= 8 {
        k1 ^= (tail[7] as u64) << 56;
    }
    if tail_len >= 7 {
        k1 ^= (tail[6] as u64) << 48;
    }
    if tail_len >= 6 {
        k1 ^= (tail[5] as u64) << 40;
    }
    if tail_len >= 5 {
        k1 ^= (tail[4] as u64) << 32;
    }
    if tail_len >= 4 {
        k1 ^= (tail[3] as u64) << 24;
    }
    if tail_len >= 3 {
        k1 ^= (tail[2] as u64) << 16;
    }
    if tail_len >= 2 {
        k1 ^= (tail[1] as u64) << 8;
    }
    if tail_len >= 1 {
        k1 ^= tail[0] as u64;
        k1 = k1.wrapping_mul(c1);
        k1 = rotl64(k1, 31);
        k1 = k1.wrapping_mul(c2);
        h1 ^= k1;
    }

    // Finalization
    h1 ^= len as u64;
    h2 ^= len as u64;

    h1 = h1.wrapping_add(h2);
    h2 = h2.wrapping_add(h1);

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 = h1.wrapping_add(h2);
    h2 = h2.wrapping_add(h1);

    let mut out = [0u8; 16];
    out[0..8].copy_from_slice(&h1.to_le_bytes());
    out[8..16].copy_from_slice(&h2.to_le_bytes());
    out
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

    u64::from_le_bytes(seed[0..8].try_into().unwrap())
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
