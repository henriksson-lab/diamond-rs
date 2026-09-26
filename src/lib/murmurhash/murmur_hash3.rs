//! Rust translation of Austin Appleby's public-domain `MurmurHash3.cpp`.

#[inline]
pub fn rotl32(x: u32, r: u32) -> u32 {
    x.rotate_left(r)
}
#[inline]
pub fn rotl64(x: u64, r: u32) -> u64 {
    x.rotate_left(r)
}
#[inline]
pub fn getblock32(bytes: &[u8]) -> u32 {
    u32::from_ne_bytes(bytes.try_into().expect("four-byte block"))
}
#[inline]
pub fn getblock64(bytes: &[u8]) -> u64 {
    u64::from_le_bytes(bytes.try_into().expect("eight-byte block"))
}
#[inline]
pub fn fmix32(mut h: u32) -> u32 {
    h ^= h >> 16;
    h = h.wrapping_mul(0x85ebca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2ae35);
    h ^ (h >> 16)
}
#[inline]
pub fn fmix64(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^ (h >> 33)
}

/// Active when upstream is built without `WITH_BLASTDB`.
pub fn murmur_hash3_x86_32(key: &[u8], seed: u32) -> u32 {
    let mut h = seed;
    let (c1, c2) = (0xcc9e2d51u32, 0x1b873593u32);
    let n = key.len() / 4;
    for block in key[..n * 4].chunks_exact(4) {
        let mut k = getblock32(block);
        k = k.wrapping_mul(c1);
        k = rotl32(k, 15);
        k = k.wrapping_mul(c2);
        h ^= k;
        h = rotl32(h, 13);
        h = h.wrapping_mul(5).wrapping_add(0xe6546b64);
    }
    let tail = &key[n * 4..];
    let mut k = 0u32;
    if tail.len() >= 3 {
        k ^= (tail[2] as u32) << 16
    }
    if tail.len() >= 2 {
        k ^= (tail[1] as u32) << 8
    }
    if !tail.is_empty() {
        k ^= tail[0] as u32;
        k = k.wrapping_mul(c1);
        k = rotl32(k, 15);
        k = k.wrapping_mul(c2);
        h ^= k
    }
    fmix32(h ^ (key.len() as u32))
}

/// Active when upstream is built without `WITH_BLASTDB`.
pub fn murmur_hash3_x86_128(key: &[u8], seed: u32) -> [u32; 4] {
    let (mut h1, mut h2, mut h3, mut h4) = (seed, seed, seed, seed);
    let (c1, c2, c3, c4) = (0x239b961bu32, 0xab0e9789u32, 0x38b34ae5u32, 0xa1e38b93u32);
    let n = key.len() / 16;
    for b in key[..n * 16].chunks_exact(16) {
        let (mut k1, mut k2, mut k3, mut k4) = (
            getblock32(&b[0..4]),
            getblock32(&b[4..8]),
            getblock32(&b[8..12]),
            getblock32(&b[12..16]),
        );
        k1 = k1.wrapping_mul(c1);
        k1 = rotl32(k1, 15).wrapping_mul(c2);
        h1 ^= k1;
        h1 = rotl32(h1, 19)
            .wrapping_add(h2)
            .wrapping_mul(5)
            .wrapping_add(0x561ccd1b);
        k2 = k2.wrapping_mul(c2);
        k2 = rotl32(k2, 16).wrapping_mul(c3);
        h2 ^= k2;
        h2 = rotl32(h2, 17)
            .wrapping_add(h3)
            .wrapping_mul(5)
            .wrapping_add(0x0bcaa747);
        k3 = k3.wrapping_mul(c3);
        k3 = rotl32(k3, 17).wrapping_mul(c4);
        h3 ^= k3;
        h3 = rotl32(h3, 15)
            .wrapping_add(h4)
            .wrapping_mul(5)
            .wrapping_add(0x96cd1c35);
        k4 = k4.wrapping_mul(c4);
        k4 = rotl32(k4, 18).wrapping_mul(c1);
        h4 ^= k4;
        h4 = rotl32(h4, 13)
            .wrapping_add(h1)
            .wrapping_mul(5)
            .wrapping_add(0x32ac3b17);
    }
    let t = &key[n * 16..];
    let (mut k1, mut k2, mut k3, mut k4) = (0u32, 0u32, 0u32, 0u32);
    for (i, &v) in t.iter().enumerate() {
        match i {
            0..=3 => k1 ^= (v as u32) << (8 * i),
            4..=7 => k2 ^= (v as u32) << (8 * (i - 4)),
            8..=11 => k3 ^= (v as u32) << (8 * (i - 8)),
            _ => k4 ^= (v as u32) << (8 * (i - 12)),
        }
    }
    if t.len() > 12 {
        k4 = k4.wrapping_mul(c4);
        k4 = rotl32(k4, 18).wrapping_mul(c1);
        h4 ^= k4
    }
    if t.len() > 8 {
        k3 = k3.wrapping_mul(c3);
        k3 = rotl32(k3, 17).wrapping_mul(c4);
        h3 ^= k3
    }
    if t.len() > 4 {
        k2 = k2.wrapping_mul(c2);
        k2 = rotl32(k2, 16).wrapping_mul(c3);
        h2 ^= k2
    }
    if !t.is_empty() {
        k1 = k1.wrapping_mul(c1);
        k1 = rotl32(k1, 15).wrapping_mul(c2);
        h1 ^= k1
    }
    let l = key.len() as u32;
    h1 ^= l;
    h2 ^= l;
    h3 ^= l;
    h4 ^= l;
    h1 = h1.wrapping_add(h2).wrapping_add(h3).wrapping_add(h4);
    h2 = h2.wrapping_add(h1);
    h3 = h3.wrapping_add(h1);
    h4 = h4.wrapping_add(h1);
    h1 = fmix32(h1);
    h2 = fmix32(h2);
    h3 = fmix32(h3);
    h4 = fmix32(h4);
    h1 = h1.wrapping_add(h2).wrapping_add(h3).wrapping_add(h4);
    h2 = h2.wrapping_add(h1);
    h3 = h3.wrapping_add(h1);
    h4 = h4.wrapping_add(h1);
    [h1, h2, h3, h4]
}

pub fn murmur_hash3_x64_128(key: &[u8], seed: [u64; 2]) -> [u64; 2] {
    let (mut h1, mut h2) = (seed[0], seed[1]);
    let (c1, c2) = (0x87c37b91114253d5u64, 0x4cf5ad432745937fu64);
    let n = key.len() / 16;
    for b in key[..n * 16].chunks_exact(16) {
        let (mut k1, mut k2) = (getblock64(&b[..8]), getblock64(&b[8..]));
        k1 = k1.wrapping_mul(c1);
        k1 = rotl64(k1, 31).wrapping_mul(c2);
        h1 ^= k1;
        h1 = rotl64(h1, 27)
            .wrapping_add(h2)
            .wrapping_mul(5)
            .wrapping_add(0x52dce729);
        k2 = k2.wrapping_mul(c2);
        k2 = rotl64(k2, 33).wrapping_mul(c1);
        h2 ^= k2;
        h2 = rotl64(h2, 31)
            .wrapping_add(h1)
            .wrapping_mul(5)
            .wrapping_add(0x38495ab5);
    }
    let t = &key[n * 16..];
    let (mut k1, mut k2) = (0u64, 0u64);
    for (i, &v) in t.iter().enumerate() {
        if i < 8 {
            k1 ^= (v as u64) << (8 * i)
        } else {
            k2 ^= (v as u64) << (8 * (i - 8))
        }
    }
    if t.len() > 8 {
        k2 = k2.wrapping_mul(c2);
        k2 = rotl64(k2, 33).wrapping_mul(c1);
        h2 ^= k2
    }
    if !t.is_empty() {
        k1 = k1.wrapping_mul(c1);
        k1 = rotl64(k1, 31).wrapping_mul(c2);
        h1 ^= k1
    }
    let l = key.len() as u64;
    h1 ^= l;
    h2 ^= l;
    h1 = h1.wrapping_add(h2);
    h2 = h2.wrapping_add(h1);
    h1 = fmix64(h1);
    h2 = fmix64(h2);
    h1 = h1.wrapping_add(h2);
    h2 = h2.wrapping_add(h1);
    [h1, h2]
}

pub fn murmur_hash3_x64_128_bytes(key: &[u8], seed: &[u8; 16]) -> [u8; 16] {
    let words = murmur_hash3_x64_128(
        key,
        [
            u64::from_ne_bytes(seed[..8].try_into().unwrap()),
            u64::from_ne_bytes(seed[8..].try_into().unwrap()),
        ],
    );
    let mut out = [0; 16];
    out[..8].copy_from_slice(&words[0].to_ne_bytes());
    out[8..].copy_from_slice(&words[1].to_ne_bytes());
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_reference_boundary_vectors() {
        let data: Vec<u8> = (0..48).map(|x| (x * 37 + 11) as u8).collect();
        let lens = [0, 1, 2, 3, 4, 7, 8, 9, 15, 16, 17, 31, 32, 33];
        let x32 = [
            0xebb6c228, 0x24cb957f, 0x299eeeea, 0x3eb13daa, 0x297b1cbc, 0x7980062b, 0xe8b63cf8,
            0x43ed2ad0, 0xb81b61fa, 0xb62dba30, 0x60c29752, 0x6882ea05, 0x025d728c, 0x0bf85430,
        ];
        let x128 = [
            [0xf7bed5a1, 0x5b576a1c, 0x5b576a1c, 0x5b576a1c],
            [0xd6a16909, 0x4ee55938, 0x4ee55938, 0x4ee55938],
            [0x1bd1769d, 0x4b0848e7, 0x4b0848e7, 0x4b0848e7],
            [0xe8806faf, 0xf705cccc, 0xf705cccc, 0xf705cccc],
            [0x093a16af, 0x0b22b3bc, 0x0b22b3bc, 0x0b22b3bc],
            [0x7038072b, 0xe6ddeef3, 0x23942965, 0x23942965],
            [0xffe8032a, 0x3d6a054a, 0xe0115e06, 0xe0115e06],
            [0xf075d7e7, 0x5e1d8b6b, 0xdd3c2c51, 0xc11c3357],
            [0xaa1b5673, 0xecd59e00, 0xc9971fc3, 0xf262c02b],
            [0xdb6a35c8, 0xd046fb59, 0xb49a54be, 0xb06ae0b7],
            [0xbe221828, 0x1a08b713, 0x5680fea5, 0xe3e11f77],
            [0x5bda3a36, 0x1afdf7c1, 0x904d8b87, 0x71769725],
            [0xba40805a, 0x64a5963b, 0x6fdd0cb0, 0xc517c6b0],
            [0x06743436, 0xc8fd8917, 0x307afdd1, 0x90317974],
        ];
        let x64 = [
            [0x03bc00795ad0f097, 0xa2c28ee76a1f820d],
            [0x5326a2889b5df090, 0x6b1978d33bde918a],
            [0x633f5b086a81d913, 0x992e66ca51f3456f],
            [0x23d0b050af8cb8a4, 0x1be98ce9e7f35b79],
            [0x2631e97ba99dc8b3, 0xae52d38b42a7e2d4],
            [0x88710f24e013f6bd, 0x44e3132a64d18626],
            [0xfa9efd933b6a35b9, 0xecabba9c1a5c4e67],
            [0xbee8ebe0a2ab0b2a, 0x653f70377e48c080],
            [0x3d18c592d8ae4bb9, 0x5f72e818c13db5c0],
            [0xccbd83b9d50f51f6, 0xbb06022209f8866c],
            [0x50466922e9ee9a3a, 0x53aeaac03d95c3f3],
            [0x53c5ab77cfb31cf6, 0xd38ef82c4dfc85f0],
            [0xa3a726a188b47921, 0x156426bcb351c96e],
            [0x6505816eedaf4566, 0x5e796e8117c625fd],
        ];
        for (i, &len) in lens.iter().enumerate() {
            assert_eq!(
                murmur_hash3_x86_32(&data[..len], 0x9747b28c),
                x32[i],
                "x86_32 len={len}"
            );
            assert_eq!(
                murmur_hash3_x86_128(&data[..len], 0x9747b28c),
                x128[i],
                "x86_128 len={len}"
            );
            assert_eq!(
                murmur_hash3_x64_128(&data[..len], [0x0123456789abcdef, 0xfedcba9876543210]),
                x64[i],
                "x64_128 len={len}"
            );
        }
    }
    #[test]
    fn empty_and_seed_endian_surface() {
        assert_eq!(murmur_hash3_x86_32(&[], 0), 0);
        assert_eq!(murmur_hash3_x86_128(&[], 0), [0; 4]);
        assert_eq!(murmur_hash3_x64_128(&[], [0, 0]), [0, 0]);
        let seed = [1, 2];
        assert_ne!(
            murmur_hash3_x64_128(b"x", seed),
            murmur_hash3_x64_128(b"x", [2, 1])
        );
    }
}
