use crate::basic::reduction::Reduction;
use crate::basic::value::{is_amino_acid, letter_mask, Letter, SEED_MASK};

/// Matches C++ `reduce_seq_generic(const Letter*, const Letter*)`.
pub fn reduce_seq_generic(seq: &[Letter], map: &[Letter]) -> [Letter; 16] {
    let mut d = [0; 16];
    for i in 0..16 {
        d[i] = map[letter_mask(seq[i]) as usize];
    }
    d
}

/// Matches C++ `reduce_seq(const Letter*, const Letter*)`.
pub fn reduce_seq(seq: &[Letter], map: &[Letter]) -> [Letter; 16] {
    assert!(seq.len() >= 16);
    assert!(map.len() >= 32);
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("ssse3") {
        // SAFETY: SSSE3 was detected and the lengths were checked above.
        return unsafe { reduce_seq_ssse3(seq, map) };
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is mandatory on AArch64 and both 16-byte
        // sequence/map table reads were checked above.
        unsafe { reduce_seq_neon(seq, map) }
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        reduce_seq_generic(seq, map)
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn reduce_seq_ssse3(seq: &[Letter], map: &[Letter]) -> [Letter; 16] {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let letters = _mm_and_si128(
        _mm_loadu_si128(seq.as_ptr().cast()),
        _mm_set1_epi8(crate::basic::value::LETTER_MASK),
    );
    let high_mask = _mm_slli_epi16(_mm_and_si128(letters, _mm_set1_epi8(0x10)), 3);
    let low_indices = _mm_or_si128(letters, high_mask);
    let high_indices = _mm_or_si128(letters, _mm_xor_si128(high_mask, _mm_set1_epi8(i8::MIN)));
    let low = _mm_shuffle_epi8(_mm_loadu_si128(map.as_ptr().cast()), low_indices);
    let high = _mm_shuffle_epi8(_mm_loadu_si128(map.as_ptr().add(16).cast()), high_indices);
    let reduced = _mm_or_si128(low, high);
    let mut out = [0; 16];
    _mm_storeu_si128(out.as_mut_ptr().cast(), reduced);
    out
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn reduce_seq_neon(seq: &[Letter], map: &[Letter]) -> [Letter; 16] {
    use std::arch::aarch64::*;

    let letters = vandq_s8(
        vld1q_s8(seq.as_ptr()),
        vdupq_n_s8(crate::basic::value::LETTER_MASK),
    );
    let indices = vreinterpretq_u8_s8(vandq_s8(letters, vdupq_n_s8(0x0f)));
    let low = vqtbl1q_s8(vld1q_s8(map.as_ptr()), indices);
    let high = vqtbl1q_s8(vld1q_s8(map.as_ptr().add(16)), indices);
    let reduced = vbslq_s8(vcgeq_s8(letters, vdupq_n_s8(16)), high, low);
    let mut out = [0; 16];
    vst1q_s8(out.as_mut_ptr(), reduced);
    out
}

/// Matches C++ `match_block_reduced(const Letter*, const Letter*, const Reduction&)`.
pub fn match_block_reduced(x: &[Letter], y: &[Letter], reduction: &Reduction) -> u32 {
    assert!(x.len() >= 16 && y.len() >= 16);
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("ssse3") {
        // SAFETY: SSSE3 was detected and both 16-byte reads are in bounds.
        return unsafe { match_block_reduced_ssse3(x, y, reduction) };
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is mandatory on AArch64 and the input lengths
        // guarantee both vector loads are in bounds.
        unsafe { match_block_reduced_neon(x, y, reduction) }
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        match_block_reduced_partial(x, y, 16, reduction)
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn match_block_reduced_ssse3(x: &[Letter], y: &[Letter], reduction: &Reduction) -> u32 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let left = reduce_seq_ssse3(x, reduction.map8());
    let right = reduce_seq_ssse3(y, reduction.map8b());
    let left = _mm_loadu_si128(left.as_ptr().cast());
    let right = _mm_loadu_si128(right.as_ptr().cast());
    _mm_movemask_epi8(_mm_cmpeq_epi8(left, right)) as u32
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn movemask_neon(value: std::arch::aarch64::int8x16_t) -> u16 {
    use std::arch::aarch64::*;

    // Extract the sign bit from every byte, then horizontally sum weighted
    // halves. Each half sums to at most 255, so the u8 reductions are exact.
    let sign = vreinterpretq_u8_s8(vshrq_n_s8(value, 7));
    let weights_data = [1u8, 2, 4, 8, 16, 32, 64, 128];
    let weights = vld1_u8(weights_data.as_ptr());
    let low = vaddv_u8(vand_u8(vget_low_u8(sign), weights)) as u16;
    let high = vaddv_u8(vand_u8(vget_high_u8(sign), weights)) as u16;
    low | (high << 8)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn match_block_reduced_neon(x: &[Letter], y: &[Letter], reduction: &Reduction) -> u32 {
    use std::arch::aarch64::*;

    let left = reduce_seq_neon(x, reduction.map8());
    let right = reduce_seq_neon(y, reduction.map8b());
    let equal = vceqq_s8(vld1q_s8(left.as_ptr()), vld1q_s8(right.as_ptr()));
    u32::from(movemask_neon(vreinterpretq_s8_u8(equal)))
}

fn match_block_reduced_partial(x: &[Letter], y: &[Letter], n: usize, reduction: &Reduction) -> u32 {
    let mut r = 0u32;
    for i in (0..n).rev() {
        r <<= 1;
        let lx = letter_mask(x[i]);
        let ly = letter_mask(y[i]);
        if !is_amino_acid(lx) || !is_amino_acid(ly) {
            continue;
        }
        if reduction.reduce(lx) == reduction.reduce(ly) {
            r |= 1;
        }
    }
    r
}

#[inline]
fn match_block_reduced_prefix(x: &[Letter], y: &[Letter], n: usize, reduction: &Reduction) -> u32 {
    if n == 16 {
        match_block_reduced(x, y, reduction)
    } else {
        match_block_reduced_partial(x, y, n, reduction)
    }
}

/// Matches C++ `reduced_match32(const Letter*, const Letter*, unsigned, const Reduction&)`.
pub fn reduced_match32(q: &[Letter], s: &[Letter], len: u32, reduction: &Reduction) -> u64 {
    let len = len as usize;
    let mut x = match_block_reduced_prefix(q, s, len.min(16), reduction) as u64;
    if len > 16 {
        x |= (match_block_reduced_prefix(&q[16..], &s[16..], (len - 16).min(16), reduction) as u64)
            << 16;
    }
    if len < 32 {
        x &= (1u64 << len) - 1;
    }
    x
}

/// Matches C++ `reduced_match(const Letter*, const Letter*, int, const Reduction&)`.
#[inline(always)]
pub fn reduced_match(q: &[Letter], s: &[Letter], len: i32, reduction: &Reduction) -> u64 {
    assert!(len <= 64);
    let len = len as usize;
    if len < 64 {
        let mask = (1u64 << len) - 1;
        let mut m = match_block_reduced_prefix(q, s, len.min(16), reduction) as u64;
        if len <= 16 {
            return m & mask;
        }
        m |= (match_block_reduced_prefix(&q[16..], &s[16..], (len - 16).min(16), reduction) as u64)
            << 16;
        if len <= 32 {
            return m & mask;
        }
        m |= (match_block_reduced_prefix(&q[32..], &s[32..], (len - 32).min(16), reduction) as u64)
            << 32;
        if len <= 48 {
            return m & mask;
        }
        m |= (match_block_reduced_prefix(&q[48..], &s[48..], len - 48, reduction) as u64) << 48;
        m & mask
    } else {
        match_block_reduced(q, s, reduction) as u64
            | ((match_block_reduced(&q[16..], &s[16..], reduction) as u64) << 16)
            | ((match_block_reduced(&q[32..], &s[32..], reduction) as u64) << 32)
            | ((match_block_reduced(&q[48..], &s[48..], reduction) as u64) << 48)
    }
}

/// Matches C++ `seed_mask(const Letter*, int)`.
#[inline(always)]
pub fn seed_mask(s: &[Letter], len: i32) -> u64 {
    assert!(len <= 64);
    assert!(len >= 0 && s.len() >= len as usize);
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("sse2") {
        // SAFETY: SSE2 is present and the kernel loads only complete chunks.
        return unsafe { seed_mask_sse2(s, len as usize) };
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is mandatory on AArch64; the kernel handles
        // only complete vector chunks and uses scalar code for the tail.
        return unsafe { seed_mask_neon(s, len as usize) };
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        let mut mask = 0u64;
        for i in 0..len as usize {
            if (s[i] & SEED_MASK) != 0 {
                mask |= 1u64 << i;
            }
        }
        mask
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn seed_mask_sse2(s: &[Letter], len: usize) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let mut mask = 0u64;
    let vector_end = len / 16 * 16;
    let seed_bit = _mm_set1_epi8(SEED_MASK);
    for offset in (0..vector_end).step_by(16) {
        let letters = _mm_loadu_si128(s.as_ptr().add(offset).cast());
        let bits = _mm_and_si128(letters, seed_bit);
        mask |= (_mm_movemask_epi8(bits) as u64) << offset;
    }
    for (offset, &letter) in s[vector_end..len].iter().enumerate() {
        if letter & SEED_MASK != 0 {
            mask |= 1u64 << (vector_end + offset);
        }
    }
    mask
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn seed_mask_neon(s: &[Letter], len: usize) -> u64 {
    use std::arch::aarch64::*;

    let mut mask = 0u64;
    let vector_end = len / 16 * 16;
    let seed_bit = vdupq_n_s8(SEED_MASK);
    for offset in (0..vector_end).step_by(16) {
        let letters = vld1q_s8(s.as_ptr().add(offset));
        let bits = vandq_s8(letters, seed_bit);
        mask |= u64::from(movemask_neon(bits)) << offset;
    }
    for (offset, &letter) in s[vector_end..len].iter().enumerate() {
        if letter & SEED_MASK != 0 {
            mask |= 1u64 << (vector_end + offset);
        }
    }
    mask
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::{MASK_LETTER, STOP_LETTER};

    #[test]
    fn test_match_block_reduced() {
        let reduction = Reduction::default_reduction();
        let x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        let mut y = x;
        let reduced = reduce_seq(&x, reduction.map8());
        for i in 0..16 {
            assert_eq!(reduced[i], reduction.map8()[x[i] as usize]);
        }
        assert_eq!(match_block_reduced(&x, &y, &reduction), 0xffff);
        y[0] = 1;
        y[1] = 11;
        y[2] = STOP_LETTER;
        let mask = match_block_reduced(&x, &y, &reduction);
        assert_eq!(mask & 1, 0);
        assert_eq!(mask & 2, 2);
        assert_eq!(mask & 4, 0);
    }

    #[test]
    fn test_reduced_match_and_seed_mask() {
        let reduction = Reduction::default_reduction();
        let q: Vec<Letter> = (0..64).map(|i| (i % 20) as Letter).collect();
        let mut s = q.clone();
        s[3] = MASK_LETTER;
        let m = reduced_match(&q, &s, 20, &reduction);
        assert_eq!(m & (1 << 0), 1);
        assert_eq!(m & (1 << 3), 0);
        assert_eq!(m >> 20, 0);
        assert_eq!(reduced_match32(&q, &s, 20, &reduction), m);

        let mut masked = vec![0i8; 64];
        masked[2] = SEED_MASK;
        masked[9] = SEED_MASK | 1;
        assert_eq!(seed_mask(&masked, 12), (1 << 2) | (1 << 9));
    }

    #[test]
    fn reduced_match_all_lengths_matches_scalar_with_masked_and_ambiguous_letters() {
        let reduction = Reduction::default_reduction();
        let q: Vec<Letter> = (0..64)
            .map(|i| {
                let letter = match i % 9 {
                    0 => STOP_LETTER,
                    1 => MASK_LETTER,
                    _ => ((i * 7 + 3) % 20) as Letter,
                };
                letter | if i % 5 == 0 { SEED_MASK } else { 0 }
            })
            .collect();
        let s: Vec<Letter> = (0..64)
            .map(|i| {
                let letter = match i % 11 {
                    0 => STOP_LETTER,
                    1 => MASK_LETTER,
                    _ => ((i * 13 + 1) % 20) as Letter,
                };
                letter | if i % 7 == 0 { SEED_MASK } else { 0 }
            })
            .collect();

        for len in 0..=64usize {
            let mut expected = 0u64;
            for offset in (0..len).step_by(16) {
                let width = (len - offset).min(16);
                expected |= u64::from(match_block_reduced_partial(
                    &q[offset..],
                    &s[offset..],
                    width,
                    &reduction,
                )) << offset;
            }
            assert_eq!(reduced_match(&q, &s, len as i32, &reduction), expected);
            if len <= 32 {
                assert_eq!(reduced_match32(&q, &s, len as u32, &reduction), expected);
            }
        }
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn x86_kernels_match_scalar_for_all_letter_codes() {
        if !std::arch::is_x86_feature_detected!("ssse3") {
            return;
        }
        let reduction = Reduction::default_reduction();
        for shift in 0..32 {
            let x: Vec<Letter> = (0..16)
                .map(|i| ((i + shift) % 32) as Letter | if i % 3 == 0 { SEED_MASK } else { 0 })
                .collect();
            let y: Vec<Letter> = (0..16)
                .map(|i| ((i * 7 + shift) % 32) as Letter | if i % 5 == 0 { SEED_MASK } else { 0 })
                .collect();
            assert_eq!(
                unsafe { match_block_reduced_ssse3(&x, &y, &reduction) },
                match_block_reduced_partial(&x, &y, 16, &reduction)
            );
            assert_eq!(
                unsafe { reduce_seq_ssse3(&x, reduction.map8()) },
                reduce_seq_generic(&x, reduction.map8())
            );
        }

        if std::arch::is_x86_feature_detected!("sse2") {
            let letters: Vec<Letter> = (0..64)
                .map(|i| {
                    if i % 3 == 0 || i % 11 == 0 {
                        SEED_MASK
                    } else {
                        4
                    }
                })
                .collect();
            for len in 0..=64 {
                let expected = letters[..len]
                    .iter()
                    .enumerate()
                    .fold(0u64, |m, (i, &x)| m | (((x & SEED_MASK != 0) as u64) << i));
                assert_eq!(unsafe { seed_mask_sse2(&letters, len) }, expected);
            }
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn neon_kernels_match_scalar_for_all_letter_codes_and_lengths() {
        let reduction = Reduction::default_reduction();
        for shift in 0..32 {
            let x: Vec<Letter> = (0..16)
                .map(|i| ((i + shift) % 32) as Letter | if i % 3 == 0 { SEED_MASK } else { 0 })
                .collect();
            let y: Vec<Letter> = (0..16)
                .map(|i| ((i * 7 + shift) % 32) as Letter | if i % 5 == 0 { SEED_MASK } else { 0 })
                .collect();
            assert_eq!(
                unsafe { match_block_reduced_neon(&x, &y, &reduction) },
                match_block_reduced_partial(&x, &y, 16, &reduction)
            );
            assert_eq!(
                unsafe { reduce_seq_neon(&x, reduction.map8()) },
                reduce_seq_generic(&x, reduction.map8())
            );
        }

        let letters: Vec<Letter> = (0..64)
            .map(|i| {
                if i % 3 == 0 || i % 11 == 0 {
                    SEED_MASK
                } else {
                    4
                }
            })
            .collect();
        for len in 0..=64 {
            let expected = letters[..len]
                .iter()
                .enumerate()
                .fold(0u64, |m, (i, &x)| m | (((x & SEED_MASK != 0) as u64) << i));
            assert_eq!(unsafe { seed_mask_neon(&letters, len) }, expected);
        }
    }
}
