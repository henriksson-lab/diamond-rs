//! Allocation-free 48-byte Hamming tile kernel.
//!
//! This translates `search/hamming/kernel.h::all_vs_all`, including its
//! four-query unroll and target-load reuse.

/// A DIAMOND Hamming fingerprint with the alignment required by the x86
/// implementations. Only the first 48 bytes participate in comparisons.
#[repr(C, align(16))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AlignedFingerprint48 {
    pub bytes: [i8; 48],
}

impl AlignedFingerprint48 {
    #[inline]
    pub const fn new(bytes: [i8; 48]) -> Self {
        Self { bytes }
    }
}

/// Compare one query tile with one target tile.
///
/// `pass_masks[q]` receives one bit per target fingerprint. Bit `t` is set
/// exactly when query `q` and target `t` have at least `threshold` equal
/// bytes. Target tiles are limited to 64 entries so the result for each query
/// fits in one word. Existing mask contents are replaced, not appended.
#[cfg(test)]
fn all_vs_all_pass_masks(
    query: &[AlignedFingerprint48],
    target: &[AlignedFingerprint48],
    threshold: u32,
    pass_masks: &mut [u64],
) {
    assert!(target.len() <= 64, "Hamming target tile exceeds 64 entries");
    assert!(
        pass_masks.len() >= query.len(),
        "Hamming output is shorter than the query tile"
    );
    pass_masks[..query.len()].fill(0);
    if query.is_empty() || target.is_empty() || threshold > 48 {
        return;
    }

    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx512bw") {
        // SAFETY: AVX-512BW is runtime-detected and masked loads consume only
        // the 48 initialized fingerprint bytes.
        unsafe { all_vs_all_pass_masks_avx512bw(query, target, threshold, pass_masks) };
        return;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 is runtime-detected. Every load consumes exactly the
        // 48 initialized bytes in an `AlignedFingerprint48`.
        unsafe { all_vs_all_pass_masks_avx2(query, target, threshold, pass_masks) };
        return;
    }

    all_vs_all_scalar(query, target, threshold, pass_masks);
}

/// AVX-512BW translation of upstream's four-query `all_vs_all` loop.
///
/// Four masked query loads remain resident in ZMM registers for the complete
/// target pass. Each 48-byte target is loaded once and compared to all four
/// queries, producing one query-major 64-bit pass mask.
#[cfg(target_arch = "x86_64")]
pub(crate) unsafe fn all_vs_all_pass_masks_avx512bw(
    query: &[AlignedFingerprint48],
    target: &[AlignedFingerprint48],
    threshold: u32,
    pass_masks: &mut [u64],
) {
    assert!(target.len() <= 64, "Hamming target tile exceeds 64 entries");
    assert!(pass_masks.len() >= query.len());
    pass_masks[..query.len()].fill(0);
    if query.is_empty() || target.is_empty() || threshold > 48 {
        return;
    }

    const VALID_48: u64 = (1u64 << 48) - 1;
    let query4 = query.len() & !3;
    let mut query_index = 0usize;
    while query_index < query4 {
        let query_ptr = query.as_ptr().add(query_index);
        let target_ptr = target.as_ptr();
        let remaining = target.len();
        let m0: u64;
        let m1: u64;
        let m2: u64;
        let m3: u64;
        std::arch::asm!(
            "kmovq k1, {valid}",
            "vmovdqu8 zmm0 {{k1}}{{z}}, [{query_ptr}]",
            "vmovdqu8 zmm1 {{k1}}{{z}}, [{query_ptr} + 48]",
            "vmovdqu8 zmm2 {{k1}}{{z}}, [{query_ptr} + 96]",
            "vmovdqu8 zmm3 {{k1}}{{z}}, [{query_ptr} + 144]",
            "xor {m0}, {m0}",
            "xor {m1}, {m1}",
            "xor {m2}, {m2}",
            "xor {m3}, {m3}",
            "xor {bit}, {bit}",
            "2:",
            "vmovdqu8 zmm4 {{k1}}{{z}}, [{target_ptr}]",
            "vpcmpeqb k2 {{k1}}, zmm0, zmm4",
            "vpcmpeqb k3 {{k1}}, zmm1, zmm4",
            "vpcmpeqb k4 {{k1}}, zmm2, zmm4",
            "vpcmpeqb k5 {{k1}}, zmm3, zmm4",
            "kmovq rax, k2",
            "popcnt rax, rax",
            "cmp eax, {threshold:e}",
            "jb 3f",
            "bts {m0}, {bit}",
            "3:",
            "kmovq rax, k3",
            "popcnt rax, rax",
            "cmp eax, {threshold:e}",
            "jb 4f",
            "bts {m1}, {bit}",
            "4:",
            "kmovq rax, k4",
            "popcnt rax, rax",
            "cmp eax, {threshold:e}",
            "jb 5f",
            "bts {m2}, {bit}",
            "5:",
            "kmovq rax, k5",
            "popcnt rax, rax",
            "cmp eax, {threshold:e}",
            "jb 6f",
            "bts {m3}, {bit}",
            "6:",
            "add {target_ptr}, 48",
            "inc {bit}",
            "dec {remaining}",
            "jnz 2b",
            valid = in(reg) VALID_48,
            query_ptr = in(reg) query_ptr,
            target_ptr = inout(reg) target_ptr => _,
            remaining = inout(reg) remaining => _,
            threshold = in(reg) threshold,
            bit = out(reg) _,
            m0 = out(reg) m0,
            m1 = out(reg) m1,
            m2 = out(reg) m2,
            m3 = out(reg) m3,
            out("rax") _,
            out("zmm0") _, out("zmm1") _, out("zmm2") _, out("zmm3") _, out("zmm4") _,
            out("k1") _, out("k2") _, out("k3") _, out("k4") _, out("k5") _,
            options(readonly, nostack),
        );
        pass_masks[query_index] = m0;
        pass_masks[query_index + 1] = m1;
        pass_masks[query_index + 2] = m2;
        pass_masks[query_index + 3] = m3;
        query_index += 4;
    }

    while query_index < query.len() {
        let query_ptr = query.as_ptr().add(query_index);
        let target_ptr = target.as_ptr();
        let remaining = target.len();
        let mask: u64;
        std::arch::asm!(
            "kmovq k1, {valid}",
            "vmovdqu8 zmm0 {{k1}}{{z}}, [{query_ptr}]",
            "xor {mask}, {mask}",
            "xor {bit}, {bit}",
            "2:",
            "vmovdqu8 zmm1 {{k1}}{{z}}, [{target_ptr}]",
            "vpcmpeqb k2 {{k1}}, zmm0, zmm1",
            "kmovq rax, k2",
            "popcnt rax, rax",
            "cmp eax, {threshold:e}",
            "jb 3f",
            "bts {mask}, {bit}",
            "3:",
            "add {target_ptr}, 48",
            "inc {bit}",
            "dec {remaining}",
            "jnz 2b",
            valid = in(reg) VALID_48,
            query_ptr = in(reg) query_ptr,
            target_ptr = inout(reg) target_ptr => _,
            remaining = inout(reg) remaining => _,
            threshold = in(reg) threshold,
            bit = out(reg) _,
            mask = out(reg) mask,
            out("rax") _,
            out("zmm0") _, out("zmm1") _,
            out("k1") _, out("k2") _,
            options(readonly, nostack),
        );
        pass_masks[query_index] = mask;
        query_index += 1;
    }
}

#[cfg(test)]
fn all_vs_all_scalar(
    query: &[AlignedFingerprint48],
    target: &[AlignedFingerprint48],
    threshold: u32,
    pass_masks: &mut [u64],
) {
    for (query_fp, mask) in query.iter().zip(pass_masks.iter_mut()) {
        let mut passed = 0u64;
        for (target_index, target_fp) in target.iter().enumerate() {
            let equal = query_fp
                .bytes
                .iter()
                .zip(&target_fp.bytes)
                .map(|(&q, &t)| u32::from(q == t))
                .sum::<u32>();
            passed |= u64::from(equal >= threshold) << target_index;
        }
        *mask = passed;
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn all_vs_all_pass_masks_avx2(
    query: &[AlignedFingerprint48],
    target: &[AlignedFingerprint48],
    threshold: u32,
    pass_masks: &mut [u64],
) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    #[inline(always)]
    unsafe fn load(fp: &AlignedFingerprint48) -> (__m256i, __m128i) {
        (
            _mm256_loadu_si256(fp.bytes.as_ptr().cast()),
            _mm_load_si128(fp.bytes.as_ptr().add(32).cast()),
        )
    }

    #[inline(always)]
    unsafe fn equal_count(
        query_lo: __m256i,
        query_hi: __m128i,
        target_lo: __m256i,
        target_hi: __m128i,
    ) -> u32 {
        let lo = _mm256_movemask_epi8(_mm256_cmpeq_epi8(query_lo, target_lo)) as u32;
        let hi = _mm_movemask_epi8(_mm_cmpeq_epi8(query_hi, target_hi)) as u32;
        lo.count_ones() + hi.count_ones()
    }

    // Match upstream's four-query unroll: query fingerprints stay resident,
    // and each target fingerprint is loaded once for four comparisons.
    let query4 = query.len() & !3;
    let mut query_index = 0usize;
    while query_index < query4 {
        let (q0_lo, q0_hi) = load(&query[query_index]);
        let (q1_lo, q1_hi) = load(&query[query_index + 1]);
        let (q2_lo, q2_hi) = load(&query[query_index + 2]);
        let (q3_lo, q3_hi) = load(&query[query_index + 3]);
        let mut m0 = 0u64;
        let mut m1 = 0u64;
        let mut m2 = 0u64;
        let mut m3 = 0u64;
        for (target_index, target_fp) in target.iter().enumerate() {
            let (target_lo, target_hi) = load(target_fp);
            let bit = 1u64 << target_index;
            if equal_count(q0_lo, q0_hi, target_lo, target_hi) >= threshold {
                m0 |= bit;
            }
            if equal_count(q1_lo, q1_hi, target_lo, target_hi) >= threshold {
                m1 |= bit;
            }
            if equal_count(q2_lo, q2_hi, target_lo, target_hi) >= threshold {
                m2 |= bit;
            }
            if equal_count(q3_lo, q3_hi, target_lo, target_hi) >= threshold {
                m3 |= bit;
            }
        }
        pass_masks[query_index] = m0;
        pass_masks[query_index + 1] = m1;
        pass_masks[query_index + 2] = m2;
        pass_masks[query_index + 3] = m3;
        query_index += 4;
    }

    while query_index < query.len() {
        let (query_lo, query_hi) = load(&query[query_index]);
        let mut mask = 0u64;
        for (target_index, target_fp) in target.iter().enumerate() {
            let (target_lo, target_hi) = load(target_fp);
            if equal_count(query_lo, query_hi, target_lo, target_hi) >= threshold {
                mask |= 1u64 << target_index;
            }
        }
        pass_masks[query_index] = mask;
        query_index += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_fingerprint(state: &mut u64) -> AlignedFingerprint48 {
        let mut bytes = [0i8; 48];
        for byte in &mut bytes {
            *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            *byte = (*state >> 56) as i8;
        }
        AlignedFingerprint48::new(bytes)
    }

    #[test]
    fn exact_threshold_edges() {
        let same = AlignedFingerprint48::new([7; 48]);
        let mut one_different = same;
        one_different.bytes[31] = 8;
        let query = [same];
        let target = [same, one_different];
        let mut masks = [u64::MAX];
        all_vs_all_pass_masks(&query, &target, 48, &mut masks);
        assert_eq!(masks, [0b01]);
        all_vs_all_pass_masks(&query, &target, 47, &mut masks);
        assert_eq!(masks, [0b11]);
        all_vs_all_pass_masks(&query, &target, 49, &mut masks);
        assert_eq!(masks, [0]);
        all_vs_all_pass_masks(&query, &target, 0, &mut masks);
        assert_eq!(masks, [0b11]);
    }

    #[test]
    fn randomized_simd_kernels_match_scalar_for_all_tile_tails() {
        let mut state = 0xd1a0_0d5e_ed12_3456;
        for query_count in 0..=17 {
            for target_count in [0, 1, 2, 3, 7, 31, 32, 33, 63, 64] {
                let query: Vec<_> = (0..query_count)
                    .map(|_| random_fingerprint(&mut state))
                    .collect();
                let target: Vec<_> = (0..target_count)
                    .map(|_| random_fingerprint(&mut state))
                    .collect();
                // Seed exact and near-exact pairs so masks are not almost
                // always zero for realistic thresholds.
                let mut target = target;
                if !query.is_empty() && !target.is_empty() {
                    target[0] = query[0];
                    if target.len() > 1 {
                        target[1] = query[query.len() - 1];
                        target[1].bytes[47] ^= 1;
                    }
                }
                for threshold in [0, 1, 12, 24, 47, 48, 49] {
                    let mut expected = vec![u64::MAX; query_count];
                    all_vs_all_scalar(&query, &target, threshold, &mut expected);
                    let mut dispatched = vec![u64::MAX; query_count];
                    all_vs_all_pass_masks(&query, &target, threshold, &mut dispatched);
                    assert_eq!(dispatched, expected);

                    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
                    if std::arch::is_x86_feature_detected!("avx2") {
                        let mut avx2 = vec![u64::MAX; query_count];
                        unsafe {
                            all_vs_all_pass_masks_avx2(&query, &target, threshold, &mut avx2)
                        };
                        assert_eq!(avx2, expected);
                    }
                    #[cfg(target_arch = "x86_64")]
                    if std::arch::is_x86_feature_detected!("avx512bw") {
                        let mut avx512 = vec![u64::MAX; query_count];
                        unsafe {
                            all_vs_all_pass_masks_avx512bw(&query, &target, threshold, &mut avx512)
                        };
                        assert_eq!(avx512, expected);
                    }
                }
            }
        }
    }

    #[test]
    fn masks_are_query_major_and_cover_bit_63() {
        let q0 = AlignedFingerprint48::new([1; 48]);
        let q1 = AlignedFingerprint48::new([2; 48]);
        let mut target = vec![AlignedFingerprint48::new([0; 48]); 64];
        target[3] = q1;
        target[63] = q0;
        let mut masks = [0; 2];
        all_vs_all_pass_masks(&[q0, q1], &target, 48, &mut masks);
        assert_eq!(masks[0], 1u64 << 63);
        assert_eq!(masks[1], 1u64 << 3);
    }

    #[test]
    fn fingerprint_layout_matches_upstream_container_stride() {
        assert_eq!(std::mem::size_of::<AlignedFingerprint48>(), 48);
        assert_eq!(std::mem::align_of::<AlignedFingerprint48>(), 16);
    }
}
