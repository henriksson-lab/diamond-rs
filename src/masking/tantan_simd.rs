//! Architecture SIMD tantan forward/backward steps.
//!
//! Direct ports of C++ `masking/tantan.cpp`: AVX2 processes eight floats per
//! register, SSE processes four, and AArch64 NEON processes four.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[cfg(target_arch = "x86_64")]
pub fn has_avx2() -> bool {
    is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma")
}

#[cfg(not(target_arch = "x86_64"))]
pub fn has_avx2() -> bool {
    false
}

/// C++ builds its four-float x86 dispatch target with SSSE3 and SSE4.1.  The
/// tantan arithmetic itself only needs older packed-float operations, but use
/// the same feature gate as that dispatch target so runtime selection mirrors
/// upstream exactly.
#[cfg(target_arch = "x86_64")]
pub fn has_sse41_ssse3() -> bool {
    is_x86_feature_detected!("sse4.1") && is_x86_feature_detected!("ssse3")
}

#[cfg(not(target_arch = "x86_64"))]
pub fn has_sse41_ssse3() -> bool {
    false
}

/// SSE horizontal sum matching `vector8_sse.h::hsum` exactly.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1,ssse3")]
unsafe fn hsum_sse(a: __m128) -> f32 {
    let shuf = _mm_shuffle_ps(a, a, 0b10_11_00_01);
    let sums = _mm_add_ps(a, shuf);
    let shuf = _mm_shuffle_ps(sums, sums, 0b01_00_11_10);
    _mm_cvtss_f32(_mm_add_ss(sums, shuf))
}

/// Four-lane x86 backend corresponding to C++'s SSE4.1/SSSE3 dispatch build.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1,ssse3")]
pub unsafe fn forward_step_sse(
    f: &mut [f32; 50],
    d: &[f32; 50],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
    f_sum_prev: f32,
) -> f32 {
    let b_old = *b;
    let vf2f = _mm_set1_ps(f2f);
    let vb_old = _mm_set1_ps(b_old);
    let mut f_sum_new = 0.0f32;
    for off in (0..48).step_by(4) {
        let vf = _mm_loadu_ps(f.as_ptr().add(off));
        let vd = _mm_loadu_ps(d.as_ptr().add(off));
        let ve = _mm_loadu_ps(e_seg.as_ptr().add(off));
        let tmp = _mm_add_ps(_mm_mul_ps(vf, vf2f), _mm_mul_ps(vb_old, vd));
        let vf_new = _mm_mul_ps(tmp, ve);
        _mm_storeu_ps(f.as_mut_ptr().add(off), vf_new);
        f_sum_new += hsum_sse(vf_new);
    }
    for off in 48..50 {
        let vf = (f[off] * f2f + b_old * d[off]) * e_seg[off];
        f[off] = vf;
        f_sum_new += vf;
    }
    *b = b_old * b2b + f_sum_prev * p_repeat_end;
    f_sum_new
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1,ssse3")]
pub unsafe fn backward_step_sse(
    f: &mut [f32; 50],
    d: &[f32; 50],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
) {
    let vf2f = _mm_set1_ps(f2f);
    let vc = _mm_set1_ps(p_repeat_end * *b);
    let mut tsum = 0.0f32;
    for off in (0..48).step_by(4) {
        let vf = _mm_mul_ps(
            _mm_loadu_ps(f.as_ptr().add(off)),
            _mm_loadu_ps(e_seg.as_ptr().add(off)),
        );
        tsum += hsum_sse(_mm_mul_ps(vf, _mm_loadu_ps(d.as_ptr().add(off))));
        _mm_storeu_ps(
            f.as_mut_ptr().add(off),
            _mm_add_ps(_mm_mul_ps(vf, vf2f), vc),
        );
    }
    for off in 48..50 {
        let vf = f[off] * e_seg[off];
        tsum += vf * d[off];
        f[off] = vf * f2f + p_repeat_end * *b;
    }
    *b = b2b * *b + tsum;
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1,ssse3")]
pub unsafe fn scale_sse(f: &mut [f32; 50], s: f32) {
    let vs = _mm_set1_ps(s);
    for off in (0..48).step_by(4) {
        let value = _mm_loadu_ps(f.as_ptr().add(off));
        _mm_storeu_ps(f.as_mut_ptr().add(off), _mm_mul_ps(value, vs));
    }
    f[48] *= s;
    f[49] *= s;
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1,ssse3")]
pub unsafe fn sum_sse(f: &[f32; 50]) -> f32 {
    let mut acc = _mm_setzero_ps();
    for off in (0..48).step_by(4) {
        acc = _mm_add_ps(acc, _mm_loadu_ps(f.as_ptr().add(off)));
    }
    hsum_sse(acc) + f[48] + f[49]
}

#[cfg(all(test, target_arch = "x86_64"))]
mod sse_tests {
    use super::*;

    fn close(left: f32, right: f32) {
        let tolerance = 2.0e-5 * left.abs().max(right.abs()).max(1.0);
        assert!((left - right).abs() <= tolerance, "{left} != {right}");
    }

    #[test]
    fn sse_forward_backward_scale_and_sum_match_scalar() {
        if !has_sse41_ssse3() {
            return;
        }
        let original = std::array::from_fn(|i| (i as f32 + 1.0) / 97.0);
        let d = std::array::from_fn(|i| (50 - i) as f32 / 701.0);
        let e: [f32; 50] = std::array::from_fn(|i| 0.73 + i as f32 / 211.0);
        let f2f = 0.951;
        let p_repeat_end = 0.049;
        let b2b = 0.993;

        let mut simd = original;
        let mut scalar = original;
        let mut simd_b = 0.625;
        let mut scalar_b = simd_b;
        let simd_sum = unsafe {
            forward_step_sse(&mut simd, &d, &e, &mut simd_b, f2f, p_repeat_end, b2b, 1.25)
        };
        let mut scalar_sum = 0.0;
        for i in 0..50 {
            scalar[i] = (scalar[i] * f2f + scalar_b * d[i]) * e[i];
            scalar_sum += scalar[i];
        }
        scalar_b = scalar_b * b2b + 1.25 * p_repeat_end;
        assert_eq!(simd, scalar);
        close(simd_sum, scalar_sum);
        assert_eq!(simd_b.to_bits(), scalar_b.to_bits());

        simd = original;
        scalar = original;
        simd_b = 0.625;
        scalar_b = simd_b;
        unsafe { backward_step_sse(&mut simd, &d, &e, &mut simd_b, f2f, p_repeat_end, b2b) };
        let mut tsum = 0.0;
        for i in 0..50 {
            let value = scalar[i] * e[i];
            tsum += value * d[i];
            scalar[i] = value * f2f + p_repeat_end * scalar_b;
        }
        scalar_b = b2b * scalar_b + tsum;
        assert_eq!(simd, scalar);
        close(simd_b, scalar_b);

        simd = original;
        unsafe { scale_sse(&mut simd, 1.75) };
        for i in 0..50 {
            assert_eq!(simd[i].to_bits(), (original[i] * 1.75).to_bits());
        }
        close(unsafe { sum_sse(&original) }, original.iter().sum());
    }
}

/// AVX2 horizontal sum: matches C++ hsum(__m256 a) exactly.
///   1. Split into two 128-bit halves, add them
///   2. Two horizontal adds
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn hsum_avx2(a: __m256) -> f32 {
    let vlow = _mm256_castps256_ps128(a);
    let vhigh = _mm256_extractf128_ps(a, 1);
    let vsum = _mm_add_ps(vlow, vhigh);
    let vsum = _mm_hadd_ps(vsum, vsum);
    let vsum = _mm_hadd_ps(vsum, vsum);
    _mm_cvtss_f32(vsum)
}

/// AVX2 forward step: matches C++ forward_step() exactly.
///
/// Processes f[0..48] with AVX2 (6 chunks of 8), then f[48..50] scalar.
///
/// Upstream calls `SIMD::fmadd` here. Its AVX2 dispatch object is compiled with
/// `-march=native`, which selects `_mm256_fmadd_ps` on an FMA-capable host.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
pub unsafe fn forward_step_avx2(
    f: &mut [f32; 50],
    d: &[f32; 50],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
    f_sum_prev: f32,
) -> f32 {
    let b_old = *b;
    let vf2f = _mm256_set1_ps(f2f);
    let vb_old = _mm256_set1_ps(b_old);
    let mut f_sum_new = 0.0f32;

    // Process 48 elements in 6 SIMD chunks
    for off in (0..48).step_by(8) {
        let vf = _mm256_loadu_ps(f.as_ptr().add(off));
        let vd = _mm256_loadu_ps(d.as_ptr().add(off));
        let ve = _mm256_loadu_ps(e_seg.as_ptr().add(off));
        let tmp = _mm256_fmadd_ps(vf, vf2f, _mm256_mul_ps(vb_old, vd));
        let vf_new = _mm256_mul_ps(tmp, ve);
        _mm256_storeu_ps(f.as_mut_ptr().add(off), vf_new);
        f_sum_new += hsum_avx2(vf_new);
    }

    // Scalar tail for elements 48, 49
    for off in 48..50 {
        let vf = f[off].mul_add(f2f, b_old * d[off]) * e_seg[off];
        f[off] = vf;
        f_sum_new += vf;
    }

    *b = b_old.mul_add(b2b, f_sum_prev * p_repeat_end);
    f_sum_new
}

/// AVX2 backward step: matches C++ backward_step() exactly.
/// See [`forward_step_avx2`] for the upstream FMA dispatch contract.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
pub unsafe fn backward_step_avx2(
    f: &mut [f32; 50],
    d: &[f32; 50],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
) {
    let vf2f = _mm256_set1_ps(f2f);
    let vc = _mm256_set1_ps(p_repeat_end * *b);
    let mut tsum = 0.0f32;

    for off in (0..48).step_by(8) {
        let vf = _mm256_loadu_ps(f.as_ptr().add(off));
        let ve = _mm256_loadu_ps(e_seg.as_ptr().add(off));
        let vd = _mm256_loadu_ps(d.as_ptr().add(off));
        let vf_e = _mm256_mul_ps(vf, ve);
        let vt = _mm256_mul_ps(vf_e, vd);
        tsum += hsum_avx2(vt);
        let vf_new = _mm256_fmadd_ps(vf_e, vf2f, vc);
        _mm256_storeu_ps(f.as_mut_ptr().add(off), vf_new);
    }

    for off in 48..50 {
        let vf = f[off] * e_seg[off];
        tsum += vf * d[off];
        f[off] = vf.mul_add(f2f, p_repeat_end * *b);
    }

    *b = b2b.mul_add(*b, tsum);
}

/// AVX2 scale: multiply all 50 elements by s (matches C++ SIMD::scale)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn scale_avx2(f: &mut [f32; 50], s: f32) {
    let vs = _mm256_set1_ps(s);
    for off in (0..48).step_by(8) {
        let v = _mm256_loadu_ps(f.as_ptr().add(off));
        _mm256_storeu_ps(f.as_mut_ptr().add(off), _mm256_mul_ps(v, vs));
    }
    f[48] *= s;
    f[49] *= s;
}

/// AVX2 sum of all 50 elements (for terminal z computation)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn sum_avx2(f: &[f32; 50]) -> f32 {
    let mut acc = _mm256_setzero_ps();
    for off in (0..48).step_by(8) {
        acc = _mm256_add_ps(acc, _mm256_loadu_ps(f.as_ptr().add(off)));
    }
    hsum_avx2(acc) + f[48] + f[49]
}

/// AArch64 NEON forward step. Upstream's NEON `fmadd` is `vfmaq_f32`, so this
/// intentionally uses fused arithmetic rather than the AVX2 mul/add sequence.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn forward_step_neon(
    f: &mut [f32; 50],
    d: &[f32; 50],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
    f_sum_prev: f32,
) -> f32 {
    use std::arch::aarch64::*;

    let b_old = *b;
    let vf2f = vdupq_n_f32(f2f);
    let vb_old = vdupq_n_f32(b_old);
    let mut f_sum_new = 0.0f32;
    for off in (0..48).step_by(4) {
        let vf = vld1q_f32(f.as_ptr().add(off));
        let vd = vld1q_f32(d.as_ptr().add(off));
        let ve = vld1q_f32(e_seg.as_ptr().add(off));
        let tmp = vfmaq_f32(vmulq_f32(vb_old, vd), vf, vf2f);
        let vf_new = vmulq_f32(tmp, ve);
        vst1q_f32(f.as_mut_ptr().add(off), vf_new);
        f_sum_new += vaddvq_f32(vf_new);
    }
    for off in 48..50 {
        let vf = (f[off] * f2f + b_old * d[off]) * e_seg[off];
        f[off] = vf;
        f_sum_new += vf;
    }
    *b = b_old * b2b + f_sum_prev * p_repeat_end;
    f_sum_new
}

/// AArch64 NEON backward step matching upstream's four-float register path.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn backward_step_neon(
    f: &mut [f32; 50],
    d: &[f32; 50],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
) {
    use std::arch::aarch64::*;

    let vf2f = vdupq_n_f32(f2f);
    let vc = vdupq_n_f32(p_repeat_end * *b);
    let mut tsum = 0.0f32;
    for off in (0..48).step_by(4) {
        let vf = vmulq_f32(
            vld1q_f32(f.as_ptr().add(off)),
            vld1q_f32(e_seg.as_ptr().add(off)),
        );
        let vt = vmulq_f32(vf, vld1q_f32(d.as_ptr().add(off)));
        tsum += vaddvq_f32(vt);
        vst1q_f32(f.as_mut_ptr().add(off), vfmaq_f32(vc, vf, vf2f));
    }
    for off in 48..50 {
        let vf = f[off] * e_seg[off];
        tsum += vf * d[off];
        f[off] = vf * f2f + p_repeat_end * *b;
    }
    *b = b2b * *b + tsum;
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn scale_neon(f: &mut [f32; 50], s: f32) {
    use std::arch::aarch64::*;

    let vs = vdupq_n_f32(s);
    for off in (0..48).step_by(4) {
        let value = vld1q_f32(f.as_ptr().add(off));
        vst1q_f32(f.as_mut_ptr().add(off), vmulq_f32(value, vs));
    }
    f[48] *= s;
    f[49] *= s;
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn sum_neon(f: &[f32; 50]) -> f32 {
    use std::arch::aarch64::*;

    let mut acc = vdupq_n_f32(0.0);
    for off in (0..48).step_by(4) {
        acc = vaddq_f32(acc, vld1q_f32(f.as_ptr().add(off)));
    }
    vaddvq_f32(acc) + f[48] + f[49]
}

#[cfg(all(test, target_arch = "aarch64"))]
mod neon_tests {
    use super::*;

    fn close(left: f32, right: f32) {
        let tolerance = 2.0e-5 * left.abs().max(right.abs()).max(1.0);
        assert!((left - right).abs() <= tolerance, "{left} != {right}");
    }

    #[test]
    fn neon_forward_backward_scale_and_sum_match_reference() {
        let original = std::array::from_fn(|i| (i as f32 + 1.0) / 97.0);
        let d = std::array::from_fn(|i| (50 - i) as f32 / 701.0);
        let e: [f32; 50] = std::array::from_fn(|i| 0.73 + i as f32 / 211.0);
        let f2f = 0.951;
        let p_repeat_end = 0.049;
        let b2b = 0.993;

        let mut forward = original;
        let mut b = 0.625;
        let b_old = b;
        let sum = unsafe {
            forward_step_neon(&mut forward, &d, &e, &mut b, f2f, p_repeat_end, b2b, 1.25)
        };
        let expected_forward: [f32; 50] = std::array::from_fn(|i| {
            let mixed = if i < 48 {
                original[i].mul_add(f2f, b_old * d[i])
            } else {
                original[i] * f2f + b_old * d[i]
            };
            mixed * e[i]
        });
        for i in 0..50 {
            assert_eq!(
                forward[i].to_bits(),
                expected_forward[i].to_bits(),
                "f[{i}]"
            );
        }
        close(sum, expected_forward.iter().sum());
        close(b, b_old * b2b + 1.25 * p_repeat_end);

        let mut backward = original;
        let mut b = 0.625;
        let b_old = b;
        unsafe {
            backward_step_neon(&mut backward, &d, &e, &mut b, f2f, p_repeat_end, b2b);
        }
        let expected_backward: [f32; 50] = std::array::from_fn(|i| {
            let value = original[i] * e[i];
            if i < 48 {
                value.mul_add(f2f, p_repeat_end * b_old)
            } else {
                value * f2f + p_repeat_end * b_old
            }
        });
        for i in 0..50 {
            assert_eq!(
                backward[i].to_bits(),
                expected_backward[i].to_bits(),
                "f[{i}]"
            );
        }
        let tsum: f32 = (0..50).map(|i| original[i] * e[i] * d[i]).sum();
        close(b, b2b * b_old + tsum);

        let mut scaled = original;
        unsafe { scale_neon(&mut scaled, 1.75) };
        for i in 0..50 {
            assert_eq!(scaled[i].to_bits(), (original[i] * 1.75).to_bits());
        }
        close(unsafe { sum_neon(&original) }, original.iter().sum());
    }
}
