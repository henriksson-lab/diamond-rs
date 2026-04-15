//! AVX2 SIMD tantan forward/backward steps.
//!
//! Direct port of C++ `masking/tantan.cpp` SIMD path using Rust `std::arch::x86_64`.
//! Processes 8 floats at a time matching the C++ AVX2 dispatch.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// Check if AVX2 + FMA are available at runtime.
#[cfg(target_arch = "x86_64")]
pub fn has_avx2_fma() -> bool {
    is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma")
}

#[cfg(not(target_arch = "x86_64"))]
pub fn has_avx2_fma() -> bool {
    false
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
/// Uses FMA: fmadd(f, f2f, b_old * d) * e_seg
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2", enable = "fma")]
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
        // tmp = fmadd(vf, vf2f, vb_old * vd) = vf * f2f + b_old * d
        let tmp = _mm256_fmadd_ps(vf, vf2f, _mm256_mul_ps(vb_old, vd));
        // vf = tmp * e_seg
        let vf_new = _mm256_mul_ps(tmp, ve);
        _mm256_storeu_ps(f.as_mut_ptr().add(off), vf_new);
        f_sum_new += hsum_avx2(vf_new);
    }

    // Scalar tail for elements 48, 49
    for off in 48..50 {
        let vf = (f[off] * f2f + b_old * d[off]) * e_seg[off];
        f[off] = vf;
        f_sum_new += vf;
    }

    *b = b_old * b2b + f_sum_prev * p_repeat_end;
    f_sum_new
}

/// AVX2 backward step: matches C++ backward_step() exactly.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2", enable = "fma")]
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
        // vf = vf * e_seg
        let vf_e = _mm256_mul_ps(vf, ve);
        // tsum += vf * d
        let vt = _mm256_mul_ps(vf_e, vd);
        tsum += hsum_avx2(vt);
        // f[off] = fmadd(vf_e, f2f, c)
        let vf_new = _mm256_fmadd_ps(vf_e, vf2f, vc);
        _mm256_storeu_ps(f.as_mut_ptr().add(off), vf_new);
    }

    for off in 48..50 {
        let vf = f[off] * e_seg[off];
        tsum += vf * d[off];
        f[off] = vf * f2f + p_repeat_end * *b;
    }

    *b = b2b * *b + tsum;
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
