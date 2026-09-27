//! Diagonal profile scanning translated from `dp/scan_diags.cpp` and
//! `dp/scan_diags.h`.

use crate::basic::value::{letter_mask, Letter};
use crate::dp::score_profile::LongScoreProfile;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::simd::{arch, Arch};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanArithmetic {
    /// C++ fallback branch: unbounded `int` accumulation.
    Scalar,
    /// C++ SSE4.1/AVX2/NEON branches: signed-saturating `int8_t`, represented
    /// as a local score in the inclusive range 0..=255.
    Saturating8,
}

fn native_arithmetic() -> ScanArithmetic {
    match arch() {
        Arch::Sse4_1 | Arch::Avx2 | Arch::Avx512 | Arch::Neon => ScanArithmetic::Saturating8,
        Arch::None | Arch::Generic => ScanArithmetic::Scalar,
    }
}

pub fn scan_diags128(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    j_begin: i32,
    j_end: i32,
    out: &mut [i32],
) {
    assert!(out.len() >= 128);
    assert!(qp.padding >= 128);
    assert!(j_begin >= 0);
    assert!(j_end >= j_begin);
    let arithmetic = native_arithmetic();
    scan_diags_fixed::<128>(qp, subject, d_begin, j_begin, j_end, out, arithmetic);
}

pub fn scan_diags64(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    j_begin: i32,
    j_end: i32,
    out: &mut [i32],
) {
    assert!(out.len() >= 64);
    assert!(qp.padding >= 128);
    assert!(j_begin >= 0);
    assert!(j_end >= j_begin);
    scan_diags64_with_arithmetic(
        qp,
        subject,
        d_begin,
        j_begin,
        j_end,
        out,
        native_arithmetic(),
    );
}

pub fn scan_diags64_with_arithmetic(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    j_begin: i32,
    j_end: i32,
    out: &mut [i32],
    arithmetic: ScanArithmetic,
) {
    scan_diags_fixed::<64>(qp, subject, d_begin, j_begin, j_end, out, arithmetic);
}

pub fn scan_diags(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    d_end: i32,
    j_begin: i32,
    j_end: i32,
    out: &mut [i32],
) {
    let selected_arch = arch();
    let avx2_branch = matches!(selected_arch, Arch::Avx2 | Arch::Avx512);
    if avx2_branch {
        assert_eq!((d_end - d_begin) % 32, 0);
    }
    let j0 = if avx2_branch {
        j_begin.max(-(d_end - 1))
    } else {
        j_begin.max(-(d_begin + 64 - 1))
    };
    scan_diags_fixed_from::<64>(qp, subject, d_begin, j0, j_end, out, native_arithmetic());
}

fn scan_diags_fixed<const LANES: usize>(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    j_begin: i32,
    j_end: i32,
    out: &mut [i32],
    arithmetic: ScanArithmetic,
) {
    let j0 = j_begin.max(-(d_begin + LANES as i32 - 1));
    scan_diags_fixed_from::<LANES>(qp, subject, d_begin, j0, j_end, out, arithmetic);
}

fn scan_diags_fixed_from<const LANES: usize>(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    j0: i32,
    j_end: i32,
    out: &mut [i32],
    arithmetic: ScanArithmetic,
) {
    assert!(out.len() >= LANES);
    assert!(j0 >= 0, "subject window begins before the sequence");
    let qlen = qp.length() as i32;
    let i0 = d_begin + j0;
    let j1 = (qlen - d_begin).min(j_end);
    assert!(j1 <= subject.len() as i32);

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if arithmetic == ScanArithmetic::Saturating8
        && LANES % 32 == 0
        && std::arch::is_x86_feature_detected!("avx2")
    {
        // SAFETY: AVX2 was detected above. LongScoreProfile guarantees the
        // padding required for each unaligned LANES-byte profile read.
        unsafe {
            scan_diags_avx2::<LANES>(qp, subject, i0, j0, j1, out);
        }
        return;
    }

    #[cfg(target_arch = "aarch64")]
    if arithmetic == ScanArithmetic::Saturating8 && LANES % 16 == 0 {
        // SAFETY: Advanced SIMD is mandatory in AArch64 and the profile owns
        // the padding needed by each unaligned 16-byte load.
        unsafe {
            scan_diags_neon::<LANES>(qp, subject, i0, j0, j1, out);
        }
        return;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if arithmetic == ScanArithmetic::Saturating8
        && LANES % 16 == 0
        && std::arch::is_x86_feature_detected!("sse4.1")
    {
        // SAFETY: SSE4.1 was detected above. LongScoreProfile guarantees the
        // padding required for each unaligned LANES-byte profile read.
        unsafe {
            scan_diags_sse41::<LANES>(qp, subject, i0, j0, j1, out);
        }
        return;
    }

    let mut v = [0i32; LANES];
    let mut max = [0i32; LANES];
    let mut i = i0;
    let mut j = j0;
    while j < j1 {
        let q = profile_get_signed(qp, subject[j as usize], i);
        for k in 0..LANES {
            let next = v[k] + q[k] as i32;
            v[k] = match arithmetic {
                ScanArithmetic::Scalar => next.max(0),
                ScanArithmetic::Saturating8 => next.clamp(0, 255),
            };
            max[k] = max[k].max(v[k]);
        }
        i += 1;
        j += 1;
    }
    out[..LANES].copy_from_slice(&max);
}

/// AArch64 NEON translation of C++'s 16-lane `ScoreVector<int8_t>` path.
/// The `SCHAR_MIN` bias turns signed saturating addition into a local score
/// accumulator in 0..=255, exactly as in the SSE4.1 implementation.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn scan_diags_neon<const LANES: usize>(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    i0: i32,
    j0: i32,
    j1: i32,
    out: &mut [i32],
) {
    use std::arch::aarch64::*;

    const MAX_VECTORS: usize = 8;
    debug_assert!(LANES % 16 == 0 && LANES <= MAX_VECTORS * 16);
    let vector_count = LANES / 16;
    let bias = vdupq_n_s8(i8::MIN);
    let mut score = [bias; MAX_VECTORS];
    let mut best = [bias; MAX_VECTORS];
    let mut i = i0;
    let mut j = j0;
    while j < j1 {
        let profile = profile_get_signed(qp, subject[j as usize], i);
        for block in 0..vector_count {
            let delta = vld1q_s8(profile.as_ptr().add(block * 16));
            score[block] = vqaddq_s8(score[block], delta);
            best[block] = vmaxq_s8(best[block], score[block]);
        }
        i += 1;
        j += 1;
    }

    let mut bytes = [i8::MIN; MAX_VECTORS * 16];
    for block in 0..vector_count {
        vst1q_s8(bytes.as_mut_ptr().add(block * 16), best[block]);
    }
    for (dst, &value) in out[..LANES].iter_mut().zip(&bytes[..LANES]) {
        *dst = i32::from(value) - i32::from(i8::MIN);
    }
}

/// AVX2 translation of C++ `scan_diags{,64,128}`. Values remain biased by
/// `SCHAR_MIN` in registers, making signed saturating addition equivalent to
/// local-alignment accumulation in the range 0..=255.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn scan_diags_avx2<const LANES: usize>(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    i0: i32,
    j0: i32,
    j1: i32,
    out: &mut [i32],
) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    const MAX_VECTORS: usize = 4;
    debug_assert!(LANES % 32 == 0 && LANES <= MAX_VECTORS * 32);
    let vector_count = LANES / 32;
    let bias = _mm256_set1_epi8(i8::MIN);
    let mut score = [bias; MAX_VECTORS];
    let mut best = [bias; MAX_VECTORS];
    let mut i = i0;
    let mut j = j0;
    while j < j1 {
        let profile = profile_get_signed(qp, subject[j as usize], i);
        for block in 0..vector_count {
            let delta = _mm256_loadu_si256(profile.as_ptr().add(block * 32).cast());
            score[block] = _mm256_adds_epi8(score[block], delta);
            best[block] = _mm256_max_epi8(best[block], score[block]);
        }
        i += 1;
        j += 1;
    }

    let mut bytes = [i8::MIN; MAX_VECTORS * 32];
    for block in 0..vector_count {
        _mm256_storeu_si256(bytes.as_mut_ptr().add(block * 32).cast(), best[block]);
    }
    for (dst, &value) in out[..LANES].iter_mut().zip(&bytes[..LANES]) {
        *dst = i32::from(value) - i32::from(i8::MIN);
    }
}

/// SSE4.1 translation of C++ `scan_diags{,64,128}`. As in the AVX2
/// implementation, `SCHAR_MIN` is the zero bias for signed saturating lanes.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
unsafe fn scan_diags_sse41<const LANES: usize>(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    i0: i32,
    j0: i32,
    j1: i32,
    out: &mut [i32],
) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    const MAX_VECTORS: usize = 8;
    debug_assert!(LANES % 16 == 0 && LANES <= MAX_VECTORS * 16);
    let vector_count = LANES / 16;
    let bias = _mm_set1_epi8(i8::MIN);
    let mut score = [bias; MAX_VECTORS];
    let mut best = [bias; MAX_VECTORS];
    let mut i = i0;
    let mut j = j0;
    while j < j1 {
        let profile = profile_get_signed(qp, subject[j as usize], i);
        for block in 0..vector_count {
            let delta = _mm_loadu_si128(profile.as_ptr().add(block * 16).cast());
            score[block] = _mm_adds_epi8(score[block], delta);
            best[block] = _mm_max_epi8(best[block], score[block]);
        }
        i += 1;
        j += 1;
    }

    let mut bytes = [i8::MIN; MAX_VECTORS * 16];
    for block in 0..vector_count {
        _mm_storeu_si128(bytes.as_mut_ptr().add(block * 16).cast(), best[block]);
    }
    for (dst, &value) in out[..LANES].iter_mut().zip(&bytes[..LANES]) {
        *dst = i32::from(value) - i32::from(i8::MIN);
    }
}

fn profile_get_signed(qp: &LongScoreProfile<i8>, letter: Letter, i: i32) -> &[i8] {
    qp.get_signed(letter_mask(letter), i as i64)
}

pub fn diag_alignment(
    scores: &[i32],
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
) -> i32 {
    let mut best = 0;
    let mut best_gap = -gap_open;
    let mut d = -1;
    for (i, &score) in scores.iter().enumerate() {
        if score < gapped_filter_diag_score {
            continue;
        }
        let i = i as i32;
        let gap_score = -gap_extend * (i - d) + best_gap;
        let mut n = score;
        if gap_score + score > best {
            best = gap_score + score;
            n = best;
        }
        if score > best {
            best = score;
            n = score;
        }
        let open_score = -gap_open + n;
        if open_score > gap_score {
            best_gap = open_score;
            d = i;
        }
    }
    best
}

pub fn diag_alignment_with_matrix(
    scores: &[i32],
    score_matrix: &ScoreMatrix,
    gapped_filter_diag_score: i32,
) -> i32 {
    diag_alignment(
        scores,
        score_matrix.gap_open(),
        score_matrix.gap_extend(),
        gapped_filter_diag_score,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;
    use crate::dp::score_profile::make_profile8;

    #[test]
    fn test_scan_diags64_scalar_profile_scores() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5];
        let subject = vec![0, 1, 2, 3];
        let qp = make_profile8(&query, None, 128, &sm);
        let mut out = vec![-1; 64];

        scan_diags64(&qp, &subject, 0, 0, subject.len() as i32, &mut out);

        let mut running = 0;
        let mut best = 0;
        for j in 0..subject.len() {
            running = (running + sm.score(query[j], subject[j])).max(0);
            best = best.max(running);
        }
        assert_eq!(out[0], best);
        assert!(out[1] >= 0);
        assert_eq!(out[63], 0);
    }

    #[test]
    fn test_scan_diags128_and_generic_match_scalar_window() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0; 12];
        let subject = vec![0; 8];
        let qp = make_profile8(&query, None, 128, &sm);
        let mut out128 = vec![0; 128];
        let mut out64 = vec![0; 64];
        let mut out_generic = vec![0; 64];

        scan_diags128(&qp, &subject, -2, 0, subject.len() as i32, &mut out128);
        scan_diags64(&qp, &subject, -2, 0, subject.len() as i32, &mut out64);
        scan_diags(
            &qp,
            &subject,
            -2,
            62,
            0,
            subject.len() as i32,
            &mut out_generic,
        );

        assert_eq!(&out128[..64], &out64[..]);
        assert_eq!(out64, out_generic);
    }

    #[test]
    fn test_diag_alignment_gap_scan() {
        let scores = [3, 0, 5, 0, 7];
        assert_eq!(diag_alignment(&scores, 2, 1, 1), 8);
        assert_eq!(diag_alignment(&scores, 2, 1, 6), 7);
    }

    #[test]
    fn scalar_and_simd_byte_arithmetic_match_cpp_saturation() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![17; 40]; // W/W is sufficiently positive to exceed 255.
        let subject = query.clone();
        let qp = make_profile8(&query, None, 128, &sm);
        let mut scalar = [0; 64];
        let mut simd = [0; 64];
        scan_diags64_with_arithmetic(
            &qp,
            &subject,
            0,
            0,
            subject.len() as i32,
            &mut scalar,
            ScanArithmetic::Scalar,
        );
        scan_diags64_with_arithmetic(
            &qp,
            &subject,
            0,
            0,
            subject.len() as i32,
            &mut simd,
            ScanArithmetic::Saturating8,
        );
        assert!(scalar[0] > 255);
        assert_eq!(simd[0], 255);
    }

    #[test]
    fn subject_letters_follow_cpp_sequence_masking() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3];
        let subject = query.clone();
        let masked: Vec<_> = subject.iter().map(|&letter| letter | SEED_MASK).collect();
        let qp = make_profile8(&query, None, 128, &sm);
        let mut plain = [0; 64];
        let mut masked_out = [0; 64];
        scan_diags64_with_arithmetic(&qp, &subject, 0, 0, 4, &mut plain, ScanArithmetic::Scalar);
        scan_diags64_with_arithmetic(
            &qp,
            &masked,
            0,
            0,
            4,
            &mut masked_out,
            ScanArithmetic::Scalar,
        );
        assert_eq!(masked_out, plain);
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn avx2_scans_match_saturating_reference() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query: Vec<_> = (0..173).map(|i| ((i * 7 + 3) % 25) as Letter).collect();
        let subject: Vec<_> = (0..91).map(|i| ((i * 13 + 2) % 25) as Letter).collect();
        let qp = make_profile8(&query, None, 128, &sm);
        for d_begin in [-40, -7, 0, 23] {
            let j0 = 0.max(-(d_begin + 128 - 1));
            let i0 = d_begin + j0;
            let j1 = (qp.length() as i32 - d_begin).min(subject.len() as i32);
            // Convert the unbounded scalar reference to the exact saturating
            // recurrence rather than merely clamping its final score.
            let mut saturating = [0; 128];
            let mut running = [0; 128];
            for j in j0..j1 {
                let row = profile_get_signed(&qp, subject[j as usize], i0 + j - j0);
                for lane in 0..128 {
                    running[lane] = (running[lane] + i32::from(row[lane])).clamp(0, 255);
                    saturating[lane] = saturating[lane].max(running[lane]);
                }
            }
            let mut actual = [0; 128];
            unsafe {
                scan_diags_avx2::<128>(&qp, &subject, i0, j0, j1, &mut actual);
            }
            assert_eq!(actual, saturating, "d_begin={d_begin}");
        }
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn sse41_scans_match_saturating_reference() {
        if !std::arch::is_x86_feature_detected!("sse4.1") {
            return;
        }
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query: Vec<_> = (0..173).map(|i| ((i * 7 + 3) % 25) as Letter).collect();
        let subject: Vec<_> = (0..91).map(|i| ((i * 13 + 2) % 25) as Letter).collect();
        let qp = make_profile8(&query, None, 128, &sm);
        for d_begin in [-40, -7, 0, 23] {
            let j0 = 0.max(-(d_begin + 128 - 1));
            let i0 = d_begin + j0;
            let j1 = (qp.length() as i32 - d_begin).min(subject.len() as i32);
            let mut saturating = [0; 128];
            let mut running = [0; 128];
            for j in j0..j1 {
                let row = profile_get_signed(&qp, subject[j as usize], i0 + j - j0);
                for lane in 0..128 {
                    running[lane] = (running[lane] + i32::from(row[lane])).clamp(0, 255);
                    saturating[lane] = saturating[lane].max(running[lane]);
                }
            }
            let mut actual = [0; 128];
            unsafe {
                scan_diags_sse41::<128>(&qp, &subject, i0, j0, j1, &mut actual);
            }
            assert_eq!(actual, saturating, "d_begin={d_begin}");
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn neon_scans_match_saturating_reference() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query: Vec<_> = (0..173).map(|i| ((i * 7 + 3) % 25) as Letter).collect();
        let subject: Vec<_> = (0..91).map(|i| ((i * 13 + 2) % 25) as Letter).collect();
        let qp = make_profile8(&query, None, 128, &sm);
        for d_begin in [-40, -7, 0, 23] {
            let j0 = 0.max(-(d_begin + 128 - 1));
            let i0 = d_begin + j0;
            let j1 = (qp.length() as i32 - d_begin).min(subject.len() as i32);
            let mut expected = [0; 128];
            let mut running = [0; 128];
            for j in j0..j1 {
                let row = profile_get_signed(&qp, subject[j as usize], i0 + j - j0);
                for lane in 0..128 {
                    running[lane] = (running[lane] + i32::from(row[lane])).clamp(0, 255);
                    expected[lane] = expected[lane].max(running[lane]);
                }
            }
            let mut actual = [0; 128];
            unsafe {
                scan_diags_neon::<128>(&qp, &subject, i0, j0, j1, &mut actual);
            }
            assert_eq!(actual, expected, "d_begin={d_begin}");
        }
    }

    #[test]
    #[should_panic]
    fn fixed_scans_require_the_cpp_output_width() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let qp = make_profile8(&[0], None, 128, &sm);
        let mut short = [0; 63];
        scan_diags64(&qp, &[0], 0, 0, 1, &mut short);
    }
}
