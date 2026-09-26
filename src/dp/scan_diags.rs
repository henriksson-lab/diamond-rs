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

    #[test]
    #[should_panic]
    fn fixed_scans_require_the_cpp_output_width() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let qp = make_profile8(&[0], None, 128, &sm);
        let mut short = [0; 63];
        scan_diags64(&qp, &[0], 0, 0, 1, &mut short);
    }
}
