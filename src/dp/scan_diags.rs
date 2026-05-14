use crate::basic::value::Letter;
use crate::dp::score_profile::LongScoreProfile;
use crate::stats::score_matrix::ScoreMatrix;

pub fn scan_diags128(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    j_begin: i32,
    j_end: i32,
    out: &mut [i32],
) {
    scan_diags_fixed::<128>(qp, subject, d_begin, j_begin, j_end, out);
}

pub fn scan_diags64(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    j_begin: i32,
    j_end: i32,
    out: &mut [i32],
) {
    scan_diags_fixed::<64>(qp, subject, d_begin, j_begin, j_end, out);
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
    let _band = d_end - d_begin;
    scan_diags_fixed::<64>(qp, subject, d_begin, j_begin, j_end, out);
}

fn scan_diags_fixed<const LANES: usize>(
    qp: &LongScoreProfile<i8>,
    subject: &[Letter],
    d_begin: i32,
    j_begin: i32,
    j_end: i32,
    out: &mut [i32],
) {
    let qlen = qp.length() as i32;
    let j0 = j_begin.max(-(d_begin + LANES as i32 - 1));
    let i0 = d_begin + j0;
    let j1 = (qlen - d_begin).min(j_end).min(subject.len() as i32);
    let lanes = LANES.min(out.len());
    let mut v = vec![0i32; lanes];
    let mut max = vec![0i32; lanes];
    let mut i = i0;
    let mut j = j0;
    while j < j1 {
        if j >= 0 {
            let q = profile_get_signed(qp, subject[j as usize], i);
            for k in 0..lanes {
                v[k] += q[k] as i32;
                v[k] = v[k].max(0);
                max[k] = max[k].max(v[k]);
            }
        }
        i += 1;
        j += 1;
    }
    out[..lanes].copy_from_slice(&max);
}

fn profile_get_signed(qp: &LongScoreProfile<i8>, letter: Letter, i: i32) -> &[i8] {
    let row = &qp.data[letter as usize];
    let pos = (i + qp.padding as i32) as usize;
    &row[pos..]
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
}
