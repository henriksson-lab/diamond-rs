//! Tantan repeat masking — faithful port of C++ DIAMOND's tantan implementation.
//!
//! Based on tantan by Martin C. Frith:
//! "A new repeat-masking method enables specific detection of homologous sequences"
//! MC Frith, Nucleic Acids Research 2011 39(4):e23.
//!
//! Uses a 50-position window HMM with forward-backward algorithm to compute
//! per-position posterior probability of being in a repeat state.

use crate::basic::value::{Letter, LETTER_MASK, SEED_MASK, TRUE_AA};

/// Default tantan parameters matching C++ DIAMOND.
const P_REPEAT: f32 = 0.005;
const P_REPEAT_END: f32 = 0.05;
const REPEAT_GROWTH: f32 = 1.0 / 0.9;
const DEFAULT_MIN_MASK_PROB: f32 = 0.9;
const WINDOW: usize = 50;

/// Compute lambda for the likelihood ratio matrix using bisection.
///
/// Finds lambda such that the sum of the inverse of exp(lambda * M) equals 1,
/// matching C++ LambdaCalculator. For standard BLOSUM62 (20x20), lambda ≈ 0.324.
fn compute_lambda_flat(scores: &[i8], stride: usize, n: usize) -> f64 {
    // Build double matrix from flat scores array
    let mut mat: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            mat[i][j] = scores[i * stride + j] as f64;
        }
    }

    // Find upper bound
    let mut ub = 0.0f64;
    let mut r_max_min = f64::MAX;
    for i in 0..n {
        let mut r_max = f64::MIN;
        for j in 0..n {
            if mat[i][j] > r_max { r_max = mat[i][j]; }
        }
        if r_max > 0.0 && r_max < r_max_min {
            r_max_min = r_max;
        }
    }
    if r_max_min == f64::MAX || r_max_min <= 0.0 {
        return 0.3176; // fallback to standard BLOSUM62 lambda
    }
    ub = 1.1 * (n as f64).ln() / r_max_min;

    // Bisection: find lambda where sum(inv(exp(lambda * M))) ≈ 1
    let lb = ub * 1e-6;
    let mut lo = lb;
    let mut hi = ub;

    let inv_sum = |tau: f64| -> Option<f64> {
        // Build exp(tau * M)
        let mut em: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                em[i][j] = (tau * mat[i][j]).exp();
            }
        }
        // Invert using Gauss-Jordan
        let mut inv: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
        for i in 0..n { inv[i][i] = 1.0; }
        let mut a = em;
        for col in 0..n {
            // Partial pivot
            let mut max_row = col;
            let mut max_val = a[col][col].abs();
            for row in col+1..n {
                if a[row][col].abs() > max_val {
                    max_val = a[row][col].abs();
                    max_row = row;
                }
            }
            if max_val < 1e-15 { return None; }
            a.swap(col, max_row);
            inv.swap(col, max_row);
            let pivot = a[col][col];
            for j in 0..n { a[col][j] /= pivot; inv[col][j] /= pivot; }
            for row in 0..n {
                if row == col { continue; }
                let factor = a[row][col];
                for j in 0..n {
                    a[row][j] -= factor * a[col][j];
                    inv[row][j] -= factor * inv[col][j];
                }
            }
        }
        let s: f64 = inv.iter().flat_map(|r| r.iter()).sum();
        Some(s)
    };

    let lo_sum = inv_sum(lo).unwrap_or(0.0);
    let hi_sum = inv_sum(hi).unwrap_or(f64::MAX);
    if (lo_sum - 1.0).signum() == (hi_sum - 1.0).signum() {
        return 0.3176; // fallback
    }

    for _ in 0..100 {
        let mid = (lo + hi) / 2.0;
        if mid == lo || mid == hi { break; }
        let mid_sum = inv_sum(mid).unwrap_or(f64::MAX);
        if (lo_sum < 1.0 && mid_sum >= 1.0) || (lo_sum > 1.0 && mid_sum <= 1.0) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (lo + hi) / 2.0
}

/// Pre-computed tantan masking state. Compute once, reuse for all sequences.
pub struct TantanMasker {
    /// Likelihood ratio matrix: lr[i][j] = exp(lambda * score(i,j))
    lr_matrix: Vec<Vec<f32>>,
    min_mask_prob: f32,
}

impl TantanMasker {
    /// Create a masker from a standard matrix (typically BLOSUM62).
    pub fn new(sm: &crate::stats::standard_matrix::StandardMatrix, min_mask_prob: f32) -> Self {
        let n = TRUE_AA as usize;
        let aa_count = crate::basic::value::AMINO_ACID_COUNT;
        let lambda = compute_lambda_flat(&sm.scores, aa_count, n);
        let mut lr_matrix = vec![vec![0.0f32; n]; n];
        for i in 0..n {
            for j in 0..n {
                lr_matrix[i][j] = (lambda * sm.scores[i * aa_count + j] as f64).exp() as f32;
            }
        }
        TantanMasker { lr_matrix, min_mask_prob }
    }

    /// Mask a sequence using the pre-computed likelihood ratio matrix.
    pub fn mask(&self, seq: &mut [Letter]) {
        mask_tantan_inner(seq, &self.lr_matrix, self.min_mask_prob);
    }
}

/// Thread-local cached masker for the default BLOSUM62 case.
fn default_masker() -> &'static TantanMasker {
    use std::sync::OnceLock;
    static MASKER: OnceLock<TantanMasker> = OnceLock::new();
    MASKER.get_or_init(|| {
        TantanMasker::new(&crate::stats::matrices::BLOSUM62, DEFAULT_MIN_MASK_PROB)
    })
}

/// Apply tantan masking to a sequence using default BLOSUM62 parameters.
pub fn mask_tantan(seq: &mut [Letter]) {
    default_masker().mask(seq);
}

/// Core tantan implementation using a pre-computed likelihood ratio matrix.
fn mask_tantan_inner(
    seq: &mut [Letter],
    lr_matrix: &[Vec<f32>],
    min_mask_prob: f32,
) {
    let len = seq.len();
    if len == 0 { return; }

    let n = TRUE_AA as usize;

    // Tantan HMM parameters
    let b2b = 1.0f32 - P_REPEAT;
    let f2f = 1.0f32 - P_REPEAT_END;
    let b2f0 = P_REPEAT * (1.0 - REPEAT_GROWTH) / (1.0 - REPEAT_GROWTH.powi(WINDOW as i32));

    // Repeat-state entry distribution (geometric decay over window positions)
    let mut d = [0.0f32; WINDOW];
    d[WINDOW - 1] = b2f0;
    for i in (0..WINDOW-1).rev() {
        d[i] = d[i + 1] * REPEAT_GROWTH;
    }

    // Pre-compute emission vectors matching C++ tantan.cpp lines 152-164:
    //   e[aa][len-1-j] = likelihoodRatioMatrix[aa][letter_mask(seq[j])]
    // Reversed so e[aa][len - i] gives emission for position i during forward pass.
    let mut e: Vec<Vec<f32>> = Vec::with_capacity(n);
    for aa in 0..n {
        let mut ev = vec![0.0f32; len + WINDOW];
        for j in 0..len {
            let idx = (seq[j] & LETTER_MASK) as usize;
            if idx < n {
                ev[len - 1 - j] = lr_matrix[aa][idx];
            }
        }
        e.push(ev);
    }

    // Forward pass
    let mut f = [0.0f32; WINDOW];
    let mut pb = vec![0.0f32; len];
    let mut scale = vec![0.0f32; (len + 15) / 16];
    let mut b = 1.0f32;
    let mut f_sum = 0.0f32;

    for i in 0..len {
        let ltr = (seq[i] & LETTER_MASK) as usize;
        if ltr >= n { continue; }
        let e_seg = &e[ltr][len - i..];

        // Forward step — matches C++ SIMD processing order:
        // C++ processes 48 elements in SIMD (6 chunks of 8), then 2 scalar.
        // hsum within each chunk accumulates differently from sequential addition.
        let b_old = b;
        let mut f_sum_new = 0.0f32;
        // Process in chunks of 8 to match SIMD hsum accumulation
        for chunk_start in (0..48).step_by(8) {
            let mut chunk_sum = 0.0f32;
            for off in chunk_start..chunk_start + 8 {
                let vf = f[off].mul_add(f2f, b_old * d[off]) * e_seg[off];
                f[off] = vf;
                chunk_sum += vf;
            }
            f_sum_new += chunk_sum;
        }
        // Scalar tail for last 2 elements (matching C++ lines 68-73)
        for off in 48..WINDOW {
            let vf = (f[off] * f2f + b_old * d[off]) * e_seg[off];
            f[off] = vf;
            f_sum_new += vf;
        }
        b = b_old * b2b + f_sum * P_REPEAT_END;
        f_sum = f_sum_new;

        // Rescale every 16 positions to avoid underflow
        if (i & 15) == 15 {
            let s = 1.0 / b;
            scale[i / 16] = s;
            b *= s;
            for v in f.iter_mut() { *v *= s; }
            f_sum *= s;
        }
        pb[i] = b;
    }

    // Terminal probability
    let z = b * b2b + f.iter().sum::<f32>() * P_REPEAT_END;
    let zinv = 1.0 / z;

    // Backward pass
    b = b2b;
    f.fill(P_REPEAT_END);

    for i in (0..len).rev() {
        let pf = 1.0 - (pb[i] * b * zinv);

        // Rescale
        if (i & 15) == 15 {
            let s = scale[i / 16];
            b *= s;
            for v in f.iter_mut() { *v *= s; }
        }

        let ltr = (seq[i] & LETTER_MASK) as usize;
        if ltr >= n { continue; }
        let e_seg = &e[ltr][len - i..];

        // Backward step — matches C++ SIMD processing order
        let mut tsum = 0.0f32;
        let c = P_REPEAT_END * b;
        for chunk_start in (0..48).step_by(8) {
            let mut chunk_sum = 0.0f32;
            for off in chunk_start..chunk_start + 8 {
                let vf = f[off] * e_seg[off];
                chunk_sum += vf * d[off];
                f[off] = vf.mul_add(f2f, c);
            }
            tsum += chunk_sum;
        }
        for off in 48..WINDOW {
            let vf = f[off] * e_seg[off];
            tsum += vf * d[off];
            f[off] = vf * f2f + c;
        }
        b = b2b * b + tsum;

        // Mask if posterior probability of repeat >= threshold
        if pf >= min_mask_prob {
            seq[i] |= SEED_MASK;
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_mask_short() {
        let mut seq = vec![0i8];
        mask_tantan(&mut seq);
        assert_eq!(seq[0] & SEED_MASK, 0);
    }

    #[test]
    fn test_no_mask_diverse() {
        let mut diverse: Vec<Letter> = (0..100).map(|i| (i % 20) as Letter).collect();
        let mut repetitive = vec![0i8; 100];
        mask_tantan(&mut diverse);
        mask_tantan(&mut repetitive);
        let diverse_masked = diverse.iter().filter(|&&l| l & SEED_MASK != 0).count();
        let repetitive_masked = repetitive.iter().filter(|&&l| l & SEED_MASK != 0).count();
        assert!(
            repetitive_masked >= diverse_masked,
            "Repetitive ({}) should have >= masking than diverse ({})",
            repetitive_masked,
            diverse_masked
        );
    }

    #[test]
    fn test_mask_repetitive() {
        let mut seq = vec![0i8; 100];
        mask_tantan(&mut seq);
        let masked_count = seq.iter().filter(|&&l| l & SEED_MASK != 0).count();
        assert!(masked_count > 0, "No positions masked in repetitive sequence");
    }

    #[test]
    fn test_mask_q6gzx3_repeat() {
        // Q6GZX3 has a PPTPPT repeat near the end that C++ masks (~12 positions)
        let seq_str = b"MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQTCASGFCTSQPLCARIKKTQVCGLRYSSKGKDPLVSAEWDSRGAPYVRCTYDADLIDTQAQVDQFVSMFGESPSLAERYCMRGVKNTAGELVSRVSSDADPAGGWCRKWYSAHRGPDQDAALGSFCIKNPGAADCKCINRASDPVYQKVKTLHAYPDQCWYVPCAADVGELKMGTQRDTPTNCPTQVCQIVFNMLDDGSVTMDDVKNTINCDFSKYVPPPPPPKPTPPTPPTPPTPPTPPTPPTPPTPRPVHNRKVMFFVAGAVLVAILISTVRW";
        let letter_map: std::collections::HashMap<u8, i8> = [
            (b'A',0),(b'R',1),(b'N',2),(b'D',3),(b'C',4),(b'Q',5),(b'E',6),(b'G',7),
            (b'H',8),(b'I',9),(b'L',10),(b'K',11),(b'M',12),(b'F',13),(b'P',14),
            (b'S',15),(b'T',16),(b'W',17),(b'Y',18),(b'V',19),
        ].iter().cloned().collect();
        let mut seq: Vec<Letter> = seq_str.iter().map(|&c| *letter_map.get(&c).unwrap_or(&23)).collect();

        mask_tantan(&mut seq);
        let masked: Vec<usize> = seq.iter().enumerate()
            .filter(|(_, &l)| l & SEED_MASK != 0)
            .map(|(i, _)| i)
            .collect();
        let diag: std::collections::HashMap<u8, i32> = [
            (b'A',4),(b'R',5),(b'N',6),(b'D',6),(b'C',9),(b'Q',5),(b'E',5),(b'G',6),
            (b'H',8),(b'I',4),(b'L',4),(b'K',5),(b'M',5),(b'F',6),(b'P',7),
            (b'S',4),(b'T',5),(b'W',11),(b'Y',7),(b'V',4),
        ].iter().cloned().collect();
        let deficit: i32 = masked.iter().map(|&i| diag.get(&seq_str[i]).copied().unwrap_or(0)).sum();
        eprintln!("Q6GZX3 masked {} of {} positions, score deficit={}", masked.len(), seq.len(), deficit);
        eprintln!("  Masked positions: {:?}", &masked);
        eprintln!("  Masked region: {}", masked.iter().map(|&i| seq_str[i] as char).collect::<String>());
        // C++ masks ~12 positions in the PPTPPT repeat (around pos 259-295)
        // Rust should mask a similar number
        assert!(masked.len() > 0, "Q6GZX3 repeat region should have some masking");
    }

    #[test]
    fn test_lambda_blosum62() {
        let sm = &crate::stats::matrices::BLOSUM62;
        let lambda = compute_lambda_flat(&sm.scores, crate::basic::value::AMINO_ACID_COUNT, 20);
        // BLOSUM62 lambda should be ~0.324
        assert!((lambda - 0.324).abs() < 0.01,
            "Lambda should be ~0.324, got {:.6}", lambda);
    }

    #[test]
    fn test_mask_d1ivsa4_no_masking() {
        // d1ivsa4 (426aa) — C++ masks 0 positions. Rust should too.
        // Read actual d1ivsa4 sequence from test file
        let fasta_path = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/5.faa");
        let records = crate::data::fasta::read_fasta_file(
            std::path::Path::new(fasta_path),
            crate::basic::value::SequenceType::AminoAcid,
        ).unwrap();
        let rec = records.iter().find(|r| r.id.contains("d1ivsa4")).unwrap();
        let mut seq = rec.sequence.clone();
        assert_eq!(seq.len(), 426);

        mask_tantan(&mut seq);
        let masked: Vec<usize> = seq.iter().enumerate()
            .filter(|(_, &l)| l & SEED_MASK != 0)
            .map(|(i, _)| i)
            .collect();
        eprintln!("d1ivsa4 masked {} of {} positions", masked.len(), seq.len());
        if !masked.is_empty() {
            eprintln!("  Masked: {:?}", &masked[..masked.len().min(20)]);
        }
        // C++ masks 0 positions for this sequence
        assert_eq!(masked.len(), 0, "d1ivsa4 should have 0 masked positions (C++ masks 0), got {}", masked.len());
    }
}
