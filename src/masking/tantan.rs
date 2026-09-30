//! Tantan repeat masking — faithful port of C++ DIAMOND's tantan implementation.
//!
//! Based on tantan by Martin C. Frith:
//! "A new repeat-masking method enables specific detection of homologous sequences"
//! MC Frith, Nucleic Acids Research 2011 39(4):e23.
//!
//! Uses a 50-position window HMM with forward-backward algorithm to compute
//! per-position posterior probability of being in a repeat state.

#[cfg(test)]
use crate::basic::value::DELIMITER_LETTER;
use crate::basic::value::{Letter, AMINO_ACID_COUNT, LETTER_MASK, MASK_LETTER, SEED_MASK, TRUE_AA};
use crate::masking::Ranges;
use crate::stats::score_matrix::ScoreMatrix;
#[cfg(all(test, target_arch = "x86_64"))]
use std::cell::Cell;
use std::cell::RefCell;

/// Default tantan parameters matching C++ DIAMOND.
const P_REPEAT: f32 = 0.005;
const P_REPEAT_END: f32 = 0.05;
const REPEAT_GROWTH: f32 = 1.0 / 0.9;
const DEFAULT_MIN_MASK_PROB: f32 = 0.9;
const WINDOW: usize = 50;
const EMISSION_ROWS: usize = AMINO_ACID_COUNT;

#[derive(Default)]
struct TantanScratch {
    emission: [Vec<f32>; EMISSION_ROWS],
    forward_background: Vec<f32>,
    scale: Vec<f32>,
}

thread_local! {
    /// Match the upstream implementation's thread-local work buffers. Tantan
    /// runs once per sequence, so allocating 32 emission matrices for every
    /// record otherwise creates substantial allocator traffic under Rayon.
    static TANTAN_SCRATCH: RefCell<TantanScratch> = RefCell::new(TantanScratch::default());
    #[cfg(all(test, target_arch = "x86_64"))]
    static FORCE_SCALAR: Cell<bool> = const { Cell::new(false) };
}

#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy, PartialEq, Eq)]
enum X86Backend {
    Avx2,
    Sse,
    Scalar,
}

#[cfg(target_arch = "x86_64")]
#[inline]
fn x86_backend() -> X86Backend {
    #[cfg(test)]
    if FORCE_SCALAR.with(Cell::get) {
        return X86Backend::Scalar;
    }
    if super::tantan_simd::has_avx2() {
        X86Backend::Avx2
    } else if super::tantan_simd::has_sse41_ssse3() {
        X86Backend::Sse
    } else {
        X86Backend::Scalar
    }
}

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
    let mut r_max_min = f64::MAX;
    for i in 0..n {
        let mut r_max = f64::MIN;
        for j in 0..n {
            if mat[i][j] > r_max {
                r_max = mat[i][j];
            }
        }
        if r_max > 0.0 && r_max < r_max_min {
            r_max_min = r_max;
        }
    }
    if r_max_min == f64::MAX || r_max_min <= 0.0 {
        return 0.3176; // fallback to standard BLOSUM62 lambda
    }
    let ub = 1.1 * (n as f64).ln() / r_max_min;

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
        for i in 0..n {
            inv[i][i] = 1.0;
        }
        let mut a = em;
        for col in 0..n {
            // Partial pivot
            let mut max_row = col;
            let mut max_val = a[col][col].abs();
            for row in col + 1..n {
                if a[row][col].abs() > max_val {
                    max_val = a[row][col].abs();
                    max_row = row;
                }
            }
            if max_val < 1e-15 {
                return None;
            }
            a.swap(col, max_row);
            inv.swap(col, max_row);
            let pivot = a[col][col];
            for j in 0..n {
                a[col][j] /= pivot;
                inv[col][j] /= pivot;
            }
            for row in 0..n {
                if row == col {
                    continue;
                }
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
        if mid == lo || mid == hi {
            break;
        }
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
        // Lambda is computed on the 20x20 standard-AA submatrix to match C++'s
        // LambdaCalculator (masking.cpp:Masking ctor uses int_matrix[20][20]).
        let n = TRUE_AA as usize;
        let aa_count = AMINO_ACID_COUNT;
        let lambda = compute_lambda_flat(&sm.scores, aa_count, n);
        // Likelihood ratio matrix covers the full 26-letter alphabet (incl.
        // B, J, Z, X, *, _). C++ Masking::Masking populates
        // `likelihoodRatioMatrixf_[i][j]` for `i,j < value_traits.alphabet_size`
        // (= 26 for amino acid). If we only fill the 20x20 sub-block, X-letter
        // positions emit 0 while C++ emits exp(lambda*score(aa,X)) ≈ 0.7 —
        // flipping mask decisions at the p_mask boundary.
        let mut lr_matrix = vec![vec![0.0f32; aa_count]; aa_count];
        for i in 0..aa_count {
            for j in 0..aa_count {
                lr_matrix[i][j] = (lambda * sm.scores[i * aa_count + j] as f64).exp() as f32;
            }
        }
        TantanMasker {
            lr_matrix,
            min_mask_prob,
        }
    }

    /// Create a masker from the active scoring matrix.
    ///
    /// Matches C++ `Masking::Masking(const ScoreMatrix&)`, which constructs
    /// tantan likelihood ratios from the current matrix rather than always
    /// using BLOSUM62.
    pub fn from_score_matrix(score_matrix: &ScoreMatrix, min_mask_prob: f32) -> Self {
        let n = TRUE_AA as usize;
        let aa_count = AMINO_ACID_COUNT;
        let mut scores = vec![0i8; aa_count * aa_count];
        for i in 0..aa_count {
            for j in 0..aa_count {
                scores[i * aa_count + j] = score_matrix.score(i as Letter, j as Letter) as i8;
            }
        }
        let lambda = compute_lambda_flat(&scores, aa_count, n);
        let mut lr_matrix = vec![vec![0.0f32; aa_count]; aa_count];
        for i in 0..aa_count {
            for j in 0..aa_count {
                let score = score_matrix.score(i as Letter, j as Letter);
                lr_matrix[i][j] = (lambda * score as f64).exp() as f32;
            }
        }
        TantanMasker {
            lr_matrix,
            min_mask_prob,
        }
    }

    /// Mask a sequence using the pre-computed likelihood ratio matrix.
    pub fn mask(&self, seq: &mut [Letter]) {
        self.mask_bit(seq);
    }

    /// C++ `Masking::operator()` path without a masking table.
    pub fn mask_hard(&self, seq: &mut [Letter]) -> Ranges {
        mask(
            seq,
            &self.lr_matrix,
            P_REPEAT,
            P_REPEAT_END,
            REPEAT_GROWTH,
            self.min_mask_prob,
            1,
        )
    }

    /// C++ `Masking::operator()` path with a masking table.
    pub fn mask_ranges(&self, seq: &mut [Letter]) -> Ranges {
        mask(
            seq,
            &self.lr_matrix,
            P_REPEAT,
            P_REPEAT_END,
            REPEAT_GROWTH,
            self.min_mask_prob,
            0,
        )
    }

    /// C++ `Masking::mask_bit` path.
    pub fn mask_bit(&self, seq: &mut [Letter]) {
        mask(
            seq,
            &self.lr_matrix,
            P_REPEAT,
            P_REPEAT_END,
            REPEAT_GROWTH,
            self.min_mask_prob,
            2,
        );
    }
}

/// Thread-local cached masker for the default BLOSUM62 case.
fn default_masker() -> &'static TantanMasker {
    use std::sync::OnceLock;
    static MASKER: OnceLock<TantanMasker> = OnceLock::new();
    MASKER
        .get_or_init(|| TantanMasker::new(&crate::stats::matrices::BLOSUM62, DEFAULT_MIN_MASK_PROB))
}

/// Apply tantan masking to a sequence using default BLOSUM62 parameters.
pub fn mask_tantan(seq: &mut [Letter]) {
    default_masker().mask(seq);
}

/// Apply tantan masking using the active score matrix.
pub fn mask_tantan_with_score_matrix(seq: &mut [Letter], score_matrix: &ScoreMatrix) {
    TantanMasker::from_score_matrix(score_matrix, DEFAULT_MIN_MASK_PROB).mask(seq);
}

/// Scalar forward step (generic fallback).
#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
fn forward_step_scalar(
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
    let mut f_sum_new = 0.0f32;
    for off in 0..50 {
        let vf = (f[off] * f2f + b_old * d[off]) * e_seg[off];
        f[off] = vf;
        f_sum_new += vf;
    }
    *b = b_old * b2b + f_sum_prev * p_repeat_end;
    f_sum_new
}

/// Scalar backward step (generic fallback).
#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
fn backward_step_scalar(
    f: &mut [f32; 50],
    d: &[f32; 50],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
) -> f32 {
    let mut tsum = 0.0f32;
    let c = p_repeat_end * *b;
    for off in 0..50 {
        let vf = f[off] * e_seg[off];
        tsum += vf * d[off];
        f[off] = vf * f2f + c;
    }
    *b = b2b * *b + tsum;
    tsum
}

/// Advance the forward probabilities by one residue.
///
/// This is the snake-case Rust mapping of C++ `forward_step`. Runtime dispatch
/// replaces C++'s compile-time `DISPATCH_ARCH` instantiations.
fn forward_step(
    f: &mut [f32; WINDOW],
    d: &[f32; WINDOW],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
    f_sum_prev: f32,
) -> f32 {
    #[cfg(target_arch = "x86_64")]
    match x86_backend() {
        // SAFETY: runtime dispatch established the required target features.
        X86Backend::Avx2 => unsafe {
            return super::tantan_simd::forward_step_avx2(
                f,
                d,
                e_seg,
                b,
                f2f,
                p_repeat_end,
                b2b,
                f_sum_prev,
            );
        },
        X86Backend::Sse => unsafe {
            return super::tantan_simd::forward_step_sse(
                f,
                d,
                e_seg,
                b,
                f2f,
                p_repeat_end,
                b2b,
                f_sum_prev,
            );
        },
        X86Backend::Scalar => {}
    }

    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: Advanced SIMD is mandatory on AArch64.
        unsafe {
            super::tantan_simd::forward_step_neon(
                f,
                d,
                e_seg,
                b,
                f2f,
                p_repeat_end,
                b2b,
                f_sum_prev,
            )
        }
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        forward_step_scalar(f, d, e_seg, b, f2f, p_repeat_end, b2b, f_sum_prev)
    }
}

/// Advance the backward probabilities by one residue and return `tsum`.
///
/// This restores the direct Rust symbol mapping for C++ `backward_step`; the
/// return value is retained even though `mask` only needs its update of `b`.
fn backward_step(
    f: &mut [f32; WINDOW],
    d: &[f32; WINDOW],
    e_seg: &[f32],
    b: &mut f32,
    f2f: f32,
    p_repeat_end: f32,
    b2b: f32,
) -> f32 {
    #[cfg(target_arch = "x86_64")]
    match x86_backend() {
        // The architecture helper performs the same update but does not expose
        // C++'s otherwise-unused `tsum`. Recover it from the defining equation.
        X86Backend::Avx2 => unsafe {
            let b_old = *b;
            super::tantan_simd::backward_step_avx2(f, d, e_seg, b, f2f, p_repeat_end, b2b);
            return *b - b2b * b_old;
        },
        X86Backend::Sse => unsafe {
            let b_old = *b;
            super::tantan_simd::backward_step_sse(f, d, e_seg, b, f2f, p_repeat_end, b2b);
            return *b - b2b * b_old;
        },
        X86Backend::Scalar => {}
    }

    #[cfg(target_arch = "aarch64")]
    {
        // The architecture helper updates `b`; recover C++'s otherwise-unused
        // `tsum` from the defining equation, as in the AVX2 dispatch.
        let b_old = *b;
        unsafe {
            super::tantan_simd::backward_step_neon(f, d, e_seg, b, f2f, p_repeat_end, b2b);
        }
        *b - b2b * b_old
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        backward_step_scalar(f, d, e_seg, b, f2f, p_repeat_end, b2b)
    }
}

/// Run the tantan forward-backward masker.
///
/// This is the direct translation of C++ `Util::tantan::mask`. `mask_mode` has
/// the same values as the original: 0 records ranges without changing the
/// sequence, 1 replaces residues with `MASK_LETTER`, and 2 sets `SEED_MASK`.
pub fn mask(
    seq: &mut [Letter],
    likelihood_ratio_matrix: &[Vec<f32>],
    p_repeat: f32,
    p_repeat_end: f32,
    repeat_growth: f32,
    p_mask: f32,
    mask_mode: i32,
) -> Ranges {
    let len = seq.len();
    if len == 0 {
        return Ranges::new();
    }

    TANTAN_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        let TantanScratch {
            emission,
            forward_background: pb,
            scale,
        } = &mut *scratch;
        let mut ranges = Ranges::new();

        let alphabet_size = AMINO_ACID_COUNT;

        // Tantan HMM parameters
        let b2b = 1.0f32 - p_repeat;
        let f2f = 1.0f32 - p_repeat_end;
        // C++ calls the floating-point overload `std::pow(float, float)` here.
        let b2f0 = p_repeat * (1.0 - repeat_growth) / (1.0 - repeat_growth.powf(WINDOW as f32));

        // Repeat-state entry distribution (geometric decay over window positions)
        let mut d = [0.0f32; WINDOW];
        d[WINDOW - 1] = b2f0;
        for i in (0..WINDOW - 1).rev() {
            d[i] = d[i + 1] * repeat_growth;
        }

        // Pre-compute emission vectors matching C++ tantan.cpp:152-164. C++ fills an
        // `e[aa]` row for every aa in 0..AMINO_ACID_COUNT (=26), and writes
        // `L[letter_mask(seq[j])]` without bounds-checking idx — meaning ambiguous
        // letters (B/J/Z/X/*/_, ids 20..25) DO contribute non-zero emissions via the
        // populated 26x26 lr_matrix block. Only DELIMITER (31) and similarly-stripped
        // values read past the 26-column initialized region (C++ UB, typically zero).
        // Normal protein records contain only the 26-letter alphabet after
        // stripping the mask bit. Validate that invariant once, then mirror
        // upstream's raw-pointer transpose without two bounds checks in each
        // of the 26 * len cells. Keep a safe zero-filling fallback for the
        // public low-level API's malformed/custom inputs.
        let all_letters_valid = seq
            .iter()
            .all(|&letter| ((letter & LETTER_MASK) as usize) < alphabet_size);
        let direct_emissions = likelihood_ratio_matrix.len() >= alphabet_size
            && likelihood_ratio_matrix[..alphabet_size]
                .iter()
                .all(|row| row.len() >= alphabet_size)
            && all_letters_valid;
        for aa in 0..alphabet_size {
            let ev = &mut emission[aa];
            ev.resize(len + WINDOW, 0.0);
            if direct_emissions {
                // SAFETY: the one-time checks above prove every source row and
                // stripped sequence letter is in bounds; resize established
                // len + WINDOW initialized destination elements.
                unsafe {
                    let source = likelihood_ratio_matrix.get_unchecked(aa).as_ptr();
                    let destination = ev.as_mut_ptr();
                    for j in 0..len {
                        let idx = (*seq.get_unchecked(j) & LETTER_MASK) as usize;
                        *destination.add(len - 1 - j) = *source.add(idx);
                    }
                    std::ptr::write_bytes(destination.add(len), 0, WINDOW);
                }
            } else {
                for j in 0..len {
                    let idx = (seq[j] & LETTER_MASK) as usize;
                    ev[len - 1 - j] = likelihood_ratio_matrix
                        .get(aa)
                        .and_then(|row| row.get(idx))
                        .copied()
                        .unwrap_or(0.0);
                }
                ev[len..len + WINDOW].fill(0.0);
            }
        }
        // The public low-level entry point accepts arbitrary encoded bytes.
        // Upstream reads outside its 26-row pointer table for those values;
        // define that malformed-input case as a zero-emission row instead.
        // Normal protein input never allocates this fallback.
        let invalid_emission = if all_letters_valid {
            Vec::new()
        } else {
            vec![0.0; len + WINDOW]
        };

        let mut f = [0.0f32; WINDOW];
        let mut d_arr = [0.0f32; WINDOW];
        d_arr.copy_from_slice(&d);
        pb.resize(len, 0.0);
        scale.resize((len + 15) / 16, 0.0);
        let mut b = 1.0f32;
        let mut f_sum = 0.0f32;

        // Forward pass
        for i in 0..len {
            let ltr = (seq[i] & LETTER_MASK) as usize;
            let e_row = if ltr < alphabet_size {
                // SAFETY: the comparison proves the row is in bounds.
                unsafe { emission.get_unchecked(ltr) }
            } else {
                &invalid_emission
            };
            let e_seg = &e_row[len - i..];

            f_sum = forward_step(&mut f, &d_arr, e_seg, &mut b, f2f, p_repeat_end, b2b, f_sum);

            // Rescale every 16 positions to avoid underflow
            if (i & 15) == 15 {
                let s = 1.0 / b;
                scale[i / 16] = s;
                b *= s;
                #[cfg(target_arch = "x86_64")]
                match x86_backend() {
                    X86Backend::Avx2 => unsafe { super::tantan_simd::scale_avx2(&mut f, s) },
                    X86Backend::Sse => unsafe { super::tantan_simd::scale_sse(&mut f, s) },
                    X86Backend::Scalar => f.iter_mut().for_each(|v| *v *= s),
                }
                #[cfg(target_arch = "aarch64")]
                unsafe {
                    super::tantan_simd::scale_neon(&mut f, s);
                }
                #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
                for v in f.iter_mut() {
                    *v *= s;
                }
                f_sum *= s;
            }
            pb[i] = b;
        }

        // Terminal probability
        let f_total = {
            #[cfg(target_arch = "x86_64")]
            match x86_backend() {
                X86Backend::Avx2 => unsafe { super::tantan_simd::sum_avx2(&f) },
                X86Backend::Sse => unsafe { super::tantan_simd::sum_sse(&f) },
                X86Backend::Scalar => f.iter().sum::<f32>(),
            }
            #[cfg(not(target_arch = "x86_64"))]
            #[cfg(not(target_arch = "aarch64"))]
            {
                f.iter().sum::<f32>()
            }
            #[cfg(target_arch = "aarch64")]
            unsafe {
                super::tantan_simd::sum_neon(&f)
            }
        };
        let z = b * b2b + f_total * p_repeat_end;
        let zinv = 1.0 / z;

        // Backward pass
        b = b2b;
        f.fill(p_repeat_end);

        for i in (0..len).rev() {
            let pf = 1.0 - (pb[i] * b * zinv);

            // Rescale
            if (i & 15) == 15 {
                let s = scale[i / 16];
                b *= s;
                #[cfg(target_arch = "x86_64")]
                match x86_backend() {
                    X86Backend::Avx2 => unsafe { super::tantan_simd::scale_avx2(&mut f, s) },
                    X86Backend::Sse => unsafe { super::tantan_simd::scale_sse(&mut f, s) },
                    X86Backend::Scalar => f.iter_mut().for_each(|v| *v *= s),
                }
                #[cfg(target_arch = "aarch64")]
                unsafe {
                    super::tantan_simd::scale_neon(&mut f, s);
                }
                #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
                for v in f.iter_mut() {
                    *v *= s;
                }
            }

            let ltr = (seq[i] & LETTER_MASK) as usize;
            let e_row = if ltr < alphabet_size {
                // SAFETY: the comparison proves the row is in bounds.
                unsafe { emission.get_unchecked(ltr) }
            } else {
                &invalid_emission
            };
            let e_seg = &e_row[len - i..];

            backward_step(&mut f, &d_arr, e_seg, &mut b, f2f, p_repeat_end, b2b);

            if pf >= p_mask {
                if mask_mode == 1 {
                    seq[i] = MASK_LETTER;
                } else if mask_mode == 2 {
                    seq[i] |= SEED_MASK;
                }
                ranges.push_front(i as i32);
            }
        }

        ranges
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_ranges(ranges: &Ranges) -> Vec<(i32, i32)> {
        let (a, b) = ranges.as_slices();
        a.iter().chain(b).copied().collect()
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn forced_scalar_mask_matches_runtime_simd() {
        if !super::super::tantan_simd::has_sse41_ssse3() {
            return;
        }
        let source: Vec<Letter> = (0..777)
            .map(|i| {
                if (180..310).contains(&i) {
                    [14, 14, 16, 14, 14, 16][i % 6]
                } else {
                    ((i * 13 + i / 7) % TRUE_AA as usize) as Letter
                }
            })
            .collect();
        let mut simd = source.clone();
        mask_tantan(&mut simd);

        struct Reset;
        impl Drop for Reset {
            fn drop(&mut self) {
                FORCE_SCALAR.with(|forced| forced.set(false));
            }
        }
        FORCE_SCALAR.with(|forced| forced.set(true));
        let _reset = Reset;
        let mut scalar = source;
        mask_tantan(&mut scalar);
        assert_eq!(simd, scalar);
    }

    #[test]
    fn test_forward_step_scalar_matches_definition() {
        let mut f = std::array::from_fn(|i| (i as f32 + 1.0) / 100.0);
        let d = std::array::from_fn(|i| (50 - i) as f32 / 1000.0);
        let e: Vec<f32> = (0..WINDOW).map(|i| 0.75 + i as f32 / 200.0).collect();
        let original_f = f;
        let mut b = 0.625;
        let b_old = b;
        let f2f = 0.95;
        let p_repeat_end = 0.05;
        let b2b = 0.995;
        let f_sum_prev = 1.25;

        let sum = forward_step_scalar(&mut f, &d, &e, &mut b, f2f, p_repeat_end, b2b, f_sum_prev);

        let mut expected_sum = 0.0f32;
        for i in 0..WINDOW {
            let expected = (original_f[i] * f2f + b_old * d[i]) * e[i];
            assert_eq!(f[i].to_bits(), expected.to_bits(), "f[{i}]");
            expected_sum += expected;
        }
        assert_eq!(sum.to_bits(), expected_sum.to_bits());
        assert_eq!(
            b.to_bits(),
            (b_old * b2b + f_sum_prev * p_repeat_end).to_bits()
        );
    }

    #[test]
    fn test_backward_step_scalar_matches_definition() {
        let mut f = std::array::from_fn(|i| (i as f32 + 1.0) / 80.0);
        let d = std::array::from_fn(|i| (50 - i) as f32 / 700.0);
        let e: Vec<f32> = (0..WINDOW).map(|i| 0.8 + i as f32 / 250.0).collect();
        let original_f = f;
        let mut b = 0.625;
        let b_old = b;
        let f2f = 0.95;
        let p_repeat_end = 0.05;
        let b2b = 0.995;

        let tsum = backward_step_scalar(&mut f, &d, &e, &mut b, f2f, p_repeat_end, b2b);

        let mut expected_sum = 0.0f32;
        for i in 0..WINDOW {
            let vf = original_f[i] * e[i];
            expected_sum += vf * d[i];
            let expected = vf * f2f + p_repeat_end * b_old;
            assert_eq!(f[i].to_bits(), expected.to_bits(), "f[{i}]");
        }
        assert_eq!(tsum.to_bits(), expected_sum.to_bits());
        assert_eq!(b.to_bits(), (b2b * b_old + expected_sum).to_bits());
    }

    #[test]
    fn test_mask_modes_and_ranges_match_cpp_contract() {
        let lr = vec![vec![1.0f32; AMINO_ACID_COUNT]; AMINO_ACID_COUNT];
        let original: Vec<Letter> = (0..64).map(|i| (i % TRUE_AA as usize) as Letter).collect();

        let mut table_only = original.clone();
        let table_ranges = mask(
            &mut table_only,
            &lr,
            P_REPEAT,
            P_REPEAT_END,
            REPEAT_GROWTH,
            -1.0,
            0,
        );
        assert_eq!(table_only, original);
        assert_eq!(flat_ranges(&table_ranges), vec![(0, 64)]);

        let mut hard = original.clone();
        let hard_ranges = mask(
            &mut hard,
            &lr,
            P_REPEAT,
            P_REPEAT_END,
            REPEAT_GROWTH,
            -1.0,
            1,
        );
        assert!(hard.iter().all(|&letter| letter == MASK_LETTER));
        assert_eq!(hard_ranges, table_ranges);

        let mut soft = original.clone();
        let soft_ranges = mask(
            &mut soft,
            &lr,
            P_REPEAT,
            P_REPEAT_END,
            REPEAT_GROWTH,
            -1.0,
            2,
        );
        assert!(soft.iter().all(|&letter| letter & SEED_MASK != 0));
        assert_eq!(soft_ranges, table_ranges);
    }

    #[test]
    fn test_mask_empty_returns_no_ranges() {
        let lr = vec![vec![1.0f32; AMINO_ACID_COUNT]; AMINO_ACID_COUNT];
        let ranges = mask(
            &mut [],
            &lr,
            P_REPEAT,
            P_REPEAT_END,
            REPEAT_GROWTH,
            DEFAULT_MIN_MASK_PROB,
            2,
        );
        assert!(ranges.is_empty());
    }

    #[test]
    fn malformed_letters_and_short_matrix_use_zero_emissions() {
        let mut seq = vec![0, DELIMITER_LETTER, 31, 1, 30, 2];
        let lr = vec![vec![1.0f32; 2]];
        let _ = mask(
            &mut seq,
            &lr,
            P_REPEAT,
            P_REPEAT_END,
            REPEAT_GROWTH,
            DEFAULT_MIN_MASK_PROB,
            0,
        );
        assert_eq!(seq, vec![0, DELIMITER_LETTER, 31, 1, 30, 2]);
    }

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
        assert!(
            masked_count > 0,
            "No positions masked in repetitive sequence"
        );
    }

    #[test]
    fn test_mask_q6gzx3_repeat() {
        // Q6GZX3 has a PPTPPT repeat near the end that C++ masks (~12 positions)
        let seq_str = b"MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQTCASGFCTSQPLCARIKKTQVCGLRYSSKGKDPLVSAEWDSRGAPYVRCTYDADLIDTQAQVDQFVSMFGESPSLAERYCMRGVKNTAGELVSRVSSDADPAGGWCRKWYSAHRGPDQDAALGSFCIKNPGAADCKCINRASDPVYQKVKTLHAYPDQCWYVPCAADVGELKMGTQRDTPTNCPTQVCQIVFNMLDDGSVTMDDVKNTINCDFSKYVPPPPPPKPTPPTPPTPPTPPTPPTPPTPPTPRPVHNRKVMFFVAGAVLVAILISTVRW";
        let letter_map: std::collections::HashMap<u8, i8> = [
            (b'A', 0),
            (b'R', 1),
            (b'N', 2),
            (b'D', 3),
            (b'C', 4),
            (b'Q', 5),
            (b'E', 6),
            (b'G', 7),
            (b'H', 8),
            (b'I', 9),
            (b'L', 10),
            (b'K', 11),
            (b'M', 12),
            (b'F', 13),
            (b'P', 14),
            (b'S', 15),
            (b'T', 16),
            (b'W', 17),
            (b'Y', 18),
            (b'V', 19),
        ]
        .iter()
        .cloned()
        .collect();
        let mut seq: Vec<Letter> = seq_str
            .iter()
            .map(|&c| *letter_map.get(&c).unwrap_or(&23))
            .collect();

        mask_tantan(&mut seq);
        let masked: Vec<usize> = seq
            .iter()
            .enumerate()
            .filter(|(_, &l)| l & SEED_MASK != 0)
            .map(|(i, _)| i)
            .collect();
        let diag: std::collections::HashMap<u8, i32> = [
            (b'A', 4),
            (b'R', 5),
            (b'N', 6),
            (b'D', 6),
            (b'C', 9),
            (b'Q', 5),
            (b'E', 5),
            (b'G', 6),
            (b'H', 8),
            (b'I', 4),
            (b'L', 4),
            (b'K', 5),
            (b'M', 5),
            (b'F', 6),
            (b'P', 7),
            (b'S', 4),
            (b'T', 5),
            (b'W', 11),
            (b'Y', 7),
            (b'V', 4),
        ]
        .iter()
        .cloned()
        .collect();
        let deficit: i32 = masked
            .iter()
            .map(|&i| diag.get(&seq_str[i]).copied().unwrap_or(0))
            .sum();
        eprintln!(
            "Q6GZX3 masked {} of {} positions, score deficit={}",
            masked.len(),
            seq.len(),
            deficit
        );
        eprintln!("  Masked positions: {:?}", &masked);
        eprintln!(
            "  Masked region: {}",
            masked
                .iter()
                .map(|&i| seq_str[i] as char)
                .collect::<String>()
        );
        // C++ masks ~12 positions in the PPTPPT repeat (around pos 259-295)
        // Rust should mask a similar number
        assert!(
            masked.len() > 0,
            "Q6GZX3 repeat region should have some masking"
        );
    }

    #[test]
    fn test_lambda_blosum62() {
        let sm = &crate::stats::matrices::BLOSUM62;
        let lambda = compute_lambda_flat(&sm.scores, crate::basic::value::AMINO_ACID_COUNT, 20);
        // BLOSUM62 lambda should be ~0.324
        assert!(
            (lambda - 0.324).abs() < 0.01,
            "Lambda should be ~0.324, got {:.6}",
            lambda
        );
    }

    #[test]
    fn test_mask_d1ivsa4_no_masking() {
        // d1ivsa4 (426aa) — C++ masks 0 positions. Rust should too.
        // Read actual d1ivsa4 sequence from test file
        let fasta_path = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/5.faa");
        let records = crate::data::fasta::read_fasta_file(
            std::path::Path::new(fasta_path),
            crate::basic::value::SequenceType::AminoAcid,
        )
        .unwrap();
        let rec = records.iter().find(|r| r.id.contains("d1ivsa4")).unwrap();
        let mut seq = rec.sequence.clone();
        assert_eq!(seq.len(), 426);

        mask_tantan(&mut seq);
        let masked: Vec<usize> = seq
            .iter()
            .enumerate()
            .filter(|(_, &l)| l & SEED_MASK != 0)
            .map(|(i, _)| i)
            .collect();
        eprintln!("d1ivsa4 masked {} of {} positions", masked.len(), seq.len());
        if !masked.is_empty() {
            eprintln!("  Masked: {:?}", &masked[..masked.len().min(20)]);
        }
        // C++ masks 0 positions for this sequence
        assert_eq!(
            masked.len(),
            0,
            "d1ivsa4 should have 0 masked positions (C++ masks 0), got {}",
            masked.len()
        );
    }
}
