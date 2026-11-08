//! Reachable complex-f32 operations used by DIAMOND's dense MCL expansion.
//!
//! Ported from Eigen 3.3.8: LU/PartialPivLU.h, Core/Redux.h,
//! Core/products/{GeneralMatrixMatrix,GeneralBlockPanelKernel,
//! TriangularSolverMatrix}.h and Core/arch/SSE/Complex.h.
//! Eigen copyright (C) 2008-2016 Gael Guennebaud and Benoit Jacob;
//! licensed under MPL-2.0 (https://www.mozilla.org/MPL/2.0/).
//! DMatrix supplies column-major storage only; its arithmetic/decompositions
//! are deliberately not used. Packet operations are expressed as scalar lanes
//! so the pinned SSE evaluation order also applies on non-x86 hosts.
use nalgebra::{Complex, DMatrix};
type C = Complex<f32>;

#[inline]
fn mul(a: C, b: C) -> C {
    C::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

/// Eigen's SSE complex packet quotient, including a real divisor converted
/// to complex by DenseBase::operator/=(const Scalar&).
pub(super) fn source_complex_packet_div(a: C, b: C) -> C {
    let numerator = C::new(a.re * b.re + a.im * b.im, -(a.re * b.im) + a.im * b.re);
    let denominator = b.re * b.re + b.im * b.im;
    C::new(numerator.re / denominator, numerator.im / denominator)
}

/// GCC's scalar complex-float division (__divsc3). The pinned libgcc
/// computes products, sum and quotient in f64, rounding only the two results
/// to f32. This differs from both num_complex's quotient and SSE packet pdiv.
/// The widened denominator covers every finite f32 without scaling.
pub(super) fn source_complex_div(a: C, b: C) -> C {
    let mut c = b.re;
    let mut d = b.im;
    let maximum = c.abs().max(d.abs());
    let denominator = (c as f64) * (c as f64) + (d as f64) * (d as f64);
    let mut x = (((a.re as f64) * (c as f64) + (a.im as f64) * (d as f64)) / denominator) as f32;
    let mut y = (((a.im as f64) * (c as f64) - (a.re as f64) * (d as f64)) / denominator) as f32;
    if x.is_nan() && y.is_nan() {
        if denominator == 0.0 && (!a.re.is_nan() || !a.im.is_nan()) {
            x = f32::INFINITY.copysign(b.re) * a.re;
            y = f32::INFINITY.copysign(b.re) * a.im;
        } else if (a.re.is_infinite() || a.im.is_infinite()) && c.is_finite() && d.is_finite() {
            let ar = if a.re.is_infinite() { 1.0f32 } else { 0.0 }.copysign(a.re);
            let ai = if a.im.is_infinite() { 1.0f32 } else { 0.0 }.copysign(a.im);
            x = f32::INFINITY * (ar * c + ai * d);
            y = f32::INFINITY * (ai * c - ar * d);
        } else if maximum.is_infinite() && a.re.is_finite() && a.im.is_finite() {
            c = if b.re.is_infinite() { 1.0f32 } else { 0.0 }.copysign(b.re);
            d = if b.im.is_infinite() { 1.0f32 } else { 0.0 }.copysign(b.im);
            x = 0.0 * (a.re * c + a.im * d);
            y = 0.0 * (a.im * c - a.re * d);
        }
    }
    C::new(x, y)
}

/// libstdc++ complex<float> real-exponent overload: positive-real shortcut,
/// otherwise exp(r*real(log(z))) * (cos(r*arg(z)), sin(r*arg(z))).
pub(super) fn source_complex_pow(z: C, power: f32) -> C {
    if z.im == 0.0 && z.re > 0.0 {
        return C::new(z.re.powf(power), 0.0);
    }
    let logarithm_real = source_log_abs(z);
    let logarithm_imag = z.im.atan2(z.re);
    let radius = (power * logarithm_real).exp();
    let angle = power * logarithm_imag;
    C::new(radius * angle.cos(), radius * angle.sin())
}

// Reachable clogf real-part reduction, GNU libc math/s_clog_template.c
// and sysdeps/ieee754/dbl-64/x2y2m1f.c, copyright FSF 1997-2022,
// LGPL-2.1-or-later. The log1p branches avoid cancellation near |z|=1.
// https://git.zx2c4.com/glibc/tree/math/s_clog_template.c
fn source_log_abs(z: C) -> f32 {
    if z.re.is_nan() || z.im.is_nan() {
        return if z.re.is_infinite() || z.im.is_infinite() {
            f32::INFINITY
        } else {
            f32::NAN
        };
    }
    let mut x = z.re.abs().max(z.im.abs());
    let mut y = z.re.abs().min(z.im.abs());
    if x == 0.0 {
        return f32::NEG_INFINITY;
    }
    let mut scale = 0;
    if x > f32::MAX / 2.0 {
        scale = -1;
        x *= 0.5;
        y = if y >= f32::MIN_POSITIVE * 2.0 {
            y * 0.5
        } else {
            0.0
        };
    } else if x < f32::MIN_POSITIVE && y < f32::MIN_POSITIVE {
        scale = 24;
        x *= 16777216.0;
        y *= 16777216.0;
    }
    if x == 1.0 && scale == 0 {
        return (y * y).ln_1p() / 2.0;
    }
    if x > 1.0 && x < 2.0 && y < 1.0 && scale == 0 {
        let mut value = (x - 1.0) * (x + 1.0);
        if y >= f32::EPSILON {
            value += y * y;
        }
        return value.ln_1p() / 2.0;
    }
    if x < 1.0 && x >= 0.5 && scale == 0 {
        if y < f32::EPSILON / 2.0 {
            return ((x - 1.0) * (x + 1.0)).ln_1p() / 2.0;
        }
        if x * x + y * y >= 0.5 {
            let (x, y) = (x as f64, y as f64);
            return (((x - 1.0) * (x + 1.0) + y * y) as f32).ln_1p() / 2.0;
        }
    }
    x.hypot(y).ln() - (scale as f32) * std::f32::consts::LN_2
}

// Eigen's single-thread complex GEMM uses two-complex packets. Accumulating
// the real and imaginary RHS contributions separately, then assembling once,
// is observably different from accumulating individual complex products.
fn packet_dot(a: impl Fn(usize) -> C, b: impl Fn(usize) -> C, count: usize) -> C {
    let (mut rr, mut ir, mut ri, mut ii) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    for k in 0..count {
        let a = a(k);
        let b = b(k);
        rr = a.re * b.re + rr;
        ir = a.im * b.re + ir;
        ri = a.re * b.im + ri;
        ii = a.im * b.im + ii;
    }
    C::new(rr - ii, ir + ri)
}
fn scalar_dot(a: impl Fn(usize) -> C, b: impl Fn(usize) -> C, count: usize, zero: bool) -> C {
    if count == 0 {
        return C::new(0.0, 0.0);
    }
    let mut value = if zero {
        C::new(0.0, 0.0)
    } else {
        mul(a(0), b(0))
    };
    for k in usize::from(!zero)..count {
        value += mul(a(k), b(k));
    }
    value
}

// Only the K split changes each coefficient's reduction. The M/N blocks are
// multiples of packet size/nr and preserve the same packet/scalar row split.
fn depth_block(depth: usize, rows: usize, cols: usize, mixed: bool) -> usize {
    let (divisor, subtract) = if mixed {
        (4 * 8 + 4 * 4, 4 * 4 * 8)
    } else {
        (2 * 8 + 4 * 8, 2 * 4 * 8)
    };
    depth_block_dimensions(depth, rows, cols, divisor, subtract)
}
fn depth_block_dimensions(
    depth: usize,
    rows: usize,
    cols: usize,
    divisor: usize,
    subtract: usize,
) -> usize {
    if depth.max(rows).max(cols) < 48 {
        return depth.max(1);
    }
    let l1 = source_l1_cache();
    let max_kc = ((l1.saturating_sub(subtract) / divisor) & !7).max(1);
    if depth <= max_kc {
        depth.max(1)
    } else if depth % max_kc == 0 {
        max_kc
    } else {
        max_kc - 8 * ((max_kc - 1 - depth % max_kc) / (8 * (depth / max_kc + 1)))
    }
}
fn source_l1_cache() -> usize {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        #[cfg(target_arch = "x86")]
        use std::arch::x86::__cpuid_count;
        #[cfg(target_arch = "x86_64")]
        use std::arch::x86_64::__cpuid_count;
        // CPUID leaf4 is the same deterministic Intel cache query used by
        // Eigen's queryCacheSizes. Unsupported leaves return no data cache.
        for index in 0..32 {
            let result = __cpuid_count(4, index);
            let kind = result.eax & 31;
            if kind == 0 {
                break;
            }
            if kind == 1 && (result.eax >> 5) & 7 == 1 {
                return (((result.ebx >> 22) + 1) as usize)
                    * ((((result.ebx >> 12) & 1023) + 1) as usize)
                    * (((result.ebx & 4095) + 1) as usize)
                    * ((result.ecx + 1) as usize);
            }
        }
    }
    32 * 1024
}

/// GeneralMatrixMatrix.h's dynamic real-f32 product, including GEBP's
/// four-accumulator scalar-row fringe. A scalar triple loop differs exactly at
/// these fringe rows, even when every operand is an ordinary finite protein
/// graph weight. This is not a special case for a particular matrix size.
pub(super) fn source_real_product(a: &DMatrix<f32>, b: &DMatrix<f32>) -> DMatrix<f32> {
    assert_eq!(a.ncols(), b.nrows());
    let (rows, cols, depth) = (a.nrows(), b.ncols(), a.ncols());
    let small = rows + cols + depth < 20 && depth > 0;
    let kc = depth_block_dimensions(depth, rows, cols, 8 * 4 + 4 * 4, 8 * 4 * 4);
    let mut result = DMatrix::zeros(rows, cols);
    for col in 0..cols {
        for row in 0..rows {
            if small {
                let mut value = a[(row, 0)] * b[(0, col)];
                for k in 1..depth {
                    value += a[(row, k)] * b[(k, col)];
                }
                result[(row, col)] = value;
            } else {
                for start in (0..depth).step_by(kc) {
                    let count = (depth - start).min(kc);
                    let mut value = 0.0f32;
                    if row >= rows / 4 * 4 && col < cols / 4 * 4 {
                        // SwappedTraits::LhsProgress==4: columns are vector
                        // lanes, and four independent K streams hide latency.
                        let end = count / 4 * 4;
                        let mut accum = [0.0f32; 4];
                        for k in (0..end).step_by(4) {
                            for lane in 0..4 {
                                accum[lane] +=
                                    a[(row, start + k + lane)] * b[(start + k + lane, col)];
                            }
                        }
                        value = (accum[0] + accum[1]) + (accum[2] + accum[3]);
                        for k in end..count {
                            value += a[(row, start + k)] * b[(start + k, col)];
                        }
                    } else {
                        for k in 0..count {
                            value += a[(row, start + k)] * b[(start + k, col)];
                        }
                    }
                    result[(row, col)] += value;
                }
            }
        }
    }
    result
}

pub(super) fn source_complex_product(a: &DMatrix<C>, b: &DMatrix<C>) -> DMatrix<C> {
    assert_eq!(a.ncols(), b.nrows());
    let (rows, cols, depth) = (a.nrows(), b.ncols(), a.ncols());
    let small = rows + cols + depth < 20 && depth > 0;
    let kc = depth_block(depth, rows, cols, false);
    let mut result = DMatrix::zeros(rows, cols);
    for col in 0..cols {
        for row in 0..rows {
            if small {
                result[(row, col)] = scalar_dot(|k| a[(row, k)], |k| b[(k, col)], depth, false);
            } else {
                for start in (0..depth).step_by(kc) {
                    let count = (depth - start).min(kc);
                    let value = if row < rows / 2 * 2 {
                        packet_dot(|k| a[(row, start + k)], |k| b[(start + k, col)], count)
                    } else {
                        scalar_dot(
                            |k| a[(row, start + k)],
                            |k| b[(start + k, col)],
                            count,
                            true,
                        )
                    };
                    result[(row, col)] += mul(value, C::new(1.0, 0.0));
                }
            }
        }
    }
    result
}

pub(super) fn source_complex_real_product(a: &DMatrix<C>, b: &DMatrix<f32>) -> DMatrix<C> {
    assert_eq!(a.ncols(), b.nrows());
    let (rows, cols, depth) = (a.nrows(), b.ncols(), a.ncols());
    let small = rows + cols + depth < 20 && depth > 0;
    let kc = depth_block(depth, rows, cols, true);
    let mut result = DMatrix::zeros(rows, cols);
    for col in 0..cols {
        for row in 0..rows {
            if small {
                let mut value = a[(row, 0)] * b[(0, col)];
                for k in 1..depth {
                    value += a[(row, k)] * b[(k, col)];
                }
                result[(row, col)] = value;
            } else {
                for start in (0..depth).step_by(kc) {
                    let mut value = C::new(0.0, 0.0);
                    for k in start..(start + kc).min(depth) {
                        value += a[(row, k)] * b[(k, col)];
                    }
                    result[(row, col)] += value;
                }
            }
        }
    }
    result
}

struct PartialLu {
    lu: DMatrix<C>,
    swaps: Vec<usize>,
    sign: f32,
}
impl PartialLu {
    fn new(matrix: &DMatrix<C>) -> Self {
        assert!(matrix.is_square());
        let mut lu = matrix.clone();
        let (swaps, transpositions) = factor(&mut lu, 0, 0, matrix.nrows(), matrix.ncols(), 256);
        Self {
            lu,
            swaps,
            sign: if transpositions % 2 == 0 { 1.0 } else { -1.0 },
        }
    }
}

fn swap_region(a: &mut DMatrix<C>, r0: usize, r1: usize, col: usize, count: usize) {
    if r0 != r1 {
        for j in col..col + count {
            a.swap((r0, j), (r1, j));
        }
    }
}
fn factor(
    a: &mut DMatrix<C>,
    row: usize,
    col: usize,
    rows: usize,
    cols: usize,
    maximum: usize,
) -> (Vec<usize>, usize) {
    let size = rows.min(cols);
    let mut swaps = vec![0; size];
    let mut transpositions = 0;
    if size <= 16 {
        for k in 0..size {
            let mut pivot = k;
            let mut biggest = a[(row + k, col + k)].re.hypot(a[(row + k, col + k)].im);
            for i in k + 1..rows {
                let value = a[(row + i, col + k)];
                let score = value.re.hypot(value.im);
                if score > biggest {
                    biggest = score;
                    pivot = i;
                }
            }
            swaps[k] = pivot;
            if biggest != 0.0 {
                swap_region(a, row + k, row + pivot, col, cols);
                transpositions += usize::from(k != pivot);
                let denominator = a[(row + k, col + k)];
                // Eigen aligns the dynamic column tail before its two-complex
                // packet loop; a prefix/suffix scalar follows GCC division.
                let begin = k + 1;
                let aligned = begin + (((col + k) * a.nrows() + row + begin) & 1);
                let end = aligned + (rows.saturating_sub(aligned) / 2) * 2;
                for i in begin..rows {
                    let numerator = a[(row + i, col + k)];
                    a[(row + i, col + k)] = if i >= aligned && i < end {
                        source_complex_packet_div(numerator, denominator)
                    } else {
                        source_complex_div(numerator, denominator)
                    };
                }
            }
            for j in k + 1..cols {
                for i in k + 1..rows {
                    let value = mul(a[(row + i, col + k)], a[(row + k, col + j)]);
                    a[(row + i, col + j)] -= value;
                }
            }
        }
    } else {
        let block = ((size / 8) / 16 * 16).max(8).min(maximum);
        for k in (0..size).step_by(block) {
            let bs = (size - k).min(block);
            let (panel, count) = factor(a, row + k, col + k, rows - k, bs, 16);
            transpositions += count;
            for i in 0..bs {
                swaps[k + i] = k + panel[i];
                swap_region(a, row + k + i, row + swaps[k + i], col, k);
            }
            if rows - k - bs != 0 {
                for i in 0..bs {
                    swap_region(
                        a,
                        row + k + i,
                        row + swaps[k + i],
                        col + k + bs,
                        cols - k - bs,
                    );
                }
                // The recursive panel uses the same four-row triangular
                // solver as inverse(), including its separate GEBP updates.
                let diagonal = a.view((row + k, col + k), (bs, bs)).into_owned();
                let mut rhs = a
                    .view((row + k, col + k + bs), (bs, cols - k - bs))
                    .into_owned();
                triangular_solve(&diagonal, &mut rhs, true);
                a.view_mut((row + k, col + k + bs), (bs, cols - k - bs))
                    .copy_from(&rhs);
                let update = source_complex_product(
                    &a.view((row + k + bs, col + k), (rows - k - bs, bs))
                        .into_owned(),
                    &a.view((row + k, col + k + bs), (bs, cols - k - bs))
                        .into_owned(),
                );
                for j in k + bs..cols {
                    for i in k + bs..rows {
                        a[(row + i, col + j)] -= update[(i - k - bs, j - k - bs)];
                    }
                }
            }
        }
    }
    (swaps, transpositions)
}

pub(super) fn source_complex_determinant(matrix: &DMatrix<C>) -> C {
    let decomposition = PartialLu::new(matrix);
    if matrix.nrows() == 0 {
        return C::new(1.0, 0.0);
    }
    let mut determinant = decomposition.lu[(0, 0)];
    for i in 1..matrix.nrows() {
        determinant = mul(determinant, decomposition.lu[(i, i)]);
    }
    mul(C::new(decomposition.sign, 0.0), determinant)
}

/// Source does not test invertibility here; the determinant gate belongs to
/// get_exp. Preserve nonfinite arithmetic for a singular matrix passed directly.
pub(super) fn source_complex_inverse(matrix: &DMatrix<C>) -> DMatrix<C> {
    let decomposition = PartialLu::new(matrix);
    let n = matrix.nrows();
    let mut output = DMatrix::<C>::identity(n, n);
    for (i, &pivot) in decomposition.swaps.iter().enumerate() {
        output.swap_rows(i, pivot);
    }
    triangular_solve(&decomposition.lu, &mut output, true);
    triangular_solve(&decomposition.lu, &mut output, false);
    output
}

fn subtract_block(
    lu: &DMatrix<C>,
    output: &mut DMatrix<C>,
    target: usize,
    count: usize,
    start: usize,
    depth: usize,
) {
    for j in 0..output.ncols() {
        for row in target..target + count {
            let value = if row - target < count / 2 * 2 {
                packet_dot(|k| lu[(row, start + k)], |k| output[(start + k, j)], depth)
            } else {
                scalar_dot(
                    |k| lu[(row, start + k)],
                    |k| output[(start + k, j)],
                    depth,
                    true,
                )
            };
            output[(row, j)] += mul(C::new(-1.0, 0.0), value);
        }
    }
}

fn triangular_solve(lu: &DMatrix<C>, output: &mut DMatrix<C>, lower: bool) {
    let n = lu.nrows();
    let kc = depth_block(n, n, output.ncols(), false);
    let mut consumed = 0;
    while consumed < n {
        let region_size = (n - consumed).min(kc);
        let region_start = if lower {
            consumed
        } else {
            n - consumed - region_size
        };
        let region_end = region_start + region_size;
        let mut solved = 0;
        while solved < region_size {
            let width = (region_size - solved).min(4);
            let start = if lower {
                region_start + solved
            } else {
                region_end - solved - width
            };
            for k in 0..width {
                let i = if lower {
                    start + k
                } else {
                    start + width - k - 1
                };
                let inverse = if lower {
                    C::new(1.0, 0.0)
                } else {
                    source_complex_div(C::new(1.0, 0.0), lu[(i, i)])
                };
                for j in 0..output.ncols() {
                    output[(i, j)] = mul(output[(i, j)], inverse);
                    let b = output[(i, j)];
                    let range = if lower {
                        i + 1..start + width
                    } else {
                        start..i
                    };
                    for row in range {
                        output[(row, j)] -= mul(b, lu[(row, i)]);
                    }
                }
            }
            let (target, count) = if lower {
                (start + width, region_end - start - width)
            } else {
                (region_start, start - region_start)
            };
            subtract_block(lu, output, target, count, start, width);
            solved += width;
        }
        let (target, count) = if lower {
            (region_end, n - region_end)
        } else {
            (0, region_start)
        };
        subtract_block(lu, output, target, count, region_start, region_size);
        consumed += region_size;
    }
}
