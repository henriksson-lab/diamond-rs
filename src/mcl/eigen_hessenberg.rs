//! Eigen 3.3.8 real Hessenberg and elementary Householder operations.
//!
//! Translated from Eigen/src/Householder/Householder.h (Benoit Jacob,
//! copyright 2010; Gael Guennebaud, copyright 2009), HouseholderSequence.h,
//! Eigenvalues/HessenbergDecomposition.h, and Core/{Redux.h,
//! products/GeneralMatrixVector.h}. Those files are licensed under MPL-2.0:
//! https://mozilla.org/MPL/2.0/ . This file is subject to that license.
//!
//! Arithmetic traversal follows the pinned GNU reference's baseline SSE
//! four-float packets, emulated in scalar Rust on every platform. Matrix
//! allocation alignment is modeled as Eigen's aligned column-major storage;
//! Rust allocation addresses therefore do not change reduction order.

use nalgebra::DMatrix;

fn reduce_packet(p: [f32; 4]) -> f32 {
    (p[0] + p[2]) + (p[1] + p[3])
}

// abs2() has no DirectAccessBit, so Redux first_default_aligned returns zero.
pub(super) fn squared_norm(v: &[f32]) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    let end = v.len() / 4 * 4;
    if end == 0 {
        let mut sum = v[0] * v[0];
        for &x in &v[1..] {
            sum += x * x;
        }
        return sum;
    }
    let packet = |i: usize| {
        [
            v[i] * v[i],
            v[i + 1] * v[i + 1],
            v[i + 2] * v[i + 2],
            v[i + 3] * v[i + 3],
        ]
    };
    let mut p = packet(0);
    if end > 4 {
        let mut q = packet(4);
        let end2 = v.len() / 8 * 8;
        for i in (8..end2).step_by(8) {
            let x = packet(i);
            let y = packet(i + 4);
            for lane in 0..4 {
                p[lane] += x[lane];
                q[lane] += y[lane];
            }
        }
        for lane in 0..4 {
            p[lane] += q[lane];
        }
        if end > end2 {
            let x = packet(end2);
            for lane in 0..4 {
                p[lane] += x[lane];
            }
        }
    }
    let mut sum = reduce_packet(p);
    for &x in &v[end..] {
        sum += x * x;
    }
    sum
}

/// Essential part, coefficient tau, and transformed leading coefficient beta.
pub(super) fn reflector(v: &[f32]) -> (Vec<f32>, f32, f32) {
    assert!(!v.is_empty());
    let tail_norm = squared_norm(&v[1..]);
    let c0 = v[0];
    if tail_norm <= f32::MIN_POSITIVE {
        (vec![0.0; v.len() - 1], 0.0, c0)
    } else {
        let mut beta = (c0 * c0 + tail_norm).sqrt();
        if c0 >= 0.0 {
            beta = -beta;
        }
        let essential = v[1..].iter().map(|&x| x / (c0 - beta)).collect();
        (essential, (beta - c0) / beta, beta)
    }
}

fn first_aligned(offset: usize, size: usize) -> usize {
    ((4 - offset % 4) % 4).min(size)
}

pub(super) fn apply_left(
    a: &mut DMatrix<f32>,
    r: usize,
    c: usize,
    rows: usize,
    cols: usize,
    e: &[f32],
    tau: f32,
) {
    apply_left_at(a, r, c, rows, cols, e, tau, 0);
}

fn apply_left_at(
    a: &mut DMatrix<f32>,
    r: usize,
    c: usize,
    rows: usize,
    cols: usize,
    e: &[f32],
    tau: f32,
    e_offset: usize,
) {
    assert_eq!(e.len(), rows.saturating_sub(1));
    if rows == 1 {
        for col in c..c + cols {
            a[(r, col)] *= 1.0 - tau;
        }
    } else if tau != 0.0 {
        // essential.adjoint()*bottom: row-major GEMV on bottom.transpose().
        let depth = rows - 1;
        let mut aligned_start = first_aligned(e_offset, depth);
        let lhs_offset = (c * a.nrows() + r + 1) % 4;
        if first_aligned(lhs_offset, depth) == depth || first_aligned(e_offset, cols) == cols {
            aligned_start = 0;
        }
        let aligned_end = if aligned_start == 0
            && (first_aligned(lhs_offset, depth) == depth || first_aligned(e_offset, cols) == cols)
        {
            0
        } else {
            aligned_start + (depth - aligned_start) / 4 * 4
        };
        let mut tmp = vec![0.0; cols];
        for col in 0..cols {
            let mut sum = 0.0;
            for i in 0..aligned_start {
                sum += a[(r + 1 + i, c + col)] * e[i];
            }
            let mut p = [0.0; 4];
            for i in (aligned_start..aligned_end).step_by(4) {
                for lane in 0..4 {
                    p[lane] += a[(r + 1 + i + lane, c + col)] * e[i + lane];
                }
            }
            if aligned_end > aligned_start {
                sum += reduce_packet(p);
            }
            for i in aligned_end..depth {
                sum += a[(r + 1 + i, c + col)] * e[i];
            }
            tmp[col] = sum + a[(r, c + col)];
        }
        for col in 0..cols {
            a[(r, c + col)] -= tau * tmp[col];
        }
        // OuterProduct evaluates the scaled left vector once, then multiplies
        // by each RHS coefficient in column-major destination order.
        let scaled: Vec<_> = e.iter().map(|&x| tau * x).collect();
        for col in 0..cols {
            for i in 0..depth {
                a[(r + 1 + i, c + col)] -= tmp[col] * scaled[i];
            }
        }
    }
}

/// Column-major real GEMV with Eigen's aligned zero destination and the
/// source matrix's physical leading dimension. Used by Householder updates
/// and EigenSolver's eigenvector backtransformation.
pub(super) fn column_product(
    a: &DMatrix<f32>,
    r: usize,
    c: usize,
    rows: usize,
    cols: usize,
    v: &[f32],
) -> Vec<f32> {
    assert_eq!(v.len(), cols);
    let n = a.nrows();
    let lhs_offset = (c * n + r) % 4;
    let lhs_alignment = first_aligned(lhs_offset, rows);
    let step = (4 - n % 4) % 4;
    let mut skip = 0;
    if lhs_alignment != rows {
        while skip < 4 && (lhs_alignment + step * skip) % 4 != 0 {
            skip += 1;
        }
        if skip == 4 {
            skip = 0;
        } else {
            skip = skip.min(cols);
        }
    }
    let aligned_end = if lhs_alignment == rows {
        0
    } else {
        rows / 4 * 4
    };
    let first_pattern = lhs_alignment != rows && step != 0 && step != 2;
    let offset1 = if first_pattern && step == 1 { 3 } else { 1 };
    let offset3 = if first_pattern && step == 1 { 1 } else { 3 };
    let bound = skip + (cols - skip) / 4 * 4;
    let mut tmp = vec![0.0; rows];
    for k in (skip..bound).step_by(4) {
        for row in 0..rows {
            let x0 = a[(r + row, c + k)] * v[k];
            let x1 = a[(r + row, c + k + offset1)] * v[k + offset1];
            let x2 = a[(r + row, c + k + 2)] * v[k + 2];
            let x3 = a[(r + row, c + k + offset3)] * v[k + offset3];
            // FirstAligned's two-packet peel uses sequential pmadd;
            // the remaining SIMD packets use balanced pair sums.
            let peeled =
                first_pattern && aligned_end >= 12 && row < ((aligned_end - 12) / 8 + 1) * 8;
            if row < aligned_end && !peeled {
                tmp[row] += (x0 + x1) + (x2 + x3);
            } else {
                tmp[row] += x0;
                tmp[row] += x1;
                tmp[row] += x2;
                tmp[row] += x3;
            }
        }
    }
    for range in [bound..cols, 0..skip] {
        for k in range {
            for row in 0..rows {
                tmp[row] += a[(r + row, c + k)] * v[k];
            }
        }
    }
    tmp
}

pub(super) fn apply_right(
    a: &mut DMatrix<f32>,
    r: usize,
    c: usize,
    rows: usize,
    cols: usize,
    e: &[f32],
    tau: f32,
) {
    assert_eq!(e.len(), cols.saturating_sub(1));
    if cols == 1 {
        for row in r..r + rows {
            a[(row, c)] *= 1.0 - tau;
        }
    } else if tau != 0.0 {
        let mut tmp = column_product(a, r, c + 1, rows, cols - 1, e);
        for row in 0..rows {
            tmp[row] += a[(r + row, c)];
        }
        for row in 0..rows {
            a[(r + row, c)] -= tau * tmp[row];
        }
        let scaled: Vec<_> = tmp.iter().map(|&x| tau * x).collect();
        for col in 0..e.len() {
            for row in 0..rows {
                a[(r + row, c + 1 + col)] -= e[col] * scaled[row];
            }
        }
    }
}

pub(super) fn hessenberg(mut input: DMatrix<f32>) -> (DMatrix<f32>, DMatrix<f32>) {
    assert_eq!(input.nrows(), input.ncols());
    let n = input.nrows();
    let mut coefficients = vec![0.0; n.saturating_sub(1)];
    for i in 0..n.saturating_sub(1) {
        let remaining = n - i - 1;
        let vector: Vec<_> = (i + 1..n).map(|row| input[(row, i)]).collect();
        let (essential, tau, beta) = reflector(&vector);
        input[(i + 1, i)] = beta;
        for (j, &x) in essential.iter().enumerate() {
            input[(i + 2 + j, i)] = x;
        }
        coefficients[i] = tau;
        apply_left_at(
            &mut input,
            i + 1,
            i + 1,
            remaining,
            remaining,
            &essential,
            tau,
            i * n + i + 2,
        );
        apply_right(&mut input, 0, i + 1, n, remaining, &essential, tau);
    }
    // HouseholderSequence::evalTo, separate destination, shift=1. It applies
    // reflectors in reverse order to successively smaller identity corners.
    let mut q = DMatrix::identity(n, n);
    for i in (0..n.saturating_sub(1)).rev() {
        let essential: Vec<_> = (i + 2..n).map(|row| input[(row, i)]).collect();
        apply_left_at(
            &mut q,
            i + 1,
            i + 1,
            n - i - 1,
            n - i - 1,
            &essential,
            coefficients[i],
            i * n + i + 2,
        );
    }
    for col in 0..n {
        for row in col + 2..n {
            input[(row, col)] = 0.0;
        }
    }
    (input, q)
}
