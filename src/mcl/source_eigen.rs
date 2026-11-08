// Numerical eigensolver translated from the recovered Eigen 3.3.8 headers.
// Copyright (C) 2008 Gael Guennebaud; 2010,2012 Jitse Niesen.
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the Mozilla Public License v. 2.0.
// https://mozilla.org/MPL/2.0/
//
// Preserve the real Hessenberg/Francis-Schur path used by Source EigenSolver,
// rather than substituting a complex Schur decomposition with different roots.
use super::eigen_complex::{source_complex_div, source_complex_packet_div};
use nalgebra::{Complex, DMatrix};

type Mat = DMatrix<f32>;
type C = Complex<f32>;

use super::eigen_hessenberg::{
    apply_left as left, apply_right as right, column_product, hessenberg, reflector,
};

fn givens(p: f32, q: f32) -> (f32, f32) {
    if q == 0. {
        (if p < 0. { -1. } else { 1. }, 0.)
    } else if p == 0. {
        (0., if q < 0. { 1. } else { -1. })
    } else if p.abs() > q.abs() {
        let t = q / p;
        let mut u = (1. + t * t).sqrt();
        if p < 0. {
            u = -u;
        }
        let c = 1. / u;
        (c, -t * c)
    } else {
        let t = p / q;
        let mut u = (1. + t * t).sqrt();
        if q < 0. {
            u = -u;
        }
        let s = -1. / u;
        (-t * s, s)
    }
}
fn rotate_rows(a: &mut Mat, p: usize, q: usize, start: usize, c: f32, s: f32) {
    for j in start..a.ncols() {
        let x = a[(p, j)];
        let y = a[(q, j)];
        a[(p, j)] = c * x + s * y;
        a[(q, j)] = -s * x + c * y;
    }
}
fn rotate_cols(a: &mut Mat, p: usize, q: usize, rows: usize, c: f32, s: f32) {
    for i in 0..rows {
        let x = a[(i, p)];
        let y = a[(i, q)];
        a[(i, p)] = c * x - s * y;
        a[(i, q)] = s * x + c * y;
    }
}
fn schur(a: &Mat) -> (Mat, Mat) {
    let n = a.nrows();
    let scale = a.iter().map(|x| x.abs()).fold(0f32, f32::max);
    if scale < f32::MIN_POSITIVE {
        return (Mat::zeros(n, n), Mat::identity(n, n));
    }
    let (t, u) = hessenberg(a / scale);
    let (mut t, mut u) = (t, u);
    let mut norm = 0.;
    for j in 0..n {
        let mut s = 0.;
        for i in 0..n.min(j + 2) {
            s += t[(i, j)].abs();
        }
        norm += s;
    }
    let zero = (norm * (f32::EPSILON * f32::EPSILON)).max(f32::MIN_POSITIVE);
    let mut end = n as isize - 1;
    let mut iteration = 0;
    let mut total = 0;
    let mut exshift = 0.;
    while end >= 0 && norm != 0. {
        let iu = end as usize;
        let mut il = iu;
        while il > 0 {
            let s = ((t[(il - 1, il - 1)].abs() + t[(il, il)].abs()) * f32::EPSILON).max(zero);
            if t[(il, il - 1)].abs() <= s {
                break;
            }
            il -= 1;
        }
        if il == iu {
            t[(iu, iu)] += exshift;
            if iu > 0 {
                t[(iu, iu - 1)] = 0.;
            }
            end -= 1;
            iteration = 0;
        } else if il + 1 == iu {
            let p = 0.5 * (t[(iu - 1, iu - 1)] - t[(iu, iu)]);
            let q = p * p + t[(iu, iu - 1)] * t[(iu - 1, iu)];
            t[(iu, iu)] += exshift;
            t[(iu - 1, iu - 1)] += exshift;
            if q >= 0. {
                let z = q.abs().sqrt();
                let (c, s) = givens(if p >= 0. { p + z } else { p - z }, t[(iu, iu - 1)]);
                rotate_rows(&mut t, iu - 1, iu, iu - 1, c, -s);
                rotate_cols(&mut t, iu - 1, iu, iu + 1, c, s);
                t[(iu, iu - 1)] = 0.;
                rotate_cols(&mut u, iu - 1, iu, n, c, s);
            }
            if iu > 1 {
                t[(iu - 1, iu - 2)] = 0.;
            }
            end -= 2;
            iteration = 0;
        } else {
            let mut shift = [
                t[(iu, iu)],
                t[(iu - 1, iu - 1)],
                t[(iu, iu - 1)] * t[(iu - 1, iu)],
            ];
            if iteration == 10 {
                exshift += shift[0];
                for i in 0..=iu {
                    t[(i, i)] -= shift[0];
                }
                let s = t[(iu, iu - 1)].abs() + t[(iu - 1, iu - 2)].abs();
                shift = [0.75 * s, 0.75 * s, (-0.4375 * s) * s];
            }
            if iteration == 30 {
                let mut s = (shift[1] - shift[0]) / 2.;
                s = s * s + shift[2];
                if s > 0. {
                    s = s.sqrt();
                    if shift[1] < shift[0] {
                        s = -s;
                    }
                    s += (shift[1] - shift[0]) / 2.;
                    s = shift[0] - shift[2] / s;
                    exshift += s;
                    for i in 0..=iu {
                        t[(i, i)] -= s;
                    }
                    shift = [0.964; 3];
                }
            }
            iteration += 1;
            total += 1;
            if total > 40 * n {
                break;
            }
            let mut im = iu - 2;
            let mut v;
            loop {
                let mm = t[(im, im)];
                let r = shift[0] - mm;
                let s = shift[1] - mm;
                v = [
                    (r * s - shift[2]) / t[(im + 1, im)] + t[(im, im + 1)],
                    t[(im + 1, im + 1)] - mm - r - s,
                    t[(im + 2, im + 1)],
                ];
                if im == il {
                    break;
                }
                let lhs = t[(im, im - 1)] * (v[1].abs() + v[2].abs());
                let rhs = v[0] * (t[(im - 1, im - 1)].abs() + mm.abs() + t[(im + 1, im + 1)].abs());
                if lhs.abs() < f32::EPSILON * rhs {
                    break;
                }
                im -= 1;
            }
            for k in im..=iu - 2 {
                let vv = if k == im {
                    v.to_vec()
                } else {
                    (k..k + 3).map(|i| t[(i, k - 1)]).collect()
                };
                let (e, tau, beta) = reflector(&vv);
                if beta != 0. {
                    if k == im && k > il {
                        t[(k, k - 1)] = -t[(k, k - 1)];
                    } else if k != im {
                        t[(k, k - 1)] = beta;
                    }
                    left(&mut t, k, k, 3, n - k, &e, tau);
                    right(&mut t, 0, k, iu.min(k + 3) + 1, 3, &e, tau);
                    right(&mut u, 0, k, n, 3, &e, tau);
                }
            }
            let vv = [t[(iu - 1, iu - 2)], t[(iu, iu - 2)]];
            let (e, tau, beta) = reflector(&vv);
            if beta != 0. {
                t[(iu - 1, iu - 2)] = beta;
                left(&mut t, iu - 1, iu - 1, 2, n - iu + 1, &e, tau);
                right(&mut t, 0, iu - 1, iu + 1, 2, &e, tau);
                right(&mut u, 0, iu - 1, n, 2, &e, tau);
            }
            for i in im + 2..=iu {
                t[(i, i - 2)] = 0.;
                if i > im + 2 {
                    t[(i, i - 3)] = 0.;
                }
            }
        }
    }
    t *= scale;
    (t, u)
}

/// EigenSolver's real-Schur eigenvalues and normalized complex eigenvectors.
pub(super) fn eigen(a: &Mat) -> (Vec<C>, DMatrix<C>) {
    let n = a.nrows();
    let (mut t, mut u) = schur(a);
    let mut vals = vec![C::new(0., 0.); n];
    let mut i = 0;
    while i < n {
        if i + 1 == n || t[(i + 1, i)] == 0. {
            vals[i] = C::new(t[(i, i)], 0.);
            i += 1;
        } else {
            let p = 0.5 * (t[(i, i)] - t[(i + 1, i + 1)]);
            let max = p.abs().max(t[(i + 1, i)].abs().max(t[(i, i + 1)].abs()));
            let t0 = t[(i + 1, i)] / max;
            let t1 = t[(i, i + 1)] / max;
            let pp = p / max;
            let z = max * (pp * pp + t0 * t1).abs().sqrt();
            vals[i] = C::new(t[(i + 1, i + 1)] + p, z);
            vals[i + 1] = C::new(t[(i + 1, i + 1)] + p, -z);
            i += 2;
        }
    }
    let mut norm = 0.;
    for j in 0..n {
        let mut sum = 0.;
        for k in j.saturating_sub(1)..n {
            sum += t[(j, k)].abs();
        }
        norm += sum;
    }
    if norm != 0. {
        let mut nn = n as isize - 1;
        while nn >= 0 {
            let k = nn as usize;
            let p = vals[k].re;
            let q = vals[k].im;
            if q == 0. {
                let (mut lastr, mut lastw) = (0., 0.);
                let mut l = k;
                t[(k, k)] = 1.;
                for i in (0..k).rev() {
                    let w = t[(i, i)] - p;
                    let mut r = 0.;
                    for j in l..=k {
                        r += t[(i, j)] * t[(j, k)];
                    }
                    if vals[i].im < 0. {
                        lastw = w;
                        lastr = r;
                    } else {
                        l = i;
                        if vals[i].im == 0. {
                            t[(i, k)] = -r / if w != 0. { w } else { f32::EPSILON * norm };
                        } else {
                            let x = t[(i, i + 1)];
                            let y = t[(i + 1, i)];
                            let d = (vals[i].re - p) * (vals[i].re - p) + vals[i].im * vals[i].im;
                            let z = (x * lastr - lastw * r) / d;
                            t[(i, k)] = z;
                            t[(i + 1, k)] = if x.abs() > lastw.abs() {
                                (-r - w * z) / x
                            } else {
                                (-lastr - y * z) / lastw
                            };
                        }
                        let z = t[(i, k)].abs();
                        if (f32::EPSILON * z) * z > 1. {
                            for j in i..n {
                                t[(j, k)] /= z;
                            }
                        }
                    }
                }
            } else if q < 0. && k > 0 {
                let (mut lastra, mut lastsa, mut lastw) = (0., 0., 0.);
                let mut l = k - 1;
                if t[(k, k - 1)].abs() > t[(k - 1, k)].abs() {
                    t[(k - 1, k - 1)] = q / t[(k, k - 1)];
                    t[(k - 1, k)] = -(t[(k, k)] - p) / t[(k, k - 1)];
                } else {
                    let z = source_complex_div(
                        C::new(0., -t[(k - 1, k)]),
                        C::new(t[(k - 1, k - 1)] - p, q),
                    );
                    t[(k - 1, k - 1)] = z.re;
                    t[(k - 1, k)] = z.im;
                }
                t[(k, k - 1)] = 0.;
                t[(k, k)] = 1.;
                for i in (0..k - 1).rev() {
                    let (mut ra, mut sa) = (0., 0.);
                    for j in l..=k {
                        ra += t[(i, j)] * t[(j, k - 1)];
                        sa += t[(i, j)] * t[(j, k)];
                    }
                    let w = t[(i, i)] - p;
                    if vals[i].im < 0. {
                        lastw = w;
                        lastra = ra;
                        lastsa = sa;
                    } else {
                        l = i;
                        if vals[i].im == 0. {
                            let z = source_complex_div(C::new(-ra, -sa), C::new(w, q));
                            t[(i, k - 1)] = z.re;
                            t[(i, k)] = z.im;
                        } else {
                            let x = t[(i, i + 1)];
                            let y = t[(i + 1, i)];
                            let mut vr = (vals[i].re - p) * (vals[i].re - p)
                                + vals[i].im * vals[i].im
                                - q * q;
                            let vi = ((vals[i].re - p) * 2.) * q;
                            if vr == 0. && vi == 0. {
                                vr = (f32::EPSILON * norm)
                                    * (w.abs() + q.abs() + x.abs() + y.abs() + lastw.abs());
                            }
                            let z = source_complex_div(
                                C::new(
                                    x * lastra - lastw * ra + q * sa,
                                    x * lastsa - lastw * sa - q * ra,
                                ),
                                C::new(vr, vi),
                            );
                            t[(i, k - 1)] = z.re;
                            t[(i, k)] = z.im;
                            if x.abs() > lastw.abs() + q.abs() {
                                t[(i + 1, k - 1)] = (-ra - w * t[(i, k - 1)] + q * t[(i, k)]) / x;
                                t[(i + 1, k)] = (-sa - w * t[(i, k)] - q * t[(i, k - 1)]) / x;
                            } else {
                                let z = source_complex_div(
                                    C::new(-lastra - y * t[(i, k - 1)], -lastsa - y * t[(i, k)]),
                                    C::new(lastw, q),
                                );
                                t[(i + 1, k - 1)] = z.re;
                                t[(i + 1, k)] = z.im;
                            }
                        }
                        let z = t[(i, k - 1)].abs().max(t[(i, k)].abs());
                        if (f32::EPSILON * z) * z > 1. {
                            for j in i..n {
                                t[(j, k - 1)] /= z;
                                t[(j, k)] /= z;
                            }
                        }
                    }
                }
                nn -= 1;
            }
            nn -= 1;
        }
        for j in (0..n).rev() {
            let coefficients = (0..=j).map(|k| t[(k, j)]).collect::<Vec<_>>();
            let v = column_product(&u, 0, 0, n, j + 1, &coefficients);
            for i in 0..n {
                u[(i, j)] = v[i];
            }
        }
    }
    let mut vectors = DMatrix::zeros(n, n);
    let mut j = 0;
    while j < n {
        if vals[j].im.abs() <= vals[j].re.abs() * (2. * f32::EPSILON) || j + 1 == n {
            for i in 0..n {
                vectors[(i, j)] = C::new(u[(i, j)], 0.);
            }
            normalize(&mut vectors, j);
        } else {
            for i in 0..n {
                vectors[(i, j)] = C::new(u[(i, j)], u[(i, j + 1)]);
                vectors[(i, j + 1)] = C::new(u[(i, j)], -u[(i, j + 1)]);
            }
            normalize(&mut vectors, j);
            normalize(&mut vectors, j + 1);
            j += 1;
        }
        j += 1;
    }
    (vals, vectors)
}
fn normalize(v: &mut DMatrix<C>, col: usize) {
    // complex abs2 has no packet implementation in Eigen3.3.8; sum scalar.
    let n = v.nrows();
    let mut sum = 0.;
    for i in 0..n {
        let z = v[(i, col)];
        sum += z.re * z.re + z.im * z.im;
    }
    if sum > 0. {
        let norm = C::new(sum.sqrt(), 0.);
        let start = ((2 - (col * n) % 2) % 2).min(n);
        let end = start + (n - start) / 2 * 2;
        for i in 0..start {
            v[(i, col)] = source_complex_div(v[(i, col)], norm);
        }
        for i in start..end {
            v[(i, col)] = source_complex_packet_div(v[(i, col)], norm);
        }
        for i in end..n {
            v[(i, col)] = source_complex_div(v[(i, col)], norm);
        }
    }
}
