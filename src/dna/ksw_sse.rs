//! Portable interpretation of KSW2's actual SSE byte-lane recurrence.
//!
//! Source: lh3/ksw2 `289609bd9e5381a13b16239d0a7703f1ff03f9ca`,
//! `ksw2_extz2_sse.c`. Copyright (c) 2018- Dana-Farber Cancer Institute;
//! Copyright (c) 2017-2018 Broad Institute, Inc. Distributed under the MIT
//! license reproduced in the sibling `ksw.rs` translation.
//!
//! This preserves the GNU engine's anti-diagonal traversal, wrapping byte
//! operations, signed comparisons, unsigned clipping, approximate maxima,
//! padded lane stores, end bonus and rotated traceback. It deliberately does
//! not replace them with an affine-gap recurrence using unbounded integers.

use super::ksw::{
    apply_zdrop, filled, push, Extz, APPROX_DROP, APPROX_MAX, EXTZ_ONLY, GENERIC_SC, NEG_INF,
    REV_CIGAR, RIGHT, SCORE_ONLY,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flavor {
    Sse2,
    Sse41,
}

fn signed_max(a: u8, b: u8) -> u8 {
    if (a as i8) > (b as i8) {
        a
    } else {
        b
    }
}

fn scalar_carry(value: u8) -> [u8; 16] {
    // cvtsi32_si128 receives an int8_t promoted to a signed int. This is
    // observable in nonstandard byte ranges, unlike a single-byte insertion.
    let mut carry = [0; 16];
    carry[..4].copy_from_slice(&(i32::from(value as i8)).to_le_bytes());
    carry
}

fn backtrack(
    trace: &[u8],
    starts: &[usize],
    ends: &[usize],
    columns: usize,
    mut target: i32,
    mut query: i32,
    reverse: bool,
) -> Result<Vec<(i32, char)>, String> {
    let mut cigar = Vec::new();
    let mut state = 0u8;
    while target >= 0 && query >= 0 {
        let row = (target + query) as usize;
        let start = *starts
            .get(row)
            .ok_or("KSW rotated traceback row outside matrix")?;
        let end = ends[row];
        let forced = if (target as usize) < start {
            Some(2)
        } else if target as usize > end {
            Some(1)
        } else {
            None
        };
        let direction = if forced.is_some() {
            0
        } else {
            let column = target as usize - start;
            if column >= columns {
                return Err("KSW rotated traceback column outside matrix".into());
            }
            *trace
                .get(row * columns + column)
                .ok_or("KSW rotated traceback outside matrix")?
        };
        if state == 0 {
            state = direction & 7;
        } else if direction >> (state + 2) & 1 == 0 {
            state = 0;
        }
        if state == 0 {
            state = direction & 7;
        }
        if let Some(forced) = forced {
            state = forced;
        }
        match state {
            0 => {
                push(&mut cigar, 'M', 1);
                target -= 1;
                query -= 1;
            }
            1 => {
                push(&mut cigar, 'D', 1);
                target -= 1;
            }
            2 => {
                push(&mut cigar, 'I', 1);
                query -= 1;
            }
            _ => return Err("Invalid KSW rotated traceback state".into()),
        }
    }
    if target >= 0 {
        push(&mut cigar, 'D', target + 1);
    }
    if query >= 0 {
        push(&mut cigar, 'I', query + 1);
    }
    if !reverse {
        cigar.reverse();
    }
    Ok(cigar)
}

#[allow(clippy::too_many_arguments)]
pub fn extz2(
    query: &[u8],
    target: &[u8],
    alphabet: usize,
    matrix: &[i8],
    gap_open: i8,
    gap_extend: i8,
    band: i32,
    zdrop: i32,
    end_bonus: i32,
    flags: u32,
    flavor: Flavor,
) -> Result<Extz, String> {
    let mut ez = Extz::default();
    if alphabet == 0 || query.is_empty() || target.is_empty() {
        return Ok(ez);
    }
    if alphabet == 1 {
        return Err("KSW SSE scoring alphabet requires at least two symbols".into());
    }
    let matrix_size = alphabet
        .checked_mul(alphabet)
        .ok_or("KSW scoring matrix size overflow")?;
    if matrix.len() < matrix_size
        || query
            .iter()
            .chain(target)
            .any(|&base| base as usize >= alphabet)
    {
        return Err("Invalid KSW scoring alphabet".into());
    }
    let qlen = i32::try_from(query.len()).map_err(|_| "KSW query exceeds Loc range")?;
    let tlen = i32::try_from(target.len()).map_err(|_| "KSW target exceeds Loc range")?;
    let rows = qlen
        .checked_add(tlen)
        .and_then(|n| n.checked_sub(1))
        .ok_or("KSW diagonal count overflow")? as usize;
    let q = i32::from(gap_open);
    let e = i32::from(gap_extend);
    let qe = q + e;
    let mut min_score = i32::from(matrix[1.min(matrix_size - 1)]);
    for &score in &matrix[1..matrix_size] {
        min_score = min_score.min(i32::from(score));
    }
    if -min_score > 2 * qe {
        return Ok(ez);
    }
    let width = if band < 0 { qlen.max(tlen) } else { band };
    let lanes = target.len().div_ceil(16) * 16;
    let query_lanes = query.len().div_ceil(16) * 16;
    let columns = query
        .len()
        .min(target.len())
        .min(width as usize + 1)
        .div_ceil(16)
        .checked_add(1)
        .and_then(|n| n.checked_mul(16))
        .ok_or("KSW trace stride overflow")?;
    let memory_size = lanes
        .checked_mul(6)
        .and_then(|n| n.checked_add(query_lanes))
        .and_then(|n| n.checked_add(16))
        .ok_or("KSW lane allocation overflow")?;
    let mut memory = filled(memory_size, 0u8)?;
    let (u, v, x, y, s, sf, qr) = (
        0,
        lanes,
        lanes * 2,
        lanes * 3,
        lanes * 4,
        lanes * 5,
        lanes * 6,
    );
    memory[sf..sf + target.len()].copy_from_slice(target);
    for (at, &base) in query.iter().rev().enumerate() {
        memory[qr + at] = base;
    }
    let with_cigar = flags & SCORE_ONLY == 0;
    let approximate = flags & APPROX_MAX != 0;
    let right = flags & RIGHT != 0;
    let mut h = filled(if approximate { 0 } else { lanes }, NEG_INF)?;
    let mut trace = if with_cigar {
        filled(
            rows.checked_mul(columns)
                .ok_or("KSW trace allocation overflow")?,
            0u8,
        )?
    } else {
        Vec::new()
    };
    let mut starts = filled(if with_cigar { rows } else { 0 }, 0usize)?;
    let mut ends = filled(if with_cigar { rows } else { 0 }, 0usize)?;
    let qe2 = (qe * 2) as u8;
    let maximum = (i32::from(matrix[0]) + qe * 2) as u8;
    let wildcard = (alphabet - 1) as u8;
    let wildcard_score = if matrix[matrix_size - 1] == 0 {
        (-e) as u8
    } else {
        matrix[matrix_size - 1] as u8
    };
    let mut last_start = -1i32;
    let mut last_end = -1i32;
    let mut h0 = 0i32;
    let mut last_h0_target = 0i32;
    for r in 0..rows as i32 {
        let start0 = 0.max(r - qlen + 1).max((r - width + 1) >> 1);
        let end0 = (tlen - 1).min(r).min((r + width) >> 1);
        if start0 > end0 {
            ez.zdropped = true;
            break;
        }
        let start = (start0 as usize / 16) * 16;
        let end = (end0 as usize + 16) / 16 * 16 - 1;
        let (first_x, first_v) = if start > 0 {
            if start as i32 - 1 >= last_start && start as i32 - 1 <= last_end {
                (memory[x + start - 1], memory[v + start - 1])
            } else {
                (0, 0)
            }
        } else {
            (0, if r != 0 { gap_open as u8 } else { 0 })
        };
        if end >= r as usize {
            memory[y + r as usize] = 0;
            memory[u + r as usize] = if r != 0 { gap_open as u8 } else { 0 };
        }
        // Preserve Source's unaligned 16-byte score stores, including the
        // padding overlap into sf. The allocation layout matches its calloc.
        if flags & GENERIC_SC == 0 {
            for t in (start0 as usize..=end0 as usize).step_by(16) {
                let mut scores = [0; 16];
                for lane in 0..16 {
                    let target_base = memory[sf + t + lane];
                    let query_at = qr as i64 + i64::from(qlen - 1 - r) + (t + lane) as i64;
                    let query_base = *memory
                        .get(
                            usize::try_from(query_at)
                                .map_err(|_| "KSW padded query offset underflow")?,
                        )
                        .ok_or("KSW padded query offset overflow")?;
                    scores[lane] = if target_base == wildcard || query_base == wildcard {
                        wildcard_score
                    } else if target_base == query_base {
                        matrix[0] as u8
                    } else {
                        matrix[1] as u8
                    };
                }
                memory[s + t..s + t + 16].copy_from_slice(&scores);
            }
        } else {
            for t in start0 as usize..=end0 as usize {
                let query_at = (qr as i64 + i64::from(qlen - 1 - r) + t as i64) as usize;
                memory[s + t] =
                    matrix[memory[sf + t] as usize * alphabet + memory[query_at] as usize] as u8;
            }
        }
        if with_cigar {
            starts[r as usize] = start;
            ends[r as usize] = end;
        }
        let mut carry_x = scalar_carry(first_x);
        let mut carry_v = scalar_carry(first_v);
        for block in (start..=end).step_by(16) {
            let old_x: [u8; 16] = memory[x + block..x + block + 16].try_into().unwrap();
            let old_v: [u8; 16] = memory[v + block..v + block + 16].try_into().unwrap();
            let old_u: [u8; 16] = memory[u + block..u + block + 16].try_into().unwrap();
            let old_y: [u8; 16] = memory[y + block..y + block + 16].try_into().unwrap();
            let scores: [u8; 16] = memory[s + block..s + block + 16].try_into().unwrap();
            for lane in 0..16 {
                let previous_x = if lane == 0 {
                    carry_x[0]
                } else {
                    old_x[lane - 1] | carry_x[lane]
                };
                let previous_v = if lane == 0 {
                    carry_v[0]
                } else {
                    old_v[lane - 1] | carry_v[lane]
                };
                let mut z = scores[lane].wrapping_add(qe2);
                let a = previous_x.wrapping_add(previous_v);
                let b = old_y[lane].wrapping_add(old_u[lane]);
                let mut direction = if if right {
                    (a as i8) >= (z as i8)
                } else {
                    (a as i8) > (z as i8)
                } {
                    1u8
                } else {
                    0
                };
                z = if flavor == Flavor::Sse41 {
                    signed_max(z, a)
                } else {
                    (if (z as i8) > 0 { z } else { 0 }).max(a)
                };
                if if right {
                    (b as i8) >= (z as i8)
                } else {
                    (b as i8) > (z as i8)
                } {
                    direction = 2;
                }
                z = z.max(b).min(maximum);
                memory[u + block + lane] = z.wrapping_sub(previous_v);
                memory[v + block + lane] = z.wrapping_sub(old_u[lane]);
                z = z.wrapping_sub(gap_open as u8);
                let a = a.wrapping_sub(z);
                let b = b.wrapping_sub(z);
                let keep_a = if right && with_cigar {
                    (a as i8) >= 0
                } else {
                    (a as i8) > 0
                };
                let keep_b = if right && with_cigar {
                    (b as i8) >= 0
                } else {
                    (b as i8) > 0
                };
                memory[x + block + lane] = if keep_a { a } else { 0 };
                memory[y + block + lane] = if keep_b { b } else { 0 };
                if with_cigar {
                    if keep_a {
                        direction |= 8;
                    }
                    if keep_b {
                        direction |= 16;
                    }
                    trace[r as usize * columns + block - start + lane] = direction;
                }
            }
            carry_x = [0; 16];
            carry_x[0] = old_x[15];
            carry_v = [0; 16];
            carry_v[0] = old_v[15];
        }
        if !approximate {
            let (mut max_h, mut max_t);
            if r > 0 {
                let last = end0 as usize;
                h[last] = if last > 0 {
                    h[last - 1] + i32::from(memory[u + last]) - qe
                } else {
                    h[last] + i32::from(memory[v + last]) - qe
                };
                max_h = h[last];
                max_t = end0;
                let end1 = start0 + (end0 - start0) / 4 * 4;
                let mut lane_h = [max_h; 4];
                let mut lane_t = [max_t; 4];
                let mut t = start0;
                while t < end1 {
                    for lane in 0..4 {
                        let at = t as usize + lane;
                        h[at] = h[at]
                            .wrapping_add(i32::from(memory[v + at]))
                            .wrapping_sub(qe);
                        if h[at] > lane_h[lane] {
                            lane_h[lane] = h[at];
                            lane_t[lane] = t;
                        }
                    }
                    t += 4;
                }
                for lane in 0..4 {
                    if max_h < lane_h[lane] {
                        max_h = lane_h[lane];
                        max_t = lane_t[lane] + lane as i32;
                    }
                }
                while t < end0 {
                    h[t as usize] += i32::from(memory[v + t as usize]) - qe;
                    if h[t as usize] > max_h {
                        max_h = h[t as usize];
                        max_t = t;
                    }
                    t += 1;
                }
            } else {
                h[0] = i32::from(memory[v]) - qe - qe;
                max_h = h[0];
                max_t = 0;
            }
            if end0 == tlen - 1 && h[end0 as usize] > ez.mte {
                ez.mte = h[end0 as usize];
                // Upstream uses its rounded lane boundary here, not end0.
                ez.mte_q = r - end as i32;
            }
            if r - start0 == qlen - 1 && h[start0 as usize] > ez.mqe {
                ez.mqe = h[start0 as usize];
                ez.mqe_t = start0;
            }
            if apply_zdrop(&mut ez, max_h, max_t, r - max_t, zdrop, e) {
                break;
            }
            if r == rows as i32 - 1 && end0 == tlen - 1 {
                ez.score = h[target.len() - 1];
            }
        } else {
            if r > 0 {
                if last_h0_target >= start0
                    && last_h0_target <= end0
                    && last_h0_target + 1 >= start0
                    && last_h0_target + 1 <= end0
                {
                    let d0 = i32::from(memory[v + last_h0_target as usize]) - qe;
                    let d1 = i32::from(memory[u + last_h0_target as usize + 1]) - qe;
                    if d0 > d1 {
                        h0 += d0;
                    } else {
                        h0 += d1;
                        last_h0_target += 1;
                    }
                } else if last_h0_target >= start0 && last_h0_target <= end0 {
                    h0 += i32::from(memory[v + last_h0_target as usize]) - qe;
                } else {
                    last_h0_target += 1;
                    h0 += i32::from(memory[u + last_h0_target as usize]) - qe;
                }
                if flags & APPROX_DROP != 0
                    && apply_zdrop(&mut ez, h0, last_h0_target, r - last_h0_target, zdrop, e)
                {
                    break;
                }
            } else {
                h0 = i32::from(memory[v]) - qe - qe;
                last_h0_target = 0;
            }
            if r == rows as i32 - 1 && end0 == tlen - 1 {
                ez.score = h0;
            }
        }
        last_start = start as i32;
        last_end = end as i32;
    }
    if with_cigar {
        let reverse = flags & REV_CIGAR != 0;
        let endpoint = if !ez.zdropped && flags & EXTZ_ONLY == 0 {
            Some((tlen - 1, qlen - 1))
        } else if !ez.zdropped && flags & EXTZ_ONLY != 0 && ez.mqe + end_bonus > ez.max {
            ez.reach_end = true;
            Some((ez.mqe_t, qlen - 1))
        } else if ez.max_t >= 0 && ez.max_q >= 0 {
            Some((ez.max_t, ez.max_q))
        } else {
            None
        };
        if let Some((target, query)) = endpoint {
            ez.cigar = backtrack(&trace, &starts, &ends, columns, target, query, reverse)?;
        }
    }
    Ok(ez)
}
