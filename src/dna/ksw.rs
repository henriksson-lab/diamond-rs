//! Scalar KSW2 extension recurrence, band traceback, flags, and Z-drop.
//!
//! Translated from lh3/ksw2 revision
//! `289609bd9e5381a13b16239d0a7703f1ff03f9ca`, `ksw2_extz.c` and `ksw2.h`.
//! DIAMOND's checked-in wrapper selects this scalar engine on macOS; its
//! Linux SSE engine has separate approximate-max/end-bonus behavior. The
//! adjacent checkout omits the dependency, but recovered Source tree
//! `69704a761fc8160e1dfb1e191383e7565d8834c9` (last present at commit
//! `08b3cbc5dc35ca1d3ab0a33eb5ee331b9a468df5`) contains byte-identical scalar
//! source/header. Its GNU SSE source differs only by MSVC feature defines;
//! the pinned DNA caller tree is unchanged since that historical snapshot.
//!
//! MIT License
//! Copyright (c) 2018- Dana-Farber Cancer Institute
//! Copyright (c) 2017-2018 Broad Institute, Inc.
//!
//! Permission is hereby granted, free of charge, to any person obtaining a
//! copy of this software and associated documentation files (the "Software"),
//! to deal in the Software without restriction, including without limitation
//! the rights to use, copy, modify, merge, publish, distribute, sublicense,
//! and/or sell copies of the Software, and to permit persons to whom the
//! Software is furnished to do so, subject to the following conditions:
//!
//! The above copyright notice and this permission notice shall be included
//! in all copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

pub const NEG_INF: i32 = -0x4000_0000;
pub const SCORE_ONLY: u32 = 0x01;
pub const RIGHT: u32 = 0x02;
pub const GENERIC_SC: u32 = 0x04;
pub const APPROX_MAX: u32 = 0x08;
pub const APPROX_DROP: u32 = 0x10;
pub const EXTZ_ONLY: u32 = 0x40;
pub const REV_CIGAR: u32 = 0x80;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Extz {
    pub max: i32,
    pub zdropped: bool,
    pub max_q: i32,
    pub max_t: i32,
    pub mqe: i32,
    pub mqe_t: i32,
    pub mte: i32,
    pub mte_q: i32,
    pub score: i32,
    pub reach_end: bool,
    /// Run-length operations in Source's MID encoding.
    pub cigar: Vec<(i32, char)>,
}

impl Default for Extz {
    fn default() -> Self {
        Self {
            max: 0,
            zdropped: false,
            max_q: -1,
            max_t: -1,
            mqe: NEG_INF,
            mqe_t: -1,
            mte: NEG_INF,
            mte_q: -1,
            score: NEG_INF,
            reach_end: false,
            cigar: Vec::new(),
        }
    }
}

pub(super) fn push(cigar: &mut Vec<(i32, char)>, op: char, length: i32) {
    if let Some((last, previous)) = cigar.last_mut() {
        if *previous == op {
            *last += length;
            return;
        }
    }
    cigar.push((length, op));
}

pub(super) fn filled<T: Clone>(length: usize, value: T) -> Result<Vec<T>, String> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|error| error.to_string())?;
    values.resize(length, value);
    Ok(values)
}

fn traceback(
    trace: &[u8],
    offsets: &[usize],
    ends: &[usize],
    columns: usize,
    mut i: i32,
    mut j: i32,
    reverse: bool,
) -> Result<Vec<(i32, char)>, String> {
    let mut cigar = Vec::new();
    let mut state = 0u8;
    while i >= 0 && j >= 0 {
        let offset = offsets[i as usize];
        if j as usize > ends[i as usize] {
            return Err("KSW traceback endpoint outside band".into());
        }
        let force = ((j as usize) < offset).then_some(2u8);
        let direction = if force.is_none() {
            let column = j as usize - offset;
            if column >= columns {
                return Err("KSW traceback endpoint outside band".into());
            }
            *trace
                .get(i as usize * columns + column)
                .ok_or("KSW traceback endpoint outside matrix")?
        } else {
            0
        };
        if state == 0 {
            state = direction & 7;
        } else if direction >> (state + 2) & 1 == 0 {
            state = 0;
        }
        if state == 0 {
            state = direction & 7;
        }
        if let Some(force) = force {
            state = force;
        }
        match state {
            0 => {
                push(&mut cigar, 'M', 1);
                i -= 1;
                j -= 1;
            }
            1 => {
                push(&mut cigar, 'D', 1);
                i -= 1;
            }
            2 => {
                push(&mut cigar, 'I', 1);
                j -= 1;
            }
            _ => return Err("Invalid KSW traceback state".into()),
        }
    }
    if i >= 0 {
        push(&mut cigar, 'D', i + 1);
    }
    if j >= 0 {
        push(&mut cigar, 'I', j + 1);
    }
    if !reverse {
        cigar.reverse();
    }
    Ok(cigar)
}

pub(super) fn apply_zdrop(
    ez: &mut Extz,
    score: i32,
    target: i32,
    query: i32,
    zdrop: i32,
    extension: i32,
) -> bool {
    if score > ez.max {
        ez.max = score;
        ez.max_t = target;
        ez.max_q = query;
    } else if target >= ez.max_t && query >= ez.max_q {
        let diagonal_distance = ((target - ez.max_t) - (query - ez.max_q)).abs();
        if zdrop >= 0 && ez.max - score > zdrop + diagonal_distance * extension {
            ez.zdropped = true;
            return true;
        }
    }
    false
}

/// KSW2 `ksw_extz` with explicit ownership and checked dimensions. Approximate
/// flags are accepted and ignored exactly as by the upstream scalar engine.
/// Empty inputs return reset state (the SIMD engine's defined behavior);
/// scalar C's invalid accesses on exhausted bands are reported as errors.
#[allow(clippy::too_many_arguments)]
pub fn extz(
    query: &[u8],
    target: &[u8],
    alphabet: usize,
    matrix: &[i8],
    gap_open: i8,
    gap_extend: i8,
    band: i32,
    zdrop: i32,
    flags: u32,
) -> Result<Extz, String> {
    let mut ez = Extz::default();
    if query.is_empty() || target.is_empty() {
        return Ok(ez);
    }
    if alphabet == 0
        || matrix.len()
            < alphabet
                .checked_mul(alphabet)
                .ok_or("KSW matrix size overflow")?
        || query
            .iter()
            .chain(target)
            .any(|&letter| letter as usize >= alphabet)
    {
        return Err("Invalid KSW scoring alphabet".into());
    }
    let qlen = i32::try_from(query.len()).map_err(|_| "KSW query exceeds Loc range")?;
    let tlen = i32::try_from(target.len()).map_err(|_| "KSW target exceeds Loc range")?;
    let width = if band < 0 { qlen.max(tlen) } else { band };
    let columns = query
        .len()
        .min((width as usize).saturating_mul(2).saturating_add(1));
    let with_cigar = flags & SCORE_ONLY == 0;
    let right = flags & RIGHT != 0;
    let mut trace = Vec::new();
    if with_cigar {
        let cells = columns
            .checked_mul(target.len())
            .ok_or("KSW traceback size overflow")?;
        trace
            .try_reserve_exact(cells)
            .map_err(|error| error.to_string())?;
        trace.resize(cells, 0);
    }
    let mut offsets = filled(if with_cigar { target.len() } else { 0 }, 0usize)?;
    let mut ends = filled(if with_cigar { target.len() } else { 0 }, 0usize)?;
    let gap = i32::from(gap_open) + i32::from(gap_extend);
    let extension = i32::from(gap_extend);
    let mut h = filled(query.len() + 1, 0i32)?;
    let mut e = filled(query.len() + 1, 0i32)?;
    e[0] = -gap - gap;
    for j in 1..=query.len() {
        if j <= width as usize {
            h[j] = -(gap + extension * (j as i32 - 1));
            e[j] = -(gap + gap + extension * j as i32);
        } else {
            h[j] = NEG_INF;
            e[j] = NEG_INF;
        }
    }
    let mut max_j = 0;
    for i in 0..tlen {
        let start = (i - width).max(0) as usize;
        let end = (i.saturating_add(width)).min(qlen - 1) as usize;
        if start > end {
            return Err("KSW extension exhausted its band".into());
        }
        let mut previous = if start > 0 {
            NEG_INF
        } else {
            -(gap + extension * i)
        };
        let mut f = if start > 0 {
            NEG_INF
        } else {
            -(gap + gap + extension * i)
        };
        let mut row_max = NEG_INF;
        if with_cigar {
            offsets[i as usize] = start;
            ends[i as usize] = end;
        }
        let profile = target[i as usize] as usize * alphabet;
        for j in start..=end {
            let diagonal = h[j] + i32::from(matrix[profile + query[j] as usize]);
            h[j] = previous;
            let prefer_diagonal = if !with_cigar || !right {
                diagonal >= e[j]
            } else {
                diagonal > e[j]
            };
            let mut direction = if prefer_diagonal { 0 } else { 1 };
            let best = diagonal.max(e[j]);
            let prefer_best = if !with_cigar || !right {
                best >= f
            } else {
                best > f
            };
            if !prefer_best {
                direction = 2;
            }
            let score = best.max(f);
            previous = score;
            // Scalar SCORE_ONLY uses the left branch's row-maximum ties.
            if if !with_cigar || !right {
                row_max <= score
            } else {
                row_max < score
            } {
                row_max = score;
                max_j = j as i32;
            }
            let open = score - gap;
            let extend_e = e[j] - extension;
            let keep_e = if right && with_cigar {
                extend_e >= open
            } else {
                extend_e > open
            };
            if keep_e {
                direction |= 0x08;
            }
            e[j] = extend_e.max(open);
            let extend_f = f - extension;
            let keep_f = if right && with_cigar {
                extend_f >= open
            } else {
                extend_f > open
            };
            if keep_f {
                direction |= 0x10;
            }
            f = extend_f.max(open);
            if with_cigar {
                trace[i as usize * columns + j - start] = direction;
            }
        }
        h[end + 1] = previous;
        e[end + 1] = NEG_INF;
        if end == query.len() - 1 && h[query.len()] > ez.mqe {
            ez.mqe = h[query.len()];
            ez.mqe_t = i;
        }
        if i == tlen - 1 {
            ez.mte = row_max;
            ez.mte_q = max_j;
        }
        if apply_zdrop(&mut ez, row_max, i, max_j, zdrop, extension) {
            break;
        }
        if i == tlen - 1 && end == query.len() - 1 {
            ez.score = h[query.len()];
        }
    }
    if with_cigar {
        if !ez.zdropped && flags & EXTZ_ONLY == 0 {
            ez.cigar = traceback(
                &trace,
                &offsets,
                &ends,
                columns,
                tlen - 1,
                qlen - 1,
                flags & REV_CIGAR != 0,
            )?;
        } else if ez.max_t >= 0 && ez.max_q >= 0 {
            ez.cigar = traceback(
                &trace,
                &offsets,
                &ends,
                columns,
                ez.max_t,
                ez.max_q,
                flags & REV_CIGAR != 0,
            )?;
        }
    }
    Ok(ez)
}
