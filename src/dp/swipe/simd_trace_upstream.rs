//! Structural AVX2 ports of upstream `dp/swipe/banded_swipe.h` traceback
//! specializations. Keep the moving-band coordinates, vector profile, rolling
//! matrix, and packed lane masks aligned with the C++ implementation.

use super::simd_trace::TraceTarget;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{Letter, LETTER_MASK, SEED_MASK};
use crate::dp::smith_waterman::SwResult;
use crate::stats::score_matrix::ScoreMatrix;
use std::arch::x86_64 as arch;

#[inline]
fn push_operation_run(operations: &mut Vec<(EditOperation, i32)>, op: EditOperation, count: i32) {
    if let Some((last_op, last_count)) = operations.last_mut() {
        if *last_op == op {
            *last_count += count;
            return;
        }
    }
    operations.push((op, count));
}

type V = arch::__m256i;

#[target_feature(enable = "avx2")]
unsafe fn standard_profile_i16(
    matrix: &[i8; 1024],
    subject: &[i8; 16],
    hard_masked: &[i8; 16],
) -> [V; 32] {
    let subject = arch::_mm_loadu_si128(subject.as_ptr().cast());
    let hard_masked = arch::_mm_loadu_si128(hard_masked.as_ptr().cast());
    let indices = arch::_mm_and_si128(subject, arch::_mm_set1_epi8(15));
    let high = arch::_mm_cmpgt_epi8(subject, arch::_mm_set1_epi8(15));
    let mut profile = [arch::_mm256_setzero_si256(); 32];
    for (query_letter, slot) in profile.iter_mut().enumerate() {
        let row = matrix.as_ptr().add(query_letter * 32);
        let low_scores = arch::_mm_loadu_si128(row.cast());
        let high_scores = arch::_mm_loadu_si128(row.add(16).cast());
        let low_scores = arch::_mm_shuffle_epi8(low_scores, indices);
        let high_scores = arch::_mm_shuffle_epi8(high_scores, indices);
        let scores = arch::_mm_blendv_epi8(low_scores, high_scores, high);
        let scores = arch::_mm_andnot_si128(hard_masked, scores);
        *slot = arch::_mm256_cvtepi8_epi16(scores);
    }
    profile
}

#[derive(Clone, Copy, Default)]
struct Trace8 {
    gap: u64,
    open: u64,
}

#[derive(Clone, Copy, Default)]
struct Trace16 {
    gap: u32,
    open: u32,
}

struct Trace16Matrix {
    packed: Vec<Trace16>,
    sparse: Vec<u8>,
    stride: usize,
}

impl Trace16Matrix {
    fn new<const SPARSE: bool>(cells: usize, lanes: usize) -> Self {
        if SPARSE {
            let stride = lanes.div_ceil(2);
            Self {
                packed: Vec::new(),
                sparse: vec![0; cells * stride],
                stride,
            }
        } else {
            Self {
                packed: vec![Trace16::default(); cells],
                sparse: Vec::new(),
                stride: 0,
            }
        }
    }

    #[inline]
    unsafe fn set<const SPARSE: bool>(
        &mut self,
        cell: usize,
        lanes: usize,
        gap_v: u32,
        gap_h: u32,
        open_v: u32,
        open_h: u32,
        active: u32,
    ) {
        if SPARSE {
            let base = cell * self.stride;
            for lane in 0..lanes {
                let bit = 1u32 << (2 * lane);
                let is_active = active & bit != 0;
                let vertical = is_active && gap_v & bit != 0;
                let horizontal = is_active && !vertical && gap_h & bit != 0;
                let state = if !is_active {
                    0
                } else if vertical {
                    2
                } else if horizontal {
                    3
                } else {
                    1
                };
                let nibble =
                    state | u8::from(open_v & bit != 0) << 2 | u8::from(open_h & bit != 0) << 3;
                let byte = self.sparse.get_unchecked_mut(base + lane / 2);
                if lane % 2 == 0 {
                    *byte = (*byte & 0xf0) | nibble;
                } else {
                    *byte = (*byte & 0x0f) | (nibble << 4);
                }
            }
        } else {
            const HMASK: u32 = 0x5555_5555;
            let active_low = active & HMASK;
            let vertical_low = (gap_v & HMASK) & active_low;
            let horizontal_low = (gap_h & HMASK) & active_low & !vertical_low;
            *self.packed.get_unchecked_mut(cell) = Trace16 {
                gap: (active_low & !vertical_low) | ((vertical_low | horizontal_low) << 1),
                open: (open_v & 0xAAAA_AAAA) | (open_h & HMASK),
            };
        }
    }

    #[inline]
    fn get<const SPARSE: bool>(&self, cell: usize, lane: usize) -> (u8, bool, bool) {
        if SPARSE {
            let byte = self.sparse[cell * self.stride + lane / 2];
            let nibble = if lane % 2 == 0 {
                byte & 0x0f
            } else {
                byte >> 4
            };
            (nibble & 3, nibble & 4 != 0, nibble & 8 != 0)
        } else {
            let low = 1u32 << (2 * lane);
            let high = 2u32 << (2 * lane);
            let state = match (
                self.packed[cell].gap & high != 0,
                self.packed[cell].gap & low != 0,
            ) {
                (false, false) => 0,
                (false, true) => 1,
                (true, false) => 2,
                (true, true) => 3,
            };
            (
                state,
                self.packed[cell].open & high != 0,
                self.packed[cell].open & low != 0,
            )
        }
    }
}

#[target_feature(enable = "avx2")]
pub(super) unsafe fn trace_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> (Vec<SwResult>, u32) {
    if query_cbs.is_empty() {
        trace_i8_impl::<false>(query, targets, matrix, query_cbs)
    } else {
        trace_i8_impl::<true>(query, targets, matrix, query_cbs)
    }
}

#[target_feature(enable = "avx2")]
unsafe fn trace_i8_impl<const HAS_CBS: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> (Vec<SwResult>, u32) {
    debug_assert_eq!(HAS_CBS, !query_cbs.is_empty());
    debug_assert!(!HAS_CBS || query_cbs.len() >= query.len());
    const LANES: usize = 32;
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    if band == 0 {
        return (vec![SwResult::default(); targets.len()], 0);
    }
    let i1 = targets
        .iter()
        .map(|target| (target.d_end - 1).max(0))
        .min()
        .unwrap_or(0);
    let i0 = i1 + 1 - band as i32;
    let mut subject_start = [0i32; LANES];
    let mut band_offset = [0usize; LANES];
    let mut columns = 0usize;
    for lane in 0..targets.len() {
        let target = &targets[lane];
        let expanded_begin = target.d_end - band as i32;
        subject_start[lane] = i1 - (target.d_end - 1);
        band_offset[lane] = (target.d_begin - expanded_begin).max(0) as usize;
        let end = ((query.len() as i32 - 1 - expanded_begin).min(target.subject.len() as i32 - 1)
            + 1)
        .max(0);
        columns = columns.max((end - subject_start[lane]).max(0) as usize);
    }
    let zero = arch::_mm256_set1_epi8(i8::MIN);
    let max_value = arch::_mm256_set1_epi8(i8::MAX);
    let mut score_row = vec![zero; band];
    let mut hgap_row = vec![zero; band + 1];
    let mut trace = vec![Trace8::default(); (columns + 1) * band];
    let mut gap_open = [0i8; LANES];
    let mut gap_extend = [0i8; LANES];
    let mut overflow_mask = 0u32;
    for lane in 0..targets.len() {
        let scale = targets[lane].matrix_scale.max(1);
        let go = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=63).contains(&go) || !(0..=63).contains(&ge) {
            overflow_mask |= 1 << lane;
        }
        gap_open[lane] = go.clamp(0, 63) as i8;
        gap_extend[lane] = ge.clamp(0, 63) as i8;
    }
    let go = arch::_mm256_loadu_si256(gap_open.as_ptr().cast());
    let ge = arch::_mm256_loadu_si256(gap_extend.as_ptr().cast());
    let mut standard_lanes = [0i8; LANES];
    for lane in 0..targets.len() {
        if targets[lane].matrix.is_none() {
            standard_lanes[lane] = -1;
        }
    }
    let standard_mask = arch::_mm256_loadu_si256(standard_lanes.as_ptr().cast());
    let row_masks: Vec<V> = (0..band)
        .map(|row| {
            let mut lanes = [0i8; LANES];
            for lane in 0..targets.len() {
                if row >= band_offset[lane] {
                    lanes[lane] = -1;
                }
            }
            arch::_mm256_loadu_si256(lanes.as_ptr().cast())
        })
        .collect();
    let mut best = [i8::MIN; LANES];
    let mut best_col = [0usize; LANES];
    let mut best_row = [0usize; LANES];

    for column in 0..columns {
        let mut subject_letter = [0usize; LANES];
        let mut active_lanes = [0i8; LANES];
        let mut hard_masked = [false; LANES];
        for lane in 0..targets.len() {
            let pos = subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                let letter = targets[lane].subject[pos as usize];
                subject_letter[lane] = (letter & LETTER_MASK) as usize;
                active_lanes[lane] = -1;
                hard_masked[lane] = letter & SEED_MASK != 0;
            }
        }
        let active = arch::_mm256_loadu_si256(active_lanes.as_ptr().cast());
        let mut profile = [zero; 32];
        for (query_letter, slot) in profile.iter_mut().enumerate() {
            let mut scores = [0i8; LANES];
            for lane in 0..targets.len() {
                if active_lanes[lane] == 0 || hard_masked[lane] {
                    continue;
                }
                scores[lane] = if let Some(adjusted) = targets[lane].matrix {
                    adjusted.scores[subject_letter[lane] * 32 + query_letter]
                } else {
                    matrix.matrix8()[query_letter * 32 + subject_letter[lane]]
                };
            }
            *slot = arch::_mm256_loadu_si256(scores.as_ptr().cast());
        }
        let moving_i0 = i0 + column as i32;
        let query_begin = moving_i0.max(0);
        let query_end = (i1 + column as i32).min(query.len() as i32 - 1) + 1;
        let mut vertical = zero;
        let mut col_best = zero;
        let mut row_counter = arch::_mm256_set1_epi8((query_begin - moving_i0) as i8);
        let mut row_max = zero;
        for q in query_begin..query_end {
            let row = (q - moving_i0) as usize;
            let cell_mask = arch::_mm256_and_si256(active, row_masks[row]);
            let query_bias = if HAS_CBS {
                *query_cbs.get_unchecked(q as usize)
            } else {
                0
            };
            let bias = arch::_mm256_and_si256(arch::_mm256_set1_epi8(query_bias), standard_mask);
            let substitution =
                arch::_mm256_adds_epi8(profile[(query[q as usize] & LETTER_MASK) as usize], bias);
            let diagonal = arch::_mm256_adds_epi8(score_row[row], substitution);
            let horizontal = hgap_row[row + 1];
            let vertical_before = vertical;
            let mut score = arch::_mm256_max_epi8(diagonal, horizontal);
            score = arch::_mm256_max_epi8(score, vertical_before);
            score = arch::_mm256_max_epi8(score, zero);
            score = arch::_mm256_blendv_epi8(zero, score, cell_mask);
            col_best = arch::_mm256_max_epi8(col_best, score);
            // Upstream VectorRowCounter records the last row equal to the
            // column maximum, not only a strictly better row.
            let at_column_max = arch::_mm256_cmpeq_epi8(col_best, score);
            row_max = arch::_mm256_blendv_epi8(row_max, row_counter, at_column_max);
            row_counter = arch::_mm256_add_epi8(row_counter, arch::_mm256_set1_epi8(1));
            let gap_v =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(score, vertical_before)) as u32;
            let gap_h =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(score, horizontal)) as u32;
            let active_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpgt_epi8(score, zero)) as u32;
            let open = arch::_mm256_subs_epi8(score, go);
            let next_horizontal =
                arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, ge), open);
            let next_vertical =
                arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical_before, ge), open);
            let open_v =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(next_vertical, open)) as u32;
            let open_h =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(next_horizontal, open)) as u32;
            // Encode four states in the existing V/H pair: 00 inactive,
            // 01 diagonal, 10 vertical, 11 horizontal. This retains exact
            // local-alignment termination without an extra active plane.
            let vertical_state = gap_v & active_bits;
            let horizontal_state = gap_h & active_bits & !vertical_state;
            let low = active_bits & !vertical_state;
            let high = vertical_state | horizontal_state;
            trace[(column + 1) * band + row] = Trace8 {
                gap: (u64::from(high) << 32) | u64::from(low),
                open: (u64::from(open_v) << 32) | u64::from(open_h),
            };
            score_row[row] = score;
            hgap_row[row] = arch::_mm256_blendv_epi8(zero, next_horizontal, cell_mask);
            vertical = arch::_mm256_blendv_epi8(zero, next_vertical, cell_mask);
            overflow_mask |=
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi8(score, max_value)) as u32;
        }
        let mut column_scores = [i8::MIN; LANES];
        let mut column_rows = [0i8; LANES];
        arch::_mm256_storeu_si256(column_scores.as_mut_ptr().cast(), col_best);
        arch::_mm256_storeu_si256(column_rows.as_mut_ptr().cast(), row_max);
        for lane in 0..targets.len() {
            if column_scores[lane] > best[lane] {
                best[lane] = column_scores[lane];
                best_col[lane] = column;
                best_row[lane] = column_rows[lane] as u8 as usize;
            }
        }
    }
    let scores: Vec<i32> = best[..targets.len()]
        .iter()
        .map(|&score| i32::from(score) - i32::from(i8::MIN))
        .collect();
    (
        finish_i8(
            query,
            targets,
            &trace,
            band,
            i0,
            &subject_start,
            &scores,
            &best_col,
            &best_row,
        ),
        overflow_mask,
    )
}

#[target_feature(enable = "avx2")]
pub(super) unsafe fn trace_i16(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> (Vec<SwResult>, u32) {
    match (targets.len() < 3, query_cbs.is_empty()) {
        (true, true) => trace_i16_impl::<true, false>(query, targets, matrix, query_cbs),
        (true, false) => trace_i16_impl::<true, true>(query, targets, matrix, query_cbs),
        (false, true) => trace_i16_impl::<false, false>(query, targets, matrix, query_cbs),
        (false, false) => trace_i16_impl::<false, true>(query, targets, matrix, query_cbs),
    }
}

#[target_feature(enable = "avx2")]
unsafe fn trace_i16_impl<const SPARSE_TRACE: bool, const HAS_CBS: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
) -> (Vec<SwResult>, u32) {
    debug_assert_eq!(HAS_CBS, !query_cbs.is_empty());
    debug_assert!(!HAS_CBS || query_cbs.len() >= query.len());
    const LANES: usize = 16;
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin).max(0) as usize)
        .max()
        .unwrap_or(0);
    if band == 0 {
        return (vec![SwResult::default(); targets.len()], 0);
    }
    let i1 = targets
        .iter()
        .map(|target| (target.d_end - 1).max(0))
        .min()
        .unwrap_or(0);
    let i0 = i1 + 1 - band as i32;
    let mut subject_start = [0i32; LANES];
    let mut band_offset = [0usize; LANES];
    let mut columns = 0usize;
    for lane in 0..targets.len() {
        let target = &targets[lane];
        let expanded_begin = target.d_end - band as i32;
        subject_start[lane] = i1 - (target.d_end - 1);
        band_offset[lane] = (target.d_begin - expanded_begin).max(0) as usize;
        let end = ((query.len() as i32 - 1 - expanded_begin).min(target.subject.len() as i32 - 1)
            + 1)
        .max(0);
        columns = columns.max((end - subject_start[lane]).max(0) as usize);
    }
    let zero = arch::_mm256_set1_epi16(i16::MIN);
    let mut score_row = vec![zero; band];
    let mut hgap_row = vec![zero; band + 1];
    let mut trace = Trace16Matrix::new::<SPARSE_TRACE>((columns + 1) * band, targets.len());
    let mut gap_open = [0i16; LANES];
    let mut gap_extend = [0i16; LANES];
    let mut overflow_mask = 0u32;
    for lane in 0..targets.len() {
        let scale = targets[lane].matrix_scale.max(1);
        let go = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=16_000).contains(&go) || !(0..=16_000).contains(&ge) {
            overflow_mask |= 1 << lane;
        }
        gap_open[lane] = go.clamp(0, 16_000) as i16;
        gap_extend[lane] = ge.clamp(0, 16_000) as i16;
    }
    let go = arch::_mm256_loadu_si256(gap_open.as_ptr().cast());
    let ge = arch::_mm256_loadu_si256(gap_extend.as_ptr().cast());
    let mut standard_lanes = [0i16; LANES];
    for lane in 0..targets.len() {
        if targets[lane].matrix.is_none() {
            standard_lanes[lane] = -1;
        }
    }
    let standard_mask = arch::_mm256_loadu_si256(standard_lanes.as_ptr().cast());
    let standard_only = targets.iter().all(|target| target.matrix.is_none());
    let row_masks: Vec<V> = (0..band)
        .map(|row| {
            let mut lanes = [0i16; LANES];
            for lane in 0..targets.len() {
                if row >= band_offset[lane] {
                    lanes[lane] = -1;
                }
            }
            arch::_mm256_loadu_si256(lanes.as_ptr().cast())
        })
        .collect();
    let mut best = [i16::MIN; LANES];
    let mut best_col = [0usize; LANES];
    let mut best_row = [0usize; LANES];

    for column in 0..columns {
        let mut subject_letter = [0usize; LANES];
        let mut subject_bytes = [0i8; LANES];
        let mut active_lanes = [0i16; LANES];
        let mut hard_masked = [false; LANES];
        let mut hard_mask_bytes = [0i8; LANES];
        for lane in 0..targets.len() {
            let pos = subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                let letter = targets[lane].subject[pos as usize];
                subject_letter[lane] = (letter & LETTER_MASK) as usize;
                subject_bytes[lane] = (letter & LETTER_MASK) as i8;
                active_lanes[lane] = -1;
                hard_masked[lane] = letter & SEED_MASK != 0;
                hard_mask_bytes[lane] = if hard_masked[lane] { -1 } else { 0 };
            }
        }
        let active = arch::_mm256_loadu_si256(active_lanes.as_ptr().cast());
        let profile = if standard_only {
            standard_profile_i16(matrix.matrix8(), &subject_bytes, &hard_mask_bytes)
        } else {
            let mut profile = [zero; 32];
            for (query_letter, slot) in profile.iter_mut().enumerate() {
                let mut scores = [0i16; LANES];
                for lane in 0..targets.len() {
                    if active_lanes[lane] == 0 || hard_masked[lane] {
                        continue;
                    }
                    scores[lane] = if let Some(adjusted) = targets[lane].matrix {
                        adjusted.scores[subject_letter[lane] * 32 + query_letter] as i16
                    } else {
                        matrix.matrix16()[query_letter * 32 + subject_letter[lane]]
                    };
                }
                *slot = arch::_mm256_loadu_si256(scores.as_ptr().cast());
            }
            profile
        };
        let moving_i0 = i0 + column as i32;
        let query_begin = moving_i0.max(0);
        let query_end = (i1 + column as i32).min(query.len() as i32 - 1) + 1;
        let mut vertical = zero;
        let mut col_best = zero;
        let mut row_counter = arch::_mm256_set1_epi16((query_begin - moving_i0) as i16);
        let mut row_max = zero;
        for q in query_begin..query_end {
            let row = (q - moving_i0) as usize;
            // The moving-band clamps q to the query and row to 0..band; the
            // profile index is masked to 0..31. These are the same pointer
            // invariants used by upstream's matrix iterator.
            let cell_mask = arch::_mm256_and_si256(active, *row_masks.get_unchecked(row));
            let q_index = q as usize;
            let query_bias = if HAS_CBS {
                *query_cbs.get_unchecked(q_index)
            } else {
                0
            };
            let bias =
                arch::_mm256_and_si256(arch::_mm256_set1_epi16(query_bias as i16), standard_mask);
            let query_letter = *query.get_unchecked(q_index) & LETTER_MASK;
            let substitution =
                arch::_mm256_adds_epi16(*profile.get_unchecked(query_letter as usize), bias);
            let diagonal = arch::_mm256_adds_epi16(*score_row.get_unchecked(row), substitution);
            let horizontal = *hgap_row.get_unchecked(row + 1);
            let vertical_before = vertical;
            let mut score = arch::_mm256_max_epi16(diagonal, horizontal);
            score = arch::_mm256_max_epi16(score, vertical_before);
            score = arch::_mm256_max_epi16(score, zero);
            score = arch::_mm256_blendv_epi8(zero, score, cell_mask);
            col_best = arch::_mm256_max_epi16(col_best, score);
            let at_column_max = arch::_mm256_cmpeq_epi16(col_best, score);
            row_max = arch::_mm256_blendv_epi8(row_max, row_counter, at_column_max);
            row_counter = arch::_mm256_add_epi16(row_counter, arch::_mm256_set1_epi16(1));
            let gap_v =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(score, vertical_before)) as u32;
            let gap_h =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(score, horizontal)) as u32;
            let active_bits =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpgt_epi16(score, zero)) as u32;
            let open = arch::_mm256_subs_epi16(score, go);
            let next_horizontal =
                arch::_mm256_max_epi16(arch::_mm256_subs_epi16(horizontal, ge), open);
            let next_vertical =
                arch::_mm256_max_epi16(arch::_mm256_subs_epi16(vertical_before, ge), open);
            let open_v =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(next_vertical, open)) as u32;
            let open_h =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(next_horizontal, open)) as u32;
            let trace_index = (column + 1) * band + row;
            trace.set::<SPARSE_TRACE>(
                trace_index,
                targets.len(),
                gap_v,
                gap_h,
                open_v,
                open_h,
                active_bits,
            );
            *score_row.get_unchecked_mut(row) = score;
            *hgap_row.get_unchecked_mut(row) =
                arch::_mm256_blendv_epi8(zero, next_horizontal, cell_mask);
            vertical = arch::_mm256_blendv_epi8(zero, next_vertical, cell_mask);
        }
        let mut column_scores = [i16::MIN; LANES];
        let mut column_rows = [0i16; LANES];
        arch::_mm256_storeu_si256(column_scores.as_mut_ptr().cast(), col_best);
        arch::_mm256_storeu_si256(column_rows.as_mut_ptr().cast(), row_max);
        for lane in 0..targets.len() {
            if column_scores[lane] > best[lane] {
                best[lane] = column_scores[lane];
                best_col[lane] = column;
                best_row[lane] = column_rows[lane] as u16 as usize;
            }
        }
    }
    for lane in 0..targets.len() {
        // `best` is a monotonic maximum, so a saturated cell remains visible
        // here.  Upstream likewise tests the final lane maximum instead of
        // materializing a saturation mask for every DP cell.
        if best[lane] == i16::MAX {
            overflow_mask |= 1 << lane;
        }
    }
    let scores: Vec<i32> = best[..targets.len()]
        .iter()
        .map(|&score| i32::from(score) - i32::from(i16::MIN))
        .collect();
    (
        finish_i16::<SPARSE_TRACE>(
            query,
            targets,
            &trace,
            band,
            i0,
            &subject_start,
            &scores,
            &best_col,
            &best_row,
        ),
        overflow_mask,
    )
}

#[allow(clippy::too_many_arguments)]
fn finish_i16<const SPARSE_TRACE: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    trace: &Trace16Matrix,
    band: usize,
    i0: i32,
    subject_start: &[i32; 16],
    scores: &[i32],
    best_col: &[usize; 16],
    best_row: &[usize; 16],
) -> Vec<SwResult> {
    targets
        .iter()
        .enumerate()
        .map(|(lane, target)| {
            if scores[lane] == 0 {
                return SwResult::default();
            }
            let mut column = best_col[lane];
            let mut row = best_row[lane];
            let mut i = i0 + column as i32 + row as i32 + 1;
            let mut j = subject_start[lane] + column as i32 + 1;
            let mut result = SwResult {
                score: scores[lane],
                query_end: i,
                subject_end: j,
                ..Default::default()
            };
            let mut operations = Vec::new();
            while i > 0 && j > 0 {
                let trace_index = (column + 1) * band + row;
                let (state, _, _) = trace.get::<SPARSE_TRACE>(trace_index, lane);
                if state == 0 {
                    break;
                }
                if state == 2 {
                    let mut count = 0;
                    loop {
                        count += 1;
                        i -= 1;
                        if row == 0 {
                            break;
                        }
                        row -= 1;
                        if i == 0 || trace.get::<SPARSE_TRACE>((column + 1) * band + row, lane).1 {
                            break;
                        }
                    }
                    push_operation_run(&mut operations, EditOperation::Insertion, count);
                    result.gap_openings += 1;
                    result.gaps += count;
                    result.length += count;
                } else if state == 3 {
                    let mut count = 0;
                    loop {
                        count += 1;
                        j -= 1;
                        if column == 0 || row + 1 >= band {
                            break;
                        }
                        column -= 1;
                        row += 1;
                        if j == 0 || trace.get::<SPARSE_TRACE>((column + 1) * band + row, lane).2 {
                            break;
                        }
                    }
                    push_operation_run(&mut operations, EditOperation::Deletion, count);
                    result.gap_openings += 1;
                    result.gaps += count;
                    result.length += count;
                } else {
                    if (query[(i - 1) as usize] & LETTER_MASK)
                        == (target.subject[(j - 1) as usize] & LETTER_MASK)
                    {
                        push_operation_run(&mut operations, EditOperation::Match, 1);
                        result.identities += 1;
                    } else {
                        push_operation_run(&mut operations, EditOperation::Substitution, 1);
                        result.mismatches += 1;
                    }
                    result.length += 1;
                    i -= 1;
                    j -= 1;
                    if column == 0 {
                        break;
                    }
                    column -= 1;
                }
            }
            result.query_begin = i;
            result.subject_begin = j;
            operations.reverse();
            result.operations = operations;
            result
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn finish_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    trace: &[Trace8],
    band: usize,
    i0: i32,
    subject_start: &[i32; 32],
    scores: &[i32],
    best_col: &[usize; 32],
    best_row: &[usize; 32],
) -> Vec<SwResult> {
    targets
        .iter()
        .enumerate()
        .map(|(lane, target)| {
            if scores[lane] == 0 {
                return SwResult::default();
            }
            let mut column = best_col[lane];
            let mut row = best_row[lane];
            let mut i = i0 + column as i32 + row as i32 + 1;
            let mut j = subject_start[lane] + column as i32 + 1;
            let mut result = SwResult {
                score: scores[lane],
                query_end: i,
                subject_end: j,
                ..Default::default()
            };
            let vmask = 1u64 << (lane + 32);
            let hmask = 1u64 << lane;
            let mut operations = Vec::new();
            while i > 0 && j > 0 {
                let trace_index = (column + 1) * band + row;
                let cell = trace[trace_index];
                let vertical_state = cell.gap & vmask != 0;
                let low_state = cell.gap & hmask != 0;
                if !vertical_state && !low_state {
                    break;
                }
                if vertical_state && !low_state {
                    let mut count = 0;
                    loop {
                        count += 1;
                        i -= 1;
                        if row == 0 {
                            break;
                        }
                        row -= 1;
                        if i == 0 || trace[(column + 1) * band + row].open & vmask != 0 {
                            break;
                        }
                    }
                    push_operation_run(&mut operations, EditOperation::Insertion, count);
                    result.gap_openings += 1;
                    result.gaps += count;
                    result.length += count;
                } else if vertical_state {
                    let mut count = 0;
                    loop {
                        count += 1;
                        j -= 1;
                        if column == 0 || row + 1 >= band {
                            break;
                        }
                        column -= 1;
                        row += 1;
                        if j == 0 || trace[(column + 1) * band + row].open & hmask != 0 {
                            break;
                        }
                    }
                    push_operation_run(&mut operations, EditOperation::Deletion, count);
                    result.gap_openings += 1;
                    result.gaps += count;
                    result.length += count;
                } else {
                    if (query[(i - 1) as usize] & LETTER_MASK)
                        == (target.subject[(j - 1) as usize] & LETTER_MASK)
                    {
                        push_operation_run(&mut operations, EditOperation::Match, 1);
                        result.identities += 1;
                    } else {
                        push_operation_run(&mut operations, EditOperation::Substitution, 1);
                        result.mismatches += 1;
                    }
                    result.length += 1;
                    i -= 1;
                    j -= 1;
                    if column == 0 {
                        break;
                    }
                    column -= 1;
                }
            }
            result.query_begin = i;
            result.subject_begin = j;
            operations.reverse();
            result.operations = operations;
            result
        })
        .collect()
}
