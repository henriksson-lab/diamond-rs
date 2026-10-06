//! Structural AVX2 ports of upstream `dp/swipe/banded_swipe.h` traceback
//! specializations. Keep the moving-band coordinates, vector profile, rolling
//! matrix, and packed lane masks aligned with the C++ implementation.

use super::simd_trace::TraceTarget;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{Letter, AMINO_ACID_COUNT, LETTER_MASK};
use crate::dp::smith_waterman::SwResult;
use crate::stats::score_matrix::ScoreMatrix;
use std::arch::x86_64 as arch;
use std::cell::RefCell;

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

#[derive(Default)]
struct TraceRowScratch {
    score8: Vec<V>,
    hgap8: Vec<V>,
    masks8: Vec<V>,
    score16: Vec<V>,
    hgap16: Vec<V>,
    masks16: Vec<V>,
    cbs: Vec<V>,
}

thread_local! {
    // Upstream's banded matrix keeps score and horizontal-gap rows in
    // thread-local MemBuffers. Trace masks remain call-local because they can
    // be much larger and retaining their peak would inflate steady-state RSS.
    static TRACE_ROW_SCRATCH: RefCell<TraceRowScratch> =
        RefCell::new(TraceRowScratch::default());
}

#[target_feature(enable = "avx2")]
unsafe fn standard_profile_i16(matrix: &[i8; 1024], subject: &[i8; 16]) -> [V; 32] {
    let subject = arch::_mm_loadu_si128(subject.as_ptr().cast());
    let indices = arch::_mm_and_si128(subject, arch::_mm_set1_epi8(15));
    let high = arch::_mm_cmpgt_epi8(subject, arch::_mm_set1_epi8(15));
    let mut profile = [arch::_mm256_setzero_si256(); 32];
    for (query_letter, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
        let row = matrix.as_ptr().add(query_letter * 32);
        let low_scores = arch::_mm_loadu_si128(row.cast());
        let high_scores = arch::_mm_loadu_si128(row.add(16).cast());
        let low_scores = arch::_mm_shuffle_epi8(low_scores, indices);
        let high_scores = arch::_mm_shuffle_epi8(high_scores, indices);
        let scores = arch::_mm_blendv_epi8(low_scores, high_scores, high);
        *slot = arch::_mm256_cvtepi8_epi16(scores);
    }
    profile
}

/// Upstream's AVX2 `SwipeProfile<int8_t>::set(SeqVector)`: preformatted
/// matrix halves and two byte shuffles replace a 32-by-lane scalar gather for
/// every target column.
#[target_feature(enable = "avx2")]
unsafe fn standard_profile_i8(
    matrix_low: &[i8; 1024],
    matrix_high: &[i8; 1024],
    subject: &[i8; 32],
) -> [V; 32] {
    let subject = arch::_mm256_and_si256(
        arch::_mm256_loadu_si256(subject.as_ptr().cast()),
        arch::_mm256_set1_epi8(LETTER_MASK as i8),
    );
    let high_mask = arch::_mm256_slli_epi16(
        arch::_mm256_and_si256(subject, arch::_mm256_set1_epi8(16)),
        3,
    );
    let low_index = arch::_mm256_or_si256(subject, high_mask);
    let high_index = arch::_mm256_or_si256(
        subject,
        arch::_mm256_xor_si256(high_mask, arch::_mm256_set1_epi8(i8::MIN)),
    );
    let mut profile = [arch::_mm256_setzero_si256(); 32];
    for (query_letter, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
        let low = arch::_mm256_loadu_si256(matrix_low.as_ptr().add(query_letter * 32).cast());
        let high = arch::_mm256_loadu_si256(matrix_high.as_ptr().add(query_letter * 32).cast());
        *slot = arch::_mm256_or_si256(
            arch::_mm256_shuffle_epi8(low, low_index),
            arch::_mm256_shuffle_epi8(high, high_index),
        );
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
    planes: Vec<u32>,
    stride: usize,
}

impl Trace16Matrix {
    fn new<const SPARSE: bool, const PLANES: bool>(cells: usize, lanes: usize) -> Self {
        if PLANES {
            Self {
                packed: Vec::new(),
                sparse: Vec::new(),
                planes: vec![0; cells],
                stride: 0,
            }
        } else if SPARSE {
            let stride = lanes.div_ceil(2);
            Self {
                packed: Vec::new(),
                sparse: vec![0; cells * stride],
                planes: Vec::new(),
                stride,
            }
        } else {
            Self {
                packed: vec![Trace16::default(); cells],
                sparse: Vec::new(),
                planes: Vec::new(),
                stride: 0,
            }
        }
    }

    #[inline]
    unsafe fn set<const SPARSE: bool, const PLANES: bool, const BMI2: bool>(
        &mut self,
        cell: usize,
        lanes: usize,
        gap_v: u32,
        gap_h: u32,
        open_v: u32,
        open_h: u32,
        active: u32,
    ) {
        if PLANES {
            debug_assert!(lanes <= 8);
            let planes = if BMI2 {
                pack_trace_planes_bmi2(gap_v, gap_h, open_v, open_h, active)
            } else {
                pack_trace_planes_portable(gap_v, gap_h, open_v, open_h, active)
            };
            *self.planes.get_unchecked_mut(cell) = planes;
        } else if SPARSE {
            let base = cell * self.stride;
            debug_assert!(lanes <= 6);
            let nibbles = if BMI2 {
                pack_trace_nibbles_bmi2(gap_v, gap_h, open_v, open_h, active)
            } else {
                pack_trace_nibbles_portable(gap_v, gap_h, open_v, open_h, active)
            };
            let destination = self.sparse.as_mut_ptr().add(base);
            match self.stride {
                1 => *destination = nibbles as u8,
                2 => std::ptr::write_unaligned(destination.cast::<u16>(), nibbles as u16),
                3 => {
                    std::ptr::write_unaligned(destination.cast::<u16>(), nibbles as u16);
                    *destination.add(2) = (nibbles >> 16) as u8;
                }
                _ => std::hint::unreachable_unchecked(),
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
    fn get<const SPARSE: bool, const PLANES: bool>(
        &self,
        cell: usize,
        lane: usize,
    ) -> (u8, bool, bool) {
        if PLANES {
            let planes = self.planes[cell];
            let bit = 1u32 << lane;
            let state = u8::from(planes & bit != 0) | (u8::from(planes & (bit << 8) != 0) << 1);
            (state, planes & (bit << 16) != 0, planes & (bit << 24) != 0)
        } else if SPARSE {
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

/// Collapse the duplicated two-bit lanes produced by `_mm256_movemask_epi8`
/// after an i16 comparison into one bit per lane.
#[inline]
fn compact_i16_movemask(mut bits: u32) -> u16 {
    bits &= 0x5555_5555;
    bits = (bits | bits >> 1) & 0x3333_3333;
    bits = (bits | bits >> 2) & 0x0f0f_0f0f;
    bits = (bits | bits >> 4) & 0x00ff_00ff;
    bits = (bits | bits >> 8) & 0x0000_ffff;
    bits as u16
}

/// Place the low eight lane bits at bit 0 of successive four-bit nibbles.
#[inline]
fn spread_nibbles(bits: u16) -> u64 {
    let mut bits = u64::from(bits & 0xff);
    bits = (bits | bits << 12) & 0x000f_000f;
    bits = (bits | bits << 6) & 0x0303_0303;
    (bits | bits << 3) & 0x1111_1111
}

#[inline]
fn pack_trace_nibbles_portable(
    gap_v: u32,
    gap_h: u32,
    open_v: u32,
    open_h: u32,
    active: u32,
) -> u64 {
    let active = compact_i16_movemask(active);
    let vertical = compact_i16_movemask(gap_v) & active;
    let horizontal = compact_i16_movemask(gap_h) & active & !vertical;
    let state_low = active & !vertical;
    let state_high = active & (vertical | horizontal);
    spread_nibbles(state_low)
        | spread_nibbles(state_high) << 1
        | spread_nibbles(compact_i16_movemask(open_v)) << 2
        | spread_nibbles(compact_i16_movemask(open_h)) << 3
}

#[target_feature(enable = "bmi2")]
unsafe fn pack_trace_nibbles_bmi2(
    gap_v: u32,
    gap_h: u32,
    open_v: u32,
    open_h: u32,
    active: u32,
) -> u64 {
    const EVEN_BYTES: u32 = 0x5555_5555;
    const NIBBLE_LOW_BITS: u64 = 0x1111_1111;
    let extract = |bits| arch::_pext_u32(bits, EVEN_BYTES) as u16;
    let active = extract(active);
    let vertical = extract(gap_v) & active;
    let horizontal = extract(gap_h) & active & !vertical;
    let state_low = active & !vertical;
    let state_high = active & (vertical | horizontal);
    arch::_pdep_u64(u64::from(state_low), NIBBLE_LOW_BITS)
        | arch::_pdep_u64(u64::from(state_high), NIBBLE_LOW_BITS) << 1
        | arch::_pdep_u64(u64::from(extract(open_v)), NIBBLE_LOW_BITS) << 2
        | arch::_pdep_u64(u64::from(extract(open_h)), NIBBLE_LOW_BITS) << 3
}

#[inline]
fn pack_trace_planes_portable(
    gap_v: u32,
    gap_h: u32,
    open_v: u32,
    open_h: u32,
    active: u32,
) -> u32 {
    let active = compact_i16_movemask(active);
    let vertical = compact_i16_movemask(gap_v) & active;
    let horizontal = compact_i16_movemask(gap_h) & active & !vertical;
    let state_low = active & !vertical;
    let state_high = active & (vertical | horizontal);
    u32::from(state_low & 0xff)
        | u32::from(state_high & 0xff) << 8
        | u32::from(compact_i16_movemask(open_v) & 0xff) << 16
        | u32::from(compact_i16_movemask(open_h) & 0xff) << 24
}

#[target_feature(enable = "bmi2")]
unsafe fn pack_trace_planes_bmi2(
    gap_v: u32,
    gap_h: u32,
    open_v: u32,
    open_h: u32,
    active: u32,
) -> u32 {
    const EVEN_BYTES: u32 = 0x5555_5555;
    let extract = |bits| arch::_pext_u32(bits, EVEN_BYTES) as u16;
    let active = extract(active);
    let vertical = extract(gap_v) & active;
    let horizontal = extract(gap_h) & active & !vertical;
    let state_low = active & !vertical;
    let state_high = active & (vertical | horizontal);
    u32::from(state_low & 0xff)
        | u32::from(state_high & 0xff) << 8
        | u32::from(extract(open_v) & 0xff) << 16
        | u32::from(extract(open_h) & 0xff) << 24
}

#[target_feature(enable = "avx2")]
pub(super) unsafe fn trace_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
    semi_global: bool,
) -> (Vec<SwResult>, u32) {
    TRACE_ROW_SCRATCH.with(|scratch| {
        if semi_global {
            trace_i8_with_scratch::<true>(
                query,
                targets,
                matrix,
                query_cbs,
                &mut scratch.borrow_mut(),
            )
        } else {
            trace_i8_with_scratch::<false>(
                query,
                targets,
                matrix,
                query_cbs,
                &mut scratch.borrow_mut(),
            )
        }
    })
}

#[target_feature(enable = "avx2")]
unsafe fn trace_i8_with_scratch<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut TraceRowScratch,
) -> (Vec<SwResult>, u32) {
    match (
        query_cbs.is_empty(),
        targets.iter().all(|target| target.matrix.is_none()),
    ) {
        (true, true) => {
            trace_i8_impl::<SEMI_GLOBAL, false, true>(query, targets, matrix, query_cbs, scratch)
        }
        (false, true) => {
            trace_i8_impl::<SEMI_GLOBAL, true, true>(query, targets, matrix, query_cbs, scratch)
        }
        (true, false) => {
            trace_i8_impl::<SEMI_GLOBAL, false, false>(query, targets, matrix, query_cbs, scratch)
        }
        (false, false) => {
            trace_i8_impl::<SEMI_GLOBAL, true, false>(query, targets, matrix, query_cbs, scratch)
        }
    }
}

#[target_feature(enable = "avx2")]
unsafe fn trace_i8_impl<const SEMI_GLOBAL: bool, const HAS_CBS: bool, const STANDARD_ONLY: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut TraceRowScratch,
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
    let delta = if SEMI_GLOBAL { 0 } else { i8::MIN };
    let zero = arch::_mm256_set1_epi8(delta);
    let TraceRowScratch {
        score8: score_row,
        hgap8: hgap_row,
        masks8: row_masks,
        cbs: cbs_vectors,
        ..
    } = scratch;
    score_row.clear();
    score_row.resize(band, zero);
    hgap_row.clear();
    hgap_row.resize(band + 1, zero);
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
    if HAS_CBS {
        cbs_vectors.clear();
        cbs_vectors.reserve(query.len());
        for &bias in query_cbs.get_unchecked(..query.len()) {
            cbs_vectors.push(arch::_mm256_and_si256(
                arch::_mm256_set1_epi8(bias),
                standard_mask,
            ));
        }
    }
    row_masks.clear();
    row_masks.reserve(band);
    for row in 0..band {
        let mut lanes = [0i8; LANES];
        for lane in 0..targets.len() {
            if row >= band_offset[lane] {
                lanes[lane] = -1;
            }
        }
        row_masks.push(arch::_mm256_loadu_si256(lanes.as_ptr().cast()));
    }
    let mut best = [delta; LANES];
    let mut best_col = [0usize; LANES];
    let mut best_row = [0usize; LANES];

    for column in 0..columns {
        let mut subject_letter = [0usize; LANES];
        let mut active_lanes = [0i8; LANES];
        for lane in 0..targets.len() {
            let pos = subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                let letter = targets[lane].subject[pos as usize];
                subject_letter[lane] = (letter & LETTER_MASK) as usize;
                active_lanes[lane] = -1;
            }
        }
        let active = arch::_mm256_loadu_si256(active_lanes.as_ptr().cast());
        let profile = if STANDARD_ONLY {
            standard_profile_i8(
                matrix.matrix8_low(),
                matrix.matrix8_high(),
                &subject_letter.map(|x| x as i8),
            )
        } else {
            let mut profile = [zero; 32];
            for (query_letter, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
                let mut scores = [0i8; LANES];
                for lane in 0..targets.len() {
                    if active_lanes[lane] == 0 {
                        continue;
                    }
                    scores[lane] = if let Some(adjusted) = targets[lane].matrix {
                        *adjusted
                            .scores
                            .get_unchecked(subject_letter[lane] * 32 + query_letter)
                    } else {
                        *matrix
                            .matrix8()
                            .get_unchecked(query_letter * 32 + subject_letter[lane])
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
        let mut row_counter = arch::_mm256_set1_epi8((query_begin - moving_i0) as i8);
        let mut row_max = zero;
        let mut row = (query_begin - moving_i0) as usize;
        let row_end = (query_end - moving_i0) as usize;
        let mut query_ptr = query.as_ptr().add(query_begin as usize);
        let mut cbs_ptr = if HAS_CBS {
            cbs_vectors.as_ptr().add(query_begin as usize)
        } else {
            std::ptr::NonNull::<V>::dangling().as_ptr()
        };
        while row < row_end {
            let cell_mask = arch::_mm256_and_si256(active, *row_masks.get_unchecked(row));
            let bias = if HAS_CBS {
                *cbs_ptr
            } else {
                arch::_mm256_setzero_si256()
            };
            let substitution = arch::_mm256_adds_epi8(
                *profile.get_unchecked((*query_ptr & LETTER_MASK) as usize),
                bias,
            );
            let diagonal = arch::_mm256_adds_epi8(*score_row.get_unchecked(row), substitution);
            let horizontal = *hgap_row.get_unchecked(row + 1);
            let vertical_before = vertical;
            let mut score = arch::_mm256_max_epi8(diagonal, horizontal);
            score = arch::_mm256_max_epi8(score, vertical_before);
            if !SEMI_GLOBAL {
                score = arch::_mm256_max_epi8(score, zero);
            }
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
            let active_bits = if SEMI_GLOBAL {
                arch::_mm256_movemask_epi8(cell_mask) as u32
            } else {
                arch::_mm256_movemask_epi8(arch::_mm256_cmpgt_epi8(score, zero)) as u32
            };
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
            *trace.get_unchecked_mut((column + 1) * band + row) = Trace8 {
                gap: (u64::from(high) << 32) | u64::from(low),
                open: (u64::from(open_v) << 32) | u64::from(open_h),
            };
            *score_row.get_unchecked_mut(row) = score;
            *hgap_row.get_unchecked_mut(row) =
                arch::_mm256_blendv_epi8(zero, next_horizontal, cell_mask);
            vertical = arch::_mm256_blendv_epi8(zero, next_vertical, cell_mask);
            query_ptr = query_ptr.add(1);
            if HAS_CBS {
                cbs_ptr = cbs_ptr.add(1);
            }
            row += 1;
        }
        let mut column_scores = [delta; LANES];
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
    for lane in 0..targets.len() {
        if best[lane] == i8::MAX {
            overflow_mask |= 1 << lane;
        }
    }
    let scores: Vec<i32> = best[..targets.len()]
        .iter()
        .map(|&score| i32::from(score) - i32::from(delta))
        .collect();
    (
        finish_i8::<SEMI_GLOBAL>(
            query,
            targets,
            matrix,
            query_cbs,
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
    semi_global: bool,
) -> (Vec<SwResult>, u32) {
    TRACE_ROW_SCRATCH.with(|scratch| {
        if semi_global {
            trace_i16_with_scratch::<true>(
                query,
                targets,
                matrix,
                query_cbs,
                &mut scratch.borrow_mut(),
            )
        } else {
            trace_i16_with_scratch::<false>(
                query,
                targets,
                matrix,
                query_cbs,
                &mut scratch.borrow_mut(),
            )
        }
    })
}

#[target_feature(enable = "avx2")]
unsafe fn trace_i16_with_scratch<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut TraceRowScratch,
) -> (Vec<SwResult>, u32) {
    if targets.len() <= 6 {
        if std::arch::is_x86_feature_detected!("bmi2") {
            trace_i16_dispatch::<SEMI_GLOBAL, true, false, true>(
                query, targets, matrix, query_cbs, scratch,
            )
        } else {
            trace_i16_dispatch::<SEMI_GLOBAL, true, false, false>(
                query, targets, matrix, query_cbs, scratch,
            )
        }
    } else if targets.len() <= 8 {
        if std::arch::is_x86_feature_detected!("bmi2") {
            trace_i16_dispatch::<SEMI_GLOBAL, true, true, true>(
                query, targets, matrix, query_cbs, scratch,
            )
        } else {
            trace_i16_dispatch::<SEMI_GLOBAL, true, true, false>(
                query, targets, matrix, query_cbs, scratch,
            )
        }
    } else {
        trace_i16_dispatch::<SEMI_GLOBAL, false, false, false>(
            query, targets, matrix, query_cbs, scratch,
        )
    }
}

#[target_feature(enable = "avx2")]
unsafe fn trace_i16_dispatch<
    const SEMI_GLOBAL: bool,
    const SPARSE_TRACE: bool,
    const PLANES: bool,
    const BMI2: bool,
>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut TraceRowScratch,
) -> (Vec<SwResult>, u32) {
    match (
        query_cbs.is_empty(),
        targets.iter().all(|target| target.matrix.is_none()),
    ) {
        (true, true) => trace_i16_impl::<SEMI_GLOBAL, SPARSE_TRACE, PLANES, BMI2, false, true>(
            query, targets, matrix, query_cbs, scratch,
        ),
        (false, true) => trace_i16_impl::<SEMI_GLOBAL, SPARSE_TRACE, PLANES, BMI2, true, true>(
            query, targets, matrix, query_cbs, scratch,
        ),
        (true, false) => trace_i16_impl::<SEMI_GLOBAL, SPARSE_TRACE, PLANES, BMI2, false, false>(
            query, targets, matrix, query_cbs, scratch,
        ),
        (false, false) => trace_i16_impl::<SEMI_GLOBAL, SPARSE_TRACE, PLANES, BMI2, true, false>(
            query, targets, matrix, query_cbs, scratch,
        ),
    }
}

#[target_feature(enable = "avx2")]
unsafe fn trace_i16_impl<
    const SEMI_GLOBAL: bool,
    const SPARSE_TRACE: bool,
    const PLANES: bool,
    const BMI2: bool,
    const HAS_CBS: bool,
    const STANDARD_ONLY: bool,
>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
    scratch: &mut TraceRowScratch,
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
    let delta = if SEMI_GLOBAL { 0 } else { i16::MIN };
    let zero = arch::_mm256_set1_epi16(delta);
    let TraceRowScratch {
        score16: score_row,
        hgap16: hgap_row,
        masks16: row_masks,
        cbs: cbs_vectors,
        ..
    } = scratch;
    score_row.clear();
    score_row.resize(band, zero);
    hgap_row.clear();
    hgap_row.resize(band + 1, zero);
    let mut trace = Trace16Matrix::new::<SPARSE_TRACE, PLANES>((columns + 1) * band, targets.len());
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
    if HAS_CBS {
        cbs_vectors.clear();
        cbs_vectors.reserve(query.len());
        for &bias in query_cbs.get_unchecked(..query.len()) {
            cbs_vectors.push(arch::_mm256_and_si256(
                arch::_mm256_set1_epi16(bias as i16),
                standard_mask,
            ));
        }
    }
    row_masks.clear();
    row_masks.reserve(band);
    for row in 0..band {
        let mut lanes = [0i16; LANES];
        for lane in 0..targets.len() {
            if row >= band_offset[lane] {
                lanes[lane] = -1;
            }
        }
        row_masks.push(arch::_mm256_loadu_si256(lanes.as_ptr().cast()));
    }
    let mut best = [delta; LANES];
    let mut best_col = [0usize; LANES];
    let mut best_row = [0usize; LANES];

    for column in 0..columns {
        let mut subject_letter = [0usize; LANES];
        let mut subject_bytes = [0i8; LANES];
        let mut active_lanes = [0i16; LANES];
        for lane in 0..targets.len() {
            let pos = subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                let letter = targets[lane].subject[pos as usize];
                subject_letter[lane] = (letter & LETTER_MASK) as usize;
                subject_bytes[lane] = (letter & LETTER_MASK) as i8;
                active_lanes[lane] = -1;
            }
        }
        let active = arch::_mm256_loadu_si256(active_lanes.as_ptr().cast());
        let profile = if STANDARD_ONLY {
            standard_profile_i16(matrix.matrix8(), &subject_bytes)
        } else {
            let mut profile = [zero; 32];
            for (query_letter, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
                let mut scores = [0i16; LANES];
                for lane in 0..targets.len() {
                    if active_lanes[lane] == 0 {
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
        let mut row = (query_begin - moving_i0) as usize;
        let row_end = (query_end - moving_i0) as usize;
        let mut query_ptr = query.as_ptr().add(query_begin as usize);
        let mut cbs_ptr = if HAS_CBS {
            cbs_vectors.as_ptr().add(query_begin as usize)
        } else {
            std::ptr::NonNull::<V>::dangling().as_ptr()
        };
        while row < row_end {
            // The moving-band clamps q to the query and row to 0..band; the
            // profile index is masked to 0..31. These are the same pointer
            // invariants used by upstream's matrix iterator.
            let cell_mask = arch::_mm256_and_si256(active, *row_masks.get_unchecked(row));
            let bias = if HAS_CBS {
                *cbs_ptr
            } else {
                arch::_mm256_setzero_si256()
            };
            let query_letter = *query_ptr & LETTER_MASK;
            let substitution =
                arch::_mm256_adds_epi16(*profile.get_unchecked(query_letter as usize), bias);
            let diagonal = arch::_mm256_adds_epi16(*score_row.get_unchecked(row), substitution);
            let horizontal = *hgap_row.get_unchecked(row + 1);
            let vertical_before = vertical;
            let mut score = arch::_mm256_max_epi16(diagonal, horizontal);
            score = arch::_mm256_max_epi16(score, vertical_before);
            if !SEMI_GLOBAL {
                score = arch::_mm256_max_epi16(score, zero);
            }
            score = arch::_mm256_blendv_epi8(zero, score, cell_mask);
            col_best = arch::_mm256_max_epi16(col_best, score);
            let at_column_max = arch::_mm256_cmpeq_epi16(col_best, score);
            row_max = arch::_mm256_blendv_epi8(row_max, row_counter, at_column_max);
            row_counter = arch::_mm256_add_epi16(row_counter, arch::_mm256_set1_epi16(1));
            let gap_v =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(score, vertical_before)) as u32;
            let gap_h =
                arch::_mm256_movemask_epi8(arch::_mm256_cmpeq_epi16(score, horizontal)) as u32;
            let active_bits = if SEMI_GLOBAL {
                arch::_mm256_movemask_epi8(cell_mask) as u32
            } else {
                arch::_mm256_movemask_epi8(arch::_mm256_cmpgt_epi16(score, zero)) as u32
            };
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
            trace.set::<SPARSE_TRACE, PLANES, BMI2>(
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
            query_ptr = query_ptr.add(1);
            if HAS_CBS {
                cbs_ptr = cbs_ptr.add(1);
            }
            row += 1;
        }
        let mut column_scores = [delta; LANES];
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
        .map(|&score| i32::from(score) - i32::from(delta))
        .collect();
    (
        finish_i16::<SEMI_GLOBAL, SPARSE_TRACE, PLANES>(
            query,
            targets,
            matrix,
            query_cbs,
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
fn finish_i16<const SEMI_GLOBAL: bool, const SPARSE_TRACE: bool, const PLANES: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
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
            // Saturated lanes are promoted to the next score-width bin and
            // upstream never enters traceback for them at this tier.
            if scores[lane] == 0 || scores[lane] == u16::MAX as i32 {
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
            let end_score = scores[lane];
            // Banded SWIPE supports target-adjusted matrices (unlike the
            // full-matrix traceback overload).  Its adjusted lanes use the
            // target's i8 substitution table, no query CBS, and scaled gaps.
            let gap_scale = if target.matrix.is_some() {
                target.matrix_scale.max(1)
            } else {
                1
            };
            let mut traceback_score = 0;
            let mut operations = Vec::new();
            while i > 0 && j > 0 && (!SEMI_GLOBAL || traceback_score < end_score) {
                let trace_index = (column + 1) * band + row;
                let (state, _, _) = trace.get::<SPARSE_TRACE, PLANES>(trace_index, lane);
                if !SEMI_GLOBAL && state == 0 {
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
                        if i == 0
                            || trace
                                .get::<SPARSE_TRACE, PLANES>((column + 1) * band + row, lane)
                                .1
                        {
                            break;
                        }
                    }
                    push_operation_run(&mut operations, EditOperation::Insertion, count);
                    result.gap_openings += 1;
                    result.gaps += count;
                    result.length += count;
                    traceback_score -=
                        (matrix.gap_open() + count * matrix.gap_extend()) * gap_scale;
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
                        if j == 0
                            || trace
                                .get::<SPARSE_TRACE, PLANES>((column + 1) * band + row, lane)
                                .2
                        {
                            break;
                        }
                    }
                    push_operation_run(&mut operations, EditOperation::Deletion, count);
                    result.gap_openings += 1;
                    result.gaps += count;
                    result.length += count;
                    traceback_score -=
                        (matrix.gap_open() + count * matrix.gap_extend()) * gap_scale;
                } else {
                    let query_index = (i - 1) as usize;
                    let query_letter = query[query_index] & LETTER_MASK;
                    let subject_letter = target.subject[(j - 1) as usize] & LETTER_MASK;
                    traceback_score += if let Some(adjusted) = target.matrix {
                        i32::from(
                            adjusted.scores[subject_letter as usize * 32 + query_letter as usize],
                        )
                    } else {
                        matrix.matrix32()[subject_letter as usize * 32 + query_letter as usize]
                            + i32::from(query_cbs.get(query_index).copied().unwrap_or(0))
                    };
                    if query_letter == subject_letter {
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
            assert_eq!(
                traceback_score, end_score,
                "traceback score does not match the endpoint score: operations={operations:?} endpoint=({i},{j})"
            );
            result.query_begin = i;
            result.subject_begin = j;
            operations.reverse();
            result.operations = operations;
            result
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn finish_i8<const SEMI_GLOBAL: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    query_cbs: &[i8],
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
            // Saturated lanes are promoted to the next score-width bin and
            // upstream never enters traceback for them at this tier.
            if scores[lane] == 0 || scores[lane] == u8::MAX as i32 {
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
            let end_score = scores[lane];
            // This is the adjusted-matrix branch of upstream banded_swipe.h,
            // not the unsupported adjusted branch in full_swipe.h.
            let gap_scale = if target.matrix.is_some() {
                target.matrix_scale.max(1)
            } else {
                1
            };
            let mut traceback_score = 0;
            let mut operations = Vec::new();
            while i > 0 && j > 0 && (!SEMI_GLOBAL || traceback_score < end_score) {
                let trace_index = (column + 1) * band + row;
                let cell = trace[trace_index];
                let vertical_state = cell.gap & vmask != 0;
                let low_state = cell.gap & hmask != 0;
                if !SEMI_GLOBAL && !vertical_state && !low_state {
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
                    traceback_score -=
                        (matrix.gap_open() + count * matrix.gap_extend()) * gap_scale;
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
                    traceback_score -=
                        (matrix.gap_open() + count * matrix.gap_extend()) * gap_scale;
                } else {
                    let query_index = (i - 1) as usize;
                    let query_letter = query[query_index] & LETTER_MASK;
                    let subject_letter = target.subject[(j - 1) as usize] & LETTER_MASK;
                    traceback_score += if let Some(adjusted) = target.matrix {
                        i32::from(
                            adjusted.scores[subject_letter as usize * 32 + query_letter as usize],
                        )
                    } else {
                        matrix.matrix32()[subject_letter as usize * 32 + query_letter as usize]
                            + i32::from(query_cbs.get(query_index).copied().unwrap_or(0))
                    };
                    if query_letter == subject_letter {
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
            assert_eq!(
                traceback_score, end_score,
                "traceback score does not match the endpoint score: operations={operations:?} endpoint=({i},{j})"
            );
            result.query_begin = i;
            result.subject_begin = j;
            operations.reverse();
            result.operations = operations;
            result
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bmi2_trace_packings_match_portable() {
        if !std::arch::is_x86_feature_detected!("bmi2") {
            return;
        }
        let mut state = 0x9e37_79b9_u32;
        let mut next = || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            // Model an i16 comparison movemask: both bytes of each lane
            // always carry the same comparison bit.
            let lanes = state & 0x5555_5555;
            lanes | lanes << 1
        };
        for _ in 0..10_000 {
            let gap_v = next();
            let gap_h = next();
            let open_v = next();
            let open_h = next();
            let active = next();
            let portable = pack_trace_planes_portable(gap_v, gap_h, open_v, open_h, active);
            let bmi2 = unsafe { pack_trace_planes_bmi2(gap_v, gap_h, open_v, open_h, active) };
            assert_eq!(bmi2, portable);
            let portable = pack_trace_nibbles_portable(gap_v, gap_h, open_v, open_h, active);
            let bmi2 = unsafe { pack_trace_nibbles_bmi2(gap_v, gap_h, open_v, open_h, active) };
            assert_eq!(bmi2, portable);
        }
    }
}
