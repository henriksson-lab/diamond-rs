//! Trace-free SWIPE score/end-coordinate kernels.
//!
//! This is the no-trace `SwipeConfig<false, VectorRowCounter<...>>` path from
//! upstream `banded_swipe.h`: byte/word lanes retain one in-place score row,
//! one in-place horizontal-gap row, and SIMD row counters.  No traceback or
//! per-cell statistics storage is allocated.

use super::simd_trace::TraceTarget;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::basic::value::AMINO_ACID_COUNT;
use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

#[cfg(target_arch = "x86")]
use std::arch::x86 as arch;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64 as arch;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Endpoint {
    pub score: i32,
    pub query_end: i32,
    pub subject_end: i32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndpointBatch {
    pub endpoints: [Endpoint; 32],
    pub overflow_mask: u32,
    pub len: usize,
}

impl EndpointBatch {
    fn empty(len: usize) -> Self {
        Self {
            endpoints: [Endpoint::default(); 32],
            overflow_mask: 0,
            len,
        }
    }
}

#[derive(Default)]
pub struct EndpointScratch {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    score8: Vec<arch::__m256i>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    hgap8: Vec<arch::__m256i>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    score16: Vec<arch::__m256i>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    hgap16: Vec<arch::__m256i>,
    i32_score: Vec<i32>,
    i32_hgap: Vec<i32>,
}

pub fn endpoint_batch_avx2_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    scratch: &mut EndpointScratch,
) -> Option<EndpointBatch> {
    if targets.is_empty()
        || targets.len() > 32
        || (!cbs.is_empty() && cbs.len() < query.len())
        || targets.iter().any(|target| target.d_end <= target.d_begin)
    {
        return None;
    }
    let max_band = if semi_global {
        i8::MAX as i32
    } else {
        u8::MAX as i32
    };
    let band = targets
        .iter()
        .map(|target| target.d_end - target.d_begin)
        .max()
        .unwrap_or(0);
    assert!(band <= max_band, "Band size exceeds row counter maximum.");
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        let has_cbs = !cbs.is_empty();
        let standard = targets.iter().all(|target| target.matrix.is_none());
        // SAFETY: AVX2 is runtime-detected; kernels size all retained rows and
        // fixed lane arrays before unchecked access.
        Some(unsafe {
            match (semi_global, has_cbs, standard) {
                (true, true, true) => {
                    endpoint_i8_impl::<true, true, true>(query, targets, matrix, cbs, scratch)
                }
                (true, true, false) => {
                    endpoint_i8_impl::<true, true, false>(query, targets, matrix, cbs, scratch)
                }
                (true, false, true) => {
                    endpoint_i8_impl::<true, false, true>(query, targets, matrix, cbs, scratch)
                }
                (true, false, false) => {
                    endpoint_i8_impl::<true, false, false>(query, targets, matrix, cbs, scratch)
                }
                (false, true, true) => {
                    endpoint_i8_impl::<false, true, true>(query, targets, matrix, cbs, scratch)
                }
                (false, true, false) => {
                    endpoint_i8_impl::<false, true, false>(query, targets, matrix, cbs, scratch)
                }
                (false, false, true) => {
                    endpoint_i8_impl::<false, false, true>(query, targets, matrix, cbs, scratch)
                }
                (false, false, false) => {
                    endpoint_i8_impl::<false, false, false>(query, targets, matrix, cbs, scratch)
                }
            }
        })
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, matrix, cbs, semi_global, scratch);
        None
    }
}

pub fn endpoint_batch_avx2_i16(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    scratch: &mut EndpointScratch,
) -> Option<EndpointBatch> {
    if targets.is_empty()
        || targets.len() > 16
        || (!cbs.is_empty() && cbs.len() < query.len())
        || targets.iter().any(|target| target.d_end <= target.d_begin)
    {
        return None;
    }
    let max_band = if semi_global {
        i16::MAX as i32
    } else {
        u16::MAX as i32
    };
    let band = targets
        .iter()
        .map(|target| target.d_end - target.d_begin)
        .max()
        .unwrap_or(0);
    assert!(band <= max_band, "Band size exceeds row counter maximum.");
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        let has_cbs = !cbs.is_empty();
        let standard = targets.iter().all(|target| target.matrix.is_none());
        // SAFETY: see the byte dispatcher above.
        Some(unsafe {
            match (semi_global, has_cbs, standard) {
                (true, true, true) => {
                    endpoint_i16_impl::<true, true, true>(query, targets, matrix, cbs, scratch)
                }
                (true, true, false) => {
                    endpoint_i16_impl::<true, true, false>(query, targets, matrix, cbs, scratch)
                }
                (true, false, true) => {
                    endpoint_i16_impl::<true, false, true>(query, targets, matrix, cbs, scratch)
                }
                (true, false, false) => {
                    endpoint_i16_impl::<true, false, false>(query, targets, matrix, cbs, scratch)
                }
                (false, true, true) => {
                    endpoint_i16_impl::<false, true, true>(query, targets, matrix, cbs, scratch)
                }
                (false, true, false) => {
                    endpoint_i16_impl::<false, true, false>(query, targets, matrix, cbs, scratch)
                }
                (false, false, true) => {
                    endpoint_i16_impl::<false, false, true>(query, targets, matrix, cbs, scratch)
                }
                (false, false, false) => {
                    endpoint_i16_impl::<false, false, false>(query, targets, matrix, cbs, scratch)
                }
            }
        })
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (query, targets, matrix, cbs, semi_global, scratch);
        None
    }
}

/// Upstream `full_swipe.h` no-trace byte tier. Unlike the banded kernel this
/// walks exactly `query.len()` rows and starts every subject lane at column 0.
pub fn endpoint_full_batch_avx2_i8(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    matrix_scale: i32,
    scratch: &mut EndpointScratch,
) -> Option<EndpointBatch> {
    if targets.is_empty() || targets.len() > 32 || (!cbs.is_empty() && cbs.len() < query.len()) {
        return None;
    }
    let max_rows = if semi_global {
        i8::MAX as usize
    } else {
        u8::MAX as usize
    };
    assert!(
        query.len() <= max_rows,
        "Query length exceeds row counter maximum."
    );
    assert!(matrix_scale == 1, "Matrix scale != 1.0 not supported.");
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        let standard = targets.iter().all(|target| target.matrix.is_none());
        // SAFETY: AVX2 was detected and fixed arrays cover all lanes.
        Some(unsafe {
            match (semi_global, !cbs.is_empty(), standard) {
                (true, true, true) => {
                    full_i8_impl::<true, true, true>(query, targets, matrix, cbs, scratch)
                }
                (true, true, false) => {
                    full_i8_impl::<true, true, false>(query, targets, matrix, cbs, scratch)
                }
                (true, false, true) => {
                    full_i8_impl::<true, false, true>(query, targets, matrix, cbs, scratch)
                }
                (true, false, false) => {
                    full_i8_impl::<true, false, false>(query, targets, matrix, cbs, scratch)
                }
                (false, true, true) => {
                    full_i8_impl::<false, true, true>(query, targets, matrix, cbs, scratch)
                }
                (false, true, false) => {
                    full_i8_impl::<false, true, false>(query, targets, matrix, cbs, scratch)
                }
                (false, false, true) => {
                    full_i8_impl::<false, false, true>(query, targets, matrix, cbs, scratch)
                }
                (false, false, false) => {
                    full_i8_impl::<false, false, false>(query, targets, matrix, cbs, scratch)
                }
            }
        })
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (
            query,
            targets,
            matrix,
            cbs,
            semi_global,
            matrix_scale,
            scratch,
        );
        None
    }
}

pub fn endpoint_full_batch_avx2_i16(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    semi_global: bool,
    matrix_scale: i32,
    scratch: &mut EndpointScratch,
) -> Option<EndpointBatch> {
    if targets.is_empty() || targets.len() > 16 || (!cbs.is_empty() && cbs.len() < query.len()) {
        return None;
    }
    let max_rows = if semi_global {
        i16::MAX as usize
    } else {
        u16::MAX as usize
    };
    assert!(
        query.len() <= max_rows,
        "Query length exceeds row counter maximum."
    );
    assert!(matrix_scale == 1, "Matrix scale != 1.0 not supported.");
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return None;
        }
        let standard = targets.iter().all(|target| target.matrix.is_none());
        // SAFETY: AVX2 was detected and fixed arrays cover all lanes.
        Some(unsafe {
            match (semi_global, !cbs.is_empty(), standard) {
                (true, true, true) => {
                    full_i16_impl::<true, true, true>(query, targets, matrix, cbs, scratch)
                }
                (true, true, false) => {
                    full_i16_impl::<true, true, false>(query, targets, matrix, cbs, scratch)
                }
                (true, false, true) => {
                    full_i16_impl::<true, false, true>(query, targets, matrix, cbs, scratch)
                }
                (true, false, false) => {
                    full_i16_impl::<true, false, false>(query, targets, matrix, cbs, scratch)
                }
                (false, true, true) => {
                    full_i16_impl::<false, true, true>(query, targets, matrix, cbs, scratch)
                }
                (false, true, false) => {
                    full_i16_impl::<false, true, false>(query, targets, matrix, cbs, scratch)
                }
                (false, false, true) => {
                    full_i16_impl::<false, false, true>(query, targets, matrix, cbs, scratch)
                }
                (false, false, false) => {
                    full_i16_impl::<false, false, false>(query, targets, matrix, cbs, scratch)
                }
            }
        })
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (
            query,
            targets,
            matrix,
            cbs,
            semi_global,
            matrix_scale,
            scratch,
        );
        None
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[derive(Clone, Copy)]
struct Geometry<const LANES: usize> {
    band: usize,
    i0: i32,
    i1: i32,
    subject_start: [i32; LANES],
    band_offset: [usize; LANES],
    columns: usize,
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn geometry<const LANES: usize>(query_len: usize, targets: &[TraceTarget<'_>]) -> Geometry<LANES> {
    let band = targets
        .iter()
        .map(|target| (target.d_end - target.d_begin) as usize)
        .max()
        .unwrap_or(0);
    let i1 = targets
        .iter()
        .map(|target| (target.d_end - 1).max(0))
        .min()
        .unwrap_or(0);
    let i0 = i1 + 1 - band as i32;
    let mut subject_start = [0; LANES];
    let mut band_offset = [0; LANES];
    let mut columns = 0;
    for (lane, target) in targets.iter().enumerate() {
        let expanded_begin = target.d_end - band as i32;
        subject_start[lane] = i1 - (target.d_end - 1);
        band_offset[lane] = (target.d_begin - expanded_begin).max(0) as usize;
        let subject_end =
            ((query_len as i32 - 1 - expanded_begin).min(target.subject.len() as i32 - 1) + 1)
                .max(0);
        columns = columns.max((subject_end - subject_start[lane]).max(0) as usize);
    }
    Geometry {
        band,
        i0,
        i1,
        subject_start,
        band_offset,
        columns,
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn standard_profile_i8(
    low: &[i8; 1024],
    high: &[i8; 1024],
    subject: &[i8; 32],
) -> [arch::__m256i; 32] {
    let subject = arch::_mm256_loadu_si256(subject.as_ptr().cast());
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
        let lo = arch::_mm256_loadu_si256(low.as_ptr().add(query_letter * 32).cast());
        let hi = arch::_mm256_loadu_si256(high.as_ptr().add(query_letter * 32).cast());
        *slot = arch::_mm256_or_si256(
            arch::_mm256_shuffle_epi8(lo, low_index),
            arch::_mm256_shuffle_epi8(hi, high_index),
        );
    }
    profile
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn full_i8_impl<const SEMI_GLOBAL: bool, const HAS_CBS: bool, const STANDARD_ONLY: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    scratch: &mut EndpointScratch,
) -> EndpointBatch {
    const LANES: usize = 32;
    let delta = if SEMI_GLOBAL { 0 } else { i8::MIN };
    let zero = arch::_mm256_set1_epi8(delta);
    scratch.score8.resize(query.len() + 1, zero);
    scratch.hgap8.resize(query.len() + 1, zero);
    scratch.score8.fill(zero);
    scratch.hgap8.fill(zero);
    let go = arch::_mm256_set1_epi8((matrix.gap_open() + matrix.gap_extend()) as i8);
    let ge = arch::_mm256_set1_epi8(matrix.gap_extend() as i8);
    let mut best = [delta; LANES];
    let mut best_query_end = [0i32; LANES];
    let mut best_subject_end = [0i32; LANES];
    let mut overflow = 0u32;
    let columns = targets
        .iter()
        .map(|target| target.subject.len())
        .max()
        .unwrap_or(0);
    for column in 0..columns {
        let mut subject = [0i8; LANES];
        let mut live_add = [i8::MIN; LANES];
        let mut live = [false; LANES];
        for lane in 0..targets.len() {
            if let Some(&letter) = targets[lane].subject.get(column) {
                subject[lane] = (letter & LETTER_MASK) as i8;
                live_add[lane] = 0;
                live[lane] = true;
            }
        }
        let live_mask = arch::_mm256_loadu_si256(live_add.as_ptr().cast());
        let profile = if STANDARD_ONLY {
            standard_profile_i8(matrix.matrix8_low(), matrix.matrix8_high(), &subject)
        } else {
            let mut profile = [zero; 32];
            for (ql, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
                let mut scores = [0i8; LANES];
                for lane in 0..targets.len() {
                    if !live[lane] {
                        continue;
                    }
                    let sl = subject[lane] as usize;
                    scores[lane] = if let Some(adjusted) = targets[lane].matrix {
                        adjusted.scores[sl * 32 + ql]
                    } else {
                        matrix.matrix8()[ql * 32 + sl]
                    };
                }
                *slot = arch::_mm256_loadu_si256(scores.as_ptr().cast());
            }
            profile
        };
        let mut diagonal = scratch.score8[0];
        scratch.score8[0] = zero;
        scratch.hgap8[0] = zero;
        let mut vertical = arch::_mm256_adds_epi8(zero, live_mask);
        let mut col_best = zero;
        let mut row_counter = arch::_mm256_set1_epi8(delta);
        let mut row_max = zero;
        for (qpos, &ql) in query.iter().enumerate() {
            let next_diagonal = scratch.score8[qpos + 1];
            let mut substitution = profile[(ql & LETTER_MASK) as usize];
            substitution = arch::_mm256_adds_epi8(substitution, live_mask);
            if HAS_CBS {
                substitution =
                    arch::_mm256_adds_epi8(substitution, arch::_mm256_set1_epi8(cbs[qpos]));
            }
            let diagonal_score = arch::_mm256_adds_epi8(diagonal, substitution);
            let horizontal = arch::_mm256_adds_epi8(scratch.hgap8[qpos + 1], live_mask);
            let mut score = arch::_mm256_max_epi8(diagonal_score, horizontal);
            score = arch::_mm256_max_epi8(score, vertical);
            col_best = arch::_mm256_max_epi8(col_best, score);
            row_max = arch::_mm256_blendv_epi8(
                row_max,
                row_counter,
                arch::_mm256_cmpeq_epi8(col_best, score),
            );
            row_counter = arch::_mm256_adds_epi8(row_counter, arch::_mm256_set1_epi8(1));
            let opened = arch::_mm256_subs_epi8(score, go);
            scratch.hgap8[qpos + 1] =
                arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, ge), opened);
            vertical = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical, ge), opened);
            scratch.score8[qpos + 1] = score;
            diagonal = next_diagonal;
        }
        let mut scores = [delta; LANES];
        let mut rows = [delta; LANES];
        arch::_mm256_storeu_si256(scores.as_mut_ptr().cast(), col_best);
        arch::_mm256_storeu_si256(rows.as_mut_ptr().cast(), row_max);
        for lane in 0..targets.len() {
            if live[lane] && scores[lane] > best[lane] {
                best[lane] = scores[lane];
                best_query_end[lane] = scores_decode_i8(rows[lane], delta) + 1;
                best_subject_end[lane] = column as i32 + 1;
            }
        }
    }
    let mut out = EndpointBatch::empty(targets.len());
    for lane in 0..targets.len() {
        out.endpoints[lane] = Endpoint {
            score: scores_decode_i8(best[lane], delta),
            query_end: best_query_end[lane],
            subject_end: best_subject_end[lane],
        };
        if best[lane] == i8::MAX {
            overflow |= 1 << lane;
        }
    }
    out.overflow_mask = overflow;
    out
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn full_i16_impl<const SEMI_GLOBAL: bool, const HAS_CBS: bool, const STANDARD_ONLY: bool>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    scratch: &mut EndpointScratch,
) -> EndpointBatch {
    const LANES: usize = 16;
    let delta = if SEMI_GLOBAL { 0 } else { i16::MIN };
    let zero = arch::_mm256_set1_epi16(delta);
    scratch.score16.resize(query.len() + 1, zero);
    scratch.hgap16.resize(query.len() + 1, zero);
    scratch.score16.fill(zero);
    scratch.hgap16.fill(zero);
    let go = arch::_mm256_set1_epi16((matrix.gap_open() + matrix.gap_extend()) as i16);
    let ge = arch::_mm256_set1_epi16(matrix.gap_extend() as i16);
    let mut best = [delta; LANES];
    let mut best_query_end = [0i32; LANES];
    let mut best_subject_end = [0i32; LANES];
    let mut overflow = 0u32;
    let columns = targets
        .iter()
        .map(|target| target.subject.len())
        .max()
        .unwrap_or(0);
    for column in 0..columns {
        let mut subject = [0usize; LANES];
        let mut live_add = [i16::MIN; LANES];
        let mut live = [false; LANES];
        for lane in 0..targets.len() {
            if let Some(&letter) = targets[lane].subject.get(column) {
                subject[lane] = (letter & LETTER_MASK) as usize;
                live_add[lane] = 0;
                live[lane] = true;
            }
        }
        let live_mask = arch::_mm256_loadu_si256(live_add.as_ptr().cast());
        let mut profile = [zero; 32];
        for (ql, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
            let mut scores = [0i16; LANES];
            for lane in 0..targets.len() {
                if !live[lane] {
                    continue;
                }
                let sl = subject[lane];
                scores[lane] = if STANDARD_ONLY {
                    matrix.matrix16()[ql * 32 + sl]
                } else if let Some(adjusted) = targets[lane].matrix {
                    adjusted.scores[sl * 32 + ql] as i16
                } else {
                    matrix.matrix16()[ql * 32 + sl]
                };
            }
            *slot = arch::_mm256_loadu_si256(scores.as_ptr().cast());
        }
        let mut diagonal = scratch.score16[0];
        scratch.score16[0] = zero;
        scratch.hgap16[0] = zero;
        let mut vertical = arch::_mm256_adds_epi16(zero, live_mask);
        let mut col_best = zero;
        let mut row_counter = arch::_mm256_set1_epi16(delta);
        let mut row_max = zero;
        for (qpos, &ql) in query.iter().enumerate() {
            let next_diagonal = scratch.score16[qpos + 1];
            let mut substitution = profile[(ql & LETTER_MASK) as usize];
            substitution = arch::_mm256_adds_epi16(substitution, live_mask);
            if HAS_CBS {
                substitution = arch::_mm256_adds_epi16(
                    substitution,
                    arch::_mm256_set1_epi16(cbs[qpos] as i16),
                );
            }
            let diagonal_score = arch::_mm256_adds_epi16(diagonal, substitution);
            let horizontal = arch::_mm256_adds_epi16(scratch.hgap16[qpos + 1], live_mask);
            let mut score = arch::_mm256_max_epi16(diagonal_score, horizontal);
            score = arch::_mm256_max_epi16(score, vertical);
            col_best = arch::_mm256_max_epi16(col_best, score);
            row_max = arch::_mm256_blendv_epi8(
                row_max,
                row_counter,
                arch::_mm256_cmpeq_epi16(col_best, score),
            );
            row_counter = arch::_mm256_adds_epi16(row_counter, arch::_mm256_set1_epi16(1));
            let opened = arch::_mm256_subs_epi16(score, go);
            scratch.hgap16[qpos + 1] =
                arch::_mm256_max_epi16(arch::_mm256_subs_epi16(horizontal, ge), opened);
            vertical = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(vertical, ge), opened);
            scratch.score16[qpos + 1] = score;
            diagonal = next_diagonal;
        }
        let mut scores = [delta; LANES];
        let mut rows = [delta; LANES];
        arch::_mm256_storeu_si256(scores.as_mut_ptr().cast(), col_best);
        arch::_mm256_storeu_si256(rows.as_mut_ptr().cast(), row_max);
        for lane in 0..targets.len() {
            if live[lane] && scores[lane] > best[lane] {
                best[lane] = scores[lane];
                best_query_end[lane] = rows[lane] as i32 - delta as i32 + 1;
                best_subject_end[lane] = column as i32 + 1;
            }
        }
    }
    let mut out = EndpointBatch::empty(targets.len());
    for lane in 0..targets.len() {
        out.endpoints[lane] = Endpoint {
            score: best[lane] as i32 - delta as i32,
            query_end: best_query_end[lane],
            subject_end: best_subject_end[lane],
        };
        if best[lane] == i16::MAX {
            overflow |= 1 << lane;
        }
    }
    out.overflow_mask = overflow;
    out
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn endpoint_i8_impl<
    const SEMI_GLOBAL: bool,
    const HAS_CBS: bool,
    const STANDARD_ONLY: bool,
>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    scratch: &mut EndpointScratch,
) -> EndpointBatch {
    const LANES: usize = 32;
    let g = geometry::<LANES>(query.len(), targets);
    let delta = if SEMI_GLOBAL { 0 } else { i8::MIN };
    let zero = arch::_mm256_set1_epi8(delta);
    scratch.score8.resize(g.band, zero);
    scratch.hgap8.resize(g.band + 1, zero);
    scratch.score8.fill(zero);
    scratch.hgap8.fill(zero);

    let mut gap_open = [0i8; LANES];
    let mut gap_extend = [0i8; LANES];
    let mut overflow = 0u32;
    let mut standard_lanes = [0i8; LANES];
    for lane in 0..targets.len() {
        let scale = if targets[lane].matrix.is_some() {
            targets[lane].matrix_scale.max(1)
        } else {
            1
        };
        let go = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=i8::MAX as i32).contains(&go) || !(0..=i8::MAX as i32).contains(&ge) {
            overflow |= 1 << lane;
        }
        gap_open[lane] = go.clamp(0, i8::MAX as i32) as i8;
        gap_extend[lane] = ge.clamp(0, i8::MAX as i32) as i8;
        if targets[lane].matrix.is_none() {
            standard_lanes[lane] = -1;
        }
    }
    let go = arch::_mm256_loadu_si256(gap_open.as_ptr().cast());
    let ge = arch::_mm256_loadu_si256(gap_extend.as_ptr().cast());
    let standard_mask = arch::_mm256_loadu_si256(standard_lanes.as_ptr().cast());

    let mut sorted_offsets = g.band_offset;
    sorted_offsets[..targets.len()].sort_unstable();
    let mut part_bounds = [0usize; 34];
    let mut part_masks = [arch::_mm256_setzero_si256(); 33];
    let mut part_count = 0;
    let mut part_begin = 0;
    for &offset in &sorted_offsets[..targets.len()] {
        let offset = offset.min(g.band);
        if offset <= part_begin || offset == g.band {
            continue;
        }
        part_bounds[part_count] = part_begin;
        let mut lanes = [i8::MIN; LANES];
        for lane in 0..targets.len() {
            if part_begin >= g.band_offset[lane] {
                lanes[lane] = 0;
            }
        }
        part_masks[part_count] = arch::_mm256_loadu_si256(lanes.as_ptr().cast());
        part_count += 1;
        part_begin = offset;
    }
    part_bounds[part_count] = part_begin;
    let mut lanes = [i8::MIN; LANES];
    for lane in 0..targets.len() {
        if part_begin >= g.band_offset[lane] {
            lanes[lane] = 0;
        }
    }
    part_masks[part_count] = arch::_mm256_loadu_si256(lanes.as_ptr().cast());
    part_count += 1;
    part_bounds[part_count] = g.band;

    let mut best = [delta; LANES];
    let mut best_col = [0usize; LANES];
    let mut best_row = [0i32; LANES];
    for column in 0..g.columns {
        let mut subject = [0i8; LANES];
        let mut live_add = [i8::MIN; LANES];
        let mut live = [false; LANES];
        for lane in 0..targets.len() {
            let pos = g.subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                subject[lane] = (targets[lane].subject[pos as usize] & LETTER_MASK) as i8;
                live_add[lane] = 0;
                live[lane] = true;
            }
        }
        let live_mask = arch::_mm256_loadu_si256(live_add.as_ptr().cast());
        let profile = if STANDARD_ONLY {
            standard_profile_i8(matrix.matrix8_low(), matrix.matrix8_high(), &subject)
        } else {
            let mut profile = [zero; 32];
            for (ql, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
                let mut scores = [0i8; LANES];
                for lane in 0..targets.len() {
                    if !live[lane] {
                        continue;
                    }
                    let sl = subject[lane] as usize;
                    scores[lane] = if let Some(adjusted) = targets[lane].matrix {
                        adjusted.scores[sl * 32 + ql]
                    } else {
                        matrix.matrix8()[ql * 32 + sl]
                    };
                }
                *slot = arch::_mm256_loadu_si256(scores.as_ptr().cast());
            }
            profile
        };
        let moving_i0 = g.i0 + column as i32;
        let query_begin = moving_i0.max(0);
        let query_end = (g.i1 + column as i32).min(query.len() as i32 - 1) + 1;
        if query_begin >= query_end {
            continue;
        }
        let active_begin = (query_begin - moving_i0) as usize;
        let active_end = (query_end - moving_i0) as usize;
        let mut vertical = zero;
        let mut col_best = zero;
        // C++ integer-promotes `zero_score() + Score(i)` and then narrows
        // into the lane. This wraps the shifted coordinate representation;
        // subsequent `VectorRowCounter::inc` operations are saturating.
        let initial_row = (delta as i32 + active_begin as i32) as i8;
        let mut row_counter = arch::_mm256_set1_epi8(initial_row);
        let mut row_max = zero;
        for part in 0..part_count {
            let row_begin = part_bounds[part].max(active_begin);
            let row_end = part_bounds[part + 1].min(active_end);
            if row_begin >= row_end {
                continue;
            }
            let cell_add = arch::_mm256_adds_epi8(part_masks[part], live_mask);
            vertical = arch::_mm256_adds_epi8(vertical, cell_add);
            for row in row_begin..row_end {
                let qpos = (moving_i0 + row as i32) as usize;
                let mut substitution = *profile.get_unchecked((query[qpos] & LETTER_MASK) as usize);
                substitution = arch::_mm256_adds_epi8(substitution, cell_add);
                if HAS_CBS {
                    let bias = arch::_mm256_and_si256(
                        arch::_mm256_set1_epi8(*cbs.get_unchecked(qpos)),
                        standard_mask,
                    );
                    substitution = arch::_mm256_adds_epi8(substitution, bias);
                }
                let diagonal =
                    arch::_mm256_adds_epi8(*scratch.score8.get_unchecked(row), substitution);
                let horizontal =
                    arch::_mm256_adds_epi8(*scratch.hgap8.get_unchecked(row + 1), cell_add);
                let mut score = arch::_mm256_max_epi8(diagonal, horizontal);
                score = arch::_mm256_max_epi8(score, vertical);
                col_best = arch::_mm256_max_epi8(col_best, score);
                row_max = arch::_mm256_blendv_epi8(
                    row_max,
                    row_counter,
                    arch::_mm256_cmpeq_epi8(col_best, score),
                );
                row_counter = arch::_mm256_adds_epi8(row_counter, arch::_mm256_set1_epi8(1));
                let opened = arch::_mm256_subs_epi8(score, go);
                let next_horizontal =
                    arch::_mm256_max_epi8(arch::_mm256_subs_epi8(horizontal, ge), opened);
                vertical = arch::_mm256_max_epi8(arch::_mm256_subs_epi8(vertical, ge), opened);
                *scratch.score8.get_unchecked_mut(row) = score;
                *scratch.hgap8.get_unchecked_mut(row) = next_horizontal;
            }
        }
        let mut scores = [delta; LANES];
        let mut rows = [delta; LANES];
        arch::_mm256_storeu_si256(scores.as_mut_ptr().cast(), col_best);
        arch::_mm256_storeu_si256(rows.as_mut_ptr().cast(), row_max);
        for lane in 0..targets.len() {
            if live[lane] && scores[lane] > best[lane] {
                best[lane] = scores[lane];
                best_col[lane] = column;
                best_row[lane] = scores_decode_i8(rows[lane], delta);
            }
        }
    }
    let mut out = EndpointBatch::empty(targets.len());
    for lane in 0..targets.len() {
        let score = scores_decode_i8(best[lane], delta);
        out.endpoints[lane] = Endpoint {
            score,
            query_end: g.i0 + best_col[lane] as i32 + best_row[lane] + 1,
            subject_end: g.subject_start[lane] + best_col[lane] as i32 + 1,
        };
        if best[lane] == i8::MAX {
            overflow |= 1 << lane;
        }
    }
    out.overflow_mask = overflow;
    out
}

#[inline]
fn scores_decode_i8(value: i8, delta: i8) -> i32 {
    value as i32 - delta as i32
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn endpoint_i16_impl<
    const SEMI_GLOBAL: bool,
    const HAS_CBS: bool,
    const STANDARD_ONLY: bool,
>(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    scratch: &mut EndpointScratch,
) -> EndpointBatch {
    const LANES: usize = 16;
    let g = geometry::<LANES>(query.len(), targets);
    let delta = if SEMI_GLOBAL { 0 } else { i16::MIN };
    let zero = arch::_mm256_set1_epi16(delta);
    scratch.score16.resize(g.band, zero);
    scratch.hgap16.resize(g.band + 1, zero);
    scratch.score16.fill(zero);
    scratch.hgap16.fill(zero);
    let mut gap_open = [0i16; LANES];
    let mut gap_extend = [0i16; LANES];
    let mut overflow = 0u32;
    let mut standard_lanes = [0i16; LANES];
    for lane in 0..targets.len() {
        let scale = if targets[lane].matrix.is_some() {
            targets[lane].matrix_scale.max(1)
        } else {
            1
        };
        let go = (matrix.gap_open() + matrix.gap_extend()).saturating_mul(scale);
        let ge = matrix.gap_extend().saturating_mul(scale);
        if !(0..=i16::MAX as i32).contains(&go) || !(0..=i16::MAX as i32).contains(&ge) {
            overflow |= 1 << lane;
        }
        gap_open[lane] = go.clamp(0, i16::MAX as i32) as i16;
        gap_extend[lane] = ge.clamp(0, i16::MAX as i32) as i16;
        if targets[lane].matrix.is_none() {
            standard_lanes[lane] = -1;
        }
    }
    let go = arch::_mm256_loadu_si256(gap_open.as_ptr().cast());
    let ge = arch::_mm256_loadu_si256(gap_extend.as_ptr().cast());
    let standard_mask = arch::_mm256_loadu_si256(standard_lanes.as_ptr().cast());
    let mut sorted_offsets = g.band_offset;
    sorted_offsets[..targets.len()].sort_unstable();
    let mut part_bounds = [0usize; 18];
    let mut part_masks = [arch::_mm256_setzero_si256(); 17];
    let mut part_count = 0;
    let mut part_begin = 0;
    for &offset in &sorted_offsets[..targets.len()] {
        let offset = offset.min(g.band);
        if offset <= part_begin || offset == g.band {
            continue;
        }
        part_bounds[part_count] = part_begin;
        let mut lanes = [i16::MIN; LANES];
        for lane in 0..targets.len() {
            if part_begin >= g.band_offset[lane] {
                lanes[lane] = 0;
            }
        }
        part_masks[part_count] = arch::_mm256_loadu_si256(lanes.as_ptr().cast());
        part_count += 1;
        part_begin = offset;
    }
    part_bounds[part_count] = part_begin;
    let mut lanes = [i16::MIN; LANES];
    for lane in 0..targets.len() {
        if part_begin >= g.band_offset[lane] {
            lanes[lane] = 0;
        }
    }
    part_masks[part_count] = arch::_mm256_loadu_si256(lanes.as_ptr().cast());
    part_count += 1;
    part_bounds[part_count] = g.band;
    let mut best = [delta; LANES];
    let mut best_col = [0usize; LANES];
    let mut best_row = [0i32; LANES];
    for column in 0..g.columns {
        let mut subject = [0usize; LANES];
        let mut live_add = [i16::MIN; LANES];
        let mut live = [false; LANES];
        for lane in 0..targets.len() {
            let pos = g.subject_start[lane] + column as i32;
            if pos >= 0 && pos < targets[lane].subject.len() as i32 {
                subject[lane] = (targets[lane].subject[pos as usize] & LETTER_MASK) as usize;
                live_add[lane] = 0;
                live[lane] = true;
            }
        }
        let live_mask = arch::_mm256_loadu_si256(live_add.as_ptr().cast());
        let mut profile = [zero; 32];
        for (ql, slot) in profile[..AMINO_ACID_COUNT].iter_mut().enumerate() {
            let mut scores = [0i16; LANES];
            for lane in 0..targets.len() {
                if !live[lane] {
                    continue;
                }
                let sl = subject[lane];
                scores[lane] = if STANDARD_ONLY {
                    matrix.matrix16()[ql * 32 + sl]
                } else if let Some(adjusted) = targets[lane].matrix {
                    adjusted.scores[sl * 32 + ql] as i16
                } else {
                    matrix.matrix16()[ql * 32 + sl]
                };
            }
            *slot = arch::_mm256_loadu_si256(scores.as_ptr().cast());
        }
        let moving_i0 = g.i0 + column as i32;
        let query_begin = moving_i0.max(0);
        let query_end = (g.i1 + column as i32).min(query.len() as i32 - 1) + 1;
        if query_begin >= query_end {
            continue;
        }
        let active_begin = (query_begin - moving_i0) as usize;
        let active_end = (query_end - moving_i0) as usize;
        let mut vertical = zero;
        let mut col_best = zero;
        let initial_row = (delta as i32 + active_begin as i32) as i16;
        let mut row_counter = arch::_mm256_set1_epi16(initial_row);
        let mut row_max = zero;
        for part in 0..part_count {
            let row_begin = part_bounds[part].max(active_begin);
            let row_end = part_bounds[part + 1].min(active_end);
            if row_begin >= row_end {
                continue;
            }
            let cell_add = arch::_mm256_adds_epi16(part_masks[part], live_mask);
            vertical = arch::_mm256_adds_epi16(vertical, cell_add);
            for row in row_begin..row_end {
                let qpos = (moving_i0 + row as i32) as usize;
                let mut substitution = *profile.get_unchecked((query[qpos] & LETTER_MASK) as usize);
                substitution = arch::_mm256_adds_epi16(substitution, cell_add);
                if HAS_CBS {
                    let bias = arch::_mm256_and_si256(
                        arch::_mm256_set1_epi16(*cbs.get_unchecked(qpos) as i16),
                        standard_mask,
                    );
                    substitution = arch::_mm256_adds_epi16(substitution, bias);
                }
                let diagonal =
                    arch::_mm256_adds_epi16(*scratch.score16.get_unchecked(row), substitution);
                let horizontal =
                    arch::_mm256_adds_epi16(*scratch.hgap16.get_unchecked(row + 1), cell_add);
                let mut score = arch::_mm256_max_epi16(diagonal, horizontal);
                score = arch::_mm256_max_epi16(score, vertical);
                col_best = arch::_mm256_max_epi16(col_best, score);
                row_max = arch::_mm256_blendv_epi8(
                    row_max,
                    row_counter,
                    arch::_mm256_cmpeq_epi16(col_best, score),
                );
                row_counter = arch::_mm256_adds_epi16(row_counter, arch::_mm256_set1_epi16(1));
                let opened = arch::_mm256_subs_epi16(score, go);
                let next_horizontal =
                    arch::_mm256_max_epi16(arch::_mm256_subs_epi16(horizontal, ge), opened);
                vertical = arch::_mm256_max_epi16(arch::_mm256_subs_epi16(vertical, ge), opened);
                *scratch.score16.get_unchecked_mut(row) = score;
                *scratch.hgap16.get_unchecked_mut(row) = next_horizontal;
            }
        }
        let mut scores = [delta; LANES];
        let mut rows = [delta; LANES];
        arch::_mm256_storeu_si256(scores.as_mut_ptr().cast(), col_best);
        arch::_mm256_storeu_si256(rows.as_mut_ptr().cast(), row_max);
        for lane in 0..targets.len() {
            if live[lane] && scores[lane] > best[lane] {
                best[lane] = scores[lane];
                best_col[lane] = column;
                best_row[lane] = rows[lane] as i32 - delta as i32;
            }
        }
    }
    let mut out = EndpointBatch::empty(targets.len());
    for lane in 0..targets.len() {
        out.endpoints[lane] = Endpoint {
            score: best[lane] as i32 - delta as i32,
            query_end: g.i0 + best_col[lane] as i32 + best_row[lane] + 1,
            subject_end: g.subject_start[lane] + best_col[lane] as i32 + 1,
        };
        if best[lane] == i16::MAX {
            overflow |= 1 << lane;
        }
    }
    out.overflow_mask = overflow;
    out
}

pub fn endpoint_batch_i32(
    query: &[Letter],
    targets: &[TraceTarget<'_>],
    matrix: &ScoreMatrix,
    cbs: &[i8],
    scratch: &mut EndpointScratch,
) -> Option<EndpointBatch> {
    if targets.is_empty() || targets.len() > 32 || (!cbs.is_empty() && cbs.len() < query.len()) {
        return None;
    }
    let mut out = EndpointBatch::empty(targets.len());
    for (lane, &target) in targets.iter().enumerate() {
        out.endpoints[lane] = endpoint_i32(query, target, matrix, cbs, scratch)?;
    }
    Some(out)
}

fn endpoint_i32(
    query: &[Letter],
    target: TraceTarget<'_>,
    matrix: &ScoreMatrix,
    cbs: &[i8],
    scratch: &mut EndpointScratch,
) -> Option<Endpoint> {
    let scale = if target.matrix.is_some() {
        target.matrix_scale.max(1)
    } else {
        1
    };
    let open = matrix
        .gap_open()
        .checked_add(matrix.gap_extend())?
        .checked_mul(scale)?;
    let extend = matrix.gap_extend().checked_mul(scale)?;
    let neg = i32::MIN / 4;
    let band = target.d_end.saturating_sub(target.d_begin).max(0) as usize;
    scratch.i32_score.resize(band, 0);
    scratch.i32_hgap.resize(band + 1, neg);
    scratch.i32_score.fill(0);
    scratch.i32_hgap.fill(neg);
    let mut best = Endpoint::default();
    for (column, &sl) in target.subject.iter().enumerate() {
        let mut vertical = neg;
        let mut col_score = 0;
        let mut col_query_end = 0;
        for row in 0..band {
            let qpos = target.d_begin + column as i32 + row as i32;
            if qpos < 0 || qpos >= query.len() as i32 {
                vertical = neg;
                continue;
            }
            let qpos = qpos as usize;
            let ql = query[qpos];
            let subst = if let Some(adjusted) = target.matrix {
                adjusted.scores[(sl & LETTER_MASK) as usize * 32 + (ql & LETTER_MASK) as usize]
                    as i32
            } else {
                matrix.score(ql & LETTER_MASK, sl & LETTER_MASK)
                    + cbs.get(qpos).copied().unwrap_or(0) as i32
            };
            let score = (scratch.i32_score[row] + subst)
                .max(scratch.i32_hgap[row + 1])
                .max(vertical)
                .max(0);
            let opened = score - open;
            scratch.i32_hgap[row] = (scratch.i32_hgap[row + 1] - extend).max(opened);
            scratch.i32_score[row] = score;
            vertical = (vertical - extend).max(opened);
            if score >= col_score {
                col_score = score;
                col_query_end = qpos as i32 + 1;
            }
        }
        if col_score > best.score {
            best = Endpoint {
                score: col_score,
                query_end: col_query_end,
                subject_end: column as i32 + 1,
            };
        }
    }
    Some(best)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;
    use crate::stats::cbs::TargetMatrix;

    fn ordinary<'a>(subject: &'a [Letter], d_begin: i32, d_end: i32) -> TraceTarget<'a> {
        TraceTarget {
            subject,
            d_begin,
            d_end,
            matrix: None,
            matrix_scale: 1,
        }
    }

    fn oracle(
        query: &[Letter],
        target: TraceTarget<'_>,
        matrix: &ScoreMatrix,
        cbs: &[i8],
    ) -> Endpoint {
        let scale = if target.matrix.is_some() {
            target.matrix_scale.max(1)
        } else {
            1
        };
        let open = (matrix.gap_open() + matrix.gap_extend()) * scale;
        let extend = matrix.gap_extend() * scale;
        let neg = i32::MIN / 4;
        let mut previous = vec![0; query.len() + 1];
        let mut previous_gap = vec![neg; query.len() + 1];
        let mut best = Endpoint::default();
        for (j, &sl) in target.subject.iter().enumerate() {
            let mut current = vec![0; query.len() + 1];
            let mut current_gap = vec![neg; query.len() + 1];
            let mut vertical = neg;
            let mut col = Endpoint::default();
            for (q, &ql) in query.iter().enumerate() {
                if (q as i32) < target.d_begin + j as i32 || (q as i32) >= target.d_end + j as i32 {
                    vertical = neg;
                    continue;
                }
                let substitution = if let Some(adjusted) = target.matrix {
                    adjusted.scores[(sl & LETTER_MASK) as usize * 32 + (ql & LETTER_MASK) as usize]
                        as i32
                } else {
                    matrix.score(ql & LETTER_MASK, sl & LETTER_MASK)
                        + cbs.get(q).copied().unwrap_or(0) as i32
                };
                let score = (previous[q] + substitution)
                    .max(previous_gap[q + 1])
                    .max(vertical)
                    .max(0);
                current[q + 1] = score;
                current_gap[q + 1] = (previous_gap[q + 1] - extend).max(score - open);
                vertical = (vertical - extend).max(score - open);
                if score >= col.score {
                    col.score = score;
                    col.query_end = q as i32 + 1;
                    col.subject_end = j as i32 + 1;
                }
            }
            if col.score > best.score {
                best = col;
            }
            previous = current;
            previous_gap = current_gap;
        }
        best
    }

    #[test]
    fn randomized_true_i8_and_i16_match_independent_oracle() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        return;
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let mut state = 0x51de_cafe_u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as u32
        };
        let mut scratch = EndpointScratch::default();
        for _ in 0..60 {
            let qlen = 8 + next() as usize % 70;
            let query: Vec<_> = (0..qlen).map(|_| (next() % 25) as Letter).collect();
            let cbs: Vec<_> = (0..qlen).map(|_| (next() % 5) as i8 - 2).collect();
            let count = 2 + next() as usize % 15;
            let subjects: Vec<Vec<Letter>> = (0..count)
                .map(|_| {
                    (0..5 + next() as usize % 70)
                        .map(|_| {
                            let letter = (next() % 25) as Letter;
                            if next() % 13 == 0 {
                                letter | SEED_MASK
                            } else {
                                letter
                            }
                        })
                        .collect()
                })
                .collect();
            let bands: Vec<_> = (0..count)
                .map(|_| {
                    let begin = (next() % 21) as i32 - 10;
                    (begin, begin + 1 + (next() % 31) as i32)
                })
                .collect();
            let targets: Vec<_> = subjects
                .iter()
                .zip(&bands)
                .map(|(subject, &(begin, end))| ordinary(subject, begin, end))
                .collect();
            let byte = endpoint_batch_avx2_i8(&query, &targets, &matrix, &cbs, false, &mut scratch)
                .unwrap();
            let word =
                endpoint_batch_avx2_i16(&query, &targets, &matrix, &cbs, false, &mut scratch)
                    .unwrap();
            for lane in 0..count {
                let expected = oracle(&query, targets[lane], &matrix, &cbs);
                if expected.score == 0 {
                    assert_eq!(word.endpoints[lane].score, 0, "word lane {lane}");
                } else {
                    assert_eq!(word.endpoints[lane], expected, "word lane {lane}");
                }
                if byte.overflow_mask & (1 << lane) == 0 {
                    if expected.score == 0 {
                        assert_eq!(byte.endpoints[lane].score, 0, "byte lane {lane}");
                    } else {
                        assert_eq!(byte.endpoints[lane], expected, "byte lane {lane}");
                    }
                }
            }
        }
    }

    #[test]
    fn mixed_adjusted_lanes_keep_profiles_gaps_and_cbs_separate() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        return;
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let adjusted = TargetMatrix::new(
            (0..1024)
                .map(|i| if i / 32 == i % 32 { 3 } else { -4 })
                .collect(),
            -4,
            3,
        );
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let a = query.clone();
        let b: Vec<_> = query.iter().rev().copied().collect();
        let targets = [
            ordinary(&a, -2, 4),
            TraceTarget {
                subject: &b,
                d_begin: 0,
                d_end: 3,
                matrix: Some(&adjusted),
                matrix_scale: 2,
            },
        ];
        let cbs = vec![7; query.len()];
        let mut scratch = EndpointScratch::default();
        for batch in [
            endpoint_batch_avx2_i8(&query, &targets, &matrix, &cbs, false, &mut scratch).unwrap(),
            endpoint_batch_avx2_i16(&query, &targets, &matrix, &cbs, false, &mut scratch).unwrap(),
        ] {
            assert_eq!(batch.overflow_mask, 0);
            for lane in 0..targets.len() {
                assert_eq!(
                    batch.endpoints[lane],
                    oracle(&query, targets[lane], &matrix, &cbs)
                );
            }
        }
    }

    #[test]
    fn shifted_and_unshifted_tiers_promote_at_upstream_limits() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        return;
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 100_000).unwrap();
        let q8 = vec![0; 40];
        let t8 = ordinary(&q8, 0, 1);
        let mut scratch = EndpointScratch::default();
        let semi8 = endpoint_batch_avx2_i8(&q8, &[t8], &matrix, &[], true, &mut scratch).unwrap();
        let local8 = endpoint_batch_avx2_i8(&q8, &[t8], &matrix, &[], false, &mut scratch).unwrap();
        assert_eq!(semi8.overflow_mask, 1);
        assert_eq!(local8.overflow_mask, 0);
        assert_eq!(local8.endpoints[0], oracle(&q8, t8, &matrix, &[]));

        let q16 = vec![0; 256];
        let bias = vec![i8::MAX; q16.len()];
        let t16 = ordinary(&q16, 0, 1);
        let semi16 =
            endpoint_batch_avx2_i16(&q16, &[t16], &matrix, &bias, true, &mut scratch).unwrap();
        let local16 =
            endpoint_batch_avx2_i16(&q16, &[t16], &matrix, &bias, false, &mut scratch).unwrap();
        assert_eq!(semi16.overflow_mask, 1);
        assert_eq!(local16.overflow_mask, 0);
        assert_eq!(local16.endpoints[0], oracle(&q16, t16, &matrix, &bias));
    }

    #[test]
    fn semi_global_does_not_add_a_per_cell_zero_clamp() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        return;
        // Pinned upstream counterexample: DELTA=0's saturating recurrence
        // selects the second query residue, whereas an added h=max(h,0)
        // selects the third.
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query = [0, 0, 0];
        let subject = [1, 1, 0];
        let target = ordinary(&subject, -2, 2);
        let mut scratch = EndpointScratch::default();
        let byte =
            endpoint_batch_avx2_i8(&query, &[target], &matrix, &[], true, &mut scratch).unwrap();
        let word =
            endpoint_batch_avx2_i16(&query, &[target], &matrix, &[], true, &mut scratch).unwrap();
        assert_eq!(
            byte.endpoints[0],
            Endpoint {
                score: 4,
                query_end: 2,
                subject_end: 3
            }
        );
        assert_eq!(word.endpoints[0], byte.endpoints[0]);
    }

    #[test]
    fn shifted_i8_row_counter_wraps_clipped_offset_above_127() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        return;
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query = vec![0, 0, 0, 0];
        let subject = vec![0, 0, 0, 0];
        let target = ordinary(&subject, -200, 1);
        let got = endpoint_batch_avx2_i8(
            &query,
            &[target],
            &matrix,
            &[],
            false,
            &mut EndpointScratch::default(),
        )
        .unwrap();
        assert_eq!(got.overflow_mask, 0);
        assert_eq!(got.endpoints[0], oracle(&query, target, &matrix, &[]));
    }

    #[test]
    fn full_matrix_kernels_use_query_rows_and_match_oracle() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        return;
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 10_000).unwrap();
        let query: Vec<_> = (0..23).map(|i| (i * 7 % 20) as Letter).collect();
        let subjects: Vec<Vec<Letter>> = (0..32)
            .map(|lane| {
                (0..7 + lane % 13)
                    .map(|i| ((i * 11 + lane) % 20) as Letter)
                    .collect()
            })
            .collect();
        let targets: Vec<_> = subjects
            .iter()
            .map(|subject| ordinary(subject, -(subject.len() as i32 - 1), query.len() as i32))
            .collect();
        let mut scratch = EndpointScratch::default();
        let byte =
            endpoint_full_batch_avx2_i8(&query, &targets, &matrix, &[], false, 1, &mut scratch)
                .unwrap();
        assert_eq!(byte.len, 32);
        for lane in 0..32 {
            assert_eq!(
                byte.endpoints[lane],
                oracle(&query, targets[lane], &matrix, &[])
            );
        }
        let word = endpoint_full_batch_avx2_i16(
            &query,
            &targets[..16],
            &matrix,
            &[],
            false,
            1,
            &mut scratch,
        )
        .unwrap();
        for lane in 0..16 {
            assert_eq!(
                word.endpoints[lane],
                oracle(&query, targets[lane], &matrix, &[])
            );
        }
    }
}
