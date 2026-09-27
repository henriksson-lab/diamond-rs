//! Scalar SWIPE wrapper pieces.
//!
//! This module ports the bottom-up, non-SIMD parts of
//! `diamond/src/dp/dp.h` and `diamond/src/dp/swipe/swipe_wrapper.cpp`: target
//! containers, value flags, bin selection, target sorting, and the sequential
//! banded swipe flow. The C++ implementation dispatches to score-vector SIMD
//! kernels; this Rust path evaluates each target with the same diagonal band
//! semantics in scalar code.

use crate::align::hsp::Hsp;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::translate::{Frame, TranslatedPosition};
use crate::basic::value::{Letter, LETTER_MASK, SEED_MASK};
use crate::data::sequence_set::SequenceSet;
use crate::dp::smith_waterman::SwResult;
use crate::stats;
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::geo;
use crate::util::interval::Interval;
use std::sync::{Arc, Mutex};

#[path = "swipe/banded_3frame_swipe.rs"]
pub mod banded_3frame_swipe;
use std::time::Instant;

#[path = "swipe/anchored_wrapper.rs"]
pub mod anchored_wrapper;

#[path = "swipe/swipe_wrapper.rs"]
pub mod swipe_wrapper;

#[path = "swipe/simd_score.rs"]
pub mod simd_score;

#[path = "swipe/simd_score8.rs"]
pub mod simd_score8;

#[path = "swipe/simd_score8_portable.rs"]
pub mod simd_score8_portable;

#[path = "swipe/simd_score_portable.rs"]
pub mod simd_score_portable;

#[path = "swipe/simd_adjusted_narrow.rs"]
mod simd_adjusted_narrow;

#[path = "swipe/simd_trace.rs"]
pub mod simd_trace;

#[path = "swipe/simd_trace_narrow.rs"]
mod simd_trace_narrow;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod simd_trace_upstream;

#[path = "swipe/simd_trace_narrow_portable.rs"]
mod simd_trace_narrow_portable;

#[path = "swipe/simd_trace_portable.rs"]
pub mod simd_trace_portable;

pub use swipe_wrapper::{
    DispatchCell, DispatchConfig, IdMaskKind, RowCounterKind, SwipeRuntimeConfig,
};

pub const BINS: usize = 6;
pub const SCORE_BINS: usize = 3;
pub const ALGO_BINS: usize = 2;
pub const BLANK_TARGET: i64 = i64::MAX;
pub const MIN_LETTERS: i32 = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Flags(pub u8);

impl Flags {
    pub const NONE: Flags = Flags(0);
    pub const PARALLEL: Flags = Flags(1);
    pub const FULL_MATRIX: Flags = Flags(2);
    pub const SEMI_GLOBAL: Flags = Flags(4);

    pub fn any(self, rhs: Flags) -> bool {
        self.0 & rhs.0 != 0
    }
}

impl std::ops::BitOr for Flags {
    type Output = Flags;

    fn bitor(self, rhs: Flags) -> Flags {
        Flags(self.0 | rhs.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct HspValues(pub u32);

impl HspValues {
    pub const NONE: HspValues = HspValues(0);
    pub const TRANSCRIPT: HspValues = HspValues(1);
    pub const QUERY_START: HspValues = HspValues(1 << 1);
    pub const QUERY_END: HspValues = HspValues(1 << 2);
    pub const TARGET_START: HspValues = HspValues(1 << 3);
    pub const TARGET_END: HspValues = HspValues(1 << 4);
    pub const IDENT: HspValues = HspValues(1 << 5);
    pub const LENGTH: HspValues = HspValues(1 << 6);
    pub const MISMATCHES: HspValues = HspValues(1 << 7);
    pub const GAP_OPENINGS: HspValues = HspValues(1 << 8);
    pub const GAPS: HspValues = HspValues(Self::IDENT.0 | Self::LENGTH.0 | Self::MISMATCHES.0);
    pub const QUERY_COORDS: HspValues = HspValues(Self::QUERY_START.0 | Self::QUERY_END.0);
    pub const TARGET_COORDS: HspValues = HspValues(Self::TARGET_START.0 | Self::TARGET_END.0);
    pub const COORDS: HspValues = HspValues(Self::QUERY_COORDS.0 | Self::TARGET_COORDS.0);

    pub fn any(self, rhs: HspValues) -> bool {
        self.0 & rhs.0 != 0
    }

    pub fn only(self, rhs: HspValues) -> bool {
        self.0 & !rhs.0 == 0
    }

    pub fn all(self, rhs: HspValues) -> bool {
        self.0 & rhs.0 == rhs.0
    }
}

impl std::ops::BitOr for HspValues {
    type Output = HspValues;

    fn bitor(self, rhs: HspValues) -> HspValues {
        HspValues(self.0 | rhs.0)
    }
}

pub fn have_coords(v: HspValues) -> bool {
    v.any(HspValues::TRANSCRIPT) || v.all(HspValues::QUERY_COORDS | HspValues::TARGET_COORDS)
}

const NO_TRACEBACK: HspValues = HspValues(
    HspValues::COORDS.0
        | HspValues::IDENT.0
        | HspValues::LENGTH.0
        | HspValues::MISMATCHES.0
        | HspValues::GAP_OPENINGS.0,
);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CarryOver {
    pub i1: i32,
    pub j1: i32,
    pub ident: i32,
    pub len: i32,
}

impl CarryOver {
    pub fn new(i1: i32, j1: i32, ident: i32, len: i32) -> Self {
        CarryOver { i1, j1, ident, len }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Anchor {
    pub query_begin: i32,
    pub query_end: i32,
    pub subject_begin: i32,
    pub subject_end: i32,
    pub d_min_left: i32,
    pub d_max_left: i32,
    pub d_min_right: i32,
    pub d_max_right: i32,
    pub prefix_score: i32,
    pub score: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DpTarget {
    /// Immutable target sequence shared across bins and alignment rounds.
    pub seq: Arc<[Letter]>,
    pub d_begin: i32,
    pub d_end: i32,
    pub cols: i32,
    pub true_target_len: i32,
    pub target_idx: i64,
    pub carry_over: CarryOver,
    pub anchor: Anchor,
    pub matrix: Option<Arc<TargetMatrix>>,
    pub matrix_scale: i32,
}

impl Default for DpTarget {
    fn default() -> Self {
        DpTarget {
            seq: Arc::from([]),
            d_begin: 0,
            d_end: 0,
            cols: 0,
            true_target_len: 0,
            target_idx: BLANK_TARGET,
            carry_over: CarryOver::default(),
            anchor: Anchor::default(),
            matrix: None,
            matrix_scale: 1,
        }
    }
}

impl DpTarget {
    pub fn banded_cols(qlen: i32, tlen: i32, d_begin: i32, d_end: i32) -> i32 {
        let pos = (d_end - 1).max(0) - (d_end - 1);
        let d0 = d_begin;
        let j1 = (qlen - 1 - d0).min(tlen - 1) + 1;
        j1 - pos
    }

    pub fn new<S>(
        seq: S,
        true_target_len: i32,
        d_begin: i32,
        d_end: i32,
        target_idx: i64,
        qlen: i32,
        carry_over: CarryOver,
        anchor: Anchor,
    ) -> Self
    where
        S: Into<Arc<[Letter]>>,
    {
        let seq = seq.into();
        let cols = Self::banded_cols(qlen, seq.len() as i32, d_begin, d_end);
        DpTarget {
            true_target_len,
            d_begin,
            d_end,
            cols,
            target_idx,
            carry_over,
            anchor,
            matrix: None,
            matrix_scale: 1,
            seq,
        }
    }

    pub fn full<S>(seq: S, true_target_len: i32, target_idx: i64, carry_over: CarryOver) -> Self
    where
        S: Into<Arc<[Letter]>>,
    {
        DpTarget {
            true_target_len,
            target_idx,
            carry_over,
            seq: seq.into(),
            ..Default::default()
        }
    }

    pub fn with_matrix(mut self, matrix: Arc<TargetMatrix>, matrix_scale: i32) -> Self {
        self.matrix = Some(matrix);
        self.matrix_scale = matrix_scale;
        self
    }

    pub fn left_i1(&self) -> i32 {
        (self.d_end - 1).max(0)
    }

    pub fn band(&self) -> i32 {
        self.d_end - self.d_begin
    }

    pub fn blank(&self) -> bool {
        self.target_idx == BLANK_TARGET
    }

    pub fn adjusted_matrix(&self) -> bool {
        self.matrix.is_some()
    }

    pub fn matrix_scale(&self) -> i32 {
        if self.adjusted_matrix() {
            self.matrix_scale
        } else {
            1
        }
    }

    pub fn cells(&self, flags: Flags, qlen: i32) -> i64 {
        if flags.any(Flags::FULL_MATRIX) {
            self.seq.len() as i64 * qlen as i64
        } else {
            self.band() as i64 * self.cols as i64
        }
    }

    pub fn extend_right(&self, qlen: i32) -> bool {
        (qlen - self.anchor.query_end).min(self.seq.len() as i32 - self.anchor.subject_end)
            >= MIN_LETTERS
    }

    pub fn extend_left(&self) -> bool {
        self.anchor.query_begin.min(self.anchor.subject_begin) >= MIN_LETTERS
    }
}

#[derive(Debug, Clone, Default)]
pub struct TargetVec {
    targets: Vec<DpTarget>,
    max_len: i32,
}

impl TargetVec {
    pub fn begin(&self) -> std::slice::Iter<'_, DpTarget> {
        self.targets.iter()
    }

    pub fn end(&self) -> std::slice::Iter<'_, DpTarget> {
        self.targets[self.targets.len()..].iter()
    }

    pub fn as_slice(&self) -> &[DpTarget] {
        &self.targets
    }

    pub fn as_mut_slice(&mut self) -> &mut [DpTarget] {
        &mut self.targets
    }

    pub fn front(&self) -> Option<&DpTarget> {
        self.targets.first()
    }

    pub fn back(&self) -> Option<&DpTarget> {
        self.targets.last()
    }

    pub fn size(&self) -> i64 {
        self.targets.len() as i64
    }

    pub fn reserve(&mut self, size: i64) {
        self.targets.reserve(size.max(0) as usize);
    }

    pub fn push_back(&mut self, target: DpTarget) {
        self.max_len = self.max_len.max(target.seq.len() as i32);
        self.targets.push(target);
    }

    pub fn push_back_vec(&mut self, targets: &TargetVec) {
        self.max_len = self.max_len.max(targets.max_len);
        self.targets.extend_from_slice(&targets.targets);
    }

    pub fn empty(&self) -> bool {
        self.targets.is_empty()
    }

    pub fn clear(&mut self) {
        self.targets.clear();
        self.max_len = 0;
    }

    pub fn max_len(&self) -> i32 {
        self.max_len
    }
}

impl std::ops::Index<usize> for TargetVec {
    type Output = DpTarget;

    fn index(&self, index: usize) -> &Self::Output {
        &self.targets[index]
    }
}

impl std::ops::IndexMut<usize> for TargetVec {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.targets[index]
    }
}

pub type Targets = [TargetVec; BINS];

pub fn targets() -> Targets {
    std::array::from_fn(|_| TargetVec::default())
}

#[derive(Clone)]
pub struct Params<'a> {
    pub query: &'a [Letter],
    pub query_id: Option<&'a str>,
    pub frame: i32,
    pub query_source_len: i32,
    pub composition_bias: Option<&'a [i8]>,
    pub flags: Flags,
    pub reverse_targets: bool,
    pub target_max_len: i32,
    pub swipe_bin: i32,
    pub v: HspValues,
    pub score_matrix: &'a ScoreMatrix,
    pub band_bin: i32,
    pub col_bin: i32,
    pub cutoff_score_8bit: i32,
    pub max_swipe_dp: i64,
    pub approx_backtrace: bool,
    pub max_evalue: f64,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub approx_min_id: f64,
    pub cbs_matrix_scale: i32,
    /// Shared sink mirroring C++ `DP::Params::stat`. `Arc<Mutex<_>>` keeps
    /// cloned reverse-pass parameters on the same accumulator.
    pub statistics: Option<Arc<Mutex<Statistics>>>,
}

impl<'a> Params<'a> {
    pub fn new(query: &'a [Letter], score_matrix: &'a ScoreMatrix) -> Self {
        Params {
            query,
            query_id: None,
            frame: 0,
            query_source_len: query.len() as i32,
            composition_bias: None,
            flags: Flags::NONE,
            reverse_targets: false,
            target_max_len: 0,
            swipe_bin: 0,
            v: HspValues::COORDS,
            score_matrix,
            band_bin: 16,
            col_bin: 16,
            cutoff_score_8bit: i8::MAX as i32,
            max_swipe_dp: 1_000_000,
            approx_backtrace: false,
            max_evalue: f64::MAX,
            query_cover: 0.0,
            subject_cover: 0.0,
            query_or_target_cover: 0.0,
            approx_min_id: 0.0,
            cbs_matrix_scale: 1,
            statistics: None,
        }
    }

    pub fn inc_stat(&self, value: StatValue, count: i64) {
        if let Some(statistics) = &self.statistics {
            statistics.lock().unwrap().inc(value, count);
        }
    }
}

pub fn sort(targets: &mut [DpTarget], band_bin: i32, col_bin: i32) {
    targets.sort_by(|a, b| {
        let i = a.left_i1();
        let j = b.left_i1();
        let b1 = a.band();
        let b2 = b.band();
        let bin_b1 = b1 / band_bin;
        let bin_b2 = b2 / band_bin;
        let t1 = a.cols;
        let t2 = b.cols;
        let bin_t1 = t1 / col_bin;
        let bin_t2 = t2 / col_bin;
        bin_b1
            .cmp(&bin_b2)
            .then_with(|| bin_t1.cmp(&bin_t2))
            .then_with(|| i.cmp(&j))
    });
}

pub fn bin_score(x: i32) -> usize {
    if x < u8::MAX as i32 {
        0
    } else if x < u16::MAX as i32 {
        1
    } else {
        2
    }
}

pub fn bin(
    v: HspValues,
    query_len: i32,
    score: i32,
    ungapped_score: i32,
    dp_size: i64,
    score_width: usize,
    mismatch_est: i32,
    cutoff_score_8bit: i32,
    max_swipe_dp: i64,
    approx_backtrace: bool,
) -> usize {
    let mut b = 0usize;
    b = b.max(bin_score(score));
    if ungapped_score > cutoff_score_8bit {
        b = b.max(1);
    }
    b = b.max(score_width);
    b = b.max(bin_score(mismatch_est));
    if v != HspValues::NONE {
        b = b.max(bin_score(query_len));
        if dp_size > max_swipe_dp {
            if v.only(NO_TRACEBACK) {
                b += SCORE_BINS;
            } else {
                b = 2;
            }
        } else if v.only(HspValues::COORDS) && !approx_backtrace {
            b += SCORE_BINS;
        }
    }
    b
}

pub fn matrix_size(query_len: i32, targets: &[DpTarget], flags: Flags, channels: i64) -> i64 {
    let mut s = 0i64;
    for target in targets {
        let cols = if flags.any(Flags::FULL_MATRIX) {
            target.seq.len() as i32
        } else {
            target.cols
        };
        let rows = if flags.any(Flags::FULL_MATRIX) {
            query_len
        } else {
            target.d_end - target.d_begin
        };
        let size = rows as i64 * cols as i64 * channels / 2;
        s = s.max(size);
    }
    s
}

pub fn reversed(v: HspValues) -> bool {
    v.only(NO_TRACEBACK)
        && v.any(
            HspValues::QUERY_START
                | HspValues::TARGET_START
                | HspValues::MISMATCHES
                | HspValues::GAP_OPENINGS,
        )
}

pub fn mismatch_est(query_len: i32, target_len: i32, aln_len: i32, v: HspValues) -> i32 {
    if !v.any(HspValues::MISMATCHES) {
        return 0;
    }
    let m = query_len.min(target_len);
    if aln_len > 0 {
        aln_len.min(m)
    } else {
        m
    }
}

fn score_only_hsp(
    target: &DpTarget,
    sw_score: i32,
    d_begin: i32,
    d_end: i32,
    p: &Params<'_>,
) -> Option<Hsp> {
    if sw_score <= 0 {
        return None;
    }
    let score = if target.adjusted_matrix() {
        sw_score
    } else {
        sw_score * p.cbs_matrix_scale
    };
    let evalue = p
        .score_matrix
        .evalue(score, p.query.len() as u32, target.true_target_len as u32);
    if evalue > p.max_evalue {
        return None;
    }
    let mut hsp = Hsp::new();
    hsp.score = score;
    hsp.bit_score = p.score_matrix.bitscore(score as f64);
    hsp.corrected_bit_score = p.score_matrix.bitscore_corrected(
        score,
        p.query.len() as u32,
        target.true_target_len as u32,
    );
    hsp.evalue = evalue;
    hsp.frame = p.frame;
    hsp.d_begin = d_begin;
    hsp.d_end = d_end;
    hsp.swipe_target = target.target_idx as i32;
    hsp.swipe_bin = p.swipe_bin;
    Some(hsp)
}

fn traceback_hsp(
    target: &DpTarget,
    target_seq: &[Letter],
    d_begin: i32,
    d_end: i32,
    sw: SwResult,
    p: &Params<'_>,
) -> Option<Hsp> {
    if sw.score <= 0 {
        return None;
    }
    let score = if target.adjusted_matrix() {
        sw.score
    } else {
        sw.score * p.cbs_matrix_scale
    };
    let evalue = p
        .score_matrix
        .evalue(score, p.query.len() as u32, target.true_target_len as u32);
    if evalue > p.max_evalue {
        return None;
    }
    let mut hsp = Hsp::new();
    hsp.backtraced = true;
    hsp.score = score;
    hsp.bit_score = p.score_matrix.bitscore(score as f64);
    hsp.corrected_bit_score = p.score_matrix.bitscore_corrected(
        score,
        p.query.len() as u32,
        target.true_target_len as u32,
    );
    hsp.evalue = evalue;
    hsp.frame = p.frame;
    if target.carry_over.i1 == 0 {
        hsp.length = sw.length;
        hsp.identities = sw.identities;
    } else {
        hsp.length = target.carry_over.len;
        hsp.identities = target.carry_over.ident;
    }
    hsp.mismatches = sw.mismatches;
    hsp.gap_openings = sw.gap_openings;
    hsp.gaps = hsp.length - hsp.identities - hsp.mismatches;
    if target.carry_over.i1 == 0 {
        hsp.query_range = Interval::new(sw.query_begin, sw.query_end);
        hsp.subject_range = Interval::new(sw.subject_begin, sw.subject_end);
        hsp.d_begin = d_begin;
        hsp.d_end = d_end;
    } else {
        let qlen = p.query.len() as i32;
        let tlen = target.seq.len() as i32;
        hsp.query_range = Interval::new(qlen - sw.query_end, target.carry_over.i1);
        hsp.subject_range = Interval::new(tlen - sw.subject_end, target.carry_over.j1);
        hsp.d_begin = -target.d_end + qlen - tlen + 1;
        hsp.d_end = -target.d_begin + qlen - tlen + 1;
    }
    hsp.query_source_range = TranslatedPosition::absolute_interval(
        TranslatedPosition::new(hsp.query_range.begin, Frame::from_index(p.frame)),
        TranslatedPosition::new(hsp.query_range.end, Frame::from_index(p.frame)),
        p.query_source_len,
        true,
    );
    hsp.subject_source_range = hsp.subject_range;
    hsp.target_seq = target.seq.clone();
    hsp.matrix = target.matrix.clone();
    hsp.swipe_target = target.target_idx as i32;
    hsp.swipe_bin = p.swipe_bin;
    let mut qi = sw.query_begin as usize;
    let mut sj = sw.subject_begin as usize;
    for (op, len) in sw.operations {
        match op {
            EditOperation::Match => {
                hsp.transcript
                    .push_with_count(EditOperation::Match, len as u32);
                hsp.positives += len;
                qi += len as usize;
                sj += len as usize;
            }
            EditOperation::Substitution => {
                for _ in 0..len {
                    let ql = p.query[qi];
                    let sl = target_seq[sj];
                    let match_score = if sl & SEED_MASK != 0 {
                        0
                    } else if let Some(matrix) = target.matrix.as_deref() {
                        matrix.scores
                            [(sl & LETTER_MASK) as usize * 32 + (ql & LETTER_MASK) as usize]
                            as i32
                    } else {
                        p.score_matrix.score(ql & LETTER_MASK, sl & LETTER_MASK)
                    };
                    if match_score > 0 {
                        hsp.positives += 1;
                    }
                    hsp.transcript
                        .push_with_letter(EditOperation::Substitution, sl);
                    qi += 1;
                    sj += 1;
                }
            }
            EditOperation::Insertion => {
                hsp.transcript
                    .push_with_count(EditOperation::Insertion, len as u32);
                qi += len as usize;
            }
            EditOperation::Deletion => {
                for _ in 0..len {
                    hsp.transcript
                        .push_with_letter(EditOperation::Deletion, target_seq[sj]);
                    sj += 1;
                }
            }
            EditOperation::FrameshiftForward | EditOperation::FrameshiftReverse => {
                for _ in 0..len {
                    hsp.transcript.push(op);
                }
            }
        }
    }
    hsp.transcript.push_terminator();
    hsp.approx_id = hsp.approx_id_percent(p.query, target_seq);
    Some(hsp)
}

pub fn dispatch_swipe(
    subject_begin: &[DpTarget],
    overflow: &mut TargetVec,
    p: &Params<'_>,
) -> Vec<Hsp> {
    let mut out = Vec::new();
    let bias = p.composition_bias.unwrap_or(&[]);
    let mut score_scratch = ScoreScratch::default();
    let mut traceback_scratch = TracebackScratch::default();

    // Target-specific composition matrices use per-lane profiles and scaled
    // gap vectors. Keep them in the same byte -> word -> i32 promotion cascade
    // as upstream, including batches which mix adjusted and ordinary lanes.
    if p.v == HspValues::NONE
        && !p.reverse_targets
        && subject_begin.iter().any(DpTarget::adjusted_matrix)
    {
        let score_bin = p.swipe_bin.max(0) as usize % SCORE_BINS;
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        let avx2 = std::arch::is_x86_feature_detected!("avx2");
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        let avx2 = false;
        let lane_width = match (avx2, score_bin) {
            (true, 0) => 32,
            (true, 1) => 16,
            (false, 0) => 16,
            (false, 1) => 8,
            _ => 8,
        };
        for chunk in subject_begin.chunks(lane_width) {
            let lane_targets: Vec<&DpTarget> = chunk
                .iter()
                .filter(|target| !target.blank() && !target.seq.is_empty())
                .collect();
            let score_targets: Vec<_> = lane_targets
                .iter()
                .map(|target| {
                    let (d_begin, d_end) = if p.flags.any(Flags::FULL_MATRIX) {
                        (-(target.seq.len() as i32 - 1), p.query.len() as i32)
                    } else {
                        (target.d_begin, target.d_end)
                    };
                    simd_trace::TraceTarget {
                        subject: target.seq.as_ref(),
                        d_begin,
                        d_end,
                        matrix: target.matrix.as_deref(),
                        matrix_scale: target.matrix_scale(),
                    }
                })
                .collect();
            let (scores, overflow_mask) = match score_bin {
                0 => simd_adjusted_narrow::score_batch_avx2_i8(
                    p.query,
                    &score_targets,
                    p.score_matrix,
                    bias,
                )
                .map(|x| (Some(x.scores), x.overflow_mask))
                .or_else(|| {
                    simd_trace_narrow_portable::score_batch_i8(
                        p.query,
                        &score_targets,
                        p.score_matrix,
                        bias,
                    )
                    .map(|x| (Some(x.scores), x.overflow_mask))
                })
                .unwrap_or((None, 0)),
                1 => simd_adjusted_narrow::score_batch_avx2_i16(
                    p.query,
                    &score_targets,
                    p.score_matrix,
                    bias,
                )
                .map(|x| (Some(x.scores), x.overflow_mask))
                .or_else(|| {
                    simd_trace_narrow_portable::score_batch_i16(
                        p.query,
                        &score_targets,
                        p.score_matrix,
                        bias,
                    )
                    .map(|x| (Some(x.scores), x.overflow_mask))
                })
                .unwrap_or((None, 0)),
                _ => (
                    simd_trace::score_adjusted_batch_avx2(
                        p.query,
                        &score_targets,
                        p.score_matrix,
                        bias,
                    )
                    .or_else(|| {
                        simd_trace_portable::score_adjusted_batch_portable(
                            p.query,
                            &score_targets,
                            p.score_matrix,
                            bias,
                        )
                    }),
                    0,
                ),
            };
            for (lane, (&target, score_target)) in
                lane_targets.iter().zip(&score_targets).enumerate()
            {
                if overflow_mask & (1 << lane) != 0 {
                    overflow.push_back(target.clone());
                    continue;
                }
                let score = scores
                    .as_ref()
                    .map(|scores| scores[lane])
                    .unwrap_or_else(|| {
                        banded_sw_cbs_score(
                            p.query,
                            score_target.subject,
                            score_target.d_begin,
                            score_target.d_end,
                            p.score_matrix,
                            bias,
                            target.matrix.as_deref(),
                            target.matrix_scale(),
                            &mut score_scratch,
                        )
                    });
                if let Some(hsp) =
                    score_only_hsp(target, score, score_target.d_begin, score_target.d_end, p)
                {
                    out.push(hsp);
                }
            }
        }
        return out;
    }

    // C++ bin 0/3 uses 32 signed-byte lanes and promotes saturated lanes to
    // the next score-width bin. Keep that flow instead of paying for i16 on
    // every low-scoring target.
    if p.v == HspValues::NONE
        && !p.reverse_targets
        && (simd_score8::available() || simd_score8_portable::available())
        && p.swipe_bin >= 0
        && p.swipe_bin as usize % SCORE_BINS == 0
        && subject_begin
            .iter()
            .all(|target| target.blank() || target.seq.is_empty() || !target.adjusted_matrix())
    {
        let mut scratch = simd_score8::Scratch8::default();
        let mut portable_scratch = simd_score8_portable::PortableScratch8::default();
        let lane_width = if simd_score8::available() { 32 } else { 16 };
        for chunk in subject_begin.chunks(lane_width) {
            let mut lanes = [simd_score::ScoreTarget {
                subject: &[],
                d_begin: 0,
                d_end: 0,
            }; 32];
            let mut lane_targets: [&DpTarget; 32] = [&subject_begin[0]; 32];
            let mut count = 0;
            for target in chunk {
                if target.blank() || target.seq.is_empty() {
                    continue;
                }
                lanes[count] = simd_score::ScoreTarget {
                    subject: target.seq.as_ref(),
                    d_begin: if p.flags.any(Flags::FULL_MATRIX) {
                        -(target.seq.len() as i32 - 1)
                    } else {
                        target.d_begin
                    },
                    d_end: if p.flags.any(Flags::FULL_MATRIX) {
                        p.query.len() as i32
                    } else {
                        target.d_end
                    },
                };
                lane_targets[count] = target;
                count += 1;
            }
            if count == 0 {
                continue;
            }
            let batch = if p.flags.any(Flags::FULL_MATRIX) {
                let mut subjects = [&[][..]; 32];
                for lane in 0..count {
                    subjects[lane] = lanes[lane].subject;
                }
                simd_score8::score_full_batch_avx2_i8(
                    p.query,
                    &subjects[..count],
                    p.score_matrix,
                    bias,
                    p.flags.any(Flags::SEMI_GLOBAL),
                    &mut scratch,
                )
                .or_else(|| {
                    simd_score8_portable::score_full_batch_portable_i8(
                        p.query,
                        &subjects[..count],
                        p.score_matrix,
                        bias,
                        p.flags.any(Flags::SEMI_GLOBAL),
                        &mut portable_scratch,
                    )
                })
            } else {
                simd_score8::score_batch_avx2_i8(
                    p.query,
                    &lanes[..count],
                    p.score_matrix,
                    bias,
                    p.flags.any(Flags::SEMI_GLOBAL),
                    &mut scratch,
                )
                .or_else(|| {
                    simd_score8_portable::score_batch_portable_i8(
                        p.query,
                        &lanes[..count],
                        p.score_matrix,
                        bias,
                        p.flags.any(Flags::SEMI_GLOBAL),
                        &mut portable_scratch,
                    )
                })
            };
            let Some(batch) = batch else {
                for &target in &lane_targets[..count] {
                    let (d_begin, d_end) = if p.flags.any(Flags::FULL_MATRIX) {
                        (-(target.seq.len() as i32 - 1), p.query.len() as i32)
                    } else {
                        (target.d_begin, target.d_end)
                    };
                    let score = banded_sw_cbs_score(
                        p.query,
                        target.seq.as_ref(),
                        d_begin,
                        d_end,
                        p.score_matrix,
                        bias,
                        None,
                        1,
                        &mut score_scratch,
                    );
                    if let Some(hsp) = score_only_hsp(target, score, d_begin, d_end, p) {
                        out.push(hsp);
                    }
                }
                continue;
            };
            for (lane, &target) in lane_targets[..count].iter().enumerate() {
                if batch.overflow_mask & (1 << lane) != 0 {
                    overflow.push_back(target.clone());
                } else if let Some(hsp) = score_only_hsp(
                    target,
                    batch.scores[lane],
                    lanes[lane].d_begin,
                    lanes[lane].d_end,
                    p,
                ) {
                    out.push(hsp);
                }
            }
        }
        return out;
    }

    // The normal first extension round is score-only. Batch ordinary forward,
    // banded targets into AVX2 lanes, retaining the scalar i32 calculation for
    // saturation and every unsupported semantic combination.
    if p.v == HspValues::NONE
        && !p.reverse_targets
        && subject_begin
            .iter()
            .all(|target| target.blank() || target.seq.is_empty() || !target.adjusted_matrix())
    {
        let mut simd_scratch = simd_score::SimdScoreScratch::default();
        let mut portable_scratch = simd_score_portable::PortableSimdScoreScratch::default();
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        let lane_width = if std::arch::is_x86_feature_detected!("avx2") {
            16
        } else {
            8
        };
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        let lane_width = 8;
        for chunk in subject_begin.chunks(lane_width) {
            let mut lanes = [simd_score::ScoreTarget {
                subject: &[],
                d_begin: 0,
                d_end: 0,
            }; 16];
            let mut lane_targets: [&DpTarget; 16] = [&subject_begin[0]; 16];
            let mut lane_count = 0usize;
            for target in chunk {
                if target.blank() || target.seq.is_empty() {
                    continue;
                }
                lanes[lane_count] = simd_score::ScoreTarget {
                    subject: target.seq.as_ref(),
                    d_begin: if p.flags.any(Flags::FULL_MATRIX) {
                        -(target.seq.len() as i32 - 1)
                    } else {
                        target.d_begin
                    },
                    d_end: if p.flags.any(Flags::FULL_MATRIX) {
                        p.query.len() as i32
                    } else {
                        target.d_end
                    },
                };
                lane_targets[lane_count] = target;
                lane_count += 1;
            }
            if lane_count == 0 {
                continue;
            }
            let mut subjects = [&[][..]; 16];
            for lane in 0..lane_count {
                subjects[lane] = lanes[lane].subject;
            }
            let batch = if p.flags.any(Flags::FULL_MATRIX) {
                simd_score::score_full_batch_avx2(
                    p.query,
                    &subjects[..lane_count],
                    p.score_matrix,
                    bias,
                    &mut simd_scratch,
                )
                .or_else(|| {
                    simd_score_portable::score_full_batch_portable_i16(
                        p.query,
                        &subjects[..lane_count],
                        p.score_matrix,
                        bias,
                        &mut portable_scratch,
                    )
                })
            } else {
                simd_score::score_batch_avx2(
                    p.query,
                    &lanes[..lane_count],
                    p.score_matrix,
                    bias,
                    &mut simd_scratch,
                )
                .or_else(|| {
                    simd_score_portable::score_batch_portable_i16(
                        p.query,
                        &lanes[..lane_count],
                        p.score_matrix,
                        bias,
                        &mut portable_scratch,
                    )
                })
            };
            let Some(batch) = batch else {
                for &target in &lane_targets[..lane_count] {
                    let (d_begin, d_end) = if p.flags.any(Flags::FULL_MATRIX) {
                        (-(target.seq.len() as i32 - 1), p.query.len() as i32)
                    } else {
                        (target.d_begin, target.d_end)
                    };
                    let score = banded_sw_cbs_score(
                        p.query,
                        target.seq.as_ref(),
                        d_begin,
                        d_end,
                        p.score_matrix,
                        bias,
                        None,
                        1,
                        &mut score_scratch,
                    );
                    if let Some(hsp) = score_only_hsp(target, score, d_begin, d_end, p) {
                        out.push(hsp);
                    }
                }
                continue;
            };
            debug_assert_eq!(batch.len, lane_count);
            for (lane, &target) in lane_targets[..lane_count].iter().enumerate() {
                let (d_begin, d_end) = if p.flags.any(Flags::FULL_MATRIX) {
                    (-(target.seq.len() as i32 - 1), p.query.len() as i32)
                } else {
                    (target.d_begin, target.d_end)
                };
                if batch.overflow_mask & (1 << lane) != 0
                    && p.swipe_bin >= 0
                    && p.swipe_bin as usize % SCORE_BINS == 1
                {
                    overflow.push_back(target.clone());
                    continue;
                }
                let score = if batch.overflow_mask & (1 << lane) != 0 {
                    banded_sw_cbs_score(
                        p.query,
                        target.seq.as_ref(),
                        d_begin,
                        d_end,
                        p.score_matrix,
                        bias,
                        None,
                        1,
                        &mut score_scratch,
                    )
                } else {
                    batch.scores[lane]
                };
                if let Some(hsp) = score_only_hsp(target, score, d_begin, d_end, p) {
                    out.push(hsp);
                }
            }
        }
        return out;
    }

    // Final-pass traceback uses the same vector recurrence as score-only
    // SWIPE, but records one compact trace byte per live lane/cell.
    if p.v != HspValues::NONE {
        let score_bin = p.swipe_bin.max(0) as usize % SCORE_BINS;
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        let avx2 = std::arch::is_x86_feature_detected!("avx2");
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        let avx2 = false;
        let lane_width = if !avx2 && simd_trace_narrow_portable::available() {
            simd_trace_narrow_portable::lane_width(score_bin)
        } else {
            match score_bin {
                0 => 32,
                1 => 16,
                _ => 8,
            }
        };
        let mut reversed_sequences: [Vec<Letter>; 32] = std::array::from_fn(|_| Vec::new());
        for chunk in subject_begin.chunks(lane_width) {
            let lane_targets: Vec<&DpTarget> = chunk
                .iter()
                .filter(|target| !target.blank() && !target.seq.is_empty())
                .collect();
            if lane_targets.is_empty() {
                continue;
            }
            if p.reverse_targets {
                for (lane, target) in lane_targets.iter().enumerate() {
                    reversed_sequences[lane].clear();
                    reversed_sequences[lane].extend(target.seq.iter().rev().copied());
                }
            }
            let trace_targets: Vec<_> = lane_targets
                .iter()
                .enumerate()
                .map(|(lane, target)| {
                    let subject = if p.reverse_targets {
                        reversed_sequences[lane].as_slice()
                    } else {
                        target.seq.as_ref()
                    };
                    let (d_begin, d_end) = if p.flags.any(Flags::FULL_MATRIX) {
                        (-(subject.len() as i32 - 1), p.query.len() as i32)
                    } else {
                        (target.d_begin, target.d_end)
                    };
                    simd_trace::TraceTarget {
                        subject,
                        d_begin,
                        d_end,
                        matrix: target.matrix.as_deref(),
                        matrix_scale: target.matrix_scale(),
                    }
                })
                .collect();
            let Some(batch) = simd_trace::trace_batch_tier_avx2(
                p.query,
                &trace_targets,
                p.score_matrix,
                bias,
                score_bin,
            )
            .or_else(|| {
                let narrow = match score_bin {
                    0 => simd_trace_narrow_portable::trace_batch_i8(
                        p.query,
                        &trace_targets,
                        p.score_matrix,
                        bias,
                    ),
                    1 => simd_trace_narrow_portable::trace_batch_i16(
                        p.query,
                        &trace_targets,
                        p.score_matrix,
                        bias,
                    ),
                    _ => None,
                };
                narrow.map(|batch| simd_trace::TraceBatch {
                    results: batch.results,
                    overflow_mask: batch.overflow_mask,
                })
            })
            .or_else(|| {
                simd_trace_portable::trace_batch_portable(
                    p.query,
                    &trace_targets,
                    p.score_matrix,
                    bias,
                )
                .map(|results| simd_trace::TraceBatch {
                    results,
                    overflow_mask: 0,
                })
            }) else {
                // Runtime fallback is handled by the scalar loop below. This
                // only occurs on a non-AVX2 host or unsupported input.
                for (&target, trace_target) in lane_targets.iter().zip(&trace_targets) {
                    let sw = banded_sw_cbs_range(
                        p.query,
                        trace_target.subject,
                        trace_target.d_begin,
                        trace_target.d_end,
                        p.score_matrix,
                        bias,
                        target.matrix.as_deref(),
                        target.matrix_scale(),
                        &mut traceback_scratch,
                    );
                    if let Some(hsp) = traceback_hsp(
                        target,
                        trace_target.subject,
                        trace_target.d_begin,
                        trace_target.d_end,
                        sw,
                        p,
                    ) {
                        out.push(hsp);
                    }
                }
                continue;
            };
            for (lane, ((&target, trace_target), sw)) in lane_targets
                .iter()
                .zip(&trace_targets)
                .zip(batch.results)
                .enumerate()
            {
                if batch.overflow_mask & (1 << lane) != 0 {
                    overflow.push_back(target.clone());
                    continue;
                }
                if let Some(hsp) = traceback_hsp(
                    target,
                    trace_target.subject,
                    trace_target.d_begin,
                    trace_target.d_end,
                    sw,
                    p,
                ) {
                    out.push(hsp);
                }
            }
        }
        return out;
    }

    for target in subject_begin {
        if target.blank() || target.seq.is_empty() {
            continue;
        }
        let mut reversed_seq;
        let target_seq: &[Letter] = if p.reverse_targets {
            reversed_seq = target.seq.to_vec();
            reversed_seq.reverse();
            reversed_seq.as_slice()
        } else {
            target.seq.as_ref()
        };
        let (d_begin, d_end) = if p.flags.any(Flags::FULL_MATRIX) {
            (-(target_seq.len() as i32 - 1), p.query.len() as i32)
        } else {
            (target.d_begin, target.d_end)
        };
        if p.v == HspValues::NONE {
            let sw_score = banded_sw_cbs_score(
                p.query,
                target_seq,
                d_begin,
                d_end,
                p.score_matrix,
                bias,
                target.matrix.as_deref(),
                target.matrix_scale(),
                &mut score_scratch,
            );
            if let Some(hsp) = score_only_hsp(target, sw_score, d_begin, d_end, p) {
                out.push(hsp);
            }
            continue;
        }
        let sw = banded_sw_cbs_range(
            p.query,
            target_seq,
            d_begin,
            d_end,
            p.score_matrix,
            bias,
            target.matrix.as_deref(),
            target.matrix_scale(),
            &mut traceback_scratch,
        );
        if let Some(hsp) = traceback_hsp(target, target_seq, d_begin, d_end, sw, p) {
            out.push(hsp);
        }
    }
    out
}

pub fn swipe_threads(
    begin: &[DpTarget],
    overflow: &mut TargetVec,
    _round: i32,
    _bin: usize,
    p: &Params<'_>,
) -> Vec<Hsp> {
    dispatch_swipe(begin, overflow, p)
}

pub fn swipe_bin(
    bin: usize,
    begin: &mut [DpTarget],
    round: i32,
    p: &mut Params<'_>,
) -> (Vec<Hsp>, TargetVec) {
    if begin.is_empty() {
        return (Vec::new(), TargetVec::default());
    }
    let timer = Instant::now();
    p.inc_stat(StatValue::SwipeTasksTotal, 1);
    let extension_stat = match bin % SCORE_BINS {
        0 => StatValue::Ext8,
        1 => StatValue::Ext16,
        _ => StatValue::Ext32,
    };
    p.inc_stat(extension_stat, begin.len() as i64);
    let cells = begin
        .iter()
        .map(|target| {
            if p.flags.any(Flags::FULL_MATRIX) {
                p.query.len() as i64 * target.seq.len() as i64
            } else {
                p.query.len() as i64 * (target.d_end - target.d_begin).max(0) as i64
            }
        })
        .sum::<i64>();
    p.inc_stat(StatValue::GrossDpCells, cells);
    p.inc_stat(StatValue::NetDpCells, cells);
    let mut overflow = TargetVec::default();
    if !p.flags.any(Flags::FULL_MATRIX) {
        sort(begin, p.band_bin, p.col_bin);
    }
    let out = swipe_threads(begin, &mut overflow, round, bin, p);
    let time_stat = if p.v.any(HspValues::TRANSCRIPT) {
        StatValue::TimeTracebackSw
    } else {
        StatValue::TimeSw
    };
    p.inc_stat(
        time_stat,
        timer.elapsed().as_micros().min(i64::MAX as u128) as i64,
    );
    (out, overflow)
}

pub fn recompute_reversed(hsps: &mut [Hsp], p: &mut Params<'_>) -> Vec<Hsp> {
    let mut dp_targets = targets();
    let qlen = p.query.len() as i32;

    for hsp in hsps.iter() {
        let qcov = hsp.query_cover(p.query_source_len);
        let tcov = hsp.subject_cover(hsp.target_seq.len() as i32);
        let qcov_filter = if p.query_or_target_cover > 0.0 {
            p.query_or_target_cover
        } else {
            p.query_cover
        };
        let tcov_filter = if p.query_or_target_cover > 0.0 {
            p.query_or_target_cover
        } else {
            p.subject_cover
        };
        let (query_min_len, subject_min_len) =
            hsp.min_range_len(qcov_filter, tcov_filter, qlen, hsp.target_seq.len() as i32);
        let qa = stats::approx_id(hsp.score, query_min_len, 0);
        let ta = stats::approx_id(hsp.score, subject_min_len, 0);
        if qcov < p.query_cover
            || tcov < p.subject_cover
            || qcov.max(tcov) < p.query_or_target_cover
            || (p.query_or_target_cover == 0.0 && qa.min(ta) < p.approx_min_id)
            || (p.query_or_target_cover > 0.0 && qa.max(ta) < p.approx_min_id)
        {
            continue;
        }
        let tlen = hsp.subject_range.end;
        if tlen <= 0 || hsp.swipe_bin < 0 {
            continue;
        }
        let reversed_seq = hsp.target_seq[..tlen as usize].to_vec();
        let band = if p.flags.any(Flags::FULL_MATRIX) {
            qlen
        } else {
            hsp.d_end - hsp.d_begin
        };
        let b = bin(
            p.v,
            band,
            hsp.score,
            0,
            i64::MAX,
            0,
            mismatch_est(hsp.query_range.end, tlen, hsp.length, p.v),
            p.cutoff_score_8bit,
            p.max_swipe_dp,
            p.approx_backtrace,
        )
        .max(hsp.swipe_bin as usize);
        debug_assert!(b >= SCORE_BINS);
        if b >= BINS {
            continue;
        }
        let carry_over = CarryOver::new(
            hsp.query_range.end,
            hsp.subject_range.end,
            hsp.identities,
            hsp.length,
        );
        let mut target = DpTarget::new(
            reversed_seq,
            hsp.target_seq.len() as i32,
            geo::rev_diag(hsp.d_end - 1, qlen, tlen),
            geo::rev_diag(hsp.d_begin, qlen, tlen) + 1,
            hsp.swipe_target as i64,
            qlen,
            carry_over,
            Anchor::default(),
        );
        if let Some(matrix) = &hsp.matrix {
            target = target.with_matrix(matrix.clone(), p.cbs_matrix_scale);
        }
        dp_targets[b].push_back(target);
    }

    let mut reversed_query = p.query.to_vec();
    reversed_query.reverse();
    let reversed_bias;
    let composition_bias = match p.composition_bias {
        Some(bias) => {
            reversed_bias = reversed_composition_bias(bias, p.query.len());
            Some(reversed_bias.as_slice())
        }
        None => None,
    };
    let mut params = Params {
        query: &reversed_query,
        composition_bias,
        reverse_targets: true,
        target_max_len: 0,
        ..p.clone()
    };
    let mut out = Vec::new();
    let mut overflow_targets: Option<Targets> = None;
    for bin_idx in SCORE_BINS..BINS {
        params.target_max_len = dp_targets[bin_idx].max_len();
        params.swipe_bin = bin_idx as i32;
        let (mut bin_out, overflow) =
            swipe_bin(bin_idx, dp_targets[bin_idx].as_mut_slice(), 1, &mut params);
        if !overflow.empty() && bin_idx + 1 < BINS {
            let targets = overflow_targets.get_or_insert_with(targets);
            for target in overflow.as_slice() {
                for hsp in hsps.iter() {
                    if hsp.swipe_target as i64 == target.target_idx {
                        let mut overflow_target = DpTarget::new(
                            hsp.target_seq.clone(),
                            hsp.target_seq.len() as i32,
                            hsp.d_begin,
                            hsp.d_end,
                            hsp.swipe_target as i64,
                            p.query.len() as i32,
                            CarryOver::default(),
                            Anchor::default(),
                        );
                        if let Some(matrix) = &hsp.matrix {
                            overflow_target =
                                overflow_target.with_matrix(matrix.clone(), p.cbs_matrix_scale);
                        }
                        targets[bin_idx + 1].push_back(overflow_target);
                    }
                }
            }
        }
        out.append(&mut bin_out);
    }
    if let Some(overflow_targets) = overflow_targets {
        out.extend(swipe(&overflow_targets, p));
    }
    out
}

fn reversed_composition_bias(bias: &[i8], query_len: usize) -> Vec<i8> {
    let mut out: Vec<i8> = bias[..query_len.min(bias.len())].to_vec();
    out.reverse();
    out.extend(std::iter::repeat_n(0, 32));
    out
}

pub fn swipe(targets: &Targets, p: &mut Params<'_>) -> Vec<Hsp> {
    let mut result: (Vec<Hsp>, TargetVec) = (Vec::new(), TargetVec::default());
    let mut out = Vec::new();
    let mut out_tmp = Vec::new();
    for algo_bin in 0..ALGO_BINS {
        for score_bin in 0..SCORE_BINS {
            let bin = algo_bin * SCORE_BINS + score_bin;
            let mut round_targets = TargetVec::default();
            round_targets.reserve(targets[bin].size() + result.1.size());
            round_targets.push_back_vec(&targets[bin]);
            round_targets.push_back_vec(&result.1);
            p.target_max_len = round_targets.max_len();
            p.swipe_bin = bin as i32;
            result = swipe_bin(bin, round_targets.as_mut_slice(), 0, p);
            if algo_bin == 0 {
                out.append(&mut result.0);
            } else {
                out_tmp.append(&mut result.0);
            }
        }
        debug_assert!(result.1.empty());
    }
    if !out_tmp.is_empty() {
        out.append(&mut recompute_reversed(&mut out_tmp, p));
    }
    out
}

pub fn swipe_set(subjects: &SequenceSet, p: &mut Params<'_>) -> Vec<Hsp> {
    let b = bin(
        p.v,
        0,
        0,
        0,
        0,
        0,
        0,
        p.cutoff_score_8bit,
        p.max_swipe_dp,
        p.approx_backtrace,
    );
    let mut round_targets = TargetVec::default();
    round_targets.reserve(subjects.len() as i64);
    for i in 0..subjects.len() {
        let seq = subjects.get(i).to_vec();
        let qlen = p.query.len() as i32;
        let tlen = seq.len() as i32;
        round_targets.push_back(DpTarget::new(
            seq,
            tlen,
            -(tlen - 1),
            qlen,
            i as i64,
            qlen,
            CarryOver::default(),
            Anchor::default(),
        ));
    }
    let (mut out, overflow) = swipe_bin(b.min(BINS - 1), round_targets.as_mut_slice(), 0, p);
    if reversed(p.v) {
        out = recompute_reversed(&mut out, p);
    }
    if b < BINS - 1 && !overflow.empty() {
        let mut targets = targets();
        targets[b + 1] = overflow;
        out.extend(swipe(&targets, p));
    }
    out
}

fn banded_sw_cbs_range(
    query: &[Letter],
    subject: &[Letter],
    d_begin: i32,
    d_end: i32,
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    target_matrix: Option<&TargetMatrix>,
    matrix_scale: i32,
    scratch: &mut TracebackScratch,
) -> SwResult {
    let qlen = query.len();
    let slen = subject.len();
    let gap_open = (score_matrix.gap_open() + score_matrix.gap_extend()) * matrix_scale;
    let gap_extend = score_matrix.gap_extend() * matrix_scale;
    let neg_inf = i32::MIN / 4;
    let use_cbs = target_matrix.is_none() && !query_cbs.is_empty();
    let rows = qlen + 1;
    let band_rows = ((d_end - d_begin + 2).max(1) as usize).min(rows);
    let col_lower = |j: usize| -> usize {
        if j == 0 {
            1
        } else {
            let target_pos = j as i32 - 1;
            (d_begin + target_pos).max(0) as usize + 1
        }
    };
    let band_offset = |i: usize, j: usize| -> Option<usize> {
        if j > slen {
            return None;
        }
        let lower = col_lower(j);
        if i < lower {
            return None;
        }
        let k = i - lower;
        if k < band_rows {
            Some(k)
        } else {
            None
        }
    };
    let band_idx =
        |i: usize, j: usize| -> Option<usize> { band_offset(i, j).map(|k| j * band_rows + k) };
    const ACTIVE: u8 = 1 << 0;
    const GAP_V: u8 = 1 << 1;
    const GAP_H: u8 = 1 << 2;
    const OPEN_V: u8 = 1 << 3;
    const OPEN_H: u8 = 1 << 4;
    let cells = band_rows * (slen + 1);
    scratch.prepare(cells, band_rows, neg_inf);
    let TracebackScratch {
        trace,
        prev_h,
        curr_h,
        prev_e,
        curr_e,
        f_col,
    } = scratch;
    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for j in 1..=slen {
        curr_h.fill(0);
        curr_e.fill(neg_inf);
        f_col.fill(neg_inf);
        let target_pos = j as i32 - 1;
        let lower = (d_begin + target_pos).max(0);
        let upper = (d_end + target_pos - 1).min(qlen as i32 - 1);
        if lower > upper {
            std::mem::swap(prev_h, curr_h);
            std::mem::swap(prev_e, curr_e);
            continue;
        }
        let lower_i = lower as usize + 1;
        let upper_i = upper as usize + 1;
        let prev_lower = col_lower(j - 1);
        let curr_col_offset = j * band_rows;
        for i in lower_i..=upper_i {
            let qpos = i - 1;
            let spos = j - 1;
            let ql = query[qpos];
            let sl = subject[spos];
            // C++ `Sequence::operator[]` strips query soft masks before
            // SWIPE scoring. A masked subject/profile lane still scores as
            // zero.
            let match_score = if sl & SEED_MASK != 0 {
                0
            } else if let Some(matrix) = target_matrix {
                matrix.scores[(sl & LETTER_MASK) as usize * 32 + (ql & LETTER_MASK) as usize] as i32
            } else {
                score_matrix.score(ql, sl)
            };
            let current_offset = i - lower_i;
            let cbs = if use_cbs { query_cbs[qpos] as i32 } else { 0 };
            let diag_score = if i > prev_lower {
                let diag_offset = i - 1 - prev_lower;
                if diag_offset < band_rows {
                    prev_h[diag_offset]
                } else {
                    0
                }
            } else {
                0
            } + match_score
                + cbs;
            let e_in = if i >= prev_lower {
                let prev_offset = i - prev_lower;
                if prev_offset < band_rows {
                    prev_e[prev_offset]
                } else {
                    neg_inf
                }
            } else {
                neg_inf
            };
            let f_in = if i > lower_i {
                f_col[current_offset - 1]
            } else {
                neg_inf
            };
            let score = diag_score.max(e_in).max(f_in).max(0);
            curr_h[current_offset] = score;

            let e_extend = e_in - gap_extend;
            let f_extend = f_in - gap_extend;
            let open = score - gap_open;
            let e_val = e_extend.max(open);
            let f_val = f_extend.max(open);
            curr_e[current_offset] = e_val;
            f_col[current_offset] = f_val;
            let trace_idx = curr_col_offset + current_offset;
            trace[trace_idx] = u8::from(score > 0) * ACTIVE
                | u8::from(score == f_in) * GAP_V
                | u8::from(score == e_in) * GAP_H
                | u8::from(f_val == open) * OPEN_V
                | u8::from(e_val == open) * OPEN_H;
            // Best-cell tie-breaking: within a column, the LATEST tied row wins
            // (matches C++ `VectorRowCounter::inc` in `cell_update.h:45-47`:
            //  `i_max = blend(i_max, i, best == current_cell)` overwrites on
            //  equality). Between columns, the earliest tied column wins
            //  (strict `>` on the column-max in `banded_swipe.h:323`). So
            //  update when the score is strictly greater, OR when it ties
            //  the current best AND we're still in the same column.
            if score > best_score || (score == best_score && j == best_j) {
                best_score = score;
                best_i = i;
                best_j = j;
            }
        }
        std::mem::swap(prev_h, curr_h);
        std::mem::swap(prev_e, curr_e);
    }

    if best_score == 0 {
        return SwResult::default();
    }

    let mut i = best_i;
    let mut j = best_j;
    let mut result = SwResult {
        score: best_score,
        query_end: i as i32,
        subject_end: j as i32,
        ..Default::default()
    };
    let mut ops = Vec::new();
    let get_trace = |i: usize, j: usize| -> u8 { band_idx(i, j).map_or(0, |idx| trace[idx]) };

    while i > 0 && j > 0 && get_trace(i, j) & ACTIVE != 0 {
        let cell_trace = get_trace(i, j);
        let ql = query[i - 1];
        let sl = subject[j - 1];
        // Match C++ TracebackVectorMatrix::walk_gap(): prefer vertical gap
        // masks, then horizontal masks, then diagonal.
        if cell_trace & GAP_V != 0 {
            let mut gap_len = 0i32;
            loop {
                gap_len += 1;
                if i == 0 {
                    break;
                }
                i -= 1;
                if get_trace(i, j) & OPEN_V != 0 || i == 0 {
                    break;
                }
            }
            ops.push((EditOperation::Insertion, gap_len));
            result.gap_openings += 1;
            result.gaps += gap_len;
            result.length += gap_len;
        } else if cell_trace & GAP_H != 0 {
            let mut gap_len = 0i32;
            loop {
                gap_len += 1;
                if j == 0 {
                    break;
                }
                j -= 1;
                if get_trace(i, j) & OPEN_H != 0 || j == 0 {
                    break;
                }
            }
            ops.push((EditOperation::Deletion, gap_len));
            result.gap_openings += 1;
            result.gaps += gap_len;
            result.length += gap_len;
        } else {
            if (ql & LETTER_MASK) == (sl & LETTER_MASK) {
                ops.push((EditOperation::Match, 1));
                result.identities += 1;
            } else {
                ops.push((EditOperation::Substitution, 1));
                result.mismatches += 1;
            }
            result.length += 1;
            i -= 1;
            j -= 1;
        }
    }

    result.query_begin = i as i32;
    result.subject_begin = j as i32;
    ops.reverse();
    result.operations = ops;
    result
}

fn banded_sw_cbs_score(
    query: &[Letter],
    subject: &[Letter],
    d_begin: i32,
    d_end: i32,
    score_matrix: &ScoreMatrix,
    query_cbs: &[i8],
    target_matrix: Option<&TargetMatrix>,
    matrix_scale: i32,
    scratch: &mut ScoreScratch,
) -> i32 {
    let qlen = query.len();
    let slen = subject.len();
    let gap_open = (score_matrix.gap_open() + score_matrix.gap_extend()) * matrix_scale;
    let gap_extend = score_matrix.gap_extend() * matrix_scale;
    let neg_inf = i32::MIN / 4;
    let use_cbs = target_matrix.is_none() && !query_cbs.is_empty();
    let rows = qlen + 1;
    let band_rows = ((d_end - d_begin + 2).max(1) as usize).min(rows);
    scratch.prepare(band_rows, neg_inf);
    let mut best_score = 0i32;
    let mut prev_begin = usize::MAX;
    let mut prev_end = 0usize;

    for j in 1..=slen {
        let target_pos = j as i32 - 1;
        let lower = (d_begin + target_pos).max(0);
        let upper = (d_end + target_pos - 1).min(qlen as i32 - 1);
        if lower > upper {
            prev_begin = usize::MAX;
            prev_end = 0;
            continue;
        }
        let lower_i = lower as usize + 1;
        let upper_i = upper as usize + 1;
        for i in lower_i..=upper_i {
            let qpos = i - 1;
            let spos = j - 1;
            let ql = query[qpos];
            let sl = subject[spos];
            let match_score = if sl & SEED_MASK != 0 {
                0
            } else if let Some(matrix) = target_matrix {
                matrix.scores[(sl & LETTER_MASK) as usize * 32 + (ql & LETTER_MASK) as usize] as i32
            } else {
                score_matrix.score(ql, sl)
            };
            let current_offset = i - lower_i;
            let cbs = if use_cbs { query_cbs[qpos] as i32 } else { 0 };
            let diag_score = if i > prev_begin {
                let diag_offset = i - 1 - prev_begin;
                if i - 1 < prev_end {
                    scratch.prev_h[diag_offset]
                } else {
                    0
                }
            } else {
                0
            } + match_score
                + cbs;
            let e_in = if i >= prev_begin {
                let prev_offset = i - prev_begin;
                if i < prev_end {
                    scratch.prev_e[prev_offset]
                } else {
                    neg_inf
                }
            } else {
                neg_inf
            };
            let f_in = if i > lower_i {
                scratch.f[current_offset - 1]
            } else {
                neg_inf
            };
            let score = diag_score.max(e_in).max(f_in).max(0);
            scratch.curr_h[current_offset] = score;
            let e_extend = e_in - gap_extend;
            let f_extend = f_in - gap_extend;
            let open = score - gap_open;
            scratch.curr_e[current_offset] = e_extend.max(open);
            scratch.f[current_offset] = f_extend.max(open);
            best_score = best_score.max(score);
        }
        std::mem::swap(&mut scratch.prev_h, &mut scratch.curr_h);
        std::mem::swap(&mut scratch.prev_e, &mut scratch.curr_e);
        prev_begin = lower_i;
        prev_end = upper_i + 1;
    }

    best_score
}

#[derive(Default)]
struct ScoreScratch {
    prev_h: Vec<i32>,
    curr_h: Vec<i32>,
    prev_e: Vec<i32>,
    curr_e: Vec<i32>,
    f: Vec<i32>,
}

#[derive(Default)]
struct TracebackScratch {
    trace: Vec<u8>,
    prev_h: Vec<i32>,
    curr_h: Vec<i32>,
    prev_e: Vec<i32>,
    curr_e: Vec<i32>,
    f_col: Vec<i32>,
}

impl TracebackScratch {
    fn prepare(&mut self, cells: usize, band_rows: usize, neg_inf: i32) {
        self.trace.clear();
        self.trace.resize(cells, 0);
        self.prev_h.clear();
        self.prev_h.resize(band_rows, 0);
        self.curr_h.clear();
        self.curr_h.resize(band_rows, 0);
        self.prev_e.clear();
        self.prev_e.resize(band_rows, neg_inf);
        self.curr_e.clear();
        self.curr_e.resize(band_rows, neg_inf);
        self.f_col.clear();
        self.f_col.resize(band_rows, neg_inf);
    }
}

impl ScoreScratch {
    fn prepare(&mut self, rows: usize, neg_inf: i32) {
        self.prev_h.clear();
        self.prev_h.resize(rows, 0);
        self.curr_h.clear();
        self.curr_h.resize(rows, 0);
        self.prev_e.clear();
        self.prev_e.resize(rows, neg_inf);
        self.curr_e.clear();
        self.curr_e.resize(rows, neg_inf);
        self.f.clear();
        self.f.resize(rows, neg_inf);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::AMINO_ACID_COUNT;

    #[test]
    fn test_dp_target_banded_cols() {
        assert_eq!(DpTarget::banded_cols(10, 10, -1, 2), 10);
        assert_eq!(DpTarget::banded_cols(10, 5, 2, 5), 5);
    }

    #[test]
    fn test_hsp_values_have_coords() {
        assert!(have_coords(HspValues::TRANSCRIPT));
        assert!(have_coords(HspValues::COORDS));
        assert!(!have_coords(HspValues::QUERY_COORDS));
    }

    #[test]
    fn test_bin_coords_promotes_to_reversed_bin() {
        let b = bin(
            HspValues::COORDS,
            100,
            40,
            0,
            100,
            0,
            0,
            127,
            i64::MAX,
            false,
        );
        assert_eq!(b, 3);
    }

    #[test]
    fn test_mismatch_est() {
        assert_eq!(mismatch_est(10, 7, 0, HspValues::MISMATCHES), 7);
        assert_eq!(mismatch_est(10, 7, 4, HspValues::MISMATCHES), 4);
        assert_eq!(mismatch_est(10, 7, 4, HspValues::COORDS), 0);
    }

    #[test]
    fn test_swipe_single_target_self_alignment() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6];
        let mut ts = targets();
        ts[0].push_back(DpTarget::new(
            query.clone(),
            query.len() as i32,
            -1,
            2,
            0,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        ));
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let statistics = Arc::new(Mutex::new(Statistics::new()));
        params.statistics = Some(statistics.clone());
        let out = swipe(&ts, &mut params);
        assert_eq!(out.len(), 1);
        assert!(out[0].score > 0);
        assert_eq!(out[0].query_range, Interval::new(0, query.len() as i32));
        assert_eq!(out[0].subject_range, Interval::new(0, query.len() as i32));
        let statistics = statistics.lock().unwrap();
        assert_eq!(statistics.get(StatValue::SwipeTasksTotal), 1);
        assert!(statistics.get(StatValue::GrossDpCells) > 0);
        assert!(statistics.get(StatValue::TimeTracebackSw) >= 0);
    }

    #[test]
    fn test_i8_saturation_promotes_target_to_i16_bin() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![17; 30];
        let target = DpTarget::new(
            query.clone(),
            query.len() as i32,
            0,
            1,
            9,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        );
        let mut direct_params = Params::new(&query, &sm);
        direct_params.swipe_bin = 0;
        direct_params.v = HspValues::NONE;
        let mut overflow = TargetVec::default();
        let direct = dispatch_swipe(std::slice::from_ref(&target), &mut overflow, &direct_params);
        assert!(
            direct.is_empty(),
            "direct={:?}, overflow={}",
            direct.iter().map(|h| h.score).collect::<Vec<_>>(),
            overflow.size()
        );
        assert_eq!(overflow.size(), 1);

        let mut bins = targets();
        bins[0].push_back(target);
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::NONE;
        let output = swipe(&bins, &mut params);
        assert_eq!(output.len(), 1);
        assert_eq!(output[0].score, 330);
        assert_eq!(output[0].swipe_target, 9);
    }

    #[test]
    fn test_traceback_i8_saturation_promotes_target_to_i16_bin() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![17; 30];
        let target = DpTarget::new(
            query.clone(),
            query.len() as i32,
            0,
            1,
            10,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        );
        let mut direct_params = Params::new(&query, &sm);
        direct_params.swipe_bin = 0;
        direct_params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let mut overflow = TargetVec::default();
        let batch = vec![target.clone(), target.clone(), target.clone()];
        let direct = dispatch_swipe(&batch, &mut overflow, &direct_params);
        assert!(direct.is_empty());
        assert_eq!(overflow.size(), 3);

        let mut bins = targets();
        bins[0].push_back_vec(&TargetVec {
            targets: batch,
            max_len: 30,
        });
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let output = swipe(&bins, &mut params);
        assert_eq!(output.len(), 3);
        assert_eq!(output[0].score, 330);
        assert_eq!(output[0].swipe_bin, 1);
        assert_eq!(output[0].query_range, Interval::new(0, 30));
        assert_eq!(output[0].subject_range, Interval::new(0, 30));
    }

    #[test]
    fn test_semi_global_i8_promotion_keeps_delta_zero_score() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![17; 20];
        let mut bins = targets();
        bins[0].push_back(DpTarget::new(
            query.clone(),
            query.len() as i32,
            0,
            1,
            5,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        ));
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::NONE;
        params.flags = Flags::SEMI_GLOBAL;
        let output = swipe(&bins, &mut params);
        assert_eq!(output.len(), 1);
        assert_eq!(output[0].score, 220);
        assert_eq!(output[0].swipe_target, 5);
        assert_eq!(output[0].swipe_bin, 1);
    }

    #[test]
    fn test_swipe_set_self_alignment() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6];
        let mut subjects = SequenceSet::new();
        subjects.push(&query);
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let out = swipe_set(&subjects, &mut params);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].swipe_target, 0);
        assert_eq!(out[0].query_range, Interval::new(0, query.len() as i32));
        assert_eq!(out[0].subject_range, Interval::new(0, query.len() as i32));
    }

    #[test]
    fn test_full_matrix_swipe_set_self_alignment() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6];
        let mut subjects = SequenceSet::new();
        subjects.push(&query);
        let mut params = Params::new(&query, &sm);
        params.flags = Flags::FULL_MATRIX;
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let out = swipe_set(&subjects, &mut params);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].swipe_target, 0);
        assert_eq!(out[0].query_range, Interval::new(0, query.len() as i32));
        assert_eq!(out[0].subject_range, Interval::new(0, query.len() as i32));
    }

    #[test]
    fn test_recompute_reversed_self_alignment() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6];
        let mut subjects = SequenceSet::new();
        subjects.push(&query);
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::COORDS;
        let out = swipe_set(&subjects, &mut params);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].swipe_target, 0);
        assert_eq!(out[0].query_range, Interval::new(0, query.len() as i32));
        assert_eq!(out[0].subject_range, Interval::new(0, query.len() as i32));
        assert!(out[0].approx_id >= 99.0);
    }

    #[test]
    fn test_adjusted_target_matrix_scores_without_cbs_or_final_rescale() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2];
        let mut scores = vec![-5i8; 32 * AMINO_ACID_COUNT];
        for i in 0..AMINO_ACID_COUNT {
            scores[i * 32 + i] = 3;
        }
        let matrix = Arc::new(TargetMatrix::new(scores, -5, 3));
        let target = DpTarget::new(
            query.clone(),
            query.len() as i32,
            -1,
            2,
            0,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        )
        .with_matrix(matrix, 2);
        assert!(target.adjusted_matrix());
        assert_eq!(target.matrix_scale(), 2);

        let cbs = [10i8; 3];
        let mut params = Params::new(&query, &sm);
        params.composition_bias = Some(&cbs);
        params.cbs_matrix_scale = 2;
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let mut overflow = TargetVec::default();
        let out = dispatch_swipe(&[target], &mut overflow, &params);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].score, 9);
        assert!(out[0].matrix.is_some());
        assert_eq!(out[0].identities, 3);
        assert_eq!(out[0].gaps, 0);
    }

    #[test]
    fn test_swipe_uses_true_target_len_for_statistics() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2];
        let true_target_len = 100;
        let target = DpTarget::new(
            query.clone(),
            true_target_len,
            -1,
            2,
            0,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        );
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let mut overflow = TargetVec::default();
        let out = dispatch_swipe(&[target], &mut overflow, &params);
        assert_eq!(out.len(), 1);
        assert_eq!(
            out[0].evalue,
            sm.evalue(out[0].score, query.len() as u32, true_target_len as u32)
        );
        assert_eq!(
            out[0].corrected_bit_score,
            sm.bitscore_corrected(out[0].score, query.len() as u32, true_target_len as u32)
        );
    }

    #[test]
    fn test_traceback_transcript_stores_substitution_subject_letter() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2];
        let subject: Vec<Letter> = vec![0, 3, 2];
        let mut scores = vec![-20i8; 32 * AMINO_ACID_COUNT];
        for i in 0..AMINO_ACID_COUNT {
            scores[i * 32 + i] = 20;
        }
        scores[3 * 32 + 1] = 7;
        let matrix = Arc::new(TargetMatrix::new(scores, -20, 20));
        let target = DpTarget::new(
            subject.clone(),
            subject.len() as i32,
            -1,
            2,
            0,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        )
        .with_matrix(matrix, 1);
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let mut overflow = TargetVec::default();
        let out = dispatch_swipe(&[target], &mut overflow, &params);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].positives, 3);

        let ops: Vec<_> = out[0].transcript.iter().collect();
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0].op, EditOperation::Match);
        assert_eq!(ops[0].count, 1);
        assert_eq!(ops[1].op, EditOperation::Substitution);
        assert_eq!(ops[1].letter, 3);
        assert_eq!(ops[2].op, EditOperation::Match);
        assert_eq!(ops[2].count, 1);
    }

    #[test]
    fn test_traceback_positives_respect_masked_subject_letter() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2];
        let subject: Vec<Letter> = vec![0, 3 | SEED_MASK, 2];
        let mut scores = vec![-20i8; 32 * AMINO_ACID_COUNT];
        for i in 0..AMINO_ACID_COUNT {
            scores[i * 32 + i] = 20;
        }
        scores[3 * 32 + 1] = 7;
        let matrix = Arc::new(TargetMatrix::new(scores, -20, 20));
        let target = DpTarget::new(
            subject,
            query.len() as i32,
            -1,
            2,
            0,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        )
        .with_matrix(matrix, 1);
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let mut overflow = TargetVec::default();
        let out = dispatch_swipe(&[target], &mut overflow, &params);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].identities, 2);
        assert_eq!(out[0].positives, 2);
    }

    #[test]
    fn test_traceback_transcript_stores_deletion_subject_letter() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1];
        let subject: Vec<Letter> = vec![0, 2, 1];
        let mut scores = vec![-20i8; 32 * AMINO_ACID_COUNT];
        for i in 0..AMINO_ACID_COUNT {
            scores[i * 32 + i] = 20;
        }
        let matrix = Arc::new(TargetMatrix::new(scores, -20, 20));
        let target = DpTarget::new(
            subject.clone(),
            subject.len() as i32,
            -2,
            2,
            0,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        )
        .with_matrix(matrix, 1);
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let mut overflow = TargetVec::default();
        let out = dispatch_swipe(&[target], &mut overflow, &params);
        assert_eq!(out.len(), 1);

        let ops: Vec<_> = out[0].transcript.iter().collect();
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0].op, EditOperation::Match);
        assert_eq!(ops[1].op, EditOperation::Deletion);
        assert_eq!(ops[1].letter, 2);
        assert_eq!(ops[2].op, EditOperation::Match);
    }

    #[test]
    fn test_swipe_exact_identity_approx_id_is_percent_identity() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2];
        let mut scores = vec![-5i8; 32 * AMINO_ACID_COUNT];
        for i in 0..AMINO_ACID_COUNT {
            scores[i * 32 + i] = 1;
        }
        let matrix = Arc::new(TargetMatrix::new(scores, -5, 1));
        let target = DpTarget::new(
            query.clone(),
            query.len() as i32,
            -5,
            1,
            0,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        )
        .with_matrix(matrix, 1);
        let mut params = Params::new(&query, &sm);
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let mut overflow = TargetVec::default();
        let out = dispatch_swipe(&[target], &mut overflow, &params);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].score, 3);
        assert_eq!(out[0].approx_id, 100.0);
        assert!(stats::approx_id(out[0].score, 3, 3) < 100.0);
    }

    #[test]
    fn test_swipe_query_source_range_uses_translated_frame() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query: Vec<Letter> = vec![0, 1, 2];
        let target = DpTarget::new(
            query.clone(),
            query.len() as i32,
            -1,
            2,
            0,
            query.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        );
        let mut params = Params::new(&query, &sm);
        params.frame = 1;
        params.query_source_len = 12;
        params.v = HspValues::TRANSCRIPT | HspValues::COORDS;
        let mut overflow = TargetVec::default();
        let out = dispatch_swipe(&[target], &mut overflow, &params);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].query_range, Interval::new(0, 3));
        assert_eq!(out[0].query_source_range, Interval::new(1, 10));
        assert_eq!(out[0].subject_source_range, out[0].subject_range);
    }

    #[test]
    fn test_reversed_composition_bias_preserves_tail_padding() {
        let mut bias = vec![1, 2, 3];
        bias.extend(std::iter::repeat_n(0, 32));
        let reversed = reversed_composition_bias(&bias, 3);
        assert_eq!(&reversed[..3], &[3, 2, 1]);
        assert_eq!(reversed.len(), 35);
        assert!(reversed[3..].iter().all(|&x| x == 0));
    }
}
