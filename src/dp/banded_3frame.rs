//! Scalar banded three-frame Smith-Waterman.
//!
//! This is a direct Rust counterpart for the core recurrence in
//! `diamond/src/dp/swipe/banded_3frame_swipe.cpp`. It evaluates three translated
//! query frames together and permits frame-shift transitions with the matrix
//! frame-shift penalty.

use crate::align::hsp::Hsp;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{Letter, LETTER_MASK};
use crate::dp::smith_waterman::SwResult;
use crate::dp::swipe::DpTarget;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::interval::Interval;

#[derive(Debug, Clone, Default)]
pub struct ThreeFrameSwResult {
    pub sw: SwResult,
    pub frame_begin: usize,
    pub frame_end: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Trace {
    Zero,
    Diag,
    ForwardShift,
    ReverseShift,
    HGapOpen,
    HGapExtend,
    VGapOpen,
    VGapExtend,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DpStat {
    pub gross_cells: usize,
    pub net_cells: usize,
}

#[derive(Debug, Clone)]
pub struct Banded3FrameSwipeMatrix {
    band: usize,
    hgap: Vec<i32>,
    score: Vec<i32>,
}

impl Banded3FrameSwipeMatrix {
    pub fn new(band: usize, _cols: usize) -> Self {
        Banded3FrameSwipeMatrix {
            band,
            hgap: vec![0; band + 3],
            score: vec![0; band + 1],
        }
    }

    pub fn begin(&mut self, offset: usize, _col: usize) -> ColumnIterator<'_> {
        ColumnIterator::new(&mut self.hgap, &mut self.score, offset)
    }

    pub fn band(&self) -> usize {
        self.band
    }
}

pub struct ColumnIterator<'a> {
    hgap: &'a mut [i32],
    score: &'a mut [i32],
    offset: usize,
    pub sm4: i32,
    pub sm3: i32,
    pub sm2: i32,
}

impl<'a> ColumnIterator<'a> {
    pub fn new(hgap: &'a mut [i32], score: &'a mut [i32], offset: usize) -> Self {
        let sm3 = score.get(offset).copied().unwrap_or(0);
        let sm2 = score.get(offset + 1).copied().unwrap_or(0);
        ColumnIterator {
            hgap,
            score,
            offset,
            sm4: 0,
            sm3,
            sm2,
        }
    }

    pub fn next(&mut self) {
        self.offset += 1;
        self.sm4 = self.sm3;
        self.sm3 = self.sm2;
        self.sm2 = self.score.get(self.offset + 1).copied().unwrap_or(0);
    }

    pub fn hgap(&self) -> i32 {
        self.hgap.get(self.offset + 3).copied().unwrap_or(0)
    }

    pub fn set_hgap(&mut self, x: i32) {
        if let Some(slot) = self.hgap.get_mut(self.offset) {
            *slot = x;
        }
    }

    pub fn set_score(&mut self, x: i32) {
        if let Some(slot) = self.score.get_mut(self.offset) {
            *slot = x;
        }
    }

    pub fn set_zero(&mut self) {
        for delta in 1..=3 {
            if let Some(slot) = self
                .offset
                .checked_sub(delta)
                .and_then(|i| self.score.get_mut(i))
            {
                *slot = 0;
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct Banded3FrameSwipeTracebackMatrix {
    band: usize,
    hgap: Vec<i32>,
    score: Vec<i32>,
    cols: usize,
}

impl Banded3FrameSwipeTracebackMatrix {
    pub fn new(band: usize, cols: usize) -> Self {
        Banded3FrameSwipeTracebackMatrix {
            band,
            hgap: vec![0; band + 3],
            score: vec![0; (band + 1) * (cols + 1)],
            cols,
        }
    }

    pub fn begin(&mut self, offset: usize, col: usize) -> TracebackColumnIterator<'_> {
        TracebackColumnIterator::new(&mut self.hgap, &mut self.score, self.band, offset, col)
    }

    pub fn band(&self) -> usize {
        self.band
    }

    pub fn traceback(
        &self,
        col: usize,
        i0: i32,
        j: i32,
        dna_len: i32,
        channel: usize,
        score: i32,
    ) -> TracebackIterator<'_> {
        let i_start = (-i0).max(0) as usize * 3;
        let i_end = self.band.min((dna_len - 2 - i0 * 3).max(0) as usize);
        let base = col.min(self.cols) * (self.band + 1);
        for i in i_start..i_end {
            let idx = base + i + channel;
            if self.score.get(idx).copied().unwrap_or(i32::MIN) == score {
                return TracebackIterator::new(
                    &self.score,
                    self.band,
                    i % 3,
                    i0 + i as i32 / 3,
                    j,
                    idx,
                );
            }
        }
        TracebackIterator::new(&self.score, self.band, 0, i0, j, base + channel)
    }
}

pub struct TracebackColumnIterator<'a> {
    hgap: &'a mut [i32],
    score: &'a mut [i32],
    band: usize,
    offset: usize,
    col: usize,
    pub sm4: i32,
    pub sm3: i32,
    pub sm2: i32,
}

impl<'a> TracebackColumnIterator<'a> {
    pub fn new(
        hgap: &'a mut [i32],
        score: &'a mut [i32],
        band: usize,
        offset: usize,
        col: usize,
    ) -> Self {
        let base = col * (band + 1);
        let sm3 = score.get(base + offset).copied().unwrap_or(0);
        let sm2 = score.get(base + offset + 1).copied().unwrap_or(0);
        TracebackColumnIterator {
            hgap,
            score,
            band,
            offset,
            col,
            sm4: 0,
            sm3,
            sm2,
        }
    }

    pub fn next(&mut self) {
        self.offset += 1;
        self.sm4 = self.sm3;
        self.sm3 = self.sm2;
        let base = self.col * (self.band + 1);
        self.sm2 = self.score.get(base + self.offset + 1).copied().unwrap_or(0);
    }

    pub fn hgap(&self) -> i32 {
        self.hgap.get(self.offset + 3).copied().unwrap_or(0)
    }

    pub fn set_hgap(&mut self, x: i32) {
        if let Some(slot) = self.hgap.get_mut(self.offset) {
            *slot = x;
        }
    }

    pub fn set_score(&mut self, x: i32) {
        let base = (self.col + 1) * (self.band + 1);
        if let Some(slot) = self.score.get_mut(base + self.offset) {
            *slot = x;
        }
    }

    pub fn set_zero(&mut self) {
        let base = (self.col + 1) * (self.band + 1);
        for delta in 1..=3 {
            if let Some(slot) = self
                .offset
                .checked_sub(delta)
                .and_then(|i| self.score.get_mut(base + i))
            {
                *slot = 0;
            }
        }
    }
}

pub struct TracebackIterator<'a> {
    band: usize,
    score: &'a [i32],
    index: usize,
    pub frame: usize,
    pub i: i32,
    pub j: i32,
}

impl<'a> TracebackIterator<'a> {
    pub fn new(score: &'a [i32], band: usize, frame: usize, i: i32, j: i32, index: usize) -> Self {
        TracebackIterator {
            band,
            score,
            index,
            frame,
            i,
            j,
        }
    }

    pub fn score(&self) -> i32 {
        self.score.get(self.index).copied().unwrap_or(0)
    }

    pub fn sm3(&self) -> i32 {
        self.score
            .get(self.index.saturating_sub(self.band + 1))
            .copied()
            .unwrap_or(0)
    }

    pub fn sm4(&self) -> i32 {
        self.score
            .get(self.index.saturating_sub(self.band + 2))
            .copied()
            .unwrap_or(0)
    }

    pub fn sm2(&self) -> i32 {
        self.score
            .get(self.index.saturating_sub(self.band))
            .copied()
            .unwrap_or(0)
    }

    pub fn walk_diagonal(&mut self) {
        self.index = self.index.saturating_sub(self.band + 1);
        self.i -= 1;
        self.j -= 1;
    }

    pub fn walk_forward_shift(&mut self) {
        self.index = self.index.saturating_sub(self.band + 2);
        self.i -= 1;
        self.j -= 1;
        if self.frame == 0 {
            self.frame = 2;
            self.i -= 1;
        } else {
            self.frame -= 1;
        }
    }

    pub fn walk_reverse_shift(&mut self) {
        self.index = self.index.saturating_sub(self.band);
        self.i -= 1;
        self.j -= 1;
        self.frame += 1;
        if self.frame == 3 {
            self.frame = 0;
            self.i += 1;
        }
    }

    pub fn walk_gap(
        &mut self,
        _d0: i32,
        _d1: i32,
        op: EditOperation,
        len: i32,
    ) -> (EditOperation, i32) {
        match op {
            EditOperation::Deletion => self.walk_hgap(len),
            EditOperation::Insertion => self.walk_vgap(len),
            _ => {}
        }
        (op, len)
    }

    pub fn walk_hgap(&mut self, len: i32) {
        self.j -= len;
        self.index = self
            .index
            .saturating_sub((len as usize).saturating_mul(self.band.saturating_sub(2).max(1)));
    }

    pub fn walk_vgap(&mut self, len: i32) {
        self.i -= len;
        self.index = self.index.saturating_sub((len as usize).saturating_mul(3));
    }
}

pub fn banded_3frame_swipe(
    query: [&[Letter]; 3],
    subject: &[Letter],
    query_anchor: usize,
    subject_anchor: usize,
    band_width: i32,
    score_matrix: &ScoreMatrix,
) -> ThreeFrameSwResult {
    let qlen = query.iter().map(|q| q.len()).max().unwrap_or(0);
    let slen = subject.len();
    if qlen == 0 || slen == 0 {
        return ThreeFrameSwResult::default();
    }

    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();
    let frame_shift = score_matrix.frame_shift();
    let neg_inf = i32::MIN / 4;
    let states = (qlen + 1) * (slen + 1) * 3;
    let idx = |frame: usize, i: usize, j: usize| -> usize { ((j * (qlen + 1) + i) * 3) + frame };

    let mut h = vec![0i32; states];
    let mut hgap = vec![neg_inf; states];
    let mut vgap = vec![neg_inf; states];
    let mut trace = vec![Trace::Zero; states];
    let mut best_score = 0i32;
    let mut best_frame = 0usize;
    let mut best_i = 0usize;
    let mut best_j = 0usize;
    let diag = query_anchor as i32 - subject_anchor as i32;

    for j in 1..=slen {
        let subject_pos = j as i32 - 1;
        let lower = (subject_pos + diag - band_width).max(0);
        let upper = (subject_pos + diag + band_width).min(qlen as i32 - 1);
        if lower > upper {
            continue;
        }

        for i in (lower as usize + 1)..=(upper as usize + 1) {
            for frame in 0..3 {
                if i > query[frame].len() {
                    continue;
                }
                let q = query[frame][i - 1] & LETTER_MASK;
                let s = subject[j - 1] & LETTER_MASK;
                let match_score = score_matrix.score(q, s);
                let current = idx(frame, i, j);

                let mut score = h[idx(frame, i - 1, j - 1)] + match_score;
                let mut source = Trace::Diag;

                let fwd_pred_frame = if frame == 0 { 2 } else { frame - 1 };
                let fwd_pred_i = if frame == 0 {
                    i.checked_sub(2)
                } else {
                    i.checked_sub(1)
                };
                if let Some(pi) = fwd_pred_i {
                    if pi <= query[fwd_pred_frame].len() {
                        let candidate =
                            h[idx(fwd_pred_frame, pi, j - 1)] + match_score - frame_shift;
                        if candidate > score {
                            score = candidate;
                            source = Trace::ForwardShift;
                        }
                    }
                }

                let rev_pred_frame = if frame == 2 { 0 } else { frame + 1 };
                let rev_pred_i = if frame == 2 {
                    Some(i)
                } else {
                    i.checked_sub(1)
                };
                if let Some(pi) = rev_pred_i {
                    if pi <= query[rev_pred_frame].len() {
                        let candidate =
                            h[idx(rev_pred_frame, pi, j - 1)] + match_score - frame_shift;
                        if candidate > score {
                            score = candidate;
                            source = Trace::ReverseShift;
                        }
                    }
                }

                let h_open = h[idx(frame, i, j - 1)] - gap_open;
                let h_extend = hgap[idx(frame, i, j - 1)] - gap_extend;
                if h_open >= h_extend {
                    hgap[current] = h_open;
                    if h_open > score {
                        score = h_open;
                        source = Trace::HGapOpen;
                    }
                } else {
                    hgap[current] = h_extend;
                    if h_extend > score {
                        score = h_extend;
                        source = Trace::HGapExtend;
                    }
                }

                let v_open = h[idx(frame, i - 1, j)] - gap_open;
                let v_extend = vgap[idx(frame, i - 1, j)] - gap_extend;
                if v_open >= v_extend {
                    vgap[current] = v_open;
                    if v_open > score {
                        score = v_open;
                        source = Trace::VGapOpen;
                    }
                } else {
                    vgap[current] = v_extend;
                    if v_extend > score {
                        score = v_extend;
                        source = Trace::VGapExtend;
                    }
                }

                if score <= 0 {
                    h[current] = 0;
                    trace[current] = Trace::Zero;
                } else {
                    h[current] = score;
                    trace[current] = source;
                    if score > best_score {
                        best_score = score;
                        best_frame = frame;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }
    }

    if best_score == 0 {
        return ThreeFrameSwResult::default();
    }

    let mut frame = best_frame;
    let mut i = best_i;
    let mut j = best_j;
    let mut ops = Vec::new();
    let mut sw = SwResult {
        score: best_score,
        query_end: best_i as i32,
        subject_end: best_j as i32,
        ..Default::default()
    };

    while i > 0 && j > 0 && h[idx(frame, i, j)] > 0 {
        match trace[idx(frame, i, j)] {
            Trace::Diag => {
                let q = query[frame][i - 1] & LETTER_MASK;
                let s = subject[j - 1] & LETTER_MASK;
                if q == s {
                    sw.identities += 1;
                    ops.push((EditOperation::Match, 1));
                } else {
                    sw.mismatches += 1;
                    ops.push((EditOperation::Substitution, 1));
                }
                sw.length += 1;
                i -= 1;
                j -= 1;
            }
            Trace::ForwardShift => {
                let q = query[frame][i - 1] & LETTER_MASK;
                let s = subject[j - 1] & LETTER_MASK;
                if q == s {
                    sw.identities += 1;
                    ops.push((EditOperation::Match, 1));
                } else {
                    sw.mismatches += 1;
                    ops.push((EditOperation::Substitution, 1));
                }
                ops.push((EditOperation::FrameshiftForward, 1));
                sw.length += 1;
                j -= 1;
                if frame == 0 {
                    frame = 2;
                    i -= 2;
                } else {
                    frame -= 1;
                    i -= 1;
                }
            }
            Trace::ReverseShift => {
                let q = query[frame][i - 1] & LETTER_MASK;
                let s = subject[j - 1] & LETTER_MASK;
                if q == s {
                    sw.identities += 1;
                    ops.push((EditOperation::Match, 1));
                } else {
                    sw.mismatches += 1;
                    ops.push((EditOperation::Substitution, 1));
                }
                ops.push((EditOperation::FrameshiftReverse, 1));
                sw.length += 1;
                j -= 1;
                if frame == 2 {
                    frame = 0;
                } else {
                    frame += 1;
                    i -= 1;
                }
            }
            Trace::HGapOpen | Trace::HGapExtend => {
                let mut gap_len = 0;
                loop {
                    gap_len += 1;
                    let t = trace[idx(frame, i, j)];
                    j -= 1;
                    if t == Trace::HGapOpen || j == 0 {
                        break;
                    }
                }
                ops.push((EditOperation::Deletion, gap_len));
                sw.gap_openings += 1;
                sw.gaps += gap_len;
                sw.length += gap_len;
            }
            Trace::VGapOpen | Trace::VGapExtend => {
                let mut gap_len = 0;
                loop {
                    gap_len += 1;
                    let t = trace[idx(frame, i, j)];
                    i -= 1;
                    if t == Trace::VGapOpen || i == 0 {
                        break;
                    }
                }
                ops.push((EditOperation::Insertion, gap_len));
                sw.gap_openings += 1;
                sw.gaps += gap_len;
                sw.length += gap_len;
            }
            Trace::Zero => break,
        }
    }

    sw.query_begin = i as i32;
    sw.subject_begin = j as i32;
    ops.reverse();
    sw.operations = ops;

    ThreeFrameSwResult {
        sw,
        frame_begin: frame,
        frame_end: best_frame,
    }
}

pub fn banded_3frame_swipe_range(
    query: [&[Letter]; 3],
    subject: &[Letter],
    d_begin: i32,
    d_end: i32,
    score_matrix: &ScoreMatrix,
) -> ThreeFrameSwResult {
    let qlen = query.iter().map(|q| q.len()).max().unwrap_or(0);
    let slen = subject.len();
    if qlen == 0 || slen == 0 || d_end <= d_begin {
        return ThreeFrameSwResult::default();
    }

    let gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
    let gap_extend = score_matrix.gap_extend();
    let frame_shift = score_matrix.frame_shift();
    let neg_inf = i32::MIN / 4;
    let states = (qlen + 1) * (slen + 1) * 3;
    let idx = |frame: usize, i: usize, j: usize| -> usize { ((j * (qlen + 1) + i) * 3) + frame };

    let mut h = vec![0i32; states];
    let mut hgap = vec![neg_inf; states];
    let mut vgap = vec![neg_inf; states];
    let mut trace = vec![Trace::Zero; states];
    let mut best_score = 0i32;
    let mut best_frame = 0usize;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for j in 1..=slen {
        let subject_pos = j as i32 - 1;
        let lower = (d_begin + subject_pos).max(0);
        let upper = (d_end + subject_pos - 1).min(qlen as i32 - 1);
        if lower > upper {
            continue;
        }

        for i in (lower as usize + 1)..=(upper as usize + 1) {
            for frame in 0..3 {
                if i > query[frame].len() {
                    continue;
                }
                let q = query[frame][i - 1] & LETTER_MASK;
                let s = subject[j - 1] & LETTER_MASK;
                let match_score = score_matrix.score(q, s);
                let current = idx(frame, i, j);

                let mut score = h[idx(frame, i - 1, j - 1)] + match_score;
                let mut source = Trace::Diag;

                let fwd_pred_frame = if frame == 0 { 2 } else { frame - 1 };
                let fwd_pred_i = if frame == 0 {
                    i.checked_sub(2)
                } else {
                    i.checked_sub(1)
                };
                if let Some(pi) = fwd_pred_i {
                    if pi <= query[fwd_pred_frame].len() {
                        let candidate =
                            h[idx(fwd_pred_frame, pi, j - 1)] + match_score - frame_shift;
                        if candidate > score {
                            score = candidate;
                            source = Trace::ForwardShift;
                        }
                    }
                }

                let rev_pred_frame = if frame == 2 { 0 } else { frame + 1 };
                let rev_pred_i = if frame == 2 {
                    Some(i)
                } else {
                    i.checked_sub(1)
                };
                if let Some(pi) = rev_pred_i {
                    if pi <= query[rev_pred_frame].len() {
                        let candidate =
                            h[idx(rev_pred_frame, pi, j - 1)] + match_score - frame_shift;
                        if candidate > score {
                            score = candidate;
                            source = Trace::ReverseShift;
                        }
                    }
                }

                let h_open = h[idx(frame, i, j - 1)] - gap_open;
                let h_extend = hgap[idx(frame, i, j - 1)] - gap_extend;
                if h_open >= h_extend {
                    hgap[current] = h_open;
                    if h_open > score {
                        score = h_open;
                        source = Trace::HGapOpen;
                    }
                } else {
                    hgap[current] = h_extend;
                    if h_extend > score {
                        score = h_extend;
                        source = Trace::HGapExtend;
                    }
                }

                let v_open = h[idx(frame, i - 1, j)] - gap_open;
                let v_extend = vgap[idx(frame, i - 1, j)] - gap_extend;
                if v_open >= v_extend {
                    vgap[current] = v_open;
                    if v_open > score {
                        score = v_open;
                        source = Trace::VGapOpen;
                    }
                } else {
                    vgap[current] = v_extend;
                    if v_extend > score {
                        score = v_extend;
                        source = Trace::VGapExtend;
                    }
                }

                if score <= 0 {
                    h[current] = 0;
                    trace[current] = Trace::Zero;
                } else {
                    h[current] = score;
                    trace[current] = source;
                    if score > best_score {
                        best_score = score;
                        best_frame = frame;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }
    }

    traceback_3frame(
        query, subject, &h, &trace, qlen, best_score, best_frame, best_i, best_j,
    )
}

pub fn traceback_3frame(
    query: [&[Letter]; 3],
    subject: &[Letter],
    h: &[i32],
    trace: &[Trace],
    qlen: usize,
    best_score: i32,
    best_frame: usize,
    best_i: usize,
    best_j: usize,
) -> ThreeFrameSwResult {
    if best_score == 0 {
        return ThreeFrameSwResult::default();
    }
    let idx = |frame: usize, i: usize, j: usize| -> usize { ((j * (qlen + 1) + i) * 3) + frame };
    let mut frame = best_frame;
    let mut i = best_i;
    let mut j = best_j;
    let mut ops = Vec::new();
    let mut sw = SwResult {
        score: best_score,
        query_end: best_i as i32,
        subject_end: best_j as i32,
        ..Default::default()
    };

    while i > 0 && j > 0 && h[idx(frame, i, j)] > 0 {
        match trace[idx(frame, i, j)] {
            Trace::Diag => {
                let q = query[frame][i - 1] & LETTER_MASK;
                let s = subject[j - 1] & LETTER_MASK;
                if q == s {
                    sw.identities += 1;
                    ops.push((EditOperation::Match, 1));
                } else {
                    sw.mismatches += 1;
                    ops.push((EditOperation::Substitution, 1));
                }
                sw.length += 1;
                i -= 1;
                j -= 1;
            }
            Trace::ForwardShift => {
                let q = query[frame][i - 1] & LETTER_MASK;
                let s = subject[j - 1] & LETTER_MASK;
                if q == s {
                    sw.identities += 1;
                    ops.push((EditOperation::Match, 1));
                } else {
                    sw.mismatches += 1;
                    ops.push((EditOperation::Substitution, 1));
                }
                ops.push((EditOperation::FrameshiftForward, 1));
                sw.length += 1;
                j -= 1;
                if frame == 0 {
                    frame = 2;
                    i -= 2;
                } else {
                    frame -= 1;
                    i -= 1;
                }
            }
            Trace::ReverseShift => {
                let q = query[frame][i - 1] & LETTER_MASK;
                let s = subject[j - 1] & LETTER_MASK;
                if q == s {
                    sw.identities += 1;
                    ops.push((EditOperation::Match, 1));
                } else {
                    sw.mismatches += 1;
                    ops.push((EditOperation::Substitution, 1));
                }
                ops.push((EditOperation::FrameshiftReverse, 1));
                sw.length += 1;
                j -= 1;
                if frame == 2 {
                    frame = 0;
                } else {
                    frame += 1;
                    i -= 1;
                }
            }
            Trace::HGapOpen | Trace::HGapExtend => {
                let mut gap_len = 0;
                loop {
                    gap_len += 1;
                    let t = trace[idx(frame, i, j)];
                    j -= 1;
                    if t == Trace::HGapOpen || j == 0 {
                        break;
                    }
                }
                ops.push((EditOperation::Deletion, gap_len));
                sw.gap_openings += 1;
                sw.gaps += gap_len;
                sw.length += gap_len;
            }
            Trace::VGapOpen | Trace::VGapExtend => {
                let mut gap_len = 0;
                loop {
                    gap_len += 1;
                    let t = trace[idx(frame, i, j)];
                    i -= 1;
                    if t == Trace::VGapOpen || i == 0 {
                        break;
                    }
                }
                ops.push((EditOperation::Insertion, gap_len));
                sw.gap_openings += 1;
                sw.gaps += gap_len;
                sw.length += gap_len;
            }
            Trace::Zero => break,
        }
    }

    sw.query_begin = i as i32;
    sw.subject_begin = j as i32;
    ops.reverse();
    sw.operations = ops;

    ThreeFrameSwResult {
        sw,
        frame_begin: frame,
        frame_end: best_frame,
    }
}

pub fn banded_3frame_swipe_targets(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    score_only: bool,
    score_matrix: &ScoreMatrix,
    stat: &mut DpStat,
    _parallel: bool,
    overflow: &mut Vec<DpTarget>,
) -> Vec<Hsp> {
    let mut out = Vec::new();
    for target in targets {
        let band = target.band().max(0) as usize;
        stat.gross_cells += band * target.cols.max(0) as usize * 3;
        if target.seq.is_empty() || target.d_end <= target.d_begin {
            continue;
        }
        let result = banded_3frame_swipe_range(
            query,
            &target.seq,
            target.d_begin,
            target.d_end,
            score_matrix,
        );
        if result.sw.score >= i16::MAX as i32 {
            overflow.push(target.clone());
            continue;
        }
        if result.sw.score <= 0 {
            continue;
        }
        stat.net_cells += result.sw.length.max(0) as usize;
        out.push(traceback_to_hsp(result, target, score_only, score_matrix));
    }
    out
}

pub fn banded_3frame_swipe_worker(
    query: [&[Letter]; 3],
    targets: &[DpTarget],
    score_only: bool,
    score_matrix: &ScoreMatrix,
    out: &mut Vec<Hsp>,
    overflow: &mut Vec<DpTarget>,
) {
    let mut stat = DpStat::default();
    out.extend(banded_3frame_swipe_targets(
        query,
        targets,
        score_only,
        score_matrix,
        &mut stat,
        true,
        overflow,
    ));
}

pub fn banded_3frame_swipe_target_range(
    query: [&[Letter]; 3],
    targets: &mut [DpTarget],
    score_matrix: &ScoreMatrix,
    stat: &mut DpStat,
    score_only: bool,
    parallel: bool,
) -> Vec<Hsp> {
    targets.sort_by_key(|target| (target.left_i1(), target.band(), target.cols));
    let mut overflow16 = Vec::new();
    let mut overflow32 = Vec::new();
    let mut out = banded_3frame_swipe_targets(
        query,
        targets,
        score_only,
        score_matrix,
        stat,
        parallel,
        &mut overflow16,
    );
    out.extend(banded_3frame_swipe_targets(
        query,
        &overflow16,
        score_only,
        score_matrix,
        stat,
        false,
        &mut overflow32,
    ));
    out
}

fn traceback_to_hsp(
    result: ThreeFrameSwResult,
    target: &DpTarget,
    score_only: bool,
    score_matrix: &ScoreMatrix,
) -> Hsp {
    let mut hsp = Hsp::new();
    hsp.backtraced = !score_only;
    hsp.score = result.sw.score;
    hsp.bit_score = score_matrix.bitscore(result.sw.score as f64);
    hsp.evalue = score_matrix.evalue(
        result.sw.score,
        result.sw.query_end.max(1) as u32,
        target.seq.len() as u32,
    );
    hsp.frame = result.frame_end as i32;
    hsp.length = result.sw.length;
    hsp.identities = result.sw.identities;
    hsp.mismatches = result.sw.mismatches;
    hsp.gap_openings = result.sw.gap_openings;
    hsp.gaps = result.sw.gaps;
    hsp.query_range = Interval::new(result.sw.query_begin, result.sw.query_end);
    hsp.subject_range = Interval::new(result.sw.subject_begin, result.sw.subject_end);
    hsp.query_source_range = hsp.query_range;
    hsp.subject_source_range = hsp.subject_range;
    hsp.target_seq = target.seq.clone();
    hsp.swipe_target = target.target_idx as i32;
    hsp.d_begin = target.d_begin;
    hsp.d_end = target.d_end;
    if !score_only {
        for (op, len) in result.sw.operations {
            for _ in 0..len {
                hsp.transcript.push(op);
            }
        }
    }
    hsp
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dp::smith_waterman::smith_waterman;
    use crate::dp::swipe::{Anchor, CarryOver, DpTarget};

    #[test]
    fn test_banded_3frame_matches_single_frame_sw() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 1000, 1, 0).unwrap();
        let q0: Vec<Letter> = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let q1: Vec<Letter> = vec![13, 13, 13, 13, 13, 13, 13];
        let q2: Vec<Letter> = vec![17, 17, 17, 17, 17, 17];
        let result = banded_3frame_swipe([&q0, &q1, &q2], &q0, 4, 4, 10, &sm);
        let sw = smith_waterman(&q0, &q0, &sm);
        assert_eq!(result.sw.score, sw.score);
        assert_eq!(result.frame_end, 0);
        assert_eq!(result.sw.identities, q0.len() as i32);
    }

    #[test]
    fn test_banded_3frame_can_use_frameshift() {
        let sm_shift = ScoreMatrix::new("blosum62", 11, 1, 1, 1, 0).unwrap();
        let sm_no_shift = ScoreMatrix::new("blosum62", 11, 1, 1000, 1, 0).unwrap();
        let q0: Vec<Letter> = vec![0, 1, 2, 3];
        let q1: Vec<Letter> = vec![0, 1, 4, 3];
        let q2: Vec<Letter> = vec![0, 1, 2, 3];
        let subject: Vec<Letter> = vec![0, 1, 4, 3];
        let shifted = banded_3frame_swipe([&q0, &q1, &q2], &subject, 2, 2, 10, &sm_shift);
        let unshifted = banded_3frame_swipe([&q0, &q1, &q2], &subject, 2, 2, 10, &sm_no_shift);
        assert!(shifted.sw.score >= unshifted.sw.score);
        assert!(shifted.sw.score > 0);
    }

    #[test]
    fn test_banded_3frame_range_matches_single_frame() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 1000, 1, 0).unwrap();
        let q0: Vec<Letter> = vec![0, 1, 2, 3, 4, 5];
        let q1: Vec<Letter> = vec![13; 5];
        let q2: Vec<Letter> = vec![17; 5];
        let result = banded_3frame_swipe_range([&q0, &q1, &q2], &q0, -1, 2, &sm);
        let sw = smith_waterman(&q0, &q0, &sm);
        assert_eq!(result.sw.score, sw.score);
        assert_eq!(result.sw.query_begin, 0);
        assert_eq!(result.sw.subject_begin, 0);
    }

    #[test]
    fn test_banded_3frame_swipe_targets_outputs_hsp() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 1000, 1, 0).unwrap();
        let q0: Vec<Letter> = vec![0, 1, 2, 3, 4, 5];
        let q1: Vec<Letter> = vec![13; 5];
        let q2: Vec<Letter> = vec![17; 5];
        let target = DpTarget::new(
            q0.clone(),
            q0.len() as i32,
            -1,
            2,
            42,
            q0.len() as i32,
            CarryOver::default(),
            Anchor::default(),
        );
        let mut stat = DpStat::default();
        let mut overflow = Vec::new();
        let out = banded_3frame_swipe_targets(
            [&q0, &q1, &q2],
            &[target],
            false,
            &sm,
            &mut stat,
            false,
            &mut overflow,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].swipe_target, 42);
        assert!(out[0].backtraced);
        assert!(overflow.is_empty());
        assert!(stat.gross_cells > 0);
    }

    #[test]
    fn test_banded_3frame_matrix_iterators() {
        let mut matrix = Banded3FrameSwipeMatrix::new(6, 2);
        assert_eq!(matrix.band(), 6);
        let mut it = matrix.begin(0, 0);
        it.set_hgap(7);
        it.set_score(11);
        assert_eq!(it.sm3, 0);
        it.next();

        let mut traceback = Banded3FrameSwipeTracebackMatrix::new(6, 2);
        assert_eq!(traceback.band(), 6);
        let mut col = traceback.begin(0, 0);
        col.set_hgap(3);
        col.set_score(5);
        col.next();
        let mut tb = traceback.traceback(1, 0, 0, 20, 0, 5);
        tb.walk_diagonal();
        tb.walk_forward_shift();
        tb.walk_reverse_shift();
        let _ = tb.walk_gap(-1, 2, EditOperation::Deletion, 1);
    }
}
