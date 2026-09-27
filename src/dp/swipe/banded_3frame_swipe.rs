//! Exact-path facade for `diamond/src/dp/swipe/banded_3frame_swipe.cpp`.
//!
//! The scalar architecture-neutral kernel remains in `dp::banded_3frame` for
//! API compatibility; this module supplies the original translated-sequence
//! dispatch boundary and re-exports its matrix/worker surface.

#[path = "banded_3frame_simd.rs"]
pub mod simd;

pub use crate::dp::banded_3frame::{
    banded_3frame_swipe_range, banded_3frame_swipe_targets, banded_3frame_swipe_worker,
    Banded3FrameSwipeMatrix, Banded3FrameSwipeTracebackMatrix, ColumnIterator, DpStat,
    ThreeFrameSwResult, Trace, TracebackColumnIterator, TracebackIterator,
};

use crate::align::hsp::Hsp;
use crate::basic::translate::{Frame, Strand, TranslatedPosition};
use crate::data::sequence_set::TranslatedSequenceView;
use crate::dp::banded_3frame::banded_3frame_swipe_target_range;
use crate::dp::swipe::DpTarget;
use crate::stats::score_matrix::ScoreMatrix;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Banded3FrameSwipeConfig {
    pub score_only: bool,
    pub parallel: bool,
    pub cbs_matrix_scale: i32,
    pub min_bit_score: f64,
    pub max_evalue: f64,
}

impl Default for Banded3FrameSwipeConfig {
    fn default() -> Self {
        Self {
            score_only: false,
            parallel: false,
            cbs_matrix_scale: 1,
            min_bit_score: 0.0,
            max_evalue: f64::MAX,
        }
    }
}

/// C++ public `banded_3frame_swipe` dispatch wrapper.
pub fn banded_3frame_swipe(
    query: TranslatedSequenceView<'_>,
    strand: Strand,
    targets: &mut [DpTarget],
    stat: &mut DpStat,
    score_matrix: &ScoreMatrix,
    config: Banded3FrameSwipeConfig,
) -> Vec<Hsp> {
    let frame_base = if strand == Strand::Forward { 0 } else { 3 };
    let frames = [
        query.frame(frame_base),
        query.frame(frame_base + 1),
        query.frame(frame_base + 2),
    ];
    let source_len = query.source().len() as i32;
    let mut hsps = banded_3frame_swipe_target_range(
        frames,
        targets,
        score_matrix,
        stat,
        config.score_only,
        config.parallel,
    );
    hsps.retain_mut(|hsp| {
        hsp.score = hsp.score.saturating_mul(config.cbs_matrix_scale);
        hsp.bit_score = score_matrix.bitscore(hsp.score as f64);
        hsp.evalue = score_matrix.evalue(
            hsp.score,
            frames[0].len() as u32,
            hsp.target_seq.len() as u32,
        );
        let offset = if config.score_only { 0 } else { hsp.frame % 3 };
        let frame = Frame::new(strand, offset);
        hsp.frame = frame.index();
        hsp.query_source_range = TranslatedPosition::absolute_interval(
            TranslatedPosition::new(hsp.query_range.begin, frame),
            TranslatedPosition::new(hsp.query_range.end, frame),
            source_len,
            true,
        );
        score_matrix.report_cutoff(
            hsp.score,
            hsp.evalue,
            config.min_bit_score,
            config.max_evalue,
        )
    });
    hsps
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::Letter;
    use crate::data::sequence_set::SequenceSet;
    use crate::dp::swipe::{Anchor, CarryOver};

    fn translated_sequence() -> SequenceSet {
        // Translation is injected directly so the test isolates strand/frame dispatch.
        let mut sequences = SequenceSet::new();
        for frame in 0..6 {
            sequences.push(&vec![frame as Letter; 6]);
        }
        sequences
    }

    #[test]
    fn wrapper_selects_reverse_frames_and_scales_scores() {
        let sequences = translated_sequence();
        let source = vec![0 as Letter; 21];
        let query = sequences.translated_seq(&source, 0, true);
        let subject = vec![3 as Letter; 6];
        let mut targets = vec![DpTarget::new(
            subject,
            6,
            -1,
            2,
            7,
            6,
            CarryOver::default(),
            Anchor::default(),
        )];
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 1000, 1, 1000).unwrap();
        let mut stat = DpStat::default();
        let hsps = banded_3frame_swipe(
            query,
            Strand::Reverse,
            &mut targets,
            &mut stat,
            &matrix,
            Banded3FrameSwipeConfig {
                cbs_matrix_scale: 2,
                ..Default::default()
            },
        );
        assert_eq!(hsps.len(), 1);
        assert!(hsps[0].score > 0 && hsps[0].score % 2 == 0);
        assert_eq!(Frame::from_index(hsps[0].frame).strand, Strand::Reverse);
        assert!(hsps[0].query_source_range.begin >= 0);
        assert!(hsps[0].query_source_range.end <= source.len() as i32);
    }

    #[test]
    fn wrapper_applies_report_cutoff_after_scaling() {
        let sequences = translated_sequence();
        let source = vec![0 as Letter; 21];
        let query = sequences.translated_seq(&source, 0, true);
        let mut targets = vec![DpTarget::new(
            vec![0 as Letter; 6],
            6,
            -1,
            2,
            7,
            6,
            CarryOver::default(),
            Anchor::default(),
        )];
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 1000, 1, 1000).unwrap();
        let mut stat = DpStat::default();
        let hsps = banded_3frame_swipe(
            query,
            Strand::Forward,
            &mut targets,
            &mut stat,
            &matrix,
            Banded3FrameSwipeConfig {
                min_bit_score: f64::MAX,
                ..Default::default()
            },
        );
        assert!(hsps.is_empty());
    }
}
