//! Translation of `diamond/src/align/full_db.cpp`.

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use crate::align::target::{GappedScoreConfig, Target};
use crate::basic::statistics::Statistics;
use crate::basic::value::{BlockId, Letter};
use crate::data::block::Block;
use crate::dp::swipe::{swipe_set, Flags, HspValues, Params};
use crate::stats::score_matrix::ScoreMatrix;

/// Align every query context against every sequence in a target block and
/// group the resulting HSPs by subject, matching C++ `full_db_align`.
pub fn full_db_align(
    query_seq: &[Vec<Letter>],
    query_cbs: &[Vec<i8>],
    flags: Flags,
    hsp_values: HspValues,
    stat: &mut Statistics,
    target_block: &Block,
    cfg: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
) -> Vec<Target> {
    let swipe_statistics = Arc::new(Mutex::new(Statistics::new()));
    let ref_seqs = target_block.seqs();
    let mut hsps = Vec::new();

    for frame in 0..cfg.query_contexts {
        let composition_bias = if cfg.comp_based_stats_hauser {
            query_cbs.get(frame).map(Vec::as_slice)
        } else {
            None
        };
        let mut params = Params::new(&query_seq[frame], score_matrix);
        params.query_id = Some("");
        params.frame = frame as i32;
        params.query_source_len = query_seq[frame].len() as i32;
        params.composition_bias = composition_bias;
        params.flags = flags | Flags::FULL_MATRIX;
        params.target_max_len = ref_seqs.max_len(0, ref_seqs.len() as BlockId);
        params.v = hsp_values;
        params.cutoff_score_8bit = cfg.cutoff_score_8bit;
        params.max_swipe_dp = cfg.max_swipe_dp;
        params.approx_backtrace = cfg.approx_backtrace;
        params.max_evalue = cfg.max_evalue;
        params.query_cover = cfg.query_cover;
        params.subject_cover = cfg.subject_cover;
        params.query_or_target_cover = cfg.query_or_target_cover;
        params.approx_min_id = cfg.approx_min_id;
        params.cbs_matrix_scale = cfg.cbs_matrix_scale;
        params.statistics = Some(swipe_statistics.clone());
        // `std::list::splice(hsp.begin(), frame_hsp, ...)` prepends every
        // context in the original, so later frames are visited first.
        let mut frame_hsps = swipe_set(ref_seqs, &mut params);
        frame_hsps.append(&mut hsps);
        hsps = frame_hsps;
    }
    *stat += &swipe_statistics.lock().unwrap();

    let mut subject_idx = BTreeMap::new();
    let mut targets = Vec::new();
    for hsp in hsps {
        let block_id = hsp.swipe_target as BlockId;
        let idx = if let Some(&idx) = subject_idx.get(&block_id) {
            idx
        } else {
            let idx = targets.len();
            subject_idx.insert(block_id, idx);
            targets.push(Target::new(
                block_id,
                ref_seqs.get(block_id as usize),
                0,
                None,
            ));
            idx
        };
        let frame = hsp.frame as usize;
        let score = hsp.score;
        let evalue = hsp.evalue;
        let hsp_frame = hsp.frame;
        let target = &mut targets[idx];
        target.hsp[frame].push(hsp);
        if score > target.filter_score {
            target.filter_evalue = evalue;
            target.filter_score = score;
            target.best_context = hsp_frame;
        }
    }

    targets
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::statistics::StatValue;
    use crate::basic::value::SequenceType;

    #[test]
    fn full_db_alignment_propagates_swipe_statistics() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5, 6];
        let mut block = Block::new();
        block
            .push_back(&query, None, None, 0, SequenceType::AminoAcid, 0, false)
            .unwrap();
        let mut statistics = Statistics::new();
        let targets = full_db_align(
            std::slice::from_ref(&query),
            &[],
            Flags::NONE,
            HspValues::COORDS,
            &mut statistics,
            &block,
            &GappedScoreConfig::default(),
            &score_matrix,
        );
        assert_eq!(targets.len(), 1);
        assert!(statistics.get(StatValue::SwipeTasksTotal) > 0);
        assert!(statistics.get(StatValue::GrossDpCells) > 0);
    }
}
