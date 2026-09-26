//! Alternative-HSP recomputation from `diamond/src/align/alt_hsp.cpp`.

use crate::align::hsp::Match;
use crate::align::target::{match_inner_culling, GappedScoreConfig};
use crate::basic::consts::MAX_CONTEXT;
use crate::basic::statistics::Statistics;
use crate::basic::value::{Letter, SUPER_HARD_MASK};
use crate::dp::swipe::{
    bin as swipe_bin, swipe, targets as make_dp_targets, CarryOver, DpTarget, Flags, HspValues,
    Params, Targets,
};
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence::is_fully_masked;
use std::sync::{Arc, Mutex};

#[derive(Debug, Clone)]
struct ActiveTarget {
    match_index: usize,
    masked_seq: [Option<Vec<Letter>>; MAX_CONTEXT as usize],
    active: u32,
}

impl ActiveTarget {
    /// C++ `ActiveTarget::ActiveTarget(vector<Match>::iterator, SequenceSet&)`.
    fn new(match_index: usize, target_match: &Match, query_contexts: usize) -> Self {
        let mut masked_seq: [Option<Vec<Letter>>; MAX_CONTEXT as usize] =
            std::array::from_fn(|_| None);
        let mut reserved = 0u32;
        for hsp in &target_match.hsps {
            let bit = 1u32 << hsp.frame;
            if reserved & bit == 0 && (hsp.frame as usize) < query_contexts {
                masked_seq[hsp.frame as usize] = Some(target_match.seq.clone());
                reserved |= bit;
            }
        }
        Self {
            match_index,
            masked_seq,
            active: 0,
        }
    }

    /// C++ `ActiveTarget` copy constructor.
    fn copy_active(&self, query_contexts: usize) -> Self {
        let mut copy = self.clone();
        for context in 0..query_contexts {
            if self.active & (1u32 << context) == 0 {
                copy.masked_seq[context] = None;
            }
        }
        copy.active = 0;
        copy
    }

    /// C++ `ActiveTarget::copy_seq`.
    fn copy_seq(&mut self, target_match: &Match) {
        for hsp in &target_match.hsps {
            let frame = hsp.frame as usize;
            if let Some(sequence) = &mut self.masked_seq[frame] {
                let begin = hsp.subject_range.begin.max(0) as usize;
                let end = hsp
                    .subject_range
                    .end
                    .max(hsp.subject_range.begin)
                    .min(sequence.len() as i32) as usize;
                sequence[begin..end].fill(SUPER_HARD_MASK);
            }
        }
    }

    /// C++ `ActiveTarget::masked`.
    fn masked(&self, context: usize) -> Option<&[Letter]> {
        self.masked_seq[context].as_deref()
    }

    /// C++ `ActiveTarget::check_fully_masked`.
    fn check_fully_masked(&mut self, query_contexts: usize) -> i32 {
        let mut count = 0;
        for context in 0..query_contexts {
            if self.active & (1u32 << context) != 0 {
                if self.masked(context).is_none_or(is_fully_masked) {
                    self.active &= !(1u32 << context);
                } else {
                    count += 1;
                }
            }
        }
        count
    }
}

type TargetVec = Vec<ActiveTarget>;

/// C++ file-local `recompute_alt_hsps` round overload.
fn recompute_alt_hsps_round(
    matches: &mut [Match],
    targets: &mut [ActiveTarget],
    query_seq: &[Vec<Letter>],
    query_cbs: &[Vec<i8>],
    query_source_len: i32,
    hsp_values: HspValues,
    config: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
    swipe_statistics: &Arc<Mutex<Statistics>>,
) -> TargetVec {
    let mut dp_targets: [Targets; MAX_CONTEXT as usize] =
        std::array::from_fn(|_| make_dp_targets());
    let query_len = query_seq[0].len() as i32;

    for (target_index, active_target) in targets.iter().enumerate() {
        let target_match = &matches[active_target.match_index];
        let dp_size = query_len as i64 * target_match.seq.len() as i64;
        let score_width = target_match
            .matrix
            .as_deref()
            .map(TargetMatrix::score_width)
            .unwrap_or(0);
        let bin = swipe_bin(
            hsp_values,
            query_len,
            0,
            0,
            dp_size,
            score_width,
            0,
            config.cutoff_score_8bit,
            config.max_swipe_dp,
            config.approx_backtrace,
        );
        for context in 0..config.query_contexts {
            if let Some(masked) = active_target.masked(context) {
                let mut dp_target = DpTarget::full(
                    masked.to_vec(),
                    masked.len() as i32,
                    target_index as i64,
                    CarryOver::default(),
                );
                if let Some(matrix) = target_match.matrix.clone() {
                    dp_target = dp_target.with_matrix(matrix, config.cbs_matrix_scale);
                }
                dp_targets[context][bin].push_back(dp_target);
            }
        }
    }

    for context in 0..config.query_contexts {
        let composition_bias = if config.comp_based_stats_hauser {
            query_cbs.get(context).map(Vec::as_slice)
        } else {
            None
        };
        let mut params = Params::new(&query_seq[context], score_matrix);
        params.query_id = Some("");
        params.frame = context as i32;
        params.query_source_len = query_source_len;
        params.composition_bias = composition_bias;
        params.flags = Flags::FULL_MATRIX;
        params.v = hsp_values;
        params.cutoff_score_8bit = config.cutoff_score_8bit;
        params.max_swipe_dp = config.max_swipe_dp;
        params.approx_backtrace = config.approx_backtrace;
        params.max_evalue = config.max_evalue;
        params.query_cover = config.query_cover;
        params.subject_cover = config.subject_cover;
        params.query_or_target_cover = config.query_or_target_cover;
        params.approx_min_id = config.approx_min_id;
        params.cbs_matrix_scale = config.cbs_matrix_scale;
        params.statistics = Some(swipe_statistics.clone());
        for hsp in swipe(&dp_targets[context], &mut params) {
            let active_index = hsp.swipe_target as usize;
            let match_index = targets[active_index].match_index;
            let begin = hsp.subject_range.begin.max(0) as usize;
            let end = hsp
                .subject_range
                .end
                .max(hsp.subject_range.begin)
                .min(matches[match_index].seq.len() as i32) as usize;
            matches[match_index].hsps.push(hsp);
            if let Some(masked) = &mut targets[active_index].masked_seq[context] {
                masked[begin..end].fill(SUPER_HARD_MASK);
            }
            targets[active_index].active |= 1u32 << context;
        }
    }

    let mut output = Vec::new();
    for active_target in targets {
        if active_target.active != 0 {
            let target_match = &mut matches[active_target.match_index];
            match_inner_culling(target_match, config.max_hsps, config.inner_culling_overlap);
            if let Some(best) = target_match.hsps.first() {
                target_match.filter_score = best.score;
                target_match.filter_evalue = best.evalue;
            }
            if active_target.check_fully_masked(config.query_contexts) > 0
                && (target_match.hsps.len() < config.max_hsps as usize || config.max_hsps == 0)
            {
                output.push(active_target.copy_active(config.query_contexts));
            }
        }
    }
    output
}

/// C++ public `recompute_alt_hsps` iterator-range overload.
pub fn recompute_alt_hsps(
    matches: &mut [Match],
    query_seq: &[Vec<Letter>],
    query_cbs: &[Vec<i8>],
    query_source_len: i32,
    hsp_values: HspValues,
    stats: &mut Statistics,
    config: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
) {
    if config.max_hsps == 1 {
        return;
    }
    let swipe_statistics = Arc::new(Mutex::new(Statistics::new()));
    let mut targets = Vec::with_capacity(matches.len());
    for (match_index, target_match) in matches.iter().enumerate() {
        let mut active_target = ActiveTarget::new(match_index, target_match, config.query_contexts);
        active_target.copy_seq(target_match);
        targets.push(active_target);
    }
    while !targets.is_empty() {
        targets = recompute_alt_hsps_round(
            matches,
            &mut targets,
            query_seq,
            query_cbs,
            query_source_len,
            hsp_values,
            config,
            score_matrix,
            &swipe_statistics,
        );
    }
    *stats += &swipe_statistics.lock().unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::util::interval::Interval;

    fn hsp(frame: i32, begin: i32, end: i32) -> Hsp {
        let mut hsp = Hsp::new();
        hsp.frame = frame;
        hsp.subject_range = Interval::new(begin, end);
        hsp
    }

    #[test]
    fn active_target_allocates_once_per_used_context_and_masks_all_hsps() {
        let mut target_match = Match::new_extension(7, &[0, 1, 2, 3, 4, 5], None, 0, 0, 1.0);
        target_match.hsps = vec![hsp(0, 0, 2), hsp(0, 4, 6), hsp(1, 2, 4)];
        let mut active = ActiveTarget::new(0, &target_match, 2);
        active.copy_seq(&target_match);
        assert_eq!(
            active.masked(0).unwrap(),
            &[
                SUPER_HARD_MASK,
                SUPER_HARD_MASK,
                2,
                3,
                SUPER_HARD_MASK,
                SUPER_HARD_MASK
            ]
        );
        assert_eq!(
            active.masked(1).unwrap(),
            &[0, 1, SUPER_HARD_MASK, SUPER_HARD_MASK, 4, 5]
        );
        assert!(active.masked(2).is_none());
    }

    #[test]
    fn fully_masked_contexts_are_deactivated_and_copy_keeps_only_active() {
        let mut target_match = Match::new_extension(7, &[0, 1], None, 0, 0, 1.0);
        target_match.hsps = vec![hsp(0, 0, 2), hsp(1, 0, 1)];
        let mut active = ActiveTarget::new(0, &target_match, 2);
        active.copy_seq(&target_match);
        active.active = 0b11;
        assert_eq!(active.check_fully_masked(2), 1);
        assert_eq!(active.active, 0b10);
        let copy = active.copy_active(2);
        assert!(copy.masked(0).is_none());
        assert_eq!(copy.masked(1).unwrap(), &[SUPER_HARD_MASK, 1]);
        assert_eq!(copy.active, 0);
    }

    #[test]
    fn max_hsps_one_is_exact_noop() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3];
        let mut target_match = Match::new_extension(5, &query, None, 0, 40, 1.0e-5);
        target_match.hsps.push(hsp(0, 0, 4));
        let before = target_match.hsps.clone();
        let mut statistics = Statistics::new();
        recompute_alt_hsps(
            std::slice::from_mut(&mut target_match),
            std::slice::from_ref(&query),
            &[],
            query.len() as i32,
            HspValues::COORDS,
            &mut statistics,
            &GappedScoreConfig::default(),
            &score_matrix,
        );
        assert_eq!(target_match.hsps.len(), before.len());
        assert_eq!(target_match.hsps[0].subject_range, before[0].subject_range);
    }

    #[test]
    fn recomputation_finds_unmasked_subject_copy() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3];
        let subject = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let mut target_match = Match::new_extension(5, &subject, None, 0, 40, 1.0e-5);
        let mut original = hsp(0, 0, 4);
        original.score = 40;
        original.evalue = 1.0e-5;
        original.query_range = Interval::new(0, 4);
        original.query_source_range = Interval::new(0, 4);
        target_match.hsps.push(original);
        let mut config = GappedScoreConfig::default();
        config.max_hsps = 2;
        let mut statistics = Statistics::new();
        recompute_alt_hsps(
            std::slice::from_mut(&mut target_match),
            std::slice::from_ref(&query),
            &[],
            query.len() as i32,
            HspValues::COORDS,
            &mut statistics,
            &config,
            &score_matrix,
        );
        assert!(target_match.hsps.len() >= 2);
        assert!(target_match
            .hsps
            .iter()
            .any(|hsp| hsp.subject_range.begin >= 4));
        assert!(statistics.get(crate::basic::statistics::StatValue::SwipeTasksTotal) > 0);
    }
}
