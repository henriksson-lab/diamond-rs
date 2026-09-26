//! Final gapped alignment pass mirrored from `diamond/src/align/gapped_final.cpp`.
//!
//! Configuration that is process-global in C++ is passed explicitly through
//! `GappedScoreConfig`, the score matrix, and the requested output values.

use super::hsp::Match;
use super::target::{
    culling_matches, match_apply_filters, match_inner_culling, recompute_alt_hsps,
    GappedScoreConfig, Target,
};
use crate::basic::consts::MAX_CONTEXT;
use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::value::Letter;
use crate::dp::swipe::{
    bin as swipe_bin, swipe, targets as make_dp_targets, Anchor as DpAnchor, CarryOver, DpTarget,
    Flags, HspValues, Params, Targets,
};
use crate::search::sensitivity::ExtensionMode;
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;

/// Matches C++ `filter_hspvalues`.
pub fn filter_hspvalues(cfg: &GappedScoreConfig) -> HspValues {
    let mut hsp_values = HspValues::NONE;
    if cfg.max_hsps != 1 {
        hsp_values = hsp_values | HspValues::QUERY_COORDS | HspValues::TARGET_COORDS;
    }
    if cfg.min_id > 0.0 {
        hsp_values = hsp_values | HspValues::IDENT | HspValues::LENGTH;
    }
    if cfg.approx_min_id > 0.0 {
        hsp_values = hsp_values | HspValues::COORDS;
    }
    if cfg.query_cover > 0.0 {
        hsp_values = hsp_values | HspValues::QUERY_COORDS;
    }
    if cfg.subject_cover > 0.0 {
        hsp_values = hsp_values | HspValues::TARGET_COORDS;
    }
    if cfg.query_or_target_cover > 0.0 {
        hsp_values = hsp_values | HspValues::COORDS;
    }
    hsp_values
}

/// Matches C++ `first_round_filter_all`.
pub fn first_round_filter_all(cfg: &GappedScoreConfig, first_round_hsp_values: HspValues) -> bool {
    if cfg.min_id > 0.0 && !first_round_hsp_values.all(HspValues::IDENT | HspValues::LENGTH) {
        return false;
    }
    if cfg.approx_min_id > 0.0 && !first_round_hsp_values.all(HspValues::COORDS) {
        return false;
    }
    if cfg.query_cover > 0.0 && !first_round_hsp_values.all(HspValues::QUERY_COORDS) {
        return false;
    }
    if cfg.subject_cover > 0.0 && !first_round_hsp_values.all(HspValues::TARGET_COORDS) {
        return false;
    }
    if cfg.query_or_target_cover > 0.0 && !first_round_hsp_values.all(HspValues::COORDS) {
        return false;
    }
    true
}

/// Matches C++ `add_dp_targets` in `gapped_final.cpp`.
pub fn add_final_dp_targets(
    target: &Target,
    target_idx: usize,
    query_seq: &[Vec<Letter>],
    dp_targets: &mut [Targets; MAX_CONTEXT as usize],
    flags: Flags,
    hsp_values: HspValues,
    cfg: &GappedScoreConfig,
) {
    assert!(
        cfg.query_contexts <= query_seq.len(),
        "final alignment query-context count exceeds query sequences"
    );
    let tlen = target.seq.len() as i32;
    let score_width = target
        .matrix
        .as_deref()
        .map(TargetMatrix::score_width)
        .unwrap_or(0);
    for frame in 0..cfg.query_contexts {
        let qlen = query_seq[frame].len() as i32;
        for hsp in &target.hsp[frame] {
            let dp_size = if flags.any(Flags::FULL_MATRIX) {
                qlen as i64 * tlen as i64
            } else {
                DpTarget::banded_cols(qlen, tlen, hsp.d_begin, hsp.d_end) as i64
                    * (hsp.d_end - hsp.d_begin) as i64
            };
            let bin = swipe_bin(
                hsp_values,
                if flags.any(Flags::FULL_MATRIX) {
                    qlen
                } else {
                    hsp.d_end - hsp.d_begin
                },
                hsp.score,
                0,
                dp_size,
                score_width,
                0,
                cfg.cutoff_score_8bit,
                cfg.max_swipe_dp,
                cfg.approx_backtrace,
            );
            let mut dp = DpTarget::new(
                target.seq.clone(),
                tlen,
                hsp.d_begin,
                hsp.d_end,
                target_idx as i64,
                qlen,
                CarryOver::default(),
                DpAnchor::default(),
            );
            if let Some(matrix) = target.matrix.clone() {
                dp = dp.with_matrix(matrix, cfg.cbs_matrix_scale);
            }
            dp_targets[frame][bin].push_back(dp);
        }
    }
}

/// Matches C++ `align(vector<Target>&, ...)` in `gapped_final.cpp`.
#[allow(clippy::too_many_arguments)]
pub fn align_targets_final(
    targets: &mut Vec<Target>,
    previous_matches: i64,
    query_seq: &[Vec<Letter>],
    query_id: &str,
    query_cbs: &[Vec<i8>],
    source_query_len: i32,
    _query_self_aln_score: f64,
    mut flags: Flags,
    first_round: HspValues,
    first_round_culling: bool,
    mode: ExtensionMode,
    stat: &mut Statistics,
    cfg: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
    output_hsp_values: HspValues,
) -> Vec<Match> {
    const MIN_STEP: i64 = 16;
    let mut matches = Vec::new();
    if targets.is_empty() {
        return matches;
    }
    assert!(
        cfg.query_contexts <= query_seq.len(),
        "final alignment query-context count exceeds query sequences"
    );
    if cfg.comp_based_stats_hauser {
        assert!(
            cfg.query_contexts <= query_cbs.len(),
            "final alignment CBS context count mismatch"
        );
    }

    let mut hsp_values = output_hsp_values;
    let copy_all = cfg.max_hsps == 1
        && first_round.all(hsp_values)
        && first_round_filter_all(cfg, first_round);
    if copy_all {
        matches.reserve(targets.len());
    }
    for target in targets.iter_mut() {
        if copy_all || target.done {
            let matrix = target.matrix.take();
            matches.push(Match::from_target_hsps(
                target.block_id,
                &target.seq,
                matrix,
                &mut target.hsp,
                target.ungapped_score,
                cfg.query_contexts,
                cfg.max_hsps,
            ));
        }
    }
    if matches.len() == targets.len() {
        for result in &mut matches {
            let subject_seq = result.seq.clone();
            match_apply_filters(
                result,
                source_query_len as u32,
                query_id,
                &query_seq[0],
                subject_seq.len() as u32,
                None,
                &subject_seq,
                cfg.min_id,
                cfg.approx_min_id,
                cfg.query_cover,
                cfg.subject_cover,
                cfg.query_or_target_cover,
                cfg.no_self_hits,
            );
        }
        return matches;
    }

    if mode == ExtensionMode::Full {
        flags = flags | Flags::FULL_MATRIX;
    }
    if mode == ExtensionMode::Global {
        flags = flags | Flags::SEMI_GLOBAL;
    }
    hsp_values = hsp_values | filter_hspvalues(cfg);

    let mut target_begin = 0usize;
    let continue_alignment = |matches_len: usize| {
        cfg.toppercent.is_some() || (matches_len as i64 + previous_matches) < cfg.max_target_seqs
    };

    while target_begin < targets.len() && continue_alignment(matches.len()) {
        let mut dp_targets: [Targets; MAX_CONTEXT as usize] =
            std::array::from_fn(|_| make_dp_targets());
        let remaining = targets.len() - target_begin;
        let step_size = if !first_round_culling && cfg.toppercent.is_none() {
            let wanted = (cfg.max_target_seqs - matches.len() as i64).max(MIN_STEP);
            let rounded = ((wanted + MIN_STEP - 1) / MIN_STEP) * MIN_STEP;
            rounded.min(remaining as i64) as usize
        } else {
            remaining
        };
        matches.reserve(step_size);
        let matches_begin = matches.len();

        for target in targets.iter_mut().skip(target_begin).take(step_size) {
            if target.done {
                continue;
            }
            add_final_dp_targets(
                target,
                matches.len(),
                query_seq,
                &mut dp_targets,
                flags,
                hsp_values,
                cfg,
            );
            let matrix = target.matrix.take();
            matches.push(Match::new_extension(
                target.block_id,
                &target.seq,
                matrix,
                target.ungapped_score,
                0,
                f64::MAX,
            ));
        }

        for frame in 0..cfg.query_contexts {
            if dp_targets[frame].iter().all(|bin| bin.size() == 0) {
                continue;
            }
            let composition_bias = cfg
                .comp_based_stats_hauser
                .then(|| query_cbs[frame].as_slice());
            let mut params = Params::new(&query_seq[frame], score_matrix);
            params.query_id = Some(query_id);
            params.frame = frame as i32;
            params.query_source_len = source_query_len;
            params.composition_bias = composition_bias;
            params.flags = flags;
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
            for hsp in swipe(&dp_targets[frame], &mut params) {
                let result = &mut matches[hsp.swipe_target as usize];
                if hsp.score > result.filter_score {
                    result.filter_evalue = hsp.evalue;
                    result.filter_score = hsp.score;
                }
                result.hsps.push(hsp);
            }
        }

        for result in matches.iter_mut().skip(matches_begin) {
            match_inner_culling(result, cfg.max_hsps, cfg.inner_culling_overlap);
            if let Some(best) = result.hsps.first() {
                result.filter_score = best.score;
                result.filter_evalue = best.evalue;
            }
            let subject_seq = result.seq.clone();
            match_apply_filters(
                result,
                source_query_len as u32,
                query_id,
                &query_seq[0],
                subject_seq.len() as u32,
                None,
                &subject_seq,
                cfg.min_id,
                cfg.approx_min_id,
                cfg.query_cover,
                cfg.subject_cover,
                cfg.query_or_target_cover,
                cfg.no_self_hits,
            );
        }
        culling_matches(&mut matches, cfg.max_target_seqs, cfg.toppercent, |score| {
            score_matrix.bitscore(score as f64)
        });
        stat.inc(StatValue::TargetHits6, step_size as i64);
        target_begin += step_size;
    }

    recompute_alt_hsps(
        &mut matches,
        query_seq,
        query_cbs,
        source_query_len,
        hsp_values,
        stat,
        cfg,
        score_matrix,
    );
    matches
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::util::interval::Interval;

    fn score_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    fn finished_target(seq: &[Letter]) -> Target {
        let mut target = Target::new(7, seq, 40, None);
        let mut hsp = Hsp::new();
        hsp.frame = 0;
        hsp.score = 40;
        hsp.evalue = 1.0e-10;
        hsp.identities = seq.len() as i32;
        hsp.length = seq.len() as i32;
        hsp.query_range = Interval::new(0, seq.len() as i32);
        hsp.query_source_range = hsp.query_range;
        hsp.subject_range = hsp.query_range;
        target.add_hit(hsp);
        target
    }

    #[test]
    fn filter_values_cover_each_upstream_requirement() {
        let mut cfg = GappedScoreConfig::default();
        cfg.max_hsps = 2;
        assert!(filter_hspvalues(&cfg).all(HspValues::COORDS));

        cfg = GappedScoreConfig::default();
        cfg.min_id = 1.0;
        assert_eq!(filter_hspvalues(&cfg), HspValues::IDENT | HspValues::LENGTH);
        assert!(!first_round_filter_all(&cfg, HspValues::IDENT));
        assert!(first_round_filter_all(
            &cfg,
            HspValues::IDENT | HspValues::LENGTH
        ));

        cfg = GappedScoreConfig::default();
        cfg.query_or_target_cover = 1.0;
        assert_eq!(filter_hspvalues(&cfg), HspValues::COORDS);
        assert!(!first_round_filter_all(&cfg, HspValues::QUERY_COORDS));
    }

    #[test]
    fn copy_all_reuses_first_round_hsp_without_running_dp() {
        let matrix = score_matrix();
        let query = vec![0, 1, 2, 3];
        let mut target = finished_target(&query);
        target.matrix = Some(std::sync::Arc::new(TargetMatrix::new(
            vec![0; 32 * crate::basic::value::AMINO_ACID_COUNT],
            0,
            0,
        )));
        let mut targets = vec![target];
        let cfg = GappedScoreConfig {
            query_contexts: 1,
            max_hsps: 1,
            max_target_seqs: 10,
            ..GappedScoreConfig::default()
        };
        let mut stat = Statistics::new();
        let matches = align_targets_final(
            &mut targets,
            0,
            std::slice::from_ref(&query),
            "query",
            &[],
            query.len() as i32,
            0.0,
            Flags::NONE,
            HspValues::COORDS,
            false,
            ExtensionMode::BandedFast,
            &mut stat,
            &cfg,
            &matrix,
            HspValues::COORDS,
        );
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].hsps.len(), 1);
        assert_eq!(matches[0].filter_score, 40);
        assert!(matches[0].matrix.is_some());
        assert!(targets[0].matrix.is_none());
        assert!(targets[0].hsp[0].is_empty());
        assert_eq!(stat.get(StatValue::TargetHits6), 0);
    }

    #[test]
    #[should_panic(expected = "final alignment CBS context count mismatch")]
    fn hauser_mode_requires_one_cbs_vector_per_context() {
        let matrix = score_matrix();
        let query = vec![0, 1, 2, 3];
        let mut targets = vec![finished_target(&query)];
        let cfg = GappedScoreConfig {
            query_contexts: 1,
            comp_based_stats_hauser: true,
            ..GappedScoreConfig::default()
        };
        let _ = align_targets_final(
            &mut targets,
            0,
            std::slice::from_ref(&query),
            "query",
            &[],
            query.len() as i32,
            0.0,
            Flags::NONE,
            HspValues::NONE,
            false,
            ExtensionMode::BandedFast,
            &mut Statistics::new(),
            &cfg,
            &matrix,
            HspValues::NONE,
        );
    }
}
