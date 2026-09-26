//! Alignment-result culling from `diamond/src/align/culling.cpp`.
//!
//! C++ reads process-global thresholds and score conversion state. Rust passes
//! those values explicitly while retaining compatibility functions used by the
//! existing alignment pipeline.

use std::cmp::Ordering;

use super::hsp::{Hsp, Match};
use super::target::Target;
use crate::basic::consts::MAX_CONTEXT;
use crate::basic::value::Letter;
use crate::data::block::Block;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HspFilterConfig {
    pub min_id: f64,
    pub approx_min_id: f64,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub no_self_hits: bool,
    /// Retains the outer C++ dispatch condition. MCL expression evaluation is
    /// unavailable in both the default upstream build and this Rust build.
    pub cluster_threshold_present: bool,
}

impl Default for HspFilterConfig {
    fn default() -> Self {
        Self {
            min_id: 0.0,
            approx_min_id: 0.0,
            query_cover: 0.0,
            subject_cover: 0.0,
            query_or_target_cover: 0.0,
            no_self_hits: false,
            cluster_threshold_present: false,
        }
    }
}

/// C++ `Extension::max_hsp_culling(list<Hsp>&)`.
pub fn max_hsp_culling(hsps: &mut Vec<Hsp>, max_hsps: u32) {
    if max_hsps > 0 && hsps.len() > max_hsps as usize {
        hsps.truncate(max_hsps as usize);
    }
}

/// C++ `Extension::inner_culling(list<Hsp>&)`.
pub fn inner_culling(hsps: &mut Vec<Hsp>, max_hsps: u32, inner_culling_overlap: f64) {
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(compare_hsps);
    if max_hsps == 1 {
        hsps.truncate(1);
        return;
    }
    let overlap = inner_culling_overlap / 100.0;
    let mut kept = Vec::with_capacity(hsps.len());
    for hsp in hsps.drain(..) {
        if !kept
            .iter()
            .any(|previous| hsp.is_enveloped_by(previous, overlap))
        {
            kept.push(hsp);
        }
    }
    *hsps = kept;
    max_hsp_culling(hsps, max_hsps);
}

/// Compatibility name retained from the earlier Rust implementation.
pub fn inner_hsp_culling(hsps: &mut Vec<Hsp>, max_hsps: u32, inner_culling_overlap: f64) {
    inner_culling(hsps, max_hsps, inner_culling_overlap);
}

fn compare_hsps(a: &Hsp, b: &Hsp) -> Ordering {
    if a.less_than_score_position(b) {
        Ordering::Less
    } else if b.less_than_score_position(a) {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

impl Target {
    pub fn inner_culling(
        &mut self,
        max_hsps: u32,
        inner_culling_overlap: f64,
        query_contexts: usize,
    ) {
        if max_hsps == 1 {
            for i in 0..MAX_CONTEXT as usize {
                if i as i32 == self.best_context {
                    self.hsp[i].sort_by(compare_hsps);
                    self.hsp[i].truncate(1);
                } else {
                    self.hsp[i].clear();
                }
            }
            return;
        }
        let mut hsps = Vec::new();
        for frame in 0..query_contexts {
            hsps.append(&mut self.hsp[frame]);
        }
        inner_culling(&mut hsps, max_hsps, inner_culling_overlap);
        for hsp in hsps {
            self.hsp[hsp.frame as usize].push(hsp);
        }
    }

    pub fn max_hsp_culling(&mut self, max_hsps: u32, query_contexts: usize) {
        for frame in 0..query_contexts {
            max_hsp_culling(&mut self.hsp[frame], max_hsps);
        }
    }
}

impl Match {
    pub fn inner_culling(&mut self, max_hsps: u32, inner_culling_overlap: f64) {
        inner_culling(&mut self.hsps, max_hsps, inner_culling_overlap);
        if let Some(best) = self.hsps.first() {
            self.filter_evalue = best.evalue;
            self.filter_score = best.score;
        }
    }

    pub fn max_hsp_culling(&mut self, max_hsps: u32) {
        max_hsp_culling(&mut self.hsps, max_hsps);
    }
}

pub fn match_inner_culling(m: &mut Match, max_hsps: u32, inner_culling_overlap: f64) {
    m.inner_culling(max_hsps, inner_culling_overlap);
}

pub fn match_max_hsp_culling(m: &mut Match, max_hsps: u32) {
    m.max_hsp_culling(max_hsps);
}

fn sort_targets(targets: &mut [Target], toppercent: Option<f64>) {
    if toppercent.is_some() {
        targets.sort_by(|a, b| {
            b.filter_score
                .cmp(&a.filter_score)
                .then_with(|| a.block_id.cmp(&b.block_id))
        });
    } else {
        targets.sort_by(|a, b| {
            a.filter_evalue
                .total_cmp(&b.filter_evalue)
                .then_with(|| b.filter_score.cmp(&a.filter_score))
                .then_with(|| a.block_id.cmp(&b.block_id))
        });
    }
}

fn output_range<F>(
    targets: &[Target],
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> usize
where
    F: FnMut(i32) -> f64,
{
    if targets.is_empty() || targets[0].filter_evalue == f64::MAX {
        return 0;
    }
    if let Some(toppercent) = toppercent {
        let cutoff = ((1.0 - toppercent / 100.0) * bitscore(targets[0].filter_score)).max(1.0);
        targets
            .iter()
            .take_while(|target| bitscore(target.filter_score) >= cutoff)
            .count()
    } else {
        finite_output_end(targets.len(), max_target_seqs, |i| targets[i].filter_evalue)
    }
}

pub fn output_range_targets<F>(
    targets: &[Target],
    max_target_seqs: i64,
    toppercent: Option<f64>,
    bitscore: F,
) -> usize
where
    F: FnMut(i32) -> f64,
{
    output_range(targets, max_target_seqs, toppercent, bitscore)
}

fn finite_output_end<F>(len: usize, max_target_seqs: i64, mut evalue: F) -> usize
where
    F: FnMut(usize) -> f64,
{
    let mut end = usize::try_from(max_target_seqs).unwrap_or(0).min(len);
    while end > 1 && evalue(end - 1) == f64::MAX {
        end -= 1;
    }
    end
}

struct TargetCullingOverload;

impl TargetCullingOverload {
    fn culling<F>(
        targets: &mut Vec<Target>,
        sort_only: bool,
        max_target_seqs: i64,
        toppercent: Option<f64>,
        mut bitscore: F,
    ) where
        F: FnMut(i32) -> f64,
    {
        sort_targets(targets, toppercent);
        if !sort_only {
            let end = output_range(targets, max_target_seqs, toppercent, &mut bitscore);
            targets.truncate(end);
        }
    }
}

pub fn culling_targets<F>(
    targets: &mut Vec<Target>,
    sort_only: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    bitscore: F,
) where
    F: FnMut(i32) -> f64,
{
    TargetCullingOverload::culling(targets, sort_only, max_target_seqs, toppercent, bitscore);
}

fn append_hits<F>(
    targets: &mut Vec<Target>,
    mut hits: Vec<Target>,
    with_culling: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> bool
where
    F: FnMut(i32) -> f64,
{
    if hits.is_empty() {
        return false;
    }
    let mut new_hits =
        toppercent.is_none() && i64::try_from(targets.len()).is_ok_and(|len| len < max_target_seqs);
    let mut append = !with_culling || new_hits;
    culling_targets(targets, append, max_target_seqs, toppercent, &mut bitscore);
    let max_score = hits.iter().map(|hit| hit.filter_score).max().unwrap_or(0);
    let min_evalue = hits
        .iter()
        .map(|hit| hit.filter_evalue)
        .fold(f64::MAX, f64::min);
    let range_end = output_range(targets, max_target_seqs, toppercent, &mut bitscore);
    if targets.is_empty()
        || (toppercent.is_none() && min_evalue <= targets[range_end - 1].filter_evalue)
        || toppercent.is_some_and(|percent| {
            max_score
                >= ((1.0 - percent / 100.0) * f64::from(targets[range_end - 1].filter_score)) as i32
        })
    {
        append = true;
        new_hits = true;
    }
    if append {
        targets.append(&mut hits);
    }
    new_hits
}

pub fn append_hits_targets<F>(
    targets: &mut Vec<Target>,
    hits: Vec<Target>,
    with_culling: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    bitscore: F,
) -> bool
where
    F: FnMut(i32) -> f64,
{
    append_hits(
        targets,
        hits,
        with_culling,
        max_target_seqs,
        toppercent,
        bitscore,
    )
}

fn sort_matches(matches: &mut [Match], toppercent: Option<f64>) {
    if toppercent.is_some() {
        matches.sort_by(|a, b| {
            b.filter_score
                .cmp(&a.filter_score)
                .then_with(|| a.target_block_id.cmp(&b.target_block_id))
        });
    } else {
        matches.sort_by(|a, b| {
            a.filter_evalue
                .total_cmp(&b.filter_evalue)
                .then_with(|| b.filter_score.cmp(&a.filter_score))
                .then_with(|| a.target_block_id.cmp(&b.target_block_id))
        });
    }
}

pub fn output_range_matches<F>(
    matches: &[Match],
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> usize
where
    F: FnMut(i32) -> f64,
{
    if matches.is_empty() || matches[0].filter_evalue == f64::MAX {
        return 0;
    }
    if let Some(toppercent) = toppercent {
        let cutoff = ((1.0 - toppercent / 100.0) * bitscore(matches[0].filter_score)).max(1.0);
        matches
            .iter()
            .take_while(|target| bitscore(target.filter_score) >= cutoff)
            .count()
    } else {
        finite_output_end(matches.len(), max_target_seqs, |i| matches[i].filter_evalue)
    }
}

struct MatchCullingOverload;

impl MatchCullingOverload {
    fn culling<F>(
        matches: &mut Vec<Match>,
        max_target_seqs: i64,
        toppercent: Option<f64>,
        mut bitscore: F,
    ) where
        F: FnMut(i32) -> f64,
    {
        sort_matches(matches, toppercent);
        let end = output_range_matches(matches, max_target_seqs, toppercent, &mut bitscore);
        matches.truncate(end);
    }
}

pub fn culling_matches<F>(
    matches: &mut Vec<Match>,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    bitscore: F,
) where
    F: FnMut(i32) -> f64,
{
    MatchCullingOverload::culling(matches, max_target_seqs, toppercent, bitscore);
}

pub fn culling_matches_with_sort_only<F>(
    matches: &mut Vec<Match>,
    sort_only: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) where
    F: FnMut(i32) -> f64,
{
    sort_matches(matches, toppercent);
    if !sort_only {
        let end = output_range_matches(matches, max_target_seqs, toppercent, &mut bitscore);
        matches.truncate(end);
    }
}

pub fn append_hits_matches<F>(
    targets: &mut Vec<Match>,
    mut hits: Vec<Match>,
    with_culling: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> bool
where
    F: FnMut(i32) -> f64,
{
    if hits.is_empty() {
        return false;
    }
    let mut new_hits =
        toppercent.is_none() && i64::try_from(targets.len()).is_ok_and(|len| len < max_target_seqs);
    let mut append = !with_culling || new_hits;
    culling_matches_with_sort_only(targets, append, max_target_seqs, toppercent, &mut bitscore);
    let max_score = hits.iter().map(|hit| hit.filter_score).max().unwrap_or(0);
    let min_evalue = hits
        .iter()
        .map(|hit| hit.filter_evalue)
        .fold(f64::MAX, f64::min);
    let range_end = output_range_matches(targets, max_target_seqs, toppercent, &mut bitscore);
    if targets.is_empty()
        || (toppercent.is_none() && min_evalue <= targets[range_end - 1].filter_evalue)
        || toppercent.is_some_and(|percent| {
            max_score
                >= ((1.0 - percent / 100.0) * f64::from(targets[range_end - 1].filter_score)) as i32
        })
    {
        append = true;
        new_hits = true;
    }
    if append {
        targets.append(&mut hits);
    }
    new_hits
}

pub fn filter_hsp(
    hsp: &Hsp,
    source_query_len: u32,
    query_title: &str,
    subject_len: u32,
    subject_title: Option<&str>,
    query_seq: &[Letter],
    subject_seq: &[Letter],
    min_id: f64,
    approx_min_id: f64,
    query_cover: f64,
    subject_cover: f64,
    query_or_target_cover: f64,
    no_self_hits: bool,
) -> bool {
    filter_hsp_with_cluster_threshold(
        hsp,
        source_query_len,
        query_title,
        subject_len,
        subject_title,
        query_seq,
        subject_seq,
        min_id,
        approx_min_id,
        query_cover,
        subject_cover,
        query_or_target_cover,
        no_self_hits,
        true,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn filter_hsp_with_cluster_threshold(
    hsp: &Hsp,
    source_query_len: u32,
    query_title: &str,
    subject_len: u32,
    subject_title: Option<&str>,
    query_seq: &[Letter],
    subject_seq: &[Letter],
    min_id: f64,
    approx_min_id: f64,
    query_cover: f64,
    subject_cover: f64,
    query_or_target_cover: f64,
    no_self_hits: bool,
    cluster_threshold_passed: bool,
) -> bool {
    let qcov = hsp.query_cover_percent(source_query_len);
    let tcov = hsp.subject_cover_percent(subject_len);
    !cluster_threshold_passed
        || hsp.id_percent() < min_id
        || (approx_min_id > 0.0 && hsp.approx_id < approx_min_id)
        || qcov < query_cover
        || tcov < subject_cover
        || (qcov < query_or_target_cover && tcov < query_or_target_cover)
        || (no_self_hits && query_seq == subject_seq && Some(query_title) == subject_title)
}

impl Match {
    #[allow(clippy::too_many_arguments)]
    pub fn apply_filters(
        &mut self,
        source_query_len: u32,
        query_title: &str,
        query_seq: &[Letter],
        subject_len: u32,
        subject_title: Option<&str>,
        subject_seq: &[Letter],
        config: &HspFilterConfig,
    ) {
        self.hsps.retain(|hsp| {
            !filter_hsp(
                hsp,
                source_query_len,
                query_title,
                subject_len,
                subject_title,
                query_seq,
                subject_seq,
                config.min_id,
                config.approx_min_id,
                config.query_cover,
                config.subject_cover,
                config.query_or_target_cover,
                config.no_self_hits,
            )
        });
        if let Some(best) = self.hsps.first() {
            self.filter_evalue = best.evalue;
            self.filter_score = best.score;
        } else {
            self.filter_evalue = f64::MAX;
            self.filter_score = 0;
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn match_apply_filters(
    target: &mut Match,
    source_query_len: u32,
    query_title: &str,
    query_seq: &[Letter],
    subject_len: u32,
    subject_title: Option<&str>,
    subject_seq: &[Letter],
    min_id: f64,
    approx_min_id: f64,
    query_cover: f64,
    subject_cover: f64,
    query_or_target_cover: f64,
    no_self_hits: bool,
) {
    target.apply_filters(
        source_query_len,
        query_title,
        query_seq,
        subject_len,
        subject_title,
        subject_seq,
        &HspFilterConfig {
            min_id,
            approx_min_id,
            query_cover,
            subject_cover,
            query_or_target_cover,
            no_self_hits,
            cluster_threshold_present: false,
        },
    );
}

#[allow(clippy::too_many_arguments)]
pub fn apply_filters_matches(
    matches: &mut [Match],
    source_query_len: u32,
    query_title: &str,
    query_seq: &[Letter],
    subject_len: u32,
    subject_title: Option<&str>,
    subject_seq: &[Letter],
    min_id: f64,
    approx_min_id: f64,
    query_cover: f64,
    subject_cover: f64,
    query_or_target_cover: f64,
    no_self_hits: bool,
) {
    let config = HspFilterConfig {
        min_id,
        approx_min_id,
        query_cover,
        subject_cover,
        query_or_target_cover,
        no_self_hits,
        cluster_threshold_present: false,
    };
    if filters_enabled(&config) {
        for target in matches {
            target.apply_filters(
                source_query_len,
                query_title,
                query_seq,
                subject_len,
                subject_title,
                subject_seq,
                &config,
            );
        }
    }
}

fn filters_enabled(config: &HspFilterConfig) -> bool {
    config.min_id > 0.0
        || config.approx_min_id > 0.0
        || config.query_cover > 0.0
        || config.subject_cover > 0.0
        || config.query_or_target_cover > 0.0
        || config.no_self_hits
        || config.cluster_threshold_present
}

/// Exact block-aware counterpart of the C++ range `apply_filters` overload.
pub fn apply_filters(
    matches: &mut [Match],
    source_query_len: u32,
    query_title: &str,
    query_seq: &[Letter],
    targets: &Block,
    config: &HspFilterConfig,
) -> Result<(), String> {
    if !filters_enabled(config) {
        return Ok(());
    }
    for target in matches {
        let id = target.target_block_id as usize;
        let subject_seq = targets.seqs().get(id);
        let subject_title = if config.no_self_hits {
            Some(
                std::str::from_utf8(targets.ids()?.get(id))
                    .map_err(|error| format!("Invalid UTF-8 target title: {error}"))?,
            )
        } else {
            None
        };
        target.apply_filters(
            source_query_len,
            query_title,
            query_seq,
            subject_seq.len() as u32,
            subject_title,
            subject_seq,
            config,
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SequenceType;
    use crate::util::interval::Interval;

    fn hsp(score: i32, evalue: f64, begin: i32, end: i32) -> Hsp {
        let mut hsp = Hsp::new();
        hsp.score = score;
        hsp.evalue = evalue;
        hsp.length = end - begin;
        hsp.identities = hsp.length;
        hsp.query_source_range = Interval::new(begin, end);
        hsp.query_range = Interval::new(begin, end);
        hsp.subject_range = Interval::new(begin, end);
        hsp
    }

    fn target(block_id: u32, score: i32, evalue: f64) -> Target {
        let mut target = Target::new(block_id, &[0, 1, 2], score, None);
        target.add_hit(hsp(score, evalue, 0, 3));
        target
    }

    fn target_match(block_id: u32, score: i32, evalue: f64) -> Match {
        let mut target = Match::new(block_id, block_id as u64);
        target.filter_score = score;
        target.filter_evalue = evalue;
        target.hsps.push(hsp(score, evalue, 0, 3));
        target
    }

    #[test]
    fn inner_culling_sorts_removes_enveloped_and_refreshes_match_cache() {
        let mut target = target_match(0, 0, f64::MAX);
        target.hsps = vec![
            hsp(80, 1.0e-8, 200, 260),
            hsp(90, 1.0e-9, 10, 90),
            hsp(100, 1.0e-10, 0, 100),
        ];
        target.inner_culling(0, 70.0);
        assert_eq!(
            target.hsps.iter().map(|h| h.score).collect::<Vec<_>>(),
            [100, 80]
        );
        assert_eq!((target.filter_score, target.filter_evalue), (100, 1.0e-10));
        target.max_hsp_culling(1);
        assert_eq!(target.hsps.len(), 1);
    }

    #[test]
    fn target_max_one_keeps_only_best_context() {
        let mut target = Target::new(3, &[0, 1], 0, None);
        let mut first = hsp(70, 1.0e-4, 0, 2);
        first.frame = 0;
        let mut best = hsp(90, 1.0e-8, 0, 2);
        best.frame = 1;
        target.add_hit(first);
        target.add_hit(best);
        target.inner_culling(1, 50.0, 2);
        assert!(target.hsp[0].is_empty());
        assert_eq!(target.hsp[1][0].score, 90);
    }

    #[test]
    fn target_and_match_culling_preserve_ties_limits_and_toppercent() {
        let mut targets = vec![
            target(2, 80, 1.0e-5),
            target(1, 100, 1.0e-20),
            target(0, 100, 1.0e-20),
        ];
        culling_targets(&mut targets, false, 2, None, f64::from);
        assert_eq!(
            targets.iter().map(|t| t.block_id).collect::<Vec<_>>(),
            [0, 1]
        );

        let mut matches = vec![
            target_match(2, 89, 1.0e-5),
            target_match(1, 90, 1.0e-4),
            target_match(0, 100, 1.0e-3),
        ];
        culling_matches(&mut matches, 9, Some(10.0), f64::from);
        assert_eq!(
            matches
                .iter()
                .map(|m| m.target_block_id)
                .collect::<Vec<_>>(),
            [0, 1]
        );
    }

    #[test]
    fn append_hits_only_retains_batches_that_can_enter_output_range() {
        let mut targets = vec![target(0, 100, 1.0e-20), target(1, 90, 1.0e-10)];
        assert!(!append_hits_targets(
            &mut targets,
            vec![target(2, 80, 1.0e-5)],
            true,
            2,
            None,
            f64::from,
        ));
        assert!(append_hits_targets(
            &mut targets,
            vec![target(3, 110, 1.0e-30)],
            true,
            2,
            None,
            f64::from,
        ));
    }

    #[test]
    fn block_aware_filters_resolve_each_matches_own_sequence_and_title() {
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2],
                Some("query"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        block
            .push_back(
                &[0, 1, 3],
                Some("other"),
                None,
                1,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let mut matches = vec![target_match(0, 100, 1.0e-10), target_match(1, 100, 1.0e-10)];
        apply_filters(
            &mut matches,
            3,
            "query",
            &[0, 1, 2],
            &block,
            &HspFilterConfig {
                no_self_hits: true,
                ..HspFilterConfig::default()
            },
        )
        .unwrap();
        assert!(matches[0].hsps.is_empty());
        assert_eq!(matches[0].filter_evalue, f64::MAX);
        assert_eq!(matches[1].hsps.len(), 1);
    }

    #[test]
    fn filter_hsp_covers_identity_coverage_combined_and_cluster_thresholds() {
        let mut alignment = hsp(80, 1.0e-6, 0, 50);
        alignment.length = 100;
        alignment.identities = 80;
        alignment.approx_id = 75.0;
        assert!(!filter_hsp(
            &alignment,
            100,
            "q",
            100,
            Some("s"),
            &[0],
            &[1],
            70.0,
            70.0,
            40.0,
            50.0,
            40.0,
            false,
        ));
        assert!(filter_hsp(
            &alignment,
            100,
            "q",
            100,
            Some("s"),
            &[0],
            &[1],
            81.0,
            70.0,
            40.0,
            50.0,
            40.0,
            false,
        ));
        assert!(filter_hsp_with_cluster_threshold(
            &alignment,
            100,
            "q",
            100,
            Some("s"),
            &[0],
            &[1],
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            false,
            false,
        ));
    }
}
