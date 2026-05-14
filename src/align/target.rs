use super::hsp::Match;
use crate::align::hsp::Hsp;
use crate::basic::value::Letter;

/// Target culling — remove redundant targets based on overlap.
///
/// When multiple targets align to the same query region, this module
/// decides which targets to keep based on score and coverage overlap.
pub struct TargetCulling {
    /// Overlap fraction threshold for considering targets redundant.
    pub overlap_threshold: f64,
}

impl TargetCulling {
    pub fn new(overlap_threshold: f64) -> Self {
        TargetCulling { overlap_threshold }
    }

    /// Apply culling to a list of matches, keeping only non-redundant ones.
    /// Matches should be pre-sorted by score (best first).
    pub fn cull(&self, matches: &mut Vec<Match>) {
        if matches.is_empty() {
            return;
        }

        // Keep track of which matches to retain
        let mut keep = vec![true; matches.len()];

        for i in 0..matches.len() {
            if !keep[i] {
                continue;
            }
            for j in (i + 1)..matches.len() {
                if !keep[j] {
                    continue;
                }
                // Check if match j is dominated by match i
                if self.is_dominated(&matches[i], &matches[j]) {
                    keep[j] = false;
                }
            }
        }

        let mut write_idx = 0;
        for (read_idx, &should_keep) in keep.iter().enumerate() {
            if should_keep {
                if write_idx != read_idx {
                    matches.swap(write_idx, read_idx);
                }
                write_idx += 1;
            }
        }
        matches.truncate(write_idx);
    }

    /// Check if target b is dominated by target a (a has higher score and overlapping HSPs).
    fn is_dominated(&self, a: &Match, b: &Match) -> bool {
        if a.top_score() <= b.top_score() {
            return false;
        }

        for hsp_b in &b.hsps {
            let mut is_covered = false;
            for hsp_a in &a.hsps {
                let overlap = hsp_a.subject_range.overlap_factor(&hsp_b.subject_range);
                if overlap >= self.overlap_threshold {
                    is_covered = true;
                    break;
                }
            }
            if !is_covered {
                return false;
            }
        }
        true
    }
}

/// Apply max-target-seqs filtering.
pub fn apply_max_target_seqs(matches: &mut Vec<Match>, max_targets: usize) {
    if matches.len() > max_targets {
        matches.truncate(max_targets);
    }
}

/// Apply top-percent filtering: keep matches within `top_percent` of the best score.
pub fn apply_top_percent(matches: &mut Vec<Match>, top_percent: f64) {
    if matches.is_empty() {
        return;
    }
    let top_score = matches[0].top_score() as f64;
    let cutoff = top_score * (1.0 - top_percent / 100.0);
    matches.retain(|m| m.top_score() as f64 >= cutoff);
}

/// Matches C++ `Extension::max_hsp_culling(list<Hsp>&)`.
pub fn max_hsp_culling(hsps: &mut Vec<Hsp>, max_hsps: u32) {
    if max_hsps > 0 && hsps.len() > max_hsps as usize {
        hsps.truncate(max_hsps as usize);
    }
}

/// Matches C++ `Extension::inner_culling(list<Hsp>&)`.
pub fn inner_hsp_culling(hsps: &mut Vec<Hsp>, max_hsps: u32, inner_culling_overlap: f64) {
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| {
        if a.less_than_score_position(b) {
            std::cmp::Ordering::Less
        } else if b.less_than_score_position(a) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });
    if max_hsps == 1 {
        hsps.truncate(1);
        return;
    }
    let overlap = inner_culling_overlap / 100.0;
    let mut i = 0;
    while i < hsps.len() {
        let enveloped = (0..i).any(|j| hsps[i].is_enveloped_by(&hsps[j], overlap));
        if enveloped {
            hsps.remove(i);
        } else {
            i += 1;
        }
    }
    if max_hsps > 0 {
        max_hsp_culling(hsps, max_hsps);
    }
}

/// Matches C++ `Match::inner_culling()`.
pub fn match_inner_culling(m: &mut Match, max_hsps: u32, inner_culling_overlap: f64) {
    inner_hsp_culling(&mut m.hsps, max_hsps, inner_culling_overlap);
}

/// Matches C++ `Match::max_hsp_culling()`.
pub fn match_max_hsp_culling(m: &mut Match, max_hsps: u32) {
    max_hsp_culling(&mut m.hsps, max_hsps);
}

/// Matches C++ `output_range()` for `Match` ranges.
pub fn output_range_matches<F>(
    matches: &[Match],
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> usize
where
    F: FnMut(i32) -> f64,
{
    if matches.is_empty() {
        return 0;
    }
    if matches[0].top_evalue() == f64::MAX {
        return 0;
    }
    if let Some(toppercent) = toppercent {
        let cutoff = ((1.0 - toppercent / 100.0) * bitscore(matches[0].top_score())).max(1.0);
        let mut i = 0;
        while i < matches.len() && bitscore(matches[i].top_score()) >= cutoff {
            i += 1;
        }
        i
    } else {
        let mut i = (max_target_seqs as usize).min(matches.len());
        while i > 1 && matches[i - 1].top_evalue() == f64::MAX {
            i -= 1;
        }
        i
    }
}

/// Matches C++ `culling(std::vector<Match>&, const Search::Config&)`.
pub fn culling_matches<F>(
    matches: &mut Vec<Match>,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) where
    F: FnMut(i32) -> f64,
{
    if toppercent.is_some() {
        matches.sort_by(|a, b| {
            b.top_score()
                .cmp(&a.top_score())
                .then_with(|| a.target_block_id.cmp(&b.target_block_id))
        });
    } else {
        matches.sort_by(|a, b| {
            a.top_evalue()
                .total_cmp(&b.top_evalue())
                .then_with(|| b.top_score().cmp(&a.top_score()))
        });
    }
    let end = output_range_matches(matches, max_target_seqs, toppercent, &mut bitscore);
    matches.truncate(end);
}

/// Matches C++ `culling(std::vector<Match>&, ..., sort_only)`.
pub fn culling_matches_with_sort_only<F>(
    matches: &mut Vec<Match>,
    sort_only: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) where
    F: FnMut(i32) -> f64,
{
    if toppercent.is_some() {
        matches.sort_by(|a, b| {
            b.top_score()
                .cmp(&a.top_score())
                .then_with(|| a.target_block_id.cmp(&b.target_block_id))
        });
    } else {
        matches.sort_by(|a, b| {
            a.top_evalue()
                .total_cmp(&b.top_evalue())
                .then_with(|| b.top_score().cmp(&a.top_score()))
        });
    }
    if !sort_only {
        let end = output_range_matches(matches, max_target_seqs, toppercent, &mut bitscore);
        matches.truncate(end);
    }
}

/// Matches C++ `append_hits(vector<Match>&, begin, end, bool, const Search::Config&)`.
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
    let mut new_hits = toppercent.is_none() && targets.len() < max_target_seqs as usize;
    let mut append = !with_culling || new_hits;

    culling_matches_with_sort_only(targets, append, max_target_seqs, toppercent, &mut bitscore);

    let mut max_score = 0;
    let mut min_evalue = f64::MAX;
    for hit in &hits {
        max_score = max_score.max(hit.top_score());
        min_evalue = min_evalue.min(hit.top_evalue());
    }

    let range_end = output_range_matches(targets, max_target_seqs, toppercent, &mut bitscore);
    if targets.is_empty()
        || range_end == 0
        || (toppercent.is_none() && min_evalue <= targets[range_end - 1].top_evalue())
        || (toppercent.is_some()
            && max_score
                >= ((1.0 - toppercent.unwrap() / 100.0) * targets[range_end - 1].top_score() as f64)
                    as i32)
    {
        append = true;
        new_hits = true;
    }

    if append {
        targets.append(&mut hits);
    }
    new_hits
}

/// Matches C++ `filter_hsp(...)` without optional MCL cluster-threshold evaluation.
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
    let qcov = hsp.query_cover_percent(source_query_len);
    let tcov = hsp.subject_cover_percent(subject_len);
    hsp.id_percent() < min_id
        || (approx_min_id > 0.0 && hsp.approx_id < approx_min_id)
        || qcov < query_cover
        || tcov < subject_cover
        || (qcov < query_or_target_cover && tcov < query_or_target_cover)
        || (no_self_hits && query_seq == subject_seq && Some(query_title) == subject_title)
}

/// Matches C++ `Match::apply_filters(...)` for an already materialized target sequence.
pub fn match_apply_filters(
    m: &mut Match,
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
    m.hsps.retain(|hsp| {
        !filter_hsp(
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
        )
    });
}

/// Matches C++ `apply_filters(vector<Match>::iterator, vector<Match>::iterator, ...)`.
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
    if min_id > 0.0
        || approx_min_id > 0.0
        || query_cover > 0.0
        || subject_cover > 0.0
        || query_or_target_cover > 0.0
        || no_self_hits
    {
        for m in matches {
            match_apply_filters(
                m,
                source_query_len,
                query_title,
                query_seq,
                subject_len,
                subject_title,
                subject_seq,
                min_id,
                approx_min_id,
                query_cover,
                subject_cover,
                query_or_target_cover,
                no_self_hits,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::basic::value::Score;
    use crate::util::interval::Interval;

    fn make_match(block_id: u32, score: Score, subject_range: (i32, i32)) -> Match {
        let mut m = Match::new(block_id, block_id as u64);
        let mut hsp = Hsp::new();
        hsp.score = score;
        hsp.evalue = 1.0 / score as f64;
        hsp.length = 100;
        hsp.identities = score.min(100);
        hsp.query_source_range = Interval::new(subject_range.0, subject_range.1);
        hsp.subject_range = Interval::new(subject_range.0, subject_range.1);
        m.hsps.push(hsp);
        m
    }

    #[test]
    fn test_max_target_seqs() {
        let mut matches = vec![
            make_match(0, 100, (0, 50)),
            make_match(1, 90, (0, 50)),
            make_match(2, 80, (0, 50)),
        ];
        apply_max_target_seqs(&mut matches, 2);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn test_top_percent() {
        let mut matches = vec![
            make_match(0, 100, (0, 50)),
            make_match(1, 95, (0, 50)),
            make_match(2, 50, (0, 50)),
        ];
        apply_top_percent(&mut matches, 10.0);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn test_inner_hsp_culling() {
        let mut a = Hsp::new();
        a.score = 100;
        a.query_source_range = Interval::new(0, 100);
        a.subject_range = Interval::new(0, 100);
        let mut b = Hsp::new();
        b.score = 90;
        b.query_source_range = Interval::new(10, 90);
        b.subject_range = Interval::new(10, 90);
        let mut c = Hsp::new();
        c.score = 80;
        c.query_source_range = Interval::new(200, 260);
        c.subject_range = Interval::new(200, 260);
        let mut hsps = vec![c, b, a];
        inner_hsp_culling(&mut hsps, 0, 70.0);
        assert_eq!(hsps.len(), 2);
        assert_eq!(hsps[0].score, 100);
        assert_eq!(hsps[1].score, 80);
    }

    #[test]
    fn test_culling_matches_evalue_and_toppercent() {
        let mut a = make_match(0, 100, (0, 50));
        a.hsps[0].evalue = 1.0e-20;
        let mut b = make_match(1, 95, (0, 50));
        b.hsps[0].evalue = 1.0e-10;
        let mut c = make_match(2, 50, (0, 50));
        c.hsps[0].evalue = 1.0e-5;
        let mut matches = vec![c.clone(), b.clone(), a.clone()];
        culling_matches(&mut matches, 2, None, |score| score as f64);
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].target_block_id, 0);
        assert_eq!(matches[1].target_block_id, 1);

        let mut matches = vec![c, b, a];
        culling_matches(&mut matches, 10, Some(10.0), |score| score as f64);
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].target_block_id, 0);
        assert_eq!(matches[1].target_block_id, 1);
    }

    #[test]
    fn test_append_hits_matches() {
        let mut a = make_match(0, 100, (0, 50));
        a.hsps[0].evalue = 1.0e-20;
        let mut b = make_match(1, 90, (0, 50));
        b.hsps[0].evalue = 1.0e-10;
        let mut targets = vec![a, b];
        let mut c = make_match(2, 80, (0, 50));
        c.hsps[0].evalue = 1.0e-5;
        assert!(!append_hits_matches(
            &mut targets,
            vec![c],
            true,
            2,
            None,
            |score| score as f64,
        ));
        assert_eq!(targets.len(), 2);

        let mut d = make_match(3, 110, (0, 50));
        d.hsps[0].evalue = 1.0e-30;
        assert!(append_hits_matches(
            &mut targets,
            vec![d],
            true,
            2,
            None,
            |score| score as f64,
        ));
        assert_eq!(targets.len(), 3);
    }

    #[test]
    fn test_filter_hsp_thresholds_and_self_hit() {
        let mut hsp = Hsp::new();
        hsp.score = 80;
        hsp.length = 100;
        hsp.identities = 80;
        hsp.approx_id = 75.0;
        hsp.query_source_range = Interval::new(0, 50);
        hsp.subject_range = Interval::new(0, 60);

        assert!(!filter_hsp(
            &hsp,
            100,
            "q",
            100,
            Some("s"),
            &[1, 2, 3],
            &[1, 2, 4],
            70.0,
            70.0,
            40.0,
            50.0,
            40.0,
            false,
        ));
        assert!(filter_hsp(
            &hsp,
            100,
            "q",
            100,
            Some("s"),
            &[1, 2, 3],
            &[1, 2, 4],
            90.0,
            70.0,
            40.0,
            50.0,
            40.0,
            false,
        ));
        assert!(filter_hsp(
            &hsp,
            100,
            "q",
            100,
            Some("q"),
            &[1, 2, 3],
            &[1, 2, 3],
            70.0,
            70.0,
            40.0,
            50.0,
            40.0,
            true,
        ));
    }

    #[test]
    fn test_apply_filters_matches() {
        let mut keep = make_match(0, 90, (0, 100));
        keep.hsps[0].length = 100;
        keep.hsps[0].identities = 90;
        let mut drop = make_match(1, 50, (0, 100));
        drop.hsps[0].length = 100;
        drop.hsps[0].identities = 50;
        let mut matches = vec![keep, drop];
        apply_filters_matches(
            &mut matches,
            100,
            "q",
            &[0, 1, 2],
            100,
            Some("s"),
            &[0, 1, 3],
            80.0,
            0.0,
            0.0,
            0.0,
            0.0,
            false,
        );
        assert_eq!(matches[0].hsps.len(), 1);
        assert!(matches[1].hsps.is_empty());
    }
}
