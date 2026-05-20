use crate::align::hsp::{Hsp, Match};
use crate::align::target::culling_matches;
use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::translate::Frame;
use crate::basic::value::{BlockId, Letter};
use crate::data::sequence_set::SequenceSet;
use crate::dp::scan_diags::{diag_alignment, scan_diags128, scan_diags64};
use crate::dp::score_profile::make_profile8;
use crate::dp::score_profile::LongScoreProfile;
use crate::dp::swipe::Flags;
use crate::dp::ungapped::DiagonalSegment;
use crate::search::hit::Hit;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::data_structures::FlatArray;
use crate::util::interval::Interval;

/// Matches C++ `Extension::SeedHit`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeedHit {
    pub i: i32,
    pub j: i32,
    pub score: i32,
    pub frame: u32,
}

impl SeedHit {
    pub fn new(i: i32, j: i32, score: i32, frame: u32) -> Self {
        Self { i, j, score, frame }
    }

    pub fn diag(&self) -> i32 {
        self.i - self.j
    }

    pub fn less_than(&self, x: &SeedHit) -> bool {
        let d1 = self.diag();
        let d2 = x.diag();
        d1 < d2 || (d1 == d2 && self.j < x.j)
    }

    pub fn query_range(&self) -> Interval {
        Interval::new(self.i, self.i + 1)
    }

    pub fn target_range(&self) -> Interval {
        Interval::new(self.j, self.j + 1)
    }

    pub fn diag_segment(&self) -> DiagonalSegment {
        DiagonalSegment::new(self.i, self.j, 1, self.score)
    }
}

impl PartialOrd for SeedHit {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SeedHit {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.less_than(other) {
            std::cmp::Ordering::Less
        } else if other.less_than(self) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

pub type ScanDiagsFn = fn(&LongScoreProfile<i8>, &[Letter], i32, i32, i32, &mut [i32]);

/// Matches C++ `Extension::TargetScore`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TargetScore {
    pub target: u32,
    pub score: u16,
}

impl TargetScore {
    pub fn less_than(&self, x: &TargetScore) -> bool {
        self.score > x.score || (self.score == x.score && self.target < x.target)
    }
}

impl PartialOrd for TargetScore {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TargetScore {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.less_than(other) {
            std::cmp::Ordering::Less
        } else if other.less_than(self) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

/// Matches C++ `Extension::SeedHitList`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedHitList {
    pub seed_hits: FlatArray<SeedHit>,
    pub target_block_ids: Vec<BlockId>,
    pub target_scores: Vec<TargetScore>,
}

impl Default for SeedHitList {
    fn default() -> Self {
        Self {
            seed_hits: FlatArray::new(),
            target_block_ids: Vec::new(),
            target_scores: Vec::new(),
        }
    }
}

/// Matches C++ `count_targets`.
pub fn count_targets(hits: &[Hit]) -> i64 {
    let Some(first) = hits.first() else {
        return 0;
    };
    let mut last = first.subject;
    let mut n = 1i64;
    for hit in &hits[1..] {
        if hit.subject != last {
            last = hit.subject;
            n += 1;
        }
    }
    n
}

/// Matches C++ `Extension::load_hits`.
pub fn load_hits(hits: &mut [Hit], ref_seqs: &SequenceSet, query_contexts: u32) -> SeedHitList {
    hits.sort_by(Hit::cmp_subject);
    let targets = count_targets(hits);
    let mut list = SeedHitList::default();
    list.seed_hits.reserve(targets as u64, hits.len() as u64);
    list.target_block_ids.reserve(targets as usize);
    list.target_scores.reserve(targets as usize);

    if hits.is_empty() {
        return list;
    }

    let total_subjects = ref_seqs.len();
    let mut target = u32::MAX;
    let mut score = 0u16;
    if (total_subjects as f64).log2() * (hits.len() as f64) < total_subjects as f64 / 10.0 {
        for hit in hits {
            let (t, local_pos) = ref_seqs.local_position(hit.subject as usize);
            let t = t as u32;
            if t != target {
                if target != u32::MAX {
                    list.target_scores.push(TargetScore {
                        target: (list.target_block_ids.len() - 1) as u32,
                        score,
                    });
                    score = 0;
                }
                list.seed_hits.next();
                target = t;
                list.target_block_ids.push(target);
            }
            list.seed_hits.push_item(SeedHit::new(
                hit.seed_offset as i32,
                local_pos as i32,
                hit.score as i32,
                hit.query % query_contexts,
            ));
            score = score.max(hit.score);
        }
    } else {
        let limits = ref_seqs.offsets();
        let mut limit_idx = 0usize;
        for hit in hits {
            let subject_offset = hit.subject as usize;
            while limits[limit_idx] <= subject_offset {
                limit_idx += 1;
            }
            let t = (limit_idx - 1) as u32;
            if t != target {
                if target != u32::MAX {
                    list.target_scores.push(TargetScore {
                        target: (list.target_block_ids.len() - 1) as u32,
                        score,
                    });
                    score = 0;
                }
                list.seed_hits.next();
                list.target_block_ids.push(t);
                target = t;
            }
            list.seed_hits.push_item(SeedHit::new(
                hit.seed_offset as i32,
                (subject_offset - limits[limit_idx - 1]) as i32,
                hit.score as i32,
                hit.query % query_contexts,
            ));
            score = score.max(hit.score);
        }
    }
    if target != u32::MAX {
        list.target_scores.push(TargetScore {
            target: (list.target_block_ids.len() - 1) as u32,
            score,
        });
    }
    list
}

/// Matches C++ `Extension::gapped_filter(const SeedHit&, ..., int band, int window, f)`.
pub fn gapped_filter_hit(
    hit: &SeedHit,
    query_profile: &[LongScoreProfile<i8>],
    target: &[Letter],
    band: i32,
    window: i32,
    f: ScanDiagsFn,
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
) -> i32 {
    let slen = target.len() as i32;
    let d = (hit.diag() - band / 2).max(-(slen - 1));
    let j0 = (hit.j - window).max(0);
    let j1 = (hit.j + window).min(slen);
    let mut scores = vec![0i32; band as usize];
    f(
        &query_profile[hit.frame as usize],
        target,
        d,
        j0,
        j1,
        &mut scores,
    );
    diag_alignment(&scores, gap_open, gap_extend, gapped_filter_diag_score)
}

/// Matches C++ `Extension::gapped_filter(begin, end, ...)`.
pub fn gapped_filter_target<F1, F2>(
    hits: &[SeedHit],
    query_profile: &[LongScoreProfile<i8>],
    target: &[Letter],
    stat: &mut Statistics,
    qlen: i32,
    cutoff_gapped1_new: F1,
    cutoff_gapped2_new: F2,
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
    query_translated: bool,
    gapped_filter_window: i32,
) -> bool
where
    F1: Fn(i32, i32) -> i32,
    F2: Fn(i32, i32) -> i32,
{
    const WINDOW1: i32 = 100;
    const MIN_STAGE2_QLEN: i32 = 100;

    let slen = target.len() as i32;
    for hit in hits {
        stat.inc(StatValue::GappedFilterHits1, 1);
        let f1 = gapped_filter_hit(
            hit,
            query_profile,
            target,
            64,
            WINDOW1,
            scan_diags64,
            gap_open,
            gap_extend,
            gapped_filter_diag_score,
        );
        if f1 > cutoff_gapped1_new(qlen, slen) {
            stat.inc(StatValue::GappedFilterHits2, 1);
            if qlen < MIN_STAGE2_QLEN && query_translated {
                return true;
            }
            let f2 = gapped_filter_hit(
                hit,
                query_profile,
                target,
                128,
                gapped_filter_window,
                scan_diags128,
                gap_open,
                gap_extend,
                gapped_filter_diag_score,
            );
            if f2 > cutoff_gapped2_new(qlen, slen) {
                return true;
            }
        }
    }
    false
}

/// Matches C++ `Extension::gapped_filter(const Sequence*, ..., FlatArray<SeedHit>::Iterator, ...)`.
pub fn gapped_filter_seed_hits<F1, F2>(
    query: &[&[Letter]],
    query_cbs: Option<&[&[i8]]>,
    seed_hits: &FlatArray<SeedHit>,
    target_block_ids: &[BlockId],
    target_seqs: &SequenceSet,
    stat: &mut Statistics,
    flags: Flags,
    cutoff_gapped1_new: F1,
    cutoff_gapped2_new: F2,
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
    query_translated: bool,
    gapped_filter_window: i32,
    score_matrix: &ScoreMatrix,
) -> (FlatArray<SeedHit>, Vec<BlockId>)
where
    F1: Fn(i32, i32) -> i32 + Copy,
    F2: Fn(i32, i32) -> i32 + Copy,
{
    let n = seed_hits.size() as usize;
    let mut hits_out = FlatArray::new();
    let mut target_ids_out = Vec::new();
    if n == 0 {
        return (hits_out, target_ids_out);
    }

    let mut query_profile = Vec::with_capacity(query.len());
    for i in 0..query.len() {
        let cbs = query_cbs.and_then(|all| all.get(i).copied());
        query_profile.push(make_profile8(query[i], cbs, 0, score_matrix));
    }

    let _ = flags;
    for i in 0..n {
        let target_block_id = target_block_ids[i];
        let target = target_seqs.get(target_block_id as usize);
        if gapped_filter_target(
            seed_hits.range(i as u64),
            &query_profile,
            target,
            stat,
            query_profile[0].length() as i32,
            cutoff_gapped1_new,
            cutoff_gapped2_new,
            gap_open,
            gap_extend,
            gapped_filter_diag_score,
            query_translated,
            gapped_filter_window,
        ) {
            target_ids_out.push(target_block_id);
            hits_out.push_back(seed_hits.range(i as u64));
        }
    }

    (hits_out, target_ids_out)
}

/// Matches C++ `seed_only_hsp`.
pub fn seed_only_hsp(hit: &SeedHit, query_source_len: u32) -> Hsp {
    let mut hsp = Hsp::new();
    hsp.seed_only = true;
    hsp.score = hit.score;
    hsp.evalue = 0.0;
    hsp.d_begin = hit.diag();
    hsp.d_end = hit.diag();
    hsp.set_begin(
        hit.i,
        hit.j,
        Frame::from_index(hit.frame as i32),
        query_source_len as i32,
    );
    hsp.set_end(
        hit.i + 1,
        hit.j + 1,
        Frame::from_index(hit.frame as i32),
        query_source_len as i32,
    );
    hsp.subject_source_range = hsp.subject_range;
    hsp
}

/// Matches C++ `seed_only_matches`.
pub fn seed_only_matches<F>(
    l: &SeedHitList,
    query_source_len: u32,
    target_seqs: Option<&SequenceSet>,
    target_oids: Option<&[u64]>,
    max_hsps: u32,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> Vec<Match>
where
    F: FnMut(i32) -> f64,
{
    let mut matches = Vec::with_capacity(l.target_scores.len());
    for target_score in &l.target_scores {
        let target_block_id = l.target_block_ids[target_score.target as usize];
        let target_oid = target_oids
            .and_then(|oids| oids.get(target_block_id as usize).copied())
            .unwrap_or(target_block_id as u64);
        let mut m = if let Some(target_seqs) = target_seqs {
            let mut m = Match::new_extension(
                target_block_id,
                target_seqs.get(target_block_id as usize),
                None,
                target_score.score as i32,
                target_score.score as i32,
                0.0,
            );
            m.target_oid = target_oid;
            m
        } else {
            let mut m = Match::new(target_block_id, target_oid);
            m.ungapped_score = target_score.score as i32;
            m.filter_score = target_score.score as i32;
            m.filter_evalue = 0.0;
            m
        };
        m.hsps = l
            .seed_hits
            .range(target_score.target as u64)
            .iter()
            .map(|hit| seed_only_hsp(hit, query_source_len))
            .collect();
        m.hsps.sort_by(|a, b| {
            if a.less_than_score_position(b) {
                std::cmp::Ordering::Less
            } else if b.less_than_score_position(a) {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });
        if max_hsps > 0 && m.hsps.len() > max_hsps as usize {
            m.hsps.truncate(max_hsps as usize);
        }
        matches.push(m);
    }
    culling_matches(&mut matches, max_target_seqs, toppercent, &mut bitscore);
    matches
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dp::score_profile::make_profile8;
    use crate::stats::score_matrix::ScoreMatrix;

    fn make_test_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn test_seed_hit_cpp_accessors_and_ordering() {
        let a = SeedHit::new(5, 3, 11, 2);
        let b = SeedHit::new(6, 3, 7, 1);
        assert_eq!(a.diag(), 2);
        assert!(a.less_than(&b));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Less);
        assert_eq!(a.query_range(), Interval::new(5, 6));
        assert_eq!(a.target_range(), Interval::new(3, 4));
        assert_eq!(a.diag_segment(), DiagonalSegment::new(5, 3, 1, 11));
    }

    #[test]
    fn test_gapped_filter_hit_scores_seed_diagonal() {
        let sm = make_test_matrix();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let target = query.clone();
        let profile = vec![make_profile8(&query, None, 128, &sm)];
        let hit = SeedHit::new(3, 3, 0, 0);

        let score = gapped_filter_hit(
            &hit,
            &profile,
            &target,
            64,
            100,
            scan_diags64,
            sm.gap_open(),
            sm.gap_extend(),
            1,
        );
        assert!(score > 0);
    }

    #[test]
    fn test_gapped_filter_target_two_stage_and_translated_short_query() {
        let sm = make_test_matrix();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let target = query.clone();
        let profile = vec![make_profile8(&query, None, 128, &sm)];
        let hits = [SeedHit::new(3, 3, 0, 0)];
        let mut stat = Statistics::new();

        let accepted = gapped_filter_target(
            &hits,
            &profile,
            &target,
            &mut stat,
            query.len() as i32,
            |_, _| -1,
            |_, _| i32::MAX,
            sm.gap_open(),
            sm.gap_extend(),
            1,
            true,
            100,
        );
        assert!(accepted);
        assert_eq!(stat.get(StatValue::GappedFilterHits1), 1);
        assert_eq!(stat.get(StatValue::GappedFilterHits2), 1);

        let mut stat = Statistics::new();
        let rejected = gapped_filter_target(
            &hits,
            &profile,
            &target,
            &mut stat,
            120,
            |_, _| -1,
            |_, _| i32::MAX,
            sm.gap_open(),
            sm.gap_extend(),
            1,
            false,
            100,
        );
        assert!(!rejected);
        assert_eq!(stat.get(StatValue::GappedFilterHits1), 1);
        assert_eq!(stat.get(StatValue::GappedFilterHits2), 1);
    }

    #[test]
    fn test_count_targets_and_load_hits_group_by_subject_block() {
        let mut ref_seqs = SequenceSet::new();
        ref_seqs.push(&[0, 1, 2, 3]);
        ref_seqs.push(&[4, 5, 6]);

        let mut hits = vec![
            Hit::with_score(5, ref_seqs.position(1, 2) as u64, 7, 30),
            Hit::with_score(4, ref_seqs.position(0, 1) as u64, 6, 20),
            Hit::with_score(2, ref_seqs.position(0, 3) as u64, 8, 40),
        ];
        assert_eq!(count_targets(&hits), 3);

        let list = load_hits(&mut hits, &ref_seqs, 3);
        assert_eq!(list.target_block_ids, vec![0, 1]);
        assert_eq!(list.seed_hits.size(), 2);
        assert_eq!(
            list.seed_hits.range(0),
            &[SeedHit::new(6, 1, 20, 1), SeedHit::new(8, 3, 40, 2),]
        );
        assert_eq!(list.seed_hits.range(1), &[SeedHit::new(7, 2, 30, 2)]);
        assert_eq!(
            list.target_scores,
            vec![
                TargetScore {
                    target: 0,
                    score: 40,
                },
                TargetScore {
                    target: 1,
                    score: 30,
                },
            ]
        );
        assert!(list.target_scores[0].less_than(&list.target_scores[1]));
        assert_eq!(
            list.target_scores[0].cmp(&list.target_scores[1]),
            std::cmp::Ordering::Less
        );
    }

    #[test]
    fn test_gapped_filter_seed_hits_filters_target_groups() {
        let sm = make_test_matrix();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut target_seqs = SequenceSet::new();
        target_seqs.push(&query);
        target_seqs.push(&[7, 6, 5, 4, 3, 2, 1, 0]);

        let mut seed_hits = FlatArray::new();
        seed_hits.push_back(&[SeedHit::new(3, 3, 0, 0)]);
        seed_hits.push_back(&[SeedHit::new(3, 3, 0, 0)]);
        let target_block_ids = vec![0, 1];
        let mut stat = Statistics::new();
        let query_refs: [&[Letter]; 1] = [&query];

        let (filtered_hits, filtered_targets) = gapped_filter_seed_hits(
            &query_refs,
            None,
            &seed_hits,
            &target_block_ids,
            &target_seqs,
            &mut stat,
            Flags::NONE,
            |_, slen| if slen == 8 { -1 } else { i32::MAX },
            |_, _| -1,
            sm.gap_open(),
            sm.gap_extend(),
            1,
            false,
            100,
            &sm,
        );

        assert_eq!(filtered_targets, vec![0, 1]);
        assert_eq!(filtered_hits.size(), 2);
        assert_eq!(filtered_hits.range(0), seed_hits.range(0));
        assert_eq!(filtered_hits.range(1), seed_hits.range(1));
        assert_eq!(stat.get(StatValue::GappedFilterHits1), 2);
        assert_eq!(stat.get(StatValue::GappedFilterHits2), 2);
    }

    #[test]
    fn test_seed_only_hsp_and_matches() {
        let hit = SeedHit::new(2, 5, 37, 0);
        let hsp = seed_only_hsp(&hit, 30);
        assert!(hsp.seed_only);
        assert_eq!(hsp.score, 37);
        assert_eq!(hsp.evalue, 0.0);
        assert_eq!(hsp.d_begin, hit.diag());
        assert_eq!(hsp.d_end, hit.diag());
        assert_eq!(hsp.query_range, Interval::new(2, 3));
        assert_eq!(hsp.subject_range, Interval::new(5, 6));
        assert_eq!(hsp.subject_source_range, hsp.subject_range);

        let mut list = SeedHitList::default();
        list.target_block_ids = vec![10, 11];
        list.target_scores = vec![
            TargetScore {
                target: 0,
                score: 37,
            },
            TargetScore {
                target: 1,
                score: 20,
            },
        ];
        list.seed_hits
            .push_back(&[SeedHit::new(2, 5, 37, 0), SeedHit::new(1, 3, 30, 0)]);
        list.seed_hits.push_back(&[SeedHit::new(4, 4, 20, 0)]);

        let matches = seed_only_matches(&list, 30, None, None, 1, 10, None, |s| s as f64);
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].target_block_id, 10);
        assert_eq!(matches[0].target_oid, 10);
        assert_eq!(matches[0].hsps.len(), 1);
        assert_eq!(matches[0].hsps[0].score, 37);
        assert!(matches[0].hsps[0].seed_only);
    }
}
