use crate::basic::packed_loc::PackedLoc;
use crate::basic::shape::Shape;
use crate::basic::value::Loc;
use crate::data::sequence_set::SequenceSet;
use crate::dp::simd_ungapped;
use crate::search::hamming::{HitField, SeedLoc as SeedLocation};
use crate::search::kmer_ranking::PackedLocId;
use crate::search::left_most::{left_most_filter, Context};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence;

/// Matches C++ `SHORT_QUERY_LEN`.
pub const SHORT_QUERY_LEN: i32 = 85;

/// Matches C++ `ungapped_cutoff`.
pub fn ungapped_cutoff<FShort, FLong>(
    query_len: i32,
    ungapped_evalue: f64,
    short_query_max_len: i32,
    short_query_ungapped_cutoff: i32,
    query_translated: bool,
    cutoff_table_short: FShort,
    cutoff_table: FLong,
) -> i32
where
    FShort: Fn(i32) -> i32,
    FLong: Fn(i32) -> i32,
{
    if ungapped_evalue == 0.0 {
        0
    } else if query_len <= short_query_max_len {
        short_query_ungapped_cutoff
    } else if query_len <= SHORT_QUERY_LEN && query_translated {
        cutoff_table_short(query_len)
    } else {
        cutoff_table(query_len)
    }
}

/// Matches C++ `ungapped_window`.
pub fn ungapped_window(
    query_len: i32,
    query_translated: bool,
    configured_ungapped_window: i32,
) -> i32 {
    if query_len <= SHORT_QUERY_LEN && query_translated {
        query_len
    } else {
        configured_ungapped_window
    }
}

/// Matches C++ `query_data(const PackedLoc q, const SequenceSet& query_seqs)`.
pub fn query_data_packed_loc(q: PackedLoc, query_seqs: &SequenceSet) -> (u32, Loc) {
    let (block_id, pos) = query_seqs.local_position(q.as_u64() as usize);
    (block_id as u32, pos as Loc)
}

/// Matches C++ `query_data(const PackedLocId q, const SequenceSet& query_seqs)`.
pub fn query_data_packed_loc_id(q: PackedLocId, query_seqs: &SequenceSet) -> (u32, Loc) {
    (
        q.block_id,
        q.loc as Loc - query_seqs.position(q.block_id as usize, 0) as Loc,
    )
}

pub struct SearchQueryOffsetConfig<'a> {
    pub score_cutoff: i32,
    pub window: i32,
    pub shape: &'a Shape,
    pub context: &'a Context<'a>,
    pub first_shape: bool,
    pub chunked: bool,
    pub hamming_filter_id: u32,
    pub self_search: bool,
    pub skip_left_most: bool,
    pub left_most_interval: i32,
    pub global_ranking_targets: bool,
    pub keep_target_id: bool,
    pub score_matrix: &'a ScoreMatrix,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SearchQueryOffsetHit {
    pub query_id: u32,
    pub subject: u64,
    pub seed_offset: Loc,
    pub score: u16,
    pub target_block_id: Option<u32>,
    pub global_ranking: bool,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SearchQueryOffsetResult {
    pub tentative_matches2: usize,
    pub tentative_matches3: usize,
    pub hits: Vec<SearchQueryOffsetHit>,
}

/// Matches C++ `search_query_offset`.
pub fn search_query_offset<SeedLoc>(
    q: SeedLoc,
    s: &[SeedLoc],
    hits: &[u32],
    query_seqs: &SequenceSet,
    ref_seqs: &SequenceSet,
    cfg: &SearchQueryOffsetConfig<'_>,
) -> SearchQueryOffsetResult
where
    SeedLoc: SeedLocation,
{
    const LANES: usize = 16;

    let query_id = q.block_id(query_seqs) as u32;
    let seed_offset = q.loc() as Loc - query_seqs.position(query_id as usize, 0) as Loc;
    let query_pos = q.loc();
    let window = cfg.window.max(0) as usize;
    let window_start = query_pos.saturating_sub(window);
    let window_left = query_pos - window_start;
    let query_window_end = (window_start + window * 2).min(query_seqs.data().len());
    let query_window = &query_seqs.data()[window_start..query_window_end];
    let query_clipped = sequence::clip(query_window, window_left as i32);
    let window_left =
        query_pos - (query_clipped.as_ptr() as usize - query_seqs.data().as_ptr() as usize);
    let window_clipped = query_clipped.len();

    let interval_mod = if cfg.left_most_interval > 0 {
        seed_offset % cfg.left_most_interval
    } else {
        window_left as Loc
    };
    let interval_overhang = (window_left as Loc - interval_mod).max(0) as usize;

    let mut out = SearchQueryOffsetResult::default();
    let mut i = 0usize;
    while i < hits.len() {
        let n = LANES.min(hits.len() - i);
        let mut subjects = Vec::with_capacity(n);
        let mut subject_indices = Vec::with_capacity(n);
        for j in 0..n {
            let subject_index = hits[i + j] as usize;
            let subject_pos = s[subject_index].loc();
            let subject_start = subject_pos.saturating_sub(window_left);
            subjects.push(&ref_seqs.data()[subject_start..]);
            subject_indices.push(subject_index);
        }

        let scores = if cfg.score_cutoff != 0 {
            simd_ungapped::window_ungapped_multi(
                query_clipped,
                &subjects,
                window_clipped,
                cfg.score_matrix,
            )
        } else {
            vec![i32::MAX; n]
        };

        for j in 0..n {
            let score = scores[j];
            if score <= cfg.score_cutoff {
                continue;
            }
            let subject_seed = s[subject_indices[j]];
            let target_block_id = subject_seed.block_id(ref_seqs) as u32;
            if cfg.self_search && target_block_id == query_id {
                continue;
            }
            out.tentative_matches2 += 1;

            let passes_left_most = if cfg.skip_left_most {
                true
            } else {
                let subject = subjects[j];
                let q_start = interval_overhang.min(query_clipped.len());
                let s_start = interval_overhang.min(subject.len());
                left_most_filter(
                    &query_clipped[q_start..],
                    &subject[s_start..],
                    window_left as i32 - interval_overhang as i32,
                    cfg.shape.length as i32,
                    cfg.context,
                    cfg.first_shape,
                    cfg.shape,
                    cfg.score_cutoff,
                    cfg.chunked,
                    cfg.hamming_filter_id,
                )
            };

            if !passes_left_most {
                continue;
            }
            out.tentative_matches3 += 1;
            out.hits.push(SearchQueryOffsetHit {
                query_id,
                subject: subject_seed.loc() as u64,
                seed_offset,
                score: score as u16,
                target_block_id: if cfg.global_ranking_targets || cfg.keep_target_id {
                    Some(target_block_id)
                } else {
                    None
                },
                global_ranking: cfg.global_ranking_targets,
            });
        }
        i += n;
    }
    out
}

/// Matches C++ `search_tile`.
pub fn search_tile<SeedLoc, F>(
    hits: &mut HitField,
    query_begin: i32,
    subject_begin: i32,
    q: &[SeedLoc],
    s: &[SeedLoc],
    mut search_query_offset: F,
) -> usize
where
    SeedLoc: Copy,
    F: FnMut(SeedLoc, &[SeedLoc], &[u32]),
{
    let query_count = hits.query_count();
    let q_begin = query_begin as usize;
    let s_begin = subject_begin as usize;
    let mut tentative_matches1 = 0usize;
    for i in 0..query_count {
        let h = hits.hits(i).to_vec();
        if h.is_empty() {
            continue;
        }
        tentative_matches1 += h.len();
        search_query_offset(q[q_begin + i], &s[s_begin..], &h);
    }
    tentative_matches1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::reduction::Reduction;
    use crate::util::algo::PatternMatcher;

    fn test_stage2_cfg<'a>(
        shape: &'a Shape,
        context: &'a Context<'a>,
        score_matrix: &'a ScoreMatrix,
    ) -> SearchQueryOffsetConfig<'a> {
        SearchQueryOffsetConfig {
            score_cutoff: 0,
            window: 4,
            shape,
            context,
            first_shape: true,
            chunked: false,
            hamming_filter_id: 0,
            self_search: false,
            skip_left_most: true,
            left_most_interval: 0,
            global_ranking_targets: false,
            keep_target_id: false,
            score_matrix,
        }
    }

    #[test]
    fn test_ungapped_cutoff_branch_order() {
        assert_eq!(ungapped_cutoff(50, 0.0, 60, 7, true, |_| 11, |_| 13), 0);
        assert_eq!(ungapped_cutoff(50, 1.0, 60, 7, true, |_| 11, |_| 13), 7);
        assert_eq!(ungapped_cutoff(80, 1.0, 60, 7, true, |_| 11, |_| 13), 11);
        assert_eq!(ungapped_cutoff(80, 1.0, 60, 7, false, |_| 11, |_| 13), 13);
        assert_eq!(ungapped_cutoff(100, 1.0, 60, 7, true, |_| 11, |_| 13), 13);
    }

    #[test]
    fn test_ungapped_window() {
        assert_eq!(ungapped_window(85, true, 40), 85);
        assert_eq!(ungapped_window(86, true, 40), 40);
        assert_eq!(ungapped_window(85, false, 40), 40);
    }

    #[test]
    fn test_query_data() {
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 1, 2, 3]);
        seqs.push(&[4, 5, 6]);

        assert_eq!(
            query_data_packed_loc(PackedLoc::new(seqs.position(1, 2) as u64), &seqs),
            (1, 2)
        );
        assert_eq!(
            query_data_packed_loc_id(PackedLocId::new(seqs.position(1, 2) as u64, 1), &seqs),
            (1, 2)
        );
    }

    #[test]
    fn test_search_tile_dispatches_nonempty_hit_rows() {
        let mut hits = HitField::new();
        hits.init(3, 8);
        hits.set(0, 2, true);
        hits.set(0, 5, true);
        hits.set(2, 1, true);

        let q = [10u32, 11, 12, 13, 14];
        let s = [20u32, 21, 22, 23, 24, 25];
        let mut calls = Vec::new();
        let count = search_tile(&mut hits, 1, 2, &q, &s, |query, subjects, h| {
            calls.push((query, subjects[0], subjects.len(), h.to_vec()));
        });

        assert_eq!(count, 3);
        assert_eq!(calls, vec![(11, 22, 4, vec![2, 5]), (13, 22, 4, vec![1])]);
    }

    #[test]
    fn test_search_tile_skips_empty_rows() {
        let mut hits = HitField::new();
        hits.init(2, 8);
        let q = [10u32, 11];
        let s = [20u32, 21];
        let mut calls = 0usize;
        let count = search_tile(&mut hits, 0, 0, &q, &s, |_query, _subjects, _h| {
            calls += 1;
        });
        assert_eq!(count, 0);
        assert_eq!(calls, 0);
    }

    #[test]
    fn test_search_query_offset_emits_hits_and_counts() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);
        let patterns = [shape.mask];
        let context = Context {
            previous_matcher: PatternMatcher::new(&patterns),
            current_matcher: PatternMatcher::new(&patterns),
            short_query_ungapped_cutoff: 0,
            seedp_mask: 0,
            reduction: &reduction,
        };
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let mut query_seqs = SequenceSet::new();
        let mut ref_seqs = SequenceSet::new();
        query_seqs.push(&[0, 1, 2, 3, 4, 5]);
        ref_seqs.push(&[0, 1, 2, 3, 4, 5]);
        ref_seqs.push(&[5, 4, 3, 2, 1, 0]);

        let q = PackedLoc::new(query_seqs.position(0, 2) as u64);
        let s = [
            PackedLoc::new(ref_seqs.position(0, 2) as u64),
            PackedLoc::new(ref_seqs.position(1, 3) as u64),
        ];
        let hits = [0u32, 1u32];
        let mut cfg = test_stage2_cfg(&shape, &context, &score_matrix);
        cfg.keep_target_id = true;

        let out = search_query_offset(q, &s, &hits, &query_seqs, &ref_seqs, &cfg);

        assert_eq!(out.tentative_matches2, 2);
        assert_eq!(out.tentative_matches3, 2);
        assert_eq!(out.hits.len(), 2);
        assert_eq!(out.hits[0].query_id, 0);
        assert_eq!(out.hits[0].seed_offset, 2);
        assert_eq!(out.hits[0].subject, ref_seqs.position(0, 2) as u64);
        assert_eq!(out.hits[0].target_block_id, Some(0));
        assert!(!out.hits[0].global_ranking);
    }

    #[test]
    fn test_search_query_offset_self_skip_and_global_ranking() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);
        let patterns = [shape.mask];
        let context = Context {
            previous_matcher: PatternMatcher::new(&patterns),
            current_matcher: PatternMatcher::new(&patterns),
            short_query_ungapped_cutoff: 0,
            seedp_mask: 0,
            reduction: &reduction,
        };
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let mut query_seqs = SequenceSet::new();
        let mut ref_seqs = SequenceSet::new();
        query_seqs.push(&[0, 1, 2, 3, 4, 5]);
        ref_seqs.push(&[0, 1, 2, 3, 4, 5]);
        ref_seqs.push(&[5, 4, 3, 2, 1, 0]);

        let q = PackedLoc::new(query_seqs.position(0, 2) as u64);
        let s = [
            PackedLoc::new(ref_seqs.position(0, 2) as u64),
            PackedLoc::new(ref_seqs.position(1, 3) as u64),
        ];
        let hits = [0u32, 1u32];
        let mut cfg = test_stage2_cfg(&shape, &context, &score_matrix);
        cfg.self_search = true;
        cfg.global_ranking_targets = true;

        let out = search_query_offset(q, &s, &hits, &query_seqs, &ref_seqs, &cfg);

        assert_eq!(out.tentative_matches2, 1);
        assert_eq!(out.tentative_matches3, 1);
        assert_eq!(out.hits.len(), 1);
        assert_eq!(out.hits[0].target_block_id, Some(1));
        assert!(out.hits[0].global_ranking);
    }
}
