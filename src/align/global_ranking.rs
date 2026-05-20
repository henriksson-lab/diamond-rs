use crate::align::gapped_filter::{load_hits, SeedHit, TargetScore};
use crate::align::hsp::Match;
use crate::align::target::{self, GappedScoreConfig};
use crate::align::ungapped::UngappedStageConfig;
use crate::basic::statistics::Statistics;
use crate::basic::value::{BlockId, Letter};
use crate::data::block::Block;
use crate::dp::swipe::{Flags, HspValues};
use crate::dp::ungapped::{ungapped_window, xdrop_ungapped};
use crate::output::intermediate::IntermediateRecord;
use crate::search::hit::Hit as SearchHit;
use crate::search::sensitivity::ExtensionMode;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::data_structures::BitVector;
use crate::util::text_buffer::TextBuffer;
use std::collections::HashMap;
use std::io::{self, Read};

pub type TargetMap = HashMap<u64, BlockId>;

/// Matches C++ `GlobalRanking::QueryList::Target`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QueryListTarget {
    pub database_id: u32,
    pub score: u16,
}

/// Matches C++ `GlobalRanking::QueryList`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct QueryList {
    pub query_block_id: u32,
    pub last_query_block_id: u32,
    pub targets: Vec<QueryListTarget>,
}

/// Matches C++ `GlobalRanking::Hit`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Hit {
    pub oid: u32,
    pub score: u16,
    pub context: u8,
}

impl Hit {
    pub fn new(oid: u32, score: u16, context: u32) -> Self {
        Self {
            oid,
            score,
            context: context as u8,
        }
    }

    pub fn from_target_id(target_id: isize) -> Self {
        Self {
            oid: target_id as u32,
            score: 0,
            context: 0,
        }
    }

    pub fn less_than(&self, x: &Hit) -> bool {
        self.score > x.score || (self.score == x.score && self.oid < x.oid)
    }

    pub fn target(&self) -> u32 {
        self.oid
    }

    pub fn cmp_oid_score(x: &Hit, y: &Hit) -> std::cmp::Ordering {
        x.oid.cmp(&y.oid).then_with(|| y.score.cmp(&x.score))
    }

    pub fn cmp_oid(x: &Hit, y: &Hit) -> bool {
        x.oid == y.oid
    }
}

impl PartialOrd for Hit {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Hit {
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

/// Matches C++ `GlobalRanking::recompute_overflow_scores`.
pub fn recompute_overflow_scores(
    hits: &[SeedHit],
    query_seq: &[Letter],
    target_seq: &[Letter],
    ungapped_window_size: usize,
    score_matrix: &ScoreMatrix,
) -> u16 {
    let mut score = 0;
    for hit in hits {
        if hit.score != u8::MAX as i32 {
            continue;
        }
        let center = hit.i.max(0) as usize;
        let query_begin = center.saturating_sub(ungapped_window_size);
        let query_end = (query_begin + ungapped_window_size * 2).min(query_seq.len());
        let window_left = center.saturating_sub(query_begin);
        let subject_center = hit.j.max(0) as usize;
        let subject_begin = subject_center.saturating_sub(window_left);
        if subject_begin >= target_seq.len() {
            continue;
        }
        let s = ungapped_window(
            &query_seq[query_begin..query_end],
            &target_seq[subject_begin..],
            query_end - query_begin,
            score_matrix,
        );
        score = score.max(s);
    }
    score.min(u16::MAX as i32) as u16
}

/// Matches C++ `GlobalRanking::ranking_list`.
pub fn ranking_list(
    _query_id: BlockId,
    target_scores: &mut [TargetScore],
    target_block_ids: &[BlockId],
    seed_hits: &crate::util::data_structures::FlatArray<SeedHit>,
    query_seq: Option<&[Letter]>,
    target_block: Option<&Block>,
    global_ranking_targets: i64,
    ungapped_window_size: usize,
    score_matrix: &ScoreMatrix,
) -> Vec<Match> {
    let mut overflows = 0usize;
    for target_score in target_scores.iter_mut() {
        if target_score.score < u8::MAX as u16 {
            break;
        }
        if target_score.score == u8::MAX as u16 {
            if let (Some(query_seq), Some(target_block)) = (query_seq, target_block) {
                let target_id = target_block_ids[target_score.target as usize];
                target_score.score = recompute_overflow_scores(
                    seed_hits.range(target_score.target as u64),
                    query_seq,
                    target_block.seqs().get(target_id as usize),
                    ungapped_window_size,
                    score_matrix,
                );
            }
            overflows += 1;
        }
    }
    if overflows > 0 {
        target_scores.sort();
    }

    let n = (global_ranking_targets.max(0) as usize).min(target_scores.len());
    let mut out = Vec::with_capacity(n);
    for target_score in target_scores.iter().take(n) {
        out.push(Match::new_extension(
            target_block_ids[target_score.target as usize],
            &[],
            None,
            target_score.score as i32,
            0,
            f64::MAX,
        ));
    }
    out
}

/// Matches C++ `GlobalRanking::write_merged_query_list_intro`.
pub fn write_merged_query_list_intro(query_id: u32, buf: &mut TextBuffer) -> usize {
    let seek_pos = buf.size();
    buf.write(query_id).write(0u32);
    seek_pos
}

/// Matches C++ `GlobalRanking::write_merged_query_list`.
pub fn write_merged_query_list(
    _record: &IntermediateRecord,
    _out: &mut TextBuffer,
    _ranking_db_filter: &mut BitVector,
    _stat: &mut Statistics,
) {
}

/// Matches C++ `GlobalRanking::finish_merged_query_list`.
pub fn finish_merged_query_list(buf: &mut TextBuffer, seek_pos: usize) {
    let n = (buf.size() - seek_pos - std::mem::size_of::<u32>() * 2) as u32;
    let begin = seek_pos + std::mem::size_of::<u32>();
    buf.data_mut()[begin..begin + std::mem::size_of::<u32>()].copy_from_slice(&n.to_ne_bytes());
}

/// Matches C++ `GlobalRanking::fetch_query_targets`.
pub fn fetch_query_targets<R: Read>(
    query_list: &mut R,
    next_query: &mut u32,
) -> io::Result<QueryList> {
    let mut out = QueryList {
        last_query_block_id: *next_query,
        ..QueryList::default()
    };
    let mut b4 = [0u8; 4];
    match query_list.read_exact(&mut b4) {
        Ok(()) => out.query_block_id = u32::from_ne_bytes(b4),
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(out),
        Err(e) => return Err(e),
    }
    *next_query = out.query_block_id + 1;
    query_list.read_exact(&mut b4)?;
    let size = u32::from_ne_bytes(b4);
    let n = size as usize / 6;
    out.targets.reserve(n);
    for _ in 0..n {
        let mut target = [0u8; 4];
        let mut score = [0u8; 2];
        query_list.read_exact(&mut target)?;
        query_list.read_exact(&mut score)?;
        out.targets.push(QueryListTarget {
            database_id: u32::from_ne_bytes(target),
            score: u16::from_ne_bytes(score),
        });
    }
    Ok(out)
}

/// Matches C++ `GlobalRanking::db_filter`.
pub fn db_filter(table: &[Hit], db_size: usize) -> BitVector {
    let mut out = BitVector::with_size(db_size);
    for hit in table {
        if hit.score != 0 {
            out.set(hit.oid as usize);
        }
    }
    out
}

/// Matches C++ `GlobalRanking::target_score`.
pub fn target_score(
    hits: &mut [SeedHit],
    query_seq: &[Vec<Letter>],
    target_seq: &[Letter],
    no_reextend: bool,
    xdrop: i32,
    score_matrix: &ScoreMatrix,
) -> (i32, u32) {
    if hits.is_empty() {
        return (0, 0);
    }
    if no_reextend {
        let mut score = hits[0].score;
        let mut context = hits[0].frame;
        for hit in &hits[1..] {
            if hit.score > score {
                score = hit.score;
                context = hit.frame;
            }
        }
        return (score, context);
    }

    hits.sort();
    let first = hits[0];
    let mut d = xdrop_ungapped(
        &query_seq[first.frame as usize],
        target_seq,
        first.i.max(0) as usize,
        first.j.max(0) as usize,
        xdrop,
        score_matrix,
    );
    let mut score = d.score;
    let mut context = first.frame;
    for hit in &hits[1..] {
        if d.diag() == hit.diag() && d.subject_end() >= hit.j {
            continue;
        }
        d = xdrop_ungapped(
            &query_seq[hit.frame as usize],
            target_seq,
            hit.i.max(0) as usize,
            hit.j.max(0) as usize,
            xdrop,
            score_matrix,
        );
        if d.score > score {
            score = d.score;
            context = hit.frame;
        }
    }
    (score, context)
}

/// Matches C++ `GlobalRanking::merge_hits`.
pub fn merge_hits(
    query: usize,
    hits: &mut Vec<Hit>,
    merged: &mut Vec<Hit>,
    ranking_table: &mut [Hit],
    global_ranking_targets: usize,
    merged_count: &mut usize,
) {
    let table_begin = query * global_ranking_targets;
    let table_end = table_begin + global_ranking_targets;
    let mut used_end = table_end;
    while used_end > table_begin && ranking_table[used_end - 1].score == 0 {
        used_end -= 1;
    }
    let count = used_end - table_begin;
    hits.extend_from_slice(&ranking_table[table_begin..used_end]);
    hits.sort_by(Hit::cmp_oid_score);
    merged.clear();
    for hit in hits.iter() {
        if merged.last().is_none_or(|last| last.oid != hit.oid) {
            merged.push(*hit);
        }
    }
    merged.sort();
    let n = global_ranking_targets.min(merged.len());
    ranking_table[table_begin..table_end].fill(Hit::default());
    ranking_table[table_begin..table_begin + n].copy_from_slice(&merged[..n]);
    *merged_count = merged_count.wrapping_add(n.wrapping_sub(count));
}

/// Matches C++ `GlobalRanking::get_query_hits`.
pub fn get_query_hits(
    seed_hits: &[SearchHit],
    target_block: &Block,
    keep_target_id: bool,
) -> Vec<Hit> {
    let mut out = Vec::new();
    if seed_hits.is_empty() {
        return out;
    }
    let target_seqs = target_block.seqs();
    let mut i = 0usize;
    while i < seed_hits.len() {
        let target = if keep_target_id {
            seed_hits[i].subject as BlockId
        } else {
            target_seqs.local_position(seed_hits[i].subject as usize).0 as BlockId
        };
        let mut score = 0u16;
        i += 1;
        while i <= seed_hits.len() {
            let same = if i < seed_hits.len() {
                let next_target = if keep_target_id {
                    seed_hits[i].subject as BlockId
                } else {
                    target_seqs.local_position(seed_hits[i].subject as usize).0 as BlockId
                };
                next_target == target
            } else {
                false
            };
            score = score.max(seed_hits[i - 1].score);
            if !same {
                break;
            }
            i += 1;
        }
        out.push(Hit::new(target_block.block_id2oid(target) as u32, score, 0));
    }
    out
}

/// Matches C++ `GlobalRanking::get_query_hits_reextend`.
pub fn get_query_hits_reextend(
    _source_query_block_id: BlockId,
    seed_hits: &mut [SearchHit],
    query_seq: &[Vec<Letter>],
    target_block: &Block,
    query_contexts: u32,
    no_reextend: bool,
    xdrop: i32,
    score_matrix: &ScoreMatrix,
) -> Vec<Hit> {
    let mut out = Vec::new();
    let mut list = load_hits(seed_hits, target_block.seqs(), query_contexts);
    for i in 0..list.target_block_ids.len() {
        let target_id = list.target_block_ids[i];
        let (score, context) = target_score(
            list.seed_hits.range_mut(i as u64),
            query_seq,
            target_block.seqs().get(target_id as usize),
            no_reextend,
            xdrop,
            score_matrix,
        );
        out.push(Hit::new(
            target_block.block_id2oid(target_id) as u32,
            score.min(u16::MAX as i32) as u16,
            context,
        ));
    }
    out
}

/// Matches C++ `GlobalRanking::extend_query(const QueryList&, ...)` up to output materialization.
pub fn extend_query_from_query_list<F1, F2>(
    query_list: &QueryList,
    db2block_id: &TargetMap,
    query_seq: &[Vec<Letter>],
    query_title: &str,
    source_query_len: i32,
    query_cbs: &[Vec<i8>],
    query_comp: &[f64; crate::basic::value::TRUE_AA as usize],
    target_block: &mut Block,
    stats: &mut Statistics,
    cfg: &GappedScoreConfig,
    ungapped_cfg: &UngappedStageConfig,
    score_matrix: &ScoreMatrix,
    output_hsp_values: HspValues,
    cutoff_gapped1_new: F1,
    cutoff_gapped2_new: F2,
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
    gapped_filter_window: i32,
) -> Vec<Match>
where
    F1: Fn(i32, i32) -> i32 + Copy,
    F2: Fn(i32, i32) -> i32 + Copy,
{
    let n = query_list.targets.len();
    let mut l = crate::align::gapped_filter::SeedHitList::default();
    l.target_block_ids.reserve(n);
    l.target_scores.reserve(n);
    l.seed_hits.reserve(n as u64, 0);
    for (i, target) in query_list.targets.iter().enumerate() {
        l.target_block_ids
            .push(*db2block_id.get(&(target.database_id as u64)).unwrap());
        l.target_scores.push(TargetScore {
            target: i as u32,
            score: target.score,
        });
        l.seed_hits.next();
        l.seed_hits
            .push_item(SeedHit::new(0, 0, target.score as i32, 0));
    }
    target::extend_seed_hit_list(
        query_list.query_block_id,
        query_seq,
        query_title,
        source_query_len,
        query_cbs,
        query_comp,
        target_block,
        stats,
        Flags::FULL_MATRIX,
        &mut l,
        ExtensionMode::Full,
        cfg,
        ungapped_cfg,
        score_matrix,
        output_hsp_values,
        cutoff_gapped1_new,
        cutoff_gapped2_new,
        gap_open,
        gap_extend,
        gapped_filter_diag_score,
        gapped_filter_window,
    )
}

/// Matches C++ `GlobalRanking::extend_query(BlockId, ...)` up to output materialization.
pub fn extend_query_from_ranking_table<F1, F2>(
    source_query_block_id: BlockId,
    ranking_table: &[Hit],
    db2block_id: &TargetMap,
    query_seq: &[Vec<Letter>],
    query_title: &str,
    source_query_len: i32,
    query_cbs: &[Vec<i8>],
    query_comp: &[f64; crate::basic::value::TRUE_AA as usize],
    target_block: &mut Block,
    stats: &mut Statistics,
    cfg: &GappedScoreConfig,
    ungapped_cfg: &UngappedStageConfig,
    score_matrix: &ScoreMatrix,
    output_hsp_values: HspValues,
    cutoff_gapped1_new: F1,
    cutoff_gapped2_new: F2,
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
    gapped_filter_window: i32,
) -> Vec<Match>
where
    F1: Fn(i32, i32) -> i32 + Copy,
    F2: Fn(i32, i32) -> i32 + Copy,
{
    let n_targets = cfg.global_ranking_targets.max(0) as usize;
    let table_begin = source_query_block_id as usize * n_targets;
    let mut table_end = table_begin + n_targets;
    while table_end > table_begin && ranking_table[table_end - 1].score == 0 {
        table_end -= 1;
    }
    if table_end == table_begin {
        return Vec::new();
    }

    let mut query_list = QueryList {
        query_block_id: source_query_block_id,
        last_query_block_id: source_query_block_id,
        targets: Vec::with_capacity(table_end - table_begin),
    };
    for hit in &ranking_table[table_begin..table_end] {
        query_list.targets.push(QueryListTarget {
            database_id: hit.oid,
            score: hit.score,
        });
    }
    extend_query_from_query_list(
        &query_list,
        db2block_id,
        query_seq,
        query_title,
        source_query_len,
        query_cbs,
        query_comp,
        target_block,
        stats,
        cfg,
        ungapped_cfg,
        score_matrix,
        output_hsp_values,
        cutoff_gapped1_new,
        cutoff_gapped2_new,
        gap_open,
        gap_extend,
        gapped_filter_diag_score,
        gapped_filter_window,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::gapped_filter::SeedHitList;

    fn sm() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn test_global_ranking_hit_order_and_grouping_predicates() {
        let mut hits = vec![Hit::new(7, 20, 1), Hit::new(3, 20, 0), Hit::new(4, 25, 2)];
        hits.sort();
        assert_eq!(
            hits.iter().map(|h| h.oid).collect::<Vec<_>>(),
            vec![4, 3, 7]
        );
        assert_eq!(Hit::from_target_id(9).target(), 9);
        assert_eq!(
            Hit::cmp_oid_score(&Hit::new(2, 8, 0), &Hit::new(2, 9, 0)),
            std::cmp::Ordering::Greater
        );
        assert!(Hit::cmp_oid(&Hit::new(2, 8, 0), &Hit::new(2, 1, 1)));
    }

    #[test]
    fn test_ranking_list_limits_and_recomputes_overflow_scores() {
        let score_matrix = sm();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut query_block = Block::new();
        query_block
            .push_back(
                &query,
                None,
                None,
                0,
                crate::basic::value::SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let mut target_block = Block::new();
        target_block
            .push_back(
                &query,
                None,
                None,
                0,
                crate::basic::value::SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();

        let mut list = SeedHitList::default();
        list.seed_hits
            .push_back(&[SeedHit::new(3, 3, u8::MAX as i32, 0)]);
        list.seed_hits.push_back(&[SeedHit::new(0, 0, 44, 0)]);
        list.target_block_ids = vec![0, 0];
        list.target_scores = vec![
            TargetScore {
                target: 0,
                score: u8::MAX as u16,
            },
            TargetScore {
                target: 1,
                score: 44,
            },
        ];

        let out = ranking_list(
            0,
            &mut list.target_scores,
            &list.target_block_ids,
            &list.seed_hits,
            Some(query_block.seqs().get(0)),
            Some(&target_block),
            1,
            4,
            &score_matrix,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].target_block_id, 0);
        assert!(out[0].ungapped_score >= 44);
    }

    #[test]
    fn test_merged_query_list_intro_finish_and_fetch() {
        let mut buf = TextBuffer::new();
        let seek = write_merged_query_list_intro(17, &mut buf);
        buf.write(9u32).write(40u16);
        finish_merged_query_list(&mut buf, seek);
        assert_eq!(u32::from_ne_bytes(buf.data()[0..4].try_into().unwrap()), 17);
        assert_eq!(u32::from_ne_bytes(buf.data()[4..8].try_into().unwrap()), 6);

        let mut next = 0u32;
        let mut input = buf.data();
        let q = fetch_query_targets(&mut input, &mut next).unwrap();
        assert_eq!(q.query_block_id, 17);
        assert_eq!(q.last_query_block_id, 0);
        assert_eq!(next, 18);
        assert_eq!(
            q.targets,
            vec![QueryListTarget {
                database_id: 9,
                score: 40,
            }]
        );

        let mut empty = &[][..];
        let q = fetch_query_targets(&mut empty, &mut next).unwrap();
        assert_eq!(q.last_query_block_id, 18);
        assert!(q.targets.is_empty());
    }

    #[test]
    fn test_db_filter_sets_ranked_targets() {
        let filter = db_filter(
            &[Hit::new(1, 5, 0), Hit::new(2, 0, 0), Hit::new(3, 6, 0)],
            5,
        );
        assert!(filter.get(1));
        assert!(filter.get(3));
        assert!(!filter.get(2));
    }

    #[test]
    fn test_target_score_no_reextend_and_xdrop_paths() {
        let score_matrix = sm();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let target = query.clone();
        let query_seq = vec![query.clone()];
        let mut hits = vec![SeedHit::new(1, 1, 10, 0), SeedHit::new(4, 4, 20, 0)];
        assert_eq!(
            target_score(&mut hits, &query_seq, &target, true, 20, &score_matrix),
            (20, 0)
        );

        let mut hits = vec![SeedHit::new(4, 4, 0, 0), SeedHit::new(1, 1, 0, 0)];
        let (score, context) =
            target_score(&mut hits, &query_seq, &target, false, 20, &score_matrix);
        assert_eq!(context, 0);
        assert!(score > 20);
        assert!(hits[0].diag() <= hits[1].diag());
    }

    #[test]
    fn test_merge_hits_updates_capped_ranking_table() {
        let mut ranking_table = vec![Hit::default(); 4];
        ranking_table[0] = Hit::new(2, 30, 0);
        ranking_table[1] = Hit::new(5, 10, 0);
        let mut hits = vec![Hit::new(5, 50, 1), Hit::new(3, 40, 0)];
        let mut merged = Vec::new();
        let mut merged_count = 0usize;
        merge_hits(
            0,
            &mut hits,
            &mut merged,
            &mut ranking_table,
            4,
            &mut merged_count,
        );
        assert_eq!(
            ranking_table[..3].iter().map(|h| h.oid).collect::<Vec<_>>(),
            vec![5, 3, 2]
        );
        assert_eq!(
            ranking_table[..3]
                .iter()
                .map(|h| h.score)
                .collect::<Vec<_>>(),
            vec![50, 40, 30]
        );
        assert_eq!(ranking_table[3], Hit::default());
        assert_eq!(merged_count, 1);
    }

    #[test]
    fn test_get_query_hits_groups_targets_and_reextends() {
        let score_matrix = sm();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut target_block = Block::new();
        target_block
            .push_back(
                &query,
                None,
                None,
                10,
                crate::basic::value::SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        target_block
            .push_back(
                &[7, 6, 5, 4, 3, 2, 1, 0],
                None,
                None,
                11,
                crate::basic::value::SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();

        let p0 = target_block.seqs().position(0, 3) as u64;
        let p1 = target_block.seqs().position(1, 2) as u64;
        let hits = vec![
            SearchHit::with_score(0, p0, 3, 12),
            SearchHit::with_score(0, p0 + 1, 4, 20),
            SearchHit::with_score(0, p1, 2, 9),
        ];
        let grouped = get_query_hits(&hits, &target_block, false);
        assert_eq!(grouped.len(), 2);
        assert_eq!(grouped[0], Hit::new(10, 20, 0));
        assert_eq!(grouped[1], Hit::new(11, 9, 0));

        let mut hits = vec![SearchHit::with_score(0, p0, 3, 0)];
        let reextended = get_query_hits_reextend(
            0,
            &mut hits,
            std::slice::from_ref(&query),
            &target_block,
            1,
            false,
            20,
            &score_matrix,
        );
        assert_eq!(reextended.len(), 1);
        assert_eq!(reextended[0].oid, 10);
        assert!(reextended[0].score > 20);
    }

    #[test]
    fn test_extend_query_from_query_list_and_ranking_table() {
        let score_matrix = sm();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut target_block = Block::new();
        target_block
            .push_back(
                &query,
                None,
                None,
                7,
                crate::basic::value::SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let mut db2block_id = TargetMap::new();
        db2block_id.insert(7, 0);
        let query_comp = crate::stats::cbs::compute_composition(&query);
        let query_list = QueryList {
            query_block_id: 0,
            last_query_block_id: 0,
            targets: vec![QueryListTarget {
                database_id: 7,
                score: 40,
            }],
        };
        let cfg = GappedScoreConfig {
            max_target_seqs: 10,
            ..GappedScoreConfig::default()
        };
        let mut stats = Statistics::new();
        let matches = extend_query_from_query_list(
            &query_list,
            &db2block_id,
            std::slice::from_ref(&query),
            "q",
            query.len() as i32,
            &[],
            &query_comp,
            &mut target_block,
            &mut stats,
            &cfg,
            &UngappedStageConfig::default(),
            &score_matrix,
            HspValues::COORDS,
            |_, _| -1,
            |_, _| -1,
            score_matrix.gap_open(),
            score_matrix.gap_extend(),
            1,
            100,
        );
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].target_block_id, 0);

        let cfg = GappedScoreConfig {
            global_ranking_targets: 2,
            max_target_seqs: 10,
            ..GappedScoreConfig::default()
        };
        let ranking_table = vec![Hit::new(7, 40, 0), Hit::default()];
        let mut stats = Statistics::new();
        let matches = extend_query_from_ranking_table(
            0,
            &ranking_table,
            &db2block_id,
            std::slice::from_ref(&query),
            "q",
            query.len() as i32,
            &[],
            &query_comp,
            &mut target_block,
            &mut stats,
            &cfg,
            &UngappedStageConfig::default(),
            &score_matrix,
            HspValues::COORDS,
            |_, _| -1,
            |_, _| -1,
            score_matrix.gap_open(),
            score_matrix.gap_extend(),
            1,
            100,
        );
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].target_block_id, 0);
    }
}
