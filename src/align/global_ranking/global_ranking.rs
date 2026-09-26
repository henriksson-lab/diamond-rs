//! Global-ranking list primitives mirrored from
//! `diamond/src/align/global_ranking/global_ranking.cpp` and its header.

use crate::align::gapped_filter::{SeedHit, TargetScore};
use crate::align::hsp::Match;
use crate::basic::statistics::Statistics;
use crate::basic::value::{BlockId, Letter, DELIMITER_LETTER};
use crate::data::block::Block;
use crate::dp::ungapped::ungapped_window;
use crate::output::intermediate::IntermediateRecord;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::data_structures::{BitVector, FlatArray};
use crate::util::text_buffer::TextBuffer;
use std::io::{self, Read};

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

fn padded_subject_window(
    target: &[Letter],
    begin: isize,
    len: usize,
) -> std::borrow::Cow<'_, [Letter]> {
    let end = begin.saturating_add(len as isize);
    if begin >= 0 && end <= target.len() as isize {
        return std::borrow::Cow::Borrowed(&target[begin as usize..end as usize]);
    }
    let mut out = vec![DELIMITER_LETTER; len];
    let source_begin = begin.max(0) as usize;
    let source_end = end.max(0).min(target.len() as isize) as usize;
    if source_begin < source_end {
        let destination_begin = (source_begin as isize - begin) as usize;
        out[destination_begin..destination_begin + source_end - source_begin]
            .copy_from_slice(&target[source_begin..source_end]);
    }
    std::borrow::Cow::Owned(out)
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
        let center = hit.i as isize;
        let query_begin = (center - ungapped_window_size as isize).max(0) as usize;
        let query_end = (center + ungapped_window_size as isize)
            .max(0)
            .min(query_seq.len() as isize) as usize;
        if query_begin >= query_end {
            continue;
        }
        let window_left = center - query_begin as isize;
        let subject_begin = hit.j as isize - window_left;
        let subject = padded_subject_window(target_seq, subject_begin, query_end - query_begin);
        let s = ungapped_window(
            &query_seq[query_begin..query_end],
            subject.as_ref(),
            query_end - query_begin,
            score_matrix,
        );
        score = score.max(s);
    }
    score.min(u16::MAX as i32) as u16
}

/// Matches C++ `GlobalRanking::ranking_list`.
#[allow(clippy::too_many_arguments)]
pub fn ranking_list(
    _query_id: BlockId,
    target_scores: &mut [TargetScore],
    target_block_ids: &[BlockId],
    seed_hits: &FlatArray<SeedHit>,
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

/// The upstream body is intentionally commented out and has no side effects.
pub fn write_merged_query_list(
    _record: &IntermediateRecord,
    _out: &mut TextBuffer,
    _ranking_db_filter: &mut BitVector,
    _stat: &mut Statistics,
) {
}

/// Matches C++ `GlobalRanking::finish_merged_query_list`.
pub fn finish_merged_query_list(buf: &mut TextBuffer, seek_pos: usize) {
    let payload_size = buf
        .size()
        .checked_sub(seek_pos)
        .and_then(|size| size.checked_sub(std::mem::size_of::<u32>() * 2))
        .expect("merged query-list header lies outside the buffer");
    let payload_size =
        u32::try_from(payload_size).expect("safe_cast: overflow (unsigned -> unsigned)");
    let begin = seek_pos + std::mem::size_of::<u32>();
    buf.data_mut()[begin..begin + std::mem::size_of::<u32>()]
        .copy_from_slice(&payload_size.to_ne_bytes());
}

/// Matches C++ `GlobalRanking::fetch_query_targets`. Rust's exclusive reader
/// and `next_query` borrows replace the upstream process-global mutex.
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
        Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => return Ok(out),
        Err(error) => return Err(error),
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

#[cfg(test)]
mod tests {
    use super::*;

    fn score_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn overflow_window_stays_centered_when_query_left_edge_is_clipped() {
        let matrix = score_matrix();
        let query = [0; 8];
        let target = [0; 8];
        let hit = SeedHit::new(1, 1, u8::MAX as i32, 0);

        let score = recompute_overflow_scores(&[hit], &query, &target, 4, &matrix);
        // C++ starts at i-window and clips at the sequence delimiter, leaving
        // [0, i+window), not a fresh two-window span starting at zero.
        let expected = ungapped_window(&query[..5], &target[..5], 5, &matrix) as u16;
        assert_eq!(score, expected);
        assert_ne!(
            score,
            ungapped_window(&query, &target, query.len(), &matrix) as u16
        );
    }

    #[test]
    fn overflow_window_preserves_target_perimeter_padding() {
        let matrix = score_matrix();
        let query = [0; 8];
        let target = [0; 8];
        let hit = SeedHit::new(1, 0, u8::MAX as i32, 0);
        let score = recompute_overflow_scores(&[hit], &query, &target, 4, &matrix);

        let padded = [DELIMITER_LETTER, 0, 0, 0, 0];
        let expected = ungapped_window(&query[..5], &padded, 5, &matrix) as u16;
        assert_eq!(score, expected);
        assert_eq!(
            recompute_overflow_scores(&[SeedHit::new(1, 0, 254, 0)], &query, &target, 4, &matrix,),
            0
        );
    }

    #[test]
    #[should_panic(expected = "merged query-list header lies outside the buffer")]
    fn finish_merged_query_list_rejects_invalid_header_position() {
        finish_merged_query_list(&mut TextBuffer::new(), 1);
    }

    #[test]
    fn fetch_query_targets_rejects_truncated_target_payload() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&3u32.to_ne_bytes());
        bytes.extend_from_slice(&6u32.to_ne_bytes());
        bytes.extend_from_slice(&9u32.to_ne_bytes());
        let mut input = bytes.as_slice();
        let mut next_query = 0;
        let error = fetch_query_targets(&mut input, &mut next_query).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
        // Upstream updates this immediately after reading the query id, before
        // attempting the size/payload reads.
        assert_eq!(next_query, 4);
    }
}
