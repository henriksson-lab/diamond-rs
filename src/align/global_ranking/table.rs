//! Translation of `diamond/src/align/global_ranking/table.cpp`.
//!
//! The table helpers were originally translated into the parent module. They
//! remain re-exported there for API compatibility; this mirrored source file
//! owns the previously missing `update_table` orchestration.

use rayon::prelude::*;

use super::Hit;
pub use super::{get_query_hits, get_query_hits_reextend, merge_hits, target_score};
use crate::basic::value::Letter;
use crate::data::block::Block;
use crate::search::hit::Hit as SearchHit;
use crate::stats::score_matrix::ScoreMatrix;

/// Explicit replacement for the C++ globals read by `update_table`:
/// `config`, `align_mode`, and the query/target members of `Search::Config`.
#[derive(Clone, Copy)]
pub struct UpdateTableConfig<'a> {
    pub query_contexts: u32,
    pub global_ranking_targets: usize,
    pub threads: usize,
    pub no_reextend: bool,
    pub xdrop: i32,
    /// Query contexts stored consecutively for each source query.
    pub query_sequences: &'a [Vec<Letter>],
    pub target_block: &'a Block,
    pub score_matrix: &'a ScoreMatrix,
}

/// C++ `Extension::GlobalRanking::update_table(Search::Config&)`.
///
/// On success, a non-empty ranking buffer is consumed just as the C++ code
/// resets `global_ranking_buffer` after processing. An empty buffer is retained
/// because the original returns before the reset in that case. The return
/// value is the C++ diagnostic `merged_count`.
pub fn update_table(
    global_ranking_buffer: &mut Option<Vec<SearchHit>>,
    ranking_table: &mut [Hit],
    config: &UpdateTableConfig<'_>,
) -> Result<usize, String> {
    if config.query_contexts == 0 {
        return Err("update_table: query_contexts must be greater than zero".to_string());
    }
    let Some(seed_hits) = global_ranking_buffer.as_mut() else {
        return Ok(0);
    };
    if seed_hits.is_empty() {
        return Ok(0);
    }

    seed_hits.sort_by(SearchHit::cmp_query_target);
    let contexts = config.query_contexts as usize;

    // C++ uses `partition_table` only to distribute complete source-query key
    // ranges. Construct those exact ranges first, then process them in the
    // configured Rayon pool. Results are merged in query order; every query
    // owns a disjoint ranking-table slice.
    let mut groups = Vec::new();
    let mut begin = 0usize;
    while begin < seed_hits.len() {
        let query = seed_hits[begin].query as usize / contexts;
        let mut end = begin + 1;
        while end < seed_hits.len() && seed_hits[end].query as usize / contexts == query {
            end += 1;
        }
        let table_end = query
            .checked_add(1)
            .and_then(|n| n.checked_mul(config.global_ranking_targets))
            .ok_or_else(|| "update_table: ranking table index overflow".to_string())?;
        if table_end > ranking_table.len() {
            return Err(format!(
                "update_table: ranking table too small for query {query}"
            ));
        }
        let query_end = query
            .checked_add(1)
            .and_then(|n| n.checked_mul(contexts))
            .ok_or_else(|| "update_table: query sequence index overflow".to_string())?;
        if query_end > config.query_sequences.len() {
            return Err(format!(
                "update_table: missing sequence contexts for query {query}"
            ));
        }
        groups.push((query, begin, end));
        begin = end;
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads.max(1))
        .build()
        .map_err(|error| format!("update_table: could not create worker pool: {error}"))?;
    let query_hits: Vec<(usize, Vec<Hit>)> = pool.install(|| {
        groups
            .par_iter()
            .map(|&(query, begin, end)| {
                let mut hits = seed_hits[begin..end].to_vec();
                let query_begin = query * contexts;
                let query_end = query_begin + contexts;
                let ranked = get_query_hits_reextend(
                    query as u32,
                    &mut hits,
                    &config.query_sequences[query_begin..query_end],
                    config.target_block,
                    config.query_contexts,
                    config.no_reextend,
                    config.xdrop,
                    config.score_matrix,
                );
                (query, ranked)
            })
            .collect()
    });

    let mut merged_count = 0usize;
    let mut merged = Vec::new();
    for (query, mut hits) in query_hits {
        merge_hits(
            query,
            &mut hits,
            &mut merged,
            ranking_table,
            config.global_ranking_targets,
            &mut merged_count,
        );
    }
    *global_ranking_buffer = None;
    Ok(merged_count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SequenceType;

    fn score_matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 100).unwrap()
    }

    fn target_block() -> Block {
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2, 3, 4, 5],
                None,
                None,
                10,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        block
            .push_back(
                &[5, 4, 3, 2, 1, 0],
                None,
                None,
                11,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        block
    }

    #[test]
    fn update_table_groups_contexts_merges_and_consumes_buffer() {
        let target = target_block();
        let matrix = score_matrix();
        let target0 = target.seqs().position(0, 2) as u64;
        let target1 = target.seqs().position(1, 3) as u64;
        // Contexts 0 and 1 belong to source query 0; contexts 2 and 3 to
        // source query 1. Deliberately unsorted to exercise the initial sort.
        let mut buffer = Some(vec![
            SearchHit::with_score(2, target0, 1, 30),
            SearchHit::with_score(1, target1, 2, 25),
            SearchHit::with_score(0, target0, 2, 20),
        ]);
        let mut table = vec![Hit::default(); 4];
        table[0] = Hit::new(11, 40, 0);
        let queries = vec![
            vec![0, 1, 2, 3, 4, 5],
            vec![5, 4, 3, 2, 1, 0],
            vec![0, 1, 2, 3, 4, 5],
            vec![5, 4, 3, 2, 1, 0],
        ];
        let config = UpdateTableConfig {
            query_contexts: 2,
            global_ranking_targets: 2,
            threads: 2,
            no_reextend: true,
            xdrop: 20,
            query_sequences: &queries,
            target_block: &target,
            score_matrix: &matrix,
        };

        assert_eq!(update_table(&mut buffer, &mut table, &config).unwrap(), 2);
        assert!(buffer.is_none());
        assert_eq!(table[0], Hit::new(11, 40, 0));
        assert_eq!(table[1], Hit::new(10, 20, 0));
        assert_eq!(table[2], Hit::new(10, 30, 0));
        assert_eq!(table[3], Hit::default());
    }

    #[test]
    fn update_table_retains_empty_buffer_like_cpp_early_return() {
        let target = target_block();
        let matrix = score_matrix();
        let queries = vec![vec![0, 1, 2, 3, 4, 5]];
        let config = UpdateTableConfig {
            query_contexts: 1,
            global_ranking_targets: 1,
            threads: 1,
            no_reextend: true,
            xdrop: 20,
            query_sequences: &queries,
            target_block: &target,
            score_matrix: &matrix,
        };
        let mut buffer = Some(Vec::new());
        let mut table = vec![Hit::default()];
        assert_eq!(update_table(&mut buffer, &mut table, &config).unwrap(), 0);
        assert_eq!(buffer, Some(Vec::new()));
    }

    #[test]
    fn update_table_validates_global_state_replacements_before_consuming() {
        let target = target_block();
        let matrix = score_matrix();
        let config = UpdateTableConfig {
            query_contexts: 0,
            global_ranking_targets: 1,
            threads: 1,
            no_reextend: true,
            xdrop: 20,
            query_sequences: &[],
            target_block: &target,
            score_matrix: &matrix,
        };
        let original = vec![SearchHit::with_score(0, 0, 0, 1)];
        let mut buffer = Some(original.clone());
        let error = update_table(&mut buffer, &mut [], &config).unwrap_err();
        assert!(error.contains("query_contexts"));
        assert_eq!(buffer, Some(original));
    }
}
