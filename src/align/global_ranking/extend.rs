//! Lifecycle orchestration from `align/global_ranking/extend.cpp`.
//!
//! Query-level alignment remains in the parent module.  Database and output
//! facilities that are global pointers in C++ are supplied by the explicit
//! [`GlobalRankingBackend`] adapter.

use std::collections::{BTreeMap, HashMap};
use std::io::Read;

use crate::basic::value::BlockId;
use crate::util::data_structures::BitVector;

use super::{fetch_query_targets, Hit, QueryList};

/// File-local C++ `db_filter`, delegated to the shared query-core facade.
pub fn db_filter(table: &[Hit], database_size: usize) -> BitVector {
    super::db_filter(table, database_size)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SequenceLoadConfig {
    pub reset_seqinfo: bool,
    pub sequences: bool,
    pub titles: bool,
    pub full_titles: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RankedTarget {
    pub block_id: BlockId,
    pub score: u16,
    pub context: u8,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct QueryOutput {
    /// `None` is the C++ null `TextBuffer*`; it still advances the reorder queue.
    pub bytes: Option<Vec<u8>>,
    pub aligned: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtendConfig {
    pub threads: usize,
    pub threads_align: usize,
    pub target_masking: bool,
    pub titles_lazy: bool,
    pub full_titles: bool,
    pub iterated: bool,
    pub current_query_block: BlockId,
    pub current_reference_block: BlockId,
    pub query_contexts: usize,
    pub query_sequence_count: usize,
    pub global_ranking_targets: usize,
    pub ranking_table: Vec<Hit>,
    pub track_aligned_queries: bool,
    pub query_aligned: Vec<bool>,
    pub iteration_query_aligned: u32,
}

impl Default for ExtendConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            threads_align: 0,
            target_masking: false,
            titles_lazy: false,
            full_titles: false,
            iterated: false,
            current_query_block: 0,
            current_reference_block: 0,
            query_contexts: 1,
            query_sequence_count: 0,
            global_ranking_targets: 0,
            ranking_table: Vec::new(),
            track_aligned_queries: false,
            query_aligned: Vec::new(),
            iteration_query_aligned: 0,
        }
    }
}

pub trait GlobalRankingBackend {
    fn database_sequence_count(&self) -> usize;
    /// Returns loaded target OIDs in target-block order.
    fn load_target_sequences(
        &mut self,
        filter: &BitVector,
        config: SequenceLoadConfig,
    ) -> Result<Vec<u64>, String>;
    fn mask_target_sequences(&mut self) -> Result<(), String>;
    fn init_dictionary_block(&mut self, begin: BlockId, count: BlockId) -> Result<(), String>;
    fn close_dictionary_block(&mut self) -> Result<(), String>;
    fn init_random_access(
        &mut self,
        query_block: BlockId,
        reference_block: BlockId,
    ) -> Result<(), String>;
    fn end_random_access(&mut self) -> Result<(), String>;
    fn align_query(
        &mut self,
        query: BlockId,
        targets: &[RankedTarget],
        iterated: bool,
    ) -> Result<QueryOutput, String>;
    fn write_output(&mut self, bytes: &[u8]) -> Result<(), String>;
    fn merge_worker_statistics(&mut self, _worker: usize) {}
    fn clear_target_sequences(&mut self);
    fn delete_merged_query_list(&mut self) -> Result<(), String> {
        Ok(())
    }
    fn finish_intermediate_file(&mut self) -> Result<(), String> {
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub(crate) struct QueryJob {
    query: BlockId,
    targets: Vec<RankedTarget>,
}

fn target_map(oids: &[u64]) -> HashMap<u64, BlockId> {
    let mut map = HashMap::with_capacity(oids.len());
    for (block_id, &oid) in oids.iter().enumerate() {
        map.insert(oid, block_id as BlockId);
    }
    map
}

fn prepare_targets<B: GlobalRankingBackend>(
    backend: &mut B,
    filter: &BitVector,
    load: SequenceLoadConfig,
    masking: bool,
) -> Result<HashMap<u64, BlockId>, String> {
    let oids = backend.load_target_sequences(filter, load)?;
    if oids.len() > BlockId::MAX as usize {
        return Err("Loaded target count exceeds BlockId".to_owned());
    }
    let map = target_map(&oids);
    if masking {
        backend.mask_target_sequences()?;
    }
    Ok(map)
}

fn extend_query_from_query_list(
    query: &QueryList,
    db2block_id: &HashMap<u64, BlockId>,
) -> Result<QueryJob, String> {
    let mut targets = Vec::with_capacity(query.targets.len());
    for target in &query.targets {
        let block_id = db2block_id
            .get(&(target.database_id as u64))
            .copied()
            .ok_or_else(|| format!("Ranked target OID was not loaded: {}", target.database_id))?;
        targets.push(RankedTarget {
            block_id,
            score: target.score,
            context: 0,
        });
    }
    Ok(QueryJob {
        query: query.query_block_id,
        targets,
    })
}

fn extend_query_from_ranking_table(
    query: BlockId,
    hits: &[Hit],
    db2block_id: &HashMap<u64, BlockId>,
) -> Result<QueryJob, String> {
    let mut targets = Vec::with_capacity(hits.len());
    for hit in hits {
        let block_id = db2block_id
            .get(&(hit.oid as u64))
            .copied()
            .ok_or_else(|| format!("Ranked target OID was not loaded: {}", hit.oid))?;
        targets.push(RankedTarget {
            block_id,
            score: hit.score,
            context: hit.context,
        });
    }
    Ok(QueryJob { query, targets })
}

fn push_ordered<B: GlobalRankingBackend>(
    backend: &mut B,
    pending: &mut BTreeMap<BlockId, Option<Vec<u8>>>,
    next_output: &mut BlockId,
    query: BlockId,
    output: Option<Vec<u8>>,
) -> Result<(), String> {
    if pending.insert(query, output).is_some() {
        return Err(format!("Duplicate query output: {query}"));
    }
    while let Some(output) = pending.remove(next_output) {
        if let Some(bytes) = output {
            backend.write_output(&bytes)?;
        }
        *next_output = next_output
            .checked_add(1)
            .ok_or_else(|| "Query output index overflow".to_owned())?;
    }
    Ok(())
}

/// C++ `align_worker`. Jobs are partitioned by worker to model the atomic
/// fetch schedule; [`push_ordered`] retains `ReorderQueue` output semantics.
pub(crate) fn align_worker<B: GlobalRankingBackend>(
    backend: &mut B,
    config: &mut ExtendConfig,
    jobs: &[QueryJob],
    worker: usize,
    worker_count: usize,
    intermediate_output: bool,
    count_new_alignments: bool,
    pending: &mut BTreeMap<BlockId, Option<Vec<u8>>>,
    next_output: &mut BlockId,
) -> Result<(), String> {
    for job in jobs.iter().skip(worker).step_by(worker_count) {
        let output = backend.align_query(job.query, &job.targets, intermediate_output)?;
        if output.aligned && config.track_aligned_queries {
            let aligned = config
                .query_aligned
                .get_mut(job.query as usize)
                .ok_or_else(|| format!("Query block ID out of range: {}", job.query))?;
            if !*aligned {
                *aligned = true;
                if count_new_alignments {
                    config.iteration_query_aligned = config
                        .iteration_query_aligned
                        .checked_add(1)
                        .ok_or_else(|| "Aligned query count overflow".to_owned())?;
                }
            }
        }
        push_ordered(backend, pending, next_output, job.query, output.bytes)?;
    }
    backend.merge_worker_statistics(worker);
    Ok(())
}

fn run_workers<B: GlobalRankingBackend>(
    backend: &mut B,
    config: &mut ExtendConfig,
    jobs: &[QueryJob],
    initial_null: &[(BlockId, Option<Vec<u8>>)],
    intermediate_output: bool,
    count_new_alignments: bool,
) -> Result<(), String> {
    let worker_count = if config.threads_align != 0 {
        config.threads_align
    } else {
        config.threads
    };
    if worker_count == 0 {
        return Ok(());
    }
    let mut pending = BTreeMap::new();
    let mut next_output = 0;
    for &(query, ref output) in initial_null {
        push_ordered(
            backend,
            &mut pending,
            &mut next_output,
            query,
            output.clone(),
        )?;
    }
    for worker in 0..worker_count {
        align_worker(
            backend,
            config,
            jobs,
            worker,
            worker_count,
            intermediate_output,
            count_new_alignments,
            &mut pending,
            &mut next_output,
        )?;
    }
    if !pending.is_empty() {
        return Err(format!("Missing query output before query {next_output}"));
    }
    Ok(())
}

/// C++ `extend(SequenceFile&, TempFile&, ...)`.
pub fn extend_merged_query_list<R: Read, B: GlobalRankingBackend>(
    backend: &mut B,
    query_list: &mut R,
    ranking_db_filter: &BitVector,
    config: &mut ExtendConfig,
) -> Result<(), String> {
    let db2block_id = prepare_targets(
        backend,
        ranking_db_filter,
        SequenceLoadConfig {
            reset_seqinfo: true,
            sequences: true,
            titles: false,
            full_titles: false,
        },
        config.target_masking,
    )?;
    let mut next_query = 0;
    let mut jobs = Vec::new();
    let mut null_outputs = Vec::new();
    loop {
        let query = fetch_query_targets(query_list, &mut next_query).map_err(|e| e.to_string())?;
        if query.targets.is_empty() {
            break;
        }
        for skipped in query.last_query_block_id..query.query_block_id {
            null_outputs.push((skipped, None));
        }
        jobs.push(extend_query_from_query_list(&query, &db2block_id)?);
    }
    run_workers(backend, config, &jobs, &null_outputs, false, false)?;
    backend.delete_merged_query_list()?;
    backend.clear_target_sequences();
    Ok(())
}

/// C++ `extend(Search::Config&, Consumer&)`.
pub fn extend<B: GlobalRankingBackend>(
    backend: &mut B,
    config: &mut ExtendConfig,
) -> Result<(), String> {
    if config.query_contexts == 0 {
        return Err("Query context count must not be zero".to_owned());
    }
    let database_count = backend.database_sequence_count();
    let filter = db_filter(&config.ranking_table, database_count);
    let db2block_id = prepare_targets(
        backend,
        &filter,
        SequenceLoadConfig {
            reset_seqinfo: true,
            sequences: true,
            titles: !config.titles_lazy,
            full_titles: config.full_titles,
        },
        config.target_masking,
    )?;
    let target_count = BlockId::try_from(db2block_id.len())
        .map_err(|_| "Loaded target count exceeds BlockId".to_owned())?;
    if config.iterated {
        config.current_reference_block = 0;
        backend.init_dictionary_block(0, target_count)?;
    } else {
        backend.init_random_access(config.current_query_block, 0)?;
    }

    let query_count = config.query_sequence_count / config.query_contexts;
    let expected_table = query_count
        .checked_mul(config.global_ranking_targets)
        .ok_or_else(|| "Ranking table size overflow".to_owned())?;
    if config.ranking_table.len() < expected_table {
        return Err("Ranking table is shorter than the configured query count".to_owned());
    }
    let mut jobs = Vec::with_capacity(query_count);
    let mut null_outputs = Vec::new();
    for query in 0..query_count {
        let begin = query * config.global_ranking_targets;
        let end = begin + config.global_ranking_targets;
        let mut used_end = end;
        while used_end > begin && config.ranking_table[used_end - 1].score == 0 {
            used_end -= 1;
        }
        if used_end == begin {
            null_outputs.push((query as BlockId, None));
            continue;
        }
        jobs.push(extend_query_from_ranking_table(
            query as BlockId,
            &config.ranking_table[begin..used_end],
            &db2block_id,
        )?);
    }
    let iterated = config.iterated;
    run_workers(backend, config, &jobs, &null_outputs, iterated, iterated)?;
    backend.clear_target_sequences();
    if config.iterated {
        backend.close_dictionary_block()?;
        backend.finish_intermediate_file()?;
    } else {
        backend.end_random_access()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct FakeBackend {
        db_count: usize,
        loaded_oids: Vec<u64>,
        load: Vec<SequenceLoadConfig>,
        loaded_filter: Vec<usize>,
        calls: Vec<String>,
        align_order: Vec<BlockId>,
        output: Vec<u8>,
    }

    impl GlobalRankingBackend for FakeBackend {
        fn database_sequence_count(&self) -> usize {
            self.db_count
        }
        fn load_target_sequences(
            &mut self,
            filter: &BitVector,
            config: SequenceLoadConfig,
        ) -> Result<Vec<u64>, String> {
            self.load.push(config);
            self.loaded_filter = (0..filter.size() as usize)
                .filter(|&i| filter.get(i))
                .collect();
            self.calls.push("load".into());
            Ok(self.loaded_oids.clone())
        }
        fn mask_target_sequences(&mut self) -> Result<(), String> {
            self.calls.push("mask".into());
            Ok(())
        }
        fn init_dictionary_block(&mut self, begin: BlockId, count: BlockId) -> Result<(), String> {
            self.calls.push(format!("dict:{begin}:{count}"));
            Ok(())
        }
        fn close_dictionary_block(&mut self) -> Result<(), String> {
            self.calls.push("close-dict".into());
            Ok(())
        }
        fn init_random_access(&mut self, query: BlockId, reference: BlockId) -> Result<(), String> {
            self.calls.push(format!("random:{query}:{reference}"));
            Ok(())
        }
        fn end_random_access(&mut self) -> Result<(), String> {
            self.calls.push("end-random".into());
            Ok(())
        }
        fn align_query(
            &mut self,
            query: BlockId,
            targets: &[RankedTarget],
            iterated: bool,
        ) -> Result<QueryOutput, String> {
            self.align_order.push(query);
            let mut bytes = vec![b'0' + query as u8, b':', b'0' + targets[0].block_id as u8];
            if iterated {
                bytes.push(b'i');
            }
            Ok(QueryOutput {
                bytes: Some(bytes),
                aligned: targets[0].score >= 10,
            })
        }
        fn write_output(&mut self, bytes: &[u8]) -> Result<(), String> {
            self.output.extend_from_slice(bytes);
            Ok(())
        }
        fn merge_worker_statistics(&mut self, worker: usize) {
            self.calls.push(format!("stats:{worker}"));
        }
        fn clear_target_sequences(&mut self) {
            self.calls.push("clear".into());
        }
        fn delete_merged_query_list(&mut self) -> Result<(), String> {
            self.calls.push("delete".into());
            Ok(())
        }
        fn finish_intermediate_file(&mut self) -> Result<(), String> {
            self.calls.push("finish".into());
            Ok(())
        }
    }

    fn query_record(query: u32, targets: &[(u32, u16)]) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&query.to_ne_bytes());
        bytes.extend_from_slice(&((targets.len() * 6) as u32).to_ne_bytes());
        for &(oid, score) in targets {
            bytes.extend_from_slice(&oid.to_ne_bytes());
            bytes.extend_from_slice(&score.to_ne_bytes());
        }
        bytes
    }

    #[test]
    fn merged_workflow_maps_oids_fills_gaps_and_deletes_input() {
        let mut input = query_record(1, &[(9, 8)]);
        input.extend(query_record(3, &[(7, 12)]));
        let mut filter = BitVector::with_size(10);
        filter.set(7);
        filter.set(9);
        let mut backend = FakeBackend {
            db_count: 10,
            loaded_oids: vec![7, 9],
            ..Default::default()
        };
        let mut config = ExtendConfig {
            threads: 2,
            target_masking: true,
            track_aligned_queries: true,
            query_aligned: vec![false; 4],
            ..Default::default()
        };
        extend_merged_query_list(&mut backend, &mut input.as_slice(), &filter, &mut config)
            .unwrap();
        assert_eq!(backend.align_order, vec![1, 3]);
        assert_eq!(backend.output, b"1:13:0");
        assert_eq!(
            backend.calls,
            vec!["load", "mask", "stats:0", "stats:1", "delete", "clear"]
        );
        assert_eq!(config.query_aligned, vec![false, false, false, true]);
        assert_eq!(
            backend.load[0],
            SequenceLoadConfig {
                reset_seqinfo: true,
                sequences: true,
                titles: false,
                full_titles: false
            }
        );
    }

    #[test]
    fn ranking_workflow_reorders_workers_and_uses_random_access() {
        let mut backend = FakeBackend {
            db_count: 8,
            loaded_oids: vec![2, 5],
            ..Default::default()
        };
        let mut config = ExtendConfig {
            threads: 2,
            titles_lazy: false,
            full_titles: true,
            current_query_block: 4,
            query_sequence_count: 3,
            global_ranking_targets: 2,
            ranking_table: vec![
                Hit::new(2, 12, 1),
                Hit::default(),
                Hit::new(5, 9, 2),
                Hit::default(),
                Hit::new(2, 11, 3),
                Hit::default(),
            ],
            track_aligned_queries: true,
            query_aligned: vec![false; 3],
            ..Default::default()
        };
        extend(&mut backend, &mut config).unwrap();
        assert_eq!(backend.loaded_filter, vec![2, 5]);
        assert_eq!(backend.align_order, vec![0, 2, 1]);
        assert_eq!(backend.output, b"0:01:12:0");
        assert!(backend.calls.contains(&"random:4:0".to_owned()));
        assert_eq!(backend.calls.last().unwrap(), "end-random");
        assert_eq!(config.query_aligned, vec![true, false, true]);
        assert!(backend.load[0].titles && backend.load[0].full_titles);
    }

    #[test]
    fn iterated_workflow_tracks_new_alignments_and_finishes_file() {
        let mut backend = FakeBackend {
            db_count: 4,
            loaded_oids: vec![1],
            ..Default::default()
        };
        let mut config = ExtendConfig {
            iterated: true,
            query_sequence_count: 2,
            global_ranking_targets: 1,
            ranking_table: vec![Hit::new(1, 10, 4), Hit::default()],
            track_aligned_queries: true,
            query_aligned: vec![false; 2],
            ..Default::default()
        };
        extend(&mut backend, &mut config).unwrap();
        assert_eq!(config.iteration_query_aligned, 1);
        assert_eq!(backend.output, b"0:0i");
        assert!(backend
            .calls
            .windows(3)
            .any(|w| w == ["clear", "close-dict", "finish"]));
    }

    #[test]
    fn ranking_workflow_rejects_missing_target_mapping() {
        let mut backend = FakeBackend {
            db_count: 4,
            loaded_oids: vec![],
            ..Default::default()
        };
        let mut config = ExtendConfig {
            query_sequence_count: 1,
            global_ranking_targets: 1,
            ranking_table: vec![Hit::new(3, 1, 0)],
            ..Default::default()
        };
        assert_eq!(
            extend(&mut backend, &mut config).unwrap_err(),
            "Ranked target OID was not loaded: 3"
        );
    }
}
