//! Top-level alignment orchestration from `diamond/src/align/align.cpp`.
//!
//! The original owns process-global hit buffers, output queues, statistics,
//! aligned-query bitsets, database random access, and worker pools.  This
//! module keeps their ordering explicit through [`AlignBackend`] and
//! [`AlignState`].

use std::ops::Range;

use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::value::{BlockId, OId};
use crate::dp::banded_3frame::DpStat;
use crate::search::hit::Hit;

pub type OutputBuffer = Vec<u8>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoadBalancing {
    QueryParallel,
    TargetParallel,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignConfig {
    pub min_task_trace_pts: usize,
    pub query_contexts: u32,
    pub swipe_all: bool,
    pub blastn: bool,
    pub frame_shift: i32,
    pub blocked_processing: bool,
    pub output_is_daa: bool,
    pub output_is_null: bool,
    pub report_unaligned: bool,
    pub query_separator: u8,
    pub track_aligned_queries: bool,
    pub track_aligned_targets: bool,
    pub load_balancing: LoadBalancing,
    pub target_sequence_count: usize,
    pub query_sequence_count: usize,
    pub threads: usize,
    pub threads_align: usize,
    pub verbosity: u32,
    pub heartbeat: bool,
    pub iterated: bool,
    pub current_query_block: BlockId,
    pub memory_limit: u64,
    pub trace_pt_fetch_size: u64,
    /// Native C++ `sizeof(Search::Hit)`, supplied explicitly because Rust's
    /// struct layout is not the packed on-disk C++ layout.
    pub hit_record_size: u64,
    pub global_ranking_targets: i64,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            min_task_trace_pts: 1,
            query_contexts: 1,
            swipe_all: false,
            blastn: false,
            frame_shift: 0,
            blocked_processing: false,
            output_is_daa: false,
            output_is_null: false,
            report_unaligned: false,
            query_separator: b'\n',
            track_aligned_queries: false,
            track_aligned_targets: false,
            load_balancing: LoadBalancing::QueryParallel,
            target_sequence_count: 0,
            query_sequence_count: 0,
            threads: 1,
            threads_align: 0,
            verbosity: 0,
            heartbeat: false,
            iterated: false,
            current_query_block: 0,
            memory_limit: 16 * 1024 * 1024 * 1024,
            trace_pt_fetch_size: u64::MAX,
            hit_record_size: 15,
            global_ranking_targets: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HitBatch {
    pub hits: Vec<Hit>,
    pub query_begin: BlockId,
    pub query_end: BlockId,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AlignmentMatches {
    pub target_block_ids: Vec<BlockId>,
}

impl AlignmentMatches {
    pub fn is_empty(&self) -> bool {
        self.target_block_ids.is_empty()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hits {
    pub query: BlockId,
    pub range: Option<Range<usize>>,
}

/// Explicit adapter for globals and heavyweight extension/output subsystems.
pub trait AlignBackend {
    fn query_memory_size(&self) -> u64;
    fn target_memory_size(&self) -> u64;
    fn allocate_hit_buffer(&mut self);
    fn load_hit_buffer(&mut self, max_bytes: u64) -> bool;
    fn retrieve_hit_buffer(&mut self) -> Result<HitBatch, String>;
    fn next_bin_size(&self) -> u64;
    fn hit_buffer_disk_size(&self) -> i64;
    fn free_hit_buffer(&mut self);

    fn init_random_access(&mut self, query_block: BlockId) -> Result<(), String>;
    fn end_random_access(&mut self) -> Result<(), String>;

    fn unaligned_query_output(&mut self, query: BlockId) -> OutputBuffer;
    /// Runs the already-translated legacy mapper/pipeline and its output pass.
    fn run_legacy_mapper(
        &mut self,
        query: BlockId,
        hits: &[Hit],
        stat: &mut Statistics,
        dp_stat: &mut DpStat,
    ) -> Result<(OutputBuffer, bool), String>;

    fn extend(
        &mut self,
        query: BlockId,
        hits: Option<&[Hit]>,
        stat: &mut Statistics,
        dp_stat: &mut DpStat,
        parallel: bool,
        dna: bool,
        global_ranking_targets: i64,
    ) -> Result<AlignmentMatches, String>;
    fn generate_output(
        &mut self,
        matches: &AlignmentMatches,
        query: BlockId,
        stat: &mut Statistics,
    ) -> Result<OutputBuffer, String>;
    fn generate_intermediate_output(
        &mut self,
        matches: &AlignmentMatches,
        query: BlockId,
    ) -> Result<OutputBuffer, String>;
    fn target_oid(&self, target_block_id: BlockId) -> OId;
    fn begin_output_batch(
        &mut self,
        _query_begin: BlockId,
        _blocked_processing: bool,
        _query_separator: u8,
    ) {
    }
    fn push_output(&mut self, query: BlockId, buffer: Option<OutputBuffer>);
    fn end_output_batch(&mut self) {}

    fn log(&mut self, _message: &str) {}
    fn elapsed_microseconds(&mut self, _stage: &'static str) -> i64 {
        0
    }
    fn worker_count(&mut self, _threads: usize, _tasks: Range<i64>) {}
    fn heartbeat_begin(&mut self, _query_end: BlockId) {}
    fn heartbeat_end(&mut self) {}
}

#[derive(Debug, Clone)]
pub struct AlignState {
    pub statistics: Statistics,
    pub dp_stat: DpStat,
    pub query_aligned: Vec<bool>,
    pub iteration_query_aligned: u64,
    pub aligned_targets: Vec<bool>,
}

impl AlignState {
    pub fn new(query_count: usize, target_count: usize) -> Self {
        Self {
            statistics: Statistics::default(),
            dp_stat: DpStat::default(),
            query_aligned: vec![false; query_count],
            iteration_query_aligned: 0,
            aligned_targets: vec![false; target_count],
        }
    }

    fn mark_query_aligned(&mut self, query: BlockId) {
        let aligned = &mut self.query_aligned[query as usize];
        if !*aligned {
            *aligned = true;
            self.iteration_query_aligned += 1;
        }
    }
}

/// C++ `make_partition`, with `config.min_task_trace_pts` and query contexts
/// made explicit.
pub fn make_partition(hits: &[Hit], min_task_trace_pts: usize, query_contexts: u32) -> Vec<i64> {
    assert!(min_task_trace_pts > 0 && query_contexts > 0);
    let mut partition = Vec::with_capacity(hits.len().div_ceil(min_task_trace_pts) + 1);
    partition.push(0);
    let mut p = 0usize;
    while p < hits.len() {
        let mut q = (p + min_task_trace_pts).min(hits.len() - 1);
        let query = hits[q].query / query_contexts;
        q += 1;
        while q < hits.len() && hits[q].query / query_contexts == query {
            q += 1;
        }
        partition.push(q as i64);
        p = q;
    }
    partition
}

pub struct HitIterator<'a> {
    pub partition: &'a [i64],
    pub parts: i64,
    pub data: &'a [Hit],
    pub query_begin: BlockId,
    pub query_end: BlockId,
    pub query_contexts: u32,
    pub swipe_all: bool,
    pub blastn: bool,
}

impl<'a> HitIterator<'a> {
    pub fn new(
        query_begin: BlockId,
        query_end: BlockId,
        data: &'a [Hit],
        partition: &'a [i64],
        query_contexts: u32,
        swipe_all: bool,
        blastn: bool,
    ) -> Self {
        Self {
            partition,
            parts: partition.len() as i64 - 1,
            data,
            query_begin,
            query_end,
            query_contexts,
            swipe_all,
            blastn,
        }
    }

    pub fn single_query(swipe_all: bool, blastn: bool) -> bool {
        swipe_all || blastn
    }

    pub fn fetch(&self, i: i64) -> Vec<Hits> {
        if Self::single_query(self.swipe_all, self.blastn) {
            return vec![Hits {
                query: i as BlockId,
                range: None,
            }];
        }
        assert!(i >= 0 && i < self.parts);
        let begin = self.partition[i as usize] as usize;
        let end = self.partition[i as usize + 1] as usize;
        let mut last_query = if begin > 0 {
            self.data[begin - 1].query / self.query_contexts + 1
        } else {
            self.query_begin
        };
        let mut result = Vec::new();
        while last_query < self.data[begin].query / self.query_contexts {
            result.push(Hits {
                query: last_query,
                range: None,
            });
            last_query += 1;
        }
        let mut cursor = begin;
        while cursor < end {
            let query = self.data[cursor].query / self.query_contexts;
            while last_query < query {
                result.push(Hits {
                    query: last_query,
                    range: None,
                });
                last_query += 1;
            }
            let group_begin = cursor;
            while cursor < end && self.data[cursor].query / self.query_contexts == query {
                cursor += 1;
            }
            result.push(Hits {
                query,
                range: Some(group_begin..cursor),
            });
            last_query = query + 1;
        }
        if i == self.parts - 1 {
            let mut query = result.last().expect("nonempty partition").query + 1;
            while query < self.query_end {
                result.push(Hits { query, range: None });
                query += 1;
            }
        }
        result
    }
}

pub fn legacy_pipeline<B: AlignBackend>(
    hit_group: &Hits,
    all_hits: &[Hit],
    cfg: &AlignConfig,
    backend: &mut B,
    state: &mut AlignState,
    stat: &mut Statistics,
    dp_stat: &mut DpStat,
) -> Result<Option<OutputBuffer>, String> {
    let hits = hit_group
        .range
        .as_ref()
        .map(|range| &all_hits[range.clone()]);
    if hits.is_none_or(<[Hit]>::is_empty) {
        if !cfg.blocked_processing && !cfg.output_is_daa && cfg.report_unaligned {
            return Ok(Some(backend.unaligned_query_output(hit_group.query)));
        }
        return Ok(None);
    }
    let (buffer, aligned) =
        backend.run_legacy_mapper(hit_group.query, hits.expect("checked above"), stat, dp_stat)?;
    if cfg.output_is_null {
        return Ok(None);
    }
    if aligned && cfg.track_aligned_queries {
        state.mark_query_aligned(hit_group.query);
    }
    Ok(Some(buffer))
}

pub fn align_worker<B: AlignBackend>(
    hit_it: &HitIterator<'_>,
    cfg: &AlignConfig,
    backend: &mut B,
    state: &mut AlignState,
    next: i64,
) -> Result<(), String> {
    let groups = hit_it.fetch(next);
    assert!(!groups.is_empty());
    let mut worker_stat = Statistics::default();
    let mut worker_dp_stat = DpStat::default();
    let parallel = cfg.swipe_all && cfg.target_sequence_count >= cfg.query_sequence_count;

    for group in groups {
        if cfg.frame_shift != 0 {
            let output = legacy_pipeline(
                &group,
                hit_it.data,
                cfg,
                backend,
                state,
                &mut worker_stat,
                &mut worker_dp_stat,
            )?;
            backend.push_output(group.query, output);
            continue;
        }
        if group.range.is_none() && !HitIterator::single_query(cfg.swipe_all, cfg.blastn) {
            backend.push_output(group.query, None);
            continue;
        }
        let hits = group
            .range
            .as_ref()
            .map(|range| &hit_it.data[range.clone()]);
        let matches = backend.extend(
            group.query,
            hits,
            &mut worker_stat,
            &mut worker_dp_stat,
            parallel,
            cfg.blastn,
            cfg.global_ranking_targets,
        )?;
        let output = if cfg.blocked_processing {
            backend.generate_intermediate_output(&matches, group.query)?
        } else {
            backend.generate_output(&matches, group.query, &mut worker_stat)?
        };
        if !matches.is_empty() && cfg.track_aligned_queries {
            state.mark_query_aligned(group.query);
        }
        if cfg.track_aligned_targets {
            for &target in &matches.target_block_ids {
                let oid = backend.target_oid(target) as usize;
                if !state.aligned_targets[oid] {
                    state.aligned_targets[oid] = true;
                }
            }
        }
        backend.push_output(group.query, Some(output));
    }
    state.statistics += &worker_stat;
    state.dp_stat.gross_cells += worker_dp_stat.gross_cells;
    state.dp_stat.net_cells += worker_dp_stat.net_cells;
    Ok(())
}

pub fn align_queries<B: AlignBackend>(
    cfg: &AlignConfig,
    backend: &mut B,
    state: &mut AlignState,
) -> Result<(), String> {
    if !cfg.blocked_processing && !cfg.iterated {
        backend.init_random_access(cfg.current_query_block)?;
    }
    let mut resident_size = backend.query_memory_size() + backend.target_memory_size();
    backend.allocate_hit_buffer();
    let initial_capacity = cfg
        .memory_limit
        .wrapping_sub(resident_size)
        .wrapping_sub(backend.next_bin_size().wrapping_mul(cfg.hit_record_size))
        .min(cfg.trace_pt_fetch_size);
    let _ = backend.load_hit_buffer(initial_capacity);
    let mut keep_going = true;

    while keep_going {
        let mut batch = backend.retrieve_hit_buffer()?;
        state.statistics.inc(
            StatValue::TimeLoadSeedHits,
            backend.elapsed_microseconds("Loading trace points"),
        );
        let next_capacity = cfg
            .memory_limit
            .wrapping_sub(resident_size)
            .wrapping_sub(backend.next_bin_size().wrapping_mul(cfg.hit_record_size))
            .min(cfg.trace_pt_fetch_size);
        keep_going = backend.load_hit_buffer(next_capacity);
        backend.log(&format!(
            "Processing {} trace points ({}).",
            batch.hits.len(),
            crate::util::string::format(
                (batch.hits.len() as u64).wrapping_mul(cfg.hit_record_size)
            )
        ));
        resident_size =
            resident_size.wrapping_add((batch.hits.len() as u64).wrapping_mul(cfg.hit_record_size));

        batch
            .hits
            .sort_unstable_by(|lhs, rhs| lhs.query.cmp(&rhs.query));
        state.statistics.inc(
            StatValue::TimeSortSeedHits,
            backend.elapsed_microseconds("Sorting trace points"),
        );
        let partition = make_partition(&batch.hits, cfg.min_task_trace_pts, cfg.query_contexts);
        let hit_it = HitIterator::new(
            batch.query_begin,
            batch.query_end,
            &batch.hits,
            &partition,
            cfg.query_contexts,
            cfg.swipe_all,
            cfg.blastn,
        );
        backend.begin_output_batch(
            batch.query_begin,
            cfg.blocked_processing,
            if cfg.blocked_processing {
                0
            } else {
                cfg.query_separator
            },
        );
        let parallel = cfg.swipe_all && cfg.target_sequence_count >= cfg.query_sequence_count;
        let threads = if cfg.load_balancing == LoadBalancing::TargetParallel || parallel {
            1
        } else if cfg.threads_align == 0 {
            cfg.threads
        } else {
            cfg.threads_align
        };
        let tasks = if cfg.swipe_all {
            batch.query_begin as i64..batch.query_end as i64
        } else {
            0..hit_it.parts
        };
        backend.worker_count(threads, tasks.clone());
        let heartbeat = cfg.verbosity >= 3
            && cfg.load_balancing == LoadBalancing::QueryParallel
            && !cfg.swipe_all
            && cfg.heartbeat;
        if heartbeat {
            backend.heartbeat_begin(batch.query_end);
        }
        for task in tasks {
            align_worker(&hit_it, cfg, backend, state, task)?;
        }
        if heartbeat {
            backend.heartbeat_end();
        }
        state.statistics.inc(
            StatValue::TimeExt,
            backend.elapsed_microseconds("Computing alignments"),
        );
        backend.end_output_batch();
    }
    state
        .statistics
        .max(StatValue::SearchTempSpace, backend.hit_buffer_disk_size());
    backend.free_hit_buffer();
    if !cfg.blocked_processing && !cfg.iterated {
        backend.end_random_access()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;

    use super::*;

    #[derive(Default)]
    struct Backend {
        batches: VecDeque<HitBatch>,
        load_results: VecDeque<bool>,
        capacities: Vec<u64>,
        outputs: Vec<(BlockId, Option<OutputBuffer>)>,
        output_batches: Vec<(BlockId, bool, u8)>,
        ended_output_batches: usize,
        extend_calls: Vec<(BlockId, bool, bool, bool, i64)>,
        legacy_calls: Vec<BlockId>,
        random_access: Vec<&'static str>,
        workers: Vec<(usize, Range<i64>)>,
        heartbeat: Vec<&'static str>,
        allocated: bool,
        freed: bool,
        next_bin: u64,
        legacy_aligned: bool,
    }

    impl AlignBackend for Backend {
        fn query_memory_size(&self) -> u64 {
            100
        }
        fn target_memory_size(&self) -> u64 {
            200
        }
        fn allocate_hit_buffer(&mut self) {
            self.allocated = true;
        }
        fn load_hit_buffer(&mut self, max_bytes: u64) -> bool {
            self.capacities.push(max_bytes);
            self.load_results.pop_front().unwrap_or(false)
        }
        fn retrieve_hit_buffer(&mut self) -> Result<HitBatch, String> {
            self.batches
                .pop_front()
                .ok_or_else(|| "no batch".to_string())
        }
        fn next_bin_size(&self) -> u64 {
            self.next_bin
        }
        fn hit_buffer_disk_size(&self) -> i64 {
            77
        }
        fn free_hit_buffer(&mut self) {
            self.freed = true;
        }
        fn init_random_access(&mut self, _query_block: BlockId) -> Result<(), String> {
            self.random_access.push("begin");
            Ok(())
        }
        fn end_random_access(&mut self) -> Result<(), String> {
            self.random_access.push("end");
            Ok(())
        }
        fn unaligned_query_output(&mut self, query: BlockId) -> OutputBuffer {
            vec![query as u8, b'U']
        }
        fn run_legacy_mapper(
            &mut self,
            query: BlockId,
            _hits: &[Hit],
            _stat: &mut Statistics,
            dp_stat: &mut DpStat,
        ) -> Result<(OutputBuffer, bool), String> {
            self.legacy_calls.push(query);
            dp_stat.net_cells += 1;
            Ok((vec![query as u8, b'L'], self.legacy_aligned))
        }
        fn extend(
            &mut self,
            query: BlockId,
            hits: Option<&[Hit]>,
            _stat: &mut Statistics,
            dp_stat: &mut DpStat,
            parallel: bool,
            dna: bool,
            global_ranking_targets: i64,
        ) -> Result<AlignmentMatches, String> {
            self.extend_calls
                .push((query, hits.is_none(), parallel, dna, global_ranking_targets));
            dp_stat.gross_cells += 2;
            Ok(AlignmentMatches {
                target_block_ids: vec![query],
            })
        }
        fn generate_output(
            &mut self,
            _matches: &AlignmentMatches,
            query: BlockId,
            _stat: &mut Statistics,
        ) -> Result<OutputBuffer, String> {
            Ok(vec![query as u8, b'O'])
        }
        fn generate_intermediate_output(
            &mut self,
            _matches: &AlignmentMatches,
            query: BlockId,
        ) -> Result<OutputBuffer, String> {
            Ok(vec![query as u8, b'I'])
        }
        fn target_oid(&self, target_block_id: BlockId) -> OId {
            target_block_id as OId
        }
        fn begin_output_batch(
            &mut self,
            query_begin: BlockId,
            blocked_processing: bool,
            query_separator: u8,
        ) {
            self.output_batches
                .push((query_begin, blocked_processing, query_separator));
        }
        fn push_output(&mut self, query: BlockId, buffer: Option<OutputBuffer>) {
            self.outputs.push((query, buffer));
        }
        fn end_output_batch(&mut self) {
            self.ended_output_batches += 1;
        }
        fn elapsed_microseconds(&mut self, _stage: &'static str) -> i64 {
            3
        }
        fn worker_count(&mut self, threads: usize, tasks: Range<i64>) {
            self.workers.push((threads, tasks));
        }
        fn heartbeat_begin(&mut self, _query_end: BlockId) {
            self.heartbeat.push("begin");
        }
        fn heartbeat_end(&mut self) {
            self.heartbeat.push("end");
        }
    }

    fn hit(query: u32, subject: u64) -> Hit {
        Hit::new(query, subject, 0)
    }

    #[test]
    fn partition_never_splits_source_query_and_fetch_fills_query_gaps() {
        let hits = vec![hit(0, 1), hit(1, 2), hit(2, 3), hit(6, 4)];
        let partition = make_partition(&hits, 2, 2);
        assert_eq!(partition, vec![0, 3, 4]);
        let iterator = HitIterator {
            partition: &partition,
            parts: 2,
            data: &hits,
            query_begin: 0,
            query_end: 5,
            query_contexts: 2,
            swipe_all: false,
            blastn: false,
        };
        assert_eq!(
            iterator.fetch(0),
            vec![
                Hits {
                    query: 0,
                    range: Some(0..2)
                },
                Hits {
                    query: 1,
                    range: Some(2..3)
                }
            ]
        );
        assert_eq!(
            iterator.fetch(1),
            vec![
                Hits {
                    query: 2,
                    range: None
                },
                Hits {
                    query: 3,
                    range: Some(3..4)
                },
                Hits {
                    query: 4,
                    range: None
                }
            ]
        );
        assert_eq!(
            HitIterator {
                swipe_all: true,
                ..iterator
            }
            .fetch(7),
            vec![Hits {
                query: 7,
                range: None
            }]
        );
    }

    #[test]
    fn legacy_pipeline_handles_unaligned_null_output_and_aligned_tracking() {
        let mut backend = Backend::default();
        let mut state = AlignState::new(3, 3);
        let mut stat = Statistics::default();
        let mut dp = DpStat::default();
        let empty = Hits {
            query: 1,
            range: None,
        };
        let cfg = AlignConfig {
            report_unaligned: true,
            ..AlignConfig::default()
        };
        assert_eq!(
            legacy_pipeline(
                &empty,
                &[],
                &cfg,
                &mut backend,
                &mut state,
                &mut stat,
                &mut dp
            )
            .unwrap(),
            Some(vec![1, b'U'])
        );

        let hits = vec![hit(0, 1)];
        let populated = Hits {
            query: 0,
            range: Some(0..1),
        };
        backend.legacy_aligned = true;
        let cfg = AlignConfig {
            track_aligned_queries: true,
            ..AlignConfig::default()
        };
        assert_eq!(
            legacy_pipeline(
                &populated,
                &hits,
                &cfg,
                &mut backend,
                &mut state,
                &mut stat,
                &mut dp
            )
            .unwrap(),
            Some(vec![0, b'L'])
        );
        assert!(state.query_aligned[0]);
        assert_eq!(state.iteration_query_aligned, 1);

        let null_cfg = AlignConfig {
            output_is_null: true,
            track_aligned_queries: true,
            ..AlignConfig::default()
        };
        assert_eq!(
            legacy_pipeline(
                &populated,
                &hits,
                &null_cfg,
                &mut backend,
                &mut state,
                &mut stat,
                &mut dp
            )
            .unwrap(),
            None
        );
    }

    #[test]
    fn align_worker_emits_missing_queries_and_forwards_global_ranking_and_tracking() {
        let hits = vec![hit(2, 1)];
        let partition = vec![0, 1];
        let iterator = HitIterator {
            partition: &partition,
            parts: 1,
            data: &hits,
            query_begin: 0,
            query_end: 3,
            query_contexts: 2,
            swipe_all: false,
            blastn: false,
        };
        let cfg = AlignConfig {
            query_contexts: 2,
            track_aligned_queries: true,
            track_aligned_targets: true,
            global_ranking_targets: 7,
            ..AlignConfig::default()
        };
        let mut backend = Backend::default();
        let mut state = AlignState::new(3, 3);
        align_worker(&iterator, &cfg, &mut backend, &mut state, 0).unwrap();
        assert_eq!(
            backend.outputs,
            vec![(0, None), (1, Some(vec![1, b'O'])), (2, None)]
        );
        assert_eq!(backend.extend_calls, vec![(1, false, false, false, 7)]);
        assert!(state.query_aligned[1]);
        assert!(state.aligned_targets[1]);
        assert_eq!(state.dp_stat.gross_cells, 2);
    }

    #[test]
    fn single_query_blastn_extends_without_seed_range_and_uses_intermediate_output() {
        let partition = vec![0];
        let iterator = HitIterator {
            partition: &partition,
            parts: 0,
            data: &[],
            query_begin: 0,
            query_end: 2,
            query_contexts: 1,
            swipe_all: false,
            blastn: true,
        };
        let cfg = AlignConfig {
            blastn: true,
            blocked_processing: true,
            ..AlignConfig::default()
        };
        let mut backend = Backend::default();
        let mut state = AlignState::new(2, 2);
        align_worker(&iterator, &cfg, &mut backend, &mut state, 1).unwrap();
        assert_eq!(backend.extend_calls, vec![(1, true, false, true, 0)]);
        assert_eq!(backend.outputs, vec![(1, Some(vec![1, b'I']))]);
    }

    #[test]
    fn align_queries_preserves_buffer_random_access_worker_heartbeat_and_statistics_lifecycle() {
        let mut backend = Backend {
            batches: VecDeque::from([HitBatch {
                hits: vec![hit(4, 2), hit(0, 1)],
                query_begin: 0,
                query_end: 3,
            }]),
            // The initial load result is intentionally ignored; the second
            // controls loop continuation, matching C++.
            load_results: VecDeque::from([true, false]),
            next_bin: 10,
            ..Backend::default()
        };
        let cfg = AlignConfig {
            min_task_trace_pts: 1,
            query_contexts: 2,
            threads: 4,
            threads_align: 2,
            verbosity: 3,
            heartbeat: true,
            memory_limit: 1000,
            trace_pt_fetch_size: 500,
            hit_record_size: 15,
            track_aligned_queries: true,
            track_aligned_targets: true,
            ..AlignConfig::default()
        };
        let mut state = AlignState::new(3, 3);
        align_queries(&cfg, &mut backend, &mut state).unwrap();
        assert!(backend.allocated && backend.freed);
        assert_eq!(backend.random_access, vec!["begin", "end"]);
        assert_eq!(backend.capacities, vec![500, 500]);
        assert_eq!(backend.workers, vec![(2, 0..1)]);
        assert_eq!(backend.output_batches, vec![(0, false, b'\n')]);
        assert_eq!(backend.ended_output_batches, 1);
        assert_eq!(backend.heartbeat, vec!["begin", "end"]);
        assert_eq!(
            backend.outputs,
            vec![
                (0, Some(vec![0, b'O'])),
                (1, None),
                (2, Some(vec![2, b'O']))
            ]
        );
        assert_eq!(state.statistics.get(StatValue::TimeLoadSeedHits), 3);
        assert_eq!(state.statistics.get(StatValue::TimeSortSeedHits), 3);
        assert_eq!(state.statistics.get(StatValue::TimeExt), 3);
        assert_eq!(state.statistics.get(StatValue::SearchTempSpace), 77);
    }
}
