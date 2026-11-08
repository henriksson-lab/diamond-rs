//! Actual grouped hit-buffer adapter for Source's alignment orchestration.
//! Owns the scoped query workers and live ordered sink; retained hit vectors
//! remain owned by one batch until every consumer finishes, then are recycled.
use std::io;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use super::align::{AlignBackend, AlignConfig, AlignState};
use crate::output::live_sink::LiveSink;
use crate::search::hit_buffer::CompactHit;

pub(crate) struct GroupedHitBatch {
    pub hits: Vec<Vec<CompactHit>>,
    /// Context indices, not source-query indices; six frames share an owner.
    pub query_begin: u32,
    pub query_end: u32,
}

/// Required concrete storage operations. The production implementation owns
/// the actual AdaptiveHitStore, including its asynchronous disk loader.
pub(crate) trait GroupedHitSource {
    fn prepare(&mut self) -> Result<(), String>;
    fn load(&mut self, max_bytes: u64) -> bool;
    fn retrieve(&mut self) -> Result<Option<GroupedHitBatch>, String>;
    fn next_bin_size(&self) -> u64;
    fn recycle(&mut self, batch: GroupedHitBatch);
    fn disk_size(&self) -> i64;
    fn finish(&mut self) -> Result<(), String>;
}

/// Borrow titles from the actual query owner; allocate only a current batch's
/// heartbeat references rather than duplicating all query-block metadata.
pub(crate) trait QueryTitle {
    fn query_title_bytes(&self) -> &[u8];
}
impl<T: AsRef<[u8]>> QueryTitle for crate::data::fasta::FastaRecord<T> {
    fn query_title_bytes(&self) -> &[u8] {
        self.id.as_ref()
    }
}

/// Completed query payloads include actual HSP/Match owners or final bytes,
/// rather than IDs into a side table. The concrete output controls buffering.
pub(crate) trait QueryOutput: Send {
    fn buffered_size(&self) -> usize;
}

/// Source make_partition over already-grouped source contexts. Does not copy
/// or flatten the retained hits just to recover their implicit query IDs.
pub(crate) fn grouped_partitions(
    hits: &[Vec<CompactHit>],
    count: usize,
    contexts: usize,
    minimum: usize,
) -> Vec<std::ops::Range<usize>> {
    assert!(contexts > 0 && minimum > 0);
    debug_assert!(hits.len() >= count.saturating_mul(contexts));
    let mut ranges = Vec::new();
    let mut begin = 0;
    let mut points = 0usize;
    for source in 0..count {
        points = points.saturating_add(
            hits[source * contexts..(source + 1) * contexts]
                .iter()
                .map(Vec::len)
                .sum::<usize>(),
        );
        if points > minimum {
            ranges.push(begin..source + 1);
            begin = source + 1;
            points = 0;
        }
    }
    if begin < count {
        if points == 0 && !ranges.is_empty() {
            ranges.last_mut().unwrap().end = count;
        } else {
            ranges.push(begin..count);
        }
    }
    ranges
}

/// Production backend. `process` is the actual DP/capture dependency and
/// `emit` writes the real completed payload. The adapter itself owns all
/// scheduling, partitioning, scratch lifetimes, heartbeat and sink draining.
pub(crate) struct NativeAlignBackend<'a, S, P, E, Scratch, O, T> {
    pub source: S,
    pub titles: &'a [T],
    pub contexts: usize,
    pub resident_bytes: u64,
    pub pool: Option<&'a rayon::ThreadPool>,
    pub serial: bool,
    pub heartbeat: bool,
    pub process: P,
    pub emit: E,
    load_elapsed: i64,
    extension_elapsed: i64,
    marker: PhantomData<fn() -> (Scratch, O)>,
}

impl<'a, S, P, E, Scratch, O, T> NativeAlignBackend<'a, S, P, E, Scratch, O, T> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        source: S,
        titles: &'a [T],
        contexts: usize,
        resident_bytes: u64,
        pool: Option<&'a rayon::ThreadPool>,
        serial: bool,
        heartbeat: bool,
        process: P,
        emit: E,
    ) -> Self {
        Self {
            source,
            titles,
            contexts,
            resident_bytes,
            pool,
            serial,
            heartbeat,
            process,
            emit,
            load_elapsed: 0,
            extension_elapsed: 0,
            marker: PhantomData,
        }
    }
}

impl<S, P, E, Scratch, O, T> AlignBackend for NativeAlignBackend<'_, S, P, E, Scratch, O, T>
where
    S: GroupedHitSource,
    T: QueryTitle + Sync,
    P: Fn(usize, &[Vec<CompactHit>], &mut Scratch) -> io::Result<O> + Sync,
    E: FnMut(usize, O) -> io::Result<()> + Send,
    Scratch: Default + Send,
    O: QueryOutput,
{
    type Batch = GroupedHitBatch;
    fn prepare_alignment(&mut self, _: &AlignConfig) -> Result<(), String> {
        if self.contexts == 0 {
            return Err("Zero alignment query contexts".into());
        }
        self.source.prepare()
    }
    fn resident_memory_size(&self) -> u64 {
        self.resident_bytes
    }
    fn load_batch(&mut self, max_bytes: u64) -> bool {
        self.source.load(max_bytes)
    }
    fn retrieve_batch(&mut self) -> Result<Option<GroupedHitBatch>, String> {
        let started = Instant::now();
        let batch = self.source.retrieve();
        self.load_elapsed = started.elapsed().as_micros() as i64;
        batch
    }
    fn next_batch_size(&self) -> u64 {
        self.source.next_bin_size()
    }
    fn batch_hit_count(&self, batch: &GroupedHitBatch) -> u64 {
        batch.hits.iter().map(|hits| hits.len() as u64).sum()
    }
    fn execute_batch(
        &mut self,
        batch: &mut GroupedHitBatch,
        cfg: &AlignConfig,
        _: &mut AlignState,
    ) -> Result<(), String> {
        let started = Instant::now();
        self.extension_elapsed = 0;
        if batch.query_begin as usize % self.contexts != 0
            || batch.query_end as usize % self.contexts != 0
        {
            return Err("Retained-hit batch splits a source query".into());
        }
        let begin = batch.query_begin as usize / self.contexts;
        let end = (batch.query_end as usize / self.contexts).min(self.titles.len());
        if end < begin || batch.hits.len() < (end - begin) * self.contexts {
            return Err("Invalid grouped alignment query range".into());
        }
        let hits = &batch.hits[..(end - begin) * self.contexts];
        if !cfg.swipe_all && hits.iter().all(Vec::is_empty) {
            return Ok(());
        }
        let count = end - begin;
        let partitions = if cfg.swipe_all {
            (0..count).map(|id| id..id + 1).collect()
        } else {
            grouped_partitions(hits, count, self.contexts, cfg.min_task_trace_pts)
        };
        let titles: Vec<_> = self.titles[begin..end]
            .iter()
            .map(QueryTitle::query_title_bytes)
            .collect();
        let pool = self.pool;
        let process = &self.process;
        let emit = &mut self.emit;
        let sink = LiveSink::new(begin, self.heartbeat, |(query, output): (usize, O)| {
            emit(query, output)
        });
        let contexts = self.contexts;
        let serial = self.serial;
        let worker_count = pool.map_or_else(
            rayon::current_num_threads,
            rayon::ThreadPool::current_num_threads,
        );
        let execute = || {
            sink.execute_bytes(begin, &titles, || {
                let process_partition =
                    |partition: usize, scratch: &mut Scratch| -> io::Result<()> {
                        let _progress = sink.task_scope();
                        for offset in partitions[partition].clone() {
                            let output = process(
                                begin + offset,
                                &hits[offset * contexts..(offset + 1) * contexts],
                                scratch,
                            )?;
                            let bytes = output.buffered_size();
                            sink.push(begin + offset, (begin + offset, output), bytes)?;
                        }
                        Ok(())
                    };
                if worker_count == 1 || serial || partitions.len() <= 1 {
                    let mut scratch = Scratch::default();
                    for index in 0..partitions.len() {
                        process_partition(index, &mut scratch)?;
                    }
                    return Ok(());
                }
                let next = AtomicUsize::new(0);
                let failed = AtomicBool::new(false);
                let error = Mutex::new(None);
                rayon::scope(|scope| {
                    for _ in 0..worker_count.min(partitions.len()) {
                        let next = &next;
                        let failed = &failed;
                        let error = &error;
                        let process_partition = &process_partition;
                        let partition_count = partitions.len();
                        scope.spawn(move |_| {
                            let mut scratch = Scratch::default();
                            while !failed.load(Ordering::Acquire) {
                                let partition = next.fetch_add(1, Ordering::Relaxed);
                                if partition >= partition_count {
                                    break;
                                }
                                if let Err(failure) = process_partition(partition, &mut scratch) {
                                    *error.lock().unwrap() = Some(failure);
                                    failed.store(true, Ordering::Release);
                                    break;
                                }
                            }
                        });
                    }
                });
                match error.into_inner().unwrap() {
                    Some(error) => Err(error),
                    None => Ok(()),
                }
            })
        };
        if let Some(pool) = pool {
            pool.install(execute)
        } else {
            execute()
        }
        .and_then(|()| sink.finish(end))
        .map_err(|error| error.to_string())?;
        self.extension_elapsed = started.elapsed().as_micros() as i64;
        Ok(())
    }
    fn recycle_batch(&mut self, batch: GroupedHitBatch) {
        self.source.recycle(batch);
    }
    fn disk_size(&self) -> i64 {
        self.source.disk_size()
    }
    fn complete_alignment(&mut self, _: &AlignConfig) -> Result<(), String> {
        self.source.finish()
    }
    fn elapsed(&mut self, stage: &'static str) -> i64 {
        match stage {
            "Loading trace points" => self.load_elapsed,
            "Computing alignments" => self.extension_elapsed,
            _ => 0,
        }
    }
    fn log_message(&mut self, message: &str) {
        let _ = crate::util::log_stream::log_stream()
            .lock()
            .unwrap()
            .write(format!("{message}\n"));
    }
}
