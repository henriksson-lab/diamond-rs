//! External-clustering alignment stage from `cluster/external/align.cpp`.

use crate::util::io::source_c_path;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use super::external::{ClusterStats, Edge, Job, PairEntryShort};
use super::pair_table::{RadixedTable, RADIX_BITS, RADIX_COUNT};
use super::seed_table::{SequenceVolumeReader, Volume, VolumedFile};
use super::storage::{BucketWriter, WorkerBucketBuffers};
use crate::basic::statistics::Statistics;
use crate::basic::value::Letter;
use crate::dp::swipe::{
    bin as swipe_bin, swipe, targets, CarryOver, DpTarget, Flags, HspValues, Params, SwipeTaskPool,
};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::hash::hash64;
use crate::util::parallel::{Atomic, TaskSet, ThreadPool};

#[derive(Debug, Clone, PartialEq)]
pub struct AlignConfig {
    pub threads: usize,
    pub approx_min_identity: f64,
    pub member_cover: f64,
    pub mutual_cover: Option<f64>,
    pub alignment_output: Option<PathBuf>,
    pub cutoff_score_8bit: i32,
    pub max_swipe_dp: i64,
    pub swipe_task_size: i64,
    pub band_bin: i32,
    pub col_bin: i32,
    pub approx_backtrace: bool,
    pub max_evalue: f64,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            approx_min_identity: 0.0,
            member_cover: 80.0,
            mutual_cover: None,
            alignment_output: None,
            cutoff_score_8bit: i8::MAX as i32,
            max_swipe_dp: 1_000_000,
            swipe_task_size: 100_000_000,
            band_bin: 24,
            col_bin: 400,
            approx_backtrace: false,
            max_evalue: f64::MAX,
        }
    }
}

#[derive(Debug)]
struct ChunkRecord {
    // Retain the actual loaded title owner along with its sequence payload.
    _title: Vec<u8>,
    sequence: Arc<[Letter]>,
}
#[derive(Debug)]
struct ChunkVolume {
    records: Vec<ChunkRecord>,
    oid_to_record: HashMap<i64, usize>,
    min_oid: i64,
    max_oid: i64,
}
#[derive(Debug)]
struct ChunkStorage {
    manifest: PathBuf,
    volumes: Vec<PathBuf>,
}
#[derive(Debug)]
struct LoadedChunkSequences {
    blocks: Vec<ChunkVolume>,
    threads: usize,
    cleanup: Option<ChunkStorage>,
    source_destructor: bool,
}
impl Drop for LoadedChunkSequences {
    fn drop(&mut self) {
        let workers = self.threads.min(self.blocks.len());
        let blocks = Mutex::new(std::mem::take(&mut self.blocks).into_iter());
        // Source deletes each real Block and its OID map in bounded workers,
        // then removes the VolumedFile. No global payload map survives this.
        std::thread::scope(|scope| {
            for _ in 0..workers {
                let blocks = &blocks;
                let spawn = std::thread::Builder::new().spawn_scoped(scope, move || loop {
                    let Some(block) = blocks.lock().unwrap().next() else {
                        break;
                    };
                    let ChunkVolume {
                        records,
                        oid_to_record,
                        ..
                    } = block;
                    drop(records);
                    drop(oid_to_record);
                });
                if spawn.is_err() && self.source_destructor {
                    // A std::thread construction failure in Source's noexcept
                    // destructor invokes its abort-only terminate handler.
                    crate::util::io::abort_source_process();
                }
            }
        });
        // A reusable library destructor cannot throw; finish any remaining
        // payloads locally if its deletion-worker creation failed.
        drop(blocks);
        if let Some(storage) = &self.cleanup {
            for path in &storage.volumes {
                let _ = fs::remove_file(source_c_path(path).as_ref());
            }
            let _ = fs::remove_file(source_c_path(&storage.manifest).as_ref());
            if let Some(directory) = storage.manifest.parent() {
                let _ = fs::remove_dir(source_c_path(directory).as_ref());
            }
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ChunkSequences {
    // Explicit in-memory library fixtures retain their former convenience
    // representation. Canonical file owners leave this map empty.
    sequences: HashMap<u64, Arc<[Letter]>>,
    loaded: Option<Arc<LoadedChunkSequences>>,
    letters: usize,
    oids: usize,
    volumes: usize,
}

impl ChunkSequences {
    pub fn from_volumes<R: SequenceVolumeReader>(
        volumes: &VolumedFile,
        reader: &mut R,
    ) -> Result<Self, String> {
        Self::load(volumes, reader, 1, None, Default::default())
    }

    fn from_manifest<R: SequenceVolumeReader>(
        path: &Path,
        reader: &mut R,
        threads: usize,
        errors: crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        let volumes = read_volumed_file(path, &errors)?;
        let cleanup = ChunkStorage {
            manifest: path.to_path_buf(),
            volumes: volumes.0.iter().map(|volume| volume.path.clone()).collect(),
        };
        Self::load(&volumes, reader, threads, Some(cleanup), errors)
    }

    fn load<R: SequenceVolumeReader>(
        volumes: &VolumedFile,
        reader: &mut R,
        threads: usize,
        cleanup: Option<ChunkStorage>,
        errors: crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        let mut blocks = Vec::with_capacity(volumes.0.len());
        if !volumes.0.is_empty() {
            if let Some(first_reader) = reader.fork_reader() {
                let slots: Arc<Vec<Mutex<Option<ChunkVolume>>>> =
                    Arc::new((0..volumes.0.len()).map(|_| Mutex::new(None)).collect());
                let paths: Arc<Vec<PathBuf>> =
                    Arc::new(volumes.0.iter().map(|volume| volume.path.clone()).collect());
                let next = Arc::new(std::sync::atomic::AtomicUsize::new(0));
                let mut first_reader = Some(first_reader);
                let mut pool = crate::util::parallel::SimpleThreadPool::new();
                for _ in 0..threads.min(paths.len()) {
                    let mut input = first_reader
                        .take()
                        .or_else(|| reader.fork_reader())
                        .ok_or("Chunk reader factory became unavailable.")?;
                    let slots = Arc::clone(&slots);
                    let paths = Arc::clone(&paths);
                    let next = Arc::clone(&next);
                    let errors = errors.clone();
                    pool.spawn(move |stop| {
                        while !stop.load(std::sync::atomic::Ordering::Relaxed) {
                            let index = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            let Some(path) = paths.get(index) else {
                                break;
                            };
                            match load_chunk_volume(path, input.as_mut()) {
                                Ok(block) => *slots[index].lock().unwrap() = Some(block),
                                Err(message) => {
                                    if errors.source_destructors() {
                                        // Source ChunkSeqs uses plain std::thread:
                                        // an input/Block exception terminates the
                                        // process, rather than reaching main.
                                        let error = errors.finish::<()>(Err(message)).unwrap_err();
                                        if let Err(error) =
                                            crate::util::io::emit_source_perror(&error)
                                        {
                                            crate::util::io::abort_source_diagnostic_error(&error);
                                        }
                                        crate::util::io::abort_source_process();
                                    }
                                    return Err(message);
                                }
                            }
                        }
                        Ok(())
                    });
                }
                pool.join_all()?;
                for slot in slots.iter() {
                    blocks.push(
                        slot.lock()
                            .unwrap()
                            .take()
                            .ok_or("Chunk volume was not loaded.")?,
                    );
                }
            } else {
                // Explicit legacy adapters without independent cursor factories
                // remain fallible and serial; actual CLI readers fork above.
                for volume in &volumes.0 {
                    blocks.push(load_chunk_volume(&volume.path, reader)?);
                }
            }
        }
        let oids = blocks.iter().map(|block| block.records.len()).sum();
        let letters = blocks
            .iter()
            .flat_map(|block| &block.records)
            .map(|record| record.sequence.len())
            .sum();
        // Attach filesystem destruction ONLY after every actual volume has
        // loaded successfully. A failing Source constructor never removes it.
        Ok(Self {
            sequences: HashMap::new(),
            loaded: Some(Arc::new(LoadedChunkSequences {
                blocks,
                threads,
                cleanup,
                source_destructor: errors.source_destructors(),
            })),
            letters,
            oids,
            volumes: volumes.0.len(),
        })
    }

    fn record(&self, oid: u64) -> Result<&Arc<[Letter]>, String> {
        if let Some(loaded) = &self.loaded {
            let signed = oid as i64;
            for block in &loaded.blocks {
                if signed < block.min_oid || signed > block.max_oid {
                    continue;
                }
                if let Some(index) = block.oid_to_record.get(&signed) {
                    return Ok(&block.records[*index].sequence);
                }
            }
            return Err("ChunkSeqs".into());
        }
        self.sequences
            .get(&oid)
            .ok_or_else(|| "ChunkSeqs".to_owned())
    }
    pub fn get(&self, oid: u64) -> Result<&[Letter], String> {
        self.record(oid).map(AsRef::as_ref)
    }
    fn get_shared(&self, oid: u64) -> Result<Arc<[Letter]>, String> {
        self.record(oid).cloned()
    }
    pub fn oids(&self) -> usize {
        if self.loaded.is_some() {
            self.oids
        } else {
            self.sequences.len()
        }
    }
    pub fn letters(&self) -> usize {
        self.letters
    }
    pub fn volumes(&self) -> usize {
        self.volumes
    }
}

fn load_chunk_volume(
    path: &Path,
    reader: &mut dyn SequenceVolumeReader,
) -> Result<ChunkVolume, String> {
    let mut records = Vec::new();
    let mut oid_to_record = HashMap::new();
    let mut min_oid = i64::MAX;
    let mut max_oid = i64::MIN;
    reader.visit_volume_bytes(path, &mut |record| {
        let oid = parse_oid(&record.id) as i64;
        min_oid = min_oid.min(oid);
        max_oid = max_oid.max(oid);
        oid_to_record.insert(oid, records.len());
        records.push(ChunkRecord {
            _title: record.id,
            sequence: record.sequence.into(),
        });
        Ok(())
    })?;
    Ok(ChunkVolume {
        records,
        oid_to_record,
        min_oid,
        max_oid,
    })
}

fn parse_oid(id: &[u8]) -> u64 {
    // ChunkSeqs stores std::atoll(ids[j]); it accepts signs, numeric prefixes,
    // empty/nonnumeric strings and arbitrary bytes after the numeric prefix.
    let mut terminated = crate::util::sequence::c_string_bytes(id).to_vec();
    terminated.push(0);
    unsafe { libc::atoll(terminated.as_ptr().cast()) as u64 }
}

fn encode_edge(edge: Edge) -> [u8; 24] {
    let mut bytes = [0u8; 24];
    bytes[0..8].copy_from_slice(&edge.rep_oid.to_ne_bytes());
    bytes[8..16].copy_from_slice(&edge.member_oid.to_ne_bytes());
    bytes[16..20].copy_from_slice(&edge.rep_len.to_ne_bytes());
    bytes[20..24].copy_from_slice(&edge.member_len.to_ne_bytes());
    bytes
}

/// C++ file-local `align_rep` using the shared full-matrix SWIPE engine.
pub fn align_rep(
    chunk_sequences: &ChunkSequences,
    pairs: &[PairEntryShort],
    score_matrix: &ScoreMatrix,
    config: &AlignConfig,
    output: &mut Vec<Edge>,
    alignment_output: Option<&mut dyn Write>,
    stats: &mut ClusterStats,
) -> Result<(), String> {
    align_rep_with_pool(
        chunk_sequences,
        pairs,
        score_matrix,
        config,
        output,
        alignment_output,
        stats,
        None,
    )
}

fn align_rep_with_pool(
    chunk_sequences: &ChunkSequences,
    pairs: &[PairEntryShort],
    score_matrix: &ScoreMatrix,
    config: &AlignConfig,
    output: &mut Vec<Edge>,
    alignment_output: Option<&mut dyn Write>,
    stats: &mut ClusterStats,
    pool: Option<&SwipeTaskPool>,
) -> Result<(), String> {
    let Some(first) = pairs.first() else {
        return Ok(());
    };
    if pairs.iter().any(|pair| pair.rep_oid != first.rep_oid) {
        return Err("Pair group contains multiple representative OIDs.".to_string());
    }
    let representative = chunk_sequences.get(first.rep_oid)?;
    let mut dp_targets = targets();
    for (index, pair) in pairs.iter().enumerate() {
        let member = chunk_sequences.get(pair.member_oid)?;
        let bin = swipe_bin(
            HspValues::COORDS,
            representative.len() as i32,
            0,
            0,
            representative.len() as i64 * member.len() as i64,
            0,
            0,
            config.cutoff_score_8bit,
            config.max_swipe_dp,
            config.approx_backtrace,
        );
        dp_targets[bin].push_back(DpTarget::full(
            chunk_sequences.get_shared(pair.member_oid)?,
            member.len() as i32,
            index as i64,
            CarryOver::default(),
        ));
        stats.extensions_computed += 1;
    }
    let swipe_statistics = Arc::new(Mutex::new(Statistics::new()));
    let mut params = Params::new(representative, score_matrix);
    params.flags = Flags::FULL_MATRIX;
    params.v = HspValues::COORDS;
    params.cutoff_score_8bit = config.cutoff_score_8bit;
    params.max_swipe_dp = config.max_swipe_dp;
    params.swipe_task_size = config.swipe_task_size;
    params.task_pool = pool.cloned();
    params.band_bin = config.band_bin;
    params.col_bin = config.col_bin;
    params.approx_backtrace = config.approx_backtrace;
    params.max_evalue = config.max_evalue;
    // External::external sets these global controls before align_rep. SWIPE
    // consumes them while selecting forward endpoints for the reverse pass,
    // before the final external edge filters below.
    if let Some(mutual_cover) = config.mutual_cover {
        params.query_cover = mutual_cover;
        params.subject_cover = mutual_cover;
    } else {
        params.query_or_target_cover = config.member_cover;
    }
    params.approx_min_id = config.approx_min_identity;
    params.statistics = Some(swipe_statistics);
    let hsps = swipe(&dp_targets, &mut params);
    stats.hits_evalue_filtered += hsps.len() as u64;
    let mut alignment_output = alignment_output;
    for hsp in hsps {
        let pair = pairs
            .get(hsp.swipe_target as usize)
            .ok_or_else(|| "Invalid SWIPE target index.".to_string())?;
        let member = chunk_sequences.get(pair.member_oid)?;
        let query_cover = hsp.query_cover_percent(representative.len() as u32);
        let subject_cover = hsp.subject_cover_percent(member.len() as u32);
        if let Some(writer) = alignment_output.as_deref_mut() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                pair.rep_oid,
                pair.member_oid,
                crate::output::format::format_default_float(query_cover),
                crate::output::format::format_default_float(subject_cover),
                crate::output::format::format_default_float(hsp.evalue)
            )
            .map_err(|error| error.to_string())?;
        }
        if hsp.approx_id_percent(representative, member) < config.approx_min_identity {
            continue;
        }
        let mut passed = false;
        if let Some(mutual_cover) = config.mutual_cover {
            if subject_cover >= mutual_cover && query_cover >= mutual_cover {
                let (rep_oid, member_oid, rep_len, member_len) = if pair.rep_oid <= pair.member_oid
                {
                    (
                        pair.rep_oid,
                        pair.member_oid,
                        representative.len() as u32,
                        member.len() as u32,
                    )
                } else {
                    (
                        pair.member_oid,
                        pair.rep_oid,
                        member.len() as u32,
                        representative.len() as u32,
                    )
                };
                output.push(Edge {
                    rep_oid,
                    member_oid,
                    rep_len,
                    member_len,
                });
                passed = true;
            }
        } else {
            if subject_cover >= config.member_cover {
                output.push(Edge {
                    rep_oid: pair.rep_oid,
                    member_oid: pair.member_oid,
                    rep_len: representative.len() as u32,
                    member_len: member.len() as u32,
                });
                passed = true;
            }
            if query_cover >= config.member_cover {
                output.push(Edge {
                    rep_oid: pair.member_oid,
                    member_oid: pair.rep_oid,
                    rep_len: member.len() as u32,
                    member_len: representative.len() as u32,
                });
                passed = true;
            }
        }
        stats.hits_filtered += u64::from(passed);
    }
    Ok(())
}

fn read_volumed_file(
    path: &Path,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<VolumedFile, String> {
    super::external::read_volumed_file_manifest_owned(path, errors)
}

#[cfg(test)]
fn visit_pair_batches(
    path: &Path,
    visitor: &mut dyn FnMut(Vec<PairEntryShort>) -> Result<(), String>,
) -> Result<(), String> {
    let errors = crate::util::io::IoErrorOwner::default();
    let mut input = errors.capture(
        crate::util::io::InputFile::with_source_profile_path(
            path,
            crate::util::io::InputFile::NO_AUTODETECT,
            65_536,
        )
        .map_err(std::io::Error::other),
    )?;
    visit_pair_batches_input(&mut input, &errors, visitor)
}

fn visit_pair_batches_input(
    input: &mut crate::util::io::InputFile,
    errors: &crate::util::io::IoErrorOwner,
    visitor: &mut dyn FnMut(Vec<PairEntryShort>) -> Result<(), String>,
) -> Result<(), String> {
    let total = errors.capture(input.file_size().map_err(std::io::Error::other))? as u64;
    let mut position = 0;
    while let Some(batch) = next_pair_batch(input, errors, total, &mut position)? {
        visitor(batch)?;
    }
    errors.capture(input.close().map_err(std::io::Error::other))?;
    Ok(())
}

fn next_pair_batch(
    input: &mut crate::util::io::InputFile,
    errors: &crate::util::io::IoErrorOwner,
    total: u64,
    position: &mut u64,
) -> Result<Option<Vec<PairEntryShort>>, String> {
    let mut header = [0; std::mem::size_of::<usize>()];
    if errors.capture(
        input
            .read_raw(&mut header[..1])
            .map_err(std::io::Error::other),
    )? == 0
    {
        return Ok(None);
    }
    errors.capture(
        input
            .read_exact(&mut header[1..])
            .map_err(std::io::Error::other),
    )?;
    let count = usize::from_ne_bytes(header);
    let bytes = count.checked_mul(16).ok_or("Pair batch size overflow.")?;
    *position += header.len() as u64;
    if bytes as u64 > total.saturating_sub(*position) {
        return Err("Short PairEntryShort batch.".into());
    }
    let mut batch = Vec::new();
    batch
        .try_reserve(count)
        .map_err(|error| error.to_string())?;
    for _ in 0..count {
        let mut pair = [0; 16];
        errors.capture(input.read_values(&mut pair).map_err(std::io::Error::other))?;
        *position += pair.len() as u64;
        batch.push(PairEntryShort {
            rep_oid: u64::from_ne_bytes(pair[..8].try_into().unwrap()),
            member_oid: u64::from_ne_bytes(pair[8..].try_into().unwrap()),
        });
    }
    Ok(Some(batch))
}

#[cfg(test)]
fn read_pair_batches(path: &Path) -> Result<Vec<Vec<PairEntryShort>>, String> {
    let mut batches = Vec::new();
    visit_pair_batches(path, &mut |batch| {
        batches.push(batch);
        Ok(())
    })?;
    Ok(batches)
}

fn alignment_output_path(prefix: &Path) -> PathBuf {
    #[cfg(unix)]
    let id = unsafe { libc::pthread_self() as usize };
    #[cfg(windows)]
    let id = unsafe { windows_sys::Win32::System::Threading::GetCurrentThreadId() as usize };
    #[cfg(not(any(unix, windows)))]
    let id = std::process::id() as usize;
    let mut path = prefix.as_os_str().to_owned();
    path.push(format!(".{id}"));
    PathBuf::from(path)
}

/// Compatibility library entry; the native command supplies its existing Arc
/// through align_owned, so no scoring matrix is copied by the production path.
pub fn align<R: SequenceVolumeReader>(
    job: &Job,
    chunk_count: usize,
    db_size: i64,
    score_matrix: &ScoreMatrix,
    config: &AlignConfig,
    reader: &mut R,
) -> Result<(RadixedTable, ClusterStats), String> {
    align_owned(
        job,
        chunk_count,
        db_size,
        Arc::new(score_matrix.clone()),
        config,
        reader,
    )
}

struct AlignmentOutput {
    file: crate::util::io::cstdio_file::OwnedCstdio,
    failed: bool,
}
impl Write for AlignmentOutput {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if !self.failed {
            self.failed = !self
                .file
                .write_text_output(bytes)
                .map_err(std::io::Error::other)?;
        }
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        if !self.failed {
            self.failed = !self
                .file
                .flush_text_output()
                .map_err(std::io::Error::other)?;
        }
        Ok(())
    }
    fn write_fmt(&mut self, arguments: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        let mut row = Vec::new();
        row.write_fmt(arguments)?;
        self.write_all(&row)?;
        self.flush() // The sole caller is Source's per-HSP std::endl row.
    }
}
impl Drop for AlignmentOutput {
    fn drop(&mut self) {
        let _ = self.file.close_file_role();
    }
}

struct PairInput {
    input: crate::util::io::InputFile,
    total: u64,
    position: u64,
}
struct AlignChunkOwner {
    sequences: ChunkSequences,
    input: Mutex<PairInput>,
    files: Arc<Mutex<BucketWriter>>,
    matrix: Arc<ScoreMatrix>,
    pool: ThreadPool,
    config: AlignConfig,
    errors: crate::util::io::IoErrorOwner,
    failure: Arc<Mutex<Option<String>>>,
    statistics: Mutex<ClusterStats>,
}
thread_local! {
    // ThreadPool workers are freshly created for each chunk. Their actual TLS
    // destructor flushes BufferArray before ThreadPool::join returns.
    static ALIGN_BUFFERS: RefCell<Option<WorkerBucketBuffers>> = const { RefCell::new(None) };
}
impl AlignChunkOwner {
    fn enqueue(self: &Arc<Self>, tasks: &Arc<TaskSet>) {
        let owner = Arc::clone(self);
        let tasks_for_task = Arc::clone(tasks);
        tasks.enqueue(move || {
            if owner.failure.lock().unwrap().is_some() {
                return;
            }
            if let Err(message) = owner.run_task(&tasks_for_task) {
                if owner.errors.source_destructors() {
                    let error = owner.errors.finish::<()>(Err(message)).unwrap_err();
                    if let Err(error) = crate::util::io::emit_source_perror(&error) {
                        crate::util::io::abort_source_diagnostic_error(&error);
                    }
                    crate::util::io::abort_source_process();
                }
                owner.failure.lock().unwrap().get_or_insert(message);
            }
        });
    }
    fn run_task(self: &Arc<Self>, tasks: &Arc<TaskSet>) -> Result<(), String> {
        ALIGN_BUFFERS.with(|slot| {
            if slot.borrow().is_none() {
                *slot.borrow_mut() = Some(WorkerBucketBuffers::new(
                    Arc::clone(&self.files),
                    RADIX_COUNT,
                    self.errors.clone(),
                    Arc::clone(&self.failure),
                ));
            }
        });
        // Source creates this per-task ofstream before reading a header, even
        // for the single EOF task. No idle worker invents an empty output file.
        let mut alignment_output = if let Some(prefix) = &self.config.alignment_output {
            let path = alignment_output_path(prefix);
            Some(AlignmentOutput {
                file: crate::util::io::cstdio_file::OwnedCstdio::open_file_role_path(&path, "a")
                    .map_err(|_| {
                        super::storage::path_error(
                            &self.errors,
                            b"Failed to open alignment output file: ",
                            &path,
                        )
                    })?,
                failed: false,
            })
        } else {
            None
        };
        let batch = {
            let mut input = self.input.lock().unwrap();
            let PairInput {
                input,
                total,
                position,
            } = &mut *input;
            next_pair_batch(input, &self.errors, *total, position)?
        };
        let Some(batch) = batch else {
            return Ok(());
        };
        // Exactly one successor is pending, in Source FIFO header/payload
        // order; it becomes runnable before this task starts expensive DP.
        self.enqueue(tasks);
        let mut stats = ClusterStats::default();
        let mut begin = 0;
        while begin < batch.len() {
            let mut end = begin + 1;
            while end < batch.len() && batch[end].rep_oid == batch[begin].rep_oid {
                end += 1;
            }
            let mut edges = Vec::new();
            let pool = SwipeTaskPool::new(
                self.pool.clone(),
                self.sequences.get_shared(batch[begin].rep_oid)?,
                Arc::clone(&self.matrix),
            );
            align_rep_with_pool(
                &self.sequences,
                &batch[begin..end],
                &self.matrix,
                &self.config,
                &mut edges,
                alignment_output.as_mut().map(|file| file as &mut dyn Write),
                &mut stats,
                Some(&pool),
            )?;
            ALIGN_BUFFERS.with(|slot| {
                let mut slot = slot.borrow_mut();
                let buffers = slot.as_mut().unwrap();
                for edge in edges {
                    let radix = (hash64(edge.member_oid) >> (64 - RADIX_BITS)) as usize;
                    buffers.write(radix, &encode_edge(edge))?;
                }
                Ok::<_, String>(())
            })?;
            begin = end;
        }
        self.statistics.lock().unwrap().add(&stats);
        Ok(())
    }
}

/// Actual native command owner, using the translated priority TaskSet and
/// sharing the already initialized scoring matrix across all worker tasks.
pub fn align_owned<R: SequenceVolumeReader>(
    job: &Job,
    chunk_count: usize,
    _db_size: i64,
    score_matrix: Arc<ScoreMatrix>,
    config: &AlignConfig,
    reader: &mut R,
) -> Result<(RadixedTable, ClusterStats), String> {
    let threads = i32::try_from(config.threads)
        .map_err(|_| "Alignment thread count exceeds signed32 range.")?;
    if threads <= 0 {
        return Err("Alignment thread count must be positive.".into());
    }
    let chunks_path = job.base_dir(None).join("chunks");
    let alignment_path = job.base_dir(None).join("alignments");
    job.capture_io(crate::util::system::mkdir_path(&alignment_path))?;
    let output = Arc::new(Mutex::new(BucketWriter::new_file_array_owned(
        &alignment_path,
        job.worker_id(),
        RADIX_COUNT,
        job.io_errors(),
    )?));
    let mut all_stats = ClusterStats::default();
    let queue = Atomic::new(chunks_path.join("align_queue"));
    let mut processed = 0;
    loop {
        let chunk = usize::try_from(queue.fetch_add_one()?).map_err(|error| error.to_string())?;
        if chunk >= chunk_count {
            break;
        }
        let chunk_path = chunks_path.join(chunk.to_string());
        let sequences = ChunkSequences::from_manifest(
            &chunk_path.join("bucket.tsv"),
            reader,
            config.threads,
            job.io_errors(),
        )?;
        job.log(&format!(
            "Computing alignments. Chunk={}/{} Volumes={} Sequences={} Letters={}",
            chunk + 1,
            chunk_count,
            sequences.volumes(),
            sequences.oids(),
            sequences.letters()
        ))?;
        let pairs_path = chunk_path.join("pairs");
        let mut input = job.open_input(&pairs_path, crate::util::io::InputFile::NO_AUTODETECT)?;
        let total = job.io(input.file_size())? as u64;
        let failure = Arc::new(Mutex::new(None));
        let pool = ThreadPool::default();
        let owner = Arc::new(AlignChunkOwner {
            sequences,
            input: Mutex::new(PairInput {
                input,
                total,
                position: 0,
            }),
            files: Arc::clone(&output),
            matrix: Arc::clone(&score_matrix),
            pool: pool.clone(),
            config: config.clone(),
            errors: job.io_errors(),
            failure,
            statistics: Mutex::new(ClusterStats::default()),
        });
        let tasks = Arc::new(TaskSet::new(&pool, 1));
        owner.enqueue(&tasks);
        pool.run(threads, true, Some(&tasks));
        pool.join();
        if let Some(error) = owner.failure.lock().unwrap().clone() {
            return Err(error);
        }
        all_stats.add(&owner.statistics.lock().unwrap());
        job.io(owner.input.lock().unwrap().input.close())?;
        fs::remove_file(source_c_path(&pairs_path).as_ref()).map_err(|error| error.to_string())?;
        drop(owner); // Actual sequence owner is freed after close/remove.
        processed += 1;
    }
    job.log_stats(&all_stats)?;
    let table = output.lock().unwrap().finish()?;
    drop(output);
    let finished = Atomic::new(alignment_path.join("finished"));
    finished.fetch_add(processed)?;
    finished.await_value(chunk_count as i64)?;
    Ok((table, all_stats))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::AMINO_ACID_ALPHABET;
    use crate::cluster::external::seed_table::SeedSequence;
    use crate::stats::score_matrix::ScoreMatrix;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    fn sequence(text: &str) -> Vec<Letter> {
        text.bytes()
            .map(|byte| {
                AMINO_ACID_ALPHABET
                    .iter()
                    .position(|&value| value == byte)
                    .unwrap() as Letter
            })
            .collect()
    }

    #[test]
    fn align_rep_emits_both_unidirectional_edges_and_exact_stats() {
        let mut sequences = ChunkSequences::default();
        sequences.sequences.insert(4, sequence("ARNDCQEGH").into());
        sequences.sequences.insert(9, sequence("ARNDCQEGH").into());
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap();
        let config = AlignConfig {
            member_cover: 80.0,
            ..AlignConfig::default()
        };
        let mut edges = Vec::new();
        let mut alignment = Vec::new();
        let mut stats = ClusterStats::default();
        align_rep(
            &sequences,
            &[PairEntryShort {
                rep_oid: 4,
                member_oid: 9,
            }],
            &matrix,
            &config,
            &mut edges,
            Some(&mut alignment),
            &mut stats,
        )
        .unwrap();
        assert_eq!(stats.extensions_computed, 1);
        assert_eq!(stats.hits_evalue_filtered, 1);
        assert_eq!(stats.hits_filtered, 1);
        assert_eq!(edges.len(), 2);
        assert_eq!((edges[0].rep_oid, edges[0].member_oid), (4, 9));
        assert_eq!((edges[1].rep_oid, edges[1].member_oid), (9, 4));
        let line = String::from_utf8(alignment).unwrap();
        let fields = line.trim_end().split('\t').collect::<Vec<_>>();
        assert_eq!(&fields[..2], &["4", "9"]);
        assert_eq!(fields.len(), 5);
        assert!(fields[2].parse::<f64>().unwrap() >= 80.0);
        assert!(fields[3].parse::<f64>().unwrap() >= 80.0);
    }

    #[test]
    fn align_rep_mutual_cover_canonicalizes_oid_order() {
        let mut sequences = ChunkSequences::default();
        sequences.sequences.insert(8, sequence("ARNDCQEGH").into());
        sequences.sequences.insert(2, sequence("ARNDCQEGH").into());
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap();
        let config = AlignConfig {
            mutual_cover: Some(80.0),
            ..AlignConfig::default()
        };
        let mut edges = Vec::new();
        let mut stats = ClusterStats::default();
        align_rep(
            &sequences,
            &[PairEntryShort {
                rep_oid: 8,
                member_oid: 2,
            }],
            &matrix,
            &config,
            &mut edges,
            None,
            &mut stats,
        )
        .unwrap();
        assert_eq!(edges.len(), 1);
        assert_eq!((edges[0].rep_oid, edges[0].member_oid), (2, 8));
        assert_eq!(stats.hits_filtered, 1);
    }

    #[test]
    fn align_rep_applies_mutual_coverage_before_reverse_alignment() {
        let mut sequences = ChunkSequences::default();
        // An unambiguous prefix match: the subject is fully covered while
        // only half of the representative is aligned. Repeated residues
        // do not distinguish equally scoring reverse endpoints.
        sequences
            .sequences
            .insert(0, sequence("ARNDCQEGHILKMFPSTWYV").into());
        sequences.sequences.insert(1, sequence("ARNDCQEGHI").into());
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap();
        let pairs = [PairEntryShort {
            rep_oid: 0,
            member_oid: 1,
        }];
        for (mutual_cover, expected_hsps) in [(None, 1), (Some(100.0), 0)] {
            let mut stats = ClusterStats::default();
            let mut edges = Vec::new();
            align_rep(
                &sequences,
                &pairs,
                &matrix,
                &AlignConfig {
                    member_cover: 100.0,
                    mutual_cover,
                    ..Default::default()
                },
                &mut edges,
                None,
                &mut stats,
            )
            .unwrap();
            assert_eq!(stats.extensions_computed, 1);
            assert_eq!(stats.hits_evalue_filtered, expected_hsps);
            assert_eq!(edges.len(), expected_hsps as usize);
        }
    }

    #[test]
    fn align_rep_applies_score_identity_bound_before_identity_shortcut() {
        let mut sequences = ChunkSequences::default();
        sequences
            .sequences
            .insert(0, sequence("AAAAAAAAAAAAAAAAAAAA").into());
        sequences
            .sequences
            .insert(1, sequence("AAAAAAAAAAAAAAAAAAAA").into());
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap();
        let mut stats = ClusterStats::default();
        let mut edges = Vec::new();
        align_rep(
            &sequences,
            &[PairEntryShort {
                rep_oid: 0,
                member_oid: 1,
            }],
            &matrix,
            &AlignConfig {
                member_cover: 100.0,
                approx_min_identity: 90.0,
                ..Default::default()
            },
            &mut edges,
            None,
            &mut stats,
        )
        .unwrap();
        // The pre-reverse score bound is 80/20*16.56+11.41 = 77.65.
        // A final identical-sequence check would instead report 100 percent.
        assert_eq!(stats.extensions_computed, 1);
        assert_eq!(stats.hits_evalue_filtered, 0);
        assert!(edges.is_empty());
    }

    #[test]
    fn pair_batches_reject_truncated_native_records() {
        let path = std::env::temp_dir().join(format!(
            "diamond-align-pairs-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ));
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&1usize.to_ne_bytes());
        bytes.extend_from_slice(&3u64.to_ne_bytes());
        fs::write(&path, bytes).unwrap();
        assert_eq!(
            read_pair_batches(&path).unwrap_err(),
            "Short PairEntryShort batch."
        );
        fs::remove_file(path).unwrap();
    }

    #[derive(Default)]
    struct FixtureReader;

    impl SequenceVolumeReader for FixtureReader {
        fn read_volume(&mut self, _path: &Path) -> Result<Vec<SeedSequence>, String> {
            Ok(vec![
                SeedSequence {
                    id: "4 representative".to_string(),
                    sequence: sequence("ARNDCQEGH"),
                },
                SeedSequence {
                    id: "9 member".to_string(),
                    sequence: sequence("ARNDCQEGH"),
                },
            ])
        }
    }

    #[derive(Clone, Default)]
    struct RealVolumeReader {
        errors: crate::util::io::IoErrorOwner,
    }
    impl SequenceVolumeReader for RealVolumeReader {
        fn fork_reader(&self) -> Option<Box<dyn SequenceVolumeReader + Send>> {
            Some(Box::new(self.clone()))
        }
        fn read_volume(&mut self, path: &Path) -> Result<Vec<SeedSequence>, String> {
            let mut records = Vec::new();
            self.visit_volume_bytes(path, &mut |record| {
                records.push(SeedSequence {
                    id: String::from_utf8(record.id).map_err(|e| e.to_string())?,
                    sequence: record.sequence,
                });
                Ok(())
            })?;
            Ok(records)
        }
        fn visit_volume_bytes(
            &mut self,
            path: &Path,
            visit: &mut dyn FnMut(super::super::seed_table::ByteSeedSequence) -> Result<(), String>,
        ) -> Result<(), String> {
            super::super::reader::visit_volume_bytes_owned(
                path,
                visit,
                &mut || Ok(()),
                &self.errors,
            )
        }
    }

    #[test]
    fn constructed_chunk_owner_cleans_only_after_last_clone_and_payload_release() {
        let root = std::env::temp_dir().join(format!(
            "diamond-chunk-owner-{}-{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&root).unwrap();
        let first = root.join("0.faa");
        let second = root.join("1.faa");
        fs::write(&first, b">4 first\nARNDCQEGH\n").unwrap();
        fs::write(&second, b">9 second\nARNDCQEGH\n").unwrap();
        let manifest = root.join("bucket.tsv");
        let mut rows = super::super::storage::manifest_line(&first, &[1]).unwrap();
        rows.extend(super::super::storage::manifest_line(&second, &[1]).unwrap());
        fs::write(&manifest, rows).unwrap();
        let errors = crate::util::io::IoErrorOwner::default();
        let mut reader = RealVolumeReader {
            errors: errors.clone(),
        };
        let owner = ChunkSequences::from_manifest(&manifest, &mut reader, 4, errors).unwrap();
        assert!(
            owner.sequences.is_empty(),
            "no global map retains canonical payloads"
        );
        assert_eq!(owner.loaded.as_ref().unwrap().blocks.len(), 2);
        assert_eq!(owner.oids(), 2);
        let weak = Arc::downgrade(owner.record(4).unwrap());
        let last = owner.clone();
        drop(owner);
        assert!(first.exists() && second.exists() && manifest.exists());
        assert!(weak.upgrade().is_some());
        drop(last);
        assert!(weak.upgrade().is_none());
        assert!(
            !root.exists(),
            "per-volume payload owners die before cache files and directory removal"
        );
    }

    #[test]
    fn failed_chunk_constructor_retains_manifest_and_original_volumes() {
        let root = std::env::temp_dir().join(format!(
            "diamond-chunk-ctor-error-{}-{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&root).unwrap();
        let first = root.join("0.faa");
        let missing = root.join("missing.faa");
        fs::write(&first, b">4 first\nARNDCQEGH\n").unwrap();
        let manifest = root.join("bucket.tsv");
        let mut rows = super::super::storage::manifest_line(&first, &[1]).unwrap();
        rows.extend(super::super::storage::manifest_line(&missing, &[1]).unwrap());
        fs::write(&manifest, rows).unwrap();
        let errors = crate::util::io::IoErrorOwner::default();
        let mut reader = RealVolumeReader {
            errors: errors.clone(),
        };
        let result = ChunkSequences::from_manifest(&manifest, &mut reader, 4, errors.clone());
        let error = errors.finish(result).unwrap_err();
        assert_eq!(
            crate::util::io::file_open_error::find(&error).unwrap().path,
            missing
        );
        assert!(first.exists() && manifest.exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn chunk_sequence_oid_uses_raw_c_numeric_prefix() {
        assert_eq!(parse_oid(b" +42\xff trailing"), 42);
        assert_eq!(parse_oid(b"-7 description"), (-7_i64) as u64);
        assert_eq!(parse_oid(b"9\0-1"), 9);
        assert_eq!(parse_oid(b"\xffnot a number"), 0);
        assert_eq!(parse_oid(b""), 0);
    }

    #[test]
    fn parallel_taskset_has_one_eof_task_and_no_idle_output_owners() {
        let root = std::env::temp_dir().join(format!(
            "diamond-align-empty-{}-{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        let job = Job::new(10, 1, 1 << 20, &root, 3).unwrap();
        let chunk = job.base_dir(None).join("chunks/0");
        fs::create_dir_all(&chunk).unwrap();
        fs::write(
            chunk.join("bucket.tsv"),
            super::super::storage::manifest_line(&chunk.join("fixture.fa"), &[2]).unwrap(),
        )
        .unwrap();
        fs::write(chunk.join("pairs"), []).unwrap();
        let matrix = Arc::new(ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap());
        let (_, stats) = align_owned(
            &job,
            1,
            10,
            matrix,
            &AlignConfig {
                threads: 4,
                alignment_output: Some(root.join("aln")),
                ..Default::default()
            },
            &mut FixtureReader,
        )
        .unwrap();
        let files: Vec<_> = fs::read_dir(&root)
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| p.file_name().unwrap().to_str().unwrap().starts_with("aln."))
            .collect();
        assert_eq!(
            files.len(),
            1,
            "only the actual EOF task opens its thread output"
        );
        assert_eq!(fs::metadata(&files[0]).unwrap().len(), 0);
        assert_eq!(stats.extensions_computed, 0);
        assert!(!chunk.join("pairs").exists());
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn parallel_pair_batches_are_consumed_once_and_flush_worker_members() {
        let root = std::env::temp_dir().join(format!(
            "diamond-align-parallel-{}-{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        let job = Job::new(10, 1, 1 << 20, &root, 3).unwrap();
        let chunk = job.base_dir(None).join("chunks/0");
        fs::create_dir_all(&chunk).unwrap();
        fs::write(
            chunk.join("bucket.tsv"),
            super::super::storage::manifest_line(&chunk.join("fixture.fa"), &[2]).unwrap(),
        )
        .unwrap();
        let mut pairs = Vec::new();
        for _ in 0..32 {
            pairs.extend_from_slice(&1usize.to_ne_bytes());
            pairs.extend_from_slice(&4u64.to_ne_bytes());
            pairs.extend_from_slice(&9u64.to_ne_bytes());
        }
        fs::write(chunk.join("pairs"), pairs).unwrap();
        let matrix = Arc::new(ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap());
        let (table, stats) = align_owned(
            &job,
            1,
            10,
            matrix,
            &AlignConfig {
                threads: 4,
                swipe_task_size: 1,
                alignment_output: Some(root.join("aln")),
                ..Default::default()
            },
            &mut FixtureReader,
        )
        .unwrap();
        assert_eq!(stats.extensions_computed, 32);
        assert_eq!(stats.hits_filtered, 32);
        let files: Vec<_> = fs::read_dir(&root)
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| p.file_name().unwrap().to_str().unwrap().starts_with("aln."))
            .collect();
        assert!(!files.is_empty() && files.len() <= 4);
        let lines: usize = files
            .iter()
            .map(|p| fs::read_to_string(p).unwrap().lines().count())
            .sum();
        assert_eq!(lines, 32);
        let records: u64 = table
            .iter()
            .map(|b| {
                super::super::storage::manifest_records(&b.path)
                    .unwrap()
                    .iter()
                    .map(|(_, n)| *n)
                    .sum::<u64>()
            })
            .sum();
        assert_eq!(records, 64); // The one HSP emits both Source directions.
        assert!(!chunk.join("pairs").exists());
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn align_consumes_chunk_pairs_and_emits_radix_manifests() {
        let root = std::env::temp_dir().join(format!(
            "diamond-external-align-{}-{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        let job = Job::new(10, 1, 1 << 20, &root, 3).unwrap();
        let chunk = job.base_dir(None).join("chunks/0");
        fs::create_dir_all(&chunk).unwrap();
        fs::write(
            chunk.join("bucket.tsv"),
            super::super::storage::manifest_line(&chunk.join("fixture.fa"), &[2]).unwrap(),
        )
        .unwrap();
        let mut pairs = Vec::new();
        pairs.extend_from_slice(&1usize.to_ne_bytes());
        pairs.extend_from_slice(&4u64.to_ne_bytes());
        pairs.extend_from_slice(&9u64.to_ne_bytes());
        fs::write(chunk.join("pairs"), pairs).unwrap();
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap();
        let (table, stats) = align(
            &job,
            1,
            10,
            &matrix,
            &AlignConfig::default(),
            &mut FixtureReader,
        )
        .unwrap();

        assert_eq!(table.len(), RADIX_COUNT);
        assert_eq!(stats.extensions_computed, 1);
        assert_eq!(stats.hits_filtered, 1);
        assert!(!chunk.join("pairs").exists());
        assert!(table.iter().all(|bucket| bucket.path.exists()));
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }
}
