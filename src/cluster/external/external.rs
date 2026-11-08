//! External-memory clustering orchestration.
//!
//! This mirrors the data contracts and orchestration boundaries in
//! `diamond/src/cluster/external/external.{h,cpp}`. The C++ implementation uses
//! cross-process atomics, radix buckets and worker pools. Rust exposes the
//! seed/radix/pair/alignment units through an explicit pipeline adapter while
//! preserving their on-disk manifests, native records, cleanup, and shared
//! job/output locks.

use crate::basic::value::OId;
use crate::masking::{MaskingAlgo, MaskingStat};
use crate::util::algo::HyperLogLog;
use crate::util::io::source_c_path;
use crate::util::io::{Compressor, InputFile, IoErrorOwner, IoResult, OutputFile};
use crate::util::parallel::Atomic;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use super::pair_table::{Bucket, RadixedTable, RADIX_BITS, RADIX_COUNT};
use super::seed_table::{SequenceVolumeReader, VolumedFile};
use crate::basic::value::AMINO_ACID_ALPHABET;
use crate::cluster::cascaded::helpers::cluster_steps_bytes;
use crate::config::Sensitivity;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ClusterStats {
    pub hits_evalue_filtered: u64,
    pub extensions_computed: u64,
    pub hits_filtered: u64,
    pub seeds_considered: u64,
    pub seeds_indexed: u64,
    pub masking_stat: MaskingStat,
}

impl ClusterStats {
    pub fn add(&mut self, other: &Self) {
        self.hits_evalue_filtered += other.hits_evalue_filtered;
        self.extensions_computed += other.extensions_computed;
        self.hits_filtered += other.hits_filtered;
        self.seeds_considered += other.seeds_considered;
        self.seeds_indexed += other.seeds_indexed;
        self.masking_stat += other.masking_stat;
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PairEntry {
    pub rep_oid: u64,
    pub member_oid: u64,
    pub rep_len: u32,
    pub member_len: u32,
}

impl Ord for PairEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.rep_oid, self.member_oid).cmp(&(other.rep_oid, other.member_oid))
    }
}

impl PartialOrd for PairEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PairEntryShort {
    pub rep_oid: u64,
    pub member_oid: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Edge {
    pub rep_oid: u64,
    pub member_oid: u64,
    pub rep_len: u32,
    pub member_len: u32,
}

impl Ord for Edge {
    fn cmp(&self, other: &Self) -> Ordering {
        self.member_oid
            .cmp(&other.member_oid)
            .then_with(|| other.rep_len.cmp(&self.rep_len))
            .then_with(|| self.rep_oid.cmp(&other.rep_oid))
    }
}

impl PartialOrd for Edge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Assignment {
    pub member_oid: u64,
    pub rep_oid: u64,
}

/// The packed 12-byte chunk-table record used by C++.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ChunkTableEntry {
    pub oid: OId,
    pub chunk: i32,
}

impl ChunkTableEntry {
    pub const ENCODED_SIZE: usize = 12;

    pub fn new(oid: OId, chunk: i32) -> Self {
        Self { oid, chunk }
    }

    pub fn key(&self) -> OId {
        self.oid
    }

    pub fn serialize(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.oid.to_ne_bytes());
        out.extend_from_slice(&self.chunk.to_ne_bytes());
    }

    pub fn deserialize(bytes: &[u8]) -> Result<Self, String> {
        if bytes.len() < Self::ENCODED_SIZE {
            return Err("Short ChunkTableEntry".to_string());
        }
        Ok(Self {
            oid: u64::from_ne_bytes(bytes[..8].try_into().unwrap()),
            chunk: i32::from_ne_bytes(bytes[8..12].try_into().unwrap()),
        })
    }
}

impl Ord for ChunkTableEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.oid, self.chunk).cmp(&(other.oid, other.chunk))
    }
}

impl PartialOrd for ChunkTableEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone)]
pub struct SizeCounter {
    pub hll: HyperLogLog,
}

impl Default for SizeCounter {
    fn default() -> Self {
        Self {
            hll: HyperLogLog::with_default_precision(),
        }
    }
}

impl SizeCounter {
    pub fn add(&mut self, oid: OId, len: u32) {
        let begin = (oid as i64) << 17;
        let end = begin + (i64::from(len) + 63) / 64;
        for value in begin..end {
            self.hll.add(value);
        }
    }
}

#[derive(Debug)]
pub struct Job {
    pub max_oid: OId,
    pub volumes: usize,
    pub mem_limit: u64,
    root_dir: PathBuf,
    worker_id: i64,
    round: i32,
    round_count: i32,
    input_count: Vec<u64>,
    start: Instant,
    io_errors: IoErrorOwner,
}

impl Job {
    pub fn new(
        max_oid: OId,
        volumes: usize,
        mem_limit: u64,
        root_dir: impl AsRef<Path>,
        worker_id: i64,
    ) -> Result<Self, String> {
        Self::new_owned(
            max_oid,
            volumes,
            mem_limit,
            root_dir,
            worker_id,
            IoErrorOwner::default(),
        )
    }

    pub fn new_owned(
        max_oid: OId,
        volumes: usize,
        mem_limit: u64,
        root_dir: impl AsRef<Path>,
        worker_id: i64,
        io_errors: IoErrorOwner,
    ) -> Result<Self, String> {
        let root_dir = root_dir.as_ref().to_path_buf();
        fs::create_dir_all(source_c_path(&(root_dir.join("round0"))).as_ref())
            .map_err(|error| error.to_string())?;
        Ok(Self {
            max_oid,
            volumes,
            mem_limit,
            root_dir,
            worker_id,
            round: 0,
            round_count: 0,
            input_count: Vec::new(),
            start: Instant::now(),
            io_errors,
        })
    }

    /// Construct a process-shared job and claim the next file-backed worker
    /// id, matching the C++ constructor.
    pub fn new_shared(
        max_oid: OId,
        volumes: usize,
        mem_limit: u64,
        root_dir: impl AsRef<Path>,
    ) -> Result<Self, String> {
        Self::new_shared_owned(
            max_oid,
            volumes,
            mem_limit,
            root_dir,
            IoErrorOwner::default(),
        )
    }

    pub fn new_shared_owned(
        max_oid: OId,
        volumes: usize,
        mem_limit: u64,
        root_dir: impl AsRef<Path>,
        io_errors: IoErrorOwner,
    ) -> Result<Self, String> {
        let root_dir = root_dir.as_ref().to_path_buf();
        fs::create_dir_all(source_c_path(&root_dir).as_ref()).map_err(|error| error.to_string())?;
        let worker_id = Atomic::new(root_dir.join("worker_id")).fetch_add_one()?;
        Self::new_owned(max_oid, volumes, mem_limit, root_dir, worker_id, io_errors)
    }

    pub fn io_errors(&self) -> IoErrorOwner {
        self.io_errors.clone()
    }

    pub fn capture_io<T>(&self, result: std::io::Result<T>) -> Result<T, String> {
        self.io_errors.capture(result)
    }

    pub fn io<T>(&self, result: IoResult<T>) -> Result<T, String> {
        self.capture_io(result.map_err(std::io::Error::other))
    }

    pub fn open_input(&self, path: &Path, flags: i32) -> Result<InputFile, String> {
        self.io(InputFile::with_source_profile_path(path, flags, 65_536))
    }

    pub fn create_output(&self, path: &Path) -> Result<OutputFile, String> {
        self.io(OutputFile::with_source_profile_path(
            path,
            Compressor::None,
            "wb",
            65_536,
        ))
    }

    pub fn read_plain_file(&self, path: &Path) -> Result<Vec<u8>, String> {
        let mut input = self.open_input(path, InputFile::NO_AUTODETECT)?;
        let mut bytes = Vec::new();
        self.capture_io(std::io::Read::read_to_end(&mut input, &mut bytes))?;
        self.io(input.close())?;
        Ok(bytes)
    }

    pub fn worker_id(&self) -> i64 {
        self.worker_id
    }

    pub fn root_dir(&self) -> &Path {
        &self.root_dir
    }

    pub fn base_dir(&self, round: Option<i32>) -> PathBuf {
        self.root_dir
            .join(format!("round{}", round.unwrap_or(self.round)))
    }

    pub fn log(&self, message: &str) -> Result<String, String> {
        let line = format!(
            "[{}, {}] {}\n",
            self.worker_id,
            self.start.elapsed().as_secs(),
            message
        );
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(source_c_path(&(self.root_dir.join("diamond_job.log"))).as_ref())
            .map_err(|error| error.to_string())?;
        file.write_all(line.as_bytes())
            .map_err(|error| error.to_string())?;
        Ok(line)
    }

    pub fn log_stats(&self, stats: &ClusterStats) -> Result<Vec<String>, String> {
        let messages = vec![
            format!(
                "Masked letters:   tantan: {}  seg: {}  motif: {}",
                stats.masking_stat.get(MaskingAlgo::Tantan),
                stats.masking_stat.get(MaskingAlgo::Seg),
                stats.masking_stat.get(MaskingAlgo::Motif)
            ),
            format!("Seeds considered: {}", stats.seeds_considered),
            format!("Seeds indexed: {}", stats.seeds_indexed),
            format!("Extensions computed: {}", stats.extensions_computed),
            format!(
                "Alignments passing e-value filter: {}",
                stats.hits_evalue_filtered
            ),
            format!("Alignments passing all filters: {}", stats.hits_filtered),
        ];
        for message in &messages {
            self.log(message)?;
        }
        Ok(messages)
    }

    pub fn next_round(&mut self) -> Result<(), String> {
        self.round += 1;
        fs::create_dir_all(source_c_path(&(self.base_dir(None))).as_ref())
            .map_err(|error| error.to_string())
    }

    pub fn round_index(&self) -> i32 {
        self.round
    }

    pub fn set_round(&mut self, input_count: u64) {
        self.input_count.push(input_count);
    }

    pub fn sparse_input_count(&self, round: usize) -> u64 {
        self.input_count[round]
    }

    pub fn set_round_count(&mut self, count: i32) {
        self.round_count = count;
    }

    pub fn last_round(&self) -> bool {
        self.round == self.round_count - 1
    }
}

pub fn base_path(file_name: &Path) -> Result<PathBuf, String> {
    let suffix = Path::new("0").join("bucket.tsv");
    if !file_name.ends_with(&suffix) {
        return Err("base_path".to_string());
    }
    Ok(file_name
        .parent()
        .and_then(Path::parent)
        .unwrap_or_else(|| Path::new(""))
        .to_path_buf())
}

pub struct ClusterChunk {
    pub id: i32,
    pairs_out: OutputFile,
    io_errors: IoErrorOwner,
    closed: bool,
    pub size: SizeCounter,
}

impl ClusterChunk {
    pub fn new(id: i32, chunks_path: &Path) -> Result<Self, String> {
        Self::new_owned(id, chunks_path, IoErrorOwner::default())
    }

    pub fn new_owned(id: i32, chunks_path: &Path, io_errors: IoErrorOwner) -> Result<Self, String> {
        let dir = chunks_path.join(id.to_string());
        fs::create_dir_all(source_c_path(&dir).as_ref()).map_err(|error| error.to_string())?;
        let path = dir.join("pairs");
        let pairs_out = io_errors.capture(
            OutputFile::with_source_profile_path(&path, Compressor::None, "wb", 65_536)
                .map_err(std::io::Error::other),
        )?;
        Ok(Self {
            id,
            pairs_out,
            io_errors,
            closed: false,
            size: SizeCounter::default(),
        })
    }

    pub fn write(
        &mut self,
        pairs: &mut Vec<PairEntryShort>,
        size: &mut SizeCounter,
    ) -> Result<(), String> {
        self.io_errors.capture(
            self.pairs_out
                .write_raw(&(pairs.len() as u64).to_ne_bytes())
                .map_err(std::io::Error::other),
        )?;
        for pair in pairs.iter() {
            self.io_errors.capture(
                self.pairs_out
                    .write_raw(&pair.rep_oid.to_ne_bytes())
                    .and_then(|_| self.pairs_out.write_raw(&pair.member_oid.to_ne_bytes()))
                    .map_err(std::io::Error::other),
            )?;
        }
        pairs.clear();
        self.size.hll.merge(&size.hll).map_err(str::to_string)?;
        size.hll = HyperLogLog::with_default_precision();
        Ok(())
    }
    fn close(&mut self) -> Result<(), String> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.io_errors
            .capture(self.pairs_out.close().map_err(std::io::Error::other))
    }
}

impl Drop for ClusterChunk {
    fn drop(&mut self) {
        let _ = self.close();
    }
}

#[derive(Debug)]
pub struct ChunkBuild {
    pub table: Vec<ChunkTableEntry>,
    pub chunks: Vec<Vec<PairEntryShort>>,
}

pub fn build_chunk_table_in_memory(
    job: &Job,
    pairs: &[PairEntry],
    max_chunk_size: i64,
) -> Result<ChunkBuild, String> {
    let mut pairs = pairs.to_vec();
    pairs.sort_unstable();
    let mut table = Vec::new();
    let mut chunks = vec![Vec::new()];
    let mut seen_in_chunk = BTreeSet::new();
    let mut estimated = 0i64;
    let mut index = 0usize;
    while index < pairs.len() {
        let rep_oid = pairs[index].rep_oid;
        let rep_len = pairs[index].rep_len;
        let mut end = index + 1;
        while end < pairs.len() && pairs[end].rep_oid == rep_oid {
            end += 1;
        }
        let mut members = BTreeMap::new();
        for pair in &pairs[index..end] {
            members.entry(pair.member_oid).or_insert(pair.member_len);
        }
        let group_size = i64::from(rep_len) + members.values().copied().map(i64::from).sum::<i64>();
        if estimated > 0 && estimated + group_size > max_chunk_size.max(1) {
            chunks.push(Vec::new());
            seen_in_chunk.clear();
            estimated = 0;
        }
        let chunk = (chunks.len() - 1) as i32;
        if seen_in_chunk.insert(rep_oid) {
            table.push(ChunkTableEntry::new(rep_oid, chunk));
        }
        for (&member_oid, &member_len) in &members {
            if seen_in_chunk.insert(member_oid) {
                table.push(ChunkTableEntry::new(member_oid, chunk));
            }
            chunks[chunk as usize].push(PairEntryShort {
                rep_oid,
                member_oid,
            });
            estimated += i64::from(member_len);
        }
        estimated += i64::from(rep_len);
        index = end;
    }
    table.sort_unstable();
    fs::create_dir_all(source_c_path(&(job.base_dir(None).join("chunk_table"))).as_ref())
        .map_err(|error| error.to_string())?;
    Ok(ChunkBuild { table, chunks })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalSequence {
    pub oid: OId,
    pub sequence: Vec<u8>,
}

pub fn build_chunks_in_memory(
    job: &Job,
    sequences: &[ExternalSequence],
    chunk_table: &[ChunkTableEntry],
    chunk_count: usize,
) -> Result<Vec<PathBuf>, String> {
    let chunks_dir = job.base_dir(None).join("chunks");
    fs::create_dir_all(source_c_path(&chunks_dir).as_ref()).map_err(|error| error.to_string())?;
    let mut wanted = BTreeMap::<OId, BTreeSet<i32>>::new();
    for entry in chunk_table {
        if entry.chunk < 0 || entry.chunk as usize >= chunk_count {
            return Err("Invalid chunk id".to_string());
        }
        wanted.entry(entry.oid).or_default().insert(entry.chunk);
    }
    let mut writers = Vec::with_capacity(chunk_count);
    let mut paths = Vec::with_capacity(chunk_count);
    for id in 0..chunk_count {
        let dir = chunks_dir.join(id.to_string());
        fs::create_dir_all(source_c_path(&dir).as_ref()).map_err(|error| error.to_string())?;
        let path = dir.join("sequences.fasta");
        writers.push(BufWriter::new(
            File::create(source_c_path(&path).as_ref()).map_err(|error| error.to_string())?,
        ));
        paths.push(path);
    }
    let mut sequences = sequences.iter().collect::<Vec<_>>();
    sequences.sort_by_key(|sequence| sequence.oid);
    for sequence in sequences {
        if let Some(chunks) = wanted.get(&sequence.oid) {
            for &chunk in chunks {
                writeln!(writers[chunk as usize], ">{0}", sequence.oid)
                    .and_then(|_| writers[chunk as usize].write_all(&sequence.sequence))
                    .and_then(|_| writers[chunk as usize].write_all(b"\n"))
                    .map_err(|error| error.to_string())?;
            }
        }
    }
    for writer in &mut writers {
        writer.flush().map_err(|error| error.to_string())?;
    }
    Ok(paths)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChunkTableConfig {
    pub threads: usize,
    /// C++ `linclust_chunk_size` after division by 64.
    pub max_chunk_size: i64,
}

impl Default for ChunkTableConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            max_chunk_size: 16 * 1024 * 1024,
        }
    }
}

/// Storage-backed C++ `build_chunk_table`.
pub fn build_chunk_table(
    job: &Job,
    pair_table: &RadixedTable,
    max_oid: OId,
    config: ChunkTableConfig,
) -> Result<(RadixedTable, usize), String> {
    if config.threads == 0 {
        return Err("Thread count must be positive.".to_owned());
    }
    let max_chunk_size = config.max_chunk_size.max(1);
    let max_processed = (max_chunk_size / config.threads as i64 / 16).clamp(1, 262_144);
    let shift = bit_length(max_oid).saturating_sub(RADIX_BITS);
    let table_dir = job.base_dir(None).join("chunk_table");
    let chunks_dir = job.base_dir(None).join("chunks");
    fs::create_dir_all(source_c_path(&table_dir).as_ref()).map_err(|error| error.to_string())?;
    fs::create_dir_all(source_c_path(&chunks_dir).as_ref()).map_err(|error| error.to_string())?;
    let mut table_buckets = super::storage::BucketWriter::new_owned(
        &table_dir,
        job.worker_id(),
        RADIX_COUNT,
        job.io_errors(),
    )?;
    let queue = Atomic::new(table_dir.join("queue"));
    let next_chunk = Atomic::new(table_dir.join("next_chunk"));
    let mut chunk_id =
        i32::try_from(next_chunk.fetch_add_one()?).map_err(|_| "Chunk id overflow.".to_owned())?;
    let mut current = ClusterChunk::new_owned(chunk_id, &chunks_dir, job.io_errors())?;
    let mut buckets_processed = 0_i64;

    loop {
        let bucket_index = usize::try_from(queue.fetch_add_one()?)
            .map_err(|_| "Negative bucket id.".to_owned())?;
        let Some(bucket) = pair_table.get(bucket_index) else {
            break;
        };
        let (mut pairs, storage) = read_pair_bucket(job, bucket)?;
        job.log(&format!(
            "Building chunk table. Bucket={}/{} Records={} Size={}",
            bucket_index + 1,
            pair_table.len(),
            pairs.len(),
            pairs.len() * 24
        ))?;
        pairs.sort_unstable();
        let mut begin = 0usize;
        while begin < pairs.len() {
            let rep_oid = pairs[begin].rep_oid;
            let rep_len = pairs[begin].rep_len;
            let mut end = begin + 1;
            while end < pairs.len() && pairs[end].rep_oid == rep_oid {
                end += 1;
            }
            push_chunk_entry(&mut table_buckets, shift, rep_oid, current.id)?;
            let mut size = SizeCounter::default();
            size.add(rep_oid, rep_len);
            let mut processed = i64::from(rep_len);
            let mut pair_buffer = Vec::new();
            let mut previous_member = None;
            for pair in &pairs[begin..end] {
                if previous_member == Some(pair.member_oid) {
                    continue;
                }
                previous_member = Some(pair.member_oid);
                push_chunk_entry(&mut table_buckets, shift, pair.member_oid, current.id)?;
                size.add(pair.member_oid, pair.member_len);
                pair_buffer.push(PairEntryShort {
                    rep_oid,
                    member_oid: pair.member_oid,
                });
                processed += i64::from(pair.member_len);
                if processed >= max_processed {
                    current.write(&mut pair_buffer, &mut size)?;
                    processed = 0;
                    if current.size.hll.estimate() >= max_chunk_size {
                        current.close()?;
                        drop(current);
                        chunk_id = i32::try_from(next_chunk.fetch_add_one()?)
                            .map_err(|_| "Chunk id overflow.".to_owned())?;
                        current = ClusterChunk::new_owned(chunk_id, &chunks_dir, job.io_errors())?;
                        push_chunk_entry(&mut table_buckets, shift, rep_oid, current.id)?;
                        size.add(rep_oid, rep_len);
                        processed += i64::from(rep_len);
                    }
                }
            }
            current.write(&mut pair_buffer, &mut size)?;
            begin = end;
        }
        if current.size.hll.estimate() >= max_chunk_size {
            current.close()?;
            drop(current);
            chunk_id = i32::try_from(next_chunk.fetch_add_one()?)
                .map_err(|_| "Chunk id overflow.".to_owned())?;
            current = ClusterChunk::new_owned(chunk_id, &chunks_dir, job.io_errors())?;
        }
        remove_manifest_storage(&bucket.path, &storage);
        buckets_processed += 1;
    }
    // Source publishes/closes the chunk table before releasing its current
    // chunk, then enters the completion-counter wait.
    let table = table_buckets.finish()?;
    current.close()?;
    drop(current);
    let finished = Atomic::new(table_dir.join("finished"));
    finished.fetch_add(buckets_processed)?;
    finished.await_value(pair_table.len() as i64)?;
    let chunk_count = usize::try_from(next_chunk.get()?)
        .map_err(|_| "Chunk count does not fit memory.".to_owned())?;
    Ok((table, chunk_count))
}

/// Storage-backed C++ `build_chunks`.
pub fn build_chunks<R: SequenceVolumeReader>(
    job: &Job,
    database: &VolumedFile,
    chunk_table: &RadixedTable,
    chunk_count: usize,
    reader: &mut R,
) -> Result<Vec<PathBuf>, String> {
    let chunks_dir = job.base_dir(None).join("chunks");
    fs::create_dir_all(source_c_path(&chunks_dir).as_ref()).map_err(|error| error.to_string())?;
    let queue = Atomic::new(chunks_dir.join("queue"));
    let mut output = super::storage::BucketWriter::new_owned(
        &chunks_dir,
        job.worker_id(),
        chunk_count,
        job.io_errors(),
    )?;
    let mut buckets_processed = 0_i64;
    let traits = crate::basic::value::ValueTraits::new(
        AMINO_ACID_ALPHABET,
        crate::basic::value::MASK_LETTER,
        b"UO-",
        crate::basic::value::SequenceType::AminoAcid,
    );
    loop {
        let bucket_index = usize::try_from(queue.fetch_add_one()?)
            .map_err(|_| "Negative bucket id.".to_owned())?;
        let Some(bucket) = chunk_table.get(bucket_index) else {
            break;
        };
        let (mut entries, storage) = read_chunk_table_bucket(job, bucket)?;
        entries.sort_unstable();
        job.log(&format!(
            "Building chunks. Bucket={}/{} Records={} Size={}",
            bucket_index + 1,
            chunk_table.len(),
            entries.len(),
            entries.len() * ChunkTableEntry::ENCODED_SIZE
        ))?;
        // Retain only one radix bucket's OID-to-chunk mapping, as in Source.
        let mut wanted = BTreeMap::<OId, BTreeSet<i32>>::new();
        for entry in entries {
            if entry.chunk < 0 || entry.chunk as usize >= chunk_count {
                return Err(format!("Invalid chunk id: {}", entry.chunk));
            }
            wanted.entry(entry.oid).or_default().insert(entry.chunk);
        }
        if let (Some((&first, _)), Some((&last, _))) =
            (wanted.first_key_value(), wanted.last_key_value())
        {
            let mut buffer = crate::util::text_buffer::TextBuffer::new();
            for volume in &database.0 {
                if volume.oid_end != 0 && (volume.oid_end <= first || volume.oid_begin > last) {
                    continue;
                }
                // Source checks file_oid < oid_end before read_seq, then uses
                // the current record's numeric title only in later rounds.
                let sequential_oid = std::cell::Cell::new(volume.oid_begin);
                let oid_end = last.wrapping_add(1);
                reader.visit_volume_bytes_while(
                    &volume.path,
                    &mut || sequential_oid.get() < oid_end,
                    &mut |record| {
                        let oid = if job.round_index() > 0 {
                            parse_oid_like_atoll(&record.id)
                        } else {
                            sequential_oid.get()
                        };
                        sequential_oid.set(oid.wrapping_add(1));
                        let Some(chunks) = wanted.get(&oid) else {
                            return Ok(());
                        };
                        buffer.clear();
                        crate::util::sequence::format(
                            &record.sequence,
                            &oid.to_string(),
                            None,
                            &mut buffer,
                            "fasta",
                            &traits,
                            160,
                        )?;
                        for &chunk in chunks {
                            // Source BufferArray's byte-string overload counts
                            // FASTA bytes, not sequences, in this volume manifest.
                            output.write_count(
                                chunk as usize,
                                buffer.data(),
                                buffer.data().len() as u64,
                            )?;
                        }
                        Ok(())
                    },
                )?;
            }
        }
        remove_manifest_storage(&bucket.path, &storage);
        buckets_processed += 1;
    }
    let table = output.finish()?;
    let outputs = table.iter().map(|bucket| bucket.path.clone()).collect();
    let finished = Atomic::new(chunks_dir.join("finished"));
    finished.fetch_add(buckets_processed)?;
    finished.await_value(chunk_table.len() as i64)?;
    Ok(outputs)
}

fn push_chunk_entry(
    buckets: &mut super::storage::BucketWriter,
    shift: u32,
    oid: OId,
    chunk: i32,
) -> Result<(), String> {
    let radix = usize::try_from(oid >> shift)
        .map_err(|_| "Chunk table radix does not fit memory.".to_owned())?;
    let mut bytes = Vec::with_capacity(ChunkTableEntry::ENCODED_SIZE);
    ChunkTableEntry::new(oid, chunk).serialize(&mut bytes);
    buckets.write(radix, &bytes)
}

fn read_pair_bucket(
    job: &Job,
    bucket: &Bucket,
) -> Result<(Vec<PairEntry>, Vec<(PathBuf, usize)>), String> {
    let storage = read_manifest_owned(&bucket.path, &job.io_errors())?;
    let mut pairs = Vec::new();
    for (path, count) in &storage {
        let mut input = job.open_input(path, 0)?;
        for _ in 0..*count {
            pairs.push(PairEntry {
                rep_oid: job.io(input.read_u64())?,
                member_oid: job.io(input.read_u64())?,
                rep_len: job.io(input.read_u32())?,
                member_len: job.io(input.read_u32())?,
            });
        }
    }
    Ok((pairs, storage))
}

fn read_chunk_table_bucket(
    job: &Job,
    bucket: &Bucket,
) -> Result<(Vec<ChunkTableEntry>, Vec<(PathBuf, usize)>), String> {
    let storage = read_manifest_owned(&bucket.path, &job.io_errors())?;
    let mut entries = Vec::new();
    for (path, count) in &storage {
        let mut input = job.open_input(path, 0)?;
        for _ in 0..*count {
            entries.push(ChunkTableEntry {
                oid: job.io(input.read_u64())?,
                chunk: job.io(input.read_i32())?,
            });
        }
    }
    Ok((entries, storage))
}

fn read_manifest(path: &Path) -> Result<Vec<(PathBuf, usize)>, String> {
    read_manifest_owned(path, &IoErrorOwner::default())
}
pub(crate) fn read_manifest_owned(
    path: &Path,
    errors: &IoErrorOwner,
) -> Result<Vec<(PathBuf, usize)>, String> {
    let bytes = fs::read(source_c_path(&(path)).as_ref()).map_err(|error| error.to_string())?;
    bytes
        .split(|&b| b == b'\n')
        .filter(|line| !line.iter().all(u8::is_ascii_whitespace))
        .map(|line| {
            let mut fields = manifest_fields(line);
            let path = fields.next().ok_or("Format error in VolumedFile")?;
            let count = parse_manifest_number(fields.next().ok_or("Format error in VolumedFile")?)?;
            Ok((
                PathBuf::from(errors.capture(
                    crate::util::io::source_bytes_os(path).map_err(std::io::Error::other),
                )?),
                count,
            ))
        })
        .collect()
}

fn manifest_fields(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(u8::is_ascii_whitespace)
        .filter(|field| !field.is_empty())
}
fn parse_manifest_number<T: std::str::FromStr>(field: &[u8]) -> Result<T, String> {
    std::str::from_utf8(field)
        .ok()
        .and_then(|field| field.parse().ok())
        .ok_or_else(|| "Format error in VolumedFile".to_owned())
}

fn open_input(path: &Path) -> Result<InputFile, String> {
    InputFile::source_new_path(path, 0).map_err(|error| error.to_string())
}

fn remove_manifest_storage(path: &Path, storage: &[(PathBuf, usize)]) {
    for (volume, _) in storage {
        let _ = fs::remove_file(source_c_path(&(volume)).as_ref());
    }
    let _ = fs::remove_file(source_c_path(&(path)).as_ref());
    if let Some(parent) = path.parent() {
        let _ = fs::remove_dir(source_c_path(&(parent)).as_ref());
    }
}

fn bit_length(value: OId) -> u32 {
    OId::BITS - value.leading_zeros()
}

fn parse_oid_like_atoll(value: impl AsRef<[u8]>) -> OId {
    let value = crate::util::sequence::c_string_bytes(value.as_ref());
    let value = &value[value
        .iter()
        .take_while(|byte| byte.is_ascii_whitespace())
        .count()..];
    let (negative, digits) = match value.first() {
        Some(b'-') => (true, &value[1..]),
        Some(b'+') => (false, &value[1..]),
        _ => (false, value),
    };
    let end = digits
        .iter()
        .position(|byte| !byte.is_ascii_digit())
        .unwrap_or(digits.len());
    let magnitude = if end == 0 {
        0
    } else {
        digits[..end]
            .iter()
            .try_fold(0_u64, |n, digit| {
                n.checked_mul(10)?.checked_add(u64::from(digit - b'0'))
            })
            .unwrap_or(0)
    };
    if negative {
        0_u64.wrapping_sub(magnitude)
    } else {
        magnitude
    }
}

#[derive(Debug)]
pub struct RoundArtifacts {
    pub chunk_table: Vec<ChunkTableEntry>,
    pub chunk_files: Vec<PathBuf>,
}

pub fn round_in_memory(
    job: &mut Job,
    pairs: &[PairEntry],
    sequences: &[ExternalSequence],
    max_chunk_size: i64,
) -> Result<RoundArtifacts, String> {
    job.set_round(sequences.len() as u64);
    job.log(&format!("Starting round {}", job.round_index()))?;
    let build = build_chunk_table_in_memory(job, pairs, max_chunk_size)?;
    let chunk_files = build_chunks_in_memory(job, sequences, &build.table, build.chunks.len())?;
    Ok(RoundArtifacts {
        chunk_table: build.table,
        chunk_files,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExternalRecordKind {
    Seed,
    Pair,
    ChunkTable,
    Edge,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExternalConfig {
    pub output_file: PathBuf,
    pub oid_output: bool,
    pub root_dir: PathBuf,
    pub memory_limit: u64,
    pub threads: usize,
    pub mutual_cover: Option<f64>,
    pub member_cover: f64,
    pub approx_min_id: f64,
    pub sensitivity: Sensitivity,
    pub min_length_ratio: f64,
    pub linclust_chunk_size: u64,
    pub cluster_steps: Vec<Vec<u8>>,
    pub file_buffer_size: usize,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
}

impl Default for ExternalConfig {
    fn default() -> Self {
        Self {
            output_file: PathBuf::new(),
            oid_output: false,
            root_dir: PathBuf::from("diamond-tmp"),
            memory_limit: 16_000_000_000,
            threads: 1,
            mutual_cover: None,
            member_cover: 80.0,
            approx_min_id: 50.0,
            sensitivity: Sensitivity::Default,
            min_length_ratio: 0.0,
            linclust_chunk_size: 10_000_000_000,
            cluster_steps: Vec::new(),
            file_buffer_size: 0,
            query_cover: 0.0,
            subject_cover: 0.0,
            query_or_target_cover: 0.0,
        }
    }
}

/// Explicit adapters for dependencies owned by the seed/pair/radix/alignment
/// translation units. Stateful implementations may retain the shared
/// `PairFileArray` between `begin_pair_table` and `finish_pair_table`.
pub trait ExternalPipeline: SequenceVolumeReader {
    fn io_errors(&self) -> IoErrorOwner {
        IoErrorOwner::default()
    }
    fn shape_count(&self, sensitivity: Sensitivity) -> usize;
    fn build_seed_table(
        &mut self,
        job: &Job,
        volumes: &VolumedFile,
        shape: usize,
    ) -> Result<RadixedTable, String>;
    fn radix_sort(
        &mut self,
        job: &Job,
        table: RadixedTable,
        kind: ExternalRecordKind,
        bits_unsorted: u32,
    ) -> Result<RadixedTable, String>;
    fn begin_pair_table(&mut self, job: &Job) -> Result<(), String>;
    fn build_pair_table(
        &mut self,
        job: &Job,
        seed_table: &RadixedTable,
        shape: usize,
        max_oid: OId,
    ) -> Result<RadixedTable, String>;
    fn finish_pair_table(&mut self) -> Result<(), String>;
    fn align(
        &mut self,
        job: &Job,
        chunk_count: usize,
        max_oid: OId,
    ) -> Result<RadixedTable, String>;
}

/// Full C++ `External::round` orchestration.
pub fn round<B: ExternalPipeline>(
    job: &mut Job,
    volumes: &VolumedFile,
    config: &mut ExternalConfig,
    backend: &mut B,
) -> Result<PathBuf, String> {
    let shape_count = backend.shape_count(config.sensitivity);
    if shape_count == 0 {
        return Err("No seed shapes configured.".to_owned());
    }
    if let Some(mutual_cover) = config.mutual_cover {
        config.min_length_ratio = if config.sensitivity < Sensitivity::Linclust40 {
            (mutual_cover / 100.0 + 0.05).min(1.0)
        } else {
            mutual_cover / 100.0 - 0.05
        };
    }
    job.log(&format!(
        "Starting round {} sensitivity {} {} shapes\n",
        job.round_index(),
        crate::run::config::sensitivity_name(config.sensitivity),
        shape_count
    ))?;
    job.set_round(volumes.0.iter().map(|volume| volume.record_count).sum());
    backend.begin_pair_table(job)?;
    let mut pair_table = RadixedTable::default();
    for shape in 0..shape_count {
        let seed_table = backend.build_seed_table(job, volumes, shape)?;
        let seed_table = backend.radix_sort(
            job,
            seed_table,
            ExternalRecordKind::Seed,
            OId::BITS - RADIX_BITS,
        )?;
        pair_table = backend.build_pair_table(job, &seed_table, shape, volumes_max_oid(volumes))?;
    }
    backend.finish_pair_table()?;
    let pair_table = backend.radix_sort(
        job,
        pair_table,
        ExternalRecordKind::Pair,
        OId::BITS - RADIX_BITS,
    )?;
    let (chunk_table, chunk_count) = build_chunk_table(
        job,
        &pair_table,
        volumes_max_oid(volumes),
        ChunkTableConfig {
            threads: config.threads,
            max_chunk_size: i64::try_from(config.linclust_chunk_size / 64).unwrap_or(i64::MAX),
        },
    )?;
    let chunk_bits = bit_length(volumes_max_oid(volumes)).saturating_sub(RADIX_BITS);
    let chunk_table =
        backend.radix_sort(job, chunk_table, ExternalRecordKind::ChunkTable, chunk_bits)?;
    build_chunks(job, volumes, &chunk_table, chunk_count, backend)?;
    let edges = backend.align(job, chunk_count, volumes_max_oid(volumes))?;
    if config.mutual_cover.is_some() {
        super::cluster::cluster_bidirectional(job, &edges, volumes)
    } else {
        let edges =
            backend.radix_sort(job, edges, ExternalRecordKind::Edge, OId::BITS - RADIX_BITS)?;
        super::cluster::cluster(job, &edges, volumes)
    }
}

/// Complete C++ `External::external` workflow. The returned `None` denotes a
/// non-winning worker at the final output lock.
pub fn external<B: ExternalPipeline>(
    config: &mut ExternalConfig,
    database: &VolumedFile,
    backend: &mut B,
) -> Result<Option<super::output::ExternalOutputSummary>, String> {
    if config.output_file.as_os_str().is_empty() {
        return Err("Option missing: output file (--out/-o)".to_owned());
    }
    if config.threads == 0 {
        return Err("Thread count must be positive.".to_owned());
    }
    config.file_buffer_size = 64 * 1024;
    let mut job = Job::new_shared_owned(
        volumes_max_oid(database),
        database.0.len(),
        config.memory_limit,
        &config.root_dir,
        backend.io_errors(),
    )?;
    if job.worker_id() == 0 {
        job.log(&format!(
            "{} coverage = {}",
            if config.mutual_cover.is_some() {
                "Bi-directional"
            } else {
                "Uni-directional"
            },
            crate::output::format::format_fixed_float(
                config.mutual_cover.unwrap_or(config.member_cover)
            )
        ))?;
        job.log(&format!(
            "Approx. id = {}",
            crate::output::format::format_fixed_float(config.approx_min_id)
        ))?;
        job.log(&format!("#Volumes = {}", database.0.len()))?;
        job.log(&format!(
            "#Sequences = {}",
            database
                .0
                .iter()
                .map(|volume| volume.record_count)
                .sum::<u64>()
        ))?;
    }
    if let Some(coverage) = config.mutual_cover {
        config.query_or_target_cover = 0.0;
        config.query_cover = coverage;
        config.subject_cover = coverage;
    } else {
        config.query_or_target_cover = config.member_cover;
        config.query_cover = 0.0;
        config.subject_cover = 0.0;
    }
    // Source raises the Windows CRT's FILE limit before creating the round's
    // bucket arrays. The CRT limit is distinct from Rust/Win32 handle limits.
    #[cfg(windows)]
    crate::util::io::cstdio_file::source_external_stdio_limit();
    let steps = if config.cluster_steps.is_empty() {
        cluster_steps_bytes(config.approx_min_id, true)
    } else {
        config.cluster_steps.clone()
    };
    job.set_round_count(steps.len() as i32);
    let mut representatives = PathBuf::new();
    for (index, step) in steps.iter().enumerate() {
        config.sensitivity = job.capture_io(sensitivity_from_step_bytes(step))?;
        let volumes = if index == 0 {
            database.clone()
        } else {
            read_volumed_file_manifest_owned(&representatives, &job.io_errors())?
        };
        representatives = round(&mut job, &volumes, config, backend)?;
        if index + 1 < steps.len() {
            job.next_round()?;
        }
    }
    let output_lock = Atomic::new(job.root_dir().join("output_lock"));
    if output_lock.fetch_add_one()? == 0 {
        let summary = super::output::output(
            &job,
            database,
            &super::output::ExternalOutputConfig {
                output_file: config.output_file.clone(),
                oid_output: config.oid_output,
                threads: config.threads,
            },
        )?;
        Ok(Some(summary))
    } else {
        Ok(None)
    }
}

fn volumes_max_oid(volumes: &VolumedFile) -> OId {
    volumes
        .0
        .iter()
        .filter_map(|volume| volume.oid_end.checked_sub(1))
        .max()
        .unwrap_or(0)
}

pub fn read_volumed_file_manifest(path: &Path) -> Result<VolumedFile, String> {
    read_volumed_file_manifest_owned(path, &IoErrorOwner::default())
}

pub fn read_volumed_file_manifest_owned(
    path: &Path,
    errors: &IoErrorOwner,
) -> Result<VolumedFile, String> {
    // Source uses ifstream here: retain its runtime diagnostic, not the
    // distinct InputFile/FileOpenException class.
    let text = match fs::read(source_c_path(&(path)).as_ref()) {
        Ok(bytes) => bytes,
        Err(_) => {
            let mut message = b"Error opening file ".to_vec();
            message.extend(errors.capture(
                crate::util::io::path_diagnostic_bytes(path).map_err(std::io::Error::other),
            )?);
            return errors.capture(Err(crate::util::io::RawDiagnosticError(message).into_io()));
        }
    };
    let mut volumes = Vec::new();
    let mut next_oid = 0_u64;
    for line in text
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty())
    {
        let fields = manifest_fields(line).collect::<Vec<_>>();
        if fields.len() < 2 {
            return Err("Format error in VolumedFile".to_owned());
        }
        let count = parse_manifest_number::<u64>(fields[1])?;
        let (begin, end) = if fields.len() >= 4 {
            (
                parse_manifest_number::<u64>(fields[2])?,
                parse_manifest_number::<u64>(fields[3])?,
            )
        } else {
            (
                next_oid,
                next_oid.checked_add(count).ok_or("Volume OID overflow")?,
            )
        };
        let (begin, end) = if begin == 0 && end == 0 {
            (
                next_oid,
                next_oid.checked_add(count).ok_or("Volume OID overflow")?,
            )
        } else {
            (begin, end)
        };
        next_oid = next_oid.checked_add(count).ok_or("Volume OID overflow")?;
        volumes.push(super::seed_table::Volume::new(
            PathBuf::from(errors.capture(
                crate::util::io::source_bytes_os(fields[0]).map_err(std::io::Error::other),
            )?),
            begin,
            end,
            count,
        ));
    }
    volumes.sort_by_key(|volume| volume.oid_begin);
    Ok(VolumedFile(volumes))
}

fn sensitivity_from_step(step: &str) -> Result<Sensitivity, String> {
    sensitivity_from_step_bytes(step.as_bytes()).map_err(|error| error.to_string())
}
fn sensitivity_from_step_bytes(step: &[u8]) -> std::io::Result<Sensitivity> {
    crate::basic::config::sensitivity_from_bytes(step.strip_suffix(b"_lin").unwrap_or(step))
}

#[cfg(test)]
mod tests {
    use super::super::pair_table::PairFileArray;
    use super::super::seed_table::{SeedSequence, Volume};
    use super::*;

    fn temp_dir(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-external-{name}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[test]
    fn packed_chunk_entry_roundtrip_and_order_are_exact() {
        let mut bytes = Vec::new();
        ChunkTableEntry::new(9, 2).serialize(&mut bytes);
        assert_eq!(bytes.len(), 12);
        assert_eq!(
            ChunkTableEntry::deserialize(&bytes).unwrap(),
            ChunkTableEntry::new(9, 2)
        );
        let mut entries = [
            ChunkTableEntry::new(9, 2),
            ChunkTableEntry::new(8, 7),
            ChunkTableEntry::new(9, 1),
        ];
        entries.sort();
        assert_eq!(
            entries,
            [
                ChunkTableEntry::new(8, 7),
                ChunkTableEntry::new(9, 1),
                ChunkTableEntry::new(9, 2)
            ]
        );
    }

    #[test]
    fn chunk_build_deduplicates_pairs_and_fasta_targets() {
        let dir = temp_dir("chunks");
        let job = Job::new(3, 1, 1024, &dir, 0).unwrap();
        let pairs = [
            PairEntry {
                rep_oid: 0,
                member_oid: 1,
                rep_len: 4,
                member_len: 3,
            },
            PairEntry {
                rep_oid: 0,
                member_oid: 1,
                rep_len: 4,
                member_len: 3,
            },
            PairEntry {
                rep_oid: 2,
                member_oid: 3,
                rep_len: 2,
                member_len: 2,
            },
        ];
        let build = build_chunk_table_in_memory(&job, &pairs, 7).unwrap();
        assert_eq!(build.chunks.len(), 2);
        assert_eq!(
            build.chunks[0],
            vec![PairEntryShort {
                rep_oid: 0,
                member_oid: 1
            }]
        );
        let sequences = [
            ExternalSequence {
                oid: 0,
                sequence: b"AAAA".to_vec(),
            },
            ExternalSequence {
                oid: 1,
                sequence: b"BBB".to_vec(),
            },
            ExternalSequence {
                oid: 2,
                sequence: b"CC".to_vec(),
            },
            ExternalSequence {
                oid: 3,
                sequence: b"DD".to_vec(),
            },
        ];
        let paths = build_chunks_in_memory(&job, &sequences, &build.table, 2).unwrap();
        assert_eq!(
            fs::read_to_string(&paths[0]).unwrap(),
            ">0\nAAAA\n>1\nBBB\n"
        );
        assert_eq!(fs::read_to_string(&paths[1]).unwrap(), ">2\nCC\n>3\nDD\n");
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn job_round_paths_logging_and_base_path_match_upstream() {
        let dir = temp_dir("job");
        let mut job = Job::new(10, 2, 4096, &dir, 7).unwrap();
        job.set_round_count(2);
        assert_eq!(job.worker_id(), 7);
        assert!(!job.last_round());
        job.log("hello").unwrap();
        assert!(fs::read_to_string(dir.join("diamond_job.log"))
            .unwrap()
            .contains("[7, 0] hello"));
        job.next_round().unwrap();
        assert!(job.last_round());
        assert!(job.base_dir(None).ends_with("round1"));
        assert_eq!(
            base_path(Path::new("tmp/table/0/bucket.tsv")).unwrap(),
            PathBuf::from("tmp/table")
        );
        assert_eq!(
            base_path(Path::new("tmp/table/bucket.tsv")).unwrap_err(),
            "base_path"
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn cluster_chunk_writes_native_count_and_pair_records() {
        let dir = temp_dir("pairs");
        let mut chunk = ClusterChunk::new(0, &dir).unwrap();
        let mut pairs = vec![PairEntryShort {
            rep_oid: 4,
            member_oid: 9,
        }];
        let mut size = SizeCounter::default();
        size.add(4, 64);
        chunk.write(&mut pairs, &mut size).unwrap();
        drop(chunk);
        let bytes = fs::read(dir.join("0/pairs")).unwrap();
        assert_eq!(bytes.len(), 8 + 16);
        assert_eq!(u64::from_ne_bytes(bytes[..8].try_into().unwrap()), 1);
        assert_eq!(u64::from_ne_bytes(bytes[8..16].try_into().unwrap()), 4);
        assert_eq!(u64::from_ne_bytes(bytes[16..24].try_into().unwrap()), 9);
        fs::remove_dir_all(dir).unwrap();
    }

    #[derive(Default)]
    struct FixtureReader {
        records: Vec<SeedSequence>,
    }

    impl SequenceVolumeReader for FixtureReader {
        fn read_volume(&mut self, _: &Path) -> Result<Vec<SeedSequence>, String> {
            Ok(self.records.clone())
        }
    }

    #[derive(Default)]
    struct RecordingPipeline {
        records: Vec<SeedSequence>,
        events: Vec<String>,
    }

    impl SequenceVolumeReader for RecordingPipeline {
        fn read_volume(&mut self, _: &Path) -> Result<Vec<SeedSequence>, String> {
            Ok(self.records.clone())
        }
    }

    impl ExternalPipeline for RecordingPipeline {
        fn shape_count(&self, _: Sensitivity) -> usize {
            1
        }

        fn build_seed_table(
            &mut self,
            job: &Job,
            _: &VolumedFile,
            shape: usize,
        ) -> Result<RadixedTable, String> {
            self.events.push(format!("seed:{shape}"));
            let path = job.base_dir(None).join(format!("seed-{shape}.tsv"));
            fs::write(&path, []).map_err(|error| error.to_string())?;
            Ok(RadixedTable(vec![Bucket::new(path, Some(0))]))
        }

        fn radix_sort(
            &mut self,
            _: &Job,
            table: RadixedTable,
            kind: ExternalRecordKind,
            bits_unsorted: u32,
        ) -> Result<RadixedTable, String> {
            self.events.push(format!("radix:{kind:?}:{bits_unsorted}"));
            Ok(table)
        }

        fn begin_pair_table(&mut self, _: &Job) -> Result<(), String> {
            self.events.push("pair:begin".into());
            Ok(())
        }

        fn build_pair_table(
            &mut self,
            job: &Job,
            _: &RadixedTable,
            shape: usize,
            max_oid: OId,
        ) -> Result<RadixedTable, String> {
            self.events.push(format!("pair:build:{shape}:{max_oid}"));
            PairFileArray::new(job.base_dir(None).join("pair-fixture"), job.worker_id())?.finish()
        }

        fn finish_pair_table(&mut self) -> Result<(), String> {
            self.events.push("pair:finish".into());
            Ok(())
        }

        fn align(
            &mut self,
            job: &Job,
            chunk_count: usize,
            max_oid: OId,
        ) -> Result<RadixedTable, String> {
            self.events.push(format!("align:{chunk_count}:{max_oid}"));
            // Model align.cpp's actual stage owner as well as its edge table;
            // cluster claims buckets through alignments/queue.
            fs::create_dir_all(job.base_dir(None).join("alignments"))
                .map_err(|error| error.to_string())?;
            let path = job.base_dir(None).join("edges.tsv");
            fs::write(&path, []).map_err(|error| error.to_string())?;
            Ok(RadixedTable(vec![Bucket::new(path, Some(0))]))
        }
    }

    #[test]
    fn chunk_builder_stops_before_bad_tail_in_original_and_numeric_oid_rounds() {
        for round in 0..=1 {
            let dir = temp_dir("chunk-read-boundary");
            let mut job = Job::new(0, 1, 1 << 20, &dir, 2).unwrap();
            if round > 0 {
                job.next_round().unwrap();
            }
            let oid = if round == 0 { 0 } else { 7 };
            let input = dir.join("volume.faa");
            fs::write(&input, b">7\nARND\n>invalid-tail\n!\n").unwrap();
            let mut table_writer = super::super::storage::BucketWriter::new_owned(
                &dir.join("chunk-table"),
                job.worker_id(),
                1,
                job.io_errors(),
            )
            .unwrap();
            let mut entry = Vec::new();
            ChunkTableEntry::new(oid, 0).serialize(&mut entry);
            table_writer.write(0, &entry).unwrap();
            let table = table_writer.finish().unwrap();
            let database = VolumedFile(vec![Volume::new(&input, 0, 10, 2)]);
            let outputs = build_chunks(
                &job,
                &database,
                &table,
                1,
                &mut super::super::seed_table::FastaVolumeReader,
            )
            .unwrap();
            let mut fasta = Vec::new();
            for manifest in outputs {
                for (path, count) in read_manifest(&manifest).unwrap() {
                    let mut source = open_input(&path).unwrap();
                    let mut data = vec![0; count];
                    source.read_values(&mut data).unwrap();
                    fasta.extend(data);
                }
            }
            assert_eq!(fasta, format!(">{oid}\nARND\n").as_bytes());
            drop(table_writer);
            drop(job);
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn storage_pipeline_deduplicates_pairs_materializes_fasta_and_cleans_inputs() {
        let dir = temp_dir("storage-pipeline");
        let job = Job::new(2, 1, 1 << 20, &dir, 4).unwrap();
        let pair_dir = dir.join("input-pairs");
        let mut pair_files = PairFileArray::new(&pair_dir, job.worker_id()).unwrap();
        for pair in [
            PairEntry {
                rep_oid: 0,
                member_oid: 1,
                rep_len: 4,
                member_len: 3,
            },
            PairEntry {
                rep_oid: 0,
                member_oid: 1,
                rep_len: 4,
                member_len: 3,
            },
            PairEntry {
                rep_oid: 0,
                member_oid: 2,
                rep_len: 4,
                member_len: 2,
            },
        ] {
            pair_files.write_msb(pair);
        }
        let pair_table = pair_files.finish().unwrap();
        let source_manifests = pair_table
            .iter()
            .map(|bucket| bucket.path.clone())
            .collect::<Vec<_>>();
        let (chunk_table, chunk_count) = build_chunk_table(
            &job,
            &pair_table,
            2,
            ChunkTableConfig {
                threads: 2,
                max_chunk_size: i64::MAX / 4,
            },
        )
        .unwrap();
        assert_eq!(chunk_count, 1);
        assert!(source_manifests.iter().all(|path| !path.exists()));

        let pair_bytes = fs::read(job.base_dir(None).join("chunks/0/pairs")).unwrap();
        assert_eq!(u64::from_ne_bytes(pair_bytes[..8].try_into().unwrap()), 2);
        assert_eq!(pair_bytes.len(), 8 + 2 * 16);

        let database = VolumedFile(vec![Volume::new("fixture.faa", 0, 3, 3)]);
        let mut reader = FixtureReader {
            records: vec![
                SeedSequence {
                    id: "a".into(),
                    sequence: vec![0; 4],
                },
                SeedSequence {
                    id: "b".into(),
                    sequence: vec![1; 3],
                },
                SeedSequence {
                    id: "c".into(),
                    sequence: vec![2; 2],
                },
            ],
        };
        let chunk_manifests = chunk_table
            .iter()
            .map(|bucket| bucket.path.clone())
            .collect::<Vec<_>>();
        let outputs =
            build_chunks(&job, &database, &chunk_table, chunk_count, &mut reader).unwrap();
        assert!(chunk_manifests.iter().all(|path| !path.exists()));
        // build_chunks returns Source VolumedFile manifests, whose counts
        // are FASTA bytes. Decode the listed compressed volume payloads.
        let volumes = read_manifest(&outputs[0]).unwrap();
        let mut fasta = Vec::new();
        for (path, count) in volumes {
            let mut input = open_input(&path).unwrap();
            let mut bytes = vec![0u8; count];
            input.read_values(&mut bytes).unwrap();
            fasta.extend(bytes);
            assert_eq!(input.read_raw(&mut [0u8; 1]).unwrap(), 0);
        }
        assert_eq!(fasta, b">0\nAAAA\n>1\nRRR\n>2\nNN\n");
        assert!(
            fs::read_to_string(job.base_dir(None).join("chunks/0/bucket.tsv"))
                .unwrap()
                .ends_with("\t21\n")
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn round_runs_the_complete_adapter_pipeline_and_writes_identity_closure() {
        let dir = temp_dir("round-pipeline");
        let mut job = Job::new(2, 1, 1 << 20, &dir, 0).unwrap();
        job.set_round_count(1);
        let volumes = VolumedFile(vec![Volume::new("fixture.faa", 0, 3, 3)]);
        let mut backend = RecordingPipeline {
            records: vec![
                SeedSequence {
                    id: "a".into(),
                    sequence: vec![0; 4],
                },
                SeedSequence {
                    id: "b".into(),
                    sequence: vec![1; 3],
                },
                SeedSequence {
                    id: "c".into(),
                    sequence: vec![2; 2],
                },
            ],
            events: Vec::new(),
        };
        let mut config = ExternalConfig {
            root_dir: dir.clone(),
            approx_min_id: 95.0,
            ..ExternalConfig::default()
        };
        assert_eq!(
            round(&mut job, &volumes, &mut config, &mut backend).unwrap(),
            PathBuf::new()
        );
        assert_eq!(
            backend.events,
            [
                "pair:begin",
                "seed:0",
                "radix:Seed:56",
                "pair:build:0:2",
                "pair:finish",
                "radix:Pair:56",
                "radix:ChunkTable:0",
                "align:1:2",
                "radix:Edge:56",
            ]
        );
        let clustering = fs::read(job.base_dir(None).join("clustering/volume0")).unwrap();
        let oids = clustering
            .chunks_exact(8)
            .map(|bytes| u64::from_ne_bytes(bytes.try_into().unwrap()))
            .collect::<Vec<_>>();
        assert_eq!(oids, [0, 1, 2]);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn external_applies_coverage_runs_rounds_and_wins_output_lock() {
        let dir = temp_dir("external-pipeline");
        let output = dir.join("clusters.tsv");
        let volumes = VolumedFile(vec![Volume::new("fixture.faa", 0, 3, 3)]);
        let mut backend = RecordingPipeline {
            records: vec![
                SeedSequence {
                    id: "a".into(),
                    sequence: vec![0; 4],
                },
                SeedSequence {
                    id: "b".into(),
                    sequence: vec![1; 3],
                },
                SeedSequence {
                    id: "c".into(),
                    sequence: vec![2; 2],
                },
            ],
            events: Vec::new(),
        };
        let mut config = ExternalConfig {
            output_file: output.clone(),
            oid_output: true,
            root_dir: dir.join("job"),
            approx_min_id: 95.0,
            member_cover: 77.0,
            ..ExternalConfig::default()
        };
        let summary = external(&mut config, &volumes, &mut backend)
            .unwrap()
            .expect("first worker owns output");
        assert_eq!(summary.cluster_count, 3);
        assert_eq!(summary.merged, [0, 1, 2]);
        assert_eq!(fs::read_to_string(output).unwrap(), "0\t0\n1\t1\n2\t2\n");
        assert_eq!(config.file_buffer_size, 64 * 1024);
        assert_eq!(config.query_or_target_cover, 77.0);
        assert_eq!(config.query_cover, 0.0);
        assert_eq!(config.subject_cover, 0.0);

        let mut missing_output = ExternalConfig::default();
        assert_eq!(
            external(&mut missing_output, &VolumedFile::default(), &mut backend).unwrap_err(),
            "Option missing: output file (--out/-o)"
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn header_ordering_and_stats_are_preserved() {
        let mut edges = [
            Edge {
                rep_oid: 9,
                member_oid: 2,
                rep_len: 100,
                member_len: 80,
            },
            Edge {
                rep_oid: 3,
                member_oid: 2,
                rep_len: 100,
                member_len: 80,
            },
            Edge {
                rep_oid: 5,
                member_oid: 2,
                rep_len: 120,
                member_len: 80,
            },
        ];
        edges.sort();
        assert_eq!(
            edges.iter().map(|edge| edge.rep_oid).collect::<Vec<_>>(),
            [5, 3, 9]
        );

        let mut total = ClusterStats::default();
        total.add(&ClusterStats {
            seeds_considered: 4,
            hits_filtered: 2,
            ..ClusterStats::default()
        });
        assert_eq!(total.seeds_considered, 4);
        assert_eq!(total.hits_filtered, 2);
    }
}
