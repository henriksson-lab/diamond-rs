//! Pair-table construction from `diamond/src/cluster/external/pair_table.cpp`
//! and its inline `build_pair_table.h` algorithms.

use crate::util::io::source_c_path;
use std::cell::RefCell;
use std::cmp::Ordering;
use std::fs;
use std::ops::{Deref, DerefMut};
use std::path::{Path, PathBuf};

use super::external::{Job, PairEntry};
use super::storage::{BucketWriter, WorkerBucketBuffers};
use crate::util::algo::Partition;
use crate::util::hash::hash64;
use crate::util::io::{CompressedBuffer, InputFile};
use crate::util::parallel::Atomic;
use crate::util::parallel::SimpleThreadPool;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering as ThreadOrdering};
use std::sync::{Arc, Mutex};

pub const RADIX_BITS: u32 = 8;
pub const RADIX_COUNT: usize = 1 << RADIX_BITS;

/// Packed C++ `SeedEntry` (8-byte seed, 8-byte signed OID, 4-byte length).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeedEntry {
    pub seed: u64,
    pub oid: i64,
    pub len: i32,
}

impl SeedEntry {
    pub const ENCODED_SIZE: usize = 20;

    pub const fn new(seed: u64, oid: i64, len: i32) -> Self {
        Self { seed, oid, len }
    }

    pub fn key(self) -> u64 {
        hash64(self.seed)
    }

    /// C++ `SeedEntry::Key::operator()` used for equal-seed merging.
    pub const fn merge_key(self) -> u64 {
        self.seed
    }

    pub fn serialize(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.seed.to_ne_bytes());
        out.extend_from_slice(&self.oid.to_ne_bytes());
        out.extend_from_slice(&self.len.to_ne_bytes());
    }

    pub fn deserialize(bytes: &[u8]) -> Result<Self, String> {
        if bytes.len() < Self::ENCODED_SIZE {
            return Err("Short SeedEntry".to_owned());
        }
        Ok(Self {
            seed: u64::from_ne_bytes(bytes[0..8].try_into().unwrap()),
            oid: i64::from_ne_bytes(bytes[8..16].try_into().unwrap()),
            len: i32::from_ne_bytes(bytes[16..20].try_into().unwrap()),
        })
    }
}

impl Ord for SeedEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.seed
            .cmp(&other.seed)
            .then_with(|| other.len.cmp(&self.len))
            .then_with(|| self.oid.cmp(&other.oid))
    }
}

impl PartialOrd for SeedEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bucket {
    pub path: PathBuf,
    pub records: Option<u64>,
}

impl Bucket {
    pub fn new(path: impl Into<PathBuf>, records: Option<u64>) -> Self {
        Self {
            path: path.into(),
            records,
        }
    }

    pub fn record_count(&self) -> Result<u64, String> {
        self.record_count_owned(&Default::default())
    }
    pub(crate) fn record_count_owned(
        &self,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<u64, String> {
        self.records.ok_or_else(|| {
            super::storage::path_error(errors, b"Record count not set for bucket: ", &self.path)
        })
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RadixedTable(pub Vec<Bucket>);

impl Deref for RadixedTable {
    type Target = [Bucket];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for RadixedTable {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl RadixedTable {
    pub fn from_manifest(path: impl AsRef<Path>) -> Result<Self, String> {
        Self::from_manifest_owned(path, &Default::default())
    }

    pub fn from_manifest_owned(
        path: impl AsRef<Path>,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        Ok(Self(
            super::storage::manifest_records_owned(path.as_ref(), errors)?
                .into_iter()
                .map(|(path, n)| Bucket::new(path, Some(n)))
                .collect(),
        ))
    }

    pub fn max_buckets(&self, mem_limit: u64, record_size: usize) -> Result<u64, String> {
        self.max_buckets_owned(mem_limit, record_size, &Default::default())
    }
    fn max_buckets_owned(
        &self,
        mem_limit: u64,
        record_size: usize,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<u64, String> {
        let mut counts = self
            .iter()
            .map(|bucket| bucket.record_count_owned(errors))
            .collect::<Result<Vec<_>, _>>()?;
        counts.sort_unstable_by(|a, b| b.cmp(a));
        let mut sum = 0_u64;
        for (index, count) in counts.iter().enumerate() {
            sum = sum.saturating_add(count.saturating_mul(record_size as u64));
            if sum >= mem_limit {
                return Ok(if index > 0 { index as u64 } else { 1 });
            }
        }
        Ok(counts.len() as u64)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairTableConfig {
    pub threads: usize,
    pub mutual_cover_present: bool,
    pub min_length_ratio: f64,
}

impl Default for PairTableConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            mutual_cover_present: false,
            min_length_ratio: 0.0,
        }
    }
}

/// FileArray-compatible output owner. Each bucket is emitted as a zlib volume
/// plus the upstream tab-separated `bucket.tsv` volume manifest.
pub struct PairFileArray {
    base_dir: PathBuf,
    worker_id: i64,
    pairs: Vec<Vec<PairEntry>>,
    bounded: Option<Arc<Mutex<BucketWriter>>>,
    legacy_buffers: RefCell<Option<WorkerBucketBuffers>>,
    failure: Arc<Mutex<Option<String>>>,
    errors: crate::util::io::IoErrorOwner,
}

impl std::fmt::Debug for PairFileArray {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        out.debug_struct("PairFileArray")
            .field("base_dir", &self.base_dir)
            .field("worker_id", &self.worker_id)
            .field("pairs", &self.pairs)
            .field("bounded", &self.bounded.is_some())
            .finish()
    }
}

impl PairFileArray {
    /// Eager in-memory compatibility adapter: creates its parent hierarchy.
    /// Actual Source FileArray workflows use the bounded owner after their
    /// explicit, single-level mkdir phase.
    pub fn new(base_dir: impl AsRef<Path>, worker_id: i64) -> Result<Self, String> {
        let base_dir = base_dir.as_ref().to_path_buf();
        for radix in 0..RADIX_COUNT {
            fs::create_dir_all(source_c_path(&(base_dir.join(radix.to_string()))).as_ref())
                .map_err(|error| error.to_string())?;
        }
        Ok(Self {
            base_dir,
            worker_id,
            pairs: vec![Vec::new(); RADIX_COUNT],
            bounded: None,
            legacy_buffers: RefCell::new(None),
            failure: Arc::new(Mutex::new(None)),
            errors: Default::default(),
        })
    }

    pub fn new_bounded(base_dir: impl AsRef<Path>, worker_id: i64) -> Result<Self, String> {
        Self::new_bounded_owned(
            base_dir,
            worker_id,
            crate::util::io::IoErrorOwner::default(),
        )
    }

    pub fn new_bounded_owned(
        base_dir: impl AsRef<Path>,
        worker_id: i64,
        errors: crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        Ok(Self {
            base_dir: base_dir.as_ref().into(),
            worker_id,
            pairs: vec![Vec::new(); RADIX_COUNT],
            bounded: Some(Arc::new(Mutex::new(BucketWriter::new_file_array_owned(
                base_dir.as_ref(),
                worker_id,
                RADIX_COUNT,
                errors.clone(),
            )?))),
            legacy_buffers: RefCell::new(None),
            failure: Arc::new(Mutex::new(None)),
            errors,
        })
    }

    pub fn write_msb(&mut self, entry: PairEntry) {
        let radix = (hash64(entry.rep_oid) >> (64 - RADIX_BITS)) as usize;
        if let Some(files) = &self.bounded {
            let buffers = self.legacy_buffers.get_mut().get_or_insert_with(|| {
                WorkerBucketBuffers::new(
                    Arc::clone(files),
                    RADIX_COUNT,
                    self.errors.clone(),
                    Arc::clone(&self.failure),
                )
            });
            if let Err(error) = buffers.write(radix, &encode_pair_entry(entry)) {
                self.failure.lock().unwrap().get_or_insert(error);
            }
        } else {
            self.pairs[radix].push(entry);
        }
    }

    fn flush_legacy(&self) -> Result<(), String> {
        let buffers = self.legacy_buffers.borrow_mut().take();
        drop(buffers);
        self.failure.lock().unwrap().clone().map_or(Ok(()), Err)
    }

    pub fn entries(&self, radix: usize) -> &[PairEntry] {
        &self.pairs[radix]
    }

    /// Source per-worker BufferArray flush, retaining the enclosing
    /// round-local FileArray owner until External::round closes it.
    fn flush_unpublished(&self) -> Result<RadixedTable, String> {
        if self.bounded.is_some() {
            self.flush_legacy()?;
            return Ok(RadixedTable(
                (0..RADIX_COUNT)
                    .map(|radix| {
                        Bucket::new(
                            self.base_dir.join(radix.to_string()).join("bucket.tsv"),
                            None,
                        )
                    })
                    .collect(),
            ));
        }
        // Preserve the older explicit in-memory writer's library behavior.
        self.finish()
    }

    pub fn finish(&self) -> Result<RadixedTable, String> {
        if let Some(writer) = &self.bounded {
            self.flush_legacy()?;
            return writer.lock().unwrap().finish();
        }
        let mut buckets = Vec::with_capacity(RADIX_COUNT);
        for radix in 0..RADIX_COUNT {
            let directory = self.base_dir.join(radix.to_string());
            let manifest = directory.join("bucket.tsv");
            let entries = &self.pairs[radix];
            if entries.is_empty() {
                fs::write(source_c_path(&manifest).as_ref(), [])
                    .map_err(|error| error.to_string())?;
            } else {
                let volume = directory.join(format!("worker_{}_volume_0", self.worker_id));
                let mut compressed = CompressedBuffer::new();
                for entry in entries {
                    compressed
                        .write(&encode_pair_entry(*entry))
                        .map_err(|error| error.to_string())?;
                }
                compressed.finish().map_err(|error| error.to_string())?;
                fs::write(source_c_path(&volume).as_ref(), compressed.data())
                    .map_err(|error| error.to_string())?;
                fs::write(
                    source_c_path(&manifest).as_ref(),
                    super::storage::manifest_line(&volume, &[entries.len() as u64])
                        .map_err(|e| e.to_string())?,
                )
                .map_err(|error| error.to_string())?;
            }
            buckets.push(Bucket::new(manifest, None));
        }
        Ok(RadixedTable(buckets))
    }
}

fn encode_pair_entry(entry: PairEntry) -> [u8; 24] {
    let mut bytes = [0_u8; 24];
    bytes[0..8].copy_from_slice(&entry.rep_oid.to_ne_bytes());
    bytes[8..16].copy_from_slice(&entry.member_oid.to_ne_bytes());
    bytes[16..20].copy_from_slice(&entry.rep_len.to_ne_bytes());
    bytes[20..24].copy_from_slice(&entry.member_len.to_ne_bytes());
    bytes
}

pub fn get_pairs_uni_cov(entries: &[SeedEntry], output: &mut PairFileArray) {
    let _ = pairs_uni_cov(entries, &mut |entry| {
        output.write_msb(entry);
        Ok(())
    });
}
pub fn get_pairs_mutual_cov(
    entries: &[SeedEntry],
    min_length_ratio: f64,
    output: &mut PairFileArray,
) {
    let _ = pairs_mutual_cov(entries, min_length_ratio, &mut |entry| {
        output.write_msb(entry);
        Ok(())
    });
}

fn pairs_uni_cov(
    entries: &[SeedEntry],
    emit: &mut dyn FnMut(PairEntry) -> Result<(), String>,
) -> Result<(), String> {
    let Some(rep) = entries.first() else {
        return Ok(());
    };
    for member in &entries[1..] {
        if rep.oid != member.oid {
            emit(PairEntry {
                rep_oid: rep.oid as u64,
                member_oid: member.oid as u64,
                rep_len: rep.len as u32,
                member_len: member.len as u32,
            })?;
        }
    }
    Ok(())
}

fn pairs_mutual_cov(
    entries: &[SeedEntry],
    min_length_ratio: f64,
    emit: &mut dyn FnMut(PairEntry) -> Result<(), String>,
) -> Result<(), String> {
    let size = entries.len();
    let mut i = 0;
    let mut j = 0;
    while i < size {
        let query_len = entries[i].len;
        let mut j1 = j;
        while j1 < size {
            let target_len = entries[j1].len;
            if f64::from(target_len) / f64::from(query_len) < min_length_ratio {
                break;
            }
            j1 += 1;
        }
        let query_position = i + (j1 - j) / 2;
        let rep = entries[query_position];
        for member in &entries[j..j1] {
            if rep.oid != member.oid {
                emit(PairEntry {
                    rep_oid: rep.oid as u64,
                    member_oid: member.oid as u64,
                    rep_len: rep.len as u32,
                    member_len: member.len as u32,
                })?;
            }
        }
        j = j1;
        if j == size {
            break;
        }
        let target_len = entries[j].len;
        while i < size {
            let query_len = entries[i].len;
            if f64::from(target_len) / f64::from(query_len) >= min_length_ratio {
                break;
            }
            i += 1;
        }
    }
    Ok(())
}

/// C++ `External::build_pair_table`: outer bucket workers and per-bucket
/// equal-key partitions share one SimpleThreadPool exception/stop owner.
pub fn build_pair_table(
    job: &Job,
    seed_table: &RadixedTable,
    shape: i32,
    max_oid: i64,
    output_files: &mut PairFileArray,
    config: PairTableConfig,
) -> Result<RadixedTable, String> {
    if config.threads == 0 || seed_table.is_empty() {
        return Err(
            "Source pair-table worker allocation requires nonzero threads and buckets.".into(),
        );
    }
    let max_concurrent =
        seed_table.max_buckets_owned(job.mem_limit, SeedEntry::ENCODED_SIZE, &job.io_errors())?;
    let concurrent_buckets = usize::try_from(max_concurrent)
        .unwrap_or(usize::MAX)
        .min(config.threads);
    if concurrent_buckets == 0 {
        return Err(
            "Source pair-table worker allocation requires a nonempty bucket budget.".into(),
        );
    }
    let bucket_workers = config.threads.div_ceil(concurrent_buckets);
    job.log(&format!("Building pair table. Concurrent buckets={concurrent_buckets} Workers per bucket={bucket_workers}"))?;
    let _shift = bit_length(max_oid as u64).saturating_sub(RADIX_BITS);
    let queue = Arc::new(Atomic::new(
        job.base_dir(None)
            .join(format!("seed_table_{shape}/build_pair_table_queue")),
    ));
    let processed = Arc::new(AtomicUsize::new(0));
    let pool = Arc::new(Mutex::new(SimpleThreadPool::new()));
    let shared_job = Arc::new(job.clone());
    let seed_table = Arc::new(seed_table.clone());
    let memory_pairs = Arc::new(Mutex::new(Vec::new()));
    let mut outer = Vec::with_capacity(concurrent_buckets);
    for _ in 0..concurrent_buckets {
        let pool_for_worker = Arc::clone(&pool);
        let queue = Arc::clone(&queue);
        let processed = Arc::clone(&processed);
        let job = Arc::clone(&shared_job);
        let table = Arc::clone(&seed_table);
        let files = output_files.bounded.clone();
        let failure = Arc::clone(&output_files.failure);
        let memory_pairs = Arc::clone(&memory_pairs);
        outer.push(pool.lock().unwrap().spawn(move |stop| {
            while !stop.load(ThreadOrdering::Relaxed) {
                let bucket_index =
                    usize::try_from(queue.fetch_add_one()?).map_err(|e| e.to_string())?;
                let Some(bucket) = table.get(bucket_index) else {
                    break;
                };
                let (mut entries, volume_paths) = read_seed_bucket(&job, bucket, config.threads)?;
                job.log(&format!(
                    "Building pair table. Bucket={}/{} Records={} Size={}",
                    bucket_index + 1,
                    table.len(),
                    crate::util::string::format(entries.len() as u64),
                    crate::util::string::format((entries.len() * SeedEntry::ENCODED_SIZE) as u64)
                ))?;
                if bucket_workers > 1 {
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(bucket_workers)
                        .build()
                        .map_err(|error| {
                            format!("Error creating pair-table sorting threads: {error}")
                        })?
                        .install(|| entries.par_sort_unstable());
                } else {
                    entries.sort_unstable();
                }
                let ranges = seed_part_ranges(&entries, bucket_workers);
                let entries = Arc::new(entries);
                let mut inner = Vec::with_capacity(ranges.len());
                for range in ranges {
                    let entries = Arc::clone(&entries);
                    let files = files.clone();
                    let failure = Arc::clone(&failure);
                    let errors = job.io_errors();
                    let memory_pairs = Arc::clone(&memory_pairs);
                    inner.push(pool_for_worker.lock().unwrap().spawn(move |stop| {
                        // A Source BufferArray belongs to every spawned input
                        // part, even when equal-key extension leaves it empty.
                        let mut buffers = files.map(|files| {
                            WorkerBucketBuffers::new(files, RADIX_COUNT, errors, failure)
                        });
                        let mut local_pairs = Vec::new();
                        let mut emit = |entry: PairEntry| {
                            if let Some(buffers) = &mut buffers {
                                let radix = (hash64(entry.rep_oid) >> (64 - RADIX_BITS)) as usize;
                                buffers.write(radix, &encode_pair_entry(entry))
                            } else {
                                local_pairs.push(entry);
                                Ok(())
                            }
                        };
                        let mut begin = range.start;
                        while !stop.load(ThreadOrdering::Relaxed) && begin < range.end {
                            let mut end = begin + 1;
                            while end < range.end
                                && entries[end].merge_key() == entries[begin].merge_key()
                            {
                                end += 1;
                            }
                            if config.mutual_cover_present {
                                pairs_mutual_cov(
                                    &entries[begin..end],
                                    config.min_length_ratio,
                                    &mut emit,
                                )?;
                            } else {
                                pairs_uni_cov(&entries[begin..end], &mut emit)?;
                            }
                            begin = end;
                        }
                        drop(emit);
                        drop(buffers); // Flush worker members before input removal.
                        memory_pairs.lock().unwrap().extend(local_pairs);
                        Ok(())
                    }));
                }
                SimpleThreadPool::join_shared(&pool_for_worker, inner)?;
                // VolumedFile::remove intentionally ignores remove/rmdir
                // results, after all input handles and bucket workers close.
                for path in volume_paths {
                    let _ = fs::remove_file(source_c_path(&path).as_ref());
                }
                let _ = fs::remove_file(source_c_path(&bucket.path).as_ref());
                if let Some(parent) = bucket.path.parent() {
                    let _ = fs::remove_dir(source_c_path(parent).as_ref());
                }
                processed.fetch_add(1, ThreadOrdering::Relaxed);
            }
            Ok(())
        }));
    }
    SimpleThreadPool::join_shared(&pool, outer)?;
    if let Some(error) = output_files.failure.lock().unwrap().clone() {
        return Err(error);
    }
    // The explicit unbounded library adapter collects logical records only;
    // actual CLI workers wrote through the retained round FileArray above.
    if output_files.bounded.is_none() {
        for entry in memory_pairs.lock().unwrap().drain(..) {
            output_files.write_msb(entry);
        }
    }
    let table = output_files.flush_unpublished()?;
    let finished = Atomic::new(
        job.base_dir(None)
            .join(format!("seed_table_{shape}/pair_table_finished")),
    );
    finished.fetch_add(processed.load(ThreadOrdering::Relaxed) as i64)?;
    finished.await_value(seed_table.len() as i64)?;
    Ok(table)
}

fn seed_part_ranges(entries: &[SeedEntry], workers: usize) -> Vec<std::ops::Range<usize>> {
    let partition = Partition::new(entries.len(), workers);
    (0..partition.parts)
        .map(|part| {
            let mut begin = partition.begin(part);
            let mut end = partition.end(part);
            if part > 0 {
                while begin < entries.len() && entries[begin - 1].key() == entries[begin].key() {
                    begin += 1;
                }
            }
            while end < entries.len() && entries[end - 1].key() == entries[end].key() {
                end += 1;
            }
            begin..end
        })
        .collect()
}

fn bit_length(value: u64) -> u32 {
    u64::BITS - value.leading_zeros()
}

/// InputBuffer<SeedEntry> reads actual disjoint volume ranges concurrently.
/// Scoped slices permit direct filling of the one Source data owner without
/// extra full-volume Vec copies or unsafe shared writes/lifetime extension.
fn read_seed_bucket(
    job: &Job,
    bucket: &Bucket,
    threads: usize,
) -> Result<(Vec<SeedEntry>, Vec<PathBuf>), String> {
    let (entries, storage) =
        super::external::read_pod_input_buffer(job, bucket, threads, 20, |_, bytes| {
            Ok(SeedEntry {
                seed: u64::from_ne_bytes(bytes[0..8].try_into().unwrap()),
                oid: i64::from_ne_bytes(bytes[8..16].try_into().unwrap()),
                len: i32::from_ne_bytes(bytes[16..20].try_into().unwrap()),
            })
        })?;
    Ok((entries, storage.into_iter().map(|(path, _)| path).collect()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    fn temp_dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond_pair_table_{name}_{}_{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, AtomicOrdering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    fn all_pairs(output: &PairFileArray) -> Vec<PairEntry> {
        let mut pairs: Vec<_> = output.pairs.iter().flatten().copied().collect();
        pairs.sort_unstable();
        pairs
    }

    #[test]
    fn seed_input_preserves_source_partial_pod_and_zero_tail() {
        let root = temp_dir("short_seed_input");
        let job = Job::new(2, 1, 1024, root.join("job"), 0).unwrap();
        let path = root.join("volume");
        let manifest = root.join("input.tsv");
        fs::write(
            &manifest,
            super::super::storage::manifest_line(&path, &[3]).unwrap(),
        )
        .unwrap();
        let payload = (1..=29).collect::<Vec<u8>>();
        fs::write(&path, &payload).unwrap();
        let (entries, _) = read_seed_bucket(&job, &Bucket::new(&manifest, None), 2).unwrap();
        let mut actual = Vec::new();
        for entry in entries {
            actual.extend_from_slice(&entry.seed.to_ne_bytes());
            actual.extend_from_slice(&entry.oid.to_ne_bytes());
            actual.extend_from_slice(&entry.len.to_ne_bytes());
        }
        let mut expected = vec![0; 60];
        expected[..payload.len()].copy_from_slice(&payload);
        assert_eq!(actual, expected);
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bounded_pair_worker_flush_does_not_publish_round_owned_volumes() {
        let root = temp_dir("round_owner_publication");
        let mut output = PairFileArray::new_bounded(&root, 7).unwrap();
        let entry = PairEntry {
            rep_oid: 4,
            member_oid: 9,
            rep_len: 100,
            member_len: 90,
        };
        output.write_msb(entry);
        let radix = (hash64(entry.rep_oid) >> (64 - RADIX_BITS)) as usize;
        let handles = output.flush_unpublished().unwrap();
        let manifest = &handles[radix].path;
        assert!(fs::read(manifest).unwrap().is_empty());
        let volume = root.join(radix.to_string()).join("worker_7_volume_0");
        // Source OutputFile retains these small compressed bytes until the
        // enclosing FileArray closes, after the worker completion barrier.
        assert_eq!(volume.metadata().unwrap().len(), 0);
        output.finish().unwrap();
        let rows = super::super::storage::manifest_records(manifest).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].1, 1);
        let mut input = InputFile::source_new(volume.to_str().unwrap(), 0).unwrap();
        let mut decoded = [0; 24];
        assert_eq!(input.read_raw(&mut decoded).unwrap(), decoded.len());
        assert_eq!(decoded.as_slice(), encode_pair_entry(entry).as_slice());
        input.close().unwrap();
        output.finish().unwrap();
        assert_eq!(
            super::super::storage::manifest_records(manifest).unwrap(),
            rows
        );
        drop(output);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn seed_layout_and_order_match_packed_header_type() {
        let entry = SeedEntry::new(0x0102_0304_0506_0708, -9, 1234);
        let mut bytes = Vec::new();
        entry.serialize(&mut bytes);
        assert_eq!(bytes.len(), SeedEntry::ENCODED_SIZE);
        assert_eq!(SeedEntry::deserialize(&bytes).unwrap(), entry);

        let mut entries = [
            SeedEntry::new(7, 2, 50),
            SeedEntry::new(7, 1, 100),
            SeedEntry::new(7, 0, 100),
            SeedEntry::new(6, 9, 1),
        ];
        entries.sort_unstable();
        assert_eq!(
            entries.iter().map(|entry| entry.oid).collect::<Vec<_>>(),
            [9, 0, 1, 2]
        );
    }

    #[test]
    fn unidirectional_uses_longest_lowest_oid_rep_and_skips_self() {
        let root = temp_dir("uni");
        let mut output = PairFileArray::new(&root, 0).unwrap();
        let entries = [
            SeedEntry::new(1, 4, 100),
            SeedEntry::new(1, 4, 90),
            SeedEntry::new(1, 8, 80),
        ];
        get_pairs_uni_cov(&entries, &mut output);
        assert_eq!(
            all_pairs(&output),
            [PairEntry {
                rep_oid: 4,
                member_oid: 8,
                rep_len: 100,
                member_len: 80,
            }]
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn mutual_cover_windows_choose_middle_representatives() {
        let root = temp_dir("mutual");
        let mut output = PairFileArray::new(&root, 0).unwrap();
        let entries = [
            SeedEntry::new(1, 0, 100),
            SeedEntry::new(1, 1, 90),
            SeedEntry::new(1, 2, 80),
            SeedEntry::new(1, 3, 40),
        ];
        get_pairs_mutual_cov(&entries, 0.8, &mut output);
        assert_eq!(
            all_pairs(&output),
            [
                PairEntry {
                    rep_oid: 1,
                    member_oid: 0,
                    rep_len: 90,
                    member_len: 100,
                },
                PairEntry {
                    rep_oid: 1,
                    member_oid: 2,
                    rep_len: 90,
                    member_len: 80,
                },
            ]
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn build_reads_compressed_seed_volume_routes_pairs_and_removes_input() {
        let root = temp_dir("build");
        let seed_dir = root.join("seed_bucket");
        fs::create_dir_all(&seed_dir).unwrap();
        let volume = seed_dir.join("volume0");
        let manifest = seed_dir.join("bucket.tsv");
        let seeds = [
            SeedEntry::new(11, 2, 80),
            SeedEntry::new(10, 3, 70),
            SeedEntry::new(10, 1, 90),
            SeedEntry::new(11, 0, 100),
        ];
        let mut compressed = CompressedBuffer::new();
        for seed in seeds {
            let mut bytes = Vec::new();
            seed.serialize(&mut bytes);
            compressed.write(&bytes).unwrap();
        }
        compressed.finish().unwrap();
        fs::write(&volume, compressed.data()).unwrap();
        fs::write(
            &manifest,
            format!("{}\t{}\n", volume.display(), seeds.len()),
        )
        .unwrap();

        let job = Job::new(3, 1, 1 << 20, root.join("job"), 5).unwrap();
        // The real seed builder creates this stage directory before the pair
        // queue is opened, even when radix manifests live in child buckets.
        fs::create_dir_all(job.base_dir(None).join("seed_table_0")).unwrap();
        let table = RadixedTable(vec![Bucket::new(&manifest, Some(seeds.len() as u64))]);
        let pair_dir = root.join("pairs");
        let mut output = PairFileArray::new(&pair_dir, job.worker_id()).unwrap();
        let result =
            build_pair_table(&job, &table, 0, 3, &mut output, PairTableConfig::default()).unwrap();

        assert_eq!(result.len(), RADIX_COUNT);
        assert_eq!(
            all_pairs(&output),
            [
                PairEntry {
                    rep_oid: 0,
                    member_oid: 2,
                    rep_len: 100,
                    member_len: 80,
                },
                PairEntry {
                    rep_oid: 1,
                    member_oid: 3,
                    rep_len: 90,
                    member_len: 70,
                },
            ]
        );
        assert!(!volume.exists());
        assert!(!manifest.exists());
        for pair in all_pairs(&output) {
            let radix = (hash64(pair.rep_oid) >> 56) as usize;
            assert!(output.entries(radix).contains(&pair));
            let line = fs::read_to_string(&result[radix].path).unwrap();
            assert!(line.ends_with("\t1\n"));
            let pair_volume = line.split_whitespace().next().unwrap();
            let mut input = InputFile::new(pair_volume, 0).unwrap();
            assert_eq!(input.read_u64().unwrap(), pair.rep_oid);
            assert_eq!(input.read_u64().unwrap(), pair.member_oid);
            assert_eq!(input.read_u32().unwrap(), pair.rep_len);
            assert_eq!(input.read_u32().unwrap(), pair.member_len);
        }
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn source_partition_keeps_equal_seed_groups_in_one_part() {
        let entries: Vec<_> = (0..8)
            .map(|oid| SeedEntry::new(if oid < 6 { 10 } else { 20 }, oid, 100))
            .collect();
        assert_eq!(seed_part_ranges(&entries, 4), vec![0..6, 6..6, 6..6, 6..8]);
        assert_eq!(
            seed_part_ranges(&[], 4),
            Vec::<std::ops::Range<usize>>::new()
        );
    }

    fn write_seed_bucket_fixture(root: &Path, number: usize) -> Bucket {
        let directory = root.join(format!("seed-bucket-{number}"));
        fs::create_dir_all(&directory).unwrap();
        let mut manifest = Vec::new();
        for part in 0..2 {
            let path = directory.join(format!("volume-{part}"));
            let mut compressed = CompressedBuffer::new();
            for offset in (part * 4..part * 4 + 4).rev() {
                let entry = SeedEntry::new(number as u64 + 10, (number * 10 + offset) as i64, 100);
                let mut bytes = Vec::new();
                entry.serialize(&mut bytes);
                compressed.write(&bytes).unwrap();
            }
            compressed.finish().unwrap();
            fs::write(&path, compressed.data()).unwrap();
            manifest.extend(super::super::storage::manifest_line(&path, &[4]).unwrap());
        }
        let path = directory.join("bucket.tsv");
        fs::write(&path, manifest).unwrap();
        Bucket::new(path, Some(8))
    }

    #[test]
    fn bounded_pair_table_runs_both_worker_levels_without_early_publication() {
        // One memory-limited bucket exercises four inner parts; four buckets
        // exercise four concurrent outer claim loops with one inner part each.
        for (bucket_count, mem_limit) in [(1, 1), (4, 1 << 20)] {
            let root = temp_dir("parallel_pair_owners");
            let table = RadixedTable(
                (0..bucket_count)
                    .map(|i| write_seed_bucket_fixture(&root, i))
                    .collect(),
            );
            let paths: Vec<_> = table.iter().map(|bucket| bucket.path.clone()).collect();
            let job = Job::new(40, 2, mem_limit, root.join("job"), 2).unwrap();
            fs::create_dir_all(job.base_dir(None).join("seed_table_0")).unwrap();
            fs::create_dir(root.join("pairs")).unwrap();
            let mut output = PairFileArray::new_bounded(root.join("pairs"), 2).unwrap();
            let result = build_pair_table(
                &job,
                &table,
                0,
                40,
                &mut output,
                PairTableConfig {
                    threads: 4,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(paths.iter().all(|path| !path.exists()));
            assert!(
                result
                    .iter()
                    .all(|bucket| fs::read(&bucket.path).unwrap().is_empty()),
                "worker BufferArrays flush but the round FileArray has not published"
            );
            output.finish().unwrap();
            let mut actual = Vec::new();
            for bucket in result.iter() {
                for (path, count) in super::super::storage::manifest_records(&bucket.path).unwrap()
                {
                    let mut input = InputFile::with_source_profile_path(&path, 0, 65_536).unwrap();
                    for _ in 0..count {
                        actual.push(PairEntry {
                            rep_oid: input.read_u64().unwrap(),
                            member_oid: input.read_u64().unwrap(),
                            rep_len: input.read_u32().unwrap(),
                            member_len: input.read_u32().unwrap(),
                        });
                    }
                    input.close().unwrap();
                }
            }
            actual.sort_unstable();
            let mut expected: Vec<_> = (0..bucket_count)
                .flat_map(|bucket| {
                    (1..8).map(move |member| PairEntry {
                        rep_oid: (bucket * 10) as u64,
                        member_oid: (bucket * 10 + member) as u64,
                        rep_len: 100,
                        member_len: 100,
                    })
                })
                .collect();
            expected.sort_unstable();
            assert_eq!(actual, expected);
            drop(output);
            drop(job);
            fs::remove_dir_all(root).unwrap();
        }
    }

    #[test]
    fn max_buckets_matches_memory_cutoff_and_missing_counts_error() {
        let table = RadixedTable(vec![
            Bucket::new("a", Some(10)),
            Bucket::new("b", Some(5)),
            Bucket::new("c", Some(1)),
        ]);
        assert_eq!(table.max_buckets(250, 20).unwrap(), 1);
        assert_eq!(table.max_buckets(201, 20).unwrap(), 1);
        assert_eq!(table.max_buckets(301, 20).unwrap(), 2);
        assert!(RadixedTable(vec![Bucket::new("x", None)])
            .max_buckets(1, 20)
            .unwrap_err()
            .contains("Record count not set"));
    }
}
