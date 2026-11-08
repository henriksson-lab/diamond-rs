//! Bounded FileArray/BufferArray storage and file-backed radix work queues.
//! Records retain the Source native-byte packed layout. Each flushed buffer
//! is a complete codec member, allowing a volume to be consumed incrementally.

use crate::util::io::source_c_path;
use std::fs;
use std::path::{Path, PathBuf};

use super::external::{base_path, ExternalRecordKind, Job};
use super::pair_table::{Bucket, RadixedTable, RADIX_BITS, RADIX_COUNT};
use crate::util::hash::hash64;
use crate::util::io::{CompressedBuffer, Compressor, IoErrorOwner, OutputFile};
use crate::util::parallel::{Atomic, FileStack};

/// Source manifests contain the path bytes verbatim; whitespace is not escaped.
pub(crate) fn manifest_line(path: &Path, counts: &[u64]) -> crate::util::io::IoResult<Vec<u8>> {
    let mut row = crate::util::io::source_os_bytes(path.as_os_str())?.into_owned();
    for count in counts {
        row.push(b'\t');
        row.extend_from_slice(count.to_string().as_bytes());
    }
    row.push(b'\n');
    Ok(row)
}

pub(crate) fn path_error(errors: &IoErrorOwner, prefix: &[u8], path: &Path) -> String {
    errors
        .capture::<()>(Err(std::io::Error::other(
            crate::util::io::IoError::path_message(prefix, path),
        )))
        .unwrap_err()
}

// FileArray concatenates separators as std::string bytes. Path::join would
// suppress the extra separator when Source's seed-table base already ends in
// one, changing stored manifest names even though both paths open the same file.
pub(crate) fn source_child_path(base: &Path, name: &str) -> PathBuf {
    let mut path = base.as_os_str().to_owned();
    path.push(std::path::MAIN_SEPARATOR.to_string());
    path.push(name);
    PathBuf::from(path)
}
fn file_array_directory(base: &Path, radix: usize) -> PathBuf {
    let mut path = source_child_path(base, &radix.to_string()).into_os_string();
    path.push(std::path::MAIN_SEPARATOR.to_string());
    PathBuf::from(path)
}
fn file_array_path(base: &Path, radix: usize, name: &str) -> PathBuf {
    let mut path = file_array_directory(base, radix).into_os_string();
    path.push(name);
    PathBuf::from(path)
}

const BUFFER_BYTES: usize = 65_536;
const MAX_VOLUME_BYTES: u64 = 1 << 30;

#[derive(Debug, Default)]
struct BucketOutput {
    compressed: Option<CompressedBuffer>,
    buffer_active: bool,
    pending_records: u64,
    file: Option<OutputFile>,
    stack: Option<FileStack>,
    volume: u64,
    volume_bytes: u64,
    volume_records: u64,
    records_total: u64,
}

#[derive(Debug)]
pub(crate) struct BucketWriter {
    directory: PathBuf,
    worker: i64,
    exclusive: bool,
    buckets: Vec<BucketOutput>,
    error: Option<String>,
    errors: IoErrorOwner,
}

impl BucketWriter {
    /// Only actual Source command owners opt into its noexcept destructors.
    /// Keep the primary exception slot separate from cleanup operations so
    /// their saved errno and emitted-perror marker cannot overwrite each other.
    fn source_cleanup<T>(
        &mut self,
        operation: impl FnOnce(&mut Self) -> Result<T, String>,
    ) -> Result<T, String> {
        if !self.errors.source_destructors() {
            return operation(self);
        }
        let primary = self.errors.clone();
        if let Err(error) = primary.emit_pending_perror() {
            crate::util::io::abort_source_diagnostic_error(&error);
        }
        let deferred = self.error.take();
        self.errors = primary.cleanup_owner();
        let result = operation(self);
        let cleanup = std::mem::replace(&mut self.errors, primary);
        match result {
            Err(message) => {
                let error = cleanup.finish::<()>(Err(message)).unwrap_err();
                if let Err(error) = crate::util::io::emit_source_perror(&error) {
                    crate::util::io::abort_source_diagnostic_error(&error);
                }
                // Source main installs a terminate handler which only aborts;
                // it never reaches its ordinary what/log exception catch.
                crate::util::io::abort_source_process();
            }
            Ok(value) => {
                if let Some(message) = deferred {
                    self.error = Some(message.clone());
                    Err(message)
                } else {
                    Ok(value)
                }
            }
        }
    }

    pub(crate) fn new(directory: &Path, worker: i64, count: usize) -> Result<Self, String> {
        Self::new_owned(directory, worker, count, IoErrorOwner::default())
    }

    pub(crate) fn new_owned(
        directory: &Path,
        worker: i64,
        count: usize,
        errors: IoErrorOwner,
    ) -> Result<Self, String> {
        Self::new_storage_owned(directory, worker, count, errors, true)
    }

    fn new_storage_owned(
        directory: &Path,
        worker: i64,
        count: usize,
        errors: IoErrorOwner,
        local_buffers: bool,
    ) -> Result<Self, String> {
        let mut buckets = Vec::with_capacity(count);
        for radix in 0..count {
            let dir = file_array_directory(directory, radix);
            fs::create_dir_all(source_c_path(&dir).as_ref()).map_err(|error| error.to_string())?;
            // FileArray opens every initial volume before its FileStack,
            // including buckets which will receive no records.
            let file = Self::open_volume(
                &file_array_path(directory, radix, &format!("worker_{worker}_volume_0")),
                &errors,
            )?;
            let stack = FileStack::try_new_owned(
                file_array_path(directory, radix, "bucket.tsv"),
                errors.clone(),
            )?;
            buckets.push(BucketOutput {
                file: Some(file),
                stack: Some(stack),
                compressed: local_buffers.then(CompressedBuffer::default),
                buffer_active: local_buffers,
                ..BucketOutput::default()
            });
        }
        Ok(Self {
            directory: directory.into(),
            worker,
            exclusive: false,
            buckets,
            error: None,
            errors,
        })
    }

    /// FileArray owner without a caller-local BufferArray. Compression belongs
    /// to the worker owners below; this owner retains only files and manifests.
    pub(crate) fn new_file_array_owned(
        directory: &Path,
        worker: i64,
        count: usize,
        errors: IoErrorOwner,
    ) -> Result<Self, String> {
        Self::new_storage_owned(directory, worker, count, errors, false)
    }

    pub(crate) fn new_exclusive_file_array_owned(
        directory: &Path,
        worker: i64,
        count: usize,
        errors: IoErrorOwner,
    ) -> Result<Self, String> {
        let mut owner = Self::new_file_array_owned(directory, worker, count, errors)?;
        owner.exclusive = true;
        Ok(owner)
    }

    pub(crate) fn buckets(&self) -> RadixedTable {
        RadixedTable(
            self.buckets
                .iter()
                .enumerate()
                .map(|(radix, bucket)| {
                    Bucket::new(
                        file_array_path(&self.directory, radix, "bucket.tsv"),
                        self.exclusive.then_some(bucket.records_total),
                    )
                })
                .collect(),
        )
    }

    fn append_worker_member(
        &mut self,
        radix: usize,
        bytes: &[u8],
        records: u64,
        errors: &IoErrorOwner,
    ) -> Result<(), String> {
        // A worker's noexcept destructor uses its cleanup exception slot; the
        // enclosing file owner must not overwrite an already pending primary.
        let previous = std::mem::replace(&mut self.errors, errors.clone());
        let result = self.append_member(radix, bytes, records);
        self.errors = previous;
        result
    }

    fn append_member(&mut self, radix: usize, bytes: &[u8], records: u64) -> Result<(), String> {
        let path = self.volume_path(radix);
        let bucket = self.buckets.get_mut(radix).ok_or("Invalid output radix")?;
        if bucket.file.is_none() {
            bucket.file = Some(Self::open_volume(&path, &self.errors)?);
        }
        self.errors.capture(
            bucket
                .file
                .as_mut()
                .unwrap()
                .write_raw(bytes)
                .map_err(std::io::Error::other),
        )?;
        bucket.volume_bytes += bytes.len() as u64;
        bucket.volume_records += records;
        if self.exclusive {
            bucket.records_total += records;
        }
        if bucket.volume_bytes >= MAX_VOLUME_BYTES {
            self.publish(radix, true)?;
            self.buckets[radix].file =
                Some(Self::open_volume(&self.volume_path(radix), &self.errors)?);
        }
        Ok(())
    }

    /// Preserve the older infallible callback API while retaining the first
    /// actual I/O failure for mandatory propagation at flush/finish.
    pub(crate) fn write_deferred(&mut self, radix: usize, bytes: &[u8]) {
        if self.error.is_none() {
            if let Err(error) = self.write(radix, bytes) {
                self.error = Some(error);
            }
        }
    }

    pub(crate) fn write(&mut self, radix: usize, bytes: &[u8]) -> Result<(), String> {
        self.write_count(radix, bytes, 1)
    }

    pub(crate) fn write_count(
        &mut self,
        radix: usize,
        bytes: &[u8],
        records: u64,
    ) -> Result<(), String> {
        if let Some(error) = &self.error {
            return Err(error.clone());
        }
        let bucket = self
            .buckets
            .get_mut(radix)
            .ok_or_else(|| "Invalid output radix".to_owned())?;
        self.errors.capture(
            bucket
                .compressed
                .as_mut()
                .expect("local BufferArray owner")
                .write(bytes)
                .map_err(std::io::Error::other),
        )?;
        bucket.buffer_active = true;
        bucket.pending_records += records;
        // BufferArray tests CompressedBuffer::size(), not the input size.
        if bucket
            .compressed
            .as_ref()
            .expect("local BufferArray owner")
            .size()
            >= BUFFER_BYTES
        {
            self.flush(radix)?;
        }
        Ok(())
    }

    fn volume_path(&self, radix: usize) -> PathBuf {
        file_array_path(
            &self.directory,
            radix,
            &format!(
                "worker_{}_volume_{}",
                self.worker, self.buckets[radix].volume
            ),
        )
    }

    fn open_volume(path: &Path, errors: &IoErrorOwner) -> Result<OutputFile, String> {
        // External::external fixes config.file_buffer_size to 64 KiB before
        // constructing any FileArray, regardless of the CLI's initial value.
        errors.capture(
            OutputFile::with_source_profile_path(path, Compressor::None, "wb", BUFFER_BYTES)
                .map_err(std::io::Error::other),
        )
    }

    fn flush(&mut self, radix: usize) -> Result<(), String> {
        if !self.buckets[radix].buffer_active {
            return Ok(());
        }
        let path = self.volume_path(radix);
        let bucket = &mut self.buckets[radix];
        if bucket.file.is_none() {
            bucket.file = Some(Self::open_volume(&path, &self.errors)?);
        }
        self.errors.capture(
            bucket
                .compressed
                .as_mut()
                .expect("local BufferArray owner")
                .finish()
                .map_err(std::io::Error::other),
        )?;
        self.errors.capture(
            bucket
                .file
                .as_mut()
                .unwrap()
                .write_raw(
                    bucket
                        .compressed
                        .as_ref()
                        .expect("local BufferArray owner")
                        .data(),
                )
                .map_err(std::io::Error::other),
        )?;
        bucket.volume_bytes += bucket
            .compressed
            .as_ref()
            .expect("local BufferArray owner")
            .size() as u64;
        bucket.volume_records += bucket.pending_records;
        if self.exclusive {
            bucket.records_total += bucket.pending_records;
        }
        bucket.pending_records = 0;
        bucket
            .compressed
            .as_mut()
            .expect("local BufferArray owner")
            .clear();
        if bucket.volume_bytes >= MAX_VOLUME_BYTES {
            self.publish(radix, true)?;
            // FileArray::write replaces a full volume immediately, even if
            // no subsequent record arrives before the owner closes.
            self.buckets[radix].file =
                Some(Self::open_volume(&self.volume_path(radix), &self.errors)?);
        }
        Ok(())
    }

    fn publish(&mut self, radix: usize, rollover: bool) -> Result<(), String> {
        let path = self.volume_path(radix);
        let bucket = &mut self.buckets[radix];
        if let Some(mut file) = bucket.file.take() {
            // FileArray::write publishes a full volume before closing it;
            // its destructor closes the last volume before publication.
            if !rollover {
                self.errors
                    .capture(file.close().map_err(std::io::Error::other))?;
            }
            if rollover || bucket.volume_records > 0 {
                // A rollover publishes even a zero-record volume; only the
                // final destructor removes an empty current file. FileStack
                // takes the Source cross-process lock around each line.
                bucket
                    .stack
                    .as_ref()
                    .unwrap()
                    .push_bytes(
                        &self.errors.capture(
                            manifest_line(&path, &[bucket.volume_records])
                                .map_err(std::io::Error::other),
                        )?,
                    )?;
            } else {
                // Source FileArray deliberately ignores remove(3)'s result.
                let _ = fs::remove_file(source_c_path(&path).as_ref());
            }
            if rollover {
                self.errors
                    .capture(file.close().map_err(std::io::Error::other))?;
            }
            bucket.volume += 1;
            bucket.volume_records = 0;
            bucket.volume_bytes = 0;
        }
        Ok(())
    }

    /// Flush and publish newly written volumes. This is repeatable: pair
    /// tables accumulate all seed shapes through the same file owner.
    pub(crate) fn finish(&mut self) -> Result<RadixedTable, String> {
        self.source_cleanup(Self::finish_inner)
    }

    fn finish_inner(&mut self) -> Result<RadixedTable, String> {
        if let Some(error) = &self.error {
            return Err(error.clone());
        }
        let mut buckets = Vec::with_capacity(self.buckets.len());
        for radix in 0..self.buckets.len() {
            self.flush(radix)?;
            self.buckets[radix].buffer_active = false;
            self.publish(radix, false)?;
            let manifest = file_array_path(&self.directory, radix, "bucket.tsv");
            // Other process workers can still append. Read aggregate counts
            // only after the stage's shared completion barrier.
            buckets.push(Bucket::new(
                manifest,
                self.exclusive.then_some(self.buckets[radix].records_total),
            ));
        }
        Ok(RadixedTable(buckets))
    }

    /// Source BufferArray destruction flushes worker records before the
    /// completion barrier, while the enclosing FileArray still owns the
    /// current volumes. Only its later destruction publishes those volumes.
    /// A full 1-GiB volume retains the normal rollover publication in flush().
    pub(crate) fn flush_unpublished(&mut self) -> Result<RadixedTable, String> {
        self.source_cleanup(Self::flush_unpublished_inner)
    }

    fn flush_unpublished_inner(&mut self) -> Result<RadixedTable, String> {
        if let Some(error) = &self.error {
            return Err(error.clone());
        }
        let mut buckets = Vec::with_capacity(self.buckets.len());
        for radix in 0..self.buckets.len() {
            // Each worker BufferArray destructor finishes even an empty
            // stream. FileArray destruction itself adds no extra member.
            if self.buckets[radix].compressed.is_some() {
                self.buckets[radix].buffer_active = true;
                self.flush(radix)?;
                self.buckets[radix].buffer_active = false;
            }
            buckets.push(Bucket::new(
                file_array_path(&self.directory, radix, "bucket.tsv"),
                None,
            ));
        }
        Ok(RadixedTable(buckets))
    }
}

impl Drop for BucketWriter {
    fn drop(&mut self) {
        // BufferArray finishes its streams before FileArray closes and
        // publishes/removes its volumes. A partial constructor never creates
        // this owner and retains its ordinary field-drop behavior.
        if self.errors.source_destructors() {
            let _ = self.source_cleanup(|writer| {
                for radix in 0..writer.buckets.len() {
                    writer.flush(radix)?;
                    writer.buckets[radix].buffer_active = false;
                    writer.publish(radix, false)?;
                }
                Ok(())
            });
            return;
        }
        // Reusable library extensions remain fallible/best effort rather
        // than inheriting the command's process-level terminate policy.
        self.errors = self.errors.cleanup_owner();
        for radix in 0..self.buckets.len() {
            let _ = self.flush(radix);
            self.buckets[radix].buffer_active = false;
            let _ = self.publish(radix, false);
        }
    }
}

pub(crate) fn manifest_records(path: &Path) -> Result<Vec<(PathBuf, u64)>, String> {
    manifest_records_owned(path, &IoErrorOwner::default())
}

pub(crate) fn manifest_records_owned(
    path: &Path,
    errors: &IoErrorOwner,
) -> Result<Vec<(PathBuf, u64)>, String> {
    Ok(
        super::external::read_volumed_file_manifest_owned(path, errors)?
            .0
            .into_iter()
            .map(|volume| (volume.path, volume.record_count))
            .collect(),
    )
}

fn record_size(kind: ExternalRecordKind) -> usize {
    match kind {
        ExternalRecordKind::Seed => 20,
        ExternalRecordKind::ChunkTable => 12,
        _ => 24,
    }
}

fn key(bytes: &[u8], kind: ExternalRecordKind) -> u64 {
    let offset = if kind == ExternalRecordKind::Edge {
        8
    } else {
        0
    };
    let value = u64::from_ne_bytes(bytes[offset..offset + 8].try_into().unwrap());
    if kind == ExternalRecordKind::ChunkTable {
        value
    } else {
        hash64(value)
    }
}

fn spill(
    job: &Job,
    bucket: Bucket,
    kind: ExternalRecordKind,
    bits: u32,
    volumes: Vec<(PathBuf, u64)>,
    count: u64,
    threads: usize,
) -> Result<Vec<Bucket>, String> {
    if count == 0 {
        // VolumedFile::remove ignores these results in Source.
        for (path, _) in volumes {
            let _ = fs::remove_file(source_c_path(&path).as_ref());
        }
        let _ = fs::remove_file(source_c_path(&bucket.path).as_ref());
        if let Some(directory) = bucket.path.parent() {
            let _ = fs::remove_dir(source_c_path(directory).as_ref());
        }
        return Ok(Vec::new());
    }
    if bits < RADIX_BITS || count.saturating_mul(record_size(kind) as u64) <= job.mem_limit {
        return Ok(vec![Bucket::new(bucket.path, Some(count))]);
    }
    let shift = bits - RADIX_BITS;
    let output = std::sync::Arc::new(std::sync::Mutex::new(
        BucketWriter::new_exclusive_file_array_owned(
            bucket.path.parent().unwrap_or(Path::new("")),
            job.worker_id(),
            RADIX_COUNT,
            job.io_errors(),
        )?,
    ));
    let next = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let failure = std::sync::Arc::new(std::sync::Mutex::new(None));
    let volumes = std::sync::Arc::new(volumes);
    let mut pool = crate::util::parallel::SimpleThreadPool::new();
    for _ in 0..threads.min(volumes.len()) {
        let job = job.clone();
        let output = std::sync::Arc::clone(&output);
        let next = std::sync::Arc::clone(&next);
        let failure = std::sync::Arc::clone(&failure);
        let volumes = std::sync::Arc::clone(&volumes);
        pool.spawn(move |stop| {
            let mut buffers =
                WorkerBucketBuffers::new(output, RADIX_COUNT, job.io_errors(), failure);
            while !stop.load(std::sync::atomic::Ordering::Relaxed) {
                let index = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let Some((path, records)) = volumes.get(index) else {
                    break;
                };
                let mut input = job.open_input(path, 0)?;
                let mut bytes = vec![0; record_size(kind)];
                for _ in 0..*records {
                    // Source radix_sort.h:38 passes sizeof(T) as a typed
                    // element count (Deserializer multiplies it by sizeof(T)).
                    // Requests can overflow its one T object when more than T
                    // bytes remain. Read safely within one record instead; the
                    // Source-undefined multi-record case gets defined Rust IO.
                    // Source ignores the count: partial/EOF reads leave the
                    // initialized prior object bytes intact within this volume.
                    let _ = job.io(input.read_raw(&mut bytes))?;
                    buffers.write(((key(&bytes, kind) >> shift) & 255) as usize, &bytes)?;
                }
                job.io(input.close())?;
            }
            // BufferArray flushes at worker destruction, before pool.join_all
            // returns and before the exclusive FileArray publishes its files.
            Ok(())
        });
    }
    pool.join_all()?;
    if let Some(error) = failure.lock().unwrap().clone() {
        return Err(error);
    }
    // Source reads exclusive cumulative counters directly from the real file
    // owner BEFORE closing it. No manifest reopen or post-close recount occurs.
    let mut timer = crate::util::log_stream::TaskTimer::with_msg("Closing the output files", 1);
    let children = output.lock().unwrap().buckets();
    output.lock().unwrap().finish()?;
    drop(output);
    timer.finish();
    job.log(&format!("Radix sorted bucket records={count}"))?;
    // One eight-bit pass only. Source retains the old input bucket here and
    // does not recursively split a skewed child which still exceeds its limit.
    Ok(children.0)
}

pub fn radix_sort(
    job: &Job,
    table: RadixedTable,
    kind: ExternalRecordKind,
    bits_unsorted: u32,
) -> Result<RadixedTable, String> {
    radix_sort_with_threads(job, table, kind, bits_unsorted, 1)
}

pub fn radix_sort_with_threads(
    job: &Job,
    table: RadixedTable,
    kind: ExternalRecordKind,
    bits_unsorted: u32,
    threads: usize,
) -> Result<RadixedTable, String> {
    if bits_unsorted < RADIX_BITS || table.is_empty() {
        return Ok(table);
    }
    let base = base_path(&table[0].path)?;
    let queue = Atomic::new(source_child_path(&base, "radix_sort_queue"));
    let out_path = source_child_path(&base, "radix_sort_out");
    let out = FileStack::try_new_owned(&out_path, job.io_errors())?;
    let mut processed = 0;
    loop {
        let index = usize::try_from(queue.fetch_add_one()?).map_err(|error| error.to_string())?;
        let Some(bucket) = table.get(index) else {
            break;
        };
        // VolumedFile snapshots the manifest once for both the diagnostic
        // and the sorting decision. Other workers may still publish files.
        let volumes = manifest_records_owned(&bucket.path, &job.io_errors())?;
        let count: u64 = volumes.iter().map(|(_, count)| count).sum();
        job.log(&format!(
            "Radix sorting. Bucket={}/{} Records={} Size={}",
            index + 1,
            table.len(),
            crate::util::string::format(count),
            crate::util::string::format(count * record_size(kind) as u64)
        ))?;
        let mut text = Vec::new();
        for result in spill(
            job,
            bucket.clone(),
            kind,
            bits_unsorted,
            volumes,
            count,
            threads,
        )? {
            text.extend_from_slice(
                &job.io(manifest_line(&result.path, &[result.records.unwrap()]))?,
            );
        }
        if !text.is_empty() {
            out.push_bytes(&text)?;
        }
        processed += 1;
    }
    let finished = Atomic::new(source_child_path(&base, "radix_sort_finished"));
    finished.fetch_add(processed)?;
    finished.await_value(table.len() as i64)?;
    RadixedTable::from_manifest_owned(out_path, &job.io_errors())
}

/// One actual Source thread-local BufferArray sharing a FileArray. Each worker
/// retains its own codec state through all tasks on that thread, then flushes
/// every radix (including empty members) before the enclosing FileArray closes.
pub(crate) struct WorkerBucketBuffers {
    files: std::sync::Arc<std::sync::Mutex<BucketWriter>>,
    buffers: Vec<CompressedBuffer>,
    records: Vec<u64>,
    errors: IoErrorOwner,
    failure: std::sync::Arc<std::sync::Mutex<Option<String>>>,
}

impl WorkerBucketBuffers {
    pub(crate) fn new(
        files: std::sync::Arc<std::sync::Mutex<BucketWriter>>,
        count: usize,
        errors: IoErrorOwner,
        failure: std::sync::Arc<std::sync::Mutex<Option<String>>>,
    ) -> Self {
        Self {
            files,
            buffers: (0..count).map(|_| CompressedBuffer::default()).collect(),
            records: vec![0; count],
            errors,
            failure,
        }
    }

    pub(crate) fn write(&mut self, radix: usize, bytes: &[u8]) -> Result<(), String> {
        self.write_count(radix, bytes, 1)
    }

    /// Source BufferArray's byte overload counts FASTA bytes, while its
    /// record overload (alignment edges) counts serialized records.
    pub(crate) fn write_count(
        &mut self,
        radix: usize,
        bytes: &[u8],
        records: u64,
    ) -> Result<(), String> {
        self.errors.capture(
            self.buffers[radix]
                .write(bytes)
                .map_err(std::io::Error::other),
        )?;
        self.records[radix] += records;
        if self.buffers[radix].size() >= BUFFER_BYTES {
            self.flush(radix, true)?;
        }
        Ok(())
    }

    fn flush(&mut self, radix: usize, reset: bool) -> Result<(), String> {
        self.errors
            .capture(self.buffers[radix].finish().map_err(std::io::Error::other))?;
        self.files.lock().unwrap().append_worker_member(
            radix,
            self.buffers[radix].data(),
            self.records[radix],
            &self.errors,
        )?;
        if reset {
            self.buffers[radix].clear();
            self.records[radix] = 0;
        }
        Ok(())
    }
}

impl Drop for WorkerBucketBuffers {
    fn drop(&mut self) {
        let source = self.errors.source_destructors();
        if source {
            if let Err(error) = self.errors.emit_pending_perror() {
                crate::util::io::abort_source_diagnostic_error(&error);
            }
            self.errors = self.errors.cleanup_owner();
        }
        for radix in 0..self.buffers.len() {
            if let Err(message) = self.flush(radix, false) {
                if source {
                    let error = self.errors.finish::<()>(Err(message)).unwrap_err();
                    if let Err(error) = crate::util::io::emit_source_perror(&error) {
                        crate::util::io::abort_source_diagnostic_error(&error);
                    }
                    crate::util::io::abort_source_process();
                }
                self.failure.lock().unwrap().get_or_insert(message);
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::InputFile;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn temp_dir() -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond-storage-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    #[cfg(unix)]
    #[test]
    fn manifest_roundtrip_retains_actual_non_utf8_volume_path() {
        use std::os::unix::ffi::OsStringExt;
        let root = temp_dir();
        let path = root.join(std::ffi::OsString::from_vec(b"volume_\xff\x80".to_vec()));
        let manifest = root.join("bucket.tsv");
        fs::write(&manifest, manifest_line(&path, &[u64::MAX]).unwrap()).unwrap();
        assert_eq!(manifest_records(&manifest).unwrap(), vec![(path, u64::MAX)]);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn file_array_preserves_source_repeated_separators_in_stored_names() {
        let root = temp_dir();
        let mut base = root.as_os_str().to_owned();
        base.push(std::path::MAIN_SEPARATOR.to_string());
        let mut files =
            BucketWriter::new_file_array_owned(&PathBuf::from(base), 3, 1, IoErrorOwner::default())
                .unwrap();
        let mut expected = root.as_os_str().to_owned();
        let sep = std::path::MAIN_SEPARATOR;
        expected.push(format!("{sep}{sep}0{sep}bucket.tsv"));
        assert_eq!(files.buckets()[0].path.as_os_str(), expected.as_os_str());
        let mut counter = root.as_os_str().to_owned();
        counter.push(format!("{sep}{sep}radix_sort_queue"));
        assert_eq!(
            source_child_path(&files.directory, "radix_sort_queue").as_os_str(),
            counter.as_os_str()
        );
        files.finish().unwrap();
        drop(files);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn output_constructor_retains_failing_volume_path_and_os_reason() {
        let root = temp_dir();
        let volume = root.join("0/worker_7_volume_0");
        fs::create_dir_all(&volume).unwrap();
        let errors = IoErrorOwner::default();
        let result = BucketWriter::new_owned(&root, 7, 1, errors.clone());
        let error = errors.finish(result.map(|_| ())).unwrap_err();
        let open = crate::util::io::file_open_error::find(&error)
            .expect("actual OutputFile constructor error must retain FileOpenError");
        assert_eq!(Path::new(&open.path), volume);
        assert!(!open.reason.is_empty());
        assert!(open.raw_os_error.is_some());
        #[cfg(unix)]
        assert_eq!(open.raw_os_error, Some(libc::EISDIR));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn scope_exit_finalizes_records_and_removes_empty_volumes_once() {
        let root = temp_dir();
        {
            let mut output = BucketWriter::new(&root, 3, 2).unwrap();
            output.write(1, b"pending record").unwrap();
        }
        assert!(!root.join("0/worker_3_volume_0").exists());
        assert!(manifest_records(&root.join("0/bucket.tsv"))
            .unwrap()
            .is_empty());
        let records = manifest_records(&root.join("1/bucket.tsv")).unwrap();
        assert_eq!(records, vec![(root.join("1/worker_3_volume_0"), 1)]);
        let mut input = InputFile::source_new(records[0].0.to_str().unwrap(), 0).unwrap();
        let mut decoded = Vec::new();
        std::io::Read::read_to_end(&mut input, &mut decoded).unwrap();
        input.close().unwrap();
        assert_eq!(decoded, b"pending record");
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn scope_exit_keeps_the_primary_input_open_error() {
        let root = temp_dir();
        let errors = IoErrorOwner::default();
        let missing = root.join("missing-input");
        let result = {
            let _output = BucketWriter::new_owned(&root, 5, 2, errors.clone()).unwrap();
            errors.capture(
                InputFile::with_source_profile(missing.to_str().unwrap(), 0, BUFFER_BYTES)
                    .map_err(std::io::Error::other),
            )
        };
        let error = errors.finish(result.map(|_| ())).unwrap_err();
        let open = crate::util::io::file_open_error::find(&error).unwrap();
        assert_eq!(Path::new(&open.path), missing);
        assert!(!root.join("0/worker_5_volume_0").exists());
        assert!(!root.join("1/worker_5_volume_0").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn file_array_owner_allocates_no_worker_codec_state() {
        let root = temp_dir();
        let mut files =
            BucketWriter::new_file_array_owned(&root, 8, 3, IoErrorOwner::default()).unwrap();
        assert!(files
            .buckets
            .iter()
            .all(|bucket| bucket.compressed.is_none()));
        assert!(files.buckets.iter().all(|bucket| !bucket.buffer_active));
        files.flush_unpublished().unwrap();
        assert!(files
            .buckets
            .iter()
            .all(|bucket| bucket.compressed.is_none()));
        files.finish().unwrap();
        drop(files);
        assert!((0..3).all(|radix| !root
            .join(radix.to_string())
            .join("worker_8_volume_0")
            .exists()));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn file_array_opens_empty_volumes_and_removes_them_at_close() {
        let root = temp_dir();
        let mut output = BucketWriter::new(&root, 7, 3).unwrap();
        for i in 0..3 {
            assert!(root.join(i.to_string()).join("worker_7_volume_0").is_file());
        }
        output.write(1, b"a record").unwrap();
        output.finish().unwrap();
        assert!(!root.join("0/worker_7_volume_0").exists());
        assert!(!root.join("2/worker_7_volume_0").exists());
        assert!(manifest_records(&root.join("0/bucket.tsv"))
            .unwrap()
            .is_empty());
        let records = manifest_records(&root.join("1/bucket.tsv")).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].1, 1);
        drop(output);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn full_volume_reopens_next_owner_before_another_record_arrives() {
        let root = temp_dir();
        let mut output = BucketWriter::new(&root, 4, 1).unwrap();
        // Exercise the rollover branch without allocating a 1-GiB fixture.
        output.buckets[0].volume_bytes = MAX_VOLUME_BYTES - 1;
        output.write(0, b"one record").unwrap();
        output.flush(0).unwrap();
        let records = manifest_records(&root.join("0/bucket.tsv")).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].1, 1);
        assert!(root.join("0/worker_4_volume_0").metadata().unwrap().len() > 0);
        assert_eq!(
            root.join("0/worker_4_volume_1").metadata().unwrap().len(),
            0
        );
        output.finish().unwrap();
        assert!(!root.join("0/worker_4_volume_1").exists());
        assert_eq!(
            manifest_records(&root.join("0/bucket.tsv")).unwrap(),
            records
        );
        drop(output);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn exclusive_file_array_retains_counts_across_rollover() {
        let root = temp_dir();
        let mut files =
            BucketWriter::new_exclusive_file_array_owned(&root, 4, 1, IoErrorOwner::default())
                .unwrap();
        let mut member = CompressedBuffer::new();
        member.write(b"first").unwrap();
        member.finish().unwrap();
        files.buckets[0].volume_bytes = MAX_VOLUME_BYTES - 1;
        files.append_member(0, member.data(), 1).unwrap();
        assert_eq!(files.buckets().0[0].records, Some(1));
        member.clear();
        member.write(b"second").unwrap();
        member.finish().unwrap();
        files.append_member(0, member.data(), 1).unwrap();
        let buckets = files.buckets();
        files.finish().unwrap();
        assert_eq!(buckets.0[0].records, Some(2));
        assert_eq!(
            manifest_records(&buckets.0[0].path)
                .unwrap()
                .iter()
                .map(|(_, n)| *n)
                .sum::<u64>(),
            2
        );
        drop(files);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn rollover_publishes_zero_record_member_before_empty_final_file_removal() {
        let root = temp_dir();
        let mut files =
            BucketWriter::new_exclusive_file_array_owned(&root, 4, 1, IoErrorOwner::default())
                .unwrap();
        let mut member = CompressedBuffer::new();
        member.finish().unwrap();
        files.buckets[0].volume_bytes = MAX_VOLUME_BYTES - 1;
        files.append_member(0, member.data(), 0).unwrap();
        assert_eq!(
            manifest_records(&root.join("0/bucket.tsv")).unwrap(),
            vec![(root.join("0/worker_4_volume_0"), 0)]
        );
        assert!(root.join("0/worker_4_volume_0").exists());
        files.finish().unwrap();
        assert!(!root.join("0/worker_4_volume_1").exists());
        drop(files);
        fs::remove_dir_all(root).unwrap();
    }

    fn seed_record(seed: u64, oid: i64) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&seed.to_ne_bytes());
        bytes.extend_from_slice(&oid.to_ne_bytes());
        bytes.extend_from_slice(&100i32.to_ne_bytes());
        bytes
    }

    #[test]
    fn oversized_radix_workers_route_single_record_volumes_and_keep_source_inputs() {
        let root = temp_dir();
        let input_dir = root.join("seeds/0");
        fs::create_dir_all(&input_dir).unwrap();
        let mut manifest = Vec::new();
        let mut expected = Vec::new();
        let mut original_paths = Vec::new();
        for index in 0..8 {
            let path = input_dir.join(format!("input-{index}"));
            let bytes = seed_record(index + 1, index as i64);
            fs::write(&path, &bytes).unwrap();
            expected.push(bytes);
            original_paths.push(path.clone());
            manifest.extend(manifest_line(&path, &[1]).unwrap());
        }
        let manifest_path = input_dir.join("bucket.tsv");
        fs::write(&manifest_path, manifest).unwrap();
        let job = Job::new(8, 8, 1, root.join("job"), 7).unwrap();
        let children = radix_sort_with_threads(
            &job,
            RadixedTable(vec![Bucket::new(&manifest_path, Some(8))]),
            ExternalRecordKind::Seed,
            RADIX_BITS,
            4,
        )
        .unwrap();
        assert_eq!(children.len(), RADIX_COUNT);
        assert_eq!(
            children
                .iter()
                .map(|child| child.records.unwrap())
                .sum::<u64>(),
            8
        );
        assert!(original_paths.iter().all(|path| path.exists()));
        assert!(
            manifest_path.exists(),
            "Source retains the oversized input after this one pass"
        );
        let mut actual = Vec::new();
        for (radix, child) in children.iter().enumerate() {
            let mut counted = 0;
            for (path, records) in manifest_records(&child.path).unwrap() {
                let mut input =
                    InputFile::with_source_profile_path(&path, 0, BUFFER_BYTES).unwrap();
                for _ in 0..records {
                    let mut record = vec![0; 20];
                    std::io::Read::read_exact(&mut input, &mut record).unwrap();
                    assert_eq!(
                        (key(&record, ExternalRecordKind::Seed) & 255) as usize,
                        radix
                    );
                    actual.push(record);
                }
                counted += records;
                input.close().unwrap();
            }
            assert_eq!(child.records, Some(counted));
        }
        actual.sort();
        expected.sort();
        assert_eq!(actual, expected);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn radix_worker_eof_retains_the_previous_record_within_one_volume() {
        let root = temp_dir();
        let input_dir = root.join("seeds/0");
        fs::create_dir_all(&input_dir).unwrap();
        let path = input_dir.join("one-physical-record");
        let bytes = seed_record(12, 4);
        fs::write(&path, &bytes).unwrap();
        let manifest = input_dir.join("bucket.tsv");
        fs::write(&manifest, manifest_line(&path, &[2]).unwrap()).unwrap();
        let job = Job::new(4, 1, 1, root.join("job"), 2).unwrap();
        let children = radix_sort_with_threads(
            &job,
            RadixedTable(vec![Bucket::new(&manifest, Some(2))]),
            ExternalRecordKind::Seed,
            RADIX_BITS,
            4,
        )
        .unwrap();
        let radix = (key(&bytes, ExternalRecordKind::Seed) & 255) as usize;
        assert_eq!(children[radix].records, Some(2));
        let rows = manifest_records(&children[radix].path).unwrap();
        let mut input = InputFile::with_source_profile_path(&rows[0].0, 0, BUFFER_BYTES).unwrap();
        let mut actual = vec![0; 40];
        std::io::Read::read_exact(&mut input, &mut actual).unwrap();
        input.close().unwrap();
        assert_eq!(actual, [bytes.clone(), bytes].concat());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn radix_worker_preserves_the_missing_input_file_error() {
        let root = temp_dir();
        let input_dir = root.join("seeds/0");
        fs::create_dir_all(&input_dir).unwrap();
        let missing = input_dir.join("missing-volume");
        let manifest = input_dir.join("bucket.tsv");
        fs::write(&manifest, manifest_line(&missing, &[2]).unwrap()).unwrap();
        let errors = IoErrorOwner::default();
        let job = Job::new_owned(4, 1, 1, root.join("job"), 2, errors.clone()).unwrap();
        let result = radix_sort_with_threads(
            &job,
            RadixedTable(vec![Bucket::new(&manifest, Some(2))]),
            ExternalRecordKind::Seed,
            RADIX_BITS,
            4,
        );
        let error = errors.finish(result).unwrap_err();
        assert_eq!(
            crate::util::io::file_open_error::find(&error).unwrap().path,
            missing
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn buffer_threshold_counts_compressed_output_instead_of_input() {
        let root = temp_dir();
        let mut output = BucketWriter::new(&root, 0, 1).unwrap();
        output.write(0, &vec![b'A'; BUFFER_BYTES * 4]).unwrap();
        assert_eq!(output.buckets[0].pending_records, 1);
        assert_eq!(output.buckets[0].volume_bytes, 0);
        output.finish().unwrap();
        let records = manifest_records(&root.join("0/bucket.tsv")).unwrap();
        assert_eq!(records.len(), 1);
        let mut input = InputFile::source_new(records[0].0.to_str().unwrap(), 0).unwrap();
        let mut bytes = vec![0; BUFFER_BYTES * 4];
        assert_eq!(input.read_raw(&mut bytes).unwrap(), bytes.len());
        assert!(bytes.iter().all(|&b| b == b'A'));
        input.close().unwrap();
        drop(output);
        fs::remove_dir_all(root).unwrap();
    }
}
