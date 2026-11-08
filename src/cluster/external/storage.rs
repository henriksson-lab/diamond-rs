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

const BUFFER_BYTES: usize = 65_536;
const MAX_VOLUME_BYTES: u64 = 1 << 30;

#[derive(Debug, Default)]
struct BucketOutput {
    compressed: CompressedBuffer,
    buffer_active: bool,
    pending_records: u64,
    file: Option<OutputFile>,
    stack: Option<FileStack>,
    volume: u64,
    volume_bytes: u64,
    volume_records: u64,
}

#[derive(Debug)]
pub(crate) struct BucketWriter {
    directory: PathBuf,
    worker: i64,
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
        let mut buckets = Vec::with_capacity(count);
        for radix in 0..count {
            let dir = directory.join(radix.to_string());
            fs::create_dir_all(source_c_path(&dir).as_ref()).map_err(|error| error.to_string())?;
            // FileArray opens every initial volume before its FileStack,
            // including buckets which will receive no records.
            let file = Self::open_volume(&dir.join(format!("worker_{worker}_volume_0")), &errors)?;
            let stack = FileStack::try_new_owned(dir.join("bucket.tsv"), errors.clone())?;
            buckets.push(BucketOutput {
                file: Some(file),
                stack: Some(stack),
                buffer_active: true,
                ..BucketOutput::default()
            });
        }
        Ok(Self {
            directory: directory.into(),
            worker,
            buckets,
            error: None,
            errors,
        })
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
                .write(bytes)
                .map_err(std::io::Error::other),
        )?;
        bucket.buffer_active = true;
        bucket.pending_records += records;
        // BufferArray tests CompressedBuffer::size(), not the input size.
        if bucket.compressed.size() >= BUFFER_BYTES {
            self.flush(radix)?;
        }
        Ok(())
    }

    fn volume_path(&self, radix: usize) -> PathBuf {
        self.directory.join(radix.to_string()).join(format!(
            "worker_{}_volume_{}",
            self.worker, self.buckets[radix].volume
        ))
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
        self.errors
            .capture(bucket.compressed.finish().map_err(std::io::Error::other))?;
        self.errors.capture(
            bucket
                .file
                .as_mut()
                .unwrap()
                .write_raw(bucket.compressed.data())
                .map_err(std::io::Error::other),
        )?;
        bucket.volume_bytes += bucket.compressed.size() as u64;
        bucket.volume_records += bucket.pending_records;
        bucket.pending_records = 0;
        bucket.compressed.clear();
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
            if bucket.volume_records > 0 {
                // FileStack takes the Source cross-process lock around each line.
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
            let manifest = self.directory.join(radix.to_string()).join("bucket.tsv");
            // Other process workers can still append. Read aggregate counts
            // only after the stage's shared completion barrier.
            buckets.push(Bucket::new(manifest, None));
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
            self.buckets[radix].buffer_active = true;
            self.flush(radix)?;
            self.buckets[radix].buffer_active = false;
            buckets.push(Bucket::new(
                self.directory.join(radix.to_string()).join("bucket.tsv"),
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
    let bytes = fs::read(source_c_path(&(path)).as_ref()).map_err(|error| error.to_string())?;
    bytes
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.iter().all(u8::is_ascii_whitespace))
        .map(|line| {
            let mut fields = line
                .split(u8::is_ascii_whitespace)
                .filter(|field| !field.is_empty());
            let path = fields.next().ok_or("Format error in VolumedFile")?;
            let count = std::str::from_utf8(fields.next().ok_or("Format error in VolumedFile")?)
                .ok()
                .and_then(|value| value.parse::<u64>().ok())
                .ok_or("Format error in VolumedFile")?;
            Ok((
                PathBuf::from(errors.capture(
                    crate::util::io::source_bytes_os(path).map_err(std::io::Error::other),
                )?),
                count,
            ))
        })
        .collect()
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
) -> Result<Vec<Bucket>, String> {
    if count == 0 {
        for (path, _) in volumes {
            fs::remove_file(source_c_path(&(path)).as_ref()).map_err(|error| error.to_string())?;
        }
        fs::remove_file(source_c_path(&bucket.path).as_ref()).map_err(|error| error.to_string())?;
        if let Some(directory) = bucket.path.parent() {
            let _ = fs::remove_dir(source_c_path(&(directory)).as_ref());
        }
        return Ok(Vec::new());
    }
    if bits < RADIX_BITS || count.saturating_mul(record_size(kind) as u64) <= job.mem_limit {
        return Ok(vec![Bucket::new(bucket.path, Some(count))]);
    }
    let shift = bits - RADIX_BITS;
    let mut writer = BucketWriter::new_owned(
        bucket.path.parent().unwrap_or(Path::new("")),
        job.worker_id(),
        RADIX_COUNT,
        job.io_errors(),
    )?;
    for (path, records) in &volumes {
        let mut input = job.open_input(path, 0)?;
        let mut bytes = vec![0; record_size(kind)];
        for _ in 0..*records {
            let mut read = 0;
            while read < bytes.len() {
                let n = job.io(input.read_raw(&mut bytes[read..]))?;
                if n == 0 {
                    return Err("Short external radix record".into());
                }
                read += n;
            }
            writer.write(((key(&bytes, kind) >> shift) & 255) as usize, &bytes)?;
        }
    }
    let mut children = writer.finish()?;
    drop(writer);
    // Source radix_cluster has an exclusive output owner, so all 256 child
    // counts (including zero) can be read without another process appending.
    for child in children.iter_mut() {
        child.records = Some(
            manifest_records_owned(&child.path, &job.io_errors())?
                .iter()
                .map(|(_, n)| n)
                .sum(),
        );
    }
    job.log(&format!("Radix sorted bucket records={count}"))?;
    // Exactly one eight-bit pass, as in radix_sort.h. Source retains the old
    // input bucket here and does not recurse even when a child still exceeds
    // mem_limit. Thus the configured memory limit remains a heuristic on
    // skewed distributions; changing pass depth can change chunk boundaries.
    Ok(children.0)
}

pub fn radix_sort(
    job: &Job,
    table: RadixedTable,
    kind: ExternalRecordKind,
    bits_unsorted: u32,
) -> Result<RadixedTable, String> {
    if bits_unsorted < RADIX_BITS || table.is_empty() {
        return Ok(table);
    }
    let base = base_path(&table[0].path)?;
    let queue = Atomic::new(base.join("radix_sort_queue"));
    let out_path = base.join("radix_sort_out");
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
        for result in spill(job, bucket.clone(), kind, bits_unsorted, volumes, count)? {
            text.extend_from_slice(
                &job.io(manifest_line(&result.path, &[result.records.unwrap()]))?,
            );
        }
        if !text.is_empty() {
            out.push_bytes(&text)?;
        }
        processed += 1;
    }
    let finished = Atomic::new(base.join("radix_sort_finished"));
    finished.fetch_add(processed)?;
    finished.await_value(table.len() as i64)?;
    RadixedTable::from_manifest_owned(out_path, &job.io_errors())
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
