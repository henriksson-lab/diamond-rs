//! External-memory clustering from `diamond/src/cluster/external/cluster.cpp`.
//!
//! File-backed queues distribute bucket and volume ownership across workers.
//! Closure ownership and completion barriers retain Source's shared job
//! protocol, native record layouts, and representative manifests.

use crate::util::io::source_c_path;
use crate::util::parallel::{Atomic, FileStack};
use std::cell::RefCell;
use std::fs;
use std::path::{Path, PathBuf};

use super::external::{read_edge_bucket, sort_input_buffer, Assignment, Edge, Job};
#[cfg(test)]
use super::pair_table::Bucket;
use super::pair_table::RadixedTable;
use super::seed_table::{Volume, VolumedFile};
use super::storage::{BucketWriter, WorkerBucketBuffers};
use crate::basic::sequence_utils;
use crate::basic::value::{OId, AMINO_ACID_ALPHABET};
use crate::util::algo::Partition;
#[cfg(test)]
use crate::util::io::CompressedBuffer;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

const EDGE_ENCODED_SIZE: usize = 24;
#[cfg(test)]
const ASSIGNMENT_ENCODED_SIZE: usize = 16;

fn representative_count(volumes: &VolumedFile) -> Result<usize, String> {
    let max_oid = volumes
        .0
        .iter()
        .filter_map(|volume| volume.oid_end.checked_sub(1))
        .max()
        .unwrap_or(0);
    usize::try_from(max_oid)
        .ok()
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| "Clustering OID range does not fit memory.".to_owned())
}

/// Source uses plain std::thread for these workers, unlike InputBuffer's
/// exception-catching SimpleThreadPool. Keep CLI termination at that boundary;
/// reusable library entry points retain their explicit fallible extension.
fn plain_workers<F>(job: &Job, count: usize, worker: F) -> Result<(), String>
where
    F: Fn(usize) -> Result<(), String> + Sync,
{
    let failure = Mutex::new(None);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for index in 0..count {
            let worker = &worker;
            let failure = &failure;
            match std::thread::Builder::new().spawn_scoped(scope, move || {
                if let Err(message) = worker(index) {
                    if job.io_errors().source_destructors() {
                        let error = job.io_errors().finish::<()>(Err(message)).unwrap_err();
                        if let Err(error) = crate::util::io::emit_source_perror(&error) {
                            crate::util::io::abort_source_diagnostic_error(&error);
                        }
                        crate::util::io::abort_source_process();
                    }
                    failure.lock().unwrap().get_or_insert(message);
                }
            }) {
                Ok(handle) => handles.push(handle),
                Err(error) => {
                    failure
                        .lock()
                        .unwrap()
                        .get_or_insert_with(|| job.capture_io::<()>(Err(error)).unwrap_err());
                    break;
                }
            }
        }
        for handle in handles {
            handle.join().unwrap();
        }
    });
    failure.into_inner().unwrap().map_or(Ok(()), Err)
}

fn identity_mapping(count: usize) -> Vec<AtomicU64> {
    (0..count as OId).map(AtomicU64::new).collect()
}

fn root_of(rep: &[AtomicU64], start: OId) -> Result<OId, String> {
    let mut root = start;
    for _ in 0..rep.len() {
        let next = rep
            .get(root as usize)
            .ok_or_else(|| format!("Representative OID out of range: {root}"))?
            .load(Ordering::Relaxed);
        if next == root {
            return Ok(root);
        }
        root = next;
    }
    Err(format!("Representative cycle contains OID {start}."))
}

/// C++ compute_closure's actual shared mapping owner. Relaxed atoms preserve
/// pointer compression without aliasing concurrent mutable slices. Source uses
/// non-atomic reads/writes and can race across partition-spanning chains; the
/// native owner defines that unsafe domain while retaining valid forest roots.
fn compute_closure_mapping(
    job: &Job,
    volumes: &VolumedFile,
    rep: &[AtomicU64],
    threads: usize,
) -> Result<(), String> {
    let expected = representative_count(volumes)?;
    if rep.len() < expected {
        return Err(format!(
            "Representative vector has {} entries, expected at least {expected}.",
            rep.len()
        ));
    }
    let parts = Partition::new(expected, threads);
    plain_workers(job, parts.parts, |part| {
        for oid in parts.begin(part)..parts.end(part) {
            let previous = rep[oid].load(Ordering::Relaxed);
            let root = root_of(rep, previous)?;
            if previous != root {
                rep[oid].store(root, Ordering::Relaxed);
            }
        }
        Ok(())
    })?;

    // Source writes closure files serially, after joining compression workers.
    let output_dir = job.base_dir(None).join("clustering");
    job.capture_io(crate::util::system::mkdir_path(&output_dir))?;
    for (index, volume) in volumes.0.iter().enumerate() {
        let begin = usize::try_from(volume.oid_begin)
            .map_err(|_| "Volume OID range does not fit memory.".to_owned())?;
        let end = usize::try_from(volume.oid_end)
            .map_err(|_| "Volume OID range does not fit memory.".to_owned())?;
        let slice = rep
            .get(begin..end)
            .ok_or_else(|| format!("Invalid OID range {begin}..{end}."))?;
        let mut bytes = Vec::with_capacity(slice.len() * std::mem::size_of::<OId>());
        for oid in slice {
            bytes.extend_from_slice(&oid.load(Ordering::Relaxed).to_ne_bytes());
        }
        let mut output = job.create_output(&output_dir.join(format!("volume{index}")))?;
        job.io(output.write_raw(&bytes))?;
        job.io(output.close())?;
    }
    Ok(())
}

/// Public compatibility surface retains the original one-worker default.
pub fn compute_closure(job: &Job, volumes: &VolumedFile, rep: &mut [OId]) -> Result<(), String> {
    compute_closure_with_threads(job, volumes, rep, 1)
}

pub fn compute_closure_with_threads(
    job: &Job,
    volumes: &VolumedFile,
    rep: &mut [OId],
    threads: usize,
) -> Result<(), String> {
    // Only this legacy borrowed-slice adapter copies values. Actual command
    // paths create their one canonical atomic mapping owner directly.
    let mapping: Vec<_> = rep.iter().copied().map(AtomicU64::new).collect();
    let result = compute_closure_mapping(job, volumes, &mapping, threads);
    for (destination, value) in rep.iter_mut().zip(&mapping) {
        *destination = value.load(Ordering::Relaxed);
    }
    result
}

pub fn compute_closure_from_assignment_file(
    job: &Job,
    assignment_file: &Path,
    volumes: &VolumedFile,
) -> Result<(), String> {
    compute_closure_from_assignment_file_with_threads(job, assignment_file, volumes, 1)
}

/// Source's assignment reader consumes complete records until physical EOF;
/// manifest record_count does not bound this loop.
pub fn compute_closure_from_assignment_file_with_threads(
    job: &Job,
    assignment_file: &Path,
    volumes: &VolumedFile,
    threads: usize,
) -> Result<(), String> {
    job.log("Computing transitive closure")?;
    let rep = identity_mapping(representative_count(volumes)?);
    let assignment_volumes =
        super::external::read_manifest_owned(assignment_file, &job.io_errors())?;
    let next = AtomicUsize::new(0);
    plain_workers(job, threads.min(assignment_volumes.len()), |thread_id| {
        loop {
            let index = next.fetch_add(1, Ordering::Relaxed);
            let Some((path, _)) = assignment_volumes.get(index) else {
                break;
            };
            job.log(&format!(
                "Reading assignments thread_id={thread_id} volume={index}"
            ))?;
            let mut input = job.open_input(path, 0)?;
            let mut bytes = [0; 16];
            while job.io(input.read_raw(&mut bytes))? == bytes.len() {
                let assignment = Assignment {
                    member_oid: OId::from_ne_bytes(bytes[..8].try_into().unwrap()),
                    rep_oid: OId::from_ne_bytes(bytes[8..].try_into().unwrap()),
                };
                if assignment.rep_oid as usize >= rep.len() {
                    return Err(format!(
                        "Assignment representative OID out of range: {}",
                        assignment.rep_oid
                    ));
                }
                let member = rep.get(assignment.member_oid as usize).ok_or_else(|| {
                    format!(
                        "Assignment member OID out of range: {}",
                        assignment.member_oid
                    )
                })?;
                // Generated member groups belong to exactly one radix bucket
                // and part. Atoms also avoid UB on conflicting custom inputs.
                member.store(assignment.rep_oid, Ordering::Relaxed);
            }
            job.io(input.close())?;
        }
        Ok(())
    })?;
    compute_closure_mapping(job, volumes, &rep, threads)?;
    remove_manifest_storage(assignment_file, &assignment_volumes)?;
    Ok(())
}

/// Library compatibility entry with an explicit single-thread context.
pub fn get_reps(job: &Job, volumes: &VolumedFile) -> Result<PathBuf, String> {
    get_reps_with_threads(job, volumes, 1)
}

/// Source get_reps creates exactly config.threads workers (not min(volumes,N));
/// every worker performs its terminal queue claim, including idle workers.
pub fn get_reps_with_threads(
    job: &Job,
    volumes: &VolumedFile,
    threads: usize,
) -> Result<PathBuf, String> {
    if job.last_round() {
        return Ok(PathBuf::new());
    }
    let base_dir = job.base_dir(None).join("reps");
    job.capture_io(crate::util::system::mkdir_path(&base_dir))?;
    let manifest_path = base_dir.join("reps.tsv");
    let manifest = std::sync::Arc::new(FileStack::try_new_owned(&manifest_path, job.io_errors())?);
    let queue = std::sync::Arc::new(Atomic::new(base_dir.join("queue")));
    let processed = std::sync::Arc::new(std::sync::atomic::AtomicI64::new(0));
    let cluster_count = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
    let volumes = std::sync::Arc::new(volumes.clone());
    let mut pool = crate::util::parallel::SimpleThreadPool::new();
    for _ in 0..threads {
        let job = job.clone();
        let volumes = std::sync::Arc::clone(&volumes);
        let manifest = std::sync::Arc::clone(&manifest);
        let queue = std::sync::Arc::clone(&queue);
        let processed = std::sync::Arc::clone(&processed);
        let cluster_count = std::sync::Arc::clone(&cluster_count);
        let base_dir = base_dir.clone();
        pool.spawn(move |stop| {
            loop {
                let index = usize::try_from(queue.fetch_add_one()?)
                    .map_err(|_| "Negative representative volume index")?;
                if stop.load(std::sync::atomic::Ordering::Relaxed) {
                    break;
                }
                let Some(volume) = volumes.0.get(index) else {
                    break;
                };
                job.log(&format!(
                    "Writing representatives. Volume={}/{} Records={}",
                    index + 1,
                    volumes.0.len(),
                    crate::util::string::format(volume.record_count)
                ))?;
                crate::util::io::console::flush_c_stdout();
                let representatives = read_clustering_volume(&job, index, volume)?;
                let path = base_dir.join(format!("{index}.faa"));
                // Source opens the actual input before its unchecked textmode
                // ofstream; exceptions are disabled and failure stays sticky.
                let output = RefCell::new(None::<crate::util::io::source_text::TextOutput>);
                let mut table_oid = volume.oid_begin;
                let mut table_index = 0usize;
                let mut next_oid = volume.oid_begin;
                let mut count = 0u64;
                super::reader::visit_volume_bytes_while_owned_with_cleanup(
                    &volume.path,
                    &mut || !stop.load(std::sync::atomic::Ordering::Relaxed),
                    &mut |record| {
                        let oid = if job.round_index() > 0 {
                            job.capture_io(parse_oid_like_atoll_bytes(&record.id))?
                        } else {
                            next_oid
                        };
                        while table_oid < oid {
                            table_oid += 1;
                            table_index += 1;
                        }
                        if table_oid != oid || table_index >= representatives.len() {
                            return Err(format!(
                                "Representative table does not contain OID {oid}."
                            ));
                        }
                        if representatives[table_index] == oid {
                            let mut bytes = format!(">{oid}\n").into_bytes();
                            bytes.extend_from_slice(
                                sequence_utils::to_string(&record.sequence, AMINO_ACID_ALPHABET)
                                    .as_bytes(),
                            );
                            bytes.push(b'\n');
                            if let Some(output) = output.borrow_mut().as_mut() {
                                output.write(&bytes);
                            }
                            count += 1;
                        }
                        next_oid = oid.checked_add(1).ok_or("Sequence OID overflow.")?;
                        table_oid += 1;
                        table_index += 1;
                        Ok(())
                    },
                    &mut || {
                        *output.borrow_mut() = Some(
                            crate::util::io::source_text::TextOutput::open_unchecked(&path, true),
                        );
                        Ok(())
                    },
                    &mut || {
                        drop(output.borrow_mut().take());
                    },
                    &job.io_errors(),
                )?;
                drop(output.into_inner()); // out.close after in.close, before publication.
                manifest.push_bytes(&job.io(super::storage::manifest_line(
                    &path,
                    &[count, volume.oid_begin, volume.oid_end],
                ))?)?;
                processed.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                cluster_count.fetch_add(count, std::sync::atomic::Ordering::Relaxed);
            }
            Ok(())
        });
    }
    pool.join_all()?;
    job.log(&format!(
        "Representatives written: {}",
        cluster_count.load(std::sync::atomic::Ordering::Relaxed)
    ))?;
    let finished = Atomic::new(base_dir.join("finished"));
    finished.fetch_add(processed.load(std::sync::atomic::Ordering::Relaxed))?;
    finished.await_value(volumes.0.len() as i64)?;
    drop(manifest);
    Ok(manifest_path)
}

/// C++ unidirectional `cluster` pass.
pub fn cluster(job: &Job, edges: &RadixedTable, volumes: &VolumedFile) -> Result<PathBuf, String> {
    cluster_with_threads(job, edges, volumes, 1)
}

pub fn cluster_with_threads(
    job: &Job,
    edges: &RadixedTable,
    volumes: &VolumedFile,
    threads: usize,
) -> Result<PathBuf, String> {
    let directory = job.base_dir(None).join("clustering");
    job.capture_io(crate::util::system::mkdir_path(&directory))?;
    let bucket_dir = directory.join("0");
    let output = Arc::new(Mutex::new(BucketWriter::new_file_array_owned(
        &directory,
        job.worker_id(),
        1,
        job.io_errors(),
    )?));
    let queue = Atomic::new(job.base_dir(None).join("alignments/queue"));
    let mut processed = 0;
    loop {
        let index = usize::try_from(queue.fetch_add_one()?)
            .map_err(|_| "Negative clustering bucket index")?;
        let Some(bucket) = edges.get(index) else {
            break;
        };
        let (mut records, storage) = read_edge_bucket(job, bucket, threads)?;
        job.log(&format!(
            "Clustering. Bucket={}/{} Records={} Size={}",
            index + 1,
            edges.len(),
            crate::util::string::format(records.len() as u64),
            crate::util::string::format((records.len() * EDGE_ENCODED_SIZE) as u64)
        ))?;
        sort_input_buffer(&mut records, threads)?;
        let ranges = edge_part_ranges(&records, threads);
        let failure = Arc::new(Mutex::new(None));
        plain_workers(job, ranges.len(), |part| {
            let mut buffer = WorkerBucketBuffers::new(
                Arc::clone(&output),
                1,
                job.io_errors(),
                Arc::clone(&failure),
            );
            let range = &ranges[part];
            let mut begin = range.start;
            while begin < range.end {
                let edge = records[begin];
                if edge.member_len < edge.rep_len
                    || (edge.member_len == edge.rep_len && edge.member_oid > edge.rep_oid)
                {
                    let mut bytes = [0; 16];
                    bytes[..8].copy_from_slice(&edge.member_oid.to_ne_bytes());
                    bytes[8..].copy_from_slice(&edge.rep_oid.to_ne_bytes());
                    buffer.write(0, &bytes)?;
                }
                begin += 1;
                while begin < range.end && records[begin].member_oid == edge.member_oid {
                    begin += 1;
                }
            }
            // Actual per-part BufferArray destruction flushes before join and
            // before the enclosing Source VolumedFile is removed.
            drop(buffer);
            Ok(())
        })?;
        if let Some(error) = failure.lock().unwrap().take() {
            return Err(error);
        }
        remove_manifest_storage(&bucket.path, &storage)?;
        processed += 1;
    }
    output.lock().unwrap().finish()?;
    drop(output);
    let manifest = bucket_dir.join("bucket.tsv");
    let finished = Atomic::new(directory.join("finished"));
    let previous = finished.fetch_add(processed)?;
    let closure_finished = Atomic::new(directory.join("closure_finished"));
    // Source elects through this completion counter alone. A late worker
    // with no claimed buckets can enter again after assignments were consumed;
    // its empty manifest then writes an identity closure. Preserve that
    // scheduling-dependent file protocol, without an additional owner lock.
    if previous + processed < edges.len() as i64 {
        closure_finished.await_value(1)?;
    } else {
        compute_closure_from_assignment_file_with_threads(job, &manifest, volumes, threads)?;
        closure_finished.fetch_add_one()?;
    }
    get_reps_with_threads(job, volumes, threads)
}

pub fn cluster_bidirectional(
    job: &Job,
    edges: &RadixedTable,
    volumes: &VolumedFile,
) -> Result<PathBuf, String> {
    cluster_bidirectional_with_threads(job, edges, volumes, 1)
}

pub fn cluster_bidirectional_with_threads(
    job: &Job,
    edges: &RadixedTable,
    volumes: &VolumedFile,
    threads: usize,
) -> Result<PathBuf, String> {
    let lock = Atomic::new(job.base_dir(None).join("cluster_bidirectional_lock"));
    let finished = Atomic::new(job.base_dir(None).join("cluster_bidirectional_finished"));
    if lock.fetch_add_one()? != 0 {
        finished.await_value(1)?;
        return get_reps_with_threads(job, volumes, threads);
    }
    job.log("Computing clustering (bi-directional coverage)")?;
    let count = representative_count(volumes)?;
    let mut storage = Vec::new();
    let mut degree = vec![0_u32; count];
    for (bucket_index, bucket) in edges.iter().enumerate() {
        let (bucket_edges, bucket_storage) = read_edge_bucket(job, bucket, threads)?;
        job.log(&format!(
            "Getting node degrees. Bucket={}/{} Records={} Size={}",
            bucket_index + 1,
            edges.len(),
            bucket_edges.len(),
            bucket_edges.len() * EDGE_ENCODED_SIZE
        ))?;
        for edge in &bucket_edges {
            *degree
                .get_mut(edge.member_oid as usize)
                .ok_or_else(|| format!("Edge member OID out of range: {}", edge.member_oid))? += 1;
            *degree.get_mut(edge.rep_oid as usize).ok_or_else(|| {
                format!("Edge representative OID out of range: {}", edge.rep_oid)
            })? += 1;
        }
        storage.push((bucket.path.clone(), bucket_storage));
    }
    let rep = identity_mapping(count);
    for (bucket_index, bucket) in edges.iter().enumerate() {
        let (bucket_edges, _) = read_edge_bucket(job, bucket, threads)?;
        job.log(&format!(
            "Assigning reps. Bucket={}/{} Records={} Size={}",
            bucket_index + 1,
            edges.len(),
            bucket_edges.len(),
            bucket_edges.len() * EDGE_ENCODED_SIZE
        ))?;
        for edge in &bucket_edges {
            assign_higher_degree(edge.member_oid, edge.rep_oid, &degree, &rep)?;
            assign_higher_degree(edge.rep_oid, edge.member_oid, &degree, &rep)?;
        }
    }
    compute_closure_mapping(job, volumes, &rep, threads)?;
    // Source publishes completion before its final VolumedFile removals.
    finished.fetch_add_one()?;
    for (manifest, bucket_storage) in storage {
        remove_manifest_storage(&manifest, &bucket_storage)?;
    }
    get_reps_with_threads(job, volumes, threads)
}

fn assign_higher_degree(
    node: OId,
    candidate: OId,
    degree: &[u32],
    rep: &[AtomicU64],
) -> Result<(), String> {
    let current = rep
        .get(node as usize)
        .ok_or_else(|| format!("Edge OID out of range: {node}"))?
        .load(Ordering::Relaxed);
    let candidate_degree = *degree
        .get(candidate as usize)
        .ok_or_else(|| format!("Edge OID out of range: {candidate}"))?;
    let current_degree = *degree
        .get(current as usize)
        .ok_or_else(|| format!("Representative OID out of range: {current}"))?;
    if candidate_degree > current_degree
        || (candidate_degree == current_degree && candidate < current)
    {
        rep[node as usize].store(candidate, Ordering::Relaxed);
    }
    Ok(())
}

fn read_clustering_volume(job: &Job, index: usize, volume: &Volume) -> Result<Vec<OId>, String> {
    let bytes = job.read_plain_file(
        &job.base_dir(None)
            .join("clustering")
            .join(format!("volume{index}")),
    )?;
    let expected = usize::try_from(volume.oid_end.saturating_sub(volume.oid_begin))
        .map_err(|_| "Volume OID range does not fit memory.".to_owned())?;
    if bytes.len() != expected * std::mem::size_of::<OId>() {
        return Err(format!(
            "Clustering volume has {} bytes, expected {}.",
            bytes.len(),
            expected * std::mem::size_of::<OId>()
        ));
    }
    Ok(bytes
        .chunks_exact(8)
        .map(|bytes| OId::from_ne_bytes(bytes.try_into().unwrap()))
        .collect())
}

#[cfg(test)]
fn write_assignments(job: &Job, assignments: &[Assignment]) -> Result<PathBuf, String> {
    let directory = job.base_dir(None).join("clustering").join("0");
    fs::create_dir_all(source_c_path(&directory).as_ref()).map_err(|error| error.to_string())?;
    let volume = directory.join(format!("worker_{}_volume_0", job.worker_id()));
    let manifest = directory.join("bucket.tsv");
    if assignments.is_empty() {
        fs::write(source_c_path(&manifest).as_ref(), []).map_err(|error| error.to_string())?;
    } else {
        let mut compressed = CompressedBuffer::new();
        for assignment in assignments {
            let mut bytes = Vec::with_capacity(ASSIGNMENT_ENCODED_SIZE);
            bytes.extend_from_slice(&assignment.member_oid.to_ne_bytes());
            bytes.extend_from_slice(&assignment.rep_oid.to_ne_bytes());
            compressed
                .write(&bytes)
                .map_err(|error| error.to_string())?;
        }
        compressed.finish().map_err(|error| error.to_string())?;
        fs::write(source_c_path(&volume).as_ref(), compressed.data())
            .map_err(|error| error.to_string())?;
        fs::write(
            source_c_path(&manifest).as_ref(),
            super::storage::manifest_line(&volume, &[assignments.len() as u64])
                .map_err(|e| e.to_string())?,
        )
        .map_err(|error| error.to_string())?;
    }
    Ok(manifest)
}

/// InputBuffer extends nominal ranges across equal Edge::key() values;
/// merge_keys subsequently groups actual member_oid values. Preserve even the
/// rare adjacent different-member hash collision's Source partition boundary.
fn edge_part_ranges(records: &[Edge], threads: usize) -> Vec<std::ops::Range<usize>> {
    let partition = Partition::new(records.len(), threads);
    let key = |index: usize| crate::util::hash::hash64(records[index].member_oid);
    (0..partition.parts)
        .map(|part| {
            let mut begin = partition.begin(part);
            let mut end = partition.end(part);
            if part > 0 {
                while begin < records.len() && key(begin - 1) == key(begin) {
                    begin += 1;
                }
            }
            while end < records.len() && key(end - 1) == key(end) {
                end += 1;
            }
            begin..end
        })
        .collect()
}

fn remove_manifest_storage(path: &Path, volumes: &[(PathBuf, usize)]) -> Result<(), String> {
    for (volume, _) in volumes {
        let _ = fs::remove_file(source_c_path(&(volume)).as_ref());
    }
    let _ = fs::remove_file(source_c_path(&(path)).as_ref());
    if let Some(directory) = path.parent() {
        let _ = fs::remove_dir(source_c_path(&(directory)).as_ref());
    }
    Ok(())
}

fn parse_oid_like_atoll(id: &str) -> Result<OId, String> {
    parse_oid_like_atoll_bytes(id.as_bytes()).map_err(|error| error.to_string())
}
fn raw_oid_error(prefix: &[u8], id: &[u8]) -> std::io::Error {
    let mut message = prefix.to_vec();
    message.extend_from_slice(id);
    crate::util::io::RawDiagnosticError(message).into_io()
}
fn parse_oid_like_atoll_bytes(id: &[u8]) -> std::io::Result<OId> {
    let bytes = crate::util::sequence::c_string_bytes(id);
    let mut index = 0;
    while bytes.get(index).is_some_and(u8::is_ascii_whitespace) {
        index += 1;
    }
    if bytes.get(index) == Some(&b'+') {
        index += 1;
    }
    let begin = index;
    let mut value = 0_u64;
    while let Some(digit) = bytes.get(index).and_then(|byte| byte.checked_sub(b'0')) {
        if digit > 9 {
            break;
        }
        value = value
            .checked_mul(10)
            .and_then(|value| value.checked_add(u64::from(digit)))
            .ok_or_else(|| raw_oid_error(b"Sequence OID overflow: ", id))?;
        index += 1;
    }
    if index == begin {
        return Err(raw_oid_error(b"Invalid sequence OID: ", id));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    fn temp_dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond_external_cluster_{name}_{}_{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    fn volumes(root: &Path, count: u64) -> VolumedFile {
        VolumedFile(vec![Volume::new(root.join("input.faa"), 0, count, count)])
    }

    fn edge_table(root: &Path, edges: &[Edge]) -> RadixedTable {
        let directory = root.join("edge_bucket");
        fs::create_dir_all(&directory).unwrap();
        let volume = directory.join("volume0");
        let manifest = directory.join("bucket.tsv");
        let mut compressed = CompressedBuffer::new();
        for edge in edges {
            let mut bytes = Vec::new();
            bytes.extend_from_slice(&edge.rep_oid.to_ne_bytes());
            bytes.extend_from_slice(&edge.member_oid.to_ne_bytes());
            bytes.extend_from_slice(&edge.rep_len.to_ne_bytes());
            bytes.extend_from_slice(&edge.member_len.to_ne_bytes());
            compressed.write(&bytes).unwrap();
        }
        compressed.finish().unwrap();
        fs::write(&volume, compressed.data()).unwrap();
        fs::write(
            &manifest,
            format!("{}\t{}\n", volume.display(), edges.len()),
        )
        .unwrap();
        RadixedTable(vec![Bucket::new(manifest, Some(edges.len() as u64))])
    }

    fn read_mapping(job: &Job, count: usize) -> Vec<OId> {
        let bytes = fs::read(job.base_dir(None).join("clustering/volume0")).unwrap();
        let values = bytes
            .chunks_exact(8)
            .map(|bytes| OId::from_ne_bytes(bytes.try_into().unwrap()))
            .collect::<Vec<_>>();
        assert_eq!(values.len(), count);
        values
    }

    #[test]
    fn compute_closure_flattens_chains_and_writes_volume_ranges() {
        let root = temp_dir("closure");
        let job = Job::new(4, 2, 1024, &root, 0).unwrap();
        let volumes = VolumedFile(vec![Volume::new("a", 0, 2, 2), Volume::new("b", 2, 5, 3)]);
        let mut rep = vec![0, 0, 1, 2, 4];
        compute_closure(&job, &volumes, &mut rep).unwrap();
        assert_eq!(rep, [0, 0, 0, 0, 4]);
        assert_eq!(
            fs::read(job.base_dir(None).join("clustering/volume0"))
                .unwrap()
                .len(),
            16
        );
        assert_eq!(
            fs::read(job.base_dir(None).join("clustering/volume1"))
                .unwrap()
                .len(),
            24
        );
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn cluster_selects_longest_then_lowest_representative_and_emits_reps() {
        let root = temp_dir("unidirectional");
        let mut job = Job::new(3, 1, 1024, &root, 7).unwrap();
        job.set_round_count(2);
        let volumes = volumes(&root, 4);
        fs::write(
            &volumes.0[0].path,
            b">zero\nARND\n>one\nARND\n>two\nARND\n>three\nARND\n",
        )
        .unwrap();
        let edges = edge_table(
            &root,
            &[
                Edge {
                    rep_oid: 2,
                    member_oid: 3,
                    rep_len: 9,
                    member_len: 4,
                },
                Edge {
                    rep_oid: 0,
                    member_oid: 3,
                    rep_len: 9,
                    member_len: 4,
                },
                Edge {
                    rep_oid: 0,
                    member_oid: 2,
                    rep_len: 9,
                    member_len: 9,
                },
                Edge {
                    rep_oid: 2,
                    member_oid: 0,
                    rep_len: 9,
                    member_len: 9,
                },
            ],
        );
        // Exercise representative selection and the storage-backed closure,
        // then use the same automatic FASTA adapter as the public workflow.
        let mut bucket_edges = read_edge_bucket(&job, &edges[0], 1).unwrap().0;
        bucket_edges.sort_unstable();
        assert_eq!(bucket_edges[0].member_oid, 0);
        let assignment_file = write_assignments(
            &job,
            &[
                Assignment {
                    member_oid: 2,
                    rep_oid: 0,
                },
                Assignment {
                    member_oid: 3,
                    rep_oid: 0,
                },
            ],
        )
        .unwrap();
        compute_closure_from_assignment_file(&job, &assignment_file, &volumes).unwrap();
        let reps = get_reps(&job, &volumes).unwrap();
        assert_eq!(read_mapping(&job, 4), [0, 1, 0, 0]);
        let manifest = fs::read_to_string(reps).unwrap();
        assert!(manifest.contains("\t2\t0\t4"));
        let fasta = fs::read_to_string(job.base_dir(None).join("reps/0.faa")).unwrap();
        let newline = if cfg!(windows) { "\r\n" } else { "\n" };
        assert!(fasta.contains(&format!(">0{newline}ARND{newline}")));
        assert!(fasta.contains(&format!(">1{newline}ARND{newline}")));
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn representatives_use_configured_workers_and_source_terminal_claims() {
        let root = temp_dir("parallel_representatives");
        let mut job = Job::new(3, 2, 1024, &root, 4).unwrap();
        job.set_round_count(2);
        let mut input = Vec::new();
        fs::create_dir_all(job.base_dir(None).join("clustering")).unwrap();
        for volume in 0..2 {
            let path = root.join(format!("input{volume}.faa"));
            fs::write(&path, b">raw_\xff\nARNDCQEGH\n>raw_\x80\nARNDCQEGH\n").unwrap();
            input.push(Volume::new(path, volume * 2, volume * 2 + 2, 2));
            let mut values = Vec::new();
            for _ in 0..2 {
                values.extend_from_slice(&(volume * 2).to_ne_bytes());
            }
            fs::write(
                job.base_dir(None)
                    .join(format!("clustering/volume{volume}")),
                values,
            )
            .unwrap();
        }
        let manifest = get_reps_with_threads(&job, &VolumedFile(input.clone()), 4).unwrap();
        assert_eq!(
            Atomic::new(job.base_dir(None).join("reps/queue"))
                .get()
                .unwrap(),
            6,
            "two successful volume claims plus four terminal worker claims"
        );
        let volumes =
            super::super::external::read_volumed_file_manifest_owned(&manifest, &job.io_errors())
                .unwrap();
        assert_eq!(volumes.0.len(), 2);
        let newline = if cfg!(windows) { "\r\n" } else { "\n" };
        for (index, volume) in volumes.0.iter().enumerate() {
            assert_eq!(volume.record_count, 1);
            assert_eq!(
                fs::read(&volume.path).unwrap(),
                format!(">{}{newline}ARNDCQEGH{newline}", index * 2).as_bytes()
            );
            assert!(
                input[index].path.exists(),
                "representative selection does not consume source inputs"
            );
        }
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn late_zero_bucket_worker_recomputes_source_identity_closure() {
        let root = temp_dir("late_closure");
        let mut first = Job::new(3, 1, 1024, &root, 0).unwrap();
        first.set_round_count(1);
        let volumes = volumes(&root, 4);
        fs::create_dir_all(first.base_dir(None).join("alignments")).unwrap();
        let edges = edge_table(
            &root,
            &[
                Edge {
                    rep_oid: 0,
                    member_oid: 2,
                    rep_len: 12,
                    member_len: 4,
                },
                Edge {
                    rep_oid: 1,
                    member_oid: 3,
                    rep_len: 11,
                    member_len: 4,
                },
            ],
        );
        cluster(&first, &edges, &volumes).unwrap();
        assert_eq!(read_mapping(&first, 4), [0, 1, 0, 1]);

        // Enter after the first worker completed, with the real shared queue
        // already exhausted and the assignment storage already consumed.
        // This is Source's reachable zero-work path, not a timing simulation.
        let mut late = Job::new(3, 1, 1024, &root, 1).unwrap();
        late.set_round_count(1);
        cluster(&late, &edges, &volumes).unwrap();
        assert_eq!(read_mapping(&late, 4), [0, 1, 2, 3]);
        assert_eq!(
            Atomic::new(late.base_dir(None).join("clustering/closure_finished"))
                .fetch_add(0)
                .unwrap(),
            2
        );
        drop(late);
        drop(first);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn cluster_end_to_end_consumes_native_edge_storage_on_last_round() {
        let root = temp_dir("cluster_e2e");
        let mut job = Job::new(3, 1, 1024, &root, 0).unwrap();
        job.set_round_count(1);
        let volumes = volumes(&root, 4);
        // Source align creates its owner directory before cluster opens the
        // shared alignments/queue. This fixture supplies edge files directly.
        fs::create_dir_all(job.base_dir(None).join("alignments")).unwrap();
        let edges = edge_table(
            &root,
            &[
                Edge {
                    rep_oid: 0,
                    member_oid: 2,
                    rep_len: 12,
                    member_len: 4,
                },
                Edge {
                    rep_oid: 1,
                    member_oid: 2,
                    rep_len: 8,
                    member_len: 4,
                },
                Edge {
                    rep_oid: 1,
                    member_oid: 3,
                    rep_len: 12,
                    member_len: 4,
                },
                Edge {
                    rep_oid: 0,
                    member_oid: 3,
                    rep_len: 12,
                    member_len: 4,
                },
            ],
        );
        assert_eq!(cluster(&job, &edges, &volumes).unwrap(), PathBuf::new());
        assert_eq!(read_mapping(&job, 4), [0, 1, 0, 0]);
        assert!(!edges[0].path.exists());
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bidirectional_prefers_highest_degree_then_lowest_oid() {
        let root = temp_dir("bidirectional");
        let mut job = Job::new(3, 1, 1024, &root, 0).unwrap();
        job.set_round_count(1);
        let volumes = volumes(&root, 4);
        let edges = edge_table(
            &root,
            &[
                Edge {
                    rep_oid: 0,
                    member_oid: 1,
                    rep_len: 4,
                    member_len: 4,
                },
                Edge {
                    rep_oid: 0,
                    member_oid: 2,
                    rep_len: 4,
                    member_len: 4,
                },
                Edge {
                    rep_oid: 3,
                    member_oid: 2,
                    rep_len: 4,
                    member_len: 4,
                },
            ],
        );
        assert_eq!(
            cluster_bidirectional(&job, &edges, &volumes).unwrap(),
            PathBuf::new()
        );
        assert_eq!(read_mapping(&job, 4), [0, 0, 0, 0]);
        assert!(!edges[0].path.exists());
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn parallel_closure_compresses_cross_partition_chains_into_ordered_volume_files() {
        let mut outputs = Vec::new();
        for threads in [1, 4] {
            let root = temp_dir("parallel_closure");
            let job = Job::new(1, 1, 1024, &root, 0).unwrap();
            let volumes = VolumedFile(vec![
                Volume::new("first", 0, 6, 6),
                Volume::new("second", 6, 12, 6),
            ]);
            // Canonical owner, not the compatibility borrowed-slice copy.
            let mapping = identity_mapping(12);
            for oid in 1..12 {
                mapping[oid].store((oid - 1) as OId, Ordering::Relaxed);
            }
            compute_closure_mapping(&job, &volumes, &mapping, threads).unwrap();
            assert!(mapping.iter().all(|oid| oid.load(Ordering::Relaxed) == 0));
            let directory = job.base_dir(None).join("clustering");
            let bytes = [
                fs::read(directory.join("volume0")).unwrap(),
                fs::read(directory.join("volume1")).unwrap(),
            ];
            assert!(bytes.iter().all(|volume| volume == &vec![0; 6 * 8]));
            outputs.push(bytes);
            drop(job);
            fs::remove_dir_all(root).unwrap();
        }
        assert_eq!(outputs[0], outputs[1]);
    }

    #[test]
    fn parallel_assignment_readers_follow_complete_records_to_eof_not_manifest_counts() {
        let root = temp_dir("assignments_eof");
        let job = Job::new(1, 1, 1024, &root, 0).unwrap();
        let volumes = volumes(&root, 6);
        let directory = root.join("assignments");
        fs::create_dir_all(&directory).unwrap();
        let manifest = directory.join("bucket.tsv");
        let mut rows = Vec::new();
        for (index, (assignments, declared)) in [
            (vec![(2u64, 0u64), (3, 2)], 1u64),
            (vec![(4, 3), (5, 4)], 99),
        ]
        .into_iter()
        .enumerate()
        {
            let path = directory.join(format!("volume{index}"));
            let mut codec = CompressedBuffer::new();
            for (member, rep) in assignments {
                codec.write(&member.to_ne_bytes()).unwrap();
                codec.write(&rep.to_ne_bytes()).unwrap();
            }
            // Source fread<T>(1) returns zero for a final partial record and
            // does not publish that incomplete Assignment into the mapping.
            codec.write(&[0xff; 5]).unwrap();
            codec.finish().unwrap();
            fs::write(&path, codec.data()).unwrap();
            rows.extend(super::super::storage::manifest_line(&path, &[declared]).unwrap());
        }
        fs::write(&manifest, rows).unwrap();
        compute_closure_from_assignment_file_with_threads(&job, &manifest, &volumes, 4).unwrap();
        assert_eq!(read_mapping(&job, 6), [0, 1, 0, 0, 0, 0]);
        assert!(!manifest.exists());
        assert!(!directory.exists());
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn member_partitions_emit_one_best_assignment_across_equal_key_boundaries() {
        let mut records = vec![
            Edge {
                rep_oid: 0,
                member_oid: 2,
                rep_len: 12,
                member_len: 4,
            },
            Edge {
                rep_oid: 1,
                member_oid: 2,
                rep_len: 10,
                member_len: 4,
            },
            Edge {
                rep_oid: 2,
                member_oid: 3,
                rep_len: 8,
                member_len: 4,
            },
            Edge {
                rep_oid: 3,
                member_oid: 5,
                rep_len: 8,
                member_len: 4,
            },
            Edge {
                rep_oid: 1,
                member_oid: 6,
                rep_len: 8,
                member_len: 4,
            },
            Edge {
                rep_oid: 6,
                member_oid: 7,
                rep_len: 8,
                member_len: 4,
            },
        ];
        records.extend((0..19).map(|i| Edge {
            rep_oid: i % 2,
            member_oid: 4,
            rep_len: 10 + i as u32,
            member_len: 4,
        }));
        records.sort_unstable();
        let ranges = edge_part_ranges(&records, 4);
        assert_eq!(ranges.len(), 4);
        assert!(ranges.iter().any(std::ops::Range::is_empty));
        let mut group_owners = std::collections::BTreeMap::new();
        for range in ranges {
            let mut previous = None;
            for edge in &records[range] {
                if previous != Some(edge.member_oid) {
                    *group_owners.entry(edge.member_oid).or_insert(0) += 1;
                    previous = Some(edge.member_oid);
                }
            }
        }
        assert_eq!(group_owners.len(), 6);
        assert!(group_owners.values().all(|&count| count == 1));
        for threads in [1, 4] {
            let root = temp_dir("member_partition");
            let mut job = Job::new(1, 1, 1024, &root, 0).unwrap();
            job.set_round_count(1);
            let volumes = volumes(&root, 8);
            fs::create_dir_all(job.base_dir(None).join("alignments")).unwrap();
            let table = edge_table(&root, &records);
            cluster_with_threads(&job, &table, &volumes, threads).unwrap();
            assert_eq!(read_mapping(&job, 8), [0, 1, 0, 0, 0, 0, 1, 1]);
            assert!(!table[0].path.exists());
            drop(job);
            fs::remove_dir_all(root).unwrap();
        }
    }

    #[test]
    fn closure_rejects_cycles_in_corrupt_assignments() {
        let root = temp_dir("cycle");
        let job = Job::new(1, 1, 1024, &root, 0).unwrap();
        let volumes = volumes(&root, 2);
        let mut rep = vec![1, 0];
        assert_eq!(
            compute_closure(&job, &volumes, &mut rep).unwrap_err(),
            "Representative cycle contains OID 1."
        );
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }
}
