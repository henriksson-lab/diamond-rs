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

use super::external::{Assignment, Edge, Job};
use super::pair_table::{Bucket, RadixedTable};
use super::seed_table::{Volume, VolumedFile};
use super::storage::BucketWriter;
use crate::basic::sequence_utils;
use crate::basic::value::{OId, AMINO_ACID_ALPHABET};
#[cfg(test)]
use crate::util::io::CompressedBuffer;

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

fn root_of(rep: &[OId], start: OId) -> Result<OId, String> {
    let mut root = start;
    for _ in 0..rep.len() {
        let next = *rep
            .get(root as usize)
            .ok_or_else(|| format!("Representative OID out of range: {root}"))?;
        if next == root {
            return Ok(root);
        }
        root = next;
    }
    Err(format!("Representative cycle contains OID {start}."))
}

/// C++ `compute_closure(Job&, const VolumedFile&, vector<uint64_t>&)`.
pub fn compute_closure(job: &Job, volumes: &VolumedFile, rep: &mut [OId]) -> Result<(), String> {
    let expected = representative_count(volumes)?;
    if rep.len() < expected {
        return Err(format!(
            "Representative vector has {} entries, expected at least {expected}.",
            rep.len()
        ));
    }
    for oid in 0..expected {
        let root = root_of(rep, rep[oid])?;
        rep[oid] = root;
    }

    let output_dir = job.base_dir(None).join("clustering");
    fs::create_dir_all(source_c_path(&output_dir).as_ref()).map_err(|error| error.to_string())?;
    for (index, volume) in volumes.0.iter().enumerate() {
        let begin = usize::try_from(volume.oid_begin)
            .map_err(|_| "Volume OID range does not fit memory.".to_owned())?;
        let end = usize::try_from(volume.oid_end)
            .map_err(|_| "Volume OID range does not fit memory.".to_owned())?;
        let slice = rep
            .get(begin..end)
            .ok_or_else(|| format!("Invalid OID range {begin}..{end}."))?;
        let mut bytes = Vec::with_capacity(slice.len() * std::mem::size_of::<OId>());
        for &oid in slice {
            bytes.extend_from_slice(&oid.to_ne_bytes());
        }
        let mut output = job.create_output(&output_dir.join(format!("volume{index}")))?;
        job.io(output.write_raw(&bytes))?;
        job.io(output.close())?;
    }
    Ok(())
}

/// C++ `compute_closure(Job&, const string&, const VolumedFile&)`.
pub fn compute_closure_from_assignment_file(
    job: &Job,
    assignment_file: &Path,
    volumes: &VolumedFile,
) -> Result<(), String> {
    job.log("Computing transitive closure")?;
    let mut rep = (0..representative_count(volumes)? as OId).collect::<Vec<_>>();
    let assignment_volumes =
        super::external::read_manifest_owned(assignment_file, &job.io_errors())?;
    for (path, count) in &assignment_volumes {
        let mut input = job.open_input(path, 0)?;
        for _ in 0..*count {
            let assignment = Assignment {
                member_oid: job.io(input.read_u64())?,
                rep_oid: job.io(input.read_u64())?,
            };
            if assignment.rep_oid as usize >= rep.len() {
                return Err(format!(
                    "Assignment representative OID out of range: {}",
                    assignment.rep_oid
                ));
            }
            let member = assignment.member_oid as usize;
            if member >= rep.len() {
                return Err(format!(
                    "Assignment member OID out of range: {}",
                    assignment.member_oid
                ));
            }
            rep[member] = assignment.rep_oid;
        }
        job.io(input.close())?;
    }
    compute_closure(job, volumes, &mut rep)?;
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
    fs::create_dir_all(source_c_path(&base_dir).as_ref()).map_err(|error| error.to_string())?;
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
    let bucket_dir = directory.join("0");
    let mut output = BucketWriter::new_owned(&directory, job.worker_id(), 1, job.io_errors())?;
    let queue = Atomic::new(job.base_dir(None).join("alignments/queue"));
    let mut processed = 0;
    loop {
        let index = usize::try_from(queue.fetch_add_one()?)
            .map_err(|_| "Negative clustering bucket index")?;
        let Some(bucket) = edges.get(index) else {
            break;
        };
        let (mut records, storage) = read_edge_bucket(job, bucket)?;
        job.log(&format!(
            "Clustering. Bucket={}/{} Records={} Size={}",
            index + 1,
            edges.len(),
            crate::util::string::format(records.len() as u64),
            crate::util::string::format((records.len() * EDGE_ENCODED_SIZE) as u64)
        ))?;
        records.sort_unstable();
        let mut begin = 0;
        while begin < records.len() {
            let edge = records[begin];
            if edge.member_len < edge.rep_len
                || (edge.member_len == edge.rep_len && edge.member_oid > edge.rep_oid)
            {
                let mut bytes = [0; 16];
                bytes[..8].copy_from_slice(&edge.member_oid.to_ne_bytes());
                bytes[8..].copy_from_slice(&edge.rep_oid.to_ne_bytes());
                output.write(0, &bytes)?;
            }
            begin += 1;
            while begin < records.len() && records[begin].member_oid == edge.member_oid {
                begin += 1;
            }
        }
        // Source's per-bucket BufferArray workers finish before file.remove,
        // retaining the enclosing FileArray across every claimed bucket.
        output.flush_unpublished()?;
        remove_manifest_storage(&bucket.path, &storage)?;
        processed += 1;
    }
    output.finish()?;
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
        compute_closure_from_assignment_file(job, &manifest, volumes)?;
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
        let (bucket_edges, bucket_storage) = read_edge_bucket(job, bucket)?;
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
    let mut rep = (0..count as OId).collect::<Vec<_>>();
    for (bucket_index, bucket) in edges.iter().enumerate() {
        let (bucket_edges, _) = read_edge_bucket(job, bucket)?;
        job.log(&format!(
            "Assigning reps. Bucket={}/{} Records={} Size={}",
            bucket_index + 1,
            edges.len(),
            bucket_edges.len(),
            bucket_edges.len() * EDGE_ENCODED_SIZE
        ))?;
        for edge in &bucket_edges {
            assign_higher_degree(edge.member_oid, edge.rep_oid, &degree, &mut rep)?;
            assign_higher_degree(edge.rep_oid, edge.member_oid, &degree, &mut rep)?;
        }
    }
    compute_closure(job, volumes, &mut rep)?;
    for (manifest, bucket_storage) in storage {
        remove_manifest_storage(&manifest, &bucket_storage)?;
    }
    finished.fetch_add_one()?;
    get_reps_with_threads(job, volumes, threads)
}

fn assign_higher_degree(
    node: OId,
    candidate: OId,
    degree: &[u32],
    rep: &mut [OId],
) -> Result<(), String> {
    let current = *rep
        .get(node as usize)
        .ok_or_else(|| format!("Edge OID out of range: {node}"))?;
    let candidate_degree = *degree
        .get(candidate as usize)
        .ok_or_else(|| format!("Edge OID out of range: {candidate}"))?;
    let current_degree = *degree
        .get(current as usize)
        .ok_or_else(|| format!("Representative OID out of range: {current}"))?;
    if candidate_degree > current_degree
        || (candidate_degree == current_degree && candidate < current)
    {
        rep[node as usize] = candidate;
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

fn read_edge_bucket(
    job: &Job,
    bucket: &Bucket,
) -> Result<(Vec<Edge>, Vec<(PathBuf, usize)>), String> {
    let storage = super::external::read_manifest_owned(&bucket.path, &job.io_errors())?;
    let mut edges = Vec::new();
    for (path, count) in &storage {
        let mut input = job.open_input(path, 0)?;
        for _ in 0..*count {
            edges.push(Edge {
                rep_oid: job.io(input.read_u64())?,
                member_oid: job.io(input.read_u64())?,
                rep_len: job.io(input.read_u32())?,
                member_len: job.io(input.read_u32())?,
            });
        }
        job.io(input.close())?;
    }
    Ok((edges, storage))
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
        let mut bucket_edges = read_edge_bucket(&job, &edges[0]).unwrap().0;
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
        fs::remove_dir_all(root).unwrap();
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
        fs::remove_dir_all(root).unwrap();
    }
}
