//! External-memory clustering from `diamond/src/cluster/external/cluster.cpp`.
//!
//! Cross-process queues in C++ only distribute ownership of buckets and
//! volumes. This port performs those same storage-backed passes
//! deterministically in one process while preserving the native record and
//! manifest formats used by the surrounding external-clustering modules.

use std::fs;
use std::path::{Path, PathBuf};

use super::external::{Assignment, Edge, Job};
use super::pair_table::{Bucket, RadixedTable};
use super::seed_table::{FastaVolumeReader, SequenceVolumeReader, Volume, VolumedFile};
use crate::basic::sequence_utils;
use crate::basic::value::{OId, AMINO_ACID_ALPHABET};
use crate::util::io::{CompressedBuffer, InputFile};

const EDGE_ENCODED_SIZE: usize = 24;
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
    fs::create_dir_all(&output_dir).map_err(|error| error.to_string())?;
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
        fs::write(output_dir.join(format!("volume{index}")), bytes)
            .map_err(|error| error.to_string())?;
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
    let assignment_volumes = read_manifest(assignment_file)?;
    for (path, count) in &assignment_volumes {
        let mut input = open_input(path)?;
        for _ in 0..*count {
            let assignment = Assignment {
                member_oid: input.read_u64().map_err(|error| error.to_string())?,
                rep_oid: input.read_u64().map_err(|error| error.to_string())?,
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
    }
    compute_closure(job, volumes, &mut rep)?;
    remove_manifest_storage(assignment_file, &assignment_volumes)?;
    Ok(())
}

/// C++ file-local `get_reps` using the normal FASTA volume adapter.
pub fn get_reps(job: &Job, volumes: &VolumedFile) -> Result<PathBuf, String> {
    let mut reader = FastaVolumeReader;
    if job.last_round() {
        return Ok(PathBuf::new());
    }
    let base_dir = job.base_dir(None).join("reps");
    fs::create_dir_all(&base_dir).map_err(|error| error.to_string())?;
    let mut manifest = String::new();
    let mut cluster_count = 0_u64;
    for (volume_index, volume) in volumes.0.iter().enumerate() {
        job.log(&format!(
            "Writing representatives. Volume={}/{} Records={}",
            volume_index + 1,
            volumes.0.len(),
            volume.record_count
        ))?;
        let representatives = read_clustering_volume(job, volume_index, volume)?;
        let records = reader.read_volume(&volume.path)?;
        let mut table_oid = volume.oid_begin;
        let mut table_index = 0usize;
        let mut next_file_oid = volume.oid_begin;
        let output_path = base_dir.join(format!("{volume_index}.faa"));
        let mut output = Vec::new();
        let mut count = 0_u64;
        for record in records {
            let file_oid = if job.round_index() > 0 {
                parse_oid_like_atoll(&record.id)?
            } else {
                next_file_oid
            };
            while table_oid < file_oid {
                table_oid += 1;
                table_index += 1;
            }
            if table_oid != file_oid || table_index >= representatives.len() {
                return Err(format!(
                    "Representative table does not contain OID {file_oid}."
                ));
            }
            if representatives[table_index] == file_oid {
                output.extend_from_slice(format!(">{file_oid}\n").as_bytes());
                output.extend_from_slice(
                    sequence_utils::to_string(&record.sequence, AMINO_ACID_ALPHABET).as_bytes(),
                );
                output.push(b'\n');
                count += 1;
            }
            next_file_oid = file_oid
                .checked_add(1)
                .ok_or_else(|| "Sequence OID overflow.".to_owned())?;
            table_oid += 1;
            table_index += 1;
        }
        fs::write(&output_path, output).map_err(|error| error.to_string())?;
        manifest.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            output_path.display(),
            count,
            volume.oid_begin,
            volume.oid_end
        ));
        cluster_count += count;
    }
    job.log(&format!("Representatives written: {cluster_count}"))?;
    let manifest_path = base_dir.join("reps.tsv");
    fs::write(&manifest_path, manifest).map_err(|error| error.to_string())?;
    Ok(manifest_path)
}

/// C++ unidirectional `cluster` pass.
pub fn cluster(job: &Job, edges: &RadixedTable, volumes: &VolumedFile) -> Result<PathBuf, String> {
    let mut assignments = Vec::new();
    for (bucket_index, bucket) in edges.iter().enumerate() {
        let (mut bucket_edges, storage) = read_edge_bucket(bucket)?;
        job.log(&format!(
            "Clustering. Bucket={}/{} Records={} Size={}",
            bucket_index + 1,
            edges.len(),
            bucket_edges.len(),
            bucket_edges.len() * EDGE_ENCODED_SIZE
        ))?;
        bucket_edges.sort_unstable();
        let mut begin = 0usize;
        while begin < bucket_edges.len() {
            let member_oid = bucket_edges[begin].member_oid;
            let edge = bucket_edges[begin];
            if edge.member_len < edge.rep_len
                || (edge.member_len == edge.rep_len && edge.member_oid > edge.rep_oid)
            {
                assignments.push(Assignment {
                    member_oid: edge.member_oid,
                    rep_oid: edge.rep_oid,
                });
            }
            begin += 1;
            while begin < bucket_edges.len() && bucket_edges[begin].member_oid == member_oid {
                begin += 1;
            }
        }
        remove_manifest_storage(&bucket.path, &storage)?;
    }
    let assignment_file = write_assignments(job, &assignments)?;
    compute_closure_from_assignment_file(job, &assignment_file, volumes)?;
    get_reps(job, volumes)
}

/// C++ `cluster_bidirectional` pass.
pub fn cluster_bidirectional(
    job: &Job,
    edges: &RadixedTable,
    volumes: &VolumedFile,
) -> Result<PathBuf, String> {
    job.log("Computing clustering (bi-directional coverage)")?;
    let count = representative_count(volumes)?;
    let mut storage = Vec::new();
    let mut degree = vec![0_u32; count];
    for (bucket_index, bucket) in edges.iter().enumerate() {
        let (bucket_edges, bucket_storage) = read_edge_bucket(bucket)?;
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
        let (bucket_edges, _) = read_edge_bucket(bucket)?;
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
    get_reps(job, volumes)
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
    let bytes = fs::read(
        job.base_dir(None)
            .join("clustering")
            .join(format!("volume{index}")),
    )
    .map_err(|error| error.to_string())?;
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

fn write_assignments(job: &Job, assignments: &[Assignment]) -> Result<PathBuf, String> {
    let directory = job.base_dir(None).join("clustering").join("0");
    fs::create_dir_all(&directory).map_err(|error| error.to_string())?;
    let volume = directory.join(format!("worker_{}_volume_0", job.worker_id()));
    let manifest = directory.join("bucket.tsv");
    if assignments.is_empty() {
        fs::write(&manifest, []).map_err(|error| error.to_string())?;
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
        fs::write(&volume, compressed.data()).map_err(|error| error.to_string())?;
        fs::write(
            &manifest,
            format!("{}\t{}\n", volume.display(), assignments.len()),
        )
        .map_err(|error| error.to_string())?;
    }
    Ok(manifest)
}

fn read_edge_bucket(bucket: &Bucket) -> Result<(Vec<Edge>, Vec<(PathBuf, usize)>), String> {
    let storage = read_manifest(&bucket.path)?;
    let mut edges = Vec::new();
    for (path, count) in &storage {
        let mut input = open_input(path)?;
        for _ in 0..*count {
            edges.push(Edge {
                rep_oid: input.read_u64().map_err(|error| error.to_string())?,
                member_oid: input.read_u64().map_err(|error| error.to_string())?,
                rep_len: input.read_u32().map_err(|error| error.to_string())?,
                member_len: input.read_u32().map_err(|error| error.to_string())?,
            });
        }
    }
    Ok((edges, storage))
}

fn read_manifest(path: &Path) -> Result<Vec<(PathBuf, usize)>, String> {
    let text = fs::read_to_string(path).map_err(|error| error.to_string())?;
    let mut volumes = Vec::new();
    for line in text.lines().filter(|line| !line.trim().is_empty()) {
        let mut fields = line.split_whitespace();
        let volume = fields
            .next()
            .ok_or_else(|| "Format error in VolumedFile".to_owned())?;
        let count = fields
            .next()
            .ok_or_else(|| "Format error in VolumedFile".to_owned())?
            .parse::<usize>()
            .map_err(|_| "Format error in VolumedFile".to_owned())?;
        volumes.push((PathBuf::from(volume), count));
    }
    Ok(volumes)
}

fn open_input(path: &Path) -> Result<InputFile, String> {
    InputFile::new(
        path.to_str()
            .ok_or_else(|| "Non-UTF-8 external-clustering volume path.".to_owned())?,
        0,
    )
    .map_err(|error| error.to_string())
}

fn remove_manifest_storage(path: &Path, volumes: &[(PathBuf, usize)]) -> Result<(), String> {
    for (volume, _) in volumes {
        let _ = fs::remove_file(volume);
    }
    let _ = fs::remove_file(path);
    if let Some(directory) = path.parent() {
        let _ = fs::remove_dir(directory);
    }
    Ok(())
}

fn parse_oid_like_atoll(id: &str) -> Result<OId, String> {
    let bytes = id.as_bytes();
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
            .ok_or_else(|| format!("Sequence OID overflow: {id}"))?;
        index += 1;
    }
    if index == begin {
        return Err(format!("Invalid sequence OID: {id}"));
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
        let mut bucket_edges = read_edge_bucket(&edges[0]).unwrap().0;
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
        assert!(fasta.contains(">0\nARND\n"));
        assert!(fasta.contains(">1\nARND\n"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn cluster_end_to_end_consumes_native_edge_storage_on_last_round() {
        let root = temp_dir("cluster_e2e");
        let mut job = Job::new(3, 1, 1024, &root, 0).unwrap();
        job.set_round_count(1);
        let volumes = volumes(&root, 4);
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
