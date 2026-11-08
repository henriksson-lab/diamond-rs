//! External-clustering output assembly from
//! `diamond/src/cluster/external/output.cpp`.
//!
//! File-backed queues retain Source's output phases and its single-pass
//! variable-record radix protocol. Accession bytes remain uninterpreted.

use crate::util::io::source_c_path;
use std::cmp::Ordering;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::{Arc, Mutex};

use super::external::Job;
use super::pair_table::{Bucket, RadixedTable, RADIX_BITS, RADIX_COUNT};
use super::seed_table::VolumedFile;
use super::storage::{source_child_path, BucketWriter, WorkerBucketBuffers};
use crate::basic::value::OId;
use crate::util::io::source_text::{SourceTextInput, TextOutput};
use crate::util::io::InputFile;
use crate::util::parallel::{Atomic, FileStack};

pub const ACC_MAPPING_NIL: OId = OId::MAX;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AccMapping {
    pub rep: OId,
    pub member: OId,
    pub rep_acc: String,
    pub member_acc: String,
}

impl AccMapping {
    pub fn new(rep: OId, member: OId, member_acc: String) -> Self {
        Self {
            rep,
            member,
            rep_acc: String::new(),
            member_acc,
        }
    }

    pub fn empty() -> Self {
        Self {
            rep: ACC_MAPPING_NIL,
            member: ACC_MAPPING_NIL,
            rep_acc: String::new(),
            member_acc: String::new(),
        }
    }

    pub const fn key(&self) -> OId {
        self.rep
    }

    pub fn serialize(&self, output: &mut Vec<u8>) {
        output.extend_from_slice(&self.rep.to_ne_bytes());
        output.extend_from_slice(&self.member.to_ne_bytes());
        output.extend_from_slice(self.rep_acc.as_bytes());
        output.push(0);
        output.extend_from_slice(self.member_acc.as_bytes());
        output.push(0);
    }

    pub fn deserialize(input: &mut InputFile) -> Result<Self, String> {
        Ok(Self {
            rep: input.read_u64().map_err(|error| error.to_string())?,
            member: input.read_u64().map_err(|error| error.to_string())?,
            rep_acc: input.read_string().map_err(|error| error.to_string())?,
            member_acc: input.read_string().map_err(|error| error.to_string())?,
        })
    }
}

impl Ord for AccMapping {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.rep, self.member).cmp(&(other.rep, other.member))
    }
}

impl PartialOrd for AccMapping {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalOutputConfig {
    pub output_file: PathBuf,
    pub oid_output: bool,
    pub threads: usize,
}

impl ExternalOutputConfig {
    pub fn new(output_file: impl Into<PathBuf>) -> Self {
        Self {
            output_file: output_file.into(),
            oid_output: false,
            threads: 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalOutputSummary {
    pub merged: Vec<OId>,
    pub cluster_count: OId,
}

pub fn read_clustering(job: &Job, round: i32) -> Result<Vec<OId>, String> {
    let len = usize::try_from(job.max_oid)
        .ok()
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| "Clustering OID range does not fit memory.".to_owned())?;
    let mut clustering = vec![0; len];
    let mut offset = 0usize;
    for volume in 0..job.volumes {
        let path = job
            .base_dir(Some(round))
            .join("clustering")
            .join(format!("volume{volume}"));
        let bytes = job.read_plain_file(&path)?;
        let record_count = bytes.len() / std::mem::size_of::<OId>();
        let end = offset
            .checked_add(record_count)
            .ok_or_else(|| "Clustering record count overflow.".to_owned())?;
        if end > clustering.len() {
            return Err(format!(
                "Clustering contains more than {} OIDs.",
                clustering.len()
            ));
        }
        for (slot, bytes) in clustering[offset..end]
            .iter_mut()
            .zip(bytes.chunks_exact(8))
        {
            *slot = OId::from_ne_bytes(bytes.try_into().unwrap());
        }
        offset = end;
    }
    Ok(clustering)
}

pub fn merge(job: &Job) -> Result<Vec<OId>, String> {
    let mut inner = read_clustering(job, job.round_index())?;
    for round in (0..job.round_index()).rev() {
        let mut outer = read_clustering(job, round)?;
        for centroid in &mut outer {
            *centroid = *inner
                .get(*centroid as usize)
                .ok_or_else(|| format!("Cluster mapping OID out of range: {centroid}"))?;
        }
        inner = outer;
    }
    Ok(inner)
}

pub fn output_oids(merged: &[OId], output_file: &Path) -> Result<OId, String> {
    output_oids_owned(merged, output_file, &Default::default())
}

fn output_oids_owned(
    merged: &[OId],
    output_file: &Path,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<OId, String> {
    let mut output = TextOutput::open(output_file, errors, b"Error opening file: ", false)?;
    let mut clusters = 0;
    for (member, &rep) in merged.iter().enumerate() {
        if rep == member as OId {
            clusters += 1;
        }
        // Source checks fopen alone; fprintf/fclose failures are ignored.
        output.write(format!("{rep}\t{member}\n").as_bytes());
    }
    Ok(clusters)
}

// The pinned GNU64 Source ABI charges two uint64_t fields and two
// std::pmr::strings (40 bytes each). This is an allocation heuristic, not the
// serialized size or Rust's String/Vec layout. The reference object performs
// the corresponding multiply by 96 in radix_sort<AccMapping>.
const SOURCE_ACC_MAPPING_SIZE: u64 = 96;

#[derive(Debug, Eq, PartialEq)]
struct SourceAccMapping {
    rep: OId,
    member: OId,
    rep_acc: Vec<u8>,
    member_acc: Vec<u8>,
}

impl Ord for SourceAccMapping {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.rep, self.member).cmp(&(other.rep, other.member))
    }
}
impl PartialOrd for SourceAccMapping {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl SourceAccMapping {
    fn serialize(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(18 + self.rep_acc.len() + self.member_acc.len());
        bytes.extend_from_slice(&self.rep.to_ne_bytes());
        bytes.extend_from_slice(&self.member.to_ne_bytes());
        bytes.extend_from_slice(&self.rep_acc);
        bytes.push(0);
        bytes.extend_from_slice(&self.member_acc);
        bytes.push(0);
        bytes
    }
    fn empty() -> Self {
        Self {
            rep: ACC_MAPPING_NIL,
            member: ACC_MAPPING_NIL,
            rep_acc: Vec::new(),
            member_acc: Vec::new(),
        }
    }
    fn deserialize(job: &Job, input: &mut InputFile) -> Result<Self, String> {
        let mut record = Self::empty();
        record.read_into(job, input)?;
        Ok(record)
    }
    fn read_into(&mut self, job: &Job, input: &mut InputFile) -> Result<(), String> {
        fn string(job: &Job, input: &mut InputFile, bytes: &mut Vec<u8>) -> Result<(), String> {
            bytes.clear();
            loop {
                let byte = job.io(input.read_value::<u8>())?;
                if byte == 0 {
                    return Ok(());
                }
                bytes.push(byte);
            }
        }
        self.rep = job.io(input.read_u64())?;
        self.member = job.io(input.read_u64())?;
        string(job, input, &mut self.rep_acc)?;
        string(job, input, &mut self.member_acc)
    }
}

/// Source std::getline removes only LF and does not throw on read failure.
/// Its round-one ifstream is unchecked; round two checks the constructor alone.
fn accession_lines(
    path: &Path,
    errors: &crate::util::io::IoErrorOwner,
    checked_open: bool,
    max_lines: Option<u64>,
    mut should_read: impl FnMut() -> bool,
    mut visit: impl FnMut(&[u8]) -> Result<(), String>,
) -> Result<(), String> {
    let mut input = if checked_open {
        SourceTextInput::open(path, errors, b"Error opening file: ")?
    } else {
        SourceTextInput::open_unchecked(path)
    };
    let mut lines = 0;
    while max_lines.map_or(true, |limit| lines < limit) && should_read() {
        let Some(line) = input.line() else { break };
        visit(&line)?;
        lines += 1;
    }
    Ok(())
}

fn accession_volumes(path: &Path) -> Result<Vec<(PathBuf, u64)>, String> {
    accession_volumes_owned(path, &Default::default())
}

fn accession_volumes_owned(
    path: &Path,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<Vec<(PathBuf, u64)>, String> {
    Ok(
        super::external::read_volumed_file_manifest_owned(path, errors)?
            .0
            .into_iter()
            .map(|volume| (volume.path, volume.record_count))
            .collect(),
    )
}

/// Borrow the actual workflow owners while retaining one local worker state
/// per Source thread. The shared stop/error slots survive every worker join.
fn scoped_output_workers(
    job: &Job,
    count: usize,
    worker: impl Fn(&AtomicBool) -> Result<(), String> + Sync,
) -> Result<(), String> {
    let stop = AtomicBool::new(false);
    let error = Mutex::new(None);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..count {
            let work = &worker;
            let stop = &stop;
            let error = &error;
            match std::thread::Builder::new().spawn_scoped(scope, move || {
                if let Err(message) = work(stop) {
                    error.lock().unwrap().get_or_insert(message);
                    stop.store(true, AtomicOrdering::Relaxed);
                }
            }) {
                Ok(handle) => handles.push(handle),
                Err(failure) => {
                    error
                        .lock()
                        .unwrap()
                        .get_or_insert_with(|| job.capture_io::<()>(Err(failure)).unwrap_err());
                    stop.store(true, AtomicOrdering::Relaxed);
                    break;
                }
            }
        }
        for handle in handles {
            if let Err(panic) = handle.join() {
                std::panic::resume_unwind(panic);
            }
        }
    });
    error.into_inner().unwrap().map_or(Ok(()), Err)
}

pub fn output_accs_round1(
    job: &Job,
    merged: &[OId],
    volumes: &VolumedFile,
) -> Result<RadixedTable, String> {
    output_accs_round1_with_threads(job, merged, volumes, 1)
}

pub fn output_accs_round1_with_threads(
    job: &Job,
    merged: &[OId],
    volumes: &VolumedFile,
    threads: usize,
) -> Result<RadixedTable, String> {
    if volumes.0.len() < job.volumes {
        return Err(format!(
            "Expected {} database volumes, found {}.",
            job.volumes,
            volumes.0.len()
        ));
    }
    let base_dir = source_child_path(&job.root_dir().join("output"), "");
    job.capture_io(crate::util::system::mkdir_path(&base_dir))?;
    let output = Arc::new(Mutex::new(BucketWriter::new_file_array_owned(
        &base_dir,
        job.worker_id(),
        RADIX_COUNT,
        job.io_errors(),
    )?));
    let failure = Arc::new(Mutex::new(None));
    let queue = Atomic::new(base_dir.join("queue_round1"));
    let shift = bit_length(job.max_oid).saturating_sub(RADIX_BITS);
    scoped_output_workers(job, threads, |stop| {
        let mut buffers = WorkerBucketBuffers::new(
            Arc::clone(&output),
            RADIX_COUNT,
            job.io_errors(),
            Arc::clone(&failure),
        );
        // Source's worker monotonic resource retains accession allocations
        // until BufferArray destruction; moving records avoids copying bytes.
        let mut arena = Vec::new();
        while !stop.load(AtomicOrdering::Relaxed) {
            let index = usize::try_from(queue.fetch_add_one()?)
                .map_err(|_| "Negative output volume index")?;
            if index >= job.volumes {
                break;
            }
            job.log(&format!(
                "Building output (round 1). Volume={}/{}",
                index + 1,
                job.volumes
            ))?;
            let path = job
                .root_dir()
                .join("accessions")
                .join(format!("{index}.txt"));
            let mut oid = volumes.0[index].oid_begin;
            accession_lines(
                &path,
                &job.io_errors(),
                false,
                None,
                || true,
                |accession| {
                    let rep = *merged
                        .get(oid as usize)
                        .ok_or_else(|| format!("Accession OID out of range: {oid}"))?;
                    let radix =
                        usize::try_from(rep >> shift).map_err(|_| "Accession radix overflow")?;
                    if radix >= RADIX_COUNT {
                        return Err("Accession radix out of range".into());
                    }
                    let entry = SourceAccMapping {
                        rep,
                        member: oid,
                        rep_acc: Vec::new(),
                        member_acc: accession.to_vec(),
                    };
                    buffers.write(radix, &entry.serialize())?;
                    arena.push(entry);
                    oid += 1;
                    Ok(())
                },
            )?;
        }
        // Source destroys the monotonic resource before BufferArray; its
        // serialized bytes are already owned by the compression buffers.
        drop(arena);
        drop(buffers);
        Ok(())
    })?;
    if let Some(error) = failure.lock().unwrap().take() {
        return Err(error);
    }
    let table = output.lock().unwrap().buckets();
    output.lock().unwrap().finish()?;
    Ok(table)
}

fn radix_cluster_accessions(
    job: &Job,
    volumes: &[(PathBuf, u64)],
    directory: &Path,
    bits_unsorted: u32,
    threads: usize,
) -> Result<RadixedTable, String> {
    let output = Arc::new(Mutex::new(BucketWriter::new_exclusive_file_array_owned(
        directory,
        job.worker_id(),
        RADIX_COUNT,
        job.io_errors(),
    )?));
    let failure = Arc::new(Mutex::new(None));
    let next = AtomicUsize::new(0);
    let shift = bits_unsorted - RADIX_BITS;
    scoped_output_workers(job, threads.min(volumes.len()), |stop| {
        let mut buffers = WorkerBucketBuffers::new(
            Arc::clone(&output),
            RADIX_COUNT,
            job.io_errors(),
            Arc::clone(&failure),
        );
        while !stop.load(AtomicOrdering::Relaxed) {
            let index = next.fetch_add(1, AtomicOrdering::Relaxed);
            let Some((path, count)) = volumes.get(index) else {
                break;
            };
            let mut input = job.open_input(path, 0)?;
            // Source reuses one non-POD object for this physical volume;
            // do not retain a separate full record for every radix write.
            let mut record = SourceAccMapping::empty();
            for _ in 0..*count {
                record.read_into(job, &mut input)?;
                let radix = ((record.rep >> shift) & (RADIX_COUNT as u64 - 1)) as usize;
                buffers.write(radix, &record.serialize())?;
            }
            drop(record);
            job.io(input.close())?;
        }
        Ok(())
    })?;
    if let Some(error) = failure.lock().unwrap().take() {
        return Err(error);
    }
    let mut timer = crate::util::log_stream::TaskTimer::with_msg("Closing the output files", 1);
    let table = output.lock().unwrap().buckets();
    output.lock().unwrap().finish()?;
    drop(output);
    timer.finish();
    job.log(&format!(
        "Radix sorted bucket records={}",
        volumes.iter().map(|(_, n)| *n).sum::<u64>()
    ))?;
    Ok(table)
}

fn radix_sort_accessions(
    job: &Job,
    table: RadixedTable,
    bits_unsorted: u32,
    threads: usize,
) -> Result<RadixedTable, String> {
    if bits_unsorted < RADIX_BITS {
        return Ok(table);
    }
    let first = table.first().ok_or("base_path")?;
    let base = super::external::base_path(&first.path)?;
    let queue = Atomic::new(source_child_path(&base, "radix_sort_queue"));
    let manifest_path = source_child_path(&base, "radix_sort_out");
    let manifest = FileStack::try_new_owned(&manifest_path, job.io_errors())?;
    let mut processed = 0_i64;
    loop {
        let index =
            usize::try_from(queue.fetch_add_one()?).map_err(|_| "Negative radix bucket index")?;
        let Some(bucket) = table.get(index) else {
            break;
        };
        let volumes = accession_volumes_owned(&bucket.path, &job.io_errors())?;
        let count: u64 = volumes.iter().map(|(_, n)| *n).sum();
        let bytes = count.wrapping_mul(SOURCE_ACC_MAPPING_SIZE);
        job.log(&format!(
            "Radix sorting. Bucket={}/{} Records={} Size={}",
            index + 1,
            table.len(),
            crate::util::string::format(count),
            crate::util::string::format(bytes)
        ))?;
        if bytes > job.mem_limit {
            let children = radix_cluster_accessions(
                job,
                &volumes,
                bucket.path.parent().ok_or("base_path")?,
                bits_unsorted,
                threads,
            )?;
            let mut text = Vec::new();
            // Source appends every child, including empty buckets. It neither
            // recurses into skewed children nor removes the old input here.
            for child in children.iter() {
                text.extend_from_slice(&job.io(super::storage::manifest_line(
                    &child.path,
                    &[child.record_count_owned(&job.io_errors())?],
                ))?);
            }
            manifest.push_bytes(&text)?;
        } else if count > 0 {
            manifest.push_bytes(&job.io(super::storage::manifest_line(&bucket.path, &[count]))?)?;
        } else {
            for (path, _) in volumes {
                let _ = fs::remove_file(source_c_path(&(path)).as_ref());
            }
            let _ = fs::remove_file(source_c_path(&bucket.path).as_ref());
            if let Some(parent) = bucket.path.parent() {
                let _ = fs::remove_dir(source_c_path(&(parent)).as_ref());
            }
        }
        processed += 1;
    }
    let finished = Atomic::new(source_child_path(&base, "radix_sort_finished"));
    finished.fetch_add(processed)?;
    finished.await_value(table.len() as i64)?;
    RadixedTable::from_manifest_owned(manifest_path, &job.io_errors())
}

fn read_acc_bucket(
    job: &Job,
    bucket: &Bucket,
    threads: usize,
) -> Result<Vec<SourceAccMapping>, String> {
    let volumes = accession_volumes_owned(&bucket.path, &job.io_errors())?;
    let next = AtomicUsize::new(0);
    let completed = Mutex::new(Vec::new());
    scoped_output_workers(job, threads.min(volumes.len()), |stop| {
        while !stop.load(AtomicOrdering::Relaxed) {
            let index = next.fetch_add(1, AtomicOrdering::Relaxed);
            let Some((path, count)) = volumes.get(index) else {
                break;
            };
            let mut input = job.open_input(path, 0)?;
            let mut records = Vec::new();
            for _ in 0..*count {
                records.push(SourceAccMapping::deserialize(job, &mut input)?);
            }
            job.io(input.close())?;
            // Source splices each completed volume's list under its mutex.
            completed.lock().unwrap().push(records);
        }
        Ok(())
    })?;
    Ok(completed
        .into_inner()
        .unwrap()
        .into_iter()
        .flatten()
        .collect())
}

pub fn output_accs(
    job: &Job,
    merged: &[OId],
    database: &VolumedFile,
    output_file: &Path,
) -> Result<OId, String> {
    output_accs_with_threads(job, merged, database, output_file, 1)
}

pub fn output_accs_with_threads(
    job: &Job,
    merged: &[OId],
    database: &VolumedFile,
    output_file: &Path,
    threads: usize,
) -> Result<OId, String> {
    let round1 = output_accs_round1_with_threads(job, merged, database, threads)?;
    let round1 = radix_sort_accessions(
        job,
        round1,
        bit_length(job.max_oid).saturating_sub(RADIX_BITS),
        threads,
    )?;
    let queue = Atomic::new(job.root_dir().join("output").join("queue_round2"));
    let digits = decimal_digits(job.max_oid);
    let cluster_count = AtomicU64::new(0);
    // Source's outer monotonic resource retains decoded accession storage
    // across buckets, through all volume workers and the final count load.
    // Rust retains the actual record owners without copying their byte data.
    let mut arena = Vec::new();
    loop {
        let index =
            usize::try_from(queue.fetch_add_one()?).map_err(|_| "Negative output bucket index")?;
        let Some(bucket) = round1.get(index) else {
            break;
        };
        let mut entries = read_acc_bucket(job, bucket, threads)?;
        job.log(&format!(
            "Building output (round 2). Bucket={}/{} Records={} Size={}",
            index + 1,
            round1.len(),
            crate::util::string::format(entries.len() as u64),
            crate::util::string::format(
                (entries.len() as u64).wrapping_mul(SOURCE_ACC_MAPPING_SIZE)
            )
        ))?;
        if entries.is_empty() {
            continue;
        }
        entries.sort();
        let oid_begin = entries[0].rep;
        let oid_end = entries
            .last()
            .unwrap()
            .rep
            .checked_add(1)
            .ok_or("Accession OID overflow")?;
        let volumes: Vec<_> = database
            .0
            .iter()
            .enumerate()
            .filter(|(_, volume)| volume.oid_end > oid_begin && volume.oid_begin < oid_end)
            .collect();
        let next = AtomicUsize::new(0);
        scoped_output_workers(job, threads.min(volumes.len()), |stop| {
            let mut my_cluster_count = 0;
            let mut position = 0;
            while !stop.load(AtomicOrdering::Relaxed) {
                let index = next.fetch_add(1, AtomicOrdering::Relaxed);
                let Some(&(volume_index, volume)) = volumes.get(index) else {
                    break;
                };
                // Source keeps one cursor per worker and one output owner per
                // claimed volume, even if it contains multiple representatives.
                let mut output: Option<TextOutput> = None;
                while position < entries.len() && entries[position].rep < volume.oid_begin {
                    position += 1;
                }
                let mut file_oid = volume.oid_begin;
                let path = job
                    .root_dir()
                    .join("accessions")
                    .join(format!("{volume_index}.txt"));
                accession_lines(
                    &path,
                    &job.io_errors(),
                    true,
                    Some(oid_end.saturating_sub(volume.oid_begin)),
                    || !stop.load(AtomicOrdering::Relaxed),
                    |accession| {
                        if position >= entries.len() {
                            file_oid += 1;
                            return Ok(());
                        }
                        if entries[position].rep > file_oid {
                            file_oid += 1;
                            return Ok(());
                        }
                        let begin = position;
                        while position < entries.len() && entries[position].rep == file_oid {
                            if output.is_none() {
                                let mut name = output_file.as_os_str().to_os_string();
                                name.push(format!(".{:0digits$}", entries[position].rep));
                                output =
                                    Some(TextOutput::open_unchecked(&PathBuf::from(name), true));
                            }
                            let file = output.as_mut().unwrap();
                            file.write(accession);
                            file.write(b"\t");
                            file.write(&entries[position].member_acc);
                            file.endl();
                            position += 1;
                        }
                        if position > begin {
                            my_cluster_count += 1;
                        }
                        file_oid += 1;
                        Ok(())
                    },
                )?;
            }
            cluster_count.fetch_add(my_cluster_count, AtomicOrdering::Relaxed);
            Ok(())
        })?;
        arena.push(entries);
    }
    Ok(cluster_count.load(AtomicOrdering::Relaxed))
}

pub fn output(
    job: &Job,
    volumes: &VolumedFile,
    config: &ExternalOutputConfig,
) -> Result<ExternalOutputSummary, String> {
    job.log("Generating output")?;
    let merged = merge(job)?;
    let cluster_count = if config.oid_output {
        output_oids_owned(&merged, &config.output_file, &job.io_errors())?
    } else {
        output_accs_with_threads(job, &merged, volumes, &config.output_file, config.threads)?
    };
    job.log(&format!("Cluster count = {cluster_count}"))?;
    cleanup_round_files(job)?;
    Ok(ExternalOutputSummary {
        merged,
        cluster_count,
    })
}

fn cleanup_round_files(job: &Job) -> Result<(), String> {
    for volume in 0..job.volumes {
        remove_if_present(
            &job.base_dir(Some(job.round_index()))
                .join("clustering")
                .join(format!("volume{volume}")),
        )?;
    }
    for round in (0..job.round_index()).rev() {
        let manifest = job.base_dir(Some(round)).join("reps").join("reps.tsv");
        let reps = super::external::read_volumed_file_manifest_owned(&manifest, &job.io_errors())?;
        for volume in reps.0 {
            remove_if_present(&volume.path)?;
        }
        remove_if_present(&manifest)?;
        if let Some(directory) = manifest.parent() {
            let _ = fs::remove_dir(source_c_path(directory).as_ref());
        }
        for volume in 0..job.volumes {
            remove_if_present(
                &job.base_dir(Some(round))
                    .join("clustering")
                    .join(format!("volume{volume}")),
            )?;
        }
    }
    Ok(())
}

fn remove_if_present(path: &Path) -> Result<(), String> {
    // C++ intentionally ignores `remove(3)` return values during cleanup.
    let _ = fs::remove_file(source_c_path(&(path)).as_ref());
    Ok(())
}

fn bit_length(value: u64) -> u32 {
    u64::BITS - value.leading_zeros()
}

fn decimal_digits(value: u64) -> usize {
    value.to_string().len()
}

#[cfg(test)]
mod tests {
    use super::super::seed_table::Volume;
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "diamond_external_output_{name}_{}_{}",
            std::process::id(),
            nonce
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    fn write_oids(path: &Path, values: &[OId]) {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).unwrap();
        }
        let mut bytes = Vec::new();
        for value in values {
            bytes.extend_from_slice(&value.to_ne_bytes());
        }
        fs::write(path, bytes).unwrap();
    }

    fn write_round(job: &Job, round: i32, volumes: &[&[OId]]) {
        for (index, values) in volumes.iter().enumerate() {
            write_oids(
                &job.base_dir(Some(round))
                    .join("clustering")
                    .join(format!("volume{index}")),
                values,
            );
        }
    }

    fn fixture_volumes(root: &Path) -> VolumedFile {
        VolumedFile(vec![
            Volume::new(root.join("db0"), 0, 3, 3),
            Volume::new(root.join("db1"), 3, 5, 2),
        ])
    }

    #[test]
    fn acc_mapping_native_serialization_sort_and_defaults_are_exact() {
        let mapping = AccMapping::new(7, 9, "member accession".to_owned());
        let mut bytes = Vec::new();
        mapping.serialize(&mut bytes);
        assert_eq!(&bytes[..8], &7_u64.to_ne_bytes());
        assert_eq!(&bytes[8..16], &9_u64.to_ne_bytes());
        assert_eq!(&bytes[16..18], &[0, b'm']);
        assert_eq!(bytes.last(), Some(&0));
        assert_eq!(AccMapping::empty().key(), OId::MAX);
        let mut values = [
            AccMapping::new(2, 9, "x".to_owned()),
            AccMapping::new(1, 8, "y".to_owned()),
            AccMapping::new(2, 3, "z".to_owned()),
        ];
        values.sort_unstable();
        assert_eq!(
            values.iter().map(|v| (v.rep, v.member)).collect::<Vec<_>>(),
            [(1, 8), (2, 3), (2, 9)]
        );
    }

    #[test]
    fn read_and_merge_compose_rounds_in_database_order() {
        let root = temp_dir("merge");
        let mut job = Job::new(4, 2, 1 << 20, &root, 0).unwrap();
        write_round(&job, 0, &[&[0, 0, 2], &[2, 4]]);
        job.next_round().unwrap();
        write_round(&job, 1, &[&[0, 0], &[1, 1, 1]]);
        assert_eq!(read_clustering(&job, 0).unwrap(), [0, 0, 2, 2, 4]);
        assert_eq!(merge(&job).unwrap(), [0, 0, 1, 1, 1]);
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn oid_output_is_complete_and_counts_self_representatives() {
        let root = temp_dir("oids");
        let output_file = root.join("clusters.tsv");
        let count = output_oids(&[0, 0, 2, 2, 4], &output_file).unwrap();
        assert_eq!(count, 3);
        assert_eq!(
            fs::read(output_file).unwrap(),
            crate::util::io::source_text::expected_text(b"0\t0\n0\t1\n2\t2\n2\t3\n4\t4\n")
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn accession_output_radixes_sorts_and_names_each_cluster() {
        let root = temp_dir("accessions");
        let job = Job::new(4, 2, 1 << 20, &root, 3).unwrap();
        fs::create_dir_all(root.join("accessions")).unwrap();
        fs::write(root.join("accessions/0.txt"), "A\nB\nC\n").unwrap();
        fs::write(root.join("accessions/1.txt"), "D\nE\n").unwrap();
        let volumes = fixture_volumes(&root);
        let prefix = root.join("clusters");
        let count = output_accs(&job, &[0, 0, 2, 2, 4], &volumes, &prefix).unwrap();
        assert_eq!(count, 3);
        assert_eq!(
            fs::read(root.join("clusters.0")).unwrap(),
            crate::util::io::source_text::expected_text(b"A\tA\nA\tB\n")
        );
        assert_eq!(
            fs::read(root.join("clusters.2")).unwrap(),
            crate::util::io::source_text::expected_text(b"C\tC\nC\tD\n")
        );
        assert_eq!(
            fs::read(root.join("clusters.4")).unwrap(),
            crate::util::io::source_text::expected_text(b"E\tE\n")
        );
        // The actual Source round-one queue is persistent. Inspect the
        // completed manifests instead of invoking that phase twice in a job.
        assert_eq!(
            (0..RADIX_COUNT)
                .map(
                    |radix| accession_volumes(&root.join(format!("output/{radix}/bucket.tsv")))
                        .unwrap()
                        .iter()
                        .map(|(_, count)| *count)
                        .sum::<u64>()
                )
                .sum::<u64>(),
            5
        );
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn same_bucket_representatives_share_one_source_volume_output_owner() {
        let root = temp_dir("same_bucket_owner");
        let count = 512;
        let job = Job::new(count - 1, 1, 1_000_000_000, &root, 0).unwrap();
        let volumes = VolumedFile(vec![Volume::new(root.join("unused"), 0, count, count)]);
        fs::create_dir_all(root.join("accessions")).unwrap();
        let ids: Vec<Vec<u8>> = (0..count)
            .map(|i| {
                let mut id = format!("A{i}").into_bytes();
                if i == 1 {
                    id.push(0xff);
                }
                id
            })
            .collect();
        let input: Vec<u8> = ids
            .iter()
            .flat_map(|id| id.iter().copied().chain([b'\n']))
            .collect();
        fs::write(root.join("accessions/0.txt"), input).unwrap();
        let merged: Vec<_> = (0..count).map(|i| i % 2).collect();
        assert_eq!(
            output_accs(&job, &merged, &volumes, &root.join("clusters")).unwrap(),
            2
        );
        assert!(!root.join("clusters.001").exists());
        let mut expected = Vec::new();
        for rep in 0..2 {
            for member in (rep..count as usize).step_by(2) {
                expected.extend_from_slice(&ids[rep]);
                expected.push(b'\t');
                expected.extend_from_slice(&ids[member]);
                expected.push(b'\n');
            }
        }
        assert_eq!(
            fs::read(root.join("clusters.000")).unwrap(),
            crate::util::io::source_text::expected_text(&expected)
        );
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn source_nonpod_memory_charge_adds_one_radix_pass_and_changes_owner_boundaries() {
        let count = 65_536_u64;
        assert_eq!(count * SOURCE_ACC_MAPPING_SIZE, 6_291_456);
        let ids: Vec<u8> = (0..count)
            .flat_map(|i| format!("A{i}\n").into_bytes())
            .collect();
        let merged: Vec<_> = (0..count).map(|i| i % 2).collect();
        for (name, memory, split) in [("large", 1_000_000_000, false), ("small", 5_000_000, true)] {
            let root = temp_dir(name);
            let job = Job::new(count - 1, 1, memory, &root, 0).unwrap();
            let volumes = VolumedFile(vec![Volume::new(root.join("unused"), 0, count, count)]);
            fs::create_dir_all(root.join("accessions")).unwrap();
            fs::write(root.join("accessions/0.txt"), &ids).unwrap();
            assert_eq!(
                output_accs(&job, &merged, &volumes, &root.join("clusters")).unwrap(),
                2
            );
            assert_eq!(root.join("clusters.00001").exists(), split);
            let first = fs::read(root.join("clusters.00000")).unwrap();
            let second = if split {
                fs::read(root.join("clusters.00001")).unwrap()
            } else {
                Vec::new()
            };
            assert_eq!(
                first.iter().filter(|&&b| b == b'\n').count()
                    + second.iter().filter(|&&b| b == b'\n').count(),
                count as usize
            );
            let manifest = fs::read_to_string(root.join("output/radix_sort_out")).unwrap();
            assert_eq!(
                manifest.lines().count(),
                if split { RADIX_COUNT } else { 1 }
            );
            assert_eq!(root.join("output/0/0/bucket.tsv").exists(), split);
            assert!(root.join("output/0/bucket.tsv").exists()); // Source retains spill input.
            drop(job);
            fs::remove_dir_all(root).unwrap();
        }
    }

    #[test]
    fn parallel_accession_workers_preserve_raw_titles_and_per_volume_output_owners() {
        let root = temp_dir("parallel_accessions");
        let job = Job::new(7, 4, 1 << 20, &root, 0).unwrap();
        fs::create_dir_all(root.join("accessions")).unwrap();
        let mut volumes = Vec::new();
        for volume in 0..4 {
            let mut text = format!("V{volume}a").into_bytes();
            text.push(0xff);
            text.extend_from_slice(format!("\nV{volume}b\n").as_bytes());
            // Use the actual Source text bytes for its accession ifstream.
            fs::write(
                root.join(format!("accessions/{volume}.txt")),
                crate::util::io::source_text::expected_text(&text),
            )
            .unwrap();
            volumes.push(Volume::new(
                root.join(format!("db{volume}")),
                volume * 2,
                volume * 2 + 2,
                2,
            ));
        }
        let prefix = root.join("clusters");
        assert_eq!(
            output_accs_with_threads(
                &job,
                &[0, 0, 2, 2, 4, 4, 6, 6],
                &VolumedFile(volumes),
                &prefix,
                4
            )
            .unwrap(),
            4
        );
        for volume in 0..4 {
            let mut rep = format!("V{volume}a").into_bytes();
            rep.push(0xff);
            let mut expected = rep.clone();
            expected.push(b'\t');
            expected.extend_from_slice(&rep);
            expected.push(b'\n');
            expected.extend_from_slice(&rep);
            expected.extend_from_slice(format!("\tV{volume}b\n").as_bytes());
            assert_eq!(
                fs::read(root.join(format!("clusters.{}", volume * 2))).unwrap(),
                crate::util::io::source_text::expected_text(&expected)
            );
        }
        // Each claimed original volume is logged once despite four readers.
        let log = fs::read_to_string(root.join("diamond_job.log")).unwrap();
        for volume in 1..=4 {
            assert_eq!(
                log.matches(&format!("Building output (round 1). Volume={volume}/4"))
                    .count(),
                1
            );
        }
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn nonpod_radix_workers_reuse_string_buffers_and_publish_exact_shared_counts() {
        let root = temp_dir("parallel_nonpod_radix");
        let job = Job::new(7, 4, 1, root.join("job"), 0).unwrap();
        let mut volumes = Vec::new();
        let mut expected = Vec::new();
        for volume in 0..4 {
            let path = root.join(format!("input{volume}"));
            let mut bytes = Vec::new();
            for record in 0..3 {
                let mapping = SourceAccMapping {
                    rep: volume + 1,
                    member: volume * 3 + record,
                    rep_acc: vec![b'R'; 64 - record as usize * 20],
                    member_acc: vec![0xfe; 48 + record as usize * 20],
                };
                let encoded = mapping.serialize();
                expected.push(encoded.clone());
                bytes.extend_from_slice(&encoded);
            }
            fs::write(&path, bytes).unwrap();
            volumes.push((path, 3));
        }
        fs::create_dir(root.join("out")).unwrap();
        let children =
            radix_cluster_accessions(&job, &volumes, &root.join("out"), RADIX_BITS, 4).unwrap();
        assert_eq!(
            children
                .iter()
                .map(|bucket| bucket.records.unwrap())
                .sum::<u64>(),
            12
        );
        let mut actual = Vec::new();
        for (radix, bucket) in children.iter().enumerate() {
            for mapping in read_acc_bucket(&job, bucket, 4).unwrap() {
                assert_eq!((mapping.rep & 255) as usize, radix);
                actual.push(mapping.serialize());
            }
        }
        actual.sort();
        expected.sort();
        assert_eq!(actual, expected);
        assert!(volumes.iter().all(|(path, _)| path.exists()));
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn cleanup_uses_checked_source_manifest_grammar_after_current_round_removal() {
        for blank in [false, true] {
            let root = temp_dir("checked_cleanup");
            let mut job = Job::new(0, 1, 1 << 20, &root, 0).unwrap();
            write_round(&job, 0, &[&[0]]);
            job.next_round().unwrap();
            write_round(&job, 1, &[&[0]]);
            let manifest = job.base_dir(Some(0)).join("reps").join("reps.tsv");
            if blank {
                fs::create_dir_all(manifest.parent().unwrap()).unwrap();
                fs::write(&manifest, b"\n").unwrap();
            }
            let error = cleanup_round_files(&job).unwrap_err();
            if blank {
                assert_eq!(error, "Format error in VolumedFile");
            } else {
                let restored = job.io_errors().finish::<()>(Err(error)).unwrap_err();
                let mut expected = b"Error opening file ".to_vec();
                expected.extend_from_slice(
                    &crate::util::io::source_os_bytes(manifest.as_os_str()).unwrap(),
                );
                assert_eq!(
                    crate::util::io::RawDiagnosticError::find(&restored),
                    Some(expected.as_slice())
                );
            }
            assert!(!job.base_dir(Some(1)).join("clustering/volume0").exists());
            assert!(job.base_dir(Some(0)).join("clustering/volume0").exists());
            drop(job);
            fs::remove_dir_all(root).unwrap();
        }
    }

    #[test]
    fn output_logs_merges_and_removes_all_round_inputs() {
        let root = temp_dir("workflow");
        let mut job = Job::new(4, 2, 1 << 20, &root, 0).unwrap();
        write_round(&job, 0, &[&[0, 0, 2], &[2, 4]]);
        let reps_dir = job.base_dir(Some(0)).join("reps");
        fs::create_dir_all(&reps_dir).unwrap();
        let rep_volume = reps_dir.join("volume0");
        fs::write(&rep_volume, b"data").unwrap();
        fs::write(
            reps_dir.join("reps.tsv"),
            format!("{}\t1\t0\t1\n", rep_volume.display()),
        )
        .unwrap();
        job.next_round().unwrap();
        write_round(&job, 1, &[&[0, 0], &[1, 1, 1]]);
        let output_file = root.join("result.tsv");
        let summary = output(
            &job,
            &VolumedFile::default(),
            &ExternalOutputConfig {
                output_file: output_file.clone(),
                oid_output: true,
                threads: 2,
            },
        )
        .unwrap();
        assert_eq!(summary.merged, [0, 0, 1, 1, 1]);
        assert_eq!(summary.cluster_count, 1);
        assert_eq!(
            fs::read(output_file).unwrap(),
            crate::util::io::source_text::expected_text(b"0\t0\n0\t1\n1\t2\n1\t3\n1\t4\n")
        );
        assert!(!rep_volume.exists());
        assert!(!reps_dir.join("reps.tsv").exists());
        for round in 0..=1 {
            for volume in 0..2 {
                assert!(!job
                    .base_dir(Some(round))
                    .join("clustering")
                    .join(format!("volume{volume}"))
                    .exists());
            }
        }
        let log = fs::read_to_string(root.join("diamond_job.log")).unwrap();
        assert!(log.contains("Generating output"));
        assert!(log.contains("Cluster count = 1"));
        drop(job);
        fs::remove_dir_all(root).unwrap();
    }
}
