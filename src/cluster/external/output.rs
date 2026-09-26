//! External-clustering output assembly from
//! `diamond/src/cluster/external/output.cpp`.
//!
//! The C++ thread pools, cross-process queues, and radix-sort workers are
//! represented by deterministic local passes over the same radix files. File
//! formats, ordering, names, round composition, and cleanup remain unchanged.

use std::cmp::Ordering;
use std::fs;
use std::path::{Path, PathBuf};

use super::external::Job;
use super::pair_table::{Bucket, RadixedTable, RADIX_BITS, RADIX_COUNT};
use super::seed_table::VolumedFile;
use crate::basic::value::OId;
use crate::util::io::{CompressedBuffer, InputFile};

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
        let bytes =
            fs::read(&path).map_err(|_| format!("Error opening file: {}", path.display()))?;
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
    let mut bytes = Vec::new();
    let mut clusters = 0;
    for (member, &rep) in merged.iter().enumerate() {
        if rep == member as OId {
            clusters += 1;
        }
        bytes.extend_from_slice(format!("{rep}\t{member}\n").as_bytes());
    }
    fs::write(output_file, bytes)
        .map_err(|_| format!("Error opening file: {}", output_file.display()))?;
    Ok(clusters)
}

pub fn output_accs_round1(
    job: &Job,
    merged: &[OId],
    volumes: &VolumedFile,
) -> Result<RadixedTable, String> {
    if volumes.0.len() < job.volumes {
        return Err(format!(
            "Expected {} database volumes, found {}.",
            job.volumes,
            volumes.0.len()
        ));
    }
    let base_dir = job.root_dir().join("output");
    fs::create_dir_all(&base_dir).map_err(|error| error.to_string())?;
    let shift = bit_length(job.max_oid).saturating_sub(RADIX_BITS);
    let mut mappings = vec![Vec::<AccMapping>::new(); RADIX_COUNT];
    for (index, volume) in volumes.0.iter().take(job.volumes).enumerate() {
        let path = job
            .root_dir()
            .join("accessions")
            .join(format!("{index}.txt"));
        let text = match fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => String::new(),
            Err(error) => return Err(error.to_string()),
        };
        let mut oid = volume.oid_begin;
        for accession in text.lines() {
            let rep = *merged
                .get(oid as usize)
                .ok_or_else(|| format!("Accession OID out of range: {oid}"))?;
            let radix = ((rep >> shift) & (RADIX_COUNT as u64 - 1)) as usize;
            mappings[radix].push(AccMapping::new(rep, oid, accession.to_owned()));
            oid += 1;
        }
    }

    let mut buckets = Vec::with_capacity(RADIX_COUNT);
    for (radix, entries) in mappings.into_iter().enumerate() {
        let directory = base_dir.join(radix.to_string());
        fs::create_dir_all(&directory).map_err(|error| error.to_string())?;
        let manifest = directory.join("bucket.tsv");
        if entries.is_empty() {
            fs::write(&manifest, []).map_err(|error| error.to_string())?;
        } else {
            let volume = directory.join(format!("worker_{}_volume_0", job.worker_id()));
            let mut compressed = CompressedBuffer::new();
            for entry in &entries {
                let mut bytes = Vec::new();
                entry.serialize(&mut bytes);
                compressed
                    .write(&bytes)
                    .map_err(|error| error.to_string())?;
            }
            compressed.finish().map_err(|error| error.to_string())?;
            fs::write(&volume, compressed.data()).map_err(|error| error.to_string())?;
            fs::write(
                &manifest,
                format!("{}\t{}\n", volume.display(), entries.len()),
            )
            .map_err(|error| error.to_string())?;
        }
        buckets.push(Bucket::new(manifest, Some(entries.len() as u64)));
    }
    Ok(RadixedTable(buckets))
}

pub fn output_accs(
    job: &Job,
    merged: &[OId],
    database: &VolumedFile,
    output_file: &Path,
) -> Result<OId, String> {
    let round1 = output_accs_round1(job, merged, database)?;
    let accession_by_oid = read_accessions(job, database)?;
    let digits = decimal_digits(job.max_oid);
    let mut cluster_count = 0;
    for bucket in round1.iter() {
        let mut entries = read_acc_bucket(bucket)?;
        if entries.is_empty() {
            continue;
        }
        entries.sort_unstable();
        let mut begin = 0;
        while begin < entries.len() {
            let rep = entries[begin].rep;
            let mut end = begin + 1;
            while end < entries.len() && entries[end].rep == rep {
                end += 1;
            }
            let rep_acc = accession_by_oid
                .get(rep as usize)
                .and_then(Option::as_deref)
                .ok_or_else(|| format!("Representative accession not found for OID: {rep}"))?;
            let path = PathBuf::from(format!("{}.{rep:0digits$}", output_file.display()));
            let mut output = Vec::new();
            for entry in &entries[begin..end] {
                output.extend_from_slice(rep_acc.as_bytes());
                output.push(b'\t');
                output.extend_from_slice(entry.member_acc.as_bytes());
                output.push(b'\n');
            }
            fs::write(&path, output)
                .map_err(|_| format!("Error opening file: {}", path.display()))?;
            cluster_count += 1;
            begin = end;
        }
    }
    Ok(cluster_count)
}

fn read_acc_bucket(bucket: &Bucket) -> Result<Vec<AccMapping>, String> {
    let manifest = fs::read_to_string(&bucket.path).map_err(|error| error.to_string())?;
    let mut entries = Vec::new();
    for line in manifest.lines().filter(|line| !line.trim().is_empty()) {
        let mut fields = line.split_whitespace();
        let path = fields
            .next()
            .ok_or_else(|| "Format error in VolumedFile".to_owned())?;
        let count = fields
            .next()
            .ok_or_else(|| "Format error in VolumedFile".to_owned())?
            .parse::<usize>()
            .map_err(|_| "Format error in VolumedFile".to_owned())?;
        let mut input = InputFile::new(path, 0).map_err(|error| error.to_string())?;
        for _ in 0..count {
            entries.push(AccMapping::deserialize(&mut input)?);
        }
    }
    Ok(entries)
}

fn read_accessions(job: &Job, database: &VolumedFile) -> Result<Vec<Option<String>>, String> {
    let mut accessions = vec![None; job.max_oid as usize + 1];
    for (index, volume) in database.0.iter().enumerate() {
        let path = job
            .root_dir()
            .join("accessions")
            .join(format!("{index}.txt"));
        let text = fs::read_to_string(&path)
            .map_err(|_| format!("Error opening file: {}", path.display()))?;
        for (offset, accession) in text.lines().enumerate() {
            let oid = volume.oid_begin as usize + offset;
            let slot = accessions
                .get_mut(oid)
                .ok_or_else(|| format!("Accession OID out of range: {oid}"))?;
            *slot = Some(accession.to_owned());
        }
    }
    Ok(accessions)
}

pub fn output(
    job: &Job,
    volumes: &VolumedFile,
    config: &ExternalOutputConfig,
) -> Result<ExternalOutputSummary, String> {
    job.log("Generating output")?;
    let merged = merge(job)?;
    let cluster_count = if config.oid_output {
        output_oids(&merged, &config.output_file)?
    } else {
        output_accs(job, &merged, volumes, &config.output_file)?
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
        if manifest.exists() {
            let text = fs::read_to_string(&manifest).map_err(|error| error.to_string())?;
            for line in text.lines().filter(|line| !line.trim().is_empty()) {
                if let Some(path) = line.split_whitespace().next() {
                    remove_if_present(Path::new(path))?;
                }
            }
            remove_if_present(&manifest)?;
            if let Some(directory) = manifest.parent() {
                let _ = fs::remove_dir(directory);
            }
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
    let _ = fs::remove_file(path);
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
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn oid_output_is_complete_and_counts_self_representatives() {
        let root = temp_dir("oids");
        let output_file = root.join("clusters.tsv");
        let count = output_oids(&[0, 0, 2, 2, 4], &output_file).unwrap();
        assert_eq!(count, 3);
        assert_eq!(
            fs::read_to_string(output_file).unwrap(),
            "0\t0\n0\t1\n2\t2\n2\t3\n4\t4\n"
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
            fs::read_to_string(root.join("clusters.0")).unwrap(),
            "A\tA\nA\tB\n"
        );
        assert_eq!(
            fs::read_to_string(root.join("clusters.2")).unwrap(),
            "C\tC\nC\tD\n"
        );
        assert_eq!(
            fs::read_to_string(root.join("clusters.4")).unwrap(),
            "E\tE\n"
        );
        let radix = output_accs_round1(&job, &[0, 0, 2, 2, 4], &volumes).unwrap();
        assert_eq!(
            radix
                .iter()
                .map(|bucket| bucket.record_count().unwrap())
                .sum::<u64>(),
            5
        );
        fs::remove_dir_all(root).unwrap();
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
            fs::read_to_string(output_file).unwrap(),
            "0\t0\n0\t1\n1\t2\n1\t3\n1\t4\n"
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
        fs::remove_dir_all(root).unwrap();
    }
}
