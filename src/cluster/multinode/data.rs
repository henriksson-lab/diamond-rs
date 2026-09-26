//! Translation of `diamond/src/cluster/multinode/data.cpp`.
//!
//! The upstream function obtains its thread count globally and consumes
//! `Job` and `SequenceFile` types owned by adjacent translation units. Rust
//! keeps those dependencies explicit through small adapter traits while
//! preserving the file-backed multi-process coordination.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Mutex;

use super::output::Volume;
use crate::basic::sequence::Sequence;
use crate::basic::value::{Letter, OId, SequenceType, AMINO_ACID_ALPHABET};
use crate::data::fasta::fasta_file::{FastaFile, FastaFileConfig};
use crate::util::parallel::{Atomic, FileStack};

/// Explicit replacement for C++ `config.threads_`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GetRepsConfig {
    pub threads: usize,
}

impl Default for GetRepsConfig {
    fn default() -> Self {
        Self { threads: 1 }
    }
}

/// Minimal `Job` surface required by this translation unit.
pub trait RepresentativeJob: Sync {
    fn last_round(&self) -> bool;
    fn base_dir(&self) -> PathBuf;
    fn round(&self) -> i32;
    fn log(&self, message: &str) -> Result<(), String>;
}

/// Format-neutral sequence record corresponding to one `read_seq` result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    pub id: String,
    pub sequence: Vec<Letter>,
}

/// Explicit adapter for C++ `SequenceFile::auto_create` and `read_seq`.
pub trait SequenceSource: Sync {
    fn read_records(&self, path: &Path) -> Result<Vec<SequenceRecord>, String>;
}

/// Concrete source used by multi-node representative volumes, which are
/// emitted as amino-acid FASTA files by the surrounding workflow.
#[derive(Debug, Clone, Copy, Default)]
pub struct FastaSequenceSource;

impl SequenceSource for FastaSequenceSource {
    fn read_records(&self, path: &Path) -> Result<Vec<SequenceRecord>, String> {
        let mut input = FastaFile::open(
            &[path.to_path_buf()],
            FastaFileConfig {
                sequence_type: SequenceType::AminoAcid,
                ..FastaFileConfig::default()
            },
        )?;
        let mut records = Vec::new();
        while let Some(record) = input.read_seq()? {
            records.push(SequenceRecord {
                id: record.id,
                sequence: record.sequence,
            });
        }
        input.close()?;
        Ok(records)
    }
}

/// Convenience entry point using the FASTA source produced by the multi-node
/// workflow.
pub fn get_reps_fasta<J: RepresentativeJob>(
    job: &J,
    volumes: &[Volume],
    config: GetRepsConfig,
) -> Result<String, String> {
    get_reps(job, volumes, config, &FastaSequenceSource)
}

/// C++ `get_reps(Job&, const VolumedFile&)`, with the global thread count and
/// `SequenceFile::auto_create` dependency made explicit.
pub fn get_reps<J, S>(
    job: &J,
    volumes: &[Volume],
    config: GetRepsConfig,
    source: &S,
) -> Result<String, String>
where
    J: RepresentativeJob,
    S: SequenceSource,
{
    if job.last_round() {
        return Ok(String::new());
    }
    if config.threads == 0 {
        return Err("Thread count must be positive.".to_string());
    }

    let base_dir = job.base_dir().join("reps");
    std::fs::create_dir_all(&base_dir).map_err(|error| error.to_string())?;
    let reps_list = FileStack::new(base_dir.join("reps.tsv"));
    let queue = Atomic::new(base_dir.join("queue"));
    let volumes_processed = AtomicU64::new(0);
    let cluster_count = AtomicU64::new(0);
    let stop = AtomicBool::new(false);
    let first_error = Mutex::new(None::<String>);

    std::thread::scope(|scope| {
        for _thread_id in 0..config.threads {
            scope.spawn(|| {
                let result = representative_worker(
                    job,
                    volumes,
                    source,
                    &base_dir,
                    &reps_list,
                    &queue,
                    &volumes_processed,
                    &cluster_count,
                    &stop,
                );
                if let Err(error) = result {
                    stop.store(true, Ordering::Relaxed);
                    let mut slot = first_error.lock().expect("error mutex poisoned");
                    if slot.is_none() {
                        *slot = Some(error);
                    }
                }
            });
        }
    });

    if let Some(error) = first_error
        .lock()
        .map_err(|error| error.to_string())?
        .take()
    {
        return Err(error);
    }

    job.log(&format!(
        "Representatives written: {}",
        cluster_count.load(Ordering::Relaxed)
    ))?;
    let finished = Atomic::new(base_dir.join("finished"));
    finished.fetch_add(volumes_processed.load(Ordering::Relaxed) as i64)?;
    finished.await_value(volumes.len() as i64)?;
    Ok(reps_list.file_name().to_string_lossy().into_owned())
}

#[allow(clippy::too_many_arguments)]
fn representative_worker<J, S>(
    job: &J,
    volumes: &[Volume],
    source: &S,
    base_dir: &Path,
    reps_list: &FileStack,
    queue: &Atomic,
    volumes_processed: &AtomicU64,
    cluster_count: &AtomicU64,
    stop: &AtomicBool,
) -> Result<(), String>
where
    J: RepresentativeJob,
    S: SequenceSource,
{
    loop {
        let volume_index = queue.fetch_add_one()?;
        if stop.load(Ordering::Relaxed) || volume_index >= volumes.len() as i64 {
            return Ok(());
        }
        let volume = &volumes[volume_index as usize];
        job.log(&format!(
            "Writing representatives. Volume={}/{} Records={}",
            volume_index + 1,
            volumes.len(),
            crate::util::string::format(volume.record_count)
        ))?;

        let id_file = job.base_dir().join(format!("rep_ids{volume_index}"));
        let rep_text = std::fs::read_to_string(&id_file)
            .map_err(|_| format!("Error opening file {}", id_file.display()))?;
        let representatives = rep_text
            .split_whitespace()
            .map(|value| {
                value
                    .parse::<OId>()
                    .map_err(|_| format!("Invalid representative oid {value}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let records = source.read_records(&volume.path)?;

        let out_file = base_dir.join(format!("{volume_index}.faa"));
        let output = File::create(&out_file)
            .map_err(|_| format!("Error opening file {}", out_file.display()))?;
        let mut output = BufWriter::new(output);
        let mut count = 0u64;
        let mut representative_index = 0usize;
        let mut file_oid = volume.oid_begin;
        for record in records {
            if stop.load(Ordering::Relaxed) {
                break;
            }
            if job.round() > 0 {
                file_oid = atoll_oid(&record.id);
            }
            if representatives.get(representative_index) == Some(&file_oid) {
                writeln!(
                    output,
                    ">{file_oid}\n{}",
                    Sequence::new(&record.sequence).to_string(AMINO_ACID_ALPHABET)
                )
                .map_err(|error| error.to_string())?;
                count += 1;
                representative_index += 1;
            }
            file_oid = file_oid.wrapping_add(1);
        }
        output.flush().map_err(|error| error.to_string())?;

        if let Some(rep) = representatives.get(representative_index) {
            return Err(format!(
                "Failed to find oid {rep} in file {}",
                volume.path.display()
            ));
        }
        reps_list.push_string(&format!(
            "{}\t{}\t{}\t{}\n",
            out_file.display(),
            count,
            volume.oid_begin,
            volume.oid_end
        ))?;
        volumes_processed.fetch_add(1, Ordering::Relaxed);
        cluster_count.fetch_add(count, Ordering::Relaxed);
    }
}

/// `std::atoll` prefix parsing used for representative FASTA IDs after round
/// zero. Invalid/non-numeric IDs map to zero like C `atoll`.
fn atoll_oid(value: &str) -> OId {
    let value = value.trim_start();
    let (negative, digits) = match value.as_bytes().first() {
        Some(b'-') => (true, &value[1..]),
        Some(b'+') => (false, &value[1..]),
        _ => (false, value),
    };
    let end = digits
        .bytes()
        .position(|byte| !byte.is_ascii_digit())
        .unwrap_or(digits.len());
    if end == 0 {
        return 0;
    }
    let magnitude = digits[..end].parse::<u64>().unwrap_or(u64::MAX);
    if negative {
        0u64.wrapping_sub(magnitude)
    } else {
        magnitude
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    struct TestDir(PathBuf);

    impl TestDir {
        fn new(label: &str) -> Self {
            let nonce = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let path = std::env::temp_dir().join(format!(
                "diamond_multinode_data_{label}_{}_{nonce}",
                std::process::id()
            ));
            std::fs::create_dir_all(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for TestDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    struct TestJob {
        base_dir: PathBuf,
        round: i32,
        last: bool,
        log: Mutex<Vec<String>>,
    }

    impl RepresentativeJob for TestJob {
        fn last_round(&self) -> bool {
            self.last
        }

        fn base_dir(&self) -> PathBuf {
            self.base_dir.clone()
        }

        fn round(&self) -> i32 {
            self.round
        }

        fn log(&self, message: &str) -> Result<(), String> {
            self.log.lock().unwrap().push(message.to_string());
            Ok(())
        }
    }

    fn job(path: &Path, round: i32, last: bool) -> TestJob {
        TestJob {
            base_dir: path.to_path_buf(),
            round,
            last,
            log: Mutex::new(Vec::new()),
        }
    }

    #[test]
    fn last_round_returns_empty_without_creating_files() {
        let dir = TestDir::new("last");
        let output = get_reps_fasta(&job(&dir.0, 2, true), &[], GetRepsConfig::default()).unwrap();
        assert!(output.is_empty());
        assert!(!dir.0.join("reps").exists());
    }

    #[test]
    fn round_zero_writes_selected_oids_and_exact_manifest() {
        let dir = TestDir::new("round0");
        let fasta = dir.0.join("input.faa");
        std::fs::write(&fasta, b">a\nARN\n>b\nDCQ\n>c\nEGH\n").unwrap();
        std::fs::write(dir.0.join("rep_ids0"), b"11\n12\n").unwrap();
        let volume = Volume::new(&fasta, 10, 13, 3);
        let test_job = job(&dir.0, 0, false);

        let manifest = get_reps_fasta(&test_job, &[volume], GetRepsConfig::default()).unwrap();
        let output = dir.0.join("reps/0.faa");
        assert_eq!(
            std::fs::read_to_string(&output).unwrap(),
            ">11\nDCQ\n>12\nEGH\n"
        );
        assert_eq!(
            std::fs::read_to_string(manifest).unwrap(),
            format!("{}\t2\t10\t13\n", output.display())
        );
        assert_eq!(
            test_job.log.lock().unwrap().last().map(String::as_str),
            Some("Representatives written: 2")
        );
    }

    #[test]
    fn later_round_uses_atoll_record_ids() {
        let dir = TestDir::new("later");
        let fasta = dir.0.join("input.faa");
        std::fs::write(&fasta, b">400 old\nARN\n>7\nDCQ\n>999\nEGH\n").unwrap();
        std::fs::write(dir.0.join("rep_ids0"), b"7\n999\n").unwrap();
        let volume = Volume::new(&fasta, 0, 3, 3);

        get_reps_fasta(&job(&dir.0, 1, false), &[volume], GetRepsConfig::default()).unwrap();
        assert_eq!(
            std::fs::read_to_string(dir.0.join("reps/0.faa")).unwrap(),
            ">7\nDCQ\n>999\nEGH\n"
        );
        assert_eq!(atoll_oid("  -2 suffix"), u64::MAX - 1);
        assert_eq!(atoll_oid("not-a-number"), 0);
    }

    #[test]
    fn missing_representative_reports_cpp_error_and_no_manifest() {
        let dir = TestDir::new("missing");
        let fasta = dir.0.join("input.faa");
        std::fs::write(&fasta, b">a\nARN\n").unwrap();
        std::fs::write(dir.0.join("rep_ids0"), b"9\n").unwrap();
        let volume = Volume::new(&fasta, 0, 1, 1);

        let error = get_reps_fasta(&job(&dir.0, 0, false), &[volume], GetRepsConfig::default())
            .unwrap_err();
        assert_eq!(
            error,
            format!("Failed to find oid 9 in file {}", fasta.display())
        );
        assert_eq!(std::fs::read(dir.0.join("reps/reps.tsv")).unwrap(), b"");
    }

    #[test]
    fn multiple_workers_claim_each_volume_once() {
        let dir = TestDir::new("workers");
        let mut volumes = Vec::new();
        for index in 0..4u64 {
            let fasta = dir.0.join(format!("input{index}.faa"));
            std::fs::write(&fasta, b">x\nARN\n").unwrap();
            std::fs::write(dir.0.join(format!("rep_ids{index}")), format!("{index}\n")).unwrap();
            volumes.push(Volume::new(fasta, index, index + 1, 1));
        }
        let manifest = get_reps_fasta(
            &job(&dir.0, 0, false),
            &volumes,
            GetRepsConfig { threads: 3 },
        )
        .unwrap();
        let manifest = std::fs::read_to_string(manifest).unwrap();
        assert_eq!(manifest.lines().count(), 4);
        for index in 0..4 {
            assert!(dir.0.join(format!("reps/{index}.faa")).exists());
        }
    }
}
