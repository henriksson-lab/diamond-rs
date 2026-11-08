//! Translation of `diamond/src/cluster/multinode/data.cpp`.
//!
//! The upstream function obtains its thread count globally and consumes
//! `Job` and `SequenceFile` types owned by adjacent translation units. Rust
//! keeps those dependencies explicit through small adapter traits while
//! preserving the file-backed multi-process coordination.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Mutex;

use super::output::Volume;
use crate::basic::sequence::Sequence;
use crate::basic::value::{Letter, OId, SequenceType, AMINO_ACID_ALPHABET};
use crate::data::block::Block;
use crate::data::fasta::fasta_file::{FastaFile, FastaFileConfig};
use crate::data::sequence_file::{ProteinInputProbe, SequenceFileFlags};
use crate::util::io::source_text::{path_error, SourceTextInput, TextOutput};
use crate::util::io::{source_c_path, IoErrorOwner};
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
    fn io_errors(&self) -> IoErrorOwner {
        IoErrorOwner::default()
    }
    fn last_round(&self) -> bool;
    fn base_dir(&self) -> PathBuf;
    fn round(&self) -> i32;
    fn log(&self, message: &str) -> Result<(), String>;
}

/// Format-neutral sequence record corresponding to one `read_seq` result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord<T = String> {
    pub id: T,
    pub sequence: Vec<Letter>,
}

pub type ByteSequenceRecord = SequenceRecord<Vec<u8>>;

/// An open sequence owner: construction, record iteration and checked close
/// remain separate Source operations. Legacy Vec adapters have no file owner.
pub trait SequenceReader {
    fn read_record(&mut self) -> Result<Option<ByteSequenceRecord>, String>;
    fn close(&mut self) -> Result<(), String>;
}
struct VecSequenceReader(std::vec::IntoIter<ByteSequenceRecord>);
impl SequenceReader for VecSequenceReader {
    fn read_record(&mut self) -> Result<Option<ByteSequenceRecord>, String> {
        Ok(self.0.next())
    }
    fn close(&mut self) -> Result<(), String> {
        Ok(())
    }
}
struct FastaSequenceReader {
    input: FastaFile,
    errors: IoErrorOwner,
}
impl SequenceReader for FastaSequenceReader {
    fn read_record(&mut self) -> Result<Option<ByteSequenceRecord>, String> {
        self.errors
            .capture(self.input.read_seq_bytes_io())
            .map(|record| {
                record.map(|record| SequenceRecord {
                    id: record.id,
                    sequence: record.sequence,
                })
            })
    }
    fn close(&mut self) -> Result<(), String> {
        self.errors.capture(self.input.close_io())
    }
}

/// Metadata and indexed operations used by the historical Source FASTA owner.
pub trait IndexedSequenceReader: SequenceReader {
    fn letters(&self) -> u64;
    fn partition(&self, max_letters: u64, start_size: u64) -> Result<(Vec<OId>, Vec<u64>), String>;
    fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String>;
    fn tell_seq(&self) -> OId;
    fn load_seqs(
        &mut self,
        max_letters: i64,
        max_seqs: OId,
        flags: SequenceFileFlags,
        ignore_warnings: bool,
    ) -> Result<Block, String> {
        let mut block = Block::new();
        let mut oid = self.tell_seq();
        let mut probe = ProteinInputProbe::new(ignore_warnings || oid != 0);
        let (mut count, mut letters) = (0u64, 0i64);
        loop {
            let Some(record) = self.read_record()? else {
                break;
            };
            if record.sequence.is_empty() {
                continue;
            }
            letters = letters.wrapping_add(
                block.push_back_bytes(
                    &record.sequence,
                    flags
                        .contains(SequenceFileFlags::TITLES)
                        .then_some(record.id.as_slice()),
                    None,
                    oid,
                    SequenceType::AminoAcid,
                    1,
                    true,
                )?,
            );
            oid = oid.wrapping_add(1);
            count = count.wrapping_add(1);
            probe.observe(&record.sequence)?;
            if if max_seqs > 0 {
                count >= max_seqs
            } else {
                letters >= max_letters
            } {
                break;
            }
        }
        probe.finish_block()?;
        Ok(block)
    }
}
struct VecIndexedSequenceReader {
    records: Vec<ByteSequenceRecord>,
    at: usize,
}
impl SequenceReader for VecIndexedSequenceReader {
    fn read_record(&mut self) -> Result<Option<ByteSequenceRecord>, String> {
        let record = self.records.get(self.at).cloned();
        self.at = self.at.saturating_add(1);
        Ok(record)
    }
    fn close(&mut self) -> Result<(), String> {
        Ok(())
    }
}
impl IndexedSequenceReader for VecIndexedSequenceReader {
    fn letters(&self) -> u64 {
        self.records.iter().fold(0u64, |n, record| {
            n.wrapping_add(record.sequence.len() as u64)
        })
    }
    fn partition(&self, max_letters: u64, start_size: u64) -> Result<(Vec<OId>, Vec<u64>), String> {
        let lengths = self
            .records
            .iter()
            .map(|record| record.sequence.len() as crate::basic::value::Loc)
            .collect::<Vec<_>>();
        Ok(crate::data::sequence_file::partition(
            &lengths,
            max_letters,
            start_size,
        ))
    }
    fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String> {
        self.at = usize::try_from(oid).map_err(|_| "OId out of bounds.")?;
        Ok(())
    }
    fn tell_seq(&self) -> OId {
        self.at as OId
    }
}
impl IndexedSequenceReader for FastaSequenceReader {
    fn letters(&self) -> u64 {
        self.input.letters()
    }
    fn partition(&self, max_letters: u64, start_size: u64) -> Result<(Vec<OId>, Vec<u64>), String> {
        let lengths = (0..self.input.sequence_count())
            .map(|oid| self.input.seq_length(oid))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(crate::data::sequence_file::partition(
            &lengths,
            max_letters,
            start_size,
        ))
    }
    fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String> {
        self.errors.capture(self.input.set_seqinfo_ptr_io(oid))
    }
    fn tell_seq(&self) -> OId {
        self.input.tell_seq()
    }
}

/// Explicit adapter for C++ `SequenceFile::auto_create` and `read_seq`.
pub trait SequenceSource: Sync {
    /// In-memory compatibility sources have no physical index to materialize.
    fn index_volume(
        &self,
        _path: &Path,
        _index: &Path,
        _errors: &IoErrorOwner,
    ) -> Result<(), String> {
        Ok(())
    }
    fn open_indexed_records(
        &self,
        path: &Path,
        _index: &Path,
        _flags: SequenceFileFlags,
        errors: &IoErrorOwner,
    ) -> Result<Box<dyn IndexedSequenceReader>, String> {
        Ok(Box::new(VecIndexedSequenceReader {
            records: self.read_records_bytes_with_errors(path, errors)?,
            at: 0,
        }))
    }
    /// Compatibility implementation for in-memory custom sources. Concrete
    /// file-backed sources override this to retain their actual input owner.
    fn open_records(
        &self,
        path: &Path,
        errors: &IoErrorOwner,
    ) -> Result<Box<dyn SequenceReader>, String> {
        Ok(Box::new(VecSequenceReader(
            self.read_records_bytes_with_errors(path, errors)?
                .into_iter(),
        )))
    }
    fn read_records_bytes_with_errors(
        &self,
        path: &Path,
        _errors: &IoErrorOwner,
    ) -> Result<Vec<ByteSequenceRecord>, String> {
        self.read_records_bytes(path)
    }
    fn read_records(&self, path: &Path) -> Result<Vec<SequenceRecord>, String>;
    /// Checked legacy adapters can supply UTF-8 records; native Source owners
    /// override this entry point to retain arbitrary std::string title bytes.
    fn read_records_bytes(&self, path: &Path) -> Result<Vec<ByteSequenceRecord>, String> {
        self.read_records(path).map(|records| {
            records
                .into_iter()
                .map(|record| SequenceRecord {
                    id: record.id.into_bytes(),
                    sequence: record.sequence,
                })
                .collect()
        })
    }
}

/// Concrete source used by multi-node representative volumes, which are
/// emitted as amino-acid FASTA files by the surrounding workflow.
#[derive(Debug, Clone, Copy, Default)]
pub struct FastaSequenceSource;

impl FastaSequenceSource {
    pub fn read_records_bytes_io(&self, path: &Path) -> std::io::Result<Vec<ByteSequenceRecord>> {
        let errors = IoErrorOwner::default();
        errors.finish(self.read_records_bytes_with_errors(path, &errors))
    }
}

impl SequenceSource for FastaSequenceSource {
    fn index_volume(&self, path: &Path, index: &Path, errors: &IoErrorOwner) -> Result<(), String> {
        errors.capture(FastaFile::index_io(path, index))
    }
    fn open_indexed_records(
        &self,
        path: &Path,
        index: &Path,
        flags: SequenceFileFlags,
        errors: &IoErrorOwner,
    ) -> Result<Box<dyn IndexedSequenceReader>, String> {
        let mut input = errors.capture(FastaFile::open_source_with_buffer(
            &[path.to_path_buf()],
            FastaFileConfig {
                flags,
                index_file: Some(index.to_path_buf()),
                sequence_type: SequenceType::AminoAcid,
            },
            65_536,
        ))?;
        input.set_io_errors(errors.clone());
        Ok(Box::new(FastaSequenceReader {
            input,
            errors: errors.clone(),
        }))
    }
    fn open_records(
        &self,
        path: &Path,
        errors: &IoErrorOwner,
    ) -> Result<Box<dyn SequenceReader>, String> {
        let mut input = errors.capture(FastaFile::open_source_with_buffer(
            &[path.to_path_buf()],
            FastaFileConfig {
                sequence_type: SequenceType::AminoAcid,
                ..Default::default()
            },
            65_536,
        ))?;
        input.set_io_errors(errors.clone());
        Ok(Box::new(FastaSequenceReader {
            input,
            errors: errors.clone(),
        }))
    }
    fn read_records(&self, path: &Path) -> Result<Vec<SequenceRecord>, String> {
        self.read_records_bytes(path)?
            .into_iter()
            .map(|record| {
                Ok(SequenceRecord {
                    id: String::from_utf8(record.id).map_err(|error| error.to_string())?,
                    sequence: record.sequence,
                })
            })
            .collect()
    }
    fn read_records_bytes(&self, path: &Path) -> Result<Vec<ByteSequenceRecord>, String> {
        self.read_records_bytes_with_errors(path, &IoErrorOwner::default())
    }
    fn read_records_bytes_with_errors(
        &self,
        path: &Path,
        errors: &IoErrorOwner,
    ) -> Result<Vec<ByteSequenceRecord>, String> {
        let mut input = self.open_records(path, errors)?;
        let mut records = Vec::new();
        while let Some(record) = input.read_record()? {
            records.push(record);
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
    get_reps_path(job, volumes, config, source)?
        .into_os_string()
        .into_string()
        .map_err(|_| "Representative path is not UTF-8; use get_reps_path".to_owned())
}

pub fn get_reps_path<J, S>(
    job: &J,
    volumes: &[Volume],
    config: GetRepsConfig,
    source: &S,
) -> Result<PathBuf, String>
where
    J: RepresentativeJob,
    S: SequenceSource,
{
    if job.last_round() {
        return Ok(PathBuf::new());
    }
    if config.threads == 0 {
        return Err("Thread count must be positive.".to_string());
    }

    let base_dir = job.base_dir().join("reps");
    std::fs::create_dir_all(source_c_path(&base_dir).as_ref())
        .map_err(|error| error.to_string())?;
    let reps_list = FileStack::try_new_owned(base_dir.join("reps.tsv"), job.io_errors())?;
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
    Ok(reps_list.file_name().to_path_buf())
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
        let mut representatives =
            SourceTextInput::open(&id_file, &job.io_errors(), b"Error opening file ")?;
        // Source opens the sequence before output but reads no sequence or
        // representative until both constructors have completed.
        let mut input = source.open_records(&volume.path, &job.io_errors())?;

        let out_file = base_dir.join(format!("{volume_index}.faa"));
        let mut output =
            TextOutput::open(&out_file, &job.io_errors(), b"Error opening file ", true)?;
        let mut count = 0u64;
        let mut representative = None;
        representatives.number_into(&mut representative);
        let mut file_oid = volume.oid_begin;
        while !stop.load(Ordering::Relaxed) {
            let Some(record) = input.read_record()? else {
                break;
            };
            if job.round() > 0 {
                file_oid = atoll_oid_bytes(&record.id);
            }
            if representative == Some(file_oid) {
                output.write(format!(">{file_oid}\n").as_bytes());
                output.write(
                    Sequence::new(&record.sequence)
                        .to_string(AMINO_ACID_ALPHABET)
                        .as_bytes(),
                );
                output.write(b"\n");
                count += 1;
                representatives.number_into(&mut representative);
            }
            file_oid = file_oid.wrapping_add(1);
        }

        if let Some(rep) = representative.filter(|_| representatives.good()) {
            // Defined num_get conversion values and sticky old values are
            // retained even when a later extraction sets failbit. Initial
            // sentry/EOF has no initialized Source value and stays None.
            return Err(path_error(
                &job.io_errors(),
                format!("Failed to find oid {rep} in file ").as_bytes(),
                &volume.path,
            ));
        }
        input.close()?;
        drop(output);
        let mut line = job.io_errors().capture(
            crate::util::io::source_os_bytes(out_file.as_os_str())
                .map(|v| v.into_owned())
                .map_err(std::io::Error::other),
        )?;
        line.extend_from_slice(
            format!("\t{count}\t{}\t{}\n", volume.oid_begin, volume.oid_end).as_bytes(),
        );
        reps_list.push_bytes(&line)?;
        volumes_processed.fetch_add(1, Ordering::Relaxed);
        cluster_count.fetch_add(count, Ordering::Relaxed);
    }
}

/// `std::atoll` prefix parsing used for representative FASTA IDs after round
/// zero. Invalid/non-numeric IDs map to zero like C `atoll`.
#[cfg(test)]
fn atoll_oid(value: &str) -> OId {
    atoll_oid_bytes(value.as_bytes())
}
fn atoll_oid_bytes(value: &[u8]) -> OId {
    let mut terminated = value.split(|&b| b == 0).next().unwrap_or_default().to_vec();
    terminated.push(0);
    unsafe { libc::strtoll(terminated.as_ptr().cast(), std::ptr::null_mut(), 10) as OId }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::source_text::expected_text;
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

    #[test]
    fn sequence_read_failure_occurs_after_output_construction() {
        let dir = TestDir::new("read-after-output");
        let path = dir.0.join("input.faa");
        std::fs::write(&path, b">bad\nA?\n").unwrap();
        std::fs::write(dir.0.join("rep_ids0"), b"0").unwrap();
        let output_path = dir.0.join("reps").join("0.faa");
        std::fs::create_dir_all(&output_path).unwrap();
        let error = get_reps_path(
            &job(&dir.0, 0, false),
            &[Volume::new(path, 0, 1, 1)],
            GetRepsConfig::default(),
            &FastaSequenceSource,
        )
        .unwrap_err();
        assert_eq!(
            error,
            format!("Error opening file {}", output_path.display())
        );
    }

    #[test]
    fn failed_rep_extraction_retains_defined_values_for_later_numeric_titles() {
        for (name, rep_ids, fasta, expected) in [
            (
                "eof",
                b"7".as_slice(),
                b">7\nARN\n>7\nDCQ\n".as_slice(),
                b">7\nARN\n>7\nDCQ\n".as_slice(),
            ),
            (
                "num-get",
                b"8x".as_slice(),
                b">8\nARN\n>0\nDCQ\n>0\nEGH\n".as_slice(),
                b">8\nARN\n>0\nDCQ\n>0\nEGH\n".as_slice(),
            ),
        ] {
            let dir = TestDir::new(name);
            let path = dir.0.join("input.faa");
            std::fs::write(&path, fasta).unwrap();
            std::fs::write(dir.0.join("rep_ids0"), rep_ids).unwrap();
            get_reps_path(
                &job(&dir.0, 1, false),
                &[Volume::new(path, 0, 3, 3)],
                GetRepsConfig::default(),
                &FastaSequenceSource,
            )
            .unwrap();
            assert_eq!(
                std::fs::read(dir.0.join("reps/0.faa")).unwrap(),
                expected_text(expected)
            );
        }
    }

    #[test]
    fn raw_title_and_os_path_representatives_use_the_byte_source() {
        let dir = TestDir::new("raw-representatives");
        #[cfg(unix)]
        let (path, _) =
            crate::util::test_os_path_or_unicode(&dir.0, b"input\xff.faa", "input-é.faa");
        #[cfg(not(unix))]
        let path = dir.0.join("input-é.faa");
        std::fs::write(&path, b">10\xff title\nARN\n").unwrap();
        let records = FastaSequenceSource.read_records_bytes_io(&path).unwrap();
        assert_eq!(records[0].id, b"10\xff title");
        assert!(FastaSequenceSource.read_records(&path).is_err());
        std::fs::write(dir.0.join("rep_ids0"), b"10").unwrap();
        let manifest = get_reps_path(
            &job(&dir.0, 1, false),
            &[Volume::new(path, 10, 11, 1)],
            GetRepsConfig { threads: 1 },
            &FastaSequenceSource,
        )
        .unwrap();
        assert_eq!(
            std::fs::read(dir.0.join("reps/0.faa")).unwrap(),
            expected_text(b">10\nARN\n")
        );
        let rows = super::super::multinode::VolumedFile::from_manifest(&manifest).unwrap();
        assert_eq!(rows.0[0].path, dir.0.join("reps/0.faa"));
        assert_eq!((rows.0[0].oid_begin, rows.0[0].oid_end), (10, 11));
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
        let output = dir.0.join("reps").join("0.faa");
        assert_eq!(
            std::fs::read(&output).unwrap(),
            expected_text(b">11\nDCQ\n>12\nEGH\n")
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
            std::fs::read(dir.0.join("reps/0.faa")).unwrap(),
            expected_text(b">7\nDCQ\n>999\nEGH\n")
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
