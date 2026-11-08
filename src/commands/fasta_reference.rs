//! Seekable native cache for a Source FASTA/FASTQ reference owner.
//! Only record offsets are retained in memory; sequences remain unmasked and
//! qualities survive random-access reloads. This is not a DIAMOND database.
use crate::util::io::source_c_path;
use std::fs::{File, OpenOptions};
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::Arc;

use super::{BlastpConfig, FastaRecord, RecordStream, TemporaryDatabase};
use crate::basic::value::SequenceType;

struct Index {
    owner: TemporaryDatabase,
    original_path: PathBuf,
    offsets: Vec<u64>,
    accessions: Option<Arc<crate::data::accession_spool::AccessionSpool>>,
    letters: u64,
}

pub(super) struct FastaReference {
    // Close the final handle before dropping the last cache owner on Windows.
    file: File,
    index: Arc<Index>,
    pub cursor: u64,
}

impl FastaReference {
    pub fn open(
        config: &BlastpConfig,
        path: PathBuf,
        sequence_type: SequenceType,
        errors: crate::util::io::IoErrorOwner,
    ) -> io::Result<Self> {
        let directory = crate::util::io::temp_file::source_temp_directory_path(
            &config.tmpdir,
            config.output.as_deref().unwrap_or(std::path::Path::new("")),
        );
        let accessions = if config.tuning.unaligned_targets.as_os_str().is_empty() {
            None
        } else {
            Some(Arc::new(
                crate::data::accession_spool::AccessionSpool::new_owned(
                    &directory,
                    config.tuning.no_unlink,
                    config.file_buffer_size,
                    errors.clone(),
                )?,
            ))
        };
        let result = (|| {
            let reserved = crate::util::io::TempFileData::init_in(false, directory, false)
                .map_err(io::Error::other)?;
            let owner = TemporaryDatabase(PathBuf::from(&reserved.name), config.tuning.no_unlink);
            let mut file = OpenOptions::new()
                .read(true)
                .write(true)
                .create(cfg!(windows))
                .truncate(cfg!(windows))
                .open(source_c_path(&owner.0).as_ref())
                .map_err(|error| crate::util::io::FileOpenError::into_io(&owner.0, error))?;
            // TempFileData owns an open path too; retain the cache independently.
            let mut retained = reserved.name.as_os_str().to_os_string();
            retained.push(".fasta-cache");
            let retained = PathBuf::from(retained);
            std::fs::rename(
                source_c_path(&owner.0).as_ref(),
                source_c_path(&retained).as_ref(),
            )?;
            let owner = TemporaryDatabase(retained, config.tuning.no_unlink);
            drop(reserved);
            let mut input = RecordStream::open_ready(
                path.clone(),
                sequence_type,
                true,
                config.file_buffer_size,
            )?;
            let mut offsets = Vec::new();
            let mut letters = 0u64;
            while let Some(record) = input.next()? {
                if let Some(spool) = &accessions {
                    spool.push_bytes(&crate::util::sequence::seqid_bytes(&record.id))?;
                }
                offsets.push(file.stream_position()?);
                let quality = record.quality.as_deref();
                file.write_all(&(record.id.len() as u64).to_le_bytes())?;
                file.write_all(&(record.sequence.len() as u64).to_le_bytes())?;
                file.write_all(
                    &quality
                        .map_or(u64::MAX, |quality| quality.len() as u64)
                        .to_le_bytes(),
                )?;
                file.write_all(record.id.as_slice())?;
                // Letter is i8; the cache stores its complete byte representation,
                // including masking bits. Bound conversion memory independently of
                // the record length and avoid reinterpreting any other type.
                for letters in record.sequence.chunks(8192) {
                    let mut bytes = [0u8; 8192];
                    for (byte, &letter) in bytes.iter_mut().zip(letters) {
                        *byte = letter as u8;
                    }
                    file.write_all(&bytes[..letters.len()])?;
                }
                if let Some(quality) = quality {
                    file.write_all(quality)?;
                }
                letters = letters
                    .checked_add(record.sequence.len() as u64)
                    .ok_or_else(|| io::Error::other("Reference letter count overflow"))?;
            }
            // Source init_read counts the real input and then seeks its owner to
            // OID zero. Preserve the failure of nonseekable stdin/FIFO references;
            // the internal cache must not silently widen that input domain.
            input.rewind_after_count()?;
            file.flush()?;
            Ok(Self {
                index: Arc::new(Index {
                    owner,
                    original_path: path,
                    offsets,
                    accessions: accessions.clone(),
                    letters,
                }),
                file,
                cursor: 0,
            })
        })();
        errors.capture_io(result)
    }

    pub fn duplicate(&self) -> io::Result<Self> {
        Ok(Self {
            file: File::open(source_c_path(&self.index.owner.0).as_ref())?,
            index: self.index.clone(),
            cursor: 0,
        })
    }
    pub fn accession_spool(&self) -> Option<Arc<crate::data::accession_spool::AccessionSpool>> {
        self.index.accessions.clone()
    }
    pub fn totals(&self) -> (u64, u64) {
        (self.index.offsets.len() as u64, self.index.letters)
    }
    pub fn file_name(&self) -> &std::path::Path {
        &self.index.original_path
    }
    /// Title-only access through the existing physical-record offset index.
    /// Final accession output must not decode or retain reference sequences.
    pub fn title(&mut self, oid: u64) -> io::Result<Vec<u8>> {
        let offset = *self.index.offsets.get(oid as usize).ok_or_else(|| {
            io::Error::other(crate::util::io::IoError::path_message(
                b"FASTA reference OID outside database ",
                self.file_name(),
            ))
        })?;
        self.file.seek(SeekFrom::Start(offset))?;
        let mut header = [0u8; 24];
        self.file.read_exact(&mut header)?;
        let count = usize::try_from(u64::from_le_bytes(header[0..8].try_into().unwrap()))
            .map_err(|_| io::Error::other("FASTA record exceeds address space"))?;
        let mut title = Vec::new();
        title.try_reserve_exact(count).map_err(io::Error::other)?;
        title.resize(count, 0);
        self.file.read_exact(&mut title)?;
        Ok(title)
    }

    pub fn record(&mut self, oid: u64) -> io::Result<FastaRecord> {
        let offset = *self.index.offsets.get(oid as usize).ok_or_else(|| {
            io::Error::other(crate::util::io::IoError::path_message(
                b"FASTA reference OID outside database ",
                self.file_name(),
            ))
        })?;
        self.file.seek(SeekFrom::Start(offset))?;
        let mut header = [0u8; 24];
        self.file.read_exact(&mut header)?;
        let title_len = u64::from_le_bytes(header[0..8].try_into().unwrap());
        let seq_len = u64::from_le_bytes(header[8..16].try_into().unwrap());
        let qual_len = u64::from_le_bytes(header[16..24].try_into().unwrap());
        fn read_bytes(file: &mut File, count: u64) -> io::Result<Vec<u8>> {
            let count = usize::try_from(count)
                .map_err(|_| io::Error::other("FASTA record exceeds address space"))?;
            let mut bytes = Vec::new();
            bytes.try_reserve_exact(count).map_err(io::Error::other)?;
            bytes.resize(count, 0);
            file.read_exact(&mut bytes)?;
            Ok(bytes)
        }
        let id = read_bytes(&mut self.file, title_len)?;
        let sequence = read_bytes(&mut self.file, seq_len)?
            .into_iter()
            .map(|byte| byte as crate::basic::value::Letter)
            .collect();
        let quality = if qual_len == u64::MAX {
            None
        } else {
            Some(read_bytes(&mut self.file, qual_len)?)
        };
        Ok(FastaRecord {
            id,
            sequence,
            quality,
        })
    }
}
