//! Native sequence-file adapter shared by the cluster command workflows.
use crate::basic::value::{Letter, OId, SequenceType};
use crate::cluster::helpers::ClusterHelperDatabase;
use crate::cluster::output::ClusterRealignDatabase;
use crate::data::blastdb::blastdb::BlastDB;
use crate::data::block::Block;
use crate::data::dmnd::DatabaseFile;
use crate::data::fasta::{
    fasta_file::{FastaFile, FastaFileConfig},
    ByteFastaRecord as FastaRecord,
};
use crate::data::sequence_file::{
    detect_type_source_with_buffer, AccessionIndex, SequenceFileConfig, SequenceFileFlags,
    SequenceFileType,
};
use std::path::PathBuf;
use std::sync::Arc;

enum Source {
    Dmnd(Box<DatabaseFile>, OId),
    Fasta(Box<FastaFile>),
    Blast(Box<BlastDB>),
    Memory {
        records: Vec<FastaRecord>,
        cursor: OId,
    },
}

impl Source {
    fn read(&mut self) -> Result<Option<(Vec<Letter>, Vec<u8>)>, String> {
        match self {
            Self::Dmnd(db, cursor) => {
                if *cursor == db.sequence_count() {
                    return Ok(None);
                }
                let record = db
                    .record_bytes_io(*cursor)
                    .map_err(|error| error.to_string())?;
                *cursor += 1;
                Ok(Some(record))
            }
            Self::Fasta(db) => Ok(db
                .read_seq_bytes_io()
                .map_err(|error| error.to_string())?
                .map(|entry| (entry.sequence, entry.id))),
            Self::Blast(db) => {
                if db.eof() {
                    Ok(None)
                } else {
                    db.read_seq_bytes().map(Some)
                }
            }
            Self::Memory { records, cursor } => {
                let record = records.get(*cursor as usize).cloned();
                if record.is_some() {
                    *cursor += 1;
                }
                Ok(record.map(|record| (record.sequence, record.id)))
            }
        }
    }

    fn read_owned(
        &mut self,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<Option<(Vec<Letter>, Vec<u8>)>, String> {
        if let Self::Dmnd(db, cursor) = self {
            if *cursor == db.sequence_count() {
                return Ok(None);
            }
            let record = errors.capture(db.record_bytes_io(*cursor))?;
            *cursor += 1;
            Ok(Some(record))
        } else if let Self::Fasta(db) = self {
            Ok(errors
                .capture(db.read_seq_bytes_io())?
                .map(|entry| (entry.sequence, entry.id)))
        } else if let Self::Blast(db) = self {
            if db.eof() {
                Ok(None)
            } else {
                errors.capture(db.source_io(|db| db.read_seq_bytes().map(Some)))
            }
        } else {
            self.read()
        }
    }

    fn seek(&mut self, oid: OId) -> Result<(), String> {
        match self {
            Self::Dmnd(db, cursor) => {
                db.set_seqinfo_ptr(oid)?;
                *cursor = oid;
                Ok(())
            }
            Self::Fasta(db) => db.set_seqinfo_ptr(oid),
            Self::Blast(db) => db.set_seqinfo_ptr(oid),
            Self::Memory { records, cursor } => {
                if oid > records.len() as OId {
                    return Err("Sequence index out of bounds.".to_owned());
                }
                *cursor = oid;
                Ok(())
            }
        }
    }

    fn tell(&self) -> OId {
        match self {
            Self::Dmnd(_, cursor) => *cursor,
            Self::Fasta(db) => db.tell_seq(),
            Self::Blast(db) => db.tell_seq(),
            Self::Memory { cursor, .. } => *cursor,
        }
    }
}

pub(crate) struct NativeDatabase {
    source: Source,
    index: AccessionIndex,
    pub(crate) count: OId,
    pub(crate) letters: u64,
    titles: bool,
    accession_spool: Option<std::sync::Arc<crate::data::accession_spool::AccessionSpool>>,
    accession_table: Option<crate::util::tsv::Table>,
    errors: crate::util::io::IoErrorOwner,
}

impl NativeDatabase {
    pub(crate) fn open(path: &str, flags: SequenceFileFlags) -> Result<Self, String> {
        Self::open_path(std::path::Path::new(path), flags)
    }
    pub(crate) fn open_with_buffer(
        path: &str,
        flags: SequenceFileFlags,
        file_buffer_size: Option<usize>,
    ) -> Result<Self, String> {
        Self::open_with_buffer_path(std::path::Path::new(path), flags, file_buffer_size)
    }
    pub(crate) fn open_source_with_buffer(
        path: &str,
        flags: SequenceFileFlags,
        file_buffer_size: Option<usize>,
        tmpdir: &std::path::Path,
        no_unlink: bool,
    ) -> std::io::Result<Self> {
        Self::open_source_with_buffer_path(
            std::path::Path::new(path),
            flags,
            file_buffer_size,
            tmpdir,
            no_unlink,
        )
    }
    pub(crate) fn open_source_with_buffer_owned(
        path: &str,
        flags: SequenceFileFlags,
        file_buffer_size: Option<usize>,
        tmpdir: &std::path::Path,
        no_unlink: bool,
        errors: &crate::util::io::IoErrorOwner,
    ) -> std::io::Result<Self> {
        Self::open_source_with_buffer_owned_path(
            std::path::Path::new(path),
            flags,
            file_buffer_size,
            tmpdir,
            no_unlink,
            errors,
        )
    }

    pub(crate) fn from_records(records: Vec<FastaRecord>) -> Result<Self, String> {
        let mut index = AccessionIndex::default();
        for (oid, record) in records.iter().enumerate() {
            index
                .add_seqid_mapping_bytes(&record.id, oid as OId)
                .map_err(|error| error.to_string())?;
        }
        Ok(Self {
            count: records.len() as OId,
            letters: records
                .iter()
                .map(|record| record.sequence.len() as u64)
                .sum(),
            source: Source::Memory { records, cursor: 0 },
            index,
            titles: true,
            accession_spool: None,
            accession_table: None,
            errors: Default::default(),
        })
    }

    pub(crate) fn open_path(
        path: &std::path::Path,
        flags: SequenceFileFlags,
    ) -> Result<Self, String> {
        Self::open_with_buffer_path(path, flags, None)
    }

    pub(crate) fn open_with_buffer_path(
        path: &std::path::Path,
        flags: SequenceFileFlags,
        file_buffer_size: Option<usize>,
    ) -> Result<Self, String> {
        Self::open_source_with_buffer_path(
            path,
            flags,
            file_buffer_size,
            std::path::Path::new(""),
            false,
        )
        .map_err(|error| error.to_string())
    }

    pub(crate) fn open_source_with_buffer_path(
        path: &std::path::Path,
        flags: SequenceFileFlags,
        file_buffer_size: Option<usize>,
        tmpdir: &std::path::Path,
        no_unlink: bool,
    ) -> std::io::Result<Self> {
        Self::open_source_with_buffer_owned_path(
            path,
            flags,
            file_buffer_size,
            tmpdir,
            no_unlink,
            &Default::default(),
        )
    }

    pub(crate) fn open_source_with_buffer_owned_path(
        path: &std::path::Path,
        flags: SequenceFileFlags,
        file_buffer_size: Option<usize>,
        tmpdir: &std::path::Path,
        no_unlink: bool,
        errors: &crate::util::io::IoErrorOwner,
    ) -> std::io::Result<Self> {
        if let Some(size) = file_buffer_size {
            Vec::<u8>::new().try_reserve_exact(size).map_err(|_| std::io::Error::other("Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage."))?;
        }
        let (kind, path) = detect_type_source_with_buffer(
            &[PathBuf::from(path)],
            flags,
            SequenceFileConfig::default(),
            file_buffer_size.unwrap_or(67_108_864),
        )?;
        // SequenceFile's base constructor creates this real TSV owner before
        // the concrete database constructor reads its sequence payload.
        let accession_spool = flags
            .contains(SequenceFileFlags::OID_TO_ACC_MAPPING)
            .then(|| {
                crate::data::accession_spool::AccessionSpool::new_owned(
                    tmpdir,
                    no_unlink,
                    file_buffer_size.unwrap_or(67_108_864),
                    errors.clone().with_source_destructors(),
                )
            })
            .transpose()?
            .map(std::sync::Arc::new);
        let result = (|| {
            let mut index = AccessionIndex::default();
            let mut fasta_totals = None;
            let mut source = match kind {
                SequenceFileType::Dmnd => Source::Dmnd(
                    // Source SequenceFile::auto_create constructs DatabaseFile
                    // on its retained InputFile, even for supplied cluster owners.
                    // Keep constructor flags unchanged; records still use SeqInfo.
                    Box::new(DatabaseFile::open_source_with_options_io(
                        &path,
                        crate::data::dmnd::DatabaseOpenOptions {
                            file_buffer_size,
                            ..Default::default()
                        },
                        false,
                    )?),
                    0,
                ),
                SequenceFileType::Fasta => {
                    let mut file = FastaFile::open_source_with_mapping(
                        &[path],
                        FastaFileConfig {
                            flags,
                            ..Default::default()
                        },
                        file_buffer_size.unwrap_or(67_108_864),
                        &mut |title, oid| {
                            index.add_seqid_mapping_bytes(title, oid)?;
                            if let Some(spool) = &accession_spool {
                                spool.push_bytes(&crate::util::sequence::seqid_bytes(title))?;
                            }
                            Ok(())
                        },
                    )?;
                    if flags.contains(SequenceFileFlags::NEED_LETTER_COUNT)
                        && flags.contains(
                            SequenceFileFlags::ACC_TO_OID_MAPPING
                                | SequenceFileFlags::OID_TO_ACC_MAPPING,
                        )
                    {
                        // Source's constructor already made this metadata and
                        // mapping pass, including its rewind. Retain that owner;
                        // a second scan changes constructor failure side effects.
                        fasta_totals = Some((file.sequence_count(), file.letters()));
                    }
                    // Later reads cross the cluster traits' String boundary.
                    // Keep their actual Source I/O diagnostic in this invocation.
                    file.set_io_errors(errors.clone());
                    Source::Fasta(Box::new(file))
                }
                SequenceFileType::Blast => Source::Blast(Box::new(
                    BlastDB::new_path_with_io_owner(&path, flags, errors)?,
                )),
                SequenceFileType::Block => {
                    unreachable!("file detection does not create an in-memory block")
                }
            };
            let (mut count, mut letters) = fasta_totals.unwrap_or_default();
            if fasta_totals.is_none() {
                while let Some((sequence, title)) = errors.finish(source.read_owned(errors))? {
                    // BLAST accession lookup remains delegated to its own decoder,
                    // but output still needs an ordinal-to-title entry per record.
                    index.add_seqid_mapping_bytes(&title, count)?;
                    if let Some(spool) = &accession_spool {
                        spool.push_bytes(&crate::util::sequence::seqid_bytes(&title))?;
                    }
                    count += 1;
                    letters += sequence.len() as u64;
                }
                errors.finish(source.seek(0))?;
            }
            Ok(Self {
                source,
                index,
                count,
                letters,
                titles: true,
                accession_spool: accession_spool.clone(),
                accession_table: None,
                errors: errors.clone(),
            })
        })();
        errors.capture_io(result)
    }

    pub(crate) fn prepare_output_accessions(&mut self, threads: usize) -> std::io::Result<()> {
        if self.accession_table.is_none() {
            if let Some(spool) = &self.accession_spool {
                self.accession_table = Some(spool.read_table(threads)?);
            }
        }
        Ok(())
    }
    /// Source `sub_db`: preserve the selected OIDs' order while streaming the
    /// original database. The translated callers supply increasing OIDs.
    pub(crate) fn selected_records(&mut self, ids: &[OId]) -> Result<Vec<FastaRecord>, String> {
        self.source.seek(0)?;
        let mut selected = Vec::with_capacity(ids.len());
        let mut next = 0;
        while next < ids.len() {
            let oid = self.source.tell();
            let Some((sequence, id)) = self.source.read_owned(&self.errors)? else {
                return Err("Subdatabase OID out of range.".to_owned());
            };
            if oid == ids[next] {
                selected.push(FastaRecord {
                    quality: None,
                    id,
                    sequence,
                });
                next += 1;
            }
        }
        Ok(selected)
    }
}

impl ClusterHelperDatabase for NativeDatabase {
    fn sequence_count(&self) -> usize {
        self.count as usize
    }
    fn accession_to_oid(&self, accession: &str) -> Result<Vec<OId>, String> {
        if let Source::Blast(db) = &self.source {
            return db.accession_to_oid(accession);
        }
        self.index.accession_to_oid(accession)
    }
    fn accession_to_oid_bytes(&self, accession: &[u8]) -> Result<Vec<OId>, String> {
        if let Source::Blast(db) = &self.source {
            return self.errors.capture(db.accession_to_oid_bytes(accession));
        }
        self.errors
            .capture(self.index.accession_to_oid_bytes(accession))
    }
    fn accession_bytes(&self, oid: OId) -> Result<&[u8], String> {
        if let Some(table) = &self.accession_table {
            if oid >= table.size() as u64 {
                return Err("Accession OID outside table.".into());
            }
            return Ok(table.record(oid as i64).string_bytes(0));
        }
        self.index.seqid_bytes(oid)
    }
    fn accession(&self, oid: OId) -> Result<&str, String> {
        if let Some(table) = &self.accession_table {
            if oid >= table.size() as u64 {
                return Err("Accession OID outside table.".into());
            }
            return std::str::from_utf8(table.record(oid as i64).string_bytes(0))
                .map_err(|error| error.to_string());
        }
        self.index.seqid(oid)
    }
}

impl ClusterRealignDatabase for NativeDatabase {
    fn sequence_count(&self) -> OId {
        self.count
    }
    fn letters(&self) -> u64 {
        self.letters
    }
    fn titles_lazy(&self) -> bool {
        false
    }
    fn set_load_flags(&mut self, _sequences: bool, titles: bool) {
        self.titles = titles;
    }
    fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String> {
        self.source.seek(oid)
    }
    fn tell_seq(&self) -> OId {
        self.source.tell()
    }
    fn seqid_bytes(&mut self, oid: OId) -> Result<Vec<u8>, String> {
        Ok(self.index.seqid_bytes(oid)?.to_vec())
    }
    fn seqid(&mut self, oid: OId) -> Result<String, String> {
        Ok(self.index.seqid(oid)?.to_owned())
    }
    fn load_seqs(&mut self, block_size: u64) -> Result<Arc<Block>, String> {
        let mut block = Block::new();
        let mut letters = 0;
        loop {
            let oid = self.source.tell();
            let Some((sequence, mut title)) = self.source.read_owned(&self.errors)? else {
                break;
            };
            crate::util::sequence::fix_title_bytes(&mut title);
            letters += block.push_back_bytes(
                &sequence,
                self.titles.then_some(title.as_slice()),
                None,
                oid,
                SequenceType::AminoAcid,
                1,
                false,
            )? as u64;
            if letters >= block_size {
                break;
            }
        }
        Ok(Arc::new(block))
    }
}

#[cfg(test)]
mod byte_tests {
    use super::*;
    #[test]
    fn raw_cluster_mapping_uses_byte_labels_and_actual_text_reader() {
        let database = NativeDatabase::from_records(vec![FastaRecord {
            id: b"q\xff".to_vec(),
            quality: None,
            sequence: vec![0, 1, 2, 3],
        }])
        .unwrap();
        let path = std::env::temp_dir().join(format!(
            "diamond-byte-cluster-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, b"q\xff\tq\xff\n").unwrap();
        let mapping = crate::cluster::helpers::read_mapping_with_source_buffer_owned::<OId, _>(
            &path,
            &database,
            false,
            64,
            &database.errors,
        )
        .unwrap();
        assert_eq!(mapping, [0]);
        assert_eq!(database.accession_bytes(0).unwrap(), b"q\xff");
        std::fs::remove_file(path).unwrap();
    }
}
