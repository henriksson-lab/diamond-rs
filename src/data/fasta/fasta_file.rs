//! Stateful translation of `data/fasta/fasta_file.{h,cpp}`.
//! Source records are tokenized from retained streams; only requested index,
//! length and accession metadata survive a counting pass. Random-access and
//! partition helpers explicitly identified below remain Rust library extensions.
#[path = "source_reader.rs"]
mod source_reader;
use crate::basic::value::{
    CharRepresentation, Letter, Loc, OId, SequenceType, AMINO_ACID_ALPHABET, MASK_LETTER,
    NUCLEOTIDE_ALPHABET,
};
use crate::data::fasta::SeqFileFormat;
use crate::data::sequence_file::{Chunk, SequenceFileFlags};
use crate::util::io::{Compressor, InputFile, IoErrorOwner, OutputFile, TempFile};
pub(crate) use source_reader::source_fastq_fields;
use source_reader::RecordReader;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

#[derive(Debug, Clone)]
pub struct FastaFileConfig {
    pub flags: SequenceFileFlags,
    pub sequence_type: SequenceType,
    pub index_file: Option<PathBuf>,
}
impl Default for FastaFileConfig {
    fn default() -> Self {
        Self {
            flags: SequenceFileFlags::NONE,
            sequence_type: SequenceType::AminoAcid,
            index_file: None,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaEntry {
    pub id: String,
    pub sequence: Vec<Letter>,
    pub quality: Option<Vec<u8>>,
    pub offset: u64,
}
/// Byte-preserving Source record. The existing String record API remains a
/// checked UTF-8 convenience and cannot represent arbitrary Source title bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaByteEntry {
    pub id: Vec<u8>,
    pub sequence: Vec<Letter>,
    pub quality: Option<Vec<u8>>,
    pub offset: u64,
}

enum Writer {
    File(OutputFile),
    Temp(TempFile),
}
impl Writer {
    fn output(&mut self) -> &mut OutputFile {
        match self {
            Self::File(file) => file,
            Self::Temp(file) => file,
        }
    }
}

pub struct FastaFile {
    paths: Vec<PathBuf>,
    files: Mutex<Vec<RecordReader>>,
    file_ptr: usize,
    oid: OId,
    format: SeqFileFormat,
    seqs: u64,
    letters: u64,
    seq_lengths: Vec<Loc>,
    index: Vec<u64>,
    partition: Vec<Chunk>,
    writer: Option<Writer>,
    sequence_type: SequenceType,
    flags: SequenceFileFlags,
    buffer_size: usize,
    errors: IoErrorOwner,
}
impl FastaFile {
    pub fn open(paths: &[PathBuf], config: FastaFileConfig) -> Result<Self, String> {
        Self::open_source_with_buffer(paths, config, crate::util::io::DEFAULT_FILE_BUFFER_SIZE)
            .map_err(|error| error.to_string())
    }
    pub fn open_source_with_buffer(
        paths: &[PathBuf],
        config: FastaFileConfig,
        buffer_size: usize,
    ) -> io::Result<Self> {
        Self::open_source_with_mapping(paths, config, buffer_size, &mut |_, _| Ok(()))
    }

    /// Bind Source SequenceFile's already-created mapping owner to the actual
    /// init_read pass. Index-based initialization never invokes this callback.
    pub fn open_source_with_mapping(
        paths: &[PathBuf],
        config: FastaFileConfig,
        buffer_size: usize,
        on_mapping: &mut dyn FnMut(&[u8], OId) -> io::Result<()>,
    ) -> io::Result<Self> {
        if paths.is_empty() {
            return Err(io::Error::other("Missing FASTA file name."));
        }
        if paths.len() > 2 {
            return Err(io::Error::other("std::exception"));
        }
        let taxon_flags = SequenceFileFlags::TAXON_MAPPING
            | SequenceFileFlags::TAXON_NODES
            | SequenceFileFlags::TAXON_RANKS
            | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES;
        if config.flags.contains(taxon_flags) {
            return Err(io::Error::other(
                "Fasta database format does not support taxonomic features.",
            ));
        }
        let first = InputFile::with_source_profile_path(&paths[0], 0, buffer_size)
            .map_err(io::Error::other)?;
        let first = RecordReader::open(first, None)?;
        let format = first.format;
        let mut files = vec![first];
        if paths.len() == 2 {
            let second = InputFile::with_source_profile_path(&paths[1], 0, buffer_size)
                .map_err(io::Error::other)?;
            // Source uses the FIRST file's tokenizer without guessing again.
            files.push(RecordReader::open(second, Some(format))?);
        }
        let mut result = Self {
            paths: paths.to_vec(),
            files: Mutex::new(files),
            file_ptr: 0,
            oid: 0,
            format,
            seqs: 0,
            letters: 0,
            seq_lengths: Vec::new(),
            index: Vec::new(),
            partition: Vec::new(),
            writer: None,
            sequence_type: config.sequence_type,
            flags: config.flags,
            buffer_size,
            errors: IoErrorOwner::default().with_source_destructors(),
        };
        // Source leaves seqs_/letters_ indeterminate without this flag. Rust
        // exposes deterministic unavailable counts (zero), never guessed totals.
        let initialize = (|| -> io::Result<()> {
            if config.flags.contains(SequenceFileFlags::NEED_LETTER_COUNT) {
                if let Some(index) = config.index_file {
                    for (offset, length) in read_index(&index)? {
                        result.index.push(offset);
                        result.seq_lengths.push(length);
                        result.seqs += 1;
                        result.letters = result.letters.wrapping_add(length as u64);
                    }
                } else {
                    while let Some(record) = result.read_seq_bytes_io()? {
                        if config.flags.contains(
                            SequenceFileFlags::ACC_TO_OID_MAPPING
                                | SequenceFileFlags::OID_TO_ACC_MAPPING,
                        ) {
                            on_mapping(&record.id, result.seqs as OId)?;
                        }
                        if config.flags.contains(SequenceFileFlags::NEED_LENGTH_LOOKUP) {
                            result.seq_lengths.push(record.sequence.len() as Loc);
                        }
                        result.seqs += 1;
                        result.letters += record.sequence.len() as u64;
                    }
                }
                result.set_seqinfo_ptr_io(0)?;
            }
            Ok(())
        })();
        if let Err(error) = initialize {
            // Source's failing low-level operation has already called perror
            // before member-destructor cleanup can itself fail/terminate. Keep
            // that ordering, and carry the marker through the outer catch.
            let error = if crate::util::io::emit_source_perror(&error).map_err(io::Error::other)? {
                crate::util::io::PerrorEmittedError::wrap(error)
            } else {
                error
            };
            return Err(error);
        }
        Ok(result)
    }
    pub(crate) fn set_io_errors(&mut self, owner: IoErrorOwner) {
        // Both FastaFile and its retained TSV File members close streams
        // from noexcept destructors in Source, including constructor unwind.
        self.errors = owner.with_source_destructors();
    }
    pub fn io_errors(&self) -> IoErrorOwner {
        self.errors.clone()
    }

    pub fn create(
        path: impl AsRef<Path>,
        overwrite: bool,
        sequence_type: SequenceType,
    ) -> Result<Self, String> {
        Self::create_source_with_buffer(
            path.as_ref(),
            overwrite,
            FastaFileConfig {
                sequence_type,
                ..Default::default()
            },
            crate::util::io::DEFAULT_FILE_BUFFER_SIZE,
        )
        .map_err(|error| error.to_string())
    }
    pub fn create_source_with_buffer(
        path: &Path,
        overwrite: bool,
        config: FastaFileConfig,
        buffer_size: usize,
    ) -> io::Result<Self> {
        let sequence_type = config.sequence_type;
        let mut writer = if path.as_os_str().is_empty() {
            Writer::Temp(TempFile::new(true).map_err(io::Error::other)?)
        } else {
            Writer::File(
                OutputFile::with_source_profile_path(
                    path,
                    Compressor::None,
                    if overwrite { "w+b" } else { "r+b" },
                    buffer_size,
                )
                .map_err(io::Error::other)?,
            )
        };
        let source_name = writer.output().file_path().to_path_buf();
        let input = InputFile::from_output_file_with_buffer_size(writer.output(), 0, buffer_size)
            .map_err(io::Error::other)?;
        let reader = RecordReader::open(input, Some(SeqFileFormat::Fasta))?;
        let mut result = Self {
            paths: vec![source_name],
            files: Mutex::new(vec![reader]),
            file_ptr: 0,
            oid: 0,
            format: SeqFileFormat::Fasta,
            seqs: 0,
            letters: 0,
            seq_lengths: Vec::new(),
            index: Vec::new(),
            partition: Vec::new(),
            writer: Some(writer),
            sequence_type,
            flags: config.flags,
            buffer_size,
            errors: IoErrorOwner::default().with_source_destructors(),
        };
        if !overwrite {
            while let Some(record) = result.read_seq_bytes_io()? {
                result.seqs += 1;
                result.letters += record.sequence.len() as u64;
            }
        }
        Ok(result)
    }
    pub fn file_count(&self) -> i64 {
        self.paths.len() as i64
    }
    pub fn files_synced_io(&self) -> io::Result<bool> {
        if self.file_ptr != 0 {
            return Ok(false);
        }
        let mut files = self
            .files
            .lock()
            .map_err(|_| io::Error::other("Poisoned FASTA input owner."))?;
        if files[0].eof && files.len() > 1 && files[1].read()?.is_some() {
            return Ok(false);
        }
        Ok(true)
    }
    pub fn files_synced(&self) -> bool {
        self.errors.capture(self.files_synced_io()).unwrap_or(false)
    }
    pub fn format(&self) -> SeqFileFormat {
        self.format
    }
    pub fn sequence_count(&self) -> u64 {
        self.seqs
    }
    pub fn letters(&self) -> u64 {
        self.letters
    }
    pub fn tell_seq(&self) -> OId {
        self.oid
    }
    pub fn eof(&self) -> bool {
        self.files
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .last()
            .is_some_and(|file| file.eof)
    }
    pub fn file_name(&self) -> Option<&Path> {
        self.paths.first().map(PathBuf::as_path)
    }
    pub fn init_seq_access(&mut self) -> Result<(), String> {
        self.set_seqinfo_ptr(0)
    }
    pub fn set_seqinfo_ptr_io(&mut self, oid: OId) -> io::Result<()> {
        if let Some(writer) = &mut self.writer {
            writer.output().rewind().map_err(io::Error::other)?;
        }
        if self.index.is_empty() {
            for file in self
                .files
                .get_mut()
                .map_err(|_| io::Error::other("Poisoned FASTA input owner."))?
            {
                file.rewind()?;
            }
            // Source deliberately preserves the current paired-file iterator.
            self.oid = 0;
            while self.oid != oid {
                // Source ignores the read result here: even EOF attempts
                // advance oid_, so seeking beyond the record count succeeds.
                self.read_seq_bytes_io()?;
            }
        } else {
            if self.paths.len() != 1 || oid as usize >= self.index.len() {
                return Err(io::Error::other("FastaFile::set_seqinfo_ptr"));
            }
            self.files
                .get_mut()
                .map_err(|_| io::Error::other("Poisoned FASTA input owner."))?[0]
                .seek(self.index[oid as usize])?;
            // Source's indexed branch does not assign oid_.
        }
        Ok(())
    }
    pub fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String> {
        let result = self.set_seqinfo_ptr_io(oid);
        self.errors.capture(result)
    }
    pub fn read_seq_bytes_io(&mut self) -> io::Result<Option<FastaByteEntry>> {
        self.oid = self.oid.wrapping_add(1);
        let files = self
            .files
            .get_mut()
            .map_err(|_| io::Error::other("Poisoned FASTA input owner."))?;
        let raw = files[self.file_ptr].read()?;
        let Some(raw) = raw else {
            return Ok(None);
        };
        let representation = encoding(self.sequence_type);
        let sequence = raw
            .sequence
            .into_iter()
            .filter(|byte| !matches!(byte, b'\n' | b'\r'))
            .map(|byte| representation.convert(byte).map_err(io::Error::other))
            .collect::<io::Result<Vec<_>>>()?;
        self.file_ptr = (self.file_ptr + 1) % files.len();
        Ok(Some(FastaByteEntry {
            id: raw.id,
            sequence,
            quality: raw.quality,
            offset: raw.offset,
        }))
    }
    pub fn read_seq_io(&mut self) -> io::Result<Option<FastaEntry>> {
        self.read_seq_bytes_io()?
            .map(|entry| {
                Ok(FastaEntry {
                    id: String::from_utf8(entry.id).map_err(io::Error::other)?,
                    sequence: entry.sequence,
                    quality: entry.quality,
                    offset: entry.offset,
                })
            })
            .transpose()
    }
    pub fn read_seq(&mut self) -> Result<Option<FastaEntry>, String> {
        let result = self.read_seq_io();
        self.errors.capture(result)
    }
    /// Rust library extension: Source seq_data on FASTA is unsupported. This
    /// scans an independent retained reader instead of caching all payloads.
    pub fn seq_data(&self, oid: OId) -> Result<Vec<Letter>, String> {
        let mut reader = Self::open_source_with_buffer(
            &self.paths,
            FastaFileConfig {
                sequence_type: self.sequence_type,
                ..Default::default()
            },
            self.buffer_size,
        )
        .map_err(|error| error.to_string())?;
        reader.set_seqinfo_ptr(oid)?;
        reader
            .read_seq_bytes_io()
            .map_err(|error| error.to_string())?
            .map(|entry| entry.sequence)
            .ok_or_else(|| "OId out of bounds.".into())
    }
    pub fn seq_length(&self, oid: OId) -> Result<Loc, String> {
        self.seq_lengths
            .get(oid as usize)
            .copied()
            .ok_or_else(|| "std::exception".into())
    }
    /// Rust offset-seek extension (Source's public seek_offset is a no-op).
    pub fn seek_offset(&mut self, offset: u64) -> Result<(), String> {
        if self.paths.len() != 1 {
            return Err("Offset seeking requires one FASTA file.".into());
        }
        let result = self
            .files
            .get_mut()
            .map_err(|_| io::Error::other("Poisoned FASTA input owner."))
            .and_then(|files| files[0].seek(offset));
        self.errors.capture(result)
    }
    /// Rust partition extension: Source FASTA partition methods are unsupported.
    pub fn create_partition_balanced(&mut self, max_letters: i64) -> Result<(), String> {
        if max_letters < 0 {
            return Err("Maximum partition letters must be non-negative.".into());
        }
        let mut reader = Self::open_source_with_buffer(
            &self.paths,
            FastaFileConfig {
                sequence_type: self.sequence_type,
                ..Default::default()
            },
            self.buffer_size,
        )
        .map_err(|error| error.to_string())?;
        let mut lengths = Vec::new();
        while let Some(record) = reader.read_seq()? {
            lengths.push(record.sequence.len());
        }
        self.partition.clear();
        let (mut begin, mut letters, mut count) = (0, 0, 0);
        for (oid, len) in lengths.iter().enumerate() {
            if count == 0 {
                begin = oid;
            }
            letters += len;
            count += 1;
            if letters > max_letters as usize || oid + 1 == lengths.len() {
                self.partition
                    .push(Chunk::new(self.partition.len() as i32, begin, count));
                letters = 0;
                count = 0;
            }
        }
        self.partition.reverse();
        Ok(())
    }
    pub fn get_n_partition_chunks(&self) -> i32 {
        self.partition.len() as i32
    }
    pub fn partitions(&self) -> &[Chunk] {
        &self.partition
    }
    /// Rust extension; Source FASTA seek_chunk is unsupported.
    pub fn seek_chunk(&mut self, chunk: &Chunk) -> Result<(), String> {
        self.set_seqinfo_ptr(chunk.offset as OId)
    }
    pub fn save_partition(&self, path: impl AsRef<Path>, annotation: &str) -> Result<(), String> {
        let mut out = BufWriter::new(File::create(path).map_err(|error| error.to_string())?);
        for chunk in &self.partition {
            write!(out, "{} {} {}", chunk.i, chunk.offset, chunk.n_seqs)
                .map_err(|error| error.to_string())?;
            if !annotation.is_empty() {
                write!(out, " {annotation}").map_err(|error| error.to_string())?;
            }
            writeln!(out).map_err(|error| error.to_string())?;
        }
        Ok(())
    }
    pub fn init_write(&mut self) -> Result<(), String> {
        let result = self
            .writer
            .as_mut()
            .ok_or_else(|| io::Error::other("FASTA file is not open for writing."))
            .and_then(|writer| {
                writer
                    .output()
                    .seek(0, std::io::SeekFrom::End(0))
                    .map_err(io::Error::other)
            });
        self.errors.capture(result)
    }

    pub fn write_seq(&mut self, sequence: &[Letter], id: &str) -> Result<(), String> {
        self.write_seq_bytes(sequence, id.as_bytes())
    }

    pub fn write_seq_bytes(&mut self, sequence: &[Letter], id: &[u8]) -> Result<(), String> {
        let alphabet = match self.sequence_type {
            SequenceType::AminoAcid => AMINO_ACID_ALPHABET,
            SequenceType::Nucleotide => NUCLEOTIDE_ALPHABET,
        };
        let mut bytes = Vec::new();
        bytes.push(b'>');
        bytes.extend_from_slice(crate::util::sequence::c_string_bytes(id));
        bytes.push(b'\n');
        for chunk in sequence.chunks(160) {
            for &letter in chunk {
                bytes.push(
                    *alphabet
                        .get((letter as u8 & 0x7f) as usize)
                        .ok_or("Invalid encoded sequence letter.")?,
                );
            }
            bytes.push(b'\n');
        }
        let result = self
            .writer
            .as_mut()
            .ok_or_else(|| io::Error::other("FASTA file is not open for writing."))
            .and_then(|writer| writer.output().write_raw(&bytes).map_err(io::Error::other));
        self.errors.capture(result)?;
        self.seqs += 1;
        self.letters += sequence.len() as u64;
        if self.flags.contains(SequenceFileFlags::NEED_LENGTH_LOOKUP) {
            self.seq_lengths.push(sequence.len() as Loc);
        }
        Ok(())
    }
    pub fn close_io(&mut self) -> io::Result<()> {
        for file in self
            .files
            .get_mut()
            .map_err(|_| io::Error::other("Poisoned FASTA input owner."))?
        {
            file.close()?;
        }
        Ok(())
    }
    pub fn close(&mut self) -> Result<(), String> {
        let result = self.close_io();
        self.errors.capture(result)
    }
    pub fn index(path: impl AsRef<Path>, destination: impl AsRef<Path>) -> Result<(), String> {
        Self::index_io(path.as_ref(), destination.as_ref()).map_err(|error| error.to_string())
    }
    pub fn index_io(path: &Path, destination: &Path) -> io::Result<()> {
        // Source opens its unchecked ofstream BEFORE checking the ifstream.
        let mut output = File::create(destination).ok().map(BufWriter::new);
        let file = File::open(path).map_err(|_| {
            io::Error::other(crate::util::io::IoError::path_message(
                b"Failed to open file: ",
                path,
            ))
        })?;
        crate::data::fasta::read_fasta_bytes(&mut BufReader::new(file), |_, sequence, offset| {
            if let Some(out) = &mut output {
                let _ = writeln!(out, "{offset}\t{}", sequence.len());
            }
            Ok(())
        })
    }
}
impl Drop for FastaFile {
    fn drop(&mut self) {
        if self.errors.source_destructors() {
            if let Err(error) = self.errors.emit_pending_perror() {
                crate::util::io::abort_source_diagnostic_error(&error);
            }
        }
        if let Err(error) = self.close_io() {
            if self.errors.source_destructors() {
                if let Err(conversion) = crate::util::io::emit_source_perror(&error) {
                    crate::util::io::abort_source_diagnostic_error(&conversion);
                }
                crate::util::io::abort_source_process();
            }
            let _ = self.errors.capture::<()>(Err(error));
        }
    }
}

fn encoding(sequence_type: SequenceType) -> CharRepresentation {
    match sequence_type {
        SequenceType::AminoAcid => {
            CharRepresentation::new(AMINO_ACID_ALPHABET, MASK_LETTER, b"UO-")
        }
        SequenceType::Nucleotide => CharRepresentation::new(NUCLEOTIDE_ALPHABET, 4, b"MRWSYKVHDBX"),
    }
}
pub fn guess_format(bytes: &[u8]) -> Result<SeqFileFormat, String> {
    match bytes.first() { Some(b'>')=>Ok(SeqFileFormat::Fasta),Some(b'@')=>Ok(SeqFileFormat::Fastq),None=>Err("Error detecting input file format. Input file seems to be empty.".into()),_=>Err("Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).".into()) }
}
fn read_index(path: &Path) -> io::Result<Vec<(u64, Loc)>> {
    let file = File::open(crate::util::io::source_c_path(path).as_ref()).map_err(|_| {
        io::Error::other(crate::util::io::IoError::path_message(
            b"Error opening file: ",
            path,
        ))
    })?;
    let mut input = BufReader::new(file);
    // C++ operator>> tokenizes whitespace, not physical lines. Keep only the
    // metadata vector, including multiple offset/length pairs on one line.
    fn token(input: &mut impl BufRead) -> io::Result<Option<Vec<u8>>> {
        let mut value = Vec::new();
        loop {
            let buffer = input.fill_buf()?;
            if buffer.is_empty() {
                return Ok((!value.is_empty()).then_some(value));
            }
            let byte = buffer[0];
            input.consume(1);
            if byte.is_ascii_whitespace() {
                if !value.is_empty() {
                    return Ok(Some(value));
                }
            } else {
                value.push(byte);
            }
        }
    }
    let mut entries = Vec::new();
    while let Some(offset) = token(&mut input)? {
        let Some(offset) = std::str::from_utf8(&offset)
            .ok()
            .and_then(|value| value.parse::<u64>().ok())
        else {
            break;
        };
        let length = token(&mut input)?.and_then(|value| {
            std::str::from_utf8(&value)
                .ok()
                .and_then(|value| value.parse::<Loc>().ok())
        });
        // Keep malformed index input at an explicit checked boundary; the
        // Source-formatted records above retain their physical byte offsets.
        let Some(length) = length else {
            return Err(io::Error::other("Invalid FASTA index"));
        };
        entries.push((offset, length));
    }
    Ok(entries)
}
#[cfg(test)]
mod tests {
    use super::*;

    fn temp(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!("diamond-rs-{name}-{}", std::process::id()))
    }

    #[test]
    fn fasta_offsets_index_and_random_access() {
        let fasta = temp("fasta-file.fa");
        let index = temp("fasta-file.fai");
        std::fs::write(&fasta, b">one full title\nARND\n>two\nCQ\n").unwrap();
        FastaFile::index(&fasta, &index).unwrap();
        assert_eq!(std::fs::read_to_string(&index).unwrap(), "0\t4\n21\t2\n");
        let mut file = FastaFile::open(
            &[fasta.clone()],
            FastaFileConfig {
                index_file: Some(index.clone()),
                flags: SequenceFileFlags::NEED_LETTER_COUNT | SequenceFileFlags::NEED_LENGTH_LOOKUP,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!((file.sequence_count(), file.letters()), (2, 6));
        file.set_seqinfo_ptr(1).unwrap();
        assert_eq!(file.tell_seq(), 0); // Source indexed seek does not assign oid_.
        let record = file.read_seq().unwrap().unwrap();
        assert_eq!(file.tell_seq(), 1);
        assert_eq!(record.id, "two");
        assert_eq!(record.sequence, vec![4, 5]);
        let _ = std::fs::remove_file(fasta);
        let _ = std::fs::remove_file(index);
    }

    #[test]
    fn fastq_multiline_sequence_reaches_source_schema_error() {
        let path = temp("fasta-file.fq");
        std::fs::write(&path, b"@read one\nAR\nND\n+\n!!\n##\n").unwrap();
        let mut file = FastaFile::open(&[path.clone()], FastaFileConfig::default()).unwrap();
        assert_eq!(file.read_seq().unwrap_err(), "Malformed FASTQ record");
        file.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn paired_files_are_round_robin_and_checked() {
        let a = temp("paired-a.fa");
        let b = temp("paired-b.fa");
        std::fs::write(&a, b">a1\nAA\n>a2\nRR\n").unwrap();
        std::fs::write(&b, b">b1\nNN\n>b2\nDD\n").unwrap();
        let mut file =
            FastaFile::open(&[a.clone(), b.clone()], FastaFileConfig::default()).unwrap();
        let ids: Vec<String> = (0..4)
            .map(|_| file.read_seq().unwrap().unwrap().id)
            .collect();
        assert_eq!(ids, ["a1", "b1", "a2", "b2"]);
        assert!(file.files_synced());
        let _ = std::fs::remove_file(a);
        let _ = std::fs::remove_file(b);
    }

    #[test]
    fn balanced_partition_matches_database_chunk_threshold() {
        let path = temp("partition.fa");
        std::fs::write(&path, b">a\nAAAA\n>b\nRRRRR\n>c\nNN\n").unwrap();
        let mut file = FastaFile::open(&[path.clone()], FastaFileConfig::default()).unwrap();
        file.create_partition_balanced(6).unwrap();
        assert_eq!(
            file.partitions(),
            &[Chunk::new(1, 2, 1), Chunk::new(0, 0, 2)]
        );
        let chunk = file.partitions()[1];
        file.seek_chunk(&chunk).unwrap();
        assert_eq!(file.read_seq().unwrap().unwrap().id, "a");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn writer_reads_same_file_before_close_and_deletes_source_destination() {
        let path = temp("writer.fa");
        let mut writer = FastaFile::create(&path, true, SequenceType::AminoAcid).unwrap();
        writer.write_seq(&[0; 161], "id title").unwrap();
        assert_eq!((writer.sequence_count(), writer.letters()), (1, 161));
        // Source out_file_->rewind flushes outer serializer before the shared
        // input reader sees the body. Closing alone does not perform this.
        writer.init_seq_access().unwrap();
        let record = writer.read_seq().unwrap().unwrap();
        assert_eq!(record.id, "id title");
        assert_eq!(record.sequence, vec![0; 161]);
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(
            bytes.split(|&byte| byte == b'\n').nth(1).unwrap().len(),
            160
        );
        writer.close().unwrap();
        assert!(!path.exists());
    }

    #[test]
    fn constructor_retains_stream_and_defers_later_sequence_error() {
        let path = temp("streaming-late-error.fa");
        std::fs::write(&path, b">first\nAA\n>second\n!\n").unwrap();
        let mut file = FastaFile::open(&[path.clone()], Default::default()).unwrap();
        assert_eq!((file.sequence_count(), file.letters()), (0, 0));
        assert_eq!(file.read_seq().unwrap().unwrap().id, "first");
        assert!(file
            .read_seq()
            .unwrap_err()
            .contains("Invalid character in sequence"));
        assert_eq!(file.tell_seq(), 2);
        file.close().unwrap();
        assert!(FastaFile::open(
            &[path.clone()],
            FastaFileConfig {
                flags: SequenceFileFlags::NEED_LETTER_COUNT,
                ..Default::default()
            }
        )
        .is_err());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn source_byte_records_do_not_require_utf8_headers() {
        let path = temp("streaming-byte-title.fa");
        std::fs::write(&path, b">\xff title\nAA\n").unwrap();
        let mut file =
            FastaFile::open_source_with_buffer(&[path.clone()], Default::default(), 64).unwrap();
        let record = file.read_seq_bytes_io().unwrap().unwrap();
        assert_eq!(record.id, b"\xff title");
        assert_eq!(record.sequence, vec![0, 0]);
        file.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn fasta_raw_eof_does_not_trigger_source_second_file_probe() {
        let a = temp("paired-probe-a.fa");
        let b = temp("paired-probe-b.fa");
        std::fs::write(&a, b">a1\nAA\n").unwrap();
        std::fs::write(&b, b">b1\nRR\n>b2\nNN\n").unwrap();
        let mut file = FastaFile::open(&[a.clone(), b.clone()], Default::default()).unwrap();
        assert_eq!(file.read_seq().unwrap().unwrap().id, "a1");
        assert_eq!(file.read_seq().unwrap().unwrap().id, "b1");
        // Source FastaTokenizer never updates TextInputFile::eof_. Even at
        // physical EOF files_synced returns true without consuming b2.
        assert!(file.files_synced_io().unwrap());
        assert!(!file.eof());
        assert!(file.read_seq().unwrap().is_none());
        assert_eq!(file.tell_seq(), 3);
        assert!(file.files_synced_io().unwrap());
        assert!(!file.eof());
        file.set_seqinfo_ptr(5).unwrap();
        assert_eq!(file.tell_seq(), 5);
        assert!(!file.eof());
        file.close().unwrap();
        std::fs::remove_file(a).unwrap();
        std::fs::remove_file(b).unwrap();
    }

    #[test]
    fn fastq_getline_eof_triggers_source_second_file_probe() {
        // Exact paired FASTQ data and operation order from the pinned Source
        // fasta-class-source-controls/proof.json (paired-fastq case).
        let a = temp("paired-probe-a.fq");
        let b = temp("paired-probe-b.fq");
        std::fs::write(&a, b"@a1\nARND\n+\nIIII\n").unwrap();
        std::fs::write(&b, b"@b1\nCQE\n+\nJJJ\n@b2\nHILK\n+\nKKKK\n").unwrap();
        let mut file = FastaFile::open(&[a.clone(), b.clone()], Default::default()).unwrap();
        assert_eq!(file.read_seq().unwrap().unwrap().id, "a1");
        assert_eq!(file.read_seq().unwrap().unwrap().id, "b1");
        assert!(file.read_seq().unwrap().is_none());
        assert_eq!(file.tell_seq(), 3);
        assert!(!file.eof());
        assert!(!file.files_synced_io().unwrap()); // Consumes b2, not its EOF.
        assert!(!file.eof());
        assert!(file.files_synced_io().unwrap()); // getline then observes EOF.
        assert!(file.eof());
        assert_eq!(file.tell_seq(), 3); // Probes do not increment FastaFile OID.
        assert!(file.read_seq().unwrap().is_none());
        assert_eq!(file.tell_seq(), 4);
        file.close().unwrap();
        std::fs::remove_file(a).unwrap();
        std::fs::remove_file(b).unwrap();
    }

    #[test]
    fn paired_second_format_is_not_guessed_at_construction() {
        let a = temp("paired-format-a.fa");
        let b = temp("paired-format-b.fq");
        std::fs::write(&a, b">a1\nAA\n").unwrap();
        std::fs::write(&b, b"@b1\nRR\n+\n!!\n").unwrap();
        let mut file = FastaFile::open(&[a.clone(), b.clone()], Default::default()).unwrap();
        assert!(file.read_seq().unwrap().is_some());
        assert_eq!(file.read_seq().unwrap_err(), "Missing fields in input line");
        file.close().unwrap();
        std::fs::remove_file(a).unwrap();
        std::fs::remove_file(b).unwrap();
    }

    #[test]
    fn readwrite_constructor_does_not_create_a_missing_file() {
        let path = temp("nonexistent-readwrite.fa");
        let _ = std::fs::remove_file(&path);
        assert!(FastaFile::create(&path, false, SequenceType::AminoAcid).is_err());
        assert!(!path.exists());
    }

    #[test]
    fn source_fastq_schema_and_multiline_quality_use_three_fields() {
        let path = temp("source-fastq-quality.fq");
        std::fs::write(&path, b"@read\nARND\n+descriptor\n!!\n##\n").unwrap();
        let mut file = FastaFile::open(&[path.clone()], Default::default()).unwrap();
        let record = file.read_seq().unwrap().unwrap();
        assert_eq!(record.sequence, vec![0, 1, 2, 3]);
        assert_eq!(record.quality.unwrap(), b"!!##");
        file.close().unwrap();
        std::fs::remove_file(path).unwrap();
        let (id, sequence, quality) = source_fastq_fields(b"@id\n@ARND\nRRRR\n+\n!!!!\n").unwrap();
        assert_eq!(id, b"id");
        assert_eq!(sequence, b"ARND");
        assert_eq!(quality, b"RRRR");
        assert_eq!(
            source_fastq_fields(b"@id\n+\n!!!!\n")
                .unwrap_err()
                .to_string(),
            "Missing fields in input line"
        );
    }

    #[test]
    fn metadata_mapping_callback_runs_in_actual_counting_pass_only() {
        let path = temp("source-mapping-phase.fa");
        std::fs::write(&path, b">one title\nAA\n>two\nRRR\n").unwrap();
        let mut observed = Vec::new();
        let mut file = FastaFile::open_source_with_mapping(
            &[path.clone()],
            FastaFileConfig {
                flags: SequenceFileFlags::NEED_LETTER_COUNT | SequenceFileFlags::OID_TO_ACC_MAPPING,
                ..Default::default()
            },
            64,
            &mut |title, oid| {
                observed.push((title.to_vec(), oid));
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(
            observed,
            vec![(b"one title".to_vec(), 0), (b"two".to_vec(), 1)]
        );
        assert_eq!(
            (file.sequence_count(), file.letters(), file.tell_seq()),
            (2, 5, 0)
        );
        file.read_seq().unwrap();
        assert_eq!(observed.len(), 2); // Ordinary reads do not republish mappings.
        file.close().unwrap();
        observed.clear();
        let mut file = FastaFile::open_source_with_mapping(
            &[path.clone()],
            FastaFileConfig {
                flags: SequenceFileFlags::OID_TO_ACC_MAPPING,
                ..Default::default()
            },
            64,
            &mut |title, oid| {
                observed.push((title.to_vec(), oid));
                Ok(())
            },
        )
        .unwrap();
        assert!(observed.is_empty());
        file.close().unwrap();
        std::fs::write(&path, b">one title\nAA\n>two\n!\n").unwrap();
        let result = FastaFile::open_source_with_mapping(
            &[path.clone()],
            FastaFileConfig {
                flags: SequenceFileFlags::NEED_LETTER_COUNT | SequenceFileFlags::OID_TO_ACC_MAPPING,
                ..Default::default()
            },
            64,
            &mut |title, oid| {
                observed.push((title.to_vec(), oid));
                Ok(())
            },
        );
        assert!(result.is_err());
        assert_eq!(observed, vec![(b"one title".to_vec(), 0)]);
        std::fs::remove_file(path).unwrap();
    }
}
