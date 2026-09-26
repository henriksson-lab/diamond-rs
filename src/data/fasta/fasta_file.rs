//! `diamond/src/data/fasta/fasta_file.{h,cpp}` translated as a concrete,
//! stateful FASTA/FASTQ sequence source.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::basic::value::{
    CharRepresentation, Letter, Loc, OId, SequenceType, AMINO_ACID_ALPHABET, MASK_LETTER,
    NUCLEOTIDE_ALPHABET,
};
use crate::data::fasta::SeqFileFormat;
use crate::data::sequence_file::{Chunk, SequenceFileFlags};

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

/// Concrete FASTA reader/writer. Multiple inputs are consumed round-robin,
/// matching C++ paired-file behavior.
pub struct FastaFile {
    paths: Vec<PathBuf>,
    files: Vec<Vec<FastaEntry>>,
    positions: Vec<usize>,
    file_ptr: usize,
    oid: OId,
    format: SeqFileFormat,
    seqs: u64,
    letters: u64,
    seq_lengths: Vec<Loc>,
    offsets: Vec<u64>,
    partition: Vec<Chunk>,
    writer: Option<BufWriter<File>>,
    sequence_type: SequenceType,
}

impl FastaFile {
    pub fn open(paths: &[PathBuf], config: FastaFileConfig) -> Result<Self, String> {
        if paths.is_empty() {
            return Err("Missing FASTA file name.".into());
        }
        if paths.len() > 2 {
            return Err("FASTA input supports at most two synchronized files.".into());
        }
        let taxon_flags = SequenceFileFlags::TAXON_MAPPING
            | SequenceFileFlags::TAXON_NODES
            | SequenceFileFlags::TAXON_RANKS
            | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES;
        if config.flags.contains(taxon_flags) {
            return Err("Fasta database format does not support taxonomic features.".into());
        }

        let mut files = Vec::with_capacity(paths.len());
        let mut format = None;
        for path in paths {
            let bytes = read_input(path)?;
            let detected = guess_format(&bytes)?;
            if format.is_some_and(|f| f != detected) {
                return Err("Synchronized sequence files must use the same format.".into());
            }
            format = Some(detected);
            files.push(parse_entries(&bytes, detected, config.sequence_type)?);
        }

        if files.len() == 2 && files[0].len() != files[1].len() {
            return Err("Synchronized sequence files contain different record counts.".into());
        }
        let seqs = files.iter().map(Vec::len).sum::<usize>() as u64;
        let letters = files
            .iter()
            .flatten()
            .map(|entry| entry.sequence.len() as u64)
            .sum();
        let mut seq_lengths = Vec::with_capacity(seqs as usize);
        let records_per_file = files.first().map_or(0, Vec::len);
        for record in 0..records_per_file {
            for entries in &files {
                seq_lengths.push(entries[record].sequence.len() as Loc);
            }
        }
        let offsets = if let Some(index) = &config.index_file {
            if paths.len() != 1 {
                return Err("A FASTA index can only be used with one input file.".into());
            }
            let indexed = read_index(index)?;
            if indexed.len() != files[0].len() {
                return Err("FASTA index record count does not match the input file.".into());
            }
            indexed.iter().map(|x| x.0).collect()
        } else if paths.len() == 1 {
            files[0].iter().map(|entry| entry.offset).collect()
        } else {
            Vec::new()
        };

        Ok(Self {
            paths: paths.to_vec(),
            positions: vec![0; files.len()],
            files,
            file_ptr: 0,
            oid: 0,
            format: format.expect("non-empty paths"),
            seqs,
            letters,
            seq_lengths,
            offsets,
            partition: Vec::new(),
            writer: None,
            sequence_type: config.sequence_type,
        })
    }

    pub fn create(
        path: impl AsRef<Path>,
        overwrite: bool,
        sequence_type: SequenceType,
    ) -> Result<Self, String> {
        let path = path.as_ref();
        let file = if overwrite {
            File::create(path)
        } else {
            std::fs::OpenOptions::new()
                .read(true)
                .append(true)
                .create(true)
                .open(path)
        }
        .map_err(|e| e.to_string())?;
        Ok(Self {
            paths: vec![path.to_path_buf()],
            files: vec![Vec::new()],
            positions: vec![0],
            file_ptr: 0,
            oid: 0,
            format: SeqFileFormat::Fasta,
            seqs: 0,
            letters: 0,
            seq_lengths: Vec::new(),
            offsets: Vec::new(),
            partition: Vec::new(),
            writer: Some(BufWriter::new(file)),
            sequence_type,
        })
    }

    pub fn file_count(&self) -> i64 {
        self.files.len() as i64
    }
    pub fn files_synced(&self) -> bool {
        self.positions.windows(2).all(|w| w[0] == w[1])
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
        self.positions
            .iter()
            .zip(&self.files)
            .all(|(p, f)| *p >= f.len())
    }
    pub fn file_name(&self) -> Option<&Path> {
        self.paths.first().map(PathBuf::as_path)
    }

    pub fn init_seq_access(&mut self) {
        self.set_seqinfo_ptr(0).expect("OID zero is valid");
    }

    pub fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String> {
        if oid > self.seqs {
            return Err("FastaFile::set_seqinfo_ptr".into());
        }
        self.positions.fill(0);
        self.file_ptr = 0;
        self.oid = 0;
        for _ in 0..oid {
            if self.read_seq()?.is_none() {
                return Err("FastaFile::set_seqinfo_ptr".into());
            }
        }
        Ok(())
    }

    pub fn read_seq(&mut self) -> Result<Option<FastaEntry>, String> {
        if self.files.is_empty() {
            return Ok(None);
        }
        let file = self.file_ptr;
        let position = self.positions[file];
        if position >= self.files[file].len() {
            return Ok(None);
        }
        let entry = self.files[file][position].clone();
        self.positions[file] += 1;
        self.file_ptr = (self.file_ptr + 1) % self.files.len();
        self.oid += 1;
        Ok(Some(entry))
    }

    pub fn seq_data(&self, oid: OId) -> Result<Vec<Letter>, String> {
        self.entry_by_oid(oid).map(|entry| entry.sequence.clone())
    }

    pub fn seq_length(&self, oid: OId) -> Result<Loc, String> {
        self.seq_lengths
            .get(oid as usize)
            .copied()
            .ok_or_else(|| "Sequence length lookup not available.".into())
    }

    pub fn seek_offset(&mut self, offset: u64) -> Result<(), String> {
        if self.files.len() != 1 {
            return Err("Offset seeking requires one FASTA file.".into());
        }
        let oid = self
            .offsets
            .binary_search(&offset)
            .map_err(|_| "FASTA offset not found.".to_string())?;
        self.set_seqinfo_ptr(oid as OId)
    }

    pub fn create_partition_balanced(&mut self, max_letters: i64) -> Result<(), String> {
        if max_letters < 0 {
            return Err("Maximum partition letters must be non-negative.".into());
        }
        self.partition.clear();
        let max_letters = max_letters as u64;
        let mut begin = 0usize;
        let mut letters = 0u64;
        let mut count = 0usize;
        for (oid, &len) in self.seq_lengths.iter().enumerate() {
            if count == 0 {
                begin = oid;
            }
            letters += len.max(0) as u64;
            count += 1;
            if letters > max_letters || oid + 1 == self.seq_lengths.len() {
                self.partition
                    .push(Chunk::new(self.partition.len() as i32, begin, count as i64));
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

    pub fn seek_chunk(&mut self, chunk: &Chunk) -> Result<(), String> {
        self.set_seqinfo_ptr(chunk.offset as OId)
    }

    pub fn save_partition(&self, path: impl AsRef<Path>, annotation: &str) -> Result<(), String> {
        let mut out = BufWriter::new(File::create(path).map_err(|e| e.to_string())?);
        for chunk in &self.partition {
            write!(out, "{} {} {}", chunk.i, chunk.offset, chunk.n_seqs)
                .map_err(|e| e.to_string())?;
            if !annotation.is_empty() {
                write!(out, " {annotation}").map_err(|e| e.to_string())?;
            }
            writeln!(out).map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    pub fn write_seq(&mut self, sequence: &[Letter], id: &str) -> Result<(), String> {
        let alphabet = match self.sequence_type {
            SequenceType::AminoAcid => AMINO_ACID_ALPHABET,
            SequenceType::Nucleotide => NUCLEOTIDE_ALPHABET,
        };
        let writer = self
            .writer
            .as_mut()
            .ok_or("FASTA file is not open for writing.")?;
        writeln!(writer, ">{id}").map_err(|e| e.to_string())?;
        for chunk in sequence.chunks(80) {
            for &letter in chunk {
                let byte = alphabet
                    .get((letter as u8 & 0x7f) as usize)
                    .ok_or("Invalid encoded sequence letter.")?;
                writer.write_all(&[*byte]).map_err(|e| e.to_string())?;
            }
            writer.write_all(b"\n").map_err(|e| e.to_string())?;
        }
        writer.flush().map_err(|e| e.to_string())?;
        self.seqs += 1;
        self.letters += sequence.len() as u64;
        self.seq_lengths.push(sequence.len() as Loc);
        Ok(())
    }

    pub fn close(&mut self) -> Result<(), String> {
        if let Some(writer) = &mut self.writer {
            writer.flush().map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    pub fn index(path: impl AsRef<Path>, destination: impl AsRef<Path>) -> Result<(), String> {
        let bytes = std::fs::read(path).map_err(|e| e.to_string())?;
        let mut output = BufWriter::new(File::create(destination).map_err(|e| e.to_string())?);
        let mut input = BufReader::new(bytes.as_slice());
        crate::data::fasta::read_fasta(&mut input, |_id, sequence, offset| {
            writeln!(output, "{offset}\t{}", sequence.len())
        })
        .map_err(|e| e.to_string())?;
        Ok(())
    }

    fn entry_by_oid(&self, oid: OId) -> Result<&FastaEntry, String> {
        if self.files.len() == 1 {
            return self.files[0]
                .get(oid as usize)
                .ok_or_else(|| "OId out of bounds.".into());
        }
        let file = oid as usize % self.files.len();
        let record = oid as usize / self.files.len();
        self.files[file]
            .get(record)
            .ok_or_else(|| "OId out of bounds.".into())
    }
}

pub fn guess_format(bytes: &[u8]) -> Result<SeqFileFormat, String> {
    match bytes.first() {
        None => Err("Error detecting input file format. Input file seems to be empty.".into()),
        Some(b'>') => Ok(SeqFileFormat::Fasta),
        Some(b'@') => Ok(SeqFileFormat::Fastq),
        _ => Err("Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).".into()),
    }
}

fn read_input(path: &Path) -> Result<Vec<u8>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut decoded = Vec::new();
        flate2::read::MultiGzDecoder::new(bytes.as_slice())
            .read_to_end(&mut decoded)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        Ok(decoded)
    } else {
        Ok(bytes)
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

fn parse_entries(
    bytes: &[u8],
    format: SeqFileFormat,
    sequence_type: SequenceType,
) -> Result<Vec<FastaEntry>, String> {
    let text = std::str::from_utf8(bytes).map_err(|e| e.to_string())?;
    let converter = encoding(sequence_type);
    match format {
        SeqFileFormat::Fasta => parse_fasta(text, &converter),
        SeqFileFormat::Fastq => parse_fastq(text, &converter),
    }
}

fn parse_fasta(text: &str, converter: &CharRepresentation) -> Result<Vec<FastaEntry>, String> {
    let mut entries = Vec::new();
    let mut current: Option<FastaEntry> = None;
    let mut offset = 0usize;
    for line_with_newline in text.split_inclusive('\n') {
        let line = line_with_newline.trim_end_matches(['\n', '\r']);
        if let Some(id) = line.strip_prefix('>') {
            if let Some(entry) = current.take() {
                if entry.sequence.is_empty() {
                    return Err("Missing fields in input line".into());
                }
                entries.push(entry);
            }
            if id.is_empty() {
                return Err(format!(
                    "FASTA format error: empty id at file offset {offset}"
                ));
            }
            current = Some(FastaEntry {
                id: id.to_string(),
                sequence: Vec::new(),
                quality: None,
                offset: offset as u64,
            });
        } else if !line.is_empty() {
            let entry = current
                .as_mut()
                .ok_or("FASTA format error: file does not start with '>'")?;
            for byte in line.bytes() {
                entry.sequence.push(converter.convert(byte)?);
            }
        }
        offset += line_with_newline.len();
    }
    if let Some(entry) = current {
        if entry.sequence.is_empty() {
            return Err("Missing fields in input line".into());
        }
        entries.push(entry);
    }
    Ok(entries)
}

fn parse_fastq(text: &str, converter: &CharRepresentation) -> Result<Vec<FastaEntry>, String> {
    let lines: Vec<&str> = text.split_inclusive('\n').collect();
    let mut entries = Vec::new();
    let mut i = 0usize;
    let mut offset = 0usize;
    while i < lines.len() {
        let start = offset;
        let header = lines[i].trim_end_matches(['\n', '\r']);
        let id = header.strip_prefix('@').ok_or("Malformed FASTQ record")?;
        offset += lines[i].len();
        i += 1;
        let mut raw_seq = String::new();
        while i < lines.len() {
            let line = lines[i].trim_end_matches(['\n', '\r']);
            offset += lines[i].len();
            i += 1;
            if line.starts_with('+') {
                break;
            }
            raw_seq.push_str(line);
        }
        if raw_seq.is_empty() {
            return Err("Malformed FASTQ record".into());
        }
        let mut quality = String::new();
        while i < lines.len() && quality.len() < raw_seq.len() {
            let line = lines[i].trim_end_matches(['\n', '\r']);
            quality.push_str(line);
            offset += lines[i].len();
            i += 1;
        }
        if quality.len() < raw_seq.len() {
            return Err("Malformed FASTQ record".into());
        }
        let sequence = raw_seq
            .bytes()
            .map(|b| converter.convert(b))
            .collect::<Result<Vec<_>, _>>()?;
        entries.push(FastaEntry {
            id: id.to_string(),
            sequence,
            quality: Some(quality.into_bytes()),
            offset: start as u64,
        });
    }
    Ok(entries)
}

fn read_index(path: &Path) -> Result<Vec<(u64, Loc)>, String> {
    std::fs::read_to_string(path)
        .map_err(|e| format!("Error opening file: {}: {e}", path.display()))?
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut fields = line.split_whitespace();
            let offset = fields
                .next()
                .ok_or("Invalid FASTA index")?
                .parse::<u64>()
                .map_err(|e| e.to_string())?;
            let len = fields
                .next()
                .ok_or("Invalid FASTA index")?
                .parse::<Loc>()
                .map_err(|e| e.to_string())?;
            Ok((offset, len))
        })
        .collect()
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
                ..Default::default()
            },
        )
        .unwrap();
        file.seek_offset(21).unwrap();
        let record = file.read_seq().unwrap().unwrap();
        assert_eq!(record.id, "two");
        assert_eq!(record.sequence, vec![4, 5]);
        let _ = std::fs::remove_file(fasta);
        let _ = std::fs::remove_file(index);
    }

    #[test]
    fn fastq_multiline_sequence_and_quality() {
        let path = temp("fasta-file.fq");
        std::fs::write(&path, b"@read one\nAR\nND\n+\n!!\n##\n").unwrap();
        let mut file = FastaFile::open(&[path.clone()], FastaFileConfig::default()).unwrap();
        let record = file.read_seq().unwrap().unwrap();
        assert_eq!(record.id, "read one");
        assert_eq!(record.sequence, vec![0, 1, 2, 3]);
        assert_eq!(record.quality.unwrap(), b"!!##");
        let _ = std::fs::remove_file(path);
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
    fn writer_updates_counts_and_roundtrips() {
        let path = temp("writer.fa");
        let mut writer = FastaFile::create(&path, true, SequenceType::AminoAcid).unwrap();
        writer.write_seq(&[0, 1, 2, 3], "id title").unwrap();
        writer.close().unwrap();
        assert_eq!((writer.sequence_count(), writer.letters()), (1, 4));
        let mut reader = FastaFile::open(&[path.clone()], FastaFileConfig::default()).unwrap();
        assert_eq!(reader.read_seq().unwrap().unwrap().id, "id title");
        let _ = std::fs::remove_file(path);
    }
}
