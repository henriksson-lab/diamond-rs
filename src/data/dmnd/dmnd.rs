//! Native DIAMOND database implementation mirroring
//! `diamond/src/data/dmnd/dmnd.{h,cpp}`.

use crate::basic::value::{Letter, Loc, OId, SequenceType, TaxId, DELIMITER_LETTER};
use crate::data::taxon_list::TaxonList;
use crate::data::taxonomy_nodes::TaxonomyNodes;
use crate::util::io::{Deserializer, VecStream};
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

/// Magic number identifying a DIAMOND database file.
pub const MAGIC_NUMBER: u64 = 0x24af8a415ee186d;

/// Current database version for protein databases.
pub const CURRENT_DB_VERSION_PROT: u32 = 3;

/// Current database version for nucleotide databases.
pub const CURRENT_DB_VERSION_NUCL: u32 = 4;

/// Minimum compatible database version.
pub const MIN_DB_VERSION: u32 = 2;

/// Minimum compatible build version.
pub const MIN_BUILD_VERSION: u32 = 74;

/// Build version stamped into .dmnd headers. Matches C++ `Const::build_version`
/// in diamond/src/basic/const.h.
pub const BUILD_VERSION: u32 = 178;

/// DIAMOND database file header (first 40 bytes).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceHeader {
    pub magic_number: u64,
    pub build: u32,
    pub db_version: u32,
    pub sequences: u64,
    pub letters: u64,
    pub pos_array_offset: u64,
}

impl Default for ReferenceHeader {
    fn default() -> Self {
        ReferenceHeader {
            magic_number: MAGIC_NUMBER,
            build: BUILD_VERSION,
            db_version: CURRENT_DB_VERSION_PROT,
            sequences: 0,
            letters: 0,
            pos_array_offset: 0,
        }
    }
}

impl ReferenceHeader {
    pub const SIZE: usize = 40;
    pub fn new() -> Self {
        Self::default()
    }

    /// Read header from a reader.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf = [0u8; Self::SIZE];
        reader.read_exact(&mut buf)?;
        Ok(ReferenceHeader {
            magic_number: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            build: u32::from_le_bytes(buf[8..12].try_into().unwrap()),
            db_version: u32::from_le_bytes(buf[12..16].try_into().unwrap()),
            sequences: u64::from_le_bytes(buf[16..24].try_into().unwrap()),
            letters: u64::from_le_bytes(buf[24..32].try_into().unwrap()),
            pos_array_offset: u64::from_le_bytes(buf[32..40].try_into().unwrap()),
        })
    }

    /// Write header to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.magic_number.to_le_bytes())?;
        writer.write_all(&self.build.to_le_bytes())?;
        writer.write_all(&self.db_version.to_le_bytes())?;
        writer.write_all(&self.sequences.to_le_bytes())?;
        writer.write_all(&self.letters.to_le_bytes())?;
        writer.write_all(&self.pos_array_offset.to_le_bytes())?;
        Ok(())
    }

    /// Validate this header is a valid DIAMOND database.
    ///
    /// Matches C++ `ReferenceHeader::init()`.
    ///
    // C++ `init()` (`diamond/src/data/dmnd/dmnd.cpp:159-164`) also rejects
    // too-old builds (pre-v74 layout differs), too-new db versions
    // (forward-incompatible writes), and the `sequences == 0` sentinel that
    // marks an incomplete database build.
    pub fn validate(&self) -> Result<(), String> {
        if self.magic_number != MAGIC_NUMBER {
            return Err("Database file is not a DIAMOND database.".into());
        }
        if self.db_version < MIN_DB_VERSION {
            return Err(format!(
                "Database version {} is not supported (minimum: {})",
                self.db_version, MIN_DB_VERSION
            ));
        }
        let max_db_version = CURRENT_DB_VERSION_PROT.max(CURRENT_DB_VERSION_NUCL);
        if self.db_version > max_db_version {
            return Err(format!(
                "Database version {} is newer than this build supports (max: {}). \
                 Please use a newer DIAMOND version.",
                self.db_version, max_db_version
            ));
        }
        if self.build < MIN_BUILD_VERSION {
            return Err(format!(
                "Database was built with build version {} which is too old (minimum: {}). \
                 Please rebuild the database.",
                self.build, MIN_BUILD_VERSION
            ));
        }
        if self.sequences == 0 {
            return Err(
                "Incomplete database file. Database building did not complete successfully.".into(),
            );
        }
        Ok(())
    }
}

/// Extended DIAMOND database header.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReferenceHeader2 {
    pub hash: [u8; 16],
    pub taxon_array_offset: u64,
    pub taxon_array_size: u64,
    pub taxon_nodes_offset: u64,
    pub taxon_names_offset: u64,
}

impl ReferenceHeader2 {
    pub const PAYLOAD_SIZE: usize = 48;
    pub fn new() -> Self {
        Self::default()
    }

    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut hash = [0u8; 16];
        reader.read_exact(&mut hash)?;
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf)?;
        Ok(ReferenceHeader2 {
            hash,
            taxon_array_offset: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            taxon_array_size: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
            taxon_nodes_offset: u64::from_le_bytes(buf[16..24].try_into().unwrap()),
            taxon_names_offset: u64::from_le_bytes(buf[24..32].try_into().unwrap()),
        })
    }

    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.hash)?;
        writer.write_all(&self.taxon_array_offset.to_le_bytes())?;
        writer.write_all(&self.taxon_array_size.to_le_bytes())?;
        writer.write_all(&self.taxon_nodes_offset.to_le_bytes())?;
        writer.write_all(&self.taxon_names_offset.to_le_bytes())?;
        Ok(())
    }

    /// Write the dynamic-record size prefix followed by the payload, matching
    /// C++ `operator<<(Serializer&, const ReferenceHeader2&)`.
    pub fn write_record_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&(Self::PAYLOAD_SIZE as u64).to_le_bytes())?;
        self.write_to(writer)
    }

    /// Read a size-delimited Header2 record and skip any future extension.
    pub fn read_record_from<R: Read + Seek>(reader: &mut R) -> io::Result<Self> {
        let mut size = [0; 8];
        reader.read_exact(&mut size)?;
        let size = u64::from_le_bytes(size);
        if size < Self::PAYLOAD_SIZE as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "ReferenceHeader2 record is too short",
            ));
        }
        let header = Self::read_from(reader)?;
        reader.seek(SeekFrom::Current(size as i64 - Self::PAYLOAD_SIZE as i64))?;
        Ok(header)
    }
}

/// Sequence info entry in the position array (trailer).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[repr(C)]
pub struct SeqInfo {
    pub pos: u64,
    pub seq_len: u32,
    pub padding: u32,
}

impl SeqInfo {
    pub const SIZE: usize = 16;

    pub fn new(pos: u64, seq_len: u32) -> Self {
        Self {
            pos,
            seq_len,
            padding: 0,
        }
    }

    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf = [0u8; 16];
        reader.read_exact(&mut buf)?;
        Ok(SeqInfo {
            pos: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            seq_len: u32::from_le_bytes(buf[8..12].try_into().unwrap()),
            padding: u32::from_le_bytes(buf[12..16].try_into().unwrap()),
        })
    }

    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.pos.to_le_bytes())?;
        writer.write_all(&self.seq_len.to_le_bytes())?;
        // C++ always serializes this reserved word as zero.
        writer.write_all(&0u32.to_le_bytes())
    }
}

/// Check if a file is a DIAMOND database.
pub fn is_diamond_db<R: Read>(reader: &mut R) -> bool {
    let mut buf = [0u8; 8];
    if reader.read_exact(&mut buf).is_err() {
        return false;
    }
    u64::from_le_bytes(buf) == MAGIC_NUMBER
}

pub fn is_diamond_db_file(path: impl AsRef<Path>) -> bool {
    let path = path.as_ref();
    if path.as_os_str().is_empty() || path == Path::new("-") {
        return false;
    }
    File::open(path)
        .map(|mut file| is_diamond_db(&mut file))
        .unwrap_or(false)
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Chunk {
    pub i: i32,
    pub offset: usize,
    pub n_seqs: i64,
}

impl Chunk {
    pub fn new(i: i32, offset: usize, n_seqs: i64) -> Self {
        Self { i, offset, n_seqs }
    }
}

pub fn to_chunk(line: &str) -> Result<Chunk, String> {
    let mut fields = line.split_whitespace();
    let i = fields
        .next()
        .ok_or_else(|| "Missing chunk index".to_string())?
        .parse()
        .map_err(|_| "Invalid chunk index".to_string())?;
    let offset = fields
        .next()
        .ok_or_else(|| "Missing chunk offset".to_string())?
        .parse()
        .map_err(|_| "Invalid chunk offset".to_string())?;
    let n_seqs = fields
        .next()
        .ok_or_else(|| "Missing chunk sequence count".to_string())?
        .parse()
        .map_err(|_| "Invalid chunk sequence count".to_string())?;
    Ok(Chunk::new(i, offset, n_seqs))
}

pub fn chunk_to_string(chunk: &Chunk) -> String {
    format!("{} {} {}", chunk.i, chunk.offset, chunk.n_seqs)
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DatabasePartition {
    pub max_letters: usize,
    pub n_seqs_total: usize,
    pub chunks: Vec<Chunk>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DatabaseOpenOptions {
    pub no_compatibility_check: bool,
    pub load_taxon_mapping: bool,
    pub load_taxon_nodes: bool,
    pub load_taxon_scientific_names: bool,
    pub load_lengths: bool,
}

/// Read-side implementation of C++ `DatabaseFile`.
#[derive(Debug)]
pub struct DatabaseFile {
    path: PathBuf,
    file: File,
    pub ref_header: ReferenceHeader,
    pub header2: ReferenceHeader2,
    pub partition: DatabasePartition,
    pos_array_cursor: OId,
    sequence_cursor: OId,
    data_start: u64,
    seq_length: Vec<Loc>,
    seqids: Vec<String>,
    taxon_list: Option<TaxonList>,
    taxon_nodes: Option<TaxonomyNodes>,
    taxon_scientific_names: Vec<String>,
}

impl DatabaseFile {
    pub const FILE_EXTENSION: &'static str = ".dmnd";

    pub fn open(path: impl AsRef<Path>) -> Result<Self, String> {
        Self::open_with_options(path, DatabaseOpenOptions::default())
    }

    /// Explicit-config counterpart of C++ `DatabaseFile::make_db`, whose
    /// inputs are process-wide `config` fields.
    pub fn make_db(
        input_files: &[&str],
        output_path: &str,
        sequence_type: SequenceType,
    ) -> io::Result<crate::data::db_builder::BuildStats> {
        crate::data::db_builder::build_db(input_files, output_path, sequence_type)
    }

    pub fn read_header<R: Read>(stream: &mut R) -> Result<ReferenceHeader, String> {
        let header = ReferenceHeader::read_from(stream).map_err(|e| e.to_string())?;
        if header.magic_number != MAGIC_NUMBER {
            Err("Database file is not a DIAMOND database.".to_string())
        } else {
            Ok(header)
        }
    }

    pub fn open_with_options(
        path: impl AsRef<Path>,
        options: DatabaseOpenOptions,
    ) -> Result<Self, String> {
        let path = resolve_database_path(path.as_ref());
        let mut file = File::open(&path).map_err(|error| error.to_string())?;
        let ref_header = ReferenceHeader::read_from(&mut file).map_err(|e| e.to_string())?;
        if ref_header.magic_number != MAGIC_NUMBER {
            return Err("Database file is not a DIAMOND database.".to_string());
        }
        if !options.no_compatibility_check {
            ref_header.validate()?;
        }
        let header2 = ReferenceHeader2::read_record_from(&mut file).map_err(|e| e.to_string())?;
        let data_start = file.stream_position().map_err(|e| e.to_string())?;
        let mut db = Self {
            path,
            file,
            ref_header,
            header2,
            partition: DatabasePartition::default(),
            pos_array_cursor: 0,
            sequence_cursor: 0,
            data_start,
            seq_length: Vec::new(),
            seqids: Vec::new(),
            taxon_list: None,
            taxon_nodes: None,
            taxon_scientific_names: Vec::new(),
        };
        db.validate_requested_metadata(options)?;
        db.load_requested_metadata(options)?;
        Ok(db)
    }

    fn validate_requested_metadata(&self, options: DatabaseOpenOptions) -> Result<(), String> {
        let mut missing = Vec::new();
        if options.load_taxon_mapping && !self.has_taxon_id_lists() {
            missing.push("taxonomy mapping information (--taxonmap option)");
        }
        if options.load_taxon_nodes && !self.has_taxon_nodes() {
            missing.push("taxonomy nodes information (--taxonnodes option)");
        }
        if options.load_taxon_scientific_names && !self.has_taxon_scientific_names() {
            missing.push("taxonomy names information (--taxonnames option)");
        }
        if missing.is_empty() {
            Ok(())
        } else {
            Err(format!(
                "Options require taxonomy information included in the database. Please use the respective options to build this information into the database when running diamond makedb: {}",
                missing.join(", ")
            ))
        }
    }

    fn load_requested_metadata(&mut self, options: DatabaseOpenOptions) -> Result<(), String> {
        if options.load_taxon_mapping {
            self.file
                .seek(SeekFrom::Start(self.header2.taxon_array_offset))
                .map_err(|e| e.to_string())?;
            let mut data = vec![0; self.header2.taxon_array_size as usize];
            self.file.read_exact(&mut data).map_err(|e| e.to_string())?;
            self.taxon_list = Some(TaxonList::new(
                data,
                self.ref_header.sequences as usize,
                self.header2.taxon_array_size as usize,
            )?);
        }
        if options.load_taxon_nodes {
            self.file
                .seek(SeekFrom::Start(self.header2.taxon_nodes_offset))
                .map_err(|e| e.to_string())?;
            let file_len = self.file.metadata().map_err(|e| e.to_string())?.len();
            let end = metadata_end(self.header2.taxon_nodes_offset, &self.header2, file_len);
            let mut data = vec![0; (end - self.header2.taxon_nodes_offset) as usize];
            self.file.read_exact(&mut data).map_err(|e| e.to_string())?;
            let mut input = Deserializer::new(VecStream::from_vec(data));
            self.taxon_nodes = Some(
                TaxonomyNodes::from_deserializer(&mut input, self.ref_header.build)
                    .map_err(|e| e.to_string())?,
            );
        }
        if options.load_taxon_scientific_names {
            self.file
                .seek(SeekFrom::Start(self.header2.taxon_names_offset))
                .map_err(|e| e.to_string())?;
            self.taxon_scientific_names = read_string_vector(&mut self.file)?;
        }
        if options.load_lengths {
            self.seq_length = (0..self.ref_header.sequences)
                .map(|oid| self.seq_info_at(oid).map(|info| info.seq_len as Loc))
                .collect::<Result<_, _>>()?;
        }
        Ok(())
    }

    pub fn seq_length_iterator(&self, start: OId) -> Result<SeqLengthIterator, String> {
        SeqLengthIterator::new(
            &self.path,
            self.ref_header.pos_array_offset,
            start,
            self.ref_header.sequences,
        )
    }

    pub fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String> {
        if oid > self.ref_header.sequences {
            return Err("Sequence index out of bounds.".to_string());
        }
        self.pos_array_cursor = oid;
        Ok(())
    }

    pub fn tell_seq(&self) -> OId {
        self.pos_array_cursor
    }

    pub fn eof(&self) -> bool {
        self.tell_seq() == self.sequence_count()
    }

    pub fn read_seqinfo(&mut self) -> Result<SeqInfo, String> {
        let info = self.seq_info_at(self.pos_array_cursor)?;
        self.pos_array_cursor += 1;
        Ok(info)
    }

    pub fn putback_seqinfo(&mut self) -> Result<(), String> {
        if self.pos_array_cursor == 0 {
            return Err("Cannot put back sequence info before the first entry.".to_string());
        }
        self.pos_array_cursor -= 1;
        Ok(())
    }

    fn seq_info_at(&mut self, oid: OId) -> Result<SeqInfo, String> {
        self.file
            .seek(SeekFrom::Start(
                self.ref_header.pos_array_offset + SeqInfo::SIZE as u64 * oid,
            ))
            .map_err(|e| e.to_string())?;
        SeqInfo::read_from(&mut self.file).map_err(|e| e.to_string())
    }

    pub fn init_seq_access(&mut self) {
        self.sequence_cursor = 0;
    }

    pub fn init_seqinfo_access(&mut self) {
        // The Rust reader seeks directly from the logical OID cursor when an
        // entry is read, so no independent buffered-file seek is required.
    }

    pub fn read_seq(&mut self) -> Result<Option<(Vec<Letter>, String)>, String> {
        if self.sequence_cursor >= self.ref_header.sequences {
            return Ok(None);
        }
        let record = self.record(self.sequence_cursor)?;
        self.sequence_cursor += 1;
        Ok(Some(record))
    }

    pub fn skip_seq(&mut self) -> Result<(), String> {
        if self.sequence_cursor >= self.ref_header.sequences {
            return Err("Unexpected end of file.".to_string());
        }
        self.record(self.sequence_cursor)?;
        self.sequence_cursor += 1;
        Ok(())
    }

    fn record(&mut self, oid: OId) -> Result<(Vec<Letter>, String), String> {
        let info = self.seq_info_at(oid)?;
        self.file
            .seek(SeekFrom::Start(info.pos))
            .map_err(|e| e.to_string())?;
        let mut framed = vec![0; info.seq_len as usize + 2];
        self.file
            .read_exact(&mut framed)
            .map_err(|e| e.to_string())?;
        if framed.first() != Some(&0xff) || framed.last() != Some(&0xff) {
            return Err("Invalid sequence framing in DIAMOND database.".to_string());
        }
        let mut id = Vec::new();
        loop {
            let mut byte = [0];
            self.file
                .read_exact(&mut byte)
                .map_err(|_| "Unexpected end of file.".to_string())?;
            if byte[0] == 0 {
                break;
            }
            id.push(byte[0]);
        }
        Ok((
            framed[1..framed.len() - 1]
                .iter()
                .map(|&letter| letter as Letter)
                .collect(),
            String::from_utf8_lossy(&id).into_owned(),
        ))
    }

    pub fn create_partition_fixed_number(&mut self, count: usize) -> Result<(), String> {
        if count == 0 {
            return Err("Partition count must be positive.".to_string());
        }
        let max_letters = (self.ref_header.letters as usize).div_ceil(count);
        self.create_partition(max_letters)
    }

    pub fn create_partition_balanced(&mut self, max_letters: i64) -> Result<(), String> {
        if max_letters < 0 {
            return Err("Maximum partition letters must be non-negative.".to_string());
        }
        self.create_partition(max_letters as usize)
    }

    pub fn create_partition(&mut self, max_letters: usize) -> Result<(), String> {
        let mut chunks = Vec::new();
        let mut letters = 0usize;
        let mut seqs = 0usize;
        let mut begin = 0usize;
        for oid in 0..self.ref_header.sequences as usize {
            let info = self.seq_info_at(oid as OId)?;
            if seqs == 0 {
                begin = oid;
            }
            letters += info.seq_len as usize;
            seqs += 1;
            if letters > max_letters || oid + 1 == self.ref_header.sequences as usize {
                chunks.push(Chunk::new(chunks.len() as i32, begin, seqs as i64));
                letters = 0;
                seqs = 0;
            }
        }
        chunks.reverse();
        self.partition = DatabasePartition {
            max_letters,
            n_seqs_total: self.ref_header.sequences as usize,
            chunks,
        };
        Ok(())
    }

    pub fn get_n_partition_chunks(&self) -> i32 {
        self.partition.chunks.len() as i32
    }

    pub fn save_partition(&self, path: impl AsRef<Path>, annotation: &str) -> Result<(), String> {
        let mut output = File::create(path).map_err(|e| e.to_string())?;
        for chunk in &self.partition.chunks {
            write!(output, "{}", chunk_to_string(chunk)).map_err(|e| e.to_string())?;
            if !annotation.is_empty() {
                write!(output, " {annotation}").map_err(|e| e.to_string())?;
            }
            writeln!(output).map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    pub fn load_partition(&mut self, path: impl AsRef<Path>) -> Result<(), String> {
        self.clear_partition();
        let input = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        for line in input.lines() {
            self.partition.chunks.push(to_chunk(line)?);
        }
        Ok(())
    }

    pub fn clear_partition(&mut self) {
        self.partition = DatabasePartition::default();
    }

    pub fn seek_chunk(&mut self, chunk: &Chunk) -> Result<(), String> {
        self.set_seqinfo_ptr(chunk.offset as OId)
    }

    pub fn seek_offset(&mut self, offset: usize) -> Result<(), String> {
        self.file
            .seek(SeekFrom::Start(offset as u64))
            .map(|_| ())
            .map_err(|e| e.to_string())
    }

    pub fn id_len(seq_info: &SeqInfo, next: &SeqInfo) -> Result<usize, String> {
        next.pos
            .checked_sub(seq_info.pos + seq_info.seq_len as u64 + 3)
            .map(|length| length as usize)
            .ok_or_else(|| "Invalid sequence-info offsets.".to_string())
    }

    pub fn sequence_count(&self) -> OId {
        self.ref_header.sequences
    }

    /// Low-level counterpart of C++ `read_seq_data`. The returned allocation
    /// includes one delimiter byte on either side of the sequence.
    pub fn read_seq_data(
        &mut self,
        len: usize,
        pos: &mut usize,
        seek: bool,
    ) -> Result<Vec<Letter>, String> {
        if seek {
            self.seek_offset(*pos)?;
        }
        let mut bytes = vec![0; len + 2];
        self.file
            .read_exact(&mut bytes)
            .map_err(|e| e.to_string())?;
        bytes[0] = DELIMITER_LETTER as u8;
        bytes[len + 1] = DELIMITER_LETTER as u8;
        Ok(bytes.into_iter().map(|byte| byte as Letter).collect())
    }

    pub fn read_id_data(&mut self, len: usize) -> Result<Vec<u8>, String> {
        let mut id = vec![0; len + 1];
        self.file.read_exact(&mut id).map_err(|e| e.to_string())?;
        Ok(id)
    }

    pub fn skip_id_data(&mut self) -> Result<(), String> {
        loop {
            let mut byte = [0];
            self.file
                .read_exact(&mut byte)
                .map_err(|_| "Unexpected end of file.".to_string())?;
            if byte[0] == 0 {
                return Ok(());
            }
        }
    }

    pub fn letters(&self) -> u64 {
        self.ref_header.letters
    }

    pub fn db_version(&self) -> i32 {
        self.ref_header.db_version as i32
    }

    pub fn program_build_version(&self) -> i32 {
        self.ref_header.build as i32
    }

    pub fn build_version(&self) -> i32 {
        self.program_build_version()
    }

    pub fn has_taxon_id_lists(&self) -> bool {
        self.header2.taxon_array_offset != 0
    }

    pub fn has_taxon_nodes(&self) -> bool {
        self.header2.taxon_nodes_offset != 0
    }

    pub fn has_taxon_scientific_names(&self) -> bool {
        self.header2.taxon_names_offset != 0
    }

    pub fn taxids(&self, oid: usize) -> Result<Vec<TaxId>, String> {
        self.taxon_list
            .as_ref()
            .ok_or_else(|| "Taxonomy mapping was not loaded.".to_string())?
            .get(oid)
    }

    pub fn rank(&self, taxid: TaxId) -> Result<i32, String> {
        Ok(self
            .taxon_nodes
            .as_ref()
            .ok_or_else(|| "Taxonomy nodes were not loaded.".to_string())?
            .rank(taxid))
    }

    pub fn get_parent(&self, taxid: TaxId) -> Result<TaxId, String> {
        self.taxon_nodes
            .as_ref()
            .ok_or_else(|| "Taxonomy nodes were not loaded.".to_string())?
            .get_parent(taxid)
    }

    pub fn max_taxid(&self) -> Result<TaxId, String> {
        Ok(self
            .taxon_nodes
            .as_ref()
            .ok_or_else(|| "Taxonomy nodes were not loaded.".to_string())?
            .max())
    }

    pub fn taxon_scientific_name(&self, taxid: TaxId) -> String {
        if taxid >= 0 {
            if let Some(name) = self.taxon_scientific_names.get(taxid as usize) {
                if !name.is_empty() {
                    return name.clone();
                }
            }
        }
        taxid.to_string()
    }

    pub fn seq_length(&mut self, oid: usize) -> Result<Loc, String> {
        if let Some(&length) = self.seq_length.get(oid) {
            return Ok(length);
        }
        Ok(self.seq_info_at(oid as OId)?.seq_len as Loc)
    }

    pub fn read_seqid_list(&mut self) -> Result<(), String> {
        self.seq_length.clear();
        self.seqids.clear();
        self.seq_length.reserve(self.sequence_count() as usize);
        self.seqids.reserve(self.sequence_count() as usize);
        for oid in 0..self.sequence_count() {
            let (sequence, mut id) = self.record(oid)?;
            crate::util::sequence::fix_title(&mut id);
            self.seq_length.push(sequence.len() as Loc);
            self.seqids.push(id);
        }
        Ok(())
    }

    pub fn seqids(&self) -> &[String] {
        &self.seqids
    }

    pub fn end_random_access(&mut self, dictionary: bool) {
        if dictionary {
            // Sequence IDs are streamed on demand in this implementation;
            // lengths are the only eager random-access dictionary we retain.
            self.seq_length.clear();
            self.seq_length.shrink_to_fit();
            self.seqids.clear();
            self.seqids.shrink_to_fit();
        }
    }

    pub fn seq_data(&mut self, _oid: usize, _dst: &mut Vec<Letter>) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn filter_by_accession(&self, _path: &str) -> Result<(), String> {
        Err("The .dmnd database format does not support filtering by accession.".to_string())
    }

    pub fn init_write(&self) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn write_seq(&self, _seq: &[Letter], _id: &str) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn file_name(&self) -> &Path {
        &self.path
    }

    pub fn file_count(&self) -> i64 {
        1
    }

    pub fn close(&mut self) {
        // `File` has no explicit close operation; dropping DatabaseFile closes
        // the handle. This method exists for parity with SequenceFile callers.
    }

    pub fn data_start(&self) -> u64 {
        self.data_start
    }

    pub fn read_sequence_with_delimiters(&mut self, oid: OId) -> Result<Vec<Letter>, String> {
        let (sequence, _) = self.record(oid)?;
        let mut output = Vec::with_capacity(sequence.len() + 2);
        output.push(DELIMITER_LETTER);
        output.extend(sequence);
        output.push(DELIMITER_LETTER);
        Ok(output)
    }
}

fn resolve_database_path(path: &Path) -> PathBuf {
    if path.exists() {
        return path.to_path_buf();
    }
    let appended = PathBuf::from(format!(
        "{}{}",
        path.display(),
        DatabaseFile::FILE_EXTENSION
    ));
    if appended.exists() {
        appended
    } else {
        path.to_path_buf()
    }
}

fn metadata_end(offset: u64, header2: &ReferenceHeader2, file_len: u64) -> u64 {
    [
        header2.taxon_array_offset,
        header2.taxon_nodes_offset,
        header2.taxon_names_offset,
        file_len,
    ]
    .into_iter()
    .filter(|&candidate| candidate > offset)
    .min()
    .unwrap_or(file_len)
}

fn read_string_vector<R: Read>(reader: &mut R) -> Result<Vec<String>, String> {
    let mut count = [0; 4];
    reader.read_exact(&mut count).map_err(|e| e.to_string())?;
    let count = u32::from_le_bytes(count) as usize;
    let mut strings = Vec::with_capacity(count);
    for _ in 0..count {
        let mut bytes = Vec::new();
        loop {
            let mut byte = [0];
            reader.read_exact(&mut byte).map_err(|e| e.to_string())?;
            if byte[0] == 0 {
                break;
            }
            bytes.push(byte[0]);
        }
        strings.push(String::from_utf8_lossy(&bytes).into_owned());
    }
    Ok(strings)
}

pub struct SeqLengthIterator {
    file: File,
    pos_array_offset: u64,
    oid: OId,
    total: OId,
    current_len: u32,
}

impl SeqLengthIterator {
    pub fn new(
        path: impl AsRef<Path>,
        pos_array_offset: u64,
        start: OId,
        total: OId,
    ) -> Result<Self, String> {
        let mut iterator = Self {
            file: File::open(path).map_err(|e| e.to_string())?,
            pos_array_offset,
            oid: start,
            total,
            current_len: 0,
        };
        if iterator.valid() {
            iterator.load_current()?;
        }
        Ok(iterator)
    }

    fn load_current(&mut self) -> Result<(), String> {
        self.file
            .seek(SeekFrom::Start(
                self.pos_array_offset + self.oid * SeqInfo::SIZE as u64,
            ))
            .map_err(|e| e.to_string())?;
        self.current_len = SeqInfo::read_from(&mut self.file)
            .map_err(|e| e.to_string())?
            .seq_len;
        Ok(())
    }

    pub fn valid(&self) -> bool {
        self.oid < self.total
    }

    pub fn oid(&self) -> OId {
        self.oid
    }

    pub fn current_len(&self) -> Option<Loc> {
        self.valid().then_some(self.current_len as Loc)
    }

    pub fn advance(&mut self) -> Result<(), String> {
        self.oid += 1;
        if self.valid() {
            self.load_current()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn fixture_path(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-dmnd-{label}-{}-{}.dmnd",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    fn write_fixture(path: &Path) {
        let mut bytes = Vec::new();
        let mut header = ReferenceHeader::new();
        header.sequences = 2;
        header.letters = 5;
        header.pos_array_offset = 110;
        header.write_to(&mut bytes).unwrap();
        ReferenceHeader2::new().write_record_to(&mut bytes).unwrap();
        assert_eq!(bytes.len(), 96);

        bytes.extend_from_slice(&[0xff, 0, 1, 2, 0xff, b'a', 0]);
        bytes.extend_from_slice(&[0xff, 3, 4, 0xff, b'b', b'b', 0]);
        SeqInfo::new(96, 3).write_to(&mut bytes).unwrap();
        SeqInfo::new(103, 2).write_to(&mut bytes).unwrap();
        SeqInfo::new(110, 0).write_to(&mut bytes).unwrap();
        std::fs::write(path, bytes).unwrap();
    }

    #[test]
    fn test_header_roundtrip() {
        let mut h = ReferenceHeader::new();
        h.sequences = 42;
        h.letters = 12345;
        h.pos_array_offset = 999;

        let mut buf = Vec::new();
        h.write_to(&mut buf).unwrap();
        assert_eq!(buf.len(), 40);

        let h2 = ReferenceHeader::read_from(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(h2.magic_number, MAGIC_NUMBER);
        assert_eq!(h2.sequences, 42);
        assert_eq!(h2.letters, 12345);
        assert_eq!(h2.pos_array_offset, 999);
    }

    #[test]
    fn test_header2_dynamic_record_and_seqinfo_binary_layout() {
        let header = ReferenceHeader2 {
            hash: [7; 16],
            taxon_array_offset: 11,
            taxon_array_size: 12,
            taxon_nodes_offset: 13,
            taxon_names_offset: 14,
        };
        let mut bytes = Vec::new();
        header.write_record_to(&mut bytes).unwrap();
        assert_eq!(bytes.len(), 8 + ReferenceHeader2::PAYLOAD_SIZE);
        assert_eq!(
            &bytes[..8],
            &(ReferenceHeader2::PAYLOAD_SIZE as u64).to_le_bytes()
        );
        assert_eq!(
            ReferenceHeader2::read_record_from(&mut Cursor::new(bytes)).unwrap(),
            header
        );

        let mut seq_info = Vec::new();
        SeqInfo {
            pos: 99,
            seq_len: 17,
            padding: u32::MAX,
        }
        .write_to(&mut seq_info)
        .unwrap();
        assert_eq!(seq_info.len(), SeqInfo::SIZE);
        assert_eq!(&seq_info[12..], &[0, 0, 0, 0]);
    }

    #[test]
    fn test_header_validation() {
        // `ReferenceHeader::new()` builds an INCOMPLETE header (sequences=0)
        // that the writer fills in. Validating it should fail, matching C++
        // which rejects `sequences == 0` as the incomplete-database sentinel.
        let h = ReferenceHeader::new();
        assert!(h.validate().is_err());

        // A complete header (sequences > 0) passes.
        let mut complete = ReferenceHeader::new();
        complete.sequences = 1;
        assert!(complete.validate().is_ok());

        // Bad magic fails.
        let mut bad = ReferenceHeader::new();
        bad.magic_number = 0;
        bad.sequences = 1;
        assert!(bad.validate().is_err());

        // Future db_version fails.
        let mut newer = ReferenceHeader::new();
        newer.sequences = 1;
        newer.db_version = 999;
        assert!(newer.validate().is_err());

        // Old build fails.
        let mut old = ReferenceHeader::new();
        old.sequences = 1;
        old.build = 1;
        assert!(old.validate().is_err());
    }

    #[test]
    fn test_is_diamond_db() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&MAGIC_NUMBER.to_le_bytes());
        assert!(is_diamond_db(&mut Cursor::new(&buf)));

        let bad_buf = vec![0u8; 8];
        assert!(!is_diamond_db(&mut Cursor::new(&bad_buf)));
    }

    #[test]
    fn test_read_real_db() {
        let db_path = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/data.dmnd");
        if let Ok(mut f) = std::fs::File::open(db_path) {
            let h = ReferenceHeader::read_from(&mut f).unwrap();
            assert_eq!(h.magic_number, MAGIC_NUMBER);
            assert!(h.validate().is_ok());
            assert!(h.sequences > 0);
        }
    }

    #[test]
    fn test_database_file_sequence_access_length_iterator_and_partition() {
        let path = fixture_path("access");
        write_fixture(&path);
        let mut db = DatabaseFile::open(&path).unwrap();
        assert_eq!(db.file_count(), 1);
        assert_eq!(db.sequence_count(), 2);
        assert_eq!(db.letters(), 5);
        assert_eq!(db.data_start(), 96);
        assert_eq!(db.read_seq().unwrap(), Some((vec![0, 1, 2], "a".into())));
        db.skip_seq().unwrap();
        assert_eq!(db.read_seq().unwrap(), None);

        db.set_seqinfo_ptr(1).unwrap();
        assert_eq!(db.tell_seq(), 1);
        assert_eq!(db.read_seqinfo().unwrap(), SeqInfo::new(103, 2));
        assert!(db.eof());
        db.putback_seqinfo().unwrap();
        assert_eq!(db.tell_seq(), 1);
        assert_eq!(db.seq_length(0).unwrap(), 3);
        assert_eq!(
            DatabaseFile::id_len(&SeqInfo::new(96, 3), &SeqInfo::new(103, 2)),
            Ok(1)
        );
        assert_eq!(
            db.read_sequence_with_delimiters(1).unwrap(),
            vec![31, 3, 4, 31]
        );
        db.read_seqid_list().unwrap();
        assert_eq!(db.seqids(), &["a".to_string(), "bb".to_string()]);
        assert_eq!(db.seq_length(1).unwrap(), 2);

        let mut lengths = db.seq_length_iterator(0).unwrap();
        assert_eq!(lengths.oid(), 0);
        assert_eq!(lengths.current_len(), Some(3));
        lengths.advance().unwrap();
        assert_eq!(lengths.current_len(), Some(2));
        lengths.advance().unwrap();
        assert!(!lengths.valid());
        assert_eq!(lengths.current_len(), None);

        db.create_partition(2).unwrap();
        assert_eq!(
            db.partition.chunks,
            vec![Chunk::new(1, 1, 1), Chunk::new(0, 0, 1)]
        );
        assert_eq!(db.get_n_partition_chunks(), 2);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_chunk_text_partition_roundtrip_and_extension_resolution() {
        assert_eq!(to_chunk("4 12 9 annotation").unwrap(), Chunk::new(4, 12, 9));
        assert_eq!(chunk_to_string(&Chunk::new(4, 12, 9)), "4 12 9");

        let path = fixture_path("partition");
        write_fixture(&path);
        let extensionless = path.with_extension("");
        let mut db = DatabaseFile::open(&extensionless).unwrap();
        db.create_partition_fixed_number(2).unwrap();
        let partition_path = fixture_path("chunks").with_extension("txt");
        db.save_partition(&partition_path, "sample").unwrap();
        let expected = db.partition.chunks.clone();
        db.clear_partition();
        db.load_partition(&partition_path).unwrap();
        assert_eq!(db.partition.chunks, expected);
        std::fs::remove_file(path).unwrap();
        std::fs::remove_file(partition_path).unwrap();
    }
}
