//! Native DIAMOND database implementation mirroring
//! `diamond/src/data/dmnd/dmnd.{h,cpp}`.

use crate::basic::value::{Letter, Loc, OId, SequenceType, TaxId, DELIMITER_LETTER};
use crate::data::taxon_list::TaxonList;
use crate::data::taxonomy_nodes::TaxonomyNodes;
use crate::util::io::{Deserializer, VecStream};
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom, Write};
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
        // Source DatabaseFile::init tests the old build/version condition
        // first, including when a header also advertises a newer version.
        if self.build < MIN_BUILD_VERSION || self.db_version < MIN_DB_VERSION {
            return Err(
                "Database was built with an older version of Diamond and is incompatible.".into(),
            );
        }
        let max_db_version = CURRENT_DB_VERSION_PROT.max(CURRENT_DB_VERSION_NUCL);
        if self.db_version > max_db_version {
            return Err(
                "Database was built with a newer version of Diamond and is incompatible.".into(),
            );
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
    /// Source EXTRA's int32 SequenceType tail; absent baseline records use 0.
    pub db_type: i32,
}

// Source prefixes sizeof(ReferenceHeader2), including C ABI tail alignment,
// while serializing only the actual fields (no padding bytes).
#[repr(C)]
struct Header2ExtraLayout {
    hash: [u8; 16],
    offsets: [u64; 4],
    db_type: i32,
}

impl ReferenceHeader2 {
    pub const PAYLOAD_SIZE: usize = if cfg!(feature = "extra") { 52 } else { 48 };
    pub const RECORD_SIZE: usize = if cfg!(feature = "extra") {
        std::mem::size_of::<Header2ExtraLayout>()
    } else {
        48
    };
    pub fn new() -> Self {
        Self::default()
    }

    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut hash = [0u8; 16];
        reader.read_exact(&mut hash)?;
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf)?;
        let mut db_type = 0;
        if cfg!(feature = "extra") {
            let mut tail = [0; 4];
            reader.read_exact(&mut tail)?;
            db_type = i32::from_le_bytes(tail);
        }
        Ok(Self {
            hash,
            taxon_array_offset: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            taxon_array_size: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
            taxon_nodes_offset: u64::from_le_bytes(buf[16..24].try_into().unwrap()),
            taxon_names_offset: u64::from_le_bytes(buf[24..32].try_into().unwrap()),
            db_type,
        })
    }
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.hash)?;
        writer.write_all(&self.taxon_array_offset.to_le_bytes())?;
        writer.write_all(&self.taxon_array_size.to_le_bytes())?;
        writer.write_all(&self.taxon_nodes_offset.to_le_bytes())?;
        writer.write_all(&self.taxon_names_offset.to_le_bytes())?;
        if cfg!(feature = "extra") {
            writer.write_all(&self.db_type.to_le_bytes())?;
        }
        Ok(())
    }
    pub fn write_record_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&(Self::RECORD_SIZE as u64).to_le_bytes())?;
        self.write_to(writer)
    }
    /// Source DynamicRecordReader fills a field with zero when it does not
    /// fit, without consuming that field's bytes; Finish consumes the remainder.
    pub fn read_record_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        fn field<R: Read, const N: usize>(reader: &mut R, left: &mut u64) -> io::Result<[u8; N]> {
            let mut bytes = [0; N];
            if *left >= N as u64 {
                reader.read_exact(&mut bytes)?;
                *left -= N as u64;
            }
            Ok(bytes)
        }
        let mut bytes = [0; 8];
        reader.read_exact(&mut bytes)?;
        let mut left = u64::from_le_bytes(bytes);
        let hash = field(reader, &mut left)?;
        let taxon_array_offset = u64::from_le_bytes(field(reader, &mut left)?);
        let taxon_array_size = u64::from_le_bytes(field(reader, &mut left)?);
        let taxon_nodes_offset = u64::from_le_bytes(field(reader, &mut left)?);
        let taxon_names_offset = u64::from_le_bytes(field(reader, &mut left)?);
        let db_type = if cfg!(feature = "extra") {
            i32::from_le_bytes(field(reader, &mut left)?)
        } else {
            0
        };
        let skipped = io::copy(&mut (&mut *reader).take(left), &mut io::sink())?;
        if skipped != left {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Unexpected end of file.",
            ));
        }
        Ok(Self {
            hash,
            taxon_array_offset,
            taxon_array_size,
            taxon_nodes_offset,
            taxon_names_offset,
            db_type,
        })
    }
    pub fn sequence_type(&self) -> Result<crate::basic::value::SequenceType, String> {
        match self.db_type {
            0 => Ok(crate::basic::value::SequenceType::AminoAcid),
            1 => Ok(crate::basic::value::SequenceType::Nucleotide),
            value => Err(format!("Invalid database sequence type: {value}")),
        }
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

/// Source `std::stoull(token)` accepts a decimal prefix, optional sign, and
/// C whitespace. It checks the unsigned magnitude before applying the sign.
fn chunk_unsigned(token: &[u8]) -> io::Result<u64> {
    let mut c_token = token.to_vec();
    c_token.push(0);
    let start = c_token.as_ptr().cast::<libc::c_char>();
    let mut end = std::ptr::null_mut();
    let value = unsafe { libc::strtoull(start, &mut end, 10) };
    if end == start.cast_mut() {
        return Err(io::Error::other("stoull"));
    }
    let consumed = unsafe { end.offset_from(start) } as usize;
    // This is the ERANGE condition for base ten, independent of the CRT's
    // thread-local errno accessor. Do not mistake a valid ULLONG_MAX for overflow.
    let mut magnitude = 0u64;
    for &byte in &c_token[..consumed] {
        if byte.is_ascii_digit() {
            magnitude = magnitude
                .checked_mul(10)
                .and_then(|n| n.checked_add(u64::from(byte - b'0')))
                .ok_or_else(|| io::Error::other("stoull"))?;
        }
    }
    Ok(value as u64)
}

/// Canonical raw record parser used by partition files and MPP work records.
/// Source splits on literal spaces, retaining empty fields (not Unicode whitespace).
pub fn to_chunk_bytes(line: &[u8]) -> io::Result<Chunk> {
    let mut fields = line.split(|&byte| byte == b' ');
    let mut next = || {
        fields.next().ok_or_else(|| {
            // Source indexes a missing vector element here; keep the safe native
            // boundary explicit rather than manufacturing a chunk from invalid reads.
            io::Error::other("Missing chunk field")
        })
    };
    let i = chunk_unsigned(next()?)? as i32;
    let offset = chunk_unsigned(next()?)? as usize;
    let n_seqs = chunk_unsigned(next()?)? as usize as i64;
    Ok(Chunk::new(i, offset, n_seqs))
}

pub fn to_chunk(line: &str) -> Result<Chunk, String> {
    to_chunk_bytes(line.as_bytes()).map_err(|error| error.to_string())
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
    /// Explicit source input-buffer capacity; None retains the library default.
    pub file_buffer_size: Option<usize>,
}

#[derive(Debug)]
enum DatabaseInput {
    Raw(File),
    Source {
        input: crate::util::io::InputFile,
        consumed: u64,
    },
}
impl Read for DatabaseInput {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Raw(file) => file.read(bytes),
            Self::Source { input, consumed } => {
                let count = input.read(bytes)?;
                *consumed += count as u64;
                Ok(count)
            }
        }
    }
}
impl Seek for DatabaseInput {
    fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
        match self {
            Self::Raw(file) => file.seek(from),
            Self::Source { input, consumed } => {
                // This position is private Rust buffering bookkeeping, never
                // Source InputFile::tell (which is invalid on decoded streams).
                if from == SeekFrom::Current(0) {
                    return Ok(*consumed);
                }
                let position = match from {
                    SeekFrom::Start(position) => position,
                    SeekFrom::Current(offset) => consumed
                        .checked_add_signed(offset)
                        .ok_or_else(|| io::Error::other("Invalid input seek."))?,
                    SeekFrom::End(offset) => (input.file_size().map_err(io::Error::other)? as u64)
                        .checked_add_signed(offset)
                        .ok_or_else(|| io::Error::other("Invalid input seek."))?,
                };
                input.seek(position as i64).map_err(io::Error::other)?;
                *consumed = position;
                Ok(position)
            }
        }
    }
}
impl DatabaseInput {
    fn source_profile(&self) -> bool {
        matches!(self, Self::Source { .. })
    }
    fn seekable(&self) -> bool {
        match self {
            Self::Raw(_) => true,
            Self::Source { input, .. } => input.seekable(),
        }
    }
    fn file_len(&mut self) -> io::Result<u64> {
        match self {
            Self::Raw(file) => Ok(file.metadata()?.len()),
            Self::Source { input, .. } => Ok(input.file_size().map_err(io::Error::other)? as u64),
        }
    }
}
trait SourceSeekPolicy {
    fn requires_absolute_seek(&self) -> bool;
}
impl SourceSeekPolicy for File {
    fn requires_absolute_seek(&self) -> bool {
        false
    }
}
impl SourceSeekPolicy for DatabaseInput {
    fn requires_absolute_seek(&self) -> bool {
        !self.seekable()
    }
}
fn seek_buffered<R: Read + Seek + SourceSeekPolicy>(
    reader: &mut BufReader<R>,
    position: u64,
) -> io::Result<()> {
    if reader.get_ref().requires_absolute_seek() {
        // Never let BufReader satisfy a Source seek from prefetched decoded
        // bytes: the actual codec owner must reject every unsupported seek.
        return reader.seek(SeekFrom::Start(position)).map(|_| ());
    }
    let current = reader.stream_position()?;
    if let Ok(offset) = i64::try_from(position as i128 - current as i128) {
        reader.seek_relative(offset)
    } else {
        reader.seek(SeekFrom::Start(position)).map(|_| ())
    }
}

/// Concrete Source command owner over InputFile's decoded stream. Sequential
/// commands preserve the constructor cursor; indexed loads perform actual
/// InputFile seeks, which compressed streams reject rather than emulate.
#[derive(Debug)]
pub(crate) struct SourceDatabaseFile {
    input: crate::util::io::InputFile,
    pub ref_header: ReferenceHeader,
}

impl SourceDatabaseFile {
    pub(crate) fn open(path: &Path, options: DatabaseOpenOptions) -> io::Result<Self> {
        let mut input = crate::util::io::InputFile::with_source_profile_path(
            path,
            crate::util::io::InputFile::BUFFERED,
            options.file_buffer_size.unwrap_or(67_108_864),
        )
        .map_err(io::Error::other)?;
        let ref_header =
            ReferenceHeader::read_from(&mut input).map_err(source_stream_read_error)?;
        if ref_header.magic_number != MAGIC_NUMBER {
            return Err(io::Error::other("Database file is not a DIAMOND database."));
        }
        // Source init returns immediately for NO_COMPATIBILITY_CHECK.
        if !options.no_compatibility_check {
            ref_header.validate().map_err(io::Error::other)?;
            ReferenceHeader2::read_record_from(&mut input).map_err(source_stream_read_error)?;
        }
        Ok(Self { input, ref_header })
    }

    pub(crate) fn read_seq(&mut self) -> io::Result<(Vec<Letter>, String)> {
        let (sequence, title) = self.read_seq_bytes()?;
        Ok((
            sequence,
            String::from_utf8(title).map_err(io::Error::other)?,
        ))
    }
    pub(crate) fn read_seq_bytes(&mut self) -> io::Result<(Vec<Letter>, Vec<u8>)> {
        let mut marker = [0];
        let _ = self.input.read(&mut marker)?;
        let mut sequence = Vec::new();
        self.input
            .read_to(&mut sequence, 0xff)
            .map_err(io::Error::other)?;
        let mut title = Vec::new();
        self.input
            .read_to(&mut title, 0)
            .map_err(io::Error::other)?;
        Ok((
            sequence
                .into_iter()
                .map(|letter| letter as Letter)
                .collect(),
            title,
        ))
    }

    /// Source load_twopass with SEQS and without TITLES, as makeindex uses it.
    /// Retain only the position metadata and one decoded sequence at a time;
    /// the caller owns the actual complete indexing Block.
    pub(crate) fn visit_sequences(
        &mut self,
        max_letters: u64,
        mut visit: impl FnMut(u64, &[Letter]) -> io::Result<()>,
    ) -> io::Result<()> {
        self.input
            .seek(self.ref_header.pos_array_offset as i64)
            .map_err(io::Error::other)?;
        let mut current = SeqInfo::read_from(&mut self.input).map_err(source_stream_read_error)?;
        let mut positions = Vec::new();
        let mut letters = 0u64;
        while current.seq_len > 0 && letters < max_letters {
            let next = SeqInfo::read_from(&mut self.input).map_err(source_stream_read_error)?;
            letters += u64::from(current.seq_len);
            positions.push(current);
            current = next;
        }
        let Some(first) = positions.first() else {
            return Ok(());
        };
        self.input
            .seek(first.pos as i64)
            .map_err(io::Error::other)?;
        let mut sequence = Vec::new();
        for (oid, info) in positions.into_iter().enumerate() {
            sequence.resize(info.seq_len as usize + 2, 0u8);
            self.input
                .read_exact(&mut sequence)
                .map_err(source_stream_read_error)?;
            // Source read_seq_data reads delimiter-framed bytes but replaces
            // the perimeter instead of validating the stored marker values.
            let letters = unsafe {
                std::slice::from_raw_parts(
                    sequence[1..sequence.len() - 1].as_ptr().cast::<Letter>(),
                    info.seq_len as usize,
                )
            };
            visit(oid as u64, letters)?;
            if !self.input.seek_forward_delim(0).map_err(io::Error::other)? {
                return Err(io::Error::other("Unexpected end of file."));
            }
        }
        Ok(())
    }
}

fn source_stream_read_error(error: io::Error) -> io::Error {
    if error.kind() == io::ErrorKind::UnexpectedEof {
        io::Error::new(io::ErrorKind::UnexpectedEof, "Unexpected end of input.")
    } else {
        error
    }
}

/// Read-side implementation of C++ `DatabaseFile`.
#[derive(Debug)]
pub struct DatabaseFile {
    path: PathBuf,
    file: BufReader<DatabaseInput>,
    index_file: Option<BufReader<File>>,
    pub ref_header: ReferenceHeader,
    pub header2: ReferenceHeader2,
    pub partition: DatabasePartition,
    pos_array_cursor: OId,
    sequence_cursor: OId,
    data_start: u64,
    seq_length: Vec<Loc>,
    seqids: Vec<String>,
    pub(crate) accession_spool:
        Option<std::sync::Arc<crate::data::accession_spool::AccessionSpool>>,
    taxon_list: Option<TaxonList>,
    taxon_nodes: Option<TaxonomyNodes>,
    taxon_scientific_names: Vec<Vec<u8>>,
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

    pub fn make_db_paths(
        input_files: &[&Path],
        output_path: &Path,
        sequence_type: SequenceType,
    ) -> io::Result<crate::data::db_builder::BuildStats> {
        crate::data::db_builder::build_db_paths(input_files, output_path, sequence_type)
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
        Self::open_with_options_io(path, options).map_err(|error| error.to_string())
    }

    /// Source command constructor preserving actual file-open failures.
    pub fn open_with_options_io(
        path: impl AsRef<Path>,
        options: DatabaseOpenOptions,
    ) -> io::Result<Self> {
        Self::open_profile(path.as_ref(), options, false, false)
    }

    /// Actual Source command owner, retaining the decoder for metadata and
    /// record operations. Indexed library constructors retain their raw profile.
    pub(crate) fn open_source_with_options_io(
        path: &Path,
        options: DatabaseOpenOptions,
        require_ranks: bool,
    ) -> io::Result<Self> {
        Self::open_profile(path, options, true, require_ranks)
    }

    fn open_profile(
        path: &Path,
        options: DatabaseOpenOptions,
        source_profile: bool,
        require_ranks: bool,
    ) -> io::Result<Self> {
        let path = resolve_database_path(path);
        let input = if source_profile {
            DatabaseInput::Source {
                input: crate::util::io::InputFile::with_source_profile_path(
                    &path,
                    crate::util::io::InputFile::BUFFERED,
                    options.file_buffer_size.unwrap_or(67_108_864),
                )
                .map_err(io::Error::other)?,
                consumed: 0,
            }
        } else {
            let file = File::open(&path)
                .map_err(|error| crate::util::io::FileOpenError::into_io(&path, error))?;
            if let Some(size) = options.file_buffer_size {
                if size < 4 && file.metadata()?.len() >= size as u64 {
                    return Err(io::Error::other("Invalid peek"));
                }
            }
            DatabaseInput::Raw(file)
        };
        // InputFile already owns Source's configured buffer. This small outer
        // window provides the historical BufRead API without duplicating it.
        let mut file = BufReader::with_capacity(
            if source_profile {
                8192
            } else {
                options.file_buffer_size.unwrap_or(8192).max(1)
            },
            input,
        );
        let source_read_error = |error: io::Error| {
            if error.kind() == io::ErrorKind::UnexpectedEof {
                io::Error::new(io::ErrorKind::UnexpectedEof, "Unexpected end of input.")
            } else {
                error
            }
        };
        let ref_header = ReferenceHeader::read_from(&mut file).map_err(source_read_error)?;
        if ref_header.magic_number != MAGIC_NUMBER {
            return Err(io::Error::other("Database file is not a DIAMOND database."));
        }
        if !options.no_compatibility_check {
            ref_header.validate().map_err(io::Error::other)?;
        }
        // Source DatabaseFile::init returns immediately after Header1 when
        // NO_COMPATIBILITY_CHECK is requested (dbinfo). In particular this
        // inspection path does not require a complete Header2 record.
        let header2 = if options.no_compatibility_check {
            ReferenceHeader2::default()
        } else {
            ReferenceHeader2::read_record_from(&mut file).map_err(source_read_error)?
        };
        let data_start = file.stream_position()?;
        // Keep index and sequence readahead independent: alternating between
        // their distant offsets must not refill a large buffer per record.
        let index_file = if !file.get_ref().seekable() {
            None
        } else {
            Some(BufReader::with_capacity(
                options.file_buffer_size.unwrap_or(8192).clamp(1, 1 << 20),
                File::open(&path)
                    .map_err(|error| crate::util::io::FileOpenError::into_io(&path, error))?,
            ))
        };
        let mut db = Self {
            path,
            file,
            index_file,
            ref_header,
            header2,
            partition: DatabasePartition::default(),
            pos_array_cursor: 0,
            sequence_cursor: 0,
            data_start,
            seq_length: Vec::new(),
            seqids: Vec::new(),
            accession_spool: None,
            taxon_list: None,
            taxon_nodes: None,
            taxon_scientific_names: Vec::new(),
        };
        db.validate_requested_metadata(options, require_ranks)
            .map_err(io::Error::other)?;
        db.load_requested_metadata(options)
            .map_err(io::Error::other)?;
        Ok(db)
    }

    fn validate_requested_metadata(
        &self,
        options: DatabaseOpenOptions,
        require_ranks: bool,
    ) -> Result<(), String> {
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
        if require_ranks && self.ref_header.build < 131 {
            missing.push("taxonomy ranks information (database needs to be built with diamond version >= 0.9.30");
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
        if self.file.get_ref().source_profile() && options.load_taxon_scientific_names {
            self.file
                .seek(SeekFrom::Start(self.header2.taxon_names_offset))
                .map_err(|e| e.to_string())?;
            self.taxon_scientific_names = read_string_vector(&mut self.file)?;
        }
        if options.load_taxon_nodes {
            self.file
                .seek(SeekFrom::Start(self.header2.taxon_nodes_offset))
                .map_err(|e| e.to_string())?;
            let file_len = self.file.get_mut().file_len().map_err(|e| e.to_string())?;
            let end = metadata_end(self.header2.taxon_nodes_offset, &self.header2, file_len);
            let mut data = vec![0; (end - self.header2.taxon_nodes_offset) as usize];
            self.file.read_exact(&mut data).map_err(|e| e.to_string())?;
            let mut input = Deserializer::new(VecStream::from_vec(data));
            self.taxon_nodes = Some(
                TaxonomyNodes::from_deserializer(&mut input, self.ref_header.build)
                    .map_err(|e| e.to_string())?,
            );
        }
        if !self.file.get_ref().source_profile() && options.load_taxon_scientific_names {
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
        self.seq_info_at_io(oid).map_err(|error| error.to_string())
    }
    fn seq_info_at_io(&mut self, oid: OId) -> io::Result<SeqInfo> {
        let position = self.ref_header.pos_array_offset + SeqInfo::SIZE as u64 * oid;
        if let Some(index_file) = &mut self.index_file {
            seek_buffered(index_file, position)?;
            SeqInfo::read_from(index_file).map_err(source_stream_read_error)
        } else {
            seek_buffered(&mut self.file, position)?;
            SeqInfo::read_from(&mut self.file).map_err(source_stream_read_error)
        }
    }

    pub fn init_seq_access(&mut self) {
        self.sequence_cursor = 0;
    }

    pub fn init_seqinfo_access(&mut self) {
        // The Rust reader seeks directly from the logical OID cursor when an
        // entry is read, so no independent buffered-file seek is required.
    }

    pub fn read_seq(&mut self) -> Result<Option<(Vec<Letter>, String)>, String> {
        self.read_seq_bytes()?
            .map(|(sequence, title)| {
                Ok((
                    sequence,
                    String::from_utf8(title).map_err(|error| error.to_string())?,
                ))
            })
            .transpose()
    }

    /// Indexed OID cursor service, independent of Source's physical stream.
    pub fn read_seq_bytes(&mut self) -> Result<Option<(Vec<Letter>, Vec<u8>)>, String> {
        if self.sequence_cursor >= self.ref_header.sequences {
            return Ok(None);
        }
        let record = self
            .record_bytes_io(self.sequence_cursor)
            .map_err(|error| error.to_string())?;
        self.sequence_cursor += 1;
        Ok(Some(record))
    }

    /// Source init_seq_access seeks its fixed C ABI header layout, which can
    /// differ from both the constructor cursor and indexed record positions.
    /// getseq deliberately does not call this before its direct first read.
    pub(crate) fn init_source_sequential(&mut self) -> io::Result<()> {
        let source_start = ReferenceHeader::SIZE as u64 + 8 + ReferenceHeader2::RECORD_SIZE as u64;
        seek_buffered(&mut self.file, source_start)
    }

    /// Source `DatabaseFile::read_seq` on the existing physical descriptor.
    /// In particular, getseq starts at the cursor left by Header2 extraction,
    /// including DynamicRecordReader::Finish's advertised-size consumption.
    /// Indexed `record`/`read_seq` callers retain their independent OID service.
    pub(crate) fn read_source_sequential(&mut self) -> io::Result<(Vec<Letter>, String)> {
        let (sequence, title) = self.read_source_sequential_bytes()?;
        Ok((
            sequence,
            String::from_utf8(title).map_err(io::Error::other)?,
        ))
    }

    /// Source std::string can contain non-UTF8 bytes, including bytes reached
    /// by its fixed physical header cursor. MCL writes these directly.
    pub(crate) fn read_source_sequential_bytes(&mut self) -> io::Result<(Vec<Letter>, Vec<u8>)> {
        use std::io::BufRead;
        let mut marker = [0u8; 1];
        // Source ignores both the leading byte value and its short-read count.
        let _ = self.file.read(&mut marker)?;
        let mut sequence = Vec::new();
        self.file.read_until(0xff, &mut sequence)?;
        if sequence.last() == Some(&0xff) {
            sequence.pop();
        }
        let mut title = Vec::new();
        self.file.read_until(0, &mut title)?;
        if title.last() == Some(&0) {
            title.pop();
        }
        Ok((
            sequence.into_iter().map(|byte| byte as Letter).collect(),
            title,
        ))
    }

    pub fn skip_seq(&mut self) -> Result<(), String> {
        if self.sequence_cursor >= self.ref_header.sequences {
            return Err("Unexpected end of file.".to_string());
        }
        self.record_bytes_io(self.sequence_cursor)
            .map_err(|error| error.to_string())?;
        self.sequence_cursor += 1;
        Ok(())
    }

    /// Source read_seqid_list, without retaining decoded sequence bytes or IDs.
    /// Scan the actual delimiter framing and preserve title repair diagnostics.
    pub(crate) fn write_accession_spool(
        &mut self,
        spool: &crate::data::accession_spool::AccessionSpool,
    ) -> io::Result<()> {
        use std::io::BufRead;
        fn until(
            reader: &mut impl BufRead,
            delimiter: u8,
            mut text: Option<&mut Vec<u8>>,
        ) -> io::Result<()> {
            loop {
                let bytes = reader.fill_buf()?;
                if bytes.is_empty() {
                    return Ok(());
                }
                let stop = bytes.iter().position(|&byte| byte == delimiter);
                let count = stop.unwrap_or(bytes.len());
                if let Some(out) = text.as_mut() {
                    out.extend_from_slice(&bytes[..count]);
                }
                reader.consume(count + usize::from(stop.is_some()));
                if stop.is_some() {
                    return Ok(());
                }
            }
        }
        // Source read_seqid_list calls init_seq_access, which seeks a fixed
        // C ABI layout offset, independently of DynamicRecordReader's consumed
        // byte count and independently of the position-array record API.
        let source_start = ReferenceHeader::SIZE as u64 + 8 + ReferenceHeader2::RECORD_SIZE as u64;
        seek_buffered(&mut self.file, source_start).map_err(io::Error::other)?;
        for _ in 0..self.sequence_count() {
            let mut marker = [0u8; 1];
            // Source read_seq ignores the byte's value and its short-read count.
            let _ = self.file.read(&mut marker)?;
            until(&mut self.file, 0xff, None)?;
            let mut title = Vec::new();
            until(&mut self.file, 0, Some(&mut title))?;
            if let Some(message) = crate::util::sequence::fix_title_bytes(&mut title) {
                crate::util::log_stream::message_stream()
                    .lock()
                    .unwrap()
                    .write(format!("Warning: {message}\n"))?;
            }
            spool.push_bytes(&crate::util::sequence::seqid_bytes(&title))?;
        }
        Ok(())
    }

    /// Read only the title of a physical OID, without allocating its sequence
    /// or changing either logical sequence cursor. Keep framing validation.
    pub fn title(&mut self, oid: OId) -> Result<String, String> {
        self.title_bytes_io(oid)
            .map_err(|error| error.to_string())
            .and_then(|id| String::from_utf8(id).map_err(|error| error.to_string()))
    }
    pub(crate) fn title_bytes_io(&mut self, oid: OId) -> io::Result<Vec<u8>> {
        let info = self.seq_info_at_io(oid)?;
        let end = info
            .pos
            .checked_add(u64::from(info.seq_len))
            .and_then(|pos| pos.checked_add(1))
            .ok_or_else(|| io::Error::other("Invalid sequence offset in DIAMOND database."))?;
        seek_buffered(&mut self.file, info.pos)?;
        let mut first = [0u8; 1];
        self.file.read_exact(&mut first)?;
        seek_buffered(&mut self.file, end)?;
        let mut last = [0u8; 1];
        self.file.read_exact(&mut last)?;
        if first != [0xff] || last != [0xff] {
            return Err(io::Error::other(
                "Invalid sequence framing in DIAMOND database.",
            ));
        }
        let mut id = Vec::new();
        loop {
            let mut byte = [0u8; 1];
            self.file.read_exact(&mut byte).map_err(|error| {
                if error.kind() == io::ErrorKind::UnexpectedEof {
                    io::Error::other("Unexpected end of file.")
                } else {
                    error
                }
            })?;
            if byte[0] == 0 {
                break;
            }
            id.push(byte[0]);
        }
        Ok(id)
    }

    /// Read an OID without changing either logical sequence cursor.
    pub fn record(&mut self, oid: OId) -> Result<(Vec<Letter>, String), String> {
        self.record_io(oid).map_err(|error| error.to_string())
    }
    pub(crate) fn record_io(&mut self, oid: OId) -> io::Result<(Vec<Letter>, String)> {
        let (sequence, id) = self.record_bytes_io(oid)?;
        Ok((sequence, String::from_utf8(id).map_err(io::Error::other)?))
    }
    pub(crate) fn record_bytes_io(&mut self, oid: OId) -> io::Result<(Vec<Letter>, Vec<u8>)> {
        let info = self.seq_info_at_io(oid)?;
        seek_buffered(&mut self.file, info.pos)?;
        let mut framed = vec![0; info.seq_len as usize + 2];
        self.file.read_exact(&mut framed)?;
        if framed.first() != Some(&0xff) || framed.last() != Some(&0xff) {
            return Err(io::Error::other(
                "Invalid sequence framing in DIAMOND database.",
            ));
        }
        let mut id = Vec::new();
        loop {
            let mut byte = [0];
            self.file.read_exact(&mut byte).map_err(|error| {
                if error.kind() == io::ErrorKind::UnexpectedEof {
                    io::Error::other("Unexpected end of file.")
                } else {
                    error
                }
            })?;
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
            id,
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

    /// Source uses an unchecked text ifstream, including silent open/read failure.
    /// Retain its real text-mode owner so Windows newline handling and raw paths
    /// do not pass through a UTF-8 decoder or a second file open.
    pub fn load_partition_io(&mut self, path: &Path) -> io::Result<()> {
        let input = crate::util::io::file::File::open_path(path, "r");
        self.clear_partition();
        let Ok(mut input) = input else {
            return Ok(());
        };
        let mut line = Vec::new();
        let mut byte = [0u8; 1];
        loop {
            let count = input.read_max(&mut byte).unwrap_or(0);
            if count == 0 {
                // getline accepts an unterminated final row only at EOF;
                // a read failure sets badbit and discards a partial row.
                if input.read_at_eof().unwrap_or(false) && !line.is_empty() {
                    self.partition.chunks.push(to_chunk_bytes(&line)?);
                }
                break;
            }
            if byte[0] == b'\n' {
                self.partition.chunks.push(to_chunk_bytes(&line)?);
                line.clear();
            } else {
                line.push(byte[0]);
            }
        }
        Ok(())
    }

    pub fn load_partition(&mut self, path: impl AsRef<Path>) -> Result<(), String> {
        self.load_partition_io(path.as_ref())
            .map_err(|error| error.to_string())
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

    /// Checked String convenience for UTF-8 taxonomy metadata.
    pub fn taxon_scientific_name(&self, taxid: TaxId) -> String {
        self.try_taxon_scientific_name(taxid)
            .expect("Raw taxonomy name is not UTF-8; use taxon_scientific_name_bytes")
    }

    pub fn try_taxon_scientific_name(
        &self,
        taxid: TaxId,
    ) -> Result<String, std::string::FromUtf8Error> {
        String::from_utf8(self.taxon_scientific_name_bytes(taxid))
    }

    pub fn taxon_scientific_name_bytes(&self, taxid: TaxId) -> Vec<u8> {
        if taxid >= 0 {
            if let Some(name) = self.taxon_scientific_names.get(taxid as usize) {
                if !name.is_empty() {
                    return name.clone();
                }
            }
        }
        taxid.to_string().into_bytes()
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
        let (sequence, _) = self
            .record_bytes_io(oid)
            .map_err(|error| error.to_string())?;
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
    let mut appended = path.as_os_str().to_os_string();
    appended.push(DatabaseFile::FILE_EXTENSION);
    let appended = PathBuf::from(appended);
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

fn read_string_vector<R: Read>(reader: &mut R) -> Result<Vec<Vec<u8>>, String> {
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
        strings.push(bytes);
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
        // Deliberately a baseline database under every build profile. This
        // verifies position-array access to legacy files; EXTRA's writer has
        // a different physical field count and advertised C ABI record size.
        bytes.extend_from_slice(&48u64.to_le_bytes());
        bytes.extend_from_slice(&[0u8; 48]);
        assert_eq!(bytes.len(), 96);

        bytes.extend_from_slice(&[0xff, 0, 1, 2, 0xff, b'a', 0]);
        bytes.extend_from_slice(&[0xff, 3, 4, 0xff, b'b', b'b', 0]);
        SeqInfo::new(96, 3).write_to(&mut bytes).unwrap();
        SeqInfo::new(103, 2).write_to(&mut bytes).unwrap();
        SeqInfo::new(110, 0).write_to(&mut bytes).unwrap();
        std::fs::write(path, bytes).unwrap();
    }

    #[test]
    fn inspection_reads_header1_without_requiring_header2() {
        let path = fixture_path("header-only-inspection");
        let header = ReferenceHeader {
            sequences: 3,
            letters: 1280,
            ..ReferenceHeader::default()
        };
        let mut bytes = Vec::new();
        header.write_to(&mut bytes).unwrap();
        std::fs::write(&path, bytes).unwrap();
        let inspection = DatabaseFile::open_with_options_io(
            &path,
            DatabaseOpenOptions {
                no_compatibility_check: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(inspection.sequence_count(), 3);
        assert_eq!(inspection.letters(), 1280);
        drop(inspection);
        let checked = DatabaseFile::open_with_options_io(&path, DatabaseOpenOptions::default())
            .err()
            .unwrap();
        assert_eq!(checked.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(checked.to_string(), "Unexpected end of input.");
        std::fs::remove_file(path).unwrap();
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
            db_type: 0,
        };
        let mut bytes = Vec::new();
        header.write_record_to(&mut bytes).unwrap();
        assert_eq!(bytes.len(), 8 + ReferenceHeader2::PAYLOAD_SIZE);
        assert_eq!(
            &bytes[..8],
            &(ReferenceHeader2::RECORD_SIZE as u64).to_le_bytes()
        );
        // These are following-stream bytes consumed by Source Finish due
        // to the ABI tail alignment advertised by the size prefix. The writer
        // itself emits no such padding.
        let finish_tail = ReferenceHeader2::RECORD_SIZE - ReferenceHeader2::PAYLOAD_SIZE;
        bytes.extend(std::iter::repeat(0xA5).take(finish_tail));
        let mut input = Cursor::new(bytes);
        assert_eq!(
            ReferenceHeader2::read_record_from(&mut input).unwrap(),
            header
        );
        assert_eq!(
            input.position(),
            (8 + ReferenceHeader2::PAYLOAD_SIZE + finish_tail) as u64
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
    fn positional_readers_use_seqinfo_offsets_despite_header2_abi_padding() {
        let path = fixture_path("positional-header2");
        let start = ReferenceHeader::SIZE as u64 + 8 + ReferenceHeader2::PAYLOAD_SIZE as u64;
        let mut header = ReferenceHeader::new();
        header.sequences = 2;
        header.letters = 5;
        header.pos_array_offset = start + 14;
        let mut bytes = Vec::new();
        header.write_to(&mut bytes).unwrap();
        ReferenceHeader2::new().write_record_to(&mut bytes).unwrap();
        assert_eq!(bytes.len() as u64, start);
        bytes.extend_from_slice(&[0xff, 0, 1, 2, 0xff, b'a', 0]);
        bytes.extend_from_slice(&[0xff, 3, 4, 0xff, b'b', b'b', 0]);
        SeqInfo::new(start, 3).write_to(&mut bytes).unwrap();
        SeqInfo::new(start + 7, 2).write_to(&mut bytes).unwrap();
        SeqInfo::new(start + 14, 0).write_to(&mut bytes).unwrap();
        std::fs::write(&path, bytes).unwrap();
        let (_, records) = crate::data::dmnd_reader::read_dmnd(&path).unwrap();
        assert_eq!(
            records
                .iter()
                .map(|r| (r.id.as_str(), r.sequence.as_slice()))
                .collect::<Vec<_>>(),
            [("a", &[0i8, 1, 2][..]), ("bb", &[3i8, 4][..])]
        );
        let mut streamed = Vec::new();
        crate::data::dmnd_reader::for_each_dmnd_record(&path, |oid, id, seq| {
            streamed.push((oid, id.to_owned(), seq.to_vec()));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            streamed,
            [
                (0, "a".to_owned(), vec![0, 1, 2]),
                (1, "bb".to_owned(), vec![3, 4])
            ]
        );
        let mut db = DatabaseFile::open(&path).unwrap();
        assert_eq!(db.record(0).unwrap(), (vec![0, 1, 2], "a".to_owned()));
        std::fs::remove_file(path).unwrap();
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
    fn partition_chunks_follow_source_unsigned_prefix_and_space_grammar() {
        assert_eq!(
            to_chunk_bytes(b"4294967295x 12\xff 9\0ignored raw\xff annotation").unwrap(),
            Chunk::new(-1, 12, 9)
        );
        assert_eq!(to_chunk_bytes(b"4\t12 9 7").unwrap(), Chunk::new(4, 9, 7));
        assert_eq!(
            to_chunk_bytes(b"-1 -1 -1").unwrap(),
            Chunk::new(-1, usize::MAX, usize::MAX as i64)
        );
        assert_eq!(
            chunk_unsigned(b"18446744073709551615suffix").unwrap(),
            u64::MAX
        );
        assert_eq!(chunk_unsigned(b"-18446744073709551615").unwrap(), 1);
        for input in [
            b"1  2 3".as_slice(),
            b"1 2 nope".as_slice(),
            b"18446744073709551616 2 3".as_slice(),
            b"1 2 -18446744073709551616".as_slice(),
        ] {
            assert_eq!(to_chunk_bytes(input).unwrap_err().to_string(), "stoull");
        }
    }

    #[test]
    fn partition_loader_preserves_raw_annotations_and_unchecked_open() {
        let path = fixture_path("raw_partition_database");
        write_fixture(&path);
        let mut db = DatabaseFile::open(&path).unwrap();
        let partition_path = fixture_path("raw_partition");
        std::fs::write(&partition_path, b"4 12 9 annotation\xff\n5 20 7\xff").unwrap();
        db.load_partition_io(&partition_path).unwrap();
        assert_eq!(
            db.partition.chunks,
            [Chunk::new(4, 12, 9), Chunk::new(5, 20, 7)]
        );
        // Source constructs the failed ifstream before clearing and silently
        // leaves an empty partition; it does not report a checked InputFile error.
        std::fs::remove_file(&partition_path).unwrap();
        db.load_partition_io(&partition_path).unwrap();
        assert_eq!(db.partition, DatabasePartition::default());
        drop(db);
        std::fs::remove_file(path).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn partition_loader_accepts_raw_os_path_and_source_c_prefix() {
        use std::os::unix::ffi::OsStringExt;
        let path = fixture_path("raw_path_partition_database");
        write_fixture(&path);
        let mut db = DatabaseFile::open(&path).unwrap();
        let directory = fixture_path("raw_path_partition_dir");
        std::fs::create_dir_all(&directory).unwrap();
        let (partition_path, _) =
            crate::util::test_os_path_or_unicode(&directory, b"partition\xff", "partition-é");
        let mut bytes = partition_path.clone().into_os_string().into_vec();
        std::fs::write(&partition_path, b"3 8 2").unwrap();
        bytes.extend_from_slice(b"\0not-part-of-open-path");
        db.load_partition_io(Path::new(&std::ffi::OsString::from_vec(bytes)))
            .unwrap();
        assert_eq!(db.partition.chunks, [Chunk::new(3, 8, 2)]);
        drop(db);
        std::fs::remove_file(path).unwrap();
        std::fs::remove_file(partition_path).unwrap();
        std::fs::remove_dir_all(directory).unwrap();
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
