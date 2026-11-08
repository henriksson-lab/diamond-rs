//! Shared sequence-file behavior from `diamond/src/data/sequence_file.{h,cpp}`.
//!
//! Concrete DMND, BLAST, FASTA, and block readers remain in their source-mirrored
//! modules.  This module contains the format-independent types and algorithms
//! which the C++ `SequenceFile` base class supplies to those readers.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

use crate::basic::value::{DictId, Letter, Loc, OId, TaxId};

pub use crate::data::blastdb::volume::{DbFilter, DecodedPackage, RawChunk, SequenceFileFlags};

pub const MAX_LINEAGE: usize = 256;
pub const DICT_EMPTY: DictId = DictId::MAX;
pub const SEQID_HEADER: &str = "seqid";

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Chunk {
    pub i: i32,
    pub offset: usize,
    pub n_seqs: i64,
}

impl Chunk {
    pub const fn new(i: i32, offset: usize, n_seqs: i64) -> Self {
        Self { i, offset, n_seqs }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceFileType {
    Dmnd,
    Blast,
    Fasta,
    Block,
}

impl std::fmt::Display for SequenceFileType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Dmnd => "Diamond database",
            Self::Blast => "BLAST database",
            Self::Fasta => "FASTA file",
            Self::Block => "",
        })
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FormatFlags(pub i32);

impl FormatFlags {
    pub const NONE: Self = Self(0);
    pub const TITLES_LAZY: Self = Self(1);
    pub const DICT_LENGTHS: Self = Self(1 << 1);
    pub const DICT_SEQIDS: Self = Self(1 << 2);
    pub const LENGTH_LOOKUP: Self = Self(1 << 3);
    pub const SEEKABLE: Self = Self(1 << 4);

    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }
}

impl std::ops::BitOr for FormatFlags {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeqInfo {
    pub pos: u64,
    pub seq_len: u32,
}

impl SeqInfo {
    pub const SIZE: usize = 16;
    pub const fn new(pos: u64, seq_len: usize) -> Self {
        Self {
            pos,
            seq_len: seq_len as u32,
        }
    }
}

/// Global-free settings used by C++ `auto_create` and `total_blocks`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SequenceFileConfig {
    pub multiprocessing: bool,
    pub target_indexed: bool,
    pub chunk_size_billions: f64,
}

impl Default for SequenceFileConfig {
    fn default() -> Self {
        Self {
            multiprocessing: false,
            target_indexed: false,
            chunk_size_billions: 2.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueryStrands {
    Both,
    Plus,
    Minus,
}

/// Six-bit translated-frame mask used by C++ `load_onepass`.
pub const fn frame_mask(strands: QueryStrands) -> i32 {
    match strands {
        QueryStrands::Both => (1 << 6) - 1,
        QueryStrands::Plus => (1 << 3) - 1,
        QueryStrands::Minus => ((1 << 3) - 1) << 3,
    }
}

pub fn total_blocks(letters: u64, config: SequenceFileConfig) -> Result<usize, String> {
    let chunk = (config.chunk_size_billions * 1e9) as u64;
    if chunk == 0 {
        return Err("Chunk size must be positive.".into());
    }
    Ok(letters.div_ceil(chunk) as usize)
}

/// Detects the same file families as C++ `SequenceFile::auto_create`, without
/// constructing a concrete reader or consulting process-wide configuration.
fn append_path_suffix(path: &std::path::Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(suffix);
    PathBuf::from(name)
}
fn source_dmnd_candidate(path: &std::path::Path) -> PathBuf {
    let candidate = append_path_suffix(path, ".dmnd");
    #[cfg(unix)]
    let has_suffix = {
        use std::os::unix::ffi::OsStrExt;
        path.as_os_str().as_bytes().ends_with(b".dmnd")
    };
    #[cfg(windows)]
    let has_suffix = {
        use std::os::windows::ffi::OsStrExt;
        path.as_os_str()
            .encode_wide()
            .collect::<Vec<_>>()
            .ends_with(&[46, 100, 109, 110, 100])
    };
    #[cfg(not(any(unix, windows)))]
    let has_suffix = path.as_os_str().as_encoded_bytes().ends_with(b".dmnd");
    if !has_suffix && candidate.exists() {
        candidate
    } else {
        path.to_path_buf()
    }
}

pub fn detect_type(
    paths: &[PathBuf],
    flags: SequenceFileFlags,
    config: SequenceFileConfig,
) -> Result<(SequenceFileType, PathBuf), String> {
    if paths.len() == 1 {
        let path = &paths[0];
        let pin = append_path_suffix(path, ".pin");
        let pal = append_path_suffix(path, ".pal");
        if pin.exists() || pal.exists() || path.extension().is_some_and(|x| x == "pal") {
            if config.multiprocessing {
                return Err("--multiprocessing is not compatible with BLAST databases.".into());
            }
            if config.target_indexed {
                return Err("--target-indexed is not compatible with BLAST databases.".into());
            }
            return Ok((SequenceFileType::Blast, path.clone()));
        }
        let dmnd = source_dmnd_candidate(path);
        if crate::data::dmnd::is_diamond_db_file(&dmnd) {
            return Ok((SequenceFileType::Dmnd, dmnd));
        }
    }
    if !flags.contains(SequenceFileFlags::NO_FASTA) && !paths.is_empty() {
        return Ok((SequenceFileType::Fasta, paths[0].clone()));
    }
    Err("Sequence file does not have a supported format.".into())
}

/// Source auto_create's actual header probe. Unlike the convenience boolean
/// detector, a missing/unreadable candidate is a typed FileOpenException;
/// only EndOfStream while reading its magic is interpreted as a non-DMND file.
pub fn detect_type_source(
    paths: &[PathBuf],
    flags: SequenceFileFlags,
    config: SequenceFileConfig,
) -> std::io::Result<(SequenceFileType, PathBuf)> {
    detect_type_source_with_buffer(paths, flags, config, 67_108_864)
}

pub fn detect_type_source_with_buffer(
    paths: &[PathBuf],
    flags: SequenceFileFlags,
    config: SequenceFileConfig,
    buffer_size: usize,
) -> std::io::Result<(SequenceFileType, PathBuf)> {
    if paths.len() == 1 {
        let path = &paths[0];
        if append_path_suffix(path, ".pin").exists()
            || append_path_suffix(path, ".pal").exists()
            || path.extension().is_some_and(|x| x == "pal")
        {
            if config.multiprocessing {
                return Err(std::io::Error::other(
                    "--multiprocessing is not compatible with BLAST databases.",
                ));
            }
            if config.target_indexed {
                return Err(std::io::Error::other(
                    "--target-indexed is not compatible with BLAST databases.",
                ));
            }
            return Ok((SequenceFileType::Blast, path.clone()));
        }
        let candidate = source_dmnd_candidate(path);
        if !candidate.as_os_str().is_empty() && candidate.as_os_str() != std::ffi::OsStr::new("-") {
            let mut input =
                crate::util::io::InputFile::with_source_profile_path(&candidate, 0, buffer_size)
                    .map_err(std::io::Error::other)?;
            let magic = match input.read_u64() {
                Ok(magic) => Some(magic),
                Err(crate::util::io::IoError::EndOfStream) => None,
                Err(error) => return Err(std::io::Error::other(error)),
            };
            input.close().map_err(std::io::Error::other)?;
            if magic == Some(crate::data::dmnd::MAGIC_NUMBER) {
                return Ok((SequenceFileType::Dmnd, candidate));
            }
        }
    }
    if !flags.contains(SequenceFileFlags::NO_FASTA) && !paths.is_empty() {
        return Ok((SequenceFileType::Fasta, paths[0].clone()));
    }
    Err(std::io::Error::other(
        "Sequence file does not have a supported format.",
    ))
}

pub fn single_oid(accession: &str, oids: &[OId]) -> Result<OId, String> {
    match oids {
        [oid] => Ok(*oid),
        [] => Err(format!("Accession not found in database: {accession}")),
        _ => Err(format!("Accession is not unique in database: {accession}")),
    }
}

#[derive(Debug, Clone, Default)]
pub struct AccessionIndex {
    acc_to_oid: HashMap<Vec<u8>, OId>,
    oid_to_acc: Vec<Vec<u8>>,
}

impl AccessionIndex {
    pub fn add_seqid_mapping(&mut self, title: &str, oid: OId) -> Result<(), String> {
        self.add_seqid_mapping_bytes(title.as_bytes(), oid)
            .map_err(|error| error.to_string())
    }
    pub fn add_seqid_mapping_bytes(&mut self, title: &[u8], oid: OId) -> std::io::Result<()> {
        let accession = crate::util::sequence::seqid_bytes(title);
        if oid != self.oid_to_acc.len() as OId {
            return Err(std::io::Error::other(
                "add_seqid_mapping: OIDs must be added consecutively",
            ));
        }
        if self.acc_to_oid.contains_key(&accession) {
            let mut message = b"Accession is not unique in database file: ".to_vec();
            message.extend_from_slice(&accession);
            return Err(crate::util::io::RawDiagnosticError(message).into_io());
        }
        self.acc_to_oid.insert(accession.clone(), oid);
        self.oid_to_acc.push(accession);
        Ok(())
    }
    pub fn accession_to_oid(&self, accession: &str) -> Result<Vec<OId>, String> {
        self.accession_to_oid_bytes(accession.as_bytes())
            .map_err(|error| error.to_string())
    }
    pub fn accession_to_oid_bytes(&self, accession: &[u8]) -> std::io::Result<Vec<OId>> {
        self.acc_to_oid
            .get(accession)
            .copied()
            .map(|oid| vec![oid])
            .ok_or_else(|| {
                let mut message = b"Accession not found in database: ".to_vec();
                message.extend_from_slice(accession);
                crate::util::io::RawDiagnosticError(message).into_io()
            })
    }
    pub fn seqid_bytes(&self, oid: OId) -> Result<&[u8], String> {
        self.oid_to_acc
            .get(oid as usize)
            .map(Vec::as_slice)
            .ok_or_else(|| "OId to accession mapping not available.".into())
    }
    pub fn seqid(&self, oid: OId) -> Result<&str, String> {
        std::str::from_utf8(self.seqid_bytes(oid)?).map_err(|error| error.to_string())
    }
}

/// C++ `partition`, including its special treatment of `start_size`.
pub fn partition(
    seq_lengths: &[Loc],
    max_block_size: u64,
    start_size: u64,
) -> (Vec<OId>, Vec<u64>) {
    let mut boundaries = vec![0];
    let mut sizes = Vec::new();
    let mut block_size = start_size;
    for (i, &length) in seq_lengths.iter().enumerate() {
        let length = length as u64;
        if block_size.wrapping_add(length) > max_block_size && block_size > 0 {
            boundaries.push(i as OId);
            sizes.push(block_size.wrapping_sub(if sizes.is_empty() { start_size } else { 0 }));
            block_size = 0;
        }
        block_size = block_size.wrapping_add(length);
    }
    boundaries.push(seq_lengths.len() as OId);
    sizes.push(block_size.wrapping_sub(if sizes.is_empty() { start_size } else { 0 }));
    (boundaries, sizes)
}

pub fn letters_filtered<F>(filter: &DbFilter, mut seq_length: F) -> Result<u64, String>
where
    F: FnMut(OId) -> Result<Loc, String>,
{
    let mut n = 0u64;
    for (oid, selected) in filter.oid_filter.iter().copied().enumerate() {
        if selected {
            n += seq_length(oid as OId)?.max(0) as u64;
        }
    }
    Ok(n)
}

/// C++ `seq_offsets`; `-1` means the selected record immediately follows the
/// previous selected record and no seek is necessary.
pub fn seq_offsets(oids: &[OId], infos: &[SeqInfo]) -> Result<Vec<i64>, String> {
    if !oids.windows(2).all(|w| w[0] <= w[1]) {
        return Err("OIds must be sorted.".into());
    }
    if oids.last().is_some_and(|&oid| oid as usize >= infos.len()) {
        return Err("OId out of bounds.".into());
    }
    Ok(oids
        .iter()
        .enumerate()
        .map(|(i, &oid)| {
            if i > 0 && oid == oids[i - 1] + 1 {
                -1
            } else {
                infos[oid as usize].pos as i64
            }
        })
        .collect())
}

pub trait TaxonomyAccess {
    fn parent(&self, taxid: TaxId) -> TaxId;
    fn rank(&self, taxid: TaxId) -> i32;
    fn max_taxid(&self) -> TaxId;
}

pub fn rank_taxid<T: TaxonomyAccess>(db: &T, mut taxid: TaxId, rank: i32) -> Result<TaxId, String> {
    for _ in 0..=64 {
        if db.rank(taxid) == rank {
            return Ok(taxid);
        }
        if taxid <= 1 {
            return Ok(0);
        }
        taxid = db.parent(taxid);
    }
    Err("Path in taxonomy too long (rank_taxid).".into())
}

pub fn rank_taxids<T: TaxonomyAccess>(
    db: &T,
    taxids: &[TaxId],
    rank: i32,
) -> Result<BTreeSet<TaxId>, String> {
    taxids.iter().map(|&t| rank_taxid(db, t, rank)).collect()
}

pub fn lineage<T: TaxonomyAccess>(db: &T, mut taxid: TaxId) -> Result<Vec<TaxId>, String> {
    let mut out = Vec::new();
    for _ in 0..=MAX_LINEAGE {
        if taxid <= 0 {
            return Ok(Vec::new());
        }
        if taxid == 1 {
            out.reverse();
            return Ok(out);
        }
        out.push(taxid);
        taxid = db.parent(taxid);
    }
    Err("Path in taxonomy too long (TaxonomyNodes::lineage).".into())
}

pub fn get_lca<T: TaxonomyAccess>(db: &T, t1: TaxId, t2: TaxId) -> Result<TaxId, String> {
    if t1 == t2 || t2 <= 0 {
        return Ok(t1);
    }
    if t1 <= 0 {
        return Ok(t2);
    }
    let mut p = t2;
    let mut ancestors = HashSet::from([p]);
    let mut reached_root = false;
    for _ in 0..MAX_LINEAGE {
        p = db.parent(p);
        if p <= 0 {
            return Ok(t1);
        }
        ancestors.insert(p);
        if p == t1 {
            return Ok(p);
        }
        if p == 1 {
            reached_root = true;
            break;
        }
    }
    if !reached_root {
        return Err("Path in taxonomy too long (get_lca).".into());
    }
    p = t1;
    for _ in 0..MAX_LINEAGE {
        if ancestors.contains(&p) {
            return Ok(p);
        }
        p = db.parent(p);
        if p <= 0 {
            return Ok(t2);
        }
    }
    Err("Path in taxonomy too long (get_lca).".into())
}

pub fn contained<T: TaxonomyAccess>(
    db: &T,
    query: TaxId,
    filter: &BTreeSet<TaxId>,
    include_invalid: bool,
) -> Result<bool, String> {
    if db.parent(query) < 0 {
        return Ok(include_invalid);
    }
    if filter.contains(&1) {
        return Ok(true);
    }
    let mut p = query;
    for _ in 0..=64 {
        if p <= 1 || filter.contains(&p) {
            return Ok(p > 1);
        }
        p = db.parent(p);
        if p <= 0 {
            return Ok(include_invalid);
        }
    }
    Err("Path in taxonomy too long (contained).".into())
}

pub fn contained_many<T: TaxonomyAccess>(
    db: &T,
    query: &[TaxId],
    filter: &BTreeSet<TaxId>,
    all: bool,
    include_invalid: bool,
) -> Result<bool, String> {
    if filter.contains(&1) {
        return Ok(true);
    }
    for &taxid in query {
        let value = contained(db, taxid, filter, include_invalid)?;
        if value != all {
            return Ok(!all);
        }
    }
    Ok(all)
}

/// Source taxonomy-list parsing, shared by real database filters and
/// provided-database workflows whose metadata access can fail afterward.
pub fn parse_taxonomy_filter(list: &str, delimiter: char) -> Result<BTreeSet<TaxId>, String> {
    parse_taxonomy_filter_bytes(list.as_bytes(), delimiter as u8)
}
pub fn parse_taxonomy_filter_bytes(list: &[u8], delimiter: u8) -> Result<BTreeSet<TaxId>, String> {
    let mut taxa = BTreeSet::new();
    if !list.is_empty() {
        let mut tokens = list.split(|byte| *byte == delimiter).peekable();
        while let Some(token) = tokens.next() {
            if token.is_empty() && tokens.peek().is_none() && list.last() == Some(&delimiter) {
                break;
            }
            let mut c = token.to_vec();
            c.push(0);
            let begin = c.as_ptr().cast();
            let mut end = std::ptr::null_mut();
            // std::stoi uses the active C locale and accepts a decimal prefix.
            // Parsing into the wider CRT type and then checking int's bounds
            // also detects overflow on LLP64 without relying on stale errno.
            let value = unsafe { libc::strtoll(begin, &mut end, 10) };
            if end == begin.cast_mut() || value < i32::MIN as i64 || value > i32::MAX as i64 {
                return Err("stoi".into());
            }
            taxa.insert(value as TaxId);
        }
    }
    if taxa.is_empty() {
        return Err("Option --taxonlist/--taxon-exclude used with empty list.".into());
    }
    if taxa.contains(&0) || taxa.contains(&1) {
        return Err(
            "Option --taxonlist/--taxon-exclude used with invalid argument (0 or 1).".into(),
        );
    }
    Ok(taxa)
}

pub fn filter_by_taxonomy<T, FT, FL>(
    db: &T,
    list: &str,
    delimiter: char,
    exclude: bool,
    sequence_count: OId,
    taxids: FT,
    seq_length: FL,
) -> Result<DbFilter, String>
where
    T: TaxonomyAccess,
    FT: FnMut(OId) -> Result<Vec<TaxId>, String>,
    FL: FnMut(OId) -> Result<Loc, String>,
{
    filter_by_taxonomy_bytes(
        db,
        list.as_bytes(),
        delimiter as u8,
        exclude,
        sequence_count,
        taxids,
        seq_length,
    )
}
pub fn filter_by_taxonomy_bytes<T, FT, FL>(
    db: &T,
    list: &[u8],
    delimiter: u8,
    exclude: bool,
    sequence_count: OId,
    mut taxids: FT,
    mut seq_length: FL,
) -> Result<DbFilter, String>
where
    T: TaxonomyAccess,
    FT: FnMut(OId) -> Result<Vec<TaxId>, String>,
    FL: FnMut(OId) -> Result<Loc, String>,
{
    let taxa = parse_taxonomy_filter_bytes(list, delimiter)?;
    let mut out = DbFilter::new(sequence_count as usize);
    for oid in 0..sequence_count {
        let c = contained_many(db, &taxids(oid)?, &taxa, exclude, exclude)?;
        if c ^ exclude {
            out.oid_filter[oid as usize] = true;
            out.letter_count += seq_length(oid)?.max(0) as u64;
        }
    }
    Ok(out)
}

#[derive(Debug, Clone)]
pub struct DictionaryEntry<T = String> {
    pub oid: OId,
    pub len: Loc,
    pub title: T,
    pub seq: Vec<Letter>,
    pub self_aln_score: f64,
}

pub type ByteDictionaryEntry = DictionaryEntry<Vec<u8>>;
impl From<DictionaryEntry> for ByteDictionaryEntry {
    fn from(entry: DictionaryEntry) -> Self {
        Self {
            oid: entry.oid,
            len: entry.len,
            title: entry.title.into_bytes(),
            seq: entry.seq,
            self_aln_score: entry.self_aln_score,
        }
    }
}
impl TryFrom<ByteDictionaryEntry> for DictionaryEntry {
    type Error = std::string::FromUtf8Error;
    fn try_from(entry: ByteDictionaryEntry) -> Result<Self, Self::Error> {
        Ok(Self {
            oid: entry.oid,
            len: entry.len,
            title: String::from_utf8(entry.title)?,
            seq: entry.seq,
            self_aln_score: entry.self_aln_score,
        })
    }
}

#[derive(Debug)]
pub struct SequenceDictionary<T = String> {
    entries: Vec<DictionaryEntry<T>>,
    block_to_dict_id: HashMap<usize, Vec<DictId>>,
}

pub type ByteSequenceDictionary = SequenceDictionary<Vec<u8>>;
impl<T> Default for SequenceDictionary<T> {
    fn default() -> Self {
        Self {
            entries: Vec::new(),
            block_to_dict_id: HashMap::new(),
        }
    }
}
impl<T> SequenceDictionary<T> {
    pub fn init_block(&mut self, block: usize, seq_count: usize, persist: bool) {
        if !persist {
            self.block_to_dict_id.clear();
        }
        self.block_to_dict_id
            .entry(block)
            .or_insert_with(|| vec![DICT_EMPTY; seq_count]);
    }

    pub fn dict_id(
        &mut self,
        block: usize,
        block_id: usize,
        entry: DictionaryEntry<T>,
    ) -> Result<DictId, String> {
        let ids = self
            .block_to_dict_id
            .get_mut(&block)
            .ok_or("Dictionary not initialized.")?;
        let slot = ids.get_mut(block_id).ok_or("Dictionary not initialized.")?;
        if *slot != DICT_EMPTY {
            return Ok(*slot);
        }
        let id = self.entries.len() as DictId;
        self.entries.push(entry);
        *slot = id;
        Ok(id)
    }

    pub fn entry(&self, dict_id: DictId) -> Result<&DictionaryEntry<T>, String> {
        usize::try_from(dict_id)
            .ok()
            .and_then(|i| self.entries.get(i))
            .ok_or_else(|| "Dictionary not loaded.".into())
    }

    pub fn free(&mut self) {
        self.entries.clear();
        self.block_to_dict_id.clear();
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Checked compatibility adapter for Source's byte-oriented FAI reader.
pub fn read_fai_file(
    path: &Path,
    seqs: i64,
    letters: i64,
    index: Option<&mut AccessionIndex>,
) -> Result<(i64, i64), String> {
    read_fai_file_io(path, seqs, letters, index).map_err(|error| error.to_string())
}

/// Source SequenceFile::read_fai_file: retain TextInputFile's physical EOF and
/// CR handling, and the old CharDelimiter tokenizer's local C-string extent.
pub fn read_fai_file_io(
    path: &Path,
    seqs: i64,
    letters: i64,
    index: Option<&mut AccessionIndex>,
) -> std::io::Result<(i64, i64)> {
    read_fai_file_with_buffer_size_io(path, seqs, letters, index, 67_108_864)
}

/// Explicit replacement for Source's global config.file_buffer_size.
pub fn read_fai_file_with_buffer_size_io(
    path: &Path,
    mut seqs: i64,
    mut letters: i64,
    mut index: Option<&mut AccessionIndex>,
    buffer_size: usize,
) -> std::io::Result<(i64, i64)> {
    let mut input =
        crate::util::io::text_input_file::ByteTextInputFile::from_source_path_bytes_with_buffer_size(path, buffer_size)
            .map_err(std::io::Error::other)?;
    loop {
        input.getline().map_err(std::io::Error::other)?;
        if input.line.is_empty() && input.eof() {
            break;
        }
        let mut tokens = crate::util::string::ByteTokenizer::new(&input.line, b'\t');
        let accession = tokens.read_bytes().map_err(std::io::Error::other)?;
        let len = tokens.read_i32().map_err(std::io::Error::other)?;
        if let Some(index) = index.as_deref_mut() {
            // The FAI reader assigns the complete token directly, replacing
            // duplicate keys. It does not call add_seqid_mapping or populate
            // the independently requested OID-to-accession vector.
            index.acc_to_oid.insert(accession, seqs as OId);
        }
        seqs += 1;
        letters += i64::from(len);
    }
    input.into_inner().close().map_err(std::io::Error::other)?;
    Ok((seqs, letters))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fai_bytes_keep_accessions_and_source_assignment_semantics() {
        let path = std::env::temp_dir().join(format!(
            "diamond-fai-bytes-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(
            &path,
            b"acc\xff\t4\0ignored\x80\nacc\xff\t7\tstuff\nother\x80\t-2\r\n",
        )
        .unwrap();
        let mut index = AccessionIndex::default();
        let result = read_fai_file_io(&path, 5, 100, Some(&mut index));
        std::fs::remove_file(&path).unwrap();
        assert_eq!(result.unwrap(), (8, 109));
        assert_eq!(index.accession_to_oid_bytes(b"acc\xff").unwrap(), vec![6]);
        assert_eq!(index.accession_to_oid_bytes(b"other\x80").unwrap(), vec![7]);
        assert!(index.seqid_bytes(5).is_err());
    }

    struct Tree {
        parent: Vec<TaxId>,
        ranks: Vec<i32>,
    }
    impl TaxonomyAccess for Tree {
        fn parent(&self, t: TaxId) -> TaxId {
            self.parent.get(t as usize).copied().unwrap_or(-1)
        }
        fn rank(&self, t: TaxId) -> i32 {
            self.ranks.get(t as usize).copied().unwrap_or(-1)
        }
        fn max_taxid(&self) -> TaxId {
            self.parent.len() as TaxId - 1
        }
    }
    fn tree() -> Tree {
        Tree {
            parent: vec![-1, 1, 1, 2, 2, 3],
            ranks: vec![-1, 0, 10, 20, 20, 30],
        }
    }

    #[test]
    fn taxonomy_filter_bytes_keep_stoi_prefix_and_int_range() {
        assert_eq!(
            parse_taxonomy_filter_bytes(b" +7\xff,9\0ignored,", b',').unwrap(),
            BTreeSet::from([7, 9])
        );
        assert_eq!(
            parse_taxonomy_filter_bytes(b"\xc2\xa07", b',').unwrap_err(),
            "stoi"
        );
        assert_eq!(
            parse_taxonomy_filter_bytes(b"2147483648", b',').unwrap_err(),
            "stoi"
        );
        assert_eq!(
            parse_taxonomy_filter_bytes(b"999999999999999999999", b',').unwrap_err(),
            "stoi"
        );
        assert_eq!(
            parse_taxonomy_filter_bytes(b"7,,", b',').unwrap_err(),
            "stoi"
        );
        assert_eq!(
            parse_taxonomy_filter_bytes(b"", b',').unwrap_err(),
            "Option --taxonlist/--taxon-exclude used with empty list."
        );
    }

    #[test]
    fn cpp_partition_preserves_initial_size_semantics() {
        assert_eq!(partition(&[4, 7, 3], 10, 2), (vec![0, 1, 3], vec![4, 10]));
        assert_eq!(partition(&[4, 3], 10, 2), (vec![0, 2], vec![7]));
        // Source converts signed Loc to uint64_t before partition arithmetic.
        assert_eq!(partition(&[-1, 2], 10, 0), (vec![0, 2], vec![1]));
        assert_eq!(partition(&[-1], 10, 2), (vec![0, 1], vec![u64::MAX]));
    }

    #[test]
    fn taxonomy_algorithms_match_base_class() {
        let db = tree();
        assert_eq!(lineage(&db, 5).unwrap(), vec![2, 3, 5]);
        assert_eq!(get_lca(&db, 5, 4).unwrap(), 2);
        assert_eq!(rank_taxid(&db, 5, 10).unwrap(), 2);
        assert_eq!(
            rank_taxids(&db, &[4, 5], 20).unwrap(),
            BTreeSet::from([3, 4])
        );
        assert!(contained(&db, 5, &BTreeSet::from([2]), false).unwrap());
        assert!(!contained_many(&db, &[5, 4], &BTreeSet::from([3]), true, false).unwrap());
    }

    #[test]
    fn taxonomy_filter_include_and_exclude() {
        let db = tree();
        let tax = [vec![5], vec![4], vec![-1]];
        let inc = filter_by_taxonomy(
            &db,
            "3",
            ',',
            false,
            3,
            |i| Ok(tax[i as usize].clone()),
            |_| Ok(10),
        )
        .unwrap();
        assert_eq!(inc.oid_filter, vec![true, false, false]);
        assert_eq!(inc.letter_count, 10);
        let exc = filter_by_taxonomy(
            &db,
            "3",
            ',',
            true,
            3,
            |i| Ok(tax[i as usize].clone()),
            |_| Ok(10),
        )
        .unwrap();
        assert_eq!(exc.oid_filter, vec![false, true, false]);
    }

    #[test]
    fn accession_index_rejects_duplicates_and_gaps() {
        let mut idx = AccessionIndex::default();
        idx.add_seqid_mapping("P12345 name", 0).unwrap();
        assert_eq!(idx.accession_to_oid("P12345").unwrap(), vec![0]);
        assert_eq!(idx.seqid(0).unwrap(), "P12345");
        assert!(idx.add_seqid_mapping("P12345 duplicate", 1).is_err());
    }

    #[test]
    fn dictionary_ids_are_stable_within_block() {
        let mut dict = SequenceDictionary::<String>::default();
        dict.init_block(7, 2, false);
        let entry = || DictionaryEntry {
            oid: 9,
            len: 2,
            title: "x".into(),
            seq: vec![1, 2],
            self_aln_score: 4.0,
        };
        assert_eq!(dict.dict_id(7, 0, entry()).unwrap(), 0);
        assert_eq!(dict.dict_id(7, 0, entry()).unwrap(), 0);
        assert_eq!(dict.entry(0).unwrap().oid, 9);
        assert_eq!(dict.len(), 1);
    }

    #[test]
    fn offsets_mark_adjacent_records() {
        let infos = [
            SeqInfo::new(10, 2),
            SeqInfo::new(20, 3),
            SeqInfo::new(40, 4),
        ];
        assert_eq!(seq_offsets(&[0, 1, 2], &infos).unwrap(), vec![10, -1, -1]);
        assert_eq!(seq_offsets(&[0, 2], &infos).unwrap(), vec![10, 40]);
    }

    #[test]
    fn explicit_chunk_configuration() {
        let cfg = SequenceFileConfig {
            chunk_size_billions: 0.000000010,
            ..Default::default()
        };
        assert_eq!(total_blocks(21, cfg).unwrap(), 3);
        assert_eq!(frame_mask(QueryStrands::Both), 0b11_1111);
        assert_eq!(frame_mask(QueryStrands::Plus), 0b00_0111);
        assert_eq!(frame_mask(QueryStrands::Minus), 0b11_1000);
    }
}

/// Source SequenceFile::load_seqs first protein-block warning. Call only for
/// the first logical block; later chunks do not repeat the ten-record probe.
const PROTEIN_DNA_ERROR: &str = "The sequences are expected to be proteins but only contain DNA letters. Use the option --ignore-warnings to proceed.";

/// Incremental first-block warning from Source's one-pass load loop. The
/// tenth all-DNA record fails before the next input record is parsed; a
/// shorter all-DNA block fails when that block finishes.
pub struct ProteinInputProbe {
    ignore: bool,
    count: usize,
    dna: usize,
}
impl ProteinInputProbe {
    pub fn new(ignore_warnings: bool) -> Self {
        Self {
            ignore: ignore_warnings,
            count: 0,
            dna: 0,
        }
    }
    pub fn observe(&mut self, sequence: &[crate::basic::value::Letter]) -> Result<(), String> {
        if self.ignore || self.count >= 10 || sequence.is_empty() {
            return Ok(());
        }
        self.count += 1;
        let traits = crate::basic::value::ValueTraits::new(
            crate::basic::value::AMINO_ACID_ALPHABET,
            crate::basic::value::MASK_LETTER,
            b" \t\r\n",
            crate::basic::value::SequenceType::AminoAcid,
        );
        if sequence.iter().all(|letter| {
            ((*letter & crate::basic::value::LETTER_MASK) as usize) < traits.alphabet.len()
        }) && crate::util::sequence::looks_like_dna(sequence, &traits)?
        {
            self.dna += 1;
            if self.dna == 10 {
                return Err(PROTEIN_DNA_ERROR.to_owned());
            }
        }
        Ok(())
    }
    pub fn finish_block(&self) -> Result<(), String> {
        if self.count != 0 && self.count == self.dna {
            Err(PROTEIN_DNA_ERROR.to_owned())
        } else {
            Ok(())
        }
    }
}

pub fn check_protein_input<'a>(
    sequences: impl IntoIterator<Item = &'a [crate::basic::value::Letter]>,
    ignore_warnings: bool,
) -> Result<(), String> {
    let mut probe = ProteinInputProbe::new(ignore_warnings);
    for sequence in sequences
        .into_iter()
        .filter(|sequence| !sequence.is_empty())
        .take(10)
    {
        probe.observe(sequence)?;
    }
    probe.finish_block()
}
