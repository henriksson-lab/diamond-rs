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
pub fn detect_type(
    paths: &[PathBuf],
    flags: SequenceFileFlags,
    config: SequenceFileConfig,
) -> Result<(SequenceFileType, PathBuf), String> {
    if paths.len() == 1 {
        let path = &paths[0];
        let path_text = path.to_string_lossy();
        let pin = PathBuf::from(format!("{path_text}.pin"));
        let pal = PathBuf::from(format!("{path_text}.pal"));
        if pin.exists() || pal.exists() || path.extension().is_some_and(|x| x == "pal") {
            if config.multiprocessing {
                return Err("--multiprocessing is not compatible with BLAST databases.".into());
            }
            if config.target_indexed {
                return Err("--target-indexed is not compatible with BLAST databases.".into());
            }
            return Ok((SequenceFileType::Blast, path.clone()));
        }
        let dmnd = if path.exists() {
            path.clone()
        } else {
            PathBuf::from(format!("{path_text}.dmnd"))
        };
        if crate::data::dmnd::is_diamond_db_file(&dmnd) {
            return Ok((SequenceFileType::Dmnd, dmnd));
        }
    }
    if !flags.contains(SequenceFileFlags::NO_FASTA) && !paths.is_empty() {
        return Ok((SequenceFileType::Fasta, paths[0].clone()));
    }
    Err("Sequence file does not have a supported format.".into())
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
    acc_to_oid: HashMap<String, OId>,
    oid_to_acc: Vec<String>,
}

impl AccessionIndex {
    pub fn add_seqid_mapping(&mut self, title: &str, oid: OId) -> Result<(), String> {
        let accession = crate::util::sequence::seqid(title);
        if oid != self.oid_to_acc.len() as OId {
            return Err("add_seqid_mapping: OIDs must be added consecutively".into());
        }
        if self.acc_to_oid.contains_key(&accession) {
            return Err(format!(
                "Accession is not unique in database file: {accession}"
            ));
        }
        self.acc_to_oid.insert(accession.clone(), oid);
        self.oid_to_acc.push(accession);
        Ok(())
    }

    pub fn accession_to_oid(&self, accession: &str) -> Result<Vec<OId>, String> {
        self.acc_to_oid
            .get(accession)
            .copied()
            .map(|x| vec![x])
            .ok_or_else(|| format!("Accession not found in database: {accession}"))
    }

    pub fn seqid(&self, oid: OId) -> Result<&str, String> {
        self.oid_to_acc
            .get(oid as usize)
            .map(String::as_str)
            .ok_or_else(|| "OId to accession mapping not available.".into())
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
        let length = length.max(0) as u64;
        if block_size + length > max_block_size && block_size > 0 {
            boundaries.push(i as OId);
            sizes.push(block_size - if sizes.is_empty() { start_size } else { 0 });
            block_size = 0;
        }
        block_size += length;
    }
    boundaries.push(seq_lengths.len() as OId);
    sizes.push(block_size - if sizes.is_empty() { start_size } else { 0 });
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

pub fn filter_by_taxonomy<T, FT, FL>(
    db: &T,
    list: &str,
    delimiter: char,
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
    let taxa: BTreeSet<TaxId> = list
        .split(delimiter)
        .filter(|x| !x.trim().is_empty())
        .map(|x| x.trim().parse::<TaxId>().map_err(|e| e.to_string()))
        .collect::<Result<_, _>>()?;
    if taxa.is_empty() {
        return Err("Option --taxonlist/--taxon-exclude used with empty list.".into());
    }
    if taxa.contains(&0) || taxa.contains(&1) {
        return Err(
            "Option --taxonlist/--taxon-exclude used with invalid argument (0 or 1).".into(),
        );
    }
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
pub struct DictionaryEntry {
    pub oid: OId,
    pub len: Loc,
    pub title: String,
    pub seq: Vec<Letter>,
    pub self_aln_score: f64,
}

#[derive(Debug, Default)]
pub struct SequenceDictionary {
    entries: Vec<DictionaryEntry>,
    block_to_dict_id: HashMap<usize, Vec<DictId>>,
}

impl SequenceDictionary {
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
        entry: DictionaryEntry,
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

    pub fn entry(&self, dict_id: DictId) -> Result<&DictionaryEntry, String> {
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

pub fn read_fai_file(
    path: &Path,
    mut seqs: i64,
    mut letters: i64,
    index: Option<&mut AccessionIndex>,
) -> Result<(i64, i64), String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut index = index;
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let accession = fields.next().ok_or("Invalid FAI record")?;
        let len: Loc = fields
            .next()
            .ok_or("Invalid FAI record")?
            .parse::<Loc>()
            .map_err(|e| e.to_string())?;
        if let Some(idx) = index.as_deref_mut() {
            idx.add_seqid_mapping(accession, seqs as OId)?;
        }
        seqs += 1;
        letters += i64::from(len);
    }
    Ok((seqs, letters))
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn cpp_partition_preserves_initial_size_semantics() {
        assert_eq!(partition(&[4, 7, 3], 10, 2), (vec![0, 1, 3], vec![4, 10]));
        assert_eq!(partition(&[4, 3], 10, 2), (vec![0, 2], vec![7]));
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
        let mut dict = SequenceDictionary::default();
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
