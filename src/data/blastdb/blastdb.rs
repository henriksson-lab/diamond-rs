//! BLAST database sequence-file implementation.
//!
//! Mirrors `diamond/src/data/blastdb/blastdb.cpp` and `blastdb.h`.

use super::pal::{absolute_parts, path_bytes, raw_error};
use crate::basic::value::{Letter, Loc, TaxId};
use crate::data::blastdb::taxdmp::{
    read_names_dmp_bytes_with_source_buffer_path, read_nodes_dmp_bytes_with_source_buffer_path,
};
use crate::data::blastdb::volume::{
    build_title_bytes, BlastDefLine, BlastVolume, ByteBlastDefLine, DbFilter, Pal, RawChunk,
    SequenceFileFlags,
};
use crate::data::taxonomy::Rank;
use crate::util::io::source_c_path;
#[cfg(test)]
use crate::util::system::PATH_SEPARATOR;
use rusqlite::{Connection, OpenFlags, OptionalExtension};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

fn path_error(prefix: &[u8], path: &Path, suffix: &[u8]) -> std::io::Error {
    let mut message = prefix.to_vec();
    let bytes = match path_bytes(path) {
        Ok(bytes) => bytes,
        Err(error) => return error,
    };
    message.extend_from_slice(&bytes);
    message.extend_from_slice(suffix);
    raw_error(message)
}

/// TextInputFile's Source byte-line contract, including codec autodetection,
/// one trailing CR removal, interior empty lines and an unterminated final line.
fn read_accession_list_source(
    file_name: &Path,
    buffer_size: usize,
) -> std::io::Result<HashMap<Vec<u8>, bool>> {
    let mut input = crate::util::io::InputFile::with_source_profile_path(file_name, 0, buffer_size)
        .map_err(std::io::Error::other)?
        .into_deserializer();
    let mut accessions = HashMap::new();
    loop {
        let mut line = Vec::new();
        let terminated = input
            .read_to(&mut line, b'\n')
            .map_err(std::io::Error::other)?;
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        if line.is_empty() && !terminated {
            break;
        }
        accessions.insert(line, false);
        if !terminated {
            break;
        }
    }
    input.close().map_err(std::io::Error::other)?;
    Ok(accessions)
}

fn match_accession_bytes(
    deflines: &[ByteBlastDefLine],
    accessions: &mut HashMap<Vec<u8>, bool>,
) -> bool {
    for defline in deflines {
        for seqid in &defline.seqids {
            if let Some(found) = accessions.get_mut(&seqid.value) {
                *found = true;
                return true;
            }
            if seqid.version.is_some() || seqid.chain.is_some() {
                let formatted = super::volume::format_seqid_bytes(seqid);
                if let Some(found) = accessions.get_mut(&formatted) {
                    *found = true;
                    return true;
                }
            }
        }
    }
    false
}

fn check_missing_accessions(
    accessions: &HashMap<Vec<u8>, bool>,
    skip_missing: bool,
) -> std::io::Result<()> {
    if !skip_missing {
        for (accession, found) in accessions {
            if !found {
                let mut message = b"Accession not found in database: ".to_vec();
                message.extend_from_slice(accession);
                message.extend_from_slice(b". Use --skip-missing-seqids to ignore.");
                return Err(crate::util::io::RawDiagnosticError(message).into_io());
            }
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BlastDbConfig {
    pub multiprocessing: bool,
}

struct SqliteConnection(Connection);

impl std::fmt::Debug for SqliteConnection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SqliteConnection").finish_non_exhaustive()
    }
}

impl SqliteConnection {
    fn open_readonly(path: impl AsRef<Path>) -> Result<Self, String> {
        Self::open(path.as_ref(), OpenFlags::SQLITE_OPEN_READ_ONLY)
    }

    fn open(path: &Path, flags: OpenFlags) -> Result<Self, String> {
        if path.as_os_str().as_encoded_bytes().contains(&0) {
            return Err("SQLite path contains a null byte".into());
        }
        Connection::open_with_flags(source_c_path(&(path)).as_ref(), flags)
            .map(Self)
            .map_err(|error| error.to_string())
    }

    fn max_taxid(&self) -> Result<TaxId, String> {
        let mut statement = self
            .0
            .prepare("SELECT max(taxid) FROM TaxidInfo;")
            .map_err(|error| format!("Failed to prepare statement: {error}"))?;
        statement
            .query_row([], |row| row.get::<_, Option<TaxId>>(0))
            .map(|value| value.unwrap_or(0))
            .map_err(|error| format!("SQLite step error: {error}"))
    }

    fn parent(&self, taxid: TaxId) -> Result<TaxId, String> {
        let mut statement = self
            .0
            .prepare("SELECT parent FROM TaxidInfo WHERE taxid = ?1 LIMIT 1;")
            .map_err(|error| format!("Failed to prepare statement: {error}"))?;
        statement
            .query_row([taxid], |row| row.get(0))
            .optional()
            .map(|value| value.unwrap_or(-1))
            .map_err(|error| format!("SQLite step error: {error}"))
    }

    #[cfg(test)]
    fn create_test_database(path: &str) -> Result<(), String> {
        let connection = Self::open(
            Path::new(path),
            OpenFlags::SQLITE_OPEN_READ_WRITE | OpenFlags::SQLITE_OPEN_CREATE,
        )?;
        connection
            .0
            .execute_batch(
                "CREATE TABLE TaxidInfo(taxid INTEGER PRIMARY KEY, parent INTEGER); \
                 INSERT INTO TaxidInfo VALUES(1,1),(3,1),(7,3),(42,7);",
            )
            .map_err(|error| error.to_string())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SeqInfo {
    pub pos: u64,
    pub seq_len: u32,
}

impl SeqInfo {
    pub const SIZE: usize = 16;

    pub fn new(pos: u64, len: usize) -> Self {
        Self {
            pos,
            seq_len: len as u32,
        }
    }
}

#[derive(Debug)]
pub struct BlastDB {
    pub(crate) accession_spool:
        Option<std::sync::Arc<crate::data::accession_spool::AccessionSpool>>,
    file_name: PathBuf,
    pal: Pal,
    taxon_db: Option<SqliteConnection>,
    taxon_mapping: Vec<(u64, TaxId)>,
    taxon_mapping_sorted: bool,
    custom_ranks: BTreeMap<Vec<u8>, i32>,
    rank_mapping: HashMap<TaxId, i32>,
    oid: u64,
    long_seqids: bool,
    flags: SequenceFileFlags,
    parent_cache: Vec<TaxId>,
    extra_names: HashMap<TaxId, Vec<u8>>,
    volumes: BTreeMap<u64, PathBuf>,
    io_errors: crate::util::io::IoErrorOwner,
    dict_oid: Vec<Vec<u64>>,
    volume: BlastVolume,
    raw_chunk_no: i32,
    seq_length: Vec<Loc>,
}

impl BlastDB {
    pub fn new(file_name: &str, flags: SequenceFileFlags) -> Result<Self, String> {
        Self::new_path(Path::new(file_name), flags).map_err(|e| e.to_string())
    }
    pub fn new_path(file_name: &Path, flags: SequenceFileFlags) -> std::io::Result<Self> {
        Self::new_path_with_io_owner(file_name, flags, &crate::util::io::IoErrorOwner::default())
    }
    pub fn new_path_with_io_owner(
        file_name: &Path,
        flags: SequenceFileFlags,
        errors: &crate::util::io::IoErrorOwner,
    ) -> std::io::Result<Self> {
        let result = Self::new_path_owned(file_name, flags, errors);
        errors.finish(result)
    }
    fn new_path_owned(
        file_name: &Path,
        flags: SequenceFileFlags,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        let pal = errors.capture(Pal::new_path(file_name))?;
        let volume = BlastVolume::new_path_owned(
            &pal.volumes[0],
            0,
            pal.oid_index[0],
            pal.oid_index[1],
            true,
            errors,
        )?;
        let mut db = Self {
            accession_spool: None,
            file_name: file_name.to_path_buf(),
            pal,
            taxon_db: None,
            taxon_mapping: Vec::new(),
            taxon_mapping_sorted: true,
            custom_ranks: BTreeMap::new(),
            rank_mapping: HashMap::new(),
            oid: 0,
            long_seqids: false,
            flags,
            parent_cache: Vec::new(),
            extra_names: HashMap::new(),
            volumes: BTreeMap::new(),
            io_errors: errors.clone(),
            dict_oid: Vec::new(),
            volume,
            raw_chunk_no: 0,
            seq_length: Vec::new(),
        };

        if db.pal.metadata.contains_key("SEQIDLIST") {
            db.flags |= SequenceFileFlags::NEED_LENGTH_LOOKUP;
        }
        if db.pal.metadata.contains_key("TAXIDLIST") {
            db.flags |= SequenceFileFlags::TAXON_MAPPING;
            db.flags |= SequenceFileFlags::TAXON_NODES;
            db.flags |=
                SequenceFileFlags::NEED_EARLY_TAXON_MAPPING | SequenceFileFlags::NEED_LENGTH_LOOKUP;
        }

        if db.flags.contains(SequenceFileFlags::TAXON_MAPPING) {
            let flags_now = db.flags;
            for oid in 0..db.pal.sequence_count {
                let deflines = db.deflines_bytes(oid, true, false, true)?;
                // Source's constructor obtains this mapping through
                // RawChunk::decode (pin.cpp): one sorted set of taxids per
                // OID, even when multiple deflines share the same taxid.
                let taxids: BTreeSet<_> = deflines
                    .into_iter()
                    .filter_map(|defline| defline.taxid)
                    .collect();
                db.taxon_mapping
                    .extend(taxids.into_iter().map(|taxid| (oid, taxid)));
            }
            db.flags = flags_now;
            db.flags &= !SequenceFileFlags::TAXON_MAPPING;
        }

        if db.flags.contains(SequenceFileFlags::NEED_LENGTH_LOOKUP) {
            db.seq_length.reserve(db.pal.sequence_count as usize);
            for i in 0..db.pal.sequence_count {
                db.open_volume(i)?;
                let volume_oid = (i - db.volume.begin) as u32;
                db.seq_length.push(db.volume.length(volume_oid));
            }
        }

        let (dbdir, _dbfile) = errors.capture(absolute_parts(file_name))?;
        if db.flags.contains(SequenceFileFlags::TAXON_RANKS) {
            let file = dbdir.join("nodes.dmp");
            if !source_c_path(&file).exists() {
                return errors.capture(Err(path_error(b"Taxonomy rank information (nodes.dmp) is missing in search path (",&dbdir,b"). Download and extract this file in the database directory: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip")));
            }
            let mut next_rank = Rank::COUNT as i32;
            let result = read_nodes_dmp_bytes_with_source_buffer_path(
                &file,
                1 << 20,
                |taxid, _parent, rank| {
                    if let Some(p) = Rank::predefined_bytes(rank) {
                        db.rank_mapping.insert(taxid, p as i32);
                    } else if let Some(&r) = db.custom_ranks.get(rank) {
                        db.rank_mapping.insert(taxid, r);
                    } else {
                        let r = next_rank;
                        next_rank += 1;
                        db.custom_ranks.insert(rank.to_vec(), r);
                        db.rank_mapping.insert(taxid, r);
                    }
                },
            );
            errors.capture(result)?;
        }

        if db.flags.contains(SequenceFileFlags::TAXON_NODES) {
            let dbpath = dbdir.join("taxonomy4blast.sqlite3");
            if !source_c_path(&dbpath).exists() {
                return errors.capture(Err(path_error(
                    b"Taxonomy database (taxonomy4blast.sqlite3) file not found in path: ",
                    &dbpath,
                    b". Make sure that the database was downloaded correctly.",
                )));
            }
            let connection =
                errors.capture(SqliteConnection::open_readonly(&dbpath).map_err(|error| {
                    path_error(
                        b"Failed to open database ",
                        &dbpath,
                        format!(": {error}").as_bytes(),
                    )
                }))?;
            let max_id = connection.max_taxid()?;
            db.parent_cache = vec![TaxId::MIN; max_id as usize + 1];
            db.taxon_db = Some(connection);
        }

        if db.flags.contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES) {
            let file = dbdir.join("names.dmp");
            if !source_c_path(&file).exists() {
                return errors.capture(Err(path_error(b"Taxonomy names information (names.dmp) is missing in search path (",&dbdir,b"). Download and extract this file in the database directory: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip")));
            }
            let result =
                read_names_dmp_bytes_with_source_buffer_path(&file, 1 << 20, |taxid, name| {
                    db.extra_names.entry(taxid).or_insert_with(|| name.to_vec());
                });
            errors.capture(result)?;
        }

        Ok(db)
    }

    pub fn new_with_config(
        file_name: &str,
        flags: SequenceFileFlags,
        config: BlastDbConfig,
    ) -> Result<Self, String> {
        let db = Self::new(file_name, flags)?;
        if config.multiprocessing {
            return Err("Multiprocessing mode is not compatible with BLAST databases.".to_string());
        }
        Ok(db)
    }

    pub fn new_path_with_config(
        file_name: &Path,
        flags: SequenceFileFlags,
        config: BlastDbConfig,
    ) -> std::io::Result<Self> {
        let db = Self::new_path(file_name, flags)?;
        if config.multiprocessing {
            return Err(std::io::Error::other(
                "Multiprocessing mode is not compatible with BLAST databases.",
            ));
        }
        Ok(db)
    }

    pub fn file_count(&self) -> i64 {
        1
    }

    pub fn open_volume(&mut self, oid: u64) -> Result<(), String> {
        if oid >= self.volume.begin && oid < self.volume.end {
            return Ok(());
        }
        let idx = self.pal.volume(oid) as usize;
        self.volume = BlastVolume::new_path_owned(
            &self.pal.volumes[idx],
            idx as i32,
            self.pal.oid_index[idx],
            self.pal.oid_index[idx + 1],
            true,
            &self.io_errors,
        )?;
        Ok(())
    }

    pub fn print_info(&self) -> String {
        String::from_utf8(
            self.print_info_bytes()
                .expect("Source name encoding failed; use print_info_bytes"),
        )
        .expect("BLAST database name is not UTF-8; use print_info_bytes")
    }
    pub fn print_info_bytes(&self) -> std::io::Result<Vec<u8>> {
        let mut out = format!(
            " (type: BLAST database, volumes: {}, sequences: {}, letters: {})\n",
            self.pal.volumes.len(),
            self.sequence_count(),
            self.letters()
        );
        if self.flags.contains(SequenceFileFlags::TAXON_RANKS) {
            if !self.custom_ranks.is_empty() {
                out.push_str(&format!(
                    "Custom taxonomic ranks in database: {}\n",
                    self.custom_ranks.len()
                ));
            }
            out.push_str(&format!(
                "Taxonomic ids assigned to ranks: {}\n",
                self.rank_mapping.len()
            ));
        }
        if self.flags.contains(SequenceFileFlags::TAXON_NODES) && !self.parent_cache.is_empty() {
            out.push_str(&format!(
                "Maximum taxid in database: {}\n",
                self.parent_cache.len() - 1
            ));
        }
        if self
            .flags
            .contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES)
        {
            out.push_str(&format!(
                "Extra taxonomic scientific names in names.dmp: {}\n",
                self.extra_names.len()
            ));
        }
        let mut bytes = b"Database: ".to_vec();
        bytes.extend_from_slice(&path_bytes(&self.file_name)?);
        bytes.extend_from_slice(out.as_bytes());
        Ok(bytes)
    }

    pub fn init_seqinfo_access(&mut self) {}

    pub fn init_seq_access(&mut self) {
        self.oid = 0;
    }

    pub fn seek_chunk(&mut self, _chunk: &Chunk) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn tell_seq(&self) -> u64 {
        self.oid
    }

    pub fn eof(&self) -> bool {
        self.oid == self.sequence_count()
    }

    pub fn read_seqinfo(&mut self) -> Result<SeqInfo, String> {
        if self.oid >= self.pal.sequence_count {
            self.oid += 1;
            return Ok(SeqInfo::new(0, 0));
        }
        let l = self.seq_length(self.oid as usize)?;
        if l == 0 {
            return Err("Database with sequence length 0 is not supported".to_string());
        }
        let out = SeqInfo::new(self.oid, l as usize);
        self.oid += 1;
        Ok(out)
    }

    pub fn putback_seqinfo(&mut self) {
        self.oid -= 1;
    }

    pub fn id_len(
        &mut self,
        seq_info: &SeqInfo,
        _seq_info_next: &SeqInfo,
    ) -> Result<usize, String> {
        self.open_volume(seq_info.pos)?;
        let volume_oid = (seq_info.pos - self.volume.begin) as u32;
        Ok(self.volume.id_len(volume_oid))
    }

    pub fn seek_offset(&mut self, _p: usize) {}

    pub fn raw_chunk(
        &mut self,
        letters: usize,
        flags: SequenceFileFlags,
    ) -> Result<RawChunk, String> {
        let mut c = self.volume.raw_chunk(letters, flags)?;
        c.no = self.raw_chunk_no;
        self.raw_chunk_no += 1;
        let oid = c.end;
        if oid < self.pal.sequence_count {
            self.open_volume(oid)?;
        }
        Ok(c)
    }

    pub fn read_seq_data(
        &mut self,
        dst: &mut Vec<Letter>,
        len: usize,
        pos: &mut usize,
        _seek: bool,
    ) -> Result<(), String> {
        self.open_volume(*pos as u64)?;
        let volume_oid = (*pos as u64 - self.volume.begin) as u32;
        let seq = self.volume.sequence(volume_oid)?;
        dst.clear();
        dst.reserve(len + 2);
        dst.push(crate::basic::value::DELIMITER_LETTER);
        dst.extend_from_slice(&seq);
        dst.push(crate::basic::value::DELIMITER_LETTER);
        *pos += 1;
        Ok(())
    }

    pub fn read_id_data(
        &mut self,
        oid: i64,
        all: bool,
        full_titles: bool,
    ) -> Result<String, String> {
        String::from_utf8(self.read_id_data_bytes(oid, all, full_titles)?)
            .map_err(|error| error.to_string())
    }

    pub fn read_id_data_bytes(
        &mut self,
        oid: i64,
        all: bool,
        full_titles: bool,
    ) -> Result<Vec<u8>, String> {
        self.fetch_seqid_bytes(oid as u64, all, full_titles)
    }

    pub fn deflines(
        &mut self,
        oid: u64,
        all: bool,
        full_titles: bool,
        taxids: bool,
    ) -> Result<Vec<BlastDefLine>, String> {
        self.deflines_bytes(oid, all, full_titles, taxids)?
            .into_iter()
            .map(ByteBlastDefLine::try_into_string)
            .collect()
    }

    pub fn deflines_bytes(
        &mut self,
        oid: u64,
        all: bool,
        full_titles: bool,
        taxids: bool,
    ) -> Result<Vec<ByteBlastDefLine>, String> {
        self.open_volume(oid)?;
        let volume_oid = (oid - self.volume.begin) as u32;
        self.volume
            .deflines_bytes(volume_oid, all, full_titles, taxids)
    }

    pub fn skip_id_data(&mut self) {}

    pub fn fetch_seqid(
        &mut self,
        oid: u64,
        all: bool,
        full_titles: bool,
    ) -> Result<String, String> {
        String::from_utf8(self.fetch_seqid_bytes(oid, all, full_titles)?)
            .map_err(|error| error.to_string())
    }

    pub fn fetch_seqid_bytes(
        &mut self,
        oid: u64,
        all: bool,
        full_titles: bool,
    ) -> Result<Vec<u8>, String> {
        self.open_volume(oid)?;
        let taxids = self.flags.contains(SequenceFileFlags::TAXON_MAPPING);
        let volume_oid = (oid - self.volume.begin) as u32;
        let deflines = self
            .volume
            .deflines_bytes(volume_oid, all, full_titles, taxids)?;
        if taxids && !self.taxon_mapping.iter().any(|&(o, _)| o == oid) {
            for i in &deflines {
                if let Some(taxid) = i.taxid {
                    self.add_taxid_mapping(&[(oid, taxid)]);
                }
            }
        }
        Ok(build_title_bytes(&deflines, b"\x01", all))
    }

    pub fn add_taxid_mapping(&mut self, taxids: &[(u64, TaxId)]) {
        if let (Some(last), Some(first)) = (self.taxon_mapping.last(), taxids.first()) {
            self.taxon_mapping_sorted &= last.0 <= first.0;
        }
        self.taxon_mapping_sorted &= taxids.windows(2).all(|pair| pair[0].0 <= pair[1].0);
        self.taxon_mapping.extend_from_slice(taxids);
    }

    /// Stream defline text without populating the optional taxonomy cache.
    /// Accession-list finalization needs titles only, even when the search
    /// owner also has TAXON_MAPPING enabled for actual matched subjects.
    pub fn title_only(&mut self, oid: u64, all: bool, full_titles: bool) -> Result<String, String> {
        String::from_utf8(self.title_only_bytes(oid, all, full_titles)?)
            .map_err(|error| error.to_string())
    }

    pub fn title_only_bytes(
        &mut self,
        oid: u64,
        all: bool,
        full_titles: bool,
    ) -> Result<Vec<u8>, String> {
        self.open_volume(oid)?;
        let volume_oid = (oid - self.volume.begin) as u32;
        let deflines = self
            .volume
            .deflines_bytes(volume_oid, all, full_titles, false)?;
        Ok(build_title_bytes(&deflines, b"\x01", all))
    }

    pub fn seqid(&mut self, oid: u64, all: bool, full_titles: bool) -> Result<String, String> {
        String::from_utf8(self.seqid_bytes(oid, all, full_titles)?)
            .map_err(|error| error.to_string())
    }

    pub fn seqid_bytes(
        &mut self,
        oid: u64,
        all: bool,
        full_titles: bool,
    ) -> Result<Vec<u8>, String> {
        self.fetch_seqid_bytes(oid, all, full_titles)
    }

    pub fn dict_seq(&mut self, dict_id: i64, _ref_block: usize) -> Result<Vec<Letter>, String> {
        let oid = self
            .dict_oid
            .first()
            .and_then(|block| usize::try_from(dict_id).ok().and_then(|id| block.get(id)))
            .copied()
            .ok_or_else(|| "Dictionary not loaded.".to_string())?;
        let mut sequence = Vec::new();
        self.seq_data(oid as usize, &mut sequence)?;
        Ok(sequence)
    }

    pub fn sequence_count(&self) -> u64 {
        self.pal.sequence_count
    }

    pub fn letters(&self) -> u64 {
        self.pal.letters
    }

    pub fn db_version(&self) -> i32 {
        self.pal.version
    }

    pub fn program_build_version(&self) -> i32 {
        0
    }

    pub fn read_seq(&mut self) -> Result<(Vec<Letter>, String), String> {
        let (sequence, title) = self.read_seq_bytes()?;
        Ok((
            sequence,
            String::from_utf8(title).map_err(|error| error.to_string())?,
        ))
    }

    /// Source deflines are byte strings; this is the canonical sequential reader.
    pub fn read_seq_bytes(&mut self) -> Result<(Vec<Letter>, Vec<u8>), String> {
        self.open_volume(self.oid)?;
        let volume_oid = (self.oid - self.volume.begin) as u32;
        let seq = self.volume.sequence(volume_oid)?;
        let deflines = self.volume.deflines_bytes(volume_oid, true, true, false)?;
        let id = build_title_bytes(&deflines, b" >", true);
        self.oid += 1;
        Ok((seq, id))
    }

    pub fn build_version(&self) -> i32 {
        0
    }

    pub fn create_partition_balanced(&mut self, _max_letters: i64) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn save_partition(
        &mut self,
        _partition_file_name: &str,
        _annotation: &str,
    ) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn get_n_partition_chunks(&mut self) -> Result<i32, String> {
        Err("Operation not supported".to_string())
    }

    pub fn set_seqinfo_ptr(&mut self, i: u64) -> Result<(), String> {
        if i != 0 {
            return Err(
                "Setting seqinfo pointer to non-zero value is not supported in BLAST databases."
                    .to_string(),
            );
        }
        self.oid = i;
        self.raw_chunk_no = 0;
        if self.volume.begin == 0 {
            self.volume.rewind()
        } else {
            self.open_volume(0)
        }
    }

    pub fn close(&mut self) {
        self.taxon_db = None;
    }

    pub fn filter_by_accession(
        &mut self,
        file_name: &str,
        skip_missing_seqids: bool,
    ) -> Result<DbFilter, String> {
        self.filter_by_accession_path(Path::new(file_name), skip_missing_seqids)
    }
    pub fn filter_by_accession_path(
        &mut self,
        file_name: &Path,
        skip_missing_seqids: bool,
    ) -> Result<DbFilter, String> {
        let mut v = DbFilter::new(self.sequence_count() as usize);
        let mut source_accessions =
            read_accession_list_source(file_name, 64 << 20).map_err(|error| error.to_string())?;
        self.set_seqinfo_ptr(0)?;
        loop {
            let chunk = self.raw_chunk(
                1_000_000_000,
                SequenceFileFlags::SEQS
                    | SequenceFileFlags::TITLES
                    | SequenceFileFlags::FULL_TITLES,
            )?;
            let pkg = chunk.decode_bytes(
                SequenceFileFlags::SEQS | SequenceFileFlags::FULL_TITLES,
                None,
                Some(&mut source_accessions),
            )?;
            for oid in pkg.oids {
                if let Some(slot) = v.oid_filter.get_mut(oid as usize) {
                    *slot = true;
                }
                v.letter_count += self.seq_length(oid as usize)? as u64;
            }
            if chunk.end >= self.sequence_count() {
                break;
            }
        }

        check_missing_accessions(&source_accessions, skip_missing_seqids)
            .map_err(|error| error.to_string())?;
        Ok(v)
    }

    /// Decode one global OID without changing the sequential logical cursor.
    /// The per-volume decoder seeks independently when random ranking reads
    /// revisit an earlier volume, retaining bounded sequence ownership.
    pub fn record_at(&mut self, oid: u64) -> Result<(Vec<Letter>, String), String> {
        let (sequence, title) = self.record_at_bytes(oid)?;
        Ok((
            sequence,
            String::from_utf8(title).map_err(|error| error.to_string())?,
        ))
    }

    pub fn record_at_bytes(&mut self, oid: u64) -> Result<(Vec<Letter>, Vec<u8>), String> {
        if oid >= self.sequence_count() {
            return Err("OID exceeds number of sequences in database".to_owned());
        }
        let mut sequence = Vec::new();
        self.seq_data(oid as usize, &mut sequence)?;
        let title = self.seqid_bytes(oid, true, true)?;
        Ok((sequence, title))
    }

    /// Accession filtering only needs deflines and indexed lengths. Avoid the
    /// old whole raw-sequence chunk allocation while constructing a PAL mask.
    pub fn filter_by_accession_metadata(
        &mut self,
        file_name: &str,
        skip_missing_seqids: bool,
    ) -> Result<DbFilter, String> {
        self.filter_by_accession_metadata_io(file_name, skip_missing_seqids, 64 << 20)
            .map_err(|error| error.to_string())
    }

    pub fn filter_by_accession_metadata_io(
        &mut self,
        file_name: &str,
        skip_missing_seqids: bool,
        buffer_size: usize,
    ) -> std::io::Result<DbFilter> {
        self.filter_by_accession_metadata_path_io(
            Path::new(file_name),
            skip_missing_seqids,
            buffer_size,
        )
    }
    pub fn filter_by_accession_metadata_path_io(
        &mut self,
        file_name: &Path,
        skip_missing_seqids: bool,
        buffer_size: usize,
    ) -> std::io::Result<DbFilter> {
        let mut accessions = read_accession_list_source(file_name, buffer_size)?;
        let mut filter = DbFilter::new(self.sequence_count() as usize);
        for oid in 0..self.sequence_count() {
            let deflines = self
                .deflines_bytes(
                    oid,
                    self.flags.contains(SequenceFileFlags::ALL_SEQIDS),
                    self.flags.contains(SequenceFileFlags::FULL_TITLES),
                    false,
                )
                .map_err(std::io::Error::other)?;
            if match_accession_bytes(&deflines, &mut accessions) {
                filter.oid_filter[oid as usize] = true;
                filter.letter_count += self
                    .seq_length(oid as usize)
                    .map_err(std::io::Error::other)?
                    .max(0) as u64;
            }
        }
        check_missing_accessions(&accessions, skip_missing_seqids)?;
        Ok(filter)
    }

    /// Apply the alias filters in the same order as Search::run. SEQIDLIST
    /// replaces TAXIDLIST when both occur in the alias, as in the source.
    pub fn alias_filter(&mut self) -> Result<Option<DbFilter>, String> {
        self.search_filter("", "", "", false)
    }

    /// Search::run filter precedence, with bounded metadata-only access.
    pub fn search_filter(
        &mut self,
        include: &str,
        exclude: &str,
        seqidlist: &str,
        skip_missing: bool,
    ) -> Result<Option<DbFilter>, String> {
        self.search_filter_io(include, exclude, seqidlist, skip_missing, 64 << 20)
            .map_err(|error| error.to_string())
    }

    /// Actual search invocation's file buffer and fatal I/O metadata owner.
    pub fn search_filter_io(
        &mut self,
        include: &str,
        exclude: &str,
        seqidlist: &str,
        skip_missing: bool,
        buffer_size: usize,
    ) -> std::io::Result<Option<DbFilter>> {
        self.search_filter_path_io(
            include,
            exclude,
            Path::new(seqidlist),
            skip_missing,
            buffer_size,
        )
    }
    pub fn search_filter_path_io(
        &mut self,
        include: &str,
        exclude: &str,
        seqidlist: &Path,
        skip_missing: bool,
        buffer_size: usize,
    ) -> std::io::Result<Option<DbFilter>> {
        self.search_filter_bytes_path_io(
            include.as_bytes(),
            exclude.as_bytes(),
            seqidlist,
            skip_missing,
            buffer_size,
        )
    }
    pub fn search_filter_bytes_path_io(
        &mut self,
        include: &[u8],
        exclude: &[u8],
        seqidlist: &Path,
        skip_missing: bool,
        buffer_size: usize,
    ) -> std::io::Result<Option<DbFilter>> {
        let errors = self.io_errors.clone();
        let result = self.search_filter_owned(
            include,
            exclude,
            seqidlist,
            skip_missing,
            buffer_size,
            &errors,
        );
        errors.finish(result)
    }

    fn search_filter_owned(
        &mut self,
        include: &[u8],
        exclude: &[u8],
        seqidlist: &Path,
        skip_missing: bool,
        buffer_size: usize,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<Option<DbFilter>, String> {
        if !include.is_empty() && !exclude.is_empty() {
            return Err("Options --taxonlist and --taxon-exclude are mutually exclusive.".into());
        }
        let explicit_taxa = !include.is_empty() || !exclude.is_empty();
        let alias_accessions = self
            .pal
            .metadata
            .get("SEQIDLIST")
            .map(|value| crate::util::io::source_bytes_os(value).map(PathBuf::from))
            .transpose();
        let alias_accessions = errors.capture(alias_accessions.map_err(std::io::Error::other))?;
        let accessions = if seqidlist.as_os_str().is_empty() {
            alias_accessions
        } else {
            Some(seqidlist.to_owned())
        };
        let taxonomy = if explicit_taxa {
            Some((
                if include.is_empty() { exclude } else { include }.to_owned(),
                b',',
                !exclude.is_empty(),
            ))
        } else if let Some(path) = self.pal.metadata.get("TAXIDLIST") {
            Some((
                {
                    let path = errors.capture(
                        crate::util::io::source_bytes_os(path)
                            .map(PathBuf::from)
                            .map_err(std::io::Error::other),
                    )?;
                    let mut input = errors
                        .capture(std::fs::File::open(source_c_path(&path).as_ref()).map_err(
                            |_| path_error(b"Cannot open TAXIDLIST file: ", &path, b""),
                        ))?;
                    let mut bytes = Vec::new();
                    // Source checks the ifstream constructor, then streams its
                    // rdbuf without checking later read failures.
                    let _ = std::io::Read::read_to_end(&mut input, &mut bytes);
                    bytes
                },
                b'\n',
                false,
            ))
        } else {
            None
        };
        let mut filter = None;
        if let Some((values, delimiter, exclude)) = taxonomy {
            // The generic taxonomy predicate uses immutable access; the source
            // parent lookup fills a cache, so keep that mutation scoped here.
            struct Access<'a>(std::cell::RefCell<&'a mut BlastDB>);
            impl crate::data::sequence_file::TaxonomyAccess for Access<'_> {
                fn parent(&self, taxid: TaxId) -> TaxId {
                    self.0.borrow_mut().get_parent(taxid)
                }
                fn rank(&self, taxid: TaxId) -> i32 {
                    self.0.borrow().rank(taxid)
                }
                fn max_taxid(&self) -> TaxId {
                    self.0.borrow().max_taxid()
                }
            }
            let count = self.sequence_count();
            let access = Access(std::cell::RefCell::new(&mut *self));
            filter = Some(crate::data::sequence_file::filter_by_taxonomy_bytes(
                &access,
                &values,
                delimiter,
                exclude,
                count,
                |oid| Ok(access.0.borrow().taxids(oid as usize)),
                |oid| access.0.borrow_mut().seq_length(oid as usize),
            )?);
        }
        if !seqidlist.as_os_str().is_empty() && self.pal.metadata.contains_key("SEQIDLIST") {
            return Err("Using --seqidlist on already filtered BLAST alias database.".into());
        }
        if let Some(path) = accessions {
            if explicit_taxa {
                return Err("--seqidlist is not compatible with taxonomy filtering.".into());
            }
            filter = Some(errors.capture(self.filter_by_accession_metadata_path_io(
                &path,
                skip_missing,
                buffer_size,
            ))?);
        }
        Ok(filter)
    }

    pub fn file_path(&self) -> &Path {
        &self.file_name
    }
    pub fn try_file_name(&self) -> std::io::Result<String> {
        self.file_name
            .to_str()
            .map(str::to_owned)
            .ok_or_else(|| std::io::Error::other("BLAST database name is not UTF-8; use file_path"))
    }
    pub fn file_name(&self) -> String {
        self.try_file_name()
            .expect("BLAST database name is not UTF-8; use file_path or try_file_name")
    }
    pub fn source_io<T>(
        &mut self,
        f: impl FnOnce(&mut Self) -> Result<T, String>,
    ) -> std::io::Result<T> {
        let errors = self.io_errors.clone();
        let result = f(self);
        errors.finish(result)
    }
    pub fn record_at_bytes_io(&mut self, oid: u64) -> std::io::Result<(Vec<Letter>, Vec<u8>)> {
        self.source_io(|db| db.record_at_bytes(oid))
    }

    pub fn taxids(&self, oid: usize) -> Vec<TaxId> {
        let mut taxids: Vec<_> = if self.taxon_mapping_sorted {
            let begin = self
                .taxon_mapping
                .partition_point(|&(entry, _)| entry < oid as u64);
            self.taxon_mapping[begin..]
                .iter()
                .take_while(|&&(entry, _)| entry == oid as u64)
                .map(|&(_, taxid)| taxid)
                .collect()
        } else {
            // Explicit external mappings may arrive in arbitrary OID order.
            self.taxon_mapping
                .iter()
                .filter_map(|&(entry, taxid)| (entry == oid as u64).then_some(taxid))
                .collect()
        };
        // The immutable reference oracle uses GNU libstdc++ 11.4:
        // unordered_multimap::emplace without a hint prepends equivalent keys,
        // and rehash preserves their relative order. The vector above stores
        // insertion order; expose the same per-OID iteration order as Source.
        taxids.reverse();
        taxids
    }

    pub fn max_taxid(&self) -> TaxId {
        self.taxon_db
            .as_ref()
            .map(|database| database.max_taxid().expect("failed to query maximum taxid"))
            .unwrap_or(0)
    }

    pub fn get_parent(&mut self, taxid: TaxId) -> TaxId {
        if taxid <= 0 {
            return taxid;
        }
        if taxid as usize >= self.parent_cache.len() {
            return -1;
        }
        let cached = self.parent_cache[taxid as usize];
        if cached != TaxId::MIN {
            return cached;
        }
        let result = self
            .taxon_db
            .as_ref()
            .expect("taxonomy cache exists without an SQLite database")
            .parent(taxid)
            .expect("failed to query parent taxid");
        self.parent_cache[taxid as usize] = result;
        result
    }

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
        self.extra_names
            .get(&taxid)
            .cloned()
            .unwrap_or_else(|| taxid.to_string().into_bytes())
    }

    pub fn seq_data(&mut self, oid: usize, dst: &mut Vec<Letter>) -> Result<(), String> {
        self.open_volume(oid as u64)?;
        let volume_oid = (oid as u64 - self.volume.begin) as u32;
        *dst = self.volume.sequence(volume_oid)?;
        Ok(())
    }

    pub fn seq_length(&mut self, oid: usize) -> Result<Loc, String> {
        if oid < self.seq_length.len() {
            Ok(self.seq_length[oid])
        } else {
            self.open_volume(oid as u64)?;
            let volume_oid = (oid as u64 - self.volume.begin) as u32;
            Ok(self.volume.length(volume_oid))
        }
    }

    pub fn end_random_access(&mut self, dictionary: bool) {
        if dictionary {
            self.dict_oid.clear();
            self.dict_oid.shrink_to_fit();
        }
    }

    pub fn accession_to_oid_bytes(&self, acc: &[u8]) -> std::io::Result<Vec<u64>> {
        let mut message = b"Accession not found in database: ".to_vec();
        message.extend_from_slice(acc);
        Err(crate::util::io::RawDiagnosticError(message).into_io())
    }

    pub fn accession_to_oid(&self, acc: &str) -> Result<Vec<u64>, String> {
        Err(format!("Accession not found in database: {acc}"))
    }

    pub fn init_write(&mut self) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn write_seq(&mut self, _seq: &[Letter], _id: &str) -> Result<(), String> {
        Err("Operation not supported".to_string())
    }

    pub fn rank(&self, taxid: TaxId) -> i32 {
        self.rank_mapping.get(&taxid).copied().unwrap_or(-1)
    }

    pub fn pal(&self) -> &Pal {
        &self.pal
    }

    pub fn long_seqids(&self) -> bool {
        self.long_seqids
    }

    pub fn volumes(&self) -> &BTreeMap<u64, PathBuf> {
        &self.volumes
    }

    pub fn flags(&self) -> SequenceFileFlags {
        self.flags
    }

    /// Supplies the OID dictionary normally assembled by the sequence-file loader.
    pub fn set_dictionary_oids(&mut self, blocks: Vec<Vec<u64>>) {
        self.dict_oid = blocks;
    }

    pub fn raw_chunk_no(&self) -> i32 {
        self.raw_chunk_no
    }
}

impl Drop for BlastDB {
    fn drop(&mut self) {
        self.close();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_be32(out: &mut Vec<u8>, value: u32) {
        out.extend_from_slice(&value.to_be_bytes());
    }

    fn write_le64(out: &mut Vec<u8>, value: u64) {
        out.extend_from_slice(&value.to_le_bytes());
    }

    fn write_pascal(out: &mut Vec<u8>, value: &str) {
        write_be32(out, value.len() as u32);
        out.extend_from_slice(value.as_bytes());
    }

    fn sample_header(title: &str, acc: &str, version: u16, taxid: u16) -> Vec<u8> {
        sample_header_bytes(title.as_bytes(), acc.as_bytes(), version, taxid)
    }

    fn sample_header_bytes(title: &[u8], acc: &[u8], version: u16, taxid: u16) -> Vec<u8> {
        let mut seqid_inner = Vec::new();
        seqid_inner.extend_from_slice(&[0xa1, (2 + acc.len()) as u8, 0x1a, acc.len() as u8]);
        seqid_inner.extend_from_slice(acc);
        seqid_inner.extend_from_slice(&[0xa3, 0x04, 0x02, 0x02]);
        seqid_inner.extend_from_slice(&version.to_be_bytes());

        let mut seqid = Vec::new();
        seqid.extend_from_slice(&[0xa4, seqid_inner.len() as u8]);
        seqid.extend_from_slice(&seqid_inner);

        let mut seqid_wrap = Vec::new();
        seqid_wrap.extend_from_slice(&[0x30, seqid.len() as u8]);
        seqid_wrap.extend_from_slice(&seqid);

        let mut defline = Vec::new();
        defline.extend_from_slice(&[0xa0, (2 + title.len()) as u8, 0x1a, title.len() as u8]);
        defline.extend_from_slice(title);
        defline.extend_from_slice(&[0xa1, seqid_wrap.len() as u8]);
        defline.extend_from_slice(&seqid_wrap);
        defline.extend_from_slice(&[0xa2, 0x04, 0x02, 0x02]);
        defline.extend_from_slice(&taxid.to_be_bytes());

        let mut defline_wrap = Vec::new();
        defline_wrap.extend_from_slice(&[0x30, defline.len() as u8]);
        defline_wrap.extend_from_slice(&defline);

        let mut root = Vec::new();
        root.extend_from_slice(&[0x30, defline_wrap.len() as u8]);
        root.extend_from_slice(&defline_wrap);
        root
    }

    fn write_volume_files(
        prefix: &std::path::Path,
        num_oids: u32,
        header_index: &[u32],
        sequence_index: &[u32],
        phr: &[u8],
        psq: &[u8],
    ) {
        let mut pin = Vec::new();
        write_be32(&mut pin, 4);
        write_be32(&mut pin, 1);
        write_pascal(&mut pin, "test title");
        write_pascal(&mut pin, "2026-05-14");
        write_be32(&mut pin, num_oids);
        write_le64(&mut pin, 5);
        write_be32(&mut pin, 4);
        for &x in header_index {
            write_be32(&mut pin, x);
        }
        for &x in sequence_index {
            write_be32(&mut pin, x);
        }

        std::fs::write(prefix.with_extension("pin"), pin).unwrap();
        std::fs::write(prefix.with_extension("phr"), phr).unwrap();
        std::fs::write(prefix.with_extension("psq"), psq).unwrap();
    }

    fn temp_dir(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "diamond-rs-db-{name}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&dir).unwrap();
        dir
    }

    #[cfg(unix)]
    #[test]
    fn raw_paths_survive_pal_volume_suffixes_accession_filter_and_diagnostics() {
        use std::os::unix::ffi::OsStringExt;
        let dir = temp_dir("raw-path");
        let (prefix, raw_names) = crate::util::test_os_path_or_unicode(&dir, b"vol\xff", "vol-é");
        let header = sample_header("title", "ACC", 1, 7);
        write_volume_files(
            &prefix,
            1,
            &[0, header.len() as u32],
            &[0, 4],
            &header,
            &[0, 1, 2, 0],
        );
        let component = |raw: &[u8], unicode: &str| {
            if raw_names {
                dir.join(std::ffi::OsString::from_vec(raw.to_vec()))
            } else {
                dir.join(unicode)
            }
        };
        let alias = component(b"alias\xfe.pal", "alias-λ.pal");
        let ids = component(b"ids\xfd", "ids-é");
        std::fs::write(&ids, b"ACC.1\n").unwrap();
        use std::os::unix::ffi::OsStrExt;
        let mut alias_bytes = b"TITLE raw\xfc\nDBLIST ".to_vec();
        alias_bytes.extend_from_slice(prefix.file_name().unwrap().as_bytes());
        alias_bytes.extend_from_slice(b"\nSEQIDLIST ");
        alias_bytes.extend_from_slice(ids.file_name().unwrap().as_bytes());
        alias_bytes.push(b'\n');
        std::fs::write(&alias, &alias_bytes).unwrap();
        let mut db = BlastDB::new_path(&alias, SequenceFileFlags::ALL).unwrap();
        assert_eq!(db.file_path(), alias);
        assert_eq!(db.pal.volumes, vec![prefix]);
        assert_eq!(db.pal.metadata["TITLE"], b"raw\xfc");
        assert_eq!(db.record_at_bytes_io(0).unwrap().0, vec![0, 20]);
        let filter = db
            .search_filter_bytes_path_io(b"", b"", Path::new(""), false, 64)
            .unwrap()
            .unwrap();
        assert_eq!(filter.oid_filter, vec![true]);
        drop(db);
        let missing = dir.join(std::ffi::OsString::from_vec(b"missing\xfb".to_vec()));
        let error = BlastDB::new_path(&missing, SequenceFileFlags::NONE).unwrap_err();
        let diagnostic = crate::util::io::RawDiagnosticError::find(&error).unwrap();
        let mut expected = b"Error opening file ".to_vec();
        expected.extend_from_slice(
            &path_bytes(&super::super::pal::append_suffix(&missing, ".phr")).unwrap(),
        );
        expected.extend_from_slice(b". ");
        assert!(diagnostic.starts_with(&expected));
        assert!(crate::util::io::file_open_error::find(&error).is_none());
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn sqlite_connection_uses_readonly_queries_and_preserves_sentinels() {
        let dir = temp_dir("sqlite");
        let populated_path = dir.join("taxonomy4blast.sqlite3");
        let populated = populated_path.to_str().unwrap();
        SqliteConnection::create_test_database(populated).unwrap();

        {
            let connection = SqliteConnection::open_readonly(populated).unwrap();
            assert_eq!(connection.max_taxid().unwrap(), 42);
            assert_eq!(connection.parent(1).unwrap(), 1);
            assert_eq!(connection.parent(3).unwrap(), 1);
            assert_eq!(connection.parent(7).unwrap(), 3);
            assert_eq!(connection.parent(42).unwrap(), 7);
            assert_eq!(connection.parent(8).unwrap(), -1);
            assert!(connection
                .0
                .execute("INSERT INTO TaxidInfo VALUES(99, 1)", [])
                .is_err());
        }

        let empty_path = dir.join("empty.sqlite3");
        {
            let connection = SqliteConnection::open(
                &empty_path,
                OpenFlags::SQLITE_OPEN_READ_WRITE | OpenFlags::SQLITE_OPEN_CREATE,
            )
            .unwrap();
            connection
                .0
                .execute(
                    "CREATE TABLE TaxidInfo(taxid INTEGER PRIMARY KEY, parent INTEGER)",
                    [],
                )
                .unwrap();
            assert_eq!(connection.max_taxid().unwrap(), 0);
            assert_eq!(connection.parent(1).unwrap(), -1);
        }

        assert_eq!(
            SqliteConnection::open_readonly("bad\0path").unwrap_err(),
            "SQLite path contains a null byte"
        );

        let wrong_schema_path = dir.join("wrong-schema.sqlite3");
        {
            let connection = Connection::open(&wrong_schema_path).unwrap();
            connection
                .execute("CREATE TABLE Other(value INTEGER)", [])
                .unwrap();
        }
        let wrong_schema =
            SqliteConnection::open_readonly(wrong_schema_path.to_str().unwrap()).unwrap();
        assert!(wrong_schema
            .max_taxid()
            .unwrap_err()
            .starts_with("Failed to prepare statement:"));

        let corrupt_path = dir.join("corrupt.sqlite3");
        std::fs::write(&corrupt_path, b"not a sqlite database").unwrap();
        let corrupt = SqliteConnection::open_readonly(corrupt_path.to_str().unwrap()).unwrap();
        assert!(corrupt
            .max_taxid()
            .unwrap_err()
            .starts_with("Failed to prepare statement:"));

        // Windows does not allow deleting a database while read-only SQLite
        // connections still hold file handles. Explicitly end those lifetimes
        // before the common cleanup below.
        drop(wrong_schema);
        drop(corrupt);
        std::fs::remove_file(populated_path).unwrap();
        std::fs::remove_file(empty_path).unwrap();
        std::fs::remove_file(wrong_schema_path).unwrap();
        std::fs::remove_file(corrupt_path).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }

    /// Validates the reader against NCBI's actual `taxdb.tar.gz` artifact.
    ///
    /// Download and extract `https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz`,
    /// then run:
    /// `DIAMOND_REAL_TAXONOMY_DB=/path/to/taxonomy4blast.sqlite3 cargo test --lib validate_real_ncbi_taxonomy4blast_database -- --ignored`
    #[test]
    #[ignore = "requires NCBI's external taxonomy4blast.sqlite3 artifact"]
    fn validate_real_ncbi_taxonomy4blast_database() {
        let path = std::env::var("DIAMOND_REAL_TAXONOMY_DB")
            .expect("DIAMOND_REAL_TAXONOMY_DB must name taxonomy4blast.sqlite3");
        let connection = SqliteConnection::open_readonly(&path).unwrap();

        // These stable NCBI taxonomy relationships exercise two distant parts
        // of the real table rather than accepting a merely compatible schema.
        assert_eq!(connection.parent(1).unwrap(), 1);
        assert_eq!(connection.parent(2).unwrap(), 131_567);
        assert_eq!(connection.parent(562).unwrap(), 561);
        assert_eq!(connection.parent(9_606).unwrap(), 9_605);
        assert!(connection.max_taxid().unwrap() > 1_000_000);
        assert_eq!(connection.parent(i32::MAX).unwrap(), -1);
    }

    /// Exercises taxonomy lookup through a BLAST database produced by NCBI
    /// `makeblastdb`, not through the synthetic volume writer used above.
    /// The database directory must also contain NCBI's
    /// `taxonomy4blast.sqlite3` from `taxdb.tar.gz`.
    #[test]
    #[ignore = "requires an external NCBI makeblastdb database and taxdb artifact"]
    fn validate_real_ncbi_blastdb_taxonomy_mapping_and_parents() {
        let prefix = std::env::var("DIAMOND_REAL_BLASTDB_PREFIX")
            .expect("DIAMOND_REAL_BLASTDB_PREFIX must name an NCBI makeblastdb prefix");
        let mut database = BlastDB::new(
            &prefix,
            SequenceFileFlags::TAXON_MAPPING | SequenceFileFlags::TAXON_NODES,
        )
        .unwrap();

        assert_eq!(database.sequence_count(), 2);
        assert_eq!(database.taxids(0), vec![9_606]);
        assert_eq!(database.taxids(1), vec![562]);
        assert_eq!(database.get_parent(9_606), 9_605);
        assert_eq!(database.get_parent(562), 561);
        assert!(database.max_taxid() > 1_000_000);
    }

    /// Compares the translated taxonomy formatters with rows emitted by the
    /// upstream DIAMOND binary from the same real NCBI BLAST database.
    /// `DIAMOND_REAL_TAXONOMY_OUTPUT` must contain outfmt columns
    /// `qseqid sseqid staxids sscinames slineages sskingdoms skingdoms sphylums`.
    #[test]
    #[ignore = "requires external NCBI taxonomy artifacts and upstream DIAMOND output"]
    fn validate_real_ncbi_taxonomy_output_fields_against_upstream() {
        use crate::data::taxonomy::{TaxonomyNode, TaxonomyTree};
        use crate::output::format::{
            print_lineage, print_rank_taxon_names, print_staxids, print_taxon_names,
        };

        let prefix = std::env::var("DIAMOND_REAL_BLASTDB_PREFIX")
            .expect("DIAMOND_REAL_BLASTDB_PREFIX must name an NCBI makeblastdb prefix");
        let expected_path = std::env::var("DIAMOND_REAL_TAXONOMY_OUTPUT")
            .expect("DIAMOND_REAL_TAXONOMY_OUTPUT must name upstream tabular output");
        let mut database = BlastDB::new(
            &prefix,
            SequenceFileFlags::TAXON_MAPPING
                | SequenceFileFlags::TAXON_NODES
                | SequenceFileFlags::TAXON_RANKS
                | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES,
        )
        .unwrap();

        let mut tree = TaxonomyTree::new();
        let mut seen = std::collections::BTreeSet::new();
        for oid in 0..database.sequence_count() as usize {
            for mut taxid in database.taxids(oid) {
                while taxid > 0 && seen.insert(taxid) {
                    let parent = database.get_parent(taxid);
                    let rank_index = database.rank(taxid);
                    let rank = usize::try_from(rank_index)
                        .ok()
                        .and_then(Rank::from_index)
                        .map_or_else(String::new, |rank| rank.name().to_owned());
                    tree.add_node(TaxonomyNode {
                        taxid,
                        parent,
                        rank,
                        name: database.taxon_scientific_name(taxid),
                    });
                    if taxid == 1 || parent == taxid {
                        break;
                    }
                    taxid = parent;
                }
            }
        }

        let mut expected_by_taxid = std::collections::BTreeMap::new();
        let expected = std::fs::read_to_string(expected_path).unwrap();
        for row in expected.lines() {
            let cells = row.split('\t').collect::<Vec<_>>();
            assert_eq!(cells.len(), 8, "unexpected upstream taxonomy row: {row}");
            expected_by_taxid.insert(cells[2].parse::<TaxId>().unwrap(), cells[2..].join("\t"));
        }

        for oid in 0..database.sequence_count() as usize {
            let taxids = database.taxids(oid);
            assert_eq!(
                taxids.len(),
                1,
                "fixture must assign one taxid per sequence"
            );
            let actual = [
                print_staxids(&taxids, false),
                print_taxon_names(taxids.iter().copied(), &tree, false),
                print_lineage(&taxids, &tree, false),
                print_rank_taxon_names(&taxids, &tree, "superkingdom", false),
                print_rank_taxon_names(&taxids, &tree, "kingdom", false),
                print_rank_taxon_names(&taxids, &tree, "phylum", false),
            ]
            .join("\t");
            assert_eq!(actual, expected_by_taxid[&taxids[0]]);
        }
    }

    /// Diagnostic wrapper-overhead benchmark. Run with:
    /// `cargo test --release --lib benchmark_sqlite_taxonomy_open_and_parent -- --ignored --nocapture`
    #[test]
    #[ignore = "release-mode diagnostic benchmark"]
    fn benchmark_sqlite_taxonomy_open_and_parent() {
        use std::hint::black_box;
        use std::time::Instant;

        const OPEN_ITERATIONS: usize = 1_000;
        const LOOKUP_ITERATIONS: usize = 100_000;
        let dir = temp_dir("sqlite-benchmark");
        let path = dir.join("taxonomy4blast.sqlite3");
        let path = path.to_str().unwrap();
        SqliteConnection::create_test_database(path).unwrap();

        let started = Instant::now();
        let mut open_checksum = 0;
        for _ in 0..OPEN_ITERATIONS {
            let connection = black_box(SqliteConnection::open_readonly(black_box(path))).unwrap();
            open_checksum += black_box(connection.max_taxid().unwrap()) as usize;
        }
        let open_elapsed = started.elapsed();

        let connection = SqliteConnection::open_readonly(path).unwrap();
        let started = Instant::now();
        let mut parent_checksum = 0i64;
        for index in 0..LOOKUP_ITERATIONS {
            let taxid = [1, 3, 7, 42, 8][index % 5];
            parent_checksum += black_box(connection.parent(black_box(taxid)).unwrap()) as i64;
        }
        let parent_elapsed = started.elapsed();

        assert_eq!(open_checksum, OPEN_ITERATIONS * 42);
        assert_eq!(parent_checksum, 220_000);
        eprintln!(
            "taxonomy SQLite: {OPEN_ITERATIONS} cold opens+max={open_elapsed:?} ({:.3} us/op); \
             {LOOKUP_ITERATIONS} parent queries={parent_elapsed:?} ({:.3} us/op)",
            open_elapsed.as_secs_f64() * 1e6 / OPEN_ITERATIONS as f64,
            parent_elapsed.as_secs_f64() * 1e6 / LOOKUP_ITERATIONS as f64,
        );

        drop(connection);
        std::fs::remove_file(path).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn opaque_deflines_survive_indexed_sequential_and_chunk_readers() {
        // Source phr.cpp assigns BER bytes directly to std::string. In
        // particular, stored strings preserve NUL; only const-char delimiters
        // in build_title use C-string termination.
        let dir = temp_dir("raw-deflines");
        let prefix = dir.join("db");
        let first = sample_header_bytes(b"description\xff\0tail", b"ACC\x80", 1, 7);
        let second = sample_header_bytes(b"second\xfe", b"ALT\xff", 2, 9);
        let mut header = vec![0x30, (first.len() + second.len() - 4) as u8];
        header.extend_from_slice(&first[2..]);
        header.extend_from_slice(&second[2..]);
        write_volume_files(
            &prefix,
            1,
            &[0, header.len() as u32],
            &[0, 4],
            &header,
            &[0, 1, 2, 0],
        );
        let mut db = BlastDB::new(prefix.to_str().unwrap(), SequenceFileFlags::NONE).unwrap();
        let expected = b"ACC\x80.1 description\xff\0tail\x01ALT\xff.2 second\xfe";
        assert_eq!(db.seqid_bytes(0, true, true).unwrap(), expected);
        assert_eq!(db.title_only_bytes(0, true, true).unwrap(), expected);
        assert_eq!(db.record_at_bytes(0).unwrap().1, expected);
        assert!(db.seqid(0, true, true).is_err());
        assert_eq!(
            db.read_seq_bytes().unwrap().1,
            b"ACC\x80.1 description\xff\0tail >ALT\xff.2 second\xfe"
        );
        db.set_seqinfo_ptr(0).unwrap();
        let flags = SequenceFileFlags::SEQS
            | SequenceFileFlags::TITLES
            | SequenceFileFlags::FULL_TITLES
            | SequenceFileFlags::ALL_SEQIDS
            | SequenceFileFlags::TAXON_MAPPING;
        let chunk = db.raw_chunk(100, flags).unwrap();
        let mut wanted = HashMap::from([(b"ALT\xff.2".to_vec(), false)]);
        let decoded = chunk.decode_bytes(flags, None, Some(&mut wanted)).unwrap();
        assert!(wanted[b"ALT\xff.2".as_slice()]);
        assert_eq!(decoded.oids, vec![0]);
        assert_eq!(decoded.ids.get(0), expected);
        assert_eq!(decoded.taxids, vec![(0, 7), (0, 9)]);
        drop(db);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn test_blastdb_read_seq_seqinfo_raw_chunk_and_filter() {
        let dir = temp_dir("basic");
        let prefix = dir.join("db");
        let h1 = sample_header("alpha", "ACC1", 1, 7);
        let h2 = sample_header("beta", "ACC2", 2, 9);
        let mut phr = Vec::new();
        phr.extend_from_slice(&h1);
        phr.extend_from_slice(&h2);
        let psq = [0, 1, 2, 0, 0, 3, 4, 0];
        write_volume_files(
            &prefix,
            2,
            &[0, h1.len() as u32, phr.len() as u32],
            &[0, 4, 8],
            &phr,
            &psq,
        );

        let mut db = BlastDB::new(prefix.to_str().unwrap(), SequenceFileFlags::NONE).unwrap();
        assert_eq!(
            BlastDB::new_with_config(
                prefix.to_str().unwrap(),
                SequenceFileFlags::NONE,
                BlastDbConfig {
                    multiprocessing: true,
                },
            )
            .unwrap_err(),
            "Multiprocessing mode is not compatible with BLAST databases."
        );
        assert_eq!(db.file_count(), 1);
        assert_eq!(db.sequence_count(), 2);
        assert_eq!(db.letters(), 5);
        assert_eq!(db.db_version(), 4);
        assert_eq!(db.program_build_version(), 0);
        assert!(!db.long_seqids());
        assert!(db.volumes().is_empty());

        let info = db.read_seqinfo().unwrap();
        assert_eq!(info, SeqInfo::new(0, 3));
        db.putback_seqinfo();
        assert_eq!(db.tell_seq(), 0);
        assert_eq!(db.id_len(&info, &SeqInfo::default()).unwrap(), h1.len());

        db.init_seq_access();
        let (seq, id) = db.read_seq().unwrap();
        assert_eq!(seq, vec![0, 20]);
        assert_eq!(id, "ACC1.1 alpha");
        assert_eq!(db.tell_seq(), 1);
        assert_eq!(db.seqid(1, true, true).unwrap(), "ACC2.2 beta");

        db.set_dictionary_oids(vec![vec![1]]);
        assert_eq!(db.dict_seq(0, 17).unwrap(), vec![4, 3]);
        assert_eq!(db.dict_seq(1, 0).unwrap_err(), "Dictionary not loaded.");
        db.end_random_access(false);
        assert_eq!(db.dict_seq(0, 0).unwrap(), vec![4, 3]);
        db.end_random_access(true);
        assert_eq!(db.dict_seq(0, 0).unwrap_err(), "Dictionary not loaded.");

        db.set_seqinfo_ptr(0).unwrap();
        let chunk = db
            .raw_chunk(
                10,
                SequenceFileFlags::SEQS
                    | SequenceFileFlags::TITLES
                    | SequenceFileFlags::TAXON_MAPPING
                    | SequenceFileFlags::FULL_TITLES,
            )
            .unwrap();
        assert_eq!(chunk.no, 0);
        assert_eq!(chunk.begin, 0);
        assert_eq!(chunk.end, 2);
        let pkg = chunk
            .decode(
                SequenceFileFlags::SEQS
                    | SequenceFileFlags::TITLES
                    | SequenceFileFlags::TAXON_MAPPING
                    | SequenceFileFlags::FULL_TITLES,
                None,
                None,
            )
            .unwrap();
        assert_eq!(pkg.ids.get(0), b"ACC1.1 alpha");
        assert_eq!(pkg.ids.get(1), b"ACC2.2 beta");
        assert_eq!(pkg.taxids, vec![(0, 7), (1, 9)]);

        let acc_file = dir.join("acc.txt");
        {
            let mut f = std::fs::File::create(&acc_file).unwrap();
            writeln!(f, "ACC2.2").unwrap();
        }
        let filter = db
            .filter_by_accession(acc_file.to_str().unwrap(), false)
            .unwrap();
        assert_eq!(filter.oid_filter, vec![false, true]);
        assert_eq!(filter.letter_count, 3);

        drop(db);
        std::fs::remove_file(prefix.with_extension("pin")).unwrap();
        std::fs::remove_file(prefix.with_extension("phr")).unwrap();
        std::fs::remove_file(prefix.with_extension("psq")).unwrap();
        std::fs::remove_file(acc_file).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn test_blastdb_taxon_names_ranks_and_unsupported_methods() {
        let dir = temp_dir("taxonomy");
        let prefix = dir.join("db");
        let h = sample_header("alpha", "ACC1", 1, 7);
        write_volume_files(&prefix, 1, &[0, h.len() as u32], &[0, 4], &h, &[0, 1, 2, 0]);
        std::fs::write(
            dir.join("names.dmp"),
            "7\t|\tAlpha species\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("nodes.dmp"),
            "7\t|\t1\t|\tspecies\t|\t\n8\t|\t1\t|\tcustom rank\t|\t\n",
        )
        .unwrap();

        let db_path = prefix.to_str().unwrap();
        let db = BlastDB::new(
            db_path,
            SequenceFileFlags::TAXON_SCIENTIFIC_NAMES | SequenceFileFlags::TAXON_RANKS,
        )
        .unwrap();
        assert_eq!(db.taxon_scientific_name(7), "Alpha species");
        assert_eq!(db.taxon_scientific_name(99), "99");
        assert_eq!(db.rank(7), Rank::Species as i32);
        assert_eq!(db.rank(8), Rank::COUNT as i32);
        assert!(db
            .print_info()
            .contains("Taxonomic ids assigned to ranks: 2"));
        assert_eq!(
            BlastDB::new(db_path, SequenceFileFlags::TAXON_NODES)
                .unwrap_err()
                .to_string(),
            format!(
                "Taxonomy database (taxonomy4blast.sqlite3) file not found in path: {}{}taxonomy4blast.sqlite3. Make sure that the database was downloaded correctly.",
                dir.to_string_lossy(),
                PATH_SEPARATOR
            )
        );

        let mut db2 = BlastDB::new(db_path, SequenceFileFlags::NONE).unwrap();
        assert_eq!(
            db2.seek_chunk(&Chunk::default()).unwrap_err(),
            "Operation not supported"
        );
        assert_eq!(
            db2.set_seqinfo_ptr(1).unwrap_err(),
            "Setting seqinfo pointer to non-zero value is not supported in BLAST databases."
        );
        assert_eq!(
            db2.accession_to_oid("missing").unwrap_err(),
            "Accession not found in database: missing"
        );
        assert_eq!(db2.max_taxid(), 0);
        assert_eq!(db2.get_parent(-3), -3);
        assert_eq!(db2.get_parent(7), -1);

        let taxonomy_db = dir.join("taxonomy4blast.sqlite3");
        SqliteConnection::create_test_database(taxonomy_db.to_str().unwrap()).unwrap();
        let mut db3 = BlastDB::new(db_path, SequenceFileFlags::TAXON_NODES).unwrap();
        assert_eq!(db3.max_taxid(), 42);
        assert_eq!(db3.get_parent(42), 7);
        assert_eq!(db3.get_parent(42), 7); // cached path
        assert_eq!(db3.get_parent(7), 3);
        assert_eq!(db3.get_parent(3), 1);
        assert_eq!(db3.get_parent(1), 1);
        assert_eq!(db3.get_parent(8), -1);
        assert!(db3.print_info().contains("Maximum taxid in database: 42"));
        db3.close();

        // Exercise the repository's synthetic BLAST volume and all three
        // taxonomy sidecars together. This is the closest local equivalent of
        // an NCBI BLAST taxonomy database without downloading the external
        // taxonomy4blast distribution.
        let mut combined = BlastDB::new(
            db_path,
            SequenceFileFlags::TAXON_MAPPING
                | SequenceFileFlags::TAXON_NODES
                | SequenceFileFlags::TAXON_RANKS
                | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES,
        )
        .unwrap();
        assert_eq!(combined.taxids(0), vec![7]);
        assert_eq!(combined.get_parent(7), 3);
        assert_eq!(combined.get_parent(3), 1);
        assert_eq!(combined.rank(7), Rank::Species as i32);
        assert_eq!(combined.taxon_scientific_name(7), "Alpha species");
        combined.close();

        drop(combined);
        drop(db3);
        drop(db2);
        drop(db);
        std::fs::remove_file(prefix.with_extension("pin")).unwrap();
        std::fs::remove_file(prefix.with_extension("phr")).unwrap();
        std::fs::remove_file(prefix.with_extension("psq")).unwrap();
        std::fs::remove_file(dir.join("names.dmp")).unwrap();
        std::fs::remove_file(dir.join("nodes.dmp")).unwrap();
        std::fs::remove_file(taxonomy_db).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn test_volume_traversal_seq_data_chunk_order_and_rewind() {
        let dir = temp_dir("volumes");
        let v1 = dir.join("vol1");
        let v2 = dir.join("vol2");
        let h1 = sample_header("first", "ONE", 1, 7);
        let h2 = sample_header("second", "TWO", 1, 9);
        write_volume_files(&v1, 1, &[0, h1.len() as u32], &[0, 4], &h1, &[0, 1, 2, 0]);
        write_volume_files(&v2, 1, &[0, h2.len() as u32], &[0, 4], &h2, &[0, 3, 4, 0]);
        let pal = dir.join("db.pal");
        std::fs::write(&pal, "DBLIST vol1 vol2\n").unwrap();

        let mut db = BlastDB::new(pal.to_str().unwrap(), SequenceFileFlags::NONE).unwrap();
        let mut sequence = Vec::new();
        db.seq_data(1, &mut sequence).unwrap();
        assert_eq!(sequence, vec![4, 3]);

        db.set_seqinfo_ptr(0).unwrap();
        let first = db.raw_chunk(100, SequenceFileFlags::SEQS).unwrap();
        assert_eq!((first.no, first.begin, first.end), (0, 0, 1));
        let second = db.raw_chunk(100, SequenceFileFlags::SEQS).unwrap();
        assert_eq!((second.no, second.begin, second.end), (1, 1, 2));
        assert_eq!(db.raw_chunk_no(), 2);

        db.set_seqinfo_ptr(0).unwrap();
        assert_eq!(db.raw_chunk_no(), 0);
        let rewound = db.raw_chunk(100, SequenceFileFlags::SEQS).unwrap();
        assert_eq!((rewound.no, rewound.begin, rewound.end), (0, 0, 1));

        drop(db);
        for prefix in [&v1, &v2] {
            std::fs::remove_file(prefix.with_extension("pin")).unwrap();
            std::fs::remove_file(prefix.with_extension("phr")).unwrap();
            std::fs::remove_file(prefix.with_extension("psq")).unwrap();
        }
        std::fs::remove_file(pal).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }

    #[test]
    fn source_accession_list_reader_decodes_gzip_and_preserves_byte_lines() {
        let dir = temp_dir("accession-byte-lines");
        let path = dir.join("ids.gz");
        let mut output = flate2::write::GzEncoder::new(
            std::fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        output.write_all(b"ONE.1\r\n\xff\r\n\nTWO.1").unwrap();
        output.finish().unwrap();
        let entries = read_accession_list_source(&path, 64).unwrap();
        assert_eq!(entries.len(), 4);
        for key in [b"ONE.1".as_slice(), b"\xff", b"", b"TWO.1"] {
            assert_eq!(entries.get(key), Some(&false));
        }
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn metadata_accession_filters_use_source_codecs_and_retain_typed_errors() {
        let dir = temp_dir("accession-source-io");
        let prefix = dir.join("db");
        let first = sample_header("", "ONE", 1, 7);
        let mut headers = first.clone();
        headers.extend_from_slice(&sample_header("", "TWO", 1, 9));
        write_volume_files(
            &prefix,
            2,
            &[0, first.len() as u32, headers.len() as u32],
            &[0, 4, 8],
            &headers,
            &[0, 1, 2, 0, 0, 3, 4, 0],
        );
        let mut database = BlastDB::new(prefix.to_str().unwrap(), SequenceFileFlags::NONE).unwrap();
        let ids = dir.join("ids.gz");
        let mut output = flate2::write::GzEncoder::new(
            std::fs::File::create(&ids).unwrap(),
            flate2::Compression::default(),
        );
        output.write_all(b"ONE.1\r\nTWO.1").unwrap();
        output.finish().unwrap();
        let filter = database
            .search_filter_io("", "", ids.to_str().unwrap(), false, 64)
            .unwrap()
            .unwrap();
        assert_eq!(filter.oid_filter, [true, true]);
        assert_eq!(filter.letter_count, 6);

        let missing = dir.join("missing.txt");
        let error = database
            .search_filter_io("", "", missing.to_str().unwrap(), false, 64)
            .unwrap_err();
        assert_eq!(
            crate::util::io::file_open_error::find(&error).unwrap().path,
            missing
        );

        let bytes = dir.join("byte-ids");
        std::fs::write(&bytes, b"\xff\r\n").unwrap();
        let error = database
            .search_filter_io("", "", bytes.to_str().unwrap(), false, 64)
            .unwrap_err();
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error).unwrap(),
            b"Accession not found in database: \xff. Use --skip-missing-seqids to ignore."
        );
        let filter = database
            .search_filter_io("", "", bytes.to_str().unwrap(), true, 64)
            .unwrap()
            .unwrap();
        assert_eq!(filter.oid_filter, [false, false]);
        drop(database);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn early_taxonomy_mapping_deduplicates_deflines_before_gnu_group_order() {
        let dir = temp_dir("early-taxonomy-deflines");
        let prefix = dir.join("db");
        let mut first = vec![0x30, 0];
        for (accession, taxid) in [("A", 9), ("B", 7), ("C", 9)] {
            let header = sample_header("", accession, 1, taxid);
            first.extend_from_slice(&header[2..]);
        }
        assert!(first.len() - 2 < 128);
        first[1] = (first.len() - 2) as u8;
        let mut headers = first.clone();
        headers.extend_from_slice(&sample_header("", "D", 1, 9));
        write_volume_files(
            &prefix,
            2,
            &[0, first.len() as u32, headers.len() as u32],
            &[0, 4, 8],
            &headers,
            &[0, 1, 2, 0, 0, 3, 4, 0],
        );
        let mut database =
            BlastDB::new(prefix.to_str().unwrap(), SequenceFileFlags::TAXON_MAPPING).unwrap();
        // RawChunk's sorted set inserts 7 then 9; GNU unordered_multimap's
        // equal-key iteration returns the most recently inserted taxid first.
        assert_eq!(database.taxids(0), [9, 7]);
        assert_eq!(database.taxids(1), [9]);
        assert!(!database.flags().contains(SequenceFileFlags::TAXON_MAPPING));
        // Explicit additions use the same insertion policy, including actual
        // duplicate pairs: Source only deduplicates within decoded packages.
        database.add_taxid_mapping(&[(0, 5), (1, 4), (0, 9)]);
        assert_eq!(database.taxids(0), [9, 5, 9, 7]);
        assert_eq!(database.taxids(1), [4, 9]);
        drop(database);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn test_early_taxonomy_mapping_is_loaded_and_flag_cleared() {
        let dir = temp_dir("early-taxonomy");
        let prefix = dir.join("db");
        let h1 = sample_header("first", "ONE", 1, 7);
        let h2 = sample_header("second", "TWO", 1, 9);
        let mut headers = h1.clone();
        headers.extend_from_slice(&h2);
        write_volume_files(
            &prefix,
            2,
            &[0, h1.len() as u32, headers.len() as u32],
            &[0, 4, 8],
            &headers,
            &[0, 1, 2, 0, 0, 3, 4, 0],
        );
        let db = BlastDB::new(prefix.to_str().unwrap(), SequenceFileFlags::TAXON_MAPPING).unwrap();
        assert_eq!(db.taxids(0), vec![7]);
        assert_eq!(db.taxids(1), vec![9]);
        assert!(!db.flags().contains(SequenceFileFlags::TAXON_MAPPING));

        drop(db);
        std::fs::remove_file(prefix.with_extension("pin")).unwrap();
        std::fs::remove_file(prefix.with_extension("phr")).unwrap();
        std::fs::remove_file(prefix.with_extension("psq")).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }
}
