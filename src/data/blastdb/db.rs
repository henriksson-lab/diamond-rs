use crate::basic::value::{Letter, Loc, TaxId};
use crate::data::blastdb::taxdmp::{read_names_dmp, read_nodes_dmp};
use crate::data::blastdb::volume::{
    build_title, BlastDefLine, BlastVolume, DbFilter, Pal, RawChunk, SequenceFileFlags,
};
use crate::data::taxonomy::Rank;
use crate::util::system::{absolute_path, exists, PATH_SEPARATOR};
use std::collections::{BTreeMap, HashMap};

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
    file_name: String,
    pal: Pal,
    taxon_mapping: Vec<(u64, TaxId)>,
    custom_ranks: BTreeMap<String, i32>,
    rank_mapping: HashMap<TaxId, i32>,
    oid: u64,
    long_seqids: bool,
    flags: SequenceFileFlags,
    parent_cache: Vec<TaxId>,
    extra_names: HashMap<TaxId, String>,
    volumes: BTreeMap<u64, String>,
    volume: BlastVolume,
    raw_chunk_no: i32,
    seq_length: Vec<Loc>,
}

impl BlastDB {
    pub fn new(file_name: &str, flags: SequenceFileFlags) -> Result<Self, String> {
        let pal = Pal::new(file_name)?;
        let volume =
            BlastVolume::new(&pal.volumes[0], 0, pal.oid_index[0], pal.oid_index[1], true)?;
        let mut db = Self {
            file_name: file_name.to_string(),
            pal,
            taxon_mapping: Vec::new(),
            custom_ranks: BTreeMap::new(),
            rank_mapping: HashMap::new(),
            oid: 0,
            long_seqids: false,
            flags,
            parent_cache: Vec::new(),
            extra_names: HashMap::new(),
            volumes: BTreeMap::new(),
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

        let (dbdir, _dbfile) = absolute_path(file_name);
        if db.flags.contains(SequenceFileFlags::TAXON_RANKS) {
            let file = format!("{dbdir}{PATH_SEPARATOR}nodes.dmp");
            if !exists(&file) {
                return Err(format!(
                    "Taxonomy rank information (nodes.dmp) is missing in search path ({dbdir}). Download and extract this file in the database directory: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"
                ));
            }
            let mut next_rank = Rank::COUNT as i32;
            read_nodes_dmp(&file, |taxid, _parent, rank| {
                if let Some(p) = Rank::predefined(rank) {
                    db.rank_mapping.insert(taxid, p as i32);
                } else if let Some(&r) = db.custom_ranks.get(rank) {
                    db.rank_mapping.insert(taxid, r);
                } else {
                    let r = next_rank;
                    next_rank += 1;
                    db.custom_ranks.insert(rank.to_string(), r);
                    db.rank_mapping.insert(taxid, r);
                }
            })
            .map_err(|e| e.to_string())?;
        }

        if db.flags.contains(SequenceFileFlags::TAXON_NODES) {
            let dbpath = format!("{dbdir}{PATH_SEPARATOR}taxonomy4blast.sqlite3");
            if !exists(&dbpath) {
                return Err(format!(
                    "Taxonomy database (taxonomy4blast.sqlite3) file not found in path: {dbpath}. Make sure that the database was downloaded correctly."
                ));
            }
            return Err(
                "SQLite taxonomy database access is not linked in this Rust crate yet".to_string(),
            );
        }

        if db.flags.contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES) {
            let file = format!("{dbdir}{PATH_SEPARATOR}names.dmp");
            if !exists(&file) {
                return Err(format!(
                    "Taxonomy names information (names.dmp) is missing in search path ({dbdir}). Download and extract this file in the database directory: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"
                ));
            }
            read_names_dmp(&file, |taxid, name| {
                db.extra_names
                    .entry(taxid)
                    .or_insert_with(|| name.to_string());
            })
            .map_err(|e| e.to_string())?;
        }

        if db.flags.contains(SequenceFileFlags::NEED_LENGTH_LOOKUP) {
            db.seq_length.reserve(db.pal.sequence_count as usize);
            for i in 0..db.pal.sequence_count {
                db.open_volume(i)?;
                let volume_oid = (i - db.volume.begin) as u32;
                db.seq_length.push(db.volume.length(volume_oid));
            }
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
        self.volume = BlastVolume::new(
            &self.pal.volumes[idx],
            idx as i32,
            self.pal.oid_index[idx],
            self.pal.oid_index[idx + 1],
            true,
        )?;
        Ok(())
    }

    pub fn print_info(&self) -> String {
        let mut out = format!(
            "Database: {} (type: BLAST database, volumes: {}, sequences: {}, letters: {})\n",
            self.file_name,
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
        out
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
        self.fetch_seqid(oid as u64, all, full_titles)
    }

    pub fn deflines(
        &mut self,
        oid: u64,
        all: bool,
        full_titles: bool,
        taxids: bool,
    ) -> Result<Vec<BlastDefLine>, String> {
        self.open_volume(oid)?;
        let volume_oid = (oid - self.volume.begin) as u32;
        self.volume.deflines(volume_oid, all, full_titles, taxids)
    }

    pub fn skip_id_data(&mut self) {}

    pub fn fetch_seqid(
        &mut self,
        oid: u64,
        all: bool,
        full_titles: bool,
    ) -> Result<String, String> {
        self.open_volume(oid)?;
        let taxids = self.flags.contains(SequenceFileFlags::TAXON_MAPPING);
        let volume_oid = (oid - self.volume.begin) as u32;
        let deflines = self.volume.deflines(volume_oid, all, full_titles, taxids)?;
        if taxids && !self.taxon_mapping.iter().any(|&(o, _)| o == oid) {
            for i in &deflines {
                if let Some(taxid) = i.taxid {
                    self.taxon_mapping.push((oid, taxid));
                }
            }
        }
        Ok(build_title(&deflines, "\x01", all))
    }

    pub fn add_taxid_mapping(&mut self, taxids: &[(u64, TaxId)]) {
        self.taxon_mapping.extend_from_slice(taxids);
    }

    pub fn seqid(&mut self, oid: u64, all: bool, full_titles: bool) -> Result<String, String> {
        self.fetch_seqid(oid, all, full_titles)
    }

    pub fn dict_seq(&mut self, _dict_id: i64, _ref_block: usize) -> Result<Vec<Letter>, String> {
        Err("Dictionary not loaded.".to_string())
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
        self.open_volume(self.oid)?;
        let volume_oid = (self.oid - self.volume.begin) as u32;
        let seq = self.volume.sequence(volume_oid)?;
        let deflines = self.volume.deflines(volume_oid, true, true, false)?;
        let id = build_title(&deflines, " >", true);
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

    pub fn close(&mut self) {}

    pub fn filter_by_accession(
        &mut self,
        file_name: &str,
        skip_missing_seqids: bool,
    ) -> Result<DbFilter, String> {
        let mut v = DbFilter::new(self.sequence_count() as usize);
        let text = std::fs::read_to_string(file_name).map_err(|e| e.to_string())?;
        let mut accs: HashMap<String, bool> =
            text.lines().map(|line| (line.to_string(), false)).collect();
        self.set_seqinfo_ptr(0)?;
        loop {
            let chunk = self.raw_chunk(
                1_000_000_000,
                SequenceFileFlags::SEQS
                    | SequenceFileFlags::TITLES
                    | SequenceFileFlags::FULL_TITLES,
            )?;
            let pkg = chunk.decode(
                SequenceFileFlags::SEQS | SequenceFileFlags::FULL_TITLES,
                None,
                Some(&mut accs),
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

        if !skip_missing_seqids {
            for (acc, found) in &accs {
                if !found {
                    return Err(format!(
                        "Accession not found in database: {acc}. Use --skip-missing-seqids to ignore."
                    ));
                }
            }
        }
        Ok(v)
    }

    pub fn file_name(&self) -> String {
        self.file_name.clone()
    }

    pub fn taxids(&self, oid: usize) -> Vec<TaxId> {
        self.taxon_mapping
            .iter()
            .filter_map(|&(o, t)| if o == oid as u64 { Some(t) } else { None })
            .collect()
    }

    pub fn max_taxid(&self) -> TaxId {
        0
    }

    pub fn get_parent(&mut self, taxid: TaxId) -> TaxId {
        if taxid <= 0 {
            return taxid;
        }
        if taxid as usize >= self.parent_cache.len() {
            return -1;
        }
        self.parent_cache[taxid as usize]
    }

    pub fn taxon_scientific_name(&self, taxid: TaxId) -> String {
        self.extra_names
            .get(&taxid)
            .cloned()
            .unwrap_or_else(|| taxid.to_string())
    }

    pub fn seq_data(&mut self, oid: usize, dst: &mut Vec<Letter>) -> Result<(), String> {
        self.open_volume(oid as u64)?;
        let volume_oid = (self.oid - self.volume.begin) as u32;
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

    pub fn end_random_access(&mut self, _dictionary: bool) {}

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

    pub fn volumes(&self) -> &BTreeMap<u64, String> {
        &self.volumes
    }

    pub fn flags(&self) -> SequenceFileFlags {
        self.flags
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
        let mut seqid_inner = Vec::new();
        seqid_inner.extend_from_slice(&[0xa1, (2 + acc.len()) as u8, 0x1a, acc.len() as u8]);
        seqid_inner.extend_from_slice(acc.as_bytes());
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
        defline.extend_from_slice(title.as_bytes());
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

        std::fs::remove_file(prefix.with_extension("pin")).unwrap();
        std::fs::remove_file(prefix.with_extension("phr")).unwrap();
        std::fs::remove_file(prefix.with_extension("psq")).unwrap();
        std::fs::remove_file(dir.join("names.dmp")).unwrap();
        std::fs::remove_file(dir.join("nodes.dmp")).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }
}
