use crate::basic::value::{Letter, Loc, TaxId, NCBI_TO_STD};
use crate::data::blastdb::asn1::{decode, Node};
use crate::data::blastdb::ber::{decode_integer, read_be32, read_le64, read_pascal_string};
use crate::data::sequence_set::{SequenceSet, StringSet};
use crate::util::io::File as DiamondFile;
use crate::util::string::ends_with;
use crate::util::system::{absolute_path, exists, is_absolute_path, PATH_SEPARATOR};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Read;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SeqId {
    pub type_: String,
    pub value: String,
    pub version: Option<i64>,
    pub chain: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BlastDefLine {
    pub title: String,
    pub seqids: Vec<SeqId>,
    pub taxid: Option<TaxId>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PinIndex {
    pub version: u32,
    pub is_protein: bool,
    pub volume_number: u32,
    pub title: String,
    pub lmdb_file: String,
    pub date: String,
    pub num_oids: u32,
    pub total_length: u64,
    pub max_length: u32,
    pub header_index: Vec<u32>,
    pub sequence_index: Vec<u32>,
    pub ambiguity_offsets_offset: usize,
    pub pin_length: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SequenceFileFlags(pub i32);

impl SequenceFileFlags {
    pub const NONE: Self = Self(0);
    pub const NO_COMPATIBILITY_CHECK: Self = Self(1);
    pub const NO_FASTA: Self = Self(1 << 1);
    pub const ALL_SEQIDS: Self = Self(1 << 2);
    pub const FULL_TITLES: Self = Self(1 << 3);
    pub const TARGET_SEQS: Self = Self(1 << 4);
    pub const SELF_ALN_SCORES: Self = Self(1 << 5);
    pub const NEED_LETTER_COUNT: Self = Self(1 << 6);
    pub const ACC_TO_OID_MAPPING: Self = Self(1 << 7);
    pub const OID_TO_ACC_MAPPING: Self = Self(1 << 8);
    pub const NEED_LENGTH_LOOKUP: Self = Self(1 << 9);
    pub const NEED_EARLY_TAXON_MAPPING: Self = Self(1 << 10);
    pub const TAXON_MAPPING: Self = Self(1 << 11);
    pub const TAXON_NODES: Self = Self(1 << 12);
    pub const TAXON_SCIENTIFIC_NAMES: Self = Self(1 << 13);
    pub const TAXON_RANKS: Self = Self(1 << 14);
    pub const SEQS: Self = Self(1 << 15);
    pub const TITLES: Self = Self(1 << 16);
    pub const QUALITY: Self = Self(1 << 17);
    pub const LAZY_MASKING: Self = Self(1 << 18);
    pub const DNA_PRESERVATION: Self = Self(1 << 19);
    pub const ALL: Self = Self(Self::SEQS.0 | Self::TITLES.0);

    pub fn contains(self, rhs: Self) -> bool {
        (self.0 & rhs.0) != 0
    }
}

impl std::ops::BitOr for SequenceFileFlags {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl std::ops::BitAnd for SequenceFileFlags {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self::Output {
        Self(self.0 & rhs.0)
    }
}

impl std::ops::BitOrAssign for SequenceFileFlags {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

impl std::ops::BitAndAssign for SequenceFileFlags {
    fn bitand_assign(&mut self, rhs: Self) {
        self.0 &= rhs.0;
    }
}

impl std::ops::Not for SequenceFileFlags {
    type Output = Self;

    fn not(self) -> Self::Output {
        Self(!self.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DbFilter {
    pub oid_filter: Vec<bool>,
    pub letter_count: u64,
}

impl DbFilter {
    pub fn new(size: usize) -> Self {
        Self {
            oid_filter: vec![false; size],
            letter_count: 0,
        }
    }

    pub fn get(&self, oid: u64) -> bool {
        self.oid_filter.get(oid as usize).copied().unwrap_or(false)
    }
}

pub struct DecodedPackage {
    pub ids: StringSet,
    pub seqs: SequenceSet,
    pub oids: Vec<u64>,
    pub taxids: Vec<(u64, TaxId)>,
    pub no: i32,
}

impl Default for DecodedPackage {
    fn default() -> Self {
        Self {
            ids: StringSet::new(),
            seqs: SequenceSet::new(),
            oids: Vec::new(),
            taxids: Vec::new(),
            no: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RawChunk {
    pub no: i32,
    pub seq_data: Vec<u8>,
    pub phr_data: Vec<u8>,
    pub seq_index: Vec<u32>,
    pub phr_index: Vec<u32>,
    pub begin: u64,
    pub end: u64,
    pub letters: usize,
}

impl RawChunk {
    pub fn empty(&self) -> bool {
        self.end <= self.begin
    }

    pub fn decode(
        &self,
        flags: SequenceFileFlags,
        filter: Option<&DbFilter>,
        mut accs: Option<&mut HashMap<String, bool>>,
    ) -> Result<DecodedPackage, String> {
        assert!(filter.is_none() || accs.is_none());
        let mut pkg = DecodedPackage {
            no: self.no,
            ..Default::default()
        };
        let n = (self.end - self.begin) as usize;
        let mut seq_ptr = 0usize;
        let mut phr_ptr = 0usize;
        let titles = flags.contains(SequenceFileFlags::TITLES);
        let seqs = flags.contains(SequenceFileFlags::SEQS);
        let taxids = flags.contains(SequenceFileFlags::TAXON_MAPPING);
        let full_titles = flags.contains(SequenceFileFlags::FULL_TITLES);
        let all_seqids = flags.contains(SequenceFileFlags::ALL_SEQIDS);
        pkg.oids.reserve(n);

        for i in 0..n {
            let oid = self.begin + i as u64;
            let mut f = filter.map(|filter| filter.get(oid)).unwrap_or(true);

            if titles || taxids || accs.is_some() {
                let lhdr = (self.phr_index[i + 1] - self.phr_index[i]) as usize;
                if f || accs.is_some() {
                    let deflines = decode_deflines(
                        &self.phr_data[phr_ptr..phr_ptr + lhdr],
                        all_seqids,
                        full_titles,
                        taxids,
                    )?;
                    if let Some(accs) = accs.as_deref_mut() {
                        f = acc_filter(&deflines, accs);
                    }
                    if f && titles {
                        let title = build_title(&deflines, "\x01", true);
                        pkg.ids.push_back(title.as_bytes());
                    }
                    if f && taxids {
                        let mut s = BTreeSet::new();
                        for j in &deflines {
                            if let Some(taxid) = j.taxid {
                                s.insert(taxid);
                            }
                        }
                        for t in s {
                            pkg.taxids.push((oid, t));
                        }
                    }
                }
                phr_ptr += lhdr;
            }

            if seqs {
                let lseq = (self.seq_index[i + 1] - self.seq_index[i]) as usize;
                if f {
                    let seq = decode_protein_sequence(&self.seq_data[seq_ptr..seq_ptr + lseq])?;
                    pkg.seqs.push(&seq);
                }
                seq_ptr += lseq;
            }

            if f {
                pkg.oids.push(oid);
            }
        }
        Ok(pkg)
    }
}

#[derive(Debug)]
pub struct BlastVolume {
    pub idx: i32,
    pub begin: u64,
    pub end: u64,
    index: PinIndex,
    phr_mapping: DiamondFile,
    psq_mapping: DiamondFile,
    seq_ptr: u32,
    hdr_ptr: u32,
}

impl BlastVolume {
    pub fn new(
        path: &str,
        idx: i32,
        begin: u64,
        end: u64,
        load_index: bool,
    ) -> Result<Self, String> {
        let phr_mapping =
            DiamondFile::open(&format!("{path}.phr"), "rb").map_err(|e| e.to_string())?;
        let psq_mapping =
            DiamondFile::open(&format!("{path}.psq"), "rb").map_err(|e| e.to_string())?;
        let mut pin = DiamondFile::open(&format!("{path}.pin"), "rb").map_err(|e| e.to_string())?;
        let index = Self::parse_pin_file(&mut pin, load_index)?;
        pin.close().map_err(|e| e.to_string())?;
        Ok(Self {
            idx,
            begin,
            end,
            index,
            phr_mapping,
            psq_mapping,
            seq_ptr: 0,
            hdr_ptr: 0,
        })
    }

    pub fn from_parts(
        index: PinIndex,
        phr_mapping: DiamondFile,
        psq_mapping: DiamondFile,
        idx: i32,
        begin: u64,
        end: u64,
    ) -> Self {
        Self {
            idx,
            begin,
            end,
            index,
            phr_mapping,
            psq_mapping,
            seq_ptr: 0,
            hdr_ptr: 0,
        }
    }

    pub fn index(&self) -> &PinIndex {
        &self.index
    }

    pub fn seq_ptr(&self) -> u32 {
        self.seq_ptr
    }

    pub fn parse_pin_file(mapping: &mut DiamondFile, load_index: bool) -> Result<PinIndex, String> {
        let mut index = PinIndex {
            version: read_be32(mapping).map_err(|e| e.to_string())?,
            ..Default::default()
        };
        if index.version != 4 && index.version != 5 {
            return Err(format!(
                "Unsupported database format version: {}",
                index.version
            ));
        }

        let seq_type_flag = read_be32(mapping).map_err(|e| e.to_string())?;
        index.is_protein = seq_type_flag == 1;

        if index.version == 5 {
            index.volume_number = read_be32(mapping).map_err(|e| e.to_string())?;
        }
        index.title = read_pascal_string(mapping).map_err(|e| e.to_string())?;

        if index.version == 5 {
            index.lmdb_file = read_pascal_string(mapping).map_err(|e| e.to_string())?;
        }

        index.date = read_pascal_string(mapping).map_err(|e| e.to_string())?;
        index.num_oids = read_be32(mapping).map_err(|e| e.to_string())?;
        index.total_length = read_le64(mapping).map_err(|e| e.to_string())?;
        index.max_length = read_be32(mapping).map_err(|e| e.to_string())?;
        if !load_index {
            return Ok(index);
        }

        let count = index.num_oids as usize + 1;
        index.header_index.reserve(count);
        index.sequence_index.reserve(count);
        for _ in 0..count {
            index
                .header_index
                .push(read_be32(mapping).map_err(|e| e.to_string())?);
        }

        for _ in 0..count {
            index
                .sequence_index
                .push(read_be32(mapping).map_err(|e| e.to_string())?);
        }

        Ok(index)
    }

    pub fn deflines(
        &mut self,
        oid: u32,
        all: bool,
        full_titles: bool,
        taxids: bool,
    ) -> Result<Vec<BlastDefLine>, String> {
        if oid >= self.index.num_oids {
            return Err("OID exceeds number of sequences in volume".to_string());
        }
        let header_offset = self.index.header_index[oid as usize] as usize;
        let next_header_offset = self.index.header_index[oid as usize + 1] as usize;
        if next_header_offset < header_offset {
            return Err("Header offsets exceed PHR file size".to_string());
        }
        let header_length = next_header_offset - header_offset;
        if oid != self.hdr_ptr {
            self.phr_mapping
                .seek(header_offset as i64, std::io::SeekFrom::Start(0))
                .map_err(|e| e.to_string())?;
        }
        self.hdr_ptr = oid + 1;
        let data = self
            .phr_mapping
            .read(header_length)
            .map_err(|e| e.to_string())?;
        decode_deflines(data, all, full_titles, taxids)
    }

    pub fn sequence(&mut self, oid: u32) -> Result<Vec<Letter>, String> {
        if oid >= self.index.num_oids {
            return Err("OID exceeds number of sequences in volume".to_string());
        }

        let start = self.index.sequence_index[oid as usize];
        let end = self.index.sequence_index[oid as usize + 1];

        if !self.index.is_protein {
            return Err("Nucleotide sequence decoding is not supported yet".to_string());
        }
        if oid != self.seq_ptr {
            self.psq_mapping
                .seek(start as i64, std::io::SeekFrom::Start(0))
                .map_err(|e| e.to_string())?;
        }
        self.seq_ptr = oid + 1;
        let data = self
            .psq_mapping
            .read((end - start) as usize)
            .map_err(|e| e.to_string())?;
        decode_protein_sequence(data)
    }

    pub fn raw_sequence(&mut self, count: u32) -> Result<Vec<u8>, String> {
        let n = (self.index.sequence_index[(self.seq_ptr + count) as usize]
            - self.index.sequence_index[self.seq_ptr as usize]) as usize;
        let mut v = vec![0u8; n];
        self.psq_mapping
            .read_exact(&mut v)
            .map_err(|e| e.to_string())?;
        self.seq_ptr += count;
        Ok(v)
    }

    pub fn raw_deflines(&mut self, count: u32) -> Result<Vec<u8>, String> {
        let n = (self.index.header_index[(self.hdr_ptr + count) as usize]
            - self.index.header_index[self.hdr_ptr as usize]) as usize;
        let mut v = vec![0u8; n];
        self.phr_mapping
            .read_exact(&mut v)
            .map_err(|e| e.to_string())?;
        self.hdr_ptr += count;
        Ok(v)
    }

    pub fn length(&self, oid: u32) -> Loc {
        length(&self.index.sequence_index, oid as usize)
    }

    pub fn id_len(&self, oid: u32) -> usize {
        id_len(&self.index.header_index, oid as usize)
    }

    pub fn raw_chunk(
        &mut self,
        letters: usize,
        flags: SequenceFileFlags,
    ) -> Result<RawChunk, String> {
        let mut begin = self.hdr_ptr;
        if !flags.contains(SequenceFileFlags::SEQS) {
            if self.seq_ptr != 0 {
                return Err("Volume::raw_chunk".to_string());
            }
        } else if !flags.contains(SequenceFileFlags::TITLES)
            && !flags.contains(SequenceFileFlags::TAXON_MAPPING)
        {
            if self.hdr_ptr != 0 {
                return Err("Volume::raw_chunk".to_string());
            }
            begin = self.seq_ptr;
        } else if self.hdr_ptr != self.seq_ptr {
            return Err(
                "Cannot read raw chunk: last accessed header and sequence OIDs do not match"
                    .to_string(),
            );
        }

        let mut end = begin;
        let mut l = 0usize;
        while end < self.index.num_oids && l < letters {
            l += self.length(end) as usize;
            end += 1;
        }
        let mut chunk = RawChunk {
            letters: l,
            begin: begin as u64 + self.begin,
            end: end as u64 + self.begin,
            ..Default::default()
        };
        let n = end - begin;
        if n == 0 {
            return Ok(chunk);
        }
        if flags.contains(SequenceFileFlags::TITLES)
            || flags.contains(SequenceFileFlags::TAXON_MAPPING)
        {
            chunk.phr_index = self.index.header_index
                [self.hdr_ptr as usize..self.hdr_ptr as usize + n as usize + 1]
                .to_vec();
            chunk.phr_data = self.raw_deflines(n)?;
        }
        if flags.contains(SequenceFileFlags::SEQS) {
            chunk.seq_index = self.index.sequence_index
                [self.seq_ptr as usize..self.seq_ptr as usize + n as usize + 1]
                .to_vec();
            chunk.seq_data = self.raw_sequence(n)?;
        }
        Ok(chunk)
    }

    pub fn rewind(&mut self) -> Result<(), String> {
        self.hdr_ptr = 0;
        self.seq_ptr = 0;
        self.phr_mapping
            .seek(0, std::io::SeekFrom::Start(0))
            .map_err(|e| e.to_string())?;
        self.psq_mapping
            .seek(0, std::io::SeekFrom::Start(0))
            .map_err(|e| e.to_string())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Pal {
    pub volumes: Vec<String>,
    pub metadata: BTreeMap<String, String>,
    pub oid_index: Vec<u64>,
    pub sequence_count: u64,
    pub letters: u64,
    pub version: i32,
}

impl Pal {
    pub fn new(path: &str) -> Result<Self, String> {
        let supported_keys: BTreeSet<&str> = [
            "TITLE",
            "MEMB_BIT",
            "SEQIDLIST",
            "NSEQ",
            "LENGTH",
            "TAXIDLIST",
        ]
        .into_iter()
        .collect();
        let (db_dir, file) = absolute_path(path);
        let mut pal = Pal::default();
        if !exists(&format!("{path}.pal")) && !ends_with(path, ".pal") {
            pal.volumes.push(format!("{db_dir}{PATH_SEPARATOR}{file}"));
        } else {
            let pal_path = if ends_with(path, ".pal") {
                path.to_string()
            } else {
                format!("{path}.pal")
            };
            let mut text = String::new();
            std::fs::File::open(&pal_path)
                .map_err(|_| format!("Unable to open PAL file: {pal_path}"))?
                .read_to_string(&mut text)
                .map_err(|e| e.to_string())?;
            for (line_number0, mut line) in text.lines().map(str::to_string).enumerate() {
                let line_number = line_number0 + 1;
                if let Some(comment) = line.find('#') {
                    line.truncate(comment);
                }
                line = trim(&line);
                if line.is_empty() {
                    continue;
                }

                let key_end = line.find([' ', '\t']).ok_or_else(|| {
                    format!("Error parsing PAL file: line {line_number} is missing a value: {line}")
                })?;
                let key = line[..key_end].to_string();
                let value = trim(&line[key_end + 1..]);
                if value.is_empty() {
                    return Err(format!(
                        "Error parsing PAL file: line {line_number} has an empty value: {line}"
                    ));
                }

                if key == "DBLIST" {
                    let vls = split_whitespace(&value);
                    if vls.is_empty() {
                        return Err(format!(
                            "Error parsing PAL file: DBLIST on line {line_number} does not list any volumes"
                        ));
                    }
                    pal.volumes.extend(vls);
                    for s in &mut pal.volumes {
                        if !is_absolute_path(s) && !s.is_empty() && !s.starts_with('"') {
                            *s = format!("{db_dir}{PATH_SEPARATOR}{s}");
                        }
                    }
                    continue;
                }

                if !supported_keys.contains(key.as_str()) {
                    return Err(format!(
                        "Error parsing PAL file: Unsupported PAL key '{key}' on line {line_number}"
                    ));
                }

                if pal.metadata.contains_key(&key) {
                    return Err(format!(
                        "Error parsing PAL file: Duplicate key '{key}' on line {line_number}"
                    ));
                }

                pal.metadata.insert(key, value);
            }
        }
        pal.sequence_count = 0;
        pal.letters = 0;
        pal.oid_index.push(0);

        let mut it = 0usize;
        while it < pal.volumes.len() {
            let volume = pal.volumes[it].clone();
            if volume.len() >= 2 && volume.starts_with('"') && volume.ends_with('"') {
                let nested = volume[1..volume.len() - 1].to_string();
                pal.volumes.remove(it);
                let nested_path = if is_absolute_path(&nested) {
                    nested
                } else {
                    format!("{db_dir}{PATH_SEPARATOR}{nested}")
                };
                let inserted = pal.recurse(&nested_path, it)?;
                it += inserted;
            } else {
                let vol = BlastVolume::new(&volume, 0, 0, 0, false)?;
                pal.sequence_count += u64::from(vol.index().num_oids);
                pal.oid_index
                    .push(u64::from(vol.index().num_oids) + *pal.oid_index.last().unwrap());
                pal.letters += vol.index().total_length;
                pal.version = vol.index().version as i32;
                it += 1;
            }
        }
        if let Some(seqidlist) = pal.metadata.get_mut("SEQIDLIST") {
            if ends_with(seqidlist, ".bsl") {
                return Err(format!(
                    "Binary SEQIDLIST files(.bsl) are not supported, use text file instead : {seqidlist}"
                ));
            }
            if !is_absolute_path(seqidlist) {
                *seqidlist = format!("{db_dir}{PATH_SEPARATOR}{seqidlist}");
            }
        }
        if let Some(taxidlist) = pal.metadata.get_mut("TAXIDLIST") {
            if !is_absolute_path(taxidlist) {
                *taxidlist = format!("{db_dir}{PATH_SEPARATOR}{taxidlist}");
            }
        }
        assert!(pal.sequence_count > 0);
        Ok(pal)
    }

    pub fn recurse(&mut self, path: &str, volume_it: usize) -> Result<usize, String> {
        let pal = Pal::new(path)?;
        let inserted = pal.volumes.len();
        self.volumes.splice(volume_it..volume_it, pal.volumes);

        for (key, value) in pal.metadata {
            if self.metadata.contains_key(&key) {
                if key == "TITLE" || key == "NSEQ" || key == "LENGTH" {
                    continue;
                }
                return Err(format!("Duplicate key '{key}' in nested PAL file: {path}"));
            } else {
                self.metadata.insert(key, value);
            }
        }
        let base = *self.oid_index.last().unwrap();
        for oid in pal.oid_index.iter().skip(1) {
            self.oid_index.push(*oid + base);
        }
        self.sequence_count += pal.sequence_count;
        self.letters += pal.letters;
        self.version = pal.version;
        Ok(inserted)
    }

    pub fn volume(&self, oid: u64) -> i32 {
        assert!(oid < self.sequence_count);
        let it = self.oid_index.partition_point(|&x| x <= oid);
        it as i32 - 1
    }
}

pub fn tag_name_from_number(num: u32) -> String {
    match num {
        0 => "local",
        1 => "gibbsq",
        2 => "gibbmt",
        3 => "giim",
        4 => "genbank",
        5 => "embl",
        6 => "pir",
        7 => "swissprot",
        8 => "patent",
        9 => "other",
        10 => "general",
        11 => "gi",
        12 => "ddbj",
        13 => "prf",
        14 => "pdb",
        15 => "tpg",
        16 => "tpe",
        17 => "tpd",
        18 => "gpipe",
        19 => "named-annot-track",
        _ => return format!("unknown-{num}"),
    }
    .to_string()
}

fn trim(text: &str) -> String {
    text.trim_matches([' ', '\t', '\r', '\n']).to_string()
}

fn split_whitespace(text: &str) -> Vec<String> {
    text.split_whitespace().map(str::to_string).collect()
}

fn acc_filter(deflines: &[BlastDefLine], accs: &mut HashMap<String, bool>) -> bool {
    for d in deflines {
        for s in &d.seqids {
            if accs.contains_key(&s.value) {
                accs.insert(s.value.clone(), true);
                return true;
            }
            if (s.version.is_some() || s.chain.is_some()) && accs.contains_key(&format_seqid(s)) {
                accs.insert(format_seqid(s), true);
                return true;
            }
        }
    }
    false
}

fn decode_seqid_into(node: &Node, seqid: &mut SeqId) {
    for n4 in &node.children {
        match n4.tag.tag_number {
            1 => {
                for n5 in &n4.children {
                    if n5.tag.tag_number == 26 {
                        seqid.value = String::from_utf8_lossy(&n5.value).into_owned();
                    }
                }
            }
            3 => {
                for n5 in &n4.children {
                    if n5.tag.tag_number == 2 {
                        seqid.version = Some(decode_integer(&n5.value));
                    }
                }
            }
            _ => {}
        }
    }
}

fn decode_seqid(node: &Node) -> SeqId {
    let mut seqid = SeqId::default();
    for n1 in &node.children {
        if n1.tag.tag_number == 16 {
            for n2 in &n1.children {
                match n2.tag.tag_number {
                    0 | 1 | 4 | 5 | 7 | 9 | 12 | 15 | 16 => {
                        seqid.type_ = tag_name_from_number(n2.tag.tag_number);
                        decode_seqid_into(n2, &mut seqid);
                        for n3 in &n2.children {
                            if n3.tag.tag_number == 16 {
                                decode_seqid_into(n3, &mut seqid);
                            }
                        }
                    }
                    14 => {
                        seqid.type_ = tag_name_from_number(n2.tag.tag_number);
                        for n3 in &n2.children {
                            if n3.tag.tag_number == 16 {
                                for n4 in &n3.children {
                                    match n4.tag.tag_number {
                                        0 => {
                                            for n5 in &n4.children {
                                                if n5.tag.tag_number == 26 {
                                                    seqid.value =
                                                        String::from_utf8_lossy(&n5.value)
                                                            .into_owned();
                                                }
                                            }
                                        }
                                        3 => {
                                            for n5 in &n4.children {
                                                if n5.tag.tag_number == 26 {
                                                    seqid.chain = Some(
                                                        String::from_utf8_lossy(&n5.value)
                                                            .into_owned(),
                                                    );
                                                }
                                            }
                                        }
                                        _ => {}
                                    }
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
    }
    seqid
}

fn decode_defline(node: &Node, full_titles: bool, taxids: bool) -> BlastDefLine {
    let mut defline = BlastDefLine::default();
    for n1 in &node.children {
        match n1.tag.tag_number {
            0 => {
                if full_titles {
                    for n2 in &n1.children {
                        if n2.tag.tag_number == 26 {
                            defline.title = String::from_utf8_lossy(&n2.value).into_owned();
                        }
                    }
                }
            }
            1 => defline.seqids.push(decode_seqid(n1)),
            2 => {
                if taxids {
                    for n2 in &n1.children {
                        if n2.tag.tag_number == 2 {
                            defline.taxid = Some(decode_integer(&n2.value) as TaxId);
                        }
                    }
                }
            }
            _ => {}
        }
    }
    defline
}

pub fn decode_deflines(
    header_data: &[u8],
    all: bool,
    full_titles: bool,
    taxids: bool,
) -> Result<Vec<BlastDefLine>, String> {
    let mut out = Vec::new();
    let nodes = decode(header_data).map_err(|e| e.to_string())?;
    if nodes.is_empty() {
        return Ok(out);
    }
    for i in &nodes[0].children {
        out.push(decode_defline(i, full_titles, taxids));
        if !all && !taxids {
            break;
        }
    }
    Ok(out)
}

pub fn format_seqid(id: &SeqId) -> String {
    if id.value.is_empty() {
        return "N/A".to_string();
    }
    let mut os = id.value.clone();
    if let Some(version) = id.version {
        os.push('.');
        os.push_str(&version.to_string());
    }
    if let Some(chain) = &id.chain {
        if !chain.is_empty() {
            os.push('_');
            os.push_str(chain);
        }
    }
    os
}

pub fn build_title(deflines: &[BlastDefLine], delimiter: &str, all: bool) -> String {
    let mut h = String::new();
    for (i, defline) in deflines.iter().enumerate() {
        if i != 0 {
            if !all {
                break;
            }
            h.push_str(delimiter);
        }
        h.push_str(&format_seqid(
            defline.seqids.first().unwrap_or(&SeqId::default()),
        ));
        h.push(' ');
        h.push_str(&defline.title);
    }
    if h.is_empty() {
        h = "N/A".to_string();
    }
    h
}

pub fn decode_protein_sequence(data: &[u8]) -> Result<Vec<Letter>, String> {
    let mut decoded = Vec::with_capacity(data.len());
    for (i, &aa) in data.iter().enumerate() {
        if aa == b'\0' {
            if i == 0 {
                continue;
            } else if i == data.len() - 1 {
                break;
            } else {
                return Err("Unexpected null terminator in sequence data".to_string());
            }
        }
        if usize::from(aa) >= NCBI_TO_STD.len() {
            return Err("Invalid amino acid code in sequence data".to_string());
        }
        decoded.push(NCBI_TO_STD[aa as usize]);
    }
    Ok(decoded)
}

pub fn length(sequence_index: &[u32], oid: usize) -> Loc {
    (sequence_index[oid + 1] - sequence_index[oid] - 1) as Loc
}

pub fn id_len(header_index: &[u32], oid: usize) -> usize {
    (header_index[oid + 1] - header_index[oid]) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_seqid_and_build_title() {
        let seqid = SeqId {
            value: "XP_001".to_string(),
            version: Some(7),
            chain: Some("A".to_string()),
            ..Default::default()
        };
        assert_eq!(format_seqid(&seqid), "XP_001.7_A");
        assert_eq!(format_seqid(&SeqId::default()), "N/A");

        let deflines = vec![
            BlastDefLine {
                title: "first protein".to_string(),
                seqids: vec![seqid],
                taxid: None,
            },
            BlastDefLine {
                title: "second protein".to_string(),
                seqids: vec![SeqId {
                    value: "YP_002".to_string(),
                    ..Default::default()
                }],
                taxid: None,
            },
        ];
        assert_eq!(
            build_title(&deflines, " | ", true),
            "XP_001.7_A first protein | YP_002 second protein"
        );
        assert_eq!(
            build_title(&deflines, " | ", false),
            "XP_001.7_A first protein"
        );
        assert_eq!(build_title(&[], " | ", true), "N/A");
    }

    #[test]
    fn test_decode_deflines_from_ber_tree() {
        let header = [
            0x30, 0x26, 0x30, 0x24, 0xa0, 0x07, 0x1a, 0x05, b't', b'i', b't', b'l', b'e', 0xa1,
            0x13, 0x30, 0x11, 0xa4, 0x0f, 0xa1, 0x07, 0x1a, 0x05, b'A', b'C', b'C', b'1', b'2',
            0xa3, 0x04, 0x02, 0x02, 0x01, 0x02, 0xa2, 0x04, 0x02, 0x02, 0x03, 0xe8,
        ];
        let deflines = decode_deflines(&header, false, true, true).unwrap();
        assert_eq!(deflines.len(), 1);
        assert_eq!(deflines[0].title, "title");
        assert_eq!(deflines[0].taxid, Some(1000));
        assert_eq!(deflines[0].seqids[0].type_, "genbank");
        assert_eq!(deflines[0].seqids[0].value, "ACC12");
        assert_eq!(deflines[0].seqids[0].version, Some(258));
    }

    #[test]
    fn test_decode_protein_sequence_cpp_null_rules() {
        assert_eq!(decode_protein_sequence(&[0, 1, 2, 0]).unwrap(), vec![0, 20]);
        assert_eq!(
            decode_protein_sequence(&[1, 0, 2]).unwrap_err(),
            "Unexpected null terminator in sequence data"
        );
        assert_eq!(
            decode_protein_sequence(&[28]).unwrap_err(),
            "Invalid amino acid code in sequence data"
        );
    }

    #[test]
    fn test_index_span_helpers() {
        assert_eq!(length(&[0, 6, 9], 0), 5);
        assert_eq!(id_len(&[10, 25, 40], 1), 15);
        assert_eq!(tag_name_from_number(19), "named-annot-track");
        assert_eq!(tag_name_from_number(77), "unknown-77");
    }

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
        version: u32,
        num_oids: u32,
        header_index: &[u32],
        sequence_index: &[u32],
        phr: &[u8],
        psq: &[u8],
    ) {
        let mut pin = Vec::new();
        write_be32(&mut pin, version);
        write_be32(&mut pin, 1);
        if version == 5 {
            write_be32(&mut pin, 9);
        }
        write_pascal(&mut pin, "test title");
        if version == 5 {
            write_pascal(&mut pin, "lmdb");
        }
        write_pascal(&mut pin, "2026-05-14");
        write_be32(&mut pin, num_oids);
        write_le64(&mut pin, 1234);
        write_be32(&mut pin, 42);
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

    fn temp_prefix(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-{name}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[test]
    fn test_parse_pin_file_version4_and_5() {
        let prefix = temp_prefix("pin");
        write_volume_files(&prefix, 4, 2, &[0, 10, 25], &[0, 4, 9], b"", b"");
        let mut pin =
            DiamondFile::open(prefix.with_extension("pin").to_str().unwrap(), "rb").unwrap();
        let index = BlastVolume::parse_pin_file(&mut pin, true).unwrap();
        assert_eq!(index.version, 4);
        assert!(index.is_protein);
        assert_eq!(index.title, "test title");
        assert_eq!(index.date, "2026-05-14");
        assert_eq!(index.num_oids, 2);
        assert_eq!(index.total_length, 1234);
        assert_eq!(index.max_length, 42);
        assert_eq!(index.header_index, vec![0, 10, 25]);
        assert_eq!(index.sequence_index, vec![0, 4, 9]);

        let prefix5 = temp_prefix("pin5");
        write_volume_files(&prefix5, 5, 1, &[0, 0], &[0, 0], b"", b"");
        let mut pin5 =
            DiamondFile::open(prefix5.with_extension("pin").to_str().unwrap(), "rb").unwrap();
        let index5 = BlastVolume::parse_pin_file(&mut pin5, false).unwrap();
        assert_eq!(index5.version, 5);
        assert_eq!(index5.volume_number, 9);
        assert_eq!(index5.lmdb_file, "lmdb");
        assert!(index5.header_index.is_empty());

        std::fs::remove_file(prefix.with_extension("pin")).unwrap();
        std::fs::remove_file(prefix.with_extension("phr")).unwrap();
        std::fs::remove_file(prefix.with_extension("psq")).unwrap();
        std::fs::remove_file(prefix5.with_extension("pin")).unwrap();
        std::fs::remove_file(prefix5.with_extension("phr")).unwrap();
        std::fs::remove_file(prefix5.with_extension("psq")).unwrap();
    }

    #[test]
    fn test_blast_volume_raw_chunk_decode_and_accession_filter() {
        let h1 = sample_header("alpha title", "ACC1", 1, 7);
        let h2 = sample_header("beta title", "ACC2", 2, 7);
        let mut phr = Vec::new();
        phr.extend_from_slice(&h1);
        phr.extend_from_slice(&h2);
        let psq = [0, 1, 2, 0, 0, 3, 4, 0];
        let prefix = temp_prefix("volume");
        write_volume_files(
            &prefix,
            4,
            2,
            &[0, h1.len() as u32, phr.len() as u32],
            &[0, 4, 8],
            &phr,
            &psq,
        );

        let mut volume = BlastVolume::new(prefix.to_str().unwrap(), 0, 10, 12, true).unwrap();
        assert_eq!(volume.id_len(0), h1.len());
        assert_eq!(volume.length(0), 3);
        assert_eq!(volume.sequence(1).unwrap(), vec![4, 3]);
        volume.rewind().unwrap();
        let chunk = volume
            .raw_chunk(
                10,
                SequenceFileFlags::SEQS
                    | SequenceFileFlags::TITLES
                    | SequenceFileFlags::TAXON_MAPPING
                    | SequenceFileFlags::FULL_TITLES,
            )
            .unwrap();
        assert_eq!(chunk.begin, 10);
        assert_eq!(chunk.end, 12);
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
        assert_eq!(pkg.oids, vec![10, 11]);
        assert_eq!(pkg.ids.get(0), b"ACC1.1 alpha title");
        assert_eq!(pkg.ids.get(1), b"ACC2.2 beta title");
        assert_eq!(pkg.seqs.get(0), &[0, 20]);
        assert_eq!(pkg.seqs.get(1), &[4, 3]);
        assert_eq!(pkg.taxids, vec![(10, 7), (11, 7)]);

        let mut accs = HashMap::from([("ACC2.2".to_string(), false)]);
        let filtered = chunk
            .decode(
                SequenceFileFlags::SEQS | SequenceFileFlags::FULL_TITLES,
                None,
                Some(&mut accs),
            )
            .unwrap();
        assert_eq!(filtered.oids, vec![11]);
        assert!(accs["ACC2.2"]);

        std::fs::remove_file(prefix.with_extension("pin")).unwrap();
        std::fs::remove_file(prefix.with_extension("phr")).unwrap();
        std::fs::remove_file(prefix.with_extension("psq")).unwrap();
    }

    #[test]
    fn test_pal_parse_volume_lookup_and_metadata_paths() {
        let dir = std::env::temp_dir().join(format!(
            "diamond-rs-pal-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&dir).unwrap();
        let v1 = dir.join("vol1");
        let v2 = dir.join("vol2");
        write_volume_files(&v1, 4, 2, &[0, 0, 0], &[0, 4, 8], b"", &[0; 8]);
        write_volume_files(&v2, 4, 3, &[0, 0, 0, 0], &[0, 4, 8, 12], b"", &[0; 12]);
        let pal_path = dir.join("db.pal");
        std::fs::write(
            &pal_path,
            "TITLE Example DB\nDBLIST vol1 vol2\nSEQIDLIST ids.txt\nTAXIDLIST taxids.txt\n",
        )
        .unwrap();

        let pal = Pal::new(pal_path.to_str().unwrap()).unwrap();
        assert_eq!(pal.sequence_count, 5);
        assert_eq!(pal.letters, 2468);
        assert_eq!(pal.oid_index, vec![0, 2, 5]);
        assert_eq!(pal.volume(0), 0);
        assert_eq!(pal.volume(2), 1);
        assert!(pal.metadata["SEQIDLIST"].ends_with("ids.txt"));
        assert!(pal.metadata["TAXIDLIST"].ends_with("taxids.txt"));

        std::fs::remove_file(v1.with_extension("pin")).unwrap();
        std::fs::remove_file(v1.with_extension("phr")).unwrap();
        std::fs::remove_file(v1.with_extension("psq")).unwrap();
        std::fs::remove_file(v2.with_extension("pin")).unwrap();
        std::fs::remove_file(v2.with_extension("phr")).unwrap();
        std::fs::remove_file(v2.with_extension("psq")).unwrap();
        std::fs::remove_file(pal_path).unwrap();
        std::fs::remove_dir(dir).unwrap();
    }
}
