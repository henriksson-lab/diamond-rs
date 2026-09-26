use crate::basic::value::{Letter, Loc, TaxId};
use crate::data::sequence_set::{SequenceSet, StringSet};
use crate::util::io::File as DiamondFile;
#[cfg(test)]
use std::collections::HashMap;

pub use super::pal::Pal;
pub use super::phr::{
    build_title, decode_deflines, format_seqid, id_len, tag_name_from_number, BlastDefLine, SeqId,
};
pub use super::pin::PinIndex;
pub use super::psq::{decode_protein_sequence, length};

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
}

#[derive(Debug)]
pub struct BlastVolume {
    pub idx: i32,
    pub begin: u64,
    pub end: u64,
    pub(super) index: PinIndex,
    pub(super) phr_mapping: DiamondFile,
    pub(super) psq_mapping: DiamondFile,
    pub(super) seq_ptr: u32,
    pub(super) hdr_ptr: u32,
}

impl BlastVolume {
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

    pub fn length(&self, oid: u32) -> Loc {
        length(&self.index.sequence_index, oid as usize)
    }
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
    fn test_deflines_reject_oid_and_descending_offsets() {
        let prefix = temp_prefix("phr-errors");
        write_volume_files(&prefix, 4, 1, &[3, 2], &[0, 1], b"abc", &[0]);
        let mut volume = BlastVolume::new(prefix.to_str().unwrap(), 0, 0, 1, true).unwrap();
        assert_eq!(
            volume.deflines(1, true, true, true).unwrap_err(),
            "OID exceeds number of sequences in volume"
        );
        assert_eq!(
            volume.deflines(0, true, true, true).unwrap_err(),
            "Header offsets exceed PHR file size"
        );
        std::fs::remove_file(prefix.with_extension("pin")).unwrap();
        std::fs::remove_file(prefix.with_extension("phr")).unwrap();
        std::fs::remove_file(prefix.with_extension("psq")).unwrap();
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
