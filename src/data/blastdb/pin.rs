//! BLAST protein index (`.pin`) parsing and indexed chunk decoding.
//!
//! Mirrors `diamond/src/data/blastdb/pin.cpp` and the corresponding declarations
//! in `volume.h`. The compatibility surface remains re-exported by `volume`.

use super::ber::{read_be32, read_le64, read_pascal_string};
use super::psq::decode_protein_sequence;
use super::volume::{
    build_title, decode_deflines, format_seqid, BlastDefLine, BlastVolume, DbFilter,
    DecodedPackage, RawChunk, SequenceFileFlags,
};
use crate::util::io::File as DiamondFile;
use std::collections::{BTreeSet, HashMap};

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

        index.is_protein = read_be32(mapping).map_err(|e| e.to_string())? == 1;
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
        let mut letter_count = 0usize;
        while end < self.index.num_oids && letter_count < letters {
            letter_count += self.length(end) as usize;
            end += 1;
        }
        let mut chunk = RawChunk {
            letters: letter_count,
            begin: u64::from(begin) + self.begin,
            end: u64::from(end) + self.begin,
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

fn acc_filter(deflines: &[BlastDefLine], accs: &mut HashMap<String, bool>) -> bool {
    for defline in deflines {
        for seqid in &defline.seqids {
            if accs.contains_key(&seqid.value) {
                accs.insert(seqid.value.clone(), true);
                return true;
            }
            if seqid.version.is_some() || seqid.chain.is_some() {
                let formatted = format_seqid(seqid);
                if accs.contains_key(&formatted) {
                    accs.insert(formatted, true);
                    return true;
                }
            }
        }
    }
    false
}

impl RawChunk {
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
            let mut selected = filter.map(|filter| filter.get(oid)).unwrap_or(true);
            if titles || taxids || accs.is_some() {
                let header_len = (self.phr_index[i + 1] - self.phr_index[i]) as usize;
                if selected || accs.is_some() {
                    let deflines = decode_deflines(
                        &self.phr_data[phr_ptr..phr_ptr + header_len],
                        all_seqids,
                        full_titles,
                        taxids,
                    )?;
                    if let Some(accs) = accs.as_deref_mut() {
                        selected = acc_filter(&deflines, accs);
                    }
                    if selected && titles {
                        pkg.ids
                            .push_back(build_title(&deflines, "\x01", true).as_bytes());
                    }
                    if selected && taxids {
                        let taxids: BTreeSet<_> = deflines
                            .iter()
                            .filter_map(|defline| defline.taxid)
                            .collect();
                        pkg.taxids
                            .extend(taxids.into_iter().map(|taxid| (oid, taxid)));
                    }
                }
                phr_ptr += header_len;
            }

            if seqs {
                let sequence_len = (self.seq_index[i + 1] - self.seq_index[i]) as usize;
                if selected {
                    let sequence =
                        decode_protein_sequence(&self.seq_data[seq_ptr..seq_ptr + sequence_len])?;
                    pkg.seqs.push(&sequence);
                }
                seq_ptr += sequence_len;
            }
            if selected {
                pkg.oids.push(oid);
            }
        }
        Ok(pkg)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accession_filter_prefers_raw_then_formatted_ids() {
        let deflines = [BlastDefLine {
            title: String::new(),
            seqids: vec![super::super::volume::SeqId {
                value: "XP_1".into(),
                version: Some(2),
                chain: Some("A".into()),
                ..Default::default()
            }],
            taxid: None,
        }];
        let mut raw = HashMap::from([("XP_1".to_string(), false)]);
        assert!(acc_filter(&deflines, &mut raw));
        assert!(raw["XP_1"]);

        let mut formatted = HashMap::from([("XP_1.2_A".to_string(), false)]);
        assert!(acc_filter(&deflines, &mut formatted));
        assert!(formatted["XP_1.2_A"]);
    }

    #[test]
    fn accession_filter_leaves_absent_ids_unchanged() {
        let deflines = [BlastDefLine {
            title: String::new(),
            seqids: vec![super::super::volume::SeqId {
                value: "XP_1".into(),
                version: Some(2),
                ..Default::default()
            }],
            taxid: None,
        }];
        let mut wanted = HashMap::from([("OTHER".to_string(), false)]);
        assert!(!acc_filter(&deflines, &mut wanted));
        assert!(!wanted["OTHER"]);
    }
}
