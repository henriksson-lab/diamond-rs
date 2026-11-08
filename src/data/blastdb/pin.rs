//! BLAST protein index (`.pin`) parsing and indexed chunk decoding.
//!
//! Mirrors `diamond/src/data/blastdb/pin.cpp` and the corresponding declarations
//! in `volume.h`. The compatibility surface remains re-exported by `volume`.

use super::ber::{read_be32, read_le64, read_pascal_bytes};
use super::psq::decode_protein_sequence;
use super::volume::{
    build_title_bytes, decode_deflines_bytes, format_seqid, format_seqid_bytes, BlastDefLine,
    BlastVolume, ByteBlastDefLine, DbFilter, DecodedPackage, RawChunk, SequenceFileFlags,
};
use crate::util::io::File as DiamondFile;
use std::collections::{BTreeSet, HashMap};

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PinIndex {
    pub version: u32,
    pub is_protein: bool,
    pub volume_number: u32,
    pub title: Vec<u8>,
    pub lmdb_file: Vec<u8>,
    pub date: Vec<u8>,
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
        Self::new_path(std::path::Path::new(path), idx, begin, end, load_index)
            .map_err(|e| e.to_string())
    }
    pub fn new_path(
        path: &std::path::Path,
        idx: i32,
        begin: u64,
        end: u64,
        load_index: bool,
    ) -> std::io::Result<Self> {
        let errors = crate::util::io::IoErrorOwner::default();
        let result = Self::new_path_owned(path, idx, begin, end, load_index, &errors);
        errors.finish(result)
    }
    pub(crate) fn new_path_owned(
        path: &std::path::Path,
        idx: i32,
        begin: u64,
        end: u64,
        load_index: bool,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        let phr_mapping = errors.capture(
            DiamondFile::open_path(&super::pal::append_suffix(path, ".phr"), "rb")
                .map_err(std::io::Error::other),
        )?;
        let psq_mapping = errors.capture(
            DiamondFile::open_path(&super::pal::append_suffix(path, ".psq"), "rb")
                .map_err(std::io::Error::other),
        )?;
        let mut pin = errors.capture(
            DiamondFile::open_path(&super::pal::append_suffix(path, ".pin"), "rb")
                .map_err(std::io::Error::other),
        )?;
        let index = Self::parse_pin_file_owned(&mut pin, load_index, errors)?;
        errors.capture(pin.close().map_err(std::io::Error::other))?;
        Ok(Self {
            idx,
            begin,
            end,
            index,
            phr_mapping,
            psq_mapping,
            seq_ptr: 0,
            hdr_ptr: 0,
            io_errors: errors.clone(),
        })
    }

    pub fn parse_pin_file(mapping: &mut DiamondFile, load_index: bool) -> Result<PinIndex, String> {
        Self::parse_pin_file_owned(
            mapping,
            load_index,
            &crate::util::io::IoErrorOwner::default(),
        )
    }
    fn parse_pin_file_owned(
        mapping: &mut DiamondFile,
        load_index: bool,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<PinIndex, String> {
        let mut index = PinIndex {
            version: errors.capture(read_be32(mapping).map_err(std::io::Error::other))?,
            ..Default::default()
        };
        if index.version != 4 && index.version != 5 {
            return Err(format!(
                "Unsupported database format version: {}",
                index.version
            ));
        }

        index.is_protein = errors.capture(read_be32(mapping).map_err(std::io::Error::other))? == 1;
        if index.version == 5 {
            index.volume_number =
                errors.capture(read_be32(mapping).map_err(std::io::Error::other))?;
        }
        index.title = errors.capture(read_pascal_bytes(mapping).map_err(std::io::Error::other))?;
        if index.version == 5 {
            index.lmdb_file =
                errors.capture(read_pascal_bytes(mapping).map_err(std::io::Error::other))?;
        }
        index.date = errors.capture(read_pascal_bytes(mapping).map_err(std::io::Error::other))?;
        index.num_oids = errors.capture(read_be32(mapping).map_err(std::io::Error::other))?;
        index.total_length = errors.capture(read_le64(mapping).map_err(std::io::Error::other))?;
        index.max_length = errors.capture(read_be32(mapping).map_err(std::io::Error::other))?;
        if !load_index {
            return Ok(index);
        }

        let count = index.num_oids as usize + 1;
        index.header_index.reserve(count);
        index.sequence_index.reserve(count);
        for _ in 0..count {
            index
                .header_index
                .push(errors.capture(read_be32(mapping).map_err(std::io::Error::other))?);
        }
        for _ in 0..count {
            index
                .sequence_index
                .push(errors.capture(read_be32(mapping).map_err(std::io::Error::other))?);
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
        self.io_errors.capture(
            self.phr_mapping
                .seek(0, std::io::SeekFrom::Start(0))
                .map_err(std::io::Error::other),
        )?;
        self.io_errors.capture(
            self.psq_mapping
                .seek(0, std::io::SeekFrom::Start(0))
                .map_err(std::io::Error::other),
        )
    }
}

pub(super) fn acc_filter(deflines: &[BlastDefLine], accs: &mut HashMap<String, bool>) -> bool {
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

pub(super) fn acc_filter_bytes(
    deflines: &[ByteBlastDefLine],
    accs: &mut HashMap<Vec<u8>, bool>,
) -> bool {
    for line in deflines {
        for id in &line.seqids {
            if let Some(found) = accs.get_mut(&id.value) {
                *found = true;
                return true;
            }
            if id.version.is_some() || id.chain.is_some() {
                if let Some(found) = accs.get_mut(&format_seqid_bytes(id)) {
                    *found = true;
                    return true;
                }
            }
        }
    }
    false
}

impl RawChunk {
    /// Compatibility for caller-provided UTF-8 accession keys. Decoded title
    /// storage has always been a byte StringSet and now stays byte-canonical.
    pub fn decode(
        &self,
        flags: SequenceFileFlags,
        filter: Option<&DbFilter>,
        accs: Option<&mut HashMap<String, bool>>,
    ) -> Result<DecodedPackage, String> {
        if let Some(accs) = accs {
            let mut bytes: HashMap<Vec<u8>, bool> = accs
                .iter()
                .map(|(key, found)| (key.as_bytes().to_vec(), *found))
                .collect();
            let result = self.decode_bytes(flags, filter, Some(&mut bytes));
            for (key, found) in accs {
                if let Some(updated) = bytes.get(key.as_bytes()) {
                    *found = *updated;
                }
            }
            result
        } else {
            self.decode_bytes(flags, filter, None)
        }
    }

    pub fn decode_bytes(
        &self,
        flags: SequenceFileFlags,
        filter: Option<&DbFilter>,
        mut accs: Option<&mut HashMap<Vec<u8>, bool>>,
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
                    let deflines = decode_deflines_bytes(
                        &self.phr_data[phr_ptr..phr_ptr + header_len],
                        all_seqids,
                        full_titles,
                        taxids,
                    )?;
                    if let Some(accs) = accs.as_deref_mut() {
                        selected = acc_filter_bytes(&deflines, accs);
                    }
                    if selected && titles {
                        pkg.ids
                            .push_back(&build_title_bytes(&deflines, b"\x01", true));
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
            title: Vec::new(),
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
            title: Vec::new(),
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
