//! BLAST protein header (`.phr`) decoding.
//!
//! Mirrors `diamond/src/data/blastdb/phr.cpp` and the related `volume.h`
//! declarations. `volume` re-exports the public types and helpers for API
//! compatibility.

use super::asn1::{decode, Node};
use super::ber::decode_integer;
use super::volume::BlastVolume;
use crate::basic::value::TaxId;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SeqId<N = String> {
    pub type_: String,
    pub value: N,
    pub version: Option<i64>,
    pub chain: Option<N>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BlastDefLine<N = String> {
    pub title: N,
    pub seqids: Vec<SeqId<N>>,
    pub taxid: Option<TaxId>,
}

pub type ByteSeqId = SeqId<Vec<u8>>;
pub type ByteBlastDefLine = BlastDefLine<Vec<u8>>;
impl ByteSeqId {
    pub fn try_into_string(self) -> Result<SeqId, String> {
        Ok(SeqId {
            type_: self.type_,
            value: String::from_utf8(self.value).map_err(|error| error.to_string())?,
            version: self.version,
            chain: self
                .chain
                .map(String::from_utf8)
                .transpose()
                .map_err(|error| error.to_string())?,
        })
    }
}
impl ByteBlastDefLine {
    pub fn try_into_string(self) -> Result<BlastDefLine, String> {
        Ok(BlastDefLine {
            title: String::from_utf8(self.title).map_err(|error| error.to_string())?,
            seqids: self
                .seqids
                .into_iter()
                .map(ByteSeqId::try_into_string)
                .collect::<Result<_, _>>()?,
            taxid: self.taxid,
        })
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

fn decode_seqid_into_bytes(node: &Node, seqid: &mut ByteSeqId) {
    for n4 in &node.children {
        match n4.tag.tag_number {
            1 => {
                for n5 in &n4.children {
                    if n5.tag.tag_number == 26 {
                        seqid.value = n5.value.clone();
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

fn decode_seqid_bytes(node: &Node) -> ByteSeqId {
    let mut seqid = ByteSeqId::default();
    for n1 in &node.children {
        if n1.tag.tag_number == 16 {
            for n2 in &n1.children {
                match n2.tag.tag_number {
                    0 | 1 | 4 | 5 | 7 | 9 | 12 | 15 | 16 => {
                        seqid.type_ = tag_name_from_number(n2.tag.tag_number);
                        decode_seqid_into_bytes(n2, &mut seqid);
                        for n3 in &n2.children {
                            if n3.tag.tag_number == 16 {
                                decode_seqid_into_bytes(n3, &mut seqid);
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
                                                    seqid.value = n5.value.clone();
                                                }
                                            }
                                        }
                                        3 => {
                                            for n5 in &n4.children {
                                                if n5.tag.tag_number == 26 {
                                                    seqid.chain = Some(n5.value.clone());
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

fn decode_defline_bytes(node: &Node, full_titles: bool, taxids: bool) -> ByteBlastDefLine {
    let mut defline = ByteBlastDefLine::default();
    for n1 in &node.children {
        match n1.tag.tag_number {
            0 => {
                if full_titles {
                    for n2 in &n1.children {
                        if n2.tag.tag_number == 26 {
                            defline.title = n2.value.clone();
                        }
                    }
                }
            }
            1 => defline.seqids.push(decode_seqid_bytes(n1)),
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

pub fn decode_deflines_bytes(
    header_data: &[u8],
    all: bool,
    full_titles: bool,
    taxids: bool,
) -> Result<Vec<ByteBlastDefLine>, String> {
    let mut out = Vec::new();
    let nodes = decode(header_data).map_err(|e| e.to_string())?;
    if nodes.is_empty() {
        return Ok(out);
    }
    for node in &nodes[0].children {
        out.push(decode_defline_bytes(node, full_titles, taxids));
        if !all && !taxids {
            break;
        }
    }
    Ok(out)
}

impl BlastVolume {
    pub fn deflines_bytes(
        &mut self,
        oid: u32,
        all: bool,
        full_titles: bool,
        taxids: bool,
    ) -> Result<Vec<ByteBlastDefLine>, String> {
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
            self.io_errors.capture(
                self.phr_mapping
                    .seek(header_offset as i64, std::io::SeekFrom::Start(0))
                    .map_err(std::io::Error::other),
            )?;
        }
        self.hdr_ptr = oid + 1;
        let data = self.io_errors.capture(
            self.phr_mapping
                .read(header_length)
                .map_err(std::io::Error::other),
        )?;
        decode_deflines_bytes(data, all, full_titles, taxids)
    }

    pub fn raw_deflines(&mut self, count: u32) -> Result<Vec<u8>, String> {
        let len = (self.index.header_index[(self.hdr_ptr + count) as usize]
            - self.index.header_index[self.hdr_ptr as usize]) as usize;
        let mut data = vec![0u8; len];
        self.io_errors.capture(
            self.phr_mapping
                .read_exact(&mut data)
                .map_err(std::io::Error::other),
        )?;
        self.hdr_ptr += count;
        Ok(data)
    }

    pub fn id_len(&self, oid: u32) -> usize {
        id_len(&self.index.header_index, oid as usize)
    }
}

pub fn id_len(header_index: &[u32], oid: usize) -> usize {
    (header_index[oid + 1] - header_index[oid]) as usize
}

pub fn format_seqid_bytes<N: AsRef<[u8]>>(id: &SeqId<N>) -> Vec<u8> {
    if id.value.as_ref().is_empty() {
        return b"N/A".to_vec();
    }
    let mut output = id.value.as_ref().to_vec();
    if let Some(version) = id.version {
        output.push(b'.');
        output.extend_from_slice(version.to_string().as_bytes());
    }
    if let Some(chain) = &id.chain {
        if !chain.as_ref().is_empty() {
            output.push(b'_');
            output.extend_from_slice(chain.as_ref());
        }
    }
    output
}
pub fn build_title_bytes<N: AsRef<[u8]>>(
    deflines: &[BlastDefLine<N>],
    delimiter: &[u8],
    all: bool,
) -> Vec<u8> {
    let mut output = Vec::new();
    let delimiter = &delimiter[..delimiter
        .iter()
        .position(|&byte| byte == 0)
        .unwrap_or(delimiter.len())];
    for (index, line) in deflines.iter().enumerate() {
        if index != 0 {
            if !all {
                break;
            }
            output.extend_from_slice(delimiter);
        }
        if let Some(id) = line.seqids.first() {
            output.extend_from_slice(&format_seqid_bytes(id));
        } else {
            output.extend_from_slice(b"N/A");
        }
        output.push(b' ');
        output.extend_from_slice(line.title.as_ref());
    }
    if output.is_empty() {
        output.extend_from_slice(b"N/A");
    }
    output
}
pub fn format_seqid(id: &SeqId) -> String {
    String::from_utf8(format_seqid_bytes(id)).expect("String SeqId must render UTF-8")
}
pub fn build_title(deflines: &[BlastDefLine], delimiter: &str, all: bool) -> String {
    String::from_utf8(build_title_bytes(deflines, delimiter.as_bytes(), all))
        .expect("String deflines must render UTF-8")
}
pub fn decode_deflines(
    data: &[u8],
    all: bool,
    full_titles: bool,
    taxids: bool,
) -> Result<Vec<BlastDefLine>, String> {
    decode_deflines_bytes(data, all, full_titles, taxids)?
        .into_iter()
        .map(ByteBlastDefLine::try_into_string)
        .collect()
}
impl BlastVolume {
    pub fn deflines(
        &mut self,
        oid: u32,
        all: bool,
        full_titles: bool,
        taxids: bool,
    ) -> Result<Vec<BlastDefLine>, String> {
        self.deflines_bytes(oid, all, full_titles, taxids)?
            .into_iter()
            .map(ByteBlastDefLine::try_into_string)
            .collect()
    }
}
#[cfg(test)]
fn decode_seqid(node: &Node) -> SeqId {
    decode_seqid_bytes(node).try_into_string().unwrap()
}
#[cfg(test)]
fn decode_defline(node: &Node, full_titles: bool, taxids: bool) -> BlastDefLine {
    decode_defline_bytes(node, full_titles, taxids)
        .try_into_string()
        .unwrap()
}

#[cfg(test)]
mod tests {
    use super::super::asn1::{Class, TagInfo};
    use super::*;

    fn node(tag_number: u32, value: &[u8], children: Vec<Node>) -> Node {
        Node {
            tag: TagInfo {
                tag_class: Class::ContextSpecific,
                constructed: !children.is_empty(),
                tag_number,
            },
            value: value.to_vec(),
            children,
        }
    }

    #[test]
    fn raw_pdb_chain_and_c_string_delimiter_preserve_source_domains() {
        let molecule = node(0, &[], vec![node(26, b"P\xff\0ID", vec![])]);
        let chain = node(3, &[], vec![node(26, b"\x80", vec![])]);
        let pdb = node(14, &[], vec![node(16, &[], vec![molecule, chain])]);
        let input = node(1, &[], vec![node(16, &[], vec![pdb])]);
        let id = decode_seqid_bytes(&input);
        assert_eq!(format_seqid_bytes(&id), b"P\xff\0ID_\x80");
        assert!(id.clone().try_into_string().is_err());
        let lines = vec![
            ByteBlastDefLine {
                title: b"T\xfe\0tail".to_vec(),
                seqids: vec![id],
                taxid: None
            };
            2
        ];
        assert_eq!(
            build_title_bytes(&lines, b"|\0ignored", true),
            b"P\xff\0ID_\x80 T\xfe\0tail|P\xff\0ID_\x80 T\xfe\0tail"
        );
    }

    #[test]
    fn pdb_seqid_decodes_molecule_and_chain() {
        let molecule = node(0, &[], vec![node(26, b"1ABC", vec![])]);
        let chain = node(3, &[], vec![node(26, b"B", vec![])]);
        let pdb = node(14, &[], vec![node(16, &[], vec![molecule, chain])]);
        let input = node(1, &[], vec![node(16, &[], vec![pdb])]);
        let id = decode_seqid(&input);
        assert_eq!(id.type_, "pdb");
        assert_eq!(id.value, "1ABC");
        assert_eq!(id.chain.as_deref(), Some("B"));
        assert_eq!(format_seqid(&id), "1ABC_B");
    }

    #[test]
    fn title_and_taxid_flags_are_independent() {
        let defline = node(
            16,
            &[],
            vec![
                node(0, &[], vec![node(26, b"description", vec![])]),
                node(2, &[], vec![node(2, &[0x2a], vec![])]),
            ],
        );
        let neither = decode_defline(&defline, false, false);
        assert_eq!(neither.title, "");
        assert_eq!(neither.taxid, None);
        let both = decode_defline(&defline, true, true);
        assert_eq!(both.title, "description");
        assert_eq!(both.taxid, Some(42));
    }

    #[test]
    fn title_formatting_preserves_empty_and_all_rules() {
        assert_eq!(build_title(&[], " | ", true), "N/A");
        let lines = [
            BlastDefLine {
                title: "first".into(),
                seqids: vec![SeqId {
                    value: "A".into(),
                    ..Default::default()
                }],
                taxid: None,
            },
            BlastDefLine {
                title: "second".into(),
                seqids: vec![SeqId {
                    value: "B".into(),
                    ..Default::default()
                }],
                taxid: None,
            },
        ];
        assert_eq!(build_title(&lines, " | ", false), "A first");
        assert_eq!(build_title(&lines, " | ", true), "A first | B second");
    }

    #[test]
    fn malformed_asn1_errors_are_propagated() {
        let error = decode_deflines(&[0x30, 0x05, 0x01], true, true, true).unwrap_err();
        assert!(!error.is_empty());
    }
}
