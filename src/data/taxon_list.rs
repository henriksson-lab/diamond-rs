use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Read};

use crate::basic::value::{OId, TaxId};
use crate::data::compact_array::CompactArray;
use crate::util::algo::{join_sorted_lists, write_varuint32};
use crate::util::sequence::{get_accession, AccessionParsing};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaxonList {
    array: CompactArray,
}

impl TaxonList {
    pub fn new(data: Vec<u8>, size: usize, data_size: usize) -> Result<Self, String> {
        Ok(Self {
            array: CompactArray::new(data, size, data_size)?,
        })
    }

    pub fn get(&self, i: usize) -> Result<Vec<i32>, String> {
        self.array.get(i)
    }

    /// C++ `TaxonList(Deserializer&, size, data_size)`.
    pub fn read_from<R: Read>(
        mut reader: R,
        size: usize,
        data_size: usize,
    ) -> Result<Self, String> {
        let mut data = vec![0; data_size];
        reader.read_exact(&mut data).map_err(|e| e.to_string())?;
        Self::new(data, size, data_size)
    }

    pub fn size(&self) -> usize {
        self.array.size()
    }

    pub fn data(&self) -> &[u8] {
        self.array.data()
    }

    pub fn build(
        acc2oid: &[(String, OId)],
        mapping_rows: &str,
        seqs: OId,
        no_parse_seqids: bool,
    ) -> Result<TaxonListBuild, String> {
        let (acc2taxid, acc_stats) = load_mapping_file(mapping_rows.as_bytes(), no_parse_seqids)?;
        let mut sorted_acc2oid = acc2oid.to_vec();
        sorted_acc2oid.sort();

        let joined = join_sorted_lists(
            &sorted_acc2oid,
            &acc2taxid,
            |row| row.0.clone(),
            |row| row.0.clone(),
            |oid, taxid| (oid.1, taxid.1),
        )?;
        let acc_matched = joined.len();
        let mut oid2taxids: BTreeMap<OId, BTreeSet<TaxId>> = BTreeMap::new();
        for (oid, taxid) in joined {
            oid2taxids.entry(oid).or_default().insert(taxid);
        }

        let mut data = Vec::new();
        let mut mapped_seqs = 0usize;
        for oid in 0..seqs {
            let mut tax_ids = oid2taxids.remove(&oid).unwrap_or_default();
            tax_ids.remove(&0);
            serialize_taxid_set(&tax_ids, &mut data);
            if !tax_ids.is_empty() {
                mapped_seqs += 1;
            }
        }
        let taxon_list = TaxonList::new(data.clone(), seqs as usize, data.len())?;
        Ok(TaxonListBuild {
            taxon_list,
            stats: TaxonListStats {
                accessions_in_database: acc2oid.len(),
                entries_in_accession_to_taxid_file: acc2taxid.len(),
                database_accessions_mapped_to_taxid: acc_matched,
                database_sequences_mapped_to_taxid: mapped_seqs,
            },
            accession_stats: acc_stats,
        })
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TaxonListStats {
    pub accessions_in_database: usize,
    pub entries_in_accession_to_taxid_file: usize,
    pub database_accessions_mapped_to_taxid: usize,
    pub database_sequences_mapped_to_taxid: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaxonListBuild {
    pub taxon_list: TaxonList,
    pub stats: TaxonListStats,
    pub accession_stats: AccessionParsing,
}

pub fn mapping_file_format(header: &str) -> Result<i32, String> {
    let mut fields = header.split('\t');
    let field1 = fields.next().unwrap_or("");
    let field2 = fields.next().unwrap_or("");
    if field1 == "accession" && field2 == "accession.version" {
        let field3 = fields.next().unwrap_or("");
        let field4 = fields.next().unwrap_or("");
        if field3 == "taxid" && field4 == "gi" && only_optional_trailing_empty(fields) {
            return Ok(0);
        }
    } else if field1 == "accession.version"
        && field2 == "taxid"
        && only_optional_trailing_empty(fields)
    {
        return Ok(1);
    }
    Err("Accession mapping file header has to be in one of these formats:\naccession\taccession.version\ttaxid\tgi\naccession.version\ttaxid".to_string())
}

fn only_optional_trailing_empty<'a>(mut fields: impl Iterator<Item = &'a str>) -> bool {
    match fields.next() {
        None => true,
        Some("") => fields.next().is_none(),
        Some(_) => false,
    }
}

pub fn load_mapping_file<R: BufRead>(
    reader: R,
    no_parse_seqids: bool,
) -> Result<(Vec<(String, TaxId)>, AccessionParsing), String> {
    let mut lines = reader.lines();
    let header = lines
        .next()
        .transpose()
        .map_err(|e| e.to_string())?
        .unwrap_or_default();
    let format = mapping_file_format(&header)?;
    let mut rows = Vec::new();
    let mut stats = AccessionParsing::default();
    let mut last = String::new();

    for (idx, line) in lines.enumerate() {
        let line_count = idx + 2;
        let line = line.map_err(|e| e.to_string())?;
        if line.is_empty() {
            break;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let (mut accession, taxid) = if format == 0 {
            if fields.len() < 3 {
                return Err(format!("Malformed input in line {line_count}"));
            }
            (fields[1].to_string(), fields[2])
        } else {
            if fields.len() < 2 {
                return Err(format!("Malformed input in line {line_count}"));
            }
            (fields[0].to_string(), fields[1])
        };
        if accession.is_empty() {
            return Err(format!("Empty accession field in line {line_count}"));
        }
        let taxid = taxid
            .trim_start_matches(|c: char| c.is_ascii_whitespace())
            .parse::<TaxId>()
            .map_err(|_| format!("Malformed input in line {line_count}"))?;
        if !no_parse_seqids {
            if let Some(i) = accession.find(":PDB=") {
                accession.truncate(i);
                stats.pdb_suffix += 1;
            }
            accession = get_accession(&accession, &mut stats);
        }
        if accession != last {
            rows.push((accession.clone(), taxid));
        }
        last = accession;
    }
    rows.sort();
    Ok((rows, stats))
}

pub fn serialize_taxid_set(tax_ids: &BTreeSet<TaxId>, out: &mut Vec<u8>) {
    write_varuint32(tax_ids.len() as u32, out);
    for &tax_id in tax_ids {
        write_varuint32(tax_id as u32, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapping_file_format() {
        assert_eq!(
            mapping_file_format("accession\taccession.version\ttaxid\tgi").unwrap(),
            0
        );
        assert_eq!(mapping_file_format("accession.version\ttaxid").unwrap(), 1);
        // C++ tokenizer accepts one delimiter immediately before its NUL.
        assert_eq!(
            mapping_file_format("accession.version\ttaxid\t").unwrap(),
            1
        );
        assert!(mapping_file_format("accession.version\ttaxid\t\t").is_err());
        assert!(mapping_file_format("accession\ttaxid").is_err());
    }

    #[test]
    fn test_load_mapping_file_formats_parse_and_dedup() {
        let data = "accession\taccession.version\ttaxid\tgi\nA\tsp|P1.2|x\t10\t0\nA\tsp|P1.2|x\t10\t0\nB\tUniRef90_Q2.1\t20\t0\nC\tfoo:PDB=1\t30\t0\n";
        let (rows, stats) = load_mapping_file(data.as_bytes(), false).unwrap();
        assert_eq!(
            rows,
            vec![
                ("P1".to_string(), 10),
                ("Q2".to_string(), 20),
                ("foo".to_string(), 30)
            ]
        );
        assert_eq!(stats.prefix_before_pipe, 2);
        assert_eq!(stats.suffix_after_pipe, 2);
        assert_eq!(stats.uniref_prefix, 1);
        assert_eq!(stats.pdb_suffix, 1);

        let data = "accession.version\ttaxid\nraw|id.1\t7\n";
        let (rows, _) = load_mapping_file(data.as_bytes(), true).unwrap();
        assert_eq!(rows, vec![("raw|id.1".to_string(), 7)]);
    }

    #[test]
    fn test_load_mapping_file_errors() {
        assert!(load_mapping_file("bad\theader\n".as_bytes(), false)
            .unwrap_err()
            .starts_with("Accession mapping file header"));
        assert_eq!(
            load_mapping_file("accession.version\ttaxid\n\t1\n".as_bytes(), false).unwrap_err(),
            "Empty accession field in line 2"
        );
        assert_eq!(
            load_mapping_file("accession.version\ttaxid\nA\tx\n".as_bytes(), false).unwrap_err(),
            "Malformed input in line 2"
        );
        assert!(load_mapping_file("".as_bytes(), false)
            .unwrap_err()
            .starts_with("Accession mapping file header"));

        let duplicate_after_sort = "accession.version\ttaxid\nA\t1\nB\t2\nA\t3\n";
        let err = TaxonList::build(
            &[("A".into(), 0), ("B".into(), 1)],
            duplicate_after_sort,
            2,
            true,
        )
        .unwrap_err();
        assert_eq!(err, "Duplicate keys: A");
    }

    #[test]
    fn test_serialize_taxid_set_and_construct_compact_array() {
        let mut set = BTreeSet::new();
        set.insert(7);
        set.insert(2);
        let mut data = Vec::new();
        serialize_taxid_set(&set, &mut data);
        assert_eq!(data, vec![5, 5, 15]);
        let list = TaxonList::new(data.clone(), 1, data.len()).unwrap();
        assert_eq!(list.get(0).unwrap(), vec![2, 7]);
        let from_reader = TaxonList::read_from(data.as_slice(), 1, data.len()).unwrap();
        assert_eq!(from_reader.get(0).unwrap(), vec![2, 7]);
        assert!(TaxonList::read_from(&data[..2], 1, data.len()).is_err());
    }

    #[test]
    fn test_taxon_list_build_join_and_stats() {
        let acc2oid = vec![
            ("P1".to_string(), 0),
            ("Q2".to_string(), 1),
            ("Q2".to_string(), 2),
            ("missing".to_string(), 3),
        ];
        let mapping = "accession.version\ttaxid\nP1\t10\nQ2\t0\nQ2\t20\n";
        let build = TaxonList::build(&acc2oid, mapping, 4, true).unwrap();
        assert_eq!(build.taxon_list.size(), 4);
        assert_eq!(build.taxon_list.get(0).unwrap(), vec![10]);
        assert_eq!(build.taxon_list.get(1).unwrap(), Vec::<i32>::new());
        assert_eq!(build.taxon_list.get(2).unwrap(), Vec::<i32>::new());
        assert_eq!(build.taxon_list.get(3).unwrap(), Vec::<i32>::new());
        assert_eq!(build.taxon_list.data(), &[3, 21, 1, 1, 1]);
        assert_eq!(
            build.stats,
            TaxonListStats {
                accessions_in_database: 4,
                entries_in_accession_to_taxid_file: 2,
                database_accessions_mapped_to_taxid: 3,
                database_sequences_mapped_to_taxid: 1,
            }
        );
    }
}
