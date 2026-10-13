use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use rayon::prelude::*;

use super::dmnd::{ReferenceHeader, ReferenceHeader2, CURRENT_DB_VERSION_PROT};
use super::fasta::{self, FastaRecord};
use crate::basic::value::OId;
use crate::basic::value::SequenceType;
use crate::data::taxon_list::TaxonList;
use crate::data::taxonomy::Taxonomy;
use crate::data::taxonomy_nodes::TaxonomyNodes;
use crate::masking::tantan::mask_tantan;
use crate::util::algo::ExternalSorter;
use crate::util::hash::murmurhash3_x64_128;
use crate::util::sequence::{accession_from_title, AccessionParsing};

/// Total size of the two on-disk reference headers (40 + 8 size prefix + 48 data).
/// Matches `out->tell()` after `*out << header << header2` in diamond/src/data/dmnd/dmnd.cpp.
const HEADER_TOTAL: u64 = 40 + 8 + 48;

/// `DatabaseFile::make_db` asks `FastaFile::load_seqs` for one billion
/// letters per block.  The loader includes the record that reaches/crosses
/// this limit and then returns the owned block for masking, writing and
/// hashing.
const MAKEDB_BLOCK_LETTERS: u64 = 1_000_000_000;

/// Build a DIAMOND database from FASTA input files.
///
/// This is the Rust implementation of the `makedb` command.
pub fn build_db(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    build_db_with_options(input_files, output_path, seq_type, &BuildOptions::default())
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BuildOptions {
    pub taxonmap: Option<String>,
    pub taxonnodes: Option<String>,
    pub taxonnames: Option<String>,
    pub no_parse_seqids: bool,
}

pub fn build_db_with_options(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
    options: &BuildOptions,
) -> io::Result<BuildStats> {
    build_db_with_block_limit(
        input_files,
        output_path,
        seq_type,
        MAKEDB_BLOCK_LETTERS,
        options,
    )
}

fn build_db_with_block_limit(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
    block_letter_limit: u64,
    options: &BuildOptions,
) -> io::Result<BuildStats> {
    debug_assert!(block_letter_limit > 0);
    let result = (|| {
        let mut output = DmndBuildWriter::create(output_path, seq_type, options.clone())?;
        let mut block = Vec::new();
        let mut block_letters = 0u64;

        for input_file in input_files {
            fasta::for_each_fasta_file(Path::new(input_file), seq_type, |record| {
                block_letters = block_letters.saturating_add(record.sequence.len() as u64);
                block.push(record);
                // Upstream's one-pass loader is a do/while: the record that
                // reaches the bound belongs to this block.
                if block_letters >= block_letter_limit {
                    output.write_block(&mut block)?;
                    block_letters = 0;
                }
                Ok(())
            })?;
        }
        output.write_block(&mut block)?;
        output.finish()
    })();

    // `DatabaseFile::make_db` removes an incomplete output on any exception.
    if result.is_err() {
        let _ = std::fs::remove_file(output_path);
    }
    result
}

/// Build a DIAMOND database from in-memory sequences.
pub fn build_db_from_records(
    records: &[FastaRecord],
    output_path: &str,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    // Keep this convenience API non-mutating.  The production `makedb` path
    // owns its records and masks them in place below, just like C++ mutates
    // each loaded Block before writing it.
    let mut records = records.to_vec();
    write_dmnd(output_path, &mut records, seq_type)
}

/// Statistics from database building.
#[derive(Debug, Clone)]
pub struct BuildStats {
    pub sequences: u64,
    pub letters: u64,
}

fn write_dmnd(
    output_path: &str,
    records: &mut Vec<FastaRecord>,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    let mut output = DmndBuildWriter::create(output_path, seq_type, BuildOptions::default())?;
    output.write_block(records)?;
    output.finish()
}

struct DmndBuildWriter {
    writer: BufWriter<std::fs::File>,
    output_path: String,
    seq_type: SequenceType,
    pos_array: Vec<(u64, u32)>,
    current_pos: u64,
    hash: [u8; 16],
    sequences: u64,
    letters: u64,
    options: BuildOptions,
    accessions: ExternalSorter<(String, OId)>,
    accession_stats: AccessionParsing,
}

impl DmndBuildWriter {
    fn create(
        output_path: &str,
        seq_type: SequenceType,
        options: BuildOptions,
    ) -> io::Result<Self> {
        let file = std::fs::File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        let mut header = ReferenceHeader::new();
        header.db_version = match seq_type {
            SequenceType::AminoAcid => CURRENT_DB_VERSION_PROT,
            SequenceType::Nucleotide => 4,
        };
        header.write_to(&mut writer)?;
        ReferenceHeader2::new().write_record_to(&mut writer)?;
        Ok(Self {
            writer,
            output_path: output_path.to_owned(),
            seq_type,
            pos_array: Vec::new(),
            current_pos: HEADER_TOTAL,
            hash: [0; 16],
            sequences: 0,
            letters: 0,
            options,
            accessions: ExternalSorter::default(),
            accession_stats: AccessionParsing::default(),
        })
    }

    fn write_block(&mut self, records: &mut Vec<FastaRecord>) -> io::Result<()> {
        if records.is_empty() {
            return Ok(());
        }
        // The upstream block is masked in place and dropped immediately after
        // this method.  No second complete encoded copy is retained.
        if self.seq_type == SequenceType::AminoAcid {
            records.par_iter_mut().for_each(|record| {
                mask_tantan(&mut record.sequence);
            });
        }
        for record in records.iter() {
            // `Letter` is i8 and DMND stores its bit pattern verbatim.  Both
            // element types have size/alignment one, and the borrow cannot outlive
            // `record.sequence`.
            let seq_bytes = unsafe {
                std::slice::from_raw_parts(
                    record.sequence.as_ptr().cast::<u8>(),
                    record.sequence.len(),
                )
            };
            // C++ `dmnd.cpp:313-314` throws "File format error: sequence of
            // length 0" — accepting empty sequences here writes a SeqInfo with
            // `seq_len == 0`, which the C++ partition logic in `dmnd.cpp:483`
            // mistakes for the sentinel (`r_next.seq_len == 0` terminates
            // partitioning), corrupting downstream reads of the same DB.
            if seq_bytes.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("File format error: sequence '{}' has length 0", record.id),
                ));
            }
            self.pos_array
                .push((self.current_pos, seq_bytes.len() as u32));

            self.writer.write_all(&[0xFF])?;
            self.writer.write_all(seq_bytes)?;
            self.writer.write_all(&[0xFF])?;
            self.writer.write_all(record.id.as_bytes())?;
            self.writer.write_all(&[0])?;
            self.current_pos += seq_bytes.len() as u64 + record.id.len() as u64 + 3;
            self.sequences += 1;
            self.letters += seq_bytes.len() as u64;

            self.hash = murmurhash3_x64_128(seq_bytes, &self.hash);
            self.hash = murmurhash3_x64_128(record.id.as_bytes(), &self.hash);

            if self.options.taxonmap.is_some() {
                for accession in accession_from_title(
                    &record.id,
                    !self.options.no_parse_seqids,
                    &mut self.accession_stats,
                ) {
                    self.accessions
                        .try_push((accession, self.sequences - 1))
                        .map_err(|error| io::Error::other(error.to_string()))?;
                }
            }
        }
        records.clear();
        Ok(())
    }

    fn finish(mut self) -> io::Result<BuildStats> {
        let pos_array_offset = self.current_pos;

        // Position array (u64 pos + u32 seq_len + u32 padding per entry).
        for &(pos, seq_len) in &self.pos_array {
            self.writer.write_all(&pos.to_le_bytes())?;
            self.writer.write_all(&seq_len.to_le_bytes())?;
            self.writer.write_all(&0u32.to_le_bytes())?;
        }
        // Sentinel: (current end-of-seqdata, 0, 0).
        self.writer.write_all(&self.current_pos.to_le_bytes())?;
        self.writer.write_all(&0u32.to_le_bytes())?;
        self.writer.write_all(&0u32.to_le_bytes())?;

        let mut header2 = ReferenceHeader2::new();
        header2.hash = self.hash;

        // Upstream calls the process-wide `taxonomy.init()` before building
        // the accession mapping. Preserve both lifetime and error order: the
        // sparse name vector remains resident through mapping/node output.
        let mut taxonomy = Taxonomy::default();
        if let Some(names_path) = self.options.taxonnames.as_deref() {
            taxonomy
                .init(names_path)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        }

        if let Some(mapping_path) = self.options.taxonmap.as_deref() {
            header2.taxon_array_offset = self.writer.stream_position()?;
            let build = TaxonList::write_from_file(
                &mut self.accessions,
                mapping_path,
                self.sequences,
                self.options.no_parse_seqids,
                &mut self.writer,
            )
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
            header2.taxon_array_size = build.bytes_written;
        }

        if let Some(nodes_path) = self.options.taxonnodes.as_deref() {
            let nodes = TaxonomyNodes::from_nodes_dmp(nodes_path)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
            header2.taxon_nodes_offset = self.writer.stream_position()?;
            nodes.write_to(&mut self.writer)?;
        }

        if self.options.taxonnames.is_some() {
            header2.taxon_names_offset = self.writer.stream_position()?;
            write_string_vector(&mut self.writer, &taxonomy.name)?;
        }

        self.writer.flush()?;
        drop(self.writer);

        // Upstream seeks to zero and rewrites both headers once all counts and
        // the chained hash are known.
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .open(&self.output_path)?;
        let mut header = ReferenceHeader::new();
        header.db_version = match self.seq_type {
            SequenceType::AminoAcid => CURRENT_DB_VERSION_PROT,
            SequenceType::Nucleotide => 4,
        };
        header.sequences = self.sequences;
        header.letters = self.letters;
        header.pos_array_offset = pos_array_offset;
        file.seek(SeekFrom::Start(0))?;
        header.write_to(&mut file)?;
        header2.write_record_to(&mut file)?;
        file.flush()?;

        Ok(BuildStats {
            sequences: self.sequences,
            letters: self.letters,
        })
    }
}

fn write_string_vector<W: Write>(writer: &mut W, values: &[String]) -> io::Result<()> {
    writer.write_all(&(values.len() as u32).to_le_bytes())?;
    for value in values {
        writer.write_all(value.as_bytes())?;
        writer.write_all(&[0])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::dmnd::{DatabaseFile, DatabaseOpenOptions, MAGIC_NUMBER};

    fn temp_path(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond_rs_{label}_{}_{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ))
    }

    #[test]
    fn test_build_stats() {
        let records = vec![
            FastaRecord {
                id: "seq1".to_string(),
                sequence: vec![0, 1, 2, 3],
            },
            FastaRecord {
                id: "seq2".to_string(),
                sequence: vec![4, 5, 6],
            },
        ];
        let original = records.clone();

        let tmp = std::env::temp_dir().join("test_diamond.dmnd");
        let stats = build_db_from_records(&records, tmp.to_str().unwrap(), SequenceType::AminoAcid)
            .unwrap();

        assert_eq!(stats.sequences, 2);
        assert_eq!(stats.letters, 7);
        for (record, original) in records.iter().zip(&original) {
            assert_eq!(record.id, original.id);
            assert_eq!(record.sequence, original.sequence);
        }

        // Verify it starts with the magic number
        let data = std::fs::read(&tmp).unwrap();
        assert!(data.len() > 8);
        let magic = u64::from_le_bytes(data[0..8].try_into().unwrap());
        assert_eq!(magic, MAGIC_NUMBER);

        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn streamed_blocks_produce_the_same_database_as_owned_records() {
        // The limit deliberately falls inside this small, varied input.  This
        // exercises the production block loop without manufacturing a huge or
        // duplicated fixture.
        let fasta_path = temp_path("streamed_blocks.faa");
        let streamed_path = temp_path("streamed_blocks.dmnd");
        let owned_path = temp_path("owned_records.dmnd");
        std::fs::write(
            &fasta_path,
            b">alpha first\nARND\n>beta\nCQE\n>gamma\nGHILK\n",
        )
        .unwrap();

        let stats = build_db_with_block_limit(
            &[fasta_path.to_str().unwrap()],
            streamed_path.to_str().unwrap(),
            SequenceType::AminoAcid,
            5,
            &BuildOptions::default(),
        )
        .unwrap();
        let records = fasta::read_fasta_file(&fasta_path, SequenceType::AminoAcid).unwrap();
        let owned_stats = build_db_from_records(
            &records,
            owned_path.to_str().unwrap(),
            SequenceType::AminoAcid,
        )
        .unwrap();

        assert_eq!(stats.sequences, 3);
        assert_eq!(stats.letters, 12);
        assert_eq!(owned_stats.sequences, stats.sequences);
        assert_eq!(owned_stats.letters, stats.letters);
        assert_eq!(
            std::fs::read(&streamed_path).unwrap(),
            std::fs::read(&owned_path).unwrap()
        );

        let _ = std::fs::remove_file(fasta_path);
        let _ = std::fs::remove_file(streamed_path);
        let _ = std::fs::remove_file(owned_path);
    }

    #[test]
    fn streaming_failure_removes_incomplete_database() {
        let fasta_path = temp_path("invalid_stream.faa");
        let output_path = temp_path("invalid_stream.dmnd");
        std::fs::write(&fasta_path, b">valid\nARND\n>invalid\nARN?\n").unwrap();

        let error = build_db_with_block_limit(
            &[fasta_path.to_str().unwrap()],
            output_path.to_str().unwrap(),
            SequenceType::AminoAcid,
            2,
            &BuildOptions::default(),
        )
        .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(!output_path.exists());

        let _ = std::fs::remove_file(fasta_path);
    }

    #[test]
    fn taxonomy_metadata_is_serialized_and_read_back_exactly() {
        let fasta_path = temp_path("taxonomy.faa");
        let mapping_path = temp_path("taxonomy.tsv");
        let nodes_path = temp_path("nodes.dmp");
        let names_path = temp_path("names.dmp");
        let output_path = temp_path("taxonomy.dmnd");
        std::fs::write(
            &fasta_path,
            b">sp|P12345.2|alpha >UniRef90_Q999.1\nMKTAYIAK\n>tr|ACC3.4|beta\nVLSPADKT\n>unmapped.7 gamma\nMPEPTIDE\n",
        )
        .unwrap();
        std::fs::write(
            &mapping_path,
            b"accession.version\ttaxid\nP12345\t10\nQ999\t20\nACC3\t30\nunused\t999\n",
        )
        .unwrap();
        std::fs::write(
            &nodes_path,
            b"1\t|\t1\t|\tno rank\t|\t\n10\t|\t1\t|\tspecies\t|\t\n20\t|\t1\t|\tgenus\t|\t\n30\t|\t20\t|\tstrain\t|\t\n",
        )
        .unwrap();
        std::fs::write(
            &names_path,
            b"1\t|\troot\t|\t\t|\tscientific name\t|\n10\t|\tAlpha species\t|\t\t|\tscientific name\t|\n10\t|\talpha\t|\t\t|\tcommon name\t|\n20\t|\tBeta genus\t|\t\t|\tscientific name\t|\n30\t|\tBeta strain\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();

        let stats = build_db_with_options(
            &[fasta_path.to_str().unwrap()],
            output_path.to_str().unwrap(),
            SequenceType::AminoAcid,
            &BuildOptions {
                taxonmap: Some(mapping_path.to_string_lossy().into_owned()),
                taxonnodes: Some(nodes_path.to_string_lossy().into_owned()),
                taxonnames: Some(names_path.to_string_lossy().into_owned()),
                no_parse_seqids: false,
            },
        )
        .unwrap();
        assert_eq!(stats.sequences, 3);
        assert_eq!(stats.letters, 24);

        let db = DatabaseFile::open_with_options(
            &output_path,
            DatabaseOpenOptions {
                load_taxon_mapping: true,
                load_taxon_nodes: true,
                load_taxon_scientific_names: true,
                ..DatabaseOpenOptions::default()
            },
        )
        .unwrap();
        assert!(db.header2.taxon_array_offset < db.header2.taxon_nodes_offset);
        assert!(db.header2.taxon_nodes_offset < db.header2.taxon_names_offset);
        assert_eq!(db.taxids(0).unwrap(), vec![10, 20]);
        assert_eq!(db.taxids(1).unwrap(), vec![30]);
        assert_eq!(db.taxids(2).unwrap(), Vec::<i32>::new());
        assert_eq!(db.get_parent(30).unwrap(), 20);
        assert_eq!(
            db.rank(10).unwrap(),
            crate::data::taxonomy::Rank::Species as i32
        );
        assert_eq!(db.taxon_scientific_name(10), "Alpha species");
        assert_eq!(db.taxon_scientific_name(11), "11");

        for path in [
            fasta_path,
            mapping_path,
            nodes_path,
            names_path,
            output_path,
        ] {
            let _ = std::fs::remove_file(path);
        }
    }

    #[test]
    fn taxonomy_failure_removes_incomplete_database() {
        let fasta_path = temp_path("taxonomy_error.faa");
        let mapping_path = temp_path("taxonomy_error.tsv");
        let output_path = temp_path("taxonomy_error.dmnd");
        std::fs::write(&fasta_path, b">P1\nMPEPTIDE\n").unwrap();
        std::fs::write(&mapping_path, b"bad\theader\nP1\t10\n").unwrap();

        let error = build_db_with_options(
            &[fasta_path.to_str().unwrap()],
            output_path.to_str().unwrap(),
            SequenceType::AminoAcid,
            &BuildOptions {
                taxonmap: Some(mapping_path.to_string_lossy().into_owned()),
                ..BuildOptions::default()
            },
        )
        .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(!output_path.exists());

        let _ = std::fs::remove_file(fasta_path);
        let _ = std::fs::remove_file(mapping_path);
    }
}
