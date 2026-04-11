use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use super::dmnd::{ReferenceHeader, ReferenceHeader2, CURRENT_DB_VERSION_PROT};
use super::fasta::{self, FastaRecord};
use crate::basic::value::SequenceType;

/// Build a DIAMOND database from FASTA input files.
///
/// This is the Rust implementation of the `makedb` command.
pub fn build_db(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    let mut sequences: Vec<FastaRecord> = Vec::new();

    for input_file in input_files {
        let path = Path::new(input_file);
        let records = fasta::read_fasta_file(path, seq_type)?;
        sequences.extend(records);
    }

    write_dmnd(output_path, &sequences, seq_type)
}

/// Build a DIAMOND database from in-memory sequences.
pub fn build_db_from_records(
    records: &[FastaRecord],
    output_path: &str,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    write_dmnd(output_path, records, seq_type)
}

/// Statistics from database building.
#[derive(Debug, Clone)]
pub struct BuildStats {
    pub sequences: u64,
    pub letters: u64,
}

fn write_dmnd(
    output_path: &str,
    records: &[FastaRecord],
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    let file = std::fs::File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    // Calculate total letters
    let total_sequences = records.len() as u64;
    let total_letters: u64 = records.iter().map(|r| r.sequence.len() as u64).sum();

    // Write header1
    let mut header = ReferenceHeader::new();
    header.sequences = total_sequences;
    header.letters = total_letters;
    header.db_version = match seq_type {
        SequenceType::AminoAcid => CURRENT_DB_VERSION_PROT,
        SequenceType::Nucleotide => 4,
    };
    header.write_to(&mut writer)?;

    // Write header2 with size prefix (matching C++ serialization format)
    let header2 = ReferenceHeader2::new();
    let h2_data_size = 48u64; // 16 hash + 4*8 offsets
    writer.write_all(&h2_data_size.to_le_bytes())?; // size prefix
    header2.write_to(&mut writer)?;

    // Write sequences
    // Format: [0xFF][sequence_data][0xFF][id\0]
    let mut pos_array: Vec<(u64, u32)> = Vec::with_capacity(records.len());
    let mut current_pos = 0u64; // relative position within sequence data

    for record in records {
        pos_array.push((current_pos, record.sequence.len() as u32));

        writer.write_all(&[0xFF])?; // delimiter
        current_pos += 1;

        for &letter in &record.sequence {
            writer.write_all(&[letter as u8])?;
            current_pos += 1;
        }

        writer.write_all(&[0xFF])?; // delimiter
        current_pos += 1;

        // Write ID as null-terminated string
        writer.write_all(record.id.as_bytes())?;
        writer.write_all(&[0])?;
        current_pos += record.id.len() as u64 + 1;
    }

    // Write position array (SeqInfo entries: pos u64 + seq_len u32 + padding u32)
    for &(pos, seq_len) in &pos_array {
        writer.write_all(&pos.to_le_bytes())?;
        writer.write_all(&seq_len.to_le_bytes())?;
        writer.write_all(&0u32.to_le_bytes())?; // padding
    }

    // Write sentinel entry
    writer.write_all(&current_pos.to_le_bytes())?;
    writer.write_all(&0u32.to_le_bytes())?;
    writer.write_all(&0u32.to_le_bytes())?;

    writer.flush()?;

    // Seek back and update header with pos_array_offset
    // Header1 is 40 bytes, header2 size prefix is 8 bytes, header2 data is 48 bytes
    let header1_size = 40u64;
    let header2_total = 8u64 + 48u64; // size prefix + data
    let pos_array_offset = header1_size + header2_total + current_pos;
    drop(writer);

    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .open(output_path)?;
    // pos_array_offset is at offset 32 in header1
    file.seek(SeekFrom::Start(32))?;
    file.write_all(&pos_array_offset.to_le_bytes())?;
    file.flush()?;

    Ok(BuildStats {
        sequences: total_sequences,
        letters: total_letters,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::dmnd::MAGIC_NUMBER;

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

        let tmp = std::env::temp_dir().join("test_diamond.dmnd");
        let stats = build_db_from_records(
            &records,
            tmp.to_str().unwrap(),
            SequenceType::AminoAcid,
        )
        .unwrap();

        assert_eq!(stats.sequences, 2);
        assert_eq!(stats.letters, 7);

        // Verify it starts with the magic number
        let data = std::fs::read(&tmp).unwrap();
        assert!(data.len() > 8);
        let magic = u64::from_le_bytes(data[0..8].try_into().unwrap());
        assert_eq!(magic, MAGIC_NUMBER);

        let _ = std::fs::remove_file(&tmp);
    }
}
