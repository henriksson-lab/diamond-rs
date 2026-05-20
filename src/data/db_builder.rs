use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use super::dmnd::{ReferenceHeader, ReferenceHeader2, CURRENT_DB_VERSION_PROT};
use super::fasta::{self, FastaRecord};
use crate::basic::value::SequenceType;
use crate::masking::tantan::mask_tantan;
use crate::util::hash::murmurhash3_x64_128;

/// Total size of the two on-disk reference headers (40 + 8 size prefix + 48 data).
/// Matches `out->tell()` after `*out << header << header2` in diamond/src/data/dmnd/dmnd.cpp.
const HEADER_TOTAL: u64 = 40 + 8 + 48;

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
    let total_sequences = records.len() as u64;

    // Apply tantan bit-masking to each amino-acid sequence (matches C++
    // mask_seqs(block->seqs(), Masking::get(), /*hard_mask=*/false, MaskingAlgo::SEG)
    // in diamond/src/data/dmnd/dmnd.cpp — with hard_mask=false the worker takes
    // the mask_bit branch which is tantan in mode 2, setting bit 7 on masked
    // positions). For nucleotides, no masking is applied.
    let mut masked: Vec<Vec<u8>> = Vec::with_capacity(records.len());
    for r in records {
        let mut seq: Vec<i8> = r.sequence.clone();
        if seq_type == SequenceType::AminoAcid {
            mask_tantan(&mut seq);
        }
        masked.push(seq.iter().map(|&l| l as u8).collect());
    }

    let total_letters: u64 = masked.iter().map(|s| s.len() as u64).sum();

    let file = std::fs::File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    // Write header1 (final pos_array_offset patched in at end).
    let mut header = ReferenceHeader::new();
    header.sequences = total_sequences;
    header.letters = total_letters;
    header.db_version = match seq_type {
        SequenceType::AminoAcid => CURRENT_DB_VERSION_PROT,
        SequenceType::Nucleotide => 4,
    };
    header.write_to(&mut writer)?;

    // Header2 with size prefix; hash is patched in at end.
    let h2_data_size = 48u64; // 16 hash + 4*8 offsets
    writer.write_all(&h2_data_size.to_le_bytes())?;
    let header2 = ReferenceHeader2::new();
    header2.write_to(&mut writer)?;

    // Sequences. C++ stores absolute file offsets in the pos array, starting at
    // out->tell() after the headers (= HEADER_TOTAL). Chain MurmurHash3_x64_128
    // over (masked_seq, id) per sequence into header2.hash (dmnd.cpp:328-329).
    let mut pos_array: Vec<(u64, u32)> = Vec::with_capacity(records.len());
    let mut current_pos = HEADER_TOTAL;
    let mut hash = [0u8; 16];

    for (record, seq_bytes) in records.iter().zip(masked.iter()) {
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
        pos_array.push((current_pos, seq_bytes.len() as u32));

        writer.write_all(&[0xFF])?;
        writer.write_all(seq_bytes)?;
        writer.write_all(&[0xFF])?;
        writer.write_all(record.id.as_bytes())?;
        writer.write_all(&[0])?;
        current_pos += seq_bytes.len() as u64 + record.id.len() as u64 + 3;

        hash = murmurhash3_x64_128(seq_bytes, &hash);
        hash = murmurhash3_x64_128(record.id.as_bytes(), &hash);
    }

    let pos_array_offset = current_pos;

    // Position array (u64 pos + u32 seq_len + u32 padding per entry).
    for &(pos, seq_len) in &pos_array {
        writer.write_all(&pos.to_le_bytes())?;
        writer.write_all(&seq_len.to_le_bytes())?;
        writer.write_all(&0u32.to_le_bytes())?;
    }
    // Sentinel: (current end-of-seqdata, 0, 0).
    writer.write_all(&current_pos.to_le_bytes())?;
    writer.write_all(&0u32.to_le_bytes())?;
    writer.write_all(&0u32.to_le_bytes())?;

    writer.flush()?;
    drop(writer);

    // Patch pos_array_offset (u64 at offset 32) and header2.hash (16 bytes at
    // offset 40 + 8 = 48). Matches C++ which seeks to 0 and rewrites both headers.
    let mut file = std::fs::OpenOptions::new().write(true).open(output_path)?;
    file.seek(SeekFrom::Start(32))?;
    file.write_all(&pos_array_offset.to_le_bytes())?;
    file.seek(SeekFrom::Start(40 + 8))?;
    file.write_all(&hash)?;
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
        let stats = build_db_from_records(&records, tmp.to_str().unwrap(), SequenceType::AminoAcid)
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
