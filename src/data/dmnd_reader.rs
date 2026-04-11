use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use crate::basic::value::Letter;
use super::dmnd::ReferenceHeader;
use super::fasta::FastaRecord;

/// Read all sequences from a DIAMOND database file.
///
/// Returns the sequences as FastaRecord objects with their IDs and Letter-encoded data.
pub fn read_dmnd(path: &Path) -> io::Result<(ReferenceHeader, Vec<FastaRecord>)> {
    let mut file = BufReader::new(File::open(path)?);

    // Read header1
    let header = ReferenceHeader::read_from(&mut file)?;
    header
        .validate()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    // Read header2 size prefix and skip header2
    let mut buf8 = [0u8; 8];
    file.read_exact(&mut buf8)?;
    let h2_size = u64::from_le_bytes(buf8);
    file.seek(SeekFrom::Current(h2_size as i64))?;

    // Read sequence data up to pos_array_offset
    let seq_data_start = file.stream_position()?;
    let seq_data_len = header.pos_array_offset - seq_data_start;
    let mut seq_data = vec![0u8; seq_data_len as usize];
    file.read_exact(&mut seq_data)?;

    // Parse sequences from the data block
    // Format: [0xFF][sequence_bytes][0xFF][id\0]
    let mut records = Vec::with_capacity(header.sequences as usize);
    let mut i = 0usize;

    while i < seq_data.len() {
        if seq_data[i] == 0xFF {
            i += 1;
            // Read sequence until next 0xFF
            let seq_start = i;
            while i < seq_data.len() && seq_data[i] != 0xFF {
                i += 1;
            }
            let seq_bytes = &seq_data[seq_start..i];
            i += 1; // skip closing 0xFF

            // Read ID until null
            let id_start = i;
            while i < seq_data.len() && seq_data[i] != 0x00 {
                i += 1;
            }
            let id = String::from_utf8_lossy(&seq_data[id_start..i]).to_string();
            i += 1; // skip null

            // Convert bytes to Letter values (they're already encoded)
            let sequence: Vec<Letter> = seq_bytes.iter().map(|&b| b as Letter).collect();

            records.push(FastaRecord { id, sequence });
        } else {
            i += 1;
        }
    }

    Ok((header, records))
}

/// Read a DIAMOND database, handling .dmnd extension.
pub fn read_dmnd_auto(path_str: &str) -> io::Result<(ReferenceHeader, Vec<FastaRecord>)> {
    let path = Path::new(path_str);
    if path.exists() {
        return read_dmnd(path);
    }
    let with_ext = path.with_extension("dmnd");
    if with_ext.exists() {
        return read_dmnd(&with_ext);
    }
    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!("Database file not found: {}", path_str),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_test_dmnd() {
        let db_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/data.dmnd"
        );
        let (header, records) = read_dmnd(Path::new(db_path)).unwrap();
        assert_eq!(header.sequences, 1);
        assert_eq!(header.letters, 426);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "d1ivsa4");
        assert_eq!(records[0].sequence.len(), 426);
    }

    #[test]
    fn test_read_galaxy_dmnd() {
        let db_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/galaxy/db.dmnd"
        );
        let (header, records) = read_dmnd(Path::new(db_path)).unwrap();
        assert!(header.sequences > 0);
        assert_eq!(records.len(), header.sequences as usize);
        for r in &records {
            assert!(!r.id.is_empty());
            assert!(!r.sequence.is_empty());
        }
    }

    #[test]
    fn test_read_auto_extension() {
        let db_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/data"
        );
        let result = read_dmnd_auto(db_path);
        assert!(result.is_ok());
        let (_, records) = result.unwrap();
        assert_eq!(records.len(), 1);
    }
}
