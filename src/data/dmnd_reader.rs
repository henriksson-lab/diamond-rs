use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use super::dmnd::{ReferenceHeader, SeqInfo};
use super::fasta::FastaRecord;
use crate::basic::value::Letter;

/// Read all sequences from a DIAMOND database file.
///
/// Returns the sequences as FastaRecord objects with their IDs and Letter-encoded data.
///
/// The reader walks the position array (SeqInfo entries) rather than byte-scanning
/// the data block for `0xFF`. C++ does the same (`dmnd.cpp:551`): for each SeqInfo
/// the sequence bytes are at `[pos+1, pos+1+seq_len)` and the id is the
/// null-terminated string at `pos+1+seq_len+1`. Byte-scanning works for the
/// current encoding (masked letters stay below 0x9F) but would break the moment
/// any encoding uses 0xFF in a payload, and it silently swallows corrupt data.
pub fn read_dmnd(path: &Path) -> io::Result<(ReferenceHeader, Vec<FastaRecord>)> {
    let mut file = BufReader::new(File::open(path)?);

    // Read header1
    let header = ReferenceHeader::read_from(&mut file)?;
    header
        .validate()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    // Read header2 size prefix and skip header2 — honor the prefix so we
    // correctly skip any future-extended header2 (e.g. C++ EXTRA builds add
    // a `db_type` int32_t making h2 52 bytes instead of 48).
    let mut buf8 = [0u8; 8];
    file.read_exact(&mut buf8)?;
    let h2_size = u64::from_le_bytes(buf8);
    file.seek(SeekFrom::Current(h2_size as i64))?;

    let seq_data_start = file.stream_position()?;
    if header.pos_array_offset < seq_data_start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "pos_array_offset is before sequence data — corrupt DMND header",
        ));
    }
    let seq_data_len = (header.pos_array_offset - seq_data_start) as usize;
    let mut seq_data = vec![0u8; seq_data_len];
    file.read_exact(&mut seq_data)?;

    // Read the position array. C++ stores `sequences + 1` entries (one sentinel
    // with `seq_len == 0` at the end). Use `SeqInfo` to extract the exact
    // sequence bounds; the bytes are `[pos+1, pos+1+seq_len)` and the id is
    // the C-string immediately following the closing `0xFF` framing byte.
    let mut records = Vec::with_capacity(header.sequences as usize);
    for _ in 0..header.sequences {
        let info = SeqInfo::read_from(&mut file)?;
        if info.pos < seq_data_start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SeqInfo.pos {} is before sequence data start {}",
                    info.pos, seq_data_start
                ),
            ));
        }
        let entry_off = (info.pos - seq_data_start) as usize;
        // Layout: [0xFF][seq_bytes; seq_len][0xFF][id; until 0x00][0x00]
        let seq_start = entry_off + 1;
        let seq_end = seq_start + info.seq_len as usize;
        if seq_end + 1 > seq_data_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SeqInfo bounds extend past sequence data — truncated DMND",
            ));
        }
        let seq_bytes = &seq_data[seq_start..seq_end];
        // Skip the closing 0xFF separator, then read the null-terminated id.
        let id_start = seq_end + 1;
        let mut id_end = id_start;
        while id_end < seq_data_len && seq_data[id_end] != 0x00 {
            id_end += 1;
        }
        let id = String::from_utf8_lossy(&seq_data[id_start..id_end]).to_string();
        // Bytes are already in Letter encoding; just bit-cast u8→i8.
        let sequence: Vec<Letter> = seq_bytes.iter().map(|&b| b as Letter).collect();
        records.push(FastaRecord { id, sequence });
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
        let db_path = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/data.dmnd");
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
        let db_path = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/data");
        let result = read_dmnd_auto(db_path);
        assert!(result.is_ok());
        let (_, records) = result.unwrap();
        assert_eq!(records.len(), 1);
    }
}
