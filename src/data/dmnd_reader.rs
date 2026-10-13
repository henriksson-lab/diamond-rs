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
    read_dmnd(&resolve_dmnd_path(path_str)?)
}

/// Visit database records in OID order without retaining the complete database.
///
/// DIAMOND's getseq path advances through the sequence-data stream and the
/// position array together. Keep separate buffered file cursors for those two
/// monotonically increasing regions so a record only remains resident until
/// the callback returns.
pub fn for_each_dmnd_record_auto<F>(path_str: &str, mut visit: F) -> io::Result<ReferenceHeader>
where
    F: FnMut(usize, &str, &[u8]) -> io::Result<()>,
{
    let path = resolve_dmnd_path(path_str)?;
    for_each_dmnd_record(&path, &mut visit)
}

/// Visit database records from an already-resolved path without retaining the
/// complete decoded record set.
pub fn for_each_dmnd_record<F>(path: &Path, mut visit: F) -> io::Result<ReferenceHeader>
where
    F: FnMut(usize, &str, &[u8]) -> io::Result<()>,
{
    let mut data = BufReader::new(File::open(&path)?);
    let mut positions = BufReader::new(File::open(path)?);

    let header = ReferenceHeader::read_from(&mut data)?;
    header
        .validate()
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;

    let mut buf8 = [0u8; 8];
    data.read_exact(&mut buf8)?;
    let h2_size = u64::from_le_bytes(buf8);
    data.seek(SeekFrom::Current(h2_size as i64))?;
    let data_start = data.stream_position()?;
    if header.pos_array_offset < data_start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "pos_array_offset is before sequence data — corrupt DMND header",
        ));
    }
    positions.seek(SeekFrom::Start(header.pos_array_offset))?;

    let mut sequence = Vec::new();
    let mut id = Vec::new();
    let mut data_position = data_start;
    for ordinal in 0..header.sequences as usize {
        let info = SeqInfo::read_from(&mut positions)?;
        if info.pos < data_start || info.pos >= header.pos_array_offset {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SeqInfo position is outside sequence data — corrupt DMND header",
            ));
        }
        if data_position != info.pos {
            data.seek(SeekFrom::Start(info.pos))?;
            data_position = info.pos;
        }

        let mut delimiter = [0u8; 1];
        data.read_exact(&mut delimiter)?;
        if delimiter[0] != 0xff {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid opening sequence delimiter in DIAMOND database",
            ));
        }
        sequence.resize(info.seq_len as usize, 0);
        data.read_exact(&mut sequence)?;
        data.read_exact(&mut delimiter)?;
        if delimiter[0] != 0xff {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid closing sequence delimiter in DIAMOND database",
            ));
        }

        id.clear();
        loop {
            let mut byte = [0u8; 1];
            data.read_exact(&mut byte)?;
            if byte[0] == 0 {
                break;
            }
            id.push(byte[0]);
        }
        let record_len = u64::from(info.seq_len)
            .checked_add(id.len() as u64)
            .and_then(|length| length.checked_add(3))
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Sequence record length overflow in DIAMOND database",
                )
            })?;
        data_position = data_position.checked_add(record_len).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Sequence record offset overflow in DIAMOND database",
            )
        })?;
        let id = String::from_utf8_lossy(&id);
        visit(ordinal, &id, &sequence)?;
    }
    Ok(header)
}

fn resolve_dmnd_path(path_str: &str) -> io::Result<std::path::PathBuf> {
    let path = Path::new(path_str);
    if path.exists() {
        return Ok(path.to_path_buf());
    }
    let with_ext = path.with_extension("dmnd");
    if with_ext.exists() {
        return Ok(with_ext);
    }
    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!("Database file not found: {}", path_str),
    ))
}
