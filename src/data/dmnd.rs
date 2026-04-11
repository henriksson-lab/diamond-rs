use std::io::{self, Read, Write};

/// Magic number identifying a DIAMOND database file.
pub const MAGIC_NUMBER: u64 = 0x24af8a415ee186d;

/// Current database version for protein databases.
pub const CURRENT_DB_VERSION_PROT: u32 = 3;

/// Current database version for nucleotide databases.
pub const CURRENT_DB_VERSION_NUCL: u32 = 4;

/// Minimum compatible database version.
pub const MIN_DB_VERSION: u32 = 2;

/// Minimum compatible build version.
pub const MIN_BUILD_VERSION: u32 = 74;

/// DIAMOND database file header (first 40 bytes).
#[derive(Debug, Clone)]
pub struct ReferenceHeader {
    pub magic_number: u64,
    pub build: u32,
    pub db_version: u32,
    pub sequences: u64,
    pub letters: u64,
    pub pos_array_offset: u64,
}

impl Default for ReferenceHeader {
    fn default() -> Self {
        ReferenceHeader {
            magic_number: MAGIC_NUMBER,
            build: 0,
            db_version: CURRENT_DB_VERSION_PROT,
            sequences: 0,
            letters: 0,
            pos_array_offset: 0,
        }
    }
}

impl ReferenceHeader {
    pub fn new() -> Self {
        Self::default()
    }

    /// Read header from a reader.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf = [0u8; 40];
        reader.read_exact(&mut buf)?;
        Ok(ReferenceHeader {
            magic_number: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            build: u32::from_le_bytes(buf[8..12].try_into().unwrap()),
            db_version: u32::from_le_bytes(buf[12..16].try_into().unwrap()),
            sequences: u64::from_le_bytes(buf[16..24].try_into().unwrap()),
            letters: u64::from_le_bytes(buf[24..32].try_into().unwrap()),
            pos_array_offset: u64::from_le_bytes(buf[32..40].try_into().unwrap()),
        })
    }

    /// Write header to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.magic_number.to_le_bytes())?;
        writer.write_all(&self.build.to_le_bytes())?;
        writer.write_all(&self.db_version.to_le_bytes())?;
        writer.write_all(&self.sequences.to_le_bytes())?;
        writer.write_all(&self.letters.to_le_bytes())?;
        writer.write_all(&self.pos_array_offset.to_le_bytes())?;
        Ok(())
    }

    /// Validate this header is a valid DIAMOND database.
    pub fn validate(&self) -> Result<(), String> {
        if self.magic_number != MAGIC_NUMBER {
            return Err("Database file is not a DIAMOND database.".into());
        }
        if self.db_version < MIN_DB_VERSION {
            return Err(format!(
                "Database version {} is not supported (minimum: {})",
                self.db_version, MIN_DB_VERSION
            ));
        }
        Ok(())
    }
}

/// Extended DIAMOND database header.
#[derive(Debug, Clone, Default)]
pub struct ReferenceHeader2 {
    pub hash: [u8; 16],
    pub taxon_array_offset: u64,
    pub taxon_array_size: u64,
    pub taxon_nodes_offset: u64,
    pub taxon_names_offset: u64,
}

impl ReferenceHeader2 {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut hash = [0u8; 16];
        reader.read_exact(&mut hash)?;
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf)?;
        Ok(ReferenceHeader2 {
            hash,
            taxon_array_offset: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            taxon_array_size: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
            taxon_nodes_offset: u64::from_le_bytes(buf[16..24].try_into().unwrap()),
            taxon_names_offset: u64::from_le_bytes(buf[24..32].try_into().unwrap()),
        })
    }

    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.hash)?;
        writer.write_all(&self.taxon_array_offset.to_le_bytes())?;
        writer.write_all(&self.taxon_array_size.to_le_bytes())?;
        writer.write_all(&self.taxon_nodes_offset.to_le_bytes())?;
        writer.write_all(&self.taxon_names_offset.to_le_bytes())?;
        Ok(())
    }
}

/// Sequence info entry in the position array (trailer).
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct SeqInfo {
    pub pos: u64,
    pub seq_len: u32,
    pub padding: u32,
}

impl SeqInfo {
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf = [0u8; 16];
        reader.read_exact(&mut buf)?;
        Ok(SeqInfo {
            pos: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            seq_len: u32::from_le_bytes(buf[8..12].try_into().unwrap()),
            padding: u32::from_le_bytes(buf[12..16].try_into().unwrap()),
        })
    }
}

/// Check if a file is a DIAMOND database.
pub fn is_diamond_db<R: Read>(reader: &mut R) -> bool {
    let mut buf = [0u8; 8];
    if reader.read_exact(&mut buf).is_err() {
        return false;
    }
    u64::from_le_bytes(buf) == MAGIC_NUMBER
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_header_roundtrip() {
        let mut h = ReferenceHeader::new();
        h.sequences = 42;
        h.letters = 12345;
        h.pos_array_offset = 999;

        let mut buf = Vec::new();
        h.write_to(&mut buf).unwrap();
        assert_eq!(buf.len(), 40);

        let h2 = ReferenceHeader::read_from(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(h2.magic_number, MAGIC_NUMBER);
        assert_eq!(h2.sequences, 42);
        assert_eq!(h2.letters, 12345);
        assert_eq!(h2.pos_array_offset, 999);
    }

    #[test]
    fn test_header_validation() {
        let h = ReferenceHeader::new();
        assert!(h.validate().is_ok());

        let mut bad = ReferenceHeader::new();
        bad.magic_number = 0;
        assert!(bad.validate().is_err());
    }

    #[test]
    fn test_is_diamond_db() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&MAGIC_NUMBER.to_le_bytes());
        assert!(is_diamond_db(&mut Cursor::new(&buf)));

        let bad_buf = vec![0u8; 8];
        assert!(!is_diamond_db(&mut Cursor::new(&bad_buf)));
    }

    #[test]
    fn test_read_real_db() {
        let db_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/data.dmnd"
        );
        if let Ok(mut f) = std::fs::File::open(db_path) {
            let h = ReferenceHeader::read_from(&mut f).unwrap();
            assert_eq!(h.magic_number, MAGIC_NUMBER);
            assert!(h.validate().is_ok());
            assert!(h.sequences > 0);
        }
    }
}
