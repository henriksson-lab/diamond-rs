use std::io::{self, Read, Write};

/// Magic number identifying a DAA (DIAMOND Alignment Archive) file.
pub const DAA_MAGIC_NUMBER: u64 = 0x3c0e53476d3ee36b;

/// Current DAA format version.
pub const DAA_VERSION: u64 = 1;

/// DAA block types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum BlockType {
    Empty = 0,
    Alignments = 1,
    RefNames = 2,
    RefLengths = 3,
}

/// DAA file header 1 (16 bytes).
#[derive(Debug, Clone, Default)]
pub struct DaaHeader1 {
    pub magic_number: u64,
    pub version: u64,
}

impl DaaHeader1 {
    pub fn new() -> Self {
        DaaHeader1 {
            magic_number: DAA_MAGIC_NUMBER,
            version: DAA_VERSION,
        }
    }

    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf = [0u8; 16];
        reader.read_exact(&mut buf)?;
        Ok(DaaHeader1 {
            magic_number: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            version: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
        })
    }

    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.magic_number.to_le_bytes())?;
        writer.write_all(&self.version.to_le_bytes())?;
        Ok(())
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.magic_number != DAA_MAGIC_NUMBER {
            return Err("File is not a DAA file.".into());
        }
        Ok(())
    }
}

/// DAA file header 2 (variable size).
#[derive(Debug, Clone)]
pub struct DaaHeader2 {
    pub diamond_build: u64,
    pub db_seqs: u64,
    pub db_seqs_used: u64,
    pub db_letters: u64,
    pub flags: u64,
    pub query_records: u64,
    pub mode: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub reward: i32,
    pub penalty: i32,
    pub reserved1: i32,
    pub reserved2: i32,
    pub reserved3: i32,
    pub k: f64,
    pub lambda: f64,
    pub evalue: f64,
    pub reserved5: f64,
    pub score_matrix: [u8; 16],
    pub block_size: [u64; 256],
    pub block_type: [u8; 256],
}

impl Default for DaaHeader2 {
    fn default() -> Self {
        DaaHeader2 {
            diamond_build: 0,
            db_seqs: 0,
            db_seqs_used: 0,
            db_letters: 0,
            flags: 0,
            query_records: 0,
            mode: 0,
            gap_open: 0,
            gap_extend: 0,
            reward: 0,
            penalty: 0,
            reserved1: 0,
            reserved2: 0,
            reserved3: 0,
            k: 0.0,
            lambda: 0.0,
            evalue: 0.0,
            reserved5: 0.0,
            score_matrix: [0; 16],
            block_size: [0; 256],
            block_type: [0; 256],
        }
    }

}

impl DaaHeader2 {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut h = Self::default();

        let mut buf8 = [0u8; 8];
        let mut buf4 = [0u8; 4];

        macro_rules! read_u64 {
            ($field:expr) => {{
                reader.read_exact(&mut buf8)?;
                $field = u64::from_le_bytes(buf8);
            }};
        }
        macro_rules! read_i32 {
            ($field:expr) => {{
                reader.read_exact(&mut buf4)?;
                $field = i32::from_le_bytes(buf4);
            }};
        }
        macro_rules! read_f64 {
            ($field:expr) => {{
                reader.read_exact(&mut buf8)?;
                $field = f64::from_le_bytes(buf8);
            }};
        }

        read_u64!(h.diamond_build);
        read_u64!(h.db_seqs);
        read_u64!(h.db_seqs_used);
        read_u64!(h.db_letters);
        read_u64!(h.flags);
        read_u64!(h.query_records);
        read_i32!(h.mode);
        read_i32!(h.gap_open);
        read_i32!(h.gap_extend);
        read_i32!(h.reward);
        read_i32!(h.penalty);
        read_i32!(h.reserved1);
        read_i32!(h.reserved2);
        read_i32!(h.reserved3);
        read_f64!(h.k);
        read_f64!(h.lambda);
        read_f64!(h.evalue);
        read_f64!(h.reserved5);
        reader.read_exact(&mut h.score_matrix)?;

        // Read block_size array (256 x u64)
        for i in 0..256 {
            reader.read_exact(&mut buf8)?;
            h.block_size[i] = u64::from_le_bytes(buf8);
        }

        // Read block_type array (256 x u8)
        reader.read_exact(&mut h.block_type)?;

        Ok(h)
    }

    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.diamond_build.to_le_bytes())?;
        writer.write_all(&self.db_seqs.to_le_bytes())?;
        writer.write_all(&self.db_seqs_used.to_le_bytes())?;
        writer.write_all(&self.db_letters.to_le_bytes())?;
        writer.write_all(&self.flags.to_le_bytes())?;
        writer.write_all(&self.query_records.to_le_bytes())?;
        writer.write_all(&self.mode.to_le_bytes())?;
        writer.write_all(&self.gap_open.to_le_bytes())?;
        writer.write_all(&self.gap_extend.to_le_bytes())?;
        writer.write_all(&self.reward.to_le_bytes())?;
        writer.write_all(&self.penalty.to_le_bytes())?;
        writer.write_all(&self.reserved1.to_le_bytes())?;
        writer.write_all(&self.reserved2.to_le_bytes())?;
        writer.write_all(&self.reserved3.to_le_bytes())?;
        writer.write_all(&self.k.to_le_bytes())?;
        writer.write_all(&self.lambda.to_le_bytes())?;
        writer.write_all(&self.evalue.to_le_bytes())?;
        writer.write_all(&self.reserved5.to_le_bytes())?;
        writer.write_all(&self.score_matrix)?;
        for &size in &self.block_size {
            writer.write_all(&size.to_le_bytes())?;
        }
        writer.write_all(&self.block_type)?;
        Ok(())
    }

    /// Get the score matrix name as a string.
    pub fn matrix_name(&self) -> String {
        let end = self
            .score_matrix
            .iter()
            .position(|&b| b == 0)
            .unwrap_or(16);
        String::from_utf8_lossy(&self.score_matrix[..end]).to_string()
    }

    /// Set the score matrix name.
    pub fn set_matrix_name(&mut self, name: &str) {
        self.score_matrix = [0; 16];
        let bytes = name.as_bytes();
        let len = bytes.len().min(15);
        self.score_matrix[..len].copy_from_slice(&bytes[..len]);
    }
}

/// Packed integer reading from a flag byte.
///
/// The flag byte encodes the byte width of subsequent packed integers.
/// Bits [1:0] = score width, [3:2] = query_begin width, [5:4] = subject_begin width.
/// Width encoding: 0 = 1 byte, 1 = 2 bytes, 2 = 4 bytes.
pub fn read_packed_int<R: Read>(reader: &mut R, width_bits: u8) -> io::Result<u32> {
    match width_bits & 0x3 {
        0 => {
            let mut buf = [0u8; 1];
            reader.read_exact(&mut buf)?;
            Ok(buf[0] as u32)
        }
        1 => {
            let mut buf = [0u8; 2];
            reader.read_exact(&mut buf)?;
            Ok(u16::from_le_bytes(buf) as u32)
        }
        2 => {
            let mut buf = [0u8; 4];
            reader.read_exact(&mut buf)?;
            Ok(u32::from_le_bytes(buf))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid packed int width",
        )),
    }
}

/// Write a packed integer using the minimum number of bytes.
pub fn write_packed_int<W: Write>(writer: &mut W, val: u32) -> io::Result<u8> {
    if val <= 0xFF {
        writer.write_all(&[val as u8])?;
        Ok(0)
    } else if val <= 0xFFFF {
        writer.write_all(&(val as u16).to_le_bytes())?;
        Ok(1)
    } else {
        writer.write_all(&val.to_le_bytes())?;
        Ok(2)
    }
}

/// Compute the flag byte for a match record.
pub fn compute_flag(score: u32, query_begin: u32, subject_begin: u32, reverse_strand: bool) -> u8 {
    let score_bits = if score <= 0xFF {
        0u8
    } else if score <= 0xFFFF {
        1
    } else {
        2
    };
    let qb_bits = if query_begin <= 0xFF {
        0u8
    } else if query_begin <= 0xFFFF {
        1
    } else {
        2
    };
    let sb_bits = if subject_begin <= 0xFF {
        0u8
    } else if subject_begin <= 0xFFFF {
        1
    } else {
        2
    };
    let strand_bit = if reverse_strand { 1u8 << 6 } else { 0 };
    score_bits | (qb_bits << 2) | (sb_bits << 4) | strand_bit
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_header1_roundtrip() {
        let h = DaaHeader1::new();
        let mut buf = Vec::new();
        h.write_to(&mut buf).unwrap();
        assert_eq!(buf.len(), 16);

        let h2 = DaaHeader1::read_from(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(h2.magic_number, DAA_MAGIC_NUMBER);
        assert_eq!(h2.version, DAA_VERSION);
    }

    #[test]
    fn test_header2_matrix_name() {
        let mut h = DaaHeader2::new();
        h.set_matrix_name("BLOSUM62");
        assert_eq!(h.matrix_name(), "BLOSUM62");
    }

    #[test]
    fn test_packed_int_roundtrip() {
        for val in [0u32, 1, 127, 255, 256, 65535, 65536, 1000000] {
            let mut buf = Vec::new();
            let width = write_packed_int(&mut buf, val).unwrap();
            let read_back = read_packed_int(&mut Cursor::new(&buf), width).unwrap();
            assert_eq!(read_back, val, "Failed for value {}", val);
        }
    }

    #[test]
    fn test_compute_flag() {
        let flag = compute_flag(50, 10, 20, false);
        assert_eq!(flag & 0x3, 0); // score fits in 1 byte
        assert_eq!((flag >> 2) & 0x3, 0); // query_begin fits in 1 byte
        assert_eq!((flag >> 4) & 0x3, 0); // subject_begin fits in 1 byte
        assert_eq!(flag & 0x40, 0); // forward strand

        let flag2 = compute_flag(300, 10, 70000, true);
        assert_eq!(flag2 & 0x3, 1); // score needs 2 bytes
        assert_eq!((flag2 >> 4) & 0x3, 2); // subject_begin needs 4 bytes
        assert_ne!(flag2 & 0x40, 0); // reverse strand
    }
}
