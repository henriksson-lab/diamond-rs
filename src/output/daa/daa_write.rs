//! Translation of `diamond/src/output/daa/daa_write.cpp`.
//!
//! The low-level implementations predate the mirrored output hierarchy and
//! live in [`crate::data::daa`].  The functions below deliberately delegate
//! to those implementations so that there is one canonical DAA wire encoder.

use std::io::{self, Seek, Write};

use crate::align::hsp::Hsp;
use crate::basic::value::{Letter, SequenceType};
use crate::data::daa::{self, DaaFile, DaaHeader2};
use crate::output::intermediate::IntermediateRecord;

/// Database access used by the `SequenceFile` overload of C++ `finish_daa`.
///
/// Titles are returned as owned strings because Rust database backends may
/// decode them lazily.  As in C++, the title is written verbatim followed by a
/// NUL byte.
pub trait DaaSequenceFile {
    fn sequence_count(&self) -> u64;
    fn dict_size(&self) -> usize;
    fn dict_title(&self, index: usize) -> String;
    fn dict_len(&self, index: usize) -> u32;
}

/// Values read from C++ globals by the `SequenceFile` overload of
/// `finish_daa` (`score_matrix`, `config`, `align_mode`, and `statistics`).
#[derive(Debug, Clone, PartialEq)]
pub struct DaaRunMetadata {
    pub db_letters: u64,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub reward: i32,
    pub penalty: i32,
    pub k: f64,
    pub lambda: f64,
    pub max_evalue: f64,
    pub matrix: String,
    pub mode: u32,
    pub aligned_queries: u64,
}

/// C++ `init_daa(OutputFile&)`.
pub fn init_daa<W: Write>(writer: &mut W) -> io::Result<()> {
    daa::init_daa(writer)
}

/// C++ `write_daa_query_record(TextBuffer&, const char*, const Sequence&)`.
///
/// `input_sequence_type` is explicit because the C++ implementation obtains
/// it from the global `align_mode`.
pub fn write_daa_query_record(
    buf: &mut Vec<u8>,
    query_name: &str,
    query: &[Letter],
    input_sequence_type: SequenceType,
) -> usize {
    daa::write_daa_query_record(buf, query_name, query, input_sequence_type)
}

/// C++ `finish_daa_query_record(TextBuffer&, size_t)`.
pub fn finish_daa_query_record(buf: &mut [u8], seek_pos: usize) {
    daa::finish_daa_query_record(buf, seek_pos)
}

/// C++ `write_daa_record(TextBuffer&, const IntermediateRecord&)`.
///
/// Rust names the overload by its argument type.
pub fn write_daa_record_intermediate(buf: &mut Vec<u8>, record: &IntermediateRecord) {
    daa::write_daa_record_intermediate(buf, record)
}

/// C++ `write_daa_record(TextBuffer&, const Hsp&, uint32_t)`.
///
/// Rust names the overload by its argument type.
pub fn write_daa_record_hsp(buf: &mut Vec<u8>, hsp: &Hsp, subject_id: u32) {
    daa::write_daa_record_hsp(buf, hsp, subject_id)
}

/// C++ `finish_daa(OutputFile&, SequenceFile&)`.
pub fn finish_daa_from_sequence_file<W, D>(
    writer: &mut W,
    db: &D,
    metadata: &DaaRunMetadata,
) -> io::Result<()>
where
    W: Write + Seek,
    D: DaaSequenceFile,
{
    let mut header = DaaHeader2::with_params(
        db.sequence_count(),
        metadata.db_letters,
        metadata.gap_open,
        metadata.gap_extend,
        metadata.reward,
        metadata.penalty,
        metadata.k,
        metadata.lambda,
        metadata.max_evalue,
        &metadata.matrix.to_ascii_lowercase(),
        metadata.mode,
    );
    header.block_type[0] = daa::BlockType::Alignments as u8;
    header.block_type[1] = daa::BlockType::RefNames as u8;
    header.block_type[2] = daa::BlockType::RefLengths as u8;

    terminate_aln_block(writer, &mut header)?;
    let count = db.dict_size();
    header.db_seqs_used = count as u64;
    header.query_records = metadata.aligned_queries;

    let mut names_size = 0u64;
    for index in 0..count {
        let title = db.dict_title(index);
        writer.write_all(title.as_bytes())?;
        writer.write_all(&[0])?;
        names_size += title.len() as u64 + 1;
    }
    header.block_size[1] = names_size;

    for index in 0..count {
        writer.write_all(&db.dict_len(index).to_ne_bytes())?;
    }
    header.block_size[2] = (count * std::mem::size_of::<u32>()) as u64;

    write_header2(writer, &header)
}

/// C++ `finish_daa(OutputFile&, DAAFile&)`.
pub fn finish_daa_from_file<W: Write + Seek>(writer: &mut W, daa_in: &DaaFile) -> io::Result<()> {
    daa::finish_daa_from_file(writer, daa_in)
}

/// C++ `finish_daa(OutputFile&, DAAFile&, StringSet, vector, int64_t)`.
pub fn finish_daa_from_refs<W: Write + Seek>(
    writer: &mut W,
    daa_in: &DaaFile,
    seq_ids: &[String],
    seq_lens: &[u32],
    query_count: i64,
) -> io::Result<()> {
    daa::finish_daa_from_refs(writer, daa_in, seq_ids, seq_lens, query_count)
}

/// C++ file-local `terminate_aln_block`.
fn terminate_aln_block<W: Write + Seek>(writer: &mut W, header: &mut DaaHeader2) -> io::Result<()> {
    daa::terminate_aln_block(writer, header)
}

/// C++ file-local `write_header2`.
fn write_header2<W: Write + Seek>(writer: &mut W, header: &DaaHeader2) -> io::Result<()> {
    daa::write_header2(writer, header)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::daa::{
        DaaHeader1, DAA_HEADER1_SIZE, DAA_HEADER2_SIZE, DAA_MAGIC_NUMBER, DAA_VERSION,
    };
    use std::io::{Cursor, SeekFrom};

    #[derive(Debug)]
    struct TestDb {
        sequence_count: u64,
        refs: Vec<(String, u32)>,
    }

    impl DaaSequenceFile for TestDb {
        fn sequence_count(&self) -> u64 {
            self.sequence_count
        }

        fn dict_size(&self) -> usize {
            self.refs.len()
        }

        fn dict_title(&self, index: usize) -> String {
            self.refs[index].0.clone()
        }

        fn dict_len(&self, index: usize) -> u32 {
            self.refs[index].1
        }
    }

    #[test]
    fn init_daa_writes_the_two_wire_headers() {
        let mut bytes = Vec::new();
        init_daa(&mut bytes).unwrap();

        assert_eq!(bytes.len(), DAA_HEADER1_SIZE + DAA_HEADER2_SIZE);
        let mut input = bytes.as_slice();
        let h1 = DaaHeader1::read_from(&mut input).unwrap();
        assert_eq!(h1.magic_number, DAA_MAGIC_NUMBER);
        assert_eq!(h1.version, DAA_VERSION);
        let h2 = DaaHeader2::read_from(&mut input).unwrap();
        let default = DaaHeader2::default();
        assert_eq!(h2.diamond_build, default.diamond_build);
        assert_eq!(h2.db_seqs, 0);
        assert_eq!(h2.block_size, [0; 256]);
        assert_eq!(h2.block_type, [0; 256]);
    }

    #[test]
    fn query_record_is_byte_exact_and_patches_its_size() {
        let mut bytes = vec![0xaa];
        let start = write_daa_query_record(
            &mut bytes,
            "query one comment",
            &[0, 1, 2],
            SequenceType::AminoAcid,
        );
        assert_eq!(start, 1);
        finish_daa_query_record(&mut bytes, start);

        // size=13, query length=3, C string `query`, flags=0, packed 0,1,2.
        assert_eq!(
            bytes,
            vec![0xaa, 13, 0, 0, 0, 3, 0, 0, 0, b'q', b'u', b'e', b'r', b'y', 0, 0, 0x20, 0x08]
        );
    }

    #[test]
    fn sequence_file_finalizer_writes_exact_blocks_and_metadata() {
        let mut output = Cursor::new(Vec::new());
        init_daa(&mut output).unwrap();
        output.write_all(&[0x11, 0x22]).unwrap();
        let db = TestDb {
            sequence_count: 9,
            refs: vec![("alpha".into(), 101), ("beta desc".into(), 202)],
        };
        let metadata = DaaRunMetadata {
            db_letters: 12_345,
            gap_open: 11,
            gap_extend: 1,
            reward: 2,
            penalty: -3,
            k: 0.041,
            lambda: 0.267,
            max_evalue: 0.001,
            matrix: "BLOSUM62".into(),
            mode: 2,
            aligned_queries: 7,
        };

        finish_daa_from_sequence_file(&mut output, &db, &metadata).unwrap();
        let bytes = output.into_inner();
        let block_start = DAA_HEADER1_SIZE + DAA_HEADER2_SIZE;
        assert_eq!(
            &bytes[block_start..block_start + 6],
            &[0x11, 0x22, 0, 0, 0, 0]
        );
        assert_eq!(
            &bytes[block_start + 6..block_start + 22],
            b"alpha\0beta desc\0"
        );
        assert_eq!(
            &bytes[block_start + 22..block_start + 26],
            &101u32.to_ne_bytes()
        );
        assert_eq!(
            &bytes[block_start + 26..block_start + 30],
            &202u32.to_ne_bytes()
        );

        let mut input = Cursor::new(bytes);
        input
            .seek(SeekFrom::Start(DAA_HEADER1_SIZE as u64))
            .unwrap();
        let header = DaaHeader2::read_from(&mut input).unwrap();
        assert_eq!(header.db_seqs, 9);
        assert_eq!(header.db_seqs_used, 2);
        assert_eq!(header.db_letters, 12_345);
        assert_eq!(header.query_records, 7);
        assert_eq!(header.matrix_name(), "blosum62");
        assert_eq!(header.block_size[0], 6);
        assert_eq!(header.block_size[1], 16);
        assert_eq!(header.block_size[2], 8);
        assert_eq!(header.block_type[..3], [1, 2, 3]);

        assert_eq!(input.into_inner().len(), block_start + 30);
    }
}
