use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;

use crate::align::hsp::{Hsp, HspContext};
use crate::basic::packed_sequence::PackedSequence;
use crate::basic::packed_transcript::{PackedOperation, PackedTranscript};
use crate::basic::translate::translate_6_frames;
use crate::basic::value::{Letter, SequenceType};
use crate::output::format::{
    write_tabular_context_row, write_tabular_context_row_json, write_tabular_query_intro,
    TabularFormat,
};
use crate::output::get_segment_flag_hsp;
use crate::output::intermediate::IntermediateRecord;
use crate::output::{edge, paf, pairwise, sam, xml};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::binary_buffer::BinaryBufferIterator;
use crate::util::interval::Interval;
use crate::util::text_buffer::find_first_of;

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
        if self.version > DAA_VERSION {
            return Err("DAA version requires later version of DIAMOND.".into());
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
            // C++ `DAA_header2()` (`daa_file.h:53-58`) initializes
            // `diamond_build(Const::build_version)`. Rust used to write 0 at
            // file offset 32, byte-diverging from any C++-produced DAA and
            // misleading downstream `view` clients keying on build version.
            diamond_build: crate::data::dmnd::BUILD_VERSION as u64,
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

    /// Matches C++ `DAA_header2::DAA_header2(db_seqs, db_letters, gap_open, gap_extend, reward, penalty, k, lambda, evalue, score_matrix, mode)`.
    #[allow(clippy::too_many_arguments)]
    pub fn with_params(
        db_seqs: u64,
        db_letters: u64,
        gap_open: i32,
        gap_extend: i32,
        reward: i32,
        penalty: i32,
        k: f64,
        lambda: f64,
        evalue: f64,
        score_matrix: &str,
        mode: u32,
    ) -> Self {
        let mut h = Self {
            db_seqs,
            db_letters,
            gap_open,
            gap_extend,
            reward,
            penalty,
            k,
            lambda,
            evalue,
            mode: mode as i32,
            ..Self::default()
        };
        h.set_matrix_name(score_matrix);
        h
    }

    /// Matches C++ `DAA_header2::DAA_header2(daa)`.
    pub fn from_daa_file(daa: &DaaFile) -> Self {
        let mut h = Self::with_params(
            daa.db_seqs(),
            daa.db_letters(),
            daa.gap_open_penalty(),
            daa.gap_extension_penalty(),
            daa.match_reward(),
            daa.mismatch_penalty(),
            daa.kappa(),
            daa.lambda(),
            daa.evalue(),
            &daa.score_matrix(),
            daa.mode(),
        );
        h.block_type[0] = BlockType::Alignments as u8;
        h.block_type[1] = BlockType::RefNames as u8;
        h.block_type[2] = BlockType::RefLengths as u8;
        h
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
        let end = self.score_matrix.iter().position(|&b| b == 0).unwrap_or(16);
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

/// Reader state for a DAA file.
pub struct DaaFile {
    reader: BufReader<File>,
    query_count: usize,
    h1: DaaHeader1,
    h2: DaaHeader2,
    ref_name: Vec<String>,
    ref_len: Vec<u32>,
}

impl DaaFile {
    /// Matches C++ `DAAFile::DAAFile(file_name)`.
    pub fn open<P: AsRef<Path>>(file_name: P) -> io::Result<Self> {
        let mut reader = BufReader::new(File::open(file_name)?);
        let h1 = DaaHeader1::read_from(&mut reader)?;
        h1.validate()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let h2 = DaaHeader2::read_from(&mut reader)?;

        if h2.block_size[0] == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid DAA file. DIAMOND run has probably not completed successfully.",
            ));
        }

        reader.seek(SeekFrom::Start(
            (DAA_HEADER1_SIZE + DAA_HEADER2_SIZE) as u64 + h2.block_size[0],
        ))?;
        let mut ref_name = Vec::with_capacity(h2.db_seqs_used as usize);
        for _ in 0..h2.db_seqs_used {
            let mut s = Vec::new();
            loop {
                let mut c = [0u8; 1];
                reader.read_exact(&mut c)?;
                if c[0] == 0 {
                    break;
                }
                s.push(c[0]);
            }
            ref_name.push(String::from_utf8_lossy(&s).to_string());
        }

        let mut ref_len = Vec::with_capacity(h2.db_seqs_used as usize);
        for _ in 0..h2.db_seqs_used {
            let mut buf = [0u8; 4];
            reader.read_exact(&mut buf)?;
            ref_len.push(u32::from_le_bytes(buf));
        }

        reader.seek(SeekFrom::Start(
            (DAA_HEADER1_SIZE + DAA_HEADER2_SIZE) as u64,
        ))?;

        Ok(Self {
            reader,
            query_count: 0,
            h1,
            h2,
            ref_name,
            ref_len,
        })
    }

    /// Matches C++ `DAAFile::diamond_build()`.
    pub fn diamond_build(&self) -> u64 {
        self.h2.diamond_build
    }

    /// Matches C++ `DAAFile::db_seqs()`.
    pub fn db_seqs(&self) -> u64 {
        self.h2.db_seqs
    }

    /// Matches C++ `DAAFile::db_seqs_used()`.
    pub fn db_seqs_used(&self) -> u64 {
        self.h2.db_seqs_used
    }

    /// Matches C++ `DAAFile::db_letters()`.
    pub fn db_letters(&self) -> u64 {
        self.h2.db_letters
    }

    /// Matches C++ `DAAFile::score_matrix()`.
    pub fn score_matrix(&self) -> String {
        self.h2.matrix_name()
    }

    /// Matches C++ `DAAFile::gap_open_penalty()`.
    pub fn gap_open_penalty(&self) -> i32 {
        self.h2.gap_open
    }

    /// Matches C++ `DAAFile::gap_extension_penalty()`.
    pub fn gap_extension_penalty(&self) -> i32 {
        self.h2.gap_extend
    }

    /// Matches C++ `DAAFile::match_reward()`.
    pub fn match_reward(&self) -> i32 {
        self.h2.reward
    }

    /// Matches C++ `DAAFile::mismatch_penalty()`.
    pub fn mismatch_penalty(&self) -> i32 {
        self.h2.penalty
    }

    /// Matches C++ `DAAFile::query_records()`.
    pub fn query_records(&self) -> u64 {
        self.h2.query_records
    }

    /// Matches C++ `DAAFile::mode()`.
    pub fn mode(&self) -> u32 {
        self.h2.mode as u32
    }

    /// Matches C++ `DAAFile::ref_name(i)`.
    pub fn ref_name(&self, i: usize) -> &str {
        &self.ref_name[i]
    }

    /// Matches C++ `DAAFile::ref_len(i)`.
    pub fn ref_len_at(&self, i: usize) -> u32 {
        self.ref_len[i]
    }

    /// Matches C++ `DAAFile::lambda()`.
    pub fn lambda(&self) -> f64 {
        self.h2.lambda
    }

    /// Matches C++ `DAAFile::kappa()`.
    pub fn kappa(&self) -> f64 {
        self.h2.k
    }

    /// Matches C++ `DAAFile::evalue()`.
    pub fn evalue(&self) -> f64 {
        self.h2.evalue
    }

    /// Matches C++ `DAAFile::block_size(i)`.
    pub fn block_size(&self, i: usize) -> u64 {
        self.h2.block_size[i]
    }

    /// Matches C++ `DAAFile::ref_len()`.
    pub fn ref_len(&self) -> &[u32] {
        &self.ref_len
    }

    /// Matches C++ `DAAFile::read_query_buffer()`.
    pub fn read_query_buffer(&mut self) -> io::Result<Option<(Vec<u8>, usize)>> {
        let mut size_buf = [0u8; 4];
        self.reader.read_exact(&mut size_buf)?;
        let size = u32::from_le_bytes(size_buf);
        if size == 0 {
            return Ok(None);
        }
        let mut buf = vec![0u8; size as usize];
        self.reader.read_exact(&mut buf)?;
        let query_num = self.query_count;
        self.query_count += 1;
        Ok(Some((buf, query_num)))
    }

    pub fn header1(&self) -> &DaaHeader1 {
        &self.h1
    }

    pub fn header2(&self) -> &DaaHeader2 {
        &self.h2
    }
}

pub const DAA_HEADER1_SIZE: usize = 16;
pub const DAA_HEADER2_SIZE: usize = 2432;
pub const DAA_ID_DELIMITERS: &str = " \x07\x08\x0c\n\r\t\x0b\x01";

/// Query record decoded from a DAA alignment block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DaaQueryRecord {
    pub query_name: String,
    pub query_num: usize,
    pub source_seq: Vec<Letter>,
    pub context: [Vec<Letter>; 6],
    raw_matches: Vec<u8>,
    mode: u32,
}

impl DaaQueryRecord {
    /// Matches C++ `DAA_query_record::DAA_query_record(file, buf, query_num)`.
    pub fn from_buffer(file: &DaaFile, buf: &[u8], query_num: usize) -> Result<Self, String> {
        let mut it = BinaryBufferIterator::new(buf);
        let query_len = it.read::<u32>()? as usize;
        let query_name = it.read_string()?;
        let flags = it.read::<u8>()?;
        let mut source_seq = Vec::new();
        let mut context: [Vec<Letter>; 6] = Default::default();

        if file.mode() == 2 {
            let byte_count = (query_len * 5).div_ceil(8);
            let packed = PackedSequence::from_raw(it.read_bytes(byte_count)?.to_vec(), false);
            context[0] = packed.unpack(5, query_len);
        } else {
            let have_n = (flags & 1) == 1;
            let bits = if have_n { 3 } else { 2 };
            let byte_count = (query_len * bits as usize).div_ceil(8);
            let packed = PackedSequence::from_raw(it.read_bytes(byte_count)?.to_vec(), have_n);
            source_seq = packed.unpack(bits, query_len);
            let dna: Vec<u8> = source_seq.iter().map(|&x| x as u8).collect();
            context = translate_6_frames(&dna).map(|frame| {
                frame
                    .into_iter()
                    .map(|letter| letter as Letter)
                    .collect::<Vec<_>>()
            });
        }

        Ok(Self {
            query_name,
            query_num,
            source_seq,
            context,
            raw_matches: it.remaining().to_vec(),
            mode: file.mode(),
        })
    }

    /// Matches C++ `DAA_query_record::begin()`.
    pub fn raw_matches(&self) -> &[u8] {
        &self.raw_matches
    }

    /// Matches C++ `DAA_query_record::raw_begin()`.
    pub fn raw_begin(&self) -> BinaryBufferIterator<'_> {
        BinaryBufferIterator::new(&self.raw_matches)
    }

    /// Matches C++ `DAA_query_record::begin(file, score_matrix)`.
    pub fn begin<'a>(
        &'a self,
        file: &'a DaaFile,
        score_matrix: &'a ScoreMatrix,
    ) -> DaaMatchIterator<'a> {
        DaaMatchIterator::new(self, file, score_matrix)
    }

    /// Matches C++ `DAA_query_record::query_len()`.
    pub fn query_len(&self) -> usize {
        if self.mode == 3 {
            self.source_seq.len()
        } else {
            self.context[0].len()
        }
    }

    pub fn query_source(&self) -> &[Letter] {
        if self.mode == 3 {
            &self.source_seq
        } else {
            &self.context[0]
        }
    }

    pub fn input_sequence_type(&self) -> SequenceType {
        if self.mode == 2 {
            SequenceType::AminoAcid
        } else {
            SequenceType::Nucleotide
        }
    }
}

/// Decoded DAA match.
#[derive(Debug, Clone)]
pub struct DaaMatch {
    pub hsp: Hsp,
    pub hsp_num: u32,
    pub hit_num: u32,
    pub subject_id: u32,
    pub subject_len: u32,
    pub subject_name: String,
}

impl DaaMatch {
    /// Matches C++ `DAA_query_record::Match::context(parent)`.
    pub fn context(&self, parent: &DaaQueryRecord) -> HspContext {
        HspContext::new(
            self.hsp.clone(),
            parent.query_num as u32,
            parent.query_num as u64,
            parent.context.to_vec(),
            parent.query_len() as i32,
            &parent.query_name,
            self.subject_id as u64,
            self.subject_len as i32,
            &self.subject_name,
            self.hit_num,
            self.hsp_num,
            Vec::new(),
            0.0,
            0.0,
        )
    }
}

/// Iterator over decoded DAA matches.
pub struct DaaMatchIterator<'a> {
    parent: &'a DaaQueryRecord,
    file: &'a DaaFile,
    score_matrix: &'a ScoreMatrix,
    it: BinaryBufferIterator<'a>,
    subject_id: u32,
    hsp_num: u32,
    hit_num: u32,
}

impl<'a> DaaMatchIterator<'a> {
    pub fn new(
        parent: &'a DaaQueryRecord,
        file: &'a DaaFile,
        score_matrix: &'a ScoreMatrix,
    ) -> Self {
        Self {
            parent,
            file,
            score_matrix,
            it: parent.raw_begin(),
            subject_id: u32::MAX,
            hsp_num: u32::MAX,
            hit_num: u32::MAX,
        }
    }
}

impl<'a> Iterator for DaaMatchIterator<'a> {
    type Item = Result<DaaMatch, String>;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.it.good() {
            return None;
        }

        let raw = match DaaRawMatch::read(&mut self.it) {
            Ok(raw) => raw,
            Err(e) => return Some(Err(e)),
        };

        let old_subject = self.subject_id;
        self.subject_id = raw.subject_id;
        if self.subject_id == old_subject {
            self.hsp_num = self.hsp_num.wrapping_add(1);
        } else {
            self.hsp_num = 0;
            self.hit_num = self.hit_num.wrapping_add(1);
        }

        if self.subject_id as usize >= self.file.ref_name.len() {
            return Some(Err("Subject id out of bounds.".to_string()));
        }

        let mut hsp = Hsp::with_backtraced(true);
        hsp.score = raw.score as i32;
        hsp.subject_range = Interval::new(raw.subject_begin as i32, raw.subject_begin as i32);
        hsp.transcript = PackedTranscript::from_bytes(&raw.transcript);

        if self.file.mode() == 3 {
            hsp.frame = if (raw.flag & (1 << 6)) == 0 {
                (raw.query_begin % 3) as i32
            } else {
                3 + ((self.parent.source_seq.len() as u32 - 1 - raw.query_begin) % 3) as i32
            };
            hsp.set_translated_query_begin(
                raw.query_begin as i32,
                self.parent.source_seq.len() as i32,
            );
        } else if self.file.mode() == 2 {
            hsp.frame = 0;
            hsp.query_range.begin = raw.query_begin as i32;
        }

        let subject_len = self.file.ref_len_at(self.subject_id as usize);
        let subject_name = self.file.ref_name(self.subject_id as usize).to_string();
        let mut parsed = DaaMatch {
            hsp,
            hsp_num: self.hsp_num,
            hit_num: self.hit_num,
            subject_id: self.subject_id,
            subject_len,
            subject_name,
        };
        let mut context = parsed.context(self.parent);
        if let Err(e) = context.parse(true, true, self.file.mode() == 3, self.score_matrix) {
            return Some(Err(e));
        }
        parsed.hsp = context.hsp();
        parsed.hsp.evalue = self.score_matrix.evalue(
            parsed.hsp.score,
            self.parent.context[0].len() as u32,
            parsed.subject_len,
        );
        parsed.hsp.bit_score = self.score_matrix.bitscore(parsed.hsp.score as f64);
        Some(Ok(parsed))
    }
}

/// Raw fields read by C++ `DAA_query_record::Match::read` before HSP parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DaaRawMatch {
    pub subject_id: u32,
    pub flag: u8,
    pub score: u32,
    pub query_begin: u32,
    pub subject_begin: u32,
    pub transcript: Vec<u8>,
}

impl DaaRawMatch {
    pub fn read(it: &mut BinaryBufferIterator<'_>) -> Result<Self, String> {
        let subject_id = it.read::<u32>()?;
        let flag = it.read::<u8>()?;
        let score = it.read_packed_u32(flag & 3)?;
        let query_begin = it.read_packed_u32((flag >> 2) & 3)?;
        let subject_begin = it.read_packed_u32((flag >> 4) & 3)?;
        let mut transcript = Vec::new();
        loop {
            let c = it.read::<u8>()?;
            transcript.push(c);
            if c == PackedOperation::terminator().code {
                break;
            }
        }
        Ok(Self {
            subject_id,
            flag,
            score,
            query_begin,
            subject_begin,
            transcript,
        })
    }
}

/// Matches C++ `copy_match_record_raw(it, out, subject_map)`.
pub fn copy_match_record_raw(
    it: &mut BinaryBufferIterator<'_>,
    out: &mut Vec<u8>,
    subject_map: &HashMap<u32, u32>,
) -> Result<(), String> {
    let r = DaaRawMatch::read(it)?;
    let mapped = subject_map
        .get(&r.subject_id)
        .ok_or_else(|| "subject id not found".to_string())?;
    out.extend_from_slice(&mapped.to_ne_bytes());
    out.push(r.flag);
    write_packed_width(out, r.score, r.flag & 3);
    write_packed_width(out, r.query_begin, (r.flag >> 2) & 3);
    write_packed_width(out, r.subject_begin, (r.flag >> 4) & 3);
    out.extend_from_slice(&r.transcript);
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DaaViewConfig {
    pub max_target_seqs: i64,
    pub toppercent: Option<f64>,
    pub forward_only: bool,
}

impl Default for DaaViewConfig {
    fn default() -> Self {
        Self {
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
        }
    }
}

/// Matches C++ `Config::output_range(n_target_seq, score, top_score, max_target_seqs, toppercent)`.
pub fn output_range(
    n_target_seq: u32,
    score: i32,
    top_score: i32,
    max_target_seqs: i64,
    toppercent: Option<f64>,
) -> bool {
    if let Some(toppercent) = toppercent {
        top_score != 0 && (1.0 - score as f64 / top_score as f64) * 100.0 <= toppercent
    } else {
        (n_target_seq as i64) < max_target_seqs
    }
}

/// Matches C++ `view_query(r, file, out, score_matrix, cfg)`.
pub fn view_query_daa(
    r: &DaaQueryRecord,
    file: &DaaFile,
    out: &mut Vec<u8>,
    score_matrix: &ScoreMatrix,
    cfg: &DaaViewConfig,
) -> Result<(), String> {
    let seek_pos = write_daa_query_record(
        out,
        &r.query_name,
        r.query_source(),
        r.input_sequence_type(),
    );
    let matches = r.begin(file, score_matrix).collect::<Result<Vec<_>, _>>()?;
    let top_score = matches.first().map(|m| m.hsp.score).unwrap_or(0);

    for m in matches {
        if m.hsp.frame > 2 && cfg.forward_only {
            continue;
        }
        if !output_range(
            m.hit_num,
            m.hsp.score,
            top_score,
            cfg.max_target_seqs,
            cfg.toppercent,
        ) {
            break;
        }
        write_daa_record_hsp(out, &m.hsp, m.subject_id);
    }
    finish_daa_query_record(out, seek_pos);
    Ok(())
}

/// Matches C++ `view_query(r, file, out, score_matrix, format, cfg, report_unaligned)`.
pub fn view_query_tabular(
    r: &DaaQueryRecord,
    file: &DaaFile,
    out: &mut Vec<u8>,
    score_matrix: &ScoreMatrix,
    format: &TabularFormat,
    cfg: &DaaViewConfig,
    report_unaligned: bool,
) -> Result<(), String> {
    if !format.is_json {
        write_tabular_query_intro(
            out,
            &r.query_name,
            r.query_len() as i32,
            r.query_source(),
            &format.fields,
            report_unaligned,
        )
        .map_err(|e| e.to_string())?;
    }

    let matches = r.begin(file, score_matrix).collect::<Result<Vec<_>, _>>()?;
    let top_score = matches.first().map(|m| m.hsp.score).unwrap_or(0);

    for m in matches {
        if m.hsp.frame > 2 && cfg.forward_only {
            continue;
        }
        if !output_range(
            m.hit_num,
            m.hsp.score,
            top_score,
            cfg.max_target_seqs,
            cfg.toppercent,
        ) {
            break;
        }
        let context = m.context(r);
        if format.is_json {
            write_tabular_context_row_json(
                out,
                &context,
                &format.fields,
                file.mode() == 3,
                false,
                0,
                0,
                &[],
                None,
                "",
                None,
            )
            .map_err(|e| e.to_string())?;
        } else {
            write_tabular_context_row(out, &context, &format.fields, file.mode() == 3, false, 0, 0)
                .map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

/// Matches C++ `view_query(r, file, out, score_matrix, cfg)`.
pub fn view_query_paf(
    r: &DaaQueryRecord,
    file: &DaaFile,
    out: &mut Vec<u8>,
    score_matrix: &ScoreMatrix,
    cfg: &DaaViewConfig,
) -> Result<(), String> {
    paf::print_query_intro(out, &r.query_name, false).map_err(|e| e.to_string())?;
    let matches = r.begin(file, score_matrix).collect::<Result<Vec<_>, _>>()?;
    let top_score = matches.first().map(|m| m.hsp.score).unwrap_or(0);

    for m in matches {
        if m.hsp.frame > 2 && cfg.forward_only {
            continue;
        }
        if !output_range(
            m.hit_num,
            m.hsp.score,
            top_score,
            cfg.max_target_seqs,
            cfg.toppercent,
        ) {
            break;
        }
        let context = m.context(r);
        paf::print_match_context(out, &context, score_matrix).map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// Matches C++ `view_query(r, file, out, score_matrix, cfg, sam_qlen_field)`.
pub fn view_query_sam(
    r: &DaaQueryRecord,
    file: &DaaFile,
    out: &mut Vec<u8>,
    score_matrix: &ScoreMatrix,
    cfg: &DaaViewConfig,
    sam_qlen_field: bool,
) -> Result<(), String> {
    sam::print_query_intro(out, &r.query_name, false).map_err(|e| e.to_string())?;
    let matches = r.begin(file, score_matrix).collect::<Result<Vec<_>, _>>()?;
    let top_score = matches.first().map(|m| m.hsp.score).unwrap_or(0);

    for m in matches {
        if m.hsp.frame > 2 && cfg.forward_only {
            continue;
        }
        if !output_range(
            m.hit_num,
            m.hsp.score,
            top_score,
            cfg.max_target_seqs,
            cfg.toppercent,
        ) {
            break;
        }
        let context = m.context(r);
        sam::print_match_context(out, &context, score_matrix, true, true, sam_qlen_field)
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// Matches C++ `view_query(r, file, out, score_matrix, cfg)`.
pub fn view_query_pairwise(
    r: &DaaQueryRecord,
    file: &DaaFile,
    out: &mut Vec<u8>,
    score_matrix: &ScoreMatrix,
    cfg: &DaaViewConfig,
) -> Result<(), String> {
    pairwise::print_query_intro(out, &r.query_name, r.query_len() as i32, false)
        .map_err(|e| e.to_string())?;
    let matches = r.begin(file, score_matrix).collect::<Result<Vec<_>, _>>()?;
    let top_score = matches.first().map(|m| m.hsp.score).unwrap_or(0);

    for m in matches {
        if m.hsp.frame > 2 && cfg.forward_only {
            continue;
        }
        if !output_range(
            m.hit_num,
            m.hsp.score,
            top_score,
            cfg.max_target_seqs,
            cfg.toppercent,
        ) {
            break;
        }
        let context = m.context(r);
        pairwise::print_match_context(out, &context, score_matrix, file.mode() == 3, false)
            .map_err(|e| e.to_string())?;
    }
    pairwise::print_query_epilog(out).map_err(|e| e.to_string())?;
    Ok(())
}

/// Matches C++ `view_query(r, file, out, score_matrix, cfg)`.
pub fn view_query_xml(
    r: &DaaQueryRecord,
    file: &DaaFile,
    out: &mut Vec<u8>,
    score_matrix: &ScoreMatrix,
    cfg: &DaaViewConfig,
) -> Result<(), String> {
    xml::print_query_intro(out, r.query_num as u64, &r.query_name, r.query_len() as i32)
        .map_err(|e| e.to_string())?;
    let matches = r.begin(file, score_matrix).collect::<Result<Vec<_>, _>>()?;
    let top_score = matches.first().map(|m| m.hsp.score).unwrap_or(0);

    for m in matches {
        if m.hsp.frame > 2 && cfg.forward_only {
            continue;
        }
        if !output_range(
            m.hit_num,
            m.hsp.score,
            top_score,
            cfg.max_target_seqs,
            cfg.toppercent,
        ) {
            break;
        }
        let context = m.context(r);
        xml::print_match_context(out, &context, score_matrix, file.mode() == 3, false, false)
            .map_err(|e| e.to_string())?;
    }
    xml::print_query_epilog(
        out,
        false,
        file.db_seqs(),
        file.db_letters(),
        score_matrix.k(),
        score_matrix.lambda(),
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

/// Matches C++ `view_query(r, file, score_matrix, cfg)`.
pub fn view_query_null(
    r: &DaaQueryRecord,
    file: &DaaFile,
    score_matrix: &ScoreMatrix,
    cfg: &DaaViewConfig,
) -> Result<(), String> {
    let matches = r.begin(file, score_matrix).collect::<Result<Vec<_>, _>>()?;
    let top_score = matches.first().map(|m| m.hsp.score).unwrap_or(0);

    for m in matches {
        if m.hsp.frame > 2 && cfg.forward_only {
            continue;
        }
        if !output_range(
            m.hit_num,
            m.hsp.score,
            top_score,
            cfg.max_target_seqs,
            cfg.toppercent,
        ) {
            break;
        }
    }
    Ok(())
}

/// Matches C++ `view_query(r, file, out, score_matrix, cfg)`.
pub fn view_query_edge(
    r: &DaaQueryRecord,
    file: &DaaFile,
    out: &mut Vec<u8>,
    score_matrix: &ScoreMatrix,
    cfg: &DaaViewConfig,
) -> Result<(), String> {
    let matches = r.begin(file, score_matrix).collect::<Result<Vec<_>, _>>()?;
    let top_score = matches.first().map(|m| m.hsp.score).unwrap_or(0);

    for m in matches {
        if m.hsp.frame > 2 && cfg.forward_only {
            continue;
        }
        if !output_range(
            m.hit_num,
            m.hsp.score,
            top_score,
            cfg.max_target_seqs,
            cfg.toppercent,
        ) {
            break;
        }
        let context = m.context(r);
        edge::print_match(out, &context).map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// Matches C++ `build_mapping(acc2oid, seq_ids, seq_lens, f)`.
pub fn build_mapping(
    acc2oid: &mut HashMap<String, u32>,
    seq_ids: &mut Vec<String>,
    seq_lens: &mut Vec<u32>,
    f: &DaaFile,
) -> HashMap<u32, u32> {
    let mut r = HashMap::new();
    for i in 0..f.db_seqs_used() as usize {
        let name = f.ref_name(i).to_string();
        let next = acc2oid.len() as u32;
        let oid = *acc2oid.entry(name.clone()).or_insert(next);
        r.insert(i as u32, oid);
        if oid == next {
            seq_ids.push(name);
            seq_lens.push(f.ref_len_at(i));
        }
    }
    r
}

/// Matches C++ `write_file(f, out, subject_map)`.
pub fn write_file<W: Write>(
    f: &mut DaaFile,
    out: &mut W,
    subject_map: &HashMap<u32, u32>,
) -> io::Result<i64> {
    let mut out_buf = Vec::new();
    let mut last_query_num = None;
    while let Some((buf, query_num)) = f.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(f, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let seek_pos = write_daa_query_record(
            &mut out_buf,
            &r.query_name,
            r.query_source(),
            r.input_sequence_type(),
        );
        let mut it = r.raw_begin();
        while it.good() {
            copy_match_record_raw(&mut it, &mut out_buf, subject_map)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        }
        finish_daa_query_record(&mut out_buf, seek_pos);
        out.write_all(&out_buf)?;
        out_buf.clear();
        last_query_num = Some(query_num);
    }
    Ok(last_query_num.map_or(0, |query_num| query_num as i64 + 1))
}

/// Matches C++ `merge_daa(input_files, output_file)`.
pub fn merge_daa_files<P: AsRef<Path>>(input_files: &[P], output_file: P) -> io::Result<i64> {
    if input_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Missing parameter: input files (--in)",
        ));
    }

    let mut files = Vec::with_capacity(input_files.len());
    let mut acc2oid = HashMap::new();
    let mut oid_maps = Vec::with_capacity(input_files.len());
    let mut seq_ids = Vec::new();
    let mut seq_lens = Vec::new();

    for file in input_files {
        let daa = DaaFile::open(file)?;
        oid_maps.push(build_mapping(
            &mut acc2oid,
            &mut seq_ids,
            &mut seq_lens,
            &daa,
        ));
        files.push(daa);
    }

    let mut out = File::create(output_file)?;
    init_daa(&mut out)?;
    let mut query_count = 0;
    for (file, subject_map) in files.iter_mut().zip(oid_maps.iter()) {
        query_count += write_file(file, &mut out, subject_map)?;
    }
    finish_daa_from_refs(&mut out, &files[0], &seq_ids, &seq_lens, query_count)?;
    Ok(query_count)
}

/// Matches C++ `write_daa_query_record(buf, query_name, query, input_sequence_type)`.
pub fn write_daa_query_record(
    buf: &mut Vec<u8>,
    query_name: &str,
    query: &[Letter],
    input_sequence_type: SequenceType,
) -> usize {
    let seek_pos = buf.len();
    buf.extend_from_slice(&0u32.to_ne_bytes());
    buf.extend_from_slice(&(query.len() as u32).to_ne_bytes());
    let name_len = find_first_of(query_name, DAA_ID_DELIMITERS);
    buf.extend_from_slice(&query_name.as_bytes()[..name_len]);
    buf.push(0);
    let s = PackedSequence::new(query, input_sequence_type);
    buf.push(if s.has_n() { 1 } else { 0 });
    buf.extend_from_slice(s.data());
    seek_pos
}

/// Matches C++ `finish_daa_query_record(buf, seek_pos)`.
pub fn finish_daa_query_record(buf: &mut [u8], seek_pos: usize) {
    let size = (buf.len() - seek_pos - std::mem::size_of::<u32>()) as u32;
    buf[seek_pos..seek_pos + std::mem::size_of::<u32>()].copy_from_slice(&size.to_ne_bytes());
}

/// Matches C++ `write_daa_record(buf, r)`.
pub fn write_daa_record_intermediate(buf: &mut Vec<u8>, r: &IntermediateRecord) {
    buf.extend_from_slice(&(r.target_dict_id as u32).to_ne_bytes());
    buf.push(r.flag);
    write_packed(buf, r.score);
    write_packed(buf, r.query_begin);
    write_packed(buf, r.subject_begin);
    for op in r.transcript.data() {
        buf.push(op.code);
    }
}

/// Matches C++ `write_daa_record(buf, hsp, subject_id)`.
pub fn write_daa_record_hsp(buf: &mut Vec<u8>, hsp: &Hsp, subject_id: u32) {
    let oriented_range = hsp.oriented_range();
    buf.extend_from_slice(&subject_id.to_ne_bytes());
    buf.push(get_segment_flag_hsp(hsp));
    write_packed(buf, hsp.score as u32);
    write_packed(buf, oriented_range.begin as u32);
    write_packed(buf, hsp.subject_range.begin as u32);
    for op in hsp.transcript.data() {
        buf.push(op.code);
    }
}

/// Matches C++ `init_daa(writer)`.
pub fn init_daa<W: Write>(writer: &mut W) -> io::Result<()> {
    DaaHeader1::new().write_to(writer)?;
    DaaHeader2::new().write_to(writer)
}

/// Matches C++ `terminate_aln_block(writer, h2)`.
pub fn terminate_aln_block<W: Write + Seek>(writer: &mut W, h2: &mut DaaHeader2) -> io::Result<()> {
    writer.write_all(&0u32.to_ne_bytes())?;
    h2.block_size[0] = writer.stream_position()? - (DAA_HEADER1_SIZE + DAA_HEADER2_SIZE) as u64;
    Ok(())
}

/// Matches C++ `write_header2(writer, h2)`.
pub fn write_header2<W: Write + Seek>(writer: &mut W, h2: &DaaHeader2) -> io::Result<()> {
    writer.seek(SeekFrom::Start(DAA_HEADER1_SIZE as u64))?;
    h2.write_to(writer)
}

/// Matches C++ `finish_daa(writer, daa_in)`.
pub fn finish_daa_from_file<W: Write + Seek>(writer: &mut W, daa_in: &DaaFile) -> io::Result<()> {
    let mut h2 = DaaHeader2::from_daa_file(daa_in);
    terminate_aln_block(writer, &mut h2)?;

    h2.db_seqs_used = daa_in.db_seqs_used();
    h2.query_records = daa_in.query_records();

    for name in &daa_in.ref_name {
        writer.write_all(name.as_bytes())?;
        writer.write_all(&[0])?;
    }
    h2.block_size[1] = daa_in.block_size(1);

    for &len in &daa_in.ref_len {
        writer.write_all(&len.to_ne_bytes())?;
    }
    h2.block_size[2] = daa_in.block_size(2);

    write_header2(writer, &h2)
}

/// Matches C++ `finish_daa(writer, daa_in, seq_ids, seq_lens, query_count)`.
pub fn finish_daa_from_refs<W: Write + Seek>(
    writer: &mut W,
    daa_in: &DaaFile,
    seq_ids: &[String],
    seq_lens: &[u32],
    query_count: i64,
) -> io::Result<()> {
    let mut h2 = DaaHeader2::from_daa_file(daa_in);
    terminate_aln_block(writer, &mut h2)?;
    h2.db_seqs_used = seq_ids.len() as u64;
    h2.query_records = query_count as u64;

    let mut s = 0u64;
    for id in seq_ids {
        writer.write_all(id.as_bytes())?;
        writer.write_all(&[0])?;
        s += id.len() as u64 + 1;
    }
    h2.block_size[1] = s;

    for &len in seq_lens {
        writer.write_all(&len.to_ne_bytes())?;
    }
    h2.block_size[2] = (seq_lens.len() * std::mem::size_of::<u32>()) as u64;

    write_header2(writer, &h2)
}

fn write_packed(out: &mut Vec<u8>, x: u32) {
    if x <= u8::MAX as u32 {
        out.push(x as u8);
    } else if x <= u16::MAX as u32 {
        out.extend_from_slice(&(x as u16).to_ne_bytes());
    } else {
        out.extend_from_slice(&x.to_ne_bytes());
    }
}

fn write_packed_width(out: &mut Vec<u8>, x: u32, width: u8) {
    match width {
        0 => out.push(x as u8),
        1 => out.extend_from_slice(&(x as u16).to_ne_bytes()),
        2 => out.extend_from_slice(&x.to_ne_bytes()),
        _ => unreachable!("invalid packed integer width"),
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
    use crate::basic::packed_sequence::PackedSequence;
    use crate::basic::packed_transcript::PackedOperation;
    use crate::basic::value::SequenceType;
    use crate::util::interval::Interval;
    use std::io::Cursor;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static DAA_TEST_FILE_ID: AtomicUsize = AtomicUsize::new(0);

    fn daa_file_for_mode(mode: i32) -> (DaaFile, std::path::PathBuf) {
        let id = DAA_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "diamond_rs_daa_record_{}_{}_{}.daa",
            std::process::id(),
            mode,
            id
        ));
        std::fs::write(&path, []).unwrap();
        let mut h2 = DaaHeader2::new();
        h2.mode = mode;
        (
            DaaFile {
                reader: BufReader::new(File::open(&path).unwrap()),
                query_count: 0,
                h1: DaaHeader1::new(),
                h2,
                ref_name: Vec::new(),
                ref_len: Vec::new(),
            },
            path,
        )
    }

    fn write_synthetic_daa(
        path: &std::path::Path,
        refs: &[(&str, u32)],
        query_name: &str,
        query: &[Letter],
        subject_id: u32,
        score: u32,
    ) {
        let (mut daa, daa_path) = daa_file_for_mode(2);
        daa.h2 = DaaHeader2::with_params(
            refs.len() as u64,
            refs.iter().map(|(_, len)| *len as u64).sum(),
            11,
            1,
            2,
            -3,
            0.7,
            1.2,
            0.001,
            "BLOSUM62",
            2,
        );

        let mut f = File::create(path).unwrap();
        init_daa(&mut f).unwrap();

        let mut buf = Vec::new();
        let seek_pos = write_daa_query_record(&mut buf, query_name, query, SequenceType::AminoAcid);
        let flag = compute_flag(score, 0, 0, false);
        buf.extend_from_slice(&subject_id.to_ne_bytes());
        buf.push(flag);
        write_packed_width(&mut buf, score, flag & 3);
        write_packed_width(&mut buf, 0, (flag >> 2) & 3);
        write_packed_width(&mut buf, 0, (flag >> 4) & 3);
        buf.push(
            PackedOperation::from_op_count(
                crate::basic::packed_transcript::EditOperation::Match,
                query.len() as u32,
            )
            .code,
        );
        buf.push(PackedOperation::terminator().code);
        finish_daa_query_record(&mut buf, seek_pos);
        f.write_all(&buf).unwrap();

        let seq_ids = refs
            .iter()
            .map(|(name, _)| (*name).to_string())
            .collect::<Vec<_>>();
        let seq_lens = refs.iter().map(|(_, len)| *len).collect::<Vec<_>>();
        finish_daa_from_refs(&mut f, &daa, &seq_ids, &seq_lens, 1).unwrap();

        std::fs::remove_file(daa_path).unwrap();
    }

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
    fn test_header1_rejects_future_version() {
        let h = DaaHeader1 {
            magic_number: DAA_MAGIC_NUMBER,
            version: DAA_VERSION + 1,
        };
        assert_eq!(
            h.validate().unwrap_err(),
            "DAA version requires later version of DIAMOND."
        );
    }

    #[test]
    fn test_header2_matrix_name() {
        let mut h = DaaHeader2::new();
        h.set_matrix_name("BLOSUM62");
        assert_eq!(h.matrix_name(), "BLOSUM62");
    }

    #[test]
    fn test_header2_with_params_and_from_daa_file() {
        let (mut daa, path) = daa_file_for_mode(2);
        daa.h2 = DaaHeader2::with_params(10, 1_000, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);

        let h = DaaHeader2::from_daa_file(&daa);
        assert_eq!(h.db_seqs, 10);
        assert_eq!(h.db_letters, 1_000);
        assert_eq!(h.gap_open, 11);
        assert_eq!(h.gap_extend, 1);
        assert_eq!(h.reward, 2);
        assert_eq!(h.penalty, -3);
        assert_eq!(h.k, 0.7);
        assert_eq!(h.lambda, 1.2);
        assert_eq!(h.evalue, 0.001);
        assert_eq!(h.matrix_name(), "BLOSUM62");
        assert_eq!(h.mode, 2);
        assert_eq!(h.block_type[0], BlockType::Alignments as u8);
        assert_eq!(h.block_type[1], BlockType::RefNames as u8);
        assert_eq!(h.block_type[2], BlockType::RefLengths as u8);

        std::fs::remove_file(path).unwrap();
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

    #[test]
    fn test_daa_file_open_reads_metadata_refs_and_query_buffers() {
        let mut h2 = DaaHeader2::new();
        h2.diamond_build = 123;
        h2.db_seqs = 4;
        h2.db_seqs_used = 2;
        h2.db_letters = 99;
        h2.query_records = 1;
        h2.mode = 2;
        h2.gap_open = 11;
        h2.gap_extend = 1;
        h2.reward = 2;
        h2.penalty = -3;
        h2.k = 0.7;
        h2.lambda = 1.2;
        h2.evalue = 0.001;
        h2.set_matrix_name("BLOSUM62");

        let query_payload = vec![7u8, 8, 9, 10];
        let align_block_size = 4 + query_payload.len() + 4;
        h2.block_size[0] = align_block_size as u64;
        h2.block_size[1] = 6;
        h2.block_size[2] = 8;
        h2.block_type[0] = BlockType::Alignments as u8;
        h2.block_type[1] = BlockType::RefNames as u8;
        h2.block_type[2] = BlockType::RefLengths as u8;

        let mut bytes = Vec::new();
        DaaHeader1::new().write_to(&mut bytes).unwrap();
        h2.write_to(&mut bytes).unwrap();
        bytes.extend_from_slice(&(query_payload.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&query_payload);
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(b"ref0\0ref1\0");
        bytes.extend_from_slice(&101u32.to_le_bytes());
        bytes.extend_from_slice(&202u32.to_le_bytes());

        let path = std::env::temp_dir().join(format!(
            "diamond_rs_daa_file_{}_{}.daa",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        std::fs::write(&path, &bytes).unwrap();

        let mut daa = DaaFile::open(&path).unwrap();
        assert_eq!(daa.header1().version, DAA_VERSION);
        assert_eq!(daa.diamond_build(), 123);
        assert_eq!(daa.db_seqs(), 4);
        assert_eq!(daa.db_seqs_used(), 2);
        assert_eq!(daa.db_letters(), 99);
        assert_eq!(daa.query_records(), 1);
        assert_eq!(daa.mode(), 2);
        assert_eq!(daa.score_matrix(), "BLOSUM62");
        assert_eq!(daa.gap_open_penalty(), 11);
        assert_eq!(daa.gap_extension_penalty(), 1);
        assert_eq!(daa.match_reward(), 2);
        assert_eq!(daa.mismatch_penalty(), -3);
        assert_eq!(daa.kappa(), 0.7);
        assert_eq!(daa.lambda(), 1.2);
        assert_eq!(daa.evalue(), 0.001);
        assert_eq!(daa.block_size(1), 6);
        assert_eq!(daa.ref_name(0), "ref0");
        assert_eq!(daa.ref_name(1), "ref1");
        assert_eq!(daa.ref_len_at(0), 101);
        assert_eq!(daa.ref_len(), &[101, 202]);
        assert_eq!(daa.read_query_buffer().unwrap(), Some((query_payload, 0)));
        assert_eq!(daa.read_query_buffer().unwrap(), None);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_daa_file_open_rejects_incomplete_run() {
        let mut bytes = Vec::new();
        DaaHeader1::new().write_to(&mut bytes).unwrap();
        DaaHeader2::new().write_to(&mut bytes).unwrap();
        let path = std::env::temp_dir().join(format!(
            "diamond_rs_daa_file_invalid_{}.daa",
            std::process::id()
        ));
        std::fs::write(&path, &bytes).unwrap();

        let err = match DaaFile::open(&path) {
            Ok(_) => panic!("expected incomplete DAA error"),
            Err(err) => err,
        };
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "Invalid DAA file. DIAMOND run has probably not completed successfully."
        );

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_daa_query_record_init_blastp() {
        let (daa, path) = daa_file_for_mode(2);
        let seq = vec![0, 1, 2, 3, 4, 5, 6];
        let packed = PackedSequence::new(&seq, SequenceType::AminoAcid);
        let raw_match = [9u8, 8, 7];

        let mut buf = Vec::new();
        buf.extend_from_slice(&(seq.len() as u32).to_ne_bytes());
        buf.extend_from_slice(b"query1\0");
        buf.push(0);
        buf.extend_from_slice(packed.data());
        buf.extend_from_slice(&raw_match);

        let r = DaaQueryRecord::from_buffer(&daa, &buf, 5).unwrap();
        assert_eq!(r.query_name, "query1");
        assert_eq!(r.query_num, 5);
        assert_eq!(r.query_len(), seq.len());
        assert_eq!(r.context[0], seq);
        assert!(r.source_seq.is_empty());
        assert_eq!(r.raw_matches(), raw_match);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_daa_query_record_init_translated_query() {
        let (daa, path) = daa_file_for_mode(3);
        let seq = vec![0, 3, 2, 4, 0, 1, 2];
        let packed = PackedSequence::new(&seq, SequenceType::Nucleotide);
        let raw_match = [1u8, 2];

        let mut buf = Vec::new();
        buf.extend_from_slice(&(seq.len() as u32).to_ne_bytes());
        buf.extend_from_slice(b"dna_query\0");
        buf.push(if packed.has_n() { 1 } else { 0 });
        buf.extend_from_slice(packed.data());
        buf.extend_from_slice(&raw_match);

        let r = DaaQueryRecord::from_buffer(&daa, &buf, 3).unwrap();
        assert_eq!(r.query_name, "dna_query");
        assert_eq!(r.query_len(), seq.len());
        assert_eq!(r.source_seq, seq);
        assert_eq!(
            r.context,
            translate_6_frames(&seq.iter().map(|&x| x as u8).collect::<Vec<_>>()).map(|frame| {
                frame
                    .into_iter()
                    .map(|letter| letter as Letter)
                    .collect::<Vec<_>>()
            })
        );
        assert_eq!(r.raw_matches(), raw_match);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_daa_raw_match_read_and_copy_match_record_raw() {
        let flag = compute_flag(300, 10, 70_000, false);
        let mut data = Vec::new();
        data.extend_from_slice(&7u32.to_ne_bytes());
        data.push(flag);
        write_packed_width(&mut data, 300, flag & 3);
        write_packed_width(&mut data, 10, (flag >> 2) & 3);
        write_packed_width(&mut data, 70_000, (flag >> 4) & 3);
        data.push(
            PackedOperation::from_op_count(
                crate::basic::packed_transcript::EditOperation::Match,
                3,
            )
            .code,
        );
        data.push(PackedOperation::terminator().code);

        let mut it = BinaryBufferIterator::new(&data);
        let raw = DaaRawMatch::read(&mut it).unwrap();
        assert_eq!(raw.subject_id, 7);
        assert_eq!(raw.score, 300);
        assert_eq!(raw.query_begin, 10);
        assert_eq!(raw.subject_begin, 70_000);
        assert_eq!(raw.transcript.len(), 2);
        assert!(!it.good());

        let mut subject_map = std::collections::HashMap::new();
        subject_map.insert(7, 11);
        let mut copied = Vec::new();
        let mut it = BinaryBufferIterator::new(&data);
        copy_match_record_raw(&mut it, &mut copied, &subject_map).unwrap();

        let mut expected = Vec::new();
        expected.extend_from_slice(&11u32.to_ne_bytes());
        expected.extend_from_slice(&data[4..]);
        assert_eq!(copied, expected);
    }

    #[test]
    fn test_daa_match_iterator_reads_blastp_matches_and_numbers_hsps() {
        let (mut daa, path) = daa_file_for_mode(2);
        daa.h2.db_letters = 1000;
        daa.h2.set_matrix_name("BLOSUM62");
        daa.h2.gap_open = 11;
        daa.h2.gap_extend = 1;
        daa.ref_name = vec!["subject0".to_string()];
        daa.ref_len = vec![100];

        let query = vec![0, 1, 2];
        let packed = PackedSequence::new(&query, SequenceType::AminoAcid);
        let mut buf = Vec::new();
        buf.extend_from_slice(&(query.len() as u32).to_ne_bytes());
        buf.extend_from_slice(b"query0\0");
        buf.push(0);
        buf.extend_from_slice(packed.data());

        for score in [42u32, 40] {
            let flag = compute_flag(score, 0, 0, false);
            buf.extend_from_slice(&0u32.to_ne_bytes());
            buf.push(flag);
            write_packed_width(&mut buf, score, flag & 3);
            write_packed_width(&mut buf, 0, (flag >> 2) & 3);
            write_packed_width(&mut buf, 0, (flag >> 4) & 3);
            buf.push(
                PackedOperation::from_op_count(
                    crate::basic::packed_transcript::EditOperation::Match,
                    query.len() as u32,
                )
                .code,
            );
            buf.push(PackedOperation::terminator().code);
        }

        let record = DaaQueryRecord::from_buffer(&daa, &buf, 7).unwrap();
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, daa.db_letters()).unwrap();
        let matches = record
            .begin(&daa, &score_matrix)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].subject_id, 0);
        assert_eq!(matches[0].subject_name, "subject0");
        assert_eq!(matches[0].hit_num, 0);
        assert_eq!(matches[0].hsp_num, 0);
        assert_eq!(matches[0].hsp.score, 42);
        assert_eq!(matches[0].hsp.frame, 0);
        assert_eq!(matches[0].hsp.query_range, Interval::new(0, 3));
        assert_eq!(matches[0].hsp.subject_range, Interval::new(0, 3));
        assert_eq!(matches[0].hsp.length, 3);
        assert_eq!(matches[0].hsp.identities, 3);
        assert!(matches[0].hsp.evalue.is_finite());
        assert!(matches[0].hsp.bit_score > 0.0);
        assert_eq!(matches[1].hit_num, 0);
        assert_eq!(matches[1].hsp_num, 1);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_daa_match_iterator_reads_translated_match() {
        let (mut daa, path) = daa_file_for_mode(3);
        daa.h2.db_letters = 1000;
        daa.h2.set_matrix_name("BLOSUM62");
        daa.h2.gap_open = 11;
        daa.h2.gap_extend = 1;
        daa.ref_name = vec!["subject0".to_string()];
        daa.ref_len = vec![100];

        let query = vec![0, 1, 2, 0, 1, 2];
        let packed = PackedSequence::new(&query, SequenceType::Nucleotide);
        let mut buf = Vec::new();
        buf.extend_from_slice(&(query.len() as u32).to_ne_bytes());
        buf.extend_from_slice(b"dna0\0");
        buf.push(if packed.has_n() { 1 } else { 0 });
        buf.extend_from_slice(packed.data());

        let flag = compute_flag(30, 0, 5, false);
        buf.extend_from_slice(&0u32.to_ne_bytes());
        buf.push(flag);
        write_packed_width(&mut buf, 30, flag & 3);
        write_packed_width(&mut buf, 0, (flag >> 2) & 3);
        write_packed_width(&mut buf, 5, (flag >> 4) & 3);
        buf.push(
            PackedOperation::from_op_count(
                crate::basic::packed_transcript::EditOperation::Match,
                2,
            )
            .code,
        );
        buf.push(PackedOperation::terminator().code);

        let record = DaaQueryRecord::from_buffer(&daa, &buf, 0).unwrap();
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, daa.db_letters()).unwrap();
        let matches = record
            .begin(&daa, &score_matrix)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].hsp.frame, 0);
        assert_eq!(matches[0].hsp.query_range, Interval::new(0, 2));
        assert_eq!(matches[0].hsp.query_source_range, Interval::new(0, 6));
        assert_eq!(matches[0].hsp.subject_range, Interval::new(5, 7));
        assert_eq!(matches[0].hsp.length, 2);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_output_range_matches_cpp_rules() {
        assert!(output_range(0, 80, 100, 1, None));
        assert!(!output_range(1, 80, 100, 1, None));
        assert!(output_range(7, 90, 100, 1, Some(10.0)));
        assert!(!output_range(7, 89, 100, 1, Some(10.0)));
        assert!(!output_range(0, 1, 0, 25, Some(10.0)));
    }

    #[test]
    fn test_view_query_daa_writes_filtered_daa_record() {
        let (mut daa, path) = daa_file_for_mode(2);
        daa.h2.db_letters = 1000;
        daa.h2.set_matrix_name("BLOSUM62");
        daa.h2.gap_open = 11;
        daa.h2.gap_extend = 1;
        daa.ref_name = vec!["subject0".to_string(), "subject1".to_string()];
        daa.ref_len = vec![100, 100];

        let query = vec![0, 1, 2];
        let packed = PackedSequence::new(&query, SequenceType::AminoAcid);
        let mut buf = Vec::new();
        buf.extend_from_slice(&(query.len() as u32).to_ne_bytes());
        buf.extend_from_slice(b"query0\0");
        buf.push(0);
        buf.extend_from_slice(packed.data());

        for (subject_id, score) in [(0u32, 50u32), (1, 49)] {
            let flag = compute_flag(score, 0, 0, false);
            buf.extend_from_slice(&subject_id.to_ne_bytes());
            buf.push(flag);
            write_packed_width(&mut buf, score, flag & 3);
            write_packed_width(&mut buf, 0, (flag >> 2) & 3);
            write_packed_width(&mut buf, 0, (flag >> 4) & 3);
            buf.push(
                PackedOperation::from_op_count(
                    crate::basic::packed_transcript::EditOperation::Match,
                    query.len() as u32,
                )
                .code,
            );
            buf.push(PackedOperation::terminator().code);
        }

        let record = DaaQueryRecord::from_buffer(&daa, &buf, 0).unwrap();
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, daa.db_letters()).unwrap();
        let mut out = Vec::new();
        view_query_daa(
            &record,
            &daa,
            &mut out,
            &score_matrix,
            &DaaViewConfig {
                max_target_seqs: 1,
                toppercent: None,
                forward_only: false,
            },
        )
        .unwrap();

        let size = u32::from_ne_bytes(out[0..4].try_into().unwrap()) as usize;
        assert_eq!(size, out.len() - 4);
        let viewed = DaaQueryRecord::from_buffer(&daa, &out[4..], 0).unwrap();
        assert_eq!(viewed.query_name, "query0");
        assert_eq!(viewed.context[0], query);
        let mut raw_it = viewed.raw_begin();
        let raw = DaaRawMatch::read(&mut raw_it).unwrap();
        assert_eq!(raw.subject_id, 0);
        assert_eq!(raw.score, 50);
        assert!(!raw_it.good());

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_build_mapping_deduplicates_reference_names() {
        let (mut daa, path) = daa_file_for_mode(2);
        daa.h2.db_seqs_used = 3;
        daa.ref_name = vec!["a".to_string(), "b".to_string(), "a".to_string()];
        daa.ref_len = vec![10, 20, 10];

        let mut acc2oid = HashMap::new();
        let mut seq_ids = Vec::new();
        let mut seq_lens = Vec::new();
        let m = build_mapping(&mut acc2oid, &mut seq_ids, &mut seq_lens, &daa);

        assert_eq!(m.get(&0), Some(&0));
        assert_eq!(m.get(&1), Some(&1));
        assert_eq!(m.get(&2), Some(&0));
        assert_eq!(seq_ids, vec!["a".to_string(), "b".to_string()]);
        assert_eq!(seq_lens, vec![10, 20]);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_write_file_remaps_subjects_and_preserves_query_record() {
        let input_path = std::env::temp_dir().join(format!(
            "diamond_rs_daa_write_file_{}_{}.daa",
            std::process::id(),
            DAA_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed)
        ));
        write_synthetic_daa(
            &input_path,
            &[("ref0", 100), ("ref1", 200)],
            "query one",
            &[0, 1, 2],
            1,
            42,
        );

        let mut daa = DaaFile::open(&input_path).unwrap();
        let mut subject_map = HashMap::new();
        subject_map.insert(1, 7);
        let mut out = Vec::new();
        assert_eq!(write_file(&mut daa, &mut out, &subject_map).unwrap(), 1);

        let size = u32::from_ne_bytes(out[0..4].try_into().unwrap()) as usize;
        assert_eq!(size, out.len() - 4);
        let record = DaaQueryRecord::from_buffer(&daa, &out[4..], 0).unwrap();
        assert_eq!(record.query_name, "query");
        let mut it = record.raw_begin();
        let raw = DaaRawMatch::read(&mut it).unwrap();
        assert_eq!(raw.subject_id, 7);
        assert_eq!(raw.score, 42);
        assert!(!it.good());

        std::fs::remove_file(input_path).unwrap();
    }

    #[test]
    fn test_merge_daa_files_deduplicates_targets_and_appends_queries() {
        let id = DAA_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir();
        let p0 = dir.join(format!(
            "diamond_rs_daa_merge_{}_{}_0.daa",
            std::process::id(),
            id
        ));
        let p1 = dir.join(format!(
            "diamond_rs_daa_merge_{}_{}_1.daa",
            std::process::id(),
            id
        ));
        let out_path = dir.join(format!(
            "diamond_rs_daa_merge_{}_{}_out.daa",
            std::process::id(),
            id
        ));

        write_synthetic_daa(&p0, &[("shared", 100)], "q0", &[0, 1, 2], 0, 30);
        write_synthetic_daa(
            &p1,
            &[("shared", 100), ("new", 200)],
            "q1",
            &[3, 4, 5],
            1,
            40,
        );

        assert_eq!(
            merge_daa_files(&[p0.clone(), p1.clone()], out_path.clone()).unwrap(),
            2
        );

        let mut merged = DaaFile::open(&out_path).unwrap();
        assert_eq!(merged.db_seqs_used(), 2);
        assert_eq!(merged.query_records(), 2);
        assert_eq!(merged.ref_name(0), "shared");
        assert_eq!(merged.ref_name(1), "new");
        assert_eq!(merged.ref_len(), &[100, 200]);

        let (buf0, qn0) = merged.read_query_buffer().unwrap().unwrap();
        let r0 = DaaQueryRecord::from_buffer(&merged, &buf0, qn0).unwrap();
        let raw0 = DaaRawMatch::read(&mut r0.raw_begin()).unwrap();
        assert_eq!(r0.query_name, "q0");
        assert_eq!(raw0.subject_id, 0);

        let (buf1, qn1) = merged.read_query_buffer().unwrap().unwrap();
        let r1 = DaaQueryRecord::from_buffer(&merged, &buf1, qn1).unwrap();
        let raw1 = DaaRawMatch::read(&mut r1.raw_begin()).unwrap();
        assert_eq!(r1.query_name, "q1");
        assert_eq!(raw1.subject_id, 1);
        assert!(merged.read_query_buffer().unwrap().is_none());

        std::fs::remove_file(p0).unwrap();
        std::fs::remove_file(p1).unwrap();
        std::fs::remove_file(out_path).unwrap();
    }

    #[test]
    fn test_write_and_finish_daa_query_record_blastp() {
        let (daa, path) = daa_file_for_mode(2);
        let seq = vec![0, 1, 2, 3, 4, 5, 6];
        let mut buf = Vec::new();

        let seek_pos = write_daa_query_record(
            &mut buf,
            "query1 description",
            &seq,
            SequenceType::AminoAcid,
        );
        finish_daa_query_record(&mut buf, seek_pos);

        let record_size = u32::from_ne_bytes(buf[0..4].try_into().unwrap()) as usize;
        assert_eq!(record_size, buf.len() - 4);
        let r = DaaQueryRecord::from_buffer(&daa, &buf[4..], 0).unwrap();
        assert_eq!(r.query_name, "query1");
        assert_eq!(r.query_len(), seq.len());
        assert_eq!(r.context[0], seq);
        assert!(r.raw_matches().is_empty());

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_write_daa_record_intermediate() {
        let mut record = IntermediateRecord::default();
        record.target_dict_id = 17;
        record.score = 300;
        record.query_begin = 10;
        record.subject_begin = 70_000;
        record.flag = compute_flag(
            record.score,
            record.query_begin,
            record.subject_begin,
            false,
        );
        record
            .transcript
            .push_with_count(crate::basic::packed_transcript::EditOperation::Match, 3);
        record.transcript.push_terminator();

        let mut buf = Vec::new();
        write_daa_record_intermediate(&mut buf, &record);

        let mut it = BinaryBufferIterator::new(&buf);
        let raw = DaaRawMatch::read(&mut it).unwrap();
        assert_eq!(raw.subject_id, 17);
        assert_eq!(raw.flag, record.flag);
        assert_eq!(raw.score, record.score);
        assert_eq!(raw.query_begin, record.query_begin);
        assert_eq!(raw.subject_begin, record.subject_begin);
        assert_eq!(raw.transcript.len(), 2);
        assert!(!it.good());
    }

    #[test]
    fn test_write_daa_record_hsp() {
        let mut hsp = Hsp::new();
        hsp.score = 1_000;
        hsp.frame = 0;
        hsp.query_source_range = Interval::new(12, 20);
        hsp.subject_range = Interval::new(70_000, 70_008);
        hsp.transcript
            .push_with_count(crate::basic::packed_transcript::EditOperation::Match, 8);
        hsp.transcript.push_terminator();

        let mut buf = Vec::new();
        write_daa_record_hsp(&mut buf, &hsp, 23);

        let mut it = BinaryBufferIterator::new(&buf);
        let raw = DaaRawMatch::read(&mut it).unwrap();
        assert_eq!(raw.subject_id, 23);
        assert_eq!(raw.flag, get_segment_flag_hsp(&hsp));
        assert_eq!(raw.score, 1_000);
        assert_eq!(raw.query_begin, 12);
        assert_eq!(raw.subject_begin, 70_000);
        assert_eq!(raw.transcript.len(), 2);
        assert!(!it.good());
    }

    #[test]
    fn test_init_and_finish_daa_from_refs() {
        let (mut daa, path) = daa_file_for_mode(2);
        daa.h2 = DaaHeader2::with_params(4, 99, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);

        let mut cur = Cursor::new(Vec::new());
        init_daa(&mut cur).unwrap();
        cur.write_all(&[1, 2, 3]).unwrap();
        finish_daa_from_refs(
            &mut cur,
            &daa,
            &["ref0".to_string(), "ref1".to_string()],
            &[101, 202],
            5,
        )
        .unwrap();

        let bytes = cur.into_inner();
        let mut header_cursor =
            Cursor::new(&bytes[DAA_HEADER1_SIZE..DAA_HEADER1_SIZE + DAA_HEADER2_SIZE]);
        let h2 = DaaHeader2::read_from(&mut header_cursor).unwrap();
        assert_eq!(h2.db_seqs, 4);
        assert_eq!(h2.db_seqs_used, 2);
        assert_eq!(h2.query_records, 5);
        assert_eq!(h2.block_size[0], 7);
        assert_eq!(h2.block_size[1], b"ref0\0ref1\0".len() as u64);
        assert_eq!(h2.block_size[2], 8);
        assert_eq!(h2.matrix_name(), "BLOSUM62");

        let ref_offset = DAA_HEADER1_SIZE + DAA_HEADER2_SIZE + h2.block_size[0] as usize;
        assert_eq!(
            &bytes[ref_offset..ref_offset + b"ref0\0ref1\0".len()],
            b"ref0\0ref1\0"
        );
        let len_offset = ref_offset + b"ref0\0ref1\0".len();
        assert_eq!(
            u32::from_ne_bytes(bytes[len_offset..len_offset + 4].try_into().unwrap()),
            101
        );
        assert_eq!(
            u32::from_ne_bytes(bytes[len_offset + 4..len_offset + 8].try_into().unwrap()),
            202
        );

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_finish_daa_from_file_copies_existing_refs() {
        let (mut daa, path) = daa_file_for_mode(2);
        daa.h2 = DaaHeader2::with_params(4, 99, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);
        daa.h2.db_seqs_used = 2;
        daa.h2.query_records = 3;
        daa.h2.block_size[1] = b"old0\0old1\0".len() as u64;
        daa.h2.block_size[2] = 8;
        daa.ref_name = vec!["old0".to_string(), "old1".to_string()];
        daa.ref_len = vec![33, 44];

        let mut cur = Cursor::new(Vec::new());
        init_daa(&mut cur).unwrap();
        finish_daa_from_file(&mut cur, &daa).unwrap();

        let bytes = cur.into_inner();
        let mut header_cursor =
            Cursor::new(&bytes[DAA_HEADER1_SIZE..DAA_HEADER1_SIZE + DAA_HEADER2_SIZE]);
        let h2 = DaaHeader2::read_from(&mut header_cursor).unwrap();
        assert_eq!(h2.db_seqs_used, 2);
        assert_eq!(h2.query_records, 3);
        assert_eq!(h2.block_size[0], 4);
        assert_eq!(h2.block_size[1], b"old0\0old1\0".len() as u64);
        assert_eq!(h2.block_size[2], 8);

        let ref_offset = DAA_HEADER1_SIZE + DAA_HEADER2_SIZE + h2.block_size[0] as usize;
        assert_eq!(
            &bytes[ref_offset..ref_offset + b"old0\0old1\0".len()],
            b"old0\0old1\0"
        );

        std::fs::remove_file(path).unwrap();
    }
}
