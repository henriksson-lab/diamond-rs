use std::collections::HashMap;
use std::io::{self, Read, Write};

use crate::align::hsp::Hsp;
use crate::basic::packed_transcript::{PackedOperation, PackedTranscript};
use crate::basic::value::{AlignMode, DictId, OId};
use crate::dp::swipe::HspValues;
use crate::output::format::FormatCode;
use crate::output::{get_segment_flag_hsp, INTERMEDIATE_RECORD_SEED_ONLY};
use crate::util::interval::Interval;

pub const FINISHED: u32 = u32::MAX;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OutputFormat {
    pub code: FormatCode,
    pub hsp_values: HspValues,
}

impl OutputFormat {
    pub fn new(code: FormatCode, hsp_values: HspValues) -> Self {
        Self { code, hsp_values }
    }

    pub fn daa() -> Self {
        Self::new(FormatCode::Daa, HspValues::TRANSCRIPT)
    }

    pub fn tabular() -> Self {
        Self::new(FormatCode::Tabular, HspValues::TRANSCRIPT)
    }

    pub fn edge() -> Self {
        Self::new(FormatCode::Tabular, HspValues::COORDS)
    }
}

#[derive(Debug, Clone)]
pub struct IntermediateRecord {
    pub query_id: u32,
    pub target_dict_id: DictId,
    pub target_oid: OId,
    pub score: u32,
    pub query_begin: u32,
    pub subject_begin: u32,
    pub query_end: u32,
    pub subject_end: u32,
    pub identities: u32,
    pub mismatches: u32,
    pub positives: u32,
    pub length: u32,
    pub gap_openings: u32,
    pub gaps: u32,
    pub evalue: f64,
    pub flag: u8,
    pub transcript: PackedTranscript,
}

impl Default for IntermediateRecord {
    fn default() -> Self {
        Self {
            query_id: 0,
            target_dict_id: 0,
            target_oid: 0,
            score: 0,
            query_begin: 0,
            subject_begin: 0,
            query_end: 0,
            subject_end: 0,
            identities: 0,
            mismatches: 0,
            positives: 0,
            length: 0,
            gap_openings: 0,
            gaps: 0,
            evalue: 0.0,
            flag: 0,
            transcript: PackedTranscript::new(),
        }
    }
}

impl IntermediateRecord {
    pub const SEED_ONLY: u8 = INTERMEDIATE_RECORD_SEED_ONLY;
    pub const FINISHED: u32 = FINISHED;

    pub fn read<R: Read>(reader: &mut R, output_format: &OutputFormat) -> io::Result<Self> {
        let mut record = Self::default();
        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8)?;
        record.target_dict_id = i64::from_ne_bytes(buf8);
        if output_format.code == FormatCode::Daa {
            reader.read_exact(&mut buf8)?;
            record.target_oid = u64::from_ne_bytes(buf8);
        }
        let mut flag = [0u8; 1];
        reader.read_exact(&mut flag)?;
        record.flag = flag[0];
        record.score = read_packed(reader, record.flag & 3)?;
        reader.read_exact(&mut buf8)?;
        record.evalue = f64::from_ne_bytes(buf8);

        if output_format.hsp_values == HspValues::NONE {
            return Ok(record);
        }

        record.query_begin = read_packed(reader, (record.flag >> 2) & 3)?;
        record.query_end = read_varint(reader)?;
        record.subject_begin = read_packed(reader, (record.flag >> 4) & 3)?;

        if record.seed_only() {
            record.subject_end = read_varint(reader)?;
            return Ok(record);
        }
        if output_format.hsp_values.any(HspValues::TRANSCRIPT) {
            let mut bytes = Vec::new();
            loop {
                let mut code = [0u8; 1];
                reader.read_exact(&mut code)?;
                bytes.push(code[0]);
                if PackedOperation::new(code[0]).is_terminator() {
                    break;
                }
            }
            record.transcript = PackedTranscript::from_bytes(&bytes);
        } else {
            record.subject_end = read_varint(reader)?;
            record.identities = read_varint(reader)?;
            record.mismatches = read_varint(reader)?;
            record.positives = read_varint(reader)?;
            record.length = read_varint(reader)?;
            record.gap_openings = read_varint(reader)?;
            record.gaps = read_varint(reader)?;
        }
        Ok(record)
    }

    pub fn frame(&self, query_source_len: i32, align_mode: i32) -> u32 {
        if align_mode == AlignMode::BLASTX {
            if (self.flag & (1 << 6)) == 0 {
                self.query_begin % 3
            } else {
                3 + ((query_source_len as u32 - 1 - self.query_begin) % 3)
            }
        } else {
            0
        }
    }

    pub fn absolute_query_range(&self) -> Interval {
        if self.query_begin < self.query_end {
            Interval::new(self.query_begin as i32, self.query_end as i32 + 1)
        } else {
            Interval::new(self.query_end as i32, self.query_begin as i32 + 1)
        }
    }

    pub fn write_query_intro(buf: &mut Vec<u8>, query_id: u32) -> usize {
        let seek_pos = buf.len();
        buf.extend_from_slice(&query_id.to_ne_bytes());
        buf.extend_from_slice(&0u32.to_ne_bytes());
        seek_pos
    }

    pub fn finish_query(buf: &mut [u8], seek_pos: usize) {
        let n = (buf.len() - seek_pos - std::mem::size_of::<u32>() * 2) as u32;
        buf[seek_pos + std::mem::size_of::<u32>()..seek_pos + std::mem::size_of::<u32>() * 2]
            .copy_from_slice(&n.to_ne_bytes());
    }

    pub fn write<W: Write>(
        writer: &mut W,
        m: &Hsp,
        _query_id: u32,
        target: DictId,
        target_oid: OId,
        output_format: &OutputFormat,
    ) -> io::Result<()> {
        let oriented_range = m.oriented_range();
        writer.write_all(&target.to_ne_bytes())?;
        if output_format.code == FormatCode::Daa {
            writer.write_all(&target_oid.to_ne_bytes())?;
        }
        writer.write_all(&[get_segment_flag_hsp(m)])?;
        write_packed(writer, m.score as u32)?;
        writer.write_all(&m.evalue.to_ne_bytes())?;
        if output_format.hsp_values == HspValues::NONE {
            return Ok(());
        }

        write_packed(writer, oriented_range.begin as u32)?;
        write_varint(writer, oriented_range.end as u32)?;
        write_packed(writer, m.subject_range.begin as u32)?;
        if m.seed_only {
            write_varint(writer, m.subject_range.end as u32)?;
            return Ok(());
        }

        if output_format.hsp_values.any(HspValues::TRANSCRIPT) {
            for op in m.transcript.data() {
                writer.write_all(&[op.code])?;
            }
        } else {
            write_varint(writer, m.subject_range.end as u32)?;
            write_varint(writer, m.identities as u32)?;
            write_varint(writer, m.mismatches as u32)?;
            write_varint(writer, m.positives as u32)?;
            write_varint(writer, m.length as u32)?;
            write_varint(writer, m.gap_openings as u32)?;
            write_varint(writer, m.gaps as u32)?;
        }
        Ok(())
    }

    pub fn write_target_score<W: Write>(
        writer: &mut W,
        target_oid: u32,
        score: i32,
    ) -> io::Result<()> {
        writer.write_all(&target_oid.to_ne_bytes())?;
        let s = score.min(u16::MAX as i32) as u16;
        writer.write_all(&s.to_ne_bytes())
    }

    pub fn finish_file<W: Write>(writer: &mut W) -> io::Result<()> {
        writer.write_all(&FINISHED.to_ne_bytes())
    }

    pub fn stats_mode(v: HspValues) -> bool {
        !v.any(HspValues::TRANSCRIPT) && v != HspValues::NONE
    }

    pub fn seed_only(&self) -> bool {
        (self.flag & Self::SEED_ONLY) != 0
    }
}

pub fn read_packed<R: Read>(reader: &mut R, length: u8) -> io::Result<u32> {
    match length {
        0 => {
            let mut b = [0u8; 1];
            reader.read_exact(&mut b)?;
            Ok(b[0] as u32)
        }
        1 => {
            let mut b = [0u8; 2];
            reader.read_exact(&mut b)?;
            Ok(u16::from_ne_bytes(b) as u32)
        }
        2 => {
            let mut b = [0u8; 4];
            reader.read_exact(&mut b)?;
            Ok(u32::from_ne_bytes(b))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid packed integer width",
        )),
    }
}

pub fn write_packed<W: Write>(writer: &mut W, x: u32) -> io::Result<()> {
    if x <= u8::MAX as u32 {
        writer.write_all(&[x as u8])
    } else if x <= u16::MAX as u32 {
        writer.write_all(&(x as u16).to_ne_bytes())
    } else {
        writer.write_all(&x.to_ne_bytes())
    }
}

pub fn read_varint<R: Read>(reader: &mut R) -> io::Result<u32> {
    let mut shift = 0;
    let mut dst = 0u32;
    loop {
        let mut b = [0u8; 1];
        reader.read_exact(&mut b)?;
        dst |= ((b[0] & 0x7f) as u32) << shift;
        if b[0] & 0x80 == 0 {
            return Ok(dst);
        }
        shift += 7;
        if shift >= 35 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid uint32 varint",
            ));
        }
    }
}

pub fn write_varint<W: Write>(writer: &mut W, mut x: u32) -> io::Result<()> {
    while x >= 0x80 {
        writer.write_all(&[((x as u8) & 0x7f) | 0x80])?;
        x >>= 7;
    }
    writer.write_all(&[x as u8])
}

pub fn copy_match_record_raw<R: Read, W: Write>(
    reader: &mut R,
    writer: &mut W,
    subject_map: &HashMap<u32, u32>,
) -> io::Result<()> {
    let mut b4 = [0u8; 4];
    reader.read_exact(&mut b4)?;
    let subject_id = u32::from_ne_bytes(b4);
    let mut flag = [0u8; 1];
    reader.read_exact(&mut flag)?;
    let score = read_packed(reader, flag[0] & 3)?;
    let query_begin = read_packed(reader, (flag[0] >> 2) & 3)?;
    let subject_begin = read_packed(reader, (flag[0] >> 4) & 3)?;
    let mapped_subject = subject_map
        .get(&subject_id)
        .copied()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "subject id missing from map"))?;
    writer.write_all(&mapped_subject.to_ne_bytes())?;
    writer.write_all(&flag)?;
    write_packed(writer, score)?;
    write_packed(writer, query_begin)?;
    write_packed(writer, subject_begin)?;
    loop {
        let mut code = [0u8; 1];
        reader.read_exact(&mut code)?;
        writer.write_all(&code)?;
        if PackedOperation::new(code[0]).is_terminator() {
            return Ok(());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::EditOperation;

    #[test]
    fn test_varint_roundtrip_boundaries() {
        for x in [0, 1, 127, 128, 16_383, 16_384, u32::MAX] {
            let mut buf = Vec::new();
            write_varint(&mut buf, x).unwrap();
            let mut cur = &buf[..];
            assert_eq!(read_varint(&mut cur).unwrap(), x);
            assert!(cur.is_empty());
        }
    }

    #[test]
    fn test_query_intro_and_finish_query() {
        let mut buf = Vec::new();
        let seek = IntermediateRecord::write_query_intro(&mut buf, 17);
        buf.extend_from_slice(&[1, 2, 3, 4, 5]);
        IntermediateRecord::finish_query(&mut buf, seek);
        assert_eq!(u32::from_ne_bytes(buf[0..4].try_into().unwrap()), 17);
        assert_eq!(u32::from_ne_bytes(buf[4..8].try_into().unwrap()), 5);
    }

    #[test]
    fn test_finish_file() {
        let mut buf = Vec::new();
        IntermediateRecord::finish_file(&mut buf).unwrap();
        assert_eq!(u32::from_ne_bytes(buf.try_into().unwrap()), FINISHED);
    }

    #[test]
    fn test_write_read_stats_record() {
        let mut hsp = Hsp::new();
        hsp.score = 300;
        hsp.evalue = 1.25e-12;
        hsp.frame = 1;
        hsp.query_source_range = Interval::new(10, 30);
        hsp.subject_range = Interval::new(700, 740);
        hsp.identities = 18;
        hsp.mismatches = 2;
        hsp.positives = 21;
        hsp.length = 30;
        hsp.gap_openings = 1;
        hsp.gaps = 3;
        let format = OutputFormat::new(
            FormatCode::Tabular,
            HspValues::COORDS
                | HspValues::IDENT
                | HspValues::MISMATCHES
                | HspValues::LENGTH
                | HspValues::GAP_OPENINGS,
        );
        let mut buf = Vec::new();
        IntermediateRecord::write(&mut buf, &hsp, 0, 42, 0, &format).unwrap();
        let mut cur = &buf[..];
        let r = IntermediateRecord::read(&mut cur, &format).unwrap();
        assert_eq!(r.target_dict_id, 42);
        assert_eq!(r.score, 300);
        assert_eq!(r.query_begin, 10);
        assert_eq!(r.query_end, 29);
        assert_eq!(r.subject_begin, 700);
        assert_eq!(r.subject_end, 740);
        assert_eq!(r.identities, 18);
        assert_eq!(r.mismatches, 2);
        assert_eq!(r.positives, 21);
        assert_eq!(r.length, 30);
        assert_eq!(r.gap_openings, 1);
        assert_eq!(r.gaps, 3);
        assert_eq!(r.evalue, hsp.evalue);
    }

    #[test]
    fn test_write_read_daa_transcript_record() {
        let mut hsp = Hsp::new();
        hsp.score = 42;
        hsp.evalue = 0.5;
        hsp.frame = 4;
        hsp.query_source_range = Interval::new(21, 3);
        hsp.subject_range = Interval::new(5, 9);
        hsp.transcript.push_with_count(EditOperation::Match, 3);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 2);
        hsp.transcript.push_terminator();
        let format = OutputFormat::daa();
        let mut buf = Vec::new();
        IntermediateRecord::write(&mut buf, &hsp, 11, -7, 1234, &format).unwrap();
        let mut cur = &buf[..];
        let r = IntermediateRecord::read(&mut cur, &format).unwrap();
        assert_eq!(r.target_dict_id, -7);
        assert_eq!(r.target_oid, 1234);
        assert_eq!(r.flag & (1 << 6), 1 << 6);
        assert_eq!(r.query_begin, 2);
        assert_eq!(r.query_end, 21);
        assert_eq!(r.subject_begin, 5);
        assert_eq!(
            r.transcript
                .data()
                .iter()
                .map(|op| op.code)
                .collect::<Vec<_>>(),
            hsp.transcript
                .data()
                .iter()
                .map(|op| op.code)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_write_read_seed_only_record() {
        let mut hsp = Hsp::new();
        hsp.seed_only = true;
        hsp.score = 7;
        hsp.evalue = 10.0;
        hsp.query_source_range = Interval::new(4, 9);
        hsp.subject_range = Interval::new(1, 6);
        let format = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let mut buf = Vec::new();
        IntermediateRecord::write(&mut buf, &hsp, 0, 2, 0, &format).unwrap();
        let mut cur = &buf[..];
        let r = IntermediateRecord::read(&mut cur, &format).unwrap();
        assert!(r.seed_only());
        assert_eq!(r.subject_end, 6);
        assert_eq!(
            r.flag & INTERMEDIATE_RECORD_SEED_ONLY,
            INTERMEDIATE_RECORD_SEED_ONLY
        );
    }

    #[test]
    fn test_frame_and_absolute_query_range() {
        let mut r = IntermediateRecord {
            query_begin: 5,
            query_end: 11,
            ..IntermediateRecord::default()
        };
        assert_eq!(r.absolute_query_range(), Interval::new(5, 12));
        assert_eq!(r.frame(100, AlignMode::BLASTP), 0);
        assert_eq!(r.frame(100, AlignMode::BLASTX), 2);
        r.flag = 1 << 6;
        r.query_begin = 8;
        r.query_end = 3;
        assert_eq!(r.absolute_query_range(), Interval::new(3, 9));
        assert_eq!(r.frame(100, AlignMode::BLASTX), 4);
    }

    #[test]
    fn test_write_target_score() {
        let mut buf = Vec::new();
        IntermediateRecord::write_target_score(&mut buf, 9, i32::MAX).unwrap();
        assert_eq!(u32::from_ne_bytes(buf[0..4].try_into().unwrap()), 9);
        assert_eq!(u16::from_ne_bytes(buf[4..6].try_into().unwrap()), u16::MAX);
    }

    #[test]
    fn test_copy_match_record_raw_remaps_subject_and_copies_transcript() {
        let subject_id = 11u32;
        let flag = 1 | (1 << 2) | (0 << 4);
        let mut input = Vec::new();
        input.extend_from_slice(&subject_id.to_ne_bytes());
        input.push(flag);
        write_packed(&mut input, 300).unwrap();
        write_packed(&mut input, 700).unwrap();
        write_packed(&mut input, 5).unwrap();
        input.push(PackedOperation::from_op_count(EditOperation::Match, 3).code);
        input.push(PackedOperation::terminator().code);
        let subject_map = HashMap::from([(11, 4)]);
        let mut output = Vec::new();
        let mut cur = &input[..];
        copy_match_record_raw(&mut cur, &mut output, &subject_map).unwrap();
        let mut out = &output[..];
        let mut b4 = [0u8; 4];
        out.read_exact(&mut b4).unwrap();
        assert_eq!(u32::from_ne_bytes(b4), 4);
        let mut copied_flag = [0u8; 1];
        out.read_exact(&mut copied_flag).unwrap();
        assert_eq!(copied_flag[0], flag);
        assert_eq!(read_packed(&mut out, copied_flag[0] & 3).unwrap(), 300);
        assert_eq!(
            read_packed(&mut out, (copied_flag[0] >> 2) & 3).unwrap(),
            700
        );
        assert_eq!(read_packed(&mut out, (copied_flag[0] >> 4) & 3).unwrap(), 5);
        assert_eq!(out, &[3, 0]);
    }
}
