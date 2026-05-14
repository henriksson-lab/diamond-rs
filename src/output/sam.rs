use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::basic::consts::VERSION_STRING;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::translate::Frame;
use crate::basic::value::{Letter, AMINO_ACID_ALPHABET, LETTER_MASK};
use crate::output::format::{format_evalue, print_title};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence::ID_DELIMITERS;

/// Data for a SAM record.
pub struct SamRecord<'a> {
    pub query_id: &'a str,
    pub subject_id: &'a str,
    pub subject_start: i32,
    pub query_seq: &'a [Letter],
    pub score: i32,
    pub evalue: f64,
    pub identities: i32,
    pub length: i32,
    pub mismatches: i32,
}

/// Write a SAM format record for a protein alignment.
pub fn write_sam_record<W: Write>(writer: &mut W, r: &SamRecord<'_>) -> io::Result<()> {
    let cigar = format!("{}M", r.length);
    let seq_str: String = r
        .query_seq
        .iter()
        .map(|&l| AMINO_ACID_ALPHABET[(l & LETTER_MASK) as usize] as char)
        .collect();

    writeln!(
        writer,
        "{}\t0\t{}\t{}\t255\t{}\t*\t0\t0\t{}\t*\tAS:i:{}\tNM:i:{}\tZR:i:{}\tZE:f:{:.2e}",
        r.query_id,
        r.subject_id,
        r.subject_start + 1,
        cigar,
        seq_str,
        r.score,
        r.mismatches,
        r.identities,
        r.evalue,
    )
}

pub fn print_md(r: &HspContext) -> String {
    let mut buf = String::new();
    let mut matches = 0u32;
    let mut del = 0u32;
    for i in r.begin_old() {
        match i.op {
            EditOperation::Match => {
                del = 0;
                matches += i.count;
            }
            EditOperation::Substitution => {
                if matches > 0 {
                    buf.push_str(&matches.to_string());
                    matches = 0;
                } else if del > 0 {
                    buf.push('0');
                    del = 0;
                }
                buf.push(AMINO_ACID_ALPHABET[i.letter as usize] as char);
            }
            EditOperation::Deletion => {
                if matches > 0 {
                    buf.push_str(&matches.to_string());
                    matches = 0;
                }
                if del == 0 {
                    buf.push('^');
                }
                buf.push(AMINO_ACID_ALPHABET[i.letter as usize] as char);
                del += 1;
            }
            _ => {}
        }
    }
    if matches > 0 {
        buf.push_str(&matches.to_string());
    }
    buf
}

pub fn print_cigar(r: &HspContext) -> String {
    let map = [0usize, 1, 2, 0, 3, 4];
    let letter = ['M', 'I', 'D', '\\', '/'];
    let mut buf = String::new();
    let mut n = 0u32;
    let mut op = 0usize;
    for i in r.begin_old() {
        let mapped = map[i.op as usize];
        if mapped == op {
            n += i.count;
        } else {
            if n > 0 {
                buf.push_str(&n.to_string());
                buf.push(letter[op]);
            }
            n = i.count;
            op = mapped;
        }
    }
    if n > 0 {
        buf.push_str(&n.to_string());
        buf.push(letter[op]);
    }
    buf
}

/// Matches C++ SamFormat::print_match(const HspContext&).
pub fn print_match_context<W: Write>(
    writer: &mut W,
    r: &HspContext,
    score_matrix: &ScoreMatrix,
    sam_qlen_field: bool,
) -> io::Result<()> {
    let query_id = r
        .query_title
        .split(|c: char| ID_DELIMITERS.contains(c))
        .next()
        .unwrap_or("");
    let target_title = print_title(&r.target_title, true, true, "<>", false);
    let frame = Frame::from_index(r.frame() as i32);
    let query_range = r.query_range();
    let query = r.query_index(r.hsp.frame);
    let begin = query_range.begin.max(0) as usize;
    let end = query_range.end.max(query_range.begin).max(0) as usize;
    let seq_str: String = query
        .get(begin..end.min(query.len()))
        .unwrap_or(&[])
        .iter()
        .map(|&l| AMINO_ACID_ALPHABET[(l & LETTER_MASK) as usize] as char)
        .collect();

    write!(
        writer,
        "{}\t0\t{}\t{}\t255\t{}\t*\t0\t0\t{}\t*\tAS:i:{}\tNM:i:{}\tZL:i:{}\tZR:i:{}\tZE:f:{}\tZI:i:{}\tZF:i:{}\tZS:i:{}\tMD:Z:{}",
        query_id,
        target_title,
        r.subject_range().begin + 1,
        print_cigar(r),
        seq_str,
        score_matrix.bitscore(r.score() as f64) as u32,
        r.length() - r.identities(),
        r.subject_len,
        r.score(),
        format_evalue(r.evalue()),
        r.identities() * 100 / r.length(),
        frame.signed_frame(),
        r.oriented_query_range().begin + 1,
        print_md(r),
    )?;
    if sam_qlen_field {
        write!(writer, "\tZQ:i:{}", r.query_source_len)?;
    }
    writeln!(writer)
}

/// Matches C++ SamFormat::print_query_intro().
pub fn print_query_intro<W: Write>(
    writer: &mut W,
    query_title: &str,
    unaligned: bool,
) -> io::Result<()> {
    if unaligned {
        let query_id = query_title
            .split(|c: char| ID_DELIMITERS.contains(c))
            .next()
            .unwrap_or("");
        writeln!(writer, "{}\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*", query_id)?;
    }
    Ok(())
}

/// Matches C++ SamFormat::print_header().
pub fn print_header<W: Write>(writer: &mut W, mode: u32, invocation: &str) -> io::Result<()> {
    let mode_str = match mode {
        2 => "BlastP",
        3 => "BlastX",
        4 => "BlastN",
        _ => "",
    };
    writeln!(writer, "@HD\tVN:1.5\tSO:query")?;
    writeln!(
        writer,
        "@PG\tPN:DIAMOND\tVN:{}\tCL:{}",
        VERSION_STRING, invocation
    )?;
    writeln!(writer, "@mm\t{}", mode_str)?;
    writeln!(writer, "@CO\t{}-like alignments", mode_str)?;
    writeln!(
        writer,
        "@CO\tReporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate"
    )?;
    Ok(())
}

/// Write SAM header.
pub fn write_sam_header<W: Write>(
    writer: &mut W,
    ref_names: &[&str],
    ref_lengths: &[i32],
) -> io::Result<()> {
    writeln!(writer, "@HD\tVN:1.5\tSO:queryname")?;
    for (name, &len) in ref_names.iter().zip(ref_lengths.iter()) {
        writeln!(writer, "@SQ\tSN:{}\tLN:{}", name, len)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::{Hsp, HspContext};
    use crate::basic::packed_transcript::EditOperation;

    #[test]
    fn test_write_sam_record() {
        let query = vec![0i8, 1, 2, 3];
        let r = SamRecord {
            query_id: "q1",
            subject_id: "s1",
            subject_start: 0,
            query_seq: &query,
            score: 50,
            evalue: 1e-10,
            identities: 4,
            length: 4,
            mismatches: 0,
        };
        let mut buf = Vec::new();
        write_sam_record(&mut buf, &r).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("q1"));
        assert!(output.contains("ARND"));
        assert!(output.contains("4M"));
        assert!(output.contains("AS:i:50"));
    }

    #[test]
    fn test_write_sam_header() {
        let mut buf = Vec::new();
        write_sam_header(&mut buf, &["seq1", "seq2"], &[100, 200]).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("@HD\tVN:1.5"));
        assert!(output.contains("@SQ\tSN:seq1\tLN:100"));
    }

    #[test]
    fn test_print_header() {
        let mut buf = Vec::new();
        print_header(&mut buf, 2, "diamond view -f sam").unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("@HD\tVN:1.5\tSO:query\n"));
        assert!(output.contains("@PG\tPN:DIAMOND\tVN:2.1.24\tCL:diamond view -f sam\n"));
        assert!(output.contains("@mm\tBlastP\n"));
    }

    #[test]
    fn test_print_md_and_cigar_from_transcript() {
        let mut hsp = Hsp::new();
        hsp.transcript.push_with_count(EditOperation::Match, 3);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 1);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 2);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 3);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 4);
        hsp.transcript.push_with_count(EditOperation::Insertion, 2);
        hsp.transcript.push(EditOperation::FrameshiftForward);
        hsp.transcript.push(EditOperation::FrameshiftReverse);

        let ctx = HspContext::new(
            hsp,
            0,
            0,
            Vec::new(),
            0,
            "",
            0,
            0,
            "",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );

        assert_eq!(print_md(&ctx), "3R^ND0C");
        assert_eq!(print_cigar(&ctx), "4M2D1M2I1\\1/");
    }

    #[test]
    fn test_print_query_intro() {
        let mut buf = Vec::new();
        print_query_intro(&mut buf, "query\x08one", true).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "query\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"
        );

        let mut buf = Vec::new();
        print_query_intro(&mut buf, "query one", false).unwrap();
        assert!(buf.is_empty());
    }

    #[test]
    fn test_print_match_context() {
        use crate::util::interval::Interval;

        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1000).unwrap();
        let mut hsp = Hsp::new();
        hsp.score = 80;
        hsp.evalue = 1e-12;
        hsp.frame = 0;
        hsp.identities = 3;
        hsp.length = 4;
        hsp.query_range = Interval::new(1, 5);
        hsp.query_source_range = Interval::new(1, 5);
        hsp.subject_range = Interval::new(20, 24);
        hsp.transcript.push_with_count(EditOperation::Match, 4);
        let ctx = HspContext::new(
            hsp,
            0,
            0,
            vec![vec![0, 1, 2, 3, 4, 5]],
            6,
            "query\x0cone",
            0,
            200,
            "target title\x01extra title",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );

        let mut buf = Vec::new();
        print_match_context(&mut buf, &ctx, &score_matrix, true).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("query\t0\ttarget title<>extra title\t21\t255\t4M"));
        assert!(output.contains("\tZL:i:200\tZR:i:80\t"));
        assert!(output.contains("\tZI:i:75\tZF:i:1\tZS:i:2\tMD:Z:4\tZQ:i:6\n"));
    }
}
