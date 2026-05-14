use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::basic::translate::{Frame, Strand};
use crate::output::format::print_title;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence::ID_DELIMITERS;

use super::format::format_evalue;

/// Data for a PAF record.
pub struct PafRecord<'a> {
    pub query_id: &'a str,
    pub query_len: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub subject_id: &'a str,
    pub subject_len: i32,
    pub subject_start: i32,
    pub subject_end: i32,
    pub identities: i32,
    pub length: i32,
    pub score: i32,
    pub evalue: f64,
}

/// Write a PAF format record.
pub fn write_paf_record<W: Write>(writer: &mut W, r: &PafRecord<'_>) -> io::Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t255\tAS:i:{}\tEV:f:{}",
        r.query_id,
        r.query_len,
        r.query_start,
        r.query_end,
        r.subject_id,
        r.subject_len,
        r.subject_start,
        r.subject_end,
        r.identities,
        r.length,
        r.score,
        format_evalue(r.evalue),
    )
}

/// Matches C++ PAFFormat::print_match(const HspContext&).
pub fn print_match_context<W: Write>(
    writer: &mut W,
    r: &HspContext,
    score_matrix: &ScoreMatrix,
) -> io::Result<()> {
    let query_id = r
        .query_title
        .split(|c: char| ID_DELIMITERS.contains(c))
        .next()
        .unwrap_or("");
    let target_title = print_title(&r.target_title, false, false, "<>", false);
    let strand = if Frame::from_index(r.frame() as i32).strand == Strand::Forward {
        '+'
    } else {
        '-'
    };

    write!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tAS:i:{}\tZR:i:{}\tZE:f:",
        query_id,
        r.query_len,
        r.query_source_range().begin,
        r.query_source_range().end - 1,
        strand,
        target_title,
        r.subject_len,
        r.subject_range().begin,
        r.subject_range().end - 1,
        r.identities(),
        r.length(),
        score_matrix.bitscore(r.score() as f64) as u32,
        r.score(),
    )?;
    write!(writer, "{}", format_evalue(r.evalue()))?;
    writeln!(writer)
}

/// Matches C++ PAFFormat::print_query_intro().
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_paf() {
        let r = PafRecord {
            query_id: "q1",
            query_len: 100,
            query_start: 0,
            query_end: 50,
            subject_id: "s1",
            subject_len: 200,
            subject_start: 10,
            subject_end: 60,
            identities: 40,
            length: 50,
            score: 100,
            evalue: 1e-10,
        };
        let mut buf = Vec::new();
        write_paf_record(&mut buf, &r).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("q1\t100\t0\t50\t+\ts1\t200\t10\t60\t40\t50\t255"));
        assert!(output.contains("AS:i:100"));
    }

    #[test]
    fn test_print_match_context() {
        use crate::align::hsp::{Hsp, HspContext};
        use crate::util::interval::Interval;

        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1000).unwrap();
        let mut hsp = Hsp::new();
        hsp.score = 80;
        hsp.evalue = 1e-12;
        hsp.frame = 0;
        hsp.identities = 9;
        hsp.length = 10;
        hsp.query_source_range = Interval::new(5, 15);
        hsp.subject_range = Interval::new(20, 30);
        let ctx = HspContext::new(
            hsp,
            0,
            0,
            Vec::new(),
            100,
            "query\x07one",
            0,
            200,
            "target title >extra title",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );

        let mut buf = Vec::new();
        print_match_context(&mut buf, &ctx, &score_matrix).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("query\t100\t5\t14\t+\ttarget\t200\t20\t29\t9\t10\t255"));
        assert!(output.contains("\tZR:i:80\tZE:f:"));
    }

    #[test]
    fn test_print_query_intro() {
        let mut buf = Vec::new();
        print_query_intro(&mut buf, "query\x0bone", true).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "query\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"
        );

        let mut buf = Vec::new();
        print_query_intro(&mut buf, "query one", false).unwrap();
        assert!(buf.is_empty());
    }
}
