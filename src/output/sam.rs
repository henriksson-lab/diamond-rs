use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::basic::consts::VERSION_STRING;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::translate::Frame;
use crate::basic::value::{Letter, AMINO_ACID_ALPHABET};
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
        .map(|&l| {
            let x = l as i32;
            let c = AMINO_ACID_ALPHABET[(x & 127) as usize] as char;
            if (x & 128) == 0 {
                c
            } else {
                c.to_ascii_lowercase()
            }
        })
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

pub fn print_md<T>(r: &HspContext<T>) -> String {
    let mut buf = String::new();
    let mut matches = 0u32;
    let mut del = 0u32;
    for i in r.begin_old() {
        match i.op {
            EditOperation::Match => {
                del = 0;
                matches = matches.wrapping_add(i.count);
            }
            EditOperation::Substitution => {
                if matches > 0 {
                    buf.push_str(&matches.to_string());
                    matches = 0;
                } else if del > 0 {
                    buf.push('0');
                    del = 0;
                }
                buf.push(r.alphabet()[i.letter as usize] as char);
            }
            EditOperation::Deletion => {
                if matches > 0 {
                    buf.push_str(&matches.to_string());
                    matches = 0;
                }
                if del == 0 {
                    buf.push('^');
                }
                buf.push(r.alphabet()[i.letter as usize] as char);
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

pub fn print_cigar<T>(r: &HspContext<T>) -> String {
    let map = [0usize, 1, 2, 0, 3, 4];
    let letter = ['M', 'I', 'D', '\\', '/'];
    let mut buf = String::new();
    let mut n = 0u32;
    let mut op = 0usize;
    for i in r.begin_old() {
        let mapped = map[i.op as usize];
        if mapped == op {
            n = n.wrapping_add(i.count);
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

/// Matches C++ `SamFormat::print_match(const HspContext&)`.
pub fn print_match_context_bytes<W: Write, T: AsRef<[u8]>>(
    writer: &mut W,
    r: &HspContext<T>,
    score_matrix: &ScoreMatrix,
    full_titles: bool,
    all_titles: bool,
    sam_qlen_field: bool,
) -> io::Result<()> {
    let query_id = crate::util::sequence::seqid_bytes(r.query_title.as_ref());
    let target_title = super::format::print_title_bytes(
        r.target_title.as_ref(),
        full_titles,
        all_titles,
        b"<>",
        false,
    );
    let sam_bit_score = score_matrix.bitscore(r.score() as f64);
    if score_matrix.source_dna_statistics().is_some() && !sam_bit_score.is_finite() {
        // Source's static blastn global yields NaN, then SAM casts it to u32.
        // Reject that actual float-to-integer UB boundary rather than emit a
        // fabricated AS score from the independent nucleotide HSP model.
        return Err(io::Error::other(
            "Nucleotide SAM AS tag converts a non-finite Source global bit score to an integer.",
        ));
    }
    let frame = Frame::from_index(r.frame() as i32);
    let query_range = r.query_range();
    let query = r.query_index(r.hsp.frame);
    if query_range.begin < 0
        || query_range.end < query_range.begin
        || query_range.end as usize > query.len()
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "query range out of bounds",
        ));
    }
    let begin = query_range.begin as usize;
    let end = query_range.end as usize;
    let seq_str: String = query
        .get(begin..end)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "query range out of bounds"))?
        .iter()
        .map(|&l| {
            let x = l as i32;
            let c = r.alphabet()[(x & 127) as usize] as char;
            if (x & 128) == 0 {
                c
            } else {
                c.to_ascii_lowercase()
            }
        })
        .collect();

    writer.write_all(&query_id)?;
    writer.write_all(b"\t0\t")?;
    writer.write_all(&target_title)?;
    write!(
        writer,
        "\t{}\t255\t{}\t*\t0\t0\t{}\t*\tAS:i:{}\tNM:i:{}\tZL:i:{}\tZR:i:{}\tZE:f:{}\tZI:i:{}\tZF:i:{}\tZS:i:{}\tMD:Z:{}",
        r.subject_range().begin + 1,
        print_cigar(r),
        seq_str,
        sam_bit_score as u32,
        r.length().wrapping_sub(r.identities()),
        r.subject_len,
        r.score(),
        format_evalue(r.evalue()),
        r.identities().wrapping_mul(100) / r.length(),
        frame.signed_frame(),
        r.oriented_query_range().begin + 1,
        print_md(r),
    )?;
    if sam_qlen_field {
        write!(writer, "\tZQ:i:{}", r.query_source_len)?;
    }
    writeln!(writer)
}

/// Matches C++ `SamFormat::print_query_intro()`.
pub fn print_query_intro_bytes<W: Write>(
    writer: &mut W,
    query_title: &[u8],
    unaligned: bool,
) -> io::Result<()> {
    if unaligned {
        let query_id = crate::util::sequence::seqid_bytes(query_title);
        writer.write_all(&query_id)?;
        writeln!(writer, "\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*")?;
    }
    Ok(())
}

/// Matches C++ `SamFormat::print_header()`.
pub fn print_header<W: Write>(writer: &mut W, mode: u32, invocation: &str) -> io::Result<()> {
    print_header_bytes(writer, mode, invocation.as_bytes())
}

pub fn print_header_bytes<W: Write>(
    writer: &mut W,
    mode: u32,
    invocation: &[u8],
) -> io::Result<()> {
    let mode_str = match mode {
        2 => "BlastP",
        3 => "BlastX",
        4 => "BlastN",
        _ => "",
    };
    writeln!(writer, "@HD\tVN:1.5\tSO:query")?;
    write!(writer, "@PG\tPN:DIAMOND\tVN:{}\tCL:", VERSION_STRING)?;
    writer.write_all(invocation)?;
    writer.write_all(b"\n")?;
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
        let query = vec![0i8, -127, 2, 3];
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
        assert!(output.contains("ArND"));
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
        assert_eq!(
            output,
            concat!(
                "@HD\tVN:1.5\tSO:query\n",
                "@PG\tPN:DIAMOND\tVN:2.1.24\tCL:diamond view -f sam\n",
                "@mm\tBlastP\n",
                "@CO\tBlastP-like alignments\n",
                "@CO\tReporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate\n",
            )
        );
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
        print_match_context(&mut buf, &ctx, &score_matrix, true, true, true).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert_eq!(
            output,
            format!(
                "query\t0\ttarget title<>extra title\t21\t255\t4M\t*\t0\t0\tRNDC\t*\tAS:i:{}\tNM:i:1\tZL:i:200\tZR:i:80\tZE:f:1.00e-12\tZI:i:75\tZF:i:1\tZS:i:2\tMD:Z:4\tZQ:i:6\n",
                score_matrix.bitscore(80.0) as u32
            )
        );
    }

    #[test]
    fn test_print_match_context_title_policy_and_range_validation() {
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
            hsp.clone(),
            0,
            0,
            vec![vec![0, -127, 2, 3, 4, 5]],
            6,
            "query one",
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
        print_match_context(&mut buf, &ctx, &score_matrix, false, false, false).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("query\t0\ttarget\t21\t255\t4M"));
        assert!(output.contains("\trNDC\t"));

        let mut bad = hsp;
        bad.query_range = Interval::new(1, 99);
        let bad_ctx = HspContext::new(
            bad,
            0,
            0,
            vec![vec![0, 1, 2]],
            3,
            "query one",
            0,
            200,
            "target title",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );
        assert_eq!(
            print_match_context(&mut Vec::new(), &bad_ctx, &score_matrix, true, true, false)
                .unwrap_err()
                .kind(),
            io::ErrorKind::InvalidData
        );
    }
}

// String compatibility entry points share the raw renderer.
pub fn print_match_context<W: Write>(
    writer: &mut W,
    r: &HspContext,
    score_matrix: &ScoreMatrix,
    full_titles: bool,
    all_titles: bool,
    sam_qlen_field: bool,
) -> io::Result<()> {
    print_match_context_bytes(
        writer,
        r,
        score_matrix,
        full_titles,
        all_titles,
        sam_qlen_field,
    )
}

pub fn print_query_intro<W: Write>(
    writer: &mut W,
    query_title: &str,
    unaligned: bool,
) -> io::Result<()> {
    print_query_intro_bytes(writer, query_title.as_bytes(), unaligned)
}
