use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::basic::translate::{Strand, TranslatedPosition};
use crate::basic::value::{Letter, AMINO_ACID_ALPHABET, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

/// Data for a pairwise alignment display.
pub struct PairwiseData<'a> {
    pub subject_id: &'a str,
    pub subject_len: i32,
    pub score: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub identities: i32,
    pub positives: i32,
    pub gaps: i32,
    pub length: i32,
    pub query_seq: &'a [Letter],
    pub subject_seq: &'a [Letter],
    pub query_start: i32,
    pub subject_start: i32,
}

/// Write a BLAST-style pairwise alignment.
pub fn write_pairwise<W: Write>(
    writer: &mut W,
    data: &PairwiseData<'_>,
    score_matrix: &ScoreMatrix,
) -> io::Result<()> {
    // Header
    writeln!(writer, ">{}", data.subject_id)?;
    writeln!(writer, "Length={}", data.subject_len)?;
    writeln!(writer)?;

    // Score line — matches C++ blast_pairwise_format.cpp:32.
    write!(
        writer,
        " Score = {} bits ({}),  Expect = ",
        crate::output::format::format_double(data.bit_score),
        data.score
    )?;
    writeln!(
        writer,
        "{}",
        crate::output::format::format_evalue(data.evalue)
    )?;

    // C++ `percentage<unsigned, unsigned>` is integer division — truncates.
    let pct_int = |n: u32, d: i32| -> i32 {
        if d <= 0 {
            0
        } else {
            (n as i32) * 100 / d
        }
    };
    let id_pct = pct_int(data.identities as u32, data.length);
    let pos_pct = pct_int(data.positives as u32, data.length);
    let gap_pct = pct_int(data.gaps as u32, data.length);
    writeln!(
        writer,
        " Identities = {}/{} ({}%), Positives = {}/{} ({}%), Gaps = {}/{} ({}%)",
        data.identities,
        data.length,
        id_pct,
        data.positives,
        data.length,
        pos_pct,
        data.gaps,
        data.length,
        gap_pct
    )?;
    writeln!(writer)?;

    // Alignment blocks (60 columns per line)
    let line_width = 60;
    let aln_len = data.query_seq.len().min(data.subject_seq.len());
    let mut qi = data.query_start;
    let mut si = data.subject_start;

    let mut pos = 0;
    while pos < aln_len {
        let end = (pos + line_width).min(aln_len);
        let block_len = end - pos;

        // Query line
        let q_start = qi + 1;
        write!(writer, "Query  {:>3}  ", q_start)?;
        for &l in &data.query_seq[pos..end] {
            writer.write_all(&[AMINO_ACID_ALPHABET[(l & LETTER_MASK) as usize]])?;
        }
        qi += block_len as i32;
        writeln!(writer, " {}", qi)?;

        // Middle line
        write!(writer, "            ")?;
        for i in pos..end {
            let ql = data.query_seq[i] & LETTER_MASK;
            let sl = data.subject_seq[i] & LETTER_MASK;
            if ql == sl {
                writer.write_all(&[AMINO_ACID_ALPHABET[ql as usize]])?;
            } else if score_matrix.score(ql, sl) > 0 {
                writer.write_all(b"+")?;
            } else {
                writer.write_all(b" ")?;
            }
        }
        writeln!(writer)?;

        // Subject line
        let s_start = si + 1;
        write!(writer, "Sbjct  {:>3}  ", s_start)?;
        for &l in &data.subject_seq[pos..end] {
            writer.write_all(&[AMINO_ACID_ALPHABET[(l & LETTER_MASK) as usize]])?;
        }
        si += block_len as i32;
        writeln!(writer, " {}", si)?;

        writeln!(writer)?;
        pos = end;
    }

    Ok(())
}

/// Matches C++ `PairwiseFormat::print_match(const HspContext&)`.
pub fn print_match_context<W: Write>(
    writer: &mut W,
    r: &HspContext,
    score_matrix: &ScoreMatrix,
    query_translated: bool,
    blastn_command: bool,
) -> io::Result<()> {
    const WIDTH: usize = 60;
    let dna_len = r.query_source_len;
    let strand = if blastn_command {
        if r.frame() != 0 {
            Strand::Forward
        } else {
            Strand::Reverse
        }
    } else if r.frame() < 3 {
        Strand::Forward
    } else {
        Strand::Reverse
    };
    // C++ `OutputFormat::print_title(out, target_title, true, true, " ")`
    // (`diamond/src/output/blast_pairwise_format.cpp:30`) tokenizes on BOTH
    // `\x01` AND the FASTA `" >"` header separator and joins surviving
    // titles with a space. A naive `replace('\x01', ' ')` only handles the
    // `\x01` case, so multi-seqid headers with `" >"` separators would
    // print verbatim instead of joined-by-space.
    let target_title = crate::output::format::print_title(&r.target_title, true, true, " ", false);
    writeln!(writer, ">{}", target_title)?;
    writeln!(writer, "Length={}", r.subject_len)?;
    writeln!(writer)?;
    write!(
        writer,
        " Score = {} bits ({}),  Expect = ",
        crate::output::format::format_double(r.bit_score()),
        r.score()
    )?;
    writeln!(
        writer,
        "{}",
        crate::output::format::format_evalue(r.evalue())
    )?;

    // C++ uses `percentage<unsigned, unsigned>(x, y) = x * 100 / y` (integer
    // division — truncates), not rounded. See diamond/src/util/util.h.
    let pct = |n: u32, d: u32| -> u32 {
        if d == 0 {
            0
        } else {
            n * 100 / d
        }
    };
    writeln!(
        writer,
        " Identities = {}/{} ({}%), Positives = {}/{} ({}%), Gaps = {}/{} ({}%)",
        r.identities(),
        r.length(),
        pct(r.identities(), r.length()),
        r.positives(),
        r.length(),
        pct(r.positives(), r.length()),
        r.gaps(),
        r.length(),
        pct(r.gaps(), r.length()),
    )?;
    if query_translated {
        writeln!(writer, " Frame = {}", r.blast_query_frame(query_translated))?;
    }
    if blastn_command {
        write!(writer, " Strand = ")?;
        if r.frame() != 0 {
            write!(writer, "Plus/")?;
        } else {
            write!(writer, "Minus/")?;
        }
        writeln!(writer, "Plus")?;
    }
    writeln!(writer)?;

    // Matches C++ blast_pairwise_format.cpp:49 exactly:
    //   digits = max(ceil(log10(subject_end)), ceil(log10(query_end)))
    // Note this gives 3 for end=1000 (log10 == 3.0), not 4 — different from
    // x.to_string().len() at exact powers of 10.
    let s_end = r.subject_range().end.max(1) as f64;
    let q_end = r.query_source_range().end.max(1) as f64;
    let digits = s_end.log10().ceil().max(q_end.log10().ceil()) as usize;
    let mut qi = r.begin();
    let mut mi = r.begin();
    let mut si = r.begin();
    while qi.good() {
        write!(
            writer,
            "Query  {:>width$}  ",
            qi.query_pos().absolute(dna_len, query_translated, false) + 1,
            width = digits
        )?;
        for _ in 0..WIDTH {
            if !qi.good() {
                break;
            }
            write!(writer, "{}", qi.query_char())?;
            qi.advance();
        }
        write!(
            writer,
            " {}",
            TranslatedPosition::oriented_position(
                qi.query_pos().in_strand(query_translated) - 1,
                strand,
                dna_len,
            ) + 1
        )?;
        writeln!(writer)?;

        for _ in 0..digits + 9 {
            write!(writer, " ")?;
        }
        for _ in 0..WIDTH {
            if !mi.good() {
                break;
            }
            write!(
                writer,
                "{}",
                mi.midline_char(score_matrix.score(mi.query(), mi.subject()))
            )?;
            mi.advance();
        }
        writeln!(writer)?;

        write!(
            writer,
            "Sbjct  {:>width$}  ",
            si.subject_pos() + 1,
            width = digits
        )?;
        for _ in 0..WIDTH {
            if !si.good() {
                break;
            }
            write!(writer, "{}", si.subject_char())?;
            si.advance();
        }
        writeln!(writer, " {}", si.subject_pos())?;
        writeln!(writer)?;
    }

    Ok(())
}

/// Write the BLAST pairwise header.
pub fn write_header<W: Write>(writer: &mut W) -> io::Result<()> {
    writeln!(writer, "BLASTP 2.3.0+")?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

/// Write the query header.
pub fn write_query_header<W: Write>(
    writer: &mut W,
    query_id: &str,
    query_len: i32,
) -> io::Result<()> {
    writeln!(writer)?;
    writeln!(writer, "Query= {}", query_id)?;
    writeln!(writer)?;
    writeln!(writer, "Length={}", query_len)?;
    writeln!(writer)?;
    Ok(())
}

/// Matches C++ `PairwiseFormat::print_query_intro()`.
pub fn print_query_intro<W: Write>(
    writer: &mut W,
    query_title: &str,
    query_len: i32,
    unaligned: bool,
) -> io::Result<()> {
    writeln!(writer, "Query= {}", query_title)?;
    writeln!(writer)?;
    writeln!(writer, "Length={}", query_len)?;
    writeln!(writer)?;
    if unaligned {
        writeln!(writer)?;
        writeln!(writer, "***** No hits found *****")?;
        writeln!(writer)?;
        writeln!(writer)?;
    }
    Ok(())
}

/// Matches C++ `PairwiseFormat::print_query_epilog()`.
pub fn print_query_epilog<W: Write>(_writer: &mut W) -> io::Result<()> {
    Ok(())
}

/// Matches C++ `PairwiseFormat::print_footer()`.
pub fn print_footer<W: Write>(_writer: &mut W) -> io::Result<()> {
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() {
        let mut buf = Vec::new();
        write_header(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("BLASTP 2.3.0+"));
    }

    #[test]
    fn test_write_pairwise() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0i8, 1, 2, 3];
        let subject = vec![0i8, 1, 2, 3];
        let data = PairwiseData {
            subject_id: "s1",
            subject_len: 4,
            score: 15,
            bit_score: 30.0,
            evalue: 1e-5,
            identities: 4,
            positives: 4,
            gaps: 0,
            length: 4,
            query_seq: &query,
            subject_seq: &subject,
            query_start: 0,
            subject_start: 0,
        };
        let mut buf = Vec::new();
        write_pairwise(&mut buf, &data, &sm).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains(">s1"));
        assert!(output.contains("Score = 30.0 bits"));
        assert!(output.contains("ARND"));
    }

    #[test]
    fn test_print_match_context() {
        use crate::align::hsp::{Hsp, HspContext};
        use crate::basic::packed_transcript::EditOperation;
        use crate::util::interval::Interval;

        let sm = ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1000).unwrap();
        let mut hsp = Hsp::new();
        hsp.score = 42;
        hsp.bit_score = 21.0;
        hsp.evalue = 1e-6;
        hsp.frame = 0;
        hsp.identities = 1;
        hsp.positives = 2;
        hsp.gaps = 1;
        hsp.length = 3;
        hsp.query_range = Interval::new(0, 3);
        hsp.query_source_range = Interval::new(0, 3);
        hsp.subject_range = Interval::new(5, 8);
        hsp.transcript.push_with_count(EditOperation::Match, 1);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 1);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 3);

        let ctx = HspContext::new(
            hsp,
            0,
            0,
            vec![vec![0, 1, 4]],
            3,
            "query",
            0,
            100,
            "subject title",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );

        let mut buf = Vec::new();
        print_match_context(&mut buf, &ctx, &sm, false, false).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains(">subject title"));
        assert!(output.contains("Length=100"));
        assert!(output.contains("Score = 21.0 bits (42),  Expect = 1.00e-06"));
        // C++ truncates: 2*100/3 = 66 (not rounded to 67).
        assert!(output.contains("Identities = 1/3 (33%), Positives = 2/3 (66%), Gaps = 1/3 (33%)"));
        assert!(output.contains("Query  1  AR- 2"));
        assert!(output.contains("A+ "));
        assert!(output.contains("Sbjct  6  ARD 8"));
    }

    #[test]
    fn test_print_query_intro_epilog_footer() {
        let mut buf = Vec::new();
        print_query_intro(&mut buf, "query one", 12, true).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("Query= query one\n\nLength=12\n\n"));
        assert!(output.contains("***** No hits found *****"));

        let mut buf = Vec::new();
        print_query_epilog(&mut buf).unwrap();
        print_footer(&mut buf).unwrap();
        assert!(buf.is_empty());
    }
}
