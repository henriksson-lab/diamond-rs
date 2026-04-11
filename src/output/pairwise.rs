use std::io::{self, Write};

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

    // Score line
    write!(
        writer,
        " Score = {} bits ({}),  Expect = ",
        data.bit_score.floor() as i64,
        data.score
    )?;
    if data.evalue == 0.0 {
        writeln!(writer, "0.0")?;
    } else {
        writeln!(writer, "{:.2e}", data.evalue)?;
    }

    // Statistics line
    let id_pct = if data.length > 0 {
        (100.0 * data.identities as f64 / data.length as f64).round() as i32
    } else {
        0
    };
    let pos_pct = if data.length > 0 {
        (100.0 * data.positives as f64 / data.length as f64).round() as i32
    } else {
        0
    };
    let gap_pct = if data.length > 0 {
        (100.0 * data.gaps as f64 / data.length as f64).round() as i32
    } else {
        0
    };
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

/// Write the BLAST pairwise header.
pub fn write_header<W: Write>(writer: &mut W) -> io::Result<()> {
    writeln!(writer, "BLASTP 2.3.0+")?;
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
        assert!(output.contains("Score = 30 bits"));
        assert!(output.contains("ARND"));
    }
}
