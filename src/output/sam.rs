use std::io::{self, Write};

use crate::basic::value::{Letter, AMINO_ACID_ALPHABET, LETTER_MASK};

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
}
