use std::io::{self, Write};

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
        r.query_id, r.query_len, r.query_start, r.query_end,
        r.subject_id, r.subject_len, r.subject_start, r.subject_end,
        r.identities, r.length, r.score, format_evalue(r.evalue),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_paf() {
        let r = PafRecord {
            query_id: "q1", query_len: 100, query_start: 0, query_end: 50,
            subject_id: "s1", subject_len: 200, subject_start: 10, subject_end: 60,
            identities: 40, length: 50, score: 100, evalue: 1e-10,
        };
        let mut buf = Vec::new();
        write_paf_record(&mut buf, &r).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("q1\t100\t0\t50\t+\ts1\t200\t10\t60\t40\t50\t255"));
        assert!(output.contains("AS:i:100"));
    }
}
