use std::io::{self, Write};

use crate::basic::value::{Letter, AMINO_ACID_ALPHABET, LETTER_MASK};
use super::format::format_evalue;

/// Data for an XML alignment hit.
pub struct XmlHit<'a> {
    pub hit_num: i32,
    pub subject_id: &'a str,
    pub subject_def: &'a str,
    pub subject_len: i32,
    pub hsps: Vec<XmlHsp<'a>>,
}

/// Data for a single HSP in XML format.
pub struct XmlHsp<'a> {
    pub num: i32,
    pub bit_score: f64,
    pub score: i32,
    pub evalue: f64,
    pub query_from: i32,
    pub query_to: i32,
    pub hit_from: i32,
    pub hit_to: i32,
    pub query_frame: i32,
    pub identity: i32,
    pub positive: i32,
    pub gaps: i32,
    pub align_len: i32,
    pub qseq: &'a [Letter],
    pub hseq: &'a [Letter],
}

/// Write BLAST XML header.
pub fn write_xml_header<W: Write>(
    writer: &mut W,
    program: &str,
    database: &str,
) -> io::Result<()> {
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">")?;
    writeln!(writer, "<BlastOutput>")?;
    writeln!(writer, "  <BlastOutput_program>{}</BlastOutput_program>", program)?;
    writeln!(writer, "  <BlastOutput_version>{} 2.3.0+</BlastOutput_version>", program)?;
    writeln!(writer, "  <BlastOutput_db>{}</BlastOutput_db>", database)?;
    writeln!(writer, "  <BlastOutput_iterations>")?;
    Ok(())
}

/// Write XML footer.
pub fn write_xml_footer<W: Write>(writer: &mut W) -> io::Result<()> {
    writeln!(writer, "  </BlastOutput_iterations>")?;
    writeln!(writer, "</BlastOutput>")?;
    Ok(())
}

/// Write a query iteration.
pub fn write_iteration_start<W: Write>(
    writer: &mut W,
    iter_num: i32,
    query_id: &str,
    query_len: i32,
) -> io::Result<()> {
    writeln!(writer, "    <Iteration>")?;
    writeln!(writer, "      <Iteration_iter-num>{}</Iteration_iter-num>", iter_num)?;
    writeln!(writer, "      <Iteration_query-def>{}</Iteration_query-def>", escape_xml(query_id))?;
    writeln!(writer, "      <Iteration_query-len>{}</Iteration_query-len>", query_len)?;
    writeln!(writer, "      <Iteration_hits>")?;
    Ok(())
}

/// Write iteration end.
pub fn write_iteration_end<W: Write>(writer: &mut W) -> io::Result<()> {
    writeln!(writer, "      </Iteration_hits>")?;
    writeln!(writer, "    </Iteration>")?;
    Ok(())
}

/// Write a single hit.
pub fn write_hit<W: Write>(writer: &mut W, hit: &XmlHit<'_>) -> io::Result<()> {
    writeln!(writer, "        <Hit>")?;
    writeln!(writer, "          <Hit_num>{}</Hit_num>", hit.hit_num)?;
    writeln!(writer, "          <Hit_id>{}</Hit_id>", escape_xml(hit.subject_id))?;
    writeln!(writer, "          <Hit_def>{}</Hit_def>", escape_xml(hit.subject_def))?;
    writeln!(writer, "          <Hit_len>{}</Hit_len>", hit.subject_len)?;
    writeln!(writer, "          <Hit_hsps>")?;

    for hsp in &hit.hsps {
        write_hsp(writer, hsp)?;
    }

    writeln!(writer, "          </Hit_hsps>")?;
    writeln!(writer, "        </Hit>")?;
    Ok(())
}

fn write_hsp<W: Write>(writer: &mut W, hsp: &XmlHsp<'_>) -> io::Result<()> {
    let qseq_str = seq_to_string(hsp.qseq);
    let hseq_str = seq_to_string(hsp.hseq);

    // Build midline
    let midline: String = hsp
        .qseq
        .iter()
        .zip(hsp.hseq.iter())
        .map(|(&q, &h)| {
            let q = q & LETTER_MASK;
            let h = h & LETTER_MASK;
            if q == h {
                AMINO_ACID_ALPHABET[q as usize] as char
            } else {
                ' '
            }
        })
        .collect();

    writeln!(writer, "            <Hsp>")?;
    writeln!(writer, "              <Hsp_num>{}</Hsp_num>", hsp.num)?;
    writeln!(writer, "              <Hsp_bit-score>{}</Hsp_bit-score>", hsp.bit_score)?;
    writeln!(writer, "              <Hsp_score>{}</Hsp_score>", hsp.score)?;
    writeln!(writer, "              <Hsp_evalue>{}</Hsp_evalue>", format_evalue(hsp.evalue))?;
    writeln!(writer, "              <Hsp_query-from>{}</Hsp_query-from>", hsp.query_from)?;
    writeln!(writer, "              <Hsp_query-to>{}</Hsp_query-to>", hsp.query_to)?;
    writeln!(writer, "              <Hsp_hit-from>{}</Hsp_hit-from>", hsp.hit_from)?;
    writeln!(writer, "              <Hsp_hit-to>{}</Hsp_hit-to>", hsp.hit_to)?;
    writeln!(writer, "              <Hsp_query-frame>{}</Hsp_query-frame>", hsp.query_frame)?;
    writeln!(writer, "              <Hsp_identity>{}</Hsp_identity>", hsp.identity)?;
    writeln!(writer, "              <Hsp_positive>{}</Hsp_positive>", hsp.positive)?;
    writeln!(writer, "              <Hsp_gaps>{}</Hsp_gaps>", hsp.gaps)?;
    writeln!(writer, "              <Hsp_align-len>{}</Hsp_align-len>", hsp.align_len)?;
    writeln!(writer, "              <Hsp_qseq>{}</Hsp_qseq>", qseq_str)?;
    writeln!(writer, "              <Hsp_hseq>{}</Hsp_hseq>", hseq_str)?;
    writeln!(writer, "              <Hsp_midline>{}</Hsp_midline>", midline)?;
    writeln!(writer, "            </Hsp>")?;
    Ok(())
}

fn seq_to_string(seq: &[Letter]) -> String {
    seq.iter()
        .map(|&l| AMINO_ACID_ALPHABET[(l & LETTER_MASK) as usize] as char)
        .collect()
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_xml_header() {
        let mut buf = Vec::new();
        write_xml_header(&mut buf, "blastp", "test.dmnd").unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("<?xml version"));
        assert!(output.contains("<BlastOutput>"));
        assert!(output.contains("blastp"));
    }

    #[test]
    fn test_escape_xml() {
        assert_eq!(escape_xml("a<b>c&d"), "a&lt;b&gt;c&amp;d");
    }

    #[test]
    fn test_write_hit() {
        let qseq = vec![0i8, 1, 2];
        let hseq = vec![0i8, 1, 3];
        let hit = XmlHit {
            hit_num: 1,
            subject_id: "s1",
            subject_def: "test subject",
            subject_len: 100,
            hsps: vec![XmlHsp {
                num: 1,
                bit_score: 50.0,
                score: 30,
                evalue: 1e-5,
                query_from: 1,
                query_to: 3,
                hit_from: 1,
                hit_to: 3,
                query_frame: 0,
                identity: 2,
                positive: 2,
                gaps: 0,
                align_len: 3,
                qseq: &qseq,
                hseq: &hseq,
            }],
        };
        let mut buf = Vec::new();
        write_hit(&mut buf, &hit).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("<Hit_id>s1</Hit_id>"));
        assert!(output.contains("<Hsp_score>30</Hsp_score>"));
    }
}
