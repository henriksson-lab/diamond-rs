use std::io::{self, Write};

use super::format::{format_double, format_evalue, print_title_xml};
use crate::align::hsp::HspContext;
use crate::basic::consts::VERSION_STRING;
use crate::basic::value::{Letter, AMINO_ACID_ALPHABET, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence::{get_accession, get_title_def, AccessionParsing};

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
pub fn write_xml_header<W: Write>(writer: &mut W, program: &str, database: &str) -> io::Result<()> {
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">")?;
    writeln!(writer, "<BlastOutput>")?;
    writeln!(
        writer,
        "  <BlastOutput_program>{}</BlastOutput_program>",
        program
    )?;
    writeln!(
        writer,
        "  <BlastOutput_version>{} 2.3.0+</BlastOutput_version>",
        program
    )?;
    writeln!(writer, "  <BlastOutput_db>{}</BlastOutput_db>", database)?;
    writeln!(writer, "  <BlastOutput_iterations>")?;
    Ok(())
}

/// Matches C++ `XMLFormat::print_header()`.
pub fn print_header<W: Write>(
    writer: &mut W,
    mode: u32,
    matrix: &str,
    gap_open: i32,
    gap_extend: i32,
    evalue: f64,
    first_query_name: &str,
    first_query_len: u32,
    database: &str,
) -> io::Result<()> {
    let mode_str = match mode {
        2 => "blastp",
        3 => "blastx",
        4 => "blastn",
        _ => "",
    };
    let first_query_cut = first_query_name
        .find('\x01')
        .unwrap_or(first_query_name.len());
    let escaped_first_query_name = escape_xml(first_query_name);
    let first_query_def =
        &escaped_first_query_name[..first_query_cut.min(escaped_first_query_name.len())];
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(
        writer,
        "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">"
    )?;
    writeln!(writer, "<BlastOutput>")?;
    writeln!(
        writer,
        "  <BlastOutput_program>{}</BlastOutput_program>",
        mode_str
    )?;
    writeln!(
        writer,
        "  <BlastOutput_version>diamond {}</BlastOutput_version>",
        VERSION_STRING
    )?;
    writeln!(
        writer,
        "  <BlastOutput_reference>Benjamin Buchfink, Xie Chao, and Daniel Huson (2015), &quot;Fast and sensitive protein alignment using DIAMOND&quot;, Nature Methods 12:59-60.</BlastOutput_reference>"
    )?;
    writeln!(writer, "  <BlastOutput_db>{}</BlastOutput_db>", database)?;
    writeln!(
        writer,
        "  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>"
    )?;
    writeln!(
        writer,
        "  <BlastOutput_query-def>{}</BlastOutput_query-def>",
        first_query_def
    )?;
    writeln!(
        writer,
        "  <BlastOutput_query-len>{}</BlastOutput_query-len>",
        first_query_len
    )?;
    writeln!(writer, "  <BlastOutput_param>")?;
    writeln!(writer, "    <Parameters>")?;
    writeln!(
        writer,
        "      <Parameters_matrix>{}</Parameters_matrix>",
        matrix
    )?;
    writeln!(
        writer,
        "      <Parameters_expect>{}</Parameters_expect>",
        evalue
    )?;
    writeln!(
        writer,
        "      <Parameters_gap-open>{}</Parameters_gap-open>",
        gap_open
    )?;
    writeln!(
        writer,
        "      <Parameters_gap-extend>{}</Parameters_gap-extend>",
        gap_extend
    )?;
    writeln!(writer, "      <Parameters_filter>F</Parameters_filter>")?;
    writeln!(writer, "    </Parameters>")?;
    writeln!(writer, "  </BlastOutput_param>")?;
    writeln!(writer, "<BlastOutput_iterations>")?;
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
    writeln!(
        writer,
        "      <Iteration_iter-num>{}</Iteration_iter-num>",
        iter_num
    )?;
    writeln!(
        writer,
        "      <Iteration_query-def>{}</Iteration_query-def>",
        escape_xml(query_id)
    )?;
    writeln!(
        writer,
        "      <Iteration_query-len>{}</Iteration_query-len>",
        query_len
    )?;
    writeln!(writer, "      <Iteration_hits>")?;
    Ok(())
}

/// Write iteration end.
pub fn write_iteration_end<W: Write>(writer: &mut W) -> io::Result<()> {
    writeln!(writer, "      </Iteration_hits>")?;
    writeln!(writer, "    </Iteration>")?;
    Ok(())
}

/// Matches C++ `XMLFormat::print_query_intro()`.
pub fn print_query_intro<W: Write>(
    writer: &mut W,
    query_oid: u64,
    query_title: &str,
    query_len: i32,
) -> io::Result<()> {
    writeln!(writer, "<Iteration>")?;
    writeln!(
        writer,
        "  <Iteration_iter-num>{}</Iteration_iter-num>",
        query_oid + 1
    )?;
    writeln!(
        writer,
        "  <Iteration_query-ID>Query_{}</Iteration_query-ID>",
        query_oid + 1
    )?;
    writeln!(
        writer,
        "  <Iteration_query-def>{}</Iteration_query-def>",
        print_title_xml(query_title, true, false, "", false)
    )?;
    writeln!(
        writer,
        "  <Iteration_query-len>{}</Iteration_query-len>",
        query_len
    )?;
    writeln!(writer, "<Iteration_hits>")?;
    Ok(())
}

/// Matches C++ `XMLFormat::print_query_epilog()`.
pub fn print_query_epilog<W: Write>(
    writer: &mut W,
    unaligned: bool,
    db_seqs: u64,
    db_letters: u64,
    kappa: f64,
    lambda: f64,
) -> io::Result<()> {
    if !unaligned {
        writeln!(writer, "  </Hit_hsps>")?;
        writeln!(writer, "</Hit>")?;
    }
    writeln!(writer, "</Iteration_hits>")?;
    writeln!(writer, "  <Iteration_stat>")?;
    writeln!(writer, "    <Statistics>")?;
    writeln!(
        writer,
        "      <Statistics_db-num>{}</Statistics_db-num>",
        db_seqs
    )?;
    writeln!(
        writer,
        "      <Statistics_db-len>{}</Statistics_db-len>",
        db_letters
    )?;
    writeln!(writer, "      <Statistics_hsp-len>0</Statistics_hsp-len>")?;
    writeln!(
        writer,
        "      <Statistics_eff-space>0</Statistics_eff-space>"
    )?;
    // C++ uses iostream default (6 significant digits via %f).
    writeln!(
        writer,
        "      <Statistics_kappa>{:.6}</Statistics_kappa>",
        kappa
    )?;
    writeln!(
        writer,
        "      <Statistics_lambda>{:.6}</Statistics_lambda>",
        lambda
    )?;
    writeln!(writer, "      <Statistics_entropy>0</Statistics_entropy>")?;
    writeln!(writer, "    </Statistics>")?;
    writeln!(writer, "  </Iteration_stat>")?;
    writeln!(writer, "</Iteration>")?;
    Ok(())
}

/// Matches C++ `XMLFormat::print_footer()`.
pub fn print_footer<W: Write>(writer: &mut W) -> io::Result<()> {
    writeln!(writer, "</BlastOutput_iterations>")?;
    write!(writer, "</BlastOutput>")?;
    Ok(())
}

/// Write a single hit.
pub fn write_hit<W: Write>(writer: &mut W, hit: &XmlHit<'_>) -> io::Result<()> {
    writeln!(writer, "        <Hit>")?;
    writeln!(writer, "          <Hit_num>{}</Hit_num>", hit.hit_num)?;
    writeln!(
        writer,
        "          <Hit_id>{}</Hit_id>",
        escape_xml(hit.subject_id)
    )?;
    writeln!(
        writer,
        "          <Hit_def>{}</Hit_def>",
        escape_xml(hit.subject_def)
    )?;
    writeln!(writer, "          <Hit_len>{}</Hit_len>", hit.subject_len)?;
    writeln!(writer, "          <Hit_hsps>")?;

    for hsp in &hit.hsps {
        write_hsp(writer, hsp)?;
    }

    writeln!(writer, "          </Hit_hsps>")?;
    writeln!(writer, "        </Hit>")?;
    Ok(())
}

/// Matches C++ `XMLFormat::print_match(const HspContext&)`.
pub fn print_match_context<W: Write>(
    writer: &mut W,
    r: &HspContext,
    score_matrix: &ScoreMatrix,
    query_translated: bool,
    xml_blord_format: bool,
    no_parse_seqids: bool,
) -> io::Result<()> {
    if r.hsp_num == 0 {
        if r.hit_num > 0 {
            writeln!(writer, "  </Hit_hsps>")?;
            writeln!(writer, "</Hit>")?;
        }

        writeln!(writer, "<Hit>")?;
        writeln!(writer, "  <Hit_num>{}</Hit_num>", r.hit_num + 1)?;
        let (id, def) = get_title_def(&r.target_title);
        if xml_blord_format {
            writeln!(writer, "  <Hit_id>gnl|BL_ORD_ID|{}</Hit_id>", r.subject_oid)?;
            writeln!(
                writer,
                "  <Hit_def>{}</Hit_def>",
                print_title_xml(&r.target_title, true, true, " &gt;", false)
            )?;
        } else {
            writeln!(writer, "  <Hit_id>{}</Hit_id>", escape_xml(&id))?;
            writeln!(
                writer,
                "  <Hit_def>{}</Hit_def>",
                print_title_xml(&def, true, true, " &gt;", false)
            )?;
        }
        let accession = if no_parse_seqids {
            id
        } else {
            get_accession(&id, &mut AccessionParsing::default())
        };
        writeln!(
            writer,
            "  <Hit_accession>{}</Hit_accession>",
            escape_xml(&accession)
        )?;
        writeln!(writer, "  <Hit_len>{}</Hit_len>", r.subject_len)?;
        writeln!(writer, "  <Hit_hsps>")?;
    }

    writeln!(writer, "    <Hsp>")?;
    writeln!(writer, "      <Hsp_num>{}</Hsp_num>", r.hsp_num + 1)?;
    writeln!(
        writer,
        "      <Hsp_bit-score>{}</Hsp_bit-score>",
        format_double(r.bit_score())
    )?;
    writeln!(writer, "      <Hsp_score>{}</Hsp_score>", r.score())?;
    writeln!(
        writer,
        "      <Hsp_evalue>{}</Hsp_evalue>",
        format_evalue(r.evalue())
    )?;
    writeln!(
        writer,
        "      <Hsp_query-from>{}</Hsp_query-from>",
        r.query_source_range().begin + 1
    )?;
    writeln!(
        writer,
        "      <Hsp_query-to>{}</Hsp_query-to>",
        r.query_source_range().end
    )?;
    writeln!(
        writer,
        "      <Hsp_hit-from>{}</Hsp_hit-from>",
        r.subject_range().begin + 1
    )?;
    writeln!(
        writer,
        "      <Hsp_hit-to>{}</Hsp_hit-to>",
        r.subject_range().end
    )?;
    writeln!(
        writer,
        "      <Hsp_query-frame>{}</Hsp_query-frame>",
        r.blast_query_frame(query_translated)
    )?;
    writeln!(writer, "      <Hsp_hit-frame>0</Hsp_hit-frame>")?;
    writeln!(
        writer,
        "      <Hsp_identity>{}</Hsp_identity>",
        r.identities()
    )?;
    writeln!(
        writer,
        "      <Hsp_positive>{}</Hsp_positive>",
        r.positives()
    )?;
    writeln!(writer, "      <Hsp_gaps>{}</Hsp_gaps>", r.gaps())?;
    writeln!(
        writer,
        "      <Hsp_align-len>{}</Hsp_align-len>",
        r.length()
    )?;
    write!(writer, "         <Hsp_qseq>")?;
    let mut it = r.begin();
    while it.good() {
        write!(writer, "{}", it.query_char())?;
        it.advance();
    }
    writeln!(writer, "</Hsp_qseq>")?;

    write!(writer, "         <Hsp_hseq>")?;
    let mut it = r.begin();
    while it.good() {
        write!(writer, "{}", it.subject_char())?;
        it.advance();
    }
    writeln!(writer, "</Hsp_hseq>")?;

    write!(writer, "      <Hsp_midline>")?;
    let mut it = r.begin();
    while it.good() {
        write!(
            writer,
            "{}",
            it.midline_char(score_matrix.score(it.query(), it.subject()))
        )?;
        it.advance();
    }
    writeln!(writer, "</Hsp_midline>")?;
    writeln!(writer, "    </Hsp>")?;
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
    writeln!(
        writer,
        "              <Hsp_bit-score>{}</Hsp_bit-score>",
        hsp.bit_score
    )?;
    writeln!(writer, "              <Hsp_score>{}</Hsp_score>", hsp.score)?;
    writeln!(
        writer,
        "              <Hsp_evalue>{}</Hsp_evalue>",
        format_evalue(hsp.evalue)
    )?;
    writeln!(
        writer,
        "              <Hsp_query-from>{}</Hsp_query-from>",
        hsp.query_from
    )?;
    writeln!(
        writer,
        "              <Hsp_query-to>{}</Hsp_query-to>",
        hsp.query_to
    )?;
    writeln!(
        writer,
        "              <Hsp_hit-from>{}</Hsp_hit-from>",
        hsp.hit_from
    )?;
    writeln!(
        writer,
        "              <Hsp_hit-to>{}</Hsp_hit-to>",
        hsp.hit_to
    )?;
    writeln!(
        writer,
        "              <Hsp_query-frame>{}</Hsp_query-frame>",
        hsp.query_frame
    )?;
    writeln!(
        writer,
        "              <Hsp_identity>{}</Hsp_identity>",
        hsp.identity
    )?;
    writeln!(
        writer,
        "              <Hsp_positive>{}</Hsp_positive>",
        hsp.positive
    )?;
    writeln!(writer, "              <Hsp_gaps>{}</Hsp_gaps>", hsp.gaps)?;
    writeln!(
        writer,
        "              <Hsp_align-len>{}</Hsp_align-len>",
        hsp.align_len
    )?;
    writeln!(writer, "              <Hsp_qseq>{}</Hsp_qseq>", qseq_str)?;
    writeln!(writer, "              <Hsp_hseq>{}</Hsp_hseq>", hseq_str)?;
    writeln!(
        writer,
        "              <Hsp_midline>{}</Hsp_midline>",
        midline
    )?;
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
    fn test_print_header() {
        let mut buf = Vec::new();
        print_header(
            &mut buf,
            2,
            "BLOSUM62",
            11,
            1,
            0.001,
            "query <one>\x01second",
            42,
            "db<&>.dmnd",
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("<BlastOutput_program>blastp</BlastOutput_program>"));
        assert!(output.contains("<BlastOutput_version>diamond 2.1.24</BlastOutput_version>"));
        assert!(output.contains("<BlastOutput_db>db<&>.dmnd</BlastOutput_db>"));
        assert!(
            output.contains("<BlastOutput_query-def>query &lt;o</BlastOutput_query-def>"),
            "{output}"
        );
        assert!(output.contains("<Parameters_gap-open>11</Parameters_gap-open>"));
        assert!(output.ends_with("<BlastOutput_iterations>\n"));
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

    #[test]
    fn test_print_match_context() {
        use crate::align::hsp::{Hsp, HspContext};
        use crate::basic::packed_transcript::EditOperation;
        use crate::util::interval::Interval;

        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1000).unwrap();
        let mut hsp = Hsp::new();
        hsp.score = 50;
        hsp.bit_score = 25.0;
        hsp.evalue = 1e-8;
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
            7,
            100,
            "sp|P0|target definition >sp|P1|other other",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );

        let mut buf = Vec::new();
        print_match_context(&mut buf, &ctx, &score_matrix, false, false, false).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("<Hit_num>1</Hit_num>"));
        assert!(output.contains("<Hit_id>sp|P0|target</Hit_id>"));
        assert!(output.contains("<Hit_def>definition &gt;sp|P1|other other</Hit_def>"));
        assert!(output.contains("<Hit_accession>P0</Hit_accession>"));
        assert!(output.contains("<Hsp_score>50</Hsp_score>"));
        assert!(output.contains("<Hsp_qseq>AR-</Hsp_qseq>"));
        assert!(output.contains("<Hsp_hseq>ARD</Hsp_hseq>"));
        assert!(output.contains("<Hsp_midline>A+ </Hsp_midline>"));
    }

    #[test]
    fn test_print_query_intro_epilog_footer() {
        let mut buf = Vec::new();
        print_query_intro(&mut buf, 2, "query <one> >second", 33).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("<Iteration_iter-num>3</Iteration_iter-num>"));
        assert!(output.contains("<Iteration_query-ID>Query_3</Iteration_query-ID>"));
        assert!(output.contains("<Iteration_query-def>query &lt;one&gt;</Iteration_query-def>"));
        assert!(output.contains("<Iteration_query-len>33</Iteration_query-len>"));

        let mut buf = Vec::new();
        print_query_epilog(&mut buf, false, 7, 99, 0.1, 0.2).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("  </Hit_hsps>\n</Hit>\n</Iteration_hits>"));
        assert!(output.contains("<Statistics_db-num>7</Statistics_db-num>"));
        assert!(output.contains("<Statistics_lambda>0.200000</Statistics_lambda>"));

        let mut buf = Vec::new();
        print_footer(&mut buf).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "</BlastOutput_iterations>\n</BlastOutput>"
        );
    }
}
