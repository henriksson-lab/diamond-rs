//! Byte-level Source output contracts. Numeric rendering remains covered by
//! the existing ASCII/UTF-8 formatter suites; these expose canonical byte flow.
use super::format::{self, FieldId};
use crate::align::hsp::{ByteHspContext, Hsp};
use crate::basic::packed_transcript::EditOperation;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::interval::Interval;

fn context() -> ByteHspContext {
    let mut hsp = Hsp::new();
    hsp.score = 30;
    hsp.bit_score = 10.0;
    hsp.evalue = 0.1;
    hsp.length = 1;
    hsp.identities = 1;
    hsp.positives = 1;
    hsp.query_range = Interval::new(0, 1);
    hsp.query_source_range = Interval::new(0, 1);
    hsp.subject_range = Interval::new(0, 1);
    hsp.subject_source_range = Interval::new(0, 1);
    hsp.transcript.push_with_count(EditOperation::Match, 1);
    ByteHspContext::new_bytes(
        hsp,
        0,
        0,
        vec![vec![0]],
        1,
        b"q\xff\0tail full".to_vec(),
        0,
        1,
        b"s\x80 desc\x01other\0drop".to_vec(),
        0,
        0,
        vec![0],
        1.0,
        1.0,
    )
}
#[test]
fn title_policies_escape_raw_bytes_and_stop_only_at_cstring_boundaries() {
    assert_eq!(
        format::print_title_bytes(b"s\xff desc\x01other\0ignored", true, true, b"<>", false),
        b"s\xff desc<>other"
    );
    assert_eq!(
        format::print_title_xml_bytes(b"x\xff<&", true, false, b"", false),
        b"x\xff&lt;&amp;"
    );
    assert_eq!(
        format::print_title_xml("caf\u{e9}<&", true, false, "", false).as_bytes(),
        b"caf\xc3\xa9&lt;&amp;"
    );
}
#[test]
fn tabular_and_json_query_titles_follow_source_c_string_extent() {
    let r = context();
    let mut out = Vec::new();
    format::write_tabular_context_row_bytes(
        &mut out,
        &r,
        &[
            FieldId::QSeqId,
            FieldId::QTitle,
            FieldId::SSeqId,
            FieldId::SAllTitles,
        ],
        false,
        false,
        0,
        0,
    )
    .unwrap();
    assert_eq!(out, b"q\xff\tq\xff\ts\x80\ts\x80 desc<>other\n");
    out.clear();
    format::write_tabular_context_row_json_bytes(
        &mut out,
        &r,
        false,
        &[FieldId::QTitle],
        false,
        false,
        0,
        0,
        &[],
        None,
        b"",
        None,
    )
    .unwrap();
    assert_eq!(out, b"\n\t{\n\t\"qtitle\":\"q\xff\"\n\t}");
}
#[test]
fn quality_slicing_is_byte_based_and_applies_source_cstring_output() {
    let r = context();
    let mut out = Vec::new();
    format::write_tabular_context_row_with_query_info_bytes(
        &mut out,
        &r,
        &[FieldId::QQual, FieldId::FullQQual],
        false,
        false,
        0,
        0,
        b"\xff\0tail",
        None,
    )
    .unwrap();
    assert_eq!(out, b"\xff\t\xff\n");
    assert_eq!(
        format::print_qqual_bytes(b"\0nonempty", Interval::new(1, 3)),
        b"*"
    );
}
#[test]
fn sam_and_paf_render_raw_ids_without_unicode_transcoding() {
    let r = context();
    let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1000).unwrap();
    let mut out = Vec::new();
    super::sam::print_match_context_bytes(&mut out, &r, &matrix, false, false, false).unwrap();
    assert!(out.starts_with(b"q\xff\t0\ts\x80\t1\t255\t"));
    out.clear();
    super::paf::print_match_context_bytes(&mut out, &r, &matrix).unwrap();
    assert!(out.starts_with(b"q\xff\t1\t0\t0\t+\ts\x80\t1\t"));
}
#[test]
fn pairwise_query_intro_uses_source_info_c_string_extent() {
    let mut out = Vec::new();
    super::pairwise::print_query_intro_bytes(&mut out, b"q\xff\0full", 1, false).unwrap();
    assert_eq!(out, b"Query= q\xff\n\nLength=1\n\n");
}
#[test]
fn xml_headers_keep_raw_database_and_source_preescape_cut_position() {
    let mut out = Vec::new();
    super::xml::print_header_bytes(
        &mut out,
        2,
        "BLOSUM62",
        11,
        1,
        0.001,
        b"<\xff\x01tail",
        1,
        b"db\xff<&",
    )
    .unwrap();
    assert!(out
        .windows(b"<BlastOutput_db>db\xff<&</BlastOutput_db>".len())
        .any(|x| x == b"<BlastOutput_db>db\xff<&</BlastOutput_db>"));
    // Source takes escaped.substr(0, original.find('\1')), not the escaped
    // position. This deliberately preserves its literal byte-count quirk.
    assert!(out
        .windows(b"<BlastOutput_query-def>&l</BlastOutput_query-def>".len())
        .any(|x| x == b"<BlastOutput_query-def>&l</BlastOutput_query-def>"));
}
#[test]
fn legacy_context_bridge_checks_utf8_instead_of_replacing_bytes() {
    assert!(context().try_to_string_context().is_err());
    assert_eq!(
        crate::util::escape_sequences::XML.escape_bytes(b"\xff<&"),
        b"\xff&lt;&amp;"
    );
}
