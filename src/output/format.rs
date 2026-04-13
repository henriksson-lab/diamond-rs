use std::io::Write;

/// Format a double value matching DIAMOND's C++ `format_double()`.
///
/// For values >= 100: floor to integer.
/// For values < 100: round to 1 decimal place, format as "X.Y".
pub fn format_double(x: f64) -> String {
    if x >= 100.0 {
        format!("{}", x.floor() as i64)
    } else if x <= -100.0 {
        format!("{}", x.ceil() as i64)
    } else {
        let rounded = (x * 10.0).round() as i64;
        let int_part = rounded / 10;
        let frac_part = (rounded % 10).abs();
        // Preserve sign for values like -0.5 where int_part rounds to 0
        if x < 0.0 && int_part == 0 {
            format!("-0.{}", frac_part)
        } else {
            format!("{}.{}", int_part, frac_part)
        }
    }
}

/// Format an e-value matching DIAMOND's C++ `print_e()`.
///
/// Zero → "0.0", otherwise `%.2e` format.
pub fn format_evalue(x: f64) -> String {
    if x == 0.0 {
        "0.0".to_string()
    } else {
        format!("{:.2e}", x)
    }
}

/// Output format codes matching the C++ enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FormatCode {
    Pairwise = 0,
    Xml = 5,
    Tabular = 6,
    Daa = 100,
    Sam = 101,
    Taxon = 102,
    Paf = 103,
    JsonFlat = 104,
}

impl FormatCode {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "0" => Some(Self::Pairwise),
            "5" | "xml" => Some(Self::Xml),
            "6" | "tab" => Some(Self::Tabular),
            "100" | "daa" => Some(Self::Daa),
            "101" | "sam" => Some(Self::Sam),
            "102" => Some(Self::Taxon),
            "103" | "paf" => Some(Self::Paf),
            "104" | "json-flat" => Some(Self::JsonFlat),
            _ => None,
        }
    }
}

/// Field IDs for tabular output (matching C++ FieldId enum).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum FieldId {
    QSeqId = 0,
    QGi = 1,
    QAcc = 2,
    QAccVer = 3,
    QLen = 4,
    SSeqId = 5,
    SAllSeqId = 6,
    SGi = 7,
    SAcc = 8,
    SAccVer = 9,
    SAllacc = 10,
    Stitle0 = 11,
    SLen = 12,
    QStart = 13,
    QEnd = 14,
    SStart = 15,
    SEnd = 16,
    QSeq = 17,
    SSeq = 18,
    EValue = 19,
    BitScore = 20,
    Score = 21,
    Length = 22,
    PIdent = 23,
    NIdent = 24,
    Mismatch = 25,
    Positive = 26,
    GapOpen = 27,
    Gaps = 28,
    PPos = 29,
    QFrameField = 30,
    QFrame = 31,
    BTop = 33,
    STaxIds = 34,
    SSciNames = 35,
    STitle = 39,
    SAllTitles = 40,
    SStrand = 41,
    QCovS = 42,
    QCovHsp = 43,
    QTitle = 45,
    FullSSeq = 48,
    QQual = 49,
    QNum = 50,
    SNum = 51,
    SCovHsp = 52,
    FullQQual = 53,
    FullQSeq = 54,
    QSeqGapped = 55,
    SSeqGapped = 56,
    QStrand = 57,
    Cigar = 58,
    SKingdoms = 59,
    SPhylums = 60,
    FullQSeqMate = 62,
    QSeqTranslated = 63,
    HspNum = 66,
    NormalizedBitscore = 67,
    NormalizedNident = 68,
    ApproxPIdent = 71,
    CorrectedBitScore = 72,
    NegEValue = 73,
}

impl FieldId {
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "qseqid" => Some(Self::QSeqId),
            "qgi" => Some(Self::QGi),
            "qacc" => Some(Self::QAcc),
            "qaccver" => Some(Self::QAccVer),
            "qlen" => Some(Self::QLen),
            "sseqid" => Some(Self::SSeqId),
            "sallseqid" => Some(Self::SAllSeqId),
            "sgi" => Some(Self::SGi),
            "sacc" => Some(Self::SAcc),
            "saccver" => Some(Self::SAccVer),
            "sallacc" => Some(Self::SAllacc),
            "slen" => Some(Self::SLen),
            "qstart" => Some(Self::QStart),
            "qend" => Some(Self::QEnd),
            "sstart" => Some(Self::SStart),
            "send" => Some(Self::SEnd),
            "qseq" => Some(Self::QSeq),
            "sseq" => Some(Self::SSeq),
            "evalue" => Some(Self::EValue),
            "bitscore" => Some(Self::BitScore),
            "score" => Some(Self::Score),
            "length" => Some(Self::Length),
            "pident" => Some(Self::PIdent),
            "nident" => Some(Self::NIdent),
            "mismatch" => Some(Self::Mismatch),
            "positive" => Some(Self::Positive),
            "gapopen" => Some(Self::GapOpen),
            "gaps" => Some(Self::Gaps),
            "ppos" => Some(Self::PPos),
            "qframe" => Some(Self::QFrame),
            "btop" => Some(Self::BTop),
            "staxids" => Some(Self::STaxIds),
            "sscinames" => Some(Self::SSciNames),
            "stitle" => Some(Self::STitle),
            "salltitles" => Some(Self::SAllTitles),
            "sstrand" => Some(Self::SStrand),
            "qcovs" => Some(Self::QCovS),
            "qcovhsp" => Some(Self::QCovHsp),
            "qtitle" => Some(Self::QTitle),
            "full_sseq" => Some(Self::FullSSeq),
            "qqual" => Some(Self::QQual),
            "qnum" => Some(Self::QNum),
            "snum" => Some(Self::SNum),
            "scovhsp" => Some(Self::SCovHsp),
            "full_qqual" => Some(Self::FullQQual),
            "full_qseq" => Some(Self::FullQSeq),
            "qseq_gapped" => Some(Self::QSeqGapped),
            "sseq_gapped" => Some(Self::SSeqGapped),
            "qstrand" => Some(Self::QStrand),
            "cigar" => Some(Self::Cigar),
            "skingdoms" => Some(Self::SKingdoms),
            "sphylums" => Some(Self::SPhylums),
            "full_qseq_mate" => Some(Self::FullQSeqMate),
            "qseq_translated" => Some(Self::QSeqTranslated),
            "hspnum" => Some(Self::HspNum),
            "nident_normalized" => Some(Self::NormalizedNident),
            "approx_pident" => Some(Self::ApproxPIdent),
            "corrected_bitscore" => Some(Self::CorrectedBitScore),
            _ => None,
        }
    }
}

/// Default BLAST tabular fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
pub const DEFAULT_TABULAR_FIELDS: &[FieldId] = &[
    FieldId::QSeqId,
    FieldId::SSeqId,
    FieldId::PIdent,
    FieldId::Length,
    FieldId::Mismatch,
    FieldId::GapOpen,
    FieldId::QStart,
    FieldId::QEnd,
    FieldId::SStart,
    FieldId::SEnd,
    FieldId::EValue,
    FieldId::BitScore,
];

/// High-Scoring Pair (HSP) — a single alignment result.
#[derive(Debug, Clone)]
pub struct Hsp {
    pub score: i32,
    pub evalue: f64,
    pub bit_score: f64,
    pub query_range: (i32, i32),
    pub subject_range: (i32, i32),
    pub query_source_range: (i32, i32),
    pub subject_source_range: (i32, i32),
    pub frame: i32,
    pub length: i32,
    pub identities: i32,
    pub mismatches: i32,
    pub positives: i32,
    pub gap_openings: i32,
    pub gaps: i32,
}

impl Hsp {
    /// Percent identity.
    pub fn pident(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            100.0 * self.identities as f64 / self.length as f64
        }
    }

    /// Percent positive.
    pub fn ppos(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            100.0 * self.positives as f64 / self.length as f64
        }
    }
}

/// Write a tab-separated row for a BLAST-6 output format.
pub fn write_tabular_row<W: Write>(
    writer: &mut W,
    query_id: &str,
    subject_id: &str,
    hsp: &Hsp,
    fields: &[FieldId],
    query_len: i32,
    subject_len: i32,
) -> std::io::Result<()> {
    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            write!(writer, "\t")?;
        }
        match field {
            FieldId::QSeqId | FieldId::QAcc | FieldId::QAccVer => write!(writer, "{}", query_id)?,
            FieldId::SSeqId | FieldId::SAcc | FieldId::SAccVer => {
                write!(writer, "{}", subject_id)?
            }
            FieldId::PIdent => write!(writer, "{}", format_double(hsp.pident()))?,
            FieldId::Length => write!(writer, "{}", hsp.length)?,
            FieldId::Mismatch => write!(writer, "{}", hsp.mismatches)?,
            FieldId::GapOpen => write!(writer, "{}", hsp.gap_openings)?,
            FieldId::QStart => write!(writer, "{}", hsp.query_source_range.0 + 1)?,
            FieldId::QEnd => write!(writer, "{}", hsp.query_source_range.1)?,
            FieldId::SStart => write!(writer, "{}", hsp.subject_source_range.0 + 1)?,
            FieldId::SEnd => write!(writer, "{}", hsp.subject_source_range.1)?,
            FieldId::EValue => write!(writer, "{}", format_evalue(hsp.evalue))?,
            FieldId::BitScore => write!(writer, "{}", format_double(hsp.bit_score))?,
            FieldId::Score => write!(writer, "{}", hsp.score)?,
            FieldId::NIdent => write!(writer, "{}", hsp.identities)?,
            FieldId::Positive => write!(writer, "{}", hsp.positives)?,
            FieldId::Gaps => write!(writer, "{}", hsp.gaps)?,
            FieldId::PPos => write!(writer, "{}", format_double(hsp.ppos()))?,
            FieldId::QLen => write!(writer, "{}", query_len)?,
            FieldId::SLen => write!(writer, "{}", subject_len)?,
            FieldId::QFrame => write!(writer, "{}", hsp.frame)?,
            FieldId::QCovHsp => {
                let cov = if query_len > 0 {
                    100.0 * (hsp.query_range.1 - hsp.query_range.0) as f64 / query_len as f64
                } else {
                    0.0
                };
                write!(writer, "{}", format_double(cov))?;
            }
            FieldId::SCovHsp => {
                let cov = if subject_len > 0 {
                    100.0 * (hsp.subject_range.1 - hsp.subject_range.0) as f64
                        / subject_len as f64
                } else {
                    0.0
                };
                write!(writer, "{}", format_double(cov))?;
            }
            _ => write!(writer, "N/A")?,
        }
    }
    writeln!(writer)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_double() {
        // >= 100: floor to integer
        assert_eq!(format_double(679.1), "679");
        assert_eq!(format_double(100.0), "100");
        assert_eq!(format_double(100.9), "100");
        // < 100: one decimal
        assert_eq!(format_double(71.7), "71.7");
        assert_eq!(format_double(99.9), "99.9");
        assert_eq!(format_double(0.0), "0.0");
        assert_eq!(format_double(50.55), "50.6"); // rounds
    }

    #[test]
    fn test_format_evalue() {
        assert_eq!(format_evalue(0.0), "0.0");
        assert_eq!(format_evalue(2.23e-247), "2.23e-247");
        assert_eq!(format_evalue(1.5e-5), "1.50e-5");
    }

    #[test]
    fn test_format_code_parse() {
        assert_eq!(FormatCode::parse("0"), Some(FormatCode::Pairwise));
        assert_eq!(FormatCode::parse("6"), Some(FormatCode::Tabular));
        assert_eq!(FormatCode::parse("tab"), Some(FormatCode::Tabular));
        assert_eq!(FormatCode::parse("xml"), Some(FormatCode::Xml));
        assert_eq!(FormatCode::parse("daa"), Some(FormatCode::Daa));
        assert_eq!(FormatCode::parse("paf"), Some(FormatCode::Paf));
        assert_eq!(FormatCode::parse("invalid"), None);
    }

    #[test]
    fn test_field_id_parse() {
        assert_eq!(FieldId::from_name("qseqid"), Some(FieldId::QSeqId));
        assert_eq!(FieldId::from_name("sseqid"), Some(FieldId::SSeqId));
        assert_eq!(FieldId::from_name("evalue"), Some(FieldId::EValue));
        assert_eq!(FieldId::from_name("bitscore"), Some(FieldId::BitScore));
        assert_eq!(FieldId::from_name("unknown"), None);
    }

    #[test]
    fn test_default_fields() {
        assert_eq!(DEFAULT_TABULAR_FIELDS.len(), 12);
    }

    #[test]
    fn test_hsp_pident() {
        let hsp = Hsp {
            score: 100,
            evalue: 1e-10,
            bit_score: 50.0,
            query_range: (0, 100),
            subject_range: (0, 100),
            query_source_range: (0, 100),
            subject_source_range: (0, 100),
            frame: 0,
            length: 100,
            identities: 80,
            mismatches: 20,
            positives: 90,
            gap_openings: 0,
            gaps: 0,
        };
        assert!((hsp.pident() - 80.0).abs() < 0.001);
        assert!((hsp.ppos() - 90.0).abs() < 0.001);
    }

    #[test]
    fn test_write_tabular() {
        let hsp = Hsp {
            score: 50,
            evalue: 1.5e-5,
            bit_score: 30.4,
            query_range: (0, 50),
            subject_range: (10, 60),
            query_source_range: (0, 50),
            subject_source_range: (10, 60),
            frame: 0,
            length: 50,
            identities: 40,
            mismatches: 10,
            positives: 45,
            gap_openings: 0,
            gaps: 0,
        };
        let mut buf = Vec::new();
        write_tabular_row(
            &mut buf,
            "query1",
            "subject1",
            &hsp,
            DEFAULT_TABULAR_FIELDS,
            100,
            200,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("query1\tsubject1\t80.0"));
        assert!(output.contains("1.50e-5") || output.contains("1.50e-05"));
    }
}
