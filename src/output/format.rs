use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{TaxId, AMINO_ACID_ALPHABET, LETTER_MASK};
use crate::data::taxonomy::TaxonomyTree;
use crate::dp::swipe::HspValues;
use crate::util::escape_sequences::XML;
use crate::util::sequence::{FASTA_HEADER_SEP, ID_DELIMITERS};

/// Matches C++ `format_double(x)`.
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

/// Matches C++ `print_e(x)`.
///
/// Zero → "0.0". Otherwise scientific with 2 fractional digits, an explicit
/// `+`/`-` exponent sign, and the exponent padded to at least 2 digits.
/// Rust's `{:.2e}` omits the sign for positive exponents and the leading zero
/// for single-digit exponents (`8.35e-9` instead of `8.35e-09`).
pub fn format_evalue(x: f64) -> String {
    if x == 0.0 {
        return "0.0".to_string();
    }
    let s = format!("{:.2e}", x);
    let bytes = s.as_bytes();
    let Some(epos) = bytes.iter().position(|&b| b == b'e') else {
        return s;
    };
    let (mantissa, exp) = s.split_at(epos);
    let exp = &exp[1..]; // strip 'e'
    let (sign, digits) = match exp.as_bytes().first() {
        Some(b'-') => ("-", &exp[1..]),
        Some(b'+') => ("+", &exp[1..]),
        _ => ("+", exp),
    };
    if digits.len() == 1 {
        format!("{}e{}0{}", mantissa, sign, digits)
    } else {
        format!("{}e{}{}", mantissa, sign, digits)
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct OutputFlags(pub u32);

impl OutputFlags {
    pub const NONE: OutputFlags = OutputFlags(0);
    pub const FULL_TITLES: OutputFlags = OutputFlags(1);
    pub const ALL_SEQIDS: OutputFlags = OutputFlags(1 << 1);
    pub const TARGET_SEQS: OutputFlags = OutputFlags(1 << 2);
    pub const SELF_ALN_SCORES: OutputFlags = OutputFlags(1 << 3);
    pub const IS_STRING: OutputFlags = OutputFlags(1 << 4);
    pub const IS_ARRAY: OutputFlags = OutputFlags(1 << 5);
    pub const SSEQID: OutputFlags = OutputFlags(1 << 6);
    pub const DEFAULT_REPORT_UNALIGNED: OutputFlags = OutputFlags(1 << 7);
    pub const NO_REALIGN: OutputFlags = OutputFlags(1 << 8);

    pub fn any(self, rhs: OutputFlags) -> bool {
        self.0 & rhs.0 != 0
    }
}

impl std::ops::BitOr for OutputFlags {
    type Output = OutputFlags;

    fn bitor(self, rhs: OutputFlags) -> OutputFlags {
        OutputFlags(self.0 | rhs.0)
    }
}

impl std::ops::BitOrAssign for OutputFlags {
    fn bitor_assign(&mut self, rhs: OutputFlags) {
        self.0 |= rhs.0;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormatKind {
    Daa,
    BlastTab,
    BlastXml,
    Sam,
    BlastPairwise,
    Null,
    Taxon,
    Paf,
    Edge,
    Json,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputCommand {
    BlastP,
    BlastX,
    View,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OutputFormatSpec {
    pub code: OutputFormatKind,
    pub query_separator: char,
    pub needs_taxon_id_lists: bool,
    pub needs_taxon_nodes: bool,
    pub needs_taxon_scientific_names: bool,
    pub needs_taxon_ranks: bool,
    pub needs_paired_end_info: bool,
    pub hsp_values: HspValues,
    pub flags: OutputFlags,
}

impl OutputFormatSpec {
    pub fn new(
        code: OutputFormatKind,
        hsp_values: HspValues,
        flags: OutputFlags,
        query_separator: char,
    ) -> Self {
        Self {
            code,
            query_separator,
            needs_taxon_id_lists: false,
            needs_taxon_nodes: false,
            needs_taxon_scientific_names: false,
            needs_taxon_ranks: false,
            needs_paired_end_info: false,
            hsp_values,
            flags,
        }
    }

    pub fn null_format() -> Self {
        Self::new(
            OutputFormatKind::Null,
            HspValues::TRANSCRIPT,
            OutputFlags::NONE,
            '\0',
        )
    }

    pub fn daa_format(salltitles: bool, sallseqid: bool) -> Self {
        let flags = OutputFlags::SSEQID
            | if salltitles {
                OutputFlags::FULL_TITLES | OutputFlags::ALL_SEQIDS
            } else if sallseqid {
                OutputFlags::ALL_SEQIDS
            } else {
                OutputFlags::NONE
            };
        Self::new(OutputFormatKind::Daa, HspValues::TRANSCRIPT, flags, '\0')
    }

    pub fn paf_format() -> Self {
        Self::new(
            OutputFormatKind::Paf,
            HspValues::TRANSCRIPT,
            OutputFlags::SSEQID | OutputFlags::DEFAULT_REPORT_UNALIGNED,
            '\0',
        )
    }

    pub fn sam_format() -> Self {
        Self::new(
            OutputFormatKind::Sam,
            HspValues::TRANSCRIPT,
            OutputFlags::SSEQID | OutputFlags::DEFAULT_REPORT_UNALIGNED,
            '\0',
        )
    }

    pub fn xml_format() -> Self {
        Self::new(
            OutputFormatKind::BlastXml,
            HspValues::TRANSCRIPT,
            OutputFlags::FULL_TITLES
                | OutputFlags::SSEQID
                | OutputFlags::DEFAULT_REPORT_UNALIGNED
                | OutputFlags::ALL_SEQIDS,
            '\0',
        )
    }

    pub fn pairwise_format() -> Self {
        Self::new(
            OutputFormatKind::BlastPairwise,
            HspValues::TRANSCRIPT,
            OutputFlags::FULL_TITLES
                | OutputFlags::SSEQID
                | OutputFlags::DEFAULT_REPORT_UNALIGNED
                | OutputFlags::ALL_SEQIDS,
            '\0',
        )
    }

    pub fn taxon_format(include_lineage: bool) -> Self {
        let mut format = Self::new(
            OutputFormatKind::Taxon,
            HspValues::NONE,
            OutputFlags::DEFAULT_REPORT_UNALIGNED,
            '\0',
        );
        format.needs_taxon_id_lists = true;
        format.needs_taxon_nodes = true;
        format.needs_taxon_scientific_names = include_lineage;
        format
    }

    pub fn edge_format() -> Self {
        Self::new(
            OutputFormatKind::Edge,
            HspValues::COORDS,
            OutputFlags::SSEQID,
            '\0',
        )
    }

    pub fn report_unaligned(&self, report_unaligned: i32) -> bool {
        report_unaligned != 0
            && (report_unaligned != -1 || self.flags.any(OutputFlags::DEFAULT_REPORT_UNALIGNED))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResolvedOutputFormat {
    Tabular(TabularFormat),
    Spec(OutputFormatSpec),
}

pub const DEFAULT_MAX_TARGET_SEQS: i64 = 25;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutputInitConfig<'a> {
    pub output_format: &'a [&'a str],
    pub daa_file: &'a str,
    pub command: OutputCommand,
    pub include_lineage: bool,
    pub frame_shift: i32,
    pub query_translated: bool,
    pub score_matrix_name: &'a str,
    pub ext: &'a str,
    pub multiprocessing: bool,
    pub global_ranking_targets: u32,
    pub toppercent: Option<f64>,
    pub max_target_seqs: Option<i64>,
    pub min_bit_score: f64,
    pub query_range_culling: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OutputInit {
    pub output_format: ResolvedOutputFormat,
    pub max_target_seqs: i64,
    pub toppercent: Option<f64>,
    pub message: String,
    pub dp_fields: HspValues,
}

/// Matches C++ `get_output_format(output_format, daa_file, command)`.
pub fn get_output_format(
    output_format: &[&str],
    daa_file: &str,
    command: OutputCommand,
    include_lineage: bool,
    frame_shift: i32,
    query_translated: bool,
    score_matrix_name: &str,
) -> Result<ResolvedOutputFormat, String> {
    if output_format.is_empty() {
        if daa_file.is_empty() || command == OutputCommand::View {
            return Ok(ResolvedOutputFormat::Tabular(TabularFormat::new(
                output_format,
                false,
                frame_shift,
                query_translated,
                score_matrix_name,
            )?));
        } else if matches!(command, OutputCommand::BlastP | OutputCommand::BlastX)
            && !daa_file.is_empty()
        {
            return Ok(ResolvedOutputFormat::Spec(OutputFormatSpec::daa_format(
                false, false,
            )));
        }
    }

    let Some(&format) = output_format.first() else {
        return Err("Invalid output format".to_string());
    };
    match format {
        "tab" | "6" => Ok(ResolvedOutputFormat::Tabular(TabularFormat::new(
            output_format,
            false,
            frame_shift,
            query_translated,
            score_matrix_name,
        )?)),
        "sam" | "101" => Ok(ResolvedOutputFormat::Spec(OutputFormatSpec::sam_format())),
        "xml" | "5" => Ok(ResolvedOutputFormat::Spec(OutputFormatSpec::xml_format())),
        "daa" | "100" => Ok(ResolvedOutputFormat::Spec(OutputFormatSpec::daa_format(
            false, false,
        ))),
        "0" => Ok(ResolvedOutputFormat::Spec(
            OutputFormatSpec::pairwise_format(),
        )),
        "null" => Ok(ResolvedOutputFormat::Spec(OutputFormatSpec::null_format())),
        "102" => Ok(ResolvedOutputFormat::Spec(OutputFormatSpec::taxon_format(
            include_lineage,
        ))),
        "paf" | "103" => Ok(ResolvedOutputFormat::Spec(OutputFormatSpec::paf_format())),
        "edge" => Ok(ResolvedOutputFormat::Spec(OutputFormatSpec::edge_format())),
        "json-flat" | "104" => Ok(ResolvedOutputFormat::Tabular(TabularFormat::new(
            output_format,
            true,
            frame_shift,
            query_translated,
            score_matrix_name,
        )?)),
        _ => Err(format!(
            "Invalid output format: {}\nAllowed values: 0,5,xml,6,tab,100,daa,101,sam,102,103,104,paf",
            format
        )),
    }
}

/// Matches C++ `init_output(int64_t& max_target_seqs)`.
pub fn init_output(cfg: &OutputInitConfig<'_>) -> Result<OutputInit, String> {
    let mut output_format = get_output_format(
        cfg.output_format,
        cfg.daa_file,
        cfg.command,
        cfg.include_lineage,
        cfg.frame_shift,
        cfg.query_translated,
        cfg.score_matrix_name,
    )?;

    if cfg.ext == "none" && cfg.frame_shift != 0 {
        return Err("Frameshift alignment does not support --ext none.".to_string());
    }
    if cfg.ext == "none"
        && matches!(
            output_format,
            ResolvedOutputFormat::Spec(OutputFormatSpec {
                code: OutputFormatKind::Daa,
                ..
            })
        )
    {
        return Err("--ext none is not supported for DAA output.".to_string());
    }

    let needs_taxonomy = match &output_format {
        ResolvedOutputFormat::Tabular(format) => {
            format.needs_taxon_id_lists
                || format.needs_taxon_nodes
                || format.needs_taxon_scientific_names
        }
        ResolvedOutputFormat::Spec(format) => {
            format.needs_taxon_id_lists
                || format.needs_taxon_nodes
                || format.needs_taxon_scientific_names
        }
    };
    if cfg.command == OutputCommand::View && needs_taxonomy {
        return Err("Taxonomy features are not supported for the DAA format.".to_string());
    }

    if matches!(
        output_format,
        ResolvedOutputFormat::Spec(OutputFormatSpec {
            code: OutputFormatKind::Daa,
            ..
        })
    ) && cfg.multiprocessing
    {
        return Err("The DAA format is not supported in multiprocessing mode.".to_string());
    }
    if matches!(
        output_format,
        ResolvedOutputFormat::Spec(OutputFormatSpec {
            code: OutputFormatKind::Daa,
            ..
        })
    ) && cfg.global_ranking_targets != 0
    {
        return Err("The DAA format is not supported in global ranking mode.".to_string());
    }

    let mut toppercent = cfg.toppercent;
    if matches!(
        output_format,
        ResolvedOutputFormat::Spec(OutputFormatSpec {
            code: OutputFormatKind::Taxon,
            ..
        })
    ) && toppercent.is_none()
        && cfg.min_bit_score == 0.0
    {
        toppercent = Some(10.0);
    }

    let max_target_seqs;
    let message;
    if let Some(top) = toppercent {
        if !(0.0..=100.0).contains(&top) {
            return Err("Allowed value range for --top is between 0.0 and 100.0".to_string());
        }
        if top == 100.0 {
            toppercent = None;
            max_target_seqs = i64::MAX;
            message = "#Target sequences to report alignments for: unlimited\n".to_string();
        } else {
            max_target_seqs = cfg.max_target_seqs.unwrap_or(DEFAULT_MAX_TARGET_SEQS);
            message = format!("Percentage range of top alignment score to report hits: {top}\n");
        }
    } else {
        let max = cfg.max_target_seqs.unwrap_or(DEFAULT_MAX_TARGET_SEQS);
        max_target_seqs = if cfg.max_target_seqs == Some(0) {
            i64::MAX
        } else {
            max
        };
        message = if max_target_seqs == i64::MAX {
            "#Target sequences to report alignments for: unlimited\n".to_string()
        } else {
            format!("#Target sequences to report alignments for: {max_target_seqs}\n")
        };
    }

    if cfg.frame_shift != 0 {
        match &mut output_format {
            ResolvedOutputFormat::Tabular(format)
                if format.hsp_values != HspValues::NONE || cfg.query_range_culling =>
            {
                format.hsp_values = HspValues::TRANSCRIPT;
            }
            ResolvedOutputFormat::Spec(format)
                if format.hsp_values != HspValues::NONE || cfg.query_range_culling =>
            {
                format.hsp_values = HspValues::TRANSCRIPT;
            }
            _ => {}
        }
    }

    let dp_fields = match &output_format {
        ResolvedOutputFormat::Tabular(format) => format.hsp_values,
        ResolvedOutputFormat::Spec(format) => format.hsp_values,
    };

    Ok(OutputInit {
        output_format,
        max_target_seqs,
        toppercent,
        message,
        dp_fields,
    })
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
    SSKingdoms = 38,
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
    Reserved1 = 74,
    Reserved2 = 75,
    SLineages = 76,
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
            "sskingdoms" => Some(Self::SSKingdoms),
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
            // C++ `blast_tab_format.cpp:104-105` registers these as
            // `normalized_bitscore` and `normalized_nident`. Accept both
            // the C++ name and the prior Rust-only `nident_normalized` spelling
            // (so existing scripts keep working) and add `normalized_bitscore`
            // which was missing from `from_name` entirely.
            "normalized_nident" | "nident_normalized" => Some(Self::NormalizedNident),
            "normalized_bitscore" => Some(Self::NormalizedBitscore),
            "approx_pident" => Some(Self::ApproxPIdent),
            "corrected_bitscore" => Some(Self::CorrectedBitScore),
            "slineages" => Some(Self::SLineages),
            "reserved1" => Some(Self::Reserved1),
            "reserved2" => Some(Self::Reserved2),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OutputField {
    pub id: FieldId,
    pub key: &'static str,
    pub clust_key: &'static str,
    pub description: &'static str,
    pub is_string: bool,
    pub is_array: bool,
}

impl FieldId {
    pub fn output_field(self) -> Option<OutputField> {
        Some(match self {
            FieldId::QSeqId => OutputField {
                id: self,
                key: "qseqid",
                clust_key: "cseqid",
                description: "Query Seq - id",
                is_string: true,
                is_array: false,
            },
            FieldId::QLen => OutputField {
                id: self,
                key: "qlen",
                clust_key: "clen",
                description: "Query sequence length",
                is_string: false,
                is_array: false,
            },
            FieldId::SSeqId => OutputField {
                id: self,
                key: "sseqid",
                clust_key: "mseqid",
                description: "Subject Seq - id",
                is_string: true,
                is_array: false,
            },
            FieldId::SAllSeqId => OutputField {
                id: self,
                key: "sallseqid",
                clust_key: "",
                description: "All subject Seq - id(s), separated by a ';'",
                is_string: false,
                is_array: true,
            },
            FieldId::SLen => OutputField {
                id: self,
                key: "slen",
                clust_key: "mlen",
                description: "Subject sequence length",
                is_string: false,
                is_array: false,
            },
            FieldId::QStart => OutputField {
                id: self,
                key: "qstart",
                clust_key: "cstart",
                description: "Start of alignment in query",
                is_string: false,
                is_array: false,
            },
            FieldId::QEnd => OutputField {
                id: self,
                key: "qend",
                clust_key: "cend",
                description: "End of alignment in query",
                is_string: false,
                is_array: false,
            },
            FieldId::SStart => OutputField {
                id: self,
                key: "sstart",
                clust_key: "mstart",
                description: "Start of alignment in subject",
                is_string: false,
                is_array: false,
            },
            FieldId::SEnd => OutputField {
                id: self,
                key: "send",
                clust_key: "mend",
                description: "End of alignment in subject",
                is_string: false,
                is_array: false,
            },
            FieldId::QSeq => OutputField {
                id: self,
                key: "qseq",
                clust_key: "",
                description: "Aligned part of query sequence",
                is_string: true,
                is_array: false,
            },
            FieldId::SSeq => OutputField {
                id: self,
                key: "sseq",
                clust_key: "",
                description: "Aligned part of subject sequence",
                is_string: true,
                is_array: false,
            },
            FieldId::EValue => OutputField {
                id: self,
                key: "evalue",
                clust_key: "evalue",
                description: "Expect value",
                is_string: false,
                is_array: false,
            },
            FieldId::BitScore => OutputField {
                id: self,
                key: "bitscore",
                clust_key: "Bitscore",
                description: "Bit score",
                is_string: false,
                is_array: false,
            },
            FieldId::Score => OutputField {
                id: self,
                key: "score",
                clust_key: "score",
                description: "Raw score",
                is_string: false,
                is_array: false,
            },
            FieldId::Length => OutputField {
                id: self,
                key: "length",
                clust_key: "length",
                description: "Alignment length",
                is_string: false,
                is_array: false,
            },
            FieldId::PIdent => OutputField {
                id: self,
                key: "pident",
                clust_key: "pident",
                description: "Percentage of identical matches",
                is_string: false,
                is_array: false,
            },
            FieldId::NIdent => OutputField {
                id: self,
                key: "nident",
                clust_key: "nident",
                description: "Number of identical matches",
                is_string: false,
                is_array: false,
            },
            FieldId::Mismatch => OutputField {
                id: self,
                key: "mismatch",
                clust_key: "mismatch",
                description: "Number of mismatches",
                is_string: false,
                is_array: false,
            },
            FieldId::Positive => OutputField {
                id: self,
                key: "positive",
                clust_key: "positive",
                description: "Number of positive - scoring matches",
                is_string: false,
                is_array: false,
            },
            FieldId::GapOpen => OutputField {
                id: self,
                key: "gapopen",
                clust_key: "gapopen",
                description: "Number of gap openings",
                is_string: false,
                is_array: false,
            },
            FieldId::Gaps => OutputField {
                id: self,
                key: "gaps",
                clust_key: "gaps",
                description: "Total number of gaps",
                is_string: false,
                is_array: false,
            },
            FieldId::PPos => OutputField {
                id: self,
                key: "ppos",
                clust_key: "ppos",
                description: "Percentage of positive - scoring matches",
                is_string: false,
                is_array: false,
            },
            FieldId::QFrame => OutputField {
                id: self,
                key: "qframe",
                clust_key: "",
                description: "Query frame",
                is_string: false,
                is_array: false,
            },
            FieldId::BTop => OutputField {
                id: self,
                key: "btop",
                clust_key: "",
                description: "Blast traceback operations (BTOP)",
                is_string: true,
                is_array: false,
            },
            FieldId::STaxIds => OutputField {
                id: self,
                key: "staxids",
                clust_key: "",
                description:
                    "Unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)",
                is_string: false,
                is_array: true,
            },
            FieldId::SSciNames => OutputField {
                id: self,
                key: "sscinames",
                clust_key: "",
                description: "Unique Subject Scientific Name(s), separated by a ';'",
                is_string: false,
                is_array: true,
            },
            FieldId::SSKingdoms => OutputField {
                id: self,
                key: "sskingdoms",
                clust_key: "",
                description: "Unique Subject Super Kingdom(s), separated by a ';'",
                is_string: false,
                is_array: true,
            },
            FieldId::STitle => OutputField {
                id: self,
                key: "stitle",
                clust_key: "",
                description: "Subject Title",
                is_string: true,
                is_array: false,
            },
            FieldId::SAllTitles => OutputField {
                id: self,
                key: "salltitles",
                clust_key: "",
                description: "All Subject Title(s), separated by a '<>'",
                is_string: false,
                is_array: true,
            },
            FieldId::QCovHsp => OutputField {
                id: self,
                key: "qcovhsp",
                clust_key: "ccovhsp",
                description: "Query coverage per HSP",
                is_string: false,
                is_array: false,
            },
            FieldId::QTitle => OutputField {
                id: self,
                key: "qtitle",
                clust_key: "",
                description: "Query title",
                is_string: true,
                is_array: false,
            },
            FieldId::FullSSeq => OutputField {
                id: self,
                key: "full_sseq",
                clust_key: "",
                description: "Subject sequence",
                is_string: true,
                is_array: false,
            },
            FieldId::QQual => OutputField {
                id: self,
                key: "qqual",
                clust_key: "",
                description: "Query quality values for the aligned part of the query",
                is_string: true,
                is_array: false,
            },
            FieldId::QNum => OutputField {
                id: self,
                key: "qnum",
                clust_key: "",
                description: "Query ordinal id",
                is_string: false,
                is_array: false,
            },
            FieldId::SNum => OutputField {
                id: self,
                key: "snum",
                clust_key: "",
                description: "Subject ordinal id",
                is_string: false,
                is_array: false,
            },
            FieldId::SCovHsp => OutputField {
                id: self,
                key: "scovhsp",
                clust_key: "mcovhsp",
                description: "Subject coverage per HSP",
                is_string: false,
                is_array: false,
            },
            FieldId::FullQQual => OutputField {
                id: self,
                key: "full_qqual",
                clust_key: "",
                description: "Query quality values",
                is_string: true,
                is_array: false,
            },
            FieldId::FullQSeq => OutputField {
                id: self,
                key: "full_qseq",
                clust_key: "",
                description: "Query sequence",
                is_string: true,
                is_array: false,
            },
            FieldId::QSeqGapped => OutputField {
                id: self,
                key: "qseq_gapped",
                clust_key: "",
                description: "Aligned part of query sequence (with gaps)",
                is_string: true,
                is_array: false,
            },
            FieldId::SSeqGapped => OutputField {
                id: self,
                key: "sseq_gapped",
                clust_key: "",
                description: "Aligned part of subject sequence (with gaps)",
                is_string: true,
                is_array: false,
            },
            FieldId::QStrand => OutputField {
                id: self,
                key: "qstrand",
                clust_key: "",
                description: "Query strand",
                is_string: true,
                is_array: false,
            },
            FieldId::Cigar => OutputField {
                id: self,
                key: "cigar",
                clust_key: "",
                description: "CIGAR string",
                is_string: true,
                is_array: false,
            },
            FieldId::SKingdoms => OutputField {
                id: self,
                key: "skingdoms",
                clust_key: "",
                description: "Unique Subject Kingdom(s), separated by a ';'",
                is_string: false,
                is_array: true,
            },
            FieldId::SPhylums => OutputField {
                id: self,
                key: "sphylums",
                clust_key: "",
                description: "Unique Subject Phylum(s), separated by a ';'",
                is_string: false,
                is_array: true,
            },
            FieldId::FullQSeqMate => OutputField {
                id: self,
                key: "full_qseq_mate",
                clust_key: "",
                description: "Query sequence of the mate",
                is_string: true,
                is_array: false,
            },
            FieldId::QSeqTranslated => OutputField {
                id: self,
                key: "qseq_translated",
                clust_key: "",
                description: "Aligned part of query sequence (translated)",
                is_string: true,
                is_array: false,
            },
            FieldId::HspNum => OutputField {
                id: self,
                key: "hspnum",
                clust_key: "",
                description: "Number of HSP within the subject",
                is_string: false,
                is_array: false,
            },
            FieldId::NormalizedBitscore => OutputField {
                id: self,
                key: "normalized_bitscore",
                clust_key: "",
                description: "Bitscore normalized by maximum self alignment score",
                is_string: false,
                is_array: false,
            },
            FieldId::NormalizedNident => OutputField {
                id: self,
                key: "normalized_nident",
                clust_key: "normalized_nident",
                description: "Number of identical matches normalized by maximum length",
                is_string: false,
                is_array: false,
            },
            FieldId::ApproxPIdent => OutputField {
                id: self,
                key: "approx_pident",
                clust_key: "approx_pident",
                description: "Approximate percentage of identical matches",
                is_string: false,
                is_array: false,
            },
            FieldId::CorrectedBitScore => OutputField {
                id: self,
                key: "corrected_bitscore",
                clust_key: "corrected_bitscore",
                description: "Bit score corrected for edge effects",
                is_string: false,
                is_array: false,
            },
            FieldId::SLineages => OutputField {
                id: self,
                key: "slineages",
                clust_key: "",
                description: "Unique Subject Lineage(s), separated by a '<>'",
                is_string: false,
                is_array: false,
            },
            _ => return None,
        })
    }

    pub fn hsp_values(self) -> HspValues {
        match self {
            FieldId::QStart => HspValues::QUERY_START,
            FieldId::QEnd => HspValues::QUERY_END,
            FieldId::SStart => HspValues::TARGET_START,
            FieldId::SEnd => HspValues::TARGET_END,
            FieldId::QSeq | FieldId::QQual | FieldId::QCovHsp => HspValues::QUERY_COORDS,
            FieldId::Length => HspValues::LENGTH,
            FieldId::PIdent | FieldId::NormalizedNident => HspValues::IDENT | HspValues::LENGTH,
            FieldId::NIdent => HspValues::IDENT,
            FieldId::Mismatch => HspValues::MISMATCHES,
            FieldId::Positive
            | FieldId::PPos
            | FieldId::SSeq
            | FieldId::BTop
            | FieldId::QSeqGapped
            | FieldId::SSeqGapped
            | FieldId::Cigar
            | FieldId::QSeqTranslated => HspValues::TRANSCRIPT,
            FieldId::GapOpen => HspValues::GAP_OPENINGS,
            FieldId::Gaps => HspValues::GAPS,
            FieldId::SCovHsp => HspValues::TARGET_COORDS,
            FieldId::ApproxPIdent => HspValues::COORDS,
            _ => HspValues::NONE,
        }
    }

    pub fn output_flags(self) -> OutputFlags {
        match self {
            FieldId::QSeqId | FieldId::QTitle => OutputFlags::IS_STRING,
            FieldId::QSeq
            | FieldId::SSeq
            | FieldId::BTop
            | FieldId::FullSSeq
            | FieldId::QQual
            | FieldId::FullQQual
            | FieldId::FullQSeq
            | FieldId::QSeqGapped
            | FieldId::SSeqGapped
            | FieldId::QStrand
            | FieldId::Cigar
            | FieldId::FullQSeqMate
            | FieldId::QSeqTranslated => OutputFlags::IS_STRING | OutputFlags::NO_REALIGN,
            FieldId::SSeqId => OutputFlags::IS_STRING | OutputFlags::SSEQID,
            FieldId::SAllSeqId => {
                OutputFlags::ALL_SEQIDS | OutputFlags::IS_ARRAY | OutputFlags::SSEQID
            }
            FieldId::STaxIds
            | FieldId::SSciNames
            | FieldId::SSKingdoms
            | FieldId::SKingdoms
            | FieldId::SPhylums => OutputFlags::IS_ARRAY | OutputFlags::NO_REALIGN,
            FieldId::STitle => {
                OutputFlags::FULL_TITLES | OutputFlags::IS_STRING | OutputFlags::SSEQID
            }
            FieldId::SAllTitles => {
                OutputFlags::ALL_SEQIDS
                    | OutputFlags::FULL_TITLES
                    | OutputFlags::IS_ARRAY
                    | OutputFlags::SSEQID
            }
            FieldId::NormalizedBitscore => OutputFlags::SELF_ALN_SCORES | OutputFlags::NO_REALIGN,
            FieldId::SLineages | FieldId::Reserved1 | FieldId::Reserved2 => OutputFlags::NO_REALIGN,
            FieldId::QFrame => OutputFlags::NO_REALIGN,
            _ => OutputFlags::NONE,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Header {
    None,
    Simple,
    Verbose,
}

impl Header {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "0" => Some(Self::None),
            "simple" => Some(Self::Simple),
            "verbose" => Some(Self::Verbose),
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Workflow {
    BlastP,
    Cluster,
    DeepClust,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TabularFormat {
    pub code: FormatCode,
    pub query_separator: char,
    pub fields: Vec<FieldId>,
    pub is_json: bool,
    pub hsp_values: HspValues,
    pub flags: OutputFlags,
    pub needs_taxon_id_lists: bool,
    pub needs_taxon_nodes: bool,
    pub needs_taxon_scientific_names: bool,
    pub needs_taxon_ranks: bool,
    pub needs_paired_end_info: bool,
    pub store_query_quality: bool,
}

impl TabularFormat {
    pub fn new(
        output_format: &[&str],
        json: bool,
        frame_shift: i32,
        query_translated: bool,
        score_matrix_name: &str,
    ) -> Result<Self, String> {
        let mut format = Self {
            code: if json {
                FormatCode::JsonFlat
            } else {
                FormatCode::Tabular
            },
            query_separator: if json { ',' } else { '\0' },
            fields: Vec::new(),
            is_json: json,
            hsp_values: HspValues::NONE,
            flags: OutputFlags::NONE,
            needs_taxon_id_lists: false,
            needs_taxon_nodes: false,
            needs_taxon_scientific_names: false,
            needs_taxon_ranks: false,
            needs_paired_end_info: false,
            store_query_quality: false,
        };
        if output_format.len() <= 1 {
            format.fields = DEFAULT_TABULAR_FIELDS.to_vec();
            format.hsp_values = if frame_shift == 0 {
                HspValues::QUERY_COORDS
                    | HspValues::TARGET_COORDS
                    | HspValues::LENGTH
                    | HspValues::IDENT
                    | HspValues::MISMATCHES
                    | HspValues::GAP_OPENINGS
            } else {
                HspValues::TRANSCRIPT
            };
            format.flags |= OutputFlags::SSEQID;
            return Ok(format);
        }
        for field in &output_format[1..] {
            let id = FieldId::from_name(field)
                .ok_or_else(|| format!("Invalid output field: {}", field))?;
            if id == FieldId::STaxIds {
                format.needs_taxon_id_lists = true;
            }
            if matches!(
                id,
                FieldId::SSciNames
                    | FieldId::SSKingdoms
                    | FieldId::SKingdoms
                    | FieldId::SPhylums
                    | FieldId::SLineages
            ) {
                format.needs_taxon_scientific_names = true;
                format.needs_taxon_id_lists = true;
            }
            if matches!(
                id,
                FieldId::SSKingdoms | FieldId::SKingdoms | FieldId::SPhylums
            ) {
                format.needs_taxon_nodes = true;
                format.needs_taxon_ranks = true;
            }
            if id == FieldId::SLineages {
                format.needs_taxon_nodes = true;
            }
            if id == FieldId::FullSSeq || id == FieldId::ApproxPIdent {
                format.flags |= OutputFlags::TARGET_SEQS;
            }
            if id == FieldId::QQual || id == FieldId::FullQQual {
                format.store_query_quality = true;
            }
            if id == FieldId::FullQSeqMate {
                format.needs_paired_end_info = true;
            }
            if id == FieldId::ApproxPIdent && score_matrix_name != "BLOSUM62" {
                return Err(
                    "Approximate identity is only supported for the BLOSUM62 scoring matrix."
                        .to_string(),
                );
            }
            if matches!(id, FieldId::FullQSeqMate | FieldId::QSeqTranslated) && !query_translated {
                return Err("Output field only supported for translated search.".to_string());
            }
            format.hsp_values = format.hsp_values | id.hsp_values();
            format.flags |= id.output_flags();
            format.fields.push(id);
        }
        Ok(format)
    }

    pub fn header_format(
        workflow: Workflow,
        output_header: Option<&[&str]>,
    ) -> Result<Header, String> {
        let cluster = matches!(workflow, Workflow::Cluster | Workflow::DeepClust);
        if workflow != Workflow::BlastP && !cluster {
            return Err("header_format".to_string());
        }
        let Some(output_header) = output_header else {
            return Ok(Header::None);
        };
        if output_header.is_empty() {
            return Ok(if cluster {
                Header::Simple
            } else {
                Header::Verbose
            });
        }
        if output_header.len() > 1 {
            return Err(format!(
                "Invalid header format: {}",
                output_header.join(" ")
            ));
        }
        let h = Header::parse(output_header[0])
            .ok_or_else(|| format!("Invalid header format: {}", output_header[0]))?;
        if h == Header::Verbose && cluster {
            return Err("Verbose header format is not supported for cluster workflow.".to_string());
        }
        Ok(h)
    }

    pub fn report_unaligned(&self, report_unaligned: i32) -> bool {
        report_unaligned != 0
            && (report_unaligned != -1 || self.flags.any(OutputFlags::DEFAULT_REPORT_UNALIGNED))
    }
}

pub fn parse_tabular_fields(fields: &[&str]) -> Result<Vec<FieldId>, String> {
    if fields.len() <= 1 {
        return Ok(DEFAULT_TABULAR_FIELDS.to_vec());
    }
    let mut out = Vec::with_capacity(fields.len() - 1);
    for field in &fields[1..] {
        match FieldId::from_name(field) {
            Some(id) => out.push(id),
            None => return Err(format!("Invalid output field: {}", field)),
        }
    }
    Ok(out)
}

/// Matches C++ `TabularFormat::output_header()`.
pub fn output_header(fields: &[FieldId], cluster: bool) -> Result<String, String> {
    let mut out = String::new();
    for (i, field) in fields.iter().enumerate() {
        let def = field
            .output_field()
            .ok_or_else(|| format!("Invalid output field: {:?}", field))?;
        let key = if cluster { def.clust_key } else { def.key };
        if cluster && key.is_empty() {
            return Err(format!(
                "Output field not supported for clustering: {}",
                def.key
            ));
        }
        if i > 0 {
            out.push('\t');
        }
        out.push_str(key);
    }
    out.push('\n');
    Ok(out)
}

/// Matches C++ `TabularFormat::print_header()`.
pub fn print_tabular_header(
    fields: &[FieldId],
    header: Header,
    is_json: bool,
    version: &str,
    invocation: &str,
) -> Result<String, String> {
    let mut out = String::new();
    match header {
        Header::Verbose => {
            out.push_str("# DIAMOND v");
            out.push_str(version);
            out.push_str(". http://github.com/bbuchfink/diamond\n");
            out.push_str("# Invocation: ");
            out.push_str(invocation);
            out.push('\n');
            out.push_str("# Fields: ");
            for (i, field) in fields.iter().enumerate() {
                let def = field
                    .output_field()
                    .ok_or_else(|| format!("Invalid output field: {:?}", field))?;
                if i > 0 {
                    out.push_str(", ");
                }
                out.push_str(def.description);
            }
            out.push('\n');
        }
        Header::Simple => out.push_str(&output_header(fields, false)?),
        Header::None => {}
    }
    if is_json {
        out.push('[');
    }
    Ok(out)
}

/// Matches C++ `TabularFormat::print_footer()`.
pub fn print_tabular_footer(is_json: bool) -> &'static str {
    if is_json {
        "\n]"
    } else {
        ""
    }
}

/// Matches C++ `OutputFormat::print_title(id, full_titles, all_titles, separator)`.
pub fn print_title(
    id: &str,
    full_titles: bool,
    all_titles: bool,
    separator: &str,
    json_array: bool,
) -> String {
    let mut out = String::new();
    let mut n = 0usize;
    let mut rest = id;
    loop {
        // C++ `StringDelimiters::scan` (`util/string/tokenizer.h:89-96`) probes
        // the delimiter array in order and returns the FIRST delimiter that
        // occurs anywhere — i.e. `\x01` always wins over ` >` if both are
        // present, regardless of position. Picking the earliest-position
        // delimiter (the obvious thing) would tokenize a multi-seqid header
        // like `"sp|A >sp|B\x01sp|C"` into 3 pieces instead of C++'s 2.
        let mut split_at = rest.len();
        let mut sep_len = 0usize;
        for sep in FASTA_HEADER_SEP {
            if let Some(i) = rest.find(sep) {
                split_at = i;
                sep_len = sep.len();
                break;
            }
        }
        let s = &rest[..split_at];
        if n > 0 {
            out.push_str(separator);
        }
        if json_array {
            out.push('"');
        }
        if full_titles {
            out.push_str(s);
        } else {
            out.push_str(
                s.split(|c: char| ID_DELIMITERS.contains(c))
                    .next()
                    .unwrap_or(""),
            );
        }
        if json_array {
            out.push('"');
        }
        n += 1;
        if !all_titles {
            break;
        }
        if split_at == rest.len() {
            break;
        }
        rest = &rest[split_at + sep_len..];
    }
    out
}

/// Matches C++ `OutputFormat::print_title(id, full_titles, all_titles, separator)`.
pub fn print_title_xml(
    id: &str,
    full_titles: bool,
    all_titles: bool,
    separator: &str,
    json_array: bool,
) -> String {
    let mut out = String::new();
    let mut n = 0usize;
    let mut rest = id;
    loop {
        // See `print_title` above: C++ scans delimiters in array order, not
        // by earliest position. `\x01` always wins over ` >` when both occur.
        let mut split_at = rest.len();
        let mut sep_len = 0usize;
        for sep in FASTA_HEADER_SEP {
            if let Some(i) = rest.find(sep) {
                split_at = i;
                sep_len = sep.len();
                break;
            }
        }
        let s = &rest[..split_at];
        if n > 0 {
            out.push_str(separator);
        }
        if json_array {
            out.push('"');
        }
        if full_titles {
            out.push_str(&XML.escape_string(s));
        } else {
            out.push_str(
                &XML.escape_string(
                    s.split(|c: char| ID_DELIMITERS.contains(c))
                        .next()
                        .unwrap_or(""),
                ),
            );
        }
        if json_array {
            out.push('"');
        }
        n += 1;
        if !all_titles {
            break;
        }
        if split_at == rest.len() {
            break;
        }
        rest = &rest[split_at + sep_len..];
    }
    out
}

pub fn print_staxids(taxids: &[TaxId], _json: bool) -> String {
    // C++ `TextBuffer::print(vector<T>, char separator)`
    // (`diamond/src/util/text_buffer.h:267-276`) takes a separator argument
    // BUT hard-codes `';'` in the body — so blast_tab_format.cpp:139's
    // `out.print(taxids, json ? ',' : ';')` always emits `;`. Mirror this
    // bug for byte-identical output rather than silently emitting "correct"
    // JSON with commas; downstream tools that parse C++ output have
    // adapted to the `;` form.
    taxids
        .iter()
        .map(|taxid| taxid.to_string())
        .collect::<Vec<_>>()
        .join(";")
}

pub fn print_taxon_names<I>(taxids: I, tree: &TaxonomyTree, json: bool) -> String
where
    I: IntoIterator<Item = TaxId>,
{
    let taxids = taxids.into_iter().collect::<Vec<_>>();
    if taxids.is_empty() {
        return "N/A".to_string();
    }
    let mut out = String::new();
    for (i, taxid) in taxids.iter().enumerate() {
        if i > 0 {
            out.push(if json { ',' } else { ';' });
        }
        if json {
            out.push('"');
        }
        out.push_str(&tree.taxon_scientific_name(*taxid));
        if json {
            out.push('"');
        }
    }
    out
}

pub fn lineage(taxid: TaxId, tree: &TaxonomyTree) -> Vec<TaxId> {
    const MAX_LINEAGE: usize = 64;
    let mut out = Vec::new();
    let mut i = taxid;
    loop {
        if i <= 0 {
            return Vec::new();
        }
        if i == 1 {
            break;
        }
        out.push(i);
        i = tree.parent(i);
        if out.len() >= MAX_LINEAGE {
            panic!("Lineage too long for taxid {}", taxid);
        }
    }
    out
}

pub fn print_lineage(taxids: &[TaxId], tree: &TaxonomyTree, json: bool) -> String {
    if taxids.is_empty() {
        return if json {
            " []".to_string()
        } else {
            "N/A".to_string()
        };
    }
    let mut lineages = std::collections::BTreeSet::new();
    for &taxid in taxids {
        let lineage = lineage(taxid, tree);
        if !lineage.is_empty() {
            lineages.insert(lineage);
        }
    }
    if lineages.is_empty() {
        return if json {
            " []".to_string()
        } else {
            "N/A".to_string()
        };
    }
    let mut out = String::new();
    if json {
        out.push_str(" [\n");
    }
    for (lineage_idx, lin) in lineages.iter().enumerate() {
        if lineage_idx > 0 {
            if json {
                out.push_str(",\n");
            } else {
                out.push_str("<>");
            }
        }
        if json {
            out.push_str("\t\t[");
        }
        for (i, taxid) in lin.iter().rev().enumerate() {
            if i > 0 {
                out.push_str(if json { ", " } else { "; " });
            }
            if json {
                out.push('"');
            }
            out.push_str(&tree.taxon_scientific_name(*taxid));
            if json {
                out.push('"');
            }
        }
        if json {
            out.push(']');
        }
    }
    if json {
        out.push_str("\n\t]");
    }
    out
}

pub fn print_rank_taxon_names(
    taxids: &[TaxId],
    tree: &TaxonomyTree,
    rank: &str,
    json: bool,
) -> String {
    print_taxon_names(tree.rank_taxids(taxids, rank), tree, json)
}

pub fn print_qqual(qual: &str, range: crate::util::interval::Interval) -> String {
    if qual.is_empty() {
        return "*".to_string();
    }
    let begin = range.begin.max(0) as usize;
    let end = range.end.max(range.begin).min(qual.len() as i32) as usize;
    qual[begin..end].to_string()
}

pub fn print_full_qqual(qual: &str) -> &str {
    if qual.is_empty() {
        "*"
    } else {
        qual
    }
}

pub fn print_full_qseq_mate(mate_seq: Option<&[crate::basic::value::Letter]>) -> String {
    let Some(mate_seq) = mate_seq else {
        return "*".to_string();
    };
    let mut out = String::new();
    for &letter in mate_seq {
        out.push(AMINO_ACID_ALPHABET[(letter & LETTER_MASK) as usize] as char);
    }
    out
}

pub fn write_tabular_context_row<W: Write>(
    writer: &mut W,
    r: &HspContext,
    fields: &[FieldId],
    query_translated: bool,
    frame_shift: bool,
    qnum_offset: u64,
    snum_offset: u64,
) -> io::Result<()> {
    for (field_idx, field) in fields.iter().enumerate() {
        if field_idx > 0 {
            write!(writer, "\t")?;
        }
        if r.seed_only() {
            match field {
                FieldId::QSeqId
                | FieldId::QLen
                | FieldId::SSeqId
                | FieldId::SAllSeqId
                | FieldId::SLen
                | FieldId::QStart
                | FieldId::QEnd
                | FieldId::SStart
                | FieldId::SEnd
                | FieldId::Score
                | FieldId::QFrame
                | FieldId::STaxIds
                | FieldId::SSciNames
                | FieldId::SSKingdoms
                | FieldId::STitle
                | FieldId::SAllTitles
                | FieldId::QTitle
                | FieldId::FullSSeq
                | FieldId::QNum
                | FieldId::SNum
                | FieldId::FullQQual
                | FieldId::FullQSeq
                | FieldId::QStrand
                | FieldId::SKingdoms
                | FieldId::SPhylums
                | FieldId::FullQSeqMate
                | FieldId::HspNum => {}
                FieldId::SLineages => {}
                _ => continue,
            }
        }
        match field {
            FieldId::QSeqId | FieldId::QAcc | FieldId::QAccVer => write!(
                writer,
                "{}",
                print_title(&r.query_title, false, false, "", false)
            )?,
            FieldId::QLen => write!(writer, "{}", r.query_len)?,
            FieldId::SSeqId | FieldId::SAcc | FieldId::SAccVer => write!(
                writer,
                "{}",
                print_title(&r.target_title, false, false, "", false)
            )?,
            FieldId::SAllSeqId | FieldId::SAllacc => write!(
                writer,
                "{}",
                print_title(&r.target_title, false, true, ";", false)
            )?,
            FieldId::SLen => write!(writer, "{}", r.subject_len)?,
            FieldId::QStart => write!(writer, "{}", r.oriented_query_range().begin + 1)?,
            FieldId::QEnd => write!(writer, "{}", r.oriented_query_range().end + 1)?,
            FieldId::SStart => write!(writer, "{}", r.subject_source_range().begin + 1)?,
            FieldId::SEnd => write!(writer, "{}", r.subject_source_range().end)?,
            FieldId::QSeq => {
                let range = r.query_source_range();
                let q = r.query_index(r.hsp.frame);
                let begin = range.begin.max(0) as usize;
                let end = range.end.max(range.begin).min(q.len() as i32) as usize;
                for &letter in &q[begin..end] {
                    write!(
                        writer,
                        "{}",
                        AMINO_ACID_ALPHABET[(letter & LETTER_MASK) as usize] as char
                    )?;
                }
            }
            FieldId::SSeq => {
                let mut it = r.begin();
                while it.good() {
                    if it.op() != EditOperation::Insertion {
                        write!(
                            writer,
                            "{}",
                            AMINO_ACID_ALPHABET[(it.subject() & LETTER_MASK) as usize] as char
                        )?;
                    }
                    it.advance();
                }
            }
            FieldId::EValue => write!(writer, "{}", format_evalue(r.evalue()))?,
            FieldId::BitScore => write!(writer, "{}", format_double(r.bit_score()))?,
            FieldId::Score => write!(writer, "{}", r.score())?,
            FieldId::Length => write!(writer, "{}", r.length())?,
            FieldId::PIdent => write!(writer, "{}", format_double(r.id_percent()))?,
            FieldId::NIdent => write!(writer, "{}", r.identities())?,
            FieldId::Mismatch => write!(writer, "{}", r.mismatches())?,
            FieldId::Positive => write!(writer, "{}", r.positives())?,
            FieldId::GapOpen => write!(writer, "{}", r.gap_openings())?,
            FieldId::Gaps => write!(writer, "{}", r.gaps())?,
            FieldId::PPos => write!(
                writer,
                "{}",
                format_double(r.positives() as f64 * 100.0 / r.length() as f64)
            )?,
            FieldId::QFrame | FieldId::QFrameField => {
                write!(writer, "{}", r.blast_query_frame(query_translated))?
            }
            FieldId::BTop => {
                let mut n_matches = 0u32;
                let mut it = r.begin();
                while it.good() {
                    match it.op() {
                        EditOperation::Match => n_matches += 1,
                        EditOperation::Substitution
                        | EditOperation::FrameshiftForward
                        | EditOperation::FrameshiftReverse => {
                            if n_matches > 0 {
                                write!(writer, "{}", n_matches)?;
                                n_matches = 0;
                            }
                            write!(writer, "{}{}", it.query_char(), it.subject_char())?;
                        }
                        EditOperation::Insertion => {
                            if n_matches > 0 {
                                write!(writer, "{}", n_matches)?;
                                n_matches = 0;
                            }
                            write!(writer, "{}-", it.query_char())?;
                        }
                        EditOperation::Deletion => {
                            if n_matches > 0 {
                                write!(writer, "{}", n_matches)?;
                                n_matches = 0;
                            }
                            write!(writer, "-{}", it.subject_char())?;
                        }
                    }
                    it.advance();
                }
                if n_matches > 0 {
                    write!(writer, "{}", n_matches)?;
                }
            }
            FieldId::STitle | FieldId::Stitle0 => write!(
                writer,
                "{}",
                print_title(&r.target_title, true, false, "<>", false)
            )?,
            FieldId::SAllTitles => write!(
                writer,
                "{}",
                print_title(&r.target_title, true, true, "<>", false)
            )?,
            FieldId::QCovHsp | FieldId::QCovS => write!(writer, "{}", format_double(r.qcovhsp()))?,
            FieldId::QTitle => write!(writer, "{}", r.query_title)?,
            FieldId::FullSSeq => {
                for &letter in &r.subject_seq {
                    write!(
                        writer,
                        "{}",
                        AMINO_ACID_ALPHABET[(letter & LETTER_MASK) as usize] as char
                    )?;
                }
            }
            FieldId::QNum => write!(writer, "{}", r.query_oid + qnum_offset)?,
            FieldId::SNum => write!(writer, "{}", r.subject_oid + snum_offset)?,
            FieldId::SCovHsp => write!(writer, "{}", format_double(r.scovhsp()))?,
            FieldId::FullQSeq => {
                for &letter in r.query_index(r.hsp.frame) {
                    write!(
                        writer,
                        "{}",
                        AMINO_ACID_ALPHABET[(letter & LETTER_MASK) as usize] as char
                    )?;
                }
            }
            FieldId::QSeqGapped => {
                let mut it = r.begin();
                while it.good() {
                    write!(writer, "{}", it.query_char())?;
                    it.advance();
                }
            }
            FieldId::SSeqGapped => {
                let mut it = r.begin();
                while it.good() {
                    write!(writer, "{}", it.subject_char())?;
                    it.advance();
                }
            }
            FieldId::QStrand => {
                if query_translated {
                    write!(
                        writer,
                        "{}",
                        if r.blast_query_frame(true) > 0 {
                            '+'
                        } else {
                            '-'
                        }
                    )?;
                } else {
                    write!(writer, "+")?;
                }
            }
            FieldId::Cigar => write!(writer, "{}", crate::output::sam::print_cigar(r))?,
            FieldId::QSeqTranslated => {
                if frame_shift {
                    let mut it = r.begin();
                    while it.good() {
                        if it.op() != EditOperation::Deletion
                            && it.op() != EditOperation::FrameshiftForward
                            && it.op() != EditOperation::FrameshiftReverse
                        {
                            write!(
                                writer,
                                "{}",
                                AMINO_ACID_ALPHABET[(it.query() & LETTER_MASK) as usize] as char
                            )?;
                        }
                        it.advance();
                    }
                } else {
                    let range = r.query_range();
                    let q = r.query_index(r.hsp.frame);
                    let begin = range.begin.max(0) as usize;
                    let end = range.end.max(range.begin).min(q.len() as i32) as usize;
                    for &letter in &q[begin..end] {
                        write!(
                            writer,
                            "{}",
                            AMINO_ACID_ALPHABET[(letter & LETTER_MASK) as usize] as char
                        )?;
                    }
                }
            }
            FieldId::ApproxPIdent => write!(writer, "{}", format_double(r.approx_id()))?,
            // C++ `info.out << r.corrected_bit_score()` (`blast_tab_format.cpp:599`)
            // routes the double through `format_double` (1-decimal under 100,
            // floor at/above 100). Writing the raw f64 diverges for any value ≥ 100.
            FieldId::CorrectedBitScore => {
                write!(writer, "{}", format_double(r.corrected_bit_score()))?
            }
            FieldId::HspNum => write!(writer, "{}", r.hsp_num)?,
            // C++ `info.out.print_d(...)` (`blast_tab_format.cpp:614,618`) uses
            // `snprintf("%lf", x)` — 6-decimal fixed. Rust's default `{}`
            // produces shortest-round-trip and diverges.
            FieldId::NormalizedBitscore => write!(
                writer,
                "{:.6}",
                r.bit_score() / r.query_self_aln_score.max(r.target_self_aln_score)
            )?,
            FieldId::NormalizedNident => write!(
                writer,
                "{:.6}",
                r.identities() as f64
                    / (r.query_index(r.hsp.frame).len() as i32).max(r.subject_len) as f64
            )?,
            FieldId::STaxIds
            | FieldId::SSciNames
            | FieldId::SSKingdoms
            | FieldId::SKingdoms
            | FieldId::SPhylums
            | FieldId::SLineages
            | FieldId::QGi
            | FieldId::SGi
            | FieldId::SStrand
            | FieldId::QQual
            | FieldId::FullQQual
            | FieldId::FullQSeqMate
            | FieldId::NegEValue
            | FieldId::Reserved1
            | FieldId::Reserved2 => write!(writer, "N/A")?,
        }
    }
    writeln!(writer)?;
    Ok(())
}

pub fn write_tabular_context_row_json<W: Write>(
    writer: &mut W,
    r: &HspContext,
    fields: &[FieldId],
    query_translated: bool,
    frame_shift: bool,
    qnum_offset: u64,
    snum_offset: u64,
    subject_taxids: &[TaxId],
    tree: Option<&TaxonomyTree>,
    query_qual: &str,
    mate_seq: Option<&[crate::basic::value::Letter]>,
) -> io::Result<()> {
    if r.hit_num != 0 {
        write!(writer, ",")?;
    }
    writeln!(writer, "\n\t{{")?;
    for (field_idx, field) in fields.iter().enumerate() {
        let output_field = field
            .output_field()
            .expect("tabular output field definition");
        let is_string = output_field.is_string;
        let is_array = output_field.is_array;
        write!(writer, "\t\"{}\":", output_field.key)?;
        if is_string {
            write!(writer, "\"")?;
        }
        if is_array {
            write!(writer, "[")?;
        }

        let seed_only_available = matches!(
            field,
            FieldId::QSeqId
                | FieldId::QLen
                | FieldId::SSeqId
                | FieldId::SAllSeqId
                | FieldId::SLen
                | FieldId::QStart
                | FieldId::QEnd
                | FieldId::SStart
                | FieldId::SEnd
                | FieldId::Score
                | FieldId::QFrame
                | FieldId::STaxIds
                | FieldId::SSciNames
                | FieldId::SSKingdoms
                | FieldId::STitle
                | FieldId::SAllTitles
                | FieldId::QTitle
                | FieldId::FullSSeq
                | FieldId::QNum
                | FieldId::SNum
                | FieldId::FullQQual
                | FieldId::FullQSeq
                | FieldId::QStrand
                | FieldId::SKingdoms
                | FieldId::SPhylums
                | FieldId::FullQSeqMate
                | FieldId::HspNum
                | FieldId::SLineages
        );
        if r.seed_only() && !seed_only_available {
            if !is_string && !is_array {
                write!(writer, "N/A")?;
            }
        } else {
            match field {
                FieldId::SAllSeqId | FieldId::SAllacc => write!(
                    writer,
                    "{}",
                    print_title(&r.target_title, false, true, ",", true)
                )?,
                FieldId::STitle | FieldId::Stitle0 => write!(
                    writer,
                    "{}",
                    print_title(&r.target_title, true, false, ",", false)
                )?,
                FieldId::SAllTitles => write!(
                    writer,
                    "{}",
                    print_title(&r.target_title, true, true, ",", true)
                )?,
                FieldId::STaxIds => write!(writer, "{}", print_staxids(subject_taxids, true))?,
                FieldId::SSciNames => {
                    if let Some(tree) = tree {
                        write!(
                            writer,
                            "{}",
                            print_taxon_names(subject_taxids.iter().copied(), tree, true)
                        )?;
                    }
                }
                FieldId::SSKingdoms => {
                    if let Some(tree) = tree {
                        write!(
                            writer,
                            "{}",
                            print_rank_taxon_names(subject_taxids, tree, "superkingdom", true)
                        )?;
                    }
                }
                FieldId::SKingdoms => {
                    if let Some(tree) = tree {
                        write!(
                            writer,
                            "{}",
                            print_rank_taxon_names(subject_taxids, tree, "kingdom", true)
                        )?;
                    }
                }
                FieldId::SPhylums => {
                    if let Some(tree) = tree {
                        write!(
                            writer,
                            "{}",
                            print_rank_taxon_names(subject_taxids, tree, "phylum", true)
                        )?;
                    }
                }
                FieldId::SLineages => {
                    if let Some(tree) = tree {
                        write!(writer, "{}", print_lineage(subject_taxids, tree, true))?;
                    }
                }
                FieldId::QQual => write!(
                    writer,
                    "{}",
                    print_qqual(query_qual, *r.query_source_range())
                )?,
                FieldId::FullQQual => write!(writer, "{}", print_full_qqual(query_qual))?,
                FieldId::FullQSeqMate => write!(writer, "{}", print_full_qseq_mate(mate_seq))?,
                _ => {
                    let mut cell = Vec::new();
                    write_tabular_context_row(
                        &mut cell,
                        r,
                        std::slice::from_ref(field),
                        query_translated,
                        frame_shift,
                        qnum_offset,
                        snum_offset,
                    )?;
                    let mut cell = String::from_utf8(cell).map_err(|err| {
                        io::Error::new(io::ErrorKind::InvalidData, err.utf8_error())
                    })?;
                    if cell.ends_with('\n') {
                        cell.pop();
                    }
                    write!(writer, "{}", cell)?;
                }
            }
        }

        if is_string {
            write!(writer, "\"")?;
        }
        if is_array {
            write!(writer, "]")?;
        }
        if field_idx < fields.len() - 1 {
            writeln!(writer, ",")?;
        } else {
            writeln!(writer)?;
        }
    }
    write!(writer, "\t}}")?;
    Ok(())
}

pub fn write_tabular_context_row_with_taxonomy<W: Write>(
    writer: &mut W,
    r: &HspContext,
    fields: &[FieldId],
    query_translated: bool,
    frame_shift: bool,
    qnum_offset: u64,
    snum_offset: u64,
    subject_taxids: &[TaxId],
    tree: &TaxonomyTree,
    json: bool,
) -> io::Result<()> {
    for (field_idx, field) in fields.iter().enumerate() {
        if field_idx > 0 {
            write!(writer, "\t")?;
        }
        match field {
            FieldId::STaxIds => write!(writer, "{}", print_staxids(subject_taxids, json))?,
            FieldId::SSciNames => write!(
                writer,
                "{}",
                print_taxon_names(subject_taxids.iter().copied(), tree, json)
            )?,
            FieldId::SSKingdoms => write!(
                writer,
                "{}",
                print_rank_taxon_names(subject_taxids, tree, "superkingdom", json)
            )?,
            FieldId::SKingdoms => write!(
                writer,
                "{}",
                print_rank_taxon_names(subject_taxids, tree, "kingdom", json)
            )?,
            FieldId::SPhylums => write!(
                writer,
                "{}",
                print_rank_taxon_names(subject_taxids, tree, "phylum", json)
            )?,
            FieldId::SLineages => write!(writer, "{}", print_lineage(subject_taxids, tree, json))?,
            _ => {
                let mut cell = Vec::new();
                write_tabular_context_row(
                    &mut cell,
                    r,
                    std::slice::from_ref(field),
                    query_translated,
                    frame_shift,
                    qnum_offset,
                    snum_offset,
                )?;
                let mut cell = String::from_utf8(cell)
                    .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.utf8_error()))?;
                if cell.ends_with('\n') {
                    cell.pop();
                }
                write!(writer, "{}", cell)?;
            }
        }
    }
    writeln!(writer)?;
    Ok(())
}

pub fn write_tabular_context_row_with_query_info<W: Write>(
    writer: &mut W,
    r: &HspContext,
    fields: &[FieldId],
    query_translated: bool,
    frame_shift: bool,
    qnum_offset: u64,
    snum_offset: u64,
    query_qual: &str,
    mate_seq: Option<&[crate::basic::value::Letter]>,
) -> io::Result<()> {
    for (field_idx, field) in fields.iter().enumerate() {
        if field_idx > 0 {
            write!(writer, "\t")?;
        }
        match field {
            FieldId::QQual => write!(
                writer,
                "{}",
                print_qqual(query_qual, *r.query_source_range())
            )?,
            FieldId::FullQQual => write!(writer, "{}", print_full_qqual(query_qual))?,
            FieldId::FullQSeqMate => write!(writer, "{}", print_full_qseq_mate(mate_seq))?,
            _ => {
                let mut cell = Vec::new();
                write_tabular_context_row(
                    &mut cell,
                    r,
                    std::slice::from_ref(field),
                    query_translated,
                    frame_shift,
                    qnum_offset,
                    snum_offset,
                )?;
                let mut cell = String::from_utf8(cell)
                    .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.utf8_error()))?;
                if cell.ends_with('\n') {
                    cell.pop();
                }
                write!(writer, "{}", cell)?;
            }
        }
    }
    writeln!(writer)?;
    Ok(())
}

pub fn write_tabular_query_intro<W: Write>(
    writer: &mut W,
    query_title: &str,
    query_len: i32,
    query_seq: &[crate::basic::value::Letter],
    fields: &[FieldId],
    report_unaligned: bool,
) -> io::Result<()> {
    if !report_unaligned {
        return Ok(());
    }
    for (field_idx, field) in fields.iter().enumerate() {
        if field_idx > 0 {
            write!(writer, "\t")?;
        }
        match field {
            FieldId::QSeqId | FieldId::QAcc | FieldId::QAccVer => write!(
                writer,
                "{}",
                print_title(query_title, false, false, "", false)
            )?,
            FieldId::QLen => write!(writer, "{}", query_len)?,
            FieldId::QFrame | FieldId::QFrameField => write!(writer, "0")?,
            FieldId::QTitle => write!(writer, "{}", query_title)?,
            FieldId::FullQSeq => {
                for &letter in query_seq {
                    write!(
                        writer,
                        "{}",
                        AMINO_ACID_ALPHABET[(letter & LETTER_MASK) as usize] as char
                    )?;
                }
            }
            FieldId::SSeqId
            | FieldId::SAllSeqId
            | FieldId::QSeq
            | FieldId::SSeq
            | FieldId::BTop
            | FieldId::STaxIds
            | FieldId::SSciNames
            | FieldId::SSKingdoms
            | FieldId::STitle
            | FieldId::SAllTitles
            | FieldId::FullSSeq
            | FieldId::QQual
            | FieldId::FullQQual
            | FieldId::QSeqGapped
            | FieldId::SSeqGapped
            | FieldId::QStrand
            | FieldId::Cigar
            | FieldId::SKingdoms
            | FieldId::SPhylums
            | FieldId::FullQSeqMate
            | FieldId::QSeqTranslated
            | FieldId::SLineages => write!(writer, "*")?,
            _ => write!(writer, "-1")?,
        }
    }
    writeln!(writer)?;
    Ok(())
}

pub fn write_tabular_query_intro_with_quality<W: Write>(
    writer: &mut W,
    query_title: &str,
    query_len: i32,
    query_seq: &[crate::basic::value::Letter],
    query_qual: &str,
    fields: &[FieldId],
    report_unaligned: bool,
) -> io::Result<()> {
    if !report_unaligned {
        return Ok(());
    }
    for (field_idx, field) in fields.iter().enumerate() {
        if field_idx > 0 {
            write!(writer, "\t")?;
        }
        match field {
            FieldId::FullQQual => write!(writer, "{}", print_full_qqual(query_qual))?,
            _ => {
                let mut cell = Vec::new();
                write_tabular_query_intro(
                    &mut cell,
                    query_title,
                    query_len,
                    query_seq,
                    std::slice::from_ref(field),
                    true,
                )?;
                let mut cell = String::from_utf8(cell)
                    .map_err(|err| io::Error::new(io::ErrorKind::InvalidData, err.utf8_error()))?;
                if cell.ends_with('\n') {
                    cell.pop();
                }
                write!(writer, "{}", cell)?;
            }
        }
    }
    writeln!(writer)?;
    Ok(())
}

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
            FieldId::QSeqId | FieldId::QAcc | FieldId::QAccVer => {
                write!(writer, "{}", print_title(query_id, false, false, "", false))?
            }
            FieldId::SSeqId | FieldId::SAcc | FieldId::SAccVer => write!(
                writer,
                "{}",
                print_title(subject_id, false, false, "", false)
            )?,
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
                    100.0 * (hsp.subject_range.1 - hsp.subject_range.0) as f64 / subject_len as f64
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
    use crate::align::hsp::{Hsp as AlignHsp, HspContext};
    use crate::basic::packed_transcript::EditOperation;
    use crate::data::taxonomy::{TaxonomyNode, TaxonomyTree};
    use crate::util::interval::Interval;

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
        // C printf %.2e pads exponent to >=2 digits, so e-5 -> e-05.
        assert_eq!(format_evalue(1.5e-5), "1.50e-05");
        assert_eq!(format_evalue(8.35e-9), "8.35e-09");
        assert_eq!(format_evalue(1.0), "1.00e+00");
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
    fn test_output_format_specs() {
        let null = OutputFormatSpec::null_format();
        assert_eq!(null.code, OutputFormatKind::Null);
        assert_eq!(null.hsp_values, HspValues::TRANSCRIPT);
        assert_eq!(null.flags, OutputFlags::NONE);

        let daa = OutputFormatSpec::daa_format(false, false);
        assert_eq!(daa.code, OutputFormatKind::Daa);
        assert!(daa.flags.any(OutputFlags::SSEQID));
        assert!(!daa.flags.any(OutputFlags::ALL_SEQIDS));
        let daa_all = OutputFormatSpec::daa_format(false, true);
        assert!(daa_all.flags.any(OutputFlags::ALL_SEQIDS));
        assert!(!daa_all.flags.any(OutputFlags::FULL_TITLES));
        let daa_titles = OutputFormatSpec::daa_format(true, false);
        assert!(daa_titles.flags.any(OutputFlags::ALL_SEQIDS));
        assert!(daa_titles.flags.any(OutputFlags::FULL_TITLES));

        let paf = OutputFormatSpec::paf_format();
        assert_eq!(paf.code, OutputFormatKind::Paf);
        assert_eq!(paf.hsp_values, HspValues::TRANSCRIPT);
        assert!(paf.flags.any(OutputFlags::SSEQID));
        assert!(paf.report_unaligned(-1));
        assert!(!paf.report_unaligned(0));
        assert!(paf.report_unaligned(1));

        let sam = OutputFormatSpec::sam_format();
        assert_eq!(sam.code, OutputFormatKind::Sam);
        assert_eq!(sam.flags, paf.flags);

        let xml = OutputFormatSpec::xml_format();
        assert_eq!(xml.code, OutputFormatKind::BlastXml);
        assert!(xml.flags.any(OutputFlags::FULL_TITLES));
        assert!(xml.flags.any(OutputFlags::ALL_SEQIDS));

        let pairwise = OutputFormatSpec::pairwise_format();
        assert_eq!(pairwise.code, OutputFormatKind::BlastPairwise);
        assert_eq!(pairwise.flags, xml.flags);

        let taxon = OutputFormatSpec::taxon_format(true);
        assert_eq!(taxon.code, OutputFormatKind::Taxon);
        assert_eq!(taxon.hsp_values, HspValues::NONE);
        assert!(taxon.needs_taxon_id_lists);
        assert!(taxon.needs_taxon_nodes);
        assert!(taxon.needs_taxon_scientific_names);

        let edge = OutputFormatSpec::edge_format();
        assert_eq!(edge.code, OutputFormatKind::Edge);
        assert_eq!(edge.hsp_values, HspValues::COORDS);
        assert!(edge.flags.any(OutputFlags::SSEQID));
    }

    #[test]
    fn test_get_output_format_selection() {
        assert!(matches!(
            get_output_format(&[], "", OutputCommand::BlastP, false, 0, true, "BLOSUM62").unwrap(),
            ResolvedOutputFormat::Tabular(_)
        ));
        assert!(matches!(
            get_output_format(
                &[],
                "out.daa",
                OutputCommand::BlastP,
                false,
                0,
                true,
                "BLOSUM62"
            )
            .unwrap(),
            ResolvedOutputFormat::Spec(OutputFormatSpec {
                code: OutputFormatKind::Daa,
                ..
            })
        ));
        assert!(matches!(
            get_output_format(
                &["json-flat", "qseqid"],
                "",
                OutputCommand::BlastP,
                false,
                0,
                true,
                "BLOSUM62"
            )
            .unwrap(),
            ResolvedOutputFormat::Tabular(TabularFormat { is_json: true, .. })
        ));
        assert!(matches!(
            get_output_format(
                &["edge"],
                "",
                OutputCommand::BlastP,
                false,
                0,
                true,
                "BLOSUM62"
            )
            .unwrap(),
            ResolvedOutputFormat::Spec(OutputFormatSpec {
                code: OutputFormatKind::Edge,
                ..
            })
        ));
        assert!(get_output_format(
            &["bad"],
            "",
            OutputCommand::BlastP,
            false,
            0,
            true,
            "BLOSUM62"
        )
        .is_err());
    }

    #[test]
    fn test_init_output_target_limits_and_toppercent() {
        let base = OutputInitConfig {
            output_format: &["tab"],
            daa_file: "",
            command: OutputCommand::BlastP,
            include_lineage: false,
            frame_shift: 0,
            query_translated: false,
            score_matrix_name: "BLOSUM62",
            ext: "",
            multiprocessing: false,
            global_ranking_targets: 0,
            toppercent: None,
            max_target_seqs: None,
            min_bit_score: 0.0,
            query_range_culling: false,
        };

        let init = init_output(&base).unwrap();
        assert_eq!(init.max_target_seqs, DEFAULT_MAX_TARGET_SEQS);
        assert_eq!(init.toppercent, None);
        assert_eq!(
            init.message,
            "#Target sequences to report alignments for: 25\n"
        );

        let mut cfg = base;
        cfg.max_target_seqs = Some(0);
        let init = init_output(&cfg).unwrap();
        assert_eq!(init.max_target_seqs, i64::MAX);
        assert_eq!(
            init.message,
            "#Target sequences to report alignments for: unlimited\n"
        );

        cfg.toppercent = Some(10.0);
        let init = init_output(&cfg).unwrap();
        assert_eq!(init.toppercent, Some(10.0));
        assert_eq!(
            init.message,
            "Percentage range of top alignment score to report hits: 10\n"
        );

        cfg.toppercent = Some(100.0);
        let init = init_output(&cfg).unwrap();
        assert_eq!(init.toppercent, None);
        assert_eq!(init.max_target_seqs, i64::MAX);

        cfg.toppercent = Some(100.1);
        assert_eq!(
            init_output(&cfg).unwrap_err(),
            "Allowed value range for --top is between 0.0 and 100.0"
        );
    }

    #[test]
    fn test_init_output_validation_and_dp_fields() {
        let base = OutputInitConfig {
            output_format: &["daa"],
            daa_file: "",
            command: OutputCommand::BlastP,
            include_lineage: false,
            frame_shift: 0,
            query_translated: false,
            score_matrix_name: "BLOSUM62",
            ext: "",
            multiprocessing: false,
            global_ranking_targets: 0,
            toppercent: None,
            max_target_seqs: None,
            min_bit_score: 0.0,
            query_range_culling: false,
        };

        let mut cfg = base;
        cfg.ext = "none";
        assert_eq!(
            init_output(&cfg).unwrap_err(),
            "--ext none is not supported for DAA output."
        );

        cfg = base;
        cfg.multiprocessing = true;
        assert_eq!(
            init_output(&cfg).unwrap_err(),
            "The DAA format is not supported in multiprocessing mode."
        );

        cfg = base;
        cfg.global_ranking_targets = 1;
        assert_eq!(
            init_output(&cfg).unwrap_err(),
            "The DAA format is not supported in global ranking mode."
        );

        cfg = base;
        cfg.output_format = &["102"];
        cfg.command = OutputCommand::View;
        assert_eq!(
            init_output(&cfg).unwrap_err(),
            "Taxonomy features are not supported for the DAA format."
        );

        cfg.command = OutputCommand::BlastP;
        let init = init_output(&cfg).unwrap();
        assert_eq!(init.toppercent, Some(10.0));

        cfg.min_bit_score = 1.0;
        let init = init_output(&cfg).unwrap();
        assert_eq!(init.toppercent, None);

        cfg.output_format = &["tab"];
        cfg.frame_shift = 15;
        cfg.ext = "none";
        assert_eq!(
            init_output(&cfg).unwrap_err(),
            "Frameshift alignment does not support --ext none."
        );

        cfg.ext = "";
        let init = init_output(&cfg).unwrap();
        assert_eq!(init.dp_fields, HspValues::TRANSCRIPT);
    }

    #[test]
    fn test_default_fields() {
        assert_eq!(DEFAULT_TABULAR_FIELDS.len(), 12);
    }

    #[test]
    fn test_tabular_field_metadata_and_parse() {
        let fields = parse_tabular_fields(&["tab", "qseqid", "sseqid", "bitscore"]).unwrap();
        assert_eq!(
            fields,
            vec![FieldId::QSeqId, FieldId::SSeqId, FieldId::BitScore]
        );
        assert_eq!(
            parse_tabular_fields(&["tab"]).unwrap(),
            DEFAULT_TABULAR_FIELDS
        );
        assert_eq!(
            parse_tabular_fields(&["tab", "bad"]).unwrap_err(),
            "Invalid output field: bad"
        );
        let def = FieldId::SAllSeqId.output_field().unwrap();
        assert_eq!(def.key, "sallseqid");
        assert!(def.is_array);
        assert!(!def.is_string);
        assert_eq!(
            FieldId::PIdent.hsp_values(),
            HspValues::IDENT | HspValues::LENGTH
        );
        assert!(FieldId::SAllSeqId
            .output_flags()
            .any(OutputFlags::ALL_SEQIDS));
        assert!(FieldId::SAllSeqId.output_flags().any(OutputFlags::SSEQID));
    }

    #[test]
    fn test_tabular_format_constructor_default_and_custom() {
        let default = TabularFormat::new(&["tab"], false, 0, false, "BLOSUM62").unwrap();
        assert_eq!(default.fields, DEFAULT_TABULAR_FIELDS);
        assert_eq!(
            default.hsp_values,
            HspValues::QUERY_COORDS
                | HspValues::TARGET_COORDS
                | HspValues::LENGTH
                | HspValues::IDENT
                | HspValues::MISMATCHES
                | HspValues::GAP_OPENINGS
        );
        assert!(default.flags.any(OutputFlags::SSEQID));
        assert!(!default.report_unaligned(-1));

        let fs = TabularFormat::new(&["tab"], false, 15, false, "BLOSUM62").unwrap();
        assert_eq!(fs.hsp_values, HspValues::TRANSCRIPT);

        let custom = TabularFormat::new(
            &[
                "tab",
                "sallseqid",
                "sscinames",
                "slineages",
                "full_sseq",
                "qqual",
            ],
            false,
            0,
            false,
            "BLOSUM62",
        )
        .unwrap();
        assert_eq!(
            custom.fields,
            vec![
                FieldId::SAllSeqId,
                FieldId::SSciNames,
                FieldId::SLineages,
                FieldId::FullSSeq,
                FieldId::QQual,
            ]
        );
        assert!(custom.needs_taxon_id_lists);
        assert!(custom.needs_taxon_scientific_names);
        assert!(custom.needs_taxon_nodes);
        assert!(custom.flags.any(OutputFlags::ALL_SEQIDS));
        assert!(custom.flags.any(OutputFlags::TARGET_SEQS));
        assert!(custom.store_query_quality);

        let json =
            TabularFormat::new(&["json-flat", "qseqid"], true, 0, false, "BLOSUM62").unwrap();
        assert_eq!(json.code, FormatCode::JsonFlat);
        assert_eq!(json.query_separator, ',');
        assert!(json.is_json);
    }

    #[test]
    fn test_tabular_format_constructor_validation_and_header_format() {
        assert_eq!(
            TabularFormat::new(&["tab", "bad"], false, 0, false, "BLOSUM62").unwrap_err(),
            "Invalid output field: bad"
        );
        assert_eq!(
            TabularFormat::new(&["tab", "approx_pident"], false, 0, false, "PAM30").unwrap_err(),
            "Approximate identity is only supported for the BLOSUM62 scoring matrix."
        );
        assert_eq!(
            TabularFormat::new(&["tab", "qseq_translated"], false, 0, false, "BLOSUM62")
                .unwrap_err(),
            "Output field only supported for translated search."
        );
        let translated =
            TabularFormat::new(&["tab", "qseq_translated"], false, 0, true, "BLOSUM62").unwrap();
        assert_eq!(translated.hsp_values, HspValues::TRANSCRIPT);

        assert_eq!(
            TabularFormat::header_format(Workflow::BlastP, None).unwrap(),
            Header::None
        );
        assert_eq!(
            TabularFormat::header_format(Workflow::BlastP, Some(&[])).unwrap(),
            Header::Verbose
        );
        assert_eq!(
            TabularFormat::header_format(Workflow::Cluster, Some(&[])).unwrap(),
            Header::Simple
        );
        assert_eq!(
            TabularFormat::header_format(Workflow::Cluster, Some(&["verbose"])).unwrap_err(),
            "Verbose header format is not supported for cluster workflow."
        );
        assert_eq!(
            TabularFormat::header_format(Workflow::BlastP, Some(&["simple", "verbose"]))
                .unwrap_err(),
            "Invalid header format: simple verbose"
        );
    }

    #[test]
    fn test_tabular_headers_and_footer() {
        let fields = [FieldId::QSeqId, FieldId::SSeqId, FieldId::BitScore];
        assert_eq!(
            output_header(&fields, false).unwrap(),
            "qseqid\tsseqid\tbitscore\n"
        );
        assert_eq!(
            output_header(&fields, true).unwrap(),
            "cseqid\tmseqid\tBitscore\n"
        );
        assert_eq!(
            output_header(&[FieldId::BTop], true).unwrap_err(),
            "Output field not supported for clustering: btop"
        );
        assert_eq!(Header::parse("simple"), Some(Header::Simple));
        assert_eq!(Header::parse("verbose"), Some(Header::Verbose));
        assert_eq!(Header::parse("0"), Some(Header::None));
        assert_eq!(print_tabular_footer(false), "");
        assert_eq!(print_tabular_footer(true), "\n]");
    }

    #[test]
    fn test_print_tabular_header() {
        let fields = [FieldId::QSeqId, FieldId::BitScore];
        assert_eq!(
            print_tabular_header(&fields, Header::Simple, false, "2.x", "diamond blastp").unwrap(),
            "qseqid\tbitscore\n"
        );
        assert_eq!(
            print_tabular_header(&fields, Header::None, true, "2.x", "diamond blastp").unwrap(),
            "["
        );
        let verbose =
            print_tabular_header(&fields, Header::Verbose, true, "2.x", "diamond blastp").unwrap();
        assert!(verbose.starts_with("# DIAMOND v2.x. http://github.com/bbuchfink/diamond\n"));
        assert!(verbose.contains("# Invocation: diamond blastp\n"));
        assert!(verbose.contains("# Fields: Query Seq - id, Bit score\n["));
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

    #[test]
    fn test_print_title_matches_diamond_token_rules() {
        assert_eq!(
            print_title("sp|A desc\x01tr|B other", false, false, "", false),
            "sp|A"
        );
        assert_eq!(
            print_title("sp|A desc\x01tr|B other", false, true, ";", false),
            "sp|A;tr|B"
        );
        assert_eq!(
            print_title("sp|A desc\x01tr|B other", true, true, "<>", false),
            "sp|A desc<>tr|B other"
        );
        assert_eq!(
            print_title("sp|A desc\x01tr|B other", false, true, ",", true),
            "\"sp|A\",\"tr|B\""
        );
        assert_eq!(
            print_title("sp|A desc >tr|B other", false, true, ";", false),
            "sp|A;tr|B"
        );
        assert_eq!(
            print_title("sp|A desc >tr|B other", true, true, "<>", false),
            "sp|A desc<>tr|B other"
        );
        assert_eq!(
            print_title_xml("sp|A <desc> >tr|B & other", true, true, " &gt;", false),
            "sp|A &lt;desc&gt; &gt;tr|B &amp; other"
        );
    }

    fn tabular_taxonomy_tree() -> TaxonomyTree {
        let mut tree = TaxonomyTree::new();
        tree.add_node(TaxonomyNode {
            taxid: 1,
            parent: 1,
            rank: "root".into(),
            name: "root".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 2,
            parent: 1,
            rank: "superkingdom".into(),
            name: "Bacteria".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 3,
            parent: 1,
            rank: "superkingdom".into(),
            name: "Archaea".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 10,
            parent: 2,
            rank: "kingdom".into(),
            name: "Pseudomonadati".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 20,
            parent: 10,
            rank: "phylum".into(),
            name: "Proteobacteria".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 30,
            parent: 3,
            rank: "phylum".into(),
            name: "Euryarchaeota".into(),
        });
        tree
    }

    #[test]
    fn test_tabular_taxonomy_print_helpers() {
        let tree = tabular_taxonomy_tree();
        assert_eq!(print_staxids(&[20, 30], false), "20;30");
        // C++ `TextBuffer::print` hard-codes `;` even for JSON, see comment
        // on `print_staxids`. Mirror this here.
        assert_eq!(print_staxids(&[20, 30], true), "20;30");
        assert_eq!(
            print_taxon_names([20, 30], &tree, false),
            "Proteobacteria;Euryarchaeota"
        );
        assert_eq!(
            print_taxon_names([20, 30], &tree, true),
            "\"Proteobacteria\",\"Euryarchaeota\""
        );
        assert_eq!(print_taxon_names([], &tree, false), "N/A");
        assert_eq!(lineage(20, &tree), vec![20, 10, 2]);
        assert_eq!(
            print_lineage(&[20, 30], &tree, false),
            "Bacteria; Pseudomonadati; Proteobacteria<>Archaea; Euryarchaeota"
        );
        assert_eq!(print_lineage(&[], &tree, false), "N/A");
        assert_eq!(print_lineage(&[], &tree, true), " []");
        assert_eq!(
            print_rank_taxon_names(&[20, 30], &tree, "superkingdom", false),
            "Bacteria;Archaea"
        );
        assert_eq!(
            print_rank_taxon_names(&[20, 30], &tree, "phylum", false),
            "Proteobacteria;Euryarchaeota"
        );
    }

    #[test]
    fn test_tabular_query_quality_and_mate_helpers() {
        assert_eq!(print_qqual("", Interval::new(1, 3)), "*");
        assert_eq!(print_qqual("abcdef", Interval::new(1, 4)), "bcd");
        assert_eq!(print_qqual("abcdef", Interval::new(-2, 99)), "abcdef");
        assert_eq!(print_full_qqual(""), "*");
        assert_eq!(print_full_qqual("abcdef"), "abcdef");
        assert_eq!(print_full_qseq_mate(None), "*");
        assert_eq!(print_full_qseq_mate(Some(&[0, 1, 2, 3])), "ARND");
    }

    #[test]
    fn test_write_tabular_context_row_from_transcript() {
        let mut hsp = AlignHsp::new();
        hsp.score = 42;
        hsp.evalue = 2.0e-5;
        hsp.bit_score = 21.5;
        hsp.corrected_bit_score = 20.25;
        hsp.approx_id = 55.5;
        hsp.frame = 0;
        hsp.length = 5;
        hsp.identities = 2;
        hsp.mismatches = 1;
        hsp.positives = 2;
        hsp.gap_openings = 2;
        hsp.gaps = 2;
        hsp.query_range = Interval::new(0, 4);
        hsp.query_source_range = Interval::new(0, 4);
        hsp.subject_range = Interval::new(0, 4);
        hsp.subject_source_range = Interval::new(0, 4);
        hsp.transcript.push_with_count(EditOperation::Match, 1);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 1);
        hsp.transcript.push_with_count(EditOperation::Insertion, 1);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 6);
        hsp.transcript.push_with_count(EditOperation::Match, 1);

        let ctx = HspContext::new(
            hsp,
            0,
            10,
            vec![vec![0, 4, 3, 6, 13]],
            5,
            "query one",
            20,
            6,
            "subj one\x01subj two",
            0,
            3,
            vec![0, 1, 6, 6, 13, 14],
            100.0,
            80.0,
        );
        let fields = [
            FieldId::QSeqId,
            FieldId::SSeqId,
            FieldId::SAllSeqId,
            FieldId::QStart,
            FieldId::QEnd,
            FieldId::SStart,
            FieldId::SEnd,
            FieldId::SSeq,
            FieldId::BTop,
            FieldId::QSeqGapped,
            FieldId::SSeqGapped,
            FieldId::Cigar,
            FieldId::QCovHsp,
            FieldId::SCovHsp,
            FieldId::QStrand,
            FieldId::QSeqTranslated,
            FieldId::HspNum,
            FieldId::NormalizedBitscore,
            FieldId::NormalizedNident,
            FieldId::CorrectedBitScore,
        ];
        let mut buf = Vec::new();
        write_tabular_context_row(&mut buf, &ctx, &fields, false, false, 1000, 2000).unwrap();
        // Match C++ field formatting:
        //   NormalizedBitscore  → `print_d` = `snprintf("%lf", x)` = 6-decimal fixed (0.215 → "0.215000")
        //   NormalizedNident    → `print_d` = "%lf" = 6-decimal fixed (1/3 → "0.333333")
        //   CorrectedBitScore   → routed through `format_double` (1-decimal under 100,
        //                          floor at/above 100). 20.25 rounds to 20.3 (half-away-from-zero).
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "query\tsubj\tsubj;subj\t1\t4\t1\t4\tAREE\t1CRD--E1\tACD-E\tAR-EE\t2M1I1D1M\t80.0\t66.7\t+\tACDE\t3\t0.215000\t0.333333\t20.3\n"
        );
    }

    #[test]
    fn test_write_tabular_context_row_with_taxonomy() {
        let tree = tabular_taxonomy_tree();
        let mut hsp = AlignHsp::new();
        hsp.score = 7;
        hsp.evalue = 1.0;
        hsp.query_source_range = Interval::new(0, 3);
        hsp.subject_source_range = Interval::new(0, 3);
        hsp.subject_range = Interval::new(0, 3);
        let ctx = HspContext::new(
            hsp,
            0,
            10,
            vec![vec![0, 1, 2]],
            3,
            "query one",
            20,
            3,
            "subj one",
            0,
            0,
            vec![0, 1, 2],
            0.0,
            0.0,
        );
        let fields = [
            FieldId::QSeqId,
            FieldId::STaxIds,
            FieldId::SSciNames,
            FieldId::SSKingdoms,
            FieldId::SKingdoms,
            FieldId::SPhylums,
            FieldId::SLineages,
        ];
        let mut buf = Vec::new();
        write_tabular_context_row_with_taxonomy(
            &mut buf,
            &ctx,
            &fields,
            false,
            false,
            0,
            0,
            &[20, 30],
            &tree,
            false,
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "query\t20;30\tProteobacteria;Euryarchaeota\tBacteria;Archaea\tPseudomonadati\tProteobacteria;Euryarchaeota\tBacteria; Pseudomonadati; Proteobacteria<>Archaea; Euryarchaeota\n"
        );
    }

    #[test]
    fn test_write_tabular_context_row_with_query_info() {
        let mut hsp = AlignHsp::new();
        hsp.score = 9;
        hsp.query_source_range = Interval::new(1, 4);
        hsp.subject_range = Interval::new(0, 3);
        hsp.subject_source_range = Interval::new(0, 3);
        let ctx = HspContext::new(
            hsp,
            0,
            10,
            vec![vec![0, 1, 2, 3, 4]],
            5,
            "query one",
            20,
            3,
            "subj one",
            0,
            0,
            vec![0, 1, 2],
            0.0,
            0.0,
        );
        let fields = [
            FieldId::QSeqId,
            FieldId::QQual,
            FieldId::FullQQual,
            FieldId::FullQSeqMate,
            FieldId::Score,
        ];
        let mut buf = Vec::new();
        write_tabular_context_row_with_query_info(
            &mut buf,
            &ctx,
            &fields,
            false,
            false,
            0,
            0,
            "abcdef",
            Some(&[0, 1, 2]),
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "query\tbcd\tabcdef\tARN\t9\n"
        );

        let mut empty = Vec::new();
        write_tabular_context_row_with_query_info(
            &mut empty,
            &ctx,
            &[FieldId::QQual, FieldId::FullQQual, FieldId::FullQSeqMate],
            false,
            false,
            0,
            0,
            "",
            None,
        )
        .unwrap();
        assert_eq!(String::from_utf8(empty).unwrap(), "*\t*\t*\n");
    }

    #[test]
    fn test_write_tabular_context_row_json() {
        let tree = tabular_taxonomy_tree();
        let mut hsp = AlignHsp::new();
        hsp.score = 9;
        hsp.query_source_range = Interval::new(1, 4);
        hsp.subject_range = Interval::new(0, 3);
        hsp.subject_source_range = Interval::new(0, 3);
        let ctx = HspContext::new(
            hsp,
            0,
            10,
            vec![vec![0, 1, 2, 3, 4]],
            5,
            "query one",
            20,
            3,
            "subj one\x01subj two",
            0,
            0,
            vec![0, 1, 2],
            0.0,
            0.0,
        );
        let fields = [
            FieldId::QSeqId,
            FieldId::SAllSeqId,
            FieldId::STaxIds,
            FieldId::SSciNames,
            FieldId::STitle,
            FieldId::QQual,
            FieldId::Score,
        ];
        let mut buf = Vec::new();
        write_tabular_context_row_json(
            &mut buf,
            &ctx,
            &fields,
            false,
            false,
            0,
            0,
            &[20],
            Some(&tree),
            "abcdef",
            None,
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "\n\t{\n\t\"qseqid\":\"query\",\n\t\"sallseqid\":[\"subj\",\"subj\"],\n\t\"staxids\":[20],\n\t\"sscinames\":[\"Proteobacteria\"],\n\t\"stitle\":\"subj one\",\n\t\"qqual\":\"bcd\",\n\t\"score\":9\n\t}"
        );
    }

    #[test]
    fn test_write_tabular_context_row_json_seed_only_blank() {
        let mut hsp = AlignHsp::new();
        hsp.seed_only = true;
        hsp.score = 9;
        hsp.query_source_range = Interval::new(0, 3);
        hsp.subject_source_range = Interval::new(0, 3);
        let ctx = HspContext::new(
            hsp,
            0,
            10,
            vec![vec![0, 1, 2]],
            3,
            "query one",
            20,
            3,
            "subj one",
            1,
            0,
            vec![0, 1, 2],
            0.0,
            0.0,
        );
        let fields = [FieldId::QSeqId, FieldId::PIdent, FieldId::QSeq];
        let mut buf = Vec::new();
        write_tabular_context_row_json(
            &mut buf,
            &ctx,
            &fields,
            false,
            false,
            0,
            0,
            &[],
            None,
            "",
            None,
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            ",\n\t{\n\t\"qseqid\":\"query\",\n\t\"pident\":N/A,\n\t\"qseq\":\"\"\n\t}"
        );
    }

    #[test]
    fn test_write_tabular_query_intro() {
        let fields = [
            FieldId::QSeqId,
            FieldId::QLen,
            FieldId::SSeqId,
            FieldId::QFrame,
            FieldId::FullQSeq,
            FieldId::EValue,
        ];
        let mut buf = Vec::new();
        write_tabular_query_intro(&mut buf, "query one", 4, &[0, 1, 2, 3], &fields, true).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "query\t4\t*\t0\tARND\t-1\n"
        );
    }

    #[test]
    fn test_write_tabular_query_intro_with_quality() {
        let fields = [
            FieldId::QSeqId,
            FieldId::FullQQual,
            FieldId::FullQSeq,
            FieldId::EValue,
        ];
        let mut buf = Vec::new();
        write_tabular_query_intro_with_quality(
            &mut buf,
            "query one",
            4,
            &[0, 1, 2, 3],
            "abcd",
            &fields,
            true,
        )
        .unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "query\tabcd\tARND\t-1\n");
    }
}
