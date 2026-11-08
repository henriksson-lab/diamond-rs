//! Concrete native search adapters for the translated transcript renderers.
use std::io::{self, Write};

use crate::align::hsp::{ByteHspContext as HspContext, Match};
use crate::basic::value::Letter;
use crate::data::block::Block;
use crate::output::blast_pairwise_format::PairwiseFormat;
use crate::output::paf_format::PafFormat;
use crate::output::xml_format::XmlFormat;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::sequence::AccessionParsing;

/// Actual reference metadata consumed by native rendering. Block merging only
/// retains matched reference titles/lengths, without fabricating sequence data.
pub(crate) trait SearchSubjects {
    fn title_bytes(&self, id: usize) -> io::Result<&[u8]>;
    fn title(&self, id: usize) -> io::Result<&str> {
        std::str::from_utf8(self.title_bytes(id)?)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
    }
    fn length(&self, id: usize) -> usize;
    fn original_oid(&self, id: usize) -> u64;
    fn sequence_count(&self) -> u64;
    fn self_aln_score(&self, _id: usize) -> Option<f64> {
        None
    }
    fn sequence(&self, _id: usize) -> Option<&[Letter]> {
        None
    }
    fn taxids(&self, _id: usize) -> &[crate::basic::value::TaxId] {
        &[]
    }
    fn taxonomy(&self) -> Option<&crate::data::taxonomy::TaxonomyTree> {
        None
    }
}

impl SearchSubjects for Block {
    fn title_bytes(&self, id: usize) -> io::Result<&[u8]> {
        Ok(self.ids().map_err(io::Error::other)?.get(id))
    }
    fn length(&self, id: usize) -> usize {
        self.seqs().length(id) as usize
    }
    fn original_oid(&self, id: usize) -> u64 {
        self.block_id2oid(id as u32)
    }
    fn sequence_count(&self) -> u64 {
        self.seqs().len() as u64
    }
    fn sequence(&self, id: usize) -> Option<&[Letter]> {
        Some(self.seqs().get(id))
    }
}

pub(crate) fn tabular_format(
    config: &super::blastp::BlastpConfig,
) -> io::Result<crate::output::blast_tab_format::ByteBlastTabFormat> {
    let mut output_format = config.outfmt.clone();
    if output_format.is_empty() {
        output_format.push(b"6".to_vec());
    }
    let json = matches!(output_format[0].as_slice(), b"104" | b"json-flat");
    output_format[0] = if json { b"104".to_vec() } else { b"6".to_vec() };
    let mut formatter = crate::output::blast_tab_format::ByteBlastTabFormat::new_bytes(
        crate::output::blast_tab_format::ByteBlastTabConfig {
            output_format,
            json,
            frame_shift: config.frame_shift,
            query_translated: config.translated_query_layout.is_some(),
            score_matrix_name: if config.tuning.matrix_file.as_os_str().is_empty() {
                config.matrix.to_ascii_uppercase()
            } else {
                b"custom".to_vec()
            },
            report_unaligned: config.report_unaligned,
            qnum_offset: config.tuning.output.qnum_offset as u64,
            snum_offset: config.tuning.output.snum_offset as u64,
            output_header: None,
            version: crate::basic::consts::VERSION_STRING.into(),
            invocation: config.sam.invocation.clone(),
            workflow: crate::output::format::Workflow::BlastP,
            ..Default::default()
        },
    )
    .map_err(io::Error::other)?;
    // Source validates the header at output start, after creating the sink.
    formatter.config.output_header = config.output_header.clone();
    Ok(formatter)
}

pub(crate) fn self_alignment_scores(
    sequences: &crate::data::sequence_set::SequenceSet,
    matrix: &ScoreMatrix,
    cbs: crate::stats::cbs::CbsMode,
) -> Vec<f64> {
    self_alignment_scores_with_window(sequences, matrix, cbs, 40)
}

pub(crate) fn self_alignment_scores_with_window(
    sequences: &crate::data::sequence_set::SequenceSet,
    matrix: &ScoreMatrix,
    cbs: crate::stats::cbs::CbsMode,
    window: i32,
) -> Vec<f64> {
    use rayon::prelude::*;
    (0..sequences.len())
        .into_par_iter()
        .map(|id| {
            let sequence = sequences.get(id);
            let bias = cbs.hauser().then(|| {
                crate::stats::cbs::hauser_correction_window(
                    sequence,
                    matrix,
                    window as u32 as usize,
                )
            });
            matrix.bitscore(crate::dp::ungapped::self_score_with_cbs(
                sequence,
                bias.as_deref(),
                matrix,
            ) as f64)
        })
        .collect()
}

/// Concrete Output::Info/HspContext adapter for the complete tabular callbacks.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_tabular_query<W: Write, S: SearchSubjects + ?Sized>(
    out: &mut W,
    format: &crate::output::blast_tab_format::ByteBlastTabFormat,
    query_id: u32,
    query_oid: u64,
    title: &[u8],
    source: &[Letter],
    frames: &[&[Letter]],
    mate: Option<&[Letter]>,
    quality: &[u8],
    query_self_aln_score: Option<f64>,
    matches: Vec<Match>,
    database: &S,
) -> io::Result<u64> {
    let unaligned = matches.iter().all(|target| target.hsps.is_empty());
    format.write_query_intro_bytes(out, title, source.len() as i32, source, quality, unaligned)?;
    let needs_self = format
        .format
        .fields
        .contains(&crate::output::format::FieldId::NormalizedBitscore);
    let mut context = HspContext {
        sequence_type: if format.config.nucleotide_query {
            crate::basic::value::SequenceType::Nucleotide
        } else {
            crate::basic::value::SequenceType::AminoAcid
        },
        query: frames.iter().map(|frame| frame.to_vec()).collect(),
        query_source: format
            .format
            .fields
            .iter()
            .any(|field| {
                matches!(
                    field,
                    crate::output::format::FieldId::QSeq | crate::output::format::FieldId::FullQSeq
                )
            })
            .then(|| source.to_vec()),
        query_source_len: source.len() as i32,
        query_len: source.len() as i32,
        query_title: crate::util::sequence::c_string_bytes(title).to_vec(),
        query_oid,
        query_id,
        query_self_aln_score: if needs_self {
            query_self_aln_score
                .ok_or_else(|| io::Error::other("Missing query self-alignment score"))?
        } else {
            0.0
        },
        ..Default::default()
    };
    let mut rows = 0;
    for (hit_num, target) in matches.into_iter().enumerate() {
        let id = target.target_block_id as usize;
        context.target_title =
            crate::util::sequence::c_string_bytes(database.title_bytes(id)?).to_vec();
        context.subject_oid = database.original_oid(id);
        context.subject_len = database.length(id) as i32;
        context.subject_seq = database.sequence(id).unwrap_or(&[]).to_vec();
        context.target_self_aln_score = if needs_self {
            database
                .self_aln_score(id)
                .ok_or_else(|| io::Error::other("Missing reference self-alignment score"))?
        } else {
            0.0
        };
        context.hit_num = hit_num as u32;
        for (hsp_num, hsp) in target.hsps.into_iter().enumerate() {
            context.hsp = hsp;
            if format
                .format
                .fields
                .contains(&crate::output::format::FieldId::ApproxPIdent)
            {
                context.hsp.approx_id = context.hsp.approx_id_percent(
                    context.query_index(context.hsp.frame),
                    &context.subject_seq,
                );
            }
            context.hsp_num = hsp_num as u32;
            format.write_match_bytes(
                out,
                &context,
                database.taxids(id),
                database.taxonomy(),
                quality,
                mate,
            )?;
            rows += 1;
        }
    }
    Ok(rows)
}

/// DAA reference dictionary title policy. `all_subject_titles` takes
/// precedence over `all_subject_ids`, matching the source DAA format flags.
#[derive(Clone, Debug, Default)]
pub struct ArchiveOptions {
    pub all_subject_titles: bool,
    pub all_subject_ids: bool,
}

impl ArchiveOptions {
    pub(crate) fn dictionary_title_bytes(&self, title: &[u8]) -> Vec<u8> {
        let flags = crate::output::format::OutputFormatSpec::daa_format(
            self.all_subject_titles,
            self.all_subject_ids,
        )
        .flags;
        if flags.any(crate::output::format::OutputFlags::FULL_TITLES) {
            // Block::dict_id copies raw loaded titles here; print_title would
            // tokenize and rewrite alternate header separators.
            title.to_owned()
        } else if flags.any(crate::output::format::OutputFlags::ALL_SEQIDS) {
            crate::util::sequence::all_seqids_bytes(title)
        } else {
            crate::util::sequence::seqid_bytes(title)
        }
    }
}

/// Search-specific SAM settings, retained across translated and blocked adapters.
#[derive(Clone, Debug, Default)]
pub struct SamOptions {
    pub invocation: Vec<u8>,
    pub all_subject_titles: bool,
    pub query_len: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SearchFormat {
    Pairwise,
    Xml,
    Paf,
    Sam,
    Null,
    Edge,
    Bin1,
    Clustering,
    Taxon { include_lineage: bool },
}

impl SearchFormat {
    pub(crate) fn from_config(config: &super::blastp::BlastpConfig) -> Option<Self> {
        let format = Self::from_name_bytes(config.outfmt.first()?)?;
        Some(match format {
            Self::Taxon { .. } => Self::Taxon {
                include_lineage: config.include_lineage,
            },
            format => format,
        })
    }

    pub(crate) fn hsp_values(self) -> crate::dp::swipe::HspValues {
        match self {
            Self::Taxon { include_lineage } => {
                crate::output::taxon_format::TaxonFormat::new(include_lineage)
                    .output_format_spec()
                    .hsp_values
            }
            Self::Edge => crate::output::format::OutputFormatSpec::edge_format().hsp_values,
            Self::Bin1 => crate::mcl::clustering_format::Bin1Format::HSP_VALUES,
            Self::Clustering => crate::dp::swipe::HspValues::NONE,
            Self::Null => {
                crate::output::output_format::NullFormat
                    .output_format_spec()
                    .hsp_values
            }
            _ => crate::dp::swipe::HspValues::TRANSCRIPT,
        }
    }

    fn report_unaligned(self, setting: i32) -> bool {
        setting != 0
            && (setting != -1
                || !matches!(
                    self,
                    Self::Null | Self::Edge | Self::Bin1 | Self::Clustering
                ))
    }

    pub(crate) fn from_name(name: &str) -> Option<Self> {
        Self::from_name_bytes(name.as_bytes())
    }
    pub(crate) fn from_name_bytes(name: &[u8]) -> Option<Self> {
        match name {
            b"0" => Some(Self::Pairwise),
            b"5" | b"xml" => Some(Self::Xml),
            b"103" | b"paf" => Some(Self::Paf),
            b"101" | b"sam" => Some(Self::Sam),
            b"null" => Some(Self::Null),
            b"edge" => Some(Self::Edge),
            b"bin1" if cfg!(feature = b"mcl") => Some(Self::Bin1),
            b"clus" if cfg!(feature = b"mcl") => Some(Self::Clustering),
            b"102" => Some(Self::Taxon {
                include_lineage: false,
            }),
            _ => None,
        }
    }

    pub(crate) fn header<W: Write>(
        self,
        out: &mut W,
        config: &super::blastp::BlastpConfig,
        matrix: &ScoreMatrix,
        first_title: &[u8],
        first_len: i32,
    ) -> io::Result<()> {
        self.header_with_query_mode(
            out,
            config,
            matrix,
            first_title,
            first_len,
            if config.translated_query_layout.is_some() {
                3
            } else {
                2
            },
        )
    }

    pub(crate) fn header_with_query_mode<W: Write>(
        self,
        out: &mut W,
        config: &super::blastp::BlastpConfig,
        matrix: &ScoreMatrix,
        first_title: &[u8],
        first_len: i32,
        query_mode: i32,
    ) -> io::Result<()> {
        match self {
            Self::Pairwise => PairwiseFormat.print_header(out),
            Self::Xml => XmlFormat.print_header_raw(
                out,
                query_mode as u32,
                &config.matrix,
                matrix.gap_open(),
                matrix.gap_extend(),
                config.max_evalue,
                first_title,
                first_len as u32,
                &crate::util::io::source_os_bytes(config.database.as_os_str())?,
            ),
            Self::Sam => crate::output::sam_format::SamFormat.print_header_bytes(
                out,
                query_mode as u32,
                &config.sam.invocation,
            ),
            Self::Paf
            | Self::Null
            | Self::Edge
            | Self::Bin1
            | Self::Clustering
            | Self::Taxon { .. } => Ok(()),
        }
    }

    pub(crate) fn footer<W: Write>(self, out: &mut W) -> io::Result<()> {
        match self {
            Self::Pairwise => PairwiseFormat.print_footer(out),
            Self::Xml => XmlFormat.print_footer(out),
            Self::Paf
            | Self::Sam
            | Self::Null
            | Self::Edge
            | Self::Bin1
            | Self::Clustering
            | Self::Taxon { .. } => Ok(()),
        }
    }

    /// Source Output::Info query intro, ordered target/HSP emission and epilog.
    /// Own each completed HSP once, retaining its transcript and source frames.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn query<S: SearchSubjects + ?Sized>(
        self,
        query_id: u32,
        query_oid: u64,
        title: &[u8],
        source_len: i32,
        frames: &[&[Letter]],
        matches: Vec<Match>,
        database: &S,
        matrix: &ScoreMatrix,
        translated: bool,
        db_letters: u64,
        report_unaligned: i32,
        sam: &SamOptions,
        controls: &super::blastp::SearchOutputOptions,
        query_self_aln_score: Option<f64>,
    ) -> io::Result<(Vec<u8>, u64)> {
        self.query_with_query_mode(
            query_id,
            query_oid,
            title,
            source_len,
            frames,
            matches,
            database,
            matrix,
            if translated { 3 } else { 2 },
            db_letters,
            report_unaligned,
            sam,
            controls,
            query_self_aln_score,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn query_with_query_mode<S: SearchSubjects + ?Sized>(
        self,
        query_id: u32,
        query_oid: u64,
        title: &[u8],
        source_len: i32,
        frames: &[&[Letter]],
        matches: Vec<Match>,
        database: &S,
        matrix: &ScoreMatrix,
        query_mode: i32,
        db_letters: u64,
        report_unaligned: i32,
        sam: &SamOptions,
        controls: &super::blastp::SearchOutputOptions,
        query_self_aln_score: Option<f64>,
    ) -> io::Result<(Vec<u8>, u64)> {
        let translated = query_mode == 3;
        let blastn = query_mode == 4;
        let mut out = Vec::new();
        let unaligned = matches.iter().all(|target| target.hsps.is_empty());
        if unaligned && !self.report_unaligned(report_unaligned) {
            return Ok((out, 0));
        }
        match self {
            Self::Pairwise => {
                PairwiseFormat.print_query_intro_bytes(&mut out, title, source_len, unaligned)?
            }
            Self::Xml => {
                XmlFormat.print_query_intro_bytes(&mut out, query_oid, title, source_len)?
            }
            Self::Paf => PafFormat.print_query_intro_bytes(&mut out, title, unaligned)?,
            Self::Sam => crate::output::sam_format::SamFormat
                .print_query_intro_bytes(&mut out, title, unaligned)?,
            Self::Bin1 => {
                crate::mcl::clustering_format::Bin1Format::print_query_intro(query_id, &mut out)
            }
            Self::Null | Self::Edge | Self::Clustering | Self::Taxon { .. } => (),
        }
        let mut context = HspContext {
            sequence_type: if blastn {
                crate::basic::value::SequenceType::Nucleotide
            } else {
                crate::basic::value::SequenceType::AminoAcid
            },
            query: if matches!(
                self,
                Self::Paf | Self::Null | Self::Edge | Self::Taxon { .. }
            ) {
                Vec::new()
            } else {
                frames.iter().map(|frame| frame.to_vec()).collect()
            },
            query_self_aln_score: query_self_aln_score.unwrap_or(0.0),
            query_source_len: source_len,
            query_len: source_len,
            query_title: crate::util::sequence::c_string_bytes(title).to_vec(),
            query_oid,
            query_id,
            ..Default::default()
        };
        let mut taxon = match self {
            Self::Taxon { include_lineage } => Some(crate::output::taxon_format::TaxonFormat::new(
                include_lineage,
            )),
            _ => None,
        };
        let mut accession_stats = AccessionParsing::default();
        let mut rows = 0;
        for (hit_num, target) in matches.into_iter().enumerate() {
            let target_id = target.target_block_id as usize;
            context.target_title =
                crate::util::sequence::c_string_bytes(database.title_bytes(target_id)?).to_vec();
            context.subject_oid = database.original_oid(target_id);
            context.subject_len = database.length(target_id) as i32;
            context.target_self_aln_score = database.self_aln_score(target_id).unwrap_or(0.0);
            context.hit_num = hit_num as u32;
            for (hsp_num, hsp) in target.hsps.into_iter().enumerate() {
                context.hsp = hsp;
                context.hsp_num = hsp_num as u32;
                match self {
                    Self::Pairwise => PairwiseFormat
                        .print_match_bytes(&mut out, &context, matrix, translated, blastn)?,
                    Self::Xml => XmlFormat.print_match_bytes(
                        &mut out,
                        &context,
                        matrix,
                        translated,
                        controls.xml_blord_format,
                        controls.no_parse_seqids,
                        &mut accession_stats,
                    )?,
                    Self::Paf => PafFormat
                        .print_match_with_command_bytes(&mut out, &context, matrix, blastn)?,
                    Self::Null => crate::output::output_format::NullFormat
                        .print_match_bytes(&mut out, &context)?,
                    Self::Edge => crate::output::edge::Edge.print_match(&mut out, &context)?,
                    Self::Bin1 => {
                        crate::mcl::clustering_format::Bin1Format::print_match(&context, &mut out)
                    }
                    Self::Clustering => controls
                        .clustering_filter
                        .as_ref()
                        .ok_or_else(|| io::Error::other("Missing clustering expression"))?
                        .format
                        .print_match(&context, translated, &mut out)
                        .map_err(io::Error::other)?,
                    Self::Taxon { .. } => taxon.as_mut().expect("taxonomy format").print_match(
                        &context,
                        database.taxids(target_id),
                        database
                            .taxonomy()
                            .ok_or_else(|| io::Error::other("Missing search taxonomy metadata"))?,
                    ),
                    Self::Sam => crate::output::sam_format::SamFormat.print_match_bytes(
                        &mut out,
                        &context,
                        matrix,
                        sam.all_subject_titles,
                        false,
                        sam.query_len,
                    )?,
                }
                rows += 1;
            }
        }
        match self {
            Self::Pairwise => PairwiseFormat.print_query_epilog(&mut out)?,
            Self::Xml => XmlFormat.print_query_epilog(
                &mut out,
                unaligned,
                database.sequence_count(),
                db_letters,
                matrix.k(),
                matrix.lambda(),
            )?,
            Self::Taxon { .. } => taxon
                .as_ref()
                .expect("taxonomy format")
                .print_query_epilog_bytes(
                    &mut out,
                    title,
                    database
                        .taxonomy()
                        .ok_or_else(|| io::Error::other("Missing search taxonomy metadata"))?,
                )?,
            Self::Paf | Self::Sam | Self::Null | Self::Edge | Self::Bin1 | Self::Clustering => (),
        }
        Ok((out, rows))
    }
}
