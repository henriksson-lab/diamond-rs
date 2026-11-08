use std::io;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

use crate::basic::value::SequenceType;
use crate::commands::blastp::{BlastpConfig, TranslatedQueryLayout};
use crate::config::Sensitivity;
use crate::data::block::Block;
use crate::data::fasta::{self, ByteFastaRecord, FastaRecord};
use crate::stats::cbs::CbsMode;

/// Configuration for a blastx run.
pub struct BlastxConfig {
    pub query_files: Vec<std::path::PathBuf>,
    pub database: std::path::PathBuf,
    pub output: Option<std::path::PathBuf>,
    pub matrix: Vec<u8>,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub max_evalue: f64,
    pub min_bit_score: f64,
    pub max_target_seqs: i64,
    /// Maximum HSPs per query/target pair; zero means unlimited.
    pub max_hsps: u32,
    /// Source --unal setting: -1 selects the format default.
    pub report_unaligned: i32,
    /// Include scientific-name lineage in taxonomy format 102.
    pub include_lineage: bool,
    /// Absence differs from a bare --header (the search workflow default).
    pub output_header: Option<Vec<Vec<u8>>>,
    /// Source compression option: empty/0=plain, 1=gzip, zstd=optional build.
    pub compression: Vec<u8>,
    pub file_buffer_size: usize,
    pub sam: crate::commands::blastp::SamOptions,
    pub archive: crate::commands::blastp::ArchiveOptions,
    pub tuning: crate::commands::blastp::SearchTuning,
    pub scoring_matrix: Option<std::sync::Arc<crate::stats::score_matrix::ScoreMatrix>>,
    pub auxiliary_output: super::query_output::AuxiliaryOutputOptions,
    pub taxonlist: Vec<u8>,
    pub taxon_exclude: Vec<u8>,
    pub seqidlist: std::path::PathBuf,
    pub skip_missing_seqids: bool,
    pub range_culling: bool,
    pub range_cover: f64,
    pub culling_overlap: f64,
    pub taxon_k: u64,
    pub query_or_subject_cover: f64,
    pub reference_taxon_rank_ids: Option<std::sync::Arc<[Vec<u32>]>>,
    pub ext_chunk_size: i64,
    pub toppercent: Option<f64>,
    pub global_ranking_targets: i64,
    /// Frameshift penalty; zero selects ordinary independent-frame extension.
    pub frame_shift: i32,
    /// Explicit source extension mode; None selects its workflow/sensitivity default.
    pub extension_mode: Option<crate::search::sensitivity::ExtensionMode>,
    /// Source sensitivity rounds; Some(empty) selects the default schedule.
    pub iterate: Option<Vec<Vec<u8>>>,
    pub min_id: f64,
    /// Minimum score-derived approximate identity percentage.
    pub approx_min_id: f64,
    /// Explicit seed-only masking, added to sensitivity-default motif masking.
    pub soft_masking: crate::masking::MaskingAlgo,
    pub threads: i32,
    /// Seed-index partitions processed per pass; zero selects the source default.
    pub index_chunks: u32,
    /// Query/reference letter budget in billions; None uses the source default.
    pub block_size: Option<f64>,
    /// Force hashed query-indexed seed filtering (--algo 1).
    pub query_indexed: bool,
    /// Source --algo ctg: dense query seed filter and contiguous encoding.
    pub contiguous_seed_mode: bool,
    pub min_query_len: i32,
    /// Exhaustively align every query against every reference without seeds.
    pub swipe_all: bool,
    pub outfmt: Vec<Vec<u8>>,
    pub sensitivity: Sensitivity,
    pub query_gencode: u32,
    pub strand: Vec<u8>,
    pub min_orf: Option<u32>,
    /// Masking mode forwarded to the downstream blastp pipeline. C++ blastx
    /// accepts `--masking` and applies the chosen masker to translated
    /// protein frames; the previous hardcoded `MaskingMode::Tantan` silently
    /// ignored user `--masking 0` / `--masking seg`.
    pub masking: crate::masking::MaskingMode,
    /// `--motif-masking` value (empty = use sensitivity-default).
    pub motif_masking: Vec<u8>,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub comp_based_stats: CbsMode,
    pub no_self_hits: bool,
    pub ungapped_xdrop_bits: f64,
    pub memory_limit: Option<u64>,
    pub tmpdir: PathBuf,
    pub target_indexed: bool,
}

/// Run blastx — translated DNA search against protein database.
///
/// Translates each DNA query into six protein contexts, uses the native seed
/// search, then ranks and extends all contexts as one source query.
pub fn run(config: &BlastxConfig) -> io::Result<()> {
    CbsMode::from_search_code(
        config.comp_based_stats as u32,
        true,
        config.global_ranking_targets,
    )
    .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
    let block_config = blastp_config(config, Vec::new());
    if let Some(budget) = super::blastp::native_block_budget(&block_config)? {
        return super::search_blocks::run_translated_files(config, block_config, budget);
    }
    let start = Instant::now();
    let (use_forward, use_reverse) = match config.strand.as_slice() {
        b"both" => (true, true),
        b"plus" => (true, false),
        b"minus" => (false, true),
        value => {
            let mut message = b"invalid value for --strand: ".to_vec();
            message.extend_from_slice(value);
            return Err(crate::util::io::RawDiagnosticError(message).into_io());
        }
    };

    // Upstream FastaFile round-robins two blastx inputs, producing
    // R1[0], R2[0], R1[1], R2[1], ... and rejects unsynchronised files.
    // This ordering is observable in output and is also what makes adjacent
    // source-query ordinals mates for `full_qseq_mate`.
    let (dna_records, paired) = load_query_records(&config.query_files, config.file_buffer_size)?;
    search_message!("Queries: {} DNA sequences", dna_records.len());

    // Upstream loads source DNA and all six protein contexts into one Block.
    // Build that same contiguous owner directly: the downstream search masks,
    // restores, fingerprints, and extends these exact sequences rather than
    // retaining six separately allocated FastaRecords plus a second arena.
    let mut query_block = Block::new();
    let dna_query_count = dna_records.len();
    let mut mate_indices = Vec::with_capacity(dna_records.len());
    let frame_mask =
        (if use_forward { 0b000111 } else { 0 }) | (if use_reverse { 0b111000 } else { 0 });

    for (dna_index, dna_rec) in dna_records.into_iter().enumerate() {
        let full_query_title = config
            .outfmt
            .first()
            .is_some_and(|name| matches!(name.as_slice(), b"0" | b"5" | b"xml"));
        let dna_title = if full_query_title {
            dna_rec.id.clone()
        } else {
            crate::util::sequence::seqid_bytes(&dna_rec.id)
        };
        if let Some(translator) = config.tuning.runtime.translator() {
            query_block.push_back_translated_with_translator_bytes(
                &dna_rec.sequence,
                Some(dna_title.as_slice()),
                dna_index as u64,
                frame_mask,
                config.min_orf.unwrap_or(0) as i32,
                translator,
            )
        } else {
            query_block.push_back_translated_with_genetic_code_bytes(
                &dna_rec.sequence,
                Some(dna_title.as_slice()),
                dna_index as u64,
                frame_mask,
                config.min_orf.unwrap_or(0) as i32,
                config.query_gencode,
            )
        }
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
        if config
            .outfmt
            .iter()
            .any(|field| matches!(field.as_slice(), b"qqual" | b"full_qqual"))
            || config.auxiliary_output.needs_quality()
            || config.tuning.store_query_quality
        {
            query_block.push_source_quality(dna_rec.quality.as_deref().unwrap_or(&[]));
        }
        mate_indices.push(paired.then_some(dna_index ^ 1));
    }
    // Reuse blastp's database loading and seed-search front end. The internal
    // translated layout makes its extension stage combine all six contexts.
    let blastp_config = blastp_config(config, mate_indices);

    let translated_frame_count = query_block.seqs().len();
    search_message!(
        "Translated {dna_query_count} DNA queries into {translated_frame_count} protein frames in {:.1}s",
        start.elapsed().as_secs_f64()
    );
    crate::commands::blastp::run_with_translated_query_block(&blastp_config, query_block)
}

pub(crate) fn blastp_config(
    config: &BlastxConfig,
    mate_indices: Vec<Option<usize>>,
) -> BlastpConfig {
    BlastpConfig {
        query_files: vec![], // not used directly — we pass translated records
        database: config.database.clone(),
        output: config.output.clone(),
        matrix: config.matrix.clone(),
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
        max_evalue: config.max_evalue,
        min_bit_score: config.min_bit_score,
        max_target_seqs: config.max_target_seqs,
        max_hsps: config.max_hsps,
        report_unaligned: config.report_unaligned,
        include_lineage: config.include_lineage,
        output_header: config.output_header.clone(),
        compression: config.compression.clone(),
        file_buffer_size: config.file_buffer_size,
        sam: config.sam.clone(),
        archive: config.archive.clone(),
        tuning: config.tuning.clone(),
        scoring_matrix: config.scoring_matrix.clone(),
        auxiliary_output: config.auxiliary_output.clone(),
        taxonlist: config.taxonlist.clone(),
        taxon_exclude: config.taxon_exclude.clone(),
        seqidlist: config.seqidlist.clone(),
        skip_missing_seqids: config.skip_missing_seqids,
        range_culling: config.range_culling,
        range_cover: config.range_cover,
        culling_overlap: config.culling_overlap,
        taxon_k: config.taxon_k,
        query_or_subject_cover: config.query_or_subject_cover,
        reference_taxon_rank_ids: config.reference_taxon_rank_ids.clone(),
        ext_chunk_size: config.ext_chunk_size,
        toppercent: config.toppercent,
        global_ranking_targets: config.global_ranking_targets,
        frame_shift: config.frame_shift,
        extension_mode: config.extension_mode,
        iterate: config.iterate.clone(),
        round_linear_target: false,
        round_query_skip: None,
        min_id: config.min_id,
        approx_min_id: config.approx_min_id,
        soft_masking: config.soft_masking,
        threads: config.threads,
        index_chunks: config.index_chunks,
        block_size: config.block_size,
        query_indexed: config.query_indexed,
        contiguous_seed_mode: config.contiguous_seed_mode,
        swipe_all: config.swipe_all,
        outfmt: config.outfmt.clone(),
        sensitivity: config.sensitivity,
        // Forward user `--masking` / `--motif-masking` to the inner blastp
        // pipeline. Previously hardcoded to `Tantan`/empty, silently dropping
        // user overrides like `--masking 0` / `seg`.
        masking: config.masking,
        motif_masking: config.motif_masking.clone(),
        min_query_len: config.min_query_len,
        query_cover: config.query_cover,
        subject_cover: config.subject_cover,
        comp_based_stats: config.comp_based_stats,
        no_self_hits: config.no_self_hits,
        ungapped_xdrop_bits: config.ungapped_xdrop_bits,
        memory_limit: config.memory_limit,
        tmpdir: config.tmpdir.clone(),
        target_indexed: config.target_indexed,
        translated_query_layout: Some(TranslatedQueryLayout { mate_indices }),
    }
}

fn load_query_records(
    paths: &[PathBuf],
    file_buffer_size: usize,
) -> io::Result<(Vec<ByteFastaRecord>, bool)> {
    let read = |path: &Path| {
        let mut records = Vec::new();
        fasta::for_each_fasta_file_with_buffer_size_bytes(
            path,
            SequenceType::Nucleotide,
            file_buffer_size,
            |record| {
                records.push(record);
                Ok(())
            },
        )?;
        Ok::<_, io::Error>(records)
    };
    match paths {
        [] => {
            crate::util::io::console::source_eprintln!(
                "Query file parameter (--query/-q) is missing. Input will be read from stdin."
            );
            Ok((read(Path::new(""))?, false))
        }
        [path] => Ok((read(path)?, false)),
        [left, right] => Ok((interleave_paired_records(read(left)?, read(right)?)?, true)),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "A maximum of 2 query files is supported in blastx mode.",
        )),
    }
}

fn interleave_paired_records<T>(
    left: Vec<FastaRecord<T>>,
    right: Vec<FastaRecord<T>>,
) -> io::Result<Vec<FastaRecord<T>>> {
    if left.len() != right.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Unequal number of sequences in paired read files.",
        ));
    }
    let mut records = Vec::with_capacity(left.len().saturating_mul(2));
    for (left, right) in left.into_iter().zip(right) {
        records.push(left);
        records.push(right);
    }
    Ok(records)
}

#[cfg(test)]
fn min_orf_len(run_len: u32, length: usize) -> u32 {
    if run_len != 0 {
        return run_len;
    }
    if length < 30 {
        1
    } else if length < 100 {
        20
    } else {
        40
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::translate;
    use crate::basic::value::Letter;
    use crate::util::sequence::find_orfs;

    fn record(id: &str, base: Letter) -> FastaRecord {
        FastaRecord {
            quality: None,
            id: id.to_owned(),
            sequence: vec![base],
        }
    }

    #[test]
    fn paired_query_records_are_round_robin_interleaved_like_upstream() {
        let records = interleave_paired_records(
            vec![record("r1/1", 0), record("r2/1", 1)],
            vec![record("r1/2", 2), record("r2/2", 3)],
        )
        .unwrap();
        assert_eq!(
            records
                .iter()
                .map(|record| record.id.as_str())
                .collect::<Vec<_>>(),
            ["r1/1", "r1/2", "r2/1", "r2/2"]
        );
    }

    #[test]
    fn paired_query_records_reject_unequal_counts_like_upstream() {
        let error = interleave_paired_records(vec![record("r1/1", 0)], Vec::new()).unwrap_err();
        assert_eq!(
            error.to_string(),
            "Unequal number of sequences in paired read files."
        );
    }

    #[test]
    fn test_blastx_translation() {
        // Test that DNA translation produces valid protein frames
        let dna = b">test\nATGCGATCGATCG\n";
        let records = fasta::read_fasta_nucleotide(&dna[..]).unwrap();
        assert_eq!(records.len(), 1);

        let dna_bytes: Vec<u8> = records[0].sequence.iter().map(|&l| l as u8).collect();
        let frames = translate::translate_6_frames(&dna_bytes);

        // Should have 6 frames
        assert_eq!(frames.len(), 6);
        // Forward frame 0: ATG CGA TCG ATC -> M R S I (4 codons from 13 bases)
        assert_eq!(frames[0].len(), 4);
    }

    #[test]
    fn test_min_orf_len_matches_cpp_defaults() {
        assert_eq!(min_orf_len(0, 29), 1);
        assert_eq!(min_orf_len(0, 30), 20);
        assert_eq!(min_orf_len(0, 99), 20);
        assert_eq!(min_orf_len(0, 100), 40);
        assert_eq!(min_orf_len(7, 100), 7);
    }

    #[test]
    fn test_min_orf_masks_short_stop_to_stop_runs() {
        let mut seq = vec![0, crate::basic::value::STOP_LETTER, 1, 2, 3];
        let min_orf = min_orf_len(3, seq.len()) as i32;
        find_orfs(&mut seq, min_orf);
        assert_eq!(
            seq,
            vec![
                crate::basic::value::MASK_LETTER,
                crate::basic::value::STOP_LETTER,
                1,
                2,
                3
            ]
        );
    }
}
