use std::io;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

use crate::basic::value::SequenceType;
use crate::commands::blastp::{BlastpConfig, TranslatedQueryLayout};
use crate::config::Sensitivity;
use crate::data::block::Block;
use crate::data::fasta::{self, FastaRecord};
use crate::stats::cbs::CbsMode;

/// Configuration for a blastx run.
pub struct BlastxConfig {
    pub query_files: Vec<String>,
    pub database: String,
    pub output: Option<String>,
    pub matrix: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub max_evalue: f64,
    pub max_target_seqs: i64,
    pub ext_chunk_size: i64,
    pub toppercent: Option<f64>,
    pub global_ranking_targets: i64,
    pub min_id: f64,
    pub threads: i32,
    pub outfmt: Vec<String>,
    pub sensitivity: Sensitivity,
    pub query_gencode: u32,
    pub strand: String,
    pub min_orf: Option<u32>,
    /// Masking mode forwarded to the downstream blastp pipeline. C++ blastx
    /// accepts `--masking` and applies the chosen masker to translated
    /// protein frames; the previous hardcoded `MaskingMode::Tantan` silently
    /// ignored user `--masking 0` / `--masking seg`.
    pub masking: crate::masking::MaskingMode,
    /// `--motif-masking` value (empty = use sensitivity-default).
    pub motif_masking: String,
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
    let start = Instant::now();
    let (use_forward, use_reverse) = match config.strand.as_str() {
        "both" => (true, true),
        "plus" => (true, false),
        "minus" => (false, true),
        value => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid value for --strand: {value}"),
            ));
        }
    };

    // Upstream FastaFile round-robins two blastx inputs, producing
    // R1[0], R2[0], R1[1], R2[1], ... and rejects unsynchronised files.
    // This ordering is observable in output and is also what makes adjacent
    // source-query ordinals mates for `full_qseq_mate`.
    let (dna_records, paired) = load_query_records(&config.query_files)?;
    eprintln!("Queries: {} DNA sequences", dna_records.len());

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
        let dna_short: String = dna_rec
            .id
            .split(|c: char| crate::util::sequence::ID_DELIMITERS.contains(c))
            .next()
            .unwrap_or("")
            .to_string();
        query_block
            .push_back_translated_with_genetic_code(
                &dna_rec.sequence,
                Some(&dna_short),
                dna_index as u64,
                frame_mask,
                config.min_orf.unwrap_or(0) as i32,
                config.query_gencode,
            )
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
        mate_indices.push(paired.then_some(dna_index ^ 1));
    }
    // Reuse blastp's database loading and seed-search front end. The internal
    // translated layout makes its extension stage combine all six contexts.
    let blastp_config = BlastpConfig {
        query_files: vec![], // not used directly — we pass translated records
        database: config.database.clone(),
        output: config.output.clone(),
        matrix: config.matrix.clone(),
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
        max_evalue: config.max_evalue,
        max_target_seqs: config.max_target_seqs,
        ext_chunk_size: config.ext_chunk_size,
        toppercent: config.toppercent,
        global_ranking_targets: config.global_ranking_targets,
        min_id: config.min_id,
        threads: config.threads,
        outfmt: config.outfmt.clone(),
        sensitivity: config.sensitivity,
        // Forward user `--masking` / `--motif-masking` to the inner blastp
        // pipeline. Previously hardcoded to `Tantan`/empty, silently dropping
        // user overrides like `--masking 0` / `seg`.
        masking: config.masking,
        motif_masking: config.motif_masking.clone(),
        min_query_len: 0,
        query_cover: config.query_cover,
        subject_cover: config.subject_cover,
        comp_based_stats: config.comp_based_stats,
        no_self_hits: config.no_self_hits,
        ungapped_xdrop_bits: config.ungapped_xdrop_bits,
        memory_limit: config.memory_limit,
        tmpdir: config.tmpdir.clone(),
        target_indexed: config.target_indexed,
        translated_query_layout: Some(TranslatedQueryLayout { mate_indices }),
    };

    let translated_frame_count = query_block.seqs().len();
    eprintln!(
        "Translated {dna_query_count} DNA queries into {translated_frame_count} protein frames in {:.1}s",
        start.elapsed().as_secs_f64()
    );
    crate::commands::blastp::run_with_translated_query_block(&blastp_config, query_block)
}

fn load_query_records(paths: &[String]) -> io::Result<(Vec<FastaRecord>, bool)> {
    match paths {
        [] => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "blastx requires a query file",
        )),
        [path] => Ok((
            fasta::read_fasta_file(Path::new(path), SequenceType::Nucleotide)?,
            false,
        )),
        [left, right] => {
            let left = fasta::read_fasta_file(Path::new(left), SequenceType::Nucleotide)?;
            let right = fasta::read_fasta_file(Path::new(right), SequenceType::Nucleotide)?;
            Ok((interleave_paired_records(left, right)?, true))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "A maximum of 2 query files is supported in blastx mode.",
        )),
    }
}

fn interleave_paired_records(
    left: Vec<FastaRecord>,
    right: Vec<FastaRecord>,
) -> io::Result<Vec<FastaRecord>> {
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
