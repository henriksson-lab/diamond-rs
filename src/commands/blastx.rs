use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

use crate::basic::translate;
use crate::basic::value::{Letter, SequenceType};
use crate::commands::blastp::{BlastpConfig, TranslatedQueryLayout, TranslatedQuerySource};
use crate::config::Sensitivity;
use crate::data::fasta::{self, FastaRecord};
use crate::stats::cbs::CbsMode;
use crate::util::sequence::find_orfs;

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
    pub memory_limit: Option<usize>,
    pub tmpdir: PathBuf,
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

    // Load query DNA sequences
    let mut dna_records = Vec::new();
    for qf in &config.query_files {
        let records = fasta::read_fasta_file(Path::new(qf), SequenceType::Nucleotide)?;
        dna_records.extend(records);
    }
    eprintln!("Queries: {} DNA sequences", dna_records.len());

    // Translate to protein in all 6 frames. Each frame becomes a private
    // protein query; the grouped extension path maps results back to the
    // source DNA query and DNA coordinates before formatting output.
    let mut protein_records = Vec::new();
    let mut translated_sources = Vec::with_capacity(dna_records.len());

    for (dna_index, dna_rec) in dna_records.iter().enumerate() {
        let dna_bytes: Vec<u8> = dna_rec.sequence.iter().map(|&l| l as u8).collect();
        let dna_len = dna_bytes.len() as i32;
        let frames =
            translate::translate_6_frames_with_genetic_code(&dna_bytes, config.query_gencode)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let dna_short: String = dna_rec
            .id
            .split(|c: char| crate::util::sequence::ID_DELIMITERS.contains(c))
            .next()
            .unwrap_or("")
            .to_string();
        translated_sources.push(TranslatedQuerySource {
            id: dna_short.clone(),
            dna_len,
        });

        for (frame_idx, frame_seq) in frames.iter().enumerate() {
            let frame = translate::Frame::from_index(frame_idx as i32);
            // Include the input ordinal in the private frame id. Two input
            // records are allowed to have the same printed FASTA id, but they
            // must still be culled independently as two source queries.
            let frame_id = format!(
                "{}_diamond_query{}_frame{}",
                dna_short,
                dna_index,
                frame.signed_frame()
            );
            let active = (frame_idx < 3 && use_forward) || (frame_idx >= 3 && use_reverse);
            let mut sequence: Vec<Letter> = if active && !frame_seq.is_empty() {
                frame_seq.iter().map(|&b| b as Letter).collect()
            } else {
                // Keep a six-record layout without giving disabled/empty
                // contexts any seedable residues.
                vec![crate::basic::value::MASK_LETTER]
            };
            if active {
                let min_orf = min_orf_len(config.min_orf.unwrap_or(0), sequence.len());
                find_orfs(&mut sequence, min_orf as i32);
            }
            protein_records.push(FastaRecord {
                id: frame_id,
                sequence,
            });
        }
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
        translated_query_layout: Some(TranslatedQueryLayout {
            sources: translated_sources,
        }),
    };

    // The native blastp front end currently accepts FASTA paths, so use a
    // temporary translated-query file for the internal handoff.
    // Write translated sequences to a temporary handoff file and run blastp.
    let temp_tag = format!(
        "{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    );
    let temp_dir = if config.tmpdir.as_os_str().is_empty() {
        std::env::temp_dir()
    } else {
        config.tmpdir.clone()
    };
    let tmp_query = temp_dir.join(format!("diamond_blastx_query_{temp_tag}.faa"));
    let translated_frame_count = protein_records.len();
    let dna_query_count = dna_records.len();
    eprintln!(
        "Translated {dna_query_count} DNA queries into {translated_frame_count} protein frames in {:.1}s",
        start.elapsed().as_secs_f64()
    );
    {
        let mut f = BufWriter::new(std::fs::File::create(&tmp_query)?);
        for rec in &protein_records {
            writeln!(f, ">{}", rec.id)?;
            for chunk in rec.sequence.chunks(60) {
                for &l in chunk {
                    let idx = (l & 0x1F) as usize;
                    if idx < crate::basic::value::AMINO_ACID_ALPHABET.len() {
                        f.write_all(&[crate::basic::value::AMINO_ACID_ALPHABET[idx]])?;
                    }
                }
                writeln!(f)?;
            }
        }
    }
    // The inner blastp reader owns the translated sequences from here. Do not
    // retain a second in-memory copy throughout seeding and extension.
    drop(protein_records);
    drop(dna_records);

    let mut bp_config = blastp_config;
    bp_config.query_files = vec![tmp_query.to_string_lossy().to_string()];
    let result = crate::commands::blastp::run(&bp_config);

    let _ = std::fs::remove_file(&tmp_query);
    result
}

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
