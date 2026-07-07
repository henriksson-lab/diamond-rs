use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::basic::translate::{self, Frame, Strand};
use crate::basic::value::{Letter, SequenceType};
use crate::commands::blastp::BlastpConfig;
use crate::config::Sensitivity;
use crate::data::fasta::{self, FastaRecord};
use crate::stats::cbs::CbsMode;
use crate::util::sequence::find_orfs;

/// Translation info for one protein frame query.
struct FrameInfo {
    original_id: String,
    frame: Frame,
    dna_len: i32,
}

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
}

/// Run blastx — translated DNA search against protein database.
///
/// Translates query DNA sequences in all 6 reading frames, then runs
/// blastp on each translated frame.
pub fn run(config: &BlastxConfig) -> io::Result<()> {
    let start = Instant::now();

    // Load query DNA sequences
    let mut dna_records = Vec::new();
    for qf in &config.query_files {
        let records = fasta::read_fasta_file(Path::new(qf), SequenceType::Nucleotide)?;
        dna_records.extend(records);
    }
    eprintln!("Queries: {} DNA sequences", dna_records.len());

    // Translate to protein in all 6 frames. Each frame becomes a synthetic
    // protein query named `{dna_id}_frame{N}`; we keep a map back to the
    // original id + frame + DNA length so we can rewrite the blastp output
    // into C++-equivalent (original id, DNA-space coords).
    let mut protein_records = Vec::new();
    let mut frame_info: HashMap<String, FrameInfo> = HashMap::new();
    let use_forward = config.strand != "minus";
    let use_reverse = config.strand != "plus";

    for dna_rec in &dna_records {
        let dna_bytes: Vec<u8> = dna_rec.sequence.iter().map(|&l| l as u8).collect();
        let dna_len = dna_bytes.len() as i32;
        let frames =
            translate::translate_6_frames_with_genetic_code(&dna_bytes, config.query_gencode)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        for (frame_idx, frame_seq) in frames.iter().enumerate() {
            // Skip frames based on strand setting
            if frame_idx < 3 && !use_forward {
                continue;
            }
            if frame_idx >= 3 && !use_reverse {
                continue;
            }

            // Skip empty frames
            if frame_seq.is_empty() {
                continue;
            }

            let frame = translate::Frame::from_index(frame_idx as i32);
            // Key `frame_info` by the *post-`print_title`* short id, since
            // the blastp inner pipeline emits column 0 via `print_title`
            // (`src/output/format.rs:1469`), which truncates at the first
            // ID delimiter (whitespace, \x01, etc.). If we keyed by the
            // raw `dna_rec.id`, queries whose header has whitespace would
            // miss the lookup and emit unchanged protein-space coords.
            let dna_short: String = dna_rec
                .id
                .split(|c: char| crate::util::sequence::ID_DELIMITERS.contains(c))
                .next()
                .unwrap_or("")
                .to_string();
            let frame_id = format!("{}_frame{}", dna_short, frame.signed_frame());
            let mut sequence: Vec<Letter> = frame_seq.iter().map(|&b| b as Letter).collect();
            let min_orf = min_orf_len(config.min_orf.unwrap_or(0), sequence.len());
            find_orfs(&mut sequence, min_orf as i32);
            frame_info.insert(
                frame_id.clone(),
                FrameInfo {
                    original_id: dna_short,
                    frame,
                    dna_len,
                },
            );
            protein_records.push(FastaRecord {
                id: frame_id,
                sequence,
            });
        }
    }
    eprintln!(
        "Translated frames: {} protein sequences",
        protein_records.len()
    );

    // Build a blastp config and delegate
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
    };

    // Run the blastp pipeline with translated sequences
    // For now, delegate to the FFI for full blastx compatibility
    // (the native blastp pipeline expects FASTA file inputs)
    eprintln!(
        "Translated {} DNA queries into {} protein frames in {:.1}s",
        dna_records.len(),
        protein_records.len(),
        start.elapsed().as_secs_f64()
    );

    // Write translated sequences to temp file and run blastp on them
    let temp_tag = format!(
        "{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    );
    let tmp_query = std::env::temp_dir().join(format!("diamond_blastx_query_{temp_tag}.faa"));
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

    // Run blastp on translated queries, redirecting its output to a temp file
    // so we can post-process: strip `_frame{N}` suffixes from query IDs and
    // convert qstart/qend from protein-space back to DNA-space.
    let final_output = blastp_config.output.clone();
    let tmp_bp_out = std::env::temp_dir().join(format!("diamond_blastx_bp_{temp_tag}.tsv"));
    let mut bp_config = blastp_config;
    bp_config.query_files = vec![tmp_query.to_string_lossy().to_string()];
    bp_config.output = Some(tmp_bp_out.to_string_lossy().to_string());
    let result = crate::commands::blastp::run(&bp_config);

    let _ = std::fs::remove_file(&tmp_query);
    result?;

    // Post-process: translate protein coords back to DNA coords.
    let bp_in = BufReader::new(std::fs::File::open(&tmp_bp_out)?);
    let mut writer: Box<dyn Write> = match final_output {
        Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    for line in bp_in.lines() {
        let line = line?;
        if line.is_empty() {
            writeln!(writer)?;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 12 {
            // Unknown format — emit as-is.
            writeln!(writer, "{}", line)?;
            continue;
        }
        let info = match frame_info.get(cols[0]) {
            Some(f) => f,
            None => {
                // No frame metadata (shouldn't happen) — emit unchanged.
                writeln!(writer, "{}", line)?;
                continue;
            }
        };
        // Default tabular: 0:qseqid 1:sseqid 2:pident 3:length 4:mismatch
        //                  5:gapopen 6:qstart 7:qend 8:sstart 9:send
        //                  10:evalue 11:bitscore
        let q_start_prot: i32 = cols[6].parse().unwrap_or(0);
        let q_end_prot: i32 = cols[7].parse().unwrap_or(0);
        let (q_start_dna, q_end_dna) =
            protein_to_dna_coords(q_start_prot, q_end_prot, info.frame, info.dna_len);
        let mut out_cols: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
        out_cols[0] = info.original_id.clone();
        out_cols[6] = q_start_dna.to_string();
        out_cols[7] = q_end_dna.to_string();
        writeln!(writer, "{}", out_cols.join("\t"))?;
    }
    writer.flush()?;
    let _ = std::fs::remove_file(&tmp_bp_out);
    Ok(())
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

/// Convert protein-space [qstart, qend] (1-based, inclusive) to the equivalent
/// DNA-space range. For forward frames the DNA range goes ascending; for reverse
/// frames it descends, matching C++ DIAMOND blastx output which reports the
/// reverse strand with qstart > qend.
fn protein_to_dna_coords(
    q_start_prot: i32,
    q_end_prot: i32,
    frame: Frame,
    dna_len: i32,
) -> (i32, i32) {
    // Each protein residue corresponds to 3 DNA bases. For forward frame f
    // (offset 0..=2), protein position i (1-based) maps to DNA bases at
    // 1-based positions [f + 3*(i-1) + 1, f + 3*i].
    match frame.strand {
        Strand::Forward => {
            let f = frame.offset;
            let dna_start = f + 3 * (q_start_prot - 1) + 1;
            let dna_end = f + 3 * q_end_prot;
            (dna_start, dna_end)
        }
        Strand::Reverse => {
            // Reverse-complement strand: protein pos i in the reversed sequence
            // corresponds to original DNA positions [L - f - 3*i + 1, L - f - 3*(i-1)].
            let f = frame.offset;
            let dna_start = dna_len - f - 3 * (q_start_prot - 1);
            let dna_end = dna_len - f - 3 * q_end_prot + 1;
            (dna_start, dna_end)
        }
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
