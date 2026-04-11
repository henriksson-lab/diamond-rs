use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::basic::translate;
use crate::basic::value::{Letter, SequenceType};
use crate::commands::blastp::BlastpConfig;
use crate::config::Sensitivity;
use crate::data::fasta::{self, FastaRecord};

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
    pub min_id: f64,
    pub threads: i32,
    pub outfmt: Vec<String>,
    pub sensitivity: Sensitivity,
    pub query_gencode: u32,
    pub strand: String,
    pub min_orf: Option<u32>,
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

    // Translate to protein in all 6 frames
    let mut protein_records = Vec::new();
    let use_forward = config.strand != "minus";
    let use_reverse = config.strand != "plus";

    for dna_rec in &dna_records {
        let dna_bytes: Vec<u8> = dna_rec.sequence.iter().map(|&l| l as u8).collect();
        let frames = translate::translate_6_frames(&dna_bytes);

        for (frame_idx, frame_seq) in frames.iter().enumerate() {
            // Skip frames based on strand setting
            if frame_idx < 3 && !use_forward {
                continue;
            }
            if frame_idx >= 3 && !use_reverse {
                continue;
            }

            // Filter by minimum ORF length
            if let Some(min_orf) = config.min_orf {
                if (frame_seq.len() as u32) < min_orf {
                    continue;
                }
            }

            // Skip empty frames
            if frame_seq.is_empty() {
                continue;
            }

            let frame = translate::Frame::from_index(frame_idx as i32);
            let frame_id = format!(
                "{}_frame{}",
                dna_rec.id,
                frame.signed_frame()
            );
            let sequence: Vec<Letter> = frame_seq.iter().map(|&b| b as Letter).collect();
            protein_records.push(FastaRecord {
                id: frame_id,
                sequence,
            });
        }
    }
    eprintln!("Translated frames: {} protein sequences", protein_records.len());

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
        min_id: config.min_id,
        threads: config.threads,
        outfmt: config.outfmt.clone(),
        sensitivity: config.sensitivity,
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
    let tmp_query = std::env::temp_dir().join("diamond_blastx_query.faa");
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

    // Run blastp on translated queries
    let mut bp_config = blastp_config;
    bp_config.query_files = vec![tmp_query.to_string_lossy().to_string()];
    let result = crate::commands::blastp::run(&bp_config);

    let _ = std::fs::remove_file(&tmp_query);
    result
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
        // Forward frame 0: ATG CGG TCG ATC -> M R S I (4 codons from 13 bases)
        assert_eq!(frames[0].len(), 4);
    }
}
