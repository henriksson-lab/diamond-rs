use std::collections::HashMap;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::basic::reduction::Reduction;
use crate::basic::shape::Shape;
use crate::basic::value::{Letter, SequenceType, DELIMITER_LETTER};
use crate::config::Sensitivity;
use crate::data::fasta;
use crate::dp::ungapped;
use crate::dp::smith_waterman;
use crate::output::format::{self, FieldId, Hsp as OutputHsp};
use crate::search::{parallel, seed_match, sensitivity};
use crate::stats::score_matrix::ScoreMatrix;

/// Configuration for a blastp run.
pub struct BlastpConfig {
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
}

/// Run a simplified blastp search.
///
/// This implements the basic pipeline:
/// 1. Load query and database sequences
/// 2. Extract seeds from both using shapes
/// 3. Find seed matches (hash join)
/// 4. Extend hits with ungapped x-drop
/// 5. Perform gapped Smith-Waterman alignment
/// 6. Filter by e-value and output
pub fn run(config: &BlastpConfig) -> io::Result<()> {
    let start = Instant::now();

    // Load database sequences - try DMND first, then FASTA
    let mut db_records = {
        let db_path = Path::new(&config.database);
        let dmnd_path = if db_path.extension().is_some_and(|e| e == "dmnd") {
            db_path.to_path_buf()
        } else {
            db_path.with_extension("dmnd")
        };

        if dmnd_path.exists() {
            let (header, records) = crate::data::dmnd_reader::read_dmnd(&dmnd_path)?;
            eprintln!(
                "Database: {} sequences, {} letters (DMND)",
                header.sequences, header.letters
            );
            records
        } else {
            let fasta_path = if db_path.extension().is_none() {
                db_path.with_extension("faa")
            } else {
                db_path.to_path_buf()
            };
            let records = fasta::read_fasta_file(&fasta_path, SequenceType::AminoAcid)?;
            eprintln!("Database: {} sequences (FASTA)", records.len());
            records
        }
    };

    // Compute total database letters for E-value normalization
    let db_letters: u64 = db_records.iter().map(|r| r.sequence.len() as u64).sum();

    // Load scoring matrix with database size
    let score_matrix = ScoreMatrix::new(
        &config.matrix,
        config.gap_open,
        config.gap_extend,
        0,
        1,
        db_letters,
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // Load query sequences
    let mut query_records = Vec::new();
    for qf in &config.query_files {
        let records = fasta::read_fasta_file(Path::new(qf), SequenceType::AminoAcid)?;
        query_records.extend(records);
    }
    eprintln!("Queries: {} sequences", query_records.len());

    // Apply tantan soft masking for seed filtering (high bit set on masked positions).
    for record in &mut db_records {
        crate::masking::mask_sequence(&mut record.sequence, crate::masking::MaskingAlgo::Tantan);
    }
    for record in &mut query_records {
        crate::masking::mask_sequence(&mut record.sequence, crate::masking::MaskingAlgo::Tantan);
    }

    // Set up seed extraction using sensitivity-appropriate shapes
    let reduction = Reduction::default_reduction();
    let shape_codes = sensitivity::get_shape_codes(config.sensitivity);
    let shapes: Vec<Shape> = shape_codes
        .iter()
        .map(|code| Shape::from_code(code, &reduction))
        .collect();
    eprintln!(
        "Sensitivity: {:?}, shapes: {} (weights: {})",
        config.sensitivity,
        shapes.len(),
        shapes.iter().map(|s| s.weight.to_string()).collect::<Vec<_>>().join(",")
    );

    // Build partitioned seed arrays and join for each shape (parallel)
    let db_seqs: Vec<&[Letter]> = db_records.iter().map(|r| r.sequence.as_slice()).collect();
    let query_seqs: Vec<&[Letter]> = query_records.iter().map(|r| r.sequence.as_slice()).collect();
    let mut all_seed_matches = Vec::new();
    for shape in &shapes {
        let matches = parallel::find_seed_matches_partitioned(&query_seqs, &db_seqs, shape, &reduction);
        all_seed_matches.extend(matches);
    }
    eprintln!("Total seed matches: {}", all_seed_matches.len());

    // Parse output format
    let fields = if config.outfmt.is_empty() || config.outfmt[0] == "6" {
        if config.outfmt.len() > 1 {
            config.outfmt[1..]
                .iter()
                .filter_map(|f| FieldId::from_name(f))
                .collect()
        } else {
            format::DEFAULT_TABULAR_FIELDS.to_vec()
        }
    } else {
        format::DEFAULT_TABULAR_FIELDS.to_vec()
    };

    // Set up output writer
    let output: Box<dyn Write> = match &config.output {
        Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut writer = output;

    let mut total_alignments = 0u64;

    // Group all seed matches by query
    let mut query_matches: HashMap<u32, Vec<&seed_match::SeedMatch>> = HashMap::new();
    for m in &all_seed_matches {
        query_matches.entry(m.query_id).or_default().push(m);
    }

    // Process each query
    for (query_idx, query_rec) in query_records.iter().enumerate() {
        let query = &query_rec.sequence;

        // Pad query with delimiters for ungapped extension
        let mut padded_query = vec![DELIMITER_LETTER];
        padded_query.extend_from_slice(query);
        padded_query.push(DELIMITER_LETTER);

        // Compute ungapped score cutoff for this query length.
        // Matches C++ CutoffTable with ungapped_evalue=10000 (default sensitivity).
        // Uses the actual database size for normalization (not 1e9) to handle
        // small databases correctly.
        let ungapped_evalue = 10000.0;
        let ungapped_cutoff = score_matrix.ungapped_cutoff_db(query.len(), ungapped_evalue);

        // CBS (composition-based statistics) correction per query position.
        // C++ default is comp-based-stats=1 (Hauser correction, window=40).
        let query_cbs = crate::stats::cbs::hauser_correction(query, &score_matrix);

        // Get seed matches for this query
        let empty_matches = Vec::new();
        let query_seed_matches = query_matches.get(&(query_idx as u32)).unwrap_or(&empty_matches);

        // Group by target
        let mut target_hits: HashMap<u32, Vec<&&seed_match::SeedMatch>> = HashMap::new();
        for m in query_seed_matches {
            target_hits.entry(m.ref_id).or_default().push(m);
        }

        // Process each target
        let mut target_results: Vec<(u32, f64, OutputHsp)> = Vec::new();

        for (&target_id, hits) in &target_hits {
            let target = &db_records[target_id as usize].sequence;

            // Ungapped extension: try all seed hits, keep best
            let mut padded_target = vec![DELIMITER_LETTER];
            padded_target.extend_from_slice(target);
            padded_target.push(DELIMITER_LETTER);

            let mut best_ungapped_score = 0i32;
            for hit in hits {
                let qa = hit.query_pos as usize + 1;
                let sa = hit.ref_pos as usize + 1;
                if qa >= padded_query.len() || sa >= padded_target.len() {
                    continue;
                }
                let ungapped = ungapped::xdrop_ungapped(
                    &padded_query,
                    &padded_target,
                    qa,
                    sa,
                    12,
                    &score_matrix,
                );
                best_ungapped_score = best_ungapped_score.max(ungapped.score);
            }

            if best_ungapped_score < ungapped_cutoff {
                continue;
            }

            // Banded Smith-Waterman with CBS and masking-aware scoring.
            // Matches C++ banded_swipe.h: band constrains alignment near seed diagonal,
            // soft-masked query positions score 0 (matching SIMD zeroing), CBS added.
            // Use the best seed hit position as the anchor for banding.
            let best_hit = hits.iter().max_by_key(|h| {
                let qa = h.query_pos as usize;
                let sa = h.ref_pos as usize;
                if qa < query.len() && sa < target.len() {
                    score_matrix.score(query[qa] & crate::basic::value::LETTER_MASK,
                                       target[sa] & crate::basic::value::LETTER_MASK)
                } else { 0 }
            }).unwrap();
            let sw = crate::dp::banded_cbs::banded_sw_cbs(
                query, target,
                best_hit.query_pos as usize,
                best_hit.ref_pos as usize,
                30, // band_width matching C++ default
                &score_matrix,
                &query_cbs,
            );
            if sw.score <= 0 {
                continue;
            }

            let evalue = score_matrix.evalue(sw.score, query.len() as u32, target.len() as u32);
            if evalue > config.max_evalue {
                continue;
            }

            let pident = if sw.length > 0 {
                100.0 * sw.identities as f64 / sw.length as f64
            } else {
                0.0
            };
            if pident < config.min_id {
                continue;
            }

            let bit_score = score_matrix.bitscore(sw.score as f64);

            let hsp = OutputHsp {
                score: sw.score,
                evalue,
                bit_score,
                query_range: (sw.query_begin, sw.query_end),
                subject_range: (sw.subject_begin, sw.subject_end),
                query_source_range: (sw.query_begin, sw.query_end),
                subject_source_range: (sw.subject_begin, sw.subject_end),
                frame: 0,
                length: sw.length,
                identities: sw.identities,
                mismatches: sw.mismatches,
                positives: sw.identities, // simplified
                gap_openings: sw.gap_openings,
                gaps: sw.gaps,
            };

            target_results.push((target_id, evalue, hsp));
        }

        // Sort by e-value
        target_results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        // Apply max target seqs
        let max_targets = config.max_target_seqs as usize;
        if target_results.len() > max_targets {
            target_results.truncate(max_targets);
        }

        // Output results
        for (target_id, _, hsp) in &target_results {
            format::write_tabular_row(
                &mut writer,
                &query_rec.id,
                &db_records[*target_id as usize].id,
                hsp,
                &fields,
                query.len() as i32,
                db_records[*target_id as usize].sequence.len() as i32,
            )?;
            total_alignments += 1;
        }
    }

    writer.flush()?;

    let elapsed = start.elapsed();
    eprintln!(
        "Reported {} alignments in {:.1}s",
        total_alignments,
        elapsed.as_secs_f64()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blastp_with_dmnd() {
        let query = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/5.faa");
        let db = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/data.dmnd");
        let output_path = std::env::temp_dir().join("test_blastp_dmnd.out");

        let config = BlastpConfig {
            query_files: vec![query.to_string()],
            database: db.to_string(),
            output: Some(output_path.to_string_lossy().to_string()),
            matrix: "blosum62".to_string(),
            gap_open: 11,
            gap_extend: 1,
            max_evalue: 0.001,
            max_target_seqs: 25,
            min_id: 0.0,
            threads: 1,
            outfmt: vec![],
            sensitivity: Sensitivity::Default,
        };

        let result = run(&config);
        assert!(result.is_ok(), "blastp with DMND failed: {:?}", result.err());

        let output = std::fs::read_to_string(&output_path).unwrap();
        assert!(!output.is_empty(), "blastp produced no output");
        let _ = std::fs::remove_file(&output_path);
    }

    #[test]
    fn test_blastp_self() {
        let fasta_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/1.faa"
        );
        let output_path = std::env::temp_dir().join("test_blastp_self.out");

        let config = BlastpConfig {
            query_files: vec![fasta_path.to_string()],
            database: fasta_path.to_string(),
            output: Some(output_path.to_string_lossy().to_string()),
            matrix: "blosum62".to_string(),
            gap_open: 11,
            gap_extend: 1,
            max_evalue: 0.001,
            max_target_seqs: 25,
            min_id: 0.0,
            threads: 1,
            outfmt: vec![],
            sensitivity: Sensitivity::Default,
        };

        let result = run(&config);
        assert!(result.is_ok(), "blastp failed: {:?}", result.err());

        // Check output file exists and has content
        let output = std::fs::read_to_string(&output_path).unwrap();
        assert!(!output.is_empty(), "blastp produced no output");
        // Self-alignment should find at least one hit
        assert!(output.lines().count() >= 1);

        let _ = std::fs::remove_file(&output_path);
    }
}
