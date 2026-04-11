use std::collections::HashMap;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::basic::reduction::Reduction;
use crate::basic::shape::Shape;
use crate::basic::value::{Letter, SequenceType};
use crate::cluster::{Cluster, ClusterResult};
use crate::config::Sensitivity;
use crate::data::fasta;
use crate::dp::smith_waterman;
use crate::search::{seed_match, sensitivity};
use crate::stats::score_matrix::ScoreMatrix;

/// Configuration for clustering.
pub struct ClusterConfig {
    pub database: String,
    pub output: String,
    pub threads: i32,
    pub member_cover: f64,
    pub approx_id: f64,
    pub sensitivity: Sensitivity,
}

/// Run greedy incremental clustering.
///
/// Algorithm:
/// 1. Sort sequences by length (longest first)
/// 2. For each sequence, search against existing centroids
/// 3. If it matches a centroid above the coverage threshold, assign it
/// 4. Otherwise, make it a new centroid
pub fn run(config: &ClusterConfig) -> io::Result<()> {
    let start = Instant::now();

    let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // Load sequences
    let db_path = Path::new(&config.database);
    let mut records = if db_path.extension().is_some_and(|e| e == "dmnd") {
        let (_, recs) = crate::data::dmnd_reader::read_dmnd(db_path)?;
        recs
    } else {
        let fasta_path = if db_path.extension().is_none() {
            db_path.with_extension("faa")
        } else {
            db_path.to_path_buf()
        };
        fasta::read_fasta_file(&fasta_path, SequenceType::AminoAcid)?
    };
    eprintln!("Sequences: {}", records.len());

    // Sort by length (longest first — these become centroids preferentially)
    records.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));

    let reduction = Reduction::default_reduction();
    let shape_codes = sensitivity::get_shape_codes(config.sensitivity);
    let shapes: Vec<Shape> = shape_codes
        .iter()
        .map(|code| Shape::from_code(code, &reduction))
        .collect();

    let cover_threshold = config.member_cover / 100.0;
    let id_threshold = config.approx_id / 100.0;

    // Greedy incremental clustering
    let mut result = ClusterResult::new();
    let mut centroid_indices: Vec<usize> = Vec::new();

    for i in 0..records.len() {
        let query = &records[i].sequence;

        // Search against existing centroids
        let mut best_centroid: Option<usize> = None;
        let mut best_score = 0i32;

        if !centroid_indices.is_empty() {
            // Build seed index from current centroids
            let centroid_seqs: Vec<&[Letter]> = centroid_indices
                .iter()
                .map(|&ci| records[ci].sequence.as_slice())
                .collect();
            let query_seqs: Vec<&[Letter]> = vec![query.as_slice()];

            for shape in &shapes {
                let index = seed_match::build_seed_index(&centroid_seqs, shape, &reduction);
                let matches = seed_match::find_seed_matches(&query_seqs, &index, shape, &reduction);

                // Group by centroid
                let mut hits_by_centroid: HashMap<u32, usize> = HashMap::new();
                for m in &matches {
                    *hits_by_centroid.entry(m.ref_id).or_default() += 1;
                }

                // Check top centroid candidates
                for (&centroid_local_id, &hit_count) in &hits_by_centroid {
                    if hit_count < 2 {
                        continue;
                    }
                    let ci = centroid_indices[centroid_local_id as usize];
                    let target = &records[ci].sequence;

                    // Quick coverage check
                    let shorter = query.len().min(target.len()) as f64;
                    let longer = query.len().max(target.len()) as f64;
                    if shorter / longer < cover_threshold {
                        continue;
                    }

                    let sw = smith_waterman::smith_waterman(query, target, &sm);
                    if sw.score > best_score {
                        let pident = if sw.length > 0 {
                            sw.identities as f64 / sw.length as f64
                        } else {
                            0.0
                        };
                        let qcov = sw.length as f64 / query.len() as f64;

                        if pident >= id_threshold && qcov >= cover_threshold {
                            best_score = sw.score;
                            best_centroid = Some(ci);
                        }
                    }
                }

                if best_centroid.is_some() {
                    break;
                }
            }
        }

        match best_centroid {
            Some(ci) => {
                // Assign to existing cluster
                for cluster in &mut result.clusters {
                    if cluster.centroid == ci as u64 {
                        cluster.members.push(i as u64);
                        break;
                    }
                }
            }
            None => {
                // Create new cluster
                centroid_indices.push(i);
                result.clusters.push(Cluster {
                    centroid: i as u64,
                    members: vec![i as u64],
                });
            }
        }

        if (i + 1) % 100 == 0 {
            eprintln!(
                "Processed {}/{} sequences, {} clusters",
                i + 1,
                records.len(),
                result.num_clusters()
            );
        }
    }

    eprintln!(
        "Clustering complete: {} sequences → {} clusters in {:.1}s",
        records.len(),
        result.num_clusters(),
        start.elapsed().as_secs_f64()
    );

    // Write output
    let file = std::fs::File::create(&config.output)?;
    let mut writer = BufWriter::new(file);
    let ids: Vec<String> = records.iter().map(|r| r.id.clone()).collect();
    crate::cluster::write_clusters(&mut writer, &result, &ids)?;
    writer.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cluster_self() {
        let db = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/1.faa");
        let out = std::env::temp_dir().join("test_cluster.out");
        let config = ClusterConfig {
            database: db.to_string(),
            output: out.to_string_lossy().to_string(),
            threads: 1,
            member_cover: 80.0,
            approx_id: 0.0,
            sensitivity: Sensitivity::Default,
        };
        let result = run(&config);
        assert!(result.is_ok(), "Cluster failed: {:?}", result.err());

        let output = std::fs::read_to_string(&out).unwrap();
        assert!(!output.is_empty());
        let _ = std::fs::remove_file(&out);
    }
}
