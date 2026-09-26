//! Command-level greedy vertex cover.
//!
//! This mirrors `diamond/src/tools/greedy_vertex_cover.cpp`. The graph
//! algorithm itself remains in [`crate::util::algo::greedy_vertex_cover`];
//! this module implements the source file's mapping/edge parsing, directional
//! coverage rules, graph construction, and output generation.

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use crate::util::algo::{greedy_vertex_cover as compute_vertex_cover, Edge};
use crate::util::data_structures::FlatArray;

/// Default used by `Cluster::DEFAULT_MEMBER_COVER` in the C++ command.
pub const DEFAULT_MEMBER_COVER: f64 = 80.0;

/// Explicit replacement for the process-wide C++ `config` object.
#[derive(Debug, Clone, PartialEq)]
pub struct GreedyVertexCoverConfig {
    pub database: PathBuf,
    pub edges: PathBuf,
    pub query_or_target_cover: f64,
    pub member_cover: Option<f64>,
    pub edge_format: String,
    pub symmetric: bool,
    pub strict_gvc: bool,
    pub no_gvc_reassign: bool,
    pub connected_component_depth: Vec<String>,
    pub centroid_out: Option<PathBuf>,
    pub output_file: Option<PathBuf>,
}

impl GreedyVertexCoverConfig {
    fn coverage_cutoff(&self) -> f64 {
        let member_cover = self.member_cover.unwrap_or(DEFAULT_MEMBER_COVER);
        // Match `std::max(a, b)` including its left-operand NaN behavior.
        if self.query_or_target_cover < member_cover {
            member_cover
        } else {
            self.query_or_target_cover
        }
    }

    fn connected_component_depth(&self) -> u64 {
        self.connected_component_depth
            .first()
            .map_or(0, |value| c_atoi(value) as u64)
    }
}

/// Generated command output, retained independently of optional output files.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GreedyVertexCoverOutput {
    pub assignments: Vec<(String, String)>,
    pub centroids: Vec<String>,
}

impl GreedyVertexCoverOutput {
    pub fn centroid_count(&self) -> usize {
        self.centroids.len()
    }

    fn assignments_text(&self) -> String {
        let mut output = String::new();
        for (centroid, member) in &self.assignments {
            output.push_str(centroid);
            output.push('\t');
            output.push_str(member);
            output.push('\n');
        }
        output
    }

    fn centroids_text(&self) -> String {
        let mut output = self.centroids.join("\n");
        if !output.is_empty() {
            output.push('\n');
        }
        output
    }
}

/// Translate C `atoi`'s command-line behavior.
fn c_atoi(value: &str) -> i32 {
    let value = value.trim_start();
    let (negative, digits) = match value.as_bytes().first() {
        Some(b'-') => (true, &value[1..]),
        Some(b'+') => (false, &value[1..]),
        _ => (false, value),
    };
    let digits = &digits[..digits
        .find(|character: char| !character.is_ascii_digit())
        .unwrap_or(digits.len())];
    let magnitude = digits.parse::<i64>().unwrap_or(0);
    let signed = if negative { -magnitude } else { magnitude };
    signed as i32
}

fn source_lines(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .map(|line| line.strip_suffix('\r').unwrap_or(line))
}

fn parse_number(value: Option<&str>, field: &str, line: usize) -> Result<f64, String> {
    value
        .ok_or_else(|| format!("Missing {field} on edge line {line}."))?
        .parse::<f64>()
        .map_err(|_| format!("Invalid {field} on edge line {line}."))
}

fn make_edge_array(mut edges: Vec<Edge<u64>>, node_count: usize) -> FlatArray<Edge<u64>, u64> {
    edges.sort_unstable_by(|left, right| {
        left.node1
            .cmp(&right.node1)
            .then_with(|| left.node2.cmp(&right.node2))
    });
    let mut limits = Vec::with_capacity(node_count + 1);
    limits.push(0);
    let mut cursor = 0usize;
    for node in 0..node_count {
        while cursor < edges.len() && edges[cursor].node1 == node as u64 {
            cursor += 1;
        }
        limits.push(cursor as u64);
    }
    FlatArray::from_limits_data(limits, edges)
}

/// Run the source workflow on in-memory text. This is the deterministic core
/// used by the file-backed command and focused tests.
pub fn greedy_vertex_cover(
    mapping_text: &str,
    edges_text: &str,
    config: &GreedyVertexCoverConfig,
) -> Result<GreedyVertexCoverOutput, String> {
    let triplets = config.edge_format == "triplet";
    if !triplets && config.symmetric {
        return Err("--symmetric requires triplet edge format".to_owned());
    }

    let mut accession_to_oid = HashMap::<String, u64>::new();
    for line in source_lines(mapping_text) {
        let accession = line.split('\t').next().unwrap_or_default();
        let next_oid = accession_to_oid.len() as u64;
        accession_to_oid
            .entry(accession.to_owned())
            .or_insert(next_oid);
    }

    let coverage = config.coverage_cutoff();
    let mut edges = Vec::<Edge<u64>>::new();
    for (line_index, line) in source_lines(edges_text).enumerate() {
        let line_number = line_index + 1;
        let mut fields = line.split('\t');
        let query = fields
            .next()
            .ok_or_else(|| format!("Missing query on edge line {line_number}."))?;
        let target = fields
            .next()
            .ok_or_else(|| format!("Missing target on edge line {line_number}."))?;
        let (query_coverage, target_coverage) = if triplets {
            (0.0, 0.0)
        } else {
            (
                parse_number(fields.next(), "query coverage", line_number)?,
                parse_number(fields.next(), "target coverage", line_number)?,
            )
        };
        let evalue = parse_number(fields.next(), "e-value", line_number)?;

        if triplets || target_coverage >= coverage || query_coverage >= coverage {
            let query_oid = *accession_to_oid.get(query).ok_or_else(|| {
                format!("Unknown query accession `{query}` on edge line {line_number}.")
            })?;
            let target_oid = *accession_to_oid.get(target).ok_or_else(|| {
                format!("Unknown target accession `{target}` on edge line {line_number}.")
            })?;
            if query_oid == target_oid {
                continue;
            }
            if triplets {
                edges.push(Edge {
                    node1: target_oid,
                    node2: query_oid,
                    weight: evalue,
                });
                if config.symmetric {
                    edges.push(Edge {
                        node1: query_oid,
                        node2: target_oid,
                        weight: evalue,
                    });
                }
            } else {
                if target_coverage >= coverage {
                    edges.push(Edge {
                        node1: query_oid,
                        node2: target_oid,
                        weight: evalue,
                    });
                }
                if query_coverage >= coverage {
                    edges.push(Edge {
                        node1: target_oid,
                        node2: query_oid,
                        weight: evalue,
                    });
                }
            }
        }
    }

    let mut edge_array = make_edge_array(edges, accession_to_oid.len());
    let representatives = compute_vertex_cover(
        &mut edge_array,
        None,
        !config.strict_gvc,
        !config.no_gvc_reassign,
        config.connected_component_depth(),
    );

    let mut accessions = vec![String::new(); accession_to_oid.len()];
    for (accession, oid) in accession_to_oid {
        accessions[oid as usize] = accession;
    }
    let mut assignments = Vec::with_capacity(accessions.len());
    let mut centroids = Vec::new();
    for (member_oid, &representative_oid) in representatives.iter().enumerate() {
        if representative_oid as usize == member_oid {
            centroids.push(accessions[member_oid].clone());
        }
        assignments.push((
            accessions[representative_oid as usize].clone(),
            accessions[member_oid].clone(),
        ));
    }
    Ok(GreedyVertexCoverOutput {
        assignments,
        centroids,
    })
}

/// File-backed adapter for [`greedy_vertex_cover`].
pub fn greedy_vertex_cover_from_files(
    config: &GreedyVertexCoverConfig,
) -> Result<GreedyVertexCoverOutput, String> {
    let mapping_text = fs::read_to_string(&config.database).map_err(|error| {
        format!(
            "Failed to read mapping file {}: {error}",
            config.database.display()
        )
    })?;
    let edges_text = fs::read_to_string(&config.edges).map_err(|error| {
        format!(
            "Failed to read edge file {}: {error}",
            config.edges.display()
        )
    })?;
    let output = greedy_vertex_cover(&mapping_text, &edges_text, config)?;
    if let Some(path) = &config.centroid_out {
        fs::write(path, output.centroids_text()).map_err(|error| {
            format!("Failed to write centroid file {}: {error}", path.display())
        })?;
    }
    if let Some(path) = &config.output_file {
        fs::write(path, output.assignments_text())
            .map_err(|error| format!("Failed to write output file {}: {error}", path.display()))?;
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> GreedyVertexCoverConfig {
        GreedyVertexCoverConfig {
            database: PathBuf::new(),
            edges: PathBuf::new(),
            query_or_target_cover: 0.0,
            member_cover: None,
            edge_format: "default".to_owned(),
            symmetric: false,
            strict_gvc: false,
            no_gvc_reassign: false,
            connected_component_depth: Vec::new(),
            centroid_out: None,
            output_file: None,
        }
    }

    #[test]
    fn symmetric_triplets_make_one_star_cluster() {
        let mut config = config();
        config.edge_format = "triplet".to_owned();
        config.symmetric = true;
        let output = greedy_vertex_cover("A\nB\nC\n", "A\tB\t1e-5\nA\tC\t2e-5\n", &config).unwrap();
        assert_eq!(output.centroids, ["A"]);
        assert_eq!(
            output.assignments,
            [
                ("A".into(), "A".into()),
                ("A".into(), "B".into()),
                ("A".into(), "C".into())
            ]
        );
    }

    #[test]
    fn default_edges_follow_both_coverage_directions() {
        let config = config();
        let target_covered =
            greedy_vertex_cover("A\nB\n", "A\tB\t10\t90\t1e-5\n", &config).unwrap();
        assert_eq!(target_covered.centroids, ["A"]);

        let query_covered = greedy_vertex_cover("A\nB\n", "A\tB\t90\t10\t1e-5\n", &config).unwrap();
        assert_eq!(query_covered.centroids, ["B"]);
    }

    #[test]
    fn validates_mode_and_accessions_and_keeps_first_mapping_occurrence() {
        let mut config = config();
        config.symmetric = true;
        assert_eq!(
            greedy_vertex_cover("A\n", "", &config).unwrap_err(),
            "--symmetric requires triplet edge format"
        );

        config.symmetric = false;
        let error = greedy_vertex_cover("A\n", "A\tB\t90\t90\t1\n", &config).unwrap_err();
        assert!(error.contains("Unknown target accession `B`"));

        let output = greedy_vertex_cover("A\tfirst\nA\tduplicate\nB\n", "", &config).unwrap();
        assert_eq!(output.assignments.len(), 2);
        assert_eq!(output.centroids, ["A", "B"]);
    }

    #[test]
    fn file_command_writes_source_ordered_outputs() {
        let root = std::env::temp_dir().join(format!(
            "diamond-gvc-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        fs::create_dir_all(&root).unwrap();
        let mut config = config();
        config.database = root.join("mapping.tsv");
        config.edges = root.join("edges.tsv");
        config.output_file = Some(root.join("clusters.tsv"));
        config.centroid_out = Some(root.join("centroids.txt"));
        fs::write(&config.database, "A\nB\n").unwrap();
        fs::write(&config.edges, "A\tB\t10\t90\t1e-5\n").unwrap();

        let output = greedy_vertex_cover_from_files(&config).unwrap();
        assert_eq!(output.centroid_count(), 1);
        assert_eq!(
            fs::read_to_string(config.output_file.as_ref().unwrap()).unwrap(),
            "A\tA\nA\tB\n"
        );
        assert_eq!(
            fs::read_to_string(config.centroid_out.as_ref().unwrap()).unwrap(),
            "A\n"
        );

        fs::remove_dir_all(root).unwrap();
    }
}
