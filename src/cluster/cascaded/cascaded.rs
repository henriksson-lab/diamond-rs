//! Cascaded clustering orchestration from `cluster/cascaded/cascaded.cpp`.
//!
//! The C++ implementation obtains its database, search engine, logging, and
//! mutable options from global state.  This module keeps the same mutations
//! explicit through [`CascadedConfig`] and [`CascadedBackend`].

use crate::basic::value::SuperBlockId;
use crate::config::{block_size, Sensitivity};
use crate::data::sequence_file::DbFilter;
use crate::util::algo::{greedy_vertex_cover, Edge};
use crate::util::data_structures::FlatArray;

use super::helpers::{
    cluster_steps_with_config, default_round_approx_id, default_round_cov, is_linclust,
    round_ccd_with_config, CallbackBidirectional, CallbackUnidirectional, CascadedHelpersConfig,
    EdgeCallback,
};
use super::Cascaded;

pub const DEFAULT_MEMORY_LIMIT: &str = "16G";
pub const CASCADED_ROUND_MAX_EVALUE: f64 = 0.001;

impl Cascaded {
    pub const fn get_description() -> &'static str {
        "Cascaded greedy vertex cover algorithm"
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphAlgo {
    GreedyVertexCover,
    LengthSorted,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CascadedSearchConfig {
    pub command: &'static str,
    pub output_format: &'static str,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub mutual_cover: bool,
    pub double_indexed: bool,
    pub max_target_seqs: i64,
    pub self_search: bool,
    pub iterate: Vec<String>,
    pub map_any: bool,
    pub lin_stage1_target: bool,
    pub db_size: u64,
    pub chunk_size: f64,
    pub lowmem: i32,
    pub sensitivity: Sensitivity,
    pub lin_stage1_query: bool,
    pub current_ref_block: i32,
    pub kmer_ranking: bool,
    pub min_length_ratio: f64,
    pub global_ranking_targets: bool,
    pub lin_stage1_combo: bool,
    pub approx_min_id: f64,
    pub max_evalue: f64,
    pub comp_based_stats: i32,
    pub hamming_ext: bool,
    pub diag_filter_cov: Option<f64>,
    pub diag_filter_id: Option<f64>,
    pub threads: i32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CascadedConfig {
    pub mutual_cover: Option<f64>,
    pub member_cover: f64,
    pub round_coverage: Vec<String>,
    pub round_approx_id: Vec<String>,
    pub helpers: CascadedHelpersConfig,
    pub approx_min_id: f64,
    pub max_evalue: f64,
    pub anchored_swipe: bool,
    pub extension_mode: String,
    pub comp_based_stats: i32,
    pub hamming_ext: bool,
    pub diag_filter_cov: Option<f64>,
    pub diag_filter_id: Option<f64>,
    pub memory_limit: u64,
    pub threads: i32,
    pub sensitivity: Sensitivity,
    pub lin_stage1_query: bool,
    pub graph_algo: GraphAlgo,
    pub weighted_gvc: bool,
    pub strict_gvc: bool,
    pub no_gvc_reassign: bool,
    pub alignment_output: Option<String>,
    // Fields mutated by `cluster`, retained because they are observable global
    // configuration in the original implementation.
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub db_size: u64,
    pub chunk_size: f64,
    pub lowmem: i32,
}

impl Default for CascadedConfig {
    fn default() -> Self {
        Self {
            mutual_cover: None,
            member_cover: 80.0,
            round_coverage: Vec::new(),
            round_approx_id: Vec::new(),
            helpers: CascadedHelpersConfig::default(),
            approx_min_id: 50.0,
            max_evalue: 0.001,
            anchored_swipe: false,
            extension_mode: String::new(),
            comp_based_stats: 1,
            hamming_ext: false,
            diag_filter_cov: None,
            diag_filter_id: None,
            memory_limit: 16_000_000_000,
            threads: 1,
            sensitivity: Sensitivity::Default,
            lin_stage1_query: false,
            graph_algo: GraphAlgo::GreedyVertexCover,
            weighted_gvc: false,
            strict_gvc: false,
            no_gvc_reassign: false,
            alignment_output: None,
            query_cover: 0.0,
            subject_cover: 0.0,
            query_or_target_cover: 0.0,
            db_size: 0,
            chunk_size: 0.0,
            lowmem: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CascadedRoundSummary {
    pub round: usize,
    pub input_sequences: u64,
    pub clusters: u64,
    pub letters: u64,
    pub edge_count: i64,
}

pub trait CascadedBackend {
    fn sequence_count(&self) -> u64;
    fn letters(&self) -> u64;
    fn letters_filtered(&self, filter: &DbFilter) -> u64;
    fn reset_statistics(&mut self);
    fn run_search(
        &mut self,
        config: &CascadedSearchConfig,
        filter: Option<&DbFilter>,
        callback: &mut dyn EdgeCallback,
    ) -> Result<(), String>;
    fn output_edges(&mut self, _path: &str, _edges: &[Edge<SuperBlockId>]) -> Result<(), String> {
        Ok(())
    }
    fn round_complete(&mut self, _summary: &CascadedRoundSummary) {}
}

/// C++ file-static `rep_bitset`.
pub fn rep_bitset(centroids: &[SuperBlockId], superset: Option<&DbFilter>) -> DbFilter {
    let mut result = DbFilter::new(centroids.len());
    for &centroid in centroids {
        let index = centroid as usize;
        if index < result.oid_filter.len()
            && superset.is_none_or(|filter| filter.get(centroid as u64))
        {
            result.oid_filter[index] = true;
        }
    }
    result
}

fn round_value(
    values: &[String],
    name: &str,
    round: usize,
    round_count: usize,
) -> Result<f64, String> {
    if values.is_empty() || round >= round_count - 1 {
        return Ok(0.0);
    }
    if values.len() >= round_count {
        return Err(format!("Too many values provided for {name}"));
    }
    let mut parsed = Vec::with_capacity(values.len());
    for value in values {
        let stripped = value.trim_start();
        if stripped.is_empty() || stripped.trim_end().len() != stripped.len() {
            return Err(format!("Invalid value provided for {name}: {value}"));
        }
        parsed.push(
            stripped
                .parse::<f64>()
                .map_err(|_| format!("Invalid value provided for {name}: {value}"))?,
        );
    }
    parsed.splice(
        0..0,
        std::iter::repeat_n(parsed[0], round_count - 1 - parsed.len()),
    );
    Ok(parsed[round])
}

fn sensitivity_from_name(name: &str) -> Result<Sensitivity, String> {
    match name {
        "faster" => Ok(Sensitivity::Faster),
        "fast" => Ok(Sensitivity::Fast),
        "default" => Ok(Sensitivity::Default),
        "linclust-40" => Ok(Sensitivity::Linclust40),
        "linclust-20" => Ok(Sensitivity::Linclust20),
        "shapes-6x10" => Ok(Sensitivity::Shapes6x10),
        "shapes-30x10" => Ok(Sensitivity::Shapes30x10),
        "mid-sensitive" => Ok(Sensitivity::MidSensitive),
        "sensitive" => Ok(Sensitivity::Sensitive),
        "more-sensitive" => Ok(Sensitivity::MoreSensitive),
        "very-sensitive" => Ok(Sensitivity::VerySensitive),
        "ultra-sensitive" => Ok(Sensitivity::UltraSensitive),
        _ => Err(format!("Invalid sensitivity level: {name}")),
    }
}

/// `Search::Config::Config()`'s protein-search `min_length_ratio` derivation.
fn derived_min_length_ratio(
    query_cover: f64,
    subject_cover: f64,
    lin_stage1_query: bool,
    sensitivity: Sensitivity,
) -> f64 {
    if query_cover >= 50.0 && query_cover == subject_cover {
        if lin_stage1_query && sensitivity < Sensitivity::Linclust40 {
            (query_cover / 100.0 + 0.05).min(0.92)
        } else {
            (query_cover / 100.0 - 0.05).max(0.0)
        }
    } else {
        0.0
    }
}

fn member_counts(mapping: &[SuperBlockId]) -> Result<Vec<SuperBlockId>, String> {
    let mut counts = vec![0u32; mapping.len()];
    for &centroid in mapping {
        let count = counts
            .get_mut(centroid as usize)
            .ok_or_else(|| format!("Centroid OID out of range: {centroid}"))?;
        *count = (*count)
            .checked_add(1)
            .ok_or_else(|| "Cluster member count overflow".to_owned())?;
    }
    Ok(counts)
}

fn len_sorted_clust(edges: &FlatArray<Edge<SuperBlockId>, SuperBlockId>) -> Vec<SuperBlockId> {
    let mut result = vec![SuperBlockId::MAX; edges.size() as usize];
    for i in 0..edges.size() {
        if result[i as usize] != SuperBlockId::MAX {
            continue;
        }
        result[i as usize] = i;
        for edge in edges.range(i) {
            if result[edge.node2 as usize] == SuperBlockId::MAX {
                result[edge.node2 as usize] = i;
            }
        }
    }
    result
}

fn make_edge_array(
    mut edges: Vec<Edge<SuperBlockId>>,
    sequence_count: SuperBlockId,
) -> FlatArray<Edge<SuperBlockId>, SuperBlockId> {
    edges.sort_unstable();
    let mut result = FlatArray::new();
    let mut begin = 0;
    for node in 0..sequence_count {
        let mut end = begin;
        while end < edges.len() && edges[end].node1 == node {
            end += 1;
        }
        result.push_back(&edges[begin..end]);
        begin = end;
    }
    result
}

/// One C++ `cluster` round.
pub fn cluster<B: CascadedBackend>(
    backend: &mut B,
    config: &mut CascadedConfig,
    filter: Option<&DbFilter>,
    member_count: Option<&[SuperBlockId]>,
    round: usize,
    round_count: usize,
) -> Result<(Vec<SuperBlockId>, i64), String> {
    backend.reset_statistics();
    let mutual = config.mutual_cover.is_some();
    let round_coverage = if config.round_coverage.is_empty() {
        default_round_cov(round_count as i32)
    } else {
        config.round_coverage.clone()
    };
    let coverage = config.mutual_cover.unwrap_or(config.member_cover);
    let round_coverage = coverage.max(round_value(
        &round_coverage,
        "--round-coverage",
        round,
        round_count,
    )?);
    if mutual {
        config.query_cover = round_coverage;
        config.subject_cover = round_coverage;
    } else {
        config.query_cover = 0.0;
        config.subject_cover = 0.0;
        config.query_or_target_cover = round_coverage;
    }
    if let Some(filter) = filter {
        config.db_size = backend.letters_filtered(filter);
    }
    (config.chunk_size, config.lowmem) = if config.lin_stage1_query && round == 0 {
        (32768.0, 1)
    } else {
        let letters = if filter.is_some() {
            config.db_size
        } else {
            backend.letters()
        };
        block_size(
            config.memory_limit.min(i64::MAX as u64) as i64,
            letters.min(i64::MAX as u64) as i64,
            config.sensitivity,
            config.lin_stage1_query,
            config.threads,
        )
    };
    let search_config = CascadedSearchConfig {
        command: "blastp",
        output_format: "edge",
        query_cover: config.query_cover,
        subject_cover: config.subject_cover,
        query_or_target_cover: config.query_or_target_cover,
        mutual_cover: config.mutual_cover.is_some(),
        double_indexed: true,
        max_target_seqs: i64::MAX,
        self_search: true,
        iterate: Vec::new(),
        map_any: false,
        lin_stage1_target: false,
        db_size: config.db_size,
        chunk_size: config.chunk_size,
        lowmem: config.lowmem,
        sensitivity: config.sensitivity,
        lin_stage1_query: config.lin_stage1_query,
        // `Search::keep_target_id` observes this value before the reference
        // loop. Upstream leaves the constructor field uninitialized; the
        // cascaded call sequence deterministically presents a nonzero value,
        // then 0, then nonzero values for this workflow. Preserve that
        // observed search state explicitly instead of relying on UB.
        current_ref_block: if round == 1 { 0 } else { 1 },
        kmer_ranking: false,
        min_length_ratio: derived_min_length_ratio(
            config.query_cover,
            config.subject_cover,
            config.lin_stage1_query,
            config.sensitivity,
        ),
        global_ranking_targets: false,
        lin_stage1_combo: false,
        approx_min_id: config.approx_min_id,
        max_evalue: config.max_evalue,
        comp_based_stats: config.comp_based_stats,
        hamming_ext: config.hamming_ext,
        diag_filter_cov: config.diag_filter_cov,
        diag_filter_id: config.diag_filter_id,
        threads: config.threads,
    };
    let mut callback: Box<dyn EdgeCallback> = if mutual {
        Box::<CallbackBidirectional>::default()
    } else {
        Box::new(CallbackUnidirectional::new(config.member_cover))
    };
    backend.run_search(&search_config, filter, callback.as_mut())?;
    let edge_count = callback.count();
    let edges: Vec<Edge<SuperBlockId>> = callback
        .edges()
        .iter()
        .map(|edge| {
            Ok(Edge::new(
                SuperBlockId::try_from(edge.node1)
                    .map_err(|_| format!("Edge OID out of range: {}", edge.node1))?,
                SuperBlockId::try_from(edge.node2)
                    .map_err(|_| format!("Edge OID out of range: {}", edge.node2))?,
                edge.weight,
            ))
        })
        .collect::<Result<_, String>>()?;
    if let Some(path) = config.alignment_output.as_deref() {
        backend.output_edges(path, &edges)?;
    }
    let sequence_count = SuperBlockId::try_from(backend.sequence_count())
        .map_err(|_| "Sequence count exceeds SuperBlockId".to_owned())?;
    if let Some(edge) = edges
        .iter()
        .find(|edge| edge.node1 >= sequence_count || edge.node2 >= sequence_count)
    {
        return Err(format!(
            "Edge OID out of database range: {} -> {}",
            edge.node1, edge.node2
        ));
    }
    let mut edge_array = make_edge_array(edges, sequence_count);
    let ccd = round_ccd_with_config(
        round as i32,
        round_count as i32,
        config.lin_stage1_query,
        &config.helpers,
    )? as SuperBlockId;
    let mapping = match config.graph_algo {
        GraphAlgo::GreedyVertexCover => greedy_vertex_cover(
            &mut edge_array,
            config.weighted_gvc.then_some(member_count).flatten(),
            !config.strict_gvc,
            !config.no_gvc_reassign,
            ccd,
        ),
        GraphAlgo::LengthSorted => len_sorted_clust(&edge_array),
    };
    Ok((mapping, edge_count))
}

/// C++ file-static `update_clustering`.
pub fn update_clustering(
    previous_filter: &DbFilter,
    previous_centroids: &[SuperBlockId],
    mut current_centroids: Vec<SuperBlockId>,
    round: usize,
) -> Result<(Vec<SuperBlockId>, DbFilter), String> {
    let filter = rep_bitset(&current_centroids, (round > 0).then_some(previous_filter));
    if round > 0 {
        if previous_centroids.len() != current_centroids.len() {
            return Err("Cascaded centroid mappings have different sizes".to_owned());
        }
        for i in 0..current_centroids.len() {
            if !previous_filter.get(i as u64) {
                let previous = previous_centroids[i] as usize;
                let mapped = *current_centroids
                    .get(previous)
                    .ok_or_else(|| format!("Previous centroid OID out of range: {previous}"))?;
                current_centroids[i] = mapped;
            }
        }
    }
    Ok((current_centroids, filter))
}

/// Complete C++ `cascaded` workflow.
pub fn cascaded<B: CascadedBackend>(
    backend: &mut B,
    config: &mut CascadedConfig,
    linear: bool,
) -> Result<Vec<SuperBlockId>, String> {
    let sequence_count = backend.sequence_count();
    if sequence_count > SuperBlockId::MAX as u64 {
        return Err(format!(
            "Workflow supports a maximum of {} input sequences.",
            SuperBlockId::MAX
        ));
    }
    let steps = cluster_steps_with_config(config.approx_min_id, linear, &config.helpers);
    let evalue_cutoff = config.max_evalue;
    let target_approx_id = config.approx_min_id;
    let anchored_swipe = config.anchored_swipe;
    let linclust = is_linclust(&steps);
    let mut filter = DbFilter::new(0);
    let mut cluster_count = sequence_count;
    let mut centroids: Vec<_> = (0..sequence_count as SuperBlockId).collect();
    if linclust {
        config.comp_based_stats = 0;
    }
    for (round, step) in steps.iter().enumerate() {
        config.lin_stage1_query = step.ends_with("_lin");
        config.anchored_swipe = anchored_swipe && (linclust || !config.lin_stage1_query);
        // This intentionally tests the saved value, as the C++ source does.
        if anchored_swipe {
            config.extension_mode = "banded-fast".to_owned();
        }
        let sensitivity = step.strip_suffix("_lin").unwrap_or(step);
        config.sensitivity = sensitivity_from_name(sensitivity)?;
        let round_approx = if config.round_approx_id.is_empty() {
            default_round_approx_id(steps.len() as i32)
        } else {
            config.round_approx_id.clone()
        };
        config.approx_min_id = target_approx_id.max(round_value(
            &round_approx,
            "--round-approx-id",
            round,
            steps.len(),
        )?);
        config.max_evalue = if round + 1 == steps.len() {
            evalue_cutoff
        } else {
            evalue_cutoff.min(CASCADED_ROUND_MAX_EVALUE)
        };
        let counts = config
            .weighted_gvc
            .then(|| member_counts(&centroids))
            .transpose()?;
        let (current, edge_count) = cluster(
            backend,
            config,
            (round > 0).then_some(&filter),
            counts.as_deref(),
            round,
            steps.len(),
        )?;
        (centroids, filter) = update_clustering(&filter, &centroids, current, round)?;
        let clusters = filter.oid_filter.iter().filter(|&&set| set).count() as u64;
        let summary = CascadedRoundSummary {
            round,
            input_sequences: cluster_count,
            clusters,
            letters: backend.letters_filtered(&filter),
            edge_count,
        };
        backend.round_complete(&summary);
        cluster_count = clusters;
    }
    Ok(centroids)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::edge::EdgeData;

    #[derive(Default)]
    struct FakeBackend {
        count: u64,
        lengths: Vec<u64>,
        rounds: Vec<Vec<EdgeData>>,
        searches: Vec<(CascadedSearchConfig, Option<Vec<bool>>)>,
        summaries: Vec<CascadedRoundSummary>,
        resets: usize,
        output: Vec<Edge<SuperBlockId>>,
    }

    fn edge_bytes(edges: &[EdgeData]) -> Vec<u8> {
        let mut bytes = Vec::new();
        for edge in edges {
            bytes.extend_from_slice(&edge.query.to_ne_bytes());
            bytes.extend_from_slice(&edge.target.to_ne_bytes());
            bytes.extend_from_slice(&edge.qcovhsp.to_ne_bytes());
            bytes.extend_from_slice(&edge.scovhsp.to_ne_bytes());
            bytes.extend_from_slice(&edge.evalue.to_ne_bytes());
        }
        bytes
    }

    impl CascadedBackend for FakeBackend {
        fn sequence_count(&self) -> u64 {
            self.count
        }
        fn letters(&self) -> u64 {
            self.lengths.iter().sum()
        }
        fn letters_filtered(&self, filter: &DbFilter) -> u64 {
            filter
                .oid_filter
                .iter()
                .zip(&self.lengths)
                .filter_map(|(&set, &n)| set.then_some(n))
                .sum()
        }
        fn reset_statistics(&mut self) {
            self.resets += 1;
        }
        fn run_search(
            &mut self,
            config: &CascadedSearchConfig,
            filter: Option<&DbFilter>,
            callback: &mut dyn EdgeCallback,
        ) -> Result<(), String> {
            self.searches
                .push((config.clone(), filter.map(|f| f.oid_filter.clone())));
            let edges = self.rounds.remove(0);
            callback.consume(&edge_bytes(&edges))
        }
        fn output_edges(
            &mut self,
            _path: &str,
            edges: &[Edge<SuperBlockId>],
        ) -> Result<(), String> {
            self.output = edges.to_vec();
            Ok(())
        }
        fn round_complete(&mut self, summary: &CascadedRoundSummary) {
            self.summaries.push(summary.clone());
        }
    }

    fn data(query: u64, target: u64, qcov: f32, scov: f32, evalue: f64) -> EdgeData {
        EdgeData {
            query,
            target,
            qcovhsp: qcov,
            scovhsp: scov,
            evalue,
        }
    }

    #[test]
    fn description_and_rep_bitset_match_cpp() {
        assert_eq!(
            Cascaded::get_description(),
            "Cascaded greedy vertex cover algorithm"
        );
        let all = rep_bitset(&[2, 2, 0, 3], None);
        assert_eq!(all.oid_filter, vec![true, false, true, true]);
        let superset = DbFilter {
            oid_filter: vec![true, false, false, true],
            letter_count: 0,
        };
        assert_eq!(
            rep_bitset(&[2, 0, 3, 3], Some(&superset)).oid_filter,
            vec![true, false, false, true]
        );
    }

    #[test]
    fn min_length_ratio_matches_search_config_derivation() {
        assert!(
            (derived_min_length_ratio(80.0, 80.0, true, Sensitivity::Fast) - 0.85).abs()
                < f64::EPSILON
        );
        assert_eq!(
            derived_min_length_ratio(90.0, 90.0, true, Sensitivity::Fast),
            0.92
        );
        assert_eq!(
            derived_min_length_ratio(80.0, 80.0, true, Sensitivity::Linclust40),
            0.75
        );
        assert_eq!(
            derived_min_length_ratio(80.0, 70.0, true, Sensitivity::Fast),
            0.0
        );
    }

    #[test]
    fn update_clustering_composes_previous_mapping() {
        let previous = DbFilter {
            oid_filter: vec![true, false, true, false],
            letter_count: 0,
        };
        let (mapping, filter) =
            update_clustering(&previous, &[0, 0, 2, 2], vec![0, 1, 0, 3], 1).unwrap();
        assert_eq!(mapping, vec![0, 0, 0, 0]);
        assert_eq!(filter.oid_filter, vec![true, false, false, false]);
    }

    #[test]
    fn cluster_selects_unidirectional_callback_and_length_order() {
        let mut backend = FakeBackend {
            count: 3,
            lengths: vec![10, 20, 30],
            rounds: vec![vec![
                data(1, 0, 90.0, 20.0, 2.0),
                data(2, 1, 10.0, 90.0, 1.0),
            ]],
            ..Default::default()
        };
        let mut config = CascadedConfig {
            graph_algo: GraphAlgo::LengthSorted,
            alignment_output: Some("edges.tsv".into()),
            hamming_ext: true,
            diag_filter_cov: Some(70.0),
            diag_filter_id: Some(80.0),
            ..Default::default()
        };
        let (mapping, count) = cluster(&mut backend, &mut config, None, None, 0, 1).unwrap();
        assert_eq!(count, 2);
        assert_eq!(
            backend.output,
            vec![Edge::new(0, 1, 2.0), Edge::new(2, 1, 1.0)]
        );
        assert_eq!(mapping, vec![0, 0, 2]);
        let search = &backend.searches[0].0;
        assert_eq!((search.command, search.output_format), ("blastp", "edge"));
        assert!(search.double_indexed && search.self_search);
        assert!(search.hamming_ext);
        assert!(!search.mutual_cover);
        assert_eq!(search.diag_filter_cov, Some(70.0));
        assert_eq!(search.diag_filter_id, Some(80.0));
        assert_eq!(
            (
                search.query_cover,
                search.subject_cover,
                search.query_or_target_cover
            ),
            (0.0, 0.0, 80.0)
        );
    }

    #[test]
    fn cascaded_runs_rounds_filters_and_restores_final_cutoffs() {
        let mut backend = FakeBackend {
            count: 3,
            lengths: vec![10, 20, 30],
            rounds: vec![vec![data(1, 0, 90.0, 90.0, 1e-5)], vec![]],
            ..Default::default()
        };
        let mut config = CascadedConfig {
            mutual_cover: Some(75.0),
            helpers: CascadedHelpersConfig {
                cluster_steps: vec!["faster_lin".into(), "default".into()],
                connected_component_depth: vec![],
            },
            max_evalue: 0.1,
            approx_min_id: 50.0,
            anchored_swipe: true,
            graph_algo: GraphAlgo::LengthSorted,
            ..Default::default()
        };
        let mapping = cascaded(&mut backend, &mut config, false).unwrap();
        assert_eq!(mapping, vec![0, 0, 2]);
        assert_eq!(backend.resets, 2);
        assert_eq!(backend.searches[0].1, None);
        assert_eq!(backend.searches[1].1, Some(vec![true, false, true]));
        assert_eq!(backend.searches[0].0.chunk_size, 32768.0);
        assert!(backend.searches[0].0.mutual_cover);
        assert_eq!(
            (
                backend.searches[0].0.query_cover,
                backend.searches[0].0.subject_cover
            ),
            (75.0, 75.0)
        );
        assert_eq!(
            backend
                .summaries
                .iter()
                .map(|s| s.clusters)
                .collect::<Vec<_>>(),
            vec![2, 2]
        );
        assert_eq!(
            backend
                .summaries
                .iter()
                .map(|s| s.edge_count)
                .collect::<Vec<_>>(),
            vec![2, 0]
        );
        assert_eq!(config.max_evalue, 0.1);
        assert_eq!(config.extension_mode, "banded-fast");
        assert!(config.anchored_swipe);
    }

    #[test]
    fn round_value_has_cpp_count_and_parse_rules() {
        assert_eq!(round_value(&["70".into()], "--x", 0, 3).unwrap(), 70.0);
        assert_eq!(round_value(&["70".into()], "--x", 2, 3).unwrap(), 0.0);
        assert!(round_value(&["70 ".into()], "--x", 0, 3).is_err());
        assert!(round_value(&["1".into(), "2".into()], "--x", 0, 2).is_err());
    }

    #[test]
    fn cascaded_rejects_sequence_counts_that_do_not_fit_superblock_ids() {
        let mut backend = FakeBackend {
            count: SuperBlockId::MAX as u64 + 1,
            ..Default::default()
        };
        let error = cascaded(&mut backend, &mut CascadedConfig::default(), false).unwrap_err();
        assert_eq!(
            error,
            format!(
                "Workflow supports a maximum of {} input sequences.",
                SuperBlockId::MAX
            )
        );
        assert_eq!(backend.resets, 0);
    }
}
