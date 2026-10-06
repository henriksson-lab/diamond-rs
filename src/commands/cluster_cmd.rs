//! Native clustering command boundary.
//!
//! The structurally translated upstream workflow lives in
//! [`crate::cluster::cascaded`]. Database, search, and output integration use
//! the translated backend traits so each command retains its upstream workflow
//! rather than substituting a shared clustering algorithm.

use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::basic::value::{OId, SequenceType, SuperBlockId};
use crate::cluster::cascaded::cascaded::{self, CascadedBackend, CascadedRoundSummary};
use crate::cluster::cascaded::wrapper::{
    BestCentroid, CascadedWrapperBackend, CascadedWrapperConfig, CentroidSearchConfig,
    LengthSortedBlock, WrapperCommand, WrapperOutput,
};
use crate::cluster::cascaded::CascadedConfig;
use crate::commands::blastp::{self, BlastpConfig};
use crate::data::fasta::{self, FastaRecord};
use crate::data::sequence_file::DbFilter;
use crate::masking::MaskingMode;
use crate::stats::cbs::CbsMode;
use crate::util::algo::Edge;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClusterWorkflow {
    Cascaded,
    LinClust,
    DeepClust,
}

/// CLI configuration retained until a concrete cascaded backend takes
/// ownership of the database and output adapters.
#[derive(Debug, Clone, PartialEq)]
pub struct ClusterConfig {
    pub workflow: ClusterWorkflow,
    pub database: String,
    pub output: String,
    pub threads: i32,
    pub member_cover: f64,
    /// `None` is important: upstream applies command-specific defaults in
    /// `Cluster::init_thresholds` (50, 90, and 0 respectively).
    pub approx_id: Option<f64>,
    pub cluster_steps: Vec<String>,
    pub alignment_output: Option<String>,
}

fn resolved_threads(threads: i32) -> i32 {
    if threads > 0 {
        threads
    } else {
        i32::try_from(
            std::thread::available_parallelism()
                .map(std::num::NonZeroUsize::get)
                .unwrap_or(1),
        )
        .unwrap_or(i32::MAX)
    }
}

/// Map CLI state into the verbatim translated cascaded wrapper.
pub fn translated_config(config: &ClusterConfig) -> CascadedWrapperConfig {
    // basic/config.cpp resolves an unset thread count before clustering setup;
    // block_size/seed partition arithmetic must never observe zero.
    let threads = resolved_threads(config.threads);
    let mut translated = CascadedWrapperConfig {
        database: Some(config.database.clone()),
        command: match config.workflow {
            ClusterWorkflow::Cascaded => WrapperCommand::Cascaded,
            ClusterWorkflow::LinClust => WrapperCommand::LinClust,
            ClusterWorkflow::DeepClust => WrapperCommand::DeepClust,
        },
        threads,
        ..CascadedWrapperConfig::default()
    };
    translated.thresholds.member_cover = Some(config.member_cover);
    translated.thresholds.approx_min_id = config.approx_id;
    translated.core.threads = threads;
    translated.core.helpers.cluster_steps = config.cluster_steps.clone();
    translated.core.alignment_output = config.alignment_output.clone();
    translated
}

#[derive(Clone)]
struct NativeBlock {
    records: Vec<FastaRecord>,
    original_oids: Vec<OId>,
}

/// Apply `Block::length_sorted` at the same linear self-search boundary as
/// `run/double_indexed.cpp`. The upstream comparison is
/// `greater<pair<Loc, BlockId>>`: length descending, then the pre-sort block
/// ID descending for equal lengths.
fn length_sort_linear_active(active: &mut [(OId, FastaRecord)], linear: bool) {
    if !linear {
        return;
    }
    active.sort_by(|(left_oid, left), (right_oid, right)| {
        right
            .sequence
            .len()
            .cmp(&left.sequence.len())
            .then_with(|| right_oid.cmp(left_oid))
    });
}

impl NativeBlock {
    fn letter_count(&self) -> u64 {
        self.records.iter().map(|r| r.sequence.len() as u64).sum()
    }
}

struct NativeDatabase {
    path: PathBuf,
    block: NativeBlock,
}

struct NativeBackend {
    output: PathBuf,
    centroids: NativeBlock,
}

impl NativeBackend {
    fn new(output: &str) -> Self {
        Self {
            output: PathBuf::from(output),
            centroids: NativeBlock {
                records: Vec::new(),
                original_oids: Vec::new(),
            },
        }
    }

    fn load(path: &str) -> Result<NativeDatabase, String> {
        let requested = Path::new(path);
        let (resolved, records) = if requested.exists() {
            if requested.extension().is_some_and(|ext| ext == "dmnd") {
                let (_, records) = crate::data::dmnd_reader::read_dmnd(requested)
                    .map_err(|error| error.to_string())?;
                (requested.to_path_buf(), records)
            } else {
                (
                    requested.to_path_buf(),
                    fasta::read_fasta_file(requested, SequenceType::AminoAcid)
                        .map_err(|error| error.to_string())?,
                )
            }
        } else {
            let mut appended = requested.as_os_str().to_owned();
            appended.push(".dmnd");
            let resolved = PathBuf::from(appended);
            let (_, records) = crate::data::dmnd_reader::read_dmnd(&resolved)
                .map_err(|error| error.to_string())?;
            (resolved, records)
        };
        let original_oids = (0..records.len() as OId).collect();
        Ok(NativeDatabase {
            path: resolved,
            block: NativeBlock {
                records,
                original_oids,
            },
        })
    }
}

impl CascadedBackend for NativeBlock {
    fn sequence_count(&self) -> u64 {
        self.records.len() as u64
    }

    fn letters(&self) -> u64 {
        self.letter_count()
    }

    fn letters_filtered(&self, filter: &DbFilter) -> u64 {
        self.records
            .iter()
            .enumerate()
            .filter(|(oid, _)| filter.get(*oid as OId))
            .map(|(_, record)| record.sequence.len() as u64)
            .sum()
    }

    fn reset_statistics(&mut self) {}

    fn run_search(
        &mut self,
        config: &cascaded::CascadedSearchConfig,
        filter: Option<&DbFilter>,
        callback: &mut dyn crate::cluster::cascaded::helpers::EdgeCallback,
    ) -> Result<(), String> {
        let mut active = self
            .records
            .iter()
            .enumerate()
            .filter(|(oid, _)| filter.is_none_or(|set| set.get(*oid as OId)))
            .map(|(oid, record)| (oid as OId, record.clone()))
            .collect::<Vec<_>>();
        length_sort_linear_active(&mut active, config.lin_stage1_query);
        let local_to_block = active.iter().map(|(oid, _)| *oid).collect::<Vec<_>>();
        let records = active
            .into_iter()
            .map(|(_, record)| record)
            .collect::<Vec<_>>();
        let target_count = records.len() as i64;
        let search = BlastpConfig {
            query_files: Vec::new(),
            database: String::new(),
            output: None,
            matrix: "blosum62".to_owned(),
            gap_open: -1,
            gap_extend: -1,
            max_evalue: config.max_evalue,
            // Upstream uses INT64_MAX to mean every target. The native
            // chunk-size arithmetic cannot add to INT64_MAX; the number of
            // loaded targets is the exactly equivalent finite upper bound.
            max_target_seqs: config.max_target_seqs.min(target_count),
            ext_chunk_size: 0,
            toppercent: None,
            global_ranking_targets: 0,
            // Cascaded search applies `approx_min_id`; exact `min_id` remains
            // unset for this workflow.
            min_id: 0.0,
            threads: config.threads,
            outfmt: vec!["6".to_owned()],
            sensitivity: config.sensitivity,
            masking: MaskingMode::None,
            motif_masking: String::new(),
            min_query_len: 0,
            query_cover: config.query_cover,
            subject_cover: config.subject_cover,
            comp_based_stats: CbsMode::parse(&config.comp_based_stats.to_string()),
            // Upstream `config.self` selects the triangular all-vs-all stage-1
            // kernel; it is independent of the `--no-self-hits` content filter.
            no_self_hits: false,
            ungapped_xdrop_bits: 12.3,
            memory_limit: None,
            tmpdir: PathBuf::new(),
            translated_query_layout: None,
        };
        let mut edges = blastp::run_edges_in_memory(
            &search,
            records.clone(),
            records,
            config.approx_min_id,
            config.lin_stage1_query,
            config.self_search,
            config.query_or_target_cover,
        )
        .map_err(|error| error.to_string())?;
        edges.sort_by_key(|edge| (edge.edge.query, edge.edge.target));
        let mut bytes = Vec::with_capacity(edges.len() * crate::output::edge::EdgeData::SIZE);
        for mut edge in edges {
            edge.edge.query = *local_to_block
                .get(edge.edge.query as usize)
                .ok_or_else(|| "Cluster query OID out of range".to_owned())?;
            edge.edge.target = *local_to_block
                .get(edge.edge.target as usize)
                .ok_or_else(|| "Cluster target OID out of range".to_owned())?;
            edge.edge
                .write(&mut bytes)
                .map_err(|error| error.to_string())?;
        }
        callback.consume(&bytes)
    }

    fn output_edges(&mut self, path: &str, edges: &[Edge<SuperBlockId>]) -> Result<(), String> {
        let file = std::fs::File::create(path).map_err(|error| error.to_string())?;
        let mut out = BufWriter::new(file);
        for edge in edges {
            let query = self
                .records
                .get(edge.node1 as usize)
                .ok_or_else(|| format!("Edge query OID out of range: {}", edge.node1))?;
            let target = self
                .records
                .get(edge.node2 as usize)
                .ok_or_else(|| format!("Edge target OID out of range: {}", edge.node2))?;
            let query_id = query
                .id
                .split_ascii_whitespace()
                .next()
                .unwrap_or(&query.id);
            let target_id = target
                .id
                .split_ascii_whitespace()
                .next()
                .unwrap_or(&target.id);
            writeln!(out, "{query_id}\t{target_id}").map_err(|error| error.to_string())?;
        }
        out.flush().map_err(|error| error.to_string())
    }

    fn round_complete(&mut self, _summary: &CascadedRoundSummary) {}
}

impl CascadedWrapperBackend for NativeBackend {
    type Database = NativeDatabase;
    type Block = NativeBlock;

    fn external(&mut self, _config: &CascadedWrapperConfig) -> Result<(), String> {
        Err("native external linclust backend is not implemented".to_owned())
    }

    fn open_database(&mut self, path: &str) -> Result<Self::Database, String> {
        Self::load(path)
    }

    fn database_name(&self, database: &Self::Database) -> String {
        database.path.to_string_lossy().into_owned()
    }

    fn database_is_blast(&self, _database: &Self::Database) -> bool {
        false
    }

    fn database_sequence_count(&self, database: &Self::Database) -> OId {
        database.block.records.len() as OId
    }

    fn database_letters(&self, database: &Self::Database) -> u64 {
        database.block.letter_count()
    }

    fn length_sort(
        &mut self,
        _database: &mut Self::Database,
        _memory_limit: i64,
        _minimizer_window: i32,
        _sketch_size: i32,
    ) -> Result<Vec<LengthSortedBlock<Self::Block>>, String> {
        Err("native cluster length-sorted streaming adapter is not implemented".to_owned())
    }

    fn block_sequence_count(&self, block: &Self::Block) -> OId {
        block.records.len() as OId
    }

    fn block_letters(&self, block: &Self::Block) -> u64 {
        block.letter_count()
    }

    fn centroid_sequence_count(&self) -> OId {
        self.centroids.records.len() as OId
    }

    fn centroid_letters(&self) -> u64 {
        self.centroids.letter_count()
    }

    fn run_centroid_search(
        &mut self,
        _block: &mut Self::Block,
        _config: &CentroidSearchConfig,
        _consumer: &mut BestCentroid,
    ) -> Result<(), String> {
        Err("native cluster centroid search adapter is not implemented".to_owned())
    }

    fn sub_database(
        &mut self,
        block: &Self::Block,
        ids: &[SuperBlockId],
    ) -> Result<Self::Block, String> {
        let mut records = Vec::with_capacity(ids.len());
        let mut original_oids = Vec::with_capacity(ids.len());
        for &id in ids {
            let index = id as usize;
            records.push(
                block
                    .records
                    .get(index)
                    .ok_or_else(|| format!("Block OID out of range: {id}"))?
                    .clone(),
            );
            original_oids.push(block.original_oids[index]);
        }
        Ok(NativeBlock {
            records,
            original_oids,
        })
    }

    fn cascaded_database(
        &mut self,
        database: &mut Self::Database,
        linear: bool,
        config: &mut CascadedConfig,
    ) -> Result<Vec<SuperBlockId>, String> {
        cascaded::cascaded(&mut database.block, config, linear)
    }

    fn cascaded_block(
        &mut self,
        block: &mut Self::Block,
        linear: bool,
        config: &mut CascadedConfig,
    ) -> Result<Vec<SuperBlockId>, String> {
        cascaded::cascaded(block, config, linear)
    }

    fn append_centroids(
        &mut self,
        block: &Self::Block,
        ids: &[SuperBlockId],
    ) -> Result<(), String> {
        let selected = self.sub_database(block, ids)?;
        self.centroids.records.extend(selected.records);
        self.centroids.original_oids.extend(selected.original_oids);
        Ok(())
    }

    fn close_block(&mut self, _block: Self::Block) -> Result<(), String> {
        Ok(())
    }

    fn write_output(
        &mut self,
        database: &Self::Database,
        output: WrapperOutput,
        oid_output: bool,
        simple_header: bool,
    ) -> Result<(), String> {
        let file = std::fs::File::create(&self.output).map_err(|error| error.to_string())?;
        let mut writer = BufWriter::new(file);
        if simple_header {
            writeln!(writer, "centroid\tmember").map_err(|error| error.to_string())?;
        }
        let mut write_pair = |centroid: OId, member: OId| -> Result<(), String> {
            if oid_output {
                writeln!(writer, "{centroid}\t{member}").map_err(|error| error.to_string())
            } else {
                let centroid_title = &database
                    .block
                    .records
                    .get(centroid as usize)
                    .ok_or_else(|| format!("Centroid OID out of range: {centroid}"))?
                    .id;
                let member_title = &database
                    .block
                    .records
                    .get(member as usize)
                    .ok_or_else(|| format!("Member OID out of range: {member}"))?
                    .id;
                let centroid_id = centroid_title
                    .split_whitespace()
                    .next()
                    .unwrap_or(centroid_title);
                let member_id = member_title
                    .split_whitespace()
                    .next()
                    .unwrap_or(member_title);
                writeln!(writer, "{centroid_id}\t{member_id}").map_err(|error| error.to_string())
            }
        };
        let mut pairs = match output {
            WrapperOutput::Dense(mapping) => mapping
                .into_iter()
                .enumerate()
                .map(|(member, centroid)| (centroid as OId, member as OId))
                .collect::<Vec<_>>(),
            WrapperOutput::Pairs(pairs) => pairs,
        };
        pairs.sort_unstable();
        for (centroid, member) in pairs {
            write_pair(centroid, member)?;
        }
        writer.flush().map_err(|error| error.to_string())
    }

    fn close_database(&mut self, _database: Self::Database) -> Result<(), String> {
        Ok(())
    }
}

/// Run native clustering.
///
/// Database ownership, block orchestration, and output use the translated
/// cascaded workflow above; it never substitutes the former one-pass greedy
/// heuristic.
pub fn run(config: &ClusterConfig) -> io::Result<()> {
    let mut translated = translated_config(config);
    let mut backend = NativeBackend::new(&config.output);
    crate::cluster::cascaded::wrapper::run(&mut backend, &mut translated).map_err(io::Error::other)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::helpers::{init_thresholds, ClusterCommand};

    fn config(workflow: ClusterWorkflow) -> ClusterConfig {
        ClusterConfig {
            workflow,
            database: "input.dmnd".to_owned(),
            output: "clusters.tsv".to_owned(),
            threads: 3,
            member_cover: 80.0,
            approx_id: None,
            cluster_steps: Vec::new(),
            alignment_output: None,
        }
    }

    #[test]
    fn linear_active_block_matches_upstream_length_sort_tie_break() {
        let record = |id: &str, len: usize| FastaRecord {
            id: id.to_owned(),
            sequence: vec![1; len],
        };
        let original = vec![
            (4, record("short", 3)),
            (7, record("equal-low", 5)),
            (9, record("equal-high", 5)),
            (12, record("long", 8)),
        ];

        let mut linear = original.clone();
        length_sort_linear_active(&mut linear, true);
        assert_eq!(
            linear.iter().map(|(oid, _)| *oid).collect::<Vec<_>>(),
            vec![12, 9, 7, 4]
        );

        let mut nonlinear = original.clone();
        length_sort_linear_active(&mut nonlinear, false);
        assert_eq!(
            nonlinear.iter().map(|(oid, _)| *oid).collect::<Vec<_>>(),
            vec![4, 7, 9, 12]
        );
    }

    #[test]
    fn workflows_map_to_distinct_upstream_command_state() {
        for (workflow, wrapper, threshold_command, default_identity) in [
            (
                ClusterWorkflow::Cascaded,
                WrapperCommand::Cascaded,
                ClusterCommand::Other,
                50.0,
            ),
            (
                ClusterWorkflow::LinClust,
                WrapperCommand::LinClust,
                ClusterCommand::LinClust,
                90.0,
            ),
            (
                ClusterWorkflow::DeepClust,
                WrapperCommand::DeepClust,
                ClusterCommand::DeepClust,
                0.0,
            ),
        ] {
            let mut translated = translated_config(&config(workflow));
            assert_eq!(translated.command, wrapper);
            translated.thresholds.command = threshold_command;
            init_thresholds(&mut translated.thresholds).unwrap();
            assert_eq!(translated.thresholds.approx_min_id, Some(default_identity));
        }
    }

    #[test]
    fn explicit_identity_is_preserved_and_native_compact_clustering_runs() {
        let mut cfg = config(ClusterWorkflow::LinClust);
        cfg.approx_id = Some(71.5);
        let translated = translated_config(&cfg);
        assert_eq!(translated.thresholds.approx_min_id, Some(71.5));

        let stem = format!("diamond-cluster-{}", std::process::id());
        let input = std::env::temp_dir().join(format!("{stem}.faa"));
        let output = std::env::temp_dir().join(format!("{stem}.tsv"));
        let sequence = "ARNDCQEGHILKMFPSTWYVARNDCQEGHILKMFPSTWYV";
        std::fs::write(&input, format!(">a\n{sequence}\n>b\n{sequence}\n")).unwrap();
        cfg.database = input.to_string_lossy().into_owned();
        cfg.output = output.to_string_lossy().into_owned();
        run(&cfg).unwrap();
        let rows = std::fs::read_to_string(&output).unwrap();
        assert_eq!(rows.lines().count(), 2);
        let centroids = rows
            .lines()
            .map(|line| line.split_once('\t').unwrap().0)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(centroids.len(), 1);
        let _ = std::fs::remove_file(input);
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn native_backend_loads_fasta_and_writes_upstream_cluster_rows() {
        let stem = format!("diamond-cluster-output-{}", std::process::id());
        let input = std::env::temp_dir().join(format!("{stem}.faa"));
        let output = std::env::temp_dir().join(format!("{stem}.tsv"));
        std::fs::write(&input, ">centroid description\nAAAA\n>member\nAAA\n").unwrap();
        let database = NativeBackend::load(input.to_str().unwrap()).unwrap();
        assert_eq!(database.block.records.len(), 2);
        assert_eq!(database.block.letter_count(), 7);

        let mut backend = NativeBackend::new(output.to_str().unwrap());
        backend
            .write_output(&database, WrapperOutput::Dense(vec![0, 0]), false, true)
            .unwrap();
        assert_eq!(
            std::fs::read_to_string(&output).unwrap(),
            "centroid\tmember\ncentroid\tcentroid\ncentroid\tmember\n"
        );
        let _ = std::fs::remove_file(input);
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn command_specific_default_identity_changes_compact_clustering() {
        let stem = format!("diamond-cluster-defaults-{}", std::process::id());
        let input = std::env::temp_dir().join(format!("{stem}.faa"));
        let original = "ARNDCQEGHILKMFPSTWYVARNDCQEGHILKMFPSTWYV";
        let changed = "VVVVVVVVHILKMFPSTWYVARNDCQEGHILKMFPSTWYV";
        std::fs::write(&input, format!(">a\n{original}\n>b\n{changed}\n")).unwrap();

        for (workflow, expected_clusters) in [
            (ClusterWorkflow::Cascaded, 1),
            (ClusterWorkflow::LinClust, 2),
            (ClusterWorkflow::DeepClust, 1),
        ] {
            let output = std::env::temp_dir().join(format!("{stem}-{workflow:?}.tsv"));
            let mut cfg = config(workflow);
            cfg.database = input.to_string_lossy().into_owned();
            cfg.output = output.to_string_lossy().into_owned();
            cfg.approx_id = None;
            run(&cfg).unwrap();
            let text = std::fs::read_to_string(&output).unwrap();
            let centroids = text
                .lines()
                .map(|line| line.split_once('\t').unwrap().0)
                .collect::<std::collections::BTreeSet<_>>();
            assert_eq!(centroids.len(), expected_clusters, "{workflow:?}: {text}");
            let _ = std::fs::remove_file(output);
        }
        let _ = std::fs::remove_file(input);
    }
}
