//! Native clustering command boundary.
//!
//! The structurally translated upstream workflow lives in
//! [`crate::cluster::cascaded`]. Database, search, and output integration use
//! the translated backend traits so each command retains its upstream workflow
//! rather than substituting a shared clustering algorithm.

use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::basic::value::{Loc, OId, SequenceType, SuperBlockId};
use crate::cluster::cascaded::cascaded::{self, CascadedBackend, CascadedRoundSummary};
use crate::cluster::cascaded::wrapper::{
    seq_mem_use, BestCentroid, CascadedWrapperBackend, CascadedWrapperConfig, CentroidSearchConfig,
    LengthSortedBlock, WrapperCommand, WrapperOutput,
};
use crate::cluster::cascaded::CascadedConfig;
use crate::cluster::cascaded::GraphAlgo;
use crate::commands::blastp::{self, BlastpConfig};
use crate::data::fasta::{self, FastaRecord};
use crate::data::sequence_file::DbFilter;
use crate::masking::{MaskingAlgo, MaskingMode};
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
    pub matrix: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub max_evalue: f64,
    pub comp_based_stats: i32,
    pub masking: MaskingMode,
    pub soft_masking: MaskingAlgo,
    pub motif_masking: String,
    pub extension_mode: String,
    pub member_cover: Option<f64>,
    pub mutual_cover: Option<f64>,
    /// `None` is important: upstream applies command-specific defaults in
    /// `Cluster::init_thresholds` (50, 90, and 0 respectively).
    pub approx_id: Option<f64>,
    pub cluster_steps: Vec<String>,
    pub connected_component_depth: Vec<String>,
    pub round_coverage: Vec<String>,
    pub round_approx_id: Vec<String>,
    pub kmer_ranking: bool,
    pub no_reassign: bool,
    pub weighted_gvc: bool,
    pub strict_gvc: bool,
    pub graph_algo: GraphAlgo,
    pub memory_limit: u64,
    pub oid_output: bool,
    pub simple_header: bool,
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
    translated.thresholds.member_cover = config.member_cover;
    translated.thresholds.mutual_cover = config.mutual_cover;
    translated.thresholds.approx_min_id = config.approx_id;
    translated.memory_limit = config.memory_limit;
    translated.oid_output = config.oid_output;
    translated.simple_header = config.simple_header;
    translated.core.threads = threads;
    translated.core.memory_limit = config.memory_limit;
    translated.core.matrix.clone_from(&config.matrix);
    translated.core.gap_open = config.gap_open;
    translated.core.gap_extend = config.gap_extend;
    translated.core.max_evalue = config.max_evalue;
    translated.core.comp_based_stats = config.comp_based_stats;
    translated.core.masking = config.masking;
    translated.core.soft_masking = config.soft_masking;
    translated
        .core
        .motif_masking
        .clone_from(&config.motif_masking);
    translated
        .core
        .extension_mode
        .clone_from(&config.extension_mode);
    translated.core.helpers.cluster_steps = config.cluster_steps.clone();
    translated.core.helpers.connected_component_depth = config.connected_component_depth.clone();
    translated.core.round_coverage = config.round_coverage.clone();
    translated.core.round_approx_id = config.round_approx_id.clone();
    translated.core.kmer_ranking = config.kmer_ranking;
    translated.core.no_gvc_reassign = config.no_reassign;
    translated.core.weighted_gvc = config.weighted_gvc;
    translated.core.strict_gvc = config.strict_gvc;
    translated.core.graph_algo = config.graph_algo;
    translated.core.alignment_output = config.alignment_output.clone();
    translated
}

#[derive(Clone)]
struct NativeBlock {
    records: Vec<FastaRecord>,
    original_oids: Vec<OId>,
    sequence_count: OId,
    letters: u64,
    source: Option<NativeBlockSource>,
    soft_tantan_masker: Arc<crate::masking::tantan::TantanMasker>,
}

#[derive(Clone)]
enum NativeBlockSource {
    Input(PathBuf),
    Spill(Arc<SpillFile>),
}

struct SpillFile {
    path: PathBuf,
}

impl Drop for SpillFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
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

/// Literal translation of `Search::keep_target_id` and the query ordering
/// predicate in upstream `run/double_indexed.cpp`.
fn keep_target_id(config: &cascaded::CascadedSearchConfig) -> bool {
    keep_target_id_values(
        config.min_length_ratio,
        config.global_ranking_targets,
        config.self_search,
        config.current_ref_block,
        config.lin_stage1_combo,
    )
}

fn keep_target_id_values(
    min_length_ratio: f64,
    global_ranking_targets: bool,
    self_search: bool,
    current_ref_block: i32,
    lin_stage1_combo: bool,
) -> bool {
    min_length_ratio != 0.0
        || global_ranking_targets
        || (self_search && current_ref_block == 0)
        || lin_stage1_combo
}

fn should_length_sort_queries(config: &cascaded::CascadedSearchConfig) -> bool {
    should_length_sort_query_values(
        keep_target_id(config),
        config.lin_stage1_query,
        config.kmer_ranking,
        config.min_length_ratio,
    )
}

fn should_length_sort_query_values(
    keep_target_id: bool,
    lin_stage1_query: bool,
    kmer_ranking: bool,
    min_length_ratio: f64,
) -> bool {
    (!keep_target_id && lin_stage1_query && !kmer_ranking) || min_length_ratio > 0.0
}

impl NativeBlock {
    fn letter_count(&self) -> u64 {
        self.letters
    }

    fn ensure_loaded(&mut self) -> Result<(), String> {
        if self.records.len() == self.sequence_count as usize {
            return Ok(());
        }
        let source = self
            .source
            .as_ref()
            .ok_or_else(|| "Sequence block has no backing store".to_owned())?;
        let (records, original_oids) = match source {
            NativeBlockSource::Input(path) => {
                let mut records = Vec::with_capacity(self.sequence_count as usize);
                visit_input_records(path, |oid, record| {
                    records.push((oid as OId, record));
                    Ok(())
                })?;
                let original_oids = records.iter().map(|(oid, _)| *oid).collect();
                let records = records.into_iter().map(|(_, record)| record).collect();
                (records, original_oids)
            }
            NativeBlockSource::Spill(file) => read_spill_records(&file.path)?,
        };
        if records.len() != self.sequence_count as usize {
            return Err("Sequence block backing store has the wrong record count".to_owned());
        }
        self.records = records;
        self.original_oids = original_oids;
        Ok(())
    }
}

struct NativeDatabase {
    path: PathBuf,
    block: NativeBlock,
    seqids: Vec<String>,
    lengths: Vec<Loc>,
}

struct NativeBackend {
    output: PathBuf,
    centroids: NativeBlock,
    /// Upstream's process-global `Masking` instance is derived from the score
    /// matrix once before `Cascaded::run`, then reused by every search round.
    soft_tantan_masker: std::sync::Arc<crate::masking::tantan::TantanMasker>,
}

fn resolve_cluster_input(path: &str) -> Result<PathBuf, String> {
    let requested = Path::new(path);
    if requested.exists() {
        return Ok(requested.to_path_buf());
    }
    let mut appended = requested.as_os_str().to_owned();
    appended.push(".dmnd");
    let resolved = PathBuf::from(appended);
    if resolved.exists() {
        Ok(resolved)
    } else {
        Err(format!("Database file not found: {path}"))
    }
}

fn visit_input_records<F>(path: &Path, mut visit: F) -> Result<(), String>
where
    F: FnMut(usize, FastaRecord) -> Result<(), String>,
{
    if path
        .extension()
        .is_some_and(|extension| extension == "dmnd")
    {
        crate::data::dmnd_reader::for_each_dmnd_record(path, |oid, id, sequence| {
            visit(
                oid,
                FastaRecord {
                    id: id.to_owned(),
                    sequence: sequence.iter().map(|&letter| letter as i8).collect(),
                },
            )
            .map_err(io::Error::other)
        })
        .map(|_| ())
        .map_err(|error| error.to_string())
    } else {
        let mut oid = 0usize;
        fasta::for_each_fasta_file(path, SequenceType::AminoAcid, |record| {
            let current = oid;
            oid += 1;
            visit(current, record).map_err(io::Error::other)
        })
        .map_err(|error| error.to_string())
    }
}

fn create_spill_file() -> Result<(Arc<SpillFile>, BufWriter<File>), String> {
    let data = crate::util::io::TempFileData::init_in(false, std::env::temp_dir(), true)
        .map_err(|error| error.to_string())?;
    let path = PathBuf::from(&data.name);
    drop(data);
    let file = OpenOptions::new()
        .write(true)
        .truncate(true)
        .open(&path)
        .map_err(|error| error.to_string())?;
    Ok((Arc::new(SpillFile { path }), BufWriter::new(file)))
}

fn write_spill_record(out: &mut BufWriter<File>, oid: OId, record: &FastaRecord) -> io::Result<()> {
    out.write_all(&oid.to_le_bytes())?;
    out.write_all(&(record.sequence.len() as u64).to_le_bytes())?;
    for &letter in &record.sequence {
        out.write_all(&[letter as u8])?;
    }
    Ok(())
}

fn read_spill_records(path: &Path) -> Result<(Vec<FastaRecord>, Vec<OId>), String> {
    let mut input = BufReader::new(File::open(path).map_err(|error| error.to_string())?);
    let mut records = Vec::new();
    let mut original_oids = Vec::new();
    loop {
        let mut oid_bytes = [0u8; 8];
        match input.read(&mut oid_bytes[..1]) {
            Ok(0) => break,
            Ok(1) => input
                .read_exact(&mut oid_bytes[1..])
                .map_err(|error| error.to_string())?,
            Ok(_) => unreachable!(),
            Err(error) => return Err(error.to_string()),
        }
        let mut length_bytes = [0u8; 8];
        input
            .read_exact(&mut length_bytes)
            .map_err(|error| error.to_string())?;
        let length = usize::try_from(u64::from_le_bytes(length_bytes))
            .map_err(|_| "Spilled sequence length exceeds address space".to_owned())?;
        let mut sequence = vec![0u8; length];
        input
            .read_exact(&mut sequence)
            .map_err(|error| error.to_string())?;
        original_oids.push(OId::from_le_bytes(oid_bytes));
        records.push(FastaRecord {
            // Upstream SequenceFile::length_sort writes the literal placeholder
            // title "X" into every temporary FastaFile.
            id: "X".to_owned(),
            sequence: sequence.into_iter().map(|letter| letter as i8).collect(),
        });
    }
    Ok((records, original_oids))
}

impl NativeBackend {
    fn new(output: &str, matrix: &str, gap_open: i32, gap_extend: i32) -> Result<Self, String> {
        let score_matrix =
            crate::stats::score_matrix::ScoreMatrix::new(matrix, gap_open, gap_extend, 0, 1, 0)?;
        let soft_tantan_masker = std::sync::Arc::new(
            crate::masking::tantan::TantanMasker::from_score_matrix(&score_matrix, 0.9),
        );
        Ok(Self {
            output: PathBuf::from(output),
            centroids: NativeBlock {
                records: Vec::new(),
                original_oids: Vec::new(),
                sequence_count: 0,
                letters: 0,
                source: None,
                soft_tantan_masker: soft_tantan_masker.clone(),
            },
            soft_tantan_masker,
        })
    }

    fn load(&self, path: &str) -> Result<NativeDatabase, String> {
        let resolved = resolve_cluster_input(path)?;
        let mut seqids = Vec::new();
        let mut lengths = Vec::new();
        let mut letters = 0u64;
        visit_input_records(&resolved, |oid, record| {
            if oid != lengths.len() {
                return Err("Input database OIDs are not contiguous".to_owned());
            }
            let length = Loc::try_from(record.sequence.len())
                .map_err(|_| "Sequence length exceeds Loc range".to_owned())?;
            letters = letters
                .checked_add(length as u64)
                .ok_or_else(|| "Database letter count overflow".to_owned())?;
            lengths.push(length);
            // `SequenceFile::add_seqid_mapping` stores `Util::Seq::seqid`,
            // not the complete FASTA title.
            seqids.push(crate::util::sequence::seqid(&record.id));
            Ok(())
        })?;
        let sequence_count = lengths.len() as OId;
        Ok(NativeDatabase {
            path: resolved.clone(),
            block: NativeBlock {
                records: Vec::new(),
                original_oids: Vec::new(),
                sequence_count,
                letters,
                source: Some(NativeBlockSource::Input(resolved)),
                soft_tantan_masker: self.soft_tantan_masker.clone(),
            },
            seqids,
            lengths,
        })
    }
}

impl CascadedBackend for NativeBlock {
    fn sequence_count(&self) -> u64 {
        self.sequence_count
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
        self.ensure_loaded()?;
        let mut active = self
            .records
            .iter()
            .enumerate()
            .filter(|(oid, _)| filter.is_none_or(|set| set.get(*oid as OId)))
            .map(|(oid, record)| (oid as OId, record.clone()))
            .collect::<Vec<_>>();
        length_sort_linear_active(&mut active, should_length_sort_queries(config));
        let local_to_block = active.iter().map(|(oid, _)| *oid).collect::<Vec<_>>();
        let records = active
            .into_iter()
            .map(|(_, record)| record)
            .collect::<Vec<_>>();
        let target_count = records.len() as i64;
        let extension_mode = if config.extension_mode.is_empty() {
            None
        } else {
            Some(crate::search::sensitivity::ExtensionMode::from_string(
                &config.extension_mode,
            )?)
        };
        let search = BlastpConfig {
            query_files: Vec::new(),
            database: String::new(),
            output: None,
            matrix: config.matrix.clone(),
            gap_open: config.gap_open,
            gap_extend: config.gap_extend,
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
            masking: config.masking,
            motif_masking: config.motif_masking.clone(),
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
            target_indexed: false,
            translated_query_layout: None,
        };
        let queries = (!config.self_search).then(|| records.clone());
        let mut edges = blastp::run_edges_in_memory(
            &search,
            records,
            queries,
            config.approx_min_id,
            config.lin_stage1_query,
            config.lin_stage1_target,
            config.self_search,
            config.query_or_target_cover,
            config.mutual_cover,
            crate::chaining::HammingExtConfig {
                hamming_ext: config.hamming_ext,
                approx_min_id: config.approx_min_id,
                query_or_target_cover: config.query_or_target_cover,
                query_cover: config.query_cover,
                subject_cover: config.subject_cover,
                max_evalue: config.max_evalue,
                diag_filter_cov: config.diag_filter_cov,
                diag_filter_id: config.diag_filter_id,
            },
            Some(self.soft_tantan_masker.as_ref()),
            config.soft_masking,
            config.freq_masking,
            (config.db_size != 0).then_some(config.db_size),
            extension_mode,
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
        self.load(path)
    }

    fn database_name(&self, database: &Self::Database) -> String {
        database.path.to_string_lossy().into_owned()
    }

    fn database_is_blast(&self, _database: &Self::Database) -> bool {
        false
    }

    fn database_sequence_count(&self, database: &Self::Database) -> OId {
        database.block.sequence_count
    }

    fn database_letters(&self, database: &Self::Database) -> u64 {
        database.block.letter_count()
    }

    fn length_sort(
        &mut self,
        database: &mut Self::Database,
        memory_limit: i64,
        minimizer_window: i32,
        sketch_size: i32,
    ) -> Result<Vec<LengthSortedBlock<Self::Block>>, String> {
        let mut lengths = database
            .lengths
            .iter()
            .copied()
            .enumerate()
            .map(|(oid, length)| (length, oid as OId))
            .collect::<Vec<_>>();
        // Exact `greater<pair<Loc,OId>>`: descending length and descending OID
        // for ties. The subsequent OID sort restores input order within each
        // temporary super-block.
        lengths.sort_unstable_by(|left, right| right.cmp(left));

        let mut block_by_oid = vec![0usize; lengths.len()];
        let mut block_counts = Vec::<OId>::new();
        let mut block_letters = Vec::<u64>::new();
        let mut size = 0i64;
        let mut sequences = 0u64;
        let mut letters_in_block = 0u64;
        let mut block = 0usize;
        for (length, oid) in &lengths {
            size = size.saturating_add(seq_mem_use(*length, 0, 1, minimizer_window, sketch_size));
            letters_in_block += *length as u64;
            sequences += 1;
            block_by_oid[*oid as usize] = block;
            if (size >= memory_limit && letters_in_block >= 1)
                || sequences >= SuperBlockId::MAX as u64
            {
                block_counts.push(sequences);
                block_letters.push(letters_in_block);
                block += 1;
                size = 0;
                sequences = 0;
                letters_in_block = 0;
            }
        }
        if size > 0 {
            block_counts.push(sequences);
            block_letters.push(letters_in_block);
        }

        let mut files = Vec::with_capacity(block_counts.len());
        let mut writers = Vec::with_capacity(block_counts.len());
        let mut mappings = (0..block_counts.len())
            .map(|index| Vec::with_capacity(block_counts[index] as usize))
            .collect::<Vec<_>>();
        for _ in 0..block_counts.len() {
            let (file, writer) = create_spill_file()?;
            files.push(file);
            writers.push(writer);
        }
        visit_input_records(&database.path, |oid, record| {
            let block = *block_by_oid
                .get(oid)
                .ok_or_else(|| format!("Input OID out of range during length sort: {oid}"))?;
            write_spill_record(&mut writers[block], oid as OId, &record)
                .map_err(|error| error.to_string())?;
            mappings[block].push(oid as OId);
            Ok(())
        })?;
        for writer in &mut writers {
            writer.flush().map_err(|error| error.to_string())?;
        }
        drop(writers);

        Ok(files
            .into_iter()
            .zip(mappings)
            .enumerate()
            .map(|(index, (file, super_block_id_to_oid))| LengthSortedBlock {
                block: NativeBlock {
                    records: Vec::new(),
                    original_oids: Vec::new(),
                    sequence_count: block_counts[index],
                    letters: block_letters[index],
                    source: Some(NativeBlockSource::Spill(file)),
                    soft_tantan_masker: self.soft_tantan_masker.clone(),
                },
                super_block_id_to_oid,
            })
            .collect())
    }

    fn block_sequence_count(&self, block: &Self::Block) -> OId {
        block.sequence_count
    }

    fn block_letters(&self, block: &Self::Block) -> u64 {
        block.letter_count()
    }

    fn centroid_sequence_count(&self) -> OId {
        self.centroids.sequence_count
    }

    fn centroid_letters(&self) -> u64 {
        self.centroids.letter_count()
    }

    fn run_centroid_search(
        &mut self,
        block: &mut Self::Block,
        config: &CentroidSearchConfig,
        consumer: &mut BestCentroid,
    ) -> Result<(), String> {
        block.ensure_loaded()?;
        self.centroids.ensure_loaded()?;
        let extension_mode = if config.extension_mode.is_empty() {
            None
        } else {
            Some(crate::search::sensitivity::ExtensionMode::from_string(
                &config.extension_mode,
            )?)
        };
        let search = BlastpConfig {
            query_files: Vec::new(),
            database: String::new(),
            output: None,
            matrix: config.matrix.clone(),
            gap_open: config.gap_open,
            gap_extend: config.gap_extend,
            max_evalue: config.max_evalue,
            max_target_seqs: config.max_target_seqs,
            ext_chunk_size: 0,
            toppercent: config.top_percent,
            global_ranking_targets: 0,
            min_id: 0.0,
            threads: config.threads,
            outfmt: vec!["6".to_owned()],
            sensitivity: config.sensitivity,
            masking: config.masking,
            motif_masking: config.motif_masking.clone(),
            min_query_len: 0,
            query_cover: config.query_cover,
            subject_cover: config.subject_cover,
            comp_based_stats: CbsMode::parse(&config.comp_based_stats.to_string()),
            no_self_hits: false,
            ungapped_xdrop_bits: 12.3,
            memory_limit: None,
            tmpdir: PathBuf::new(),
            target_indexed: false,
            translated_query_layout: None,
        };
        let mut targets = self
            .centroids
            .records
            .iter()
            .cloned()
            .enumerate()
            .collect::<Vec<_>>();
        if config.lin_stage1_target {
            targets.sort_unstable_by(|(left_oid, left), (right_oid, right)| {
                right
                    .sequence
                    .len()
                    .cmp(&left.sequence.len())
                    .then_with(|| right_oid.cmp(left_oid))
            });
        }
        let target_to_centroid = targets
            .iter()
            .map(|(centroid, _)| *centroid as OId)
            .collect::<Vec<_>>();
        let target_records = targets
            .into_iter()
            .map(|(_, record)| record)
            .collect::<Vec<_>>();
        let edges = blastp::run_edges_in_memory(
            &search,
            target_records,
            Some(block.records.clone()),
            config.approx_min_id,
            config.lin_stage1_query,
            config.lin_stage1_target,
            config.self_search,
            config.query_or_target_cover,
            config.query_cover > 0.0 && config.subject_cover > 0.0,
            crate::chaining::HammingExtConfig {
                hamming_ext: config.hamming_ext,
                approx_min_id: config.approx_min_id,
                query_or_target_cover: config.query_or_target_cover,
                query_cover: config.query_cover,
                subject_cover: config.subject_cover,
                max_evalue: config.max_evalue,
                diag_filter_cov: config.diag_filter_cov,
                diag_filter_id: config.diag_filter_id,
            },
            Some(self.soft_tantan_masker.as_ref()),
            config.soft_masking,
            config.freq_masking,
            (config.db_size != 0).then_some(config.db_size),
            extension_mode,
        )
        .map_err(|error| error.to_string())?;
        let mut bytes = Vec::with_capacity(edges.len() * crate::output::edge::EdgeData::SIZE);
        for mut edge in edges {
            edge.edge.target = *target_to_centroid
                .get(edge.edge.target as usize)
                .ok_or_else(|| "Centroid search target OID out of range".to_owned())?;
            edge.edge
                .write(&mut bytes)
                .map_err(|error| error.to_string())?;
        }
        consumer.consume(&bytes)
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
            sequence_count: ids.len() as OId,
            letters: ids
                .iter()
                .map(|&id| block.records[id as usize].sequence.len() as u64)
                .sum(),
            source: None,
            soft_tantan_masker: self.soft_tantan_masker.clone(),
        })
    }

    fn cascaded_database(
        &mut self,
        database: &mut Self::Database,
        linear: bool,
        config: &mut CascadedConfig,
    ) -> Result<Vec<SuperBlockId>, String> {
        database.block.ensure_loaded()?;
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
        self.centroids.sequence_count = self.centroids.records.len() as OId;
        self.centroids.letters = self
            .centroids
            .records
            .iter()
            .map(|record| record.sequence.len() as u64)
            .sum();
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
                let centroid_title = database
                    .seqids
                    .get(centroid as usize)
                    .ok_or_else(|| format!("Centroid OID out of range: {centroid}"))?;
                let member_title = database
                    .seqids
                    .get(member as usize)
                    .ok_or_else(|| format!("Member OID out of range: {member}"))?;
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
    let mut backend = NativeBackend::new(
        &config.output,
        &config.matrix,
        config.gap_open,
        config.gap_extend,
    )
    .map_err(io::Error::other)?;
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
            matrix: "blosum62".to_owned(),
            gap_open: -1,
            gap_extend: -1,
            max_evalue: 0.001,
            comp_based_stats: 1,
            masking: MaskingMode::None,
            soft_masking: MaskingAlgo::Tantan,
            motif_masking: String::new(),
            extension_mode: String::new(),
            member_cover: None,
            mutual_cover: None,
            approx_id: None,
            cluster_steps: Vec::new(),
            connected_component_depth: Vec::new(),
            round_coverage: Vec::new(),
            round_approx_id: Vec::new(),
            kmer_ranking: false,
            no_reassign: false,
            weighted_gvc: false,
            strict_gvc: false,
            graph_algo: GraphAlgo::GreedyVertexCover,
            memory_limit: 16_000_000_000,
            oid_output: false,
            simple_header: false,
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
    fn query_sort_predicate_preserves_self_search_target_ids() {
        assert!(!keep_target_id_values(0.0, false, true, 1, false));
        assert!(keep_target_id_values(0.0, false, true, 0, false));
        assert!(keep_target_id_values(0.1, false, false, 1, false));
        assert!(keep_target_id_values(0.0, true, false, 1, false));
        assert!(keep_target_id_values(0.0, false, false, 1, true));
        assert!(should_length_sort_query_values(false, true, false, 0.0));
        assert!(!should_length_sort_query_values(true, true, false, 0.0));
        assert!(!should_length_sort_query_values(false, true, true, 0.0));
        assert!(!should_length_sort_query_values(false, false, false, 0.0));
        assert!(should_length_sort_query_values(true, false, true, 0.5));
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
    fn translated_config_forwards_cluster_controls_without_defaulting_member_cover_early() {
        let mut cfg = config(ClusterWorkflow::Cascaded);
        cfg.mutual_cover = Some(77.0);
        cfg.cluster_steps = vec!["faster_lin".into(), "default".into()];
        cfg.connected_component_depth = vec!["2".into()];
        cfg.round_coverage = vec!["81".into(), "82".into()];
        cfg.round_approx_id = vec!["55".into(), "60".into()];
        cfg.kmer_ranking = true;
        cfg.no_reassign = true;
        cfg.weighted_gvc = true;
        cfg.strict_gvc = true;
        cfg.graph_algo = GraphAlgo::LengthSorted;
        cfg.memory_limit = 3_000_000_000;
        cfg.matrix = "pam250".to_owned();
        cfg.gap_open = 15;
        cfg.gap_extend = 3;
        cfg.oid_output = true;
        cfg.simple_header = true;

        let translated = translated_config(&cfg);
        assert_eq!(translated.thresholds.member_cover, None);
        assert_eq!(translated.thresholds.mutual_cover, Some(77.0));
        assert_eq!(translated.memory_limit, 3_000_000_000);
        assert_eq!(translated.core.memory_limit, 3_000_000_000);
        assert_eq!(translated.core.matrix, "pam250");
        assert_eq!(translated.core.gap_open, 15);
        assert_eq!(translated.core.gap_extend, 3);
        assert_eq!(translated.core.helpers.cluster_steps, cfg.cluster_steps);
        assert_eq!(
            translated.core.helpers.connected_component_depth,
            cfg.connected_component_depth
        );
        assert_eq!(translated.core.round_coverage, cfg.round_coverage);
        assert_eq!(translated.core.round_approx_id, cfg.round_approx_id);
        assert!(translated.core.kmer_ranking);
        assert!(translated.core.no_gvc_reassign);
        assert!(translated.core.weighted_gvc);
        assert!(translated.core.strict_gvc);
        assert_eq!(translated.core.graph_algo, GraphAlgo::LengthSorted);
        assert!(translated.oid_output);
        assert!(translated.simple_header);
    }

    #[test]
    fn native_backend_loads_fasta_and_writes_upstream_cluster_rows() {
        let stem = format!("diamond-cluster-output-{}", std::process::id());
        let input = std::env::temp_dir().join(format!("{stem}.faa"));
        let output = std::env::temp_dir().join(format!("{stem}.tsv"));
        std::fs::write(&input, ">centroid description\nAAAA\n>member\nAAA\n").unwrap();
        let mut backend = NativeBackend::new(output.to_str().unwrap(), "blosum62", -1, -1).unwrap();
        let mut database = backend.load(input.to_str().unwrap()).unwrap();
        assert_eq!(database.block.sequence_count, 2);
        assert!(database.block.records.is_empty());
        database.block.ensure_loaded().unwrap();
        assert_eq!(database.block.records.len(), 2);
        assert_eq!(database.block.letter_count(), 7);

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
    fn native_length_sort_uses_two_pass_spills_and_cpp_block_order() {
        let stem = format!(
            "diamond-cluster-length-sort-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        );
        let input = std::env::temp_dir().join(format!("{stem}.faa"));
        let output = std::env::temp_dir().join(format!("{stem}.tsv"));
        std::fs::write(
            &input,
            ">zero\nAAAA\n>one\nAAAAAAAA\n>two\nAAAAAA\n>three\nAAAAAAAA\n>four\nAAA\n",
        )
        .unwrap();
        let mut backend = NativeBackend::new(output.to_str().unwrap(), "blosum62", -1, -1).unwrap();
        let mut database = backend.load(input.to_str().unwrap()).unwrap();
        assert!(database.block.records.is_empty());

        // Faster uses sketch_size=21. A 360-byte budget forces multiple
        // super-blocks while retaining more than one sequence in a block.
        let mut blocks = backend.length_sort(&mut database, 360, 0, 21).unwrap();
        assert_eq!(blocks.len(), 3);
        assert_eq!(blocks[0].super_block_id_to_oid, [1, 3]);
        assert_eq!(blocks[1].super_block_id_to_oid, [0, 2]);
        assert_eq!(blocks[2].super_block_id_to_oid, [4]);
        for item in &blocks {
            assert!(item.block.records.is_empty());
        }
        for item in &mut blocks {
            item.block.ensure_loaded().unwrap();
            assert_eq!(item.block.records.len(), item.super_block_id_to_oid.len());
            assert!(item.block.records.iter().all(|record| record.id == "X"));
            assert_eq!(item.block.original_oids, item.super_block_id_to_oid);
        }

        let spill_paths = blocks
            .iter()
            .map(|item| match item.block.source.as_ref().unwrap() {
                NativeBlockSource::Spill(file) => file.path.clone(),
                NativeBlockSource::Input(_) => panic!("expected spill-backed super-block"),
            })
            .collect::<Vec<_>>();
        assert!(spill_paths.iter().all(|path| path.exists()));
        drop(blocks);
        assert!(spill_paths.iter().all(|path| !path.exists()));
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
