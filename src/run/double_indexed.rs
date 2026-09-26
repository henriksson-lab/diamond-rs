//! Double-indexed search orchestration from `diamond/src/run/double_indexed.cpp`.
//!
//! The original translation unit coordinates process globals, sequence-file
//! implementations, temporary consumers, multiprocessing stacks, masking,
//! seed search, alignment, and output.  This port owns the decision making and
//! iteration locally while exposing those concrete services through
//! [`DoubleIndexedBackend`].

use std::cmp::max;
use std::path::{Path, PathBuf};

use crate::config::Sensitivity;
use crate::data::flags::SeedEncoding;
use crate::data::sequence_file::{SequenceFileFlags, SequenceFileType};
use crate::masking::MaskingAlgo;

use super::config::{Algo, Round};

pub const MAX_INDEX_QUERY_SIZE: u64 = 32 * 1024 * 1024;
pub const MAX_HASH_SET_SIZE: usize = 8 * 1024 * 1024;
pub const MIN_QUERY_INDEXED_DB_SIZE: u64 = 256 * 1024 * 1024;

pub const STACK_ALIGN_TODO: &str = "align_todo";
pub const STACK_ALIGN_WIP: &str = "align_wip";
pub const STACK_ALIGN_DONE: &str = "align_done";
pub const STACK_JOIN_TODO: &str = "join_todo";
pub const STACK_JOIN_WIP: &str = "join_wip";
pub const STACK_JOIN_REDO: &str = "join_redo";
pub const STACK_JOIN_DONE: &str = "join_done";

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct OutputFlags(pub u32);

impl OutputFlags {
    pub const SSEQID: Self = Self(1 << 0);
    pub const FULL_TITLES: Self = Self(1 << 1);
    pub const ALL_SEQIDS: Self = Self(1 << 2);
    pub const TARGET_SEQS: Self = Self(1 << 3);
    pub const SELF_ALN_SCORES: Self = Self(1 << 4);

    pub const fn contains(self, flag: Self) -> bool {
        self.0 & flag.0 != 0
    }
}

impl std::ops::BitOr for OutputFlags {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceBlock {
    pub sequence_count: usize,
    pub source_sequence_count: usize,
    pub letters: u64,
    pub raw_bytes: u64,
    pub oid_begin: u64,
    pub oid_end: u64,
    pub has_ids: bool,
}

impl SequenceBlock {
    pub const fn empty() -> Self {
        Self {
            sequence_count: 0,
            source_sequence_count: 0,
            letters: 0,
            raw_bytes: 0,
            oid_begin: 0,
            oid_end: 0,
            has_ids: false,
        }
    }

    pub const fn is_empty(&self) -> bool {
        self.sequence_count == 0
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatabaseInfo {
    pub file_type: SequenceFileType,
    pub sequence_count: u64,
    pub letters: u64,
    pub total_blocks: usize,
    pub file_name: String,
    pub alias_taxidlist: Option<String>,
    pub alias_seqidlist: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuerySeedMode {
    None,
    Hashed,
    Contiguous,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistogramFilter {
    None,
    HashedQuerySeeds,
    ContiguousQuerySeeds,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputTarget {
    Master,
    Temporary,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BufferPlan {
    pub keep_target_id: bool,
    pub target_histogram_entries: usize,
    pub query_histogram_entries: Option<usize>,
    pub index_chunks: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub enum WorkflowOp {
    LogRss,
    SortReferenceByLength,
    SortQueryByLength,
    PreserveUnmaskedReference,
    MaskReference(MaskingAlgo),
    MaskQuery(MaskingAlgo),
    ComputeReferenceSelfAlignment,
    ComputeQuerySelfAlignment,
    InitDictionary {
        query_block: usize,
        reference_block: usize,
    },
    InitDictionaryBlock {
        reference_block: usize,
        sequence_count: usize,
        persist: bool,
    },
    CloseDictionaryBlock {
        persist: bool,
    },
    InitHitBuffer,
    InitGlobalRankingBuffer,
    AllocateGlobalRankingTable {
        entries: usize,
    },
    BuildReferenceHistogram(HistogramFilter),
    BuildQueryHistogram {
        hashed: bool,
    },
    AllocateBuffers(BufferPlan),
    FreeBuffers,
    LoadTargetSeedIndex(String),
    SearchShape {
        shape: usize,
        query_iteration: usize,
    },
    UpdateGlobalRanking,
    FinishHitBuffer,
    ClearQueryMasking,
    OpenTemporaryOutput {
        path: Option<PathBuf>,
    },
    AlignQueries(OutputTarget),
    FinishIntermediate,
    CloseTemporaryOutput,
    BuildHashedQuerySeedSet,
    BuildContiguousQuerySeedSet,
    SetQuerySkipFromAligned,
    BuildDnaReferenceIndex,
    ClearQuerySeedSets,
    SetDatabaseFlags(SequenceFileFlags),
    MultiprocessingSearchBegin {
        query_block: usize,
        reference_block: usize,
    },
    MultiprocessingSearchEnd {
        query_block: usize,
        reference_block: usize,
    },
    WriteSkippedSelfBlock {
        query_block: usize,
        reference_block: usize,
    },
    ExtendGlobalRanking(OutputTarget),
    SetupSearch(Sensitivity),
    JoinBlocks {
        reference_blocks: usize,
        multiprocessing: bool,
    },
    WriteUnalignedQueries,
    WriteAlignedQueries,
    RecoverMultiprocessing {
        max_query_chunks: usize,
    },
    InitMultiprocessing {
        directory: PathBuf,
        chunk_letters: u64,
    },
    SavePartition {
        query_block: usize,
        path: PathBuf,
    },
    OpenOutput,
    InitDaa,
    PrintHeader,
    PrintFooter,
    FinishDaa {
        query_blocks: usize,
    },
    FinalizeOutput,
    CloseAuxiliaryOutputs,
    CloseQueryInput,
    WriteUnalignedTargets(String),
    CloseDatabase,
    FreeConfig,
    PrintDatabaseInfo,
    ApplyTaxonomyFilter {
        values: String,
        exclude: bool,
        separator: char,
    },
    ApplyAccessionFilter(String),
    SetScoreMatrixDatabaseLetters(u64),
    ConfigureNucleotideScoring,
    StatisticsReset,
    StatisticsPrint,
    Message(String),
}

pub trait DoubleIndexedBackend {
    fn l3_cache_size(&self) -> usize;
    fn database_info(&mut self, flags: SequenceFileFlags) -> Result<DatabaseInfo, String>;
    fn filtered_database_totals(&self) -> Option<(u64, u64)>;
    fn load_query_block(
        &mut self,
        max_letters: u64,
        self_search: bool,
        oid_offset: u64,
        flags: SequenceFileFlags,
    ) -> Result<Option<SequenceBlock>, String>;
    fn rewind_query(&mut self) -> Result<(), String>;
    fn load_reference_block(
        &mut self,
        max_letters: u64,
        chunk: Option<usize>,
        oid_offset: u64,
    ) -> Result<Option<SequenceBlock>, String>;
    fn reference_partition_chunks(&self) -> usize;
    fn query_seed_table_size(&mut self, mode: QuerySeedMode) -> Result<usize, String>;
    fn histogram_entries(&self, query: bool) -> usize;
    fn seed_partition_bits(&self, shape_weight: u32, threads: usize, chunks: u32) -> u32;
    fn multiprocessing_chunks(&mut self, query_block: usize) -> Result<Vec<usize>, String>;
    fn take_iteration_aligned(&mut self) -> u64;
    fn stop_exists(&self) -> bool;
    fn remove_stop(&mut self) -> Result<(), String>;
    fn perform(&mut self, operation: WorkflowOp) -> Result<(), String>;
}

#[derive(Debug, Clone, PartialEq)]
pub struct DoubleIndexedConfig {
    pub parallel_tmpdir: PathBuf,
    pub output_file: String,
    pub unaligned_targets: String,
    pub taxonlist: String,
    pub taxon_exclude: String,
    pub seqidlist: String,
    pub query_file_missing: bool,
    pub self_search: bool,
    pub target_indexed: bool,
    pub multiprocessing: bool,
    pub mp_recover: bool,
    pub mp_init: bool,
    pub mp_self: bool,
    pub mp_query_chunk: Option<usize>,
    pub global_ranking_targets: usize,
    pub hit_keep_target_id: bool,
    pub kmer_ranking: bool,
    pub target_exclusively_owned: bool,
    pub swipe_all: bool,
    pub blastn: bool,
    pub no_self_hits: bool,
    pub store_query_quality: bool,
    pub lin_stage1_query: bool,
    pub lin_stage1_target: bool,
    pub lin_stage1_combo: bool,
    pub frame_shift: i32,
    pub threads: usize,
    pub shape_count: usize,
    pub shape_weight: u32,
    pub chunk_size_billions: Option<f64>,
    pub block_size_letters: u64,
    pub algo: Algo,
    pub sensitivity: Vec<Round>,
    pub supports_query_indexed: Vec<Sensitivity>,
    pub index_chunks: Option<u32>,
    pub seed_encoding: SeedEncoding,
    pub query_masking: MaskingAlgo,
    pub target_masking: MaskingAlgo,
    pub soft_masking: MaskingAlgo,
    pub min_length_ratio: f64,
    pub minimizer_window: i32,
    pub sketch_size: i32,
    pub gapped_filter_evalue: f64,
    pub gapped_filter_evalue1: f64,
    pub output_flags: OutputFlags,
    pub output_titles_lazy: bool,
    pub output_is_daa: bool,
    pub output_needs_taxon_id_lists: bool,
    pub output_needs_taxon_nodes: bool,
    pub output_needs_taxon_scientific_names: bool,
    pub output_needs_taxon_ranks: bool,
    pub taxon_culling: bool,
    pub track_aligned_queries: bool,
    pub query_contexts: usize,
    pub db_size: Option<u64>,
}

impl Default for DoubleIndexedConfig {
    fn default() -> Self {
        Self {
            parallel_tmpdir: PathBuf::from("."),
            output_file: String::new(),
            unaligned_targets: String::new(),
            taxonlist: String::new(),
            taxon_exclude: String::new(),
            seqidlist: String::new(),
            query_file_missing: false,
            self_search: false,
            target_indexed: false,
            multiprocessing: false,
            mp_recover: false,
            mp_init: false,
            mp_self: false,
            mp_query_chunk: None,
            global_ranking_targets: 0,
            hit_keep_target_id: false,
            kmer_ranking: false,
            target_exclusively_owned: true,
            swipe_all: false,
            blastn: false,
            no_self_hits: false,
            store_query_quality: false,
            lin_stage1_query: false,
            lin_stage1_target: false,
            lin_stage1_combo: false,
            frame_shift: 0,
            threads: 1,
            shape_count: 1,
            shape_weight: 0,
            chunk_size_billions: None,
            block_size_letters: 2_000_000_000,
            algo: Algo::Auto,
            sensitivity: vec![Round::new(Sensitivity::Default, false)],
            supports_query_indexed: vec![
                Sensitivity::Faster,
                Sensitivity::Fast,
                Sensitivity::Default,
            ],
            index_chunks: None,
            seed_encoding: SeedEncoding::SpacedFactor,
            query_masking: MaskingAlgo::None,
            target_masking: MaskingAlgo::None,
            soft_masking: MaskingAlgo::None,
            min_length_ratio: 0.0,
            minimizer_window: 0,
            sketch_size: 0,
            gapped_filter_evalue: 0.0,
            gapped_filter_evalue1: 0.0,
            output_flags: OutputFlags::default(),
            output_titles_lazy: false,
            output_is_daa: false,
            output_needs_taxon_id_lists: false,
            output_needs_taxon_nodes: false,
            output_needs_taxon_scientific_names: false,
            output_needs_taxon_ranks: false,
            taxon_culling: false,
            track_aligned_queries: false,
            query_contexts: 1,
            db_size: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DoubleIndexedState {
    pub config: DoubleIndexedConfig,
    pub database: Option<DatabaseInfo>,
    pub query: Option<SequenceBlock>,
    pub target: Option<SequenceBlock>,
    pub current_query_block: usize,
    pub current_reference_block: usize,
    pub blocked_processing: bool,
    pub lazy_masking: bool,
    pub lin_stage1_target: bool,
    pub query_seed_mode: QuerySeedMode,
    pub seed_partition_bits: u32,
    pub query_aligned: Vec<bool>,
    pub aligned_total: u64,
    pub cutoff_gapped1: Option<f64>,
    pub cutoff_gapped2: Option<f64>,
}

impl DoubleIndexedState {
    pub fn new(config: DoubleIndexedConfig) -> Self {
        Self {
            lin_stage1_target: config.lin_stage1_target,
            config,
            database: None,
            query: None,
            target: None,
            current_query_block: 0,
            current_reference_block: 0,
            blocked_processing: false,
            lazy_masking: false,
            query_seed_mode: QuerySeedMode::None,
            seed_partition_bits: 0,
            query_aligned: Vec::new(),
            aligned_total: 0,
            cutoff_gapped1: None,
            cutoff_gapped2: None,
        }
    }
}

pub fn use_query_index(table_size: usize, l3_cache_size: usize) -> bool {
    table_size <= max(MAX_HASH_SET_SIZE, l3_cache_size)
}

pub fn get_ref_part_file_name(
    parallel_tmpdir: &Path,
    prefix: &str,
    query: usize,
    suffix: &str,
) -> PathBuf {
    let suffix = if suffix.is_empty() {
        String::new()
    } else {
        format!("{suffix}_")
    };
    parallel_tmpdir.join(format!("{prefix}_{suffix}{query}"))
}

pub fn get_ref_block_tmpfile_name(parallel_tmpdir: &Path, query: usize, block: usize) -> PathBuf {
    parallel_tmpdir.join(format!("ref_block_{query}_{block}"))
}

pub fn alloc_buffers<B: DoubleIndexedBackend>(
    state: &DoubleIndexedState,
    backend: &B,
) -> BufferPlan {
    BufferPlan {
        keep_target_id: keep_target_id(state),
        target_histogram_entries: backend.histogram_entries(false),
        query_histogram_entries: (!state.config.target_indexed)
            .then(|| backend.histogram_entries(true)),
        index_chunks: state.config.index_chunks.unwrap_or(1),
    }
}

pub fn keep_target_id(state: &DoubleIndexedState) -> bool {
    state.config.hit_keep_target_id
        || state.config.min_length_ratio != 0.0
        || state.config.global_ranking_targets > 0
        || (state.config.self_search && state.current_reference_block == 0)
        || state.config.lin_stage1_combo
}

pub fn run_ref_chunk<B: DoubleIndexedBackend>(
    query_iteration: usize,
    state: &mut DoubleIndexedState,
    backend: &mut B,
) -> Result<(), String> {
    backend.perform(WorkflowOp::LogRss)?;
    let query = state.query.as_ref().ok_or("Query block is not loaded.")?;
    let target = state
        .target
        .as_ref()
        .ok_or("Reference block is not loaded.")?;
    if (state.lin_stage1_target || state.config.min_length_ratio > 0.0)
        && !state.config.kmer_ranking
        && state.config.target_exclusively_owned
    {
        backend.perform(WorkflowOp::SortReferenceByLength)?;
    }
    if state.config.output_flags.contains(OutputFlags::TARGET_SEQS) {
        backend.perform(WorkflowOp::PreserveUnmaskedReference)?;
    }
    if state.config.target_masking != MaskingAlgo::None && !state.lazy_masking {
        backend.perform(WorkflowOp::MaskReference(state.config.target_masking))?;
    }
    if state
        .config
        .output_flags
        .contains(OutputFlags::SELF_ALN_SCORES)
    {
        backend.perform(WorkflowOp::ComputeReferenceSelfAlignment)?;
    }
    let persist_dictionary = state.config.output_is_daa || state.config.sensitivity.len() > 1;
    if ((state.blocked_processing || state.config.output_is_daa)
        && state.config.global_ranking_targets == 0)
        || state.config.sensitivity.len() > 1
    {
        if state.config.multiprocessing
            || (state.current_reference_block == 0
                && (!state.config.output_is_daa || state.current_query_block == 0)
                && query_iteration == 0)
        {
            backend.perform(WorkflowOp::InitDictionary {
                query_block: state.current_query_block,
                reference_block: state.current_reference_block,
            })?;
        }
        if state.config.global_ranking_targets == 0 {
            backend.perform(WorkflowOp::InitDictionaryBlock {
                reference_block: state.current_reference_block,
                sequence_count: target.sequence_count,
                persist: persist_dictionary,
            })?;
        }
    }
    backend.perform(if state.config.global_ranking_targets > 0 {
        WorkflowOp::InitGlobalRankingBuffer
    } else {
        WorkflowOp::InitHitBuffer
    })?;

    if !state.config.swipe_all {
        let filter = match state.query_seed_mode {
            QuerySeedMode::Contiguous => HistogramFilter::ContiguousQuerySeeds,
            QuerySeedMode::Hashed => HistogramFilter::HashedQuerySeeds,
            QuerySeedMode::None => HistogramFilter::None,
        };
        backend.perform(WorkflowOp::BuildReferenceHistogram(filter))?;
        backend.perform(WorkflowOp::AllocateBuffers(alloc_buffers(state, backend)))?;
        if state.config.target_indexed {
            let file = state
                .database
                .as_ref()
                .map(|db| format!("{}.seed_idx", db.file_name))
                .unwrap_or_else(|| ".seed_idx".to_owned());
            backend.perform(WorkflowOp::LoadTargetSeedIndex(file))?;
        }
        if !state.config.blastn {
            for shape in 0..state.config.shape_count {
                if state.config.global_ranking_targets > 0 {
                    backend.perform(WorkflowOp::InitGlobalRankingBuffer)?;
                }
                backend.perform(WorkflowOp::SearchShape {
                    shape,
                    query_iteration,
                })?;
                if state.config.global_ranking_targets > 0 {
                    backend.perform(WorkflowOp::UpdateGlobalRanking)?;
                }
            }
            if state.config.global_ranking_targets == 0 {
                backend.perform(WorkflowOp::FinishHitBuffer)?;
            }
        } else {
            backend.perform(WorkflowOp::BuildDnaReferenceIndex)?;
        }
        backend.perform(WorkflowOp::FreeBuffers)?;
        backend.perform(WorkflowOp::ClearQueryMasking)?;
        backend.perform(WorkflowOp::LogRss)?;
    }

    let temporary = (state.blocked_processing || state.config.sensitivity.len() > 1)
        && state.config.global_ranking_targets == 0;
    if temporary {
        backend.perform(WorkflowOp::OpenTemporaryOutput {
            path: state.config.multiprocessing.then(|| {
                get_ref_block_tmpfile_name(
                    &state.config.parallel_tmpdir,
                    state.current_query_block,
                    state.current_reference_block,
                )
            }),
        })?;
    }
    if state.config.global_ranking_targets == 0 {
        backend.perform(WorkflowOp::AlignQueries(if temporary {
            OutputTarget::Temporary
        } else {
            OutputTarget::Master
        }))?;
    }
    if temporary {
        backend.perform(WorkflowOp::FinishIntermediate)?;
    }
    state.target = None;
    backend.perform(WorkflowOp::CloseDictionaryBlock {
        persist: persist_dictionary,
    })?;
    let _ = query;
    Ok(())
}

pub fn run_query_iteration<B: DoubleIndexedBackend>(
    query_iteration: usize,
    state: &mut DoubleIndexedState,
    backend: &mut B,
) -> Result<(), String> {
    let (query_sequence_count, query_source_count, query_letters, query_oid_end) = {
        let query = state.query.as_ref().ok_or("Query block is not loaded.")?;
        (
            query.sequence_count,
            query.source_sequence_count,
            query.letters,
            query.oid_end,
        )
    };
    if query_iteration > 0 && state.query_aligned.len() != query_source_count {
        state.query_aligned.resize(query_source_count, false);
    }
    if query_iteration > 0 {
        backend.perform(WorkflowOp::SetQuerySkipFromAligned)?;
    }
    let sensitivity = state
        .config
        .sensitivity
        .get(query_iteration)
        .ok_or("Query iteration is out of range.")?
        .sensitivity;
    if state.config.algo == Algo::Auto
        && (!state.config.supports_query_indexed.contains(&sensitivity)
            || query_letters > MAX_INDEX_QUERY_SIZE
            || state.database.as_ref().map_or(0, |db| db.letters) < MIN_QUERY_INDEXED_DB_SIZE
            || state.config.target_indexed
            || state.config.swipe_all
            || state.config.minimizer_window != 0
            || state.config.sketch_size != 0)
    {
        state.config.algo = Algo::DoubleIndexed;
    }
    if matches!(state.config.algo, Algo::Auto | Algo::QueryIndexed) {
        backend.perform(WorkflowOp::BuildHashedQuerySeedSet)?;
        let table_size = backend.query_seed_table_size(QuerySeedMode::Hashed)?;
        if state.config.algo == Algo::Auto && !use_query_index(table_size, backend.l3_cache_size())
        {
            state.config.algo = Algo::DoubleIndexed;
            state.query_seed_mode = QuerySeedMode::None;
            backend.perform(WorkflowOp::ClearQuerySeedSets)?;
        } else {
            state.config.algo = Algo::QueryIndexed;
            state.query_seed_mode = QuerySeedMode::Hashed;
            state.config.seed_encoding = SeedEncoding::Hashed;
        }
    }
    if state.config.algo == Algo::CtgSeed {
        backend.perform(WorkflowOp::BuildContiguousQuerySeedSet)?;
        let _ = backend.query_seed_table_size(QuerySeedMode::Contiguous)?;
        state.query_seed_mode = QuerySeedMode::Contiguous;
        state.config.seed_encoding = SeedEncoding::Contiguous;
    }
    let index_chunks = *state.config.index_chunks.get_or_insert_with(|| {
        if state.config.algo == Algo::DoubleIndexed {
            sensitivity_index_chunks(sensitivity)
        } else {
            1
        }
    });
    state.seed_partition_bits = backend.seed_partition_bits(
        state.config.shape_weight,
        state.config.threads,
        index_chunks,
    );
    state.lazy_masking = state.config.algo != Algo::DoubleIndexed
        && state.config.target_masking != MaskingAlgo::None
        && state.config.frame_shift == 0;
    if !state.config.blastn && state.config.gapped_filter_evalue != 0.0 {
        state.cutoff_gapped1 = Some(state.config.gapped_filter_evalue1);
        state.cutoff_gapped2 = Some(state.config.gapped_filter_evalue);
    }
    if state.current_query_block == 0 && query_iteration == 0 {
        backend.perform(WorkflowOp::Message(format!(
            "Algorithm: {:?}",
            state.config.algo
        )))?;
    }
    if state.config.global_ranking_targets > 0 {
        let contexts = state.config.query_contexts.max(1);
        backend.perform(WorkflowOp::AllocateGlobalRankingTable {
            entries: query_sequence_count * state.config.global_ranking_targets / contexts,
        })?;
    }
    if !state.config.swipe_all && !state.config.target_indexed {
        backend.perform(WorkflowOp::BuildQueryHistogram {
            hashed: state.query_seed_mode == QuerySeedMode::Hashed,
        })?;
    }
    backend.perform(WorkflowOp::LogRss)?;
    let mut db_flags = SequenceFileFlags::SEQS;
    if (!state.config.output_titles_lazy && state.config.output_flags.contains(OutputFlags::SSEQID))
        || state.config.no_self_hits
    {
        db_flags |= SequenceFileFlags::TITLES;
    }
    if state.lazy_masking {
        db_flags |= SequenceFileFlags::LAZY_MASKING;
    }
    if state.config.output_flags.contains(OutputFlags::FULL_TITLES) {
        db_flags |= SequenceFileFlags::FULL_TITLES;
    }
    if state.config.output_flags.contains(OutputFlags::ALL_SEQIDS) {
        db_flags |= SequenceFileFlags::ALL_SEQIDS;
    }
    backend.perform(WorkflowOp::SetDatabaseFlags(db_flags))?;

    if state.config.multiprocessing {
        for block in backend.multiprocessing_chunks(state.current_query_block)? {
            backend.perform(WorkflowOp::MultiprocessingSearchBegin {
                query_block: state.current_query_block,
                reference_block: block,
            })?;
            state.current_reference_block = block;
            if !state.config.mp_self || block >= state.current_query_block {
                state.target = backend.load_reference_block(0, Some(block), 0)?;
                if state
                    .target
                    .as_ref()
                    .is_some_and(|target| !target.is_empty())
                {
                    state.blocked_processing = true;
                    run_ref_chunk(query_iteration, state, backend)?;
                }
            } else {
                backend.perform(WorkflowOp::WriteSkippedSelfBlock {
                    query_block: state.current_query_block,
                    reference_block: block,
                })?;
            }
            backend.perform(WorkflowOp::CloseTemporaryOutput)?;
            backend.perform(WorkflowOp::MultiprocessingSearchEnd {
                query_block: state.current_query_block,
                reference_block: block,
            })?;
        }
    } else {
        let mut offset = if state.config.self_search && !state.config.lin_stage1_query {
            query_oid_end
        } else {
            0
        };
        state.current_reference_block = 0;
        loop {
            if state.config.self_search
                && ((state.config.lin_stage1_query
                    && state.current_reference_block == state.current_query_block)
                    || (!state.config.lin_stage1_query && state.current_reference_block == 0))
            {
                state.target = state.query.clone();
                if state.config.lin_stage1_query {
                    offset = query_oid_end;
                }
            } else {
                state.target =
                    backend.load_reference_block(state.config.block_size_letters, None, offset)?;
            }
            let Some(target) = state.target.as_ref() else {
                break;
            };
            if target.is_empty() {
                break;
            }
            if state.current_reference_block == 0 {
                state.blocked_processing = state.config.global_ranking_targets > 0
                    || (target.sequence_count as u64)
                        < state.database.as_ref().map_or(0, |db| db.sequence_count);
            }
            run_ref_chunk(query_iteration, state, backend)?;
            state.current_reference_block += 1;
        }
    }
    backend.perform(WorkflowOp::ClearQuerySeedSets)?;
    state.query_seed_mode = QuerySeedMode::None;
    if state.config.global_ranking_targets > 0 {
        let temporary = state.config.sensitivity.len() > 1;
        if temporary {
            backend.perform(WorkflowOp::OpenTemporaryOutput { path: None })?;
        }
        backend.perform(WorkflowOp::ExtendGlobalRanking(if temporary {
            OutputTarget::Temporary
        } else {
            OutputTarget::Master
        }))?;
    }
    Ok(())
}

pub fn run_query_chunk<B: DoubleIndexedBackend>(
    state: &mut DoubleIndexedState,
    backend: &mut B,
) -> Result<(), String> {
    let source_count = state
        .query
        .as_ref()
        .ok_or("Query block is not loaded.")?
        .source_sequence_count as u64;
    if state.config.track_aligned_queries {
        state.query_aligned = vec![false; source_count as usize];
    }
    if state
        .config
        .output_flags
        .contains(OutputFlags::SELF_ALN_SCORES)
    {
        backend.perform(WorkflowOp::ComputeQuerySelfAlignment)?;
    }
    backend.perform(WorkflowOp::LogRss)?;
    state.aligned_total = 0;
    for query_iteration in 0..state.config.sensitivity.len() {
        if state.aligned_total >= source_count {
            break;
        }
        let round = state.config.sensitivity[query_iteration];
        state.lin_stage1_target = state.config.lin_stage1_target || round.linearize;
        let linearization_count = usize::from(state.lin_stage1_target)
            + usize::from(state.config.lin_stage1_query)
            + usize::from(state.config.lin_stage1_combo);
        if linearization_count > 1 {
            return Err(
                "Multiple linearization options are not allowed to be used together".into(),
            );
        }
        backend.perform(WorkflowOp::SetupSearch(round.sensitivity))?;
        run_query_iteration(query_iteration, state, backend)?;
        if state.config.sensitivity.len() > 1 {
            state.aligned_total += backend.take_iteration_aligned();
        }
    }
    backend.perform(WorkflowOp::LogRss)?;
    if state.blocked_processing
        || state.config.multiprocessing
        || state.config.sensitivity.len() > 1
    {
        backend.perform(WorkflowOp::JoinBlocks {
            reference_blocks: if state.config.multiprocessing {
                backend.reference_partition_chunks()
            } else {
                state.current_reference_block
            },
            multiprocessing: state.config.multiprocessing,
        })?;
    }
    if state.config.track_aligned_queries {
        backend.perform(WorkflowOp::WriteUnalignedQueries)?;
        backend.perform(WorkflowOp::WriteAlignedQueries)?;
    }
    state.query = None;
    Ok(())
}

pub fn master_thread<B: DoubleIndexedBackend>(
    state: &mut DoubleIndexedState,
    backend: &mut B,
) -> Result<(), String> {
    backend.perform(WorkflowOp::LogRss)?;
    if state.config.multiprocessing && state.config.mp_recover {
        backend.perform(WorkflowOp::RecoverMultiprocessing {
            max_query_chunks: 65_536,
        })?;
        if backend.stop_exists() {
            backend.remove_stop()?;
        }
        return Ok(());
    }
    if state.config.multiprocessing {
        backend.perform(WorkflowOp::InitMultiprocessing {
            directory: state.config.parallel_tmpdir.clone(),
            chunk_letters: (state.config.chunk_size_billions.unwrap_or(0.0) * 1e9) as u64,
        })?;
    }
    let mut query_flags = SequenceFileFlags::ALL;
    if state.config.output_flags.contains(OutputFlags::FULL_TITLES) {
        query_flags |= SequenceFileFlags::FULL_TITLES;
    }
    if state.config.output_flags.contains(OutputFlags::ALL_SEQIDS) {
        query_flags |= SequenceFileFlags::ALL_SEQIDS;
    }
    if state.config.store_query_quality {
        query_flags |= SequenceFileFlags::QUALITY;
    }
    if state.config.multiprocessing && state.config.mp_init {
        let mut block_count = 0usize;
        let mut offset = 0u64;
        loop {
            let block = backend.load_query_block(
                (state.config.chunk_size_billions.unwrap_or(0.0) * 1e9) as u64,
                state.config.self_search,
                offset,
                query_flags,
            )?;
            let Some(block) = block else { break };
            if block.is_empty() {
                break;
            }
            offset = block.oid_end;
            block_count += 1;
        }
        backend.rewind_query()?;
        for query_block in 0..block_count {
            backend.perform(WorkflowOp::SavePartition {
                query_block,
                path: get_ref_part_file_name(
                    &state.config.parallel_tmpdir,
                    STACK_ALIGN_TODO,
                    query_block,
                    "",
                ),
            })?;
        }
        return Ok(());
    }
    if state.config.query_file_missing && !state.config.self_search {
        backend.perform(WorkflowOp::Message(
            "Query file parameter (--query/-q) is missing. Input will be read from stdin."
                .to_owned(),
        ))?;
    }
    backend.perform(WorkflowOp::OpenOutput)?;
    if state.config.output_is_daa {
        backend.perform(WorkflowOp::InitDaa)?;
    }
    let mut offset = 0u64;
    state.current_query_block = 0;
    loop {
        let block = backend.load_query_block(
            state.config.block_size_letters,
            state.config.self_search,
            offset,
            query_flags,
        )?;
        let Some(block) = block else { break };
        if block.is_empty() {
            break;
        }
        offset = block.oid_end;
        state.query = Some(block);
        if state
            .config
            .mp_query_chunk
            .is_some_and(|selected| selected != state.current_query_block)
        {
            state.current_query_block += 1;
            continue;
        }
        if (!keep_target_id(state) && state.config.lin_stage1_query && !state.config.kmer_ranking)
            || state.config.min_length_ratio > 0.0
        {
            backend.perform(WorkflowOp::SortQueryByLength)?;
        }
        if state.current_query_block == 0
            && !state.config.output_is_daa
            && state.query.as_ref().is_some_and(|query| query.has_ids)
        {
            backend.perform(WorkflowOp::PrintHeader)?;
        }
        if state.config.query_masking != MaskingAlgo::None {
            backend.perform(WorkflowOp::MaskQuery(state.config.query_masking))?;
        }
        run_query_chunk(state, backend)?;
        if backend.stop_exists() {
            break;
        }
        state.current_query_block += 1;
    }
    if !state.config.self_search {
        backend.perform(WorkflowOp::CloseQueryInput)?;
    }
    if state.config.output_is_daa {
        backend.perform(WorkflowOp::FinishDaa {
            query_blocks: state.current_query_block,
        })?;
    } else {
        backend.perform(WorkflowOp::PrintFooter)?;
    }
    backend.perform(WorkflowOp::FinalizeOutput)?;
    backend.perform(WorkflowOp::CloseAuxiliaryOutputs)?;
    if !state.config.unaligned_targets.is_empty() {
        backend.perform(WorkflowOp::WriteUnalignedTargets(
            state.config.unaligned_targets.clone(),
        ))?;
    }
    backend.perform(WorkflowOp::CloseDatabase)?;
    backend.perform(WorkflowOp::FreeConfig)?;
    backend.perform(WorkflowOp::StatisticsPrint)?;
    Ok(())
}

pub fn run<B: DoubleIndexedBackend>(
    mut config: DoubleIndexedConfig,
    backend: &mut B,
) -> Result<DoubleIndexedState, String> {
    if config.chunk_size_billions.is_none() {
        let final_sensitivity = config
            .sensitivity
            .last()
            .map_or(Sensitivity::Default, |round| round.sensitivity);
        config.chunk_size_billions = Some(if final_sensitivity >= Sensitivity::VerySensitive {
            0.4
        } else {
            2.0
        });
    }
    backend.perform(WorkflowOp::Message(format!(
        "Temporary directory: {}",
        config.parallel_tmpdir.display()
    )))?;
    backend.perform(WorkflowOp::StatisticsReset)?;
    let taxon_filter = !config.taxonlist.is_empty() || !config.taxon_exclude.is_empty();
    let mut flags = SequenceFileFlags::NEED_LETTER_COUNT;
    if config.output_needs_taxon_id_lists || taxon_filter || config.taxon_culling {
        flags |= SequenceFileFlags::TAXON_MAPPING;
    }
    if config.output_needs_taxon_nodes || taxon_filter || config.taxon_culling {
        flags |= SequenceFileFlags::TAXON_NODES;
    }
    if config.output_needs_taxon_scientific_names {
        flags |= SequenceFileFlags::TAXON_SCIENTIFIC_NAMES;
    }
    if config.output_needs_taxon_ranks || config.taxon_culling {
        flags |= SequenceFileFlags::TAXON_RANKS;
    }
    if config.output_flags.contains(OutputFlags::ALL_SEQIDS) {
        flags |= SequenceFileFlags::ALL_SEQIDS;
    }
    if config.output_flags.contains(OutputFlags::FULL_TITLES) || config.no_self_hits {
        flags |= SequenceFileFlags::FULL_TITLES;
    }
    if config.output_flags.contains(OutputFlags::TARGET_SEQS) {
        flags |= SequenceFileFlags::TARGET_SEQS;
    }
    if config.output_flags.contains(OutputFlags::SELF_ALN_SCORES) {
        flags |= SequenceFileFlags::SELF_ALN_SCORES;
    }
    if !config.unaligned_targets.is_empty() {
        flags |= SequenceFileFlags::OID_TO_ACC_MAPPING;
    }
    if taxon_filter {
        flags |=
            SequenceFileFlags::NEED_EARLY_TAXON_MAPPING | SequenceFileFlags::NEED_LENGTH_LOOKUP;
    }
    if !config.seqidlist.is_empty() {
        flags |= SequenceFileFlags::NEED_LENGTH_LOOKUP;
    }
    let database = backend.database_info(flags)?;
    if config.multiprocessing && database.file_type == SequenceFileType::Fasta {
        return Err("Multiprocessing mode is not compatible with FASTA databases.".into());
    }
    backend.perform(WorkflowOp::PrintDatabaseInfo)?;
    let alias_taxfilter = database.alias_taxidlist.clone();
    if taxon_filter {
        if !config.taxonlist.is_empty() && !config.taxon_exclude.is_empty() {
            return Err("Options --taxonlist and --taxon-exclude are mutually exclusive.".into());
        }
        backend.perform(WorkflowOp::ApplyTaxonomyFilter {
            values: if config.taxonlist.is_empty() {
                config.taxon_exclude.clone()
            } else {
                config.taxonlist.clone()
            },
            exclude: !config.taxon_exclude.is_empty(),
            separator: ',',
        })?;
    } else if let Some(path) = alias_taxfilter {
        backend.perform(WorkflowOp::ApplyTaxonomyFilter {
            values: path,
            exclude: false,
            separator: '\n',
        })?;
    }
    let mut seqidlist = config.seqidlist.clone();
    if let Some(alias) = &database.alias_seqidlist {
        if !seqidlist.is_empty() {
            return Err("Using --seqidlist on already filtered BLAST alias database.".into());
        }
        seqidlist.clone_from(alias);
    }
    if !seqidlist.is_empty() {
        if taxon_filter {
            return Err("--seqidlist is not compatible with taxonomy filtering.".into());
        }
        backend.perform(WorkflowOp::ApplyAccessionFilter(seqidlist))?;
    }
    let filtered_letters = backend
        .filtered_database_totals()
        .map(|(_, letters)| letters)
        .filter(|letters| *letters != 0);
    backend.perform(WorkflowOp::SetScoreMatrixDatabaseLetters(
        config
            .db_size
            .or(filtered_letters)
            .unwrap_or(database.letters),
    ))?;
    if config.blastn {
        backend.perform(WorkflowOp::ConfigureNucleotideScoring)?;
    }
    let mut state = DoubleIndexedState::new(config);
    state.database = Some(database);
    master_thread(&mut state, backend)?;
    backend.perform(WorkflowOp::LogRss)?;
    Ok(state)
}

fn sensitivity_index_chunks(sensitivity: Sensitivity) -> u32 {
    match sensitivity {
        Sensitivity::Faster | Sensitivity::Fast => 1,
        Sensitivity::Default | Sensitivity::Linclust40 | Sensitivity::Linclust20 => 2,
        Sensitivity::Shapes6x10 | Sensitivity::Shapes30x10 | Sensitivity::MidSensitive => 3,
        Sensitivity::Sensitive | Sensitivity::MoreSensitive => 4,
        Sensitivity::VerySensitive | Sensitivity::UltraSensitive => 8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::VecDeque;

    struct Backend {
        events: Vec<WorkflowOp>,
        queries: VecDeque<SequenceBlock>,
        references: VecDeque<SequenceBlock>,
        iterations: VecDeque<u64>,
        seed_table_size: usize,
        l3: usize,
        stop: bool,
        database: DatabaseInfo,
    }

    impl Default for Backend {
        fn default() -> Self {
            Self {
                events: Vec::new(),
                queries: VecDeque::new(),
                references: VecDeque::new(),
                iterations: VecDeque::new(),
                seed_table_size: 1,
                l3: 1,
                stop: false,
                database: DatabaseInfo {
                    file_type: SequenceFileType::Dmnd,
                    sequence_count: 2,
                    letters: MIN_QUERY_INDEXED_DB_SIZE,
                    total_blocks: 1,
                    file_name: "db.dmnd".into(),
                    alias_taxidlist: None,
                    alias_seqidlist: None,
                },
            }
        }
    }

    impl DoubleIndexedBackend for Backend {
        fn l3_cache_size(&self) -> usize {
            self.l3
        }

        fn database_info(&mut self, _: SequenceFileFlags) -> Result<DatabaseInfo, String> {
            Ok(self.database.clone())
        }

        fn filtered_database_totals(&self) -> Option<(u64, u64)> {
            None
        }

        fn load_query_block(
            &mut self,
            _: u64,
            _: bool,
            _: u64,
            _: SequenceFileFlags,
        ) -> Result<Option<SequenceBlock>, String> {
            Ok(self.queries.pop_front())
        }

        fn rewind_query(&mut self) -> Result<(), String> {
            Ok(())
        }

        fn load_reference_block(
            &mut self,
            _: u64,
            _: Option<usize>,
            _: u64,
        ) -> Result<Option<SequenceBlock>, String> {
            Ok(self.references.pop_front())
        }

        fn reference_partition_chunks(&self) -> usize {
            self.database.total_blocks
        }

        fn query_seed_table_size(&mut self, _: QuerySeedMode) -> Result<usize, String> {
            Ok(self.seed_table_size)
        }

        fn histogram_entries(&self, query: bool) -> usize {
            if query {
                11
            } else {
                17
            }
        }

        fn seed_partition_bits(&self, weight: u32, threads: usize, chunks: u32) -> u32 {
            weight + threads as u32 + chunks
        }

        fn multiprocessing_chunks(&mut self, _: usize) -> Result<Vec<usize>, String> {
            Ok((0..self.database.total_blocks).collect())
        }

        fn take_iteration_aligned(&mut self) -> u64 {
            self.iterations.pop_front().unwrap_or(0)
        }

        fn stop_exists(&self) -> bool {
            self.stop
        }

        fn remove_stop(&mut self) -> Result<(), String> {
            self.stop = false;
            Ok(())
        }

        fn perform(&mut self, operation: WorkflowOp) -> Result<(), String> {
            self.events.push(operation);
            Ok(())
        }
    }

    fn block(begin: u64, count: usize) -> SequenceBlock {
        SequenceBlock {
            sequence_count: count,
            source_sequence_count: count,
            letters: count as u64 * 100,
            raw_bytes: count as u64 * 120,
            oid_begin: begin,
            oid_end: begin + count as u64,
            has_ids: true,
        }
    }

    #[test]
    fn thresholds_and_file_names_match_upstream() {
        assert!(use_query_index(MAX_HASH_SET_SIZE, 1));
        assert!(use_query_index(20, 20));
        assert!(!use_query_index(MAX_HASH_SET_SIZE + 1, 1));
        let root = Path::new("tmp");
        assert_eq!(
            get_ref_part_file_name(root, STACK_ALIGN_TODO, 7, "redo"),
            PathBuf::from("tmp/align_todo_redo_7")
        );
        assert_eq!(
            get_ref_block_tmpfile_name(root, 3, 9),
            PathBuf::from("tmp/ref_block_3_9")
        );
    }

    #[test]
    fn query_iteration_selects_query_index_and_runs_reference_lifecycle() {
        let mut config = DoubleIndexedConfig {
            algo: Algo::Auto,
            target_masking: MaskingAlgo::Tantan,
            shape_count: 2,
            shape_weight: 5,
            threads: 3,
            output_flags: OutputFlags::TARGET_SEQS | OutputFlags::SELF_ALN_SCORES,
            ..DoubleIndexedConfig::default()
        };
        config.index_chunks = None;
        let mut state = DoubleIndexedState::new(config);
        state.database = Some(Backend::default().database);
        state.query = Some(block(0, 2));
        let mut backend = Backend {
            seed_table_size: MAX_HASH_SET_SIZE,
            references: VecDeque::from([block(0, 2)]),
            ..Backend::default()
        };
        run_query_iteration(0, &mut state, &mut backend).unwrap();
        assert_eq!(state.config.algo, Algo::QueryIndexed);
        assert_eq!(state.config.seed_encoding, SeedEncoding::Hashed);
        assert!(state.lazy_masking);
        assert_eq!(state.seed_partition_bits, 9);
        assert!(backend
            .events
            .contains(&WorkflowOp::BuildHashedQuerySeedSet));
        assert!(backend
            .events
            .contains(&WorkflowOp::BuildReferenceHistogram(
                HistogramFilter::HashedQuerySeeds
            )));
        assert_eq!(
            backend
                .events
                .iter()
                .filter(|event| matches!(event, WorkflowOp::SearchShape { .. }))
                .count(),
            2
        );
        assert!(backend.events.contains(&WorkflowOp::FinishHitBuffer));
        assert!(backend
            .events
            .contains(&WorkflowOp::AlignQueries(OutputTarget::Master)));
    }

    #[test]
    fn oversized_query_forces_double_indexed_and_gapped_cutoffs() {
        let config = DoubleIndexedConfig {
            algo: Algo::Auto,
            gapped_filter_evalue: 0.01,
            gapped_filter_evalue1: 0.1,
            ..DoubleIndexedConfig::default()
        };
        let mut state = DoubleIndexedState::new(config);
        state.database = Some(Backend::default().database);
        let mut query = block(0, 1);
        query.letters = MAX_INDEX_QUERY_SIZE + 1;
        state.query = Some(query);
        let mut backend = Backend {
            references: VecDeque::from([block(0, 2)]),
            ..Backend::default()
        };
        run_query_iteration(0, &mut state, &mut backend).unwrap();
        assert_eq!(state.config.algo, Algo::DoubleIndexed);
        assert_eq!(state.config.index_chunks, Some(2));
        assert!(!state.lazy_masking);
        assert_eq!(state.cutoff_gapped1, Some(0.1));
        assert_eq!(state.cutoff_gapped2, Some(0.01));
        assert!(!backend
            .events
            .contains(&WorkflowOp::BuildHashedQuerySeedSet));
    }

    #[test]
    fn iterated_global_ranking_extends_then_joins_once() {
        let config = DoubleIndexedConfig {
            algo: Algo::DoubleIndexed,
            global_ranking_targets: 4,
            query_contexts: 2,
            sensitivity: vec![
                Round::new(Sensitivity::Fast, false),
                Round::new(Sensitivity::Default, false),
            ],
            ..DoubleIndexedConfig::default()
        };
        let mut state = DoubleIndexedState::new(config);
        state.database = Some(Backend::default().database);
        state.query = Some(block(0, 2));
        let mut backend = Backend {
            references: VecDeque::from([block(0, 2), SequenceBlock::empty(), block(0, 2)]),
            iterations: VecDeque::from([1, 1]),
            ..Backend::default()
        };
        run_query_chunk(&mut state, &mut backend).unwrap();
        assert_eq!(state.aligned_total, 2);
        assert_eq!(
            backend
                .events
                .iter()
                .filter(|event| matches!(event, WorkflowOp::ExtendGlobalRanking(_)))
                .count(),
            2
        );
        assert!(backend.events.contains(&WorkflowOp::JoinBlocks {
            reference_blocks: 1,
            multiprocessing: false,
        }));
        assert!(backend
            .events
            .contains(&WorkflowOp::AllocateGlobalRankingTable { entries: 4 }));
    }

    #[test]
    fn top_level_run_applies_filters_and_finalizes_output() {
        let config = DoubleIndexedConfig {
            parallel_tmpdir: PathBuf::from("work"),
            algo: Algo::DoubleIndexed,
            query_masking: MaskingAlgo::Tantan,
            taxonlist: "2,9606".into(),
            output_flags: OutputFlags::FULL_TITLES | OutputFlags::ALL_SEQIDS,
            track_aligned_queries: true,
            ..DoubleIndexedConfig::default()
        };
        let mut backend = Backend {
            queries: VecDeque::from([block(0, 2)]),
            references: VecDeque::from([block(0, 2)]),
            ..Backend::default()
        };
        let state = run(config, &mut backend).unwrap();
        assert_eq!(state.config.chunk_size_billions, Some(2.0));
        assert!(backend.events.contains(&WorkflowOp::ApplyTaxonomyFilter {
            values: "2,9606".into(),
            exclude: false,
            separator: ',',
        }));
        assert!(backend
            .events
            .contains(&WorkflowOp::MaskQuery(MaskingAlgo::Tantan)));
        assert!(backend.events.contains(&WorkflowOp::PrintHeader));
        assert!(backend.events.contains(&WorkflowOp::PrintFooter));
        assert!(backend.events.contains(&WorkflowOp::FinalizeOutput));
        assert!(backend.events.contains(&WorkflowOp::WriteUnalignedQueries));
        assert!(backend.events.contains(&WorkflowOp::WriteAlignedQueries));
        assert!(backend.events.contains(&WorkflowOp::CloseDatabase));
    }

    #[test]
    fn multiprocessing_recovery_removes_stop_and_returns_without_output() {
        let config = DoubleIndexedConfig {
            multiprocessing: true,
            mp_recover: true,
            ..DoubleIndexedConfig::default()
        };
        let mut state = DoubleIndexedState::new(config);
        let mut backend = Backend {
            stop: true,
            ..Backend::default()
        };
        master_thread(&mut state, &mut backend).unwrap();
        assert!(!backend.stop);
        assert_eq!(
            backend.events,
            [
                WorkflowOp::LogRss,
                WorkflowOp::RecoverMultiprocessing {
                    max_query_chunks: 65_536
                }
            ]
        );
    }
}
