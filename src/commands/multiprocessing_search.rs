//! Real native MPP adapter for the Source reference-block work protocol.
//! Nested under search_blocks so query/reference owners and normal output
//! renderers are reused rather than reconstructed through a second pipeline.

use crate::util::io::source_c_path;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use std::sync::Arc;

use super::{QueryChunk, QueryInput, ReferenceDictionary, ReferenceInfo, ReferenceSource};
use crate::commands::blastp::{BlastpConfig, NativeSearchWriter};
use crate::commands::search_output::SearchSubjects;
use crate::data::persistent_dictionary::{
    dictionary_file_name, load_dictionary_bytes as load_dictionary, DictionaryLayout,
    DictionaryWriter,
};
use crate::data::sequence_file::ByteDictionaryEntry as DictionaryEntry;
use crate::output::intermediate::{IntermediateRecord, OutputFormat};
use crate::run::double_indexed::{
    get_ref_block_tmpfile_name, get_ref_part_file_name, STACK_ALIGN_TODO,
};
use crate::run::multiprocessing::{stop_exists, MultiprocessingQueues};

fn needs_self_scores(config: &BlastpConfig) -> bool {
    config
        .outfmt
        .iter()
        .any(|field| field == b"normalized_bitscore")
        || config
            .tuning
            .output
            .clustering_filter
            .as_ref()
            .is_some_and(|filter| {
                filter
                    .format
                    .flags
                    .any(crate::output::format::OutputFlags::SELF_ALN_SCORES)
            })
}

fn dictionary_layout(config: &BlastpConfig, writer: &NativeSearchWriter) -> DictionaryLayout {
    DictionaryLayout {
        lengths: true,
        titles: true,
        sequences: writer.needs_reference_sequences(),
        self_scores: needs_self_scores(config),
    }
}

fn save_dictionary(
    path: &Path,
    dictionary: &ReferenceDictionary,
    layout: DictionaryLayout,
    buffer: usize,
) -> io::Result<()> {
    let mut output = DictionaryWriter::create(path, layout, buffer)?;
    output.init_block(0, dictionary.entries.len(), false);
    for (id, entry) in dictionary.entries.iter().enumerate() {
        output.dict_id_bytes(
            0,
            id,
            &DictionaryEntry {
                oid: entry.oid,
                len: entry.len as i32,
                title: entry.title.clone(),
                seq: entry.sequence.as_deref().unwrap_or(&[]).to_vec(),
                self_aln_score: entry.self_aln_score.unwrap_or(0.0),
            },
        )?;
    }
    output.close_block(false)?;
    output.finish()
}

/// Translate local per-file dictionary IDs into the actual joined metadata
/// arena while retaining every other serialized record byte unchanged.
fn remap_payload(payload: &mut [u8], ids: &[usize], format: &OutputFormat) -> io::Result<()> {
    let mut position = 0;
    while position < payload.len() {
        let (local, bytes) = {
            let mut input = io::Cursor::new(&payload[position..]);
            let record = IntermediateRecord::read(&mut input, format)?;
            (record.target_dict_id, input.position() as usize)
        };
        let id = usize::try_from(local)
            .ok()
            .and_then(|id| ids.get(id))
            .ok_or_else(|| io::Error::other("Dictionary not loaded."))?;
        payload[position..position + 8].copy_from_slice(&(*id as i64).to_ne_bytes());
        position += bytes;
    }
    Ok(())
}

pub(super) fn sort_join_only_queries(
    chunk: &mut QueryChunk,
    config: &BlastpConfig,
) -> io::Result<Vec<u64>> {
    if config.tuning.query_pre_sorted {
        return Ok(Vec::new());
    }
    let ratio = if let Some(ratio) = config.tuning.search_min_length_ratio {
        ratio
    } else {
        crate::run::config::derive_min_length_ratio(
            &crate::run::config::ConfigInput {
                min_length_ratio: config.tuning.min_length_ratio,
                query_cover: config.query_cover,
                subject_cover: config.subject_cover,
                query_translated: config.translated_query_layout.is_some(),
                lin_stage1_query: config.tuning.setup.lin_stage1_query,
                ..Default::default()
            },
            &[crate::run::config::Round::new(config.sensitivity, false)],
        )
        .map_err(io::Error::other)?
    };
    super::prepare_master_query_chunk(chunk, config, config.tuning.search_block.reference, ratio)
}

/// Even a worker claiming no reference jobs prepares its query owner before
/// joining: Source masks queries and computes their self scores unconditionally.
pub(super) fn prepare_join_only_queries(
    chunk: &mut QueryChunk,
    config: &BlastpConfig,
    matrix: &crate::stats::score_matrix::ScoreMatrix,
    masker: Option<&crate::masking::tantan::TantanMasker>,
) {
    // The caller may share the TANTAN kernel for soft seed masking even
    // when Source Search::Config disables query hard masking (clustering).
    // Kernel availability alone must not replace residues with X here.
    // Source SEG masking applies only to targets, not queries.
    let masker = masker.filter(|_| config.masking == crate::masking::MaskingMode::Tantan);
    let prepare = |sequence: &mut [crate::basic::value::Letter]| {
        crate::masking::remove_bit_mask(sequence);
        if let Some(masker) = masker {
            masker.mask_hard(sequence);
        }
    };
    if let Some(block) = chunk.translated.as_mut() {
        for frame in 0..block.seqs().len() {
            prepare(block.seqs_mut().get_mut(frame));
        }
        if needs_self_scores(config) {
            chunk.self_aln_scores =
                crate::commands::search_output::self_alignment_scores_with_window(
                    block.seqs(),
                    matrix,
                    config.comp_based_stats,
                    config.tuning.cbs_window,
                );
        }
    } else {
        for record in &mut chunk.records {
            prepare(&mut record.sequence);
        }
        if needs_self_scores(config) {
            let mut sequences = crate::data::sequence_set::SequenceSet::new();
            for record in &chunk.records {
                sequences.push(&record.sequence);
            }
            chunk.self_aln_scores =
                crate::commands::search_output::self_alignment_scores_with_window(
                    &sequences,
                    matrix,
                    config.comp_based_stats,
                    config.tuning.cbs_window,
                );
        }
    }
}

pub(super) fn run(
    config: &BlastpConfig,
    references: ReferenceSource,
    input: QueryInput<'_>,
    budget: u64,
) -> io::Result<()> {
    run_retaining_references(config, references, input, budget).map(drop)
}

pub(super) fn run_retaining_references(
    config: &BlastpConfig,
    references: ReferenceSource,
    input: QueryInput<'_>,
    budget: u64,
) -> io::Result<ReferenceSource> {
    run_impl(config, references, input, budget, None, 0.0, None, None)
}

/// The supplied clustering consumer is bypassed by Source MPP, but the
/// ordinary queue, partition, extension, dictionary and chunk output work remains.
pub(super) fn run_cluster(
    config: &BlastpConfig,
    references: ReferenceSource,
    input: QueryInput<'_>,
    budget: u64,
    options: crate::commands::blastp::ClusterSearchOptions<'_>,
    owner: &super::cluster_search::ClusterSearchOwner,
) -> io::Result<()> {
    run_impl(
        config,
        references,
        input,
        budget,
        None,
        0.0,
        Some(options),
        Some(owner),
    )
    .map(drop)
}

pub(super) fn run_dna(
    config: &BlastpConfig,
    references: ReferenceSource,
    input: QueryInput<'_>,
    budget: u64,
    dna: &crate::dna::native::DnaOptions,
    dna_min_length_ratio: f64,
) -> io::Result<()> {
    run_impl(
        config,
        references,
        input,
        budget,
        Some(dna),
        dna_min_length_ratio,
        None,
        None,
    )
    .map(drop)
}

fn run_impl(
    config: &BlastpConfig,
    mut references: ReferenceSource,
    mut input: QueryInput<'_>,
    budget: u64,
    dna: Option<&crate::dna::native::DnaOptions>,
    dna_min_length_ratio: f64,
    cluster: Option<crate::commands::blastp::ClusterSearchOptions<'_>>,
    cluster_owner: Option<&super::cluster_search::ClusterSearchOwner>,
) -> io::Result<ReferenceSource> {
    crate::data::accession_spool::source_scope(references.accession_spool(), || {
        let mut source_config = config.clone();
        if dna.is_none() {
            crate::commands::blastp::begin_algorithm_scope(&mut source_config)?;
        }
        let config = &source_config;
        let mp = &config.tuning.multiprocessing;
        if mp.directory.as_os_str().is_empty() {
            return Err(io::Error::other(
                "--multiprocessing requires setting --parallel-tmpdir",
            ));
        }
        if config.iterate.is_some() {
            return Err(io::Error::other(
                "Iterated search is not compatible with --multiprocessing.",
            ));
        }
        if !config.taxonlist.is_empty() || !config.taxon_exclude.is_empty() {
            return Err(io::Error::other(
                "Multiprocessing mode is not compatible with database filtering.",
            ));
        }
        match &references {
            ReferenceSource::Database {
                _temporary: None, ..
            } => (),
            ReferenceSource::Blast { .. } => {
                return Err(io::Error::other(
                    "--multiprocessing is not compatible with BLAST databases.",
                ))
            }
            _ => {
                return Err(io::Error::other(
                    "Multiprocessing mode is not compatible with FASTA databases.",
                ))
            }
        }
        let errors = references.source_errors();
        if mp.recover {
            let _ = crate::util::system::log_rss();
            return crate::run::multiprocessing::recover_owned(&mp.directory, &errors)
                .map(|()| references);
        }
        let queues = MultiprocessingQueues::new_owned(&mp.directory, errors)?;
        let partition = if let ReferenceSource::Database { file, .. } = &mut references {
            file.create_partition_balanced(i64::try_from(budget).map_err(io::Error::other)?)
                .map_err(io::Error::other)?;
            file.partition.chunks.clone()
        } else {
            unreachable!()
        };
        let check_query_protein = input.checks_protein_input();
        if mp.initialize {
            let _ = crate::util::system::log_rss();
            let mut query = 0;
            loop {
                let chunk = input.chunk(budget)?;
                if chunk.count() == 0 {
                    break;
                }
                if dna.is_none() && check_query_protein && query == 0 && chunk.translated.is_none()
                {
                    crate::data::sequence_file::check_protein_input(
                        chunk
                            .records
                            .iter()
                            .map(|record| record.sequence.as_slice()),
                        config.tuning.ignore_warnings,
                    )
                    .map_err(io::Error::other)?;
                }
                query += 1;
            }
            input.rewind_after_count()?;
            for query in 0..query {
                if let ReferenceSource::Database { file, .. } = &references {
                    file.save_partition(
                        get_ref_part_file_name(&mp.directory, STACK_ALIGN_TODO, query, ""),
                        &format!("# query_chunk={query}"),
                    )
                    .map_err(io::Error::other)?;
                }
            }
            return Ok(references);
        }
        let (sequence_count, full_letters) = references.totals();
        let unaligned_targets = super::UnalignedTargets::new(config, &mut references)?;
        let dna_statistics = dna
            .map(|options| {
                crate::dna::build_score::Blastn_Score::new(
                    options.reward,
                    options.penalty,
                    if config.gap_open == -1 {
                        5
                    } else {
                        config.gap_open
                    },
                    if config.gap_extend == -1 {
                        2
                    } else {
                        config.gap_extend
                    },
                    full_letters,
                    i64::try_from(sequence_count).map_err(|_| {
                        "Database sequence count exceeds nucleotide statistics domain".to_owned()
                    })?,
                )
            })
            .transpose()
            .map_err(io::Error::other)?;
        let matrix = Arc::new(if let Some(statistics) = dna_statistics {
            let mut matrix =
                crate::stats::score_matrix::ScoreMatrix::source_blastn_metadata(statistics)
                    .map_err(io::Error::other)?;
            matrix.set_db_letters(if config.tuning.db_size != 0 {
                config.tuning.db_size
            } else {
                references.scoring_letters()
            });
            matrix
        } else {
            crate::commands::blastp::native_score_matrix(config, references.scoring_letters())?
        });
        let masker = (dna.is_none() && config.masking == crate::masking::MaskingMode::Tantan)
            .then(|| crate::commands::blastp::native_tantan_masker(config, &matrix));
        let query_mode = if dna.is_some() {
            4
        } else if config.translated_query_layout.is_some() {
            3
        } else {
            2
        };
        let min_length_ratio = if dna.is_some() {
            dna_min_length_ratio
        } else {
            config.tuning.search_min_length_ratio.unwrap_or(
                crate::run::config::derive_min_length_ratio(
                    &crate::run::config::ConfigInput {
                        min_length_ratio: config.tuning.min_length_ratio,
                        query_cover: config.query_cover,
                        subject_cover: config.subject_cover,
                        query_translated: config.translated_query_layout.is_some(),
                        lin_stage1_query: config.tuning.setup.lin_stage1_query,
                        ..Default::default()
                    },
                    &[crate::run::config::Round::new(config.sensitivity, false)],
                )
                .map_err(io::Error::other)?,
            )
        };
        crate::run::double_indexed::owned::run(MppBackend {
            config,
            references,
            input,
            budget,
            dna,
            dna_min_length_ratio,
            dna_statistics,
            cluster,
            cluster_owner,
            queues,
            partition,
            sequence_count,
            full_letters,
            matrix,
            masker,
            query_mode,
            output: None,
            unaligned_targets,
            min_length_ratio,
            check_query_protein,
            query_begin: 0,
            source_reference_block: 0,
        })
    })
}

/// Real loaded query owner and the Source work-stack handles for its lifetime.
struct MppQuery {
    chunk: QueryChunk,
    ordinal: usize,
    begin: u64,
    master_oids: Vec<u64>,
    output_oids: Vec<u64>,
    aligned: Vec<bool>,
    prepared: bool,
    master_prepared: bool,
    block_config: BlastpConfig,
    work: Option<crate::run::multiprocessing::QueryQueues>,
}

/// A claimed queue record stays in WIP until data/dictionary finalization
/// succeeds. Dropping a failed owner never acknowledges or reclaims the job.
struct ClaimedReference {
    record: String,
    chunk: crate::data::dmnd::Chunk,
}

/// Output resources are constructed by begin_master only after query input
/// construction succeeds, matching Source's error/truncation ordering.
struct MppOutput {
    main: NativeSearchWriter,
    format: OutputFormat,
    layout: DictionaryLayout,
    auxiliary: crate::commands::query_output::QueryOutputs,
}

struct MppBackend<'cfg, 'input, 'options> {
    config: &'cfg BlastpConfig,
    references: ReferenceSource,
    input: QueryInput<'input>,
    budget: u64,
    dna: Option<&'cfg crate::dna::native::DnaOptions>,
    dna_min_length_ratio: f64,
    dna_statistics: Option<crate::dna::build_score::Blastn_Score>,
    cluster: Option<crate::commands::blastp::ClusterSearchOptions<'options>>,
    cluster_owner: Option<&'cfg super::cluster_search::ClusterSearchOwner>,
    queues: MultiprocessingQueues,
    partition: Vec<crate::data::dmnd::Chunk>,
    sequence_count: u64,
    full_letters: u64,
    matrix: Arc<crate::stats::score_matrix::ScoreMatrix>,
    masker: Option<Arc<crate::masking::tantan::TantanMasker>>,
    query_mode: i32,
    output: Option<MppOutput>,
    unaligned_targets: Option<super::UnalignedTargets>,
    min_length_ratio: f64,
    check_query_protein: bool,
    query_begin: u64,
    source_reference_block: u32,
}

impl crate::run::double_indexed::owned::DoubleIndexedBackend for MppBackend<'_, '_, '_> {
    type Query = MppQuery;
    type Round = ();
    type Reference = ClaimedReference;
    type Output = ReferenceSource;

    fn capture_source_error(&self, error: io::Error) -> io::Error {
        let errors = self.references.source_errors();
        errors.capture_io::<()>(Err(error)).unwrap_err()
    }
    fn database_sequence_count(&self) -> u64 {
        self.sequence_count
    }
    fn sensitivity_round_count(&self) -> usize {
        1
    }
    fn begin_master(&mut self) -> io::Result<()> {
        self.input.prepare_for_master()?;
        let config = self.config;
        let mut main =
            NativeSearchWriter::new_with_query_mode(config, self.sequence_count, self.query_mode)?;
        main.set_database_letters(self.full_letters);
        let format = main.intermediate_format();
        let layout = dictionary_layout(config, &main);
        let auxiliary = crate::commands::query_output::QueryOutputs::new(
            &config.auxiliary_output,
            if self.dna.is_some() || config.translated_query_layout.is_some() {
                crate::basic::value::SequenceType::Nucleotide
            } else {
                crate::basic::value::SequenceType::AminoAcid
            },
            config.file_buffer_size,
        )?;
        self.output = Some(MppOutput {
            main,
            format,
            layout,
            auxiliary,
        });
        Ok(())
    }
    fn load_query(&mut self, ordinal: usize) -> io::Result<Option<MppQuery>> {
        let chunk = self.input.chunk(self.budget)?;
        let count = chunk.count();
        if count == 0 {
            return Ok(None);
        }
        if self.dna.is_none()
            && self.check_query_protein
            && ordinal == 0
            && chunk.translated.is_none()
        {
            crate::data::sequence_file::check_protein_input(
                chunk
                    .records
                    .iter()
                    .map(|record| record.sequence.as_slice()),
                self.config.tuning.ignore_warnings,
            )
            .map_err(io::Error::other)?;
        }
        let begin = self.query_begin;
        self.query_begin += count as u64;
        Ok(Some(MppQuery {
            chunk,
            ordinal,
            begin,
            master_oids: Vec::new(),
            output_oids: Vec::new(),
            aligned: vec![false; count],
            prepared: false,
            master_prepared: self.dna.is_none(),
            block_config: self.config.clone(),
            work: None,
        }))
    }
    fn query_is_selected(&self, ordinal: usize, _: &MppQuery) -> bool {
        let selected = self.config.tuning.multiprocessing.query_chunk;
        selected < 0 || selected as usize == ordinal
    }
    fn prepare_query(
        &mut self,
        query: &mut MppQuery,
        _: usize,
        previous_reference: usize,
    ) -> io::Result<()> {
        let config = self.config;
        query.master_oids = super::prepare_master_query_chunk(
            &mut query.chunk,
            config,
            u32::try_from(previous_reference).map_err(io::Error::other)?,
            self.min_length_ratio,
        )?;
        query.output_oids = query.master_oids.clone();
        if query.ordinal == 0 {
            let first = query.chunk.query(0)?;
            self.output
                .as_mut()
                .expect("master output initialized")
                .main
                .header_bytes(&self.matrix, first.title, first.source.len())?;
        }
        if self.dna.is_some() {
            if usize::from(config.tuning.setup.lin_stage1_target)
                + usize::from(config.tuning.setup.lin_stage1_query)
                + usize::from(config.tuning.setup.lin_stage1_combo)
                > 1
            {
                return Err(io::Error::other(
                    "Multiple linearization options are not allowed to be used together",
                ));
            }
            if config.target_indexed || config.query_indexed {
                return Err(io::Error::other("Unsupported reduction."));
            }
            if config.contiguous_seed_mode {
                return Err(io::Error::other("Unsupported contiguous seed."));
            }
        }
        query.work = Some(self.queues.query(query.ordinal)?);
        query.block_config.tuning.query_pre_sorted = true;
        query.block_config.tuning.search_min_length_ratio = Some(self.min_length_ratio);
        query.block_config.tuning.search_block.query = u32::try_from(query.ordinal)
            .map_err(|_| io::Error::other("Query block index exceeds uint32"))?;
        query.block_config.tuning.search_block.self_search =
            config.tuning.self_search || self.cluster.is_some_and(|options| options.self_search);
        query.block_config.scoring_matrix = Some(self.matrix.clone());
        if let Some(layout) = &mut query.block_config.translated_query_layout {
            layout.mate_indices = query.chunk.mates.clone();
        }
        if query.master_prepared {
            prepare_join_only_queries(
                &mut query.chunk,
                &query.block_config,
                &self.matrix,
                self.masker.as_deref(),
            );
        }
        Ok(())
    }
    fn query_has_pending_alignments(&self, query: &MppQuery) -> bool {
        query.chunk.count() != 0
    }
    fn begin_round(&mut self, query: &mut MppQuery, _: usize) -> io::Result<()> {
        if self.dna.is_none() {
            self.config.tuning.algorithm.begin_round()?;
            super::prepare_unsearched_query_chunk(
                &query.block_config,
                &mut query.chunk,
                &self.matrix,
                self.masker.as_deref(),
                false,
                self.references.scoring_letters(),
            )?;
        }
        Ok(())
    }
    fn load_reference(
        &mut self,
        query: &mut MppQuery,
        _: &mut (),
        _: usize,
    ) -> io::Result<Option<ClaimedReference>> {
        Ok(self
            .queues
            .claim_reference(query.work.as_ref().expect("query queues prepared"))?
            .map(|(record, chunk)| ClaimedReference { record, chunk }))
    }
    fn run_reference(
        &mut self,
        query: &mut MppQuery,
        _: &mut (),
        claimed: ClaimedReference,
    ) -> io::Result<crate::run::double_indexed::owned::ReferenceFlow> {
        let config = self.config;
        let mp = &config.tuning.multiprocessing;
        let references = &mut self.references;
        let queues = &self.queues;
        let partition = &self.partition;
        let dna = self.dna;
        let dna_statistics = self.dna_statistics;
        let dna_min_length_ratio = self.dna_min_length_ratio;
        let cluster = self.cluster;
        let cluster_owner = self.cluster_owner;
        let matrix = &self.matrix;
        let masker = &self.masker;
        let resources = self.output.as_ref().expect("master output initialized");
        let layout = resources.layout;
        let format = resources.format;
        let sequence_count = self.sequence_count;
        let unaligned_targets = &self.unaligned_targets;
        let query_block = query.ordinal;
        let master_prepared = query.master_prepared;
        let MppQuery {
            chunk,
            master_oids,
            output_oids,
            prepared,
            aligned,
            block_config,
            work,
            ..
        } = query;
        let work = work.as_ref().expect("query queues prepared");
        let ClaimedReference {
            record,
            chunk: reference,
        } = claimed;
        let reference_id = usize::try_from(reference.i).map_err(io::Error::other)?;
        block_config.tuning.search_block.reference = u32::try_from(reference_id)
            .map_err(|_| io::Error::other("Reference block index exceeds uint32"))?;
        self.source_reference_block = block_config.tuning.search_block.reference;
        let path = get_ref_block_tmpfile_name(&mp.directory, query_block, reference_id);
        let mut output = io::BufWriter::with_capacity(
            config.file_buffer_size,
            File::create(source_c_path(&path).as_ref())?,
        );
        let mut dictionary = ReferenceDictionary::new(sequence_count, layout.sequences);
        if !mp.self_compare || reference_id >= query_block {
            let count = u64::try_from(reference.n_seqs).map_err(io::Error::other)?;
            let mut records = Vec::new();
            let mut oids = Vec::new();
            for oid in reference.offset as u64..reference.offset as u64 + count {
                if let ReferenceSource::Database {
                    filter: Some(filter),
                    ..
                } = &references
                {
                    if !filter.get(oid) {
                        continue;
                    }
                }
                records.push(references.record(oid)?);
                oids.push(oid);
            }
            block_config.reference_taxon_rank_ids =
                if config.taxon_k != 0 && cluster_owner.is_none() {
                    Some(references.species_rank_ids(&oids)?.into())
                } else {
                    None
                };
            if dna.is_some()
                && (config.tuning.setup.lin_stage1_target || dna_min_length_ratio > 0.0)
                && !config.tuning.kmer_ranking
            {
                crate::dna::native::sort_records_by_length(&mut records, &mut oids)
                    .map_err(io::Error::other)?;
            }
            let mut capture = if let (Some(options), Some(statistics)) = (dna, dna_statistics) {
                let _ = crate::util::system::log_rss();
                let mut sequences = crate::data::sequence_set::LetterStringSet::new();
                for record in &records {
                    sequences.push_back(&record.sequence);
                }
                sequences.finish_reserve();
                let threads = if config.threads > 0 {
                    config.threads
                } else {
                    crate::util::system::hardware_concurrency() as i32
                };
                let engine = crate::dna::native::NativeDnaBlock::new_with_seed_shapes_bytes(
                    sequences,
                    statistics,
                    config.sensitivity,
                    options.clone(),
                    config.tuning.setup.minimizer_window,
                    config.min_query_len,
                    threads,
                    config.index_chunks.max(1) as i32,
                    config.max_evalue,
                    &config.tuning.setup.shape_mask,
                    config.tuning.setup.shapes,
                )
                .map_err(io::Error::other)?;
                let _ = crate::util::system::log_rss();
                let _ = crate::util::system::log_rss();
                let queries = std::mem::take(&mut chunk.records);
                let matches = queries
                    .iter()
                    .map(|query| engine.search(&query.sequence).map_err(io::Error::other))
                    .collect::<io::Result<Vec<_>>>()?;
                let query_self_aln_scores = queries
                    .iter()
                    .map(|query| {
                        matrix.bitscore(
                            query
                                .sequence
                                .iter()
                                .filter(|&&l| (l & crate::basic::value::LETTER_MASK) < 4)
                                .count() as f64
                                * f64::from(options.reward),
                        )
                    })
                    .collect();
                let target_self_aln_scores = records
                    .iter()
                    .map(|record| {
                        matrix.bitscore(
                            record
                                .sequence
                                .iter()
                                .filter(|&&l| (l & crate::basic::value::LETTER_MASK) < 4)
                                .count() as f64
                                * f64::from(options.reward),
                        )
                    })
                    .collect();
                crate::commands::blastp::NativeSearchCapture {
                    queries,
                    matches,
                    query_self_aln_scores,
                    target_self_aln_scores,
                    ..Default::default()
                }
            } else if let Some(options) = cluster {
                crate::commands::blastp::search_cluster_block(
                    &block_config,
                    records.clone(),
                    Some(std::mem::take(&mut chunk.records)),
                    references.scoring_letters(),
                    *prepared || master_prepared,
                    options,
                )?
            } else {
                crate::commands::blastp::search_block(
                    &block_config,
                    records.clone(),
                    std::mem::take(&mut chunk.records),
                    chunk.translated.take(),
                    references.scoring_letters(),
                    *prepared || master_prepared,
                    masker.as_deref(),
                    None,
                )?
            };
            chunk.records = std::mem::take(&mut capture.queries);
            chunk.translated = capture.translated.take();
            if let Some(tracker) = &unaligned_targets {
                tracker.mark_capture(&capture, &oids)?;
            }
            if !*prepared {
                if dna.is_none() && !capture.query_oids.is_empty() {
                    *output_oids = if master_oids.is_empty() {
                        capture.query_oids.clone()
                    } else {
                        capture
                            .query_oids
                            .iter()
                            .map(|&local| master_oids[local as usize])
                            .collect()
                    };
                }
                if !capture.query_self_aln_scores.is_empty() {
                    chunk.self_aln_scores = std::mem::take(&mut capture.query_self_aln_scores);
                }
            }
            for (query, matches) in capture.matches.iter().enumerate() {
                if matches.iter().any(|target| !target.hsps.is_empty()) {
                    aligned[query] = true;
                }
            }
            dictionary.write_block(&mut capture, &records, &oids, &format, &mut output)?;
            *prepared = true;
        } else {
            IntermediateRecord::finish_file(&mut output)?;
        }
        output.flush()?;
        drop(output);
        save_dictionary(
            &dictionary_file_name(&mp.directory, query_block, reference_id),
            &dictionary,
            layout,
            config.file_buffer_size,
        )?;
        queues.finish_reference(&work, &record, reference, partition.len())?;
        let _ = crate::util::system::log_rss();

        Ok(crate::run::double_indexed::owned::ReferenceFlow::Continue)
    }
    fn finish_round(&mut self, query: &mut MppQuery, _: (), _: usize) -> io::Result<()> {
        let matrix = &self.matrix;
        let masker = &self.masker;
        let references = &self.references;
        let dna = self.dna;

        if !query.prepared {
            if let Some(options) = dna {
                query.chunk.self_aln_scores = query
                    .chunk
                    .records
                    .iter()
                    .map(|query| {
                        matrix.bitscore(
                            query
                                .sequence
                                .iter()
                                .filter(|&&l| (l & crate::basic::value::LETTER_MASK) < 4)
                                .count() as f64
                                * f64::from(options.reward),
                        )
                    })
                    .collect();
            } else {
                if let Some(oids) = super::prepare_unsearched_query_chunk(
                    &query.block_config,
                    &mut query.chunk,
                    &matrix,
                    masker.as_deref(),
                    !query.master_prepared,
                    references.scoring_letters(),
                )? {
                    query.output_oids = oids;
                }
            }
        }

        Ok(())
    }
    fn log_after_reference_scan(&self) -> bool {
        false
    }
    fn finish_query(
        &mut self,
        mut query: MppQuery,
    ) -> io::Result<crate::run::double_indexed::owned::QueryCompletion> {
        let config = self.config;
        let mp = &config.tuning.multiprocessing;
        let sequence_count = self.sequence_count;
        let full_letters = self.full_letters;
        let query_mode = self.query_mode;
        let resources = self.output.as_mut().expect("master output initialized");
        let layout = resources.layout;
        let format = resources.format;
        let partition = &self.partition;
        let matrix = &self.matrix;
        let queues = &self.queues;
        let references = &mut self.references;
        let cluster_owner = self.cluster_owner;
        let auxiliary = &mut resources.auxiliary;
        let query_count = query.chunk.count();
        let query_block = query.ordinal;
        let query_begin = query.begin;
        let MppQuery {
            chunk,
            output_oids,
            aligned,
            work,
            ..
        } = &mut query;
        let work = work.as_ref().expect("query queues prepared");
        if let Some(token) = queues.claim_join(&work)? {
            self.source_reference_block = u32::try_from(partition.len())
                .map_err(|_| io::Error::other("Reference block count exceeds uint32"))?;
            let mut dictionary = ReferenceDictionary::new(sequence_count, layout.sequences);
            let mut mappings = Vec::new();
            let mut files = Vec::new();
            for reference in 0..partition.len() {
                let path = dictionary_file_name(&mp.directory, query_block, reference);
                let entries = load_dictionary(&path, layout, config.file_buffer_size)?;
                std::fs::remove_file(source_c_path(&(path)).as_ref())?;
                let mut ids = Vec::new();
                for entry in entries {
                    let id = dictionary.entries.len();
                    ids.push(id);
                    dictionary.by_oid.insert(entry.oid, id);
                    dictionary.entries.push(ReferenceInfo {
                        oid: entry.oid,
                        title: entry.title,
                        len: entry.len as usize,
                        sequence: layout.sequences.then(|| entry.seq.into()),
                        self_aln_score: layout.self_scores.then_some(entry.self_aln_score),
                        taxids: Vec::new(),
                    });
                }
                mappings.push(ids);
                files.push(io::BufReader::with_capacity(
                    config.file_buffer_size,
                    File::open(
                        source_c_path(
                            &(get_ref_block_tmpfile_name(&mp.directory, query_block, reference)),
                        )
                        .as_ref(),
                    )?,
                ));
            }
            references.load_matched_taxonomy(&mut dictionary)?;
            let mut output_config = config.clone();
            let mut output_path = config
                .output
                .as_deref()
                .unwrap_or(std::path::Path::new(""))
                .as_os_str()
                .to_os_string();
            output_path.push(format!("_{query_block}"));
            output_config.output = Some(output_path.into());
            let mut output = NativeSearchWriter::new_with_query_mode(
                &output_config,
                sequence_count,
                query_mode,
            )?;
            output.set_database_letters(full_letters);
            for entry in &dictionary.entries {
                output.reserve_archive_subject_bytes(entry.oid, &entry.title, entry.len);
            }
            let mut fetcher = crate::output::join_blocks::StreamJoinFetcher::from_readers(files)?;
            let mut exhausted = false;
            let mut last_query = None;
            let (renderer, mut sink) = output.joined_parts(matrix);
            let render = |query_id: usize, matches, aligned| {
                let query = chunk.query(query_id)?;
                renderer.render_query_bytes(
                    query_id as u32,
                    query.oid.unwrap_or(
                        query_begin
                            + output_oids
                                .get(query_id)
                                .copied()
                                .unwrap_or(query_id as u64),
                    ),
                    query.title,
                    query.source,
                    &query.frames,
                    query.mate,
                    query.quality,
                    query.self_aln_score,
                    matches,
                    &dictionary,
                    matrix,
                    aligned,
                )
            };
            crate::output::join_blocks::run_join_workers(
                config.threads.max(1) as usize,
                || {
                    let payload = crate::output::join_blocks::fetch_join_payload(
                        &mut fetcher,
                        &mut exhausted,
                    )?;
                    if let Some(payload) = &payload {
                        last_query = Some(payload.query_id);
                    }
                    Ok(payload)
                },
                |mut payload| {
                    for (bytes, mapping) in payload.buffers.iter_mut().zip(&mappings) {
                        remap_payload(bytes, mapping, &format)?;
                    }
                    if config.taxon_k != 0 {
                        if let Some(owner) = cluster_owner {
                            if let Some(bytes) =
                                payload.buffers.iter().find(|bytes| !bytes.is_empty())
                            {
                                let mut input = std::io::Cursor::new(bytes);
                                let record = IntermediateRecord::read(&mut input, &format)?;
                                return Err(owner.taxonomy_access(
                                    dictionary.original_oid(record.target_dict_id as usize),
                                ));
                            }
                        }
                    }
                    let mut output = crate::commands::blastp::NativeRenderedQuery::join_task();
                    if config.report_unaligned != 0
                        && format.code != crate::output::format::FormatCode::Daa
                    {
                        for id in payload.unaligned_from..payload.query_id {
                            output.append(render(id as usize, Vec::new(), false)?);
                        }
                    }
                    let query = chunk.query(payload.query_id as usize)?;
                    let frames: Vec<_> = query.frames.iter().map(|frame| frame.to_vec()).collect();
                    let matches = super::merge_query(
                        &payload.buffers,
                        &frames,
                        query.source.len(),
                        &dictionary,
                        matrix,
                        config,
                        format,
                    )?;
                    output.append(render(payload.query_id as usize, matches, true)?);
                    Ok(output)
                },
                |output| sink.write(output),
            )?;
            fetcher.finish();
            if config.report_unaligned != 0 && format.code != crate::output::format::FormatCode::Daa
            {
                for id in last_query.map_or(0, |id| id + 1)..query_count as u32 {
                    sink.write(render(id as usize, Vec::new(), false)?)?;
                }
            }
            output.finish_multiprocessing_chunk()?;
            queues.finish_join(&work, &token)?;
            for reference in 0..partition.len() {
                std::fs::remove_file(
                    source_c_path(
                        &(get_ref_block_tmpfile_name(&mp.directory, query_block, reference)),
                    )
                    .as_ref(),
                )?;
            }
        }
        for pass in [false, true] {
            for query in 0..query_count {
                if aligned[query] == pass {
                    let query = chunk.query(query)?;
                    auxiliary.write_query_bytes(
                        pass,
                        query.title,
                        query.source,
                        (!query.quality.is_empty()).then_some(query.quality),
                    )?;
                }
            }
        }

        Ok(crate::run::double_indexed::owned::QueryCompletion {
            reference_counter: self.source_reference_block as usize,
        })
    }
    fn stop_exists(&self) -> bool {
        stop_exists()
    }
    fn finish(mut self) -> io::Result<ReferenceSource> {
        crate::data::accession_spool::source_scope(self.references.accession_spool(), || {
            let output = self.output.take().expect("master output initialized");
            output.auxiliary.finish()?;
            output.main.finish(&self.matrix)?;
            if let Some(tracker) = self.unaligned_targets {
                tracker.finish_spooled()?;
            }
            Ok(self.references)
        })
    }
}
