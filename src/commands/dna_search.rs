//! Nucleotide orchestration inside the shared block-reader/merger owner.
use super::*;
use crate::basic::value::SequenceType;
use crate::data::fasta::ByteFastaRecord;
use crate::dna::build_score::Blastn_Score;
use crate::dna::native::NativeDnaBlock;

pub(crate) fn run(options: &crate::commands::blastn::BlastnConfig, budget: u64) -> io::Result<()> {
    crate::basic::statistics::with_search_statistics_prepared(
        !options.search.tuning.multiprocessing.enabled
            || (!options.search.tuning.multiprocessing.initialize
                && !options.search.tuning.multiprocessing.recover),
        || crate::commands::source_config::begin_search_temp(&options.search),
        || {
            if !options.search.tuning.self_search && options.search.query_files.len() > 1 {
                return Err(io::Error::other("--query/-q has more than one argument."));
            }
            let normalized =
                crate::commands::blastp::normalize_output_config_for_dna(&options.search)?;
            let config = normalized.as_ref();
            let errors = crate::util::io::IoErrorOwner::default();
            let search_setup =
                errors.finish(crate::run::config::Config::from_input_with_io_owner(
                    &crate::run::config::ConfigInput {
                        blastn: true,
                        gap_open: config.gap_open,
                        gap_extend: config.gap_extend,
                        iterate: config.iterate.clone(),
                        self_search: config.tuning.self_search,
                        sensitivity: config.sensitivity,
                        multiprocessing: config.tuning.multiprocessing.enabled,
                        global_ranking_targets: config.global_ranking_targets,
                        target_indexed: config.target_indexed,
                        frame_shift: config.frame_shift,
                        algo: if config.contiguous_seed_mode {
                            crate::run::config::Algo::CtgSeed
                        } else if config.query_indexed {
                            crate::run::config::Algo::QueryIndexed
                        } else {
                            crate::run::config::Algo::Auto
                        },
                        minimizer_window: config.tuning.setup.minimizer_window.unwrap_or(0),
                        freq_masking: config.tuning.setup.freq_masking,
                        seed_cut: config.tuning.setup.seed_cut,
                        freq_sd: config.tuning.setup.freq_sd.unwrap_or(0.0),
                        query_cover: config.query_cover,
                        subject_cover: config.subject_cover,
                        min_length_ratio: config.tuning.min_length_ratio,
                        lin_stage1_query: config.tuning.setup.lin_stage1_query,
                        lin_stage1_target: config.tuning.setup.lin_stage1_target,

                        ..Default::default()
                    },
                    &errors,
                ))?;
            if let Some(message) = &search_setup.iteration_message {
                search_message!("{message}");
            }
            let timing = crate::dna::native::NativeDnaTiming::new();
            let _timing_scope = timing.enter();
            // Source --swipe leaves dna_ref_index null, while ordinary DNA extension
            // dereferences it. Reject that UB boundary instead of silently running a
            // seeded search. Global ranking bypasses align_queries and stays empty.
            if config.swipe_all && config.global_ranking_targets == 0 {
                return Err(io::Error::other(
                    "Exhaustive nucleotide search (--swipe) is not supported.",
                ));
            }
            let references =
                ReferenceSource::open_with_sequence_type(config, SequenceType::Nucleotide)?;
            crate::data::accession_spool::source_scope(references.accession_spool(), || {
                let (sequence_count, full_letters) = references.totals();
                let selected_sequence_count = references.selected_sequence_count();
                let metadata_db_letters = if config.tuning.db_size != 0 {
                    config.tuning.db_size
                } else {
                    references.scoring_letters()
                };
                let statistics = Blastn_Score::new(
                    options.dna.reward,
                    options.dna.penalty,
                    search_setup.gap_open,
                    search_setup.gap_extend,
                    full_letters,
                    i64::try_from(sequence_count).map_err(io::Error::other)?,
                )
                .map_err(io::Error::other)?;
                // The static Source global is zero-initialized, not an automatic object.
                // Its metadata is separate from real Blastn_Score HSP statistics.
                let mut matrix =
                    ScoreMatrix::source_blastn_metadata(statistics).map_err(io::Error::other)?;
                // Source's separate global matrix uses override/filtered letters, while
                // Blastn_Score's HSP model retains cfg.db_letters from the full database.
                matrix.set_db_letters(metadata_db_letters);
                if config.tuning.multiprocessing.enabled {
                    if config.tuning.self_search {
                        let input = QueryInput::Protein(ProteinQueries::Database(Box::new(
                            references.query_owner(config)?,
                        )));
                        return super::multiprocessing::run_dna(
                            config,
                            references,
                            input,
                            budget,
                            &options.dna,
                            search_setup.min_length_ratio,
                        );
                    }
                    let path = config.query_files.first().cloned().unwrap_or_default();
                    let keep_quality = config
                        .outfmt
                        .iter()
                        .any(|field| matches!(field.as_slice(), b"qqual" | b"full_qqual"))
                        || config.auxiliary_output.needs_quality()
                        || config.tuning.store_query_quality;
                    let stream = RecordStream::open(
                        path.into(),
                        SequenceType::Nucleotide,
                        keep_quality,
                        config.file_buffer_size,
                    );
                    return super::multiprocessing::run_dna(
                        config,
                        references,
                        QueryInput::Protein(ProteinQueries::Stream(stream)),
                        budget,
                        &options.dna,
                        search_setup.min_length_ratio,
                    );
                }
                crate::run::double_indexed::owned::run(DnaBackend {
                    options,
                    config,
                    budget,
                    references,
                    sequence_count,
                    full_letters,
                    selected_sequence_count,
                    statistics,
                    matrix,
                    setup: search_setup,
                    writer: None,
                    auxiliary: None,
                    tracker: None,
                    input: None,
                    query_begin: 0,
                    skipped_source_count: None,
                })
            })
        },
    )
}

use crate::run::double_indexed::owned::{DoubleIndexedBackend, QueryCompletion, ReferenceFlow};

/// Actual nucleotide source-query owner across all rounds and references.
struct DnaQuery {
    records: Vec<ByteFastaRecord>,
    oids: Vec<u64>,
    oid_end: u64,
    ordinal: usize,
    source_begin: u64,
    dictionary: Option<ReferenceDictionary>,
    files: Vec<IntermediateFile>,
    direct: Option<Vec<Vec<Match>>>,
    aligned: Vec<bool>,
    reference_counter: usize,
}
struct DnaRound {
    sensitivity: crate::run::config::Round,
    linear_target: bool,
    previously_aligned: usize,
    ranking: Option<Vec<crate::align::global_ranking::Hit>>,
    reference_counter: usize,
}
/// Concrete ordinary DNA services. MPP returns through its actual lazy queue
/// owner before this backend is constructed; no jobs are fabricated here.
struct DnaBackend<'a> {
    options: &'a crate::commands::blastn::BlastnConfig,
    config: &'a BlastpConfig,
    budget: u64,
    references: ReferenceSource,
    sequence_count: u64,
    full_letters: u64,
    selected_sequence_count: u64,
    statistics: Blastn_Score,
    matrix: ScoreMatrix,
    setup: crate::run::config::Config,
    writer: Option<NativeSearchWriter>,
    auxiliary: Option<crate::commands::query_output::QueryOutputs>,
    tracker: Option<UnalignedTargets>,
    input: Option<QueryInput<'a>>,
    query_begin: u64,
    // A deselected block has no finish_query callback; advance its real
    // ordinal before the next read, preserving Source's skip lifecycle.
    skipped_source_count: Option<usize>,
}
impl DoubleIndexedBackend for DnaBackend<'_> {
    type Query = DnaQuery;
    type Round = DnaRound;
    type Reference = ReferenceChunk;
    type Output = ();
    fn capture_source_error(&self, error: io::Error) -> io::Error {
        let errors = self.references.source_errors();
        errors.capture_io::<()>(Err(error)).unwrap_err()
    }
    fn database_sequence_count(&self) -> u64 {
        self.sequence_count
    }
    fn sensitivity_round_count(&self) -> usize {
        self.setup.sensitivity.len()
    }
    fn begin_master(&mut self) -> io::Result<()> {
        let config = self.config;
        let sequence_count = self.sequence_count;
        let full_letters = self.full_letters;
        let references = &mut self.references;
        let path = config.query_files.first().cloned().unwrap_or_default();
        let keep_quality = config
            .outfmt
            .iter()
            .any(|field| matches!(field.as_slice(), b"qqual" | b"full_qqual"))
            || config.auxiliary_output.needs_quality()
            || config.tuning.store_query_quality;
        let mut input = if config.tuning.self_search {
            QueryInput::Protein(ProteinQueries::Database(Box::new(
                references.query_owner(config)?,
            )))
        } else {
            QueryInput::Protein(ProteinQueries::Stream(RecordStream::open(
                path.into(),
                SequenceType::Nucleotide,
                keep_quality,
                config.file_buffer_size,
            )))
        };

        input.prepare_for_master()?;
        let mut writer = NativeSearchWriter::new_with_query_mode(config, sequence_count, 4)?;
        writer.set_database_letters(full_letters);
        let auxiliary = crate::commands::query_output::QueryOutputs::new(
            &config.auxiliary_output,
            SequenceType::Nucleotide,
            config.file_buffer_size,
        )?;
        let tracker = UnalignedTargets::new(config, references)?;

        self.writer = Some(writer);
        self.auxiliary = Some(auxiliary);
        self.tracker = tracker;
        self.input = Some(input);
        Ok(())
    }
    fn load_query(&mut self, ordinal: usize) -> io::Result<Option<DnaQuery>> {
        if let Some(count) = self.skipped_source_count.take() {
            self.query_begin = self
                .query_begin
                .checked_add(count as u64)
                .ok_or_else(|| io::Error::other("Query ordinal overflow"))?;
        }
        let chunk = self
            .input
            .as_mut()
            .ok_or_else(|| io::Error::other("DNA query input is not initialized"))?
            .chunk(self.budget)?;
        let oid_end = chunk
            .protein_oids
            .iter()
            .copied()
            .max()
            .map_or(0, |oid| oid + 1);
        let records = chunk.records;
        if records.is_empty() {
            return Ok(None);
        }
        let oids = if self.config.tuning.self_search {
            chunk.protein_oids
        } else {
            (0..records.len() as u64).collect()
        };
        let aligned = vec![false; records.len()];
        self.skipped_source_count = Some(records.len());
        Ok(Some(DnaQuery {
            records,
            oids,
            oid_end,
            ordinal,
            source_begin: self.query_begin,
            dictionary: None,
            files: Vec::new(),
            direct: None,
            aligned,
            reference_counter: 0,
        }))
    }
    fn query_is_selected(&self, ordinal: usize, _: &DnaQuery) -> bool {
        self.config.tuning.multiprocessing.query_chunk < 0
            || self.config.tuning.multiprocessing.query_chunk as usize == ordinal
    }
    fn prepare_query(&mut self, query: &mut DnaQuery, _: usize, _: usize) -> io::Result<()> {
        self.skipped_source_count = None;
        let config = self.config;
        let search_setup = &self.setup;
        let queries = &mut query.records;
        let query_oids = &mut query.oids;
        let query_begin = query.source_begin;
        let matrix = &self.matrix;
        let sequence_count = self.sequence_count;
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| io::Error::other("DNA output is not initialized"))?;
        let keep_ids = config.tuning.self_search
            || search_setup.min_length_ratio != 0.0
            || config.global_ranking_targets != 0
            || config.tuning.setup.lin_stage1_combo;
        if search_setup.min_length_ratio > 0.0
            || (config.tuning.setup.lin_stage1_query && !keep_ids && !config.tuning.kmer_ranking)
        {
            crate::dna::native::sort_records_by_length(queries, query_oids)
                .map_err(io::Error::other)?;
        }
        if query_begin == 0 {
            writer.header_bytes(&matrix, &queries[0].id, queries[0].sequence.len())?;
        }
        validate_seed_encoding(config)?;
        let dictionary =
            ReferenceDictionary::new(sequence_count, writer.needs_reference_sequences());

        query.dictionary = Some(dictionary);
        Ok(())
    }
    fn query_has_pending_alignments(&self, query: &DnaQuery) -> bool {
        query.aligned.iter().any(|&aligned| !aligned)
    }
    fn begin_round(&mut self, query: &mut DnaQuery, ordinal: usize) -> io::Result<DnaRound> {
        let config = self.config;
        let round = &self.setup.sensitivity[ordinal];
        let queries = &query.records;
        let linear_target = config.tuning.setup.lin_stage1_target || round.linearize;
        if usize::from(linear_target)
            + usize::from(config.tuning.setup.lin_stage1_query)
            + usize::from(config.tuning.setup.lin_stage1_combo)
            > 1
        {
            return Err(io::Error::other(
                "Multiple linearization options are not allowed to be used together",
            ));
        }
        // Source allocates this table for blastn too, but its DNA branch
        // bypasses every search_shape/update_table call. Retain that real
        // empty-ranking lifecycle rather than invent a nucleotide ranker.
        let ranking = if config.global_ranking_targets != 0 {
            let count = usize::try_from(config.global_ranking_targets)
                .ok()
                .and_then(|n| queries.len().checked_mul(n))
                .ok_or_else(|| {
                    io::Error::other("cannot create std::vector larger than max_size()")
                })?;
            let mut table = Vec::new();
            table.try_reserve_exact(count).map_err(|_| {
                io::Error::other("cannot create std::vector larger than max_size()")
            })?;
            table.resize(count, crate::align::global_ranking::Hit::default());
            Some(table)
        } else {
            None
        };

        self.references.rewind();
        if config.tuning.self_search && !config.tuning.setup.lin_stage1_query {
            self.references.seek_oid(query.oid_end);
        }
        Ok(DnaRound {
            sensitivity: *round,
            linear_target,
            previously_aligned: query.aligned.iter().filter(|&&aligned| aligned).count(),
            ranking,
            reference_counter: 0,
        })
    }
    fn load_reference(
        &mut self,
        query: &mut DnaQuery,
        round: &mut DnaRound,
        cursor: usize,
    ) -> io::Result<Option<ReferenceChunk>> {
        let config = self.config;
        let query_block = query.ordinal;
        let query_oid_end = query.oid_end;
        let queries = &query.records;
        let query_oids = &query.oids;
        let references = &mut self.references;
        let budget = self.budget;
        let alias = config.tuning.self_search
            && ((!config.tuning.setup.lin_stage1_query && cursor == 0)
                || (config.tuning.setup.lin_stage1_query && cursor == query_block));
        let reference = if alias {
            if config.tuning.setup.lin_stage1_query {
                references.seek_oid(query_oid_end);
            }
            ReferenceChunk {
                records: queries.clone(),
                oids: query_oids.clone(),
            }
        } else {
            references.chunk(budget)?
        };
        round.reference_counter = cursor
            .checked_add(1)
            .ok_or_else(|| io::Error::other("Reference block index overflow"))?;
        if reference.records.is_empty() {
            return Ok(None);
        }

        Ok(Some(reference))
    }
    fn run_reference(
        &mut self,
        query: &mut DnaQuery,
        round_state: &mut DnaRound,
        reference: ReferenceChunk,
    ) -> io::Result<ReferenceFlow> {
        let ReferenceChunk {
            mut records,
            mut oids,
        } = reference;
        let config = self.config;
        let options = self.options;
        let search_setup = &self.setup;
        let round = &round_state.sensitivity;
        let linear_target = round_state.linear_target;
        let ranking = &round_state.ranking;
        let queries = &query.records;
        let aligned = &mut query.aligned;
        let dictionary = query
            .dictionary
            .as_mut()
            .ok_or_else(|| io::Error::other("DNA dictionary is not initialized"))?;
        let files = &mut query.files;
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| io::Error::other("DNA output is not initialized"))?;
        let format = writer.intermediate_format();
        let tracker = &self.tracker;
        let references = &mut self.references;
        let matrix = &self.matrix;
        let statistics = self.statistics;
        let selected_sequence_count = self.selected_sequence_count;
        let iterated = self.setup.sensitivity.len() > 1;
        let threads = if config.threads > 0 {
            config.threads
        } else {
            crate::util::system::hardware_concurrency() as i32
        };
        let chunks = if config.index_chunks > 0 {
            config.index_chunks as i32
        } else {
            1
        };
        let minimizer_window = config.tuning.setup.minimizer_window;
        let _ = crate::util::system::log_rss();
        if (linear_target || search_setup.min_length_ratio > 0.0) && !config.tuning.kmer_ranking {
            crate::dna::native::sort_records_by_length(&mut records, &mut oids)
                .map_err(io::Error::other)?;
        }
        let engine = if config.swipe_all && ranking.is_some() {
            None
        } else {
            let mut target_set = crate::data::sequence_set::LetterStringSet::new();
            for record in &records {
                target_set.push_back(&record.sequence);
            }
            target_set.finish_reserve();
            let engine = NativeDnaBlock::new_with_seed_shapes_bytes(
                target_set,
                statistics,
                round.sensitivity,
                options.dna.clone(),
                minimizer_window,
                config.min_query_len,
                threads,
                chunks,
                config.max_evalue,
                &config.tuning.setup.shape_mask,
                config.tuning.setup.shapes,
            )
            .map_err(io::Error::other)?;
            // Source samples after the DNA reference index is
            // built and after its query-buffer cleanup stage.
            let _ = crate::util::system::log_rss();
            let _ = crate::util::system::log_rss();
            Some(engine)
        };
        // Global-ranking Source skips align_queries completely. For
        // ordinary DNA, HitIterator bypasses query_skip in later rounds.
        let matches = if ranking.is_some() {
            (0..queries.len()).map(|_| Vec::new()).collect()
        } else {
            queries
                .iter()
                .map(|query| {
                    engine
                        .as_ref()
                        .expect("ordinary DNA reference index")
                        .search(&query.sequence)
                        .map_err(io::Error::other)
                })
                .collect::<io::Result<Vec<_>>>()?
        };
        let target_self_aln_scores = records
            .iter()
            .map(|record| {
                matrix.bitscore(
                    record
                        .sequence
                        .iter()
                        .filter(|&&l| (l & crate::basic::value::LETTER_MASK) < 4)
                        .count() as f64
                        * f64::from(options.dna.reward),
                )
            })
            .collect();
        let mut capture = NativeSearchCapture {
            matches,
            target_self_aln_scores,
            ..Default::default()
        };
        if let Some(tracker) = &tracker {
            tracker.mark_capture(&capture, &oids)?;
        }
        // Source tracks the aligned union for stopping/auxiliary
        // output, while DNA extension itself ignores query_skip.
        if !iterated && records.len() as u64 == selected_sequence_count {
            for (id, matches) in capture.matches.iter().enumerate() {
                aligned[id] |= !matches.is_empty();
            }
            dictionary.remap_direct(&mut capture, &records, &oids);
            references.load_matched_taxonomy(dictionary)?;
            query.direct = Some(capture.matches);
            return Ok(ReferenceFlow::FinishRound);
        }
        let (file, mut output) = IntermediateFile::new(config)?;
        for (id, matches) in capture.matches.iter().enumerate() {
            aligned[id] |= !matches.is_empty();
        }
        dictionary.write_block(&mut capture, &records, &oids, &format, &mut output)?;
        references.load_matched_taxonomy(dictionary)?;
        for subject in &dictionary.entries {
            writer.reserve_archive_subject_bytes(subject.oid, &subject.title, subject.len);
        }
        std::io::Write::flush(&mut output)?;
        drop(output);
        files.push(file);
        Ok(ReferenceFlow::Continue)
    }
    fn finish_round(&mut self, query: &mut DnaQuery, round: DnaRound, _: usize) -> io::Result<()> {
        let ranking = round.ranking;
        let sequence_count = self.sequence_count;
        let aligned = &query.aligned;
        let queries = &query.records;
        let previously_aligned = round.previously_aligned;
        let iterated = self.setup.sensitivity.len() > 1;
        if let Some(table) = ranking {
            // The actual global-ranking filter has no selected OIDs;
            // therefore its target-load/DP worklists are empty.
            let filter = crate::align::global_ranking::db_filter(
                &table,
                usize::try_from(sequence_count).map_err(io::Error::other)?,
            );
            debug_assert_eq!(filter.one_count(), 0);
        }
        if iterated {
            let total = aligned.iter().filter(|&&value| value).count();
            search_message!(
                "Aligned {}/{} queries in this iteration, {}/{} total.",
                total - previously_aligned,
                queries.len(),
                total,
                queries.len()
            );
        }
        query.reference_counter = round.reference_counter;
        Ok(())
    }
    fn finish_query(&mut self, query: DnaQuery) -> io::Result<QueryCompletion> {
        let DnaQuery {
            records: queries,
            oids: query_oids,
            source_begin: query_begin,
            dictionary,
            files,
            direct,
            aligned,
            reference_counter,
            ..
        } = query;
        let dictionary =
            dictionary.ok_or_else(|| io::Error::other("DNA dictionary is not initialized"))?;
        let config = self.config;
        let options = self.options;
        let matrix = &self.matrix;
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| io::Error::other("DNA output is not initialized"))?;
        let auxiliary = self
            .auxiliary
            .as_mut()
            .ok_or_else(|| io::Error::other("DNA auxiliary output is not initialized"))?;
        let format = writer.intermediate_format();
        let joined_output = direct.is_none() && !files.is_empty();
        let mut fetcher = if direct.is_none() && !files.is_empty() {
            Some(crate::output::join_blocks::StreamJoinFetcher::from_readers(
                files
                    .iter()
                    .map(IntermediateFile::reader)
                    .collect::<io::Result<Vec<_>>>()?,
            )?)
        } else {
            None
        };
        if joined_output {
            let mut fetcher = fetcher.take().unwrap();
            let mut exhausted = false;
            let mut trailing_from = 0usize;
            let (renderer, mut sink) = writer.joined_parts(matrix);
            let render = |id: usize, matches: Vec<crate::align::hsp::Match>, aligned: bool| {
                let query = &queries[id];
                // Only context zero owns untranslated DNA sequence data.
                let frames = vec![
                    query.sequence.clone(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                ];
                let frame_refs = frames.iter().map(Vec::as_slice).collect::<Vec<_>>();
                let self_score = query
                    .sequence
                    .iter()
                    .filter(|&&l| (l & crate::basic::value::LETTER_MASK) < 4)
                    .count() as f64
                    * f64::from(options.dna.reward);
                renderer.render_query_bytes(
                    id as u32,
                    if config.tuning.self_search {
                        query_oids[id]
                    } else {
                        query_begin + query_oids[id]
                    },
                    &query.id,
                    &query.sequence,
                    &frame_refs,
                    None,
                    query.quality.as_deref().unwrap_or(&[]),
                    Some(matrix.bitscore(self_score)),
                    matches,
                    &dictionary,
                    matrix,
                    aligned,
                )
            };
            let report_unaligned = config.report_unaligned != 0
                && !matches!(
                    config.outfmt.first().map(Vec::as_slice),
                    Some(b"100" | b"daab")
                );
            crate::output::join_blocks::run_join_workers(
                config.threads as usize,
                || {
                    let payload = crate::output::join_blocks::fetch_join_payload(
                        &mut fetcher,
                        &mut exhausted,
                    )?;
                    if let Some(payload) = &payload {
                        trailing_from = payload.query_id as usize + 1;
                    }
                    Ok(payload)
                },
                |payload| {
                    let id = payload.query_id as usize;
                    let query = &queries[id];
                    let frames = vec![
                        query.sequence.clone(),
                        Vec::new(),
                        Vec::new(),
                        Vec::new(),
                        Vec::new(),
                        Vec::new(),
                    ];
                    let matches = merge_query(
                        &payload.buffers,
                        &frames,
                        query.sequence.len(),
                        &dictionary,
                        matrix,
                        config,
                        format,
                    )?;
                    let mut task = crate::commands::blastp::NativeRenderedQuery::join_task();
                    if report_unaligned {
                        for missing in payload.unaligned_from as usize..id {
                            task.append(render(missing, Vec::new(), false)?);
                        }
                    }
                    task.append(render(id, matches, true)?);
                    Ok(task)
                },
                |rendered| sink.write(rendered),
            )?;
            fetcher.finish();
            if report_unaligned {
                for missing in trailing_from..queries.len() {
                    sink.write(render(missing, Vec::new(), false)?)?;
                }
            }
        } else {
            let mut direct = direct.map(Vec::into_iter);
            for (id, query) in queries.iter().enumerate() {
                // Source TranslatedSequence(Sequence) initializes frame0 and
                // default-constructs five empty (length0,null data) Sequence views.
                let frames = vec![
                    query.sequence.clone(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                ];
                let matches = if let Some(direct) = direct.as_mut() {
                    direct.next().unwrap_or_default()
                } else {
                    Vec::new()
                };
                let frame_refs = frames.iter().map(Vec::as_slice).collect::<Vec<_>>();
                let self_score = query
                    .sequence
                    .iter()
                    .filter(|&&l| (l & crate::basic::value::LETTER_MASK) < 4)
                    .count() as f64
                    * f64::from(options.dna.reward);
                writer.write_query_bytes(
                    id as u32,
                    if config.tuning.self_search {
                        query_oids[id]
                    } else {
                        query_begin + query_oids[id]
                    },
                    &query.id,
                    &query.sequence,
                    &frame_refs,
                    None,
                    query.quality.as_deref().unwrap_or(&[]),
                    Some(matrix.bitscore(self_score)),
                    matches,
                    &dictionary,
                    &matrix,
                )?;
            }
        }
        for pass in [false, true] {
            for (id, query) in queries.iter().enumerate() {
                if aligned[id] == pass {
                    auxiliary.write_query_bytes(
                        pass,
                        &query.id,
                        &query.sequence,
                        query.quality.as_deref(),
                    )?;
                }
            }
        }
        self.query_begin = query_begin
            .checked_add(queries.len() as u64)
            .ok_or_else(|| io::Error::other("Query ordinal overflow"))?;

        Ok(QueryCompletion { reference_counter })
    }
    fn stop_exists(&self) -> bool {
        crate::run::multiprocessing::stop_exists()
    }
    fn finish(self) -> io::Result<()> {
        crate::data::accession_spool::source_scope(self.references.accession_spool(), || {
            self.writer
                .ok_or_else(|| io::Error::other("DNA output is not initialized"))?
                .finish(&self.matrix)?;
            self.auxiliary
                .ok_or_else(|| io::Error::other("DNA auxiliary output is not initialized"))?
                .finish()?;
            if let Some(tracker) = self.tracker {
                tracker.finish_spooled()?;
            }
            Ok(())
        })
    }
}

/// HASHED dispatch in Source enum_seeds.h supports only four-bit reductions;
/// DNA has two bits. CONTIGUOUS likewise supports only 4-bit shapes5/6/7.
pub(super) fn validate_seed_encoding(config: &BlastpConfig) -> io::Result<()> {
    if config.swipe_all {
        return Ok(());
    }
    if config.target_indexed || config.query_indexed {
        return Err(io::Error::other("Unsupported reduction."));
    }
    if config.contiguous_seed_mode {
        return Err(io::Error::other("Unsupported contiguous seed."));
    }
    Ok(())
}
