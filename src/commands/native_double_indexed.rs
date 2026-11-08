//! Concrete owner services for the translated double-indexed orchestrator.
use super::*;
use crate::run::config::Round;
use crate::run::double_indexed::owned::{
    self, DoubleIndexedBackend, QueryCompletion, ReferenceFlow,
};
use std::sync::Arc;

struct NativeBackend<'a> {
    config: BlastpConfig,
    references: ReferenceSource,
    input: QueryInput<'a>,
    budget: u64,
    rounds: Vec<Round>,
    iterated: bool,
    min_length_ratio: f64,
    sequence_count: u64,
    full_database_letters: u64,
    selected_sequence_count: u64,
    letters: u64,
    matrix: ScoreMatrix,
    prepared_matrix: Arc<ScoreMatrix>,
    masker: Option<Arc<crate::masking::tantan::TantanMasker>>,
    writer: Option<NativeSearchWriter>,
    format: Option<OutputFormat>,
    auxiliary: Option<crate::commands::query_output::QueryOutputs>,
    unaligned_targets: Option<UnalignedTargets>,
    query_begin: u64,
    check_query_protein: bool,
}

struct QueryOwner {
    chunk: QueryChunk,
    query_count: usize,
    query_oid_end: u64,
    query_block: usize,
    query_begin: u64,
    block_config: BlastpConfig,
    files: Vec<IntermediateFile>,
    dictionary: ReferenceDictionary,
    prepared: bool,
    query_prepared: bool,
    master_oids: Vec<u64>,
    output_oids: Vec<u64>,
    direct: Option<(Vec<Vec<Match>>, Vec<bool>)>,
    aligned: Vec<bool>,
    source_reference_block: u32,
}
struct RoundOwner {
    ranking: Option<crate::commands::blastp::NativeRankingState>,
    previously_aligned: usize,
    first_reference_chunk: bool,
}

pub(super) fn run(
    config: &BlastpConfig,
    mut references: ReferenceSource,
    input: QueryInput<'_>,
    budget: u64,
    writer: Option<NativeSearchWriter>,
    rounds: Vec<Round>,
    min_length_ratio: f64,
) -> io::Result<ReferenceSource> {
    crate::data::accession_spool::source_scope(references.accession_spool(), || {
        let (sequence_count, full_database_letters) = references.totals();
        let selected_sequence_count = references.selected_sequence_count();
        let unaligned_targets = UnalignedTargets::new(config, &mut references)?;
        let letters = references.scoring_letters();
        let matrix = super::super::blastp::native_score_matrix(config, letters)?;
        let prepared_matrix = Arc::new(matrix.clone());
        let masker = (config.masking == crate::masking::MaskingMode::Tantan)
            .then(|| super::super::blastp::native_tantan_masker(config, &matrix));
        let check_query_protein = input.checks_protein_input();
        owned::run(NativeBackend {
            config: config.clone(),
            references,
            input,
            budget,
            iterated: rounds.len() > 1,
            rounds,
            min_length_ratio,
            sequence_count,
            full_database_letters,
            selected_sequence_count,
            letters,
            matrix,
            prepared_matrix,
            masker,
            writer,
            format: None,
            auxiliary: None,
            unaligned_targets,
            query_begin: 0,
            check_query_protein,
        })
    })
}

impl DoubleIndexedBackend for NativeBackend<'_> {
    type Query = QueryOwner;
    type Round = RoundOwner;
    type Reference = ReferenceChunk;
    type Output = ReferenceSource;
    fn capture_source_error(&self, error: io::Error) -> io::Error {
        let errors = self.references.source_errors();
        errors.capture_io::<()>(Err(error)).unwrap_err()
    }
    fn database_sequence_count(&self) -> u64 {
        self.sequence_count
    }
    fn sensitivity_round_count(&self) -> usize {
        self.rounds.len()
    }
    fn begin_master(&mut self) -> io::Result<()> {
        self.input.prepare_for_master()?;
        if self.writer.is_none() {
            self.writer = Some(NativeSearchWriter::new(&self.config, self.sequence_count)?);
        }
        let writer = self.writer.as_mut().unwrap();
        writer.set_database_letters(self.full_database_letters);
        self.format = Some(writer.intermediate_format());
        self.auxiliary = Some(crate::commands::query_output::QueryOutputs::new(
            &self.config.auxiliary_output,
            if self.config.translated_query_layout.is_some() {
                crate::basic::value::SequenceType::Nucleotide
            } else {
                crate::basic::value::SequenceType::AminoAcid
            },
            self.config.file_buffer_size,
        )?);
        Ok(())
    }
    fn load_query(&mut self, ordinal: usize) -> io::Result<Option<QueryOwner>> {
        let chunk = self.input.chunk(self.budget)?;
        let query_count = chunk.count();
        if query_count == 0 {
            return Ok(None);
        }
        let query_oid_end = chunk
            .protein_oids
            .iter()
            .copied()
            .max()
            .map_or(0, |oid| oid + 1);
        if self.check_query_protein && ordinal == 0 && chunk.translated.is_none() {
            crate::data::sequence_file::check_protein_input(
                chunk
                    .records
                    .iter()
                    .map(|record| record.sequence.as_slice()),
                self.config.tuning.ignore_warnings,
            )
            .map_err(io::Error::other)?;
        }
        let query_begin = self.query_begin;
        self.query_begin = self
            .query_begin
            .checked_add(query_count as u64)
            .ok_or_else(|| io::Error::other("Query ordinal overflow"))?;
        let dictionary = ReferenceDictionary::new(
            self.sequence_count,
            self.writer.as_ref().unwrap().needs_reference_sequences(),
        );
        Ok(Some(QueryOwner {
            chunk,
            query_count,
            query_oid_end,
            query_block: ordinal,
            query_begin,
            block_config: self.config.clone(),
            files: Vec::new(),
            dictionary,
            prepared: false,
            query_prepared: false,
            master_oids: Vec::new(),
            output_oids: Vec::new(),
            direct: None,
            aligned: vec![false; query_count],
            source_reference_block: 0,
        }))
    }
    fn query_is_selected(&self, ordinal: usize, _query: &QueryOwner) -> bool {
        self.config.tuning.multiprocessing.query_chunk < 0
            || self.config.tuning.multiprocessing.query_chunk as usize == ordinal
    }
    fn prepare_query(
        &mut self,
        q: &mut QueryOwner,
        _ordinal: usize,
        previous_reference: usize,
    ) -> io::Result<()> {
        let config = &self.config;
        let matrix = &self.matrix;
        let min_length_ratio = self.min_length_ratio;
        let prepared_matrix = &self.prepared_matrix;
        let masker = self.masker.as_deref();
        let writer = self.writer.as_mut().unwrap();
        q.master_oids = prepare_master_query_chunk(
            &mut q.chunk,
            config,
            u32::try_from(previous_reference).map_err(io::Error::other)?,
            min_length_ratio,
        )?;
        if q.query_block == 0 {
            let query = q.chunk.query(0)?;
            writer.header_bytes(&matrix, query.title, query.source.len())?;
        }
        q.block_config = config.clone();
        q.block_config.tuning.query_pre_sorted = true;
        q.block_config.tuning.search_min_length_ratio = Some(min_length_ratio);
        q.block_config.tuning.search_block.query = u32::try_from(q.query_block)
            .map_err(|_| io::Error::other("Query block index exceeds uint32"))?;
        q.block_config.tuning.search_block.self_search = config.tuning.self_search;
        q.block_config.scoring_matrix = Some(prepared_matrix.clone());
        if let Some(layout) = &mut q.block_config.translated_query_layout {
            layout.mate_indices = q.chunk.mates.clone();
        }
        // Source masks the query owner after its header, before loading any
        // reference block. Later captures reuse this prepared sequence arena.
        multiprocessing::prepare_join_only_queries(
            &mut q.chunk,
            &q.block_config,
            &matrix,
            masker.as_deref(),
        );

        q.output_oids = q.master_oids.clone();
        q.query_prepared = true;
        Ok(())
    }
    fn query_has_pending_alignments(&self, q: &QueryOwner) -> bool {
        q.aligned.iter().any(|&value| !value)
    }
    fn begin_round(&mut self, q: &mut QueryOwner, ordinal: usize) -> io::Result<RoundOwner> {
        let config = &self.config;
        let settings = self.rounds[ordinal];
        let matrix = &self.matrix;
        let letters = self.letters;
        let masker = self.masker.as_deref();
        let previously_aligned = q.aligned.iter().filter(|&&value| value).count();
        q.block_config.sensitivity = settings.sensitivity;
        q.block_config.tuning.algorithm.begin_round()?;
        q.block_config.round_linear_target = settings.linearize;
        let mut skip = crate::util::data_structures::BitVector::with_size(q.query_count);
        for (query, &value) in q.aligned.iter().enumerate() {
            if value {
                skip.set(query);
            }
        }
        q.block_config.round_query_skip = Some(std::sync::Arc::new(skip));
        // Real query setup/AUTO seed-cache construction precedes reference
        // I/O for each sensitivity iteration, as in run_query_iteration.
        prepare_unsearched_query_chunk(
            &q.block_config,
            &mut q.chunk,
            &matrix,
            masker.as_deref(),
            false,
            letters,
        )?;
        let ranking = if config.global_ranking_targets != 0 {
            let count = usize::try_from(config.global_ranking_targets)
                .ok()
                .and_then(|count| q.query_count.checked_mul(count))
                .ok_or_else(|| {
                    io::Error::other("cannot create std::vector larger than max_size()")
                })?;
            let mut table = Vec::new();
            table.try_reserve_exact(count).map_err(|_| {
                io::Error::other("cannot create std::vector larger than max_size()")
            })?;
            table.resize(count, crate::align::global_ranking::Hit::default());
            Some(super::super::blastp::NativeRankingState {
                table,
                reference_begin: 0,
                reference_oids: None,
            })
        } else {
            None
        };

        self.references.rewind();
        if config.tuning.self_search && !config.tuning.setup.lin_stage1_query {
            self.references.seek_oid(q.query_oid_end);
        }
        Ok(RoundOwner {
            ranking,
            previously_aligned,
            first_reference_chunk: true,
        })
    }
    fn load_reference(
        &mut self,
        q: &mut QueryOwner,
        round: &mut RoundOwner,
        cursor: usize,
    ) -> io::Result<Option<ReferenceChunk>> {
        let config = &self.config;
        let references = &mut self.references;
        let budget = self.budget;
        q.block_config.tuning.search_block.reference =
            u32::try_from(cursor).map_err(io::Error::other)?;
        q.source_reference_block = q.block_config.tuning.search_block.reference;
        let alias = config.tuning.self_search
            && ((!config.tuning.setup.lin_stage1_query
                && q.block_config.tuning.search_block.reference == 0)
                || (config.tuning.setup.lin_stage1_query
                    && q.block_config.tuning.search_block.reference
                        == q.block_config.tuning.search_block.query));
        let ReferenceChunk { records, oids } = if alias {
            if config.tuning.setup.lin_stage1_query {
                references.seek_oid(q.query_oid_end);
            }
            ReferenceChunk {
                records: q.chunk.records.clone(),
                oids: q.chunk.protein_oids.clone(),
            }
        } else {
            references.chunk(budget)?
        };
        if records.is_empty() {
            return Ok(None);
        }
        if round.first_reference_chunk
            && !alias
            && references.checks_protein_input()
            && oids.first() == Some(&0)
        {
            crate::data::sequence_file::check_protein_input(
                records.iter().map(|record| record.sequence.as_slice()),
                config.tuning.ignore_warnings,
            )
            .map_err(io::Error::other)?;
        }
        round.first_reference_chunk = false;

        Ok(Some(ReferenceChunk { records, oids }))
    }
    fn run_reference(
        &mut self,
        q: &mut QueryOwner,
        round: &mut RoundOwner,
        reference: ReferenceChunk,
    ) -> io::Result<ReferenceFlow> {
        let config = &self.config;
        let references = &mut self.references;
        let letters = self.letters;
        let iterated = self.iterated;
        let selected_sequence_count = self.selected_sequence_count;
        let masker = self.masker.as_deref();
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| io::Error::other("Search output is not initialized"))?;
        let format = self
            .format
            .ok_or_else(|| io::Error::other("Search format is not initialized"))?;
        let unaligned_targets = self.unaligned_targets.as_ref();

        let ReferenceChunk { records, oids } = reference;
        // The metadata owner outlives serialization of this actual block.
        if let Some(state) = round.ranking.as_mut() {
            state.reference_begin = oids.first().copied().unwrap_or(0);
            state.reference_oids = Some(oids.clone().into());
        }
        q.block_config.reference_taxon_rank_ids = if config.taxon_k != 0 {
            Some(references.species_rank_ids(&oids)?.into())
        } else {
            None
        };
        let mut capture = super::super::blastp::search_block(
            &q.block_config,
            records.clone(),
            std::mem::take(&mut q.chunk.records),
            q.chunk.translated.take(),
            letters,
            q.query_prepared,
            masker.as_deref(),
            round.ranking.take(),
        )?;
        q.chunk.records = std::mem::take(&mut capture.queries);
        q.chunk.translated = capture.translated.take();
        q.query_prepared = true;
        if !q.prepared {
            if !q.chunk.protein_oids.is_empty() {
                q.chunk.protein_oids = capture
                    .query_oids
                    .iter()
                    .map(|&local| q.chunk.protein_oids[local as usize])
                    .collect();
            }
            // Prepared Source query owners already hold these scores. A
            // capture omits recomputation and must not erase that owner data.
            if !capture.query_self_aln_scores.is_empty() {
                q.chunk.self_aln_scores = std::mem::take(&mut capture.query_self_aln_scores);
            }
            if !capture.query_oids.is_empty() {
                q.output_oids = capture
                    .query_oids
                    .iter()
                    .map(|&local| q.master_oids[local as usize])
                    .collect();
            }
        }
        if config.global_ranking_targets != 0 {
            round.ranking = capture.ranking.take();
            q.prepared = true;
            return Ok(ReferenceFlow::Continue);
        }
        if let Some(tracker) = &unaligned_targets {
            tracker.mark_capture(&capture, &oids)?;
        }
        if iterated || config.auxiliary_output.enabled() {
            for (query, matches) in capture.matches.iter().enumerate() {
                if matches.iter().any(|target| !target.hsps.is_empty()) {
                    q.aligned[query] = true;
                }
            }
        }
        if !iterated && records.len() as u64 == selected_sequence_count {
            q.dictionary.remap_direct(&mut capture, &records, &oids);
            references.load_matched_taxonomy(&mut q.dictionary)?;
            q.direct = Some((capture.matches, capture.seeded));
            q.prepared = true;
            // Source advances the loop counter before its EOF load.
            q.source_reference_block = q
                .block_config
                .tuning
                .search_block
                .reference
                .checked_add(1)
                .ok_or_else(|| io::Error::other("Reference block index exceeds uint32"))?;
            return Ok(ReferenceFlow::FinishRound);
        }
        let (file, mut output) = IntermediateFile::new(config)?;
        q.dictionary
            .write_block(&mut capture, &records, &oids, &format, &mut output)?;
        references.load_matched_taxonomy(&mut q.dictionary)?;
        for subject in &q.dictionary.entries {
            writer.reserve_archive_subject_bytes(subject.oid, &subject.title, subject.len);
        }
        std::io::Write::flush(&mut output)?;
        drop(output);
        q.files.push(file);
        q.prepared = true;

        Ok(ReferenceFlow::Continue)
    }
    fn finish_round(
        &mut self,
        q: &mut QueryOwner,
        round: RoundOwner,
        loaded_references: usize,
    ) -> io::Result<()> {
        let references = &mut self.references;
        let matrix = &self.matrix;
        let letters = self.letters;
        let iterated = self.iterated;
        let masker = self.masker.as_deref();
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| io::Error::other("Search output is not initialized"))?;
        let unaligned_targets = self.unaligned_targets.as_ref();

        q.source_reference_block = u32::try_from(loaded_references).map_err(io::Error::other)?;
        if !q.prepared {
            if let Some(oids) = prepare_unsearched_query_chunk(
                &q.block_config,
                &mut q.chunk,
                &matrix,
                masker.as_deref(),
                !q.query_prepared,
                letters,
            )? {
                q.output_oids = oids;
            }
            q.query_prepared = true;
        }
        // Empty target loads leave an empty ranking table; query reporting
        // still runs below, and the Source iteration message still occurs.
        if let Some(state) = round.ranking.filter(|_| q.prepared) {
            let intermediate = global_ranking::extend(
                &q.block_config,
                references,
                &q.chunk,
                writer,
                &matrix,
                &state.table,
                q.query_begin,
                &q.output_oids,
                iterated.then_some(&mut q.dictionary),
                &mut q.aligned,
                unaligned_targets.as_ref().map(|tracker| &tracker.aligned),
            )?;
            if let Some(file) = intermediate {
                q.files.push(file);
            }
            if iterated {
                q.source_reference_block = 0;
            }
        }
        if iterated {
            let total = q.aligned.iter().filter(|&&value| value).count();
            search_message!(
                "Aligned {}/{} queries in this iteration, {}/{} total.",
                total - round.previously_aligned,
                q.query_count,
                total,
                q.query_count
            );
        }

        Ok(())
    }
    fn finish_query(&mut self, mut q: QueryOwner) -> io::Result<QueryCompletion> {
        let config = &self.config;
        let matrix = &self.matrix;
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| io::Error::other("Search output is not initialized"))?;
        let format = self
            .format
            .ok_or_else(|| io::Error::other("Search format is not initialized"))?;

        let auxiliary = self.auxiliary.as_mut().unwrap();
        if !q.prepared {
            for query_id in 0..q.query_count {
                let query = q.chunk.query(query_id)?;
                writer.write_query_bytes(
                    query_id as u32,
                    query.oid.unwrap_or_else(|| {
                        q.query_begin
                            + q.output_oids
                                .get(query_id)
                                .copied()
                                .unwrap_or(query_id as u64)
                    }),
                    query.title,
                    query.source,
                    &query.frames,
                    query.mate,
                    query.quality,
                    query.self_aln_score,
                    Vec::new(),
                    &q.dictionary,
                    &matrix,
                )?;
            }
        } else if let Some((matches, seeded)) = q.direct.take() {
            for (query_id, matches) in matches.into_iter().enumerate() {
                if !seeded[query_id] {
                    continue;
                }
                let query = q.chunk.query(query_id)?;
                writer.write_query_bytes(
                    query_id as u32,
                    query.oid.unwrap_or_else(|| {
                        q.query_begin
                            + q.output_oids
                                .get(query_id)
                                .copied()
                                .unwrap_or(query_id as u64)
                    }),
                    query.title,
                    query.source,
                    &query.frames,
                    query.mate,
                    query.quality,
                    query.self_aln_score,
                    matches,
                    &q.dictionary,
                    &matrix,
                )?;
            }
        } else if !q.files.is_empty() {
            let readers = q
                .files
                .iter()
                .map(IntermediateFile::reader)
                .collect::<io::Result<Vec<_>>>()?;
            let mut fetcher = crate::output::join_blocks::StreamJoinFetcher::from_readers(readers)?;
            let mut exhausted = false;
            let mut last_query = None;
            let (renderer, mut sink) = writer.joined_parts(matrix);
            let render = |query_id: usize, matches, aligned| {
                let query = q.chunk.query(query_id)?;
                renderer.render_query_bytes(
                    query_id as u32,
                    query.oid.unwrap_or_else(|| {
                        q.query_begin
                            + q.output_oids
                                .get(query_id)
                                .copied()
                                .unwrap_or(query_id as u64)
                    }),
                    query.title,
                    query.source,
                    &query.frames,
                    query.mate,
                    query.quality,
                    query.self_aln_score,
                    matches,
                    &q.dictionary,
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
                |payload| {
                    let mut output = crate::commands::blastp::NativeRenderedQuery::join_task();
                    if config.report_unaligned != 0
                        && format.code != crate::output::format::FormatCode::Daa
                    {
                        for id in payload.unaligned_from..payload.query_id {
                            output.append(render(id as usize, Vec::new(), false)?);
                        }
                    }
                    let query = q.chunk.query(payload.query_id as usize)?;
                    let frames: Vec<Vec<Letter>> =
                        query.frames.iter().map(|frame| frame.to_vec()).collect();
                    let matches = merge_query(
                        &payload.buffers,
                        &frames,
                        query.source.len(),
                        &q.dictionary,
                        matrix,
                        &q.block_config,
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
                for id in last_query.map_or(0, |id| id + 1)..q.query_count as u32 {
                    sink.write(render(id as usize, Vec::new(), false)?)?;
                }
            }
        }
        for pass in [false, true] {
            for (query_id, &is_aligned) in q.aligned.iter().enumerate() {
                if is_aligned != pass {
                    continue;
                }
                let query = q.chunk.query(query_id)?;
                auxiliary.write_query_bytes(
                    is_aligned,
                    query.title,
                    query.source,
                    (!query.quality.is_empty()).then_some(query.quality),
                )?;
            }
        }

        Ok(QueryCompletion {
            reference_counter: q.source_reference_block as usize,
        })
    }
    fn stop_exists(&self) -> bool {
        crate::run::multiprocessing::stop_exists()
    }
    fn finish(mut self) -> io::Result<ReferenceSource> {
        crate::data::accession_spool::source_scope(self.references.accession_spool(), || {
            self.writer.take().unwrap().finish(&self.matrix)?;
            self.auxiliary.take().unwrap().finish()?;
            if let Some(tracker) = self.unaligned_targets.take() {
                tracker.finish_spooled()?;
            }
            Ok(self.references)
        })
    }
}
