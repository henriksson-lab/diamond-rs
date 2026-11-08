//! Actual provided-owner clustering services for the shared Source orchestrator.
use super::*;
use crate::run::config::Round;
use crate::run::double_indexed::owned::{
    self, DoubleIndexedBackend, QueryCompletion, ReferenceFlow,
};
use crate::stats::score_matrix::ScoreMatrix;
use std::sync::{Arc, Mutex};

struct Backend<'a> {
    config: BlastpConfig,
    database: &'a [FastaRecord],
    input: &'a [FastaRecord],
    owner: &'a ClusterSearchOwner,
    options: ClusterSearchOptions<'a>,
    rounds: Vec<Round>,
    iterated: bool,
    min_length_ratio: f64,
    budget: u64,
    database_letters: u64,
    query_cursor: usize,
    matrix: ScoreMatrix,
    fallback_masker: Option<Arc<crate::masking::tantan::TantanMasker>>,
    auxiliary: Option<crate::commands::query_output::QueryOutputs>,
    unaligned_targets: Option<super::super::UnalignedTargets>,
    global_aligned: Option<Mutex<crate::util::data_structures::BitVector>>,
    edges: Vec<InMemorySearchEdge>,
}
struct QueryOwner {
    query_range: Range<usize>,
    query_block: usize,
    records: Vec<FastaRecord>,
    query_oids: Vec<u64>,
    self_scores: Vec<f64>,
    config: BlastpConfig,
    prepared: bool,
    query_prepared: bool,
    dictionary: ReferenceDictionary,
    files: Vec<IntermediateFile>,
    direct: Option<Vec<Vec<crate::align::hsp::Match>>>,
    aligned: Vec<bool>,
    source_reference_block: u32,
}
struct RoundOwner<'a> {
    ranking: Option<blastp::NativeRankingState>,
    options: ClusterSearchOptions<'a>,
    previously_aligned: usize,
    reference_cursor: usize,
}
struct ReferenceOwner {
    targets: Vec<FastaRecord>,
    target_oids: Vec<u64>,
    alias: bool,
}
fn next_range(records: &[FastaRecord], begin: usize, budget: u64) -> Option<Range<usize>> {
    if begin >= records.len() {
        return None;
    }
    let mut end = begin;
    let mut letters = 0u64;
    while end < records.len() && (end == begin || letters < budget) {
        letters = letters.saturating_add(records[end].sequence.len() as u64);
        end += 1;
    }
    Some(begin..end)
}

pub(super) fn run<'a>(
    config: &BlastpConfig,
    database: &'a [FastaRecord],
    input: &'a [FastaRecord],
    budget: u64,
    database_letters: u64,
    options: ClusterSearchOptions<'a>,
    owner: &'a ClusterSearchOwner,
    rounds: Vec<Round>,
    min_length_ratio: f64,
) -> io::Result<Vec<InMemorySearchEdge>> {
    let matrix = blastp::native_score_matrix(config, database_letters)?;
    let mut config = config.clone();
    config.tuning.search_min_length_ratio = Some(min_length_ratio);
    blastp::begin_algorithm_scope(&mut config)?;
    config.approx_min_id = options.approx_min_id;
    config.query_or_subject_cover = options.query_or_target_cover;
    config.block_size = None;
    config.scoring_matrix = Some(Arc::new(matrix.clone()));
    config.tuning.query_pre_sorted = true;
    let fallback_masker = (options.masker.is_none()
        && config.masking == crate::masking::MaskingMode::Tantan)
        .then(|| blastp::native_tantan_masker(&config, &matrix));
    let unaligned_targets = if let Some(accessions) = &owner.accessions {
        super::super::UnalignedTargets::from_accessions(&config, owner.sequence_count, |oid| {
            accessions
                .get(oid as usize)
                .cloned()
                .ok_or_else(|| io::Error::other("Target tracking OID outside original owner."))
        })?
    } else {
        // Source's provided spill/subdatabase owner was constructed without
        // OID_TO_ACC_MAPPING. Preserve its missing-spool finalization boundary.
        super::super::UnalignedTargets::with_count(&config, owner.sequence_count)
    };
    let global_aligned =
        (unaligned_targets.is_some() && config.global_ranking_targets != 0).then(|| {
            Mutex::new(crate::util::data_structures::BitVector::with_size(
                database.len(),
            ))
        });
    crate::data::accession_spool::source_scope(
        unaligned_targets
            .as_ref()
            .and_then(|tracker| tracker.accession_spool.clone()),
        || {
            owned::run(Backend {
                config,
                database,
                input,
                owner,
                options,
                iterated: rounds.len() > 1,
                rounds,
                min_length_ratio,
                budget,
                database_letters,
                query_cursor: 0,
                matrix,
                fallback_masker,
                auxiliary: None,
                unaligned_targets,
                global_aligned,
                edges: Vec::new(),
            })
        },
    )
}

impl<'a> DoubleIndexedBackend for Backend<'a> {
    type Query = QueryOwner;
    type Round = RoundOwner<'a>;
    type Reference = ReferenceOwner;
    type Output = Vec<InMemorySearchEdge>;
    fn capture_source_error(&self, error: io::Error) -> io::Error {
        let errors = self
            .unaligned_targets
            .as_ref()
            .and_then(|tracker| tracker.accession_spool.as_ref())
            .map_or_else(Default::default, |spool| spool.io_errors());
        errors.capture_io::<()>(Err(error)).unwrap_err()
    }
    fn database_sequence_count(&self) -> u64 {
        self.owner.sequence_count as u64
    }
    fn sensitivity_round_count(&self) -> usize {
        self.rounds.len()
    }
    fn initial_reference_counter(&self) -> usize {
        self.config.tuning.search_block.reference as usize
    }
    fn begin_master(&mut self) -> io::Result<()> {
        self.auxiliary = Some(crate::commands::query_output::QueryOutputs::new(
            &self.config.auxiliary_output,
            crate::basic::value::SequenceType::AminoAcid,
            self.config.file_buffer_size,
        )?);
        Ok(())
    }
    fn load_query(&mut self, ordinal: usize) -> io::Result<Option<QueryOwner>> {
        let Some(query_range) = next_range(self.input, self.query_cursor, self.budget) else {
            return Ok(None);
        };
        self.query_cursor = query_range.end;
        let records = self.input[query_range.clone()].to_vec();
        let n = records.len();
        Ok(Some(QueryOwner {
            query_oids: (query_range.start as u64..query_range.end as u64).collect(),
            query_range,
            query_block: ordinal,
            records,
            self_scores: Vec::new(),
            config: self.config.clone(),
            prepared: false,
            query_prepared: false,
            dictionary: ReferenceDictionary::new(self.owner.sequence_count as u64, false),
            files: Vec::new(),
            direct: None,
            aligned: vec![false; n],
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
        q.config.tuning.search_block.query =
            u32::try_from(q.query_block).map_err(io::Error::other)?;
        q.config.tuning.search_block.reference =
            u32::try_from(previous_reference).map_err(io::Error::other)?;
        q.config.tuning.search_block.self_search = self.options.self_search;
        let mut chunk = super::super::QueryChunk {
            records: std::mem::take(&mut q.records),
            protein_oids: std::mem::take(&mut q.query_oids),
            self_aln_scores: Vec::new(),
            translated: None,
            mates: Vec::new(),
        };
        super::super::prepare_master_query_chunk(
            &mut chunk,
            &q.config,
            u32::try_from(previous_reference).map_err(io::Error::other)?,
            self.min_length_ratio,
        )?;
        let masker = self.options.masker.or(self.fallback_masker.as_deref());
        super::super::multiprocessing::prepare_join_only_queries(
            &mut chunk,
            &q.config,
            &self.matrix,
            masker,
        );
        q.records = chunk.records;
        q.query_oids = chunk.protein_oids;
        q.self_scores = chunk.self_aln_scores;
        q.query_prepared = true;
        Ok(())
    }
    fn query_has_pending_alignments(&self, q: &QueryOwner) -> bool {
        q.aligned.iter().any(|&value| !value)
    }
    fn begin_round(&mut self, q: &mut QueryOwner, ordinal: usize) -> io::Result<RoundOwner<'a>> {
        let settings = self.rounds[ordinal];
        q.config.sensitivity = settings.sensitivity;
        q.config.round_linear_target = settings.linearize;
        q.config.tuning.algorithm.begin_round()?;
        let mut skip = crate::util::data_structures::BitVector::with_size(q.records.len());
        for (id, &aligned) in q.aligned.iter().enumerate() {
            if aligned {
                skip.set(id);
            }
        }
        q.config.round_query_skip = Some(Arc::new(skip));
        let mut options = self.options;
        options.linear_target |= settings.linearize;
        // Build the real query cache with the same explicit overrides used by
        // the reference kernel, before the first reference is loaded.
        let mut cache_config = q.config.clone();
        cache_config.tuning.setup.freq_masking = options.freq_masking;
        cache_config.tuning.setup.lin_stage1_query |= options.linear_query;
        cache_config.round_linear_target |= options.linear_target;
        let mut chunk = super::super::QueryChunk {
            records: std::mem::take(&mut q.records),
            protein_oids: q.query_oids.clone(),
            self_aln_scores: q.self_scores.clone(),
            translated: None,
            mates: Vec::new(),
        };
        super::super::prepare_query_cache_with_overrides(
            &cache_config,
            &mut chunk,
            &self.matrix,
            options.masker.or(self.fallback_masker.as_deref()),
            false,
            self.database_letters,
            options.soft_masking,
            options.freq_masking,
        )?;
        q.records = chunk.records;
        let ranking = if q.config.global_ranking_targets != 0 {
            let n = usize::try_from(q.config.global_ranking_targets)
                .ok()
                .and_then(|width| q.records.len().checked_mul(width))
                .ok_or_else(|| {
                    io::Error::other("cannot create std::vector larger than max_size()")
                })?;
            let mut table = Vec::new();
            table.try_reserve_exact(n).map_err(|_| {
                io::Error::other("cannot create std::vector larger than max_size()")
            })?;
            table.resize(n, crate::align::global_ranking::Hit::default());
            Some(blastp::NativeRankingState {
                table,
                reference_begin: 0,
                reference_oids: None,
            })
        } else {
            None
        };
        Ok(RoundOwner {
            ranking,
            options,
            previously_aligned: q.aligned.iter().filter(|&&value| value).count(),
            reference_cursor: if q.config.tuning.self_search
                && !q.config.tuning.setup.lin_stage1_query
            {
                q.query_range.end
            } else {
                0
            },
        })
    }
    fn load_reference(
        &mut self,
        q: &mut QueryOwner,
        round: &mut RoundOwner<'a>,
        cursor: usize,
    ) -> io::Result<Option<ReferenceOwner>> {
        q.config.tuning.search_block = SearchBlockContext {
            query: u32::try_from(q.query_block).map_err(io::Error::other)?,
            reference: u32::try_from(cursor).map_err(io::Error::other)?,
            self_search: self.options.self_search,
        };
        q.source_reference_block = q.config.tuning.search_block.reference;
        let alias = q.config.tuning.self_search
            && ((!q.config.tuning.setup.lin_stage1_query && cursor == 0)
                || (q.config.tuning.setup.lin_stage1_query && cursor == q.query_block));
        if alias {
            if q.config.tuning.setup.lin_stage1_query {
                round.reference_cursor = q.query_range.end;
            }
            Ok(Some(ReferenceOwner {
                targets: q.records.clone(),
                target_oids: q.query_oids.clone(),
                alias: true,
            }))
        } else {
            let Some(range) = next_range(self.database, round.reference_cursor, self.budget) else {
                return Ok(None);
            };
            round.reference_cursor = range.end;
            Ok(Some(ReferenceOwner {
                targets: self.database[range.clone()].to_vec(),
                target_oids: (range.start as u64..range.end as u64).collect(),
                alias: false,
            }))
        }
    }
    fn run_reference(
        &mut self,
        q: &mut QueryOwner,
        round: &mut RoundOwner<'a>,
        reference: ReferenceOwner,
    ) -> io::Result<ReferenceFlow> {
        let ReferenceOwner {
            targets,
            target_oids,
            alias,
        } = reference;
        let database = self.database;
        let database_letters = self.database_letters;
        let iterated = self.iterated;
        let owner = self.owner;
        let unaligned_targets = self.unaligned_targets.as_ref();
        let format = OutputFormat::edge();
        if let Some(state) = &mut round.ranking {
            state.reference_begin = target_oids.first().copied().unwrap_or(0);
            state.reference_oids = Some(target_oids.clone().into());
        }
        let mut capture = blastp::search_cluster_block_ranked(
            &q.config,
            targets.clone(),
            if alias {
                None
            } else {
                Some(std::mem::take(&mut q.records))
            },
            database_letters,
            q.query_prepared,
            round.options,
            round.ranking.take(),
        )?;
        if !q.prepared && !alias {
            q.query_oids = capture
                .query_oids
                .iter()
                .map(|&local| q.query_oids[local as usize])
                .collect();
        }
        if !q.prepared {
            if !capture.query_self_aln_scores.is_empty() {
                q.self_scores = std::mem::take(&mut capture.query_self_aln_scores);
            }
        }
        q.records = std::mem::take(&mut capture.queries);
        q.prepared = true;
        q.query_prepared = true;
        if let Some(tracker) = &unaligned_targets {
            let physical_oids = target_oids
                .iter()
                .map(|&oid| {
                    if owner.accessions.is_some() {
                        owner.original_oids[oid as usize]
                    } else {
                        oid
                    }
                })
                .collect::<Vec<_>>();
            tracker.mark_capture(&capture, &physical_oids)?;
        }
        if q.config.global_ranking_targets != 0 {
            round.ranking = capture.ranking.take();
            return Ok(ReferenceFlow::Continue);
        }
        if iterated || q.config.auxiliary_output.enabled() {
            for (id, matches) in capture.matches.iter().enumerate() {
                q.aligned[id] |= matches.iter().any(|target| !target.hsps.is_empty());
            }
        }
        if !iterated && targets.len() == database.len() {
            q.dictionary
                .remap_direct(&mut capture, &targets, &target_oids);
            q.direct = Some(capture.matches);
            return Ok(ReferenceFlow::FinishRound);
        }
        let (file, mut writer) = IntermediateFile::new(&q.config)?;
        q.dictionary
            .write_block(&mut capture, &targets, &target_oids, &format, &mut writer)?;
        writer.flush()?;
        drop(writer);
        q.files.push(file);

        Ok(ReferenceFlow::Continue)
    }
    fn finish_round(
        &mut self,
        q: &mut QueryOwner,
        round: RoundOwner<'a>,
        loaded_references: usize,
    ) -> io::Result<()> {
        let database = self.database;
        let database_letters = self.database_letters;
        let options = self.options;
        let iterated = self.iterated;
        let matrix = &self.matrix;
        let edges = &mut self.edges;
        let global_aligned = &self.global_aligned;
        q.config.tuning.search_block.reference =
            u32::try_from(loaded_references).map_err(io::Error::other)?;
        q.source_reference_block = q.config.tuning.search_block.reference;
        if !q.prepared {
            let mut chunk = super::super::QueryChunk {
                protein_oids: q.query_oids.clone(),
                self_aln_scores: q.self_scores.clone(),
                records: std::mem::take(&mut q.records),
                translated: None,
                mates: Vec::new(),
            };
            let fallback_masker = (options.masker.is_none()
                && q.config.masking == crate::masking::MaskingMode::Tantan)
                .then(|| blastp::native_tantan_masker(&q.config, &matrix));
            if let Some(order) = super::super::prepare_unsearched_query_chunk(
                &q.config,
                &mut chunk,
                &matrix,
                options.masker.or(fallback_masker.as_deref()),
                !q.query_prepared,
                database_letters,
            )? {
                q.query_oids = order
                    .iter()
                    .map(|&local| q.query_oids[local as usize])
                    .collect();
            }
            q.self_scores = chunk.self_aln_scores;
            q.records = chunk.records;
            q.query_prepared = true;
        }
        if let Some(state) = round.ranking {
            let chunk = super::super::QueryChunk {
                protein_oids: q.query_oids.clone(),
                self_aln_scores: q.self_scores.clone(),
                records: std::mem::take(&mut q.records),
                translated: None,
                mates: Vec::new(),
            };
            let mut refs = super::super::global_ranking::SliceReferences {
                records: database,
                physical_totals: (self.owner.sequence_count as u64, self.owner.letters),
            };
            let mut edge_config = q.config.clone();
            edge_config.outfmt = vec!["edge".into()];
            let mut writer = blastp::NativeSearchWriter::new_with_consumer(
                &edge_config,
                self.owner.sequence_count as u64,
                Box::new(io::sink()),
            )?;
            if iterated {
                if let Some(file) = super::super::global_ranking::extend(
                    &edge_config,
                    &mut refs,
                    &chunk,
                    &mut writer,
                    &matrix,
                    &state.table,
                    q.query_range.start as u64,
                    &q.query_oids,
                    Some(&mut q.dictionary),
                    &mut q.aligned,
                    global_aligned.as_ref(),
                )? {
                    q.files.push(file);
                }
            } else {
                edges.extend(super::super::global_ranking::extend_captured(
                    &edge_config,
                    &mut refs,
                    &chunk,
                    &mut writer,
                    &matrix,
                    &state.table,
                    q.query_range.start as u64,
                    &q.query_oids,
                    &mut q.aligned,
                    global_aligned.as_ref(),
                )?);
            }
            q.records = chunk.records;
            if iterated {
                q.source_reference_block = 0;
            }
        }
        if iterated {
            let total = q.aligned.iter().filter(|&&value| value).count();
            search_message!(
                "Aligned {}/{} queries in this iteration, {}/{} total.",
                total - round.previously_aligned,
                q.records.len(),
                total,
                q.records.len()
            );
        }

        Ok(())
    }
    fn finish_query(&mut self, mut q: QueryOwner) -> io::Result<QueryCompletion> {
        let auxiliary = self.auxiliary.as_mut().unwrap();
        if q.config.global_ranking_targets != 0 && !self.iterated {
            write_auxiliary_queries(auxiliary, &q.records, &q.aligned)?;
            return Ok(QueryCompletion {
                reference_counter: q.source_reference_block as usize,
            });
        }
        let database = self.database;
        let owner = self.owner;
        let matrix = &self.matrix;
        let edges = &mut self.edges;
        let format = OutputFormat::edge();
        let mut merged = if let Some(matches) = q.direct.take() {
            for query in &matches {
                crate::basic::statistics::merge_statistics(&blastp::generated_output_statistics(
                    query,
                ));
            }
            matches
        } else {
            let readers = q
                .files
                .iter()
                .map(IntermediateFile::reader)
                .collect::<io::Result<Vec<_>>>()?;
            let mut fetcher = crate::output::join_blocks::StreamJoinFetcher::from_readers(readers)?;
            let mut result = Vec::with_capacity(q.records.len());
            for (query_id, query) in q.records.iter().enumerate() {
                if fetcher.next() == query_id as u32 {
                    fetcher.fetch_next()?;
                    if q.config.taxon_k != 0 && fetcher.buf.iter().any(|buffer| !buffer.is_empty())
                    {
                        // Source asks db->taxids(target_oid) before culling
                        // the first actual joined target, never on direct/no-hit output.
                        let mut input = std::io::Cursor::new(
                            fetcher
                                .buf
                                .iter()
                                .find(|buffer| !buffer.is_empty())
                                .unwrap(),
                        );
                        let record = crate::output::intermediate::IntermediateRecord::read(
                            &mut input, &format,
                        )?;
                        let local = q.dictionary.original_oid(record.target_dict_id as usize);
                        let oid = owner
                            .original_oids
                            .get(local as usize)
                            .copied()
                            .unwrap_or(local);
                        return Err(owner.taxonomy_access(oid));
                    }
                    result.push(merge_query(
                        &fetcher.buf,
                        &[query.sequence.clone()],
                        query.sequence.len(),
                        &q.dictionary,
                        &matrix,
                        &q.config,
                        format,
                    )?);
                } else {
                    result.push(Vec::new());
                }
            }
            fetcher.finish();
            result
        };
        for (query_id, matches) in merged.iter_mut().enumerate() {
            let query = &q.records[query_id].sequence;
            for target in matches {
                let Some(hsp) = target.hsps.first() else {
                    continue;
                };
                let target_id = target.target_block_id as usize;
                let subject_oid = q.dictionary.original_oid(target_id);
                let subject = database
                    .get(subject_oid as usize)
                    .map_or(&[][..], |record| record.sequence.as_slice());
                edges.push(InMemorySearchEdge {
                    edge: crate::output::edge::EdgeData {
                        query: q.query_oids[query_id],
                        target: q.dictionary.original_oid(target_id),
                        qcovhsp: hsp.query_cover_percent(query.len() as u32) as f32,
                        scovhsp: hsp.subject_cover_percent(q.dictionary.length(target_id) as u32)
                            as f32,
                        evalue: hsp.corrected_bit_score,
                    },
                    approx_id: hsp.approx_id_percent(query, subject),
                    identity: hsp.id_percent(),
                });
            }
        }

        write_auxiliary_queries(auxiliary, &q.records, &q.aligned)?;
        Ok(QueryCompletion {
            reference_counter: q.source_reference_block as usize,
        })
    }
    fn stop_exists(&self) -> bool {
        crate::run::multiprocessing::stop_exists()
    }
    fn finish(mut self) -> io::Result<Vec<InMemorySearchEdge>> {
        crate::data::accession_spool::source_scope(
            self.unaligned_targets
                .as_ref()
                .and_then(|tracker| tracker.accession_spool.clone()),
            || {
                self.auxiliary.take().unwrap().finish()?;
                if let Some(tracker) = self.unaligned_targets {
                    if let Some(compact) = self.global_aligned {
                        let compact = compact
                            .into_inner()
                            .map_err(|_| io::Error::other("Target tracking lock poisoned."))?;
                        let mut full = tracker
                            .aligned
                            .lock()
                            .map_err(|_| io::Error::other("Target tracking lock poisoned."))?;
                        for local in 0..compact.size() {
                            if compact.get(local as usize) {
                                let oid = if self.owner.accessions.is_some() {
                                    self.owner.original_oids[local as usize]
                                } else {
                                    local
                                };
                                if oid >= full.size() {
                                    return Err(io::Error::other(
                                        "Target tracking OID outside database.",
                                    ));
                                }
                                full.set(oid as usize);
                            }
                        }
                    }
                    tracker.finish_spooled()?;
                }
                Ok(self.edges)
            },
        )
    }
}
