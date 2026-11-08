//! Concrete native database, parallel-DP and output adapters for global ranking.
//! The enclosing block driver retains the table across the reference scan;
//! this phase reloads only its ranked OID union and performs full extension.

use crate::output::intermediate::{IntermediateRecord, OutputFormat};
use std::io::{self, Write};

use rayon::prelude::*;

use super::{QueryChunk, ReferenceDictionary, ReferenceSource};
use crate::align::gapped_filter::{SeedHit, SeedHitList, TargetScore};
use crate::align::global_ranking::{
    ExtendConfig, GlobalRankingBackend, Hit, QueryOutput, RankedTarget, SequenceLoadConfig,
};
use crate::align::target::GappedScoreConfig;
use crate::align::ungapped::UngappedStageConfig;
use crate::basic::statistics::Statistics;
use crate::basic::value::{BlockId, SequenceType, LETTER_MASK, MASK_LETTER, SEED_MASK};
use crate::commands::blastp::{BlastpConfig, NativeSearchWriter};
use crate::commands::search_output::SearchFormat;
use crate::data::block::Block;
use crate::dp::swipe::{Flags, HspValues};
use crate::masking::MaskingMode;
use crate::search::sensitivity::ExtensionMode;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::data_structures::BitVector;

pub(super) trait RankingReferences {
    fn totals(&self) -> (u64, u64);
    /// Original physical owner totals for Output::Info; compact rank indices
    /// and scoring matrix database-size overrides remain independent.
    fn metadata_totals(&self) -> (u64, u64) {
        self.totals()
    }
    fn rewind(&mut self);
    fn init_seq_access(&mut self);
    fn checks_protein_input(&self) -> bool;
    fn record(&mut self, oid: u64) -> io::Result<crate::data::fasta::ByteFastaRecord>;
    fn load_matched_taxonomy(&mut self, dictionary: &mut ReferenceDictionary) -> io::Result<()>;
}

impl RankingReferences for ReferenceSource {
    fn checks_protein_input(&self) -> bool {
        ReferenceSource::checks_protein_input(self)
    }
    fn totals(&self) -> (u64, u64) {
        ReferenceSource::totals(self)
    }
    fn rewind(&mut self) {
        ReferenceSource::rewind(self);
    }
    fn init_seq_access(&mut self) {
        match self {
            ReferenceSource::Database { file, .. } => {
                file.init_seq_access();
                file.init_seqinfo_access();
            }
            ReferenceSource::Blast { file, .. } => {
                file.init_seq_access();
                file.init_seqinfo_access();
            }
            ReferenceSource::Fasta(_) | ReferenceSource::Memory { .. } => (),
        }
    }
    fn record(&mut self, oid: u64) -> io::Result<crate::data::fasta::ByteFastaRecord> {
        ReferenceSource::record(self, oid)
    }
    fn load_matched_taxonomy(&mut self, dictionary: &mut ReferenceDictionary) -> io::Result<()> {
        ReferenceSource::load_matched_taxonomy(self, dictionary)
    }
}

/// Borrow the actual clustering database owner; only the ranked OID union is
/// copied into the DP target block, as in Source's selective database load.
pub(super) struct SliceReferences<'a> {
    pub records: &'a [crate::data::fasta::ByteFastaRecord],
    pub physical_totals: (u64, u64),
}
impl RankingReferences for SliceReferences<'_> {
    fn metadata_totals(&self) -> (u64, u64) {
        self.physical_totals
    }
    fn checks_protein_input(&self) -> bool {
        false
    }
    fn totals(&self) -> (u64, u64) {
        (
            self.records.len() as u64,
            self.records
                .iter()
                .map(|record| record.sequence.len() as u64)
                .sum(),
        )
    }
    fn rewind(&mut self) {}
    fn init_seq_access(&mut self) {}
    fn record(&mut self, oid: u64) -> io::Result<crate::data::fasta::ByteFastaRecord> {
        self.records
            .get(oid as usize)
            .cloned()
            .ok_or_else(|| io::Error::other("Ranked target OID outside database."))
    }
    fn load_matched_taxonomy(&mut self, _dictionary: &mut ReferenceDictionary) -> io::Result<()> {
        Ok(())
    }
}

struct NativeRankingBackend<'a> {
    io_errors: crate::util::io::IoErrorOwner,
    database_totals: (u64, u64),
    physical_totals: (u64, u64),
    edge_capture:
        Option<&'a std::sync::Mutex<Vec<Vec<crate::commands::blastp::InMemorySearchEdge>>>>,
    ordered_next: BlockId,
    aligned_targets: Option<&'a std::sync::Mutex<crate::util::data_structures::BitVector>>,
    config: &'a BlastpConfig,
    references: &'a mut dyn RankingReferences,
    queries: &'a QueryChunk,
    writer: &'a mut NativeSearchWriter,
    matrix: &'a ScoreMatrix,
    query_begin: u64,
    query_oids: &'a [u64],
    targets: Block,
    dictionary: ReferenceDictionary,
    format: Option<SearchFormat>,
    tabular: Option<crate::output::blast_tab_format::ByteBlastTabFormat>,
    pool: rayon::ThreadPool,
    random_access: Option<(BlockId, BlockId)>,
    persistent_dictionary: Option<&'a mut ReferenceDictionary>,
    intermediate_ids: Vec<u32>,
    intermediate_format: Option<OutputFormat>,
    intermediate_output: Option<std::io::BufWriter<std::fs::File>>,
    intermediate_finished: bool,
}

/// Immutable per-query worker state; no file writer or database cursor crosses
/// the parallel boundary. Results return to the translated ordered sink.
struct Worker<'a> {
    edge_capture:
        Option<&'a std::sync::Mutex<Vec<Vec<crate::commands::blastp::InMemorySearchEdge>>>>,
    aligned_targets: Option<&'a std::sync::Mutex<crate::util::data_structures::BitVector>>,
    config: &'a BlastpConfig,
    queries: &'a QueryChunk,
    database_letters: u64,
    matrix: &'a ScoreMatrix,
    targets: &'a Block,
    dictionary: &'a ReferenceDictionary,
    format: Option<SearchFormat>,
    tabular: Option<&'a crate::output::blast_tab_format::ByteBlastTabFormat>,
    query_begin: u64,
    query_oids: &'a [u64],
    intermediate: Option<(&'a [u32], OutputFormat)>,
}

impl Worker<'_> {
    fn align(&self, query_id: BlockId, targets: &[RankedTarget]) -> io::Result<QueryOutput> {
        let query = self.queries.query(query_id as usize)?;
        let contexts = query.frames.len();
        let mut seeds = SeedHitList::default();
        seeds.seed_hits.reserve(targets.len() as u64, 0);
        for (i, target) in targets.iter().enumerate() {
            seeds.target_block_ids.push(target.block_id);
            seeds.target_scores.push(TargetScore {
                target: i as u32,
                score: target.score,
            });
            seeds.seed_hits.next();
            seeds.seed_hits.push_item(SeedHit::new(
                0,
                0,
                target.score as i32,
                target.context as u32,
            ));
        }
        let cbs: Vec<Vec<i8>> = if self.config.comp_based_stats.hauser() {
            query
                .frames
                .iter()
                .map(|frame| {
                    crate::stats::cbs::hauser_correction_window(
                        frame,
                        self.matrix,
                        self.config.tuning.cbs_window as u32 as usize,
                    )
                })
                .collect()
        } else {
            Vec::new()
        };
        let composition = crate::stats::cbs::compute_composition(query.frames[0]);
        let cfg = GappedScoreConfig {
            anchored_swipe: self.config.tuning.setup.anchored_swipe,
            log_query: self.config.tuning.log_query,
            target_parallel_verbosity: self.config.tuning.target_parallel_verbosity,
            threads_align: self.config.tuning.threads_align,
            swipe_error: Some(Default::default()),
            add_self_aln: self.config.tuning.add_self_aln,
            self_: self.config.tuning.self_search,
            current_query_block: self.config.tuning.search_block.query,
            current_ref_block: if self.intermediate.is_some() {
                0
            } else {
                self.config.tuning.search_block.reference
            },
            clustering_filter: self.config.tuning.output.clustering_filter.clone(),
            clustering_query_self_score: query.self_aln_score,
            clustering_target_self_scores: self
                .config
                .tuning
                .output
                .clustering_filter
                .as_ref()
                .map(|_| {
                    self.dictionary
                        .entries
                        .iter()
                        .map(|entry| entry.self_aln_score.unwrap_or(0.0))
                        .collect::<Vec<_>>()
                        .into()
                }),
            cluster_threshold_present: self
                .config
                .tuning
                .output
                .clustering_filter
                .as_ref()
                .is_some_and(|filter| filter.minimum.is_some()),
            cbs_matrix_scale: self.config.tuning.matrix_scale,
            matrix_adjust_tolerance: self.config.tuning.cbs_err_tolerance,
            matrix_adjust_max_iterations: self.config.tuning.cbs_it_limit,
            query_contexts: contexts,
            query_translated: contexts == 6,
            swipe_task_size: self.config.tuning.swipe_task_size,
            transcript_len_estimate: self.config.tuning.transcript_len_estimate,
            band_bin: self.config.tuning.band_bin,
            col_bin: self.config.tuning.col_bin,
            narrow_band_cov: self.config.tuning.narrow_band_cov,
            narrow_band_factor: self.config.tuning.narrow_band_factor,
            min_band_overlap: self.config.tuning.min_band_overlap,
            cutoff_score_8bit: self.config.tuning.cutoff_score_8bit,
            max_swipe_dp: self.config.tuning.max_swipe_dp,
            approx_backtrace: self.config.tuning.approx_backtrace,
            tantan_min_mask_prob: self.config.tuning.tantan_min_mask_prob,
            cbs_window: self.config.tuning.cbs_window,
            tantan_masker: self.config.tuning.runtime.tantan(),
            padding: self.config.tuning.padding.unwrap_or(0),
            no_ranking: self.config.tuning.no_ranking,
            mapany: self.config.tuning.mapany,
            target_hard_cap: (self.config.tuning.target_hard_cap != 0)
                .then_some(self.config.tuning.target_hard_cap),
            ranking_score_drop_factor: self.config.tuning.ranking_score_drop_factor,
            ranking_cutoff_bitscore: self.config.tuning.ranking_cutoff_bitscore,
            inner_culling_overlap: self.config.culling_overlap,
            query_or_target_cover: self.config.query_or_subject_cover,
            comp_based_stats_hauser: self.config.comp_based_stats.hauser(),
            query_cover: self.config.query_cover,
            subject_cover: self.config.subject_cover,
            no_self_hits: self.config.no_self_hits,
            max_evalue: self.config.max_evalue,
            min_bit_score: self.config.min_bit_score,
            min_id: self.config.min_id,
            approx_min_id: self.config.approx_min_id,
            max_target_seqs: self.config.max_target_seqs,
            max_hsps: self.config.max_hsps,
            ext_chunk_size: self.config.ext_chunk_size,
            toppercent: self.config.toppercent,
            global_ranking_targets: self.config.global_ranking_targets,
            sensitivity: self.config.sensitivity,
            ..Default::default()
        };
        let ungapped = UngappedStageConfig {
            anchored_swipe: self.config.tuning.setup.anchored_swipe,
            hamming_ext: crate::commands::blastp::native_hamming_config(
                self.config,
                self.config.approx_min_id,
                self.config.query_or_subject_cover,
            ),
            mutual_cover: self.config.tuning.mutual_cover.is_some(),
            lin_stage1_query: self.config.tuning.setup.lin_stage1_query,
            lin_stage1_target: self.config.tuning.setup.lin_stage1_target
                || self.config.round_linear_target,
            log_extend: self.config.tuning.log_extend,
            chaining: self.config.tuning.chaining,
            cbs_matrix_scale: self.config.tuning.matrix_scale,
            query_contexts: contexts,
            query_translated: contexts == 6,
            comp_based_stats: self.config.comp_based_stats,
            matrix_adjust_thresholds: self.config.tuning.matrix_adjust_thresholds,
            matrix_adjust_tolerance: self.config.tuning.cbs_err_tolerance,
            matrix_adjust_max_iterations: self.config.tuning.cbs_it_limit,
            xdrop: self.matrix.rawscore_int(self.config.ungapped_xdrop_bits),
            ..Default::default()
        };
        let values = self.format.map_or_else(
            || {
                self.tabular
                    .unwrap()
                    .format
                    .fields
                    .iter()
                    .fold(HspValues::NONE, |values, field| values | field.hsp_values())
            },
            SearchFormat::hsp_values,
        );
        let values = values
            | self
                .config
                .tuning
                .output
                .clustering_filter
                .as_ref()
                .map_or(HspValues::NONE, |filter| filter.format.hsp_values);
        let mut statistics = Statistics::new();
        let mut matches = crate::align::target::extend_seed_hit_list(
            query_id,
            &query.frames,
            query.title,
            query.source.len() as i32,
            &cbs,
            &composition,
            self.targets,
            &mut statistics,
            Flags::FULL_MATRIX,
            &mut seeds,
            ExtensionMode::Full,
            &cfg,
            &ungapped,
            self.matrix,
            values,
            |_, _| -1,
            |_, _| -1,
            self.matrix.gap_open(),
            self.matrix.gap_extend(),
            self.matrix.rawscore_int(12.0),
            200,
        );
        cfg.check_swipe_error()?;
        crate::basic::statistics::merge_statistics(&statistics);
        crate::commands::blastp::validate_added_self_targets(&matches, self.targets.seqs().len())?;
        if let Some(filter) = &self.config.tuning.output.clustering_filter {
            filter.check().map_err(io::Error::other)?;
        }
        if self.config.no_self_hits {
            matches.retain(|target| {
                let id = target.target_block_id as usize;
                let sequence = self.targets.seqs().get(id);
                target.is_synthetic_self_alignment()
                    || query.title != self.dictionary.entries[id].title
                    || sequence.len() != query.frames[0].len()
                    || sequence
                        .iter()
                        .zip(query.frames[0])
                        .any(|(&a, &b)| a & LETTER_MASK != b & LETTER_MASK)
            });
        }
        if let Some(tracker) = self.aligned_targets {
            let mut aligned = tracker
                .lock()
                .map_err(|_| io::Error::other("Target tracking lock poisoned."))?;
            for target in &matches {
                let oid = self.dictionary.entries[target.target_block_id as usize].oid;
                if oid >= aligned.size() {
                    return Err(io::Error::other("Target tracking OID outside database."));
                }
                aligned.set(oid as usize);
            }
        }
        let output_statistics = crate::commands::blastp::generated_output_statistics(&matches);
        let aligned = matches.iter().any(|target| !target.hsps.is_empty());
        if let Some(capture) = self.edge_capture {
            let oid = query.oid.unwrap_or_else(|| {
                self.query_begin
                    + self
                        .query_oids
                        .get(query_id as usize)
                        .copied()
                        .unwrap_or(query_id as u64)
            });
            let mut edges = Vec::with_capacity(matches.len());
            for target in &matches {
                let Some(hsp) = target.hsps.first() else {
                    continue;
                };
                let id = target.target_block_id as usize;
                let subject = self.targets.seqs().get(id);
                edges.push(crate::commands::blastp::InMemorySearchEdge {
                    edge: crate::output::edge::EdgeData {
                        query: oid,
                        target: self.dictionary.entries[id].oid,
                        qcovhsp: hsp.query_cover_percent(query.source.len() as u32) as f32,
                        scovhsp: hsp.subject_cover_percent(subject.len() as u32) as f32,
                        evalue: hsp.corrected_bit_score,
                    },
                    approx_id: hsp.approx_id_percent(query.frames[0], subject),
                    identity: hsp.id_percent(),
                });
            }
            capture
                .lock()
                .map_err(|_| io::Error::other("Ranked edge capture lock poisoned."))?
                [query_id as usize] = edges;
            crate::basic::statistics::merge_statistics(&output_statistics);
            return Ok(QueryOutput {
                bytes: Some(Vec::new()),
                aligned,
            });
        }
        if let Some((dictionary_ids, format)) = self.intermediate {
            let mut records = Vec::new();
            for target in matches {
                let block_id = target.target_block_id as usize;
                for hsp in target.hsps {
                    IntermediateRecord::write(
                        &mut records,
                        &hsp,
                        query_id,
                        dictionary_ids[block_id] as i64,
                        self.dictionary.entries[block_id].oid,
                        &format,
                    )?;
                }
            }
            let mut bytes = Vec::new();
            if !records.is_empty() {
                bytes.write_all(&query_id.to_ne_bytes())?;
                let size = u32::try_from(records.len())
                    .map_err(|_| io::Error::other("intermediate query exceeds uint32"))?;
                bytes.write_all(&size.to_ne_bytes())?;
                bytes.extend(records);
            }
            return Ok(QueryOutput {
                bytes: Some(bytes),
                aligned,
            });
        }
        let oid = query.oid.unwrap_or_else(|| {
            self.query_begin
                + self
                    .query_oids
                    .get(query_id as usize)
                    .copied()
                    .unwrap_or(query_id as u64)
        });
        let bytes = if let Some(format) = self.format {
            format
                .query(
                    query_id,
                    oid,
                    query.title,
                    query.source.len() as i32,
                    &query.frames,
                    matches,
                    self.dictionary,
                    self.matrix,
                    contexts == 6,
                    self.database_letters,
                    self.config.report_unaligned,
                    &self.config.sam,
                    &self.config.tuning.output,
                    query.self_aln_score,
                )?
                .0
        } else {
            let mut bytes = Vec::new();
            crate::commands::search_output::write_tabular_query(
                &mut bytes,
                self.tabular.unwrap(),
                query_id,
                oid,
                query.title,
                query.source,
                &query.frames,
                query.mate,
                query.quality,
                query.self_aln_score,
                matches,
                self.dictionary,
            )?;
            bytes
        };
        crate::basic::statistics::merge_statistics(&output_statistics);
        Ok(QueryOutput {
            bytes: Some(bytes),
            aligned,
        })
    }
}

impl NativeRankingBackend<'_> {
    fn worker(&self) -> Worker<'_> {
        Worker {
            edge_capture: self.edge_capture,
            aligned_targets: self.aligned_targets,
            database_letters: self.physical_totals.1,
            config: self.config,
            queries: self.queries,
            matrix: self.matrix,
            targets: &self.targets,
            dictionary: &self.dictionary,
            format: self.format,
            tabular: self.tabular.as_ref(),
            query_begin: self.query_begin,
            query_oids: self.query_oids,
            intermediate: self
                .intermediate_format
                .map(|format| (self.intermediate_ids.as_slice(), format)),
        }
    }

    fn self_scores(&mut self) {
        if self
            .config
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
            || self
                .config
                .outfmt
                .iter()
                .any(|field| field == b"normalized_bitscore")
        {
            let scores = self.pool.install(|| {
                crate::commands::search_output::self_alignment_scores_with_window(
                    self.targets.seqs(),
                    self.matrix,
                    self.config.comp_based_stats,
                    self.config.tuning.cbs_window,
                )
            });
            for (entry, score) in self.dictionary.entries.iter_mut().zip(scores) {
                entry.self_aln_score = Some(score);
            }
            if let Some(dictionary) = self.persistent_dictionary.as_mut() {
                for entry in &self.dictionary.entries {
                    if let Some(&id) = dictionary.by_oid.get(&entry.oid) {
                        dictionary.entries[id].self_aln_score = entry.self_aln_score;
                    }
                }
            }
        }
    }
}

impl GlobalRankingBackend for NativeRankingBackend<'_> {
    fn database_sequence_count(&self) -> usize {
        self.database_totals.0 as usize
    }

    fn load_target_sequences(
        &mut self,
        filter: &BitVector,
        load: SequenceLoadConfig,
    ) -> Result<Vec<u64>, String> {
        if !load.reset_seqinfo || !load.sequences || self.random_access.is_some() {
            return Err("Invalid ranked database loading lifecycle".into());
        }
        self.references.rewind();
        let mut oids = Vec::new();
        for oid in 0..self.database_totals.0 {
            if !filter.get(oid as usize) {
                continue;
            }
            let mut record = self
                .references
                .record(oid)
                .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
            crate::masking::remove_bit_mask(&mut record.sequence);
            self.dictionary.insert(oid, &record, None, false);
            self.targets.push_back_bytes(
                &record.sequence,
                Some(&record.id),
                None,
                oid,
                SequenceType::AminoAcid,
                0,
                false,
            )?;
            oids.push(oid);
        }
        if self.references.checks_protein_input() {
            crate::data::sequence_file::check_protein_input(
                (0..self.targets.seqs().len()).map(|id| self.targets.seqs().get(id)),
                self.config.tuning.ignore_warnings,
            )?;
        }
        self.references
            .load_matched_taxonomy(&mut self.dictionary)
            .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
        if self.config.masking == MaskingMode::None {
            self.self_scores();
        }
        if let Some(dictionary) = self.persistent_dictionary.as_mut() {
            for entry in &self.dictionary.entries {
                let record = crate::data::fasta::ByteFastaRecord {
                    id: entry.title.clone(),
                    sequence: self
                        .targets
                        .seqs()
                        .get(self.intermediate_ids.len())
                        .to_vec(),
                    quality: None,
                };
                let id = dictionary.insert(entry.oid, &record, entry.self_aln_score, false);
                dictionary.entries[id].taxids = entry.taxids.clone();
                self.intermediate_ids.push(id as u32);
            }
            self.references
                .load_matched_taxonomy(dictionary)
                .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
        }
        Ok(oids)
    }

    fn mask_target_sequences(&mut self) -> Result<(), String> {
        let tantan = (self.config.masking == MaskingMode::Tantan)
            .then(|| crate::commands::blastp::native_tantan_masker(self.config, self.matrix));
        for id in 0..self.targets.seqs().len() {
            let sequence = self.targets.seqs_mut().get_mut(id);
            match self.config.masking {
                MaskingMode::Tantan => {
                    tantan.as_ref().unwrap().mask(sequence);
                    for letter in sequence {
                        if *letter & SEED_MASK != 0 {
                            *letter = MASK_LETTER;
                        }
                    }
                }
                MaskingMode::BlastSeg => {
                    crate::masking::mask_sequence(sequence, crate::masking::MaskingAlgo::Seg)
                }
                MaskingMode::None => (),
            }
        }
        self.self_scores();
        Ok(())
    }

    fn init_random_access(&mut self, query: BlockId, reference: BlockId) -> Result<(), String> {
        if self.random_access.is_some() {
            return Err("Ranked random access already active".into());
        }
        self.references.init_seq_access();
        // Native OID reads use independently owned buffered handles rather than
        // the source's global random-access dictionaries. This state guards
        // the corresponding output lifetime; metadata is already materialized.
        self.random_access = Some((query, reference));
        Ok(())
    }

    fn end_random_access(&mut self) -> Result<(), String> {
        self.random_access
            .take()
            .ok_or_else(|| "Ranked random access was not active".to_owned())?;
        Ok(())
    }

    fn align_query(
        &mut self,
        query: BlockId,
        targets: &[RankedTarget],
        iterated: bool,
    ) -> Result<QueryOutput, String> {
        if iterated != self.intermediate_format.is_some() || self.random_access.is_none() {
            return Err("Invalid ranked alignment lifecycle".into());
        }
        self.worker()
            .align(query, targets)
            .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())
    }

    fn align_batch(
        &mut self,
        queries: &[(BlockId, &[RankedTarget])],
        iterated: bool,
    ) -> Result<Option<Vec<QueryOutput>>, String> {
        if iterated != self.intermediate_format.is_some() || self.random_access.is_none() {
            return Err("Invalid ranked alignment lifecycle".into());
        }
        let worker = self.worker();
        self.pool
            .install(|| {
                queries
                    .par_iter()
                    .map(|&(query, targets)| {
                        worker
                            .align(query, targets)
                            .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .map(Some)
    }

    fn align_batch_ordered(
        &mut self,
        queries: &[(BlockId, &[RankedTarget])],
        skipped: &[(BlockId, Option<Vec<u8>>)],
        iterated: bool,
    ) -> Result<Option<Vec<bool>>, String> {
        if iterated != self.intermediate_format.is_some() || self.random_access.is_none() {
            return Err("Invalid ranked alignment lifecycle".into());
        }
        let Some(&(last, _)) = queries.last() else {
            return Ok(Some(Vec::new()));
        };
        let begin = self.ordered_next;
        let mut end = last.checked_add(1).ok_or("Ranked query index overflow")?;
        while skipped
            .binary_search_by_key(&end, |(query, _)| *query)
            .is_ok()
        {
            end = end.checked_add(1).ok_or("Ranked query index overflow")?;
        }
        let worker = Worker {
            edge_capture: self.edge_capture,
            aligned_targets: self.aligned_targets,
            database_letters: self.physical_totals.1,
            config: self.config,
            queries: self.queries,
            matrix: self.matrix,
            targets: &self.targets,
            dictionary: &self.dictionary,
            format: self.format,
            tabular: self.tabular.as_ref(),
            query_begin: self.query_begin,
            query_oids: self.query_oids,
            intermediate: self
                .intermediate_format
                .map(|format| (self.intermediate_ids.as_slice(), format)),
        };
        let heartbeat = self.config.tuning.heartbeat
            && !self.config.swipe_all
            && crate::util::log_stream::log_stream()
                .lock()
                .unwrap()
                .to_cout();
        let titles = (begin..end)
            .map(|query| self.queries.query(query as usize).map(|query| query.title))
            .collect::<io::Result<Vec<_>>>()
            .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
        let mut aligned = vec![false; queries.len()];
        let intermediate = &mut self.intermediate_output;
        let writer = &mut self.writer;
        let sink = crate::output::live_sink::LiveSink::new(
            begin as usize,
            heartbeat,
            |(index, output): (Option<usize>, QueryOutput)| {
                if let Some(index) = index {
                    aligned[index] = output.aligned;
                }
                if let Some(bytes) = output.bytes {
                    if iterated {
                        intermediate
                            .as_mut()
                            .ok_or_else(|| {
                                io::Error::other("Ranked intermediate stream is absent")
                            })?
                            .write_all(&bytes)?;
                    } else {
                        writer.write_ranked_bytes(&bytes)?;
                    }
                }
                Ok(())
            },
        );
        let first_skipped = skipped.partition_point(|(query, _)| *query < begin);
        for (query, bytes) in skipped[first_skipped..]
            .iter()
            .take_while(|(query, _)| *query < end)
        {
            let size = bytes.as_ref().map_or(0, Vec::capacity);
            sink.push(
                *query as usize,
                (
                    None,
                    QueryOutput {
                        bytes: bytes.clone(),
                        aligned: false,
                    },
                ),
                size,
            )
            .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
        }
        sink.execute_bytes(begin as usize, &titles, || {
            self.pool.install(|| {
                queries
                    .par_iter()
                    .enumerate()
                    .try_for_each(|(index, &(query, targets))| {
                        let _progress_scope = sink.task_scope();
                        let output = worker.align(query, targets)?;
                        let size = output.bytes.as_ref().map_or(0, Vec::capacity);
                        sink.push(query as usize, (Some(index), output), size)
                    })
            })
        })
        .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
        sink.finish(end as usize)
            .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
        self.ordered_next = end;
        Ok(Some(aligned))
    }

    fn write_output(&mut self, bytes: &[u8]) -> Result<(), String> {
        if self.random_access.is_none() {
            return Err("Ranked output outside random-access lifetime".into());
        }
        if self.intermediate_format.is_some() {
            self.intermediate_output
                .as_mut()
                .ok_or_else(|| "Ranked intermediate stream is absent".to_owned())?
                .write_all(bytes)
                .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())
        } else {
            self.writer
                .write_ranked_bytes(bytes)
                .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())
        }
    }

    fn clear_target_sequences(&mut self) {
        self.targets = Block::new();
        self.dictionary.entries.clear();
        self.dictionary.by_oid.clear();
    }

    fn init_dictionary_block(&mut self, query: BlockId, reference: BlockId) -> Result<(), String> {
        if self.persistent_dictionary.is_none() || self.intermediate_format.is_none() {
            return Err("Iterated ranked dictionary owner is absent".into());
        }
        self.init_random_access(query, reference)
    }
    fn close_dictionary_block(&mut self) -> Result<(), String> {
        self.end_random_access()
    }
    fn delete_merged_query_list(&mut self) -> Result<(), String> {
        // The active Source driver uses the persistent ranking table. Its
        // alternative merged-list writer and extend call in output/join_blocks.cpp
        // are commented out; this backend therefore never owns such a file.
        Err("Table-based global ranking does not own a merged query list".into())
    }
    fn finish_intermediate_file(&mut self) -> Result<(), String> {
        if self.intermediate_format.is_none()
            || self.random_access.is_some()
            || self.intermediate_finished
        {
            return Err("Invalid ranked intermediate finish lifecycle".into());
        }
        let output = self
            .intermediate_output
            .as_mut()
            .ok_or_else(|| "Ranked intermediate stream is absent".to_owned())?;
        IntermediateRecord::finish_file(output)
            .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
        output
            .flush()
            .map_err(|error| self.io_errors.capture::<()>(Err(error)).unwrap_err())?;
        self.intermediate_finished = true;
        Ok(())
    }
}

pub(super) fn extend(
    config: &BlastpConfig,
    references: &mut dyn RankingReferences,
    queries: &QueryChunk,
    writer: &mut NativeSearchWriter,
    matrix: &ScoreMatrix,
    table: &[Hit],
    query_begin: u64,
    query_oids: &[u64],
    dictionary: Option<&mut ReferenceDictionary>,
    aligned: &mut Vec<bool>,
    aligned_targets: Option<&std::sync::Mutex<crate::util::data_structures::BitVector>>,
) -> io::Result<Option<super::IntermediateFile>> {
    extend_impl(
        config,
        references,
        queries,
        writer,
        matrix,
        table,
        query_begin,
        query_oids,
        dictionary,
        aligned,
        aligned_targets,
        None,
    )
}

pub(super) fn extend_captured(
    config: &BlastpConfig,
    references: &mut dyn RankingReferences,
    queries: &QueryChunk,
    writer: &mut NativeSearchWriter,
    matrix: &ScoreMatrix,
    table: &[Hit],
    query_begin: u64,
    query_oids: &[u64],
    aligned: &mut Vec<bool>,
    aligned_targets: Option<&std::sync::Mutex<crate::util::data_structures::BitVector>>,
) -> io::Result<Vec<crate::commands::blastp::InMemorySearchEdge>> {
    let capture = std::sync::Mutex::new(vec![Vec::new(); queries.count()]);
    extend_impl(
        config,
        references,
        queries,
        writer,
        matrix,
        table,
        query_begin,
        query_oids,
        None,
        aligned,
        aligned_targets,
        Some(&capture),
    )?;
    Ok(capture
        .into_inner()
        .map_err(|_| io::Error::other("Ranked edge capture lock poisoned."))?
        .into_iter()
        .flatten()
        .collect())
}

fn extend_impl(
    config: &BlastpConfig,
    references: &mut dyn RankingReferences,
    queries: &QueryChunk,
    writer: &mut NativeSearchWriter,
    matrix: &ScoreMatrix,
    table: &[Hit],
    query_begin: u64,
    query_oids: &[u64],
    dictionary: Option<&mut ReferenceDictionary>,
    aligned: &mut Vec<bool>,
    aligned_targets: Option<&std::sync::Mutex<crate::util::data_structures::BitVector>>,
    edge_capture: Option<&std::sync::Mutex<Vec<Vec<crate::commands::blastp::InMemorySearchEdge>>>>,
) -> io::Result<Option<super::IntermediateFile>> {
    let normalized = crate::commands::blastp::normalize_clustering_filter(config)?;
    let config = normalized.as_ref();
    let threads = if config.tuning.threads_align != 0 {
        config.tuning.threads_align as usize
    } else if config.threads > 0 {
        config.threads as usize
    } else {
        rayon::current_num_threads().max(1)
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;
    let format = SearchFormat::from_config(config);
    let tabular = if format.is_none() {
        Some(crate::commands::search_output::tabular_format(config)?)
    } else {
        None
    };
    let iterated = dictionary.is_some();
    let intermediate_format = iterated.then(|| writer.intermediate_format());
    let intermediate = if iterated {
        Some(super::IntermediateFile::new(config)?)
    } else {
        None
    };
    let (file, output) =
        intermediate.map_or((None, None), |(file, output)| (Some(file), Some(output)));
    let io_errors = crate::util::io::IoErrorOwner::default();
    let mut backend = NativeRankingBackend {
        io_errors: io_errors.clone(),
        database_totals: references.totals(),
        physical_totals: references.metadata_totals(),
        edge_capture,
        ordered_next: 0,
        aligned_targets,
        dictionary: ReferenceDictionary::new(
            references.metadata_totals().0,
            writer.needs_reference_sequences(),
        ),
        config,
        references,
        queries,
        writer,
        matrix,
        query_begin,
        query_oids,
        targets: Block::new(),
        format,
        tabular,
        pool,
        random_access: None,
        persistent_dictionary: dictionary,
        intermediate_ids: Vec::new(),
        intermediate_format,
        intermediate_output: output,
        intermediate_finished: false,
    };
    let contexts = if config.translated_query_layout.is_some() {
        6
    } else {
        1
    };
    let mut extension = ExtendConfig {
        current_query_block: config.tuning.search_block.query,
        current_reference_block: config.tuning.search_block.reference,
        threads,
        target_masking: config.masking != MaskingMode::None,
        full_titles: true,
        query_contexts: contexts,
        query_sequence_count: queries.count() * contexts,
        global_ranking_targets: config.global_ranking_targets as usize,
        ranking_table: table.to_vec(),
        iterated,
        track_aligned_queries: iterated || config.auxiliary_output.enabled(),
        query_aligned: aligned.clone(),
        ..Default::default()
    };
    io_errors.finish(crate::align::global_ranking::extend(
        &mut backend,
        &mut extension,
    ))?;
    *aligned = extension.query_aligned;
    drop(backend);
    Ok(file)
}
