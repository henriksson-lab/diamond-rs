//! Concrete adapters for the legacy translated-query extension pipeline.
//!
//! The source uses global Blocks and score state. This adapter borrows their
//! native owners and returns normal `Match` records at the same output/culling
//! boundary, so command-wide writers retain headers, dictionaries and streams.

use std::collections::HashMap;
use std::io;
use std::sync::Arc;

use super::banded_swipe_pipeline::{BandedSwipeBackend, Pipeline, PipelineConfig, RunConfig};
use super::query_mapper::{
    GenerateOutputConfig, MatchContext, OutputMode, QueryIntro, QueryMapper, QueryMapperBackend,
    QueryMapperConfig, QueryMapperOutput,
};
use crate::align::hsp::{Hsp, Match};
use crate::basic::statistics::Statistics;
use crate::basic::translate::Strand;
use crate::basic::value::{BlockId, Letter, OId};
use crate::data::block::Block;
use crate::data::sequence_set::SequenceSet;
use crate::dp::banded_3frame::{banded_3frame_swipe_strand_with_controls_and_bins, DpStat};
use crate::dp::swipe::DpTarget;
use crate::dp::ungapped::DiagonalSegment;
use crate::search::hit::Hit;
use crate::stats::score_matrix::ScoreMatrix;

#[derive(Debug, Clone)]
pub struct NativeConfig {
    pub cbs_matrix_scale: i32,
    pub band_bin: i32,
    pub col_bin: i32,
    pub swipe_chunk_size: u32,
    pub transcript_len_estimate: f64,
    pub mapper: QueryMapperConfig,
    pub pipeline: PipelineConfig,
    pub max_target_seqs: i64,
    pub max_evalue: f64,
    pub min_bit_score: f64,
    pub ungapped_xdrop_bits: f64,
    pub target_parallel: bool,
    pub report_unaligned: bool,
    /// Source null format runs DP but omits final output/culling and aligned tracking.
    pub discard_output: bool,
}

impl Default for NativeConfig {
    fn default() -> Self {
        Self {
            cbs_matrix_scale: 1,
            band_bin: 24,
            col_bin: 400,
            swipe_chunk_size: 256,
            transcript_len_estimate: 1.0,
            mapper: QueryMapperConfig {
                query_contexts: 6,
                query_translated: true,
                ..Default::default()
            },
            pipeline: PipelineConfig {
                frame_shift: true,
                ..Default::default()
            },
            max_target_seqs: 25,
            max_evalue: 0.001,
            min_bit_score: 0.0,
            ungapped_xdrop_bits: 12.3,
            target_parallel: false,
            report_unaligned: false,
            discard_output: false,
        }
    }
}

/// A real Block/score/DP adapter, reusable by CLI and typed native searches.
/// Query contexts must already be translated, masked and restored for extension.
pub struct NativeBackend<'a> {
    queries: &'a Block,
    query_sequences: &'a SequenceSet,
    targets: &'a Block,
    matrix: &'a ScoreMatrix,
    config: NativeConfig,
    query_titles: Vec<&'a [u8]>,
    target_titles: Vec<&'a [u8]>,
    target_rank_ids: Option<&'a [Vec<u32>]>,
    pool: rayon::ThreadPool,
}

impl<'a> NativeBackend<'a> {
    pub fn new(
        queries: &'a Block,
        targets: &'a Block,
        matrix: &'a ScoreMatrix,
        config: NativeConfig,
        target_rank_ids: Option<&'a [Vec<u32>]>,
    ) -> io::Result<Self> {
        Self::new_with_query_sequences(
            queries,
            queries.seqs(),
            targets,
            matrix,
            config,
            target_rank_ids,
        )
    }

    /// Borrow the restored arena while its metadata Block temporarily has its
    /// sequences taken by the native seed pipeline; no translation or copy.
    pub fn new_with_query_sequences(
        queries: &'a Block,
        query_sequences: &'a SequenceSet,
        targets: &'a Block,
        matrix: &'a ScoreMatrix,
        mut config: NativeConfig,
        target_rank_ids: Option<&'a [Vec<u32>]>,
    ) -> io::Result<Self> {
        if config.mapper.query_contexts != 6
            || !config.mapper.query_translated
            || !config.pipeline.frame_shift
            || matrix.frame_shift() == 0
        {
            return Err(io::Error::other("Frameshift extension requires six translated contexts and a nonzero matrix penalty."));
        }
        let query_count = queries.source_seqs().len();
        if query_sequences.len() != query_count.saturating_mul(6)
            || queries.oid_count() as usize != query_count
            || targets.oid_count() as usize != targets.seqs().len()
        {
            return Err(io::Error::other("Invalid frameshift Block layout."));
        }
        if config.max_target_seqs <= 0 {
            return Err(io::Error::other(
                "Frameshift target limit must be positive.",
            ));
        }
        if (config.mapper.taxon_k != 0 && target_rank_ids.is_none())
            || target_rank_ids.is_some_and(|ids| ids.len() != targets.seqs().len())
        {
            return Err(io::Error::other(
                "Frameshift taxonomy ranks are absent or have the wrong target count.",
            ));
        }
        let titles = |block: &'a Block, count: usize| -> io::Result<Vec<&'a [u8]>> {
            if count == 0 {
                return Ok(Vec::new());
            }
            let ids = block.ids().map_err(io::Error::other)?;
            if ids.size() as usize != count {
                return Err(io::Error::other(
                    "Frameshift Block titles have the wrong sequence count.",
                ));
            }
            Ok((0..count).map(|i| ids.get(i)).collect())
        };
        let query_titles = titles(queries, query_count)?;
        let target_titles = titles(targets, targets.seqs().len())?;
        let threads = if config.pipeline.threads == 0 {
            rayon::current_num_threads().max(1)
        } else {
            config.pipeline.threads
        };
        config.pipeline.threads = threads;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(io::Error::other)?;
        Ok(Self {
            queries,
            query_sequences,
            targets,
            matrix,
            config,
            query_titles,
            target_titles,
            target_rank_ids,
            pool,
        })
    }

    /// Run Source's target loading/ranking, optional score-only culling, exact
    /// three-frame traceback, inner culling and final report filtering.
    pub fn extend_query(
        &self,
        query_id: BlockId,
        hits: Vec<Hit>,
        statistics: &mut Statistics,
        dp_stat: &mut DpStat,
    ) -> io::Result<NativeQueryResult> {
        self.extend_query_impl(query_id, hits, statistics, dp_stat, true)
    }

    /// Collect ordinary records for the command's final writer. The mapper's
    /// temporary output must not count final matches before all blocks/rounds
    /// have been joined; DP/extension counters still belong to this worker.
    pub(crate) fn extend_query_for_capture(
        &self,
        query_id: BlockId,
        hits: Vec<Hit>,
        statistics: &mut Statistics,
        dp_stat: &mut DpStat,
    ) -> io::Result<NativeQueryResult> {
        self.extend_query_impl(query_id, hits, statistics, dp_stat, false)
    }

    fn extend_query_impl(
        &self,
        query_id: BlockId,
        hits: Vec<Hit>,
        statistics: &mut Statistics,
        dp_stat: &mut DpStat,
        count_output: bool,
    ) -> io::Result<NativeQueryResult> {
        if query_id as usize >= self.queries.source_seqs().len()
            || hits.iter().any(|hit| hit.query / 6 != query_id)
        {
            return Err(io::Error::other(
                "Frameshift hits do not belong to the source query.",
            ));
        }
        for hit in &hits {
            let position = usize::try_from(hit.subject).map_err(|_| {
                io::Error::other("Frameshift seed position exceeds native address space.")
            })?;
            if position < self.targets.seqs().offsets()[0]
                || position >= self.targets.seqs().raw_len()
            {
                return Err(io::Error::other(
                    "Frameshift seed position is outside its sequence.",
                ));
            }
            let (subject, position) = self.targets.seqs().local_position(hit.subject as usize);
            let frame = hit.query % 6;
            if subject >= self.targets.seqs().len()
                || position >= self.targets.seqs().get(subject).len()
                || hit.seed_offset as usize >= self.query_seq(query_id, frame).len()
            {
                return Err(io::Error::other(
                    "Frameshift seed position is outside its sequence.",
                ));
            }
        }
        let mut mapper =
            QueryMapper::new(query_id as usize, hits, self, self.config.mapper.clone());
        mapper.target_parallel = self.config.target_parallel;
        mapper.init();
        let mut pipeline = Pipeline::new(mapper, dp_stat, self.config.pipeline);
        pipeline.run(
            statistics,
            RunConfig {
                max_target_seqs: self.config.max_target_seqs,
            },
        );
        if self.config.discard_output {
            return Ok(NativeQueryResult {
                matches: Vec::new(),
                aligned: false,
                report_unaligned: false,
            });
        }
        let mut output = CollectedOutput::new(self.targets, self.config.report_unaligned);
        let mut output_statistics = Statistics::new();
        let aligned = pipeline.mapper.generate_output(
            &mut output,
            if count_output {
                statistics
            } else {
                &mut output_statistics
            },
            GenerateOutputConfig {
                max_target_seqs: self.config.max_target_seqs,
                blocked_processing: false,
                current_ref_block: 0,
            },
        );
        Ok(NativeQueryResult {
            matches: output.matches,
            aligned,
            report_unaligned: output.unaligned,
        })
    }
}

impl QueryMapperBackend<[u8]> for NativeBackend<'_> {
    fn source_query_len(&self, query: u32) -> u32 {
        self.query_sequences
            .reverse_translated_len(query as usize * 6) as u32
    }
    fn query_seq(&self, query: u32, frame: u32) -> &[Letter] {
        self.query_sequences
            .get(query as usize * 6 + frame as usize)
    }
    fn query_source_seq(&self, query: u32) -> &[Letter] {
        self.queries.source_seqs().get(query as usize)
    }
    fn query_title(&self, query: u32) -> &[u8] {
        crate::util::sequence::c_string_bytes(self.query_titles[query as usize])
    }
    fn query_oid(&self, query: u32) -> OId {
        self.queries.block_id2oid(query)
    }
    fn local_target_position(&self, position: u64) -> (usize, usize) {
        self.targets.seqs().local_position(position as usize)
    }
    fn target_seq(&self, subject: BlockId) -> &[Letter] {
        self.targets.seqs().get(subject as usize)
    }
    fn target_oid(&self, subject: BlockId) -> OId {
        self.targets.block_id2oid(subject)
    }
    fn target_taxon_rank_ids(&self, oid: OId) -> Vec<u32> {
        let id = self.targets.oid2block_id(oid).expect("loaded target OID");
        self.target_rank_ids
            .map_or_else(Vec::new, |ranks| ranks[id as usize].clone())
    }
    fn target_title(&self, subject: BlockId, _: OId) -> Vec<u8> {
        crate::util::sequence::c_string_bytes(self.target_titles[subject as usize]).to_vec()
    }
    fn target_dict_id(&self, _: usize, subject: BlockId) -> usize {
        subject as usize
    }
    fn score_matrix(&self) -> &ScoreMatrix {
        self.matrix
    }
    fn xdrop_ungapped(
        &self,
        query: &[Letter],
        subject: &[Letter],
        qa: usize,
        sa: usize,
    ) -> DiagonalSegment {
        crate::dp::ungapped::xdrop_ungapped(
            query,
            subject,
            qa,
            sa,
            self.matrix.rawscore_int(self.config.ungapped_xdrop_bits),
            self.matrix,
        )
    }
    fn report_cutoff(&self, score: i32, evalue: f64) -> bool {
        self.matrix.report_cutoff(
            score,
            evalue,
            self.config.min_bit_score,
            self.config.max_evalue,
        )
    }
    fn bitscore(&self, score: i32) -> f64 {
        self.matrix.bitscore(score as f64)
    }
}

impl BandedSwipeBackend<[u8]> for NativeBackend<'_> {
    fn banded_3frame_swipe(
        &self,
        query: u32,
        strand: Strand,
        targets: &mut [DpTarget],
        stat: &mut DpStat,
        score_only: bool,
        target_parallel: bool,
    ) -> Vec<Hsp> {
        let offset = if strand == Strand::Forward { 0 } else { 3 };
        let frames = [
            self.query_seq(query, offset),
            self.query_seq(query, offset + 1),
            self.query_seq(query, offset + 2),
        ];
        let mut hsps = self.pool.install(|| {
            banded_3frame_swipe_strand_with_controls_and_bins(
                frames,
                strand,
                self.query_source_seq(query).len() as i32,
                targets,
                self.matrix,
                stat,
                score_only,
                target_parallel,
                self.config.swipe_chunk_size,
                self.config.transcript_len_estimate * self.config.cbs_matrix_scale as f64,
                self.config.band_bin,
                self.config.col_bin,
            )
        });
        hsps.retain_mut(|hsp| {
            hsp.score = hsp.score.saturating_mul(self.config.cbs_matrix_scale);
            hsp.bit_score = self.matrix.bitscore(hsp.score as f64);
            hsp.evalue = self.matrix.evalue(
                hsp.score,
                frames[0].len() as u32,
                hsp.target_seq.len() as u32,
            );
            self.report_cutoff(hsp.score, hsp.evalue)
        });
        hsps
    }
}

#[derive(Debug)]
pub struct NativeQueryResult {
    pub matches: Vec<Match>,
    pub aligned: bool,
    pub report_unaligned: bool,
}

/// Collect the exact translated output boundary into the native writer's
/// ordinary records. The OID map preserves loaded-block IDs across culling.
struct CollectedOutput {
    matches: Vec<Match>,
    target_ids: HashMap<OId, BlockId>,
    report_unaligned: bool,
    unaligned: bool,
    active_query: Option<(BlockId, usize)>,
}
impl CollectedOutput {
    fn new(targets: &Block, report_unaligned: bool) -> Self {
        Self {
            matches: Vec::new(),
            target_ids: (0..targets.seqs().len())
                .map(|id| (targets.block_id2oid(id as BlockId), id as BlockId))
                .collect(),
            report_unaligned,
            unaligned: false,
            active_query: None,
        }
    }
}
impl QueryMapperOutput<[u8]> for CollectedOutput {
    fn mode(&self) -> OutputMode {
        OutputMode::Regular
    }
    fn report_unaligned(&self) -> bool {
        self.report_unaligned
    }
    fn begin_query(&mut self, intro: QueryIntro<'_, [u8]>) -> usize {
        assert!(self.active_query.is_none());
        let begin = self.matches.len();
        self.active_query = Some((intro.query_id, begin));
        self.unaligned = intro.unaligned;
        begin
    }
    fn write_match(&mut self, context: MatchContext<'_, [u8]>) {
        assert_eq!(
            self.active_query.map(|state| state.0),
            Some(context.query_id)
        );
        if context.hsp_index == 0 {
            let block_id = self.target_ids[&context.database_id];
            let mut target = Match::new(block_id, context.database_id);
            target.seq = Arc::from(context.subject);
            target.filter_score = context.hsp.score;
            target.filter_evalue = context.hsp.evalue;
            self.matches.push(target);
        }
        let target = self.matches.last_mut().expect("query target started");
        assert_eq!(target.target_oid, context.database_id);
        target.hsps.push(context.hsp.clone());
    }
    fn finish_query(&mut self, begin: usize) {
        let (_, expected) = self.active_query.take().expect("query output started");
        assert_eq!(begin, expected);
    }
}

#[cfg(test)]
mod byte_title_tests {
    use super::super::query_mapper::Target;
    use super::*;
    use crate::basic::value::SequenceType;
    use crate::util::interval::Interval;

    #[test]
    fn actual_frameshift_owner_forwards_byte_titles_through_mapper_output() {
        let mut queries = Block::new();
        queries
            .push_back_translated_with_genetic_code_bytes(
                &[0, 3, 2, 0, 3, 2, 0, 3, 2, 0, 3, 2],
                Some(b"q\xff\0stored"),
                0,
                0b111111,
                1,
                1,
            )
            .unwrap();
        let mut targets = Block::new();
        targets
            .push_back_bytes(
                &[0, 1, 2, 3],
                Some(b"s\xfe\0stored"),
                None,
                7,
                SequenceType::AminoAcid,
                1,
                false,
            )
            .unwrap();
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 15, 1, 0).unwrap();
        let mut config = NativeConfig::default();
        config.pipeline.threads = 1;
        let backend = NativeBackend::new(&queries, &targets, &matrix, config, None).unwrap();
        assert_eq!(queries.ids().unwrap().get(0), b"q\xff\0stored");
        assert_eq!(backend.query_title(0), b"q\xff");
        assert_eq!(backend.target_title(0, 7), b"s\xfe");
        let mut mapper = QueryMapper::new(0, vec![], &backend, backend.config.mapper.clone());
        let mut target = Target::new(0, 0, vec![0, 1, 2, 3], vec![]);
        target.filter_score = 100;
        target.filter_evalue = 1e-20;
        let mut hsp = Hsp::default();
        hsp.score = 100;
        hsp.evalue = 1e-20;
        hsp.length = 4;
        hsp.identities = 4;
        hsp.query_range = Interval::new(0, 4);
        hsp.query_source_range = Interval::new(0, 12);
        hsp.subject_range = Interval::new(0, 4);
        target.hsps.push(hsp);
        mapper.targets.push(target);
        #[derive(Default)]
        struct Capture {
            query: Vec<u8>,
            subject: Vec<u8>,
        }
        impl QueryMapperOutput<[u8]> for Capture {
            fn mode(&self) -> OutputMode {
                OutputMode::Regular
            }
            fn begin_query(&mut self, intro: QueryIntro<'_, [u8]>) -> usize {
                self.query = intro.title.to_vec();
                0
            }
            fn write_match(&mut self, context: MatchContext<'_, [u8]>) {
                assert_eq!(context.query_title, self.query);
                self.subject = context.target_title.to_vec();
            }
            fn finish_query(&mut self, _: usize) {}
        }
        let mut capture = Capture::default();
        assert!(mapper.generate_output(
            &mut capture,
            &mut Statistics::new(),
            GenerateOutputConfig {
                max_target_seqs: 1,
                blocked_processing: false,
                current_ref_block: 0,
            }
        ));
        assert_eq!(capture.query, b"q\xff");
        assert_eq!(capture.subject, b"s\xfe");
    }
}
