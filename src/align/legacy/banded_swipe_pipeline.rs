//! Legacy banded-SWIPE extension pipeline.
//!
//! Mirrors `align/legacy/banded_swipe_pipeline.cpp`. Database/query ownership
//! lives in [`QueryMapperBackend`]; the old global translated query, score
//! matrix and SWIPE dispatcher are represented by [`BandedSwipeBackend`].

use crate::align::hsp::Hsp;
use crate::align::legacy::query_mapper::{
    QueryMapper, QueryMapperBackend, SeedHit, Target as QueryTarget,
};
use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::translate::{Frame, Strand, TranslatedPosition};
use crate::dp::banded_3frame::DpStat;
use crate::dp::swipe::{Anchor, CarryOver, DpTarget};
use crate::util::interval::{Interval, IntervalPartition};

/// Explicit replacement for the C++ `translated_query`, score globals, and
/// `banded_3frame_swipe` dispatch.
pub trait BandedSwipeBackend: QueryMapperBackend {
    fn log(&self, _message: &str) {}

    fn banded_3frame_swipe(
        &self,
        query_id: u32,
        strand: Strand,
        targets: &mut [DpTarget],
        dp_stat: &mut DpStat,
        score_only: bool,
        target_parallel: bool,
    ) -> Vec<Hsp>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PipelineConfig {
    /// `None` is C++'s blank option; `run` installs the default 32.
    pub padding: Option<i32>,
    /// `None` is C++'s `-1` sentinel and selects 0.4.
    pub rank_ratio: Option<f64>,
    /// `None` is C++'s `-1.0` sentinel and selects 1000.
    pub rank_factor: Option<f64>,
    pub query_range_culling: bool,
    pub frame_shift: bool,
    pub threads: usize,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            padding: None,
            rank_ratio: None,
            rank_factor: None,
            query_range_culling: false,
            frame_shift: false,
            threads: 1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunConfig {
    pub max_target_seqs: i64,
}

/// Methods supplied by the derived C++ `ExtensionPipeline::BandedSwipe::Target`.
/// They operate on the base query-mapper target with all former globals passed
/// explicitly.
impl QueryTarget {
    pub fn ungapped_stage(&mut self, hits: &[SeedHit]) {
        let first = hits.first().expect("target must contain a seed hit");
        self.top_hit = first.clone();
        for hit in &hits[1..] {
            if hit.ungapped.score > self.top_hit.ungapped.score {
                self.top_hit = hit.clone();
            }
        }
        self.filter_score = self.top_hit.ungapped.score;
    }

    pub fn ungapped_query_range(&self, query_dna_len: i32, query_translated: bool) -> Interval {
        let frame = Frame::from_index(self.top_hit.frame as i32);
        let i0 = (self.top_hit.query_pos as i32 - self.top_hit.subject_pos as i32).max(0);
        let i1 = (self.top_hit.query_pos as i32 + self.subject.len() as i32
            - self.top_hit.subject_pos as i32)
            .min(frame.length(query_dna_len));
        TranslatedPosition::absolute_interval(
            TranslatedPosition::new(i0, frame),
            TranslatedPosition::new(i1, frame),
            query_dna_len,
            query_translated,
        )
    }

    pub fn add_strand(
        &self,
        hits: &[SeedHit],
        output: &mut Vec<DpTarget>,
        target_idx: i32,
        padding: i32,
        query_len: i32,
    ) {
        let Some(first) = hits.first() else {
            return;
        };
        let d_min = -(self.subject.len() as i32 - 1);
        let d_max = query_len - 1;
        let mut d0 = (first.diagonal() - padding).max(d_min);
        let mut d1 = (first.diagonal() + padding).min(d_max);
        for hit in &hits[1..] {
            if hit.diagonal() - d1 <= padding {
                d1 = (hit.diagonal() + padding).min(d_max);
            } else {
                output.push(DpTarget::new(
                    self.subject.clone(),
                    self.subject.len() as i32,
                    d0,
                    d1,
                    target_idx as i64,
                    query_len,
                    CarryOver::default(),
                    Anchor::default(),
                ));
                d0 = (hit.diagonal() - padding).max(d_min);
                d1 = (hit.diagonal() + padding).min(d_max);
            }
        }
        output.push(DpTarget::new(
            self.subject.clone(),
            self.subject.len() as i32,
            d0,
            d1,
            target_idx as i64,
            query_len,
            CarryOver::default(),
            Anchor::default(),
        ));
    }

    pub fn add(
        &self,
        hits: &mut [SeedHit],
        forward: &mut Vec<DpTarget>,
        reverse: &mut Vec<DpTarget>,
        target_idx: i32,
        padding: i32,
        query_len: i32,
        target_parallel: bool,
    ) {
        if target_parallel {
            hits.sort_by(SeedHit::compare_diag_strand);
        } else {
            hits.sort_by(SeedHit::compare_diag_strand2);
        }
        let reverse_begin = hits
            .iter()
            .position(|hit| hit.strand() == Strand::Reverse)
            .unwrap_or(hits.len());
        self.add_strand(
            &hits[..reverse_begin],
            forward,
            target_idx,
            padding,
            query_len,
        );
        self.add_strand(
            &hits[reverse_begin..],
            reverse,
            target_idx,
            padding,
            query_len,
        );
    }

    pub fn set_filter_score(&mut self) {
        self.filter_score = 0;
        self.filter_evalue = f64::MAX;
        for hsp in &self.hsps {
            self.filter_score = self.filter_score.max(hsp.score);
            self.filter_evalue = self.filter_evalue.min(hsp.evalue);
        }
    }

    pub fn reset(&mut self) {
        self.hsps.clear();
    }

    pub fn finish(&mut self, source_query_len: i32, query_translated: bool, frame_shift: bool) {
        self.inner_culling();
        if frame_shift {
            return;
        }
        for hsp in &mut self.hsps {
            let frame = Frame::from_index(hsp.frame);
            hsp.query_source_range = TranslatedPosition::absolute_interval(
                TranslatedPosition::new(hsp.query_range.begin, frame),
                TranslatedPosition::new(hsp.query_range.end, frame),
                source_query_len,
                query_translated,
            );
        }
    }

    pub fn is_outranked_partition(
        &self,
        partition: &IntervalPartition,
        source_query_len: i32,
        query_translated: bool,
        rank_ratio: f64,
        toppercent: Option<f64>,
        query_range_cover: f64,
    ) -> bool {
        let range = self.ungapped_query_range(source_query_len, query_translated);
        let covered = if let Some(toppercent) = toppercent {
            let min_score =
                (self.filter_score as f64 / rank_ratio / (1.0 - toppercent / 100.0)) as i32;
            partition.covered_max_score(range, min_score)
        } else {
            let min_score = (self.filter_score as f64 / rank_ratio) as i32;
            partition.covered_min_score(range, min_score)
        };
        covered as f64 / range.length() as f64 * 100.0 >= query_range_cover
    }
}

pub struct Pipeline<'a, 'd, B: BandedSwipeBackend> {
    pub mapper: QueryMapper<'a, B>,
    pub dp_stat: &'d mut DpStat,
    pub config: PipelineConfig,
}

impl<'a, 'd, B: BandedSwipeBackend> Pipeline<'a, 'd, B> {
    pub fn new(
        mapper: QueryMapper<'a, B>,
        dp_stat: &'d mut DpStat,
        config: PipelineConfig,
    ) -> Self {
        Self {
            mapper,
            dp_stat,
            config,
        }
    }

    pub fn target(&mut self, i: usize) -> &mut QueryTarget {
        &mut self.mapper.targets[i]
    }

    pub fn range_ranking(&mut self, max_target_seqs: i64) {
        let rank_ratio = self.config.rank_ratio.unwrap_or(0.4);
        self.mapper.targets.sort_by(QueryTarget::compare_score);
        let mut partition = IntervalPartition::new(max_target_seqs);
        let mut i = 0usize;
        while i < self.mapper.targets.len() {
            if self.mapper.targets[i].is_outranked_partition(
                &partition,
                self.mapper.source_query_len as i32,
                self.mapper.config.query_translated,
                rank_ratio,
                self.mapper.config.toppercent,
                self.mapper.config.culling.query_range_cover,
            ) {
                self.mapper.targets.remove(i);
            } else {
                let range = self.mapper.targets[i].ungapped_query_range(
                    self.mapper.source_query_len as i32,
                    self.mapper.config.query_translated,
                );
                partition.insert(range, self.mapper.targets[i].filter_score);
                i += 1;
            }
        }
    }

    pub fn run_swipe(&mut self, score_only: bool) {
        let mut forward = Vec::new();
        let mut reverse = Vec::new();
        let query_len = self.mapper.query_seq(0).len() as i32;
        let padding = self.config.padding.unwrap_or(32);
        let target_parallel = self.mapper.target_parallel;
        for i in 0..self.mapper.targets.len() {
            let begin = self.mapper.targets[i].begin;
            let end = self.mapper.targets[i].end;
            self.mapper.targets[i].add(
                &mut self.mapper.seed_hits[begin..end],
                &mut forward,
                &mut reverse,
                i as i32,
                padding,
                query_len,
                target_parallel,
            );
        }
        let backend = self.mapper.backend();
        let mut hsps = backend.banded_3frame_swipe(
            self.mapper.query_id,
            Strand::Forward,
            &mut forward,
            self.dp_stat,
            score_only,
            target_parallel,
        );
        hsps.extend(backend.banded_3frame_swipe(
            self.mapper.query_id,
            Strand::Reverse,
            &mut reverse,
            self.dp_stat,
            score_only,
            target_parallel,
        ));
        for hsp in hsps {
            let target = hsp.swipe_target as usize;
            self.mapper.targets[target].hsps.push(hsp);
        }
    }

    /// Deterministic equivalent of the C++ 64-target atomic work allocator.
    /// Per-worker maxima are merged identically after all targets are visited.
    pub fn build_ranking_worker(
        targets: &[QueryTarget],
        workers: usize,
        interval_count: usize,
    ) -> Vec<Vec<i32>> {
        assert!(workers > 0);
        let mut intervals = vec![vec![0; interval_count]; workers];
        for (i, chunk) in targets.chunks(64).enumerate() {
            let worker = i % workers;
            for target in chunk {
                target.add_ranges(&mut intervals[worker]);
            }
        }
        intervals
    }

    pub fn run(&mut self, stat: &mut Statistics, cfg: RunConfig) {
        if self.config.padding.is_none() {
            self.config.padding = Some(32);
        }
        if self.mapper.n_targets() == 0 {
            return;
        }
        stat.inc(StatValue::TargetHits0, self.mapper.n_targets());

        if !self.mapper.target_parallel {
            for target in &mut self.mapper.targets {
                let begin = target.begin;
                let end = target.end;
                target.ungapped_stage(&self.mapper.seed_hits[begin..end]);
            }
            if !self.config.query_range_culling {
                self.mapper.rank_targets(
                    self.config.rank_ratio.unwrap_or(0.4),
                    self.config.rank_factor.unwrap_or(1e3),
                    cfg.max_target_seqs,
                );
            } else {
                self.range_ranking(cfg.max_target_seqs);
            }
        } else {
            self.mapper.backend().log(&format!(
                "Query: {}; Seed hits: {}; Targets: {}",
                self.mapper.query_id,
                self.mapper.seed_hits.len(),
                self.mapper.n_targets()
            ));
        }

        if self.mapper.n_targets() > cfg.max_target_seqs || self.mapper.config.toppercent.is_some()
        {
            stat.inc(StatValue::TargetHits1, self.mapper.n_targets());
            self.run_swipe(true);

            if self.mapper.target_parallel {
                let interval_count =
                    (self.mapper.source_query_len as usize + QueryTarget::INTERVAL as usize - 1)
                        / QueryTarget::INTERVAL as usize;
                let mut intervals = Self::build_ranking_worker(
                    &self.mapper.targets,
                    self.config.threads,
                    interval_count,
                );
                for worker in 1..intervals.len() {
                    for i in 0..interval_count {
                        intervals[0][i] = intervals[0][i].max(intervals[worker][i]);
                    }
                }
                if let Some(toppercent) = self.mapper.config.toppercent {
                    let threshold = 1.0 - toppercent / 100.0;
                    self.mapper
                        .targets
                        .retain(|target| !target.is_outranked(&intervals[0], threshold));
                }
                self.mapper.backend().log(&format!(
                    "Targets after score-only ranking: {}",
                    self.mapper.targets.len()
                ));
            } else {
                for target in &mut self.mapper.targets {
                    target.set_filter_score();
                }
                self.mapper.score_only_culling(cfg.max_target_seqs);
            }
        }

        stat.inc(StatValue::TargetHits2, self.mapper.n_targets());
        for target in &mut self.mapper.targets {
            target.reset();
        }
        self.run_swipe(false);
        for target in &mut self.mapper.targets {
            target.finish(
                self.mapper.source_query_len as i32,
                self.mapper.config.query_translated,
                self.config.frame_shift,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;

    use super::*;
    use crate::align::legacy::query_mapper::QueryMapperConfig;
    use crate::basic::value::{Letter, OId};
    use crate::dp::ungapped::DiagonalSegment;
    use crate::search::hit::Hit;
    use crate::stats::score_matrix::ScoreMatrix;

    struct Backend {
        matrix: ScoreMatrix,
        query: Vec<Vec<Letter>>,
        targets: Vec<Vec<Letter>>,
        calls: RefCell<Vec<(Strand, bool, bool, usize)>>,
    }

    impl Backend {
        fn new() -> Self {
            Self {
                matrix: ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap(),
                query: (0..6).map(|frame| vec![frame as Letter; 20]).collect(),
                targets: vec![vec![0; 20], vec![1; 20], vec![2; 20]],
                calls: RefCell::new(Vec::new()),
            }
        }
    }

    impl QueryMapperBackend for Backend {
        fn source_query_len(&self, _query_id: u32) -> u32 {
            60
        }
        fn query_seq(&self, _query_id: u32, frame: u32) -> &[Letter] {
            &self.query[frame as usize]
        }
        fn query_source_seq(&self, _query_id: u32) -> &[Letter] {
            &self.query[0]
        }
        fn query_title(&self, _query_id: u32) -> &str {
            "query"
        }
        fn query_oid(&self, query_id: u32) -> OId {
            query_id as OId
        }
        fn local_target_position(&self, packed_position: u64) -> (usize, usize) {
            (
                (packed_position / 100) as usize,
                (packed_position % 100) as usize,
            )
        }
        fn target_seq(&self, subject_id: u32) -> &[Letter] {
            &self.targets[subject_id as usize]
        }
        fn target_oid(&self, subject_id: u32) -> OId {
            subject_id as OId
        }
        fn target_taxon_rank_ids(&self, _oid: OId) -> Vec<u32> {
            Vec::new()
        }
        fn target_title(&self, subject_id: u32, _oid: OId) -> String {
            subject_id.to_string()
        }
        fn target_dict_id(&self, _current_ref_block: usize, subject_id: u32) -> usize {
            subject_id as usize
        }
        fn score_matrix(&self) -> &ScoreMatrix {
            &self.matrix
        }
        fn xdrop_ungapped(
            &self,
            _query: &[Letter],
            _subject: &[Letter],
            query_pos: usize,
            subject_pos: usize,
        ) -> DiagonalSegment {
            let subject = subject_pos / 10;
            DiagonalSegment::new(
                query_pos as i32,
                subject_pos as i32,
                4,
                100 - subject as i32 * 10,
            )
        }
        fn report_cutoff(&self, _score: i32, _evalue: f64) -> bool {
            true
        }
        fn bitscore(&self, score: i32) -> f64 {
            score as f64
        }
    }

    impl BandedSwipeBackend for Backend {
        fn banded_3frame_swipe(
            &self,
            _query_id: u32,
            strand: Strand,
            targets: &mut [DpTarget],
            dp_stat: &mut DpStat,
            score_only: bool,
            target_parallel: bool,
        ) -> Vec<Hsp> {
            self.calls
                .borrow_mut()
                .push((strand, score_only, target_parallel, targets.len()));
            dp_stat.gross_cells += targets.len();
            targets
                .iter()
                .map(|target| {
                    let mut hsp = Hsp::new();
                    hsp.swipe_target = target.target_idx as i32;
                    hsp.score = if score_only {
                        100 - target.target_idx as i32 * 50
                    } else {
                        80 - target.target_idx as i32 * 10
                    };
                    hsp.evalue = 1.0 / hsp.score as f64;
                    hsp.frame = if strand == Strand::Forward { 0 } else { 3 };
                    hsp.query_range = Interval::new(1, 6);
                    hsp.query_source_range = Interval::new(0, 30);
                    hsp.subject_range = Interval::new(0, 5);
                    hsp.length = 5;
                    hsp.identities = 5;
                    hsp
                })
                .collect()
        }
    }

    fn mapper_config(toppercent: Option<f64>) -> QueryMapperConfig {
        QueryMapperConfig {
            query_contexts: 6,
            query_translated: true,
            toppercent,
            ..QueryMapperConfig::default()
        }
    }

    #[test]
    fn derived_target_stages_select_top_hit_build_bands_and_finish_coordinates() {
        let hits = vec![
            SeedHit::new(0, 0, 10, 10, DiagonalSegment::new(10, 10, 4, 30)),
            SeedHit::new(0, 0, 12, 13, DiagonalSegment::new(13, 12, 4, 50)),
            SeedHit::new(4, 0, 8, 15, DiagonalSegment::new(15, 8, 4, 40)),
        ];
        let mut target = QueryTarget::new(0, 0, vec![0; 20], vec![]);
        target.ungapped_stage(&hits);
        assert_eq!(target.top_hit.ungapped.score, 50);
        assert_eq!(target.filter_score, 50);

        let mut sorted = hits;
        let mut forward = Vec::new();
        let mut reverse = Vec::new();
        target.add(&mut sorted, &mut forward, &mut reverse, 7, 2, 20, false);
        assert_eq!((forward.len(), reverse.len()), (1, 1));
        assert_eq!(forward[0].target_idx, 7);
        assert_eq!((forward[0].d_begin, forward[0].d_end), (-2, 3));

        target.hsps = vec![{
            let mut h = Hsp::new();
            h.score = 60;
            h.evalue = 0.1;
            h.frame = 0;
            h.query_range = Interval::new(2, 5);
            h.query_source_range = Interval::new(0, 0);
            h
        }];
        target.set_filter_score();
        assert_eq!((target.filter_score, target.filter_evalue), (60, 0.1));
        target.finish(60, true, false);
        assert_eq!(target.hsps[0].query_source_range, Interval::new(6, 15));
        target.reset();
        assert!(target.hsps.is_empty());
    }

    #[test]
    fn range_ranking_uses_partition_capacity_min_score_and_stable_score_order() {
        let backend = Backend::new();
        let mut mapper = QueryMapper::new(0, vec![], &backend, mapper_config(None));
        mapper.config.culling.query_range_cover = 100.0;
        let mut stat = DpStat::default();
        let mut pipeline = Pipeline::new(
            mapper,
            &mut stat,
            PipelineConfig {
                rank_ratio: Some(1.0),
                query_range_culling: true,
                ..PipelineConfig::default()
            },
        );
        for (subject, score) in [(0, 100), (1, 90)] {
            let mut target = QueryTarget::new(0, subject, vec![0; 20], vec![]);
            target.top_hit = SeedHit::new(0, subject, 0, 0, DiagonalSegment::new(0, 0, 4, score));
            target.filter_score = score;
            pipeline.mapper.targets.push(target);
        }
        pipeline.range_ranking(1);
        assert_eq!(pipeline.mapper.targets.len(), 1);
        assert_eq!(pipeline.mapper.targets[0].filter_score, 100);
    }

    #[test]
    fn nonparallel_run_performs_score_only_then_traceback_and_updates_stats() {
        let backend = Backend::new();
        let hits = vec![Hit::new(0, 0, 5), Hit::new(0, 110, 5)];
        let mut mapper = QueryMapper::new(0, hits, &backend, mapper_config(None));
        mapper.init();
        let mut dp_stat = DpStat::default();
        let mut pipeline = Pipeline::new(mapper, &mut dp_stat, PipelineConfig::default());
        let mut statistics = Statistics::default();
        pipeline.run(&mut statistics, RunConfig { max_target_seqs: 1 });
        assert_eq!(pipeline.config.padding, Some(32));
        assert_eq!(pipeline.mapper.targets.len(), 1);
        assert_eq!(pipeline.mapper.targets[0].hsps.len(), 1);
        assert_eq!(statistics.get(StatValue::TargetHits0), 2);
        assert_eq!(statistics.get(StatValue::TargetHits1), 2);
        assert_eq!(statistics.get(StatValue::TargetHits2), 1);
        let calls = backend.calls.borrow();
        assert_eq!(calls.len(), 4);
        assert!(calls[0].1 && calls[1].1);
        assert!(!calls[2].1 && !calls[3].1);
    }

    #[test]
    fn target_parallel_toppercent_removes_only_targets_outranked_in_every_bin() {
        let backend = Backend::new();
        let hits = vec![Hit::new(0, 0, 5), Hit::new(0, 100, 5)];
        let mut mapper = QueryMapper::new(0, hits, &backend, mapper_config(Some(10.0)));
        mapper.target_parallel = true;
        mapper.init();
        let mut dp_stat = DpStat::default();
        let mut pipeline = Pipeline::new(
            mapper,
            &mut dp_stat,
            PipelineConfig {
                threads: 2,
                ..PipelineConfig::default()
            },
        );
        let mut statistics = Statistics::default();
        pipeline.run(
            &mut statistics,
            RunConfig {
                max_target_seqs: 10,
            },
        );
        assert_eq!(pipeline.mapper.targets.len(), 1);
        assert!(backend.calls.borrow().iter().all(|call| call.2));
        assert_eq!(statistics.get(StatValue::TargetHits2), 1);
    }

    #[test]
    fn empty_pipeline_returns_before_statistics_and_swipe() {
        let backend = Backend::new();
        let mapper = QueryMapper::new(0, vec![], &backend, mapper_config(None));
        let mut dp_stat = DpStat::default();
        let mut pipeline = Pipeline::new(mapper, &mut dp_stat, PipelineConfig::default());
        let mut statistics = Statistics::default();
        pipeline.run(&mut statistics, RunConfig { max_target_seqs: 1 });
        assert_eq!(pipeline.config.padding, Some(32));
        assert_eq!(statistics.get(StatValue::TargetHits0), 0);
        assert!(backend.calls.borrow().is_empty());
    }
}
