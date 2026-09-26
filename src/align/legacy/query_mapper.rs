//! Legacy query-to-target mapping orchestration.
//!
//! This is the Rust counterpart of `align/legacy/query_mapper.{h,cpp}`.  The
//! C++ implementation reaches into process-global query, database, scoring,
//! and output objects.  Those dependencies are explicit traits here, while
//! preserving the ordering, culling, filtering, and output lifecycle.

use std::cmp::Ordering;

use crate::align::hsp::Hsp;
use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::translate::{Frame, Strand, TranslatedPosition};
use crate::basic::value::{BlockId, Letter, OId};
use crate::dp::ungapped::DiagonalSegment;
use crate::output::target_culling::{
    CullingResult, CullingTarget, TargetCulling, TargetCullingConfig,
};
use crate::search::hit::Hit;
use crate::stats::hauser_correction::HauserCorrection;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::geo::DiagonalSegmentT;
use crate::util::hsp::ApproxHsp;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SeedHit {
    pub frame: u32,
    pub subject: u32,
    pub subject_pos: u32,
    pub query_pos: u32,
    pub ungapped: DiagonalSegment,
    pub prefix_score: u32,
}

impl SeedHit {
    pub fn new(
        frame: u32,
        subject: u32,
        subject_pos: u32,
        query_pos: u32,
        ungapped: DiagonalSegment,
    ) -> Self {
        Self {
            frame,
            subject,
            subject_pos,
            query_pos,
            prefix_score: ungapped.score as u32,
            ungapped,
        }
    }

    pub fn diagonal(&self) -> i32 {
        self.query_pos as i32 - self.subject_pos as i32
    }

    /// C++ `SeedHit::operator<` (higher score sorts first).
    pub fn compare_score(&self, rhs: &Self) -> Ordering {
        rhs.ungapped.score.cmp(&self.ungapped.score)
    }

    pub fn is_enveloped(&self, hsps: &[Hsp], dna_len: i32, query_translated: bool) -> bool {
        let d = self.diagonal_segment();
        hsps.iter()
            .any(|hsp| hsp.envelopes(&d, dna_len, query_translated))
    }

    pub fn diagonal_segment(&self) -> DiagonalSegmentT {
        DiagonalSegmentT::from_diagonal_segment(
            &self.ungapped,
            Frame::from_index(self.frame as i32),
        )
    }

    pub fn query_source_range(
        &self,
        dna_len: i32,
        query_translated: bool,
    ) -> crate::util::interval::Interval {
        self.diagonal_segment()
            .query_absolute_range(dna_len, query_translated)
    }

    pub fn strand(&self) -> Strand {
        Frame::from_index(self.frame as i32).strand
    }

    pub fn compare_pos(x: &Self, y: &Self) -> Ordering {
        x.ungapped.subject_end().cmp(&y.ungapped.subject_end())
    }

    pub fn compare_diag(x: &Self, y: &Self) -> Ordering {
        x.frame
            .cmp(&y.frame)
            .then_with(|| x.diagonal().cmp(&y.diagonal()))
            .then_with(|| x.ungapped.j.cmp(&y.ungapped.j))
    }

    pub fn compare_diag_strand(x: &Self, y: &Self) -> Ordering {
        (x.strand() as u8)
            .cmp(&(y.strand() as u8))
            .then_with(|| x.diagonal().cmp(&y.diagonal()))
            .then_with(|| x.ungapped.j.cmp(&y.ungapped.j))
    }

    pub fn compare_diag_strand2(x: &Self, y: &Self) -> Ordering {
        (x.strand() as u8)
            .cmp(&(y.strand() as u8))
            .then_with(|| x.diagonal().cmp(&y.diagonal()))
            .then_with(|| x.subject_pos.cmp(&y.subject_pos))
    }

    pub fn frame_of(x: &Self) -> u32 {
        x.frame
    }
}

#[derive(Debug, Clone)]
pub struct Target {
    pub subject_block_id: BlockId,
    pub subject: Vec<Letter>,
    pub filter_score: i32,
    pub filter_evalue: f64,
    pub filter_time: f32,
    pub outranked: bool,
    pub begin: usize,
    pub end: usize,
    pub hsps: Vec<Hsp>,
    pub ts: Vec<ApproxHsp>,
    pub top_hit: SeedHit,
    pub taxon_rank_ids: Vec<u32>,
}

impl Target {
    pub const INTERVAL: i32 = 64;

    pub fn filter_only(filter_score: i32, filter_evalue: f64) -> Self {
        Self {
            subject_block_id: 0,
            subject: Vec::new(),
            filter_score,
            filter_evalue,
            filter_time: 0.0,
            outranked: false,
            begin: 0,
            end: 0,
            hsps: Vec::new(),
            ts: Vec::new(),
            top_hit: SeedHit::default(),
            taxon_rank_ids: Vec::new(),
        }
    }

    pub fn new(
        begin: usize,
        subject_block_id: BlockId,
        subject: Vec<Letter>,
        mut taxon_rank_ids: Vec<u32>,
    ) -> Self {
        taxon_rank_ids.sort_unstable();
        taxon_rank_ids.dedup();
        Self {
            subject_block_id,
            subject,
            filter_score: 0,
            filter_evalue: f64::MAX,
            filter_time: 0.0,
            outranked: false,
            begin,
            end: 0,
            hsps: Vec::new(),
            ts: Vec::new(),
            top_hit: SeedHit::default(),
            taxon_rank_ids,
        }
    }

    pub fn compare_evalue(lhs: &Self, rhs: &Self) -> Ordering {
        if lhs.filter_evalue < rhs.filter_evalue {
            Ordering::Less
        } else if rhs.filter_evalue < lhs.filter_evalue {
            Ordering::Greater
        } else if lhs.filter_evalue == rhs.filter_evalue {
            Self::compare_score(lhs, rhs)
        } else {
            // C++ `<` returns false in both directions for NaN, so stable_sort
            // retains the input order rather than applying the score tie-break.
            Ordering::Equal
        }
    }

    pub fn compare_score(lhs: &Self, rhs: &Self) -> Ordering {
        rhs.filter_score
            .cmp(&lhs.filter_score)
            .then_with(|| lhs.subject_block_id.cmp(&rhs.subject_block_id))
    }

    pub fn fill_source_ranges(&mut self, query_len: usize, query_translated: bool) {
        for hsp in &mut self.ts {
            hsp.query_source_range = TranslatedPosition::absolute_interval(
                TranslatedPosition::new(hsp.query_range.begin, Frame::from_index(hsp.frame)),
                TranslatedPosition::new(hsp.query_range.end, Frame::from_index(hsp.frame)),
                query_len as i32,
                query_translated,
            );
        }
    }

    pub fn envelopes(&self, t: &ApproxHsp, p: f64) -> bool {
        self.ts
            .iter()
            .any(|hsp| t.query_source_range.overlap_factor(&hsp.query_source_range) >= p)
    }

    pub fn is_enveloped(&self, target: &Self, p: f64) -> bool {
        self.ts.iter().all(|hsp| target.envelopes(hsp, p))
    }

    pub fn is_enveloped_by_any(&self, targets: &[Target], p: f64, min_score: i32) -> bool {
        targets
            .iter()
            .any(|target| self.is_enveloped(target, p) && target.filter_score >= min_score)
    }

    pub fn add_ranges(&self, bins: &mut [i32]) {
        for hsp in &self.hsps {
            let i0 = hsp.query_source_range.begin / Self::INTERVAL;
            let i1 = (hsp.query_source_range.end / Self::INTERVAL).min(bins.len() as i32 - 1);
            for i in i0..=i1 {
                bins[i as usize] = bins[i as usize].max(hsp.score);
            }
        }
    }

    pub fn is_outranked(&self, bins: &[i32], threshold: f64) -> bool {
        for hsp in &self.hsps {
            let i0 = hsp.query_source_range.begin / Self::INTERVAL;
            let i1 = (hsp.query_source_range.end / Self::INTERVAL).min(bins.len() as i32 - 1);
            for i in i0..=i1 {
                if hsp.score as f64 >= bins[i as usize] as f64 * threshold {
                    return false;
                }
            }
        }
        true
    }

    pub fn inner_culling(&mut self) {
        self.hsps.sort_by(|a, b| {
            if a.less_than_score_position(b) {
                Ordering::Less
            } else if b.less_than_score_position(a) {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        });
        if let Some(first) = self.hsps.first() {
            self.filter_score = first.score;
            self.filter_evalue = first.evalue;
        } else {
            self.filter_score = 0;
            self.filter_evalue = f64::MAX;
        }
        let mut i = 0;
        while i < self.hsps.len() {
            let enveloped = self.hsps[..i]
                .iter()
                .any(|prior| self.hsps[i].query_range_enveloped_by(prior, 0.5));
            if enveloped {
                self.hsps.remove(i);
            } else {
                i += 1;
            }
        }
    }

    pub fn apply_filters(&mut self, dna_len: i32, subject_len: i32, filters: HspFilters<'_>) {
        let _query_title = filters.query_title; // C++ accepts but does not use this parameter.
        self.hsps.retain(|hsp| {
            !(hsp.id_percent() < filters.min_id
                || hsp.query_cover_percent(dna_len as u32) < filters.query_cover
                || hsp.subject_cover_percent(subject_len as u32) < filters.subject_cover)
        });
    }
}

#[derive(Debug, Clone, Copy)]
pub struct HspFilters<'a> {
    pub min_id: f64,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_title: &'a str,
}

#[derive(Debug, Clone)]
pub struct QueryMapperConfig {
    pub query_contexts: u32,
    pub query_translated: bool,
    pub log_query: bool,
    pub hauser: bool,
    pub cbs_window: usize,
    pub taxon_k: u32,
    pub toppercent: Option<f64>,
    pub min_id: f64,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub max_hsps: u32,
    pub culling: TargetCullingConfig,
}

impl Default for QueryMapperConfig {
    fn default() -> Self {
        Self {
            query_contexts: 1,
            query_translated: false,
            log_query: false,
            hauser: false,
            cbs_window: 40,
            taxon_k: 0,
            toppercent: None,
            min_id: 0.0,
            query_cover: 0.0,
            subject_cover: 0.0,
            max_hsps: 0,
            culling: TargetCullingConfig::default(),
        }
    }
}

/// Explicit adapter for the C++ query/target/database and score globals.
pub trait QueryMapperBackend {
    fn source_query_len(&self, query_id: u32) -> u32;
    fn query_seq(&self, query_id: u32, frame: u32) -> &[Letter];
    fn query_source_seq(&self, query_id: u32) -> &[Letter];
    fn query_title(&self, query_id: u32) -> &str;
    fn query_oid(&self, query_id: u32) -> OId;
    fn local_target_position(&self, packed_position: u64) -> (usize, usize);
    fn target_seq(&self, subject_id: BlockId) -> &[Letter];
    fn target_oid(&self, subject_id: BlockId) -> OId;
    fn target_taxon_rank_ids(&self, oid: OId) -> Vec<u32>;
    fn target_title(&self, subject_id: BlockId, oid: OId) -> String;
    fn target_dict_id(&self, current_ref_block: usize, subject_id: BlockId) -> usize;
    fn score_matrix(&self) -> &ScoreMatrix;
    fn xdrop_ungapped(
        &self,
        query: &[Letter],
        subject: &[Letter],
        query_pos: usize,
        subject_pos: usize,
    ) -> DiagonalSegment;
    fn report_cutoff(&self, score: i32, evalue: f64) -> bool;
    fn bitscore(&self, score: i32) -> f64;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputMode {
    Regular,
    Daa,
    Intermediate,
}

pub struct QueryIntro<'a> {
    pub mode: OutputMode,
    pub unaligned: bool,
    pub query_id: u32,
    pub query_oid: OId,
    pub title: &'a str,
    pub source_seq: &'a [Letter],
}

pub struct MatchContext<'a> {
    pub mode: OutputMode,
    pub hsp: &'a Hsp,
    pub query_id: u32,
    pub query_oid: OId,
    pub query_title: &'a str,
    /// The context selected from C++ `TranslatedSequence` by `hsp.frame`.
    pub query: &'a [Letter],
    pub database_id: OId,
    pub subject_len: u32,
    pub target_title: &'a str,
    pub target_index: u32,
    pub hsp_index: u32,
    pub subject: &'a [Letter],
    pub dict_id: usize,
}

/// Explicit adapter for regular, DAA, and intermediate output writers.
pub trait QueryMapperOutput {
    fn mode(&self) -> OutputMode;
    fn report_unaligned(&self) -> bool {
        false
    }
    fn begin_query(&mut self, intro: QueryIntro<'_>) -> usize;
    fn write_match(&mut self, context: MatchContext<'_>);
    fn finish_query(&mut self, seek_pos: usize);
}

#[derive(Debug, Clone, Copy)]
pub struct GenerateOutputConfig {
    pub max_target_seqs: i64,
    pub blocked_processing: bool,
    pub current_ref_block: usize,
}

pub struct QueryMapper<'a, B: QueryMapperBackend> {
    pub source_hits: Vec<Hit>,
    pub query_id: u32,
    pub targets_finished: u32,
    pub next_target: u32,
    pub source_query_len: u32,
    pub unaligned_from: u32,
    pub targets: Vec<Target>,
    pub seed_hits: Vec<SeedHit>,
    pub query_cb: Vec<HauserCorrection>,
    pub target_parallel: bool,
    pub config: QueryMapperConfig,
    backend: &'a B,
}

impl<'a, B: QueryMapperBackend> QueryMapper<'a, B> {
    pub fn new(
        query_id: usize,
        source_hits: Vec<Hit>,
        backend: &'a B,
        config: QueryMapperConfig,
    ) -> Self {
        Self {
            source_hits,
            query_id: query_id as u32,
            targets_finished: 0,
            next_target: 0,
            source_query_len: backend.source_query_len(query_id as u32),
            unaligned_from: 0,
            targets: Vec::new(),
            seed_hits: Vec::new(),
            query_cb: Vec::new(),
            target_parallel: false,
            config,
            backend,
        }
    }

    pub fn init(&mut self) {
        if self.config.log_query {
            eprintln!(
                "Query = {}\t{}",
                self.backend.query_title(self.query_id),
                self.query_id
            );
        }
        if self.config.hauser {
            for frame in 0..self.config.query_contexts {
                self.query_cb.push(HauserCorrection::new(
                    self.query_seq(frame),
                    self.backend.score_matrix(),
                    self.config.cbs_window,
                ));
            }
        }
        let count = self.count_targets();
        self.targets.reserve(count);
        if count != 0 {
            self.load_targets();
        }
        debug_assert_eq!(self.targets.len(), count);
    }

    fn count_targets(&mut self) -> usize {
        self.source_hits.sort_by(Hit::cmp_subject);
        let mut subject_id = usize::MAX;
        let mut n_subject = 0usize;
        for hit in &self.source_hits {
            let (local_subject, subject_pos) = self.backend.local_target_position(hit.subject);
            let frame = hit.query % self.config.query_contexts;
            let diagonal = if self.target_parallel {
                DiagonalSegment::default()
            } else {
                self.backend.xdrop_ungapped(
                    self.backend.query_seq(self.query_id, frame),
                    self.backend.target_seq(local_subject as BlockId),
                    hit.seed_offset as usize,
                    subject_pos,
                )
            };
            if self.target_parallel || diagonal.score > 0 {
                if local_subject != subject_id {
                    subject_id = local_subject;
                    n_subject += 1;
                }
                self.seed_hits.push(SeedHit::new(
                    frame,
                    local_subject as u32,
                    subject_pos as u32,
                    hit.seed_offset,
                    diagonal,
                ));
            }
        }
        n_subject
    }

    fn load_targets(&mut self) {
        let mut subject_id = u32::MAX;
        for (i, seed_hit) in self.seed_hits.iter().enumerate() {
            if seed_hit.subject != subject_id {
                if let Some(previous) = self.targets.last_mut() {
                    previous.end = i;
                }
                let oid = self.backend.target_oid(seed_hit.subject);
                let taxons = if self.config.taxon_k != 0 {
                    self.backend.target_taxon_rank_ids(oid)
                } else {
                    Vec::new()
                };
                self.targets.push(Target::new(
                    i,
                    seed_hit.subject,
                    self.backend.target_seq(seed_hit.subject).to_vec(),
                    taxons,
                ));
                subject_id = seed_hit.subject;
            }
        }
        self.targets.last_mut().expect("nonempty seed hits").end = self.seed_hits.len();
    }

    pub fn rank_targets(&mut self, ratio: f64, factor: f64, max_target_seqs: i64) {
        if self.config.taxon_k != 0 && self.config.toppercent.is_none() {
            return;
        }
        self.targets.sort_by(Target::compare_score);
        assert!(!self.targets.is_empty() && max_target_seqs > 0);
        let score = if let Some(toppercent) = self.config.toppercent {
            (self.targets[0].filter_score as f64 * (1.0 - toppercent / 100.0) * ratio) as i32
        } else {
            let min_idx = (self.targets.len() as i64).min(max_target_seqs) as usize;
            (self.targets[min_idx - 1].filter_score as f64 * ratio) as i32
        };
        let cap = if self.config.toppercent.is_some() || max_target_seqs == i64::MAX {
            i64::MAX
        } else {
            (max_target_seqs as f64 * factor) as i64
        };
        let keep = self
            .targets
            .iter()
            .enumerate()
            .position(|(i, target)| target.filter_score < score || i as i64 >= cap)
            .unwrap_or(self.targets.len());
        self.targets.truncate(keep);
    }

    fn culling_target(target: &Target) -> CullingTarget<'_> {
        CullingTarget {
            filter_score: target.filter_score,
            taxon_rank_ids: &target.taxon_rank_ids,
            hsps: &target.hsps,
        }
    }

    fn culling_config(&self) -> TargetCullingConfig {
        TargetCullingConfig {
            taxon_k: self.config.taxon_k,
            toppercent: self.config.toppercent,
            ..self.config.culling
        }
    }

    pub fn score_only_culling(&mut self, max_target_seqs: i64) {
        if self.config.toppercent.is_none() {
            self.targets.sort_by(Target::compare_evalue);
        } else {
            self.targets.sort_by(Target::compare_score);
        }
        let culling_config = self.culling_config();
        let mut culling = TargetCulling::from_config(max_target_seqs, &culling_config);
        let mut i = 0usize;
        while i < self.targets.len() {
            if !self
                .backend
                .report_cutoff(self.targets[i].filter_score, self.targets[i].filter_evalue)
            {
                break;
            }
            let (code, coverage) = culling.cull_target(
                Self::culling_target(&self.targets[i]),
                &culling_config,
                |score| self.backend.bitscore(score),
            );
            match code {
                CullingResult::Finished => break,
                CullingResult::Next => {
                    self.targets.remove(i);
                }
                CullingResult::Include => {
                    if coverage < 0.1 {
                        culling.add_target(
                            Self::culling_target(&self.targets[i]),
                            &culling_config,
                            |score| self.backend.bitscore(score),
                        );
                    }
                    i += 1;
                }
            }
        }
        self.targets.truncate(i);
    }

    pub fn generate_output<O: QueryMapperOutput>(
        &mut self,
        output: &mut O,
        stat: &mut Statistics,
        cfg: GenerateOutputConfig,
    ) -> bool {
        if self.config.toppercent.is_none() {
            self.targets.sort_by(Target::compare_evalue);
        } else {
            self.targets.sort_by(Target::compare_score);
        }
        let culling_config = self.culling_config();
        let mut culling = TargetCulling::from_config(cfg.max_target_seqs, &culling_config);
        let backend = self.backend;
        let query_title = backend.query_title(self.query_id);
        let query_source_seq = backend.query_source_seq(self.query_id);
        let query_oid = backend.query_oid(self.query_id);
        let mode = if cfg.blocked_processing {
            OutputMode::Intermediate
        } else {
            output.mode()
        };
        let mut n_hsp = 0u32;
        let mut n_target_seq = 0u32;
        let mut seek_pos = 0usize;

        for target in &mut self.targets {
            let subject_id = target.subject_block_id;
            let database_id = backend.target_oid(subject_id);
            let target_title = if cfg.blocked_processing {
                String::new()
            } else {
                backend.target_title(subject_id, database_id)
            };
            // Intermediate output resolves the dictionary entry before
            // filtering, exactly where the C++ branch does. DAA resolves it
            // separately for every emitted HSP below.
            let intermediate_dict_id = if cfg.blocked_processing {
                Some(backend.target_dict_id(cfg.current_ref_block, subject_id))
            } else {
                None
            };
            let subject_len = backend.target_seq(subject_id).len() as u32;
            target.apply_filters(
                self.source_query_len as i32,
                subject_len as i32,
                HspFilters {
                    min_id: self.config.min_id,
                    query_cover: self.config.query_cover,
                    subject_cover: self.config.subject_cover,
                    query_title,
                },
            );
            if target.hsps.is_empty() {
                continue;
            }
            match culling
                .cull_target(Self::culling_target(target), &culling_config, |score| {
                    backend.bitscore(score)
                })
                .0
            {
                CullingResult::Next => continue,
                CullingResult::Finished => break,
                CullingResult::Include => {}
            }
            culling.add_target(Self::culling_target(target), &culling_config, |score| {
                backend.bitscore(score)
            });
            let mut hit_hsps = 0u32;
            for hsp in &target.hsps {
                if self.config.max_hsps > 0 && hit_hsps >= self.config.max_hsps {
                    break;
                }
                if n_hsp == 0 {
                    seek_pos = output.begin_query(QueryIntro {
                        mode,
                        unaligned: false,
                        query_id: self.query_id,
                        query_oid,
                        title: query_title,
                        source_seq: query_source_seq,
                    });
                }
                output.write_match(MatchContext {
                    mode,
                    hsp,
                    query_id: self.query_id,
                    query_oid,
                    query_title,
                    query: backend.query_seq(self.query_id, hsp.frame as u32),
                    database_id,
                    subject_len,
                    target_title: &target_title,
                    target_index: n_target_seq,
                    hsp_index: hit_hsps,
                    subject: &target.subject,
                    dict_id: match mode {
                        OutputMode::Intermediate => intermediate_dict_id.unwrap(),
                        OutputMode::Daa => {
                            backend.target_dict_id(cfg.current_ref_block, subject_id)
                        }
                        OutputMode::Regular => 0,
                    },
                });
                n_hsp += 1;
                hit_hsps += 1;
            }
            n_target_seq += 1;
        }

        if n_hsp > 0 {
            output.finish_query(seek_pos);
        } else if !cfg.blocked_processing && mode != OutputMode::Daa && output.report_unaligned() {
            let seek = output.begin_query(QueryIntro {
                mode,
                unaligned: true,
                query_id: self.query_id,
                query_oid,
                title: query_title,
                source_seq: query_source_seq,
            });
            output.finish_query(seek);
        }
        if !cfg.blocked_processing {
            stat.inc(StatValue::Matches, n_hsp as i64);
            stat.inc(StatValue::Pairwise, n_target_seq as i64);
            if n_hsp > 0 {
                stat.inc(StatValue::Aligned, 1);
            }
        }
        n_hsp > 0
    }

    pub fn n_targets(&self) -> i64 {
        self.targets.len() as i64
    }
    pub fn finished(&self) -> bool {
        self.targets_finished as usize == self.targets.len()
    }
    pub fn query_seq(&self, frame: u32) -> &[Letter] {
        self.backend.query_seq(self.query_id, frame)
    }
    pub fn query_source_seq(&self) -> &[Letter] {
        self.backend.query_source_seq(self.query_id)
    }
    pub fn backend(&self) -> &B {
        self.backend
    }
    pub fn fill_source_ranges(&mut self) {
        for target in &mut self.targets {
            target.fill_source_ranges(self.source_query_len as usize, self.config.query_translated);
        }
    }
}

/// Rust counterpart for the abstract C++ `QueryMapper::run` slot.
pub trait QueryMapperRun {
    fn run(&mut self, stat: &mut Statistics, cfg: GenerateOutputConfig);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::interval::Interval;

    struct Backend {
        matrix: ScoreMatrix,
        queries: Vec<Vec<Letter>>,
        targets: Vec<Vec<Letter>>,
    }

    impl Backend {
        fn new() -> Self {
            Self {
                matrix: ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap(),
                queries: vec![vec![0; 30], vec![1; 30]],
                targets: vec![vec![0; 40], vec![1; 50], vec![2; 60]],
            }
        }
    }

    impl QueryMapperBackend for Backend {
        fn source_query_len(&self, _query_id: u32) -> u32 {
            30
        }
        fn query_seq(&self, _query_id: u32, frame: u32) -> &[Letter] {
            &self.queries[frame as usize]
        }
        fn query_source_seq(&self, _query_id: u32) -> &[Letter] {
            &self.queries[0]
        }
        fn query_title(&self, _query_id: u32) -> &str {
            "query"
        }
        fn query_oid(&self, query_id: u32) -> OId {
            query_id as OId + 100
        }
        fn local_target_position(&self, packed: u64) -> (usize, usize) {
            ((packed / 100) as usize, (packed % 100) as usize)
        }
        fn target_seq(&self, subject_id: BlockId) -> &[Letter] {
            &self.targets[subject_id as usize]
        }
        fn target_oid(&self, subject_id: BlockId) -> OId {
            subject_id as OId + 1000
        }
        fn target_taxon_rank_ids(&self, oid: OId) -> Vec<u32> {
            vec![(oid - 900) as u32]
        }
        fn target_title(&self, subject_id: BlockId, _oid: OId) -> String {
            format!("target-{subject_id}")
        }
        fn target_dict_id(&self, block: usize, subject_id: BlockId) -> usize {
            block * 10 + subject_id as usize
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
            DiagonalSegment::new(
                query_pos as i32,
                subject_pos as i32,
                5,
                if subject_pos == 9 { 0 } else { 20 },
            )
        }
        fn report_cutoff(&self, score: i32, evalue: f64) -> bool {
            score >= 20 && evalue <= 10.0
        }
        fn bitscore(&self, score: i32) -> f64 {
            score as f64
        }
    }

    fn hsp(score: i32, begin: i32, end: i32) -> Hsp {
        let mut h = Hsp::new();
        h.score = score;
        h.evalue = 1.0 / score as f64;
        h.query_source_range = Interval::new(begin, end);
        h.query_range = h.query_source_range;
        h.subject_range = Interval::new(0, end - begin);
        h.length = end - begin;
        h.identities = h.length;
        h
    }

    #[test]
    fn seed_hit_geometry_and_comparators_match_header_inline_methods() {
        let forward = SeedHit::new(0, 2, 3, 8, DiagonalSegment::new(8, 3, 5, 40));
        let reverse = SeedHit::new(4, 2, 1, 9, DiagonalSegment::new(9, 1, 2, 30));
        assert_eq!(forward.diagonal(), 5);
        assert_eq!(forward.strand(), Strand::Forward);
        assert_eq!(reverse.strand(), Strand::Reverse);
        assert_eq!(SeedHit::frame_of(&reverse), 4);
        assert_eq!(forward.compare_score(&reverse), Ordering::Less);
        assert_eq!(
            SeedHit::compare_diag_strand(&forward, &reverse),
            Ordering::Less
        );
        assert_eq!(forward.query_source_range(30, false), Interval::new(8, 13));
        assert!(forward.is_enveloped(&[hsp(50, 8, 13)], 30, false));
    }

    #[test]
    fn target_enveloping_bins_inner_culling_and_filters_are_faithful() {
        let mut target = Target::new(0, 1, vec![], vec![]);
        target.ts = vec![ApproxHsp::from_query_source_range(Interval::new(0, 10))];
        let mut covering = Target::new(0, 2, vec![], vec![]);
        covering.ts = vec![ApproxHsp::from_query_source_range(Interval::new(0, 20))];
        assert!(target.is_enveloped(&covering, 1.0));
        assert!(target.is_enveloped_by_any(&[covering], 1.0, 0));

        target.hsps = vec![hsp(100, 0, 100), hsp(80, 25, 75), hsp(70, 128, 160)];
        target.inner_culling();
        assert_eq!(
            target.hsps.iter().map(|h| h.score).collect::<Vec<_>>(),
            vec![100, 70]
        );
        assert_eq!(target.filter_score, 100);
        let mut bins = vec![0; 4];
        target.add_ranges(&mut bins);
        assert_eq!(bins, vec![100, 100, 70, 0]);
        assert!(!target.is_outranked(&bins, 1.0));
        target.apply_filters(
            200,
            200,
            HspFilters {
                min_id: 90.0,
                query_cover: 25.0,
                subject_cover: 10.0,
                query_title: "ignored",
            },
        );
        assert_eq!(target.hsps.len(), 1);
    }

    #[test]
    fn init_sorts_hits_drops_failed_extensions_and_loads_group_bounds() {
        let backend = Backend::new();
        let hits = vec![Hit::new(0, 109, 4), Hit::new(1, 5, 3), Hit::new(0, 102, 2)];
        let mut config = QueryMapperConfig::default();
        config.query_contexts = 2;
        config.taxon_k = 1;
        let mut mapper = QueryMapper::new(0, hits, &backend, config);
        mapper.init();
        assert_eq!(mapper.seed_hits.len(), 2);
        assert_eq!(mapper.targets.len(), 2);
        assert_eq!((mapper.targets[0].begin, mapper.targets[0].end), (0, 1));
        assert_eq!((mapper.targets[1].begin, mapper.targets[1].end), (1, 2));
        assert_eq!(mapper.targets[1].taxon_rank_ids, vec![101]);

        let hits = vec![Hit::new(0, 109, 4)];
        let mut mapper = QueryMapper::new(0, hits, &backend, QueryMapperConfig::default());
        mapper.target_parallel = true;
        mapper.init();
        assert_eq!(mapper.seed_hits[0].ungapped, DiagonalSegment::default());
    }

    #[test]
    fn rank_and_score_only_culling_follow_score_evalue_and_cutoff_rules() {
        let backend = Backend::new();
        let mut mapper = QueryMapper::new(0, vec![], &backend, QueryMapperConfig::default());
        mapper.targets = vec![
            Target::filter_only(80, 0.3),
            Target::filter_only(100, 0.2),
            Target::filter_only(90, 0.1),
        ];
        mapper.rank_targets(1.0, 2.0, 2);
        assert_eq!(
            mapper
                .targets
                .iter()
                .map(|t| t.filter_score)
                .collect::<Vec<_>>(),
            vec![100, 90]
        );
        mapper.score_only_culling(1);
        assert_eq!(mapper.targets.len(), 1);
        assert_eq!(mapper.targets[0].filter_score, 90); // e-value order precedes score without top-percent.

        let mut mapper = QueryMapper::new(0, vec![], &backend, QueryMapperConfig::default());
        mapper.targets.push(Target::filter_only(10, 0.01));
        mapper.score_only_culling(1);
        assert!(mapper.targets.is_empty()); // A failed report cutoff terminates and erases the tail.
    }

    #[derive(Default)]
    struct Output {
        mode: Option<OutputMode>,
        intros: Vec<OutputMode>,
        matches: Vec<(OutputMode, u32, u32, usize)>,
        finishes: usize,
        unaligned: bool,
    }

    impl QueryMapperOutput for Output {
        fn mode(&self) -> OutputMode {
            self.mode.unwrap_or(OutputMode::Regular)
        }
        fn report_unaligned(&self) -> bool {
            self.unaligned
        }
        fn begin_query(&mut self, intro: QueryIntro<'_>) -> usize {
            self.intros.push(intro.mode);
            17
        }
        fn write_match(&mut self, context: MatchContext<'_>) {
            self.matches.push((
                context.mode,
                context.target_index,
                context.hsp_index,
                context.dict_id,
            ));
        }
        fn finish_query(&mut self, seek_pos: usize) {
            assert_eq!(seek_pos, 17);
            self.finishes += 1;
        }
    }

    #[test]
    fn generate_output_preserves_lazy_intro_limits_modes_and_statistics() {
        let backend = Backend::new();
        let mut config = QueryMapperConfig::default();
        config.max_hsps = 1;
        let mut mapper = QueryMapper::new(0, vec![], &backend, config);
        let mut a = Target::new(0, 0, backend.targets[0].clone(), vec![]);
        a.filter_score = 100;
        a.filter_evalue = 0.01;
        a.hsps = vec![hsp(100, 0, 20), hsp(90, 20, 30)];
        let mut b = Target::new(0, 1, backend.targets[1].clone(), vec![]);
        b.filter_score = 90;
        b.filter_evalue = 0.02;
        b.hsps = vec![hsp(90, 0, 20)];
        mapper.targets = vec![a, b];
        let mut output = Output::default();
        let mut stat = Statistics::default();
        assert!(mapper.generate_output(
            &mut output,
            &mut stat,
            GenerateOutputConfig {
                max_target_seqs: 2,
                blocked_processing: false,
                current_ref_block: 3
            }
        ));
        assert_eq!(output.intros, vec![OutputMode::Regular]);
        assert_eq!(output.matches.len(), 2);
        assert_eq!(output.finishes, 1);
        assert_eq!(stat.get(StatValue::Matches), 2);
        assert_eq!(stat.get(StatValue::Pairwise), 2);
        assert_eq!(stat.get(StatValue::Aligned), 1);

        let mut output = Output::default();
        let mut stat = Statistics::default();
        assert!(mapper.generate_output(
            &mut output,
            &mut stat,
            GenerateOutputConfig {
                max_target_seqs: 2,
                blocked_processing: true,
                current_ref_block: 3
            }
        ));
        assert_eq!(output.intros, vec![OutputMode::Intermediate]);
        assert_eq!(output.matches[0].3, 30);
        assert_eq!(stat.get(StatValue::Matches), 0);
    }

    #[test]
    fn regular_unaligned_output_emits_empty_query_but_daa_does_not() {
        let backend = Backend::new();
        let mut mapper = QueryMapper::new(0, vec![], &backend, QueryMapperConfig::default());
        let mut regular = Output {
            unaligned: true,
            ..Output::default()
        };
        assert!(!mapper.generate_output(
            &mut regular,
            &mut Statistics::default(),
            GenerateOutputConfig {
                max_target_seqs: 10,
                blocked_processing: false,
                current_ref_block: 0
            }
        ));
        assert_eq!((regular.intros.len(), regular.finishes), (1, 1));
        let mut daa = Output {
            mode: Some(OutputMode::Daa),
            unaligned: true,
            ..Output::default()
        };
        assert!(!mapper.generate_output(
            &mut daa,
            &mut Statistics::default(),
            GenerateOutputConfig {
                max_target_seqs: 10,
                blocked_processing: false,
                current_ref_block: 0
            }
        ));
        assert!(daa.intros.is_empty());
    }
}
