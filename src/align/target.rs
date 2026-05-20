use super::hsp::Match;
use crate::align::gapped_filter::{
    gapped_filter_seed_hits, load_hits, seed_only_matches, SeedHit, SeedHitList,
};
use crate::align::global_ranking;
use crate::align::hsp::Hsp;
use crate::align::ungapped::{ungapped_stage_seed_hits, UngappedStageConfig};
use crate::basic::consts::MAX_CONTEXT;
use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::value::{BlockId, Letter, SUPER_HARD_MASK, TRUE_AA};
use crate::config::Sensitivity;
use crate::data::block::Block;
use crate::dp::anchored::{anchored_swipe, AnchoredSwipeConfig};
use crate::dp::swipe::{
    bin as swipe_bin, have_coords, swipe, swipe_set, targets as make_dp_targets,
    Anchor as DpAnchor, CarryOver, DpTarget, Flags, HspValues, Params, Targets,
};
use crate::masking::{mask_sequence, MaskingAlgo};
use crate::search::hit::Hit;
use crate::search::sensitivity::ExtensionMode;
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::data_structures::FlatArray;
use crate::util::hsp::{Anchor, ApproxHsp};
use crate::util::sequence::is_fully_masked;
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct Target {
    pub block_id: BlockId,
    pub seq: Vec<Letter>,
    pub filter_score: i32,
    pub filter_evalue: f64,
    pub best_context: i32,
    pub ungapped_score: i32,
    pub hsp: [Vec<Hsp>; MAX_CONTEXT as usize],
    pub matrix: Option<Arc<TargetMatrix>>,
    pub done: bool,
}

impl Target {
    /// Matches C++ `Extension::Target::Target`.
    pub fn new(
        block_id: BlockId,
        seq: &[Letter],
        ungapped_score: i32,
        matrix: Option<Arc<TargetMatrix>>,
    ) -> Self {
        Self {
            block_id,
            seq: seq.to_vec(),
            filter_score: 0,
            filter_evalue: f64::MAX,
            best_context: 0,
            ungapped_score,
            hsp: std::array::from_fn(|_| Vec::new()),
            matrix,
            done: false,
        }
    }

    /// Matches C++ `Target::add_hit(Hsp&&)`.
    pub fn add_hit(&mut self, hsp: Hsp) {
        if hsp.evalue < self.filter_evalue {
            self.filter_evalue = hsp.evalue;
            self.filter_score = hsp.score;
            self.best_context = hsp.frame;
        }
        self.hsp[hsp.frame as usize].push(hsp);
    }

    /// Matches C++ `Target::add_hit(list<Hsp>&, iterator)`.
    pub fn add_hit_from_vec(&mut self, hsps: &mut Vec<Hsp>, index: usize) {
        let hsp = hsps.remove(index);
        let frame = hsp.frame as usize;
        self.hsp[frame].push(hsp);
        let last = self.hsp[frame].last().unwrap();
        if last.score > self.filter_score {
            self.filter_evalue = last.evalue;
            self.filter_score = last.score;
            self.best_context = last.frame;
        }
    }

    /// Matches C++ `Target::add_hit(const ApproxHsp&, Loc)`.
    pub fn add_approx_hit(&mut self, h: &ApproxHsp, qlen: i32, score_matrix: &ScoreMatrix) {
        let hsp = Hsp::from_approx_hsp(h, qlen, self.seq.len() as i32, score_matrix);
        let frame = h.frame as usize;
        if h.evalue < self.filter_evalue {
            self.filter_evalue = h.evalue;
            self.filter_score = h.score;
        }
        self.hsp[frame].push(hsp);
        self.done = true;
    }

    pub fn comp_evalue(t: &Target, u: &Target) -> bool {
        t.filter_evalue < u.filter_evalue
            || (t.filter_evalue == u.filter_evalue && Self::comp_score(t, u))
    }

    pub fn comp_score(t: &Target, u: &Target) -> bool {
        t.filter_score > u.filter_score
            || (t.filter_score == u.filter_score && t.block_id < u.block_id)
    }

    /// Matches C++ `Target::inner_culling`.
    pub fn inner_culling(
        &mut self,
        max_hsps: u32,
        inner_culling_overlap: f64,
        query_contexts: usize,
    ) {
        if max_hsps == 1 {
            for i in 0..MAX_CONTEXT as usize {
                if i as i32 == self.best_context {
                    self.hsp[i].sort_by(|a, b| {
                        if a.less_than_score_position(b) {
                            std::cmp::Ordering::Less
                        } else if b.less_than_score_position(a) {
                            std::cmp::Ordering::Greater
                        } else {
                            std::cmp::Ordering::Equal
                        }
                    });
                    self.hsp[i].truncate(1);
                } else {
                    self.hsp[i].clear();
                }
            }
            return;
        }
        let mut hsps = Vec::new();
        for frame in 0..query_contexts {
            hsps.append(&mut self.hsp[frame]);
        }
        inner_hsp_culling(&mut hsps, max_hsps, inner_culling_overlap);
        for hsp in hsps {
            self.hsp[hsp.frame as usize].push(hsp);
        }
    }

    /// Matches C++ `Target::max_hsp_culling`.
    pub fn max_hsp_culling(&mut self, max_hsps: u32, query_contexts: usize) {
        for frame in 0..query_contexts {
            max_hsp_culling(&mut self.hsp[frame], max_hsps);
        }
    }
}

/// Matches C++ `output_range()` for `Target` ranges.
pub fn output_range_targets<F>(
    targets: &[Target],
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> usize
where
    F: FnMut(i32) -> f64,
{
    if targets.is_empty() || targets[0].filter_evalue == f64::MAX {
        return 0;
    }
    if let Some(toppercent) = toppercent {
        let cutoff = ((1.0 - toppercent / 100.0) * bitscore(targets[0].filter_score)).max(1.0);
        let mut i = 0;
        while i < targets.len() && bitscore(targets[i].filter_score) >= cutoff {
            i += 1;
        }
        i
    } else {
        let mut i = (max_target_seqs as usize).min(targets.len());
        while i > 1 && targets[i - 1].filter_evalue == f64::MAX {
            i -= 1;
        }
        i
    }
}

/// Matches C++ `culling(std::vector<Target>&, bool, const Search::Config&)`.
pub fn culling_targets<F>(
    targets: &mut Vec<Target>,
    sort_only: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) where
    F: FnMut(i32) -> f64,
{
    if toppercent.is_some() {
        targets.sort_by(|a, b| {
            b.filter_score
                .cmp(&a.filter_score)
                .then_with(|| a.block_id.cmp(&b.block_id))
        });
    } else {
        targets.sort_by(|a, b| {
            a.filter_evalue
                .total_cmp(&b.filter_evalue)
                .then_with(|| b.filter_score.cmp(&a.filter_score))
                .then_with(|| a.block_id.cmp(&b.block_id))
        });
    }
    if !sort_only {
        let end = output_range_targets(targets, max_target_seqs, toppercent, &mut bitscore);
        targets.truncate(end);
    }
}

/// Matches C++ `append_hits(vector<Target>&, begin, end, bool, const Search::Config&)`.
pub fn append_hits_targets<F>(
    targets: &mut Vec<Target>,
    mut hits: Vec<Target>,
    with_culling: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> bool
where
    F: FnMut(i32) -> f64,
{
    if hits.is_empty() {
        return false;
    }
    let mut new_hits = toppercent.is_none() && targets.len() < max_target_seqs as usize;
    let mut append = !with_culling || new_hits;

    culling_targets(targets, append, max_target_seqs, toppercent, &mut bitscore);

    let mut max_score = 0;
    let mut min_evalue = f64::MAX;
    for hit in &hits {
        max_score = max_score.max(hit.filter_score);
        min_evalue = min_evalue.min(hit.filter_evalue);
    }

    let range_end = output_range_targets(targets, max_target_seqs, toppercent, &mut bitscore);
    if targets.is_empty()
        || range_end == 0
        || (toppercent.is_none() && min_evalue <= targets[range_end - 1].filter_evalue)
        || (toppercent.is_some()
            && max_score
                >= ((1.0 - toppercent.unwrap() / 100.0)
                    * targets[range_end - 1].filter_score as f64) as i32)
    {
        append = true;
        new_hits = true;
    }

    if append {
        targets.append(&mut hits);
    }
    new_hits
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GappedScoreConfig {
    pub padding: i32,
    pub narrow_band_cov: f64,
    pub narrow_band_factor: f64,
    pub min_band_overlap: f64,
    pub anchored_swipe: bool,
    pub query_contexts: usize,
    pub cutoff_score_8bit: i32,
    pub max_swipe_dp: i64,
    pub approx_backtrace: bool,
    pub max_hsps: u32,
    pub inner_culling_overlap: f64,
    pub comp_based_stats_hauser: bool,
    pub comp_based_stats_matrix_adjust: bool,
    pub sensitivity: Sensitivity,
    pub max_evalue: f64,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub approx_min_id: f64,
    pub cbs_matrix_scale: i32,
    pub min_id: f64,
    pub no_self_hits: bool,
    pub max_target_seqs: i64,
    pub toppercent: Option<f64>,
    pub no_ranking: bool,
    pub global_ranking_targets: i64,
    pub ext_chunk_size: i64,
    pub mapany: bool,
    pub target_hard_cap: Option<i64>,
    pub ranking_score_drop_factor: f64,
    pub ranking_cutoff_bitscore: f64,
    pub cluster_threshold_present: bool,
    pub add_self_aln: bool,
    pub self_: bool,
    pub current_ref_block: BlockId,
    pub current_query_block: BlockId,
    pub lazy_masking: bool,
    pub target_masking: MaskingAlgo,
    pub gapped_filter_evalue: f64,
    pub query_translated: bool,
    pub swipe_all: bool,
    pub lin_stage1_query: bool,
    pub min_bit_score: f64,
    pub max_alignments: i64,
    pub ungapped_window: usize,
}

impl Default for GappedScoreConfig {
    fn default() -> Self {
        Self {
            padding: 0,
            narrow_band_cov: 0.0,
            narrow_band_factor: 1.0,
            min_band_overlap: 0.0,
            anchored_swipe: false,
            query_contexts: 1,
            // C++ default is 240 (`diamond/src/basic/config.cpp:552`).
            // The earlier `i8::MAX` (127) caused 16-bit dispatch to kick in
            // for any ungapped score >= 128, leaving the faster 8-bit kernel
            // unused. Bit-exact, since the 16-bit kernel computes the same
            // score; speed-only.
            cutoff_score_8bit: 240,
            max_swipe_dp: 1_000_000,
            approx_backtrace: false,
            max_hsps: 1,
            inner_culling_overlap: 50.0,
            comp_based_stats_hauser: false,
            comp_based_stats_matrix_adjust: false,
            sensitivity: Sensitivity::Default,
            max_evalue: f64::MAX,
            query_cover: 0.0,
            subject_cover: 0.0,
            query_or_target_cover: 0.0,
            approx_min_id: 0.0,
            cbs_matrix_scale: 1,
            min_id: 0.0,
            no_self_hits: false,
            max_target_seqs: 25,
            toppercent: None,
            no_ranking: false,
            global_ranking_targets: 0,
            ext_chunk_size: 0,
            mapany: false,
            target_hard_cap: None,
            // C++ defaults (`diamond/src/basic/config.cpp:579,581`): the
            // chunked-ranking termination guards. Rust used 0.9 / 0.0 which
            // both terminate chunking later than C++, producing extra
            // (later-evicted) work but identical results — fixed for parity
            // with C++ chunk counts.
            ranking_score_drop_factor: 0.95,
            ranking_cutoff_bitscore: 25.0,
            cluster_threshold_present: false,
            add_self_aln: false,
            self_: false,
            current_ref_block: 0,
            current_query_block: 0,
            lazy_masking: false,
            target_masking: MaskingAlgo::None,
            gapped_filter_evalue: 0.0,
            query_translated: false,
            swipe_all: false,
            lin_stage1_query: false,
            min_bit_score: 0.0,
            max_alignments: i64::MAX,
            ungapped_window: 48,
        }
    }
}

/// Matches C++ `ranking_chunk_size`.
pub fn ranking_chunk_size(
    target_count: i64,
    ref_letters: i64,
    max_target_seqs: i64,
    cfg: &GappedScoreConfig,
) -> i64 {
    const MAX_CHUNK_SIZE: i64 = 400;
    const MIN_CHUNK_SIZE: i64 = 128;
    const MAPANY_CHUNK_SIZE: i64 = 16;

    if cfg.no_ranking || cfg.global_ranking_targets > 0 {
        return target_count;
    }
    if cfg.ext_chunk_size > 0 {
        return cfg.ext_chunk_size;
    }
    if cfg.mapany {
        return MAPANY_CHUNK_SIZE;
    }
    let default_letters = if cfg.sensitivity >= Sensitivity::VerySensitive {
        800.0 * 1.0e6
    } else {
        2.0 * 1.0e9
    };
    let block_mult = ((ref_letters as f64 / default_letters).round() as i64).max(1);
    if cfg.toppercent.is_some() {
        return MIN_CHUNK_SIZE * block_mult;
    }
    let multiple = ((max_target_seqs + 31) / 32) * 32;
    let size = MIN_CHUNK_SIZE.max(multiple.min(MAX_CHUNK_SIZE)) * block_mult;
    cfg.target_hard_cap.map_or(size, |cap| size.min(cap))
}

/// Matches C++ `have_filters`.
pub fn have_filters(cfg: &GappedScoreConfig) -> bool {
    cfg.min_id > 0.0
        || cfg.approx_min_id > 0.0
        || cfg.query_cover > 0.0
        || cfg.subject_cover > 0.0
        || cfg.query_or_target_cover > 0.0
}

/// Matches C++ `first_round_hspv`.
pub fn first_round_hspv(cfg: &GappedScoreConfig, output_hsp_values: HspValues) -> HspValues {
    let mut first_round = HspValues::NONE;
    if cfg.min_id > 0.0 {
        first_round = first_round | HspValues::IDENT | HspValues::LENGTH;
    }
    if cfg.query_cover > 0.0 {
        first_round = first_round | HspValues::QUERY_COORDS;
    }
    if cfg.subject_cover > 0.0 {
        first_round = first_round | HspValues::TARGET_COORDS;
    }
    if cfg.cluster_threshold_present {
        first_round = first_round | output_hsp_values;
    }
    first_round
}

/// Matches C++ `ranking_terminate`.
pub fn ranking_terminate(
    new_hits: bool,
    last_tail_score: i32,
    tail_score: i32,
    targets_processed: i64,
    targets_aligned: i64,
    cfg: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
) -> bool {
    if cfg
        .target_hard_cap
        .is_some_and(|cap| targets_processed >= cap)
    {
        return true;
    }
    if cfg.mapany && cfg.toppercent.is_none() && targets_aligned > 0 {
        return true;
    }
    !new_hits
        && (last_tail_score == 0
            || tail_score as f64 / last_tail_score as f64 <= cfg.ranking_score_drop_factor
            || score_matrix.bitscore(tail_score as f64) < cfg.ranking_cutoff_bitscore)
}

/// Matches C++ `add_self_aln`.
pub fn add_self_aln(cfg: &GappedScoreConfig) -> bool {
    cfg.add_self_aln
        && ((cfg.self_ && cfg.current_ref_block == 0)
            || (!cfg.self_ && cfg.current_query_block == cfg.current_ref_block))
}

/// Matches C++ `lazy_masking`.
pub fn lazy_masking(target_block_ids: &[BlockId], targets: &mut Block, algo: MaskingAlgo) -> usize {
    if algo == MaskingAlgo::None {
        return 0;
    }
    let mut seq = Vec::new();
    let mut n = 0usize;
    for &block_id in target_block_ids {
        if targets.fetch_seq_if_unmasked(block_id as usize, &mut seq) {
            mask_sequence(&mut seq, algo);
            targets.write_masked_seq(block_id as usize, &seq);
            n += 1;
        }
    }
    n
}

/// Matches C++ `extend_chunk`.
pub fn extend_chunk<F1, F2>(
    query_seq: &[Vec<Letter>],
    query_id: Option<&str>,
    source_query_len: i32,
    query_cbs: &[Vec<i8>],
    query_comp: &[f64; TRUE_AA as usize],
    seed_hits: &FlatArray<SeedHit>,
    target_block_ids: &[BlockId],
    target_block: &mut Block,
    stat: &mut Statistics,
    flags: Flags,
    hsp_values: HspValues,
    mode: ExtensionMode,
    cfg: &GappedScoreConfig,
    ungapped_cfg: &UngappedStageConfig,
    score_matrix: &ScoreMatrix,
    cutoff_gapped1_new: F1,
    cutoff_gapped2_new: F2,
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
    gapped_filter_window: i32,
) -> Vec<Target>
where
    F1: Fn(i32, i32) -> i32 + Copy,
    F2: Fn(i32, i32) -> i32 + Copy,
{
    const GAPPED_FILTER_MIN_QLEN: i32 = 85;
    let n = seed_hits.data_size() as i64;
    stat.inc(StatValue::TargetHits2, n);

    if cfg.lazy_masking && cfg.global_ranking_targets == 0 {
        let masked = lazy_masking(target_block_ids, target_block, cfg.target_masking);
        stat.inc(StatValue::MaskedLazy, masked as i64);
    }

    let gf_seed_hits;
    let gf_target_block_ids;
    let (seed_hits, target_block_ids) = if cfg.gapped_filter_evalue > 0.0
        && cfg.global_ranking_targets == 0
        && (!cfg.query_translated || query_seq[0].len() as i32 >= GAPPED_FILTER_MIN_QLEN)
    {
        let query_refs: Vec<&[Letter]> = query_seq.iter().map(Vec::as_slice).collect();
        let query_cbs_refs: Vec<&[i8]> = query_cbs.iter().map(Vec::as_slice).collect();
        let query_cbs_arg = if query_cbs_refs.is_empty() {
            None
        } else {
            Some(query_cbs_refs.as_slice())
        };
        let (filtered_hits, filtered_ids) = gapped_filter_seed_hits(
            &query_refs,
            query_cbs_arg,
            seed_hits,
            target_block_ids,
            target_block.seqs(),
            stat,
            flags,
            cutoff_gapped1_new,
            cutoff_gapped2_new,
            gap_open,
            gap_extend,
            gapped_filter_diag_score,
            cfg.query_translated,
            gapped_filter_window,
            score_matrix,
        );
        gf_seed_hits = filtered_hits;
        gf_target_block_ids = filtered_ids;
        (&gf_seed_hits, gf_target_block_ids.as_slice())
    } else {
        (seed_hits, target_block_ids)
    };

    stat.inc(StatValue::TargetHits3, seed_hits.data_size() as i64);

    let mut targets = ungapped_stage_seed_hits(
        query_seq,
        query_cbs,
        query_comp,
        seed_hits,
        target_block_ids,
        stat,
        target_block,
        mode,
        ungapped_cfg,
        score_matrix,
    );
    align_work_targets(
        &mut targets,
        query_seq,
        query_id,
        query_cbs,
        source_query_len,
        flags,
        hsp_values,
        mode,
        stat,
        cfg,
        score_matrix,
    )
}

/// Matches C++ `extend(BlockId, const Search::Config&, Statistics&, DP::Flags, SeedHitList&)`.
pub fn extend_seed_hit_list<F1, F2>(
    query_id: BlockId,
    query_seq: &[Vec<Letter>],
    query_title: &str,
    source_query_len: i32,
    query_cbs: &[Vec<i8>],
    query_comp: &[f64; TRUE_AA as usize],
    target_block: &mut Block,
    stat: &mut Statistics,
    flags: Flags,
    seed_hit_list: &mut SeedHitList,
    mode: ExtensionMode,
    cfg: &GappedScoreConfig,
    ungapped_cfg: &UngappedStageConfig,
    score_matrix: &ScoreMatrix,
    output_hsp_values: HspValues,
    cutoff_gapped1_new: F1,
    cutoff_gapped2_new: F2,
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
    gapped_filter_window: i32,
) -> Vec<Match>
where
    F1: Fn(i32, i32) -> i32 + Copy,
    F2: Fn(i32, i32) -> i32 + Copy,
{
    const UNIFIED_TARGET_LEN: u32 = 50;
    let query_len = query_seq.first().map_or(0, Vec::len) as u32;
    let self_aln_score = if target_block.has_self_aln() {
        target_block.self_aln_score(query_id as i64)
    } else {
        0.0
    };
    let target_count = seed_hit_list.target_block_ids.len();

    if mode == ExtensionMode::None {
        seed_hit_list.target_scores.sort();
        return seed_only_matches(
            seed_hit_list,
            source_query_len as u32,
            Some(target_block.seqs()),
            None,
            cfg.max_hsps,
            cfg.max_target_seqs,
            cfg.toppercent,
            |score| score_matrix.bitscore(score as f64),
        );
    }

    if target_count == 0 {
        let mut matches = if cfg.swipe_all {
            let mut aligned_targets = full_db_align(
                query_seq,
                query_cbs,
                flags,
                HspValues::NONE,
                stat,
                target_block,
                cfg,
                score_matrix,
            );
            culling_targets(
                &mut aligned_targets,
                false,
                cfg.max_target_seqs,
                cfg.toppercent,
                |score| score_matrix.bitscore(score as f64),
            );
            align_targets_final(
                &mut aligned_targets,
                0,
                query_seq,
                query_title,
                query_cbs,
                source_query_len,
                self_aln_score,
                flags,
                HspValues::NONE,
                false,
                mode,
                stat,
                cfg,
                score_matrix,
                output_hsp_values,
            )
        } else {
            Vec::new()
        };
        if add_self_aln(cfg) && !matches.iter().any(|m| m.target_block_id == query_id) {
            matches.push(Match::self_match(query_id, &query_seq[0]));
        }
        culling_matches(&mut matches, cfg.max_target_seqs, cfg.toppercent, |score| {
            score_matrix.bitscore(score as f64)
        });
        return matches;
    }

    let chunk_size = ranking_chunk_size(
        target_count as i64,
        target_block.seqs().letters() as i64,
        cfg.max_target_seqs,
        cfg,
    )
    .max(1) as usize;

    let mut i0 = 0usize;
    let mut i1 = chunk_size.min(seed_hit_list.target_scores.len());
    if cfg.toppercent.is_none()
        && cfg.min_bit_score == 0.0
        && (i1 - i0) < cfg.max_target_seqs as usize
        && (cfg.ext_chunk_size == 0 || cfg.lin_stage1_query)
    {
        while i1 < seed_hit_list.target_scores.len()
            && score_matrix.evalue(
                seed_hit_list.target_scores[i1].score as i32,
                query_len,
                UNIFIED_TARGET_LEN,
            ) <= cfg.max_evalue
        {
            i1 += 16.min(seed_hit_list.target_scores.len() - i1);
        }
    }

    let first_round = if cfg.anchored_swipe
        || (cfg.sensitivity <= Sensitivity::Shapes30x10
            && output_hsp_values.only(HspValues::COORDS)
            && cfg.toppercent.is_none()
            && cfg.max_target_seqs as usize >= seed_hit_list.target_scores.len())
    {
        HspValues::COORDS
    } else {
        HspValues::NONE
    };
    let first_round_culling = !have_filters(cfg) || cfg.toppercent.is_some();
    let mut tail_score = 0i32;
    let mut seed_hits_chunk = FlatArray::new();
    let mut target_block_ids_chunk = Vec::new();
    let mut matches = Vec::new();

    while i0 < seed_hit_list.target_scores.len() {
        let mut aligned_targets = Vec::new();
        let mut new_hits;
        let mut round_new_hits: bool;

        loop {
            let current_chunk_size = i1 - i0;
            let multi_chunk = current_chunk_size < seed_hit_list.target_scores.len();
            if multi_chunk {
                target_block_ids_chunk.clear();
                seed_hits_chunk.clear();
                let data_size: u64 = seed_hit_list.target_scores[i0..i1]
                    .iter()
                    .map(|s| seed_hit_list.seed_hits.count(s.target as u64))
                    .sum();
                seed_hits_chunk.reserve(current_chunk_size as u64, data_size);
                target_block_ids_chunk.reserve(current_chunk_size);
                for target_score in &seed_hit_list.target_scores[i0..i1] {
                    target_block_ids_chunk
                        .push(seed_hit_list.target_block_ids[target_score.target as usize]);
                    seed_hits_chunk
                        .push_back(seed_hit_list.seed_hits.range(target_score.target as u64));
                }
            }

            let v = extend_chunk(
                query_seq,
                Some(query_title),
                source_query_len,
                query_cbs,
                query_comp,
                if multi_chunk {
                    &seed_hits_chunk
                } else {
                    &seed_hit_list.seed_hits
                },
                if multi_chunk {
                    &target_block_ids_chunk
                } else {
                    &seed_hit_list.target_block_ids
                },
                target_block,
                stat,
                flags,
                first_round,
                mode,
                cfg,
                ungapped_cfg,
                score_matrix,
                cutoff_gapped1_new,
                cutoff_gapped2_new,
                gap_open,
                gap_extend,
                gapped_filter_diag_score,
                gapped_filter_window,
            );

            stat.inc(StatValue::TargetHits4, v.len() as i64);
            new_hits = !v.is_empty();
            round_new_hits = new_hits;
            if multi_chunk {
                new_hits = append_hits_targets(
                    &mut aligned_targets,
                    v,
                    first_round_culling,
                    cfg.max_target_seqs,
                    cfg.toppercent,
                    |score| score_matrix.bitscore(score as f64),
                );
            } else {
                aligned_targets = v;
            }

            i0 = i1;
            i1 += chunk_size.min(seed_hit_list.target_scores.len() - i1);
            let previous_tail_score = tail_score;
            if new_hits && i1 > 0 {
                tail_score = seed_hit_list.target_scores[i1 - 1].score as i32;
            }
            if i0 >= seed_hit_list.target_scores.len()
                || ranking_terminate(
                    new_hits,
                    previous_tail_score,
                    seed_hit_list.target_scores[i1.saturating_sub(1)].score as i32,
                    i1 as i64,
                    aligned_targets.len() as i64,
                    cfg,
                    score_matrix,
                )
            {
                break;
            }
        }

        if cfg.swipe_all {
            aligned_targets = full_db_align(
                query_seq,
                query_cbs,
                flags,
                HspValues::NONE,
                stat,
                target_block,
                cfg,
                score_matrix,
            );
        }

        culling_targets(
            &mut aligned_targets,
            !first_round_culling,
            cfg.max_target_seqs,
            cfg.toppercent,
            |score| score_matrix.bitscore(score as f64),
        );
        stat.inc(StatValue::TargetHits5, aligned_targets.len() as i64);

        let mut round_matches = align_targets_final(
            &mut aligned_targets,
            matches.len() as i64,
            query_seq,
            query_title,
            query_cbs,
            source_query_len,
            self_aln_score,
            flags,
            first_round,
            first_round_culling,
            mode,
            stat,
            cfg,
            score_matrix,
            output_hsp_values,
        );
        matches.append(&mut round_matches);

        if cfg.toppercent.is_some()
            || matches.len() as i64 >= cfg.max_target_seqs
            || i0 >= seed_hit_list.target_scores.len()
            || !round_new_hits
            || (cfg.mapany && !matches.is_empty())
        {
            break;
        }
    }

    if add_self_aln(cfg) && !matches.iter().any(|m| m.target_block_id == query_id) {
        matches.push(Match::self_match(query_id, &query_seq[0]));
    }
    culling_matches(&mut matches, cfg.max_target_seqs, cfg.toppercent, |score| {
        score_matrix.bitscore(score as f64)
    });
    matches
}

/// Matches C++ `extend(BlockId, Search::Hit*, Search::Hit*, ...)`.
pub fn extend<F1, F2, G>(
    query_id: BlockId,
    hits: &mut [Hit],
    query_seq: &[Vec<Letter>],
    query_title: &str,
    source_query_len: i32,
    query_cbs: &[Vec<i8>],
    query_comp: &[f64; TRUE_AA as usize],
    target_block: &mut Block,
    stat: &mut Statistics,
    flags: Flags,
    mode: ExtensionMode,
    cfg: &GappedScoreConfig,
    ungapped_cfg: &UngappedStageConfig,
    score_matrix: &ScoreMatrix,
    output_hsp_values: HspValues,
    cutoff_gapped1_new: F1,
    cutoff_gapped2_new: F2,
    gap_open: i32,
    gap_extend: i32,
    gapped_filter_diag_score: i32,
    gapped_filter_window: i32,
    global_ranking_cb: Option<G>,
) -> Vec<Match>
where
    F1: Fn(i32, i32) -> i32 + Copy,
    F2: Fn(i32, i32) -> i32 + Copy,
    G: FnOnce(BlockId, &SeedHitList) -> Vec<Match>,
{
    let mut list = load_hits(hits, target_block.seqs(), cfg.query_contexts as u32);
    stat.inc(StatValue::TargetHits0, list.target_block_ids.len() as i64);

    let target_count = list.target_block_ids.len() as i64;
    if target_count == 0 && !cfg.swipe_all {
        if mode == ExtensionMode::None {
            return Vec::new();
        }
        if add_self_aln(cfg) {
            return vec![Match::self_match(query_id, &query_seq[0])];
        }
        return Vec::new();
    }

    let chunk_size = ranking_chunk_size(
        target_count,
        target_block.seqs().letters() as i64,
        cfg.max_target_seqs,
        cfg,
    );
    if chunk_size < target_count || cfg.global_ranking_targets > 0 {
        list.target_scores.sort();
        if cfg.global_ranking_targets > 0 {
            if let Some(global_ranking_cb) = global_ranking_cb {
                return global_ranking_cb(query_id, &list);
            }
            return global_ranking::ranking_list(
                query_id,
                &mut list.target_scores,
                &list.target_block_ids,
                &list.seed_hits,
                query_seq.first().map(Vec::as_slice),
                Some(target_block),
                cfg.global_ranking_targets,
                cfg.ungapped_window,
                score_matrix,
            );
        }
    }

    extend_seed_hit_list(
        query_id,
        query_seq,
        query_title,
        source_query_len,
        query_cbs,
        query_comp,
        target_block,
        stat,
        flags,
        &mut list,
        mode,
        cfg,
        ungapped_cfg,
        score_matrix,
        output_hsp_values,
        cutoff_gapped1_new,
        cutoff_gapped2_new,
        gap_open,
        gap_extend,
        gapped_filter_diag_score,
        gapped_filter_window,
    )
}

/// Matches C++ `filter_hspvalues`.
pub fn filter_hspvalues(cfg: &GappedScoreConfig) -> HspValues {
    let mut hsp_values = HspValues::NONE;
    if cfg.max_hsps != 1 {
        hsp_values = hsp_values | HspValues::QUERY_COORDS | HspValues::TARGET_COORDS;
    }
    if cfg.min_id > 0.0 {
        hsp_values = hsp_values | HspValues::IDENT | HspValues::LENGTH;
    }
    if cfg.approx_min_id > 0.0 {
        hsp_values = hsp_values | HspValues::COORDS;
    }
    if cfg.query_cover > 0.0 {
        hsp_values = hsp_values | HspValues::QUERY_COORDS;
    }
    if cfg.subject_cover > 0.0 {
        hsp_values = hsp_values | HspValues::TARGET_COORDS;
    }
    if cfg.query_or_target_cover > 0.0 {
        hsp_values = hsp_values | HspValues::COORDS;
    }
    hsp_values
}

/// Matches C++ `first_round_filter_all`.
pub fn first_round_filter_all(cfg: &GappedScoreConfig, first_round_hsp_values: HspValues) -> bool {
    if cfg.min_id > 0.0 && !first_round_hsp_values.all(HspValues::IDENT | HspValues::LENGTH) {
        return false;
    }
    if cfg.approx_min_id > 0.0 && !first_round_hsp_values.all(HspValues::COORDS) {
        return false;
    }
    if cfg.query_cover > 0.0 && !first_round_hsp_values.all(HspValues::QUERY_COORDS) {
        return false;
    }
    if cfg.subject_cover > 0.0 && !first_round_hsp_values.all(HspValues::TARGET_COORDS) {
        return false;
    }
    if cfg.query_or_target_cover > 0.0 && !first_round_hsp_values.all(HspValues::COORDS) {
        return false;
    }
    true
}

/// Matches C++ `Extension::band`.
pub fn band(len: i32, mode: ExtensionMode, padding: i32) -> i32 {
    if padding > 0 {
        return padding;
    }
    if mode == ExtensionMode::BandedFast {
        if len < 50 {
            12
        } else if len < 100 {
            16
        } else if len < 250 {
            30
        } else if len < 350 {
            40
        } else {
            64
        }
    } else if len < 50 {
        15
    } else if len < 100 {
        20
    } else if len < 150 {
        30
    } else if len < 200 {
        50
    } else if len < 250 {
        60
    } else if len < 350 {
        100
    } else if len < 500 {
        120
    } else {
        150
    }
}

/// Matches C++ `hsp_band`.
pub fn hsp_band(
    base_band: i32,
    qlen: i32,
    tlen: i32,
    hsp: &ApproxHsp,
    narrow_band_cov: f64,
    narrow_band_factor: f64,
) -> i32 {
    if narrow_band_cov == 0.0 {
        return base_band;
    }
    if hsp.query_range.length() as f64 / qlen as f64 >= narrow_band_cov
        || hsp.subject_range.length() as f64 / tlen as f64 >= narrow_band_cov
    {
        ((hsp.d_max - hsp.d_min) as f64 * narrow_band_factor) as i32
    } else {
        base_band
    }
}

/// Matches C++ `add_dp_targets`.
pub fn add_dp_targets(
    target: &crate::align::ungapped::WorkTarget,
    target_idx: BlockId,
    matrix: Option<Arc<TargetMatrix>>,
    query_seq: &[Vec<Letter>],
    dp_targets: &mut [Targets; MAX_CONTEXT as usize],
    hsp_values: HspValues,
    mode: ExtensionMode,
    cfg: &GappedScoreConfig,
) {
    let base_band = band(query_seq[0].len() as i32, mode, cfg.padding);
    let slen = target.seq.len() as i32;
    let score_width = matrix
        .as_deref()
        .map(TargetMatrix::score_width)
        .unwrap_or(0);
    for frame in 0..cfg.query_contexts {
        let qlen = query_seq[frame].len() as i32;
        if mode == ExtensionMode::Full {
            if target.ungapped_score[frame] == 0 {
                continue;
            }
            let b = swipe_bin(
                hsp_values,
                qlen,
                0,
                target.ungapped_score[frame],
                qlen as i64 * slen as i64,
                score_width,
                0,
                cfg.cutoff_score_8bit,
                cfg.max_swipe_dp,
                cfg.approx_backtrace,
            );
            let mut dp = DpTarget::new(
                target.seq.clone(),
                slen,
                0,
                0,
                target_idx as i64,
                qlen,
                CarryOver::default(),
                DpAnchor::default(),
            );
            if let Some(matrix) = matrix.clone() {
                dp = dp.with_matrix(matrix, 1);
            }
            dp_targets[frame][b].push_back(dp);
            continue;
        }
        if target.hsp[frame].is_empty() {
            continue;
        }

        let mut d0 = i32::MAX;
        let mut d1 = i32::MIN;
        let mut score = 0;
        let mut anchor = Anchor::default();

        for hsp in &target.hsp[frame] {
            let b = hsp_band(
                base_band,
                qlen,
                slen,
                hsp,
                cfg.narrow_band_cov,
                cfg.narrow_band_factor,
            );
            let b0 = (hsp.d_min - b).max(-(slen - 1));
            let b1 = (hsp.d_max + 1 + b).min(qlen);
            let overlap = crate::util::interval::Interval::new(d0, d1)
                .intersect(&crate::util::interval::Interval::new(b0, b1))
                .length() as f64;
            if d0 != i32::MAX
                && ((overlap / (d1 - d0) as f64 > cfg.min_band_overlap)
                    || (overlap / (b1 - b0) as f64 > cfg.min_band_overlap))
                && !cfg.anchored_swipe
            {
                d0 = d0.min(b0);
                d1 = d1.max(b1);
                score = score.max(hsp.score);
                if hsp.max_diag.segment.score > anchor.segment.score {
                    anchor = hsp.max_diag.clone();
                }
            } else {
                if d0 != i32::MAX {
                    let dp_size =
                        DpTarget::banded_cols(qlen, slen, d0, d1) as i64 * (d1 - d0) as i64;
                    let bin = swipe_bin(
                        hsp_values,
                        d1 - d0,
                        0,
                        score,
                        dp_size,
                        score_width,
                        0,
                        cfg.cutoff_score_8bit,
                        cfg.max_swipe_dp,
                        cfg.approx_backtrace,
                    );
                    let mut dp = DpTarget::new(
                        target.seq.clone(),
                        slen,
                        d0,
                        d1,
                        target_idx as i64,
                        qlen,
                        CarryOver::default(),
                        DpAnchor {
                            query_begin: anchor.segment.i,
                            query_end: anchor.segment.query_end(),
                            subject_begin: anchor.segment.j,
                            subject_end: anchor.segment.subject_end(),
                            d_min_left: anchor.d_min_left,
                            d_max_left: anchor.d_max_left,
                            d_min_right: anchor.d_min_right,
                            d_max_right: anchor.d_max_right,
                            prefix_score: anchor.prefix_score,
                            score: anchor.segment.score,
                        },
                    );
                    if let Some(matrix) = matrix.clone() {
                        dp = dp.with_matrix(matrix, 1);
                    }
                    dp_targets[frame][bin].push_back(dp);
                }
                d0 = b0;
                d1 = b1;
                score = hsp.score;
                anchor = hsp.max_diag.clone();
            }
        }
        let dp_size = DpTarget::banded_cols(qlen, slen, d0, d1) as i64 * (d1 - d0) as i64;
        let bin = swipe_bin(
            hsp_values,
            d1 - d0,
            0,
            score,
            dp_size,
            score_width,
            0,
            cfg.cutoff_score_8bit,
            cfg.max_swipe_dp,
            cfg.approx_backtrace,
        );
        let mut dp = DpTarget::new(
            target.seq.clone(),
            slen,
            d0,
            d1,
            target_idx as i64,
            qlen,
            CarryOver::default(),
            DpAnchor {
                query_begin: anchor.segment.i,
                query_end: anchor.segment.query_end(),
                subject_begin: anchor.segment.j,
                subject_end: anchor.segment.subject_end(),
                d_min_left: anchor.d_min_left,
                d_max_left: anchor.d_max_left,
                d_min_right: anchor.d_min_right,
                d_max_right: anchor.d_max_right,
                prefix_score: anchor.prefix_score,
                score: anchor.segment.score,
            },
        );
        if let Some(matrix) = matrix.clone() {
            dp = dp.with_matrix(matrix, 1);
        }
        dp_targets[frame][bin].push_back(dp);
    }
}

/// Matches C++ `align(vector<WorkTarget>&, const Sequence*, ...)` in `gapped_score.cpp`.
pub fn align_work_targets(
    targets: &mut Vec<crate::align::ungapped::WorkTarget>,
    query_seq: &[Vec<Letter>],
    query_id: Option<&str>,
    query_cbs: &[Vec<i8>],
    source_query_len: i32,
    mut flags: Flags,
    hsp_values: HspValues,
    mode: ExtensionMode,
    stat: &mut Statistics,
    cfg: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
) -> Vec<Target> {
    let mut dp_targets: [Targets; MAX_CONTEXT as usize] =
        std::array::from_fn(|_| make_dp_targets());
    let mut r = Vec::new();
    if targets.is_empty() {
        return r;
    }
    r.reserve(targets.len());
    let mut cbs_targets = 0;

    for i in 0..targets.len() {
        r.push(Target::new(
            targets[i].block_id,
            &targets[i].seq,
            targets[i].ungapped_score[0],
            targets[i].matrix.clone(),
        ));
        if targets[i].done {
            debug_assert_eq!(targets[i].hsp[0].len(), 1);
            debug_assert_eq!(cfg.query_contexts, 1);
            let frame = targets[i].hsp[0][0].frame as usize;
            r.last_mut().unwrap().add_approx_hit(
                &targets[i].hsp[0][0],
                query_seq[frame].len() as i32,
                score_matrix,
            );
        } else {
            add_dp_targets(
                &targets[i],
                i as BlockId,
                r.last().unwrap().matrix.clone(),
                query_seq,
                &mut dp_targets,
                hsp_values,
                mode,
                cfg,
            );
        }
        if targets[i].matrix.is_some() {
            cbs_targets += 1;
        }
    }
    stat.inc(StatValue::TargetHits3Cbs, cbs_targets);

    if mode == ExtensionMode::Full {
        flags = flags | Flags::FULL_MATRIX;
    }
    if mode == ExtensionMode::Global {
        flags = flags | Flags::SEMI_GLOBAL;
    }

    for frame in 0..cfg.query_contexts {
        let n: i64 = dp_targets[frame].iter().map(|v| v.size()).sum();
        if n == 0 {
            continue;
        }
        let cbs = if cfg.comp_based_stats_hauser {
            query_cbs.get(frame).map(Vec::as_slice)
        } else {
            None
        };
        let mut hsps = if cfg.anchored_swipe {
            let mut acfg = AnchoredSwipeConfig::new(&query_seq[frame], score_matrix);
            acfg.query_cbs = cbs;
            acfg.recompute_adjusted = cfg.comp_based_stats_matrix_adjust;
            acfg.sensitivity = cfg.sensitivity;
            acfg.query_cover = cfg.query_cover;
            acfg.subject_cover = cfg.subject_cover;
            acfg.query_or_target_cover = cfg.query_or_target_cover;
            acfg.max_evalue = cfg.max_evalue;
            anchored_swipe(&mut dp_targets[frame], &acfg)
        } else {
            let mut params = Params::new(&query_seq[frame], score_matrix);
            params.query_id = query_id;
            params.frame = frame as i32;
            params.query_source_len = source_query_len;
            params.composition_bias = cbs;
            params.flags = flags;
            params.v = hsp_values;
            params.cutoff_score_8bit = cfg.cutoff_score_8bit;
            params.max_swipe_dp = cfg.max_swipe_dp;
            params.approx_backtrace = cfg.approx_backtrace;
            params.max_evalue = cfg.max_evalue;
            params.query_cover = cfg.query_cover;
            params.subject_cover = cfg.subject_cover;
            params.query_or_target_cover = cfg.query_or_target_cover;
            params.approx_min_id = cfg.approx_min_id;
            params.cbs_matrix_scale = cfg.cbs_matrix_scale;
            swipe(&dp_targets[frame], &mut params)
        };
        while !hsps.is_empty() {
            let idx = hsps[0].swipe_target as usize;
            r[idx].add_hit_from_vec(&mut hsps, 0);
        }
    }

    let mut r2 = Vec::with_capacity(r.len());
    for mut target in r {
        if target.filter_evalue != f64::MAX {
            if cfg.max_hsps == 1 || have_coords(hsp_values) {
                target.inner_culling(cfg.max_hsps, cfg.inner_culling_overlap, cfg.query_contexts);
            }
            r2.push(target);
        }
    }
    r2
}

/// Matches C++ `full_db_align`.
pub fn full_db_align(
    query_seq: &[Vec<Letter>],
    query_cbs: &[Vec<i8>],
    flags: Flags,
    hsp_values: HspValues,
    stat: &mut Statistics,
    target_block: &Block,
    cfg: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
) -> Vec<Target> {
    let _ = stat;
    let ref_seqs = target_block.seqs();
    let mut hsp = Vec::new();

    for frame in 0..cfg.query_contexts {
        let cbs = if cfg.comp_based_stats_hauser {
            query_cbs.get(frame).map(Vec::as_slice)
        } else {
            None
        };
        let mut params = Params::new(&query_seq[frame], score_matrix);
        params.query_id = Some("");
        params.frame = frame as i32;
        params.query_source_len = query_seq[frame].len() as i32;
        params.composition_bias = cbs;
        params.flags = flags | Flags::FULL_MATRIX;
        params.target_max_len = ref_seqs.max_len(0, ref_seqs.len() as BlockId);
        params.v = hsp_values;
        params.cutoff_score_8bit = cfg.cutoff_score_8bit;
        params.max_swipe_dp = cfg.max_swipe_dp;
        params.approx_backtrace = cfg.approx_backtrace;
        params.max_evalue = cfg.max_evalue;
        params.query_cover = cfg.query_cover;
        params.subject_cover = cfg.subject_cover;
        params.query_or_target_cover = cfg.query_or_target_cover;
        params.approx_min_id = cfg.approx_min_id;
        params.cbs_matrix_scale = cfg.cbs_matrix_scale;
        hsp.append(&mut swipe_set(ref_seqs, &mut params));
    }

    let mut subject_idx = BTreeMap::new();
    let mut r = Vec::new();
    while !hsp.is_empty() {
        let block_id = hsp[0].swipe_target as BlockId;
        let idx = if let Some(&idx) = subject_idx.get(&block_id) {
            idx
        } else {
            let idx = r.len();
            subject_idx.insert(block_id, idx);
            r.push(Target::new(
                block_id,
                ref_seqs.get(block_id as usize),
                0,
                None,
            ));
            idx
        };
        r[idx].add_hit_from_vec(&mut hsp, 0);
    }

    r
}

/// Matches C++ `add_dp_targets` in `gapped_final.cpp`.
pub fn add_final_dp_targets(
    target: &Target,
    target_idx: usize,
    query_seq: &[Vec<Letter>],
    dp_targets: &mut [Targets; MAX_CONTEXT as usize],
    flags: Flags,
    hsp_values: HspValues,
    cfg: &GappedScoreConfig,
) {
    let tlen = target.seq.len() as i32;
    let score_width = target
        .matrix
        .as_deref()
        .map(TargetMatrix::score_width)
        .unwrap_or(0);
    for frame in 0..cfg.query_contexts {
        let qlen = query_seq[frame].len() as i32;
        for hsp in &target.hsp[frame] {
            let dp_size = if flags.any(Flags::FULL_MATRIX) {
                qlen as i64 * tlen as i64
            } else {
                DpTarget::banded_cols(qlen, tlen, hsp.d_begin, hsp.d_end) as i64
                    * (hsp.d_end - hsp.d_begin) as i64
            };
            let bin = swipe_bin(
                hsp_values,
                if flags.any(Flags::FULL_MATRIX) {
                    qlen
                } else {
                    hsp.d_end - hsp.d_begin
                },
                hsp.score,
                0,
                dp_size,
                score_width,
                0,
                cfg.cutoff_score_8bit,
                cfg.max_swipe_dp,
                cfg.approx_backtrace,
            );
            let mut dp = DpTarget::new(
                target.seq.clone(),
                tlen,
                hsp.d_begin,
                hsp.d_end,
                target_idx as i64,
                qlen,
                CarryOver::default(),
                DpAnchor::default(),
            );
            if let Some(matrix) = target.matrix.clone() {
                dp = dp.with_matrix(matrix, cfg.cbs_matrix_scale);
            }
            dp_targets[frame][bin].push_back(dp);
        }
    }
}

/// Matches C++ `align(vector<Target>&, ...)` in `gapped_final.cpp`.
pub fn align_targets_final(
    targets: &mut Vec<Target>,
    previous_matches: i64,
    query_seq: &[Vec<Letter>],
    query_id: &str,
    query_cbs: &[Vec<i8>],
    source_query_len: i32,
    _query_self_aln_score: f64,
    mut flags: Flags,
    first_round: HspValues,
    first_round_culling: bool,
    mode: ExtensionMode,
    stat: &mut Statistics,
    cfg: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
    output_hsp_values: HspValues,
) -> Vec<Match> {
    const MIN_STEP: i64 = 16;
    let mut r = Vec::new();
    if targets.is_empty() {
        return r;
    }

    let mut hsp_values = output_hsp_values;
    let copy_all = cfg.max_hsps == 1
        && first_round.all(hsp_values)
        && first_round_filter_all(cfg, first_round);
    if copy_all {
        r.reserve(targets.len());
    }
    for target in targets.iter_mut() {
        if copy_all || target.done {
            r.push(Match::from_target_hsps(
                target.block_id,
                &target.seq,
                target.matrix.clone(),
                &mut target.hsp,
                target.ungapped_score,
                cfg.query_contexts,
                cfg.max_hsps,
            ));
        }
    }
    if r.len() == targets.len() {
        for m in &mut r {
            let subject_seq = m.seq.clone();
            match_apply_filters(
                m,
                source_query_len as u32,
                query_id,
                &query_seq[0],
                subject_seq.len() as u32,
                None,
                &subject_seq,
                cfg.min_id,
                cfg.approx_min_id,
                cfg.query_cover,
                cfg.subject_cover,
                cfg.query_or_target_cover,
                cfg.no_self_hits,
            );
        }
        return r;
    }

    if mode == ExtensionMode::Full {
        flags = flags | Flags::FULL_MATRIX;
    }
    if mode == ExtensionMode::Global {
        flags = flags | Flags::SEMI_GLOBAL;
    }
    hsp_values = hsp_values | filter_hspvalues(cfg);

    let mut it = 0usize;
    let goon = |matches_len: usize| {
        cfg.toppercent.is_some() || (matches_len as i64 + previous_matches) < cfg.max_target_seqs
    };

    while it < targets.len() && goon(r.len()) {
        let mut dp_targets: [Targets; MAX_CONTEXT as usize] =
            std::array::from_fn(|_| make_dp_targets());
        let remaining = targets.len() - it;
        let step_size = if !first_round_culling && cfg.toppercent.is_none() {
            let want = (cfg.max_target_seqs - r.len() as i64).max(MIN_STEP);
            let step = ((want + MIN_STEP - 1) / MIN_STEP) * MIN_STEP;
            step.min(remaining as i64) as usize
        } else {
            remaining
        };
        r.reserve(step_size);
        let matches_begin = r.len();

        for target in targets.iter().skip(it).take(step_size) {
            if target.done {
                continue;
            }
            add_final_dp_targets(
                target,
                r.len(),
                query_seq,
                &mut dp_targets,
                flags,
                hsp_values,
                cfg,
            );
            r.push(Match::new_extension(
                target.block_id,
                &target.seq,
                target.matrix.clone(),
                target.ungapped_score,
                0,
                f64::MAX,
            ));
        }

        for frame in 0..cfg.query_contexts {
            let n: i64 = dp_targets[frame].iter().map(|v| v.size()).sum();
            if n == 0 {
                continue;
            }
            let cbs = if cfg.comp_based_stats_hauser {
                query_cbs.get(frame).map(Vec::as_slice)
            } else {
                None
            };
            let mut params = Params::new(&query_seq[frame], score_matrix);
            params.query_id = Some(query_id);
            params.frame = frame as i32;
            params.query_source_len = source_query_len;
            params.composition_bias = cbs;
            params.flags = flags;
            params.v = hsp_values;
            params.cutoff_score_8bit = cfg.cutoff_score_8bit;
            params.max_swipe_dp = cfg.max_swipe_dp;
            params.approx_backtrace = cfg.approx_backtrace;
            params.max_evalue = cfg.max_evalue;
            params.query_cover = cfg.query_cover;
            params.subject_cover = cfg.subject_cover;
            params.query_or_target_cover = cfg.query_or_target_cover;
            params.approx_min_id = cfg.approx_min_id;
            params.cbs_matrix_scale = cfg.cbs_matrix_scale;
            let mut hsps = swipe(&dp_targets[frame], &mut params);
            while !hsps.is_empty() {
                let idx = hsps[0].swipe_target as usize;
                r[idx].add_hit_from_vec(&mut hsps, 0);
            }
        }

        for m in r.iter_mut().skip(matches_begin) {
            match_inner_culling(m, cfg.max_hsps, cfg.inner_culling_overlap);
            // C++ `culling.cpp:84-86` reads `hsp.front()` after a sort by
            // `Hsp::operator<` (score DESC, d_begin ASC, qpos ASC).
            // `match_inner_culling` sorts via `less_than_score_position`
            // which is the same comparator, so `hsps[0]` is the C++-front
            // element. Rust's `Iterator::max_by` returns the LAST equal-max
            // element — on score ties this picks the highest d_begin/qpos
            // (worst by C++ ordering) and inverts the tie-break, perturbing
            // downstream culling and output row order.
            if let Some(best) = m.hsps.first() {
                m.filter_score = best.score;
                m.filter_evalue = best.evalue;
            }
            let subject_seq = m.seq.clone();
            match_apply_filters(
                m,
                source_query_len as u32,
                query_id,
                &query_seq[0],
                subject_seq.len() as u32,
                None,
                &subject_seq,
                cfg.min_id,
                cfg.approx_min_id,
                cfg.query_cover,
                cfg.subject_cover,
                cfg.query_or_target_cover,
                cfg.no_self_hits,
            );
        }
        culling_matches(&mut r, cfg.max_target_seqs, cfg.toppercent, |score| {
            score_matrix.bitscore(score as f64)
        });
        stat.inc(StatValue::TargetHits6, step_size as i64);
        it += step_size;
    }

    recompute_alt_hsps(
        &mut r,
        query_seq,
        query_cbs,
        source_query_len,
        hsp_values,
        stat,
        cfg,
        score_matrix,
    );
    r
}

#[derive(Debug, Clone)]
struct ActiveTarget {
    match_index: usize,
    masked_seq: [Option<Vec<Letter>>; MAX_CONTEXT as usize],
    active: u32,
}

impl ActiveTarget {
    /// Matches C++ `ActiveTarget::ActiveTarget`.
    fn new(match_index: usize, m: &Match, query_contexts: usize) -> Self {
        let mut masked_seq: [Option<Vec<Letter>>; MAX_CONTEXT as usize] =
            std::array::from_fn(|_| None);
        let mut reserved = 0u32;
        for h in &m.hsps {
            let bit = 1u32 << h.frame;
            if reserved & bit == 0 && (h.frame as usize) < query_contexts {
                masked_seq[h.frame as usize] = Some(m.seq.clone());
                reserved |= bit;
            }
        }
        Self {
            match_index,
            masked_seq,
            active: 0,
        }
    }

    /// Matches C++ `ActiveTarget::copy_seq`.
    fn copy_seq(&mut self, m: &Match) {
        for h in &m.hsps {
            let frame = h.frame as usize;
            if let Some(seq) = &mut self.masked_seq[frame] {
                let begin = h.subject_range.begin.max(0) as usize;
                let end = h
                    .subject_range
                    .end
                    .max(h.subject_range.begin)
                    .min(seq.len() as i32) as usize;
                seq[begin..end].fill(SUPER_HARD_MASK);
            }
        }
    }

    /// Matches C++ `ActiveTarget::check_fully_masked`.
    fn check_fully_masked(&mut self, query_contexts: usize) -> i32 {
        let mut n = 0;
        for i in 0..query_contexts {
            if self.active & (1u32 << i) != 0 {
                if self.masked_seq[i].as_deref().is_none_or(is_fully_masked) {
                    self.active &= !(1u32 << i);
                } else {
                    n += 1;
                }
            }
        }
        n
    }
}

/// Matches C++ `recompute_alt_hsps(TargetVec&, ...)`.
fn recompute_alt_hsps_round(
    matches: &mut [Match],
    targets: &mut [ActiveTarget],
    query_seq: &[Vec<Letter>],
    query_cbs: &[Vec<i8>],
    query_source_len: i32,
    hsp_values: HspValues,
    cfg: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
) -> Vec<ActiveTarget> {
    let mut dp_targets: [Targets; MAX_CONTEXT as usize] =
        std::array::from_fn(|_| make_dp_targets());
    let qlen = query_seq[0].len() as i32;

    for (target_idx, active) in targets.iter().enumerate() {
        let m = &matches[active.match_index];
        let dp_size = qlen as i64 * m.seq.len() as i64;
        let score_width = m
            .matrix
            .as_deref()
            .map(TargetMatrix::score_width)
            .unwrap_or(0);
        let bin = swipe_bin(
            hsp_values,
            qlen,
            0,
            0,
            dp_size,
            score_width,
            0,
            cfg.cutoff_score_8bit,
            cfg.max_swipe_dp,
            cfg.approx_backtrace,
        );
        for context in 0..cfg.query_contexts {
            if let Some(masked) = &active.masked_seq[context] {
                let mut dp = DpTarget::full(
                    masked.clone(),
                    masked.len() as i32,
                    target_idx as i64,
                    CarryOver::default(),
                );
                if let Some(matrix) = m.matrix.clone() {
                    dp = dp.with_matrix(matrix, cfg.cbs_matrix_scale);
                }
                dp_targets[context][bin].push_back(dp);
            }
        }
    }

    for context in 0..cfg.query_contexts {
        let cbs = if cfg.comp_based_stats_hauser {
            query_cbs.get(context).map(Vec::as_slice)
        } else {
            None
        };
        let mut params = Params::new(&query_seq[context], score_matrix);
        params.query_id = Some("");
        params.frame = context as i32;
        params.query_source_len = query_source_len;
        params.composition_bias = cbs;
        params.flags = Flags::FULL_MATRIX;
        params.v = hsp_values;
        params.cutoff_score_8bit = cfg.cutoff_score_8bit;
        params.max_swipe_dp = cfg.max_swipe_dp;
        params.approx_backtrace = cfg.approx_backtrace;
        params.max_evalue = cfg.max_evalue;
        params.query_cover = cfg.query_cover;
        params.subject_cover = cfg.subject_cover;
        params.query_or_target_cover = cfg.query_or_target_cover;
        params.approx_min_id = cfg.approx_min_id;
        params.cbs_matrix_scale = cfg.cbs_matrix_scale;
        let mut hsps = swipe(&dp_targets[context], &mut params);
        while !hsps.is_empty() {
            let active_idx = hsps[0].swipe_target as usize;
            let match_idx = targets[active_idx].match_index;
            let hsp = hsps.remove(0);
            let begin = hsp.subject_range.begin.max(0) as usize;
            let end = hsp
                .subject_range
                .end
                .max(hsp.subject_range.begin)
                .min(matches[match_idx].seq.len() as i32) as usize;
            matches[match_idx].hsps.push(hsp);
            if let Some(masked) = &mut targets[active_idx].masked_seq[context] {
                masked[begin..end].fill(SUPER_HARD_MASK);
            }
            targets[active_idx].active |= 1u32 << context;
        }
    }

    let mut out = Vec::new();
    for active in targets {
        if active.active != 0 {
            let m = &mut matches[active.match_index];
            match_inner_culling(m, cfg.max_hsps, cfg.inner_culling_overlap);
            // See comment above on `max_by` vs `first()` tie-break ordering.
            if let Some(best) = m.hsps.first() {
                m.filter_score = best.score;
                m.filter_evalue = best.evalue;
            }
            if active.check_fully_masked(cfg.query_contexts) > 0
                && (m.hsps.len() < cfg.max_hsps as usize || cfg.max_hsps == 0)
            {
                let mut next = active.clone();
                for i in 0..cfg.query_contexts {
                    if next.active & (1u32 << i) == 0 {
                        next.masked_seq[i] = None;
                    }
                }
                next.active = 0;
                out.push(next);
            }
        }
    }
    out
}

/// Matches C++ `recompute_alt_hsps(vector<Match>::iterator, ...)`.
pub fn recompute_alt_hsps(
    matches: &mut [Match],
    query_seq: &[Vec<Letter>],
    query_cbs: &[Vec<i8>],
    query_source_len: i32,
    hsp_values: HspValues,
    _stats: &mut Statistics,
    cfg: &GappedScoreConfig,
    score_matrix: &ScoreMatrix,
) {
    if cfg.max_hsps == 1 {
        return;
    }
    let mut targets = Vec::with_capacity(matches.len());
    for i in 0..matches.len() {
        let mut active = ActiveTarget::new(i, &matches[i], cfg.query_contexts);
        active.copy_seq(&matches[i]);
        targets.push(active);
    }
    while !targets.is_empty() {
        targets = recompute_alt_hsps_round(
            matches,
            &mut targets,
            query_seq,
            query_cbs,
            query_source_len,
            hsp_values,
            cfg,
            score_matrix,
        );
    }
}

/// Target culling — remove redundant targets based on overlap.
///
/// When multiple targets align to the same query region, this module
/// decides which targets to keep based on score and coverage overlap.
pub struct TargetCulling {
    /// Overlap fraction threshold for considering targets redundant.
    pub overlap_threshold: f64,
}

impl TargetCulling {
    pub fn new(overlap_threshold: f64) -> Self {
        TargetCulling { overlap_threshold }
    }

    /// Apply culling to a list of matches, keeping only non-redundant ones.
    /// Matches should be pre-sorted by score (best first).
    pub fn cull(&self, matches: &mut Vec<Match>) {
        if matches.is_empty() {
            return;
        }

        // Keep track of which matches to retain
        let mut keep = vec![true; matches.len()];

        for i in 0..matches.len() {
            if !keep[i] {
                continue;
            }
            for j in (i + 1)..matches.len() {
                if !keep[j] {
                    continue;
                }
                // Check if match j is dominated by match i
                if self.is_dominated(&matches[i], &matches[j]) {
                    keep[j] = false;
                }
            }
        }

        let mut write_idx = 0;
        for (read_idx, &should_keep) in keep.iter().enumerate() {
            if should_keep {
                if write_idx != read_idx {
                    matches.swap(write_idx, read_idx);
                }
                write_idx += 1;
            }
        }
        matches.truncate(write_idx);
    }

    /// Check if target b is dominated by target a (a has higher score and overlapping HSPs).
    fn is_dominated(&self, a: &Match, b: &Match) -> bool {
        if a.top_score() <= b.top_score() {
            return false;
        }

        for hsp_b in &b.hsps {
            let mut is_covered = false;
            for hsp_a in &a.hsps {
                let overlap = hsp_a.subject_range.overlap_factor(&hsp_b.subject_range);
                if overlap >= self.overlap_threshold {
                    is_covered = true;
                    break;
                }
            }
            if !is_covered {
                return false;
            }
        }
        true
    }
}

/// Apply max-target-seqs filtering.
pub fn apply_max_target_seqs(matches: &mut Vec<Match>, max_targets: usize) {
    if matches.len() > max_targets {
        matches.truncate(max_targets);
    }
}

/// Apply top-percent filtering: keep matches within `top_percent` of the best score.
pub fn apply_top_percent(matches: &mut Vec<Match>, top_percent: f64) {
    if matches.is_empty() {
        return;
    }
    let top_score = matches[0].top_score() as f64;
    let cutoff = top_score * (1.0 - top_percent / 100.0);
    matches.retain(|m| m.top_score() as f64 >= cutoff);
}

/// Matches C++ `Extension::max_hsp_culling(list<Hsp>&)`.
pub fn max_hsp_culling(hsps: &mut Vec<Hsp>, max_hsps: u32) {
    if max_hsps > 0 && hsps.len() > max_hsps as usize {
        hsps.truncate(max_hsps as usize);
    }
}

/// Matches C++ `Extension::inner_culling(list<Hsp>&)`.
pub fn inner_hsp_culling(hsps: &mut Vec<Hsp>, max_hsps: u32, inner_culling_overlap: f64) {
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| {
        if a.less_than_score_position(b) {
            std::cmp::Ordering::Less
        } else if b.less_than_score_position(a) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });
    if max_hsps == 1 {
        hsps.truncate(1);
        return;
    }
    let overlap = inner_culling_overlap / 100.0;
    let mut i = 0;
    while i < hsps.len() {
        let enveloped = (0..i).any(|j| hsps[i].is_enveloped_by(&hsps[j], overlap));
        if enveloped {
            hsps.remove(i);
        } else {
            i += 1;
        }
    }
    if max_hsps > 0 {
        max_hsp_culling(hsps, max_hsps);
    }
}

/// Matches C++ `Match::inner_culling()`.
pub fn match_inner_culling(m: &mut Match, max_hsps: u32, inner_culling_overlap: f64) {
    inner_hsp_culling(&mut m.hsps, max_hsps, inner_culling_overlap);
}

/// Matches C++ `Match::max_hsp_culling()`.
pub fn match_max_hsp_culling(m: &mut Match, max_hsps: u32) {
    max_hsp_culling(&mut m.hsps, max_hsps);
}

/// Matches C++ `output_range()` for `Match` ranges.
pub fn output_range_matches<F>(
    matches: &[Match],
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> usize
where
    F: FnMut(i32) -> f64,
{
    if matches.is_empty() {
        return 0;
    }
    // Use cached `filter_evalue` / `filter_score` to match C++
    // (`diamond/src/align/culling.cpp:201`, `extend.h:58-63`). After
    // `match_inner_culling` the cached values come from `hsps[0]` (the
    // max-score HSP under `less_than_score_position` ordering); the
    // alternative `top_evalue()`/`top_score()` scan over all HSPs and on
    // CBS-adjusted alignments can pick a lower-score HSP whose evalue is
    // smaller, which perturbs ordering and the empty-row guard relative to
    // C++.
    if matches[0].filter_evalue == f64::MAX {
        return 0;
    }
    if let Some(toppercent) = toppercent {
        let cutoff = ((1.0 - toppercent / 100.0) * bitscore(matches[0].filter_score)).max(1.0);
        let mut i = 0;
        while i < matches.len() && bitscore(matches[i].filter_score) >= cutoff {
            i += 1;
        }
        i
    } else {
        let mut i = (max_target_seqs as usize).min(matches.len());
        while i > 1 && matches[i - 1].filter_evalue == f64::MAX {
            i -= 1;
        }
        i
    }
}

/// Matches C++ `culling(std::vector<Match>&, const Search::Config&)`.
pub fn culling_matches<F>(
    matches: &mut Vec<Match>,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) where
    F: FnMut(i32) -> f64,
{
    // Match C++ `Match::cmp_evalue` (`extend.h:58-63`): full chain is
    // `evalue ASC → score DESC → target_block_id ASC`. The block_id tiebreak
    // is required for deterministic order on equal-evalue-equal-score targets.
    // Sort by cached filter fields to match C++ `Match::cmp_evalue` /
    // `cmp_score` in `extend.h:58-72`. See `output_range_matches` above.
    if toppercent.is_some() {
        matches.sort_by(|a, b| {
            b.filter_score
                .cmp(&a.filter_score)
                .then_with(|| a.target_block_id.cmp(&b.target_block_id))
        });
    } else {
        matches.sort_by(|a, b| {
            a.filter_evalue
                .total_cmp(&b.filter_evalue)
                .then_with(|| b.filter_score.cmp(&a.filter_score))
                .then_with(|| a.target_block_id.cmp(&b.target_block_id))
        });
    }
    let end = output_range_matches(matches, max_target_seqs, toppercent, &mut bitscore);
    matches.truncate(end);
}

/// Matches C++ `culling(std::vector<Match>&, ..., sort_only)`.
pub fn culling_matches_with_sort_only<F>(
    matches: &mut Vec<Match>,
    sort_only: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) where
    F: FnMut(i32) -> f64,
{
    // See `culling_matches`: sort by cached filter fields, not `top_*`.
    if toppercent.is_some() {
        matches.sort_by(|a, b| {
            b.filter_score
                .cmp(&a.filter_score)
                .then_with(|| a.target_block_id.cmp(&b.target_block_id))
        });
    } else {
        matches.sort_by(|a, b| {
            a.filter_evalue
                .total_cmp(&b.filter_evalue)
                .then_with(|| b.filter_score.cmp(&a.filter_score))
                .then_with(|| a.target_block_id.cmp(&b.target_block_id))
        });
    }
    if !sort_only {
        let end = output_range_matches(matches, max_target_seqs, toppercent, &mut bitscore);
        matches.truncate(end);
    }
}

/// Matches C++ `append_hits(vector<Match>&, begin, end, bool, const Search::Config&)`.
pub fn append_hits_matches<F>(
    targets: &mut Vec<Match>,
    mut hits: Vec<Match>,
    with_culling: bool,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    mut bitscore: F,
) -> bool
where
    F: FnMut(i32) -> f64,
{
    if hits.is_empty() {
        return false;
    }
    let mut new_hits = toppercent.is_none() && targets.len() < max_target_seqs as usize;
    let mut append = !with_culling || new_hits;

    culling_matches_with_sort_only(targets, append, max_target_seqs, toppercent, &mut bitscore);

    // Compute the batch's best filter_score / filter_evalue for the
    // "did this batch contribute?" check. Match C++ `append_hits`
    // (`culling.cpp:127-136`), which uses cached filter fields.
    let mut max_score = 0;
    let mut min_evalue = f64::MAX;
    for hit in &hits {
        max_score = max_score.max(hit.filter_score);
        min_evalue = min_evalue.min(hit.filter_evalue);
    }

    let range_end = output_range_matches(targets, max_target_seqs, toppercent, &mut bitscore);
    if targets.is_empty()
        || range_end == 0
        || (toppercent.is_none() && min_evalue <= targets[range_end - 1].filter_evalue)
        || (toppercent.is_some()
            && max_score
                >= ((1.0 - toppercent.unwrap() / 100.0)
                    * targets[range_end - 1].filter_score as f64) as i32)
    {
        append = true;
        new_hits = true;
    }

    if append {
        targets.append(&mut hits);
    }
    new_hits
}

/// Matches C++ `filter_hsp(...)` without optional MCL cluster-threshold evaluation.
pub fn filter_hsp(
    hsp: &Hsp,
    source_query_len: u32,
    query_title: &str,
    subject_len: u32,
    subject_title: Option<&str>,
    query_seq: &[Letter],
    subject_seq: &[Letter],
    min_id: f64,
    approx_min_id: f64,
    query_cover: f64,
    subject_cover: f64,
    query_or_target_cover: f64,
    no_self_hits: bool,
) -> bool {
    let qcov = hsp.query_cover_percent(source_query_len);
    let tcov = hsp.subject_cover_percent(subject_len);
    hsp.id_percent() < min_id
        || (approx_min_id > 0.0 && hsp.approx_id < approx_min_id)
        || qcov < query_cover
        || tcov < subject_cover
        || (qcov < query_or_target_cover && tcov < query_or_target_cover)
        || (no_self_hits && query_seq == subject_seq && Some(query_title) == subject_title)
}

/// Matches C++ `Match::apply_filters(...)` for an already materialized target sequence.
pub fn match_apply_filters(
    m: &mut Match,
    source_query_len: u32,
    query_title: &str,
    query_seq: &[Letter],
    subject_len: u32,
    subject_title: Option<&str>,
    subject_seq: &[Letter],
    min_id: f64,
    approx_min_id: f64,
    query_cover: f64,
    subject_cover: f64,
    query_or_target_cover: f64,
    no_self_hits: bool,
) {
    m.hsps.retain(|hsp| {
        !filter_hsp(
            hsp,
            source_query_len,
            query_title,
            subject_len,
            subject_title,
            query_seq,
            subject_seq,
            min_id,
            approx_min_id,
            query_cover,
            subject_cover,
            query_or_target_cover,
            no_self_hits,
        )
    });
    // See comment above on `max_by` vs `first()` tie-break ordering — this
    // path runs post-`inner_culling`, where `hsps[0]` is the C++-front
    // element after sorting by score DESC, d_begin ASC, qpos ASC.
    if let Some(best) = m.hsps.first() {
        m.filter_score = best.score;
        m.filter_evalue = best.evalue;
    } else {
        m.filter_score = 0;
        m.filter_evalue = f64::MAX;
    }
}

/// Matches C++ `apply_filters(vector<Match>::iterator, vector<Match>::iterator, ...)`.
pub fn apply_filters_matches(
    matches: &mut [Match],
    source_query_len: u32,
    query_title: &str,
    query_seq: &[Letter],
    subject_len: u32,
    subject_title: Option<&str>,
    subject_seq: &[Letter],
    min_id: f64,
    approx_min_id: f64,
    query_cover: f64,
    subject_cover: f64,
    query_or_target_cover: f64,
    no_self_hits: bool,
) {
    if min_id > 0.0
        || approx_min_id > 0.0
        || query_cover > 0.0
        || subject_cover > 0.0
        || query_or_target_cover > 0.0
        || no_self_hits
    {
        for m in matches {
            match_apply_filters(
                m,
                source_query_len,
                query_title,
                query_seq,
                subject_len,
                subject_title,
                subject_seq,
                min_id,
                approx_min_id,
                query_cover,
                subject_cover,
                query_or_target_cover,
                no_self_hits,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::basic::value::{Score, SequenceType};
    use crate::dp::ungapped::DiagonalSegment;
    use crate::stats::score_matrix::ScoreMatrix;
    use crate::util::hsp::Anchor;
    use crate::util::interval::Interval;

    fn make_match(block_id: u32, score: Score, subject_range: (i32, i32)) -> Match {
        let mut m = Match::new(block_id, block_id as u64);
        let mut hsp = Hsp::new();
        hsp.score = score;
        hsp.evalue = 1.0 / score as f64;
        hsp.length = 100;
        hsp.identities = score.min(100);
        hsp.query_source_range = Interval::new(subject_range.0, subject_range.1);
        hsp.subject_range = Interval::new(subject_range.0, subject_range.1);
        m.hsps.push(hsp);
        // Mirror what `match_inner_culling` does in production: cache the
        // best HSP's (score, evalue) into `filter_score`/`filter_evalue`,
        // since downstream culling now reads those fields directly.
        if let Some(best) = m.hsps.first() {
            m.filter_score = best.score;
            m.filter_evalue = best.evalue;
        }
        m
    }

    fn make_target(block_id: BlockId, score: Score, evalue: f64) -> Target {
        let mut t = Target::new(block_id, &[0, 1, 2, 3, 4, 5], score, None);
        let mut hsp = Hsp::new();
        hsp.score = score;
        hsp.evalue = evalue;
        hsp.frame = 0;
        hsp.length = 10;
        hsp.identities = 10;
        hsp.query_source_range = Interval::new(0, 10);
        hsp.query_range = Interval::new(0, 10);
        hsp.subject_range = Interval::new(0, 10);
        t.add_hit(hsp);
        t
    }

    #[test]
    fn test_max_target_seqs() {
        let mut matches = vec![
            make_match(0, 100, (0, 50)),
            make_match(1, 90, (0, 50)),
            make_match(2, 80, (0, 50)),
        ];
        apply_max_target_seqs(&mut matches, 2);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn test_top_percent() {
        let mut matches = vec![
            make_match(0, 100, (0, 50)),
            make_match(1, 95, (0, 50)),
            make_match(2, 50, (0, 50)),
        ];
        apply_top_percent(&mut matches, 10.0);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn test_inner_hsp_culling() {
        let mut a = Hsp::new();
        a.score = 100;
        a.query_source_range = Interval::new(0, 100);
        a.subject_range = Interval::new(0, 100);
        let mut b = Hsp::new();
        b.score = 90;
        b.query_source_range = Interval::new(10, 90);
        b.subject_range = Interval::new(10, 90);
        let mut c = Hsp::new();
        c.score = 80;
        c.query_source_range = Interval::new(200, 260);
        c.subject_range = Interval::new(200, 260);
        let mut hsps = vec![c, b, a];
        inner_hsp_culling(&mut hsps, 0, 70.0);
        assert_eq!(hsps.len(), 2);
        assert_eq!(hsps[0].score, 100);
        assert_eq!(hsps[1].score, 80);
    }

    #[test]
    fn test_culling_matches_evalue_and_toppercent() {
        let mut a = make_match(0, 100, (0, 50));
        a.hsps[0].evalue = 1.0e-20;
        let mut b = make_match(1, 95, (0, 50));
        b.hsps[0].evalue = 1.0e-10;
        let mut c = make_match(2, 50, (0, 50));
        c.hsps[0].evalue = 1.0e-5;
        let mut matches = vec![c.clone(), b.clone(), a.clone()];
        culling_matches(&mut matches, 2, None, |score| score as f64);
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].target_block_id, 0);
        assert_eq!(matches[1].target_block_id, 1);

        let mut matches = vec![c, b, a];
        culling_matches(&mut matches, 10, Some(10.0), |score| score as f64);
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].target_block_id, 0);
        assert_eq!(matches[1].target_block_id, 1);
    }

    #[test]
    fn test_append_hits_matches() {
        let mut a = make_match(0, 100, (0, 50));
        a.hsps[0].evalue = 1.0e-20;
        let mut b = make_match(1, 90, (0, 50));
        b.hsps[0].evalue = 1.0e-10;
        let mut targets = vec![a, b];
        let mut c = make_match(2, 80, (0, 50));
        c.hsps[0].evalue = 1.0e-5;
        assert!(!append_hits_matches(
            &mut targets,
            vec![c],
            true,
            2,
            None,
            |score| score as f64,
        ));
        assert_eq!(targets.len(), 2);

        let mut d = make_match(3, 110, (0, 50));
        d.hsps[0].evalue = 1.0e-30;
        assert!(append_hits_matches(
            &mut targets,
            vec![d],
            true,
            2,
            None,
            |score| score as f64,
        ));
        assert_eq!(targets.len(), 3);
    }

    #[test]
    fn test_filter_hsp_thresholds_and_self_hit() {
        let mut hsp = Hsp::new();
        hsp.score = 80;
        hsp.length = 100;
        hsp.identities = 80;
        hsp.approx_id = 75.0;
        hsp.query_source_range = Interval::new(0, 50);
        hsp.subject_range = Interval::new(0, 60);

        assert!(!filter_hsp(
            &hsp,
            100,
            "q",
            100,
            Some("s"),
            &[1, 2, 3],
            &[1, 2, 4],
            70.0,
            70.0,
            40.0,
            50.0,
            40.0,
            false,
        ));
        assert!(filter_hsp(
            &hsp,
            100,
            "q",
            100,
            Some("s"),
            &[1, 2, 3],
            &[1, 2, 4],
            90.0,
            70.0,
            40.0,
            50.0,
            40.0,
            false,
        ));
        assert!(filter_hsp(
            &hsp,
            100,
            "q",
            100,
            Some("q"),
            &[1, 2, 3],
            &[1, 2, 3],
            70.0,
            70.0,
            40.0,
            50.0,
            40.0,
            true,
        ));
    }

    #[test]
    fn test_apply_filters_matches() {
        let mut keep = make_match(0, 90, (0, 100));
        keep.hsps[0].length = 100;
        keep.hsps[0].identities = 90;
        let mut drop = make_match(1, 50, (0, 100));
        drop.hsps[0].length = 100;
        drop.hsps[0].identities = 50;
        let mut matches = vec![keep, drop];
        apply_filters_matches(
            &mut matches,
            100,
            "q",
            &[0, 1, 2],
            100,
            Some("s"),
            &[0, 1, 3],
            80.0,
            0.0,
            0.0,
            0.0,
            0.0,
            false,
        );
        assert_eq!(matches[0].hsps.len(), 1);
        assert!(matches[1].hsps.is_empty());
    }

    #[test]
    fn test_target_add_hit_and_inner_culling() {
        let mut target = Target::new(7, &[0, 1, 2, 3], 42, None);
        let mut h1 = Hsp::new();
        h1.score = 70;
        h1.evalue = 1.0e-5;
        h1.frame = 0;
        h1.query_source_range = Interval::new(0, 20);
        h1.query_range = Interval::new(0, 20);
        h1.subject_range = Interval::new(0, 20);
        let mut h2 = h1.clone();
        h2.score = 90;
        h2.evalue = 1.0e-9;
        h2.frame = 1;
        target.add_hit(h1);
        target.add_hit(h2);
        assert_eq!(target.filter_score, 90);
        assert_eq!(target.best_context, 1);

        target.inner_culling(1, 70.0, 2);
        assert!(target.hsp[0].is_empty());
        assert_eq!(target.hsp[1].len(), 1);
        assert_eq!(target.hsp[1][0].score, 90);
    }

    #[test]
    fn test_hsp_from_approx_and_target_add_approx_hit() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let approx = ApproxHsp::from_parts(
            -2,
            3,
            55,
            0,
            Interval::new(1, 6),
            Interval::new(2, 7),
            Anchor::new(DiagonalSegment::new(1, 2, 5, 55), -2, 3, -2, 3, 55),
            1.0e-6,
        );
        let hsp = Hsp::from_approx_hsp(&approx, 10, 12, &sm);
        assert!(hsp.backtraced);
        assert_eq!(hsp.score, 55);
        assert_eq!(hsp.length, 5);
        assert_eq!(hsp.identities, 5);
        assert_eq!(hsp.subject_source_range, Interval::new(2, 7));

        let mut target = Target::new(3, &[0, 1, 2, 3, 4, 5, 6], 0, None);
        target.add_approx_hit(&approx, 10, &sm);
        assert!(target.done);
        assert_eq!(target.hsp[0].len(), 1);
        assert_eq!(target.filter_score, 55);
    }

    #[test]
    fn test_culling_and_append_targets() {
        let mut targets = vec![
            make_target(2, 70, 1.0e-4),
            make_target(0, 100, 1.0e-20),
            make_target(1, 90, 1.0e-10),
        ];
        culling_targets(&mut targets, false, 2, None, |score| score as f64);
        assert_eq!(targets.len(), 2);
        assert_eq!(targets[0].block_id, 0);
        assert_eq!(targets[1].block_id, 1);

        let mut low = make_target(4, 60, 1.0e-3);
        assert!(!append_hits_targets(
            &mut targets,
            vec![low.clone()],
            true,
            2,
            None,
            |score| score as f64,
        ));
        low.filter_score = 120;
        low.filter_evalue = 1.0e-30;
        assert!(append_hits_targets(
            &mut targets,
            vec![low],
            true,
            2,
            None,
            |score| score as f64,
        ));
    }

    #[test]
    fn test_band_and_hsp_band_thresholds() {
        assert_eq!(band(20, ExtensionMode::BandedFast, 0), 12);
        assert_eq!(band(120, ExtensionMode::BandedFast, 0), 30);
        assert_eq!(band(120, ExtensionMode::BandedSlow, 0), 30);
        assert_eq!(band(700, ExtensionMode::BandedSlow, 0), 150);
        assert_eq!(band(700, ExtensionMode::BandedSlow, 23), 23);

        let h = ApproxHsp::from_parts(
            -5,
            5,
            40,
            0,
            Interval::new(0, 80),
            Interval::new(0, 20),
            Anchor::default(),
            f64::MAX,
        );
        assert_eq!(hsp_band(50, 100, 100, &h, 0.0, 0.5), 50);
        assert_eq!(hsp_band(50, 100, 100, &h, 0.7, 0.5), 5);
    }

    #[test]
    fn test_add_dp_targets_full_and_banded() {
        let mut work = crate::align::ungapped::WorkTarget {
            block_id: 0,
            seq: vec![0, 1, 2, 3, 4, 5, 6, 7],
            ungapped_score: [0; MAX_CONTEXT as usize],
            hsp: std::array::from_fn(|_| Vec::new()),
            matrix: None,
            done: false,
        };
        work.ungapped_score[0] = 33;
        let mut dp_targets: [Targets; MAX_CONTEXT as usize] =
            std::array::from_fn(|_| crate::dp::swipe::targets());
        add_dp_targets(
            &work,
            9,
            None,
            &[vec![0, 1, 2, 3, 4, 5, 6, 7]],
            &mut dp_targets,
            HspValues::NONE,
            ExtensionMode::Full,
            &GappedScoreConfig::default(),
        );
        let total: i64 = dp_targets[0].iter().map(|v| v.size()).sum();
        assert_eq!(total, 1);
        assert_eq!(dp_targets[0][0].as_slice()[0].target_idx, 9);

        let mut work = crate::align::ungapped::WorkTarget {
            block_id: 0,
            seq: vec![0, 1, 2, 3, 4, 5, 6, 7],
            ungapped_score: [0; MAX_CONTEXT as usize],
            hsp: std::array::from_fn(|_| Vec::new()),
            matrix: None,
            done: false,
        };
        work.hsp[0].push(ApproxHsp::from_parts(
            -1,
            1,
            45,
            0,
            Interval::new(1, 5),
            Interval::new(1, 5),
            Anchor::new(DiagonalSegment::new(1, 1, 4, 45), -1, 1, -1, 1, 45),
            f64::MAX,
        ));
        let mut dp_targets: [Targets; MAX_CONTEXT as usize] =
            std::array::from_fn(|_| crate::dp::swipe::targets());
        add_dp_targets(
            &work,
            4,
            None,
            &[vec![0, 1, 2, 3, 4, 5, 6, 7]],
            &mut dp_targets,
            HspValues::NONE,
            ExtensionMode::BandedFast,
            &GappedScoreConfig::default(),
        );
        let total: i64 = dp_targets[0].iter().map(|v| v.size()).sum();
        assert_eq!(total, 1);
        let target = dp_targets[0]
            .iter()
            .find_map(|v| v.as_slice().first())
            .unwrap();
        assert_eq!(target.target_idx, 4);
        assert!(target.d_begin <= -1);
        assert!(target.d_end >= 2);
        assert_eq!(target.anchor.score, 45);
    }

    #[test]
    fn test_align_work_targets_done_and_swipe_paths() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut done = crate::align::ungapped::WorkTarget {
            block_id: 3,
            seq: query.clone(),
            ungapped_score: [0; MAX_CONTEXT as usize],
            hsp: std::array::from_fn(|_| Vec::new()),
            matrix: None,
            done: true,
        };
        done.hsp[0].push(ApproxHsp::from_parts(
            0,
            0,
            50,
            0,
            Interval::new(0, 8),
            Interval::new(0, 8),
            Anchor::new(DiagonalSegment::new(0, 0, 8, 50), 0, 0, 0, 0, 50),
            1.0e-9,
        ));
        let mut targets = vec![done];
        let mut stat = Statistics::new();
        let out = align_work_targets(
            &mut targets,
            std::slice::from_ref(&query),
            Some("q"),
            &[],
            query.len() as i32,
            Flags::NONE,
            HspValues::COORDS,
            ExtensionMode::BandedFast,
            &mut stat,
            &GappedScoreConfig::default(),
            &sm,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].block_id, 3);
        assert_eq!(out[0].filter_score, 50);
        assert_eq!(out[0].hsp[0].len(), 1);

        let mut work = crate::align::ungapped::WorkTarget {
            block_id: 4,
            seq: query.clone(),
            ungapped_score: [0; MAX_CONTEXT as usize],
            hsp: std::array::from_fn(|_| Vec::new()),
            matrix: None,
            done: false,
        };
        work.hsp[0].push(ApproxHsp::from_parts(
            0,
            0,
            40,
            0,
            Interval::new(0, 8),
            Interval::new(0, 8),
            Anchor::new(DiagonalSegment::new(0, 0, 8, 40), 0, 0, 0, 0, 40),
            f64::MAX,
        ));
        let mut targets = vec![work];
        let mut stat = Statistics::new();
        let out = align_work_targets(
            &mut targets,
            std::slice::from_ref(&query),
            Some("q"),
            &[],
            query.len() as i32,
            Flags::NONE,
            HspValues::COORDS,
            ExtensionMode::BandedFast,
            &mut stat,
            &GappedScoreConfig::default(),
            &sm,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].block_id, 4);
        assert!(out[0].filter_score > 0);
        assert_eq!(out[0].hsp[0].len(), 1);
    }

    #[test]
    fn test_full_db_align_groups_hsps_by_target() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut block = Block::new();
        block
            .push_back(&query, None, None, 0, SequenceType::AminoAcid, 0, false)
            .unwrap();
        block
            .push_back(
                &[7, 6, 5, 4, 3, 2, 1, 0],
                None,
                None,
                1,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();

        let mut cfg = GappedScoreConfig::default();
        cfg.query_contexts = 1;
        let mut stat = Statistics::new();
        let out = full_db_align(
            std::slice::from_ref(&query),
            &[],
            Flags::NONE,
            HspValues::COORDS,
            &mut stat,
            &block,
            &cfg,
            &sm,
        );

        assert!(!out.is_empty());
        assert!(out.iter().any(|target| target.block_id == 0));
        for target in &out {
            let hsp_count: usize = target.hsp.iter().map(Vec::len).sum();
            assert!(hsp_count > 0);
        }
    }

    #[test]
    fn test_filter_hspvalues_and_first_round_filter_all() {
        let mut cfg = GappedScoreConfig::default();
        cfg.max_hsps = 2;
        cfg.min_id = 90.0;
        cfg.query_cover = 80.0;
        cfg.subject_cover = 70.0;
        let v = filter_hspvalues(&cfg);
        assert!(v.all(HspValues::QUERY_COORDS));
        assert!(v.all(HspValues::TARGET_COORDS));
        assert!(v.all(HspValues::IDENT | HspValues::LENGTH));
        assert!(!first_round_filter_all(&cfg, HspValues::COORDS));
        assert!(first_round_filter_all(
            &cfg,
            HspValues::COORDS | HspValues::IDENT | HspValues::LENGTH
        ));
    }

    #[test]
    fn test_add_final_dp_targets_and_align_targets_final() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut target = Target::new(11, &query, 50, None);
        let mut hsp = Hsp::new();
        hsp.score = 50;
        hsp.evalue = 1.0e-8;
        hsp.frame = 0;
        hsp.d_begin = 0;
        hsp.d_end = 1;
        hsp.query_range = Interval::new(0, 8);
        hsp.query_source_range = Interval::new(0, 8);
        hsp.subject_range = Interval::new(0, 8);
        target.add_hit(hsp);

        let mut cfg = GappedScoreConfig::default();
        cfg.max_target_seqs = 10;
        let mut dp_targets: [Targets; MAX_CONTEXT as usize] =
            std::array::from_fn(|_| crate::dp::swipe::targets());
        add_final_dp_targets(
            &target,
            0,
            std::slice::from_ref(&query),
            &mut dp_targets,
            Flags::NONE,
            HspValues::COORDS,
            &cfg,
        );
        let total: i64 = dp_targets[0].iter().map(|v| v.size()).sum();
        assert_eq!(total, 1);

        let mut targets = vec![target];
        let mut stat = Statistics::new();
        let out = align_targets_final(
            &mut targets,
            0,
            std::slice::from_ref(&query),
            "q",
            &[],
            query.len() as i32,
            0.0,
            Flags::NONE,
            HspValues::NONE,
            false,
            ExtensionMode::BandedFast,
            &mut stat,
            &cfg,
            &sm,
            HspValues::COORDS,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].target_block_id, 11);
        assert!(!out[0].hsps.is_empty());
        assert!(out[0].filter_score > 0);
    }

    #[test]
    fn test_recompute_alt_hsps_finds_unmasked_copy() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3];
        let subject = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let mut m = Match::new_extension(5, &subject, None, 0, 40, 1.0e-5);
        let mut hsp = Hsp::new();
        hsp.score = 40;
        hsp.evalue = 1.0e-5;
        hsp.frame = 0;
        hsp.query_range = Interval::new(0, 4);
        hsp.query_source_range = Interval::new(0, 4);
        hsp.subject_range = Interval::new(0, 4);
        m.hsps.push(hsp);
        let mut matches = vec![m];
        let mut cfg = GappedScoreConfig::default();
        cfg.max_hsps = 2;
        let mut stat = Statistics::new();
        recompute_alt_hsps(
            &mut matches,
            std::slice::from_ref(&query),
            &[],
            query.len() as i32,
            HspValues::COORDS,
            &mut stat,
            &cfg,
            &sm,
        );
        assert!(matches[0].hsps.len() >= 2);
        assert!(matches[0].hsps.iter().any(|h| h.subject_range.begin >= 4));
    }

    #[test]
    fn test_extend_decision_helpers_match_cpp_rules() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let mut cfg = GappedScoreConfig::default();
        assert_eq!(ranking_chunk_size(1_000, 10_000, 25, &cfg), 128);
        cfg.ext_chunk_size = 33;
        assert_eq!(ranking_chunk_size(1_000, 10_000, 25, &cfg), 33);
        cfg.ext_chunk_size = 0;
        cfg.mapany = true;
        assert_eq!(ranking_chunk_size(1_000, 10_000, 25, &cfg), 16);
        cfg.mapany = false;
        cfg.toppercent = Some(10.0);
        assert_eq!(ranking_chunk_size(1_000, 10_000, 25, &cfg), 128);

        assert!(!have_filters(&cfg));
        cfg.min_id = 90.0;
        assert!(have_filters(&cfg));
        assert!(first_round_hspv(&cfg, HspValues::COORDS).all(HspValues::IDENT | HspValues::LENGTH));
        cfg.min_id = 0.0;
        cfg.cluster_threshold_present = true;
        assert!(first_round_hspv(&cfg, HspValues::TRANSCRIPT).all(HspValues::TRANSCRIPT));

        cfg.toppercent = None;
        cfg.target_hard_cap = Some(10);
        assert!(ranking_terminate(true, 100, 99, 10, 0, &cfg, &sm));
        cfg.target_hard_cap = None;
        cfg.mapany = true;
        assert!(ranking_terminate(true, 100, 99, 1, 1, &cfg, &sm));
        cfg.mapany = false;
        cfg.ranking_score_drop_factor = 0.9;
        assert!(ranking_terminate(false, 100, 80, 1, 0, &cfg, &sm));

        cfg.add_self_aln = true;
        cfg.self_ = false;
        cfg.current_query_block = 2;
        cfg.current_ref_block = 2;
        assert!(add_self_aln(&cfg));
        cfg.current_ref_block = 1;
        assert!(!add_self_aln(&cfg));

        let self_match = Match::self_match(7, &[0, 1, 2]);
        assert_eq!(self_match.target_block_id, 7);
        assert_eq!(self_match.filter_score, i32::MAX);
        assert_eq!(self_match.hsps[0].subject_range, Interval::new(0, 3));
    }

    #[test]
    fn test_lazy_masking_counts_unmasked_targets() {
        let mut block = Block::new();
        block
            .push_back(&[0, 1, 2], None, None, 0, SequenceType::AminoAcid, 0, false)
            .unwrap();
        block
            .push_back(&[3, 4, 5], None, None, 1, SequenceType::AminoAcid, 0, false)
            .unwrap();
        assert_eq!(lazy_masking(&[0, 1], &mut block, MaskingAlgo::None), 0);
        assert_eq!(lazy_masking(&[0, 1], &mut block, MaskingAlgo::Seg), 2);
        assert_eq!(lazy_masking(&[0, 1], &mut block, MaskingAlgo::Seg), 0);
    }

    #[test]
    fn test_extend_chunk_runs_ungapped_and_alignment_stages() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut block = Block::new();
        block
            .push_back(&query, None, None, 0, SequenceType::AminoAcid, 0, false)
            .unwrap();
        let mut seed_hits = FlatArray::new();
        seed_hits.push_back(&[SeedHit::new(3, 3, 40, 0)]);
        let target_ids = vec![0];
        let query_comp = crate::stats::cbs::compute_composition(&query);
        let mut stat = Statistics::new();
        let out = extend_chunk(
            std::slice::from_ref(&query),
            Some("q"),
            query.len() as i32,
            &[],
            &query_comp,
            &seed_hits,
            &target_ids,
            &mut block,
            &mut stat,
            Flags::NONE,
            HspValues::COORDS,
            ExtensionMode::BandedFast,
            &GappedScoreConfig::default(),
            &UngappedStageConfig::default(),
            &sm,
            |_, _| -1,
            |_, _| -1,
            sm.gap_open(),
            sm.gap_extend(),
            1,
            100,
        );
        assert_eq!(stat.get(StatValue::TargetHits2), 1);
        assert_eq!(stat.get(StatValue::TargetHits3), 1);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].block_id, 0);
        assert!(out[0].filter_score > 0);
    }

    #[test]
    fn test_extend_seed_hit_list_seed_only_and_gapped_paths() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut block = Block::new();
        block
            .push_back(&query, None, None, 0, SequenceType::AminoAcid, 0, false)
            .unwrap();

        let mut list = SeedHitList::default();
        list.seed_hits.push_back(&[SeedHit::new(3, 3, 40, 0)]);
        list.target_block_ids.push(0);
        list.target_scores
            .push(crate::align::gapped_filter::TargetScore {
                target: 0,
                score: 40,
            });

        let query_comp = crate::stats::cbs::compute_composition(&query);
        let mut stat = Statistics::new();
        let seed_only = extend_seed_hit_list(
            0,
            std::slice::from_ref(&query),
            "q",
            query.len() as i32,
            &[],
            &query_comp,
            &mut block,
            &mut stat,
            Flags::NONE,
            &mut list.clone(),
            ExtensionMode::None,
            &GappedScoreConfig::default(),
            &UngappedStageConfig::default(),
            &sm,
            HspValues::COORDS,
            |_, _| -1,
            |_, _| -1,
            sm.gap_open(),
            sm.gap_extend(),
            1,
            100,
        );
        assert_eq!(seed_only.len(), 1);
        assert_eq!(seed_only[0].seq, query);
        assert_eq!(seed_only[0].filter_score, 40);
        assert!(seed_only[0].hsps[0].seed_only);

        let mut stat = Statistics::new();
        let gapped = extend_seed_hit_list(
            0,
            std::slice::from_ref(&query),
            "q",
            query.len() as i32,
            &[],
            &query_comp,
            &mut block,
            &mut stat,
            Flags::NONE,
            &mut list,
            ExtensionMode::BandedFast,
            &GappedScoreConfig::default(),
            &UngappedStageConfig::default(),
            &sm,
            HspValues::COORDS,
            |_, _| -1,
            |_, _| -1,
            sm.gap_open(),
            sm.gap_extend(),
            1,
            100,
        );
        assert_eq!(stat.get(StatValue::TargetHits4), 1);
        assert_eq!(stat.get(StatValue::TargetHits5), 1);
        assert_eq!(gapped.len(), 1);
        assert_eq!(gapped[0].target_block_id, 0);
        assert!(gapped[0].filter_score > 0);

        let subject_pos = block.seqs().position(0, 3) as u64;
        let mut hits = vec![Hit::with_score(0, subject_pos, 3, 40)];
        let mut stat = Statistics::new();
        let from_hits = extend(
            0,
            &mut hits,
            std::slice::from_ref(&query),
            "q",
            query.len() as i32,
            &[],
            &query_comp,
            &mut block,
            &mut stat,
            Flags::NONE,
            ExtensionMode::BandedFast,
            &GappedScoreConfig::default(),
            &UngappedStageConfig::default(),
            &sm,
            HspValues::COORDS,
            |_, _| -1,
            |_, _| -1,
            sm.gap_open(),
            sm.gap_extend(),
            1,
            100,
            Option::<fn(BlockId, &SeedHitList) -> Vec<Match>>::None,
        );
        assert_eq!(stat.get(StatValue::TargetHits0), 1);
        assert_eq!(from_hits.len(), 1);
        assert_eq!(from_hits[0].target_block_id, 0);

        let mut hits = vec![Hit::with_score(0, subject_pos, 3, 40)];
        let mut cfg = GappedScoreConfig::default();
        cfg.global_ranking_targets = 1;
        let mut stat = Statistics::new();
        let ranked = extend(
            0,
            &mut hits,
            std::slice::from_ref(&query),
            "q",
            query.len() as i32,
            &[],
            &query_comp,
            &mut block,
            &mut stat,
            Flags::NONE,
            ExtensionMode::Full,
            &cfg,
            &UngappedStageConfig::default(),
            &sm,
            HspValues::COORDS,
            |_, _| -1,
            |_, _| -1,
            sm.gap_open(),
            sm.gap_extend(),
            1,
            100,
            Option::<fn(BlockId, &SeedHitList) -> Vec<Match>>::None,
        );
        assert_eq!(ranked.len(), 1);
        assert_eq!(ranked[0].target_block_id, 0);
        assert_eq!(ranked[0].ungapped_score, 40);
    }
}
