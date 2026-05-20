use crate::basic::packed_loc::PackedLoc;
use crate::basic::value::{letter_mask, Letter, DELIMITER_LETTER, LETTER_MASK};
use crate::data::sequence_set::SequenceSet;
use crate::search::kmer_ranking::{KmerRanking, PackedLocId};
use crate::util::math::bit_length;

pub type Container = Vec<[Letter; 48]>;

pub trait SeedLoc: Copy {
    fn loc(self) -> usize;
    fn block_id(self, seqs: &SequenceSet) -> usize;
}

impl SeedLoc for usize {
    fn loc(self) -> usize {
        self
    }

    fn block_id(self, seqs: &SequenceSet) -> usize {
        seqs.local_position(self).0
    }
}

impl SeedLoc for PackedLoc {
    fn loc(self) -> usize {
        self.as_u64() as usize
    }

    fn block_id(self, seqs: &SequenceSet) -> usize {
        seqs.local_position(self.as_u64() as usize).0
    }
}

impl SeedLoc for PackedLocId {
    fn loc(self) -> usize {
        self.loc as usize
    }

    fn block_id(self, _seqs: &SequenceSet) -> usize {
        self.block_id as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage1KernelKind {
    Stage1,
    Stage1Self,
    Stage1QueryLin,
    Stage1QueryLinRanked,
    Stage1TargetLin,
    Stage1LongestComboLin,
    Stage1MutualCov,
    Stage1SelfMutualCov,
    Stage1MutualCovQueryLin,
    Stage1MutualCovTargetLin,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Stage1DispatchConfig {
    pub lin_stage1_combo: bool,
    pub lin_stage1_query: bool,
    pub lin_stage1_target: bool,
    pub min_length_ratio: f64,
    pub self_search: bool,
    pub current_ref_block: u32,
    pub global_ranking_targets: bool,
    pub hit_keep_target_id: bool,
}

/// Matches C++ `stage1_dispatch(const Search::Config*, PackedLocId)`.
pub fn stage1_dispatch_packed_loc_id(cfg: &Stage1DispatchConfig) -> Stage1KernelKind {
    if cfg.lin_stage1_combo {
        return Stage1KernelKind::Stage1LongestComboLin;
    }
    if cfg.lin_stage1_query {
        return if cfg.min_length_ratio > 0.0 {
            Stage1KernelKind::Stage1MutualCovQueryLin
        } else {
            Stage1KernelKind::Stage1QueryLinRanked
        };
    }
    if cfg.lin_stage1_target {
        return if cfg.min_length_ratio > 0.0 {
            Stage1KernelKind::Stage1MutualCovTargetLin
        } else {
            Stage1KernelKind::Stage1TargetLin
        };
    }
    if cfg.min_length_ratio > 0.0 {
        return if cfg.self_search && cfg.current_ref_block == 0 {
            Stage1KernelKind::Stage1SelfMutualCov
        } else {
            Stage1KernelKind::Stage1MutualCov
        };
    }
    if cfg.self_search && cfg.current_ref_block == 0 {
        return Stage1KernelKind::Stage1Self;
    }
    Stage1KernelKind::Stage1
}

/// Matches C++ `stage1_dispatch(const Search::Config*, PackedLoc)`.
pub fn stage1_dispatch_packed_loc(cfg: &Stage1DispatchConfig) -> Stage1KernelKind {
    if cfg.lin_stage1_query {
        Stage1KernelKind::Stage1QueryLin
    } else if cfg.lin_stage1_target {
        Stage1KernelKind::Stage1TargetLin
    } else if cfg.self_search && cfg.current_ref_block == 0 {
        Stage1KernelKind::Stage1Self
    } else {
        Stage1KernelKind::Stage1
    }
}

/// Matches C++ `keep_target_id(const Search::Config*)`.
pub fn keep_target_id(cfg: &Stage1DispatchConfig) -> bool {
    cfg.hit_keep_target_id
        || cfg.min_length_ratio != 0.0
        || cfg.global_ranking_targets
        || (cfg.self_search && cfg.current_ref_block == 0)
        || cfg.lin_stage1_combo
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct HitField {
    shift: usize,
    word_shift: usize,
    words_per_query: usize,
    data: Vec<u64>,
    hits: Vec<u32>,
}

impl HitField {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn init(&mut self, query_count: usize, max_target: usize) {
        self.shift = bit_length(max_target as u64).max(8) as usize;
        self.word_shift = self.shift - 6;
        self.words_per_query = 1usize << self.word_shift;
        self.data.clear();
        self.data.resize(query_count * self.words_per_query, 0);
        self.hits.clear();
    }

    pub fn set(&mut self, query: u32, target: u32, v: bool) {
        let bit = ((query as usize) << self.shift) | target as usize;
        let w = bit >> 6;
        let m = 1u64 << (bit & 63);
        let word = &mut self.data[w];
        *word ^= ((0u64.wrapping_sub(v as u64)) ^ *word) & m;
    }

    pub fn hits(&mut self, query: usize) -> &[u32] {
        self.hits.clear();
        let base_w = query << self.word_shift;
        for off in 0..self.words_per_query {
            let mut w = self.data[base_w + off];
            while w != 0 {
                let tz = w.trailing_zeros();
                self.hits.push(((off << 6) | tz as usize) as u32);
                w &= w - 1;
            }
        }
        &self.hits
    }

    pub fn query_count(&self) -> usize {
        self.data.len() >> self.word_shift
    }

    pub fn shift(&self) -> usize {
        self.shift
    }

    pub fn words_per_query(&self) -> usize {
        self.words_per_query
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FingerPrint {
    pub r: [Letter; 48],
}

impl Default for FingerPrint {
    fn default() -> Self {
        Self { r: [0; 48] }
    }
}

impl FingerPrint {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_array(a: &[Letter; 48]) -> Self {
        Self { r: *a }
    }

    pub fn from_seq_center(q: &[Letter], center: usize) -> Self {
        let mut r = [0; 48];
        Self::load(q, center, &mut r);
        Self { r }
    }

    pub fn load(q: &[Letter], center: usize, dst: &mut [Letter; 48]) {
        for (i, out) in dst.iter_mut().enumerate() {
            let pos = center as isize + i as isize - 16;
            *out = if pos < 0 || pos as usize >= q.len() {
                DELIMITER_LETTER
            } else {
                letter_mask(q[pos as usize]) & LETTER_MASK
            };
        }
    }

    pub fn load_masked(q: &[Letter], center: usize, dst: &mut [Letter; 48]) {
        Self::load(q, center, dst);
    }

    pub fn match_count(&self, rhs: &FingerPrint) -> u32 {
        let mut n = 0;
        for i in 0..48 {
            if self.r[i] == rhs.r[i] {
                n += 1;
            }
        }
        n
    }
}

pub fn all_vs_all(
    a: &[[Letter; 48]],
    b: &[[Letter; 48]],
    out: &mut HitField,
    hamming_filter_id: u32,
) {
    let na = a.len() as u32;
    let nb = b.len() as u32;
    let na2 = na & !3;
    let mut i = 0u32;
    while i < na2 {
        let e1 = FingerPrint::from_array(&a[i as usize]);
        let e2 = FingerPrint::from_array(&a[i as usize + 1]);
        let e3 = FingerPrint::from_array(&a[i as usize + 2]);
        let e4 = FingerPrint::from_array(&a[i as usize + 3]);
        for j in 0..nb {
            let fb = FingerPrint::from_array(&b[j as usize]);
            out.set(i, j, e1.match_count(&fb) >= hamming_filter_id);
            out.set(i + 1, j, e2.match_count(&fb) >= hamming_filter_id);
            out.set(i + 2, j, e3.match_count(&fb) >= hamming_filter_id);
            out.set(i + 3, j, e4.match_count(&fb) >= hamming_filter_id);
        }
        i += 4;
    }
    while i < na {
        let e = FingerPrint::from_array(&a[i as usize]);
        for j in 0..nb {
            out.set(
                i,
                j,
                e.match_count(&FingerPrint::from_array(&b[j as usize])) >= hamming_filter_id,
            );
        }
        i += 1;
    }
}

pub fn all_vs_all_self(a: &[[Letter; 48]], out: &mut HitField, hamming_filter_id: u32) {
    let na = a.len() as u32;
    let na2 = na & !3;
    let mut i = 0u32;
    while i < na2 {
        let e1 = FingerPrint::from_array(&a[i as usize]);
        let e2 = FingerPrint::from_array(&a[i as usize + 1]);
        let e3 = FingerPrint::from_array(&a[i as usize + 2]);
        let e4 = FingerPrint::from_array(&a[i as usize + 3]);
        out.set(i, i + 1, e1.match_count(&e2) >= hamming_filter_id);
        out.set(i, i + 2, e1.match_count(&e3) >= hamming_filter_id);
        out.set(i, i + 3, e1.match_count(&e4) >= hamming_filter_id);
        out.set(i + 1, i + 2, e2.match_count(&e3) >= hamming_filter_id);
        out.set(i + 1, i + 3, e2.match_count(&e4) >= hamming_filter_id);
        out.set(i + 2, i + 3, e3.match_count(&e4) >= hamming_filter_id);
        for j in i + 4..na {
            let fa = FingerPrint::from_array(&a[j as usize]);
            out.set(i, j, e1.match_count(&fa) >= hamming_filter_id);
            out.set(i + 1, j, e2.match_count(&fa) >= hamming_filter_id);
            out.set(i + 2, j, e3.match_count(&fa) >= hamming_filter_id);
            out.set(i + 3, j, e4.match_count(&fa) >= hamming_filter_id);
        }
        i += 4;
    }
    while i < na {
        let e = FingerPrint::from_array(&a[i as usize]);
        for j in i + 1..na {
            out.set(
                i,
                j,
                e.match_count(&FingerPrint::from_array(&a[j as usize])) >= hamming_filter_id,
            );
        }
        i += 1;
    }
}

pub fn all_vs_all_mutual_cov(
    query_lengths: &[usize],
    target_lengths: &[usize],
    a: &[[Letter; 48]],
    b: &[[Letter; 48]],
    out: &mut HitField,
    hamming_filter_id: u32,
    min_length_ratio: f64,
) {
    let na = a.len() as u32;
    let nb = b.len() as u32;
    let mut j0 = 0u32;
    let mut j1 = 0u32;
    for i in 0..na {
        let e = FingerPrint::from_array(&a[i as usize]);
        let qlen = query_lengths[i as usize] as f64;
        while j0 < nb {
            let tlen = target_lengths[j0 as usize] as f64;
            let lr = qlen / tlen;
            if lr >= min_length_ratio {
                break;
            }
            j0 += 1;
        }
        j1 = j1.max(j0);
        for j in j0..j1 {
            out.set(
                i,
                j,
                e.match_count(&FingerPrint::from_array(&b[j as usize])) >= hamming_filter_id,
            );
        }
        while j1 < nb {
            let tlen = target_lengths[j1 as usize] as f64;
            let lr = tlen / qlen;
            if lr < min_length_ratio {
                break;
            }
            out.set(
                i,
                j1,
                e.match_count(&FingerPrint::from_array(&b[j1 as usize])) >= hamming_filter_id,
            );
            j1 += 1;
        }
    }
}

pub fn all_vs_all_self_mutual_cov(
    lengths: &[usize],
    a: &[[Letter; 48]],
    out: &mut HitField,
    hamming_filter_id: u32,
    min_length_ratio: f64,
) -> u64 {
    let na = a.len() as u32;
    let mut seed_hits = 0u64;
    for i in 0..na {
        let e = FingerPrint::from_array(&a[i as usize]);
        let qlen = lengths[i as usize] as f64;
        for j in i + 1..na {
            let tlen = lengths[j as usize] as f64;
            if tlen / qlen < min_length_ratio {
                break;
            }
            seed_hits += 1;
            out.set(
                i,
                j,
                e.match_count(&FingerPrint::from_array(&a[j as usize])) >= hamming_filter_id,
            );
        }
    }
    seed_hits
}

pub fn load_fps<L: SeedLoc>(seed_locs: &[L], seqs: &SequenceSet, out: &mut Container) {
    out.clear();
    out.resize(seed_locs.len(), [0; 48]);
    for (dst, &p) in out.iter_mut().zip(seed_locs) {
        FingerPrint::load(seqs.data(), p.loc(), dst);
    }
}

/// Matches C++ `stage1(...)`.
pub fn stage1<L, F>(
    q: &[L],
    s: &[L],
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    tile_size: usize,
    hamming_filter_id: u32,
    mut search_tile: F,
) -> u64
where
    L: SeedLoc,
    F: FnMut(&mut HitField, usize, usize, &[L], &[L]),
{
    assert!(tile_size > 0);
    let mut vq = Container::new();
    let mut vs = Container::new();
    load_fps(s, target_seqs, &mut vs);
    load_fps(q, query_seqs, &mut vq);
    let mut hits = HitField::new();
    let qs = vq.len();
    let ss = vs.len();
    let mut i = 0usize;
    while i < qs {
        let mut j = 0usize;
        while j < ss {
            let tq = tile_size.min(qs - i);
            let ts = tile_size.min(ss - j);
            hits.init(tq, ts);
            all_vs_all(&vq[i..i + tq], &vs[j..j + ts], &mut hits, hamming_filter_id);
            search_tile(&mut hits, i, j, q, s);
            j += tile_size;
        }
        i += tile_size;
    }
    (q.len() * s.len()) as u64
}

/// Matches C++ `stage1_self(...)`.
pub fn stage1_self<L, F>(
    s: &[L],
    target_seqs: &SequenceSet,
    tile_size: usize,
    hamming_filter_id: u32,
    mut search_tile: F,
) -> u64
where
    L: SeedLoc,
    F: FnMut(&mut HitField, usize, usize, &[L], &[L]),
{
    assert!(tile_size > 0);
    let mut vs = Container::new();
    load_fps(s, target_seqs, &mut vs);
    let mut hits = HitField::new();
    let ss = vs.len();
    let mut i = 0usize;
    while i < ss {
        let tq = tile_size.min(ss - i);
        hits.init(tq, tq);
        all_vs_all_self(&vs[i..i + tq], &mut hits, hamming_filter_id);
        search_tile(&mut hits, i, i, s, s);
        let mut j = i + tile_size;
        while j < ss {
            let ts = tile_size.min(ss - j);
            hits.init(tq, ts);
            all_vs_all(&vs[i..i + tq], &vs[j..j + ts], &mut hits, hamming_filter_id);
            search_tile(&mut hits, i, j, s, s);
            j += tile_size;
        }
        i += tile_size;
    }
    (s.len() * (s.len().saturating_sub(1)) / 2) as u64
}

/// Matches C++ `stage1_query_lin(...)`.
pub fn stage1_query_lin<F>(
    q: &[PackedLoc],
    s: &[PackedLoc],
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    tile_size: usize,
    hamming_filter_id: u32,
    mut search_tile: F,
) -> u64
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLoc], &[PackedLoc]),
{
    assert!(tile_size > 0);
    let mut vq = Container::new();
    let mut vs = Container::new();
    load_fps(&q[..1], query_seqs, &mut vq);
    load_fps(s, target_seqs, &mut vs);
    let mut hits = HitField::new();
    let ss = vs.len();
    let mut j = 0usize;
    while j < ss {
        let ts = tile_size.min(ss - j);
        hits.init(1, ts);
        all_vs_all(&vq, &vs[j..j + ts], &mut hits, hamming_filter_id);
        search_tile(&mut hits, 0, j, q, s);
        j += tile_size;
    }
    s.len() as u64
}

/// Matches C++ `stage1_query_lin_ranked(...)`.
pub fn stage1_query_lin_ranked<F>(
    q: &[PackedLocId],
    s: &[PackedLocId],
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    ranking: &KmerRanking,
    tile_size: usize,
    hamming_filter_id: u32,
    mut search_tile: F,
) -> u64
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLocId], &[PackedLocId]),
{
    assert!(tile_size > 0);
    let pivot = ranking.highest_ranking(q) as usize;
    let mut vq = Container::new();
    let mut vs = Container::new();
    load_fps(&q[pivot..pivot + 1], query_seqs, &mut vq);
    load_fps(s, target_seqs, &mut vs);
    let mut hits = HitField::new();
    let ss = vs.len();
    let mut j = 0usize;
    while j < ss {
        let ts = tile_size.min(ss - j);
        hits.init(1, ts);
        all_vs_all(&vq, &vs[j..j + ts], &mut hits, hamming_filter_id);
        search_tile(&mut hits, pivot, j, q, s);
        j += tile_size;
    }
    s.len() as u64
}

/// Matches C++ `stage1_target_lin(...)`.
pub fn stage1_target_lin<L, F>(
    q: &[L],
    s: &[L],
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    tile_size: usize,
    hamming_filter_id: u32,
    mut search_tile: F,
) -> u64
where
    L: SeedLoc,
    F: FnMut(&mut HitField, usize, usize, &[L], &[L]),
{
    assert!(tile_size > 0);
    let mut vq = Container::new();
    let mut vs = Container::new();
    load_fps(q, query_seqs, &mut vq);
    load_fps(&s[..1], target_seqs, &mut vs);
    let mut hits = HitField::new();
    let qs = vq.len();
    let mut j = 0usize;
    while j < qs {
        let tq = tile_size.min(qs - j);
        hits.init(tq, 1);
        all_vs_all(&vq[j..j + tq], &vs, &mut hits, hamming_filter_id);
        search_tile(&mut hits, j, 0, q, s);
        j += tile_size;
    }
    q.len() as u64
}

/// Matches C++ `stage1_longest_combo_lin(...)`.
pub fn stage1_longest_combo_lin<F>(
    q: &[PackedLocId],
    s: &[PackedLocId],
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    tile_size: usize,
    hamming_filter_id: u32,
    mut search_tile: F,
) -> u64
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLocId], &[PackedLocId]),
{
    assert!(tile_size > 0);
    let mut query_pivot = true;
    let mut pivot = 0usize;
    let mut pivot_len = query_seqs.length(q[0].block_id as usize);
    for (i, loc) in q.iter().enumerate().skip(1) {
        let len = query_seqs.length(loc.block_id as usize);
        if len > pivot_len {
            pivot = i;
            pivot_len = len;
        }
    }
    for (i, loc) in s.iter().enumerate() {
        let len = target_seqs.length(loc.block_id as usize);
        if len > pivot_len {
            query_pivot = false;
            pivot = i;
            pivot_len = len;
        }
    }
    if query_pivot {
        let mut vq = Container::new();
        let mut vs = Container::new();
        load_fps(&q[pivot..pivot + 1], query_seqs, &mut vq);
        load_fps(s, target_seqs, &mut vs);
        let mut hits = HitField::new();
        let ss = vs.len();
        let mut j = 0usize;
        while j < ss {
            let ts = tile_size.min(ss - j);
            hits.init(1, ts);
            all_vs_all(&vq, &vs[j..j + ts], &mut hits, hamming_filter_id);
            search_tile(&mut hits, pivot, j, q, s);
            j += tile_size;
        }
        s.len() as u64
    } else {
        let mut vq = Container::new();
        let mut vs = Container::new();
        load_fps(q, query_seqs, &mut vq);
        load_fps(&s[pivot..pivot + 1], target_seqs, &mut vs);
        let mut hits = HitField::new();
        let qs = vq.len();
        let mut j = 0usize;
        while j < qs {
            let tq = tile_size.min(qs - j);
            hits.init(tq, 1);
            all_vs_all(&vq[j..j + tq], &vs, &mut hits, hamming_filter_id);
            search_tile(&mut hits, j, pivot, q, s);
            j += tile_size;
        }
        q.len() as u64
    }
}

/// Matches C++ `stage1_mutual_cov(...)`.
pub fn stage1_mutual_cov<F>(
    q: &[PackedLocId],
    s: &[PackedLocId],
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    tile_size: usize,
    hamming_filter_id: u32,
    min_length_ratio: f64,
    mut search_tile: F,
) -> u64
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLocId], &[PackedLocId]),
{
    assert!(tile_size > 0);
    let mut vq = Container::new();
    let mut vs = Container::new();
    load_fps(s, target_seqs, &mut vs);
    load_fps(q, query_seqs, &mut vq);
    let mut hits = HitField::new();
    let qs = vq.len();
    let ss = vs.len();
    let mut i = 0usize;
    while i < qs {
        let mut j = 0usize;
        while j < ss {
            let tq = tile_size.min(qs - i);
            let ts = tile_size.min(ss - j);
            hits.init(tq, ts);
            let q_len: Vec<usize> = q[i..i + tq]
                .iter()
                .map(|x| query_seqs.length(x.block_id as usize) as usize)
                .collect();
            let s_len: Vec<usize> = s[j..j + ts]
                .iter()
                .map(|x| target_seqs.length(x.block_id as usize) as usize)
                .collect();
            all_vs_all_mutual_cov(
                &q_len,
                &s_len,
                &vq[i..i + tq],
                &vs[j..j + ts],
                &mut hits,
                hamming_filter_id,
                min_length_ratio,
            );
            search_tile(&mut hits, i, j, q, s);
            j += tile_size;
        }
        i += tile_size;
    }
    (q.len() * s.len()) as u64
}

/// Matches C++ `stage1_self_mutual_cov(...)`.
pub fn stage1_self_mutual_cov<F>(
    s: &[PackedLocId],
    target_seqs: &SequenceSet,
    hamming_filter_id: u32,
    min_length_ratio: f64,
    mut search_tile: F,
) -> u64
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLocId], &[PackedLocId]),
{
    let mut vs = Container::new();
    load_fps(s, target_seqs, &mut vs);
    let mut hits = HitField::new();
    hits.init(vs.len(), vs.len());
    let lengths: Vec<usize> = s
        .iter()
        .map(|x| target_seqs.length(x.block_id as usize) as usize)
        .collect();
    let seed_hits = all_vs_all_self_mutual_cov(
        &lengths,
        &vs,
        &mut hits,
        hamming_filter_id,
        min_length_ratio,
    );
    search_tile(&mut hits, 0, 0, s, s);
    seed_hits
}

/// Matches C++ `stage1_mutual_cov_query_lin(...)`.
pub fn stage1_mutual_cov_query_lin<F>(
    q: &[PackedLocId],
    s: &[PackedLocId],
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    hamming_filter_id: u32,
    min_length_ratio: f64,
    self_search: bool,
    mut search_tile: F,
) -> u64
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLocId], &[PackedLocId]),
{
    let mut vq = Container::new();
    let mut vs = Container::new();
    load_fps(s, target_seqs, &mut vs);
    load_fps(q, query_seqs, &mut vq);
    let mut hits = HitField::new();
    let qs = vq.len();
    let ss = vs.len();
    let mut j = 0usize;
    let mut seed_hits = 0u64;
    let mut i = 0usize;
    while i < qs {
        let qlen = query_seqs.length(q[i].block_id as usize) as f64;
        let mut j1 = j;
        while j1 < ss {
            let tlen = target_seqs.length(s[j1].block_id as usize) as f64;
            if tlen / qlen < min_length_ratio {
                break;
            }
            j1 += 1;
        }
        let ts = j1 - j;
        hits.init(1, ts);
        let qpos = if self_search { i + (j1 - j) / 2 } else { i };
        all_vs_all(
            &vq[qpos..qpos + 1],
            &vs[j..j + ts],
            &mut hits,
            hamming_filter_id,
        );
        search_tile(&mut hits, qpos, j, q, s);
        seed_hits += ts as u64;
        j = j1;
        if j == ss {
            break;
        }
        let tlen = target_seqs.length(s[j].block_id as usize) as f64;
        while i < qs {
            let qlen = query_seqs.length(q[i].block_id as usize) as f64;
            if tlen / qlen >= min_length_ratio {
                break;
            }
            i += 1;
        }
    }
    seed_hits
}

/// Matches C++ `stage1_mutual_cov_target_lin(...)`.
pub fn stage1_mutual_cov_target_lin<F>(
    q: &[PackedLocId],
    s: &[PackedLocId],
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    hamming_filter_id: u32,
    min_length_ratio: f64,
    mut search_tile: F,
) -> u64
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLocId], &[PackedLocId]),
{
    let mut vq = Container::new();
    let mut vs = Container::new();
    load_fps(s, target_seqs, &mut vs);
    load_fps(q, query_seqs, &mut vq);
    let mut hits = HitField::new();
    let qs = vq.len();
    let ss = vs.len();
    let mut i = 0usize;
    let mut seed_hits = 0u64;
    let mut j = 0usize;
    while j < ss {
        let tlen = target_seqs.length(s[j].block_id as usize) as f64;
        let mut i1 = i;
        while i1 < qs {
            let qlen = query_seqs.length(q[i1].block_id as usize) as f64;
            if qlen / tlen < min_length_ratio {
                break;
            }
            i1 += 1;
        }
        let tq = i1 - i;
        hits.init(tq, 1);
        all_vs_all(&vq[i..i + tq], &vs[j..j + 1], &mut hits, hamming_filter_id);
        search_tile(&mut hits, i, j, q, s);
        seed_hits += tq as u64;
        i = i1;
        if i == qs {
            break;
        }
        let qlen = query_seqs.length(q[i].block_id as usize) as f64;
        while j < ss {
            let tlen = target_seqs.length(s[j].block_id as usize) as f64;
            if qlen / tlen >= min_length_ratio {
                break;
            }
            j += 1;
        }
    }
    seed_hits
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seqset(lengths: &[usize]) -> SequenceSet {
        let mut seqs = SequenceSet::new();
        for &len in lengths {
            let seq: Vec<Letter> = (0..len).map(|i| (i % 20) as Letter).collect();
            seqs.push(&seq);
        }
        seqs
    }

    fn fp(v: Letter) -> [Letter; 48] {
        [v; 48]
    }

    #[test]
    fn test_hit_field_set_hits_and_query_count() {
        let mut field = HitField::new();
        field.init(3, 100);
        assert_eq!(field.shift(), 8);
        assert_eq!(field.words_per_query(), 4);
        assert_eq!(field.query_count(), 3);
        field.set(1, 2, true);
        field.set(1, 99, true);
        field.set(1, 2, false);
        field.set(2, 64, true);
        assert_eq!(field.hits(0), &[]);
        assert_eq!(field.hits(1), &[99]);
        assert_eq!(field.hits(2), &[64]);
    }

    #[test]
    fn test_fingerprint_load_and_match() {
        let seq: Vec<Letter> = (0..80).map(|i| (i % 20) as Letter).collect();
        let mut a = [0; 48];
        FingerPrint::load(&seq, 20, &mut a);
        assert_eq!(a[0], seq[4]);
        assert_eq!(a[47], seq[51]);
        let fp1 = FingerPrint::from_array(&a);
        let fp2 = FingerPrint::from_seq_center(&seq, 20);
        assert_eq!(fp1.match_count(&fp2), 48);
        a[0] = -1;
        assert_eq!(FingerPrint::from_array(&a).match_count(&fp2), 47);
    }

    #[test]
    fn test_fingerprint_load_masks_and_pads_boundaries() {
        let mut seq: Vec<Letter> = (0..20).map(|i| (i % 20) as Letter).collect();
        seq[0] = -128 | 5;
        let mut a = [0; 48];
        FingerPrint::load(&seq, 0, &mut a);
        assert!(a[..16].iter().all(|&x| x == DELIMITER_LETTER));
        assert_eq!(a[16], 5);
        assert!(a[36..].iter().all(|&x| x == DELIMITER_LETTER));
    }

    #[test]
    fn test_all_vs_all_sets_matching_targets() {
        let a = [fp(1), fp(2), fp(3), fp(4), fp(5)];
        let mut b = [fp(1), fp(9), fp(5)];
        b[1][0] = 2;
        let mut field = HitField::new();
        field.init(a.len(), b.len());
        all_vs_all(&a, &b, &mut field, 48);
        assert_eq!(field.hits(0), &[0]);
        assert_eq!(field.hits(1), &[]);
        assert_eq!(field.hits(4), &[2]);
        all_vs_all(&a, &b, &mut field, 1);
        assert_eq!(field.hits(1), &[1]);
    }

    #[test]
    fn test_all_vs_all_self_sets_upper_triangle() {
        let a = [fp(1), fp(2), fp(1), fp(3), fp(1)];
        let mut field = HitField::new();
        field.init(a.len(), a.len());
        all_vs_all_self(&a, &mut field, 48);
        assert_eq!(field.hits(0), &[2, 4]);
        assert_eq!(field.hits(1), &[]);
        assert_eq!(field.hits(2), &[4]);
    }

    #[test]
    fn test_all_vs_all_mutual_cov_windows_by_length_ratio() {
        let q = [fp(1), fp(2), fp(3)];
        let s = [fp(1), fp(2), fp(3), fp(4)];
        let q_len = [100, 80, 50];
        let s_len = [120, 100, 80, 40];
        let mut field = HitField::new();
        field.init(q.len(), s.len());
        all_vs_all_mutual_cov(&q_len, &s_len, &q, &s, &mut field, 48, 0.75);
        assert_eq!(field.hits(0), &[0]);
        assert_eq!(field.hits(1), &[1]);
        assert_eq!(field.hits(2), &[]);
    }

    #[test]
    fn test_all_vs_all_self_mutual_cov_counts_considered_pairs() {
        let a = [fp(1), fp(2), fp(1), fp(4)];
        let lengths = [100, 90, 75, 40];
        let mut field = HitField::new();
        field.init(a.len(), a.len());
        let seed_hits = all_vs_all_self_mutual_cov(&lengths, &a, &mut field, 48, 0.8);
        assert_eq!(seed_hits, 2);
        assert_eq!(field.hits(0), &[]);
        assert_eq!(field.hits(1), &[]);
    }

    #[test]
    fn test_load_fps() {
        let mut seqs = SequenceSet::new();
        let seq: Vec<Letter> = (0..100).map(|i| (i % 23) as Letter).collect();
        seqs.push(&seq);
        let p0 = seqs.position(0, 20);
        let p1 = seqs.position(0, 30);
        let mut out = Vec::new();
        load_fps(&[p0, p1], &seqs, &mut out);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0][0], letter_mask(seqs.data()[p0 - 16]) & LETTER_MASK);
        assert_eq!(out[1][47], seqs.data()[p1 + 31]);
    }

    #[test]
    fn test_stage1_tiles_all_pairs() {
        let query = seqset(&[80]);
        let target = seqset(&[80]);
        let q = [
            query.position(0, 20),
            query.position(0, 30),
            query.position(0, 40),
        ];
        let s = [target.position(0, 20), target.position(0, 30)];
        let mut calls = Vec::new();
        let seed_hits = stage1(&q, &s, &query, &target, 2, 48, |hits, i, j, _, _| {
            calls.push((i, j, hits.query_count()));
        });
        assert_eq!(seed_hits, 6);
        assert_eq!(calls, vec![(0, 0, 2), (2, 0, 1)]);
    }

    #[test]
    fn test_stage1_self_tiles_diagonal_and_off_diagonal() {
        let target = seqset(&[90]);
        let s = [
            target.position(0, 20),
            target.position(0, 30),
            target.position(0, 40),
        ];
        let mut calls = Vec::new();
        let seed_hits = stage1_self(&s, &target, 2, 48, |hits, i, j, _, _| {
            calls.push((i, j, hits.query_count()));
        });
        assert_eq!(seed_hits, 3);
        assert_eq!(calls, vec![(0, 0, 2), (0, 2, 2), (2, 2, 1)]);
    }

    #[test]
    fn test_stage1_query_and_target_linear_offsets() {
        let query = seqset(&[90]);
        let target = seqset(&[90]);
        let q = [
            PackedLoc::new(query.position(0, 20) as u64),
            PackedLoc::new(query.position(0, 30) as u64),
        ];
        let s = [
            PackedLoc::new(target.position(0, 20) as u64),
            PackedLoc::new(target.position(0, 30) as u64),
            PackedLoc::new(target.position(0, 40) as u64),
        ];
        let mut query_calls = Vec::new();
        assert_eq!(
            stage1_query_lin(&q, &s, &query, &target, 2, 48, |_, i, j, _, _| {
                query_calls.push((i, j));
            }),
            3
        );
        assert_eq!(query_calls, vec![(0, 0), (0, 2)]);

        let mut target_calls = Vec::new();
        assert_eq!(
            stage1_target_lin(&q, &s, &query, &target, 2, 48, |_, i, j, _, _| {
                target_calls.push((i, j));
            }),
            2
        );
        assert_eq!(target_calls, vec![(0, 0)]);
    }

    #[test]
    fn test_stage1_query_lin_ranked_uses_highest_ranking() {
        let query = seqset(&[90, 110, 100]);
        let target = seqset(&[90]);
        let q = [
            PackedLocId::new(query.position(0, 20) as u64, 0),
            PackedLocId::new(query.position(1, 20) as u64, 1),
            PackedLocId::new(query.position(2, 20) as u64, 2),
        ];
        let s = [PackedLocId::new(target.position(0, 20) as u64, 0)];
        let ranking = KmerRanking::from_queries(&query);
        let mut calls = Vec::new();
        let seed_hits =
            stage1_query_lin_ranked(&q, &s, &query, &target, &ranking, 4, 48, |_, i, j, _, _| {
                calls.push((i, j));
            });
        assert_eq!(seed_hits, 1);
        assert_eq!(calls, vec![(1, 0)]);
    }

    #[test]
    fn test_stage1_longest_combo_selects_target_pivot() {
        let query = seqset(&[80, 90]);
        let target = seqset(&[85, 120]);
        let q = [
            PackedLocId::new(query.position(0, 20) as u64, 0),
            PackedLocId::new(query.position(1, 20) as u64, 1),
        ];
        let s = [
            PackedLocId::new(target.position(0, 20) as u64, 0),
            PackedLocId::new(target.position(1, 20) as u64, 1),
        ];
        let mut calls = Vec::new();
        let seed_hits =
            stage1_longest_combo_lin(&q, &s, &query, &target, 4, 48, |_, i, j, _, _| {
                calls.push((i, j));
            });
        assert_eq!(seed_hits, 2);
        assert_eq!(calls, vec![(0, 1)]);
    }

    #[test]
    fn test_stage1_self_mutual_cov_counts_filtered_pairs() {
        let target = seqset(&[100, 90, 70]);
        let s = [
            PackedLocId::new(target.position(0, 20) as u64, 0),
            PackedLocId::new(target.position(1, 20) as u64, 1),
            PackedLocId::new(target.position(2, 20) as u64, 2),
        ];
        let mut calls = Vec::new();
        let seed_hits = stage1_self_mutual_cov(&s, &target, 48, 0.8, |hits, i, j, _, _| {
            calls.push((i, j, hits.query_count()));
        });
        assert_eq!(seed_hits, 1);
        assert_eq!(calls, vec![(0, 0, 3)]);
    }

    #[test]
    fn test_stage1_mutual_cov_linear_windows() {
        let query = seqset(&[100, 80, 60]);
        let target = seqset(&[120, 100, 80, 60]);
        let q: Vec<_> = (0..3)
            .map(|i| PackedLocId::new(query.position(i, 20) as u64, i as u32))
            .collect();
        let s: Vec<_> = (0..4)
            .map(|i| PackedLocId::new(target.position(i, 20) as u64, i as u32))
            .collect();
        let mut query_calls = Vec::new();
        let query_hits = stage1_mutual_cov_query_lin(
            &q,
            &s,
            &query,
            &target,
            48,
            0.75,
            false,
            |hits, i, j, _, _| query_calls.push((i, j, hits.words_per_query())),
        );
        assert_eq!(query_hits, 4);
        assert_eq!(query_calls, vec![(0, 0, 4), (1, 3, 4)]);

        let mut target_calls = Vec::new();
        let target_hits =
            stage1_mutual_cov_target_lin(&q, &s, &query, &target, 48, 0.75, |hits, i, j, _, _| {
                target_calls.push((i, j, hits.query_count()))
            });
        assert_eq!(target_hits, 3);
        assert_eq!(target_calls, vec![(0, 0, 1), (1, 1, 1), (2, 2, 1)]);
    }

    #[test]
    fn test_stage1_dispatch_packed_loc_id_branches() {
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig {
                lin_stage1_combo: true,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1LongestComboLin
        );
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig {
                lin_stage1_query: true,
                min_length_ratio: 0.5,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1MutualCovQueryLin
        );
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig {
                lin_stage1_query: true,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1QueryLinRanked
        );
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig {
                lin_stage1_target: true,
                min_length_ratio: 0.5,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1MutualCovTargetLin
        );
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig {
                lin_stage1_target: true,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1TargetLin
        );
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig {
                min_length_ratio: 0.5,
                self_search: true,
                current_ref_block: 0,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1SelfMutualCov
        );
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig {
                min_length_ratio: 0.5,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1MutualCov
        );
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig {
                self_search: true,
                current_ref_block: 0,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1Self
        );
        assert_eq!(
            stage1_dispatch_packed_loc_id(&Stage1DispatchConfig::default()),
            Stage1KernelKind::Stage1
        );
    }

    #[test]
    fn test_stage1_dispatch_packed_loc_and_keep_target_id() {
        assert_eq!(
            stage1_dispatch_packed_loc(&Stage1DispatchConfig {
                lin_stage1_query: true,
                min_length_ratio: 0.5,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1QueryLin
        );
        assert_eq!(
            stage1_dispatch_packed_loc(&Stage1DispatchConfig {
                lin_stage1_target: true,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1TargetLin
        );
        assert_eq!(
            stage1_dispatch_packed_loc(&Stage1DispatchConfig {
                self_search: true,
                current_ref_block: 0,
                ..Stage1DispatchConfig::default()
            }),
            Stage1KernelKind::Stage1Self
        );
        assert_eq!(
            stage1_dispatch_packed_loc(&Stage1DispatchConfig::default()),
            Stage1KernelKind::Stage1
        );

        assert!(!keep_target_id(&Stage1DispatchConfig::default()));
        assert!(keep_target_id(&Stage1DispatchConfig {
            hit_keep_target_id: true,
            ..Stage1DispatchConfig::default()
        }));
        assert!(keep_target_id(&Stage1DispatchConfig {
            min_length_ratio: 0.1,
            ..Stage1DispatchConfig::default()
        }));
        assert!(keep_target_id(&Stage1DispatchConfig {
            global_ranking_targets: true,
            ..Stage1DispatchConfig::default()
        }));
        assert!(keep_target_id(&Stage1DispatchConfig {
            self_search: true,
            current_ref_block: 0,
            ..Stage1DispatchConfig::default()
        }));
        assert!(keep_target_id(&Stage1DispatchConfig {
            lin_stage1_combo: true,
            ..Stage1DispatchConfig::default()
        }));
    }
}
