use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::Letter;
use crate::dp::swipe::DpTarget;
use crate::dp::ungapped::{score_range, DiagonalSegment};
use crate::stats;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::hsp::{Anchor, ApproxHsp};
use crate::util::interval::Interval;
use std::collections::BTreeMap;
use std::io::{self, Write};

pub const SPACE_PENALTY: f64 = 0.1;
pub const LINK_PADDING: i32 = 10;
pub const REVERSE_LINK_MIN_OVERHANG: i32 = 10;
pub const DEFAULT_CHAINING_RANGE_COVER: usize = 8;
pub const DEFAULT_CHAINING_STACKED_HSP_RATIO: f64 = 0.5;
pub const END: usize = usize::MAX;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DiagonalNode {
    pub segment: DiagonalSegment,
    pub link_idx: i32,
    pub prefix_score: i32,
    pub path_max: i32,
    pub path_min: i32,
}

impl DiagonalNode {
    pub const ESTIMATE: i32 = 0;
    pub const FINISHED: i32 = 1;

    pub fn new(query_pos: i32, subject_pos: i32, len: i32, score: i32, link_idx: i32) -> Self {
        Self::from_segment_with_link(
            DiagonalSegment::new(query_pos, subject_pos, len, score),
            link_idx,
        )
    }

    pub fn from_segment(segment: DiagonalSegment) -> Self {
        Self::from_segment_with_link(segment, -1)
    }

    pub fn from_segment_with_link(segment: DiagonalSegment, link_idx: i32) -> Self {
        let score = segment.score;
        Self {
            segment,
            link_idx,
            prefix_score: score,
            path_max: score,
            path_min: score,
        }
    }

    pub fn deactivate(&mut self) {
        self.link_idx = 0;
    }

    pub fn reset(&mut self) {
        self.link_idx = -1;
        self.prefix_score = self.segment.score;
        self.path_max = self.segment.score;
        self.path_min = self.segment.score;
    }

    pub fn is_maximum(&self) -> bool {
        self.path_max == self.prefix_score
    }

    pub fn rel_score(&self) -> i32 {
        if self.prefix_score == self.path_max {
            self.prefix_score
        } else {
            self.prefix_score - self.path_min
        }
    }

    pub fn cmp_prefix_score(x: &DiagonalNode, y: &DiagonalNode) -> bool {
        x.prefix_score > y.prefix_score
    }

    pub fn cmp_rel_score(x: &DiagonalNode, y: &DiagonalNode) -> bool {
        x.rel_score() > y.rel_score()
    }
}

impl std::ops::Deref for DiagonalNode {
    type Target = DiagonalSegment;

    fn deref(&self) -> &Self::Target {
        &self.segment
    }
}

impl std::ops::DerefMut for DiagonalNode {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.segment
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Edge {
    pub prefix_score: i32,
    pub path_max: i32,
    pub j: i32,
    pub path_min: i32,
    pub prefix_score_begin: i32,
    pub node_in: u32,
    pub node_out: u32,
}

impl Edge {
    pub fn new(
        prefix_score: i32,
        path_max: i32,
        j: i32,
        node_in: u32,
        node_out: u32,
        path_min: i32,
        prefix_score_begin: i32,
    ) -> Self {
        Self {
            prefix_score,
            path_max,
            j,
            path_min,
            prefix_score_begin,
            node_in,
            node_out,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DiagGraph {
    pub nodes: Vec<DiagonalNode>,
    pub edges: Vec<Edge>,
}

impl DiagGraph {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn init(&mut self) {
        self.nodes.clear();
        self.edges.clear();
    }

    pub fn init_node(&mut self, node: usize) -> Result<(), String> {
        if self.edges.len() >= i32::MAX as usize {
            return Err("Too many edges.".to_string());
        }
        self.nodes[node].link_idx = self.edges.len() as i32;
        Ok(())
    }

    pub fn load(&mut self, segments: &[DiagonalSegment]) {
        let mut d = i32::MIN;
        let mut max_j_end = d;
        for segment in segments {
            let d2 = segment.diag();
            if d2 != d {
                d = d2;
                self.nodes.push(DiagonalNode::from_segment(segment.clone()));
                max_j_end = self.nodes.last().unwrap().subject_end();
            } else if max_j_end < segment.j {
                self.nodes.push(DiagonalNode::from_segment(segment.clone()));
                max_j_end = max_j_end.max(self.nodes.last().unwrap().subject_end());
            }
        }
    }

    pub fn clear_edges(&mut self) {
        self.edges.clear();
        for node in &mut self.nodes {
            node.deactivate();
        }
    }

    pub fn sort(&mut self) {
        self.nodes
            .sort_by(|x, y| x.j.cmp(&y.j).then_with(|| x.i.cmp(&y.i)));
    }

    pub fn prune(&mut self) {
        let mut finished = Vec::new();
        let mut window: Vec<DiagonalNode> = Vec::new();
        for d in &self.nodes {
            let mut n = 0usize;
            let mut i = 0usize;
            while i < window.len() {
                if window[i].subject_end() > d.j {
                    if window[i].score >= d.score
                        && window[i].j <= d.j
                        && window[i].subject_end() >= d.subject_end()
                    {
                        n += 1;
                    }
                    i += 1;
                } else {
                    finished.push(window.remove(i));
                }
            }
            if n <= DEFAULT_CHAINING_RANGE_COVER {
                window.push(d.clone());
            }
        }
        finished.extend(window);
        self.nodes = finished;
    }

    pub fn add_edge(&mut self, edge: Edge) -> usize {
        for j in edge.node_in as usize + 1..self.nodes.len() {
            if self.nodes[j].link_idx == -1 {
                break;
            }
            self.nodes[j].link_idx += 1;
        }
        let d = &mut self.nodes[edge.node_in as usize];
        if edge.prefix_score > d.prefix_score {
            d.prefix_score = edge.prefix_score;
            d.path_max = edge.path_max;
            d.path_min = edge.path_min;
        }
        let insert = d.link_idx as usize;
        d.link_idx += 1;
        self.edges.insert(insert, edge);
        insert
    }

    pub fn get_edge(&self, node: usize, j: i32) -> Option<usize> {
        let d = &self.nodes[node];
        if d.score == 0 {
            return Some((d.link_idx - 1) as usize);
        }
        if self.edges.is_empty() || d.link_idx <= 0 {
            return None;
        }
        let mut max_score = d.score;
        let mut max_i = None;
        let mut i = d.link_idx - 1;
        while i >= 0 {
            let edge = &self.edges[i as usize];
            if edge.node_in as usize != node {
                break;
            }
            if edge.j < j && edge.prefix_score > max_score {
                max_i = Some(i as usize);
                max_score = edge.prefix_score;
            }
            i -= 1;
        }
        max_i
    }

    pub fn prefix_score(&self, node: usize, j: i32) -> (i32, i32, i32) {
        if let Some(i) = self.get_edge(node, j) {
            let edge = &self.edges[i];
            (
                self.nodes[node].score.max(edge.prefix_score),
                self.nodes[node].score.max(edge.path_max),
                edge.path_min,
            )
        } else {
            let score = self.nodes[node].score;
            (score, score, score)
        }
    }

    pub fn top_node(&self) -> usize {
        let mut top_score = 0;
        let mut top_node = END;
        for (k, node) in self.nodes.iter().enumerate() {
            if node.prefix_score > top_score {
                top_node = k;
                top_score = node.prefix_score;
            }
        }
        top_node
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Link {
    pub subject_pos1: i32,
    pub query_pos1: i32,
    pub subject_pos2: i32,
    pub query_pos2: i32,
    pub score1: i32,
    pub score2: i32,
}

impl Default for Link {
    fn default() -> Self {
        Self {
            subject_pos1: -1,
            query_pos1: 0,
            subject_pos2: 0,
            query_pos2: 0,
            score1: 0,
            score2: 0,
        }
    }
}

impl Link {
    pub fn new(_target: u32, query_pos: i32, subject_pos: i32, score1: i32, score2: i32) -> Self {
        Self {
            subject_pos1: subject_pos,
            query_pos1: query_pos,
            subject_pos2: 0,
            query_pos2: 0,
            score1,
            score2,
        }
    }

    pub fn transpose(&mut self) -> &mut Self {
        std::mem::swap(&mut self.subject_pos1, &mut self.query_pos1);
        std::mem::swap(&mut self.subject_pos2, &mut self.query_pos2);
        self
    }

    pub fn reset(&mut self) {
        self.subject_pos1 = -1;
        self.score1 = 0;
        self.score2 = 0;
    }
}

pub fn get_hgap_link(
    d1: &DiagonalSegment,
    d2: &DiagonalSegment,
    query: &[Letter],
    subject: &[Letter],
    l: &mut Link,
    padding: i32,
    score_matrix: &ScoreMatrix,
) -> i32 {
    let d = d1.diag() - d2.diag();
    let j2_end =
        d2.j.max(d1.subject_last() + d + 1 + padding)
            .min(d2.subject_last());
    let mut j1;
    let space;
    if d1.subject_last() < d2.j - d - 1 {
        j1 = d1.subject_last();
        space = true;
    } else {
        j1 = (d2.j - d - 1 - padding).max(d1.j);
        space = false;
    }
    let mut j2 = j1 + d + 1;
    let mut i1 = d1.i + (j1 - d1.j);
    let mut i2 = i1 + 1;
    if j2 > d2.subject_last() {
        l.reset();
        return i32::MIN;
    }
    let mut score1 = 0;
    let mut score2 = score_range(
        query,
        subject,
        i2 as usize,
        j2 as usize,
        d2.j as usize,
        score_matrix,
    ) + d2.score
        - score_range(
            query,
            subject,
            d2.i as usize,
            d2.j as usize,
            j2 as usize,
            score_matrix,
        );
    let mut max_score = i32::MIN;
    loop {
        if score1 + score2 > max_score {
            max_score = score1 + score2;
            l.query_pos1 = i1;
            l.subject_pos1 = j1;
            l.query_pos2 = i2;
            l.subject_pos2 = j2;
            l.score1 = score1;
            l.score2 = score2;
        }
        score2 -= score_matrix.score(query[i2 as usize], subject[j2 as usize]);
        i1 += 1;
        i2 += 1;
        j1 += 1;
        j2 += 1;
        if j2 > j2_end {
            break;
        }
        score1 += score_matrix.score(query[i1 as usize], subject[j1 as usize]);
    }
    let j1_end = j2_end - d;
    if space {
        l.score1 += d1.score;
    } else {
        l.score1 += d1.score
            - score_range(
                query,
                subject,
                (d1.diag() + j1_end) as usize,
                j1_end as usize,
                d1.subject_end() as usize,
                score_matrix,
            )
            + score_range(
                query,
                subject,
                d1.query_end() as usize,
                d1.subject_end() as usize,
                j1_end as usize,
                score_matrix,
            )
            - score1;
    }
    max_score
}

pub fn get_vgap_link(
    d1: &DiagonalSegment,
    d2: &DiagonalSegment,
    query: &[Letter],
    subject: &[Letter],
    l: &mut Link,
    padding: i32,
    score_matrix: &ScoreMatrix,
) -> i32 {
    let s = get_hgap_link(
        &d1.transpose(),
        &d2.transpose(),
        subject,
        query,
        l,
        padding,
        score_matrix,
    );
    l.transpose();
    s
}

pub fn get_link(
    d1: &DiagonalSegment,
    d2: &DiagonalSegment,
    query: &[Letter],
    subject: &[Letter],
    l: &mut Link,
    padding: i32,
    score_matrix: &ScoreMatrix,
) -> i32 {
    if d1.diag() < d2.diag() {
        get_vgap_link(d1, d2, query, subject, l, padding, score_matrix)
    } else {
        get_hgap_link(d1, d2, query, subject, l, padding, score_matrix)
    }
}

pub fn merge_score(h1: &ApproxHsp, h2: &ApproxHsp) -> i32 {
    const GAP_PENALTY: f64 = 0.5;
    let gq = h2.query_range.begin - h1.query_range.end;
    let gt = h2.subject_range.begin - h1.subject_range.end;
    if gq < 0 || gt < 0 {
        return 0;
    }
    let s = h1.score + h2.score;
    if gq > gt {
        (s as f64 - gq as f64 * GAP_PENALTY - gt as f64 * SPACE_PENALTY) as i32
    } else {
        (s as f64 - gt as f64 * GAP_PENALTY - gq as f64 * SPACE_PENALTY) as i32
    }
}

pub fn merge(h1: &ApproxHsp, h2: &ApproxHsp) -> ApproxHsp {
    let mut h = ApproxHsp::new(h1.frame as u32, 0);
    h.d_max = h1.d_max.max(h2.d_max);
    h.d_min = h1.d_min.min(h2.d_min);
    h.query_range = Interval::new(h1.query_range.begin, h2.query_range.end);
    h.query_source_range = h.query_range;
    h.subject_range = Interval::new(h1.subject_range.begin, h2.subject_range.end);
    h.score = merge_score(h1, h2);
    h.evalue = 0.0;
    if h1.max_diag.segment.score > h2.max_diag.segment.score {
        h.max_diag = h1.max_diag.clone();
        h.max_diag.d_max_right = h.max_diag.d_max_right.max(h2.d_max);
        h.max_diag.d_min_right = h.max_diag.d_min_right.min(h2.d_min);
    } else {
        h.max_diag = h2.max_diag.clone();
        h.max_diag.d_max_left = h.max_diag.d_max_left.max(h1.d_max);
        h.max_diag.d_min_left = h.max_diag.d_min_left.min(h1.d_min);
    }
    h
}

pub fn merge_hsps(hsps: &mut Vec<ApproxHsp>) {
    let mut it = 0usize;
    while it < hsps.len() {
        let mut it2 = it + 1;
        while it2 < hsps.len() {
            let score12 = merge_score(&hsps[it], &hsps[it2]);
            if score12 > hsps[it].score.max(hsps[it2].score) {
                hsps[it] = merge(&hsps[it], &hsps[it2]);
                hsps.remove(it2);
            } else {
                let score21 = merge_score(&hsps[it2], &hsps[it]);
                if score21 > hsps[it].score.max(hsps[it2].score) {
                    hsps[it] = merge(&hsps[it2], &hsps[it]);
                    hsps.remove(it2);
                } else {
                    it2 += 1;
                }
            }
        }
        it += 1;
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HammingExtConfig {
    pub hamming_ext: bool,
    pub approx_min_id: f64,
    pub query_or_target_cover: f64,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub max_evalue: f64,
    pub diag_filter_cov: Option<f64>,
    pub diag_filter_id: Option<f64>,
}

impl Default for HammingExtConfig {
    fn default() -> Self {
        Self {
            hamming_ext: false,
            approx_min_id: 0.0,
            query_or_target_cover: 0.0,
            query_cover: 0.0,
            subject_cover: 0.0,
            max_evalue: f64::MAX,
            diag_filter_cov: None,
            diag_filter_id: None,
        }
    }
}

pub fn find_aln(
    segments: &mut [DiagonalSegment],
    qlen: i32,
    tlen: i32,
    config: &HammingExtConfig,
    score_matrix: &ScoreMatrix,
) -> ApproxHsp {
    for it in segments {
        let ev = score_matrix.evalue(it.score, qlen as u32, tlen as u32);
        if (it.id_percent() >= config.approx_min_id
            || stats::approx_id(it.score, it.len, it.len) >= config.approx_min_id)
            && ((config.query_or_target_cover > 0.0
                && it.cov_percent(qlen).max(it.cov_percent(tlen)) >= config.query_or_target_cover)
                || (config.query_or_target_cover == 0.0
                    && it.cov_percent(qlen) >= config.query_cover
                    && it.cov_percent(tlen) >= config.subject_cover))
            && ev <= config.max_evalue
        {
            return ApproxHsp::from_parts(
                0,
                0,
                it.score,
                0,
                it.query_range(),
                it.subject_range(),
                Anchor::from_diagonal_segment(it.clone()),
                ev,
            );
        }
    }
    ApproxHsp::new(0, 0)
}

pub fn filter(
    segments: &mut [DiagonalSegment],
    qlen: i32,
    tlen: i32,
    use_cov_filter: bool,
    config: &HammingExtConfig,
) -> ApproxHsp {
    const TOLERANCE_FACTOR: f64 = 1.1;
    const ID_MIN_COV: f64 = 80.0;
    segments.sort_by(|x, y| y.score.cmp(&x.score));
    let mut ident = 0;
    let mut len = 0;
    let qtol = (qlen as f64 * TOLERANCE_FACTOR) as i32;
    let ttol = (tlen as f64 * TOLERANCE_FACTOR) as i32;
    for it in segments {
        if len + it.len > qtol || len + it.len > ttol {
            continue;
        }
        ident += it.identities;
        len += it.len;
    }
    let qcov = len as f64 / qlen as f64 * 100.0;
    let tcov = len as f64 / tlen as f64 * 100.0;
    if let Some(diag_filter_cov) = config.diag_filter_cov {
        if use_cov_filter
            && ((config.query_or_target_cover > 0.0 && qcov.max(tcov) < diag_filter_cov)
                || (config.query_cover > 0.0 && qcov < diag_filter_cov)
                || (config.subject_cover > 0.0 && tcov < diag_filter_cov))
        {
            return ApproxHsp::new(0, -1);
        }
    }
    if let Some(diag_filter_id) = config.diag_filter_id {
        if qcov.max(tcov) >= ID_MIN_COV
            && len > 0
            && ident as f64 / len as f64 * 100.0 < diag_filter_id
        {
            return ApproxHsp::new(0, -1);
        }
    }
    ApproxHsp::new(0, 0)
}

pub fn hamming_ext(
    segments: &mut [DiagonalSegment],
    qlen: i32,
    tlen: i32,
    use_cov_filter: bool,
    config: &HammingExtConfig,
    score_matrix: &ScoreMatrix,
) -> ApproxHsp {
    if config.hamming_ext {
        let h = find_aln(segments, qlen, tlen, config, score_matrix);
        if h.score > 0 {
            return h;
        }
    }
    if config.diag_filter_cov.is_some() || config.diag_filter_id.is_some() {
        return filter(segments, qlen, tlen, use_cov_filter, config);
    }
    ApproxHsp::new(0, 0)
}

pub fn disjoint_approx_hsp(hsps: &[ApproxHsp], t: &ApproxHsp, cutoff: i32) -> bool {
    for begin in hsps {
        let ot = t.subject_range.overlap_factor(&begin.subject_range);
        let oq = t.query_range.overlap_factor(&begin.query_range);
        if (1.0 - ot.min(oq)) * t.score as f64 / begin.score as f64
            >= DEFAULT_CHAINING_STACKED_HSP_RATIO
        {
            continue;
        }
        if (1.0 - ot.max(oq)) * (t.score as f64) < cutoff as f64 {
            return false;
        }
    }
    true
}

pub fn disjoint_diagonal_segment(hsps: &[ApproxHsp], d: &DiagonalSegment, cutoff: i32) -> bool {
    for begin in hsps {
        let ot = d.subject_range().overlap_factor(&begin.subject_range);
        let oq = d.query_range().overlap_factor(&begin.query_range);
        if (1.0 - ot.min(oq)) * d.score as f64 / begin.score as f64
            >= DEFAULT_CHAINING_STACKED_HSP_RATIO
        {
            continue;
        }
        if (1.0 - ot.max(oq)) * (d.score as f64) < cutoff as f64 {
            return false;
        }
    }
    true
}

#[derive(Clone)]
pub struct Aligner<'a> {
    pub query: &'a [Letter],
    pub subject: &'a [Letter],
    pub log: bool,
    pub frame: u32,
    pub diags: DiagGraph,
    pub window: BTreeMap<i32, usize>,
    pub score_matrix: &'a ScoreMatrix,
}

impl<'a> Aligner<'a> {
    pub fn new(
        query: &'a [Letter],
        subject: &'a [Letter],
        log: bool,
        frame: u32,
        score_matrix: &'a ScoreMatrix,
    ) -> Self {
        Self {
            query,
            subject,
            log,
            frame,
            diags: DiagGraph::new(),
            window: BTreeMap::new(),
            score_matrix,
        }
    }

    pub fn get_approximate_link(
        &mut self,
        d_idx: usize,
        e_idx: usize,
        space_penalty: f64,
        _max_i: i32,
    ) -> i32 {
        let d = self.diags.nodes[d_idx].clone();
        let e = self.diags.nodes[e_idx].clone();
        let shift = d.diag() - e.diag();
        let gap_score = if shift != 0 {
            -self.score_matrix.gap_open() - shift.abs() * self.score_matrix.gap_extend()
        } else {
            0
        };
        let space = if shift > 0 {
            d.j - e.subject_last()
        } else {
            d.i - e.query_last()
        };
        let mut prefix_score = 0;
        let mut _link_score = 0;
        let mut link_j = 0;
        let mut path_max = 0;
        let mut path_min = 0;
        let mut prefix_score_begin = 0;
        if space <= 0 || space_penalty == 0.0 {
            if let Some(edge_idx) = self.diags.get_edge(d_idx, d.j) {
                let edge = &self.diags.edges[edge_idx];
                if edge.prefix_score > e.prefix_score + gap_score + d.score {
                    return 0;
                }
            }
            let mut link = Link::default();
            if get_link(
                &e.segment,
                &d.segment,
                self.query,
                self.subject,
                &mut link,
                LINK_PADDING,
                self.score_matrix,
            ) > 0
            {
                let diff1 = e.score - link.score1;
                let prefix = self.diags.prefix_score(e_idx, link.subject_pos1);
                let prefix_e = prefix.0;
                path_max = prefix.1;
                path_min = prefix.2;
                prefix_score = prefix_e - diff1 + gap_score + link.score2;
                if let Some(edge_idx) = self.diags.get_edge(d_idx, link.subject_pos2) {
                    if self.diags.edges[edge_idx].prefix_score > prefix_score {
                        return 0;
                    }
                }
                prefix_score_begin = prefix_score - link.score2;
                path_min = path_min.min(prefix_score - link.score2);
                if prefix_e == path_max {
                    path_max -= diff1;
                }
                _link_score = link.score1 + link.score2 + gap_score;
                link_j = link.subject_pos2;
            }
        } else {
            prefix_score = e.prefix_score + gap_score
                - (space_penalty * (space - 1).max(0) as f64) as i32
                + d.score;
            if let Some(edge_idx) = self.diags.get_edge(d_idx, d.j) {
                if self.diags.edges[edge_idx].prefix_score > prefix_score {
                    return 0;
                }
            }
            prefix_score_begin = prefix_score - d.score;
            path_max = e.path_max;
            path_min = e.path_min;
            path_min = path_min.min(prefix_score - d.score);
            _link_score = e.score + d.score + gap_score;
            link_j = d.j;
        }

        if prefix_score > d.score {
            path_max = path_max.max(prefix_score);
            self.diags.add_edge(Edge::new(
                prefix_score,
                path_max,
                link_j,
                d_idx as u32,
                e_idx as u32,
                if prefix_score == path_max {
                    prefix_score
                } else {
                    path_min
                },
                prefix_score_begin,
            ));
        }
        prefix_score
    }

    pub fn forward_pass<I>(&mut self, iter: I, init: bool, space_penalty: f64)
    where
        I: IntoIterator<Item = usize>,
    {
        self.window.clear();

        for node in iter {
            if init {
                let _ = self.diags.init_node(node);
            }
            let d = self.diags.nodes[node].clone();
            let dd = d.diag();
            self.window.entry(dd).or_insert(node);

            let mut max_j = 0;
            let left_keys: Vec<i32> = self.window.range(..dd).map(|(k, _)| *k).rev().collect();
            for key in left_keys {
                let Some(&e_node) = self.window.get(&key) else {
                    continue;
                };
                let e = self.diags.nodes[e_node].clone();
                if e.prefix_score - (space_penalty * (d.j - e.subject_end()).max(0) as f64) as i32
                    <= 0
                {
                    self.window.remove(&key);
                    continue;
                }
                if e.subject_end() < max_j {
                    continue;
                }
                self.get_approximate_link(node, e_node, space_penalty, max_j);
                max_j = max_j.max(d.j.min(e.subject_end()));
                if e.subject_end() - (d.subject_end() - (e.diag() - d.diag()).min(0))
                    >= REVERSE_LINK_MIN_OVERHANG
                {
                    self.get_approximate_link(e_node, node, space_penalty, max_j);
                }
            }

            let mut max_i = 0;
            let right_keys: Vec<i32> = self.window.range(dd..).map(|(k, _)| *k).collect();
            for key in right_keys {
                let Some(&e_node) = self.window.get(&key) else {
                    continue;
                };
                if e_node == node {
                    continue;
                }
                let e = self.diags.nodes[e_node].clone();
                if e.prefix_score - (space_penalty * (d.j - e.subject_end()).max(0) as f64) as i32
                    <= 0
                {
                    self.window.remove(&key);
                    continue;
                }
                if e.query_end() < max_i {
                    continue;
                }
                self.get_approximate_link(node, e_node, space_penalty, max_i);
                if e.i < d.i {
                    max_i = max_i.max(e.query_end().min(d.i));
                }
                if e.subject_end() - (d.subject_end() - (e.diag() - d.diag()).min(0))
                    >= REVERSE_LINK_MIN_OVERHANG
                {
                    self.get_approximate_link(e_node, node, space_penalty, max_i);
                }
            }
            self.window.insert(dd, node);
        }
    }

    pub fn backtrace_old(
        &self,
        node: usize,
        j_end: i32,
        t: &mut ApproxHsp,
        score_max: i32,
        mut score_min: i32,
        max_shift: i32,
        next: &mut u32,
    ) -> bool {
        let d = &self.diags.nodes[node];
        let edge_idx = self.diags.get_edge(node, j_end);
        let mut at_end = edge_idx.is_none();
        let prefix_score = edge_idx
            .map(|idx| self.diags.edges[idx].prefix_score)
            .unwrap_or(d.score);
        if prefix_score > score_max {
            return false;
        }

        score_min = score_min.min(
            edge_idx
                .map(|idx| self.diags.edges[idx].prefix_score_begin)
                .unwrap_or(0),
        );

        let mut j;
        if let Some(idx) = edge_idx {
            let f = &self.diags.edges[idx];
            let e = &self.diags.nodes[f.node_out as usize];
            let shift = d.diag() - e.diag();
            j = f.j;

            if shift.abs() <= max_shift {
                let bt = self.backtrace_old(
                    f.node_out as usize,
                    if shift > 0 { j } else { j + shift },
                    t,
                    score_max,
                    score_min,
                    max_shift,
                    next,
                );
                if !bt {
                    if f.prefix_score_begin > score_min {
                        return false;
                    } else {
                        at_end = true;
                    }
                }
            } else {
                *next = f.node_out;
                at_end = true;
            }
        } else {
            j = d.j;
        }

        if at_end {
            t.query_range.begin = d.i;
            t.subject_range.begin = d.j;
            t.score = score_max - score_min;
            j = d.j;
        } else if let Some(idx) = edge_idx {
            j = self.diags.edges[idx].j;
        }

        let dd = d.diag();
        t.d_max = t.d_max.max(dd);
        t.d_min = t.d_min.min(dd);
        if d.score > t.max_diag.segment.score {
            t.max_diag = crate::util::hsp::Anchor::from_diagonal_segment(d.segment.clone());
            t.max_diag.prefix_score = prefix_score;
            t.max_diag.d_max_left = t.max_diag.d_max_right.max(t.max_diag.d_max_left).max(dd);
            t.max_diag.d_min_left = t.max_diag.d_min_right.min(t.max_diag.d_min_left).min(dd);
            t.max_diag.d_max_right = dd;
            t.max_diag.d_min_right = dd;
        } else {
            t.max_diag.d_max_right = t.max_diag.d_max_right.max(dd);
            t.max_diag.d_min_right = t.max_diag.d_min_right.min(dd);
        }

        let _ = j;
        true
    }

    pub fn backtrace(
        &self,
        top_node: usize,
        t: &mut ApproxHsp,
        max_shift: i32,
        next: &mut u32,
        max_j: i32,
    ) {
        let mut traits = ApproxHsp::new(self.frame, 0);
        if top_node != END {
            let d = &self.diags.nodes[top_node];
            traits.subject_range.end = d.subject_end();
            traits.query_range.end = d.query_end();
            let score_min = d.prefix_score;
            self.backtrace_old(
                top_node,
                d.subject_end().min(max_j),
                &mut traits,
                d.prefix_score,
                score_min,
                max_shift,
                next,
            );
        } else {
            traits.score = 0;
        }
        *t = traits;
    }

    pub fn backtrace_top(
        &self,
        top_node: usize,
        ts: &mut Vec<ApproxHsp>,
        t_begin: &mut usize,
        cutoff: i32,
        max_shift: i32,
    ) -> i32 {
        let mut top_node = top_node;
        let mut max_score = 0;
        let mut max_j = self.subject.len() as i32;
        loop {
            let mut t = ApproxHsp::new(self.frame, 0);
            let mut next = u32::MAX;
            self.backtrace(top_node, &mut t, max_shift, &mut next, max_j);
            if t.score > 0 {
                max_j = t.subject_range.begin;
            }
            if t.score >= cutoff && disjoint_approx_hsp(&ts[*t_begin..], &t, cutoff) {
                let was_end = *t_begin == ts.len();
                max_score = max_score.max(t.score);
                ts.push(t);
                if was_end {
                    *t_begin = ts.len() - 1;
                }
            }
            if next == u32::MAX {
                break;
            }
            top_node = next as usize;
        }
        max_score
    }

    pub fn backtrace_all(&self, ts: &mut Vec<ApproxHsp>, cutoff: i32, max_shift: i32) -> i32 {
        let mut top_nodes: Vec<usize> = self
            .diags
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, d)| (d.rel_score() >= cutoff).then_some(i))
            .collect();
        top_nodes.sort_by(|&x, &y| {
            self.diags.nodes[y]
                .rel_score()
                .cmp(&self.diags.nodes[x].rel_score())
        });

        let mut max_score = 0;
        let mut t_begin = ts.len();
        for node in top_nodes {
            if disjoint_diagonal_segment(&ts[t_begin..], &self.diags.nodes[node].segment, cutoff) {
                max_score =
                    max_score.max(self.backtrace_top(node, ts, &mut t_begin, cutoff, max_shift));
            }
        }
        max_score
    }

    pub fn run(
        &mut self,
        ts: &mut Vec<ApproxHsp>,
        space_penalty: f64,
        cutoff: i32,
        max_shift: i32,
    ) -> i32 {
        self.diags.sort();
        self.diags.prune();
        self.forward_pass(0..self.diags.nodes.len(), true, space_penalty);
        self.backtrace_all(ts, cutoff, max_shift)
    }

    pub fn run_segments(
        &mut self,
        ts: &mut Vec<ApproxHsp>,
        segments: &[DiagonalSegment],
        band: i32,
    ) -> i32 {
        self.diags.init();
        self.diags.load(segments);
        self.run(ts, SPACE_PENALTY, 19, band)
    }
}

pub fn run(
    query: &[Letter],
    subject: &[Letter],
    segments: &[DiagonalSegment],
    log: bool,
    frame: u32,
    band: i32,
    score_matrix: &ScoreMatrix,
    no_chaining_merge_hsps: bool,
) -> (i32, Vec<ApproxHsp>) {
    if segments.len() == 1 {
        let d = segments[0].diag();
        let anchor = Anchor::new(segments[0].clone(), d, d, d, d, segments[0].score);
        return (
            segments[0].score,
            vec![ApproxHsp::from_parts(
                d,
                d,
                segments[0].score,
                frame as i32,
                segments[0].query_range(),
                segments[0].subject_range(),
                anchor,
                f64::MAX,
            )],
        );
    }
    let mut aligner = Aligner::new(query, subject, log, frame, score_matrix);
    let mut ts = Vec::new();
    let score = aligner.run_segments(&mut ts, segments, band);
    if !no_chaining_merge_hsps {
        merge_hsps(&mut ts);
    }
    (score, ts)
}

pub fn run_targets(_query: &[Letter], _targets: &[DpTarget]) -> Vec<crate::align::hsp::Hsp> {
    Vec::new()
}

pub fn print_diag<W: Write>(
    i0: i32,
    j0: i32,
    l: i32,
    score: i32,
    diags: &DiagGraph,
    query: &[Letter],
    subject: &[Letter],
    score_matrix: &ScoreMatrix,
    out: &mut W,
) -> io::Result<()> {
    let ds = DiagonalSegment::new(i0, j0, l, 0);
    let mut n = 0usize;
    for (idx, d) in diags.nodes.iter().enumerate() {
        if d.intersect(&ds).len > 0 {
            if d.score == 0 {
                continue;
            }
            let diff = score_range(
                query,
                subject,
                d.query_end() as usize,
                d.subject_end() as usize,
                (j0 + l) as usize,
                score_matrix,
            );
            if n > 0 {
                write!(out, "(")?;
            }
            let prefix_score = score
                + score_range(
                    query,
                    subject,
                    (i0 + l) as usize,
                    (j0 + l) as usize,
                    d.subject_end() as usize,
                    score_matrix,
                )
                - diff.min(0);
            let prefix_score2 = diags.prefix_score(idx, j0 + l).0;
            write!(
                out,
                "Diag n={} i={} j={} len={} prefix_score={} prefix_score2={}",
                idx, i0, j0, l, prefix_score, prefix_score2
            )?;
            if n > 0 {
                write!(out, ")")?;
            }
            writeln!(out)?;
            n += 1;
        }
    }
    if n == 0 {
        writeln!(
            out,
            "Diag n=x i={} j={} len={} prefix_score={}",
            i0, j0, l, score
        )?;
    }
    Ok(())
}

pub fn smith_waterman<W: Write>(
    query: &[Letter],
    subject: &[Letter],
    diags: &DiagGraph,
    score_matrix: &ScoreMatrix,
    out: &mut W,
) -> io::Result<()> {
    let hsp = crate::dp::smith_waterman::smith_waterman(query, subject, score_matrix);
    let mut qpos = hsp.query_begin;
    let mut spos = hsp.subject_begin;
    let mut i0 = -1;
    let mut j0 = -1;
    let mut l = 0;
    let mut score = 0;

    for (op, count) in hsp.operations {
        match op {
            EditOperation::Match | EditOperation::Substitution => {
                for _ in 0..count {
                    if i0 < 0 {
                        i0 = qpos;
                        j0 = spos;
                        l = 0;
                    }
                    score += score_matrix.score(query[qpos as usize], subject[spos as usize]);
                    l += 1;
                    qpos += 1;
                    spos += 1;
                }
            }
            EditOperation::Deletion => {
                for _ in 0..count {
                    if i0 >= 0 {
                        print_diag(i0, j0, l, score, diags, query, subject, score_matrix, out)?;
                        score -= score_matrix.gap_open() + score_matrix.gap_extend();
                        i0 = -1;
                        j0 = -1;
                    } else {
                        score -= score_matrix.gap_extend();
                    }
                    spos += 1;
                }
            }
            EditOperation::Insertion => {
                for _ in 0..count {
                    if i0 >= 0 {
                        print_diag(i0, j0, l, score, diags, query, subject, score_matrix, out)?;
                        score -= score_matrix.gap_open() + score_matrix.gap_extend();
                        i0 = -1;
                        j0 = -1;
                    } else {
                        score -= score_matrix.gap_extend();
                    }
                    qpos += 1;
                }
            }
            EditOperation::FrameshiftForward | EditOperation::FrameshiftReverse => {}
        }
    }
    print_diag(i0, j0, l, score, diags, query, subject, score_matrix, out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::score_matrix::ScoreMatrix;
    use crate::util::hsp::Anchor;

    #[test]
    fn test_diagonal_node_methods() {
        let segment = DiagonalSegment::new(2, 5, 4, 20);
        let mut node = DiagonalNode::from_segment(segment.clone());
        assert_eq!(node.segment, segment);
        assert!(node.is_maximum());
        assert_eq!(node.rel_score(), 20);
        node.path_max = 25;
        node.path_min = 12;
        node.prefix_score = 18;
        assert_eq!(node.rel_score(), 6);
        node.deactivate();
        assert_eq!(node.link_idx, 0);
        node.reset();
        assert_eq!(node.link_idx, -1);
        assert_eq!(node.prefix_score, 20);
    }

    #[test]
    fn test_diag_graph_load_sort_edges_and_top() {
        let segments = vec![
            DiagonalSegment::new(0, 0, 4, 10),
            DiagonalSegment::new(2, 2, 3, 8),
            DiagonalSegment::new(7, 4, 2, 20),
            DiagonalSegment::new(8, 6, 2, 15),
        ];
        let mut graph = DiagGraph::new();
        graph.load(&segments);
        assert_eq!(graph.nodes.len(), 3);
        graph.sort();
        assert!(graph.nodes.windows(2).all(|w| w[0].j <= w[1].j));
        graph.init_node(1).unwrap();
        graph.add_edge(Edge::new(30, 30, graph.nodes[1].j + 1, 1, 0, 10, 20));
        assert_eq!(graph.top_node(), 1);
        assert_eq!(graph.prefix_score(1, graph.nodes[1].j + 2), (30, 30, 10));
        graph.clear_edges();
        assert!(graph.edges.is_empty());
        assert!(graph.nodes.iter().all(|n| n.link_idx == 0));
    }

    #[test]
    fn test_diag_graph_prune_default_cover_keeps_noncovered() {
        let mut graph = DiagGraph::new();
        for k in 0..10 {
            graph
                .nodes
                .push(DiagonalNode::from_segment(DiagonalSegment::new(
                    k,
                    k,
                    10 - k,
                    10 + k,
                )));
        }
        graph.prune();
        assert!(!graph.nodes.is_empty());
        assert!(graph.nodes.len() <= 10);
    }

    #[test]
    fn test_link_transpose_reset_and_gap_links() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let query = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let subject = query.clone();
        let d1 = DiagonalSegment::new(0, 0, 4, sm.score(0, 0) * 4);
        let d2 = DiagonalSegment::new(5, 4, 4, sm.score(0, 0) * 4);
        let mut link = Link::default();
        let score = get_link(&d1, &d2, &query, &subject, &mut link, LINK_PADDING, &sm);
        assert!(score > i32::MIN);
        assert!(link.subject_pos1 >= 0);
        let old_subject = link.subject_pos1;
        let old_query = link.query_pos1;
        link.transpose();
        assert_eq!(link.subject_pos1, old_query);
        assert_eq!(link.query_pos1, old_subject);
        link.reset();
        assert_eq!(link.subject_pos1, -1);
        assert_eq!(link.score1, 0);
    }

    #[test]
    fn test_merge_score_merge_and_merge_hsps() {
        let h1 = ApproxHsp::from_parts(
            0,
            0,
            30,
            1,
            Interval::new(0, 10),
            Interval::new(0, 10),
            Anchor::new(DiagonalSegment::new(0, 0, 10, 30), 0, 0, 0, 0, 30),
            1.0,
        );
        let h2 = ApproxHsp::from_parts(
            1,
            1,
            25,
            1,
            Interval::new(12, 20),
            Interval::new(11, 19),
            Anchor::new(DiagonalSegment::new(12, 11, 8, 25), 1, 1, 1, 1, 25),
            1.0,
        );
        assert_eq!(merge_score(&h1, &h2), 53);
        assert_eq!(merge_score(&h2, &h1), 0);

        let merged = merge(&h1, &h2);
        assert_eq!(merged.d_min, 0);
        assert_eq!(merged.d_max, 1);
        assert_eq!(merged.query_range, Interval::new(0, 20));
        assert_eq!(merged.subject_range, Interval::new(0, 19));
        assert_eq!(merged.score, 53);
        assert_eq!(merged.max_diag.d_max_right, 1);

        let mut hsps = vec![h1, h2];
        merge_hsps(&mut hsps);
        assert_eq!(hsps.len(), 1);
        assert_eq!(hsps[0].score, 53);
    }

    #[test]
    fn test_aligner_forward_pass_adds_spaced_link() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0; 16];
        let subject = vec![0; 16];
        let score = sm.score(0, 0) * 4;
        let segments = vec![
            DiagonalSegment::new(0, 0, 4, score),
            DiagonalSegment::new(5, 5, 4, score),
        ];

        let mut aligner = Aligner::new(&query, &subject, false, 0, &sm);
        aligner.diags.load(&segments);
        aligner.forward_pass(0..aligner.diags.nodes.len(), true, SPACE_PENALTY);

        assert_eq!(aligner.diags.edges.len(), 1);
        assert_eq!(aligner.diags.edges[0].node_in, 1);
        assert_eq!(aligner.diags.edges[0].node_out, 0);
        assert!(aligner.diags.nodes[1].prefix_score > score);
    }

    #[test]
    fn test_backtrace_all_collects_approx_hsp() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0; 16];
        let subject = vec![0; 16];
        let score = sm.score(0, 0) * 4;
        let segments = vec![
            DiagonalSegment::new(0, 0, 4, score),
            DiagonalSegment::new(5, 5, 4, score),
        ];

        let mut aligner = Aligner::new(&query, &subject, false, 3, &sm);
        aligner.diags.load(&segments);
        aligner.forward_pass(0..aligner.diags.nodes.len(), true, SPACE_PENALTY);

        let mut ts = Vec::new();
        let max_score = aligner.backtrace_all(&mut ts, 1, 100);
        assert_eq!(max_score, aligner.diags.nodes[1].prefix_score);
        assert_eq!(ts.len(), 1);
        assert_eq!(ts[0].frame, 3);
        assert_eq!(ts[0].query_range, Interval::new(0, 9));
        assert_eq!(ts[0].subject_range, Interval::new(0, 9));
        assert!(disjoint_approx_hsp(&[], &ts[0], 1));
        assert!(!disjoint_diagonal_segment(&ts, &segments[0], score));
    }

    #[test]
    fn test_hamming_ext_find_aln_and_filters() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let mut segments = vec![DiagonalSegment::with_identities(0, 0, 10, 60, 9)];
        let config = HammingExtConfig {
            hamming_ext: true,
            approx_min_id: 80.0,
            query_cover: 50.0,
            subject_cover: 50.0,
            max_evalue: f64::MAX,
            ..HammingExtConfig::default()
        };

        let h = hamming_ext(&mut segments, 20, 20, true, &config, &sm);
        assert_eq!(h.score, 60);
        assert_eq!(h.query_range, Interval::new(0, 10));
        assert_eq!(h.subject_range, Interval::new(0, 10));

        let mut low_identity = vec![
            DiagonalSegment::with_identities(0, 0, 6, 10, 1),
            DiagonalSegment::with_identities(6, 6, 6, 10, 1),
        ];
        let reject_id = HammingExtConfig {
            diag_filter_id: Some(50.0),
            ..HammingExtConfig::default()
        };
        let h = hamming_ext(&mut low_identity, 12, 12, true, &reject_id, &sm);
        assert_eq!(h.score, -1);

        let mut short = vec![DiagonalSegment::with_identities(0, 0, 2, 10, 2)];
        let reject_cov = HammingExtConfig {
            query_cover: 80.0,
            subject_cover: 80.0,
            diag_filter_cov: Some(80.0),
            ..HammingExtConfig::default()
        };
        let h = hamming_ext(&mut short, 10, 10, true, &reject_cov, &sm);
        assert_eq!(h.score, -1);
    }

    #[test]
    fn test_run_single_and_multi_segment() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0; 16];
        let subject = vec![0; 16];
        let score = sm.score(0, 0) * 4;
        let single = vec![DiagonalSegment::new(2, 1, 4, score)];

        let (single_score, single_hsps) = run(&query, &subject, &single, false, 2, 100, &sm, false);
        assert_eq!(single_score, score);
        assert_eq!(single_hsps.len(), 1);
        assert_eq!(single_hsps[0].frame, 2);
        assert_eq!(single_hsps[0].d_min, single[0].diag());
        assert_eq!(single_hsps[0].d_max, single[0].diag());

        let segments = vec![
            DiagonalSegment::new(0, 0, 4, score),
            DiagonalSegment::new(5, 5, 4, score),
        ];
        let (multi_score, multi_hsps) = run(&query, &subject, &segments, false, 1, 100, &sm, true);
        assert!(multi_score > score);
        assert_eq!(multi_hsps.len(), 1);
        assert_eq!(multi_hsps[0].frame, 1);
    }

    #[test]
    fn test_print_diag_and_smith_waterman_debug_output() {
        let sm = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0; 8];
        let subject = vec![0; 8];
        let score = sm.score(0, 0) * 4;
        let mut graph = DiagGraph::new();
        graph.load(&[DiagonalSegment::new(0, 0, 4, score)]);

        let mut out = Vec::new();
        print_diag(0, 0, 4, score, &graph, &query, &subject, &sm, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("Diag n=0 i=0 j=0 len=4"));
        assert!(text.contains("prefix_score="));

        let mut out = Vec::new();
        smith_waterman(&query[..4], &subject[..4], &graph, &sm, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("Diag n=0 i=0 j=0 len=4"));
    }
}
