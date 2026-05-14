use crate::basic::value::{BlockId, Loc};
use crate::dna::seed_set_dna::SeedMatch;

pub const MIN_OVERLAP_PERCENTAGE_SECONDARY: f64 = 0.5;

#[derive(Debug, Clone)]
pub struct ChainingParameters {
    pub max_dist_x: i32,
    pub max_dist_y: i32,
    pub band_width: i32,
    pub max_skip: i32,
    pub max_iterations: i32,
    pub map_percentage_target: f32,
    pub min_chain_score: i32,
    pub chain_pen_gap: f32,
    pub chain_pen_skip: f32,
    pub max_overlap_extension: f32,
    pub best_hsp_only: bool,
}

impl ChainingParameters {
    pub fn new(gap: f32, skip: f32, min_chain_score: i32, max_overlap_extension: f32) -> Self {
        Self {
            max_dist_x: 1000,
            max_dist_y: 1000,
            band_width: 300,
            max_skip: 25,
            max_iterations: 3000,
            map_percentage_target: 0.99,
            min_chain_score,
            chain_pen_gap: gap,
            chain_pen_skip: skip,
            max_overlap_extension,
            best_hsp_only: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnchorData {
    pub predecessor_anchor: Vec<i64>,
    pub best_score_anchor: Vec<i32>,
    pub peak_score_anchor: Vec<i32>,
    pub pre_predecessor_anchor: Vec<i32>,
    pub anchor_used: Vec<bool>,
}

impl AnchorData {
    pub fn new(n: usize) -> Self {
        Self {
            predecessor_anchor: vec![0; n],
            best_score_anchor: vec![0; n],
            peak_score_anchor: vec![0; n],
            pre_predecessor_anchor: vec![0; n],
            anchor_used: vec![false; n],
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chain {
    pub chain_score: i32,
    pub target_id: BlockId,
    pub mapping_quality: u8,
    pub is_primary: u8,
    pub reverse: bool,
    pub anchors: Vec<Anchor>,
}

impl Chain {
    pub fn new(rev: bool) -> Self {
        Self {
            chain_score: 0,
            target_id: 0,
            mapping_quality: 0,
            is_primary: 0,
            reverse: rev,
            anchors: Vec::new(),
        }
    }

    pub fn overlap_in_query(&self, other_chain: &Chain) -> i32 {
        self.anchors[0].i.min(other_chain.anchors[0].i)
            - self.anchors[self.anchors.len() - 1]
                .i_start()
                .max(other_chain.anchors[other_chain.anchors.len() - 1].i_start())
    }

    pub fn overlap_in_target(&self, other_chain: &Chain) -> i32 {
        self.anchors[0].j.min(other_chain.anchors[0].j)
            - self.anchors[self.anchors.len() - 1]
                .j_start()
                .max(other_chain.anchors[other_chain.anchors.len() - 1].j_start())
    }

    pub fn compute_mapping_quality(&mut self, score_secondary_chain: i32) {
        let score_ratio = f64::from(score_secondary_chain) / f64::from(self.chain_score);
        let quality_score = 40.0
            * (1.0 - score_ratio)
            * 1.0_f64.min(self.anchors.len() as f64 / 10.0)
            * f64::from(self.chain_score).ln();
        self.mapping_quality = (quality_score * 60.0 / 312.0) as u8;
    }
}

impl PartialOrd for Chain {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Chain {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.chain_score.cmp(&other.chain_score)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anchor {
    pub i: Loc,
    pub j: Loc,
    pub span: i32,
}

impl Anchor {
    pub fn new(i: Loc, j: Loc, span: i32) -> Self {
        Self { i, j, span }
    }

    pub fn i_start(&self) -> Loc {
        self.i - self.span
    }

    pub fn j_start(&self) -> Loc {
        self.j - self.span
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScoreIndexPair<T1, T2> {
    pub score: T1,
    pub index: T2,
}

impl<T1, T2> ScoreIndexPair<T1, T2>
where
    T1: Default,
    T2: Default,
{
    pub fn new() -> Self {
        Self {
            score: T1::default(),
            index: T2::default(),
        }
    }
}

impl<T1, T2> ScoreIndexPair<T1, T2> {
    pub fn from_values(first: T1, second: T2) -> Self {
        Self {
            score: first,
            index: second,
        }
    }
}

impl<T1, T2> PartialOrd for ScoreIndexPair<T1, T2>
where
    T1: Ord,
    T2: Ord,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<T1, T2> Ord for ScoreIndexPair<T1, T2>
where
    T1: Ord,
    T2: Ord,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score
            .cmp(&other.score)
            .then_with(|| self.index.cmp(&other.index))
    }
}

pub fn only_keep_best_chains_per_target(chains: &mut Vec<Chain>, chain_cutoff_percentage: f32) {
    if chains.is_empty() {
        return;
    }

    chains.sort_by(|a, b| {
        if a.target_id == b.target_id {
            b.chain_score.cmp(&a.chain_score)
        } else {
            a.target_id.cmp(&b.target_id)
        }
    });

    let mut last_target_id = chains[0].target_id;
    chains[0].is_primary = 1;
    let mut cutoff = (chain_cutoff_percentage * chains[0].chain_score as f32) as i32;

    for i in 1..chains.len() {
        if chains[i].chain_score >= cutoff && chains[i].target_id == last_target_id {
            chains[i].is_primary = 1;
        } else if chains[i].target_id != last_target_id {
            last_target_id = chains[i].target_id;
            chains[i].is_primary = 1;
            cutoff = (chain_cutoff_percentage * chains[i].chain_score as f32) as i32;
        }
    }

    chains.sort_by(|a, b| b.is_primary.cmp(&a.is_primary));

    let first_secondary = chains
        .iter()
        .position(|chain| chain.is_primary == 0)
        .unwrap_or(chains.len());
    chains.truncate(first_secondary);
}

pub fn detect_primary_chains(chains: &mut [Chain]) {
    if chains.is_empty() {
        return;
    }

    let mut best_secondary_score_of_primary = vec![0; chains.len()];
    let mut primary_chain_indices = vec![0usize];

    let mut chain_span = Vec::with_capacity(chains.len());
    for chain in chains.iter() {
        chain_span.push(chain.anchors[0].i - chain.anchors[chain.anchors.len() - 1].i_start());
    }

    for index_chain in 1..chains.len() {
        let mut is_primary = true;

        for &index_primary_chain in &primary_chain_indices {
            let overlap_length = chains[index_chain].overlap_in_query(&chains[index_primary_chain]);
            if overlap_length < 1 {
                continue;
            }

            let overlap_percentage = f64::from(overlap_length)
                / f64::from(chain_span[index_chain].min(chain_span[index_primary_chain]));

            if overlap_percentage >= MIN_OVERLAP_PERCENTAGE_SECONDARY {
                is_primary = false;
                best_secondary_score_of_primary[index_primary_chain] =
                    best_secondary_score_of_primary[index_primary_chain]
                        .max(chains[index_chain].chain_score);
            }
        }
        if is_primary {
            primary_chain_indices.push(index_chain);
        }
    }

    for i in primary_chain_indices {
        chains[i].compute_mapping_quality(best_secondary_score_of_primary[i]);
    }
}

pub fn find_chain_start(
    max_drop: i32,
    score_end: u64,
    index_end: u64,
    anchor_data: &AnchorData,
) -> i64 {
    let mut index = index_end as i64;
    let mut index_max_score = index;
    let mut max_score = 0;
    if index < 0 || anchor_data.anchor_used[index as usize] {
        return index;
    }

    loop {
        index = anchor_data.predecessor_anchor[index as usize];
        let score_difference_anchors = if index < 0 {
            score_end as i32
        } else {
            score_end as i32 - anchor_data.best_score_anchor[index as usize]
        };
        if score_difference_anchors > max_score {
            max_score = score_difference_anchors;
            index_max_score = index;
        } else if max_score - score_difference_anchors > max_drop {
            break;
        }
        if index <= -1 || anchor_data.anchor_used[index as usize] {
            break;
        }
    }
    index_max_score
}

pub fn chain_backtrack(
    anchor_data: &mut AnchorData,
    min_chain_score: i32,
    max_drop: i32,
    seed_match_begin: &[SeedMatch],
    is_reverse: bool,
    best_hsp_only: bool,
) -> Vec<Chain> {
    let mut potential_chain_ends = Vec::new();
    for i in 0..anchor_data.best_score_anchor.len() {
        if anchor_data.best_score_anchor[i] >= min_chain_score {
            potential_chain_ends.push(ScoreIndexPair::from_values(
                anchor_data.best_score_anchor[i] as u64,
                i as u64,
            ));
        }
    }

    if potential_chain_ends.is_empty() {
        return Vec::new();
    }

    if best_hsp_only {
        potential_chain_ends.sort_by(|a, b| b.score.cmp(&a.score));
        if potential_chain_ends.len() > 4 {
            potential_chain_ends.truncate(4);
        }
    }

    potential_chain_ends.sort();

    let mut detected_chains = Vec::new();
    anchor_data
        .anchor_used
        .resize(anchor_data.anchor_used.len(), false);
    anchor_data.anchor_used.fill(false);

    for k in (0..potential_chain_ends.len()).rev() {
        if anchor_data.anchor_used[potential_chain_ends[k].index as usize] {
            continue;
        }

        let chain_start_index = find_chain_start(
            max_drop,
            potential_chain_ends[k].score,
            potential_chain_ends[k].index,
            anchor_data,
        );

        let mut chain = Chain::new(is_reverse);
        let mut idx_chain_end = potential_chain_ends[k].index as i64;
        while idx_chain_end != chain_start_index {
            anchor_data.anchor_used[idx_chain_end as usize] = true;
            let seed = seed_match_begin[idx_chain_end as usize];
            chain
                .anchors
                .push(Anchor::new(seed.i(), seed.j(), seed.ungapped_score()));
            idx_chain_end = anchor_data.predecessor_anchor[idx_chain_end as usize];
        }

        let score = if idx_chain_end < 0 {
            potential_chain_ends[k].score as i32
        } else {
            potential_chain_ends[k].score as i32
                - anchor_data.best_score_anchor[idx_chain_end as usize]
        };

        if score >= min_chain_score && !chain.anchors.is_empty() {
            chain.target_id = seed_match_begin[0].id();
            chain.chain_score = score;
            detected_chains.push(chain);
        }
    }
    detected_chains
}

pub fn compute_score(
    second_match: &SeedMatch,
    first_match: &SeedMatch,
    chaining_parameters: &ChainingParameters,
) -> i32 {
    let distance_query = second_match.i_start() - first_match.i();
    let distance_query_end_to_end = second_match.i() - first_match.i();
    if distance_query_end_to_end < 1 || distance_query > chaining_parameters.max_dist_x {
        return i32::MIN;
    }

    let distance_target = second_match.j_start() - first_match.j();
    let distance_target_end_to_end = second_match.j() - first_match.j();
    if distance_target_end_to_end == 0 || distance_target > chaining_parameters.max_dist_y {
        return i32::MIN;
    }

    let distance_diagonal = if distance_target_end_to_end > distance_query_end_to_end {
        distance_target_end_to_end - distance_query_end_to_end
    } else {
        distance_query_end_to_end - distance_target_end_to_end
    };

    if distance_diagonal > chaining_parameters.band_width {
        return i32::MIN;
    }

    let distance_skip = distance_target.abs().min(distance_query.abs());
    let distance_gap_end_to_end = distance_target_end_to_end.min(distance_query_end_to_end);

    let mut score = second_match.ungapped_score().min(distance_gap_end_to_end);
    if distance_diagonal != 0 {
        let lin_pen = chaining_parameters.chain_pen_gap * distance_diagonal as f32
            + chaining_parameters.chain_pen_skip * distance_skip as f32;
        let log_pen = if distance_diagonal > 0 {
            log2_approximate((distance_diagonal + 1) as f32)
        } else {
            0.0
        };
        score -= (lin_pen + 0.5 * log_pen) as i32;
    }

    score
}

pub fn chaining_dynamic_program(
    chaining_parameters: &ChainingParameters,
    seed_match_begin: &[SeedMatch],
    is_reverse: bool,
) -> Vec<Chain> {
    let total_number_matches = seed_match_begin.len();
    let mut anchor_data = AnchorData::new(total_number_matches);
    let mut total_max_score = 0;

    let mut max_score_index = -1i64;
    for index_second_match in 0..total_number_matches {
        let mut index_predecessor = -1i64;
        let mut max_score = seed_match_begin[index_second_match].ungapped_score();
        let mut n_skip = 0;

        let mut start = 0usize;
        while start < index_second_match
            && seed_match_begin[index_second_match].j_start()
                > seed_match_begin[start].j() + chaining_parameters.max_dist_x
        {
            start += 1;
        }

        let max_iterations_start =
            index_second_match.saturating_sub(chaining_parameters.max_iterations as usize);
        start = start.max(max_iterations_start);

        let mut index_first_match = index_second_match as i64 - 1;
        while index_first_match >= start as i64 {
            let score = compute_score(
                &seed_match_begin[index_second_match],
                &seed_match_begin[index_first_match as usize],
                chaining_parameters,
            );
            if score != i32::MIN {
                let score = score + anchor_data.best_score_anchor[index_first_match as usize];

                if score > max_score {
                    max_score = score;
                    index_predecessor = index_first_match;
                    if n_skip > 0 {
                        n_skip -= 1;
                    }
                } else if anchor_data.pre_predecessor_anchor[index_first_match as usize]
                    == index_second_match as i32
                {
                    n_skip += 1;
                    if n_skip > chaining_parameters.max_skip {
                        break;
                    }
                }

                if anchor_data.predecessor_anchor[index_first_match as usize] > -1 {
                    let pre = anchor_data.predecessor_anchor[index_first_match as usize] as usize;
                    anchor_data.pre_predecessor_anchor[pre] = index_second_match as i32;
                }
            }
            index_first_match -= 1;
        }
        let end_j = index_first_match;
        if max_score_index < 0
            || seed_match_begin[index_second_match].j_start()
                - seed_match_begin[max_score_index as usize].j()
                > chaining_parameters.max_dist_x
        {
            let mut max = i32::MIN;
            max_score_index = -1;
            index_first_match = index_second_match as i64 - 1;
            while index_first_match >= start as i64 {
                if max < anchor_data.best_score_anchor[index_first_match as usize] {
                    max = anchor_data.best_score_anchor[index_first_match as usize];
                    max_score_index = index_first_match;
                }
                index_first_match -= 1;
            }
        }

        if max_score_index >= 0 && max_score_index < end_j {
            let score_extend_max = compute_score(
                &seed_match_begin[index_second_match],
                &seed_match_begin[max_score_index as usize],
                chaining_parameters,
            );
            if score_extend_max != i32::MIN
                && max_score
                    < score_extend_max + anchor_data.best_score_anchor[max_score_index as usize]
            {
                max_score =
                    score_extend_max + anchor_data.best_score_anchor[max_score_index as usize];
                index_predecessor = max_score_index;
            }
        }

        anchor_data.best_score_anchor[index_second_match] = max_score;
        anchor_data.predecessor_anchor[index_second_match] = index_predecessor;
        anchor_data.peak_score_anchor[index_second_match] = if index_predecessor > -1
            && anchor_data.peak_score_anchor[index_predecessor as usize] > max_score
        {
            anchor_data.peak_score_anchor[index_predecessor as usize]
        } else {
            max_score
        };
        if max_score_index < 0
            || anchor_data.best_score_anchor[max_score_index as usize]
                < anchor_data.best_score_anchor[index_second_match]
        {
            max_score_index = index_second_match as i64;
        }
        total_max_score = total_max_score.max(max_score);
    }

    chain_backtrack(
        &mut anchor_data,
        chaining_parameters.min_chain_score,
        chaining_parameters.band_width,
        seed_match_begin,
        is_reverse,
        chaining_parameters.best_hsp_only,
    )
}

pub fn log2_approximate(x: f32) -> f32 {
    let mut i = x.to_bits();
    let mut log_2 = ((i >> 23) & 255) as f32 - 128.0;
    i &= !(255 << 23);
    i += 127 << 23;
    let f = f32::from_bits(i);
    log_2 += (-0.34484843 * f + 2.02466578) * f - 0.67487759;
    log_2
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seed(i: Loc, id: BlockId, j: Loc, score: i32) -> SeedMatch {
        let mut s = SeedMatch::new(i, id, j, score);
        s.score(score);
        s
    }

    fn chain(target_id: BlockId, score: i32, anchors: &[(Loc, Loc, i32)]) -> Chain {
        let mut c = Chain::new(false);
        c.target_id = target_id;
        c.chain_score = score;
        c.anchors = anchors
            .iter()
            .map(|&(i, j, span)| Anchor::new(i, j, span))
            .collect();
        c
    }

    #[test]
    fn test_seed_match_accessors_and_ordering() {
        let mut s = SeedMatch::new(17, 3, 29, 11);
        assert_eq!(s.ungapped_score(), 11);
        s.score(7);
        assert_eq!(s.i(), 17);
        assert_eq!(s.j(), 29);
        assert_eq!(s.id(), 3);
        assert_eq!(s.i_start(), 10);
        assert_eq!(s.j_start(), 22);

        assert!(SeedMatch::new(1, 1, 2, 4) > SeedMatch::new(1, 2, 2, 4));
        assert!(SeedMatch::new(1, 1, 2, 8) > SeedMatch::new(1, 1, 2, 4));
    }

    #[test]
    fn test_only_keep_best_chains_per_target() {
        let mut chains = vec![
            chain(2, 50, &[(50, 50, 10)]),
            chain(1, 40, &[(40, 40, 10)]),
            chain(1, 20, &[(20, 20, 10)]),
            chain(2, 35, &[(35, 35, 10)]),
            chain(3, 8, &[(8, 8, 4)]),
        ];

        only_keep_best_chains_per_target(&mut chains, 0.75);

        let mut kept = chains
            .iter()
            .map(|c| (c.target_id, c.chain_score, c.is_primary))
            .collect::<Vec<_>>();
        kept.sort();
        assert_eq!(kept, vec![(1, 40, 1), (2, 50, 1), (3, 8, 1)]);
    }

    #[test]
    fn test_detect_primary_chains_sets_mapping_quality() {
        let mut chains = vec![
            chain(1, 100, &[(100, 100, 10), (60, 60, 10)]),
            chain(1, 80, &[(95, 120, 10), (65, 90, 10)]),
            chain(1, 70, &[(20, 200, 5), (5, 185, 5)]),
        ];

        detect_primary_chains(&mut chains);

        assert_eq!(chains[0].mapping_quality, 1);
        assert_eq!(chains[1].mapping_quality, 0);
        assert!(chains[2].mapping_quality > 0);
    }

    #[test]
    fn test_compute_score_rejects_and_penalizes_like_cpp() {
        let params = ChainingParameters::new(0.5, 0.1, 1, 0.0);
        let first = seed(10, 0, 10, 5);
        let second = seed(20, 0, 21, 6);
        assert_eq!(compute_score(&second, &first, &params), 5);

        let mut far = params.clone();
        far.band_width = 0;
        assert_eq!(compute_score(&second, &first, &far), i32::MIN);

        let behind = seed(9, 0, 20, 4);
        assert_eq!(compute_score(&behind, &first, &params), i32::MIN);
    }

    #[test]
    fn test_chain_backtrack_builds_chain_from_predecessors() {
        let seeds = vec![seed(10, 7, 10, 5), seed(20, 7, 20, 5), seed(30, 7, 30, 5)];
        let mut data = AnchorData::new(3);
        data.best_score_anchor = vec![5, 10, 15];
        data.predecessor_anchor = vec![-1, 0, 1];

        let chains = chain_backtrack(&mut data, 10, 300, &seeds, true, false);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].target_id, 7);
        assert_eq!(chains[0].chain_score, 15);
        assert!(chains[0].reverse);
        assert_eq!(chains[0].anchors.len(), 3);
    }

    #[test]
    fn test_chaining_dynamic_program_simple_diagonal_chain() {
        let params = ChainingParameters::new(0.5, 0.1, 10, 0.0);
        let seeds = vec![seed(10, 4, 10, 5), seed(20, 4, 20, 5), seed(30, 4, 30, 5)];

        let chains = chaining_dynamic_program(&params, &seeds, false);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].target_id, 4);
        assert_eq!(chains[0].chain_score, 15);
        assert_eq!(chains[0].anchors.len(), 3);
    }
}
