use crate::align::hsp::{Hsp, Match};
use crate::basic::value::{BlockId, Letter};
use crate::dna::alignment::{
    build_hsp_from_cigar, compute_ksw_cigar, compute_wfa_cigar, AlignmentStatus, Cigar, DnaScoring,
    KswMode, WFA_BAND_EXTENSION,
};
use crate::dna::build_score::Blastn_Score;
use crate::dna::chain::{
    chaining_dynamic_program, detect_primary_chains, Anchor, Chain, ChainingParameters,
};
use crate::dna::extension::DnaExtensionAlgo;
use crate::dna::extension_seed_matches::merge_and_extend_seeds;
use crate::dna::seed_set_dna::SeedMatch;
use crate::util::interval::Interval;
use std::sync::Arc;

#[derive(Debug, Clone, Copy)]
pub struct ChainExtensionConfig<'a> {
    pub score_builder: &'a Blastn_Score,
    pub algorithm: DnaExtensionAlgo,
    pub zdrop_extension: i32,
    pub zdrop_global: i32,
    pub band_extension: i32,
    pub band_global: i32,
    pub max_evalue: f64,
    pub best_hsp_only: bool,
}

impl ChainExtensionConfig<'_> {
    fn scoring(self) -> DnaScoring {
        DnaScoring {
            reward: self.score_builder.reward(),
            penalty: self.score_builder.penalty(),
            gap_open: self.score_builder.gap_open(),
            gap_extend: self.score_builder.gap_extend(),
        }
    }

    fn extension_band(self) -> i32 {
        if self.algorithm == DnaExtensionAlgo::Wfa {
            WFA_BAND_EXTENSION
        } else {
            self.band_extension
        }
    }
}

pub fn compute_residue_matches_of_chain(anchors: &[Anchor], _kmer_size: i32) -> i32 {
    let mut total_matching_residues = anchors[anchors.len() - 1].span;

    for i in (1..anchors.len()).rev() {
        total_matching_residues += anchors[i - 1]
            .span
            .min(anchors[i - 1].i - anchors[i].i)
            .min(anchors[i - 1].j - anchors[i].j);
    }

    total_matching_residues
}

pub fn build_map_hsp(
    target_block_id: BlockId,
    target: &[Letter],
    begin_chain: &[Chain],
    kmer_size: i32,
) -> Match {
    let mut m = Match::new(target_block_id, target_block_id as u64);

    for chain in begin_chain {
        let mut map_hsp = Hsp::new();

        map_hsp.query_range.begin = chain.anchors[chain.anchors.len() - 1].i_start();
        map_hsp.subject_range.begin = chain.anchors[chain.anchors.len() - 1].j_start();
        map_hsp.query_range.end = chain.anchors[0].i;
        map_hsp.subject_range.end = chain.anchors[0].j;

        map_hsp.identities = compute_residue_matches_of_chain(&chain.anchors, kmer_size);
        map_hsp.length = (chain.anchors[0].i - chain.anchors[chain.anchors.len() - 1].i_start())
            .max(chain.anchors[0].j - chain.anchors[chain.anchors.len() - 1].j_start());
        map_hsp.mapping_quality = chain.mapping_quality;
        map_hsp.n_anchors = chain.anchors.len() as i32;

        map_hsp.transcript.push_terminator();
        map_hsp.target_seq = Arc::from(target);
        map_hsp.query_source_range = map_hsp.query_range;
        map_hsp.subject_source_range = if chain.reverse {
            Interval::new(map_hsp.subject_range.end, map_hsp.subject_range.begin)
        } else {
            Interval::new(map_hsp.subject_range.begin, map_hsp.subject_range.end)
        };
        map_hsp.frame = chain.reverse as i32 + 2;

        m.hsps.push(map_hsp);
    }

    m
}

pub fn compute_chains(
    mut seed_hits: Vec<SeedMatch>,
    query: &[Letter],
    targets: &[Vec<Letter>],
    kmer_size: i32,
    is_reverse: bool,
    chaining_parameters: &ChainingParameters,
) -> Vec<Chain> {
    if seed_hits.is_empty() {
        return Vec::new();
    }

    seed_hits = merge_and_extend_seeds(&mut seed_hits, query, targets, kmer_size);

    seed_hits.sort_by(|a, b| {
        if a.id() == b.id() {
            a.j().cmp(&b.j())
        } else {
            a.id().cmp(&b.id())
        }
    });

    let mut chains = Vec::new();
    let mut begin = 0usize;
    while begin < seed_hits.len() {
        let target_id = seed_hits[begin].id();
        let mut end = begin + 1;
        while end < seed_hits.len() && seed_hits[end].id() == target_id {
            end += 1;
        }

        let mut new_chains =
            chaining_dynamic_program(chaining_parameters, &seed_hits[begin..end], is_reverse);
        chains.append(&mut new_chains);

        begin = end;
    }

    chains
}

fn align_segment(
    cfg: ChainExtensionConfig<'_>,
    query: &[Letter],
    target: &[Letter],
    extension: &mut Cigar,
    left: bool,
    global: bool,
    band: i32,
) -> Result<AlignmentStatus, String> {
    let query_bytes: Vec<u8> = query.iter().map(|&letter| letter as u8).collect();
    let target_bytes: Vec<u8> = target.iter().map(|&letter| letter as u8).collect();
    if cfg.algorithm == DnaExtensionAlgo::Wfa {
        let query = std::str::from_utf8(&query_bytes)
            .map_err(|_| "WFA query contains a non-UTF-8 nucleotide code".to_owned())?;
        let target = std::str::from_utf8(&target_bytes)
            .map_err(|_| "WFA target contains a non-UTF-8 nucleotide code".to_owned())?;
        compute_wfa_cigar(cfg.scoring(), query, extension, left, global, target, band)
    } else {
        Ok(compute_ksw_cigar(
            &target_bytes,
            &query_bytes,
            cfg.scoring(),
            if left {
                KswMode::Left
            } else if global {
                KswMode::Global
            } else {
                KswMode::Right
            },
            extension,
            if global {
                cfg.zdrop_global
            } else {
                cfg.zdrop_extension
            },
            band,
        ))
    }
}

#[allow(clippy::too_many_arguments)]
pub fn extend_new_at_peak(
    cfg: ChainExtensionConfig<'_>,
    extension: &mut Cigar,
    chain: &Chain,
    query: &[Letter],
    target: &[Letter],
    _target_block_id: BlockId,
    start_i: i32,
    start_j: i32,
    now_anchor_index: i32,
) -> Result<(Hsp, i32), String> {
    let peak_score_cigar_index = extension.peak_score_cigar_index as usize;
    extension.get_cigar_data().truncate(peak_score_cigar_index);
    extension.score = extension.peak_score;

    let peak_anchor = chain
        .anchors
        .get(extension.peak_score_anchor_index as usize)
        .ok_or_else(|| "peak anchor index is outside the chain".to_owned())?;
    let query_begin = usize::try_from(peak_anchor.i)
        .map_err(|_| "negative query anchor coordinate".to_owned())?;
    let target_begin = usize::try_from(peak_anchor.j)
        .map_err(|_| "negative target anchor coordinate".to_owned())?;
    let query_right = query
        .get(query_begin..)
        .ok_or_else(|| "query anchor coordinate is outside the sequence".to_owned())?;
    let target_end = target.len().min(
        target_begin
            .saturating_add(cfg.band_extension.max(0) as usize)
            .saturating_add(query_right.len())
            .min(target_begin.saturating_add(query_right.len().saturating_mul(2))),
    );
    let target_right = target
        .get(target_begin..target_end)
        .ok_or_else(|| "target anchor coordinate is outside the sequence".to_owned())?;
    align_segment(
        cfg,
        query_right,
        target_right,
        extension,
        false,
        false,
        cfg.extension_band(),
    )?;

    Ok((
        build_hsp_from_cigar(
            extension,
            target,
            query,
            start_i,
            start_j,
            chain.reverse,
            cfg.score_builder,
            cfg.max_evalue,
        ),
        now_anchor_index,
    ))
}

pub fn extend_between_anchors(
    cfg: ChainExtensionConfig<'_>,
    target_block_id: BlockId,
    chain: &Chain,
    query: &[Letter],
    target: &[Letter],
    mut anchor_idx: i32,
) -> Result<(Hsp, i32), String> {
    let anchor = chain
        .anchors
        .get(anchor_idx as usize)
        .ok_or_else(|| "anchor index is outside the chain".to_owned())?;
    let start_i = anchor.i_start();
    let start_j = anchor.j_start();
    let mut extension = Cigar::with_reserve(anchor_idx.max(0) as usize * 3);

    if start_i > 0 && start_j > 0 {
        let (mut prev_i, mut prev_j) = (0, 0);
        if (anchor_idx as usize) + 1 < chain.anchors.len() {
            let previous = chain.anchors[anchor_idx as usize + 1];
            if previous.i <= start_i && previous.j <= start_j {
                prev_i = previous.i;
                prev_j = previous.j;
            }
        }
        let mut query_left = query[prev_i as usize..start_i as usize].to_vec();
        query_left.reverse();
        let target_left_begin = prev_j.max(
            (start_j - query_left.len() as i32 - cfg.band_extension)
                .max(start_j - query_left.len() as i32 * 2),
        );
        let mut target_left = target[target_left_begin as usize..start_j as usize].to_vec();
        target_left.reverse();
        align_segment(
            cfg,
            &query_left,
            &target_left,
            &mut extension,
            true,
            false,
            cfg.extension_band(),
        )?;
    } else {
        extension.set_max_values(-1, -1);
    }

    let mut anchor_distance_query = i32::MAX;
    let mut anchor_distance_target = i32::MAX;
    while anchor_idx > 0 {
        let current = chain.anchors[anchor_idx as usize];
        let next = chain.anchors[anchor_idx as usize - 1];
        if anchor_distance_query > current.span && anchor_distance_target > current.span {
            extension.extend_cigar_op(current.span as u32, 'M');
            extension.score += current.span * cfg.score_builder.reward();
        }
        if extension.score > extension.peak_score {
            extension.peak_score = extension.score;
            extension.peak_score_cigar_index = extension.get_cigar_data_const().len() as i32;
            extension.peak_score_anchor_index = anchor_idx;
        }

        anchor_distance_query = next.i - current.i;
        anchor_distance_target = next.j - current.j;
        if anchor_distance_query > next.span && anchor_distance_target > next.span {
            let query_right = &query[current.i as usize..next.i_start() as usize];
            let target_right = &target[current.j as usize..next.j_start() as usize];
            let alignment_band = (query_right.len() as i32 - target_right.len() as i32).abs()
                + cfg
                    .band_global
                    .min(query_right.len().min(target_right.len()) as i32 / 2);
            let status = align_segment(
                cfg,
                query_right,
                target_right,
                &mut extension,
                false,
                true,
                alignment_band,
            )?;
            if matches!(
                status,
                AlignmentStatus::Dropped | AlignmentStatus::NegativeScore
            ) {
                if extension.score >= extension.peak_score {
                    return Ok((
                        build_hsp_from_cigar(
                            &extension,
                            target,
                            query,
                            start_i,
                            start_j,
                            chain.reverse,
                            cfg.score_builder,
                            cfg.max_evalue,
                        ),
                        anchor_idx - 1,
                    ));
                }
                return extend_new_at_peak(
                    cfg,
                    &mut extension,
                    chain,
                    query,
                    target,
                    target_block_id,
                    start_i,
                    start_j,
                    anchor_idx - 1,
                );
            }
        } else if anchor_distance_query <= next.span
            && anchor_distance_query < anchor_distance_target
        {
            let gaps = anchor_distance_target - anchor_distance_query;
            extension.extend_cigar_op(gaps as u32, 'D');
            extension.score -= gaps * cfg.score_builder.gap_extend() + cfg.score_builder.gap_open();
            extension.extend_cigar_op(anchor_distance_query as u32, 'M');
            extension.score += anchor_distance_query * cfg.score_builder.reward();
        } else if anchor_distance_target <= next.span
            && anchor_distance_target < anchor_distance_query
        {
            let gaps = anchor_distance_query - anchor_distance_target;
            extension.extend_cigar_op(gaps as u32, 'I');
            extension.score -= gaps * cfg.score_builder.gap_extend() + cfg.score_builder.gap_open();
            extension.extend_cigar_op(anchor_distance_target as u32, 'M');
            extension.score += anchor_distance_target * cfg.score_builder.reward();
        } else {
            return Err("Error in chaining extension: No anchor overlap case matched".to_owned());
        }
        anchor_idx -= 1;
    }

    let last = chain.anchors[0];
    if anchor_distance_query > last.span && anchor_distance_target > last.span {
        extension.extend_cigar_op(last.span as u32, 'M');
        extension.score += last.span * cfg.score_builder.reward();
    }
    let query_right = &query[last.i as usize..];
    let target_begin = last.j as usize;
    let target_end = target.len().min(
        target_begin
            .saturating_add(cfg.band_extension.max(0) as usize)
            .saturating_add(query_right.len())
            .min(target_begin.saturating_add(query_right.len().saturating_mul(2))),
    );
    align_segment(
        cfg,
        query_right,
        &target[target_begin..target_end],
        &mut extension,
        false,
        false,
        cfg.extension_band(),
    )?;
    if extension.score >= extension.peak_score {
        Ok((
            build_hsp_from_cigar(
                &extension,
                target,
                query,
                start_i,
                start_j,
                chain.reverse,
                cfg.score_builder,
                cfg.max_evalue,
            ),
            anchor_idx - 1,
        ))
    } else {
        extend_new_at_peak(
            cfg,
            &mut extension,
            chain,
            query,
            target,
            target_block_id,
            start_i,
            start_j,
            anchor_idx - 1,
        )
    }
}

pub fn extend_chains(
    cfg: ChainExtensionConfig<'_>,
    target_block_id: BlockId,
    chains: &[Chain],
    query: &[Letter],
    query_reverse: &[Letter],
    target: &[Letter],
    chaining_parameters: &ChainingParameters,
) -> Result<Match, String> {
    let mut m = Match::new(target_block_id, target_block_id as u64);
    let mut extended_ranges: Vec<(i32, i32, i32, i32)> = Vec::new();
    for chain in chains {
        let first = chain.anchors[0];
        let last = chain.anchors[chain.anchors.len() - 1];
        let query_range = first.i - last.i_start();
        let target_range = first.j - last.j_start();
        if extended_ranges.iter().any(|&(qb, qe, tb, te)| {
            let query_overlap = first.i.min(qe) - last.i_start().max(qb);
            let target_overlap = first.j.min(te) - last.j_start().max(tb);
            query_overlap > (chaining_parameters.max_overlap_extension * query_range as f32) as i32
                && target_overlap
                    > (chaining_parameters.max_overlap_extension * target_range as f32) as i32
        }) {
            continue;
        }

        let selected_query = if chain.reverse { query_reverse } else { query };
        let mut anchor_idx = chain.anchors.len() as i32 - 1;
        loop {
            let (hsp, current_index) = extend_between_anchors(
                cfg,
                target_block_id,
                chain,
                selected_query,
                target,
                anchor_idx,
            )?;
            if hsp.evalue < cfg.max_evalue {
                extended_ranges.push((
                    hsp.query_range.begin,
                    hsp.query_range.end,
                    hsp.subject_range.begin,
                    hsp.subject_range.end,
                ));
                m.hsps.push(hsp);
            }
            anchor_idx = current_index;
            if anchor_idx < 0 {
                break;
            }
        }
    }
    if cfg.best_hsp_only {
        if let Some(best) = m.hsps.iter().max_by_key(|hsp| hsp.score).cloned() {
            m.hsps.clear();
            m.hsps.push(best);
        }
    }
    Ok(m)
}

pub fn chaining_and_extension(
    mut chains: Vec<Chain>,
    targets: &[Vec<Letter>],
    chain_fraction_align: f32,
    chaining_out: bool,
    kmer_size: i32,
    query: &[Letter],
    query_reverse: &[Letter],
    chaining_parameters: &ChainingParameters,
    extension_config: ChainExtensionConfig<'_>,
) -> Result<Vec<Match>, String> {
    let mut matches = Vec::new();

    if chains.is_empty() {
        return Ok(matches);
    }

    chains.sort_by(|a, b| b.cmp(a));

    if chaining_out {
        detect_primary_chains(&mut chains);
    }

    let map_score_threshold = (chains[0].chain_score as f32 * chain_fraction_align) as i32;
    let lower_bound = chains
        .iter()
        .position(|chain| chain.chain_score < map_score_threshold)
        .unwrap_or(chains.len());
    chains.truncate(lower_bound);

    chains.sort_by(|a, b| {
        if a.target_id == b.target_id {
            b.chain_score.cmp(&a.chain_score)
        } else {
            a.target_id.cmp(&b.target_id)
        }
    });

    let mut begin = 0usize;
    while begin < chains.len() {
        let target_id = chains[begin].target_id;
        let mut end = begin + 1;
        while end < chains.len() && chains[end].target_id == target_id {
            end += 1;
        }

        let m = if chaining_out {
            build_map_hsp(
                target_id,
                &targets[target_id as usize],
                &chains[begin..end],
                kmer_size,
            )
        } else {
            extend_chains(
                extension_config,
                target_id,
                &chains[begin..end],
                query,
                query_reverse,
                &targets[target_id as usize],
                chaining_parameters,
            )?
        };

        if !m.hsps.is_empty() {
            matches.push(m);
        }

        begin = end;
    }

    Ok(matches)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn score_builder() -> Blastn_Score {
        Blastn_Score::new_with_parameters(2, -3, 5, 1, 1_000, 2, 0.625, 0.41, 1.1)
    }

    fn extension_config(score_builder: &Blastn_Score) -> ChainExtensionConfig<'_> {
        ChainExtensionConfig {
            score_builder,
            algorithm: DnaExtensionAlgo::Ksw,
            zdrop_extension: 40,
            zdrop_global: 100,
            band_extension: 40,
            band_global: 30,
            max_evalue: f64::INFINITY,
            best_hsp_only: false,
        }
    }

    fn dna(s: &[u8]) -> Vec<Letter> {
        s.iter().map(|&x| x as Letter).collect()
    }

    fn chain(reverse: bool, anchors: &[(i32, i32, i32)]) -> Chain {
        let mut c = Chain::new(reverse);
        c.mapping_quality = 17;
        c.anchors = anchors
            .iter()
            .map(|&(i, j, span)| Anchor::new(i, j, span))
            .collect();
        c
    }

    #[test]
    fn test_compute_residue_matches_of_chain() {
        let anchors = vec![
            Anchor::new(40, 42, 10),
            Anchor::new(25, 28, 8),
            Anchor::new(10, 12, 7),
        ];

        assert_eq!(compute_residue_matches_of_chain(&anchors, 10), 25);
    }

    #[test]
    fn test_build_map_hsp_forward_and_reverse() {
        let target = dna(b"AACCGGTT");
        let chains = vec![
            chain(false, &[(40, 42, 10), (25, 28, 8), (10, 12, 7)]),
            chain(true, &[(30, 35, 6), (20, 23, 5)]),
        ];

        let m = build_map_hsp(3, &target, &chains, 10);

        assert_eq!(m.target_block_id, 3);
        assert_eq!(m.hsps.len(), 2);
        assert_eq!(m.hsps[0].query_range, Interval::new(3, 40));
        assert_eq!(m.hsps[0].subject_range, Interval::new(5, 42));
        assert_eq!(m.hsps[0].identities, 25);
        assert_eq!(m.hsps[0].length, 37);
        assert_eq!(m.hsps[0].mapping_quality, 17);
        assert_eq!(m.hsps[0].n_anchors, 3);
        assert_eq!(m.hsps[0].subject_source_range, Interval::new(5, 42));
        assert_eq!(m.hsps[0].frame, 2);
        assert!(m.hsps[0].transcript.data().last().unwrap().is_terminator());

        assert_eq!(m.hsps[1].subject_source_range, Interval::new(35, 18));
        assert_eq!(m.hsps[1].frame, 3);
    }

    fn seed(i: i32, id: BlockId, j: i32, score: i32) -> SeedMatch {
        let mut s = SeedMatch::new(i, id, j, score);
        s.score(score);
        s
    }

    #[test]
    fn test_compute_chains_groups_seed_hits_by_target() {
        let query = dna(b"AACCGGTTAACCGGTT");
        let targets = vec![dna(b"AACCGGTTAACCGGTT"), dna(b"TTCCAACCGGTTAACC")];
        let seed_hits = vec![
            seed(2, 1, 6, 3),
            seed(8, 0, 8, 3),
            seed(2, 0, 2, 3),
            seed(12, 0, 12, 3),
        ];
        let params = ChainingParameters::new(0.5, 0.1, 6, 0.5);

        let mut chains = compute_chains(seed_hits, &query, &targets, 3, false, &params);
        chains.sort_by_key(|c| c.target_id);

        assert_eq!(chains.len(), 2);
        assert_eq!(chains[0].target_id, 0);
        assert!(chains[0].chain_score >= 6);
        assert_eq!(chains[1].target_id, 1);
    }

    #[test]
    fn test_chaining_and_extension_chaining_out_builds_grouped_map_matches() {
        let targets = vec![dna(b"AACCGGTT"), dna(b"TTCCAACC")];
        let mut c0 = chain(false, &[(40, 42, 10), (25, 28, 8), (10, 12, 7)]);
        c0.target_id = 0;
        c0.chain_score = 100;
        let mut c1 = chain(true, &[(30, 35, 6), (20, 23, 5)]);
        c1.target_id = 1;
        c1.chain_score = 80;
        let mut filtered = chain(false, &[(10, 10, 4)]);
        filtered.target_id = 0;
        filtered.chain_score = 10;

        let scoring = score_builder();
        let params = ChainingParameters::new(0.5, 0.1, 6, 0.5);
        let matches = chaining_and_extension(
            vec![filtered, c1, c0],
            &targets,
            0.0,
            true,
            10,
            &[],
            &[],
            &params,
            extension_config(&scoring),
        )
        .unwrap();

        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].target_block_id, 0);
        assert_eq!(matches[0].hsps.len(), 2);
        assert_eq!(matches[0].hsps[0].n_anchors, 3);
        assert_eq!(matches[1].target_block_id, 1);
        assert_eq!(matches[1].hsps[0].frame, 3);
    }

    #[test]
    fn test_extend_between_anchors_aligns_single_anchor_to_both_ends() {
        let query = dna(b"ACGTACGT");
        let target = query.clone();
        let mut c = chain(false, &[(4, 4, 4)]);
        c.target_id = 0;
        let scoring = score_builder();

        let (hsp, next) =
            extend_between_anchors(extension_config(&scoring), 0, &c, &query, &target, 0).unwrap();

        assert_eq!(next, -1);
        assert_eq!(hsp.query_range, Interval::new(0, 8));
        assert_eq!(hsp.subject_range, Interval::new(0, 8));
        assert_eq!(hsp.score, 16);
        assert_eq!(hsp.identities, 8);
        assert!(hsp.transcript.data().last().unwrap().is_terminator());
    }

    #[test]
    fn test_chaining_and_extension_non_mapping_mode_returns_alignment() {
        let query = dna(b"ACGTACGT");
        let targets = vec![query.clone()];
        let mut c = chain(false, &[(4, 4, 4)]);
        c.target_id = 0;
        c.chain_score = 20;
        let scoring = score_builder();
        let params = ChainingParameters::new(0.5, 0.1, 6, 0.5);

        let matches = chaining_and_extension(
            vec![c],
            &targets,
            0.0,
            false,
            4,
            &query,
            &query,
            &params,
            extension_config(&scoring),
        )
        .unwrap();

        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].hsps.len(), 1);
        assert_eq!(matches[0].hsps[0].score, 16);
    }

    #[test]
    fn test_extend_between_anchors_globally_aligns_inter_anchor_gap() {
        let query = dna(b"AACCGGTTAA");
        let target = dna(b"AACCGGGTTAA");
        let mut c = chain(false, &[(8, 9, 2), (4, 4, 2)]);
        c.target_id = 0;
        let scoring = score_builder();

        let (hsp, next) =
            extend_between_anchors(extension_config(&scoring), 0, &c, &query, &target, 1).unwrap();

        assert_eq!(next, -1);
        assert_eq!(hsp.query_range, Interval::new(0, 10));
        assert_eq!(hsp.subject_range, Interval::new(0, 11));
        assert_eq!(hsp.gaps, 1);
        assert_eq!(hsp.gap_openings, 1);
    }

    #[test]
    fn test_extend_chains_suppresses_overlapping_later_chain() {
        let query = dna(b"ACGTACGT");
        let target = query.clone();
        let mut first = chain(false, &[(4, 4, 4)]);
        first.target_id = 0;
        let second = first.clone();
        let scoring = score_builder();
        let params = ChainingParameters::new(0.5, 0.1, 6, 0.5);

        let m = extend_chains(
            extension_config(&scoring),
            0,
            &[first, second],
            &query,
            &query,
            &target,
            &params,
        )
        .unwrap();

        assert_eq!(m.hsps.len(), 1);
    }
}
