use crate::align::hsp::{Hsp, Match};
use crate::basic::value::{BlockId, Letter};
use crate::dna::chain::{
    chaining_dynamic_program, detect_primary_chains, Anchor, Chain, ChainingParameters,
};
use crate::dna::extension_seed_matches::merge_and_extend_seeds;
use crate::dna::seed_set_dna::SeedMatch;
use crate::util::interval::Interval;

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
        map_hsp.target_seq = target.to_vec();
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

pub fn chaining_and_extension(
    mut chains: Vec<Chain>,
    targets: &[Vec<Letter>],
    chain_fraction_align: i32,
    chaining_out: bool,
    kmer_size: i32,
) -> Vec<Match> {
    let mut matches = Vec::new();

    if chains.is_empty() {
        return matches;
    }

    chains.sort_by(|a, b| b.cmp(a));

    if chaining_out {
        detect_primary_chains(&mut chains);
    }

    let map_score_threshold = chains[0].chain_score * chain_fraction_align;
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
            Match::new(target_id, target_id as u64)
        };

        if !m.hsps.is_empty() {
            matches.push(m);
        }

        begin = end;
    }

    matches
}

#[cfg(test)]
mod tests {
    use super::*;

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

        let matches = chaining_and_extension(vec![filtered, c1, c0], &targets, 0, true, 10);

        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].target_block_id, 0);
        assert_eq!(matches[0].hsps.len(), 2);
        assert_eq!(matches[0].hsps[0].n_anchors, 3);
        assert_eq!(matches[1].target_block_id, 1);
        assert_eq!(matches[1].hsps[0].frame, 3);
    }
}
