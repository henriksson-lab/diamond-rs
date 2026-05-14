use crate::basic::value::Letter;
use crate::dna::seed_set_dna::SeedMatch;

pub fn extend_seeds_ungapped(
    seed_match_begin: &[SeedMatch],
    query: &[Letter],
    target: &[Letter],
    kmer_size: i32,
) -> Vec<SeedMatch> {
    let mut new_hits: Vec<SeedMatch> = Vec::with_capacity(seed_match_begin.len());
    let mut prev_diagonal = i32::MIN;

    let query_length = query.len() as i32;
    let target_length = target.len() as i32;
    for hit in seed_match_begin {
        let diagonal = hit.i() - hit.j();

        if diagonal == prev_diagonal && hit.i() <= new_hits[new_hits.len() - 1].i() {
            continue;
        }

        let mut left_i = hit.i() - 1;
        let mut left_j = hit.j() - 1;
        while left_i > -1 && left_j > -1 && query[left_i as usize] == target[left_j as usize] {
            left_i -= 1;
            left_j -= 1;
        }

        let mut right_i = hit.i() + kmer_size;
        let mut right_j = hit.j() + kmer_size;
        while right_i < query_length
            && right_j < target_length
            && query[right_i as usize] == target[right_j as usize]
        {
            right_i += 1;
            right_j += 1;
        }

        let mut new_hit = SeedMatch::new(right_i, hit.id(), right_j, kmer_size);
        new_hit.score(right_i - left_i - 1);
        new_hits.push(new_hit);

        prev_diagonal = diagonal;
    }

    new_hits
}

pub fn compare_by_target_and_diagonal(a: &SeedMatch, b: &SeedMatch) -> bool {
    if a.id() != b.id() {
        return a.id() < b.id();
    }
    if (a.i() - a.j()) != (b.i() - b.j()) {
        return (a.i() - a.j()) < (b.i() - b.j());
    }
    a.i() < b.i()
}

pub fn merge_and_extend_seeds(
    seed_hits: &mut [SeedMatch],
    query: &[Letter],
    targets: &[Vec<Letter>],
    kmer_size: i32,
) -> Vec<SeedMatch> {
    seed_hits.sort_by(|a, b| {
        if compare_by_target_and_diagonal(a, b) {
            std::cmp::Ordering::Less
        } else if compare_by_target_and_diagonal(b, a) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });

    let mut new_seed_hits = Vec::with_capacity(seed_hits.len());
    let mut begin = 0usize;
    while begin < seed_hits.len() {
        let target_id = seed_hits[begin].id();
        let mut end = begin + 1;
        while end < seed_hits.len() && seed_hits[end].id() == target_id {
            end += 1;
        }

        let new_hits = extend_seeds_ungapped(
            &seed_hits[begin..end],
            query,
            &targets[target_id as usize],
            kmer_size,
        );
        new_seed_hits.extend(new_hits);
        begin = end;
    }

    new_seed_hits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::BlockId;

    fn seed(i: i32, id: BlockId, j: i32, score: i32) -> SeedMatch {
        let mut s = SeedMatch::new(i, id, j, score);
        s.score(score);
        s
    }

    #[test]
    fn test_extend_seeds_ungapped_extends_left_and_right() {
        let query = b"AACCGGTT".iter().map(|&x| x as Letter).collect::<Vec<_>>();
        let target = b"TACCGGTA".iter().map(|&x| x as Letter).collect::<Vec<_>>();
        let hits = vec![seed(2, 0, 2, 3)];

        let extended = extend_seeds_ungapped(&hits, &query, &target, 3);

        assert_eq!(extended.len(), 1);
        assert_eq!(extended[0].i(), 7);
        assert_eq!(extended[0].j(), 7);
        assert_eq!(extended[0].ungapped_score(), 6);
    }

    #[test]
    fn test_extend_seeds_ungapped_skips_contained_same_diagonal_hits() {
        let query = b"AACCGGTT".iter().map(|&x| x as Letter).collect::<Vec<_>>();
        let target = b"AACCGGTT".iter().map(|&x| x as Letter).collect::<Vec<_>>();
        let hits = vec![seed(1, 0, 1, 2), seed(3, 0, 3, 2)];

        let extended = extend_seeds_ungapped(&hits, &query, &target, 2);

        assert_eq!(extended.len(), 1);
        assert_eq!(extended[0].i(), 8);
        assert_eq!(extended[0].j(), 8);
        assert_eq!(extended[0].ungapped_score(), 8);
    }

    #[test]
    fn test_compare_by_target_and_diagonal() {
        assert!(compare_by_target_and_diagonal(
            &seed(10, 1, 9, 2),
            &seed(10, 2, 5, 2)
        ));
        assert!(compare_by_target_and_diagonal(
            &seed(8, 1, 10, 2),
            &seed(10, 1, 9, 2)
        ));
        assert!(compare_by_target_and_diagonal(
            &seed(5, 1, 4, 2),
            &seed(8, 1, 7, 2)
        ));
    }

    #[test]
    fn test_merge_and_extend_seeds_groups_by_target_after_sorting() {
        let query = b"AACCGGTT".iter().map(|&x| x as Letter).collect::<Vec<_>>();
        let targets = vec![
            b"AACCGGTT".iter().map(|&x| x as Letter).collect::<Vec<_>>(),
            b"TACCGGTA".iter().map(|&x| x as Letter).collect::<Vec<_>>(),
        ];
        let mut hits = vec![seed(2, 1, 2, 3), seed(3, 0, 3, 2), seed(1, 0, 1, 2)];

        let extended = merge_and_extend_seeds(&mut hits, &query, &targets, 2);

        assert_eq!(extended.len(), 2);
        assert_eq!(extended[0].id(), 0);
        assert_eq!(extended[0].ungapped_score(), 8);
        assert_eq!(extended[1].id(), 1);
        assert_eq!(extended[1].ungapped_score(), 6);
    }
}
