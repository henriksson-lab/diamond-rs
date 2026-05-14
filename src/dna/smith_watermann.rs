use crate::align::hsp::{Hsp, Match};
use crate::basic::value::{BlockId, Letter};

pub fn scoring_function(match_score: i32, mismatch: i32, first: Letter, second: Letter) -> i32 {
    if first == second {
        match_score
    } else {
        mismatch
    }
}

pub fn dynamic_programm(
    target: &[Letter],
    query: &[Letter],
    match_score: i32,
    mismatch: i32,
    gapopen: i32,
    gapextend: i32,
) -> i32 {
    let mut dp_matrix = vec![vec![0; target.len() + 1]; query.len() + 1];
    let mut hgap = vec![0; target.len() + 1];

    let mut max = 0;
    for i in 1..dp_matrix.len() {
        let mut vgap = 0;
        for j in 1..dp_matrix[0].len() {
            let mut s = dp_matrix[i - 1][j - 1]
                + scoring_function(match_score, mismatch, query[i - 1], target[j - 1]);
            s = s.max(hgap[j]).max(vgap).max(0);
            dp_matrix[i][j] = s;
            vgap -= gapextend;
            hgap[j] -= gapextend;
            let open = s - gapopen - gapextend;
            vgap = vgap.max(open);
            hgap[j] = hgap[j].max(open);

            if dp_matrix[i][j] > max {
                max = dp_matrix[i][j];
            }
        }
    }

    max
}

pub fn local_alignment(
    target_seqs: &[Vec<Letter>],
    query_sequence: &[Letter],
    query_id: BlockId,
    match_score: i32,
    mismatch: i32,
    gapopen: i32,
    gapextend: i32,
) -> Vec<Match> {
    let mut matches = Vec::new();

    for (i, target_sequence) in target_seqs.iter().enumerate() {
        let mut m = Match::new(query_id, i as u64);
        let score = dynamic_programm(
            target_sequence,
            query_sequence,
            match_score,
            mismatch,
            gapopen,
            gapextend,
        );
        m.hsps.push(Hsp::with_score(false, score, 0));
        matches.push(m);
    }

    matches
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(s: &[u8]) -> Vec<Letter> {
        s.iter().map(|&x| x as Letter).collect()
    }

    #[test]
    fn test_scoring_function() {
        assert_eq!(scoring_function(2, -3, b'A' as Letter, b'A' as Letter), 2);
        assert_eq!(scoring_function(2, -3, b'A' as Letter, b'C' as Letter), -3);
    }

    #[test]
    fn test_dynamic_programm_identical_and_mismatch() {
        let query = dna(b"ACGT");
        let target = dna(b"ACGT");
        assert_eq!(dynamic_programm(&target, &query, 2, -3, 5, 2), 8);

        let target_mismatch = dna(b"AGGT");
        assert_eq!(dynamic_programm(&target_mismatch, &query, 2, -3, 5, 2), 4);
    }

    #[test]
    fn test_dynamic_programm_gap_extension_recurrence() {
        let query = dna(b"AAAA");
        let target = dna(b"AAATAA");
        assert_eq!(dynamic_programm(&target, &query, 2, -3, 5, 1), 6);
    }

    #[test]
    fn test_local_alignment_builds_matches() {
        let query = dna(b"ACGT");
        let targets = vec![dna(b"ACGT"), dna(b"TTTT")];
        let matches = local_alignment(&targets, &query, 9, 2, -3, 5, 2);

        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].target_block_id, 9);
        assert_eq!(matches[0].target_oid, 0);
        assert_eq!(matches[0].hsps[0].score, 8);
        assert_eq!(matches[1].target_oid, 1);
        assert_eq!(matches[1].hsps[0].score, 2);
    }
}
