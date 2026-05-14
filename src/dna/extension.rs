use crate::align::hsp::Hsp;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::Letter;
use crate::dna::seed_set_dna::SeedMatch;
use crate::util::interval::Interval;

pub const KSW2_END_BONUS: i32 = 5;
pub const KSW2_BAND: i32 = 40;
pub const WFA_CUTOFF_STEPS: i32 = 10;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DnaExtensionAlgo {
    Ksw,
    Wfa,
}

impl DnaExtensionAlgo {
    pub fn to_string(self) -> &'static str {
        match self {
            Self::Ksw => "ksw",
            Self::Wfa => "wfa",
        }
    }

    pub fn from_string(s: &str) -> Option<Self> {
        match s {
            "ksw" => Some(Self::Ksw),
            "wfa" => Some(Self::Wfa),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExtendedSeed {
    pub i_min_extended: i32,
    pub i_max_extended: i32,
    pub j_min_extended: i32,
    pub j_max_extended: i32,
    pub length: i32,
}

impl ExtendedSeed {
    pub fn new(i_min: i32, i_max: i32, j_min: i32, j_max: i32) -> Self {
        Self {
            i_min_extended: i_min,
            i_max_extended: i_max,
            j_min_extended: j_min,
            j_max_extended: j_max,
            length: i_max - i_min,
        }
    }
}

pub fn intersection(hit: &SeedMatch, extended: &[ExtendedSeed]) -> bool {
    if extended.is_empty() {
        return false;
    }

    extended.iter().any(|s| {
        hit.i_start() >= s.i_min_extended
            && hit.i() <= s.i_max_extended
            && hit.j_start() >= s.j_min_extended
            && hit.j() <= s.j_max_extended
    })
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CigarShort {
    pub cigar_data: Vec<(i32, char)>,
    score: i32,
    max_query: i32,
    max_target: i32,
}

impl CigarShort {
    pub fn score(&self) -> i32 {
        self.score
    }

    pub fn max_query(&self) -> i32 {
        self.max_query
    }

    pub fn max_target(&self) -> i32 {
        self.max_target
    }
}

impl std::ops::Add for CigarShort {
    type Output = CigarShort;

    fn add(mut self, other: CigarShort) -> CigarShort {
        self.cigar_data.extend(other.cigar_data);
        self.score += other.score;
        self
    }
}

pub fn cigar_to_hsp_seed_match(
    target: &[Letter],
    query: &[Letter],
    hit: &SeedMatch,
    out: &mut Hsp,
    reverse: bool,
) {
    let mut pattern_pos = hit.i_start();
    let mut text_pos = hit.j_start();
    out.query_range.begin = pattern_pos;
    out.subject_range.begin = text_pos;

    for _ in 0..hit.ungapped_score() {
        out.push_match(target[text_pos as usize], query[pattern_pos as usize], true);
        pattern_pos += 1;
        text_pos += 1;
    }

    out.query_range.end = pattern_pos;
    out.subject_range.end = text_pos;
    out.transcript.push_terminator();
    out.target_seq = target.to_vec();
    out.query_source_range = out.query_range;
    out.subject_source_range = if reverse {
        Interval::new(out.subject_range.end, out.subject_range.begin)
    } else {
        Interval::new(out.subject_range.begin, out.subject_range.end)
    };
    out.frame = reverse as i32;
}

pub fn cigar_to_hsp(
    cigar: &CigarShort,
    target: &[Letter],
    query: &[Letter],
    pos_i: i32,
    pos_j: i32,
    out: &mut Hsp,
    reverse: bool,
) {
    let mut pattern_pos = pos_i - cigar.max_query() - 1;
    let mut text_pos = pos_j - cigar.max_target() - 1;
    out.query_range.begin = pattern_pos;
    out.subject_range.begin = text_pos;

    for operation in &cigar.cigar_data {
        match operation.1 {
            'M' | '=' | 'X' => {
                for _ in 0..operation.0 {
                    out.push_match(target[text_pos as usize], query[pattern_pos as usize], true);
                    pattern_pos += 1;
                    text_pos += 1;
                }
            }
            'D' => {
                let end = (text_pos + operation.0) as usize;
                out.push_gap(EditOperation::Deletion, operation.0, &target[..end]);
                text_pos += operation.0;
            }
            'I' => {
                let end = (pattern_pos + operation.0).max(0) as usize;
                out.push_gap(EditOperation::Insertion, operation.0, &query[..end]);
                pattern_pos += operation.0;
            }
            _ => {}
        }
    }

    out.query_range.end = pattern_pos;
    out.subject_range.end = text_pos;
    out.transcript.push_terminator();
    out.target_seq = target.to_vec();
    out.query_source_range = out.query_range;
    out.subject_source_range = if reverse {
        Interval::new(out.subject_range.end, out.subject_range.begin)
    } else {
        Interval::new(out.subject_range.begin, out.subject_range.end)
    };
    out.frame = reverse as i32 + 2;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::EditOperation;
    use crate::basic::value::BlockId;

    fn dna(s: &[u8]) -> Vec<Letter> {
        s.iter().map(|&x| x as Letter).collect()
    }

    fn seed(i: i32, id: BlockId, j: i32, score: i32) -> SeedMatch {
        let mut s = SeedMatch::new(i, id, j, score);
        s.score(score);
        s
    }

    #[test]
    fn test_dna_extension_algo_parse_and_display() {
        assert_eq!(
            DnaExtensionAlgo::from_string("ksw"),
            Some(DnaExtensionAlgo::Ksw)
        );
        assert_eq!(
            DnaExtensionAlgo::from_string("wfa"),
            Some(DnaExtensionAlgo::Wfa)
        );
        assert_eq!(DnaExtensionAlgo::from_string("bad"), None);
        assert_eq!(DnaExtensionAlgo::Wfa.to_string(), "wfa");
    }

    #[test]
    fn test_extended_seed_and_intersection() {
        let ext = vec![ExtendedSeed::new(5, 20, 10, 25)];
        assert_eq!(ext[0].length, 15);
        assert!(intersection(&seed(15, 0, 20, 5), &ext));
        assert!(!intersection(&seed(25, 0, 30, 5), &ext));
        assert!(!intersection(&seed(15, 0, 20, 5), &[]));
    }

    #[test]
    fn test_cigar_short_add_and_accessors() {
        let mut left = CigarShort::default();
        left.cigar_data.push((2, '='));
        left.score = 4;
        left.max_query = 1;
        left.max_target = 2;

        let mut right = CigarShort::default();
        right.cigar_data.push((1, 'I'));
        right.score = -3;

        let combined = left + right;
        assert_eq!(combined.cigar_data, vec![(2, '='), (1, 'I')]);
        assert_eq!(combined.score(), 1);
        assert_eq!(combined.max_query(), 1);
        assert_eq!(combined.max_target(), 2);
    }

    #[test]
    fn test_cigar_to_hsp_seed_match() {
        let target = dna(b"AACCGG");
        let query = dna(b"AACAGG");
        let hit = seed(4, 0, 4, 4);
        let mut hsp = Hsp::new();

        cigar_to_hsp_seed_match(&target, &query, &hit, &mut hsp, false);

        assert_eq!(hsp.query_range, Interval::new(0, 4));
        assert_eq!(hsp.subject_range, Interval::new(0, 4));
        assert_eq!(hsp.subject_source_range, Interval::new(0, 4));
        assert_eq!(hsp.frame, 0);
        assert_eq!(hsp.length, 4);
        assert_eq!(hsp.identities, 3);
        assert_eq!(hsp.mismatches, 1);
    }

    #[test]
    fn test_cigar_to_hsp_short() {
        let target = dna(b"AACCGGTT");
        let query = dna(b"AACAGGTT");
        let mut cigar = CigarShort::default();
        cigar.cigar_data = vec![(3, '='), (1, 'I'), (2, 'M'), (1, 'D')];
        let mut hsp = Hsp::new();

        cigar_to_hsp(&cigar, &target, &query, 1, 1, &mut hsp, true);

        assert_eq!(hsp.query_range, Interval::new(0, 6));
        assert_eq!(hsp.subject_range, Interval::new(0, 6));
        assert_eq!(hsp.subject_source_range, Interval::new(6, 0));
        assert_eq!(hsp.frame, 3);
        assert_eq!(hsp.length, 7);
        assert_eq!(hsp.gaps, 2);
        let ops = hsp.transcript.iter().collect::<Vec<_>>();
        assert_eq!(ops[0].op, EditOperation::Match);
        assert!(hsp.transcript.data().last().unwrap().is_terminator());
    }
}
