pub const KSW2_ZDROP_EXTENSION: i32 = 40;
pub const KSW2_ZDROP_BETWEEN_ANCHORS: i32 = 100;
pub const KSW2_BAND_EXTENSION: i32 = 40;
pub const KSW2_BAND_GLOBAL: i32 = 30;
pub const WFA_BAND_EXTENSION: i32 = 20;
pub const WFA_ZDROP_EXTENSION: i32 = 100;
pub const WFA_ZDROP_GLOBAL: i32 = 500;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cigar {
    query_extension_distance: i32,
    target_extension_distance: i32,
    cigar_data: Vec<(i32, char)>,
    pub score: i32,
    pub peak_score: i32,
    pub peak_score_cigar_index: i32,
    pub peak_score_anchor_index: i32,
}

impl Default for Cigar {
    fn default() -> Self {
        Self::new()
    }
}

impl Cigar {
    pub fn new() -> Self {
        Self {
            query_extension_distance: 0,
            target_extension_distance: 0,
            cigar_data: Vec::new(),
            score: 0,
            peak_score: 0,
            peak_score_cigar_index: 0,
            peak_score_anchor_index: 0,
        }
    }

    pub fn with_reserve(reserve_size: usize) -> Self {
        let mut cigar = Self::new();
        cigar.cigar_data.reserve(reserve_size);
        cigar
    }

    pub fn reserve_cigar_space(&mut self, reserve_size: usize) {
        self.cigar_data.reserve(reserve_size);
    }

    pub fn extend_cigar(&mut self, other_vector: &[(i32, char)]) {
        self.cigar_data.extend_from_slice(other_vector);
    }

    pub fn extend_cigar_op(&mut self, length: u32, cigar_operation: char) {
        self.cigar_data.push((length as i32, cigar_operation));
    }

    pub fn query_extension_distance(&self) -> i32 {
        self.query_extension_distance
    }

    pub fn target_extension_distance(&self) -> i32 {
        self.target_extension_distance
    }

    pub fn set_max_values(&mut self, query_start: i32, target_start: i32) {
        self.query_extension_distance = query_start;
        self.target_extension_distance = target_start;
    }

    pub fn get_cigar_data_const(&self) -> &[(i32, char)] {
        &self.cigar_data
    }

    pub fn get_cigar_data(&mut self) -> &mut Vec<(i32, char)> {
        &mut self.cigar_data
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum AlignmentStatus {
    NotDropped = 0,
    Dropped = 1,
    NegativeScore = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DnaScoring {
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

pub fn compute_wfa_cigar_from_string(
    scoring: DnaScoring,
    cigar: &str,
    extension: &mut Cigar,
    left: bool,
) -> Result<AlignmentStatus, String> {
    let mut cigar_data = Vec::new();
    let mut max_query = -1;
    let mut max_target = -1;
    let mut steps = 0;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            steps = steps * 10 + (c as i32 - '0' as i32);
            continue;
        }
        cigar_data.push((steps, c));
        match c {
            '=' => {
                extension.score += steps * scoring.reward;
                max_query += steps;
                max_target += steps;
            }
            'X' => {
                extension.score += steps * scoring.penalty;
                max_query += steps;
                max_target += steps;
            }
            'I' => {
                extension.score -= scoring.gap_open + (steps * scoring.gap_extend);
                max_query += steps;
            }
            'D' => {
                extension.score -= scoring.gap_open + (steps * scoring.gap_extend);
                max_target += steps;
            }
            _ => return Err(format!("WFA Cigar_short: Invalid Cigar_short Symbol {}", c)),
        }

        steps = 0;
    }

    if left {
        cigar_data.reverse();
        extension.set_max_values(max_query, max_target);
    } else if extension.score < 1 {
        return Ok(AlignmentStatus::NegativeScore);
    }
    extension.extend_cigar(&cigar_data);

    Ok(AlignmentStatus::NotDropped)
}

pub fn build_hsp_from_cigar(
    cigar: &Cigar,
    target: &[crate::basic::value::Letter],
    query: &[crate::basic::value::Letter],
    first_anchor_i: i32,
    first_anchor_j: i32,
    is_reverse: bool,
) -> crate::align::hsp::Hsp {
    use crate::basic::packed_transcript::EditOperation;
    use crate::util::interval::Interval;

    let mut align_hsp = crate::align::hsp::Hsp::new();

    let mut query_pos = first_anchor_i - cigar.query_extension_distance() - 1;
    let mut target_pos = first_anchor_j - cigar.target_extension_distance() - 1;
    align_hsp.query_range.begin = query_pos;
    align_hsp.subject_range.begin = target_pos;

    for operation in cigar.get_cigar_data_const() {
        match operation.1 {
            'M' | '=' | 'X' => {
                for _ in 0..operation.0 {
                    align_hsp.push_match(
                        target[target_pos as usize],
                        query[query_pos as usize],
                        true,
                    );
                    target_pos += 1;
                    query_pos += 1;
                }
            }
            'D' => {
                let end = (target_pos + operation.0) as usize;
                align_hsp.push_gap(EditOperation::Deletion, operation.0, &target[..end]);
                target_pos += operation.0;
            }
            'I' => {
                let end = (query_pos + operation.0).max(0) as usize;
                align_hsp.push_gap(EditOperation::Insertion, operation.0, &query[..end]);
                query_pos += operation.0;
            }
            _ => {}
        }
    }

    align_hsp.score = cigar.score;
    align_hsp.query_range.end = query_pos;
    align_hsp.subject_range.end = target_pos;
    align_hsp.transcript.push_terminator();
    align_hsp.target_seq = target.to_vec();
    align_hsp.query_source_range = align_hsp.query_range;
    align_hsp.subject_source_range = if is_reverse {
        Interval::new(align_hsp.subject_range.end, align_hsp.subject_range.begin)
    } else {
        Interval::new(align_hsp.subject_range.begin, align_hsp.subject_range.end)
    };
    align_hsp.frame = is_reverse as i32;

    align_hsp
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::EditOperation;
    use crate::basic::value::Letter;
    use crate::util::interval::Interval;

    fn dna(s: &[u8]) -> Vec<Letter> {
        s.iter().map(|&x| x as Letter).collect()
    }

    #[test]
    fn test_cigar_methods() {
        let mut c = Cigar::with_reserve(8);
        c.reserve_cigar_space(16);
        c.extend_cigar(&[(2, 'M'), (1, 'I')]);
        c.extend_cigar_op(3, 'D');
        c.set_max_values(4, 5);

        assert_eq!(c.query_extension_distance(), 4);
        assert_eq!(c.target_extension_distance(), 5);
        assert_eq!(c.get_cigar_data_const(), &[(2, 'M'), (1, 'I'), (3, 'D')]);
        c.get_cigar_data().push((1, 'X'));
        assert_eq!(c.get_cigar_data_const()[3], (1, 'X'));
    }

    #[test]
    fn test_build_hsp_from_cigar() {
        let target = dna(b"AACCGGTT");
        let query = dna(b"AACAGGTT");
        let mut cigar = Cigar::new();
        cigar.set_max_values(0, 0);
        cigar.extend_cigar(&[(3, 'M'), (1, 'I'), (2, 'M'), (1, 'D')]);
        cigar.score = 12;

        let hsp = build_hsp_from_cigar(&cigar, &target, &query, 1, 1, true);

        assert_eq!(hsp.score, 12);
        assert_eq!(hsp.query_range, Interval::new(0, 6));
        assert_eq!(hsp.subject_range, Interval::new(0, 6));
        assert_eq!(hsp.query_source_range, Interval::new(0, 6));
        assert_eq!(hsp.subject_source_range, Interval::new(6, 0));
        assert_eq!(hsp.frame, 1);
        assert_eq!(hsp.length, 7);
        assert_eq!(hsp.gaps, 2);
        assert_eq!(hsp.gap_openings, 2);
        let ops = hsp.transcript.iter().collect::<Vec<_>>();
        assert_eq!(ops[0].op, EditOperation::Match);
        assert!(hsp.transcript.data().last().unwrap().is_terminator());
    }

    #[test]
    fn test_compute_wfa_cigar_from_string_scores_and_reverses_left_extension() {
        let scoring = DnaScoring {
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 1,
        };
        let mut cigar = Cigar::new();

        let status = compute_wfa_cigar_from_string(scoring, "3=2X4I5D", &mut cigar, true).unwrap();

        assert_eq!(status, AlignmentStatus::NotDropped);
        assert_eq!(cigar.score, 6 - 6 - 9 - 10);
        assert_eq!(cigar.query_extension_distance(), 8);
        assert_eq!(cigar.target_extension_distance(), 9);
        assert_eq!(
            cigar.get_cigar_data_const(),
            &[(5, 'D'), (4, 'I'), (2, 'X'), (3, '=')]
        );
    }

    #[test]
    fn test_compute_wfa_cigar_from_string_negative_and_invalid() {
        let scoring = DnaScoring {
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 1,
        };
        let mut cigar = Cigar::new();

        let status = compute_wfa_cigar_from_string(scoring, "1X", &mut cigar, false).unwrap();
        assert_eq!(status, AlignmentStatus::NegativeScore);
        assert!(cigar.get_cigar_data_const().is_empty());

        let err =
            compute_wfa_cigar_from_string(scoring, "1Z", &mut Cigar::new(), false).unwrap_err();
        assert!(err.contains("Invalid Cigar_short Symbol Z"));
    }
}
