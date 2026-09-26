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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KswMode {
    Left,
    Right,
    Global,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TraceState {
    Match,
    Insertion,
    Deletion,
}

/// Safe scalar counterpart of the affine-gap engines used by KSW2 and WFA.
/// It starts at (0, 0); extension mode ends at the best-scoring cell, while
/// global mode ends at the bottom-right cell.
fn affine_cigar(
    query: &[u8],
    target: &[u8],
    scoring: DnaScoring,
    global: bool,
    band: Option<i32>,
) -> (i32, usize, usize, Vec<(i32, char)>) {
    let cols = target.len() + 1;
    let cells = (query.len() + 1) * cols;
    let neg = i32::MIN / 4;
    let mut m = vec![neg; cells];
    let mut ins = vec![neg; cells];
    let mut del = vec![neg; cells];
    let mut tm = vec![TraceState::Match; cells];
    let mut ti = vec![TraceState::Match; cells];
    let mut td = vec![TraceState::Match; cells];
    m[0] = 0;
    let index = |i: usize, j: usize| i * cols + j;

    for i in 0..=query.len() {
        let j_begin = band
            .map(|width| i.saturating_sub(width.max(0) as usize))
            .unwrap_or(0);
        let j_end = band
            .map(|width| target.len().min(i.saturating_add(width.max(0) as usize)))
            .unwrap_or(target.len());
        for j in j_begin..=j_end {
            if i == 0 && j == 0 {
                continue;
            }
            let at = index(i, j);
            if i > 0 {
                let prev = index(i - 1, j);
                let open_m = m[prev].saturating_sub(scoring.gap_open + scoring.gap_extend);
                let open_d = del[prev].saturating_sub(scoring.gap_open + scoring.gap_extend);
                let extend = ins[prev].saturating_sub(scoring.gap_extend);
                let (score, state) = [
                    (open_m, TraceState::Match),
                    (open_d, TraceState::Deletion),
                    (extend, TraceState::Insertion),
                ]
                .into_iter()
                .max_by_key(|x| x.0)
                .unwrap();
                ins[at] = score;
                ti[at] = state;
            }
            if j > 0 {
                let prev = index(i, j - 1);
                let open_m = m[prev].saturating_sub(scoring.gap_open + scoring.gap_extend);
                let open_i = ins[prev].saturating_sub(scoring.gap_open + scoring.gap_extend);
                let extend = del[prev].saturating_sub(scoring.gap_extend);
                let (score, state) = [
                    (open_m, TraceState::Match),
                    (open_i, TraceState::Insertion),
                    (extend, TraceState::Deletion),
                ]
                .into_iter()
                .max_by_key(|x| x.0)
                .unwrap();
                del[at] = score;
                td[at] = state;
            }
            if i > 0 && j > 0 {
                let prev = index(i - 1, j - 1);
                let (base, state) = [
                    (m[prev], TraceState::Match),
                    (ins[prev], TraceState::Insertion),
                    (del[prev], TraceState::Deletion),
                ]
                .into_iter()
                .max_by_key(|x| x.0)
                .unwrap();
                m[at] = base.saturating_add(if query[i - 1] == target[j - 1] {
                    scoring.reward
                } else {
                    scoring.penalty
                });
                tm[at] = state;
            }
        }
    }

    let mut best = if global {
        let at = index(query.len(), target.len());
        [
            (m[at], TraceState::Match),
            (ins[at], TraceState::Insertion),
            (del[at], TraceState::Deletion),
        ]
        .into_iter()
        .max_by_key(|x| x.0)
        .map(|(score, state)| (score, query.len(), target.len(), state))
        .unwrap()
    } else {
        (0, 0, 0, TraceState::Match)
    };
    if !global {
        for i in 0..=query.len() {
            for j in 0..=target.len() {
                let at = index(i, j);
                for (score, state) in [
                    (m[at], TraceState::Match),
                    (ins[at], TraceState::Insertion),
                    (del[at], TraceState::Deletion),
                ] {
                    if score > best.0 {
                        best = (score, i, j, state);
                    }
                }
            }
        }
    }

    let (score, end_i, end_j, mut state) = best;
    let (mut i, mut j) = (end_i, end_j);
    let mut raw = Vec::new();
    while i > 0 || j > 0 {
        let at = index(i, j);
        match state {
            TraceState::Match if i > 0 && j > 0 => {
                raw.push(if query[i - 1] == target[j - 1] {
                    '='
                } else {
                    'X'
                });
                state = tm[at];
                i -= 1;
                j -= 1;
            }
            TraceState::Insertion if i > 0 => {
                raw.push('I');
                state = ti[at];
                i -= 1;
            }
            TraceState::Deletion if j > 0 => {
                raw.push('D');
                state = td[at];
                j -= 1;
            }
            _ => break,
        }
    }
    raw.reverse();
    let mut cigar: Vec<(i32, char)> = Vec::new();
    for op in raw {
        if let Some(last) = cigar.last_mut().filter(|last| last.1 == op) {
            last.0 += 1;
        } else {
            cigar.push((1, op));
        }
    }
    (score, end_i, end_j, cigar)
}

fn cigar_string(cigar: &[(i32, char)]) -> String {
    cigar.iter().map(|(n, op)| format!("{n}{op}")).collect()
}

pub fn compute_wfa_extension(
    query_sequence: &str,
    target_sequence: &str,
    scoring: DnaScoring,
    _band: i32,
) -> (AlignmentStatus, String) {
    let (_, _, _, cigar) = affine_cigar(
        query_sequence.as_bytes(),
        target_sequence.as_bytes(),
        scoring,
        false,
        None,
    );
    (AlignmentStatus::NotDropped, cigar_string(&cigar))
}

pub fn compute_wfa_global(
    query_sequence: &str,
    target_sequence: &str,
    scoring: DnaScoring,
    _band: i32,
) -> (AlignmentStatus, String) {
    let (_, _, _, cigar) = affine_cigar(
        query_sequence.as_bytes(),
        target_sequence.as_bytes(),
        scoring,
        true,
        None,
    );
    (AlignmentStatus::NotDropped, cigar_string(&cigar))
}

pub fn compute_wfa_cigar(
    scoring: DnaScoring,
    query_sequence: &str,
    extension: &mut Cigar,
    left: bool,
    global: bool,
    target_sequence: &str,
    band: i32,
) -> Result<AlignmentStatus, String> {
    let (status, cigar) = if global {
        compute_wfa_global(query_sequence, target_sequence, scoring, band)
    } else {
        compute_wfa_extension(query_sequence, target_sequence, scoring, band)
    };
    compute_wfa_cigar_from_string(scoring, &cigar, extension, left)?;
    Ok(status)
}

pub fn compute_ksw_cigar(
    target_sequence: &[u8],
    query_sequence: &[u8],
    scoring: DnaScoring,
    mode: KswMode,
    extension: &mut Cigar,
    _zdrop: i32,
    band: i32,
) -> AlignmentStatus {
    let global = mode == KswMode::Global;
    let (score, query_end, target_end, cigar) =
        affine_cigar(query_sequence, target_sequence, scoring, global, Some(band));
    extension.score += score;
    if mode == KswMode::Left {
        extension.set_max_values(query_end as i32 - 1, target_end as i32 - 1);
    } else if extension.score < 1 {
        return AlignmentStatus::NegativeScore;
    }
    for (length, op) in cigar {
        extension.extend_cigar_op(
            length as u32,
            if matches!(op, '=' | 'X') { 'M' } else { op },
        );
    }
    AlignmentStatus::NotDropped
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
    score_builder: &crate::dna::build_score::Blastn_Score,
    max_evalue: f64,
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
    align_hsp.bit_score = score_builder.blast_bit_score(align_hsp.score);
    align_hsp.evalue = score_builder.blast_e_value(align_hsp.score, query.len() as i32);
    if align_hsp.evalue >= max_evalue {
        return align_hsp;
    }
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

        let score_builder = crate::dna::build_score::Blastn_Score::new_with_parameters(
            2, -3, 5, 1, 1000, 2, 0.625, 0.41, 1.1,
        );
        let hsp = build_hsp_from_cigar(
            &cigar,
            &target,
            &query,
            1,
            1,
            true,
            &score_builder,
            f64::INFINITY,
        );

        assert_eq!(hsp.score, 12);
        assert_eq!(hsp.query_range, Interval::new(0, 6));
        assert_eq!(hsp.subject_range, Interval::new(0, 6));
        assert_eq!(hsp.query_source_range, Interval::new(0, 6));
        assert_eq!(hsp.subject_source_range, Interval::new(6, 0));
        assert_eq!(hsp.frame, 1);
        assert_eq!(hsp.bit_score, score_builder.blast_bit_score(12));
        assert_eq!(hsp.evalue, score_builder.blast_e_value(12, 8));
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

    #[test]
    fn test_wfa_global_and_extension_entry_points() {
        let scoring = DnaScoring {
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 1,
        };
        let (_, global) = compute_wfa_global("ACGT", "ACCT", scoring, 20);
        assert_eq!(global, "2=1X1=");

        let (_, extension) = compute_wfa_extension("ACGTTT", "ACGTAA", scoring, 20);
        assert_eq!(extension, "4=");

        let mut parsed = Cigar::new();
        assert_eq!(
            compute_wfa_cigar(scoring, "ACGT", &mut parsed, false, true, "ACCT", 20).unwrap(),
            AlignmentStatus::NotDropped
        );
        assert_eq!(parsed.score, 3);
        assert_eq!(
            parsed.get_cigar_data_const(),
            &[(2, '='), (1, 'X'), (1, '=')]
        );
    }

    #[test]
    fn test_ksw_entry_point_modes_and_negative_score() {
        let scoring = DnaScoring {
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 1,
        };
        let mut left = Cigar::new();
        assert_eq!(
            compute_ksw_cigar(
                b"ACGTAA",
                b"ACGTTT",
                scoring,
                KswMode::Left,
                &mut left,
                40,
                40
            ),
            AlignmentStatus::NotDropped
        );
        assert_eq!(left.score, 8);
        assert_eq!(
            (
                left.query_extension_distance(),
                left.target_extension_distance()
            ),
            (3, 3)
        );
        assert_eq!(left.get_cigar_data_const(), &[(4, 'M')]);

        let mut right = Cigar::new();
        assert_eq!(
            compute_ksw_cigar(b"T", b"A", scoring, KswMode::Right, &mut right, 40, 40),
            AlignmentStatus::NegativeScore
        );
        assert!(right.get_cigar_data_const().is_empty());
    }
}
