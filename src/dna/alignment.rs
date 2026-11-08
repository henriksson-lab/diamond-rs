pub const KSW2_ZDROP_EXTENSION: i32 = 40;
use std::sync::Arc;

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

/// Typed compatibility API; checked production callers use the `try_` entry.
/// The Source band setters are commented out, so `band` has no active effect.
pub fn compute_wfa_extension(
    query_sequence: &str,
    target_sequence: &str,
    scoring: DnaScoring,
    band: i32,
) -> (AlignmentStatus, String) {
    try_compute_wfa_extension(query_sequence, target_sequence, scoring, band)
        .unwrap_or_else(|_| (AlignmentStatus::Dropped, String::new()))
}

pub fn try_compute_wfa_extension(
    query_sequence: &str,
    target_sequence: &str,
    scoring: DnaScoring,
    _band: i32,
) -> Result<(AlignmentStatus, String), String> {
    wavefront_alignment(query_sequence, target_sequence, scoring, true)
}

pub fn compute_wfa_global(
    query_sequence: &str,
    target_sequence: &str,
    scoring: DnaScoring,
    band: i32,
) -> (AlignmentStatus, String) {
    try_compute_wfa_global(query_sequence, target_sequence, scoring, band)
        .unwrap_or_else(|_| (AlignmentStatus::Dropped, String::new()))
}

pub fn try_compute_wfa_global(
    query_sequence: &str,
    target_sequence: &str,
    scoring: DnaScoring,
    _band: i32,
) -> Result<(AlignmentStatus, String), String> {
    wavefront_alignment(query_sequence, target_sequence, scoring, false)
}

fn wavefront_alignment(
    query: &str,
    target: &str,
    scoring: DnaScoring,
    extension: bool,
) -> Result<(AlignmentStatus, String), String> {
    // Source passes target as WFA pattern and query as text. Its wavefront
    // match cost is zero; reward is used later while parsing the CIGAR.
    let result = super::wavefront::align(
        target.as_bytes(),
        query.as_bytes(),
        scoring
            .penalty
            .checked_neg()
            .ok_or("WFA mismatch penalty overflow")?,
        scoring.gap_open,
        scoring.gap_extend,
        extension,
    )?;
    let status = if result.status == 0 {
        AlignmentStatus::NotDropped
    } else {
        AlignmentStatus::Dropped
    };
    Ok((status, result.cigar_string()))
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
        try_compute_wfa_global(query_sequence, target_sequence, scoring, band)?
    } else {
        try_compute_wfa_extension(query_sequence, target_sequence, scoring, band)?
    };
    let parsed_status = compute_wfa_cigar_from_string(scoring, &cigar, extension, left)?;
    // Source returns NEGATIVE_SCORE before appending a non-left CIGAR,
    // independently of the underlying wavefront completion status.
    if parsed_status == AlignmentStatus::NegativeScore {
        return Ok(parsed_status);
    }
    Ok(status)
}

/// Compatibility entry point accepting either Source's numeric bases or the
/// older typed Rust API's ASCII DNA convenience. Production adapters should
/// use `try_compute_ksw_cigar` for checked, encoded input and error propagation.
pub fn compute_ksw_cigar(
    target_sequence: &[u8],
    query_sequence: &[u8],
    scoring: DnaScoring,
    mode: KswMode,
    extension: &mut Cigar,
    zdrop: i32,
    band: i32,
) -> AlignmentStatus {
    let target = encoded_or_ascii(target_sequence);
    let query = encoded_or_ascii(query_sequence);
    try_compute_ksw_cigar(&target, &query, scoring, mode, extension, zdrop, band)
        .unwrap_or(AlignmentStatus::NegativeScore)
}

pub(crate) fn encoded_or_ascii(sequence: &[u8]) -> Vec<u8> {
    if sequence.iter().any(|&base| base >= 5)
        && sequence.iter().all(|base| b"ACGTNacgtn".contains(base))
    {
        sequence
            .iter()
            .map(|base| match base.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect()
    } else {
        sequence.to_vec()
    }
}

pub(crate) fn ksw_extension_result(
    target: &[u8],
    query: &[u8],
    scoring: DnaScoring,
    mode: KswMode,
    zdrop: i32,
    band: i32,
) -> Result<super::ksw::Extz, String> {
    let reward = scoring.reward as i8;
    let penalty = if scoring.penalty < 0 {
        scoring.penalty
    } else {
        scoring.penalty.wrapping_neg()
    } as i8;
    let mut matrix = [0i8; 25];
    for i in 0..4 {
        for j in 0..4 {
            matrix[i * 5 + j] = if i == j { reward } else { penalty };
        }
    }
    let mut flags = super::ksw::APPROX_MAX | super::ksw::APPROX_DROP;
    if mode != KswMode::Global {
        flags |= super::ksw::EXTZ_ONLY;
    }
    if mode == KswMode::Left {
        flags |= super::ksw::REV_CIGAR;
    }
    if cfg!(target_vendor = "apple") {
        // DIAMOND's __APPLE__ branch calls the scalar engine.
        super::ksw::extz(
            query,
            target,
            5,
            &matrix,
            scoring.gap_open as i8,
            scoring.gap_extend as i8,
            band,
            zdrop,
            flags,
        )
    } else {
        // Interpret the GNU SSE engine's byte operations on every Rust
        // target. The source preprocessor selects these two intrinsic paths;
        // an explicit low-level Flavor remains available for conformance.
        let flavor = if cfg!(target_feature = "sse4.1") {
            super::ksw_sse::Flavor::Sse41
        } else {
            super::ksw_sse::Flavor::Sse2
        };
        super::ksw_sse::extz2(
            query,
            target,
            5,
            &matrix,
            scoring.gap_open as i8,
            scoring.gap_extend as i8,
            band,
            zdrop,
            5,
            flags,
            flavor,
        )
    }
}

/// Source `compute_ksw_cigar`, selecting its actual scalar or GNU state machine.
/// Inputs are encoded A/C/G/T/N (0..4). Unlike the compatibility wrapper,
/// malformed alphabets, exhausted bands, and allocation errors propagate.
pub fn try_compute_ksw_cigar(
    target_sequence: &[u8],
    query_sequence: &[u8],
    scoring: DnaScoring,
    mode: KswMode,
    extension: &mut Cigar,
    zdrop: i32,
    band: i32,
) -> Result<AlignmentStatus, String> {
    let result = ksw_extension_result(target_sequence, query_sequence, scoring, mode, zdrop, band)?;
    extension.score += if mode != KswMode::Global || result.zdropped {
        result.max
    } else {
        result.score
    };
    if mode == KswMode::Left {
        extension.set_max_values(result.max_q, result.max_t);
    } else if extension.score < 1 {
        return Ok(AlignmentStatus::NegativeScore);
    }
    for (length, op) in result.cigar {
        extension.extend_cigar_op(length as u32, op);
    }
    Ok(if result.zdropped {
        AlignmentStatus::Dropped
    } else {
        AlignmentStatus::NotDropped
    })
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
    let build_started = std::time::Instant::now();
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
    align_hsp.target_seq = Arc::from(target);
    align_hsp.query_source_range = align_hsp.query_range;
    align_hsp.subject_source_range = if is_reverse {
        Interval::new(align_hsp.subject_range.end, align_hsp.subject_range.begin)
    } else {
        Interval::new(align_hsp.subject_range.begin, align_hsp.subject_range.end)
    };
    align_hsp.frame = is_reverse as i32;

    super::native::record_time(2, build_started.elapsed());
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
