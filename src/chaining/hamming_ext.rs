//! Approximate Hamming extension and diagonal pre-filters.
//!
//! This mirrors `diamond/src/chaining/hamming_ext.cpp`. Configuration that is
//! global in C++ is explicit here, keeping the numerical and ordering behavior
//! testable without process-global state.

use crate::dp::ungapped::DiagonalSegment;
use crate::stats;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::hsp::{Anchor, ApproxHsp};
use crate::util::misc::safe_cast;

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

/// Return the first diagonal segment satisfying identity, coverage, and E-value filters.
pub fn find_aln(
    segments: &mut [DiagonalSegment],
    qlen: i32,
    tlen: i32,
    config: &HammingExtConfig,
    score_matrix: &ScoreMatrix,
) -> ApproxHsp {
    for segment in segments {
        let evalue = score_matrix.evalue(segment.score, qlen as u32, tlen as u32);
        if (segment.id_percent() >= config.approx_min_id
            || stats::approx_id(segment.score, segment.len, segment.len) >= config.approx_min_id)
            && ((config.query_or_target_cover > 0.0
                && segment.cov_percent(qlen).max(segment.cov_percent(tlen))
                    >= config.query_or_target_cover)
                || (config.query_or_target_cover == 0.0
                    && segment.cov_percent(qlen) >= config.query_cover
                    && segment.cov_percent(tlen) >= config.subject_cover))
            && evalue <= config.max_evalue
        {
            return ApproxHsp::from_parts(
                0,
                0,
                segment.score,
                0,
                segment.query_range(),
                segment.subject_range(),
                Anchor::from_diagonal_segment(segment.clone()),
                evalue,
            );
        }
    }
    ApproxHsp::new(0, 0)
}

/// Apply aggregate coverage and identity filters to score-sorted segments.
pub fn filter(
    segments: &mut [DiagonalSegment],
    qlen: i32,
    tlen: i32,
    use_cov_filter: bool,
    config: &HammingExtConfig,
) -> ApproxHsp {
    const TOLERANCE_FACTOR: f64 = 1.1;
    const ID_MIN_COV: f64 = 80.0;

    // `slice::sort_by` is stable, matching C++ `stable_sort`.
    segments.sort_by(|x, y| y.score.cmp(&x.score));
    let mut identities = 0;
    let mut len = 0;
    // C++ uses `safe_cast<Loc>` here.  In particular, it reports an
    // out-of-range sequence length instead of applying Rust's saturating
    // float-to-integer cast.
    let query_tolerance = safe_cast::<i32, f64>(qlen as f64 * TOLERANCE_FACTOR)
        .unwrap_or_else(|message| panic!("{message}"));
    let target_tolerance = safe_cast::<i32, f64>(tlen as f64 * TOLERANCE_FACTOR)
        .unwrap_or_else(|message| panic!("{message}"));
    for segment in segments {
        if len + segment.len > query_tolerance || len + segment.len > target_tolerance {
            continue;
        }
        identities += segment.identities;
        len += segment.len;
    }
    let query_coverage = len as f64 / qlen as f64 * 100.0;
    let target_coverage = len as f64 / tlen as f64 * 100.0;
    if let Some(diag_filter_cov) = config.diag_filter_cov {
        if use_cov_filter
            && ((config.query_or_target_cover > 0.0
                && query_coverage.max(target_coverage) < diag_filter_cov)
                || (config.query_cover > 0.0 && query_coverage < diag_filter_cov)
                || (config.subject_cover > 0.0 && target_coverage < diag_filter_cov))
        {
            return ApproxHsp::new(0, -1);
        }
    }
    if let Some(diag_filter_id) = config.diag_filter_id {
        // The C++ 0/0 result is NaN and therefore does not compare less than
        // the cutoff. The explicit guard preserves that result without a
        // floating-point invalid operation.
        if query_coverage.max(target_coverage) >= ID_MIN_COV
            && len > 0
            && identities as f64 / len as f64 * 100.0 < diag_filter_id
        {
            return ApproxHsp::new(0, -1);
        }
    }
    ApproxHsp::new(0, 0)
}

/// Try approximate extension first, then run configured diagonal filters.
pub fn hamming_ext(
    segments: &mut [DiagonalSegment],
    qlen: i32,
    tlen: i32,
    use_cov_filter: bool,
    config: &HammingExtConfig,
    score_matrix: &ScoreMatrix,
) -> ApproxHsp {
    if config.hamming_ext {
        let hsp = find_aln(segments, qlen, tlen, config, score_matrix);
        if hsp.score > 0 {
            return hsp;
        }
    }
    if config.diag_filter_cov.is_some() || config.diag_filter_id.is_some() {
        return filter(segments, qlen, tlen, use_cov_filter, config);
    }
    ApproxHsp::new(0, 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::interval::Interval;

    fn score_matrix() -> ScoreMatrix {
        ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap()
    }

    #[test]
    fn find_aln_returns_first_passing_segment_without_reordering() {
        let mut segments = vec![
            DiagonalSegment::with_identities(2, 3, 8, 40, 8),
            DiagonalSegment::with_identities(4, 5, 10, 80, 10),
        ];
        let original = segments.clone();
        let config = HammingExtConfig {
            hamming_ext: true,
            approx_min_id: 90.0,
            query_or_target_cover: 40.0,
            ..HammingExtConfig::default()
        };

        let hsp = hamming_ext(&mut segments, 20, 40, true, &config, &score_matrix());
        assert_eq!(hsp.score, 40);
        assert_eq!(hsp.query_range, Interval::new(2, 10));
        assert_eq!(hsp.subject_range, Interval::new(3, 11));
        assert_eq!(segments, original);
    }

    #[test]
    fn filter_sorts_stably_and_enforces_tolerance_and_identity() {
        let mut segments = vec![
            DiagonalSegment::with_identities(0, 0, 5, 20, 0),
            DiagonalSegment::with_identities(5, 5, 6, 30, 6),
            DiagonalSegment::with_identities(11, 11, 1, 10, 1),
        ];
        let config = HammingExtConfig {
            diag_filter_id: Some(60.0),
            ..HammingExtConfig::default()
        };

        // The first two score-sorted lengths exactly fill floor(10 * 1.1).
        // The final length is skipped, leaving 6/11 identities (< 60%).
        let hsp = filter(&mut segments, 10, 10, true, &config);
        assert_eq!(hsp.score, -1);
        assert_eq!(
            segments.iter().map(|s| s.score).collect::<Vec<_>>(),
            [30, 20, 10]
        );
    }

    #[test]
    fn use_cov_filter_only_gates_coverage_not_identity() {
        let config = HammingExtConfig {
            query_cover: 90.0,
            subject_cover: 90.0,
            diag_filter_cov: Some(90.0),
            ..HammingExtConfig::default()
        };
        let segment = DiagonalSegment::with_identities(0, 0, 2, 10, 2);

        assert_eq!(
            filter(&mut [segment.clone()], 10, 10, true, &config).score,
            -1
        );
        assert_eq!(filter(&mut [segment], 10, 10, false, &config).score, 0);
    }

    #[test]
    fn find_aln_honors_approximate_identity_and_separate_coverages() {
        // The recorded identity is deliberately too low, while the score
        // estimate is high enough: C++ accepts when either identity estimate
        // passes.  Separate query/subject thresholds must both pass.
        let segment = DiagonalSegment::with_identities(3, 7, 8, 40, 0);
        let config = HammingExtConfig {
            approx_min_id: 80.0,
            query_cover: 40.0,
            subject_cover: 20.0,
            ..HammingExtConfig::default()
        };
        let hsp = find_aln(&mut [segment.clone()], 20, 40, &config, &score_matrix());
        assert_eq!(hsp.score, 40);
        assert_eq!(hsp.max_diag.segment, segment);

        let reject_subject = HammingExtConfig {
            subject_cover: 21.0,
            ..config
        };
        assert_eq!(
            find_aln(
                &mut [DiagonalSegment::with_identities(3, 7, 8, 40, 0)],
                20,
                40,
                &reject_subject,
                &score_matrix(),
            )
            .score,
            0
        );
    }

    #[test]
    fn query_or_target_coverage_uses_the_better_side() {
        let config = HammingExtConfig {
            approx_min_id: 0.0,
            query_or_target_cover: 50.0,
            // These are ignored while the OR threshold is enabled.
            query_cover: 100.0,
            subject_cover: 100.0,
            ..HammingExtConfig::default()
        };
        let hsp = find_aln(
            &mut [DiagonalSegment::with_identities(0, 0, 5, 20, 5)],
            10,
            100,
            &config,
            &score_matrix(),
        );
        assert_eq!(hsp.score, 20);
    }

    #[test]
    fn evalue_rejection_falls_through_to_diagonal_filter() {
        let config = HammingExtConfig {
            hamming_ext: true,
            approx_min_id: 0.0,
            query_cover: 1.0,
            max_evalue: -1.0,
            diag_filter_cov: Some(90.0),
            ..HammingExtConfig::default()
        };
        let mut segments = [DiagonalSegment::with_identities(0, 0, 5, 40, 5)];

        // E-values are non-negative, so approximate extension cannot return
        // this segment. `hamming_ext` must then execute the configured filter.
        assert_eq!(
            hamming_ext(&mut segments, 10, 10, true, &config, &score_matrix()).score,
            -1
        );
    }

    #[test]
    fn score_sort_is_stable_for_equal_segments() {
        let mut segments = [
            DiagonalSegment::new(1, 1, 1, 10),
            DiagonalSegment::new(2, 2, 1, 20),
            DiagonalSegment::new(3, 3, 1, 10),
        ];
        let _ = filter(&mut segments, 10, 10, true, &HammingExtConfig::default());
        assert_eq!(
            segments.iter().map(|segment| segment.i).collect::<Vec<_>>(),
            [2, 1, 3]
        );
    }

    #[test]
    #[should_panic(expected = "safe_cast: out of range (float -> int)")]
    fn filter_rejects_out_of_range_tolerance_like_safe_cast() {
        let _ = filter(&mut [], i32::MAX, 1, true, &HammingExtConfig::default());
    }
}
