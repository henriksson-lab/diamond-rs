use crate::basic::packed_transcript::PackedTranscript;
use crate::basic::value::{Letter, Score};
use crate::util::interval::Interval;

/// High-Scoring Pair — a single alignment between query and subject.
///
/// This is the core alignment result type, containing the alignment
/// coordinates, scores, statistics, and optionally the alignment transcript.
#[derive(Debug, Clone)]
pub struct Hsp {
    /// Whether this HSP has a full backtrace (alignment transcript).
    pub backtraced: bool,
    /// Whether this is a seed-only hit (no extension).
    pub seed_only: bool,
    /// Raw alignment score.
    pub score: Score,
    /// E-value.
    pub evalue: f64,
    /// Bit score.
    pub bit_score: f64,
    /// Corrected bit score (composition-adjusted).
    pub corrected_bit_score: f64,
    /// Approximate identity.
    pub approx_id: f64,
    /// Reading frame (0 for blastp, 0-5 for blastx).
    pub frame: i32,
    /// Alignment length.
    pub length: i32,
    /// Number of identical positions.
    pub identities: i32,
    /// Number of mismatches.
    pub mismatches: i32,
    /// Number of positive-scoring positions.
    pub positives: i32,
    /// Number of gap openings.
    pub gap_openings: i32,
    /// Number of gaps (total gap positions).
    pub gaps: i32,
    /// Query alignment range (in query coordinates).
    pub query_range: Interval,
    /// Subject alignment range.
    pub subject_range: Interval,
    /// Query range in source coordinates (for blastx: nucleotide coords).
    pub query_source_range: Interval,
    /// Subject range in source coordinates.
    pub subject_source_range: Interval,
    /// Target sequence data (for output).
    pub target_seq: Vec<Letter>,
    /// Alignment transcript (CIGAR-like edit operations).
    pub transcript: PackedTranscript,
    /// Internal: swipe target index.
    pub swipe_target: i32,
    /// Internal: diagonal band begin.
    pub d_begin: i32,
    /// Internal: diagonal band end.
    pub d_end: i32,
}

impl Default for Hsp {
    fn default() -> Self {
        Hsp {
            backtraced: false,
            seed_only: false,
            score: 0,
            evalue: f64::MAX,
            bit_score: 0.0,
            corrected_bit_score: 0.0,
            approx_id: 0.0,
            frame: 0,
            length: 0,
            identities: 0,
            mismatches: 0,
            positives: 0,
            gap_openings: 0,
            gaps: 0,
            query_range: Interval::default(),
            subject_range: Interval::default(),
            query_source_range: Interval::default(),
            subject_source_range: Interval::default(),
            target_seq: Vec::new(),
            transcript: PackedTranscript::new(),
            swipe_target: 0,
            d_begin: 0,
            d_end: 0,
        }
    }
}

impl Hsp {
    pub fn new() -> Self {
        Self::default()
    }

    /// Percent identity.
    pub fn pident(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            100.0 * self.identities as f64 / self.length as f64
        }
    }

    /// Percent positive.
    pub fn ppos(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            100.0 * self.positives as f64 / self.length as f64
        }
    }

    /// Query coverage of this HSP.
    pub fn query_cover(&self, query_len: i32) -> f64 {
        if query_len == 0 {
            0.0
        } else {
            100.0 * self.query_range.length() as f64 / query_len as f64
        }
    }

    /// Subject coverage of this HSP.
    pub fn subject_cover(&self, subject_len: i32) -> f64 {
        if subject_len == 0 {
            0.0
        } else {
            100.0 * self.subject_range.length() as f64 / subject_len as f64
        }
    }

    /// Whether this HSP passes the given filters.
    pub fn passes_filters(
        &self,
        max_evalue: f64,
        min_id: f64,
        min_query_cover: f64,
        min_subject_cover: f64,
        query_len: i32,
        subject_len: i32,
    ) -> bool {
        self.evalue <= max_evalue
            && self.pident() >= min_id
            && self.query_cover(query_len) >= min_query_cover
            && self.subject_cover(subject_len) >= min_subject_cover
    }
}

/// A match between a query and a target, containing one or more HSPs.
#[derive(Debug, Clone)]
pub struct Match {
    /// Target sequence block ID.
    pub target_block_id: u32,
    /// Target sequence original ID.
    pub target_oid: u64,
    /// List of HSPs for this query-target pair.
    pub hsps: Vec<Hsp>,
}

impl Match {
    pub fn new(target_block_id: u32, target_oid: u64) -> Self {
        Match {
            target_block_id,
            target_oid,
            hsps: Vec::new(),
        }
    }

    /// Best score among all HSPs.
    pub fn top_score(&self) -> Score {
        self.hsps.iter().map(|h| h.score).max().unwrap_or(0)
    }

    /// Best e-value among all HSPs.
    pub fn top_evalue(&self) -> f64 {
        self.hsps
            .iter()
            .map(|h| h.evalue)
            .fold(f64::MAX, f64::min)
    }

    /// Sort HSPs by score descending.
    pub fn sort_by_score(&mut self) {
        self.hsps.sort_by(|a, b| b.score.cmp(&a.score));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hsp_defaults() {
        let hsp = Hsp::new();
        assert_eq!(hsp.score, 0);
        assert_eq!(hsp.evalue, f64::MAX);
        assert!(!hsp.backtraced);
    }

    #[test]
    fn test_hsp_pident() {
        let mut hsp = Hsp::new();
        hsp.length = 100;
        hsp.identities = 80;
        assert!((hsp.pident() - 80.0).abs() < 0.001);
    }

    #[test]
    fn test_hsp_filters() {
        let mut hsp = Hsp::new();
        hsp.evalue = 1e-10;
        hsp.length = 100;
        hsp.identities = 80;
        hsp.query_range = Interval::new(0, 80);
        hsp.subject_range = Interval::new(0, 80);

        assert!(hsp.passes_filters(0.001, 50.0, 50.0, 0.0, 100, 200));
        assert!(!hsp.passes_filters(1e-20, 50.0, 0.0, 0.0, 100, 200));
    }

    #[test]
    fn test_match_top_score() {
        let mut m = Match::new(0, 0);
        let mut h1 = Hsp::new();
        h1.score = 50;
        let mut h2 = Hsp::new();
        h2.score = 100;
        m.hsps.push(h1);
        m.hsps.push(h2);
        assert_eq!(m.top_score(), 100);
    }
}
