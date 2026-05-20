use crate::basic::consts::MAX_CONTEXT;
use crate::basic::packed_transcript::{EditOperation, PackedTranscript};
use crate::basic::translate::{Frame, Strand, TranslatedPosition};
use crate::basic::value::{BlockId, Letter, OId, Score, AMINO_ACID_ALPHABET, LETTER_MASK};
use crate::stats;
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::geo::DiagonalSegmentT;
use crate::util::hsp::ApproxHsp;
use crate::util::interval::Interval;
use std::sync::Arc;

pub fn normalized_range(pos: u32, len: i32, strand: Strand) -> Interval {
    if strand == Strand::Forward {
        Interval::new(pos as i32, pos as i32 + len)
    } else {
        Interval::new(pos as i32 + 1 + len, pos as i32 + 1)
    }
}

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
    /// Mapping quality for mapping-style output.
    pub mapping_quality: u8,
    /// Number of anchors used for mapping-style output.
    pub n_anchors: i32,
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
    /// Composition-adjusted target matrix used for this HSP, when present.
    pub matrix: Option<Arc<TargetMatrix>>,
    /// Alignment transcript (CIGAR-like edit operations).
    pub transcript: PackedTranscript,
    /// Internal: swipe target index.
    pub swipe_target: i32,
    /// Internal: SWIPE bin that produced the HSP.
    pub swipe_bin: i32,
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
            mapping_quality: 0,
            n_anchors: 0,
            query_range: Interval::default(),
            subject_range: Interval::default(),
            query_source_range: Interval::default(),
            subject_source_range: Interval::default(),
            target_seq: Vec::new(),
            matrix: None,
            transcript: PackedTranscript::new(),
            swipe_target: 0,
            swipe_bin: 0,
            d_begin: 0,
            d_end: 0,
        }
    }
}

impl Hsp {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_backtraced(backtraced: bool) -> Self {
        Hsp {
            backtraced,
            ..Self::default()
        }
    }

    pub fn with_score(backtraced: bool, score: Score, swipe_target: i32) -> Self {
        Hsp {
            backtraced,
            score,
            swipe_target,
            ..Self::default()
        }
    }

    /// Matches C++ `Hsp::Hsp(const ApproxHsp&, Loc, Loc)`.
    pub fn from_approx_hsp(
        h: &ApproxHsp,
        query_len: i32,
        target_len: i32,
        score_matrix: &ScoreMatrix,
    ) -> Self {
        Hsp {
            backtraced: true,
            score: h.score,
            frame: h.frame,
            length: h.query_range.length(),
            identities: h.query_range.length(),
            positives: h.query_range.length(),
            query_source_range: h.query_range,
            query_range: h.query_range,
            subject_source_range: h.subject_range,
            subject_range: h.subject_range,
            evalue: h.evalue,
            bit_score: score_matrix.bitscore(h.score as f64),
            corrected_bit_score: score_matrix.bitscore_corrected(
                h.score,
                query_len as u32,
                target_len as u32,
            ),
            approx_id: 100.0,
            ..Self::default()
        }
    }

    /// Matches C++ `Hsp::clear()`.
    pub fn clear(&mut self) {
        self.seed_only = false;
        self.score = 0;
        self.frame = 0;
        self.length = 0;
        self.identities = 0;
        self.mismatches = 0;
        self.positives = 0;
        self.gap_openings = 0;
        self.gaps = 0;
        self.mapping_quality = 0;
        self.n_anchors = 0;
        self.transcript.clear();
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

    /// Matches C++ `Hsp::diagonal_bounds()`.
    pub fn diagonal_bounds(&self) -> (i32, i32) {
        let mut d0 = i32::MAX;
        let mut d1 = i32::MIN;
        let mut it = self.begin();
        while it.good() {
            let d = it.query_pos.translated - it.subject_pos;
            d0 = d0.min(d);
            d1 = d1.max(d);
            it.advance();
        }
        (d0, d1)
    }

    pub fn oriented_range(&self) -> Interval {
        if self.frame < 3 {
            Interval::new(
                self.query_source_range.begin,
                self.query_source_range.end - 1,
            )
        } else {
            Interval::new(
                self.query_source_range.end - 1,
                self.query_source_range.begin,
            )
        }
    }

    pub fn set_translated_query_begin(&mut self, oriented_query_begin: i32, dna_len: i32) {
        let f = if self.frame <= 2 {
            self.frame + 1
        } else {
            2 - self.frame
        };
        if f > 0 {
            self.query_range.begin = (oriented_query_begin - (f - 1)) / 3;
        } else {
            self.query_range.begin = (dna_len + f - oriented_query_begin) / 3;
        }
    }

    pub fn set_translated_query_end(&mut self, oriented_query_end: i32, dna_len: i32) {
        let f = if self.frame <= 2 {
            self.frame + 1
        } else {
            2 - self.frame
        };
        if f > 0 {
            self.query_range.end = (oriented_query_end - 2 - (f - 1)) / 3 + 1;
        } else {
            self.query_range.end = (dna_len + f - oriented_query_end) / 3 + 1;
        }
    }

    pub fn blast_query_frame(&self, query_translated: bool) -> i32 {
        if query_translated {
            if self.frame <= 2 {
                self.frame + 1
            } else {
                2 - self.frame
            }
        } else {
            0
        }
    }

    pub fn begin(&self) -> HspIterator<'_> {
        HspIterator::new(self)
    }

    pub fn cmp_evalue(a: &Hsp, b: &Hsp) -> bool {
        a.evalue < b.evalue || (a.evalue == b.evalue && a.score > b.score)
    }

    pub fn less_than_score_position(&self, rhs: &Hsp) -> bool {
        self.score > rhs.score
            || (self.score == rhs.score
                && (self.d_begin < rhs.d_begin
                    || (self.d_begin == rhs.d_begin
                        && self.query_source_range.begin < rhs.query_source_range.begin)))
    }

    pub fn partial_score(&self, hsp: &Hsp) -> i32 {
        let overlap = self.subject_range.overlap_factor(&hsp.subject_range).max(
            self.query_source_range
                .overlap_factor(&hsp.query_source_range),
        );
        ((1.0 - overlap) * self.score as f64) as i32
    }

    pub fn id_percent(&self) -> f64 {
        self.identities as f64 * 100.0 / self.length as f64
    }

    pub fn query_cover_percent(&self, query_source_len: u32) -> f64 {
        self.query_source_range.length() as f64 * 100.0 / query_source_len as f64
    }

    pub fn subject_cover_percent(&self, subject_len: u32) -> f64 {
        self.subject_range.length() as f64 * 100.0 / subject_len as f64
    }

    pub fn envelopes(&self, d: &DiagonalSegmentT, dna_len: i32, query_translated: bool) -> bool {
        self.query_source_range
            .contains(&d.query_absolute_range(dna_len, query_translated))
            || self.subject_range.contains(&d.subject_range())
    }

    pub fn is_weakly_enveloped(&self, hsp: &Hsp) -> bool {
        const OVERLAP_FACTOR: f64 = 0.9;
        self.score <= hsp.score
            && self.subject_range.overlap_factor(&hsp.subject_range) >= OVERLAP_FACTOR
            && self.query_range.overlap_factor(&hsp.query_range) >= OVERLAP_FACTOR
    }

    pub fn is_weakly_enveloped_by<'a>(
        &self,
        hsps: impl IntoIterator<Item = &'a Hsp>,
        cutoff: i32,
    ) -> bool {
        hsps.into_iter().any(|hsp| self.partial_score(hsp) < cutoff)
    }

    pub fn is_enveloped_by(&self, hsp: &Hsp, p: f64) -> bool {
        self.query_source_range
            .overlap_factor(&hsp.query_source_range)
            >= p
            || self.subject_range.overlap_factor(&hsp.subject_range) >= p
    }

    pub fn is_enveloped_by_any<'a>(&self, hsps: impl IntoIterator<Item = &'a Hsp>, p: f64) -> bool {
        hsps.into_iter().any(|hsp| self.is_enveloped_by(hsp, p))
    }

    pub fn query_range_enveloped_by(&self, hsp: &Hsp, p: f64) -> bool {
        self.query_source_range
            .overlap_factor(&hsp.query_source_range)
            >= p
    }

    pub fn query_range_enveloped_by_any<'a>(
        &self,
        hsps: impl IntoIterator<Item = &'a Hsp>,
        p: f64,
    ) -> bool {
        hsps.into_iter()
            .any(|hsp| self.query_range_enveloped_by(hsp, p))
    }

    pub fn is_identity(&self, query: &[Letter], target: &[Letter]) -> bool {
        let qb = self.query_range.begin as usize;
        let qe = self.query_range.end as usize;
        let sb = self.subject_range.begin as usize;
        let se = self.subject_range.end as usize;
        qe <= query.len()
            && se <= target.len()
            && self.query_range.length() == self.subject_range.length()
            && query[qb..qe]
                .iter()
                .zip(&target[sb..se])
                .all(|(&x, &y)| (x & LETTER_MASK) == (y & LETTER_MASK))
    }

    pub fn approx_id_percent(&self, query: &[Letter], target: &[Letter]) -> f64 {
        if self.is_identity(query, target) {
            100.0
        } else {
            stats::approx_id(
                self.score,
                self.query_range.length(),
                self.subject_range.length(),
            )
        }
    }

    pub fn min_range_len(&self, qcov: f64, tcov: f64, qlen: i32, tlen: i32) -> (i32, i32) {
        (
            (self.query_range.end
                - (self.query_range.end as f64 - qcov * qlen as f64 / 100.0).floor() as i32)
                .max(0),
            (self.subject_range.end
                - (self.subject_range.end as f64 - tcov * tlen as f64 / 100.0).floor() as i32)
                .max(0),
        )
    }

    /// Matches C++ `Hsp::push_match(Letter, Letter, bool)`.
    pub fn push_match(&mut self, q: Letter, s: Letter, positive: bool) {
        if q == s {
            self.transcript.push_with_count(EditOperation::Match, 1_u32);
            self.identities += 1;
            self.positives += 1;
        } else {
            self.transcript
                .push_with_letter(EditOperation::Substitution, s);
            self.mismatches += 1;
            if positive {
                self.positives += 1;
            }
        }
        self.length += 1;
    }

    /// Matches C++ `Hsp::push_gap(EditOperation, int, const Sequence&)`.
    ///
    /// For deletions, `subject` is the valid subject prefix ending at the C++
    /// pointer position; letters are emitted in the same `subject[-i]` order.
    pub fn push_gap(&mut self, op: EditOperation, length: i32, subject: &[Letter]) {
        self.gap_openings += 1;
        self.length += length;
        self.gaps += length;
        if op == EditOperation::Insertion {
            self.transcript
                .push_with_count(EditOperation::Insertion, length as u32);
        } else {
            let length = length as usize;
            assert!(length <= subject.len());
            for i in 0..length {
                self.transcript
                    .push_with_letter(EditOperation::Deletion, subject[subject.len() - 1 - i]);
            }
        }
    }

    /// Matches C++ `Hsp::push_back(const DiagonalSegmentT&, ...)`.
    pub fn push_back_segment(
        &mut self,
        d: &DiagonalSegmentT,
        query: &[Letter],
        subject: &[Letter],
        reversed: bool,
        score_matrix: &ScoreMatrix,
    ) {
        assert!(d.len >= 0);
        if reversed {
            for offset in 0..d.len {
                let i = (d.query_last().translated - offset) as usize;
                let j = (d.subject_last() - offset) as usize;
                let ls = subject[j];
                let lq = query[i];
                self.push_match(lq, ls, score_matrix.score(ls, lq) > 0);
            }
        } else {
            for offset in 0..d.len {
                let i = (d.i.translated + offset) as usize;
                let j = (d.j + offset) as usize;
                let ls = subject[j];
                let lq = query[i];
                self.push_match(lq, ls, score_matrix.score(ls, lq) > 0);
            }
        }
    }

    /// Matches C++ `Hsp::splice(const DiagonalSegmentT&, const DiagonalSegmentT&, ...)`.
    pub fn splice_segment(
        &mut self,
        a: &DiagonalSegmentT,
        b: &DiagonalSegmentT,
        subject: &[Letter],
        reversed: bool,
    ) {
        let mut i0 = a.query_last();
        let j0 = a.subject_last();
        let fs = i0.frame_shift(&b.i);
        if fs == 1 {
            i0.shift_forward();
            self.transcript.push(EditOperation::FrameshiftForward);
        } else if fs == -1 {
            i0.shift_back();
            self.transcript.push(EditOperation::FrameshiftReverse);
        }
        let d0 = i0.translated - j0;
        let d1 = b.diag();
        if d1 > d0 {
            self.transcript
                .push_with_count(EditOperation::Insertion, (d1 - d0) as u32);
        } else if d1 < d0 {
            let begin = (j0 + 1) as usize;
            let end = b.j as usize;
            assert!(begin <= end && end <= subject.len());
            if reversed {
                for &letter in subject[begin..end].iter().rev() {
                    self.transcript
                        .push_with_letter(EditOperation::Deletion, letter);
                }
            } else {
                for &letter in &subject[begin..end] {
                    self.transcript
                        .push_with_letter(EditOperation::Deletion, letter);
                }
            }
        }
        let shift = (d1 - d0).abs();
        if shift > 0 {
            self.length += shift;
            self.gap_openings += 1;
            self.gaps += shift;
        }
    }

    /// Matches C++ `Hsp::set_begin(int, int, Frame, int)`.
    pub fn set_begin(&mut self, i: i32, j: i32, frame: Frame, dna_len: i32) {
        self.subject_range.begin = j;
        self.query_range.begin = i;
        self.frame = frame.index();
        let absolute = TranslatedPosition::new(i, frame).absolute(dna_len, true, false);
        if frame.strand == Strand::Forward {
            self.query_source_range.begin = absolute;
        } else {
            self.query_source_range.end = absolute + 1;
        }
    }

    /// Matches C++ `Hsp::set_end(int, int, Frame, int)`.
    pub fn set_end(&mut self, i: i32, j: i32, frame: Frame, dna_len: i32) {
        self.subject_range.end = j;
        self.query_range.end = i;
        let absolute = TranslatedPosition::new(i, frame).absolute(dna_len, true, false);
        if frame.strand == Strand::Forward {
            self.query_source_range.end = absolute;
        } else {
            self.query_source_range.begin = absolute + 1;
        }
    }

    /// Matches C++ `Hsp::set_begin(const DiagonalSegmentT&, int)`.
    pub fn set_begin_segment(&mut self, d: &DiagonalSegmentT, dna_len: i32) {
        self.subject_range.begin = d.j;
        self.query_range.begin = d.i.translated;
        self.frame = d.i.frame.index();
        let absolute = d.i.absolute(dna_len, true, false);
        if d.i.frame.strand == Strand::Forward {
            self.query_source_range.begin = absolute;
        } else {
            self.query_source_range.end = absolute + 1;
        }
    }

    /// Matches C++ `Hsp::set_end(const DiagonalSegmentT&, int)`.
    pub fn set_end_segment(&mut self, d: &DiagonalSegmentT, dna_len: i32) {
        self.subject_range.end = d.subject_end();
        self.query_range.end = d.query_end().translated;
        let absolute = d.query_end().absolute(dna_len, true, false);
        if d.i.frame.strand == Strand::Forward {
            self.query_source_range.end = absolute;
        } else {
            self.query_source_range.begin = absolute + 1;
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

/// Matches C++ `Hsp::Iterator::Iterator(const Hsp&)`.
pub struct HspIterator<'a> {
    pub query_pos: TranslatedPosition,
    pub subject_pos: i32,
    ptr: &'a [crate::basic::packed_transcript::PackedOperation],
    ptr_idx: usize,
    count: u32,
}

impl<'a> HspIterator<'a> {
    pub fn new(parent: &'a Hsp) -> Self {
        let ptr = parent.transcript.data();
        let count = ptr.first().map(|op| op.count()).unwrap_or(0);
        HspIterator {
            query_pos: TranslatedPosition::new(
                parent.query_range.begin,
                Frame::from_index(parent.frame),
            ),
            subject_pos: parent.subject_range.begin,
            ptr,
            ptr_idx: 0,
            count,
        }
    }

    pub fn good(&self) -> bool {
        self.ptr_idx < self.ptr.len() && !self.ptr[self.ptr_idx].is_terminator()
    }

    pub fn op(&self) -> EditOperation {
        self.ptr[self.ptr_idx].op()
    }

    pub fn letter(&self) -> Letter {
        self.ptr[self.ptr_idx].letter()
    }

    pub fn advance(&mut self) -> &mut Self {
        match self.op() {
            EditOperation::Deletion => self.subject_pos += 1,
            EditOperation::Insertion => self.query_pos.translated += 1,
            EditOperation::Match | EditOperation::Substitution => {
                self.query_pos.translated += 1;
                self.subject_pos += 1;
            }
            EditOperation::FrameshiftForward => self.query_pos.shift_forward(),
            EditOperation::FrameshiftReverse => self.query_pos.shift_back(),
        }
        self.count -= 1;
        if self.count == 0 {
            self.ptr_idx += 1;
            self.count = self.ptr.get(self.ptr_idx).map(|op| op.count()).unwrap_or(0);
        }
        self
    }

    pub fn to_next_match(&mut self) -> bool {
        loop {
            self.advance();
            if !self.good()
                || (self.op() != EditOperation::Deletion && self.op() != EditOperation::Insertion)
            {
                return self.good();
            }
        }
    }
}

/// Matches C++ `HspContext::HspContext()`.
#[derive(Debug, Clone)]
pub struct HspContext {
    pub query: Vec<Vec<Letter>>,
    pub query_source_len: i32,
    pub query_title: String,
    pub target_title: String,
    pub query_id: BlockId,
    pub query_oid: OId,
    pub subject_oid: OId,
    pub query_len: i32,
    pub subject_len: i32,
    pub hit_num: u32,
    pub hsp_num: u32,
    pub query_self_aln_score: f64,
    pub target_self_aln_score: f64,
    pub subject_seq: Vec<Letter>,
    pub hsp: Hsp,
}

impl Default for HspContext {
    fn default() -> Self {
        HspContext {
            query: Vec::new(),
            query_source_len: 0,
            query_title: String::new(),
            target_title: String::new(),
            query_id: 0,
            query_oid: 0,
            subject_oid: 0,
            query_len: 0,
            subject_len: 0,
            hit_num: 0,
            hsp_num: 0,
            query_self_aln_score: 0.0,
            target_self_aln_score: 0.0,
            subject_seq: Vec::new(),
            hsp: Hsp::new(),
        }
    }
}

impl HspContext {
    pub fn new(
        hsp: Hsp,
        query_id: BlockId,
        query_oid: OId,
        query: Vec<Vec<Letter>>,
        query_source_len: i32,
        query_title: impl Into<String>,
        subject_oid: OId,
        subject_len: i32,
        subject_title: impl Into<String>,
        hit_num: u32,
        hsp_num: u32,
        subject_seq: Vec<Letter>,
        query_self_aln_score: f64,
        target_self_aln_score: f64,
    ) -> Self {
        HspContext {
            query,
            query_source_len,
            query_title: query_title.into(),
            target_title: subject_title.into(),
            query_id,
            query_oid,
            subject_oid,
            query_len: query_source_len,
            subject_len,
            hit_num,
            hsp_num,
            query_self_aln_score,
            target_self_aln_score,
            subject_seq,
            hsp,
        }
    }

    pub fn begin(&self) -> HspContextIterator<'_> {
        HspContextIterator::new(self)
    }

    pub fn begin_old(&self) -> crate::basic::packed_transcript::TranscriptIterator<'_> {
        self.hsp.transcript.iter()
    }

    pub fn score(&self) -> u32 {
        self.hsp.score as u32
    }

    pub fn evalue(&self) -> f64 {
        self.hsp.evalue
    }

    pub fn bit_score(&self) -> f64 {
        self.hsp.bit_score
    }

    pub fn corrected_bit_score(&self) -> f64 {
        self.hsp.corrected_bit_score
    }

    pub fn approx_id(&self) -> f64 {
        self.hsp.approx_id
    }

    pub fn frame(&self) -> u32 {
        self.hsp.frame as u32
    }

    pub fn length(&self) -> u32 {
        self.hsp.length as u32
    }

    pub fn identities(&self) -> u32 {
        self.hsp.identities as u32
    }

    pub fn mismatches(&self) -> u32 {
        self.hsp.mismatches as u32
    }

    pub fn positives(&self) -> u32 {
        self.hsp.positives as u32
    }

    pub fn gap_openings(&self) -> u32 {
        self.hsp.gap_openings as u32
    }

    pub fn gaps(&self) -> u32 {
        self.hsp.gaps as u32
    }

    pub fn approx_id_percent(&self) -> f64 {
        self.hsp
            .approx_id_percent(self.query_index(self.hsp.frame), &self.subject_seq)
    }

    pub fn qcovhsp(&self) -> f64 {
        self.query_source_range().length() as f64 * 100.0 / self.query_len as f64
    }

    pub fn scovhsp(&self) -> f64 {
        self.subject_range().length() as f64 * 100.0 / self.subject_len as f64
    }

    pub fn id_percent(&self) -> f64 {
        self.identities() as f64 * 100.0 / self.length() as f64
    }

    pub fn query_source_range(&self) -> &Interval {
        &self.hsp.query_source_range
    }

    pub fn subject_source_range(&self) -> &Interval {
        &self.hsp.subject_source_range
    }

    pub fn query_range(&self) -> &Interval {
        &self.hsp.query_range
    }

    pub fn subject_range(&self) -> &Interval {
        &self.hsp.subject_range
    }

    pub fn oriented_query_range(&self) -> Interval {
        self.hsp.oriented_range()
    }

    pub fn blast_query_frame(&self, query_translated: bool) -> i32 {
        self.hsp.blast_query_frame(query_translated)
    }

    pub fn transcript(&self) -> PackedTranscript {
        self.hsp.transcript.clone()
    }

    pub fn seed_only(&self) -> bool {
        self.hsp.seed_only
    }

    pub fn hsp(&self) -> Hsp {
        self.hsp.clone()
    }

    pub fn query_index(&self, frame: i32) -> &[Letter] {
        &self.query[frame as usize]
    }

    /// Matches C++ `HspContext::parse(bool, bool, bool, const ScoreMatrix&)`.
    pub fn parse(
        &mut self,
        need_transcript: bool,
        view_command: bool,
        query_translated: bool,
        score_matrix: &ScoreMatrix,
    ) -> Result<&mut Self, String> {
        if self.hsp.seed_only || (!need_transcript && !view_command) {
            self.hsp.query_source_range = TranslatedPosition::absolute_interval(
                TranslatedPosition::new(
                    self.hsp.query_range.begin,
                    Frame::from_index(self.hsp.frame),
                ),
                TranslatedPosition::new(
                    self.hsp.query_range.end,
                    Frame::from_index(self.hsp.frame),
                ),
                self.query_source_len,
                query_translated,
            );
            self.hsp.subject_source_range = self.hsp.subject_range;
            if !self.hsp.seed_only && !self.subject_seq.is_empty() {
                self.hsp.approx_id = self
                    .hsp
                    .approx_id_percent(self.query_index(self.hsp.frame), &self.subject_seq);
            }
            return Ok(self);
        }

        let begin_query_pos = TranslatedPosition::new(
            self.hsp.query_range.begin,
            Frame::from_index(self.hsp.frame),
        );
        let mut length = 0;
        let mut identities = 0;
        let mut mismatches = 0;
        let mut gap_openings = 0;
        let mut positives = 0;
        let mut gaps = 0;
        let mut d = 0u32;
        let mut it = self.hsp.begin();

        while it.good() {
            let query_pos = it.query_pos;
            length += 1;
            let frame = query_pos.frame.index() as usize;
            if frame >= self.query.len()
                || query_pos.translated < 0
                || query_pos.translated as usize >= self.query[frame].len()
            {
                return Err("Query sequence index out of bounds.".to_string());
            }
            match it.op() {
                EditOperation::Match => {
                    identities += 1;
                    positives += 1;
                    d = 0;
                }
                EditOperation::Substitution => {
                    mismatches += 1;
                    let q = self.query[frame][query_pos.translated as usize];
                    let s = it.letter();
                    if score_matrix.score(q, s) > 0 {
                        positives += 1;
                    }
                    d = 0;
                }
                EditOperation::Insertion | EditOperation::Deletion => {
                    if d == 0 {
                        gap_openings += 1;
                    }
                    d += 1;
                    gaps += 1;
                }
                _ => {}
            }
            it.advance();
        }

        let query_pos = it.query_pos;
        let subject_pos = it.subject_pos;
        self.hsp.length = length;
        self.hsp.identities = identities;
        self.hsp.mismatches = mismatches;
        self.hsp.gap_openings = gap_openings;
        self.hsp.positives = positives;
        self.hsp.gaps = gaps;
        self.hsp.query_range.end = query_pos.translated;
        self.hsp.subject_range.end = subject_pos;
        self.hsp.subject_source_range = self.hsp.subject_range;
        self.hsp.query_source_range = TranslatedPosition::absolute_interval(
            begin_query_pos,
            query_pos,
            self.query_source_len,
            query_translated,
        );
        if !self.subject_seq.is_empty() {
            self.hsp.approx_id = self
                .hsp
                .approx_id_percent(self.query_index(self.hsp.frame), &self.subject_seq);
        }

        Ok(self)
    }
}

impl PartialEq for HspContext {
    fn eq(&self, other: &Self) -> bool {
        self.query_oid == other.query_oid
    }
}

impl Eq for HspContext {}

impl PartialOrd for HspContext {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HspContext {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.query_oid.cmp(&other.query_oid)
    }
}

/// Matches C++ `HspContext::Iterator::Iterator(const HspContext&)`.
pub struct HspContextIterator<'a> {
    parent: &'a HspContext,
    inner: HspIterator<'a>,
}

impl<'a> HspContextIterator<'a> {
    pub fn new(parent: &'a HspContext) -> Self {
        HspContextIterator {
            parent,
            inner: parent.hsp.begin(),
        }
    }

    pub fn good(&self) -> bool {
        self.inner.good()
    }

    pub fn advance(&mut self) -> &mut Self {
        self.inner.advance();
        self
    }

    pub fn op(&self) -> EditOperation {
        self.inner.op()
    }

    pub fn query_pos(&self) -> TranslatedPosition {
        self.inner.query_pos
    }

    pub fn subject_pos(&self) -> i32 {
        self.inner.subject_pos
    }

    pub fn query(&self) -> Letter {
        let frame = self.inner.query_pos.frame.index() as usize;
        self.parent.query[frame][self.inner.query_pos.translated as usize]
    }

    pub fn subject(&self) -> Letter {
        match self.op() {
            EditOperation::Substitution | EditOperation::Deletion => self.inner.letter(),
            _ => self.query(),
        }
    }

    pub fn query_char(&self) -> char {
        match self.op() {
            EditOperation::Deletion => '-',
            EditOperation::FrameshiftForward => '\\',
            EditOperation::FrameshiftReverse => '/',
            _ => AMINO_ACID_ALPHABET[(self.query() & LETTER_MASK) as usize] as char,
        }
    }

    pub fn subject_char(&self) -> char {
        match self.op() {
            EditOperation::Insertion
            | EditOperation::FrameshiftForward
            | EditOperation::FrameshiftReverse => '-',
            _ => AMINO_ACID_ALPHABET[(self.subject() & LETTER_MASK) as usize] as char,
        }
    }

    pub fn midline_char(&self, score: i32) -> char {
        match self.op() {
            EditOperation::Match => {
                AMINO_ACID_ALPHABET[(self.query() & LETTER_MASK) as usize] as char
            }
            EditOperation::Substitution => {
                if score > 0 {
                    '+'
                } else {
                    ' '
                }
            }
            _ => ' ',
        }
    }
}

/// A match between a query and a target, containing one or more HSPs.
#[derive(Debug, Clone)]
pub struct Match {
    /// Target sequence block ID.
    pub target_block_id: u32,
    /// Target sequence original ID.
    pub target_oid: u64,
    pub seq: Vec<Letter>,
    pub matrix: Option<Arc<TargetMatrix>>,
    pub filter_score: Score,
    pub filter_evalue: f64,
    pub ungapped_score: Score,
    /// List of HSPs for this query-target pair.
    pub hsps: Vec<Hsp>,
}

impl Match {
    pub fn new(target_block_id: u32, target_oid: u64) -> Self {
        Match {
            target_block_id,
            target_oid,
            seq: Vec::new(),
            matrix: None,
            filter_score: 0,
            filter_evalue: f64::MAX,
            ungapped_score: 0,
            hsps: Vec::new(),
        }
    }

    /// Matches C++ `Extension::Match::Match(BlockId, Sequence, TargetMatrix, Score, Score, double)`.
    pub fn new_extension(
        target_block_id: u32,
        seq: &[Letter],
        matrix: Option<Arc<TargetMatrix>>,
        ungapped_score: Score,
        filter_score: Score,
        filter_evalue: f64,
    ) -> Self {
        Match {
            target_block_id,
            target_oid: 0,
            seq: seq.to_vec(),
            matrix,
            filter_score,
            filter_evalue,
            ungapped_score,
            hsps: Vec::new(),
        }
    }

    /// Matches C++ `Match::add_hit(list<Hsp>&, iterator)`.
    pub fn add_hit_from_vec(&mut self, hsps: &mut Vec<Hsp>, index: usize) {
        let hsp = hsps.remove(index);
        self.hsps.push(hsp);
        let last = self.hsps.last().unwrap();
        if last.score > self.filter_score {
            self.filter_evalue = last.evalue;
            self.filter_score = last.score;
        }
    }

    /// Matches C++ `Match::Match(BlockId, Sequence, TargetMatrix, array<list<Hsp>>, int)`.
    pub fn from_target_hsps(
        target_block_id: u32,
        seq: &[Letter],
        matrix: Option<Arc<TargetMatrix>>,
        hsps: &mut [Vec<Hsp>; MAX_CONTEXT as usize],
        ungapped_score: Score,
        query_contexts: usize,
        max_hsps: u32,
    ) -> Self {
        if max_hsps != 1 {
            panic!("Match::Match max_hsps != 1.");
        }
        let mut m = Match::new_extension(target_block_id, seq, matrix, ungapped_score, 0, f64::MAX);
        for frame_hsps in hsps.iter_mut().take(query_contexts) {
            m.hsps.append(frame_hsps);
        }
        if m.hsps.is_empty() {
            panic!("Match::Match hsp.empty()");
        }
        m.hsps.sort_by(|a, b| {
            if a.less_than_score_position(b) {
                std::cmp::Ordering::Less
            } else if b.less_than_score_position(a) {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });
        m.hsps.truncate(1);
        m.filter_evalue = m.hsps[0].evalue;
        m.filter_score = m.hsps[0].score;
        m
    }

    /// Matches C++ `Match::self_match`.
    pub fn self_match(query_id: BlockId, query_seq: &[Letter]) -> Self {
        let mut m = Match::new_extension(query_id, query_seq, None, 0, Score::MAX, 0.0);
        let mut hsp = Hsp::new();
        hsp.evalue = 0.0;
        hsp.score = Score::MAX;
        hsp.bit_score = f64::MAX;
        hsp.query_range = Interval::new(0, query_seq.len() as i32);
        hsp.query_source_range = Interval::new(0, query_seq.len() as i32);
        hsp.subject_range = Interval::new(0, query_seq.len() as i32);
        m.hsps.push(hsp);
        m
    }

    /// Best score among all HSPs.
    pub fn top_score(&self) -> Score {
        if self.hsps.is_empty() {
            self.filter_score
        } else {
            self.hsps.iter().map(|h| h.score).max().unwrap_or(0)
        }
    }

    /// Best e-value among all HSPs.
    pub fn top_evalue(&self) -> f64 {
        if self.hsps.is_empty() {
            self.filter_evalue
        } else {
            self.hsps.iter().map(|h| h.evalue).fold(f64::MAX, f64::min)
        }
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

        let hsp = Hsp::with_backtraced(true);
        assert!(hsp.backtraced);
        assert_eq!(hsp.score, 0);

        let hsp = Hsp::with_score(true, 42, 7);
        assert!(hsp.backtraced);
        assert_eq!(hsp.score, 42);
        assert_eq!(hsp.swipe_target, 7);
    }

    #[test]
    fn test_normalized_range() {
        assert_eq!(
            normalized_range(10, 5, Strand::Forward),
            Interval::new(10, 15)
        );
        assert_eq!(
            normalized_range(10, 5, Strand::Reverse),
            Interval::new(16, 11)
        );
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
    fn test_hsp_envelope_and_identity_methods() {
        let mut h1 = Hsp::new();
        h1.score = 80;
        h1.query_range = Interval::new(10, 30);
        h1.query_source_range = Interval::new(10, 30);
        h1.subject_range = Interval::new(20, 40);

        let mut h2 = Hsp::new();
        h2.score = 100;
        h2.query_range = Interval::new(9, 31);
        h2.query_source_range = Interval::new(9, 31);
        h2.subject_range = Interval::new(19, 41);

        assert!(h1.is_weakly_enveloped(&h2));
        assert!(h1.is_enveloped_by(&h2, 0.9));
        assert!(h1.query_range_enveloped_by(&h2, 0.9));
        assert!(h1.is_weakly_enveloped_by([&h2], 20));
        assert_eq!(h1.min_range_len(50.0, 25.0, 100, 200), (50, 50));

        let query = vec![0, 1, 2, 3];
        let target = vec![0, 1, 2, 3];
        let mut ident = Hsp::new();
        ident.query_range = Interval::new(0, 4);
        ident.subject_range = Interval::new(0, 4);
        assert!(ident.is_identity(&query, &target));
        assert_eq!(ident.approx_id_percent(&query, &target), 100.0);
    }

    #[test]
    fn test_hsp_diagonal_bounds_and_frame_accessors() {
        let mut hsp = Hsp::new();
        hsp.query_range = Interval::new(10, 0);
        hsp.subject_range = Interval::new(5, 0);
        hsp.transcript.push_with_count(EditOperation::Match, 2);
        hsp.transcript.push_with_count(EditOperation::Insertion, 2);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 8);

        assert_eq!(hsp.diagonal_bounds(), (5, 7));

        hsp.frame = 1;
        hsp.query_source_range = Interval::new(7, 16);
        assert_eq!(hsp.oriented_range(), Interval::new(7, 15));
        assert_eq!(hsp.blast_query_frame(true), 2);
        assert_eq!(hsp.blast_query_frame(false), 0);
        hsp.set_translated_query_begin(7, 100);
        hsp.set_translated_query_end(16, 100);
        assert_eq!(hsp.query_range, Interval::new(2, 5));

        hsp.frame = 4;
        hsp.query_source_range = Interval::new(86, 95);
        assert_eq!(hsp.oriented_range(), Interval::new(94, 86));
        assert_eq!(hsp.blast_query_frame(true), -2);
        hsp.set_translated_query_begin(94, 100);
        hsp.set_translated_query_end(86, 100);
        assert_eq!(hsp.query_range, Interval::new(1, 5));
    }

    #[test]
    fn test_hsp_iterator_positions_and_to_next_match() {
        let mut hsp = Hsp::new();
        hsp.query_range = Interval::new(10, 0);
        hsp.subject_range = Interval::new(20, 0);
        hsp.transcript.push_with_count(EditOperation::Insertion, 2);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 7);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 8);
        hsp.transcript.push(EditOperation::FrameshiftForward);
        hsp.transcript.push_with_count(EditOperation::Match, 1);

        let mut it = hsp.begin();
        assert!(it.good());
        assert_eq!(it.op(), EditOperation::Insertion);
        assert_eq!(it.query_pos.translated, 10);
        assert_eq!(it.subject_pos, 20);

        it.advance();
        assert_eq!(it.op(), EditOperation::Insertion);
        assert_eq!(it.query_pos.translated, 11);
        assert_eq!(it.subject_pos, 20);

        assert!(it.to_next_match());
        assert_eq!(it.op(), EditOperation::Substitution);
        assert_eq!(it.query_pos.translated, 12);
        assert_eq!(it.subject_pos, 21);

        it.advance();
        assert_eq!(it.op(), EditOperation::FrameshiftForward);
        let before_shift = it.query_pos;
        it.advance();
        assert_eq!(it.op(), EditOperation::Match);
        assert_ne!(it.query_pos.frame, before_shift.frame);
        assert_eq!(it.subject_pos, 22);

        it.advance();
        assert!(!it.good());
    }

    #[test]
    fn test_hsp_context_accessors_and_ordering() {
        let mut hsp = Hsp::new();
        hsp.score = 80;
        hsp.evalue = 1e-20;
        hsp.bit_score = 40.0;
        hsp.corrected_bit_score = 38.0;
        hsp.approx_id = 91.0;
        hsp.frame = 0;
        hsp.length = 20;
        hsp.identities = 15;
        hsp.mismatches = 3;
        hsp.positives = 16;
        hsp.gap_openings = 1;
        hsp.gaps = 2;
        hsp.query_range = Interval::new(0, 20);
        hsp.query_source_range = Interval::new(10, 30);
        hsp.subject_range = Interval::new(50, 70);
        hsp.subject_source_range = hsp.subject_range;

        let ctx = HspContext::new(
            hsp.clone(),
            3,
            11,
            vec![vec![0; 20]],
            100,
            "query",
            22,
            200,
            "subject",
            5,
            6,
            vec![0; 20],
            44.0,
            55.0,
        );

        assert_eq!(ctx.score(), 80);
        assert_eq!(ctx.evalue(), 1e-20);
        assert_eq!(ctx.bit_score(), 40.0);
        assert_eq!(ctx.corrected_bit_score(), 38.0);
        assert_eq!(ctx.approx_id(), 91.0);
        assert_eq!(ctx.frame(), 0);
        assert_eq!(ctx.length(), 20);
        assert_eq!(ctx.identities(), 15);
        assert_eq!(ctx.mismatches(), 3);
        assert_eq!(ctx.positives(), 16);
        assert_eq!(ctx.gap_openings(), 1);
        assert_eq!(ctx.gaps(), 2);
        assert_eq!(ctx.qcovhsp(), 20.0);
        assert_eq!(ctx.scovhsp(), 10.0);
        assert_eq!(ctx.id_percent(), 75.0);
        assert_eq!(*ctx.query_source_range(), Interval::new(10, 30));
        assert_eq!(*ctx.subject_source_range(), Interval::new(50, 70));
        assert_eq!(*ctx.query_range(), Interval::new(0, 20));
        assert_eq!(*ctx.subject_range(), Interval::new(50, 70));
        assert_eq!(ctx.oriented_query_range(), Interval::new(10, 29));
        assert!(!ctx.seed_only());
        assert_eq!(ctx.hsp().score, hsp.score);

        let mut later = ctx.clone();
        later.query_oid = 12;
        assert!(ctx < later);
    }

    #[test]
    fn test_hsp_context_iterator_letters_and_chars() {
        let mut hsp = Hsp::new();
        hsp.query_range = Interval::new(0, 0);
        hsp.subject_range = Interval::new(0, 0);
        hsp.transcript.push_with_count(EditOperation::Match, 1);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 2);
        hsp.transcript.push_with_count(EditOperation::Insertion, 1);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 3);
        hsp.transcript.push(EditOperation::FrameshiftForward);

        let ctx = HspContext::new(
            hsp,
            0,
            0,
            vec![vec![0, 1, 4]],
            9,
            "",
            0,
            10,
            "",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );

        let mut it = ctx.begin();
        assert!(it.good());
        assert_eq!(it.op(), EditOperation::Match);
        assert_eq!(it.query(), 0);
        assert_eq!(it.subject(), 0);
        assert_eq!(it.query_char(), 'A');
        assert_eq!(it.subject_char(), 'A');
        assert_eq!(it.midline_char(0), 'A');

        it.advance();
        assert_eq!(it.op(), EditOperation::Substitution);
        assert_eq!(it.query(), 1);
        assert_eq!(it.subject(), 2);
        assert_eq!(it.query_char(), 'R');
        assert_eq!(it.subject_char(), 'N');
        assert_eq!(it.midline_char(1), '+');
        assert_eq!(it.midline_char(0), ' ');

        it.advance();
        assert_eq!(it.op(), EditOperation::Insertion);
        assert_eq!(it.query_char(), 'C');
        assert_eq!(it.subject_char(), '-');

        it.advance();
        assert_eq!(it.op(), EditOperation::Deletion);
        assert_eq!(it.query_char(), '-');
        assert_eq!(it.subject_char(), 'D');

        it.advance();
        assert_eq!(it.op(), EditOperation::FrameshiftForward);
        assert_eq!(it.query_char(), '\\');
        assert_eq!(it.subject_char(), '-');
    }

    #[test]
    fn test_hsp_context_parse_transcript() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1000).unwrap();
        let mut hsp = Hsp::new();
        hsp.frame = 0;
        hsp.query_range = Interval::new(0, 0);
        hsp.subject_range = Interval::new(0, 0);
        hsp.transcript.push_with_count(EditOperation::Match, 1);
        hsp.transcript
            .push_with_letter(EditOperation::Substitution, 2);
        hsp.transcript.push_with_count(EditOperation::Insertion, 1);
        hsp.transcript.push_with_letter(EditOperation::Deletion, 3);

        let mut ctx = HspContext::new(
            hsp,
            0,
            0,
            vec![vec![0, 1, 2, 3]],
            3,
            "",
            0,
            10,
            "",
            0,
            0,
            vec![0, 2, 4],
            0.0,
            0.0,
        );

        ctx.parse(true, false, false, &score_matrix).unwrap();
        assert_eq!(ctx.length(), 4);
        assert_eq!(ctx.identities(), 1);
        assert_eq!(ctx.mismatches(), 1);
        assert_eq!(ctx.gap_openings(), 1);
        assert_eq!(ctx.gaps(), 2);
        assert_eq!(ctx.query_range().end, 3);
        assert_eq!(ctx.subject_range().end, 3);
        assert_eq!(*ctx.subject_source_range(), Interval::new(0, 3));
        assert_eq!(*ctx.query_source_range(), Interval::new(0, 3));
    }

    #[test]
    fn test_hsp_context_parse_without_transcript_and_bounds_error() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1000).unwrap();
        let mut hsp = Hsp::new();
        hsp.frame = 0;
        hsp.query_range = Interval::new(1, 4);
        hsp.subject_range = Interval::new(2, 7);

        let mut ctx = HspContext::new(
            hsp,
            0,
            0,
            vec![vec![0, 1, 2, 3]],
            10,
            "",
            0,
            10,
            "",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );
        ctx.parse(false, false, false, &score_matrix).unwrap();
        assert_eq!(*ctx.query_source_range(), Interval::new(1, 4));
        assert_eq!(*ctx.subject_source_range(), Interval::new(2, 7));

        ctx.hsp.query_range = Interval::new(4, 0);
        ctx.hsp.transcript.clear();
        ctx.hsp.transcript.push_with_count(EditOperation::Match, 1);
        assert_eq!(
            ctx.parse(true, false, false, &score_matrix).unwrap_err(),
            "Query sequence index out of bounds."
        );
    }

    #[test]
    fn test_hsp_comparators() {
        let mut a = Hsp::new();
        a.evalue = 1e-20;
        a.score = 80;
        a.d_begin = 2;
        a.query_source_range = Interval::new(5, 10);

        let mut b = Hsp::new();
        b.evalue = 1e-10;
        b.score = 70;
        b.d_begin = 3;
        b.query_source_range = Interval::new(3, 8);

        assert!(Hsp::cmp_evalue(&a, &b));
        assert!(a.less_than_score_position(&b));

        b.evalue = a.evalue;
        b.score = 90;
        assert!(Hsp::cmp_evalue(&b, &a));
    }

    #[test]
    fn test_hsp_push_back_segment_forward_and_reversed() {
        let score_matrix = ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1000).unwrap();
        let frame = Frame::new(Strand::Forward, 0);
        let d = DiagonalSegmentT::new(TranslatedPosition::new(1, frame), 2, 3, 0);
        let query = vec![0, 1, 2, 3, 4];
        let subject = vec![0, 0, 1, 7, 3, 4];

        let mut forward = Hsp::new();
        forward.push_back_segment(&d, &query, &subject, false, &score_matrix);
        assert_eq!(forward.length, 3);
        assert_eq!(forward.identities, 2);
        assert_eq!(forward.mismatches, 1);
        assert_eq!(forward.positives, 2);
        let ops: Vec<_> = forward.transcript.iter().collect();
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0].op, EditOperation::Match);
        assert_eq!(ops[0].count, 1);
        assert_eq!(ops[1].op, EditOperation::Substitution);
        assert_eq!(ops[1].letter, 7);
        assert_eq!(ops[2].op, EditOperation::Match);
        assert_eq!(ops[2].count, 1);

        let mut reversed = Hsp::new();
        reversed.push_back_segment(&d, &query, &subject, true, &score_matrix);
        assert_eq!(reversed.length, 3);
        assert_eq!(reversed.identities, 2);
        assert_eq!(reversed.mismatches, 1);
        let ops: Vec<_> = reversed.transcript.iter().collect();
        assert_eq!(ops[0].op, EditOperation::Match);
        assert_eq!(ops[0].count, 1);
        assert_eq!(ops[1].op, EditOperation::Substitution);
        assert_eq!(ops[1].letter, 7);
        assert_eq!(ops[2].op, EditOperation::Match);
        assert_eq!(ops[2].count, 1);
    }

    #[test]
    fn test_hsp_splice_segment_insertion_and_deletion() {
        let frame = Frame::new(Strand::Forward, 0);
        let a = DiagonalSegmentT::new(TranslatedPosition::new(0, frame), 0, 2, 0);

        let mut insertion = Hsp::new();
        let b = DiagonalSegmentT::new(TranslatedPosition::new(4, frame), 2, 2, 0);
        insertion.splice_segment(&a, &b, &[0, 1, 2, 3, 4], false);
        assert_eq!(insertion.length, 2);
        assert_eq!(insertion.gap_openings, 1);
        assert_eq!(insertion.gaps, 2);
        let ops: Vec<_> = insertion.transcript.iter().collect();
        assert_eq!(ops.len(), 1);
        assert_eq!(ops[0].op, EditOperation::Insertion);
        assert_eq!(ops[0].count, 2);

        let mut deletion = Hsp::new();
        let b = DiagonalSegmentT::new(TranslatedPosition::new(2, frame), 4, 2, 0);
        deletion.splice_segment(&a, &b, &[0, 1, 8, 9, 4], false);
        assert_eq!(deletion.length, 2);
        assert_eq!(deletion.gap_openings, 1);
        assert_eq!(deletion.gaps, 2);
        let ops: Vec<_> = deletion.transcript.iter().collect();
        assert_eq!(ops.len(), 2);
        assert_eq!(ops[0].op, EditOperation::Deletion);
        assert_eq!(ops[0].letter, 8);
        assert_eq!(ops[1].op, EditOperation::Deletion);
        assert_eq!(ops[1].letter, 9);

        let mut reversed_deletion = Hsp::new();
        reversed_deletion.splice_segment(&a, &b, &[0, 1, 8, 9, 4], true);
        let ops: Vec<_> = reversed_deletion.transcript.iter().collect();
        assert_eq!(ops[0].letter, 9);
        assert_eq!(ops[1].letter, 8);
    }

    #[test]
    fn test_hsp_splice_segment_frameshift() {
        let a_frame = Frame::new(Strand::Forward, 0);
        let b_frame = Frame::new(Strand::Forward, 1);
        let a = DiagonalSegmentT::new(TranslatedPosition::new(0, a_frame), 0, 1, 0);
        let b = DiagonalSegmentT::new(TranslatedPosition::new(0, b_frame), 0, 1, 0);

        let mut hsp = Hsp::new();
        hsp.splice_segment(&a, &b, &[0, 1, 2], false);
        let ops: Vec<_> = hsp.transcript.iter().collect();
        assert_eq!(ops.len(), 1);
        assert_eq!(ops[0].op, EditOperation::FrameshiftForward);
        assert_eq!(hsp.length, 0);
        assert_eq!(hsp.gaps, 0);
    }

    #[test]
    fn test_hsp_envelopes_diagonal_segment_t() {
        let frame = Frame::new(Strand::Forward, 1);
        let mut hsp = Hsp::new();
        hsp.query_source_range = Interval::new(5, 30);
        hsp.subject_range = Interval::new(40, 60);

        let query_enveloped = DiagonalSegmentT::new(TranslatedPosition::new(2, frame), 70, 4, 80);
        assert!(hsp.envelopes(&query_enveloped, 100, true));

        let subject_enveloped =
            DiagonalSegmentT::new(TranslatedPosition::new(20, frame), 45, 5, 80);
        assert!(hsp.envelopes(&subject_enveloped, 100, true));

        let outside = DiagonalSegmentT::new(TranslatedPosition::new(20, frame), 70, 5, 80);
        assert!(!hsp.envelopes(&outside, 100, true));
    }

    #[test]
    fn test_hsp_set_begin_and_end_forward_frame() {
        let frame = Frame::new(Strand::Forward, 1);
        let mut hsp = Hsp::new();

        hsp.set_begin(2, 20, frame, 100);
        hsp.set_end(5, 30, frame, 100);

        assert_eq!(hsp.frame, 1);
        assert_eq!(hsp.query_range, Interval::new(2, 5));
        assert_eq!(hsp.subject_range, Interval::new(20, 30));
        assert_eq!(hsp.query_source_range, Interval::new(7, 16));
    }

    #[test]
    fn test_hsp_set_begin_and_end_reverse_frame() {
        let frame = Frame::new(Strand::Reverse, 2);
        let mut hsp = Hsp::new();

        hsp.set_begin(1, 20, frame, 100);
        hsp.set_end(4, 30, frame, 100);

        assert_eq!(hsp.frame, 5);
        assert_eq!(hsp.query_range, Interval::new(1, 4));
        assert_eq!(hsp.subject_range, Interval::new(20, 30));
        assert_eq!(hsp.query_source_range, Interval::new(86, 95));
    }

    #[test]
    fn test_hsp_set_begin_and_end_segment_overloads() {
        let frame = Frame::new(Strand::Forward, 1);
        let d = DiagonalSegmentT::new(TranslatedPosition::new(2, frame), 20, 3, 30);
        let mut hsp = Hsp::new();

        hsp.set_begin_segment(&d, 100);
        hsp.set_end_segment(&d, 100);

        assert_eq!(hsp.frame, 1);
        assert_eq!(hsp.query_range, Interval::new(2, 5));
        assert_eq!(hsp.subject_range, Interval::new(20, 23));
        assert_eq!(hsp.query_source_range, Interval::new(7, 16));
    }

    #[test]
    fn test_hsp_push_match_and_gap_methods() {
        let mut hsp = Hsp::new();

        hsp.push_match(4, 4, false);
        hsp.push_match(4, 7, true);
        hsp.push_gap(EditOperation::Insertion, 3, &[]);
        hsp.push_gap(EditOperation::Deletion, 2, &[10, 11, 12]);

        assert_eq!(hsp.length, 7);
        assert_eq!(hsp.identities, 1);
        assert_eq!(hsp.mismatches, 1);
        assert_eq!(hsp.positives, 2);
        assert_eq!(hsp.gap_openings, 2);
        assert_eq!(hsp.gaps, 5);

        let ops: Vec<_> = hsp.transcript.iter().collect();
        assert_eq!(ops.len(), 5);
        assert_eq!(ops[0].op, EditOperation::Match);
        assert_eq!(ops[0].count, 1);
        assert_eq!(ops[1].op, EditOperation::Substitution);
        assert_eq!(ops[1].letter, 7);
        assert_eq!(ops[2].op, EditOperation::Insertion);
        assert_eq!(ops[2].count, 3);
        assert_eq!(ops[3].op, EditOperation::Deletion);
        assert_eq!(ops[3].letter, 12);
        assert_eq!(ops[4].op, EditOperation::Deletion);
        assert_eq!(ops[4].letter, 11);
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
