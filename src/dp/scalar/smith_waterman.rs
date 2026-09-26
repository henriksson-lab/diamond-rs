//! Scalar Smith–Waterman facade for `dp/scalar/smith_waterman.cpp` and the
//! accompanying `scalar.h`, `traceback.h`, and `needleman_wunsch.h` helpers.
//!
//! The established `dp::smith_waterman` module remains the compatibility API;
//! this mirrored module re-exports it and supplies the original HSP-oriented
//! entry point plus the header-only semiglobal routine.

use crate::align::hsp::Hsp;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{Letter, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::interval::Interval;

pub use crate::dp::smith_waterman::{
    needleman_wunsch, smith_waterman, smith_waterman_cbs, SwResult,
};

/// C++ `smith_waterman(Sequence, Sequence, Hsp&)` with the score matrix made
/// explicit instead of read from process-global state.
pub fn smith_waterman_into_hsp(
    query: &[Letter],
    subject: &[Letter],
    out: &mut Hsp,
    score_matrix: &ScoreMatrix,
) {
    let result = smith_waterman(query, subject, score_matrix);

    out.clear();
    out.score = result.score;
    out.length = result.length;
    out.identities = result.identities;
    out.query_range = Interval::new(result.query_begin, result.query_end);
    out.subject_range = Interval::new(result.subject_begin, result.subject_end);
    out.query_source_range = out.query_range;

    let mut query_pos = result.query_begin as usize;
    let mut subject_pos = result.subject_begin as usize;
    for (op, count) in result.operations {
        match op {
            EditOperation::Match => {
                out.transcript
                    .push_with_count(EditOperation::Match, count as u32);
                query_pos += count as usize;
                subject_pos += count as usize;
            }
            EditOperation::Substitution => {
                for _ in 0..count {
                    out.transcript.push_with_letter(
                        EditOperation::Substitution,
                        subject[subject_pos] & LETTER_MASK,
                    );
                    query_pos += 1;
                    subject_pos += 1;
                }
            }
            EditOperation::Deletion => {
                for _ in 0..count {
                    out.transcript.push_with_letter(
                        EditOperation::Deletion,
                        subject[subject_pos] & LETTER_MASK,
                    );
                    subject_pos += 1;
                }
            }
            EditOperation::Insertion => {
                out.transcript
                    .push_with_count(EditOperation::Insertion, count as u32);
                query_pos += count as usize;
            }
            EditOperation::FrameshiftForward | EditOperation::FrameshiftReverse => {
                unreachable!("scalar protein Smith-Waterman cannot emit frameshifts")
            }
        }
    }
    debug_assert_eq!(query_pos, result.query_end as usize);
    debug_assert_eq!(subject_pos, result.subject_end as usize);
    out.transcript.push_terminator();
}

/// Header-only C++ `nw_semiglobal`: align the complete query to the best
/// target prefix/end position using affine gaps.
pub fn nw_semiglobal(query: &[Letter], target: &[Letter], score_matrix: &ScoreMatrix) -> i32 {
    let rows = query.len() + 1;
    let cols = target.len() + 1;
    let index = |i: usize, j: usize| j * rows + i;
    let mut matrix = vec![0i32; rows * cols];
    let mut hgap = vec![0i32; rows];

    for i in 1..=query.len() {
        hgap[i] = -score_matrix.gap_open() - i as i32 * score_matrix.gap_extend();
        matrix[index(i, 0)] = hgap[i] - score_matrix.gap_open() - score_matrix.gap_extend();
    }

    let mut best = i32::MIN;
    for j in 1..=target.len() {
        let mut vgap = -score_matrix.gap_open() - score_matrix.gap_extend();
        for i in 1..=query.len() {
            let mut score = matrix[index(i - 1, j - 1)]
                + crate::dp::smith_waterman::score_letters(
                    query[i - 1],
                    target[j - 1],
                    score_matrix,
                );
            score = score.max(vgap).max(hgap[i]);
            matrix[index(i, j)] = score;

            vgap -= score_matrix.gap_extend();
            hgap[i] -= score_matrix.gap_extend();
            let open = score - score_matrix.gap_open() - score_matrix.gap_extend();
            vgap = vgap.max(open);
            hgap[i] = hgap[i].max(open);
        }
        best = best.max(matrix[index(query.len(), j)]);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;

    fn matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn hsp_entry_point_writes_exact_ranges_and_terminated_transcript() {
        let sm = matrix();
        let query = [0, 1, 2, 3];
        let subject = [0, 1, 4, 3];
        let mut out = Hsp::default();
        smith_waterman_into_hsp(&query, &subject, &mut out, &sm);

        assert!(out.score > 0);
        assert_eq!(out.query_source_range, out.query_range);
        assert_eq!(out.length, out.query_range.length());
        assert!(out.transcript.data().last().unwrap().is_terminator());
        let substitutions: Vec<_> = out
            .transcript
            .data()
            .iter()
            .filter(|op| op.op() == EditOperation::Substitution)
            .collect();
        assert_eq!(substitutions.len(), 1);
        assert_eq!(substitutions[0].letter(), 4);
    }

    #[test]
    fn empty_local_alignment_still_has_cpp_terminator() {
        let sm = matrix();
        let mut out = Hsp::default();
        smith_waterman_into_hsp(&[0], &[17], &mut out, &sm);
        assert_eq!(out.score, 0);
        assert_eq!(out.transcript.data().len(), 1);
        assert!(out.transcript.data()[0].is_terminator());
    }

    #[test]
    fn scalar_sequence_mask_matches_cpp_sequence_operator() {
        let sm = matrix();
        let plain = smith_waterman(&[0, 0], &[0, 0], &sm);
        let masked = smith_waterman(&[0 | SEED_MASK, 0 | SEED_MASK], &[0, 0], &sm);
        assert_eq!(masked.score, plain.score);
        assert_eq!(masked.identities, plain.identities);
    }

    #[test]
    fn semiglobal_selects_best_target_end() {
        let sm = matrix();
        let query = [0, 1, 2];
        let exact = nw_semiglobal(&query, &[17, 0, 1, 2, 17], &sm);
        let global = needleman_wunsch(&query, &[17, 0, 1, 2, 17], &sm).0;
        assert!(exact > global);
        assert_eq!(exact, sm.score(0, 0) + sm.score(1, 1) + sm.score(2, 2));
    }
}
