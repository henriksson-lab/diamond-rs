//! Smith-Waterman diagnostics for the chaining graph.
//!
//! This mirrors `diamond/src/chaining/smith_waterman.cpp`. The C++ routines
//! write to `std::cout`; Rust accepts an explicit writer and propagates I/O
//! errors while preserving the emitted bytes.

use std::io::{self, Write};

use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{letter_mask, Letter};
use crate::dp::ungapped::{score_range, DiagonalSegment};
use crate::stats::score_matrix::ScoreMatrix;

use super::DiagGraph;

pub fn print_diag<W: Write>(
    query_begin: i32,
    subject_begin: i32,
    len: i32,
    score: i32,
    diags: &DiagGraph,
    query: &[Letter],
    subject: &[Letter],
    score_matrix: &ScoreMatrix,
    out: &mut W,
) -> io::Result<()> {
    let segment = DiagonalSegment::new(query_begin, subject_begin, len, 0);
    let mut count = 0usize;
    for (idx, diagonal) in diags.nodes.iter().enumerate() {
        if diagonal.intersect(&segment).len <= 0 || diagonal.score == 0 {
            continue;
        }
        let diff = score_range(
            query,
            subject,
            diagonal.query_end() as usize,
            diagonal.subject_end() as usize,
            (subject_begin + len) as usize,
            score_matrix,
        );
        if count > 0 {
            write!(out, "(")?;
        }
        let prefix_score = score
            + score_range(
                query,
                subject,
                (query_begin + len) as usize,
                (subject_begin + len) as usize,
                diagonal.subject_end() as usize,
                score_matrix,
            )
            - diff.min(0);
        let graph_prefix_score = diags.prefix_score(idx, subject_begin + len).0;
        write!(
            out,
            "Diag n={} i={} j={} len={} prefix_score={} prefix_score2={}",
            idx, query_begin, subject_begin, len, prefix_score, graph_prefix_score
        )?;
        if count > 0 {
            write!(out, ")")?;
        }
        writeln!(out)?;
        count += 1;
    }
    if count == 0 {
        writeln!(
            out,
            "Diag n=x i={} j={} len={} prefix_score={}",
            query_begin, subject_begin, len, score
        )?;
    }
    Ok(())
}

pub fn smith_waterman<W: Write>(
    query: &[Letter],
    subject: &[Letter],
    diags: &DiagGraph,
    score_matrix: &ScoreMatrix,
    out: &mut W,
) -> io::Result<()> {
    let hsp = crate::dp::smith_waterman::smith_waterman(query, subject, score_matrix);
    let mut query_pos = hsp.query_begin;
    let mut subject_pos = hsp.subject_begin;
    let mut query_begin = -1;
    let mut subject_begin = -1;
    let mut len = 0;
    let mut score = 0;

    for (operation, count) in hsp.operations {
        match operation {
            EditOperation::Match | EditOperation::Substitution => {
                for _ in 0..count {
                    if query_begin < 0 {
                        query_begin = query_pos;
                        subject_begin = subject_pos;
                        len = 0;
                    }
                    score += score_matrix.score(
                        letter_mask(query[query_pos as usize]),
                        letter_mask(subject[subject_pos as usize]),
                    );
                    len += 1;
                    query_pos += 1;
                    subject_pos += 1;
                }
            }
            EditOperation::Deletion => {
                for _ in 0..count {
                    if query_begin >= 0 {
                        print_diag(
                            query_begin,
                            subject_begin,
                            len,
                            score,
                            diags,
                            query,
                            subject,
                            score_matrix,
                            out,
                        )?;
                        score -= score_matrix.gap_open() + score_matrix.gap_extend();
                        query_begin = -1;
                        subject_begin = -1;
                    } else {
                        score -= score_matrix.gap_extend();
                    }
                    subject_pos += 1;
                }
            }
            EditOperation::Insertion => {
                for _ in 0..count {
                    if query_begin >= 0 {
                        print_diag(
                            query_begin,
                            subject_begin,
                            len,
                            score,
                            diags,
                            query,
                            subject,
                            score_matrix,
                            out,
                        )?;
                        score -= score_matrix.gap_open() + score_matrix.gap_extend();
                        query_begin = -1;
                        subject_begin = -1;
                    } else {
                        score -= score_matrix.gap_extend();
                    }
                    query_pos += 1;
                }
            }
            EditOperation::FrameshiftForward | EditOperation::FrameshiftReverse => {}
        }
    }
    print_diag(
        query_begin,
        subject_begin,
        len,
        score,
        diags,
        query,
        subject,
        score_matrix,
        out,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;

    fn matrix() -> ScoreMatrix {
        ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap()
    }

    #[test]
    fn perfect_alignment_diagnostic_is_byte_exact() {
        let matrix = matrix();
        let query = vec![0; 4];
        let subject = vec![0; 4];
        let score = matrix.score(0, 0) * 4;
        let mut graph = DiagGraph::new();
        graph.load(&[DiagonalSegment::new(0, 0, 4, score)]);
        let mut out = Vec::new();

        smith_waterman(&query, &subject, &graph, &matrix, &mut out).unwrap();

        assert_eq!(
            out,
            format!("Diag n=0 i=0 j=0 len=4 prefix_score={score} prefix_score2={score}\n")
                .as_bytes()
        );
    }

    #[test]
    fn masked_letters_are_rescored_like_cpp_sequence_access() {
        let matrix = matrix();
        let mut query = vec![0; 4];
        let mut subject = vec![0; 4];
        query[1] |= SEED_MASK;
        subject[2] |= SEED_MASK;
        let score = matrix.score(0, 0) * 4;
        let mut graph = DiagGraph::new();
        graph.load(&[DiagonalSegment::new(0, 0, 4, score)]);
        let mut out = Vec::new();

        smith_waterman(&query, &subject, &graph, &matrix, &mut out).unwrap();

        assert_eq!(
            out,
            format!("Diag n=0 i=0 j=0 len=4 prefix_score={score} prefix_score2={score}\n")
                .as_bytes()
        );
    }

    #[test]
    fn print_diag_fallback_is_byte_exact() {
        let matrix = matrix();
        let query = vec![0; 3];
        let subject = vec![0; 3];
        let graph = DiagGraph::new();
        let mut out = Vec::new();

        print_diag(1, 1, 2, 7, &graph, &query, &subject, &matrix, &mut out).unwrap();

        assert_eq!(out, b"Diag n=x i=1 j=1 len=2 prefix_score=7\n");
    }

    #[test]
    fn deletion_splits_diagonals_and_accumulates_gap_cost_exactly() {
        let matrix = matrix();
        let query = vec![0; 20];
        let mut subject = vec![0; 10];
        subject.extend_from_slice(&[17; 6]); // W: A/W mismatches make one gap optimal.
        subject.extend_from_slice(&[0; 10]);
        let graph = DiagGraph::new();
        let mut out = Vec::new();

        smith_waterman(&query, &subject, &graph, &matrix, &mut out).unwrap();

        assert_eq!(
            out,
            b"Diag n=x i=0 j=0 len=10 prefix_score=40\n\
Diag n=x i=10 j=16 len=10 prefix_score=63\n"
        );
    }
}
