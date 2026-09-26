//! Chaining graph backtrace.
//!
//! This mirrors `diamond/src/chaining/backtrace.cpp`. C++ overloads are given
//! descriptive snake-case suffixes, while the pre-existing approximate-only
//! Rust methods remain compatibility wrappers.

use crate::align::hsp::Hsp;
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::letter_mask;
use crate::dp::ungapped::DiagonalSegment;
use crate::util::hsp::{Anchor, ApproxHsp};

use super::{Aligner, DEFAULT_CHAINING_STACKED_HSP_RATIO, END};

pub fn disjoint_approx_hsp_with_ratio(
    hsps: &[ApproxHsp],
    candidate: &ApproxHsp,
    cutoff: i32,
    stacked_hsp_ratio: f64,
) -> bool {
    for previous in hsps {
        let target_overlap = candidate
            .subject_range
            .overlap_factor(&previous.subject_range);
        let query_overlap = candidate.query_range.overlap_factor(&previous.query_range);
        if (1.0 - target_overlap.min(query_overlap)) * candidate.score as f64
            / previous.score as f64
            >= stacked_hsp_ratio
        {
            continue;
        }
        if (1.0 - target_overlap.max(query_overlap)) * (candidate.score as f64) < cutoff as f64 {
            return false;
        }
    }
    true
}

pub fn disjoint_approx_hsp(hsps: &[ApproxHsp], candidate: &ApproxHsp, cutoff: i32) -> bool {
    disjoint_approx_hsp_with_ratio(hsps, candidate, cutoff, DEFAULT_CHAINING_STACKED_HSP_RATIO)
}

pub fn disjoint_diagonal_segment_with_ratio(
    hsps: &[ApproxHsp],
    candidate: &DiagonalSegment,
    cutoff: i32,
    stacked_hsp_ratio: f64,
) -> bool {
    for previous in hsps {
        let target_overlap = candidate
            .subject_range()
            .overlap_factor(&previous.subject_range);
        let query_overlap = candidate
            .query_range()
            .overlap_factor(&previous.query_range);
        if (1.0 - target_overlap.min(query_overlap)) * candidate.score as f64
            / previous.score as f64
            >= stacked_hsp_ratio
        {
            continue;
        }
        if (1.0 - target_overlap.max(query_overlap)) * (candidate.score as f64) < cutoff as f64 {
            return false;
        }
    }
    true
}

pub fn disjoint_diagonal_segment(
    hsps: &[ApproxHsp],
    candidate: &DiagonalSegment,
    cutoff: i32,
) -> bool {
    disjoint_diagonal_segment_with_ratio(
        hsps,
        candidate,
        cutoff,
        DEFAULT_CHAINING_STACKED_HSP_RATIO,
    )
}

#[derive(Clone, Copy)]
struct BacktraceNode {
    node: usize,
    score_min: i32,
    subject_end: i32,
}

impl Aligner<'_> {
    /// Recursive C++ `Aligner::backtrace_old`, including optional transcript output.
    #[allow(clippy::too_many_arguments)]
    pub fn backtrace_old_with_hsp(
        &self,
        node: usize,
        subject_end: i32,
        mut out: Option<&mut Hsp>,
        traits: &mut ApproxHsp,
        score_max: i32,
        mut score_min: i32,
        max_shift: i32,
        next: &mut u32,
    ) -> bool {
        let diagonal = &self.diags.nodes[node];
        let edge_idx = self.diags.get_edge(node, subject_end);
        let mut at_end = edge_idx.is_none();
        let prefix_score = edge_idx
            .map(|idx| self.diags.edges[idx].prefix_score)
            .unwrap_or(diagonal.score);
        if prefix_score > score_max {
            return false;
        }

        score_min = score_min.min(
            edge_idx
                .map(|idx| self.diags.edges[idx].prefix_score_begin)
                .unwrap_or(0),
        );

        let mut subject_pos = edge_idx
            .map(|idx| self.diags.edges[idx].j)
            .unwrap_or(diagonal.j);
        if let Some(idx) = edge_idx {
            let edge = &self.diags.edges[idx];
            let predecessor = &self.diags.nodes[edge.node_out as usize];
            let shift = diagonal.diag() - predecessor.diag();

            if shift.abs() <= max_shift {
                let traced = self.backtrace_old_with_hsp(
                    edge.node_out as usize,
                    if shift > 0 {
                        subject_pos
                    } else {
                        subject_pos + shift
                    },
                    out.as_deref_mut(),
                    traits,
                    score_max,
                    score_min,
                    max_shift,
                    next,
                );
                if !traced {
                    if edge.prefix_score_begin > score_min {
                        return false;
                    }
                    at_end = true;
                }
            } else {
                *next = edge.node_out;
                at_end = true;
            }
        }

        if at_end {
            if let Some(hsp) = out.as_deref_mut() {
                hsp.query_range.begin = diagonal.i;
                hsp.subject_range.begin = diagonal.j;
                hsp.score = score_max - score_min;
            }
            traits.query_range.begin = diagonal.i;
            traits.subject_range.begin = diagonal.j;
            traits.score = score_max - score_min;
            subject_pos = diagonal.j;
        } else if let (Some(idx), Some(hsp)) = (edge_idx, out.as_deref_mut()) {
            let edge = &self.diags.edges[idx];
            let predecessor = &self.diags.nodes[edge.node_out as usize];
            let shift = diagonal.diag() - predecessor.diag();
            if shift > 0 {
                hsp.transcript
                    .push_with_count(EditOperation::Insertion, shift as u32);
                hsp.length += shift;
            } else if shift < 0 {
                for j in edge.j + shift..edge.j {
                    hsp.transcript.push_with_letter(
                        EditOperation::Deletion,
                        letter_mask(self.subject[j as usize]),
                    );
                    hsp.length += 1;
                }
            }
        }

        let diag = diagonal.diag();
        traits.d_max = traits.d_max.max(diag);
        traits.d_min = traits.d_min.min(diag);
        if diagonal.score > traits.max_diag.segment.score {
            traits.max_diag = Anchor::from_diagonal_segment(diagonal.segment.clone());
            traits.max_diag.prefix_score = prefix_score;
            traits.max_diag.d_max_left = traits
                .max_diag
                .d_max_right
                .max(traits.max_diag.d_max_left)
                .max(diag);
            traits.max_diag.d_min_left = traits
                .max_diag
                .d_min_right
                .min(traits.max_diag.d_min_left)
                .min(diag);
            traits.max_diag.d_max_right = diag;
            traits.max_diag.d_min_right = diag;
        } else {
            traits.max_diag.d_max_right = traits.max_diag.d_max_right.max(diag);
            traits.max_diag.d_min_right = traits.max_diag.d_min_right.min(diag);
        }

        if let Some(hsp) = out.as_deref_mut() {
            while subject_pos < subject_end {
                let subject_letter = letter_mask(self.subject[subject_pos as usize]);
                let query_letter = letter_mask(self.query[(diag + subject_pos) as usize]);
                if subject_letter == query_letter {
                    hsp.transcript.push(EditOperation::Match);
                    hsp.identities += 1;
                } else {
                    hsp.transcript
                        .push_with_letter(EditOperation::Substitution, subject_letter);
                }
                hsp.length += 1;
                subject_pos += 1;
            }
        }
        true
    }

    /// Approximate-only compatibility wrapper for the former Rust API.
    pub fn backtrace_old(
        &self,
        node: usize,
        subject_end: i32,
        traits: &mut ApproxHsp,
        score_max: i32,
        score_min: i32,
        max_shift: i32,
        next: &mut u32,
    ) -> bool {
        self.backtrace_old_with_hsp(
            node,
            subject_end,
            None,
            traits,
            score_max,
            score_min,
            max_shift,
            next,
        )
    }

    /// Iterative C++ `Aligner::backtrace(node, j_end, ...)` overload.
    #[allow(clippy::too_many_arguments)]
    pub fn backtrace_iterative(
        &self,
        node: usize,
        subject_end: i32,
        mut out: Option<&mut Hsp>,
        traits: &mut ApproxHsp,
        score_max: i32,
        score_min: i32,
        max_shift: i32,
        next: &mut u32,
    ) {
        let mut nodes = vec![BacktraceNode {
            node,
            score_min,
            subject_end,
        }];
        let mut returned = false;
        let mut return_value = false;

        while let Some(current) = nodes.last().copied() {
            let diagonal = &self.diags.nodes[current.node];
            let edge_idx = self.diags.get_edge(current.node, current.subject_end);
            let mut at_end = edge_idx.is_none();
            let prefix_score = edge_idx
                .map(|idx| self.diags.edges[idx].prefix_score)
                .unwrap_or(diagonal.score);
            if !returned && prefix_score > score_max {
                returned = true;
                return_value = false;
                nodes.pop();
                continue;
            }

            let score_min = current.score_min.min(
                edge_idx
                    .map(|idx| self.diags.edges[idx].prefix_score_begin)
                    .unwrap_or(0),
            );
            nodes.last_mut().unwrap().score_min = score_min;

            if let Some(idx) = edge_idx {
                let edge = &self.diags.edges[idx];
                let predecessor = &self.diags.nodes[edge.node_out as usize];
                let shift = diagonal.diag() - predecessor.diag();
                if shift.abs() <= max_shift {
                    if !returned {
                        nodes.push(BacktraceNode {
                            node: edge.node_out as usize,
                            score_min,
                            subject_end: if shift > 0 { edge.j } else { edge.j + shift },
                        });
                        continue;
                    }
                    if !return_value {
                        if edge.prefix_score_begin > score_min {
                            return_value = false;
                            nodes.pop();
                            continue;
                        }
                        at_end = true;
                    }
                } else {
                    *next = edge.node_out;
                    at_end = true;
                }
            }

            if at_end {
                if let Some(hsp) = out.as_deref_mut() {
                    hsp.query_range.begin = diagonal.i;
                    hsp.subject_range.begin = diagonal.j;
                    hsp.score = score_max - score_min;
                }
                traits.query_range.begin = diagonal.i;
                traits.subject_range.begin = diagonal.j;
                traits.score = score_max - score_min;
            } else if let (Some(idx), Some(hsp)) = (edge_idx, out.as_deref_mut()) {
                let edge = &self.diags.edges[idx];
                let predecessor = &self.diags.nodes[edge.node_out as usize];
                let shift = diagonal.diag() - predecessor.diag();
                if shift > 0 {
                    hsp.transcript
                        .push_with_count(EditOperation::Insertion, shift as u32);
                    hsp.length += shift;
                } else if shift < 0 {
                    for j in edge.j + shift..edge.j {
                        hsp.transcript.push_with_letter(
                            EditOperation::Deletion,
                            letter_mask(self.subject[j as usize]),
                        );
                        hsp.length += 1;
                    }
                }
            }

            let diag = diagonal.diag();
            traits.d_max = traits.d_max.max(diag);
            traits.d_min = traits.d_min.min(diag);
            if let Some(hsp) = out.as_deref_mut() {
                let mut subject_pos = edge_idx
                    .filter(|_| !at_end)
                    .map(|idx| self.diags.edges[idx].j)
                    .unwrap_or(diagonal.j);
                while subject_pos < current.subject_end {
                    let subject_letter = letter_mask(self.subject[subject_pos as usize]);
                    let query_letter = letter_mask(self.query[(diag + subject_pos) as usize]);
                    if subject_letter == query_letter {
                        hsp.transcript.push(EditOperation::Match);
                        hsp.identities += 1;
                    } else {
                        hsp.transcript
                            .push_with_letter(EditOperation::Substitution, subject_letter);
                    }
                    hsp.length += 1;
                    subject_pos += 1;
                }
            }
            returned = true;
            return_value = true;
            nodes.pop();
        }
    }

    /// C++ top-node overload with optional full transcript output.
    pub fn backtrace_with_hsp(
        &self,
        top_node: usize,
        mut out: Option<&mut Hsp>,
        traits: &mut ApproxHsp,
        max_shift: i32,
        next: &mut u32,
        max_subject_end: i32,
    ) {
        let mut result = ApproxHsp::new(self.frame, 0);
        if top_node != END {
            let diagonal = &self.diags.nodes[top_node];
            if let Some(hsp) = out.as_deref_mut() {
                hsp.transcript.clear();
                hsp.query_range.end = diagonal.query_end();
                hsp.subject_range.end = diagonal.subject_end();
            }
            result.subject_range.end = diagonal.subject_end();
            result.query_range.end = diagonal.query_end();
            self.backtrace_old_with_hsp(
                top_node,
                diagonal.subject_end().min(max_subject_end),
                out.as_deref_mut(),
                &mut result,
                diagonal.prefix_score,
                diagonal.prefix_score,
                max_shift,
                next,
            );
        } else {
            result.score = 0;
            if let Some(hsp) = out.as_deref_mut() {
                hsp.score = 0;
            }
        }
        if let Some(hsp) = out.as_deref_mut() {
            hsp.transcript.push_terminator();
        }
        *traits = result;
    }

    pub fn backtrace(
        &self,
        top_node: usize,
        traits: &mut ApproxHsp,
        max_shift: i32,
        next: &mut u32,
        max_subject_end: i32,
    ) {
        self.backtrace_with_hsp(top_node, None, traits, max_shift, next, max_subject_end);
    }

    /// C++ per-top-node overload, including its optional debug HSP list.
    #[allow(clippy::too_many_arguments)]
    pub fn backtrace_top_with_hsps(
        &self,
        top_node: usize,
        hsps: &mut Vec<Hsp>,
        traits: &mut Vec<ApproxHsp>,
        traits_begin: &mut usize,
        cutoff: i32,
        max_shift: i32,
        stacked_hsp_ratio: f64,
    ) -> i32 {
        let mut top_node = top_node;
        let mut max_score = 0;
        let mut max_subject_end = self.subject.len() as i32;
        loop {
            let mut hsp = self.log.then(|| Hsp::with_backtraced(true));
            let mut candidate = ApproxHsp::new(self.frame, 0);
            let mut next = u32::MAX;
            self.backtrace_with_hsp(
                top_node,
                hsp.as_mut(),
                &mut candidate,
                max_shift,
                &mut next,
                max_subject_end,
            );
            if candidate.score > 0 {
                max_subject_end = candidate.subject_range.begin;
            }
            if candidate.score >= cutoff
                && disjoint_approx_hsp_with_ratio(
                    &traits[*traits_begin..],
                    &candidate,
                    cutoff,
                    stacked_hsp_ratio,
                )
            {
                let was_end = *traits_begin == traits.len();
                max_score = max_score.max(candidate.score);
                traits.push(candidate);
                if was_end {
                    *traits_begin = traits.len() - 1;
                }
                if let Some(hsp) = hsp {
                    hsps.push(hsp);
                }
            }
            if next == u32::MAX {
                break;
            }
            top_node = next as usize;
        }
        max_score
    }

    pub fn backtrace_top(
        &self,
        top_node: usize,
        traits: &mut Vec<ApproxHsp>,
        traits_begin: &mut usize,
        cutoff: i32,
        max_shift: i32,
    ) -> i32 {
        self.backtrace_top_with_hsps(
            top_node,
            &mut Vec::new(),
            traits,
            traits_begin,
            cutoff,
            max_shift,
            DEFAULT_CHAINING_STACKED_HSP_RATIO,
        )
    }

    /// C++ all-candidates overload with explicit global-config replacement.
    pub fn backtrace_all_with_hsps(
        &self,
        hsps: &mut Vec<Hsp>,
        traits: &mut Vec<ApproxHsp>,
        cutoff: i32,
        max_shift: i32,
        stacked_hsp_ratio: f64,
    ) -> i32 {
        let mut top_nodes: Vec<usize> = self
            .diags
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, diagonal)| (diagonal.rel_score() >= cutoff).then_some(i))
            .collect();
        top_nodes.sort_by(|&x, &y| {
            self.diags.nodes[y]
                .rel_score()
                .cmp(&self.diags.nodes[x].rel_score())
        });

        let mut max_score = 0;
        let mut traits_begin = traits.len();
        for node in top_nodes {
            if disjoint_diagonal_segment_with_ratio(
                &traits[traits_begin..],
                &self.diags.nodes[node].segment,
                cutoff,
                stacked_hsp_ratio,
            ) {
                max_score = max_score.max(self.backtrace_top_with_hsps(
                    node,
                    hsps,
                    traits,
                    &mut traits_begin,
                    cutoff,
                    max_shift,
                    stacked_hsp_ratio,
                ));
            }
        }
        max_score
    }

    pub fn backtrace_all(&self, traits: &mut Vec<ApproxHsp>, cutoff: i32, max_shift: i32) -> i32 {
        self.backtrace_all_with_hsps(
            &mut Vec::new(),
            traits,
            cutoff,
            max_shift,
            DEFAULT_CHAINING_STACKED_HSP_RATIO,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::EditOperation;
    use crate::basic::value::SEED_MASK;
    use crate::stats::score_matrix::ScoreMatrix;
    use crate::util::interval::Interval;

    fn aligner<'a>(
        query: &'a [i8],
        subject: &'a [i8],
        score_matrix: &'a ScoreMatrix,
    ) -> Aligner<'a> {
        let score = score_matrix.score(0, 0) * 4;
        let segments = [
            DiagonalSegment::new(0, 0, 4, score),
            DiagonalSegment::new(5, 5, 4, score),
        ];
        let mut aligner = Aligner::new(query, subject, true, 2, score_matrix);
        aligner.diags.load(&segments);
        aligner.forward_pass(
            0..aligner.diags.nodes.len(),
            true,
            super::super::SPACE_PENALTY,
        );
        aligner
    }

    #[test]
    fn transcript_backtrace_preserves_ranges_operations_and_terminator() {
        let matrix = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0; 9];
        let mut subject = vec![0; 9];
        subject[4] = 1;
        let aligner = aligner(&query, &subject, &matrix);
        let top = aligner.diags.top_node();
        let mut hsp = Hsp::with_backtraced(true);
        let mut traits = ApproxHsp::new(2, 0);
        let mut next = u32::MAX;

        aligner.backtrace_with_hsp(
            top,
            Some(&mut hsp),
            &mut traits,
            100,
            &mut next,
            subject.len() as i32,
        );

        assert_eq!(next, u32::MAX);
        assert_eq!(hsp.query_range, Interval::new(0, 9));
        assert_eq!(hsp.subject_range, Interval::new(0, 9));
        assert_eq!(hsp.length, 9);
        assert_eq!(hsp.identities, 8);
        let operations: Vec<_> = hsp
            .transcript
            .iter()
            .map(|operation| (operation.op, operation.count, operation.letter))
            .collect();
        assert_eq!(
            operations,
            vec![
                (EditOperation::Match, 4, 0),
                (EditOperation::Substitution, 1, 1),
                (EditOperation::Match, 4, 0),
            ]
        );
        assert!(hsp.transcript.data().last().unwrap().is_terminator());
        assert_eq!(traits.query_range, hsp.query_range);
        assert_eq!(traits.subject_range, hsp.subject_range);
        assert_eq!(traits.score, hsp.score);
    }

    #[test]
    fn iterative_and_recursive_backtraces_agree_on_core_result() {
        let matrix = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let query = vec![0; 9];
        let subject = vec![0; 9];
        let aligner = aligner(&query, &subject, &matrix);
        let top = aligner.diags.top_node();
        let diagonal = &aligner.diags.nodes[top];
        let mut recursive_hsp = Hsp::with_backtraced(true);
        let mut recursive = ApproxHsp::new(2, 0);
        let mut recursive_next = u32::MAX;
        assert!(aligner.backtrace_old_with_hsp(
            top,
            diagonal.subject_end(),
            Some(&mut recursive_hsp),
            &mut recursive,
            diagonal.prefix_score,
            diagonal.prefix_score,
            100,
            &mut recursive_next,
        ));

        let mut iterative_hsp = Hsp::with_backtraced(true);
        let mut iterative = ApproxHsp::new(2, 0);
        let mut iterative_next = u32::MAX;
        aligner.backtrace_iterative(
            top,
            diagonal.subject_end(),
            Some(&mut iterative_hsp),
            &mut iterative,
            diagonal.prefix_score,
            diagonal.prefix_score,
            100,
            &mut iterative_next,
        );

        assert_eq!(iterative.score, recursive.score);
        assert_eq!(iterative.query_range.begin, recursive.query_range.begin);
        assert_eq!(iterative.subject_range.begin, recursive.subject_range.begin);
        assert_eq!(iterative.d_min, recursive.d_min);
        assert_eq!(iterative.d_max, recursive.d_max);
        assert_eq!(iterative_hsp.length, recursive_hsp.length);
        assert_eq!(iterative_hsp.identities, recursive_hsp.identities);
        assert_eq!(
            iterative_hsp.transcript.data(),
            recursive_hsp.transcript.data()
        );
    }

    #[test]
    fn transcript_backtrace_emits_insertion_and_deletion_shifts() {
        let matrix = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let cases = [
            (
                vec![0; 10],
                vec![0; 9],
                DiagonalSegment::new(6, 5, 4, matrix.score(0, 0) * 4),
                EditOperation::Insertion,
            ),
            (
                vec![0; 9],
                vec![0; 10],
                DiagonalSegment::new(5, 6, 4, matrix.score(0, 0) * 4),
                EditOperation::Deletion,
            ),
        ];

        for (query, subject, second, expected_gap) in cases {
            let first = DiagonalSegment::new(0, 0, 4, matrix.score(0, 0) * 4);
            let mut aligner = Aligner::new(&query, &subject, true, 0, &matrix);
            aligner.diags.load(&[first, second]);
            aligner.forward_pass(
                0..aligner.diags.nodes.len(),
                true,
                super::super::SPACE_PENALTY,
            );
            assert_eq!(aligner.diags.edges.len(), 1);

            let mut hsp = Hsp::with_backtraced(true);
            let mut traits = ApproxHsp::new(0, 0);
            let mut next = u32::MAX;
            aligner.backtrace_with_hsp(
                aligner.diags.top_node(),
                Some(&mut hsp),
                &mut traits,
                10,
                &mut next,
                subject.len() as i32,
            );
            let operations: Vec<_> = hsp
                .transcript
                .iter()
                .map(|operation| operation.op)
                .collect();
            assert!(operations.contains(&expected_gap));
            assert_eq!(hsp.length, 10);
        }
    }

    #[test]
    fn masked_letters_compare_and_serialize_as_unmasked() {
        let matrix = ScoreMatrix::new("BLOSUM62", -1, -1, -1, 1, 0).unwrap();
        let mut query = vec![0; 9];
        let mut subject = vec![0; 9];
        query[2] |= SEED_MASK;
        subject[3] |= SEED_MASK;
        let aligner = aligner(&query, &subject, &matrix);
        let top = aligner.diags.top_node();
        let diagonal = &aligner.diags.nodes[top];

        let mut recursive_hsp = Hsp::with_backtraced(true);
        let mut recursive = ApproxHsp::new(0, 0);
        let mut next = u32::MAX;
        assert!(aligner.backtrace_old_with_hsp(
            top,
            diagonal.subject_end(),
            Some(&mut recursive_hsp),
            &mut recursive,
            diagonal.prefix_score,
            diagonal.prefix_score,
            100,
            &mut next,
        ));
        assert_eq!(recursive_hsp.identities, 9);
        assert_eq!(
            recursive_hsp
                .transcript
                .iter()
                .map(|operation| operation.op)
                .collect::<Vec<_>>(),
            [EditOperation::Match]
        );

        let mut iterative_hsp = Hsp::with_backtraced(true);
        let mut iterative = ApproxHsp::new(0, 0);
        aligner.backtrace_iterative(
            top,
            diagonal.subject_end(),
            Some(&mut iterative_hsp),
            &mut iterative,
            diagonal.prefix_score,
            diagonal.prefix_score,
            100,
            &mut next,
        );
        assert_eq!(iterative_hsp.identities, 9);

        // Deletion payloads also come through Sequence::operator[] in C++.
        let mut deletion_subject = vec![0; 10];
        deletion_subject[5] |= SEED_MASK;
        let deletion_query = vec![0; 9];
        let first = DiagonalSegment::new(0, 0, 4, matrix.score(0, 0) * 4);
        let second = DiagonalSegment::new(5, 6, 4, matrix.score(0, 0) * 4);
        let mut deletion_aligner =
            Aligner::new(&deletion_query, &deletion_subject, true, 0, &matrix);
        deletion_aligner.diags.load(&[first, second]);
        deletion_aligner.forward_pass(
            0..deletion_aligner.diags.nodes.len(),
            true,
            super::super::SPACE_PENALTY,
        );
        let mut deletion_hsp = Hsp::with_backtraced(true);
        let mut deletion_traits = ApproxHsp::new(0, 0);
        deletion_aligner.backtrace_with_hsp(
            deletion_aligner.diags.top_node(),
            Some(&mut deletion_hsp),
            &mut deletion_traits,
            10,
            &mut next,
            deletion_subject.len() as i32,
        );
        let deletion = deletion_hsp
            .transcript
            .iter()
            .find(|operation| operation.op == EditOperation::Deletion)
            .unwrap();
        assert_eq!(deletion.letter, 0);
    }

    #[test]
    fn stacked_hsp_ratio_is_explicit_and_configurable() {
        let previous = ApproxHsp::from_parts(
            0,
            0,
            100,
            0,
            Interval::new(0, 10),
            Interval::new(0, 10),
            Anchor::default(),
            0.0,
        );
        let candidate = ApproxHsp::from_parts(
            0,
            0,
            100,
            0,
            Interval::new(5, 15),
            Interval::new(5, 15),
            Anchor::default(),
            0.0,
        );
        assert!(disjoint_approx_hsp_with_ratio(
            &[previous.clone()],
            &candidate,
            60,
            0.5
        ));
        assert!(!disjoint_approx_hsp_with_ratio(
            &[previous],
            &candidate,
            60,
            0.6
        ));
    }
}
