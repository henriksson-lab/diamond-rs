//! Left-most filtering for a subject window already proven to be wholly inside
//! one database sequence.
//!
//! The regular implementation must search for delimiter letters at both ends
//! of every candidate window.  Its caller can often prove from the sequence
//! index that no delimiter is present; this variant preserves the filtering
//! logic while eliding that redundant scan.

use crate::basic::shape::Shape;
use crate::basic::value::Letter;
use crate::data::seed_histogram::SeedPartitionRange;
use crate::search::hamming::match_positions;
use crate::search::left_most::Context;
use crate::search::sse_dist::{reduced_match, seed_mask};

#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub fn left_most_filter_with_range_unclipped(
    query: &[Letter],
    query_start: usize,
    query_len: usize,
    subject: &[Letter],
    subject_start: usize,
    seed_offset: i32,
    seed_len: i32,
    context: &Context,
    first_shape: bool,
    shape: &Shape,
    score_cutoff: i32,
    chunked: bool,
    hamming_filter_id: u32,
    current_range: Option<SeedPartitionRange>,
) -> bool {
    const WINDOW_LEFT: i32 = 16;
    const WINDOW_RIGHT: i32 = 32;

    let d = (seed_offset - WINDOW_LEFT).max(0) as usize;
    let window_left = seed_offset.min(WINDOW_LEFT).max(0) as usize;
    debug_assert!(query_start + query_len <= query.len());
    debug_assert!(subject_start + query_len <= subject.len());
    let q_base = query_start + d;
    let s_base = subject_start + d;
    let q = &query[q_base..query_start + query_len];
    let s = &subject[s_base..subject_start + query_len];
    let window = (query_len - d).min(window_left + 1 + WINDOW_RIGHT as usize);

    let match_mask = reduced_match(q, s, window as i32, context.reduction);
    let query_seed_mask = !seed_mask(q, window as i32);

    let len_left = window_left as u32 + seed_len as u32 - 1;
    let match_mask_left = (((1u64 << len_left) - 1) & match_mask) as u32;
    let query_mask_left = (((1u64 << len_left) - 1) & query_seed_mask) as u32;
    let left_hit = context.current_matcher.hit(match_mask_left, len_left) & query_mask_left;

    if first_shape && !chunked {
        if left_hit == 0 {
            return true;
        }
        return !verify_hits_backed(
            left_hit,
            query,
            q_base,
            subject,
            s_base,
            0,
            score_cutoff,
            true,
            match_mask_left,
            shape,
            chunked,
            hamming_filter_id,
            context.seedp_mask,
            context.reduction,
            current_range,
        );
    }

    let len_right = window as u32 - window_left as u32 - 1;
    let match_mask_right = (match_mask >> (window_left + 1)) as u32;
    let query_mask_right = (query_seed_mask >> (window_left + 1)) as u32;
    let right_matcher = if chunked {
        &context.current_matcher
    } else {
        &context.previous_matcher
    };
    let right_hit = right_matcher.hit(match_mask_right, len_right) & query_mask_right;

    if left_hit == 0 && right_hit == 0 {
        return true;
    }

    (left_hit == 0
        || !verify_hits_backed(
            left_hit,
            query,
            q_base,
            subject,
            s_base,
            0,
            score_cutoff,
            true,
            match_mask_left,
            shape,
            chunked,
            hamming_filter_id,
            context.seedp_mask,
            context.reduction,
            current_range,
        ))
        && (right_hit == 0
            || !verify_hits_backed(
                right_hit,
                query,
                q_base,
                subject,
                s_base,
                window_left + 1,
                score_cutoff,
                false,
                match_mask_right,
                shape,
                chunked,
                hamming_filter_id,
                context.seedp_mask,
                context.reduction,
                current_range,
            ))
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn verify_hits_backed(
    mut mask: u32,
    query: &[Letter],
    query_base: usize,
    subject: &[Letter],
    subject_base: usize,
    offset: usize,
    _score_cutoff: i32,
    left: bool,
    match_mask: u32,
    shape: &Shape,
    chunked: bool,
    hamming_filter_id: u32,
    seedp_mask: crate::basic::seed::PackedSeed,
    reduction: &crate::basic::reduction::Reduction,
    current_range: Option<SeedPartitionRange>,
) -> bool {
    let mut shift = 0usize;
    while mask != 0 {
        let i = mask.trailing_zeros() as usize;
        let relative_center = offset + i + shift;
        let query_center = query_base + relative_center;
        let subject_center = subject_base + relative_center;
        let mut eligible = true;
        if chunked && (shape.mask & (match_mask >> (i + shift))) == shape.mask {
            if subject_center + shape.length as usize > subject.len() {
                eligible = false;
            } else {
                // SAFETY: the fit check covers every selected shape position.
                let partition = unsafe {
                    shape.seed_partition_unchecked(
                        &subject[subject_center..],
                        reduction,
                        seedp_mask,
                    )
                };
                if let Some(partition) = partition {
                    if let Some(range) = current_range {
                        eligible = (left && range.lower_or_equal(partition))
                            || (!left && range.lower(partition));
                    } else {
                        let range = crate::data::seed_histogram::CURRENT_RANGE.lock().unwrap();
                        eligible = (left && range.lower_or_equal(partition))
                            || (!left && range.lower(partition));
                    }
                } else {
                    eligible = false;
                }
            }
        }
        if eligible
            && match_positions(query, query_center, subject, subject_center) >= hamming_filter_id
        {
            return true;
        }
        mask >>= i + 1;
        shift += i + 1;
    }
    false
}
