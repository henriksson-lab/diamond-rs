use crate::basic::reduction::Reduction;
use crate::basic::seed::{seed_partition, PackedSeed};
use crate::basic::shape::Shape;
use crate::basic::value::Letter;
use crate::data::seed_histogram::CURRENT_RANGE;
use crate::search::hamming::FingerPrint;
use crate::search::sse_dist::{reduced_match, seed_mask};
use crate::util::algo::PatternMatcher;
use crate::util::sequence;

#[derive(Clone)]
pub struct Context<'a> {
    pub previous_matcher: PatternMatcher,
    pub current_matcher: PatternMatcher,
    pub short_query_ungapped_cutoff: i32,
    pub seedp_mask: PackedSeed,
    pub reduction: &'a Reduction,
}

/// Matches C++ `verify_hit`.
pub fn verify_hit(
    q: &[Letter],
    s: &[Letter],
    center: usize,
    _score_cutoff: i32,
    left: bool,
    match_mask: u32,
    shape: &Shape,
    chunked: bool,
    hamming_filter_id: u32,
    seedp_mask: PackedSeed,
    reduction: &Reduction,
) -> bool {
    if chunked && (shape.mask & match_mask) == shape.mask {
        let Some(seed) = shape.set_seed(&s[center..], reduction) else {
            return false;
        };
        let p = seed_partition(seed, seedp_mask);
        let current_range = CURRENT_RANGE.lock().unwrap();
        if left && !current_range.lower_or_equal(p) {
            return false;
        }
        if !left && !current_range.lower(p) {
            return false;
        }
    }
    if center < 16 || center + 32 > q.len() || center + 32 > s.len() {
        return false;
    }
    let fq = FingerPrint::from_seq_center(q, center);
    let fs = FingerPrint::from_seq_center(s, center);
    fq.match_count(&fs) >= hamming_filter_id
}

/// Matches C++ `verify_hits`.
pub fn verify_hits(
    mut mask: u32,
    q: &[Letter],
    s: &[Letter],
    offset: usize,
    score_cutoff: i32,
    left: bool,
    match_mask: u32,
    shape: &Shape,
    chunked: bool,
    hamming_filter_id: u32,
    seedp_mask: PackedSeed,
    reduction: &Reduction,
) -> bool {
    let mut shift = 0usize;
    while mask != 0 {
        let i = mask.trailing_zeros() as usize;
        let center = offset + i + shift;
        if verify_hit(
            q,
            s,
            center,
            score_cutoff,
            left,
            match_mask >> (i + shift),
            shape,
            chunked,
            hamming_filter_id,
            seedp_mask,
            reduction,
        ) {
            return true;
        }
        mask >>= i + 1;
        shift += i + 1;
    }
    false
}

/// Matches C++ `left_most_filter`.
pub fn left_most_filter(
    query: &[Letter],
    subject: &[Letter],
    seed_offset: i32,
    seed_len: i32,
    context: &Context,
    first_shape: bool,
    shape: &Shape,
    score_cutoff: i32,
    chunked: bool,
    hamming_filter_id: u32,
) -> bool {
    const WINDOW_LEFT: i32 = 16;
    const WINDOW_RIGHT: i32 = 32;

    let d = (seed_offset - WINDOW_LEFT).max(0) as usize;
    let mut window_left = seed_offset.min(WINDOW_LEFT).max(0) as usize;
    let q0 = &query[d..];
    let s0 = &subject[d..];
    let mut window = query.len() - d;
    window = window.min(window_left + 1 + WINDOW_RIGHT as usize);

    let subject_clipped = sequence::clip(&s0[..window.min(s0.len())], window_left as i32);
    window = window.saturating_sub(s0[..window.min(s0.len())].len() - subject_clipped.len());

    let clipped_start = subject_clipped.as_ptr() as usize - s0.as_ptr() as usize;
    let clipped_start = clipped_start / std::mem::size_of::<Letter>();
    let q = &q0[clipped_start..];
    let s = &s0[clipped_start..];
    window_left -= clipped_start;
    window -= clipped_start;

    let match_mask = reduced_match(q, s, window as i32, context.reduction);
    let query_seed_mask = !seed_mask(q, window as i32);

    let len_left = window_left as u32 + seed_len as u32 - 1;
    let match_mask_left = (((1u64 << len_left) - 1) & match_mask) as u32;
    let query_mask_left = (((1u64 << len_left) - 1) & query_seed_mask) as u32;

    let left_hit = context.current_matcher.hit(match_mask_left, len_left) & query_mask_left;

    if first_shape && !chunked {
        return left_hit == 0
            || !verify_hits(
                left_hit,
                q,
                s,
                0,
                score_cutoff,
                true,
                match_mask_left,
                shape,
                chunked,
                hamming_filter_id,
                context.seedp_mask,
                context.reduction,
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

    (left_hit == 0
        || !verify_hits(
            left_hit,
            q,
            s,
            0,
            score_cutoff,
            true,
            match_mask_left,
            shape,
            chunked,
            hamming_filter_id,
            context.seedp_mask,
            context.reduction,
        ))
        && (right_hit == 0
            || !verify_hits(
                right_hit,
                q,
                s,
                window_left + 1,
                score_cutoff,
                false,
                match_mask_right,
                shape,
                chunked,
                hamming_filter_id,
                context.seedp_mask,
                context.reduction,
            ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::seed::seedp_mask;
    use crate::basic::shape::Shape;
    use crate::data::seed_histogram::SeedPartitionRange;

    fn context<'a>(shape: &Shape, reduction: &'a Reduction) -> Context<'a> {
        Context {
            previous_matcher: PatternMatcher::new(&[shape.mask]),
            current_matcher: PatternMatcher::new(&[shape.mask]),
            short_query_ungapped_cutoff: 0,
            seedp_mask: seedp_mask(8),
            reduction,
        }
    }

    #[test]
    fn test_verify_hit_hamming_filter() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);
        let q: Vec<Letter> = (0..80).map(|i| (i % 20) as Letter).collect();
        let mut s = q.clone();
        assert!(verify_hit(
            &q,
            &s,
            24,
            0,
            true,
            shape.mask,
            &shape,
            false,
            48,
            seedp_mask(8),
            &reduction,
        ));
        s[10] = 9;
        assert!(!verify_hit(
            &q,
            &s,
            24,
            0,
            true,
            shape.mask,
            &shape,
            false,
            48,
            seedp_mask(8),
            &reduction,
        ));
        assert!(verify_hit(
            &q,
            &s,
            24,
            0,
            true,
            shape.mask,
            &shape,
            false,
            47,
            seedp_mask(8),
            &reduction,
        ));
    }

    #[test]
    fn test_verify_hits_scans_set_bits() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);
        let q: Vec<Letter> = (0..96).map(|i| (i % 20) as Letter).collect();
        let mut s = q.clone();
        s[0] = 9;
        assert!(verify_hits(
            0b101,
            &q,
            &s,
            16,
            0,
            true,
            shape.mask,
            &shape,
            false,
            48,
            seedp_mask(8),
            &reduction,
        ));
    }

    #[test]
    fn test_verify_hit_chunked_partition_range() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);
        let q: Vec<Letter> = (0..80).map(|i| (i % 20) as Letter).collect();
        let s = q.clone();
        let seed = shape.set_seed(&s[24..], &reduction).unwrap();
        let p = seed_partition(seed, seedp_mask(8));
        {
            let mut r = CURRENT_RANGE.lock().unwrap();
            *r = SeedPartitionRange::with_bounds(0, p);
        }
        assert!(!verify_hit(
            &q,
            &s,
            24,
            0,
            true,
            shape.mask,
            &shape,
            true,
            48,
            seedp_mask(8),
            &reduction,
        ));
        {
            let mut r = CURRENT_RANGE.lock().unwrap();
            *r = SeedPartitionRange::with_bounds(p, p + 1);
        }
        assert!(verify_hit(
            &q,
            &s,
            24,
            0,
            true,
            shape.mask,
            &shape,
            true,
            48,
            seedp_mask(8),
            &reduction,
        ));
        *CURRENT_RANGE.lock().unwrap() = SeedPartitionRange::new();
    }

    #[test]
    fn test_left_most_filter_accepts_current_leftmost_seed() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);
        let context = context(&shape, &reduction);
        let q: Vec<Letter> = (0..96).map(|i| (i % 20) as Letter).collect();
        let s = q.clone();
        assert!(left_most_filter(
            &q,
            &s,
            32,
            shape.length,
            &context,
            true,
            &shape,
            0,
            false,
            48,
        ));
    }

    #[test]
    fn test_left_most_filter_rejects_prior_identical_seed() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111", &reduction);
        let context = context(&shape, &reduction);
        let q: Vec<Letter> = (0..96).map(|i| (i % 20) as Letter).collect();
        let s = q.clone();
        assert!(!left_most_filter(
            &q,
            &s,
            32,
            shape.length,
            &context,
            false,
            &shape,
            0,
            false,
            48,
        ));
    }
}
