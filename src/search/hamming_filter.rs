//! Stage-1 hamming filter for seed matches.
//!
//! Ports C++ `search/hamming/kernel.h:all_vs_all`: for each seed hit, count
//! exact-letter matches over a 48-position window (16 letters before the seed
//! through 31 after). If the count is at least `hamming_filter_id` the hit
//! survives. The high bit (soft-mask) is stripped before comparing — this
//! mirrors C++'s `letter_mask` on the fingerprint load path.
//!
//! In the C++ pipeline this filter eliminates roughly 90% of seed hits before
//! ungapped extension and is the primary mechanism that keeps DIAMOND's
//! default-mode output as selective as it is.
use crate::basic::value::{Letter, LETTER_MASK};
use crate::search::seed_match::SeedMatch;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

/// 48-letter fingerprint window matching C++ `FingerPrint::load`:
/// 16 letters before the seed anchor and 32 letters from the anchor onward.
const FP_BEFORE: usize = 16;
const FP_AFTER: usize = 32;
const FP_LEN: usize = FP_BEFORE + FP_AFTER;

/// Count exact-letter matches over the 48-position window around (q_pos, r_pos).
///
/// C++ `FingerPrint::load` (`search/hamming/finger_print.h:55-93`) unconditionally
/// reads 48 bytes via SIMD from `q-16` and `t-16` regardless of sequence
/// boundaries — the `SequenceSet` buffer holds `DELIMITER_LETTER` padding on
/// both sides, so out-of-range reads return the same byte on both sides and
/// the `match()` count effectively treats those positions as matches.
///
/// Rust holds raw `&[Letter]` slices with no padding, so we must simulate:
/// when **both** positions are out of range (before-start or past-end), the
/// padding-vs-padding equality counts as a match; when **only one** side is
/// out of range, the real-vs-padding comparison never matches and we don't
/// count it. The previous "skip if either side OOR" undercounted matches by
/// up to ~16 near sequence boundaries and silently dropped hits C++ keeps.
#[inline]
fn fingerprint_match(query: &[Letter], target: &[Letter], q_pos: usize, r_pos: usize) -> u32 {
    let mut count = 0u32;
    let qlen = query.len() as isize;
    let tlen = target.len() as isize;
    for i in 0..FP_LEN {
        let q_idx = q_pos as isize + i as isize - FP_BEFORE as isize;
        let r_idx = r_pos as isize + i as isize - FP_BEFORE as isize;
        let q_oor = q_idx < 0 || q_idx >= qlen;
        let r_oor = r_idx < 0 || r_idx >= tlen;
        if q_oor && r_oor {
            // Both sides read padding (DELIMITER_LETTER on both) → match.
            count += 1;
            continue;
        }
        if q_oor || r_oor {
            // One side is padding, the other is a real letter — never matches.
            continue;
        }
        if (query[q_idx as usize] & LETTER_MASK) == (target[r_idx as usize] & LETTER_MASK) {
            count += 1;
        }
    }
    count
}

/// Apply the stage-1 hamming filter in parallel. Returns the surviving matches.
pub fn apply_hamming_filter(
    matches: Vec<SeedMatch>,
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    hamming_filter_id: u32,
) -> Vec<SeedMatch> {
    matches
        .par_iter()
        .filter_map(|m| {
            let q = query_seqs[m.query_id as usize];
            let t = ref_seqs[m.ref_id as usize];
            if fingerprint_match(q, t, m.query_pos as usize, m.ref_pos as usize)
                >= hamming_filter_id
            {
                Some(*m)
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_window_full_match() {
        let q: Vec<Letter> = (0..64).map(|i| (i % 20) as Letter).collect();
        let r = q.clone();
        // Anchor in the middle so the whole 48-letter window is in-range.
        assert_eq!(fingerprint_match(&q, &r, 32, 32), FP_LEN as u32);
    }

    #[test]
    fn boundary_positions_zero_pad() {
        let q: Vec<Letter> = (0..40).map(|i| (i % 20) as Letter).collect();
        let r = q.clone();
        // q_pos=0 means positions 0..32 are real (matching) and positions
        // -16..0 are out-of-range on BOTH sides — now treated as padding-
        // vs-padding matches (DELIMITER vs DELIMITER) per C++. So the count
        // is the full 48: 16 phantom matches + 32 real matches.
        let n = fingerprint_match(&q, &r, 0, 0);
        assert_eq!(n, FP_LEN as u32);
    }

    #[test]
    fn one_sided_oor_does_not_match() {
        // q_pos=0, r_pos=32: for i=0..15, q_idx=-16..-1 (OOR) and r_idx=16..31
        // (in-range real letters). Padding-vs-real never matches. For i=16..47,
        // both sides are in-range but q[0..32] = i%20 and r[32..64] = (i+32)%20
        // — offset by 32 = 12 mod 20, so they never align letter-for-letter.
        // Total expected: 0.
        let q: Vec<Letter> = (0..40).map(|i| (i % 20) as Letter).collect();
        let r: Vec<Letter> = (0..80).map(|i| (i % 20) as Letter).collect();
        let n = fingerprint_match(&q, &r, 0, 32);
        assert_eq!(n, 0);
    }

    #[test]
    fn mismatching_window_returns_zero() {
        let q: Vec<Letter> = (0..64).map(|i| (i % 20) as Letter).collect();
        let r: Vec<Letter> = (0..64).map(|i| ((i + 1) % 20) as Letter).collect();
        assert_eq!(fingerprint_match(&q, &r, 32, 32), 0);
    }

    #[test]
    fn high_bit_is_stripped() {
        // soft-masked positions still count as matches when underlying letter agrees.
        let q: Vec<Letter> = (0..64).map(|i| (i % 20) as Letter).collect();
        let mut r = q.clone();
        for x in r.iter_mut() {
            *x |= 0x80u8 as Letter;
        }
        assert_eq!(fingerprint_match(&q, &r, 32, 32), FP_LEN as u32);
    }
}
