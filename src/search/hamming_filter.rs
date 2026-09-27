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
use crate::basic::value::{Letter, DELIMITER_LETTER, LETTER_MASK};
use crate::search::seed_match::SeedMatch;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::cell::RefCell;

/// 48-letter fingerprint window matching C++ `FingerPrint::load`:
/// 16 letters before the seed anchor and 32 letters from the anchor onward.
const FP_BEFORE: usize = 16;
const FP_AFTER: usize = 32;
const FP_LEN: usize = FP_BEFORE + FP_AFTER;

type CachedFingerprint = ((u32, u32), [Letter; FP_LEN]);

#[derive(Default)]
struct FingerprintScratch {
    query: Vec<CachedFingerprint>,
    target: Vec<CachedFingerprint>,
}

thread_local! {
    // C++ keeps `vq`/`vs` in each worker's WorkSet. Do the same here so the
    // many small seed partitions do not allocate two vectors apiece.
    static FINGERPRINT_SCRATCH: RefCell<FingerprintScratch> =
        RefCell::new(FingerprintScratch::default());
}

#[inline]
fn fingerprint_equal_count(query: &[Letter; FP_LEN], target: &[Letter; FP_LEN]) -> u32 {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: AVX2 was detected and both fixed arrays contain all 48
            // bytes loaded by the kernel.
            return unsafe { fingerprint_equal_count_avx2(query, target) };
        }
        if std::arch::is_x86_feature_detected!("sse2") {
            // SAFETY: SSE2 was detected and both fixed arrays contain all 48
            // bytes loaded by the kernel.
            return unsafe { fingerprint_equal_count_sse2(query, target) };
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        // NEON is mandatory for AArch64, and both arrays contain all 48 bytes
        // consumed by the three vector loads.
        return unsafe { fingerprint_equal_count_neon(query, target) };
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        query
            .iter()
            .zip(target)
            .map(|(&q, &t)| u32::from(q == t))
            .sum()
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn fingerprint_equal_count_neon(query: &[Letter; FP_LEN], target: &[Letter; FP_LEN]) -> u32 {
    use std::arch::aarch64::*;

    let ones = vdupq_n_u8(1);
    let mut sum = vdupq_n_u16(0);
    for offset in [0, 16, 32] {
        let q = vld1q_u8(query.as_ptr().add(offset).cast());
        let t = vld1q_u8(target.as_ptr().add(offset).cast());
        let equal = vandq_u8(vceqq_u8(q, t), ones);
        sum = vpadalq_u8(sum, equal);
    }
    u32::from(vaddvq_u16(sum))
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn fingerprint_equal_count_avx2(query: &[Letter; FP_LEN], target: &[Letter; FP_LEN]) -> u32 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let q0 = _mm256_loadu_si256(query.as_ptr().cast());
    let t0 = _mm256_loadu_si256(target.as_ptr().cast());
    let q1 = _mm_loadu_si128(query.as_ptr().add(32).cast());
    let t1 = _mm_loadu_si128(target.as_ptr().add(32).cast());
    let m0 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(q0, t0)) as u32;
    let m1 = _mm_movemask_epi8(_mm_cmpeq_epi8(q1, t1)) as u32;
    m0.count_ones() + m1.count_ones()
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn fingerprint_equal_count_sse2(query: &[Letter; FP_LEN], target: &[Letter; FP_LEN]) -> u32 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let mut count = 0;
    for offset in [0, 16, 32] {
        let q = _mm_loadu_si128(query.as_ptr().add(offset).cast());
        let t = _mm_loadu_si128(target.as_ptr().add(offset).cast());
        count += (_mm_movemask_epi8(_mm_cmpeq_epi8(q, t)) as u32).count_ones();
    }
    count
}

/// Count exact-letter matches over the 48-position window around (q_pos, r_pos).
///
/// C++ `FingerPrint::load` (`search/hamming/finger_print.h:55-93`) unconditionally
/// reads 48 bytes via SIMD from `q-16` and `t-16` regardless of sequence
/// boundaries. C++ `SequenceSet` places one delimiter between records rather
/// than filling the whole SIMD overread with delimiters, so positions outside
/// both logical records must not be treated as synthetic identities.
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
        if q_oor || r_oor {
            continue;
        }
        let q_letter = query[q_idx as usize] & LETTER_MASK;
        let r_letter = target[r_idx as usize] & LETTER_MASK;
        if q_letter == r_letter {
            count += 1;
        }
    }
    count
}

/// Apply the stage-1 hamming filter in parallel. Returns the surviving matches.
///
/// The predicate is evaluated into a one-byte-per-hit keep mask, then the
/// input allocation is compacted in place.  In particular, do not express
/// this as `par_iter().filter_map(...).collect()`: `SeedMatch` is 32 bytes and
/// that formulation keeps a second, survivor-sized match allocation alive
/// alongside the raw join output at the seed stage's memory high-water mark.
/// `Vec::retain` is stable, so shape/partition join emission order is kept.
pub fn apply_hamming_filter(
    mut matches: Vec<SeedMatch>,
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    hamming_filter_id: u32,
) -> Vec<SeedMatch> {
    let keep: Vec<u8> = matches
        .par_iter()
        .map(|m| {
            let q = query_seqs[m.query_id as usize];
            let t = ref_seqs[m.ref_id as usize];
            u8::from(
                fingerprint_match(q, t, m.query_pos as usize, m.ref_pos as usize)
                    >= hamming_filter_id,
            )
        })
        .collect();

    let mut i = 0usize;
    matches.retain(|_| {
        let retain = keep[i] != 0;
        i += 1;
        retain
    });
    matches
}

#[inline]
fn sequence_set_letter(
    seqs: &[&[Letter]],
    restores: &[Vec<(usize, Letter)>],
    mut seq_id: usize,
    mut pos: isize,
) -> Letter {
    loop {
        let len = seqs[seq_id].len() as isize;
        if pos < 0 {
            if pos == -1 || seq_id == 0 {
                return DELIMITER_LETTER;
            }
            seq_id -= 1;
            pos += seqs[seq_id].len() as isize + 1;
        } else if pos >= len {
            if pos == len || seq_id + 1 == seqs.len() {
                return DELIMITER_LETTER;
            }
            pos -= len + 1;
            seq_id += 1;
        } else {
            let pos = pos as usize;
            return restores[seq_id]
                .binary_search_by_key(&pos, |&(p, _)| p)
                .map_or(seqs[seq_id][pos] & LETTER_MASK, |i| {
                    restores[seq_id][i].1 & LETTER_MASK
                });
        }
    }
}

pub(crate) fn retain_hamming_filter_sequence_set(
    matches: &mut Vec<SeedMatch>,
    query_seqs: &[&[Letter]],
    query_restores: &[Vec<(usize, Letter)>],
    ref_seqs: &[&[Letter]],
    ref_restores: &[Vec<(usize, Letter)>],
    hamming_filter_id: u32,
) {
    const TILE_SIZE: usize = 64;

    #[inline]
    fn load(
        seqs: &[&[Letter]],
        restores: &[Vec<(usize, Letter)>],
        seq_id: u32,
        pos: u32,
    ) -> [Letter; FP_LEN] {
        std::array::from_fn(|i| {
            sequence_set_letter(
                seqs,
                restores,
                seq_id as usize,
                pos as isize + i as isize - FP_BEFORE as isize,
            )
        })
    }

    // `sort_merge_seed_matches*` emits one contiguous cross product for each
    // shared seed key, in q-major order. Mirror C++ `load_fps/all_vs_all`:
    // materialize every q/r fingerprint once per natural seed group and reuse
    // it for all comparisons. The two vectors retain capacity across groups,
    // so the hot path performs no per-hit allocation and its auxiliary RSS is
    // O(nq + nr), rather than O(nq * nr).
    let mut scratch = FINGERPRINT_SCRATCH.with(|slot| std::mem::take(&mut *slot.borrow_mut()));
    let FingerprintScratch {
        query: query_fps,
        target: target_fps,
    } = &mut scratch;
    let mut write = 0usize;
    let mut group_begin = 0usize;
    while group_begin < matches.len() {
        let seed = matches[group_begin].seed;
        let shape_id = matches[group_begin].shape_id;
        let mut group_end = group_begin + 1;
        while group_end < matches.len()
            && matches[group_end].seed == seed
            && matches[group_end].shape_id == shape_id
        {
            group_end += 1;
        }

        if group_end == group_begin + 1 {
            let m = matches[group_begin];
            let query_fp = load(query_seqs, query_restores, m.query_id, m.query_pos);
            let target_fp = load(ref_seqs, ref_restores, m.ref_id, m.ref_pos);
            if fingerprint_equal_count(&query_fp, &target_fp) >= hamming_filter_id {
                matches[write] = m;
                write += 1;
            }
            group_begin = group_end;
            continue;
        }

        query_fps.clear();
        target_fps.clear();
        let first_query = (
            matches[group_begin].query_id,
            matches[group_begin].query_pos,
        );
        let mut last_query = None;
        for m in &matches[group_begin..group_end] {
            let query = (m.query_id, m.query_pos);
            if last_query != Some(query) {
                query_fps.push((query, load(query_seqs, query_restores, query.0, query.1)));
                last_query = Some(query);
            }
            if query == first_query {
                let target = (m.ref_id, m.ref_pos);
                target_fps.push((target, load(ref_seqs, ref_restores, target.0, target.1)));
            }
        }

        // The production join guarantees this rectangular q-major layout.
        // Check that invariant in debug builds without paying another full
        // cross-product walk in optimized production builds.
        debug_assert_eq!(
            query_fps.len().saturating_mul(target_fps.len()),
            group_end - group_begin
        );
        debug_assert!(query_fps.iter().enumerate().all(|(qi, (query, _))| {
            target_fps.iter().enumerate().all(|(ri, (target, _))| {
                let m = matches[group_begin + qi * target_fps.len() + ri];
                (m.query_id, m.query_pos) == *query && (m.ref_id, m.ref_pos) == *target
            })
        }));

        let mut read = group_begin;
        for (_, query_fp) in query_fps.iter() {
            // Tiling keeps the active target fingerprints in L1 while
            // retaining q-major order, hence stable survivor order.
            for target_tile in target_fps.chunks(TILE_SIZE) {
                for (_, target_fp) in target_tile {
                    let m = matches[read];
                    read += 1;
                    if fingerprint_equal_count(query_fp, target_fp) >= hamming_filter_id {
                        matches[write] = m;
                        write += 1;
                    }
                }
            }
        }
        group_begin = group_end;
    }
    matches.truncate(write);
    FINGERPRINT_SCRATCH.with(|slot| *slot.borrow_mut() = scratch);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simd_fingerprint_count_matches_scalar_for_all_bytes() {
        let mut query = [0; FP_LEN];
        let mut target = [0; FP_LEN];
        for shift in 0..64 {
            for i in 0..FP_LEN {
                query[i] = ((i * 17 + shift) & 31) as Letter;
                target[i] = ((i * 11 + shift * 3) & 31) as Letter;
                if (i + shift) % 7 == 0 {
                    target[i] = query[i];
                }
            }
            let expected = query.iter().zip(&target).filter(|(q, t)| q == t).count() as u32;
            assert_eq!(fingerprint_equal_count(&query, &target), expected);
        }
    }

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
        // q_pos=0 means positions 0..32 are real and positions -16..0 are
        // outside both logical records. Only the 32 real positions count.
        let n = fingerprint_match(&q, &r, 0, 0);
        assert_eq!(n, FP_AFTER as u32);
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

    #[test]
    fn filter_compacts_in_place_and_preserves_order() {
        let query_match = vec![1; 64];
        let query_miss = vec![2; 64];
        let target = vec![1; 64];
        let query_seqs = vec![query_match.as_slice(), query_miss.as_slice()];
        let ref_seqs = vec![target.as_slice()];
        let matches = vec![
            SeedMatch {
                query_id: 0,
                query_pos: 24,
                ref_id: 0,
                ref_pos: 24,
                seed: 11,
                shape_id: 0,
            },
            SeedMatch {
                query_id: 1,
                query_pos: 24,
                ref_id: 0,
                ref_pos: 24,
                seed: 22,
                shape_id: 0,
            },
            SeedMatch {
                query_id: 0,
                query_pos: 25,
                ref_id: 0,
                ref_pos: 25,
                seed: 33,
                shape_id: 0,
            },
        ];
        let allocation = matches.as_ptr();

        let filtered = apply_hamming_filter(matches, &query_seqs, &ref_seqs, FP_LEN as u32);

        assert_eq!(filtered.as_ptr(), allocation);
        assert_eq!(
            filtered.iter().map(|m| m.seed).collect::<Vec<_>>(),
            [11, 33]
        );
    }

    #[test]
    fn sequence_set_lookup_matches_cpp_delimiters_and_restores() {
        let a = vec![1, 2, 3];
        let b = vec![4, 23, 6];
        let c = vec![7, 8, 9];
        let seqs = vec![a.as_slice(), b.as_slice(), c.as_slice()];
        let restores = vec![Vec::new(), vec![(1, 5)], Vec::new()];

        assert_eq!(
            sequence_set_letter(&seqs, &restores, 1, -1),
            DELIMITER_LETTER
        );
        assert_eq!(sequence_set_letter(&seqs, &restores, 1, -2), 3);
        assert_eq!(sequence_set_letter(&seqs, &restores, 1, 1), 5);
        assert_eq!(
            sequence_set_letter(&seqs, &restores, 1, 3),
            DELIMITER_LETTER
        );
        assert_eq!(sequence_set_letter(&seqs, &restores, 1, 4), 7);
        assert_eq!(
            sequence_set_letter(&seqs, &restores, 0, -2),
            DELIMITER_LETTER
        );
    }
}
