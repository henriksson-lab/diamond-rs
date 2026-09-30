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
use crate::search::hamming::FingerPrint;
use crate::search::hamming_all_vs_all::AlignedFingerprint48;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::search::hamming_all_vs_all::{
    all_vs_all_pass_masks_avx2, all_vs_all_pass_masks_avx512bw,
};
use crate::search::seed_match::SeedMatch;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::cell::RefCell;

/// 48-letter fingerprint window matching C++ `FingerPrint::load`:
/// 16 letters before the seed anchor and 32 letters from the anchor onward.
const FP_BEFORE: usize = 16;
const FP_AFTER: usize = 32;
const FP_LEN: usize = FP_BEFORE + FP_AFTER;

type CachedFingerprint = ((u32, u32), [Letter; FP_LEN]);

/// C++ stores fingerprints as a separate `vector<array<char, 48>>`. Keeping
/// locations out of this hot array gives the same 48-byte stride (rather than
/// 56 bytes for `(location, fingerprint)`) and guarantees the alignment used
/// by its SSE/AVX loads.
#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
pub(crate) const FP_SCALAR: u8 = 0;
pub(crate) const FP_SSE2: u8 = 1;
pub(crate) const FP_AVX2: u8 = 2;
pub(crate) const FP_AVX512BW: u8 = 3;
pub(crate) const FP_NEON: u8 = 4;

#[derive(Default)]
struct FingerprintScratch {
    query: Vec<CachedFingerprint>,
    target: Vec<CachedFingerprint>,
}

#[derive(Default)]
struct StreamingFingerprintScratch {
    query: Vec<AlignedFingerprint48>,
    target: Vec<AlignedFingerprint48>,
}

thread_local! {
    // C++ keeps `vq`/`vs` in each worker's WorkSet. Do the same here so the
    // many small seed partitions do not allocate two vectors apiece.
    static FINGERPRINT_SCRATCH: RefCell<FingerprintScratch> =
        RefCell::new(FingerprintScratch::default());
    static STREAMING_FINGERPRINT_SCRATCH: RefCell<StreamingFingerprintScratch> =
        RefCell::new(StreamingFingerprintScratch::default());
}

#[inline]
fn fingerprint_equal_count(query: &[Letter; FP_LEN], target: &[Letter; FP_LEN]) -> u32 {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx512bw") {
            // SAFETY: AVX-512BW was detected; the masked loads touch exactly
            // the 48 initialized bytes in each fingerprint.
            return unsafe { fingerprint_equal_count_avx512(query, target) };
        }
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

/// Compile-time-selected fingerprint comparator. Runtime CPU dispatch happens
/// at the seed-group boundary so the all-vs-all inner loop contains no feature
/// test. Keeping this specialization local avoids cloning the much larger
/// stage-1 pipeline for every supported ISA.
#[inline(always)]
fn fingerprint_equal_count_for<const KERNEL: u8>(
    query: &[Letter; FP_LEN],
    target: &[Letter; FP_LEN],
) -> u32 {
    match KERNEL {
        FP_AVX512BW => {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            // SAFETY: this instantiation is selected only after AVX-512BW
            // detection at the group boundary.
            unsafe {
                return fingerprint_equal_count_avx512(query, target);
            }
        }
        FP_AVX2 => {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            // SAFETY: this instantiation is selected only after AVX2
            // detection at the group boundary.
            unsafe {
                return fingerprint_equal_count_avx2(query, target);
            }
        }
        FP_SSE2 => {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            // SAFETY: this instantiation is selected only after SSE2
            // detection at the group boundary.
            unsafe {
                return fingerprint_equal_count_sse2(query, target);
            }
        }
        FP_NEON => {
            #[cfg(target_arch = "aarch64")]
            // SAFETY: NEON is mandatory on AArch64.
            unsafe {
                return fingerprint_equal_count_neon(query, target);
            }
        }
        _ => {}
    }
    query
        .iter()
        .zip(target)
        .map(|(&q, &t)| u32::from(q == t))
        .sum()
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
unsafe fn fingerprint_equal_count_avx512(
    query: &[Letter; FP_LEN],
    target: &[Letter; FP_LEN],
) -> u32 {
    let valid = (1u64 << FP_LEN) - 1;
    let bits: u64;
    std::arch::asm!(
        "kmovq k1, {valid}",
        "vmovdqu8 zmm0 {{k1}}{{z}}, [{query}]",
        "vmovdqu8 zmm1 {{k1}}{{z}}, [{target}]",
        "vpcmpeqb k2 {{k1}}, zmm0, zmm1",
        "kmovq {bits}, k2",
        valid = in(reg) valid,
        query = in(reg) query.as_ptr(),
        target = in(reg) target.as_ptr(),
        bits = lateout(reg) bits,
        out("zmm0") _,
        out("zmm1") _,
        out("k1") _,
        out("k2") _,
        options(readonly, nostack, preserves_flags),
    );
    bits.count_ones()
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
fn sequence_set_letter(seqs: &[&[Letter]], mut seq_id: usize, mut pos: isize) -> Letter {
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
            return seqs[seq_id][pos as usize] & LETTER_MASK;
        }
    }
}

#[inline]
fn load_sequence_set_fingerprint(seqs: &[&[Letter]], seq_id: u32, pos: u32) -> [Letter; FP_LEN] {
    let seq = seqs[seq_id as usize];
    let center = pos as usize;
    if center >= FP_BEFORE && center + FP_AFTER <= seq.len() {
        // The production caller supplies the original, pre-motif sequence.
        // Consequently every interior window is contiguous and can follow
        // C++ `FingerPrint::load` directly without consulting a sparse motif
        // restoration table.
        return FingerPrint::from_seq_center(seq, center).r;
    }
    std::array::from_fn(|i| {
        sequence_set_letter(
            seqs,
            seq_id as usize,
            pos as isize + i as isize - FP_BEFORE as isize,
        )
    })
}

// Compatibility path for the materialized stage-1 helpers. Production blastp
// uses `visit_hamming_group` with original sequence backing and never enters
// this sparse-restoration path.
#[inline]
fn sequence_set_letter_restored(
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

#[inline]
fn load_sequence_set_fingerprint_restored(
    seqs: &[&[Letter]],
    restores: &[Vec<(usize, Letter)>],
    seq_id: u32,
    pos: u32,
) -> [Letter; FP_LEN] {
    let seq = seqs[seq_id as usize];
    let restore = &restores[seq_id as usize];
    let center = pos as usize;
    if center >= FP_BEFORE && center + FP_AFTER <= seq.len() {
        let begin = center - FP_BEFORE;
        let end = center + FP_AFTER;
        let first_restore = restore.partition_point(|&(restore_pos, _)| restore_pos < begin);
        if restore
            .get(first_restore)
            .is_none_or(|&(restore_pos, _)| restore_pos >= end)
        {
            return FingerPrint::from_seq_center(seq, center).r;
        }
    }
    std::array::from_fn(|i| {
        sequence_set_letter_restored(
            seqs,
            restores,
            seq_id as usize,
            pos as isize + i as isize - FP_BEFORE as isize,
        )
    })
}

/// Visit a joined seed group's passing q×r pairs without materializing that
/// cross product. Fingerprints are loaded once per location and reused, as in
/// C++ `all_vs_all`.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn visit_hamming_group<F>(
    query_locs: &[(u32, u32)],
    target_locs: &[(u32, u32)],
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    hamming_filter_id: u32,
    visit: F,
) where
    F: FnMut((u32, u32), (u32, u32), usize),
{
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx512bw") {
            return visit_hamming_group_for::<FP_AVX512BW, F>(
                query_locs,
                target_locs,
                query_seqs,
                ref_seqs,
                hamming_filter_id,
                visit,
            );
        }
        if std::arch::is_x86_feature_detected!("avx2") {
            return visit_hamming_group_for::<FP_AVX2, F>(
                query_locs,
                target_locs,
                query_seqs,
                ref_seqs,
                hamming_filter_id,
                visit,
            );
        }
        if std::arch::is_x86_feature_detected!("sse2") {
            return visit_hamming_group_for::<FP_SSE2, F>(
                query_locs,
                target_locs,
                query_seqs,
                ref_seqs,
                hamming_filter_id,
                visit,
            );
        }
    }
    #[cfg(target_arch = "aarch64")]
    return visit_hamming_group_for::<FP_NEON, F>(
        query_locs,
        target_locs,
        query_seqs,
        ref_seqs,
        hamming_filter_id,
        visit,
    );
    #[cfg(not(target_arch = "aarch64"))]
    visit_hamming_group_for::<FP_SCALAR, F>(
        query_locs,
        target_locs,
        query_seqs,
        ref_seqs,
        hamming_filter_id,
        visit,
    )
}

#[inline]
pub(crate) fn visit_hamming_group_for<const KERNEL: u8, F>(
    query_locs: &[(u32, u32)],
    target_locs: &[(u32, u32)],
    query_seqs: &[&[Letter]],
    ref_seqs: &[&[Letter]],
    hamming_filter_id: u32,
    mut visit: F,
) where
    F: FnMut((u32, u32), (u32, u32), usize),
{
    const TILE_SIZE: usize = 64;
    let mut scratch =
        STREAMING_FINGERPRINT_SCRATCH.with(|slot| std::mem::take(&mut *slot.borrow_mut()));
    scratch.query.clear();
    scratch.target.clear();
    scratch.query.extend(query_locs.iter().map(|&(id, pos)| {
        AlignedFingerprint48::new(load_sequence_set_fingerprint(query_seqs, id, pos))
    }));
    scratch.target.extend(target_locs.iter().map(|&(id, pos)| {
        AlignedFingerprint48::new(load_sequence_set_fingerprint(ref_seqs, id, pos))
    }));
    let mut pass_masks = [0u64; TILE_SIZE];
    for (query_loc_tile, query_fp_tile) in query_locs
        .chunks(TILE_SIZE)
        .zip(scratch.query.chunks(TILE_SIZE))
    {
        for (target_tile_index, (target_loc_tile, target_fp_tile)) in target_locs
            .chunks(TILE_SIZE)
            .zip(scratch.target.chunks(TILE_SIZE))
            .enumerate()
        {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            if KERNEL == FP_AVX512BW {
                // SAFETY: this specialization is selected only after
                // AVX-512BW detection at the group boundary. Masked loads
                // consume exactly the 48 initialized fingerprint bytes.
                unsafe {
                    all_vs_all_pass_masks_avx512bw(
                        query_fp_tile,
                        target_fp_tile,
                        hamming_filter_id,
                        &mut pass_masks[..query_fp_tile.len()],
                    );
                }
            } else if KERNEL == FP_AVX2 {
                // SAFETY: AVX2 was detected at the group boundary. Inputs
                // contain complete fingerprints, target tiles have at most 64
                // rows, and the output covers every query in this tile.
                unsafe {
                    all_vs_all_pass_masks_avx2(
                        query_fp_tile,
                        target_fp_tile,
                        hamming_filter_id,
                        &mut pass_masks[..query_fp_tile.len()],
                    );
                }
            } else {
                pass_masks[..query_fp_tile.len()].fill(0);
                for (query_index, query_fp) in query_fp_tile.iter().enumerate() {
                    for (target_index, target_fp) in target_fp_tile.iter().enumerate() {
                        if fingerprint_equal_count_for::<KERNEL>(&query_fp.bytes, &target_fp.bytes)
                            >= hamming_filter_id
                        {
                            pass_masks[query_index] |= 1u64 << target_index;
                        }
                    }
                }
            }
            #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
            {
                pass_masks[..query_fp_tile.len()].fill(0);
                for (query_index, query_fp) in query_fp_tile.iter().enumerate() {
                    for (target_index, target_fp) in target_fp_tile.iter().enumerate() {
                        if fingerprint_equal_count_for::<KERNEL>(&query_fp.bytes, &target_fp.bytes)
                            >= hamming_filter_id
                        {
                            pass_masks[query_index] |= 1u64 << target_index;
                        }
                    }
                }
            }
            // Upstream writes a HitField and consumes it query-major after the
            // comparison kernel. These masks reproduce that ordering without
            // allocating a cross-product-sized match vector.
            for (query_index, &query_loc) in query_loc_tile.iter().enumerate() {
                let mut targets = pass_masks[query_index];
                while targets != 0 {
                    let target_index = targets.trailing_zeros() as usize;
                    visit(
                        query_loc,
                        target_loc_tile[target_index],
                        target_tile_index * TILE_SIZE + target_index,
                    );
                    targets &= targets - 1;
                }
            }
        }
    }
    STREAMING_FINGERPRINT_SCRATCH.with(|slot| *slot.borrow_mut() = scratch);
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
            let query_fp = load_sequence_set_fingerprint_restored(
                query_seqs,
                query_restores,
                m.query_id,
                m.query_pos,
            );
            let target_fp =
                load_sequence_set_fingerprint_restored(ref_seqs, ref_restores, m.ref_id, m.ref_pos);
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
                query_fps.push((
                    query,
                    load_sequence_set_fingerprint_restored(
                        query_seqs,
                        query_restores,
                        query.0,
                        query.1,
                    ),
                ));
                last_query = Some(query);
            }
            if query == first_query {
                let target = (m.ref_id, m.ref_pos);
                target_fps.push((
                    target,
                    load_sequence_set_fingerprint_restored(
                        ref_seqs,
                        ref_restores,
                        target.0,
                        target.1,
                    ),
                ));
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
    fn streaming_group_matches_materialized_filter() {
        let query_a: Vec<Letter> = (0..96).map(|i| (i % 20) as Letter).collect();
        let query_b: Vec<Letter> = (0..96).map(|i| ((i + 3) % 20) as Letter).collect();
        let target_a = query_a.clone();
        let target_b: Vec<Letter> = (0..96).map(|i| ((i + 7) % 20) as Letter).collect();
        let queries: Vec<&[Letter]> = vec![&query_a, &query_b];
        let targets: Vec<&[Letter]> = vec![&target_a, &target_b];
        let query_locs = vec![(0, 32), (1, 40)];
        let target_locs = vec![(0, 32), (1, 40)];
        let mut expected = Vec::new();
        for &query in &query_locs {
            for &target in &target_locs {
                if fingerprint_match(
                    queries[query.0 as usize],
                    targets[target.0 as usize],
                    query.1 as usize,
                    target.1 as usize,
                ) >= 20
                {
                    expected.push((query, target));
                }
            }
        }
        let mut actual = Vec::new();
        visit_hamming_group(
            &query_locs,
            &target_locs,
            &queries,
            &targets,
            20,
            |query, target, _| actual.push((query, target)),
        );
        assert_eq!(actual, expected);
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
            sequence_set_letter_restored(&seqs, &restores, 1, -1),
            DELIMITER_LETTER
        );
        assert_eq!(sequence_set_letter_restored(&seqs, &restores, 1, -2), 3);
        assert_eq!(sequence_set_letter_restored(&seqs, &restores, 1, 1), 5);
        assert_eq!(sequence_set_letter(&seqs, 1, 1), 23);
        assert_eq!(
            sequence_set_letter_restored(&seqs, &restores, 1, 3),
            DELIMITER_LETTER
        );
        assert_eq!(sequence_set_letter_restored(&seqs, &restores, 1, 4), 7);
        assert_eq!(
            sequence_set_letter_restored(&seqs, &restores, 0, -2),
            DELIMITER_LETTER
        );
    }
}
