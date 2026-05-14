use crate::basic::packed_loc::PackedLoc;
use crate::basic::reduction::Reduction;
use crate::basic::shape::Shape;
use crate::basic::value::{letter_mask, Letter, SEED_MASK, TRUE_AA};
use crate::data::flags::{PackedLocId, SeedEncoding};
use crate::data::seed_histogram::SeedPartitionRange;
use crate::data::sequence_set::SequenceSet;
use crate::util::data_structures::DoubleArray;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeedStats {
    pub good_seed_positions: usize,
    pub low_complexity_seeds: usize,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MaskSeedsStats {
    pub seed_count: usize,
    pub masked_seed_count: usize,
    pub query_count: usize,
    pub target_count: usize,
}

pub trait SeedLocBytes {
    const ELEM_SIZE: usize;
    fn pos_from_bytes(bytes: &[u8]) -> u64;
}

impl SeedLocBytes for PackedLoc {
    const ELEM_SIZE: usize = 5;

    fn pos_from_bytes(bytes: &[u8]) -> u64 {
        let high = bytes[0] as u64;
        let low = u32::from_ne_bytes(bytes[1..5].try_into().unwrap()) as u64;
        (high << 32) | low
    }
}

impl SeedLocBytes for PackedLocId {
    const ELEM_SIZE: usize = 9;

    fn pos_from_bytes(bytes: &[u8]) -> u64 {
        let high = bytes[0] as u64;
        let low = u32::from_ne_bytes(bytes[1..5].try_into().unwrap()) as u64;
        (high << 32) | low
    }
}

/// Matches C++ `Search::seed_is_complex`.
pub fn seed_is_complex(seq: &[Letter], shape: &Shape, cut: f64, reduction: &Reduction) -> bool {
    let mut count = [0u32; TRUE_AA as usize];
    for i in 0..shape.weight as usize {
        let l = letter_mask(seq[shape.positions[i] as usize]);
        if l >= TRUE_AA as Letter {
            return false;
        }
        count[reduction.reduce(l) as usize] += 1;
    }
    let mut entropy = ln_factorial(shape.weight as u32);
    for i in 0..reduction.size() as usize {
        entropy -= ln_factorial(count[i]);
    }
    entropy >= cut
}

/// Matches C++ `Search::seed_is_complex_unreduced`.
pub fn seed_is_complex_unreduced(
    seq: &mut [Letter],
    shape: &Shape,
    cut: f64,
    mask_seeds: bool,
    stats: &mut SeedStats,
) -> bool {
    let mut count = [0u32; TRUE_AA as usize];
    for i in 0..shape.weight as usize {
        let l = letter_mask(seq[shape.positions[i] as usize]);
        if l >= TRUE_AA as Letter {
            if mask_seeds {
                seq[0] |= SEED_MASK;
            }
            return false;
        }
        count[l as usize] += 1;
    }
    stats.good_seed_positions += 1;
    let mut entropy = ln_factorial(shape.weight as u32);
    for i in 0..TRUE_AA as usize {
        entropy -= ln_factorial(count[i]);
    }
    if entropy < cut {
        if mask_seeds {
            seq[0] |= SEED_MASK;
        }
        stats.low_complexity_seeds += 1;
        return false;
    }
    true
}

/// Matches C++ `Search::mask_seeds`.
pub fn mask_seeds<SeedLoc>(
    shape: &Shape,
    range: &SeedPartitionRange,
    query_seed_hits: &mut [DoubleArray],
    ref_seed_hits: &mut [DoubleArray],
    query_seqs: &mut SequenceSet,
    seed_encoding: SeedEncoding,
    cut: f64,
    reduction: &Reduction,
) -> MaskSeedsStats
where
    SeedLoc: SeedLocBytes,
{
    if seed_encoding != SeedEncoding::SpacedFactor {
        return MaskSeedsStats::default();
    }

    let partition_count = range.size() as usize;
    let mut stats = MaskSeedsStats::default();
    for p in 0..partition_count {
        let mut q_it = query_seed_hits[p].begin_mut_with_elem_size(SeedLoc::ELEM_SIZE);
        let mut s_it = ref_seed_hits[p].begin_mut_with_elem_size(SeedLoc::ELEM_SIZE);
        while q_it.good() && s_it.good() {
            stats.seed_count += 1;
            let query_count = q_it.count() as usize;
            let target_count = s_it.count() as usize;
            let query_pos = {
                let query_hits = q_it.range_bytes(SeedLoc::ELEM_SIZE);
                SeedLoc::pos_from_bytes(&query_hits[..SeedLoc::ELEM_SIZE])
            };
            if !seed_is_complex(query_seqs.data_at(query_pos), shape, cut, reduction) {
                stats.masked_seed_count += 1;
                stats.query_count += query_count;
                stats.target_count += target_count;
                let query_hits = q_it.range_bytes(SeedLoc::ELEM_SIZE).to_vec();
                for hit in query_hits.chunks_exact(SeedLoc::ELEM_SIZE) {
                    let pos = SeedLoc::pos_from_bytes(hit);
                    query_seqs.data_mut_at(pos)[0] |= SEED_MASK;
                }
                q_it.erase(SeedLoc::ELEM_SIZE);
                s_it.erase(SeedLoc::ELEM_SIZE);
            } else {
                q_it.increment(SeedLoc::ELEM_SIZE);
                s_it.increment(SeedLoc::ELEM_SIZE);
            }
        }
    }
    stats
}

pub fn ln_factorial(n: u32) -> f64 {
    (2..=n).map(|i| (i as f64).ln()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ln_factorial() {
        assert_eq!(ln_factorial(0), 0.0);
        assert_eq!(ln_factorial(1), 0.0);
        assert!((ln_factorial(4) - (24.0f64).ln()).abs() < 1.0e-12);
    }

    #[test]
    fn test_seed_is_complex_reduced() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);
        assert!(!seed_is_complex(&[0, 0, 0, 0], &shape, 1.0, &reduction));
        assert!(seed_is_complex(&[0, 1, 3, 4], &shape, 1.0, &reduction));
        assert!(!seed_is_complex(&[0, 1, 3, 24], &shape, 0.0, &reduction));
    }

    #[test]
    fn test_seed_is_complex_unreduced_stats_and_mask() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);
        let mut stats = SeedStats::default();
        let mut low = [0, 0, 0, 0];
        assert!(!seed_is_complex_unreduced(
            &mut low, &shape, 1.0, true, &mut stats
        ));
        assert_eq!(stats.good_seed_positions, 1);
        assert_eq!(stats.low_complexity_seeds, 1);
        assert_ne!(low[0] & SEED_MASK, 0);

        let mut good = [0, 1, 2, 3];
        assert!(seed_is_complex_unreduced(
            &mut good, &shape, 1.0, true, &mut stats
        ));
        assert_eq!(stats.good_seed_positions, 2);

        let mut invalid = [0, 1, 2, 24];
        assert!(!seed_is_complex_unreduced(
            &mut invalid,
            &shape,
            0.0,
            true,
            &mut stats
        ));
        assert_ne!(invalid[0] & SEED_MASK, 0);
        assert_eq!(stats.good_seed_positions, 2);
    }

    #[test]
    fn test_mask_seeds_packed_loc_erases_low_complexity_join_group() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 0, 0, 0, 1, 2, 3, 4]);
        let low_pos = seqs.position(0, 0) as u64;
        let good_pos = seqs.position(0, 4) as u64;

        let mut query_bytes = Vec::new();
        query_bytes.extend_from_slice(&1u32.to_ne_bytes());
        let low = PackedLoc::new(low_pos);
        query_bytes.push(low.high);
        query_bytes.extend_from_slice(&low.low.to_ne_bytes());
        query_bytes.extend_from_slice(&1u32.to_ne_bytes());
        let good = PackedLoc::new(good_pos);
        query_bytes.push(good.high);
        query_bytes.extend_from_slice(&good.low.to_ne_bytes());

        let mut ref_bytes = Vec::new();
        ref_bytes.extend_from_slice(&2u32.to_ne_bytes());
        for pos in [20, 25] {
            let loc = PackedLoc::new(pos);
            ref_bytes.push(loc.high);
            ref_bytes.extend_from_slice(&loc.low.to_ne_bytes());
        }
        ref_bytes.extend_from_slice(&1u32.to_ne_bytes());
        let loc = PackedLoc::new(30);
        ref_bytes.push(loc.high);
        ref_bytes.extend_from_slice(&loc.low.to_ne_bytes());

        let mut query_hits = [DoubleArray::from_bytes(query_bytes, 18)];
        let mut ref_hits = [DoubleArray::from_bytes(ref_bytes, 23)];
        let stats = mask_seeds::<PackedLoc>(
            &shape,
            &SeedPartitionRange::with_bounds(0, 1),
            &mut query_hits,
            &mut ref_hits,
            &mut seqs,
            SeedEncoding::SpacedFactor,
            1.0,
            &reduction,
        );

        assert_eq!(
            stats,
            MaskSeedsStats {
                seed_count: 2,
                masked_seed_count: 1,
                query_count: 1,
                target_count: 2,
            }
        );
        assert_ne!(seqs.data_at(low_pos)[0] & SEED_MASK, 0);
        assert_eq!(seqs.data_at(good_pos)[0] & SEED_MASK, 0);
        assert_eq!(query_hits[0].begin_with_elem_size(5).count(), 1);
        assert_eq!(ref_hits[0].begin_with_elem_size(5).count(), 1);
    }

    #[test]
    fn test_mask_seeds_packed_loc_id_uses_nine_byte_entries() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("1111", &reduction);
        let mut seqs = SequenceSet::new();
        seqs.push(&[0, 0, 0, 0]);
        let pos = seqs.position(0, 0) as u64;
        let loc = PackedLoc::new(pos);

        let mut query_bytes = Vec::new();
        query_bytes.extend_from_slice(&1u32.to_ne_bytes());
        query_bytes.push(loc.high);
        query_bytes.extend_from_slice(&loc.low.to_ne_bytes());
        query_bytes.extend_from_slice(&7u32.to_ne_bytes());

        let mut ref_bytes = Vec::new();
        ref_bytes.extend_from_slice(&1u32.to_ne_bytes());
        ref_bytes.push(loc.high);
        ref_bytes.extend_from_slice(&loc.low.to_ne_bytes());
        ref_bytes.extend_from_slice(&11u32.to_ne_bytes());

        let mut query_hits = [DoubleArray::from_bytes(query_bytes, 13)];
        let mut ref_hits = [DoubleArray::from_bytes(ref_bytes, 13)];
        let stats = mask_seeds::<PackedLocId>(
            &shape,
            &SeedPartitionRange::with_bounds(0, 1),
            &mut query_hits,
            &mut ref_hits,
            &mut seqs,
            SeedEncoding::SpacedFactor,
            1.0,
            &reduction,
        );

        assert_eq!(stats.masked_seed_count, 1);
        assert!(!query_hits[0].begin_with_elem_size(9).good());
        assert!(!ref_hits[0].begin_with_elem_size(9).good());
        assert_ne!(seqs.data_at(pos)[0] & SEED_MASK, 0);
    }
}
