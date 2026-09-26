use crate::basic::packed_loc::PackedLoc;
use crate::basic::value::letter_mask;
use crate::data::flags::PackedLocId;
use crate::data::seed_histogram::SeedPartitionRange;
use crate::data::sequence_set::SequenceSet;
use crate::search::seed_complexity::SeedLocBytes;
use crate::util::data_structures::DoubleArray;
use crate::util::misc::Sd;
use rayon::prelude::*;

pub struct FrequentSeeds;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrequentSeedsConfig {
    pub freq_sd: f64,
    pub threads: usize,
}

impl Default for FrequentSeedsConfig {
    fn default() -> Self {
        Self {
            freq_sd: 0.0,
            threads: 1,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FrequentSeedsBuildStats {
    pub ref_mean: f64,
    pub ref_sd: f64,
    pub query_mean: f64,
    pub query_sd: f64,
    pub ref_max_n: u32,
    pub query_max_n: u32,
    pub masked_positions: u32,
}

impl FrequentSeeds {
    pub const HASH_TABLE_FACTOR: f64 = 1.3;

    /// Matches C++ `compute_sd(range, query_seed_hits, ref_seed_hits)`.
    pub fn compute_sd<SeedLoc>(
        range: &SeedPartitionRange,
        query_seed_hits: &[DoubleArray],
        ref_seed_hits: &[DoubleArray],
    ) -> (Vec<Sd>, Vec<Sd>)
    where
        SeedLoc: SeedLocBytes,
    {
        let partition_count = range.size() as usize;
        let mut ref_out = vec![Sd::new(); partition_count];
        let mut query_out = vec![Sd::new(); partition_count];
        for p in 0..partition_count {
            let mut ref_sd = Sd::new();
            let mut query_sd = Sd::new();
            let mut q_it = query_seed_hits[p].begin_with_elem_size(SeedLoc::ELEM_SIZE);
            let mut s_it = ref_seed_hits[p].begin_with_elem_size(SeedLoc::ELEM_SIZE);
            while q_it.good() && s_it.good() {
                query_sd.add(q_it.count() as f64);
                ref_sd.add(s_it.count() as f64);
                q_it.increment(SeedLoc::ELEM_SIZE);
                s_it.increment(SeedLoc::ELEM_SIZE);
            }
            ref_out[p] = ref_sd;
            query_out[p] = query_sd;
        }
        (ref_out, query_out)
    }

    /// Matches C++ `FrequentSeeds::build_worker(rel_seedp, query_seed_hits, ref_seed_hits, query_seqs, ref_max_n, query_max_n, counts)`.
    pub fn build_worker<SeedLoc>(
        rel_seedp: usize,
        query_seed_hits: &mut [DoubleArray],
        ref_seed_hits: &mut [DoubleArray],
        query_seqs: &mut SequenceSet,
        ref_max_n: u32,
        query_max_n: u32,
        counts: &mut [u32],
    ) where
        SeedLoc: SeedLocBytes,
    {
        let (n, positions) = filter_partition::<SeedLoc>(
            &mut query_seed_hits[rel_seedp],
            &mut ref_seed_hits[rel_seedp],
            ref_max_n,
            query_max_n,
        );
        for pos in positions {
            query_seqs.data_mut_at(pos)[0] |= crate::basic::value::SEED_MASK;
        }
        counts[rel_seedp] = n;
    }

    /// C++ `build` with thread count and frequency cutoff made explicit.
    pub fn build_with_config<SeedLoc>(
        _sid: u32,
        range: &SeedPartitionRange,
        query_seed_hits: &mut [DoubleArray],
        ref_seed_hits: &mut [DoubleArray],
        query_seqs: &mut SequenceSet,
        config: FrequentSeedsConfig,
    ) -> Result<FrequentSeedsBuildStats, String>
    where
        SeedLoc: SeedLocBytes + Send + Sync,
    {
        let partition_count = range.size() as usize;
        if config.threads == 0 {
            return Err("Frequent seed thread count must be positive.".into());
        }
        if query_seed_hits.len() < partition_count || ref_seed_hits.len() < partition_count {
            return Err("Seed hit arrays do not cover the partition range.".into());
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()
            .map_err(|e| e.to_string())?;

        let groups: Vec<(Sd, Sd)> = pool.install(|| {
            (0..partition_count)
                .into_par_iter()
                .map(|p| partition_sd::<SeedLoc>(&query_seed_hits[p], &ref_seed_hits[p]))
                .collect()
        });
        let ref_sds = groups.iter().map(|(r, _)| *r).collect::<Vec<_>>();
        let query_sds = groups.iter().map(|(_, q)| *q).collect::<Vec<_>>();
        let ref_sd = Sd::from_groups(&ref_sds);
        let query_sd = Sd::from_groups(&query_sds);
        let ref_max_n = (ref_sd.mean() + config.freq_sd * ref_sd.sd()) as u32;
        let query_max_n = (query_sd.mean() + config.freq_sd * query_sd.sd()) as u32;

        let filtered: Vec<(u32, Vec<u64>)> = pool.install(|| {
            query_seed_hits[..partition_count]
                .par_iter_mut()
                .zip(ref_seed_hits[..partition_count].par_iter_mut())
                .map(|(query, reference)| {
                    filter_partition::<SeedLoc>(query, reference, ref_max_n, query_max_n)
                })
                .collect()
        });
        let mut masked_positions = 0u32;
        for (count, positions) in filtered {
            masked_positions = masked_positions.wrapping_add(count);
            // Applying the collected masks serially avoids the C++ data race
            // when the same query position occurs in multiple partitions.
            for pos in positions {
                query_seqs.data_mut_at(pos)[0] |= crate::basic::value::SEED_MASK;
            }
        }
        Ok(FrequentSeedsBuildStats {
            ref_mean: ref_sd.mean(),
            ref_sd: ref_sd.sd(),
            query_mean: query_sd.mean(),
            query_sd: query_sd.sd(),
            ref_max_n,
            query_max_n,
            masked_positions,
        })
    }

    /// Matches C++ `FrequentSeeds::build(sid, range, query_seed_hits, ref_seed_hits, query_seqs, freq_sd)`.
    pub fn build<SeedLoc>(
        _sid: u32,
        range: &SeedPartitionRange,
        query_seed_hits: &mut [DoubleArray],
        ref_seed_hits: &mut [DoubleArray],
        query_seqs: &mut SequenceSet,
        freq_sd: f64,
    ) -> FrequentSeedsBuildStats
    where
        SeedLoc: SeedLocBytes + Send + Sync,
    {
        Self::build_with_config::<SeedLoc>(
            _sid,
            range,
            query_seed_hits,
            ref_seed_hits,
            query_seqs,
            FrequentSeedsConfig {
                freq_sd,
                threads: 1,
            },
        )
        .expect("single-thread frequent seed configuration is valid")
    }

    pub fn build_packed_loc(
        sid: u32,
        range: &SeedPartitionRange,
        query_seed_hits: &mut [DoubleArray],
        ref_seed_hits: &mut [DoubleArray],
        query_seqs: &mut SequenceSet,
        freq_sd: f64,
    ) -> FrequentSeedsBuildStats {
        Self::build::<PackedLoc>(
            sid,
            range,
            query_seed_hits,
            ref_seed_hits,
            query_seqs,
            freq_sd,
        )
    }

    pub fn build_packed_loc_id(
        sid: u32,
        range: &SeedPartitionRange,
        query_seed_hits: &mut [DoubleArray],
        ref_seed_hits: &mut [DoubleArray],
        query_seqs: &mut SequenceSet,
        freq_sd: f64,
    ) -> FrequentSeedsBuildStats {
        Self::build::<PackedLocId>(
            sid,
            range,
            query_seed_hits,
            ref_seed_hits,
            query_seqs,
            freq_sd,
        )
    }

    /// Matches C++ `FrequentSeeds::clear_masking(seqs)`.
    pub fn clear_masking(seqs: &mut SequenceSet) {
        for i in 0..seqs.len() {
            let seq = seqs.get_mut(i);
            for letter in seq {
                *letter = letter_mask(*letter);
            }
        }
    }
}

fn partition_sd<SeedLoc>(query: &DoubleArray, reference: &DoubleArray) -> (Sd, Sd)
where
    SeedLoc: SeedLocBytes,
{
    let mut ref_sd = Sd::new();
    let mut query_sd = Sd::new();
    let mut q_it = query.begin_with_elem_size(SeedLoc::ELEM_SIZE);
    let mut s_it = reference.begin_with_elem_size(SeedLoc::ELEM_SIZE);
    while q_it.good() && s_it.good() {
        query_sd.add(q_it.count() as f64);
        ref_sd.add(s_it.count() as f64);
        q_it.increment(SeedLoc::ELEM_SIZE);
        s_it.increment(SeedLoc::ELEM_SIZE);
    }
    (ref_sd, query_sd)
}

fn filter_partition<SeedLoc>(
    query: &mut DoubleArray,
    reference: &mut DoubleArray,
    ref_max_n: u32,
    query_max_n: u32,
) -> (u32, Vec<u64>)
where
    SeedLoc: SeedLocBytes,
{
    let mut count = 0u32;
    let mut positions = Vec::new();
    let mut q_it = query.begin_mut_with_elem_size(SeedLoc::ELEM_SIZE);
    let mut s_it = reference.begin_mut_with_elem_size(SeedLoc::ELEM_SIZE);
    while q_it.good() && s_it.good() {
        if s_it.count() > ref_max_n || q_it.count() > query_max_n {
            count = count.wrapping_add(s_it.count());
            positions.extend(
                q_it.range_bytes(SeedLoc::ELEM_SIZE)
                    .chunks_exact(SeedLoc::ELEM_SIZE)
                    .map(SeedLoc::pos_from_bytes),
            );
            q_it.erase(SeedLoc::ELEM_SIZE);
            s_it.erase(SeedLoc::ELEM_SIZE);
        } else {
            q_it.increment(SeedLoc::ELEM_SIZE);
            s_it.increment(SeedLoc::ELEM_SIZE);
        }
    }
    (count, positions)
}

pub static FREQUENT_SEEDS: FrequentSeeds = FrequentSeeds;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;

    fn push_packed_loc_group(out: &mut Vec<u8>, locs: &[PackedLoc]) {
        out.extend_from_slice(&(locs.len() as u32).to_ne_bytes());
        for loc in locs {
            out.push(loc.high);
            let low = loc.low;
            out.extend_from_slice(&low.to_ne_bytes());
        }
    }

    fn push_packed_loc_id_group(out: &mut Vec<u8>, locs: &[(PackedLoc, u32)]) {
        out.extend_from_slice(&(locs.len() as u32).to_ne_bytes());
        for (loc, block_id) in locs {
            out.push(loc.high);
            let low = loc.low;
            out.extend_from_slice(&low.to_ne_bytes());
            out.extend_from_slice(&block_id.to_ne_bytes());
        }
    }

    #[test]
    fn test_clear_masking() {
        let mut seqs = SequenceSet::new();
        seqs.push(&[SEED_MASK | 1, 2, SEED_MASK | 3]);
        seqs.push(&[4, SEED_MASK | 5]);

        FrequentSeeds::clear_masking(&mut seqs);

        assert_eq!(seqs.get(0), &[1, 2, 3]);
        assert_eq!(seqs.get(1), &[4, 5]);
        assert_eq!(FrequentSeeds::HASH_TABLE_FACTOR, 1.3);
        let _ = FREQUENT_SEEDS;
    }

    #[test]
    fn test_build_packed_loc_masks_frequent_join_groups() {
        let mut seqs = SequenceSet::new();
        seqs.push(&[1, 2, 3, 4, 5, 6, 7, 8]);
        let p0 = PackedLoc::new(seqs.position(0, 1) as u64);
        let p1 = PackedLoc::new(seqs.position(0, 3) as u64);

        let mut query_bytes = Vec::new();
        push_packed_loc_group(&mut query_bytes, &[p0]);
        push_packed_loc_group(&mut query_bytes, &[p1]);

        let mut ref_bytes = Vec::new();
        push_packed_loc_group(
            &mut ref_bytes,
            &[
                PackedLoc::new(10),
                PackedLoc::new(11),
                PackedLoc::new(12),
                PackedLoc::new(13),
            ],
        );
        push_packed_loc_group(&mut ref_bytes, &[PackedLoc::new(20)]);

        let mut query_hits = [DoubleArray::from_bytes(query_bytes, 18)];
        let mut ref_hits = [DoubleArray::from_bytes(ref_bytes, 33)];
        let stats = FrequentSeeds::build_packed_loc(
            0,
            &SeedPartitionRange::with_bounds(0, 1),
            &mut query_hits,
            &mut ref_hits,
            &mut seqs,
            0.0,
        );

        assert_eq!(stats.ref_max_n, 2);
        assert_eq!(stats.query_max_n, 1);
        assert_eq!(stats.masked_positions, 4);
        assert_ne!(seqs.data_at(p0.as_u64())[0] & SEED_MASK, 0);
        assert_eq!(seqs.data_at(p1.as_u64())[0] & SEED_MASK, 0);
        assert_eq!(query_hits[0].begin_with_elem_size(5).count(), 1);
        assert_eq!(ref_hits[0].begin_with_elem_size(5).count(), 1);
    }

    #[test]
    fn test_build_packed_loc_id_masks_nine_byte_groups() {
        let mut seqs = SequenceSet::new();
        seqs.push(&[1, 2, 3, 4]);
        let p0 = PackedLoc::new(seqs.position(0, 1) as u64);

        let mut query_bytes = Vec::new();
        push_packed_loc_id_group(&mut query_bytes, &[(p0, 7), (PackedLoc::new(3), 8)]);

        let mut ref_bytes = Vec::new();
        push_packed_loc_id_group(&mut ref_bytes, &[(PackedLoc::new(10), 1)]);

        let mut query_hits = [DoubleArray::from_bytes(query_bytes, 22)];
        let mut ref_hits = [DoubleArray::from_bytes(ref_bytes, 13)];
        let stats = FrequentSeeds::build_packed_loc_id(
            0,
            &SeedPartitionRange::with_bounds(0, 1),
            &mut query_hits,
            &mut ref_hits,
            &mut seqs,
            0.0,
        );

        assert_eq!(stats.query_max_n, 2);
        assert_eq!(stats.ref_max_n, 1);
        assert_eq!(stats.masked_positions, 0);
        assert_eq!(seqs.data_at(p0.as_u64())[0] & SEED_MASK, 0);
        assert!(query_hits[0].begin_with_elem_size(9).good());
        assert!(ref_hits[0].begin_with_elem_size(9).good());
    }

    fn run_parallel_case(
        threads: usize,
    ) -> (
        FrequentSeedsBuildStats,
        Vec<DoubleArray>,
        Vec<DoubleArray>,
        SequenceSet,
    ) {
        let mut seqs = SequenceSet::new();
        seqs.push(&[1, 2, 3, 4]);
        let pos = PackedLoc::new(seqs.position(0, 1) as u64);

        let mut q0 = Vec::new();
        push_packed_loc_group(&mut q0, &[pos]);
        let mut q1 = Vec::new();
        push_packed_loc_group(&mut q1, &[pos, pos, pos]);
        let mut r0 = Vec::new();
        push_packed_loc_group(
            &mut r0,
            &[
                PackedLoc::new(10),
                PackedLoc::new(11),
                PackedLoc::new(12),
                PackedLoc::new(13),
            ],
        );
        let mut r1 = Vec::new();
        push_packed_loc_group(&mut r1, &[PackedLoc::new(20)]);
        let mut query = vec![
            DoubleArray::from_bytes(q0, 9),
            DoubleArray::from_bytes(q1, 19),
        ];
        let mut reference = vec![
            DoubleArray::from_bytes(r0, 24),
            DoubleArray::from_bytes(r1, 9),
        ];
        let stats = FrequentSeeds::build_with_config::<PackedLoc>(
            0,
            &SeedPartitionRange::with_bounds(0, 2),
            &mut query,
            &mut reference,
            &mut seqs,
            FrequentSeedsConfig {
                freq_sd: 0.0,
                threads,
            },
        )
        .unwrap();
        (stats, query, reference, seqs)
    }

    #[test]
    fn test_explicit_parallel_build_matches_single_thread_statistics_and_masks() {
        let single = run_parallel_case(1);
        let parallel = run_parallel_case(2);
        assert_eq!(single.0, parallel.0);
        assert_eq!(single.1, parallel.1);
        assert_eq!(single.2, parallel.2);
        assert_eq!(single.3.get(0), parallel.3.get(0));

        let stats = single.0;
        assert!((stats.ref_mean - 2.5).abs() < 1e-12);
        assert!((stats.ref_sd - 3.0f64.sqrt()).abs() < 1e-12);
        assert!((stats.query_mean - 2.0).abs() < 1e-12);
        assert!((stats.query_sd - (4.0f64 / 3.0).sqrt()).abs() < 1e-12);
        assert_eq!((stats.ref_max_n, stats.query_max_n), (2, 2));
        assert_eq!(stats.masked_positions, 5);
        assert_ne!(single.3.get(0)[1] & SEED_MASK, 0);
        assert!(!single.1[0].begin_with_elem_size(5).good());
        assert!(!single.1[1].begin_with_elem_size(5).good());
    }

    #[test]
    fn test_explicit_parallel_config_validation() {
        let mut query = Vec::new();
        let mut reference = Vec::new();
        let mut seqs = SequenceSet::new();
        assert_eq!(
            FrequentSeeds::build_with_config::<PackedLoc>(
                0,
                &SeedPartitionRange::with_bounds(0, 1),
                &mut query,
                &mut reference,
                &mut seqs,
                FrequentSeedsConfig {
                    freq_sd: 0.0,
                    threads: 0,
                },
            )
            .unwrap_err(),
            "Frequent seed thread count must be positive."
        );
    }
}
