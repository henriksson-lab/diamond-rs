use crate::basic::reduction::Reduction;
use crate::basic::seed::{seed_partition, seedp_count, seedp_mask, PackedSeed, SeedPartition};
use crate::basic::shape_config::ShapeConfig;
use crate::basic::value::Loc;
use crate::data::block::Block;
use crate::data::enum_seeds::{
    enum_seeds_block, EnumSeedsCallback, EnumSeedsContext, EnumSeedsFilter,
};
use crate::data::flags::EnumCfg;
use crate::util::algo::Partition;
use std::sync::Mutex;

pub type ShapeHistogram = Vec<Vec<u32>>;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeedPartitionRange {
    begin: SeedPartition,
    end: SeedPartition,
}

impl SeedPartitionRange {
    pub const fn new() -> Self {
        Self { begin: 0, end: 0 }
    }

    pub const fn with_bounds(begin: SeedPartition, end: SeedPartition) -> Self {
        Self { begin, end }
    }

    /// Matches C++ `SeedPartitionRange::contains(i)`.
    pub fn contains(&self, i: SeedPartition) -> bool {
        i >= self.begin && i < self.end
    }

    /// Matches C++ `SeedPartitionRange::begin()`.
    pub fn begin(&self) -> SeedPartition {
        self.begin
    }

    /// Matches C++ `SeedPartitionRange::end()`.
    pub fn end(&self) -> SeedPartition {
        self.end
    }

    /// Matches C++ `SeedPartitionRange::lower(i)`.
    pub fn lower(&self, i: SeedPartition) -> bool {
        i < self.begin
    }

    /// Matches C++ `SeedPartitionRange::lower_or_equal(i)`.
    pub fn lower_or_equal(&self, i: SeedPartition) -> bool {
        i < self.end
    }

    /// Matches C++ `SeedPartitionRange::size()`.
    pub fn size(&self) -> SeedPartition {
        self.end - self.begin
    }
}

pub static CURRENT_RANGE: Mutex<SeedPartitionRange> = Mutex::new(SeedPartitionRange::new());

/// Matches C++ `partition_size(hst, p)`.
pub fn partition_size(hst: &ShapeHistogram, p: usize) -> usize {
    let mut s = 0usize;
    for row in hst {
        s += row[p] as usize;
    }
    s
}

/// Matches C++ `hst_size(hst, range)`.
pub fn hst_size(hst: &ShapeHistogram, range: &SeedPartitionRange) -> usize {
    let mut s = 0usize;
    for i in range.begin()..range.end() {
        s += partition_size(hst, i as usize);
    }
    s
}

pub struct SeedHistogram {
    p: Vec<u32>,
    data: Vec<ShapeHistogram>,
}

impl SeedHistogram {
    /// Matches C++ `SeedHistogram::SeedHistogram()`.
    pub fn new() -> Self {
        Self {
            p: Vec::new(),
            data: Vec::new(),
        }
    }

    /// Matches C++ `SeedHistogram::SeedHistogram(seqs, serial, filter, enum_cfg, seedp_bits, threads, shapes, reduction, min_query_len, query_contexts)`.
    pub fn from_block<Filter>(
        seqs: &mut Block,
        serial: bool,
        filter: &Filter,
        enum_cfg: &EnumCfg<'_>,
        seedp_bits: i32,
        threads: u32,
        shapes: &ShapeConfig,
        reduction: &Reduction,
        min_query_len: Loc,
        query_contexts: usize,
    ) -> Result<Self, String>
    where
        Filter: EnumSeedsFilter,
    {
        let p = seqs.seqs().partition(threads, false, false);
        let seedp = seedp_count(seedp_bits) as usize;
        let n_shapes = shapes.count() as usize;
        let mut callbacks: Vec<SeedHistogramCallback> = (0..p.len().saturating_sub(1))
            .map(|_| SeedHistogramCallback::new(seedp_bits, n_shapes, seedp))
            .collect();

        let mut cfg = EnumCfg {
            partition: Some(&p),
            shape_begin: enum_cfg.shape_begin,
            shape_end: enum_cfg.shape_end,
            code: enum_cfg.code,
            skip: enum_cfg.skip,
            filter_masked_seeds: enum_cfg.filter_masked_seeds,
            mask_seeds: enum_cfg.mask_seeds,
            seed_cut: enum_cfg.seed_cut,
            soft_masking: enum_cfg.soft_masking,
            minimizer_window: enum_cfg.minimizer_window,
            filter_low_complexity_seeds: enum_cfg.filter_low_complexity_seeds,
            mask_low_complexity_seeds: enum_cfg.mask_low_complexity_seeds,
            sketch_size: enum_cfg.sketch_size,
        };
        let ctx = EnumSeedsContext {
            shapes,
            reduction,
            min_query_len,
            query_contexts,
        };

        if serial {
            for s in 0..shapes.count() {
                cfg.shape_begin = s;
                cfg.shape_end = s + 1;
                enum_seeds_block(seqs, &mut callbacks, filter, &cfg, &ctx)?;
            }
        } else {
            cfg.shape_begin = 0;
            cfg.shape_end = shapes.count();
            enum_seeds_block(seqs, &mut callbacks, filter, &cfg, &ctx)?;
        }

        let mut data = vec![vec![vec![0u32; seedp]; p.len().saturating_sub(1)]; n_shapes];
        for (seqp, cb) in callbacks.into_iter().enumerate() {
            for (shape, row) in cb.rows.into_iter().enumerate() {
                data[shape][seqp] = row;
            }
        }

        Ok(Self { p, data })
    }

    /// Matches C++ `SeedHistogram::get(sid)`.
    pub fn get(&self, sid: usize) -> &ShapeHistogram {
        &self.data[sid]
    }

    /// Matches C++ `SeedHistogram::max_chunk_size(index_chunks)`.
    pub fn max_chunk_size(&self, index_chunks: i32) -> usize {
        let mut max = 0usize;
        let p = Partition::new(self.seedp() as usize, index_chunks as usize);
        for shape in 0..self.data.len() {
            for chunk in 0..p.parts {
                max = max.max(hst_size(
                    &self.data[shape],
                    &SeedPartitionRange::with_bounds(
                        p.begin(chunk) as SeedPartition,
                        p.end(chunk) as SeedPartition,
                    ),
                ));
            }
        }
        max
    }

    /// Matches C++ `SeedHistogram::partition()`.
    pub fn partition(&self) -> &Vec<u32> {
        &self.p
    }

    /// Matches C++ `SeedHistogram::seedp()`.
    pub fn seedp(&self) -> i32 {
        self.data[0][0].len() as i32
    }
}

struct SeedHistogramCallback {
    seedp_mask: PackedSeed,
    rows: Vec<Vec<u32>>,
}

impl SeedHistogramCallback {
    fn new(seedp_bits: i32, n_shapes: usize, seedp: usize) -> Self {
        Self {
            seedp_mask: seedp_mask(seedp_bits),
            rows: vec![vec![0; seedp]; n_shapes],
        }
    }
}

impl EnumSeedsCallback for SeedHistogramCallback {
    fn call(&mut self, seed: u64, _pos: usize, _block_id: usize, shape: i32) -> bool {
        self.rows[shape as usize][seed_partition(seed, self.seedp_mask) as usize] += 1;
        true
    }
}

impl Default for SeedHistogram {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::reduction::Reduction;
    use crate::basic::shape_config::ShapeConfig;
    use crate::basic::value::SequenceType;
    use crate::data::block::Block;
    use crate::data::flags::{EnumCfg, SeedEncoding, NO_FILTER};
    use crate::masking::MaskingAlgo;

    #[test]
    fn test_seed_partition_range() {
        let r = SeedPartitionRange::with_bounds(3, 8);
        assert!(r.contains(3));
        assert!(r.contains(7));
        assert!(!r.contains(8));
        assert!(r.lower(2));
        assert!(!r.lower(3));
        assert!(r.lower_or_equal(7));
        assert!(!r.lower_or_equal(8));
        assert_eq!(r.begin(), 3);
        assert_eq!(r.end(), 8);
        assert_eq!(r.size(), 5);
    }

    #[test]
    fn test_partition_and_histogram_size() {
        let hst = vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8]];
        assert_eq!(partition_size(&hst, 0), 6);
        assert_eq!(partition_size(&hst, 2), 10);
        assert_eq!(hst_size(&hst, &SeedPartitionRange::with_bounds(1, 3)), 18);
    }

    #[test]
    fn test_seed_histogram_accessors() {
        let h = SeedHistogram {
            p: vec![0, 2, 4],
            data: vec![
                vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8]],
                vec![vec![2, 0, 1, 0], vec![1, 1, 1, 1]],
            ],
        };
        assert_eq!(h.partition(), &vec![0, 2, 4]);
        assert_eq!(h.seedp(), 4);
        assert_eq!(h.get(1)[0], vec![2, 0, 1, 0]);
        assert_eq!(h.max_chunk_size(2), 22);
    }

    #[test]
    fn test_current_range() {
        let mut r = CURRENT_RANGE.lock().unwrap();
        *r = SeedPartitionRange::with_bounds(1, 2);
        assert_eq!(r.size(), 1);
        *r = SeedPartitionRange::new();
    }

    #[test]
    fn test_seed_histogram_from_block_spaced_factor() {
        let reduction = Reduction::default_reduction();
        let shapes =
            ShapeConfig::from_codes(&["11".to_string(), "101".to_string()], 0, &reduction).unwrap();
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2, 3],
                Some("s0"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        block
            .push_back(
                &[1, 2, 3],
                Some("s1"),
                None,
                1,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let cfg = EnumCfg {
            partition: None,
            shape_begin: 0,
            shape_end: shapes.count(),
            code: SeedEncoding::SpacedFactor,
            skip: None,
            filter_masked_seeds: false,
            mask_seeds: false,
            seed_cut: 0.0,
            soft_masking: MaskingAlgo::None,
            minimizer_window: 0,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: 0,
        };

        let h = SeedHistogram::from_block(
            &mut block, false, &NO_FILTER, &cfg, 2, 1, &shapes, &reduction, 0, 1,
        )
        .unwrap();

        assert_eq!(h.seedp(), 4);
        assert_eq!(h.partition(), &vec![0, 2]);
        assert_eq!(h.get(0), &vec![vec![2, 1, 2, 0]]);
        assert_eq!(h.get(1), &vec![vec![2, 0, 1, 0]]);
        assert_eq!(h.max_chunk_size(2), 3);
    }

    #[test]
    fn test_seed_histogram_from_block_serial_matches_parallel() {
        let reduction = Reduction::default_reduction();
        let shapes =
            ShapeConfig::from_codes(&["11".to_string(), "101".to_string()], 0, &reduction).unwrap();
        let mut parallel_block = Block::new();
        let mut serial_block = Block::new();
        for block in [&mut parallel_block, &mut serial_block] {
            block
                .push_back(
                    &[0, 1, 2, 3],
                    Some("s0"),
                    None,
                    0,
                    SequenceType::AminoAcid,
                    0,
                    false,
                )
                .unwrap();
            block
                .push_back(
                    &[1, 2, 3],
                    Some("s1"),
                    None,
                    1,
                    SequenceType::AminoAcid,
                    0,
                    false,
                )
                .unwrap();
        }
        let cfg = EnumCfg {
            partition: None,
            shape_begin: 0,
            shape_end: shapes.count(),
            code: SeedEncoding::SpacedFactor,
            skip: None,
            filter_masked_seeds: false,
            mask_seeds: false,
            seed_cut: 0.0,
            soft_masking: MaskingAlgo::None,
            minimizer_window: 0,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: 0,
        };

        let parallel = SeedHistogram::from_block(
            &mut parallel_block,
            false,
            &NO_FILTER,
            &cfg,
            2,
            1,
            &shapes,
            &reduction,
            0,
            1,
        )
        .unwrap();
        let serial = SeedHistogram::from_block(
            &mut serial_block,
            true,
            &NO_FILTER,
            &cfg,
            2,
            1,
            &shapes,
            &reduction,
            0,
            1,
        )
        .unwrap();

        assert_eq!(serial.partition(), parallel.partition());
        assert_eq!(serial.get(0), parallel.get(0));
        assert_eq!(serial.get(1), parallel.get(1));
    }
}
