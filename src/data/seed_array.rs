use crate::basic::packed_loc::PackedLoc;
use crate::basic::seed::{
    seed_partition, seed_partition_offset, seedp_mask, PackedSeed, SeedOffset,
};
use crate::data::block::Block;
use crate::data::enum_seeds::{enum_seeds_block, EnumSeedsCallback, EnumSeedsContext};
use crate::data::flags::{EnumCfg, NoFilter, PackedLocId, SeedEncoding};
use crate::data::seed_histogram::{partition_size, SeedPartitionRange, ShapeHistogram};
use crate::data::seed_set::{HashedSeedSet, SeedSet};
use crate::search::seed_complexity::SeedStats;
use crate::util::algo::{HashJoinRecord, HashJoinValue};

pub trait SeedLocation: Copy + Default + Eq {
    fn make_seed_loc(pos: PackedLoc, block_id: u32) -> Self;
}

impl SeedLocation for PackedLoc {
    fn make_seed_loc(pos: PackedLoc, _block_id: u32) -> Self {
        pos
    }
}

impl SeedLocation for PackedLocId {
    fn make_seed_loc(pos: PackedLoc, block_id: u32) -> Self {
        PackedLocId::with_block_id(pos, block_id)
    }
}

#[repr(C, packed)]
#[derive(Clone, Copy, Eq)]
pub struct SeedArrayEntry<SeedLoc: SeedLocation> {
    pub key: SeedOffset,
    pub value: SeedLoc,
}

impl<SeedLoc: SeedLocation> SeedArrayEntry<SeedLoc> {
    pub fn new(key: SeedOffset, value: SeedLoc) -> Self {
        Self { key, value }
    }

    pub fn with_block_id(key: SeedOffset, pos: PackedLoc, block_id: u32) -> Self {
        Self {
            key,
            value: SeedLoc::make_seed_loc(pos, block_id),
        }
    }

    pub fn get_key(&self) -> u32 {
        self.key
    }
}

impl HashJoinValue for PackedLoc {
    fn write_ne_bytes(self, out: &mut Vec<u8>) {
        out.push(self.high);
        out.extend_from_slice(&self.low.to_ne_bytes());
    }
}

impl HashJoinValue for PackedLocId {
    fn write_ne_bytes(self, out: &mut Vec<u8>) {
        let pos = self.pos;
        pos.write_ne_bytes(out);
        let block_id = self.block_id;
        out.extend_from_slice(&block_id.to_ne_bytes());
    }
}

impl<SeedLoc> HashJoinRecord for SeedArrayEntry<SeedLoc>
where
    SeedLoc: SeedLocation + HashJoinValue,
{
    type Key = SeedOffset;
    type Value = SeedLoc;

    fn key(&self) -> Self::Key {
        self.key
    }

    fn value(&self) -> Self::Value {
        self.value
    }
}

impl<SeedLoc: SeedLocation> Default for SeedArrayEntry<SeedLoc> {
    fn default() -> Self {
        Self {
            key: 0,
            value: SeedLoc::default(),
        }
    }
}

impl<SeedLoc: SeedLocation> PartialEq for SeedArrayEntry<SeedLoc> {
    fn eq(&self, rhs: &Self) -> bool {
        let key = self.key;
        let value = self.value;
        let rhs_key = rhs.key;
        let rhs_value = rhs.value;
        key == rhs_key && value == rhs_value
    }
}

impl<SeedLoc: SeedLocation + std::fmt::Debug> std::fmt::Debug for SeedArrayEntry<SeedLoc> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let key = self.key;
        let value = self.value;
        f.debug_struct("SeedArrayEntry")
            .field("key", &key)
            .field("value", &value)
            .finish()
    }
}

impl<SeedLoc: SeedLocation> PartialOrd for SeedArrayEntry<SeedLoc> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<SeedLoc: SeedLocation> Ord for SeedArrayEntry<SeedLoc> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let key = self.key;
        let other_key = other.key;
        key.cmp(&other_key)
    }
}

#[derive(Debug, Clone)]
pub struct SeedArray<SeedLoc: SeedLocation> {
    key_bits: i32,
    data: Option<Vec<SeedArrayEntry<SeedLoc>>>,
    begin: Vec<i64>,
    entries: Vec<Vec<SeedArrayEntry<SeedLoc>>>,
    stats: SeedStats,
}

impl<SeedLoc: SeedLocation> SeedArray<SeedLoc> {
    pub fn alloc_buffer(
        hst: &crate::data::seed_histogram::SeedHistogram,
        index_chunks: i32,
    ) -> usize {
        std::mem::size_of::<SeedArrayEntry<SeedLoc>>() * hst.max_chunk_size(index_chunks)
    }

    pub fn from_histogram<Filter>(
        seqs: &mut Block,
        hst: &ShapeHistogram,
        range: &SeedPartitionRange,
        seedp_bits: i32,
        filter: &Filter,
        enum_cfg: &EnumCfg<'_>,
        ctx: &EnumSeedsContext<'_>,
    ) -> Result<Self, String>
    where
        Filter: crate::data::enum_seeds::EnumSeedsFilter,
    {
        if enum_cfg.shape_end - enum_cfg.shape_begin > 1 {
            return Err("SeedArray construction for >1 shape.".to_string());
        }
        let key_bits = seed_bits(enum_cfg.code, seedp_bits, ctx)?;
        let mut begin = Vec::with_capacity(range.size() as usize + 1);
        begin.push(0);
        for i in range.begin()..range.end() {
            let next = *begin.last().unwrap() + partition_size(hst, i as usize) as i64;
            begin.push(next);
        }

        let mut data = vec![SeedArrayEntry::default(); *begin.last().unwrap() as usize];
        let mut iterators = build_iterators(&begin, hst, range);
        let partition = enum_cfg
            .partition
            .ok_or_else(|| "EnumCfg::partition is required.".to_string())?;
        let mut callbacks: Vec<BuildCallback<SeedLoc>> = (0..partition.len() - 1)
            .map(|_| BuildCallback::new(*range, &iterators, seedp_bits))
            .collect();
        let stats = enum_seeds_block(seqs, &mut callbacks, filter, enum_cfg, ctx)?;

        for cb in callbacks {
            for (shape, row) in cb.entries.into_iter().enumerate() {
                for (part, entries) in row.into_iter().enumerate() {
                    let mut dst = iterators[shape][part];
                    for entry in entries {
                        data[dst] = entry;
                        dst += 1;
                    }
                    iterators[shape][part] = dst;
                }
            }
        }

        Ok(Self {
            key_bits,
            data: Some(data),
            begin,
            entries: Vec::new(),
            stats,
        })
    }

    pub fn one_pass<Filter>(
        seqs: &mut Block,
        range: &SeedPartitionRange,
        seedp_bits: i32,
        filter: &Filter,
        enum_cfg: &EnumCfg<'_>,
        ctx: &EnumSeedsContext<'_>,
        n_parts: u32,
    ) -> Result<Self, String>
    where
        Filter: crate::data::enum_seeds::EnumSeedsFilter,
    {
        if enum_cfg.shape_end - enum_cfg.shape_begin > 1 {
            return Err("SeedArray construction for >1 shape.".to_string());
        }
        let key_bits = seed_bits(enum_cfg.code, seedp_bits, ctx)?;
        let seq_partition = seqs.seqs().partition(n_parts.max(1), false, false);
        let cfg = EnumCfg {
            partition: Some(&seq_partition),
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
        let mut callbacks: Vec<OnePassBuildCallback<SeedLoc>> = (0..seq_partition.len() - 1)
            .map(|_| OnePassBuildCallback::new(*range, seedp_bits))
            .collect();
        let stats = enum_seeds_block(seqs, &mut callbacks, filter, &cfg, ctx)?;

        let mut counts = vec![0usize; range.size() as usize];
        for cb in &callbacks {
            for i in 0..range.size() as usize {
                counts[i] += cb.it.out[i].len();
            }
        }

        let mut entries = Vec::with_capacity(range.size() as usize);
        for i in 0..range.size() as usize {
            let mut dst = Vec::with_capacity(counts[i]);
            for cb in &mut callbacks {
                dst.append(&mut cb.it.out[i]);
            }
            entries.push(dst);
        }

        Ok(Self {
            key_bits,
            data: None,
            begin: Vec::new(),
            entries,
            stats,
        })
    }

    pub fn begin(&self, i: usize) -> &[SeedArrayEntry<SeedLoc>] {
        if let Some(data) = &self.data {
            &data[self.begin[i] as usize..self.begin[i + 1] as usize]
        } else {
            &self.entries[i]
        }
    }

    pub fn begin_mut(&mut self, i: usize) -> &mut [SeedArrayEntry<SeedLoc>] {
        if let Some(data) = &mut self.data {
            &mut data[self.begin[i] as usize..self.begin[i + 1] as usize]
        } else {
            &mut self.entries[i]
        }
    }

    pub fn size_partition(&self, i: usize) -> usize {
        if self.data.is_some() {
            (self.begin[i + 1] - self.begin[i]) as usize
        } else {
            self.entries[i].len()
        }
    }

    pub fn size(&self) -> usize {
        if self.data.is_some() {
            *self.begin.last().unwrap_or(&0) as usize
        } else {
            self.entries.iter().map(Vec::len).sum()
        }
    }

    pub fn stats(&self) -> &SeedStats {
        &self.stats
    }

    pub fn key_bits(&self) -> i32 {
        self.key_bits
    }
}

pub fn seed_bits(
    code: SeedEncoding,
    seedp_bits: i32,
    ctx: &EnumSeedsContext<'_>,
) -> Result<i32, String> {
    match code {
        SeedEncoding::Hashed => Ok((std::mem::size_of::<SeedOffset>() * 8) as i32),
        SeedEncoding::SpacedFactor => Ok((ctx.shapes.get(0).weight as f64
            * ctx.reduction.bit_size_exact())
        .ceil() as i32
            - seedp_bits),
        SeedEncoding::Contiguous => {
            Ok(ctx.shapes.get(0).length * ctx.reduction.bit_size() - seedp_bits)
        }
    }
}

struct BufferedWriter<SeedLoc: SeedLocation> {
    seedp_mask: PackedSeed,
    seedp_bits: i32,
    range: SeedPartitionRange,
    out: Vec<Vec<SeedArrayEntry<SeedLoc>>>,
}

impl<SeedLoc: SeedLocation> BufferedWriter<SeedLoc> {
    fn new(range: SeedPartitionRange, seedp_bits: i32) -> Self {
        Self {
            seedp_mask: seedp_mask(seedp_bits),
            seedp_bits,
            range,
            out: vec![Vec::new(); range.size() as usize],
        }
    }

    fn push(&mut self, key: PackedSeed, value: usize, block_id: u32) {
        let p = seed_partition(key, self.seedp_mask);
        if self.range.contains(p) {
            let d = (p - self.range.begin()) as usize;
            self.out[d].push(SeedArrayEntry::with_block_id(
                seed_partition_offset(key, self.seedp_bits as u64),
                PackedLoc::new(value as u64),
                block_id,
            ));
        }
    }
}

fn build_iterators(
    begin: &[i64],
    hst: &ShapeHistogram,
    range: &SeedPartitionRange,
) -> Vec<Vec<usize>> {
    let mut iterators = vec![Vec::with_capacity(range.size() as usize); hst.len()];
    for i in 0..range.size() as usize {
        iterators[0].push(begin[i] as usize);
    }
    for i in 1..hst.len() {
        for j in range.begin()..range.end() {
            let d = (j - range.begin()) as usize;
            let prev = iterators[i - 1][d];
            iterators[i].push(prev + hst[i - 1][j as usize] as usize);
        }
    }
    iterators
}

struct BuildCallback<SeedLoc: SeedLocation> {
    range: SeedPartitionRange,
    it: BufferedWriter<SeedLoc>,
    entries: Vec<Vec<Vec<SeedArrayEntry<SeedLoc>>>>,
}

impl<SeedLoc: SeedLocation> BuildCallback<SeedLoc> {
    fn new(range: SeedPartitionRange, iterators: &[Vec<usize>], seedp_bits: i32) -> Self {
        Self {
            range,
            it: BufferedWriter::new(range, seedp_bits),
            entries: vec![vec![Vec::new(); range.size() as usize]; iterators.len()],
        }
    }
}

impl<SeedLoc: SeedLocation> EnumSeedsCallback for BuildCallback<SeedLoc> {
    fn call(&mut self, seed: u64, pos: usize, block_id: usize, shape: i32) -> bool {
        self.it.push(seed, pos, block_id as u32);
        let p = seed_partition(seed, self.it.seedp_mask);
        if self.range.contains(p) {
            let d = (p - self.range.begin()) as usize;
            self.entries[shape as usize][d].push(SeedArrayEntry::with_block_id(
                seed_partition_offset(seed, self.it.seedp_bits as u64),
                PackedLoc::new(pos as u64),
                block_id as u32,
            ));
        }
        true
    }
}

struct OnePassBufferedWriter<SeedLoc: SeedLocation> {
    seedp_bits: i32,
    out: Vec<Vec<SeedArrayEntry<SeedLoc>>>,
    range: SeedPartitionRange,
}

impl<SeedLoc: SeedLocation> OnePassBufferedWriter<SeedLoc> {
    fn new(range: SeedPartitionRange, seedp_bits: i32) -> Self {
        Self {
            seedp_bits,
            out: vec![Vec::new(); range.size() as usize],
            range,
        }
    }

    fn push(&mut self, key: PackedSeed, value: usize) {
        let p = seed_partition(key, seedp_mask(self.seedp_bits));
        if self.range.contains(p) {
            let d = (p - self.range.begin()) as usize;
            self.out[d].push(SeedArrayEntry::new(
                seed_partition_offset(key, self.seedp_bits as u64),
                SeedLoc::make_seed_loc(PackedLoc::new(value as u64), 0),
            ));
        }
    }
}

struct OnePassBuildCallback<SeedLoc: SeedLocation> {
    it: OnePassBufferedWriter<SeedLoc>,
}

impl<SeedLoc: SeedLocation> OnePassBuildCallback<SeedLoc> {
    fn new(range: SeedPartitionRange, seedp_bits: i32) -> Self {
        Self {
            it: OnePassBufferedWriter::new(range, seedp_bits),
        }
    }
}

impl<SeedLoc: SeedLocation> EnumSeedsCallback for OnePassBuildCallback<SeedLoc> {
    fn call(&mut self, seed: u64, pos: usize, _block_id: usize, _shape: i32) -> bool {
        self.it.push(seed, pos);
        true
    }
}

impl crate::data::enum_seeds::EnumSeedsFilter for SeedSet {
    fn contains(&self, key: u64, shape_id: i32) -> bool {
        self.contains(key, shape_id as u64)
    }
}

impl crate::data::enum_seeds::EnumSeedsFilter for HashedSeedSet {
    fn contains(&self, key: u64, shape_id: i32) -> bool {
        self.contains(key, shape_id as u64)
    }
}

impl crate::data::enum_seeds::EnumSeedsFilter for &NoFilter {
    fn contains(&self, key: u64, shape_id: i32) -> bool {
        (*self).contains(key, shape_id as u64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::reduction::Reduction;
    use crate::basic::shape_config::ShapeConfig;
    use crate::basic::value::SequenceType;
    use crate::data::flags::NO_FILTER;
    use crate::data::seed_histogram::SeedPartitionRange;
    use crate::masking::MaskingAlgo;

    fn enum_cfg<'a>(partition: Option<&'a Vec<u32>>, code: SeedEncoding) -> EnumCfg<'a> {
        EnumCfg {
            partition,
            shape_begin: 0,
            shape_end: 1,
            code,
            skip: None,
            filter_masked_seeds: false,
            mask_seeds: false,
            seed_cut: 0.0,
            soft_masking: MaskingAlgo::None,
            minimizer_window: 0,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: 0,
        }
    }

    fn block() -> Block {
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2, 3],
                Some("q0"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        block
            .push_back(
                &[4, 5, 6],
                Some("q1"),
                None,
                1,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        block
    }

    #[test]
    fn test_seed_array_entry_layout_make_loc_and_order() {
        assert_eq!(
            std::mem::size_of::<SeedArrayEntry<PackedLoc>>(),
            std::mem::size_of::<SeedOffset>() + std::mem::size_of::<PackedLoc>()
        );
        assert_eq!(
            std::mem::size_of::<SeedArrayEntry<PackedLocId>>(),
            std::mem::size_of::<SeedOffset>() + std::mem::size_of::<PackedLocId>()
        );
        let a = SeedArrayEntry::<PackedLoc>::with_block_id(2, PackedLoc::new(10), 99);
        assert_eq!(a.get_key(), 2);
        assert_eq!(a.value.as_u64(), 10);
        let b = SeedArrayEntry::<PackedLocId>::with_block_id(1, PackedLoc::new(7), 42);
        let b_value = b.value;
        let b_block_id = b_value.block_id;
        assert_eq!(b_value.as_u64(), 7);
        assert_eq!(b_block_id, 42);
        assert!(b < SeedArrayEntry::<PackedLocId>::with_block_id(2, PackedLoc::new(0), 0));
    }

    #[test]
    fn test_seed_bits_matches_cpp_cases() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["101".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        assert_eq!(seed_bits(SeedEncoding::Hashed, 4, &ctx).unwrap(), 32);
        assert_eq!(seed_bits(SeedEncoding::SpacedFactor, 1, &ctx).unwrap(), 6);
        assert_eq!(seed_bits(SeedEncoding::Contiguous, 1, &ctx).unwrap(), 11);
    }

    #[test]
    fn test_seed_array_one_pass_builds_partitioned_entries() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["11".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut block = block();
        let range = SeedPartitionRange::with_bounds(0, 4);
        let cfg = enum_cfg(None, SeedEncoding::SpacedFactor);
        let sa = SeedArray::<PackedLoc>::one_pass(&mut block, &range, 2, &NO_FILTER, &cfg, &ctx, 2)
            .unwrap();
        assert_eq!(sa.key_bits(), 5);
        assert_eq!(sa.size(), 5);
        let total: usize = (0..range.size() as usize)
            .map(|i| sa.size_partition(i))
            .sum();
        assert_eq!(total, sa.size());
        assert_eq!(sa.stats(), &SeedStats::default());
    }

    #[test]
    fn test_seed_array_from_histogram_flat_storage() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["11".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut block = block();
        let range = SeedPartitionRange::with_bounds(0, 4);
        let partition = vec![0, 1, 2];
        let cfg = enum_cfg(Some(&partition), SeedEncoding::SpacedFactor);
        let hst = vec![vec![2, 1, 1, 1]];
        let sa = SeedArray::<PackedLocId>::from_histogram(
            &mut block, &hst, &range, 2, &NO_FILTER, &cfg, &ctx,
        )
        .unwrap();
        assert_eq!(sa.size(), 5);
        assert_eq!(sa.size_partition(0), 2);
        assert_eq!(sa.begin(0).len(), 2);
        assert!(sa.begin(0).iter().any(|e| e.value.block_id == 0));
        assert!(sa.begin(0).iter().any(|e| e.value.block_id == 1));
    }

    #[test]
    fn test_seed_array_rejects_multiple_shapes() {
        let reduction = Reduction::default_reduction();
        let shapes =
            ShapeConfig::from_codes(&["11".to_string(), "101".to_string()], 0, &reduction).unwrap();
        let ctx = EnumSeedsContext {
            shapes: &shapes,
            reduction: &reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        let mut block = block();
        let range = SeedPartitionRange::with_bounds(0, 4);
        let partition = vec![0, 2];
        let mut cfg = enum_cfg(Some(&partition), SeedEncoding::SpacedFactor);
        cfg.shape_end = 2;
        let err = SeedArray::<PackedLoc>::from_histogram(
            &mut block,
            &vec![vec![0; 4], vec![0; 4]],
            &range,
            2,
            &NO_FILTER,
            &cfg,
            &ctx,
        )
        .unwrap_err();
        assert_eq!(err, "SeedArray construction for >1 shape.");
    }
}
