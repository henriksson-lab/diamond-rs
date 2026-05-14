use crate::basic::reduction::Reduction;
use crate::basic::shape_config::ShapeConfig;
use crate::data::block::Block;
use crate::data::enum_seeds::{enum_seeds_block, EnumSeedsCallback};
use crate::data::flags::{EnumCfg, SeedEncoding, NO_FILTER};
use crate::masking::MaskingAlgo;
use crate::util::data_structures::{HashSet as DiamondHashSet, Identity, Modulo2};
use crate::util::math::next_pow2;

pub type Table = DiamondHashSet<Modulo2, Identity>;

pub const SEED_INDEX_MAGIC_NUMBER: u64 = 0x2d6ba306ecbf6aba;
pub const SEED_INDEX_VERSION: u32 = 0;
pub const SEED_INDEX_HEADER_SIZE: usize = 16;
pub const HASH_TABLE_FACTOR: f64 = 1.25;

#[derive(Debug, Clone)]
pub struct SeedSet {
    data: Vec<bool>,
    coverage: f64,
}

impl SeedSet {
    /// Matches C++ `SeedSet::SeedSet`.
    pub fn from_block(
        seqs: &mut Block,
        max_coverage: f64,
        skip: Option<&Vec<bool>>,
        seed_cut: f64,
        soft_masking: MaskingAlgo,
        shapes: &ShapeConfig,
        reduction: &Reduction,
    ) -> Result<Self, String> {
        if !shapes.get(0).contiguous() {
            return Err("Contiguous seed required.".to_string());
        }
        let data_size = (1u64 << reduction.bit_size()).pow(shapes.get(0).length as u32) as usize;
        let max_coverage =
            (max_coverage * (reduction.size() as f64).powi(shapes.get(0).length)) as usize;
        let mut cb = vec![SeedSetCallback::new(vec![false; data_size], max_coverage)];
        let partition = seqs.seqs().partition(1, false, false);
        let cfg = EnumCfg {
            partition: Some(&partition),
            shape_begin: 0,
            shape_end: 1,
            code: SeedEncoding::Contiguous,
            skip,
            filter_masked_seeds: true,
            mask_seeds: false,
            seed_cut,
            soft_masking,
            minimizer_window: 0,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: 0,
        };
        let ctx = crate::data::enum_seeds::EnumSeedsContext {
            shapes,
            reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        enum_seeds_block(seqs, &mut cb, &NO_FILTER, &cfg, &ctx)?;
        let cb = cb.pop().unwrap();
        let coverage = cb.coverage as f64 / (reduction.size() as f64).powi(shapes.get(0).length);
        Ok(Self {
            data: cb.data,
            coverage,
        })
    }

    /// Matches C++ `SeedSet::contains`.
    pub fn contains(&self, key: u64, _shape: u64) -> bool {
        self.data[key as usize]
    }

    /// Matches C++ `SeedSet::coverage`.
    pub fn coverage(&self) -> f64 {
        self.coverage
    }
}

#[derive(Debug, Clone)]
pub struct HashedSeedSet {
    data: Vec<Table>,
}

impl HashedSeedSet {
    /// Matches C++ `HashedSeedSet::HashedSeedSet(Block&, ...)`.
    pub fn from_block(
        seqs: &mut Block,
        skip: Option<&Vec<bool>>,
        seed_cut: f64,
        soft_masking: MaskingAlgo,
        shapes: &ShapeConfig,
        reduction: &Reduction,
    ) -> Result<Self, String> {
        let shape_count = shapes.count() as usize;
        let initial_capacity = next_pow2(seqs.seqs().letters() as f64 * HASH_TABLE_FACTOR);
        let mut cb = vec![HashedSeedSetCallback::new(
            (0..shape_count)
                .map(|_| Table::new(initial_capacity))
                .collect(),
        )];
        let partition = seqs.seqs().partition(1, false, false);
        let cfg = EnumCfg {
            partition: Some(&partition),
            shape_begin: 0,
            shape_end: shapes.count(),
            code: SeedEncoding::Hashed,
            skip,
            filter_masked_seeds: false,
            mask_seeds: false,
            seed_cut,
            soft_masking,
            minimizer_window: 0,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: 0,
        };
        let ctx = crate::data::enum_seeds::EnumSeedsContext {
            shapes,
            reduction,
            min_query_len: 0,
            query_contexts: 1,
        };
        enum_seeds_block(seqs, &mut cb, &NO_FILTER, &cfg, &ctx)?;
        let sizes: Vec<usize> = cb[0].data.iter().map(Table::load).collect();
        cb[0].data = sizes
            .iter()
            .map(|&size| Table::new(next_pow2(size as f64 * HASH_TABLE_FACTOR)))
            .collect();
        enum_seeds_block(seqs, &mut cb, &NO_FILTER, &cfg, &ctx)?;
        for table in &mut cb[0].data {
            table.finish();
        }
        Ok(Self {
            data: cb.pop().unwrap().data,
        })
    }

    /// Matches C++ `HashedSeedSet::HashedSeedSet(const string&)`.
    pub fn from_index_file<P: AsRef<std::path::Path>>(
        index_file: P,
        shapes: &ShapeConfig,
    ) -> Result<Self, String> {
        let bytes = std::fs::read(index_file).map_err(|e| e.to_string())?;
        Self::from_index_bytes(&bytes, shapes)
    }

    /// Matches C++ `HashedSeedSet::HashedSeedSet(const string&)` after mmap.
    pub fn from_index_bytes(buf: &[u8], shapes: &ShapeConfig) -> Result<Self, String> {
        if buf.len() < SEED_INDEX_HEADER_SIZE {
            return Err("Invalid seed index file.".to_string());
        }
        let magic = u64::from_ne_bytes(buf[0..8].try_into().unwrap());
        if magic != SEED_INDEX_MAGIC_NUMBER {
            return Err("Invalid seed index file.".to_string());
        }
        let version = u32::from_ne_bytes(buf[8..12].try_into().unwrap());
        if version != SEED_INDEX_VERSION {
            return Err("Invalid seed index file version.".to_string());
        }
        let shape_count = i32::from_ne_bytes(buf[12..16].try_into().unwrap());
        if shape_count != shapes.count() {
            return Err("Index has a different number of shapes.".to_string());
        }

        let shape_count_usize = shape_count as usize;
        let sizes_begin = SEED_INDEX_HEADER_SIZE;
        let sizes_end = sizes_begin + std::mem::size_of::<usize>() * shape_count_usize;
        if buf.len() < sizes_end {
            return Err("Invalid seed index file.".to_string());
        }
        let mut sizes = Vec::with_capacity(shape_count_usize);
        for i in 0..shape_count_usize {
            let begin = sizes_begin + i * std::mem::size_of::<usize>();
            let end = begin + std::mem::size_of::<usize>();
            sizes.push(usize::from_ne_bytes(buf[begin..end].try_into().unwrap()));
        }

        let mut data = Vec::with_capacity(shape_count_usize);
        let mut data_ptr = sizes_end;
        for size in sizes {
            let end = data_ptr + size + Table::PADDING;
            if buf.len() < end {
                return Err("Invalid seed index file.".to_string());
            }
            data.push(Table::from_data(buf[data_ptr..end].to_vec(), size));
            data_ptr = end;
        }

        Ok(Self { data })
    }

    /// Matches C++ `HashedSeedSet::contains`.
    pub fn contains(&self, key: u64, shape: u64) -> bool {
        self.data[shape as usize].contains(key)
    }

    /// Matches C++ `HashedSeedSet::table`.
    pub fn table(&self, i: usize) -> &Table {
        &self.data[i]
    }

    /// Matches C++ `HashedSeedSet::max_table_size`.
    pub fn max_table_size(&self) -> usize {
        self.data.iter().map(Table::size).max().unwrap()
    }
}

struct SeedSetCallback {
    coverage: usize,
    max_coverage: usize,
    data: Vec<bool>,
}

impl SeedSetCallback {
    fn new(data: Vec<bool>, max_coverage: usize) -> Self {
        Self {
            coverage: 0,
            max_coverage,
            data,
        }
    }
}

impl EnumSeedsCallback for SeedSetCallback {
    fn call(&mut self, seed: u64, _pos: usize, _block_id: usize, _shape: i32) -> bool {
        if !self.data[seed as usize] {
            self.data[seed as usize] = true;
            self.coverage += 1;
            if self.coverage > self.max_coverage {
                return false;
            }
        }
        true
    }
}

struct HashedSeedSetCallback {
    data: Vec<Table>,
}

impl HashedSeedSetCallback {
    fn new(data: Vec<Table>) -> Self {
        Self { data }
    }
}

impl EnumSeedsCallback for HashedSeedSetCallback {
    fn call(&mut self, seed: u64, _pos: usize, _block_id: usize, shape: i32) -> bool {
        self.data[shape as usize].insert(seed);
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::seed_iterator::HashedSeedIterator;
    use crate::basic::value::{SequenceType, MASK_LETTER};

    #[test]
    fn test_seed_index_constants() {
        assert_eq!(SEED_INDEX_MAGIC_NUMBER, 0x2d6ba306ecbf6aba);
        assert_eq!(SEED_INDEX_VERSION, 0);
        assert_eq!(SEED_INDEX_HEADER_SIZE, 16);
        assert_eq!(HASH_TABLE_FACTOR, 1.25);
    }

    #[test]
    fn test_seed_set_accessors() {
        let seed_set = SeedSet {
            data: vec![false, true, false, true],
            coverage: 0.5,
        };
        assert!(!seed_set.contains(0, 9));
        assert!(seed_set.contains(1, 9));
        assert!(seed_set.contains(3, 0));
        assert_eq!(seed_set.coverage(), 0.5);
    }

    #[test]
    fn test_hashed_seed_set_accessors() {
        let mut a = Table::new(8);
        a.insert(3);
        a.insert(7);
        a.finish();
        let mut b = Table::new(8);
        b.insert(11);
        b.finish();
        let hashed = HashedSeedSet { data: vec![a, b] };

        assert!(hashed.contains(3, 0));
        assert!(!hashed.contains(3, 1));
        assert!(hashed.contains(11, 1));
        assert_eq!(hashed.table(0).load(), 2);
        assert_eq!(hashed.max_table_size(), 8);
    }

    #[test]
    fn test_seed_set_from_block_contiguous_and_coverage() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["11111".to_string()], 0, &reduction).unwrap();
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2, 3, 4, 5],
                Some("q1"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let seed_set = SeedSet::from_block(
            &mut block,
            1.0,
            None,
            0.0,
            MaskingAlgo::None,
            &shapes,
            &reduction,
        )
        .unwrap();
        let key0 = shapes[0]
            .set_seed_shifted(block.seqs().get(0), &reduction)
            .unwrap();
        let key1 = shapes[0]
            .set_seed_shifted(&block.seqs().get(0)[1..], &reduction)
            .unwrap();
        assert!(seed_set.contains(key0, 0));
        assert!(seed_set.contains(key1, 0));
        assert!(seed_set.coverage() > 0.0);
    }

    #[test]
    fn test_seed_set_from_block_rejects_spaced_shape() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["101".to_string()], 0, &reduction).unwrap();
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2],
                Some("q1"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let err = SeedSet::from_block(
            &mut block,
            1.0,
            None,
            0.0,
            MaskingAlgo::None,
            &shapes,
            &reduction,
        )
        .unwrap_err();
        assert_eq!(err, "Contiguous seed required.");
    }

    #[test]
    fn test_hashed_seed_set_from_block_multiple_shapes() {
        let reduction = Reduction::default_reduction();
        let shapes =
            ShapeConfig::from_codes(&["111".to_string(), "1011".to_string()], 0, &reduction)
                .unwrap();
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2, 3, MASK_LETTER],
                Some("q1"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let hashed = HashedSeedSet::from_block(
            &mut block,
            None,
            0.0,
            MaskingAlgo::None,
            &shapes,
            &reduction,
        )
        .unwrap();
        assert_eq!(hashed.data.len(), 2);
        let seq = block.seqs().get(0).to_vec();
        let it0 =
            HashedSeedIterator::<4>::new(&seq, block.seqs().length(0), shapes.get(0), &reduction);
        let it1 =
            HashedSeedIterator::<4>::new(&seq, block.seqs().length(0), shapes.get(1), &reduction);
        let key0 = it0.get();
        let key1 = it1.get();
        assert!(hashed.contains(key0, 0));
        assert!(hashed.contains(key1, 1));
        assert!(hashed.max_table_size() > 0);
    }

    #[test]
    fn test_hashed_seed_set_from_index_bytes() {
        let reduction = Reduction::default_reduction();
        let shapes =
            ShapeConfig::from_codes(&["111".to_string(), "1011".to_string()], 0, &reduction)
                .unwrap();
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2, 3, 4, 5],
                Some("q1"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        let hashed = HashedSeedSet::from_block(
            &mut block,
            None,
            0.0,
            MaskingAlgo::None,
            &shapes,
            &reduction,
        )
        .unwrap();

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&SEED_INDEX_MAGIC_NUMBER.to_ne_bytes());
        bytes.extend_from_slice(&SEED_INDEX_VERSION.to_ne_bytes());
        bytes.extend_from_slice(&(shapes.count() as u32).to_ne_bytes());
        for i in 0..shapes.count() as usize {
            bytes.extend_from_slice(&hashed.table(i).size().to_ne_bytes());
        }
        for i in 0..shapes.count() as usize {
            bytes.extend_from_slice(hashed.table(i).data());
        }

        let parsed = HashedSeedSet::from_index_bytes(&bytes, &shapes).unwrap();
        assert_eq!(parsed.max_table_size(), hashed.max_table_size());
        let seq = block.seqs().get(0).to_vec();
        let it0 =
            HashedSeedIterator::<4>::new(&seq, block.seqs().length(0), shapes.get(0), &reduction);
        let it1 =
            HashedSeedIterator::<4>::new(&seq, block.seqs().length(0), shapes.get(1), &reduction);
        assert!(parsed.contains(it0.get(), 0));
        assert!(parsed.contains(it1.get(), 1));
    }

    #[test]
    fn test_hashed_seed_set_from_index_bytes_rejects_invalid_headers() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        assert_eq!(
            HashedSeedSet::from_index_bytes(&[], &shapes).unwrap_err(),
            "Invalid seed index file."
        );

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&SEED_INDEX_MAGIC_NUMBER.to_ne_bytes());
        bytes.extend_from_slice(&(SEED_INDEX_VERSION + 1).to_ne_bytes());
        bytes.extend_from_slice(&(shapes.count() as u32).to_ne_bytes());
        assert_eq!(
            HashedSeedSet::from_index_bytes(&bytes, &shapes).unwrap_err(),
            "Invalid seed index file version."
        );
    }
}
