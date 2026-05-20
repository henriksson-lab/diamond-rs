use crate::basic::packed_loc::PackedLoc;
use crate::basic::sequence::Sequence;
use crate::basic::value::{BlockId, Loc, OId};
use crate::masking::MaskingAlgo;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeedEncoding {
    SpacedFactor,
    Hashed,
    Contiguous,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct NoFilter;

impl NoFilter {
    /// Matches C++ `NoFilter::contains(seed, shape)`.
    pub fn contains(&self, _seed: u64, _shape: u64) -> bool {
        true
    }
}

pub static NO_FILTER: NoFilter = NoFilter;

#[repr(C, packed)]
#[derive(Clone, Copy, Default, Eq)]
pub struct PackedLocId {
    pub pos: PackedLoc,
    pub block_id: u32,
}

impl PackedLocId {
    pub fn new(pos: PackedLoc) -> Self {
        Self { pos, block_id: 0 }
    }

    pub fn with_block_id(pos: PackedLoc, block_id: u32) -> Self {
        Self { pos, block_id }
    }

    /// Matches C++ `PackedLocId::operator uint64_t()`.
    pub fn as_u64(self) -> u64 {
        self.pos.as_u64()
    }
}

impl PartialEq for PackedLocId {
    fn eq(&self, rhs: &Self) -> bool {
        let pos = self.pos;
        let block_id = self.block_id;
        let rhs_pos = rhs.pos;
        let rhs_block_id = rhs.block_id;
        pos == rhs_pos && block_id == rhs_block_id
    }
}

impl std::fmt::Debug for PackedLocId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let pos = self.pos;
        let block_id = self.block_id;
        f.debug_struct("PackedLocId")
            .field("pos", &pos)
            .field("block_id", &block_id)
            .finish()
    }
}

impl From<PackedLocId> for u64 {
    fn from(value: PackedLocId) -> Self {
        value.as_u64()
    }
}

/// Matches C++ `block_id(i)`.
pub fn block_id_packed_loc_id(i: PackedLocId) -> u32 {
    i.block_id
}

/// Matches C++ `block_id(i)`.
pub fn block_id_packed_loc(_i: PackedLoc) -> u32 {
    panic!("Unsupported");
}

pub struct EnumCfg<'a> {
    pub partition: Option<&'a Vec<u32>>,
    pub shape_begin: i32,
    pub shape_end: i32,
    pub code: SeedEncoding,
    pub skip: Option<&'a Vec<bool>>,
    pub filter_masked_seeds: bool,
    pub mask_seeds: bool,
    pub seed_cut: f64,
    pub soft_masking: MaskingAlgo,
    pub minimizer_window: Loc,
    pub filter_low_complexity_seeds: bool,
    pub mask_low_complexity_seeds: bool,
    pub sketch_size: Loc,
}

pub struct SeqInfo<'a> {
    pub block_id: BlockId,
    pub oid: OId,
    pub title: Option<&'a str>,
    pub qual: Option<&'a str>,
    pub len: Loc,
    pub source_seq: Sequence<'a>,
    pub mate_seq: Sequence<'a>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_filter_contains() {
        assert!(NO_FILTER.contains(123, 7));
        assert!(NoFilter.contains(0, 0));
    }

    #[test]
    fn test_packed_loc_id_layout_and_accessors() {
        assert_eq!(std::mem::size_of::<PackedLocId>(), 9);
        assert_eq!(std::mem::align_of::<PackedLocId>(), 1);

        let p = PackedLoc::new(0x12_3456_789a);
        let x = PackedLocId::with_block_id(p, 42);
        assert_eq!(x.as_u64(), 0x12_3456_789a);
        assert_eq!(u64::from(x), 0x12_3456_789a);
        assert_eq!(block_id_packed_loc_id(x), 42);
        assert_eq!(PackedLocId::new(p).as_u64(), p.as_u64());
    }

    #[test]
    #[should_panic(expected = "Unsupported")]
    fn test_block_id_packed_loc_panics() {
        let _ = block_id_packed_loc(PackedLoc::new(1));
    }

    #[test]
    fn test_enum_cfg_and_seq_info_fields() {
        let partition = vec![0, 3];
        let skip = vec![false, true];
        let cfg = EnumCfg {
            partition: Some(&partition),
            shape_begin: 0,
            shape_end: 1,
            code: SeedEncoding::Contiguous,
            skip: Some(&skip),
            filter_masked_seeds: true,
            mask_seeds: false,
            seed_cut: 0.8,
            soft_masking: MaskingAlgo::None,
            minimizer_window: 0,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: 0,
        };
        assert_eq!(cfg.partition.unwrap()[1], 3);
        assert_eq!(cfg.skip.unwrap()[1], true);
        assert_eq!(cfg.code, SeedEncoding::Contiguous);

        let source = [0, 1, 2];
        let mate = [3, 4];
        let info = SeqInfo {
            block_id: 7,
            oid: 9,
            title: Some("title"),
            qual: None,
            len: 3,
            source_seq: Sequence::new(&source),
            mate_seq: Sequence::new(&mate),
        };
        assert_eq!(info.source_seq.length(), 3);
        assert_eq!(info.mate_seq.length(), 2);
    }
}
