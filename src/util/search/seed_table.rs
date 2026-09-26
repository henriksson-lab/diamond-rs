//! Direct translation of `util/search/seed_table.cpp` and
//! `util/search/seed_table.h`.
//!
//! The upstream implementation is deliberately incomplete: neighbor marking
//! is a no-op and construction only records the six-letter table size. This
//! module preserves those semantics exactly instead of speculating about the
//! intended lookup layout.

use crate::basic::reduction::Reduction;
use crate::basic::sequence::Sequence;
use crate::basic::shape::Shape;
use crate::basic::shape_config::ShapeConfig;
use crate::basic::value::Letter;

pub const WORD_SIZE: usize = 6;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedTable {
    size: usize,
    lookup: Vec<u32>,
    ptr: Vec<u32>,
}

impl SeedTable {
    /// Matches the currently empty-body C++ `SeedTable::SeedTable`.
    pub fn new(seq: Sequence<'_>, reduction: &Reduction, shapes: &ShapeConfig) -> Self {
        let _ = (seq, shapes);
        Self {
            size: (reduction.size() as usize).pow(WORD_SIZE as u32),
            lookup: Vec::new(),
            ptr: Vec::new(),
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn lookup(&self) -> &[u32] {
        &self.lookup
    }

    pub fn ptr(&self) -> &[u32] {
        &self.ptr
    }
}

/// Matches the no-op shape overload of C++ `seed_neighbors`.
pub fn seed_neighbors_shape(seed: &[Letter], shape: &Shape, out: &mut [bool]) {
    let _ = (seed, shape, out);
}

/// Matches the shape-config overload: allocate an all-false vector, invoke the
/// no-op helper once per shape, and return it.
pub fn seed_neighbors(seed: &[Letter], shapes: &ShapeConfig, size: usize) -> Vec<bool> {
    let mut neighbors = vec![false; size];
    for shape in 0..shapes.count() as usize {
        seed_neighbors_shape(seed, shapes.get(shape), &mut neighbors);
    }
    neighbors
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::AMINO_ACID_ALPHABET;

    #[test]
    fn neighbor_helpers_preserve_upstream_no_op_behavior() {
        let reduction = Reduction::default_reduction();
        let shapes =
            ShapeConfig::from_codes(&["111".to_string(), "1011".to_string()], 0, &reduction)
                .unwrap();
        let seed = [0, 1, 2, 3];
        assert_eq!(seed_neighbors(&seed, &shapes, 5), vec![false; 5]);

        let mut existing = [true, false, true];
        seed_neighbors_shape(&seed, shapes.get(0), &mut existing);
        assert_eq!(existing, [true, false, true]);
    }

    #[test]
    fn constructor_only_sets_reduction_power_size() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let residues = [0, 1, 2, 3];
        let table = SeedTable::new(Sequence::new(&residues), &reduction, &shapes);
        assert_eq!(table.size(), 10usize.pow(WORD_SIZE as u32));
        assert!(table.lookup().is_empty());
        assert!(table.ptr().is_empty());

        let binary = Reduction::new("A RNDCQEGHILKMFPSTWYV", AMINO_ACID_ALPHABET);
        let empty_shapes = ShapeConfig::new();
        let binary_table = SeedTable::new(Sequence::new(&[]), &binary, &empty_shapes);
        assert_eq!(binary_table.size(), 2usize.pow(WORD_SIZE as u32));
        assert!(binary_table.lookup().is_empty());
        assert!(binary_table.ptr().is_empty());
    }

    #[test]
    fn compatibility_module_reexports_the_mirrored_type() {
        let reduction = Reduction::default_reduction();
        let table = crate::util::seed_table::SeedTable::new(
            Sequence::new(&[]),
            &reduction,
            &ShapeConfig::new(),
        );
        assert_eq!(table.size(), 10usize.pow(WORD_SIZE as u32));
    }
}
