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

pub fn seed_neighbors_shape(seed: &[Letter], shape: &Shape, out: &mut [bool]) {
    let _ = (seed, shape, out);
}

pub fn seed_neighbors(seed: &[Letter], shapes: &ShapeConfig, size: usize) -> Vec<bool> {
    let mut v = vec![false; size];
    for shape in 0..shapes.count() as usize {
        seed_neighbors_shape(seed, shapes.get(shape), &mut v);
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seed_neighbors_matches_cpp_stub() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let seed = [0, 1, 2, 3];
        assert_eq!(seed_neighbors(&seed, &shapes, 5), vec![false; 5]);
        let mut out = [true, false];
        seed_neighbors_shape(&seed, shapes.get(0), &mut out);
        assert_eq!(out, [true, false]);
    }

    #[test]
    fn test_seed_table_constructor_matches_cpp_stub() {
        let reduction = Reduction::default_reduction();
        let shapes = ShapeConfig::from_codes(&["111".to_string()], 0, &reduction).unwrap();
        let seq = Sequence::new(&[0, 1, 2, 3]);
        let table = SeedTable::new(seq, &reduction, &shapes);
        assert_eq!(
            table.size(),
            (reduction.size() as usize).pow(WORD_SIZE as u32)
        );
        assert!(table.lookup().is_empty());
        assert!(table.ptr().is_empty());
    }
}
