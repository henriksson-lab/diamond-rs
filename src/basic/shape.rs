use super::reduction::Reduction;
use super::seed::{PackedSeed, MAX_SEED_WEIGHT};
use super::value::{is_amino_acid, Letter, LETTER_MASK, MASK_LETTER, STOP_LETTER, DELIMITER_LETTER};

/// A spaced seed shape/pattern used for seed extraction.
///
/// Pattern is a string of '1' and '0' characters, where '1' indicates
/// a position that contributes to the seed value.
#[derive(Default)]
pub struct Shape {
    pub length: i32,
    pub weight: i32,
    pub positions: [u32; MAX_SEED_WEIGHT],
    pub d: u32,
    pub mask: u32,
    pub rev_mask: u32,
    pub long_mask: u64,
}

impl Shape {
    /// Create an empty shape.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a shape from a pattern code string (e.g., "111101011").
    pub fn from_code(code: &str, reduction: &Reduction) -> Self {
        let b = reduction.bit_size() as u64;
        let mut shape = Shape::new();

        for (i, ch) in code.bytes().enumerate() {
            shape.rev_mask <<= 1;
            shape.long_mask <<= b;
            if ch == b'1' {
                assert!((shape.weight as usize) < MAX_SEED_WEIGHT);
                shape.positions[shape.weight as usize] = i as u32;
                shape.weight += 1;
                shape.mask |= 1 << i;
                shape.rev_mask |= 1;
                shape.long_mask |= (1u64 << b) - 1;
            }
        }
        shape.length = code.len() as i32;
        if shape.weight >= 2 {
            shape.d = shape.positions[(shape.weight / 2 - 1) as usize];
        }
        shape
    }

    /// Extract a packed seed from a sequence at the current position.
    /// Returns None if any position contains a non-amino-acid letter.
    #[inline]
    pub fn set_seed(&self, seq: &[Letter], reduction: &Reduction) -> Option<PackedSeed> {
        let mut s: PackedSeed = 0;
        for i in 0..self.weight as usize {
            let l = seq[self.positions[i] as usize] & LETTER_MASK;
            if !is_amino_acid(l) {
                return None;
            }
            let r = reduction.reduce(l);
            s *= reduction.size() as u64;
            s += r as u64;
        }
        Some(s)
    }

    /// Extract a packed seed using bit-shifting.
    #[inline]
    pub fn set_seed_shifted(&self, seq: &[Letter], reduction: &Reduction) -> Option<PackedSeed> {
        let mut s: PackedSeed = 0;
        let b = reduction.bit_size() as u64;
        for i in 0..self.weight as usize {
            let l = seq[self.positions[i] as usize] & LETTER_MASK;
            if l == MASK_LETTER || l == DELIMITER_LETTER || l == STOP_LETTER {
                return None;
            }
            let r = reduction.reduce(l);
            s <<= b;
            s |= r as u64;
        }
        Some(s)
    }

    /// Extract a packed seed from a pre-reduced sequence.
    #[inline]
    pub fn set_seed_reduced(&self, seq: &[Letter], reduction: &Reduction) -> Option<PackedSeed> {
        let mut s: PackedSeed = 0;
        for i in 0..self.weight as usize {
            let l = seq[self.positions[i] as usize] & LETTER_MASK;
            if l == MASK_LETTER {
                return None;
            }
            s *= reduction.size() as u64;
            s += l as u64;
        }
        Some(s)
    }

    /// Whether the pattern is contiguous (no gaps).
    pub fn contiguous(&self) -> bool {
        self.length == self.weight
    }

    /// Number of bits needed to represent the seed value.
    pub fn bit_length(&self, reduction: &Reduction) -> i32 {
        let max_val = (reduction.size() as i64).pow(self.weight as u32) - 1;
        if max_val <= 0 {
            0
        } else {
            64 - max_val.leading_zeros() as i32
        }
    }
}

impl std::fmt::Display for Shape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.length {
            if self.mask & (1 << i) != 0 {
                write!(f, "1")?;
            } else {
                write!(f, "0")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_contiguous() {
        let r = Reduction::default_reduction();
        let s = Shape::from_code("1111", &r);
        assert!(s.contiguous());
        assert_eq!(s.length, 4);
        assert_eq!(s.weight, 4);
    }

    #[test]
    fn test_shape_spaced() {
        let r = Reduction::default_reduction();
        let s = Shape::from_code("110101", &r);
        assert!(!s.contiguous());
        assert_eq!(s.length, 6);
        assert_eq!(s.weight, 4);
        assert_eq!(s.positions[0], 0);
        assert_eq!(s.positions[1], 1);
        assert_eq!(s.positions[2], 3);
        assert_eq!(s.positions[3], 5);
    }

    #[test]
    fn test_shape_display() {
        let r = Reduction::default_reduction();
        let s = Shape::from_code("110101", &r);
        assert_eq!(format!("{}", s), "110101");
    }

    #[test]
    fn test_shape_set_seed() {
        let r = Reduction::default_reduction();
        let s = Shape::from_code("111", &r);
        // Sequence: A(0), R(1), N(2)
        let seq = vec![0i8, 1, 2];
        let seed = s.set_seed(&seq, &r);
        assert!(seed.is_some());
    }

    #[test]
    fn test_shape_set_seed_masked() {
        let r = Reduction::default_reduction();
        let s = Shape::from_code("111", &r);
        // Sequence with mask letter should return None
        let seq = vec![0i8, MASK_LETTER, 2];
        let seed = s.set_seed(&seq, &r);
        assert!(seed.is_none());
    }
}
