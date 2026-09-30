use super::reduction::Reduction;
use super::seed::{PackedSeed, SeedPartition, MAX_SEED_WEIGHT};

/// `ln(n!)` for small n. The Hauser/SEG-style entropy formula uses this in a
/// tight loop. Keep this as a literal table instead of a lazily initialized
/// `OnceLock`: seed enumeration calls it for every candidate seed, so even the
/// fast initialized path shows up in real blastp profiles.
#[inline]
fn ln_factorial(n: u32) -> f64 {
    const TABLE: [f64; 32] = [
        0.0,
        0.0,
        0.6931471805599453,
        1.791759469228055,
        3.1780538303479453,
        4.787491742782046,
        6.579251212010101,
        8.525161361065415,
        10.60460290274525,
        12.80182748008147,
        15.104412573075518,
        17.502307845873887,
        19.98721449566189,
        22.552163853123425,
        25.191221182738683,
        27.89927138384089,
        30.671860106080675,
        33.50507345013689,
        36.39544520803305,
        39.339884187199495,
        42.335616460753485,
        45.38013889847691,
        48.47118135183523,
        51.60667556776438,
        54.784729398112326,
        58.00360522298052,
        61.261701761002,
        64.55753862700634,
        67.88974313718154,
        71.257038967168,
        74.65823634883017,
        78.0922235533153,
    ];
    if let Some(&x) = TABLE.get(n as usize) {
        x
    } else {
        (2..=n).map(|i| (i as f64).ln()).sum()
    }
}
use super::value::{
    is_amino_acid, Letter, DELIMITER_LETTER, LETTER_MASK, MASK_LETTER, STOP_LETTER,
};

/// A spaced seed shape/pattern used for seed extraction.
///
/// Pattern is a string of '1' and '0' characters, where '1' indicates
/// a position that contributes to the seed value.
#[derive(Clone, Copy, Debug, Default)]
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
        // All production protein-search shapes currently have ten selected
        // positions under the Murphy-10 reduction.  Spell that case at its
        // fixed width, just as the native compiler does for DIAMOND's hot seed
        // iterator: this removes the loop-carried base-10 multiply chain and
        // combines two reduction lookups in one small L1-resident table.
        // Keep the generic path for custom reductions and test shapes.
        if self.weight == 10 && reduction.size() == 10 {
            debug_assert!(seq.len() >= self.length as usize);
            macro_rules! letter {
                ($index:expr) => {{
                    // SAFETY: Shape::from_code records selected positions
                    // inside `length`, and callers provide a complete window.
                    let position = unsafe { *self.positions.get_unchecked($index) as usize };
                    let raw = unsafe { *seq.get_unchecked(position) };
                    if raw & crate::basic::value::SEED_MASK != 0 {
                        return None;
                    }
                    raw & LETTER_MASK
                }};
            }
            let p0 = reduction.reduce_pair10(letter!(0), letter!(1));
            let p1 = reduction.reduce_pair10(letter!(2), letter!(3));
            let p2 = reduction.reduce_pair10(letter!(4), letter!(5));
            let p3 = reduction.reduce_pair10(letter!(6), letter!(7));
            let p4 = reduction.reduce_pair10(letter!(8), letter!(9));
            if p0 == u16::MAX
                || p1 == u16::MAX
                || p2 == u16::MAX
                || p3 == u16::MAX
                || p4 == u16::MAX
            {
                return None;
            }
            return Some(
                p0 as u64 * 100_000_000
                    + p1 as u64 * 1_000_000
                    + p2 as u64 * 10_000
                    + p3 as u64 * 100
                    + p4 as u64,
            );
        }
        let mut s: PackedSeed = 0;
        for i in 0..self.weight as usize {
            let raw = seq[self.positions[i] as usize];
            // Match C++'s `soft_masking` behavior: any soft-masked position
            // (high bit set, from tantan) makes the entire seed invalid for
            // indexing. Mirrors what `EnumCfg::soft_masking` triggers in
            // C++ DIAMOND seed enumeration.
            if (raw & crate::basic::value::SEED_MASK) != 0 {
                return None;
            }
            let l = raw & LETTER_MASK;
            if !is_amino_acid(l) {
                return None;
            }
            let r = reduction.reduce(l);
            s *= reduction.size() as u64;
            s += r as u64;
        }
        Some(s)
    }

    /// Extract a seed after the caller has proved that the complete shape fits.
    ///
    /// # Safety
    /// `seq.len()` must be at least `self.length`, and `self` must retain the
    /// position/weight invariants established by [`Shape::from_code`].
    #[inline(always)]
    pub unsafe fn set_seed_unchecked(
        &self,
        seq: &[Letter],
        reduction: &Reduction,
    ) -> Option<PackedSeed> {
        debug_assert!(seq.len() >= self.length as usize);
        debug_assert!((self.weight as usize) <= MAX_SEED_WEIGHT);
        let mut seed: PackedSeed = 0;
        let size = reduction.size() as u64;
        for i in 0..self.weight as usize {
            // SAFETY: both unchecked indices are covered by this method's
            // contract and Shape's constructor invariants.
            let position = unsafe { *self.positions.get_unchecked(i) as usize };
            let raw = unsafe { *seq.get_unchecked(position) };
            // This path mirrors C++ `Shape::set_seed` under `SEQ_MASK`: soft
            // masking is stripped before the amino-acid validity check. Seed
            // enumeration uses the checked method above, which deliberately
            // rejects soft-masked positions.
            let letter = raw & LETTER_MASK;
            if !is_amino_acid(letter) {
                return None;
            }
            seed = seed * size + reduction.reduce(letter) as u64;
        }
        Some(seed)
    }

    /// Compute only the partition bits needed by chunked left-most checking.
    ///
    /// The production protein shapes have weight 10 with the Murphy-10
    /// reduction. Spell that case at its fixed width so the compiler can
    /// schedule the independent constant-weight terms instead of retaining a
    /// ten-step multiply dependency chain. The result is algebraically
    /// identical to `set_seed_unchecked(...) & partition_mask`.
    ///
    /// # Safety
    /// The requirements are identical to [`Shape::set_seed_unchecked`].
    #[inline(always)]
    pub unsafe fn seed_partition_unchecked(
        &self,
        seq: &[Letter],
        reduction: &Reduction,
        partition_mask: PackedSeed,
    ) -> Option<SeedPartition> {
        if self.weight != 10 || reduction.size() != 10 {
            return self
                .set_seed_unchecked(seq, reduction)
                .map(|seed| (seed & partition_mask) as SeedPartition);
        }
        debug_assert!(seq.len() >= self.length as usize);
        macro_rules! letter {
            ($index:expr) => {{
                let position = *self.positions.get_unchecked($index) as usize;
                *seq.get_unchecked(position) & LETTER_MASK
            }};
        }
        let l0 = letter!(0);
        let l1 = letter!(1);
        let l2 = letter!(2);
        let l3 = letter!(3);
        let l4 = letter!(4);
        let l5 = letter!(5);
        let l6 = letter!(6);
        let l7 = letter!(7);
        let l8 = letter!(8);
        let l9 = letter!(9);
        let p0 = reduction.reduce_pair10(l0, l1);
        let p1 = reduction.reduce_pair10(l2, l3);
        let p2 = reduction.reduce_pair10(l4, l5);
        let p3 = reduction.reduce_pair10(l6, l7);
        let p4 = reduction.reduce_pair10(l8, l9);
        if p0 == u16::MAX || p1 == u16::MAX || p2 == u16::MAX || p3 == u16::MAX || p4 == u16::MAX {
            return None;
        }
        let seed = p0 as u64 * 100_000_000
            + p1 as u64 * 1_000_000
            + p2 as u64 * 10_000
            + p3 as u64 * 100
            + p4 as u64;
        Some((seed & partition_mask) as SeedPartition)
    }

    /// Fused seed extraction + complexity check. Returns the seed value if
    /// the seed is valid AND its REDUCED-alphabet composition has entropy
    /// ≥ `cut`. Iterates the spaced positions once, accumulating both the
    /// seed value and a REDUCED-AA histogram, then computes entropy from
    /// the histogram. Replaces a `seed_is_complex` + `set_seed` pair that
    /// walked the same spaced positions twice.
    ///
    /// Matches C++ `Search::seed_is_complex(seed, shape, cut)`.
    /// the entropy is over the Murphy-10 reduced alphabet, NOT the full
    /// 20-letter unreduced one. The reduced filter is strictly more
    /// aggressive on tandem repeats (e.g. "QQQNNQLT"-style sequences) —
    /// using unreduced counts here lets through seeds C++ would drop in
    /// its post-join `mask_seeds` step, producing Rust-only hits.
    #[inline]
    pub fn set_seed_with_complexity(
        &self,
        seq: &[Letter],
        reduction: &Reduction,
        complexity_cut: f64,
    ) -> Option<PackedSeed> {
        use crate::basic::value::TRUE_AA;
        let mut s: PackedSeed = 0;
        let red_size = reduction.size() as usize;
        // Reduction has at most 10 classes (Murphy-10). 16 is a safe upper bound.
        let mut counts = [0u32; 16];
        let weight = self.weight as usize;
        let size = reduction.size() as u64;
        for i in 0..weight {
            let raw = seq[self.positions[i] as usize];
            if (raw & crate::basic::value::SEED_MASK) != 0 {
                return None;
            }
            let l = raw & LETTER_MASK;
            if (l as i32) >= TRUE_AA {
                return None;
            }
            let r = reduction.reduce(l);
            counts[r as usize] += 1;
            s = s * size + r as u64;
        }
        if complexity_cut > 0.0 {
            let mut entropy = ln_factorial(weight as u32);
            for &c in counts[..red_size].iter() {
                if c > 0 {
                    entropy -= ln_factorial(c);
                }
            }
            if entropy < complexity_cut {
                return None;
            }
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
    #[inline(always)]
    pub fn set_seed_reduced(&self, seq: &[Letter], reduction: &Reduction) -> Option<PackedSeed> {
        if reduction.size() == 10 && matches!(self.weight, 7 | 8 | 10) {
            debug_assert!(seq.len() >= self.length as usize);
            macro_rules! digit {
                ($index:expr) => {{
                    // SAFETY: selected positions are bounded by Shape::length.
                    let position = unsafe { *self.positions.get_unchecked($index) as usize };
                    let letter = unsafe { *seq.get_unchecked(position) } & LETTER_MASK;
                    if letter == MASK_LETTER {
                        return None;
                    }
                    letter as u64
                }};
            }
            return match self.weight {
                7 => {
                    let p0 = digit!(0) * 10 + digit!(1);
                    let p1 = digit!(2) * 10 + digit!(3);
                    let p2 = digit!(4) * 10 + digit!(5);
                    Some((p0 * 10_000 + p1 * 100 + p2) * 10 + digit!(6))
                }
                8 => {
                    let p0 = digit!(0) * 10 + digit!(1);
                    let p1 = digit!(2) * 10 + digit!(3);
                    let p2 = digit!(4) * 10 + digit!(5);
                    let p3 = digit!(6) * 10 + digit!(7);
                    Some(p0 * 1_000_000 + p1 * 10_000 + p2 * 100 + p3)
                }
                10 => {
                    let p0 = digit!(0) * 10 + digit!(1);
                    let p1 = digit!(2) * 10 + digit!(3);
                    let p2 = digit!(4) * 10 + digit!(5);
                    let p3 = digit!(6) * 10 + digit!(7);
                    let p4 = digit!(8) * 10 + digit!(9);
                    Some(p0 * 100_000_000 + p1 * 1_000_000 + p2 * 10_000 + p3 * 100 + p4)
                }
                _ => unreachable!(),
            };
        }
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

    #[test]
    fn specialized_weight10_partition_matches_full_seed() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("111101110111", &reduction);
        assert_eq!(shape.weight, 10);
        let mut state = 0x5eed_10u64;
        for _ in 0..512 {
            let sequence: Vec<Letter> = (0..shape.length)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let letter = ((state >> 32) & 31) as Letter;
                    if state & 7 == 0 {
                        letter | crate::basic::value::SEED_MASK
                    } else {
                        letter
                    }
                })
                .collect();
            let full = unsafe { shape.set_seed_unchecked(&sequence, &reduction) };
            for bits in [4, 8, 10, 16] {
                let mask = (1u64 << bits) - 1;
                let expected = full.map(|seed| (seed & mask) as SeedPartition);
                let actual = unsafe { shape.seed_partition_unchecked(&sequence, &reduction, mask) };
                assert_eq!(actual, expected);
            }
        }
    }
}
