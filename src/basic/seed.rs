//! Expanded seed operations translated from `diamond/src/basic/basic.cpp` and
//! `diamond/src/basic/seed.h`.  The C++ implementation reads `seed_weight` and
//! the score matrix from process-wide configuration; the Rust API takes both
//! explicitly so callers cannot accidentally enumerate with stale globals.

use super::value::Letter;
use crate::stats::score_matrix::ScoreMatrix;

/// Packed seed representation (uint64_t in C++).
pub type PackedSeed = u64;

/// Seed offset type for indexing.
pub type SeedOffset = u32;

/// Seed partition type.
pub type SeedPartition = u32;

/// Maximum seed weight (number of positions in a spaced seed pattern).
pub const MAX_SEED_WEIGHT: usize = 32;

#[inline]
pub fn seedp_mask(seedp_bits: i32) -> PackedSeed {
    (1u64 << seedp_bits) - 1
}

#[inline]
pub fn seedp_count(seedp_bits: i32) -> PackedSeed {
    1u64 << seedp_bits
}

#[inline]
pub fn seed_partition(s: PackedSeed, mask: PackedSeed) -> SeedPartition {
    (s & mask) as SeedPartition
}

#[inline]
pub fn seed_partition_offset(s: PackedSeed, seedp_bits: PackedSeed) -> SeedOffset {
    (s >> seedp_bits) as SeedOffset
}

/// A seed as an array of letter values (expanded form).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Seed {
    data: [Letter; MAX_SEED_WEIGHT],
}

impl Seed {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn get(&self, i: usize) -> Letter {
        self.data[i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, val: Letter) {
        self.data[i] = val;
    }

    /// Convert to packed representation using alphabet size 20.
    pub fn packed(&self, weight: usize) -> u64 {
        assert!(weight <= MAX_SEED_WEIGHT, "seed weight exceeds storage");
        let mut s = 0u64;
        for i in 0..weight {
            s *= 20;
            s += self.data[i] as u64;
        }
        s
    }

    /// Score this seed against `rhs` over the configured seed weight.
    ///
    /// This is C++ `Seed::score`, with its implicit `config.seed_weight` and
    /// global `score_matrix` made explicit.
    pub fn score(&self, rhs: &Seed, weight: usize, score_matrix: &ScoreMatrix) -> i32 {
        assert!(weight <= MAX_SEED_WEIGHT, "seed weight exceeds storage");
        (0..weight)
            .map(|i| score_matrix.score(self.data[i], rhs.data[i]))
            .sum()
    }

    /// Enumerate all true-amino-acid seeds scoring at least `threshold`.
    ///
    /// Results have the same base-20/lexicographic order as the recursive C++
    /// implementation. `self` is not modified, unlike the temporarily-mutated
    /// C++ receiver.
    pub fn enum_neighborhood(
        &self,
        threshold: i32,
        weight: usize,
        score_matrix: &ScoreMatrix,
    ) -> Vec<Seed> {
        assert!(weight > 0, "seed weight must be positive");
        assert!(weight <= MAX_SEED_WEIGHT, "seed weight exceeds storage");

        let mut candidate = self.clone();
        let mut out = Vec::new();
        let initial_score = self.score(self, weight, score_matrix);
        self.enum_neighborhood_at(
            0,
            threshold,
            weight,
            initial_score,
            score_matrix,
            &mut candidate,
            &mut out,
        );
        out
    }

    #[allow(clippy::too_many_arguments)]
    fn enum_neighborhood_at(
        &self,
        pos: usize,
        threshold: i32,
        weight: usize,
        score: i32,
        score_matrix: &ScoreMatrix,
        candidate: &mut Seed,
        out: &mut Vec<Seed>,
    ) {
        let original = self.data[pos];
        let score_without_original = score - score_matrix.score(original, original);
        for replacement in 0..20 {
            let new_score =
                score_without_original + score_matrix.score(original, replacement as Letter);
            candidate.data[pos] = replacement as Letter;
            if new_score >= threshold {
                if pos + 1 < weight {
                    self.enum_neighborhood_at(
                        pos + 1,
                        threshold,
                        weight,
                        new_score,
                        score_matrix,
                        candidate,
                        out,
                    );
                } else {
                    out.push(candidate.clone());
                }
            }
        }
        candidate.data[pos] = original;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blosum62() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, -1, 1, 1_000).unwrap()
    }

    #[test]
    fn test_seedp_mask() {
        assert_eq!(seedp_mask(4), 0xF);
        assert_eq!(seedp_mask(8), 0xFF);
    }

    #[test]
    fn test_seed_partition() {
        let s: PackedSeed = 0xABCD;
        assert_eq!(seed_partition(s, 0xFF), 0xCD);
    }

    #[test]
    fn test_seed_packed() {
        let mut seed = Seed::new();
        seed.set(0, 1);
        seed.set(1, 2);
        // packed = 1*20 + 2 = 22
        assert_eq!(seed.packed(2), 22);
    }

    #[test]
    fn test_seed_score_uses_explicit_weight_and_matrix() {
        let matrix = blosum62();
        let mut a = Seed::new();
        let mut b = Seed::new();
        a.set(0, 0); // A
        a.set(1, 1); // R
        b.set(0, 0); // A
        b.set(1, 2); // N
        assert_eq!(a.score(&b, 1, &matrix), matrix.score(0, 0));
        assert_eq!(
            a.score(&b, 2, &matrix),
            matrix.score(0, 0) + matrix.score(1, 2)
        );
    }

    #[test]
    fn test_enum_neighborhood_matches_exhaustive_base20_order() {
        let matrix = blosum62();
        let mut seed = Seed::new();
        seed.set(0, 0); // A
        seed.set(1, 1); // R
        let threshold = 5;

        let actual = seed.enum_neighborhood(threshold, 2, &matrix);
        let mut expected = Vec::new();
        for a in 0..20 {
            for b in 0..20 {
                if matrix.score(0, a) + matrix.score(1, b) >= threshold {
                    let mut candidate = Seed::new();
                    candidate.set(0, a);
                    candidate.set(1, b);
                    expected.push(candidate);
                }
            }
        }
        assert_eq!(actual, expected);
        assert!(actual.contains(&seed));
        assert_eq!(seed.get(0), 0);
        assert_eq!(seed.get(1), 1);
    }
}
