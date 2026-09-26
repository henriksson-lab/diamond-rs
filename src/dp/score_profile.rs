//! Long score profiles.
//!
//! Direct Rust counterpart of `diamond/src/dp/score_profile.h` for scalar
//! profile construction.

use crate::basic::value::{letter_mask, Letter, AMINO_ACID_COUNT, TRUE_AA};
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::simd::Arch;

pub const DEFAULT_PADDING: usize = 128;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LongScoreProfile<Score: Copy + Default> {
    pub data: Vec<Vec<Score>>,
    pub padding: usize,
}

impl<Score: Copy + Default> LongScoreProfile<Score> {
    pub const DEFAULT_PADDING: usize = DEFAULT_PADDING;

    /// Matches C++ `LongScoreProfile::LongScoreProfile`.
    pub fn new<P>(padding: P) -> Self
    where
        P: TryInto<i64>,
    {
        let padding = match padding.try_into() {
            Ok(padding) => padding,
            Err(_) => panic!("score-profile padding does not fit int64"),
        };
        LongScoreProfile {
            data: vec![Vec::new(); AMINO_ACID_COUNT],
            padding: padding.max(Self::DEFAULT_PADDING as i64) as usize,
        }
    }

    /// Matches C++ `LongScoreProfile::length`.
    pub fn length(&self) -> usize {
        self.data
            .first()
            .map(|row| row.len().saturating_sub(2 * self.padding))
            .unwrap_or(0)
    }

    /// Compatibility accessor using an unsigned profile offset.
    pub fn get(&self, letter: Letter, i: usize) -> &[Score] {
        self.get_signed(letter, i as i64)
    }

    /// Matches C++ `LongScoreProfile::get(Letter, int)`, including access to
    /// left-padding cells through a negative offset.
    pub fn get_signed(&self, letter: Letter, i: i64) -> &[Score] {
        let row = &self.data[letter as usize];
        let pos = i + self.padding as i64;
        assert!(pos >= 0, "score-profile offset precedes left padding");
        &row[pos as usize..]
    }

    /// Matches C++ `LongScoreProfile::pointers`.
    pub fn pointers(&self, offset: usize) -> Vec<&[Score]> {
        self.pointers_signed(offset as i64)
    }

    /// Signed-offset counterpart of C++ `pointers(int)`.
    pub fn pointers_signed(&self, offset: i64) -> Vec<&[Score]> {
        let mut v = Vec::with_capacity(AMINO_ACID_COUNT);
        for letter in 0..AMINO_ACID_COUNT {
            v.push(self.get_signed(letter as Letter, offset));
        }
        v
    }

    /// Matches C++ `LongScoreProfile::reverse` (including both padding areas).
    pub fn reverse(&self) -> Self {
        let mut r = self.clone();
        for row in &mut r.data {
            row.reverse();
        }
        r
    }
}

/// Rust counterpart of C++'s templated score-matrix `make_profile` overload.
///
/// C++ has an AVX2-only specialization inside this function. Its byte-vector
/// `operator+=` saturates and skips CBS for ambiguous-letter rows. All other
/// dispatched architectures use the scalar loop, whose `int8_t +=` narrows by
/// wrapping on DIAMOND's supported compilers and applies CBS to every row.
fn make_profile<Score>(
    seq: &[Letter],
    cbs: Option<&[i8]>,
    padding: usize,
    matrix: &ScoreMatrix,
) -> LongScoreProfile<Score>
where
    Score: Copy + Default + From<i8>,
{
    // Standard C++ builds leave WITH_AVX512 disabled, so an AVX-512-capable
    // host is reported as AVX2 and takes this specialization. Rust's global
    // detector reports the hardware capability independently of build flags,
    // hence both Rust variants map to the C++ AVX2 path here.
    let avx2_specialization = matches!(crate::util::simd::arch(), Arch::Avx2 | Arch::Avx512);
    make_profile_for_dispatch(seq, cbs, padding, matrix, avx2_specialization)
}

fn make_profile_for_dispatch<Score>(
    seq: &[Letter],
    cbs: Option<&[i8]>,
    padding: usize,
    matrix: &ScoreMatrix,
    avx2_specialization: bool,
) -> LongScoreProfile<Score>
where
    Score: Copy + Default + From<i8>,
{
    if let Some(cbs) = cbs {
        assert!(
            cbs.len() >= seq.len(),
            "CBS vector must cover the complete sequence"
        );
    }

    let mut profile = LongScoreProfile::new(padding);
    let len = seq.len() + 2 * profile.padding;
    for row in &mut profile.data {
        row.resize(len, Score::from(-1));
    }

    for letter in 0..AMINO_ACID_COUNT {
        let scores = &matrix.matrix8()[letter << 5..(letter + 1) << 5];
        for (i, &subject) in seq.iter().enumerate() {
            let mut score = scores[letter_mask(subject) as usize];
            if let Some(cbs) = cbs {
                if !avx2_specialization || letter < TRUE_AA as usize {
                    score = if avx2_specialization {
                        score.saturating_add(cbs[i])
                    } else {
                        score.wrapping_add(cbs[i])
                    };
                }
            }
            profile.data[letter][i + profile.padding] = Score::from(score);
        }
    }
    profile
}

pub fn make_profile8(
    seq: &[Letter],
    cbs: Option<&[i8]>,
    padding: usize,
    matrix: &ScoreMatrix,
) -> LongScoreProfile<i8> {
    make_profile(seq, cbs, padding, matrix)
}

pub fn make_profile16(
    seq: &[Letter],
    cbs: Option<&[i8]>,
    padding: usize,
    matrix: &ScoreMatrix,
) -> LongScoreProfile<i16> {
    make_profile(seq, cbs, padding, matrix)
}

/// Compatibility name retained for existing Rust callers.
pub fn make_profile16_target_matrix(
    seq: &[Letter],
    matrix: &TargetMatrix,
    padding: usize,
) -> LongScoreProfile<i16> {
    make_profile16_from_target_matrix(seq, matrix, padding)
}

/// Rust mapping of the `Stats::TargetMatrix` overload of C++ `make_profile16`.
pub fn make_profile16_from_target_matrix(
    seq: &[Letter],
    matrix: &TargetMatrix,
    padding: usize,
) -> LongScoreProfile<i16> {
    make_profile_from_target_matrix(seq, matrix, padding)
}

/// Rust counterpart of C++'s templated target-matrix `make_profile` overload.
fn make_profile_from_target_matrix<Score>(
    seq: &[Letter],
    matrix: &TargetMatrix,
    padding: usize,
) -> LongScoreProfile<Score>
where
    Score: Copy + Default + From<i8>,
{
    let mut profile = LongScoreProfile::new(padding);
    let len = seq.len() + 2 * profile.padding;
    for row in &mut profile.data {
        row.resize(len, Score::from(-1));
    }
    for letter in 0..AMINO_ACID_COUNT {
        let row = &matrix.scores[letter << 5..(letter + 1) << 5];
        for (i, &subject) in seq.iter().enumerate() {
            profile.data[letter][i + profile.padding] =
                Score::from(row[letter_mask(subject) as usize]);
        }
    }
    profile
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SEED_MASK;

    #[test]
    fn test_long_score_profile_length_and_get() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq: Vec<Letter> = vec![0, 1, 2, 3];
        let profile = make_profile16(&seq, None, 4, &sm);
        assert_eq!(profile.padding, DEFAULT_PADDING);
        assert_eq!(profile.length(), seq.len());
        assert_eq!(profile.get(0, 0)[0], sm.score(0, seq[0]) as i16);
        assert_eq!(profile.get_signed(0, -1)[0], -1);
        assert_eq!(profile.get_signed(0, seq.len() as i64)[0], -1);
        let pointers = profile.pointers_signed(-1);
        assert_eq!(pointers.len(), AMINO_ACID_COUNT);
        assert!(pointers.iter().all(|row| row[0] == -1));

        let negative_padding = LongScoreProfile::<i8>::new(-50_i64);
        assert_eq!(
            negative_padding.padding,
            LongScoreProfile::<i8>::DEFAULT_PADDING
        );
    }

    #[test]
    fn test_long_score_profile_cbs_and_reverse() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq: Vec<Letter> = vec![0, 1, 2];
        let cbs = vec![1, -2, 3];
        let profile = make_profile8(&seq, Some(&cbs), 128, &sm);
        assert_eq!(profile.get(0, 0)[0], (sm.score(0, seq[0]) + 1) as i8);
        let rev = profile.reverse();
        assert_eq!(rev.data[0].first(), profile.data[0].last());
        assert_eq!(rev.get(0, 0)[0], profile.get(0, seq.len() - 1)[0]);
    }

    #[test]
    fn test_cbs_byte_arithmetic_matches_both_cpp_dispatch_paths() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq = vec![SEED_MASK]; // masked A; C++ Sequence/ScoreVector strips the bit
        let cbs = vec![127];

        for avx2 in [false, true] {
            let p8 = make_profile_for_dispatch::<i8>(&seq, Some(&cbs), 0, &sm, avx2);
            let p16 = make_profile_for_dispatch::<i16>(&seq, Some(&cbs), 0, &sm, avx2);
            for letter in 0..AMINO_ACID_COUNT {
                let base = sm.matrix8()[letter << 5];
                let expected = if avx2 && letter >= TRUE_AA as usize {
                    base
                } else if avx2 {
                    base.saturating_add(cbs[0])
                } else {
                    base.wrapping_add(cbs[0])
                };
                assert_eq!(p8.get(letter as Letter, 0)[0], expected, "row {letter}");
                // C++ performs byte arithmetic before `store_expanded` converts
                // the AVX2 result to int16_t.
                assert_eq!(
                    p16.get(letter as Letter, 0)[0],
                    expected as i16,
                    "row {letter}"
                );
            }
        }

        let runtime = make_profile16(&seq, Some(&cbs), 0, &sm);
        let uses_avx2 = matches!(crate::util::simd::arch(), Arch::Avx2 | Arch::Avx512);
        let forced = make_profile_for_dispatch::<i16>(&seq, Some(&cbs), 0, &sm, uses_avx2);
        assert_eq!(runtime, forced);
    }

    #[test]
    fn test_empty_profile_has_cpp_padding_and_zero_logical_length() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let profile = make_profile8(&[], None, 0, &sm);
        assert_eq!(profile.length(), 0);
        assert_eq!(profile.data.len(), AMINO_ACID_COUNT);
        assert!(profile
            .data
            .iter()
            .all(|row| row.len() == 2 * DEFAULT_PADDING && row.iter().all(|&x| x == -1)));
    }

    #[test]
    fn test_make_profile16_target_matrix() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let q = crate::stats::cbs::compute_composition(&[0, 1, 2, 3, 4]);
        let t = crate::stats::cbs::compute_composition(&[0, 1, 1, 2, 3]);
        let matrix = TargetMatrix::from_hauser_global(&q, &t, &sm);
        let seq: Vec<Letter> = vec![0, 1, 2];
        let profile = make_profile16_target_matrix(&seq, &matrix, 4);
        let canonical = make_profile16_from_target_matrix(&seq, &matrix, 4);
        assert_eq!(profile.length(), seq.len());
        assert_eq!(profile.get(0, 0)[0], matrix.scores[0] as i16);
        assert_eq!(profile, canonical);

        let masked = vec![seq[0] | SEED_MASK];
        let masked_profile = make_profile16_from_target_matrix(&masked, &matrix, 4);
        assert_eq!(
            masked_profile.get(0, 0)[0],
            matrix.scores[seq[0] as usize] as i16
        );
    }
}
