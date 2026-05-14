//! Long score profiles.
//!
//! Direct Rust counterpart of `diamond/src/dp/score_profile.h` for scalar
//! profile construction.

use crate::basic::value::{Letter, AMINO_ACID_COUNT};
use crate::stats::cbs::TargetMatrix;
use crate::stats::score_matrix::ScoreMatrix;

pub const DEFAULT_PADDING: usize = 128;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LongScoreProfile<Score: Copy + Default> {
    pub data: Vec<Vec<Score>>,
    pub padding: usize,
}

impl<Score: Copy + Default> LongScoreProfile<Score> {
    pub fn new(padding: usize) -> Self {
        LongScoreProfile {
            data: vec![Vec::new(); AMINO_ACID_COUNT],
            padding: padding.max(DEFAULT_PADDING),
        }
    }

    pub fn length(&self) -> usize {
        self.data
            .first()
            .map(|row| row.len().saturating_sub(2 * self.padding))
            .unwrap_or(0)
    }

    pub fn get(&self, letter: Letter, i: usize) -> &[Score] {
        let row = &self.data[letter as usize];
        &row[i + self.padding..]
    }

    pub fn pointers(&self, offset: usize) -> Vec<&[Score]> {
        let mut v = Vec::with_capacity(AMINO_ACID_COUNT);
        for letter in 0..AMINO_ACID_COUNT {
            v.push(self.get(letter as Letter, offset));
        }
        v
    }

    pub fn reverse(&self) -> Self {
        let mut r = self.clone();
        for row in &mut r.data {
            row.reverse();
        }
        r
    }
}

pub fn make_profile8(
    seq: &[Letter],
    cbs: Option<&[i8]>,
    padding: usize,
    matrix: &ScoreMatrix,
) -> LongScoreProfile<i8> {
    let mut profile = LongScoreProfile::new(padding);
    let len = seq.len() + 2 * profile.padding;
    for row in &mut profile.data {
        row.resize(len, 0);
    }
    for letter in 0..AMINO_ACID_COUNT {
        for (i, &subject) in seq.iter().enumerate() {
            let bias = cbs.and_then(|b| b.get(i)).copied().unwrap_or(0) as i32;
            let score = matrix.score(letter as Letter, subject) + bias;
            profile.data[letter][i + profile.padding] =
                score.clamp(i8::MIN as i32, i8::MAX as i32) as i8;
        }
    }
    profile
}

pub fn make_profile16(
    seq: &[Letter],
    cbs: Option<&[i8]>,
    padding: usize,
    matrix: &ScoreMatrix,
) -> LongScoreProfile<i16> {
    let mut profile = LongScoreProfile::new(padding);
    let len = seq.len() + 2 * profile.padding;
    for row in &mut profile.data {
        row.resize(len, 0);
    }
    for letter in 0..AMINO_ACID_COUNT {
        for (i, &subject) in seq.iter().enumerate() {
            let bias = cbs.and_then(|b| b.get(i)).copied().unwrap_or(0) as i32;
            let score = matrix.score(letter as Letter, subject) + bias;
            profile.data[letter][i + profile.padding] =
                score.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
        }
    }
    profile
}

pub fn make_profile16_target_matrix(
    seq: &[Letter],
    matrix: &TargetMatrix,
    padding: usize,
) -> LongScoreProfile<i16> {
    let mut profile = LongScoreProfile::new(padding);
    let len = seq.len() + 2 * profile.padding;
    for row in &mut profile.data {
        row.resize(len, 0);
    }
    for letter in 0..AMINO_ACID_COUNT {
        let row = &matrix.scores[letter << 5..];
        for (i, &subject) in seq.iter().enumerate() {
            profile.data[letter][i + profile.padding] = row[subject as usize] as i16;
        }
    }
    profile
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_long_score_profile_length_and_get() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let seq: Vec<Letter> = vec![0, 1, 2, 3];
        let profile = make_profile16(&seq, None, 4, &sm);
        assert_eq!(profile.padding, DEFAULT_PADDING);
        assert_eq!(profile.length(), seq.len());
        assert_eq!(profile.get(0, 0)[0], sm.score(0, seq[0]) as i16);
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
    }

    #[test]
    fn test_make_profile16_target_matrix() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let q = crate::stats::cbs::compute_composition(&[0, 1, 2, 3, 4]);
        let t = crate::stats::cbs::compute_composition(&[0, 1, 1, 2, 3]);
        let matrix = TargetMatrix::from_hauser_global(&q, &t, &sm);
        let seq: Vec<Letter> = vec![0, 1, 2];
        let profile = make_profile16_target_matrix(&seq, &matrix, 4);
        assert_eq!(profile.length(), seq.len());
        assert_eq!(profile.get(0, 0)[0], matrix.scores[0] as i16);
    }
}
