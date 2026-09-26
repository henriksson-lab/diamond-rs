#![allow(non_snake_case)]

use std::sync::Mutex;

pub const SEED: i64 = 31415;
const R_OFF: usize = 12;
const STATE_LEN: usize = 33;

#[derive(Debug, Clone)]
struct RandomState {
    state: [i64; STATE_LEN],
    rJ: usize,
    rK: usize,
}

impl RandomState {
    const fn new() -> Self {
        Self {
            state: [
                0xd53f1852, 0xdfc78b83, 0x4f256096, 0x0e643df7, 0x82c359bf, 0xc7794dfa, 0xd5e9ffaa,
                0x2c8cb64a, 0x2f07b334, 0xad5a7eb5, 0x96dc0cde, 0x6fc24589, 0xa5853646, 0xe71576e2,
                0x0dae30df, 0xb09ce711, 0x5e56ef87, 0x4b4b0082, 0x6f4f340e, 0xc5bb17e8, 0xd788d765,
                0x67498087, 0x9d7aba26, 0x261351d4, 0x411ee7ea, 0x0393a263, 0x2c5a5835, 0xc115fcd8,
                0x25e9132c, 0xd0c6e906, 0xc2bc5b2d, 0x6c065c98, 0x6e37bd55,
            ],
            rJ: R_OFF,
            rK: STATE_LEN - 1,
        }
    }

    fn number(&mut self) -> i64 {
        let r = self.state[self.rK].wrapping_add(self.state[self.rJ]);
        self.state[self.rK] = r;

        // C post-decrements both pointers before applying its wrap tests.
        if self.rJ == 0 {
            self.rJ = STATE_LEN - 1;
            self.rK -= 1;
        } else if self.rK == 0 {
            self.rK = STATE_LEN - 1;
            self.rJ -= 1;
        } else {
            self.rJ -= 1;
            self.rK -= 1;
        }
        (r >> 1) & 0x7fff_ffff
    }
}

static RANDOM_STATE: Mutex<RandomState> = Mutex::new(RandomState::new());

#[cfg(test)]
pub(crate) static TEST_RANDOM_LOCK: Mutex<()> = Mutex::new(());

pub fn seed(x: i64) {
    let mut random = RANDOM_STATE.lock().unwrap();
    random.state[0] = x;

    for i in 1..STATE_LEN {
        random.state[i] = 1103515245i64
            .wrapping_mul(random.state[i - 1])
            .wrapping_add(12345);
    }

    random.rJ = R_OFF;
    random.rK = STATE_LEN - 1;

    for _ in 0..(10 * STATE_LEN) {
        random.number();
    }
}

pub fn number() -> i64 {
    let mut random = RANDOM_STATE.lock().unwrap();
    random.number()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seeded_sequence_is_repeatable() {
        let _guard = TEST_RANDOM_LOCK.lock().unwrap();
        seed(SEED);
        let first: Vec<i64> = (0..8).map(|_| number()).collect();
        seed(SEED);
        let second: Vec<i64> = (0..8).map(|_| number()).collect();
        assert_eq!(first, second);
        assert_eq!(
            first,
            vec![
                1253310707, 1103551984, 1095634718, 2007891671, 930113008, 2138746150, 1613398876,
                128803596,
            ]
        );
    }

    #[test]
    fn test_default_state_number_range() {
        let _guard = TEST_RANDOM_LOCK.lock().unwrap();
        let value = number();
        assert!((0..=0x7fffffff).contains(&value));
    }

    #[test]
    fn negative_and_extreme_seeds_retain_signed_long_wrapping() {
        let _guard = TEST_RANDOM_LOCK.lock().unwrap();
        seed(-1);
        let negative: Vec<_> = (0..8).map(|_| number()).collect();
        assert_eq!(
            negative,
            [
                463_757_167,
                1_727_264_988,
                105_497_114,
                1_446_351_531,
                938_629_108,
                1_746_641_206,
                2_062_592_476,
                1_337_070_476,
            ]
        );
        seed(i64::MAX);
        assert_eq!((0..8).map(|_| number()).collect::<Vec<_>>(), negative);
    }
}
