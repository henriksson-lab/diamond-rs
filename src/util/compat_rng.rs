//! Deterministic process-global pseudo-random generator.
//!
//! Upstream DIAMOND calls C `rand`/`srand`, whose algorithm and `RAND_MAX`
//! differ between C runtimes. This module specifies the algorithm in Rust
//! instead: the 31-word additive generator used by GNU libc, including its
//! seeding and warm-up procedure. It preserves our established GNU/Linux
//! parity vectors while producing exactly the same sequence on every target.
//!
//! The state remains process-global and serialized, matching the observable
//! upstream API. Callers that reseed or draw concurrently still impose a
//! shared ordering and should not expect per-thread streams.

use std::sync::{Mutex, OnceLock};

pub const RAND_MAX: i32 = i32::MAX;

const DEGREE: usize = 31;
const SEPARATION: usize = 3;

struct DeterministicRand {
    state: [u32; DEGREE],
    front: usize,
    rear: usize,
}

impl DeterministicRand {
    fn seeded(seed: u32) -> Self {
        let mut rng = Self {
            state: [0; DEGREE],
            front: SEPARATION,
            rear: 0,
        };
        rng.reseed(seed);
        rng
    }

    fn reseed(&mut self, seed: u32) {
        // glibc stores the unsigned API seed in its signed 32-bit state before
        // running the Schrage form of the Park-Miller recurrence. Spell out
        // that conversion so seeds above i32::MAX are deterministic too.
        let mut word = i64::from((if seed == 0 { 1 } else { seed }) as i32);
        self.state[0] = word as u32;
        for slot in &mut self.state[1..] {
            let high = word / 127_773;
            let low = word % 127_773;
            word = 16_807 * low - 2_836 * high;
            if word < 0 {
                word += 2_147_483_647;
            }
            *slot = word as u32;
        }
        self.front = SEPARATION;
        self.rear = 0;
        for _ in 0..(DEGREE * 10) {
            let _ = self.next_i31();
        }
    }

    #[inline]
    fn next_i31(&mut self) -> i32 {
        let value = self.state[self.front].wrapping_add(self.state[self.rear]);
        self.state[self.front] = value;
        self.front += 1;
        if self.front == DEGREE {
            self.front = 0;
        }
        self.rear += 1;
        if self.rear == DEGREE {
            self.rear = 0;
        }
        (value >> 1) as i32
    }
}

fn global() -> &'static Mutex<DeterministicRand> {
    static RNG: OnceLock<Mutex<DeterministicRand>> = OnceLock::new();
    RNG.get_or_init(|| Mutex::new(DeterministicRand::seeded(1)))
}

#[inline]
pub fn next_i31() -> i32 {
    global()
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .next_i31()
}

#[inline]
pub fn seed(seed: u32) {
    global()
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .reseed(seed);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_vectors_are_fully_specified() {
        let mut rng = DeterministicRand::seeded(1);
        assert_eq!(
            std::array::from_fn::<_, 8, _>(|_| rng.next_i31()),
            [
                1_804_289_383,
                846_930_886,
                1_681_692_777,
                1_714_636_915,
                1_957_747_793,
                424_238_335,
                719_885_386,
                1_649_760_492,
            ]
        );
        rng.reseed(12_345);
        assert_eq!(
            std::array::from_fn::<_, 8, _>(|_| rng.next_i31()),
            [
                383_100_999,
                858_300_821,
                357_768_173,
                455_528_251,
                133_005_921,
                116_285_904,
                591_987_137,
                102_557_902,
            ]
        );
    }

    #[test]
    fn zero_seed_is_seed_one() {
        let mut zero_rng = DeterministicRand::seeded(0);
        let zero = std::array::from_fn::<_, 16, _>(|_| zero_rng.next_i31());
        let mut one_rng = DeterministicRand::seeded(1);
        let one = std::array::from_fn::<_, 16, _>(|_| one_rng.next_i31());
        assert_eq!(zero, one);
    }

    #[test]
    fn full_u32_seed_domain_is_defined() {
        let max_first = DeterministicRand::seeded(u32::MAX).next_i31();
        assert_eq!(DeterministicRand::seeded(u32::MAX).next_i31(), max_first);
        let high_bit_first = DeterministicRand::seeded(0x8000_0000).next_i31();
        assert_eq!(
            DeterministicRand::seeded(0x8000_0000).next_i31(),
            high_bit_first
        );
        assert_ne!(max_first, high_bit_first);
    }
}
