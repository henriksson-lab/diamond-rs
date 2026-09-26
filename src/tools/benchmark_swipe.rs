//! Portable counterpart of the architecture-specific SWIPE cell benchmark.

pub const ROW_LEN: usize = 256;
pub const PROFILE_ROWS: usize = 32;
pub const CHANNELS: usize = 16;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Cell {
    pub score: [i8; CHANNELS],
    pub identities: [u16; CHANNELS],
    pub len: [u16; CHANNELS],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BenchmarkResult {
    pub cell_updates: u64,
    pub checksum: u64,
}

/// One row of the local-alignment recurrence exercised by the original SIMD
/// microbenchmark. Rust uses explicit initialized gap state instead of the
/// source benchmark's default-constructed vector temporaries.
pub fn update_row(
    query: &[u8; ROW_LEN],
    diagonal_cell: &mut [Cell; ROW_LEN],
    horizontal_gap: &mut [[i8; CHANNELS]; ROW_LEN],
    profile: &[[i8; CHANNELS]; PROFILE_ROWS],
) {
    let mut vertical_gap = [i8::MIN; CHANNELS];
    for i in 0..ROW_LEN {
        let old_diagonal = diagonal_cell[i];
        let mut next = Cell::default();
        for lane in 0..CHANNELS {
            let substitution =
                old_diagonal.score[lane].saturating_add(profile[query[i] as usize][lane]);
            let score = substitution
                .max(horizontal_gap[i][lane])
                .max(vertical_gap[lane])
                .max(0);
            next.score[lane] = score;
            next.len[lane] = old_diagonal.len[lane].saturating_add(1);
            next.identities[lane] = old_diagonal.identities[lane]
                .saturating_add((profile[query[i] as usize][lane] > 0) as u16);

            // The benchmark only measures update throughput. Fixed penalties
            // make the otherwise uninitialized C++ temporaries deterministic.
            horizontal_gap[i][lane] = horizontal_gap[i][lane]
                .saturating_sub(1)
                .max(score.saturating_sub(5));
            vertical_gap[lane] = vertical_gap[lane]
                .saturating_sub(1)
                .max(score.saturating_sub(5));
        }
        diagonal_cell[i] = next;
    }
}

fn next_random(state: &mut u64) -> u32 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1);
    (*state >> 32) as u32
}

/// Run the SWIPE microbenchmark workload and return an observable checksum so
/// the optimizer cannot discard it. The source's one-million iteration count
/// is supplied by callers rather than hard-wired into tests.
pub fn swipe_cell_update(iterations: usize, seed: u64) -> BenchmarkResult {
    let mut random = seed;
    let mut query = [0u8; ROW_LEN];
    for value in &mut query {
        *value = (next_random(&mut random) % PROFILE_ROWS as u32) as u8;
    }
    let mut profile = [[0i8; CHANNELS]; PROFILE_ROWS];
    for row in &mut profile {
        for value in row {
            *value = (next_random(&mut random) % 20) as i8 - 10;
        }
    }
    let mut diagonal = [Cell::default(); ROW_LEN];
    let mut horizontal = [[i8::MIN; CHANNELS]; ROW_LEN];
    for _ in 0..iterations {
        update_row(&query, &mut diagonal, &mut horizontal, &profile);
    }
    let first = diagonal[0];
    let checksum = first
        .score
        .iter()
        .map(|&value| value as i64 as u64)
        .chain(first.identities.iter().map(|&value| value as u64))
        .chain(first.len.iter().map(|&value| value as u64))
        .fold(0u64, |sum, value| sum.wrapping_add(value));
    BenchmarkResult {
        cell_updates: iterations as u64 * ROW_LEN as u64 * CHANNELS as u64,
        checksum,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_update_is_initialized_and_saturating() {
        let query = [0; ROW_LEN];
        let mut diagonal = [Cell::default(); ROW_LEN];
        let mut gaps = [[i8::MIN; CHANNELS]; ROW_LEN];
        let mut profile = [[0; CHANNELS]; PROFILE_ROWS];
        profile[0].fill(4);
        update_row(&query, &mut diagonal, &mut gaps, &profile);
        assert_eq!(diagonal[0].score, [4; CHANNELS]);
        assert_eq!(diagonal[0].identities, [1; CHANNELS]);
        assert_eq!(diagonal[0].len, [1; CHANNELS]);
    }

    #[test]
    fn benchmark_has_exact_work_count_and_deterministic_checksum() {
        let first = swipe_cell_update(8, 7);
        let second = swipe_cell_update(8, 7);
        assert_eq!(first, second);
        assert_eq!(first.cell_updates, 8 * 256 * 16);
        assert_ne!(first.checksum, 0);
    }
}
