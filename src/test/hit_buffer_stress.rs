//! Rust counterpart of `diamond/src/test/hit_buffer_stress.cpp`.

use crate::search::hit_buffer::{HitBuffer, HitBufferMode};

const BIN_COUNT: usize = 32;
const QUERIES_PER_BIN: usize = 1_000;
const QUERY_COUNT: usize = BIN_COUNT * QUERIES_PER_BIN;
const HITS_PER_QUERY: usize = 11;
const TARGET_LEN: u64 = 50_000;

fn hit_fp(query: u32, subject: u64, seed_offset: u32, score: u16) -> u64 {
    (query as u64)
        .wrapping_mul(2_654_435_769)
        .wrapping_add(subject.wrapping_mul(2_246_822_519))
        .wrapping_add((seed_offset as u64).wrapping_mul(3_266_489_917))
        .wrapping_add((score as u64).wrapping_mul(668_265_261))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HitBufferTestResult {
    pub passed: bool,
    pub expected_hits: usize,
    pub actual_hits: usize,
    pub expected_checksum: u64,
    pub actual_checksum: u64,
}

fn run_single_mode_with_size(
    membuf_mode: bool,
    bin_count: usize,
    queries_per_bin: usize,
) -> Result<HitBufferTestResult, String> {
    let query_count = bin_count * queries_per_bin;
    let key_partition: Vec<u32> = (1..=bin_count)
        .map(|bin| (bin * queries_per_bin) as u32)
        .collect();
    let mode = if membuf_mode {
        HitBufferMode::Memory
    } else {
        HitBufferMode::Disk
    };
    let mut buffer = HitBuffer::with_limits(
        key_partition,
        std::env::temp_dir(),
        false,
        1,
        query_count as u32,
        TARGET_LEN + 1,
        mode,
    )?;

    // The upstream stress test has one writer per bin. The translated buffer
    // intentionally performs packet I/O synchronously, so threads prepare the
    // same independent per-bin streams and ownership is merged into Writer in
    // deterministic bin order.
    let mut generators = Vec::with_capacity(bin_count);
    for bin in 0..bin_count {
        generators.push(std::thread::spawn(move || {
            let mut rows = Vec::with_capacity(queries_per_bin * HITS_PER_QUERY);
            let mut checksum = 0_u64;
            let query_begin = bin * queries_per_bin;
            let query_end = query_begin + queries_per_bin;
            for query in query_begin..query_end {
                let seed_offset = (query % 64) as u32;
                for hit in 0..HITS_PER_QUERY {
                    let subject = ((query * HITS_PER_QUERY + hit) as u64 % TARGET_LEN) + 1;
                    let score = ((query * 7 + hit * 13) % 65_534 + 1) as u16;
                    rows.push((query as u32, subject, seed_offset, score));
                    checksum =
                        checksum.wrapping_add(hit_fp(query as u32, subject, seed_offset, score));
                }
            }
            (rows, checksum)
        }));
    }

    let mut expected_checksum = 0_u64;
    for generator in generators {
        let (rows, checksum) = generator
            .join()
            .map_err(|_| "hit-buffer writer generator panicked".to_string())?;
        expected_checksum = expected_checksum.wrapping_add(checksum);
        let mut writer = buffer.writer();
        let mut previous_query = None;
        for (query, subject, seed_offset, score) in rows {
            if previous_query != Some(query) {
                writer.new_query(query, seed_offset)?;
                previous_query = Some(query);
            }
            writer.write(query, subject, score, 0)?;
        }
    }

    buffer.try_finish_writing()?;
    buffer.alloc_buffer();
    let mut actual_hits = 0_usize;
    let mut actual_checksum = 0_u64;
    while buffer.load(usize::MAX) {
        let Some((hits, _, _)) = buffer.try_retrieve()? else {
            break;
        };
        for hit in hits {
            actual_checksum = actual_checksum.wrapping_add(hit_fp(
                hit.query,
                hit.subject,
                hit.seed_offset,
                hit.score,
            ));
        }
        actual_hits += hits.len();
    }
    buffer.free_buffer();

    let expected_hits = query_count * HITS_PER_QUERY;
    Ok(HitBufferTestResult {
        passed: actual_hits == expected_hits && actual_checksum == expected_checksum,
        expected_hits,
        actual_hits,
        expected_checksum,
        actual_checksum,
    })
}

fn run_single_mode(membuf_mode: bool) -> Result<HitBufferTestResult, String> {
    debug_assert_eq!(QUERY_COUNT, BIN_COUNT * QUERIES_PER_BIN);
    run_single_mode_with_size(membuf_mode, BIN_COUNT, QUERIES_PER_BIN)
}

pub fn run_hit_buffer_stress_test() -> i32 {
    println!("\nHitBuffer stress test");
    println!("=====================");
    println!(
        "Threads = {}",
        std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)
    );
    let mut failures = 0;
    for membuf in [true, false] {
        let mode = if membuf { "in-memory (membuf)" } else { "disk" };
        print!("  Mode: {mode} ... ");
        match run_single_mode(membuf) {
            Ok(result) if result.passed => {
                println!("PASSED ({} hits, checksum ok)", result.actual_hits)
            }
            Ok(result) => {
                println!(
                    "FAILED\n    expected hits={} actual={}\n    expected checksum={} actual={}",
                    result.expected_hits,
                    result.actual_hits,
                    result.expected_checksum,
                    result.actual_checksum
                );
                failures += 1;
            }
            Err(error) => {
                println!("EXCEPTION: {error}");
                failures += 1;
            }
        }
    }
    println!("  Result: {}/2 passed", 2 - failures);
    println!("=====================");
    failures
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_storage_modes_preserve_every_hit_and_checksum() {
        for memory in [true, false] {
            let result = run_single_mode_with_size(memory, 4, 100).unwrap();
            assert!(result.passed, "{result:?}");
            assert_eq!(result.expected_hits, 4 * 100 * HITS_PER_QUERY);
        }
    }

    #[test]
    fn fingerprint_wraps_like_cpp_uint64() {
        assert_eq!(hit_fp(0, 0, 0, 0), 0);
        assert_eq!(
            hit_fp(u32::MAX, u64::MAX, u32::MAX, u16::MAX),
            6_983_481_896_302_944_870
        );
    }
}
