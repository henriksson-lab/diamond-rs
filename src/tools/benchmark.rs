//! Portable benchmark driver translated from `diamond/src/tools/benchmark.cpp`.
//!
//! Architecture-specific kernels are explicit dependencies. This keeps the
//! source workloads executable without process globals or compile-time SIMD
//! dispatch and makes their work counts testable.

use super::benchmark_swipe::{swipe_cell_update, BenchmarkResult as CellBenchmarkResult};
use std::collections::BTreeMap;

pub const S1: &[u8] = b"mpeeeysefkelilqkelhvvyalshvcgqdrtllasillriflhekleslllctlndreismedeattlfrattlastlmeqymkatatqfvhhalkdsilkimeskqscelspskleknedvntnlthllnilselvekifmaseilpptlryiygclqksvqhkwptnttmrtrvvsgfvflrlicpailnprmfniisdspspiaartlilvaksvqnlanlvefgakepymegvnpfiksnkhrmimfldelgnvpelpdttehsrtdlsrdlaalheicvahsdelrtlsnergaqqhvlkkllaitellqqkqnqyt";
pub const S2: &[u8] = b"erlvelvtmmgdqgelpiamalanvvpcsqwdelarvlvtlfdsrhllyqllwnmfskeveladsmqtlfrgnslaskimtfcfkvygatylqklldpllrivitssdwqhvsfevdptrlepsesleenqrnllqmtekffhaiissssefppqlrsvchclyqvvsqrfpqnsigavgsamflrfinpaivspyeagildkkpppiierglklmskilqsianhvlftkeehmrpfndfvksnfdaarrffldiasdcptsdavnhslsfisdgnvlalhrllwnnqekigqylssnrdhkavgrrpfdkmatllaylgppe";
pub const S3: &[u8] = b"ttfgrcavksnqagggtrshdwwpcqlrldvlrqfqpsqnplggdfdyaeafqsldyeavkkdiaalmtesqdwwpadfgnygglfvrmawhsagtyramdgrggggmgqqrfaplnswpdnqnldkarrliwpikqkygnkiswadlmlltgnvalenmgfktlgfgggradtwqsdeavywgaettfvpqgndvrynnsvdinaradklekplaathmgliyvnpegpngtpdpaasakdireafgrmgmndtetvaliagghafgkthgavkgsnigpapeaadlgmqglgwhnsvgdgngpnqmtsgleviwtktptkwsngyleslinnnwtlvespagahqweavngtvdypdpfdktkfrkatmltsdlalindpeylkisqrwlehpeeladafakawfkllhrdlgpttrylgpevp";
pub const S4: &[u8] = b"lvhvasvekgrsyedfqkvynaialklreddeydnyigygpvlvrlawhisgtwdkhdntggsyggtyrfkkefndpsnaglqngfkflepihkefpwissgdlfslggvtavqemqgpkipwrcgrvdtpedttpdngrlpdadkdagyvrtffqrlnmndrevvalmgahalgkthlknsgyegpggaannvftnefylnllnedwklekndanneqwdsksgymmlptdysliqdpkylsivkeyandqdkffkdfskafekllengitfpkdapspfifktleeqgl";

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SwipeMode {
    Int8,
    Int16,
    Statistics,
    MatrixAdjust,
    Cbs,
    Traceback,
    BandedCbs,
    Banded,
    BandedTraceback,
    Anchored,
    AnchoredPipeline,
}

pub trait BenchmarkKernels {
    fn ungapped_window(&mut self, query: &[u8], target: &[u8], len: usize) -> i32;
    fn score_shuffle(&mut self, letter: usize, sequence: &[u8]) -> u64;
    fn run_swipe(&mut self, query: &[u8], target: &[u8], mode: SwipeMode) -> u64;
    fn diagonal_scores(&mut self, query: &[u8], target: &[u8]) -> [i32; 128];
    fn evalue_norm(&mut self, score: i32, query_len: usize) -> f64;
    fn calculate_evalue(&mut self, score: i32, query_len: usize, subject_len: usize) -> f64;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkloadResult {
    pub operations: u64,
    pub checksum: u64,
}

impl WorkloadResult {
    const fn new(operations: u64, checksum: u64) -> Self {
        Self {
            operations,
            checksum,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BenchmarkConfig {
    pub iterations: usize,
    pub threads: usize,
    pub simd_channels: usize,
}

impl Default for BenchmarkConfig {
    fn default() -> Self {
        Self {
            iterations: 1,
            threads: 1,
            simd_channels: 16,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BenchmarkReport {
    pub workloads: BTreeMap<&'static str, WorkloadResult>,
    pub cell_update: Option<CellBenchmarkResult>,
}

/// Exact scalar transpose including the source's left zero-padding when
/// fewer than `WIDTH` input rows are supplied.
pub fn transpose_scalar<const WIDTH: usize>(data: &[&[i8]], n: usize, output: &mut [i8]) {
    assert!(n <= WIDTH && data.len() >= n && output.len() >= WIDTH * WIDTH);
    let padding = WIDTH - n;
    for x in 0..padding {
        for y in 0..WIDTH {
            output[y * WIDTH + x] = 0;
        }
    }
    for x in padding..WIDTH {
        let input = data[x + n - WIDTH];
        assert!(input.len() >= WIDTH);
        for y in 0..WIDTH {
            output[y * WIDTH + x] = input[y];
        }
    }
}

/// The active source function is empty; its former stress body is commented.
pub fn hit_buffer() {}

pub fn benchmark_ungapped<K: BenchmarkKernels>(
    kernels: &mut K,
    query: &[u8],
    target: &[u8],
    iterations: usize,
) -> WorkloadResult {
    assert!(query.len() >= 64 && target.len() >= 64);
    let mut checksum = 0u64;
    for _ in 0..iterations {
        checksum = checksum.wrapping_add(kernels.ungapped_window(query, target, 64) as i64 as u64);
    }
    WorkloadResult::new(iterations as u64 * 64, checksum)
}

pub fn benchmark_ssse3_shuffle<K: BenchmarkKernels>(
    kernels: &mut K,
    sequence: &[u8],
    iterations: usize,
    channels: usize,
) -> WorkloadResult {
    let mut checksum = 0u64;
    for iteration in 0..iterations {
        checksum = checksum.wrapping_add(kernels.score_shuffle(iteration & 15, sequence));
    }
    WorkloadResult::new(iterations as u64 * channels as u64, checksum)
}

/// Both SIMD calls in the vendored source are commented out, so this records
/// the measured loop work without inventing a kernel invocation.
pub fn benchmark_ungapped_sse(iterations: usize, lanes: usize) -> WorkloadResult {
    WorkloadResult::new(iterations as u64 * lanes as u64 * 64, 0)
}

pub fn benchmark_transpose(iterations: usize, width: usize) -> WorkloadResult {
    assert!(matches!(width, 16 | 32));
    let mut input = vec![0i8; width * width];
    for (index, value) in input.iter_mut().enumerate() {
        *value = index as i8;
    }
    let mut output = vec![0i8; width * width];
    for _ in 0..iterations {
        let rows: Vec<_> = input.chunks_exact(width).collect();
        if width == 16 {
            transpose_scalar::<16>(&rows, 16, &mut output);
        } else {
            transpose_scalar::<32>(&rows, 32, &mut output);
        }
        input[0] = output[0];
    }
    let checksum = output
        .iter()
        .fold(0u64, |sum, &value| sum.wrapping_add(value as i64 as u64));
    WorkloadResult::new(iterations as u64 * width as u64 * width as u64, checksum)
}

fn repeat_swipe<K: BenchmarkKernels>(
    kernels: &mut K,
    query: &[u8],
    target: &[u8],
    mode: SwipeMode,
    iterations: usize,
    lanes: usize,
    cells_per_lane: usize,
) -> WorkloadResult {
    let mut checksum = 0u64;
    for _ in 0..iterations {
        checksum = checksum.wrapping_add(kernels.run_swipe(query, target, mode));
    }
    WorkloadResult::new(
        iterations as u64 * cells_per_lane as u64 * lanes as u64,
        checksum,
    )
}

pub fn mt_swipe<K: BenchmarkKernels>(
    kernels: &mut K,
    query: &[u8],
    target: &[u8],
    config: BenchmarkConfig,
) -> WorkloadResult {
    let query_len = query.len().min(255);
    // The injected adapter is called sequentially; its explicit thread count
    // retains the source workload without sharing `&mut K` unsafely.
    repeat_swipe(
        kernels,
        &query[..query_len],
        target,
        SwipeMode::Int8,
        config.iterations * config.threads,
        config.simd_channels,
        query_len * target.len(),
    )
}

pub fn swipe<K: BenchmarkKernels>(
    kernels: &mut K,
    query: &[u8],
    target: &[u8],
    config: BenchmarkConfig,
) -> BTreeMap<SwipeMode, WorkloadResult> {
    let query = &query[..query.len().min(255)];
    [
        SwipeMode::Int8,
        SwipeMode::Int16,
        SwipeMode::Statistics,
        SwipeMode::MatrixAdjust,
        SwipeMode::Cbs,
        SwipeMode::Traceback,
    ]
    .into_iter()
    .map(|mode| {
        (
            mode,
            repeat_swipe(
                kernels,
                query,
                target,
                mode,
                config.iterations,
                config.simd_channels,
                query.len() * target.len(),
            ),
        )
    })
    .collect()
}

pub fn banded_swipe<K: BenchmarkKernels>(
    kernels: &mut K,
    query: &[u8],
    target: &[u8],
    iterations: usize,
) -> BTreeMap<SwipeMode, WorkloadResult> {
    [
        SwipeMode::BandedCbs,
        SwipeMode::Banded,
        SwipeMode::BandedTraceback,
    ]
    .into_iter()
    .map(|mode| {
        (
            mode,
            repeat_swipe(
                kernels,
                query,
                target,
                mode,
                iterations,
                16,
                query.len() * 65,
            ),
        )
    })
    .collect()
}

pub fn anchored_swipe<K: BenchmarkKernels>(
    kernels: &mut K,
    query: &[u8],
    target: &[u8],
    iterations: usize,
) -> BTreeMap<SwipeMode, WorkloadResult> {
    let query = &query[..query.len().min(128)];
    let target = &target[..target.len().min(128)];
    [SwipeMode::Anchored, SwipeMode::AnchoredPipeline]
        .into_iter()
        .map(|mode| {
            (
                mode,
                repeat_swipe(kernels, query, target, mode, iterations, 16, 128 * 64),
            )
        })
        .collect()
}

pub fn diag_scores<K: BenchmarkKernels>(
    kernels: &mut K,
    query: &[u8],
    target: &[u8],
    iterations: usize,
) -> WorkloadResult {
    let mut checksum = 0u64;
    for iteration in 0..iterations {
        let scores = kernels.diagonal_scores(query, target);
        checksum = checksum.wrapping_add(scores[iteration & 127] as i64 as u64);
    }
    WorkloadResult::new(iterations as u64 * target.len() as u64 * 128, checksum)
}

pub fn evalue<K: BenchmarkKernels>(kernels: &mut K, iterations: usize) -> [WorkloadResult; 2] {
    let mut normalized = 0.0f64;
    let mut alp = 0.0f64;
    for score in 0..iterations {
        normalized += kernels.evalue_norm(score as i32, 300);
        alp += kernels.calculate_evalue(300, 300, 300);
    }
    [
        WorkloadResult::new(iterations as u64, normalized.to_bits()),
        WorkloadResult::new(iterations as u64, alp.to_bits()),
    ]
}

/// The active source benchmark only clears this matrix; both optimizer calls
/// are commented out.
pub fn matrix_adjust(_query: &[u8], _target: &[u8], iterations: usize) -> Vec<f64> {
    let mut matrix = vec![0.0; 20 * 20];
    for _ in 0..iterations {
        matrix.fill(0.0);
    }
    matrix
}

pub fn benchmark<K: BenchmarkKernels>(
    kind: Option<&str>,
    kernels: &mut K,
    config: BenchmarkConfig,
) -> BenchmarkReport {
    if kind == Some("swipe") {
        return BenchmarkReport {
            workloads: BTreeMap::new(),
            cell_update: Some(swipe_cell_update(config.iterations, 1)),
        };
    }

    let mut workloads = BTreeMap::new();
    let _ = matrix_adjust(S1, S2, config.iterations);
    for (mode, result) in anchored_swipe(kernels, S1, S2, config.iterations) {
        workloads.insert(mode_name(mode), result);
    }
    for (mode, result) in swipe(kernels, S3, S4, config) {
        workloads.insert(mode_name(mode), result);
    }
    workloads.insert(
        "diagonal_scores",
        diag_scores(kernels, S1, S2, config.iterations),
    );
    for (mode, result) in banded_swipe(kernels, S1, S2, config.iterations) {
        workloads.insert(mode_name(mode), result);
    }
    let [normal, alp] = evalue(kernels, config.iterations);
    workloads.insert("evalue", normal);
    workloads.insert("evalue_alp", alp);
    workloads.insert(
        "ungapped",
        benchmark_ungapped(kernels, &S1[34..], &S2[33..], config.iterations),
    );
    workloads.insert(
        "shuffle",
        benchmark_ssse3_shuffle(kernels, S1, config.iterations, config.simd_channels),
    );
    workloads.insert(
        "ungapped_sse",
        benchmark_ungapped_sse(config.iterations, 16),
    );
    workloads.insert("transpose16", benchmark_transpose(config.iterations, 16));
    workloads.insert("transpose32", benchmark_transpose(config.iterations, 32));
    BenchmarkReport {
        workloads,
        cell_update: None,
    }
}

const fn mode_name(mode: SwipeMode) -> &'static str {
    match mode {
        SwipeMode::Int8 => "swipe_int8",
        SwipeMode::Int16 => "swipe_int16",
        SwipeMode::Statistics => "swipe_statistics",
        SwipeMode::MatrixAdjust => "swipe_matrix_adjust",
        SwipeMode::Cbs => "swipe_cbs",
        SwipeMode::Traceback => "swipe_traceback",
        SwipeMode::BandedCbs => "banded_cbs",
        SwipeMode::Banded => "banded",
        SwipeMode::BandedTraceback => "banded_traceback",
        SwipeMode::Anchored => "anchored",
        SwipeMode::AnchoredPipeline => "anchored_pipeline",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct Kernels {
        calls: usize,
    }

    impl BenchmarkKernels for Kernels {
        fn ungapped_window(&mut self, query: &[u8], target: &[u8], len: usize) -> i32 {
            self.calls += 1;
            query[..len]
                .iter()
                .zip(&target[..len])
                .map(|(a, b)| i32::from(a == b))
                .sum()
        }

        fn score_shuffle(&mut self, letter: usize, sequence: &[u8]) -> u64 {
            self.calls += 1;
            u64::from(sequence[letter])
        }

        fn run_swipe(&mut self, query: &[u8], target: &[u8], mode: SwipeMode) -> u64 {
            self.calls += 1;
            query.len() as u64 + target.len() as u64 + mode as u64
        }

        fn diagonal_scores(&mut self, query: &[u8], target: &[u8]) -> [i32; 128] {
            self.calls += 1;
            [query.len() as i32 - target.len() as i32; 128]
        }

        fn evalue_norm(&mut self, score: i32, query_len: usize) -> f64 {
            self.calls += 1;
            score as f64 / query_len as f64
        }

        fn calculate_evalue(&mut self, score: i32, query_len: usize, subject_len: usize) -> f64 {
            self.calls += 1;
            score as f64 / (query_len * subject_len) as f64
        }
    }

    #[test]
    fn scalar_transpose_pads_and_transposes() {
        let rows = [[1i8, 2, 3, 4], [5, 6, 7, 8]];
        let refs = rows.iter().map(|row| row.as_slice()).collect::<Vec<_>>();
        let mut output = [99i8; 16];
        transpose_scalar::<4>(&refs, 2, &mut output);
        assert_eq!(output, [0, 0, 1, 5, 0, 0, 2, 6, 0, 0, 3, 7, 0, 0, 4, 8]);
    }

    #[test]
    fn individual_workloads_have_exact_source_cell_counts() {
        let mut kernels = Kernels::default();
        let ungapped = benchmark_ungapped(&mut kernels, S1, S2, 3);
        assert_eq!(ungapped.operations, 3 * 64);
        let banded = banded_swipe(&mut kernels, S1, S2, 2);
        assert_eq!(
            banded[&SwipeMode::Banded].operations,
            2 * S1.len() as u64 * 65 * 16
        );
        let diagonal = diag_scores(&mut kernels, S1, S2, 2);
        assert_eq!(diagonal.operations, 2 * S2.len() as u64 * 128);
    }

    #[test]
    fn driver_covers_active_source_benchmarks_and_swipe_shortcut() {
        let mut kernels = Kernels::default();
        let report = benchmark(None, &mut kernels, BenchmarkConfig::default());
        assert!(report.workloads.contains_key("swipe_int8"));
        assert!(report.workloads.contains_key("anchored_pipeline"));
        assert!(report.workloads.contains_key("transpose32"));
        assert!(report.cell_update.is_none());

        let shortcut = benchmark(Some("swipe"), &mut kernels, BenchmarkConfig::default());
        assert!(shortcut.workloads.is_empty());
        assert_eq!(shortcut.cell_update.unwrap().cell_updates, 256 * 16);
    }

    #[test]
    fn commented_source_paths_remain_observable_noops() {
        hit_buffer();
        assert_eq!(
            benchmark_ungapped_sse(5, 32),
            WorkloadResult::new(5 * 32 * 64, 0)
        );
        assert_eq!(matrix_adjust(S1, S2, 3), vec![0.0; 400]);
    }
}
