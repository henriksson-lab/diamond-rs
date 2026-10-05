# SIMD translation closure

Audited against the production dispatch object list in `diamond/CMakeLists.txt`.
The Rust implementation has native x86/x86-64 and AArch64 paths plus scalar
fallbacks on other stable-Rust targets.

| Upstream production path | Rust implementation | Native tiers |
| --- | --- | --- |
| Main SWIPE score/trace | `src/dp/swipe/simd_score*.rs`, `simd_trace*.rs`, `simd_adjusted_narrow.rs` | AVX2 32xi8, 16xi16, 8xi32; SSE/NEON 16xi8, 8xi16 and trace; per-lane adjusted matrices and gap scales |
| Anchored SWIPE | `src/dp/anchored.rs` | AVX2 8xi32, matching the upstream AVX2-only dispatch |
| Three-frame SWIPE | `src/dp/swipe/banded_3frame_simd.rs` | AVX2 16xi16, SSE2/NEON 8xi16 score; scalar traceback as upstream |
| Score profiles | `src/dp/score_profile.rs` | AVX2 and SSSE3 |
| Diagonal scan | `src/dp/scan_diags.rs` | AVX2, SSE4.1, NEON |
| Ungapped windows | `src/dp/ungapped_simd.rs` | AVX2, SSE4.1, NEON, including SIMD input transposes |
| Tantan | `src/masking/tantan_simd.rs` | AVX2, SSE4.1/SSSE3, NEON |
| Hamming stage/fingerprints | `src/search/hamming.rs`, `hamming_filter.rs` | AVX2, SSE2, NEON |
| Reduced seed distance | `src/search/sse_dist.rs` | SSSE3/SSE2, NEON |
| Matrix adjustment | `src/stats/target_freq_simd.rs` | AVX2/SSE2 |
| Byte transposes | `src/util/simd.rs` | AVX2 32x32, SSE2/NEON 16x16 |
| Hash fingerprints / bit counts | `src/util/data_structures.rs` | SSE2 and NEON |

The SWIPE score-width cascade preserves per-lane saturation promotion from
i8 to i16 to i32. Adjusted and ordinary targets may share a batch; their
per-lane matrices, gap scales, and composition-bias semantics remain distinct.
Tests cover banded/full matrices, semi-global scoring, reverse inputs, adjusted
matrices, masked letters, composition bias, tie ordering, traceback operations,
and saturation.

## Deliberate exclusions

- The vendored build can create AVX-512 objects, but its production dispatch
  macros declare and select only SSE4.1, AVX2, and NEON. AVX-512 therefore has
  no reachable upstream dispatch case.
- SIMD in `tools/benchmark*.cpp` is benchmark-only. Commented shape/seed and
  score-profile blocks are not production code.
- ARMv7 scalar fallback cross-compiles. ARMv7 NEON cannot be expressed on
  stable Rust today because `std::arch::arm` NEON intrinsics and runtime NEON
  detection are unstable (`rust-lang/rust` issues 111800 and 111190).
  AArch64 NEON is fully compiled and checked.

## Validation baseline

- Release suite: 1,650 passed, 0 failed, 1 ignored microbenchmark; integration
  and doc tests passed. External large-data tests remain explicitly ignored
  when their optional fixtures are absent.
- AArch64: `cargo check --target aarch64-unknown-linux-gnu --offline --all-targets`.
- ARMv7 scalar fallback: `cargo zigbuild --target armv7-unknown-linux-gnueabihf --lib --offline`.
- Small realistic blastp fixture: byte-identical output, SHA-256
  `ed82209f28b8b6327673e00ac4c82cf4da74c72fb250a03b792afd755a4f26d8`.
- Isolated seven-run medians on the closure host: one thread C++ 0.66 s /
  20,592 KiB, Rust 0.94 s / 16,156 KiB; four threads C++ 0.32 s / 28,896 KiB,
  Rust 0.49 s / 24,376 KiB. Rust peak RSS is 21.5% lower at one thread and
  15.6% lower at four threads. This small fixture does not yet meet the speed
  target, so the benchmark records that gap rather than claiming closure on
  performance.

Use `scripts/compare_cpp_rust.sh` to reproduce speed, peak RSS, and byte-parity
measurements. The script alternates blastp execution order and reports medians.
