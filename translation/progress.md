# Translation progress

Last refreshed: 2026-09-26

The authoritative completion record is `file_map.csv`: each row represents an
original source file that has been manually audited, not merely name-matched.
All 198 C/C++ implementation files (189 `.cpp`, 9 `.c`) now have completed
audit rows: 180 translated, 17 parity-tested, and one replaced. This includes
all 140 core files and all 58 implementation files under the formerly excluded
`lib`, `contrib`, `test`, and `tools` areas.

## Current verification

- `cargo check --all-targets`: passed
- `cargo fmt --all --check`: passed
- Focused final-source tests: 18 passed
- `cargo test --release --lib -- --test-threads=1`: 1,599 passed
- `cargo test --release --tests -- --test-threads=1`: 1,618 passed,
  26 fixture-dependent real-data comparisons ignored
- `git diff --check`: passed
- Native `makeidx` seed index: byte-identical to the C++ executable on the
  repository test database
- `scripts/compare_cpp_rust.sh`: exercised end to end on the bundled
  389-sequence fixture; it correctly detects the current four-extra-hit parity
  gap while retaining a unified diff and time/RSS measurements

## CCC snapshot

| snapshot | Rust functions | raw missing | partial/stub candidates |
|---|---:|---:|---:|
| before source-file batches | 3,769 | 1,023 | 37 |
| wave 2 | 3,836 | 1,000 | 40 |
| wave 3 | 3,943 | 955 | 41 |
| wave 4 | 4,095 | 922 | 46 |
| wave 5 | 4,190 | 909 | 48 |
| wave 6 | 4,245 | 904 | 48 |
| wave 7 | 4,482 | 879 | 48 |
| wave 8 | 4,783 | 828 | 48 |
| wave 9 | 4,817 | 819 | 48 |
| wave 10 | 4,896 | 793 | 50 |
| wave 11 | 5,319 | 712 | 69 |
| final source-file wave | 5,937 | 534 | 82 |

The raw counters are deliberately not percentages. They count overloads,
constructors, destructors, templates, vendored C/C++ routines, and ambiguous
common names independently, while Rust often represents those with traits,
generic functions, ownership, or suffixed overload names. Conversely, a name
match does not prove semantic parity. Completion therefore requires a
source-file audit and focused tests recorded in `file_map.csv`.

Current generated reports are `/tmp/diamond-rust-final.json` and
`/tmp/diamond-cpp-final.json`; refresh instructions are in `README.md`.

## Source files not yet complete

None. End-to-end integration and parity work remains; source-file audit
completion does not imply every translated workflow is selected by the CLI.

## Important completeness caveats

- The default Cargo build has no C++ FFI feature. Unmatched commands and
  unsupported option combinations fall through `src/main.rs` to `run_legacy`,
  which returns an error in that build. Native `blastp` still routes DAA output,
  SEG masking, `--max-hsps != 1`, composition modes above 1, and global ranking
  to that fallback; `blastx` additionally routes frameshifting and `--swipe`.
- `src/ffi/ffi.rs` now translates the FFI dispatcher as a facade over the
  audited Rust `run::main`; `src/ffi.rs::run_cpp` remains an optional C++
  conformance bridge when the `ffi` feature is enabled.
- SEG and motif masking are now implemented and wired. The SEG implementation
  produced range-for-range identical output to C++ on 610 deterministic cases.
- Translation-unit coverage is ahead of end-to-end integration. In particular,
  the native `cluster`/`linclust`/`deepclust` CLI currently calls the simplified
  greedy implementation in `commands/cluster_cmd.rs`, not the translated
  cascaded, external, or multinode workflows. `commands/blastp.rs` likewise
  describes its native path as a simplified pipeline even though it now reuses
  many of the lower-level translated stages.
- The 58 implementation files under upstream `lib`, `contrib`, `test`, and
  `tools` now have explicit translate/replace classifications in `file_map.csv`.
- `TODO.md` contains useful historical performance notes but is stale in at
  least one important respect: banded CBS-aware alignment is now implemented
  and used by `dp::swipe`; its old “port banded swipe” next step should not be
  treated as the current inventory.
- The enabled repeat-sequence C++ comparison permits an 11-bit-score delta, and
  the larger 25/100/250/500-query dynamic comparisons are ignored unless the
  `/tmp/bench_*` fixtures are installed. The green default suite therefore does
  not establish full result parity across real workloads.
