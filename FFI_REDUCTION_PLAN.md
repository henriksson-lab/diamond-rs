# Native ABI and libc reduction plan

Status: implementation complete; external macOS/Windows runtime gates open  
Started: 2026-10-06

## Goal

Reduce handwritten platform ABI declarations and unsafe descriptor ownership while preserving exact DIAMOND behavior, output parity, speed, and RSS.

The practical target is **no handwritten operating-system/libc ABI declarations in application modules**, not “no operating-system ABI calls.” Rust's standard library and portable system crates ultimately use those APIs too. Native SQLite, pure-Rust zstd, and the optional C++ conformance oracle are tracked separately below.

## Non-negotiable constraints

- [x] Preserve byte-for-byte output parity unless a format's documented contract is decoded equivalence rather than identical compressed bytes.
- [x] Establish the GNU raw/downstream reference vectors, then intentionally
  standardize that sequence across targets instead of preserving divergent C
  runtime generators.
- [x] Do not remove GNU `malloc_trim` until realistic RSS measurements prove it unnecessary.
- [x] Keep the optional upstream C++ FFI disabled by default and out of production benchmarks.
- [x] Make ownership explicit with RAII; do not exchange one set of raw descriptors or pointers for another.
- [x] Benchmark natural, nonduplicated data with repeated alternating runs. Check one and multiple threads and the relevant `16G`, `0G`, and adaptive-spill memory modes.
- [x] Treat every observable behavior change as a separate, documented decision. This project is pursuing a faithful translation, not an approximation.

## Baseline inventory

At the start of this plan, `src/` contains 17 handwritten `unsafe extern` blocks. The important groups are:

- `src/util/sequence.rs`: libc `memchr` in the stage-2 `clip` hot path.
- `src/util/io.rs`, `util/io/temp_file.rs`, `input_file.rs`, `output_file.rs`, and `file_sink.rs`: `dup`, `close`, raw descriptor ownership, and `FromRawFd`/`IntoRawFd` transitions.
- `src/util/system/system.rs`: Unix `open`, `mmap`, `munmap`, `close`, Linux `sysconf`, and Windows console APIs.
- `src/util/system/get_rss.rs`: Unix `getrusage`, Linux page size, macOS Mach task information, and Windows process-memory APIs.
- `src/util/io/output_file.rs`: `posix_fadvise` with handwritten constants.
- `src/stats/*` and `src/lib/tantan/*`: deterministic Rust RNG in `compat_rng`;
  no C RNG calls remain.
- `src/commands/blastp.rs`: GNU `malloc_trim`, currently used to reduce RSS after large releases/spills.
- `src/data/blastdb/blastdb.rs`: raw `libsqlite3-sys` handles and statements.
- `src/ffi.rs` and `build.rs`: optional upstream DIAMOND C++ conformance adapter (`ffi` feature, default off).
- `zstd`: safe Rust-facing API backed by native `zstd-sys`.

### Baseline tasks

- [x] Record `rg` inventory of all `extern "C"`, `extern "system"`, `link_name`, and raw-FD conversions.
- [x] Record default and `--features ffi` dependency trees separately.
- [x] Capture current full test results on Linux.
- [x] Capture representative search/clustering parity, wall time, and peak RSS before rewrites.
- [x] Add a CI check that rejects new handwritten ABI blocks outside a small allowlist.

## Rewrite checkpoint evidence (historical)

The following focused checks and full all-targets run passed on Linux x86-64
GNU during the rewrite. Some source-audit and adjacent-upstream integration
suites named below were subsequently removed so routine tests are
self-contained. The current self-contained all-target suite passes 1,697
library tests with 6 ignored, together with its binary and integration tests.
These records do not replace the unchecked native cross-platform gates below.

Dependency inventory was captured with `cargo tree -e normal --no-default-features` and `cargo tree -e all --features ffi`. The normal runtime graph contains the maintained bindings/wrappers introduced here (`libc`, `memmap2`, `nix`, `rustix`, `rusqlite`) plus bundled native `libsqlite3-sys` and pure-Rust `zstd-pure-rs`. `zstd-sys` is no longer in the graph. `libc::memchr` is an evidence-backed hot-path exception explained in Phase 1. `cmake 0.1.58` appears only in the all-edge `ffi` graph as the build dependency for the optional C++ oracle; enabling `ffi` does not add another Rust runtime dependency.

- `cargo test --lib clip_matches_scalar_reference_for_edge_layouts`: 1 passed, 1,743 filtered out.
- `cargo test --lib util::io::temp_file::tests`: 6 passed and 1,736 filtered out, covering linked/unlinked ownership, shared clone cursors, ownership transfer, close behavior, and truncation semantics.
- `cargo test --lib util::system::`: 13 passed, 1,731 filtered out, including the isolated touched-memory RSS child test.
- `cargo test --lib data::blastdb::blastdb::tests`: 5 passed, 4 external/diagnostic tests ignored, and 1,735 filtered out.
- `cargo test --test native_abi_boundary`: 1 passed.
- `cargo test --test compat_rng_vectors`: 1 passed and 1 subprocess helper ignored; only GNU libc vectors are established.
- Final `cargo check --all-targets` passed on Linux x86-64 GNU; `cargo check --target x86_64-pc-windows-gnu --all-targets` also passed without warnings.
- Final `cargo test --all-targets -- --test-threads=1` on Linux x86-64 GNU: library tests 1,736 passed and 6 ignored; binary tests 7 passed; `equivalence` 8 passed; `pipeline` 6 passed; the remaining integration suites passed, with data-dependent diagnostics ignored where their external fixtures were unavailable.
- `cargo test --features ffi --test regression -- --test-threads=1`: 6 passed, including the upstream built-in 20/20 command matrix.

### Representative checkpoint-to-rewrite benchmark

The clean tracked checkpoint `c954916f8c48e945c8f69f025478053568396551` was extracted with `git archive` into `/tmp` and compared with the current FFI-rewritten worktree using `rustc 1.98.1`, identical release profiles, and the current lockfile's compatible shared dependency versions. The checkpoint's `open` declaration was already variadic, so no source patch was required. The two frozen benchmark binaries have SHA-256 hashes `294bc1cdffe171c822e1bfbe0f2255f67f819ede2efd5cc328d060edbbcb40ef` (checkpoint) and `31a188f7178546c1871476811eead5fd379f6cb2641f4ebc92dadf1573101f19` (rewrite).

Each cell below is the mean of four release-mode `/usr/bin/time` runs. Runs alternated checkpoint/rewrite order and reversed that order on even repeats. The quick search used 300 natural, nonduplicated AMR protein queries against a 4,000-sequence, 1,358,898-letter database. The adaptive case used 2,000 queries against the same database because the quick fixture does not accumulate enough hits for a genuine delayed spill. Every `64M` run logged the memory-limit transition and read back 21.1–21.3 MiB of compressed spill data. The clustering comparison used a natural 500-sequence, 179,068-letter AMR database. Other work was present on the shared host, so small timing differences should be interpreted with the repeated alternating design in mind.

| Workload | Memory | Threads | Checkpoint wall | Rewrite wall | Wall ratio | Checkpoint RSS | Rewrite RSS | RSS ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `blastp`, 300 queries | `16G` | 1 | 3.720 s | 3.825 s | 1.028 | 28,017 KiB | 27,958 KiB | 0.998 |
| `blastp`, 300 queries | `16G` | 4 | 1.098 s | 1.130 s | 1.030 | 30,915 KiB | 31,269 KiB | 1.011 |
| `blastp`, 300 queries | `0G` | 1 | 3.733 s | 3.878 s | 1.039 | 18,122 KiB | 17,916 KiB | 0.989 |
| `blastp`, 300 queries | `0G` | 4 | 1.170 s | 1.203 s | 1.028 | 25,133 KiB | 25,401 KiB | 1.011 |
| `blastp`, 2,000 queries | delayed `64M` spill | 1 | 20.875 s | 21.580 s | 1.034 | 69,398 KiB | 69,503 KiB | 1.002 |
| `blastp`, 2,000 queries | delayed `64M` spill | 4 | 5.943 s | 6.390 s | 1.075 | 67,815 KiB | 71,409 KiB | 1.053 |
| `cluster`, 500 sequences | n/a | 1 | 2.520 s | 2.488 s | 0.987 | 27,048 KiB | 26,972 KiB | 0.997 |
| `cluster`, 500 sequences | n/a | 4 | 1.125 s | 1.165 s | 1.036 | 39,442 KiB | 37,944 KiB | 0.962 |

All 32 quick-search outputs were byte-identical (SHA-256 `00a88bb1b3c48c69bcc11b352aa4c14e36104d4d2359b1653c426c126fdd88c6`), all 16 delayed-spill outputs were byte-identical (`d324e3efeaa48cff4885dfb4ad9c880a699b5f9a8c372b97f84629e2f25ccfc9`), and all 16 cluster outputs were byte-identical (`9fd51c2dd37fc0a946afab68508fced0b8ef246a68f9c1123aea21691216f86c`). Raw commands, timings, RSS records, output files, hashes, host metadata, and spill diagnostics are retained under `.tmp/ffi-prepost/`.

Parity passed for this representative gate, but the first rewrite was consistently 2.8–3.9% slower in the quick cells and 3.4–7.5% slower in the delayed-spill cells. The four-thread delayed-spill cell also used 5.3% more peak RSS. Profiling attributed the repeatable speed loss to the portable `memchr` replacement, as detailed in Phase 1. That implementation and binary were rejected rather than treating the regression as acceptable noise.

After restoring upstream's exact C-runtime `memchr` kernel through the maintained `libc` crate binding, the entire matrix was repeated from scratch with four balanced alternating runs per cell. The fixed binary has SHA-256 `93155e5cf506d559902784ab244cb00afde72aef2737f09bf611dc2db2ca77e6`.

| Workload | Memory | Threads | Checkpoint wall | Fixed rewrite wall | Wall ratio | Checkpoint RSS | Fixed rewrite RSS | RSS ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `blastp`, 300 queries | `16G` | 1 | 3.593 s | 3.590 s | 0.999 | 28,041 KiB | 27,910 KiB | 0.995 |
| `blastp`, 300 queries | `16G` | 4 | 1.033 s | 1.040 s | 1.007 | 31,187 KiB | 31,184 KiB | 1.000 |
| `blastp`, 300 queries | `0G` | 1 | 3.498 s | 3.528 s | 1.009 | 18,471 KiB | 17,899 KiB | 0.969 |
| `blastp`, 300 queries | `0G` | 4 | 1.123 s | 1.098 s | 0.978 | 26,205 KiB | 24,525 KiB | 0.936 |
| `blastp`, 2,000 queries | delayed `64M` spill | 1 | 19.683 s | 19.620 s | 0.997 | 69,360 KiB | 69,426 KiB | 1.001 |
| `blastp`, 2,000 queries | delayed `64M` spill | 4 | 5.578 s | 5.565 s | 0.998 | 68,305 KiB | 69,302 KiB | 1.015 |
| `cluster`, 500 sequences | n/a | 1 | 2.403 s | 2.380 s | 0.991 | 26,992 KiB | 27,150 KiB | 1.006 |
| `cluster`, 500 sequences | n/a | 4 | 1.138 s | 1.133 s | 0.996 | 38,681 KiB | 39,466 KiB | 1.020 |

All 64 fixed-matrix outputs remained byte-identical with the checkpoint and retained the three hashes above. Every adaptive run again performed a genuine delayed transition and read back 21.1–21.3 MiB of compressed spill data. Search wall time is now within 0.9% of the checkpoint except for a 2.2% improvement in the forced-disk four-thread cell; delayed-spill time is within 0.3%. RSS differences are within run-to-run variability except for a 3.1–6.4% improvement in forced-disk mode. The representative parity, speed, and RSS acceptance gate therefore passes after the correction. Raw fixed-run artifacts and the exact driver are retained under `.tmp/ffi-prepost/fixed-results/` and `.tmp/ffi-prepost/run_fixed.sh`.

## Parallel execution model

The work should proceed in waves. Within a wave, workers own disjoint modules. The integrator alone edits `Cargo.toml` and `Cargo.lock`, preventing dependency and lockfile conflicts.

### Wave 1 — independent low-risk work

- Worker A: `memchr` replacement and hot-path tests.
- Worker B: temp-file/File ownership rewrite across `util/io`.
- Worker C: platform RSS/console bindings and runtime tests.
- Integrator: add agreed dependencies, maintain the inventory, run the complete suite, and merge each lane separately.

The mmap rewrite may run beside the descriptor rewrite because it is isolated in `util/system/system.rs`. The `fadvise` rewrite must follow or belong to the descriptor worker because both touch `output_file.rs`.

### Wave 2 — independent safe wrappers

- Worker A: mmap RAII replacement.
- Worker B: typed `fadvise` and Linux cache discovery.
- Worker C: safe SQLite wrapper while retaining the SQLite engine.

### Wave 3 — parity-sensitive investigation

- One worker captures per-platform RNG reference vectors.
- One worker runs `malloc_trim` A/B RSS experiments.
- One worker audits zstd interoperability and the remaining C++ oracle routes.

Workers must not combine unrelated behavior changes. Each lane should produce a focused commit that can be reverted independently.

## Phase 1 — replace the handwritten `memchr` ABI

- [x] Add `libc` as a direct dependency for its maintained, target-specific `memchr` declaration.
- [x] Replace the handwritten declaration with `libc::memchr`, retaining upstream's exact C-runtime kernel and pointer arithmetic.
- [x] Add exhaustive comparison with a scalar reference for empty input, all anchors, delimiters before/at/after the anchor, and adjacent delimiters.
- [x] Benchmark delimiter-free, dense-delimiter, and near-anchor buffers in release mode.
- [x] Run full search output parity and realistic speed/RSS checks.

The initial implementation used the portable `memchr` crate. Although its reproducible diagnostic was much faster than a scalar loop, `perf` on the natural delayed-spill workload showed that this was not faithful enough for the actual target: the crate's AVX2 search/dispatch symbols consumed a combined 3.18% of samples, versus 1.36% in glibc `__memchr_evex` at the checkpoint. The intermediate rewrite used 23,869 ms task-clock versus 22,411 ms at the checkpoint and caused the 2.8–7.5% regression recorded above. This implementation was rejected.

The corrected implementation calls the same C-runtime kernel as upstream through the maintained `libc` crate declaration. The handwritten ABI is still gone, but the optimized target runtime remains selectable by glibc's IFUNC machinery. A new profile measured 1.39% in `__memchr_evex` and 22,431 ms task-clock, matching the checkpoint's 1.36% and 22,411 ms. It also retired slightly fewer instructions than the checkpoint in those profile runs (180.85 billion versus 181.01 billion). Profile data, reports, counter output, and exact result files are under `.tmp/ffi-prepost/profile/`.

The reproducible ignored diagnostic `benchmark_clip_memchr_layouts` compares the corrected `clip` with its scalar reference over 20,000 iterations of 64 KiB buffers. On this Linux x86-64 host in release mode, C-runtime `memchr`/scalar elapsed-time ratios were 0.040 for delimiter-free data, 0.634 for a delimiter every 17 bytes, and 0.014 for delimiters adjacent to the anchor. Checksums were equal in every layout. Command: `cargo test --release --lib benchmark_clip_memchr_layouts -- --ignored --nocapture`.

Acceptance: no handwritten `memchr` ABI remains and the stage-2 hot path has no material regression. Calling the maintained `libc` binding is an explicit evidence-backed performance exception to preferring a pure-Rust wrapper; it preserves the exact upstream kernel on every target C runtime without maintaining an ABI declaration locally.

## Phase 2 — make temp-file ownership safe

Preferred design: store an owning `std::fs::File` in `TempFileData`, because every in-tree consumer ultimately needs a `File`. Use `Option<File>` only where ownership must be taken.

- [x] Add characterization tests for linked and unlinked temp files before changing the API.
- [x] Test that clones share the file cursor, matching Unix `dup`/`File::try_clone` behavior.
- [x] Replace `TempFileData.fd: i32`, the `-1` sentinel, `IntoRawFd`, and its manual `Drop` with an owning `File`.
- [x] Add private `try_clone_file` and, if required, `take_file` operations.
- [x] Replace `InputFile::from_temp_file_data` and `OutputFile::from_temp_file_data` duplication with `File::try_clone`.
- [x] Replace unsafe `FileSink::from_fd` with safe `FileSink::from_file`.
- [x] Replace the exported Linux-numbered `posix_flags` helper with portable typed mode parsing over the exact accepted `fopen` spellings.
- [x] Replace the remaining `File::temporary` raw-FD conversion.
- [x] Remove the module-level `dup` declaration only after every consumer has migrated.
- [x] Preserve the current behavior that an already-open temporary file is not retruncated when wrapped.
- [x] Preserve existing error text where it is user-visible.

Acceptance: no raw descriptor ownership or manual `close` remains in the temp-file path; linked and unlinked files remain readable/writable through owning files and clones; ownership-transfer, shared-cursor, close, and truncation behavior tests pass. Descriptor counting is intentionally not part of this gate: no leak was observed, and platform-specific process-resource enumeration would add test-only complexity without diagnosing a known defect.

## Phase 3 — replace the raw mmap interface

The current `mmap_file -> (*mut u8, usize, i32)` and `unmap_file` API has no production caller; only its unit test uses it.

- [x] Decide whether to remove the unused public API or retain it as an owning `MappedFile`/`memmap2::Mmap` API.
- [x] Open through `std::fs::File::open`.
- [x] Contain the mapping's unavoidable unsafe construction in one documented function.
- [x] Return an owning mapping that dereferences to `[u8]`; remove manual unmap and close operations.
- [x] Preserve and test current zero-length-file behavior: empty files still return an error.
- [x] Intentionally extend the owning `memmap2::Mmap` API to Windows; this removes the former platform-specific unsupported branch without exposing raw platform handles.

Acceptance: `open`, `mmap`, `munmap`, and `close` declarations disappear from `system.rs`, and no raw mapping pointer escapes.

## Phase 4 — typed file advice and cache discovery

- [x] Replace handwritten `posix_fadvise` with `rustix::fs::fadvise` and `Advice::WillNeed`.
- [x] Continue ignoring the result initially, matching current observable behavior.
- [x] Document that the current expression `SEQUENTIAL | WILLNEED` numerically resolves to `WILLNEED`; issuing two calls would be an intentional behavior/performance change requiring its own benchmark.
- [x] Replace hardcoded Linux `_SC_LEVEL3_CACHE_SIZE = 194` with `/sys/devices/system/cpu/cpu*/cache/index*` discovery.
- [x] Define L3 semantics as the maximum shared unified/data L3 domain, not a sum across CPUs.
- [x] Add parser fixtures for malformed sizes, overflow, multiple CPUs, and multiple sockets.
- [x] Return zero on unsupported platforms until a validated native implementation exists.

Acceptance: no `sysconf` or `posix_fadvise` declarations remain.

## Phase 5 — portable process and console APIs

Current RSS is correctness-sensitive because search spilling consults it. Peak RSS is also part of our benchmark contract.

- [x] Keep Linux `/proc/self/statm`, but replace the page-size `sysconf` call with `rustix::param::page_size()`.
- [x] Replace the untyped `[c_long; 64]` `getrusage` scratch buffer with a maintained typed binding. Evaluate `nix::sys::resource`; if its dependency cost is unacceptable, isolate `libc::{rusage, getrusage}` behind one private module.
- [x] Replace handwritten macOS Mach layouts/constants with maintained `mach2` bindings, or another narrow maintained wrapper selected during implementation.
- [x] Replace Windows process memory declarations with `windows-sys` and `K32GetProcessMemoryInfo`.
- [x] Initialize Windows `PROCESS_MEMORY_COUNTERS.cb` to its structure size; the current code leaves it zero.
- [x] Replace Windows console declarations with `windows-sys`, preserving the upstream behavior of addressing stdout even when `err` is true.
- [x] Add native runtime tests asserting nonzero current/peak RSS and checking that touching a 16–64 MiB allocation raises current RSS within broad scheduling/allocator tolerances.
- [x] Test realistic Linux delayed-spill behavior and selected-limit output
  equivalence before and after the ABI rewrite.
- [ ] Bracket the exact pre/post-rewrite transition boundary with repeated
  alternating runs, and validate native macOS/Windows spill decisions.
  - [x] Linux: establish the scheduler-robust pre/post bracket.
  - [ ] Validate native macOS and Windows spill decisions.
- [x] Add a target-native current-tree spill-conformance driver and CI matrix
  for Linux, macOS arm64/x86-64, and Windows x86-64.
- [ ] Retain successful native CI artifacts for every target before treating
  the cross-platform spill-decision requirement as complete.

Linux post-rewrite behavior has been characterized on natural, nonduplicated AMR data. With 2,000 protein queries against 4,000 references, `blastp` at 64 MiB passed its initial in-memory check and later migrated accumulated hits to disk at both one and four threads. The default 16 GB, forced `0G`, and adaptive 64 MiB outputs were byte-identical (39,909 rows; SHA-256 `d324e3efeaa48cff4885dfb4ad9c880a699b5f9a8c372b97f84629e2f25ccfc9`). Peak RSS was 98,720/99,248 KiB at 16 GB, 38,504/48,300 KiB at `0G`, and 69,484/69,940 KiB at 64 MiB (one/four threads).

The natural 250 kb *E. coli* `blastx` fixture also produced byte-identical 16 GB, `0G`, and 16 MiB outputs, but it retained only 491 hits. Its small-limit transition therefore occurred at the first RSS check; it did not accumulate the 8 MiB required for a later check. Limits from 64–88 MiB did not spill. This fixture is recorded as insufficient evidence for a genuine delayed `blastx` spill, not as a pass; the larger natural SRR1703350/Swiss-Prot matrix in Phase 8 subsequently establishes that Linux case. Raw timings, RSS logs, outputs, and transition diagnostics are retained under `.tmp/spill-boundary/results/`. Cross-platform and pre/post-boundary equivalence remain open, so the checklist item stays unchecked.

The Linux pre/post boundary is now bracketed on the natural 2,000-query/4,000-reference `blastp` fixture at four threads. Four alternating runs per binary and limit, reversing checkpoint/rewrite order on even repeats, gave 8/8 spills at `96M` and 0/8 spills at `104M`. Every `96M` run performed a genuine delayed transition and read back 20.6--20.7 MiB of compressed spill data. The closer nominal limits are not stable process-RSS boundaries on this workload: at `99M`, the checkpoint spilled in 1/4 runs and the rewrite in 2/4; at `100M`, each spilled in 2/4. That mixed result is retained rather than overfit to a single scheduler/allocator sample. All 32 outputs were byte-identical (39,909 rows; SHA-256 `d324e3efeaa48cff4885dfb4ad9c880a699b5f9a8c372b97f84629e2f25ccfc9`) and all runs exited zero. Thus `96M`/`104M` is the narrowest repeated scheduler-robust Linux bracket established by this pass, and the ABI rewrite does not change the spill decision at either stable side. Timing is recorded but is not used as boundary evidence because other work was present on the host. The exact driver, summaries, full diagnostics, hashes, timings, and provenance are retained under `.tmp/spill-boundary/run_prepost_boundary.sh`, `.tmp/spill-boundary/prepost-boundary-results/`, and `.tmp/spill-boundary/prepost-robust-boundary-results/`. Native macOS and Windows validation remains open, so the parent checklist item remains unchecked.

`scripts/run_native_spill_conformance.sh` now supplies the portable
current-tree gate. It downloads the immutable NCBI AMRFinderPlus
`2025-12-03.1` protein release and verifies its SHA-256 before selecting 4,000
reference and 2,000 query records evenly across the source. The selected FASTA
hashes are pinned too, and the selector rejects duplicate identifiers. On each
native runner the driver discovers the RSS transition instead of reusing a
Linux limit: it establishes the startup-spill and no-spill sides, searches for
a run that migrates already-retained hits and reads a nonzero compressed
payload, moves into the interior of the observed delayed-spill interval, then
requires two repeated delayed spills. Two `16G` runs must remain in memory.
Every calibration and gate output must be byte-identical to the first no-spill
output. Commands, stdout, stderr, output files, hashes, selected limit, binary
hash, fixture hashes, and host provenance are uploaded by the
`native-spill-conformance` CI job for Ubuntu, macOS arm64, macOS x86-64, and
Windows x86-64. The same driver is also run against the statically linked
Linux musl release binary. Hashing selects `sha256sum` or the stock macOS
`shasum -a 256`, and artifact-directory checks use Bash rather than relying on
GNU `find`. A `conformance-passed.txt` marker containing hashes of the summary
and output-hash manifest is written only after both delayed-spill repeats
pass; output paths in the hash manifest are artifact-relative. The driver
validates the complete relocatable artifact before reporting success, and the
same validator checks target metadata, vector lengths, and C++/Rust raw-vector
identity for RNG artifacts before upload. Partial diagnostic artifacts are
therefore not completion evidence.

The driver itself passed end to end on Linux in a local release run. Its
target-native calibration selected `64M`; both gate repeats performed genuine
delayed spills, migrating 2,084,609 and 1,539,946 retained hits and reading
20.9 and 21.0 MiB. Both `16G` runs and all 13 calibration/gate runs produced
the same output SHA-256,
`4b7f0c4ca2d566f4efa0edf8bbe5659ac2ce350e0110e90564e81db5911ceaeb`.
This validates the driver mechanics only. It is not macOS/Windows evidence,
and the CI job compares current-tree spill and no-spill modes rather than the
historical checkpoint binary. The earlier Windows checkpoint has a known
zero-initialized `PROCESS_MEMORY_COUNTERS.cb`, so treating its RSS decisions as
a portable oracle would preserve a known defect. The pre/post and native
macOS/Windows checkboxes therefore remain open until native artifacts are
available and reviewed.

After the portable hash and evidence-marker hardening, the complete driver was
rerun locally against the pinned natural fixture. It selected `65M`; the two
gate runs migrated 2,063,585 and 2,084,609 retained hits and read 20.9 MiB
each. Both `16G` controls and all fourteen outputs were byte-identical with
SHA-256 `6922e5b64db2d4d10217d94e3ee7e9c060d08b324b1d0cf5eef79fc5088f59ab`.
The generated success marker and its summary/output-manifest hashes were
verified. The standalone artifact validator accepted the relocatable complete
spill and GNU RNG captures and rejected marker-corrupted copies of each. This
reconfirms the current Linux driver but does not close the unchecked native
macOS/Windows gates.

Acceptance: no handwritten Mach or Windows ABI layouts/declarations remain; Linux, macOS, and Windows native runtime tests pass; current and peak RSS retain their platform units and semantics.

## Phase 6 — safe SQLite boundary, same engine

This phase reduces local unsafe code, not the native dependency.

- [x] Replace direct `libsqlite3-sys` connection/statement ownership with `rusqlite` or an equivalently maintained safe API.
- [x] Retain bundled SQLite initially for reproducible Windows/macOS/Linux builds.
- [x] Preserve readonly behavior, missing-parent result, maximum taxid, cache behavior, schema errors, and corrupt-database errors.
- [x] Validate against a real `taxonomy4blast.sqlite3` database and taxonomy-aware output fields.
- [x] Measure cold open and repeated parent lookups.
- [x] Correct README taxonomy claims based on the functionality actually validated.

Acceptance: application code contains no raw SQLite pointers or manual finalize/close calls, with exact query-result parity.

The ignored release diagnostic `benchmark_sqlite_taxonomy_open_and_parent` uses the same schema and query strings as production. On this Linux host, 1,000 read-only opens plus `max(taxid)` took 95.948 ms (95.948 us/op), and 100,000 prepared parent lookups took 1.222 s (12.225 us/op). It checks result checksums and can be reproduced with `cargo test --release --lib benchmark_sqlite_taxonomy_open_and_parent -- --ignored --nocapture`.

The external-data gate was run on 2026-10-06 against NCBI's actual `taxdb.tar.gz` (`dd77625e0a2073257f7aa4adaac8e99357ad96947006abbe2cea539be2f15ad8`) and `new_taxdump.tar.gz` (`64fff11538f8a25d5cf497417ce7b0347c937fa56be476a1cceb463ced4eb892`). NCBI BLAST+ `makeblastdb` 2.12.0+ produced a two-record v5 protein database with taxids 9606 and 562. `validate_real_ncbi_taxonomy4blast_database` verified stable parent rows and the real maximum taxid; `validate_real_ncbi_blastdb_taxonomy_mapping_and_parents` verified the upstream-generated volume, embedded taxid mapping, and SQLite parent lookup together. Finally, upstream DIAMOND 2.1.24 generated `staxids`, `sscinames`, `slineages`, `sskingdoms`, `skingdoms`, and `sphylums` rows from those same files, and `validate_real_ncbi_taxonomy_output_fields_against_upstream` matched every translated field byte-for-byte. This exposed and fixed one fidelity bug: C++ retains the taxid `0` sentinel when a requested rank is absent, while `TaxonomyTree::rank_taxids` had filtered it out.

These tests are ignored because the current NCBI inputs total hundreds of MiB. Their environment variables and exact invocation contract are documented beside each test in `src/data/blastdb/blastdb.rs`. This validates the translated reader and formatters; it does not claim that the native `blastp`/`blastx` command path now carries taxonomy metadata. That separate command-integration gap remains in `GAPS.md` and the README.

## Phase 7 — deterministic cross-platform RNG

- [x] Inventory every translated C `rand`/`srand` consumer and preserve draw order.
- [x] Establish raw GNU/Linux vectors for seeds `1` and `12345`.
- [x] Replace the target C runtime with a specified Rust implementation of the GNU 31-word additive generator.
- [x] Fix `RAND_MAX` at `i32::MAX` on every target.
- [x] Snapshot `standard_normal`, TANTAN lambda/probabilities, and `initParametersWithErrors` in isolated processes.
- [x] Execute the same raw and downstream suite against glibc and musl.

The active contract is target-independent: `compat_rng` uses only `std`,
retains the established GNU sequence, and serializes one process-global state.
This intentionally differs from upstream DIAMOND on C runtimes whose
`rand` algorithm or `RAND_MAX` differs. Concurrent consumers still share
draw order, matching the observable global-state model.

The former target-native C++ RNG oracle, its source fixture, its musl driver,
and its CI matrix were removed on 2026-10-07. They tested the old
platform-runtime contract and would now reject the intended deterministic
behavior. `tests/compat_rng_vectors.rs` is self-contained and covers both raw
integers and downstream numerical results.

The upstream `run/tools.cpp::random_seqs` utility remains a separate command
surface gap: the translated tools facade delegates selection to
`ToolsBackend::random_index`, and no production backend implements it.

Acceptance: the generator sequence and normalization range are identical on
all supported targets, while the established GNU/Linux downstream results
remain unchanged.

## Phase 8 — retain and measure GNU `malloc_trim`

- [x] Keep `malloc_trim` behind `cfg(all(target_os = "linux", target_env = "gnu"))` and one documented safe wrapper, using the maintained `libc` declaration rather than a handwritten ABI block.
- [x] Add an A/B feature or test hook without changing the production default.
- [x] Compare the default adaptive-memory policy and forced disk `0G` for `blastp`, using one and four threads.
- [x] Require byte parity and report two-repeat mean wall time and peak RSS for the scoped A/B workload.
- [x] Decide from the measured evidence: retain `malloc_trim`; the tested allocator/lifetime design does not reach equivalent RSS when trimming is disabled.
- [x] Exercise a nonzero adaptive limit that spills in the middle of a `blastp` run at one and four threads.
- [x] Repeat the trim/no-trim `blastx` A/B matrix for default `16G`, forced disk `0G`, and a documented `16M` startup transition at one and four threads.
- [x] Find a natural, nonduplicated `blastx` workload that produces a genuine nonzero delayed spill after retained hits accumulate, then repeat its trim/no-trim A/B.

The nonzero delayed-spill case is established for `blastp` at 64 MiB and both one/four threads, as described in Phase 5. The initial small-fixture `blastx` A/B below has a `16M` transition at the first RSS check and reads back 0.0 MiB, so it is explicitly not treated as delayed-spill evidence. The later SRR1703350/Swiss-Prot matrix closes that requirement with a measured nonzero in-memory migration.

### `malloc_trim` A/B evidence

The Linux x86-64 GNU fixture used 300 protein queries against the existing 4,000-sequence, 1,357,780-letter DMND database. Two release binaries were built in separate target directories: the production default and `--features disable-malloc-trim`. Each cell ran twice with trim/no-trim order reversed:

```text
diamond blastp \
  -q .tmp/validation-fixtures/amr-search/query-300.faa \
  -d .tmp/bench-tmp/diamond-compare.H22XYK/rust-db.dmnd \
  -o RESULT.tsv --threads {1,4} --tmpdir .tmp/malloc-trim-ab/scratch \
  [--memory-limit 0G]
```

All 16 result files were byte-identical, with SHA-256 `1b589d3ab5e7f1779707dcc3eced0ef381131ca266c0f9ff0870fb528cdfa776`.

| Memory mode | Threads | Trim | Mean time | Mean peak RSS |
|---|---:|---|---:|---:|
| Default adaptive policy | 1 | enabled | 4.090 s | 25,558 KiB |
| Default adaptive policy | 1 | disabled | 3.910 s | 26,000 KiB |
| Default adaptive policy | 4 | enabled | 1.500 s | 29,002 KiB |
| Default adaptive policy | 4 | disabled | 1.560 s | 35,116 KiB |
| Forced disk `0G` | 1 | enabled | 4.050 s | 17,188 KiB |
| Forced disk `0G` | 1 | disabled | 4.165 s | 17,766 KiB |
| Forced disk `0G` | 4 | enabled | 1.235 s | 23,846 KiB |
| Forced disk `0G` | 4 | disabled | 1.285 s | 28,714 KiB |

The default-policy fixture stayed below its 16 GB limit, so it exercised the adaptive in-memory path without spilling. The `0G` cases entered disk mode immediately. This initial A/B did not exercise a nonzero mid-run spill; the later 64 MiB boundary matrix in Phase 5 establishes that behavior for `blastp`. With only two repetitions these timings are diagnostic, but the four-thread RSS reductions are large enough that the evidence supports retaining `malloc_trim`; it does not support removal.

### `blastx` trim/no-trim A/B evidence

The natural 250 kb *E. coli* DNA fixture was searched against the same 4,000-protein database using the existing production-default and `disable-malloc-trim` release binaries. Each cell ran twice with order reversed on repeat two. All 24 result files were byte-identical, with SHA-256 `dc8c6af6abc52a0a506fcd945d2c8fa7288dab1514b602b12362c24ceb7d62a1`.

| Memory mode | Threads | Trim | Mean time | Mean peak RSS |
|---|---:|---|---:|---:|
| In-memory `16G` | 1 | enabled | 1.055 s | 39,848 KiB |
| In-memory `16G` | 1 | disabled | 1.260 s | 43,052 KiB |
| In-memory `16G` | 4 | enabled | 0.590 s | 39,716 KiB |
| In-memory `16G` | 4 | disabled | 0.590 s | 42,918 KiB |
| Forced disk `0G` | 1 | enabled | 1.045 s | 39,412 KiB |
| Forced disk `0G` | 1 | disabled | 1.055 s | 43,134 KiB |
| Forced disk `0G` | 4 | enabled | 0.610 s | 40,042 KiB |
| Forced disk `0G` | 4 | disabled | 0.590 s | 43,402 KiB |
| Startup transition `16M` | 1 | enabled | 1.030 s | 39,538 KiB |
| Startup transition `16M` | 1 | disabled | 1.025 s | 43,010 KiB |
| Startup transition `16M` | 4 | enabled | 0.610 s | 40,058 KiB |
| Startup transition `16M` | 4 | disabled | 0.630 s | 43,360 KiB |

The timings are diagnostic because this fixture completes in 0.6–1.5 seconds and one no-trim `16G` one-thread sample was visibly interrupted. The RSS signal is nevertheless consistent in every cell: trimming lowered mean peak RSS by 7.5–8.6%, while byte output stayed exact. Both `0G` and `16M` logged the memory-limit transition at their first check and read back 0.0 MiB of compressed spill data; `16G` did not spill. Consequently this matrix strengthens the decision to retain trimming but does **not by itself** satisfy the delayed-spill `blastx` item. Exact commands, alternating order, timings, outputs, diagnostics, hashes, and binary/input provenance are retained under `.tmp/malloc-trim-ab/blastx-results/` and `.tmp/malloc-trim-ab/run_blastx.sh`.

### Delayed-spill natural `blastx` evidence

The remaining delayed-spill case uses both existing 5,000-record mates from public paired-end run SRR1703350 and the existing 574,627-sequence, 208,482,574-letter Swiss-Prot database. These are 10,000 distinct natural DNA reads; no sequence was copied, repeated, or synthesized to enlarge the workload. Default-sensitivity search retained 764,440 stage-2 hits and reported 79,648 alignments. A `1160M` calibration stayed in memory, while `1150M` passed the first RSS check and later migrated 526,641--547,543 already-retained hits before reading 5.1 MiB of compressed spill data. The transition is therefore a genuine nonzero delayed spill rather than the earlier startup approximation.

Production-default and `disable-malloc-trim` release binaries each ran twice at one and four threads, with their order reversed on repeat two:

```text
diamond blastx \
  --query subset_1.fastq --query subset_2.fastq \
  --db sprot_rust.dmnd --out RESULT.tsv \
  --threads {1,4} --memory-limit 1150M --tmpdir RUN_TMP
```

| Threads | Trim | Mean time | Mean peak RSS |
|---:|---|---:|---:|
| 1 | enabled | 99.370 s | 1,128,454 KiB |
| 1 | disabled | 97.340 s | 1,130,104 KiB |
| 4 | enabled | 31.945 s | 1,132,032 KiB |
| 4 | disabled | 31.685 s | 1,134,700 KiB |

All eight outputs were byte-identical (SHA-256 `41583b9fca760e27aa63f73633dfc43dc7e3194fb43e7e9fb9bc011f12183efd`). Every run migrated more than 526,000 retained hits and read back 5.1 MiB. At this scale the immutable Swiss-Prot database dominates RSS, so trimming saves a smaller 1,650 KiB (0.15%) at one thread and 2,668 KiB (0.24%) at four threads. It costs 2.1% and 0.8% wall time respectively in these two diagnostic repeats. This does not overturn the retain decision: it shows that the large fixed database footprint dilutes the allocator-release benefit seen consistently in the smaller complete matrix. Exact timings, RSS, stderr diagnostics, outputs, hashes, binary/input checksums, and the alternating driver are retained under `.tmp/malloc-trim-ab/blastx-delayed/results/` and `.tmp/malloc-trim-ab/run_blastx_delayed.sh`.

The retained call now goes through `libc::malloc_trim`. This is the same glibc symbol and behavior measured above, but it removes the final handwritten production libc declaration. Acceptance: retain the call unless measurements justify removal. Changing the global allocator is a separate project.

## Deferred native dependencies and oracle

These are not blockers for eliminating handwritten libc/platform ABI code:

- [x] Replace the zstd C backend with `zstd-pure-rs` after validating C++↔Rust interoperability, truncated/concatenated/empty frames, throughput, peak RSS, and compressed size.
- [x] Keep bundled SQLite through the safe-wrapper phase. A pure-Rust SQL replacement would require its own format/correctness evaluation.
- [x] Keep upstream `diamond_core` behind the default-off `ffi` feature as a conformance oracle until all native routes and the validation matrix are complete.
- [x] Keep default CI and benchmarks free of `--features ffi`; remove the
  adjacent-upstream oracle from routine CI and use the external comparison
  harness for parity work.
- [x] Inventory all `--legacy` routes before proposing removal of the C++ adapter.

### C++ oracle route inventory (historical checkpoint)

This inventory describes the ABI-reduction checkpoint, not current command
coverage. Native SEG multiple-HSP and standard composition-mode search integration were completed in
November 2026; current evidence and remaining routes are tracked in
`TRANSLATION_COMPLETION.md` and the executable's routing code.

The default Ubuntu/macOS/Windows CI matrix runs `cargo build --all-targets` and
`cargo test --all-targets` without `ffi`. macOS is represented by both
`macos-latest` (ARM64) and `macos-15-intel` (x86-64). Source-audit tests and
the CI-side upstream checkout were removed on 2026-10-07; parity comparisons
now belong to the external benchmark harness rather than the repository test
graph. The optional `ffi` feature remains available for manual diagnosis only
from a source checkout containing the adjacent `diamond/` tree. That tree is
intentionally absent from packaged crates. Packaged all-feature builds compile
with the C++ adapter unavailable; production legacy routes still return an
error. Removing those production routes remains an active requirement.

The shared native parser/dispatcher is `src/commands/dispatch.rs`; the small
`src/main.rs` executable bridge routes its untranslated results to `run_legacy`
in the following cases. Public `diamond::run` uses only native dispatch and
returns an explicit error for those results. Native parser/command failures
return statuses without exiting the embedding process. The bridge remains a
production dependency until its valid routes are translated. Without the
non-Windows `ffi` build this route reports an error rather than changing engines.

Remaining executable routes include:

- An explicit `--legacy` on commands with a native arm.
- Any command not matched by a native arm. Native arms currently cover `version`, help, `dbinfo`, `makedb`, `getseq`, `makeidx`, supported `blastp`/`blastx`, `cluster`/`linclust`/`deepclust`, `merge-daa`, supported `view` formats, `greedy-vertex-cover`, `realign`, supported `reassign`/`recluster`, `prepdb`, `roc`, `benchmark`, and `test`.
- `view` formats other than pairwise (`0`), XML (`5`/`xml`), tabular (`6`/`tab`), DAA (`100`/`daa`), SAM (`101`/`sam`), PAF (`103`/`paf`), flat JSON (`104`/`json-flat`), `null`, and `edge`.
- `blastp`/`blastx` with an unknown option; a non-native output format or invalid output field; NCBI BLAST input; nonzero `--global-ranking`; contiguous/combined indexing; or remaining controls such as `--soft-masking`, `--approx-id`, `--file-buffer-size`, `--header`, `--compress`, or `--iterate`. Blocked taxonomy fields now use the concrete native metadata adapter. `--unal` is connected for transcript, tabular, and flat JSON formats; explicit DAA unaligned handling remains open. All standard tabular/JSON fields and dynamic taxonomy ranks use the native context adapter. See the live route predicates and `TRANSLATION_COMPLETION.md` for the full evolving inventory.
- SEG masking, unlimited/multiple HSPs, all standard CBS modes with source blastx restrictions, `--index-chunks`, explicit `--algo 0/1`, and `--query-parallel-limit` validation now run natively. Explicit `--block-size` uses actual protein/translated file blocks and disk intermediate merging; the source default budgets and supplied translated-block API are now connected for native FASTA/DMND search. BLAST/clustering blocking and scale proof remain open. `diamond test` now passes all twenty source cases natively, including blocking.
- Additionally for `blastx`: nonzero `--frameshift`, `--swipe`, or any `--min-query-len` request.

### zstd backend evidence

The interoperability probe compiled the checked-in upstream implementation itself, rather than substituting the `zstd` command-line tool:

```text
c++ -std=c++17 -O3 -DNDEBUG -ffunction-sections -fdata-sections \
  -Idiamond/src cpp_probe.cpp diamond/src/util/io/zstd_stream.cpp \
  -Wl,--gc-sections -lzstd -pthread -o cpp_probe
```

The probe linked upstream `ZstdSink` and `zstd_decompress` to system libzstd 1.4.8. The initial Rust probe used `zstd-sys` 1.5.7; the current adapter uses `zstd-pure-rs` 0.2.0 with the same translated `ZstdSink`/`zstd_decompress` interfaces. C++-encoded data decoded byte-for-byte in Rust and Rust-encoded data decoded byte-for-byte in C++ for an empty input, the natural 4,695,584-byte AMRProt FASTA, and concatenated two-frame forms of that FASTA. Both whole-file decoders rejected frames truncated by three bytes. The AMRProt input and each decoded output had SHA-256 `41d5ebf4f807c9590f27de7dd23c9f037afb4d8efd613790c54847fb1fc2896b`. A fixed frame emitted by the upstream C++ stream is also retained as a normal Rust whole/stream decoder unit test.

The pure-Rust migration was rechecked on Linux with that natural AMRProt fixture. Both cross-decode directions were byte-exact, and the 694,480-byte pure-Rust frame was byte-identical to the frame previously emitted through `zstd-sys` 1.5.7. Version 0.1.2's higher-level `ZSTD_decompressStream` staged a complete compressed and decoded frame, which could make RSS scale with large DIAMOND spill files. Version 0.2.0 contains the bounded block-streaming fix, including reusable staging, partial-input cursors, window-limit checks, magicless streams, concatenated/skippable frames, and dictionary lifecycle handling; all 1,089 library tests pass. DIAMOND has upgraded to the published 0.2.0 crate but temporarily retains its already-validated low-level block adapter. Migrating that adapter to the repaired high-level API is a separate cleanup requiring another spill/RSS gate.

The fix was committed as `zstd-pure-rs` commit
`7f34eae2a0ea7a328f47311b5c11893f8795630e` (`streaming`) and released in
0.2.0. The release uses 0.2 rather than 0.1.3 because the repair added state to
the public `ZSTD_DCtx`; that type is now non-exhaustive so future internal
decoder-state additions do not repeatedly break downstream construction.

Ten focused DIAMOND stream tests passed, including empty, concatenated, truncated, rewind, small-output, small-read, bounded staging, and large temporary-file cases; the 22 general I/O tests and three TSV tests also passed. The full target-native spill gate then passed on the pinned natural AMRProt selection: all no-spill, startup-spill, calibration, and two repeated delayed-spill outputs had SHA-256 `6922e5b64db2d4d10217d94e3ee7e9c060d08b324b1d0cf5eef79fc5088f59ab`. Both selected 65 MiB delayed-spill runs migrated 2,063,585 retained hits and read 20.9 MiB of compressed spill data. Artifacts are under `/tmp/diamond-zstd-pure-spill-20261006` on the validation host.

Three release probe runs gave pure-Rust encode wall times of 0.03–0.05 s at 15,360 KiB RSS versus 0.05–0.06 s and 13,760–14,376 KiB for upstream C++/libzstd. Decode timings were too short and noisy for a speed claim; peak RSS was 7,360–7,680 KiB for Rust and 10,240 KiB for C++. This is a focused regression screen, not broad benchmarking.

A modest release-mode comparison used the natural 11,416,585-byte `calibrate-blastp-64m-t4.tsv` search result. Seven runs per implementation alternated order; the table reports medians from `/usr/bin/time` and includes file I/O. Compression compares each implementation's stream encoder. Decompression deliberately crosses implementations (C++ decodes the Rust frame and Rust decodes the C++ frame). The source and decoded outputs shared SHA-256 `d324e3efeaa48cff4885dfb4ad9c880a699b5f9a8c372b97f84629e2f25ccfc9`.

| Operation | Implementation | Median time | Median peak RSS | Compressed size |
|---|---|---:|---:|---:|
| Encode | Upstream C++/libzstd 1.4.8 | 0.09 s | 20,392 KiB | 682,564 bytes (5.9787%) |
| Encode | Rust wrapper/zstd-sys 1.5.7 | 0.05 s | 19,520 KiB | 681,412 bytes (5.9686%) |
| Cross-decode | Upstream C++ decoding Rust | 0.11 s | 16,960 KiB | — |
| Cross-decode | Rust decoding C++ | 0.12 s | 16,320 KiB | — |

The workload is intentionally modest and sub-second, so these timings are a regression screen rather than a broad performance claim. This table is the pre-migration native-backend baseline; the focused pure-Rust measurements above are the migration checkpoint. Compressed-byte identity is not required across different zstd implementations, although the migrated encoder did reproduce the prior Rust frame exactly on the natural AMRProt fixture; decoded-byte identity remains the format contract.

## Validation matrix for every wave

- [x] `cargo fmt --all -- --check`
- [x] `cargo check --all-targets`
- [x] Full `cargo test --lib` and integration tests on Linux.
- [ ] Native tests, not build-only checks, on Ubuntu, macOS, and Windows CI.
- [x] Linux x86-64 GNU runtime coverage.
- [x] Linux musl build and runtime coverage where available.
- [ ] macOS ARM64 runtime plus x86-64 compile/runtime where available.
- [ ] Windows x86-64 MSVC runtime coverage.
- [x] Selected 32-bit Linux/Windows compile checks to catch layout and `usize` assumptions.
- [x] Exact output parity on the affected representative search and clustering commands.
- [x] Repeated release-mode wall-time and peak-RSS comparison on realistic data before and after the completed rewrite wave.
- [x] Update this checklist, `GAPS.md`, validation results, and stable README benchmark tables as evidence is produced.

### Portability CI configured, awaiting run evidence

The workflow keeps three portability checks separate from the native
Ubuntu/macOS/Windows test matrix:

- `linux-musl` builds every target for `x86_64-unknown-linux-musl`, executes
  the Rust tests, and runs the same target-native spill gate as the primary OS
  matrix. The obsolete C-runtime RNG oracle is no longer part of this job.
- `linux-32-bit` installs `gcc-multilib` and runs `cargo check --target i686-unknown-linux-gnu --all-targets`.
- `windows-32-bit` runs `cargo check --target i686-pc-windows-msvc --all-targets` on a native Windows runner, allowing bundled SQLite/zstd C build scripts to use the matching MSVC environment.

None of these jobs enables `ffi` or downloads the upstream source tree.
Selected 32-bit compile coverage is checked above from the successful local
GNU target checks described below. The independently configured native i686
MSVC job and the overall native-CI completion gate remain unproven until
GitHub Actions completes successfully; local workflow parsing and installed-
target checks do not prove native MSVC or runtime gates.

Linux musl was additionally exercised locally, rather than inferred from the
workflow definition. With Rust's `x86_64-unknown-linux-musl` standard library
and Zig 0.17 as the musl C/linker toolchain, `cargo zigbuild --offline --target
x86_64-unknown-linux-musl --all-targets` produced statically linked binaries.
Every generated normal test harness was then executed natively on this x86-64
Linux host. The library result was 1,712 passed and 6 ignored; the binary and
normal integration suites also passed (7, 8, 1, 1, 6, and 4 tests in their
nonempty harnesses). Data-dependent diagnostic suites remained ignored by
their existing contracts. The historical target-C-runtime RNG comparison was
later superseded by the deterministic Rust generator described in Phase 7;
its script and fixture are no longer part of the repository. Build products are retained under
`.tmp/zig-target/`. This closes the local
musl build/runtime item; the configured GitHub Actions job must still run
successfully before CI itself can be claimed green.

The selected 32-bit gate was also run locally after installing the Rust i686
standard libraries. A full Zig-backed
`i686-unknown-linux-gnu --all-targets` build and
`cargo check --target i686-pc-windows-gnu --all-targets` both pass. This gate
found and fixed three genuine assumptions: the structural AVX2 traceback
module imported the x86-64 intrinsic namespace and used the x86-64-only
`_pdep_u64`; AVX-512 assembly used 64-bit general registers on 32-bit x86; and
the 16 GB search ceiling was stored as `usize`. The x86-64 kernels remain
unchanged; 32-bit x86 now uses its existing portable traceback tier and AVX2
Hamming path, while the AVX-512 path is correctly x86-64-only. Search memory
limits are now represented as `u64`, matching the CLI domain without changing
RSS or hit-buffer sizes. The attempted 32-bit Linux executable cannot be used
as runtime evidence in this sandbox (execution exited under its syscall
policy), so only compile coverage is claimed. Native i686 MSVC CI remains an
independent configured gate.

The current workflow has five self-contained jobs: the four-OS native build
matrix, the four-OS spill matrix, Linux musl, Linux i686, and Windows i686.
The native matrix explicitly covers both hosted macOS architectures.
`CARGO_TARGET_DIR=.tmp/portability-check cargo check --target
x86_64-pc-windows-gnu --all-targets` also passed with the installed GNU Windows
cross-toolchain. That check does not substitute for native MSVC or macOS runs.

The full formatting gate initially exposed a rustfmt 1.8.0 nontermination in the pre-existing mechanically translated `src/lib/blast/blastn_score.rs`. Every source file was checked individually to isolate it. The already-formatted 1,600-line unit is now explicitly skipped at its module declaration, allowing `cargo fmt --all -- --check` to complete successfully while every changed source and test remains covered.

### Post-rewrite broad parity rerun

The worktree was rebuilt in release mode after the FFI rewrite and compared
directly with upstream DIAMOND 2.1.24.  The broad-matrix Rust binary has
SHA-256 `0aa98a5adcf5fff86a33af3d6a3a18fdd2b5676f8c2f0bfcef13a1bcbfd3df07`.
The exact invocation was:

```text
cargo build --release --bin diamond
env VALIDATION_AREAS='sensitivity clustering scoring view' \
  RESULT_DIR=.tmp/ffi-final-parity-current \
  PROTEIN_REFERENCE_FASTA=.tmp/validation-fixtures/amr-search/reference-1000.faa \
  PROTEIN_QUERY_FASTA=.tmp/validation-fixtures/amr-search/query-300.faa \
  CLUSTER_FASTA=.tmp/validation-fixtures/amr-cluster/reference-500.faa \
  NUCLEOTIDE_QUERY_FASTA=.tmp/validation-fixtures/ecoli-250k.fna \
  RUST_BIN="$PWD/target/release/diamond" \
  CPP_BIN="$PWD/diamond/build/diamond" \
  bash scripts/run_validation_plan.sh
```

This was a correctness matrix, not a stable timing benchmark.  It retained the
exact command, stdout, stderr, output, timing, RSS, and comparator result for
each case.  The natural protein fixtures contain 1,000 references/339,123
residues and 300 queries/102,101 residues; clustering uses 500 natural AMR
proteins/179,068 residues; blastx uses one natural 250,000-base genomic
region.  Fixture hashes are in the retained manifest.

| Area | Cases | Result against upstream |
|---|---:|---|
| Search sensitivity | 40 | 40 pass with byte-identical tabular output: every preset from true default through ultra-sensitive for blastp and blastx on quick and natural data, plus the selected four-thread and `--memory-limit 0G` cases |
| Scoring | 50 | 50 pass: 42 accepted cases have byte-identical output across all eight supported matrices/default and alternate gaps; 8 rejected-gap cases match upstream rejection semantics |
| Clustering, compact fixture | 6 | 6 pass for `cluster`, `linclust`, and `deepclust`, default and `--approx-id 90 --member-cover 80` |
| Clustering, natural fixture | 12 | 2 pass (`deepclust` default at one and four threads); 10 fail with genuine centroid/member differences |
| DAA production and view | 56 | 56 pass; flat JSON is byte-identical, including pinned upstream 2.1.24's missing inter-query separators |

The aggregate is 154 passes, 10 failures, and 0
blocked cases.  The natural clustering failures are not ordering artifacts:
the comparator canonicalizes centroid/member pairs by FASTA identifier.  For
example, default `cluster` has the same coarse counts in both implementations
(169 clusters, 500 members, largest cluster 48, 116 singletons), but different
canonical pair hashes.  At `--approx-id 90 --member-cover 80`, upstream has
254 clusters/195 singletons while Rust has 260/201.  The ten failing cases and
their canonical hashes are recorded in each `compare.log`.

Broad-matrix evidence is retained under `.tmp/ffi-final-parity-current/`, with
the fixture and toolchain identity in `manifest.tsv` and the machine-readable
result set in `summary.tsv`.  The byte-parity correction touched only the DAA
viewer. Its rebuilt binary has SHA-256
`82117f711cae2e6c21caf6932911a793ea60317d20f398aa371aa88130977122`;
the complete 56-case focused rerun is retained under
`.tmp/ffi-final-view-byte-parity/`. Together these reruns close the previously
missing natural DAA/view validation pass and reconfirm search and scoring
after the FFI changes. At that stage the broad exact-parity completion
criterion was not met because ten natural clustering cases remained real
fidelity gaps; the later follow-ups below close those cases and rerun the
entire matrix.

A focused clustering rerun on 2026-10-06, after restoring upstream's
cluster-only Hamming-extension configuration flow, improved the natural
fixture from 2/12 to 7/12 passing cases. The complete clustering area produced
13 passes and 5 failures across 6 compact and 12 natural cases, with no blocked
cases. Natural `cluster --approx-id 90 --member-cover 80`, `linclust` default
at one and four threads, `linclust --approx-id 90 --member-cover 80`, and
`deepclust --approx-id 90 --member-cover 80` now pass. The five remaining
natural failures are `cluster` default at one and four threads and the three
commands at `--approx-id 50 --member-cover 60`. Commands, outputs, canonical
comparisons, timing, and RSS are retained under
`.tmp/ffi-cluster-hamming-parity/`. A single-round natural
`--cluster-steps faster_lin` characterization is byte-identical (923
alignments, 195 clusters, 500 pairs, 141 singletons; SHA-256
`63defdc92c491b5e117e3dc9447237ff458ea568faffb9f43e545f06b1ab1e85`
for both TSV files). At that diagnostic stage exact clustering parity remained
open; the later round-order and endpoint-kernel follow-up below closes it for
the 18-case clustering matrix.

A fresh release rebuild containing the faithful double-precision coverage
calculation was then used to rerun the five failures plus a passing control.
The control still passed and all five failures remained the same single
42-member centroid swap. Round 1 now agrees in aggregate (for default: 923
HSPs, 1,831 directed edges, and 195 clusters). The first measured divergence
is round 2: upstream/Rust HSP counts are 30/28 for default, 31/30 for the
identity-50 `cluster`/`deepclust` path, and 36/34 for identity-50 `linclust`.
This narrows the remaining work to exact round-2 active-OID, seed-hit, and HSP
comparison; it does not justify changing thresholds or scores.

A subsequent faithful round-2 and endpoint-kernel translation closes that
diagnostic gap. Rust had re-sorted every linear round by sequence length,
whereas upstream preserves active OID order after the first round when
`keep_target_id` is set. Restoring that ordering and the literal mutual-cover
predicate removed the centroid/member differences. The no-trace endpoint path
now follows upstream's 32-lane i8, 16-lane i16, and scalar i32 promotion
cascade, with moving strict bands, SIMD row counters, per-lane adjusted
matrices/gap scaling, and dedicated full-matrix kernels.

Two fresh clustering matrices pass 18/18 exact comparisons against upstream:
all six compact and all twelve natural `cluster`, `linclust`, and `deepclust`
cases, including one/four threads and the identity-50/coverage-60 cases. The
first is retained under `.tmp/cluster-endpoint-structural-final/`; its repeat
is under `.tmp/cluster-endpoint-structural-repeat2/`. The structural SIMD
rewrite reduced aggregate Rust time from the retained scalar baseline's
37.68 s to 30.91 s and 26.83 s in the two runs (1.22x and 1.40x), while maximum
RSS was 9,280 and 8,760 KiB versus the baseline's 8,688 KiB. These are
shared-host diagnostic timings, not publication-grade benchmark medians.

The final rebuilt binary was then validated from scratch in three independent
area runs. Search sensitivity passed 40/40 cases, scoring passed 50/50,
clustering passed 18/18, and DAA production/view passed 56/56: 164/164 exact
comparisons, with no failures, blocked cases, or expected-difference labels.
The retained results are under
`.tmp/ffi-final-parity-fixed-search-scoring/`,
`.tmp/ffi-final-parity-fixed-clustering/`, and
`.tmp/ffi-final-parity-fixed-view/`. This closes the broad exact-parity
criterion for the explicitly validated matrix; it is not a claim about the
separate unchecked command/mode inventory in `GAPS.md`.

## Completion criteria

- [x] No handwritten OS/libc ABI declarations remain outside an explicit, reviewed parity-exception module.
- [x] No public raw descriptor or mmap pointer ownership remains.
- [ ] Linux, macOS, and Windows CI run relevant tests successfully.
- [ ] Search spilling still obeys the configured memory limit on every supported platform.
- [x] Exact parity remains established for the validated search, clustering, scoring, and DAA/view matrices.
- [x] Speed and RSS have not regressed in the completed representative rewrite matrix; every retained native call has a written parity or performance justification.


## Native SAM search checkpoint (2026-11-07)

`blastp`/`blastx` output `101`/`sam` no longer routes to C++ merely for its
format. Native SAM also carries `--unal`, `--sam-query-len` and SAM
`--salltitles` through real blocked search contexts. Source fixture coverage,
including the direct/intermediate title distinction, is recorded in
`tests/data/sam_search/` and `TRANSLATION_COMPLETION.md`. Advanced search
controls can still select production fallback; the overall FFI removal gate
remains open.


## Native edge/null search checkpoint (2026-11-07)

Search output `edge` and `null` now use concrete native callbacks with the
source coordinate/transcript DP requirements. Binary edges retain original
IDs and exact corrected-score bytes through block merging. The source
fixture matrix and two discovered numerical/translated-coverage corrections
are recorded in `tests/data/edge_null/` and `TRANSLATION_COMPLETION.md`.
Taxonomy format 102 and advanced search controls remain production fallback
work; the overall FFI-removal gate is not complete.

## Native taxonomy search checkpoint (2026-11-07)

`blastp`/`blastx -f 102` and `--include-lineage` now use the concrete native
classifier and actual metadata. All standard source search formats are
connected; advanced controls can still select production fallback. Classifier
top defaults, explicit top/score overrides, conditional metadata requirements,
and the source title-sorting undefined-behavior exception are documented with
bounded evidence in `tests/data/taxon_search/` and `TRANSLATION_COMPLETION.md`.
The whole-program FFI-removal gate remains open.


## Native search output controls checkpoint (2026-11-07)

Search headers, plain/gzip compression, filename suffix policy and file buffer
sizes now use the native OutputFile/Serializer lifecycle. The buffer setting
also reaches query/reference input, DMND readers and block intermediates.
Stock zlib through flate2 replaces zlib-rs to preserve source gzip bytes;
this C codec dependency is independent of the optional DIAMOND C++ adapter.
Lazy bounded scratch avoids duplicate large output allocations. Source option
protocol, header validation timing and error file side effects are preserved.
See `tests/data/output_options/README.md` for 348 CLI comparisons, eight exact
fixed-invocation library comparisons, large-stream checks and 21 error cases.
Invocation/XML-path substitutions are explicitly distinguished from raw-byte
comparisons. Advanced search fallback routes, optional modes and the outstanding
LambdaCalculator.cc inventory audit still prevent whole-program completion.


## Parallel search-control and TANTAN checkpoint (2026-11-07)

Approximate identity and explicit soft masking now connect to native CLI,
blocked search and supplied-query APIs, including hashed/translated mask
lifecycles. Focused tests pass 152 direct source comparisons and 12 separately
classified compatibility checks for an original SEG negative-index UB boundary.
Source variant outcomes are 17 successes and one further crash; expected
compatibility snapshots are corroborated by two independently completing variants.
The active LambdaCalculator numerical audit fixes decimal rounding and LU/solve
FMA, with 36 exact state/probability/RNG vectors and 19 formatter probes. Its .cc
row and explicit helper mappings complete the 199-file inventory. Configured
buffers now preflight impossible capacities safely. Details and limitations
are in TRANSLATION_COMPLETION.md and the self-contained fixture READMEs.
Full application integration and removal of production fallback remain open;
199 mapped files do not establish end-to-end completion.


## Code-first integration checkpoint (2026-11-07)

Native source now connects global ranking, iteration, exhaustive SWIPE,
legacy frameshift, ordinary extension modes, contiguous seeds, signed
context-aware minimum lengths, block-streamed BLAST references, baseline
registry preflight and getseq's real parser/source adapter. One shared
all-target Cargo check and native debug executable build pass. Broad
runtime/parity/platform/package checks are intentionally deferred per user
priority; these paths are code-connected, not newly parity-proven.
See TRANSLATION_COMPLETION.md for concrete phases, source-UB compatibility
choices and remaining requirements. Production fallback still exists for
other valid untranslated modes; the end-to-end goal remains active.
