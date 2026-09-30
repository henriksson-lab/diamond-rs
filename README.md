# diamond-rs

Rust port of the [DIAMOND](https://github.com/bbuchfink/diamond) protein sequence aligner.

DIAMOND is a high-performance sequence aligner for protein and translated DNA searches, designed for big sequence data analysis. This crate provides both a CLI binary and a library API.

**This crate is under translation. Do not use it. Do not trust any text below**

* 2026-09-30: More optimization
* 2026-09-28: Ondisk blastp mode, but also faster inmem mode (not in original). Parity on broader datasets
* 2006-09-27: Closing the gaps in translation, optimization. Some more work left
* 2026-08-01: CI added. Full translation is blocked until BLAST is translated
* 2026-07-07: New audit; state to be checked

## This is an LLM-mediated faithful (hopefully) translation, not the original code!

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this page](https://rewrites.bio/)
but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back to the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* Do not trust the benchmarks on this page. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* Check the original Github pages for information about the package. This README is kept sparse on purpose. It is not meant to be the primary source of information



## Status

This project is an ongoing port of the DIAMOND C++ codebase to Rust. Currently:

- **CLI**: `blastp` and `blastx` run natively in Rust by default; C++ FFI fallback is only built for non-Windows conformance testing with `--features ffi`
- **Native Rust commands**: `blastp`, `blastx`, `makedb`, `dbinfo`, `getseq`, `version`, `help`
- **Parallel**: Seed search uses rayon for multi-threaded processing
- **SIMD**: SSE4.1/AVX2 dynamic-programming kernels and runtime-selected AVX-512 search workers on supported x86-64 hosts
- **Library API**: Core types, scoring matrices, DP kernels, FASTA parsing, and seed search
- **Tests**: More than 1,600 passing library tests plus CLI and integration suites, including the C++ regression inventory and native-vs-FFI equivalence
- **Not yet translated**: SQLite-backed taxonomy lookup for NCBI BLAST databases (`taxonomy4blast.sqlite3`), used by taxonomy-aware output fields such as `slineages`, `sskingdoms`, `skingdoms`, and `sphylums`

## Building

### Prerequisites

- Rust 1.70+
- Default native Rust build: no CMake or C++ toolchain required
- Optional non-Windows FFI test build: CMake 2.6+, a C++ compiler, zlib, SQLite3, and pthreads. SQLite3 is currently required by the vendored C++ build, but SQLite-backed BLAST taxonomy lookup has not yet been translated into native Rust.

For the optional FFI build on Ubuntu/Debian:
```bash
sudo apt-get install g++ cmake zlib1g-dev libsqlite3-dev
```

### Build

```bash
cargo build --release

# Host-optimized benchmark build (matches upstream's -march=native build)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Build the non-Windows C++ FFI backend for conformance testing
cargo build --features ffi
```

### Test

```bash
# Run all tests (single-threaded for FFI safety)
cargo test --release -- --test-threads=1

# Run tests that compare against the vendored C++ FFI backend
cargo test --release --features ffi -- --test-threads=1

# Run just the 20 regression tests
cargo run --release --features ffi -- test
```

## CLI Usage

```bash
# Build a database from FASTA
diamond makedb --in reference.faa -d reference

# Protein-protein alignment
diamond blastp -q query.faa -d reference -o results.txt

# Translated DNA-protein alignment
diamond blastx -q reads.fna -d reference -o results.txt

# View database info
diamond dbinfo -d reference

# Custom output format
diamond blastp -q query.faa -d reference -o results.txt \
    -f 6 qseqid sseqid pident length evalue bitscore
```

## Library Usage

```rust
use diamond::prelude::*;

// Parse FASTA sequences
let records = diamond::data::fasta::read_fasta_amino_acid(
    b">query\nMKTAYIAKQRQISFVKSHFSRQLE\n" as &[u8]
).unwrap();

// Create a scoring matrix
let score_matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();

// Run Smith-Waterman alignment
let query = &records[0].sequence;
let result = diamond::dp::smith_waterman::smith_waterman(query, query, &score_matrix);
println!("Self-alignment score: {}", result.score);
println!("Identity: {}/{}", result.identities, result.length);

// Ungapped extension
let diag = diamond::dp::ungapped::xdrop_ungapped(
    query, query, 5, 5, 12, &score_matrix
);
println!("Ungapped score: {}", diag.score);

// Seed extraction and matching
let reduction = diamond::basic::reduction::Reduction::default_reduction();
let shape = diamond::basic::shape::Shape::from_code("111111", &reduction);
let seeds = diamond::search::seed_match::extract_seeds(query, &shape, &reduction);
println!("Seeds extracted: {}", seeds.len());

// E-value calculation
let evalue = score_matrix.evalue(result.score, query.len() as u32, query.len() as u32);
let bitscore = score_matrix.bitscore(result.score as f64);
println!("E-value: {:.2e}, Bit score: {:.1}", evalue, bitscore);
```

## Benchmarks

Original benchmark baseline: vendored upstream DIAMOND from `https://github.com/bbuchfink/diamond.git`, commit `1d162b4fefb5` (`v2.1.24-2-g1d162b4f-dirty`).

### Default adaptive mode (`--memory-limit 16G`)

Native Rust `blastp` and `blastx` default to a **16 GB soft process-RSS
ceiling**. Below that ceiling, retained hits stay in memory for speed. When the
ceiling is reached, all hits retained so far are migrated to DIAMOND's
compressed temporary-bin format, and every later hit remains disk-backed. In
other words, the 16 GB default still spills; it is not a requirement to have
16 GB available. Input data, indices, an active batch, allocator bookkeeping,
and one active disk bin also consume memory, so this is a spill trigger rather
than an OS-enforced hard cap and peak RSS can overshoot it.

This default deliberately differs from upstream DIAMOND `blastp`, whose hit
buffer is always disk-backed. The value follows upstream's 16 GB clustering
memory default. Use `--memory-limit 0G` to start directly in disk mode for the
closest memory-policy comparison with upstream, or set another nonzero limit
to spill earlier or later.

The following table uses the same upstream baseline and pinned CPUs as the
forced-disk table below. Rust was invoked without `--memory-limit`, exercising
the actual 16 GB CLI default. Times and peak resident set sizes are medians;
the speed ratio is C++ time divided by Rust time, while the RSS ratio is Rust
divided by C++.

| Dataset | Threads | Runs | C++ time | Rust time | Speed ratio | C++ RSS | Rust RSS | RSS ratio | Byte parity |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| Human proteins, 2,000 ref / 1,000 query | 1 | 5 | 0.88 s | 1.24 s | 0.710x | 33,648 KiB | 22,216 KiB | 0.660x | PASS |
| Human proteins, 2,000 ref / 1,000 query | 4 | 5 | 0.60 s | 0.46 s | 1.304x | 37,512 KiB | 31,388 KiB | 0.837x | PASS |
| AMRFinderPlus, 1,000 ref / 500 query | 1 | 5 | 1.60 s | 1.54 s | 1.039x | 24,452 KiB | 15,900 KiB | 0.650x | PASS |
| AMRFinderPlus, 1,000 ref / 500 query | 4 | 5 | 0.55 s | 0.47 s | 1.170x | 30,440 KiB | 18,876 KiB | 0.620x | PASS |
| AMRFinderPlus stress, 4,000 ref / 2,000 query | 1 | 5 | 18.50 s | 16.39 s | 1.129x | 54,028 KiB | 105,468 KiB | 1.952x | PASS |
| AMRFinderPlus stress, 4,000 ref / 2,000 query | 4 | 7 | 4.97 s | 4.62 s | 1.076x | 68,176 KiB | 111,060 KiB | 1.629x | PASS |

The human workload uses 1,120,714 reference and 546,969 query residues. Its
reference, query, and output SHA-256 hashes are respectively
`b10d71548bf9be8318c40ed3a47b4ee5c1628b58bd22002b9bd69bf2e140c6bd`,
`eed1277c3133ce2c88c35c42f32e392d2fccd2528a90fe0ada114f088217c244`,
and `da8d5858689815795c9121cc45994e446d8f92696c455d6e9357f72b19374172`.

The independent AMR workload uses 346,146 reference and 170,358 query
residues. Its source is the NCBI AMRFinderPlus `2026-03-24.1` `AMRProt.fa`
snapshot (SHA-256
`6b5b02061f2a3132e516f951f296ae766381e875b09dbe12867c1c4693646be4`).
Its reference, query, and byte-identical output SHA-256 hashes are respectively
`a59a1d3f96d14fe7ad814f4635b65d804a16e6d2e1a43c29288bb905dde04da0`,
`fa0c0f05486af5ad5849c8f52ebe455a2c1c6571a8497a930782b4be96c729e9`,
and `51847a3bde8172949ad6d934af235a3163a040a51c681bbe476f8e92f126eb71`.

### Upstream-compatible forced-disk mode (`--memory-limit 0G`)

Upstream DIAMOND's normal `blastp` hit buffer is disk-backed; it does not call
the external BLASTP program. Upstream v2.1.24 rejects `--memory-limit` for this
workflow, so the apples-to-apples commands below use normal upstream `blastp`
and Rust `blastp --memory-limit 0G`. Zero enters disk mode before stage 1.

These are medians of alternating-order runs on independent, non-synthetic
protein datasets. Both implementations were built for the host CPU. The C++
build uses `-march=native`; the Rust benchmark build uses
`-C target-cpu=native` but keeps general code at AVX2 on x86. Rust then
runtime-selects its explicitly specialized AVX-512 search workers. This keeps
the long score/trace loops at the same AVX2 tier as upstream instead of letting
LLVM's extended EVEX register allocation lower their sustained clock. Rust's
release profile uses ThinLTO with one codegen unit. The same pinned CPUs were
used for each C++/Rust pair.

| Dataset | Threads | Runs | C++ time | Rust time | Speed ratio | C++ RSS | Rust RSS | RSS ratio | Byte parity |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| Human proteins, 2,000 ref / 1,000 query | 1 | 5 | 0.88 s | 1.34 s | 0.657x | 33,648 KiB | 22,212 KiB | 0.660x | PASS |
| Human proteins, 2,000 ref / 1,000 query | 4 | 5 | 0.60 s | 0.60 s | 1.000x | 37,512 KiB | 30,312 KiB | 0.808x | PASS |
| AMRFinderPlus, 1,000 ref / 500 query | 1 | 5 | 1.60 s | 1.76 s | 0.909x | 24,452 KiB | 12,356 KiB | 0.505x | PASS |
| AMRFinderPlus, 1,000 ref / 500 query | 4 | 5 | 0.55 s | 0.56 s | 0.982x | 30,440 KiB | 19,172 KiB | 0.630x | PASS |
| AMRFinderPlus stress, 4,000 ref / 2,000 query | 1 | 5 | 18.50 s | 18.77 s | 0.986x | 54,028 KiB | 42,836 KiB | 0.793x | PASS |
| AMRFinderPlus stress, 4,000 ref / 2,000 query | 4 | 7 | 4.97 s | 5.08 s | 0.978x | 68,176 KiB | 53,728 KiB | 0.788x | PASS |

The stress case evaluates 81.7 million raw seed pairs without materializing
the cross product. Both implementations emit 40,297 lines with SHA-256
`6ea16d64c71f833970aa8f045adbf0cd22a35c2706c69cebbc883f0643911653`.
Rust uses less RSS in every row. It reaches wall-time parity on the four-thread
human fixture and is within 2.2% on both four-thread AMR fixtures. On the
heavier stress fixture Rust is 1.4% slower at one thread and 2.2% slower at
four threads, while using 20.7% and 21.2% less RSS respectively. Fixed-cost
latency is still visible on the smaller one-thread inputs. Direct paired
hardware-counter runs show that isolating AVX-512 to the search worker reduces
four-thread task-clock by about 2% versus allowing LLVM to use AVX-512 registers
throughout; its one-thread cost is about 0.5%. On the final stress build Rust
executes 2.6% fewer instructions and 1.7% fewer cycles than C++ at one thread.
Hoisting the Hamming ISA choice into a const-generic worker removes another
1.55% of Rust's instructions; matching upstream's final-best saturation check
in the i16 trace kernel removes 0.3% more.

Native `blastp` and `blastx` can cap the fast in-memory retained-hit path and
spill to DIAMOND's compressed temporary-bin format. Use zero to select disk
mode immediately:

```bash
diamond blastp ... --memory-limit 0G --tmpdir /path/to/fast/scratch
```

Nonzero values provide the adaptive in-memory mode and use DIAMOND's decimal
`K`, `M`, `G`, or `T` suffixes. This is a soft process-RSS ceiling, not an
OS-enforced limit: the
loaded query/database, seed arrays, allocator overhead, and one active query
bin still have to fit in RAM. Once the measured RSS approaches the ceiling,
existing hits are migrated to disk and all subsequent hits stay disk-backed.
Temporary files are removed when the search finishes, including error exits.
Small hit sets should normally remain in memory unless `0G` was requested.

Other jobs were active on the benchmark host. Alternating implementation order,
CPU affinity, and medians reduce bias, but do not eliminate it. Stress-workload
ranges were C++/Rust 17.51–21.78 s / 18.17–22.87 s at one thread and
4.94–5.03 s / 4.97–5.18 s at four threads. Raw spread should be considered
alongside the medians, especially where an external job appeared mid-run.

The benchmark driver and full fixture details are in
`scripts/compare_real_cpp_rust.sh` and
`translation/optimization_benchmark.md`.

### Bundled smoke benchmark

The table below uses bundled test data, measured with the same native-CPU
release build, one thread, and `/usr/bin/time` (median of 5 runs):

| Operation | Input | C++ Original | Rust (Native) | Speedup | Peak RSS Ratio (Rust/C++) |
|-----------|-------|--------------|---------------|---------|---------------------------|
| `blastp` | `5.faa` (389 queries) vs `data.dmnd` | 0.59s, 20.8 MiB | 1.00s, 11.2 MiB | 0.590x | **0.540x** |
| `makedb` | `data.faa` (389 sequences) | 0.03s, 12.2 MiB | 0.02s, 4.4 MiB | **1.500x** | **0.359x** |

Speedup is C++ wall time divided by Rust wall time. Peak RSS ratio is Rust peak resident memory divided by C++ peak resident memory; lower is better.

Reproduce and refresh the comparison with:

```bash
scripts/compare_cpp_rust.sh
```

The harness builds both release binaries, uses the bundled 389-sequence
protein dataset, reports median wall time and peak RSS for `makedb` and
`blastp`, and compares a stable 12-column blastp output byte-for-byte. It exits
nonzero and retains a unified diff when parity fails. Quick smoke runs can use
`REPETITIONS=1`; `REFERENCE_FASTA`, `QUERY_FASTA`, `THREADS`, `RUST_BIN`, and
`CPP_BIN` can be overridden for later scaling experiments. To benchmark
existing equivalent databases without rebuilding, set both `CPP_DB` and
`RUST_DB` to their database paths (with or without `.dmnd`). Set
`RUST_MEMORY_LIMIT=0G` for the upstream-compatible forced-disk comparison.

The broader matrix exercises both the 16 GB adaptive default and immediate
disk mode, a 64 MB mid-run spill threshold, gzip input, sensitivity presets,
masking and composition-based statistics, a non-default score matrix,
representative identity/coverage/top-score filters, an extended tabular output
schema, and native `blastx`. The translated-search cases additionally cover
very/ultra-sensitive search, CBS and masking disabled, each strand separately,
genetic code 11, a fixed minimum ORF, and translated-query output fields. It
runs every mode at one and four threads,
alternates C++/Rust execution order, checks output bytes, and writes one TSV
row per comparison with time, RSS, and both ratios:

```bash
REFERENCE_FASTA=/path/to/real-reference.faa \
PROTEIN_QUERY_FASTA=/path/to/disjoint-real-query.faa \
NUCLEOTIDE_QUERY_FASTA=/path/to/real-genomic-region.fna \
BENCH_CPUSET_1=4 BENCH_CPUSET_4=4-7 \
  scripts/compare_mode_matrix.sh
```

`MATRIX_MODES`, `MATRIX_THREADS`, `REPETITIONS`, and `RESULT_DIR` can narrow a
diagnostic run. The default is five repetitions; use at least that many for a
published comparison, particularly on a shared host. `SEARCH_ARGS`,
`OUTFMT_FIELDS`, and `COMMAND=blastp|blastx` expose the same controls in the
single-comparison driver. The nucleotide input should be a natural sequence
set or contiguous genomic region—repeating a small fixture makes timing longer
without adding realistic seeds, alignments, or failure modes.

### Search-mode parity matrix

The following is the five-run matrix on the bundled 389-record protein corpus
(97,484 residues on each side), pinned to CPU 10 or CPUs 10–13. This compact
fixture is intended to exercise modes and exact output, not to represent
large-dataset throughput. `Speed` is C++ time / Rust time and `RSS` is Rust /
C++; higher speed and lower RSS are better. Every row is byte-identical.

| Mode | Threads | C++ time | Rust time | Speed | C++ RSS | Rust RSS | RSS ratio | Parity |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| Default 16G | 1 | 0.33 s | 0.45 s | 0.733x | 23,352 KiB | 12,484 KiB | 0.535x | PASS |
| Default 16G | 4 | 0.36 s | 0.35 s | 1.029x | 28,940 KiB | 19,412 KiB | 0.671x | PASS |
| 64M ceiling | 1 | 0.31 s | 0.42 s | 0.738x | 23,292 KiB | 12,484 KiB | 0.536x | PASS |
| 64M ceiling | 4 | 0.25 s | 0.16 s | 1.563x | 28,816 KiB | 16,768 KiB | 0.582x | PASS |
| Forced disk, 0G | 1 | 0.34 s | 0.45 s | 0.756x | 23,348 KiB | 12,524 KiB | 0.536x | PASS |
| Forced disk, 0G | 4 | 0.44 s | 0.21 s | 2.095x | 28,620 KiB | 18,392 KiB | 0.643x | PASS |
| gzip input | 1 | 0.33 s | 0.44 s | 0.750x | 23,420 KiB | 12,524 KiB | 0.535x | PASS |
| gzip input | 4 | 0.31 s | 0.18 s | 1.722x | 28,884 KiB | 16,156 KiB | 0.559x | PASS |
| `--faster` | 1 | 0.22 s | 0.22 s | 1.000x | 18,868 KiB | 8,000 KiB | 0.424x | PASS |
| `--faster` | 4 | 0.28 s | 0.10 s | 2.800x | 28,936 KiB | 8,616 KiB | 0.298x | PASS |
| `--more-sensitive` | 1 | 1.61 s | 2.17 s | 0.742x | 29,484 KiB | 16,592 KiB | 0.563x | PASS |
| `--more-sensitive` | 4 | 0.94 s | 0.82 s | 1.146x | 37,540 KiB | 20,656 KiB | 0.550x | PASS |
| `--ultra-sensitive` | 1 | 3.30 s | 5.40 s | 0.611x | 32,976 KiB | 19,400 KiB | 0.588x | PASS |
| `--ultra-sensitive` | 4 | 1.70 s | 2.30 s | 0.739x | 51,860 KiB | 31,020 KiB | 0.598x | PASS |
| CBS disabled | 1 | 0.27 s | 0.38 s | 0.711x | 21,768 KiB | 12,364 KiB | 0.568x | PASS |
| CBS disabled | 4 | 0.19 s | 0.15 s | 1.267x | 28,940 KiB | 15,716 KiB | 0.543x | PASS |
| Masking disabled | 1 | 0.27 s | 0.37 s | 0.730x | 23,004 KiB | 12,140 KiB | 0.528x | PASS |
| Masking disabled | 4 | 0.21 s | 0.14 s | 1.500x | 27,312 KiB | 14,292 KiB | 0.523x | PASS |
| BLOSUM45 | 1 | 0.28 s | 0.38 s | 0.737x | 22,392 KiB | 12,696 KiB | 0.567x | PASS |
| BLOSUM45 | 4 | 0.22 s | 0.16 s | 1.375x | 28,688 KiB | 15,308 KiB | 0.534x | PASS |
| Identity 70% | 1 | 0.28 s | 0.40 s | 0.700x | 22,328 KiB | 12,524 KiB | 0.561x | PASS |
| Identity 70% | 4 | 0.21 s | 0.16 s | 1.313x | 28,620 KiB | 15,744 KiB | 0.550x | PASS |
| Query/subject cover 50% | 1 | 0.27 s | 0.38 s | 0.711x | 21,972 KiB | 12,064 KiB | 0.549x | PASS |
| Query/subject cover 50% | 4 | 0.21 s | 0.16 s | 1.313x | 29,608 KiB | 13,444 KiB | 0.454x | PASS |
| Top 10% | 1 | 0.27 s | 0.36 s | 0.750x | 18,488 KiB | 7,364 KiB | 0.398x | PASS |
| Top 10% | 4 | 0.19 s | 0.14 s | 1.357x | 28,936 KiB | 8,252 KiB | 0.285x | PASS |
| Extended tabular fields | 1 | 0.28 s | 0.38 s | 0.737x | 23,012 KiB | 12,844 KiB | 0.558x | PASS |
| Extended tabular fields | 4 | 0.21 s | 0.16 s | 1.313x | 28,940 KiB | 17,688 KiB | 0.611x | PASS |

The 64 MB setting does not spill on that compact input. A separate three-run
AMRFinderPlus stress comparison (4,000 reference and 2,000 query proteins;
1,357,780 and 679,382 residues) crosses the limit during stage 1. Each run
logged migration of existing in-memory hits, retained all later hits on disk,
and read 21.9 MiB of compressed spill data. All 40,297 output rows remained
byte-identical (SHA-256 `6ea16d64c71f833970aa8f045adbf0cd22a35c2706c69cebbc883f0643911653`).

| Adaptive mid-run spill | Threads | Runs | C++ time | Rust time | Speed | C++ RSS | Rust RSS | RSS ratio | Parity |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| `--memory-limit 64M` | 1 | 3 | 20.10 s | 19.19 s | 1.047x | 54,024 KiB | 62,592 KiB | 1.159x | PASS |
| `--memory-limit 64M` | 4 | 3 | 5.12 s | 5.85 s | 0.875x | 65,824 KiB | 66,996 KiB | 1.018x | PASS |

### Native blastx comparison

The translated-search comparison uses one natural contiguous 250,000 nt
region from *E. coli* K-12 MG1655 against the independent 4,000-protein AMR
reference above, with `--sensitive`. It is not made longer by repeating reads.
These are medians of five alternating runs. All four
cases emit the same 25 rows with SHA-256
`05b35eb5a1a6f89425ae58fe120c12eff0dbe19c3ff9c3d1944bc7e339124bb6`.

| Rust memory mode | Threads | Runs | C++ time | Rust time | Speed | C++ RSS | Rust RSS | RSS ratio | Parity |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| Default 16G | 1 | 5 | 5.22 s | 2.09 s | 2.498x | 77,192 KiB | 80,204 KiB | 1.039x | PASS |
| Default 16G | 4 | 5 | 2.51 s | 0.99 s | 2.535x | 80,308 KiB | 81,332 KiB | 1.013x | PASS |
| Forced disk, 0G | 1 | 5 | 5.17 s | 1.74 s | 2.971x | 77,308 KiB | 80,776 KiB | 1.045x | PASS |
| Forced disk, 0G | 4 | 5 | 2.61 s | 0.98 s | 2.663x | 80,756 KiB | 81,536 KiB | 1.010x | PASS |

The blastx path keeps each DNA query's six translated contexts together, then
culls globally and traces only the retained targets, matching upstream rather
than applying `-k` independently to every frame. Seed arrays are built,
sorted, and joined once per shape; their parallel fill uses exact per-worker
partition ranges instead of one contended atomic increment per seed. For long
translated queries, Tantan masking reuses one workspace instead of retaining
one large workspace per Rayon worker. Rust is substantially faster on this
short fixture; its peak RSS is within 1--5% of C++. The `0G` rows exercise the same
compressed temporary-bin spill path used after the default 16 GB soft limit
is reached.

As a separate multi-query coverage check, a natural 42-contig bacterial
assembly (2,476,164 nt; SHA-256
`b1be888ba95d0d9daaa6844063b8b1c1b3074d8a133e99b5e208e34b3fe061ce`)
was searched against the same 4,000-protein database with `--very-sensitive`.
These are three-run alternating medians. The host had other active jobs, so
the implementations alternate first-run order and the counter check below is
CPU-pinned.

| Rust memory mode | Threads | Runs | C++ time | Rust time | Speed ratio | C++ RSS | Rust RSS | RSS ratio | Byte parity |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| Default 16G | 1 | 3 | 8.93 s | 8.70 s | 1.026x | 347,680 KiB | 170,532 KiB | 0.490x | PASS (531 rows) |
| Default 16G | 4 | 3 | 6.63 s | 3.44 s | 1.927x | 359,140 KiB | 226,888 KiB | 0.632x | PASS (531 rows) |
| Forced disk, 0G | 1 | 3 | 8.95 s | 8.84 s | 1.012x | 347,872 KiB | 170,184 KiB | 0.489x | PASS (531 rows) |
| Forced disk, 0G | 4 | 3 | 6.67 s | 3.42 s | 1.950x | 363,280 KiB | 236,892 KiB | 0.652x | PASS (531 rows) |

The very-sensitive output is byte-identical (SHA-256
`48b11ee5f4932930510d21d46c1843026fa4e6322a330721fa434c987e9d04f2`).
The four-thread Rust path is 1.93–1.95x faster and uses 35–37% less RSS on this
larger workload. The one-thread path is 1–3% faster while using about 51% less
RSS. Compact i16 traceback uses four-bit nibbles for batches of one to six
targets and four byte-wide bit planes for batches of seven or eight. The
layout choice is hoisted above the DP loop with const generics; the bit-plane
case removes BMI2 deposit work only where both layouts occupy the same four
bytes per cell. Score, horizontal-gap, and row-mask buffers retain capacity per
alignment worker, matching upstream's thread-local `MemBuffer` lifetime, while
the much larger traceback allocation remains call-local. Translated workers
also retain their flattened-hit and six Hauser-CBS buffers across source
queries. Stage 1 reuses worker-local match batches and compact sequence-offset
locators, and updates the raw-match atomic only once per partition. The AVX2
i8 and i16 score recurrences process two cells per loop iteration, reducing
the hot path's branch count without changing its operation order. Forced-disk
translated search coalesces adjacent decoded bins into a
bounded 64 MiB extension window. Hits are still written to and read from the
compressed spill files; the window restores worker utilization without
duplicating retained hit records.

A final pinned one-thread hardware-counter run on this fixture retired 68.68
billion instructions and 27.18 billion cycles. The comparable upstream run
retired 72.15 billion instructions and 25.73 billion cycles; the repeated
alternating wall-clock results above therefore remain the primary end-to-end
measurement. Repeated validation covered default 16G and forced-disk 0G at one
and four threads on both natural blastx fixtures; all 32 outputs matched the
hashes above. Wall-time measurements were repeated because this host was
shared. An upstream-dense traceback experiment raised Rust RSS to about 352
MiB without improving time, so the compact hybrid layout was retained.
One-run post-change parity checks also passed with CBS disabled, masking
disabled, plus-only and minus-only strands, genetic code 11, and the extended
translated-query output schema.

For optimization work, use the larger real-sequence harness rather than
duplicating the bundled records:

```bash
REAL_PROTEIN_FASTA=/path/to/reviewed-proteins.fasta \
  scripts/compare_real_cpp_rust.sh
```

It deterministically selects disjoint reference/query records, defaults to
2,000 references and 1,000 queries. Use four threads for throughput comparison
or one thread for a hardware-counter profiling run. See
`translation/optimization_benchmark.md` for fixture hashes and calibration.

The bundled quick fixture, human real-sequence fixture, and both AMR fixtures
currently have byte-identical C++ and Rust output. Every harness run checks
parity, exits nonzero on a mismatch, and retains its raw outputs and unified
diff. When built on a non-Windows target with `--features ffi`, the `--legacy`
flag falls back to C++ FFI for conformance testing.

## Architecture

```
src/
  basic/      - Core types: Letter, Sequence, Seed, Shape, Reduction
  stats/      - Scoring matrices (BLOSUM/PAM), E-value computation
  data/       - File formats: FASTA, DMND database, DAA archive
  dp/         - Dynamic programming: ungapped x-drop, Smith-Waterman, banded DP, SIMD (SSE4.1/AVX2)
  masking/    - Tantan repeat masking
  search/     - Seed extraction, hash join, hit buffer
  align/      - HSP, Match, target culling
  output/     - Output formats: tabular, pairwise, XML, SAM, PAF
  commands/   - CLI commands (dbinfo, makedb, blastp)
  cluster/    - Clustering types
  config.rs   - CLI definition with clap
```

### SIMD Support

The DP kernels and search workers use `std::arch` intrinsics with runtime detection:
- **SSE4.1**: 16-way parallel ungapped scoring
- **AVX2**: 32-way parallel ungapped scoring
- **AVX-512BW**: Specialized seed-search workers when the host supports them
- **Scalar fallback**: Works on all platforms

## Citation

When using DIAMOND in published research, please cite:

> Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life
> scale using DIAMOND", *Nature Methods* **18**, 366-368 (2021).
> [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)

If you use our translation, we recommend that you also cite the precise version you use. If you link to [crates.io](http://crates.io), you can cite the version number;
but if you link to our Git repository, for reproducibility, it is better that you provide the URL to the repository and the git hash (Github lists it high up on the page as 7 letters, under the Code button, e.g. '21751cd')

In addition, we appreciate if you cite the paper below describing the translation approach. If for some reason you struggle with journal citation limits, please prioritizing citing the original software over our translation paper.

> Johan Henriksson. Static analysis-guided agentic AI translation enables Rust as a full stack bioinformatics language. arXiv:2608.13029, 2026. https://doi.org/10.48550/arXiv.2608.13029

## License

Apache-2.0
