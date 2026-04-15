# diamond-rs

Rust port of the [DIAMOND](https://github.com/bbuchfink/diamond) protein sequence aligner.

DIAMOND is a high-performance sequence aligner for protein and translated DNA searches, designed for big sequence data analysis. This crate provides both a CLI binary and a library API.



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

- **CLI**: `blastp` and `blastx` run natively in Rust by default; `--legacy` for C++ FFI
- **Native Rust commands**: `blastp`, `blastx`, `makedb`, `dbinfo`, `getseq`, `version`, `help`
- **Parallel**: Seed search uses rayon for multi-threaded processing
- **SIMD**: SSE4.1/AVX2 vectorized ungapped scoring
- **Library API**: Core types, scoring matrices, DP kernels, FASTA parsing, and seed search
- **Tests**: 206 tests including all 20 C++ regression tests + native-vs-FFI equivalence

## Building

### Prerequisites

- Rust 1.70+
- CMake 2.6+
- C++ compiler (GCC or Clang)
- zlib, SQLite3, pthreads

On Ubuntu/Debian:
```bash
sudo apt-get install g++ cmake zlib1g-dev libsqlite3-dev
```

### Build

```bash
cargo build --release

# With native CPU optimizations (recommended for benchmarks)
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

### Test

```bash
# Run all tests (single-threaded for FFI safety)
cargo test --release -- --test-threads=1

# Run just the 20 regression tests
cargo run --release -- test
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

Measured with `RUSTFLAGS="-C target-cpu=native"` on test data (single-threaded, 389 queries vs 1 reference, 3 runs each):

| Operation | C++ Original | Rust (Native) | Speedup |
|-----------|-------------|---------------|---------|
| `blastp` (389 queries) | 0.14s | 0.06s | **2.3x** |
| `makedb` | 0.02s | <0.01s | **>2x** |
| `test` (20 regressions) | 28.7s | 28.8s (FFI) | 1.0x |

The native Rust pipeline produces matching scores, coordinates, and bit scores. The `--legacy` flag falls back to C++ FFI for full bit-exact compatibility.

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

The DP kernels use `std::arch` intrinsics with runtime detection:
- **SSE4.1**: 16-way parallel ungapped scoring
- **AVX2**: 32-way parallel ungapped scoring
- **Scalar fallback**: Works on all platforms

## Citation

When using DIAMOND in published research, please cite:

> Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life
> scale using DIAMOND", *Nature Methods* **18**, 366-368 (2021).
> [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)

## License

Apache-2.0
