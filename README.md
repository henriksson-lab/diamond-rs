# diamond-rs

Rust port of the [DIAMOND](https://github.com/bbuchfink/diamond) protein sequence aligner.

DIAMOND is a high-performance sequence aligner for protein and translated DNA searches, designed for big sequence data analysis. This crate provides both a CLI binary and a library API.

This is a translation of the original code and not the authoritative implementation. This code should generate bitwise
equal output to the original. Please report any deviations

The aim of this project is to increase performance, especially by providing this code through a type-safe library interface.
The code can also be compiled to be used for webassembly.


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
