# Optimization benchmark

The optimization workload uses distinct, unmodified records from a real
protein collection. It does not repeat or synthesize sequences. The fixture
preparer interleaves two disjoint samples from one source FASTA so reference
and query inputs cover the same portion of the source distribution without
sharing identifiers.

The calibrated default is:

- 2,000 reference proteins, 1,120,714 residues
- 1,000 query proteins, 546,969 residues
- four threads and five alternating-order measurements for throughput comparison
- one thread for a stable hardware-counter profiling run

Prepare and run it with:

```bash
REAL_PROTEIN_FASTA=/path/to/reviewed-proteins.fasta \
  scripts/compare_real_cpp_rust.sh
```

`REFERENCE_COUNT`, `QUERY_COUNT`, `THREADS`, `REPETITIONS`, and `FIXTURE_DIR`
are configurable. `SKIP_BUILD=1` avoids rebuilding during short optimization
loops. On x86-64 the harness defaults to `-C target-cpu=x86-64-v3`: a native
build on the development host introduced AVX-512 instructions into the AVX2
kernel and was slower. The wrapper prints source and derived-fixture SHA-256
hashes, sequence and residue counts, median wall time, peak RSS, timing spread,
and exact output parity. A parity failure retains the raw outputs, timing table,
and unified diff and returns a nonzero status.

## Calibrated Swiss-Prot snapshot

The calibration used 20,431 reviewed human Swiss-Prot records with decompressed
source SHA-256
`b205ed9970fee62413e2f192c6fb82daab86e4a5897ef521f907ecbfcd7b3e63`.
The derived fixture hashes are:

- reference: `b10d71548bf9be8318c40ed3a47b4ee5c1628b58bd22002b9bd69bf2e140c6bd`
- query: `eed1277c3133ce2c88c35c42f32e392d2fccd2528a90fe0ada114f088217c244`

After structurally retranslating the AVX2 SWIPE score and traceback kernels,
the final five-run four-thread medians were C++ 0.80 s / 36,308 KiB and Rust
0.76 s / 30,940 KiB. Rust is 1.053x faster and uses 14.8% less peak RSS. The
single-thread medians were C++ 1.55 s / 30,064 KiB and Rust 2.11 s /
23,156 KiB (0.735x speed and 0.770x RSS ratios).

Output is byte-identical, with SHA-256
`da8d5858689815795c9121cc45994e446d8f92696c455d6e9357f72b19374172`.

## Independent AMR stress fixture

An independent run uses NCBI AMRFinderPlus `2026-03-24.1` `AMRProt.fa`, not a
resample of the human source. The source SHA-256 is
`6b5b02061f2a3132e516f951f296ae766381e875b09dbe12867c1c4693646be4`.
The 1,000-reference / 500-query fixture contains 346,146 and 170,358 residues;
its hashes are:

- reference: `a59a1d3f96d14fe7ad814f4635b65d804a16e6d2e1a43c29288bb905dde04da0`
- query: `fa0c0f05486af5ad5849c8f52ebe455a2c1c6571a8497a930782b4be96c729e9`

Final alternating-order medians after the AMR parity fixes are:

| Threads | C++ time / RSS | Rust time / RSS | C++/Rust speed | Rust/C++ RSS | Parity |
|---:|---:|---:|---:|---:|:---:|
| 1 | 2.25 s / 21,812 KiB | 4.81 s / 16,372 KiB | 0.468x | 0.751x | PASS |
| 4 | 1.00 s / 29,404 KiB | 2.07 s / 23,760 KiB | 0.483x | 0.808x | PASS |

A 4,000-reference / 2,000-query diagnostic evaluates 81.7 million raw seed
pairs. The streaming Rust pipeline completed in 15.11 s / 173,948 KiB, versus
32.47 s / 4,526,948 KiB before the fix (2.15x faster, 96.2% less peak RSS).
C++ measured 5.89 s / 66,612 KiB. Both outputs match byte-for-byte: 40,297
lines with SHA-256
`6ea16d64c71f833970aa8f045adbf0cd22a35c2706c69cebbc883f0643911653`.
The Rust join now walks each matching seed group through the Hamming and
left-most filters in batches of at most 32 instead of allocating its q×r pair
cross product.

Other jobs were active on the benchmark host. The harness alternates C++/Rust
order and reports medians. Final AMR timing ranges were 2.19–2.35 s /
4.66–4.95 s at one thread and 0.95–1.07 s / 2.05–2.48 s at four threads.

## Profiling

Generate repeated hardware counters and flat profiles for both implementations
with:

```bash
REAL_PROTEIN_FASTA=/path/to/reviewed-proteins.fasta \
  scripts/profile_real_cpp_rust.sh
```

The default uses one thread, three `perf stat` repetitions, and one sampled
`perf record` per implementation. It retains the counter files, `perf.data`,
flat reports, outputs, and databases under `FIXTURE_DIR/profiles`, and verifies
byte parity. `THREADS` and `STAT_REPETITIONS` are configurable.

The pre-retranslation profile localized the runtime gap to SWIPE dynamic
programming:

| Rust symbol group | Rust cycles |
|---|---:|
| AVX2 traceback | 41.69% |
| AVX2 16-bit score | 17.10% |
| AVX2 8-bit score | 16.13% |
| Seed array build/count | 12.17% |

The three SWIPE kernels total 74.92% of Rust samples. The corresponding C++
SWIPE variants total about 33.8% of its samples. In the matched three-run
counter measurement Rust executed 30.09 billion instructions in 3.624 s,
versus C++ at 6.63 billion in 1.606 s. Heap allocator symbols were approximately
0.1% in the earlier profile and below the 0.5% reporting threshold in the final
profile, so allocation is not the active speed bottleneck.

That profile led to a direct port of upstream's common moving band,
vector-resident score profile, rolling score/gap rows, last-equal row-counter
rule, and packed traceback masks. Sparse one- and two-lane i16 batches use a
four-bit-per-live-lane trace representation because a fixed 16-lane matrix
caused a measured 19.2 MiB allocation for one two-lane batch. Dense batches
retain the upstream fixed-width mask representation.

The post-retranslation five-run one-thread counter pass measured Rust at
1.728 s (SD 0.021 s) and C++ at 1.478 s (SD 0.054 s). Rust executed 13.49
billion instructions versus 6.63 billion for C++. Its remaining sampled
hotspots were i16 traceback (25.48%), seed construction (20.66%), i16 scoring
(10.15%), i8 traceback (7.42%), and i8 scoring (6.75%). The single-thread
result is still about 17% behind upstream, while the four-thread throughput
run above is 1.68x faster; future tuning should therefore start with i16
traceback and seed construction rather than allocation changes.

An earlier calibration that sampled only the leading source region exposed two
differing target selections despite equal row counts. That biased selector was
rejected for performance work in favor of full-source sampling, but the result
confirms why real distinct sequences—not cloned inputs—are needed for separate
correctness stress fixtures.
