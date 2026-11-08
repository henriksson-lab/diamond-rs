# Standard benchmark evidence

The pinned baseline is DIAMOND 2.1.24 build 178, source commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, with AVX2 enabled and EXTRA,
WITH_DNA, WITH_MCL, and WITH_AVX512 disabled. `probe.json` records hashes.

`cpp-default.stdout` and `cpp-default.stderr` are from the unmodified pinned
executable running `diamond benchmark` at its full iteration counts. Exit
status is **1**: after the diagonal phase, the banded target iterator rejects
the absent adjusted matrix. Timing numbers are observations, not golden
performance requirements. The warning depends on the execution date.
`invalid-options.json` contains unmodified executable diagnostics.

The kernel values in `probe.json` have a different evidence scope: a copied
`benchmark.cpp` was instrumented to print otherwise unreported results and
reduce loop counts to one. Only its AVX2 benchmark object was replaced when
linking against available source-build objects. This instrumented executable
is **not** the immutable oracle, and those existing objects do not reproduce
the immutable executable byte for byte. The native reduced-workload tests
assert these diagnostic values; full command comparison independently uses
the immutable executable output.

The Rust implementation executes actual native kernels. It omits the C++
benchmark's writes through destroyed bin-zero entries and out-of-bounds
`scores[i & 128]` read; neither has a defined observable result. This does not
establish optional EXTRA benchmark coverage, heterogeneous anchored-target
parity, score overflow parity, or equivalent performance.
