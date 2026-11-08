# Native binary edge and null search fixtures

Oracle: unchanged DIAMOND 2.1.24 build 178, source commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, executable SHA-256
`d8d08608b96706db62e44036fc979f9a658e72f59fc19d8727546daa616c42d8`.
The C++ tree/executable were neither edited nor rebuilt for this checkpoint.

The 84 cases in `cases.tsv` use `../tabular_context/` inputs. They cover both
query types, indexing algorithms, direct/multiple reference blocks, default/0/1
unaligned settings, protein CBS 0–5, translated CBS 0/1, TANTAN/SEG masking,
identity/coverage/top-score/bit-score filters and multiple HSPs. Each runs at
one thread with in-memory hit buffers and four threads with forced-disk buffers:
168 exact file comparisons. Three also repeat through `diamond::run`.

Twelve precision cases use 60 deterministically generated variant sequences
and their DNA encodings, including reverse frames. Eight cases check 180 binary
records apiece across query types, algorithms and block sizes. Four more check
120 records apiece with equal 50% query/subject coverage, exercising length
sorting and original sequence IDs. They also caught a translated reverse-pass
coverage bug: the adapter divided amino-acid length by DNA length instead of
using the source-coordinate interval. Both native reverse-pass adapters now
use `query_cover_percent`, as the source does. Each runs at one/four threads, yielding 24
file comparisons and 3,840 compared binary records. Separate supplied protein
and translated-query APIs add four exact edge comparisons. Thus the new checks
cover 196 source-output comparisons plus three dispatcher repeats and two
malformed-FASTQ diagnostics. The null cases require observed seed-search work;
malformed input also proves null does not skip input processing.

All runtime tests build their own databases in Rust and need no C++ checkout.
The source manifests preserve commands/statuses/hashes; `provenance.json`
records input/output identities. Raw binary files are compared without float
rounding or tolerance. The 32-byte little-endian record layout in this oracle
is two u64 original IDs, two f32 coverages, then one f64 corrected bit score
(the source field is misleadingly named `evalue`). The implementation uses
native-endian writes as source does; supported CI platforms are little-endian.

Numerical compatibility: the pinned executable was compiled with
`-march=native` and contracts multiply/add expressions in corrected bit scores
and the ALP log-area calculation. Rust now uses explicit `mul_add` at those
operations to reproduce that rounding. This fixes previously hidden ULP
mismatches exposed by binary output; decimal tabular fields masked them.
`AlignmentEvaluer::log_area` also passes subject/query dimensions in that order,
which the native score-matrix adapter now preserves. This is evidence against
the identified oracle build, not a claim that every C++ compiler/build emits
identical floating-point bytes.

`null` retains the source's TRANSCRIPT DP requirements and inherited no-op
callbacks. `edge` requests COORDS, does not need query sequence copies, and
serializes actual HSP values/original IDs after block merging. Neither emits
unaligned rows or headers. Advanced search controls and taxonomy format 102
remain separate translation/integration work. Whole-program console text,
large-scale performance and all numerical inputs are outside this bounded
fixture matrix.
