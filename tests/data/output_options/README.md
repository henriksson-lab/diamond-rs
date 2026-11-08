# Native search output controls

Oracle: unchanged DIAMOND 2.1.24 build 178, source commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, executable SHA-256
`d8d08608b96706db62e44036fc979f9a658e72f59fc19d8727546daa616c42d8`.
The upstream tree and binary were neither edited nor rebuilt.

`cases.json` records the exact source arguments, invocations, status and hashes
for 174 successful runs. `cases.tsv` drives self-contained native tests at one
and four threads, in-memory and forced-disk hit buffers: 348 comparisons.
The inputs reuse `../tabular_context/` and `../search_blocks/taxonomy/`.
The tests build databases using native Rust, with no C++ dependency.

Coverage includes protein/translated queries; direct/multiple reference blocks;
all standard text/binary search formats; tabular/JSON absent, bare, zero, simple,
verbose and repeated headers; dynamic taxonomy fields; plain/gzip output;
file-buffer sizes 4, 7, 4096, 65536 and 2097152; source numeric-prefix parsing;
automatic gzip/DAA suffixes and --no-auto-append. Non-tabular formats ignore
header values, including invalid ones, as source does. Bare/repeated --header
values preserve source vector-option semantics. Version headers use 2.1.24.

All decoded payloads match. For invocation-bearing headers, only the invocation
field is replaced with the exact native argv. XML database paths are replaced
with the actual native database path. 252 metadata-free CLI comparisons
also compare the original complete wire bytes, including gzip headers/trailers.
Three cases repeat through the in-process dispatcher.

`api-source.json` and `api-invocations.tsv` retain four additional source runs.
Eight supplied-query library comparisons (protein/translated, tabular/SAM,
one/four threads) supply the identical invocation string explicitly and compare
complete original gzip bytes with verbose headers; no metadata is substituted.
The typed API treats its output path as explicit; automatic suffixes are CLI policy.

Twenty-three malformed/unsupported option checks (24 on 32-bit targets) cover statuses, diagnostics and
empty-file creation side effects. `errors-source.json`,
`buffer-errors-source.json` and `compression-errors-probe.json` retain source
observations. Source builds without Zstd reject `--compress zstd` after creating
an empty file; native search preserves that baseline despite the separate
low-level Zstd codec. Even `--compress 0` is incompatible with DAA. Buffers 0--3
fail the source input peek. Negative and impossible positive sizes produce the source allocation-error
status/message through checked conversion and a fallible capacity preflight.
`oversized-buffer-source.json` preserves bounded source observations and their
resource limits. The preflight does not guarantee later allocations under
changing memory pressure. The typed search API validates impossible capacities too.

`large-source.json` records a separate source search over 1200 generated
2000-residue sequences, requesting full_qseq. The deterministic RNG parameters
and output hashes are retained. `large.gzip` contains 2,401,200 decoded bytes
and 1,456,039 compressed bytes. The test regenerates and verifies the decoded
payload and compares complete gzip streams at 65536/67108864 buffer sizes.
The compressed data exceeds the native bounded 1 MiB codec scratch space.

Production search connects the translated OutputFile/Serializer/compressor
lifecycle. The configured buffer size also reaches FASTA/FASTQ/gzip inputs,
DMND data/index readers, FASTA reference database construction and intermediate
block files. Independent buffered data/index handles retain read windows across
alternating random accesses. Codec scratch allocation is bounded and lazy;
serializer flushing forwards directly instead of copying into a second buffer.

`flate2` now uses stock C zlib through libz-sys: zlib-rs produced different
compressed bytes for identical decoded data. This is a compression dependency,
not the DIAMOND C++ adapter, and requires no adjacent source tree or CMake.
Gzip OS metadata is pinned to Unix (3), matching this source oracle on every
native target. These fixtures establish bounded Linux-source compatibility;
new macOS package jobs still require hosted execution. They do not establish
whole-program parity, performance parity, stdout terminal behavior, view-command
compression, optional upstream build modes, or completion of the translation.
