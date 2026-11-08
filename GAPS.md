# Validation gaps

This checklist tracks native DIAMOND modes that are not yet covered by a
realistic upstream-C++ comparison of output parity, speed, and peak RSS.
Passing Rust-only unit tests does not close an item unless the item explicitly
asks only for internal coverage.

The final 2026-10-06 matrix passed all 164 cases in its defined search,
scoring, clustering, and DAA/view scope. No failure from that matrix is known
to remain. Unchecked items below are untested or untranslated surface area,
not known failures unless the text explicitly says otherwise.

## Active translation batch (2026-11-08)

- Raw OS arguments/paths and GVC/TSV byte content are integrated. The retained
  18 raw-OS, six GVC and thirteen TSV controls match Source; their broader
  combination/platform scope remains qualified.
- The latest batch implements same-pool nested DP tasks, explicit scoped query
  pools, parallel external volume readers/buffers, retained multinode streams,
  and Source's formatted manifest extraction. It passes isolated default,
  optional and Windows compilation plus 204 focused tests (one SIMD case ignored).
- The retained multinode indexed partition/load_seqs adapter is now concrete.
  The next batch adds parallel seed/pair/chunk workers and retained FASTA index
  streams; its frozen package passes three compile profiles and 142 focused tests.
  A subsequent frozen batch implements oversized radix spill parallelism,
  parallel chunk input loading/sorting and the concrete retained multinode
  search/GVC backend. It also connects actual intermediate generation and
  corrects seed-log formatting. Its frozen package passes three compile
  profiles and 112 focused tests; the final formatter-only package separately
  passes default compilation and scoped multinode/formatter checks. Fresh
  Source comparisons and broader command/profile interactions remain gates.
- Oversized POD spill retains a documented Source stack-overflow boundary:
  radix_sort.h:38 supplies sizeof(T) as a typed element count. Rust safely
  decodes one record and preserves the zero-initialized/previous-record bytes
  on EOF; Source overflow behavior is not a parity claim.
- The immutable d4fac90 release passes three natural 9,704-record controls
  at one/four threads with one cross-process worker: exact mapping, streams,
  alignment row multisets and counters. Full artifact comparison reopened
  concrete chunk destructor, representative-worker and literal path gaps.
  Those fixes plus bounded parallel blocked/accession output now pass three
  compile profiles and 101 focused tests. New release comparisons against the
  preserved Source captures verify all three code fixes. Directory modes
  remain different in c3ea396; the subsequent mkdir/log package passes three
  compile profiles and 79 focused tests but awaits new-release comparisons.
  External Edge loading/sort, assignment workers/EOF and closure compression
  now pass three compile profiles and 83 focused tests. Realign default-thread
  follow-up passes default compilation and three integration controls. The
  earlier 8965365 release closes all 17 directory modes and all three constructor
  controls; the new clustering/POD release still needs upstream comparisons.
  Broader final-snapshot parity and full artifact qualification remain open.
- Hosted CI verification, the Source-shared two-worker hazard, and final parity
  gates remain open. No complete-translation claim applies.

## Current implementation and remaining evidence (2026-11-08)

Baseline commands, enabled optional adapters, external clustering, streamed BLAST
references and shared search/output orchestration are connected in Rust. There
is no automatic production C++ fallback. Earlier adapter-missing statements
below are historical and superseded; unchecked validation items remain evidence
gates, not proof that their code is absent. Current results and explicit Source
compatibility boundaries are in `TRANSLATION_COMPLETION.md` and
`translation/progress.md`. The overall translation goal is unachieved. The user resumed translation work while hosted CI is monitored.

The latest external two-worker gate remains unresolved: a late unconditional
empty-chunk allocation can disagree with another worker's captured chunk count.
Pinned Source has the same scheduling hazard; its preserved two-worker and
one-worker runs also produce different mappings. Current storage-owner and BLAST metadata corrections pass their focused
checks. The corrective codecs match all 68 direct Source snapshots; the recorded
one-worker OID and qualified default-accession CLI controls also pass. Later
ownership/console snapshots retain their focused evidence scope. See the 2026-11-08 follow-up in `TRANSLATION_COMPLETION.md` for the exact
snapshot scopes; earlier selected passes do not close these newer findings.

Final-snapshot integration/profile/package checks, hosted platform runtime/RSS
results, unavailable historical data and natural production-scale boundaries
remain open. A successful compile, inventory count or bounded comparison does
not close those gates.

The real 1.1-billion-residue default makedb comparison found six mask-bit
differences despite successful builds and identical counts. Rust's runtime FMA
selection did not match Source's default compiler profile. The kernel correction
and strict original-record tests now pass; the corrected full default-masking
database is byte-identical to the preserved Source output (1,231,665,708 bytes).
All four crossed dbinfo calls also match. This closes that specific scale case;
the separate production taxonomy-spill case also passes as recorded below;
platform and final-snapshot runtime gates remain open. The corrected frozen snapshot also passes 1,876 default tests and all 20
release hashes; these precede the pending sorter-ownership batch. See
`TRANSLATION_COMPLETION.md` and
`.tmp/native-post-translation-20261107/scale-data/makedb-current/`.

## Optimization rotation rule

Do not optimize one command indefinitely. Work in measured rounds across four
lanes: search (`blastp`/`blastx`), clustering, database/index utilities, and
DAA/view. In each round, first require exact upstream parity on natural,
nonduplicated input whose steady-state work takes roughly 3--5 seconds or
longer. Then test at most one profile-backed hypothesis in a lane before
rotating, unless that hypothesis yields a clear faithful improvement. Record
runtime and peak RSS together, retain only measured improvements, and never
trade correctness or upstream semantics for a benchmark result.

Default external accession output now has its variable-record radix stage and
Source per-volume output-stream scope connected. The 9,704-row natural default
workflow matches all 170 filenames/bytes and diagnostics against the explicit
GNU 12 output-component Source profile. Original GNU 11's allocator-backed list
sort abort remains separately qualified. POSIX FileSink stdio write/close timing
and the external destructor policy pass their direct Source comparisons,
including actual CLI ENOSPC failure and primary-error/cleanup-error ordering
followed by SIGABRT. Input/shared-FILE controls and all eleven retained FASTA
class controls also match unchanged Source exactly. Subsequent spool lifetime,
console, filename and signal fixes pass the recorded native compile profiles
and focused checks. Windows runtime/ACP diagnostic bytes, further raw-byte title
integration and final-snapshot broader workflows remain open. These selected
controls do not establish full platform or command-combination parity.

## Raw argv and filename evidence

The executable now uses OS-native arguments and paths through the connected
parser, input/output owners, diagnostics and invocation metadata. The 18 retained
raw-OS CLI controls match Source exactly, including full DMND/DAA artifacts and
the unknown-output-format diagnostic. There is no lossy argument encoding.

Checked String library adapters remain available. The retained multinode library
and partition-file follow-ups have their own code scope and do not invalidate
or broaden those CLI captures. Windows ACP versus wide-path support and Source's
XML escape-table undefined behavior retain their explicit qualifications.

## Priority order

- [x] Test true-default and every sensitivity preset for `blastp` and
  `blastx`.
- [x] Compare `cluster`, `linclust`, and `deepclust` against upstream.
- [x] Compare every supported scoring matrix.
- [x] Exercise adaptive mid-run spilling in `blastx`.
- [x] Compare every native `view` format.
- [x] Compare `makeidx`, `merge-daa`, and `getseq` on successful nontrivial
  workflows. All three generators now produce byte-identical output; native
  persistent-index consumption is also exact on the natural scout below. The
  refreshed `merge-daa` case combines two contiguous, disjoint 1,500-query
  natural shards into a byte-identical 509 MiB archive. Three-pair medians are
  6.02/4.50 s upstream/Rust and 210,416/9,656 KiB, so that lane is closed
  without profiling.
- [x] Compare compressed natural `makedb` and restore upstream thread behavior.
  On the 89 MiB gzip Swiss-Prot input, profiling showed TANTAN—not inflate—as
  the old bottleneck. Parallel in-place masking plus the correct hardware-
  concurrency default improves exact old/new Rust medians from 17.49 to
  4.85 s. Fresh upstream/Rust medians are 6.38/4.10 s with RSS
  552,272/379,412 KiB; explicit `-p1` remains exact and serial.

## Portability

The long-term target is native support on Linux, macOS, and Windows. New
low-level code must use the target platform's actual ABI and must not expose
Linux-only declarations on other targets.

- [ ] Keep routine CI builds and tests green on Linux, macOS, and Windows.
- [x] Add target-specific compile checks for OS-gated FFI and filesystem code
  where the CI toolchains permit them. CI now covers i686 GNU/MSVC in addition
  to native OS jobs; local i686 Linux and Windows-GNU builds pass.
- [ ] Replace explicit unsupported platform branches with native
  implementations when the affected feature is part of a supported workflow.

### Native-ABI reduction follow-up

The implementation and detailed evidence are tracked in `FFI_REDUCTION_PLAN.md`.
Handwritten OS ABI declarations have been replaced with maintained portable
crates, and raw temporary-file/mmap ownership has been replaced with RAII.
The parity-sensitive `memchr` and glibc `malloc_trim` calls use target-specific
declarations maintained by the `libc` crate. C `rand`/`srand` has been replaced
with a fully specified Rust implementation of the established GNU sequence, so
the generator and `RAND_MAX` are identical on every target. The optional
default-off C++ conformance adapter is the sole remaining handwritten native
ABI boundary. It is a source-checkout-only diagnostic: the published package
intentionally excludes the upstream C++ tree. Its build script detects that
omission, so packaged `--all-features` builds succeed with the adapter
explicitly unavailable instead of looking for absent sources.
The default published configuration has now been exercised directly:
`cargo package --locked --allow-dirty` passes its verification build, and the
unpacked package passes `cargo test --locked --all-targets` with 1,742 library
tests plus all runnable binary/integration targets. The archive contains no
`diamond/` tree or C/C++ sources. `rusqlite` resolves through bundled
`libsqlite3-sys`, and
the packaged executable has no dynamic `libsqlite3` dependency. Enabling the
checkout-only `ffi` feature in the package emits a diagnostic warning and
leaves `--legacy` unavailable; packaged all-features builds do not require the
excluded tree. Ubuntu CI now repeats both the isolated `cargo package --locked`
verification and the unpacked package's complete all-target test suite, so a
new compile-time dependency on checkout-only files fails before release.
The remaining validation gaps are deliberately not approximated:

- [ ] Observe the updated native RSS tests on actual macOS ARM64 and Windows
  MSVC CI runners; cross-compilation alone does not validate RSS units.
  The current source additionally passes full linked cross-builds for macOS
  ARM64/x86-64, Linux musl, Linux i686, and Windows GNU x86-64. Windows
  standard streams now duplicate the inherited handle instead of opening
  `CONIN$`/`CONOUT$`, preserving pipes and redirected CI output. Actual hosted
  MSVC execution remains the authoritative open gate.
- [x] Add Linux musl runtime and selected 32-bit compile coverage. Local musl
  build/runtime and selected i686 Linux/Windows compile evidence are recorded
  in `FFI_REDUCTION_PLAN.md`; native CI evidence remains a separate open gate.
- [x] Replace target-dependent C `rand`/`srand` with a deterministic Rust
  generator and run raw plus downstream statistical vectors on every target.
  The retained sequence is byte-for-byte compatible with the established GNU
  vectors; divergence from platform-specific upstream C runtimes is explicit.
- [x] Validate `rusqlite` parent lookup with an external, real
  `taxonomy4blast.sqlite3` and compare translated taxonomy field formatting
  byte-for-byte with upstream output from the same NCBI artifacts.
- [x] Wire the translated DMND taxonomy metadata into native `blastp` and
  `blastx`. Search loads only the mapping/nodes/names sections required by the
  requested fields and emits exact `staxids`, names, lineage, kingdom, and
  phylum output. Direct blastp/blastx comparisons against C++- and Rust-built
  databases are byte-identical; ordinary searches retain the lean reader and
  do not allocate an empty taxid vector per subject.
  - [x] Extend the native `makedb` configuration boundary with `--taxonmap`,
    `--taxonnodes`, `--taxonnames`, and `--no-parse-seqids`; do not parse these
    as inert generic options.
  - [x] During each existing one-billion-letter FASTA block, extract every
    accession from each full title with `accession_from_title`, pair it with
    the global OID, and release the sequence block as today. Upstream does
    this only when `--taxonmap` is present and preserves the accession parsing
    counters.
  - [x] Make the translated `ExternalSorter` genuinely external before using
    it here. Upstream spills sorted runs at roughly 2 GiB and k-way merges
    them; Rust now uses the same threshold/accounting, portable temporary run
    files, the same four-byte run-order indices, a parallel unstable run sort,
    and a heap k-way merge.
  - [x] Translate `TaxonList::build` at its streaming boundaries: accept only
    the two upstream mapping headers, suppress consecutive duplicate mapping
    accessions before sorting, sorted-join accession/OID with accession/taxid,
    externally sort by OID, remove taxid 0, deduplicate taxids per sequence,
    and write one varuint count followed by sorted varuint taxids for every
    OID including unmapped OIDs. Set `taxon_array_offset` before these raw
    bytes and `taxon_array_size` afterwards; there is no extra container
    header.
  - [x] Serialize `--taxonnodes` as upstream does: `u32` parent-array length,
    one `i32` parent per taxid, then one rank byte per entry. Serialize
    `--taxonnames` as `u32` vector length followed by one NUL-terminated
    scientific-name string per taxid, including empty slots. Set the three
    `ReferenceHeader2` offsets in mapping/nodes/names order and rewrite both
    headers only after the trailer completes.
  - [x] Preserve failure cleanup and validate byte-for-byte databases for
    mapping-only, nodes-only, names-only, all-three, parsed/raw seqids,
    multiple accessions per title, unmapped records, taxid 0, and duplicate
    mapping keys. The compact fixtures are byte-identical in all five option
    combinations. A three-pair natural 150,000-record `--taxonmap` comparison
    is also exact: 4.39/4.46 seconds and 126,172/105,600 KiB upstream/Rust.
    Each implementation opens the other's database. The production threshold
    now passes too: a verified official mapping projection contains 122,563,424
    accepted rows and 2,147,483,664 bytes. Both unchanged executables produce
    identical 295,634,030-byte databases, report all four identical taxonomy
    counts, and create two actual unlinked full-bucket runs with identical
    sampled prefixes and sizes. Source/native wall time is 93.13/283.91 s and
    RSS 2,709,548/3,806,048 KiB in this single shared-host scout; no performance
    improvement is claimed. This proves the accession-to-taxid spill, not a
    matched-OID sorter spill or a scan of every temporary-file row. Evidence:
    `.tmp/native-post-translation-20261107/scale-data/taxonomy-spill-current/`.
    Failure cleanup follows Source's sequence-loop boundary; later taxonomy
    errors retain the incomplete database. Rust zeros Source's uninitialized
    initial Header1 position field, with every defined byte compared strictly.
- [ ] Verify the enabled native NCBI BLAST reference path while retaining its
    documented Source compatibility boundary. Historical diagnosis: the pinned
    upstream `load_parallel`
    path invokes `StringSetBase::append` on an empty destination; that removes
    its leading 256-byte perimeter while retaining offsets based at 256, then
    reads an uninitialised/out-of-bounds perimeter from the temporary decoded
    package. On the real two-record fixture this produces blank subject IDs,
    corrupted scores, and traceback failures despite `.psq` decoding exactly
    to the source FASTA. Rust's random-access reader produces the intended
    sequences and taxonomy, but that is not byte parity with the upstream UB,
    so native safe decoding is explicitly qualified instead of claimed as
    parity with undefined Source bytes. Native streaming/filter/taxonomy/OID
    routing is now connected; broader defined-path conformance remains open.
- [x] Find a natural `blastx` workload whose retained-hit buffer grows enough
  to cross a nonzero memory limit after the initial RSS check.
  - [x] The distinct-read SRR1703350/Swiss-Prot fixture crossed 1150 MiB after
    accumulating 526,641--547,543 retained hits and read back 5.1 MiB of spill
    data. Default, forced-disk, and delayed-spill outputs were byte-identical.

## Native commands

- [x] Compare `cluster` on realistic data: parity, speed, and peak RSS.
- [x] Compare `linclust` on realistic data: parity, speed, and peak RSS.
- [x] Compare `deepclust` on realistic data: parity, speed, and peak RSS.
  - The three commands now select distinct translated workflows. Exact compact
    parity is established for all six tested default/identity cases: `cluster`,
    `linclust`, and `deepclust`, each at its command default and at
    `--approx-id 90 --member-cover 80`. The realistic 500-protein matrix now
    also passes all 12 cases, including default at one/four threads and the
    identity-50/coverage-60 cases. Two independent post-rewrite matrices are
    18/18 exact overall; artifacts are under
    `.tmp/cluster-endpoint-structural-final/` and
    `.tmp/cluster-endpoint-structural-repeat2/`. The second run took 26.83 s
    in Rust versus 37.68 s for the retained scalar-Rust baseline, with 8,760
    versus 8,688 KiB maximum RSS. The earlier 13/18 result under
    `.tmp/ffi-cluster-hamming-parity/` remains useful diagnostic history: its
    round-2 discrepancy was fixed by faithfully preserving active OID order
    after round 1 and by translating the upstream no-trace endpoint kernels.
- [x] Compare `makeidx` output and performance with upstream. The first pass
  found a real difference because Rust's `Block::soft_mask(Motif)` did no motif
  masking. Restoring upstream's masking-table operation makes both 2.1 MiB
  indexes byte-identical (`9eb38148...`). The post-fix scouts were 0.21
  s/20,800 KiB C++ and 0.25 s/22,480 KiB Rust. A later steady-state fixture of
  150,000 distinct Swiss-Prot records produces identical 129 MiB indexes
  (`dbe81bcd...`) but initially measured 8.00 s/394,968 KiB upstream versus
  10.55 s/437,128 KiB Rust. Removing the
  apparent per-shape sequence copies changed a controlled median only from
  9.70 to 9.59 seconds with flat RSS, so that experiment was reverted and the
  profile was followed instead. Restoring upstream's threaded soft-masking
  stage and `--threads` handling improves the Rust default median by 15.3%.
  Streaming DMND records directly into the block then removes the second
  retained decoded database representation and cuts Rust RSS by 20.2% in a
  controlled A/B. A fixed-invariant hashed-enumerator specialization preserves
  both passes and exact insertion order while hoisting the production `makeidx`
  no-skip/no-filter decisions. Three exact Rust A/B pairs improve 7.79 to 7.45
  seconds at flat RSS. A later isolated p1 refresh measured
  8.93 s/394,972 KiB upstream versus 9.57 s/349,760 KiB Rust (speed 0.933x,
  RSS ratio 0.886x). Exact `perf` profiles showed that seed enumeration differs
  by only 0.85% in absolute cycles; motif masking accounted for about 57% of
  the remaining cycle gap. Rust was executing hardware division for the fixed
  `% 20^7` rolling-kmer modulus while upstream strength-reduced the same
  template constant. A const-specialized identity iterator preserves exact
  output and improves the Rust p1 median from 9.52 to 8.88 seconds at flat
  RSS, closing the fresh 8.93-second upstream gap. A fresh default-thread
  comparison is 7.45 s/394,968 KiB upstream versus 7.61 s/347,236 KiB Rust;
  the remaining 2.1% is below the material threshold and the RSS gap is closed.
  A post-streaming profile put 70.5% of Rust cycles in hashed-seed table
  construction before the specialization, matching upstream's 68.6% dominant
  phase. An unsafe raw-pointer probe experiment changed five-run medians only from 10.61 to
  10.52 s with mixed pairwise direction and flat RSS, so it was reverted.
  A separate upstream-shaped unchecked sequence load in
  `HashedSeedIterator` was also exact but changed five-run medians from
  10.49 to 10.51 s with flat RSS; it too was reverted. Do not revisit bounds
  checks without stronger evidence.
- [x] Finish exact persistent-index search parity. Native `blastp` and `blastx`
  now consume `--target-indexed -c1` through the translated hashed-seed stage 0
  and the regular Hamming/ungapped/extension pipeline. On the natural 75-query
  by 7,000-target scout, both implementations enumerate exactly 28,375 and
  28,108 query seeds and 2,286,238 and 2,265,280 reference seeds for the two
  shapes. Both then produce 1,683,263 raw pairs, 76,067 left-most survivors,
  and 1,241 byte-identical rows (SHA-256 `100681ab...`). The parity fix was a
  faithful hash-join correction: extracted bits select the initial table
  bucket, but the complete 32-bit seed remains the equality key. One release
  observation was 1.16 s/44,524 KiB C++ versus 1.06 s/67,228 KiB Rust. The
  remaining RSS gap is a separate optimization item; Rust currently owns the
  2.0 MiB seed-index read while upstream mmaps it. Artifacts are under
  `.tmp/target-indexed-check/`.
  The subsequent direct-fill translation brought forced-disk RSS to parity,
  and a bounded zero-copy adapter pass now lets stage-2 consume the original
  borrowed packed locations plus HitField indices. Five exact old/new Rust
  pairs improve median runtime from 4.83 to 4.06 s (15.9%) with effectively
  flat RSS. Fresh complete forced-disk C++/Rust medians are 3.94/4.01 s and
  45,504/46,404 KiB. The indexed-search experiment is closed while work
  rotates to other commands; artifacts are under
  `.tmp/indexed-zero-copy-agent/` and `.tmp/indexed-zero-copy-final/`.
- [x] Add a successful, nontrivial `merge-daa` comparison. The original 30 MiB
  two-input scout deliberately repeated one query set and remains only a
  correctness check. The accepted fixture instead splits 3,000 natural queries
  into two contiguous, disjoint 1,500-record shards. Both implementations
  deduplicate the per-file reference lists to 6,615 targets, append 2,894 query
  records, and produce the same 532,968,582-byte SHA-256 `7c168e79...`.
  Three-pair medians are 6.02 s/210,416 KiB upstream and 4.50 s/9,656 KiB Rust;
  cross-view output is also exact at 64,204 rows and SHA-256 `dde6d4f5...`.
  Evidence is under `.tmp/merge-daa-two-natural-20261012/`.
- [x] Compare `getseq` output, speed, and peak RSS with upstream. All 7,000
  records and selected `--seq 1 3` output are byte-identical. The selected
  scout was 0.03 s/9,740 KiB C++ versus 0.02 s/9,600 KiB Rust. This comparison
  also fixed Rust's former seqid interpretation: upstream uses 1-based numeric
  sequence ordinals and accepts a space-separated list.
  A later 574,627-record Swiss-Prot rotation found that Rust materialized the
  complete database. A faithful sequential reader now streams sequence data
  alongside the position table. Full and selected outputs are byte-identical;
  RSS falls from 600,960 to 3,200 KiB (upstream 137,740 KiB). A bounded
  formatting follow-up also restores upstream's reusable per-record
  `TextBuffer` boundary. Uncontended Rust improves from 1.69--1.72 to
  1.05--1.13 s; three fresh, exact pairs under shared-host contention give
  medians of 3.48 s upstream and 3.03 s Rust, with every pair favoring Rust.
- [x] Match standard upstream's rejection of non-`getseq` controls instead of
  silently ignoring them. The apparent `--reverse` and `--hard-masked` paths
  belong to upstream's disabled `EXTRA` option group and are not exposed by a
  normal DIAMOND 2.1.24 build. Rust now matches upstream stderr and exit status
  for those flags, `--oid-list`, `--chunk-size`, `--query`/`-q`, and
  `--block-size`/`-b`. Standard `getseq` functionality is `-d`, `-o`, and
  `--seq`, all already translated.
- [x] Stream `makedb` in upstream-sized blocks for inputs beyond one billion
  letters. Rust now mirrors `FastaFile::load_seqs(1e9)`: the record that
  reaches/crosses the limit completes the block, which is masked in place,
  written, hashed, and dropped while only the 16-byte-per-sequence position
  trailer remains resident. A small-limit regression exercises multiple
  blocks without a fabricated huge fixture and produces the same DMND bytes as
  the owned-record path. Natural 3.4-million- and 59.4-million-letter inputs are
  byte-identical to upstream; the latter measures 4.44 s/126,084 KiB upstream
  and 4.23 s/96,000 KiB Rust at one thread.
- [x] Validate peak RSS across the real one-billion-letter `makedb` block
  boundary. The checksummed UniParc 2026_03 prefix contains 2,862,293 distinct
  complete records and 1,100,000,118 residues. Corrected native default masking
  produces all 1,231,665,708 database bytes identically to Source. One-thread
  scouts measured Source 96.39 s/1,252,152 KiB and native 90.96 s/1,403,520 KiB.
  These are single shared-host measurements, not a stable performance claim.
  Source logs prove two nonempty chunks; native phase events are not exposed,
  so its chunk count is not claimed as an independent runtime observation.
  Commands, limits, complete record/layout checks and crossed reads are under
  `.tmp/native-post-translation-20261107/scale-data/makedb-current/fixed-proof.json`.
- [x] Measure `dbinfo` speed and peak RSS. Essential output is identical; the
  startup-scale scout was 0.01 s/9,740 KiB C++ versus 0.00 s/2,560 KiB Rust.

## DAA viewing

Scouting evidence from 2026-10-04 and the producer follow-up from 2026-10-05
are summarized in `VALIDATION_RESULTS.md`. Native Rust now produces both
blastp and blastx DAA files. Compact compatibility-square checks render
byte-identical tabular alignments with either producer and either viewer.

Use realistic DAA input and compare native output with upstream for every
format:

- [x] DAA (`100`/`daa`) when viewing an upstream-produced archive.
- [x] Tabular (`6`/`tab`).
- [x] Flat JSON (`104`/`json-flat`), including byte-exact comparison of
  multi-query blastp archives from both producers.
- [x] PAF (`103`/`paf`).
- [x] SAM (`101`/`sam`).
- [x] Pairwise (`0`).
- [x] XML (`5`/`xml`).
- [x] Null output.
- [x] Edge output.
- [x] Scout speed and peak RSS for the practically important view formats.
- [x] Exercise unlimited `view -k 0` on a genuinely nonduplicated, large
  natural archive. The 509 MiB archive renders 1,678,045 byte-identical rows;
  this exposed and fixed Rust's former zero-results interpretation of `-k 0`.
  A later alternating null-output scout was stable: 1.77--1.80 s C++ versus
  6.85--6.86 s Rust, while Rust used only about 6.5 MiB RSS versus 821--863
  MiB C++. This is now a distinct speed/RSS tradeoff to optimize, not an
  unresolved timing artifact. `perf` attributes 36.5% of Rust cycles to
  `HspContext::parse`, then 10.8% to escaped title output and 5.8% to
  transcript-buffer delimiter scanning. Preserve streaming RSS while reducing
  repeated transcript parsing. Profile and timings are under
  `.tmp/cross-tool-rotation/`.
- [x] Restore upstream's ordered parallel DAA tabular-view execution boundary.
  Rust previously ignored `-p/--threads`, so the apparent default-view gap
  compared a serial Rust run with 40 upstream workers. Rust now uses a bounded,
  ordered producer/consumer pipeline capped at 32 resident queries and retains
  a direct one-thread streaming path.
  On the same natural archive, repeated `/dev/null` medians are 0.70/0.56 s
  C++/Rust at one thread, 0.44/0.27 s at four threads, and 1.11/0.81 s at the
  40-thread host default. Rust RSS ratios are 0.055x, 0.083x, and 0.091x.
  All modes retain the exact 458 MiB output hash. Artifacts are under
  `.tmp/view-ordered-workers/` and `.tmp/portfolio-view-*`.
  A follow-up tested upstream-sized queued windows at four threads. Windows of
  128 and 384 records improved contended runtime by about 24%, but raised RSS
  from roughly 33 MiB to 92--247 MiB, so both were rejected and fully reverted.
  The accepted 32-resident-record pipeline eliminates the batch barrier without
  retaining those wider windows.

Pinned upstream DIAMOND 2.1.24 has a viewer-specific bug: its `ViewWriter`
does not insert `OutputFormat::query_separator`, so multi-query `json-flat`
archives are malformed at query boundaries. The Rust viewer deliberately
matches those bytes rather than silently repairing the upstream format. Rust
also matches upstream's rejection of both `view --unal 0` and `view --unal 1`.

- [x] Native blastp DAA production readable by both C++ and Rust viewers.
- [x] Native blastx DAA production readable by both C++ and Rust viewers.
  On the natural 1 Mb translated-search workload, producer medians are
  4.62/3.62 s upstream/Rust and 268,764/188,564 KiB. The four-way producer/
  viewer square renders byte-identical TSV, and matched direct output confirms
  DAA serialization is not a performance deficit.
- [x] Repeat the complete producer/viewer and format matrix on the larger
  natural fixtures using the updated validation runner (56/56 exact or
  format-specific parity checks pass in `.tmp/ffi-final-view-byte-parity/`).

The viewer now has both synthetic Rust tests and realistic byte-for-byte
upstream comparisons. Multi-query flat JSON is intentionally tested as exact
upstream bytes even though those bytes are not valid JSON.

Remaining DAA/view work found by the 2026-10-12 breadth audit:

- [x] Apply upstream's ordered worker boundary to DAA, PAF, SAM, pairwise,
  XML, null, and edge output, not only tabular/flat JSON. The principal defect
  was also semantic at the internal boundary: Rust collected every decoded
  match before applying `-k 25`, while upstream streams and stops at the output
  range. The translation now uses 32-query tasks, at most `3 * threads`
  resident ordered tasks, and streaming match iterators. On the natural 509 MiB
  archive, three-run exact null medians are 0.54/0.26 seconds upstream/Rust at
  one thread and 0.46/0.21 at four threads; Rust RSS is 0.314x/0.207x. PAF is
  likewise exact and 1.68x/2.00x faster. Evidence is in
  `.tmp/view-null-workers-20261012/`.
- [x] Translate `--qnum-offset` and `--snum-offset` through native tabular and
  flat-JSON viewing. The parser mirrors upstream's `atoll` input behavior and
  its signed-`int64_t` to unsigned-`OId` conversion, including modulo-2^64
  output for negative offsets. Direct tabular and flat-JSON probes on the
  natural archive are byte exact; evidence is in
  `.tmp/view-offsets-20261012/`.
- [x] Give native `view` the same strict option/value validation as upstream:
  reject unknown switches and malformed numbers, enforce the `--top`/`-k`
  conflict and the 0--100 `--top` range, and validate `--header`. Direct probes
  now match upstream output or rejection for joined long options, C-style
  numeric prefixes, repeated options, filtering, headers, and formats.
- [x] Stream native `view` output to stdout when `-o` is absent. Tabular and
  redirected binary DAA stdout are byte-identical to upstream; the latter
  retains the format's required seek/finalization behavior.
- [x] Match upstream's `view` help and verbosity/quiet presentation. The
  native help text is byte-identical apart from upstream's expired-build ANSI
  reset, stdout output is quiet by default, file output uses normal verbosity,
  and `--verbose`, `--quiet`, and `--log` now select the upstream metadata,
  selection-limit, and debug diagnostics. Timer durations and the
  date-dependent expired-build warning are intentionally not reproduced.
- [x] Reject `--compress` and `--sam-query-len` for `view`, matching upstream's
  workflow parser (in addition to the existing `--unal` rejection), instead
  of accepting and silently ignoring them.

## Search sensitivity

### `blastp`

- [x] `--fast`.
- [x] `--mid-sensitive`.
- [x] `--sensitive`.
- [x] `--very-sensitive`.
- [x] `blastx --min-orf` on a steady-state natural workload. The 2026-10-12
  one-thread rotation used a contiguous 1 Mb *E. coli* query against 1,000
  natural AMR proteins. All six outputs were byte-identical; Rust was 4.34 s
  and 188,568 KiB at the median versus 5.05 s and 266,664 KiB upstream.

All listed natural-data presets match. The compact `--mid-sensitive` boundary
case is also byte-identical after restoring upstream's contiguous
`SequenceSet` fingerprint semantics.

The final matrix covers default and all seven named presets on compact and
natural inputs, with the representative multithreaded spot checks described in
`VALIDATION_RESULTS.md`.

A fresh post-optimization `--very-sensitive` rotation closes the last stale
timing deficit. On 300 disjoint natural queries by 1,000 targets, p1 medians are
3.99/3.96 seconds upstream/Rust with Rust at 0.730x RSS. A separate 1,200-query
p4 fixture is exact at 24,403 rows and 4.97/4.97 seconds, with Rust at 0.788x
RSS. Artifacts are under `.tmp/blastp-very-sensitive-refresh-20261012/`.

### Cluster option fidelity refresh (2026-10-12)

- [x] Match upstream's last-occurrence-wins parsing for repeated
  `--cluster-steps`; Rust previously appended repeated occurrences.
- [x] Reject native `--aln-out` without `--parallel-tmpdir` with the same
  error as upstream instead of emitting a non-upstream two-column file.
- [x] Wire the translated cluster controls that the manual CLI previously
  ignored: `--mutual-cover`, `--connected-component-depth`, `--no-reassign`,
  `--round-coverage`, `--round-approx-id`, `--kmer-ranking`, `--weighted-gvc`,
  `--strict-gvc`, `--graph-algo`, and `-M`/`--memory-limit`. Explicit and
  mutual cover now remain distinguishable so their upstream conflict is
  rejected identically.
- [x] Wire the already-translated `--oid-output` and simple `--header` cluster
  output controls instead of silently ignoring them. Unsupported verbose
  headers are rejected explicitly.
- [x] Exercise a genuinely nondefault three-round `--cluster-steps` pipeline
  on 4,000 distinct AMR proteins. Three alternating one-thread runs produced
  the same 4,000-line output byte-for-byte. Rust's median was 2.29 s /
  17,304 KiB; upstream's was 2.49 s / 23,160 KiB.
- [x] Exercise `--mutual-cover 80` on the same natural workload. This exposed
  and fixed two search-boundary errors: the aliased self-search block was
  length-sorted twice, and the derived minimum-length ratio was not forwarded
  to stage-1 dispatch. Every round's alignment and cluster counts now match;
  final output is byte-identical. Rust's median was 7.35 s / 18,740 KiB versus
  7.19 s / 23,156 KiB upstream.
- [x] Forward `--matrix`, `--gapopen`, and `--gapextend` through clustering
  instead of silently ignoring them. This was a demonstrated parity failure,
  not only an option audit: on the natural 500-protein fixture, default output
  was exact (`5f945c2b...`), but upstream `--matrix pam250` changed its output
  to `59d25e29...` while Rust emitted the unchanged default hash. The score
  state now passes through every native boundary; default, PAM250, and
  BLOSUM62 with accepted alternate gaps 9/2 are byte-identical to upstream.
  Upstream creates
  the global `ScoreMatrix` once in `basic/config.cpp` from the requested
  matrix/gaps before `Cascaded::run`; each round retains it and only updates
  its database-letter normalization. The minimal native path is therefore
  `parse_native_cluster_config` -> `ClusterConfig` -> `CascadedWrapperConfig`/
  `CascadedConfig` -> the `BlastpConfig` constructed by
  `NativeBlock::run_search`. Round handling continues to cap intermediate
  e-values at 0.001, disable CBS for Linclust step sets, and apply the existing
  sensitivity/extension mutations exactly as upstream does. Every config
  handoff is unit-tested and the former 500-protein PAM250 failure is retained
  as an end-to-end artifact. Profiling the initially exact but slow path found
  63.77% of sampled Rust cycles in TANTAN's `LambdaCalculator::inverse_at`:
  Rust reconstructed the immutable soft-masker once per cascaded round while
  upstream constructs its global `Masking` object once before
  `Cascaded::run`. The native backend now shares one command-lifetime masker
  across all round blocks. Three alternating exact-output pairs measure
  1.93 s/28,112 KiB upstream versus 1.86 s/7,680 KiB Rust (medians), closing
  the 6.89 s pre-fix runtime gap while retaining much lower Rust RSS.
- [x] Forward the remaining ordinary cluster search controls that upstream
  permits: `--evalue`, `--comp-based-stats`, `--masking`, `--soft-masking`,
  `--motif-masking`, and `--ext`. All six now pass through every cascaded
  round with upstream defaults and mutations intact. The 4,000-protein
  `--evalue 1e-20` acceptance is exact in all six alternating runs at median
  5.38/5.29 seconds and 25,992/17,816 KiB upstream/Rust. Independent natural
  controls for CBS 0, hard TANTAN, soft masking 0, motif masking 0, and full
  extension are also byte-identical.
- [ ] Verify the implemented external/multiprocessing cluster pipeline for
  `--parallel-tmpdir --aln-out`, including its real chunk/job/storage and richer
  alignment artifacts. Final-profile, multi-worker and scale parity remain
  open; this is no longer an absent native implementation.
- [ ] Verify split-input `linclust` within the documented Source UB policy.
  The following historical 2026-10-12 evidence must not be reclassified as an
  unimplemented adapter or a universally defined upstream result. The
  length-sorted streaming adapter forms the same natural
  4,000-protein blocks as upstream (block 1: 3,619 sequences and 1,302,291
  letters). The four second-block mismatches came from dropping upstream's
  `lin_stage1_target` state at the native stage-1 boundary: Rust evaluated the
  complete query/target seed-group cross product instead of selecting the
  first length-sorted target pivot with `stage1_target_lin`. The translated
  dispatcher now selects that exact kernel and also carries the wrapper's
  `freq_masking` state, which disables sensitivity-derived motif masking as in
  upstream. Three alternating `linclust -p1 -M 5M` pairs are byte-identical
  (SHA-256 `254ff1333d763d0a9c86ef94e4eae6550224e4fd66d3fe7d57b3fe7be8bfaa8f`):
  median upstream/Rust time is 1.59/1.28 s and RSS is 25,460/16,320 KiB
  (speed 1.242x, RSS ratio 0.641x). Evidence:
  `.tmp/cluster-lowmem-agent/README.md`. This closes the original 4/4,000
  target-selection gap, but a larger 10,000-protein, three-super-block `-M 5M`
  control remains non-identical (4,238 sorted lines differ). Its first-block
  cascaded round already reports 717 clusters in Rust versus 719 upstream,
  before any centroid-vs-later-block search, so it is a separate cascaded
  search fidelity gap and remains open. The earliest cause is now isolated:
  upstream reads `Search::Config::current_ref_block` before initializing it.
  On the 4,000-record/two-block process it length-sorts the first linear query
  block, while on the 10,000-record/three-block process it skips that sort.
  Rust deliberately keeps this state defined and therefore cannot reproduce
  both outputs by a verbatim source translation; an input-size heuristic was
  rejected. Forcing the observed 10,000-record order reduces the delta to
  three centroid assignments and makes every reported alignment/cluster count
  exact, but breaks 4,000-record parity, so that experiment was not retained.

### `blastx`

- [x] True default sensitivity with no sensitivity flag.
- [x] `--faster`.
- [x] `--fast`.
- [x] `--mid-sensitive`.
- [x] `--more-sensitive`.

The final matrix tests true-default `blastx` without a sensitivity flag and all
seven named presets. Historical harness aliases that added `--sensitive` are
not used as evidence for the true-default case.

## Scoring matrices

BLOSUM62 and BLOSUM45 have end-to-end comparison coverage. Add realistic
parity, speed, and RSS comparisons for:

- [x] BLOSUM50, including default composition statistics, statistics disabled,
  and an alternate 13/3 gap scheme.
- [x] BLOSUM80.
- [x] BLOSUM90.
- [x] PAM30.
- [x] PAM70.
- [x] PAM250.
- [x] Representative non-default gap-open and gap-extension combinations for
  each matrix were exercised.

## Search options and filters

- [ ] Reconstruct and rerun the historical 500-query mouse scoring fixture.
  The original checkpoint differed from its C++ snapshot in 137 scoring rows;
  the current binary passes later natural scoring matrices, but this exact
  corpus is unavailable and has not been revalidated.

- [x] `--max-target-seqs` with a non-default value (`-k 5`). The natural
  75-by-1,000 AMR scout is byte-identical (202 rows, SHA-256 `50471250...`).
- [x] Unlimited `--max-target-seqs 0` for native `blastp` and `view`.
  Native `blastp` emits the same 1,999 rows and SHA-256 `2934068a...` at one
  and four threads on the 75-by-1,000 AMR fixture; the larger viewer evidence
  is recorded above.
- [x] Strict non-default `--evalue` (`1e-20`): byte-identical on the same
  natural fixture (533 rows, SHA-256 `da63e217...`).
- [x] `--no-self-hits`. A 500-by-itself natural search first exposed that the
  translated final-filter path lacked target titles. The fixed path reproduces
  upstream's post-ranking placement and its exact 6,973 rows (SHA-256
  `3bba2f79...`), including the fact that a removed self hit is not backfilled.
- [x] Custom `--xdrop 20`: byte-identical (549 rows, SHA-256 `77046cda...`).
- [x] Non-default `--ext-chunk-size 16`: byte-identical (549 rows, SHA-256
  `77046cda...`).
- [x] `blastp --min-query-len`. Pinned upstream applies this globally while
  enumerating both query and reference seed arrays, despite the option name.
  Rust now preserves that observable behavior; the natural 300-by-7,000
  comparison is byte-identical at 5,941 rows (SHA-256 `1b7c7264...`).
- [x] Motif masking. Explicit on/off comparisons on 300 natural queries by
  7,000 targets are byte-identical to upstream. The switch changes six
  alignments (6,461 versus 6,467 rows), proving both paths were exercised;
  hashes and scout measurements are recorded in `VALIDATION_RESULTS.md`.
- [x] Combinations of filters rather than only isolated identity, coverage,
  and top-score thresholds. A natural 300-by-7,000 run combining identity,
  query/subject coverage, top-score, and strict e-value emits the same 39,516
  rows as upstream (SHA-256 `45aae572...`). Upstream leave-one-out controls
  prove identity and query coverage independently remove real alignments.

## Memory modes

- [x] Make a realistic `blastx` run cross a nonzero intermediate memory limit
  during stage 1, verifying migration of existing hits and continued
  disk-backed operation.
- [x] Compare parity, speed, and peak RSS for that adaptive `blastx` spill at
  one and multiple threads.

Default 16G and immediate `0G` spilling are measured for `blastx`. Adaptive
64M mid-run spilling is measured for `blastp`; the larger natural
SRR1703350/Swiss-Prot workload establishes delayed `blastx` spilling at
1150M. The current three-pair `blastp --sensitive -p4` refresh is exact and
confirms actual delayed migration/readback. After absolute-location filtering
and the shared worker histogram, upstream/Rust medians are 27.59/27.41 seconds
and 92,716/71,148 KiB: Rust is at runtime parity while using 23.3% less RSS.
Artifacts are under `.tmp/delayed-spill-hist-refresh/`; the immediately prior
profile under `/tmp/diamond-spill-refresh/` confirms spill compression and file
I/O remain below 0.01% of CPU.

## Output fields

Add byte-parity coverage for native tabular fields not included in the
current extended-output comparison:

- [x] `positive`.
- [x] `ppos`.
- [x] `qcovhsp`.
- [x] `scovhsp`.
- [x] Query accession variants (`qacc`, `qaccver`). Pinned upstream DIAMOND
  rejects these BLAST aliases as invalid fields; Rust now rejects them too.
- [x] Subject accession variants (`sacc`, `saccver`). Pinned upstream likewise
  does not register these names, and Rust no longer accepts them as aliases.
- [x] Test useful combinations of these fields rather than only one field at
  a time. A combined natural 300-by-7,000 search emits 6,461 exact rows with
  SHA-256 `f49a90dd73d42a038a539e3a24270f09ac1d91fc14da2d3fc8831db85603d9fe`.

## Input handling

- [x] Multiple query files in one invocation. Native Rust now matches pinned
  upstream's synchronized paired-end semantics: R1/R2 records are interleaved,
  unequal counts and more than two files are rejected, and `full_qseq_mate`
  prints the opposite ordinal mate without duplicating sequence storage. On
  5,000 genuine SRR1703350 pairs against 7,000 AMR proteins, all one- and
  four-thread outputs contain the same 3,081 rows and SHA-256 `23c16303...`.
  The translated query path now follows upstream's ownership model: one
  authoritative `Block` owns the source DNA and all six translated contexts,
  masking/restoration occurs in place, and extension borrows slices from that
  storage. The restoration path also distinguishes `SEED_MASK` from
  `LETTER_MASK`, preserving propagated seed masks exactly. Focused
  motif-disabled multi-shape and target-indexed checks remain byte-identical
  at 25 rows, with SHA-256 `c4b1cad2e798...` and `7890ad7a0ab1...`
  respectively.
  Three balanced one-thread runs measure 6.38 s / 47,784 KiB upstream,
  6.07 s / 36,656 KiB for Rust at the default 16G limit (speed 1.051x, RSS
  0.767x), and 6.49 s / 37,408 KiB for Rust at forced-disk 0G (speed 0.983x,
  RSS 0.783x). Thus the former dual-ownership RSS excess is closed in both
  memory modes while exact 3,081-row output is retained.
  After the retained absolute-location and faithful all-shape worker-histogram
  translations, fresh four-thread medians are 2.64 s/62,824 KiB upstream
  versus 2.30 s/53,580 KiB Rust (speed 1.148x, RSS ratio 0.853x). An earlier
  worker-histogram attempt changed the output hash and regressed to 5.54 s; it
  was reverted rather than reused. The accepted implementation preserves
  upstream worker sequence ranges and exact per-partition seed order.
- [x] Multiple reference arguments for `makedb`. Pinned upstream rejects two
  values in one `--in` occurrence and treats repeated `--in` occurrences as
  assignment (the last wins), rather than concatenation. Rust now preserves
  both behaviors; its last-input database is byte-identical to upstream.
- [x] Compressed `blastx` query input. The existing contiguous natural 1 Mbp
  *E. coli* query was gzip-compressed without changing its records and searched
  at `--mid-sensitive`; upstream and Rust emit the same 25 rows and SHA-256
  `c4b1cad2...`, also matching the uncompressed-input result.
- [x] Compressed multi-file input combinations. On 5,000 natural SRR1703350
  pairs, gzip/plain, plain/gzip, and gzip/gzip `blastx` inputs all match the
  all-plain baseline and upstream exactly: 3,081 rows, 1,111,055 bytes, SHA-256
  `23c16303...`. Repeated mixed gzip/plain `makedb --in` occurrences preserve
  upstream's last-occurrence-wins behavior in both orders and produce
  byte-identical databases. Supplying two mixed values in one occurrence is
  rejected by both programs with exit 1 and the same core error. Artifacts are
  under `.tmp/compressed-multifile-parity/`.

## Completion criteria

For each performance-relevant item:

- [ ] Use natural, nonduplicated data.
- [ ] Run both one-thread and multithreaded comparisons where applicable.
- [ ] Alternate C++/Rust execution order and repeat runs on the shared host.
- [ ] Record median wall time and peak RSS, plus C++/Rust ratios.
- [ ] Check byte-identical output for every paired run, or document the exact
  semantic reason byte identity is inapplicable.
- [ ] Test both the default 16G Rust memory limit and `0G` forced spilling for
  search modes where retained hits are involved.
- [ ] Add stable results to the README and remove the corresponding checklist
  item.

## Current performance portfolio (2026-10-12)

- [x] Recheck a startup-independent one-thread translated search on a natural
  contiguous genomic input. After the full-key hash-join correction, three
  exact 1 Mb `blastx --mid-sensitive` pairs give 4.49 s C++ versus 3.64 s
  Rust (1.234x), with Rust RSS at 0.703x.
- [x] Recheck all three clustering commands after the shared stage-1 changes.
  The compact/natural matrix passes 18/18; clustering is not the current
  worst measured runtime offender.
- [x] Exercise DAA production plus every supported view format after the
  search changes. The format matrix remains exact, and large pairwise view is
  at runtime parity with much lower Rust RSS.
- [x] Close the four-thread `blastp --sensitive --memory-limit 0G` stage-1
  gap. Matching upstream's whole-worker AVX2 dispatch, rather than mixing in
  an independently selected AVX-512 fingerprint kernel, produced a previous
  near-parity series. Matched profiles then exposed a translation mismatch:
  Rust materialized a transient 32-byte `SeedMatch` for every survivor, while
  upstream passes compact target indices from `HitField` directly into
  `search_query_offset`. Restoring a compact per-query boundary improves the
  current exact medians from 28.75/34.02 s to 27.44/31.32 s C++/Rust (0.876x).
  Rust RSS remains lower at 70,832 versus 89,452 KiB. A final profile found
  that Rust scalarized every partial left-most-filter tail while upstream
  performs a complete 16-byte SIMD load from delimiter-padded sequence
  storage and masks the unused lanes. Restoring that exact access pattern and
  fusing the whole-window SSSE3 reduction closes the gap: current exact
  medians are 27.36 s C++ versus 27.85 s Rust (0.982x), with 94,120 versus
  69,588 KiB RSS (0.739x).
  A faithful single-query stage-2 boundary was also restored, but six combined
  alternating pairs show it is performance-neutral within host variance
  (approximately 28.43 s C++ versus 29.53 s Rust). A compact-target scratch
  rewrite regressed to 31.25 s and was removed.
- [x] Repeat the natural one-thread `blastp --ultra-sensitive` row after the
  shared hash-join fix. Three exact pairs give 4.62 s C++ versus 2.94 s Rust
  (1.571x), with Rust RSS at 0.547x. The later heavier 300-by-1,000 acceptance
  row is the current comparison: 16.07/15.45 seconds and 29,872/22,208 KiB.
- [x] Reduce `makedb` peak RSS on the 574,627-record natural Swiss-Prot
  fixture. Rust had retained the original records plus a second complete
  masked copy; upstream masks its owned block in place. Restoring that lifetime
  reduces the Rust median from 562,768 to 350,400 KiB. Fresh exact C++/Rust
  medians are 16.34/15.83 s and 402,172/350,720 KiB (speed 1.032x, RSS
  0.872x).
- [x] Keep optimization passes portfolio-based: after any shared search-core
  rewrite, rerun at least one clustering case and one DAA/view case before
  further search-only tuning. After the padded left-most-filter rewrite, the
  4,000-record `cluster` control and 24 MiB DAA tabular-view control retain
  their established exact SHA-256 values. The 509 MiB viewer was then
  profiled, a non-beneficial pointer-iterator experiment was removed, and work
  rotated through 10,000-record `linclust` plus natural BLOSUM50/PAM250
  scoring controls. Both scoring outputs are byte-identical; details are
  recorded in `VALIDATION_RESULTS.md`.
- [x] Separate persistent-index RSS caused by Rust's 16G in-memory default from
  the translated indexed path itself. On 300 natural queries by 7,000 targets,
  the pre-allocation-fix three-pair checkpoint gave 4.03 s/46,760 KiB
  upstream versus 4.64 s/82,908 KiB Rust. `--memory-limit 0G` lowered Rust to
  68,360 KiB but left a measured 4.69 s versus 3.95 s gap. That rotation was
  bounded while clustering and DAA/view received their turns; artifacts are
  under `.tmp/rotation-target-indexed-compact/` and
  `.tmp/rotation-target-indexed-0g/`.
  A later faithful allocation fix removes the complete duplicate reference
  seed array: forced-disk Rust now uses 46,404 KiB versus 46,260 KiB upstream
  (1.003x), while default Rust uses 64,316 KiB. Exact parity is unchanged;
  forced-disk runtime remains slower at 4.88 versus 3.95 seconds. Updated
  artifacts are under `.tmp/indexed-rss-agent/`.
- [x] Profile and optimize native tabular `view` on a larger archive without
  regressing the already-fast pairwise formatter. The merged 30 MiB archive
  decodes to an exact 133,214-line table. Streaming matches like upstream and
  replacing per-row generic title search reduced three exact Rust runs to a
  0.67 s median from the initial 0.92 s scout; upstream's median was 0.16 s
  (speed 0.239x). Rust RSS remains much lower at 5,120 versus 114,844 KiB.
  `perf` now leaves traceback parsing as the largest cost. Repeat on a 3--5
  second workload before accepting a deeper context/traceback rewrite. The
  subsequent 509 MiB natural archive provides that workload: C++ is
  1.77--1.80 s at 841--883 MiB RSS. Direct borrowed-query HSP parsing improves
  Rust from 6.85--6.86 s to 6.44--6.60 s at about 6 MiB RSS.
  `HspContext::parse` accounts for 36.5% of sampled Rust cycles; preserve
  streaming instead of buying speed by materializing the archive. Direct HSP
  parsing and a reusable per-query output context are now implemented. The
  latter improves a controlled three-pair Rust A/B median from 7.59 s to 6.99
  s and RSS from 6,508 to 6,252 KiB while retaining the exact 458 MiB output.
  The later ordered-worker translation supersedes the serial/default
  comparison: repeated compute-focused medians are documented above, while
  the direct one-thread path preserves the low streaming RSS.
- [x] Check the apparent PAM250 scoring gap on a startup-independent workload.
  It disappears as useful work grows: at 7,000 references and 3,000 disjoint
  natural queries, three exact pairs give 2.30 s C++ versus 2.17 s Rust
  (1.060x), with Rust RSS at 0.493x. The small 1,000-by-300 row was dominated
  by fixed/runtime overhead and is not an optimization target.
- [x] Repeat `deepclust --approx-id 50 --member-cover 60`. Five exact natural
  pairs give median times of 1.60 s C++ and 1.59 s Rust, with Rust RSS at
  6,708 versus 27,036 KiB. The original single-run gap was shared-host noise.
- [x] Scale the clustering rotation beyond the 500-protein scout. On 4,000
  distinct AMRFinderPlus proteins, profiling found that Rust rebuilt the
  immutable query-length `KmerRanking` for every joined seed bucket instead of
  once per shape as upstream does. Restoring that lifetime reduces one-thread
  `cluster` from about 9.0 s to a 6.00 s median versus 5.56 s upstream. The
  four-thread medians are 2.78 s Rust versus 2.95 s upstream. All outputs are
  byte-identical, and Rust RSS ratios are 0.691x and 0.276x respectively.
  Natural 4,000-record `linclust` and `deepclust` controls also remain exact.
  A later bounded seed-extraction pass translated upstream's prevalidated raw
  iterator for the production weight-12/base-10 shape. Five alternating
  old/new Rust pairs retained the exact output and reduced the median from
  5.57 to 5.34 s (4.1%), with no added allocation. This experiment is closed;
  fresh complete C++/Rust medians are 5.48/5.36 s at one thread and
  3.03/2.63 s at four threads. Further work rotates to a different workflow.
  A separate three-pair default one-thread `deepclust` refresh is also exact:
  medians are 9.42 s/28,152 KiB upstream and 9.12 s/18,688 KiB Rust. The Rust
  timing spread was 7.67--9.56 s under shared-host contention, so only parity
  and the lower 0.664x RSS ratio are considered firm conclusions.
- [x] Calibrate and optimize `linclust` on larger natural, nonduplicated input.
  The translated libstdc++ introsort processed every sketch candidate even
  though faster mode retains only the bottom 21. An exact cutoff selection now
  sorts the small retained set and falls back to the original exact introsort
  when equal hashes straddle the boundary. On 7,000 records Rust improves from
  2.88 s to a 2.14 s median. On the 10,000-record acceptance input, the exact
  pre-i16-profile one-thread medians were 2.80 s C++ versus 3.14 s Rust, and
  four-thread medians were 1.42 s C++ versus 1.19 s Rust. The output SHA is
  identical in every run.
  A faithful follow-up translated upstream's dynamic refill of completed AVX2
  endpoint lanes, but the 10,000-record median changed only from 3.14 to
  3.13 s; the experiment was fully reverted. The remaining one-thread gap is
  inside the i16 endpoint kernel/profile construction, not lane scheduling.
  That profile construction has now been translated: AVX2 byte-matrix shuffles
  replace Rust's scalar indexed loads. Exact old/new Rust medians improve from
  3.48 to 3.10 s, and fresh upstream/Rust medians are 2.97/3.11 s at one
  thread and 1.37/1.02 s at four threads, with essentially unchanged Rust RSS.
  Finally, the in-memory clustering search now mirrors upstream's aliased
  query/target block instead of cloning all records for self-search. Five
  exact old/new runs improve from 3.13 s/38,756 KiB to 2.61 s/33,980 KiB.
  Fresh upstream/Rust medians are 2.93/2.79 s and 35,932/34,104 KiB, so the
  one-thread speed and RSS deficits are closed. Further optimization rotates
  to other commands rather than continuing to tune `linclust`.
  A separate balanced `deepclust` refresh through the same backend is exact
  and measures 9.21/8.90 s upstream/Rust with RSS 27,912/17,152 KiB; its
  whole-command lane is closed as well.
- [x] Make native `makedb --in /dev/stdin` accept non-seekable input like
  upstream. Compression detection now peeks through `BufRead::fill_buf`
  without consuming or seeking. Plain and gzip pipes both pass; the streamed
  10,000-record DMND is byte-identical to upstream, and its `linclust` output
  has the same SHA-256 as the timed acceptance runs.

The post-rebuild cross-command checkpoint is retained under
`.tmp/validation-cross-tool-natural-checkpoint/`: 164/164 rows pass exact
parity across search sensitivities, clustering, scoring, DAA production, and
view.  This is the required rotation checkpoint before another focused
optimization pass.

- [x] Remove the standard-search `Stage1Candidate` survivor materialization.
  Rust now passes ordered `u32` target indices and borrows the original target
  arrays through stage 2, matching upstream. Three exact A/B pairs improve
  16.53 to 15.70 s and reduce RSS 28,924 to 28,056 KiB. Fresh exact
  upstream/Rust medians are 16.07/15.45 s and 29,872/22,208 KiB, closing the
  larger natural 300-by-1,000 one-thread `blastp --ultra-sensitive` lane.
- [x] Translate upstream's in-memory hit-bin layout for the Rust 16G fast path.
  On target-indexed 300-by-7,000 p1, Rust 0G matches ordinary upstream RSS
  (46,404 versus 46,732 KiB), while Rust 16G uses 63,912 KiB. A like-for-like
  upstream `--hit-membuf` control uses 58,848 KiB. Rust currently retains
  895,237 hits in 300 geometrically growing aligned vectors; upstream uses 16
  coarse pre-reserved packed bins. Rust now mirrors that 15-byte, 16-bin
  ownership layout and consumes bins in stable query order. Exact old/new
  medians are 4.17/4.18 s while RSS falls 63,520 to 56,776 KiB (10.6%); this
  is below the like-for-like upstream `--hit-membuf` median of 58,848 KiB.
  Forced 0G spilling remains exact at 46,404 KiB.
- [x] Test hoisting the left-most fingerprint ISA choice into the existing
  stage-1 dispatch. The refreshed 64M delayed-spill workload is exact at 28.58/30.72 s
  upstream/Rust and 92,944/77,012 KiB; spill writing is only 0.10% of Rust CPU.
  Rust currently rechecks AVX2 support and crosses an out-of-line wrapper for
  every alternative seed even though its caller is already the dispatched
  AVX2 filter. A faithful const-specialized scout was exact but slowed the
  median from 5.38 to 5.57 s (3.5%), so it was reverted. Do not retry dispatch
  hoisting alone; the residual stage-1 gap needs a different translation
  hypothesis.

### Cross-command stage-1 rotation — 2026-10-12

- [x] Keep upstream-style absolute `SequenceSet` positions through fingerprint
  rejection and recover record-local coordinates only for surviving hits.
  The natural 300-by-4,000 sensitive p4/64M scout is exact in all six runs and
  improves the Rust median from 6.03 to 5.73 seconds (5.0%) while lowering RSS
  from 38,320 to 37,780 KiB. Exact acceptance controls improve paired p4
  `blastx` from 3.16 to 3.10 seconds, p1 ultra-sensitive `blastp` from 16.02
  to 15.72 seconds, and the sampled p4 `linclust` run from 1.15 to 1.09
  seconds. A separate contiguous p4 `blastx` control changed from 6.47 to 6.56
  seconds (1.4% slower) with flat RSS, so the portfolio favors retention.
  Artifacts are under `.tmp/absolute-loc-scout/`.
- [x] Reject the fused motif-mask lifecycle experiment for `makeidx`. It
  retained the exact 134,217,792-byte index and SHA-256 `dbe81bcd...`, but
  regressed from 8.05 to 9.64 seconds with flat RSS. The candidate was fully
  reverted; artifacts remain under `.tmp/makeidx-fused-motif/`.
- [x] Rotate through unrelated commands before further search tuning. Fresh
  exact p1 controls find no material deficit: `cluster` is 5.46/5.34 seconds
  upstream/Rust with Rust at 0.646x RSS; `linclust` is 2.81/2.93 seconds with
  Rust at 0.919x RSS; unlimited tabular `view` is 10.23/7.57 seconds with Rust
  at 0.036x RSS; and pairwise `view` is 1.10/0.99 seconds with Rust at 0.059x
  RSS. Clustering and DAA/view are closed for this pass.
- [x] Translate upstream's reusable multithread seed histogram rather than
  recounting every query and reference seed inside each of four index chunks.
  The histogram retains worker-by-shape-by-partition counts, and its prefix
  ranges let parallel fills preserve exact record and seed order. On paired
  natural p4 `blastx`, all six Rust A/B outputs retain SHA-256 `23c16303...`;
  median runtime improves 3.04 to 2.30 seconds (24.3%) while RSS changes only
  52,996 to 53,580 KiB. The independent sensitive p4 control improves 5.51 to
  5.06 seconds and remains exact. Artifacts are under
  `.tmp/blastx-histogram-reuse-20261012/`.
- [x] Specialize the fixed production `makeidx` hashed enumerator. The generic
  skip/filter path remains available, while the no-skip/no-filter path hoists
  invariant decisions outside the seed loop. All six 129 MiB indexes retain
  SHA-256 `dbe81bcd...`; Rust improves 7.79 to 7.45 seconds with flat RSS.
  At that checkpoint upstream was 7.08 seconds, leaving a 5.2% gap. The later
  const-modulus masking translation documented above closes the p1 gap and
  leaves only 2.1% in the default-thread row. Artifacts are under
  `.tmp/makeidx-specialization-accept/` and
  `.tmp/makeidx-const-modulus-20261012/`.
