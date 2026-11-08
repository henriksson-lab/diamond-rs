# Translation progress

Historical checkpoint: 2026-09-26

Current inventory correction (2026-11-07): the five removed test modules have
been restored with self-contained corpus literals and focused tests. Their
file-map rows describe the current implementation and its limits. A concrete
native regression backend and `test` CLI now pass all twenty source hashes,
including actual blocking and query-indexed search. Explicit protein and DNA
file blocking uses native intermediate files and the block joiner. Source
default budgets, supplied translated blocks and blocked taxonomy are connected;
BLAST/clustering adapters and scale checks remain unfinished. The internal debug
logging path is explicit and source standard-test CLI rejects `--log`. The historical
198-file count below also omitted `lib/tantan/LambdaCalculator.cc`, so it does
not establish full inventory or end-to-end completion. Use
`../TRANSLATION_COMPLETION.md` for active integration and parity work.

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

At the historical checkpoint: none. The current inventory is missing the Rust
counterparts for upstream `test/queue.cpp`, `test/test_cases.cpp`,
`test/hit_buffer_stress.cpp`, `test/data.cpp`, and `test/test.cpp`.
End-to-end integration and parity work also remains; source-file audit
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

## 2026-11-07 native output-context checkpoint

The concrete tabular/flat JSON backend is connected to protein and translated
search, including blocked merging. Frozen source fixtures cover all 52 standard
static fields and 48 dynamic taxonomy rank fields, original sequence/quality
payloads, phase-specific self-alignment scores, masking/indexing differences,
unaligned callbacks and the source's invalid JSON separators/quoting. The focused
five-test suite passes 216 byte comparisons plus ten protocol/error checks.
See `tests/data/tabular_context/README.md` and `TRANSLATION_COMPLETION.md` for
bounded coverage and compatibility decisions. Earlier caveats above are
historical; use the live dispatcher for current fallback decisions.

The inventory audit currently finds 199 original implementation files and 198
mapped files, with `diamond/src/lib/tantan/LambdaCalculator.cc` still unmapped.
An existing Rust module is not sufficient evidence to close that audit. Overall
CLI integration/parity and removal of production FFI remain active work.

## 2026-11-07 native SAM search checkpoint

`SamFormat` now serves native `blastp`/`blastx` output `101`/`sam`, with actual
HSP transcripts, invocation headers, unaligned callbacks, `--sam-query-len`
and SAM `--salltitles`. Both direct and intermediate output paths are connected;
source dictionary title truncation is preserved. Debug evidence covers 168
frozen source CLI comparisons, three in-process dispatcher calls, four flag
errors, and two supplied translated-query block comparisons. See
`tests/data/sam_search/README.md` for exact scope. Other production fallback
routes and the unmapped `LambdaCalculator.cc` audit remain open.

## 2026-11-07 native edge/null checkpoint

The concrete binary Edge and NullFormat callbacks are connected to native
protein/translated search with their source coordinate/transcript DP requirements.
Original OIDs survive merging and coverage-based length sorting. The fixture
matrix covers 196 exact output comparisons, three in-process dispatcher repeats
and two malformed-query diagnostics. It exposed corrected-score FMA rounding
and translated reverse-pass coverage bugs; both have focused passing regressions.
See `tests/data/edge_null/README.md` for scope and compatibility decisions.
The taxonomy-format/advanced-mode fallback routes and inventory audit remain open.

## 2026-11-07 native taxonomy classification checkpoint

Format 102 and --include-lineage now use the concrete native per-query
classifier with actual metadata, conditional name loading, NONE DP requirements
and source top-score defaults. Explicit top/score overrides work through CLI
and supplied-query APIs; unlimited-target chunk rounding now uses the translated
overflow-safe helper. Focused validation covers 284 direct source-output
comparisons, three dispatcher repeats and ten metadata/option errors. Eight
additional checks validate a defined native result against safe source variants
for a documented title-sorting source-UB boundary; they are counted separately
from direct parity. See tests/data/taxon_search/README.md for evidence and scope.
All standard baseline search formats are connected. Advanced fallback routes,
optional modes and the LambdaCalculator.cc inventory audit remain open.


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

## Controls and EXTRA code-first batch (2026-11-07)

Shared native search tuning/scoring, database filtering, range/species culling,
aligned/unaligned query output, recursive clustering configuration and Source
I/O controls are now connected. Default-off `extra` connects six more tools,
including the Source smith-waterman unknown-DNA-matrix error. CLI repeated
scalar/vector options now use Source last-assignment semantics. External
linclust and multiprocess persistence/workflow adapters are in progress.
Broad runtime/parity testing remains deferred per user instruction; see
TRANSLATION_COMPLETION.md for exact behavior and Source-UB compatibility scope.

Combined `cargo check --offline --locked --all-targets --features extra` passes
for this batch (16.40 seconds); compile log and source hashes are retained in
`.tmp/native-options-20261107/`. The second EXTRA batch also connects split,
hashseqs, listseeds, wc, cut and the Source deprecated roc-id error. This check
compiles existing tests without executing them; it is not conformance evidence.

## External/distributed and concrete startup batch (2026-11-07)

External linclust and multiprocess DMND search now connect their real native
file/queue/dictionary/search/join pipelines. Optional commands use concrete
translated Config dependencies; find-shapes and Source main-switch fallthroughs
are connected. Shared all-target check with `extra` passes (16.39 seconds),
retained with source hashes in `.tmp/native-options-20261107/distributed-*`.
This establishes compilation only. Runtime/parity/recovery and platform/package
checks remain deferred; the full translation goal is still active.

## Native I/O, logging and remaining control connections (2026-11-07)

Ordinary clustering and recursive adapters forward configured buffers,
temporary directories, no-unlink, warnings and full search configuration.
Actual DP limits/backtrace settings and matrix-based TANTAN thresholds reach
ordinary/global/external/realign paths. Long-read defaults, search logging and
linear target search are connected. Shared all-target check with `extra`
passes in 15.37 seconds (`io-logging-compile.log`). Later output-offset/XML,
Source-ignored option and clustering ranking-budget edits await the next shared
check. These are code connections; runtime parity is still deliberately deferred.

The optional DNA dependency audit found missing KSW/WFA source dependencies in
the adjacent checkout. Its existing Rust wrappers also use surrogate algorithms
that omit z-drop/band semantics. Actual KSW recurrence translation is in progress
against a separately pinned canonical primary source; no DNA completion or
pinned-build parity claim is made for those absent dependencies.

Linear query/combo stage dispatch, real extension selection, output ID controls,
Source-ignored option parsing, cluster ranking/partition settings and actual
scalar KSW dependency wrappers now pass the combined all-target `extra` check
(12.43 seconds; `.tmp/native-options-20261107/linear-output-*`). DNA SSE/WFA/CLI,
parallel decode, clustering-specific block orchestration and final conformance
remain active translation work. Broad tests remain deferred per instruction.
The same code also passes the default-feature all-target check (16.65 seconds;
`default-compile.log`). Neither check executes tests or requires C++ FFI.

## Actual loader/clustering and database controls batch (2026-11-07)

Actual cluster query/reference blocks and disk joins, ordered parallel BLAST
raw/decode and hit-packet pipelines, frameshift ranking ratio, pre-join target
tracking/accession output, makedb runtime/stdin/warning controls and GNU KSW
state-machine port now pass the combined `extra` all-target check (14.97s).
Logs/source snapshot are in `.tmp/native-options-20261107/loader-cluster-*`.
Source-defined EDGE tag and native invalid-value errors replace two remaining
misrepresentations. Broad tests remain deferred. Heartbeat live counters, WFA,
DNA CLI/reporting and final parity/platform/package gates remain open.

Header/writer/common parsing and native-only production routing pass the
`extra` all-target check (10.51s). The actual DNA command behind default-off
`dna` passes its all-target check (20.42s). Source ABI, default scoring and
DNA quirks/undefined boundaries are detailed in TRANSLATION_COMPLETION.md;
no runtime parity or dependency-fork identity is claimed. Stdin integration is
prepared; AUTO, DNA MPP/timers, optional MCL and final conformance remain active.

## AUTO, DNA multiprocessing and concrete MCL batch (2026-11-07)

The combined all-target `dna,mcl` compile passes (20.11s; log
`.tmp/native-options-20261107/mcl-auto-getseq-compile-fixed.log`). This connects
command-owned AUTO/cache selection across query rounds, empty-reference and
MPP no-work paths, and fixes the direct single-reference prepared-state
transition. DNA now uses the shared real MPP queues/dictionaries/join lifecycle
and actual timing counters. Optional default-off `mcl` connects CLI parsing,
native streamed search, bounded graph checkpoints, component processing and
expression thresholds before HSP culling, retaining pre-lazy target self scores.
Source's MCL MPP chunk-output/empty-graph behavior is preserved explicitly.

EXTRA getseq now implements title/OID selection, reverse/hard-masked output,
letter budgets and actual file buffering; its FASTA source is a bounded pull
reader. Source-specific title/default-output and unterminated-list quirks are
retained. Shared startup now maps getseq and resolves explicit zero threads
to hardware concurrency. The default-feature all-target compile also passes
(17.08s; `mcl-auto-getseq-default-compile.log`). No broad tests were run. Callgraph translation,
DNA dependency/setup provenance and final parity/platform/package gates remain
open; this compile is integration evidence, not a parity claim.

## Concrete lifecycle and workflow batch (2026-11-07)

The next combined all-target `dna,mcl` check passes in 16.15s
(`.tmp/native-options-20261107/lifecycle-workflow-compile-fixed.log`). It covers
Source DNA sorting/OID permutations, first-shape seed configuration, asymmetric
reference reduction, full-database HSP statistics and scaled WFA metadata;
MCL's actual worker/cluster-ID scheduling, component spill deduplication and
sparse numerical paths; per-invocation search pools and actual stage2 windows;
shared reassign/recluster tuning and realign scoring/DP/buffering; and configured
optional-tool IO. MCL's zero-worker/progress Source boundaries use explicitly
defined native behavior; its IDs retain Source scheduling/container sensitivity.

Makedb now prepares compression/format readers before destination creation,
forces Source's 4 MiB capacity and carries its actual temporary directory into
all taxonomy sort relations. Getseq uses the same prepared input boundary.
The concrete regression backend reuses production configuration/logging and
EXTRA stdout mode. No runtime suite was executed. Subsequent generic native
dispatcher integration, bounded TSV tool input, DAA IO controls and Source ABI
spill accounting are ongoing and are not covered by this completed check.

## Concrete native entry point and bounded tools (2026-11-07)

The actual CLI/library path now prepares typed commands, enters the Rust FFI
facade and translated generic main dispatcher, and executes concrete native
backend methods. Invocation-owned logging and command locks are released on
both success and error, allowing backend reuse. EXTRA preparation is separate
from actual tool execution. This is no longer an unused trait or a callback
into the old dispatcher. The all-target `dna,mcl` check passes in 12.24s and
default-feature check in 16.37s (`native-entry-tsv-daa-*-compile*.log`).

The same batch connects actual buffered DAA input/output/finalization, bounded
parallel TSV parsing/projection, Source input-open logging order and zero-buffer
diagnostics. ALP matrix/vector helpers now use six-significant-digit primitive
float output and formatted numeric-prefix extraction, retaining Source's
different machine/symmetric input failure checks. Taxonomy spill thresholds
account for the C++ string/pair ABI instead of Rust object layouts. No runtime
suite was run. Subsequent clustered global/iteration orchestration, full
workflow configuration, mapping IO and Source FASTA codec work remain outside
these compile results until the next aggregate check.

## Clustered rounds, full workflow controls and Source codecs (2026-11-07)

The next all-target `dna,mcl` check passes in 15.02s
(`cluster-rounds-codecs-compile-fixed.log`). Reassign/recluster now use the real
common search configuration, reapply Source's workflow-specific overrides and
original-database scoring size, and submit all queries to one edge search.
That concrete driver owns sensitivity rounds, original OIDs, pre-join aligned
query state, ordinary intermediate merging and actual global-ranking tables/DP.
Global extension borrows the existing database and copies only the ranked OID
union; direct typed edge capture occurs at the real HSP boundary in query order.
Live ranking/load-balancing controls now reach their actual native consumers.

Mapping input now streams configured/compressed lines with Source's count and
tokenizer behavior; mapping output stages one textual row. FASTA prepared and
streamed inputs reuse the translated codec stack, including raw zlib and mixed
members. Native Source paths preserve regular-file-only compression detection;
typed generic-reader convenience retains explicit autodetection. Source FASTA
paths enforce the pinned zstd-off profile; other command IO profile boundaries
remain a subsequent integration task. Low-level Rust codecs stay available.

No runtime suite was executed. The broader shared BasicConfig backend,
add-self-alignment context, remaining command-profile IO, final feature parity
and platform/package/performance gates remain open.

## Shared Source runtime, database owners and command codec profile (2026-11-07)

The all-target `dna,mcl,zstd` check passes in 20.65s and the default-feature
all-target check in 20.43s (`source-runtime-owner-compile-fixed.log` and
`source-runtime-owner-default-compile.log` under `.tmp/native-options-20261107`).
This batch connects the real shared BasicConfig startup to initialized score
matrix, masking and translation resources, preserves Source logging transitions
across regression constructors, and routes benchmark execution through those
resources. Command IO now applies the explicit default-off Source zstd profile;
low-level Rust codecs remain available independently. Makdb opens and guesses
its input before truncating its destination and forwards temporary directories
to the actual taxonomy sorters.

Clustering now retains the actual original or temporary sequence-file owner,
original OIDs, active filters, real MPP partition/chunk execution and complete
common search configuration. Self-alignment contexts reach ordinary, global,
iterated and MPP consumers. Source-invalid owner metadata access is reported at
the actual access boundary rather than silently loading a different resource.
The original C++ tree, objects and executable were not rebuilt or modified.

These are compile results, not a new runtime parity claim. Subsequent CBS,
alignment scheduling, diagnostics, chaining, DP reservation and radix/hash join
translations remain outside this evidence until the next shared check.

## Actual join recursion and DP/search control integration (2026-11-07)

The combined all-target `dna,mcl,zstd` check passes in 17.00s
(`source-join-dp-controls-compile-fixed.log`). The first attempt identified only
a configuration equality derive involving the new shared DP failure channel;
its wrapper now compares channel identity and preserves cloned error propagation.
The successful batch includes actual ascending radix join recursion, reusable
scratch arenas, stable linear group scatter, packed value widths, hash-table
failure propagation and all six join controls reaching real stage0 partitions.
Nondefault controls select that path; default optimized seed-join equivalence
remains a general parity gate. Source's buffered-radix implementation drops
partial eight-record buffers and may read uninitialized gaps; native bounded
tail scatter is an explicit compatibility decision, not literal parity there.

The same batch connects whole-query trace-point task scheduling, `palign`, actual
DP PARALLEL batch workers and fallible thread-pool failures before output;
Source query/extension diagnostics and target-parallel timers; all live chaining
controls including backtrace ratios; Source band bins and stable target order;
minimum/narrow band controls; CBS solver limits; transcript reservations and
atomic frameshift chunk scheduling with serial overflow retry after worker joins.
Hardware thread defaults now follow the host C++ count rather than Rust quota
availability. Forty-six fields with no Source consumer (including no constructor
consumer) retain their accepted no-effect behavior after canonical registration
checks. Kmer ranking reaches the generic path, including its real PackedLoc
unsupported boundary. No broad runtime suite was executed.

Code translation is still incomplete: the audit identifies live CBS window,
automatic motif cap, explicit length ratio, ordinary diagonal-filter controls
and complete self-query database ownership as the next concrete integration
batch. Final parity, feature-profile, Apple execution, package and performance
evidence remain open; the goal is active.

The same batch also passes the default-feature all-target compile in 16.42s
(`source-join-dp-controls-default-compile.log`). This checks cfg-specific build
integration only; it does not upgrade runtime or platform parity claims.

## Self-query ownership and remaining live search controls (2026-11-07)

The combined all-target `dna,mcl,zstd` check passes in 9.19s and the default
all-target check in 9.20s (`source-self-cbs-minratio-*-compile.log`). This batch
connects database-owned self queries for protein and DNA, real filtered physical
OIDs, reference seek/alias scheduling, MPP query ownership and Source's iterate
restriction. Translated self-search reaches a documented missing Source DNA
context diagnostic before query/stdin access. CBS window and motif caps now
reach ordinary, translated, global, self-score, workflow, benchmark and actual
seed-block consumers. Explicit length ratio reaches the Source setup factory
and remains fixed across sensitivity rounds before ID retention, length sorting
and real stage-1 filtering. Diagonal/Hamming controls preserve explicit-zero
presence and Source's one-context/mutual-coverage restrictions. Actual aln-out
setup and clustering threshold expression/requirements are connected, with
Source-invalid output-class lookup diagnosed only at the real HSP boundary.

The routing audit also repairs six implemented controls missing from the native
allowlist, final short-alias recognition, and vector grouping for numeric-looking
values such as `--query -1`. Source-only unused and other-workflow controls remain
explicitly classified; original registration checks still govern availability.
No runtime suite was executed, and compile success does not prove parity.

A subsequent FASTA reference-owner audit found that native seekability was
implemented by converting FASTA to a masked DMND, while Source retains an
unmasked FastaFile owner. That concrete owner correction is underway and is not
covered by the compile results above. A proposed makedb SEG replacement was
rejected after following the complete Source call graph: hard_mask=false ignores
the SEG argument and invokes TANTAN bit masking, which native makedb already
implements. The audit records this Source quirk rather than changing algorithms.

## Concrete FASTA reference owner and complete search routing classification (2026-11-07)

The FASTA owner batch passes all-target `dna,mcl,zstd` compilation in 6.96s
(`source-fasta-owner-compile-fixed.log`). The first attempt found only signed
Letter/raw-byte conversion errors in the new spool; bounded byte conversion
preserves the actual letter/mask representation. Native FASTA/FASTQ references
now retain an explicit owner, original filename, raw records, qualities,
physical OIDs/counts, independent self-query cursors and real random access.
Only offset/index metadata is retained in memory. Automatic masked-DMND
conversion has been removed; actual FASTA taxonomy/accession/MPP errors and
input rewind failures preserve the Source file type and supported input domain.
Global ranking reloads through that same owner rather than a substituted DB.

The static search-routing classification now accounts for all 316 Source
registrations: 203 common native routes and 113 explicitly classified no-effect
registrations for ordinary protein/translated search. The canonical original
parser still enforces profile/workflow permissions and argument counts. This is
routing classification, not proof of every function/command, numeric result or
output. Full callgraph review specifically confirms that Source opens FASTA
without invoking makedb, so database-build taxonomy paths have no search effect.
The audit artifact is `source-search-routing-classification.json` under the same
temporary evidence directory. No runtime suite was executed.

The FASTA owner batch also passes the default-feature all-target compile in
7.11s (`source-fasta-owner-default-compile.log`). No runtime tests were added or
run in this code-first batch.

## Master search, historical DNA dependencies, and statistics integration (2026-11-07)

Further Source callgraph audits connected real master-query preparation, one-time
length sorting and physical OIDs, header/masking order, reference iteration,
self/global-self distinctions, selected reference counts, iteration announcements,
and stop-sentinel handling. Seed annotation clearing now occurs at the actual
post-shape boundary. Dynamic seed partition widths and large arenas select the
existing generic stage0 path; the ordinary optimized path retains the Source
extension hit ordering. Global ranking's equal-key sorts now use the translated
GNU introsort movement rather than Rust stable sorting. These changes compile
with all targets and `dna,mcl,zstd` (8.49s master batch; 7.66s subsequent sort/WFA
batch), and defaults (7.95s and 8.00s respectively). Logs use the
`source-master-stage0-ref-*` and `source-sort-wfa-cluster-*` prefixes under
`.tmp/native-options-20261107/`.

The previously unresolved historical DNA dependencies were recovered read-only
from the original repository's Git history. KSW tree
`69704a761fc8160e1dfb1e191383e7565d8834c9` and WFA v2.3 tree
`3dc90a492f4a2397114e0c51bc168f1471dd398a` are preserved with blob IDs, SHA256
hashes and identical pinned DNA caller-tree evidence in
`.tmp/native-dna-20261107/recovered-source/provenance.json`. Active WFA source
comparison corrected pruned-wavefront bounds in recurrence/backtrace and the
negative-score CIGAR parser status. Canonical KSW differs only in inactive MSVC
macro definitions for the supported GNU/Apple paths. Dependency provenance is
therefore recovered; execution parity of the optional DNA profile remains open.
No original Source files, objects or executable were rebuilt or changed.

The same audit corrected prior overly broad undefined-state assumptions. Source's
static DNA ScoreMatrix metadata is zero-initialized, and default translated
contexts are defined empty sequences. Native search now retains those defined
states, uses actual DNA HSP statistics only at Source's DNA reconstruction branch,
and diagnoses invalid character access at the real access boundary. Undefined
float-to-integer and unassigned automatic-field cases remain documented safe
boundaries. Source float formatting, signed NaN spelling, fixed/default precision
and signed remainder formatting were also connected to actual output paths.
The DNA defined-state batch compiles in 10.01s; float formatting in 6.55s and
9.80s with optional/default features. These are compile results, not output parity.

Actual seed, extension, traceback, join and output statistics now fold once into
the process statistics object. Final reports reset once at the real Search::run
boundary and print after successful finalization; recursive blocks/rounds do not
reset it. Joined output uses the joiner's actual counts and avoids a second fold
in the final writer. Source's separate three-frame DP statistics remain separate.
Search temporary space records the real maximum before buffers are freed. The
statistics batch passes all-target `dna,mcl,zstd` compilation in 7.96s
(`source-statistics-compile-fixed.log`); default compilation also passes
(`source-statistics-default-compile.log`). The initial batch failed only on a
letter-count/address-sized reporting interface and was corrected before the pass.

The default-off C++ oracle now accepts `DIAMOND_FFI_SOURCE_DIR` for an isolated
source copy and follows explicit Cargo profile toggles. It never supplies a native
production fallback. Missing oracle sources leave the optional adapter unavailable;
no C++ build is required for ordinary use. Safe feature compilation uses an empty
oracle override rather than touching the adjacent user-owned Source tree.
Broad runtime testing remains deferred while these concrete translation audits
continue. Hosted Apple, complete runtime/profile/package and performance gates
remain open. Generic orchestration model tests are not evidence that production
uses DoubleIndexedBackend: the concrete native search drivers must be audited
independently, as this batch does.

### DP controls, iteration seed cache, and direct-output followup

The next independent worker batch corrected three additional concrete paths:
final traceback and full-database SWIPE now receive configured sort bins;
full-database/alternate-HSP work receives transcript allocation controls; and
alternate HSPs receive actual alignment-thread controls and the shared DP error
sink. Source evidence is `.tmp/native-statistics-20261107/dp-output-controls.txt`.
Direct native fallback render/edge/DAA paths now fold output statistics at actual
successful generation/encoding, leaving intermediate captures to the downstream
writer/join. This closes the output path omitted from the preceding batch.

The query-index/AUTO seed-cache builder now activates its precomputed union mask
for both hash-table passes, including iterations with a present all-false skip
mask. The specialized unfiltered path also applies/removes masking separately
per pass and honors an actual NONE selector. Previously the skip-mask path could
enumerate unmasked starts while the no-skip specialization accidentally masked
regardless of NONE. DAA view's metadata/selection messages and XML statistics now
use the shared Source floating-point formatters, including non-finite spellings.
The generic orchestration module documentation explicitly distinguishes its test
adapter from actual production search drivers.

All-target `dna,mcl,zstd` compilation passes in 7.22s
(`source-dp-seed-cache-compile.log`). Earlier statistics default compilation passes
in 8.26s and the safe explicit `ffi,dna,mcl,zstd` unavailable-oracle profile in
11.01s. The latter used the empty `DIAMOND_FFI_SOURCE_DIR` override and did not
build C++. No broad runtime suite was run in this followup.

The DP/seed-cache followup also passes default-feature all-target compilation in
6.61s (`source-dp-seed-cache-default-compile.log`). Goal remains active; these
checks establish compilation only.

## Viewer, frameshift, and cluster stream batch (2026-11-07)

Native DAA viewing now performs Source's synchronous first-record render before
printing the format header, then starts the actual bounded 32-record worker
queue at the second record. Empty archives avoid worker startup. Format/header
validation and metadata printing follow their real Source phases. Tabular uses
the shared bounded worker implementation; failed worker creation closes the
sender so already-started receivers can terminate. Source's literal zero query
OID in DAA view contexts was checked and retained.

Three-frame DP no longer invents gross/net cell counters: the corresponding
Source increments are commented out. Its defined zero-score admission and empty
scalar traceback now retain Source's initial maximum column and exact first-cell
coordinates, subject to the real report cutoff. Invalid/unwritten Source buffer
cases retain an explicit bounded native compatibility boundary. Proof is under
`.tmp/native-statistics-20261107/three-frame-{counter-correction,zero-score}.txt`.

Cluster realignment now writes actual per-block-pair serialized HspContext temp
streams and merges their heads with Source's heap operations instead of retaining
all member-block contexts in memory. The private codec follows basic/match.h,
including DNA feature fields, omitted query sequence/transcript fields, and
fresh deserialized contexts. Native realign/recluster forward tmpdir, file buffer,
no-unlink lifecycle, transcript reservation, alignment controls and real deferred
SWIPE errors. EXTRA-only permitted realignment coverage/approx-ID filters reach
the actual DP predicate; Source's explicit DBL_MAX realignment E-value remains.
Successful realignment reports existing statistics without resetting or adding a
search timer. Public convenience wrappers remain available.

Actual native search drivers now emit OS RSS readings at matched master/query/
iteration/reference/seed-buffer boundaries through Source's debug stream. Search
cleanup/report readings are shared at the real lifecycle boundary, including
successful MPP initialization/recovery. Real hit retrieval and completed
extension/output batches now charge their actual elapsed time; no sort timing
is fabricated when grouped storage avoids that sort. Source host CPU defaults
replace the remaining quota-aware native fallback defaults.

The translated typed roc-id routine now promotes parsed f32 identity to f64 for
Source's multiplication before signed truncation/division, preserves defined
negative bin-zero inputs, and preserves Source's uint8 family-index wrapping.
Only genuinely invalid Source conversions/array indices are diagnosed. The CLI
still returns Source's deprecation error. Evidence:
`.tmp/native-tools-20261107/rocid-numeric-audit.txt`.

The combined batch passes all-target `dna,mcl,zstd` compilation in 11.30s and
default compilation in 11.40s (`source-view-cluster-frameshift-{,default-}compile.log`
under `.tmp/native-options-20261107/`). A redundant new Read import was removed
after these checks; subsequent checks will cover that cleanup. Broad runtime
suites remain deferred during concrete translation work.

Structural audits pass for all 199 implementation-file mappings and for current
command registrations: 21 standard /45 EXTRA commands, 299 standard /316 DNA
option registrations, with no inventory/name discrepancies. The registry artifact
needed only its Rust config hash refreshed; original Source hashes/registrations
were unchanged. These audits prove inventory scope, not runtime functionality.

The backend inventory at this checkpoint finds 15 adapter traits, 12 with native
implementations. AlignBackend and DoubleIndexedBackend still have only test
adapters; Source's multinode invocation is commented out. Production alignment
behavior maps to real grouped native drivers, but that does not resolve the
adapter integration requirement. The coherent alignment integration now underway
will move actual grouped batch loading, prefetch/recycle and scoped output
execution into a real adapter, retaining owned Matches rather than a target-ID
side table. See `.tmp/native-align-driver-20261107/source-production-map.txt` and
`.tmp/native-backends-20261107/inventory.json` for the pre-refactor checkpoint.

## Real production alignment and outer owner adapters (2026-11-07)

AlignBackend is now implemented by the actual grouped native alignment backend.
The generic align_queries driver owns real preload/retrieve/prefetch, memory
accounting, batch execution/recycling and teardown. New align/native.rs owns
actual grouped vectors, bounded scoped workers, worker scratch, heartbeat and
LiveSink ordering. Actual QueryOutput owns complete Match/HSP payloads or final
bytes; no target-ID side table substitutes for them. The independent command
loading/executor loops were removed. The flat-record API remains available as
an explicit RecordAlignBackend compatibility layer. Provenance and hashes are
under `.tmp/native-align-driver-20261107/production-integration.txt` and
`integration-hashes.json`. All-target combined/default checks pass in
15.63s/16.07s. The first attempt found one migrated diagnostic needing explicit
io::Error-to-String conversion, corrected before the passing check.

Production double-indexed orchestration now lives in
run/double_indexed_owned.rs, exported through double_indexed::owned. Its concrete
associated owners and required phases replace the old independent master/query/
round/reference loops in ordinary protein, blastx, DNA and MPP search. Actual
NativeBackend, DnaBackend and MppBackend own the loaded query/reference/input,
real scoring/index/extension dependencies, dictionaries, temporary files,
alignment state, joins and output. Lazy MPP Reference owners retain exact claimed
WIP records; data/dictionaries finalize before acknowledgment. Real join tokens,
auxiliary output, stop checks and source counters remain actual owner operations.
Initialization/recovery remain their real preflight early-return paths. The old
metadata/WorkflowOp API is explicitly a compatibility model, rather than an
unimplemented production adapter. Provided-owner clustering extraction is the
next still-active integration batch, not silently counted complete here.

Real query construction/detection now precedes output creation for all these
normal paths. Missing or malformed query construction therefore cannot truncate
the output first. MPP constructs its normal output resources during begin_master;
recovery/initialization retain their distinct Source lifecycle. Corrected Source
callgraph evidence confirms MPP rejects global ranking, iteration and taxonomy
filtering before executing their unsupported combinations; SWIPE reaches actual
exhaustive extension. The audit artifact is
`.tmp/native-options-20261107/mpp-global-swipe-callgraph.json`.

The same combined batch includes MCL's pinned GNU component/member ordering,
backed by the oracle's GNU 11.4 compiler records and libstdc++ prime/insertion/
rehash policy. Physical component IDs follow Source first-visited members; an
already-remapped component entrypoint prevents applying that permutation twice.
The batch passes all-target dna,mcl,zstd in 15.59s and defaults in 15.45s
(`source-owned-search-mcl-{,default-}compile.log`; the optional passing attempt
has the `-fixed` suffix). Initial errors were one missing query-owner lifetime
and two moved query-count spellings; they were corrected before these passes.
No broad runtime suite or C++ build was run.

A further DNA audit separates precise Source failures from native safety choices:
ordinary Source DNA starts hit-buffer writers without finish_writing, schedules
zero hit partitions, and can terminate when unique_ptr reset encounters the active
writer destructor exception. Separately, Dna::Index retains a borrowed ref_buffer
that the caller frees, and its destructor later uses delete[] on that buffer.
The native owner keeps real NativeDnaBlock algorithms with safe buffer ownership;
no dummy hits or pointer faults are introduced. Source's startup invalid numeric
cast and these lifetime defects are explicit compatibility boundaries. Optional
DNA literal command parity remains unproven, alongside the outstanding profile,
runtime, package/platform and performance gates.

## Provided-owner clustering and actual accession-spool completion (2026-11-07)

Provided-owner clustering now has the fourth real double-indexed backend in
commands/cluster_double_indexed.rs. Shared owned orchestration replaces its
independent query/round/reference loops and pre-materialized range vectors.
Actual aliases, global/local self distinctions, query sorting/OIDs/quality loss,
round masks/cache/ranking, intermediates, callbacks and edges remain owned
services. Query-cache setup receives the same explicit masking/frequency/linear
overrides as the actual reference kernel. Honest physical owner letters now
accompany physical row counts from NativeBlock or the complete typed input;
global SliceReferences retains separate compact lookup/index and physical
reporting domains. Source's existing forced filtered scoring override remains
unchanged and was rechecked end to end.

Unaligned-target tracking now uses actual packed BitVector storage. A real
TempFile-backed one-String TSV accession spool is constructed/populated at the
Source owner initialization/read boundary, then Arc-shared into the tracker.
There is no eager full-ID cache during search or final selected-only live lookup.
Finalization opens its output first, parses the complete actual spool table and
then emits rows, including the all-aligned case. Source's empty BLAST mapping
and absent provided-owner mapping remain distinct explicit safe boundaries.
Concrete buffer-aware file constructors and StreamingFile.from_input provide the
real descriptor/TSV services; ordinary typed defaults remain unchanged. This
supersedes the intermediate live-lookup malformed-input timing limitation in the
earlier tracker artifact. Proof is under
`.tmp/native-unaligned-targets-20261107/spool-provenance.txt`.

The spill audit confirmed literal X titles, discarded qualities and separate
original OIDs were already present. It corrected raw annotation-bit retention to
match Source's SEQ_MASK FASTA print/read roundtrip. seqid now uses the existing
Source byte scanner, including C-string NUL termination.

All-target combined dna,mcl,zstd compilation passes in 10.42s and defaults in
10.55s (`source-cluster-owner-accessions-{,default-}compile.log`). A redundant
mutable self binding was removed before the default pass. No broad runtime tests
or C++ builds were performed during these translation batches.

The refreshed qualified backend inventory identifies 14 active production
contracts with native implementations, two explicit compatibility models and
one upstream-inactive multinode invocation. The owned double-indexed contract
has four actual driver implementations. This closes the known disconnected
orchestration work in the code inventory; counts and compilation do not prove
runtime parity. Fresh runtime/profile/package/platform/performance checks now
remain the next completion phase, and any behavior gaps they reveal still
require translation fixes. The goal remains active and unachieved.


### Runtime-discovered translation fixes and Source numerical dependencies (2026-11-07)

After the known owner/orchestrator inventory closed, the first default library
batch passed 1772 cases with 15 failures and 6 ignored cases. Batched corrections
to stale Source-incompatible fixtures and native code reduced this to 1785 passed,
2 failed and 6 ignored. The integration batch passed all 20 upstream built-in
regression hashes. Eight integration targets still failed; these are recorded in
`.tmp/native-post-translation-20261107/default-integration-targets.log` and must
not be represented as a passing suite.

Actual code defects identified in this phase include shared self-search seed
ownership, hardmasking during a soft-only clustering pass, clearing an already
computed master query self-score table, and CLI stdout ownership/stdio tie
semantics. The corresponding native fixes are being batched before rebuilding.
Fixture corrections preserve exact biological output and Source stream bytes;
ANSI sequences are not indiscriminately removed. Source probes distinguish a
buffered reset suffix on help output from the reset discarded when a POSIX
OutputFile closes descriptor 1. Embedded library calls retain their host stdout.

The dense fractional MCL path now translates the recovered Eigen 3.3.8
Hessenberg/Francis-Schur eigensolver, complex arithmetic and partial-pivot LU,
products and reductions. Sparse duplicate ordering and tiny-value pruning also
follow Source. A standalone current-source replay matches all 54 exact numerical
fixtures, including complex fractional expansion through size 17. This result
covers numerical stages, not complete MCL CLI parity. The exact pinned full MCL
CLI cannot compile because its caller still uses removed transport APIs; an
isolated numerical harness invokes unchanged Source matrix functions instead.
All C++ work used /tmp copies and left the adjacent user-owned tree unchanged.

These results expose and fix code gaps; they do not complete the goal. Shared
Cargo integration of this latest batch, the remaining runtime failures, optional
profiles, package/platform execution and performance gates remain open.


The follow-up code batch also corrects actual MCL checkpoint component/block
ordering, independently of the later work-size sort. Realign/reassign/recluster
and wrapper output adapters now use Source OutputFile ownership; the wrapper
opens output before either search branch. GVC retains its Source ofstream/no-
output semantics and MCL retains its Source cout/ofstream row flushing.
All-target native dna,mcl,zstd compilation passes in 8.81s, recorded in
`.tmp/native-post-translation-20261107/current-batch-combined-check.log`.
The current frozen code hashes are in current-batch-source-hashes.json; default
library and affected integration targets are being checked as one batch.


The affected-target rerun is complete: default library 1786 passed / 1 failed /
6 ignored; GVC 2/2, MCL numerical 1/1 (54 exact cases), realign 3/3, reassign
3/3, search output 5/5, tabular context 5/5 and view stdout 2/2 pass. Recluster
has 2 passes and one recursive representative mismatch. Logs are in
current-batch-targeted-tests.log. This is a selected batch, not a complete suite.

Source debugging identifies two further actual stage-0/stage-2 integration gaps:
clustering's explicit softmask override was omitted from the authoritative seed
Block tables, and generic borrowed partition filtering lacks Source's same-
physical-record self-hit rejection before tentative-match counters. The initial
softmask table correction compiles with dna,mcl,zstd but the focused LinClust
execution exposes a histogram/enumeration mask-state assertion; this remains an
open code defect, not a relaxed expected output. A coherent follow-up correction
is underway. The recursive representative expectation remains unchanged too.
No production fallback has been reintroduced and the full goal stays active.


### Seed annotations and actual generic stage-2 self filtering (2026-11-07)

Source Reduction::reduce_seq reads biological letters through Sequence's
annotation-stripping access. Native reduction previously indexed raw annotated
bytes, producing different partition keys between the histogram and fill phases
once real softmask tables were connected. Reduction now masks annotations while
retaining biological X/stop values. The actual borrowed partition callback also
rejects same physical-record self hits after the score predicate and before
Source tentative-match counters, using the actual first-reference-block self
context. These are production-path corrections, not relaxed histogram or output
assertions.

The combined dna,mcl,zstd all-target build passes in 7.42s. All eight cluster
command unit tests pass, including both previously failing biological partitions;
the new reduction invariant and actual native regression logging test also pass.
The recursive integration sweep now passes its original representative mismatch
and later reference variants, but still fails at natural approx90 p1/16G after
58.86s; this remains an open parity/code investigation. Exact commands, statuses
and timings are in masking-reduction-focused-results.json under the current
post-translation checkpoint. No broad green suite was repeated.

Package-list inspection confirms the new numerical modules, frozen vectors and
startup fixture helper are packaged, and no adjacent C++ tree is included. This
is list evidence only. Isolated-package CI now additionally executes native
features dna,mcl,zstd on Linux and both explicit macOS architectures, preserving
default tests and the all-feature absent-oracle build. Hosted results remain
unverified; workflow changes alone do not prove platform parity.


### Existing clustering compatibility state reconnected (2026-11-07)

The remaining natural approx90 mismatch traced to a concrete owner-state loss:
the cascaded facade already supplies an explicit reference-counter compatibility
value for Source's uninitialized first-read/sort boundary, but the newly shared
provided-owner driver discarded it. The driver now receives that existing value
through an initial_reference_counter hook, defaulting to zero for ordinary/DNA/
MPP adapters, then continues retaining actual traversal counters. No fresh
round heuristic or altered biological expectation was introduced. The Source
undefined-state qualification remains explicit; this is diagnostic-profile
compatibility, not proof of a defined upstream first-counter value.

Combined dna,mcl,zstd all-target compilation passes in 10.17s. All three
recursive-clustering integration tests now pass, with the full matrix completing
in 82.29s (`provided-initial-reference-{combined-check,recluster-test}.log`).
The earlier default integration failures are now corrected, including the
separately passing actual regression-log fixture; no complete rerun is claimed.

An isolated packaged snapshot builds both native dna,mcl,zstd and exact
EXTRA-only profiles without adjacent C++ sources. The archive hash, executable
hashes and exact staging paths are recorded in native-package-location.json.
This snapshot predates the final reference-counter hook. A 44-case optional
Source tool capture/comparison (37 defined comparisons) exposed concrete cut
stdout-owner and file-open-diagnostic gaps, now under correction. Signal and
nontermination boundaries are classified separately, not silently passed.


### Actual file-error transport and temporary-owner lifecycle (2026-11-07)

The Source FileSource/FileSink openers print perror before throwing their
FileOpenException; main catches that class without another message. Native
openers now capture typed path/OS diagnostics, preserve them through actual
String-returning adapters via invocation-owned IoErrorOwner and through the
cloneable RunError bridge, and report the dedicated protocol once. This is
actual operation metadata, without message matching or an extra open probe.
Makedb final-header rewriting retains the original descriptor. Source-specific
GVC inputs now use InputFile codecs; validation precedes opening, and its
ofstream output failure remains silent as Source defines. Dbinfo/makeidx
receive their actual EXTRA file-buffer consumer; makeidx input extension
fallback follows the existing Source helper.

Source Search::run always evaluates TempFile::get_temp_dir even when quiet.
That genuine temporary owner is now connected before statistics reset, within
total timing. Explicit tmpdir wins, otherwise Source's output directory is used,
otherwise the empty/current directory. Existing query/reference, accession,
intermediate and graph owners use the policy rather than invented probes.
Unix uses real mkstemp and captured temporary-file diagnostics. Windows generates
a name and defers physical creation to its actual constructors, preserving the
FileOpen error class there. Context streams retain the original Unix descriptor
and their actual create/write/flush/merge errors. Original clustering accession
metadata uses a real constructor-time TSV and final typed table. Its final
explicit deletion is independent of --no-unlink, as Source's Tsv::File close
requires; a native paired case confirms identical rows and an empty tempdir.

The coherent batch compiles for default and dna,mcl,zstd, including all targets
on x86_64-pc-windows-gnu (36.19s); this is compilation, not Windows runtime proof.
Paired checks pass all 20 standard missing-input cases and all six GVC codec/
ordering/output-failure cases byte for byte. All 14 output-open cases pass:
ten exact streams and four with only the precisely validated six-character
mkstemp suffix differing. Of 37 defined EXTRA captures, 35 are byte-exact; the
other two differ solely in Source Release versus Rust Debug assertions. Signals
and nontermination are separate classifications, not passing parity cases.
Evidence is under standard-file-errors/, extra-tools-source/, gvc-codecs/,
output-open-errors-fixed/ and accession-no-unlink-{source,native}.json.

### Optional CBS validation and physical DMND record boundaries (2026-11-07)

Optional-profile library execution identified an actual CBS count error: Source's
enum explicitly reassigns 5 after 6, so its implicit COUNT is 6. Both profiles
reject CLI modes 6/7. Three native guards now agree; all six bounded native
comparisons against pinned baseline/EXTRA capture identical diagnostics and
statuses. Parser fixtures now exercise EXTRA's permitted workflow bypass while
retaining strict arity and testing target-indexed constraints at actual startup.

The exact Source header harness confirms baseline Header2 sizeof/payload 48/48
and EXTRA64 56/52. Serialization prefixes sizeof while writing only the fields.
Tests now distinguish this from the physical sequence boundary. Actual bulk and
streamed positional readers use the first SeqInfo.pos instead of treating the
advertised record size as physical payload length. The Source-specific accession
API separately retains its literal fixed init_seq_access offset; indexed reads
remain independent. An earlier suspicion that Source int reads omit byte
accounting was disproved by the exact harness: native accounting was already
correct and no such production change remains. Evidence is in
header2-int-proof/provenance.json and feature-parser-audit/.

After the corrected batch, default library tests pass 1792 cases with six ignored;
native dna,mcl,zstd library tests pass 1791 with six ignored. Combined/default and
Windows GNU all-target checks pass; the latest exact commands/statuses are in
source-feature-count-header-results.json. Current release compilation is now
underway for fresh core parity/performance checks. The goal remains active;
complete integration/profile/platform/package/performance evidence is still
required.

The fresh exact pinned baseline C++ oracle independently passes all 20 built-in
regression cases, status 0 in 26.47s. Its auxiliary queue tests pass and disk
HitBuffer checksum passes; in-memory auxiliary stress reports bad_alloc under
the explicitly imposed 4 GiB address-space cap. This qualification is retained
in baseline-builtin-test/results.json. All isolated source hashes and adjacent
Source witnesses remain unchanged.


### Sequential getseq, trace-free SWIPE and real reservation policy (2026-11-07)

Fresh matching EXTRA databases are byte-identical between Source and native.
Roundtrip checks exposed two further defined Source flows. Getseq reads from the
actual post-header descriptor cursor rather than the positional API; because
EXTRA64 advertises 56 header bytes but serializes 52, this affects its first
record. CLI getseq now uses the real sequential decoder/cursor, without
truncating indexed output or padding serialized files. Indexed library services
remain independent. Source/native getseq (including selected, reversed,
hardmasked and chunked output), dbinfo and SWIPE on both databases match exactly
in extra-dmnd-roundtrip/verified-physical-reader/proof.json.

The coordinate-only default-sensitivity SWIPE optimization declares COORDS but
actually runs NONE and copies initialized DummyRowCounter HSPs. Source's query
interval is [0,1], with a real subject endpoint; adding statistics/transcripts,
very-sensitive or top-percent changes selects real coordinates. Native now
preserves the actual declaration and the trace-free DP specialization, including
Source coverage filtering, with no output rewriting. Actual coordinate/stat/
transcript paths retain their real ranges. Independent baseline and EXTRA
captures cover all six branches (12 Source runs); the new unit covers a long
query without dummy row-counter limits and a target endpoint before target end.

A native release run under the same 4 GiB bound as Source then aborted on a
999999990-byte allocation before its first regression case. The automatic
16 GB soft RSS ceiling had incorrectly become an eager reservation in every
hit bin. Automatic mode now grows with actual hits and retains actual
RSS-triggered disk migration. Explicit --hit-membuf connects Source's fallible
1 GiB-per-bin reservation and persistent memory selection; native 0G still
forces disk. A typed BadAllocError preserves the actual allocation class through
the runtime bridge without mistaking an OS ENOMEM for std::bad_alloc.

Default library tests pass 1795 with six ignored; dna,mcl,zstd pass 1794 with six
ignored, and combined/Windows GNU checks and default release compilation pass.
The corrected native release completes all 20 built-in cases with status 0 in
21.326s under the same 4 GiB cap and leaves no files. Its auxiliary memory stress
passes where Source's explicitly bounded memory reservation fails, so these
whole-command times are not an equal-work performance comparison. The Source
and native core regression hashes pass independently. Raw proof is under
native-release-builtin-test-fixed/; the prior native abort is retained under
native-release-builtin-test/ rather than erased.

Explicit --hit-membuf failure also matches Source quiet streams/status without
aborting. The --log variant identified terminal/log ordering: Source emits the
terminal catch diagnostic before its log while streams are still live. Native
now reports at that actual owner boundary, retains the original typed error plus
reported state, and avoids a duplicate outer diagnostic. Generic adapters keep
their supplied-stream reporting. This latest narrow ordering correction is in
shared verification; do not extrapolate earlier evidence to its result.


### Current stable parity checkpoint (2026-11-07)

The terminal/log ownership correction compiles for combined native features and
Windows GNU; the four generic-main and two Rust-FFI-facade tests pass. Actual
NativeBackend reports the Source catch diagnostic before its still-live log,
then marks the typed cause as already reported so outer CLI/library frontends
do not duplicate it. Explicit source-memory allocation failure now matches the
quiet Source streams exactly and the logged allocation hint/error ordering with
one logfile error. A further Source ostream evaluation detail was retained:
`Temporary directory: ` is inserted before the real resolver, with the suffix
only after success. Its logged and quiet failure probes pass the complete
protocol with only the exact six-character temporary nonce qualified.

The final paired regressions pass all 12 SWIPE branches (both Source profiles),
all 20 standard input-open errors and all six GVC cases. Source/native EXTRA
makedb bytes match (525 bytes); both .seed_idx files also match (1088 bytes with
two 512-byte shape tables). Target-indexed results match for both executables
using either index. Static Source tracing confirms that indexing uses SeqInfo
positions, unlike getseq's direct post-header sequential path. Proof paths:
extra-dmnd-roundtrip/swipe-dummy-proof/native-comparison.json,
extra-dmnd-roundtrip/verified-physical-reader/proof.json,
extra-dmnd-roundtrip/verified-index/proof.json and
output-open-errors-prefix-fixed/results.json.

All source workers are frozen; no live parent Cargo sessions remain. Current
code hashes are current-final-code-hashes.json. The 199-file inventory and
316-name registry checks pass. Canonical registry regeneration changed only the
Rust Config hash; original Source hashes and registrations are unchanged.
No current all-target integration/package/hosted-OS/performance completion claim
is made from these selected results. The goal remains active and unachieved.


### External clustering owner and BLAST metadata follow-up (2026-11-08)

The active external `linclust --parallel-tmpdir` protocol now uses Source's
versioned job directory. A real two-process comparison found further translation
gaps: premature pair-volume publication, an extra closure-owner guard, missing
pre-reverse SWIPE filters, and different chunk-owner close order. These have
been corrected from the concrete Source call graph. The subsequent native
one-worker control matches all 9,704 mapping rows and both diagnostic streams.
Its 26,946 alignment rows have identical values; their raw serialization exposed
a further six-significant-digit formatting gap, now corrected in the writer.
The first corrected snapshot's raw alignment comparison remains recorded as a
failure; numeric diagnostic normalization is not a wire-parity pass.

The corrected two-worker run still exposes an initial empty-chunk race: one
worker sees 12 chunks, while a late worker unconditionally allocates a thirteenth
empty chunk and waits for its completion. The same unconditional allocation and
polling order exist in pinned Source. The preserved Source two-worker run also
has schedule-dependent late closure behavior, producing a different mapping
from its own one-worker run. These are defined file-protocol scheduling hazards;
no new barrier, artificial delay or lazy allocation has been substituted for
Source behavior. The actual native two-worker gate remains unresolved.

A further concrete owner audit found that native FileArray storage lazily
opened volumes, bypassed Source OutputFile staging, and flushed at uncompressed
rather than compressed byte counts. The current correction eagerly retains
volume and FileStack owners, uses the actual 64-KiB external stream buffer,
restores rollover/final-close order, and logs Source's radix claims. Focused
verification of this new batch is pending. The optional CompressedBuffer profile
also requires Source's actual Zstd stream rather than unconditional zlib.

The BLAST constructor now deduplicates taxids across a sequence's deflines,
matching RawChunk's set, and exposes the pinned GNU reference's equivalent-key
iteration order. Explicit mapping additions still retain duplicates. The
multi-defline fixture and compiler/header provenance are recorded in
`blastdb-defined-mapping-audit/proof.json`; runtime verification is pending.
Detailed external evidence is in `external-two-worker/` under
`.tmp/native-post-translation-20261107/`. These changes do not establish complete
external, BLAST, optional-profile or hosted-platform parity. The goal is active.


The combined follow-up batch is now implemented and frozen for focused checks.
Clustering assignments also use the actual FileArray/BufferArray owner instead
of a separate hardcoded zlib encoder. Optional internal streams select Zstd
through the native codec. The actual BLAST CLI accession filter now uses
Source-profile byte-line input, including gzip decoding, while TAXIDLIST keeps
its distinct raw-file semantics. Typed input errors and arbitrary missing
accession bytes reach the terminal and log unchanged, including Source's
`what()` NUL terminator behavior; established String-returning library entry
points remain available. Gzip/byte/error fixtures and a raw diagnostic boundary
test accompany this code. No new runtime pass is claimed before the combined
checks complete.

The two-worker source-order replay reached its original bound and was reaped
with statuses `[-15, 0]` after 602 seconds. The timeout and strict mapping and
alignment mismatches remain preserved in
`external-two-worker/comparison-source-order-two-worker.json`.
Actual Source compressed-buffer probes are captured once for both profiles:
the default probe links the unchanged Source object, and the Zstd probe compiles
only byte-identical isolated copies of the wrapper against existing libzstd
1.4.8. Native comparison against these captures is pending; the native codec
backend's different version is explicitly recorded.

Direct comparison of the actual compressed-buffer APIs found a remaining gap in
that frozen package despite its focused tests passing. All default final bytes
match, but 12 intermediate snapshots differ because the zlib convenience
encoder buffers output differently from Source's direct deflate calls. This
changes the size used by BufferArray's flush threshold. The corrective port now
needs direct streaming calls and Source's retained output-buffer growth policy.
The optional profile matches 33 of 34 byte snapshots and every size index; its
empty-stream frame header differs between library implementations. That last
case is under a separate protocol audit. The first comparisons and failure
artifacts remain in `compressed-buffer-profile/native-proof.json`; they are not
reported as full profile parity. The completed 5a package build does not close
these newer findings.

The empty-Zstd protocol audit rules out a C-library version difference: both
1.4.8 and 1.5.7 emit the same canonical empty frame. The pure-Rust endStream
wrapper bypasses generic finalization's implicit zero-size pledge. The
correction applies that actual encoder parameter transition to never-fed
CompressedBuffer and ZstdSink streams. Default compression now uses direct
streaming calls with Source's buffer growth and reuse policy. Both compile
profiles and the focused default codec/empty-sink checks pass in package
`34ce17f2…`; strict recovery of all Source snapshots and the fresh runtime
controls are pending. Earlier failed captures remain unchanged.

Strict recovery now passes all 68 original Source compressed-buffer snapshots:
34/34 in each profile, requiring exact intermediate prefix bytes, stage sizes
and final bytes. The correction uses direct streaming calls and the encoder's
actual empty-size parameter transition. Evidence is in
`compressed-buffer-profile/native-corrective-proof.json` against package
`34ce17f2…`; the earlier 22/34 and 33/34 failures remain unchanged. No Source
capture was regenerated. Both compile profiles and the 12 focused codec checks
pass; fresh external CLI controls remain pending. This closes these specific
wrapper defects and does not claim complete external or optional-profile parity.

The corrective optimized binary (`38baedec…`, package `34ce17f2…`) now passes the
actual one-worker external CLI control strictly: all 9,704 mapping rows, the
26,946 raw alignment rows as a multiset, and both diagnostic streams match the
preserved Source run. No numeric normalization is used. Volume claims are
unique, chunk counts match and there is one output winner. See
`external-two-worker/comparison-prefix-one-worker.json`. The new two-worker
control is running against the unchanged original Source capture; its result
is still open. An additional invocation-owned external I/O error transport
audit is in progress. The full goal remains active.

The external I/O audit confirmed another concrete gap with a physical directory
in place of the first pair-table output file. Both executables fail immediately,
but Source emits its dedicated FileOpen diagnostic with the saved OS reason;
the current native adapter flattens that cause into a generic String error.
The before capture is `external-io-owner-audit/proof-before.json`. Invocation-owned
transport is now being connected to the existing outer error owner, preserving
public String APIs and Source's distinct constructor classes. A separate reader
ordering correction retains the actual input owner before creating accession
output. This new batch is in progress and has not been compiled; earlier success
proofs do not establish its error-path parity.

The current error-path audit also found incomplete default accession-output
integration. The OID-output control above does not cover this path. Native
accession output skips Source's final non-POD AccMapping radix stage and opens
one file per representative; Source retains one output stream per source volume
inside each radix bucket. That changes file names and contents when groups
share an owner. The missing radix stage and exact owner scope are being ported;
this is required native integration, not an evidence-only gap. Source's unchecked
ofstream and ignored final OID fprintf/fclose failures are also being retained.
The first shared I/O compile caught a missing Read import; its failed package
is preserved, and the expanded batch has not yet been compiled.


The expanded external batch now passes both all-target compile profiles and
45 focused external tests. Four actual I/O recoveries match the original Source
captures exactly: obstructed output, missing input, empty input and invalid first
byte. Status, stdout, stderr, and every artifact name and defined byte match,
including job logs and removal of 512 empty volume files. Evidence:
`external-io-owner-audit/proof-current.json`.

Default accession output is now concretely integrated: streamed byte lines,
variable-length NUL-terminated records, Source's final non-POD radix queues and
single eight-bit spill pass, the pinned GNU 96-byte record charge, and one output
stream per volume inside a radix bucket. Four output-API controls match payloads,
statuses, group counts and phase counters, including 65,536 records whose output
changes from one file to two at the Source charged-memory boundary, and ignored
final OID write failure. See `external-accession-output/proof-current.json`.

These accession controls have an explicit Source toolchain qualification. The
original GNU 11 object aborts in `std::pmr::list::sort`; a minimal allocator-backed
list reproduces the abort with GNU 11 and succeeds with the installed GNU 12
headers. The working comparison compiles only byte-identical `output.cpp` with
GNU 12; every other baseline Source object and original build remains unchanged.
Its object/probe/CLI identities, headers and the reproducer are recorded in
`external-accession-output/source-gcc12-build.json`. Original GNU 11 failures
remain preserved and are not relabelled as runtime parity passes.

The complete natural default non-OID external CLI also matches this qualified
working Source profile: all 170 output filenames and every byte match, covering
9,704 rows; all 26,946 raw alignment rows, streams, activity totals and queue
counters match. Source/native exited zero in 43.02/23.04 seconds in this single
shared-host run. Evidence: `external-two-worker/comparison-accession-cli.json`.
This is a new default-mode workflow, not a retry of the two-worker comparison.

The original two-worker gate still fails at Source's shared late empty-chunk
hazard. Separate open work remains in FileSink's C stdio buffering, saved write
and close diagnostics, and Source's fatal noexcept-destructor behavior. Native
best-effort Drop is not claimed equivalent to Source abort on cleanup failure.
The Windows maintained bindings also need a separate CRT compatibility plan;
none of these remaining requirements is closed by the selected successful cases.
The goal remains active.


### 2026-11-08 POSIX output stdio and fatal cleanup checkpoint

Frozen package `b2f8e941…` contains actual libc FILE ownership/buffering rather
than raw Rust writes, with saved operation-specific diagnostics through the
native CLI and library bridges. Both native all-target compile profiles and
160 focused I/O/database-builder/external checks pass. Its immutable default
release is `882623719d04237e…`; later input/Windows edits are a separate batch.

Three direct FileSink controls match unchanged Source objects exactly: a tiny
`/dev/full` write succeeds then close fails; a larger write fails but close
succeeds; ordinary tell/physical-size/seek/rewind/final-byte behavior matches.
The API driver explicitly emits the saved diagnostics and does not establish
CLI reporting by itself. Separately, actual `getseq -o /dev/full --log` matches
Source's ordinary close failure and log error line; CLI stderr comparison
substitutes only the exact executable invocation line, with full streams kept.

Actual external CLI controls verify the Source destructor policy: both normal
cleanup failure and missing-input followed by cleanup failure abort with
SIGABRT. Status, stdout, stderr and every artifact match Source exactly
(1,547/1,545 entries). The primary open perror precedes the cleanup ENOSPC
perror, and the ordinary main catch is never reached. Source uses its actual
empty compressed member during unwind; the harness invents no writes.
A genuinely readerless stdout pipe also gives the same SIGPIPE termination
as Source. These new Native recoveries do not rerun the preserved Source cases.

Evidence: `.tmp/native-post-translation-20261107/stdio-source-plan/controls/`
(`native-current-linked1/proof.json`, `cli-getseq-native-proof.json`,
`broken-pipe-current.json`) and `fatal-cli/native-comparison.json`.

The next translation batch addresses actual fread/feof/rewind/close input
behavior and shared FILE state across output/input aliases. New unchanged-object
Source controls capture directory read errors, nonseekable stdin and shared
EOF/cursor effects under `stdio-source-plan/input-controls/`. Windows CRT
bindings/ownership are being implemented separately; no Windows runtime parity
is claimed. The two-worker scheduling gate and other open completion evidence
remain open. The full goal is active.


### 2026-11-08 retained FASTA, spool and console followthrough

The retained FASTA class now reads records from actual InputFile owners rather
than caching complete files. Metadata flags control its count/index pass;
provided mapping callbacks run during that constructor pass, including before a
later record fails. NativeDatabase uses this class for its CLI FASTA path and
reuses the resulting counts/mapping pass. The previous whole-record Memory path
remains a Rust extension only for explicitly supplied records.

Eleven direct public API controls match preserved unchanged Source objects
exactly in status/stdout/stderr, including paired FASTQ's consuming sync probe,
index/count/cursor behavior, permissive FASTQ quality fields, and named writer
flush/deletion boundaries. Source raw FASTA tokenization leaves TextInputFile's
EOF flag false; the first frozen native snapshot matched only 7/11 cases, then
the literal EOF correction recovered all eleven without rerunning Source.
Evidence: `fasta-class-source-controls/native-7655-before/proof.json` and
`native-c27-eof-recovery/proof.json`. Counts without Source's metadata flag are
indeterminate upstream; Rust exposes zero as an explicit safe unavailable state.

The spool now creates and retains both the actual TempFile writer and input
alias at construction. Destruction closes the input alias, discarding pending
outer serializer data exactly as Source TSV File does. Explicit reading rewinds
the writer then the same input. Source-mode close failure preserves the primary
perror before cleanup perror and abort, without proceeding to deletion. Shared
IO error transport and provided backend hooks retain actual typed errors before
constructor, search-loop and consuming-finalizer owners drop. Generic library
construction remains a distinct nonfatal extension.

Actual system FILE bindings now carry console output on GNU/musl/Darwin and
Windows. Cout/cerr failbits, tied stdout flush, stderr unit-buffering and untied
perror remain distinct from C printf. Generic caller-supplied writers retain
their existing contract. Proven command adapters and message-stream call order
are connected; further low-level direct output sites remain qualified.

Frozen package `7655ac01…` passes three all-target compile profiles, 191 focused
units and nine cluster adapter integration tests. The subsequent `c27b8b29…`
package passes the same compile profiles plus 194 focused units (the empty
owned-driver filter is not a coverage claim); its default debug library supplies
the eleven exact FASTA recoveries. The subsequent `34a90749…` package adds wide
Windows filename APIs/stat-layout assertions, binary-only inherited-SIGPIPE
capture and message routing/order corrections. It also passes the three compile
profiles and focused checks. These are default/optional Linux and Windows GNU
cross-compilation, not hosted Windows runtime or full Apple crate validation.

An actual pre-fix Darwin module reproduces E0617 for both Apple architectures;
the corrected C variadic mode promotion passes both targeted module checks.
Evidence: `stdio-source-plan/apple-primitive-check/proof.json`.

The executable now preserves its inherited SIGPIPE action across Rust startup;
the library installs no initializer. Native-only recovery matches both original
Source controls: SIG_DFL gives signal 13 with empty stderr; inherited SIG_IGN
returns zero with empty stderr. The previous unconditional default restoration
failed the ignored-parent control. Evidence:
`stdio-source-plan/inherited-sigpipe/proof-current.json`.

Wide Windows CRT opens/stat preserve mode/ownership and avoid UTF-8-as-ACP
filename corruption. Non-ACP-representable Unicode paths are an explicit Rust
extension; Source ACP argv/diagnostic-byte equivalence, actual Windows runtime
and MSVC ARM64 header/link validation remain unproven. Bare FastaFile base-mapping
ownership outside provided adapters, full raw-byte title propagation, remaining
console consumers, the two-worker scheduling gate and final-snapshot broader
workflow/performance checks remain open. These checkpoints do not establish
full translation completion.


### Raw byte title ownership batch (2026-11-08; verification continues)

Native FASTA/FASTQ records, DMND/BLAST deflines, search dictionaries and query
views, DAA records, output renderers, auxiliary sequence output, clustering
callbacks and EXTRA tools now use byte title/quality cores. Checked String
convenience APIs remain; Source C-string truncation is applied locally at its
actual consumers. This includes raw accession joins and persisted MPP
reference dictionaries. Alignment and taxonomy numeric algorithms are unchanged.

Frozen isolated package `5a0c1e993b6b1ae19f2e03fe969fa62f980207c3e35e451af4b7e47d98ea187b`
contains 1,954 files, exact copied inputs and no adjacent Source tree. Linux
all-target default/optional (`dna,mcl,zstd`) and Windows GNU optional cross-checks
pass. Focused filters report FASTA 28 passed, DAA 49 passed and BLAST 48 passed /
4 ignored. The byte filter reports 61 passed / 1 failed: its new taxonomy fixture
expected LEB encoding `[1,7]`, while the unchanged production serializer and
pinned Source encode singleton taxid 7 as `[3,15]` (`(x<<1)|1`). The workspace
fixture is corrected with an additional decode assertion; the original failure
remains recorded and is not relabeled as a pass. These filters overlap.

The same immutable package builds the default binary/library. Binary SHA-256 is
`d85215ce9e88186d12ce362a2a0e53f240e0ea1d879297dec79e3ccf1c5fe544`.
Evidence: `.tmp/native-post-translation-20261107/raw-title-integration-corrective-package/`.
Actual upstream raw-title CLI comparisons are underway. A separate subsequent
batch is completing Block library byte APIs and taxonomy scientific-name byte
ownership/rendering; it is not covered by this package. Neither these focused
checks nor the title implementation establishes full translation completion.


The first actual raw-title CLI control captured 20 commands per implementation:
18 match status, stdout, stderr and every output artifact exactly, including
complete successful DMND/DAA files, DAA merge and raw FASTQ auxiliary output.
Two stdout comparisons differed. The defined embedded-NUL query-title case
exposed missing C-string boundaries in native context construction/rendering;
those fixes are in the next package, awaiting native-only comparison against
the retained Source capture. The non-ASCII XML case reaches Source undefined
behavior: `escape_sequences.h` allocates `seq_[UCHAR_MAX][7]` (255 rows) and
indexes it with `(size_t)char`. In the verified GCC 11 signed-char profile,
0x80–0xff index outside the table; even unsigned-char 0xff exceeds it. Rust
keeps a safe 256-entry identity/escape table. Its differing bytes are an explicit
safe extension for this undefined Source domain, not a parity pass or an
emulation of accidental reads. Exact compiler/header/offset evidence and all
original captures are retained under `raw-title-cli-controls/`, including
`xml-source-ub-evidence.json`. Nine distinct ASCII-only NUL format controls have
also been captured once from Source, isolating its defined boundaries. No
original Source case was retried or normalized.


The next isolated package, `bfd84cccd22f99b4b29267cf7d91b3226f76a5f9a8fc9c56cc7b440cff119fed`,
passes default/optional Linux and Windows GNU all-target checks. Its focused
filters pass: byte 66, taxonomy 14, database builder 9, Block 10, taxdmp 6 and
DAA 49 (overlapping filters, not a unique-test total). The corrected varint
fixture passes here. Raw scientific names/ranks now retain bytes through Source
text readers, DMND serialization, tree metadata and actual taxonomy output;
Block byte views/dictionaries and HSP C-string construction boundaries are
connected. Evidence is in `raw-taxonomy-library-package/`.

Immutable binary `ab47e5202a05e022d179cc5fd4aa41fd16ca24007a76433e99607ce2d86dc3c7`
then matches all 10 retained NUL controls exactly, including all nine format
controls. Source was not rerun. Five new taxonomy controls also match: four
successful commands compare every byte, including the complete 500-byte DMND
SHA-256 `f0dcfd26a7ea5a36380766c9cc7367ad9d151bde90d7fd977618b8fb685f2365`;
the invalid raw-rank control matches status 1 and raw stderr exactly, with only
its partial file's Source-indeterminate `pos_array_offset` field `[32,40)`
qualified. Successful files exclude nothing. Evidence:
`raw-title-cli-controls/native-bfd84ccc-nul-recovery/proof.json` and
`raw-taxonomy-cli-controls/first-bfd84ccc/proof.json`.

A subsequent code audit identified a remaining UTF-8 title gate in the actual
frameshift capture backend (`align/legacy/native.rs` and its QueryMapper
contracts). That byte-contract migration is in progress and is not covered by
the package above. Full translation and final-snapshot parity remain open.


## Native backend and parallel storage continuation (2026-11-08)

The retained multinode library now has a concrete native search/GVC backend,
with real indexed FASTA owners, retained initialized configuration and one typed
I/O error owner. Source's commented CLI call remains inactive. External radix
spill and POD chunk inputs now use real parallel readers and shared output
owners. The native intermediate writer calls the translated
`generate_intermediate_output` with the actual dictionary and original OIDs;
the duplicate framing loop is removed. Seed progress uses Source's numeric
formatter. These edits await the combined frozen-package checks below.

The oversized POD spill path has an explicit Source undefined-behavior boundary:
`diamond/src/cluster/external/radix_sort.h:38` calls
`in.read(&x, sizeof(T))`; inherited `Deserializer::read(T*, size_t)` at
`diamond/src/util/io/deserializer.h:118` interprets the second argument as an
element count and requests `sizeof(T)^2` bytes into one stack object. Native
code reads at most one record into its retained buffer, propagating real I/O
errors and preserving the previous bytes on short reads/EOF. All four actual
record types have zeroing Source constructors, so initial buffer state is also
retained. Source's volume workers and file-array ownership are preserved. No stack overwrite is emulated and no byte-parity claim
applies to that undefined Source domain. Defined paths and final-snapshot
upstream comparisons still require evidence; this does not complete the goal.

A read-only dispatch audit again finds no automatic production C++ call. The
only executable oracle branch requires explicit `--legacy`, the default-off
`ffi` feature, linked oracle availability and a non-Windows target. Public
native dispatch never enters that branch. Evidence and input hashes are retained
under `.tmp/translation-native-backends-spill-20261108/dispatch-audit.json`.


Native backend/storage verification (2026-11-08): frozen package
`e95142bf7ac1de960cef6e6ca17d5b9bd3f8a219aefa1982050b714b987d1f64`
passes isolated default, optional dna/mcl/zstd and Windows GNU all-target checks.
It also passes 112 focused tests: 69 external, 22 multinode, eight alignment
output, ten intermediate and three native search-block integration cases. The
new multinode integration uses real startup, raw-title FASTA index/read owners,
BLASTP and GVC; it checks exact assignment bytes and the shared initialized
matrix owner. The error control verifies the exact typed Source index diagnostic.
Parallel loader, spill routing, EOF reuse, rollover counters and error-owner
controls pass. Source was not run or rebuilt. All packaged bytes stayed frozen.

A final narrowly scoped follow-up uses the existing Source fixed-float formatter
for multinode startup values, preserving signed/lowercase nonfinite spellings
accepted by Source atof. Its isolated default all-target and multinode/formatter
checks are recorded separately; the earlier three-profile result applies to
`e95142bf...`, before this formatting-only follow-up. Full upstream parity,
final-snapshot integration/profile/platform/RSS gates and the Source-shared
external two-process-worker scheduling boundary remain open. No goal-completion
claim is made.

Final formatter-only package
`c8fe6a374e452bb9876d3575ce529a6f2eaecf86df4684ccfe33d3cfc5061622`
passes default all-target compilation and 51 focused tests (22 multinode,
29 formatter). Workspace production/test inputs still match that frozen package
exactly. The previous three-profile package remains the evidence for optional
and Windows compilation before this small common-formatter substitution.


## Natural parallel controls and reopened ownership gaps (2026-11-08)

Immutable release `d4fac90` (package `c8fe6a37...`, binary
`c1f400dc57bf994231f47b94b183374d8cb0c1dcddea9e1ab16c5d29e48f3db1`)
passes three natural 9,704-record / 32-volume / two-round OID-output controls:
legacy one-thread reuse, one-thread task-size-one and four-thread task-size-one.
All cases use one cross-process worker.
All mapping bytes, raw stdout/stderr, raw alignment row multisets, counters and
coverage match their pinned Source captures. Exactly two new Source cases were
captured once; the earlier legacy Source capture was reused. Every process is
terminal and reaped, and every generated artifact difference is retained.
Evidence: `.tmp/native-post-translation-20261107/external-parallel-controls/d4fac90-proof.json`.

These are defined-output passes, not complete artifact parity. Full comparisons
found retained chunk cache files absent in Source, a serial representative queue
(Source p4 terminal count36 versus native33), and missing repeated separators
in serialized seed paths. Those concrete code gaps are being translated in the
next batch. Elapsed seconds Source/native were 44.04/46.28, 39.68/46.80 and
36.21/41.14. Peak RSS KiB were 133016/190200, 172740/191016 and 322100/306056.
These are measured points, not a performance parity claim. Thread-name/file-order
and elapsed-time differences remain separately visible, not normalized away.

The next code batch also replaces serial joins with bounded parallel fetch,
merge/cull and rendering in ordinary, DNA, multiprocessing and supplied-owner
clustering. Only ordered commits serialize. Readonly prepared DAA dictionaries
avoid full per-worker copies; explicit fetched state preserves Source's empty
DAA record after all targets are culled. CLI join-worker errors use Source's
raw stderr/log catch and CRT exit1; embedded callers preserve typed errors.
External accession construction, non-POD radix/loading and final volume output
now use configured workers and retained Source text/stream owners. Prior-round
cleanup always parses the required manifest, including Source error grammar.
These changes await their combined frozen-package verification.

The non-POD input allocator has a narrow Source race boundary: input_buffer.h
workers allocate both pmr list nodes and strings through one shared monotonic
resource before the splice mutex. GNU11 memory_resource::do_allocate mutates
its current buffer/available fields without synchronization. Concurrent
allocation schedules can race; this does not invalidate every multithread run
or excuse unrelated record/order differences. Rust uses independent owned
records and retains the actual decoded storage through the enclosing phase.

Output callback traits in align/output.rs are library reference surfaces, not
production backend implementations. Source generate_output has concrete
counterparts in NativeSearchWriter/NativeQueryRenderer and the real format
printers, preserving physical OIDs/raw title, sequence, self-score, taxonomy,
mate and quality owners. Direct final statistics are counted once; intermediate
emission does not count them, and the final join owns its counters. Legacy
capture feeds these real owners rather than a target-ID side table. The shared
generate_intermediate_output is invoked by the real block writer. No missing
workflow was identified merely from test-only callback trait implementations;
this is an explicit alternate Source-function mapping, not a trait-integration
claim. Full translation/parity completion remains unproven.


Combined parallel/ownership checkpoint (2026-11-08): frozen package
`c2d5604704f484c6efb3baa2968ca04ab8752d19679212cb02f82bd67f3095ce`
passes isolated default, dna/mcl/zstd and Windows GNU all-target compilation.
All 101 selected tests pass: external77, join11, fetched-renderer1,
worker-exit2, native blocked-search3 and supplied-block library7. The worker-exit
parent launches the guarded probe subprocess to verify raw stderr/log bytes,
status 1 and no outer object unwind; ordinary embedded calls return normally.
Ownership tests verify final-clone payload release before cache removal,
constructor-failure file retention, configured representative terminal claims,
literal path bytes, shared radix counts and Source cleanup error grammar.
Native blocked outputs and supplied protein/translated contexts pass their
existing exact format controls. Package and workspace source/test bytes match.

The first package `21a3367f...` passed default compilation but stopped at an
optional DNA module qualification error. The corrected absolute module path is
the only source difference before the passing batch; the failure capture is
retained. No new Source case was run during these checks. A new release and
Native-only comparisons against the preserved Source controls remain necessary
to prove the three natural artifact fixes. Final integration/profile/platform
parity and the broader goal remain open.

Additional code follow-up: FileArray still uses recursive create_dir_all where
Source uses its single-directory mkdir contract (including ignored EEXIST).
This can change the early error class for a blocked radix directory; normal
artifact checks do not establish constructor-failure parity. Caller directory
creation and Source's Windows ignored-result behavior need a coherent callee/
caller translation before that error-domain gate can close.


## Directory/log callee and caller checkpoint (2026-11-08)

Frozen package `78ff05d0fdf467968094368701f8a3cfeff7c2517034bb78bfa4e0e63c666400`
passes isolated default, dna/mcl/zstd and Windows GNU all-target checks, 78
external tests and the focused single-directory test (79 selected tests total).
Source-style single-level mkdir uses mode 0755, ignores POSIX EEXIST (including
regular files), and preserves Windows's ignored CreateDirectory result. Actual
callers create directories in Source order; convenience-only recursive APIs
remain explicit. Job retains its real FileStack log owner, opened after root/
round0 and before worker-ID claim. Optional unchecked diagnostic stream failures
stay silent; actual FileStack errors propagate. Windows fixture cleanup drops
live Jobs after assertions, preserving Source handle ownership/share modes.

The prior package `620f030c...` passed all three compile profiles but found eight
runtime failures: one missing Source clustering-directory creation affected
three flows; five direct fixtures had relied on implicit recursive parents.
The corrective batch adds that real caller step and explicit fixture parents,
with all output assertions retained. Failure evidence stays separate.

The immutable c3ea396 release passes all 24 defined checks across the three
natural controls, using preserved Source captures only. All three identified
artifact defects are fixed: consumed chunks disappear, representative terminal
claims match33/36, and forced-control seed manifests match byte-for-byte with
literal repeated separators. Remaining classified differences are 17 directory
modes (addressed by this new mkdir code but not yet release-verified), thread
files/row order, elapsed log time, and legacy work-root paths. Evidence is
`external-parallel-controls/captures/c3ea396-first/recovery.json`.

Three new short Source constructor controls were captured once with binary `17a...`
and all input hashes unchanged: missing parent, existing versioned-root file,
and existing seed radix-directory file. All stop at the intended Source stage
with status 1. They preserve raw streams/argv/artifacts under
`external-mkdir-controls/source-proof.json`; Native comparison awaits this
checkpoint's immutable release. Source was not rebuilt or retried.

The next concrete code batch is external clustering: Edge InputBuffer loading/
sorting, per-member assignment workers, assignment-volume reads to EOF rather
than manifest count, and parallel closure compression. Source keeps the
bidirectional degree/assignment loops serial; those must remain serial. Exact
plan/hashes are retained under `external-cluster-parallel-audit/plan.json`.
A Source closure race can occur when one part compresses a link another reads;
Rust must use a real safe shared mapping owner, not unsafe aliases or a stale
snapshot. This qualifies those race schedules, not all multithreaded inputs.
Full translation and final parity remain open.


## External clustering workers and Source POD EOF checkpoint (2026-11-08)

Frozen package `da7c86366c2839e1546810aefa1ff6c655ece82c6fe0f00b799feb8d3f2df74f`
passes isolated default, dna/mcl/zstd and Windows GNU all-target checks plus 83
external tests. Actual Edge loading and sorting, per-member assignment workers,
assignment-volume reads and closure compression use configured workers. The
canonical mapping is one shared AtomicU64 vector initialized from identity;
valid forest roots are compressed in place without unsafe mutable aliases or
an obsolete snapshot. Source's cross-partition non-atomic compression race stays
explicitly qualified; bidirectional degree/assignment loops remain serial.
Assignment input consumes complete 16-byte records to EOF, ignoring declared
counts and a final partial record just as Source does.

The shared POD input owner now preserves Source's zero-filled unconsumed records
and actual partial-record bytes when bulk input ends early, propagating real
read/close failures. Pair 24 / Chunk 12 / Edge 24 / Seed 20 layouts decode native-endian
bytes through bounded per-worker scratch; there is no second whole packed-input
copy. Source hash-key partition boundaries and per-worker BufferArray lifetimes
are retained. New controls exercise partial inputs, cross-partition chains,
incorrect assignment manifest counts and member groups crossing nominal parts.

A subsequent one-file realign correction uses startup's resolved thread count;
omitted/-p0 no longer adopts RAYON_NUM_THREADS over Source hardware concurrency.
Package `0132764a620e5c424e034bd63ffcb37d7d5f670971a07907bcb4f192808fc437`
passes default all-target compilation and three existing realign integration
controls. The earlier three-profile/83-test result precedes only this small
thread-budget substitution. Current production/test inputs match that final
frozen package exactly. Runtime/platform/upstream gates remain open.

The immutable `8965365` release passes all 24 defined checks across the three
natural controls; all 17 directory-mode differences are closed in each case.
Remaining raw differences are classified thread file/row ordering, elapsed
log time and old legacy work-root paths. All three constructor controls also
match status 1, raw stdout/stderr and every artifact exactly against Source.
No Source reruns occurred. Evidence: `external-parallel-controls/captures/8965365-first/recovery.json`
and `external-mkdir-controls/native-8965365/comparison.json`. Those release
proofs precede this new clustering/POD/realign code; fresh Native comparisons
and broader final-snapshot command/profile coverage still remain necessary.


## CI run 42514413206 repair (2026-11-08)

The pushed `038f53c` run exposed one default-test failure on each macOS host and
eight on Windows. Windows text-position queries could return -1 for CRT LF-only
cookies. Native incorrectly set the input failbit there, unlike Source tellg,
which discarded later FASTA records. The input owner now returns the failed
position while preserving line extraction; a literal-LF regression retains all
five original lines. This preserves Source cookie quirks rather than inventing
new offsets. Random-access/index fixtures use platform text input and assert
exact physical positions, output bytes and the original full record counts.

Other fixtures now build paths component by component for native separators.
The FileArray regular-file fixture distinguishes Darwin's earlier ENOTDIR
mkdir failure from Linux's ignored EEXIST and Windows's ignored directory result;
all branches retain exact error-stage/path assertions. Production directory
semantics and archive/manifest bytes are unchanged.

Frozen package `221e26e3a15c539e269e0a0c20b85f86d03a620509d46d593f4ea3d7ea95f0ca`
passes isolated default, dna/mcl/zstd and Windows GNU all-target checks plus 125
focused tests (external 83, multinode 22, FASTA 15 and text 5). Current source/test
bytes match the package. These local checks do not prove hosted macOS/Windows
runtime success; the fixes require a pushed rerun. The original Linux job's
optional-feature step is still running and is monitored without restarting it.
Translation/parity completion remains unproven; CI repair is the current priority.
