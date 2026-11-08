# End-to-end translation completion

Goal established 2026-11-06; not achieved. The user resumed translation work
on 2026-11-08 while the next hosted CI run is monitored.

Completion requires native Rust implementations integrated into the executable
and library, verified parity against pinned upstream DIAMOND, and no production
command or option depending on C++ FFI. C++ may remain only as an explicit
testing/conformance adapter behind a default-off Cargo feature. File audit
coverage and green tests alone do not establish application completeness.

Use the source-to-source method in `translation/README.md`: translate callees
before callers, preserve source structure and identifiable function mappings,
audit every function, test it, integrate its concrete dependencies, and compare
its observable behavior with upstream. An abstract backend implemented only
by test doubles leaves its production integration incomplete.

## Current integration checkpoint (2026-11-08)

The raw OS argument/path, byte-option, GVC and TSV integration batch is now
committed. Isolated default and optional-feature checks, Windows cross-checks,
focused I/O tests and all optional integration targets passed in their recorded
snapshots, including all 20 frozen regression hashes. All 18 raw-OS CLI controls
match the preserved Source captures exactly; the six GVC and thirteen TSV
controls also match. These are bounded comparisons, not universal parity.

CI follow-ups `cdbd43e` and `437fc9e` repair workflow option projection, independent
zstd-frame reset, platform fixture assumptions and Windows test handle
inheritance. The workflow now has seven jobs covering four native hosts, musl,
and two 32-bit targets. Both hosted macOS spill checks passed in run 42472717649;
verification of the latest CI repair remains pending.

Resumed code work addresses raw partition-file parsing, retained multinode
library byte/path semantics, and external clustering's volume read boundary
and parallel alignment ownership. The frozen corrective package `5c816812…`
passes isolated default/optional/Windows all-target checks and 86 focused tests.
The activated task pool also now preserves Source's completion-notification
mutex; the FileArray owner avoids allocating worker-local codec state.
Multinode remains a retained library translation;
Source's inactive CLI call is not made into a new command. The shared external
two-worker empty-chunk scheduling hazard remains an explicit open boundary.
The subsequent frozen package `b46aac0c…` adds actual same-pool nested DP,
parallel volume building, retained multinode input/output ownership, and shared
formatted text parsing. Null-pool dispatch and explicit query-pool dispatch are
separate; single-worker pooled calls retain Source task counts and completed
ranges merge under its mutex. The source-compatible byte manifest parser now
serves the connected external stages as well as the retained library.

All three native compile profiles and 204 focused tests pass (one SIMD case
ignored). Windows text-file fixture expectations were corrected against the
CRT output in hosted run 42478857080; its full macOS ARM64 job passed. Fresh
Source parallel-branch comparisons remain pending. The retained multinode
prepare_inputs indexed partition/load_seqs adapter is still incomplete, and no
full parallel, platform or translation-completion claim is made.

## Current coverage (2026-11-08)

All baseline commands and the enabled EXTRA, DNA and MCL adapters have concrete
native routes. Search controls, output formats, streamed BLAST references,
external clustering and provided-owner orchestration are connected. Production
execution never falls back to C++; only explicit `--legacy` uses the optional,
default-off conformance adapter. This is implementation coverage, not proof of
all command/option combinations or optional-profile parity.

The arithmetic-corrected frozen package (`d7144c14…`) passed 1,876 tests
with no failures, plus all 20 independent release regression hashes. Its full
billion-residue database also matches Source byte-for-byte. Earlier EXTRA
roundtrips, crossed reads and strict I/O-error comparisons retain their recorded
snapshot scope. Subsequent sorter ownership, statistics and cleanup corrections
pass their focused CLI checks and the production taxonomy-spill comparison.
The subsequent full 186a default suite passed 1,880 tests; the optional profile
recovered 1,879 passes after five Source-backed profile fixture corrections.
Later view/serializer changes have targeted evidence recorded below. Hosted
macOS/Windows checks and remaining historical-data/profile gates stay open.
Source undefined behavior and optional Source build failures retain their
explicit qualifications in the dated checkpoints.

The subsequent external-storage and BLAST metadata batch passes both native
compile profiles and 66 focused checks. Direct Source probes then exposed
intermediate compression-buffer differences; the corrective package now matches
all 68 Source snapshots exactly and passes its 12 focused codec checks. The actual one-worker OID workflow now matches the original Source, and the
default accession workflow matches all 170 files and diagnostics against the
explicit GNU 12 output-component Source profile. Original GNU 11 aborts in that
component's allocator-backed list sort. The two-worker run still timed out at
the Source-shared empty-chunk race. The POSIX FileSink stdio/fatal-cleanup
translation passes its frozen-package verification and direct Source controls.
Input/shared-stream and retained FASTA class controls also pass; the latest
spool lifetime, console and filename changes pass native compile profiles.
Platform runtime and further compatibility gates remain open. These results do not establish complete external
or platform parity.

The sections below are historical checkpoints. Statements that an adapter was
missing, automatically routed to C++, or awaiting a later fix describe that
checkpoint; later entries and this current coverage supersede those claims.

The latest constructor/ownership audit found further behavior gaps despite
that route coverage: `makeidx` omitted Header2 initialization, `dbinfo` wrongly
required Header2 in upstream's Header1-only inspection mode, and MCL eagerly
cached titles and propagated output-stream failures that upstream suppresses.
These paths now follow the concrete upstream owners and operation order. The
combined native `dna,mcl,zstd` all-target check passes; final packaged runtime
verification of this batch remains open. Evidence and before/after scope are
under `.tmp/native-post-translation-20261107/`, including
`retained-owner-final-batch-check.json`, `mcl-owner-output/provenance.json`, and
`makeidx-seg-loader/header-order/constructor-protocol-before.json`.

Ten focused constructor cases now match upstream exactly, including status and
both output streams (`final-dmnd-constructor/proof.json`). These include
upstream's old/new database-version diagnostics and precedence.
MCL now uses the retained owner's physical sequential cursor and preserves raw
title bytes through accession scanning and output. An unchanged Source object
probe verifies the cursor, delimiter, EOF and raw-byte behavior; full Source
MCL CLI parity remains qualified by its obsolete, unbuildable API calls. The
latest version-diagnostic and raw-title changes pass the refreshed native
all-target check. The actual MCL CLI now preserves the probe's raw accession
bytes exactly (`mcl-owner-output/native-raw-cli-after.json`); the retained
before capture demonstrates the previous lossy conversion. Refreshed isolated
package checks remain scoped by their recorded source hashes. The 19 default
integration targets affected by disk exhaustion recovered with 55 passing
tests; that run's binary predates these latest constructor/MCL changes, so it
does not replace their focused current-code evidence.

The additional `makeidx` SEG-bit clearing mirrors upstream's loader state.
Existing seed consumers already ignored those bits, and the annotated fixture
produced identical index bytes before this change. It is not evidence of a
previous index-output defect; see `makeidx-seg-loader/consumer-proof.json`.

A later input-owner audit confirmed a compressed-DMND integration gap:
upstream accepts gzip databases for Header1 inspection and sequential getseq,
and rejects random indexing at the actual compressed-stream seek. Native raw
file adapters previously rejected those inputs earlier. The concrete streaming
owner now reaches dbinfo/getseq/makeidx and reference search, including Source
metadata validation/loading order. Eleven compressed-input cases match Source
exactly, including selected records, small buffers, actual unsupported seeks,
matrix/metadata error precedence and output creation. Current blocked-search
and taxonomy workflows also pass. Evidence is in
`compressed-dmnd/current/proof.json` and
`current-package-final/compressed-results.json`.

The actual scale test uses 2,862,293 distinct public protein records containing
1,100,000,118 residues, independently checked against the complete, checksummed
UniParc 2026_03 gzip object. Both default makedb commands succeeded, but six
soft-mask bits differed and therefore changed the header hash. Source's default
AVX2 object uses separate multiply/add; Rust incorrectly enabled FMA from CPU
capabilities irrespective of the compiler profile. The corrected kernel follows
the Source build contract. All six original records match both corresponding
Source arithmetic profiles exactly, and 15 focused tests pass. The corrected
release completed the actual full native rerun with default masking: its entire
1,231,665,708-byte database is byte-identical to the preserved Source database,
and all four crossed dbinfo calls match. This closes this particular large
makedb comparison, not the separate taxonomy-spill, optional-profile or hosted
platform gates. Evidence is in `scale-data/completion-proof.json` and
`scale-data/makedb-current/fixed-proof.json`.

The corrected release also matches all eleven natural AMR matrix cases against
preserved Source captures: all eight sensitivity presets, the four-thread
control, accepted BLOSUM50 scoring and invalid-gap diagnostics. The previous
BLOSUM50 run differed in 25 alignments; the arithmetic correction closes that
failure over the complete 3,000-query input. Evidence is in
`natural-current-release/tantan-fixed/results.json`; this is the arithmetic
snapshot, before the subsequent sorter/statistics changes.

The original 5,000-pair SRR1703350 inputs were recovered from the sibling SKESA
fixture named in the historical command logs. Both exact recorded input hashes
match. Current default Source/native blastx runs at one and four threads each
produce the same 3,081 rows and the historical output hash, including mate
sequence output. Evidence is in `paired-current-release/manifest.json` and
`paired-current-release/proof.json`; this does not reconstruct the separate
missing 500-query mouse inputs.

Preparing the real taxonomy-spill test uncovered a further transport gap:
native sorter runs used named-file reopen operations, while Source transfers
the actual temporary descriptor to its reader. The correction now retains the
real descriptor, unlink policy, buffer size and typed errors, including empty
final runs. Actual taxonomy counts now reach the final statistics table after
output closure. Both native all-target checks and 13 focused tests pass; the
production-threshold native comparison now passes against the completed
Source capture and verified mapping projection, as detailed below. Taxonomy
integer parsing also preserves Source's signed saturation, narrowing cast and malformed-row error precedence. Fresh Source
CLI captures and check results are in `sorter-owner-cli/` and
`sorter-owner-final-batch-results.json`.

The first rebuilt CLI passed five focused cases but exposed an incorrect error
cleanup boundary: Source deletes output only for sequence-loop failures and
retains partial databases for later taxonomy/finalization failures. Rust's
handler had covered both phases. The corrected boundary passes the combined
all-target compile check and 11 focused builder/sorter tests. The library's established `InvalidData`
classification is preserved for validation errors while captured typed I/O
errors keep their original identity. Release `a3275472` passes all six focused
CLI checks. Source also leaves Header1's initial
`pos_array_offset` uninitialized; on failed finalization only, Rust keeps its
defined zero value. Comparisons require exact file presence, length and every
defined byte, while recording this specific eight-byte Source qualification.
Successful databases retain whole-file byte equality.
The complete official 2025_03 mapping archive now matches its pinned MD5
(12,624,381,784 bytes; SHA-256
`87d427395d1e0dd38b248ed8dda557b7d19695b5c6c77717e1fbfbefb5e81932`).
The projection is independently verified at 122,563,424 accepted rows and
2,147,483,664 bytes, including full gzip CRC/EOF. The actual paired spill gate
passes: both unchanged executables report all four identical taxonomy counts
and produce byte-identical 295,634,030-byte databases (SHA-256
`aa48be8e92be68bb8cb6068a03d879069ce38c2413006363d0494646ffe88780`).
Two full-bucket unlinked runs per process have identical sizes and sampled
prefix bytes; temporary directories are empty afterward. Source/native wall
is 93.13/283.91 s and RSS 2,709,548/3,806,048 KiB in this single shared-host
scout. This closes accession-to-taxid production spilling, not matched-OID
spilling, whole temporary-file scans or performance parity. Full raw identities
and qualifications are in `scale-data/taxonomy-spill-current/`.

The final frozen package (`186a046d…`, release `a3275472…`) passes both native
compile profiles, the isolated all-feature/all-target compile without C++,
11 builder/sorter tests and all six focused CLI comparisons. Its complete
default suite passed 1,880 tests with no failures. The `dna,mcl,zstd` profile
completed all 33 targets; five baseline-only fixture expectations were corrected
against exact EXTRA captures, and focused recovery gives 1,879 passing checks,
33 ignored and one release-covered skip. Original failure logs remain intact.
These broad results precede the later view/serializer changes.

The later EXTRA view gap is now connected: compression, SAM/XML controls,
real Search::Config/output construction and taxonomy-name loading use their
actual runtime owners. Public ViewConfig/run APIs remain compatible, and
baseline workflow restrictions stay intact. Source comparisons also corrected
the DAA TSV mate callback (`*` for absent mate) and database-serializer error
ownership: real tell/seek/full-chunk commits survive, while pending bytes are
discarded on finalization failure. Typed names/nodes errors retain their OS
reason, and the actual post-mapping tell commits metadata before nodes loading.

`view-extra-compression/consolidated-proof.json` records 35 distinct comparisons
across two explicit snapshots: 27 retained passing codec/constructor/XML/SAM
cases and eight fresh mate/buffer/taxonomy-error cases on package `87f9421f…`.
Successful gzip bytes match wholly. Failed DAA gzip archives additionally
qualify Source's uninitialized decoded bytes [24,144); failed DMND files qualify
only initial Header1 bytes [32,40). Native zeros those fields and every defined
byte is checked. Both current compile profiles and 20 focused view/builder/nodes
tests pass. The current isolated package also passes default public-dispatch/
stdout tests (4/4), the all-feature/all-target compile with C++ absent, and the
Windows GNU native-feature cross-check. This is cross-compilation, not hosted
Windows execution. These selected results do not establish every remaining
workflow interaction.

The next external-job audit found a queue-layout mismatch: Source places job
files beneath `parallel_tmpdir/diamond-tmp-2.1.24`, while the native CLI adapter
passed the parent directory directly. The adapter now adds Source's versioned
suffix; the library Job API keeps its explicit-root contract. The actual
multi-worker integration comparison remains pending.

## Required coverage

- [x] Inventory the full pinned upstream command and option registry, including
  optional build configurations. Optional upstream features must be accounted
  for explicitly rather than silently dropped from the completion claim.
  The source inventory audit now includes `.cc`: the historical 198-file map
  omitted `lib/tantan/LambdaCalculator.cc`. Its completed numerical/function
  audit now brings the file inventory to 199/199; command/option completeness
  remains a separate unfinished requirement.
- [x] Connect all standard commands to concrete native implementations.
  Current routing covers all 21 baseline registrations. Actual production
  adapters exist for all 14 active backend contracts; two retained metadata/
  record interfaces are explicit compatibility models and Source's multinode
  invocation is inactive. The ordinary routes never call C++. This checks
  connection only: the separate behavior-verification requirement below
  remains open. Historical command
  checkpoints are retained below, with later entries superseding old gaps.
- [ ] Verify the connected native search controls end to end: database/query
  blocking, algorithm selection, hard/soft masking, composition modes, global
  ranking, iteration, score thresholds and tuning controls. Protein/translated
  default budgets, supplied blocks, streamed BLAST and clustering owners now
  reach concrete implementations; applicable scale and interaction proof remain.
- [ ] Verify search output formats and fields, headers, compression, and
  auxiliary sequence output. Tabular, flat JSON, DAA, pairwise, XML, PAF, SAM,
  edge/null and taxonomy output have native consumers. Registration or a
  successful bounded case alone does not establish every output combination.
- [ ] Verify the integrated translated-search frameshift/SWIPE and query-length
  controls across relevant input, output, blocking and thread modes.
- [ ] Retain and verify NCBI database and split-input clustering compatibility
  boundaries. Native BLAST references are enabled with real metadata/filter/OID
  handling; Source's perimeter UB is documented rather than reproduced. The
  uninitialized Source clustering counter remains a qualified diagnostic
  compatibility policy, without input-size heuristics or universal byte claims.
- [x] Connect translated library dispatchers to concrete backends and ensure
  the public command-running API does not require C++. Native RunBackend owns
  typed preparation/execution; native AlignBackend and owned double-indexed
  adapters own actual query/reference/output lifetimes for protein, blastx,
  DNA, multiprocessing and provided-owner clustering. Default and explicit
  dna,mcl,zstd all-target builds compile. The library preserves the host stdout
  descriptor; only the executable enables Source process-descriptor ownership.
  This checks concrete integration, not complete runtime/platform parity.
- [x] Remove automatic C++ execution: the executable invokes the default-off
  oracle only for explicit `--legacy`; unconnected production routes return
  native errors. This does not complete native coverage or parity. Unknown
  commands/options receive native registration diagnostics.
- [ ] Verify every completed command/option against the pinned upstream build
  with corresponding build features, retaining commands, outputs, statuses,
  and diffs. Use natural data for workflow/performance comparisons and focused
  fixtures for algorithm boundaries. Cover applicable thread and memory modes.
- [ ] Keep native Linux/macOS/Windows CI green and verify isolated packaged
  builds without the adjacent C++ tree. Cross-compilation is not runtime proof.
- [ ] Reconcile README, file audit, TODO, and validation claims with the final
  native coverage and evidence before declaring completion.

## Verified native additions

2026-11-06: `greedy-vertex-cover` now reaches the translated graph implementation
using the translated upstream command-line parser. Eight retained upstream
snapshots cover directional/default and triplet edges, symmetry, 50/80/95%
member coverage, reassignment, and connected-component depths 1/2.
`tests/native_gvc.rs` runs without C++ sources or the FFI feature and also checks
option cardinality and the symmetry constraint.

External comparison evidence is under `.tmp/native-gvc-parity-20261106/`:
96 deterministic graph cases have byte-identical cluster and centroid files.
A natural graph generated by upstream self-search of 500 distinct AMR proteins
contains 20,187 edges; six GVC modes at one/four threads give 12/12 byte-identical
comparisons. The natural FASTA SHA-256 is
`60caa8e27deaee62141bd7ccb79eaa717bc086243178451d09f6873c5672c04f`.
Raw commands and results for that graph are retained in `natural/commands.json`
and `natural/results.json`. This evidence establishes the tested graph modes,
not completion of other commands, performance parity, all malformed-input
diagnostics, or all graph I/O/scale boundaries.

2026-11-06: CI package verification now includes native macOS ARM64 and Intel,
with compiler-host assertions and isolated all-target/all-features builds.
The new hosted jobs still need to run. Existing macOS checkout build/test jobs
were already present; they must not be described as newly introduced.

Local validation: `cargo test --offline --all-targets` passes (1,736 library
tests passed with six ignored; all runnable CLI/integration targets passed).
Actionlint and workflow shell/YAML checks pass. Changed Rust files pass
rustfmt and `git diff --check` passes. Workspace-wide `cargo fmt --all --
--check` reports an existing formatting difference in the Windows
`UnlockFileEx` call in `src/util/parallel/filestack.rs`; that file was not
changed by this work.

2026-11-07: removed the automatic C++ route for `blastp`/`blastx --masking seg`.
The existing native masking path applies SEG to targets and leaves queries
unmasked, as upstream does. Forty natural tabular comparisons cover protein
default/fast/sensitive and translated default/mid-sensitive modes, CBS 0/1,
one/four threads, and Rust 16G/0G memory modes. All are byte-identical.
Four natural DAA comparisons (both commands and CBS modes) are byte-identical
after the wire fixes below. Artifacts, commands, hashes, and the tabular driver
are under `.tmp/native-seg-parity-20261106/`; use `daa-fixed-results.json` for
the post-fix archive results.

The archive comparisons exposed broader native DAA defects, now fixed:

- `Hsp::push_match` writes one count-1 operation per matching cell. Native
  traceback had combined these into runs, preserving alignments but changing
  serialized bytes. Long insertion chunks also need the backward traceback's
  remainder-first order. A focused 127-cell gap regression covers that boundary.
- Default DAA output stores the first subject ID rather than its full title.
- Upstream serializes `Config::reward`/`penalty`, which its parser never
  initializes. These are distinct from `match_reward`/`mismatch_penalty` and
  unused in protein/translated search. Rust now deliberately writes stable
  zeroes, matching the pinned executable's observed header values. This is an
  explicit deterministic compatibility choice at an upstream undefined-value
  boundary, not a claim that every upstream build emits the same unused bytes.

Eight self-contained upstream DAA snapshots and eight tabular snapshots are
retained in `tests/data/seg/`, including masked/unmasked controls that produce
different alignments. The integration test checks tabular output at 16G/0G
and exact archives for both commands with CBS 0/1. The full native all-target
suite passes with 1,737 library tests passed and six ignored; raw output is
retained as `all-targets.log` in the same artifact directory.

The developer-only `scripts/audit_translation_inventory.py` verifies that
recorded modules still exist and reports unmapped C/C++ implementation files.
At that checkpoint it failed for the five removed test modules and omitted
`.cc` audit entry; later checkpoints restore those modules and the entry. It is not a parity proof and adds no upstream-source
dependency to Cargo builds or tests.

2026-11-07: native `blastp` and `blastx` now support `--max-hsps`, including
zero for unlimited HSPs. The limit crosses both command configuration boundaries
and both extension configurations. Tabular output emits every retained HSP;
DAA already serializes every retained HSP. Cluster-internal searches retain
their upstream one-HSP configuration. The CLI enforces the upstream scalar
parameter count and unsigned numeric conversion before selecting a search
route; malformed options no longer fall through to C++.

Evidence is under `.tmp/native-hsps-parity-20261107/`:

- `results.json`: 32/32 byte-identical natural tabular comparisons. Protein
  more-sensitive and translated mid-sensitive modes cover max-HSP 0/2, CBS
  0/1, one/four threads, and Rust 16G/0G memory modes. The protein scout gains
  two genuine additional rows over the single-HSP control (1077 vs 1075 at
  CBS 0; 1073 vs 1071 at CBS 1). The translated natural scout retains 25 rows;
  actual forward/reverse multiple-HSP expansion is covered by the compact
  regression rather than inferred from that unchanged natural count.
- `daa-results.json`: 8/8 natural archives are byte-identical for both commands,
  CBS modes, and max-HSP 0/2.
- `seg-combined-results.json`: 16/16 tabular/archive combinations also pass
  with SEG masking, checking the interaction with the previous native addition.
- `field-mode-results.json`: 24/24 default, ID-only, score-only, and
  statistics-only output comparisons pass on compact and reduced regression
  inputs. The requested fields must select upstream's actual DP requirements.
- `final-provenance.json` records input/binary hashes, compiler version, and
  upstream commit/build features. Existing local C++ changes were inspected:
  source changes are whitespace/comment encoding and one header include;
  CMake adds the optional static adapter without replacing executable logic.
- `all-targets-final.log`: all runnable native tests pass, including 1740
  library tests (six ignored), 19 binary tests, and the integration targets.
  The ignored external-data tests remain outside this validation claim.

The comparisons exposed and fixed three integration defects:

- Narrow SIMD traceback used moving-band row offsets for full-matrix bounds.
  Padding for a long target could wrap a byte counter and produce a negative
  query endpoint. Full-matrix bounds now use query-row counters and promote
  score widths when those counters cannot represent the query. Focused tests
  exercise long-target padding and the counter-width boundary.
- Full-matrix HSPs must retain their constructor's zero diagonal bounds.
  Artificial full-range diagonals changed equal-score culling precedence and
  displaced a valid earlier banded alignment. All native HSP constructors now
  preserve the upstream key, with a focused regression and retained AMR case.
- Tabular DP requirements now come from `FieldId::hsp_values`, including
  transcript requirements for positive counts. Direct-output DP statistics
  remain distinct from the transcript parsing used by upstream's joined-block
  and DAA-view paths; reparsing direct output would change gap-opening counts.

`tests/native_multi_hsp.rs` uses self-contained upstream snapshots from
`tests/data/multi_hsp/`. Reordered domains produce three disjoint HSPs, and
translated input includes forward/reverse-complement records. Limits 0/1/2/3
exercise actual row-count changes, memory/disk processing, and exact archives.
A reduced AMR regression retains the equal-score alignment and weaker disjoint
HSP statistics. Prior failing comparisons and the initial panic are retained
in the artifact directory; `before-output-integration/` is historical evidence,
not the final result set.


2026-11-07: standard protein-search composition modes 2–5 now run natively.
Modes 0/1 remain available for translated searches; the standard upstream
rejections for translated modes 2–5, mode 6+, and global ranking combined
with adjusted scoring are preserved before opening inputs. The public native
blastp/blastx entry points enforce the same mode restrictions. This covers the
standard build surface, not the experimental EXTRA mode or hidden CBS tuning
options, which remain part of the wider option inventory.

Connecting these modes exposed four integration errors: conditional selection
received background scores instead of frequencies; optimizer defaults disabled
iterations instead of using tolerance 1e-8 and 2000 iterations; scaling used
gapped rather than ungapped lambda; and mode 5 reconstructed its starting
matrix from integer scores instead of using the separately rounded matrix
built from frequency ratios. The existing translated matrix constructor and
numerical routines now supply those dependencies.

`tests/native_cbs.rs` retains upstream tabular/DAA snapshots for all six standard
modes, one/unlimited HSPs on reordered domains, and a reduced natural CopS/
histidine-kinase case that catches mode-5 starting-matrix errors. Native tests
exercise both 16G and forced-disk 0G output. The reduced natural case was
necessary: the compact case alone passed despite that error.

Artifacts are under `.tmp/native-cbs-parity-20261107/`. `matrix-results.json`
contains 32/32 exact compact comparisons across eight matrices and modes 2–5.
`natural100-results.json` contains 4/4 exact more-sensitive unmasked comparisons
of 100 natural AMR queries against 1000 references (1078 rows per mode).
`before-scaled-matrix-fix/` preserves the failing mode-5 comparisons;
`provenance.json` identifies the corrected debug/release binaries and inputs.
`all-targets-final.log` records the passing full native suite: 1740 library
tests with six ignored, 19 binary tests, and all runnable integration targets,
including the three new composition tests. Changed Rust formatting and
`git diff --check` pass. No FFI feature was enabled for these Rust runs.

`combined-results.json` additionally records 16/16 exact natural comparisons
for modes 2–5 with tantan/SEG: tabular at four threads with forced-disk Rust
storage and DAA at one thread in memory. The tabular comparisons contain
1064 rows with tantan and 1076 with SEG. DAA equality compares whole archive
bytes; the artifact driver’s newline counts have no row meaning for DAA.
These comparisons establish the listed scenarios rather than universal
composition, masking, matrix, or full-translation parity.


2026-11-07: `realign` now has a concrete native command route. The source-
mirrored workflow uses native FASTA/DMND readers, the translated cluster parser,
full-matrix SWIPE, tabular format validation and rendering, and a native thread
pool. BLAST databases follow the pinned upstream accession-lookup rejection
policy rather than routing through FFI. Database loading preserves OIDs while
using separate position and sequence cursors.

This integration exposed disconnected behavior in the existing translated code:

- Full-matrix output follows asynchronous SIMD lane completion and refill,
  including early score saturation and promotion. Native batched computations
  now replay those completion events; score-only full matrices also preserve
  the byte-to-word promotion instead of bypassing it with a wider calculation.
- Split-block output must merge stream heads with the upstream heap policy.
  The former flattened stable sort and linear merge changed equal-centroid
  member ordering. `merge_sorted_files` now explicitly follows the pinned
  Linux heap operations, with an independent five-stream C++ snapshot test.
  Upstream's comparator is non-strict for equal keys, so matching that chosen
  executable is an explicit compatibility policy; other C++ standard libraries
  are not claimed to produce the same ordering.
- Upstream sends serialized/deserialized HSP contexts to the formatter. That
  record omits query sequence data and corrected bitscores. Native callbacks
  now reflect this boundary, including the resulting normalized-identity
  denominator. Upstream also leaves `hsp_num` uninitialized on this path;
  native Rust deliberately retains zero instead of reproducing the observed
  stack value 32767. This field is checked against the stated policy rather
  than counted as strict upstream byte parity.

`tests/native_realign.rs` retains 24 upstream snapshots for compact reordered
and 43-target lane-refill fixtures, plus two extra-field snapshots. It checks
96 core output combinations across FASTA/DMND, whole/split blocks, modes 0–5
on the compact fixture, and one/four threads. Extra-field tests cover titles,
coverage, signed ordinal offsets, normalized identity and corrected scores.
Invalid format/workflow/mode diagnostics are checked before input I/O.
On AVX2 the snapshots compare whole bytes. Other SIMD widths compare complete
rows because upstream's lane width changes completion order; exact ARM output
ordering still needs comparison against a corresponding upstream executable.

External evidence is in `.tmp/native-realign-parity-20261107/`:
`results.json` has 88/88 exact compact/natural cases; `dynamic-results.json`
has 24/24 exact lane-refill cases; `extras-results.json` has 6/6 exact defined-
field comparisons at four threads. The natural fixture contains 100 AMR
proteins and 100 mappings to ten centroids. Input/binary hashes and the
undefined-field policy are in `provenance.json`. Earlier diffs are retained in
`historical-diffs/` and `before-intermediate-context-fix/`.
The completed initial-debug cases include a 1K split natural run; broader
natural coverage uses 32K, and lane-refill coverage uses 16K, to avoid the
upstream thread-pool heartbeat delay on hundreds of tiny block pairs.

The final full native all-target suite passes with 1741 library tests and six
ignored, 19 binary tests, and all runnable integration targets, including
three realignment tests. See `all-targets-final.log`. Changed Rust files pass
rustfmt and `git diff --check` passes. No Rust run enabled the FFI feature.

Remaining `realign` completion work includes exact no-argument documentation
and global help registration, comprehensive malformed-input and warning policy,
scale/memory behavior of the translated in-memory intermediate contexts versus
upstream temporary files, and platform-specific upstream conformance. The
native route and the bounded comparisons do not establish full command-surface
or whole-translation completion.

`final-native-results.json` rechecks all 118 comparisons with the final release
binary against the retained upstream outputs, recording both SHA-256 digests
and each executed native command. All 118 are exact. The final binary and
relevant source hashes are recorded in `provenance.json`; this recheck includes
the intermediate-context corrections rather than relying on the earlier sweep.


2026-11-07: `reassign` now has a concrete native database/search/output adapter
for its connected command options. Its source-mirrored workflow initializes
thresholds, reads the complete clustering, extracts member/representative
subdatabases in OID order, performs native protein search, consumes typed edge
records, updates membership, and emits the translated clustering output.
The database adapter is shared with `realign`; Rust does not load the adjacent
C++ tree for these paths. The initialized masking, identity and diagonal-filter
settings now cross the previously disconnected `ReassignSearchConfig` boundary.

The integration exposed an additional native search bug: equal coverage cutoffs
length-sort both sides, but in-memory edges reported the sorted indices rather
than original subdatabase OIDs. Edge-mode searches now retain and restore both
permutations. The asymmetric coverage fixture changes representatives between
member-only and mutual coverage, detecting both query and target OID mistakes.
Configured tantan maskers also now supply the same threshold for hard and soft
masking, as the upstream singleton does.

The command has 28 registered native options. Eighteen advanced options still
route through the existing legacy adapter via the explicit
`commands::reassign::requires_legacy` list. `option-inventory.json` records the
source registry and that partition. This is route/parse coverage, not proof of
all command-surface behavior: exact no-argument documentation, verbosity/logging,
remaining tuning dependencies, and exhaustive malformed-input behavior remain
part of the full completion goal.

Upstream defects prevent counting ordinary nonempty command runs as unmodified
parity. The pinned executable fails with `Block::ids()` or a segmentation fault:
its temporary member database lacks the sequence/title flags expected by search.
Its inline `update_clustering` also checks an unsigned NIL sentinel with `>= 0`,
then can index outside the centroid vector. Native Rust supplies concrete query
records and retains the existing checked NIL policy: an unmapped member keeps
its previous assignment. An isolated diagnostic executable repairs just those
two call sites. The pinned source tree/executable were not modified.

Soft SEG exposed a further upstream undefined-input boundary. The diagnostic
crashed in vendor `s_ShiftWin1`; its 128-entry alphabet arrays are indexed by
signed sequence bytes. Instrumentation observed high seed-mask bits in those
buffers. A separate diagnostic normalizes residues with `LETTER_MASK` before
SEG, matching the existing native masking path. The original failure, a
signal-handler backtrace, observed tagged-byte counts, and the additional patch
are retained. The SEG-normalized executable is kept separate from the ordinary
two-repair reassignment diagnostic.

Evidence is under `.tmp/native-reassign-parity-20261107/`:

- `diagnostic-results.json`: 40/40 exact output comparisons against the
  explicitly repaired two-call-site diagnostic. Compact and natural inputs
  cover FASTA/DMND, member/mutual coverage, 90% approximate identity, CBS 0/1/4,
  SEG masking, PAM250, full extension, and numeric OID output at four threads.
- `tuning-results.json`: 31/32 equal outputs with further connected tuning
  controls. The remaining entry records the diagnostic's soft-SEG crash and a
  successful native run; it is not counted as parity.
- `seg-normalized-results.json`: 4/4 exact soft-SEG outputs against the separate
  normalization diagnostic, with FASTA/DMND and one/four threads.
- `unmodified-singleton-results.json`: 16/16 exact output/status comparisons
  with the unmodified pinned executable. These include both database formats,
  header/no-header, one/four threads and accession/numeric output. Singleton
  inputs contain no members to search and avoid the broken integration.
- `unmodified-diagnostic-results.json` and `unmodified-enum-errors.json`: eight
  configuration/enum error comparisons match the unmodified executable.
- `final-native-results.json`: the final release binary rechecks all 91
  qualified output comparisons against the retained snapshots. It labels
  unmodified and repaired oracles separately and records hashes and commands.
  `provenance.json`, build scripts/commands and patch files identify every
  diagnostic repair, binary, relevant source file and input.

`tests/native_reassign.rs` retains these explicit oracle boundaries: 56 native
output combinations use documented diagnostic snapshots, plus unmodified
singleton output and configuration-error checks. Both input formats, one/four
threads and 16G/forced-disk 0G paths are exercised. The natural data consists
of 36 distinct SHV/TEM/CTX-M proteins, six representatives and 30 deliberately
misassigned members. PAM250's no-hit result verifies retained assignments.

The final `cargo test --offline --all-targets --quiet` run passes with 1741
library tests and six ignored, 19 binary tests, and all runnable integration
targets, including three reassignment tests. See `all-targets-final.log`.
Changed Rust files pass rustfmt and `git diff --check` passes. All Rust runs
used the default feature set without FFI. This establishes the listed native
paths and compatibility decisions; the full translation goal remains active.

2026-11-07: `recluster` now has a concrete native CLI/database/realignment/
search/cascaded/output backend. The translated recursive workflow realigns the
current assignments, searches unaligned members against existing centroids,
clusters still-unmapped sequences with the translated cascade, recursively
refines that subset, and merges its local OIDs into the original mapping.
Initialized identity/masking/diagonal thresholds cross the search adapter;
mutable coverage, effective database size and cascade scoring settings survive
recursive calls. Intermediate subdatabase and mapback ownership is released
before recursion where the source releases it.

Natural data exposed an integration defect absent from the compact examples:
upstream assigns an empty-but-present `config.iterate`, enabling its automatic
sensitivity schedule. A direct final-sensitivity search produced different
assignments. The concrete shared search adapter now uses translated
`run::config::Config` schedules and removes queries with earlier output from
later rounds. Reassignment retains its absent-iterate behavior. Native search
also now length-sorts linearized targets independently of query coverage,
retaining original OIDs. Shared translated sensitivity enum parsing supplies
the source's sorted permitted-value diagnostics, and the absent shapes-6x10
iteration schedule returns the pinned build's `map::at` error rather than a
Rust panic.

Evidence is under `.tmp/native-recluster-parity-20261107/`:

- `sweep-results.json` and `final-native-results.json`: 84/84 exact outputs,
  replayed with the final release binary. Sixteen compare with the unmodified
  executable for already-valid clusters; 68 compare with an explicitly
  repaired diagnostic. FASTA/DMND, one/four threads, 16G/32K memory budgets,
  identity and mutual coverage, CBS 0/1/4, SEG, full extension, PAM250,
  custom cascade steps, round cutoffs, reassignment controls and OID output
  are covered. The natural 36-protein fixture catches iterated-search changes;
  the compact fixture forces cascading, recursion and original-OID merging.
- `final-error-results.json`: seven status/error comparisons match the
  unmodified executable, including invalid workflow options, incompatible
  coverage cutoffs and invalid/missing sensitivity schedules.
- `reassign-regression-results.json`: the final release binary retains all
  91 previously qualified reassignment comparisons after sharing the adapter.
- `initial-observations.json` retains the unmodified command's recursive-input
  failure. Like reassignment, its temporary query file lacks sequence/ID load
  flags, and its unsigned unmapped sentinel is checked with `>= 0`.
  `oracle-repairs.patch` records exactly those two repairs in an isolated
  recluster object. The original source tree and executable remain untouched.
  Diagnostic outputs are never counted as unmodified-upstream parity.
- `provenance.json`, build/comparison/replay drivers, commands and stderr files
  record the pinned build, inputs, diagnostic patch and source/binary hashes.

`tests/native_recluster.rs` retains three native integration tests: 45 output
combinations against qualified diagnostic snapshots, three already-valid
unmodified output checks and five configuration errors. Its data and explicit
oracle provenance are included in `tests/data/recluster/`. The expensive
PAM250 debug optimizer is exercised once in the retained test; the release
conformance sweep covers all four database/thread/memory combinations.

The adapter registers 33 options; 14 advanced options still route to the
explicit legacy adapter. `option-inventory.json` records this parse/route
partition, not complete semantics. Remaining command work includes propagating
all parsed tuning controls through cascaded rounds (including extension chunk
size and temporary-directory selection), exhaustive validation/help/logging,
file-backed subdatabase ownership instead of decoded in-memory subsets,
source block/index scheduling, bounded memory at scale and platform-specific
upstream conformance. The native route and these bounded results do not
establish complete command or project parity. The translation goal is active.

The final `cargo test --offline --all-targets --quiet` run passes: 1741 library
checks (six ignored), 19 binary checks and every runnable integration target,
including the three new reclustering tests and existing reassignment tests.
See `all-targets-final.log`. Changed Rust files pass rustfmt and
`git diff --check` passes. The final 84 output, seven error and 91 reassignment
regression comparisons all pass with the release binary. Recorded source,
fixture and executable hashes were revalidated after testing. All Rust runs
used default features without FFI; no completion claim is made for the
remaining dependencies listed above.

2026-11-07: native `blastp` and `blastx` now connect completed alignments
directly to the existing translated pairwise (`0`), XML (`5`/`xml`) and PAF
(`103`/`paf`) renderers. These paths preserve transcripts, target/HSP ordering,
source frames, query OIDs, full titles where required, XML database statistics,
and the header/query/footer lifecycle. They do not serialize through DAA or
invoke C++. Remaining unsupported search options retain their explicit legacy
routes; this is not full search-command completion.

The integration fixed three additional mismatches. Blastx previously discarded
full query titles needed by pairwise/XML. XML header escaping used the original
title length as a limit even without a multi-title delimiter, truncating an
expanded escaped title; the source's byte-offset policy for an actual delimiter
is retained safely. XML's E-value header now uses the source iostream's six
significant digits instead of Rust Display's full decimal form. ROC-ID's
existing default-float formatter now shares that implementation, including
notation selection after rounding across an exponent boundary.

The native hit store now uses translated `SequenceSet::partition_with_contexts`
length-weighted boundaries, reduced to source-query boundaries for blastx and
expanded back to context indices for storage. Disk migration retains the same
partition instead of rebuilding count-based boundaries. Repeated empty final
bins need no physical buffer. A 44-query heterogeneous fixture exposed another
upstream dispatch condition: ordinary `align_worker` skips queries without any
retained seed hits, even inside a populated bin. Such queries do not emit an
unaligned record. Seeded queries whose alignments fail final filters can emit
one. The native per-query guard now mirrors that source distinction. `--unal`
is parsed with source scalar cardinality/conversion for the three new formats;
its unsupported tabular/DAA paths remain explicitly legacy-routed.

Evidence under `.tmp/native-search-output-20261107/`:

- `final-native-results.json`: 101/101 complete output files match the
  unmodified pinned executable. The 96-case matrix covers both commands,
  pairwise/XML/PAF, composition 0/1, one/unlimited HSPs, coverage rejection,
  explicit/default unaligned policy, one/four threads, compact/natural inputs,
  and native 16G/forced-disk 0G storage. Five additional cases cover format
  aliases, database path forms and E-value header formatting. XML paths are
  identical between executables; these are byte comparisons, not normalized
  XML comparisons. Native-only memory settings are omitted from C++ commands.
- `bin-final-results.json`: 12/12 complete outputs match for 44 alternating
  long/short/unseeded protein or DNA queries. These cross heterogeneous bin
  boundaries and retain gaps in original XML query ordinals. The initial
  mismatches and both executables' outputs remain available separately.
- `matrix-results.json`: 12/12 release outputs match for composition modes
  2–5 on the natural CopS/histidine-kinase scoring fixture, for all three new
  formats. With the eight translated tie-order checks below, these give 133 exact
  output comparisons.
- `tie-final-results.json`: 8/8 complete XML/PAF/tabular/DAA outputs match
  for long translated exact matches with zero E-values and equal scores.
  This exposed a pre-existing missing `!query_translated` condition in the
  native minimum-length-ratio optimization. It incorrectly length-sorted
  translated targets and reversed tied hits. Native sorting now follows the
  derived ratio, with the source's protein-only condition. One/four threads
  and native memory/disk paths retain the original target order.
- `reassign-regression-results.json` and
  `recluster-regression-results.json`: all 91 and 84 earlier qualified release
  comparisons remain exact after changing hit partitioning. Their original
  unmodified/explicitly repaired oracle distinctions remain in force.
- `command-inventory.json` separates the 21 standard source registrations from
  the 24 `EXTRA`-only registrations. An existing CLI arm is not backend or
  option completion: the native `test` arm still errors, and `benchmark` still
  needs its concrete kernel adapter. Standard `roc` is a deprecation error in
  `run/main.cpp`; its old ROC algorithm is not the active CLI behavior.
- `deprecated-results.json`: native `roc` and `prepdb` validate the complete
  standard no-option workflow registry, including hidden/disabled names and
  short aliases. Of 644 status/stdout/error comparisons, 642 match raw and
  all 644 match after removing ANSI color controls. The two no-argument
  upstream runs leave a reset on stdout from their old-build warning; native
  console output retains its existing absence of that warning artifact.
  General console/age-warning/color policy is still open, not claimed as raw
  byte parity. Invalid flags no longer require FFI to receive source errors.

`tests/native_search_output.rs` retains compact and natural output snapshots,
heterogeneous-bin records, format aliases, both commands, HSP/coverage/masking
settings, unaligned policy, E-value formatting and native deprecated workflow
validation. `tests/data/search_output/README.md` records their provenance.
Retained XML snapshots replace only the variable database filename with an
explicit placeholder; the external comparisons use unmodified whole files.
`provenance.json` records the final binaries, source/fixture hashes and build
settings. No Rust comparison enables FFI. The full translation goal is active.

Validation for this addition includes a passing full all-target run, followed
by final-state library/binary checks and all five native output integration
tests after the translated tie correction: 1741 library checks pass (six
ignored), 19 binary checks pass, and all five focused integration tests pass.
See `all-targets-verified.log`, `library-cli-verified.log` and
`output-tests-verified.log`. The final release replays retain all 133 output
comparisons and all 175 earlier qualified clustering comparisons. Changed Rust
files pass rustfmt and `git diff --check` passes. These checks do not close the
remaining benchmark, built-in regression runner, library-dispatch, optional
build, search-option, other format, memory/scale or console-policy requirements.


2026-11-07: restored upstream regression modules and native index partitioning
---------------------------------------------------------------------------

The removed `src/test/{data,test_cases,test,queue,hit_buffer_stress}.rs`
translations are restored without `include_str!`/`include_bytes!` access to the
adjacent checkout. The ordered corpus is embedded as Rust literals: 389 records,
97484 amino-acid letters. All twenty source command lines and reference hashes
remain in source order. Ten focused module tests pass. File-map entries now
reflect these implementations and their limits rather than missing paths.
The developer inventory audit finds no missing Rust module; its remaining
failure is the previously omitted `lib/tantan/LambdaCalculator.cc` inventory row.
This inventory result alone does not prove end-to-end completion.

`BlastpConfig`/`BlastxConfig::index_chunks` is passed into the actual native
seed-array partition loop and chunk-aware left-most filtering. Zero selects
the sensitivity default (one for persistent target-indexed search). Native CLI
parsing preserves attached `-c1`, separate short/long values, repeated scalar
replacement, source integer conversion, and parameter-count errors. Existing
internal clustering configurations retain auto selection. This removes the
legacy route for native search index partitioning.

The pinned source registers `query_parallel_limit` but has no read of that
field anywhere in its C/C++ sources. Its case named "target-parallel" therefore
runs the same ordinary search as the more-sensitive case. Native CLI validation
now accepts that scalar and preserves its parser errors without inventing a
scheduling effect. This does not establish target-parallel DP implementation.

`blastp::run_with_records` connects already decoded database/query records to
the existing native seed/extension/output pipeline. `NativeRegressionBackend`
uses this boundary with source-order in-memory SequenceFiles, rewind/statistics
reset, temporary output ownership, and a per-case Rayon pool. It retains the
source runner's empty database-header name. Blocked search, query-indexed
search, and debug logging return explicit unfinished-dependency errors; none is
silently ignored. The full `test` CLI remains disabled. The hit-buffer stress
translation concurrently generates independent per-bin streams then ingests
them synchronously; it does not prove simultaneous packet-writer behavior.

Evidence is under `.tmp/native-regression-20261107/`. The unmodified pinned C++
`diamond test` passes all 20 source hashes. The initial corpus matrix established 18 byte-exact whole-file native/C++ CLI
comparisons; `corpus-results.json` now includes the nineteenth case described
below. The initial subset covered:
default/threaded, three sensitivities, unlimited HSPs, the source's inactive
query-parallel-limit scalar, CBS 0/2/3/4, top/target limits, high E-value,
BLOSUM50, and pairwise/XML/PAF. These commands include the
real native `-c1` integration; no Rust comparison enables FFI.
`tests/native_regression.rs` additionally passes all 18 supported frozen source
hashes through the concrete in-process backend, including the empty XML
database-header name and a separate worker pool per case. The two unfinished
dependencies are tested as explicit errors. Its unoptimized run passes in
324.27 seconds and its release run in 8.62 seconds; no C++ executable is read
by these Cargo tests.

`partition-results.json` records 40 exact comparisons across `-c0/1/2/4/8`,
blastp/blastx, default/very-sensitive, and native memory/disk hit storage. The
source command omits the Rust-only `--memory-limit` extension and otherwise
uses the same arguments. All 58 comparison pairs return zero and produce
identical full files. Library/CLI
validation passes 1751 library tests (6 ignored) and 20 binary tests; scoped
formatting and all-target compilation checks pass. These checks include the
restored corpus and modules. `cargo package --offline --locked --allow-dirty`
packages 650 files and successfully compiles the isolated crate without the
adjacent C++ tree. The isolated package also passes all-target/all-feature
compilation with the C++ tree absent; its diagnostic feature produces the
expected missing-adapter warning. All five existing native transcript-output
regression tests pass after these changes. These are Linux checks; the added
macOS ARM64/Intel CI jobs await a hosted run. `provenance.json` records source,
binary, package, and evidence hashes. The complete translation goal remains
active.


2026-11-07: explicit query-indexed search integrated
-------------------------------------------------

`--algo 1` now selects the real translated hashed query seed set and stage-0
reference filter in native blastp/blastx. `--algo 0` explicitly selects the
existing double-indexed pipeline. Scalar argument validation preserves source
cardinality errors and repeated-value replacement. The native configs carry
this choice into protein and translated query owners; clustering retains its
previous ordinary algorithm. Persistent target indexing continues to use the
same stage-0 implementation with the opposite seed filter.

Query-indexed reference seeding reads unmasked residues. The native ownership
boundary records actual retained-hit target IDs and prepares their configured
hard masks after all shapes finish, before immutable parallel extension. It
reuses the same active-matrix tantan masker used for the query. This hoists the
source's per-ranking-chunk pure masking operation to the retained-candidate
boundary; it can mask candidates that source ranking would never extend.
Match/sequence behavior is verified below; exact lazy-work counters and scale
performance of this scheduling adaptation are not established. No mask is
substituted for the required unmasked seed enumeration. Final query SEED_MASK
annotations move back to the authoritative protein arena before extension.

The real in-process regression backend now passes all 19 supported frozen
source hashes, including the full 389-protein query-indexed case. The remaining
source case uses database/query blocks; regression debug logging also remains
explicitly unfinished. The full CLI runner remains disabled. Broader algorithm
completion still requires automatic query-index selection for qualifying large
databases, contiguous-seed mode, combined query/persistent target indexing,
query-indexed self/linear clustering, and high-partition/scale validation;
these are not covered by the explicit
algorithm-choice evidence here. Combined indexing is rejected explicitly by
the native API and remains an explicit CLI legacy route.

Evidence uses the same pinned unmodified C++ binary. Current corpus comparisons
cover 19 whole-file equalities. The compact and natural-protein query-indexed matrices each have 120 exact
output files and six matching translated-CBS rejection statuses/error messages
(240 full-file equalities and 12 matching errors combined). They span both
query modes, no/tantan/SEG masking, c1/default/c2/c4, unlimited HSPs, CBS variants,
BLOSUM50, all standard search sensitivities including faster/sketching, and
pairwise/XML/PAF. Query/reference histograms receive the same query-indexed
sketch/minimizer settings as enumeration.
Native comparisons use 0G disk hit storage while C++ uses its standard storage.
Four ordinary persistent-index regressions also match for blastp/blastx and
one/four threads, protecting the shared state refactor. The six compact
upstream tabular snapshots in `tests/data/query_indexed/` also drive 18 native
CLI comparisons in `tests/native_regression.rs` across both query types, three
masks, three chunk counts, and memory/disk storage. These Cargo tests need no
C++ checkout; fixture provenance includes the oracle digest and actual source
arguments. Library checks pass
1751 tests (six ignored) and the binary passes 22 tests. Further final-build
and natural-protein comparisons are recorded with the evidence files under
`.tmp/native-regression-20261107/`; the overall goal remains active.


Source registry follow-up: native search now also accepts `--linclust-20`,
`--shapes-6x10`, and `--shapes-30x10`, selecting their existing translated
shape/trait/extension implementations. Previously search routed them to FFI
and `makeidx` silently selected the default sensitivity. Source-compatible
sensitivity selection now rejects mutually exclusive switches and scalar
values supplied to these boolean options; repeated identical switches
remain allowed. The ten registered switches plus default are covered by native
selection tests. The shared selector is used by blastp, blastx, and makeidx.

The expanded two-dataset query-indexed matrix covers all eleven source modes,
with 240 exact files and twelve matched rejections. Another 36 exact files
exercise the three newly exposed sensitivities through double-indexed search.
All eleven makeidx modes produce byte-identical indexes, and all 22 associated
persistent-index searches (both query modes) match the pinned source. Source
index bytes for the three previously ignored flags are frozen alongside the
six tabular fixtures; a native Cargo integration test regenerates those three
indexes and compares every byte. Evidence is in `sensitivity-double*-results.json`
and `sensitivity-index-results.json` under the current evidence directory.


Final package validation for this batch: the 660-file package was extracted
into a separate directory and passed `cargo check --all-targets --all-features`
with no adjacent C++ tree and an isolated target directory. The expected
unavailable diagnostic-adapter warning is recorded in
`query-index-isolated-all-features.log`. The four native regression integration
tests pass in release mode, covering all 19 supported hashes, 18 CLI snapshots,
three exact generated indexes, and explicit blocking rejection. The development
cache was explicitly rebuilt after an earlier package verification left stale
package-relative dependency metadata; diagnostic runs using the older backend
are not counted as validation. Subsequent package checks use isolated targets.
The full translation goal remains active; blocking is the next core dependency.


Persistent-index configuration now validates the already parsed final `-c`
value, instead of re-reading its first raw spelling. This preserves source
behavior for repeated values, leading zeros, and integer prefixes. Fourteen
source/native comparisons pass: eight complete files and six matching
status/error rejections. These cases and every final matrix row record the
actual C++ and Rust binary digests. The final release integration run passes
all four tests; scoped formatting, ABI confinement, 1751 library tests (six
ignored), 22 binary tests, and the refreshed isolated all-feature/all-target
package compilation also pass. `provenance.json` records the final source and
artifact hashes. No final parity runtime enabled FFI.


2026-11-07: native blocked search and complete built-in regression command
------------------------------------------------------------------------

This checkpoint supersedes the earlier 18/19-case and disabled-CLI statements
above. `diamond test` now executes all twenty original cases through the real
native backend and returns success with 20/20 hashes. Its translated command
parser permits the source standard-test option groups, rejects general/search
options unavailable to that command, validates `--algo`, and preserves the
outer configuration reset behavior. `test --bootstrap` matches all twenty
source hashes, returns status 1, and retains twenty output files. Internal
`NativeRegressionBackend::run_search --log` remains unfinished, but the source
standard `test` command itself rejects `--log`; it is not a silent no-op.
The hit-buffer stress scheduling limitation recorded above still applies.

Explicit `blastp`/`blastx -b` now searches actual query and reference blocks.
Protein queries use the existing native FASTA/FASTQ/compression parser with
one-record backpressure; reference DMND records are read incrementally and
rewound for each query block. FASTA references are spooled through native
`makedb` to a private temporary DMND. Translated input streams one or two files
in round-robin order, keeps complete pairs, and budgets the translated ORF
letters returned by `Block::push_back_translated_with_genetic_code`, as source
does. The authoritative six-context query arena and its preprocessing state
survive reference blocks. The same active-matrix Tantan masker is reused.

Completed HSPs pass through the real intermediate codec into temporary files.
The translated block joiner reads one query payload per reference block,
merges source heap heads, performs global culling, reconstructs field-dependent
HSP values and transcripts, and renders tabular/pairwise/XML/PAF/DAA output.
Matched reference metadata is retained without fabricated sequences. When the
reference fits one actual block, source's direct output path is preserved.
Original query/target OIDs survive coverage sorting; ordinary mutual-coverage
search now reaches the already translated stage-1 kernel. Reverse DNA tabular
coordinates use strand orientation, and DAA dictionary IDs are reserved during
reference-block serialization before global culling.

`tests/native_search_blocks.rs` retains fifteen unmodified upstream snapshots
covering both query types, all five output formats, direct archives, reverse
DNA coordinates, sorted query IDs, and mutual-coverage seed filtering. Each
runs at one/four threads and in-memory/forced-disk hit-buffer modes. Snapshot
commands, hashes and XML-path-only normalization are recorded in
`tests/data/search_blocks/provenance.json`. These tests require no C++ checkout.
Current development evidence lives under `.tmp/native-blocks-20261107/`;
earlier debug comparison logs are exploratory and may precede the fixes above.

Limits at that checkpoint (default budgets, taxonomy and supplied-block
adapters are superseded by the following checkpoint): absent `block_size`
still used the monolithic route. Blocked taxonomy output, BLAST database input,
clustering edge output, and supplied translated blocks were not connected.
Those paths stayed legacy or returned explicit unsupported errors. The adapter
currently captures one complete query/reference block's matches before
serialization, so bounded high-scale output queues and performance parity
are not established. Automatic/contiguous/combined indexing, advanced search
options, benchmark, public library dispatch, and other goal requirements above
remain unfinished; production FFI dependency removal is not complete.

Final scoped evidence for this checkpoint: the release binary SHA-256 is
`d3a7a06e10109fd1b1059f057332ae814aa80eea672329aab2bd19ab257895bb`.
All 156 full-file comparisons match the immutable pinned source binary:
60 compact/natural protein outputs across three block budgets and five formats;
50 advanced outputs covering query indexing, matrix adjustment, top/k culling,
and mutual coverage; 30 translated outputs; and 16 paired-read outputs covering
both strands, alternate genetic code, minimum ORF length and full mate sequence.
DAA comparisons hash the actual suffixed `.daa` file, never a missing path.
`provenance.json` records the final implementation and result digests; every
result row records both executable digests. Final native `test` returns 0 with
20/20; final bootstrap returns 1, matches the twenty source hashes, and retains
twenty outputs. Release regression/search-output/block suites pass, including
an explicit unequal paired-file rejection. Library tests pass 1,753 with six
ignored; binary parser tests pass 23. These are Linux runtime checks; the newly
added macOS packaged-crate jobs still await their hosted run.


2026-11-07: default budgets, supplied translated blocks, and blocked taxonomy
---------------------------------------------------------------------------

Normal native FASTA/DMND searches now use the same real block driver even when
`-b` is absent: two billion letters below very-sensitive and 400 million for
very/ultra-sensitive, matching `run/double_indexed.cpp::run`. Explicit zero
selects the same source default. The concrete loader budgets complete protein
records or the actual translated ORF letters and retains complete read pairs.
This connects automatic traversal rather than treating the budget as a no-op.
The pre-existing eager native BLAST library adapter remains only for absent
`block_size`; explicit blocked BLAST input and its separate source
compatibility policy remain open. CLI BLAST search stays legacy.

`Block` retains each translated loader return value, so supplied contexts can
be divided with the exact source ORF budget without inferring it from masked
residues or retranslating with a different genetic code. Source DNA, all six
original contexts, annotations, titles, original OIDs and mate relations pass
through the actual native reference blocks and intermediate merger. Supplied
mate layouts are validated before database I/O. Frozen source tests cover all
five formats plus paired genetic-code-4/min-ORF/mate-sequence output at one/four
threads and both hit-buffer memory modes. No C++ files are read by these tests.

Blocked taxonomy fields no longer route to FFI. The DMND reader validates and
loads exactly the mapping, names, nodes and ranks requested by the fields.
Matched dictionary entries retain their actual taxids and the renderer tree
contains assigned nodes and their ancestors. Original target OIDs survive
length sorting and intermediate serialization. The direct and merged paths
use the same taxonomy rendering helpers as ordinary search. In-memory protein
record APIs lacking database metadata still reject taxonomy explicitly.

Scoped release evidence is under `.tmp/native-block-adapters-20261107/`.
The initial exploratory taxonomy matrix contained 24 comparisons. The final
matrix expands to 64 full-file equalities across both query types, default/
three explicit budgets, mask/no-mask, explicit algorithms 0/1, one/four threads,
and 16G/forced-disk hit buffers. Another 88 default-budget tabular/DAA files
match across all eleven CLI sensitivity modes and both query types/mask modes.
Partial-metadata checks add eight exact outputs and 24 matching source error
statuses/messages. The minimal mapping-only and mapping-plus-names databases
exercise field-specific component loading. Taxonomy database producer bytes
also match. Four unmodified source taxonomy snapshots and the source metadata
inputs are retained under `tests/data/search_blocks/taxonomy/` with provenance.

This closes those specific block adapters, not the full translation goal.
BLAST blocked input, clustering edge blocking, scale/performance proof, remaining
search options and formats, benchmark, public command-running library dispatch,
optional build configurations and automatic production FFI routes remain open.

Final binary for this adapter checkpoint: SHA-256
`04696a6e5d4f5ee0e6b6aec5b796068305e516fd4b59fa9cda0224562bdffa7c`.
The final matrices contain 172 byte-identical output files and 24 matching
error statuses/messages; the additional twelve outputs exercise gzip queries,
gzip references and both compressed together in both query modes, tabular/DAA.
Each row records both executable digests; `provenance.json` records result,
driver and implementation hashes. The corrected default sensitivity cases
omit any sensitivity flag, as the source requires. Earlier exploratory logs
are not final proof.

`cargo test --offline --release --all-targets` passes: 1,754 library tests with
six ignored, 23 binary parser tests, and every runnable integration target.
The debug library suite passes 1,755 with six ignored (one debug-only test).
Scoped rustfmt and `git diff --check` pass. The 691-file packaged crate also
passes isolated all-target/all-feature checking without the C++ tree. These
remain Linux checks; hosted macOS/Windows validation and broader performance
parity are not inferred from them.


2026-11-07: shared native executable and public library entry
-----------------------------------------------------------

`diamond::run` is now an unconditional native library function. It uses the
same parsed configurations and concrete native command functions as the
executable, moved from `src/main.rs` to `src/commands/dispatch.rs`. The parser's
23 existing tests move with that code. This replaces the old feature-gated
public API that always called C++ or immediately returned unavailable.

CLI command failures and parser failures return statuses/IO errors through the
native entry rather than terminating the embedding process. A child-process
regression calls nine malformed/untranslated cases and then successfully calls
version in the same process; it asserts that no FFI-unavailable diagnostic is
used. Another in-process test builds a real database/index, searches protein
and DNA queries at one/four threads with 16G/forced-disk buffers against frozen
source outputs, and writes/views a DAA archive. The library example now compiles
without a C++ feature and is covered by doctests.

Search, database/index, clustering and reassign commands executed by this
adapter receive per-invocation Rayon pools; repeated calls do not inherit the
first call's worker count. Core search and makedb avoid creating a redundant
global pool while executing inside a supplied pool. Commands with their own
worker-pool adapters retain those adapters.

The native library dispatcher contains no C++ call. The temporary automatic
legacy bridge remains exclusively in the executable and is still an open
production dependency. `diamond::run` reports untranslated routes explicitly;
it does not make those features implemented. Benchmark and the remaining
advanced/optional controls are not connected by this entry-point change.
The separate generic source `run::main::RunBackend` workflow still requires
its full concrete-backend/registry integration; this shared CLI adapter is not
proof of that abstract API's completion. Source-level fatal ALP matrix/vector
serialization helpers also remain separate from the connected command routes.

Final validation of this dispatcher checkpoint: release executable SHA-256
`d49508ef42af197df2e35a34686f42794589c72ba3f278e94ba5e5d52f20be58`.
All 152 final-executable comparisons (88 default-budget and 64 taxonomy files)
match the unmodified pinned source binary. Result/implementation hashes are
in `.tmp/native-dispatch-20261107/provenance.json`. The preceding checkpoint's
result ledgers remain byte-identical to their original recorded digests.

`cargo test --offline --release --all-targets` passes with 1,777 library tests
and six ignored; the executable's 23 parser tests now reside in that library
count. Every runnable integration target passes, including all twenty source
regression hashes and both new public-library tests. Two doctests pass. The
693-file isolated package passes all-target/all-feature checking and actually
executes the public-library workflow/error tests with all features enabled,
without the C++ tree. Scoped formatting and diff whitespace checks pass.
These are Linux checks; this entry change does not supply hosted platform or
benchmark/performance parity evidence.

## Standard benchmark backend and anchored padding (2026-11-07)

The standard `benchmark` command and `--benchmark` alias now have a concrete
native command/backend reached by both CLI and library dispatch. They execute
the source iteration counts over real matrix, anchored, SWIPE and diagonal
kernels. Reduced-count tests call the same backend; no test-double backend is
used to establish this integration. Standard benchmark options are rejected
before entering the workloads, matching six pinned executable diagnostics.

The pinned standard C++ executable itself fails with status 1 after its
diagonal phase: the earlier SWIPE phase leaves matrix adjustment enabled, but
the banded target iterator receives no adjusted matrices. Native execution
preserves that boundary after executing the preceding workloads. Full-count
stdout matches byte for byte, and stderr matches after replacing only timing
values. This is compatibility with the source failure, not a successful
completion of all benchmark phases or performance equivalence.

An instrumented copy of the pinned benchmark body, reduced to one iteration
and linked against available source-build objects, exposed an anchored score
mismatch: 79 in the previous native kernel versus 80 in that diagnostic.
The source prepends a masked target column and processes complete thirteen
column blocks using query-profile padding. Scalar and AVX2 native kernels now
retain those cells, including leading gaps. The reduced native workloads
match the diagnostic's anchored score 80, empty anchored pipeline result,
six SWIPE phases yielding 32 identical-target HSPs at score 124, and diagonal
endpoints 16/11. The diagnostic is not the immutable oracle; this distinction
and source/binary hashes are retained in `tests/data/benchmark/README.md` and
`probe.json`.

The C++ benchmark contains writes through destroyed bin-zero entries and an
out-of-bounds volatile score read. Native execution omits those undefined
accesses while retaining the actual kernel work and defined outputs.
Heterogeneous anchored lane grouping/refill, saturation, broader performance
fidelity, optional EXTRA benchmarks, and the remaining advanced command routes
remain open. The observed native anchored phase is substantially slower than
the recorded C++ run. The full translation goal is still active.

Evidence for this checkpoint is under `.tmp/native-benchmark-20261107/`.
The first full-count native release executable has SHA-256
`009f5a1362dd0d8cdfe21dd996bb74b1dda31f9dd4fd326ed97348c5715ab785`.
Its command comparison records exact stdout, timing-normalized stderr,
matching status 1, and all six exact invalid-option diagnostics. The debug
library suite passed 1,780 tests with six ignored before the additional
benchmark option-validation test, which separately passed. The final release
all-targets suite passes 1,780 library tests with six ignored and every runnable
integration target, including all twenty frozen upstream regression hashes.
The 700-file isolated package passes all-target/all-feature checking and both
real-kernel benchmark tests with all features enabled and no C++ source tree.
`provenance.json` records implementation, package, executable, and result hashes.
Scoped rustfmt, workflow actionlint, and diff whitespace checks pass.
The source inventory audit still reports the known missing
`LambdaCalculator.cc` map/audit entry; it is not counted as passing completion.

## Native bit-score cutoff integration (2026-11-07)

`blastp` and `blastx` now parse and forward `--min-score` to the existing native
`ScoreMatrix::report_cutoff` implementation. Nonzero thresholds replace the
reporting E-value cutoff; zero retains E-value filtering. The setting reaches
initial gapped scoring, final traceback, full-database alignment, and
alternative-HSP recomputation. Both SWIPE result passes and their reverse
parameter copies preserve it. `GappedScoreConfig::min_bit_score` also disables
the E-value-based ranking shortcut as in source. Protein and translated block
adapters clone the same complete configuration. Anchored-wrapper filtering
retains the pinned source's separate E-value behavior.

The new self-contained `native_min_score` integration suite compares 48 frozen
unchanged-executable cases across score-only, coordinate, and transcript
requirements. Protein CBS modes 0–5 and translated modes 0/1 cover zero,
retaining, fractional-boundary, partial, and rejecting thresholds against
`-e 1e-200`. It also checks one/four threads, actual query/reference blocks,
forced-disk hit buffers, and four byte-exact DAA outputs. Fourteen numeric
conversion cases match source and execute through `diamond::run`, including
hexadecimal prefixes, negative thresholds, NaN, infinity, and repeated scalar
options. Six malformed-argument diagnostics are frozen from the source.

Debug integration validation passes all three tests after fixing propagation
into final, full-database and alternative-HSP passes. The final release
all-targets run passes 1,780 library tests with six ignored and every runnable
integration target, including the new cutoff tests and all twenty upstream
regression hashes. The 756-file isolated archive passes all-target/all-feature
checking and all three cutoff integration tests with all features enabled and
no C++ tree. Twelve direct final-release pairwise/XML/PAF comparisons also
match the unchanged executable byte for byte. Ninety-six additional direct
source field-set comparisons verify that projecting the transcript fixtures
to score-only/coordinate columns is valid for these inputs. Earlier failing
diagnostic logs are retained separately from final evidence.

The final release executable SHA-256 is
`443765986c3c08d3db7288f057595c62027deddaac6bb0cba9973117a1f46faa`.
Fixtures and their source/input hashes are under `tests/data/min_score/`;
working evidence and complete provenance are under
`.tmp/native-min-score-20261107/`. Formatting and whitespace checks pass.
Remaining advanced/optional routes and the full translation completion
requirements remain open.

## Native tabular and flat JSON search contexts (2026-11-07)

The concrete `BlastTabFormat` callbacks now receive actual native search HSPs,
original query DNA, translated frames, original reference sequences, full titles,
mate sequences, FASTQ quality, taxonomy and source-phase self-alignment scores.
This replaces compact summary formatting in both eager and blocked protein and
translated search. Block dictionaries retain original reference sequences only
for fields that need them, and carry normalized-score scalars separately. Source
reference masking order differs between double-indexed and query-indexed search;
the self-alignment tables are captured at those distinct phases. Quality is
retained only when requested by the output fields.

All 52 standard registered tabular fields and 48 dynamic taxonomy rank fields
are connected. Dynamic IDs match source registration order 77–124. Flat JSON
`104`/`json-flat` and tabular `--unal` use the same native formatter. Source
missing separators and inconsistent taxonomy quoting are preserved deliberately;
the resulting flat JSON is not generally valid JSON. Native FASTQ parsing retains
full titles, accepts multiline quality and source-defined short full-quality
records, and rejects multiline sequence records. For C++ out-of-allocation
`qqual` slicing, Rust explicitly clamps both bounds rather than invoking undefined
behavior; no byte-parity claim is made for that source undefined read.

The self-contained five-test `native_tabular_context` suite passes 216 frozen
source output comparisons plus ten protocol/error cases. It covers protein CBS
0–5, translated CBS 0/1, FASTA/FASTQ, paired translated input, one/four threads,
small query/reference blocks, forced-disk hit buffers, TANTAN/SEG, both indexing
algorithms, multiple HSPs, taxonomy and unaligned callbacks. Field registry and
input/output hashes are recorded in `tests/data/tabular_context/`. Small failed
output files match the source by abandoning unflushed buffers; large failed-file
prefixes and general CLI diagnostic priority remain unverified. The oracle
executable remains unchanged. Final release and isolated-package verification
are recorded in `.tmp/native-tabular-context-20261107/`.

This closes the output-field routing gap for this registry, not the overall
translation goal. Other advanced controls, output formats and production fallback
routes remain open; FFI has not yet been reduced to testing-only use.

Final validation for this checkpoint: the release all-targets command exits zero,
with 1,780 library tests passing (six ignored) and all runnable integration targets
passing, including the five context tests and all twenty upstream regression
hashes. The 897-file isolated archive passes all-target/all-feature checking and
all five context integration tests with all features enabled and no `diamond/`
tree. Scoped formatting and whitespace checks pass. The final release executable
SHA-256 is `b339b261b671a662058a1c3fc7ab1da80fb97671363d18902f981742358839c6`;
source, fixture, log and package hashes are in the checkpoint's `provenance.json`.
New macOS ARM64/Intel package CI jobs are configured; their hosted execution is
still pending. The inventory audit remains open on `LambdaCalculator.cc`.


## Native SAM search integration (2026-11-07)

Search format `101`/`sam` now selects the concrete `SamFormat` backend for
protein and translated native search. It receives actual source frames and
HSP transcripts through both direct and intermediate output adapters. Headers
retain argv invocation and source mode. Public `SamOptions`, carried by both
search configurations and through translated/block clones, supplies invocation,
full subject titles and the optional source query-length tag. The CLI connects
`--salltitles` for SAM and `--sam-query-len`, with source-compatible rejection
of values on these boolean switches. Other-format full-title handling remains
separate work and continues to route explicitly to the legacy adapter.

Source block dictionaries discard SAM full titles even with `--salltitles`,
while direct reference blocks keep their loaded titles. The native SAM
intermediate format and dictionary insertion retain this distinction. The
existing CIGAR, MD, aligned amino-acid sequence, oriented source-coordinate,
frame, identity and score callbacks are used directly rather than reconstructed
from compact summary rows. Unaligned records follow the source format default
and block-dependent seed handling.

The self-contained debug suites pass all 168 frozen source CLI comparisons,
three additional in-process dispatcher calls, four malformed-flag diagnostics,
and two supplied translated-query block comparisons. Protein CBS 0–5,
translated CBS 0/1, both indexing algorithms, FASTA/FASTQ, TANTAN/SEG,
one/four threads and in-memory/forced-disk buffers are covered. Only the
invocation header field is substituted with the exact argv expected for each
native run; all other SAM bytes compare directly. Fixtures and hashes are in
`tests/data/sam_search/`. Final release/package evidence is retained separately
under `.tmp/native-sam-20261107/`. Other output formats/advanced controls,
frameshift search and the remaining production FFI dependencies are open.


Final validation for SAM integration: release `--all-targets` exits zero with
1,780 library tests passing, six ignored, and all runnable integration targets
passing, including all twenty upstream regression hashes and the existing
216-comparison tabular/JSON suite. The 987-file isolated archive passes
all-target/all-feature checking and the SAM/supplied-query integration suites
with all features enabled and no C++ tree. Dispatcher unit tests (23), scoped
formatting and whitespace checks pass. The release executable SHA-256 is
`08e44a55383983a4e052b204e81218115038e8a179a1458bb832a56d45100aa7`. Source, fixture, package and
log hashes are retained in `.tmp/native-sam-20261107/provenance.json`.
The unchanged oracle hash and the open `LambdaCalculator.cc` inventory gap
are recorded separately. Overall translation and production FFI removal remain
incomplete.


## Native edge and null search integration (2026-11-07)

`blastp`/`blastx -f edge` and `-f null` now select concrete native callbacks.
The format adapter has explicit DP requirements: edge requests COORDS; null
retains the source TRANSCRIPT pipeline and inherited no-op output methods.
The binary edge callback receives actual HSP coordinate coverages, original
query/reference OIDs and corrected bit scores after direct or intermediate
block processing. Edge/null avoid query-frame copies in output callbacks.
Native entry points and the supplied protein/translated-query APIs use these
same backends. Explicit unaligned settings are accepted without emitting rows
for these formats. Existing transcript formats retain their own defaults.

Raw binary comparisons exposed corrected-score ULP differences previously
hidden by decimal formatting. The unchanged `-march=native` source executable
contracts multiply/add operations in the final corrected score and ALP
log-area calculation. Explicit Rust `mul_add` now preserves that rounding,
and the score-matrix adapter preserves the source evaluator's subject/query
log-area dimension order. This is a compatibility decision against the pinned
build, not a claim that different C++ compiler flags yield identical doubles.

A 60-query varied-score/length/frame probe additionally exposed translated
coverage rejection in both native reverse-pass adapters. They used aligned
amino-acid length divided by DNA query length. Both now use the existing
source-coordinate `query_cover_percent` method, matching source. The 50% equal
coverage regression now retains the same 120 records per source case; the
unfiltered cases retain 180. Original IDs are checked under length sorting.

Fixtures in `tests/data/edge_null/` cover 196 byte-exact source-output
comparisons: 168 CLI matrix cases, 24 varied-score/sorting cases and four
supplied-query API cases. Three dispatcher repeats and two malformed-query
diagnostics supplement this matrix. Null runs also require observed seed work
rather than accepting a trivial early empty-output return. Focused checks
pass; final release/package verification is retained separately under
`.tmp/native-edge-null-20261107/`. Taxonomy format 102, advanced search
controls, other production fallback dependencies and the remaining inventory
audit are still open.


Final validation for this checkpoint: release all-targets exits zero, with
1,780 library tests passing (six ignored) and all runnable integrations passing.
This includes the new edge/null tests, supplied-query APIs, 216 tabular/JSON
comparisons, SAM comparisons and all twenty upstream regression hashes. The
1,093-file isolated archive passes all-target/all-feature checking and all eight
edge/null and supplied-query integration tests with all features enabled and
no C++ tree. Scoped formatting and whitespace checks pass. Release executable
SHA-256: `48b1b96cbff7d377884bfab672f8fc0d823f3a129c54498e80f594997372e8c7`. Source, fixture,
package, disassembly and final log hashes are retained in the checkpoint's
`provenance.json`. Earlier failed/partial probes are retained separately.
The known inventory gap and remaining production FFI routes remain open.

## Native taxonomy classification integration (2026-11-07)

Search format `102` now invokes the concrete TaxonFormat accumulator for
protein and translated queries. Actual matched subject taxid lists, node
parents and conditional scientific names survive direct and intermediate
block output. A new accumulator is created per query, preserving minimum
mapped E-value/LCA state and source unaligned reporting. Its DP requirements
are NONE unless filters need additional values; output avoids query sequence
copies. Both CLI and supplied protein/translated-query APIs use this path.
`--include-lineage` is a native boolean setting and names are required only
when it is used; rank metadata is not requested for classification.

The source classifier implicitly enables a 10% top-score range when neither
an explicit range nor a nonzero bit-score cutoff is supplied. A shared native
normalizer applies that once, outside captured block-core calls. Explicit
`--top 100` clears the range and selects unlimited targets; nonzero
`--min-score` suppresses the default. CLI top values use the existing C++-style
float parser, validate range and reject explicit --top/-k combinations. A
previous debug overflow in ranking-chunk rounding for unlimited targets is
fixed by calling the existing translated overflow-safe `make_multiple` helper,
as the source does. The weak-reference fixture distinguishes these behaviors
at the classification level and through both supplied-query APIs.

Focused tests pass 284 direct source-output comparisons, three dispatcher
repeats and ten metadata/option errors. Eight additional checks cover a source
undefined-behavior exception, counted separately: original protein taxonomy
runs with equal 50% coverage can index an unloaded reference-title array in
`Block::length_sorted`, producing allocator aborts/signals or vector errors.
The native path retains valid metadata and computes the defined filter/LCA
result. For these fixtures only, source `--no-self-hits` forces title loading
without excluding any query/subject pair, providing a successful reference
for that native result. Original failing commands/diagnostics and alternative
source commands are preserved in `tests/data/taxon_search/`. These eight
checks do not establish byte parity with an aborted source invocation or
prove that the workaround is equivalent on arbitrary data. No C++ source
patch or rebuilt oracle was used; ptrace was unavailable in this environment.

The standard source search-output registry is now connected natively. Advanced
search controls, optional modes, the unmapped `LambdaCalculator.cc` audit and
remaining production fallback routes remain open. Final release/package
verification is retained under `.tmp/native-taxon-20261107/`.


Final validation for taxonomy integration: release all-targets exits zero with
1,780 library tests passing (six ignored) and every runnable integration target
passing. This includes all twenty upstream regression hashes and the prior
SAM/edge/null/tabular/JSON suites after shared top normalization and unlimited
chunk-rounding changes. The 1,251-file isolated archive passes all-target/
all-feature checking and all nine taxonomy/supplied-query integration tests
with all features enabled and no C++ tree. Scoped formatting and whitespace
checks pass. The release executable SHA-256 is
`7218dad816e5c212eb9ac8f65d7e24727bf2f0e5919ed3c531b5d12ac7f675ae`. Source, fixture, package and
final-log hashes are in `.tmp/native-taxon-20261107/provenance.json`; source
crash probes are distinguished explicitly from direct passing conformance.
The overall translation goal remains active.


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


The output-controls checkpoint passed `cargo test --offline --locked --release
--all-targets`: 1780 library tests passed (six ignored), along with all integration
suites and the twenty upstream regression hashes. The isolated 1614-entry crate
passed an offline all-target/all-feature check without the adjacent C++ tree.
The immutable source executable hash and exact source/fixture snapshot are in
`.tmp/native-output-options-20261107/provenance.json`. Later parallel edits are
tracked as a separate batch; this evidence does not validate those later changes.
MacOS default stock-zlib bitstream compatibility remains a hosted test gate;
linking different system zlib versions is not proof of identical compressed bytes.

The isolated packaged runtime suites also passed with all features: seven native
block-library tests and three output-option tests, including 348 CLI comparisons.
This verifies that this checkpoint needs no adjacent C++ checkout at runtime.


## Parallel numerical and search-control batch (2026-11-07)

`--approx-id` and `--soft-masking` now flow through the typed protein/translated
configs, ordinary dispatch, blocked/captured search and supplied-query APIs.
Approximate thresholds reach Hamming cutoff selection, score/coordinate DP
requirements, reverse-pass rejection and final approximate-identity filtering;
explicit clustering overrides retain precedence. Soft masks combine automatic
motif with explicit SEG/TANTAN using the actual score-matrix Masking object.
Prepared tables reach hashed Blocks and authoritative translated contexts,
restore residues before extension, and preserve seed annotations across chunks.
The source option spellings/count/repetition behavior remains explicit.

Focused validation passes 136 direct CLI comparisons plus 16 supplied-query
comparisons, nine parser cases and the impossible typed buffer check. TANTAN
changes emitted outputs for both query types and both algorithms; custom
BLOSUM45 masking and threshold-driven coordinate requests are exercised.
Six original source configurations crash in SEG; bounded gdb catches signed
byte 0x8a becoming index -118 into 128-entry alphabet tables and then corrupting
composition indexing. Source hashes remain unchanged. Of 18 safe-variant probes,
17 complete and one larger-block variant also crashes. All 6 algorithm-0
snapshots agree with independently completing softnone variants; 12 native
checks use those corroborated results as a defined compatibility choice.
They are not direct parity with the crashing invocations. Larger-block variants
that change approximate-identity fields are explicitly rejected as references.
See tests/data/search_controls and the source-seg-crash evidence beneath
.tmp/native-search-controls-20261107. Existing reassign/recluster and supplied
block-library focused regressions pass; source's repaired diagnostic cases
retain their prior separate classification.

The outstanding LambdaCalculator.cc row is now audited and mapped. Its active
Rust implementation had arithmetic decimal rounding and uncontracted LU/solve
products. Six-significant-digit formatting/parsing and explicit FMA match all 36
immutable-object state vectors, full probability bits and RNG continuations;
19 formatter probes include halfway values, subnormals and signed special
values. Seven baseline states differed, including a different valid lambda
and RNG consumption. Cargo focused tests now pass without C++ at runtime.
Allocation/deletion/swap/transpose/sum mappings and defined invalid-input
safety decisions are recorded in tests/data/tantan_lambda/README.md and the
file map. The inventory script passes 199/199 with no missing modules or rows;
this is an inventory result, not an application-completion claim. Exact target
math-library behavior on macOS/Windows remains a hosted gate.

Configured buffer conversion is checked; impossible capacities receive a
fallible allocation preflight in CLI and typed search before stdio construction.
Negative/maximum/saturated values retain source allocation diagnostics and
create no output. This does not guarantee subsequent allocations under changing
memory pressure. The output-option protocol suite passes 23 cases on 64-bit
(24 on 32-bit), with bounded original probes retained in the fixtures.

The full combined release/package checks are recorded separately after the
frozen focused snapshot. Global ranking, iteration, title/archive policies,
translated minimum-length/frameshift/swipe controls, BLAST-database source-UB
policy, optional build modes and other registered options remain open.
Production automatic fallback still exists; the overall goal remains active.


The final combined `cargo test --offline --locked --release --all-targets`
passes: 1781 library tests (six ignored), all integration targets and all 20
built-in source regression hashes. The isolated 1719-file crate passes a fresh
all-target/all-feature check with no adjacent C++ tree. Its all-feature runtime
checks pass seven block-library, two output-option, four search-control and one
lambda integration tests. The 348 output-option CLI comparisons additionally
pass in the full release run; that long test is deliberately omitted from the
isolated focused rerun. Exact source/package/binary/fixture/log hashes and
platform limitations are recorded in
`.tmp/native-search-controls-20261107/provenance.json`. Scoped formatting,
whitespace and the 199-file inventory audit pass. This closes this batch only;
the end-to-end goal and production FFI-removal requirement remain active.


## Native archive dictionary policy and registration inventory (2026-11-07)

ArchiveOptions connects the translated DAA flags and raw-title/all-ID/first-ID
selection at all three native dictionary insertion sites. Typed protein and
translated configurations and block/spill adapters retain it. CLI salltitles,
sallseqid and DAA unal0/1 no longer route to C++. Source full-title precedence,
actual alternate-header separators and aligned-only DAA records are preserved.
The focused archive suite passes 156 CLI comparisons, 48 archive/view outputs,
24 supplied-query archive comparisons and four protocol errors. Existing seven
block-library and two SAM regression tests also pass. The eager BLAST reader
still has its separate source-UB gate; dictionary code is updated there without
claiming new eager-path runtime parity. See tests/data/archive_titles.

The new reproducible registry extraction enumerates all source registrations,
not just exposed help names: baseline 21 commands/299 options (161 enabled,
138 hidden disabled); WITH_DNA adds 17 options; EXTRA adds 24 commands and enables
hidden fields while disabling workflow checks. All 316 names match the Rust
registry, with no unknown arities. 836 bounded baseline parser/doc probes pass.
Defaults are registration expressions; post-store normalization is not inferred.
The final short map includes disabled seed-weight shadowing enabled window's
-w alias. Alternate compiled variants are modeled statically, not runtime-proved.
This closes enumeration gaps but does not prove native option semantics;
optional build modes and remaining native integration gates stay open.

Global-ranking prerequisites now preserve context in table-based final extension
and source uint16 score narrowing (65536->0), with an isolated original-header
probe and 21 passing focused helper tests. The concrete native ranking lifecycle
is being integrated separately; helper tests alone are not CLI parity evidence.


## Code-first native integration checkpoint (2026-11-07)

Per the user's instruction, the current phase prioritizes code translation and
concrete backend integration. Broad runtime/parity/package/platform runs are
explicitly deferred until that integration is connected. Compilation is not
being treated as proof of parity or completion.

The following are now connected in source and pass the shared all-target Cargo
check plus a native debug executable build:

- Global ranking: per-shape/per-reference persistent ungapped ranking tables,
  original OIDs and translated contexts, ranked-union reload, actual parallel
  full-matrix extension and ordered native output. The source uint16 narrowing
  and source-required output restrictions remain explicit.
- Iteration: source schedules including real linear rounds, stable source-query
  skip snapshots, actual extension-based aligned tracking, reference rewinds,
  persistent dictionaries/intermediate streams and one final block join.
  Global-ranking iteration writes real intermediate records.
- SWIPE: exhaustive every-query/every-target extension without synthetic seed
  hits, source score/final passes and filter-aware culling. Its interactions
  with iteration, global ranking and soft/motif settings follow source branches.
- Frameshift: concrete legacy mapper and banded pipeline backends, actual
  matrix penalty, grouped three-frame DP and overflow handling, strand/source
  coordinates/transcripts, legacy empty-hit/null behavior, and native output.
  Banded scoring now uses packed band storage and rolling affine states rather
  than allocating full query-by-subject matrices; traceback retains the actual
  affine-gap path and distinct start/end frame positions.
- Ordinary extension modes: typed and CLI selection of banded-fast/banded-slow,
  full, global and none reaches the existing concrete kernels. Explicit
  clustering overrides retain precedence. Seed-only protein coordinates and
  filter bypass match the corresponding source path.
- Contiguous seeds: actual dense bitset/query filtering and contiguous encoding,
  source shape/reduction setup and partitioning, contexts/minimum/skip support;
  this is not hashed mode relabeled as ctg.
- Minimum query length: signed source parsing and global context-aware filtering
  through counting/filling/sketch/histogram/hashed/contiguous paths on both sides.
  The original's incomplete six-context reference-group out-of-bounds case is
  handled safely by excluding incomplete positive-threshold groups.
- BLAST references: bounded native record chunks/rewinds/random reads, volumes
  and PAL masks with sparse original OIDs, conditional taxonomy/title metadata,
  filtered scoring letters versus full reporting totals, and native CLI routes.
  The previous source StringSet memory-corruption boundary uses safe decoded
  Rust blocks; no direct parity is claimed against a crashing/corrupt process.
- CLI protocol: a self-contained baseline registration preflight uses the
  translated parser for real arity/workflow/disabled/short-shadow behavior.
  Validation does not pretend the corresponding option backend is implemented.
  Native convenience flags remain explicit; original command-prefix invocation
  metadata is preserved. Getseq now has its actual typed baseline parser,
  last-wins --seq semantics and streamed DMND/BLAST/FASTA sources.

Evidence: `.tmp/native-global-ranking-20261107/expanded-native-compile.log`,
`native-code-build.log` and `code-connection-snapshot.json`. The pinned C++
executable remains unchanged. Runtime correctness of these new integration
paths is not yet established by conformance runs. Earlier snapshots retain
only the bounded checks actually recorded for them.

Remaining work includes advanced search/culling/filter/scoring/output controls,
query/target-index combinations, clustering controls, optional build modes and
EXTRA commands, and the final parity/platform/package and performance gates.
Automatic production legacy routing still exists for untranslated valid modes;
removing it without replacing their functionality would not satisfy the goal.
The full original objective remains active and unachieved.

## Search controls and optional tools integration (2026-11-07)

This is a code-first checkpoint, with broad conformance deliberately deferred.
Native search now consumes the shared Source setup controls rather than merely
accepting their names: custom scoring matrices, stop score and database-size
calibration, shapes, seed/frequency/sketch/bin settings, configurable Stage0
and Hamming tiles, short/long ungapped filters, gapped filter thresholds/windows,
ranking, bands, anchored extension, DP task budgets and matrix-adjust thresholds.
A command-local matrix snapshot avoids recalibrating custom matrices per block.

Native DMND/BLAST filtering preserves original OIDs and implements taxonomy and
sequence-ID selection with Source precedence and missing-ID policy. Range and
species culling reach actual target/intermediate joins. Auxiliary aligned and
unaligned query files use real sequence/quality data and alignment tracking
across reference blocks and iteration rounds. Retained temporary files and
first-block protein-input warnings now consume their corresponding controls.

Reassignment and recursive reclustering forward the shared controls through
subdatabases, masking, search and realignment; ordinary clustering also receives
the shared tuning and matrix cache. All cascade budget calls consume the real
no-block-size-limit setting. Source linclust-chunk-size has no runtime read in
reassignment or recursive reclustering: its only read is external linclust's
chunk table (external.cpp:156); recursive recluster.cpp:128 calls cascaded
directly. Treating it as a no-op here preserves Source behavior. External
linclust and multiprocess search remain separate unfinished adapters.

Default-off Cargo feature `extra` selects Source's EXTRA registration profile.
The native tools backend connects random-seqs, info, mask, fastq2fasta, reverse
and smith-waterman. The last command faithfully reaches Source's error:
ScoreMatrix("DNA") is unknown in this checkout; a DNA matrix is not invented.
Source commented-out input reads in several other tools leave infinite loops;
the native readers terminate at EOF. This defined compatibility behavior is
explicitly distinct from literal nontermination parity. GNU random arithmetic
uses explicit wrapping where Source signed overflow is undefined.

CLI scalar and vector selection now preserves Source's last-occurrence
assignment behavior, including query files and output field lists. This closes
an implementation gap in the old convenience parsing helpers.

These additions do not establish runtime parity or complete the goal. Optional
commands, external/distributed backends and final platform/package/conformance
gates remain open; production fallback removal must follow their replacement.

The combined `cargo check --offline --locked --all-targets --features extra`
passes in 16.40 seconds; retained evidence is
`.tmp/native-options-20261107/compile.log` and `code-snapshot.json`. The same
snapshot includes native split/hashseqs/listseeds/wc/cut and deprecated roc-id
routing, using real generic database readers, masking/seed enumeration, hashing,
compressed chunk outputs and TSV processing. Existing tests were type-checked,
not executed. No new broad runtime or parity results are claimed.

## Distributed and optional startup code batch (2026-11-07)

External linclust now has a concrete native pipeline, rather than a trait-only
port. Its manifest reader, seed/pair/radix/chunk stages, real DP, closure and
representative output connect to the CLI's parallel-tmpdir workflow. Shared
queues, stage barriers and closure ownership coordinate workers; compressed
writers and volume/pair readers stream data. Radix partitioning follows Source's
single eight-bit spill, including retained parent files and empty children.
Skewed child buckets may exceed the Source memory-budget heuristic; recursively
partitioning them would change intermediate order and chunk boundaries, so the
native translation deliberately retains Source behavior. Runtime parity and
multihost/recovery evidence are still deferred, not established by this code.

Multiprocess search now connects real native DMND partition loads, OS-locked
work/recovery queues, named block streams and query output chunks. Its on-disk
dictionaries use Source's OID/length/title/optional-sequence/f64 records with
implicit per-block IDs; joins reload and remap these dictionaries. This batch
requires a shared compile check and later conformance/recovery runs before any
claim of verified distributed behavior.

Optional commands now execute the translated basic/config startup through a
concrete owned backend: scoring, masking, CBS mode, genetic-code initialization,
temporary-directory handling and string-enum validation have real dependencies.
Source EXTRA's seventh CBS mode and named matrix scale are exposed explicitly.
Mask/listseeds use the initialized matrix rather than a hardcoded BLOSUM62
masker. Generic config now forwards the long-read-adjusted frameshift penalty
into its scoring callback. Source-declared commands without a run/main switch
arm return Source's status 1 natively after startup; they are not assigned an
invented algorithm. The registered fallthroughs are sort, dbstat, seed-stat,
simulate-seqs, upgma, upgmamc, compute-medoids, mutate, model-seqs and
make-seed-table. Optional find-shapes now reopens and streams its actual input
for each greedy selection round and emits Source's progress/pattern output.
Its exhausted zero-pattern formatter is undefined in Source; native code reports
an explicit error at that boundary. Full optional startup diagnostics/help and
runtime parity remain validation requirements.

Shared `cargo check --offline --locked --all-targets --features extra` passes
for the distributed/startup batch (16.39 seconds). Evidence:
`.tmp/native-options-20261107/distributed-compile.log` and
`distributed-snapshot.json`. This is type-checking evidence only: broad runtime,
source parity, crash recovery and platform/package checks remain deferred.

## Linear/output and actual KSW dependency batch (2026-11-07)

Source linear query/combo modes now reach the real Stage0 packed-location
WorkSet dispatch and its existing ranking tables, mutual-coverage and longest
combo kernels. Source keep-target-ID and query sorting rules are carried into
this path. Actual setup extension selection now controls linear DP, so
linclust-banded-ext is no longer overridden by an unconditional full-matrix
choice; explicit clustering overrides retain precedence. The original blastx
query-linear length-sort path can index six-context sequence data using
single-context ID/OID metadata; native code keeps source-query groups intact
at this undefined Source boundary, without a direct parity claim there.

Search qnum/snum offsets (including unsigned wrapping of signed offsets), XML
BL_ORD IDs and raw XML accessions now reach real formatters, including global
ranking. Source-ignored search options no longer require C++: deprecated
constructor-local scalars and test/ROC/clustering-only fields retain typed
native parsing with no invented search effect. The explicit legacy diagnostic
switch now bypasses all native routes consistently. Ordinary clustering forwards
actual ext-chunk-size ranking and Source-calculated index partition counts.

The DNA scalar KSW wrapper now uses an actual KSW recurrence with bands, peak
coordinates, z-drop, traceback and Source flags instead of a generic affine-DP
substitute. Canonical upstream revision and hashes are retained under
`.tmp/native-dna-20261107/`; the dependency's MIT license is retained in code.
The adjacent DIAMOND checkout lacks the referenced KSW/WFA dependency trees, so
identity with its absent dependency revision is unproven. GNU's SSE variant,
WFA and native blastn remain unfinished; the scalar addition is not a claim of
DNA CLI parity. Source zero-pattern formatting and other previously documented
undefined cases retain explicit defined-native boundaries.

Shared `cargo check --offline --locked --all-targets --features extra` passes
in 12.43 seconds (`linear-output-compile.log` and
`linear-output-snapshot.json` under `.tmp/native-options-20261107/`). This
compiles tests without running them. No broad runtime, parity, packaging or
performance result is implied.
The same source also passes `cargo check --offline --locked --all-targets`
with default features in 16.65 seconds; `default-compile.log` retains that
result. Both checks keep the C++ feature off and execute no runtime tests.

## Actual block/loader/database and GNU KSW code batch (2026-11-07)

Clustering's native cascaded and centroid adapters now schedule real
query/reference blocks using Source letter budgets, self-search alias rules,
length sorting and original OID composition. Results pass through existing
on-disk intermediates, reference dictionaries and the actual merge/culling
logic. Whole-database sorting and the former one-shot in-memory search are no
longer substituted for this workflow. Edge intermediates use Source's actual
EDGE=9 tag and coordinate flags; the CLI still does not accept numeric format 9.

Native BLAST loading now uses a bounded rawchunk producer, parallel decoders
and ordered consumer, with actual load-threads/minichunk controls and Source
raw-letter chunk boundaries. Compressed hit packets use the same ordered worker
pipeline. Source-invalid zero/negative worker situations that can deadlock are
handled with a defined positive worker count. Legacy frameshift rank-ratio
and unaligned-targets are connected to real pre-join extension matches, OIDs
and a duplicate-preserving accession arena. BLAST's absent Source accession
map is an explicit undefined-source compatibility boundary: native code uses
actual BLAST record accessions there, not an invented source-parity result.

The `makedb` command now forwards actual buffers, masking matrix and
probability, ignore-warnings, sequence type, suffix policy and log/quiet state;
missing --in uses the actual stdin reader. Local thread pools replace mutation
of Rayon's global pool. The builder's streaming protein warning rejects the
tenth DNA-only record before parsing the next record, or a shorter DNA-only
first block at block end. Its Source soft-mask worker is TANTAN bit masking
even though make_db passes the SEG flag: mask_seqs ignores that algorithm when
hard_mask=false. The translation retains this actual Source branch instead of
changing database bytes to SEG. Legacy in-memory builder convenience APIs keep
their existing warning policy; native command runtime uses Source policy.

The GNU KSW engine now follows the real anti-diagonal byte-difference state
machine, signed/unsigned lane operations, saturation, approximate maximum path,
z-drop, end bonus and rotated traceback; Apple selects Source's scalar variant.
A scalar affine DP is no longer used as a substitute for the GNU KSW path.
Missing DIAMOND dependency revision, WFA and native DNA reporting/CLI remain
unresolved; canonical-source provenance alone does not prove pinned-fork parity.

Shared `cargo check --offline --locked --all-targets --features extra` passes
in 14.97 seconds. Evidence is `.tmp/native-options-20261107/loader-cluster-*`.
This is compile evidence, with broad runtime/parity testing still deferred.
Heartbeat's live output/scheduling counters remain an active translation task.
Invalid output/algorithm values and invalid target-index combinations now reach
native validation rather than the C++ diagnostic fallback.

## Native routing, shared configuration and DNA command (2026-11-07)

Automatic production C++ execution has been removed. The executable calls the
optional comparison oracle only for an explicit `--legacy` request; normal
unconnected requests fail natively and still count as unfinished work. The
Cargo `ffi` feature remains default-off and needs the adjacent source checkout.
This satisfies the routing requirement, not whole-program parity/completion.

Source-wide command logging, real debug log sinks and restoration, protein-mode
strand/genetic-code validation, shared ORF parsing, scoring scale and additional
DP/legacy settings are connected. Generic dbinfo now opens DMND or BLAST with
Source's NO_FASTA/NO_COMPATIBILITY_CHECK behavior. Database detection now honors
Source's existing .dmnd candidate preference even when the original path exists.
A common parsed search-config factory replaces the duplicated protein literal
and supplies the DNA command without inventing a second option registry.

Header2 now implements Source EXTRA's db_type tail and C ABI declared record
size separately from its physical serialized payload: on 64-bit targets Source
writes prefix56 then52 payload bytes, without the four padding bytes. The writer
uses actual physical offsets; the reader retains Source's declared-record cursor
and implements DynamicRecordReader's per-field zero fill and remainder reads.
Baseline48-byte headers remain unchanged. This is source-code correspondence;
feature-specific database/runtime parity is still pending.

Default-off `dna` implies `extra` and selects EXTRA+WITH_DNA registration. The
native command owns real nucleotide reference blocks, actual minimizer indexing,
both-strand KSW/WFA and fallible chained extension, validated Blastn_Score,
streamed query/reference chunks, ordinary intermediate joining, taxonomy/OIDs,
auxiliary output and nucleotide rendering. Actual Source run/config resolves
DNA gap defaults to5/2 (with reward2/penalty-3), correcting the earlier incomplete
constructor-only audit; raw unsupported explicit combinations still fail.

Source DNA behavior retained includes reference-only min-query-length filtering,
re-search of already-aligned queries during partial iteration, direct target
order, blocked frame-context loss, archive reward/penalty0/0 from distinct
unused config fields, and global-ranking's defined empty table/result lifecycle.
Global scoring/context reads left uninitialized by Source are instead backed
by actual DNA statistics/owned contexts, explicitly classified as a safe native
boundary; no direct pinned-fork DNA parity is claimed. Canonical missing KSW/WFA
dependencies remain separately pinned, with DIAMOND fork identity unresolved.

Shared all-target checks pass for the header/writer/common configuration batch
with `extra` (10.51s) and the command/engine batch with `dna` (20.42s), retained
under `.tmp/native-options-20261107/{header-writer,dna-command}-*`. These are
compilation results, with tests not executed. Query stdin, AUTO selection,
DNA MPP/timers, optional MCL and remaining Source callgraph findings continue;
all final conformance/platform/package/performance gates remain open.

### Further concrete callgraph integration (2026-11-07)

The AUTO/DNA-MPP/MCL/getseq batch passes all-target `dna,mcl` (20.11s) and
default-feature (17.08s) checks. The later lifecycle/workflow batch passes
all-target `dna,mcl` (16.15s). Logs are under
`.tmp/native-options-20261107/`; `translation/progress.md` identifies the exact
batch contents. These are compile checks only; broad tests remain deferred
while translation continues, as requested. In particular, neither optional
feature parity nor final hosted macOS/package/performance gates are proven.

Source-call audits continue to distinguish real gaps from dead alternatives:
merged-list global ranking's only Source invocation is commented out, and the
packed-location k-mer constructor's unsupported result is Source behavior.
The native Source-like entry-point facade still needs a concrete prepared
backend used by actual dispatch; that refactor is underway. Prepared input,
real shared tuning, live DNA/MCL paths and bounded tool dependencies are being
translated before expanding runtime comparisons.

### Current code-first checkpoint (2026-11-07, subsequent shared batches)

The earlier entry-point/backend refactor mentioned above is now connected to
actual prepared native commands, shared initialized Source resources, regression
sessions and real benchmark/tool dependencies. Production dispatch has no
automatic C++ route; the explicit default-off FFI oracle remains a separate
conformance option. This is not a claim of complete translation or parity.

Further code batches connect actual radix/hash recursion, task scheduling,
parallel DP, chaining/band controls, frameshift chunking/reservations, self-query
DB ownership, CBS window/motif caps, minimum length ratios, Hamming filters and
clustering thresholds. Exact compile log names/results are in
`translation/progress.md`; broad runtime suites remain deferred by request.
FASTA owner fidelity remains active code work, and full command/option callgraph,
feature-profile parity, hosted Apple/package execution and performance gates
remain open. File/name coverage and successful compilation do not close them.


### Subsequent native integration checkpoint (2026-11-07)

The earlier FASTA-owner work is now implemented and compiled. Concrete master
query/reference iteration, seed/global-ranking sorting, DNA defined-state/output
semantics and real worker/join/output statistics have since been connected.
`translation/progress.md` records batch scope and compile evidence. These changes
remain subject to runtime parity verification; the goal is active and incomplete.

The historical DNA dependency identity gate described in older checkpoints is
superseded by read-only recovery of the actual vendored KSW/WFA trees and an
identical pinned caller tree from Git history. The manifest and active-path
comparison are in `.tmp/native-dna-20261107/`. Optional DNA runtime parity is
still open. Prior broad claims of uninitialized global DNA score metadata or
all nonzero-frame accesses being undefined are corrected: static metadata is
zero-initialized and default sequence contexts are defined empty; only actual
invalid accesses/conversions require a safe compatibility boundary.

Production CLI/library execution has no automatic FFI route. The explicit oracle
feature remains off by default and can use an isolated source directory. No
original adjacent Source file/object/binary was rebuilt. All-target compilation
is a development check, not completion evidence. Full command/function runtime
parity, optional profiles, hosted macOS/package execution and performance gates
remain required. DoubleIndexedBackend is an orchestration model with a test
adapter; concrete production drivers, rather than facade/name coverage, determine
native integration status.


### Viewer/cluster/frameshift checkpoint (2026-11-07)

Further real native code now preserves DAA viewer first-record/header ordering,
Source three-frame zero-score behavior, serialized cluster realignment streams,
real DP failure/control propagation and final realignment reports. Actual search
RSS and hit-loading/extension timing effects are connected. The combined/default
all-target batch compiles (11.30s/11.40s); exact evidence and scope are in
translation/progress.md. Broad runtime testing remains deferred while translation
continues.

The registry inventory requirement is now independently verified against the
current tree, including both registration preprocessor dimensions (EXTRA and
WITH_DNA); other feature flags introduce no Config registration branches. This
closes inventory only. Command/option behavior and parity remain open. The
199-file module audit likewise passes and is structural evidence only.

The backend audit explicitly identifies AlignBackend and DoubleIndexedBackend
as still lacking production adapters; Source's multinode call is inactive.
Actual concrete driver behavior and generic model tests are distinguished.
A coherent production alignment adapter refactor is underway, preserving real
HSP ownership, bounded grouped hit storage and the scoped ordered executor.
Completion remains unproven until this integration and all required runtime,
feature-profile, hosted platform/package and performance gates pass.

### Actual production adapters connected (2026-11-07)

The earlier disconnected alignment backend is superseded by a real grouped
native AlignBackend with actual bounded loading/execution/output owners. Shared
owner-based double-indexed orchestration now drives ordinary protein/blastx,
DNA and MPP search through concrete native adapters; its old metadata interface
remains explicitly a compatibility model. Both combined-feature and default
all-target batches compile; exact scope/results are in translation/progress.md.
Provided-owner clustering is still being moved onto that shared orchestrator.

The optional DNA audit additionally establishes Source unfinished-writer failure
and dangling/mismatched ref-buffer ownership. Native safe ownership and real
translated dependencies are retained with explicit compatibility classification;
this does not establish successful Source DNA CLI parity. The full completion
goal remains active, including all required runtime/profile/platform/package and
performance evidence.

### Known orchestration and accession-owner code gaps closed (2026-11-07)

Provided-owner clustering now uses the actual shared owner orchestrator too.
Its physical reporting metadata, compact lookup domain and filtered scoring
metadata are kept distinct. Unaligned targets use packed bits and a real disk
accession spool through Source's constructor/final-table lifecycle; the interim
live-lookup limitation is superseded. Combined/default all-target compilation
passes (10.42s/10.55s). The qualified backend inventory now has native adapters
for all 14 active contracts; two retained interfaces are explicit compatibility
models and Source multinode invocation is inactive.

These are concrete code/integration and compile results. They close the known
adapter inventory gaps, not the required runtime/profile/package/platform or
performance gates. Fresh parity work is the next phase; the full goal remains
active and completion is unproven.


### Runtime-discovered fixes and exact MCL numerical translation (2026-11-07)

The first integration batch passes all 20 upstream built-in regression hashes,
but eight integration targets failed and the second default library run still
has two failures (1785 passed / 2 failed / 6 ignored). Real shared-search masking,
query self-score ownership and CLI stdout defects are fixed in the pending
shared build; stale fixture assumptions are corrected against Source probes.
Numerical MCL dependencies now use translated Eigen 3.3.8 arithmetic, with 54
standalone exact Source vectors passing. This is numerical-stage evidence only:
full pinned MCL CLI fails to build against removed Source transport APIs.

The latest fixes still need shared Cargo integration. Remaining runtime,
optional-profile, package/platform and performance requirements stay open.
The goal remains active and unachieved; details and evidence are recorded in
translation/progress.md and .tmp/native-post-translation-20261107/.


### Seed-table/self-hit corrections and platform coverage (2026-11-07)

The actual seed-table override, annotated reduction and generic stage-2 self-hit
corrections compile together with dna,mcl,zstd (7.42s). All eight cluster command
unit cases, the reduction invariant and native regression logging now pass.
Recursive clustering advances through the original failed fixture but still has
a natural approx90 byte mismatch; assertions are unchanged and the goal remains
active. Current detailed commands/statuses are in translation/progress.md and
.tmp/native-post-translation-20261107/masking-reduction-focused-results.json.

CI now runs optional native-feature tests from the isolated package on both
Apple architectures as well as Linux. No hosted green result is claimed.


### Recursive matrix and isolated native profiles (2026-11-07)

All three recursive-clustering integration tests now pass (82.29s), after
reconnecting the facade's already documented compatibility counter to its actual
provided-owner driver. Source's uninitialized first-read qualification remains;
no new heuristic or weakened expectation was used. Combined compilation passes
in 10.17s. Isolated packaged native combined and EXTRA-only profiles build without
C++ sources. Fresh EXTRA captures exposed cut stdout and file-open diagnostics
still requiring correction. Completion remains unproven and the goal active.


### Current I/O and optional-profile checkpoint (2026-11-07)

Actual typed file-open/temporary errors, Source directory policy and genuine
pre-reset temporary initialization are connected through the native drivers.
GVC compression/ordering/ostream behavior and clustering accession-table
ownership/cleanup now follow Source. Paired input/output-error and GVC cases
pass with only explicitly qualified random temp suffixes or build-profile
assertions. The CBS CLI count and positional EXTRA DMND reader issues are fixed.

Default library tests pass 1792 with six ignored; dna,mcl,zstd library tests pass
1791 with six ignored. Default/combined/Windows GNU all-target compilation passes.
The fresh exact baseline Source independently passes all 20 regression cases;
its bounded auxiliary-memory failure remains qualified. Fresh native release,
full integration/profile and hosted platform/package/performance verification
remain open. Full details and limits are in translation/progress.md.


### Fresh native release core parity (2026-11-07)

The corrected default native release passes all 20 built-in regression cases
under the same 4 GiB bound as the fresh pinned Source (status 0; no leftover
files). This follows real fixes to automatic reservation policy, actual Source
sequential getseq and the trace-free SWIPE DummyRowCounter path. Current library
runs pass 1795 default and 1794 native-feature tests, each with six ignored.
Exact database bytes and bounded EXTRA roundtrip branches also match. Scope and
auxiliary-memory qualifications are retained in translation/progress.md.

A final terminal/log ordering correction is undergoing focused checks. Broader
integration, optional-profile and hosted platform/package/performance gates stay
open; the full goal remains active and completion unproven.


### Stable selected parity checkpoint (2026-11-07)

All 12 Source SWIPE comparisons and the 26 existing standard-input/GVC cases
pass. EXTRA databases and seed indexes are byte-identical, with matching crossed
getseq/dbinfo/SWIPE/target-indexed results. Logged allocation and temporary-file
failure ordering is corrected and verified, including the real pre-call stream
prefix. Inventory/registry checks pass; Source hashes remain unchanged.

These are selected bounded results. Broader current integration, package,
hosted platform and performance gates remain open. The goal stays active.


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


Worker/index translation checkpoint (2026-11-08): actual external seed, pair and
chunk workers retain their independent buffers, shared file owners and Source
cleanup order. The retained multinode prepare_inputs path now uses concrete
FASTA indexing, indexed reads and Source partition/load behavior. Source text
streams preserve physical offsets, numeric extraction and CRT output semantics.

Frozen corrective package
8f676da4b2d853f7fc4636d6a6a72bdd99507ae532f76d83ea70d0977d8f96e2
passes isolated default, optional dna/mcl/zstd and Windows GNU all-target checks,
plus 142 focused external, multinode, FASTA, sequence, text and scheduler tests.
The earlier package's one failing fixture is retained: Source may emit an empty
worker batch header. The corrected fixture checks every header, exact pair
payload and complete file length; production code did not change between those
two frozen packages. Inventory audit passes 199/199; this is not parity proof.

Hosted run 42478857080 (81214dc) finished with six of seven jobs passing,
including both macOS hosts. Its eight Windows text fixture failures are fixed
in local commit 894bfbb; hosted verification of that commit remains pending.
No new Source comparisons were run for this batch. Remaining code includes
parallel oversized radix spill, generic parallel chunk loading/sorting, the
concrete retained multinode search/GVC backend and Source seed-log number
formatting. Source's inactive multinode CLI remains inactive. The full
translation/parity goal stays open.


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
