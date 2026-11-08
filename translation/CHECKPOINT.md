# Translation checkpoint — 2026-11-08

The full native translation/parity goal is active and unachieved. Production
command dispatch is Rust; the C++ oracle requires an explicit default-off `ffi`
feature and explicit selection. The adjacent Source tree stays immutable.

Latest committed code: `d4fac90`. Its frozen backend/storage batch passes three
compile profiles and 112 focused tests. The subsequent formatter-only package
passes default compilation and 51 scoped tests. These are named snapshot checks,
not proof of final parity.

The immutable `d4fac90` release matches three natural 9,704-record, two-round
OID-output controls at one/four threads against pinned Source: exact mapping,
raw streams, alignment row multiset and stage counters. All cases use one
cross-process worker. Full artifact comparison found chunk cleanup,
representative-worker and literal-path gaps; their fixes and bounded parallel
blocked/accession output are now being translated in the worktree. The combined frozen package now passes default/optional/Windows compilation
and 101 focused tests; its new release parity remains pending. Source comparison proofs are
retained under `.tmp/native-post-translation-20261107/external-parallel-controls/`.

Latest hosted run 42478857080 still tests `81214dc`: six of seven jobs passed,
including both macOS hosts. Its eight Windows text expectations are fixed in
local commit `894bfbb`; hosted verification awaits a push. Final-snapshot
integration/profile/platform/RSS checks, broader upstream parity and Source's
shared two-cross-process-worker scheduling hazard remain open.

Build products, temporary captures, runtime logs and the immutable C++ checkout
are excluded from commits. The entries below retain earlier checkpoint scopes;
later evidence supersedes their historical pause or incomplete-batch status.

CI repair follow-through (2026-11-08): addressed failures from hosted runs
42464375688 and 42466492609. Internal cluster/reassign searches now validate the
original workflow arguments before applying their trusted edge-format override.
Independent zstd frames reset repeat offsets and entropy state. Fixtures now
respect enabled compression features, generated ABI declarations, Windows reader
lifetimes, and byte-exact checkout contents. Spill calibration requires two
consecutive genuine delayed spills at one limit within a bounded search, while
search failures and byte mismatches remain fatal.

Local verification passed: isolated default/optional all-target checks, Windows
GNU all-target cross-check, default library and combined integration coverage,
optional library and all integration targets (including all 20 frozen regression
hashes), both profiles of the corrected compression fixture, and an all-features
package build without the C++ tree. The final reader-only fixture batch passed
120 default I/O tests, 119 optional I/O tests, and Windows cross-compilation.
All 18 raw-OS CLI controls matched preserved C++ results exactly without rerunning
C++; bounded calibration and error-propagation shell checks passed. Fresh hosted
macOS/Windows runtime results remain pending a push and CI rerun.

Follow-up for run 42472717649: both hosted macOS spill checks passed. Corrected
the numerical fixture's platform scope: all 54 GNU-profile vectors remain exact,
while three libm/NaN-sensitive cases also have independent mathematical checks
that run on every platform. Console-text expectations now use Windows CRT line
endings; archive/data bytes remain exact. A Windows-only guard prevents archive
subprocess tests from inheriting the parallel library test's open CRT handles.

At the user's request, CI now has seven jobs instead of fourteen: one packaged
native matrix for Linux, both macOS architectures, and Windows, plus musl and two
32-bit checks. Default tests, all-feature builds without C++, and spill checks
remain on native hosts; the full optional-feature runtime suite runs on Linux.
Local isolated default/optional/Windows compilation, protocol, numerical,
archive-library, stdout/help, and deprecated-command checks passed. Workflow YAML
and shell syntax checks passed. Fresh hosted verification of this batch remains
pending; the broader translation goal stays paused.

The user resumed translation work while monitoring the pushed CI run on
2026-11-08. The ten-minute check of run 42477158649 found Windows cutoff tests
blocked by inherited CRT handles. Commit 6e18720 serializes the Windows test
harness and removes the archive-only guard; internal parallel application cases
remain enabled. Other hosted jobs were still running at the last check.

The resumed translation batch adds raw unsigned-prefix partition parsing,
retained multinode library byte/path and text-stream semantics, a pre-read OID
boundary for external volumes, and Source TaskSet alignment dispatch with
per-worker codec buffers. The concrete command shares its scoring owner across
those tasks. Source-equivalent locking closes a task-completion lost-wakeup
window; the shared file-array owner allocates no worker codec state.

Corrective package 5c816812dceb517e2933f01288a0379034bc0325b7c1d13050a56c29edb30fe2
passed isolated default/optional/Windows all-target checks and 86 focused tests
(53 external, 21 multinode, 11 DMND, one successor-chain regression). The original
failed snapshot remains preserved. These checks do not establish fresh C++
parity for the new parallel branch. Nested DP pool reuse, parallel volume
building, retained multinode streaming ownership, and the Source-shared
two-worker scheduling hazard remain open; the translation goal is not achieved.

Active-goal continuation: nested external DP now uses priority-zero tasks on the
actual parent pool, including one-worker execution. Arc input/matrix/target
owners survive queued and reverse-pass work; null-pool calls dispatch directly.
Ordinary query and anchored callers explicitly share their configured pool;
completed DP ranges merge under the Source completion mutex. External volume
building uses independent readers and worker codec buffers with atomic claims.
The retained multinode library now retains real sequence/text owners and Source
formatted-extraction state. Its indexed prepare_inputs partition adapter remains
an explicit code gap; no inactive CLI route has been invented.

The shared text stream utility moved to util/io/source_text.rs and also serves
connected external manifests, preserving numeric prefixes, optional-field
fail-state values, blank-row errors and raw paths. Hosted Windows run
42478857080 exposed eight LF fixture assumptions after the CRT translation;
expected text bytes are now projected without normalizing actual/binary data.
The same pushed run's full macOS ARM64 job passed.

Frozen package b46aac0cd3a62ba07d37ae37ae828c6b06f43cbdb24094805c06662982393415
passes default/optional/Windows all-target checks and 204 focused tests with one
SIMD case ignored. The first batch's null-pool counter fixture failure is retained
separately. Fresh C++ parallel/profile parity and remaining completion gates stay
open. No original Source tree or algorithm object was rebuilt or modified.


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
status1 and no outer object unwind; ordinary embedded calls return normally.
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
