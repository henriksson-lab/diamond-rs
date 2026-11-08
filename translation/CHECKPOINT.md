# Translation checkpoint — 2026-11-08

This is an in-progress checkpoint committed at the user’s request. The full translation/parity goal remains open.

Production dispatch uses Rust; the C++ conformance oracle is available only behind the explicit, default-off `ffi` feature. Native command integration, raw sequence-title/taxonomy handling, Source I/O behavior, and macOS CI configuration are included in this checkpoint.

The current atomic migration of OS arguments, paths, byte-valued options, and their CLI/library callers has **not been compiled or tested as a combined snapshot**. Individual edited batches passed formatting and whitespace checks. Earlier validation recorded in the completion/progress documents applies to its named snapshots, not automatically to this checkpoint. macOS CI has not yet run on a hosted runner.

Next work: finish central search caller and test-fixture integration, compile the combined batch, then verify remaining parity and platform gaps. The external two-worker scheduling case remains unresolved. Other Source undefined-behavior and unavailable-platform qualifications remain documented in GAPS.md and TRANSLATION_COMPLETION.md.

Build outputs, temporary validation workspaces, runtime logs, and the adjacent immutable C++ source checkout are excluded from this commit.

The user paused the broader translation goal on 2026-11-08 and prioritized the
CI/CLI fixes. Current verification processes and all earlier Source captures
are retained; the pause is not a completion claim.

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
