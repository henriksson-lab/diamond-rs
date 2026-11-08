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
