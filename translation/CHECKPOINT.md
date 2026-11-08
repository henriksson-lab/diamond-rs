# Translation checkpoint — 2026-11-08

This is an in-progress checkpoint committed at the user’s request. The full translation/parity goal remains open.

Production dispatch uses Rust; the C++ conformance oracle is available only behind the explicit, default-off `ffi` feature. Native command integration, raw sequence-title/taxonomy handling, Source I/O behavior, and macOS CI configuration are included in this checkpoint.

The current atomic migration of OS arguments, paths, byte-valued options, and their CLI/library callers has **not been compiled or tested as a combined snapshot**. Individual edited batches passed formatting and whitespace checks. Earlier validation recorded in the completion/progress documents applies to its named snapshots, not automatically to this checkpoint. macOS CI has not yet run on a hosted runner.

Next work: finish central search caller and test-fixture integration, compile the combined batch, then verify remaining parity and platform gaps. The external two-worker scheduling case remains unresolved. Other Source undefined-behavior and unavailable-platform qualifications remain documented in GAPS.md and TRANSLATION_COMPLETION.md.

Build outputs, temporary validation workspaces, runtime logs, and the adjacent immutable C++ source checkout are excluded from this commit.
