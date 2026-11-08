# Remaining work

This file lists current work only. Historical failures and completed fixes are
recorded in `VALIDATION_RESULTS.md` and `FFI_REDUCTION_PLAN.md`.

Current coverage (2026-11-07): baseline commands and enabled optional adapters,
including external clustering and streamed BLAST references, have concrete Rust
implementations. C++ is only an explicit default-off conformance adapter.
Current selected parity results are recorded in `TRANSLATION_COMPLETION.md` and
`translation/progress.md`; the full completion goal remains active. Final
integration/profile/package/platform and scale verification is still required.

The historical 2026-10-06 comparison matrix passed all 164 tested cases against
upstream DIAMOND 2.1.24: 40 search-sensitivity cases, 50 scoring cases, 18
clustering cases, and 56 DAA producer/viewer cases. No failure from that matrix
is known to remain in the current code.

That result is bounded evidence, not a claim that the translation is bug-free
or that every DIAMOND mode has been translated. The active work is:

## Native coverage and fidelity

- [x] Compare `makeidx` output and persistent-index search behavior with
  upstream on natural data. Index bytes and native consumption are exact;
  current speed/RSS measurements are tracked in the README.
- [x] Add a successful nontrivial `merge-daa` comparison. Two disjoint natural
  1,500-query archives now exercise reference deduplication and query appending;
  the 509 MiB merged archive and cross-view TSV are byte-identical to upstream.
- [x] Compare `getseq` output and measure `getseq`/`dbinfo` time and RSS. Full
  and selected-record output parity plus utility timings are recorded in
  `VALIDATION_RESULTS.md`. Standard upstream rejects the apparent hidden
  `--reverse`/`--hard-masked` controls; native Rust now rejects those and the
  other non-workflow options with matching diagnostics rather than silently
  ignoring them.
- [x] Extend the ordered DAA-view worker boundary beyond tabular/flat JSON to
  DAA, PAF, SAM, pairwise, XML, null, and edge output. The faithful streaming
  translation closes the 509 MiB null-output deficit: Rust now takes 0.26/0.21
  seconds at one/four threads versus 0.54/0.46 upstream, with lower RSS.
- [x] Implement upstream's one-billion-letter `makedb` block lifecycle. Rust
  now masks, writes, hashes, and releases each owned block before reading the
  next while retaining only the position trailer. A real natural >1B boundary
  RSS run remains an evidence gate, not an implementation gap.
- [x] Wire DMND taxonomy metadata into native `blastp` and `blastx`.
  Taxonomy-enabled `makedb`, selective metadata loading, and every supported
  taxonomy output field now match upstream. The production-threshold natural
  spill evidence and NCBI BLAST-database UB policy remain separate gates in
  `GAPS.md`.
- [x] Validate the enumerated native search options, output fields, multi-file
  inputs, and compressed-input combinations listed in `GAPS.md`. Taxonomy-aware
  NCBI database routing and the unavailable historical mouse corpus remain
  separate compatibility/verification gates below.
- [ ] Reconstruct the historical 500-query mouse fixture and rerun
  `tests/divergent_blastp_rows.rs`. The old checkpoint had 137 scoring-field
  differences; later natural search/scoring matrices pass, but this exact
  corpus has not been rerun and must not be described as fixed yet.

## Known user-visible limitations

- [x] Retain the repository-only C++ `ffi`/`--legacy` adapter as a diagnostic.
  It remains disabled by default and absent from routine benchmarks. Packaged
  all-features builds detect the intentionally absent upstream source tree and
  leave the adapter explicitly unavailable instead of failing to compile.
- [ ] Decide whether to preserve or repair upstream DIAMOND 2.1.24's malformed
  multi-query `json-flat` output. Rust currently reproduces those bytes for
  exact parity, so this is a known inherited behavior rather than a Rust/C++
  difference.
- [ ] Complete native macOS and Windows runtime validation, including RSS units
  and memory-limit spill decisions. Cross-compilation is not runtime evidence.
- [ ] Verify the connected external `linclust --parallel-tmpdir --aln-out`
  workflow at the final snapshot. Its actual native seed/pair/chunk/job and
  alignment-artifact pipeline is implemented; broader Source parity and
  multi-worker/scale evidence remain required.
- [x] Match upstream's wrapper-only target-linear Hamming/stage-1 pruning in
  the 4,000-protein split-input `linclust`. The missing `lin_stage1_target`
  handoff now dispatches the exact
  target-linear kernel instead of the full seed-group cross product. All three
  repeated 4,000-protein forced-split outputs are byte-identical; see
  `.tmp/cluster-lowmem-agent/README.md`.
- [ ] Retain the Source compatibility qualification for the historical
  10,000-protein forced-split `linclust` comparison. The following diagnosis
  records the 2026-10-12 evidence, not a missing external/native adapter. The
  first super-block already produces 717 clusters in Rust versus 719 upstream,
  so this does not originate in the now-exact later-block target-linear
  centroid search. The divergence is caused by upstream reading its
  uninitialized `current_ref_block` before deciding whether to length-sort the
  first linear query block: the 4,000-record run sorts, while the 10,000-record
  run does not. Do not emulate that process-layout accident with an input-size
  heuristic; either establish an upstream-defined rule/fix or document the UB
  compatibility boundary explicitly.
- [x] Forward cluster `--matrix`, `--gapopen`, and `--gapextend` through the
  translated per-round search. Default, PAM250, and alternate-gap natural
  outputs now match upstream exactly.
- [x] Forward the remaining standard cluster search controls. `--evalue`,
  `--comp-based-stats`, hard/soft/motif masking, and `--ext` now propagate
  through every round and pass exact natural controls; the repeated per-round
  TANTAN construction remains fixed.

## Packaging

- [x] Track `Cargo.lock` in version control; this repository ships a CLI
  application and now has reproducible dependency resolution.
- [x] Keep the `ffi` feature buildable but explicitly unavailable to package
  consumers. The build script only links the adapter when the adjacent
  upstream tree exists; packaged `--all-features` builds do not require or
  silently substitute for that intentionally excluded source tree.

Performance measurements in the README are checkpoints, not promises. Profile
again only after establishing parity for a newly covered mode; do not optimize
unequal work.
