# Focused parity and performance-scouting plan

Status: completed 2026-10-06. The final rerun passed 164/164 cases. This file
is retained as the executed protocol; results and artifact paths are in
`VALIDATION_RESULTS.md`. Work outside this plan remains in `GAPS.md`.

This plan covers search sensitivity presets, clustering, scoring matrices, and
DAA/view interoperability. Its purpose is to find correctness gaps and identify
the slowest or most memory-hungry modes for later profiling. It is deliberately
not a publication-quality benchmark campaign.

## Progress board

Keep this section as the high-level tracker. Detailed cases are checklists in
the sections below; observations and raw-result paths belong in
`VALIDATION_RESULTS.md`.

| Area | Compact parity | Natural-data parity | 1-thread scout | Multithread spot check | Worst offender recorded |
|---|:---:|:---:|:---:|:---:|:---:|
| Search sensitivities | [x] | [x] | [x] | [x] | [x] |
| Clustering | [x] | [x] | [x] | [x] | [x] |
| Scoring and gaps | [x] | [x] | [x] | [x] | [x] |
| DAA production/view | [x] | [x] | [x] | N/A | [x] |

## Suggested execution order

Work through these as four small, independently resumable batches. Check off
the detailed rows below and put measurements/diffs in `VALIDATION_RESULTS.md`.

1. **SENS** — all sensitivity presets on the quick fixtures; then only the
   unresolved or worst-looking cases on natural data.
2. **SCORE** — matrix defaults, one accepted alternate gap pair, and one
   rejected pair per matrix.
3. **CLUSTER** — the three clustering commands on quick data, followed by the
   family-rich natural fixture. Compare memberships, not just counts.
4. **DAA** — build one blastp and one blastx archive with each producer, then
   run the producer/viewer compatibility square and all view formats.

A batch is complete when parity status, one elapsed-time/RSS observation, raw
artifact paths, and any follow-up offender are recorded. A parity failure does
not block recording the rest of that batch, but it does block interpreting its
timing as an optimization result.

### Run budget

- Use the correctness run itself to collect one indicative wall time and peak
  RSS; do not add dedicated benchmark repetitions.
- Run the broad matrix with one thread. Use four threads only for the small
  spot checks explicitly named below.
- Repeat a pair only when shared-host load makes the result implausible. Note
  both observations rather than selecting the more favorable one.
- Stop performance analysis when outputs are not comparable. Record the
  correctness gap first; optimizing unequal work is misleading.
- Do not profile or optimize during this pass. At the end of each area, record
  at most the two worst time ratios and two worst RSS ratios for later work.

### Deferred-optimization note format

Add one entry per offender to `VALIDATION_RESULTS.md`:

```text
- [ ] AREA / CASE: parity status; Rust/C++ time ratio; Rust/C++ RSS ratio;
      likely phase if known; raw artifact directory; suggested next probe.
```

## Outcome and limits

For every case, compare the native Rust command with the bundled, pinned
upstream C++ DIAMOND revision. Record correctness first, followed by one
indicative wall-time and peak-RSS observation. Do not optimize while running
this matrix; add a note and a follow-up item when a mode is a clear offender.

The broad matrix runs at one thread. Run a four-thread spot check only for the
representative cases identified below. One run per implementation is enough
for scouting. If host load makes a result implausible, repeat that pair once
with reversed execution order and retain both observations. Do not present
these measurements as stable benchmarks.

## Common setup

- [x] Build C++ and Rust in release mode for the host CPU and record the exact
  Git commit, compiler versions, and command lines.
- [x] Use the same input FASTA for both `makedb` invocations, but let each
  implementation consume the database it produced.
- [x] Keep thread count, temporary directory filesystem, output fields, and all
  explicit options identical.
- [x] Capture exit status, stdout/stderr, output SHA-256, elapsed wall time, and
  maximum RSS using `/usr/bin/time -v` or the existing comparison driver.
- [x] Save raw commands, logs, outputs, and diffs below one dated results
  directory. Do not keep only the summary.
- [x] Use byte comparison where the format is deterministic. If byte comparison
  is inappropriate, use the format-specific semantic comparison described
  below and record the exact reason.
- [x] Mark each case `PASS`, `FAIL`, `EXPECTED DIFFERENCE`, or `BLOCKED`.

### Data sets

Use complete natural records; never lengthen a workload by duplicating them.

1. **Quick protein set:** the bundled 389-record corpus for fast iteration.
2. **Real protein set:** a disjoint reference/query split produced by
   `scripts/prepare_real_benchmark.sh`. Prefer the documented AMRFinderPlus
   snapshot or another versioned public protein corpus. Aim for about 3--5
   seconds in the slowest search preset, but accept shorter runs here because
   timing is diagnostic only.
3. **Clustering set:** a natural protein collection containing genuine
   homologous families, near duplicates, fragments, and singletons. It must not
   be the same sequences copied repeatedly. A few thousand records is enough.
4. **Translated set:** one natural genomic region or small assembly and the
   independent real protein reference used for the existing blastx checks.

Create a small manifest in the results directory containing source URL/version,
record and residue counts, and SHA-256 for every input. Record enough provenance
to recreate any untracked large fixture.

Use this stable tabular schema for search comparisons unless a test explicitly
requires another format:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

## 1. Search sensitivity presets

Run every preset once on the quick protein set to catch parsing and dispatch
errors. Run the missing presets below on the real set at one thread. The true
default cases must contain **no sensitivity flag**.

### blastp

- [x] no flag (default; control)
- [x] `--faster` (existing control)
- [x] `--fast`
- [x] `--mid-sensitive`
- [x] `--sensitive`
- [x] `--more-sensitive` (existing control)
- [x] `--very-sensitive`
- [x] `--ultra-sensitive` (existing control)

### blastx

- [x] no flag (true default; fix/avoid the historical `blastx-default` harness
  alias that silently adds `--sensitive`)
- [x] `--faster`
- [x] `--fast`
- [x] `--mid-sensitive`
- [x] `--sensitive` (existing control)
- [x] `--more-sensitive`
- [x] `--very-sensitive` (existing control)
- [x] `--ultra-sensitive` (existing control)

For each result, require byte-identical tabular output. On failure, retain the
first differing query plus both complete query result groups; a matching row
count is not parity. Also record raw-hit/output-row counts where available, as
these help distinguish seed-search differences from ranking or formatting
differences.

Four-thread spot checks:

- [x] blastp true default and the slowest Rust/C++ ratio found above.
- [x] blastx true default and the slowest Rust/C++ ratio found above.
- [x] For those four cases only, repeat Rust with its default 16 GB soft limit
  and `--memory-limit 0G`; upstream remains in its normal disk-backed mode.

Add a follow-up optimization note when Rust is more than about 20% slower or
uses more than about 25% extra RSS. These are triage thresholds, not release
criteria.

## 2. Clustering

Clustering was high risk at the start of this plan because the native
dispatcher sent `cluster`, `linclust`, and `deepclust` into one Rust
implementation. The completed translation now uses distinct workflows; the
historical warning below explains why membership, rather than counts, was the
acceptance test.

Run on both the bundled quick set and the natural clustering set:

- [x] `cluster` with defaults.
- [x] `linclust` with defaults.
- [x] `deepclust` with defaults.
- [x] Each command with `--approx-id 90 --member-cover 80`.
- [x] One lower-identity case, selected after confirming that both CLIs accept
  the same explicit options.

Compare clustering semantically as well as textually:

1. Validate that every input ID occurs exactly once as a member.
2. Parse each line as `(centroid, member)` and sort the pairs before comparison
   so harmless line ordering is ignored.
3. Require the same centroid/member pairs for strict parity. Do not hide a
   changed centroid or cluster membership by comparing only cluster sizes.
4. Record cluster count, singleton count, largest-cluster size, and a hash of
   the sorted pairs to make failures easy to characterize.

Run one-thread scouting measurements for all cases. Run four threads only for
the three default commands on the natural set. If parity fails substantially,
stop performance interpretation for that case—the algorithms are not yet doing
equivalent work.

## 3. Scoring matrices and gap penalties

BLOSUM62 is the control and BLOSUM45 already has end-to-end coverage. Run the
following with blastp default sensitivity on both the quick and real protein
sets:

- [x] BLOSUM62 with its default gap penalties (control).
- [x] BLOSUM45 with a known accepted gap pair (control).
- [x] BLOSUM50.
- [x] BLOSUM80.
- [x] BLOSUM90.
- [x] PAM30.
- [x] PAM70.
- [x] PAM250.

For each non-control matrix:

- [x] First omit `--gapopen/--gapextend` and confirm Rust selects the same
  defaults as upstream.
- [x] Obtain the upstream-accepted gap pairs from the bundled source/runtime
  validation rather than guessing them.
- [x] Test one additional accepted non-default gap pair.
- [x] Test one rejected pair and compare nonzero exit status and the essential
  diagnostic, without requiring byte-identical version banners or paths.

Require byte-identical output for successful searches, including score,
bit-score, and e-value fields. When a mismatch occurs, rerun the first differing
query with `score`, `positive`, `ppos`, `qseq`, and `sseq` added so matrix lookup,
traceback, and statistical-parameter errors can be separated.

Use one thread for the full matrix. Four-thread timing/RSS is needed only for
the slowest observed matrix and one PAM matrix. Record an optimization note if
a matrix-specific path is materially worse than BLOSUM62 after accounting for
different hit counts.

## 4. DAA production and `view`

Create one nontrivial blastp DAA and one blastx DAA containing multiple queries,
multiple hits, gapped alignments, and at least one query without a reported hit.
If native Rust search cannot currently produce DAA output, record the producer
rows as `BLOCKED` and open that wiring gap; do not silently use `--legacy` as the
Rust producer. The viewer half can still proceed with a C++-produced DAA.
Exercise a compatibility square rather than testing only files produced by the
same implementation:

| Producer | Viewer | Required |
|---|---|:---:|
| C++ | C++ | baseline |
| C++ | Rust | yes |
| Rust | C++ | yes |
| Rust | Rust | yes |

For both DAA inputs/viewers, test:

- [x] DAA (`100` / `daa`) round trip.
- [x] Tabular (`6` / `tab`) with the stable 12 fields above.
- [x] Flat JSON (`104` / `json-flat`).
- [x] PAF (`103` / `paf`).
- [x] SAM (`101` / `sam`).
- [x] Pairwise (`0`).
- [x] XML (`5` / `xml`).
- [x] Null output.
- [x] Edge output.

Comparison rules:

- Tabular, PAF, pairwise, and edge: try exact bytes first. If ordering alone
  differs, compare sorted complete records and retain the ordering diff.
- JSON produced directly by search: parse and compare values structurally;
  also retain the textual diff. For DAA `view`, compare exact bytes: pinned
  upstream DIAMOND 2.1.24 omits commas between per-query buffers, and the
  verbatim-translation contract requires preserving that malformed output.
- SAM: compare records exactly after removing only demonstrably volatile
  program/version header fields. Do not reorder alignments.
- XML: parse and compare the document structure after removing only
  demonstrably volatile version/invocation fields.
- Null: require a zero-byte output and successful status.
- DAA: compare raw bytes and header metadata first. Regardless of raw-byte
  equality, require both viewers to render equivalent tabular results from both
  producers. A DAA that only its own producer can read is a failure.

Also exercise view-time filters on one C++-produced and one Rust-produced DAA:

- [x] default `-k` behavior.
- [x] `--max-target-seqs 1`.
- [x] one `--top` value.
- [x] `--forwardonly` for blastx.
- [x] unaligned reporting/rejection where the workflow contract permits it.

Collect one-thread wall time and RSS once for each format on the same DAA. Run
four threads only if the view command actually supports threaded execution and
does different work. Flag the two worst Rust/C++ time ratios and two worst RSS
ratios for later profiling; do not tune them during this pass.

## Result notes template

The executable harness is `scripts/run_validation_plan.sh`. A full scouting
run with prepared natural fixtures is:

```bash
PROTEIN_REFERENCE_FASTA=/path/to/reference.faa \
PROTEIN_QUERY_FASTA=/path/to/query.faa \
CLUSTER_FASTA=/path/to/family-rich-proteins.faa \
NUCLEOTIDE_QUERY_FASTA=/path/to/natural-region.fna \
PROTEIN_SOURCE_URL=https://example/source \
PROTEIN_SOURCE_VERSION=version \
DNA_SOURCE_URL=https://example/genome \
RESULT_DIR=.tmp/validation-run \
  scripts/run_validation_plan.sh
```

Use `VALIDATION_AREAS="sensitivity scoring"` to select sections. The harness
returns nonzero when it finds a parity failure, which is expected while gaps
remain; inspect `summary.tsv`, `notes.md`, and the retained `artifacts/`
directory. `scripts/validation_compare.py` contains the cluster, JSON, SAM,
XML, sorted-record, and exact-byte comparators.

Append one row per case to a TSV or Markdown table using these fields:

```text
date area dataset command threads cpp_seconds rust_seconds speed_ratio_cpp_over_rust cpp_rss_kib rust_rss_kib rss_ratio_rust_over_cpp parity status artifact_dir notes
```

For each area, finish with a short note containing:

- the correctness failures, ordered by severity;
- the slowest two Rust/C++ ratios;
- the largest two Rust/C++ RSS ratios;
- whether extra work appears to be algorithmic, I/O, allocation/RSS, or output
  formatting;
- the exact saved artifact directory needed to reproduce it.

## Completion checklist

- [x] Every sensitivity preset has a parity result and true default was not
  confused with `--sensitive`.
- [x] All three clustering commands were compared independently and cluster
  membership—not only counts—was checked.
- [x] Every supported BLOSUM/PAM matrix has default-gap parity plus one explicit
  accepted gap-pair comparison.
- [x] C++ and Rust DAA files are readable by both viewers, and every native view
  format has a recorded parity result.
- [x] Indicative time/RSS observations and raw artifacts are retained.
- [x] Worst offenders are listed as later profiling/optimization tasks without
  treating single-run timings as stable benchmark claims.
