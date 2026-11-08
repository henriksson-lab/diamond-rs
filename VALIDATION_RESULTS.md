# Validation scouting results

## Current status

The latest complete rerun is the **Final combined exact-parity rerun —
2026-10-06** near the end of this file. It passed 164/164 cases with zero
failures or blocked cases: search sensitivity 40/40, scoring 50/50, clustering
18/18, and DAA production/view 56/56. Earlier failure tables are retained as
diagnostic history and must not be read as the current state.

This result applies only to the explicitly tested matrix. Open command,
option, input, portability, and packaging coverage is tracked in `GAPS.md` and
`TODO.md`; no untested surface is claimed to have parity.

Date: 2026-10-04  
Rust commit: `be03674` plus the current validation-harness worktree changes  
Upstream: bundled DIAMOND 2.1.24 build

These are single-run scouting observations, not publication-quality
benchmarks. They are intended to preserve parity failures and identify the
first areas to profile. Raw commands, stderr, outputs, diffs, timings, and
manifests are retained under `.tmp/validation-real-*`.

## Fixtures

No records were duplicated or synthesized.

| Role | Natural source | Records / span | SHA-256 |
|---|---|---:|---|
| Protein source | NCBI AMRFinderPlus database 2026-08-07.1 `AMRProt.fa` | source snapshot | `41d5ebf4f807c9590f27de7dd23c9f037afb4d8efd613790c54847fb1fc2896b` |
| Search reference | evenly sampled AMR proteins | 1,000 | `45bf7cbdc841a7418565a26279696b7f812f2c841e1253acde19fe50717235c7` |
| Search query | disjoint, evenly sampled AMR proteins | 300 | `3e6b70ccb70c8f0ebee02718296f8899383a91d01f4954f263a8e2e104e58029` |
| Clustering | evenly sampled AMR proteins | 500 | `60caa8e27deaee62141bd7ccb79eaa717bc086243178451d09f6873c5672c04f` |
| Translated query | contiguous bases 1--250,000 of *E. coli* K-12 MG1655, GCF_000005845.2 | 250,000 nt | `6e1dfc494ad02df445e81c61236bfe4b8903805b877135f64b4aebb54ea7554d` |

The compact matrix also used the bundled 389-record protein corpus. Input
paths, counts, residue totals, host details, compiler versions, and hashes are
in each run's `manifest.tsv` and `sources.tsv`.

## Initial scouting result (superseded)

| Area | Pass | Fail | Blocked | Primary follow-up |
|---|---:|---:|---:|---|
| Sensitivity | 39 | 1 | 0 | blastp `--mid-sensitive` misses two reciprocal compact-fixture hits |
| Scoring | 41 | 7 | 0 | PAM250 seed/search behavior; BLOSUM50 traceback/statistics on natural data |
| Clustering | 1 | 17 | 0 | replace the shared greedy implementation with faithful command-specific algorithms |
| DAA/view | 22 | 3 | 2 | native DAA producer, `--unal` behavior, and valid multi-query JSON |
| **Total** | **103** | **28** | **2** | |

## Search sensitivity

All eight blastp and all eight blastx presets were compared at one thread on
both compact and natural data. True defaults used no sensitivity flag. Default
and ultra-sensitive also received four-thread scouting checks, including Rust
default memory behavior and `--memory-limit 0G`.

- Every natural-data sensitivity case is byte-identical.
- Compact blastp `--mid-sensitive` emits 1,336 Rust rows versus 1,338 upstream
  rows. Rust misses reciprocal `d3t6ka_`/`d3crna1` alignments.
- Natural blastp ranges from 0.995x to 2.229x C++/Rust speed. The closest case
  is `--very-sensitive`: C++ 4.06 s, Rust 4.08 s.
- Natural blastx is 1.410x--4.237x faster in this scouting pass. Its largest
  observed RSS excess is small: 1.042x at `--mid-sensitive`.
- The four forced-disk spot checks are byte-identical. They are retained in
  `.tmp/validation-real-sensitivity/summary.tsv`.

Historical classification: the remaining parity problem at this checkpoint
was sensitivity-specific search coverage, not formatting or retained-hit
storage. The final matrix closes it.

## Scoring matrices

Default gaps, one accepted alternate gap pair, and a rejected `1/1` pair were
tested for BLOSUM45/50/62/80/90 and PAM30/70/250. Rejected-pair status and
essential diagnostics agree.

- BLOSUM62, BLOSUM45, BLOSUM80, BLOSUM90, PAM30, and PAM70 are byte-identical
  on compact and natural data for all successful cases.
- BLOSUM50 passes the compact corpus but differs on natural data with the same
  row count. The first difference changes one 453-residue alignment from 286
  to 287 mismatches, with corresponding identity, score, and e-value changes.
- PAM250 is a major correctness failure. On the compact corpus upstream emits
  8 rows and Rust emits 709, indicating matrix-dependent seed/filter behavior
  rather than a small rendering difference. It also fails on natural data.
- Successful natural scoring cases use less Rust RSS and are faster in these
  one-off observations; correctness remains the gating concern.
- Four-thread spot checks pass: BLOSUM62 is 2.333x faster at 0.499x RSS and
  PAM70 is 2.467x faster at 0.513x RSS. These single observations are stored in
  `.tmp/validation-real-scoring-t4/summary.tsv`.

Classification: inspect matrix-dependent reduction/seed filtering for PAM250,
then BLOSUM50 traceback and statistical parameter application.

## Clustering

`cluster`, `linclust`, and `deepclust` were independently compared at defaults,
90% identity/80% member coverage, and 50% identity/60% member coverage. Default
commands were also run with four threads. Comparisons validate every FASTA ID
occurs exactly once and then compare sorted `(centroid, member)` pairs.

Only compact `linclust --approx-id 90 --member-cover 80` matches. All natural
cases differ. At natural defaults:

| Command | Upstream clusters | Rust clusters | Upstream/Rust speed | Rust/upstream RSS |
|---|---:|---:|---:|---:|
| `cluster` | 169 | 143 | 0.404x | 1.516x |
| `linclust` | 261 | 143 | 0.082x | 2.654x |
| `deepclust` | 118 | 143 | 0.768x | 1.191x |

The identical 143-cluster Rust result across all three commands confirmed that
they shared one algorithm at this checkpoint despite upstream's distinct
workflows.
Four threads do not materially accelerate Rust because the implementation is
effectively serial. Rebuilding the centroid seed index during incremental
clustering is the likely dominant performance and allocation problem.

Classification: this is an algorithm-fidelity problem first and an
algorithmic-allocation/performance problem second. Do not micro-optimize the
current shared implementation.

### Clustering translation follow-up — 2026-10-06

The table above is retained as historical scouting evidence for the former
shared greedy implementation. The native commands now use separate translated
cascaded workflows. On the bundled 389-record compact corpus, all six current
single-thread comparisons are byte-identical to upstream:

| Command | Case | Clusters | Output SHA-256 |
|---|---|---:|---|
| `cluster` | default | 385 | `43ca4fa73047dc65a3ffcfb76f9d92584fc1589a1b62193d71a7737ef43fa847` |
| `cluster` | `--approx-id 90 --member-cover 80` | 388 | `95e6a84f8084016de0dbc04cf7821fe7200b36df1cfed8861c21c44fedeed72c` |
| `linclust` | default | 389 | `cb8a631c410cca130f98268a841b5c7c14c74e2bbc15fdff4d4286c8d2f8b96d` |
| `linclust` | `--approx-id 90 --member-cover 80` | 389 | `cb8a631c410cca130f98268a841b5c7c14c74e2bbc15fdff4d4286c8d2f8b96d` |
| `deepclust` | default | 175 | `5b0090356cdc3f0b18c8ea8a26bb7b374f082fc5e36b49e5eb15848e166efa91` |
| `deepclust` | `--approx-id 90 --member-cover 80` | 388 | `95e6a84f8084016de0dbc04cf7821fe7200b36df1cfed8861c21c44fedeed72c` |

The last compact divergence was caused by active-block ordering in linear
self-search rounds. Upstream sorts by descending `(sequence length, prior
BlockId)` before seed enumeration; carrying that exact permutation through the
local-to-original OID map restored the direction-sensitive CBS alignment and
default `deepclust` parity. Natural-data parity and performance/RSS remain to
be rerun against the translated workflows; no result from the historical table
should be treated as a current benchmark.

## DAA and view

Upstream-produced blastp and blastx DAA files were read by both viewers. DAA,
tabular, PAF, SAM, pairwise, XML, null, edge, `-k 1`, `--top 10`, and blastx
`--forwardonly` outputs match under their format-specific comparators. Raw DAA
round trips are byte-identical.

- Blastp `json-flat` output is byte-identical but invalid JSON in both tools:
  separators are missing between some query groups. Structural parsing catches
  this; byte comparison alone would not.
- Upstream rejects `view --unal 1` for this workflow, while Rust accepts it.
  Blastp and blastx therefore fail status parity.
- At this checkpoint, native Rust blastp/blastx could not produce DAA output.
  The reverse producer/viewer half of the compatibility square was explicitly
  `BLOCKED`; the harness does not substitute `--legacy`.
- View operations are too short for meaningful performance claims. Pairwise,
  XML, and edge had the worst observed Rust time ratios and should be revisited
  only with a larger DAA after correctness work.

Historical classification: DAA production was missing functionality; JSON and
`--unal` were output/CLI semantics. The later follow-up below closes these
parity gaps. Viewer RSS was not an offender in this scout.

## Historical optimization order (completed or superseded)

1. Retranslate clustering as three faithful workflows; its current results are
   neither comparable nor performant.
2. Fix PAM250 seed/filter parity.
3. Fix the BLOSUM50 alignment discrepancy and compact `--mid-sensitive` misses.
4. Wire native DAA production and complete the producer/viewer square.
5. Correct or explicitly document upstream-compatible malformed JSON and
   reject unsupported `view --unal` consistently.
6. Only after parity, profile clustering and the slowest larger-DAA view
   formats with repeated measurements.

## Reproduction

The command form is documented in `VALIDATION_PLAN.md`. The actual outputs for
this pass are:

- `.tmp/validation-real-sensitivity/{summary.tsv,notes.md,artifacts/}`
- `.tmp/validation-real-scoring/{summary.tsv,notes.md,artifacts/}`
- `.tmp/validation-real-scoring-t4/{summary.tsv,notes.md,artifacts/}`
- `.tmp/validation-real-clustering/{summary.tsv,notes.md,artifacts/}`
- `.tmp/validation-real-view/{summary.tsv,notes.md,artifacts/}`

The runner deliberately exits nonzero when parity failures are present.

## Native-ABI reduction follow-up — 2026-10-06

The libc/platform cleanup was regression-tested against clean Rust checkpoint
`c954916` with four alternating runs per cell on natural AMRFinderPlus data.
The first portable-`memchr` build was rejected after profiling explained a
2.8–7.5% search slowdown: its AVX2 dispatch/search path replaced glibc's EVEX
`memchr`. Calling the same runtime kernel through the maintained `libc` crate
removed the handwritten declaration and restored performance.

The final rewrite is within 0.9% of the checkpoint in seven of eight search
cells and 2.2% faster in the remaining forced-disk four-thread cell. Genuine
delayed 64 MiB spill runs are within 0.3%; representative clustering is within
1.0%. RSS is effectively neutral except for a 3.1–6.4% reduction in forced
disk mode. All 64 final paired outputs are byte-identical. Detailed timing/RSS
tables, hashes, `perf` attribution, and retained artifacts are in
`FFI_REDUCTION_PLAN.md` and `.tmp/ffi-prepost/`.

The remaining Linux spill-decision boundary was tested with four alternating
runs per checkpoint/rewrite binary and limit on the natural 2,000-query
fixture. Both binaries spilled in all eight `96M` runs and neither spilled in
any of the eight `104M` runs; all delayed spills read back 20.6--20.7 MiB.
Closer `99M` and `100M` trials were deliberately retained as instability
evidence: decisions were mixed for both binaries as process RSS varied between
samples. All 32 outputs exited zero and had SHA-256
`d324e3efeaa48cff4885dfb4ad9c880a699b5f9a8c372b97f84629e2f25ccfc9`.
This closes the Linux half at the narrowest scheduler-robust bracket
established by the pass; native macOS and Windows spill-decision validation
remains open. Reproduction artifacts are under `.tmp/spill-boundary/`.

A current-binary performance refresh keeps this path visible in the wider
portfolio. Three paired `blastp --sensitive -p4` runs on 2,000 natural queries
against 4,000 natural targets give medians of 27.53 s/92,200 KiB for upstream's
normal disk path and 29.59 s/79,944 KiB for Rust at `--memory-limit 64M`
(speed 0.930x, RSS ratio 0.867x). All six 41,556-row outputs have SHA-256
`18c53b2d8afc218482cd9836e778f530964ab3ee39c3f27f130a3bbbacca63ad`.
Each Rust run explicitly migrated 1,594,528 retained hits and read 37.5 MiB of
compressed spill data. The 7.5% runtime gap is retained as a future bounded
optimization target; it does not reopen the already-closed spill correctness
work. Artifacts are under `.tmp/spill64-refresh/runs/`.

## Correctness follow-up — 2026-10-05

These checks were focused correctness validations rather than a new benchmark
pass; the historical scouting counts and timings above remain unchanged.

- PAM250 now follows upstream Tantan lambda-calculation failure semantics
  instead of using a BLOSUM fallback. Default compact output is byte-identical
  at 8 rows; natural default output and an accepted `15/3` gap variant are also
  byte-identical.
- Pinned upstream's DAA `ViewWriter` does not apply the flat-JSON query
  separator between per-query buffers. The Rust translation now preserves
  that malformed multi-query byte stream exactly instead of repairing it.
- Rust and upstream now both reject `view --unal 0` and `view --unal 1` with
  status 1 and the essential diagnostic `Option is not permitted for this
  workflow: unal`.
- Native Rust DAA production is implemented for blastp and blastx. The compact
  blastp archive contains 731 alignments and the blastx archive contains 5;
  both C++ and Rust viewers render byte-identical tabular output for archives
  from either producer. The blastp render hash is
  `ed82209f28b8b6327673e00ac4c82cf4da74c72fb250a03b792afd755a4f26d8`;
  blastx is
  `c9e71aaf60f902be5c84c0c55278fb3461dc85705fc0c1bba1d96b27bd8b3fc8`.
- The validation runner now creates both producer archives and runs both
  viewers/formats instead of recording Rust DAA production as blocked.

The earlier release-mode quick matrix is retained in
`.tmp/validation-daa-final/`. Its two JSON differences led to the upstream
source audit above. The later natural-fixture rerun in
`.tmp/ffi-final-view-byte-parity/` records 56 passes, 0 differences, 0
failures, and 0 blocked cases across both producers, both viewers, every
listed view format/filter, and `--unal 0/1` rejection. All four blastp
producer/viewer flat-JSON files have SHA-256
`acd49b19d4cece6bea7e59a5ee9274cb097616f58d2883d19b22bbe8ee3f3ffd`;
all four blastx files have
`0117a6f8e5d4463dda2c6b3aefa07971455e22b4a4045aec196b54e485693f9d`.

No DAA/view fidelity case remains open in the current matrix. The malformed
flat JSON is an upstream behavior, covered by a focused byte-exact regression
test and retained artifacts rather than classified as an intentional Rust
difference.

## Post-FFI broad correctness rerun — 2026-10-06

The release binary rebuilt from the current FFI-reduction worktree (SHA-256
`0aa98a5adcf5fff86a33af3d6a3a18fdd2b5676f8c2f0bfcef13a1bcbfd3df07`)
was compared with upstream DIAMOND 2.1.24 using the full validation runner and
the retained natural AMR/genomic fixtures.  The exact invocation and fixture
sizes/hashes are recorded in `FFI_REDUCTION_PLAN.md` and
`.tmp/ffi-final-parity-current/manifest.tsv`.

The original broad matrix produced 152 passes, 2 known differences, 10 failures, and no
blocked cases:

- all 40 search-sensitivity cases passed with byte-identical output;
- all 50 scoring cases passed (42 exact accepted outputs and 8 matching
  rejected-gap cases);
- all 6 compact clustering cases passed, but only 2 of 12 natural clustering
  cases passed;
- all DAA producers and 52 of 54 producer/viewer format combinations passed
  in that run. A source audit then established that pinned upstream's
  malformed blastp flat JSON is the literal behavior to translate, and the
  focused rerun passes all 56 DAA/view cases byte-for-byte or under the stated
  format-specific comparator.

The ten natural clustering failures are real membership/centroid differences,
not line ordering: the comparator resolves FASTA identifiers and sorts pairs
before hashing.  Therefore this run does not support an overall exact-parity
claim.  It does establish current-source parity for the complete tested search
and scoring matrices, compact clustering, natural default `deepclust` at one
and four threads. At that stage, after the DAA/view rerun, only the ten natural
clustering cases remained parity failures. Commands and full artifacts are
retained under `.tmp/ffi-final-parity-current/` and
`.tmp/ffi-final-view-byte-parity/`.

### Cluster Hamming-extension follow-up — 2026-10-06

Upstream enables its cluster-only Hamming extension when approximate identity
is at least 50%. Rust previously computed that option in the clustering
wrapper but did not propagate it into the ungapped stage, causing qualifying
partial HSPs to proceed through full Smith-Waterman extension. After threading
the exact identity, coverage, E-value, diagonal-filter, mutual-cover, and
linear-stage configuration, the natural single `faster_lin` round is
byte-identical: 923 alignments, 195 clusters, 500 pairs, 141 singletons, and
TSV SHA-256
`63defdc92c491b5e117e3dc9447237ff458ea568faffb9f43e545f06b1ab1e85`.

The complete clustering-area rerun under
`.tmp/ffi-cluster-hamming-parity/` has 13 passes, 5 failures, and no blocked
cases. All 6 compact cases pass; natural clustering improves from 2/12 to
7/12 passes. The remaining natural failures are `cluster` default at one and
four threads plus `cluster`, `linclust`, and `deepclust` at
`--approx-id 50 --member-cover 60`. This is measured progress, not a complete
clustering-parity claim.

The five failures were repeated after rebuilding with the corrected
double-precision coverage calculation; an identity-90 control still passes,
while the failures retain the same 42-member centroid swap. Default round 1
now matches upstream's aggregate 923 HSPs, 1,831 directed edges, and 195
clusters. The earliest measured count divergence is round 2: upstream/Rust
HSP counts are 30/28 for default, 31/30 for identity-50
`cluster`/`deepclust`, and 36/34 for identity-50 `linclust`. Further work must
compare exact round-2 active OIDs, seed hits, and HSP pairs.

### Round-2 ordering and structural endpoint follow-up — 2026-10-06

The remaining clustering differences were caused by translation gaps, not a
score or threshold tolerance. Upstream preserves active OID order after the
first linear round when `keep_target_id` is enabled; Rust had length-sorted
every round. The endpoint-only path also lacked upstream's structural
no-trace SIMD kernel. Restoring the exact ordering/predicates and translating
the 32-lane i8, 16-lane i16, and scalar i32 promotion cascade closes the
measured gap.

Two fresh complete clustering matrices pass 18/18 exact comparisons: all six
compact cases and all twelve natural `cluster`, `linclust`, and `deepclust`
cases. The retained result sets are
`.tmp/cluster-endpoint-structural-final/summary.tsv` and
`.tmp/cluster-endpoint-structural-repeat2/summary.tsv`. The structural SIMD
matrix took 30.91 s and 26.83 s in the two runs, versus 37.68 s for the
retained scalar-Rust baseline (1.22x and 1.40x faster). Maximum Rust RSS was
9,280 and 8,760 KiB versus 8,688 KiB for that baseline. Timing was collected
on a shared host and is diagnostic; exact output parity is the acceptance
result.

## Final combined exact-parity rerun — 2026-10-06

After integrating the structural endpoint kernels and exact semi-global
narrow score/trace behavior, a first broad run exposed one integration
regression: traced HSP range starts had been zeroed even though their
transcripts retained the correct traceback origin. Restoring
`sw.query_begin`/`sw.subject_begin` and the local zero-cell traceback boundary
made the saved 731-alignment quick blastp oracle byte-identical again.

The rebuilt binary then passed every fresh area matrix:

| Area | Result | Retained evidence |
|---|---:|---|
| Search sensitivity | 40/40 exact | `.tmp/ffi-final-parity-fixed-search-scoring/` |
| Scoring | 50/50 exact | `.tmp/ffi-final-parity-fixed-search-scoring/` |
| Clustering | 18/18 exact | `.tmp/ffi-final-parity-fixed-clustering/` |
| DAA production/view | 56/56 exact | `.tmp/ffi-final-parity-fixed-view/` |

The aggregate is 164/164 passes, with zero failures, blocked cases, or
expected-difference classifications. Both final DAA producer comparisons are
exact after decoding: blastp has SHA-256
`b6d56b102f9e396ea8ac27883456ec49608b01a4f0c1d42334b41a8f855a323e`
and blastx has
`63ba7f46a8927e5d6aadd6181b7a29d1186515bbc4d828b9c4b27b9f3a07183e`
for both upstream and Rust TSV output.

## CI and performance portfolio follow-up — 2026-10-11

The current worktree passes `cargo test --all-targets --no-fail-fast`: 1,707
library tests pass, 6 diagnostics are ignored, and every binary/integration
target passes. CI-specific fixes include AArch64 NEON transpose ordering,
Windows-safe test paths, Windows-safe locked `FileStack` I/O, and platform
correct unlink/path expectations. Native Windows and macOS CI must still be
observed after publication; local cross-checks cannot validate their runtime
filesystem and RSS behavior.

Performance work was deliberately split across three command families:

- Search: one histogram per shape is now reused across one-thread index
  chunks, matching upstream. A five-run natural 1 Mb `blastx
  --mid-sensitive` comparison is within 1.1% of upstream and uses 29.7% less
  RSS. The heavier forced-disk `blastp --sensitive` comparison remains open:
  three alternating four-thread pairs give 26.44 s / 97,832 KiB for C++ and
  32.60 s / 72,644 KiB for Rust, with exact SHA-256
  `18c53b2d8afc218482cd9836e778f530964ab3ee39c3f27f130a3bbbacca63ad`.
- Clustering: the post-change compact/natural matrix passes 18/18. Realistic
  rows are broadly at parity or faster and use substantially less RSS; the
  retained scout is `.tmp/perf-portfolio-post-tile-cluster/`.
- DAA/view: all checked formats remain exact. On the 23 MiB DAA producing a
  101 MiB pairwise report, five alternating runs give 2.73 s C++ versus 2.71
  s Rust while Rust uses 5,440 KiB versus 204,868 KiB. Generic per-character
  formatting was replaced by direct byte writes, matching the upstream text
  buffer behavior.

The sensitive-search profile places the remaining runtime difference in the
stage-1 worker. Restoring upstream's 1024-target Hamming tile and query-major
stage-2 handoff improves Rust while retaining exact bytes, but does not close
the gap. This result is intentionally documented as unfinished rather than
treated as parity because default in-memory mode is faster.

### Cross-command rebuild checkpoint

After rejecting a const-generic stage-1 dispatch experiment that regressed
the large sensitive-search workload, the restored release binary was checked
across the complete natural portfolio. The retained result set is
`.tmp/validation-cross-tool-natural-checkpoint/summary.tsv`: 164/164 exact
passes, with no failures or blocked cases. It covers all search sensitivities,
`cluster`/`linclust`/`deepclust`, all supported scoring matrices and selected
gap penalties, both DAA producers, and every supported `view` format.

The rotation initially found two candidates outside the already-known large
sensitive-search case, but startup-independent follow-up cleared both. On
7,000 references and 3,000 disjoint natural queries, PAM250 is exact and
faster in Rust across three pairs (2.17 s versus 2.30 s, speed ratio 1.060x,
RSS ratio 0.493x). Five exact repetitions of natural `deepclust --approx-id
50 --member-cover 60` give medians of 1.59 s Rust and 1.60 s C++, with RSS
6,708 versus 27,036 KiB. The smaller apparent deficits were startup/shared-host
effects and are not current optimization targets.

### Whole-worker AVX2 dispatch follow-up

Fresh paired profiles of the 716,941,027-raw-pair sensitive workload showed
that upstream executes `Search::ARCH_AVX2::stage1`, while Rust independently
selected an AVX-512 fingerprint sub-kernel before returning to AVX2 stage 2.
Production Hamming dispatch now selects the upstream AVX2 worker boundary as
a unit. The AVX-512 kernel remains covered by direct SIMD-vs-scalar tests.

Three alternating `blastp --sensitive --threads 4` pairs with Rust forced to
`--memory-limit 0G` are byte-identical at SHA-256
`18c53b2d8afc218482cd9836e778f530964ab3ee39c3f27f130a3bbbacca63ad`.
Medians improve from the previous 26.44/32.60 s C++/Rust comparison to
28.42/29.48 s under the current shared-host load (speed ratio 0.964x), with
RSS 91,476/68,960 KiB (ratio 0.754x). The first isolated Rust diagnostic was
29.20 s. Cross-tool acceptance after this change passes 18/18 clustering and
56/56 DAA/view comparisons in
`.tmp/validation-avx2-dispatch-cross-tool/`.

The next fidelity pass removed a redundant query-boundary scan inside stage
2. The streaming caller already emits exactly one query per callback, just as
upstream `search_query_offset` receives one query plus its target-index list.
Keeping the existing `SeedMatch` scratch layout while consuming that contract
directly is performance-neutral within shared-host variance: across the two
three-pair series, the combined medians are about 28.43 s C++ and 29.53 s
Rust. The more invasive compact-target scratch experiment was rejected after
an exact but slower 31.25 s full run. Post-change cross-tool validation again
passes 18/18 clustering and 56/56 DAA/view rows under
`.tmp/validation-single-query-stage2-cross-tool/`.

### Native utility-command rotation — 2026-10-11

The next pass deliberately moved away from the search kernel and exercised
successful natural-data workflows for `getseq`, `makeidx`, `merge-daa`, and
`dbinfo`. Artifacts are retained under `.tmp/validation-misc-commands/`.

- `getseq` all-record output is byte-identical on the 7,000-sequence database.
  A selected `--seq 1 3` comparison first exposed that Rust incorrectly read
  a single seqid while upstream reads space-separated, 1-based numeric
  ordinals. After translating that behavior, both selected outputs have
  SHA-256 `9ce34dc247907c81aeaf3ad2c555c16aee0e622ec275cb2e444eb4f01127d706`.
  The selected startup scout was 0.03 s/9,740 KiB upstream and 0.02 s/9,600
  KiB Rust.
- `merge-daa` merged the natural upstream and Rust archives into the same
  byte-for-byte 30 MiB DAA, SHA-256
  `ed76c6a1ad46871198bb09c1ec0d41fbf3b341b9bdf24cd76d853dd9d2463178`.
  One scout measured 0.18 s/72,312 KiB upstream and 0.21 s/8,000 KiB Rust.
  Both viewers then emitted the same 133,214-line tabular file. Because the
  inputs are results for the same query set, this is a repeated-query
  merge-correctness scout rather than a natural optimization fixture. Its
  tabular view nevertheless exposed a runtime target: 0.19 s upstream versus
  0.92 s Rust,
  while RSS was 115,348 versus 4,800 KiB. This needs a longer repeated workload
  before optimization acceptance.
- `dbinfo` reports the same database type, format/build, sequence count, and
  letter count. Its startup-scale scout was 0.01 s/9,740 KiB upstream versus
  0.00 s/2,560 KiB Rust.
- The initial `makeidx` comparison found different deterministic 2.1 MiB
  index bytes and an upstream cross-consumption difference of 64,208 versus
  64,266 alignments. The cause was a direct translation omission:
  `Block::soft_mask(Motif)` did not apply upstream's motif masking-table
  operation. After restoring it, both complete indexes have SHA-256
  `9eb381487e6af429f93d38061dc0a595f5df943108aa717f6105c51174c11597`.
  The post-fix construction scouts were 0.21 s/20,800 KiB upstream and 0.25
  s/22,480 KiB Rust. Native consumption was added and made exact in the
  follow-up below.

### Persistent-index consumption scout — 2026-10-11

Native search now loads and uses the byte-identical `.seed_idx` artifact for
`--target-indexed -c1`; it no longer routes that mode away from Rust. The
translation uses hashed seed encoding, one-pass target-set-filtered query
enumeration, the translated hash join, and the normal stage-1 Hamming and
downstream extension path. Per-shape query/reference seed counts exactly match
upstream (28,375/2,286,238 and 28,108/2,265,280).

The first scout exposed one extra weak alignment and 4,150,617 Rust raw pairs
versus 1,683,263 upstream pairs. The translated hash table incorrectly stored
the extracted bucket bits as the equality key; upstream uses those bits only
for initial placement and retains the full 32-bit seed for equality. After the
faithful correction, both implementations produce 1,683,263 raw pairs, 76,067
left-most survivors, and 1,241 byte-identical rows (`100681ab...`). The release
observation was 1.16 s/44,524 KiB upstream and 1.06 s/67,228 KiB Rust. The RSS
gap remains worth investigating; Rust's owned 2.0 MiB index read contributes,
whereas upstream mmaps the index. Complete artifacts are in
`.tmp/target-indexed-check/`.

To keep this work from narrowing onto one search option, the same build was
then checked on natural `deepclust` data. The 500 centroid/member pairs remain
exact after sorting (SHA-256 `46c21546...`); the post-hash-fix observation was
3.34 s and 32,444 KiB upstream versus 2.69 s and 7,492 KiB Rust. Artifacts are under
`.tmp/cross-tool-rotation/`.

The same rotation measured the already exact 509 MiB natural DAA viewer case
with output directed to `/dev/null`, alternating implementation order. C++ was
1.77--1.80 s at 841,192--883,284 KiB RSS; Rust was 6.85--6.86 s at
6,516--6,520 KiB. A 199 Hz cycle profile assigns 36.5% to transcript parsing
(`HspContext::parse`), 10.8% to escaped-title emission, and 5.8% to scanning
the packed transcript delimiter. The next viewer optimization should retain
the streaming memory advantage and remove repeated parse work rather than
materializing the archive like upstream.

### Shared hash-join correction and portfolio rerun — 2026-10-11

The translated open-addressing join had incorrectly retained the extracted
bucket fragment as its equality key. Upstream uses that fragment only to pick
the initial slot and stores and compares the complete 32-bit seed key. Rust
now does the same, including upstream's ceiling when sizing the table. A
collision regression covers distinct full keys sharing an initial bucket.

The release build was then rotated across command families rather than tuned
against a single executable path:

- Three alternating natural 1 Mb `blastx --mid-sensitive` pairs are exact at
  SHA-256 `c4b1cad2e798d3535b5f88ab4f57595f4212c0d33a0fa6c1f5ccafabcabfa3ca`.
  Medians are 4.95 s/268,808 KiB C++ and 4.18 s/188,392 KiB Rust.
- Three alternating natural `blastp --ultra-sensitive` pairs are exact at
  SHA-256 `2a2570ec2fe235c626b344dfc45e174d16db0305499535c1d02b8701518740a6`.
  Medians are 4.62 s/23,944 KiB C++ and 2.94 s/13,108 KiB Rust.
- Three alternating four-thread forced-disk `blastp --sensitive` pairs are
  exact at SHA-256
  `18c53b2d8afc218482cd9836e778f530964ab3ee39c3f27f130a3bbbacca63ad`.
  Medians are 28.75 s/88,108 KiB C++ and 34.02 s/70,936 KiB Rust. This remains
  an open stage-1 performance gap; the fresh result supersedes the earlier
  near-parity timing on the shared host.
- Natural 500-record `deepclust` remains exact after sorting at SHA-256
  `46c21546f09971d49798a6dbb40463f75a624d921b08553d913b5b44e13ea64c`;
  the scout was 3.34 s/32,444 KiB C++ and 2.69 s/7,492 KiB Rust.
- The 509 MiB DAA tabular-view case remains exact while exposing a separate
  streaming tradeoff: 1.77--1.80 s at 841,192--883,284 KiB C++ versus
  6.44--6.60 s at 5,972--6,508 KiB Rust after eliminating the per-match
  temporary context used only for transcript parsing. A full 458 MiB output
  still has SHA-256 `cf0162c04c022c78f110c8f6b9f073f9253343ac032f60f935044d2314b72dd8`.

Artifacts are retained under `.tmp/hashfix-blastx-mid/`,
`.tmp/hashfix-blastp-ultra/`, `.tmp/hashfix-blastp-forced/`,
`.tmp/cross-tool-rotation/`, and `.tmp/target-indexed-check/`.

Matched 99 Hz profiles of the forced-disk case are retained as
`rust-current.perf.data` and `cpp-current.perf.data` in
`.tmp/hashfix-blastp-forced/`. Rust attributes 55.6% self cycles to
`filter_partition_to_hits_avx2`; upstream attributes 68.7% to its fused
`ARCH_AVX2::stage1`. Rust additionally pays separately visible Hamming and
visitor costs while moving each survivor through a transient 32-byte
`SeedMatch`. Upstream's `HitField` instead gives `search_query_offset` compact
target indices directly. This is the next faithful stage-boundary correction;
the evidence does not support tuning spill I/O or changing SIMD dispatch.

The compact boundary is now implemented for ordinary search. Each callback
carries one shared query coordinate and 16-byte target candidates in upstream
tile order; the materialized `SeedMatch` adapter remains available for other
callers and persistent-index search retains its distinct stage-0 path. Three
alternating forced-disk pairs are byte-identical at SHA-256
`18c53b2d8afc218482cd9836e778f530964ab3ee39c3f27f130a3bbbacca63ad`.
Medians improve from 28.75/34.02 s to 27.44/31.32 s C++/Rust (speed 0.876x),
with RSS 89,452/70,832 KiB (ratio 0.792x). Artifacts are under
`.tmp/compact-stage1-forced/`.

The same release preserves the exact established hashes for natural
one-thread `blastp --ultra-sensitive` and `blastx --mid-sensitive`. A natural
`deepclust --approx-id 50 --member-cover 60` rotation also matches upstream
after sorting at SHA-256
`52fc1f944a6e477bff16b086e4d84df6dba11a3a03241370fc49474b3fc6d34c`;
Rust took 1.48 s and 7,392 KiB. The 509 MiB DAA viewer completed independently
at 7.58 s and 6,508 KiB; this single contended observation is a regression
gate, not a replacement for the repeated viewer benchmark.

A follow-up experiment hoisted reduced-match and fingerprint CPU dispatch into
the AVX2 worker with const-generic specialization. Although all parity gates
passed, the heavy Rust screen regressed from the accepted 31.32 s median to
33.15 s, so the experiment was removed. Rebuilding the compact-only source
restored an exact 31.27 s/70,244 KiB screen. The next search pass should use
the retained post-compact profile rather than reintroducing that specialization.

### Reusable DAA tabular context — 2026-10-11

The tabular viewer still called `DaaMatch::context` for every row after direct
HSP parsing had completed. That constructed an owned `HspContext` by cloning
all query frames, both titles, and the packed transcript for each of 1,678,045
alignments. Upstream advances one reusable match context. Rust now owns query
state once per query record, moves each parsed HSP into that context, and
reuses the target-title allocation.

Three alternating Rust-to-Rust pairs on the 509 MiB natural archive give a
7.59 s/6,508 KiB baseline median and a 6.99 s/6,252 KiB reusable-context
median (1.086x faster, RSS ratio 0.961x). A complete 458 MiB output remains
byte-identical to upstream at SHA-256
`cf0162c04c022c78f110c8f6b9f073f9253343ac032f60f935044d2314b72dd8`.
The final rebuilt release binary completed the same archive-to-`/dev/null`
sanity gate in 7.17 s at 6,564 KiB RSS. The full Linux suite (now 1,707 library
tests), Windows GNU and AArch64 cross-checks, formatting check, diff check, and
offline package construction all pass. Viewer optimization is closed for this
rotation; the next profiling pass must use a different command family.
Artifacts are under `.tmp/view-reused-context/`.

### Clustering ranking lifetime — 2026-10-11

The portfolio rotation expanded default `cluster` from 500 to 4,000 distinct
AMRFinderPlus proteins. Two initial exact Rust runs took 8.98 and 9.02 s,
versus 5.58 and 5.33 s upstream, while Rust retained the lower RSS. Both logs
localized the gap to the 15-shape Linclust40 round (5.5 s Rust versus about
2.9 s upstream), and paired `perf` profiles attributed 31.1% of Rust cycles to
the query-length `KmerRanking` path.

The translation constructed that ranking inside every joined seed bucket.
Upstream `search_shape` constructs it once per shape and passes the same
immutable object to all partition work sets. Restoring that lifetime reduces
the Linclust40 round to about 3.1 s. Three alternating one-thread pairs give
5.56 s/27,448 KiB upstream and 6.00 s/18,956 KiB Rust (speed 0.927x, RSS
0.691x). Three sequential four-thread pairs give 2.95 s/74,528 KiB upstream
and 2.78 s/20,556 KiB Rust (speed 1.061x, RSS 0.276x). Every output has SHA-256
`f79ce9267234ce989cf1c8d53f6ae566dea04263295c042065be0c2374cc3f4b`.

Natural 4,000-record controls remain byte-identical for `linclust` and
`deepclust`. The single-run observations were 1.23/1.64 s and 9.46/9.80 s
upstream/Rust respectively. The short `linclust` result is a candidate for a
larger calibrated pass, not an accepted benchmark conclusion. A one-thread
`blastp --ultra-sensitive` search and a 216,087-row DAA tabular view also
remain byte-identical, satisfying the cross-command rotation gate. Raw
timings, outputs, and profiles are under
`.tmp/rotation-cluster-4000/`.

### Exact bottom-N sketch selection — 2026-10-11

The next portfolio rotation scaled `linclust` to 7,000 distinct natural AMR
proteins. Its initial exact comparison was 2.01 s/28,492 KiB upstream versus
2.88 s/28,548 KiB Rust. Matched profiles measured 9.07B sampled Rust cycles
versus 5.46B upstream; Rust sketch construction plus its hand-translated
libstdc++ introsort accounted for about 3.33B cycles versus 1.51B upstream.

Faster mode retains only the bottom 21 hashes. Rust now selects the exact
cutoff hash first and sorts the retained set. If equal hashes straddle the
cutoff, it falls back to the existing libstdc++-exact introsort because its
unstable permutation determines which positions survive. Focused tests cover
both the unique-cutoff fast path and the split-tie fallback. A post-change
profile reduces total sampled cycles to 6.55B and removes the full introsort
from the significant hotspots. On 7,000 records the exact Rust median is 2.14
s, down from 2.88 s, versus 2.05 s upstream. Four-thread medians are 0.81 s
Rust and 1.07 s upstream.

The acceptance fixture combines the disjoint 7,000-reference and 3,000-query
natural AMR sets, with no duplicate identifiers. Three exact one-thread pairs
give 2.80 s/36,400 KiB upstream and 3.14 s/39,072 KiB Rust (speed 0.892x, RSS
1.073x). Three exact four-thread pairs give 1.42 s/57,792 KiB upstream and
1.19 s/40,036 KiB Rust (speed 1.193x, RSS 0.693x). Every output has SHA-256
`a6ea9d7874b51508bc61a3cd5e2794bace72170facf76f01d8d918f64136b48c`.
Natural `cluster`, `deepclust`, and `blastp --faster` controls remain
byte-identical. Artifacts and before/after profiles are under
`.tmp/rotation-linclust-7000/` and `.tmp/rotation-linclust-10000/`.

While constructing the heavier input, native `makedb --in /dev/stdin` exposed
an `Illegal seek`: gzip detection read two bytes and rewound the input. It now
peeks through `BufRead::fill_buf` and remains forward-only. A plain 10,000-
record pipe produces a DMND byte-identical to upstream (SHA-256
`c3fd018c658e3618abc7172951fcb2322077af0bb0f555690a5dc91f4ffe6956`),
and `linclust` from that database reproduces the acceptance output SHA above.
An actual gzip pipe also produces a database byte-identical to an uncompressed
file input. The timed results remain controlled: both binaries read the same
upstream-built DMND.

### Padded left-most SIMD blocks — 2026-10-11

The current four-thread forced-disk profile still placed 59.3% of Rust cycles
in the compact stage-1 filter. The pointer-based upstream implementation loads
the final partial `reduced_match` and `seed_mask` block as a complete 16-byte
SIMD vector from `SequenceSet` delimiter/perimeter storage, then masks unused
lanes. Rust instead scalarized every partial block. The backed left-most path
now verifies that the padded bytes exist, performs the same full loads, and
falls back to the bounded implementation only at a genuine buffer edge. Its
x86 SSSE3 implementation reduces the complete window in one call, matching
upstream's inlined ISA-level operation rather than redispatching each block.

The first matched sampled run reduced total Rust cycles from approximately
419.9B to 397.1B and reduced the compact-filter share from 59.3% to 53.7%.
Three subsequent alternating exact pairs measured 27.36 s/94,120 KiB upstream
and 27.85 s/69,588 KiB Rust (speed 0.982x, RSS 0.739x). Every output has
SHA-256
`18c53b2d8afc218482cd9836e778f530964ab3ee39c3f27f130a3bbbacca63ad`.
The 4,000-record `cluster` control retains SHA-256
`f79ce9267234ce989cf1c8d53f6ae566dea04263295c042065be0c2374cc3f4b`,
and the 24 MiB DAA tabular-view control retains SHA-256
`723b772a2b4aac2caf08a981c8b44b8d51bc57ae853aadabc0b655141fda8325`.
Artifacts and both before/after profiles are under
`.tmp/rotation-forced-search-current/`.

The first tabular-view profile attributed 25.9% of cycles to
`HspContext::parse`, with additional work in retaining all decoded matches and
generic string-pattern searches for every title. Matching upstream's streaming
iterator and scanning the fixed ASCII delimiters directly reduced three exact
Rust runs to a 0.67 s median; three interleaved upstream runs had a 0.16 s
median. RSS medians were 5,120 and 114,844 KiB respectively. This is a real
improvement from the 0.92 s initial Rust scout, but the sub-second fixture is
too short for an acceptance benchmark and the remaining 4.2x gap still needs
a 3--5 second natural workload. The complete quick producer/viewer matrix
passes 56/56 after the changes under
`.tmp/validation-view-streaming-title/`.

The subsequent steady-state viewer fixture is genuinely natural and
nonduplicated: 3,000 distinct AMRFinderPlus queries were searched against
7,000 distinct targets with upstream `blastp -k 0`, producing a 509 MiB DAA.
Upstream `view -k 0` emitted 1,678,045 rows (458 MiB), while Rust initially
emitted nothing. Normalizing the documented zero sentinel to unlimited at the
CLI boundary fixes both search and viewing. The two large TSVs now share
SHA-256 `cf0162c04c022c78f110c8f6b9f073f9253343ac032f60f935044d2314b72dd8`.
The hot DAA iterator was then translated more literally: it borrows the packed
transcript through parsing rather than growing one byte at a time and copying
it again. A low-load Rust scout improved from about 8.28 s to 8.16 s and used
6,512 KiB peak RSS, but the host became heavily contended. Unpinned upstream
runs ranged from 1.87 to 6.24 s and pinned upstream runs from 15.04 to 23.09 s;
pinned Rust ranged from 10.16 to 13.95 s. These raw values establish neither a
stable speed ratio nor an accepted optimization. Output remained exact in all
runs. A profile of the natural workload attributes about 30.5% of Rust samples
to `HspContext::parse`, 13.7% to raw-match decoding, and 8.5% to title escaping,
so future viewer work has named targets without monopolizing the broader port.
After the borrowed-transcript change, the complete producer/viewer regression
matrix passes all 56 rows under `.tmp/validation-view-direct-transcript/`.

The same `-k 0` normalization was checked independently in native `blastp` on
75 distinct queries against 1,000 targets. One- and four-thread C++ and Rust
outputs all contain 1,999 rows with SHA-256
`2934068a5b0fb70006a5495e3d0eabe1df3a51034161d1fa89ca96eded6eff99`.
The startup-scale scouts were 0.61 s/17,236 KiB C++ versus 0.28 s/8,192 KiB
Rust at one thread, and 0.38 s/28,896 KiB versus 0.13 s/11,308 KiB at four
threads; they are correctness scouts, not steady-state benchmark claims.

A separate natural search-option rotation then checked five independent
controls rather than remaining on viewer optimization. `-k 5` is exact at 202
rows (`50471250...`), `-e 1e-20` at 533 rows (`da63e217...`), `--xdrop 20` at
549 rows (`77046cda...`), and `--ext-chunk-size 16` at the same 549 exact rows.
The initial 500-by-itself `--no-self-hits` comparison exposed a real omission:
the translated final HSP filter had no target title and retained self matches.
The faithful fix compares the full titles and biological residues after
ranking/extension, ignoring internal mask annotation bits. It intentionally
does not backfill a discarded self target, matching upstream's stage boundary.
The corrected outputs contain 6,973 identical rows with SHA-256
`3bba2f79fe1c65eb3e0f6a4c03d05abb133d959cb1e9714f75ac63b6364b3b6d`.

Persistent-index construction and native consumption are now byte-exact. Both
implementations load `.seed_idx`, build the reference seed array, and build a
one-pass query seed array filtered by the stored target seed set before the
hash join. Rust therefore follows the same stage boundary rather than
approximating target-indexed search with the ordinary path.

### Bounded optimization rotation and scoring controls — 2026-10-11

Optimization passes now have an explicit stop rule: make at most one or two
profile-supported faithful changes in one command family, retain only changes
with a repeatable whole-command benefit, then rotate to another family. Two
experiments were removed under that rule. Translating the DAA HSP iterator to
upstream-style raw pointers preserved the 458 MiB tabular output exactly but
changed the large viewer median only from about 7.02 to 7.00 seconds. Replacing
the motif table's randomized hasher with an identity hasher also preserved
output, but hashbrown's power-of-two table exposed poorer packed-kmer bucket
distribution: `mask_motifs` grew from roughly 3.4% to 6.8% of sampled cycles
and one-thread `linclust` remained at 3.11--3.14 seconds. Both experiments were
reverted rather than accumulating unproven unsafe code or special cases.

The next rotation exercised the historically fragile BLOSUM50 and PAM250
scoring paths with 300 natural AMR queries against 7,000 distinct AMR targets.
Both outputs are byte-identical. The 6,460-row BLOSUM50 scout took 5.97 seconds
and 36,320 KiB upstream versus 5.00 seconds and 42,240 KiB Rust (speed 1.194x,
RSS ratio 1.163x), SHA-256
`56de01a50ac30860b11857b4bbb228a6d114a934fc7636079dfb2aa05758de89`.
The 69-row PAM250 correctness control took 1.84 seconds and 31,596 KiB upstream
versus 1.87 seconds and 15,824 KiB Rust, SHA-256
`fd52217edeedc1ca2841ead2900a79957558fd7389ec10154592baa4c85face8`.
These are single alternating scouts, not stable benchmark-table claims.
Artifacts are under `.tmp/rotation-motif-hasher/`,
`.tmp/view-pointer-iterator/`, and `.tmp/rotation-scoring-current/`.

A parallel clustering audit found that upstream dynamically refills completed
AVX2 endpoint lanes while Rust uses fixed-width batches. Translating that
behavior preserved the exact 10,000-record `linclust` output but changed the
three-pair Rust median only from 3.14 to 3.13 seconds. The experiment was
removed under the bounded stop rule; future work should compare the generated
i16 endpoint/profile kernel itself rather than revisit lane scheduling.

### Target-indexed memory-mode checkpoint — 2026-10-11

The next bounded search rotation scaled persistent-index consumption to 300
natural queries against 7,000 targets. Three default-memory pairs are exact
at SHA-256 `289cd6758fdee6b47ae7e0e9cda9d728e7383a7dc2763af15525d17021f6da0e`:
4.03 s/46,760 KiB upstream versus 4.64 s/82,908 KiB Rust. Three forced-disk
Rust runs (`--memory-limit 0G`) retain the same hash and give median results of
3.95 s/46,260 KiB upstream versus 4.69 s/68,360 KiB Rust. Thus the documented
16G in-memory default accounts for about 15 MiB of this RSS difference, but
does not explain the remaining 22 MiB or the 16--19% running-time gap.

The target-indexed adapter now passes compact stage-1 candidates directly to
the ordinary query-major ungapped filter, matching upstream's
`search_query_offset` boundary instead of materializing a second full
`SeedMatch` representation. Exact parity was unchanged; the whole-command
timing change was only about one percent and is not claimed as a speed win.
The bounded stop rule therefore moves optimization to clustering and DAA/view
before any further indexed-search change. Artifacts are under
`.tmp/rotation-target-indexed-compact/` and
`.tmp/rotation-target-indexed-0g/`.

A subsequent allocation audit found the remaining RSS cause: Rust accumulated
the complete 2.28-million-entry reference seed array in partitioned vectors
and then copied it into its final preallocated array. Upstream uses 16-entry
buffers to write directly into disjoint histogram-derived slices. Translating
that lifetime reduces forced-disk Rust RSS from 68,360 to 46,404 KiB versus
46,260 KiB upstream (RSS ratio 1.003x). Default-memory Rust falls from 82,908
to 64,316 KiB; its remaining excess is the documented in-memory hit store.
All repeated outputs retain the same exact hash. Forced-disk runtime is now
4.88 s versus 3.95 s upstream, so this closes RSS rather than speed. Artifacts
are under `.tmp/indexed-rss-agent/`.

A later bounded speed pass removed another adapter-only representation. The
stage-1 HitField indices now address the original borrowed packed target
locations directly, matching upstream's `JoinIterator` boundary instead of
copying query locations, target locations, and `Stage1Candidate` records.
Five alternating and reverse-order old/new Rust pairs in forced-disk mode all
emit 6,462 rows with SHA-256
`289cd6758fdee6b47ae7e0e9cda9d728e7383a7dc2763af15525d17021f6da0e`.
Median runtime improves from 4.83 to 4.06 seconds (15.9%; all five pairs),
while median RSS changes from 46,400 to 46,720 KiB. The packed borrow is
guarded by size/alignment/range assertions and focused layout, offset, and
ordering tests. This experiment is closed; results are under
`.tmp/indexed-zero-copy-agent/`.

A fresh post-integration three-pair comparison measures the complete current
implementations rather than old/new Rust binaries. Forced-disk medians are
3.94 s/45,504 KiB upstream and 4.01 s/46,404 KiB Rust (speed 0.983x, RSS
1.020x). Default 16G Rust takes 4.00 s/63,972 KiB against the same upstream
baseline; the extra memory is the documented in-memory retained-hit store.
All nine outputs have the same SHA-256. Artifacts are under
`.tmp/indexed-zero-copy-final/`.

### Ordered parallel DAA viewing — 2026-10-11

The large-view profile was based on an invalid execution-mode comparison:
upstream defaulted to 40 ordered workers on this host, while Rust accepted but
ignored `-p/--threads` and always rendered serially. Rust now follows the
upstream boundary with ordered 32-query batches on a bounded Rayon pool;
explicit one-thread operation retains direct streaming so it does not pay the
batch memory cost.

The complete 458 MiB tabular output has SHA-256
`cf0162c04c022c78f110c8f6b9f073f9253343ac032f60f935044d2314b72dd8`
at one thread, four threads, and the 40-thread host default. With output sent
to `/dev/null` to isolate decode/render running cost, three alternating pairs
give these medians:

| Threads | C++ time | Rust time | Speed | C++ RSS | Rust RSS | RSS ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.70 s | 0.56 s | 1.250x | 111,416 KiB | 6,172 KiB | 0.055x |
| 4 | 0.44 s | 0.27 s | 1.630x | 169,580 KiB | 14,144 KiB | 0.083x |
| 40 (default) | 1.11 s | 0.81 s | 1.370x | 459,716 KiB | 41,728 KiB | 0.091x |

The earlier 6.4--7.0 second Rust observations remain useful history for the
serial parser improvements, but no longer describe the default execution
mode. File-output timings are strongly affected by storage/page-cache state;
the table therefore reports repeated null-output running cost, while parity
was separately established from complete files. Artifacts are under
`.tmp/view-ordered-workers/`.

A bounded four-thread follow-up confirmed that scheduling, rather than a
different transcript kernel, was the remaining gap: both implementations spent
roughly 43--45% of sampled work in DAA match/transcript parsing. Increasing
Rust's barrier window from 32 records to 128 or 384 improved contended runtime
by about 24%, but increased RSS from about 33 MiB to 92--247 MiB. Both exact
experiments were reverted; artifacts are under `.tmp/view-p4-agent/`.

The accepted replacement is an ordered producer/consumer pipeline with a hard
cap of 32 resident records across queued, running, completed, and pending
work. Slots are refilled after each ordered write, eliminating the batch
barrier without widening the memory window. The complete 479,998,113-byte
output remains exact at SHA-256
`cf0162c04c022c78f110c8f6b9f073f9253343ac032f60f935044d2314b72dd8`.
Three fresh alternating upstream/Rust four-thread null-output pairs give
medians of 0.44 s/169,580 KiB upstream and 0.27 s/14,144 KiB Rust. An
independent alternating old/new Rust comparison under heavier host contention
improved 3.70 to 1.96 seconds and reduced RSS from 32,836 to 25,620 KiB.
Artifacts are under `.tmp/portfolio-view-*`.

### AVX2 i16 score-profile construction — 2026-10-12

The remaining one-thread `linclust` profile showed Rust building every i16
substitution-profile column with 32 by 16 scalar indexed loads. Upstream's
`ScoreVector<int16_t>` instead loads the biased byte matrix and uses AVX2
`pshufb`, unpack/mask, and bias subtraction. The direct translation is now
used by both full and banded standard-matrix i16 endpoint kernels; adjusted
matrices retain their scalar construction path.

All natural 10,000-record outputs retain SHA-256
`a6ea9d7874b51508bc61a3cd5e2794bace72170facf76f01d8d918f64136b48c`.
Alternating old/new Rust medians improve from 3.48 to 3.10 seconds (10.9%)
with flat RSS. Fresh three-pair one-thread medians are 2.97 s/36,452 KiB
upstream and 3.11 s/38,440 KiB Rust. Four-thread medians are
1.37 s/59,588 KiB upstream and 1.02 s/40,360 KiB Rust. Artifacts are under
`.tmp/rotation-linclust-i16-profile/`.

### Combined tabular scoring and coverage fields — 2026-10-12

A natural 300-query by 7,000-target correctness rotation combined `positive`,
`ppos`, `qcovhsp`, and `scovhsp` with both sequence identifiers. Upstream and
Rust emit the same 6,461 rows, SHA-256
`f49a90dd73d42a038a539e3a24270f09ac1d91fc14da2d3fc8831db85603d9fe`.
The correctness scout took 5.70 s/33,936 KiB upstream and 4.81 s/42,488 KiB
Rust; it is not promoted to the repeated benchmark table.

This pass also corrected an overextension in the Rust field parser. Pinned
DIAMOND 2.1.24 rejects BLAST aliases such as `qacc`, `qaccver`, `sacc`, and
`saccver`; they are enum placeholders but are absent from upstream's registered
tabular field map. Rust now rejects those names (and the other absent BLAST-only
aliases) instead of silently formatting them as `qseqid`/`sseqid`. Unknown
field names stay on the native parser so default builds report the field error
rather than attempting an unavailable optional C++ fallback. Artifacts are
under `.tmp/rotation-output-fields/`.

### Global minimum-length seed enumeration — 2026-10-12

The natural `blastp --min-query-len 200` rotation initially found five extra
Rust alignments. The upstream seed enumerator reads `config.min_query_len` for
every `SequenceSet`, so the nominally query-specific option also removes short
reference sequences. Rust had applied it only to query seed arrays. The native
single-thread, partition-range, sketch, count-reuse, and low-complexity paths
now all preserve upstream's global behavior, with a focused regression for
short queries and short references.

The corrected 300-query by 7,000-target run is byte-identical: 5,941 rows and
SHA-256 `1b7c7264343fd187dc0b562c355dbff28213741d707a16ad46c4b58c8fb481a8`.
The single correctness scout took 5.81 s/31,052 KiB upstream and
4.92 s/42,752 KiB Rust. It is not promoted to the repeated benchmark table.
Artifacts are under `.tmp/rotation-min-query-len/`.

### Portfolio rotation profiles — 2026-10-12

A bounded clustering rotation confirmed one-thread `cluster` as the largest
remaining clustering gap: an exact natural 4,000-record scout took
5.41 s/27,248 KiB upstream and 6.02 s/18,312 KiB Rust. `deepclust` was only
3.1% slower and the current 10,000-record `linclust` comparison is 4.7% slower.
Matched profiles attribute about 4.80 billion Rust cycles versus 3.55 billion
upstream cycles to full i16 endpoint SWIPE. The structural candidate is
upstream's `AsyncTargetBuffer`, which refills completed SIMD lanes immediately;
Rust still waits for the longest sequence in each fixed 16-target batch.
Because an earlier linclust refill trial barely moved runtime, any cluster-
specific translation retains a strict accept/revert threshold. Artifacts are
under `.tmp/portfolio-cluster-current/`.

The separate natural one-thread `blastx --mid-sensitive` check is not a current
performance gap: exact output took 4.54 s/268,716 KiB upstream and
3.75 s/188,252 KiB Rust, corroborating the repeated README row. Its Rust
profile nevertheless records seed-array construction as the largest hotspot
(24.5%); upstream's 16-entry partition-local buffered writer remains a
faithful candidate for a later mode where measurements justify it. The trial
was reverted rather than retained without validation. Artifacts are under
`.tmp/portfolio-blastx-mid/`.

The same contiguous natural 1 Mbp *E. coli* query was gzip-compressed to close
the compressed-`blastx` input gap without duplicating or otherwise synthesizing
records. At one thread and `--mid-sensitive`, both implementations emit the
same 25 rows and SHA-256
`c4b1cad2e798d3535b5f88ab4f57595f4212c0d33a0fa6c1f5ccafabcabfa3ca`,
which is also the uncompressed-input hash. The correctness scout took
4.58 s/268,660 KiB upstream and 3.76 s/188,064 KiB Rust. Artifacts are
`.tmp/compressed-blastx-*` and the compressed input is under
`.tmp/perf-portfolio-baseline/`.

Explicit motif masking was checked in both states on the natural 300-query by
7,000-target protein fixture. With `--motif-masking 1`, both implementations
emit 6,461 rows and SHA-256
`c21b52d3b766dadea854060aa1cdacc27cb20a41bd53d7eb2160bb9fec868cce`;
the scouts were 5.98 s/33,924 KiB upstream and 4.95 s/42,492 KiB Rust. With
`--motif-masking 0`, both emit 6,467 rows and SHA-256
`c790a10df97d4a546c25dddc473b23a7ba33ae69b90401a7be8311cc74ed71fc`;
the scouts were 5.77 s/34,924 KiB and 4.86 s/42,932 KiB. The six-alignment
difference proves the option changed real search behavior rather than merely
passing through an inert parser path. Artifacts and exact commands are under
`.tmp/rotation-motif-masking/`.

A separate natural 300-by-7,000 combined-filter case used `--id 50
--query-cover 80 --subject-cover 80 --top 20 --evalue 1e-50`. Upstream and
Rust are byte-identical at 39,516 rows, SHA-256
`45aae57232d485a21784f7b0692ebbef5c49e81c2cf79f10d3ffbdd8930601c0`.
The correctness scout took 6.25 s/42,832 KiB upstream and 5.74 s/42,688 KiB
Rust. Upstream leave-one-out controls establish that multiple filters are
active: omitting identity adds 654 rows, while omitting query coverage adds 10.
Exact commands and artifacts are under `.tmp/rotation-combined-filters/`.

### Multi-argument `makedb` semantics — 2026-10-12

Pinned upstream does not support concatenating reference files through
`makedb --in`: one occurrence followed by two paths is rejected with `Too many
arguments provided for option --in.` Repeating the option is assignment, not
append: `--in part-0.faa --in part-1.faa` builds only `part-1.faa`. Rust had
incorrectly accumulated repeated occurrences and then rejected them. Its CLI
parser now preserves last-occurrence assignment while retaining the upstream
multi-value rejection.

The natural 7,000-target FASTA was split at its 3,500th record into two
disjoint files whose concatenation is byte-identical to the original. Both
implementations' repeated-option commands build the same 3,500-sequence DMND,
SHA-256 `c7c3adbbe6f949a8d7caf4c507096483bd4dab8f6d1623053df5d4a92530f7e0`.
The equivalent last-input release scouts took 0.13 s/15,936 KiB upstream and
0.10 s/7,040 KiB Rust.
Artifacts and exact commands are under `.tmp/rotation-makedb-multi-input/`.

### Rejected bounded optimization scouts — 2026-10-12

Two source-faithful lifetime/scheduling candidates retained exact output but
failed the bounded performance gate and were fully reverted:

- Reusing target-indexed stage-1 fingerprint and hit-field storage like
  upstream's long-lived `WorkSet` produced 5.02 s/46,408 KiB at `0G`, versus
  accepted current runs of 4.85--4.90 s/46,288--46,720 KiB.
- Translating full-i16 dynamic lane refill specifically for the natural
  4,000-record one-thread `cluster` case produced 6.11 s/18,316 KiB, versus
  the accepted current 6.02 s/18,312 KiB.

These results prevent repeated tuning of candidates that do not improve the
whole command. Artifacts are under `.tmp/indexed-scratch-agent/` and
`.tmp/cluster-refill-agent/`.

The next bounded clustering hypothesis did pass the acceptance gate. Upstream
extracts reduced production seeds through a raw iterator after validating the
shape and reduction once. Rust now has an equivalent prevalidated weight-12,
base-10 path rather than rebuilding checked slices and indexing the shape for
every seed. Across five alternating old/new Rust pairs on the natural 4,000-
record case, every output had SHA-256
`f79ce9267234ce989cf1c8d53f6ae566dea04263295c042065be0c2374cc3f4b`.
Median time improved from 5.57 to 5.34 seconds (4.1%; the candidate won four
of five pairs); median RSS was 18,636 versus 18,948 KiB with overlapping
ranges and no new allocation. Hardware counters independently showed 6.3%
fewer cycles, 7.1% fewer instructions, and 19.2% fewer branches. A 2,000-case
randomized equivalence test includes masked residues. This closes the bounded
cluster experiment; the optimization rotation moves to other workflows.
Artifacts are under `.tmp/cluster-seed12-agent/`.

Fresh complete C++/Rust comparisons after integration confirm the whole-command
effect. At one thread, three-run medians are 5.48 s/25,820 KiB upstream and
5.36 s/18,632 KiB Rust. At four threads they are 3.03 s/80,184 KiB and
2.63 s/20,552 KiB. Every output has the same SHA-256 above. Artifacts are
under `.tmp/cluster-seed12-final/`.

### Utility-command portfolio rotation — 2026-10-12

A 574,627-record, 208,482,574-residue natural Swiss-Prot fixture exposed a
different command-level memory problem. Native `getseq` materialized every ID
and sequence, while upstream walks the sequence-data and position streams.
The translated streaming reader keeps only the current record. Its full
286,605,314-byte FASTA and a selected first/middle/last-record output are
byte-identical to upstream (SHA-256 `3c15b06c...` and `317bd165...`). Rust RSS
falls from 600,960 to 3,200 KiB and the streaming-only median improves from
1.82 to 1.71 seconds. A second bounded pass then restored upstream's reusable
`TextBuffer` boundary: Rust formats a complete record once instead of calling
through `dyn Write` per residue. Uncontended Rust runs improve to 1.05--1.13
seconds. In a fresh three-pair series affected by shared-host contention, all
pairs still favor Rust and medians are 3.48 s/137,740 KiB upstream versus
3.03 s/3,200 KiB Rust. Artifacts are under `.tmp/getseq-format-agent/` and
`.tmp/utility-final/`.

On the same input, `makedb` produces a byte-identical database (SHA-256
`fd2b4de8...`). Rust formerly kept both the full original record set and a
second fully masked copy. Upstream masks its owned `Block` in place; restoring
that lifetime reduces the Rust median from 562,768 to 350,400 KiB and improves
old/new Rust time from 14.89 to 14.59 seconds. Fresh alternating complete-run
medians are 16.34 s/402,172 KiB upstream and 15.83 s/350,720 KiB Rust (speed
1.032x, RSS 0.872x). Artifacts are under `.tmp/makedb-rss-agent/` and
`.tmp/utility-final/`.
`dbinfo` fields match except for upstream's cosmetic leading ANSI reset, and a
larger `merge-daa` benchmark was not fabricated from duplicated query sets.
The source fixture remains under `.tmp/utility-portfolio/swissprot/`.

The next distinct one-thread `linclust` hypothesis was rejected rather than
extended. Profiles assign 45.18% of samples to `full_i16_impl`, but reusing an
upstream-style batch-lifetime score profile changed the exact 10,000-record
median from 2.94 to 3.15 seconds (7.1% slower) with flat RSS. All six outputs
retained SHA-256
`a6ea9d7874b51508bc61a3cd5e2794bace72170facf76f01d8d918f64136b48c`;
the source was fully reverted. Artifacts are under
`.tmp/linclust-profile-lifetime-agent/`.

The portfolio then rotated to a fresh default one-thread `deepclust`
whole-command check on 4,000 distinct AMRFinderPlus proteins. Three paired
runs give medians of 9.42 s/28,152 KiB upstream and 9.12 s/18,688 KiB Rust
(speed 1.033x, RSS ratio 0.664x). All six outputs contain 4,000 rows and have
SHA-256
`671fb44c4c5729ed732386fe3a916e9ded661e66b7f6a2e6c2f35d54262bdf99`.
Rust ranged from 7.67 to 9.56 seconds on the shared host, so this closes the
parity/RSS check without claiming a stable small timing difference. Artifacts
are under `.tmp/deepclust-refresh-20261012/`.

A separate one-thread `linclust` pass revisited motif-table hashing after the
surrounding seed-path changes. Upstream hashes the encoded `u64` k-mer
directly, while Rust's general-purpose `HashSet` used keyed SipHash. Five
paired Rust A/B runs retain SHA-256
`a6ea9d7874b51508bc61a3cd5e2794bace72170facf76f01d8d918f64136b48c`
and improve the median from 2.90 to 2.86 seconds (1.4%) with flat RSS. Fresh
three-run whole-command medians are 2.79 s/36,092 KiB upstream and
2.86 s/39,080 KiB Rust. This newer result supersedes the pre-seed-path
experiment recorded above. Artifacts are under `.tmp/linclust-motif-hash/`.

The next bounded `linclust` correction removed a larger translation mismatch
at the command boundary. Upstream aliases its target block to its query block
for linear self-search, whereas Rust cloned all 10,000 records into distinct
query and target owners. The in-memory edge-search API now represents that
case explicitly and shares the sorted, masked sequence storage through seed
enumeration; ordinary two-input search is unchanged. Five controlled old/new
Rust runs retain SHA-256
`a6ea9d7874b51508bc61a3cd5e2794bace72170facf76f01d8d918f64136b48c`
and improve the median from 3.13 s/38,756 KiB to 2.61 s/33,980 KiB (16.6%
faster, 12.3% less RSS). A fresh three-pair whole-command comparison is also
exact and gives 2.93 s/35,932 KiB upstream versus 2.79 s/34,104 KiB Rust
(speed 1.050x, RSS ratio 0.949x). This closes the one-thread `linclust` lane;
artifacts are under `.tmp/linclust-shared-self/`.

The same shared self-search backend is exercised independently by all five
`deepclust` rounds. A balanced three-pair refresh on 4,000 natural proteins is
exact at 4,000 rows with SHA-256
`671fb44c4c5729ed732386fe3a916e9ded661e66b7f6a2e6c2f35d54262bdf99`.
Its medians are 9.21 s/27,912 KiB upstream and 8.90 s/17,152 KiB Rust (speed
1.035x, RSS ratio 0.615x), so that separate command has no current deficit to
chase. Artifacts are under `.tmp/deepclust-audit-20261012/`.

The distinct utility/index lane scaled `makeidx` to 150,000 nonduplicated
Swiss-Prot records (59,368,545 letters). Both implementations emit the same
134,217,792-byte index with SHA-256
`dbe81bcde14fd8792554d893ef9dfdfc51d2f493974bec777f899312a3ea33b7`.
The discovery baseline gave medians of 8.00 s/394,968 KiB upstream and
10.55 s/437,128 KiB Rust (speed 0.758x, RSS ratio 1.107x). A faithful
const-generic experiment
borrowed each sequence for the non-mutating hashed-seed path instead of
cloning it once per shape, but controlled old/new medians were only
9.70/9.59 seconds with flat RSS. The change was reverted rather than treating
the visible allocation count as the cause.

Profiling instead showed the dominant missing translation: upstream runs its
per-sequence soft-masking stage on `config.threads_` workers, while Rust ran it
serially and ignored `makeidx --threads`. Rust now discovers independent motif
ranges on a bounded Rayon pool and applies them in deterministic block order.
The explicit one-thread medians are 10.09 s/394,972 KiB upstream and
11.07 s/437,124 KiB Rust (speed 0.911x, RSS ratio 1.107x). With the host's
default 40 threads they are 8.21 s/394,972 KiB and 8.94 s/434,744 KiB (speed
0.918x, RSS ratio 1.101x), improving Rust's default median by 15.3%. Both modes
retain the exact index hash above. Evidence is under
`.tmp/makeidx-portfolio-agent/`, `.tmp/makeidx-borrow-agent/`, and
`.tmp/makeidx-current/`; accepted threaded results are under
`.tmp/makeidx-threaded/`.

The next bounded `makeidx` pass addressed RSS rather than extending masking
tuning. Rust decoded the entire DMND into owned `FastaRecord`s and then copied
them into `Block`, retaining both representations while building the index;
upstream streams directly into its block. The translated path now uses the
same record visitor as `getseq` and bit-casts each borrowed byte slice to the
alignment-one `Letter` representation only for the callback lifetime. Three
old/new default-thread pairs retain the index hash above, reduce median RSS
from 434,716 to 346,872 KiB (20.2%), and leave time effectively flat at
9.08/9.02 seconds. Fresh complete default medians are 8.09 s/394,968 KiB
upstream and 8.71 s/346,892 KiB Rust (speed 0.929x, RSS ratio 0.878x);
one-thread medians are 9.50 s/394,972 KiB and 10.50 s/349,440 KiB (speed
0.905x, RSS ratio 0.885x). Artifacts are under `.tmp/makeidx-streaming/`.

The later fixed-invariant hashed-enumerator translation keeps the generic
skip/filter path but hoists `makeidx`'s known no-skip/no-filter decisions out
of the seed loop. Three exact Rust A/B pairs improve 7.79 to 7.45 seconds at
flat RSS. A fresh three-run upstream median is 7.08 seconds/394,968 KiB versus
Rust 7.45 seconds/347,248 KiB (speed 0.950x, RSS ratio 0.879x). All outputs are
134,217,792 bytes with SHA-256 `dbe81bcde...`; artifacts are under
`.tmp/makeidx-specialization-accept/`.

A later fresh one-thread portfolio confirms that search is not the remaining
material regression. Natural `blastx --mid-sensitive` is exact at 25 rows and
SHA-256 `c4b1cad2...`; medians are 4.49 s/268,576 KiB upstream and
3.64 s/188,920 KiB Rust (speed 1.234x, RSS ratio 0.703x). Natural `blastp
--ultra-sensitive` is exact at 5,210 rows and SHA-256 `23fb031b...`; medians
are 16.07 s/29,872 KiB upstream and 15.45 s/22,208 KiB Rust (speed 1.040x,
RSS ratio 0.743x).

The same rotation initially left p1 `makeidx` as the sole material deficit: exact
134,217,792-byte indexes with SHA-256 `dbe81bcde...` measure
8.93 s/394,972 KiB upstream and 9.57 s/349,760 KiB Rust (speed 0.933x, RSS
ratio 0.886x). Isolated exact-output `perf` runs record 28.394 billion cycles
upstream and 30.318 billion in Rust. Seed enumeration is already essentially
at parity in absolute cycles (20.989 versus 21.167 billion, +0.85%). Motif
masking costs 6.337 billion cycles in Rust versus 5.247 billion upstream and
explains about 57% of the total gap; setup/streaming accounts for smaller
pieces. This rules out more speculative hash-loop tuning and selects masking
as the next bounded faithful-translation target. Portfolio artifacts are under
`.tmp/p1-portfolio-refresh/`; raw profiles are under
`/tmp/diamond-makeidx-profile/`.

The selected masking fix specializes only the identity-reduced k-mer advance,
making the fixed `% 20^7` modulus visible to LLVM just as upstream's template
does. Assembly now contains reciprocal multiplication rather than hardware
`div`. Three alternating exact baseline/candidate pairs improve Rust from
9.52 to 8.88 seconds (6.7%) at flat 349,760 KiB RSS, closing the fresh
8.93-second upstream p1 gap. A fresh three-pair default-thread comparison is
7.45 s/394,968 KiB upstream versus 7.61 s/347,236 KiB Rust (speed 0.979x,
RSS ratio 0.879x), below the material threshold. All twelve indexes retain the
same 134,217,792-byte SHA-256 above. Evidence is under
`.tmp/makeidx-const-modulus-20261012/`.

### CI portability audit — 2026-10-12

GitHub Actions run `38143645156` on the checkpoint commit passed Ubuntu,
Intel macOS, Linux musl, both 32-bit checks, and all five native delayed-spill
jobs. Its remaining Windows and macOS ARM test failures predate the current
worktree fixes. The current source contains no upstream `include_str!`/
`include_bytes!` dependencies and no hand-declared runtime `open` symbol;
SQLite resolves through bundled `libsqlite3-sys`.

Full linked cross-builds pass for macOS ARM64/x86-64, Linux musl, Linux i686,
and Windows GNU x86-64; Windows 32-bit source also passes `cargo check`.
Windows standard input/output now duplicate the process's inherited handles
instead of opening `CONIN$`/`CONOUT$`, so pipes, redirection, and headless CI
use the same stream selected by the parent process. The local full test suite,
package verification build, formatting check, and diff check pass. A hosted
Windows MSVC/macOS ARM rerun remains required after the worktree is committed;
cross-compilation is not claimed as native runtime evidence.

A fresh one-thread profile after streaming attributes 70.5% of Rust cycles to
the same random-probe hashed-seed construction that accounts for 68.6% of the
upstream profile; input streaming is below 0.5%. Mirroring upstream raw-pointer
loads by removing Rust bounds checks from the hash-table probe was tested as a
bounded unsafe A/B. Five exact pairs kept index SHA
`dbe81bcde14fd8792554d893ef9dfdfc51d2f493974bec777f899312a3ea33b7`,
but medians changed only from 10.61 to 10.52 seconds with mixed pairwise
direction and flat RSS. The unsafe change was reverted. Evidence is under
`.tmp/makeidx-table-probe/`; further work needs a stronger translation
hypothesis rather than additional probe-loop tweaking.

The later portfolio rotation independently tested the remaining pointer/load
difference in `HashedSeedIterator`: Rust validated the sequence boundary once
and used upstream-shaped unchecked loads inside the rolling-seed loop. All
runs retained index SHA-256 `dbe81bcde14fd8792554d893ef9dfdfc51d2f493974bec777f899312a3ea33b7`.
Five alternated p1 runs measured 10.49 s/349,760 KiB for the frozen baseline
and 10.51 s/349,440 KiB for the candidate. The change was fully reverted;
the residual `makeidx` gap is not caused by these bounds checks.

Two other lanes were refreshed before continuing index work. One-thread
`blastx --mid-sensitive` on a contiguous 1,000,000-nt *E. coli* query remains
exact at 25 rows (SHA-256 `c4b1cad2...`) and measures 4.61 s/268,892 KiB
upstream versus 3.73 s/188,252 KiB Rust. Pairwise `view -f 0` on the natural
3,000-query by 7,000-target DAA remains exact at 101,270,326 bytes (SHA-256
`00fb5456...`) and five-run medians are 1.11 s/102,720 KiB upstream versus
1.04 s/5,120 KiB Rust. Those lanes are closed rather than profiled further;
artifacts are under `.tmp/blastx-mid-refresh-20261012/` and
`.tmp/daa-view-refresh-20261012/`.

The one-thread `blastp --ultra-sensitive` lane was also enlarged so startup no
longer dominates: 300 distinct AMRFinderPlus queries (102,101 residues) against
1,000 distinct targets (339,123 residues). All six alternating outputs are
exact at 5,210 rows and SHA-256
`23fb031bc14fa3f994394616e812fdd71a19b7976cba35f97e8e0d41a46a7f58`.
Medians are 17.66 s/29,048 KiB upstream and 18.35 s/28,988 KiB Rust (speed
0.962x, RSS ratio 0.998x), although individual paired speed ratios span
0.898x--1.014x under shared-host noise. Profiles localize the remaining
bounded discrepancy to standard stage-1 survivor handling: upstream passes
compact target indices plus the original packed locations into stage 2, while
Rust materializes a 16-byte `Stage1Candidate` per survivor before walking it
again. The standard path now retains ordered `u32` target indices and borrows
the original target arrays through stage 2, matching upstream's boundary.
Three frozen baseline/candidate pairs are exact and improve 16.53 to 15.70 s
(5.0%) while reducing RSS from 28,924 to 28,056 KiB (3.0%). A fresh balanced
upstream/candidate series gives 16.16 s/30,056 KiB upstream versus
15.59 s/28,344 KiB Rust (speed 1.037x, RSS ratio 0.943x); all six outputs are
5,210 rows with the same SHA above. Evidence is under
`.tmp/blastp-ultra-300-refresh/` and `.tmp/blastp-ultra-zero-copy/`.

The target-indexed memory lane was refreshed separately on 300 natural AMR
queries against 7,000 targets. All nine ordinary C++/Rust-0G/Rust-16G outputs
are exact at 6,462 rows and SHA-256
`289cd6758fdee6b47ae7e0e9cda9d728e7383a7dc2763af15525d17021f6da0e`.
Medians are 4.28 s/46,732 KiB upstream, 4.34 s/46,404 KiB Rust-0G, and
4.34 s/63,912 KiB Rust-16G. Forced spill therefore has RSS parity; the default
fast path retains more memory by design. A like-for-like upstream `--hit-membuf`
control measures 4.30 s/58,848 KiB, narrowing the residual Rust excess to
5,064 KiB. Upstream stores its 895,237 retained hits in 16 coarse pre-reserved
packed bins, whereas Rust grows 300 independently allocated, aligned per-query
vectors. Translating that packed coarse-bin layout is the bounded default-RSS
hypothesis. Artifacts are under `.tmp/indexed-rss-refresh-20261012/`.

That layout is now translated: the in-memory path stores upstream-shaped
15-byte hits in 16 coarse pre-reserved bins and materializes one bin at a time
in stable query order. Three exact frozen old/new pairs retain SHA-256
`289cd6758fdee6b47ae7e0e9cda9d728e7383a7dc2763af15525d17021f6da0e`,
reduce median RSS from 63,520 to 56,776 KiB (10.6%), and leave time flat at
4.17/4.18 seconds. The prior like-for-like upstream `--hit-membuf` median is
4.30 s/58,848 KiB, so Rust is modestly better on both measures. A forced-0G
guard remains exact at 4.26 s/46,404 KiB. Evidence is under
`.tmp/indexed-rss-layout-agent-20261012/`.

Archive utilities were rotated in independently with a valid one-input
`merge-daa` rewrite of the existing 532,968,582-byte natural all-target DAA
(2,894 queries, 6,615 used targets). Combining the other local large archive
would duplicate all 2,894 query IDs, so it was deliberately not used as an
artificial multi-input benchmark. Every input/output and cross-tool hash is
`e0371a564e9165a5b1e0c6fdebbd3abd7d6da1c7ca065a739b8c519dd67de23b`.
Balanced three-pair medians under variable shared-host load are
6.34 s/143,824 KiB upstream and 5.51 s/8,212 KiB Rust (speed 1.151x, RSS ratio
0.057x); Rust won each paired time. No optimization profile was taken because
there is no Rust deficit. Artifacts are under
`.tmp/merge-daa-natural-refresh-20261012/`.

The defining multi-input path was then exercised without duplicating data.
The natural 3,000-query corpus was split into contiguous, disjoint 1,500-query
shards and upstream produced both DAAs against the same 7,000-target database.
The inputs are 172 and 338 MiB with distinct SHA-256 values. Both mergers map
6,496 and 5,664 per-file targets to 6,615 unique output targets and append
2,894 query records. All six merged outputs are 532,968,582 bytes with SHA-256
`7c168e7931d0289e7426a7571e302c35c30ef94dbe17b805853f5552f959c600`.
Alternating three-pair medians are 6.02 s/210,416 KiB upstream and
4.50 s/9,656 KiB Rust (speed 1.338x, RSS ratio 0.046x), although filesystem
contention produced wide 1.11--8.79 second ranges. Cross-view validation is
also exact: upstream on the Rust archive and Rust on the upstream archive emit
64,204 rows, 18,951,353 bytes, and SHA-256 `dde6d4f5...`. Evidence is under
`.tmp/merge-daa-two-natural-20261012/`.

The search rotation independently refreshed the last stale sensitivity
deficit. Natural `blastp --very-sensitive -p1` on 300 disjoint queries by
1,000 targets is exact at 5,183 rows and SHA-256 `eca7ad25...`; medians are
3.99 s/25,252 KiB upstream and 3.96 s/18,444 KiB Rust (speed 1.008x, RSS ratio
0.730x). The missing multithread cell uses 1,200 contiguous natural queries
outside the 1,000-record target shard so startup does not dominate. Its exact
24,403-row SHA-256 is `5962e399...`; p4 medians are 4.97 s/51,484 KiB upstream
and 4.97 s/40,572 KiB Rust. This closes the sensitivity lane and rotates work
away from search. Evidence is under
`.tmp/blastp-very-sensitive-refresh-20261012/`.

The four-thread delayed-spill `blastp --sensitive` lane was refreshed on the
natural 2,000-query by 4,000-target workload. All six outputs retain SHA-256
`18c53b2d8afc218482cd9836e778f530964ab3ee39c3f27f130a3bbbacca63ad`;
every current Rust run migrates about 2.22 million hits and reloads 37.5 MiB of
compressed spill. After absolute-location filtering and reusable worker
histograms, medians are 27.59 s/92,716 KiB upstream and 27.41 s/71,148 KiB
Rust (speed 1.007x, RSS ratio 0.767x). The immediately preceding profile put
spill compression and file I/O below 0.01% of CPU; the remaining work is
stage-1 filtering, not the disk path. Current artifacts are under
`.tmp/delayed-spill-hist-refresh/`; earlier diagnostic history remains under
`.tmp/spill64-audit-20261012/`.

That ISA-hoisting hypothesis was then tested and rejected rather than inferred
from the profile. On a shorter natural 300-by-4,000 `--sensitive -p4` scout at
64M, all six outputs retain SHA-256
`a60992426ffcf94e909913edef548d328d412fbd4a9f549322ac6fe57476534f`,
but baseline/candidate medians regress from 5.38 to 5.57 seconds (3.5%) with
overlapping RSS. The specialization was fully reverted and the expensive
30-second acceptance workload was not rerun. Future stage-1 work needs a
different hypothesis than merely removing the inner feature check.

The database-construction lane then tested the same 574,627-record natural
Swiss-Prot input through its real 93,392,974-byte gzip representation. Both
tools emit the established 293,486,822-byte DMND with SHA-256
`fd2b4de8e95dda1fcb901956aa2e88cf712142764b68a9076825c753c285be98`.
Profiling put about 83% of the old Rust runtime in TANTAN masking and only 5.4%
in gzip inflate. Upstream partitions independent sequences over its configured
workers; Rust both masked serially and accidentally defaulted `makedb` to one
thread. Rust now masks owned records with Rayon while preserving record order,
and the CLI default uses host parallelism while explicit `-p1` stays serial.
Three exact old/new default-thread pairs improve the Rust median from 17.49 to
4.85 seconds (3.606x) with RSS rising from 350,720 to 379,548 KiB for worker
scratch. Fresh balanced medians are 6.38 s/552,272 KiB upstream and
4.10 s/379,412 KiB Rust (speed 1.556x, RSS ratio 0.687x). The explicit-p1
candidate remains exact at 14.71 s/350,720 KiB. Artifacts are under
`.tmp/makedb-gzip-swissprot-refresh-20261012/`.

Native DAA production was measured separately on the contiguous 1,000,000-nt
*E. coli* `blastx --mid-sensitive -p1` workload. Producer medians are
4.62 s/268,764 KiB upstream and 3.62 s/188,564 KiB Rust (speed 1.276x, RSS
ratio 0.702x). The archives are deterministic but use different valid
encodings: upstream writes 270,419 bytes (`e51bdf4d...`) and Rust 266,099 bytes
(`f9271624...`). Both upstream and Rust viewers render both archives to the
same 25-row, 4,135-byte TSV with SHA-256
`c4b1cad2e798d3535b5f88ab4f57595f4212c0d33a0fa6c1f5ccafabcabfa3ca`.
A fresh direct-TSV control is 4.58/3.65 s upstream/Rust, showing that DAA
serialization adds no Rust deficit. This output-mode lane is closed without a
source change.

### Two-file `blastx` audit — 2026-10-12

Source inspection established that upstream treats two query files as
synchronized paired-end input rather than concatenated files. It interleaves
records by ordinal, requires equal record counts, accepts positional pairing
without title matching, searches/ranks each mate independently, and exposes
the opposite sequence through `full_qseq_mate`. Rust now implements that same
round-robin load and validation boundary. Mate relationships are stored as
indices into the existing source-query vector, so the fix does not duplicate
read data. This pass also fixed a latent formatter error: mate bases were
previously rendered through the amino-acid alphabet rather than upstream's
nucleotide `input_value_traits`.

Validation uses the unmodified 5,000-pair SRR1703350 FASTQ subset (4,715,168
bytes total) against 7,000 natural AMR proteins. The R1/R2 SHA-256 values are
`1b6d931d6b004dee5a379f8cb636afaee770c2ed4322002f264c049f6ff5666b` and
`f65b796f713840cfbc62c5991c18884ede9cc4e358c460076c9b0e0cdb6dd078`.
Output explicitly includes `qseqid full_qseq_mate sseqid bitscore`, so mere
concatenation cannot pass. Every one- and four-thread output has 3,081 rows and
SHA-256 `23c16303d397882e2090cfa3fc5fd1ffc083740a88f8495e71d4fdcba1a8508b`.

The accepted memory pass removes the Rust-only translated-FASTA round trip,
builds the extension copy directly in contiguous `SequenceSet` storage, and
does not allocate 60,000 empty TANTAN-span vectors when the clustering-only
override is disabled. Refreshed balanced three-run medians are
6.09 s/47,664 KiB upstream versus 5.64 s/49,716 KiB Rust at one thread (speed
1.080x, RSS ratio 1.043x). The final four-thread medians are 2.64 s/62,824 KiB
upstream versus 2.30 s/53,580 KiB Rust (speed 1.148x, RSS ratio 0.853x).

A matched `perf` run attributed 40.94% of Rust cycles to the per-chunk seed
count pass and 37.65% to fill, while upstream amortizes its histogram before
the chunk loop. An early worker-histogram candidate changed the output hash to
`ce94f0c2...` and regressed to 5.54 seconds, so it was fully reverted. The
accepted implementation instead ports upstream's worker-by-shape-by-partition
counts and prefix ranges, preserving exact fill order. Six Rust A/B outputs
retain SHA-256 `23c16303...`; runtime improves 3.04 to 2.30 seconds with RSS
changing 52,996 to 53,580 KiB. Artifacts are under
`.tmp/blastx-histogram-reuse-20261012/`.

### NCBI BLAST database search audit — 2026-10-12

The translated reader decodes `.tmp/real-taxonomy/db2.psq` byte-for-byte to
the two source FASTA sequences (297 and 339 residues), and all six taxonomy
columns match upstream for each taxid. The route is not accepted as parity.
The pinned C++ search emits blank subject IDs and scores `1324/179/1369/531`,
whereas both Rust and upstream's FASTA route emit the correct
`1480/1409/1706/1406`; requesting traceback-dependent fields makes upstream
fail with `Traceback error`.

Source isolation found undefined behavior in the upstream BLAST database
loader. `load_parallel` appends a decoded temporary `StringSet` into an empty
destination. `StringSetBase::append` removes the destination's only 256-byte
leading perimeter, but preserves limits whose first offset is 256, and also
assumes a trailing 256-byte source perimeter that the decoder has not
materialized. Consequently the search addresses shifted and out-of-bounds
sequence/title data. A diagnostic rebuild that materialized the source
perimeter changed the already-corrupt result, confirming that output depends
on the invalid memory contents; the probe was then fully reverted and the
original binary output restored byte-for-byte. Reproducing this behavior in
Rust would be neither deterministic nor a faithful portable translation.
Taxonomy-aware native routing therefore remains disabled and is tracked in
`GAPS.md` rather than being reported as semantic parity.

### Breadth checkpoint after absolute-location stage 1 — 2026-10-12

The retained stage-1 translation preserves upstream-style absolute sequence
locations through Hamming rejection and delays record/local-coordinate decoding
until a hit survives. The primary natural sensitive p4/64M scout remains exact
and improves 6.03 to 5.73 seconds with RSS falling 38,320 to 37,780 KiB.

Acceptance was deliberately cross-command. Paired natural p4 `blastx` retains
SHA-256 `23c16303...` in a new three-pair series: upstream/Rust medians are
2.64/2.30 seconds and 62,824/53,580 KiB after the reusable-histogram
translation. One-thread ultra-sensitive `blastp`
is exact at 16.02/15.72 seconds in the frozen-baseline A/B, with equal RSS.
The exact `linclust` control also improves in its sampled pair. The only
negative control is a separate contiguous p4 `blastx` case at 6.47/6.56
seconds (1.4%), with flat RSS; this is below the material threshold.

The same rotation rejected an exact but 20% slower fused-motif `makeidx`
candidate and fully reverted it. Independent current-release checks found no
optimization target in clustering or DAA/view: exact `cluster` is 2.2% faster
than upstream with 35.4% lower RSS; exact `linclust` is about 4% slower with
8.1% lower RSS; large unlimited tabular view is 1.35x faster at 0.036x RSS;
and pairwise view is 1.11x faster at 0.059x RSS.

The final all-shape worker histogram is a direct translation of upstream's
`SeedHistogram`/`build_iterators` boundary: one count traversal supplies stable
worker ranges for every index chunk. It removes three redundant count scans
without atomics or a merge copy. A separate sensitive p4 control remains exact
and improves 5.51 to 5.06 seconds. The heavier delayed-spill acceptance row is
also exact and reaches 27.41 seconds versus 27.59 upstream while using 0.767x
RSS.

### Search and clustering feature rotation — 2026-10-12

The missing steady-state `--min-orf` search check used one contiguous
1,000,000-nt *E. coli* query against 1,000 natural AMR targets with
`blastx --mid-sensitive --min-orf 30 --threads 1`. Three balanced alternating
pairs produced the same 25 rows and SHA-256
`a9bbeabb81fea0d921290ec6b866c64575b749c38db2a7fc457c2046d773d9ff`.
Upstream medians were 5.05 seconds / 266,664 KiB; Rust medians were 4.34
seconds / 188,568 KiB. This closes the earlier one-second, four-thread-only
coverage gap without identifying an optimization deficit.

An explicit clustering pipeline then exercised
`--cluster-steps faster_lin fast_lin linclust-40_lin default` on 4,000 distinct
AMR proteins. All six 4,000-line cluster files were byte-identical with
SHA-256 `f79ce9267234ce989cf1c8d53f6ae566dea04263295c042065be0c2374cc3f4b`.
Upstream medians were 5.49 seconds / 23,484 KiB and Rust medians were 5.43
seconds / 16,956 KiB. Parser inspection also found and fixed a fidelity bug:
upstream replaces a vector option on each occurrence, whereas Rust had
appended repeated `--cluster-steps` groups.

The attempted native `--aln-out` comparison correctly failed upstream with
`Alignment output to file is not supported without --parallel-tmpdir.` Rust
had accepted that combination and emitted a two-column approximation. Native
Rust now rejects it with the upstream error. The valid upstream external
`linclust --parallel-tmpdir --aln-out` route remains explicitly untranslated
and is recorded as a gap rather than being approximated. Detailed run records
are under `.tmp/blastx-min-orf-refresh-20261012/` and
`.tmp/cluster-features-refresh-20261012/`.

The subsequent option-wiring pass found that native clustering silently
ignored several already-translated controls. They are now forwarded through
the command boundary. A natural `--mutual-cover 80` run initially disagreed
because Rust length-sorted the aliased self-search block twice and derived the
minimum-length ratio without forwarding it to stage-1 dispatch. Removing the
second sort and selecting upstream's `stage1_self_mutual_cov` kernel restores
every per-round alignment/cluster count and the final 4,000-line output
(SHA-256 `f8419aebdf960702016839e12597b13bb4ef7e38894d47a1e24fb93bbeea091f`).
Three alternating medians are 7.19 s / 23,156 KiB upstream and 7.35 s /
18,740 KiB Rust. A genuinely nondefault step vector (`faster_lin fast_lin
default`) is also exact (SHA-256
`afbbf91e675f5db64bef459e4350a8fbf23521f9619ab32bf0b016311bf9dc95`),
at 2.49 s / 23,160 KiB upstream versus 2.29 s / 17,304 KiB Rust.

At this checkpoint, the paired one-thread `blastx --mid-sensitive` RSS excess persisted under
forced-disk `0G`, ruling out the adaptive in-memory hit store: exact medians
are 6.03 s / 47,784 KiB upstream and 5.69 s / 49,844 KiB Rust. A direct
frame-allocation reuse experiment remained exact and improved time slightly,
but raised RSS to a 52,828 KiB median, so it was fully reverted. The remaining
1.043x RSS ratio is attributed to the known dual ownership of per-frame record
vectors and contiguous stage-2 sequence storage; fixing it requires making the
contiguous block authoritative through mask/remove-mask, as upstream does.
The later single-owner follow-up below closes this historical gap.
Artifacts are under `.tmp/cluster-option-wiring-20261012/` and
`.tmp/blastx-paired-rss-isolation-20261012/`.

### Streaming makedb and non-tabular DAA workers — 2026-10-12

The database builder now follows upstream's bounded `load_seqs(1e9)` lifetime:
the record that reaches the one-billion-letter threshold completes its block,
then that owned block is masked in place, written, included in the chained
hash, and released. Only the 16-byte-per-sequence position trailer remains
until final header rewriting. A small-limit multi-block regression avoids an
artificial huge fixture and is byte-identical to the owned-record path. On
150,000 natural Swiss-Prot records (59.37 million letters), upstream and Rust
produce the same DMND SHA-256; upstream measured 4.44 s / 126,084 KiB and Rust
4.23 s / 96,000 KiB. A natural input crossing one billion letters remains the
boundary-scale RSS evidence gate.

The DAA breadth audit found that non-tabular formats still collected every
decoded match before applying the default `-k 25` and ignored the translated
ordered worker boundary. The faithful replacement streams the match iterator,
stops at the output range, groups 32 queries per task, keeps at most
`3 * threads` ordered tasks resident, and covers DAA, PAF, SAM, pairwise, XML,
null, and edge output. On the natural 509 MiB archive, three alternating exact
null-output medians are 0.54/0.26 s upstream/Rust at one thread and 0.46/0.21 s
at four threads, with Rust RSS ratios 0.314x and 0.207x. Exact PAF medians are
0.69/0.41 s and 0.42/0.21 s. Natural DAA, PAF, pairwise, XML, and edge hashes
match at one and four threads; SAM bodies match after excluding the literal
invocation stored in `@PG CL:`. Evidence is under
`.tmp/view-null-workers-20261012/`.

### Heavy search rotation refresh — 2026-10-12

Three distinct natural, nonduplicated lanes were rerun as three alternating
C++/Rust pairs after the utility and viewer changes. All six outputs in each
lane are byte-identical. One-thread `blastx --mid-sensitive` on the contiguous
1 Mb *E. coli* query is 4.52/3.47 s upstream/Rust with Rust at 0.703x RSS.
One-thread `blastp --ultra-sensitive` on 300 by 1,000 distinct AMR proteins is
17.04/15.88 s with Rust at 0.718x RSS. Four-thread
`blastp --very-sensitive` on 1,200 by 1,000 proteins is 4.60/4.52 s with Rust
at 0.837x RSS. No lane exposes a repeatable Rust runtime regression, so this
rotation deliberately made no search-kernel change. Raw evidence is under
`.tmp/search-heavy-refresh-20261012/`.
# Native view ordinal offsets (2026-10-12)

Native `view --qnum-offset/--snum-offset` now preserves upstream's unusual but
observable conversion boundary: command-line values are parsed as signed
`int64_t`, then converted to unsigned `OId` during addition. Consequently a
query ordinal of zero with `--qnum-offset -1` prints
`18446744073709551615`. Rust uses explicit wrapping addition so debug and
release builds agree.

Direct upstream/Rust comparisons on `.tmp/native-cpp-blastp-quick.daa` were
byte exact for both tabular (`-f 6 qnum snum`, offsets `-1`/`20`) and flat JSON
(`-f 104 qnum snum`, offsets `10`/`-2`). SHA-256 was respectively
`458a9c95a796f12e16b0df277a7d9ef9d41ffd04284c8ad705cf21c91263844c`
and `89012510a487702311a49f04ca111afa10aa52ada98b385fdeefa8ddfcc26935`.
Focused CLI-parser and renderer tests pass. Artifacts are in
`.tmp/view-offsets-20261012/`.

### Strict view boundary and stdout — 2026-10-12

Native `view` now applies upstream-style option grouping rather than scanning
known flags independently. Unknown options, forbidden joined long forms,
invalid field counts, the `--top`/`-k` conflict, out-of-range `--top`, invalid
headers, and unsupported formats produce the same success/failure decisions.
C-compatible numeric prefixes, repeated-option last-wins behavior, `--top 100`
as unlimited, bare `--header` as verbose, and the version/invocation header
are preserved. Direct top/offset/header comparisons are byte-identical.

Omitting `-o` now uses stdout. Direct natural-fixture tabular stdout is exact
(SHA-256 `fefe587d...`), as is binary DAA stdout redirected to a seekable file
(`f316eb18...`).

### View help and console presentation — 2026-10-13

Native `view` now prints the pinned upstream option and output-field
documentation, including the same spacing and final blank line. Direct help
comparison is byte-identical after excluding the ANSI reset caused solely by
upstream's date-dependent expired-build warning. The missing-DAA diagnostic
also uses upstream's exact wording.

Console routing now follows upstream: stdout output is quiet by default, file
output uses normal verbosity, `--quiet` suppresses normal messages,
`--verbose` adds the archive build version, and `--log` overrides `--quiet`
and adds invocation/configuration diagnostics. Selection-limit messages for
the default, `--top`, and unlimited `-k 0` cases are emitted in the same
position after archive metadata. Timer-duration lines and the expired-build
warning remain intentionally outside the comparison. The self-contained
`view_stdout` integration tests and six focused `commands::view` unit tests
pass without reading vendored upstream source files.

The metadata path also reconstructs `ScoreMatrix` before printing Lambda and
K, matching upstream instead of echoing the raw DAA header fields. A synthetic
archive whose stored values deliberately differ verifies the upstream-derived
`BLOSUM62 0.267/0.041` presentation exactly.

### Taxonomy-enabled makedb — 2026-10-12

The translated builder now extracts accessions while each bounded FASTA block
is resident, externally sorts accession relations, performs the sorted taxid
join, streams the compact per-OID taxid array, and writes nodes/ranks and
scientific-name trailers before rewriting both headers. Temporary runs use the
upstream 2 GiB allocation threshold, portable files, in-place unstable run
sorting, and heap merging. Incomplete databases and abandoned runs are
removed on error.

Small fixtures are byte-identical to upstream for parsed mapping-only,
`--no-parse-seqids` mapping-only, nodes-only, names-only, and all-three input.
The all-three SHA-256 is
`9e5d89bbbf2c5840d7aa98b9f4d1657883891fea49dc86af6432855ccf877a9d`.

The acceptance workload uses 150,000 contiguous, nonduplicated Swiss-Prot
records (59,368,545 letters) and derives its real accession/taxid mapping from
the unchanged UniProt headers. Three alternating one-thread pairs are exact,
SHA-256 `61d7d6db2f5b7bd26c1916abb098c031ae2400c3d80c25870559caa319d277bc`.
Median upstream/Rust measurements are 4.39/4.46 seconds and 126,172/105,600
KiB (speed 0.984x, RSS ratio 0.837x). Each implementation opens the other's
database. Evidence and commands are under `.tmp/taxonomy-natural-150k/`; the
fixture generator is `scripts/prepare_taxonomy_benchmark.py`.

### Native taxonomy search output — 2026-10-12

Native DMND search now derives `SequenceFileFlags` from the requested output
fields and loads only the required mapping, nodes/ranks, and names sections.
It reconstructs the same ancestor/rank/name tree used by upstream and supplies
per-subject taxids to both protein and translated search. Ordinary searches
retain the lean reader rather than allocating an empty 24-byte taxid vector
for every subject.

Direct C++/Rust output is byte-identical across `qseqid`, `sseqid`, `staxids`,
`sscinames`, `sskingdoms`, `skingdoms`, `sphylums`, `slineages`, and `score`.
Blastp SHA-256 is
`d5143565a6fb4cea25f31362a42881f095aaec97f2b3d96aa002984296720e54`;
blastx SHA-256 is
`cd6740dab8deea280425198265d2e9111069666aa44f552cc8ddb36356e8b41b`.
The Rust-built database gives the same blastp result, and missing mapping data
fails in both implementations. Evidence is under
`.tmp/taxonomy-search-parity/`.

### Remaining ordinary cluster controls — 2026-10-12

`--evalue`, `--comp-based-stats`, `--masking`, `--soft-masking`,
`--motif-masking`, and `--ext` now pass through the CLI, wrapper, cascaded
state, every round, and the in-memory search backend. Defaults remain hard
masking none, soft TANTAN, CBS mode 1, E-value 0.001, and sensitivity-derived
motif/extension behavior. Non-final rounds cap a user E-value at 0.001 while
the final round restores the requested cutoff.

On 4,000 distinct natural AMR proteins, `cluster --evalue 1e-20 -p1` is exact
in all six alternating runs (SHA-256
`cebeeb7ad240379fb8e06a336e4b4e02e819447116f65b320e38b69b7275de05`).
Median upstream/Rust measurements are 5.38/5.29 seconds and 25,992/17,816
KiB (speed 1.017x, RSS ratio 0.685x). Independent 500-protein controls for
CBS 0, hard TANTAN, soft masking 0, motif masking 0, and full extension are
also byte-identical. Evidence is in `.tmp/cluster-controls-accept/`.

### Packaged all-features policy — 2026-10-13

The optional C++ `ffi` feature remains usable as a conformance diagnostic from
a repository checkout containing `diamond/`. The build script now exposes the
adapter only when that source tree is present. The intentionally source-free
crate archive therefore compiles with `--all-features` and reports `--legacy`
as unavailable instead of failing during CMake configuration.

Validated with an offline `cargo package --locked --allow-dirty`, followed by
`cargo check --locked --all-features --all-targets --offline` inside the
expanded 388-file, 6.1 MiB package. A packaged `cargo run --features ffi --
blastp --legacy` emitted the expected unavailable-adapter diagnostic.

### Forced-split linclust target-linear parity — 2026-10-13

The low-memory centroid search now carries upstream's `lin_stage1_target`
state into stage 1. It therefore dispatches `stage1_target_lin` and selects the
first occurrence from the length-sorted target seed group, rather than
evaluating Rust's former full query/target cross product. The split wrapper's
`freq_masking` state is also propagated into motif-mask selection.

On 4,000 distinct natural AMR proteins, `linclust -p1 -M 5M` forms the same
two super-blocks in both implementations. Three alternating pairs all produce
the same 4,000-row file with SHA-256
`254ff1333d763d0a9c86ef94e4eae6550224e4fd66d3fe7d57b3fe7be8bfaa8f`.
Median upstream/Rust measurements are 1.59/1.28 seconds and 25,460/16,320 KiB
(speed 1.242x, RSS ratio 0.641x). Evidence and per-run files are under
`.tmp/cluster-lowmem-agent/`.

### Current CI fixes and isolated-package tests — 2026-10-13

The remaining hosted Windows x86-64 test failures were test-path portability
bugs rather than production linker failures. Database-builder tests no longer
embed Rust test names containing `::` in Windows filenames, and the multinode
FASTA test constructs `reps/0.faa` from path components instead of a literal
slash-containing component. The focused suites pass, as does a fresh
`x86_64-pc-windows-gnu --all-targets` cross-check.

The generated, expanded package now passes the complete offline
`cargo test --locked --all-targets` run: 1,736 library tests pass, six are
ignored, and every runnable binary/integration target passes. This verifies
that routine tests no longer depend on the checkout-only upstream source tree.

### AVX2 ungapped padded-tail translation — 2026-10-13

The upstream AVX2 ungapped kernel transposes a complete 32-byte vector for the
final partial block of its fixed 48-byte window, relying on `SequenceSet`
perimeter padding. Rust instead copied the final 16 bytes into 32 temporary
lane buffers. The translated kernel now uses the upstream full-width path when
the borrowed slices expose that padding and retains the copied-tail fallback
for generic exact-length callers. A direct regression proves both paths return
the same scores.

One bounded natural `blastp --ultra-sensitive -p1` A/B remained byte-identical
to upstream at 5,210 rows (SHA-256 `23fb031b...`) and changed 16.65 to 16.57
seconds while reducing RSS from 20,900 to 20,584 KiB. This single contended
pair is evidence against regression, not a stable performance claim. Raw
evidence is under `.tmp/blastp-padded-tail-results/`.

### Paired blastx single-owner translation — 2026-10-13

The translated-search query lifetime now matches upstream: a single
authoritative `Block` owns each query's source DNA and six translated contexts,
masking and restoration mutate that storage in place, and extension consumes
borrowed slices instead of a duplicate owned frame set. The restoration logic
also uses the distinct `SEED_MASK` and `LETTER_MASK` semantics, so propagated
seed masks survive the soft/motif mask cycle. The motif-disabled natural
multi-shape control is byte-identical at 25 rows with SHA-256
`c4b1cad2e798d3535b5f88ab4f57595f4212c0d33a0fa6c1f5ccafabcabfa3ca`.
The target-indexed control is likewise byte-identical at 25 rows with SHA-256
`7890ad7a0ab136ad2e7d2255e7221d5d3741c854bf4eedc6cddf93338224fc18`,
guarding both the repeated mask lifecycle and the alternate search path.

Three balanced repetitions used 5,000 natural SRR1703350 read pairs against
7,000 distinct AMR proteins. All nine upstream, Rust-default, and Rust-0G
outputs contain 3,081 rows and have SHA-256
`23c16303d397882e2090cfa3fc5fd1ffc083740a88f8495e71d4fdcba1a8508b`.
Medians are 6.38 s / 47,784 KiB for upstream, 6.07 s / 36,656 KiB for Rust's
default 16G mode (speed 1.051x, RSS ratio 0.767x), and 6.49 s / 37,408 KiB for
Rust forced-disk 0G (speed 0.983x, RSS ratio 0.783x). The single-owner rewrite
therefore removes the former paired-blastx RSS excess in both memory modes,
while forced spilling remains within 2% of upstream runtime on a shared host.
Commands, full hashes, inputs, load notes, and raw measurements are in
`.tmp/blastx-paired-single-owner-20261013-1154/`.

