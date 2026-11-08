Upstream DIAMOND 2.1.24 snapshots, pinned commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, generated on Linux AVX2 on
2026-11-07. Native integration tests require no C++ tree or FFI.

The compact fixture combines the reordered-domain query and two references
from `../multi_hsp/`. The dynamic fixture uses its single domain as centroid
with 42 deterministic complex flanking sequences (Python Random seed 712).
It exercises more than one SIMD vector of targets, lane refills, byte-score
saturation and promotion, and variable subject completion order. These are
algorithm regressions, not natural-data performance fixtures.

Run `realign -d INPUT --clusters CASE.clusters.tsv --header -p1
--comp-based-stats CBS -M MEMORY -o OUTPUT`. Input can be the FASTA or a
database built with upstream `makedb`. Default output has nine fields; `stats`
uses the fields in `tests/native_realign.rs`, and `ids` uses
`--outfmt 6 qseqid sseqid score`. Snapshots cover CBS 0/1 and whole/split
database loading. Upstream accepts 2–5 but its realignment DP enables Hauser
correction only for code 1, so the compact test also checks those codes against
the mode-0 snapshots.

AVX2 hosts compare exact bytes, including member ordering. Other hosts compare
complete rows: upstream's asynchronous subject lanes have different widths
on SSE/NEON and therefore can finish in a different order. These fixtures do
not establish exact output-order parity against an ARM upstream executable.

The upstream file-merging comparator treats equal centroid keys as less than
each other. Rust deliberately uses explicit heap operations matching the
pinned Linux executable, with a separate five-stream equal-key regression.
This fixes the otherwise different split-block ordering; it is not a claim
that an invalid comparator has portable ordering across C++ standard libraries.

Extra-field snapshots omit the input/output header and request signed ordinal
ID offsets, full titles, coverage, corrected bitscores, and normalized identity.
Upstream serializes and deserializes the HSP contexts before invoking the
formatter: query sequence data is absent and corrected bitscores reset to zero.
Consequently normalized identity uses the subject length in this workflow.
The native callback now preserves that intermediate-record boundary.
Upstream does not initialize `hsp_num` on this path (the local executable
happened to emit 32767); Rust deliberately emits its deterministic zero default.
The integration test checks that choice separately from exact upstream fields.
