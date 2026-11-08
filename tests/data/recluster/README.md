Native recursive clustering regressions for DIAMOND 2.1.24, pinned commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, generated 2026-11-07.
Cargo tests use these snapshots without C++ sources or the FFI feature.

The `reference` fixture extends the reassignment fixture with two related
180-residue proteins (Python Random seed 731, every fourth residue changed in
the second protein) and an unrelated 150-residue singleton. All three start
in the wrong clusters. They fail initial realignment and centroid search;
the workflow cascades the unmapped subset, recursively refines it and maps
the resulting local OIDs back to the original database. At 50% identity the
new family joins; at 90% its proteins remain separate. The original compact
records also exercise reassignment and asymmetric coverage.

The natural input and initial clustering are shared with `../reassign/`:
36 distinct SHV/TEM/CTX-M proteins with 30 deliberately misplaced members.
Unlike `reassign`, upstream `recluster` enables automatic iterated search.
Its early linear rounds choose different assignments from a direct search
at the final sensitivity. The natural snapshots detect that difference.

Snapshots named `diagnostic.tsv` come from an isolated executable with the
two repairs in `oracle-repairs.patch`: load the temporary query sequences
and IDs, and compare the unsigned unmapped sentinel correctly. The original
command otherwise fails with `Block::ids()` or invalid memory access before
completing ordinary recursive inputs. Only the recluster object was replaced;
all other objects come from the pinned build. These snapshots are explicitly
qualified diagnostic conformance, not unmodified-upstream parity.

`valid.upstream.tsv` comes from the unmodified executable. Its input is
`../reassign/reference.faa` and `valid.clusters.tsv`; every member already
qualifies against its assigned centroid, so no query search is needed.
Configuration errors in `tests/native_recluster.rs` also match the unmodified
executable. The native implementation retains checked mappings instead of
reproducing invalid upstream memory accesses.

Run `recluster -d INPUT --clusters CLUSTERS --header -p1 -o OUTPUT`, adding
the variant options in `tests/native_recluster.rs`. Tests cover FASTA and
DMND inputs, one/four threads, and 16G/32K memory budgets. Build the DMND
input using `makedb --in INPUT -d DB -p1`.

The reproducible comparison/build drivers, commands, input/output digests,
unmodified failures and final native replay evidence are retained under
`.tmp/native-recluster-parity-20261107/` in the development checkout.
