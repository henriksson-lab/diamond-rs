Native reassignment regressions for pinned DIAMOND 2.1.24 commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, generated 2026-11-07.
Cargo tests use the included data and require no C++ sources or FFI.

Nonempty-search snapshots are named `diagnostic.tsv` deliberately: the
unmodified upstream workflow fails with `Block::ids()` or a segmentation fault
because its temporary member database lacks the sequence/title load flags.
It also tests an unsigned sentinel with `>= 0` and can index before the centroid
array for an unmapped member. The included `oracle-repairs.patch` records the
two diagnostic-only call-site repairs. Only that reassignment object was
recompiled; the executable uses the pinned build's remaining objects.
These snapshots establish conformance to the explicitly repaired diagnostic,
not byte parity with the unmodified command. Rust keeps an unmapped member's
previous assignment instead of reproducing invalid C++ memory access.

`reference.faa` uses Python Random seed 518 to create two related representatives
and four misplaced members. One representative extends an exact member by 70
residues, while the other has mismatches. Member-only and mutual coverage must
select different representatives. Members of different lengths also catch
lost original OIDs when native search length-sorts queries and targets.
The PAM250 case finds no qualifying hits and checks unchanged assignments.
The natural fixture contains 36 distinct SHV, TEM, and CTX-M AMR proteins,
with six representatives and 30 deliberately misplaced members.

Run `reassign -d INPUT --clusters CLUSTERS --header -p4 -o OUTPUT` with
the additional options listed in `tests/native_reassign.rs`. Build a DMND
input with `makedb --in INPUT -d db -p1`. Native tests cover FASTA/DMND,
one/four threads, and 16G/forced-disk 0G storage.

`singletons.upstream.tsv` is from the unmodified upstream executable. It
contains no non-representative members and therefore avoids the broken search
integration and sentinel path. Configuration-error snapshots in the tests
were also verified against the unmodified executable.

Soft-SEG snapshots use one additional explicit diagnostic repair, recorded in
`oracle-seg-normalization.patch`: clear residue annotation bits before calling
vendor SEG. A signal-handler backtrace of the two-repair diagnostic reached
`s_ShiftWin1`, whose 128-entry alphabet is indexed directly by signed sequence
bytes; the diagnostic also observed seed-tagged input bytes. Normalization
eliminates this invalid indexing and yields matching one/four-thread output.
Rust already normalizes residues on this path. This extra repair is separate
from the ordinary reassignment diagnostic's two workflow fixes.
