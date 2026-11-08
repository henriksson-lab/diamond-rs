Native search-output snapshots from DIAMOND 2.1.24, pinned commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, generated 2026-11-07.
These are outputs from the unmodified executable, not repaired diagnostics.
Cargo tests need neither the upstream checkout nor FFI.

The compact inputs extend `../multi_hsp/` with XML-special title characters
and a low-complexity query that produces no retained trace points. The domain
query produces several disjoint HSPs. DNA inputs cover both translated strands.
The coverage variants reject alignments after seeding, distinguishing emitted
unaligned-query records from entirely empty source hit bins, which upstream
never schedules. The ordinary alignment worker also skips individual queries
with no retained seed hits inside nonempty bins. `--unal 0/1` changes the former and does not change the latter.

Natural query inputs are eight full-length SHV proteins drawn from
`../reassign/natural.faa`; the DNA input encodes the first four on both strands.
That shared 36-protein FASTA is the natural reference. The natural comparisons
exercise many target hits and source ordering, not only self matches.

For a snapshot, build the relevant reference with `makedb --in REFERENCE
-d db -p1`, then run `blastp` or `blastx -q QUERY -d db.dmnd -p1 --masking 0
-f FORMAT`, adding the options in `tests/native_search_output.rs`.
Pairwise uses `0`, XML uses `xml`/`5`, and PAF uses `paf`/`103`.
XML snapshots replace only the runtime `<BlastOutput_db>` path with
`@DATABASE@`. All other output bytes remain unchanged. The external comparison
uses the same real path in both commands and compares every byte directly.

The development checkout retains the full 96-case command/output comparison
and final release replay under `.tmp/native-search-output-20261107/`.
It covers both commands, CBS 0/1, one/unlimited HSPs, coverage filters,
unaligned reporting, one/four threads, and native 16G/forced-disk 0G storage.
The source and native CLI differ on permitting the native memory-limit
extension for these searches; that option is supplied only to the native
comparison command. The upstream executable uses its ordinary buffering.

`small-evalue.xml` additionally checks iostream six-significant-digit header
formatting at `--evalue 1e-10`. Escaped first-query titles detect the previous
header truncation, and blastx titles detect loss of full query descriptions.

`bins-query.faa`/`bins-query.fna` contain 44 alternating long, short and
unseeded queries. This exceeds the 16-bin default and forces heterogeneous
length-weighted boundaries. Its snapshots prove that unseeded queries remain
absent even beside seeded queries, while original XML query ordinals retain
those gaps. Both memory and disk output are checked.

`tie-reference.faa` and `tie-query.fna` use Python Random seed 5238 to create
an exact 1000-residue query and two references (1000 and 1500 residues). Both
hits have raw score 5748 and zero E-value. At mutual 50% coverage, source
blastx must preserve the shorter-first input order: the automatically derived
protein length-ratio optimization does not apply to translated search.
XML/PAF snapshots retain that ordering; the external sweep also checks tabular
and DAA bytes at one/four threads. This is a defined tie, not an oracle repair.
