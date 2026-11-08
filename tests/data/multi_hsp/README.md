Pinned upstream DIAMOND 2.1.24 snapshots, generated 2026-11-07 from baseline
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`.

The synthetic query contains three complex domains separated by unknown
residues. The principal target contains the same domains in reverse order,
forcing three disjoint HSPs rather than one collinear alignment. A second
target is a single-domain control. The DNA queries encode the same protein
on the forward and reverse-complement strands; unknown residues use NNN codons.
This is an algorithm fixture, not a benchmark.

Build with `makedb --in reference.faa -d db -p1`. Run `blastp -q query.faa` or
`blastx -q query.fna` with `-d db -p1 --masking 0 --max-hsps MAX
--comp-based-stats CBS`, using the fields in `tests/native_multi_hsp.rs` for
tabular snapshots or `--outfmt 100` for archives. MAX is 0/1/2/3 and CBS is 0/1.
Each query produces four rows for unlimited/three HSPs, three rows for two HSPs,
and two rows for one HSP, proving that the limit changes the result. The default DAA
dictionary contains IDs only; unused nucleotide header fields are zero in the
pinned comparison executable (see `../seg/README.md` for the compatibility
decision around upstream's uninitialized fields).

The `reduced-*` regression uses a complete AMR NarA query (QHA94808.1) and
CplR target (AJP12342.1), plus a CesB query segment (ABD14712.1 residues
811--1241) and ClbB target segment (AMQ58400.1 residues 2031--2449) with a
synthetic complex anchor. The segments retain flanking context from the
natural 100-query/1000-reference scout. It checks that an equal-score
full-matrix alternative does not displace the earlier banded alignment, and
checks transcript-derived statistics for a weaker disjoint HSP. These are
constructed regression inputs, not a natural performance workload.
Generate `reduced.tsv`/`reduced.daa` with the same fields and `blastp
--more-sensitive --comp-based-stats 1 --max-hsps 0 -p1`.
