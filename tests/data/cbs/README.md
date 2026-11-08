Upstream DIAMOND 2.1.24 snapshots, pinned commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, generated 2026-11-07.
The tests build native databases and need no C++ sources or FFI.

`compact` uses the reordered-domain FASTA files in `../multi_hsp/`. Snapshots
cover composition modes 0–5 with one/unlimited HSPs, masking disabled, and
default sensitivity. Matrix adjustment changes the scores relative to mode 0.

`scaled` uses the included natural AMR sequences: CopS query AAA25804.1 and
three histidine-kinase subjects from the AMR reference corpus. Searches use
`--more-sensitive --masking 0`, composition modes 0–5, and one HSP. This case
exposes the distinction between the integer scoring matrix and the separately
rounded frequency-ratio matrix used by mode 5 to estimate its lambda ratio.
Rebuilding the latter as integer scores times scale changes the output.

Run upstream `makedb --in REFERENCE -d db -p1`, then `blastp -q QUERY -d db
-p1 --masking 0 --comp-based-stats MODE --max-hsps MAXIMUM`, adding
`--more-sensitive` for `scaled`. Output uses `--outfmt 100` for DAA or the
fields in `tests/native_cbs.rs` for tabular snapshots. Archive headers retain
the pinned executable's zero values for the unused reward/penalty slots,
as documented in `../seg/README.md`.
