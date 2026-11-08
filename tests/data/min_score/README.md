# Bit-score cutoff conformance

Fixtures come from the unchanged DIAMOND 2.1.24 build 178 executable pinned
in `provenance.json`. Inputs are `tests/data/multi_hsp/{reference.faa,
query.faa,query.fna}`. The source database is constructed from that reference.
No C++ checkout is read by the Rust tests.

All cases use `--masking 0 --max-hsps 0 -e 1e-200`. Threshold zero retains
the E-value cutoff and produces no alignments. Nonzero `--min-score` overrides
that cutoff using unrounded bit scores. Thresholds 100, 170, 174.5, 175 and
180 exercise retained, partially retained, and empty results. Protein cases
cover CBS modes 0–5; translated cases cover the permitted modes 0 and 1.
`source-results.json` records commands, output hashes, and row counts.

TSV fields are `qseqid sseqid score bitscore evalue qstart qend sstart send
positive`; `positive` requires transcript computation. Four DAA fixtures also
verify serialized output. Native integration tests additionally project the
TSV fields to exercise score-only and coordinate passes, query/reference
blocking, one/four threads, and in-memory/forced-disk hit buffers. Numeric
conversion cases check source C-style prefixes, hexadecimal floats, NaN,
infinity, negative values, and repeated scalar options through the public
in-process native runner.

`source-field-projections.json` records 96 direct C++ runs with score-only or
coordinate field sets, confirming the projected expectations for these inputs.
Its qualification records removal of a leading terminal-color reset if present;
no TSV values are normalized.
