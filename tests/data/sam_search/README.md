# Native SAM search fixtures

Frozen from the unchanged DIAMOND 2.1.24 build 178 executable, source commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, executable SHA-256
`d8d08608b96706db62e44036fc979f9a658e72f59fc19d8727546daa616c42d8`.
The C++ tree and executable were not edited or rebuilt.

The 84 cases in `cases.tsv` use the query and three-reference fixtures in
`../tabular_context/`. `cases.json` preserves complete source arguments,
status and file hashes. They cover protein/translated search, double/query
indexing, one/multiple reference blocks, default/0/1/2 unaligned settings,
full subject titles, query-length tags, numeric/text format names, protein
CBS modes 0–5, translated CBS 0/1, FASTA/FASTQ, and TANTAN/SEG masking.
Each case runs with one thread/in-memory hit buffers and four threads/forced
disk hit buffers: 168 comparisons. Four further checks compare malformed
boolean-option diagnostics against `invalid-flags-source.json`.

The complete SAM byte stream is compared after replacing only the source
header's invocation text with the exact native invocation expected from
argv. This field necessarily records different executable, input, database,
output paths, thread counts and memory settings in each run. Every other
header and record byte is unchanged. Three cases also execute the public
in-process `diamond::run` dispatcher. A separate supplied translated-query
block API test compares case 041 at one/four threads (two more comparisons).
All runtime tests are self-contained and build databases natively.

Compatibility details retained from source:

- SAM contains aligned amino-acid sequence even for translated search, with
  CIGAR/MD, oriented source coordinates, frame and score tags from actual HSPs.
- For `--salltitles`, a direct reference block retains its full loaded titles.
  The intermediate block dictionary retains sequence IDs because the SAM
  format lacks FULL_TITLES. Tiny reference blocks therefore produce shorter
  SAM reference names than the single-block case.
- Unaligned output follows the format default and source block/seed behavior;
  the fixture set does not imply that every seedless query is reported in
  every upstream block arrangement.
- `--sam-query-len` and `--salltitles` are boolean switches. Both appended values
  and equals syntax are rejected as in the source. Full-title handling for
  other search formats remains separate integration work.

This is bounded search-format evidence, not proof of whole-program parity.
SAM frameshift search and other advanced search controls remain open.
