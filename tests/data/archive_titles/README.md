# Native archive title and ID policy

Oracle: unchanged DIAMOND 2.1.24 build 178, source commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, executable SHA-256
`d8d08608b96706db62e44036fc979f9a658e72f59fc19d8727546daa616c42d8`.
No original source, object or executable was edited or rebuilt.

`cases.json` records 78 normally completing original invocations and hashes.
The self-contained integration test runs each at one/four threads with
in-memory/forced-disk hit buffers: 156 complete CLI comparisons. DAA outputs
also undergo 48 native archive/view comparisons. Supplied protein records,
protein query files with a decoded reference Block, and supplied translated
Blocks add 24 complete archive comparisons. Four checks cover boolean option
count/equals-syntax errors. Tests create their databases natively and do not
require the C++ tree or adapter at runtime.

The reference FASTA contains descriptions, XML characters, actual byte 0x01
(four occurrences), byte 0x02 and a space-before-greater-than title separator.
Default archive dictionaries retain the first ID; --sallseqid retains all IDs;
--salltitles copies the raw loaded full title including control separators.
Full titles take precedence when both switches are set. Raw archive sizes
and bytes distinguish all three policies for both query types.

The typed ArchiveOptions is forwarded through BlastxConfig/BlastpConfig.
Its helper uses the translated DAA format flags and sequence-title functions.
All dictionary insertion paths share it: per-reference reservation, late query
insertion, and eager final dictionary construction. Ordinary file/supplied-query
tests exercise the blocked adapter; this does not claim separate eager BLAST
input end-to-end coverage, whose source-UB policy remains unfinished.

Other formats accept these switches with source-specific effects. SAM uses
salltitles but ignores sallseqid. Tabular/JSON fields named salltitles/sallseqid
request their own metadata independently of the switches. Pairwise, XML, PAF,
edge and null accept the switches as no-ops. Only the actual invocation header
and XML database path are substituted with each native command's expected
metadata where those values necessarily differ. DAA byte comparisons have
no substitutions.

The input contains unaligned queries: source protein archives retain one of
two queries, translated archives two of three. Explicit --unal 0/1 still writes
aligned-only DAA query records; the native archive writer preserves that
behavior. `source-noops.json` independently records other-format and DAA-unal
byte equality. This closes these native routing gaps, not the full translation.
