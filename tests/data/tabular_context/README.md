# Native tabular and flat JSON context fixtures

Oracle: unchanged DIAMOND 2.1.24 build 178, source commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, executable SHA-256
`d8d08608b96706db62e44036fc979f9a658e72f59fc19d8727546daa616c42d8`.
The executable and source tree were neither edited nor rebuilt for these cases.
Source result manifests retain arguments, statuses and output hashes. Runtime
Rust tests build their own databases and need no C++ checkout or executable.

`native_tabular_context` contains five integration tests. The 66 TSV cases in
`cases.tsv`, four dynamic-rank cases, 22 flat JSON cases in `json-cases.tsv`,
and 16 masking/algorithm cases each run with one thread/in-memory hit buffers
and four threads/forced-disk hit buffers: 216 output comparisons. Ten further
protocol/error checks cover short quality, malformed FASTQ, invalid fields,
and unsupported unaligned callbacks. These are bounded fixtures, not proof
of parity for all inputs or CLI error priorities.

The field registry coverage is all 52 standard statically registered source
fields and all 48 dynamically registered taxonomy rank fields. Optional
DP_STAT fields are outside this oracle build. Full query/reference sequence,
translated DNA and amino acids, gapped transcripts, identifiers, full titles,
mate data, quality, taxonomy and self-alignment scores flow from actual search
contexts through query/reference block merging. FASTQ tests include multiline
quality; source rejects multiline sequence records. CBS modes 0–5 are covered
for proteins, and 0/1 for translated queries. Tiny blocks exercise merging;
`--unal 1` cases retain the source's block-dependent handling of seedless queries.

The masked-reference fixtures distinguish original display sequences from
masked extension sequences. Self-alignment scores are captured at the source
phase: after eager reference masking in double-indexed search, before lazy
masking in query-indexed search. Their normalized scores can differ between
algorithms, so each has its own frozen oracle output.

Compatibility decisions:

- Flat JSON compares raw bytes. Source omits separators between queries and
  between multiple HSPs of the first target, leaves string content unescaped,
  and has inconsistent taxonomy quoting/separators. These behaviors are
  deliberately preserved; these fixtures are not generally valid JSON.
- The short-quality source fixture requests only `full_qqual`, for which the
  source output is defined. Slicing `qqual` beyond the quality allocation in
  C++ is undefined; Rust explicitly clamps both slice bounds. Unit tests
  check this safe behavior without claiming byte parity for undefined reads.
- Four unaligned callback errors preserve the source's empty small output
  file by discarding unflushed buffered bytes on failure. Partial-file bytes
  after large failed writes remain unverified.
- Invalid-field diagnostics use valid query/database arguments. General
  missing-input diagnostic priority remains outside this fixture matrix.

`provenance.json` records oracle identity and hashes of all included evidence
and inputs. Historical exploratory failures remain in `.tmp/` separately from
final verification logs; they are not passing evidence.
