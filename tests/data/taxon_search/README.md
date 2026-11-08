# Native taxonomy search fixtures

Oracle: unchanged DIAMOND 2.1.24 build 178, source commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, executable SHA-256
`d8d08608b96706db62e44036fc979f9a658e72f59fc19d8727546daa616c42d8`.
The C++ tree and executable were not edited or rebuilt.

The 124 successful cases in `cases.tsv` cover protein/translated queries,
both indexing algorithms, direct/multiple query/reference blocks, default/0/1
unaligned settings, optional lineage, mapped/unmapped/multiple-taxid subjects,
leaf/root LCA, absent scientific names when unused, target/score/E-value/
coverage filtering and one/unlimited HSPs. They reuse query inputs from
`../tabular_context/` and full database inputs from `../search_blocks/taxonomy/`.
Twelve more cases use `weak.faa` to distinguish the classifier's default 10%
top-score range from explicit `--top 100` and nonzero `--min-score`, which
suppresses the default. Each runs with one thread/in-memory hit buffers and
four threads/forced-disk hit buffers: 272 exact CLI output comparisons. Three
cases also repeat through the public in-process dispatcher. Twelve supplied
protein/translated-query API comparisons check the same default/overrides.
Thus 284 comparisons are byte exact against normally completing source runs.
Ten additional checks cover missing metadata and malformed/incompatible options.
All runtime tests build their own databases natively; C++ is not required.

The concrete TaxonFormat accumulator is created once per query, receives real
subject taxid lists/nodes/scientific names and HSP E-values, and emits the source
LCA/minimum-mapped-E-value classification. DP requests NONE unless filtering
needs further values. Names are required only with `--include-lineage`; no
rank metadata is requested. The native configuration normalizer applies the
classifier's implicit `--top 10` once, normalizes explicit `--top 100` to
unlimited targets and retains explicit bit-score behavior for CLI and APIs.
The ranking chunk size uses the existing translated overflow-safe
`make_multiple`, as source does, including the unlimited-target sentinel.

## Source undefined behavior and the native compatibility decision

Four original protein cases (equal 50% query/subject coverage, with/without
lineage, direct/multiple blocks) do not have a defined source result. The pinned
executable aborted with allocator errors/signals; repeats also produced a
vector-length exception or segmentation fault. `source-aborts.json` retains
the original commands/statuses/diagnostics. They are not counted as passing
byte parity against those original invocations.

Source `run/double_indexed.cpp` leaves reference titles unloaded when format
102 lacks SSEQID. Equal coverage enables reference length sorting.
`Block::length_sorted` unconditionally indexes/copies `ids_` despite that array
being empty. This static source evidence and the title-loading probe identify
the invalid-access boundary. Adding `--no-self-hits` forces title loading and
avoids it. The fixture's query/subject names and sequences are distinct, so
this flag does not exclude records here. `source-safe-coverage.json` retains
those successful alternative commands and their frozen `.safe.tsv` outputs.

Rust retains valid reference metadata through sorting and performs the defined
filter/LCA calculation. Eight native comparisons (four original cases, one/four
threads) check that result against the safe source variants. These are explicit
source-UB compatibility checks, separate from the 284 direct parity comparisons.
The native result is not an emulation of allocator corruption. This exception
is specific to the recorded fixtures; it is not universal proof that adding
`--no-self-hits` leaves arbitrary searches unchanged. A debugger trace could
not be collected because the environment disallowed ptrace; no C++ patch or
instrumented rebuild was used.

`cases.json`, `weak-source.json` and `errors-source.json` retain source evidence.
`provenance.json` records fixture/input/source-file hashes. Whole-program CLI
error priority, every numerical input, optional upstream build features and
large-scale performance remain outside this bounded matrix.
