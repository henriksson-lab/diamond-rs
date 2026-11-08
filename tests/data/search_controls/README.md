# Native approximate identity and soft masking

Oracle: unchanged DIAMOND 2.1.24 build 178, source commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`, executable SHA-256
`d8d08608b96706db62e44036fc979f9a658e72f59fc19d8727546daa616c42d8`.
The original tree and executable were not edited or rebuilt.

`cases.json` records exact source invocations and output hashes; `cases.tsv`
drives 68 normally completing source cases at one/four threads and
in-memory/forced-disk hit buffers: 136 direct CLI comparisons. Sixteen supplied
protein/translated-query API comparisons use the same source snapshots.
Nine parser-error checks and a typed impossible-buffer check are separate.
The tests build reference databases natively and need no C++ at runtime.

Coverage includes approximate thresholds around the source 50/90 Hamming
cutoff boundaries, actual approximate identity filtering and thresholds >100;
score-only fields requiring filter-driven coordinates; soft none/SEG/TANTAN;
automatic motif combinations; hard masking; double/query-indexed algorithms;
direct/multiple query/reference blocks; detailed tabular, score-only, DAA and
binary edge output; repeated options and malformed values. BLOSUM45 soft TANTAN
cases use score-only fields: source disallows approx_pident output with that
matrix but allows approximate identity filtering itself.

TANTAN changes complete emitted outputs for both query types and both
algorithms: baseline 000/012/024/036 differ from 004/016/028/040. Threshold 101
score-only cases are actually empty. These distinguish connected behavior
from no-op option acceptance. Existing cluster/reassign overrides retain
precedence while ordinary search reads the typed configuration.

Combined automatic motif and explicit SEG/TANTAN flags flow to the translated
Masking object built from the actual score matrix. Installed tables reach
hashed query/reference Blocks and ordinary enumeration. Translated contexts
are masked in the authoritative sequence arena, restored before extension,
and prepared tables are reused across shapes/chunks rather than accumulated.
Explicit `none` retains automatic motif masking. The source accepted strings
are exactly 0, none, seg and tantan; other spellings are rejected.

## Source SEG undefined behavior

Six original translated query-indexed SEG cases crash; their exact arguments,
statuses and diagnostics are preserved in `source-failures.json`. They are
excluded from direct parity. Of 18 variant probes (larger query block, algo0, softnone for each original),
17 complete and one larger-block variant also crashes;
`source-compatibility.json` records every status.
Twelve additional native checks compare the original tiny-block configuration
against source algo0 snapshots that agree byte-for-byte with independently
successful softnone snapshots; see `compatibility-selection.json`. The larger
query-block snapshots also change approximate-identity fields in some cases,
so they are not used as equivalent references. Those checks establish
a defined compatibility choice, not byte parity with a crashing process.

The bounded debugger investigation caught original query byte 0x8a (residue 10
plus SEED_MASK 0x80) being sign-extended to -118 before indexing SEG's 128-entry
alphaflag/alphaindex tables. The invalid lookup then produced index 16843009
for a 20-entry composition array. Repeated soft masking of a prepared but empty
SEG table exposes this boundary during later seed enumeration. Native code
tracks prepared-empty tables explicitly and strips seed-mask bits before
SEG alphabet conversion, preserving valid sequence data and memory safety.

Debugger commands/output, numbered original source excerpts and equal
before/after source/binary hashes are retained under
`.tmp/native-search-controls-20261107/source-seg-crash/`. No original patch,
source rebuild or core dump was used. This fixture matrix is bounded evidence;
other advanced modes and production FFI routes remain unfinished.
