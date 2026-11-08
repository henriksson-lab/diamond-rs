# TANTAN LambdaCalculator original-object fixtures

These fixtures audit `diamond/src/lib/tantan/LambdaCalculator.cc`, the calculator
called by the active native TantanMasker. It differs from the unrelated
`diamond/src/masking/lambda.cpp` helper.

The oracle is the original, unchanged GNU/Linux object linked into DIAMOND
2.1.24 build 178, source commit `1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`.
Its object SHA-256 is
`2c85f6fc7952073f71423672284e8087712d6c25f7574e0c977e58c878cee4ca`.
Only a new external driver was compiled and linked against that existing
object. No upstream source, object or executable was edited or rebuilt.
`source-provenance.json` records exact compiler arguments, immutable hashes,
fixture hashes and the driver hash. Derivation diagnostics remain under
`.tmp/native-lambda-20261107/` in the working checkout.

`matrices.tsv` supplies 36 matrix/seed/state-transition cases: all eight
registered standard matrices, asymmetric/transposed and simple DNA matrices,
zero-axis and singular cases, positive/negative/empty invalid states, reseeding
and continuation. `states.tsv` retains bad state, full lambda bits, both complete
probability-vector bit arrays and the next GNU RNG draw. One sequential native
integration test owns every RNG-consuming case; the RNG is process-global.
Reset behavior is checked after every calculation. Runtime tests require no C++.

Baseline Rust differed in seven states. Explicit fused subtract-products in
LU and both triangular solves match the pinned object's instructions and fix
all 36 states in standalone diagnostics. A zero-axis case previously chose a
different valid lambda and consumed a different number of random draws;
a loose tolerance would have concealed that behavioral difference.

`rounding.tsv` retains 19 exact snprintf("%g")/atof probes: halfway neighbors,
exponent carry and notation boundaries, subnormal/1e-300, signed zero, infinity
and NaN. The native helper uses six-significant-digit scientific formatting
followed by decimal parsing, avoiding arithmetic round's different tie rule
and scale overflow. The focused unit test preserves all output bit patterns.

The constructor/accessors/reset/calculate methods map to LambdaCalculator's
Rust methods. `roundToFewDigits` maps to `round_to_few_digits`; `max_index`,
`lu_pivoting`, `calculate_inv_sum`, `binary_search`, `calculate_lambda` and
`check_lambda` retain their names; `find_ub` maps to `find_upper_bound` and
`solvep` to `solve_p`. Matrix allocation/deletion is initialized owned Vec
storage plus RAII. Matrix/vector swaps use slice swaps; transpose is inlined
in `invert`; ordered summation is inside `calculate_inv_sum`. `inverse_at`
factors owned matrix construction/inversion for reuse by these routines.
`calculate_flat_i8` is the safe native score-matrix adapter used by masking.

Rust rejects malformed dimensions and nonfinite probabilities safely, and frees
failed-path allocations. Source assumes valid pointers and leaks some failed
allocations; these are defined safety decisions, not omitted source functions.
The deterministic GNU-compatible RNG remains an explicit portability choice.

Fused rounding is specified on every ISA, but exp/ln still depend on target
math libraries. Exact GNU/Linux vectors are retained on every platform without
skips/tolerance; macOS/Windows bit parity remains a hosted evidence gate.
Standalone derivation is complete; coordinated Cargo and pipeline validation
is recorded separately by the parent batch.
