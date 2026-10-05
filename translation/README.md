# Translation audit

This directory records the source-to-source audit against the local DIAMOND
baseline at commit `1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`.

## Conventions

- Work bottom-up from the call graph and finish one original source file per
  batch.
- Mirror the original directory hierarchy below `src/` when adding modules.
- Translate C++ function and field names to Rust `snake_case`; record names
  that CCC cannot infer in `ccc_mapping.toml`.
- A file is `translated` only after every function has an explicit Rust
  counterpart and focused tests pass. `parity-tested` additionally requires a
  direct comparison with the C++ implementation or its serialized output.
- Constructors, destructors, operators, templates, and RAII may map to Rust
  traits or ownership rather than one syntactic function. Record these in the
  file map instead of adding artificial functions merely to satisfy matching.

## Refreshing the CCC audit

```text
ccc-rs analyze src -l rust --recurse -o /tmp/diamond-rust.json
ccc-rs analyze diamond/src -l cpp --recurse -o /tmp/diamond-cpp.json
ccc-rs order /tmp/diamond-cpp.json --strict -o /tmp/diamond-order.csv
ccc-rs order-annotate /tmp/diamond-order.csv \
  --source /tmp/diamond-cpp.json \
  --rust /tmp/diamond-rust.json \
  --mapping translation/ccc_mapping.toml \
  -o /tmp/diamond-order-annotated.csv
ccc-rs missing /tmp/diamond-rust.json /tmp/diamond-cpp.json \
  --mapping translation/ccc_mapping.toml
```

CCC output is a guide, not completion proof. Name-only matching is ambiguous
for overloaded functions and common Rust names such as `new`, `from`, and
`drop`; validate those against the source file and tests.

