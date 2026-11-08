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

First check that the recorded files still exist and that the inventory includes
all implementation extensions:

```text
python3 scripts/audit_translation_inventory.py
```

This developer-only check requires the adjacent upstream checkout. It returns
nonzero for missing modules or incomplete inventory/audit entries; it does not
prove function-level correspondence or runtime parity and is not part of
packaged Cargo builds.

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


## Command and option registration audit

`command_options.json` inventories every registration in the pinned Config
constructor, including disabled entries and EXTRA/WITH_DNA branches. Refresh
or check it with the developer-only script:

```text
python3 scripts/audit_command_options.py --check translation/command_options.json
python3 scripts/audit_command_options.py --probe .tmp/registry-probes.json
```

The baseline has 21 commands and 299 registered options (161 enabled; 138
hidden/disabled). WITH_DNA adds 17 registrations. EXTRA adds 24 commands,
enables hidden options and removes workflow-membership checks. All 316 names
occur in Rust's registration inventory, but that is name coverage only.
The artifact records storage types/arities, explicit default/min-count
expressions, command groups and final short-option shadowing. Post-store
normalization and actual native execution require separate integration audits.
Only the existing baseline executable is probed; alternate feature variants
are statically modeled and remain runtime conformance gates. Probes abort at
an invalid parser sentinel before executing commands, under resource limits.
