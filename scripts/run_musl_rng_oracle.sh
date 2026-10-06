#!/usr/bin/env bash
# Build and run the pinned RNG conformance oracle against musl. This keeps the
# target C++ probe and Rust test on the same C runtime.
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
if [[ ${SKIP_MUSL_BUILD:-0} == 1 ]]; then
    default_target_dir="$repo_dir/target"
else
    default_target_dir="$repo_dir/.tmp/zig-target"
fi
target_dir=${CARGO_TARGET_DIR:-"$default_target_dir"}
cache_root=${ZIG_CACHE_ROOT:-"$repo_dir/.tmp/zig-cache"}
capture=${DIAMOND_RNG_ORACLE_CAPTURE:-"$repo_dir/.tmp/rng-oracle-linux-musl.txt"}
zig_bin=${ZIG:-zig}

mkdir -p -- "$target_dir" "$cache_root/global" "$cache_root/local" "$cache_root/tmp" "$(dirname -- "$capture")"

export CARGO_TARGET_DIR="$target_dir"
export CARGO_ZIGBUILD_CACHE_DIR="$cache_root/cargo-zigbuild"
export CARGO_ZIGBUILD_ZIG_PATH="$zig_bin"
export ZIG="$zig_bin"
export ZIG_GLOBAL_CACHE_DIR="$cache_root/global"
export ZIG_LOCAL_CACHE_DIR="$cache_root/local"
export TMPDIR="$cache_root/tmp"

if [[ ${SKIP_MUSL_BUILD:-0} != 1 ]]; then
    cargo zigbuild --offline --target x86_64-unknown-linux-musl --test compat_rng_oracle
fi

test_dir="$target_dir/x86_64-unknown-linux-musl/debug/deps"
test_bin=$(find "$test_dir" -maxdepth 1 -type f -name 'compat_rng_oracle-*' -perm -111 \
    -printf '%T@\t%p\n' | sort -nr | head -n 1 | cut -f 2-)
if [[ -z "$test_bin" ]]; then
    echo "error: musl RNG oracle test binary not found under $test_dir" >&2
    exit 1
fi

CXX="$repo_dir/scripts/zig-cxx-musl.sh" \
DIAMOND_RNG_ORACLE_CAPTURE="$capture" \
"$test_bin" --ignored --exact target_c_runtime_and_downstream_paths_match_pinned_cpp \
    --test-threads=1 --nocapture

sha256sum "$capture"
