#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source_fasta=${REAL_PROTEIN_FASTA:-}
fixture_dir=${FIXTURE_DIR:-/tmp/diamond-real-benchmark}
reference_count=${REFERENCE_COUNT:-2000}
query_count=${QUERY_COUNT:-1000}
threads=${THREADS:-1}
stat_repetitions=${STAT_REPETITIONS:-3}
rust_bin=${RUST_BIN:-"$repo_dir/target/release/diamond"}
cpp_bin=${CPP_BIN:-"$repo_dir/diamond/build/diamond"}

if [[ -z "$source_fasta" ]]; then
    echo "error: set REAL_PROTEIN_FASTA to a diverse, uncompressed protein FASTA" >&2
    exit 2
fi
if ! command -v perf >/dev/null 2>&1; then
    echo "error: perf is required" >&2
    exit 2
fi
if [[ ! "$threads" =~ ^[1-9][0-9]*$ || ! "$stat_repetitions" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: THREADS and STAT_REPETITIONS must be positive integers" >&2
    exit 2
fi

REFERENCE_COUNT=$reference_count QUERY_COUNT=$query_count \
    "$repo_dir/scripts/prepare_real_benchmark.sh" "$source_fasta" "$fixture_dir"

if [[ ${SKIP_BUILD:-0} != 1 ]]; then
    benchmark_rustflags=${RUSTFLAGS:-}
    if [[ -z "$benchmark_rustflags" ]]; then
        benchmark_arch=$(uname -m)
        if [[ "$benchmark_arch" == x86_64 ]]; then
            benchmark_rustflags='-C target-cpu=x86-64-v3'
        else
            benchmark_rustflags='-C target-cpu=native'
        fi
    fi
    RUSTFLAGS="$benchmark_rustflags" cargo build \
        --manifest-path "$repo_dir/Cargo.toml" --release --offline
    if [[ ! -x "$cpp_bin" ]]; then
        cmake -S "$repo_dir/diamond" -B "$repo_dir/diamond/build" \
            -DCMAKE_BUILD_TYPE=Release
        cmake --build "$repo_dir/diamond/build" -j "$threads"
    fi
fi

if [[ ! -x "$rust_bin" || ! -x "$cpp_bin" ]]; then
    echo "error: Rust or C++ DIAMOND binary is missing" >&2
    exit 2
fi

profile_root=${PROFILE_DIR:-"$fixture_dir/profiles"}
mkdir -p -- "$profile_root"
work_dir=$(mktemp -d "$profile_root/run.XXXXXX")
reference="$fixture_dir/reference-${reference_count}.faa"
query="$fixture_dir/query-${query_count}.faa"
rust_db="$work_dir/rust-db"
cpp_db="$work_dir/cpp-db"
outfmt=(6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)

"$rust_bin" makedb --in "$reference" --db "$rust_db" --threads "$threads" \
    >"$work_dir/rust-makedb.stdout" 2>"$work_dir/rust-makedb.stderr"
"$cpp_bin" makedb --in "$reference" --db "$cpp_db" --threads "$threads" \
    >"$work_dir/cpp-makedb.stdout" 2>"$work_dir/cpp-makedb.stderr"

profile_one() {
    local implementation=$1
    local binary=$2
    local database=$3
    local output="$work_dir/$implementation.tsv"
    local command=("$binary" blastp -q "$query" -d "$database" -o "$output"
        --threads "$threads" --quiet --outfmt "${outfmt[@]}")

    perf stat -r "$stat_repetitions" \
        -e task-clock,cycles,instructions,cache-references,cache-misses,branches,branch-misses \
        -o "$work_dir/$implementation-stat.txt" -- "${command[@]}"
    perf record -F 997 -g --call-graph dwarf,16384 \
        -o "$work_dir/$implementation-perf.data" -- "${command[@]}"
    perf report --stdio --no-children --call-graph none --sort symbol \
        --percent-limit 0.5 -i "$work_dir/$implementation-perf.data" \
        >"$work_dir/$implementation-report.txt"
}

profile_one rust "$rust_bin" "$rust_db"
profile_one cpp "$cpp_bin" "$cpp_db"

if ! cmp -s -- "$work_dir/rust.tsv" "$work_dir/cpp.tsv"; then
    diff -u -- "$work_dir/cpp.tsv" "$work_dir/rust.tsv" >"$work_dir/parity.diff" || true
    echo "error: output parity failed; see $work_dir/parity.diff" >&2
    exit 1
fi

echo "Profile artifacts: $work_dir"
echo "Output parity: PASS"
echo "Rust counters: $work_dir/rust-stat.txt"
echo "Rust flat profile: $work_dir/rust-report.txt"
echo "C++ counters: $work_dir/cpp-stat.txt"
echo "C++ flat profile: $work_dir/cpp-report.txt"
