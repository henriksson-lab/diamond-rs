#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
rust_bin=${RUST_BIN:-"$repo_dir/target/release/diamond"}
cpp_bin=${CPP_BIN:-"$repo_dir/diamond/build/diamond"}
reference=${REFERENCE_FASTA:-"$repo_dir/diamond/src/test/data.faa"}
query=${QUERY_FASTA:-"$repo_dir/diamond/src/test/5.faa"}
repetitions=${REPETITIONS:-3}
threads=${THREADS:-1}
keep_work=${KEEP_WORK:-0}

if [[ ! -x /usr/bin/time ]]; then
    echo "error: /usr/bin/time is required for peak-RSS measurement" >&2
    exit 2
fi
if [[ ! -s "$reference" || ! -s "$query" ]]; then
    echo "error: reference and query FASTA files must be nonempty" >&2
    exit 2
fi
if [[ ! "$repetitions" =~ ^[1-9][0-9]*$ || ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: REPETITIONS and THREADS must be positive integers" >&2
    exit 2
fi

if [[ ${SKIP_BUILD:-0} != 1 ]]; then
    echo "Building Rust release binary..." >&2
    benchmark_rustflags=${RUSTFLAGS:-}
    if [[ -z "$benchmark_rustflags" ]]; then
        benchmark_arch=$(uname -m)
        if [[ "$benchmark_arch" == x86_64 ]]; then
            # Keep LLVM from introducing AVX-512 into the explicitly AVX2
            # kernels on AVX-512 hosts. That code is slower for this workload.
            benchmark_rustflags='-C target-cpu=x86-64-v3'
        else
            benchmark_rustflags='-C target-cpu=native'
        fi
    fi
    RUSTFLAGS="$benchmark_rustflags" cargo build \
        --manifest-path "$repo_dir/Cargo.toml" --release --offline
    if [[ ! -x "$cpp_bin" ]]; then
        echo "Building vendored C++ release binary..." >&2
        cmake -S "$repo_dir/diamond" -B "$repo_dir/diamond/build" \
            -DCMAKE_BUILD_TYPE=Release
        cmake --build "$repo_dir/diamond/build" -j "$threads"
    fi
fi

if [[ ! -x "$rust_bin" || ! -x "$cpp_bin" ]]; then
    echo "error: Rust or C++ DIAMOND binary is missing" >&2
    exit 2
fi

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/diamond-compare.XXXXXX")
cleanup() {
    if [[ "$keep_work" == 1 ]]; then
        echo "Kept benchmark artifacts in $work_dir" >&2
    else
        rm -rf -- "$work_dir"
    fi
}
trap cleanup EXIT

fasta_stats() {
    awk '
        /^>/ { ++records; next }
        { gsub(/[[:space:]]/, ""); residues += length }
        END { printf "%d\t%d\n", records + 0, residues + 0 }
    ' "$1"
}

IFS=$'\t' read -r reference_records reference_residues < <(fasta_stats "$reference")
IFS=$'\t' read -r query_records query_residues < <(fasta_stats "$query")
printf 'Workload: reference=%s records/%s residues; query=%s records/%s residues; threads=%s; runs=%s\n' \
    "$reference_records" "$reference_residues" "$query_records" "$query_residues" \
    "$threads" "$repetitions"
printf 'Reference SHA-256: '
sha256sum "$reference" | awk '{ print $1 }'
printf 'Query SHA-256: '
sha256sum "$query" | awk '{ print $1 }'

metrics="$work_dir/metrics.tsv"
printf 'implementation\toperation\trun\tseconds\tpeak_rss_kib\n' > "$metrics"

measure() {
    local implementation=$1
    local operation=$2
    local run=$3
    shift 3
    local measurement="$work_dir/time.txt"
    /usr/bin/time -f '%e\t%M' -o "$measurement" "$@" \
        >"$work_dir/${implementation}-${operation}-${run}.stdout" \
        2>"$work_dir/${implementation}-${operation}-${run}.stderr"
    local seconds rss
    IFS=$'\t' read -r seconds rss < "$measurement"
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$implementation" "$operation" "$run" "$seconds" "$rss" >> "$metrics"
}

cpp_db="$work_dir/cpp-db"
rust_db="$work_dir/rust-db"
cpp_output="$work_dir/cpp.tsv"
rust_output="$work_dir/rust.tsv"
outfmt=(6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)

for run in $(seq 1 "$repetitions"); do
    rm -f -- "$cpp_db.dmnd" "$rust_db.dmnd"
    measure cpp makedb "$run" "$cpp_bin" makedb --in "$reference" --db "$cpp_db" \
        --threads "$threads"
    measure rust makedb "$run" "$rust_bin" makedb --in "$reference" --db "$rust_db" \
        --threads "$threads"
done

# Build fresh databases outside the blastp measurements.
"$cpp_bin" makedb --in "$reference" --db "$cpp_db" --threads "$threads" \
    >"$work_dir/cpp-makedb.stdout" 2>"$work_dir/cpp-makedb.stderr"
"$rust_bin" makedb --in "$reference" --db "$rust_db" --threads "$threads" \
    >"$work_dir/rust-makedb.stdout" 2>"$work_dir/rust-makedb.stderr"

for run in $(seq 1 "$repetitions"); do
    if ((run % 2)); then
        implementations=(cpp rust)
    else
        implementations=(rust cpp)
    fi
    for implementation in "${implementations[@]}"; do
        if [[ "$implementation" == cpp ]]; then
            measure cpp blastp "$run" "$cpp_bin" blastp -q "$query" -d "$cpp_db" \
                -o "$cpp_output" --threads "$threads" --outfmt "${outfmt[@]}"
        else
            measure rust blastp "$run" "$rust_bin" blastp -q "$query" -d "$rust_db" \
                -o "$rust_output" --threads "$threads" --outfmt "${outfmt[@]}"
        fi
    done
done

parity=FAIL
if cmp -s -- "$cpp_output" "$rust_output"; then
    parity=PASS
else
    diff -u -- "$cpp_output" "$rust_output" > "$work_dir/parity.diff" || true
fi

median() {
    local implementation=$1
    local operation=$2
    local column=$3
    awk -F '\t' -v impl="$implementation" -v op="$operation" -v col="$column" \
        '$1 == impl && $2 == op { print $col }' "$metrics" \
        | sort -n \
        | awk '{ values[NR] = $1 } END { if (NR % 2) print values[(NR + 1) / 2]; else print (values[NR / 2] + values[NR / 2 + 1]) / 2 }'
}

spread() {
    local implementation=$1
    local operation=$2
    awk -F '\t' -v impl="$implementation" -v op="$operation" '
        $1 == impl && $2 == op {
            ++n
            sum += $4
            sumsq += $4 * $4
            if (n == 1 || $4 < min) min = $4
            if (n == 1 || $4 > max) max = $4
        }
        END {
            mean = sum / n
            variance = sumsq / n - mean * mean
            if (variance < 0) variance = 0
            printf "mean %.3f s, SD %.3f s, range %.2f-%.2f s", mean, sqrt(variance), min, max
        }
    ' "$metrics"
}

printf '\n%-8s %-8s %12s %14s\n' implementation operation seconds peak_rss_kib
for operation in makedb blastp; do
    for implementation in cpp rust; do
        printf '%-8s %-8s %12s %14s\n' \
            "$implementation" "$operation" \
            "$(median "$implementation" "$operation" 4)" \
            "$(median "$implementation" "$operation" 5)"
    done
done

cpp_seconds=$(median cpp blastp 4)
rust_seconds=$(median rust blastp 4)
cpp_rss=$(median cpp blastp 5)
rust_rss=$(median rust blastp 5)
awk -v cpp="$cpp_seconds" -v rust="$rust_seconds" \
    'BEGIN { printf "\nblastp speedup (C++/Rust): %.3fx\n", cpp / rust }'
awk -v cpp="$cpp_rss" -v rust="$rust_rss" \
    'BEGIN { printf "blastp RSS ratio (Rust/C++): %.3fx\n", rust / cpp }'
printf 'C++ blastp timing: %s\n' "$(spread cpp blastp)"
printf 'Rust blastp timing: %s\n' "$(spread rust blastp)"
printf 'blastp byte parity: %s\n' "$parity"
printf 'C++ SHA-256: '
sha256sum "$cpp_output" | awk '{print $1}'
printf 'Rust SHA-256: '
sha256sum "$rust_output" | awk '{print $1}'
if [[ "$keep_work" == 1 || "$parity" != PASS ]]; then
    printf 'Raw metrics: %s\n' "$metrics"
else
    echo 'Set KEEP_WORK=1 to retain raw metrics and outputs.'
fi

if [[ "$parity" != PASS ]]; then
    echo "parity diff: $work_dir/parity.diff" >&2
    keep_work=1
    exit 1
fi
