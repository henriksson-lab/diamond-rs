#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
rust_bin=${RUST_BIN:-"$repo_dir/target/release/diamond"}
cpp_bin=${CPP_BIN:-"$repo_dir/diamond/build/diamond"}
reference=${REFERENCE_FASTA:-"$repo_dir/diamond/src/test/data.faa"}
query=${QUERY_FASTA:-"$repo_dir/diamond/src/test/5.faa"}
cpp_db_override=${CPP_DB:-}
rust_db_override=${RUST_DB:-}
repetitions=${REPETITIONS:-3}
threads=${THREADS:-1}
keep_work=${KEEP_WORK:-0}
benchmark_label=${BENCH_LABEL:-custom}
summary_file=${SUMMARY_FILE:-}
search_command=${COMMAND:-blastp}
search_args_text=${SEARCH_ARGS:-}
outfmt_text=${OUTFMT_FIELDS:-"qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"}
rust_memory_limit=${RUST_MEMORY_LIMIT:-}
benchmark_cpuset=${BENCH_CPUSET:-}
benchmark_prefix=()
search_args=()
outfmt=()

if [[ "$search_command" != blastp && "$search_command" != blastx ]]; then
    echo "error: COMMAND must be blastp or blastx" >&2
    exit 2
fi
if [[ -n "$search_args_text" ]]; then
    # Benchmark flags and values do not contain whitespace. Keeping this as a
    # plain word list avoids eval and makes the exact invocation visible.
    read -r -a search_args <<< "$search_args_text"
fi
read -r -a outfmt <<< "$outfmt_text"
if ((${#outfmt[@]} == 0)); then
    echo "error: OUTFMT_FIELDS must contain at least one tabular field" >&2
    exit 2
fi

if [[ -n "$benchmark_cpuset" ]]; then
    if ! command -v taskset >/dev/null 2>&1; then
        echo "error: BENCH_CPUSET requires taskset" >&2
        exit 2
    fi
    benchmark_prefix=(taskset --cpu-list "$benchmark_cpuset")
fi

if [[ ! -x /usr/bin/time ]]; then
    echo "error: /usr/bin/time is required for peak-RSS measurement" >&2
    exit 2
fi
if [[ -n "$cpp_db_override" && -z "$rust_db_override" ]] || \
    [[ -z "$cpp_db_override" && -n "$rust_db_override" ]]; then
    echo "error: CPP_DB and RUST_DB must either both be set or both be unset" >&2
    exit 2
fi
use_prebuilt_db=0
if [[ -n "$cpp_db_override" ]]; then
    use_prebuilt_db=1
    cpp_db_override=${cpp_db_override%.dmnd}
    rust_db_override=${rust_db_override%.dmnd}
    if [[ ! -s "$cpp_db_override.dmnd" || ! -s "$rust_db_override.dmnd" ]]; then
        echo "error: CPP_DB and RUST_DB must name nonempty DIAMOND databases" >&2
        exit 2
    fi
elif [[ ! -s "$reference" ]]; then
    echo "error: reference FASTA must be nonempty" >&2
    exit 2
fi
if [[ ! -s "$query" ]]; then
    echo "error: query FASTA must be nonempty" >&2
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
        # Match the vendored C++ release build's native-CPU policy. Individual
        # kernels still perform their own runtime feature dispatch.
        benchmark_rustflags='-C target-cpu=native'
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
    local fasta=$1
    local magic
    magic=$(od -An -tx1 -N2 -- "$fasta" | tr -d ' \n')
    local reader=(cat -- "$fasta")
    if [[ "$magic" == 1f8b ]]; then
        reader=(gzip -cd -- "$fasta")
    fi
    "${reader[@]}" | awk '
        /^>/ { ++records; next }
        { gsub(/[[:space:]]/, ""); residues += length }
        END { printf "%d\t%d\n", records + 0, residues + 0 }
    '
}

IFS=$'\t' read -r query_records query_residues < <(fasta_stats "$query")
if ((use_prebuilt_db)); then
    printf 'Workload: prebuilt databases; query=%s records/%s residues; threads=%s; runs=%s\n' \
        "$query_records" "$query_residues" "$threads" "$repetitions"
    printf 'C++ database: %s.dmnd\nRust database: %s.dmnd\n' \
        "$cpp_db_override" "$rust_db_override"
else
    IFS=$'\t' read -r reference_records reference_residues < <(fasta_stats "$reference")
    printf 'Workload: reference=%s records/%s residues; query=%s records/%s residues; threads=%s; runs=%s\n' \
        "$reference_records" "$reference_residues" "$query_records" "$query_residues" \
        "$threads" "$repetitions"
fi
printf 'Search mode: %s' "$search_command"
if ((${#search_args[@]})); then
    printf ' %q' "${search_args[@]}"
fi
printf '\n'
printf 'Output fields:'
printf ' %q' "${outfmt[@]}"
printf '\n'
if [[ -n "$benchmark_cpuset" ]]; then
    printf 'CPU affinity: %s (applied identically to C++ and Rust)\n' "$benchmark_cpuset"
fi
if [[ -n "$rust_memory_limit" ]]; then
    printf 'Rust hit-buffer memory limit: %s (upstream blastp is disk-backed)\n' \
        "$rust_memory_limit"
else
    printf 'Rust hit-buffer memory limit: 16G (native default; adaptive spill)\n'
fi
if ((!use_prebuilt_db)); then
    printf 'Reference SHA-256: '
    sha256sum "$reference" | awk '{ print $1 }'
fi
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
    local stdout="$work_dir/${implementation}-${operation}-${run}.stdout"
    local stderr="$work_dir/${implementation}-${operation}-${run}.stderr"
    if ! /usr/bin/time -f '%e\t%M' -o "$measurement" "${benchmark_prefix[@]}" "$@" \
        >"$stdout" 2>"$stderr"; then
        keep_work=1
        printf 'error: %s %s run %s failed; stderr follows:\n' \
            "$implementation" "$operation" "$run" >&2
        sed -n '1,160p' "$stderr" >&2
        return 1
    fi
    local seconds rss
    IFS=$'\t' read -r seconds rss < "$measurement"
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$implementation" "$operation" "$run" "$seconds" "$rss" >> "$metrics"
}

if ((use_prebuilt_db)); then
    cpp_db=$cpp_db_override
    rust_db=$rust_db_override
else
    cpp_db="$work_dir/cpp-db"
    rust_db="$work_dir/rust-db"
fi
cpp_output=
rust_output=
outfmt=(6 "${outfmt[@]}")
rust_search_extra=()
if [[ -n "$rust_memory_limit" ]]; then
    rust_search_extra+=(--memory-limit "$rust_memory_limit" --tmpdir "$work_dir")
fi

parity=PASS
expected_output_sha256=${EXPECTED_OUTPUT_SHA256:-}
if ((!use_prebuilt_db)); then
    for run in $(seq 1 "$repetitions"); do
        rm -f -- "$cpp_db.dmnd" "$rust_db.dmnd"
        measure cpp makedb "$run" "$cpp_bin" makedb --in "$reference" --db "$cpp_db" \
            --threads "$threads"
        measure rust makedb "$run" "$rust_bin" makedb --in "$reference" --db "$rust_db" \
            --threads "$threads"
    done

    # Build fresh databases outside the search measurements.
    "$cpp_bin" makedb --in "$reference" --db "$cpp_db" --threads "$threads" \
        >"$work_dir/cpp-makedb.stdout" 2>"$work_dir/cpp-makedb.stderr"
    "$rust_bin" makedb --in "$reference" --db "$rust_db" --threads "$threads" \
        >"$work_dir/rust-makedb.stdout" 2>"$work_dir/rust-makedb.stderr"
fi

for run in $(seq 1 "$repetitions"); do
    cpp_output="$work_dir/cpp-${run}.tsv"
    rust_output="$work_dir/rust-${run}.tsv"
    if ((run % 2)); then
        implementations=(cpp rust)
    else
        implementations=(rust cpp)
    fi
    for implementation in "${implementations[@]}"; do
        if [[ "$implementation" == cpp ]]; then
            measure cpp "$search_command" "$run" "$cpp_bin" "$search_command" \
                -q "$query" -d "$cpp_db" -o "$cpp_output" --threads "$threads" \
                "${search_args[@]}" --outfmt "${outfmt[@]}"
        else
            measure rust "$search_command" "$run" "$rust_bin" "$search_command" \
                -q "$query" -d "$rust_db" -o "$rust_output" --threads "$threads" \
                "${rust_search_extra[@]}" "${search_args[@]}" --outfmt "${outfmt[@]}"
        fi
    done
    if ! cmp -s -- "$cpp_output" "$rust_output"; then
        parity=FAIL
        diff -u -- "$cpp_output" "$rust_output" > "$work_dir/parity-run-${run}.diff" || true
    fi
    if [[ -n "$expected_output_sha256" ]]; then
        cpp_sha=$(sha256sum "$cpp_output" | awk '{print $1}')
        rust_sha=$(sha256sum "$rust_output" | awk '{print $1}')
        if [[ "$cpp_sha" != "$expected_output_sha256" || "$rust_sha" != "$expected_output_sha256" ]]; then
            parity=FAIL
            printf 'error: run %s output hash differs from EXPECTED_OUTPUT_SHA256=%s (C++=%s, Rust=%s)\n' \
                "$run" "$expected_output_sha256" "$cpp_sha" "$rust_sha" >&2
        fi
    fi
done

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
operations=("$search_command")
if ((!use_prebuilt_db)); then
    operations=(makedb "$search_command")
fi
for operation in "${operations[@]}"; do
    for implementation in cpp rust; do
        printf '%-8s %-8s %12s %14s\n' \
            "$implementation" "$operation" \
            "$(median "$implementation" "$operation" 4)" \
            "$(median "$implementation" "$operation" 5)"
    done
done

cpp_seconds=$(median cpp "$search_command" 4)
rust_seconds=$(median rust "$search_command" 4)
cpp_rss=$(median cpp "$search_command" 5)
rust_rss=$(median rust "$search_command" 5)
awk -v cpp="$cpp_seconds" -v rust="$rust_seconds" \
    -v command="$search_command" \
    'BEGIN { printf "\n%s speedup (C++/Rust): %.3fx\n", command, cpp / rust }'
awk -v cpp="$cpp_rss" -v rust="$rust_rss" \
    -v command="$search_command" \
    'BEGIN { printf "%s RSS ratio (Rust/C++): %.3fx\n", command, rust / cpp }'
printf 'C++ %s timing: %s\n' "$search_command" "$(spread cpp "$search_command")"
printf 'Rust %s timing: %s\n' "$search_command" "$(spread rust "$search_command")"
printf '%s byte parity across all %s run pairs: %s\n' \
    "$search_command" "$repetitions" "$parity"
printf 'C++ SHA-256: '
sha256sum "$cpp_output" | awk '{print $1}'
printf 'Rust SHA-256: '
sha256sum "$rust_output" | awk '{print $1}'
if [[ -n "$summary_file" ]]; then
    speed_ratio=$(awk -v cpp="$cpp_seconds" -v rust="$rust_seconds" \
        'BEGIN { printf "%.6f", cpp / rust }')
    rss_ratio=$(awk -v cpp="$cpp_rss" -v rust="$rust_rss" \
        'BEGIN { printf "%.6f", rust / cpp }')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$benchmark_label" "$search_command" "$threads" "$repetitions" \
        "${rust_memory_limit:-16G-default}" "$search_args_text" \
        "$cpp_seconds" "$rust_seconds" "$speed_ratio" "$cpp_rss" "$rust_rss" \
        "$rss_ratio" "$parity" >> "$summary_file"
fi
if [[ "$keep_work" == 1 || "$parity" != PASS ]]; then
    printf 'Raw metrics: %s\n' "$metrics"
else
    echo 'Set KEEP_WORK=1 to retain raw metrics and outputs.'
fi

if [[ "$parity" != PASS ]]; then
    echo "parity artifacts: $work_dir/parity-run-*.diff" >&2
    keep_work=1
    exit 1
fi
