#!/usr/bin/env bash
# Target-native conformance gate for the RSS-controlled delayed-spill path.
#
# The numeric RSS boundary is deliberately discovered on each runner. Process
# RSS includes allocator, loader, and target-runtime state, so a limit copied
# from another operating system (or even another runner image) is not useful
# conformance evidence.
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=lib/portable.sh
source "$root/scripts/lib/portable.sh"
artifact_dir=${1:-"$root/target/native-spill-conformance"}
source_url=${SPILL_SOURCE_URL:-https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/2025-12-03.1/AMRProt.fa}
source_sha256=${SPILL_SOURCE_SHA256:-b99bfcc6eda5c74a094cc37a5d26d5d06a2dd57037bd00f493efe465dfd84bae}
reference_sha256=abfd1896a97cf8dd3c6fba33965f0be9dc3ce6d1e7184e8109e239ce53e2563b
query_sha256=09f9b80a050e5cc4df32e073a8e485cb288238c084c9970d85ead1425ddcc889
threads=${SPILL_THREADS:-4}
max_limit_mib=${SPILL_MAX_LIMIT_MIB:-2048}

if [[ -n "${DIAMOND_BIN:-}" ]]; then
    diamond_bin=$DIAMOND_BIN
elif [[ -x "$root/target/release/diamond" ]]; then
    diamond_bin="$root/target/release/diamond"
elif [[ -x "$root/target/release/diamond.exe" ]]; then
    diamond_bin="$root/target/release/diamond.exe"
else
    echo "error: release diamond executable not found; run cargo build --release --bin diamond" >&2
    exit 2
fi

for tool in awk cmp curl grep sed tail; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "error: required tool is unavailable: $tool" >&2
        exit 2
    fi
done
if [[ ! "$threads" =~ ^[1-9][0-9]*$ || ! "$max_limit_mib" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: SPILL_THREADS and SPILL_MAX_LIMIT_MIB must be positive integers" >&2
    exit 2
fi
if [[ -e "$artifact_dir" ]]; then
    shopt -s nullglob dotglob
    artifact_entries=("$artifact_dir"/*)
    shopt -u nullglob dotglob
    if ((${#artifact_entries[@]} != 0)); then
        echo "error: artifact directory must be absent or empty: $artifact_dir" >&2
        exit 2
    fi
fi

mkdir -p "$artifact_dir/fixture" "$artifact_dir/runs" "$artifact_dir/spill-tmp"
source_fasta="$artifact_dir/fixture/AMRProt-2025-12-03.1.fa"
if [[ -n "${SPILL_SOURCE_FASTA:-}" ]]; then
    cp "$SPILL_SOURCE_FASTA" "$source_fasta"
else
    curl --fail --location --retry 3 --retry-all-errors \
        --output "$source_fasta" "$source_url"
fi

actual_source_sha=$(sha256_file "$source_fasta")
if [[ "$actual_source_sha" != "$source_sha256" ]]; then
    echo "error: source FASTA checksum mismatch: expected $source_sha256, got $actual_source_sha" >&2
    exit 1
fi

REFERENCE_COUNT=4000 QUERY_COUNT=2000 \
    "$root/scripts/prepare_real_benchmark.sh" "$source_fasta" "$artifact_dir/fixture" \
    >"$artifact_dir/fixture/selection.log"
reference="$artifact_dir/fixture/reference-4000.faa"
query="$artifact_dir/fixture/query-2000.faa"
actual_reference_sha=$(sha256_file "$reference")
actual_query_sha=$(sha256_file "$query")
if [[ "$actual_reference_sha" != "$reference_sha256" || "$actual_query_sha" != "$query_sha256" ]]; then
    echo "error: deterministic fixture selection differs from the pinned fixture" >&2
    echo "reference expected=$reference_sha256 actual=$actual_reference_sha" >&2
    echo "query expected=$query_sha256 actual=$actual_query_sha" >&2
    exit 1
fi

db="$artifact_dir/fixture/amr-4000.dmnd"
"$diamond_bin" makedb --in "$reference" --db "$db" --threads "$threads" \
    >"$artifact_dir/makedb.stdout" 2>"$artifact_dir/makedb.stderr"

{
    printf 'source_url=%s\n' "$source_url"
    printf 'source_sha256=%s\n' "$actual_source_sha"
    printf 'reference_sha256=%s\n' "$actual_reference_sha"
    printf 'query_sha256=%s\n' "$actual_query_sha"
    printf 'diamond_sha256=%s\n' "$(sha256_file "$diamond_bin")"
    printf 'threads=%s\n' "$threads"
    printf 'repository_revision=%s\n' "${GITHUB_SHA:-$(git -C "$root" rev-parse HEAD 2>/dev/null || printf unknown)}"
    printf 'github_run_id=%s\n' "${GITHUB_RUN_ID:-local}"
    printf 'github_run_attempt=%s\n' "${GITHUB_RUN_ATTEMPT:-local}"
    printf 'runner_os=%s\n' "${RUNNER_OS:-unknown}"
    printf 'runner_arch=%s\n' "${RUNNER_ARCH:-unknown}"
    printf 'uname='; uname -a || true
    printf 'rustc='; rustc --version || true
} >"$artifact_dir/provenance.txt"
printf 'case\tlimit_mib\tclassification\tmigrated_hits\tread_mib\toutput_sha256\n' \
    >"$artifact_dir/summary.tsv"

baseline_output=""
RUN_CLASS=""
RUN_MIGRATED=""
RUN_READ_MIB=""
run_serial=0

run_case() {
    local label=$1
    local limit_mib=$2
    local output stderr stdout output_sha spill_lines read_lines
    run_serial=$((run_serial + 1))
    output="$artifact_dir/runs/$(printf '%03d' "$run_serial")-${label}-${limit_mib}M.tsv"
    stdout="${output%.tsv}.stdout"
    stderr="${output%.tsv}.stderr"
    {
        printf '%q ' "$diamond_bin" blastp -q "$query" -d "$db" -o "$output" \
            --threads "$threads" --memory-limit "${limit_mib}M" \
            --tmpdir "$artifact_dir/spill-tmp"
        printf '\n'
    } >"${output%.tsv}.command"
    "$diamond_bin" blastp -q "$query" -d "$db" -o "$output" \
        --threads "$threads" --memory-limit "${limit_mib}M" \
        --tmpdir "$artifact_dir/spill-tmp" >"$stdout" 2>"$stderr"
    if [[ ! -s "$output" ]]; then
        echo "error: search produced an empty result: $label at ${limit_mib}M" >&2
        exit 1
    fi
    if [[ -n "$baseline_output" ]] && ! cmp -s "$baseline_output" "$output"; then
        echo "error: byte output differs from no-spill baseline: $label at ${limit_mib}M" >&2
        cmp "$baseline_output" "$output" >"${output%.tsv}.cmp" 2>&1 || true
        exit 1
    fi

    spill_lines=$(grep -c 'Hit buffer: RSS reached --memory-limit; spilling [0-9][0-9]* retained hits' "$stderr" || true)
    read_lines=$(grep -c 'Hit buffer: read [0-9][0-9.]* MiB of compressed spill data' "$stderr" || true)
    if ((spill_lines == 0 && read_lines == 0)); then
        RUN_CLASS=no-spill
        RUN_MIGRATED=0
        RUN_READ_MIB=0
    elif ((spill_lines == 1 && read_lines == 1)); then
        RUN_MIGRATED=$(sed -n 's/.*spilling \([0-9][0-9]*\) retained hits.*/\1/p' "$stderr" | tail -n 1)
        RUN_READ_MIB=$(sed -n 's/.*read \([0-9][0-9.]*\) MiB of compressed spill data.*/\1/p' "$stderr" | tail -n 1)
        if [[ -z "$RUN_MIGRATED" || -z "$RUN_READ_MIB" ]]; then
            echo "error: could not parse spill diagnostics for $label" >&2
            exit 1
        fi
        if ((RUN_MIGRATED > 0)) && awk -v value="$RUN_READ_MIB" 'BEGIN { exit !(value > 0.0) }'; then
            RUN_CLASS=delayed-spill
        else
            RUN_CLASS=startup-spill
        fi
    else
        echo "error: incomplete or repeated spill diagnostics for $label: spill=$spill_lines read=$read_lines" >&2
        exit 1
    fi
    output_sha=$(sha256_file "$output")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$limit_mib" "$RUN_CLASS" "$RUN_MIGRATED" "$RUN_READ_MIB" "$output_sha" \
        >>"$artifact_dir/summary.tsv"
    printf '%-22s limit=%4sM class=%-13s migrated=%s read=%s MiB sha256=%s\n' \
        "$label" "$limit_mib" "$RUN_CLASS" "$RUN_MIGRATED" "$RUN_READ_MIB" "$output_sha"
}

# A limit orders of magnitude beyond this fixture's footprint is the stable
# no-spill side. Run it twice so platform RSS reporting must consistently keep
# the search in memory, and use its bytes as the parity oracle for every probe.
run_case no-spill-a 16384
if [[ "$RUN_CLASS" != no-spill ]]; then
    echo "error: 16 GiB calibration run unexpectedly spilled" >&2
    exit 1
fi
baseline_output="$artifact_dir/runs/001-no-spill-a-16384M.tsv"
run_case no-spill-b 16384
if [[ "$RUN_CLASS" != no-spill ]]; then
    echo "error: repeated 16 GiB run unexpectedly spilled" >&2
    exit 1
fi

# Establish a startup-spill lower side, then double until this target's RSS
# implementation reports a no-spill upper side. Remember the first genuine
# delayed transition encountered. This calibration avoids baking a Linux RSS
# number into macOS or Windows CI.
run_case calibration-low 1
if [[ "$RUN_CLASS" != startup-spill ]]; then
    echo "error: 1 MiB did not exercise the expected immediate spill side" >&2
    exit 1
fi
startup_limit=1
delayed_limit=0
no_spill_limit=0
probe_limit=32
while ((probe_limit <= max_limit_mib)); do
    run_case calibration "$probe_limit"
    case "$RUN_CLASS" in
        startup-spill) startup_limit=$probe_limit ;;
        delayed-spill)
            if ((delayed_limit == 0)); then
                delayed_limit=$probe_limit
            fi
            ;;
        no-spill)
            no_spill_limit=$probe_limit
            break
            ;;
    esac
    probe_limit=$((probe_limit * 2))
done
if ((no_spill_limit == 0)); then
    echo "error: no target-native no-spill upper bound found through ${max_limit_mib}M" >&2
    exit 1
fi

# If doubling stepped over the delayed window, bisect the startup/no-spill
# bracket until an integer-MiB delayed limit is observed.
if ((delayed_limit == 0)); then
    lower=$startup_limit
    upper=$no_spill_limit
    while ((upper - lower > 1)); do
        middle=$(((lower + upper) / 2))
        run_case calibration-bisect "$middle"
        case "$RUN_CLASS" in
            startup-spill) lower=$middle ;;
            delayed-spill) delayed_limit=$middle; break ;;
            no-spill) upper=$middle ;;
        esac
    done
fi
if ((delayed_limit == 0)); then
    echo "error: target-native calibration found no genuine delayed-spill limit" >&2
    exit 1
fi

# Move away from the startup edge by locating the upper end of the observed
# delayed interval. A classification reversal is evidence of an unstable RSS
# environment, so fail instead of silently selecting a convenient threshold.
lower=$delayed_limit
upper=$no_spill_limit
while ((upper - lower > 1)); do
    middle=$(((lower + upper) / 2))
    run_case calibration-upper "$middle"
    case "$RUN_CLASS" in
        delayed-spill) lower=$middle ;;
        no-spill) upper=$middle ;;
        startup-spill)
            echo "error: RSS classification reversed while locating the delayed-spill interval" >&2
            exit 1
            ;;
    esac
done
selected_limit=$(((delayed_limit + lower) / 2))
printf 'selected_limit_mib=%s\n' "$selected_limit" >>"$artifact_dir/provenance.txt"

# Repeat the selected interior point. Both runs must migrate already-retained
# hits, read a nonzero compressed payload, and remain byte-identical to both
# no-spill runs. This is the actual platform conformance gate; calibration
# probes are retained as diagnostic evidence.
for repeat in 1 2; do
    run_case "delayed-gate-${repeat}" "$selected_limit"
    if [[ "$RUN_CLASS" != delayed-spill ]]; then
        echo "error: selected ${selected_limit}M boundary was not a stable delayed spill" >&2
        exit 1
    fi
done

for output in "$artifact_dir"/runs/*.tsv; do
    relative_output="runs/${output##*/}"
    printf '%s  %s\n' "$(sha256_file "$output")" "$relative_output"
done >"$artifact_dir/output-hashes.sha256"
{
    printf 'status=passed\n'
    printf 'summary_sha256=%s\n' "$(sha256_file "$artifact_dir/summary.tsv")"
    printf 'output_hashes_sha256=%s\n' "$(sha256_file "$artifact_dir/output-hashes.sha256")"
} >"$artifact_dir/conformance-passed.txt"
bash "$root/scripts/validate_native_conformance_artifact.sh" spill "$artifact_dir"
echo "native spill conformance passed; artifacts: $artifact_dir"
