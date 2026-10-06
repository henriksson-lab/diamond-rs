#!/usr/bin/env bash
# Validate downloaded native spill or RNG conformance evidence.
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=lib/portable.sh
source "$script_dir/lib/portable.sh"

usage() {
    echo "usage: $0 spill ARTIFACT_DIR | rng CAPTURE_FILE" >&2
    exit 2
}

if (($# != 2)); then
    usage
fi

mode=$1
artifact=$2

validate_spill() {
    local dir=$1
    local marker="$dir/conformance-passed.txt"
    local summary="$dir/summary.tsv"
    local output_hashes="$dir/output-hashes.sha256"
    local provenance="$dir/provenance.txt"
    local expected actual relative rows manifests summary_output_hash

    for required in "$marker" "$summary" "$output_hashes" "$provenance"; do
        if [[ ! -s "$required" ]]; then
            echo "error: missing spill evidence file: $required" >&2
            exit 1
        fi
    done
    if ! grep -qx 'status=passed' "$marker"; then
        echo "error: spill success marker is absent or invalid" >&2
        exit 1
    fi

    expected=$(sed -n 's/^summary_sha256=//p' "$marker")
    actual=$(sha256_file "$summary")
    if [[ -z "$expected" || "$actual" != "$expected" ]]; then
        echo "error: spill summary hash mismatch" >&2
        exit 1
    fi
    expected=$(sed -n 's/^output_hashes_sha256=//p' "$marker")
    actual=$(sha256_file "$output_hashes")
    if [[ -z "$expected" || "$actual" != "$expected" ]]; then
        echo "error: spill output-manifest hash mismatch" >&2
        exit 1
    fi

    awk -F '\t' '
        NR == 1 {
            if ($0 != "case\tlimit_mib\tclassification\tmigrated_hits\tread_mib\toutput_sha256")
                exit 1
            next
        }
        {
            ++rows
            if (length($6) != 64) exit 1
            if (first_hash == "") first_hash = $6
            if ($6 != first_hash) exit 1
        }
        $1 == "no-spill-a" || $1 == "no-spill-b" {
            if ($3 != "no-spill" || $4 != 0 || $5 != 0) exit 1
            ++no_spill
        }
        $1 == "delayed-gate-1" || $1 == "delayed-gate-2" {
            if ($3 != "delayed-spill" || $4 !~ /^[0-9]+$/ || $4 <= 0 || $5 <= 0.0) exit 1
            ++delayed
        }
        END {
            if (rows == 0 || no_spill != 2 || delayed != 2) exit 1
        }
    ' "$summary" || {
        echo "error: spill summary invariants failed" >&2
        exit 1
    }

    rows=$(awk 'END {print NR - 1}' "$summary")
    summary_output_hash=$(awk -F '\t' 'NR == 2 {print $6}' "$summary")
    manifests=0
    while read -r expected relative; do
        relative=${relative#\*}
        if [[ "$relative" == /* || "$relative" == *..* || "$relative" != runs/*.tsv ]]; then
            echo "error: non-relocatable output path in manifest: $relative" >&2
            exit 1
        fi
        actual=$(sha256_file "$dir/$relative")
        if [[ "$actual" != "$expected" ]]; then
            echo "error: spill output hash mismatch: $relative" >&2
            exit 1
        fi
        if [[ "$expected" != "$summary_output_hash" ]]; then
            echo "error: spill manifest and summary parity hashes differ: $relative" >&2
            exit 1
        fi
        manifests=$((manifests + 1))
    done <"$output_hashes"
    if ((manifests != rows)); then
        echo "error: spill manifest has $manifests outputs for $rows summary rows" >&2
        exit 1
    fi

    echo "validated spill conformance artifact: $dir"
}

validate_rng() {
    local file=$1
    local cpp rust label fields

    if [[ ! -s "$file" ]]; then
        echo "error: missing RNG capture: $file" >&2
        exit 1
    fi
    for label in target-os target-arch target-env; do
        if ! grep -Eq "^${label} [^[:space:]]+\$" "$file"; then
            echo "error: missing RNG target metadata: $label" >&2
            exit 1
        fi
    done
    if ! grep -qx 'status passed' "$file"; then
        echo "error: RNG success marker is absent or invalid" >&2
        exit 1
    fi
    for label in rand-1 rand-12345 rust-rand-1 rust-rand-12345; do
        fields=$(awk -v label="$label" '$1 == label {print NF}' "$file")
        if [[ "$fields" != 101 ]]; then
            echo "error: RNG vector $label has ${fields:-0} fields, expected 101" >&2
            exit 1
        fi
    done
    for label in 1 12345; do
        cpp=$(awk -v label="rand-$label" '$1 == label {$1 = ""; sub(/^ /, ""); print}' "$file")
        rust=$(awk -v label="rust-rand-$label" '$1 == label {$1 = ""; sub(/^ /, ""); print}' "$file")
        if [[ "$cpp" != "$rust" ]]; then
            echo "error: C++ and Rust RNG vectors differ for seed $label" >&2
            exit 1
        fi
    done
    for specification in 'standard-normal 17' 'tantan 42' 'evaluator 41'; do
        label=${specification% *}
        expected=${specification##* }
        fields=$(awk -v label="$label" '$1 == label {print NF}' "$file")
        if [[ "$fields" != "$expected" ]]; then
            echo "error: RNG downstream vector $label has ${fields:-0} fields, expected $expected" >&2
            exit 1
        fi
    done

    echo "validated RNG conformance artifact: $file"
}

case "$mode" in
    spill) validate_spill "$artifact" ;;
    rng) validate_rng "$artifact" ;;
    *) usage ;;
esac
