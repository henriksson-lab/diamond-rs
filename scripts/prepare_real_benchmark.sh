#!/usr/bin/env bash
set -euo pipefail

if (($# < 1 || $# > 2)); then
    echo "usage: $0 SOURCE_FASTA [OUTPUT_DIR]" >&2
    exit 2
fi

source_fasta=$1
output_dir=${2:-/tmp/diamond-real-benchmark}
reference_count=${REFERENCE_COUNT:-2000}
query_count=${QUERY_COUNT:-300}

if [[ ! -s "$source_fasta" ]]; then
    echo "error: source FASTA is missing or empty: $source_fasta" >&2
    exit 2
fi
if [[ ! "$reference_count" =~ ^[1-9][0-9]*$ || ! "$query_count" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: REFERENCE_COUNT and QUERY_COUNT must be positive integers" >&2
    exit 2
fi

mkdir -p -- "$output_dir"
reference="$output_dir/reference-${reference_count}.faa"
query="$output_dir/query-${query_count}.faa"
: > "$reference"
: > "$query"

# Select evenly across the full source, then interleave reference/query
# assignment within that sample. Every output record is complete and
# unmodified; no sequence is duplicated or synthesized.
source_count=$(awk '/^>/ { ++records } END { print records + 0 }' "$source_fasta")
sample_count=$((reference_count + query_count))
if ((source_count < sample_count)); then
    echo "error: source contains $source_count records but $sample_count are required" >&2
    exit 1
fi

awk -v reference="$reference" -v query="$query" \
    -v reference_limit="$reference_count" -v query_limit="$query_count" \
    -v source_count="$source_count" -v sample_limit="$sample_count" '
    BEGIN {
        next_selected = int(0.5 * source_count / sample_limit) + 1
    }
    /^>/ {
        ++seen
        destination = ""
        if (seen == next_selected) {
            ++selected
            desired_queries = int(selected * query_limit / sample_limit)
            if (desired_queries > queries) {
                destination = query
                ++queries
            } else {
                destination = reference
                ++references
            }
            next_selected = int((selected + 0.5) * source_count / sample_limit) + 1
        } else {
            destination = ""
        }
    }
    destination != "" { print > destination }
    END {
        if (references != reference_limit || queries != query_limit) {
            printf "error: source has too few records (selected %d/%d reference, %d/%d query)\n", \
                references, reference_limit, queries, query_limit > "/dev/stderr"
            exit 1
        }
    }
    ' "$source_fasta"

if ! awk '
    /^>/ {
        id = $1
        sub(/^>/, "", id)
        if (seen[id]++)
            duplicate = 1
    }
    END { exit duplicate }
    ' "$reference" "$query"
then
    echo "error: selected records contain duplicate FASTA identifiers" >&2
    exit 1
fi

reference_sha=$(sha256sum "$reference" | awk '{ print $1 }')
query_sha=$(sha256sum "$query" | awk '{ print $1 }')
source_sha=$(sha256sum "$source_fasta" | awk '{ print $1 }')
printf 'source=%s sha256=%s\n' "$source_fasta" "$source_sha"
printf 'reference=%s records=%s sha256=%s\n' "$reference" "$reference_count" "$reference_sha"
printf 'query=%s records=%s sha256=%s\n' "$query" "$query_count" "$query_sha"
