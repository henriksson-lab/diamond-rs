#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source_fasta=${REAL_PROTEIN_FASTA:-}
fixture_dir=${FIXTURE_DIR:-/tmp/diamond-real-benchmark}
reference_count=${REFERENCE_COUNT:-2000}
query_count=${QUERY_COUNT:-1000}

if [[ -z "$source_fasta" ]]; then
    echo "error: set REAL_PROTEIN_FASTA to a diverse, uncompressed protein FASTA" >&2
    echo "example: REAL_PROTEIN_FASTA=/path/to/uniprot_sprot_human.fasta $0" >&2
    exit 2
fi

REFERENCE_COUNT=$reference_count QUERY_COUNT=$query_count \
    "$repo_dir/scripts/prepare_real_benchmark.sh" "$source_fasta" "$fixture_dir"

REFERENCE_FASTA="$fixture_dir/reference-${reference_count}.faa" \
QUERY_FASTA="$fixture_dir/query-${query_count}.faa" \
THREADS=${THREADS:-4} REPETITIONS=${REPETITIONS:-5} KEEP_WORK=${KEEP_WORK:-1} \
    "$repo_dir/scripts/compare_cpp_rust.sh"

