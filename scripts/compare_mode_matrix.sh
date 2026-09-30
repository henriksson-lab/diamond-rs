#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
driver="$repo_dir/scripts/compare_cpp_rust.sh"
reference=${REFERENCE_FASTA:-"$repo_dir/diamond/src/test/data.faa"}
protein_query=${PROTEIN_QUERY_FASTA:-${QUERY_FASTA:-"$repo_dir/diamond/src/test/5.faa"}}
nucleotide_query=${NUCLEOTIDE_QUERY_FASTA:-}
repetitions=${REPETITIONS:-5}
matrix_threads=${MATRIX_THREADS:-"1 4"}
matrix_modes=${MATRIX_MODES:-"default16 spill64m disk0 gzip-input faster more ultra cbs0 mask0 blosum45 filter-id filter-cover filter-top output-extra blastx-default blastx-disk blastx-very blastx-ultra blastx-cbs0 blastx-mask0 blastx-strand-plus blastx-strand-minus blastx-gencode11 blastx-min-orf blastx-output-extra"}
result_dir=${RESULT_DIR:-"$repo_dir/.tmp/benchmark-matrix"}
summary="$result_dir/summary.tsv"

mkdir -p -- "$result_dir"
printf 'label\tcommand\tthreads\truns\trust_memory_limit\tsearch_args\tcpp_seconds\trust_seconds\tspeed_ratio\tcpp_rss_kib\trust_rss_kib\trss_ratio\tparity\n' > "$summary"

run_mode() {
    local label=$1
    local command=$2
    local memory_limit=$3
    local query=$4
    local search_args=$5
    local threads=$6
    local outfmt=$7
    local cpuset=''
    if [[ "$threads" == 1 ]]; then
        cpuset=${BENCH_CPUSET_1:-}
    elif [[ "$threads" == 4 ]]; then
        cpuset=${BENCH_CPUSET_4:-}
    fi

    printf '\n=== %s, %s thread(s) ===\n' "$label" "$threads"
    BENCH_LABEL="$label" COMMAND="$command" SEARCH_ARGS="$search_args" \
        RUST_MEMORY_LIMIT="$memory_limit" REFERENCE_FASTA="$reference" \
        QUERY_FASTA="$query" THREADS="$threads" REPETITIONS="$repetitions" \
        OUTFMT_FIELDS="$outfmt" \
        BENCH_CPUSET="$cpuset" SUMMARY_FILE="$summary" SKIP_BUILD=1 KEEP_WORK=0 \
        "$driver" | tee "$result_dir/${label}-t${threads}.log"
}

failures=0
for mode in $matrix_modes; do
    command=blastp
    mode_reference=$reference
    query=$protein_query
    memory_limit=0G
    search_args=''
    outfmt='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    case "$mode" in
        default16) memory_limit='' ;;
        spill64m) memory_limit=64M ;;
        disk0) ;;
        gzip-input)
            mkdir -p -- "$result_dir/compressed-input"
            gzip -c -- "$reference" > "$result_dir/compressed-input/reference.faa.gz"
            gzip -c -- "$protein_query" > "$result_dir/compressed-input/query.faa.gz"
            mode_reference="$result_dir/compressed-input/reference.faa.gz"
            query="$result_dir/compressed-input/query.faa.gz"
            ;;
        faster) search_args='--faster' ;;
        more) search_args='--more-sensitive' ;;
        ultra) search_args='--ultra-sensitive' ;;
        cbs0) search_args='--comp-based-stats 0' ;;
        mask0) search_args='--masking 0' ;;
        blosum45) search_args='--matrix BLOSUM45 --gapopen 14 --gapextend 2' ;;
        filter-id) search_args='--id 70' ;;
        filter-cover) search_args='--query-cover 50 --subject-cover 50' ;;
        filter-top) search_args='--top 10' ;;
        output-extra)
            outfmt='qseqid qlen sseqid slen pident nident length mismatch gapopen gaps qstart qend sstart send evalue bitscore score'
            ;;
        blastx-default | blastx-disk | blastx-very | blastx-very-default | \
        blastx-very-disk | blastx-ultra | \
        blastx-cbs0 | blastx-mask0 | blastx-strand-plus | \
        blastx-strand-minus | blastx-gencode11 | blastx-min-orf | \
        blastx-output-extra)
            if [[ -z "$nucleotide_query" ]]; then
                echo "Skipping $mode: set NUCLEOTIDE_QUERY_FASTA to a real DNA FASTA." >&2
                continue
            fi
            command=blastx
            query=$nucleotide_query
            # A natural 250 kb region reaches the 3--5 second diagnostic
            # range at this preset and exercises several translated-search
            # shapes rather than measuring mostly fixed startup cost.
            search_args='--sensitive'
            [[ "$mode" == blastx-default || "$mode" == blastx-very-default ]] && memory_limit=''
            case "$mode" in
                blastx-very | blastx-very-default | blastx-very-disk)
                    search_args='--very-sensitive'
                    ;;
                blastx-ultra) search_args='--ultra-sensitive' ;;
                blastx-cbs0) search_args='--sensitive --comp-based-stats 0' ;;
                blastx-mask0) search_args='--sensitive --masking 0' ;;
                blastx-strand-plus) search_args='--sensitive --strand plus' ;;
                blastx-strand-minus) search_args='--sensitive --strand minus' ;;
                blastx-gencode11) search_args='--sensitive --query-gencode 11' ;;
                blastx-min-orf) search_args='--sensitive --min-orf 30' ;;
                blastx-output-extra)
                    outfmt='qseqid qlen sseqid slen pident nident length mismatch gapopen gaps qstart qend sstart send qframe evalue bitscore score'
                    ;;
            esac
            ;;
        *)
            echo "error: unknown MATRIX_MODES entry: $mode" >&2
            exit 2
            ;;
    esac

    for threads in $matrix_threads; do
        reference=$mode_reference
        if ! run_mode "$mode" "$command" "$memory_limit" "$query" "$search_args" "$threads" "$outfmt"; then
            failures=$((failures + 1))
        fi
    done
    reference=${REFERENCE_FASTA:-"$repo_dir/diamond/src/test/data.faa"}
done

printf '\nMatrix summary: %s\n' "$summary"
if ((failures)); then
    echo "error: $failures matrix run(s) failed (including parity failures)" >&2
    exit 1
fi
