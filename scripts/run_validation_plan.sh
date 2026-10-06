#!/usr/bin/env bash
# Execute the focused checks in VALIDATION_PLAN.md.
set -uo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cpp_bin=${CPP_BIN:-"$repo_dir/diamond/build/diamond"}
rust_bin=${RUST_BIN:-"$repo_dir/target/release/diamond"}
result_dir=${RESULT_DIR:-"$repo_dir/.tmp/validation-$(date +%Y%m%d-%H%M%S)"}
areas=${VALIDATION_AREAS:-"sensitivity clustering scoring view"}
quick_ref=${QUICK_REFERENCE_FASTA:-"$repo_dir/diamond/src/test/data.faa"}
quick_query=${QUICK_QUERY_FASTA:-"$repo_dir/diamond/src/test/5.faa"}
real_ref=${PROTEIN_REFERENCE_FASTA:-}
real_query=${PROTEIN_QUERY_FASTA:-}
cluster_fasta=${CLUSTER_FASTA:-${PROTEIN_REFERENCE_FASTA:-}}
dna_query=${NUCLEOTIDE_QUERY_FASTA:-}
compare_py="$repo_dir/scripts/validation_compare.py"
outfmt=(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)

mkdir -p -- "$result_dir" "$result_dir/artifacts" "$result_dir/db" "$result_dir/tmp"
summary="$result_dir/summary.tsv"
notes="$result_dir/notes.md"
manifest="$result_dir/manifest.tsv"
printf 'date\tarea\tdataset\tcommand\tthreads\tcpp_seconds\trust_seconds\tspeed_ratio_cpp_over_rust\tcpp_rss_kib\trust_rss_kib\trss_ratio_rust_over_cpp\tparity\tstatus\tartifact_dir\tnotes\n' > "$summary"
printf '# Validation scouting notes\n\nRaw, single-run timings below are diagnostic and are not stable benchmark claims.\n\n' > "$notes"
printf 'role\tpath\trecords\tresidues\tsha256\n' > "$manifest"
printf 'source\tvalue\n' > "$result_dir/sources.tsv"
printf 'protein\t%s\n' "${PROTEIN_SOURCE_URL:-unspecified}" >> "$result_dir/sources.tsv"
printf 'protein_version\t%s\n' "${PROTEIN_SOURCE_VERSION:-unspecified}" >> "$result_dir/sources.tsv"
printf 'dna\t%s\n' "${DNA_SOURCE_URL:-unspecified}" >> "$result_dir/sources.tsv"

contains_area() {
    [[ " $areas " == *" $1 "* ]]
}

fasta_manifest() {
    local role=$1 path=$2 records residues hash
    if [[ ! -s "$path" ]]; then
        return
    fi
    records=$(awk '/^>/ { n++ } END { print n + 0 }' "$path")
    residues=$(awk '!/^>/ { gsub(/[[:space:]]/, ""); n += length } END { print n + 0 }' "$path")
    hash=$(sha256sum "$path" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\t%s\n' "$role" "$path" "$records" "$residues" "$hash" >> "$manifest"
}

fasta_manifest quick_reference "$quick_ref"
fasta_manifest quick_query "$quick_query"
[[ -n "$real_ref" ]] && fasta_manifest real_reference "$real_ref"
[[ -n "$real_query" ]] && fasta_manifest real_query "$real_query"
[[ -n "$cluster_fasta" ]] && fasta_manifest clustering "$cluster_fasta"
[[ -n "$dna_query" ]] && fasta_manifest translated_query "$dna_query"
{
    printf 'git_commit\t%s\n' "$(git -C "$repo_dir" rev-parse HEAD)"
    printf 'cpp_version\t%s\n' "$($cpp_bin version 2>&1 | head -n 1)"
    printf 'rust_version\t%s\n' "$($rust_bin version 2>&1 | head -n 1)"
    printf 'rustc\t%s\n' "$(rustc --version)"
    printf 'host\t%s\n' "$(uname -a)"
} >> "$manifest"

if [[ ! -x "$cpp_bin" ]]; then
    echo "error: C++ binary is missing: $cpp_bin" >&2
    exit 2
fi
if [[ ! -x "$rust_bin" ]]; then
    echo "error: Rust binary is missing: $rust_bin (run cargo build --release)" >&2
    exit 2
fi

run_timed() {
    local prefix=$1
    shift
    local timing="${prefix}.time" stdout="${prefix}.stdout" stderr="${prefix}.stderr"
    printf '%q ' "$@" > "${prefix}.command"
    printf '\n' >> "${prefix}.command"
    /usr/bin/time -f '%e\t%M' -o "$timing" "$@" > "$stdout" 2> "$stderr"
    RUN_STATUS=$?
    if [[ -s "$timing" ]]; then
        IFS=$'\t' read -r RUN_SECONDS RUN_RSS < <(tail -n 1 "$timing")
    else
        RUN_SECONDS=NA
        RUN_RSS=NA
    fi
}

ratio() {
    local numerator=$1 denominator=$2
    if [[ "$numerator" == NA || "$denominator" == NA ]]; then
        printf 'NA'
    else
        awk -v a="$numerator" -v b="$denominator" 'BEGIN { if (b == 0) print "NA"; else printf "%.3f", a / b }'
    fi
}

record_blocked() {
    local area=$1 dataset=$2 command=$3 threads=$4 why=$5 artifact=${6:-}
    printf '%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNOT_RUN\tBLOCKED\t%s\t%s\n' \
        "$(date -Iseconds)" "$area" "$dataset" "$command" "$threads" "$artifact" "$why" >> "$summary"
}

record_pair() {
    local area=$1 dataset=$2 command=$3 threads=$4 artifact=$5 comparator=$6 cpp_out=$7 rust_out=$8 note=${9:-}
    local parity=FAIL status=FAIL compare_log="$artifact/compare.log"
    if ((CPP_STATUS == 0 && RUST_STATUS == 0)); then
        local compare_command=(python3 "$compare_py" "$comparator" "$cpp_out" "$rust_out")
        if [[ "$comparator" == cluster && -n "${CLUSTER_COMPARE_FASTA:-}" ]]; then
            compare_command+=(--fasta "$CLUSTER_COMPARE_FASTA")
        fi
        if "${compare_command[@]}" > "$compare_log" 2>&1; then
            parity=PASS
            status=PASS
        fi
    else
        printf 'C++ exit=%s; Rust exit=%s\n' "$CPP_STATUS" "$RUST_STATUS" > "$compare_log"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(date -Iseconds)" "$area" "$dataset" "$command" "$threads" \
        "$CPP_SECONDS" "$RUST_SECONDS" "$(ratio "$CPP_SECONDS" "$RUST_SECONDS")" \
        "$CPP_RSS" "$RUST_RSS" "$(ratio "$RUST_RSS" "$CPP_RSS")" \
        "$parity" "$status" "$artifact" "$note" >> "$summary"
    [[ "$status" == PASS ]]
}

make_databases() {
    local label=$1 fasta=$2
    local dir="$result_dir/db/$label"
    mkdir -p -- "$dir"
    run_timed "$dir/cpp-makedb" "$cpp_bin" makedb --in "$fasta" -d "$dir/cpp"
    local cpp_status=$RUN_STATUS
    run_timed "$dir/rust-makedb" "$rust_bin" makedb --in "$fasta" -d "$dir/rust"
    local rust_status=$RUN_STATUS
    if ((cpp_status != 0 || rust_status != 0)); then
        echo "error: makedb failed for $label; see $dir" >&2
        return 1
    fi
    printf '%s\t%s\n' "$dir/cpp" "$dir/rust"
}

quick_dbs=$(make_databases quick "$quick_ref") || exit 1
IFS=$'\t' read -r quick_cpp_db quick_rust_db <<< "$quick_dbs"
real_cpp_db=
real_rust_db=
if [[ -n "$real_ref" && -n "$real_query" && -s "$real_ref" && -s "$real_query" ]]; then
    real_dbs=$(make_databases real "$real_ref") || exit 1
    IFS=$'\t' read -r real_cpp_db real_rust_db <<< "$real_dbs"
fi

run_search_case() {
    local dataset=$1 workflow=$2 label=$3 query=$4 cpp_db=$5 rust_db=$6 threads=$7
    shift 7
    local artifact="$result_dir/artifacts/sensitivity-${dataset}-${workflow}-${label}-t${threads}"
    mkdir -p -- "$artifact"
    run_timed "$artifact/cpp" "$cpp_bin" "$workflow" -d "$cpp_db" -q "$query" -o "$artifact/cpp.tsv" \
        --threads "$threads" --outfmt 6 "${outfmt[@]}" "$@"
    CPP_STATUS=$RUN_STATUS CPP_SECONDS=$RUN_SECONDS CPP_RSS=$RUN_RSS
    run_timed "$artifact/rust" "$rust_bin" "$workflow" -d "$rust_db" -q "$query" -o "$artifact/rust.tsv" \
        --threads "$threads" --outfmt 6 "${outfmt[@]}" "$@"
    RUST_STATUS=$RUN_STATUS RUST_SECONDS=$RUN_SECONDS RUST_RSS=$RUN_RSS
    record_pair sensitivity "$dataset" "$workflow:$label" "$threads" "$artifact" exact "$artifact/cpp.tsv" "$artifact/rust.tsv" || true
}

run_search_disk_case() {
    local dataset=$1 workflow=$2 label=$3 query=$4 cpp_db=$5 rust_db=$6 threads=$7
    shift 7
    local artifact="$result_dir/artifacts/sensitivity-${dataset}-${workflow}-${label}-disk0-t${threads}"
    mkdir -p -- "$artifact"
    run_timed "$artifact/cpp" "$cpp_bin" "$workflow" -d "$cpp_db" -q "$query" -o "$artifact/cpp.tsv" \
        --threads "$threads" --outfmt 6 "${outfmt[@]}" "$@"
    CPP_STATUS=$RUN_STATUS CPP_SECONDS=$RUN_SECONDS CPP_RSS=$RUN_RSS
    run_timed "$artifact/rust" "$rust_bin" "$workflow" -d "$rust_db" -q "$query" -o "$artifact/rust.tsv" \
        --threads "$threads" --memory-limit 0G --outfmt 6 "${outfmt[@]}" "$@"
    RUST_STATUS=$RUN_STATUS RUST_SECONDS=$RUN_SECONDS RUST_RSS=$RUN_RSS
    record_pair sensitivity "$dataset" "$workflow:$label:disk0" "$threads" "$artifact" exact "$artifact/cpp.tsv" "$artifact/rust.tsv" || true
}

run_sensitivities_for() {
    local dataset=$1 workflow=$2 query=$3 cpp_db=$4 rust_db=$5 preset
    for preset in default faster fast mid-sensitive sensitive more-sensitive very-sensitive ultra-sensitive; do
        if [[ "$preset" == default ]]; then
            run_search_case "$dataset" "$workflow" true-default "$query" "$cpp_db" "$rust_db" 1
        else
            run_search_case "$dataset" "$workflow" "$preset" "$query" "$cpp_db" "$rust_db" 1 "--$preset"
        fi
    done
}

if contains_area sensitivity; then
    printf '## Search sensitivity\n\n' >> "$notes"
    run_sensitivities_for quick blastp "$quick_query" "$quick_cpp_db" "$quick_rust_db"
    run_sensitivities_for quick blastx "${dna_query:-$repo_dir/diamond/src/test/galaxy/nucleotide.fasta}" "$quick_cpp_db" "$quick_rust_db"
    if [[ -n "$real_cpp_db" ]]; then
        run_sensitivities_for real blastp "$real_query" "$real_cpp_db" "$real_rust_db"
        if [[ -n "$dna_query" && -s "$dna_query" ]]; then
            run_sensitivities_for real blastx "$dna_query" "$real_cpp_db" "$real_rust_db"
        else
            record_blocked sensitivity real blastx 1 'Set NUCLEOTIDE_QUERY_FASTA to a natural genomic region or assembly.'
        fi
    else
        record_blocked sensitivity real blastp 1 'Set PROTEIN_REFERENCE_FASTA and PROTEIN_QUERY_FASTA to a disjoint natural fixture.'
        record_blocked sensitivity real blastx 1 'Set real protein fixtures and NUCLEOTIDE_QUERY_FASTA.'
    fi
    # Limited multithread scouting: control and deliberately expensive preset.
    run_search_case quick blastp true-default "$quick_query" "$quick_cpp_db" "$quick_rust_db" 4
    run_search_case quick blastp ultra-sensitive "$quick_query" "$quick_cpp_db" "$quick_rust_db" 4 --ultra-sensitive
    run_search_case quick blastx true-default "${dna_query:-$repo_dir/diamond/src/test/galaxy/nucleotide.fasta}" "$quick_cpp_db" "$quick_rust_db" 4
    run_search_case quick blastx ultra-sensitive "${dna_query:-$repo_dir/diamond/src/test/galaxy/nucleotide.fasta}" "$quick_cpp_db" "$quick_rust_db" 4 --ultra-sensitive
    run_search_disk_case quick blastp true-default "$quick_query" "$quick_cpp_db" "$quick_rust_db" 4
    run_search_disk_case quick blastp ultra-sensitive "$quick_query" "$quick_cpp_db" "$quick_rust_db" 4 --ultra-sensitive
    run_search_disk_case quick blastx true-default "${dna_query:-$repo_dir/diamond/src/test/galaxy/nucleotide.fasta}" "$quick_cpp_db" "$quick_rust_db" 4
    run_search_disk_case quick blastx ultra-sensitive "${dna_query:-$repo_dir/diamond/src/test/galaxy/nucleotide.fasta}" "$quick_cpp_db" "$quick_rust_db" 4 --ultra-sensitive
fi

run_cluster_case() {
    local dataset=$1 command=$2 label=$3 fasta=$4 cpp_db=$5 rust_db=$6 threads=$7
    shift 7
    local artifact="$result_dir/artifacts/clustering-${dataset}-${command}-${label}-t${threads}"
    mkdir -p -- "$artifact"
    run_timed "$artifact/cpp" "$cpp_bin" "$command" -d "${cpp_db}.dmnd" -o "$artifact/cpp.tsv" --threads "$threads" "$@"
    CPP_STATUS=$RUN_STATUS CPP_SECONDS=$RUN_SECONDS CPP_RSS=$RUN_RSS
    run_timed "$artifact/rust" "$rust_bin" "$command" -d "${rust_db}.dmnd" -o "$artifact/rust.tsv" --threads "$threads" "$@"
    RUST_STATUS=$RUN_STATUS RUST_SECONDS=$RUN_SECONDS RUST_RSS=$RUN_RSS
    CLUSTER_COMPARE_FASTA=$fasta
    record_pair clustering "$dataset" "$command:$label" "$threads" "$artifact" cluster "$artifact/cpp.tsv" "$artifact/rust.tsv" || true
    unset CLUSTER_COMPARE_FASTA
}

if contains_area clustering; then
    printf '\n## Clustering\n\n' >> "$notes"
    for command in cluster linclust deepclust; do
        run_cluster_case quick "$command" default "$quick_ref" "$quick_cpp_db" "$quick_rust_db" 1
        run_cluster_case quick "$command" id90-cover80 "$quick_ref" "$quick_cpp_db" "$quick_rust_db" 1 --approx-id 90 --member-cover 80
    done
    if [[ -n "$cluster_fasta" && -s "$cluster_fasta" ]]; then
        cluster_dbs=$(make_databases cluster "$cluster_fasta") || exit 1
        IFS=$'\t' read -r cluster_cpp_db cluster_rust_db <<< "$cluster_dbs"
        for command in cluster linclust deepclust; do
            run_cluster_case real "$command" default "$cluster_fasta" "$cluster_cpp_db" "$cluster_rust_db" 1
            run_cluster_case real "$command" id90-cover80 "$cluster_fasta" "$cluster_cpp_db" "$cluster_rust_db" 1 --approx-id 90 --member-cover 80
            run_cluster_case real "$command" id50-cover60 "$cluster_fasta" "$cluster_cpp_db" "$cluster_rust_db" 1 --approx-id 50 --member-cover 60
            run_cluster_case real "$command" default "$cluster_fasta" "$cluster_cpp_db" "$cluster_rust_db" 4
        done
    else
        for command in cluster linclust deepclust; do
            record_blocked clustering real "$command" 1 'Set CLUSTER_FASTA to a natural family-rich protein collection.'
        done
    fi
fi

matrix_spec() {
    case "$1" in
        BLOSUM45) printf '14 2 13 3' ;;
        BLOSUM50) printf '13 2 13 3' ;;
        BLOSUM62) printf '11 1 11 2' ;;
        BLOSUM80) printf '10 1 25 2' ;;
        BLOSUM90) printf '10 1 9 2' ;;
        PAM30) printf '9 1 7 2' ;;
        PAM70) printf '10 1 8 2' ;;
        PAM250) printf '14 2 15 3' ;;
    esac
}

run_scoring_case() {
    local dataset=$1 matrix=$2 label=$3 query=$4 cpp_db=$5 rust_db=$6
    shift 6
    local artifact="$result_dir/artifacts/scoring-${dataset}-${matrix}-${label}"
    mkdir -p -- "$artifact"
    run_timed "$artifact/cpp" "$cpp_bin" blastp -d "$cpp_db" -q "$query" -o "$artifact/cpp.tsv" \
        --threads 1 --outfmt 6 "${outfmt[@]}" --matrix "$matrix" "$@"
    CPP_STATUS=$RUN_STATUS CPP_SECONDS=$RUN_SECONDS CPP_RSS=$RUN_RSS
    run_timed "$artifact/rust" "$rust_bin" blastp -d "$rust_db" -q "$query" -o "$artifact/rust.tsv" \
        --threads 1 --outfmt 6 "${outfmt[@]}" --matrix "$matrix" "$@"
    RUST_STATUS=$RUN_STATUS RUST_SECONDS=$RUN_SECONDS RUST_RSS=$RUN_RSS
    record_pair scoring "$dataset" "$matrix:$label" 1 "$artifact" exact "$artifact/cpp.tsv" "$artifact/rust.tsv" || true
}

run_rejected_gap() {
    local matrix=$1
    local artifact="$result_dir/artifacts/scoring-quick-${matrix}-rejected-gap"
    mkdir -p -- "$artifact"
    run_timed "$artifact/cpp" "$cpp_bin" blastp -d "$quick_cpp_db" -q "$quick_query" -o "$artifact/cpp.tsv" --threads 1 --matrix "$matrix" --gapopen 1 --gapextend 1
    CPP_STATUS=$RUN_STATUS CPP_SECONDS=$RUN_SECONDS CPP_RSS=$RUN_RSS
    run_timed "$artifact/rust" "$rust_bin" blastp -d "$quick_rust_db" -q "$quick_query" -o "$artifact/rust.tsv" --threads 1 --matrix "$matrix" --gapopen 1 --gapextend 1
    RUST_STATUS=$RUN_STATUS RUST_SECONDS=$RUN_SECONDS RUST_RSS=$RUN_RSS
    local parity=FAIL status=FAIL
    if ((CPP_STATUS != 0 && RUST_STATUS != 0)) && python3 "$compare_py" errors "$artifact/cpp.stderr" "$artifact/rust.stderr" > "$artifact/compare.log" 2>&1; then
        parity=PASS status=PASS
    fi
    printf '%s\tscoring\tquick\t%s:rejected-gap\t1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\trejection parity\n' \
        "$(date -Iseconds)" "$matrix" "$CPP_SECONDS" "$RUST_SECONDS" "$(ratio "$CPP_SECONDS" "$RUST_SECONDS")" \
        "$CPP_RSS" "$RUST_RSS" "$(ratio "$RUST_RSS" "$CPP_RSS")" "$parity" "$status" "$artifact" >> "$summary"
}

run_scoring_spot_case() {
    local matrix=$1 query=$2 cpp_db=$3 rust_db=$4
    local artifact="$result_dir/artifacts/scoring-real-${matrix}-default-t4"
    mkdir -p -- "$artifact"
    run_timed "$artifact/cpp" "$cpp_bin" blastp -d "$cpp_db" -q "$query" -o "$artifact/cpp.tsv" \
        --threads 4 --outfmt 6 "${outfmt[@]}" --matrix "$matrix"
    CPP_STATUS=$RUN_STATUS CPP_SECONDS=$RUN_SECONDS CPP_RSS=$RUN_RSS
    run_timed "$artifact/rust" "$rust_bin" blastp -d "$rust_db" -q "$query" -o "$artifact/rust.tsv" \
        --threads 4 --outfmt 6 "${outfmt[@]}" --matrix "$matrix"
    RUST_STATUS=$RUN_STATUS RUST_SECONDS=$RUN_SECONDS RUST_RSS=$RUN_RSS
    record_pair scoring real "$matrix:default:t4" 4 "$artifact" exact "$artifact/cpp.tsv" "$artifact/rust.tsv" || true
}

if contains_area scoring; then
    printf '\n## Scoring matrices\n\n' >> "$notes"
    if [[ "${VALIDATION_SCORING_SPOTS_ONLY:-0}" != 1 ]]; then
        for matrix in BLOSUM62 BLOSUM45 BLOSUM50 BLOSUM80 BLOSUM90 PAM30 PAM70 PAM250; do
            read -r default_go default_ge alternate_go alternate_ge <<< "$(matrix_spec "$matrix")"
            run_scoring_case quick "$matrix" default "$quick_query" "$quick_cpp_db" "$quick_rust_db"
            run_scoring_case quick "$matrix" explicit-default "$quick_query" "$quick_cpp_db" "$quick_rust_db" --gapopen "$default_go" --gapextend "$default_ge"
            run_scoring_case quick "$matrix" alternate "$quick_query" "$quick_cpp_db" "$quick_rust_db" --gapopen "$alternate_go" --gapextend "$alternate_ge"
            run_rejected_gap "$matrix"
            if [[ -n "$real_cpp_db" ]]; then
                run_scoring_case real "$matrix" default "$real_query" "$real_cpp_db" "$real_rust_db"
                run_scoring_case real "$matrix" alternate "$real_query" "$real_cpp_db" "$real_rust_db" --gapopen "$alternate_go" --gapextend "$alternate_ge"
            else
                record_blocked scoring real "$matrix" 1 'Set PROTEIN_REFERENCE_FASTA and PROTEIN_QUERY_FASTA.'
            fi
        done
    fi
    if [[ -n "$real_cpp_db" ]]; then
        # Minimal multithread scouting: the lowest observed natural-data speed
        # ratio in the baseline pass and one representative PAM matrix.
        run_scoring_spot_case BLOSUM62 "$real_query" "$real_cpp_db" "$real_rust_db"
        run_scoring_spot_case PAM70 "$real_query" "$real_cpp_db" "$real_rust_db"
    fi
fi

run_view_case() {
    local dataset=$1 label=$2 format=$3 comparator=$4 daa=$5
    shift 5
    local safe_label=${label//[^A-Za-z0-9_-]/_}
    local artifact="$result_dir/artifacts/view-${dataset}-${safe_label}"
    mkdir -p -- "$artifact"
    run_timed "$artifact/cpp" "$cpp_bin" view -a "$daa" -o "$artifact/cpp.out" --outfmt "$format" "$@"
    CPP_STATUS=$RUN_STATUS CPP_SECONDS=$RUN_SECONDS CPP_RSS=$RUN_RSS
    run_timed "$artifact/rust" "$rust_bin" view -a "$daa" -o "$artifact/rust.out" --outfmt "$format" "$@"
    RUST_STATUS=$RUN_STATUS RUST_SECONDS=$RUN_SECONDS RUST_RSS=$RUN_RSS
    record_pair view "$dataset" "view:$label" 1 "$artifact" "$comparator" "$artifact/cpp.out" "$artifact/rust.out" || true
}

run_view_unal_rejection() {
    local dataset=$1 value=$2 daa=$3
    local artifact="$result_dir/artifacts/view-${dataset}-unal-${value}"
    local diagnostic='Option is not permitted for this workflow: unal'
    local parity=FAIL status=FAIL
    mkdir -p -- "$artifact"
    run_timed "$artifact/cpp" "$cpp_bin" view -a "$daa" -o "$artifact/cpp.out" --unal "$value"
    CPP_STATUS=$RUN_STATUS CPP_SECONDS=$RUN_SECONDS CPP_RSS=$RUN_RSS
    run_timed "$artifact/rust" "$rust_bin" view -a "$daa" -o "$artifact/rust.out" --unal "$value"
    RUST_STATUS=$RUN_STATUS RUST_SECONDS=$RUN_SECONDS RUST_RSS=$RUN_RSS
    if ((CPP_STATUS != 0 && RUST_STATUS != 0)) \
        && grep -Fq "$diagnostic" "$artifact/cpp.stderr" \
        && grep -Fq "$diagnostic" "$artifact/rust.stderr"; then
        parity=PASS
        status=PASS
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(date -Iseconds)" view "$dataset" "view:unal-$value-rejected" 1 \
        "$CPP_SECONDS" "$RUST_SECONDS" "$(ratio "$CPP_SECONDS" "$RUST_SECONDS")" \
        "$CPP_RSS" "$RUST_RSS" "$(ratio "$RUST_RSS" "$CPP_RSS")" \
        "$parity" "$status" "$artifact" 'both implementations must reject this option' >> "$summary"
}

run_daa_view_suite() {
    local dataset=$1 workflow=$2 query=$3 cpp_db=$4 rust_db=$5
    local daa_artifact="$result_dir/artifacts/daa-${dataset}-producers"
    mkdir -p -- "$daa_artifact"
    run_timed "$daa_artifact/cpp-${workflow}" "$cpp_bin" "$workflow" -d "$cpp_db" -q "$query" -a "$daa_artifact/archive" --threads 1
    local producer_cpp_status=$RUN_STATUS producer_cpp_seconds=$RUN_SECONDS producer_cpp_rss=$RUN_RSS
    run_timed "$daa_artifact/rust-${workflow}" "$rust_bin" "$workflow" -d "$rust_db" -q "$query" -a "$daa_artifact/rust-archive" --threads 1
    local producer_rust_status=$RUN_STATUS producer_rust_seconds=$RUN_SECONDS producer_rust_rss=$RUN_RSS
    if ((producer_cpp_status != 0)) || [[ ! -s "$daa_artifact/archive.daa" ]]; then
        record_blocked view "$dataset" cpp-daa-producer 1 'Upstream failed to produce the baseline DAA.' "$daa_artifact"
        return
    fi
    if ((producer_rust_status != 0)) || [[ ! -s "$daa_artifact/rust-archive.daa" ]]; then
        record_blocked view "$dataset" rust-daa-producer 1 'Native Rust failed to produce its DAA.' "$daa_artifact"
        return
    fi
    local cpp_daa="$daa_artifact/archive.daa"
    local rust_daa="$daa_artifact/rust-archive.daa"

    # Compare producer semantics through the native Rust viewer. Raw archives
    # need not be byte-identical because equivalent packed transcripts can use
    # different run lengths; rendered alignments must be identical.
    "$rust_bin" view -a "$cpp_daa" -o "$daa_artifact/cpp-producer.tsv" --outfmt 6 "${outfmt[@]}" >"$daa_artifact/cpp-render.stdout" 2>"$daa_artifact/cpp-render.stderr"
    local cpp_render_status=$?
    "$rust_bin" view -a "$rust_daa" -o "$daa_artifact/rust-producer.tsv" --outfmt 6 "${outfmt[@]}" >"$daa_artifact/rust-render.stdout" 2>"$daa_artifact/rust-render.stderr"
    local rust_render_status=$?
    CPP_STATUS=$producer_cpp_status CPP_SECONDS=$producer_cpp_seconds CPP_RSS=$producer_cpp_rss
    RUST_STATUS=$producer_rust_status RUST_SECONDS=$producer_rust_seconds RUST_RSS=$producer_rust_rss
    if ((cpp_render_status != 0 || rust_render_status != 0)); then
        CPP_STATUS=1
        RUST_STATUS=1
    fi
    record_pair view "$dataset" "$workflow:daa-producer" 1 "$daa_artifact" exact "$daa_artifact/cpp-producer.tsv" "$daa_artifact/rust-producer.tsv" || true

    local producer daa
    for producer in cpp rust; do
        if [[ "$producer" == cpp ]]; then daa=$cpp_daa; else daa=$rust_daa; fi
        run_view_case "$dataset-$producer" daa 100 exact "$daa"
        run_view_case "$dataset-$producer" tab 6 exact "$daa" "${outfmt[@]}"
        # Pinned upstream view output is malformed for a multi-query DAA: its
        # ViewWriter omits OutputFormat::query_separator between query buffers.
        # Translation fidelity therefore requires exact bytes, not a repaired
        # JSON document.
        run_view_case "$dataset-$producer" json-flat 104 exact "$daa"
        run_view_case "$dataset-$producer" paf 103 sorted-lines "$daa"
        run_view_case "$dataset-$producer" sam 101 sam "$daa"
        run_view_case "$dataset-$producer" pairwise 0 exact "$daa"
        run_view_case "$dataset-$producer" xml 5 xml "$daa"
        run_view_case "$dataset-$producer" null null exact "$daa"
        run_view_case "$dataset-$producer" edge edge exact "$daa"
        run_view_case "$dataset-$producer" tab-k1 6 exact "$daa" "${outfmt[@]}" --max-target-seqs 1
        run_view_case "$dataset-$producer" tab-top10 6 exact "$daa" "${outfmt[@]}" --top 10
        run_view_unal_rejection "$dataset-$producer" 0 "$daa"
        run_view_unal_rejection "$dataset-$producer" 1 "$daa"
        if [[ "$workflow" == blastx ]]; then
            run_view_case "$dataset-$producer" tab-forward 6 exact "$daa" "${outfmt[@]}" --forwardonly
        fi
    done
}

if contains_area view; then
    printf '\n## DAA and view\n\n' >> "$notes"
    if [[ -n "$real_cpp_db" ]]; then
        run_daa_view_suite blastp-daa blastp "$real_query" "$real_cpp_db" "$real_rust_db"
        run_daa_view_suite blastx-daa blastx "${dna_query:-$repo_dir/diamond/src/test/galaxy/nucleotide.fasta}" "$real_cpp_db" "$real_rust_db"
    else
        run_daa_view_suite blastp-daa blastp "$quick_query" "$quick_cpp_db" "$quick_rust_db"
        run_daa_view_suite blastx-daa blastx "$repo_dir/diamond/src/test/galaxy/nucleotide.fasta" "$repo_dir/diamond/src/test/galaxy/db.dmnd" "$repo_dir/diamond/src/test/galaxy/db.dmnd"
    fi
fi

python3 - "$summary" "$notes" <<'PY'
import csv, sys
summary, notes = sys.argv[1:]
with open(summary, newline='') as stream:
    rows = list(csv.DictReader(stream, delimiter='\t'))
fails = [r for r in rows if r['status'] == 'FAIL']
blocked = [r for r in rows if r['status'] == 'BLOCKED']
timed = [r for r in rows if r['speed_ratio_cpp_over_rust'] not in ('NA', '') and r['rss_ratio_rust_over_cpp'] not in ('NA', '')]
slow = sorted(timed, key=lambda r: float(r['speed_ratio_cpp_over_rust']))[:2]
rss = sorted(timed, key=lambda r: float(r['rss_ratio_rust_over_cpp']), reverse=True)[:2]
with open(notes, 'a') as out:
    out.write(f"\n## Triage summary\n\n- PASS: {sum(r['status'] == 'PASS' for r in rows)}\n")
    out.write(f"- EXPECTED DIFFERENCE: {sum(r['status'] == 'EXPECTED_DIFFERENCE' for r in rows)}\n")
    out.write(f"- FAIL: {len(fails)}\n- BLOCKED: {len(blocked)}\n")
    out.write("- Slowest observed: " + ", ".join(f"{r['area']}/{r['command']} ({r['speed_ratio_cpp_over_rust']}x)" for r in slow) + "\n")
    out.write("- Largest RSS ratios: " + ", ".join(f"{r['area']}/{r['command']} ({r['rss_ratio_rust_over_cpp']}x)" for r in rss) + "\n")
    if fails:
        out.write("\n### Correctness failures\n\n")
        for row in fails:
            out.write(f"- `{row['area']}/{row['dataset']}/{row['command']}`: `{row['artifact_dir']}`\n")
    if blocked:
        out.write("\n### Blocked coverage\n\n")
        for row in blocked:
            out.write(f"- `{row['area']}/{row['dataset']}/{row['command']}`: {row['notes']}\n")
PY

echo "Validation summary: $summary"
echo "Triage notes:       $notes"
echo "Fixture manifest:   $manifest"
if awk -F '\t' 'NR > 1 && $13 == "FAIL" { found=1 } END { exit !found }' "$summary"; then
    exit 1
fi
