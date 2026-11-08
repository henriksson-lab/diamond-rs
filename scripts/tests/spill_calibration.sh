#!/usr/bin/env bash
# Deterministic control-flow checks; no search binary, network, or Cargo.
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
source "$root/scripts/lib/spill_calibration.sh"

calls=0
limits=""
run_case() {
    calls=$((calls + 1))
    limits="$limits $2"
    case "$calls" in
        1) RUN_CLASS=delayed-spill ;;
        2|3) RUN_CLASS=no-spill ;;
        4|5) RUN_CLASS=delayed-spill ;;
        *) echo 'unexpected extra probe' >&2; return 9 ;;
    esac
}
select_stable_delayed_limit 196 192 128
[[ "$calls" == 5 && "$limits" == ' 196 196 196 192 192' ]]
[[ "$SELECTED_LIMIT" == 192 ]]
[[ "$GATE_SAMPLE_1" == validation-2-1 && "$GATE_SAMPLE_2" == validation-2-2 ]]
printf 'PASS noisy calibration retains failed samples and requires a consecutive same-limit pair\n'

calls=0
run_case() {
    calls=$((calls + 1))
    RUN_CLASS=no-spill
}
if select_stable_delayed_limit 196 192 128; then
    echo 'FAIL: absent delayed transition was accepted' >&2
    exit 1
fi
[[ "$calls" == 18 && "$SELECTED_LIMIT" == 0 ]]
printf 'PASS no genuine delayed transition fails after 18 bounded samples\n'

calls=0
run_case() {
    calls=$((calls + 1))
    return 17 # Stand-in for a search/diagnostic/output-integrity failure.
}
status=0
select_stable_delayed_limit 196 192 128 || status=$?
[[ "$status" == 17 && "$calls" == 1 ]]
printf 'PASS physical/search/integrity failures propagate immediately without retry\n'

# A checked caller disables errexit throughout the callee. Test the production
# run_case itself with a failing fake search that leaves plausible old output.
# Keeping a later successful command reproduces the original masking hazard.
case_definition=$(sed -n '/^run_case() {/,/^}/p' "$root/scripts/run_native_spill_conformance.sh")
eval "$case_definition"
artifact_dir=$(mktemp -d)
trap 'rm -rf "$artifact_dir"' EXIT
mkdir -p "$artifact_dir/runs" "$artifact_dir/spill-tmp"
cat >"$artifact_dir/failing-search" <<'FAKE'
#!/usr/bin/env bash
while (($#)); do
    if [[ "$1" == -o ]]; then
        shift
        printf 'plausible matching data\n' >"$1"
        break
    fi
    shift
done
exit 23
FAKE
chmod +x "$artifact_dir/failing-search"
diamond_bin="$artifact_dir/failing-search"
query=unused
db=unused
threads=4
run_serial=0
baseline_output=""
status=0
run_case fake-search-failure 192 || status=$?
[[ "$status" == 23 ]]
[[ ! -e "$artifact_dir/summary.tsv" ]]
printf 'PASS actual run_case preserves search status23 under ORcontext despite nonempty output\n'

cat >"$artifact_dir/successful-search" <<'FAKE'
#!/usr/bin/env bash
while (($#)); do
    if [[ "$1" == -o ]]; then
        shift
        printf 'plausible matching data\n' >"$1"
        break
    fi
    shift
done
exit 0
FAKE
chmod +x "$artifact_dir/successful-search"
diamond_bin="$artifact_dir/successful-search"
sha256_file() { return 29; }
status=0
run_case fake-hash-failure 192 || status=$?
[[ "$status" == 29 ]]
[[ ! -e "$artifact_dir/summary.tsv" ]]
printf 'PASS actual run_case preserves hash failure29 under ORcontext\n'
