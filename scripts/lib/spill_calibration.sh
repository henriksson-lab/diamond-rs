#!/usr/bin/env bash
# run_case must retain every sample and enforce diagnostic/output integrity.
# Only RSS classifications may vary between otherwise equivalent processes.
select_stable_delayed_limit() {
    local preferred=$1 observed=$2 startup=$3
    local candidate=$preferred previous=-1 attempt repeat consecutive label
    local first_label=""
    SELECTED_LIMIT=0
    GATE_SAMPLE_1=""
    GATE_SAMPLE_2=""
    # Six candidates, at most three samples each. Start with the calibrated
    # interior and the known delayed point, then approach the startup side.
    for attempt in 1 2 3 4 5 6; do
        if ((attempt == 2)); then
            candidate=$observed
        elif ((attempt > 2)); then
            candidate=$(((candidate + startup) / 2))
        fi
        if ((candidate <= startup || candidate == previous)); then
            continue
        fi
        previous=$candidate
        consecutive=0
        for repeat in 1 2 3; do
            label="validation-${attempt}-${repeat}"
            run_case "$label" "$candidate" || return $?
            if [[ "$RUN_CLASS" == delayed-spill ]]; then
                consecutive=$((consecutive + 1))
                if ((consecutive == 1)); then
                    first_label=$label
                else
                    SELECTED_LIMIT=$candidate
                    GATE_SAMPLE_1=$first_label
                    GATE_SAMPLE_2=$label
                    return 0
                fi
            else
                consecutive=0
            fi
        done
    done
    echo 'error: no same-limit consecutive delayed-spill pair within bounded calibration' >&2
    return 1
}
