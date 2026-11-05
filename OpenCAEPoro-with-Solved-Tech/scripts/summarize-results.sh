#!/usr/bin/env bash

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

field()
{
    local file="$1"
    local key="$2"
    awk -v key="${key}" '
        {
            prefix = key "="
            for (i = 1; i <= NF; ++i) {
                if (index($i, prefix) == 1) {
                    split($i, parts, "=")
                    value = parts[2]
                }
            }
        }
        END { if (value != "") print value }
    ' "${file}"
}

printf 'method,r,exit,wall_s,e2e_s,rows,nnz,matrix_generation_s,common_s,solver_preprocess_s,setup_s,solve_s,iterations,true_relative_residual\n'
shopt -s nullglob
for meta in "${RESULT_DIR}"/*.meta; do
    method=
    factor=
    exit_code=
    wall_seconds=
    log=
    source "${meta}"
    log_path="${log}"
    if [[ -z "${log_path}" || ! -f "${log_path}" ]]; then
        continue
    fi
    rows="$(field "${log_path}" rows)"
    nnz="$(field "${log_path}" nnz)"
    setup_seconds="$(field "${log_path}" setup_seconds)"
    solve_seconds="$(field "${log_path}" solve_seconds)"
    iterations="$(field "${log_path}" iterations)"
    residual="$(field "${log_path}" relative_residual)"
    common_seconds="$(field "${log_path}" common_total_seconds)"
    matrix_seconds="$(field "${log_path}" matrix_generation_seconds)"
    preprocess_seconds="$(field "${log_path}" solver_preprocess_seconds)"
    e2e_seconds="$(field "${log_path}" OCP_SIMULATOR_E2E_SECONDS)"
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "${method}" "${factor}" "${exit_code}" "${wall_seconds}" \
        "${e2e_seconds}" "${rows}" "${nnz}" "${matrix_seconds}" \
        "${common_seconds}" "${preprocess_seconds}" "${setup_seconds}" \
        "${solve_seconds}" "${iterations}" "${residual}"
done
