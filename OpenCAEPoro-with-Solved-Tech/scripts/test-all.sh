#!/usr/bin/env bash

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
mkdir -p "${LOG_DIR}" "${RESULT_DIR}"
export OCP_PACKAGE_ROOT="${PACKAGE_ROOT}"

for executable in \
    "${BUILD_ROOT}/hypre-cpu/main/testOpenCAEPoro" \
    "${BUILD_ROOT}/hypre-gpu/main/testOpenCAEPoro" \
    "${BUILD_ROOT}/linglongamg/main/testOpenCAEPoro"; do
    if [[ ! -x "${executable}" ]]; then
        printf 'Executable not found or not executable: %s\n' "${executable}" >&2
        exit 1
    fi
done

sbatch --job-name=ocp-spe10-hypre-cpu-r1-r4 \
    --array=1-4 --chdir="${PACKAGE_ROOT}" --export=ALL \
    --output="${LOG_DIR}/hypre-cpu-r1-r4-%A_%a.out" \
    "${SCRIPT_DIR}/test-hypre-cpu.sbatch"
sbatch --job-name=ocp-spe10-hypre-cpu-r5 \
    --array=5 --cpus-per-task="${OCP_CPU_R5_THREADS:-8}" --mem=116G \
    --chdir="${PACKAGE_ROOT}" \
    --export="ALL,OMP_NUM_THREADS=${OCP_CPU_R5_THREADS:-8}" \
    --output="${LOG_DIR}/hypre-cpu-r5-%A_%a.out" \
    "${SCRIPT_DIR}/test-hypre-cpu.sbatch"
sbatch --chdir="${PACKAGE_ROOT}" --export=ALL \
    --output="${LOG_DIR}/hypre-gpu-%A_%a.out" \
    "${SCRIPT_DIR}/test-hypre-gpu.sbatch"
sbatch --chdir="${PACKAGE_ROOT}" --export=ALL \
    --output="${LOG_DIR}/linglongamg-%A_%a.out" \
    "${SCRIPT_DIR}/test-linglongamg.sbatch"

printf 'Submitted HYPRE CPU, HYPRE GPU and LinglongAMG r1-r5 tests.\n'
