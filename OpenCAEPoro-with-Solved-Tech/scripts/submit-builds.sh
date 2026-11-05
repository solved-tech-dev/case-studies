#!/usr/bin/env bash

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
mkdir -p "${LOG_DIR}" "${BUILD_ROOT}" "${RESULT_DIR}"
export OCP_PACKAGE_ROOT="${PACKAGE_ROOT}"

if [[ ! -f "${SOURCE_DIR}/CMakeLists.txt" ]]; then
    printf 'Prepared source tree not found: %s\n' "${SOURCE_DIR}" >&2
    printf 'Run %s/scripts/prepare-source.sh first.\n' "${PACKAGE_ROOT}" >&2
    exit 1
fi

sbatch --job-name=ocp-build-hypre-cpu \
    --chdir="${PACKAGE_ROOT}" \
    --export=ALL \
    --output="${LOG_DIR}/build-hypre-cpu-%j.out" \
    "${SCRIPT_DIR}/build-hypre-cpu.sbatch"
sbatch --job-name=ocp-build-hypre-gpu \
    --chdir="${PACKAGE_ROOT}" \
    --export=ALL \
    --output="${LOG_DIR}/build-hypre-gpu-%j.out" \
    "${SCRIPT_DIR}/build-hypre-gpu.sbatch"
sbatch --job-name=ocp-build-linglongamg \
    --chdir="${PACKAGE_ROOT}" \
    --export=ALL \
    --output="${LOG_DIR}/build-linglongamg-%j.out" \
    "${SCRIPT_DIR}/build-linglongamg.sbatch"
