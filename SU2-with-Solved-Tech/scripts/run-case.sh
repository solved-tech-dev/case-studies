#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
require_slurm
method="${1:?method is required}"
case_name="${2:?case name is required}"
require_env XIANZHIKRYLOV_ROOT CUDA_ROOT
require_command "${PYTHON_BIN}"
require_file "${BUILD_ROOT}/xianzhikrylov/install/bin/SU2_CFD"
case "${method}" in
    cpu) backend=cpu ;;
    xianzhikrylov) backend=gpu ;;
    *) printf 'Unknown method: %s\n' "${method}" >&2; exit 1 ;;
esac
mkdir -p "${LOG_DIR}" "${RESULT_DIR}"
log_path="${3:-${LOG_DIR}/${method}-${case_name}-${SLURM_JOB_ID}.log}"
output_path="${RESULT_DIR}/${method}-${case_name}-${SLURM_JOB_ID}"
if [[ -e "${output_path}" || -L "${output_path}" ]]; then
    printf 'Output directory already exists: %s\n' "${output_path}" >&2
    exit 1
fi
mkdir -p "$(dirname "${log_path}")"
export LD_LIBRARY_PATH="${XIANZHIKRYLOV_ROOT}/lib:${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH:-}"
verification_args=(--verify)
if [[ "${SU2_TIMING:-0}" == "1" ]]; then
    verification_args=(--no-verify)
fi
set -o noclobber
srun --ntasks=1 --cpu-bind=cores "${PYTHON_BIN}" "${SCRIPT_DIR}/run_cases.py" \
    --backend "${backend}" --cases "${case_name}" \
    --binary "${BUILD_ROOT}/xianzhikrylov/install/bin/SU2_CFD" \
    --threads "${SU2_THREADS:-${SLURM_CPUS_PER_TASK:-8}}" \
    --outer-iterations "${SU2_ITER:-20}" "${verification_args[@]}" \
    --output "${output_path}" > "${log_path}" 2>&1
