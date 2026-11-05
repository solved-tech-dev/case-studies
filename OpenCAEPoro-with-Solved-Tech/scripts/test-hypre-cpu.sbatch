#!/usr/bin/env bash
#SBATCH --job-name=ocp-spe10-hypre-cpu
#SBATCH --partition=cpu
#SBATCH --nodelist=cpu9462
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=96G
#SBATCH --time=00:20:00
#SBATCH --array=1-5%1

set -euo pipefail
SCRIPT_DIR="${OCP_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/scripts"
source "${SCRIPT_DIR}/common.sh"
factor="${SLURM_ARRAY_TASK_ID:-1}"
log_path="${OCP_LOG_PATH:-${LOG_DIR}/hypre-cpu-r${factor}-${SLURM_JOB_ID}.log}"
"${SCRIPT_DIR}/run-case.sh" hypre-cpu "${factor}" \
    "${BUILD_ROOT}/hypre-cpu/main/testOpenCAEPoro" "${log_path}"
