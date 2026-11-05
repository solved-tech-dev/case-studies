#!/usr/bin/env bash
#SBATCH --job-name=ocp-spe10-hypre-gpu
#SBATCH --partition=a100
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=00:20:00
#SBATCH --array=1-5%1

set -euo pipefail
SCRIPT_DIR="${OCP_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/scripts"
source "${SCRIPT_DIR}/common.sh"
factor="${SLURM_ARRAY_TASK_ID:-1}"
log_path="${OCP_LOG_PATH:-${LOG_DIR}/hypre-gpu-r${factor}-${SLURM_JOB_ID}.log}"
"${SCRIPT_DIR}/run-case.sh" hypre-gpu "${factor}" \
    "${BUILD_ROOT}/hypre-gpu/main/testOpenCAEPoro" "${log_path}"
