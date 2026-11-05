#!/usr/bin/env bash
#SBATCH --job-name=ocp-build-hypre-cpu
#SBATCH --partition=cpu
#SBATCH --nodelist=cpu9462
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:20:00

set -euo pipefail
SCRIPT_DIR="${OCP_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/scripts"
BUILD_JOBS=16 "${SCRIPT_DIR}/build-hypre-cpu.sh"
