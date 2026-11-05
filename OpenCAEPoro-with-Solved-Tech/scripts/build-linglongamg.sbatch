#!/usr/bin/env bash
#SBATCH --job-name=ocp-build-linglongamg
#SBATCH --partition=a100
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=00:20:00

set -euo pipefail
SCRIPT_DIR="${OCP_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/scripts"
BUILD_JOBS=8 "${SCRIPT_DIR}/build-linglongamg.sh"
