#!/usr/bin/env bash
#SBATCH --job-name=su2-build-xianzhikrylov
#SBATCH --partition=a100
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=00:30:00
set -euo pipefail
SCRIPT_DIR="${SU2_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/scripts"
bash "${SCRIPT_DIR}/build-xianzhikrylov.sh"
