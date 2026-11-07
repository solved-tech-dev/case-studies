#!/usr/bin/env bash
#SBATCH --job-name=su2-xianzhikrylov
#SBATCH --partition=a100
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=00:20:00
#SBATCH --array=1-5%1
set -euo pipefail
SCRIPT_DIR="${SU2_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/scripts"
cases=(Cavity Kenics IncTurbulentFlatPlate IncNACA Urban)
index="${SLURM_ARRAY_TASK_ID:-1}"
bash "${SCRIPT_DIR}/run-case.sh" xianzhikrylov "${cases[$((index - 1))]}"
