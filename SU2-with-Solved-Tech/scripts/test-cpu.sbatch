#!/usr/bin/env bash
#SBATCH --job-name=su2-cpu9462
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --nodelist=cpu9462
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=110G
#SBATCH --time=00:20:00
#SBATCH --array=1-5%1
set -euo pipefail
SCRIPT_DIR="${SU2_PACKAGE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/scripts"
cases=(Cavity Kenics IncTurbulentFlatPlate IncNACA Urban)
index="${SLURM_ARRAY_TASK_ID:-1}"
bash "${SCRIPT_DIR}/run-case.sh" cpu "${cases[$((index - 1))]}"
