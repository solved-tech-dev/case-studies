#!/usr/bin/env bash
#SBATCH --job-name=occ-openfoam
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
exec bash scripts/run-case.sh "$@"
