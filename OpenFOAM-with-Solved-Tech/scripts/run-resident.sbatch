#!/usr/bin/env bash
#SBATCH --job-name=openfoam-resident
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
set -euo pipefail
# Submit from the package root after sourcing OpenFOAM and exporting SDK paths.
cd "${SLURM_SUBMIT_DIR:?}"
export RESIDENT_CPU_THREADS="${SLURM_CPUS_PER_TASK:-20}"
exec bash scripts/run-resident.sh "$@"
