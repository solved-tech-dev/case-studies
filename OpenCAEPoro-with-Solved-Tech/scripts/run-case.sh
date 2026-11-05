#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
require_slurm

method="${1:?method is required}"
factor="${2:?refine factor is required}"
executable="${3:?executable is required}"
log_path="${4:?log path is required}"

if [[ ! "${factor}" =~ ^[1-5]$ ]]; then
    printf 'Refine factor must be one of 1,2,3,4,5.\n' >&2
    exit 1
fi
require_file "${executable}"
require_file "${SOURCE_DIR}/data/spe10/spe10.data"
mkdir -p "$(dirname "${log_path}")" "${RESULT_DIR}"
status_path="${RESULT_DIR}/${method}-r${factor}-${SLURM_JOB_ID}.status"
meta_path="${RESULT_DIR}/${method}-r${factor}-${SLURM_JOB_ID}.meta"

case "${method}" in
    hypre-cpu)
        require_env CPU_MPI_ROOT
        activate_mpi "${CPU_MPI_ROOT}"
        export MPICH_CC=gcc
        export MPICH_CXX=g++
        export OMP_NUM_THREADS="${OMP_NUM_THREADS:-64}"
        export OCP_LINEAR_SOLVER=hypre
        ;;
    hypre-gpu)
        require_env GPU_MPI_ROOT CUDA_ROOT OPENBLAS_LIB_DIR
        activate_mpi "${GPU_MPI_ROOT}"
        export PATH="${CUDA_ROOT}/bin:${PATH}"
        export LD_LIBRARY_PATH="${OPENBLAS_LIB_DIR}:${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"
        export OMPI_CC=gcc
        export OMPI_CXX=g++
        export OMPI_MCA_opal_cuda_support=false
        export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
        export OCP_LINEAR_SOLVER=hypre
        ;;
    linglongamg)
        require_env GPU_MPI_ROOT CUDA_ROOT LINGLONGAMG_ROOT OPENBLAS_LIB_DIR
        activate_mpi "${GPU_MPI_ROOT}"
        export PATH="${CUDA_ROOT}/bin:${PATH}"
        export LD_LIBRARY_PATH="${LINGLONGAMG_ROOT}/lib:${OPENBLAS_LIB_DIR}:${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"
        export OMPI_CC=gcc
        export OMPI_CXX=g++
        export OMPI_MCA_opal_cuda_support=false
        export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
        export OCP_LINEAR_SOLVER=linglongamg
        export OCP_LINGLONGAMG_CONFIG="${OCP_LINGLONGAMG_CONFIG:-${PACKAGE_ROOT}/examples/linglongamg.json}"
        require_file "${OCP_LINGLONGAMG_CONFIG}"
        ;;
    *)
        printf 'Unknown method: %s\n' "${method}" >&2
        exit 1
        ;;
esac

export OCP_SPE10_REFINE_FACTOR="${factor}"
export OCP_SINGLE_SOLVE_BENCHMARK=1
export OCP_HYPRE_GMRES_RESTART="${OCP_HYPRE_GMRES_RESTART:-10}"
export OCP_HYPRE_INTERP_MAX_ELEMENTS="${OCP_HYPRE_INTERP_MAX_ELEMENTS:-4}"
export OCP_HYPRE_STRONG_THRESHOLD="${OCP_HYPRE_STRONG_THRESHOLD:-0.25}"
export OCP_HYPRE_AGGRESSIVE_LEVELS="${OCP_HYPRE_AGGRESSIVE_LEVELS:-0}"
export OCP_HYPRE_RAP2="${OCP_HYPRE_RAP2:-1}"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_STACKSIZE="${OMP_STACKSIZE:-256M}"
export OPENBLAS_NUM_THREADS=1

if [[ "${factor}" == "5" ]]; then
    export OCP_SPE10_R5_USE_R4_TIME_RATIO=1
else
    unset OCP_SPE10_R5_USE_R4_TIME_RATIO
fi

{
    printf 'method=%s\n' "${method}"
    printf 'factor=%s\n' "${factor}"
    printf 'executable=%s\n' "${executable}"
    printf 'hostname=%s\n' "$(hostname)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "${meta_path}"

cd "${SOURCE_DIR}"
started_epoch="$(date +%s)"
rc=1
set +e
srun --ntasks=1 --cpu-bind=cores \
    /usr/bin/time -f 'OCP_SIMULATOR_E2E_SECONDS=%e MAX_RSS_KIB=%M EXIT_STATUS=%x' \
    "${executable}" data/spe10/spe10.data \
    method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1 > "${log_path}" 2>&1
rc=$?
set -e
finished_epoch="$(date +%s)"

{
    printf 'exit=%s\n' "${rc}"
    printf 'wall_seconds=%s\n' "$((finished_epoch - started_epoch))"
    printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'log=%s\n' "${log_path}"
} >> "${meta_path}"
printf '%s\n' "${rc}" > "${status_path}"
exit "${rc}"
