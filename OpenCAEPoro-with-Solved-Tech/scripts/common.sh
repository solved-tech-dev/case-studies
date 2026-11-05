#!/usr/bin/env bash

set -euo pipefail

PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="${OCP_WORK_DIR:-${PACKAGE_ROOT}/work}"
SOURCE_DIR="${OCP_SOURCE_DIR:-${WORK_DIR}/OpenCAEPoroX}"
BUILD_ROOT="${OCP_BUILD_ROOT:-${PACKAGE_ROOT}/build}"
RESULT_DIR="${OCP_RESULT_DIR:-${PACKAGE_ROOT}/results}"
LOG_DIR="${OCP_LOG_DIR:-${PACKAGE_ROOT}/logs}"

CPU_MPI_ROOT="${CPU_MPI_ROOT:-}"
GPU_MPI_ROOT="${GPU_MPI_ROOT:-}"
CMAKE_BIN="${CMAKE_BIN:-cmake}"
CUDA_ROOT="${CUDA_ROOT:-}"
METIS_ROOT="${METIS_ROOT:-}"
PARMETIS_ROOT="${PARMETIS_ROOT:-}"

HYPRE_CPU_CONFIG="${HYPRE_CPU_CONFIG:-}"
HYPRE_GPU_CONFIG="${HYPRE_GPU_CONFIG:-}"
LINGLONGAMG_ROOT="${LINGLONGAMG_ROOT:-}"

OPENBLAS_ROOT="${OPENBLAS_ROOT:-}"
OPENBLAS_LIB_DIR="${OPENBLAS_LIB_DIR:-${OPENBLAS_ROOT:+${OPENBLAS_ROOT}/lib}}"
BLAS_LIB="${BLAS_LIB:-${OPENBLAS_LIB_DIR:+${OPENBLAS_LIB_DIR}/libopenblas.so}}"
LAPACK_LIB="${LAPACK_LIB:-${OPENBLAS_LIB_DIR:+${OPENBLAS_LIB_DIR}/libopenblas.so}}"
CPU_BLAS_LIB="${CPU_BLAS_LIB:-}"
CPU_LAPACK_LIB="${CPU_LAPACK_LIB:-}"
PARTITION_INCLUDE_FLAGS="${PARTITION_INCLUDE_FLAGS:-}"

require_env()
{
    local name
    for name in "$@"; do
        if [[ -z "${!name:-}" ]]; then
            printf 'Required environment variable is not set: %s\n' "${name}" >&2
            exit 1
        fi
    done
}

require_slurm()
{
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        printf 'This script performs compute work and must run through sbatch/salloc/srun.\n' >&2
        exit 1
    fi
}

require_file()
{
    if [[ ! -e "$1" ]]; then
        printf 'Required file not found: %s\n' "$1" >&2
        exit 1
    fi
}

require_command()
{
    local command_or_path="$1"
    if [[ "${command_or_path}" == */* ]]; then
        require_file "${command_or_path}"
    elif ! command -v "${command_or_path}" >/dev/null 2>&1; then
        printf 'Required command not found: %s\n' "${command_or_path}" >&2
        exit 1
    fi
}

activate_mpi()
{
    local mpi_root="$1"
    require_file "${mpi_root}/bin/mpicc"
    require_file "${mpi_root}/bin/mpicxx"
    export PATH="${mpi_root}/bin:${PATH}"
    export LD_LIBRARY_PATH="${mpi_root}/lib:${LD_LIBRARY_PATH:-}"
}

hypre_config_path()
{
    local configured="$1"
    if [[ -f "${configured}/HYPREConfig.cmake" ]]; then
        printf '%s\n' "$configured"
        return
    fi
    if [[ -f "${configured}/lib/cmake/HYPRE/HYPREConfig.cmake" ]]; then
        printf '%s\n' "${configured}/lib/cmake/HYPRE"
        return
    fi
    printf 'HYPREConfig.cmake not found under %s\n' "$configured" >&2
    exit 1
}
