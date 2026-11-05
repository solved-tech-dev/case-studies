#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
require_slurm
require_env CPU_MPI_ROOT HYPRE_CPU_CONFIG METIS_ROOT PARMETIS_ROOT \
    CPU_BLAS_LIB CPU_LAPACK_LIB
require_file "${SOURCE_DIR}/CMakeLists.txt"
activate_mpi "${CPU_MPI_ROOT}"
require_command "${CMAKE_BIN}"
hypre_config="$(hypre_config_path "${HYPRE_CPU_CONFIG}")"
require_file "${METIS_ROOT}/include/metis.h"
require_file "${PARMETIS_ROOT}/include/parmetis.h"
require_file "${CPU_BLAS_LIB}"
require_file "${CPU_LAPACK_LIB}"
PARTITION_INCLUDE_FLAGS="${PARTITION_INCLUDE_FLAGS:--I${METIS_ROOT}/include -I${PARMETIS_ROOT}/include}"

export MPICH_CC=gcc
export MPICH_CXX=g++
export OMP_NUM_THREADS=1

build_dir="${BUILD_ROOT}/hypre-cpu"
"${CMAKE_BIN}" -S "${SOURCE_DIR}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_EXE_LINKER_FLAGS=-no-pie \
    -DCMAKE_C_COMPILER="${CPU_MPI_ROOT}/bin/mpicc" \
    -DCMAKE_CXX_COMPILER="${CPU_MPI_ROOT}/bin/mpicxx" \
    -DCMAKE_CXX_FLAGS="${PARTITION_INCLUDE_FLAGS}" \
    -DOCP_USE_HYPRE=ON \
    -DHYPRE_DIR="${hypre_config}" \
    -DOCP_USE_OPENMP=ON \
    -DOCP_USE_MPI=ON \
    -DOCP_USE_METIS=OFF \
    -DOCP_USE_PARMETIS=OFF \
    -DBLAS_LIBRARIES="${CPU_BLAS_LIB}" \
    -DLAPACK_LIBRARIES="${CPU_LAPACK_LIB}"
"${CMAKE_BIN}" --build "${build_dir}" --parallel "${BUILD_JOBS:-16}"
