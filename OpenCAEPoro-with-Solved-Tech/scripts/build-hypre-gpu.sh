#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
require_slurm
require_env GPU_MPI_ROOT CUDA_ROOT HYPRE_GPU_CONFIG METIS_ROOT PARMETIS_ROOT \
    BLAS_LIB LAPACK_LIB
require_file "${SOURCE_DIR}/CMakeLists.txt"
activate_mpi "${GPU_MPI_ROOT}"
export PATH="${CUDA_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"
require_command "${CMAKE_BIN}"
require_file "${BLAS_LIB}"
require_file "${LAPACK_LIB}"
require_file "${METIS_ROOT}/include/metis.h"
require_file "${PARMETIS_ROOT}/include/parmetis.h"
hypre_config="$(hypre_config_path "${HYPRE_GPU_CONFIG}")"
PARTITION_INCLUDE_FLAGS="${PARTITION_INCLUDE_FLAGS:--I${METIS_ROOT}/include -I${PARMETIS_ROOT}/include}"

export OMPI_CC=gcc
export OMPI_CXX=g++
export OMPI_MCA_opal_cuda_support=false
export OMP_NUM_THREADS=1

build_dir="${BUILD_ROOT}/hypre-gpu"
"${CMAKE_BIN}" -S "${SOURCE_DIR}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_EXE_LINKER_FLAGS=-no-pie \
    -DCMAKE_C_COMPILER="${GPU_MPI_ROOT}/bin/mpicc" \
    -DCMAKE_CXX_COMPILER="${GPU_MPI_ROOT}/bin/mpicxx" \
    -DCMAKE_CXX_FLAGS="${PARTITION_INCLUDE_FLAGS}" \
    -DOCP_USE_HYPRE=ON \
    -DHYPRE_DIR="${hypre_config}" \
    -DCUDAToolkit_ROOT="${CUDA_ROOT}" \
    -DOCP_USE_OPENMP=ON \
    -DOCP_USE_MPI=ON \
    -DOCP_USE_METIS=OFF \
    -DOCP_USE_PARMETIS=OFF \
    -DBLAS_LIBRARIES="${BLAS_LIB}" \
    -DLAPACK_LIBRARIES="${LAPACK_LIB}"
"${CMAKE_BIN}" --build "${build_dir}" --parallel "${BUILD_JOBS:-8}"
