#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
require_slurm
require_env GPU_MPI_ROOT CUDA_ROOT LINGLONGAMG_ROOT METIS_ROOT PARMETIS_ROOT \
    BLAS_LIB LAPACK_LIB
require_file "${SOURCE_DIR}/CMakeLists.txt"
activate_mpi "${GPU_MPI_ROOT}"
export PATH="${CUDA_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${LINGLONGAMG_ROOT}/lib:${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"
require_command "${CMAKE_BIN}"
require_file "${LINGLONGAMG_ROOT}/lib/linglongamg.h"
require_file "${LINGLONGAMG_ROOT}/lib/liblinglongamg.so"
require_file "${BLAS_LIB}"
require_file "${LAPACK_LIB}"
require_file "${METIS_ROOT}/include/metis.h"
require_file "${PARMETIS_ROOT}/include/parmetis.h"
PARTITION_INCLUDE_FLAGS="${PARTITION_INCLUDE_FLAGS:--I${METIS_ROOT}/include -I${PARMETIS_ROOT}/include}"

export OMPI_CC=gcc
export OMPI_CXX=g++
export OMPI_MCA_opal_cuda_support=false
export OMP_NUM_THREADS=1

build_dir="${BUILD_ROOT}/linglongamg"
"${CMAKE_BIN}" -S "${SOURCE_DIR}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_EXE_LINKER_FLAGS=-no-pie \
    -DCMAKE_C_COMPILER="${GPU_MPI_ROOT}/bin/mpicc" \
    -DCMAKE_CXX_COMPILER="${GPU_MPI_ROOT}/bin/mpicxx" \
    -DCMAKE_CXX_FLAGS="${PARTITION_INCLUDE_FLAGS}" \
    -DOCP_USE_LINGLONGAMG=ON \
    -DLINGLONGAMG_DIR="${LINGLONGAMG_ROOT}" \
    -DLINGLONGAMG_INCLUDE_DIR="${LINGLONGAMG_ROOT}/lib" \
    -DCUDAToolkit_ROOT="${CUDA_ROOT}" \
    -DOCP_USE_OPENMP=ON \
    -DOCP_USE_MPI=ON \
    -DOCP_USE_METIS=OFF \
    -DOCP_USE_PARMETIS=OFF \
    -DBLAS_LIBRARIES="${BLAS_LIB}" \
    -DLAPACK_LIBRARIES="${LAPACK_LIB}"
"${CMAKE_BIN}" --build "${build_dir}" --parallel "${BUILD_JOBS:-8}"
