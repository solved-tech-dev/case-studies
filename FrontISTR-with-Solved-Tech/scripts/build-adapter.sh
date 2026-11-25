#!/usr/bin/env bash
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
: "${PANGULU_ROOT:?Set PANGULU_ROOT to the directory containing include and lib}"
cuda=${CUDA_HOME:-/usr/local/cuda}
mkdir -p "$root/build"
"${CXX:-g++}" -O2 -std=c++17 -fPIC -shared -fopenmp \
    -I"$PANGULU_ROOT/include" -I"${CUDA_INCLUDE:-$cuda/include}" -I"${CUDA_NVCC_INCLUDE:-$cuda/include}" \
    "$root/adapter/bridge.cpp" -L"$PANGULU_ROOT/lib" -L"${CUDA_LIB:-$cuda/lib64}" \
    -Wl,-soname,libfrontistr_pangulu.so -Wl,--no-as-needed -lpangulu -lpangulu_platform -l:libcudart.so.12 -ldl \
    -o "$root/build/libfrontistr_pangulu.so"
