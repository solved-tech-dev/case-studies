#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
: "${FEBIO_SRC:?Set FEBIO_SRC to the matching FEBio source or SDK include root (contains FECore)}"
: "${FEBIO_LIB:?Set FEBIO_LIB to the directory containing libfecore.so}"
: "${PANGULU_ROOT:?Set PANGULU_ROOT to the PanguLU 26.46 installation}"
cuda=${CUDA_HOME:-/usr/local/cuda}
mkdir -p "$root/build"
"${CXX:-g++}" -std=c++17 -O2 -shared -fPIC -fopenmp -DLINUX \
    -I"$FEBIO_SRC" -I"$PANGULU_ROOT/include" -I"${CUDA_INCLUDE:-$cuda/include}" \
    -I"${CUDA_NVCC_INCLUDE:-$cuda/include}" "$root/adapter/febio_gpu.cpp" \
    -L"$FEBIO_LIB" -lfecore -L"${CUDA_LIB:-$cuda/lib64}" -l:libcudart.so.12 \
    -L"$PANGULU_ROOT/lib" -Wl,--no-as-needed -lpangulu -lpangulu_platform -ldl \
    -o "$root/build/pangulu.so"
