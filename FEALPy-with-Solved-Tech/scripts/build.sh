#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
: "${PANGULU_ROOT:?Set PANGULU_ROOT to the installed PanguLU 26.46 SDK directory}"
cuda=${CUDA_HOME:-/usr/local/cuda}
mkdir -p "$root/build"
for precision in real complex; do
    flags=()
    if [[ $precision == complex ]]; then flags+=(-DFEALPY_COMPLEX); fi
    "${CXX:-g++}" -std=c++17 -O2 -fPIC -shared "${flags[@]}" \
        -I"$PANGULU_ROOT/include" -I"${CUDA_INCLUDE:-$cuda/include}" \
        -I"${CUDA_NVCC_INCLUDE:-$cuda/include}" "$root/adapter/gpu.cpp" \
        -L"$PANGULU_ROOT/lib" -L"${CUDA_LIB:-$cuda/lib64}" -Wl,--no-as-needed \
        -lpangulu -lpangulu_platform -l:libcudart.so.12 \
        -o "$root/build/libfealpy_pangulu_$precision.so"
done
