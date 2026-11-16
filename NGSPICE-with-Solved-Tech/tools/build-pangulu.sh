#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
set -euo pipefail

if (( $# != 1 )); then
    printf 'Usage: %s BUILD_DIR\n' "$0" >&2
    exit 2
fi
: "${PANGULU_ROOT:?Set PANGULU_ROOT to your PanguLU 26.46 SDK directory}"
: "${CUDA_HOME:?Set CUDA_HOME to your CUDA 12.9 toolkit directory}"
root=$(cd "$(dirname "$0")/.." && pwd)
sdk=$(cd "$PANGULU_ROOT" && pwd)
cuda_headers=$(cd "${CUDA_RUNTIME_INCLUDE:-$CUDA_HOME/include}" && pwd)
cuda_lib=$(cd "${CUDA_RUNTIME_LIB_DIR:-$CUDA_HOME/lib64}" && pwd)
mkdir -p "$1"
build=$(cd "$1" && pwd)

include_flags=(-I"$sdk/include" -I"$cuda_headers")
if [[ -n ${CUDA_COMPILER_INCLUDE:-} ]]; then
    include_flags+=(-I"$CUDA_COMPILER_INCLUDE")
fi
"${CXX:-g++}" -std=c++17 -O3 -fopenmp "${include_flags[@]}" "$root/adapter/replay.cpp" \
    -L"$sdk/lib" -L"$cuda_lib" -Wl,-rpath,"$sdk/lib" -Wl,-rpath,"$cuda_lib" \
    -Wl,--no-as-needed -lpangulu -lpangulu_platform -Wl,--as-needed -l:libcudart.so.12 \
    -o "$build/pangulu-replay"
printf 'PanguLU replay executable: %s\n' "$build/pangulu-replay"
