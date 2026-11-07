#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
require_slurm
require_env XIANZHIKRYLOV_ROOT CUDA_ROOT
require_file "${SOURCE_DIR}/meson.build"
require_command "${PYTHON_BIN}"
require_command "${NINJA_BIN}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/build_xianzhikrylov.py" \
    --source-dir "${SOURCE_DIR}" --workspace "${BUILD_ROOT}/xianzhikrylov" \
    --sdk-root "${XIANZHIKRYLOV_ROOT}" --cuda-root "${CUDA_ROOT}" \
    --cc "${CC:-gcc}" --cxx "${CXX:-g++}" --ninja "${NINJA_BIN}" \
    --jobs "${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-8}}"
