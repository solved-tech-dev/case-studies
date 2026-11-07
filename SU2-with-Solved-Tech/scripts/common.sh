#!/usr/bin/env bash
set -euo pipefail
PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="${SU2_WORK_DIR:-${PACKAGE_ROOT}/work}"
SOURCE_DIR="${SU2_SOURCE_DIR:-${WORK_DIR}/SU2}"
BUILD_ROOT="${SU2_BUILD_ROOT:-${PACKAGE_ROOT}/build}"
RESULT_DIR="${SU2_RESULT_DIR:-${PACKAGE_ROOT}/results}"
LOG_DIR="${SU2_LOG_DIR:-${PACKAGE_ROOT}/logs}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
NINJA_BIN="${NINJA_BIN:-ninja}"
XIANZHIKRYLOV_ROOT="${XIANZHIKRYLOV_ROOT:-}"
CUDA_ROOT="${CUDA_ROOT:-}"

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

