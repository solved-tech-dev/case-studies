#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
UPSTREAM_DIR="${1:-${OCP_UPSTREAM_DIR:-}}"
TARGET_DIR="${2:-${OCP_SOURCE_DIR:-${PACKAGE_ROOT}/work/OpenCAEPoroX}}"

if [[ -z "${UPSTREAM_DIR}" ]]; then
    printf 'Usage: %s <OpenCAEPoroX-source-directory> [output-directory]\n' "$0" >&2
    exit 1
fi
if [[ ! -f "${UPSTREAM_DIR}/CMakeLists.txt" ]]; then
    printf 'Not an OpenCAEPoroX source tree: %s\n' "${UPSTREAM_DIR}" >&2
    exit 1
fi
if [[ -e "${TARGET_DIR}" ]]; then
    printf 'Target already exists; refusing to overwrite it: %s\n' "${TARGET_DIR}" >&2
    exit 1
fi

mkdir -p "$(dirname "${TARGET_DIR}")"
if command -v rsync >/dev/null 2>&1; then
    rsync -a \
        --exclude='.git' \
        --exclude='build*' \
        --exclude='lib' \
        "${UPSTREAM_DIR}/" "${TARGET_DIR}/"
else
    cp -a "${UPSTREAM_DIR}" "${TARGET_DIR}"
fi
cp -a "${PACKAGE_ROOT}/overlay/." "${TARGET_DIR}/"

printf 'Prepared source tree: %s\n' "${TARGET_DIR}"
