#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
UPSTREAM_DIR="${1:-${SU2_UPSTREAM_DIR:-}}"
TARGET_DIR="${2:-${SOURCE_DIR}}"
if [[ -z "${UPSTREAM_DIR}" ]]; then
    printf 'Usage: %s <SU2-8.5.0-source-directory> [output-directory]\n' "$0" >&2
    exit 1
fi
require_file "${UPSTREAM_DIR}/meson.build"
if ! grep -q "version:.*8.5.0" "${UPSTREAM_DIR}/meson.build"; then
    printf 'Expected SU2 8.5.0 source: %s\n' "${UPSTREAM_DIR}" >&2
    exit 1
fi
require_command "${PYTHON_BIN}"
"${PYTHON_BIN}" - "${UPSTREAM_DIR}" "${TARGET_DIR}" <<'PY'
from pathlib import Path
import sys
source, target = (Path(value).resolve() for value in sys.argv[1:])
if target == source or source in target.parents:
    raise SystemExit('The prepared-source output must be outside the upstream source tree.')
PY
mkdir -p "$(dirname "${TARGET_DIR}")"
lock_dir="${TARGET_DIR}.prepare.lock"
if ! mkdir "${lock_dir}"; then
    printf 'Another source preparation owns %s\n' "${lock_dir}" >&2
    exit 1
fi
staging=''
trap '
    if [[ -n "${staging}" ]]; then rm -rf "${staging}"; fi
    rmdir "${lock_dir}"
' EXIT
staging="$(mktemp -d "$(dirname "${TARGET_DIR}")/.su2-prepare.XXXXXX")"
if [[ -e "${TARGET_DIR}" || -L "${TARGET_DIR}" ]]; then
    printf 'Target already exists: %s\n' "${TARGET_DIR}" >&2
    exit 1
fi
if command -v rsync >/dev/null 2>&1; then
    rsync -aL --exclude='.git' --exclude='/build*/' --exclude='__pycache__/' \
        "${UPSTREAM_DIR}/" "${staging}/"
else
    cp -RL "${UPSTREAM_DIR}/." "${staging}/"
    rm -rf "${staging}/.git"
fi
cp -R "${PACKAGE_ROOT}/overlay/." "${staging}/"
mv "${staging}" "${TARGET_DIR}"
printf 'Prepared source tree: %s\n' "${TARGET_DIR}"
