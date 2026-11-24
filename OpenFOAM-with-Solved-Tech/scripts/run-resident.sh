#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/resident-common.sh"
require_openfoam
sdk_environment
resident_environment
exec python3 "$package_dir/scripts/run_resident.py" "$@"
