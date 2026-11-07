#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sbatch "${SCRIPT_DIR}/build-xianzhikrylov.sbatch"
