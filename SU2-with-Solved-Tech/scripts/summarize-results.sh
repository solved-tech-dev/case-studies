#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"
"${PYTHON_BIN}" - "${RESULT_DIR}" <<'PY'
import csv
from pathlib import Path
import sys

fields = ('case', 'backend', 'threads', 'outer_iterations', 'returncode',
          'application_wall_seconds', 'whole_linear_seconds', 'gpu_execute_seconds',
          'gpu_prepare_seconds', 'gpu_solve_seconds', 'linear_iterations',
          'mean_reported_residual')
writer = csv.DictWriter(sys.stdout, fieldnames=fields, extrasaction='ignore')
writer.writeheader()
for path in sorted(Path(sys.argv[1]).glob('*/summary.csv')):
    with path.open(newline='') as stream:
        writer.writerows(csv.DictReader(stream))
PY
