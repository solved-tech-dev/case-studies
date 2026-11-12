#!/usr/bin/env bash
set -euo pipefail
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
build=$(cd -- "${PHENGLEI_BUILD_DIR:-"$package/build/linglong"}" && pwd)
sdk=$(cd -- "${LINGLONGAMG_ROOT:?Set LINGLONGAMG_ROOT}" && pwd)
test_binary="$build/PHengLEIv3d0/CFD/linglong_session_test"
audit_library="$build/lib/CFD/liblinglong_session_audit.so"
test -x "$test_binary" && test -f "$audit_library" || {
  printf 'Build with PHENGLEI_LINGLONG_TESTS=ON first.\n' >&2
  exit 2
}
output=${1:?Specify a new test output directory}
mkdir -p -- "$(dirname -- "$output")"
mkdir -- "$output"
output=$(cd -- "$output" && pwd)
python3 - "$package/examples/linglongamg_config.json" "$output/config.json" <<'PY'
import json
import sys
from pathlib import Path

config = json.loads(Path(sys.argv[1]).read_text())
config["solve"].update(tolerance=1e-10, max_iters=200, print_level=0)
Path(sys.argv[2]).write_text(json.dumps(config, indent=2) + "\n")
PY
export PHENGLEI_LINGLONG_LIBRARY="$audit_library"
export LINGLONG_AUDIT_REAL_LIBRARY="$sdk/lib/liblinglongamg.so"
export LINGLONG_AUDIT_FILE="$output/reuse-api.log"
unset LINGLONG_AUDIT_DENY_QUERY
"$test_binary" "$output/config.json" > "$output/reuse.log" 2>&1
export LINGLONG_AUDIT_FILE="$output/denied-api.log"
export LINGLONG_AUDIT_DENY_QUERY=1
"$test_binary" "$output/config.json" --deny-query > "$output/denied.log" 2>&1
printf 'Session tests passed. Logs in %s\n' "$output"
