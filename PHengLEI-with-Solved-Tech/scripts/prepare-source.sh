#!/usr/bin/env bash
set -euo pipefail
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source_dir=$(cd -- "${PHENGLEI_UPSTREAM_DIR:?Set PHENGLEI_UPSTREAM_DIR to the authorized source tree}" && pwd)
target=${PHENGLEI_SOURCE_DIR:-"$package/work/PHengLEI"}
commit=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["commit"])' "$package/upstream.lock.json")

mkdir -p -- "$(dirname -- "$target")"
mkdir -- "$target"
created=$target
trap 'test -z "$created" || rm -rf -- "$created"' EXIT
if git -C "$source_dir" cat-file -e "$commit^{commit}" 2>/dev/null; then
  git -C "$source_dir" archive "$commit" | tar -xf - -C "$target"
else
  (cd -- "$source_dir" && tar --exclude='./.git' -cf - .) | tar -xf - -C "$target"
fi
python3 - "$target" "$package/upstream.lock.json" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

root = Path(sys.argv[1])
with open(sys.argv[2]) as file:
    expected = json.load(file)["base_file_sha256"]
for name, digest in expected.items():
    path = root / name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        sys.exit(f"Source does not match the pinned baseline: {name}")
PY
(
  cd -- "$target"
  # The case-study package may itself be in a Git checkout.
  export GIT_CEILING_DIRECTORIES="$(dirname -- "$(pwd -P)")"
  git apply --check "$package/patches/0001-linglong-pressure.patch"
  git apply "$package/patches/0001-linglong-pressure.patch"
)
cp -R "$package/overlay/." "$target/"
eigen_rel=PHengLEI/3rdparty/eigen
eigen=${EIGEN_ROOT:-"$source_dir/$eigen_rel"}
if [[ ! -f "$target/$eigen_rel/Eigen/Eigen" && -f "$eigen/Eigen/Eigen" ]]; then
  mkdir -p -- "$target/$eigen_rel"
  cp -R "$eigen/." "$target/$eigen_rel/"
elif [[ -n "${EIGEN_ROOT:-}" && ! -f "$eigen/Eigen/Eigen" ]]; then
  printf 'EIGEN_ROOT must contain Eigen/Eigen: %s\n' "$eigen" >&2
  exit 2
fi
created=
printf 'Prepared %s\n' "$target"
