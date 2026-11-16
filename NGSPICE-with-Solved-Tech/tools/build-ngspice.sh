#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
set -euo pipefail

if (( $# < 1 || $# > 2 )); then
    printf 'Usage: %s BUILD_DIR [LOCAL_UPSTREAM_REPOSITORY]\n' "$0" >&2
    exit 2
fi
root=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$1"
build=$(cd "$1" && pwd)
read -r repository tag commit < <(python3 -c '
import json, sys
with open(sys.argv[1]) as source:
    lock = json.load(source)
print(lock["repository"], lock["tag"], lock["commit"])
' "$root/upstream.lock.json")
source="$build/ngspice-src"
if [[ -e "$source" ]]; then
    printf 'Use a new build directory. Source already exists at %s\n' "$source" >&2
    exit 2
fi
if (( $# == 2 )); then
    mkdir "$source"
    git -C "$2" archive "$commit" | tar -x -C "$source"
else
    git clone --quiet --depth 1 --branch "$tag" "$repository" "$source"
    if [[ $(git -C "$source" rev-parse HEAD) != "$commit" ]]; then
        printf 'The upstream tag does not match upstream.lock.json\n' >&2
        exit 2
    fi
fi

cd "$source"
git apply "$root/integration/ngspice/capture.patch"
cp "$root/integration/ngspice/capture.c" src/maths/ni/solvedtech-capture.c
./autogen.sh
mkdir "$build/ngspice-build"
cd "$build/ngspice-build"
CFLAGS=${CFLAGS:--O2 -std=gnu11} "$source/configure" \
    --without-x --without-readline --disable-debug --disable-xspice --disable-cider \
    --disable-osdi --with-ngshared=no --enable-klu
make -j "${JOBS:-2}"
printf 'ngspice capture executable: %s\n' "$build/ngspice-build/src/ngspice"
