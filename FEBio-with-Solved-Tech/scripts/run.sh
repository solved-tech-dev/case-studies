#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
case "${1:-}" in
    bp10) name=bp10-12 ;;
    tu01) name=tu01-8 ;;
    mi06) name=mi06-refined ;;
    *) echo "Usage: $0 {bp10|tu01|mi06} NEW_OUTPUT_DIRECTORY" >&2; exit 2 ;;
esac
[[ $# == 2 ]] || { echo 'Supply a new output directory.' >&2; exit 2; }
deck="$root/cases/$1/$name.feb"
plugin="$root/build/pangulu.so"
[[ -r "$deck" && -r "$plugin" ]] || { echo 'Prepare the case and build the plugin first.' >&2; exit 1; }
exe=$(command -v "${FEBIO_BIN:-febio4}")
mkdir -- "$2"
cd -- "$2"
"$exe" -noconfig -import "$plugin" -i "$deck" -o application.log -p physical.xplt >console.log 2>&1
grep -q 'N O R M A L   T E R M I N A T I O N' application.log
grep -q 'FEBIO_RESIDUAL ' console.log
test -s physical.xplt
