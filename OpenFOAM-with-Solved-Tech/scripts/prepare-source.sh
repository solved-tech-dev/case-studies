#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

[[ $# -le 1 ]] || die "Usage: $0 [source-work-directory]"
destination=${1:-"$package_dir/work/source"}
[[ ! -e $destination && ! -L $destination ]] || die "Path already exists: $destination"
mkdir -p -- "$(dirname -- "$destination")"
cp -R -- "$package_dir/overlay" "$destination"
printf 'Source directory: %s\n' "$destination"
