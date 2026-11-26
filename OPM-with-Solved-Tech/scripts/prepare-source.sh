#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: prepare-source.sh --route=<linglongamg|xianzhikrylov|cpr> [OPTIONS]

Prepare a patched opm-simulators working tree in work/source-<route> by copying
OPM_SOURCE_ROOT, applying overlay/patches/opm-2026.04-solved-tech.patch, and
installing only the route-specific overlay source files.

Options:
  --route=ROUTE         Required route: linglongamg, xianzhikrylov, or cpr
  --source-root=DIR     Upstream opm-simulators 2026.04 tree (default: $OPM_SOURCE_ROOT)
  --dest=DIR            Destination directory (default: work/source-<route>)
  --force               Remove existing destination directory before preparing
  -h, --help            Show this help message and exit
EOF
}

route=""
source_root="${OPM_SOURCE_ROOT:-}"
dest=""
force=0

for arg in "$@"; do
    case "$arg" in
        --route=*) route="${arg#--route=}" ;;
        --source-root=*) source_root="${arg#--source-root=}" ;;
        --dest=*) dest="${arg#--dest=}" ;;
        --force) force=1 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $arg (run with --help)" ;;
    esac
done

validate_route "$route"
[[ -n $source_root ]] || die "Set OPM_SOURCE_ROOT or pass --source-root=DIR."
[[ -r "$source_root/CMakeLists.txt" && -r "$source_root/opm/simulators/linalg/ISTLSolver.cpp" ]] ||
    die "Invalid opm-simulators source tree: $source_root"

dest=${dest:-"$work_root/source-$route"}
if [[ -e $dest ]]; then
    if [[ $force -eq 1 ]]; then
        rm -rf -- "$dest"
    else
        die "Destination already exists: $dest (pass --force to replace)"
    fi
fi

mkdir -p -- "$(dirname -- "$dest")"
cp -a -- "$source_root" "$dest"

patch_file="$package_dir/overlay/patches/opm-2026.04-solved-tech.patch"
[[ -r $patch_file ]] || die "Missing patch file: $patch_file"

if command -v git >/dev/null 2>&1 && git -C "$dest" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git -C "$dest" apply --check "$patch_file" || die "Patch check failed on $dest"
    git -C "$dest" apply "$patch_file"
else
    patch -d "$dest" -p1 --dry-run < "$patch_file" >/dev/null || die "Patch dry-run failed on $dest"
    patch -d "$dest" -p1 < "$patch_file"
fi

case "$route" in
    linglongamg)
        cp -R -- "$package_dir/overlay/linglongamg/opm/." "$dest/opm/"
        ;;
    xianzhikrylov)
        cp -R -- "$package_dir/overlay/common/opm/." "$dest/opm/"
        cp -R -- "$package_dir/overlay/xianzhikrylov/opm/." "$dest/opm/"
        ;;
    cpr)
        cp -R -- "$package_dir/overlay/common/opm/." "$dest/opm/"
        cp -R -- "$package_dir/overlay/cpr/opm/." "$dest/opm/"
        cp -- "$package_dir/overlay/linglongamg/opm/simulators/linalg/linglongamg_sdk.h" \
            "$dest/opm/simulators/linalg/linglongamg_sdk.h"
        ;;
esac

printf '%s\n' "$route" > "$dest/.opm-solved-tech-route"
info "Prepared $route source tree at: $dest"
