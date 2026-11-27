#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: build.sh [OPTIONS]

Build OPM Flow for the configured route and verify backend link isolation.

Options:
  --route=ROUTE         Route to build (default: work/active-route from configure.sh)
  --build-dir=DIR       CMake build dir (default: work/build-<route>)
  --target=TARGET       CMake target (default: flow)
  --jobs=N              Parallel build jobs (default: BUILD_JOBS or nproc)
  -h, --help            Show this help message and exit
EOF
}

route=""
build_dir=""
target="flow"
jobs="${BUILD_JOBS:-$(nproc 2>/dev/null || echo 4)}"

for arg in "$@"; do
    case "$arg" in
        --route=*) route="${arg#--route=}" ;;
        --build-dir=*) build_dir="${arg#--build-dir=}" ;;
        --target=*) target="${arg#--target=}" ;;
        --jobs=*) jobs="${arg#--jobs=}" ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $arg (run with --help)" ;;
    esac
done

if [[ -z $route ]]; then
    [[ -r "$work_root/active-route" ]] ||
        die "No active route recorded. Run scripts/configure.sh --route=<linglongamg|xianzhikrylov> first."
    route=$(cat "$work_root/active-route")
fi
validate_route "$route"

build_dir=${build_dir:-"$work_root/build-$route"}
[[ -d $build_dir && -r "$build_dir/route.env" ]] ||
    die "Build directory $build_dir is not configured. Run scripts/configure.sh --route=$route first."

# shellcheck disable=SC1091
source "$build_dir/route.env"
[[ ${ROUTE:-} == "$route" ]] ||
    die "Mismatch: $build_dir/route.env records ROUTE=${ROUTE:-unset}, expected $route."

cmake --build "$build_dir" --target "$target" --parallel "$jobs"

flow_bin="$build_dir/bin/flow"
if [[ -x $flow_bin ]] && command -v ldd >/dev/null 2>&1; then
    ldd_out=$(ldd "$flow_bin" 2>&1 || true)
    case "$route" in
        linglongamg)
            if grep -q "libXianzhiKrylov" <<< "$ldd_out"; then
                die "Isolation violation: linglongamg binary links libXianzhiKrylov!"
            fi
            info "Verified $flow_bin does NOT link libXianzhiKrylov."
            ;;
        xianzhikrylov)
            if grep -q "liblinglongamg" <<< "$ldd_out"; then
                die "Isolation violation: xianzhikrylov binary links liblinglongamg!"
            fi
            grep -q "libXianzhiKrylov" <<< "$ldd_out" ||
                die "Expected $flow_bin to link libXianzhiKrylov."
            info "Verified $flow_bin links libXianzhiKrylov and does NOT link liblinglongamg."
            ;;
    esac
fi

info "Build completed for route=$route ($flow_bin)"
