#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: check-env.sh --route=<linglongamg|xianzhikrylov|cpr> [OPTIONS]

Validate required toolchains, SDK headers/libraries, symbols, and license files
for the selected OPM GPU route without requiring unrelated backends.

Options:
  --route=ROUTE         Required route: linglongamg, xianzhikrylov, or cpr
  --require-opm         Also verify OPM_SOURCE_ROOT points to opm-simulators
  --skip-license-check  Skip checking SDK license files (for static/dry-run checks)
  -h, --help            Show this help message and exit
EOF
}

route=""
require_opm=0
skip_license=0

for arg in "$@"; do
    case "$arg" in
        --route=*) route="${arg#--route=}" ;;
        --require-opm) require_opm=1 ;;
        --skip-license-check) skip_license=1 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $arg (run with --help)" ;;
    esac
done

validate_route "$route"

if [[ $require_opm -eq 1 ]]; then
    : "${OPM_SOURCE_ROOT:?Set OPM_SOURCE_ROOT to your upstream opm-simulators 2026.04 source tree.}"
    [[ -r "$OPM_SOURCE_ROOT/CMakeLists.txt" && -r "$OPM_SOURCE_ROOT/opm/simulators/linalg/ISTLSolver.cpp" ]] ||
        die "OPM_SOURCE_ROOT does not look like an opm-simulators source tree: $OPM_SOURCE_ROOT"
fi

: "${CUDA_HOME:?Set CUDA_HOME to the installed CUDA toolkit root.}"
[[ -x "$CUDA_HOME/bin/nvcc" ]] || die "Cannot find nvcc under CUDA_HOME=$CUDA_HOME"

check_linglongamg() {
    : "${LINGLONGAMG_ROOT:?Set LINGLONGAMG_ROOT to the LinglongAMG SDK root for route=$route.}"
    local inc lib lic
    inc=$(find_linglong_header "$LINGLONGAMG_ROOT") ||
        die "Cannot find linglongamg.h under LINGLONGAMG_ROOT=$LINGLONGAMG_ROOT (or LINGLONGAMG_INCLUDE_DIR)."
    lib=$(find_linglong_library "$LINGLONGAMG_ROOT") ||
        die "Cannot find liblinglongamg.so under LINGLONGAMG_ROOT=$LINGLONGAMG_ROOT (or LINGLONGAMG_LIBRARY)."

    local syms
    syms=$(nm -D --defined-only "$lib" 2>/dev/null || true)
    if [[ -n $syms ]]; then
        for sym in \
            linglongamg_create_from_config \
            linglongamg_read_csr \
            linglongamg_csr_destroy \
            linglongamg_setup \
            linglongamg_reset \
            linglongamg_destroy; do
            grep -Eq "[[:space:]]${sym}(@@[^[:space:]]+)?$" <<< "$syms" ||
                die "LinglongAMG library $lib does not export required symbol: $sym"
        done
        if [[ $route == cpr ]]; then
            if ! grep -Eq "[[:space:]]linglongamg_apply_cycle(@@[^[:space:]]+)?$" <<< "$syms"; then
                grep -Eq "[[:space:]]linglongamg_solve(@@[^[:space:]]+)?$" <<< "$syms" &&
                    grep -Eq "[[:space:]]linglongamg_get_solve_info(@@[^[:space:]]+)?$" <<< "$syms" ||
                    die "CPR route requires linglongamg_apply_cycle or linglongamg_solve + linglongamg_get_solve_info in $lib"
            fi
        fi
    fi

    if [[ $skip_license -eq 0 ]]; then
        lic=$(find_license_file "${LINGLONGAMG_LICENSE:-}" "$LINGLONGAMG_ROOT") ||
            die "Cannot find LinglongAMG license file (set LINGLONGAMG_LICENSE or place license.key under LINGLONGAMG_ROOT)."
        info "LinglongAMG license: $lic"
    fi
    info "LinglongAMG include: $inc"
    info "LinglongAMG library: $lib"
}

check_xianzhikrylov() {
    : "${XIANZHIKRYLOV_ROOT:?Set XIANZHIKRYLOV_ROOT to the XianzhiKrylov SDK root for route=$route.}"
    local inc lib lic
    inc=$(find_xianzhi_include "$XIANZHIKRYLOV_ROOT") ||
        die "Cannot find XianzhiKrylov_CprComponents.h and XianzhiKrylov_MethodOptions.h under XIANZHIKRYLOV_ROOT=$XIANZHIKRYLOV_ROOT."
    lib=$(find_xianzhi_library "$XIANZHIKRYLOV_ROOT") ||
        die "Cannot find libXianzhiKrylov.so under XIANZHIKRYLOV_ROOT=$XIANZHIKRYLOV_ROOT (or XIANZHIKRYLOV_LIBRARY)."

    local syms
    syms=$(nm -D --defined-only "$lib" 2>/dev/null || true)
    if [[ -n $syms ]]; then
        for sym in \
            xk_internal_cprw_krylov_create \
            xk_internal_cprw_krylov_solve \
            xk_internal_cprw_krylov_destroy \
            xk_internal_cprw_ilu_create \
            xk_internal_cprw_ilu_prepare \
            xk_internal_cprw_ilu_apply \
            xk_internal_cprw_ilu_destroy; do
            grep -Eq "[[:space:]]${sym}(@@[^[:space:]]+)?$" <<< "$syms" ||
                die "XianzhiKrylov library $lib does not export required symbol: $sym"
        done
    fi

    if [[ $skip_license -eq 0 ]]; then
        lic=$(find_license_file "${XIANZHIKRYLOV_LICENSE:-}" "$XIANZHIKRYLOV_ROOT") ||
            die "Cannot find XianzhiKrylov license file (set XIANZHIKRYLOV_LICENSE or place license.key under XIANZHIKRYLOV_ROOT)."
        info "XianzhiKrylov license: $lic"
    fi
    info "XianzhiKrylov include: $inc"
    info "XianzhiKrylov library: $lib"
}

case "$route" in
    linglongamg)
        check_linglongamg
        info "Route isolation: XianzhiKrylov is NOT required and will NOT be linked."
        ;;
    xianzhikrylov)
        check_xianzhikrylov
        info "Route isolation: LinglongAMG is NOT required and will NOT be linked."
        ;;
    cpr)
        check_linglongamg
        check_xianzhikrylov
        info "Route cpr verified both LinglongAMG and XianzhiKrylov SDKs."
        ;;
esac

info "Environment check passed for route=$route (CUDA_HOME=$CUDA_HOME)."
