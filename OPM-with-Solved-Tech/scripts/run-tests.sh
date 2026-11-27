#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: run-tests.sh [OPTIONS]

Run public repository verification checks. Without SDKs installed, use
--static-only to run all host static checks (syntax, configs, --help, invalid
route errors, no-SDK preflight, public tree hygiene, and optional patch-apply
dry-run).

Options:
  --static-only         Run host static checks only (no CUDA/SDK/OPM build)
  --opm-source=DIR      Optional clean opm-simulators tree for patch-apply test
  -h, --help            Show this help message and exit
EOF
}

static_only=0
opm_source="${OPM_SOURCE_ROOT:-}"

for arg in "$@"; do
    case "$arg" in
        --static-only) static_only=1 ;;
        --opm-source=*) opm_source="${arg#--opm-source=}" ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $arg (run with --help)" ;;
    esac
done

info "1/6 Checking shell script syntax..."
for s in "$package_dir"/scripts/*.sh; do
    bash -n "$s"
done

info "2/6 Validating example configs, manifests, and day-1 summary..."
python3 "$package_dir/tests/config/validate_configs.py"
python3 "$package_dir/tests/norne-day1/check-day1-summary.py"

info "3/6 Checking --help on all scripts..."
for s in "$package_dir"/scripts/*.sh; do
    bash "$s" --help >/dev/null
done

info "4/6 Checking invalid-route error handling and no-SDK preflight..."
for bad_route in invalid cpr; do
    if bash "$package_dir/scripts/configure.sh" --route="$bad_route" >/dev/null 2>&1; then
        die "configure.sh --route=$bad_route unexpectedly succeeded"
    fi
done
for r in linglongamg xianzhikrylov; do
    if env -u LINGLONGAMG_ROOT -u XIANZHIKRYLOV_ROOT -u CUDA_HOME \
        bash "$package_dir/scripts/check-env.sh" --route="$r" >/dev/null 2>&1; then
        die "check-env.sh --route=$r unexpectedly succeeded without SDK env vars"
    fi
done

info "5/6 Checking public tree hygiene..."
bash "$package_dir/scripts/verify-public-tree.sh"

if [[ -n $opm_source && -r "$opm_source/CMakeLists.txt" ]]; then
    info "6/6 Testing patch + route-isolated overlay preparation against $opm_source..."
    tmp_base=$(mktemp -d)
    trap 'rm -rf -- "$tmp_base"' EXIT
    for r in linglongamg xianzhikrylov; do
        OPM_SOLVED_TECH_WORK_ROOT="$tmp_base" \
            bash "$package_dir/scripts/prepare-source.sh" \
                --route="$r" --source-root="$opm_source" --dest="$tmp_base/src-$r" >/dev/null
        case "$r" in
            linglongamg)
                [[ -r "$tmp_base/src-$r/opm/simulators/linalg/LinglongAMGPreconditioner.hpp" ]] ||
                    die "Missing LinglongAMGPreconditioner.hpp for route=linglongamg"
                [[ ! -e "$tmp_base/src-$r/opm/simulators/linalg/NorneXianzhi.hpp" ]] ||
                    die "Unexpected NorneXianzhi.hpp in route=linglongamg"
                ;;
            xianzhikrylov)
                [[ -r "$tmp_base/src-$r/opm/simulators/linalg/NorneXianzhi.hpp" ]] ||
                    die "Missing NorneXianzhi.hpp for route=xianzhikrylov"
                [[ ! -e "$tmp_base/src-$r/opm/simulators/linalg/LinglongAMGPreconditioner.hpp" ]] ||
                    die "Unexpected LinglongAMGPreconditioner.hpp in route=xianzhikrylov"
                ;;
        esac
    done
    info "Patch + route-isolation dry-run passed for both routes."
else
    info "6/6 Skipped patch-apply test (OPM_SOURCE_ROOT / --opm-source not provided)."
fi

if [[ $static_only -eq 1 ]]; then
    info "Static verification completed (GPU smoke tests and Norne day-1 runs skipped by --static-only)."
    exit 0
fi

info "Checking GPU smoke / runtime prerequisites..."
for r in linglongamg xianzhikrylov; do
    if ! bash "$package_dir/scripts/check-env.sh" --route="$r"; then
        die "Runtime checks require CUDA + LinglongAMG + XianzhiKrylov SDKs (or run with --static-only)."
    fi
done
info "All runtime SDK prerequisites detected; build smoke binaries via tests/smoke/README.md."
