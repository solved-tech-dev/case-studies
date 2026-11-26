#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: run-case.sh --route=<linglongamg|xianzhikrylov|cpr> --deck=FILE [OPTIONS]

Run an OPM Flow case using the selected route configuration and runtime library
dispatch.

Options:
  --route=ROUTE         Required route: linglongamg, xianzhikrylov, or cpr
  --deck=FILE           Path to Eclipse .DATA deck
  --config=FILE         Solver JSON config (default: examples/configs/<route>.json)
  --flow=BIN            Path to flow binary (default: work/build-<route>/bin/flow)
  --output-dir=DIR      Output directory (default: work/run-<route>)
  --no-output           Pass --enable-ecl-output=false (recommended for timing runs)
  --profile             Enable per-system stderr timing logs
  --true-residual       Enable independent true-residual verification
  -h, --help            Show this help message and exit
EOF
}

route=""
deck=""
config=""
flow_bin=""
output_dir=""
no_output=0
profile=0
true_residual=0

for arg in "$@"; do
    case "$arg" in
        --route=*) route="${arg#--route=}" ;;
        --deck=*) deck="${arg#--deck=}" ;;
        --config=*) config="${arg#--config=}" ;;
        --flow=*) flow_bin="${arg#--flow=}" ;;
        --output-dir=*) output_dir="${arg#--output-dir=}" ;;
        --no-output) no_output=1 ;;
        --profile) profile=1 ;;
        --true-residual) true_residual=1 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $arg (run with --help)" ;;
    esac
done

validate_route "$route"
[[ -n $deck && -r $deck ]] || die "Specify a readable deck via --deck=FILE."

config=${config:-"$package_dir/examples/configs/$route.json"}
[[ -r $config ]] || die "Config file not found: $config"

flow_bin=${flow_bin:-"$work_root/build-$route/bin/flow"}
[[ -x $flow_bin ]] || die "Flow executable not found: $flow_bin (run configure.sh and build.sh first)"

output_dir=${output_dir:-"$work_root/run-$route"}
mkdir -p -- "$output_dir"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OPM_SOLVED_TECH_ROUTE="$route"

if [[ $profile -eq 1 ]]; then
    export OPM_SOLVED_TECH_PROFILE=1
fi
if [[ $true_residual -eq 1 ]]; then
    export OPM_SOLVED_TECH_TRUE_RESIDUAL=1
    export NORNE_GPU_CPRW_AUDIT=1
fi

case "$route" in
    linglongamg)
        : "${LINGLONGAMG_ROOT:?Set LINGLONGAMG_ROOT before running route=linglongamg.}"
        ll_lib=$(find_linglong_library "$LINGLONGAMG_ROOT") ||
            die "Cannot find liblinglongamg.so under LINGLONGAMG_ROOT=$LINGLONGAMG_ROOT."
        export LINGLONGAMG_OPM_LIBRARY="${LINGLONGAMG_OPM_LIBRARY:-$ll_lib}"
        export LINGLONGAMG_OPM_CONFIG="${LINGLONGAMG_OPM_CONFIG:-$package_dir/examples/configs/linglong-gmres-amg.json}"
        export LD_LIBRARY_PATH="$(dirname -- "$ll_lib")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
        ;;
    xianzhikrylov)
        : "${XIANZHIKRYLOV_ROOT:?Set XIANZHIKRYLOV_ROOT before running route=xianzhikrylov.}"
        xk_lib=$(find_xianzhi_library "$XIANZHIKRYLOV_ROOT") ||
            die "Cannot find libXianzhiKrylov.so under XIANZHIKRYLOV_ROOT=$XIANZHIKRYLOV_ROOT."
        export LD_LIBRARY_PATH="$(dirname -- "$xk_lib")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
        ;;
    cpr)
        : "${LINGLONGAMG_ROOT:?Set LINGLONGAMG_ROOT before running route=cpr.}"
        : "${XIANZHIKRYLOV_ROOT:?Set XIANZHIKRYLOV_ROOT before running route=cpr.}"
        ll_lib=$(find_linglong_library "$LINGLONGAMG_ROOT") ||
            die "Cannot find liblinglongamg.so under LINGLONGAMG_ROOT=$LINGLONGAMG_ROOT."
        xk_lib=$(find_xianzhi_library "$XIANZHIKRYLOV_ROOT") ||
            die "Cannot find libXianzhiKrylov.so under XIANZHIKRYLOV_ROOT=$XIANZHIKRYLOV_ROOT."
        export NORNE_AMG_LIBRARY="${NORNE_AMG_LIBRARY:-$ll_lib}"
        export NORNE_AMG_CONFIG="${NORNE_AMG_CONFIG:-$package_dir/examples/configs/linglong-vcycle.json}"
        export NORNE_AMG_ALLOW_SINGLE_STEP="${NORNE_AMG_ALLOW_SINGLE_STEP:-1}"
        export LD_LIBRARY_PATH="$(dirname -- "$xk_lib"):$(dirname -- "$ll_lib")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
        ;;
esac

flow_args=(
    "$deck"
    "--output-dir=$output_dir"
    "--linear-solver-configuration-file=$config"
    "--matrix-add-well-contributions=false"
    "--enable-opm-rst-file=false"
)
if [[ $no_output -eq 1 ]]; then
    flow_args+=("--enable-ecl-output=false")
fi

info "Running route=$route: $flow_bin ${flow_args[*]}"
exec "$flow_bin" "${flow_args[@]}"
