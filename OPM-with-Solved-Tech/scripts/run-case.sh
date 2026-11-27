#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: run-case.sh --route=<linglongamg|xianzhikrylov> --deck=FILE [OPTIONS]

Run an OPM Flow case using the selected route configuration and runtime library
dispatch.

Options:
  --route=ROUTE         Required route: linglongamg or xianzhikrylov
  --deck=FILE           Path to Eclipse .DATA deck
  --config=FILE         OPM linear solver JSON (default: examples/configs/<route>.json)
  --xk-config=FILE      XianzhiKrylov JSON config (route=xianzhikrylov; default: examples/configs/xk-bsr3-bicgstab-ilu0.json)
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
xk_config=""
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
        --xk-config=*) xk_config="${arg#--xk-config=}" ;;
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
fi

well_flag="--matrix-add-well-contributions=false"

case "$route" in
    linglongamg)
        : "${LINGLONGAMG_ROOT:?Set LINGLONGAMG_ROOT before running route=linglongamg.}"
        ll_lib=$(find_linglong_library "$LINGLONGAMG_ROOT") ||
            die "Cannot find liblinglongamg.so under LINGLONGAMG_ROOT=$LINGLONGAMG_ROOT."
        export LINGLONGAMG_OPM_LIBRARY="${LINGLONGAMG_OPM_LIBRARY:-$ll_lib}"
        export LINGLONGAMG_OPM_CONFIG="${LINGLONGAMG_OPM_CONFIG:-$package_dir/examples/configs/linglong-gmres-amg.json}"
        export LD_LIBRARY_PATH="$(dirname -- "$ll_lib")${CUDA_HOME:+:$CUDA_HOME/lib64}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
        ;;
    xianzhikrylov)
        : "${XIANZHIKRYLOV_ROOT:?Set XIANZHIKRYLOV_ROOT before running route=xianzhikrylov.}"
        xk_lib=$(find_xianzhi_library "$XIANZHIKRYLOV_ROOT") ||
            die "Cannot find libXianzhiKrylov.so under XIANZHIKRYLOV_ROOT=$XIANZHIKRYLOV_ROOT."
        xk_cfg=${xk_config:-"${XIANZHIKRYLOV_OPM_CONFIG:-$package_dir/examples/configs/xk-bsr3-bicgstab-ilu0.json}"}
        [[ -r $xk_cfg ]] || die "XianzhiKrylov config file not found: $xk_cfg"
        export XIANZHIKRYLOV_OPM_CONFIG="$xk_cfg"
        export LD_LIBRARY_PATH="$(dirname -- "$xk_lib")${CUDA_HOME:+:$CUDA_HOME/lib64}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
        well_flag="--matrix-add-well-contributions=true"
        ;;
esac

flow_args=(
    "$deck"
    "--output-dir=$output_dir"
    "--linear-solver=$config"
    "$well_flag"
    "--enable-opm-rst-file=false"
)
if [[ $no_output -eq 1 ]]; then
    flow_args+=("--enable-ecl-output=false")
fi

info "Running route=$route: $flow_bin ${flow_args[*]}"
exec "$flow_bin" "${flow_args[@]}"
