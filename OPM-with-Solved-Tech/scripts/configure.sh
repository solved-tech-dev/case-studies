#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: configure.sh --route=<linglongamg|xianzhikrylov> [OPTIONS]

Validate SDK dependencies and configure a CMake/Ninja build for the selected
OPM route. Fails fast when required backends or headers are missing; never
cross-links unrelated backends.

Options:
  --route=ROUTE         Required route: linglongamg or xianzhikrylov
  --source-dir=DIR      Prepared OPM source dir (default: work/source-<route>)
  --build-dir=DIR       CMake build dir (default: work/build-<route>)
  --prefix-path=PATHS   Extra CMAKE_PREFIX_PATH entries (semicolon-separated)
  --cuda-arch=ARCH      CUDA architecture (default: 120 for RTX 5090 / SM120)
  --skip-license-check  Skip checking SDK license files during preflight
  -h, --help            Show this help message and exit
EOF
}

route=""
source_dir=""
build_dir=""
prefix_path="${CMAKE_PREFIX_PATH:-}"
cuda_arch="${OPM_CUDA_ARCH:-120}"
skip_license=0

for arg in "$@"; do
    case "$arg" in
        --route=*) route="${arg#--route=}" ;;
        --source-dir=*) source_dir="${arg#--source-dir=}" ;;
        --build-dir=*) build_dir="${arg#--build-dir=}" ;;
        --prefix-path=*) prefix_path="${arg#--prefix-path=}" ;;
        --cuda-arch=*) cuda_arch="${arg#--cuda-arch=}" ;;
        --skip-license-check) skip_license=1 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $arg (run with --help)" ;;
    esac
done

validate_route "$route"

check_args=("--route=$route")
[[ $skip_license -eq 1 ]] && check_args+=("--skip-license-check")
bash "$package_dir/scripts/check-env.sh" "${check_args[@]}"

source_dir=${source_dir:-"$work_root/source-$route"}
build_dir=${build_dir:-"$work_root/build-$route"}

if [[ ! -d $source_dir ]]; then
    : "${OPM_SOURCE_ROOT:?Set OPM_SOURCE_ROOT or run scripts/prepare-source.sh --route=$route first.}"
    bash "$package_dir/scripts/prepare-source.sh" --route="$route" --dest="$source_dir"
fi

if [[ -r "$source_dir/.opm-solved-tech-route" ]]; then
    prepared_route=$(cat "$source_dir/.opm-solved-tech-route")
    [[ $prepared_route == "$route" ]] ||
        die "Source directory $source_dir was prepared for route=$prepared_route, not route=$route."
fi

cmake_args=(
    -S "$source_dir"
    -B "$build_dir"
    -G Ninja
    -DCMAKE_BUILD_TYPE=Release
    -DOPM_SOLVED_TECH_ROUTE="$route"
    -DUSE_GPU_BRIDGE=OFF
    -DUSE_HYPRE=OFF
    -DUSE_AMGX=OFF
    -DUSE_OPENCL=OFF
    -DBUILD_TESTING=OFF
    -DBUILD_FLOW_VARIANTS=OFF
    -DOPM_ENABLE_PYTHON=OFF
    -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc"
    -DCUDAToolkit_ROOT="$CUDA_HOME"
    -DCUDA_TOOLKIT_ROOT_DIR="$CUDA_HOME"
    -DCMAKE_CUDA_ARCHITECTURES="$cuda_arch"
    -DCUDA_ARCH_LIST="$cuda_arch"
)

if [[ -n $prefix_path ]]; then
    cmake_args+=(-DCMAKE_PREFIX_PATH="$prefix_path;$CUDA_HOME")
else
    cmake_args+=(-DCMAKE_PREFIX_PATH="$CUDA_HOME")
fi

ll_inc=""
ll_lib=""
xk_inc=""
xk_lib=""

case "$route" in
    linglongamg)
        ll_inc=$(find_linglong_header "$LINGLONGAMG_ROOT")
        ll_lib=$(find_linglong_library "$LINGLONGAMG_ROOT")
        cmake_args+=(-DLINGLONGAMG_INCLUDE_DIR="$ll_inc")
        ;;
    xianzhikrylov)
        xk_inc=$(find_xianzhi_include "$XIANZHIKRYLOV_ROOT")
        xk_lib=$(find_xianzhi_library "$XIANZHIKRYLOV_ROOT")
        cmake_args+=(
            -DXIANZHIKRYLOV_INCLUDE_DIR="$xk_inc"
            -DXIANZHIKRYLOV_LIBRARY="$xk_lib"
        )
        ;;
esac

mkdir -p -- "$build_dir"
cmake --fresh "${cmake_args[@]}"

cat > "$build_dir/route.env" <<EOF
ROUTE=$route
SOURCE_DIR=$source_dir
BUILD_DIR=$build_dir
CUDA_HOME=$CUDA_HOME
LINGLONGAMG_INCLUDE_DIR=$ll_inc
LINGLONGAMG_LIBRARY=$ll_lib
XIANZHIKRYLOV_INCLUDE_DIR=$xk_inc
XIANZHIKRYLOV_LIBRARY=$xk_lib
EOF

mkdir -p -- "$work_root"
printf '%s\n' "$route" > "$work_root/active-route"
info "Configured route=$route in $build_dir"
