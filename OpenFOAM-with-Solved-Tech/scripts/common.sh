#!/usr/bin/env bash

package_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
case_template="$package_dir/examples/occDrivAerStaticMesh"

die() {
    printf '%s\n' "$*" >&2
    exit 1
}

require_openfoam() {
    [[ ${WM_PROJECT_VERSION:-} == v2412 ]] ||
        die "Source the OpenFOAM v2412 environment before running this script."
    [[ ${WM_OPTIONS:-} == *DPInt32* ]] ||
        die "This package requires the OpenFOAM DPInt32 build."
    export FOAM_FILEHANDLER=uncollated
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
}

require_backend() {
    case "$1" in
        linglongamg|xianzhikrylov) ;;
        *) die "Backend must be linglongamg or xianzhikrylov." ;;
    esac
}

sdk_environment() {
    : "${CUDA_HOME:?Set CUDA_HOME to the installed CUDA toolkit root.}"
    export CUDA_INCLUDE_DIR="${CUDA_INCLUDE_DIR:-$CUDA_HOME/include}"
    export CUDA_LIB_DIR="${CUDA_LIB_DIR:-$CUDA_HOME/lib64}"
    case "$1" in
        linglongamg)
            : "${LINGLONGAMG_SDK:?Set LINGLONGAMG_SDK to the installed SDK root.}"
            export LINGLONGAMG_INCLUDE_DIR="${LINGLONGAMG_INCLUDE_DIR:-$LINGLONGAMG_SDK/lib}"
            sdk_lib="$LINGLONGAMG_SDK/lib"
            sdk_license="${LINGLONGAMG_LICENSE:-$LINGLONGAMG_SDK/license.key}"
            executable=timedSimpleFoam
            adapter=libFoamPressure.so
            ;;
        xianzhikrylov)
            : "${XIANZHIKRYLOV_SDK:?Set XIANZHIKRYLOV_SDK to the installed SDK root.}"
            export XIANZHIKRYLOV_LIB_DIR="${XIANZHIKRYLOV_LIB_DIR:-$XIANZHIKRYLOV_SDK/lib}"
            sdk_lib="$XIANZHIKRYLOV_LIB_DIR"
            sdk_license="${XIANZHIKRYLOV_LICENSE:-$XIANZHIKRYLOV_SDK/license.key}"
            executable=benchSimpleFoam
            adapter=libBenchLinear.so
            ;;
    esac
    export LD_LIBRARY_PATH="$sdk_lib:$CUDA_LIB_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
}
