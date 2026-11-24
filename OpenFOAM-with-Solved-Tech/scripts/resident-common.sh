#!/usr/bin/env bash

package_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)

die() { printf '%s\n' "$*" >&2; exit 1; }

require_openfoam() {
    [[ ${WM_PROJECT_VERSION:-} == v2412 && ${WM_OPTIONS:-} == *DPInt32* ]] ||
        die "Source an OpenFOAM v2412 DPInt32 environment."
    export FOAM_FILEHANDLER=uncollated
}

sdk_environment() {
    : "${PREPARATION_SDK:?Set PREPARATION_SDK to the resident-r1 preparation directory.}"
    : "${XIANZHIKRYLOV_SDK:?Set XIANZHIKRYLOV_SDK to the matching resident-r1 core directory.}"
    : "${RESIDENT_CUDA:?Set RESIDENT_CUDA to the CUDA 13.1 toolkit directory.}"
    export PREPARATION_SDK="$(cd -- "$PREPARATION_SDK" && pwd -P)"
    export XIANZHIKRYLOV_SDK="$(cd -- "$XIANZHIKRYLOV_SDK" && pwd -P)"
    export RESIDENT_CUDA="$(cd -- "$RESIDENT_CUDA" && pwd -P)"
    [[ -r $PREPARATION_SDK/include/XianzhiPreparation2647.hpp &&
       -r $PREPARATION_SDK/lib/libLinearSystemPreparation.so &&
       -r $XIANZHIKRYLOV_SDK/include/XianzhiKrylovAPI.h &&
       -r $XIANZHIKRYLOV_SDK/lib/libXianzhiKrylov.so &&
       -r $RESIDENT_CUDA/include/cuda_runtime_api.h ]] || die "Incomplete resident-r1 SDK/toolkit."
    export RESIDENT_BUILD="${RESIDENT_BUILD:-$package_dir/work/resident/build}"
    mkdir -p -- "$RESIDENT_BUILD"
    export RESIDENT_BUILD="$(cd -- "$RESIDENT_BUILD" && pwd -P)"
    export LD_LIBRARY_PATH="$PREPARATION_SDK/lib:$XIANZHIKRYLOV_SDK/lib:$RESIDENT_CUDA/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
}

resident_environment() {
    : "${LINGLONGAMG_SDK:?Set LINGLONGAMG_SDK to the matching LinglongAMG installation.}"
    export LINGLONGAMG_SDK="$(cd -- "$LINGLONGAMG_SDK" && pwd -P)"
    [[ -r $LINGLONGAMG_SDK/lib/liblinglongamg.so ]] || die "Missing lib/liblinglongamg.so."
    export LD_LIBRARY_PATH="$LINGLONGAMG_SDK/lib:$LD_LIBRARY_PATH"
    export OMP_NUM_THREADS="${RESIDENT_CPU_THREADS:-20}" OPENBLAS_NUM_THREADS=1
    export XK_PRIVATE_DEVICE_VECTOR_BORROW=1 XK_PRIVATE_BICGSTAB_COMPACT_VECTORS=1
    export XK_PRIVATE_BICGSTAB_TRUE_RESTART=1 XK_PRIVATE_AMG_RELEASE_AFTER_SOLVE=1
    export XK_PRIVATE_GS_CSR_ONLY=1 XK_PRIVATE_GS_RELEASE_VECTORS=1
    export XK_PRIVATE_GS_COMPACT_PARTITIONS=1 XK_PRIVATE_GS_RELEASE_WORKSPACE=1
    export XK_PRIVATE_RESIDENT_GS_DIRECT=1 XK_PRIVATE_RESIDENT_GS_BATCH=1
    export RESIDENT_GPU_CELL_GEOMETRY=1 RESIDENT_GPU_BOUNDARY_GEOMETRY=1
    export RESIDENT_GPU_INITIAL_FLUX=1
    unset XK_PRIVATE_GS_SINGLE_PARTITION RESIDENT_OUTER_STEPS RESIDENT_STATE_DUMPS
    unset XK_PRIVATE_PRESSURE_CAPTURE RESIDENT_CUPTI_TRANSFER_LOG RESIDENT_CUPTI_KERNELS
    unset XK_PRIVATE_RESIDENT_PROFILE XK_SDK_PROFILE XK_PRIVATE_RESIDENT_SYSTEM_FINGERPRINT
    unset XK_PRIVATE_RESIDENT_POINTER_TRACE XK_PRIVATE_RESIDENT_STORAGE_LOG XK_PRIVATE_RESIDENT_MEMORY_LOG
}
