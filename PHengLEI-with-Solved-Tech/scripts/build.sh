#!/usr/bin/env bash
set -euo pipefail
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
mode=${1:-linglong}
if (($#)); then shift; fi
source_dir=$(cd -- "${PHENGLEI_SOURCE_DIR:-"$package/work/PHengLEI"}" && pwd)
if [[ ! -f "$source_dir/PHengLEI/3rdparty/eigen/Eigen/Eigen" ]]; then
  printf 'Eigen headers are missing. Prepare a new source copy with EIGEN_ROOT set.\n' >&2
  exit 2
fi
build_dir=${PHENGLEI_BUILD_DIR:-"$package/build/$mode"}
cmake=${CMAKE_BIN:-cmake}
hypre=$(cd -- "${HYPRE_ROOT:?Set HYPRE_ROOT to the HYPRE CPU installation}" && pwd)
compiler=${MPICXX:-mpicxx}
linglong=OFF
tests=OFF
extra=()
case "$mode" in
  linglong)
    linglong=ON
    tests=${PHENGLEI_LINGLONGAMG_TESTS:-OFF}
    sdk=$(cd -- "${LINGLONGAMG_ROOT:?Set LINGLONGAMG_ROOT to the SDK}" && pwd)
    cuda=$(cd -- "${CUDA_ROOT:?Set CUDA_ROOT to CUDA 13.x}" && pwd)
    extra=(-DLINGLONGAMG_ROOT="$sdk" -DCUDAToolkit_ROOT="$cuda")
    ;;
  hypre) ;;
  *) printf 'Usage: %s [linglong|hypre] [additional CMake options]\n' "$0" >&2; exit 2 ;;
esac
mkdir -p -- "$build_dir"
build_dir=$(cd -- "$build_dir" && pwd)
"$cmake" -S "$source_dir/PHengLEI" -B "$build_dir" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$compiler" \
  -DBUILD_API=ON -DBUILD_CFD=ON -DBUILD_MESH=ON -DBUILD_IncomSolver=ON \
  -DBUILD_GPU=OFF -DBUILD_GMRESSolver=OFF -DBUILD_SpecDiffHybSolver=OFF \
  -DBUILD_DemoSolver=OFF -DBUILD_LBMSolverMPI=OFF -DBUILD_LBMSolverOMP=OFF \
  -DUSE_HypreLib=ON -DUSE_UnapLib=OFF -DUSE_YhamgLib=OFF -DUSE_DOUBLE_PRECISION=ON \
  -DUSE_LinglongAMG="$linglong" -DPHENGLEI_LINGLONGAMG_TESTS="$tests" \
  -DPHENGLEI_HYPRE_ROOT="$hypre" \
  -DBUILD_LIBRARY_TYPE=Static -DPHENGLEI_LIBRARY_OUTPUT_ROOT="$build_dir/lib" \
  -DPHENGLEI_EXECUTABLE_NAME=PHengLEIv3d0 "${extra[@]}" "$@"
"$cmake" --build "$build_dir" --parallel "${BUILD_JOBS:-4}"
printf 'Executable %s\n' "$build_dir/PHengLEIv3d0/PHengLEIv3d0"
