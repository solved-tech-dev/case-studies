#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

[[ $# -ge 1 && $# -le 2 ]] || die "Usage: $0 BACKEND [source-work-directory]"
backend=$1
require_backend "$backend"
require_openfoam
sdk_environment "$backend"

source_root=${2:-"$package_dir/work/source"}
[[ -d $source_root/$backend ]] || die "Run scripts/prepare-source.sh first."
build_root=$(cd -- "$source_root/$backend" && pwd -P)
mkdir -p -- "$build_root/bin" "$build_root/lib"
mpi_cxx=${MPI_CXX:-mpicxx}
SOLVED_MPI_CFLAGS=$("$mpi_cxx" --showme:compile)
SOLVED_MPI_LIBS=$("$mpi_cxx" --showme:link)
export SOLVED_MPI_CFLAGS SOLVED_MPI_LIBS
export LD_LIBRARY_PATH="$build_root/lib:$LD_LIBRARY_PATH"
jobs=${BUILD_JOBS:-1}

case "$backend" in
    linglongamg)
        [[ -r $LINGLONGAMG_INCLUDE_DIR/linglongamg.h && -r $sdk_lib/liblinglongamg.so ]] ||
            die "Cannot find linglongamg.h or liblinglongamg.so in the configured SDK."
        export FOAM_PRESSURE_ROOT="$build_root"
        "${CXX:-g++}" -std=c++17 -O2 -fPIC -shared \
            "$build_root/adapter/BenchIoGuard.C" -o "$build_root/lib/libBenchIoGuard.so" -ldl
        ;;
    xianzhikrylov)
        [[ -r $XIANZHIKRYLOV_SDK/include/XianzhiKrylovAPI.h && -r $sdk_lib/libXianzhiKrylov.so ]] ||
            die "Cannot find XianzhiKrylovAPI.h or libXianzhiKrylov.so in the configured SDK."
        exported_symbols=$(nm -D --defined-only "$sdk_lib/libXianzhiKrylov.so")
        for symbol in xk_bench_set_tolerances xk_internal_get_last_stats; do
            grep -Eq "[[:space:]]${symbol}(@@[^[:space:]]+)?$" <<< "$exported_symbols" ||
                die "The selected XianzhiKrylov library does not export $symbol."
        done
        export BENCH_GPU_ROOT="$build_root"
        ;;
esac

(cd -- "$build_root/adapter" && wmake -j "$jobs" libso)
(cd -- "$build_root/application" && wmake -j "$jobs")
printf 'Executable: %s/bin/%s\n' "$build_root" "$executable"
