#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
: "${FRONTISTR_SRC:?Set FRONTISTR_SRC to a patched FrontISTR 5.9 checkout}"
cmake=${CMAKE:-cmake}
"$cmake" -S "$FRONTISTR_SRC" -B "$root/build/frontistr" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS=-ffree-line-length-none \
    -DWITH_MPI=ON -DWITH_OPENMP=ON -DWITH_MUMPS=ON -DWITH_MKL=OFF \
    -DWITH_METIS=OFF -DWITH_ML=OFF -DWITH_REFINER=OFF -DWITH_REVOCAP=OFF -DWITH_NETCDF=OFF \
    -DFRONTISTR_PANGULU_BRIDGE="$root/build/libfrontistr_pangulu.so" "$@"
"$cmake" --build "$root/build/frontistr" --target fistr1 --parallel "${BUILD_JOBS:-4}"
