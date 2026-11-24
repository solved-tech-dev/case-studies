#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/resident-common.sh"
require_openfoam
sdk_environment
command -v wmake >/dev/null || die "wmake is unavailable."
command -v mpicc >/dev/null || die "Use the MPI environment matching OpenFOAM."
export RESIDENT_MPI_CFLAGS="${RESIDENT_MPI_CFLAGS:-$(mpicc --showme:compile)}"
export RESIDENT_MPI_LIBS="${RESIDENT_MPI_LIBS:-$(mpicc --showme:link)}"
mkdir -p -- "$RESIDENT_BUILD/source" "$RESIDENT_BUILD/bin"
cp -R -- "$package_dir/overlay/resident/application/." "$RESIDENT_BUILD/source/"
(
    cd -- "$RESIDENT_BUILD/source"
    wmake -j "${BUILD_JOBS:-2}"
)
printf 'Built %s\n' "$RESIDENT_BUILD/bin/rawResident"
