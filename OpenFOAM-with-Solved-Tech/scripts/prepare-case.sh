#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

[[ $# -eq 2 ]] || die "Usage: $0 POLYMESH_DIRECTORY PREPARED_CASE_DIRECTORY"
require_openfoam
mesh=$1
destination=$2
[[ ! -e $destination && ! -L $destination ]] || die "Path already exists: $destination"
for field in points faces owner neighbour boundary; do
    [[ -s $mesh/$field || -s $mesh/$field.gz ]] ||
        die "Missing polyMesh/$field in $mesh."
done

mkdir -p -- "$(dirname -- "$destination")"
staged=$(mktemp -d "$(dirname -- "$destination")/.occ-preparation.XXXXXXXX")
printf 'Preparing OCC case in %s\n' "$staged"
cp -R -- "$case_template/0.orig" "$staged/0"
cp -R -- "$case_template/constant" "$case_template/system" "$staged/"
cp -a --reflink=auto -- "$mesh" "$staged/constant/polyMesh"

mpi_run=${MPI_RUN:-mpirun}
(
    cd -- "$staged"
    decomposePar -constant > log.decomposePar 2>&1
    "$mpi_run" -np 20 renumberMesh -constant -overwrite -parallel > log.renumberMesh 2>&1
    "$mpi_run" -np 20 potentialFoam -initialiseUBCs -parallel > log.potentialFoam 2>&1
    "$mpi_run" -np 20 applyBoundaryLayer -ybl 0.0450244 -parallel > log.applyBoundaryLayer 2>&1
)
mv -T -- "$staged" "$destination"
printf 'Prepared case: %s\n' "$destination"
