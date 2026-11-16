#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

[[ $# -ge 3 && $# -le 4 ]] ||
    die "Usage: $0 BACKEND PREPARED_CASE_DIRECTORY RUN_DIRECTORY [SIMPLE_STEPS]"
backend=$1
require_backend "$backend"
require_openfoam
sdk_environment "$backend"
prepared=$(cd -- "$2" && pwd -P)
destination=$3
steps=${4:-5}
[[ $steps =~ ^[1-9][0-9]*$ ]] || die "SIMPLE_STEPS must be a positive integer."
[[ ! -e $destination && ! -L $destination ]] || die "Path already exists: $destination"
[[ -d $prepared/constant ]] || die "Missing constant directory in $prepared."
for ((rank=0; rank<20; rank++)); do
    for directory in 0 constant/polyMesh; do
        [[ -d $prepared/processor$rank/$directory ]] ||
            die "Missing processor$rank/$directory in $prepared."
    done
done

source_root=${SOLVED_OPENFOAM_SOURCE:-"$package_dir/work/source"}
build_root=$(cd -- "$source_root/$backend" && pwd -P)
[[ -x $build_root/bin/$executable && -r $build_root/lib/$adapter ]] ||
    die "Run scripts/build.sh $backend before starting the case."
[[ -r $sdk_license ]] || die "Install the SDK license and set its path before running."
sdk_license=$(realpath -- "$sdk_license")

mkdir -p -- "$(dirname -- "$destination")"
mkdir -- "$destination"
run_dir=$(cd -- "$destination" && pwd -P)
case_dir="$run_dir/case"
runtime="$run_dir/runtime"
mkdir -p -- "$case_dir" "$runtime/bin" "$run_dir/results"
chmod 700 "$runtime"
cp -R -- "$case_template/system" "$case_dir/"
ln -s -- "$prepared/constant" "$case_dir/constant"
for ((rank=0; rank<20; rank++)); do
    processor="processor$rank"
    mkdir -- "$case_dir/$processor"
    ln -s -- "$prepared/$processor/0" "$prepared/$processor/constant" "$case_dir/$processor/"
done
ln -s -- "$sdk_license" "$runtime/license.key"
cp -p -- "$build_root/bin/$executable" "$runtime/bin/"
cp -- "$package_dir/examples/$backend.json" "$run_dir/solver.json"

ordering_args=()
if [[ $backend == linglongamg && -n ${LINGLONGAMG_ORDERING_DIR:-} ]]; then
    ordering_args=(--ordering "$LINGLONGAMG_ORDERING_DIR")
fi
python3 "$package_dir/scripts/configure-case.py" \
    "$backend" "$case_dir" "$run_dir/solver.json" "$steps" "${ordering_args[@]}"

export LD_LIBRARY_PATH="$build_root/lib:$LD_LIBRARY_PATH"
command_prefix=(env)
case "$backend" in
    linglongamg)
        export FOAM_PRESSURE_RESULTS="$run_dir/results"
        export FOAM_PRESSURE_FROZEN_REPEATS=1 FOAM_PRESSURE_HIERARCHY_REUSE=0
        export CUDA_CACHE_DISABLE=1
        command_prefix+=("LD_PRELOAD=$build_root/lib/libBenchIoGuard.so")
        ;;
    xianzhikrylov)
        export XK_BENCH_RESULTS="$run_dir/results"
        unset XK_CPU_RESIDUAL_TARGETS
        ;;
esac
unset FOAM_BENCH_CHECKPOINT_STEP FOAM_BENCH_STOP_AFTER_CHECKPOINT FOAM_BENCH_DIAGNOSTIC_FIELD_DUMP

printf 'Starting %s. Log: %s/log.openfoam\n' "$backend" "$run_dir"
(
    cd -- "$runtime"
    "${MPI_RUN:-mpirun}" -np 20 "${command_prefix[@]}" \
        "$runtime/bin/$executable" -case "$case_dir" -parallel -noFunctionObjects
) > "$run_dir/log.openfoam" 2>&1
python3 - "$run_dir/results/main.json" "$steps" <<'PY'
import json, sys
from pathlib import Path

completed = json.loads(Path(sys.argv[1]).read_text())["completed_steps"]
if completed != int(sys.argv[2]):
    raise SystemExit(f"OpenFOAM completed {completed} steps, expected {sys.argv[2]}.")
PY
printf 'Completed %s. Run directory: %s\n' "$backend" "$run_dir"
