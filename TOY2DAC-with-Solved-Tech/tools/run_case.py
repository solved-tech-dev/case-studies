#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Prepare or run one Marmousi 3D case."""
import argparse
import array
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "upstream/TOY2DAC-2.6/run_marmousi_template/vp_Marmousi_exact"
GRIDS = json.loads((ROOT / "cases/grids.json").read_text())
PHASES = ("analysis", "numeric", "solve")


def validate(grid, output):
    stdout = (output / "stdout.log").read_text()
    match = re.search(r"^MARMOUSI3D_RESIDUAL n=(\d+) nnz=(\d+) relative_l2=(\S+)$", stdout, re.M)
    if not match or (int(match[1]), int(match[2])) != (grid["n"], grid["nnz"]):
        raise ValueError("missing or inconsistent matrix and residual record")
    residual = float(match[3])
    if not 0 <= residual <= 1e-8:
        raise ValueError(f"relative residual check failed ({residual})")
    records = re.findall(r"^PANGULU_PHASE phase=(\w+) seconds=(\S+) status=(-?\d+)$", stdout, re.M)
    if tuple(row[0] for row in records) != PHASES or any(row[2] != "0" for row in records):
        raise ValueError("missing or unsuccessful native solver phases")
    seconds = {phase: float(elapsed) for phase, elapsed, _ in records}
    if any(not math.isfinite(t) or t < 0 for t in seconds.values()):
        raise ValueError("invalid native phase timings")
    outputs = ("field_complex.bin", "field.vtk", "model.vtk", "midplane_y.csv", "metadata.txt")
    for name in outputs:
        path = output / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing or empty output {name}")
    expected_shape = tuple(n + 1 for n in grid["element_grid_xyz"])
    field = output / "field_complex.bin"
    if field.stat().st_size != 12 + 16 * grid["n"]:
        raise ValueError("incomplete three-component complex field")
    with field.open("rb") as stream:
        if struct.unpack("<3i", stream.read(12)) != expected_shape:
            raise ValueError("field grid does not match the case")
        values = array.array("d")
        values.fromfile(stream, 2 * grid["n"])
        if sys.byteorder == "big":
            values.byteswap()
        if not all(map(math.isfinite, values)):
            raise ValueError("non-finite displacement field")
    result = {"n": grid["n"], "nnz": grid["nnz"], "nrhs": grid["nrhs"], "relative_l2": residual,
              "node_grid_xyz": expected_shape, "field_finite": True, "phase_s": seconds,
              "core_s": sum(seconds.values())}
    (output / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{grid['n']:,} equations, relative residual {residual:.3g}; outputs in {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=GRIDS)
    parser.add_argument("--output", type=Path, required=True, help="new run directory")
    parser.add_argument("--build", type=Path, default=ROOT / "_work/build")
    parser.add_argument("--prepare-only", action="store_true", help="copy the inputs without building or running")
    args = parser.parse_args()
    grid, output, build = GRIDS[args.case], args.output.resolve(), args.build.resolve()
    binary = build / "bin/marmousi3d_pangulu"
    if not args.prepare_only:
        if not binary.is_file() or not (build / "runtime.json").is_file():
            parser.error(f"build the application first (expected {binary})")
        runtime = json.loads((build / "runtime.json").read_text())["library_paths"]
    output.mkdir(parents=True)
    (output / "case.json").write_text(json.dumps(dict(grid, case=args.case), indent=2) + "\n")
    model = output / MODEL.name
    shutil.copy2(MODEL, model)
    if args.prepare_only:
        print(output)
        return
    nx, ny, nz = grid["element_grid_xyz"]
    command = [str(binary), "--nx", str(nx), "--ny", str(ny), "--nz", str(nz),
               "--spacing", str(grid["spacing_m"]), "--frequency", str(grid["frequency_hz"]),
               "--marmousi", str(model), "--output", str(output)]
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = ":".join(filter(None, [*runtime, env.get("LD_LIBRARY_PATH")]))
    with (output / "stdout.log").open("w") as out, (output / "stderr.log").open("w") as err:
        result = subprocess.run(command, cwd=output, env=env, stdout=out, stderr=err)
    if result.returncode:
        raise SystemExit(f"application exited with {result.returncode}; see {output / 'stdout.log'} and stderr.log")
    try:
        validate(grid, output)
    except (ValueError, EOFError, OSError) as error:
        raise SystemExit(f"validation failed: {error}") from error


if __name__ == "__main__":
    main()
