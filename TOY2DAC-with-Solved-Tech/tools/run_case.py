#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Generate synthetic observations and run one Marmousi 3D inversion case."""
import argparse
import array
from collections import Counter
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "upstream/TOY2DAC-2.6/run_marmousi_template/vp_Marmousi_exact"
GRIDS = json.loads((ROOT / "cases/fwi-grids.json").read_text())
NATIVE = ("analysis", "numeric", "solve")
CORE = {"initial_gradient", "iteration_gradient", "line_search"}


def fields(line):
    return dict(token.split("=", 1) for token in line.split()[1:] if "=" in token)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def check(condition, reason):
    if not condition:
        raise ValueError(reason)


def binary_values(path, count):
    values = array.array("d")
    check(path.stat().st_size == count * values.itemsize, f"unexpected size of {path.name}")
    values.frombytes(path.read_bytes())
    if sys.byteorder == "big":
        values.byteswap()
    check(all(map(math.isfinite, values)), f"nonfinite value in {path.name}")
    return values


def scan(path, grid):
    log = path.read_text()
    systems = re.split(r"(?=^FWI_SYSTEM )", log, flags=re.M)[1:]
    check(bool(systems), f"{path.name}: no linear systems")
    counts = Counter()
    timings = {}
    residuals = []
    for block in systems:
        system = fields(block.splitlines()[0])
        stage = system["stage"]
        check(int(system["nrhs"]) == grid["nrhs"], f"{stage}: wrong RHS count")
        configs = re.findall(r"^PANGULU_CONFIG .+$", block, re.M)
        check(len(configs) == 1, f"{stage}: missing solver configuration")
        config = fields(configs[0])
        check((int(config["n"]), int(config["nnz"])) == (grid["n"], grid["nnz"]),
              f"{stage}: unexpected matrix dimensions")
        check(config["version"] == "26.46" and config["precision"] == "complex128",
              f"{stage}: unexpected solver or precision")
        native = re.findall(r"^PANGULU_PHASE phase=(\w+) seconds=(\S+) status=(-?\d+)$", block, re.M)
        check(tuple(p for p, _, _ in native) == NATIVE, f"{stage}: incomplete native phases")
        phase_s = timings.setdefault(stage, dict.fromkeys(NATIVE, 0.0))
        for phase, duration, status in native:
            seconds = float(duration)
            check(status == "0" and math.isfinite(seconds) and seconds >= 0, f"{stage}: failed {phase}")
            phase_s[phase] += seconds
        matches = re.findall(r"^FWI_RESIDUAL .*columns=1 max_relative=(\S+)$", block, re.M)
        check(len(matches) == 1 and 0 <= float(matches[0]) <= 1e-8,
              f"{stage}: original COO residual check failed")
        residuals.append(float(matches[0]))
        counts[stage] += 1
    return log, timings, counts, max(residuals)


def physical_errors(grid, output):
    ex, ey, ez = grid["element_grid_xyz"]
    nx, ny, nz = ex + 1, ey + 1, ez + 1
    points = nx * ny * nz
    truth = binary_values(output / "observations/vp-true.bin", points)
    initial = binary_values(output / "observations/vp-initial.bin", points)
    recovered = binary_values(output / "vp-final.bin", points)
    interior = [((z * ny + y) * nx + x) for z in range(1, nz-1)
                for y in range(1, ny-1) for x in range(1, nx-1)]
    norm = math.fsum(truth[i] ** 2 for i in interior)
    check(norm > 0 and all(1200 - 1e-8 <= v <= 5000 + 1e-8 for v in recovered),
          "invalid final velocity model")
    vp_errors = [math.sqrt(math.fsum((v[i]-truth[i])**2 for i in interior) / norm)
                 for v in (initial, recovered)]

    fitted, heldout = [], []
    for y in (3, 5, 7, 9, 11):
        for x in range(5, 96, 2):
            rec = (y-3)//2*46 + (x-5)//2
            if (x-49)**2 + (y-7)**2 >= 25:
                group = heldout if ((x-5)//2+2*y) % 5 == 0 else fitted
                group.extend(range(6*rec, 6*rec+6))
    wave_squares = [[], []]
    for j in range(len(grid["frequencies_hz"])):
        values = 6 * grid["receivers"]
        obs = binary_values(output / "observations" / f"data-{j}.bin", values)
        pred = binary_values(output / f"predicted-{j}.bin", values)
        for group, dst in zip((fitted, heldout), wave_squares):
            energy = math.fsum(obs[i] ** 2 for i in group)
            check(energy > 0, "empty observation energy")
            dst.append(math.fsum((pred[i]-obs[i])**2 for i in group) / energy)
    return vp_errors, [math.sqrt(sum(v)/len(v)) for v in wave_squares]


def validate(grid, output, wall):
    observations, _, obs_counts, obs_residual = scan(output / "generate.stdout.log", grid)
    check("FWI_OBSERVATIONS_COMPLETE" in observations and obs_counts == {"synthetic_data": 2},
          "missing synthetic observations")
    log, phases, systems, residual = scan(output / "invert.stdout.log", grid)
    iterations = grid["iterations"]
    check(f"FWI_COMPLETE accepted_iterations={iterations} requested_iterations={iterations}" in log,
          "the inversion did not complete")
    matched = re.search(r"^FWI_MATCHED_FORWARD relative_receiver_l2=(\S+)$", log, re.M)
    check(matched is not None and 0 <= float(matched[1]) <= 1e-8, "receiver wavefield check failed")
    gradients = [float(v) for v in re.findall(r"^FWI_GRAD_CHECK .*relative_error=(\S+)$", log, re.M)]
    check(len(gradients) == 3 and all(math.isfinite(x) and x >= 0 for x in gradients)
          and gradients[-1] < 1e-3, "directional derivative check failed")
    rows = [{name: float(value) for name, value in fields(line).items()}
            for line in log.splitlines() if line.startswith("FWI_ITER ")]
    check([row["step"] for row in rows] == list(range(iterations+1)), "invalid iteration sequence")
    check(all(math.isfinite(v) for row in rows for v in row.values()) and
          all(b["objective"] <= a["objective"] for a, b in zip(rows, rows[1:])),
          "invalid objective history")
    vp_errors, waveform = physical_errors(grid, output)
    close = lambda a, b: math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-12)
    check(close(vp_errors[0], rows[0]["vp_relative"]) and close(vp_errors[1], rows[-1]["vp_relative"]) and
          close(waveform[0], rows[-1]["waveform_relative"]) and close(waveform[1], rows[-1]["heldout_relative"]),
          "recomputed output errors disagree with the application")
    core_phase_s = {phase: sum(times[phase] for name, times in phases.items() if name in CORE)
                    for phase in NATIVE}
    result = dict(accepted_iterations=iterations, max_original_coo_residual=max(obs_residual,residual),
                  matched_forward_receiver_l2=float(matched[1]), gradient_relative_errors=gradients,
                  initial_and_final_vp_l2=vp_errors, final_fit_and_heldout_l2=waveform, iterations=rows,
                  native_phase_s_by_stage=phases, core_phase_s=core_phase_s, core_s=sum(core_phase_s.values()),
                  core_systems=sum(n for name,n in systems.items() if name in CORE),
                  linear_systems_by_stage=systems, wall_s=wall)
    save(output / "validation.json", result)
    print(f"{iterations} updates, residual {result['max_original_coo_residual']:.3g}, "
          f"solver time {result['core_s']:.3f} s; outputs in {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=GRIDS)
    parser.add_argument("--output", type=Path, required=True, help="new run directory")
    parser.add_argument("--build", type=Path, default=ROOT / "_work/build")
    parser.add_argument("--prepare-only", action="store_true", help="copy inputs without running")
    args = parser.parse_args()
    grid, output, build = GRIDS[args.case], args.output.resolve(), args.build.resolve()
    binary = build / "bin" / f"marmousi3d_fwi_h{grid['spacing_m']}"
    if not args.prepare_only:
        if not binary.is_file() or not (build / "runtime.json").is_file():
            parser.error(f"build the application first (expected {binary})")
        runtime = json.loads((build / "runtime.json").read_text())["library_paths"]
    output.mkdir(parents=True)
    save(output / "case.json", dict(grid, case=args.case))
    model = output / MODEL.name
    shutil.copy2(MODEL, model)
    if args.prepare_only:
        print(output)
        return
    obs = output / "observations"
    obs.mkdir()
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = ":".join(filter(None, [*runtime, env.get("LD_LIBRARY_PATH")]))
    wall = {}
    for operation in ("generate", "invert"):
        start = time.perf_counter()
        with (output / f"{operation}.stdout.log").open("w") as out, (output / f"{operation}.stderr.log").open("w") as err:
            result = subprocess.run([str(binary), operation, str(model), str(obs)], cwd=output, env=env,
                                    stdout=out, stderr=err)
        wall[operation] = time.perf_counter() - start
        if result.returncode:
            raise SystemExit(f"{operation} exited with {result.returncode}; see {output / (operation + '.stderr.log')}")
    try:
        validate(grid, output, wall)
    except (ValueError, KeyError, OSError) as error:
        raise SystemExit(f"validation failed: {error}") from error


if __name__ == "__main__":
    main()
