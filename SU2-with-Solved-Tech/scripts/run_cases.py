#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Run the five SU2 cases with the CPU or XianzhiKrylov backend."""
import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parent.parent
CASES = ("Cavity", "Kenics", "IncTurbulentFlatPlate", "IncNACA", "Urban")
LINEAR_SOLVER_ITER = {
    "Cavity": 30, "Kenics": 20, "IncTurbulentFlatPlate": 20, "IncNACA": 10, "Urban": 10,
}

# Preserved from the integration archive's run_a100_compare.py. These select
# implementation paths in the matching SDK build, alongside its public API.
BLOCK_FASTPATH = {
    "SU2_CSR_RECBLOCK_GPU_ILU_MODE": "level_cta_map32_sync",
    "SU2_FGMRES_CGS_PASSES": "2",
    "SU2_FGMRES_BATCHED_CGS": "1",
    "XK_BSR_GPU_SYMBOLIC": "1",
    "XK_BSR_JACOBI_FIXED": "1",
    "XK_BSR_JACOBI_COOPERATIVE": "1",
    "XK_CSR_CGS_SINGLE_READBACK": "1",
    "XK_CSR_ILU_DIAGMAP": "1",
}
FASTPATH_ENV = {
    "IncTurbulentFlatPlate": {**BLOCK_FASTPATH, "SU2_GPU_REORDER": "1"},
    "IncNACA": {**BLOCK_FASTPATH, "SU2_GPU_REORDER": "1"},
    "Urban": {
        **BLOCK_FASTPATH, "XK_MSR_CACHED_BLOCK3": "1",
        "XK_MSR_COMPACT_WORKSPACE": "1", "XK_MSR_COMPACT_PAIR": "1",
        "XK_MSR_SUBWARP": "1", "XK_MSR_OWNED_INDEX_CACHE": "6",
        "XK_LDU_CHAIN_PLAN": "1", "XK_LDU_NUMERIC_CHAIN": "1",
        "XK_LDU_APPLY_TAGS": "1",
    },
}


def patch_config(source, destination, changes):
    lines = []
    remaining = dict(changes)
    for line in source.read_text().splitlines():
        match = re.match(r"^\s*([A-Z0-9_]+)\s*=", line)
        if match and match.group(1) in changes:
            key = match.group(1)
            lines.append(f"{key}= {changes[key]}")
            remaining.pop(key, None)
        else:
            lines.append(line)
    lines += [f"{key}= {value}" for key, value in remaining.items()]
    destination.write_text("\n".join(lines) + "\n")


def read_calls(path):
    records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    if not records:
        raise ValueError(f"No linear-call records in {path}")
    return records


def run_case(case, args, output):
    source = ROOT / "examples" / case
    output.mkdir()
    for mesh in source.glob("*.su2.gz"):
        with gzip.open(mesh, "rb") as compressed, (output / mesh.stem).open("wb") as expanded:
            shutil.copyfileobj(compressed, expanded)
    changes = {
        "ITER": args.outer_iterations, "CONV_STARTITER": args.outer_iterations + 1,
        "RESTART_SOL": "NO", "LINEAR_SOLVER": "FGMRES", "LINEAR_SOLVER_PREC": "ILU",
        "LINEAR_SOLVER_ILU_FILL_IN": 0, "ENABLE_CUDA": "YES" if args.backend == "gpu" else "NO",
        "LINEAR_SOLVER_ITER": LINEAR_SOLVER_ITER[case],
    }
    if args.inner_iterations is not None:
        changes["LINEAR_SOLVER_ITER"] = args.inner_iterations
    patch_config(source / "run_100.cfg", output / "run.cfg", changes)
    env = os.environ.copy()
    env.pop("SU2_XK_LEVEL", None)
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_PROC_BIND="close", OMP_PLACES="cores",
               SU2_XIANZHI_BACKEND=args.backend, SU2_XIANZHI_OUTPUT_DIR=str(output),
               SU2_XIANZHI_VERIFY="1" if args.verify else "0",
               SU2_XK_NUMERIC_REFRESH=str(args.numeric_refresh))
    if args.backend == "gpu":
        env["SU2_XK_ZERO"] = "1"
        env.update(FASTPATH_ENV.get(case, {}))
    recorded_env = {key: value for key, value in env.items()
                    if key.startswith(("SU2_", "XK_", "OMP_"))}
    # Keep environment records specific to solver controls. Licences and the
    # caller's other environment variables are not copied into result files.
    recorded_env = {key: value for key, value in recorded_env.items()
                    if not any(word in key.upper() for word in ("LICENSE", "TOKEN", "PASSWORD", "SECRET"))}
    (output / "solver-environment.json").write_text(json.dumps(recorded_env, indent=2) + "\n")
    started = time.monotonic()
    with (output / "application.log").open("w") as log:
        try:
            result = subprocess.run([str(args.binary), "run.cfg"], cwd=output, env=env,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout)
            returncode = result.returncode
        except subprocess.TimeoutExpired:
            returncode = 124
            log.write(f"\nRunner timeout after {args.timeout} seconds.\n")
    row = {
        "case": case, "backend": args.backend, "threads": args.threads,
        "outer_iterations": args.outer_iterations, "numeric_refresh": args.numeric_refresh,
        "verification": args.verify, "returncode": returncode,
        "application_wall_seconds": time.monotonic() - started,
        "whole_linear_seconds": None, "gpu_execute_seconds": None,
        "gpu_prepare_seconds": None, "gpu_solve_seconds": None,
        "linear_calls": 0, "linear_iterations": None, "mean_reported_residual": None,
    }
    try:
        calls = read_calls(output / "whole-linear-calls.jsonl")
        row.update(whole_linear_seconds=sum(call["wall_seconds"] for call in calls),
                   linear_calls=len(calls),
                   linear_iterations=sum(call["iterations"] for call in calls),
                   mean_reported_residual=sum(call["residual"] for call in calls) / len(calls))
        if args.backend == "gpu":
            gpu_calls = read_calls(output / "linear-solves.jsonl")
            row["gpu_execute_seconds"] = sum(call["execute_seconds"] for call in gpu_calls)
            row["gpu_prepare_seconds"] = sum(call["analysis_seconds"] + call["numeric_prepare_seconds"]
                                             for call in gpu_calls)
            row["gpu_solve_seconds"] = sum(call["solve_seconds"] for call in gpu_calls)
            if not math.isfinite(row["gpu_execute_seconds"]):
                raise ValueError("Non-finite GPU timing")
        if not all(math.isfinite(row[key]) for key in ("whole_linear_seconds", "mean_reported_residual")):
            raise ValueError("Non-finite linear-call statistics")
    except (OSError, KeyError, ValueError) as error:
        row["statistics_error"] = str(error)
        row["returncode"] = returncode or 1
        row.update(whole_linear_seconds=None, gpu_execute_seconds=None,
                   gpu_prepare_seconds=None, gpu_solve_seconds=None,
                   linear_calls=0, linear_iterations=None, mean_reported_residual=None)
    (output / "summary.json").write_text(json.dumps(row, indent=2, allow_nan=False) + "\n")
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("cpu", "gpu"), required=True)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--binary", type=Path, default=ROOT / "build/xianzhikrylov/install/bin/SU2_CFD")
    parser.add_argument("--output", type=Path, required=True, help="New output directory.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--outer-iterations", type=int, default=20)
    parser.add_argument("--inner-iterations", type=int, help="Override the selected cases' latest-validation budgets.")
    parser.add_argument("--numeric-refresh", type=int, default=1,
                        help="GPU preconditioner rebuild interval. 1 rebuilds every call, 0 reuses factors.")
    parser.add_argument("--verify", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--timeout", type=float, default=900)
    args = parser.parse_args()
    if min(args.threads, args.outer_iterations, args.timeout) <= 0 or args.numeric_refresh < 0:
        parser.error("Threads, iterations and timeout must be positive, refresh interval must be nonnegative.")
    if args.inner_iterations is not None and args.inner_iterations < 1:
        parser.error("--inner-iterations must be positive.")
    if len(args.cases) != len(set(args.cases)):
        parser.error("Select each case once.")
    args.binary = args.binary.resolve()
    if not args.binary.is_file() or not os.access(args.binary, os.X_OK):
        parser.error(f"Executable not found: {args.binary}")
    args.output = args.output.resolve()
    if args.output.exists():
        parser.error("The output directory already exists. Choose a new --output.")
    for case in args.cases:
        if not (ROOT / "examples" / case / "run_100.cfg").is_file():
            parser.error(f"Missing configuration for {case}.")
        if not list((ROOT / "examples" / case).glob("*.su2.gz")):
            parser.error(f"Missing compressed mesh for {case}.")
    args.output.mkdir(parents=True)
    rows = []
    for case in args.cases:
        row = run_case(case, args, args.output / case)
        rows.append(row)
        (args.output / "summary.json").write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
        print(json.dumps(row, ensure_ascii=False, allow_nan=False), flush=True)
    fields = [key for key in rows[0] if key != "statistics_error"]
    with (args.output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return 0 if all(row["returncode"] == 0 for row in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
