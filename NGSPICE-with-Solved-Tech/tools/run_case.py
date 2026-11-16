#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Capture and replay the original ngspice linear systems for one published case."""
import argparse
import filecmp
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]
CASES = json.loads((ROOT / "cases/manifest.json").read_text())["cases"]
CAPTURE = re.compile(
    r"^NGSPICE_CAPTURE call=(\d+) n=(\d+) nnz=(\d+) graph=(\d+) values=(\d+) residual=(\S+) time=(\S+)$",
    re.MULTILINE,
)
NATIVE = re.compile(r"^NATIVE call=(\d+) phase=(\w+) seconds=(\S+) status=(\d+)$", re.MULTILINE)
VERIFY = re.compile(r"^VERIFY call=(\d+) residual=(\S+) klu_relative_l2=(\S+)$", re.MULTILINE)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def run_logged(command, directory, name, environment):
    start = time.perf_counter()
    log = directory / name
    with log.open("x") as stream:
        result = subprocess.run(command, cwd=directory, env=environment,
                                stdout=stream, stderr=subprocess.STDOUT)
    seconds = time.perf_counter() - start
    require(result.returncode == 0, f"process exited {result.returncode}; see {log}")
    return log.read_text(errors="replace"), seconds


def check_raw(path):
    """Verify that the real binary operating-point plot contains one full record."""
    with path.open("rb") as stream:
        meta = {}
        for line in stream:
            if line == b"Binary:\n":
                entries = int(meta.get(b"no. variables", b"0"))
                require(meta.get(b"flags") == b"real" and
                        meta.get(b"no. points") == b"1" and
                        meta.get(b"plotname") == b"Operating Point" and
                        entries >= 2 and path.stat().st_size - stream.tell() == 8 * entries,
                        "incomplete operating-point output")
                return
            if not line.startswith(b"\t") and b":" in line:
                key, value = line.split(b":", 1)
                meta[key.lower()] = value.strip()
    raise RuntimeError("missing binary operating-point output")


def read_capture(directory, case):
    log = (directory / "ngspice.log").read_text(errors="replace")
    require(not re.search(r"analysis .*aborted|timestep too small|^Error:|^CAPTURE_ERROR",
                          log, re.IGNORECASE | re.MULTILINE), "ngspice reported an error")
    require("Using KLU as Direct Linear Solver" in log, "expected native KLU")
    check_raw(directory / "field.raw")
    stems = []
    records = []
    for i, row in enumerate(CAPTURE.findall(log), start=1):
        call, n, nz, graph, values = map(int, row[:5])
        residual, simulation_time = map(float, row[5:])
        require((call, n, nz, graph, values) ==
                (i, case["matrix_order"], case["csr_entries"], 1, 1),
                "unexpected system dimensions or factor reuse sequence")
        require(simulation_time == 0 and math.isfinite(residual) and 0 <= residual <= 1e-8,
                "invalid KLU operating-point reference")
        stem = directory / "systems" / f"call-{i:04d}"
        for suffix, size in ((".b.f64", n * 8), (".x.f64", n * 8)):
            require(Path(str(stem) + suffix).stat().st_size == size, f"invalid {stem}{suffix}")
        if i == 1:
            for suffix, size in ((".rp.i64", (n + 1) * 8), (".ci.i32", nz * 4), (".a.f64", nz * 8)):
                require(Path(str(stem) + suffix).stat().st_size == size, f"invalid {stem}{suffix}")
        else:
            require(filecmp.cmp(str(stems[0]) + ".b.f64", str(stem) + ".b.f64", shallow=False),
                    "RHS differs from the first native solve")
        stems.append(stem)
        records.append(dict(call=call, matrix_order=n, csr_entries=nz,
                            graph_from=graph, coefficients_from=values, klu_residual=residual))
    require(len(stems) == case["calls"], "unexpected KLU solve count")
    return stems, records


def save_result(path, payload):
    with path.open("x") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")


def capture(args, case):
    with args.netlist.open("rb") as deck:
        require(hashlib.file_digest(deck, "sha256").hexdigest() == case["netlist_sha256"],
                "netlist does not match the published case")
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True)
    systems = directory / "systems"
    systems.mkdir()
    env = dict(os.environ, NGSPICE_CAPTURE=str(systems), OMP_NUM_THREADS="1",
               OMP_DYNAMIC="FALSE", SPICE_ASCIIRAWFILE="0")
    _, wall = run_logged([args.ngspice.resolve(), "-n", "-b", "-r", "field.raw",
                          args.netlist.resolve()], directory, "ngspice.log", env)
    _, records = read_capture(directory, case)
    save_result(directory / "capture.json", dict(case=args.case, wall_s=wall, calls=records))


def replay(args, case):
    capture_dir = args.capture_dir.resolve()
    marker = json.loads((capture_dir / "capture.json").read_text())
    require(marker["case"] == args.case, "capture directory belongs to another case")
    stems, records = read_capture(capture_dir, case)
    require(marker["calls"] == records, "capture metadata does not describe the native log")
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.host_threads), OMP_DYNAMIC="FALSE")
    order, entries = case["matrix_order"], case["csr_entries"]
    cmd = [args.pangulu_replay.resolve(), str(order), str(entries), *stems]
    log, wall = run_logged(cmd, directory, "pangulu.log", env)
    phases = [
        dict(call=int(call), phase=phase, seconds=float(seconds), status=int(status))
        for call, phase, seconds, status in NATIVE.findall(log)
    ]
    expected = [(1, "analysis"), (1, "numeric")]
    expected += [(call, "solve") for call in range(1, case["calls"] + 1)]
    require([(p["call"], p["phase"]) for p in phases] == expected and
            all(p["status"] == 0 and math.isfinite(p["seconds"]) and p["seconds"] >= 0 for p in phases),
            "unexpected native phase sequence")
    checks = [
        dict(call=int(call), residual=float(residual), klu_relative_l2=float(error))
        for call, residual, error in VERIFY.findall(log)
    ]
    require([v["call"] for v in checks] == list(range(1, case["calls"] + 1)) and
            all(0 <= v["residual"] <= 1e-8 and 0 <= v["klu_relative_l2"] <= 1e-6 for v in checks) and
            f"REPLAY_COMPLETE calls={case['calls']}" in log,
            "incomplete or failed solution validation")
    result = dict(case=args.case, matrix_order=order, csr_entries=entries, host_threads=args.host_threads,
                  native_solver_s=sum(p["seconds"] for p in phases), replay_wall_s=wall,
                  phases=phases, validation=checks)
    gpu = re.search(r"^GPU name=(.+) runtime=(\d+)$", log, re.MULTILINE)
    require(gpu is not None, "missing GPU runtime identity")
    result.update(gpu=gpu[1], cuda_runtime=int(gpu[2]))
    save_result(directory / "result.json", result)
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    native = modes.add_parser("capture")
    native.add_argument("--ngspice", type=Path, required=True)
    native.add_argument("--netlist", type=Path, required=True)
    gpu = modes.add_parser("replay")
    gpu.add_argument("--pangulu-replay", type=Path, required=True)
    gpu.add_argument("--capture-dir", type=Path, required=True)
    gpu.add_argument("--host-threads", type=int, default=1, help="PanguLU host OpenMP threads")
    for subcommand in (native, gpu):
        subcommand.add_argument("case", choices=CASES)
        subcommand.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "capture":
        capture(args, CASES[args.case])
    else:
        require(args.host_threads > 0, "--host-threads must be positive")
        replay(args, CASES[args.case])


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        raise SystemExit(str(error)) from error
