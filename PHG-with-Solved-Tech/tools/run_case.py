#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Run a PHG case and check the solver, assembled size, and native outputs."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=sorted(p.parent.name for p in (ROOT / "cases").glob("*/case.json")))
    parser.add_argument("--workspace", default="_work/clean")
    parser.add_argument("--profile", choices=("benchmark", "connection-test"), default="benchmark")
    parser.add_argument("--backend", choices=("pangulu", "test-reference"), default="pangulu")
    parser.add_argument("--output", required=True)
    parser.add_argument("--solver-threads", type=int, default=40)
    parser.add_argument("--vtk", action="store_true")
    parser.add_argument("--capture-linear", action="store_true")
    parser.add_argument("--timeout", type=int, help="Optional application limit in seconds")
    args = parser.parse_args()
    if args.backend != "pangulu" and args.profile != "connection-test":
        parser.error("Reference backends are available for connection-test only")
    if args.solver_threads <= 0 or args.timeout is not None and args.timeout <= 0:
        parser.error("Thread count and timeout must be positive")
    if args.capture_linear and args.backend != "pangulu":
        parser.error("--capture-linear requires the PanguLU backend")
    return args


def main():
    args = arguments()
    config = json.loads((ROOT / "cases" / args.case / "case.json").read_text())
    profile = config["profiles"][args.profile]
    workspace, output = Path(args.workspace).resolve(), Path(args.output).resolve()
    subprocess.run([sys.executable, str(ROOT / "tools/verify_sources.py"), "--workspace", str(workspace)], check=True)
    executable = workspace / "bin" / config["program"]
    env = dict(os.environ)
    for key in ("LD_PRELOAD", "PHG_REPRO_VTK", "PHG_REPRO_LINEAR_DIR", "PHG_PANGULU_BRIDGE"):
        env.pop(key, None)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PHG_REPRO_FIELDS=str(output / "fields"))
    if args.vtk:
        env["PHG_REPRO_VTK"] = "1"
    if args.backend == "pangulu":
        env.update(PHG_PANGULU_BRIDGE=str(workspace / "bin/libphg_pangulu.so"),
                   PHG_PANGULU_THREADS=str(args.solver_threads), PHG_PANGULU_REUSE_LEVEL=config["reuse"])
    elif args.backend == "test-reference":
        env["PHG_PANGULU_BRIDGE"] = str(workspace / "bin/libtest_reference.so")
    if "PHG_PANGULU_BRIDGE" in env and not Path(env["PHG_PANGULU_BRIDGE"]).is_file():
        raise SystemExit("Build the selected backend before running a case")
    output.mkdir(parents=True, exist_ok=False)
    (output / "fields").mkdir()
    if args.capture_linear:
        (output / "linear").mkdir()
        env["PHG_REPRO_LINEAR_DIR"] = str(output / "linear")

    command = [str(executable), "-solver", "pangulu", "-solver_symmetry", "unsym", "-vtk_precision", "double"]
    if config["program"] == "solid":
        tree = ET.parse(workspace / "solid" / config["input"])
        root = tree.getroot()
        for key, value in profile["xml_updates"].items():
            tag, attribute = key.split(".")
            node = root if tag == "model" else root.find(tag)
            if node is None or attribute not in node.attrib:
                raise ValueError(f"Unknown PHG-Solid XML attribute {key}")
            node.set(attribute, value)
        tree.write(output / "input.xml", encoding="utf-8", xml_declaration=True)
        shutil.copyfile(workspace / "solid" / config["mesh"], output / "cantilever.mesh")
        command.extend(["-f", "input.xml"])
    else:
        command.extend(["-mesh_file", str(workspace / "phg-0.9.7" / config["mesh"]),
                        *config["common_args"], *profile["args"]])

    start = time.monotonic()
    with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
        process = subprocess.Popen(command, cwd=output, env=env, stdout=stdout, stderr=stderr,
                                   start_new_session=True)
        timed_out = False
        try:
            rc = process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                rc = process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                rc = process.wait()
    elapsed = time.monotonic() - start
    stderr = (output / "stderr.log").read_text(errors="replace")
    stdout = (output / "stdout.log").read_text(errors="replace")
    residuals = [float(x) for x in re.findall(r"relative_residual=([\d.eE+-]+)", stderr)]
    systems = [list(map(int, row)) for row in re.findall(r"PHG_PANGULU_RESULT call=\d+ n=(\d+) nnz=(\d+)", stderr)]
    phases = [dict(call=int(call), phase=name, seconds=float(seconds), status=int(status), cuda_status=int(cuda))
              for call, name, seconds, status, cuda in re.findall(
                  r"PHG_PANGULU_PHASE call=(\d+) phase=(\w+) seconds=([\d.eE+-]+) status=(\d+) cuda_status=(\d+)", stderr)]
    fields = {}
    for meta in sorted((output / "fields").glob("*.meta"), key=lambda p: int(p.stem.split("-")[-1])):
        item = dict(line.split("=", 1) for line in meta.read_text().splitlines())
        data = meta.with_suffix(".f64")
        count = int(item["count"])
        if item["name"] in config["fields"]:
            fields[item["name"]] = dict(file=str(data.relative_to(output)), count=count,
                                        mesh_leaves=int(item["mesh_leaves"]),
                                        finite=item["finite"] == "1" and item["scalar_bytes"] == "8"
                                        and data.stat().st_size == 8 * count)
    vtk = list(output.glob("interface-*.vtk")) if config["program"] == "interface" else list((output / "fields").glob("*.vtk"))
    issues = []
    if rc or timed_out:
        issues.append("application failed or timed out")
    if set(fields) != set(config["fields"]) or not all(item["finite"] for item in fields.values()):
        issues.append("native fields missing or invalid")
    if (args.vtk or config["program"] == "interface") and not any(p.stat().st_size for p in vtk):
        issues.append("native VTK missing")
    if not residuals or any(not 0 <= v <= 1e-8 for v in residuals):
        issues.append("original-system residual missing or exceeds 1e-8")
    if args.backend == "pangulu":
        expected_phases = [(i, name) for i in range(1, len(systems) + 1)
                           for name in ("analysis", "numeric_factorisation", "solve")]
        if ("PHG_REPRO_BACKEND PanguLU public API 26.44" not in stdout or not systems
                or len(residuals) != len(systems)
                or [(p["call"], p["phase"]) for p in phases] != expected_phases
                or any(p["status"] or p["cuda_status"] or not math.isfinite(p["seconds"]) for p in phases)):
            issues.append("PanguLU identity or native phase sequence invalid")
        if "systems" in profile and systems != profile["systems"]:
            issues.append("original CSR sizes differ from the case profile")
        if systems and sum(f["count"] for f in fields.values()) != systems[-1][0]:
            issues.append("native field size does not match the final matrix order")
    elif args.backend == "test-reference" and "PHG_REPRO_BACKEND TEST_ONLY_CPU_REFERENCE" not in stdout:
        issues.append("CPU test bridge identity missing")

    result = dict(case=args.case, profile=args.profile, backend=args.backend, passed=not issues, issues=issues,
                  exit_code=rc, timed_out=timed_out, systems=systems, residuals=residuals,
                  native_phases=phases, native_core_s=sum(p["seconds"] for p in phases) if phases and not issues else None,
                  application_wall_s=elapsed, fields=fields, vtk=[str(p.relative_to(output)) for p in vtk])
    (output / "run-report.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 1 if issues else 0


if __name__ == "__main__":
    sys.exit(main())
