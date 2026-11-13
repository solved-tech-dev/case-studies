#!/usr/bin/env python3
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
from shutil import copy2, copytree, rmtree
import subprocess
import time

PACKAGE = Path(__file__).resolve().parents[1]

def set_parameter(content, kind, name, value):
    pattern = rb"(?m)^([ \t]*" + kind + rb"[ \t]+" + name + rb"[ \t]*=[ \t]*)[^;\r\n]+;"
    content, count = re.subn(pattern, lambda match: match[1] + value + b";", content)
    if count > 1:
        raise ValueError(f"Duplicate parameter {name.decode()}")
    if count == 0:
        content += b"\n" + kind + b" " + name + b" = " + value + b";\n"
    return content

def stage(run, inputs, manifest, args):
    run.mkdir(parents=True, exist_ok=False)
    try:
        copytree(inputs / "bin", run / "bin")
        (run / "grid").mkdir()
        (run / "results").mkdir()
        for entry in manifest["grid"]:
            name = entry["file"]
            dest = run / "grid" / name.removesuffix(".gz")
            digest = hashlib.sha256()
            with gzip.open(inputs / "grid" / name, "rb") as compressed, dest.open("wb") as output:
                while data := compressed.read(1024 * 1024):
                    digest.update(data)
                    output.write(data)
            if digest.hexdigest() != entry["uncompressed_sha256"] or dest.stat().st_size != entry["uncompressed_bytes"]:
                raise ValueError(f"Grid checksum mismatch: {name}")

        path = run / "bin" / "cfd_para_incompressible.hypara"
        data = path.read_bytes()
        for key in (b"maxSimuStep", b"intervalStepFlow", b"intervalStepPlot"):
            data = set_parameter(data, b"int", key, str(args.steps).encode())
        data = set_parameter(data, b"int", b"mathLibType", b"0")
        data = set_parameter(data, b"int", b"useLinglongPressure", b"1" if args.backend == "linglong" else b"0")
        if args.backend == "linglong":
            copy2(PACKAGE / "examples" / "linglongamg_config.json", run / "bin")
            data = set_parameter(data, b"int", b"linglongDevice", str(args.device).encode())
            data = set_parameter(data, b"string", b"linglongConfigFile", b'"./bin/linglongamg_config.json"')
            data = set_parameter(data, b"string", b"linglongLogFile", b'"./results/linglong_pressure.csv"')
        path.write_bytes(data)
    except BaseException:
        rmtree(run)
        raise

def main():
    parser = argparse.ArgumentParser(description="Stage and run one single-zone PHengLEI case.")
    parser.add_argument("case", choices=("H11", "H12", "H13"))
    parser.add_argument("backend", choices=("linglong", "hypre"))
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.steps < 1 or args.device < 0:
        parser.error("--steps must be positive and --device must be nonnegative")
    inputs = PACKAGE / "examples" / args.case
    manifest = json.loads((inputs / "case.json").read_text())
    run = args.run_dir.resolve()
    exe = Path(os.environ.get("PHENGLEI_EXE", PACKAGE / "build" / args.backend / "PHengLEIv3d0" / "PHengLEIv3d0")).resolve()
    if not args.prepare_only and not (exe.is_file() and os.access(exe, os.X_OK)):
        parser.error(f"Build first or set PHENGLEI_EXE to an executable: {exe}")
    stage(run, inputs, manifest, args)
    summary = {"case": args.case, "backend": args.backend, "requested_steps": args.steps, "cells": manifest["cells"]}
    if args.prepare_only:
        summary["state"] = "prepared"
        result = 0
    else:
        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "1")
        started = time.monotonic()
        with (run / "solver.log").open("wb") as log:
            result = subprocess.run([str(exe)], cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        summary.update(
            state="exited", exit_code=result, wall_s=time.monotonic() - started,
            omp_threads=env["OMP_NUM_THREADS"],
        )
    (run / "run.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(run)
    return result if result >= 0 else 128 - result

if __name__ == "__main__":
    raise SystemExit(main())
