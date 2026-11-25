#!/usr/bin/env python3
"""Run one published FrontISTR model with PanguLU in a new directory."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

from validate import validate

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("conrod", "hinge", "heat"))
    parser.add_argument("output", type=Path, help="new run directory")
    parser.add_argument("--fistr1", default="fistr1", help="patched fistr1 executable")
    args = parser.parse_args()
    exe = shutil.which(args.fistr1)
    if exe is None:
        parser.error(f"fistr1 executable not found: {args.fistr1}")
    exe = str(Path(exe).resolve())
    case = ROOT / "cases" / args.case
    meta = json.loads((case / "model.json").read_text())
    mesh = meta["basename"] + ".msh"
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    digest = hashlib.sha256()
    with gzip.open(case / (mesh + ".gz"), "rb") as src, (output / mesh).open("wb") as dst:
        for chunk in iter(lambda: src.read(1024 * 1024), b""):
            digest.update(chunk)
            dst.write(chunk)
    if digest.hexdigest() != meta["mesh_sha256"]:
        raise ValueError("Unpacked mesh differs from the reference")
    for name in (meta["basename"] + ".cnt", "hecmw_ctrl.dat"):
        shutil.copyfile(case / name, output / name)
    env = dict(os.environ, FRONTISTR_SOLVER="pangulu")
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("FRONTISTR_PANGULU_THREADS", "40")
    with (output / "application.log").open("w") as log:
        subprocess.run([exe], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(json.dumps(validate(meta, output), indent=2))


if __name__ == "__main__":
    main()
