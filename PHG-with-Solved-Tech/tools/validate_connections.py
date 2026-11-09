#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Compare all three small native application fields."""
import argparse
import array
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--workspace", default="_work/clean")
    p.add_argument("--output", required=True)
    p.add_argument("--candidate", choices=("test-reference", "pangulu"), default="test-reference")
    args = p.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    reference = "test-reference"
    summary = []
    for case_file in sorted((ROOT / "cases").glob("*/case.json")):
        case = case_file.parent.name
        reports = {}
        dirs = {}
        for backend in dict.fromkeys((reference, args.candidate)):
            dirs[backend] = output / f"{case}-{backend}"
            command = [sys.executable, str(ROOT / "tools/run_case.py"), case, "--workspace", args.workspace,
                       "--profile", "connection-test", "--backend", backend,
                       "--output", str(dirs[backend]), "--timeout", "120", "--vtk"]
            if backend == "pangulu":
                command.extend(["--solver-threads", "4"])
            result = subprocess.run(command, capture_output=True, text=True)
            (output / f"{case}-{backend}.log").write_text(result.stdout + result.stderr)
            result.check_returncode()
            reports[backend] = json.loads((dirs[backend] / "run-report.json").read_text())
        errors = {}
        for name, expected in (reports[reference]["fields"].items() if args.candidate=="pangulu" else ()):
            actual = reports[args.candidate]["fields"][name]
            assert (expected["count"], expected["mesh_leaves"]) == (actual["count"], actual["mesh_leaves"]), name
            u = array.array("d", ()); v = array.array("d", ())
            u.frombytes((dirs[reference] / expected["file"]).read_bytes())
            v.frombytes((dirs[args.candidate] / actual["file"]).read_bytes())
            rel = math.sqrt(sum((x-y)**2 for x,y in zip(u,v))) / max(math.sqrt(sum(x*x for x in u)),1e-300)
            assert math.isfinite(rel) and rel <= 1e-6, (case, name, rel)
            errors[name] = rel
        summary.append(dict(case=case, field_relative_l2=errors,
                            max_original_residual=max(reports[args.candidate]["residuals"])))
        print(case, "PASS", errors, flush=True)
    (output / "connection-validation.json").write_text(
        json.dumps(dict(candidate=args.candidate, reference=reference, cases=summary), indent=2) + "\n")


if __name__ == "__main__":
    main()
