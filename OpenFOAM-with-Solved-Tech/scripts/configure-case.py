#!/usr/bin/env python3
"""Select equations for one solver route and set its SIMPLE step count."""

import argparse
import json
from pathlib import Path
import re


def field_span(text, name):
    entry = re.search(rf"(?m)^    {re.escape(name)}\s*\n    \{{", text)
    if entry is None:
        raise ValueError(f"Missing solver entry for {name}")
    start = text.index("{", entry.start())
    depth = 0
    for end in range(start, len(text)):
        depth += (text[end] == "{") - (text[end] == "}")
        if depth == 0:
            return start, end + 1
    raise ValueError(f"Unclosed solver entry for {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("backend", choices=("linglongamg", "xianzhikrylov"))
    parser.add_argument("case", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("steps", type=int)
    parser.add_argument("--ordering", type=Path)
    args = parser.parse_args()
    if args.steps < 1:
        parser.error("steps must be positive")
    if not args.config.is_file():
        parser.error(f"missing SDK configuration: {args.config}")
    if args.ordering is not None:
        if args.backend != "linglongamg":
            parser.error("--ordering applies to linglongamg")
        for rank in range(20):
            file = args.ordering / f"processor{rank}.bin"
            if not file.is_file() or not file.stat().st_size:
                parser.error(f"missing ordering data: {file}")

    fv_solution = args.case / "system" / "fvSolution"
    text = fv_solution.read_text()
    config = json.dumps(str(args.config.resolve()), ensure_ascii=False)
    if args.backend == "linglongamg":
        executable = "timedSimpleFoam"
        start, end = field_span(text, "p")
        native_pressure = text[start:end]
        phi_start, phi_end = field_span(text, "Phi")
        text = text[:phi_start] + native_pressure + text[phi_end:]
        ordering = ""
        if args.ordering is not None:
            pattern = json.dumps(str(args.ordering.resolve() / "processor%r.bin"), ensure_ascii=False)
            ordering = f"        permutationFile {pattern};\n"
        pressure = (
            "{\n"
            "        solver          FoamPressure;\n"
            "        backend         linglongamg;\n"
            f"        configFile      {config};\n"
            f"{ordering}"
            "        device          0;\n"
            "        tolerance       0;\n"
            "        relTol          0;\n"
            "        log             0;\n"
            "    }"
        )
        text = text[:start] + pressure + text[end:]
    else:
        executable = "benchSimpleFoam"
        transport = (
            "{\n"
            "        solver          BenchLinear;\n"
            "        backend         xianzhi;\n"
            f"        configFile      {config};\n"
            "        tolerance       1e-8;\n"
            "        relTol          0;\n"
            "        maxIter         2000;\n"
            "        log             0;\n"
            "    }"
        )
        for field in ("U", "k", "omega"):
            start, end = field_span(text, field)
            text = text[:start] + transport + text[end:]
    fv_solution.write_text(text)

    control = args.case / "system" / "controlDict"
    text = re.sub(r"(?m)^application\s+\w+;", f"application     {executable};", control.read_text())
    text = re.sub(r"(?m)^endTime\s+\d+;", f"endTime         {args.steps};", text)
    control.write_text(text)


if __name__ == "__main__":
    main()
