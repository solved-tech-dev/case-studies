#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
import sys

FIELDS = (
    "case", "backend", "requested_steps", "cells", "exit_code", "wall_s", "pressure_calls",
    "pressure_calls_converged", "pressure_iterations", "max_relative_rhs",
    "adapter_setup_ms", "adapter_solve_ms", "sdk_setup_ms", "sdk_solve_ms",
)

def read_summary(directory):
    summary = json.loads((directory / "run.json").read_text())
    pressure = directory / "results" / "linglongamg_pressure.csv"
    if not pressure.is_file():
        pressure = directory / "results" / "linglong_pressure.csv"
    if pressure.is_file():
        with pressure.open(newline="") as source:
            rows = list(csv.DictReader(source))
        if rows:
            summary.update(
                pressure_calls=len(rows),
                pressure_calls_converged=sum(int(row["converged"]) for row in rows),
                pressure_iterations=sum(int(row["iterations"]) for row in rows),
                max_relative_rhs=max(float(row["relative_rhs"]) for row in rows),
                adapter_setup_ms=sum(float(row["setup_ms"]) for row in rows),
                adapter_solve_ms=sum(float(row["solve_ms"]) for row in rows),
                sdk_setup_ms=sum(float(row["sdk_setup_ms"]) for row in rows),
                sdk_solve_ms=sum(float(row["sdk_solve_ms"]) for row in rows),
            )
    return summary

def main():
    parser = argparse.ArgumentParser(description="Summarize PHengLEI run folders as CSV.")
    parser.add_argument("runs", nargs="+", type=Path)
    args = parser.parse_args()
    writer = csv.DictWriter(sys.stdout, FIELDS, extrasaction="ignore")
    writer.writeheader()
    for directory in args.runs:
        writer.writerow(read_summary(directory))

if __name__ == "__main__":
    main()
