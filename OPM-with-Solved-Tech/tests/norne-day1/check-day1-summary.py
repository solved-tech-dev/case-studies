#!/usr/bin/env python3
"""Validate structural integrity of examples/norne/reference-day1-summary.json."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
summary_path = ROOT / "examples" / "norne" / "reference-day1-summary.json"
data = json.loads(summary_path.read_text(encoding="utf-8"))

assert data["active_cells"] == 2843196
assert data["main_scalar_unknowns"] == 8529588
assert set(data["routes"].keys()) == {"linglongamg", "xianzhikrylov"}

ll = data["routes"]["linglongamg"]["day1"]["linglongamg"]
assert ll["linear_total_s"] == round(ll["linear_setup_s"] + ll["linear_solve_s"], 2)

xk = data["routes"]["xianzhikrylov"]["public_sdk_rtx5090_verification"]
assert xk["smoke_rel1"] < 1e-12
assert xk["max_converged_true_relative_residual"] < 0.005

print("check-day1-summary: OK (2 routes verified against recorded reference summary)")
