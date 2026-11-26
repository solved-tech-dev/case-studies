#!/usr/bin/env python3
"""Validate structural integrity of examples/norne/reference-day1-summary.json."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
summary_path = ROOT / "examples" / "norne" / "reference-day1-summary.json"
data = json.loads(summary_path.read_text(encoding="utf-8"))

assert data["active_cells"] == 2843196
assert data["main_scalar_unknowns"] == 8529588
assert set(data["routes"].keys()) == {"linglongamg", "xianzhikrylov", "cpr"}

ll = data["routes"]["linglongamg"]["day1"]["linglongamg"]
assert ll["linear_total_s"] == round(ll["linear_setup_s"] + ll["linear_solve_s"], 2)

xk = data["routes"]["xianzhikrylov"]["day1_no_output_timing"]
assert xk["xianzhikrylov_bsr3_ilu0"]["max_true_relative_residual"] < 0.005

cpr = data["routes"]["cpr"]["day1_paired_no_output"]
assert cpr["physical_summary_exact_match"]["FPR"] == 271.699219
assert cpr["gpu_cpr"]["linear_s"] < cpr["cpu_1rank_hypre_cpu"]["linear_s"]

print("check-day1-summary: OK (3 routes verified against recorded reference summary)")
