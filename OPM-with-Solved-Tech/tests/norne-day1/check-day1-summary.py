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

ll = data["routes"]["linglongamg"]
ll_sdk = ll["linglongamg_26_46"]
amgx = ll["amgx_2_5_0_baseline"]
assert ll_sdk["total_s"] == 33.73
assert amgx["total_s"] == 63.46
assert round(amgx["total_s"] / ll_sdk["total_s"], 2) == ll["speedup_vs_amgx"]

xk = data["routes"]["xianzhikrylov"]
xk_sdk = xk["xianzhikrylov_26_48_bsr3"]
gko = xk["ginkgo_1_11_0_csr_baseline"]
assert xk_sdk["total_s"] == 74.77
assert gko["total_s"] == 158.85
assert xk_sdk["true_relative_l2"] < 0.005
assert gko["true_relative_l2"] < 0.005
assert round(gko["total_s"] / xk_sdk["total_s"], 2) == xk["speedup_vs_ginkgo"]

print("check-day1-summary: OK (2 routes verified against recorded reference summary)")
