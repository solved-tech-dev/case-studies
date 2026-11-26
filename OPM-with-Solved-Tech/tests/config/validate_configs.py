#!/usr/bin/env python3
"""Validate public JSON configurations, manifests, and absence of private paths."""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PRIVATE_RE = re.compile(
    r"(?:^|[\s\"'=:,;])(" + "|".join([
        "/ho" + "me/",
        "/mnt/" + "NAS",
        "/wo" + "rk2/",
        "/opt/" + "spack/",
        "/var/" + "tmp/",
    ]) + r")"
)

REQUIRED_JSON = [
    "examples/configs/linglongamg.json",
    "examples/configs/linglong-gmres-amg.json",
    "examples/configs/xianzhikrylov.json",
    "examples/configs/cpr.json",
    "examples/configs/linglong-vcycle.json",
    "examples/norne/reference-day1-summary.json",
    "manifests/source.json",
    "manifests/sdk-contract.json",
]

for rel in REQUIRED_JSON:
    p = ROOT / rel
    if not p.is_file():
        raise SystemExit(f"Missing JSON file: {rel}")
    raw = p.read_text(encoding="utf-8")
    if PRIVATE_RE.search(raw):
        raise SystemExit(f"Private path found in {rel}")
    data = json.loads(raw)
    if rel == "examples/configs/linglongamg.json":
        assert data["solver"] == "bicgstab"
        assert data["preconditioner"]["type"] == "cprw"
        coarse_pc = data["preconditioner"]["coarsesolver"]["preconditioner"]
        assert coarse_pc["type"] == "linglongamg"
        assert not Path(coarse_pc["config_file"]).is_absolute()
    elif rel == "examples/configs/xianzhikrylov.json":
        assert data["solver"] == "bicgstab"
        assert data["preconditioner"]["type"] == "paroverilu0"
        assert data["preconditioner"]["relaxation"] == 1.0
    elif rel == "examples/configs/cpr.json":
        assert data["solver"] == "bicgstab"
        assert data["preconditioner"]["type"] == "cprw"
        assert data["preconditioner"]["weight_type"] == "trueimpes"
        assert data["preconditioner"]["add_wells"] is True
        assert data["preconditioner"]["use_well_weights"] is False
    elif rel == "examples/configs/linglong-vcycle.json":
        assert data["solve"]["method"] == "amg"
        assert data["solve"]["cycle"] == "v"
        assert data["solve"]["max_iters"] == 1
    elif rel == "manifests/source.json":
        assert data["public_repository_policy"]["sdk_binaries_included"] is False
        assert data["public_repository_policy"]["licenses_included"] is False
    elif rel == "manifests/sdk-contract.json":
        assert set(data["routes"].keys()) == {"linglongamg", "xianzhikrylov", "cpr"}

print(f"validate_configs: OK ({len(REQUIRED_JSON)} JSON files verified)")
