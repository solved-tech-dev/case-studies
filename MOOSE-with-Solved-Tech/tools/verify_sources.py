#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Verify the preserved upstream files and the release's bilingual entry points."""

import hashlib, json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
lock = json.loads((root / "upstream.lock.json").read_text())
for row in lock["files"]:
    path = root / row["file"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        raise SystemExit("Upstream input changed: " + row["file"])
for name in [
    "README",
    "docs/BUILD-AND-RUN",
    "docs/CASES",
    "docs/INTEGRATION",
    "docs/VALIDATION",
]:
    for suffix in [".md", ".zh-CN.md"]:
        if not (root / (name + suffix)).is_file():
            raise SystemExit("Missing documentation: " + name + suffix)
print(
    f'{len(lock["files"])} original inputs/meshes intact; English and Chinese entry points present.'
)
