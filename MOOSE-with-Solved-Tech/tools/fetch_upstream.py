#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Fetch the pinned MOOSE source into a new directory; leave existing installs alone."""

import argparse, json, subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument("destination")
a = p.parse_args()
dest = Path(a.destination).resolve()
lock = json.loads((root / "upstream.lock.json").read_text())
if dest.exists():
    p.error("Destination exists; use a new directory")
subprocess.run(
    [
        "git",
        "clone",
        "--filter=blob:none",
        "--no-checkout",
        lock["repository"],
        str(dest),
    ],
    check=True,
)
subprocess.run(
    ["git", "-C", str(dest), "checkout", "--detach", lock["commit"]], check=True
)
print(dest)
print(
    "Next: build/install compatible PETSc and libMesh using the MOOSE instructions. See docs/BUILD-AND-RUN.md."
)
