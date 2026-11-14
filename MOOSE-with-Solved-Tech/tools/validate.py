#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Check complete fields and residuals, without comparing execution times."""

import argparse, json, re
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file


def exodus(path):
    with netcdf_file(path, "r", mmap=False) as f:
        v = f.variables
        labels = [
            b"".join(row).split(b"\0", 1)[0].decode().strip()
            for row in v["name_nod_var"][:]
        ]
        fields = {
            label: v[f"vals_nod_var{i}"][:].copy() for i, label in enumerate(labels, 1)
        }
        # Validate auxiliary/element fields too, not just the headline variable.
        for name, var in v.items():
            if name.startswith("vals_elem_var") or name == "vals_glo_var":
                fields[name] = var[:].copy()
        coords = (
            v["coord"][:].copy()
            if "coord" in v
            else np.array(
                [v["coord" + axis][:] for axis in "xyz" if "coord" + axis in v]
            )
        )
        return v["time_whole"][:].copy(), fields, coords


def postprocessors(folder, cfg, expected_times):
    names = cfg["case_config"].get("postprocessors", [])
    if not names:
        return {}
    data = np.genfromtxt(folder / "physics.csv", delimiter=",", names=True)
    data = np.atleast_1d(data)
    if (
        "time" not in data.dtype.names
        or len(data) != len(expected_times)
        or not np.allclose(data["time"], expected_times, rtol=0, atol=1e-12)
    ):
        raise ValueError("Incomplete postprocessor time history")
    if any(name not in data.dtype.names for name in names):
        raise ValueError("Missing required physical postprocessor")
    values = {name: data[name] for name in data.dtype.names if name != "time"}
    if not all(np.isfinite(a).all() for a in values.values()):
        raise ValueError("Nonfinite physical postprocessor")
    return values


def validate(folder, reference=None):
    folder = Path(folder)
    cfg = json.loads((folder / "run.json").read_text())
    text = (folder / "application.log").read_text(errors="replace")
    if re.search(
        r"DIVERGED|Solve Did NOT Converge|Aborting as solve did not converge",
        text,
        re.I,
    ):
        raise ValueError("Application did not converge")
    files = list(folder.glob("physics*.e"))
    if len(files) != 1:
        raise ValueError(f"Expected one complete Exodus output; found {files}")
    times, fields, coords = exodus(files[0])
    expected = cfg["case_config"]["expected_output_times"]
    if times.shape != (len(expected),) or not np.allclose(
        times, expected, rtol=0, atol=1e-12
    ):
        raise ValueError("Incomplete physical time history")
    for name in cfg["case_config"]["fields"]:
        if name not in fields:
            raise ValueError("Missing field " + name)
    if not all(np.isfinite(a).all() for a in fields.values()):
        raise ValueError("Nonfinite physical field")
    pp = postprocessors(folder, cfg, times)
    checks = [
        float(x) for x in re.findall(r"MOOSE_PANGULU_CHECK .*?relative_l2=(\S+)", text)
    ]
    phases = re.findall(
        r"MOOSE_PANGULU_NATIVE phase=(\S+) seconds=(\S+) status=(\d+)", text
    )
    if cfg["solver"] == "pangulu":
        if "MOOSE_PANGULU_PROVIDER version=26.46 abi=2" not in text or not checks:
            raise ValueError("Missing PanguLU provider or residual evidence")
        if not all(np.isfinite(x) and x <= 1e-8 for x in checks):
            raise ValueError("Residual failed")
        if not phases or any(int(s) for _, _, s in phases):
            raise ValueError("Native phase failed")
    report = dict(
        completed=True,
        fields_finite=True,
        physical_times=times.tolist(),
        field_names=list(fields),
        postprocessor_names=list(pp),
        linear_rhs_checked=len(checks),
        max_relative_residual=max(checks, default=None),
        native_stage_seconds={},
    )
    for name, t, _ in phases:
        report["native_stage_seconds"][name] = report["native_stage_seconds"].get(
            name, 0
        ) + float(t)
    report["native_core_seconds"] = (
        sum(report["native_stage_seconds"].values()) if phases else None
    )
    if cfg["case"] == "hex20-heat":
        target = 200 * coords[0] + 100 * coords[1] + 200 * coords[2]
        error = float(
            np.linalg.norm(fields["temp"][-1] - target) / np.linalg.norm(target)
        )
        if not np.isfinite(error) or error > 1e-6:
            raise ValueError("Analytic temperature check failed")
        report["analytic_temperature_relative_l2"] = error
    if reference:
        ref_files = list(Path(reference).glob("physics*.e"))
        if len(ref_files) != 1:
            raise ValueError("Reference Exodus file missing or ambiguous")
        rt, rf, rc = exodus(ref_files[0])
        if (
            not np.array_equal(times, rt)
            or set(fields) != set(rf)
            or not np.array_equal(coords, rc)
        ):
            raise ValueError("Reference geometry, times or fields differ")
        errors = {}
        absolute_errors = {}
        for name, value in fields.items():
            if value.shape != rf[name].shape:
                raise ValueError("Field shape mismatch")
            denominator = np.linalg.norm(rf[name])
            difference = float(np.linalg.norm(value - rf[name]))
            # Relative error alone is undefined or misleading for a physical
            # component that is zero up to roundoff. Use an absolute RMS floor
            # in model units while retaining both measured errors in the report.
            bound = 1e-12 * np.sqrt(value.size) + 1e-6 * denominator
            if not np.isfinite(difference) or difference > bound:
                raise ValueError(f"Field mismatch: {name} absolute_l2={difference}")
            errors[name] = difference / denominator if denominator else None
            absolute_errors[name] = {
                "absolute_l2": difference,
                "absolute_rms": difference / np.sqrt(value.size),
                "max_absolute": float(np.max(np.abs(value - rf[name]))),
                "accepted_l2_bound": float(bound),
            }
        report["reference_field_relative_l2"] = errors
        report["reference_field_absolute_errors"] = absolute_errors
        report["reference_field_tolerances"] = {
            "relative_l2": 1e-6,
            "absolute_rms": 1e-12,
        }
        rp = postprocessors(Path(reference), cfg, times)
        if set(pp) != set(rp):
            raise ValueError("Reference postprocessor set differs")
        pp_errors = {}
        for name, value in pp.items():
            difference = float(np.linalg.norm(value - rp[name]))
            scale = float(np.linalg.norm(rp[name]))
            # Zero-valued constraints and symmetry components need an absolute
            # tolerance as well as a relative one, in the model's own units.
            if not np.isfinite(difference) or difference > 1e-12 + 1e-6 * scale:
                raise ValueError(f"Postprocessor mismatch: {name} {difference}")
            pp_errors[name] = {
                "absolute_l2": difference,
                "relative_l2": difference / scale if scale else None,
            }
        report["reference_postprocessor_errors"] = pp_errors
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder")
    p.add_argument("--reference")
    a = p.parse_args()
    print(json.dumps(validate(a.folder, a.reference), indent=2, allow_nan=False))
