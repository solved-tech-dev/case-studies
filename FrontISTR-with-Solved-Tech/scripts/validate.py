"""Check solve history and final FrontISTR v2 ASCII nodal result files."""
import json
import math
from pathlib import Path
import re


def field_norms(model, folder):
    path = folder / model["final_result"]
    fields = model["fields"]
    with path.open() as stream:
        for line in stream:
            if line.strip() == "*data":
                break
        counts = list(map(int, stream.readline().split()))
        if counts != [model["nodes"], model["elements"]]:
            raise ValueError(f"Unexpected mesh counts in {path}: {counts}")
        if list(map(int, stream.readline().split())) != [len(fields), 0]:
            raise ValueError("Unexpected node/element field counts")
        widths = list(map(int, stream.readline().split()))
        names = [stream.readline().strip() for _ in widths]
        if list(zip(names, widths)) != [(f["name"], f["width"]) for f in fields]:
            raise ValueError("Unexpected field names or widths")
        sums = [0.0] * len(fields)
        for node in range(1, model["nodes"] + 1):
            if int(stream.readline()) != node:
                raise ValueError(f"Expected node {node}")
            row = []
            while len(row) < sum(widths):
                line = stream.readline()
                if not line:
                    raise ValueError(f"Truncated field at node {node}")
                row.extend(map(float, line.split()))
            if len(row) != sum(widths) or not all(map(math.isfinite, row)):
                raise ValueError(f"Invalid field values at node {node}")
            offset = 0
            for i, width in enumerate(widths):
                sums[i] += sum(x*x for x in row[offset:offset+width])
                offset += width
        if stream.read().strip():
            raise ValueError(f"Unexpected trailing records in {path}")
    return [math.sqrt(s) for s in sums]


def validate(model, folder):
    log = (folder / "application.log").read_text(errors="replace")
    if "FrontISTR Completed" not in log:
        raise ValueError("FrontISTR did not report completion")
    records = re.findall(r"^FRONTISTR_SOLVE call=(\d+) n=(\d+) reuse_analysis=(\d) reuse_factor=(\d) residual=(\S+)$",
                         log, re.MULTILINE)
    if len(records) != model["calls"]:
        raise ValueError("Unexpected number of PanguLU solves")
    residuals = []
    for index, (call, size, reuse_analysis, reuse_factor, raw) in enumerate(records):
        if (int(call), int(size), int(reuse_analysis), int(reuse_factor)) != (
                index + 1, model["equations"], int(index > 0), 0):
            raise ValueError(f"Unexpected system size or reuse at solve {call}")
        residual = float(raw)
        if not math.isfinite(residual) or not 0 <= residual <= 1e-8:
            raise ValueError(f"Residual check failed at solve {call}: {residual}")
        residuals.append(residual)
    norms = field_norms(model, folder)
    for field, value in zip(model["fields"], norms):
        if not math.isclose(value, field["reference_l2"], rel_tol=1e-6, abs_tol=1e-12):
            raise ValueError(f"Unexpected final field norm for {field['name']}: {value}")
    return {"case": model["case"], "calls": len(records), "max_residual": max(residuals),
            "final_field_l2": {f["name"]: v for f, v in zip(model["fields"], norms)}}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("conrod", "hinge", "heat"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    print(json.dumps(validate(json.loads((root / "cases" / args.case / "model.json").read_text()), args.output), indent=2))
