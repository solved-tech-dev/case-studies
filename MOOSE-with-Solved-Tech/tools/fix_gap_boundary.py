# SPDX-License-Identifier: AGPL-3.0-only
"""Restore the seed mesh's one-sided contact surfaces after Gmsh conversion."""

import numpy as np
from scipy.io import netcdf_file


def fix(source, target, factor):
    with netcdf_file(source, "r", mmap=False) as f:
        active = []
        offset = 0
        for i, block in enumerate(f.variables["eb_prop1"][:], 1):
            count = f.dimensions[f"num_el_in_blk{i}"]
            if block in (1, 3):
                active.extend(range(offset + 1, offset + count + 1))
            offset += count
        masks = {}
        report = {}
        for i, boundary in enumerate(f.variables["ss_prop1"][:], 1):
            old = f.variables[f"elem_ss{i}"][:]
            mask = np.isin(old, active)
            if not mask.all():
                masks[f"num_side_ss{i}"] = mask
            if mask.sum() != 10 * factor:
                raise ValueError("Contact boundary count changed unexpectedly")
            report[int(boundary)] = dict(before=len(old), after=int(mask.sum()))
        with netcdf_file(target, "w", version=2) as out:
            for key, value in f._attributes.items():
                setattr(out, key, value)
            for name, size in sorted(
                f.dimensions.items(), key=lambda item: item[1] is not None
            ):
                out.createDimension(
                    name, int(masks[name].sum()) if name in masks else size
                )
            for name, var in f.variables.items():
                dest = out.createVariable(name, var.typecode(), var.dimensions)
                for key, value in var._attributes.items():
                    setattr(dest, key, value)
                values = var[:].copy()
                for axis, dim in enumerate(var.dimensions):
                    if dim in masks:
                        values = np.compress(masks[dim], values, axis=axis)
                dest[:] = values
    return report
