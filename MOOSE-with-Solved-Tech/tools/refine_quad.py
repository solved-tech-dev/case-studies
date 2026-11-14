#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Subdivide the original QUAD4 faces without surface reprojection."""

import argparse, json
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file


def names(var):
    return [b"".join(row).split(b"\0", 1)[0].decode().strip() for row in var[:]]


def refine(source, target, factor):
    if factor < 1:
        raise ValueError("Subdivision factor must be positive")
    with netcdf_file(source, "r", mmap=False) as f:
        v = f.variables
        n = f.dimensions["num_nodes"]
        coords = np.zeros((n, 3))
        if "coord" in v:
            coords[:, : v["coord"].shape[0]] = v["coord"][:].T
        else:
            for j, axis in enumerate("xyz"):
                if "coord" + axis in v:
                    coords[:, j] = v["coord" + axis][:]
        blocks = [int(x) for x in v["eb_prop1"][:]]
        parents = []
        for i, block in enumerate(blocks, 1):
            q = v[f"connect{i}"][:]
            if q.shape[1] != 4:
                raise ValueError("Only QUAD4 seed meshes are supported")
            parents.extend((block, tuple(int(x) - 1 for x in row)) for row in q)
        edges = {}
        point_sets = {}
        boundary_names = {}
        if "ns_prop1" in v:
            labels = names(v["ns_names"])
            for i, b in enumerate(v["ns_prop1"][:], 1):
                b = int(b)
                ns = set(int(x) - 1 for x in v[f"node_ns{i}"][:])
                boundary_names[b] = labels[i - 1]
                if len(ns) == 1:
                    point_sets[b] = next(iter(ns))
                for _, q in parents:
                    for a, c in zip(q, q[1:] + q[:1]):
                        if a in ns and c in ns:
                            edges[tuple(sorted((a, c)))] = b
        if "ss_prop1" in v:
            for i, b in enumerate(v["ss_prop1"][:], 1):
                b = int(b)
                boundary_names[b] = str(b)
                for elem, side in zip(v[f"elem_ss{i}"][:], v[f"side_ss{i}"][:]):
                    q = parents[int(elem) - 1][1]
                    s = int(side) - 1
                    edges[tuple(sorted((q[s], q[(s + 1) % 4])))] = b
    nodes = []
    lookup = {}
    quads = []
    lines = {}

    def node_id(x):
        key = tuple(np.round(x, 12))
        if key not in lookup:
            lookup[key] = len(nodes) + 1
            nodes.append(tuple(x))
        return lookup[key]

    for block, q in parents:
        a, b, c, e = coords[list(q)]
        ids = np.empty((factor + 1, factor + 1), dtype=np.int64)
        for j in range(factor + 1):
            t = j / factor
            for i in range(factor + 1):
                s = i / factor
                ids[j, i] = node_id(
                    (1 - s) * (1 - t) * a
                    + s * (1 - t) * b
                    + s * t * c
                    + (1 - s) * t * e
                )
        for j in range(factor):
            for i in range(factor):
                quads.append(
                    (
                        block,
                        int(ids[j, i]),
                        int(ids[j, i + 1]),
                        int(ids[j + 1, i + 1]),
                        int(ids[j + 1, i]),
                    )
                )
        for side, seq in enumerate(
            [ids[0, :], ids[:, -1], ids[-1, ::-1], ids[::-1, 0]]
        ):
            edge = tuple(sorted((q[side], q[(side + 1) % 4])))
            if edge in edges:
                for x, y in zip(seq[:-1], seq[1:]):
                    lines[tuple(sorted((int(x), int(y))))] = edges[edge]
    points = {b: node_id(coords[p]) for b, p in point_sets.items()}
    with Path(target).open("x") as out:
        out.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        physical = [(2, b, str(b)) for b in blocks] + [
            (0 if b in points else 1, b, n or str(b)) for b, n in boundary_names.items()
        ]
        out.write(f"$PhysicalNames\n{len(physical)}\n")
        for dim, b, label in physical:
            out.write(f'{dim} {b} "{label}"\n')
        out.write(f"$EndPhysicalNames\n$Nodes\n{len(nodes)}\n")
        for i, x in enumerate(nodes, 1):
            out.write(f"{i} {x[0]:.17g} {x[1]:.17g} {x[2]:.17g}\n")
        out.write(f"$EndNodes\n$Elements\n{len(quads)+len(lines)+len(points)}\n")
        i = 0
        for b, x in points.items():
            i += 1
            out.write(f"{i} 15 2 {b} {b} {x}\n")
        for (x, y), b in lines.items():
            i += 1
            out.write(f"{i} 1 2 {b} {b} {x} {y}\n")
        for b, a, c, e, g in quads:
            i += 1
            out.write(f"{i} 3 2 {b} {b} {a} {c} {e} {g}\n")
        out.write("$EndElements\n")
    return dict(
        seed_nodes=n,
        seed_elements=len(parents),
        factor=factor,
        nodes=len(nodes),
        elements=len(quads),
        boundary_names=boundary_names,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source")
    p.add_argument("target")
    p.add_argument("--factor", type=int, required=True)
    a = p.parse_args()
    print(json.dumps(refine(a.source, a.target, a.factor), indent=2))
