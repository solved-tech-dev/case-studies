"""Exact parent-element h-refinement of the native biphasic quarter cylinder.

Keep rounded original geometry, ALL z-displacement constraints (plane strain),
pressure boundary sets and original oriented inner pressure facets. No analytic
circle replacement, material/control/plot changes or unreported new loading.
"""

import copy
import itertools
import math
import xml.etree.ElementTree as ET


def nid(i, j, k, nr, nt):
    return 1 + i + (nr + 1) * (j + (nt + 1) * k)


def elements(nr, nt, nz):
    return {
        1
        + i
        + nr
        * (j + nt * k): [
            nid(i, j, k, nr, nt),
            nid(i + 1, j, k, nr, nt),
            nid(i + 1, j + 1, k, nr, nt),
            nid(i, j + 1, k, nr, nt),
            nid(i, j, k + 1, nr, nt),
            nid(i + 1, j, k + 1, nr, nt),
            nid(i + 1, j + 1, k + 1, nr, nt),
            nid(i, j + 1, k + 1, nr, nt),
        ]
        for k in range(nz)
        for j in range(nt)
        for i in range(nr)
    }


def sets(nr, nt, nz):
    return {
        "FixedDisplacement01": {
            nid(i, nt, k, nr, nt) for k in range(nz + 1) for i in range(nr + 1)
        },
        "FixedDisplacement02": {
            nid(i, 0, k, nr, nt) for k in range(nz + 1) for i in range(nr + 1)
        },
        "FixedDisplacement03": set(range(1, (nr + 1) * (nt + 1) * (nz + 1) + 1)),
        "FixedFluidPressure01": {
            nid(nr, j, k, nr, nt) for k in range(nz + 1) for j in range(nt + 1)
        },
        "PrescribedBC06": {
            nid(0, j, k, nr, nt) for k in range(nz + 1) for j in range(nt + 1)
        },
    }


def facets(nr, nt, nz):
    return {
        1
        + j
        + nt
        * k: [
            nid(0, j + 1, k, nr, nt),
            nid(0, j, k, nr, nt),
            nid(0, j, k + 1, nr, nt),
            nid(0, j + 1, k + 1, nr, nt),
        ]
        for k in range(nz)
        for j in range(nt)
    }


def build(source, output, f):
    assert isinstance(f, int) and 1 <= f <= 12
    old = ET.fromstring(source.read_bytes().decode("iso-8859-1"))
    assert old.find("Module").get("type") == "biphasic"
    assert old.find("Control/solver/linear_solver") is None
    ET.SubElement(old.find("Control/solver"), "linear_solver", type="febio_gpu_lu")
    assert (
        old.findtext("Control/time_steps") == "10"
        and old.findtext("Control/step_size") == "0.1"
    )
    assert old.find("Material/material/solid").get("type") == "Holmes-Mow"
    assert old.find("Material/material/permeability").get("type") == "perm-Holmes-Mow"
    mesh = old.find("Mesh")
    nodes = mesh.find("Nodes")
    cells = mesh.find("Elements")
    assert cells.get("type") == "hex8" and cells.get("name") == "Part1"
    points = {int(n.get("id")): tuple(map(float, n.text.split(","))) for n in nodes}
    assert set(points) == set(range(1, 379)) and all(
        math.isfinite(x) for p in points.values() for x in p
    )
    assert {
        int(e.get("id")): list(map(int, e.text.split(","))) for e in cells
    } == elements(20, 8, 1)
    assert {
        s.get("name"): set(map(int, s.text.split(","))) for s in mesh.findall("NodeSet")
    } == sets(20, 8, 1)
    assert (
        len(mesh.findall("Surface")) == 1
        and mesh.find("Surface").get("name") == "PressureLoad1"
    )
    assert all(e.tag == "quad4" for e in mesh.find("Surface"))
    assert {
        int(e.get("id")): list(map(int, e.text.split(",")))
        for e in mesh.find("Surface")
    } == facets(20, 8, 1)
    assert (
        len(mesh.findall("ElementSet")) == 1
        and mesh.find("ElementSet").get("name") == "unnamed_set01"
    )
    assert list(map(int, mesh.find("ElementSet").text.split(","))) == [2]
    assert not any(
        "unnamed_set01" in ET.tostring(s, encoding="unicode")
        for s in old
        if s.tag != "Mesh"
    )
    # Exact extrusion and parent-map positivity over each complete reference
    # square. A bilinear2D determinant is affine in u,v, so corner minima bound
    # the whole square. The z-map has unit positive derivative.
    parent_dets = []
    for j in range(9):
        for i in range(21):
            p, q = points[nid(i, j, 0, 20, 8)], points[nid(i, j, 1, 20, 8)]
            assert p[:2] == q[:2] and p[2] == 0 and q[2] == 1
    for j in range(8):
        for i in range(20):
            q = [
                points[nid(i + di, j + dj, 0, 20, 8)]
                for di, dj in [(0, 0), (1, 0), (1, 1), (0, 1)]
            ]
            for u, v in itertools.product([0.0, 1.0], repeat=2):
                du = [
                    (1 - v) * (q[1][c] - q[0][c]) + v * (q[2][c] - q[3][c])
                    for c in range(2)
                ]
                dv = [
                    (1 - u) * (q[3][c] - q[0][c]) + u * (q[2][c] - q[1][c])
                    for c in range(2)
                ]
                parent_dets.append(du[0] * dv[1] - du[1] * dv[0])
    assert min(parent_dets) > 0
    nr, nt, nz = 20 * f, 8 * f, f
    refined = {}

    def parent_position(index, maximum):
        cell = min(index // f, maximum - 1)
        return cell, (index - cell * f) / f

    for k in range(nz + 1):
        iz, w = parent_position(k, 1)
        for j in range(nt + 1):
            it, v = parent_position(j, 8)
            for i in range(nr + 1):
                ir, u = parent_position(i, 20)
                if i % f == j % f == k % f == 0:
                    value = points[nid(i // f, j // f, k // f, 20, 8)]
                else:
                    value = tuple(
                        sum(
                            points[nid(ir + di, it + dj, iz + dk, 20, 8)][c]
                            * (
                                (u if di else 1 - u)
                                * (v if dj else 1 - v)
                                * (w if dk else 1 - w)
                            )
                            for di, dj, dk in itertools.product([0, 1], repeat=3)
                        )
                        for c in range(3)
                    )
                refined[nid(i, j, k, nr, nt)] = value
    assert all(
        refined[nid(i * f, j * f, k * f, nr, nt)] == points[nid(i, j, k, 20, 8)]
        for k in range(2)
        for j in range(9)
        for i in range(21)
    )
    out = copy.deepcopy(old)
    m = out.find("Mesh")
    n = m.find("Nodes")
    e = m.find("Elements")
    for item in list(n):
        n.remove(item)
    for item in list(e):
        e.remove(item)
    for key, p in sorted(refined.items()):
        ET.SubElement(n, "node", id=str(key)).text = ",".join(
            format(x, ".17g") for x in p
        )
    new_cells = elements(nr, nt, nz)
    for key, c in sorted(new_cells.items()):
        ET.SubElement(e, "elem", id=str(key)).text = ",".join(map(str, c))
    boundaries = sets(nr, nt, nz)
    for s in m.findall("NodeSet"):
        s.text = ",".join(map(str, sorted(boundaries[s.get("name")])))
    surface = m.find("Surface")
    for item in list(surface):
        surface.remove(item)
    for key, c in sorted(facets(nr, nt, nz).items()):
        ET.SubElement(surface, "quad4", id=str(key)).text = ",".join(map(str, c))
    children = [
        1 + i + nr * (j + nt * k)
        for k in range(f)
        for j in range(f)
        for i in range(f, 2 * f)
    ]
    assert len(children) == f**3
    m.find("ElementSet").text = ",".join(map(str, children))
    unchanged = {
        a.tag: ET.tostring(a) == ET.tostring(b)
        for a, b in zip(old, out)
        if a.tag != "Mesh"
    }
    assert all(unchanged.values()) and len(old) == len(out)
    # Native small measured1032 verifies the free-equation estimate: all uz
    # fixed, symmetry ux/uy fixed; zero outer p fixed, prescribed inner p retained.
    estimate = (
        3 * len(refined)
        - len(boundaries["FixedDisplacement01"])
        - len(boundaries["FixedDisplacement02"])
        - len(boundaries["FixedFluidPressure01"])
    )
    if f == 1:
        assert estimate == 1032
    ET.ElementTree(out).write(output, encoding="iso-8859-1", xml_declaration=True)
