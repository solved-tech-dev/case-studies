"""Refine selected axis-aligned HEX8 verification models in parent coordinates."""

from itertools import product
import xml.etree.ElementTree as ET
import numpy as np

hex8 = np.array(
    [
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    ]
)


def numbers(text):
    return [int(i.strip()) for i in text.split(",") if i.strip()]


def build(source, output):
    factors = np.asarray((45, 45, 1))
    tree = ET.parse(source)
    mesh = tree.find("Mesh")
    assert tree.find("MeshData") is None and set(x.tag for x in mesh) <= {
        "Nodes",
        "Elements",
        "NodeSet",
        "Surface",
        "ElementSet",
    }
    preserved = [ET.tostring(x) for x in tree.getroot() if x.tag != "Mesh"]
    assert len(mesh.findall("Elements")) == len(mesh.findall("Nodes")) == 1
    nodes, elems = mesh.find("Nodes"), mesh.find("Elements")
    assert elems.get("type") == "hex8"
    xyz = np.array([list(map(float, a.text.split(","))) for a in nodes])
    node_ids = [int(a.get("id")) for a in nodes]
    axes = [np.unique(xyz[:, i]) for i in range(3)]
    old_dims = np.array([len(a) - 1 for a in axes])
    ids = np.arange(1, 1 + len(nodes)).reshape(tuple(old_dims + 1))
    ijk = {
        idx: np.array([np.searchsorted(a, x) for a, x in zip(axes, coord)])
        for idx, coord in zip(node_ids, xyz)
    }
    assert len(nodes) == len({tuple(i) for i in ijk.values()})
    for idx, coord in ijk.items():
        ids[tuple(coord)] = idx
    grids = [
        np.r_[
            [(u + (v - u) * j / f) for u, v in zip(a, a[1:]) for j in range(f)], a[-1]
        ]
        for a, f in zip(axes, factors)
    ]
    for new, old, f in zip(grids, axes, factors):
        assert np.array_equal(new[::f], old) and np.all(np.diff(new) > 0)
    dims = old_dims * factors
    new_ids = np.arange(1, np.prod(dims + 1) + 1).reshape(tuple(dims + 1))
    new_eids = np.arange(1, np.prod(dims) + 1).reshape(tuple(dims))
    origins = {}
    for e in elems:
        conn = numbers(e.text)
        origin = ijk[conn[0]]
        assert np.array_equal([ijk[n] for n in conn], origin + hex8)
        origins[int(e.get("id"))] = origin
    assert (
        len(origins) == len({tuple(v) for v in origins.values()}) == np.prod(old_dims)
    )
    for block in mesh.findall("NodeSet"):
        original = numbers(block.text)
        selected = set(original)
        mask = np.isin(ids, list(selected))
        if mask.all():
            chosen = new_ids.ravel()
        else:
            coarse = np.zeros_like(ids, dtype=bool)
            fine = np.zeros_like(new_ids, dtype=bool)
            for axis, edge in product(range(3), (0, -1)):
                sl = [slice(None)] * 3
                sl[axis] = edge
                sl = tuple(sl)
                if mask[sl].all():
                    coarse[sl] = True
                    fine[sl] = True
            if not np.array_equal(mask, coarse):
                # mi06 combines a fixed face with individual corner nodes.
                # The thickness mesh is unchanged; keep those nodal constraints at
                # their original positions, without expanding them to another face.
                assert factors[2] == 1
                extra = selected - set(ids[coarse])
                assert extra and extra <= set(
                    ids[:: old_dims[0], :: old_dims[1], :].ravel()
                )
                chosen = np.union1d(
                    new_ids[fine],
                    [int(new_ids[tuple(ijk[n] * factors)]) for n in extra],
                )
            else:
                chosen = new_ids[fine]
        block.text = ",".join(map(str, chosen))
    for surface in mesh.findall("Surface"):
        old = list(surface)
        for face in old:
            surface.remove(face)
        count = 0
        for face in old:
            assert face.tag == "quad4"
            a, b, c, d = [ijk[n] for n in numbers(face.text)]
            u, v = b - a, d - a
            assert (
                np.abs(u).sum() == np.abs(v).sum() == 1
                and u @ v == 0
                and np.array_equal(c, a + u + v)
            )
            for i, j in product(
                range(int(abs(u) @ factors)), range(int(abs(v) @ factors))
            ):
                base = a * factors + i * u + j * v
                quad = [new_ids[tuple(base + t)] for t in (0, u, u + v, v)]
                count += 1
                ET.SubElement(surface, "quad4", id=str(count)).text = ",".join(
                    map(str, quad)
                )
    for block in mesh.findall("ElementSet"):
        children = []
        for n in numbers(block.text):
            base = origins[n] * factors
            sl = tuple(slice(int(k), int(k + f)) for k, f in zip(base, factors))
            children.extend(new_eids[sl].ravel())
        block.text = ",".join(map(str, children))
    del nodes[:], elems[:]
    for idx, coord in enumerate(product(*grids), 1):
        ET.SubElement(nodes, "node", id=str(idx)).text = ",".join(
            format(x, ".17g") for x in coord
        )
    offsets = np.array([new_ids[tuple(row)] - 1 for row in hex8])
    for eid, base in enumerate(new_ids[:-1, :-1, :-1].ravel(), 1):
        ET.SubElement(elems, "elem", id=str(eid)).text = ",".join(
            map(str, base + offsets)
        )
    assert preserved == [ET.tostring(x) for x in tree.getroot() if x.tag != "Mesh"]
    ET.SubElement(tree.find("Control/solver"), "linear_solver", type="febio_gpu_lu")
    tree.write(output, encoding="iso-8859-1", xml_declaration=True)
