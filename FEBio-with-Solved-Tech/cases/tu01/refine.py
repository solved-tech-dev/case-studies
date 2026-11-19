"""Exact original Cartesian HEX8 subdivision and boundary-plane preservation."""

import copy
import xml.etree.ElementTree as E


def node(i, j, k, nx, ny, nz):
    return 1 + (i * (ny + 1) + j) * (nz + 1) + k


def cells(nx, ny, nz):
    return {
        1
        + (i * ny + j) * nz
        + k: [
            node(i, j, k, nx, ny, nz),
            node(i + 1, j, k, nx, ny, nz),
            node(i + 1, j + 1, k, nx, ny, nz),
            node(i, j + 1, k, nx, ny, nz),
            node(i, j, k + 1, nx, ny, nz),
            node(i + 1, j, k + 1, nx, ny, nz),
            node(i + 1, j + 1, k + 1, nx, ny, nz),
            node(i, j + 1, k + 1, nx, ny, nz),
        ]
        for i in range(nx)
        for j in range(ny)
        for k in range(nz)
    }


def build(source, output, factor):
    assert 1 <= factor <= 8
    old = E.parse(source).getroot()
    assert old.attrib == {"version": "4.0"}
    assert old.find("Module").get("type") == "solid"
    assert old.find("Material/material").get("type") == "Mooney-Rivlin"
    assert [old.findtext("Material/material/" + k) for k in ["c1", "c2", "k"]] == [
        "20",
        "0",
        "100",
    ]
    assert (
        old.findtext("Control/time_steps") == "10"
        and old.findtext("Control/step_size") == "0.1"
    )
    m = old.find("Mesh")
    n = m.find("Nodes")
    e = m.find("Elements")
    assert e.get("type") == "hex8" and e.get("name") == "Part13"
    points = {int(x.get("id")): tuple(map(float, x.text.split(","))) for x in n}
    axes = [sorted({p[k] for p in points.values()}) for k in range(3)]
    assert list(map(len, axes)) == [21, 6, 6]
    assert points == {
        node(i, j, k, 20, 5, 5): (x, y, z)
        for i, x in enumerate(axes[0])
        for j, y in enumerate(axes[1])
        for k, z in enumerate(axes[2])
    }
    originalcells = {int(x.get("id")): list(map(int, x.text.split(","))) for x in e}
    assert originalcells == cells(20, 5, 5)

    def boundaries(p):
        ends = {n for n, q in p.items() if q[0] in [axes[0][0], axes[0][-1]]}
        return {
            "FixedDisplacement01": ends,
            "FixedDisplacement02": ends,
            "FixedDisplacement03": {n for n, q in p.items() if q[0] == axes[0][0]},
            "PrescribedBC243": {n for n, q in p.items() if q[0] == axes[0][-1]},
        }

    assert {
        s.get("name"): set(map(int, s.text.split(","))) for s in m.findall("NodeSet")
    } == boundaries(points)
    assert (
        len(m.findall("ElementSet")) == 1 and m.find("ElementSet").text.strip() == "240"
    )
    assert not any(
        "unnamed_set01" in E.tostring(s, encoding="unicode")
        for s in old
        if s.tag != "Mesh"
    )
    out = copy.deepcopy(old)
    solver = out.find("Control/solver")
    assert solver.find("linear_solver") is None
    E.SubElement(solver, "linear_solver", type="febio_gpu_lu")
    hooked = copy.deepcopy(out)
    a = [
        [
            x + (y - x) * j / factor
            for x, y in zip(axis, axis[1:])
            for j in range(factor)
        ]
        + [axis[-1]]
        for axis in axes
    ]
    nx, ny, nz = 20 * factor, 5 * factor, 5 * factor
    p = {
        node(i, j, k, nx, ny, nz): (x, y, z)
        for i, x in enumerate(a[0])
        for j, y in enumerate(a[1])
        for k, z in enumerate(a[2])
    }
    assert all(
        p[node(i * factor, j * factor, k * factor, nx, ny, nz)]
        == points[node(i, j, k, 20, 5, 5)]
        for i in range(21)
        for j in range(6)
        for k in range(6)
    )
    mesh = out.find("Mesh")
    nodes = mesh.find("Nodes")
    elems = mesh.find("Elements")
    for x in list(nodes):
        nodes.remove(x)
    for x in list(elems):
        elems.remove(x)
    for key, q in p.items():
        E.SubElement(nodes, "node", id=str(key)).text = ",".join(
            format(x, ".17g") for x in q
        )
    for key, q in cells(nx, ny, nz).items():
        E.SubElement(elems, "elem", id=str(key)).text = ",".join(map(str, q))
    bs = boundaries(p)
    for s in mesh.findall("NodeSet"):
        s.text = ",".join(map(str, sorted(bs[s.get("name")])))
    originalindex = 239
    i = originalindex // 25
    j = (originalindex % 25) // 5
    k = originalindex % 5
    children = [
        1 + ((i * factor + di) * ny + j * factor + dj) * nz + k * factor + dk
        for di in range(factor)
        for dj in range(factor)
        for dk in range(factor)
    ]
    mesh.find("ElementSet").text = ",".join(map(str, children))
    unchanged = {
        x.tag: E.tostring(x) == E.tostring(y)
        for x, y in zip(hooked, out)
        if x.tag != "Mesh"
    }
    assert all(unchanged.values())
    E.ElementTree(out).write(output, encoding="iso-8859-1", xml_declaration=True)
