#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Generate the five original ngspice operating-point netlists."""
import argparse
import hashlib
import json
from pathlib import Path


def tile_network(nx, ny, width, height, scheme):
    yield f"* {scheme} power distribution network with resistive device loads\n"
    yield ".options klu\nVDD supply 0 DC 1.05\n"
    for y in range(ny):
        for x in range(nx):
            here = f"n{x}_{y}"
            tile_x, tile_y = x // width, y // height
            local_x, local_y = x % width, y % height
            if scheme == "interleaved":
                h = local_x < width-1 or local_y in (3, height//2)
                v = local_y < height-1 or local_x in (4, width//2, width-1)
                v = v and (local_x % (2 + tile_x % 3) == 0 or local_x in (4, width//2, width-1))
            else:
                h = local_x < width-1 or local_y % 12 in (0, 1)
                v = (local_y < height-1 or local_x % 12 in (0, 1)) and (x+y//height) % 2 == 0
            if h and x+1 < nx:
                yield f"RH{here} {here} n{x+1}_{y} {0.045+0.015*(y%4):.3f}\n"
            if v and y+1 < ny:
                yield f"RV{here} {here} n{x}_{y+1} {0.050+0.020*(x%3):.3f}\n"
            if local_y == 0 and local_x == 0:
                yield f"RF{here} supply {here} 0.008\n"
            yield f"RL{here} {here} 0 {50000+500*((17*x+23*y)%97)}\n"
    yield ".op\n.save all\n.end\n"


def strapped_network(name, nx, ny, diagonal):
    yield f"* {name} DC power distribution with resistive load and supply paths\n"
    yield ".options klu\nVDD vdd 0 DC 1.0\n"
    for y in range(ny):
        for x in range(nx):
            here = f"n{x}_{y}_0"
            if x+1 < nx:
                yield f"RH{here} {here} n{x+1}_{y}_0 {0.04+0.02*(y%3):.3f}\n"
            if y+1 < ny and x%4==0:
                yield f"RV{here} {here} n{x}_{y+1}_0 {0.05+0.015*(x%4):.3f}\n"
            if y%13==0 and x+17<nx:
                yield f"RS{here} {here} n{x+17}_{y}_0 0.12\n"
            if diagonal and x+1<nx and y+1<ny and (x+y)%3==0:
                yield f"RD{here} {here} n{x+1}_{y+1}_0 0.18\n"
            if x%32==0 and y%32==0:
                yield f"RF{here} vdd {here} 0.04\n"
            conductance=(x*19+y*37)%113
            yield f"RL{here} {here} 0 {100000+1000*conductance}\n"
    yield ".op\n.save all\n.end\n"


def two_rail(nx,ny):
    yield "* Two coupled rails with resistive differential loads\n"
    yield ".options klu\nVDD vdd 0 DC 1.1\n"
    for y in range(ny):
        for x in range(nx):
            here=f"n{x}_{y}"
            if x+1<nx:
                yield f"RPH{here} p{here} pn{x+1}_{y} 0.075\n"
                yield f"RNH{here} m{here} mn{x+1}_{y} 0.055\n"
            if y+1<ny and x%2==0:
                yield f"RPV{here} p{here} pn{x}_{y+1} 0.085\n"
            if y+1<ny and x%3==0:
                yield f"RNV{here} m{here} mn{x}_{y+1} 0.065\n"
            # Resistor and independent current sink model the DC load.
            yield f"RL{here} p{here} m{here} {750+(7*x+13*y)%1500}\n"
            yield f"IL{here} p{here} m{here} DC {5+(11*x+17*y)%29}u\n"
            if x%12==0 and y%12==0:
                yield f"RBP{here} vdd p{here} 0.03\n"
                yield f"RBN{here} m{here} 0 0.025\n"
    yield ".op\n.save all\n.end\n"


def main():
    root = Path(__file__).resolve().parents[1]
    cases = json.loads((root / "cases/manifest.json").read_text())["cases"]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=[*cases, "all"])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for name in cases if args.case == "all" else [args.case]:
        case = cases[name]
        nx, ny = case["grid"]
        if "tile" in case:
            lines = tile_network(nx, ny, *case["tile"], case["topology"])
        elif case["topology"] == "dual-rail":
            lines = two_rail(nx, ny)
        else:
            lines = strapped_network(name, nx, ny, case["topology"] == "diagonal")
        path = args.output_dir / (name + ".cir")
        with path.open("x", encoding="ascii", newline="\n") as stream:
            stream.writelines(lines)
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != case["netlist_sha256"]:
            raise RuntimeError(f"{path}: generated netlist differs from the tested input")
        print(path, digest, flush=True)


if __name__ == "__main__":
    main()
