#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Run a complete native MOOSE case using PanguLU or a PETSc reference solver."""

import argparse, json, os, subprocess, time
from pathlib import Path
from refine_quad import refine
from fix_gap_boundary import fix

ROOT = Path(__file__).resolve().parents[1]


def run(cmd, folder, label):
    (folder / (label + "-command.json")).write_text(json.dumps(cmd, indent=2) + "\n")
    with (folder / (label + ".log")).open("w") as log:
        result = subprocess.run(cmd, cwd=folder, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(
            f'{label} exited {result.returncode}; see {folder/(label+".log")}'
        )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "case", choices=sorted(x.name for x in (ROOT / "cases").iterdir() if x.is_dir())
    )
    p.add_argument("--profile", choices=["small", "million"], default="small")
    p.add_argument("--exe", required=True)
    p.add_argument("--plugin", default=str(ROOT / "build/libmoose_pangulu.so"))
    p.add_argument("--solver", choices=["pangulu", "reference"], default="pangulu")
    p.add_argument(
        "--reference-solver",
        default="petsc",
        help="A real double LU solver available in your PETSc, e.g. petsc or mumps",
    )
    p.add_argument("--output", required=True)
    p.add_argument("--prepare-only", action="store_true")
    a = p.parse_args()
    exe = Path(a.exe).expanduser().resolve()
    folder = Path(a.output).expanduser().resolve()
    if not exe.is_file():
        p.error("MOOSE executable does not exist")
    if (
        a.solver == "pangulu"
        and not a.prepare_only
        and not Path(a.plugin).resolve().is_file()
    ):
        p.error("Build the PanguLU plugin first")
    folder.mkdir(parents=True, exist_ok=False)
    case_dir = ROOT / "cases" / a.case
    case = json.loads((case_dir / "case.json").read_text())
    grid = case["grids"][a.profile]
    upstream = case_dir / "upstream"
    inp = upstream / case["input"]
    args = []
    # These are the same mesh, coupling and output changes used for the enlarged cases.
    # The original input stays byte-for-byte intact in upstream/.
    if a.case == "thermomech":
        args = [
            f'Mesh/generated/nx={grid["nx"]}',
            f'Mesh/generated/ny={grid["ny"]}',
            "Physics/SolidMechanics/QuasiStatic/all/temperature=T",
            "Physics/SolidMechanics/QuasiStatic/all/use_finite_deform_jacobian=true",
            "Physics/SolidMechanics/QuasiStatic/all/eigenstrain_names=thermal_expansion",
        ]
    elif a.case == "annular-torsion":
        args = [
            f'Mesh/ring/nr={grid["nr"]}',
            f'Mesh/ring/nt={grid["nt"]}',
            f'Mesh/extrude/num_layers={grid["layers"]}',
            "Physics/SolidMechanics/QuasiStatic/all/generate_output=vonmises_stress stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz",
        ]
    elif a.case == "mug":
        args = [
            f'Mesh/file={upstream/"mug.e"}',
            f'Mesh/uniform_refine={grid["refine"]}',
        ]
    elif a.case in ("pinched-cylinder", "gap-heat"):
        seed = upstream / (
            "cyl_sym_10x10.e" if a.case == "pinched-cylinder" else "meshed_gap.e"
        )
        report = refine(seed, folder / "refined.msh", grid["factor"])
        (folder / "mesh.json").write_text(json.dumps(report, indent=2) + "\n")
        key = "Mesh/mesh/file" if a.case == "pinched-cylinder" else "Mesh/fmesh/file"
        run(
            [
                str(exe),
                "-i",
                str(inp),
                f'{key}={folder/"refined.msh"}',
                "--mesh-only",
                str(folder / "mesh.e"),
            ],
            folder,
            "mesh",
        )
        mesh = folder / "mesh.e"
        if a.case == "gap-heat":
            mesh = folder / "mesh-active-boundaries.e"
            report = fix(folder / "mesh.e", mesh, grid["factor"])
            (folder / "contact-boundaries.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )
        args = [f"{key}={mesh}"]
        if a.case == "gap-heat" and a.profile == "small":
            args.append(f"Executioner/nl_max_its={grid['nl_max_its']}")
    else:
        n = grid["cells"]
        mesh_input = folder / "mesh.i"
        mesh_input.write_text(f"""[Mesh]
  [cube]
    type = GeneratedMeshGenerator
    dim = 3
    nx = {n}
    ny = {n}
    nz = {n}
    elem_type = HEX20
    subdomain_ids = '1'
  []
  [exterior]
    type = RenameBoundaryGenerator
    input = cube
    old_boundary = 'left right bottom top front back'
    new_boundary = '10 10 10 10 10 10'
  []
[]
""")
        run(
            [str(exe), "-i", str(mesh_input), "--mesh-only", str(folder / "mesh.e")],
            folder,
            "mesh",
        )
        args = [f'Mesh/file={folder/"mesh.e"}']
    options = folder / "direct.i"
    solver = "pangulu" if a.solver == "pangulu" else a.reference_solver
    # SMP/full supplies the actual full Jacobian, not a diagonal preconditioner.
    options.write_text(f"""[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]
[Executioner]
  solve_type = NEWTON
  petsc_options = ''
  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'preonly lu {solver}'
[]
[Outputs]
  exodus = true
  csv = true
  execute_on = 'INITIAL TIMESTEP_END'
[]
""")
    cmd = [
        str(exe),
        "-i",
        str(inp),
        str(options),
        "--n-threads=1",
        "--disable-thread-safe-output",
    ] + args
    if a.case == "gap-heat":
        cmd += ["Outputs/exodus=false"]  # original named Exodus output is retained
    cmd += [
        f'Outputs/file_base={folder/"physics"}',
        "-snes_converged_reason",
        "-snes_monitor",
        "-ksp_converged_reason",
    ]
    if a.solver == "pangulu":
        cmd += ["-dll_prepend", str(Path(a.plugin).resolve())]
    manifest = dict(
        case=a.case,
        profile=a.profile,
        solver=solver,
        executable=str(exe),
        command=cmd,
        case_config=case,
    )
    (folder / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if a.prepare_only:
        print(folder)
        return
    start = time.monotonic()
    run(cmd, folder, "application")
    from validate import validate

    elapsed = time.monotonic() - start
    report = validate(folder)
    report["application_wall_seconds"] = elapsed
    (folder / "validation.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
