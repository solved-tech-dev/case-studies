"""Prescribed displacement on a tetrahedral unit cube (FEALPy Exp0005)."""
import sys
from pathlib import Path

import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.fem.linear_elasticity_lfem_model import LinearElasticityLFEMModel

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from adapter.driver import arguments, solve, validate_field


def main():
    args = arguments(__doc__)
    bm.set_backend("numpy")
    model = LinearElasticityLFEMModel()
    model.set_pde(5)
    nx, ny, nz = (8, 8, 8) if args.small else (1024, 17, 17)
    model.mesh = model.pde.init_mesh["uniform_tet"](nx=nx, ny=ny, nz=nz)
    model.set_space_degree(p=1)
    matrix, load = model.linear_system(model.mesh, model.p)
    matrix, load = model.apply_bc(matrix, load)
    model.uh[:] = solve(matrix, load)
    points = model.space.interpolation_points()
    exact = np.asarray(model.pde.displacement(points)).ravel()
    error = np.linalg.norm(model.uh.array - exact) / np.linalg.norm(exact)
    validate_field(model.uh, error, "relative_nodal_displacement_error", 0.005, args.output)


if __name__ == "__main__":
    main()
