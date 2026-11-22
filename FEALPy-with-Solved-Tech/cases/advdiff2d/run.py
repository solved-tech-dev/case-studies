"""Quadratic triangular advection–diffusion–reaction model (FEALPy Exp0001)."""
import sys
from pathlib import Path

import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.fem import (
    BilinearForm, DirichletBC, LinearForm, ScalarConvectionIntegrator,
    ScalarDiffusionIntegrator, ScalarMassIntegrator, ScalarSourceIntegrator,
)
from fealpy.functionspace import LagrangeFESpace
from fealpy.model.diffusion_convection_reaction.exp0001 import Exp0001

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from adapter.driver import arguments, solve, validate_field


def main():
    args = arguments(__doc__)
    bm.set_backend("numpy")
    n = 8 if args.small else 500
    pde = Exp0001()
    mesh = pde.init_mesh["uniform_tri"](nx=n, ny=n)
    space = LagrangeFESpace(mesh, p=2)
    form = BilinearForm(space)
    form.add_integrator(ScalarDiffusionIntegrator(coef=pde.diffusion_coef().reshape(1, 1, 2, 2), q=5))
    form.add_integrator(ScalarMassIntegrator(coef=float(np.asarray(pde.reaction_coef()).item()), q=5))
    form.add_integrator(ScalarConvectionIntegrator(coef=pde.convection_coef().reshape(1, 1, 2), q=5))
    source = LinearForm(space)
    source.add_integrator(ScalarSourceIntegrator(source=pde.source, q=5))
    bc = DirichletBC(space, gd=pde.dirichlet, threshold=pde.is_dirichlet_boundary)
    matrix, load = bc.apply(form.assembly(), source.assembly())
    field = space.function()
    field[:] = solve(matrix, load)
    exact = np.asarray(pde.solution(space.interpolation_points())).ravel()
    error = np.linalg.norm(field.array - exact) / np.linalg.norm(exact)
    validate_field(field, error, "relative_nodal_scalar_error", 0.005, args.output)


if __name__ == "__main__":
    main()
