"""Frequency-domain Maxwell field with Robin boundary data (FEALPy curlcurl example 1)."""
import sys
from pathlib import Path

from fealpy.backend import backend_manager as bm
from fealpy.fem.curlcurl_lfem_model import CurlCurlLFEMModel

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from adapter.driver import arguments, solve, validate_field


def main():
    args = arguments(__doc__)
    bm.set_backend("numpy")
    n = 6 if args.small else 316
    model = CurlCurlLFEMModel(dict(
        pde=1, nx=n, ny=n, space_degree=1, init_mesh="uniform_tri",
        solver="direct", method="standard", wave_number=1.0, gamma=-0.07,
        pbar_log=False, log_level="WARNING",
    ))
    matrix, load = model.linear_system["standard"](model.mesh, model.p)
    field = model.space.function(dtype=bm.complex128)
    field[:] = solve(matrix, load)
    error = float(model.mesh.error(model.pde.solution, field.value))
    validate_field(field, error, "electric_field_L2_error", 0.05, args.output)


if __name__ == "__main__":
    main()
