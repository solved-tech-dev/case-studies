# FEALPy with PanguLU

[简体中文](README.zh-CN.md) · [Build and run](docs/BUILD-AND-RUN.md) · [Solver integration](docs/INTEGRATION.md)

Compute displacement, a frequency-domain electric field, and a steady transported scalar with FEALPy and PanguLU 26.46. Each example builds its mesh and sparse system from the original FEALPy model. The adapter performs one general LU factorization and one single-RHS solve on a CUDA GPU. Each script checks both the original equations and the finite-element field.

| Case and instructions | Discretization | Equations | Stored CSR entries |
| --- | --- | ---: | ---: |
| [3D linear elasticity](cases/elastic3d/README.md) | Linear tetrahedra, three displacement components | 996,300 | 33,228,054 |
| [2D frequency-domain Maxwell](cases/maxwell2d/README.md) | First-family Nédélec elements, complex coefficients | 999,824 | 11,585,824 |
| [2D advection–diffusion–reaction](cases/advdiff2d/README.md) | Quadratic Lagrange triangles | 1,002,001 | 11,456,053 |

Follow the build guide to install the matching FEALPy revision and build the real and complex adapters. Run any case with `--small` for a quick check or use its default mesh to reproduce the system in the table. The original FEALPy checkout needs no edits. [Sources and licensing](NOTICE.md)
