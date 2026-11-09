# PHG cases

[简体中文](CASES.zh-CN.md)

## Adaptive solid stress

PHG-Solid's `examples/adaptive.xml` computes the static response of a loaded cantilever using first-order tetrahedral displacement elements. The upstream input defines the fixed end and the surface load. After an initial uniform mesh refinement, the application uses a residual-based error estimator to decide where to refine the mesh locally. The configuration sets a budget of 1,000,000 degrees of freedom. The application checks that budget after completing each assembly and solve, so the final matrix has order 1,268,499.

| Adaptive step | Matrix order | Stored CSR entries |
| ---: | ---: | ---: |
| 1 | 421,443 | 14,740,623 |
| 2 | 452,952 | 15,857,580 |
| 3 | 530,697 | 18,733,887 |
| 4 | 677,124 | 25,224,914 |
| 5 | 895,794 | 33,405,310 |
| 6 | 1,268,499 | 47,802,493 |

## Linear static cantilever with quadratic elements

PHG-Solid's `examples/static.xml` uses quadratic tetrahedral displacement elements. The configuration starts from the upstream `cantilever.mesh` and applies 13 levels of initial refinement. It assembles and solves one system on the resulting fixed mesh. The matrix has order 833,091 and 61,575,691 stored CSR entries.

## Spherical material interface (DG1)

PHG's `examples/interface.c` defines a spherical interface of radius 0.3 centered at `(0.5, 0.5, 0.5)` in a unit cube. The diffusion coefficient is 0.5 inside the sphere and 2.0 outside. The upstream program uses an extended finite element formulation to handle jumps in the solution and flux at the interface. This configuration uses the tetrahedral mesh `test/cube.dat` supplied with PHG, DG1 elements, and 14 levels of initial refinement. It assembles and solves one system on the resulting fixed mesh. The matrix has order 786,432 and 8,625,216 stored CSR entries.

The physical parameters come from the official upstream archives. See `cases/<case>/case.json` for the settings applied at run time. The CSR entry counts include explicit zeros.
