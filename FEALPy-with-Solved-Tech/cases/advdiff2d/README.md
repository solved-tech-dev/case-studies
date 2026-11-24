# 2D steady advection–diffusion–reaction

[简体中文](README.zh-CN.md) · [Build and run](../../docs/BUILD-AND-RUN.md) · [Entry point](run.py)

On the unit square the diffusion tensor is diag(2, 3), the velocity is (1, −1) and the reaction coefficient is 4. Homogeneous Dirichlet boundary data and a manufactured source give u = sin(πx) sin(πy). The advection term produces a nonsymmetric matrix. [Upstream model](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/model/diffusion_convection_reaction/exp0001.py)

The 500 × 500 grid is split into 500,000 triangles. Quadratic Lagrange elements have degrees of freedom at vertices and edge midpoints. The resulting matrix has 1,002,001 rows and 11,456,053 CSR entries. The small mesh uses an 8 × 8 partition and gives 289 equations.

## Run and interpret the output

From the package root:

```bash
python cases/advdiff2d/run.py --small
python cases/advdiff2d/run.py --output output/advdiff2d.npy
```

FEALPy assembles diffusion, mass, convection and load using quadrature parameter q = 5, then applies the Dirichlet conditions. The saved float64 vector gives the scalar at the interpolation points of the P2 space. Its relative nodal error against the analytical field must be ≤ 0.005. Each mesh requires one RHS. The full CSR relative residual ‖Ax − b‖₂ / ‖b‖₂ must be ≤ 10⁻⁸ before field validation and export.
