# 2D frequency-domain electric field

[简体中文](README.zh-CN.md) · [Build and run](../../docs/BUILD-AND-RUN.md) · [Entry point](run.py)

The unit-square model solves curl curl E − k²E = f at k = 1. A complex Robin term on the boundary represents an impedance condition. The volume and boundary sources are taken from FEALPy curlcurl example 1, with analytical components x y (1 − x)(1 − y) and sin(πx) sin(πy). [Upstream model](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/model/curlcurl/exp0001.py)

A 316 × 316 partition produces 199,712 triangles. The standard first-family Nédélec formulation (FEALPy parameter p = 1) gives 999,824 unknowns and 11,585,824 CSR entries. The small mesh uses a 6 × 6 partition and gives 384 equations.

## Run and interpret the output

From the package root:

```bash
python cases/maxwell2d/run.py --small
python cases/maxwell2d/run.py --output output/maxwell2d.npy
```

The matrix and right-hand side are complex128. The saved complex vector contains finite-element coefficients. Reconstruct the field with model.space.function(dtype=bm.complex128) and evaluate its value on cells with the original basis. FEALPy integrates the absolute L² error against the analytical vector field. The accepted bound is 0.05. Each mesh requires one RHS. The full CSR relative residual ‖Ax − b‖₂ / ‖b‖₂ must be ≤ 10⁻⁸ before field validation and export.
