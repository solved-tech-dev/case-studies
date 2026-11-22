# 3D linear elasticity

[简体中文](README.zh-CN.md) · [Build and run](../../docs/BUILD-AND-RUN.md) · [Entry point](run.py)

A homogeneous isotropic unit cube has zero body force and a prescribed linear displacement on all six faces. The displacement combines extension and shear. FEALPy linear elasticity example 5 supplies the boundary data, the analytical three-component displacement and the Lamé constants λ = μ = 4 × 10⁵. [Upstream model](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/model/linear_elasticity/exp0005.py)

The unit cube is partitioned into 1,024 × 17 × 17 cells and split into 1,775,616 tetrahedra. Linear elements give 332,100 nodes, 996,300 displacement unknowns and 33,228,054 CSR entries. The refinement is directional. The small mesh uses 8 × 8 × 8 subdivisions and gives 2,187 equations.

## Run and interpret the output

From the package root:

```bash
python cases/elastic3d/run.py --small
python cases/elastic3d/run.py --output output/elastic3d.npy
```

The model assembles the stiffness and load with linear_system() and applies displacement data with apply_bc(). The saved float64 vector contains interleaved x, y and z displacement components at each node. Reshape it to (number_of_nodes, 3) to recover the nodal displacement. The relative nodal error against the analytical displacement must be ≤ 0.005. Each mesh requires one RHS. The full CSR relative residual ‖Ax − b‖₂ / ‖b‖₂ must be ≤ 10⁻⁸ before field validation and export.
