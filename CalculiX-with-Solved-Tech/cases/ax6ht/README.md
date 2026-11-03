# Ax6ht

[English](README.md) | [简体中文](README.zh-CN.md)

Axisymmetric steady heat transfer with quadratic CAX6 triangles and conductivity 50. Nodes 1 and 46 have prescribed temperatures 0 and 100. Point-constraint localization depends on mesh refinement.

## Model parameters

| Item | Value |
|---|---|
| Element type | CAX6 |
| Mesh dimension | 2D axisymmetric |
| Subdivision | 125 |
| Nodes | 1,002,501 |
| Elements | 500,000 |
| Equations n | 1,002,499 |
| Complete nnz | 11,509,973 |
| Linear systems per application | 2 |

## Files

- [Refined input](input/model.inp.gz): gzip-compressed complete CalculiX model; decompressed by the runner.
- [Parameters](parameters.json): mesh, solver configuration, and source information.
- [Physics cards](physics-cards.inp): material, boundary, load, step, and output cards for reference.
- [Upstream input files](source/): CalculiX 2.23 verification example; [source distribution](https://www.dhondt.de/ccx_2.23.test.tar.bz2).

## Run

From `CalculiX-with-Solved-Tech/`, after [building and configuring the executable](../../docs/BUILD-AND-RUN.md):

```bash
mkdir -p _work/runs
python3 tools/run_case.py ax6ht \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/ax6ht
```

Add `--profile small` and choose another output directory to use the upstream mesh. The adapter uses L0 ordering and L2 structural reuse, updating coefficients and refactorizing on graph-unchanged calls.

See [Examples](../../docs/CASES.md) for the physical interpretation and [Checking a calculation](../../docs/VALIDATION.md) for numerical checks.
