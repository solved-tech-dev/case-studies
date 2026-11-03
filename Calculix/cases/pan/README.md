# Pan

[English](README.md) | [简体中文](README.zh-CN.md)

An axisymmetric CAX6 iron-pan model with conductivity 50, bottom flux 31500, air/food convection coefficients 10/150, and ambient temperature 293.15. Heat balance uses full-revolution 2πr weighting.

## Model parameters

| Item | Value |
|---|---|
| Element type | CAX6 |
| Mesh dimension | 2D axisymmetric |
| Subdivision | 16 |
| Nodes | 890,630 |
| Elements | 442,880 |
| Equations n | 890,630 |
| Complete nnz | 10,206,762 |
| Linear systems per application | 2 |

## Files

- [Refined input](input/model.inp.gz): gzip-compressed complete CalculiX model; decompressed by the runner.
- [Parameters](parameters.json): mesh, solver configuration, and source information.
- [Physics cards](physics-cards.inp): material, boundary, load, step, and output cards for reference.
- [Upstream input files](source/): CalculiX community Launcher example; [source distribution](https://github.com/calculix/examples/tree/a079253398c05f96b019e8c07ce076e370ddeb75/launcher/Therm/pan).

## Run

From `Calculix/`, after [building and configuring the executable](../../docs/BUILD-AND-RUN.md):

```bash
mkdir -p _work/runs
python3 tools/run_case.py pan \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/pan
```

Add `--profile small` and choose another output directory to use the upstream mesh. The adapter uses L0 ordering and L2 structural reuse, updating coefficients and refactorizing on graph-unchanged calls.

See [Examples](../../docs/CASES.md) for the physical interpretation and [Checking a calculation](../../docs/VALIDATION.md) for numerical checks.
