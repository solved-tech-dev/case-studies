# Capacitor

[English](README.md) | [简体中文](README.zh-CN.md)

A three-dimensional C3D20 coaxial sector using scalar diffusion as a potential analogy. The coefficient is 8.8541878176e-12, with inner and outer values2 and 1.

## Model parameters

| Item | Value |
|---|---|
| Element type | C3D20 |
| Mesh dimension | 3D |
| Subdivision | 3 × 16 × 400 |
| Nodes | 969,687 |
| Elements | 192,000 |
| Equations n | 952,073 |
| Complete nnz | 47,134,131 |
| Linear systems per application | 1 |

## Files

- [Refined input](input/model.inp.gz): gzip-compressed complete CalculiX model; decompressed by the runner.
- [Parameters](parameters.json): mesh, solver configuration, and source information.
- [Physics cards](physics-cards.inp): material, boundary, load, step, and output cards for reference.
- [Upstream input files](source/): CalculiX 2.23 verification example; [source distribution](https://www.dhondt.de/ccx_2.23.test.tar.bz2).

## Run

From `CalculiX-with-Solved-Tech/`, after [building and configuring the executable](../../docs/BUILD-AND-RUN.md):

```bash
mkdir -p _work/runs
python3 tools/run_case.py capacitor \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/capacitor
```

Add `--profile small` and choose another output directory to use the upstream mesh. The adapter uses L0 ordering and L2 structural reuse, updating coefficients and refactorizing on graph-unchanged calls.

See [Examples](../../docs/CASES.md) for the physical interpretation and [Checking a calculation](../../docs/VALIDATION.md) for numerical checks.
