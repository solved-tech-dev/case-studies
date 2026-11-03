# Dam

[English](README.md) | [简体中文](README.zh-CN.md)

Planar orthotropic diffusion with CPS8 elements, representing a groundwater/scalar-diffusion analogy. The two prescribed boundary values are 8.5344 and 3.9624.

## Model parameters

| Item | Value |
|---|---|
| Element type | CPS8 |
| Mesh dimension | 2D |
| Subdivision | 32 |
| Nodes | 985,756 |
| Elements | 327,680 |
| Equations n | 984,474 |
| Complete nnz | 15,390,186 |
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
python3 tools/run_case.py dam \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/dam
```

Add `--profile small` and choose another output directory to use the upstream mesh. The adapter uses L0 ordering and L2 structural reuse, updating coefficients and refactorizing on graph-unchanged calls.

See [Examples](../../docs/CASES.md) for the physical interpretation and [Checking a calculation](../../docs/VALIDATION.md) for numerical checks.
