# Oneel20fi

[English](README.md) | [简体中文](README.zh-CN.md)

A unit-cube heat-transfer model with C3D20R elements, conductivity 50, one face at temperature 0, and convection to ambient −12 with coefficient 10. The analytic steady solution is T=−2x, with flux (100,0,0).

## Model parameters

| Item | Value |
|---|---|
| Element type | C3D20R |
| Mesh dimension | 3D |
| Subdivision | 44 |
| Nodes | 358,425 |
| Elements | 85,184 |
| Equations n | 352,440 |
| Complete nnz | 19,849,804 |
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
python3 tools/run_case.py oneel20fi \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/oneel20fi
```

Add `--profile small` and choose another output directory to use the upstream mesh. The adapter uses L0 ordering and L2 structural reuse, updating coefficients and refactorizing on graph-unchanged calls.

See [Examples](../../docs/CASES.md) for the physical interpretation and [Checking a calculation](../../docs/VALIDATION.md) for numerical checks.
