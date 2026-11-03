# CalculiX + PanguLU

[English](README.md) | [简体中文](README.zh-CN.md)

Run five CalculiX finite element examples with PanguLU: **Dam, Ax6ht, Capacitor, Oneel20fi, and Pan**.

The integration connects CalculiX 2.23 to PanguLU at the linear-system solver interface. CalculiX handles finite element assembly, materials, boundary conditions, Newton iterations, and result output; the adapter passes the assembled system to PanguLU and returns the solution.

PanguLU is a solver product from 解得科技. To purchase PanguLU and obtain a license, visit the [company website](https://solved-tech.com/). The PanguLU SDK and license are installed separately.

The adapter, project tools, tests, and documentation authored by 解得科技 use [AGPL-3.0-only](LICENSE). CalculiX and upstream examples retain their respective copyrights and licenses; see [License notice](LICENSE-NOTICE.md).

## Documentation

- [Build and run](docs/BUILD-AND-RUN.md): dependencies, compilation, and example commands.
- [Integration](docs/INTEGRATION.md): solver entry points, data conversion, and structural reuse.
- [Examples](docs/CASES.md): physical models and mesh parameters.
- [Checking a calculation](docs/VALIDATION.md): residual checks and physical-field validation.
- [Code organization](docs/OWNERSHIP.md): CalculiX components and integration components.
- [License notice](LICENSE-NOTICE.md).

## Directory layout

```text
adapter/                  CalculiX–PanguLU adapter
cases/<case>/source/      Upstream example inputs and include files
cases/<case>/input/       Refined input model, compressed with gzip
cases/<case>/parameters.json  Mesh parameters and system dimensions
provenance/upstream/      CalculiX 2.23 source archive and license
patches/                 Optional Newton-stage timing patch
tools/                   Build, run, validation, and mesh-generation tools
tests/                   Input-preparation and coefficient-update checks
docs/                    English and Chinese documentation
```

## Quick start

From the `case-studies` repository root, enter `Calculix/` and follow the [build instructions](docs/BUILD-AND-RUN.md). Once the executable and runtime environment are configured:

```bash
mkdir -p _work/runs
python3 tools/run_case.py dam \
  --binary _work/build/ccx_pangulu \
  --output _work/runs/dam
```

Replace `dam` with `ax6ht`, `capacitor`, `oneel20fi`, or `pan` to run another example. Results are written to the selected output directory.

## Supported configuration

The examples use real, symmetric, unshifted steady-state scalar systems on Linux with one process and one NVIDIA GPU. The adapter uses PanguLU 26.42, ABI2, CUDA 12.9, L0 ordering, and L2 structural reuse. Graph-unchanged calls refresh the matrix coefficients and perform numerical refactorization.

The provided models cover planar diffusion, axisymmetric heat transfer, and three-dimensional scalar problems. See the example descriptions for the geometry and physical interpretation of each model.
