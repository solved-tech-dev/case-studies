# PHG with PanguLU

[English](README.md) · [简体中文](README.zh-CN.md)

Run three native finite-element cases with PanguLU 26.44. PHG and PHG-Solid handle the mesh, material model, boundary conditions, assembly, and output. The added `pangulu` backend passes the assembled double-precision system to PanguLU and returns the solution to the original application.

| Case | Upstream entry point | Matrix order | Linear solves |
| --- | --- | ---: | ---: |
| `solid-adaptive-stress` | PHG-Solid `examples/adaptive.xml` | 421,443 → 1,268,499 | 6 |
| `solid-static-p2` | PHG-Solid `examples/static.xml` | 833,091 | 1 |
| `interface-sphere-dg1` | PHG `examples/interface.c` | 786,432 | 1 |

Each `cases/<case>/case.json` records the parameters and matrix sizes. The setup tool obtains the pinned upstream sources and meshes, and the runner writes the selected configuration to its own output directory.

Follow [Build and run](docs/BUILD-AND-RUN.md) to install the dependencies and build the two native executables. Then run a case from the project root.

```bash
python3 tools/run_case.py solid-static-p2 --output _work/runs/static-p2 --vtk
```

PanguLU requires a separately obtained SDK and a license for the target GPU. [Solved Tech](https://solved-tech.com/)

[Case details](docs/CASES.md) · [Integration and outputs](docs/INTEGRATION.md) · [License and attribution](LICENSE-NOTICE.md)
