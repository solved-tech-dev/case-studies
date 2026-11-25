# FrontISTR with PanguLU

[简体中文](README.zh-CN.md) · [Build and run](docs/BUILD.md) · [Sources and licenses](NOTICE.md)

Solve the linear systems assembled by FrontISTR 5.9 with PanguLU 26.48 on a CUDA GPU. FrontISTR retains the material models, boundary conditions and native field output. The package provides the exact refined meshes used in the three examples.

| Case | Nodes | Elements | Equations | Final field |
|---|---:|---:|---:|---|
| [Connecting rod under prescribed displacement](cases/conrod) | 333,340 | 226,614 quadratic tetrahedra | 1,000,020 | [Displacement](cases/conrod/field.png), nodal stresses |
| [Elastic hinge under nodal load](cases/hinge) | 333,340 | 227,819 quadratic tetrahedra | 1,000,020 | [Displacement](cases/hinge/field.png), nodal stresses |
| [Heat conduction in an aluminium block](cases/heat) | 914,550 | 868,320 hexahedra | 914,550 | [Temperature](cases/heat/field.png) |

The [Fortran patch](patches/frontistr-pangulu.patch) dispatches to the GPU at the direct-solver entry used by the examples. The [COO bridge](adapter/bridge.cpp) expands symmetric storage, merges duplicates and checks matrix structure and values for reuse. The [GPU backend](adapter/gpu_backend.inc) builds CSC storage and applies PanguLU analysis, factorization or refactorization, and solve. The supported path uses one MPI process, one CUDA GPU, real double precision and one right-hand side per call.

The examples default to 4 OpenMP threads in FrontISTR and 40 host threads during each PanguLU call, matching the A100 PCIe 40 GB measurements with the CUDA 12.9 SDK. Set `FRONTISTR_PANGULU_THREADS` to the number of solver threads available to your process.

Each run verifies the unpacked mesh, original-equation residuals, solve count, reuse decisions and final nodal field norms. The heat example performs six linear solves with one structural analysis. See the [build guide](docs/BUILD.md) for the source revision, configuration and commands.
