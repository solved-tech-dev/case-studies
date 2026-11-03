# Code organization and attribution

[English](OWNERSHIP.md) | [简体中文](OWNERSHIP.zh-CN.md) · [Home](../README.md)

## Components

| Component | Implementation and location |
|---|---|
| Finite elements, materials, assembly, constraints, Newton control, field output | CalculiX 2.23; source archive under `provenance/upstream/` |
| Upstream example inputs | Four CalculiX verification examples and the community Pan example; `cases/<case>/source/` |
| Refined meshes and boundary mapping | Project mesh-preparation tool, `tools/generate_mesh.py` |
| PanguLU solver-entry implementation | `adapter/ccx_pangulu.cpp` |
| CalculiX entry declarations | `adapter/spooles_entry.h` |
| Ordering, analysis, numerical factorization, and solve algorithms | External PanguLU library |
| Build, run, residual, and field checks | `tools/` and `tests/` |
| Newton-stage timing | Optional `patches/nonlingeo-timing.patch` |

## Build configuration

The build tool compiles CalculiX 2.23 with the PanguLU entry implementation, using 32-bit ITG and LP64 BLAS/LAPACK. It enables the CalculiX SPOOLES call branch, ARPACK, MATRIXSTORAGE, and NETWORKOUT.

The application-side `spooles.h` declarations are supplied by the adapter header. The `spooles.c` solver wrapper is replaced by the adapter object. Constraint-processing helpers use the external SPOOLES library, and `dgesv` is supplied by external LAPACK. The executable is linked from the main object and a CalculiX source archive.

The optional timing patch adds measurements around Newton stages.

## Model preparation

`source/` contains the upstream input models. `input/model.inp.gz` contains the refined model used by the runner. Refinement subdivides the parent elements, maps sets and boundary faces to child elements, and enables full-field output.

Material parameters, physical loads, boundary values, and step controls are specified in each case's input and `physics-cards.inp`. In Ax6ht, the `STEP,NLGEOM` setting and prescribed temperatures at nodes 1 and 46 are part of the model.

## Attribution and licenses

CalculiX and upstream examples retain their authorship and license notices. The adapter and preparation tools are project components from 解得科技. PanguLU is installed and licensed separately.

See [License notice](../LICENSE-NOTICE.md) and the source URLs in each case's `parameters.json`.

Project additions by 解得科技 are licensed under AGPL-3.0-only; see [License notice](../LICENSE-NOTICE.md).
