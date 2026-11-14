# MOOSE with Solved Tech

[中文说明](README.zh-CN.md)

Run six MOOSE finite-element applications with PanguLU 26.46 through a native PETSc direct-solver plugin. The original physics objects remain in MOOSE; this package provides the solver connection, original case inputs and seed meshes, reproducible grid preparation, and complete application runners.

| Case | What it computes | Enlarged matrix rows |
|---|---|---:|
| `thermomech` | Temperature, thermal expansion and displacement over five time steps | 1,002,501 |
| `annular-torsion` | Displacement and stress in a hollow cylinder under torque | 1,005,312 |
| `mug` | Steady scalar diffusion through a cup-shaped domain | 1,349,545 |
| `pinched-cylinder` | Displacements and rotations of a thin cylindrical shell under point loading | 1,017,005 |
| `gap-heat` | Heat transfer between two solids across a conducting gap | 982,802 |
| `hex20-heat` | A quadratic-element heat-conduction patch problem | 988,281 |

These sizes describe the supplied `million` profiles, with equation numbering from the validated source configuration. A `small` profile is also provided for every case. Both use the same material, load and time-control definitions. The HEX20 enlarged mesh is a structured unit cube replacing the original distorted seed mesh; see the case documentation.

## Start here

1. Read [build and run](docs/BUILD-AND-RUN.md) for the MOOSE, PETSc and PanguLU prerequisites and exact commands.
2. Read [cases and grid preparation](docs/CASES.md) for physical meaning and the changes applied to upstream inputs.
3. Read [integration](docs/INTEGRATION.md) for native callbacks, exact-zero handling, reuse and timing definitions.
4. Read [validation](docs/VALIDATION.md) for residuals, physical fields and reference comparisons.

The main code entry is `integration/petsc/solver_pangulu.cpp`; `adapter/pangulu_solver.hpp` wraps the PanguLU 26.46 API. `app/` registers the original Solid Mechanics and Heat Transfer modules. `tools/run_case.py` launches full MOOSE applications, not matrix-replay benchmarks.

The current plugin supports **one MPI rank, one NVIDIA GPU, real FP64 sequential AIJ matrices and general LU**. Independent RHS blocks supplied through PETSc `MatMatSolve` use the native PanguLU multiple-vector API. Newton-dependent RHS are solved in their original application sequence.

This source distribution contains no measured benchmark tables, performance figures, historical simulation results, SDKs or license keys. Running a case produces physical output and reports its measured phase times; there are no expected-time thresholds.

See [license and attribution](LICENSE-NOTICE.md). Upstream input files are byte-for-byte preserved and identified in `upstream.lock.json`.
