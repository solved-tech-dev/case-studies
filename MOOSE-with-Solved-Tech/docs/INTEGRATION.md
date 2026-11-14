# PanguLU integration and reuse

## Native application path

`MOOSE → libMesh → PETSc KSP/PC LU → integration/petsc/solver_pangulu.cpp → adapter/pangulu_solver.hpp → PanguLU 26.46`

`app/src/base/MooseSolvedApp.C` registers original modules without replacing physics kernels. Loading `libmoose_pangulu.so` through `-dll_prepend` invokes `PetscDLLibraryRegister_moose_pangulu()`, registering the `pangulu` general-LU solver for sequential AIJ matrices. Selecting `preonly lu pangulu` returns each solution through PETSc to MOOSE; iteration control and physical output remain native application work.

| File | Responsibility |
|---|---|
| `integration/petsc/solver_pangulu.cpp` | Register factor callbacks, accept matrices and single/block RHS, return solutions and errors |
| `adapter/sparse.hpp` | Preserve original full CSR, map values into CSC, independently evaluate the original-system residual |
| `adapter/pangulu_solver.hpp` | PanguLU 26.46 ABI 2 API, device transfer, phases, reuse and native timing |
| `app/` | Register original Solid Mechanics, Heat Transfer and required Ray Tracing modules |
| `tools/run_case.py` | Prepare grids, add full-Jacobian options, run and check complete applications |

## Matrix and reuse rules

The plugin uses real FP64 general LU. It neither selects Cholesky nor extracts a matrix triangle. CSR-to-CSC conversion preserves every stored entry, including explicit zeros, without perturbing values or dropping zeros. Each numeric callback checks the dimensions, complete row pointers and column indices, and compares all coefficients exactly.

- Changed graph: release the old context, analyze and factor the new graph.
- Same graph, changed values: retain analysis, update every coefficient and refactor numerically.
- Same graph and all coefficients equal: retain numerical factors for subsequent RHS.

The SDK uses `PANGULU_REORDERING_METHOD=L1` and `PANGULU_REUSE_LEVEL=L2`. The cache belongs to its source PETSc matrix, survives recreation of its factor object, and is released with the source matrix. Separate source matrices own separate caches.

## RHS and error handling

`MatSolve` uses a one-column multiple-vector descriptor. `MatMatSolve` packs the independent dense RHS block supplied by PETSc and executes one native block SOLVE. The plugin does not invent future Newton RHS or batch application steps that depend on previous solutions.

After each RHS, an independent CPU calculation evaluates `||Ax-b||₂/||b||₂` using the full original CSR and RHS. The default finite-result threshold is `1e-8`; zero RHS uses absolute residual. Evidence time is separate. Nonzero API statuses, nonfinite results and excessive residuals propagate through PETSc errors; no silent solver fallback is used.

## Timing semantics

`MOOSE_PANGULU_NATIVE` records synchronized, non-overlapping native phases: `analysis_order_symbolic` (including the API's analysis/preprocessing work), `numeric`, `refactor_numeric`, and `solve`. Combined analysis is not split into fictitious separate ordering and symbolic measurements. Their sum is the solver core time.

Conversion, allocation, transfer, independent residuals, MOOSE assembly, nonlinear control and output are outside that sum. `application_wall_seconds` measures the complete application subprocess, already including its internal work; it must not be added to native time. Validation has no historical-time threshold.

## Supported configurations

Linux x86_64, PanguLU 26.46 ABI 2, a CUDA 12.9-compatible runtime, and real double-precision PETSc. This adapter supports one MPI rank, sequential AIJ and one GPU. It uses PETSc factor callbacks and internal matrix structures: compile against the same PETSc headers, libraries and MPI as MOOSE. Development checks target PETSc 3.25.4. Distributed matrices, complex arithmetic and transpose solves are not implemented.

PETSc's registration API is documented in [MatSolverTypeRegister](https://petsc.org/release/manualpages/Mat/MatSolverTypeRegister/).
