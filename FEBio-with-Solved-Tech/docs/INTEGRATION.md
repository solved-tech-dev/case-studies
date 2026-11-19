# How the FEBio plugin works

[简体中文](INTEGRATION.zh-CN.md)

The model selects `febio_gpu_lu` in `Control/solver/linear_solver`. The `PluginInitialize` entry in `FEBio-with-Solved-Tech/adapter/febio_gpu.cpp` registers the `FebioGpuLU` class. FEBio loads `FEBio-with-Solved-Tech/build/pangulu.so` via `-import`.

| FEBio entry | Action |
| --- | --- |
| CreateSparseMatrix / SetSparseMatrix | Provide or accept zero-based full CSR storage for assembly by FEBio. |
| PreProcess / Factor | Compare row pointers and column indices exactly. Analyze and factor a new pattern. For an identical pattern, refresh every coefficient and refactor if values changed. Reuse the numeric factors when coefficient bytes also match. |
| BackSolve | Solve the current FEBio right-hand side, return the solution and check the residual using the original CSR and a copy of the right-hand side. |
| Destroy | Mark the system unready. Factor checks the next assembled pattern before reuse. Solver destruction releases the backend resources. |

`FEBio-with-Solved-Tech/adapter/gpu_backend.inc` converts CSR to CSC, preserving explicit zeros and a map back to the original coefficients. It uses double-precision general LU, L0 ordering and L2 reuse through the native `pangulu_execute` phases. FEBio supplies one right-hand side at a time through BackSolve, so the vector descriptor has one column.
